require 'yaml'
include Math

class CalcMoment #名前MomentMethodは使われてる？
  #元の定数違う，多分間違い  BOLTZ = 1.380658e-16
  BOLTZ = 1.38064852e-16
  PLANCK = 6.626e-27
  HBAR = PLANCK/(2.0*PI)
  AVOGADO = 6.023e23

  def initialize(structure)
    puts "Hi,moment"
    @structure=structure
    puts "structure = #{@structure}"
    select
    @structure == "vasp" ? calc_moment_vasp : calc_moment
  end

  def check(*val)
    p val.join(", ")
  end

  def select(file='POTCAR')
    src = YAML.load_file(file)
    #src = YAML.load_file(POTCAR)
    p @@potential=case src[:type]
    #ここの微分式POTCARのほうにいれたい．メソッドそのままPOTCARに入れるのもありかも．．
    when 'lj_jindo'
      DiffLjJindo.new(src)
    else
      p "*********************missed potcar*******************"
      exit(0)
    end
  end

  def calc_moment
    p "start calc_moment"
    a0 = calc_a0(@@potential.m, @@potential.n, @@potential.r0)
    p "a0="
    check(a0)
    for i in 1..10
      temp = 100*(i-1)
      temp = 10 if i==1
      p theta = BOLTZ*temp
      gap = Array.new(6){ Array.new(201) }
      aa1,aa2=[],[]
      aa1[0]=a0
      comp_gap = 0
      for change in 1..5 #width to displace
        for same in 1..200
          gap[change][same] = (same-1)*(10**(-(9+change))).to_f #start 1e-10
          aa1[change] = aa1[change-1] + gap[change][same]
          aa2[change] = aa1[change]*sqrt(2)
          a1, a2 = aa1[change], aa2[change]
          k = calc_k(a1,a2,@structure)
          atom_mass = @@potential.atom_mass/AVOGADO
          omega = sqrt(k/atom_mass)
          x = HBAR*omega/(2.0*theta)
          gamma = calc_gamma(a1,a2,@structure)
          gt_k2 = gamma*theta/k**2
          u0 = calc_u0(a1,a2)
          psi0 = calc_psi0(x,theta)
          large_a = calc_large_a(x, gt_k2)
          y0 = calc_y0(k, gamma, theta, large_a)
          case @structure
          when "fitting_test","vasp"
            psi_nonli = calc_psi_linear(k, x, gamma, theta)#とりあえず変数はそのままnonli
          else
            psi_nonli = calc_psi_nonli(k, x, gamma,theta) #nonlinear fcc thesis(1988)p516(19)
          end
          total_gap = comp_gap + gap[change][same]
          break if y0 - total_gap < 0
        end
        aa1[change] = aa1[change] - (10**(-(9+change))).to_f
        aa2[change] = aa1[change]*sqrt(2)
        comp_gap += gap[change][same-1]
        a1_cal = a0+comp_gap
        a2_cal = sqrt(2)*a1_cal
      end
      #ここから出力
      a1_cal = a1_cal*1.0e8
      a2_cal = a2_cal*1.0e8
      puts "T(K), a1, a2_cal(A), k, g="
      check(temp, a1_cal, a2_cal, k, gamma)
      puts "temp, u0, harmonic free, free="
      check(temp, u0, psi0, psi_nonli)
      #eVに単位変換
      u0_ev = ev_from_kelvin(u0)
      psi0_ev = ev_from_erg(psi0)
      psi_nonli_ev = ev_from_erg(psi_nonli)
      psi = u0_ev + psi0_ev + psi_nonli_ev
      puts "u0_ev, psi0_ev, psi_nonli_ev, psi, large_a"
      check(u0_ev, psi0_ev, psi_nonli_ev, psi, large_a)
      puts "\n"
    end
  end

  def calc_moment_vasp
    p "calc_moment_vasp"
    #それぞれの定数を再定義，もともとはergだがJにするために．
    @boltz = 1.38064852e-23
    @planck = 6.626e-34
    @hbar = @planck/(2.0*PI)
    @avogadro = 6.023e23

    a0 = 2.5713342796904244e-10#vasp Cu
    p "a0"
    check(a0)
    for i in 1..10
      temp = 100*(i-1)
      temp = 10 if i==1
      theta = @boltz*temp
      gap = Array.new(6){ Array.new(201) }
      aa1,aa2=[],[]
      aa1[0]=a0
      comp_gap = 0
      for change in 1..5 #width to displace
        for same in 1..900
          gap[change][same] = (same-1)*(10**(-(11+change))).to_f#start 1e-12
          aa1[change] = aa1[change-1] + gap[change][same]
          aa2[change] = aa1[change]*sqrt(2)
          a1, a2 = aa1[change], aa2[change]
          a1 = a1*1e10#angstromに
          #p "k"
          k = (192.3313091-69.33663132*a1+135.4613676*(a1-2.5713342796904244)**2+88.20445888*(a1-2.5713342796904244)**3)
          k = k*1.60218e-19#evからJに変換
          atom_mass = 63.5460/@avogadro
          omega = sqrt(k/atom_mass)
          x = @hbar*omega/(2.0*theta)
          #p "gamma"
          gamma = (-1089.896157+529.2267532*a1)/6
          gamma = gamma*1.60218e-19#evからJに変換
          gt_k2 = gamma*theta/k**2
          ###################ここまでok
          #u0 = calc_u0(a1,a2)
          psi0 = calc_psi0(x,theta)
          large_a = calc_large_a(x, gt_k2)
          y0 = calc_y0(k, gamma, theta, large_a)
          psi_linear = calc_psi_linear(k, x, gamma, theta)#とりあえず変数はそのままnonli
          total_gap = comp_gap + gap[change][same]
          break if y0 - total_gap < 0
          a1 = a1*1e10
        end
        aa1[change] = aa1[change] - (10**(-(11+change))).to_f
        aa2[change] = aa1[change]*sqrt(2)
        comp_gap += gap[change][same-1]
        a1_cal = a0+comp_gap
        a2_cal = sqrt(2)*a1_cal
      end
      a1_cal = a1_cal*1.0e10
      a2_cal = a2_cal*1.0e10
      puts "T(K), a1, a2_cal(A), k, g="
      check(temp, a1_cal, a2_cal, k, gamma)
      puts "temp, u0, harmonic free, free="
      check(temp, "u0", "psi0", "psi_nonli")
      puts "\n"
    end

    p a1_cal

  end

  def calc_a0(m, n, r0)
    an=12.0 + 6.0/sqrt(2)**m
    am=12.0 + 6.0/sqrt(2)**n
    return exp( log(r0) + log(an/am)/(m-n) )
  end

  def calc_k(a1,a2,structure)
    case structure
    when "fitting_test"
      k = (5.692097474e6-2.141523420e14*a1+5.413711083e22*(a1-2.512603723e-8)**2-2.352757090e30*(a1-2.512603723e-8)**3)
    when "vasp"
      #0番目のポテンシャルのみだから1/2倍しないでいい．とりあえずevからergに変換しとく
      a1=a1*1e10
      k = (192.3313091-69.33663132*a1+135.4613676*(a1-2.5713342796904244)**2+88.20445888*(a1-2.5713342796904244)**3)*1.60218e-12
    else
    k = 2.0*@@potential.de2dr2(a1) + 4.0*@@potential.dedr(a1)/a1 + @@potential.de2dr2(a2) + 2.0*@@potential.dedr(a2)/a2;
    k =k*BOLTZ
    end
    return k
  end

  def calc_gamma(a1,a2,structure)
    case structure
    when "jindofcc"
      diff_u4_a1 = 2.0*@@potential.de4dr4(a1) + 12.0*@@potential.de3dr3(a1)/a1 - 42.0*@@potential.de2dr2(a1)/(a1*a1) + 42.0*@@potential.dedr(a1)/(a1*a1*a1)
      diff_u4_a2 = 2.0*@@potential.de4dr4(a2) + 12.0*@@potential.de2dr2(a2)/(a2*a2) - 12.0*@@potential.dedr(a2)/(a2*a2*a2)
      diff_x2y2_a1 = @@potential.de4dr4(a1) + 2.0*@@potential.de3dr3(a1)/a1 - 8.0*@@potential.de2dr2(a1)/(a1*a1) + 8.0*@@potential.dedr(a1)/(a1*a1*a1)
      diff_x2y2_a2 = 4.0*@@potential.de3dr3(a2)/a2 - 11.0*@@potential.de2dr2(a2)/(a2*a2) + 11.0*@@potential.dedr(a2)/(a2*a2*a2)
    when "sakakifcc"
      diff_u4_a1 = 2.0*@@potential.de4dr4(a1) + 12.0*@@potential.de3dr3(a1)/a1 - 6.0*@@potential.de2dr2(a1)/(a1*a1) + 6.0*@@potential.dedr(a1)/(a1*a1*a1)
      diff_u4_a2 = 2.0*@@potential.de4dr4(a2) + 12.0*@@potential.de2dr2(a2)/(a2*a2) - 12.0*@@potential.dedr(a2)/(a2*a2*a2)
      diff_x2y2_a1 = @@potential.de4dr4(a1) + 2.0*@@potential.de3dr3(a1)/a1 + 3.0*@@potential.de2dr2(a1)/(a1*a1) - 3.0*@@potential.dedr(a1)/(a1*a1*a1)
      diff_x2y2_a2 = 4.0*@@potential.de3dr3(a2)/a2 - 6.0*@@potential.de2dr2(a2)/(a2*a2) + 6.0*@@potential.dedr(a2)/(a2*a2*a2)
    when "fitting_test"
      return (4.629669951e23-1.411654254e31*a1)/6.0
    when "vasp"
      a1=a1*1e10
      return (-1089.896157+529.2267532*a1)/6.0*1.60218e-12
    else
      p "*******************missed structure**********************"
      exit(0)
    end
    @gamma1 = (1.0/48)*(diff_u4_a1 + diff_u4_a2)*BOLTZ
    @gamma2 = (6.0/48)*(diff_x2y2_a1 + diff_x2y2_a2)*BOLTZ
    return gamma = 4.0*(@gamma1+@gamma2)#p516(18)
    #return 4*@gamma1
  end

  def calc_u0(a1,a2)
    #最近接１２，第２近接６はfccだから．．．POTCARにその情報を入れたほうがいいかも
    z1=12
    z2=6
    return u0 = 0.5*(z1*@@potential.energy(a1) + z2*@@potential.energy(a2))
  end

  def calc_psi0(x, theta)
    arg0 = 1.0 - exp(-2.0*x)
    return psai0 = 3.0*theta*(x+log(arg0))
  end

  def calc_large_a(x, gt_k2)
    xcothx = x/tanh(x)
    xcothx = x*(exp(x)+exp(-x))/(exp(x)-exp(-x))
    xcothx2 = xcothx**2
    xcothx3 = xcothx**3
    xcothx4 = xcothx**4
    xcothx5 = xcothx**5
    small_a1 = 1.0 + xcothx/2.0
    small_a2 = 13.0/3.0 + 47.0*xcothx/6.0 + 23.0*xcothx2/6.0 + xcothx3/2.0
    small_a3 = -(25.0/3.0 + 121.0*xcothx/6.0 + 50.0*xcothx2/3.0 + 16.0*xcothx3/3.0 + xcothx4/2.0)
    small_a4 = 43.0/3.0 + 93.0*xcothx/2.0 + 169.0*xcothx2/3.0 + 83.0*xcothx3/3.0 + 22.0*xcothx4/3.0 + xcothx5/2.0
    return large_a = small_a1 + small_a2*gt_k2**2 + small_a3*gt_k2**3 + small_a4*gt_k2**4
  end

  def calc_y0(k, gamma, theta, large_a)
    return sqrt( 2.0*gamma*theta**2*large_a/(3.0*k**3))
  end

  def calc_psi_nonli(k, x, gamma,theta)
    fac_gam = @gamma1*@gamma1 + 2.0*@gamma1*@gamma2
    @gamma1
    @gamma2
    xcothx = x/tanh(x)
    #xcothx = x*(exp(x)+exp(-x))/(exp(x)-exp(-x))
    xcothx2 = xcothx**2
    fac1 = 1.0 + xcothx/2.0
    fac2 = 1.0 + xcothx
    term1 = (theta/k)**2 * (@gamma2*xcothx2 - 2.0*@gamma1*fac1/3.0)
    term2 = (4.0*@gamma2**2 / 3.0)*xcothx*fac1-2.0*fac_gam*fac1*fac2
    term2 = (2.0*theta**3/k**4)*term2
    return 3.0*(term1+term2)
  end
  def calc_psi_linear(k, x, gamma,theta)#追加で実装　論文(1988)p515
    xcothx = x/tanh(x)
    #xcothx = x*(exp(x)+exp(-x))/(exp(x)-exp(-x))
    fac1 = 1.0 + xcothx/2.0
    fac2 = 1.0 + xcothx
    term1 = ((gamma*theta**2)/(6.0*k**2))*fac1
    term2 = ((gamma**2*theta**3)/(4.0*k**4))*fac1*fac2
    return -term1-term2
  end

  def ev_from_kelvin(a)
    return a*8.617385e-05
  end

  def ev_from_erg(a)
    return a*6.2415064e11
  end
end

class MomentPlot < CalcMoment
  def initialize(structure, range)
    puts "Hi,MomentPlot"
    @structure=structure
    puts "structure = #{@structure}"
    select
    #2.5e-8から2.6e-8
    min = 2.6e-8 - (1e-9)*range
    max = 2.5e-8 + (1e-9)*range
    stepno = (max - min)/50
    moment_plot( min, max, stepno)
  end

  def moment_plot(min, max, stepno)
    puts "start moment_check"
    make_gp
      #make k date ．ベタ書き
    tmp=""
    min.step(max, stepno) do |r|
      tmp << r.to_s
      tmp << " "
      tmp << calc_k(r, r*sqrt(2), @structure).to_s
      tmp << " "
      tmp << (2.0*@@potential.de2dr2(r)*BOLTZ).to_s
      tmp << " "
      tmp << (4.0*@@potential.dedr(r)/r*BOLTZ).to_s
      tmp << " "
      tmp << (@@potential.de2dr2(r*sqrt(2))*BOLTZ).to_s
      tmp << " "
      tmp << (2.0*@@potential.dedr(r*sqrt(2))/(r*sqrt(2))*BOLTZ).to_s
      tmp << "\n"
    #  printf("%e %.10e\n", r, calc_k(r, r*sqrt(2)))
      File.open("./k.txt",'w'){|io| io.print(tmp)}
    end
    #make gamma data
    tmp=""
    min.step(max, stepno) do |r|
      tmp << r.to_s
      tmp << " "
      tmp << calc_gamma(r, r*sqrt(2), @structure).to_s
      tmp << " "
      tmp << (@gamma1*4).to_s
      tmp << " "
      tmp << (@gamma2*4).to_s
      tmp << "\n"
      #printf("%e %.10e\n", r, calc_k(r, r*sqrt(2)))
      File.open("./gamma.txt",'w'){|io| io.print(tmp)}
    end
    system "gnuplot plot_k.gp"
    system "gnuplot plot_gamma.gp"
  end
  def make_gp #gnuplot用，ベタ書き，それぞれのポテンシャルのフォルダの中に保存される．
    File.open("./plot_k.gp",'w'){|io|
      io.print("
      set format x \"%.2e\"\n
      set xlabel \"a1\"
      set ylabel \"erg\"
      plot \"./k.txt\" using 1:2 title \"k\" w lp\n
      replot \"./k.txt\" using 1:3 title \"term 1,    2*de2dr2(a1)\" w lp\n
      replot \"./k.txt\" using 1:4 title \"term 2,    4*dedr(a1)/a1\" w lp\n
      replot \"./k.txt\" using 1:5 title \"term 3,    de2dr2(a2)\" w lp\n
      replot \"./k.txt\" using 1:6 title \"term 4,    2*dedr(a2)/a2\" w lp\n
      ")}

    File.open("./plot_gamma.gp",'w'){|io|
      io.print("
      set format x \"%.2e\"\n
      set xlabel \"a1\"
      set ylabel \"erg\"
      plot \"./gamma.txt\" using 1:2 title \"gamma\" w lp\n
      replot \"./gamma.txt\" using 1:3 title \"xxxx\" w lp\n
      replot \"./gamma.txt\" using 1:4 title \"xxyy\" w lp\n
      ")}
  end
end
#lj.rbにclassLJ_jindoがあるけど，後の切り分けのために新たに作成.理想はPOTCARの中
class DiffLjJindo
  attr_reader :m, :n, :r0, :atom_mass
  def initialize(src)
    @@d0,@@m,@@n,@@r0=src[:d0],src[:m],src[:n],src[:r0_no_angstrom]
    @m, @n, @r0, @atom_mass=src[:m],src[:n], src[:r0_no_angstrom], src[:atom_mass_no_avogado]#計算のために参照できるようにしとく
  end

  def energy(r)
    ene=-@@d0*((r/@@r0)**(-@@n)*@@m-(r/@@r0)**(-@@m)*@@n)/(@@m-@@n)
  end

  def dedr(r)
    dedr=-@@d0*@@m*@@n / ((@@n-@@m)*r) * ((@@r0/r)**@@n - (@@r0/r)**@@m)
  end

  def de2dr2(r)
    de2dr2 = @@d0*@@m*@@n / ((@@n-@@m)*r**2) * ((@@n+1)*(@@r0/r)**@@n - (@@m+1)*(@@r0/r)**@@m)
  end

  def de3dr3(r)
    de3dr3 = @@d0*@@m*@@n / ((@@n-@@m)*r**3) * ((@@m+1)*(@@m+2)*(@@r0/r)**@@m - (@@n+1)*(@@n+2)*(@@r0/r)**@@n)
  end

  def de4dr4(r)
    de4dr4 = @@d0*@@m*@@n / ((@@n-@@m)*r**4) * ((@@n+1)*(@@n+2)*(@@n+3)*(@@r0/r)**@@n - (@@m+1)*(@@m+2)*(@@m+3)*(@@r0/r)**@@m)
  end
end
