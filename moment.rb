#Cu,Auで値一致，Agは負の平方根取るからだめ，fortranの方もNaNになってる
include Math

def sqrt(n)
  Math::sqrt(n)
end

def log(n)
  Math::log(n)
end

def exp(n)
  Math::exp(n)
end

def check(*val)
  p val.join(", ")
end

def diff_sign?(n1, n2)
  return true if n1*n2 < 0
  return false
end

def deriv(rr) #動作確認済み
  rh = R0/rr
  r2 = rr*rr
  r3 = rr**3
  r4 = rr**4
  vnp1 = VN+1
  vnp2 = VN+2
  vnp3 = VN+3
  vmp1 = VM+1
  vmp2 = VM+2
  vmp3 = VM+3
  $p_der[1] = -D0*VM*VN / ((VN-VM)*rr) * (rh**VN - rh**VM) #potential 1次微分
  $p_der[2] = D0*VM*VN / ((VN-VM)*r2) * (vnp1*rh**VN - vmp1*rh**VM)
  $p_der[3] = D0*VM*VN / ((VN-VM)*r3) * (vmp1*vmp2*rh**VM - vnp1*vnp2*rh**VN)
  $p_der[4] = D0*VM*VN / ((VN-VM)*r4) * (vnp1*vnp2*vnp3*rh**VN - vmp1*vmp2*vmp3*rh**VM)
end

$p_der=[]
der_u4 = []
der_x2y2 = []
diff = []

PI = Math::PI
BOLTZ = 1.380658e-16
PLANCK = 6.626e-27
HBAR = PLANCK/(2.0*PI)
AVOGADO = 6.023e23

material = "Cu"
#原子固有の定数
case material
when "Cu" then
  puts "Calculation of Cu"
  D0 = 4125.70
  VN = 9.00
  VM = 5.50
  R0 = 2.5487e-8 #0Kにおける最近接原子間距離
  ALAT = 3.6153e-8 #使われてない
  ATOM_MASS = 63.5460/AVOGADO #Cuの原子量が63.546
when "Ag" then
  puts "Calculation of Ag"
  D0 = 3325.6
  VN = 11.5
  VM = 5.5
  R0 = 2.876e-8
  ALAT = 4.0856e-8
  ATOM_MASS = 107.8682/AVOGADO 
when "Au" then
  puts "Calculation of Au"
  D0 = 4683.0
  VN = 10.5
  VM = 5.5
  R0 = 2.8751e-8
  ALAT = 4.0783e-8
  ATOM_MASS = 196.96654/AVOGADO 
end

a1 = sqrt(2)*ALAT/2.0#使われてない
a2 = ALAT#使われてない
fac_fcc2 = sqrt(2)
fac_fcc3 = sqrt(3)
an = 12.0 + 6.0/fac_fcc2**VM
am = 12.0 + 6.0/fac_fcc2**VN
a0 = exp( log(R0) + log(an/am)/(VM-VN) )#今後の計算に関係してくる．．最後にa0+y0_calしてるから，最近接原子の初期値っぽい．．

p "a0, R0, an, am, D0"
check(a0, R0, an, am, D0)

p "VM, VN, D0, a1, a2"
check(VM, VN, D0, a1, a2)

for it in 1..10
  temp = 100*(it-1);
  temp = 10 if it==1
  theta = BOLTZ*temp
  p "temp(K)=#{temp}"

  y0_cal = 0.0
  aa1 = []
  aa2 = []
  aa1[0] = R0 #今回は初期値a0
  aa1[0] = a0
  y0_tenta = Array.new(6){ Array.new(1000) } #要素数とりあえず合わしとく
  pot_der = Array.new(10){ Array.new(10) } #とりあえず多めに宣言
  #ループの始まり
  for kaisu in 1..5 do
    for jy0 in 1..200 do
      y0_tenta[kaisu][jy0] = (jy0-1)*(10**(-(9+kaisu))).to_f #0,1e-10, 2e-10, 3e-10.....0, 1e-11, 2e-11,....
#      y0_tenta[kaisu][jy0] = (jy0-1)*(10**(-(11+kaisu))).to_f #0,1e-10, 2e-10, 3e-10.....0, 1e-11, 2e-11,....
      #原子の距離に対応するaa1,aa2を少しずつ増やして回してる
      #距離増やして力の釣り合う距離を探してる？，０を越えたら符号が変わるからそこでループ終了して，少し小さい値でもう一度を繰り返してる
      aa1[kaisu] = aa1[kaisu-1] + y0_tenta[kaisu][jy0]
      aa2[kaisu] = aa1[kaisu]*sqrt(2)
      a1 = aa1[kaisu]
      a2 = aa2[kaisu]

      deriv(a1)
      for i in 1..4 do
        pot_der[i][1] = BOLTZ*$p_der[i]
      end

      deriv(a2)
      for i in 1..4 do
        pot_der[i][2] = BOLTZ*$p_der[i]
      end

      vkfcc = 2.0*pot_der[2][1] + 4.0*pot_der[1][1]/a1 + pot_der[2][2] + 2.0*pot_der[1][2]/a2;
      #            p vkfcc/ATOM_MASS#平方根正負チェック
      omega = sqrt(vkfcc/ATOM_MASS)
      x = HBAR*omega/(2.0*theta)
      
      structure = "Jindo_fcc"
      case structure
      when "Jindo_fcc"
      der_u4[1] = 2.0*pot_der[4][1] + 12.0*pot_der[3][1]/a1 - 42.0*pot_der[2][1]/(a1*a1) + 42.0*pot_der[1][1]/(a1*a1*a1)
      der_u4[2] = 2.0*pot_der[4][2] + 12.0*pot_der[2][2]/(a2*a2) - 12.0*pot_der[1][2]/(a2*a2*a2)
      der_x2y2[1] = pot_der[4][1] + 2.0*pot_der[3][1]/a1 - 8.0*pot_der[2][1]/(a1*a1) + 8.0*pot_der[1][1]/(a1*a1*a1)
      der_x2y2[2] = 4.0*pot_der[3][2]/a2 - 11.0*pot_der[2][2]/(a2*a2) + 11.0*pot_der[1][2]/(a2*a2*a2)
      when "Sakaki_fcc"
      der_u4[1] = 2.0*pot_der[4][1] + 12.0*pot_der[3][1]/a1 - 6.0*pot_der[2][1]/(a1*a1) + 6.0*pot_der[1][1]/(a1*a1*a1)
      der_u4[2] = 2.0*pot_der[4][2] + 12.0*pot_der[2][2]/(a2*a2) - 12.0*pot_der[1][2]/(a2*a2*a2)
      der_x2y2[1] = pot_der[4][1] + 2.0*pot_der[3][1]/a1 + 3.0*pot_der[2][1]/(a1*a1) - 3.0*pot_der[1][1]/(a1*a1*a1)
      der_x2y2[2] = 4.0*pot_der[3][2]/a2 - 6.0*pot_der[2][2]/(a2*a2) + 6.0*pot_der[1][2]/(a2*a2*a2)
      end
      gamma1 = (1.0/48)*(der_u4[1] + der_u4[2])#p516(18)
      gamma2 = (6.0/48)*(der_x2y2[1] + der_x2y2[2])
      gamma = 4.0*(gamma1+gamma2)
      gtbyk2 = gamma*theta/vkfcc**2

      #Calculation of Free-Energy of Crystal
      r0a1 = R0/a1
      r0a2 = R0/a2
      z1 = 12.0
      z2 = 6.0

      pote_a1 = D0*(VM*r0a1**VN - VN*r0a1**VM)/(VN-VM)#ポテンシャルの計算，a1が距離に対応
      pote_a2 = D0*(VM*r0a2**VN - VN*r0a2**VM)/(VN-VM)#a2が距離に対応
      u0 = 0.5*(z1*pote_a1 + z2*pote_a2) #ポテンシャル
      arg0 = 1.0 - exp(-2.0*x)
      psai0 = 3.0*theta*(x+log(arg0))#thesis p515,fcc is p516

      xcthx = x*(exp(x)+exp(-x))/(exp(x)-exp(-x)) #x cothx
      xcthx2 = xcthx**2
      xcthx3 = xcthx**3
      xcthx4 = xcthx**4
      xcthx5 = xcthx**5

      fac1 = 1.0 + xcthx/2.0
      fac2 = 1.0 + xcthx
      small_a1 = fac1
      small_a2 = 13.0/3.0 + 47.0*xcthx/6.0 + 23.0*xcthx2/6.0 + xcthx3/2.0
      small_a3 = 25.0/3.0 + 121.0*xcthx/6.0 + 50.0*xcthx2/3.0 + 16.0*xcthx3/3.0 + xcthx4/2.0
      small_a3 = -small_a3
      small_a4 = 43.0/3.0 + 93.0*xcthx/2.0 + 169.0*xcthx2/3.0 + 83.0*xcthx3/3.0 + 22.0*xcthx4/3.0 + xcthx5/2.0
      alarge = small_a1 + small_a2*gtbyk2**2 + small_a3*gtbyk2**3 + small_a4*gtbyk2**4

      y0 = 2.0*gamma*theta**2*alarge/(3.0*vkfcc**3) #外力がない場合の最近接原子間距離?熱膨張のみ
      y0 = sqrt(y0)
      fac_gam = gamma1*gamma1 + 2.0*gamma1*gamma2

      term1 = (theta/vkfcc)**2 * (gamma2*xcthx2 - 2.0*gamma1*fac1/3.0)
      term2 = (4.0*gamma2**2 / 3.0)*xcthx*fac1-2.0*fac_gam*fac1*fac2
      term2 = (2.0*theta**3/vkfcc**4)*term2

      psai_nonli = 3.0*(term1+term2) #nonlinear? thesis(1988)p516(19)
      psai = u0 + psai0 + psai_nonli #使ってない
      #      "psai=#{psai}"
      y0_comp = y0_cal + y0_tenta[kaisu][jy0]
      diff[jy0] = y0 - y0_comp #0以下とかの方がよさそう．．
      #      p "jy0=#{jy0} ##############################################################################"
      break if diff_sign?(diff[jy0], diff[jy0-1]) if jy0 > 1
    end
    #ここまでjy0
    aa1[kaisu] = aa1[kaisu] + y0_tenta[kaisu][jy0-1] - y0_tenta[kaisu][jy0] if jy0 >1 #一応jy0 > 1つけとく
    aa2[kaisu] = aa1[kaisu]*sqrt(2)

    y0_cal = y0_cal + y0_tenta[kaisu][jy0-1]
    a1_cal = a0+y0_cal
    a2_cal = sqrt(2)*a1_cal
    # p "kaisu=#{kaisu} ##############################################################################"
    #ここまでkaisu
  end
  a1_cal = a1_cal*1.0e8
  a2_cal = a2_cal*1.0e8
  harmonic = u0+psai0#これ使ってない
  puts "T(K), a1, a2_cal(A), k, g="
  check(temp, a1_cal, a2_cal, vkfcc, gamma)
  puts "temp, u0, harmonic free, free="
  check(temp, u0, psai0, psai_nonli)
  #eVに単位変換
  u0_ev = u0*8.617385e-05
  psai0_ev = psai0*6.2415064e11
  psai_nonli_ev = psai_nonli*6.2415064e11
  psai = u0_ev + psai0_ev + psai_nonli_ev
  puts "u0_ev, psai0_ev, psai_nonli_ev, psai"
  check(u0_ev, psai0_ev, psai_nonli_ev, psai)
  puts"\n"
end
