      implicit real*8 (a-h,o-z)
      common/blk1/ vm,vn,d0,r0,a1,a2
      dimension pot_der1(3),pot_der2(3),pot_der3(3),pot_der4(3)
     * ,der_u4(3),der_x2y2(3),diff(1000),y0_tenta(5,1000)
     * ,aa1(0:20),aa2(0:20)
c****** Calculation of Thermal Expansion Coefficients of
c     alat, rnn(nearest neighbour distance), and rnnn
c     will be determined self-consistently, should not given
c     as input parameters
c***** Atomic Mass from Rika Nenpiyou (1991)
c***** Lattice Constant from X-Ray Diffraction B.D. Cullity
c**   /2000.7.29/
      open(10,file='vhung.fcc.data',status='unknown')
c
c      write(6,*) '/home/jindo/vietnam/vhung.fcc.f'
c      write(10,*) '/home/jindo/vietnam/vhung.fcc.f'
c
      write(6,*) 'Cal. of Temp Dependence of Lattice Const'
      write(10,*) 'Cal. of Temp Dependence of Lattice Const'
      pi=2.d0*dasin(1.d0)
      s3=dsqrt(3.d0)
      s2=dsqrt(2.d0)
      boltz_con=1.380658d-16
      planck=6.626d-27
      hbar=planck/(2.d0*pi)
      avogado=6.023d+23
	mat=2
      go to (1,2,3), mat
    1 continue
      write(6,*) 'Calculation of Ag'
      write(10,*) 'Calculation of Ag'
      d0=3325.6
      vn=11.5d0
      vm=5.5d0
      r0=2.8760d-08
      alat=4.0856d-08
      atom_mass=107.8682d0/avogado
      go to 10
    2 continue
      write(6,*) 'Calculation of Cu'
      write(10,*) 'Calculation of Cu'
      d0=4125.7d0
      vn=9.0d0
      vm=5.5d0
      r0=2.5487d-08
      alat=3.6153d-08
      atom_mass=63.546d0/avogado
      go to 10
    3 continue
      write(6,*) 'Calculation of Au'
      write(10,*) 'Calculation of Au'
      d0=4683.d0
      vn=10.5d0
      vm=5.5d0
      r0=2.8751d-08
      alat=4.0783d-08
      atom_mass=196.96654d0/avogado
   10 continue
c   
      r0_zero=r0
      a1=s2*alat/2.d0
      a2=alat
c     
      fac_fcc2=dsqrt(2.d0)
      fac_fcc3=dsqrt(3.d0)
c      
c      an=12.d0+6.d0/fac_fcc2**vm+24.d0/fac_fcc3**vm
c      am=12.d0+6.d0/fac_fcc2**vn+24.d0/fac_fcc3**vn
c
      an=12.d0+6.d0/fac_fcc2**vm
      am=12.d0+6.d0/fac_fcc2**vn
c
      dammi=dlog(r0)+dlog(an/am)/(vm-vn)
      a0=dexp(dammi)
c    
      write(6,315) a0,r0,an,am,d0
      write(10,315) a0,r0,an,am,d0
  315 format(1h ,'a0,r0,an,am,d0=',5d14.6)
c
c***** a0=nearest-neighbour distance at 0K
c
      write(6,*) 'below a1 and a2 from experiments'
      write(10,*) 'below a1 and a2 from experiments'
      write(6,300) vm,vn,d0,a1,a2
      write(10,300) vm,vn,d0,a1,a2
  300 format(1h ,'vm,vn,d0,a1,a2=',2f7.2,f9.3,3d14.7)
c
      itest=2
      if(itest.ge.2) go to 100 
c***** Calculation of potential energy versus radius
      imax=51
      do 20 i=1,imax
      ri=0.02d0*dfloat(i)+0.84d0
      r=ri*a1
      rh=r0/r
      pot=d0/(vn-vm)*(vm*rh**vn-vn*rh**vm) 
      write(6,301) i,ri,r,pot
  301 format(1h ,'i,ri,r,pot=',i3,f8.4,f9.5,f10.5)
   20 continue
  100 continue
ccccccccccccccccccccccccccccccccccccccccccccc
cc    Calculation of k=mw**2   and r(gamma)=Vu van Hung
      write(6,*) 'Calcul. of k=mw**2 and r(gamma)'
c 
      itmax=10
c      itmax=1
      do 120 it=1,itmax
      temp=100.d0*dfloat(it-1)
      if(it.eq.1) temp=10.d0
      teta=boltz_con*temp
c      write(6,*) 'temp(K)=', temp
	y0_cal=0.d0
       aa1(0)=r0
       aa1(0)=a0
      do 200 kaisu=1,5
      do 210 jy0=1,200
      if(kaisu.eq.1) y0_tenta(kaisu,jy0)=dfloat(jy0-1)*1.0d-10
      if(kaisu.eq.2) y0_tenta(kaisu,jy0)=dfloat(jy0-1)*1.0d-11
      if(kaisu.eq.3) y0_tenta(kaisu,jy0)=dfloat(jy0-1)*1.0d-12
      if(kaisu.eq.4) y0_tenta(kaisu,jy0)=dfloat(jy0-1)*1.0d-13
      if(kaisu.eq.5) y0_tenta(kaisu,jy0)=dfloat(jy0-1)*1.0d-14
ccccccccccccccccccccccccccccccccccccccccccccccccc
      aa1(kaisu)=aa1(kaisu-1)+y0_tenta(kaisu,jy0)
      aa2(kaisu)=aa1(kaisu)*s2
      a1=aa1(kaisu)
      a2=aa2(kaisu)
ccccccccccccccccccccccccccccccccccccccccccccccccc
      do 101 k=1,2
      if(k.eq.1) rr=a1
      if(k.eq.2) rr=a2
      call deriv(rr,pt_der1,pt_der2,pt_der3,pt_der4)
      pot_der1(k)=pt_der1*1.380658d-16
      pot_der2(k)=pt_der2*1.380658d-16
      pot_der3(k)=pt_der3*1.380658d-16
      pot_der4(k)=pt_der4*1.380658d-16
  101 continue
c**  
      vkfcc=2.d0*pot_der2(1)  + 4.d0*pot_der1(1)/a1 
     *          +pot_der2(2)  +2.d0*pot_der1(2)/a2
      omega=dsqrt(vkfcc/atom_mass)
c      write(6,*) 'test=', dsqrt(vkfcc/atom_mass)!平方根の中身チェック
      x=hbar*omega/(2.d0*teta)
c      write(6,*) 'vkfcc=', vkfcc
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      der_u4(1)=2.d0*pot_der4(1)+12.d0*pot_der3(1)/a1
     *-42.d0*pot_der2(1)/(a1*a1)+42.d0*pot_der1(1)/(a1*a1*a1)
      der_u4(2)=2.d0*pot_der4(2)+12.d0*pot_der2(2)/(a2*a2)
     *-12.d0*pot_der1(2)/(a2*a2*a2)
c*  
      der_x2y2(1)=pot_der4(1)+2.d0*pot_der3(1)/a1
     *-8.d0*pot_der2(1)/(a1*a1)+8.d0*pot_der1(1)/(a1*a1*a1)
      der_x2y2(2)=4.d0*pot_der3(2)/a2-11.d0*pot_der2(2)/(a2*a2)
     *+11.d0*pot_der1(2)/(a2*a2*a2)
c*
      gamma1=(1.d0/48.d0)*(der_u4(1) + der_u4(2) )
      gamma2=(6.d0/48.d0)*(der_x2y2(1)+der_x2y2(2) )
      gamma=4.d0*( gamma1+gamma2 )
      gtbyk2=gamma*teta/vkfcc**2
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccc Calculation of Free-Energy of Crystal ccccccccccccccccccc
      r0a1=r0/a1
      r0a2=r0/a2
ccccccccccccccc Below z1 and z2 added on 99.3.16.
      z1=12.d0
      z2=6.d0
      pote_a1=d0*( vm*r0a1**vn - vn*r0a1**vm )/(vn-vm)
      pote_a2=d0*( vm*r0a2**vn - vn*r0a2**vm )/(vn-vm)
      U0=0.5d0*( z1*pote_a1+z2*pote_a2 )
      arg0=1.d0-dexp(-2.d0*x)
      psai0=3.d0*teta*( x + dlog(arg0) )
      xcthx=x*( dexp(x)+dexp(-x) )/( dexp(x)-dexp(-x) )
      xcthx2=xcthx**2
      xcthx3=xcthx**3
      xcthx4=xcthx**4
      xcthx5=xcthx**5
      fac1=1.d0+xcthx/2.d0
      fac2=1.d0+xcthx
      small_a1=fac1
      small_a2=13.d0/3.d0+47.d0*xcthx/6.d0+23.d0*xcthx2/6.d0 
     *        +xcthx3/2.d0
      small_a3=25.d0/3.d0+121.d0*xcthx/6.d0+50.d0*xcthx2/3.d0
     *        +16.d0*xcthx3/3.d0+xcthx4/2.d0
      small_a3=-small_a3
      small_a4=43.d0/3.d0+93.d0*xcthx/2.d0+169.d0*xcthx2/3.d0
     *        +83.d0*xcthx3/3.d0+22.d0*xcthx4/3.d0+xcthx5/2.d0
      alarge=small_a1+small_a2*gtbyk2**2+small_a3*gtbyk2**3
     *      +small_a4*gtbyk2**4
c      y0=2.d0*gamma*teta*teta*alarge/(3.d0*vkbcc**3)
c***** below corrected on 99.3.18. **********
      y0=2.d0*gamma*teta*teta*alarge/(3.d0*vkfcc**3)
      y0=dsqrt(y0)
      fac_gam=gamma1*gamma1+2.d0*gamma1*gamma2
c      term1=(teta/vkbcc)**2*(gamma2*xcthx2-2.d0*gamma1*fac1/3.d0)
c     99.3.18 below two vkbcc corrected
      term1=(teta/vkfcc)**2*(gamma2*xcthx2-2.d0*gamma1*fac1/3.d0)
      term2=(4.d0*gamma2**2/3.d0)*xcthx*fac1-2.d0*fac_gam*fac1*fac2
      term2=(2.d0*teta**3/vkfcc**4)*term2
c      term2=(2.d0*teta**3/vkbcc**4)*term2
      psai_nonli=3.d0*( term1+term2 )
      psai=U0 + psai0 + psai_nonli
c    
c      diff(jy0)=y0-y0_tenta(kaisu,jy0)
      y0_comp=y0_cal + y0_tenta(kaisu,jy0)
ccccccccccc
      diff(jy0)=y0-y0_comp
      jy0m1=jy0-1
      if(jy0m1.lt.1) jy0m1=1
      diff_sign=dsign(1.d0,diff(jy0))-dsign(1.d0,diff(jy0m1))
      diff_sign=dabs( diff_sign )
cccccc above gamma1 and gamma2 for bcc structure ******************
c      write(6,302) rr,pot_der1(1),pot_der2(1),pot_der1(2),pot_der2(2)
c      write(6,303) rr,pot_der3(1),pot_der4(1),pot_der3(2),pot_der4(2)
c      write(6,307) rr,alarge,y0,y0_tenta(kaisu,jy0),diff
c      write(6,308) a1,a2,y0,y0_tenta(kaisu,jy0),diff(jy0),diff_sign
      if(diff_sign.gt.1.5d0) go to 205
  210 continue
  205 continue
c      write(6,308) a1,a2,y0,y0_tenta(kaisu,jy0),diff(jy0),diff_sign
      aa1(kaisu)=aa1(kaisu)+y0_tenta(kaisu,jy0-1)-y0_tenta(kaisu,jy0)
      aa2(kaisu)=aa1(kaisu)*s2
      y0_cal=y0_cal + y0_tenta(kaisu,jy0-1)
      a1_cal=a0+y0_cal
      a2_cal=dsqrt(2.d0)*a1_cal
c      write(6,309)  kaisu,jy0,y0_cal,a1_cal
  200 continue
      a1_cal=a1_cal*1.0d+08
      a2_cal=a2_cal*1.0d+08
      hamonic=U0 + psai0
      write(6,310) temp,a1_cal,a2_cal,vkfcc,gamma
c      write(10,310) temp,a1_cal,a2_cal,vkfcc,gamma
      write(6,312) temp,U0,psai0,psai_nonli
  120 continue
c
  302 format(1h ,'rr,pot_der1,2=',d10.5,4d15.7)
  303 format(1h ,'rr,pot_der3,4=',d10.5,4d15.7)
  304 format(1h ,'rr,vkfcc,der_u4,der_x2y2(1)=',d10.5,5d15.7)
  305 format(1h ,'rr,vkfcc,der_u4,der_x2y2(2)=',d10.5,5d15.7)
  306 format(1h ,'rr,vkfcc,gamma1,gamma2=',d10.5,5d15.7)
  307 format(1h ,'rr,alarge,y0,y0_t,dif=',d10.5,7d15.7)
  308 format(1h ,'a1,a2,y0,y0_t,dif,dif_sign=',10d15.7)
  309 format(1h ,'kaisu,jy0,y0_cal,a1_cal=',2i3,10d15.7)
  310 format(1h ,'     T(K),a1,a2_cal(A),k,g=',f10.2,2f12.5,2d15.5)
  312 format(1h ,'temp,U0,harmonic free,free=',f10.2,5d15.5)
      stop
      end
      subroutine deriv(rr,pot_der1,pot_der2,pot_der3,pot_der4)
      implicit real*8 (a-h,o-z)
      common/blk1/ vm,vn,d0,r0,a1,a2
      rh=r0/rr
      r2=rr*rr
      r3=rr**3
      r4=rr**4
      vnp1=vn+1
      vmp1=vm+1
      vnp2=vn+2
      vmp2=vm+2
      vnp3=vn+3
      vmp3=vm+3
      pot_der1=-d0*vm*vn/((vn-vm)*rr)*(rh**vn-rh**vm)
      pot_der2=d0*vm*vn/((vn-vm)*r2)*(vnp1*rh**vn-vmp1*rh**vm)
      pot_der3=d0*vm*vn/((vn-vm)*r3)*(vmp1*vmp2*rh**vm
     *                               -vnp1*vnp2*rh**vn)
      pot_der4=d0*vm*vn/((vn-vm)*r4)*(vnp1*vnp2*vnp3*rh**vn
     *                               -vmp1*vmp2*vmp3*rh**vm)
      return
      end
