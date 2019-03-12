module chiral_constants 
  use constants, only : pi  

  ! chiral order definition
  INTEGER, PARAMETER :: LO   = 0
  ! all contributions vanish at order 1
  ! due to parity and time-reversal  invariance
  INTEGER, PARAMETER :: NLO  = 2
  INTEGER, PARAMETER :: NNLO = 3
  INTEGER, PARAMETER :: N3LO = 4
  
  real*8, parameter :: c1 =  -0.9186395287347203d0
  real*8, parameter :: c2 =   0.d0 
  real*8, parameter :: c3 =  -0.3888687492763241d0
  real*8, parameter :: c4 =   0.4310327160829740d0
  real*8, parameter :: gA = 1.29d0 
  real*8, parameter :: gA2 = gA*gA 
  real*8, parameter :: gA4 = gA2*gA2 
  real*8, parameter :: fpi = 92.4d0 
  real*8, parameter :: fpi2 = fpi*fpi
  real*8, parameter :: fpi4 = fpi2*fpi2
  real*8, parameter :: twopi = pi*2.d0
  real*8, parameter :: pi2 = pi*pi
  real*8 :: const(21) 
  real*8, parameter :: proton_mass = 938.272d0 
  real*8, parameter :: neutron_mass = 939.5653d0
  real*8 :: mnuc(-1:1), mpi(-1:2)
  REAL*8 :: mnuc2(-1:1)   ! nucleon mass squared
  REAL*8 :: mpi2(-1:2)    ! pion mass squared
  REAL*8 :: mpi4(-1:2)    ! mpi2 squared
  REAL*8 :: mpi5(-1:2)    ! mpi^5
  REAL*8 :: twompi(-1:2)  ! two times pion mass
  REAL*8 :: fourmpi2(-1:2)! four times pion mass squared
  COMPLEX*16 :: sigma_x(2,2), sigma_y(2,2), sigma_z(2,2)
  REAL*8, PARAMETER :: lambda = 500.D0 
  REAL*8, PARAMETER :: sfr = 700.d0 
  REAL*8, PARAMETER :: sfr2 = sfr*sfr
  REAL*8 :: sfr_heavyside(-1:2) ! THETA(sfr-twompi)
  REAL*8 :: CS(-1:1), CT(-1:1), cnlo(1:7) 
  REAL*8 :: CE, CD, Econst, Dconst 
  
CONTAINS 
  
  subroutine init_chp_constants 
    implicit none   
    double precision :: c1s0(-1:1), c3s1(-1:1), cnlo_pw(1:7)
    integer :: i 
    
    mnuc(-1) = proton_mass
    mnuc(0)  = 0.5d0*(proton_mass+neutron_mass)
    mnuc(+1) = neutron_mass
  
    mpi(-1)  = 139.5702d0
    mpi(0)   = 134.9766d0
    mpi(+1)  = 139.5702d0
    mpi(+2)  = ( mpi(-1)+mpi(0)+mpi(+1) ) /3.d0  

    mnuc2(:)   = mnuc*mnuc 
    mpi2(:)    = mpi*mpi   
    mpi4(:)    = mpi2*mpi2
    mpi5(:)    = mpi4*mpi 
    twompi(:)  = 2.0D0*mpi 
    fourmpi2(:)= 4.0D0*mpi2
  
 
    const = 0.0D0
    !1: 1/(2pi)^3
    const(1) = 1.0D0/(twopi**3)
    !2: gA^2/(4*fpi^2)
    const(2) = gA2/(4.0D0*fpi2)
    !3: 384pi^2*fpi^4
    const(3) = 384.0D0*pi2*fpi4
    !4: 5gA^4-4gA^2-1
    const(4) = 5.0D0*gA4-4.0D0*gA2-1.0D0
    !5: 23gA^4-10gA^2-1
    const(5) = 23.0D0*gA4-10.0D0*gA2 - 1.0D0
    !6: 48gA^4
    const(6) = 48.0D0*gA4
    !7: 3gA^4
    const(7) = 3.0D0*gA4
    !8: 64pi^2fpi^4
    const(8) = 64.0D0*pi2*fpi4
    !9: 3gA^2/16pifpi^4
    const(9) = 3.0D0*gA2/(16.0D0*pi*fpi4)
    !10: ga^2/16
    const(10) = ga2/16.0D0
    !11: 2.0D0*(2c1-c3)
    const(11) = 2.0D0*(2.0D0*c1-c3)
    !12 : const(7)/256pifpi^4
    const(12) = const(7)/(256.0D0*pi*fpi4)
    !13: gA^2/(128pifpi^4)
    const(13) = gA2/(128.0D0*pi*fpi4)
    !14: gA^4/(128pifpi^4)
    const(14) = gA4/(128.0D0*pi*fpi4)
    !15: 3gA^4/(512pifpi^4)
    const(15) = 3.0D0*gA4/(512.0D0*pi*fpi4)
    !16: gA2/(32pifpi^4)
    const(16) = gA2/(32.0D0*pi*fpi4)
    !17: gA2/8
    const(17) = gA2/8.0D0
    !18: gA4/(256pifpi^4)
    const(18) = gA4/(256.0D0*pi*fpi4)
    !19: 3gA4/(32pifpi^4)
    const(19) = 3.0D0*gA4/(32.0D0*pi*fpi4)
    !20: const(16)*(1-gA2)
    const(20) = const(16)*(1.0D0-gA2)
    !21: gA/(8*fpi^2)
    const(21) = gA/(8.0D0*fpi2)
     
    
    sigma_x = 0.d0
    sigma_y = 0.d0
    sigma_z = 0.d0 
    
    sigma_x(1,2) = 1.d0 
    sigma_x(2,1) = 1.d0 
    
    sigma_z(1,1) = 1.d0 
    sigma_z(2,2) = -1.d0 
    
    sigma_y(1,2) = dcmplx(0.d0,-1.d0) 
    sigma_y(2,1) = dcmplx(0.d0, 1.d0) 
    
    DO i=-1,2
       IF ( (sfr-twompi(i)) <  0.0D0 ) sfr_heavyside(i) = 0.0D0
       IF ( (sfr-twompi(i)) >= 0.0D0 ) sfr_heavyside(i) = 1.0D0
    END DO
    

    !
    ! leading order contacts in PW
    !
    c1s0(-1) = -0.1513660372031080E+00
    c1s0(0)  = -0.1521410882366787E+00
    c1s0(1)  = -0.1517647459006913E+00

    c3s1(-1) = -0.1584341766228121E+00
    c3s1(0)  = -0.1584341766228121E+00
    c3s1(1)  = -0.1584341766228121E+00
    
    ! the LO CIBcontacts are input in units of 10^4/GeV^2
    ! this = 10^-2/MeV^2 
    c1s0 = c1s0*0.01d0 
    c3s1 = c3s1*0.01d0 
    
    CS = (c1s0 + 3.d0*c3s1) /16.d0/pi
    CT = (c3s1 - c1s0) /16.d0/pi
    
    !
    ! next-to-leading order contacts in PW
    !
    cnlo_pw = 0.d0 
    cnlo_pw(1) =  0.2404021944134705E+01
    cnlo_pw(2) =  0.1263390763475578E+01
    cnlo_pw(3) =  0.4170455420556486E+00
    cnlo_pw(4) = -0.7826584999752046E+00
    cnlo_pw(5) =  0.9283846626623043E+00
    cnlo_pw(6) =  0.6181414190474576E+00
    cnlo_pw(7) = -0.6778085114063558E+00
    
    ! NLO
    ! the NLO contacts are input in units of 10^4/GeV^4
    ! this = 10^-8/MeV^4
    DO i=1,7
       CNLO_PW(i) = CNLO_PW(i) * 1.E-08
    END DO
    
!!$    cnlo = 0.d0
!!$    cnlo(1) = 1.d0* 1.D-08/64./pi

    cnlo(1) = (-5.d0*cnlo_pw(7)+6.d0*cnlo_pw(5)-3.d0*cnlo_pw(4)-3.d0*cnlo_pw(3)-cnlo_pw(2)+2.d0*cnlo_pw(1))/(64.d0*pi)
    cnlo(2) = ( 5.d0*cnlo_pw(7)+6.d0*cnlo_pw(5)+3.d0*cnlo_pw(4)+3.d0*cnlo_pw(3)+cnlo_pw(2)+2.d0*cnlo_pw(1))/(16.d0*pi)
    cnlo(3) = -( 2.d0*cnlo_pw(7)-2.d0*dsqrt(2.d0)*cnlo_pw(6)-2.d0*cnlo_pw(5)-3.d0*cnlo_pw(3)+cnlo_pw(2)+2.d0*cnlo_pw(1))/(64.d0*pi)
    cnlo(4) = -(-2.d0*cnlo_pw(7)-2.d0*dsqrt(2.d0)*cnlo_pw(6)-2.d0*cnlo_pw(5)+3.d0*cnlo_pw(3)-cnlo_pw(2)+2.d0*cnlo_pw(1))/(16.d0*pi)
    cnlo(5) = -(-5.d0*cnlo_pw(7)+3.d0*cnlo_pw(4)+2.d0*cnlo_pw(2))/(16.d0*pi)
    cnlo(6) = ( cnlo_pw(7)-6.d0*dsqrt(2.d0)*cnlo_pw(6)-3.d0*cnlo_pw(4)+2.d0*cnlo_pw(2))/(64.d0*pi)
    cnlo(7) = -(cnlo_pw(7)+6.d0*dsqrt(2.d0)*cnlo_pw(6)-3.d0*cnlo_pw(4)+2.d0*cnlo_pw(2))/(16.d0*pi) 
    

    !
    ! 3NF CONTACT parameters, fit from Petr and Sofia to H3/He3 and H3 half life. 
    !  
    !   cD = -0.39 +/- 0.07, cE=-0.398 +0.015/-0.016

    cD = -0.39d0 
    cE = -0.398d0
    Econst = cE/fpi4/lambda 
    Dconst = cD/fpi2/lambda
    
    ! -----------------------------------------------------
  end subroutine init_chp_constants
  
  ! static one pion exchange, [1] eq 4.5
  ! without isospin structure
  ! q2  : momentum transfer squared
  ! impi: determines which mpi2 to use, 
  FUNCTION chp_one_pion_exchange(q2, impi) RESULT(res)
  
    IMPLICIT NONE 
    REAL*8 , INTENT(IN) :: q2
    INTEGER, INTENT(IN) :: impi
    REAL*8 :: res
    
    res = 0.0D0
    res = -1.0D0 * const(2)/ (q2 + mpi2(impi))
    
  END FUNCTION chp_one_pion_exchange
  

  ! NLO function Eq 4.9
  ! q2  : momentum transfer squared
  ! impi: determines which mpi2 to use, 
  FUNCTION chp_NLO_two_pion_exchange_Wc(q2, L, w, impi) RESULT(res)
    
    REAL*8 , INTENT(IN) :: q2
    REAL*8 , INTENT(IN) :: L
    REAL*8 , INTENT(IN) :: w
    INTEGER, INTENT(IN) :: impi
    REAL*8 :: res
    
    res = 0.0D0
    res = -L *(fourmpi2(impi)*const(4) + q2*const(5) + const(6)*mpi4(impi)/(w*w))/const(3) 
!    write(6,*) L,fourmpi2(impi),const(4)
    
  END FUNCTION chp_NLO_two_pion_exchange_Wc
  
  ! NLO function Eq 4.10
  ! q2  : momentum transfer squared
  ! impi: determines which mpi2 to use, 
  FUNCTION chp_NLO_two_pion_exchange_Vs(q2,L,impi) RESULT(res)
    
    REAL*8 , INTENT(IN) :: q2
    REAL*8 , INTENT(IN) :: L
    INTEGER, INTENT(IN) :: impi
    REAL*8 :: res
    
    res = 0.0D0
    res = +const(7)*L*q2/const(8)
    
  END FUNCTION chp_NLO_two_pion_exchange_Vs

  ! NLO Vt function Eq 4.10
  ! q2  : momentum transfer squared
  ! impi: determines which mpi2 to use, 
  FUNCTION chp_NLO_two_pion_exchange_VT(L,impi) RESULT(res)
    
    REAL*8 , INTENT(IN) :: L
    INTEGER, INTENT(IN) :: impi
    REAL*8 :: res
    
    res = 0.0D0 
    res = -const(7)*L/const(8)
    
  END FUNCTION chp_NLO_two_pion_exchange_VT

  ! NLO loop function w [1] Eq 4.12 (DR and SFR)
  ! q2  : momentum transfer squared
  ! impi: determines which mpi2 to use, 
  FUNCTION chp_NLO_two_pion_exchange_loop_w(q2, impi) RESULT(res)
    
    REAL*8 , INTENT(IN) :: q2
    INTEGER, INTENT(IN) :: impi
    REAL*8 :: res
    
    res = 0.0D0
    res = DSQRT(fourmpi2(impi) + q2)


    
  END FUNCTION chp_NLO_two_pion_exchange_loop_w

  ! NLO SFR loop function s [2] Eq 2.16
  ! impi: determines which mpi2 to use, 
  FUNCTION chp_NLO_sfr_two_pion_exchange_loop_s(impi) RESULT(res)
    
    INTEGER, INTENT(IN) :: impi
    REAL*8 :: res
    
    res = 0.0D0
    res = DSQRT(sfr2 - fourmpi2(impi))

  END FUNCTION chp_NLO_sfr_two_pion_exchange_loop_s
  
  ! NLO SFR loop function L [2] Eq 2.16
  ! q   : momentum transfer
  ! q2  : momentum transfer squared
  ! w,s : SFR loop functions
  ! impi: determines which mpi2 to use, 
  FUNCTION chp_NLO_sfr_two_pion_exchange_loop_L(q, q2, w, s, impi) RESULT(res)
    
    REAL*8 , INTENT(INOUT) :: q
    REAL*8 , INTENT(INOUT) :: q2
    REAL*8 , INTENT(IN) :: w
    REAL*8 , INTENT(IN) :: s
    INTEGER, INTENT(IN) :: impi
    REAL*8 :: res, eps 
    
    eps = 1.e-8 
    res = 0.0D0
    IF (sfr_heavyside(impi) == 0.0D0) return
    
    if ( q == 0.d0 ) then  
       q = q + eps 
       q2 = q*q 
    end if
    res = w * dlog( (sfr2*w*w+q2*s*s+2.0D0*sfr*q*w*s)/( fourmpi2(impi)*(sfr2+q2) ) )/(2.0D0*q)
    !write(6,'(3g)') q, q2, res !(sfr-twompi(i)), sfr, sfr2, q,q2 ! (sfr2*w*w+q2*s*s+2.0D0*sfr*q*w*s), ( fourmpi2(impi)*(sfr2+q2) )
    !write(6,*) 0.5d0*( 2.d0*s/sfr/dsqrt(fourmpi2(impi)) - sfr2*fourmpi2(impi)) * w 

  END FUNCTION chp_NLO_sfr_two_pion_exchange_loop_L



  FUNCTION chp_sigma_dot_sigma(ms1,ms2,ms3,ms4) RESULT(res)
    
    IMPLICIT NONE 
    INTEGER, INTENT(IN) :: ms1,ms2,ms3,ms4
    REAL*8 :: res
    complex*16 :: res1
    COMPLEX*16 :: chi1(2), chi2(2), chi3(2), chi4(2)
    INTEGER :: i1 
    
    res = 0.0D0
    chi1 = 0.d0; chi2 = 0.d0; chi3 = 0.d0; chi4 = 0.d0
    i1 = nint(1.5-0.5*ms1) 
    chi1(i1) = 1.d0 
    i1 = nint(1.5-0.5*ms2) 
    chi2(i1) = 1.d0 
    i1 = nint(1.5-0.5*ms3) 
    chi3(i1) = 1.d0 
    i1 = nint(1.5-0.5*ms4) 
    chi4(i1) = 1.d0 
    
 
    res1 = dot_product(chi1, matmul( sigma_x, chi3)) * dot_product(chi2, matmul( sigma_x, chi4)) &
         + dot_product(chi1, matmul( sigma_y, chi3)) * dot_product(chi2, matmul( sigma_y, chi4)) &
         + dot_product(chi1, matmul( sigma_z, chi3)) * dot_product(chi2, matmul( sigma_z, chi4)) 
    res = res1
    
    
  END FUNCTION chp_sigma_dot_sigma

  


!!$
!!$  FUNCTION chp_sigma_dot_sigma(ms1,ms2,ms3,ms4) RESULT(res)
!!$    
!!$    IMPLICIT NONE 
!!$    INTEGER, INTENT(IN) :: ms1,ms2,ms3,ms4
!!$    REAL*8 :: res, spin 
!!$    COMPLEX*16 :: res1, jx1, jx2, jy1, jy2, jz1, jz2 
!!$    REAL*8 :: m1,m2,m3,m4
!!$    INTEGER :: i1 
!!$    
!!$    m1 = 0.5d0*dble(ms1); m2 = 0.5d0*dble(ms2)
!!$    m3 = 0.5d0*dble(ms3); m4 = 0.5d0*dble(ms4)
!!$    
!!$    jx1 = 0.d0
!!$    jx2 = 0.d0 
!!$    jy1 = 0.d0
!!$    jy2 = 0.d0 
!!$    jz1 = 0.d0
!!$    jz2 = 0.d0 
!!$    
!!$    spin = 0.5d0 
!!$    res = 0.0D0
!!$    if ( ms1 == ms3 ) THEN
!!$       jx1 = 0.d0
!!$       jy1 = 0.d0 
!!$       jz1 = 2.d0*m1
!!$    end if
!!$    if ( ms1 == ms3 + 2 ) then 
!!$       jx1 = sqrt( (spin+m3+1.d0)*(spin-m3) )
!!$       jy1 = -dcmplx(0.d0,1.d0)*sqrt( (spin+m3+1.d0)*(spin-m3) )
!!$       jz1 = 0.d0
!!$    end if
!!$    if ( ms1 == ms3 - 2 ) then 
!!$       jx1 = sqrt( (spin-m3+1.d0)*(spin+m3) )
!!$       jy1 = dcmplx(0.d0,1.d0)*sqrt( (spin-m3+1.d0)*(spin+m3) )
!!$       jz1 = 0.d0
!!$    end if
!!$    
!!$    if ( ms2 == ms4 ) THEN
!!$       jx2 = 0.d0
!!$       jy2 = 0.d0 
!!$       jz2 = 2.d0*m2 
!!$    end if
!!$    if ( ms2 == ms4 + 2 ) then 
!!$       jx2 = sqrt( (spin+m4+1.d0)*(spin-m4) )
!!$       jy2 = -dcmplx(0.d0,1.d0)*sqrt( (spin+m4+1.d0)*(spin-m4) )
!!$       jz2 = 0.d0
!!$    end if
!!$    if ( ms2 == ms4 - 2 ) then 
!!$       jx2 = sqrt( (spin-m4+1.d0)*(spin+m4) )
!!$       jy2 = dcmplx(0.d0,1.d0)*sqrt( (spin-m4+1.d0)*(spin+m4) )
!!$       jz2 = 0.d0
!!$    end if
!!$            
!!$    res1 = jx1*jx2 + jy1*jy2 + jz1*jz2
!!$    res = res1
!!$    
!!$  END FUNCTION chp_sigma_dot_sigma

  
  FUNCTION chp_tau_dot_tau(ms1,ms2,ms3,ms4) RESULT(res)
    
    IMPLICIT NONE 
    INTEGER, INTENT(IN) :: ms1,ms2,ms3,ms4
    REAL*8 :: res
    complex*16 :: res1
    COMPLEX*16 :: chi1(2), chi2(2), chi3(2), chi4(2)
    INTEGER :: i1 

    res = 0.0D0
    chi1 = 0.d0; chi2 = 0.d0; chi3 = 0.d0; chi4 = 0.d0
    i1 = nint(1.5-0.5*ms1) 
    chi1(i1) = 1.d0 
    i1 = nint(1.5-0.5*ms2) 
    chi2(i1) = 1.d0 
    i1 = nint(1.5-0.5*ms3) 
    chi3(i1) = 1.d0 
    i1 = nint(1.5-0.5*ms4) 
    chi4(i1) = 1.d0 
    
 
    res1 = dot_product(chi1, matmul( sigma_x, chi3)) * dot_product(chi2, matmul( sigma_x, chi4)) &
         + dot_product(chi1, matmul( sigma_y, chi3)) * dot_product(chi2, matmul( sigma_y, chi4)) &
         + dot_product(chi1, matmul( sigma_z, chi3)) * dot_product(chi2, matmul( sigma_z, chi4)) 
    
    res = res1
        
  END FUNCTION chp_tau_dot_tau

  
!!$  FUNCTION chp_tau_dot_tau(ms1,ms2,ms3,ms4) RESULT(res)
!!$    
!!$    IMPLICIT NONE 
!!$    INTEGER, INTENT(IN) :: ms1,ms2,ms3,ms4
!!$    REAL*8 :: res, spin 
!!$    COMPLEX*16 :: res1, jx1, jx2, jy1, jy2, jz1, jz2 
!!$    REAL*8 :: m1,m2,m3,m4
!!$    INTEGER :: i1 
!!$    
!!$    m1 = 0.5d0*dble(ms1); m2 = 0.5d0*dble(ms2)
!!$    m3 = 0.5d0*dble(ms3); m4 = 0.5d0*dble(ms4)
!!$    
!!$    jx1 = 0.d0
!!$    jx2 = 0.d0 
!!$    jy1 = 0.d0
!!$    jy2 = 0.d0 
!!$    jz1 = 0.d0
!!$    jz2 = 0.d0 
!!$    
!!$    spin = 0.5d0 
!!$    res = 0.0D0
!!$    if ( ms1 == ms3 ) THEN
!!$       jx1 = 0.d0
!!$       jy1 = 0.d0 
!!$       jz1 = 2.d0*m1
!!$    end if
!!$    if ( ms1 == ms3 + 2 ) then 
!!$       jx1 = sqrt( (spin+m3+1.d0)*(spin-m3) )
!!$       jy1 = -dcmplx(0.d0,1.d0)*sqrt( (spin+m3+1.d0)*(spin-m3) )
!!$       jz1 = 0.d0
!!$    end if
!!$    if ( ms1 == ms3 - 2 ) then 
!!$       jx1 = sqrt( (spin-m3+1.d0)*(spin+m3) )
!!$       jy1 = dcmplx(0.d0,1.d0)*sqrt( (spin-m3+1.d0)*(spin+m3) )
!!$       jz1 = 0.d0
!!$    end if
!!$    
!!$    if ( ms2 == ms4 ) THEN
!!$       jx2 = 0.d0
!!$       jy2 = 0.d0 
!!$       jz2 = 2.d0*m2 
!!$    end if
!!$    if ( ms2 == ms4 + 2 ) then 
!!$       jx2 = sqrt( (spin+m4+1.d0)*(spin-m4) )
!!$       jy2 = -dcmplx(0.d0,1.d0)*sqrt( (spin+m4+1.d0)*(spin-m4) )
!!$       jz2 = 0.d0
!!$    end if
!!$    if ( ms2 == ms4 - 2 ) then 
!!$       jx2 = sqrt( (spin-m4+1.d0)*(spin+m4) )
!!$       jy2 = dcmplx(0.d0,1.d0)*sqrt( (spin-m4+1.d0)*(spin+m4) )
!!$       jz2 = 0.d0
!!$    end if
!!$            
!!$    res1 = jx1*jx2 + jy1*jy2 + jz1*jz2
!!$    res = res1
!!$    
!!$  END FUNCTION chp_tau_dot_tau
  
  
  FUNCTION chp_sigma1_dot_q_sigma2_dot_q(ms1,ms2,ms3,ms4,q) RESULT(res)

    IMPLICIT NONE 
    INTEGER, INTENT(IN) :: ms1,ms2,ms3,ms4
    REAL*8, INTENT(IN) :: q(3) 
    REAL*8 :: res
    COMPLEX*16 :: res1, res2 
    COMPLEX*16 :: chi1(2), chi2(2), chi3(2), chi4(2), mat(2,2) 
    INTEGER :: i1 
    
    chi1 = 0.d0; chi2 = 0.d0; chi3 = 0.d0; chi4 = 0.d0
    i1 = nint(1.5-0.5*ms1) 
    chi1(i1) = 1.d0 
    i1 = nint(1.5-0.5*ms2) 
    chi2(i1) = 1.d0 
    i1 = nint(1.5-0.5*ms3) 
    chi3(i1) = 1.d0 
    i1 = nint(1.5-0.5*ms4) 
    chi4(i1) = 1.d0 
    
    mat = q(1)*sigma_x+q(2)*sigma_y+q(3)*sigma_z 
    res1 = dot_product(chi1, matmul( mat, chi3))
    res2 = dot_product(chi2, matmul( mat, chi4))
    
    res = res1 * res2  

    
  END FUNCTION chp_sigma1_dot_q_sigma2_dot_q

!!$
!!$  FUNCTION chp_sigma1_dot_q_sigma2_dot_q(ms1,ms2,ms3,ms4,q) RESULT(res)
!!$
!!$    IMPLICIT NONE 
!!$    INTEGER, INTENT(IN) :: ms1,ms2,ms3,ms4
!!$    REAL*8, INTENT(IN) :: q(3) 
!!$    REAL*8 :: res, spin 
!!$    COMPLEX*16 :: res1,res2, jx1, jx2, jy1, jy2, jz1, jz2 
!!$    REAL*8 :: m1,m2,m3,m4
!!$    INTEGER :: i1 
!!$    
!!$    m1 = 0.5d0*dble(ms1); m2 = 0.5d0*dble(ms2)
!!$    m3 = 0.5d0*dble(ms3); m4 = 0.5d0*dble(ms4)
!!$    
!!$    jx1 = 0.d0
!!$    jx2 = 0.d0 
!!$    jy1 = 0.d0
!!$    jy2 = 0.d0 
!!$    jz1 = 0.d0
!!$    jz2 = 0.d0 
!!$    
!!$    spin = 0.5d0 
!!$    res = 0.0D0
!!$    if ( ms1 == ms3 ) THEN
!!$       jx1 = 0.d0
!!$       jy1 = 0.d0 
!!$       jz1 = 2.d0*m1
!!$    end if
!!$    if ( ms1 == ms3 + 2 ) then 
!!$       jx1 = sqrt( (spin+m3+1.d0)*(spin-m3) )
!!$       jy1 = -dcmplx(0.d0,1.d0)*sqrt( (spin+m3+1.d0)*(spin-m3) )
!!$       jz1 = 0.d0
!!$    end if
!!$    if ( ms1 == ms3 - 2 ) then 
!!$       jx1 = sqrt( (spin-m3+1.d0)*(spin+m3) )
!!$       jy1 = dcmplx(0.d0,1.d0)*sqrt( (spin-m3+1.d0)*(spin+m3) )
!!$       jz1 = 0.d0
!!$    end if
!!$    
!!$    if ( ms2 == ms4 ) THEN
!!$       jx2 = 0.d0
!!$       jy2 = 0.d0 
!!$       jz2 = 2.d0*m2 
!!$    end if
!!$    if ( ms2 == ms4 + 2 ) then 
!!$       jx2 = sqrt( (spin+m4+1.d0)*(spin-m4) )
!!$       jy2 = -dcmplx(0.d0,1.d0)*sqrt( (spin+m4+1.d0)*(spin-m4) )
!!$       jz2 = 0.d0
!!$    end if
!!$    if ( ms2 == ms4 - 2 ) then 
!!$       jx2 = sqrt( (spin-m4+1.d0)*(spin+m4) )
!!$       jy2 = dcmplx(0.d0,1.d0)*sqrt( (spin-m4+1.d0)*(spin+m4) )
!!$       jz2 = 0.d0
!!$    end if
!!$            
!!$    res1 = jx1*q(1) + jy1*q(2) + jz1*q(3)
!!$    res2 = jx2*q(1) + jy2*q(2) + jz2*q(3)
!!$    res = res1*res2 
!!$    
!!$  end FUNCTION chp_sigma1_dot_q_sigma2_dot_q
  
  FUNCTION chp_spin_dot_qxk(ms1,ms2,ms3,ms4,qxk) RESULT(res)

    IMPLICIT NONE 

    INTEGER, INTENT(IN) :: ms1,ms2,ms3,ms4
    REAL*8, INTENT(IN) :: qxk(3) 
    REAL*8 :: res, spin 
    COMPLEX*16 :: res1,res2, jx1, jx2, jy1, jy2, jz1, jz2 
    REAL*8 :: m1,m2,m3,m4
    INTEGER :: i1 
    
    m1 = 0.5d0*dble(ms1); m2 = 0.5d0*dble(ms2)
    m3 = 0.5d0*dble(ms3); m4 = 0.5d0*dble(ms4)
    
    jx1 = 0.d0
    jx2 = 0.d0 
    jy1 = 0.d0
    jy2 = 0.d0 
    jz1 = 0.d0
    jz2 = 0.d0 
    
    spin = 0.5d0 
    res = 0.0D0
    if ( ms1 == ms3 ) THEN
       jx1 = 0.d0
       jy1 = 0.d0 
       jz1 = 2.d0*m1
    end if
    if ( ms1 == ms3 + 2 ) then 
       jx1 = sqrt( (spin+m3+1.d0)*(spin-m3) )
       jy1 = -dcmplx(0.d0,1.d0)*sqrt( (spin+m3+1.d0)*(spin-m3) )
       jz1 = 0.d0
    end if
    if ( ms1 == ms3 - 2 ) then 
       jx1 = sqrt( (spin-m3+1.d0)*(spin+m3) )
       jy1 = dcmplx(0.d0,1.d0)*sqrt( (spin-m3+1.d0)*(spin+m3) )
       jz1 = 0.d0
    end if
    
    if ( ms2 == ms4 ) THEN
       jx2 = 0.d0
       jy2 = 0.d0 
       jz2 = 2.d0*m2 
    end if
    if ( ms2 == ms4 + 2 ) then 
       jx2 = sqrt( (spin+m4+1.d0)*(spin-m4) )
       jy2 = -dcmplx(0.d0,1.d0)*sqrt( (spin+m4+1.d0)*(spin-m4) )
       jz2 = 0.d0
    end if
    if ( ms2 == ms4 - 2 ) then 
       jx2 = sqrt( (spin-m4+1.d0)*(spin+m4) )
       jy2 = dcmplx(0.d0,1.d0)*sqrt( (spin-m4+1.d0)*(spin+m4) )
       jz2 = 0.d0
    end if
            
    !res1 = -dcmplx(0.d0,0.5d0) * ( (jx1+jx2)*qxk(1) + (jy1+jy2)*qxk(2) + (jz1+jz2)*qxk(3) )
    res1 = -( (jx1+jx2)*qxk(1) + (jy1+jy2)*qxk(2) + (jz1+jz2)*qxk(3) )
    
    res = res1 
    
    
  END FUNCTION chp_spin_dot_qxk
  
  
  ! regulator, eq 4.63 
  ! pfinal : final momentum
  ! pinit  : initial momentum
  ! n      : cutoff order 
  ! LAMBDA is accessed from the chp constant chp_lambda

  FUNCTION freg(pfinal, pinit, n) RESULT(res)
    
    REAL*8 , INTENT(IN) :: pfinal(3), pinit(3)
    REAL*8 , INTENT(IN) :: n
    REAL*8 :: res,exponent, p2, pp2 
    
    res = 0.0D0
    p2 = sum( pfinal*pfinal) 
    pp2 = sum( pinit*pinit )
    
    exponent = (p2**n/lambda**(2.d0*n) + &
         pp2**n/lambda**(2.0D0*n) )
    
    res = dexp(-exponent)
    
  END FUNCTION freg
  
  ! local 3nf regulator, eq 11 in P. Navratil's paper
  ! pfinal : final momentum
  ! pinit  : initial momentum
  ! n      : cutoff order 
  ! LAMBDA is accessed from the chp constant chp_lambda
  FUNCTION freg_3nflocal(q2) RESULT(res)
    
    REAL*8 , INTENT(IN) :: q2
    REAL*8 :: res,exponent
    
    res = 0.0D0
    exponent = (q2/lambda**2)**2    
    res = dexp(-exponent)
    
  END FUNCTION freg_3nflocal
  
  
  
  FUNCTION chp_3nf_Dterm(ms1,ms2,ms3,ms4,q) RESULT(res)
    
    IMPLICIT NONE 
    INTEGER, INTENT(IN) :: ms1,ms2,ms3,ms4
    REAL*8, INTENT(IN) :: q(3) 
    REAL*8 :: res, q2
    COMPLEX*16 :: res1
    COMPLEX*16 :: chi1(2), chi2(2), chi3(2), chi4(2), mat1_x(2,2) , mat1_y(2,2) , mat1_z(2,2), mat(2,2)  
    INTEGER :: i1 
    
    chi1 = 0.d0; chi2 = 0.d0; chi3 = 0.d0; chi4 = 0.d0
    i1 = nint(1.5-0.5*ms1) 
    chi1(i1) = 1.d0 
    i1 = nint(1.5-0.5*ms2) 
    chi2(i1) = 1.d0 
    i1 = nint(1.5-0.5*ms3) 
    chi3(i1) = 1.d0 
    i1 = nint(1.5-0.5*ms4) 
    chi4(i1) = 1.d0 
    
    q2 = sum( q*q )
    mat =  q(1)*sigma_x+q(2)*sigma_y+q(3)*sigma_z
    mat1_x = matmul(sigma_x, mat)
    mat1_y = matmul(sigma_y, mat)
    mat1_z = matmul(sigma_z, mat)
    
    
    res1 = dot_product(chi1, matmul( mat1_x, chi3))*dot_product(chi2, matmul( sigma_x, chi4)) + & 
         dot_product(chi1, matmul( mat1_y, chi3))*dot_product(chi2, matmul( sigma_y, chi4)) + & 
         dot_product(chi1, matmul( mat1_z, chi3))*dot_product(chi2, matmul( sigma_z, chi4)) 
    
    res = -Dconst*const(21)*res1/( q2 + mpi2(2))
    
  END FUNCTION chp_3nf_Dterm

  
  FUNCTION chp_sigma1_dot_q_ope(ms1,ms2,q) RESULT(res)

    IMPLICIT NONE 
    INTEGER, INTENT(IN) :: ms1,ms2
    REAL*8, INTENT(IN) :: q(3) 
    REAL*8 :: res, q2
    COMPLEX*16 :: res1 
    COMPLEX*16 :: chi1(2), chi2(2), mat(2,2) 
    INTEGER :: i1 
    
    chi1 = 0.d0; chi2 = 0.d0
    i1 = nint(1.5-0.5*ms1) 
    chi1(i1) = 1.d0 
    i1 = nint(1.5-0.5*ms2) 
    chi2(i1) = 1.d0 
    
    mat = q(1)*sigma_x+q(2)*sigma_y+q(3)*sigma_z 
    res1 = dot_product(chi1, matmul( mat, chi2))
    
    q2 = sum( q*q )
    res = res1 /( q2 + mpi2(2))

    
  END FUNCTION chp_sigma1_dot_q_ope
  
  
  FUNCTION chp_sigma1_dot_qXq(ms1,ms2,q2,q3) RESULT(res)

    IMPLICIT NONE 
    INTEGER, INTENT(IN) :: ms1,ms2
    REAL*8, INTENT(IN) :: q2(3), q3(3) 
    REAL*8 :: res, qxq(3)
    COMPLEX*16 :: chi1(2), chi2(2), mat(2,2) 
    INTEGER :: i1 
    
    qxq(1) = q2(2)*q3(3)-q2(3)*q3(2)
    qxq(2) = q2(3)*q3(1)-q2(1)*q3(3)
    qxq(3) = q2(1)*q3(2)-q2(2)*q3(1)
    
    chi1 = 0.d0; chi2 = 0.d0
    i1 = nint(1.5-0.5*ms1) 
    chi1(i1) = 1.d0 
    i1 = nint(1.5-0.5*ms2) 
    chi2(i1) = 1.d0 
    
    mat = qxq(1)*sigma_x+qxq(2)*sigma_y+qxq(3)*sigma_z 
    res = dot_product(chi1, matmul( mat, chi2))
        
  END FUNCTION chp_sigma1_dot_qXq
  
  FUNCTION chp_tau1_dot_tauXtau(t1,t2,t3,t4,t5,t6) RESULT(res)

    IMPLICIT NONE 
    INTEGER, INTENT(IN) :: t1,t2,t3,t4,t5,t6
    REAL*8 :: res, qxq(3)
    COMPLEX*16 :: chi1(2), chi2(2), chi3(2), chi4(2), chi5(2), chi6(2), & 
         mat(2,2), facx, facy, facz  
    INTEGER :: i1 
    
    chi1 = 0.d0; chi2 = 0.d0; chi3 = 0.d0; chi4 = 0.d0; chi5 = 0.d0; chi6 = 0.d0
    i1 = nint(1.5-0.5*t1) 
    chi1(i1) = 1.d0 
    i1 = nint(1.5-0.5*t2) 
    chi2(i1) = 1.d0 
    i1 = nint(1.5-0.5*t3) 
    chi3(i1) = 1.d0 
    i1 = nint(1.5-0.5*t4) 
    chi4(i1) = 1.d0 
    i1 = nint(1.5-0.5*t5) 
    chi5(i1) = 1.d0 
    i1 = nint(1.5-0.5*t6) 
    chi6(i1) = 1.d0 
    
    facx = dot_product(chi3,matmul(sigma_y,chi5))*dot_product(chi3,matmul(sigma_z,chi6))  - &
         dot_product(chi3,matmul(sigma_z,chi5))*dot_product(chi3,matmul(sigma_y,chi6))  

    facy = dot_product(chi3,matmul(sigma_z,chi5))*dot_product(chi3,matmul(sigma_x,chi6))  - &
         dot_product(chi3,matmul(sigma_x,chi5))*dot_product(chi3,matmul(sigma_z,chi6))  

    facz = dot_product(chi3,matmul(sigma_x,chi5))*dot_product(chi3,matmul(sigma_y,chi6))  - &
         dot_product(chi3,matmul(sigma_y,chi5))*dot_product(chi3,matmul(sigma_x,chi6))  


    mat = facx*sigma_x+facy*sigma_y+facz*sigma_z 
    res = dot_product(chi1, matmul( mat, chi2))
        
  END FUNCTION chp_tau1_dot_tauXtau
  
  

end module chiral_constants



module chiral_potentials 

  implicit none 
  
  TYPE, PUBLIC :: chp_int_type
     INTEGER          :: val
     CHARACTER(LEN=12):: name
     LOGICAL          :: set
  END TYPE chp_int_type

  TYPE, PUBLIC :: chp_real_type
     REAL*8           :: val
     CHARACTER(LEN=12):: name
     LOGICAL          :: set
  END TYPE chp_real_type

  TYPE(chp_int_type), PRIVATE :: chp_chiral_order ! LO, NLO, NNLO, N3LO
  TYPE(chp_real_type), PRIVATE :: chp_regcut
  
  
contains 
  
  real*8 function chiral_pot(p,q,r,s)
    USE single_particle_orbits
    USE constants
    use chiral_constants
    use ang_mom_functions !, only: tjs 
    
    implicit none 
    INTEGER :: p,q,r,s, m1,m2,m3,m4, spin, iph, t1,t2,t3,t4, Tiso
    INTEGER :: nx1, ny1, nz1, nx2, ny2, nz2, nx3, ny3, nz3, nx4, ny4, nz4 
    REAL*8 :: k1(3), k2(3), k3(3), k4(3), kmean(3) 
    REAL*8 :: qtrans(3), prel(3), pprel(3), qxk(3)
    REAL*8 :: q2, p2, kmean2, qabs, pp2, vdir, vexc, cg1, cg2, sum1
    REAL*8 :: delta, nucleon_mass, relativity_factor, dij
    REAL*8 :: VLO, VNLO, VNNLO, VN3LO 
    ! LIST OF OPERATOR STRUCTURES
    REAL*8 :: Vc
    REAL*8 :: Wc
    REAL*8 :: Vs
    REAL*8 :: Ws
    REAL*8 :: VLS
    REAL*8 :: WLS
    REAL*8 :: VsigL
    REAL*8 :: WsigL
    REAL*8 :: VT
    REAL*8 :: WT
    REAL*8 :: Vsigk
    REAL*8 :: Wsigk
    
    REAL*8 :: loop_w
    REAL*8 :: loop_s
    REAL*8 :: loop_L
    
    !
    ! set chiral order 
    !
    chp_chiral_order%val = LO 
    ! set the regulator cutoffs depending on the chiral order
    IF (chp_chiral_order%val == LO   ) chp_regcut%val = 2.0D0
    IF (chp_chiral_order%val == NLO  ) chp_regcut%val = 2.0D0
    IF (chp_chiral_order%val == NNLO ) chp_regcut%val = 3.0D0
    IF (chp_chiral_order%val == N3LO ) chp_regcut%val = 2.0D0
    

    Vc      = 0.0D0
    Wc      = 0.0D0
    Vs      = 0.0D0
    Ws      = 0.0D0
    VLS     = 0.0D0
    WLS     = 0.0D0
    VsigL   = 0.0D0
    WsigL   = 0.0D0
    VT      = 0.0D0
    WT      = 0.0D0
    Vsigk   = 0.0D0
    Wsigk   = 0.0D0
  
    vlo = 0.0D0 ; vnlo = 0.0D0 ; vnnlo = 0.0D0 ; vn3lo = 0.0D0

  
    nx1 = all_orbit%nx(p)
    ny1 = all_orbit%ny(p)
    nz1 = all_orbit%nz(p)
    nx2 = all_orbit%nx(q)
    ny2 = all_orbit%ny(q)
    nz2 = all_orbit%nz(q)
    nx3 = all_orbit%nx(r)
    ny3 = all_orbit%ny(r)
    nz3 = all_orbit%nz(r)
    nx4 = all_orbit%nx(s)
    ny4 = all_orbit%ny(s)
    nz4 = all_orbit%nz(s)
    ! 
    ! Conservation of linear momentum
    !
    if ( nx1 + nx2 /= nx3 + nx4 ) return 
    if ( ny1 + ny2 /= ny3 + ny4 ) return 
    if ( nz1 + nz2 /= nz3 + nz4 ) return 

    k1(1) = all_orbit%kx(p)
    k1(2) = all_orbit%ky(p)
    k1(3) = all_orbit%kz(p)
    k2(1) = all_orbit%kx(q)
    k2(2) = all_orbit%ky(q)
    k2(3) = all_orbit%kz(q)
    k3(1) = all_orbit%kx(r)
    k3(2) = all_orbit%ky(r)
    k3(3) = all_orbit%kz(r)
    k4(1) = all_orbit%kx(s)
    k4(2) = all_orbit%ky(s)
    k4(3) = all_orbit%kz(s)
    
    ! momenta in MeV
    k1 = k1*hbarc
    k2 = k2*hbarc
    k3 = k3*hbarc
    k4 = k4*hbarc

    ! 
    ! conservation of spin and isospin 
    !
    if ( all_orbit%szp(p) + all_orbit%szp(q) /= all_orbit%szp(r) + all_orbit%szp(s) ) return 
    if ( all_orbit%itzp(p) + all_orbit%itzp(q) /= all_orbit%itzp(r) + all_orbit%itzp(s) ) return 
  
    m1 = all_orbit%szp(p) 
    m2 = all_orbit%szp(q) 
    m3 = all_orbit%szp(r) 
    m4 = all_orbit%szp(s) 
  
    t1 = all_orbit%itzp(p) 
    t2 = all_orbit%itzp(q) 
    t3 = all_orbit%itzp(r) 
    t4 = all_orbit%itzp(s) 
  
    ! 
    ! RELATIVE MOMENTA <prel |v| pprel > 
    ! 
    prel  = 0.5d0*(k1-k2)
    pprel = 0.5d0*(k3-k4)
    p2  = sum(prel*prel) 
    pp2 = sum(pprel*pprel)
    !
    ! AVERAGE MOMENTA kav = 1/2(prel+pprel)
    !
    kmean = 0.5d0*( prel + pprel ) 
    kmean2 = sum(kmean*kmean)

    nucleon_mass = mnuc((t1+t2)/2)
    relativity_factor = nucleon_mass/ & 
         sqrt( sqrt(( nucleon_mass**2 + p2) * (nucleon_mass**2 + pp2 )) )
    
    !
    ! momentum transfer 
    !
    qtrans = prel - pprel
    q2 = sum(qtrans*qtrans) 
    qabs = dsqrt(q2)
    
    !
    !  cross product between momentum transfer and average momenta q X k 
    !
    qxk(1) = qtrans(2)*kmean(3)-qtrans(3)*kmean(2) 
    qxk(2) = qtrans(3)*kmean(1)-qtrans(1)*kmean(3) 
    qxk(3) = qtrans(1)*kmean(2)-qtrans(2)*kmean(1) 
    
    
    ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    !
    ! LEADING ORDER 
    ! 
    ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    
    !
    ! leading order one-pion exchange 
    ! 
    vlo = chp_sigma1_dot_q_sigma2_dot_q(m1,m2,m3,m4,qtrans)*chp_tau_dot_tau(t1,t2,t3,t4)*chp_one_pion_exchange(q2, 2)
    
!!$    WT = 0.d0 
!!$    DO Tiso = 0, 2, 2
!!$       if ( Tiso == 0 ) sum1 = -3.d0 
!!$       if ( Tiso == 2 ) sum1 =  1.d0 
!!$       
!!$       cg1 = tjs(1,1,Tiso,t1,t2,-(t1+t2))*iph( (t1+t2)/2 )*sqrt(Tiso+1.d0)
!!$       cg2 = tjs(1,1,Tiso,t3,t4,-(t3+t4))*iph( (t3+t4)/2 )*sqrt(Tiso+1.d0)
!!$       WT = WT + cg1*cg2*sum1  
!!$       
!!$    END DO
    
    !vlo = WT*chp_sigma1_dot_q_sigma2_dot_q(m1,m2,m3,m4,qtrans)*chp_one_pion_exchange(q2, 2)!*delta(t1,t3)*delta(t2,t4)!*delta(m1,m3)*delta(m2,m4)/3.
    
!!$    WT = 0.d0 
!!$    ! LO CIB one-pion exchanges
!!$    ! [1] Eq. 4.77-4.79
!!$    IF (t1 + t2 == 0) THEN 
!!$       WT = - chp_one_pion_exchange(q2, 0) 
!!$       DO Tiso = 0, 2, 2
!!$          cg1 = tjs(1,1,Tiso,t1,t2,-(t1+t2))*iph( (t1+t2)/2 )*sqrt(Tiso+1.d0)
!!$          cg2 = tjs(1,1,Tiso,t3,t4,-(t3+t4))*iph( (t3+t4)/2 )*sqrt(Tiso+1.d0)
!!$          WT = WT + cg1*cg2*iph(Tiso/2+1)*2.0D0*chp_one_pion_exchange(q2,1)
!!$       END DO
!!$    END IF
!!$    IF (t1 + t2 /= 0) THEN
!!$       WT = chp_one_pion_exchange(q2, 0) *chp_tau_dot_tau(t1,t2,t3,t4)
!!$    END IF
!!$    

    
    !vlo = chp_sigma1_dot_q_sigma2_dot_q(m1,m2,m3,m4,qtrans)*WT 
        
    !vlo=vlo +
    !
    ! leading order contacts 
    !
    vlo = vlo + ( CS((t1+t2)/2)*delta(m1,m3)*delta(m2,m4)*delta(t1,t3)*delta(t2,t4) + & 
         CT((t1+t2)/2) * chp_sigma_dot_sigma(m1,m2,m3,m4) *delta(t1,t3)*delta(t2,t4) )  *10.d0
    

    
    ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    !
    ! NEXT TO LEADING ORDER 
    ! 
    ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    
    !
    ! NLO pion-exchange 
    ! 
    ! irreducible, or non-polynomial, NLO two-pion exchanges
    
    loop_w = chp_NLO_two_pion_exchange_loop_w(q2,2)
    loop_s = chp_NLO_sfr_two_pion_exchange_loop_s(2)
    loop_L = chp_NLO_sfr_two_pion_exchange_loop_L(qabs, q2, loop_w, loop_s, 2)
    
    

    Wc = chp_tau_dot_tau(t1,t2,t3,t4)*delta(m1,m3)*delta(m2,m4)* & 
         chp_NLO_two_pion_exchange_Wc(q2, loop_L, loop_w, 2)
    Vs = chp_NLO_two_pion_exchange_Vs(q2, loop_L, 2)* &
         chp_sigma_dot_sigma(m1,m2,m3,m4)*delta(t1,t3)*delta(t2,t4) 
    VT = chp_NLO_two_pion_exchange_VT(loop_L, 2)* & 
         chp_sigma1_dot_q_sigma2_dot_q(m1,m2,m3,m4,qtrans)*delta(t1,t3)*delta(t2,t4) 
    
    vnlo = ( WC + VS + VT ) 
    !
    ! next-to-leading order contacts 
    !
    vnlo = vnlo + cnlo(1)*q2*delta(m1,m3)*delta(m2,m4)*delta(t1,t3)*delta(t2,t4) + & 
         cnlo(2)*kmean2*delta(m1,m3)*delta(m2,m4)*delta(t1,t3)*delta(t2,t4) + & 
         ( cnlo(3)*q2 + cnlo(4)*kmean2 )*chp_sigma_dot_sigma(m1,m2,m3,m4)*delta(t1,t3)*delta(t2,t4) + & 
         cnlo(5)*chp_spin_dot_qxk(m1,m2,m3,m4,qxk) + & 
         cnlo(6)*chp_sigma1_dot_q_sigma2_dot_q(m1,m2,m3,m4,qtrans)*delta(t1,t3)*delta(t2,t4) + & 
         cnlo(7)*chp_sigma1_dot_q_sigma2_dot_q(m1,m2,m3,m4,kmean)*delta(t1,t3)*delta(t2,t4) 
    
    !
    ! sum up all orders 
    !
    !
    ! regulator and relativity factor 
    !
    !vdir = (vlo+vnlo+vnnlo+vn3lo) * relativity_factor * freg(prel, pprel, chp_regcut%val) 
    
    vdir = (vlo + vnlo) * relativity_factor * freg(prel, pprel, chp_regcut%val) 
    
    !
    ! SUBTRACT TWO-BODY CENTER-OF-MASS KINETIC ENERGY
    !
    chiral_pot = hbarc**3 * (vdir)/volume - sum(k1*k2)*delta(p,r) * delta(q,s)/p_mass/below_ef! + 0.001* ( 1.5d0*delta(t1,t3)*delta(t2,t4) & 
!         + 2.d0*chp_tau_dot_tau(t1,t2,t3,t4))
        
    
  end function chiral_pot

  real*8 function chiral_3nf(p,q,r,s,t,u)
    USE single_particle_orbits
    USE constants
    use chiral_constants
    use ang_mom_functions, only : tjs 
    
    implicit none 
    INTEGER :: p,q,r,s,t,u, m1,m2,m3,m4,m5,m6, spin, iph, t1,t2,t3,t4,t5,t6, Tiso
    INTEGER :: nx1, ny1, nz1, nx2, ny2, nz2, nx3, ny3, nz3, nx4, ny4, nz4, nx5, ny5, nz5 , nx6, ny6, nz6 
    REAL*8 :: k1(3), k2(3), k3(3), k4(3), k5(3), k6(3), kmean(3) 
    REAL*8 :: qtrans1(3), qtrans2(3), qtrans3(3), prel(3), pprel(3)
    REAL*8 :: q22, q32, q12
    REAL*8 :: delta, nucleon_mass
    REAL*8 :: VNNLO
    
    vnnlo = 0.0D0 
    
  
    nx1 = all_orbit%nx(p)
    ny1 = all_orbit%ny(p)
    nz1 = all_orbit%nz(p)
    nx2 = all_orbit%nx(q)
    ny2 = all_orbit%ny(q)
    nz2 = all_orbit%nz(q)
    nx3 = all_orbit%nx(r)
    ny3 = all_orbit%ny(r)
    nz3 = all_orbit%nz(r)
    nx4 = all_orbit%nx(s)
    ny4 = all_orbit%ny(s)
    nz4 = all_orbit%nz(s)
    nx5 = all_orbit%nx(t)
    ny5 = all_orbit%ny(t)
    nz5 = all_orbit%nz(t)
    nx6 = all_orbit%nx(u)
    ny6 = all_orbit%ny(u)
    nz6 = all_orbit%nz(u)
    ! 
    ! Conservation of linear momentum
    !
    if ( nx1 + nx2 + nx3 /= nx4 + nx5 + nx6 ) return 
    if ( ny1 + ny2 + ny3 /= ny4 + ny5 + ny6 ) return 
    if ( nz1 + nz2 + nz3 /= nz4 + nz5 + nz6 ) return 
    

    k1(1) = all_orbit%kx(p)
    k1(2) = all_orbit%ky(p)
    k1(3) = all_orbit%kz(p)
    k2(1) = all_orbit%kx(q)
    k2(2) = all_orbit%ky(q)
    k2(3) = all_orbit%kz(q)
    k3(1) = all_orbit%kx(r)
    k3(2) = all_orbit%ky(r)
    k3(3) = all_orbit%kz(r)
    k4(1) = all_orbit%kx(s)
    k4(2) = all_orbit%ky(s)
    k4(3) = all_orbit%kz(s)
    k5(1) = all_orbit%kx(t)
    k5(2) = all_orbit%ky(t)
    k5(3) = all_orbit%kz(t)
    k6(1) = all_orbit%kx(u)
    k6(2) = all_orbit%ky(u)
    k6(3) = all_orbit%kz(u)
    
    ! momenta in MeV
    k1 = k1*hbarc
    k2 = k2*hbarc
    k3 = k3*hbarc
    k4 = k4*hbarc
    k5 = k5*hbarc
    k6 = k6*hbarc

    ! 
    ! conservation of spin and isospin 
    !
    if ( all_orbit%szp(p) + all_orbit%szp(q)+ all_orbit%szp(r) /= & 
         all_orbit%szp(s) + all_orbit%szp(t) + all_orbit%szp(u) ) return 
    if ( all_orbit%itzp(p) + all_orbit%itzp(q)+ all_orbit%itzp(r) /= & 
         all_orbit%itzp(s) + all_orbit%itzp(t)+ all_orbit%itzp(u) ) return 
    
    m1 = all_orbit%szp(p) 
    m2 = all_orbit%szp(q) 
    m3 = all_orbit%szp(r) 
    m4 = all_orbit%szp(s) 
    m5 = all_orbit%szp(s) 
    m6 = all_orbit%szp(s) 
  
    t1 = all_orbit%itzp(p) 
    t2 = all_orbit%itzp(q) 
    t3 = all_orbit%itzp(r) 
    t4 = all_orbit%itzp(s) 
    t5 = all_orbit%itzp(t) 
    t6 = all_orbit%itzp(u) 
  
    !
    ! Momentum transfer of particle 2 and 3 
    !
    qtrans1 = k4-k1
    qtrans2 = k5-k2
    qtrans3 = k6-k3
    
    Q12 = sum(qtrans1*qtrans1)
    Q22 = sum(qtrans2*qtrans2)
    Q32 = sum(qtrans3*qtrans3)
    
    !
    ! NNLO cE contact term 
    !
    vnnlo = Econst * ( chp_tau_dot_tau(t1,t4,t2,t5)*delta(m1,m4)*delta(m2,m5)*delta(r,u) * & 
         freg_3nflocal(q22)*freg_3nflocal(q12) + & 
         chp_tau_dot_tau(t1,t4,t3,t6)*delta(m1,m4)*delta(m3,m6)*delta(q,t) * & 
         freg_3nflocal(q32)*freg_3nflocal(q12) + & 
         chp_tau_dot_tau(t2,t5,t3,t6)*delta(m2,m5)*delta(m3,m6)*delta(p,s) * & 
         freg_3nflocal(q22)*freg_3nflocal(q32))
    
    !
    ! NNLO cD contact term 
    !
    vnnlo = vnnlo + & 
         chp_3nf_Dterm(m1,m2,m4,m5,qtrans1)*chp_tau_dot_tau(t1,t2,t4,t5) * & 
         freg_3nflocal(q22)*freg_3nflocal(q12) + &
         chp_3nf_Dterm(m1,m3,m4,m6,qtrans1)*chp_tau_dot_tau(t1,t3,t4,t6) * & 
         freg_3nflocal(q32)*freg_3nflocal(q12) + &
         chp_3nf_Dterm(m1,m2,m4,m5,qtrans2)*chp_tau_dot_tau(t1,t2,t4,t5) * & 
         freg_3nflocal(q22)*freg_3nflocal(q12) + &
         chp_3nf_Dterm(m2,m3,m5,m6,qtrans2)*chp_tau_dot_tau(t2,t3,t5,t6) * & 
         freg_3nflocal(q22)*freg_3nflocal(q32) + &
         chp_3nf_Dterm(m1,m3,m4,m6,qtrans3)*chp_tau_dot_tau(t1,t3,t4,t6) * & 
         freg_3nflocal(q32)*freg_3nflocal(q12) + &
         chp_3nf_Dterm(m2,m3,m5,m6,qtrans3)*chp_tau_dot_tau(t2,t3,t5,t6) * & 
         freg_3nflocal(q22)*freg_3nflocal(q32) 
  
    !
    ! NNLO two-pion exchange
    !
    vnnlo = vnnlo - c1 *const(2)*(4.d0*mpi2(2)**2/fpi2) * ( & 
         chp_sigma1_dot_q_ope(m2,m5,qtrans2) * chp_sigma1_dot_q_ope(m3,m6,qtrans3) * &
         chp_tau_dot_tau(t2,t5,t3,t6)*delta(p,s) * freg_3nflocal(q22)*freg_3nflocal(q32) +  & 
         chp_sigma1_dot_q_ope(m2,m5,qtrans2) * chp_sigma1_dot_q_ope(m1,m4,qtrans1) * &
         chp_tau_dot_tau(t2,t5,t1,t4)*delta(r,u) * freg_3nflocal(q22)*freg_3nflocal(q12) +  & 
         chp_sigma1_dot_q_ope(m1,m4,qtrans1) * chp_sigma1_dot_q_ope(m3,m6,qtrans3) * &
         chp_tau_dot_tau(t1,t4,t3,t6)*delta(q,t) * freg_3nflocal(q12)*freg_3nflocal(q32) )
    
    vnnlo = vnnlo + c3 *const(2)*(2.d0/fpi2) * ( & 
         chp_sigma1_dot_q_ope(m2,m5,qtrans2) * chp_sigma1_dot_q_ope(m3,m6,qtrans3)* & 
         chp_tau_dot_tau(t2,t5,t3,t6)*delta(p,s)*sum(qtrans2*qtrans3) * freg_3nflocal(q22)*freg_3nflocal(q32) +  & 
         chp_sigma1_dot_q_ope(m2,m5,qtrans2) * chp_sigma1_dot_q_ope(m1,m4,qtrans1)* & 
         chp_tau_dot_tau(t2,t5,t1,t4)*delta(r,u)*sum(qtrans2*qtrans1) * freg_3nflocal(q22)*freg_3nflocal(q12) +  & 
         chp_sigma1_dot_q_ope(m1,m4,qtrans1) * chp_sigma1_dot_q_ope(m3,m6,qtrans3)* & 
         chp_tau_dot_tau(t1,t4,t3,t6)*delta(q,t)*sum(qtrans1*qtrans3) * freg_3nflocal(q12)*freg_3nflocal(q32) )
    
    vnnlo = vnnlo + c4 *const(2)*(1.d0/fpi2) * ( &
         chp_sigma1_dot_q_ope(m2,m5,qtrans2) * chp_sigma1_dot_q_ope(m3,m6,qtrans3) * & 
         chp_tau1_dot_tauXtau(t1,t2,t3,t4,t5,t6) *chp_sigma1_dot_qXq(m1,m4,qtrans2,qtrans3) * & 
         freg_3nflocal(q22)*freg_3nflocal(q32) + &  
         chp_sigma1_dot_q_ope(m2,m5,qtrans2) * chp_sigma1_dot_q_ope(m1,m4,qtrans1) * & 
         chp_tau1_dot_tauXtau(t3,t2,t1,t6,t5,t4) *chp_sigma1_dot_qXq(m3,m6,qtrans2,qtrans1) * & 
         freg_3nflocal(q22)*freg_3nflocal(q12) + &  
         chp_sigma1_dot_q_ope(m1,m4,qtrans1) * chp_sigma1_dot_q_ope(m3,m6,qtrans3) * & 
         chp_tau1_dot_tauXtau(t2,t1,t3,t5,t4,t6) *chp_sigma1_dot_qXq(m2,m5,qtrans1,qtrans3) * & 
         freg_3nflocal(q12)*freg_3nflocal(q32) )
    
    
    chiral_3nf = hbarc**6*vnnlo/volume**2
    
  end function chiral_3nf
  !
  ! Exact deuteron energy -2.18266175748238
  !
  real*8 function vmom_minnesota(p,q,r,s)
    USE single_particle_orbits
    USE constants
    use chiral_constants
  
    implicit none 
    INTEGER :: p,q,r,s, m1,m2,m3,m4, spin, iph, t1,t2,t3,t4, Tiso
    INTEGER :: nx1, ny1, nz1, nx2, ny2, nz2, nx3, ny3, nz3, nx4, ny4, nz4 
    REAL*8 :: k1(3), k2(3), k3(3), k4(3), kmean(3) 
    REAL*8 :: qtrans(3), prel(3), pprel(3)
    REAL*8 :: q2, p2, qabs, pp2, vdir, vexc, cg1, cg2 
    REAL*8 :: delta, nucleon_mass, relativity_factor
    REAL*8 :: v0r, v0t, v0s, kr,kt,ks, vr,vt,vs,vcentral, vsigma, spin_exc1, spin_exc2
  
    nx1 = all_orbit%nx(p)
    ny1 = all_orbit%ny(p)
    nz1 = all_orbit%nz(p)
    nx2 = all_orbit%nx(q)
    ny2 = all_orbit%ny(q)
    nz2 = all_orbit%nz(q)
    nx3 = all_orbit%nx(r)
    ny3 = all_orbit%ny(r)
    nz3 = all_orbit%nz(r)
    nx4 = all_orbit%nx(s)
    ny4 = all_orbit%ny(s)
    nz4 = all_orbit%nz(s)
    ! 
    ! Conservation of linear momentum
    !
    if ( nx1 + nx2 /= nx3 + nx4 ) return 
    if ( ny1 + ny2 /= ny3 + ny4 ) return 
    if ( nz1 + nz2 /= nz3 + nz4 ) return 
  
    k1(1) = all_orbit%kx(p)
    k1(2) = all_orbit%ky(p)
    k1(3) = all_orbit%kz(p)
    k2(1) = all_orbit%kx(q)
    k2(2) = all_orbit%ky(q)
    k2(3) = all_orbit%kz(q)
    k3(1) = all_orbit%kx(r)
    k3(2) = all_orbit%ky(r)
    k3(3) = all_orbit%kz(r)
    k4(1) = all_orbit%kx(s)
    k4(2) = all_orbit%ky(s)
    k4(3) = all_orbit%kz(s)
  
  
    ! 
    ! conservation of spin and isospin 
    !
    if ( all_orbit%szp(p) + all_orbit%szp(q) /= all_orbit%szp(r) + all_orbit%szp(s) ) return 
    if ( all_orbit%itzp(p) + all_orbit%itzp(q) /= all_orbit%itzp(r) + all_orbit%itzp(s) ) return 
  
    m1 = all_orbit%szp(p) 
    m2 = all_orbit%szp(q) 
    m3 = all_orbit%szp(r) 
    m4 = all_orbit%szp(s) 
  
    t1 = all_orbit%itzp(p) 
    t2 = all_orbit%itzp(q) 
    t3 = all_orbit%itzp(r) 
    t4 = all_orbit%itzp(s) 
  
    ! 
    ! RELATIVE MOMENTA <prel |v| pprel > 
    ! 
    prel  = 0.5d0*(k1-k2)
    pprel = 0.5d0*(k3-k4)
    !
    ! momentum transfer 
    !
    qtrans = prel - pprel
    q2 = sum(qtrans*qtrans) 
    
    
    v0r=200.0  ! MeV
    v0t=300.d0 !178.0  ! MeV
    v0s=300.d0 !91.85  ! MeV
    kr=1.487  ! fm**-2
    kt=0.639  ! fm**-2
    ks=0.465  ! fm**-2
    
    
    ! r-space 
    !vr=v0r*exp(-kr*rr**2)
    !vt=-v0t*exp(-kt*rr**2)
    !vs=-v0s*exp(-ks*rr**2)
    
    vr =  v0r * pi**1.5d0 * exp(-q2/(4.d0*kr) )/ (kr**1.5d0) 
    vt = -v0t * pi**1.5d0 * exp(-q2/(4.d0*kt) )/ (kt**1.5d0)
    vs = -v0s * pi**1.5d0 * exp(-q2/(4.d0*ks) )/ (ks**1.5d0)
    
    vcentral=0.25d0*(vr+vs)*delta(m1,m3)*delta(m2,m4)*delta(t1,t3)*delta(t2,t4) 
    vsigma=-0.25d0*(vr+vs)*chp_sigma_dot_sigma(m1,m2,m3,m4)*delta(t1,t3)*delta(t2,t4)
!!$


    vr = vr * (3.d0*delta(m1,m3)*delta(m2,m4)*delta(t1,t3)*delta(t2,t4) - & 
         chp_tau_dot_tau(t1,t2,t3,t4)*delta(m1,m3)*delta(m2,m4) - & 
         chp_sigma_dot_sigma(m1,m2,m3,m4)*delta(t1,t3)*delta(t2,t4) - & 
         chp_tau_dot_tau(t1,t2,t3,t4)*chp_sigma_dot_sigma(m1,m2,m3,m4) )/8.d0
    
    vs = vs * (3.d0*delta(m1,m3)*delta(m2,m4)*delta(t1,t3)*delta(t2,t4) + & 
         chp_tau_dot_tau(t1,t2,t3,t4)*delta(m1,m3)*delta(m2,m4) - & 
         3.d0*chp_sigma_dot_sigma(m1,m2,m3,m4)*delta(t1,t3)*delta(t2,t4) - & 
         chp_tau_dot_tau(t1,t2,t3,t4)*chp_sigma_dot_sigma(m1,m2,m3,m4) )/16.d0
    
    vt = vt * (3.d0*delta(m1,m3)*delta(m2,m4)*delta(t1,t3)*delta(t2,t4) - & 
         3.d0*chp_tau_dot_tau(t1,t2,t3,t4)*delta(m1,m3)*delta(m2,m4) + & 
         chp_sigma_dot_sigma(m1,m2,m3,m4)*delta(t1,t3)*delta(t2,t4) - & 
         chp_tau_dot_tau(t1,t2,t3,t4)*chp_sigma_dot_sigma(m1,m2,m3,m4) )/16.d0
    
!!$
    
    !vcentral(r_12)+vsigma(r_12)*sigma_1.sigma_2.
    
    vdir = (vs+vr+vt)
    !vdir = (vcentral + vsigma)
    
    !if ( t1 + t1 == 0 ) then
    !!   vdir = vdir
    !else 
    !   vdir = vdir / 2.d0
    !end if
    
    !
    ! Subtract-body part of CoM kinetic energy
    !
    vmom_minnesota = vdir/volume - sum(k1*k2)*delta(p,r) * delta(q,s) /p_mass/below_ef
    
    
  end function vmom_minnesota
  
  real*8 function vmom_minnesota_new(p,q,r,s)
    USE single_particle_orbits
    USE constants
    use chiral_constants
  
    implicit none 
    INTEGER :: p,q,r,s, m1,m2,m3,m4, spin, iph,i, t1,t2,t3,t4, Tiso
    INTEGER :: nx1, ny1, nz1, nx2, ny2, nz2, nx3, ny3, nz3, nx4, ny4, nz4 
    REAL*8 :: k1(3), k2(3), k3(3), k4(3), kmean(3) 
    REAL*8 :: qtrans(3), prel(3), pprel(3), ss ,tt, pref_i, plusST, minusST
    REAL*8 :: q2, p2, qabs, pp2, vdir, vexc, cg1, cg2 
    REAL*8 :: delta, nucleon_mass, relativity_factor
    REAL*8 :: v0r, v0t, v0s, kr,kt,ks, vr,vt,vcentral, vsigma, spin_exc1, spin_exc2
    real(8):: mu(3),VC(3),WC(3),VS(3),VCtild(3)
    real(8):: V_i(3),mu_i(3),W_i(3), M_i(3), B_i(3),H_i(3)
    
    
    
    nx1 = all_orbit%nx(p)
    ny1 = all_orbit%ny(p)
    nz1 = all_orbit%nz(p)
    nx2 = all_orbit%nx(q)
    ny2 = all_orbit%ny(q)
    nz2 = all_orbit%nz(q)
    nx3 = all_orbit%nx(r)
    ny3 = all_orbit%ny(r)
    nz3 = all_orbit%nz(r)
    nx4 = all_orbit%nx(s)
    ny4 = all_orbit%ny(s)
    nz4 = all_orbit%nz(s)
    ! 
    ! Conservation of linear momentum
    !
    if ( nx1 + nx2 /= nx3 + nx4 ) return 
    if ( ny1 + ny2 /= ny3 + ny4 ) return 
    if ( nz1 + nz2 /= nz3 + nz4 ) return 
  
    k1(1) = all_orbit%kx(p)
    k1(2) = all_orbit%ky(p)
    k1(3) = all_orbit%kz(p)
    k2(1) = all_orbit%kx(q)
    k2(2) = all_orbit%ky(q)
    k2(3) = all_orbit%kz(q)
    k3(1) = all_orbit%kx(r)
    k3(2) = all_orbit%ky(r)
    k3(3) = all_orbit%kz(r)
    k4(1) = all_orbit%kx(s)
    k4(2) = all_orbit%ky(s)
    k4(3) = all_orbit%kz(s)
  
  
    ! 
    ! conservation of spin and isospin 
    !
    if ( all_orbit%szp(p) + all_orbit%szp(q) /= all_orbit%szp(r) + all_orbit%szp(s) ) return 
    if ( all_orbit%itzp(p) + all_orbit%itzp(q) /= all_orbit%itzp(r) + all_orbit%itzp(s) ) return 
  
    m1 = all_orbit%szp(p) 
    m2 = all_orbit%szp(q) 
    m3 = all_orbit%szp(r) 
    m4 = all_orbit%szp(s) 
  
    t1 = all_orbit%itzp(p) 
    t2 = all_orbit%itzp(q) 
    t3 = all_orbit%itzp(r) 
    t4 = all_orbit%itzp(s) 
  
    ! 
    ! RELATIVE MOMENTA <prel |v| pprel > 
    ! 
    prel  = 0.5d0*(k1-k2)
    pprel = 0.5d0*(k3-k4)
    !
    ! momentum transfer 
    !
    qtrans = prel - pprel
    q2 = sum(qtrans*qtrans) 
    
    
    mu_i = (/1.487,0.639,0.465/)
    V_i = (/200.,-300.,-300./)
    W_i = (/.5,.25,.25/)
    M_i = (/.5,.25,.25/)
    B_i = (/0.,.25,-.25/)
    H_i = (/0.,.25,-.25/)

    DO i = 1, 3
       VC(i) = V_i(i)*(W_i(i)+ .5d0*(B_i(i)-H_i(i)))
       VS(i) = V_i(i)*B_i(i)/2.d0
       WC(i) = -V_i(i)*H_i(i)/2.d0
       VCtild(i) = V_i(i)*M_i(i)
    ENDDO
    
    
    ! r-space 
    !vr=v0r*exp(-kr*rr**2)
    !vt=-v0t*exp(-kt*rr**2)
    !vs=-v0s*exp(-ks*rr**2)
    
    tt = chp_tau_dot_tau(t1,t2,t3,t4)*delta(m1,m3)*delta(m2,m4)
    ss = chp_sigma_dot_sigma(m1,m2,m3,m4)*delta(t1,t3)*delta(t2,t4)
    vdir = 0.d0 
    DO i = 1, 3
       pref_i = 2./pi * dexp(-q2/4/mu_i(i))*(pi/mu_i(i))**1.5 /(8.*pi)
       !pref_i = dexp(-q2/4/mu_i(i))*(pi/mu_i(i))**1.5
       minusST = VCtild(i)*delta(t1,t3)*delta(t2,t4)*delta(m1,m3)*delta(m2,m4)
       
       plusST = VC(i)*delta(t1,t3)*delta(t2,t4)*delta(m1,m3)*delta(m2,m4) + WC(i)*tt + VS(i)*ss
       vdir = vdir + pref_i*(plusST+minusST)
    ENDDO
    
    vdir = vdir/volume
    
    !
    ! Subtract-body part of CoM kinetic energy
    !
    vmom_minnesota_new = vdir - sum(k1*k2)*delta(p,r) * delta(q,s) /p_mass/below_ef
    
    
  end function vmom_minnesota_new


end module chiral_potentials

  


real*8 function vmom_yukawa_rel(p,q)
  USE single_particle_orbits
  USE constants

  implicit none 
  INTEGER :: p,q,m1,m2,m3,m4, spin, iph, t1,t2,t3,t4
  INTEGER :: nx1, ny1, nz1, nx2, ny2, nz2, nx3, ny3, nz3, nx4, ny4, nz4 
  REAL*8 :: kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4 
  REAL*8 :: q2, vdir, vexc, V0, mu
  REAL*8 :: delta, va, vb, mu_a, mu_b

  ! v is in units MeVfm 
  Va = -4.7*hbarc ! -3.179d0 * hbarc  
  Vb =  7.291d0 * hbarc 
  ! mu is in units of fm-2
  mu_a = 305.86d0/hbarc
  mu_b = 613.69d0/hbarc 
  
  vmom_yukawa_rel = 0.d0
  kx1 = all_orbit%kx(p)
  ky1 = all_orbit%ky(p)
  kz1 = all_orbit%kz(p)
  kx2 = all_orbit%kx(q)
  ky2 = all_orbit%ky(q)
  kz2 = all_orbit%kz(q)
  
  
  q2 = ( (kx1-kx2)**2 + (ky1-ky2)**2 + (kz1-kz2)**2) 
  vdir = (va/ ( q2 + mu_a**2) +vb/ ( q2 + mu_b**2) )
  
  vmom_yukawa_rel = 4.*pi*(vdir)/(volume) 
  

end function vmom_yukawa_rel





real*8 function vrspace_minnesota(rr)
  USE single_particle_orbits
  USE constants

  implicit none 
  REAL*8 :: rr  
  REAL*8 :: v0r, v0t, v0s, kr,kt,ks, vr,vt,vs,vcentral, vsigma
  REAL*8 :: delta, va, vb, mu_a, mu_b
  
  v0r=200.0  ! MeV
  v0t=178.0  ! MeV
  v0s=91.85  ! MeV
  kr=1.487  ! fm**-2
  kt=0.639  ! fm**-2
  ks=0.465  ! fm**-2
  vr=v0r*exp(-kr*rr**2)
  vt=-v0t*exp(-kt*rr**2)
  vs=-v0s*exp(-ks*rr**2)
  vcentral=0.25d0*(vr+vs)
  vsigma=-vcentral

  vrspace_minnesota = vcentral 

end function vrspace_minnesota



real*8 function vrspace_yukawa(r)
  USE single_particle_orbits
  USE constants

  implicit none 
  REAL*8 :: r  
  REAL*8 :: q2, vdir, vexc, V0, mu
  REAL*8 :: delta, va, vb, mu_a, mu_b

  
  ! v is in units MeVfm 
  Va = -3.179d0 * hbarc  
  Vb =  7.291d0 * hbarc 
  ! mu is in units of fm-2
  mu_a = 305.86d0/hbarc
  mu_b = 613.69d0/hbarc 
  
  vrspace_yukawa = va*exp(-mu_a*r)/r + vb*exp(-mu_b*r)/r 
  
end function vrspace_yukawa


real*8 function vmom_yukawa(p,q,r,s)
  USE single_particle_orbits
  USE constants

  implicit none 
  INTEGER :: p,q,r,s, m1,m2,m3,m4, spin, iph, t1,t2,t3,t4
  INTEGER :: nx1, ny1, nz1, nx2, ny2, nz2, nx3, ny3, nz3, nx4, ny4, nz4 
  REAL*8 :: kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4 
  REAL*8 :: q2, vdir, vexc, V0, mu
  REAL*8 :: delta, va, vb, mu_a, mu_b
  
  ! v is in units MeVfm 
  Va = -4.7*hbarc ! -3.179d0 * hbarc  
  Vb =  7.291d0 * hbarc 
  ! mu is in units of fm-2
  mu_a = 305.86d0/hbarc
  mu_b = 613.69d0/hbarc 
  
  vmom_yukawa = 0.d0
  nx1 = all_orbit%nx(p)
  ny1 = all_orbit%ny(p)
  nz1 = all_orbit%nz(p)
  nx2 = all_orbit%nx(q)
  ny2 = all_orbit%ny(q)
  nz2 = all_orbit%nz(q)
  nx3 = all_orbit%nx(r)
  ny3 = all_orbit%ny(r)
  nz3 = all_orbit%nz(r)
  nx4 = all_orbit%nx(s)
  ny4 = all_orbit%ny(s)
  nz4 = all_orbit%nz(s)
  ! 
  ! Conservation of linear momentum
  !
  if ( nx1 + nx2 /= nx3 + nx4 ) return 
  if ( ny1 + ny2 /= ny3 + ny4 ) return 
  if ( nz1 + nz2 /= nz3 + nz4 ) return 

  kx1 = all_orbit%kx(p)
  ky1 = all_orbit%ky(p)
  kz1 = all_orbit%kz(p)
  kx2 = all_orbit%kx(q)
  ky2 = all_orbit%ky(q)
  kz2 = all_orbit%kz(q)
  kx3 = all_orbit%kx(r)
  ky3 = all_orbit%ky(r)
  kz3 = all_orbit%kz(r)
  kx4 = all_orbit%kx(s)
  ky4 = all_orbit%ky(s)
  kz4 = all_orbit%kz(s)

  ! 
  ! conservation of spin and isospin 
  !
  if ( all_orbit%szp(p) + all_orbit%szp(q) /= all_orbit%szp(r) + all_orbit%szp(s) ) return 
  if ( all_orbit%itzp(p) + all_orbit%itzp(q) /= all_orbit%itzp(r) + all_orbit%itzp(s) ) return 
  
  m1 = all_orbit%szp(p) 
  m2 = all_orbit%szp(q) 
  m3 = all_orbit%szp(r) 
  m4 = all_orbit%szp(s) 
  
  t1 = all_orbit%itzp(p) 
  t2 = all_orbit%itzp(q) 
  t3 = all_orbit%itzp(r) 
  t4 = all_orbit%itzp(s) 
  
  
  q2 = 0.25d0*( (kx1-kx2-kx3+kx4)**2 + (ky1-ky2-ky3+ky4)**2 + (kz1-kz2-kz3+kz4)**2) 
  vdir = (va/ ( q2 + mu_a**2) +vb/ ( q2 + mu_b**2) )*delta(m1,m3)*delta(m2,m4)*delta(t1,t3)*delta(t2,t4)
  
  q2 = 0.25d0*( (kx1-kx2-kx4+kx3)**2 + (ky1-ky2-ky4+ky3)**2 + (kz1-kz2-kz4+kz3)**2)
  vexc =  (va/ ( q2 + mu_a**2) +vb/ ( q2 + mu_b**2) ) *delta(m1,m4)*delta(m2,m3)*delta(t1,t4)*delta(t2,t3)
  
  
  vmom_yukawa = 4.*pi*(vdir - vexc)/volume 
  
    
  vdir = ( kx1*kx2 + ky1*ky2 + kz1*kz2 )*delta(p,r) * delta(q,s) 
  vexc = ( kx1*kx2 + ky1*ky2 + kz1*kz2 )*delta(p,s) * delta(q,r) 
  
  vmom_yukawa = vmom_yukawa - (vdir - vexc)*hbarc**2/p_mass/below_ef
  

end function vmom_yukawa


real*8 function kin_energy(p,q,r,s)
  USE single_particle_orbits
  USE constants
  use one_body_operators
  use t2_storage
  use configurations
  USE KSPACE 

  implicit none 
  INTEGER :: p,q,r,s, m1,m2,m3,m4, spin, iph
  INTEGER :: nx1, ny1, nz1, nx2, ny2, nz2, nx3, ny3, nz3, nx4, ny4, nz4 
  REAL*8 :: q2, vdir, vexc, V0, mu
  REAL*8 :: delta
  
  
  kin_energy = 0.d0
  nx1 = all_orbit%nx(p)
  ny1 = all_orbit%ny(p)
  nz1 = all_orbit%nz(p)
  nx2 = all_orbit%nx(q)
  ny2 = all_orbit%ny(q)
  nz2 = all_orbit%nz(q)
  nx3 = all_orbit%nx(r)
  ny3 = all_orbit%ny(r)
  nz3 = all_orbit%nz(r)
  nx4 = all_orbit%nx(s)
  ny4 = all_orbit%ny(s)
  nz4 = all_orbit%nz(s)
  
  ! 
  ! Conservation of linear momentum
  !
  if ( nx1 + nx2 /= nx3 + nx4 ) return 
  if ( ny1 + ny2 /= ny3 + ny4 ) return 
  if ( nz1 + nz2 /= nz3 + nz4 ) return 
  ! 
  ! conservation of spin and isospin 
  !
  if ( all_orbit%szp(p) + all_orbit%szp(q) /= all_orbit%szp(r) + all_orbit%szp(s) ) return 
  if ( all_orbit%itzp(p) + all_orbit%itzp(q) /= all_orbit%itzp(r) + all_orbit%itzp(s) ) return 
  
  m1 = all_orbit%szp(p) 
  m2 = all_orbit%szp(q) 
  m3 = all_orbit%szp(r) 
  m4 = all_orbit%szp(s) 
  
  
  vdir = tkin(p,r)*delta(q,s)+tkin(q,s)*delta(p,r) 
  vexc = tkin(p,s)*delta(q,r)+tkin(q,r)*delta(p,s)

  kin_energy = vdir - vexc
  
end function kin_energy

