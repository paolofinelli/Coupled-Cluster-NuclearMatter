module chiral_tables 
  complex*16, allocatable :: sigma_dot_q_tab(:,:,:,:,:), sigma_dot_qxq_tab(:,:,:,:,:)
  complex*16, allocatable :: sigma_dot_q_ope_tab(:,:,:,:,:)
  complex*16 :: tau1_dot_tauXtau_tab(-1:1,-1:1,-1:1,-1:1,-1:1,-1:1)
  complex*16 :: tau_dot_tau_tab(-1:1,-1:1,-1:1,-1:1)
  real*8 :: delta_tab(-1:1,-1:1) 
  
end module chiral_tables

module chiral_constants 
  use constants, only : pi  
  
  ! basic mesh info
  TYPE, PUBLIC :: chp_mesh_info
     INTEGER  :: amount
     REAL(8) :: xmin, xmax
  END TYPE chp_mesh_info
  
  ! GAUSS-LEGENDRE MESH POINTS
  ! Gauss-Legendre mesh point x, corresponding integration weight w and corresponding x*x*w-value
  TYPE, PUBLIC :: chp_gauleg_mesh_point
     SEQUENCE
     REAL(8) :: x, w, xxw
  END TYPE chp_gauleg_mesh_point
  
  ! mesh points and weights in momentum space
  TYPE, PUBLIC :: chp_gauleg_mesh
     TYPE(chp_mesh_info)                                    :: info
     TYPE(chp_gauleg_mesh_point), DIMENSION(:), ALLOCATABLE :: pnt
  END TYPE chp_gauleg_mesh
  
  ! chiral order definition
  INTEGER, PARAMETER :: LO   = 0
  ! all contributions vanish at order 1
  ! due to parity and time-reversal  invariance
  INTEGER, PARAMETER :: NLO  = 2
  INTEGER, PARAMETER :: NNLO = 3
  INTEGER, PARAMETER :: N3LO = 4
  
  real*8, parameter :: gA =  1.29000000000000004D+00  !1.29d0 
  real*8, parameter :: gA2 = gA*gA 
  real*8, parameter :: gA4 = gA2*gA2 
  real*8, parameter :: fpi = 9.24000000000000057D+01 !92.4d0 
  real*8, parameter :: fpi2 = fpi*fpi
  real*8, parameter :: fpi4 = fpi2*fpi2
  real*8, parameter :: fpi_inv = 1.d0/fpi
  real*8, parameter :: twopi = pi*2.d0
  real*8, parameter :: pi2 = pi*pi
  real*8 :: const(100) 
  real*8, parameter :: proton_mass  = 9.38272046000000046D+02  !938.272d0 
  real*8, parameter :: neutron_mass = 9.39565379000000007D+02  !939.5653d0
  real*8 :: mnuc(-1:1), mpi(-1:2)
  real*8 :: mnuc_inv(-1:1)! nucleon mass inversed
  REAL*8 :: mnuc2(-1:1)   ! nucleon mass squared
  REAL*8 :: mpi2(-1:2)    ! pion mass squared
  REAL*8 :: mpi3(-1:2)    ! pion mass cubed
  REAL*8 :: mpi4(-1:2)    ! mpi2 squared
  REAL*8 :: mpi5(-1:2)    ! mpi^5
  REAL*8 :: twompi(-1:2)  ! two times pion mass
  REAL*8 :: fourmpi2(-1:2)! four times pion mass squared
  COMPLEX*16 :: sigma_x(2,2), sigma_y(2,2), sigma_z(2,2)
  REAL*8 :: c1,c2,c3,c4
  REAL*8 :: C1_3NF, C2_3NF, C3_3NF, C4_3NF
  REAL*8, PARAMETER :: lambda = 450.D0 
  REAL*8, PARAMETER :: sfr = 700.d0 
  REAL*8, PARAMETER :: sfr2 = sfr*sfr
  CHARACTER(LEN=2), PARAMETER :: chp_itope = 'EM'
  REAL*8 :: sfr_heavyside(-1:2) ! THETA(sfr-twompi)
  REAL*8 :: CS(-1:1), CT(-1:1), cnlo(1:7), cn3lo(1:15) 
  REAL*8 :: CE, CD, Econst, Dconst, lec_c1, lec_c2, lec_c3, lec_c4   
  REAL*8 :: lec_c1_3NF, lec_c2_3NF, lec_c3_3NF, lec_c4_3NF, z_low   

! x and z are integration variables. Some of the expressions for the
! 2PE 2-loop diagrams contains integrals that must be solved numerically,
! The integration variables are x and z = 2*m_\pi / \mu using the notation
! from [1] and [2]. The number of integration points used here are
! enough to get highly accurate results.
! The integral over x is independent of q, and thus need only be done once
! for each system (np, pp, nn).
  INTEGER,PARAMETER :: nof_2PE_2loop_int_x_points = 50
  INTEGER,PARAMETER :: nof_2PE_2loop_int_z_points = 30
  ! variables for calculating optional 2PE 2loop diagrams
  LOGICAL :: chp_2PE_2loop_int_VTS_data_set
  LOGICAL :: chp_2PE_2loop_int_WC_DR_data_set
  LOGICAL :: chp_2PE_2loop_int_WC_SFR_data_set
  ! Integration meshs for non-analytical integrals in some of the 2PE 2-loop diagrams
  TYPE(chp_gauleg_mesh) :: chp_2PE_2loop_int_mesh_x
  TYPE(chp_gauleg_mesh) :: chp_2PE_2loop_int_mesh_z
  TYPE(chp_gauleg_mesh) :: chp_2PE_2loop_int_mesh_z_VC_SFR
  ! Factors in the x-integrations that need only be calculated once
  REAL(8) :: chp_2PE_2loop_int_VTS_x_fact(nof_2PE_2loop_int_z_points)
  REAL(8) :: chp_2PE_2loop_int_WC_DR_x_fact(nof_2PE_2loop_int_z_points)
  REAL(8) :: chp_2PE_2loop_int_WC_SFR_x_fact(nof_2PE_2loop_int_z_points)
  
  
  !! Value taken from Andreas
  !! These are the scale-independent LECs
  
  real*8,PARAMETER :: d1_plus_d2    = 6.22349655930370016D+00*1.0D-6
  real*8,PARAMETER :: d3            =-5.30654650657754967D+00*1.0D-6
  real*8,PARAMETER :: d5            =-4.63631448771476995D-01*1.0D-6
  real*8,PARAMETER :: d14_minus_d15 =-1.09953639023018006D+01*1.0D-6
 
  
CONTAINS 
  
  subroutine init_chp_constants 
    use parallel 
    implicit none   
    double precision :: c1s0(-1:1), c3s1(-1:1), cnlo_pw(1:7), cn3lo_pw(1:15)
    real*8 :: lambdachi
    integer :: i, row
    REAL*8 :: N3LOPW2ST(15,15), N3LOST2PW(15,15)
    
    mnuc(-1) = proton_mass
    mnuc(0)  = 0.5d0*(proton_mass+neutron_mass)
    mnuc(+1) = neutron_mass
    
    mnuc_inv = 1.D0 / mnuc
    
    mpi(-1)  = 1.39570179999999993D+02 !139.5702d0
    mpi(0)   = 1.34976599999999991D+02 !134.9766d0
    mpi(+1)  = 1.39570179999999993D+02 !139.5702d0
    mpi(+2)  = ( mpi(-1)+mpi(0)+mpi(+1) ) /3.d0  
    
    
    mnuc2(:)   = mnuc*mnuc 
    mpi2(:)    = mpi*mpi   
    mpi3(:)    = mpi*mpi2
    mpi4(:)    = mpi2*mpi2
    mpi5(:)    = mpi4*mpi 
    twompi(:)  = 2.0D0*mpi 
    fourmpi2(:)= 4.0D0*mpi2
    
    chp_2PE_2loop_int_VTS_data_set = .FALSE.
    chp_2PE_2loop_int_WC_DR_data_set  = .FALSE.
    chp_2PE_2loop_int_WC_SFR_data_set  = .FALSE.
    
    
!!$    !
!!$    ! N2LO_opt(500) LECs 
!!$    ! 
!!$    LEC_c1 =  -0.9186395287347203d0
!!$    LEC_c2 =   0.d0 
!!$    LEC_c3 =  -3.888687492763241d0
!!$    LEC_c4 =   4.310327160829740d0
!!$    LEC_c1_3NF =  LEC_c1 !-0.9186395287347203d0
!!$    LEC_c2_3NF =   0.d0 
!!$    LEC_c3_3NF =  LEC_c3 !-3.888687492763241d0
!!$    LEC_c4_3NF =  LEC_c4 ! 4.310327160829740d0
!!$    
!!$    !
!!$    ! leading order contacts in PW
!!$    !
!!$    c1s0(-1) = -0.1513660372031080D+00
!!$    c1s0(0)  = -0.1521410882366787D+00
!!$    c1s0(1)  = -0.1517647459006913D+00
!!$    
!!$    c3s1(-1) = -0.1584341766228121D+00
!!$    c3s1(0)  = -0.1584341766228121D+00
!!$    c3s1(1)  = -0.1584341766228121D+00
!!$    
!!$    !
!!$    ! next-to-leading order contacts in PW
!!$    !
!!$    cnlo_pw(1) =  0.2404021944134705D+01
!!$    cnlo_pw(2) =  0.1263390763475578D+01
!!$    cnlo_pw(3) =  0.4170455420556486D+00
!!$    cnlo_pw(4) = -0.7826584999752046D+00
!!$    cnlo_pw(5) =  0.9283846626623043D+00
!!$    cnlo_pw(6) =  0.6181414190474576D+00
!!$    cnlo_pw(7) = -0.6778085114063558D+00
!!$
!!$    
    
    
    
    LEC_c1 =  -6.85008592144712991D-01  !-1.399355229598470D0    
    LEC_c2 =   2.99405406565804988D+00   !1.736403902935290D0    
    LEC_c3 =  -4.11681871693645984D+00  !-4.608490948168820D0    
    LEC_c4 =   5.35080249059179014D+00   ! 3.704035249000000D0    
    LEC_c1_3NF =  LEC_c1 !-0.9186395287347203d0
    LEC_c2_3NF =   0.d0 
    LEC_c3_3NF =  LEC_c3 !-3.888687492763241d0
    LEC_c4_3NF =  LEC_c4 ! 4.310327160829740d0
    
    
    !
    ! leading order contacts in PW
    !
    c1s0(-1) = -1.54345070623573999D-01 ! Ct_1S0pp   -0.158407847000000D0  OG
    c1s0(0)  = -1.55450952594823999D-01 ! Ct_1S0np   -0.160117163000000D0  
    c1s0(1)  = -1.55002481436744999D-01 ! Ct_1S0nn   -0.159470077000000D0      
    c3s1(-1) = -1.61531977709016006D-01 ! Ct_3S1pp  -0.180121483808290D0
    c3s1(0)  = -1.61531977709016006D-01 ! Ct_3S1np  -0.180121483808290D0
    c3s1(1)  = -1.61531977709016006D-01 ! Ct_3S1nn  -0.180121483808290D0  
    
    
    !
    ! next-to-leading order contacts in PW
    !
    cnlo_pw(1) = 2.43124677048033000D+00 ! C_1S0        2.039737161000000D0 ! C_1S0    OG
    cnlo_pw(2) = 1.25740974180182996D+00 ! C_3P0        1.112798744000000D0 ! C_3P0       
    cnlo_pw(3) = 3.79234207699529979D-01 ! C_1P1        0.065408409052280D0 ! C_1P1       
    cnlo_pw(4) = -8.10456470917395921D-01 ! C_3P1      -0.611643165572000D0 ! C_3P1       
    cnlo_pw(5) = 6.75058823794391971D-01 ! C_3S1       1.123342269000000D0 ! C_3S1       
    cnlo_pw(6) = 6.79668131973288014D-01 ! C_3S1-3D1   0.570361534000000D0 ! C_3S1-3D1   
    cnlo_pw(7) = -6.40320460232033040D-01 ! C_3P2      -0.616100345000000D0 ! C_3P2       
    

    cn3lo_pw(1) =  -6.80955215444071982D+00 ! Dh_1S0      -9.053098975549741D0 ! Dh_1S0       OG
    cn3lo_pw(2) =  -2.83771264573585000D+01 ! D_1S0      -41.174147949508580D0 ! D_1S0        
    cn3lo_pw(3) =   3.18984094348776015D+00 ! D_3P0         0.975566993000000D0 ! D_3P0        
    cn3lo_pw(4) =   1.42720692941284000D+00 ! D_1P1         4.992115255686630D0 ! D_1P1        
    cn3lo_pw(5) =   7.67533624169093009D+00 ! D_3P1         9.662526669000000D0 ! D_3P1        
    cn3lo_pw(6) =   3.85425727120557005D+00 ! Dh_3S1        8.241947272000001D0 ! Dh_3S1       
    cn3lo_pw(7) =  -3.26672846982458012D+01 ! D_3S1       -5.530061984562440D0 ! D_3S1        
    cn3lo_pw(8) =  -2.43147613205035018D+00 ! D_3D1        0.696673967313840D0 ! D_3D1        
    cn3lo_pw(9) =   4.03300295736832015D-01 ! Dh_3S1-3D1    6.208007914000000D0 ! Dh_3S1-3D1   
    cn3lo_pw(10) = -5.37768260640681994D+00 ! D_3S1-3D1    0.690083154537830D0 ! D_3S1-3D1    
    cn3lo_pw(11) = -4.39661753305750036D+00 ! D_1D2       -5.458404347000000D0 ! D_1D2        
    cn3lo_pw(12) = -1.72535934926913015D+00 ! D_3D2       -1.005291586000000D0 ! D_3D2        
    cn3lo_pw(13) =  6.72814338560019998D+00 ! D_3P2         7.707987498161421D0 ! D_3P2        
    cn3lo_pw(14) = -5.56020678597204054D-01 ! D_3P2-3F2   -1.189722077358910D0 ! D_3P2-3F2    
    cn3lo_pw(15) = -2.89647473116136034D+00 ! D_3D3        2.997513967341410D0 ! D_3D3        
    





!!$    !
!!$    ! N3LO_EM(500) LECs  (N3LO Entem&Machleidt)
!!$    ! 
!!$    LEC_c1 =  -0.810000000000000D0
!!$    LEC_c2 =   2.800000000000000D0
!!$    LEC_c3 =  -3.200000000000000D0
!!$    LEC_c4 =   5.400000000000000D0
!!$    LEC_c1_3NF =  LEC_c1 !-0.9186395287347203d0
!!$    LEC_c2_3NF =   0.d0 
!!$    LEC_c3_3NF =  LEC_c3 !-3.888687492763241d0
!!$    LEC_c4_3NF =  LEC_c4 ! 4.310327160829740d0
!!$    
!!$    
!!$
!!$    !
!!$    ! leading order contacts in PW
!!$    !
!!$    c1s0(-1) = -0.145286000000000D0 
!!$    c1s0(0)  = -0.147167000000000D0 
!!$    c1s0(1)  = -0.146285000000000D0 
!!$    
!!$    c3s1(-1) = -0.118972496000000D0 
!!$    c3s1(0)  = -0.118972496000000D0 
!!$    c3s1(1)  = -0.118972496000000D0 
!!$    
!!$    
!!$    !
!!$    ! next-to-leading order contacts in PW
!!$    !
!!$    cnlo_pw(1) =  2.380000000000000D0 ! C_1S0    
!!$    cnlo_pw(2) =  1.487000000000000D0 ! C_3P0
!!$    cnlo_pw(3) =  0.656000000000000D0 ! C_1P1
!!$    cnlo_pw(4) = -0.630000000000000D0 ! C_3P1
!!$    cnlo_pw(5) =  0.760000000000000D0 ! C_3S1
!!$    cnlo_pw(6) =  0.826000000000000D0 ! C_3S1-3D1
!!$    cnlo_pw(7) = -0.538000000000000D0 ! C_3P2
!!$    
!!$    !
!!$    ! N3LO contacts in PW
!!$    !
!!$    cn3lo_pw(1) =     -2.545000000000000D0 ! Dh_1S0
!!$    cn3lo_pw(2) =    -16.000000000000000D0 ! D_1S0
!!$    cn3lo_pw(3) =      0.245000000000000D0 ! D_3P0
!!$    cn3lo_pw(4) =      5.250000000000000D0 ! D_1P1
!!$    cn3lo_pw(5) =      2.350000000000000D0 ! D_3P1
!!$    cn3lo_pw(6) =      7.000000000000001D0 ! Dh_3S1
!!$    cn3lo_pw(7) =      6.550000000000000D0 ! D_3S1
!!$    cn3lo_pw(8) =     -2.800000000000000D0 ! D_3D1
!!$    cn3lo_pw(9) =      2.250000000000000D0 ! Dh_3S1-3D1
!!$    cn3lo_pw(10) =     6.610000000000000D0 ! D_3S1-3D1
!!$    cn3lo_pw(11) =    -1.770000000000000D0 ! D_1D2
!!$    cn3lo_pw(12) =    -1.460000000000000D0 ! D_3D2
!!$    cn3lo_pw(13) =     2.295000000000000D0 ! D_3P2
!!$    cn3lo_pw(14) =    -0.465000000000000D0 ! D_3P2-3F2
!!$    cn3lo_pw(15) =     5.660000000000000D0 ! D_3D3
!!$    
    

!!$    LEC_c1 =  -1.121521199632590D0
!!$    LEC_c2 =   0.d0
!!$    LEC_c3 =  -3.925005856486820D0
!!$    LEC_c4 =   3.765687158585920D0
!!$    LEC_c1_3NF =  LEC_c1 !-0.9186395287347203d0
!!$    LEC_c2_3NF =   0.d0 
!!$    LEC_c3_3NF =  LEC_c3 !-3.888687492763241d0
!!$    LEC_c4_3NF =  LEC_c4 ! 4.310327160829740d0
!!$    
!!$
!!$    !
!!$    ! leading order contacts in PW
!!$    !
!!$    c1s0(-1) = -0.158149379370110D0 ! Ct_1S0pp
!!$    c1s0(0)  = -0.159822449578320D0 ! Ct_1S0np
!!$    c1s0(1)  = -0.159150268280180D0 ! Ct_1S0nn
!!$    
!!$    c3s1(-1) = -0.177674364499000D0 ! Ct_3S1pp
!!$    c3s1(0)  = -0.177674364499000D0 ! Ct_3S1np
!!$    c3s1(1)  = -0.177674364499000D0 ! Ct_3S1nn
!!$    
!!$    !
!!$    ! next-to-leading order contacts in PW
!!$    !
!!$    cnlo_pw(1) = 2.539367785050380D0 
!!$    cnlo_pw(2) = 1.398365591876140D0 
!!$    cnlo_pw(3) = 0.555958765133350D0 
!!$    cnlo_pw(4) =-1.136095263327820D0 
!!$    cnlo_pw(5) = 1.002892673483510D0
!!$    cnlo_pw(6) = 0.600716048335960D0
!!$    cnlo_pw(7) =-0.802300295338460D0
    
    
    !
    ! Get right units of LECs 
    !
    c1 = LEC_c1*1.0D-3
    c2 = LEC_c2*1.0D-3
    c3 = LEC_c3*1.0D-3
    c4 = LEC_c4*1.0D-3

    c1_3NF = LEC_c1_3NF*1.0D-3
    c2_3NF = LEC_c2_3NF*1.0D-3
    c3_3NF = LEC_c3_3NF*1.0D-3
    c4_3NF = LEC_c4_3NF*1.0D-3
    
    
    ! NLO
    ! the NLO contacts are input in units of 10^4/GeV^4
    ! this = 10^-8/MeV^4
    DO i=1,7
       CNLO_PW(i) = CNLO_PW(i) * 1.D-08
    END DO
    
    ! the LO CIBcontacts are input in units of 10^4/GeV^2
    ! this = 10^-2/MeV^2 
    c1s0 = c1s0*0.01d0 
    c3s1 = c3s1*0.01d0 
    
    ! See Eq. 4.39 p. 26
    CS = (c1s0 + 3.d0*c3s1) /16.d0/pi
    CT = (c3s1 - c1s0) /16.d0/pi
    
    
    
    ! See Eq. 4.42 p. 26 
    cnlo(1) = (-5.d0*cnlo_pw(7)+6.d0*cnlo_pw(5)-3.d0*cnlo_pw(4)-3.d0*cnlo_pw(3)-cnlo_pw(2)+2.d0*cnlo_pw(1))/(64.d0*pi)
    cnlo(2) = ( 5.d0*cnlo_pw(7)+6.d0*cnlo_pw(5)+3.d0*cnlo_pw(4)+3.d0*cnlo_pw(3)+cnlo_pw(2)+2.d0*cnlo_pw(1))/(16.d0*pi)
    cnlo(3) = -( 2.d0*cnlo_pw(7)-2.d0*dsqrt(2.d0)*cnlo_pw(6)-2.d0*cnlo_pw(5)-3.d0*cnlo_pw(3)+cnlo_pw(2)+2.d0*cnlo_pw(1))/(64.d0*pi)
    cnlo(4) = -(-2.d0*cnlo_pw(7)-2.d0*dsqrt(2.d0)*cnlo_pw(6)-2.d0*cnlo_pw(5)+3.d0*cnlo_pw(3)-cnlo_pw(2)+2.d0*cnlo_pw(1))/(16.d0*pi)
    cnlo(5) = -(-5.d0*cnlo_pw(7)+3.d0*cnlo_pw(4)+2.d0*cnlo_pw(2))/(16.d0*pi)
    cnlo(6) = ( cnlo_pw(7)-6.d0*dsqrt(2.d0)*cnlo_pw(6)-3.d0*cnlo_pw(4)+2.d0*cnlo_pw(2))/(64.d0*pi)
    cnlo(7) = -(cnlo_pw(7)+6.d0*dsqrt(2.d0)*cnlo_pw(6)-3.d0*cnlo_pw(4)+2.d0*cnlo_pw(2))/(16.d0*pi) 
    
    
    if ( iam == 0 ) WRITE(6,"(A12,F30.16)") 'C1', cnlo(1)* 1.D8
    if ( iam == 0 ) WRITE(6,"(A12,F30.16)") 'C2', cnlo(2)* 1.D8
    if ( iam == 0 ) WRITE(6,"(A12,F30.16)") 'C3', cnlo(3)* 1.D8
    if ( iam == 0 ) WRITE(6,"(A12,F30.16)") 'C4', cnlo(4)* 1.D8
    if ( iam == 0 ) WRITE(6,"(A12,F30.16)") 'C5', cnlo(5)* 1.D8
    if ( iam == 0 ) WRITE(6,"(A12,F30.16)") 'C6', cnlo(6)* 1.D8
    if ( iam == 0 ) WRITE(6,"(A12,F30.16)") 'C7', cnlo(7)* 1.D8
    

    !
    ! Add N3LO contact parameters 
    ! See Eq. E.2 p. 69
    ! 
    ! @N3LO       PW
    ! DN3LO( 1) : D1S0t     D1
    ! DN3LO( 2) : D1S0      D2
    ! DN3LO( 3) : D3P0      D3
    ! DN3LO( 4) : D1P1      D4
    ! DN3LO( 5) : D3P1      D5
    ! DN3LO( 6) : D3S1t     D6
    ! DN3LO( 7) : D3S1      D7
    ! DN3LO( 8) : D3D1      D8
    ! DN3LO( 9) : D3S1-3D1t D9
    ! DN3LO(10) : D3S1-3D1  D10
    ! DN3LO(11) : D1D2      D11
    ! DN3LO(12) : D3D2      D12
    ! DN3LO(13) : D3P2      D13
    ! DN3LO(14) : D3P2-3F2  D14
    ! DN3LO(15) : D3D3      D15
    
    N3LOST2PW = transpose(reshape((/ &
         1D0, 1D0/16D0, 0.25D0, 0D0, -3D0, -3D0/16D0, -0.75D0, 0D0, 0D0, 0D0, -1D0, -0.25D0, -0.25D0, -1D0/16D0, 0D0, &
         10D0/3D0, 5D0/24D0, 1D0/6D0, 2D0/3D0, -10D0, -0.625D0, -0.5D0, -2D0, 0D0, 0D0, -10D0/3D0, -1D0/6D0, -1D0/6D0, -5D0/24D0, -2D0/3D0, &
         -4D0/3D0, 1D0/12D0, 0D0, 0D0, -4D0/3D0, 1D0/12D0, 0D0, 0D0, -2D0/3D0, -1D0/6D0, 8D0/3D0, 1D0/3D0, -1D0/3D0, -1D0/6D0, 0D0, &
         -4D0/3D0, 1D0/12D0, 0D0, 0D0, 4D0, -0.25D0, 0D0, 0D0, 0D0, 0D0, 4D0/3D0, 0D0, 0D0, -1D0/12D0, 0D0, &
         -4D0/3D0, 1D0/12D0, 0D0, 0D0, -4D0/3D0, 1D0/12D0, 0D0, 0D0, -1D0/3D0, -1D0/12D0, -2D0, -1D0/6D0, 1D0/6D0, 1D0/8D0, 0D0, &
         1D0, 1D0/16D0, 0.25D0, 0D0, 1D0, 1D0/16D0, 0.25D0, 0D0, 0D0, 0D0, 1D0/3D0, 1D0/12D0, 1D0/12D0, 1D0/48D0, 0D0, &
         10D0/3D0, 5D0/24D0, 1D0/6D0, 2D0/3D0, 10D0/3D0, 5D0/24D0, 1D0/6D0, 2D0/3D0, 0D0, 0D0, 10D0/9D0, 1D0/18D0, 1D0/18D0, 5D0/72D0, 2D0/9D0, &
         8D0/15D0, 1D0/30D0, -2D0/15D0, -2D0/15D0, 8D0/15D0, 1D0/30D0, -2D0/15D0, -2D0/15D0, 2D0/5D0, -0.1D0, -4D0/9D0, 1D0/9D0, 1D0/9D0, -1D0/36D0, -16D0/45D0, &
         0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, -DSQRT(2D0)*2D0/3D0, -DSQRT(2D0)/6D0, -DSQRT(2D0)/6D0, -DSQRT(2D0)/24D0, 0D0, &
         0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, -DSQRT(2D0)*14D0/9D0, DSQRT(2D0)/18D0, DSQRT(2D0)/18D0, -DSQRT(2D0)*7D0/72D0, DSQRT(2D0)*2D0/9D0, &
         8D0/15D0, 1D0/30D0, -2D0/15D0, -2D0/15D0, -8D0/5D0, -0.1D0, 0.4D0, 0.4D0, 0D0, 0D0, -8D0/15D0, 2D0/15D0, 2D0/15D0, -1D0/30D0, 2D0/15D0, &
         8D0/15D0, 1D0/30D0, -2D0/15D0, -2D0/15D0, 8D0/15D0, 1D0/30D0, -2D0/15D0, -2D0/15D0, 2D0/15D0, -1D0/30D0, 0.8D0, -0.2D0, -0.2D0, 0.05D0, 4D0/15D0, &
         -4D0/3D0, 1D0/12D0, 0D0, 0D0, -4D0/3D0, 1D0/12D0, 0D0, 0D0, 1D0/3D0, 1D0/12D0, -2D0/15D0, 1D0/30D0, -1D0/30D0, 1D0/120D0, 0D0, &
         0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, DSQRT(6D0)*4D0/15D0, -DSQRT(6D0)/15D0, DSQRT(6D0)/15D0, -DSQRT(6D0)/60D0, 0D0, &
         8D0/15D0, 1D0/30D0, -2D0/15D0, -2D0/15D0, 8D0/15D0, 1D0/30D0, -2D0/15D0, -2D0/15D0, -4D0/15D0, 1D0/15D0, 0D0, 0D0, 0D0, 0D0, -2D0/15D0 &
         /), shape(N3LOST2PW)))
    N3LOST2PW = 4.d0*pi * N3LOST2PW
    N3LOPW2ST = N3LOST2PW
    CALL chp_matinv(N3LOPW2ST, 15)
    
    ! the N3LO contacts are input in units of 10^4/GeV^6
    ! this = 10^(4-3*6)/MeV^6 = 10^-14/MeV^6
    cN3LO_pw  = cN3LO_pw * 1.D-14
    
    DO row=1,15
       cn3lo(row) = SUM(N3LOPW2ST(row,:)*cn3lo_pw(:))
    END DO
    
    
    
    if ( iam == 0 ) WRITE(6,"(A12,F30.16)")   'D1' , cn3lo( 1)* 1.D14
    if ( iam == 0 ) WRITE(6,"(A12,F30.16)")   'D2' , cn3lo( 2)* 1.D14
    if ( iam == 0 ) WRITE(6,"(A12,F30.16)")   'D3' , cn3lo( 3)* 1.D14
    if ( iam == 0 ) WRITE(6,"(A12,F30.16)")   'D4' , cn3lo( 4)* 1.D14
    if ( iam == 0 ) WRITE(6,"(A12,F30.16)")   'D5' , cn3lo( 5)* 1.D14
    if ( iam == 0 ) WRITE(6,"(A12,F30.16)")   'D6' , cn3lo( 6)* 1.D14
    if ( iam == 0 ) WRITE(6,"(A12,F30.16)")   'D7' , cn3lo( 7)* 1.D14
    if ( iam == 0 ) WRITE(6,"(A12,F30.16)")   'D8' , cn3lo( 8)* 1.D14
    if ( iam == 0 ) WRITE(6,"(A12,F30.16)")   'D9' , cn3lo( 9)* 1.D14
    if ( iam == 0 ) WRITE(6,"(A12,F30.16)")   'D10', cn3lo(10)* 1.D14
    if ( iam == 0 ) WRITE(6,"(A12,F30.16)")   'D11', cn3lo(11)* 1.D14
    if ( iam == 0 ) WRITE(6,"(A12,F30.16)")   'D12', cn3lo(12)* 1.D14
    if ( iam == 0 ) WRITE(6,"(A12,F30.16)")   'D13', cn3lo(13)* 1.D14
    if ( iam == 0 ) WRITE(6,"(A12,F30.16)")   'D14', cn3lo(14)* 1.D14
    if ( iam == 0 ) WRITE(6,"(A12,F30.16)")   'D15', cn3lo(15)* 1.D14
    
    
    ! NNLO 3NF constants
    !
    !
    !read(5,*);read(5,*) cE, cD 
    
    

    !cE = -0.398
    !cD = -0.39
    lambdaChi = 700 ! MeV
    Econst = cE/(fpi4*lambdaChi)
    Dconst = cD/(fpi2*lambdaChi)
    
    !
    ! This is for SFR
    !
    DO i=-1,2
       IF ( (sfr-twompi(i)) <  0.0D0 ) sfr_heavyside(i) = 0.0D0
       IF ( (sfr-twompi(i)) >= 0.0D0 ) sfr_heavyside(i) = 1.0D0
    END DO
    
    z_low = 2 * mpi(2) / sfr
    IF (z_low > 1) z_low = 1
    
    
    ! z = 2*m_\pi / \mu, using the notation of [1] and [2] for the 2PE 2-loop diagrams
    ! In DR, the integration is done from z = 0 and in SFR from z = 2 * m_\pi / \Lambda_{SFR}
    CALL chp_setup_2PE_2loop_int_data(nof_2PE_2loop_int_x_points, nof_2PE_2loop_int_z_points, z_low)
    
    
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
    !22: 1/(pi^2*fpi^4)
    const(22) = 1d0/(pi2*fpi4)


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
    

    
    

  end subroutine init_chp_constants
 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      This routine calculates gauss-legendre mesh points and weights
!      INPUT:
!      mesh%info%xmin     : lower limit of the integration interval
!      mesh%info%xmax     : upper limit ---------- "" -------------
!      mesh%info%amount   : the desired number of mesh points
!      OUTPUT:
!      mesh%pnt(:)%x      : gauss-legendre mesh points on the interval (x1,x2)
!      mesh%pnt(:)%w      : the corresponding weights
!      FROM               : Numerical recipes
!      F90 version        : M. Hjorth-Jensen
!      Object interface   : M. Kartamyshev
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE chp_setup_gauleg_mesh (mesh)
    !use chiral_constants
    implicit none
    TYPE(chp_gauleg_mesh), INTENT(INOUT)        :: mesh
    INTEGER                                 :: i, j, m, n
    REAL(8)                                :: x1, x2
    REAL(8), DIMENSION(:), ALLOCATABLE     :: x, w
    REAL(8)                                :: p1,p2,p3,pp,xl,xm,z,z1
    REAL(8), PARAMETER                     :: EPS = 3.D-14
    
    ALLOCATE(x(mesh%info%amount))
    ALLOCATE(w(mesh%info%amount))
    
! allocate points and weights storages
!    CALL chp_destroy_gauleg_mesh (mesh)
    IF (ALLOCATED(mesh%pnt) ) DEALLOCATE (mesh%pnt)
    IF (mesh%info%amount <= 0 .OR. mesh%info%xmin > mesh%info%xmax) THEN
       WRITE(*,*) ': incorrect mesh info', mesh%info ; STOP
    ENDIF
    
    ALLOCATE( mesh%pnt( 1:mesh%info%amount ) )
    mesh%pnt(:) = chp_gauleg_mesh_point(0.0D0, 0.0D0, 0.0D0)
    
    ! set values of local variables
    x1 = mesh%info%xmin ; x2 = mesh%info%xmax; n = mesh%info%amount
    
    m=(n+1)/2
    xm=0.5D0*(x2+x1)
    xl=0.5D0*(x2-x1)
    DO i=1,m
       z1=0.0D0
       z=COS(pi*(i - 0.25D0)/(n + 0.5D0))
       DO WHILE ( ABS(z-z1) > EPS)
          p1=1.0D0
          p2=0.0D0
          DO j=1,n
             p3=p2
             p2=p1
             p1=((2.0D0*j-1.0D0)*z*p2-(j-1.0D0)*p3)/j
          END DO
          pp=n*(z*p1-p2)/(z*z-1.0D0)
          z1=z
          z=z-p1/pp
       END DO
       x(i)=xm-xl*z
       x(n+1-i)=xm+xl*z
       w(i)=2.0D0*xl/((1.0D0-z*z)*pp*pp)
       w(n+1-i)=w(i)
    ENDDO
    
! set return values
    mesh%pnt(:)%x = x(:) ; mesh%pnt(:)%w = w(:) ; mesh%pnt(:)%xxw = x(:) * x(:) * w(:)

    DEALLOCATE(w)
    DEALLOCATE(x)

  END SUBROUTINE chp_setup_gauleg_mesh



  !
  !     Given an NxN matrix A(N,N), this routine replaces it by the LU
  !     decomposed one, where the matrix elements are stored in the same
  !     matrix A. The array indx is  an output vector which records the row
  !     permutation effected by the partial pivoting. d is the determinant
  !
  SUBROUTINE chp_lu_decompose(a,n,indx,d)
    
    IMPLICIT NONE
    INTEGER :: n, i, j, k, imax
    REAL(8) :: sum , tiny, aamax, dum, d
    REAL(8), DIMENSION(n,n) :: a
    INTEGER, DIMENSION(n) :: indx
    REAL(8), ALLOCATABLE :: vv(:)
    
    tiny=1.0e-20
    ALLOCATE ( vv(n) )
    D=1.
    DO i=1,n
       aamax=0.
       DO j=1,n
          IF (ABS(a(i,j)) > aamax) aamax=ABS(a(i,j))
       ENDDO
       !     Zero is the largest element
       IF (aamax == 0.) STOP 'Singular matrix.'
       !     No nonzero largest element
       vv(i)=1./aamax
    ENDDO
    !     loop over columns
    DO j=1,n
       !     solves equation 2.3.12 except for i=j of Numerical Recipes
       IF (j > 1) THEN
          DO i=1,j-1
             sum=a(i,j)
             IF (i > 1)THEN
                DO k=1,i-1
                   sum=sum-a(i,k)*a(k,j)
                ENDDO
                a(i,j)=sum
             ENDIF
          ENDDO
       ENDIF
       !    start searching for largest pivot element
       aamax=0.
       DO i=j,n
          sum=a(i,j)
          IF (j > 1)THEN
             DO k=1,j-1
                sum=sum-a(i,k)*a(k,j)
             ENDDO
             a(i,j)=sum
          ENDIF
          dum=vv(i)*ABS(sum)
          IF (dum >= aamax) THEN
             imax=i
             aamax=dum
          ENDIF
       ENDDO
       !    interchange of rows
       IF (j /= imax)THEN
          DO k=1,n
             dum=a(imax,k)
             a(imax,k)=a(j,k)
             a(j,k)=dum
          ENDDO
          !    change of parity for determinant
          d=-d
          vv(imax)=vv(j)
       ENDIF
       indx(j)=imax
       IF(j /= n) THEN
          IF(a(j,j) == 0.) a(j,j)=tiny
          dum=1./a(j,j)
          DO i=j+1,n
             a(i,j)=a(i,j)*dum
          ENDDO
       ENDIF
       !    set up determinant
       d=d*a(j,j)
    ENDDO
    IF(a(n,n) == 0.)  a(n,n)=tiny
    DEALLOCATE ( vv)
    
  END SUBROUTINE chp_lu_decompose
  
  !     Solves set of linear equations Ax=b, A is input as an LU decompomsed
  !     matrix and indx keeps track of the permutations of the rows. b is input
  !     as the right-hand side vector b and returns the solution x. A, n and indx
  !     are not modified by this routine. This function takes into that b can contain
  !     many zeros and is therefore suitable for matrix inversion
  
  
  SUBROUTINE chp_lu_linear_equation(a,n,indx,b)
    
    IMPLICIT NONE
    INTEGER :: n, ii, ll, i, j
    REAL(8) :: sum
    REAL(8), DIMENSION(n,n) :: a
    REAL(8), DIMENSION(n) :: b
    INTEGER, DIMENSION(n) :: indx
    
    ii=0
    !     First we solve equation 2.3.6 of numerical recipes
    DO i=1,n
       ll=indx(i)
       sum=b(ll)
       b(ll)=b(i)
       IF (ii /= 0)THEN
          DO j=ii,i-1
             sum=sum-a(i,j)*b(j)
          ENDDO
       ELSEIF (sum /= 0.) THEN
          ii=i
       ENDIF
       b(i)=sum
    ENDDO
    !     then we solve equation 2.3.7
    DO i=n,1,-1
       sum=b(i)
       IF (i < n) THEN
          DO j=i+1,n
             sum=sum-a(i,j)*b(j)
          ENDDO
       ENDIF
       !     store a component of the solution x in the same place as b
       b(i)=sum/a(i,i)
    ENDDO
    
  END SUBROUTINE chp_lu_linear_equation
  
   !            Routines to do mtx inversion, from Numerical
  !            Recepies, Teukolsky et al. Routines included
  !            below are MATINV, LUDCMP and LUBKSB. See chap 2
  !            of Numerical Recipes for further details
  !            Recoded in FORTRAN 90 by M. Hjorth-Jensen
  !
  SUBROUTINE chp_matinv(a,n)
    
    IMPLICIT NONE
    
    INTEGER, INTENT(IN) :: n
    INTEGER :: i, j
    REAL(8), DIMENSION(n,n), INTENT(INOUT)  :: a
    REAL(8), ALLOCATABLE :: y(:,:)
    REAL(8) :: d
    INTEGER, ALLOCATABLE :: indx(:)
    
    ALLOCATE (y( n, n))  ; ALLOCATE ( indx (n))
    y=0.
    !     setup identity matrix
    DO i=1,n
       y(i,i)=1.
    ENDDO
    !     LU decompose the matrix just once
    CALL  chp_lu_decompose(a,n,indx,d)
    
    !     Find inverse by columns
    DO j=1,n
       CALL chp_lu_linear_equation(a,n,indx,y(:,j))
    ENDDO
    !     The original matrix a was destroyed, now we equate it with the inverse y
    a=y
    
    DEALLOCATE ( y ); DEALLOCATE ( indx )
    
  END SUBROUTINE chp_matinv
  






  subroutine init_cEcD
    
    real*8 :: lambdachi
    lambdaChi = 700 ! MeV
    Econst = cE/(fpi4*lambdaChi)
    Dconst = cD/(fpi2*lambdaChi)


  end subroutine init_cEcD

! z = 2*m_\pi / \mu, using the notation of [1] and [2] for the 2PE 2-loop diagrams
! In DR, the integration is done from z_low = 0 and in SFR from z_low = 2 * m_\pi / \Lambda_{SFR}
! The expressions for some of the 2PE 2-loop diagrams contains integrals that must be solved numerically
! This function sets up the integration meshs used in the numerical integrations
  SUBROUTINE chp_setup_2PE_2loop_int_data(nof_x_points, nof_z_points, z_low)
    INTEGER, INTENT(IN) :: nof_x_points
    INTEGER, INTENT(IN) :: nof_z_points
    REAL(8), INTENT(IN) :: z_low

    chp_2PE_2loop_int_mesh_x%info = chp_mesh_info(nof_x_points, 0.0D0, 1.0D0)
    CALL chp_setup_gauleg_mesh(chp_2PE_2loop_int_mesh_x)

    chp_2PE_2loop_int_mesh_z%info = chp_mesh_info(nof_z_points, z_low, 1.0D0)
    CALL chp_setup_gauleg_mesh(chp_2PE_2loop_int_mesh_z)

! TODO: Magic number: 0.95. The reason for using 0.95 is that the z integral is a bit
! tricky near z = 1, therefore I do numerical integration up to 0.95, then I use the
! integral of the taylor expansion around z=1 from 0.95 to 1.0 to get the last part.
! This is needed to get accurate results.
    chp_2PE_2loop_int_mesh_z_VC_SFR%info = chp_mesh_info(nof_z_points, z_low, 0.95D0)
    CALL chp_setup_gauleg_mesh(chp_2PE_2loop_int_mesh_z_VC_SFR)

    chp_2PE_2loop_int_VTS_data_set     = .false.
    !    chp_2PE_2loop_int_WC_DR_data_set   = .false.
    chp_2PE_2loop_int_WC_SFR_data_set  = .false.
  END SUBROUTINE


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
  
  
  ! NLO function [1] 4.9 OR [2] 2.14 (W_C part) OR [3] B1 OR [4] 6
! q2  : momentum transfer squared
! impi: determines which mpi2 to use,
  FUNCTION chp_two_pion_exchange_1loop_0_Wc(q2, L, w, impi) RESULT(res)
    
    REAL(8) , INTENT(IN) :: q2
    REAL(8) , INTENT(IN) :: L
    REAL(8) , INTENT(IN) :: w
    INTEGER, INTENT(IN) :: impi
    REAL(8) :: res
    
    res = -L *(fourmpi2(impi)*const(4) + q2*const(5) + const(6)*mpi4(impi)/(w*w))/const(3)
    
  END FUNCTION chp_two_pion_exchange_1loop_0_Wc
  
! NLO function [1] 4.10 OR [2] 2.14 (V_S part) OR [3] B2 OR [4] 6
! q2  : momentum transfer squared
! impi: determines which mpi2 to use,
  FUNCTION chp_two_pion_exchange_1loop_0_Vs(q2,L,impi) RESULT(res)
    
    REAL(8) , INTENT(IN) :: q2
    REAL(8) , INTENT(IN) :: L
    INTEGER, INTENT(IN) :: impi
    REAL(8) :: res
    
    res = const(7)*L*q2/const(8)
    
  END FUNCTION chp_two_pion_exchange_1loop_0_Vs

! NLO function [1] 4.10 OR [2] 2.14 (V_T part) OR [3] B2 OR [4] 6
! q2  : momentum transfer squared
! impi: determines which mpi2 to use,
  FUNCTION chp_two_pion_exchange_1loop_0_VT(L,impi) RESULT(res)
    
    REAL(8) , INTENT(IN) :: L
    INTEGER, INTENT(IN) :: impi
    REAL(8) :: res
    
    res = -const(7)*L/const(8)
    
  END FUNCTION chp_two_pion_exchange_1loop_0_VT


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
    res = const(7)*L*q2/const(8)
    
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
    
    eps = 1.0D-8 
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
  
   ! NNLO loop function wtilde SQUARED [1] Eq 4.20 (DR)
  ! q2  : momentum transfer squared
  ! impi: determines which mpi2 to use, 
  FUNCTION chp_NNLO_two_pion_exchange_loop_wtilde2(q2, impi) RESULT(res)
    
    REAL*8 , INTENT(IN) :: q2
    INTEGER, INTENT(IN) :: impi
    REAL*8 :: res
    
    res = 0.0D0
    res = 2.0D0*mpi2(impi) + q2
    
  END FUNCTION chp_NNLO_two_pion_exchange_loop_wtilde2

  ! NNLO loop function wtilde [2] Eq 2.17 (SFR)
  ! q   : momentum transfer
  ! q2  : momentum transfer squared
  ! impi: determines which mpi to use, 
  FUNCTION chp_NNLO_sfr_two_pion_exchange_loop_A(q, q2, impi) RESULT(res)
    
    REAL*8 , INTENT(INOUT) :: q
    REAL*8 , INTENT(INOUT) :: q2
    INTEGER, INTENT(IN) :: impi
    REAL*8 :: res
    REAL*8 :: eps
    
    res = 0.0D0
    eps = 1.0D-8
    IF (sfr_heavyside(impi) == 0.0D0) return
    
    if ( q == 0.d0 ) then  
       q = q + eps 
       q2 = q*q 
    end if
    
    res = datan( q*(sfr-twompi(impi) )/(q2 + sfr*twompi(impi) ) )/(2.0D0*q)
    
  END FUNCTION chp_NNLO_sfr_two_pion_exchange_loop_A
  
  ! NNLO function [1] Eq 4.13 (ci part) OR [2] Eq 2.15 (V_C part) OR [3] Eq C1
! q2    : momentum transfer squared
! A     : NNLO loop function A [1]Eq. 4.19
! wt2   : NNLO loop function wtilde Squared(!) ([1] Eq. 4.20)
! impi  : determines which mpi2 to use
  FUNCTION chp_two_pion_exchange_1loop_d_Vc(q2, A, wt2, impi) RESULT(res)
    
    REAL(8) , INTENT(IN) :: q2
    REAL(8) , INTENT(IN) :: A
    REAL(8) , INTENT(IN) :: wt2
    INTEGER, INTENT(IN) :: impi
    REAL(8) :: res
    
    res = const(9)*(- &
            (mpi2(impi)*const(11) - q2*c3) *wt2*A)

  END FUNCTION chp_two_pion_exchange_1loop_d_Vc

! NNLO/N3LO function [1] Eq 4.13 (M_N part) and Eq 4.21 OR [2] Eq 2.23 (V_C part) OR [3] Eq D7
! q2    : momentum transfer squared
! w     : Eq 4.12
! A     : NNLO loop function A Eq. 4.19
! wt2   : NNLO loop function wtilde Squared(!) (Eq. 4.20)
! impi  : determines which mpi2 to use
! chp_itope;  EM or KW, choose which iterated 2pe 'model' to use
! this is set in the module header
  FUNCTION chp_two_pion_exchange_1loop_r_Vc(q2, w, A, wt2, impi, imnuc) RESULT(res)
    
    REAL(8) , INTENT(IN) :: q2
    REAL(8) , INTENT(IN) :: w
    REAL(8) , INTENT(IN) :: A
    REAL(8) , INTENT(IN) :: wt2
    INTEGER, INTENT(IN) :: impi
    INTEGER, INTENT(IN) :: imnuc
    REAL(8) :: res
    
!    elseif(chiral_2PE_1loop_r_mode%val == chiral_2PE_1loop_r_mode_EM2015) then
! [3] D7 OR [4] 19+22
       res = 3*const(14)*(mpi5(impi)/(2*w*w) + (2*mpi2(impi) + q2) * (q2 - mpi2(impi))*A) / mnuc(imnuc)

  END FUNCTION chp_two_pion_exchange_1loop_r_Vc
  
! NNLO function Eq [1] 4.14 and 4.22
! q2    : momentum transfer squared
! w     : Eq 4.12
! A     : NNLO loop function A Eq. 4.19
! wt2   : NNLO loop function wtilde Squared(!) (Eq. 4.20)
! impi  : determines which mpi2 to use
! chp_itope;  EM or KW, choose which iterated 2pe 'model' to use
! this is set in the module header
  FUNCTION chp_two_pion_exchange_1loop_r_Wc(q2, w, A, wt2, impi, imnuc) RESULT(res)
    
    REAL(8) , INTENT(IN) :: q2
    REAL(8) , INTENT(IN) :: w
    REAL(8) , INTENT(IN) :: A
    REAL(8) , INTENT(IN) :: wt2
    INTEGER, INTENT(IN) :: impi
    INTEGER, INTENT(IN) :: imnuc
    REAL(8) :: res
    

![3] D8 = [4] 19 + 22
       res = 2.0D0*const(13) * (3.0D0*gA2*mpi5(impi)/(2.0D0*w*w) + &
             (gA2*(3.0D0*mpi2(impi) + 2.0D0*q2) - 2.0D0*mpi2(impi) - q2)*&
             (2.0D0*mpi2(impi) + q2)*A)/mnuc(imnuc)

    
  END FUNCTION chp_two_pion_exchange_1loop_r_Wc
  
! NNLO function Eq [1] 4.15 and 4.23
! w     : Eq 4.12
! A     : NNLO loop function A Eq. 4.19
! wt2   : NNLO loop function wtilde Squared(!) (Eq. 4.20)
! impi  : determines which mpi2 to use
! chp_itope;  EM or KW, choose which iterated 2pe 'model' to use
! this is set in the module header
  FUNCTION chp_two_pion_exchange_1loop_r_VT(q2, w, A, wt2, impi, imnuc) RESULT(res)

    REAL(8) , INTENT(IN) :: q2
    REAL(8) , INTENT(IN) :: w
    REAL(8) , INTENT(IN) :: A
    REAL(8) , INTENT(IN) :: wt2
    INTEGER, INTENT(IN) :: impi
    INTEGER, INTENT(IN) :: imnuc
    REAL(8) :: res
    
![3] D9 = [4] 19 + 22
       res = 3.0D0*const(18)*(5.0D0*mpi2(impi) + 2.0D0*q2)*A/mnuc(imnuc)
    
  END FUNCTION chp_two_pion_exchange_1loop_r_VT

  ! NNLO function Eq [1] 4.15 and 4.23
! q2    : momentum transfer squared
! w     : Eq 4.12
! A     : NNLO loop function A Eq. 4.19
! wt2   : NNLO loop function wtilde Squared(!) (Eq. 4.20)
! impi  : determines which mpi2 to use
! chp_itope;  EM or KW, choose which iterated 2pe 'model' to use
! this is set in the module header
  FUNCTION chp_two_pion_exchange_1loop_r_Vs(q2, w, A, wt2, impi, imnuc) RESULT(res)
    
    REAL(8) , INTENT(IN) :: q2
    REAL(8) , INTENT(IN) :: w
    REAL(8) , INTENT(IN) :: A
    REAL(8) , INTENT(IN) :: wt2
    INTEGER, INTENT(IN) :: impi
    INTEGER, INTENT(IN) :: imnuc
    REAL(8) :: res
    
    res = -q2 * chp_two_pion_exchange_1loop_r_VT(q2, w, A, wt2, impi, imnuc)

  END FUNCTION chp_two_pion_exchange_1loop_r_Vs
  
! NNLO function Eq [1] 4.16
! q2    : momentum transfer squared
! w     : Eq 4.12
! A     : NNLO loop function A Eq. 4.19
! impi  : determines which mpi2 to use
! chp_itope;  EM or KW, choose which iterated 2pe 'model' to use
! this is set in the module header
  FUNCTION chp_two_pion_exchange_1loop_d_WT(w, A) RESULT(res)
    
    REAL(8) , INTENT(IN) :: w
    REAL(8) , INTENT(IN) :: A
    REAL(8) :: res
    
! [1] 4.16
    res = -1.0D0*const(16)*A*( c4*w*w )
    
  END FUNCTION chp_two_pion_exchange_1loop_d_WT
  
  
! NNLO function Eq [1] 4.16 and 4.24
! q2    : momentum transfer squared
! w     : Eq 4.12
! A     : NNLO loop function A Eq. 4.19
! impi  : determines which mpi2 to use
! chp_itope;  EM or KW, choose which iterated 2pe 'model' to use
! this is set in the module header
  FUNCTION chp_two_pion_exchange_1loop_r_WT(q2, w, A, impi, imnuc) RESULT(res)
    
    REAL(8) , INTENT(IN) :: q2
    REAL(8) , INTENT(IN) :: w
    REAL(8) , INTENT(IN) :: A
    INTEGER, INTENT(IN) :: impi
    INTEGER, INTENT(IN) :: imnuc
    REAL(8) :: res
    
! [3] D10 = [4] 19 + 22
       res = const(13)*(gA2*(3.0D0*mpi2(impi) + q2) - w**2)*A/mnuc(imnuc)
 
    
  END FUNCTION chp_two_pion_exchange_1loop_r_WT

! NNLO function Eq [1] 4.16
! q2    : momentum transfer squared
! w     : Eq 4.12
! A     : NNLO loop function A Eq. 4.19
! impi  : determines which mpi2 to use
! chp_itope;  EM or KW, choose which iterated 2pe 'model' to use
! this is set in the module header
  FUNCTION chp_two_pion_exchange_1loop_d_Ws(q2, w, A) RESULT(res)
    
    REAL(8) , INTENT(IN) :: q2
    REAL(8) , INTENT(IN) :: w
    REAL(8) , INTENT(IN) :: A
    REAL(8) :: res
    
    res = -q2 * chp_two_pion_exchange_1loop_d_WT(w, A)
    
  END FUNCTION chp_two_pion_exchange_1loop_d_Ws
  
  ! NNLO function Eq [1] 4.16 and 4.24
! q2    : momentum transfer squared
! w     : Eq 4.12
! A     : NNLO loop function A Eq. 4.19
! impi  : determines which mpi2 to use
! chp_itope;  EM or KW, choose which iterated 2pe 'model' to use
! this is set in the module header
  FUNCTION chp_two_pion_exchange_1loop_r_Ws(q2, w, A, impi, imnuc) RESULT(res)
    
    REAL(8) , INTENT(IN) :: q2
    REAL(8) , INTENT(IN) :: w
    REAL(8) , INTENT(IN) :: A
    INTEGER, INTENT(IN) :: impi
    INTEGER, INTENT(IN) :: imnuc
    REAL(8) :: res
    
    res = -q2 * chp_two_pion_exchange_1loop_r_WT(q2, w, A, impi, imnuc)
    
  END FUNCTION chp_two_pion_exchange_1loop_r_Ws

! NNLO function Eq [1] 4.17
! A     : NNLO loop function A Eq. 4.19
! wt2   : NNLO loop function wtilde Squared(!) (Eq. 4.20)
! impi  : determines which mpi2 to use
  FUNCTION chp_two_pion_exchange_1loop_r_VLS(A, wt2, imnuc) RESULT(res)
    
    REAL(8) , INTENT(IN) :: A
    REAL(8) , INTENT(IN) :: wt2
    INTEGER, INTENT(IN) :: imnuc
    REAL(8) :: res
    
    res = const(19) * wt2*A/mnuc(imnuc)
             
  END FUNCTION chp_two_pion_exchange_1loop_r_VLS

! NNLO function Eq [1] 4.18
! w   : NNLO loop function wtilde Squared(!) (Eq. 4.20)
! A     : NNLO loop function A Eq. 4.19
! impi  : determines which mpi2 to use
  FUNCTION chp_two_pion_exchange_1loop_r_WLS(w, A, imnuc) RESULT(res)
    
    REAL(8) , INTENT(IN) :: w
    REAL(8) , INTENT(IN) :: A
    INTEGER, INTENT(IN) :: imnuc
    REAL(8) :: res
    
    res = const(20)*w*w*A/mnuc(imnuc)
             
  END FUNCTION chp_two_pion_exchange_1loop_r_WLS
  
  
  
  ! NNLO function Eq [1] 4.13 and 4.21
  ! q2    : momentum transfer squared
  ! w     : Eq 4.12
  ! A     : NNLO loop function A Eq. 4.19
  ! wt2   : NNLO loop function wtilde Squared(!) (Eq. 4.20)
  ! impi  : determines which mpi2 to use
  ! chp_itope;  EM or KW, choose which iterated 2pe 'model' to use
  ! this is set in the module header
  FUNCTION chp_NNLO_two_pion_exchange_Vc(q2, w, A, wt2, impi) RESULT(res)
    
    REAL*8 , INTENT(IN) :: q2
    REAL*8 , INTENT(IN) :: w
    REAL*8 , INTENT(IN) :: A
    REAL*8 , INTENT(IN) :: wt2
    INTEGER, INTENT(IN) :: impi
    REAL*8 :: res
    
!!$    res = 0.0D0
!!$    res = const(9)*(const(10)*mpi5(impi)/(mnuc(0)*w*w) - &
!!$         (mpi2(impi)*const(11) - q2*c3 - q2*3.0D0*const(10)/mnuc(0) ) *wt2*A)
!!$    
!!$    ! IF twopion exchange in BbS , i.e. EM format
!!$    ! add correction to the NNLO central term
!!$    IF (chp_itope == 'EM') THEN
!!$       res = res - const(12)*(mpi(impi)*w*w+wt2*wt2*A)/mnuc(0)
!!$    END IF

!  new power counting scheme
    !from 1loop_d_Vc
    res = const(9)*(-(mpi2(impi)*const(11) - q2*c3) *wt2*A)
    !from 1loop_r_Vc
     res = res + 3*const(14)*(mpi5(impi)/(2*w*w) + (2*mpi2(impi) + q2) * (q2 - mpi2(impi))*A) / mnuc(0)
     
  END FUNCTION chp_NNLO_two_pion_exchange_Vc


  ! NNLO function Eq [1] 4.14 and 4.22
  ! q2    : momentum transfer squared
  ! w     : Eq 4.12
  ! A     : NNLO loop function A Eq. 4.19
  ! wt2   : NNLO loop function wtilde Squared(!) (Eq. 4.20)
  ! impi  : determines which mpi2 to use
  ! chp_itope;  EM or KW, choose which iterated 2pe 'model' to use
  ! this is set in the module header
  FUNCTION chp_NNLO_two_pion_exchange_Wc(q2, w, A, wt2, impi) RESULT(res)
    
    REAL*8 , INTENT(IN) :: q2
    REAL*8 , INTENT(IN) :: w
    REAL*8 , INTENT(IN) :: A
    REAL*8 , INTENT(IN) :: wt2
    INTEGER, INTENT(IN) :: impi
    REAL*8 :: res
    
!!$    res = 0.0D0
!!$    res = const(13) * (3.0D0*gA2*mpi5(impi)/(w*w) - & 
!!$         (fourmpi2(impi) + 2.0D0*q2 - gA2*(fourmpi2(impi)+3.0D0*q2))*wt2*A)/mnuc(0)
!!$             
!!$    ! IF twopion exchange in BbS , i.e. EM format
!!$    ! add correction to the NNLO central term
!!$    IF (chp_itope == 'EM') THEN
!!$       res = res + const(14)*(mpi(impi)*w*w + wt2*wt2*A)/mnuc(0)
!!$    END IF
    
   !from loop_r_Wc 
    res =   2.0D0*const(13) * (3.0D0*gA2*mpi5(impi)/(2.0D0*w*w) + &
             (gA2*(3.0D0*mpi2(impi) + 2.0D0*q2) - 2.0D0*mpi2(impi) - q2)*&
             (2.0D0*mpi2(impi) + q2)*A)/mnuc(0)

  END FUNCTION chp_NNLO_two_pion_exchange_Wc
  
  ! NNLO function Eq [1] 4.15 and 4.23
  ! w     : Eq 4.12
  ! A     : NNLO loop function A Eq. 4.19
  ! wt2   : NNLO loop function wtilde Squared(!) (Eq. 4.20)
  ! impi  : determines which mpi2 to use
  ! chp_itope;  EM or KW, choose which iterated 2pe 'model' to use
  ! this is set in the module header
  FUNCTION chp_NNLO_two_pion_exchange_VT(w, A, wt2, impi) RESULT(res)
    
    REAL*8 , INTENT(IN) :: w
    REAL*8 , INTENT(IN) :: A
    REAL*8 , INTENT(IN) :: wt2
    INTEGER, INTENT(IN) :: impi
    REAL*8 :: res
!!$    
!!$    res = 0.0D0
!!$    res = 3.0D0*const(15)*wt2*A/mnuc(0)
!!$             
!!$    ! IF twopion exchange in BbS , i.e. EM format
!!$    ! add correction to the NNLO central term
!!$    IF (chp_itope == 'EM') THEN
!!$       res = res + const(15)*(mpi(impi) + w*w*A )/mnuc(0)
!!$    END IF
    
    !from  1loop_r_VT
    ![3] D9 = [4] 19 + 22
    res = 3.0D0*const(18)*(5.0D0*mpi2(impi) + 2.0D0*q2)*A/mnuc(0)
       
       
  END FUNCTION chp_NNLO_two_pion_exchange_VT

  ! NNLO function Eq [1] 4.15 and 4.23
  ! q2    : momentum transfer squared
  ! w     : Eq 4.12
  ! A     : NNLO loop function A Eq. 4.19
  ! wt2   : NNLO loop function wtilde Squared(!) (Eq. 4.20)
  ! impi  : determines which mpi2 to use
  ! chp_itope;  EM or KW, choose which iterated 2pe 'model' to use
  ! this is set in the module header
  FUNCTION chp_NNLO_two_pion_exchange_Vs(q2, w, A, wt2, impi) RESULT(res)
    
    REAL*8 , INTENT(IN) :: q2
    REAL*8 , INTENT(IN) :: w
    REAL*8 , INTENT(IN) :: A
    REAL*8 , INTENT(IN) :: wt2
    INTEGER, INTENT(IN) :: impi
    REAL*8 :: res
    
!!$    res = 0.0D0
!!$    res = -3.0D0*q2*const(15)*wt2*A/mnuc(0)
!!$
!!$    ! IF twopion exchange in BbS , i.e. EM format
!!$    ! add correction to the NNLO central term
!!$    IF (chp_itope == 'EM') THEN
!!$       res = res - 1.0D0*const(15)*q2*(mpi(impi) + w*w*A)/mnuc(0)
!!$    END IF 

    res = -q2*3.0D0*const(18)*(5.0D0*mpi2(impi) + 2.0D0*q2)*A/mnuc(0)
    
    !WRITE(*,*)q2
    !WRITE(*,*)const(15)
    !WRITE(*,*)wt2
    !WRITE(*,*)A
    !WRITE(*,*)mnuc(0)
    !WRITE(*,*)const(15)
    !WRITE(*,*)mpi(impi)
    !WRITE(*,*)w*w
    !WRITE(*,*)res
    
  END FUNCTION chp_NNLO_two_pion_exchange_Vs
  
  ! NNLO function Eq [1] 4.16 and 4.24
  ! q2    : momentum transfer squared
  ! w     : Eq 4.12
  ! A     : NNLO loop function A Eq. 4.19
  ! impi  : determines which mpi2 to use
  ! chp_itope;  EM or KW, choose which iterated 2pe 'model' to use
  ! this is set in the module header
  FUNCTION chp_NNLO_two_pion_exchange_WT(q2, w, A, impi) RESULT(res)
    
    REAL*8 , INTENT(IN) :: q2
    REAL*8 , INTENT(IN) :: w
    REAL*8 , INTENT(IN) :: A
    INTEGER, INTENT(IN) :: impi
    REAL*8 :: res
!!$    
!!$    res = 0.0D0
!!$    res = -1.0D0*const(16)*A*( (c4 + 1.0D0/(4.0D0*mnuc(0)))*w*w - &
!!$         const(17)*(10.0D0*mpi2(impi) + 3.0D0*q2)/mnuc(0))
!!$             
!!$    ! IF twopion exchange in BbS , i.e. EM format
!!$    ! add correction to the NNLO central term
!!$    IF (chp_itope == 'EM') THEN
!!$       res = res - const(18)*(mpi(impi) + w*w*A)/mnuc(0)
!!$    END IF
    
    res =  -1.0D0*const(16)*A*( c4*w*w )
    ! [3] D10 = [4] 19 + 22
    res = res + const(13)*(gA2*(3.0D0*mpi2(impi) + q2) - w**2)*A/mnuc(0)
    
    
  END FUNCTION chp_NNLO_two_pion_exchange_WT 

  ! NNLO function Eq [1] 4.16 and 4.24
  ! q2    : momentum transfer squared
  ! w     : Eq 4.12
  ! A     : NNLO loop function A Eq. 4.19
  ! impi  : determines which mpi2 to use
  ! chp_itope;  EM or KW, choose which iterated 2pe 'model' to use
  ! this is set in the module header
  FUNCTION chp_NNLO_two_pion_exchange_Ws(q2, w, A, impi) RESULT(res)
    
    REAL*8 , INTENT(IN) :: q2
    REAL*8 , INTENT(IN) :: w
    REAL*8 , INTENT(IN) :: A
    INTEGER, INTENT(IN) :: impi
    REAL*8 :: res
    
!!$    res = 0.0D0
!!$    res = q2*const(16)*A*( (c4+1.0D0/(4.0D0*mnuc(0)))*w*w - & 
!!$         const(17)*(10.0D0*mpi2(impi) +3.0D0*q2)/mnuc(0))
!!$             
!!$    ! IF twopion exchange in BbS , i.e. EM format
!!$    ! add correction to the NNLO central term
!!$    IF (chp_itope == 'EM') THEN
!!$       res = res + const(18)*q2*(mpi(impi) + w*w*A)/mnuc(0)
!!$    END IF
    
    res =  q2*1.0D0*const(16)*A*( c4*w*w )
    ! [3] D10 = [4] 19 + 22
    res = res -q2*const(13)*(gA2*(3.0D0*mpi2(impi) + q2) - w**2)*A/mnuc(0)
    !res = -q2*res
    
    
  END FUNCTION chp_NNLO_two_pion_exchange_Ws

  ! NNLO function Eq [1] 4.17
  ! A     : NNLO loop function A Eq. 4.19
  ! wt2   : NNLO loop function wtilde Squared(!) (Eq. 4.20)
  ! impi  : determines which mpi2 to use
  FUNCTION chp_NNLO_two_pion_exchange_VLS(A, wt2, impi) RESULT(res)
    
    REAL*8 , INTENT(IN) :: A
    REAL*8 , INTENT(IN) :: wt2
    INTEGER, INTENT(IN) :: impi
    REAL*8 :: res
    
    res = 0.0D0
    res = const(19) * wt2*A/mnuc(0)
             
  END FUNCTION chp_NNLO_two_pion_exchange_VLS

  ! NNLO function Eq [1] 4.18
  ! w   : NNLO loop function wtilde Squared(!) (Eq. 4.20)
  ! A     : NNLO loop function A Eq. 4.19
  ! impi  : determines which mpi2 to use
  FUNCTION chp_NNLO_two_pion_exchange_WLS(w, A, impi) RESULT(res)
    
    REAL*8 , INTENT(IN) :: w
    REAL*8 , INTENT(IN) :: A
    INTEGER, INTENT(IN) :: impi
    REAL*8 :: res
    
    res = 0.0D0
    res = const(20)*w*w*A/mnuc(0)
             
  END FUNCTION chp_NNLO_two_pion_exchange_WLS


  ! N3LO two pion exchange functions
  ! k2    : mean momentum squared
  ! q2    : momentum transfer squared
  ! w     : Eq 4.12
  ! L     : L loop function (depends on what regularization is used)
  ! wt2   : NNLO loop function wtilde Squared(!) (Eq. 4.20)
  ! A     : NNLO loop function A Eq. 4.19
  ! impi  : determines which mpi2 to use
  
  ! [1]Eq. D.1
  FUNCTION chp_N3LO_2PE_Vc_1loop_ci2 (w, L, wt2, impi) RESULT(res)
    REAL(8) , INTENT(IN) :: w
    REAL(8) , INTENT(IN) :: L
    REAL(8) , INTENT(IN) :: wt2
    INTEGER, INTENT(IN) :: impi
    REAL(8) :: res
    
    res = c2 * (1D0/6D0) * w * w + c3 * wt2 - c1 * fourmpi2(impi)
    res = ((3D0 / 16D0) * const(22)) * L * (res*res + c2*c2*w*w*w*W*(1D0/45D0))
  END FUNCTION chp_N3LO_2PE_Vc_1loop_ci2

  ! [1]Eq. D.4
  FUNCTION chp_N3LO_2PE_Vc_1loop_ci1 ( q2, w, L, impi, imnuc) RESULT(res)
    REAL(8) , INTENT(IN) :: q2
    REAL(8) , INTENT(IN) :: w
    REAL(8) , INTENT(IN) :: L
    INTEGER, INTENT(IN) :: impi
    INTEGER, INTENT(IN) :: imnuc
    REAL(8) :: res
    
    res = (-gA2 * const(22) * mnuc_inv(imnuc) * (1D0 / 32D0)) * L * ( &
              (c2 - 6D0*c3) * q2 * q2 + 4D0 * (6D0*c1 + c2 - 3D0*c3)*q2*mpi2(impi) + &
              6D0 * (c2 - 2D0*c3)*mpi4(impi) + &
              24D0*(2D0*c1 + c3)*mpi4(impi)*mpi2(impi)/(w*w) )
  END FUNCTION chp_N3LO_2PE_Vc_1loop_ci1
  
  ! [1]Eq. D.9
  FUNCTION chp_N3LO_2PE_Vc_1loop_ci0 (q2, w, L, impi, imnuc) RESULT(res)
    REAL(8) , INTENT(IN) :: q2
    REAL(8) , INTENT(IN) :: w
    REAL(8) , INTENT(IN) :: L
    INTEGER, INTENT(IN) :: impi
    INTEGER, INTENT(IN) :: imnuc
    REAL(8) :: res
    
    res = -gA2*gA2*const(22)*mnuc_inv(imnuc)*mnuc_inv(imnuc)*(1D0/32D0) * ( &
              L*(2D0*mpi4(impi)*mpi4(impi)/(w*w*w*w) + 8D0*mpi3(impi)*mpi3(impi)/(w*w) - &
                  q2*q2 - 2D0*mpi4(impi) ) + mpi3(impi)*mpi3(impi)*0.5D0/(w*w) )
  END FUNCTION chp_N3LO_2PE_Vc_1loop_ci0

  ! [1]Eq. D.18
  FUNCTION chp_N3LO_2PE_Vc_2loop(q2, wt2, A, impi) RESULT(res)
    REAL(8) , INTENT(IN) :: q2
    REAL(8) , INTENT(IN) :: wt2
    REAL(8) , INTENT(IN) :: A
    INTEGER, INTENT(IN) :: impi
    REAL(8) :: res
    
    res = 3D0*gA2*wt2*A*const(22)*const(2)*(1D0/256D0)*( &
              (mpi2(impi) + 2D0*q2)*(2D0*mpi(impi) + wt2*A) + &
              4D0*gA2*mpi(impi)*wt2)
  END FUNCTION chp_N3LO_2PE_Vc_2loop     
  
  
!!!!!   XXXXXXXXX   New Stuff XXXXXXXXXX 
  ! [1]Eq. D.5
  FUNCTION chp_N3LO_2PE_Wc_1loop_ci1 (q2, w, L, impi, imnuc) RESULT(res)
    REAL(8) , INTENT(IN) :: q2
    REAL(8) , INTENT(IN) :: w
    REAL(8) , INTENT(IN) :: L
    INTEGER, INTENT(IN) :: impi
    INTEGER, INTENT(IN) :: imnuc
    REAL(8) :: res
    
    res = -c4*q2*L*const(22)*mnuc_inv(imnuc)*(1D0/192D0) * ( &
              gA2*(8D0*mpi2(impi) + 5D0*q2) + w*w)
  END FUNCTION chp_N3LO_2PE_Wc_1loop_ci1

 ! [1]Eq. D.10
  FUNCTION chp_N3LO_2PE_Wc_1loop_ci0 (k2,q2, w, L, impi, imnuc) RESULT(res)
    REAL(8) , INTENT(IN) :: k2
    REAL(8) , INTENT(IN) :: q2
    REAL(8) , INTENT(IN) :: w
    REAL(8) , INTENT(IN) :: L
    INTEGER, INTENT(IN) :: impi
    INTEGER, INTENT(IN) :: imnuc
    REAL(8) :: res
    
    res = -const(22)*mnuc_inv(imnuc)*mnuc_inv(imnuc)*(1D0/768D0)*( &
              L*(8D0*gA2*(1.5D0*q2*q2 + 3D0*mpi2(impi)*q2 + 3D0*mpi4(impi) - &
                      6D0*mpi3(impi)*mpi3(impi)/(w*w) - k2*(8D0*mpi2(impi) + 5D0*q2)) + &
                  4D0*gA2*gA2*(k2*(20D0*mpi2(impi) + 7D0*q2 - 16D0*mpi4(impi)/(w*w)) + &
                      16D0*mpi4(impi)*mpi4(impi)/(w*w*w*w) + &
                      12D0*mpi3(impi)*mpi3(impi)/(w*w) - 4D0*mpi4(impi)*q2/(w*w) - &
                      5D0*q2*q2 - 6D0*mpi2(impi)*q2 - 6D0*mpi4(impi)) - 4D0*k2*w*w) + &
              16D0*gA2*gA2*mpi3(impi)*mpi3(impi)/(w*w))
  END FUNCTION chp_N3LO_2PE_Wc_1loop_ci0

  ! [1]Eq. D.20
  FUNCTION chp_N3LO_2PE_Wc_2loop_a (q2, w, L, wt2,impi) RESULT(res)
    REAL(8) , INTENT(IN) :: q2
    REAL(8) , INTENT(IN) :: w
    REAL(8) , INTENT(IN) :: L
    REAL(8) , INTENT(IN) :: wt2
    INTEGER, INTENT(IN) :: impi
    REAL(8) :: res
    
    ! const(22) = 1 / (pi^2 * f_pi^4)
    !!  di's are the scale-independent versions
    res = L*const(22)*fpi_inv*fpi_inv*pi_inv*pi_inv*(1D0/18432D0)*( &
              192D0*pi2*fpi2*w*w*d3*(2D0*gA2*wt2 - 0.6D0*(gA2 - 1D0)*w*w) + &
              (6D0*gA2*wt2 - (gA2 - 1D0)*w*w)*(384D0*pi2*fpi2*(wt2*d1_plus_d2 + &
                      4D0*mpi2(impi)*d5) + L*(4D0*mpi2(impi)*(1D0+2D0*gA2) + &
                          q2*(1D0 + 5D0*gA2)) - (q2*(1D0/3D0)*(5D0 + 13D0*gA2) + &
                              8D0*mpi2(impi)*(1D0 + 2D0*gA2))))
  END FUNCTION chp_N3LO_2PE_Wc_2loop_a   

  !!  According to Machleidt Appendix D loop(b) contributions are negligible for scattering
  
  ! Helper function needed by chp_calculate_2PE_2loop_int_*_data functions
  FUNCTION chp_N3LO_2PE_2loop_int_helper(x) RESULT(res)
    REAL(8), INTENT(IN) :: x
    REAL(8) :: x2, x4, x2_inv, A
    REAL(8) :: res

    IF(x < 0.05) THEN
      ! Use Taylor expansion around 0 with terms up to x^6 for small x
      x2 = x*x
      x4 = x2*x2

      res = (((-8.0D0/315.0D0)*x4*x2 + (2.0D0/35.0D0)*x4) - 0.2*x2) - 4.0D0/3.0D0
    ELSE
      x2_inv = 1.0D0 / (x*x)
      A = sqrt(1 + x2_inv)

      res = x2_inv - (1+x2_inv)*A * log(x * (1 + A))
    END IF
  END FUNCTION chp_N3LO_2PE_2loop_int_helper
!!$
!!$  ! [1]Eq. D.21 and D.22
!!$  SUBROUTINE chp_calculate_2PE_2loop_int_WC_data
!!$    INTEGER :: z_nr, x_nr, i, j
!!$    REAL(8) :: fact, C1
!!$    REAL(8) :: res, z, z2, zt, D1
!!$    REAL(8) :: w, x, y, y2, y_z
!!$
!!$    IF(chp_2PE_2loop_int_WC_data_set) return
!!$
!!$    z_nr = chp_2PE_2loop_int_mesh_z%info%amount
!!$    x_nr = chp_2PE_2loop_int_mesh_x%info%amount
!!$    fact = 1.0D0 / (2048.0D0 * pi2 * pi2)
!!$    C1 = 2.0D0 * gA2 * gA2 / 3.0D0
!!$
!!$!$OMP parallel do schedule(static) default(none) shared(z_nr, x_nr, fact, C1, chp_2PE_2loop_int_mesh_z, gA2, chp_2PE_2loop_int_mesh_x, chp_2PE_2loop_int_WC_x_fact) private(i, j, res, z, z2, zt, D1, w, x, y, y2, y_z)
!!$    DO i = 1, z_nr
!!$      res = 0
!!$      z = chp_2PE_2loop_int_mesh_z%pnt(i)%x
!!$      z2 = z * z
!!$      zt = sqrt(1 - z2)
!!$      D1 = gA2 * (2 - z2)
!!$
!!$      DO j = 1, x_nr
!!$        w = chp_2PE_2loop_int_mesh_x%pnt(j)%w
!!$        x = chp_2PE_2loop_int_mesh_x%pnt(j)%x
!!$        y = x * zt
!!$        y2 = y * y
!!$        y_z = y / z
!!$
!!$        res = res + fact * w * zt * z * ((gA2 - 1) * y2 - D1) * (-y2 + 2*y * sqrt(z2 + y2) * log(y_z + sqrt(1 + y_z*y_z)) + C1 * (2 - y2 - z2) * (chp_N3LO_2PE_2loop_int_helper(y_z) + 5.0D0/6.0D0))
!!$      END DO
!!$
!!$      chp_2PE_2loop_int_WC_x_fact(i) = res
!!$    END DO
!!$!$OMP end parallel do
!!$
!!$    chp_2PE_2loop_int_WC_data_set = .TRUE.
!!$  END SUBROUTINE chp_calculate_2PE_2loop_int_WC_data
!!$
!!$  FUNCTION chp_N3LO_2PE_Wc_2loop_b(q2,impi) RESULT(res)
!!$    REAL(8) , INTENT(IN) :: q2
!!$    INTEGER, INTENT(IN) :: impi
!!$    REAL(8) :: res
!!$    INTEGER :: i
!!$    REAL(8) :: fact, a1
!!$    
!!$    if(chp_use_2PE_2loop_int%val == 0) then
!!$      res = 0
!!$      return
!!$    end if
!!$
!!$    call chp_calculate_2PE_2loop_int_WC_data
!!$
!!$    fact = fpi_inv**6
!!$    a1 = 4.0D0 * mpi2(impi)
!!$
!!$!$OMP parallel do schedule(static) default(none) shared(q2, impi, res, fact, a1, chp_2PE_2loop_int_mesh_z, chp_2PE_2loop_int_WC_x_fact) private(i)
!!$    DO i = 1, nof_theta_int_points
!!$      res(i) = fact * q2(i)**3 * sum(chp_2PE_2loop_int_mesh_z%pnt(:)%w * chp_2PE_2loop_int_WC_x_fact(:) / (a1 + chp_2PE_2loop_int_mesh_z%pnt(:)%x**2 * q2(i)))
!!$    END DO
!!$!$OMP end parallel do
!!$  END FUNCTION chp_N3LO_2PE_Wc_2loop_b   

! [1]Eq. D.7
  FUNCTION chp_N3LO_2PE_VLS_1loop_ci1(w, L,impi, imnuc) RESULT(res)
    REAL(8) , INTENT(IN) :: w
    REAL(8) , INTENT(IN) :: L
    INTEGER, INTENT(IN) :: impi
    INTEGER, INTENT(IN) :: imnuc
    REAL(8) :: res
    
    res = (c2*gA2*const(22)*mnuc_inv(imnuc)*(1D0/8D0))*w*w*L
  END FUNCTION chp_N3LO_2PE_VLS_1loop_ci1

  ! [1]Eq. D.13
  FUNCTION chp_N3LO_2PE_VLS_1loop_ci0(q2, w, L,impi, imnuc) RESULT(res)
    REAL(8) , INTENT(IN) :: q2
    REAL(8) , INTENT(IN) :: w
    REAL(8) , INTENT(IN) :: L
    INTEGER, INTENT(IN) :: impi
    INTEGER, INTENT(IN) :: imnuc
    REAL(8) :: res
    
    res = (gA2*gA2*const(22)*0.25D0*mnuc_inv(imnuc)*mnuc_inv(imnuc))*L*( &
              (11D0/32D0)*q2 + mpi4(impi)/(w*w))
  END FUNCTION chp_N3LO_2PE_VLS_1loop_ci0


  ! [1]Eq. D.8
  FUNCTION chp_N3LO_2PE_WLS_1loop_ci1(q2, w, L,impi, imnuc) RESULT(res)
    REAL(8) , INTENT(IN) :: q2
    REAL(8) , INTENT(IN) :: w
    REAL(8) , INTENT(IN) :: L
    INTEGER, INTENT(IN) :: impi
    INTEGER, INTENT(IN) :: imnuc
    REAL(8) :: res
    
    res = (-c4*const(22)*mnuc_inv(imnuc)*(1D0/48D0))*L*( &
              gA2*(8D0*mpi2(impi) + 5D0*q2) + w*w)
  END FUNCTION chp_N3LO_2PE_WLS_1loop_ci1
!!$
  ! [1]Eq. D.14
  FUNCTION chp_N3LO_2PE_WLS_1loop_ci0(q2, w, L,impi, imnuc) RESULT(res)
    REAL(8) , INTENT(IN) :: q2
    REAL(8) , INTENT(IN) :: w
    REAL(8) , INTENT(IN) :: L
    INTEGER, INTENT(IN) :: impi
    INTEGER, INTENT(IN) :: imnuc
    REAL(8) :: res
    
    res = (const(22)*mnuc_inv(imnuc)*mnuc_inv(imnuc)*(1D0/256D0))*L*( &
              16D0*gA2*(mpi2(impi) + 0.375D0*q2) + (4D0/3D0)*gA2*gA2*( &
                  4D0*mpi4(impi)/(w*w) - (11D0/4D0)*q2 - 9D0*mpi2(impi)) - w*w)
  END FUNCTION chp_N3LO_2PE_WLS_1loop_ci0

  ! [1]Eq. D.11
  FUNCTION chp_N3LO_2PE_VT_1loop_ci0 (k2, q2, w, L,impi, imnuc) RESULT(res)
    REAL(8) , INTENT(IN) :: k2
    REAL(8) , INTENT(IN) :: q2
    REAL(8) , INTENT(IN) :: w
    REAL(8) , INTENT(IN) :: L
    INTEGER, INTENT(IN) :: impi
    INTEGER, INTENT(IN) :: imnuc
    REAL(8) :: res
    
    res = (gA2*gA2*const(22)*mnuc_inv(imnuc)*mnuc_inv(imnuc)*(1D0/32D0))*L*( &
              k2 + 0.625D0*q2 + mpi4(impi)/(w*w))
  END FUNCTION chp_N3LO_2PE_VT_1loop_ci0

  ! [1]Eq. D.24
  FUNCTION chp_N3LO_2PE_VT_2loop_a (w, L,impi) RESULT(res)
    REAL(8) , INTENT(IN) :: w
    REAL(8) , INTENT(IN) :: L
    INTEGER, INTENT(IN) :: impi
    REAL(8) :: res
    !!! Check again d14_minus_d15 !!!!
    res = (-gA2*const(22)*(1D0/32D0)*d14_minus_d15)*w*w*L
  END FUNCTION chp_N3LO_2PE_VT_2loop_a   
  ! [1]Eq. D.2
  FUNCTION chp_N3LO_2PE_WT_1loop_ci2 (w, L,impi) RESULT(res)
    REAL(8) , INTENT(IN) :: w
    REAL(8) , INTENT(IN) :: L
    INTEGER, INTENT(IN) :: impi
    REAL(8) :: res
    
    res = (c4*c4*const(22)*(1D0/96D0))*w*w*L
  END FUNCTION chp_N3LO_2PE_WT_1loop_ci2

  ! [1]Eq. D.6
  FUNCTION chp_N3LO_2PE_WT_1loop_ci1 (q2,w, L,impi, imnuc) RESULT(res)
    REAL(8) , INTENT(IN) :: q2
    REAL(8) , INTENT(IN) :: w
    REAL(8) , INTENT(IN) :: L
    INTEGER, INTENT(IN) :: impi
    INTEGER, INTENT(IN) :: imnuc
    REAL(8) :: res
    
    res = (-c4*const(22)*mnuc_inv(imnuc)*(1D0/192D0))*L*( &
              gA2*(16D0*mpi2(impi) + 7D0*q2) - w*w)
  END FUNCTION chp_N3LO_2PE_WT_1loop_ci1

  ! [1]Eq. D.12
  FUNCTION chp_N3LO_2PE_WT_1loop_ci0 (q2, w, L, impi, imnuc) RESULT(res)
    REAL(8) , INTENT(IN) :: q2
    REAL(8) , INTENT(IN) :: w
    REAL(8) , INTENT(IN) :: L
    INTEGER, INTENT(IN) :: impi
    INTEGER, INTENT(IN) :: imnuc
    REAL(8) :: res
    
    res = (mnuc_inv(imnuc)*mnuc_inv(imnuc)*const(22)*(1D0/1536D0))*L*( &
              4D0*gA2*gA2*(7D0*mpi2(impi) + (17D0/4D0)*q2 + 4D0*mpi4(impi)/(w*w)) - &
              32D0*gA2*(mpi2(impi) + (7D0/16D0)*q2) + w*w)
  END FUNCTION chp_N3LO_2PE_WT_1loop_ci0

  ! [1]Eq. D.27
  FUNCTION chp_N3LO_2PE_WT_2loop  (w, A, impi) RESULT(res)
    REAL(8) , INTENT(IN) :: w
    REAL(8) , INTENT(IN) :: A
    INTEGER, INTENT(IN) :: impi
    REAL(8) :: res
    
    res = (gA2*gA2*const(22)*fpi_inv*fpi_inv*(1D0/2048D0))*w*w*A*( &
              w*w*A + 2D0*mpi(impi)*(1D0 + 2D0*gA2))
  END FUNCTION chp_N3LO_2PE_WT_2loop     

  ! [1]Eq. D.11
  FUNCTION chp_N3LO_2PE_Vs_1loop_ci0 (k2, q2, w, L,impi, imnuc) RESULT(res)
    REAL(8) , INTENT(IN) :: k2
    REAL(8) , INTENT(IN) :: q2
    REAL(8) , INTENT(IN) :: w
    REAL(8) , INTENT(IN) :: L
    INTEGER, INTENT(IN) :: impi
    INTEGER, INTENT(IN) :: imnuc
    REAL(8) :: res
    
    res = -q2 * chp_N3LO_2PE_VT_1loop_ci0(k2, q2, w, L, impi, imnuc)
  END FUNCTION chp_N3LO_2PE_Vs_1loop_ci0

  ! [1]Eq. D.24
  FUNCTION chp_N3LO_2PE_Vs_2loop_a   ( q2, w, L, impi) RESULT(res)
    REAL(8) , INTENT(IN) :: q2
    REAL(8) , INTENT(IN) :: w
    REAL(8) , INTENT(IN) :: L
    INTEGER, INTENT(IN) :: impi
    REAL(8) :: res
    
    res = -q2 * chp_N3LO_2PE_VT_2loop_a(w, L, impi)
  END FUNCTION chp_N3LO_2PE_Vs_2loop_a   

!!$  ! [1]Eq. D.25
!!$  FUNCTION chp_N3LO_2PE_Vs_2loop_b (q2,impi) RESULT(res)
!!$    REAL(8) , INTENT(IN) :: q2
!!$    INTEGER, INTENT(IN) :: impi
!!$    REAL(8) :: res
!!$    
!!$    res = -q2 * chp_N3LO_2PE_VT_2loop_b(q2, impi)
!!$  END FUNCTION chp_N3LO_2PE_Vs_2loop_b   

 ! [1]Eq. D.2
  FUNCTION chp_N3LO_2PE_Ws_1loop_ci2 (q2, w, L,impi) RESULT(res)
    REAL(8) , INTENT(IN) :: q2
    REAL(8) , INTENT(IN) :: w
    REAL(8) , INTENT(IN) :: L
    INTEGER, INTENT(IN) :: impi
    REAL(8) :: res
    
    res = -q2 * chp_N3LO_2PE_WT_1loop_ci2(w, L, impi)
  END FUNCTION chp_N3LO_2PE_Ws_1loop_ci2

  ! [1]Eq. D.6
  FUNCTION chp_N3LO_2PE_Ws_1loop_ci1 (q2, w, L,impi, imnuc) RESULT(res)
    REAL(8) , INTENT(IN) :: q2
    REAL(8) , INTENT(IN) :: w
    REAL(8) , INTENT(IN) :: L
    INTEGER, INTENT(IN) :: impi
    INTEGER, INTENT(IN) :: imnuc
    REAL(8) :: res
    
    res = -q2 * chp_N3LO_2PE_WT_1loop_ci1(q2, w, L, impi, imnuc)
  END FUNCTION chp_N3LO_2PE_Ws_1loop_ci1

  ! [1]Eq. D.12
  FUNCTION chp_N3LO_2PE_Ws_1loop_ci0 (q2, w, L,impi, imnuc) RESULT(res)
    REAL(8) , INTENT(IN) :: q2
    REAL(8) , INTENT(IN) :: w
    REAL(8) , INTENT(IN) :: L
    INTEGER, INTENT(IN) :: impi
    INTEGER, INTENT(IN) :: imnuc
    REAL(8) :: res
    
    res = -q2 * chp_N3LO_2PE_WT_1loop_ci0(q2, w, L, impi, imnuc)
  END FUNCTION chp_N3LO_2PE_Ws_1loop_ci0
!!$
  ! [1]Eq. D.27
  FUNCTION chp_N3LO_2PE_Ws_2loop(q2, w,A, impi) RESULT(res)
    REAL(8) , INTENT(IN) :: q2
    REAL(8) , INTENT(IN) :: w
    REAL(8) , INTENT(IN) :: A
    INTEGER, INTENT(IN) :: impi
    REAL(8) :: res
    
    res = -q2 * chp_N3LO_2PE_WT_2loop(w, A, impi)
  END FUNCTION chp_N3LO_2PE_Ws_2loop     

  ! [1]Eq. D.15
  FUNCTION chp_N3LO_2PE_VsigL_1loop_ci0(L,impi, imnuc) RESULT(res)
    REAL(8) , INTENT(IN) :: L
    INTEGER, INTENT(IN) :: impi
    INTEGER, INTENT(IN) :: imnuc
    REAL(8) :: res
    
    res = (gA2*gA2*const(22)*mnuc_inv(imnuc)*mnuc_inv(imnuc)*(1D0/32D0))*L
  END FUNCTION chp_N3LO_2PE_VsigL_1loop_ci0

!
! xxxxxxxxxx SFR 2PE 2Loop  xxxxxxxxxx
!
  ! [2]Eq. 2.19 and 2.20, the part for V_c(q)
! Due to the large derivative of the z-integrand w.r.t z close to z=1, the integral
! is done numerically up to 0.95, then the contribution from 0.95 to 1.0 as a function
! of q^2 / m_\pi^2 is calculated from a taylor expansion around z = 1.
! The relative error of the approximate expression for the integral from 0.95 to 1.0
! is smaller than 2.6e-10 for all values of q^2 / m_\pi^2
  FUNCTION chp_N3LO_2PE_Vc_2loop_SFR( q2, impi) RESULT(res)
    REAL(8) , INTENT(INOUT) :: q2
    INTEGER, INTENT(IN) :: impi
    REAL(8) :: res
    INTEGER :: i, j
    INTEGER :: z_nr
    REAL(8) :: fact
    REAL(8) :: a1, s
    REAL(8) :: s2
    REAL(8) :: z, z2, zt
! Magic numbers, I got these from a Taylor expansion calculated in Mathematica
    INTEGER, PARAMETER :: appr_order = 6
    REAL(8), PARAMETER, DIMENSION(0:appr_order) :: c = (/ 3593.92733718761d0, 5431.46429291511d0, 3420.6675650193d0, 1149.1163208075d0, 217.17174525573d0, 21.89317093303d0, 0.9197592690185d0 /)
    REAL(8) :: L
    REAL(8) :: x, m_pi_inv
    REAL(8) :: d1, d2, d3, d4, d5
    
!    if(chiral_2PE_2loop_int%val == 0) then
!      res = 0
!      return
!    end if

    z_nr = chp_2PE_2loop_int_mesh_z_VC_SFR%info%amount
    L = chp_2PE_2loop_int_mesh_z_VC_SFR%info%xmin
    
    
    fact = fpi_inv**6 * gA2*gA2 * 3.0d0 / (pi2*2.0d0**15)
    
    a1 = 4.0D0 * mpi2(impi)
    m_pi_inv = 1 / mpi(impi)
    d1 = (L-1)*32*(1+4*gA2)*mpi2(impi)*mpi2(impi);
    d2 = -8*(L-1)*(-29 + L + L*L + 4*(-11 + L + L*L)*gA2)*mpi2(impi)/3;
    d3 = 2*(L-1)*(193 + 172*gA2 + L*(L+1)*(-47 - 68*gA2 + 3*L*L*(1 + 4*gA2)))/15;
    d4 = 32*mpi(impi)*(1 + 4*gA2)*mpi2(impi);
    d5 = 64*mpi(impi)*(1 + gA2);
    
! Calculate the part from 0.95 to 1.0
      x = q2 * m_pi_inv**2
      s = c(appr_order)

      DO j = appr_order - 1, 0, -1
        s = x * s + c(j)
      END DO

      s = s * m_pi_inv**2 * (1/(4+x))**(appr_order+1)
! Calculate the rest as a numerical integration from z_low to 0.95
      DO j = 1, z_nr
        z = chp_2PE_2loop_int_mesh_z_VC_SFR%pnt(j)%x
        
        z2 = z*z
        zt = sqrt(1 - z2)
        s = s + chp_2PE_2loop_int_mesh_z_VC_SFR%pnt(j)%w * (2-z2)*(z2-8)*(z2-2)*z*0.5d0*log((1+z)/(1-z)) / (a1 + z2 * q2)
      END DO
      s = s * q2**3
      x = 0.5*sqrt(x)
      s2 = s + d1 + q2*(d2 + q2*d3) + (0.5*a1 + q2) * (d4 + d5*q2) * (atan(x) - atan(L*x)) / sqrt(q2)
      res = fact * s2 


  END FUNCTION chp_N3LO_2PE_Vc_2loop_SFR
  

! [2]Eq. 2.19 and 2.20, the part for W_C(q)
! Calculates all q-independent stuff in the z-integral (z = 2*m_\pi / \mu), i.e. it need only be done once.
  SUBROUTINE chp_calculate_2PE_2loop_int_WC_SFR_data
    INTEGER :: z_nr, x_nr, i, j
    REAL(8) :: c1, c2, cz1, cz2
    REAL(8) :: z, z2, zt
    REAL(8) :: res
    REAL(8) :: w, x

    IF(chp_2PE_2loop_int_WC_SFR_data_set) return

    z_nr = chp_2PE_2loop_int_mesh_z%info%amount
    x_nr = chp_2PE_2loop_int_mesh_x%info%amount
    c1 = -gA2**3
    c2 = gA2*gA2*(2*gA2-1)

    DO i = 1, z_nr
      res = 0
      z = chp_2PE_2loop_int_mesh_z%pnt(i)%x
      z2 = z * z
      zt = sqrt(1 - z2)
      cz1 = c1 * (2-z2)**2
      cz2 = c2 * (2-z2)*(1-z2)

      DO j = 1, x_nr
        w = chp_2PE_2loop_int_mesh_x%pnt(j)%w
        x = chp_2PE_2loop_int_mesh_x%pnt(j)%x

        res = res + w * chp_N3LO_2PE_2loop_int_helper(x*zt/z) * (cz1 + cz2*x*x)
      END DO

      chp_2PE_2loop_int_WC_SFR_x_fact(i) = res
      
   end do
  

    chp_2PE_2loop_int_WC_SFR_data_set = .TRUE.
  END SUBROUTINE chp_calculate_2PE_2loop_int_WC_SFR_data

  
! [2]Eq. 2.19 and 2.20, the part for W_c(q), LEC independent part
  FUNCTION chp_N3LO_2PE_Wc_2loop_SFR_a(q2,impi) RESULT(res)
    REAL(8) , INTENT(IN) :: q2
    INTEGER, INTENT(IN) :: impi
    REAL(8) :: res
    INTEGER :: i, j
    INTEGER :: z_nr
    REAL(8) :: fact
    REAL(8) :: a1, a2, a3, a4, b1, b2, b3
    REAL(8) :: d1
    REAL(8) :: z, z2, z4, zt
    REAL(8) :: s
    
!!$    if(chiral_2PE_2loop_int%val == 0) then
!!$      res = 0
!!$      return
!!$    end if
    
    call chp_calculate_2PE_2loop_int_WC_SFR_data
    
    z_nr = chp_2PE_2loop_int_mesh_z%info%amount
    
    fact = fpi_inv**6 / (pi2*pi2*3.0d0*2.0d0**10)
    d1 = 4.0D0 * mpi2(impi)
    a1 = -(1+2*gA2)**2/3.0d0
    a2 = 1 + 6*gA2 + 8*gA2*gA2
    a3 = -gA2*(8+15*gA2)
    a4 = (-4 + 29*gA2 + 122*gA2*gA2 + 3*gA2**3) / 15.0d0
    b1 = -(171 + 2*gA2*(1+gA2)*(327+49*gA2)) / 450.0d0
    b2 = (-73 + 1748*gA2 + 2549*gA2*gA2 + 726*gA2**3) / 450.0d0
    b3 = -(-64 + 389*gA2 + 1782*gA2*gA2 + 1093*gA2**3) / 450.0d0
    
    
    !DO i = 1, nof_theta_int_points
    s = 0
    DO j = 1, z_nr
       z = chp_2PE_2loop_int_mesh_z%pnt(j)%x
       z2 = z**2
       z4 = z2**2
       zt = sqrt(1-z2)
       s = s + chp_2PE_2loop_int_mesh_z%pnt(j)%w * z / (d1 + z2 * q2) * &
            ((a1*z4*z2 + a2*z4 + a3*z2 + a4)*log((1+zt)/z) + &
            (chp_2PE_2loop_int_WC_SFR_x_fact(j) + b1*z4 + b2*z2 + b3)*zt)
    END DO
    res = fact * q2**3 * s
    ! END DO

    
  END FUNCTION chp_N3LO_2PE_Wc_2loop_SFR_a

  ! [2]Eq. 2.19 and 2.20, the part for W_c(q), LEC dependent part
  FUNCTION chp_N3LO_2PE_Wc_2loop_SFR_b(q2,impi) RESULT(res)
    REAL(8) , INTENT(INOUT) :: q2
    INTEGER, INTENT(IN) :: impi
    REAL(8) :: res
    INTEGER :: i, j
    INTEGER :: z_nr
    REAL(8) :: fact1, fact2, fact3
    REAL(8) :: z, z2, zt
    REAL(8) :: s1, s2, s3
    REAL(8) :: d1, tmp
    
!!$    if(chiral_2PE_2loop_int%val == 0) then
!!$      res = 0
!!$      return
!!$    end if


    call chp_calculate_2PE_2loop_int_WC_SFR_data

    z_nr = chp_2PE_2loop_int_mesh_z%info%amount

    d1 = 4.0D0 * mpi2(impi)
    fact1 = -d1_plus_d2 * fpi_inv**4 / (pi2 * 96 )
    fact2 = -d3         * fpi_inv**4 / (pi2 * 480)
    fact3 =  d5         * fpi_inv**4 / (pi2 * 48 )


!    DO i = 1, nof_theta_int_points
      s1 = 0
      s2 = 0
      s3 = 0
      DO j = 1, z_nr
        z = chp_2PE_2loop_int_mesh_z%pnt(j)%x
        z2 = z**2
        zt = sqrt(1-z2)
        tmp = chp_2PE_2loop_int_mesh_z%pnt(j)%w * z * zt / (d1 + z2 * q2)
        s1 = s1 + tmp * (2-z2)*(  (gA2-1)*(1-z2) - 3*gA2*(2-z2))
        s2 = s2 + tmp * (1-z2)*(3*(gA2-1)*(1-z2) - 5*gA2*(2-z2))
        s3 = s3 + tmp * z2    *(  (gA2-1)*(1-z2) - 3*gA2*(2-z2))
      END DO
      res = q2**3 * (fact1 * s1 + fact2 * s2 + fact3 * s3)
 !   END DO

  END FUNCTION chp_N3LO_2PE_Wc_2loop_SFR_b


 ! [2]Eq. 2.19 and 2.20, the part for V_T(q), LEC dependent part
  FUNCTION chp_N3LO_2PE_VT_2loop_SFR_a ( q2,impi) RESULT(res)
    REAL(8) , INTENT(IN) :: q2
    INTEGER, INTENT(IN) :: impi
    REAL(8) :: res
    INTEGER :: i, j
    INTEGER :: z_nr
    REAL(8) :: fact, s
    REAL(8) :: a1
    REAL(8) :: z, z2, zt
    
!!$    if(chiral_2PE_2loop_int%val == 0) then
!!$      res = 0
!!$      return
!!$    end if

    z_nr = chp_2PE_2loop_int_mesh_z%info%amount

    fact = -d14_minus_d15 * fpi_inv**4 * gA2 / (pi2*32.0d0)
    a1 = 4.0D0 * mpi2(impi)


!    DO i = 1, nof_theta_int_points
      s = 0
      DO j = 1, z_nr
        z = chp_2PE_2loop_int_mesh_z%pnt(j)%x
        z2 = z*z
        zt = sqrt(1 - z2)
        s = s + chp_2PE_2loop_int_mesh_z%pnt(j)%w * z*zt*(1-z2) / (a1 + z2 * q2)
      END DO
      res = fact * q2**2 * s
 !   END DO

  END FUNCTION chp_N3LO_2PE_VT_2loop_SFR_a

! [1]Eq. D.25 (LEC independent part), identical with [2]Eq. 2.19 and 2.20, part for V_TS (LEC independent part).
! Only the integration limits differ in [1] (DR) and [2] (SFR)
! Calculates all q-independent stuff in the z-integral (z = 2*m_\pi / \mu) for the LEC independent part, 
! i.e. it need only be done once.

  SUBROUTINE chp_calculate_2PE_2loop_int_VTS_data
    INTEGER :: z_nr, x_nr, i, j
    REAL(8) :: fact, z_fact
    REAL(8) :: res, z, zt
    REAL(8) :: w, x

    IF(chp_2PE_2loop_int_VTS_data_set) return

    z_nr = chp_2PE_2loop_int_mesh_z%info%amount
    x_nr = chp_2PE_2loop_int_mesh_x%info%amount
    fact = -1.0D0 / (1024.0D0 * pi2 * pi2)


    DO i = 1, z_nr
      res = 0
      z = chp_2PE_2loop_int_mesh_z%pnt(i)%x
      zt = sqrt(1 - z*z)
      z_fact = fact * chp_2PE_2loop_int_mesh_z%pnt(i)%w * z * zt * (1-z*z)

      DO j = 1, x_nr
        w = chp_2PE_2loop_int_mesh_x%pnt(j)%w
        x = chp_2PE_2loop_int_mesh_x%pnt(j)%x

        res = res + z_fact * w * (1 - x*x) * (chp_N3LO_2PE_2loop_int_helper(x*zt/z) - 1.0D0/6.0D0)
      END DO

      chp_2PE_2loop_int_VTS_x_fact(i) = res
   END DO


    chp_2PE_2loop_int_VTS_data_set = .TRUE.
  END SUBROUTINE chp_calculate_2PE_2loop_int_VTS_data
 
! [1]Eq. D.25 (LEC independent part), identical with [2]Eq. 2.19 and 2.20, part for V_T(q) (LEC independent part).
  FUNCTION chp_N3LO_2PE_VT_2loop_b(q2,impi) RESULT(res)
    REAL(8) , INTENT(IN) :: q2
    INTEGER, INTENT(IN) :: impi
    REAL(8) :: res
    INTEGER :: i, j
    REAL(8) :: fact
    REAL(8) :: a1
    INTEGER :: z_nr
    
!!$    if(chiral_2PE_2loop_int%val == 0) then
!!$      res = 0
!!$      return
!!$    end if

    call chp_calculate_2PE_2loop_int_VTS_data

    z_nr = chp_2PE_2loop_int_mesh_z%info%amount
    fact = (gA * fpi_inv)**6
    a1 = 4.0D0 * mpi2(impi)


!    DO i = 1, nof_theta_int_points
      res = 0
      DO j = 1, z_nr
        res = res + chp_2PE_2loop_int_VTS_x_fact(j) / (a1 + chp_2PE_2loop_int_mesh_z%pnt(j)%x**2 * q2)
      END DO
      res = res * fact * q2**2
 !   END DO

  END FUNCTION chp_N3LO_2PE_VT_2loop_b   
  
  ! [2]Eq. 2.19 and 2.20, par for W_T(q)
  FUNCTION chp_N3LO_2PE_WT_2loop_SFR(q2, impi) RESULT(res)
    REAL(8) , INTENT(IN) :: q2
    INTEGER, INTENT(IN) :: impi
    REAL(8) :: res
    INTEGER :: i, j
    INTEGER :: z_nr
    REAL(8) :: fact
    REAL(8) :: a1
    REAL(8) :: s
    REAL(8) :: z, z2, zt
    
!!$    if(chiral_2PE_2loop_int%val == 0) then
!!$      res = 0
!!$      return
!!$    end if

    z_nr = chp_2PE_2loop_int_mesh_z%info%amount

    fact = -fpi_inv**6 * gA2 * gA2 / (pi2*4096.0d0)
    a1 = 4.0D0 * mpi2(impi)

!    DO i = 1, nof_theta_int_points
      s = 0
      DO j = 1, z_nr
        z = chp_2PE_2loop_int_mesh_z%pnt(j)%x
        z2 = z*z
        zt = sqrt(1 - z2)
        s = s + chp_2PE_2loop_int_mesh_z%pnt(j)%w * (1-z2) * ((z2-1)*z*0.5d0*log((1+z)/(1-z)) + (1+2*gA2)*z2) / (a1 + z2 * q2)
      END DO
      res = fact * q2**2 * s
 !   END DO

  END FUNCTION chp_N3LO_2PE_WT_2loop_SFR

! [2]Eq. 2.19 and 2.20, the part for V_S(q), LEC dependent part
  FUNCTION chp_N3LO_2PE_Vs_2loop_SFR_a(q2, impi) RESULT(res)
    REAL(8) , INTENT(IN) :: q2
    INTEGER, INTENT(IN) :: impi
    REAL(8) :: res
    
    res = -q2 * chp_N3LO_2PE_VT_2loop_SFR_a(q2, impi)
    
  END FUNCTION chp_N3LO_2PE_Vs_2loop_SFR_a
  
! [1]Eq. D.25 (LEC independent part), identical with [2]Eq. 2.19 and 2.20, part for V_S(q) (LEC independent part).
  FUNCTION chp_N3LO_2PE_Vs_2loop_b(q2, impi) RESULT(res)
    REAL(8) , INTENT(IN) :: q2
    INTEGER, INTENT(IN) :: impi
    REAL(8) :: res
    
    res = -q2 * chp_N3LO_2PE_VT_2loop_b(q2, impi)
    
  END FUNCTION chp_N3LO_2PE_Vs_2loop_b   
  
  
! [2]Eq. 2.19 and 2.20, par for W_S(q)
  FUNCTION chp_N3LO_2PE_Ws_2loop_SFR(    q2, impi) RESULT(res)
    REAL(8) , INTENT(IN) :: q2
    INTEGER, INTENT(IN) :: impi
    REAL(8) :: res
    
    res = -q2 * chp_N3LO_2PE_WT_2loop_SFR(q2, impi)
  END FUNCTION chp_N3LO_2PE_Ws_2loop_SFR


  
! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!
!  End of the N3LO 2PE 2LOOP DIAGRAMS
!
! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXxXXX
  
  FUNCTION chp_sigma_dot_sigma_mtx(ms1,ms2,ms3,ms4) RESULT(res)
    
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

    !if ( abs( aimag(res1)) > 1.e-10 ) write(6,*)  'yep', res1
    res = res1
    
    
  END FUNCTION chp_sigma_dot_sigma_mtx

  FUNCTION chp_sigma_dot_sigma(ms1,ms2,ms3,ms4) RESULT(res)
    
    IMPLICIT NONE 
    INTEGER, INTENT(IN) :: ms1,ms2,ms3,ms4
    REAL*8 :: res, spin 
    COMPLEX*16 :: res1, jx1, jx2, jy1, jy2, jz1, jz2 
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
            
    res1 = jx1*jx2 + jy1*jy2 + jz1*jz2
    res = res1
    
  END FUNCTION chp_sigma_dot_sigma

  FUNCTION chp_tau_dot_tau_mtx(ms1,ms2,ms3,ms4) RESULT(res)
    
    IMPLICIT NONE 
    INTEGER, INTENT(IN) :: ms1,ms2,ms3,ms4
    REAL*8 :: res, delta
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
    
    !if ( abs( aimag(res1)) > 1.e-10 ) write(6,*)  'yep', res1
    
  END FUNCTION chp_tau_dot_tau_mtx
  
  
  FUNCTION chp_tau_dot_tau(ms1,ms2,ms3,ms4) RESULT(res)
    
    IMPLICIT NONE 
    INTEGER, INTENT(IN) :: ms1,ms2,ms3,ms4
    REAL*8 :: res, spin 
    COMPLEX*16 :: res1, jx1, jx2, jy1, jy2, jz1, jz2 
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
            
    res1 = jx1*jx2 + jy1*jy2 + jz1*jz2
    res = res1
    
  END FUNCTION chp_tau_dot_tau
  
  FUNCTION chp_sigma_dot_q_mtx(ms1,ms3,q) RESULT(res)

    IMPLICIT NONE 
    INTEGER, INTENT(IN) :: ms1, ms3
    REAL*8, INTENT(IN) :: q(3) 
    COMPLEX*16 :: res
    COMPLEX*16 :: res1, res2 
    COMPLEX*16 :: chi1(2), chi3(2), mat(2,2) 
    INTEGER :: i1 
    
    chi1 = 0.d0; chi3 = 0.d0
    i1 = nint(1.5-0.5*ms1) 
    chi1(i1) = 1.d0 
    i1 = nint(1.5-0.5*ms3) 
    chi3(i1) = 1.d0 
    
    mat = q(1)*sigma_x+q(2)*sigma_y+q(3)*sigma_z 
    res1 = dot_product(chi1, matmul( mat, chi3))
    res = res1 
    
  end FUNCTION chp_sigma_dot_q_mtx

  
  FUNCTION chp_sigma1_dot_q_sigma2_dot_q_mtx(ms1,ms2,ms3,ms4,q) RESULT(res)
    
    IMPLICIT NONE 
    INTEGER, INTENT(IN) :: ms1,ms2,ms3,ms4
    REAL*8, INTENT(IN) :: q(3) 
    COMPLEX*16 :: res
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
    !if ( abs( aimag(res1*res2) ) > 1.e-10 ) write(6,*) res1*res2 
    

  END FUNCTION chp_sigma1_dot_q_sigma2_dot_q_mtx

  FUNCTION chp_sigma1_dot_q_sigma2_dot_q(ms1,ms2,ms3,ms4,q) RESULT(res)

    IMPLICIT NONE 
    INTEGER, INTENT(IN) :: ms1,ms2,ms3,ms4
    REAL*8, INTENT(IN) :: q(3) 
    COMPLEX*16 :: res
    COMPLEX*16 :: res1,res2, jx1, jx2, jy1, jy2, jz1, jz2 
    REAL*8 :: m1,m2,m3,m4, spin 
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
            
    res1 = jx1*q(1) + jy1*q(2) + jz1*q(3)
    res2 = jx2*q(1) + jy2*q(2) + jz2*q(3)
    res = res1*res2 
    
  end FUNCTION chp_sigma1_dot_q_sigma2_dot_q

  FUNCTION chp_spin_dot_qxk_mtx(ms1,ms2,ms3,ms4,qxk) RESULT(res)

    IMPLICIT NONE 
    INTEGER, INTENT(IN) :: ms1,ms2,ms3,ms4
    REAL*8, INTENT(IN) :: qxk(3)
    real*8 :: delta 
    COMPLEX*16 :: res
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
    
    mat =-dcmplx(0.d0,1.d0) * 0.5D0*(qxk(1)*sigma_x+qxk(2)*sigma_y+qxk(3)*sigma_z)
    
    res1 = dot_product(chi1, matmul( mat, chi3))*delta(ms2,ms4)
    res2 = dot_product(chi2, matmul( mat, chi4))*delta(ms1,ms3)
    
    res = (res1+res2)
    
    
  END FUNCTION chp_spin_dot_qxk_mtx
  
  
  FUNCTION chp_spin_dot_qxk(ms1,ms2,ms3,ms4,qxk) RESULT(res)

    IMPLICIT NONE 

    INTEGER, INTENT(IN) :: ms1,ms2,ms3,ms4
    REAL*8, INTENT(IN) :: qxk(3) 
    REAL*8 ::  spin, delta, d1, d2
    COMPLEX*16 :: res, res1,res2, jx1, jx2, jy1, jy2, jz1, jz2 
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
    
    d1 = delta(ms1,ms3)
    d2 = delta(ms2,ms4)
    res1 = -dcmplx(0.d0,0.5d0) * ( (jx1*d2+jx2*d1)*qxk(1) + (jy1*d2+jy2*d1)*qxk(2) + (jz1*d2+jz2*d1)*qxk(3) )
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
    !exponent = (q2/lambda**2)**4    
    exponent = (q2/lambda**2)**2    
    !exponent = (q2/400.d0**2)**2    
    !res =  exp(-exponent)
    res = 1.d0 
    
  END FUNCTION freg_3nflocal


  ! non-local 3nf regulator, eq 3 in http://arxiv.org/abs/1304.2212
  ! pfinal : final momentum
  ! pinit  : initial momentum
  ! n      : cutoff order 
  ! LAMBDA is accessed from the chp constant chp_lambda
  FUNCTION freg_3nfnonlocal(k1,k2,k3) RESULT(res)
    
    use constants, only : hbarc 
    REAL*8 , INTENT(IN) :: k1(3), k2(3), k3(3)
    real*8 :: p(3), q(3) 
    REAL*8 :: res,exponent, q2, llambda, cutoff   
    INTEGER :: nexp 
    
    p = 0.5d0*(k1-k2)
    q = (2./3.)*(k3-0.5d0*(k1+k2))
    
    cutoff = 394.d0
    llambda = 4.d0 * cutoff**2 
    ! use power nexp=3 to match with Kyle's 3NF non-local HOSC matrix elements 
    nexp = 4
    res = 0.0D0
    q2 = 4.d0*dot_product(p,p)+ 3.d0*dot_product(q,q)
    exponent = (q2/llambda)**nexp
    res = dexp(-exponent)
    
    
  END FUNCTION freg_3nfnonlocal

  
  FUNCTION chp_3nf_Dterm(ms1,ms2,ms3,ms4,q) RESULT(res)
    
    IMPLICIT NONE 
    INTEGER, INTENT(IN) :: ms1,ms2,ms3,ms4
    REAL*8, INTENT(IN) :: q(3) 
    REAL*8 :: q2
    COMPLEX*16 :: res, res1
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
    REAL*8 :: q2
    COMPLEX*16 :: res, res1 
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
  !
  ! New function Aug. 8, 2013
  ! 
  FUNCTION chp_sigma_dot_qxk_mtx(ms1,ms3,q,k) RESULT(res)
    
    IMPLICIT NONE 
    INTEGER, INTENT(IN) :: ms1,ms3
    REAL*8, INTENT(IN) :: q(3), k(3)
    REAL*8 :: qxk(3)
    real*8 :: delta 
    COMPLEX*16 :: res
    COMPLEX*16 :: res1
    COMPLEX*16 :: chi1(2), chi3(2), mat(2,2) 
    INTEGER :: i1 
    
    qxk(1) = q(2)*k(3)-q(3)*k(2) 
    qxk(2) = q(3)*k(1)-q(1)*k(3) 
    qxk(3) = q(1)*k(2)-q(2)*k(1) 
    
    
    chi1 = 0.d0; chi3 = 0.d0
    i1 = nint(1.5-0.5*ms1) 
    chi1(i1) = 1.d0 
    i1 = nint(1.5-0.5*ms3) 
    chi3(i1) = 1.d0 
    
    mat = (qxk(1)*sigma_x+qxk(2)*sigma_y+qxk(3)*sigma_z)
    
    res1 = dot_product(chi1, matmul( mat, chi3))
    res = (res1)
    
    
  END FUNCTION chp_sigma_dot_qxk_mtx
  
  
  FUNCTION chp_sigma1_dot_qXq(ms1,ms2,q2,q3) RESULT(res)

    IMPLICIT NONE 
    INTEGER, INTENT(IN) :: ms1,ms2
    REAL*8, INTENT(IN) :: q2(3), q3(3) 
    REAL*8 :: qxq(3)
    COMPLEX*16 :: chi1(2), chi2(2), mat(2,2), res 
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
    REAL*8 :: qxq(3)
    COMPLEX*16 :: chi1(2), chi2(2), chi3(2), chi4(2), chi5(2), chi6(2), & 
         mat(2,2), facx, facy, facz,res   
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
    
    ! andreas: t2.t3 = (t2yt3z-t2zt3y, t2zt33-t2xt3z, t2xt3y-t2yt3x)
    ! facx = chi2[sigma_y]chi5 * chi3[sigma_z]chi6
    ! etc
    !
    ! ... therefore I have changed chi3 to chi2 in the first factor
    ! of both terms in the entire vector

    !facx = dot_product(chi3,matmul(sigma_y,chi5))*dot_product(chi3,matmul(sigma_z,chi6))  - &
    !     dot_product(chi3,matmul(sigma_z,chi5))*dot_product(chi3,matmul(sigma_y,chi6))  
    !
    !facy = dot_product(chi3,matmul(sigma_z,chi5))*dot_product(chi3,matmul(sigma_x,chi6))  - &
    !     dot_product(chi3,matmul(sigma_x,chi5))*dot_product(chi3,matmul(sigma_z,chi6))  
    !
    !facz = dot_product(chi3,matmul(sigma_x,chi5))*dot_product(chi3,matmul(sigma_y,chi6))  - &
    !     dot_product(chi3,matmul(sigma_y,chi5))*dot_product(chi3,matmul(sigma_x,chi6))  

    facx = dot_product(chi2,matmul(sigma_y,chi5))*dot_product(chi3,matmul(sigma_z,chi6))  - &
         dot_product(chi2,matmul(sigma_z,chi5))*dot_product(chi3,matmul(sigma_y,chi6))  

    facy = dot_product(chi2,matmul(sigma_z,chi5))*dot_product(chi3,matmul(sigma_x,chi6))  - &
         dot_product(chi2,matmul(sigma_x,chi5))*dot_product(chi3,matmul(sigma_z,chi6))  

    facz = dot_product(chi2,matmul(sigma_x,chi5))*dot_product(chi3,matmul(sigma_y,chi6))  - &
         dot_product(chi2,matmul(sigma_y,chi5))*dot_product(chi3,matmul(sigma_x,chi6))  
    
    ! andreas: t1.fac = chi1[t1.fac]chi4
    !
    !
    
    mat = facx*sigma_x+facy*sigma_y+facz*sigma_z 
    res = dot_product(chi1, matmul( mat, chi4))
    
  END FUNCTION chp_tau1_dot_tauXtau
  
  !NEW FUNCTION: we need (tau x tau).tau and not what is coded above [namely tau . (tau x tau)]
  
  FUNCTION chp_tau1Xtau2_dot_tau3(t1,t2,t3,t4,t5,t6) RESULT(res)

    IMPLICIT NONE 
    INTEGER, INTENT(IN) :: t1,t2,t3,t4,t5,t6
    complex*16 :: res
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
    
    facx = dot_product(chi1,matmul(sigma_y,chi4))*dot_product(chi2,matmul(sigma_z,chi5))  - &
         dot_product(chi1,matmul(sigma_z,chi4))*dot_product(chi2,matmul(sigma_y,chi5))  

    facy = dot_product(chi1,matmul(sigma_z,chi4))*dot_product(chi2,matmul(sigma_x,chi5))  - &
         dot_product(chi1,matmul(sigma_x,chi4))*dot_product(chi2,matmul(sigma_z,chi5))  

    facz = dot_product(chi1,matmul(sigma_x,chi4))*dot_product(chi2,matmul(sigma_y,chi5))  - &
         dot_product(chi1,matmul(sigma_y,chi4))*dot_product(chi2,matmul(sigma_x,chi5))  


    mat = facx*sigma_x+facy*sigma_y+facz*sigma_z 
    res = dot_product(chi3, matmul( mat, chi6))
        
  END FUNCTION chp_tau1Xtau2_dot_tau3

  

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
  TYPE(chp_real_type), PRIVATE :: chp_regcut, chp_regcut_nnlo

  
  
contains 
  
  complex*16 function chiral_pot(p,q,r,s)
    USE single_particle_orbits
    USE constants
    use chiral_constants
    use ang_mom_functions, only : tjs 
    
    implicit none 
    INTEGER, INTENT(IN) :: p,q,r,s 
    INTEGER :: m1,m2,m3,m4, spin, iph, t1,t2,t3,t4, Tiso
    INTEGER :: nx1, ny1, nz1, nx2, ny2, nz2, nx3, ny3, nz3, nx4, ny4, nz4, mt1, mt2, mt3, mt4 
    
    REAL*8 :: k1(3), k2(3), k3(3), k4(3), kmean(3) 
    REAL*8 :: qtrans(3), prel(3), pprel(3), qxk(3)
    REAL*8 :: q2, p2, kmean2, qabs, pp2, cg1, cg2 
    REAL*8 :: delta, nucleon_mass, relativity_factor
    
    COMPLEX*16 :: vlo, vnlo, vnnlo, vn3lo, vdir
    COMPLEX*16 :: cont_lo, cont_nlo, cont_n3lo, v1pe_tiso(0:1)
    ! LIST OF OPERATOR STRUCTURES
    COMPLEX*16 :: Vc
    COMPLEX*16 :: Wc
    COMPLEX*16 :: Vs
    COMPLEX*16 :: Ws
    COMPLEX*16 :: VLS
    COMPLEX*16 :: WLS
    COMPLEX*16 :: VsigL
    COMPLEX*16 :: WsigL
    COMPLEX*16 :: VT
    COMPLEX*16 :: WT
    COMPLEX*16 :: Vsigk
    COMPLEX*16 :: Wsigk
    COMPLEX*16 :: term_CS
    COMPLEX*16 :: term_CT
    COMPLEX*16 :: term_C(1:7)
    COMPLEX*16 :: term_D(1:15)  !!! ADDED will hold N3LO contact terms -> cont_n3lo
    REAL*8 :: loop_w
    REAL*8 :: loop_wtilde2
    REAL*8 :: loop_s
    REAL*8 :: loop_L
    REAL*8 :: loop_A
    real*8 :: sigma_dot_sigma,qxk2
    INTEGER :: imnuc_2PE
    
    ! specify chiral order
    chp_chiral_order%val = N3LO 
    
    ! See Eq. 4.7 p. 17
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
  
    term_C = 0.0D0
    term_D = 0.0D0  !! initialize N3LO contact terms to zero

    vlo = 0.0D0 ; vnlo = 0.0D0 ; vnnlo = 0.0D0 ; vn3lo = 0.0D0
    
    cont_lo = 0.0D0 ; cont_nlo = 0.0D0 ; cont_n3lo = 0.0D0

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
    ! conservation of isospin 
    !
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
    
    IF (chp_chiral_order%val == LO ) THEN
       
       ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
       !
       ! LEADING ORDER 
       ! 
       ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
       
       chp_regcut%val = 2.0D0
 
!NEW
       ! relativistic 1/m^2 correction to static 1PE
       ! (1-(p^2+pp^2)/2mN^2 + ...)
       ! multiplied inside chp_one_pion_exchange
       ! Formally, this correction should only enter
       ! the static 1PE at N3LO.
       !relcorr_1PE = (1.0D0-((pp2 + p2)/(2.0D0*mnuc(0)**2)))
       
       !
       ! leading order one-pion exchange 
       ! 
       WT = chp_sigma1_dot_q_sigma2_dot_q_mtx(m1,m2,m3,m4,qtrans)*&
            chp_one_pion_exchange(q2, 2)*chp_tau_dot_tau_mtx(t1,t2,t3,t4)
    
       
       vlo = WT
       
       
       !
       ! leading order contacts 
       !
       term_CS = CS(0)*delta(m1,m3)*delta(m2,m4)*delta(t1,t3)*delta(t2,t4)
       term_CT = CT(0)*chp_sigma_dot_sigma_mtx(m1,m2,m3,m4)*delta(t1,t3)*delta(t2,t4) 
       cont_lo = term_CS + term_CT
       
    END IF
    
    IF (chp_chiral_order%val >= NLO) THEN
       
       ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
       !
       ! NEXT TO LEADING ORDER 
       ! 
       ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
       
       chp_regcut%val = 2.0D0
       
       !
       ! NLO pion-exchange 
       ! 
       ! irreducible, or non-polynomial, NLO two-pion exchanges
       
       loop_w = chp_NLO_two_pion_exchange_loop_w(q2,2)
       loop_s = chp_NLO_sfr_two_pion_exchange_loop_s(2)
       loop_L = chp_NLO_sfr_two_pion_exchange_loop_L(qabs, q2, loop_w, loop_s, 2)
       
       ! LO CIB one-pion exchanges
       ! [1] Eq. 4.77-4.79
       
       !
       ! Note: the LO contribution is modified at NLO
       ! to include CIB effects. Therefore, should you go beyond
       ! LO, make sure to not include the vlo contribution 
       ! in the sum defining vdir.
       !
       
       !vlo = chp_sigma1_dot_q_sigma2_dot_q_mtx(m1,m2,m3,m4,qtrans)*chp_tau_dot_tau_mtx(t1,t2,t3,t4)*chp_one_pion_exchange(q2, 2)
       
       WT = 0.d0 
       IF (t1 + t2 == 0) THEN 
          
          WT = 0.d0 
          DO Tiso = 0, 2, 2
             cg1 = tjs(1,1,Tiso,t1,t2,-(t1+t2))*iph( (t1+t2)/2 )*sqrt(Tiso+1.d0)
             cg2 = tjs(1,1,Tiso,t3,t4,-(t3+t4))*iph( (t3+t4)/2 )*sqrt(Tiso+1.d0)
             
             WT = WT + cg1*cg2*(-chp_one_pion_exchange(q2, 0) + iph(Tiso/2+1)*2.0D0*chp_one_pion_exchange(q2,1)) 
          end DO
          
          
       END IF
       IF (t1 + t2 /= 0) THEN
          WT = chp_one_pion_exchange(q2, 0) 
       END IF
       
!!$       IF (t1+t2 == 0) THEN
!!$          WT = -1.0D0*chp_one_pion_exchange(q2, 0) + &
!!$               chp_minus_power(CHN%T+1)*2.0D0*chp_one_pion_exchange(q2, 1)
!!$       ELSE
!!$          WT = chp_one_pion_exchange(q2, relcorr_1PE, 0)
!!$       END IF
!!$    ELSE
!!$       ! static OPE contribution
!!$       ! using average pion mass only
!!$       
!!$       WT = iso(T)*chp_one_pion_exchange(q2, relcorr_1PE, 2)
!!$    END IF
       
       
       vlo = WT*chp_sigma1_dot_q_sigma2_dot_q_mtx(m1,m2,m3,m4,qtrans)
       
       
       !
       ! leading order CIB contacts 
       !
       term_CS = CS((t1+t2)/2)*delta(m1,m3)*delta(m2,m4)*delta(t1,t3)*delta(t2,t4)
       term_CT = CT((t1+t2)/2) * chp_sigma_dot_sigma_mtx(m1,m2,m3,m4) *delta(t1,t3)*delta(t2,t4) 
       cont_lo = term_CS + term_CT
       !cont_lo = 4.*cont_lo
       
!!$       Wc = chp_tau_dot_tau_mtx(t1,t2,t3,t4)*delta(m1,m3)*delta(m2,m4)* & 
!!$            chp_NLO_two_pion_exchange_Wc(q2, loop_L, loop_w, 2)
!!$       
!!$       Vs = chp_NLO_two_pion_exchange_Vs(q2, loop_L, 2)* &
!!$            chp_sigma_dot_sigma_mtx(m1,m2,m3,m4)*delta(t1,t3)*delta(t2,t4) 
!!$       
!!$       VT = chp_NLO_two_pion_exchange_VT(loop_L, 2)* & 
!!$            chp_sigma1_dot_q_sigma2_dot_q_mtx(m1,m2,m3,m4,qtrans)*delta(t1,t3)*delta(t2,t4) 

       Wc = chp_two_pion_exchange_1loop_0_Wc(q2, loop_L, loop_w, 2)*&
            chp_tau_dot_tau_mtx(t1,t2,t3,t4)*delta(m1,m3)*delta(m2,m4)
       Vs = chp_two_pion_exchange_1loop_0_Vs(q2, loop_L, 2)* &                                                
            chp_sigma_dot_sigma_mtx(m1,m2,m3,m4)*delta(t1,t3)*delta(t2,t4)
       VT = chp_two_pion_exchange_1loop_0_VT(loop_L, 2)*&
            chp_sigma1_dot_q_sigma2_dot_q_mtx(m1,m2,m3,m4,qtrans)*delta(t1,t3)*delta(t2,t4) 
       
       vnlo = ( WT + WC + VS + VT ) 
       
       !
       ! next-to-leading order contacts 
       !
       term_C(1) = cnlo(1)*q2*delta(m1,m3)*delta(m2,m4)*delta(t1,t3)*delta(t2,t4) 
       
       term_C(2) = cnlo(2)*kmean2*delta(m1,m3)*delta(m2,m4)*delta(t1,t3)*delta(t2,t4)
       
       term_C(3) = cnlo(3)*q2*chp_sigma_dot_sigma_mtx(m1,m2,m3,m4)*delta(t1,t3)*delta(t2,t4)
       
       term_C(4) = cnlo(4)*kmean2*chp_sigma_dot_sigma_mtx(m1,m2,m3,m4)*delta(t1,t3)*delta(t2,t4)

       term_C(5) = cnlo(5)*chp_spin_dot_qxk_mtx(m1,m2,m3,m4,qxk)*delta(t1,t3)*delta(t2,t4)
       
       term_C(6) = cnlo(6)*chp_sigma1_dot_q_sigma2_dot_q_mtx(m1,m2,m3,m4,qtrans)*delta(t1,t3)*delta(t2,t4)

       term_C(7) = cnlo(7)*chp_sigma1_dot_q_sigma2_dot_q_mtx(m1,m2,m3,m4,kmean)*delta(t1,t3)*delta(t2,t4) 
       
       cont_nlo = SUM(term_C)
       
    END IF
    
    IF (chp_chiral_order%val >= NNLO) THEN
       
       ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
       !
       ! NEXT TO NEXT TO LEADING ORDER 
       ! 
       ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
       
       chp_regcut_nnlo%val = 3.0D0
       
       ! NNLO loop functions
       
       loop_wtilde2 = chp_NNLO_two_pion_exchange_loop_wtilde2(q2,2)
       loop_A = chp_NNLO_sfr_two_pion_exchange_loop_A(qabs, q2, 2)
       
       ! chp_itope type is set in the header of this function
!!$       Vc  = chp_NNLO_two_pion_exchange_Vc (q2, loop_w, loop_A, loop_wtilde2, 2)*&
!!$            delta(m1,m3)*delta(m2,m4)*delta(t1,t3)*delta(t2,t4)
!!$       
!!$       Wc  = chp_NNLO_two_pion_exchange_Wc (q2, loop_w, loop_A, loop_wtilde2, 2)*&
!!$            chp_tau_dot_tau_mtx(t1,t2,t3,t4)*delta(m1,m3)*delta(m2,m4)
!!$       
!!$       VLS = chp_NNLO_two_pion_exchange_VLS(loop_A, loop_wtilde2, 2)*&
!!$            chp_spin_dot_qxk_mtx(m1,m2,m3,m4,qxk)*delta(t1,t3)*delta(t2,t4)
!!$       
!!$       WLS = chp_NNLO_two_pion_exchange_WLS(loop_w, loop_A, 2)*&
!!$            chp_spin_dot_qxk_mtx(m1,m2,m3,m4,qxk)*chp_tau_dot_tau_mtx(t1,t2,t3,t4)
!!$       
!!$       VT  = chp_NNLO_two_pion_exchange_VT (loop_w, loop_A, loop_wtilde2, 2)*&
!!$            chp_sigma1_dot_q_sigma2_dot_q_mtx(m1,m2,m3,m4,qtrans)*delta(t1,t3)*delta(t2,t4) 
!!$       
!!$       WT  = chp_NNLO_two_pion_exchange_WT (q2, loop_w, loop_A, 2)*&
!!$            chp_sigma1_dot_q_sigma2_dot_q_mtx(m1,m2,m3,m4,qtrans)*chp_tau_dot_tau_mtx(t1,t2,t3,t4)
!!$       
!!$       Vs  = chp_NNLO_two_pion_exchange_Vs (q2, loop_w, loop_A, loop_wtilde2, 2)*&
!!$            chp_sigma_dot_sigma_mtx(m1,m2,m3,m4)*delta(t1,t3)*delta(t2,t4) 
!!$       
!!$       Ws  = chp_NNLO_two_pion_exchange_Ws (q2, loop_w, loop_A, 2)*&
!!$            chp_sigma_dot_sigma_mtx(m1,m2,m3,m4)*chp_tau_dot_tau_mtx(t1,t2,t3,t4)
       
       Vc  = (chp_two_pion_exchange_1loop_d_Vc (q2,loop_A, loop_wtilde2, 2)+&
             chp_two_pion_exchange_1loop_r_Vc (q2, loop_w, loop_A, loop_wtilde2, 2, 0))*&
             delta(m1,m3)*delta(m2,m4)*delta(t1,t3)*delta(t2,t4)
       
       WT  = (chp_two_pion_exchange_1loop_d_WT (    loop_w, loop_A )+&
            chp_two_pion_exchange_1loop_r_WT (q2, loop_w, loop_A, 2, 0))*&
             chp_sigma1_dot_q_sigma2_dot_q_mtx(m1,m2,m3,m4,qtrans)*chp_tau_dot_tau_mtx(t1,t2,t3,t4)
            
       Ws  = (chp_two_pion_exchange_1loop_d_Ws (q2, loop_w, loop_A )+&
            chp_two_pion_exchange_1loop_r_Ws (q2, loop_w, loop_A, 2, 0))*&
             chp_sigma_dot_sigma_mtx(m1,m2,m3,m4)*chp_tau_dot_tau_mtx(t1,t2,t3,t4)
      
       Wc  = chp_two_pion_exchange_1loop_r_Wc (q2, loop_w, loop_A, loop_wtilde2, 2, 0)*&
             chp_tau_dot_tau_mtx(t1,t2,t3,t4)*delta(m1,m3)*delta(m2,m4)
            
       VLS = chp_two_pion_exchange_1loop_r_VLS( loop_A, loop_wtilde2, 0)*&
            chp_spin_dot_qxk_mtx(m1,m2,m3,m4,qxk)*delta(t1,t3)*delta(t2,t4)

       WLS = chp_two_pion_exchange_1loop_r_WLS( loop_w, loop_A , 0)*&
            chp_spin_dot_qxk_mtx(m1,m2,m3,m4,qxk)*chp_tau_dot_tau_mtx(t1,t2,t3,t4)

       VT  = chp_two_pion_exchange_1loop_r_VT (q2, loop_w, loop_A, loop_wtilde2, 2, 0)*&
            chp_sigma1_dot_q_sigma2_dot_q_mtx(m1,m2,m3,m4,qtrans)*delta(t1,t3)*delta(t2,t4) 

       Vs  = chp_two_pion_exchange_1loop_r_Vs (q2, loop_w, loop_A, loop_wtilde2, 2, 0)*&
            chp_sigma_dot_sigma_mtx(m1,m2,m3,m4)*delta(t1,t3)*delta(t2,t4) 
     
       
       vnnlo = (Vc + Wc + VLS + WLS + VT + WT + Vs + Ws)
       
    END IF
    
    IF (chp_chiral_order%val >= N3LO) THEN
       
       ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
       !
       ! NEXT TO NEXT TO NEXT TO LEADING ORDER 
       ! 
       ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
       
       chp_regcut_nnlo%val = 3.0D0
       
       vn3lo = 0.0D0
       cont_n3lo = 0.0D0
       ! Note that N3LO(EM) uses isospin dependent nucleon mass. For now we use the average mass.
       imnuc_2PE = 0
       
       Vc  = (chp_N3LO_2PE_Vc_1loop_ci2( loop_w, loop_L, loop_wtilde2, 2) + &
            chp_N3LO_2PE_Vc_2loop_SFR(q2,2))*& 
            delta(m1,m3)*delta(m2,m4)*delta(t1,t3)*delta(t2,t4)
       
       
!!$       !! XXXXXXX  New Stuff XXXXXXX  !!!!
       Wc  = chp_tau_dot_tau_mtx(t1,t2,t3,t4)*delta(m1,m3)*delta(m2,m4) *&
            ( chp_N3LO_2PE_Wc_2loop_SFR_a(q2,2) +&
            chp_N3LO_2PE_Wc_2loop_SFR_b(q2,2) )
       
       VLS = 0.d0
       WLS = 0.d0
         
       VT  = chp_sigma1_dot_q_sigma2_dot_q_mtx(m1,m2,m3,m4,qtrans)*delta(t1,t3)*delta(t2,t4)*&
            (chp_N3LO_2PE_VT_2loop_SFR_a(q2,2)+ chp_N3LO_2PE_VT_2loop_b(q2,2))
       
       WT  = chp_sigma1_dot_q_sigma2_dot_q_mtx(m1,m2,m3,m4,qtrans)*chp_tau_dot_tau_mtx(t1,t2,t3,t4)*&
            (chp_N3LO_2PE_WT_1loop_ci2 (loop_w, loop_L,2) +&
            chp_N3LO_2PE_WT_2loop_SFR(q2,2) )
       
       Vs  = chp_sigma_dot_sigma_mtx(m1,m2,m3,m4)*delta(t1,t3)*delta(t2,t4)*&
            (chp_N3LO_2PE_Vs_2loop_SFR_a(q2,2) + chp_N3LO_2PE_Vs_2loop_b(q2,2))
       
       Ws  =  chp_sigma_dot_sigma_mtx(m1,m2,m3,m4)*chp_tau_dot_tau_mtx(t1,t2,t3,t4)*&
            (chp_N3LO_2PE_Ws_1loop_ci2 (q2, loop_w, loop_L,2) +& 
            chp_N3LO_2PE_Ws_2loop_SFR(q2,2))
       
       VsigL = 0.d0
         
       !vn3lo = Vc
       vn3lo = Vc +Wc + VLS + WLS+ VT+ WT+ Vs + Ws + VsigL
       
         
         !
         ! N3LO Contact terms Eq 4.43 
         !
         sigma_dot_sigma = chp_sigma_dot_sigma_mtx(m1,m2,m3,m4)
         qxk2 = dot_product(qxk,qxk)
         term_D(1) = cn3lo(1)*q2*q2*delta(m1,m3)*delta(m2,m4)*delta(t1,t3)*delta(t2,t4) 
         term_D(2) = cn3lo(2)*kmean2*kmean2*delta(m1,m3)*delta(m2,m4)*delta(t1,t3)*delta(t2,t4) 
         term_D(3) = cn3lo(3)*q2*kmean2*delta(m1,m3)*delta(m2,m4)*delta(t1,t3)*delta(t2,t4) 
         term_D(4) = cn3lo(4)*qxk2*delta(m1,m3)*delta(m2,m4)*delta(t1,t3)*delta(t2,t4) 
         term_D(5) = cn3lo(5)*q2*q2*sigma_dot_sigma*delta(t1,t3)*delta(t2,t4)
         term_D(6) = cn3lo(6)*kmean2*kmean2*sigma_dot_sigma*delta(t1,t3)*delta(t2,t4)
         term_D(7) = cn3lo(7)*q2*kmean2*sigma_dot_sigma*delta(t1,t3)*delta(t2,t4)
         term_D(8) = cn3lo(8)*qxk2*sigma_dot_sigma*delta(t1,t3)*delta(t2,t4)
         term_D(9) = cn3lo(9)*q2*chp_spin_dot_qxk_mtx(m1,m2,m3,m4,qxk)*delta(t1,t3)*delta(t2,t4)
         term_D(10) = cn3lo(10)*kmean2*chp_spin_dot_qxk_mtx(m1,m2,m3,m4,qxk)*delta(t1,t3)*delta(t2,t4)
         term_D(11) = cn3lo(11)*q2*chp_sigma1_dot_q_sigma2_dot_q_mtx(m1,m2,m3,m4,qtrans)*delta(t1,t3)*delta(t2,t4)
         term_D(12) = cn3lo(12)*kmean2*chp_sigma1_dot_q_sigma2_dot_q_mtx(m1,m2,m3,m4,qtrans)*delta(t1,t3)*delta(t2,t4)
         term_D(13) = cn3lo(13)*q2*chp_sigma1_dot_q_sigma2_dot_q_mtx(m1,m2,m3,m4,kmean)*delta(t1,t3)*delta(t2,t4)
         term_D(14) = cn3lo(14)*kmean2*chp_sigma1_dot_q_sigma2_dot_q_mtx(m1,m2,m3,m4,kmean)*delta(t1,t3)*delta(t2,t4)
         term_D(15) = cn3lo(15)*chp_sigma1_dot_q_sigma2_dot_q_mtx(m1,m2,m3,m4,qxk)*delta(t1,t3)*delta(t2,t4)

         cont_n3lo = SUM(term_D)

    END IF
    
    ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    !
    ! DONE WITH CHIRAL ORDER CONTRIBUTIONS
    !
    ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    
    ! sum up all orders 
    ! regulator and relativity factor 
    IF (chp_chiral_order%val >= NNLO) THEN
       vdir = (vlo + cont_lo + vnlo + cont_nlo + vnnlo + cont_n3lo + vn3lo) * & 
          relativity_factor * freg(prel, pprel, chp_regcut_nnlo%val)
    ELSE
       vdir = (vlo + cont_lo + vnlo + cont_nlo + vnnlo + cont_n3lo + vn3lo) * &
            relativity_factor * freg(prel, pprel, chp_regcut%val)
    ENDIF
    
    
    !vdir = (vlo + cont_lo + cont_nlo + vnlo) * & 
    !     relativity_factor * freg(prel, pprel, chp_regcut_nnlo%val)
    
    !vdir = (vlo + cont_lo + vnlo + cont_nlo + vnnlo + vn3lo + cont_n3lo )
    
    chiral_pot =  hbarc**3 * (vdir)/volume 
    ! if(real(chiral_pot)>1) print*,chiral_pot

    
    !vdir = ( vlo + cont_lo + cont_nlo + vnlo + vnnlo )
    !* & 
    !     relativity_factor * freg(prel, pprel, chp_regcut%val)
    !chiral_pot =  hbarc**3 * (vdir) ! /(2.*pi)**3 
    
    
    !
    ! SUBTRACT TWO-BODY CENTER-OF-MASS KINETIC ENERGY
    !
    !vdir = sum(k1*k2)*delta(p,r) * delta(q,s) 
    !chiral_pot  = chiral_pot - vdir/p_mass/below_ef
    
  end function chiral_pot


  

   
  real*8 function intrinsic_kinetic(p,q,r,s)
    USE single_particle_orbits
    USE constants
    use chiral_constants
    use ang_mom_functions, only : tjs 
    
    implicit none 
    integer, intent(in) :: p,q,r,s
    INTEGER :: m1,m2,m3,m4, spin, iph, t1,t2,t3,t4, Tiso
    INTEGER :: nx1, ny1, nz1, nx2, ny2, nz2, nx3, ny3, nz3, nx4, ny4, nz4 
    REAL*8 :: k1(3), k2(3), k3(3), k4(3), kmean(3) 
    REAL*8 :: qtrans(3), prel(3), pprel(3), qxk(3)
    REAL*8 :: q2, p2, kmean2, qabs, pp2, vdir, vexc, cg1, cg2 
    REAL*8 :: delta, nucleon_mass, relativity_factor
    REAL*8 :: vlo, vnlo, vnnlo, vn3lo
    REAL*8 :: cont_lo, cont_nlo, cont_n3lo
    
    intrinsic_kinetic  = 0.d0 
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
!!$    if ( nx1 + nx2 /= nx3 + nx4 ) return 
!!$    if ( ny1 + ny2 /= ny3 + ny4 ) return 
!!$    if ( nz1 + nz2 /= nz3 + nz4 ) return 

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
    
    if ( dot_product(k1+k2-k3-k4,k1+k2-k3-k4).gt.1.d-4) return
    ! momenta in MeV
    k1 = k1*hbarc
    k2 = k2*hbarc
    
    !
    ! SUBTRACT TWO-BODY CENTER-OF-MASS KINETIC ENERGY
    !
    intrinsic_kinetic  = -sum(k1*k2)*delta(p,r) * delta(q,s) /p_mass/below_ef
    
  end function intrinsic_kinetic
  
  
  complex*16 function chiral_3nf_asym(i,j,p,k,l,q)
    USE single_particle_orbits
    USE constants
    use chiral_constants
    
    implicit none 
    INTEGER, intent(in) :: i,j,p,k,l,q
    
    chiral_3nf_asym =  ( chiral_3nf(i,j,p,k,l,q) - chiral_3nf(i,j,p,l,k,q) &
         - chiral_3nf(i,j,p,q,l,k) - chiral_3nf(i,j,p,k,q,l) &
         + chiral_3nf(i,j,p,l,q,k) + chiral_3nf(i,j,p,q,k,l) )
    
  end function chiral_3nf_asym
  
  complex*16 function chiral_3nf(p,q,r,s,t,u)
    USE single_particle_orbits
    USE constants
    use chiral_constants
    use chiral_tables 
    use ang_mom_functions, only : tjs 
    
    implicit none 
    INTEGER :: p,q,r,s,t,u, m1,m2,m3,m4,m5,m6, spin, iph, t1,t2,t3,t4,t5,t6, Tiso
    INTEGER :: nx1, ny1, nz1, nx2, ny2, nz2, nx3, ny3, nz3, nx4, ny4, nz4, nx5, ny5, nz5 , nx6, ny6, nz6 
    REAL*8 :: k1(3), k2(3), k3(3), k4(3), k5(3), k6(3), kmean(3), q1xq2(3), q2xq3(3), q1xq3(3)
    REAL*8 :: qtrans1(3), qtrans2(3), qtrans3(3), prel(3), pprel(3)
    REAL*8 :: q22, q32, q12, freg12, freg22, freg32
    REAL*8 :: delta, nucleon_mass
    complex*16 :: VNNLO, vcd, vce, v2pe, v3nf_sum, v3nf_sum1, v3nf_sum2, v3nf_sum3
    complex*16 :: sigma_dot_q252, sigma_dot_q142, sigma_dot_q251, sigma_dot_q141, sigma_dot_q_ope252, sigma_dot_q_ope141 
    complex*16 :: sigma_dot_q363, sigma_dot_q143, sigma_dot_q361, sigma_dot_q_ope363
    complex*16 :: sigma_dot_q253, sigma_dot_q362
    
    chiral_3nf=cmplx(0.d0)
    vce = 0.d0
    vcd = 0.d0
    v2pe = 0.d0

    ! 
    ! conservation of isospin 
    !
    if ( all_orbit%itzp(p) + all_orbit%itzp(q)+ all_orbit%itzp(r) /= & 
         all_orbit%itzp(s) + all_orbit%itzp(t)+ all_orbit%itzp(u) ) return 

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

    
    m1 = all_orbit%szp(p) 
    m2 = all_orbit%szp(q) 
    m3 = all_orbit%szp(r) 
    m4 = all_orbit%szp(s) 
    m5 = all_orbit%szp(t) 
    m6 = all_orbit%szp(u) 
  
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
    
    ! andreas: I dont really understand the naming convention here (Q12,Q22,Q32)
    ! But I kept it.
    !

    Q12 = sum(qtrans1*qtrans1)
    Q22 = sum(qtrans2*qtrans2)
    Q32 = sum(qtrans3*qtrans3)


    !
    !  cross product between momentum transfer and average momenta q X k 
    !
    q1xq2(1) = qtrans1(2)*qtrans2(3)-qtrans1(3)*qtrans2(2) 
    q1xq2(2) = qtrans1(3)*qtrans2(1)-qtrans1(1)*qtrans2(3) 
    q1xq2(3) = qtrans1(1)*qtrans2(2)-qtrans1(2)*qtrans2(1) 
    
    q1xq3(1) = qtrans1(2)*qtrans3(3)-qtrans1(3)*qtrans3(2) 
    q1xq3(2) = qtrans1(3)*qtrans3(1)-qtrans1(1)*qtrans3(3) 
    q1xq3(3) = qtrans1(1)*qtrans3(2)-qtrans1(2)*qtrans3(1) 
    
    q2xq3(1) = qtrans2(2)*qtrans3(3)-qtrans2(3)*qtrans3(2) 
    q2xq3(2) = qtrans2(3)*qtrans3(1)-qtrans2(1)*qtrans3(3) 
    q2xq3(3) = qtrans2(1)*qtrans3(2)-qtrans2(2)*qtrans3(1) 
    
    
    !
    ! NNLO cE contact term 
    !

    !
    ! compute all terms that are depending on freg(q12), freg(q22) 
    !
    
    ! factor 1/2 omitted because only three permutations are summed (123 and 213 are equal contributions etc.)
    ! 123 / 213 
    
!!$    sigma_dot_q252 = chp_sigma_dot_q_mtx(m2,m5,qtrans2)
!!$    sigma_dot_q142 = chp_sigma_dot_q_mtx(m1,m4,qtrans2) 
!!$    sigma_dot_q251 = chp_sigma_dot_q_mtx(m2,m5,qtrans1)
!!$    sigma_dot_q141 = chp_sigma_dot_q_mtx(m1,m4,qtrans1) 
!!$    sigma_dot_q363 = chp_sigma_dot_q_mtx(m3,m6,qtrans3)
!!$    sigma_dot_q143 = chp_sigma_dot_q_mtx(m1,m4,qtrans3) 
!!$    sigma_dot_q361 = chp_sigma_dot_q_mtx(m3,m6,qtrans1)
!!$    sigma_dot_q253 = chp_sigma_dot_q_mtx(m2,m5,qtrans3)
!!$    sigma_dot_q362 = chp_sigma_dot_q_mtx(m3,m6,qtrans2)
!!$    
!!$   
!!$    
    sigma_dot_q252 = sigma_dot_q_tab(m2,m5,nx5-nx2,ny5-ny2,nz5-nz2)
    sigma_dot_q142 = sigma_dot_q_tab(m1,m4,nx5-nx2,ny5-ny2,nz5-nz2)
    sigma_dot_q362 = sigma_dot_q_tab(m3,m6,nx5-nx2,ny5-ny2,nz5-nz2)
    sigma_dot_q251 = sigma_dot_q_tab(m2,m5,nx4-nx1,ny4-ny1,nz4-nz1) 
    sigma_dot_q141 = sigma_dot_q_tab(m1,m4,nx4-nx1,ny4-ny1,nz4-nz1) 
    sigma_dot_q361 = sigma_dot_q_tab(m3,m6,nx4-nx1,ny4-ny1,nz4-nz1) 
    sigma_dot_q363 = sigma_dot_q_tab(m3,m6,nx6-nx3,ny6-ny3,nz6-nz3)
    sigma_dot_q143 = sigma_dot_q_tab(m1,m4,nx6-nx3,ny6-ny3,nz6-nz3)
    sigma_dot_q253 = sigma_dot_q_tab(m2,m5,nx6-nx3,ny6-ny3,nz6-nz3)
    
!!$    
!!$    sigma_dot_q_ope252 = chp_sigma1_dot_q_ope(m2,m5,qtrans2)  
!!$    sigma_dot_q_ope141 = chp_sigma1_dot_q_ope(m1,m4,qtrans1)
!!$    sigma_dot_q_ope363 = chp_sigma1_dot_q_ope(m3,m6,qtrans3)  
    
    sigma_dot_q_ope252 = sigma_dot_q_ope_tab(m2,m5,nx5-nx2,ny5-ny2,nz5-nz2)
    sigma_dot_q_ope141 = sigma_dot_q_ope_tab(m1,m4,nx4-nx1,ny4-ny1,nz4-nz1) 
    sigma_dot_q_ope363 = sigma_dot_q_ope_tab(m3,m6,nx6-nx3,ny6-ny3,nz6-nz3)

    
    
    v3nf_sum1 = & 
         ! factor out [tau.tau * delta_tab(r,u)] from vce and vcd  
         tau_dot_tau_tab(t1,t2,t4,t5)*delta_tab(m3,m6)*delta_tab(t3,t6)*( &   
         ! vce terms 
         Econst * delta_tab(m1,m4)*delta_tab(m2,m5) + & 
         ! vcd terms 
         ! 123
         - Dconst * const(21) * ( & 
         sigma_dot_q252 * & 
         sigma_dot_q142/( q22 + mpi2(2) ) + & 
         ! 213
         sigma_dot_q141 * & 
         sigma_dot_q251/( q12 + mpi2(2) ) ) + & 
         sigma_dot_q_ope141 * sigma_dot_q_ope252 * ( - c1_3nf *const(2)*(4.d0*mpi2(2)/fpi2) +  & 
         ! 123 
         c3_3nf *const(2)*(2.d0/fpi2) * sum(qtrans1*qtrans2) ) ) ! + & 
    
    v3nf_sum1 = v3nf_sum1 + c4_3nf *const(2)*(1.d0/fpi2) * 0.5d0* sigma_dot_q_ope141 * sigma_dot_q_ope252 * & 
         ( & 
         !123
         chp_sigma_dot_q_mtx(m3,m6,q1xq2) * ( tau1_dot_tauXtau_tab(t3,t1,t2,t6,t4,t5) - & 
         ! 213 
         tau1_dot_tauXtau_tab(t3,t2,t1,t6,t5,t4) ) ) 
    
    
    
    ! 132/312 
    v3nf_sum2 = & 
         ! factor out [tau.tau * delta_tab(r,u)] from vce and vcd  
         tau_dot_tau_tab(t1,t3,t4,t6)*delta_tab(m2,m5)*delta_tab(t2,t5)* ( &   
         ! vce terms 
         Econst * delta_tab(m1,m4)*delta_tab(m3,m6) + & 
         ! vcd terms 
         ! 123
         - Dconst * const(21) * ( & 
         sigma_dot_q363 * & 
         sigma_dot_q143/( q32 + mpi2(2) ) + & 
         ! 213
         sigma_dot_q141 * & 
         sigma_dot_q361/( q12 + mpi2(2) ) ) + & 
         sigma_dot_q_ope141 * sigma_dot_q_ope363 * ( - c1_3nf *const(2)*(4.d0*mpi2(2)/fpi2) +  & 
         ! 123 
         c3_3nf *const(2)*(2.d0/fpi2) * sum(qtrans1*qtrans3) ) ) ! + & 

    v3nf_sum2 = v3nf_sum2 + c4_3nf *const(2)*(1.d0/fpi2) * 0.5d0* sigma_dot_q_ope141 * sigma_dot_q_ope363 * & 
         ( & 
         !123
         chp_sigma_dot_q_mtx(m2,m5,q1xq3) * ( tau1_dot_tauXtau_tab(t2,t1,t3,t5,t4,t6) - & 
         ! 213 
         tau1_dot_tauXtau_tab(t2,t3,t1,t5,t6,t4) ) ) 

        
    ! 321/231 
    v3nf_sum3 = & 
         ! factor out [tau.tau * delta_tab(r,u)] from vce and vcd  
         tau_dot_tau_tab(t2,t3,t5,t6)*delta_tab(m1,m4)*delta_tab(t1,t4)* ( &   
         ! vce terms 
         Econst * delta_tab(m2,m5)*delta_tab(m3,m6) + & 
         ! vcd terms 
         ! 123
         - Dconst * const(21) * ( & 
         sigma_dot_q363 * & 
         sigma_dot_q253/( q32 + mpi2(2) ) + & 
         ! 213
         sigma_dot_q252 * & 
         sigma_dot_q362/( q22 + mpi2(2) ) ) + & 
         sigma_dot_q_ope252 * sigma_dot_q_ope363 * ( - c1_3nf *const(2)*(4.d0*mpi2(2)/fpi2) +  & 
         ! 123 
         c3_3nf *const(2)*(2.d0/fpi2) * sum(qtrans2*qtrans3) ) ) !+ & 
    
    v3nf_sum3 = v3nf_sum3 + c4_3nf *const(2)*(1.d0/fpi2) * 0.5d0* sigma_dot_q_ope252 * sigma_dot_q_ope363 * & 
         ( & 
         !123
         chp_sigma_dot_q_mtx(m1,m4,q2xq3) * ( tau1_dot_tauXtau_tab(t1,t2,t3,t4,t5,t6) - & 
         ! 213 
         tau1_dot_tauXtau_tab(t1,t3,t2,t4,t6,t5) ) ) 
    
    
    freg12 = freg_3nflocal(q12)
    freg22 = freg_3nflocal(q22)
    freg32 = freg_3nflocal(q32)
    
    v3nf_sum = freg12*freg22 * v3nf_sum1 + & 
         freg12*freg32 * v3nf_sum2 + & 
         freg22*freg32 * v3nf_sum3 
    chiral_3nf = freg_3nfnonlocal(k1,k2,k3)*freg_3nfnonlocal(k4,k5,k6)*hbarc**6*(v3nf_sum)/volume**2
    !chiral_3nf = hbarc**6*(v3nf_sum)/volume**2
    
  end function chiral_3nf

  
  complex*16 function chiral_3nf_asym2(kk1,m1,t1,kk2,m2,t2,kk3,m3,t3,kk4,m4,t4,kk5,m5,t5,kk6,m6,t6)
    USE single_particle_orbits
    USE constants
    use chiral_constants
    
    implicit none 
    INTEGER, intent(in) :: m1,m2,m3,m4,m5,m6, t1,t2,t3,t4,t5,t6
    REAL*8, intent(in) :: kk1(3), kk2(3), kk3(3), kk4(3), kk5(3), kk6(3)
    
    chiral_3nf_asym2 =  ( & 
         + chiral_3nf2(kk1,m1,t1,kk2,m2,t2,kk3,m3,t3, kk4,m4,t4,kk5,m5,t5,kk6,m6,t6) & 
         - chiral_3nf2(kk1,m1,t1,kk2,m2,t2,kk3,m3,t3, kk5,m5,t5,kk4,m4,t4,kk6,m6,t6) & 
         - chiral_3nf2(kk1,m1,t1,kk2,m2,t2,kk3,m3,t3, kk6,m6,t6,kk5,m5,t5,kk4,m4,t4) & 
         - chiral_3nf2(kk1,m1,t1,kk2,m2,t2,kk3,m3,t3, kk4,m4,t4,kk6,m6,t6,kk5,m5,t5) & 
         + chiral_3nf2(kk1,m1,t1,kk2,m2,t2,kk3,m3,t3, kk5,m5,t5,kk6,m6,t6,kk4,m4,t4) & 
         + chiral_3nf2(kk1,m1,t1,kk2,m2,t2,kk3,m3,t3, kk6,m6,t6,kk4,m4,t4,kk5,m5,t5) )
    
  end function chiral_3nf_asym2


  complex*16 function chiral_3nf2(kk1,m1,t1,kk2,m2,t2,kk3,m3,t3,kk4,m4,t4,kk5,m5,t5,kk6,m6,t6)
    USE single_particle_orbits
    USE constants
    use chiral_constants
    use ang_mom_functions, only : tjs 
    
    implicit none 
    INTEGER, intent(in) :: m1,m2,m3,m4,m5,m6, t1,t2,t3,t4,t5,t6
    REAL*8, intent(in) :: kk1(3), kk2(3), kk3(3), kk4(3), kk5(3), kk6(3)
    INTEGER :: spin, iph, Tiso
    REAL*8 :: k1(3), k2(3), k3(3), k4(3), k5(3), k6(3), kmean(3) 
    REAL*8 :: qtrans1(3), qtrans2(3), qtrans3(3), prel(3), pprel(3)
    REAL*8 :: q22, q32, q12
    REAL*8 :: delta, nucleon_mass
    complex*16 :: VNNLO, vcd, vce, v2pe, v3nf_sum 
    complex*16 :: sigma_dot_q252, sigma_dot_q142, sigma_dot_q251, sigma_dot_q141, & 
         sigma_dot_q_ope252, sigma_dot_q_ope141 
    complex*16 :: sigma_dot_q363, sigma_dot_q143, sigma_dot_q361, sigma_dot_q_ope363
    complex*16 :: sigma_dot_q253, sigma_dot_q362


    
    chiral_3nf2=cmplx(0.d0)
    vce = 0.d0
    vcd = 0.d0
    v2pe = 0.d0
    
        
    ! momenta in MeV
    k1 = kk1*hbarc
    k2 = kk2*hbarc
    k3 = kk3*hbarc
    k4 = kk4*hbarc
    k5 = kk5*hbarc
    k6 = kk6*hbarc
    
    
    !
    ! Momentum transfer of particle 2 and 3 
    !
    qtrans1 = k4-k1
    qtrans2 = k5-k2
    qtrans3 = k6-k3
    
    ! andreas: I dont really understand the naming convention here (Q12,Q22,Q32)
    ! But I kept it.
    !

    Q12 = sum(qtrans1*qtrans1)
    Q22 = sum(qtrans2*qtrans2)
    Q32 = sum(qtrans3*qtrans3)
    
    !
    ! NNLO cE contact term 
    !
    
    !
    ! compute all terms that are depending on freg(q12), freg(q22) 
    !
    
    
    ! factor 1/2 omitted because only three permutations are summed (123 and 213 are equal contributions etc.)
    ! 123 / 213 
    
    sigma_dot_q252 = chp_sigma_dot_q_mtx(m2,m5,qtrans2)
    sigma_dot_q142 = chp_sigma_dot_q_mtx(m1,m4,qtrans2) 
    sigma_dot_q251 = chp_sigma_dot_q_mtx(m2,m5,qtrans1)
    sigma_dot_q141 = chp_sigma_dot_q_mtx(m1,m4,qtrans1) 
    sigma_dot_q_ope252 = chp_sigma1_dot_q_ope(m2,m5,qtrans2)  
    sigma_dot_q_ope141 = chp_sigma1_dot_q_ope(m1,m4,qtrans1)
    sigma_dot_q363 = chp_sigma_dot_q_mtx(m3,m6,qtrans3)
    sigma_dot_q143 = chp_sigma_dot_q_mtx(m1,m4,qtrans3) 
    sigma_dot_q361 = chp_sigma_dot_q_mtx(m3,m6,qtrans1)
    sigma_dot_q_ope363 = chp_sigma1_dot_q_ope(m3,m6,qtrans3)  
    sigma_dot_q253 = chp_sigma_dot_q_mtx(m2,m5,qtrans3)
    sigma_dot_q362 = chp_sigma_dot_q_mtx(m3,m6,qtrans2)
    
    v3nf_sum = freg_3nflocal(q12)*freg_3nflocal(q22) * ( & 
         ! factor out [tau.tau * delta(r,u)] from vce and vcd  
         chp_tau_dot_tau_mtx(t1,t2,t4,t5)*delta(m3,m6)*delta(t3,t6)* ( &   
         ! vce terms 
         Econst * delta(m1,m4)*delta(m2,m5) + & 
         ! vcd terms 
         ! 123
         - Dconst * const(21) * ( & 
         sigma_dot_q252 * & 
         sigma_dot_q142/( q22 + mpi2(2) ) + & 
         ! 213
         sigma_dot_q141 * & 
         sigma_dot_q251/( q12 + mpi2(2) ) ) + & 
         sigma_dot_q_ope141 * sigma_dot_q_ope252 * ( - c1_3nf *const(2)*(4.d0*mpi2(2)/fpi2) +  & 
         ! 123 
         c3_3nf *const(2)*(2.d0/fpi2) * dot_product(qtrans1,qtrans2) ) ) + & 
         c4_3nf *const(2)*(1.d0/fpi2) * 0.5d0* sigma_dot_q_ope141 * sigma_dot_q_ope252 * & 
         ( & 
         !123
         chp_sigma_dot_qxk_mtx(m3,m6,qtrans1,qtrans2) * chp_tau1_dot_tauXtau(t3,t1,t2,t6,t4,t5) + & 
         ! 213 
         chp_sigma_dot_qxk_mtx(m3,m6,qtrans2,qtrans1) * chp_tau1_dot_tauXtau(t3,t2,t1,t6,t5,t4) )  & 
         )

     
    
    ! 132/312 
    v3nf_sum = v3nf_sum + freg_3nflocal(q12)*freg_3nflocal(q32) * ( & 
         ! factor out [tau.tau * delta(r,u)] from vce and vcd  
         chp_tau_dot_tau_mtx(t1,t3,t4,t6)*delta(m2,m5)*delta(t2,t5)* ( &   
         ! vce terms 
         Econst * delta(m1,m4)*delta(m3,m6) + & 
         ! vcd terms 
         ! 123
         - Dconst * const(21) * ( & 
         sigma_dot_q363 * & 
         sigma_dot_q143/( q32 + mpi2(2) ) + & 
         ! 213
         sigma_dot_q141 * & 
         sigma_dot_q361/( q12 + mpi2(2) ) ) + & 
         sigma_dot_q_ope141 * sigma_dot_q_ope363 * ( - c1_3nf *const(2)*(4.d0*mpi2(2)/fpi2) +  & 
         ! 123 
         c3_3nf *const(2)*(2.d0/fpi2) * dot_product(qtrans1,qtrans3) ) ) + & 
         c4_3nf *const(2)*(1.d0/fpi2) * 0.5d0* sigma_dot_q_ope141 * sigma_dot_q_ope363 * & 
         ( & 
         !123
         chp_sigma_dot_qxk_mtx(m2,m5,qtrans1,qtrans3) * chp_tau1_dot_tauXtau(t2,t1,t3,t5,t4,t6) + & 
         ! 213 
         chp_sigma_dot_qxk_mtx(m2,m5,qtrans3,qtrans1) * chp_tau1_dot_tauXtau(t2,t3,t1,t5,t6,t4) )  & 
         )
    
        
    ! 321/231 
    v3nf_sum = v3nf_sum + freg_3nflocal(q22)*freg_3nflocal(q32) * ( & 
         ! factor out [tau.tau * delta(r,u)] from vce and vcd  
         chp_tau_dot_tau_mtx(t2,t3,t5,t6)*delta(m1,m4)*delta(t1,t4)*( &   
         ! vce terms 
         Econst * delta(m2,m5)*delta(m3,m6) + & 
         ! vcd terms 
         ! 123
         - Dconst * const(21) * ( & 
         sigma_dot_q363 * & 
         sigma_dot_q253/( q32 + mpi2(2) ) + & 
         ! 213
         sigma_dot_q252 * & 
         sigma_dot_q362/( q22 + mpi2(2) ) ) + & 
         sigma_dot_q_ope252 * sigma_dot_q_ope363 * ( - c1_3nf *const(2)*(4.d0*mpi2(2)/fpi2) +  & 
         ! 123 
         c3_3nf *const(2)*(2.d0/fpi2) * dot_product(qtrans2,qtrans3) ) ) + & 
         c4_3nf *const(2)*(1.d0/fpi2) * 0.5d0* sigma_dot_q_ope252 * sigma_dot_q_ope363 * & 
         ( & 
         !123
         chp_sigma_dot_qxk_mtx(m1,m4,qtrans2,qtrans3) * chp_tau1_dot_tauXtau(t1,t2,t3,t4,t5,t6) + & 
         ! 213 
         chp_sigma_dot_qxk_mtx(m1,m4,qtrans3,qtrans2) * chp_tau1_dot_tauXtau(t1,t3,t2,t4,t6,t5) )  & 
         )
    
    chiral_3nf2 = freg_3nfnonlocal(k1,k2,k3)*freg_3nfnonlocal(k4,k5,k6)*hbarc**6*(v3nf_sum)
    
  end function chiral_3nf2


  complex*16 function chiral_3nf_fregnonlocal(p,q,r,s,t,u)
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
    complex*16 :: VNNLO, vcd, vce, v2pe, v3nf_sum 
    complex*16 :: sigma_dot_q252, sigma_dot_q142, sigma_dot_q251, sigma_dot_q141, & 
         sigma_dot_q_ope252, sigma_dot_q_ope141 
    complex*16 :: sigma_dot_q363, sigma_dot_q143, sigma_dot_q361, sigma_dot_q_ope363
    complex*16 :: sigma_dot_q253, sigma_dot_q362
    
    chiral_3nf_fregnonlocal=cmplx(0.d0)
    vce = 0.d0
    vcd = 0.d0
    v2pe = 0.d0
    
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
    
    ! 
    ! conservation of isospin 
    !
    if ( all_orbit%itzp(p) + all_orbit%itzp(q)+ all_orbit%itzp(r) /= & 
         all_orbit%itzp(s) + all_orbit%itzp(t)+ all_orbit%itzp(u) ) return 


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

    
    m1 = all_orbit%szp(p) 
    m2 = all_orbit%szp(q) 
    m3 = all_orbit%szp(r) 
    m4 = all_orbit%szp(s) 
    m5 = all_orbit%szp(t) 
    m6 = all_orbit%szp(u) 
  
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
    
    ! andreas: I dont really understand the naming convention here (Q12,Q22,Q32)
    ! But I kept it.
    !

    Q12 = sum(qtrans1*qtrans1)
    Q22 = sum(qtrans2*qtrans2)
    Q32 = sum(qtrans3*qtrans3)
    
    !
    ! NNLO cE contact term 
    !
    
    !
    ! compute all terms that are depending on freg(q12), freg(q22) 
    !
    
    
    ! factor 1/2 omitted because only three permutations are summed (123 and 213 are equal contributions etc.)
    ! 123 / 213 
    
    !c1_3nf = 0.d0
    !c2_3nf = 0.d0
    !c3_3nf = 0.d0
    !c4_3nf = 0.d0 
    
    !Dconst = 0.d0
    !Econst = 0.d0
    
    sigma_dot_q252 = chp_sigma_dot_q_mtx(m2,m5,qtrans2)
    sigma_dot_q142 = chp_sigma_dot_q_mtx(m1,m4,qtrans2) 
    sigma_dot_q251 = chp_sigma_dot_q_mtx(m2,m5,qtrans1)
    sigma_dot_q141 = chp_sigma_dot_q_mtx(m1,m4,qtrans1) 
    sigma_dot_q_ope252 = chp_sigma1_dot_q_ope(m2,m5,qtrans2)  
    sigma_dot_q_ope141 = chp_sigma1_dot_q_ope(m1,m4,qtrans1)
    sigma_dot_q363 = chp_sigma_dot_q_mtx(m3,m6,qtrans3)
    sigma_dot_q143 = chp_sigma_dot_q_mtx(m1,m4,qtrans3) 
    sigma_dot_q361 = chp_sigma_dot_q_mtx(m3,m6,qtrans1)
    sigma_dot_q_ope363 = chp_sigma1_dot_q_ope(m3,m6,qtrans3)  
    sigma_dot_q253 = chp_sigma_dot_q_mtx(m2,m5,qtrans3)
    sigma_dot_q362 = chp_sigma_dot_q_mtx(m3,m6,qtrans2)
    
    v3nf_sum = ( & 
         ! factor out [tau.tau * delta(r,u)] from vce and vcd  
         chp_tau_dot_tau_mtx(t1,t2,t4,t5)*delta(m3,m6)* ( &   
         ! vce terms 
         Econst * delta(m1,m4)*delta(m2,m5) + & 
         ! vcd terms 
         ! 123
         - Dconst * const(21) * ( & 
         sigma_dot_q252 * & 
         sigma_dot_q142/( q22 + mpi2(2) ) + & 
         ! 213
         sigma_dot_q141 * & 
         sigma_dot_q251/( q12 + mpi2(2) ) ) + & 
         sigma_dot_q_ope141 * sigma_dot_q_ope252 * ( - c1_3nf *const(2)*(4.d0*mpi2(2)/fpi2) +  & 
         ! 123 
         c3_3nf *const(2)*(2.d0/fpi2) * dot_product(qtrans1,qtrans2) ) ) + & 
         c4_3nf *const(2)*(1.d0/fpi2) * 0.5d0* sigma_dot_q_ope141 * sigma_dot_q_ope252 * & 
         ( & 
         !123
         chp_sigma_dot_qxk_mtx(m3,m6,qtrans1,qtrans2) * chp_tau1_dot_tauXtau(t3,t1,t2,t6,t4,t5) + & 
         ! 213 
         chp_sigma_dot_qxk_mtx(m3,m6,qtrans2,qtrans1) * chp_tau1_dot_tauXtau(t3,t2,t1,t6,t5,t4) )  & 
         )

!!$    v3nf_sum = freg_3nflocal(q12)*freg_3nflocal(q22) * ( & 
!!$         ! factor out [tau.tau * delta(r,u)] from vce and vcd  
!!$         chp_tau_dot_tau_mtx(t1,t2,t4,t5)*delta(r,u)* ( &   
!!$         ! vce terms 
!!$         Econst * delta(m1,m4)*delta(m2,m5) + & 
!!$         ! vcd terms 
!!$         ! 123
!!$         - Dconst * const(21) * ( & 
!!$         chp_sigma_dot_q_mtx(m2,m5,qtrans2) * & 
!!$         chp_sigma_dot_q_mtx(m1,m4,qtrans2)/( q22 + mpi2(2) ) + & 
!!$         ! 213
!!$         chp_sigma_dot_q_mtx(m1,m4,qtrans1) * & 
!!$         chp_sigma_dot_q_mtx(m2,m5,qtrans1)/( q12 + mpi2(2) ) ) + & 
!!$         chp_sigma1_dot_q_ope(m1,m4,qtrans1) * chp_sigma1_dot_q_ope(m2,m5,qtrans2) * ( - c1_3nf *const(2)*(4.d0*mpi2(2)/fpi2) +  & 
!!$         ! 123 
!!$         c3_3nf *const(2)*(2.d0/fpi2) * dot_product(qtrans1,qtrans2) ) ) + & 
!!$         c4_3nf *const(2)*(1.d0/fpi2) * 0.5d0* chp_sigma1_dot_q_ope(m1,m4,qtrans1) * chp_sigma1_dot_q_ope(m2,m5,qtrans2) * & 
!!$         ( & 
!!$         !123
!!$         chp_sigma_dot_qxk_mtx(m3,m6,qtrans1,qtrans2) * chp_tau1_dot_tauXtau(t3,t1,t2,t6,t4,t5) + & 
!!$         ! 213 
!!$         chp_sigma_dot_qxk_mtx(m3,m6,qtrans2,qtrans1) * chp_tau1_dot_tauXtau(t3,t2,t1,t6,t5,t4) )  & 
!!$         )
!!$    

    
     
    
    ! 132/312 
    v3nf_sum = v3nf_sum + ( & 
         ! factor out [tau.tau * delta(r,u)] from vce and vcd  
         chp_tau_dot_tau_mtx(t1,t3,t4,t6)*delta(m2,m5)* ( &   
         ! vce terms 
         Econst * delta(m1,m4)*delta(m3,m6) + & 
         ! vcd terms 
         ! 123
         - Dconst * const(21) * ( & 
         sigma_dot_q363 * & 
         sigma_dot_q143/( q32 + mpi2(2) ) + & 
         ! 213
         sigma_dot_q141 * & 
         sigma_dot_q361/( q12 + mpi2(2) ) ) + & 
         sigma_dot_q_ope141 * sigma_dot_q_ope363 * ( - c1_3nf *const(2)*(4.d0*mpi2(2)/fpi2) +  & 
         ! 123 
         c3_3nf *const(2)*(2.d0/fpi2) * dot_product(qtrans1,qtrans3) ) ) + & 
         c4_3nf *const(2)*(1.d0/fpi2) * 0.5d0* sigma_dot_q_ope141 * sigma_dot_q_ope363 * & 
         ( & 
         !123
         chp_sigma_dot_qxk_mtx(m2,m5,qtrans1,qtrans3) * chp_tau1_dot_tauXtau(t2,t1,t3,t5,t4,t6) + & 
         ! 213 
         chp_sigma_dot_qxk_mtx(m2,m5,qtrans3,qtrans1) * chp_tau1_dot_tauXtau(t2,t3,t1,t5,t6,t4) )  & 
         )
    
!!$    ! 132/312 
!!$    v3nf_sum = v3nf_sum + freg_3nflocal(q12)*freg_3nflocal(q32) * ( & 
!!$         ! factor out [tau.tau * delta(r,u)] from vce and vcd  
!!$         chp_tau_dot_tau_mtx(t1,t3,t4,t6)*delta(q,t)* ( &   
!!$         ! vce terms 
!!$         Econst * delta(m1,m4)*delta(m3,m6) + & 
!!$         ! vcd terms 
!!$         ! 123
!!$         - Dconst * const(21) * ( & 
!!$         chp_sigma_dot_q_mtx(m3,m6,qtrans3) * & 
!!$         chp_sigma_dot_q_mtx(m1,m4,qtrans3)/( q32 + mpi2(2) ) + & 
!!$         ! 213
!!$         chp_sigma_dot_q_mtx(m1,m4,qtrans1) * & 
!!$         chp_sigma_dot_q_mtx(m3,m6,qtrans1)/( q12 + mpi2(2) ) ) + & 
!!$         chp_sigma1_dot_q_ope(m1,m4,qtrans1) * chp_sigma1_dot_q_ope(m3,m6,qtrans3) * ( - c1_3nf *const(2)*(4.d0*mpi2(2)/fpi2) +  & 
!!$         ! 123 
!!$         c3_3nf *const(2)*(2.d0/fpi2) * dot_product(qtrans1,qtrans3) ) ) + & 
!!$         c4_3nf *const(2)*(1.d0/fpi2) * 0.5d0* chp_sigma1_dot_q_ope(m1,m4,qtrans1) * chp_sigma1_dot_q_ope(m3,m6,qtrans3) * & 
!!$         ( & 
!!$         !123
!!$         chp_sigma_dot_qxk_mtx(m2,m5,qtrans1,qtrans3) * chp_tau1_dot_tauXtau(t2,t1,t3,t5,t4,t6) + & 
!!$         ! 213 
!!$         chp_sigma_dot_qxk_mtx(m2,m5,qtrans3,qtrans1) * chp_tau1_dot_tauXtau(t2,t3,t1,t5,t6,t4) )  & 
!!$         )

    
        
    ! 321/231 
    v3nf_sum = v3nf_sum + ( & 
         ! factor out [tau.tau * delta(r,u)] from vce and vcd  
         chp_tau_dot_tau_mtx(t2,t3,t5,t6)*delta(m1,m4)* ( &   
         ! vce terms 
         Econst * delta(m2,m5)*delta(m3,m6) + & 
         ! vcd terms 
         ! 123
         - Dconst * const(21) * ( & 
         sigma_dot_q363 * & 
         sigma_dot_q253/( q32 + mpi2(2) ) + & 
         ! 213
         sigma_dot_q252 * & 
         sigma_dot_q362/( q22 + mpi2(2) ) ) + & 
         sigma_dot_q_ope252 * sigma_dot_q_ope363 * ( - c1_3nf *const(2)*(4.d0*mpi2(2)/fpi2) +  & 
         ! 123 
         c3_3nf *const(2)*(2.d0/fpi2) * dot_product(qtrans2,qtrans3) ) ) + & 
         c4_3nf *const(2)*(1.d0/fpi2) * 0.5d0* sigma_dot_q_ope252 * sigma_dot_q_ope363 * & 
         ( & 
         !123
         chp_sigma_dot_qxk_mtx(m1,m4,qtrans2,qtrans3) * chp_tau1_dot_tauXtau(t1,t2,t3,t4,t5,t6) + & 
         ! 213 
         chp_sigma_dot_qxk_mtx(m1,m4,qtrans3,qtrans2) * chp_tau1_dot_tauXtau(t1,t3,t2,t4,t6,t5) )  & 
         )
!!$
!!$    v3nf_sum = v3nf_sum + freg_3nflocal(q22)*freg_3nflocal(q32) * ( & 
!!$         ! factor out [tau.tau * delta(r,u)] from vce and vcd  
!!$         chp_tau_dot_tau_mtx(t2,t3,t5,t6)*delta(p,s)* ( &   
!!$         ! vce terms 
!!$         Econst * delta(m2,m5)*delta(m3,m6) + & 
!!$         ! vcd terms 
!!$         ! 123
!!$         - Dconst * const(21) * ( & 
!!$         chp_sigma_dot_q_mtx(m3,m6,qtrans3) * & 
!!$         chp_sigma_dot_q_mtx(m2,m5,qtrans3)/( q32 + mpi2(2) ) + & 
!!$         ! 213
!!$         chp_sigma_dot_q_mtx(m2,m5,qtrans2) * & 
!!$         chp_sigma_dot_q_mtx(m3,m6,qtrans2)/( q22 + mpi2(2) ) ) + & 
!!$         chp_sigma1_dot_q_ope(m2,m5,qtrans2) * chp_sigma1_dot_q_ope(m3,m6,qtrans3) * ( - c1_3nf *const(2)*(4.d0*mpi2(2)/fpi2) +  & 
!!$         ! 123 
!!$         c3_3nf *const(2)*(2.d0/fpi2) * dot_product(qtrans2,qtrans3) ) ) + & 
!!$         c4_3nf *const(2)*(1.d0/fpi2) * 0.5d0* chp_sigma1_dot_q_ope(m2,m5,qtrans2) * chp_sigma1_dot_q_ope(m3,m6,qtrans3) * & 
!!$         ( & 
!!$         !123
!!$         chp_sigma_dot_qxk_mtx(m1,m4,qtrans2,qtrans3) * chp_tau1_dot_tauXtau(t1,t2,t3,t4,t5,t6) + & 
!!$         ! 213 
!!$         chp_sigma_dot_qxk_mtx(m1,m4,qtrans3,qtrans2) * chp_tau1_dot_tauXtau(t1,t3,t2,t4,t6,t5) )  & 
!!$         )
!!$    
!!$  
    
    chiral_3nf_fregnonlocal = freg_3nfnonlocal(k1,k2,k3)*freg_3nfnonlocal(k4,k5,k6) * & 
         hbarc**6*(v3nf_sum)/volume**2
    
  end function chiral_3nf_fregnonlocal

  complex*16 function chiral_3nf_contact(p,q,r,s,t,u)
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
    complex*16 :: VNNLO, vcd, vce, v2pe, v3nf_sum 
    complex*16 :: sigma_dot_q252, sigma_dot_q142, sigma_dot_q251, sigma_dot_q141, & 
         sigma_dot_q_ope252, sigma_dot_q_ope141 
    complex*16 :: sigma_dot_q363, sigma_dot_q143, sigma_dot_q361, sigma_dot_q_ope363
    complex*16 :: sigma_dot_q253, sigma_dot_q362
    
    chiral_3nf_contact=cmplx(0.d0)
    vce = 0.d0
    vcd = 0.d0
    v2pe = 0.d0
    
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
    
    ! 
    ! conservation of isospin 
    !
    if ( all_orbit%itzp(p) + all_orbit%itzp(q)+ all_orbit%itzp(r) /= & 
         all_orbit%itzp(s) + all_orbit%itzp(t)+ all_orbit%itzp(u) ) return 


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

    
    m1 = all_orbit%szp(p) 
    m2 = all_orbit%szp(q) 
    m3 = all_orbit%szp(r) 
    m4 = all_orbit%szp(s) 
    m5 = all_orbit%szp(t) 
    m6 = all_orbit%szp(u) 
  
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
    
    ! andreas: I dont really understand the naming convention here (Q12,Q22,Q32)
    ! But I kept it.
    !

    Q12 = sum(qtrans1*qtrans1)
    Q22 = sum(qtrans2*qtrans2)
    Q32 = sum(qtrans3*qtrans3)
    
    !
    ! NNLO cE contact term 
    !
    
    !
    ! compute all terms that are depending on freg(q12), freg(q22) 
    !
    
    
    ! factor 1/2 omitted because only three permutations are summed (123 and 213 are equal contributions etc.)
    ! 123 / 213 
    
    !Dconst = 0.d0
    !Econst = 0.d0
    
    v3nf_sum = 0.d0 
    v3nf_sum = ( & 
         ! factor out [tau.tau * delta(r,u)] from vce and vcd  
         delta(m3,m6) * delta(t3,t6) * ( &   
         ! vce terms 
         Econst * delta(m1,m4)*delta(m2,m5) ) )
    
    
    ! 132/312 
    v3nf_sum = v3nf_sum +  ( & 
         ! factor out [tau.tau * delta(r,u)] from vce and vcd  
         delta(m2,m5) * delta(t2,t5)* ( &   
         ! vce terms 
         Econst * delta(m1,m4)*delta(m3,m6) ) )
    
    ! 321/231 
    v3nf_sum = v3nf_sum +  ( & 
         ! factor out [tau.tau * delta(r,u)] from vce and vcd  
         delta(m1,m4) * delta(t1,t4) * ( &   
         ! vce terms 
         Econst * delta(m2,m5)*delta(m3,m6) ) )       
    
    !chiral_3nf_contact =  freg_3nfnonlocal(k1,k2,k3) *freg_3nfnonlocal(k4,k5,k6) * & 
    !     hbarc**6*(v3nf_sum)/volume**2
    chiral_3nf_contact =  hbarc**6*(v3nf_sum)/volume**2
    
  end function chiral_3nf_contact

  complex*16 function chiral_3nf_old(p,q,r,s,t,u)
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
    complex*16 :: VNNLO, vcd, vce, v2pe
    
    chiral_3nf_old=cmplx(0.d0)
    vce = 0.d0
    vcd = 0.d0
    v2pe = 0.d0
    
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
    
    ! 
    ! conservation of isospin 
    !
    if ( all_orbit%itzp(p) + all_orbit%itzp(q)+ all_orbit%itzp(r) /= & 
         all_orbit%itzp(s) + all_orbit%itzp(t)+ all_orbit%itzp(u) ) return 


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

    
    m1 = all_orbit%szp(p) 
    m2 = all_orbit%szp(q) 
    m3 = all_orbit%szp(r) 
    m4 = all_orbit%szp(s) 
    m5 = all_orbit%szp(t) 
    m6 = all_orbit%szp(u) 
  
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
    
    ! andreas: I dont really understand the naming convention here (Q12,Q22,Q32)
    ! But I kept it.
    !

    Q12 = sum(qtrans1*qtrans1)
    Q22 = sum(qtrans2*qtrans2)
    Q32 = sum(qtrans3*qtrans3)
    
    !
    ! NNLO cE contact term 
    !
    
    ! andreas: defined Econst
    !
    ! factor 1/2 omitted because only three permutations are summed (123 and 213 are equal contributions etc.)
    vce = Econst * ( chp_tau_dot_tau_mtx(t1,t2,t4,t5)*delta(m1,m4)*delta(m2,m5)*delta(m3,m6)*delta(t3,t6) * & 
         freg_3nflocal(q12)*freg_3nflocal(q22) + & 
         chp_tau_dot_tau_mtx(t1,t3,t4,t6)*delta(m1,m4)*delta(m3,m6)*delta(m2,m5)*delta(t2,t5) * & 
         freg_3nflocal(q12)*freg_3nflocal(q32) + & 
         chp_tau_dot_tau_mtx(t2,t3,t5,t6)*delta(m2,m5)*delta(m3,m6)*delta(m1,m4)*delta(t1,t4) * & 
         freg_3nflocal(q22)*freg_3nflocal(q32) ) 
    
    !
    ! NNLO cD term 
    !
    vcd = - Dconst * const(21) * &
         !123
         ( freg_3nflocal(q12)*freg_3nflocal(q22) * & 
         chp_tau_dot_tau_mtx(t1,t2,t4,t5) * chp_sigma_dot_q_mtx(m2,m5,qtrans2) * & 
         chp_sigma_dot_q_mtx(m1,m4,qtrans2) * delta(m3,m6)*delta(t3,t6)/( q22 + mpi2(2) ) + & 
         !213
         freg_3nflocal(q12)*freg_3nflocal(q22) * & 
         chp_tau_dot_tau_mtx(t2,t1,t5,t4) * chp_sigma_dot_q_mtx(m1,m4,qtrans1) * & 
         chp_sigma_dot_q_mtx(m2,m5,qtrans1) * delta(m3,m6)*delta(t3,t6)/( q12 + mpi2(2) ) + & 
         !132
         freg_3nflocal(q12)*freg_3nflocal(q32) * & 
         chp_tau_dot_tau_mtx(t1,t3,t4,t6) * chp_sigma_dot_q_mtx(m3,m6,qtrans3) * & 
         chp_sigma_dot_q_mtx(m1,m4,qtrans3) * delta(m2,m5)*delta(t2,t5)/( q32 + mpi2(2) ) + & 
         !312
         freg_3nflocal(q12)*freg_3nflocal(q32) * & 
         chp_tau_dot_tau_mtx(t3,t1,t6,t4) * chp_sigma_dot_q_mtx(m3,m6,qtrans1) * & 
         chp_sigma_dot_q_mtx(m1,m4,qtrans1) * delta(m2,m5)*delta(t2,t5)/( q12 + mpi2(2) ) + & 
         !231
         freg_3nflocal(q22)*freg_3nflocal(q32) * & 
         chp_tau_dot_tau_mtx(t2,t3,t5,t6) * chp_sigma_dot_q_mtx(m3,m6,qtrans3) * & 
         chp_sigma_dot_q_mtx(m2,m5,qtrans3) * delta(m1,m4)*delta(t1,t4)/( q32 + mpi2(2) ) + & 
         !321
         freg_3nflocal(q22)*freg_3nflocal(q32) * & 
         chp_tau_dot_tau_mtx(t3,t2,t6,t5) * chp_sigma_dot_q_mtx(m3,m6,qtrans2) * & 
         chp_sigma_dot_q_mtx(m2,m5,qtrans2) * delta(m1,m4)*delta(t1,t4)/( q22 + mpi2(2) ) )

    
    !
    ! NNLO two-pion exchange
    !
    
    v2pe = - c1_3nf *const(2)*(4.d0*mpi2(2)/fpi2) * ( & !only three permutations, omitted factor 1/2
         ! 123 
         chp_sigma1_dot_q_ope(m1,m4,qtrans1) * chp_sigma1_dot_q_ope(m2,m5,qtrans2) *  &
         chp_tau_dot_tau_mtx(t1,t2,t4,t5)*delta(m3,m6)*delta(t3,t6) * freg_3nflocal(q22)*freg_3nflocal(q12) +  & 
         ! 231 
         chp_sigma1_dot_q_ope(m2,m5,qtrans2) * chp_sigma1_dot_q_ope(m3,m6,qtrans3) *  &
         chp_tau_dot_tau_mtx(t2,t3,t5,t6)*delta(m1,m4)*delta(t1,t4) * freg_3nflocal(q22)*freg_3nflocal(q32) +  & 
         ! 312
         chp_sigma1_dot_q_ope(m3,m6,qtrans3) * chp_sigma1_dot_q_ope(m1,m4,qtrans1) * &
         chp_tau_dot_tau_mtx(t1,t3,t4,t6)*delta(m2,m5)*delta(t2,t5) * freg_3nflocal(q12)*freg_3nflocal(q32))

    v2pe = v2pe + c3_3nf *const(2)*(2.d0/fpi2) * ( & !only three permutations, omitted factor 1/2
         ! 123 
         chp_sigma1_dot_q_ope(m1,m4,qtrans1) * chp_sigma1_dot_q_ope(m2,m5,qtrans2) *  &
         chp_tau_dot_tau_mtx(t1,t2,t4,t5)* dot_product(qtrans1,qtrans2) * & 
         delta(m3,m6)*delta(t3,t6) * freg_3nflocal(q22)*freg_3nflocal(q12) +  & 
         ! 231 
         chp_sigma1_dot_q_ope(m2,m5,qtrans2) * chp_sigma1_dot_q_ope(m3,m6,qtrans3) *  &
         chp_tau_dot_tau_mtx(t2,t3,t5,t6)* dot_product(qtrans2,qtrans3) * & 
         delta(m1,m4)*delta(t1,t4) * freg_3nflocal(q22)*freg_3nflocal(q32) +  & 
         ! 312
         chp_sigma1_dot_q_ope(m3,m6,qtrans3) * chp_sigma1_dot_q_ope(m1,m4,qtrans1) * &
         chp_tau_dot_tau_mtx(t1,t3,t4,t6)* dot_product(qtrans1,qtrans3) * & 
         delta(m2,m5)*delta(t2,t5) * freg_3nflocal(q12)*freg_3nflocal(q32))
    
    v2pe = v2pe + c4_3nf *const(2)*(1.d0/fpi2) * 0.5d0*( &  !the factor 1/2 is there because all 6 permutations follow
         !123
         chp_sigma1_dot_q_ope(m1,m4,qtrans1) * chp_sigma1_dot_q_ope(m2,m5,qtrans2) * & 
         chp_sigma_dot_qxk_mtx(m3,m6,qtrans1,qtrans2) * chp_tau1_dot_tauXtau(t3,t1,t2,t6,t4,t5) * & 
         freg_3nflocal(q12)*freg_3nflocal(q22) + & 
         !231
         chp_sigma1_dot_q_ope(m2,m5,qtrans2) * chp_sigma1_dot_q_ope(m3,m6,qtrans3) * & 
         chp_sigma_dot_qxk_mtx(m1,m4,qtrans2,qtrans3) * chp_tau1_dot_tauXtau(t1,t2,t3,t4,t5,t6) * & 
         freg_3nflocal(q22)*freg_3nflocal(q32) + & 
         !312
         chp_sigma1_dot_q_ope(m3,m6,qtrans3) * chp_sigma1_dot_q_ope(m1,m4,qtrans1) * & 
         chp_sigma_dot_qxk_mtx(m2,m5,qtrans3,qtrans1) * chp_tau1_dot_tauXtau(t2,t3,t1,t5,t6,t4) * & 
         freg_3nflocal(q32)*freg_3nflocal(q12) + & 
         ! 213 
         chp_sigma1_dot_q_ope(m2,m5,qtrans2) * chp_sigma1_dot_q_ope(m1,m4,qtrans1) * & 
         chp_sigma_dot_qxk_mtx(m3,m6,qtrans2,qtrans1) * chp_tau1_dot_tauXtau(t3,t2,t1,t6,t5,t4) * & 
         freg_3nflocal(q22)*freg_3nflocal(q12) + &  
         ! 132
         chp_sigma1_dot_q_ope(m1,m4,qtrans1) * chp_sigma1_dot_q_ope(m3,m6,qtrans3) * & 
         chp_sigma_dot_qxk_mtx(m2,m5,qtrans1,qtrans3) * chp_tau1_dot_tauXtau(t2,t1,t3,t5,t4,t6) * & 
         freg_3nflocal(q12)*freg_3nflocal(q32) + & 
         ! 321
         chp_sigma1_dot_q_ope(m3,m6,qtrans3) * chp_sigma1_dot_q_ope(m2,m5,qtrans2) * & 
         chp_sigma_dot_qxk_mtx(m1,m4,qtrans3,qtrans2) * chp_tau1_dot_tauXtau(t1,t3,t2,t4,t6,t5) * &
         freg_3nflocal(q32)*freg_3nflocal(q22) )
    
    chiral_3nf_old = hbarc**6*(vce+vcd+v2pe)/volume**2
    
  end function chiral_3nf_old
  
  

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
    REAL*8 :: v0r, v0t, v0s, kr,kt,ks, vr,vt,vs,vs_ex,vcentral, vsigma, spin_exc1, spin_exc2
  
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
    !if ( all_orbit%szp(p) + all_orbit%szp(q) /= all_orbit%szp(r) + all_orbit%szp(s) ) return 
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
    v0t=178.0  ! MeV
    v0s=91.85  ! MeV
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
    
    !vr = 0.5*vr 
    !vs = 4.*vs 
    
    vcentral=0.25d0*(vr+vs)*delta(m1,m3)*delta(m2,m4)*delta(t1,t3)*delta(t2,t4) 
    vsigma=-0.25d0*(vr+vs)*chp_sigma_dot_sigma(m1,m2,m3,m4)*delta(t1,t3)*delta(t2,t4)
    
    
    
    !vcentral(r_12)+vsigma(r_12)*sigma_1.sigma_2.
    !vdir = vcentral + vsigma
    !vdir = 10.*vcentral + 4.*vsigma
    
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
    
    
    !vdir = vs*delta(m1,m3)*delta(m2,m4)*delta(t1,t3)*delta(t2,t4) !(vs+vr+vt)
    vdir = vs+vr+vt
    !vmom_minnesota 
    
    !
    ! subtract two-body part of CoM kinetic energy
    !
    
    vmom_minnesota = vdir/volume ! - sum(k1*k2)*delta(p,r) * delta(q,s) *hbarc**2/p_mass/below_ef
    

  end function vmom_minnesota


  real*8 function vmom_minnesota2(k1,m1,t1,k2,m2,t2,k3,m3,t3,k4,m4,t4)
    USE single_particle_orbits
    USE constants
    use chiral_constants
  
    implicit none 
    INTEGER, intent(in) :: m1,m2,m3,m4, t1,t2,t3,t4
    INTEGER :: nx1, ny1, nz1, nx2, ny2, nz2, nx3, ny3, nz3, nx4, ny4, nz4 
    REAL*8, intent(in) :: k1(3), k2(3), k3(3), k4(3)
    REAL*8 :: qtrans(3), prel(3), pprel(3)
    REAL*8 :: q2, p2, qabs, pp2, vdir, vexc, cg1, cg2 
    REAL*8 :: delta, nucleon_mass, relativity_factor
    REAL*8 :: v0r, v0t, v0s, kr,kt,ks, vr,vt,vs,vs_ex,vcentral, vsigma, spin_exc1, spin_exc2
    
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
    v0t=178.0  ! MeV
    v0s=91.85  ! MeV
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
    
    vdir = vs+vr+vt
        
    !
    ! subtract two-body part of CoM kinetic energy
    !
    
    vmom_minnesota2 = vdir ! - sum(k1*k2)*delta(p,r) * delta(q,s) *hbarc**2/p_mass/below_ef
    

  end function vmom_minnesota2

  
  complex*16 function chiral_pot2(kk1,m1,t1,kk2,m2,t2,kk3,m3,t3,kk4,m4,t4)
    USE single_particle_orbits
    USE constants
    use chiral_constants
    use ang_mom_functions, only : tjs 
    
    implicit none 
    INTEGER, intent(in) :: m1,m2,m3,m4, t1,t2,t3,t4
    REAL*8, intent(in) :: kk1(3), kk2(3), kk3(3), kk4(3)
    INTEGER :: spin, iph, Tiso
    INTEGER :: nx1, ny1, nz1, nx2, ny2, nz2, nx3, ny3, nz3, nx4, ny4, nz4 

!!$    REAL*8 :: kmean(3),  k1(3), k2(3), k3(3), k4(3)
!!$    REAL*8 :: qtrans(3), prel(3), pprel(3), qxk(3)
!!$    REAL*8 :: q2, p2, kmean2, qabs, pp2
!!$    complex*16 :: vdir, vexc, cg1, cg2 
!!$    REAL*8 :: delta, nucleon_mass, relativity_factor
!!$    

    REAL*8 :: k1(3), k2(3), k3(3), k4(3), kmean(3) 
    REAL*8 :: qtrans(3), prel(3), pprel(3), qxk(3)
    REAL*8 :: q2, p2, kmean2, qabs, pp2, cg1, cg2 
    REAL*8 :: delta, nucleon_mass, relativity_factor
    
    COMPLEX*16 :: vlo, vnlo, vnnlo, vn3lo, vdir
    COMPLEX*16 :: cont_lo, cont_nlo, cont_n3lo
    ! LIST OF OPERATOR STRUCTURES
    COMPLEX*16 :: Vc
    COMPLEX*16 :: Wc
    COMPLEX*16 :: Vs
    COMPLEX*16 :: Ws
    COMPLEX*16 :: VLS
    COMPLEX*16 :: WLS
    COMPLEX*16 :: VsigL
    COMPLEX*16 :: WsigL
    COMPLEX*16 :: VT
    COMPLEX*16 :: WT
    COMPLEX*16 :: Vsigk
    COMPLEX*16 :: Wsigk
    COMPLEX*16 :: term_CS
    COMPLEX*16 :: term_CT
    COMPLEX*16 :: term_C(1:7)
    REAL*8 :: loop_w
    REAL*8 :: loop_wtilde2
    REAL*8 :: loop_s
    REAL*8 :: loop_L
    REAL*8 :: loop_A
    
        
    
    chiral_pot2 = 0.d0 
    
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
  
    term_C = 0.0D0
    
    vlo = 0.0D0 ; vnlo = 0.0D0 ; vnnlo = 0.0D0 ; vn3lo = 0.0D0
    
    cont_lo = 0.0D0 ; cont_nlo = 0.0D0 ; cont_n3lo = 0.0D0
    
    ! momenta in MeV
    k1 = kk1*hbarc
    k2 = kk2*hbarc
    k3 = kk3*hbarc
    k4 = kk4*hbarc
        
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
    
    ! LO CIB one-pion exchanges
    ! [1] Eq. 4.77-4.79
    
    !
    ! Note: the LO contribution is modified at NLO
    ! to include CIB effects. Therefore, should you go beyond
    ! LO, make sure to not include the vlo contribution 
    ! in the sum defining vdir.
    !
    
    !vlo = chp_sigma1_dot_q_sigma2_dot_q(m1,m2,m3,m4,qtrans)*chp_tau_dot_tau(t1,t2,t3,t4)*chp_one_pion_exchange(q2, 2)
    
    WT = 0.d0 
    IF (t1 + t2 == 0) THEN 
       
       WT = 0.d0 
       DO Tiso = 0, 2, 2
          cg1 = tjs(1,1,Tiso,t1,t2,-(t1+t2))*iph( (t1+t2)/2 )*sqrt(Tiso+1.d0)
          cg2 = tjs(1,1,Tiso,t3,t4,-(t3+t4))*iph( (t3+t4)/2 )*sqrt(Tiso+1.d0)
          
          WT = WT + cg1*cg2*(-chp_one_pion_exchange(q2, 0) + iph(Tiso/2+1)*2.0D0*chp_one_pion_exchange(q2,1)) 
       end DO
       
       
    END IF
    IF (t1 + t2 /= 0) THEN
       WT = chp_one_pion_exchange(q2, 0) 
    END IF
    !! should rel_corr_1PE be here ?
    vlo = WT*chp_sigma1_dot_q_sigma2_dot_q_mtx(m1,m2,m3,m4,qtrans)

    
    !
    ! leading order CIB contacts 
    !
    term_CS = CS((t1+t2)/2)*delta(m1,m3)*delta(m2,m4)*delta(t1,t3)*delta(t2,t4)
    term_CT = CT((t1+t2)/2) * chp_sigma_dot_sigma_mtx(m1,m2,m3,m4) *delta(t1,t3)*delta(t2,t4) 
    cont_lo = term_CS + term_CT
    !cont_lo = 4.*cont_lo
    
    Wc = chp_tau_dot_tau_mtx(t1,t2,t3,t4)*delta(m1,m3)*delta(m2,m4)* & 
         chp_NLO_two_pion_exchange_Wc(q2, loop_L, loop_w, 2)
    
    Vs = chp_NLO_two_pion_exchange_Vs(q2, loop_L, 2)* &
         chp_sigma_dot_sigma_mtx(m1,m2,m3,m4)*delta(t1,t3)*delta(t2,t4) 
    
    VT = chp_NLO_two_pion_exchange_VT(loop_L, 2)* & 
         chp_sigma1_dot_q_sigma2_dot_q_mtx(m1,m2,m3,m4,qtrans)*delta(t1,t3)*delta(t2,t4) 
    
    vnlo = ( WT + WC + VS + VT ) 
    
    !
    ! next-to-leading order contacts 
    !
    term_C(1) = cnlo(1)*q2*delta(m1,m3)*delta(m2,m4)*delta(t1,t3)*delta(t2,t4) 
    
    term_C(2) = cnlo(2)*kmean2*delta(m1,m3)*delta(m2,m4)*delta(t1,t3)*delta(t2,t4)
    
    term_C(3) = cnlo(3)*q2*chp_sigma_dot_sigma_mtx(m1,m2,m3,m4)*delta(t1,t3)*delta(t2,t4)
    
    term_C(4) = cnlo(4)*kmean2*chp_sigma_dot_sigma_mtx(m1,m2,m3,m4)*delta(t1,t3)*delta(t2,t4)
    
    term_C(5) = cnlo(5)*chp_spin_dot_qxk_mtx(m1,m2,m3,m4,qxk)*delta(t1,t3)*delta(t2,t4)
    
    term_C(6) = cnlo(6)*chp_sigma1_dot_q_sigma2_dot_q_mtx(m1,m2,m3,m4,qtrans)*delta(t1,t3)*delta(t2,t4)
    
    term_C(7) = cnlo(7)*chp_sigma1_dot_q_sigma2_dot_q_mtx(m1,m2,m3,m4,kmean)*delta(t1,t3)*delta(t2,t4) 
    
    cont_nlo = SUM(term_C)
       
       
       
    ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    !
    ! NEXT TO NEXT TO LEADING ORDER 
    ! 
    ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    
    chp_regcut%val = 3.0D0
    
    ! NNLO loop functions
    
    loop_wtilde2 = chp_NNLO_two_pion_exchange_loop_wtilde2(q2,2)
    loop_A = chp_NNLO_sfr_two_pion_exchange_loop_A(qabs, q2, 2)
    
    ! chp_itope type is set in the header of this function
    Vc  = chp_NNLO_two_pion_exchange_Vc (q2, loop_w, loop_A, loop_wtilde2, 2)*&
         delta(m1,m3)*delta(m2,m4)*delta(t1,t3)*delta(t2,t4)
    
    Wc  = chp_NNLO_two_pion_exchange_Wc (q2, loop_w, loop_A, loop_wtilde2, 2)*&
         chp_tau_dot_tau_mtx(t1,t2,t3,t4)*delta(m1,m3)*delta(m2,m4)
    
    VLS = chp_NNLO_two_pion_exchange_VLS(loop_A, loop_wtilde2, 2)*&
         chp_spin_dot_qxk_mtx(m1,m2,m3,m4,qxk)*delta(t1,t3)*delta(t2,t4)
    
    WLS = chp_NNLO_two_pion_exchange_WLS(loop_w, loop_A, 2)*&
         chp_spin_dot_qxk_mtx(m1,m2,m3,m4,qxk)*chp_tau_dot_tau_mtx(t1,t2,t3,t4)
    
    VT  = chp_NNLO_two_pion_exchange_VT (loop_w, loop_A, loop_wtilde2, 2)*&
         chp_sigma1_dot_q_sigma2_dot_q_mtx(m1,m2,m3,m4,qtrans)*delta(t1,t3)*delta(t2,t4) 
    
    WT  = chp_NNLO_two_pion_exchange_WT (q2, loop_w, loop_A, 2)*&
         chp_sigma1_dot_q_sigma2_dot_q_mtx(m1,m2,m3,m4,qtrans)*chp_tau_dot_tau_mtx(t1,t2,t3,t4)
    
    Vs  = chp_NNLO_two_pion_exchange_Vs (q2, loop_w, loop_A, loop_wtilde2, 2)*&
         chp_sigma_dot_sigma_mtx(m1,m2,m3,m4)*delta(t1,t3)*delta(t2,t4) 
    
    Ws  = chp_NNLO_two_pion_exchange_Ws (q2, loop_w, loop_A, 2)*&
         chp_sigma_dot_sigma_mtx(m1,m2,m3,m4)*chp_tau_dot_tau_mtx(t1,t2,t3,t4)
    
    vnnlo = (Vc + Wc + VLS + WLS + VT + WT + Vs + Ws)
    
    
    ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    !
    ! DONE WITH CHIRAL ORDER CONTRIBUTIONS
    !
    ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    
    ! sum up all orders 
    ! regulator and relativity factor 
    
    vdir = (vlo + cont_lo + vnlo + cont_nlo + vnnlo + cont_n3lo + vn3lo) * &
         relativity_factor * freg(prel, pprel, 3.d0 )

    chiral_pot2 =  hbarc**3 * (vdir) 
    
  end function chiral_pot2


  
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
  
  !vmom_yukawa = vmom_yukawa - (vdir - vexc)*hbarc**2/p_mass/below_ef
  vmom_yukawa = vmom_yukawa - (vdir)*hbarc**2/p_mass/below_ef

end function vmom_yukawa



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
  REAL*8 :: delta, s2, t2, va, vb, mu_a, mu_b
  INTEGER :: lang, spin, ispin
  
  spin = 0
  ispin = 0 
  lang = 0 

  v0r=200.0  ! MeV
  v0t=178.0  ! MeV
  v0s=91.85  ! MeV
  kr=1.487  ! fm**-2
  kt=0.639  ! fm**-2
  ks=0.465  ! fm**-2
  vr=v0r*exp(-kr*rr**2)
  vt=-v0t*exp(-kt*rr**2)
  vs=-v0s*exp(-ks*rr**2)
  
  !vr = vr 
  !vs = 3.*vs 
  
  vcentral=0.25d0*(vr+vs)
  vsigma=-vcentral
  
  ispin = 0 
  if ( mod(lang+spin, 2) == 0 ) ispin = 1
  !write(6,*) lang, spin, ispin,  mod(lang+spin, 2)
  
  if ( spin == 1 ) s2 = -3.d0
  if ( spin == 0 ) s2 =  1.d0
  if ( ispin == 1 ) t2 = -3.d0
  if ( ispin == 0 ) t2 =  1.d0
  
  vr = vr * (3.d0 - t2 - s2 - t2*s2 )/8.d0
  vs = vs * (3.d0 + t2 - 3.d0*s2 - t2*s2 )/16.d0
  vt = vt * (3.d0 - 3.d0*t2 + s2 - t2*s2 )/16.d0

  vrspace_minnesota = vr + vs + vt !vcentral + vsigma

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

