!             Coupled Cluster Nuclear Matter Project
!
!             Author:   Thomas Papenbrock and Gaute Hagen 
!             ADDRESS:  Physics Division, Oak Ridge National Laboratory
!             E-MAIL:   hageng@ornl.gov, tpapenbr@utk.edu
!             LANGUAGE: F90/F95, MPI
!             LAST UPGRADE : May 2012
PROGRAM ccm_kspace 
  use parallel
  USE ang_mom_functions
  use chiral_constants
  use single_particle_orbits
  USE constants

  implicit none
  real*8  :: factor, startwtime , endwtime, x, diff, e00, einf
  complex*16 :: e0_av, mbpt2_av, eccsd_av, eccsdt_av
  integer :: nx, ny, nz, i , nxx, nyy, nzz 
  integer :: spec_points(10,10,10)
  
  call mpi_init(ierror)
  
  call mpi_comm_size(mpi_comm_world,num_procs,ierror)
  call mpi_comm_rank(mpi_comm_world,iam,ierror)
  
  
  open(5,file='ccm_in', access="SEQUENTIAL")
  master=0
  
  
  !call diagonalize_h_osc_basis
  !stop
  startwtime = MPI_WTIME()
  
  read(5,*);read(5,*)  cE, cD 
  CALL commons_to_angmom
  call init_chp_constants 
  
  IF (IAM == 0 ) WRITE(6,*) 
  IF (IAM == 0 ) WRITE(6,"(A12)") '3NF LECs:'
  IF (IAM == 0 ) WRITE(6,"(A12,F30.16)") 'C1', LEC_c1_3NF 
  IF (IAM == 0 ) WRITE(6,"(A12,F30.16)") 'C3', LEC_c3_3NF 
  IF (IAM == 0 ) WRITE(6,"(A12,F30.16)") 'C4', LEC_c4_3NF
  IF (IAM == 0 ) WRITE(6,"(A12,F30.16)") 'cE', cE
  IF (IAM == 0 ) WRITE(6,"(A12,F30.16)") 'cD', cD
  IF (IAM == 0 ) WRITE(6,*) 
  call allocate_sp_data 
  ! 
  ! setup special points from 10^3 twist angles in [0:pi] and minimizing: diff = abs(E0_inf + E0)/N 
  ! This is for nuclear matter and N2LO_opt (with CIB)
  spec_points = 0 
  if ( kf == 1.d0  ) spec_points(1,4,4) = 1 
  if ( kf == 1.1d0 ) spec_points(1,2,4) = 1 
  if ( kf == 1.2d0 ) spec_points(2,2,4) = 1 
  if ( kf == 1.3d0 ) spec_points(1,3,4) = 1 
  if ( kf == 1.4d0 ) spec_points(4,4,7) = 1 
  if ( kf == 1.5d0 ) spec_points(2,3,5) = 1 
  if ( kf == 1.6d0 ) spec_points(1,5,10) = 1 
  if ( kf == 1.7d0 ) spec_points(4,6,8) = 1 
  if ( kf == 1.8d0 ) spec_points(4,4,8) = 1 
  if ( kf == 1.9d0 ) spec_points(2,7,8) = 1 
  if ( kf == 2.d0  ) spec_points(3,7,10) = 1 
  

  
  !call plot_3nf  !compute_hf_inf
  !goto 10 
  
  select case( boundary_conditions )
  case('TABC') 
     e0_av =0.d0 
     mbpt2_av = 0.d0 
     eccsd_av = 0.d0 
     eccsdt_av = 0.d0 
     
     diff = 100.d0 
     do nx = 1, size(twist_angle, 2)
        do ny = nx, size(twist_angle, 2)
           do nz = ny, size(twist_angle, 2)
              
              !if ( spec_points(nx,ny,nz) == 0 ) cycle 
              
              !twist_angle = 0.d0 
              CALL setup_sp_data(nx,ny,nz)
              !call setup_channel_structures
              !call setup_ph_channel_structures
              call normal_ordered_hamiltonian
              !call setup_t_amplitudes
              !CALL ccd_iter 
              if ( nx == ny .and. nx == nz ) factor = 1.d0
              if ( nx == ny .and. nx /= nz ) factor = 3.d0 
              if ( nx == nz .and. nx /= ny ) factor = 3.d0
              if ( ny == nz .and. ny /= nx ) factor = 3.d0 
              if ( nx /= ny .and. nx /= nz .and. ny /= nz ) factor = 6.d0 
           
              e0_av = e0_av + factor * wxx(nx)*wxx(ny)*wxx(nz)*e0/(1.d0*pi)**3
              mbpt2_av = mbpt2_av + factor * wxx(nx)*wxx(ny)*wxx(nz)*mbpt2/(1.d0*pi)**3
              eccsd_av = eccsd_av + factor * wxx(nx)*wxx(ny)*wxx(nz)*eccsd/(1.d0*pi)**3

              e0_av = e0 
!!$           mbpt2_av = mbpt2
!!$           eccsd_av = eccsd
              einf = -2.16599648
           
              ! for nn only 
              !einf = -5.996597
              !einf = -7.9968459
              !einf = -10.100443
              !einf = -12.182785
              !einf = -14.087980
              !einf = -15.630038
              !einf = -16.596422
              !einf = -16.754645
              !einf = -15.862494
              !einf = -13.682065
              !if ( abs(e0/below_ef - einf ) < diff ) then 
              !   diff = abs(e0/below_ef - einf ) 
              !   nxx = nx 
              !   nyy = ny 
              !   nzz = nz 
              !   e00 = e0/below_ef 
              !end if
              
              call deallocate_structures
           
           
           end do
        end do
     end do
     
     !  if (iam == 0 ) write(6,*) diff, nxx, nyy, nzz 
     !  if ( iam == 0 ) write(6,*) e00
     
     !10 continue  
  
  
  case( 'PBC' )
     
     twist_angle = 0.d0 
     CALL setup_sp_data(1,1,1)
     call precalc_chp_functions
     
     !call compute_v3nf_memory
     call setup_channel_structures
     call setup_ph_channel_structures
     call normal_ordered_hamiltonian
     
     !call mbpt2_3nf
     !goto 10
     
     call setup_t_amplitudes
     CALL ccd_iter 
     !call deallocate_structures
     
     e0_av = e0
     eccsd_av = eccsd
     mbpt2_av = mbpt2 
     eccsdt_av = eccsdt
  end select
  

!10 continue  
  open(1137,file='mydata.dat',status='old',action='write',form='formatted',position="append")
  if ( iam == 0 ) write(6,*)
  if ( iam == 0 ) write(6,*)
  if ( iam == 0 ) write(6,*) '---------------------------------------------------------------------'
  if ( iam == 0 ) write(6,*) 'Boxsize', lx, 'Nmax', nmax
  if ( iam == 0 ) write(6,*) 'A, Density, Kf', below_ef, dens, kf   
  if ( iam == 0 ) write(6,*) 'Twist averaged E0/N', e0_av/below_ef 
  if ( iam == 0 ) write(6,*) 'Twist averaged mbpt2', mbpt2_av/below_ef 
  if ( iam == 0 ) write(6,*) 'Twist averaged CCSD', eccsd_av/below_ef
  if ( iam == 0 ) write(6,*) 'kinetic_TF', 0.3d0* (hbarc**2*kf**2/p_mass)*below_ef  
  if ( iam == 0 ) write(6,*) dens, e0_av/below_ef, eccsd_av
  if ( iam == 0 ) write(6,*)
  if ( iam == 0 ) write(6,*) 'E/N', real(kf), real(dens), real(e0_av/below_ef), & 
       real(mbpt2_av/below_ef), real(eccsd_av/below_ef), real(eccsdt_av/below_ef)

  if ( iam == 0 ) write(1137,'(5f10.3)') cE, cD, real(kf), real(dens), real(eccsd_av/below_ef)

  !'(a,1x,6(f10.5,2x))'
  !write(6,*) below_ef, abs(1.d0  + (e0_av/below_ef)/16.671932d0)
  !write(6,*) below_ef, abs(1.d0  + (e0_av/below_ef)/21.479137d0 ) 
  !write(6,*) below_ef, abs(1.d0  + (e0_av/below_ef)/17.68343315d0)

  if ( iam == 0 ) write(6,*) '---------------------------------------------------------------------'
  endwtime = MPI_WTIME()
  if ( iam == 0 ) write(6,*) 'Total execution time for CCD code', endwtime - startwtime

  call mpi_finalize(ierror)

END PROGRAM ccm_kspace


!
! HF energy in the Thermodynamic Limit. 
!
SUBROUTINE plot_3nf
  USE PARALLEL
  USE single_particle_orbits
  USE constants
  use one_body_operators
  use wave_functions
  use chiral_potentials
    
  IMPLICIT NONE
  real*8, allocatable :: klab(:), theta(:), phi(:), wklab(:), wtheta(:), wphi(:)
  real*8, allocatable :: kx(:), wkx(:), ky(:), wky(:), kz(:), wkz(:)
  integer :: kk1,kk2,kk3,kk4,kk5, i, j
  integer :: nk, nphi, ntheta,ii,jj, kk, i1,j1,k1,i2,j2,k2, i3,j3,k3,tz1, tz2,tz3, sz1, sz2,sz3, mm ,nn, nconfs 
  real*8 :: p1(3), p2(3), p3(3), p4(3), p5(3), p6(3), p(3), q(3)
  complex*16 :: vint, sum1, hf_ene, hf_sum, tmp_sum, hf_3nf 
  
  nk = 50
  allocate( wklab(nk), klab(nk) )
  allocate( wtheta(100), theta(100) )
  
  call gauss_legendre(kf,10.d0,klab,wklab,nk)
  call gauss_legendre(0.d0,pi,theta, wtheta, 100)
  
  tz1 =   1
  tz2 =  -1
  tz3 =  -1
  
  sz1 =  -1
  sz2 =   1
  sz3 =  -1
  
  do i = 1, nk 
     do j = 1, 100
        
        p1(1) = 0.
        p1(2) = 0.
        p1(3) = 0.
     
        p2(1) = 0.
        p2(2) = 0.
        p2(3) = kf
     
        p3(1) = 0.
        p3(2) = kf
        p3(3) = 0.
     
        p4(1) = 0. 
        p4(2) = 0. 
        p4(3) = klab(i) 
        
        p5(1) = klab(i)*sin(theta(j))
        p5(2) = 0.d0 
        p5(3) = klab(i)*cos(theta(j))
        
        p6 = p1 + p2 + p3 - p4 -p5
        vint = ( chiral_3nf_asym2(p1,sz1,tz1,p2,sz2,tz2,p3,sz3,tz3, p4,sz1,tz1,p5,sz2,tz2,p6,sz3,tz3)) 
        
        p = 0.5d0*(p4-p5)
        q = (2./3.)*(p6-0.5d0*(p4+p5))
        
        write(61,'(4(g20.10,1x))') sqrt(dot_product(p,p)), SQRT(dot_product(q,q)), dble(vint), aimag(vint)
     end do
  end do
  
  
end SUBROUTINE plot_3nf

! HF energy in the Thermodynamic Limit. 
!
SUBROUTINE compute_hf_inf
  USE PARALLEL
  USE single_particle_orbits
  USE constants
  use one_body_operators
  use wave_functions
  use chiral_potentials
    
  IMPLICIT NONE
  real*8, allocatable :: klab(:), theta(:), phi(:), wklab(:), wtheta(:), wphi(:)
  real*8, allocatable :: kx(:), wkx(:), ky(:), wky(:), kz(:), wkz(:)
  integer :: kk1,kk2,kk3,kk4,kk5, i
  integer :: nk, nphi, ntheta,ii,jj, kk, i1,j1,k1,i2,j2,k2, i3,j3,k3,tz1, tz2,tz3, sz1, sz2,sz3, mm ,nn, nconfs 
  real*8 :: p1(3), p2(3), p3(3)
  complex*16 :: vint, sum1, hf_ene, hf_sum, tmp_sum, hf_3nf 
  integer, allocatable :: spherical_confs(:) 
  integer ::  nconfs_tot, number_mtxel_iam, diff, n
  integer, allocatable :: nconf_low(:), nconf_high(:)
  
  nk = 6
  nphi = 18
  ntheta = 18
  
  allocate( wklab(nk), klab(nk) )
  allocate( theta(ntheta), wtheta(ntheta) )
  allocate( phi(nphi), wphi(nphi) )
  
  call gauss_legendre(0.d0,pi,theta,wtheta,ntheta)
  call gauss_legendre(0.d0,kf,klab,wklab,nk)
  call gauss_legendre(0.d0,2.d0*pi,phi,wphi,nphi)
  
  
  
  
  nconfs = 0 
  do tz1 = -1, 1, 2 
     do sz1 = -1, 1, 2 
        do i1 = 1, nk
           do j1 = 1, ntheta
              do k1 = 1, nphi
              
                 nconfs= nconfs + 1 
              end do
           end do
        end do
     end do
  end do
  allocate( spherical_confs(5*nconfs) )
  
  nconfs = 0 
  spherical_confs = 0
  do tz1 = -1, 1, 2
     do sz1 = -1, 1, 2 
        do i1 = 1, nk
           do j1 = 1, ntheta
              do k1 = 1, nphi
              
                 nconfs= nconfs + 1 
                 kk1 = nconfs*5 
                 kk2 = nconfs*5 - 1
                 kk3 = nconfs*5 - 2
                 kk4 = nconfs*5 - 3
                 kk5 = nconfs*5 - 4
                 spherical_confs(kk1) = tz1 
                 spherical_confs(kk2) = sz1
                 spherical_confs(kk3) = i1
                 spherical_confs(kk4) = j1
                 spherical_confs(kk5) = k1
                 
              end do
           end do
        end do
     end do
  end do

  allocate( nconf_low(num_procs), nconf_high(num_procs) )
  nconf_low = 0
  nconf_high = 0

  !
  ! find number of configs on each process
  !
  nconfs_tot = nconfs
  number_mtxel_iam = floor(real(nconfs_tot/num_procs ))
  !if ( iam == 0 ) write(6,*) number_mtxel_iam 
  diff = nconfs - num_procs * number_mtxel_iam 
  
  !if ( iam == 0 ) write(6,*) diff
  !if ( iam == 0 ) write(6,*) 'Total number of mtx-elements', bra_confs
  !if ( iam == 0 ) write(6,*) 'Total number of mtx-elements', number_mtxel_iam*num_procs + diff
  
  !if ( iam == num_procs - 1 ) nconfs = number_mtxel_iam + diff
  !if ( iam /= num_procs - 1 ) nconfs = number_mtxel_iam  
  
  n = 1 
  do i = 1, num_procs 
     if ( i < num_procs ) then
        nconf_low(i) = n 
        nconf_high(i) = n + number_mtxel_iam - 1
     else
        nconf_low(i) = n 
        nconf_high(i) = n + number_mtxel_iam - 1 + diff
     end if
     n = n + number_mtxel_iam 
  end do
  
  do i = 1, num_procs
     if ( iam == 0 ) write(6,*) 'iam',i,'lower/upper', nconf_low(i), nconf_high(i)
  end do
  
  if ( iam == 0 ) write(6,*) nconfs 
  
  
  hf_sum = 0.d0 
  do ii = nconf_low(iam+1), nconf_high(iam+1) !1, nconfs 
     
     kk1 = ii*5 
     kk2 = ii*5 - 1
     kk3 = ii*5 - 2
     kk4 = ii*5 - 3
     kk5 = ii*5 - 4
     tz1 = spherical_confs(kk1)
     sz1 = spherical_confs(kk2)
     i1 = spherical_confs(kk3) 
     j1 = spherical_confs(kk4) 
     k1 = spherical_confs(kk5) 
     
     p1(1) = klab(i1)*sin(theta(j1))*cos(phi(k1))  
     p1(2) = klab(i1)*sin(theta(j1))*sin(phi(k1))  
     p1(3) = klab(i1)*cos(theta(j1))
     
     tmp_sum = 0.d0 
     !$omp parallel default(shared) private(jj,kk1,kk2,kk3,kk4,kk5,tz2,sz2,i2,j2,k2,p2, vint) 
     !$omp do schedule(dynamic), reduction(+:tmp_sum)
     do jj = ii+1, nconfs 
        
        kk1 = jj*5 
        kk2 = jj*5 - 1
        kk3 = jj*5 - 2
        kk4 = jj*5 - 3
        kk5 = jj*5 - 4
        tz2 = spherical_confs(kk1)
        sz2 = spherical_confs(kk2)
        i2 = spherical_confs(kk3) 
        j2 = spherical_confs(kk4) 
        k2 = spherical_confs(kk5) 
        
        p2(1) = klab(i2)*sin(theta(j2))*cos(phi(k2))  
        p2(2) = klab(i2)*sin(theta(j2))*sin(phi(k2))  
        p2(3) = klab(i2)*cos(theta(j2))
                          
        !vint = ( vmom_minnesota2(p1,sz1,tz1,p2,sz2,tz2, p1,sz1,tz1,p2,sz2,tz2) -  &
        !     vmom_minnesota2(p1,sz1,tz1,p2,sz2,tz2, p2,sz2,tz2,p1,sz1,tz1) )
        !write(6,*) p1, p2 
        vint = ( chiral_pot2(p1,sz1,tz1,p2,sz2,tz2, p1,sz1,tz1,p2,sz2,tz2) -  &
             chiral_pot2(p1,sz1,tz1,p2,sz2,tz2, p2,sz2,tz2,p1,sz1,tz1) )
        !write(6,*) vint 
        tmp_sum  = tmp_sum + wklab(i1)*wklab(i2)*wtheta(j1)*wtheta(j2)*wphi(k1)*wphi(k2) * & 
             klab(i1)**2*klab(i2)**2*sin(theta(j1))*sin(theta(j2)) * vint 
        
     end do
     !$omp end do
     !$omp end parallel

     hf_sum = hf_sum + tmp_sum 

  end do
  
  call mpi_reduce(hf_sum,hf_ene,1,mpi_complex16,mpi_sum,master, &
       mpi_comm_world,ierror)
  
  hf_ene = hf_ene/(2.*pi)**6/dens
  if ( iam == 0 ) write(6,*) 'neutron matter', dens, real(hf_ene), real(hf_ene + 0.3d0* (2.d0*hbarc**2*kf**2/p_mass))
  if ( iam == 0 ) write(6,*) 'sym. nucl. matter', dens, real(hf_ene), real(hf_ene + 0.3d0* (hbarc**2*kf**2/p_mass))

  return 
  
  hf_sum = 0.d0 
  !
  ! 3nf hf contribution 
  !
  do ii = nconf_low(iam+1), nconf_high(iam+1) !1, nconfs 
     
     kk1 = ii*5 
     kk2 = ii*5 - 1
     kk3 = ii*5 - 2
     kk4 = ii*5 - 3
     kk5 = ii*5 - 4
     tz1 = spherical_confs(kk1)
     sz1 = spherical_confs(kk2)
     i1 = spherical_confs(kk3) 
     j1 = spherical_confs(kk4) 
     k1 = spherical_confs(kk5) 
     
     p1(1) = klab(i1)*sin(theta(j1))*cos(phi(k1))  
     p1(2) = klab(i1)*sin(theta(j1))*sin(phi(k1))  
     p1(3) = klab(i1)*cos(theta(j1))
     
     tmp_sum = 0.d0 
     !$omp parallel default(shared) private(jj,kk1,kk2,kk3,kk4,kk5,tz2,sz2,i2,j2,k2,p2, & 
     !$omp kk, tz3,sz3,i3,j3,k3,p3,vint )
     !$omp do schedule(dynamic), reduction(+:tmp_sum)
     do jj = ii+1, nconfs 
        
        kk1 = jj*5 
        kk2 = jj*5 - 1
        kk3 = jj*5 - 2
        kk4 = jj*5 - 3
        kk5 = jj*5 - 4
        tz2 = spherical_confs(kk1)
        sz2 = spherical_confs(kk2)
        i2 = spherical_confs(kk3) 
        j2 = spherical_confs(kk4) 
        k2 = spherical_confs(kk5) 
        
        p2(1) = klab(i2)*sin(theta(j2))*cos(phi(k2))  
        p2(2) = klab(i2)*sin(theta(j2))*sin(phi(k2))  
        p2(3) = klab(i2)*cos(theta(j2))
        
        do kk = jj+1, nconfs 
        
           kk1 = kk*5 
           kk2 = kk*5 - 1
           kk3 = kk*5 - 2
           kk4 = kk*5 - 3
           kk5 = kk*5 - 4
           tz3 = spherical_confs(kk1)
           sz3 = spherical_confs(kk2)
           i3 = spherical_confs(kk3) 
           j3 = spherical_confs(kk4) 
           k3 = spherical_confs(kk5) 
        
           p3(1) = klab(i3)*sin(theta(j3))*cos(phi(k3))  
           p3(2) = klab(i3)*sin(theta(j3))*sin(phi(k3))  
           p3(3) = klab(i3)*cos(theta(j3))
           
           vint = ( chiral_3nf_asym2(p1,sz1,tz1,p2,sz2,tz2,p3,sz3,tz3, p1,sz1,tz1,p2,sz2,tz2,p3,sz3,tz3)) 
           
           tmp_sum  = tmp_sum + wklab(i1)*wklab(i2)*wklab(i3)*wtheta(j1)*wtheta(j2)*wtheta(j3)*wphi(k1)*wphi(k2)*wphi(k3) * & 
                klab(i1)**2*klab(i2)**2*klab(i3)**2*sin(theta(j1))*sin(theta(j2))*sin(theta(j3)) * vint 
        
        end do
     end do
     !$omp end do
     !$omp end parallel
     hf_sum = hf_sum + tmp_sum 
  end do
  
  call mpi_reduce(hf_sum,hf_3nf,1,mpi_complex16,mpi_sum,master, &
       mpi_comm_world,ierror)
  
  hf_ene = hf_ene + hf_3nf/(2.d0*pi)**9/dens
  
  write(6,*) 'iam', iam, 'done'
  if ( iam == 0 ) write(6,*) 'neutron matter', dens, real(hf_ene), real(hf_ene + 0.3d0* (2.d0*hbarc**2*kf**2/p_mass))
  if ( iam == 0 ) write(6,*) 'sym. nucl. matter', dens, real(hf_ene), real(hf_ene + 0.3d0* (hbarc**2*kf**2/p_mass))
  
  
  
end SUBROUTINE compute_hf_inf



subroutine deallocate_structures
  USE single_particle_orbits
  USE constants
  use one_body_operators
  use t2_storage
  use configurations
  USE KSPACE 
  use CHIRAL_POTENTIALS

  IMPLICIT NONE
  INTEGER :: i, channel, channel1, channel2, channel3, channel4, channel6, channel5, local_ch_id , nx,ny,nz,tz 
  
  do channel   = 1, channels%number_pppp_confs 
     nx = channels%pppp_quantum_numbers(channel*4)
     ny = channels%pppp_quantum_numbers(channel*4-1)
     nz = channels%pppp_quantum_numbers(channel*4-2)
     tz = channels%pppp_quantum_numbers(channel*4-3)

     !
     ! check that pppp channel is restricted to hhpp channels. 
     !
     channel2 = locate_channel(3,tz, nx, ny, nz)
     if ( channel2 == 0 ) cycle 
     if ( check_my_channel(channel2) == 0 ) cycle
     local_ch_id = channel2
     
     deallocate( vnn_pppp(channel2)%val ) 
  end do
  deallocate( vnn_pppp )
  
  do i = 1, size( lookup_hhhh_configs, 2 ) 
     DEALLOCATE( lookup_hhhh_configs(1,i)%ival ) 
  end do
  do i = 1, size( lookup_hhhp_configs, 2 )
     DEALLOCATE( lookup_hhhp_configs(1,i)%ival ) 
     DEALLOCATE( lookup_hhhp_configs(2,i)%ival )
  end do
  do i = 1, size( lookup_hhpp_configs, 2 )
     DEALLOCATE( lookup_hhpp_configs(1,i)%ival ) 
     DEALLOCATE( lookup_hhpp_configs(2,i)%ival )
  end do
  do i = 1, size( lookup_pppp_configs, 2 ) 
     DEALLOCATE( lookup_pppp_configs(1,i)%ival ) 
  end do
  do i = 1, size( lookup_hphp_configs, 2 ) 
     DEALLOCATE( lookup_hphp_configs(1,i)%ival ) 
  end do
  do i = 1, size( lookup_ph_phhp_configs, 2) 
     DEALLOCATE( lookup_ph_phhp_configs(1,i)%ival )
     DEALLOCATE( lookup_ph_phhp_configs(2,i)%ival )
  end do
  do i = 1, size( lookup_ph_hphp_configs, 2) 
     DEALLOCATE( lookup_ph_hphp_configs(1,i)%ival )
  end do
  do i = 1, size( lookup_ph_hpph_configs, 2) 
     DEALLOCATE( lookup_ph_hpph_configs(1,i)%ival )
     DEALLOCATE( lookup_ph_hpph_configs(2,i)%ival )
  end do
  
  DEALLOCATE( lookup_ph_hpph_configs )
  DEALLOCATE( lookup_ph_hphp_configs )
  DEALLOCATE( lookup_ph_phhp_configs )
  deallocate(  lookup_hhhh_configs ) 
  deallocate(  lookup_hhhp_configs ) 
  deallocate(  lookup_hhpp_configs ) 
  deallocate(  lookup_hphp_configs ) 
  deallocate(  lookup_pppp_configs ) 
  deallocate( locate_t2channel ) 
  deallocate( locate_channel ) 
  deallocate( locate_ph_channel )
  
  deallocate( ph_ph_phhp%ival )
  deallocate( hp_ph_phhp%ival ) 
  deallocate( hp_ph_hphp%ival )
  deallocate( hp_ph_hpph%ival ) 
  deallocate( ph_ph_hpph%ival ) 
  
  deallocate( hh_hhhh%ival ) 
  deallocate( hh_hhhp%ival ) 
  deallocate( hh_hhpp%ival ) 
  deallocate( hp_hphp%ival ) 
!  deallocate( hp_hppp%ival ) 
  deallocate( hp_hhhp%ival ) 
  deallocate( pp_hhpp%ival ) 
!  deallocate( pp_hppp%ival ) 
  deallocate( pp_pppp%ival ) 
  DEALLOCATE( channels%hhhh_quantum_numbers ) 
  DEALLOCATE( channels%hhhp_quantum_numbers ) 
  DEALLOCATE( channels%hhpp_quantum_numbers ) 
  DEALLOCATE( channels%hphp_quantum_numbers ) 
  DEALLOCATE( channels%pppp_quantum_numbers ) 
  
  DEALLOCATE( channels%ph_phhp_quantum_numbers )
  DEALLOCATE( channels%ph_hphp_quantum_numbers )
  DEALLOCATE( channels%ph_hpph_quantum_numbers )
  
!!$  do i = 1, size( vnn_ph_hpph )
!!$     deallocate( vnn_ph_hpph(i)%val )
!!$  end do
  do i = 1, size( vnn_hhhh ) 
     deallocate( vnn_hhhh(i)%val )
  end do
  do i = 1, size( vnn_hhhp ) 
     deallocate( vnn_hhhp(i)%val )
  end do
  do i = 1, size( vnn_hhpp ) 
     deallocate( vnn_hhpp(i)%val )
  end do
  do i = 1, size( vnn_hphp ) 
     deallocate( vnn_hphp(i)%val )
  end do
  deallocate( vnn_hhhh ) 
  deallocate( vnn_hhhp ) 
  deallocate( vnn_hhpp ) 
  deallocate( vnn_hphp ) 
!!$  deallocate( vnn_ph_hpph )
  
  do i = 1, size( t2_ccm )
     deallocate( t2_ccm(i)%val )
     deallocate( t2_ccm_eqn(i)%val )
  end do
  
  deallocate( t2_ccm, t2_ccm_eqn )
  
  do i = 1, size( t2_ph_ccm) 
     deallocate( t2_ph_ccm(i)%val )
     deallocate( t2_ph_ccm_eqn(i)%val )
  end do
  deallocate( t2_ph_ccm, t2_ph_ccm_eqn )

  deallocate( mapping ) 
  deallocate( check_my_channel ) 
  deallocate( my_channel_low, my_channel_high )
  deallocate( fully_stored_channel )
  deallocate( mapping_hhpp ) 
  deallocate( check_my_channel_hhpp ) 
  deallocate( my_hhpp_channel_low, my_hhpp_channel_high )
  deallocate( mapping_hphp ) 
  deallocate( check_my_channel_hphp ) 
  deallocate( mapping_ph_hphp )
  deallocate( check_my_channel_ph_hphp )
  deallocate( my_ph_hphp_channel_low, my_ph_hphp_channel_high) 
  
  DEALLOCATE( orbitof )


  deallocate( tkin, fock_mtx )
  
end subroutine deallocate_structures

!
! read in the single-particle data
! and set up the single-particle orbitals; 
! allocate sp space.
!
SUBROUTINE allocate_sp_data
  USE PARALLEL
  USE single_particle_orbits
  USE constants
  use one_body_operators
  use wave_functions
  use chiral_potentials
    
  IMPLICIT NONE
  INTEGER :: i,j,k, ii, jj, a, nn, b, c,d, ndim, i1,j1, n
  INTEGER :: nx,ny,nz,sz,tzp, ntwist,  nkx, nky, nkz, p, q, r, s
  REAL*8 ::  sum1, vmom_yukawa_rel, input, vrspace_yukawa, vrspace_minnesota, kx,ky,kz,rr
  CHARACTER(len=100) :: calc_type, input_type 
  
  !     Read number of particles 
  read(5,*);read(5,*) below_ef
  !read(5,*);read(5,*) Lx, Ly, Lz, nmax
  read(5,*);read(5,*) calc_type, input_type 
  read(5,*);read(5,*) boundary_conditions 
  read(5,*);read(5,*) input, ntwist, nmax
  read(5,*);read(5,*) cc_approx
  read(5,*);read(5,*) tnf_switch, tnf_approx
  
  if ( iam == 0 ) write(6,'(2(a,1x))') 'Boundary Conditions', boundary_conditions 
  if ( iam == 0 ) write(6,'(2(a,1x))') 'CC approx', cc_approx
  if ( iam == 0 ) write(6,*) 'TNF:', tnf_switch, 'TNF-approx:', tnf_approx

  szmin = -1 
  szmax =  1
  
  allocate( twist_angle(3, ntwist))
  twist_angle = 0.d0 
  allocate( xx(ntwist), wxx(ntwist) ) 
  xx = 0.d0; wxx = 0.d0 
  call gauss_legendre(0.d0,pi,xx,wxx,ntwist)
  
  
  sum1 = 0. 
  do i = 1, ntwist
     twist_angle(:,i) = xx(i) 
     !twist_angle(:,i) = dble(i-1.d0)*pi/(ntwist-1)
     !write(6,*) twist_angle(:,i) !, xx(i), wxx(i)
  end do
  !stop
  
  select case( calc_type ) 
  case('pnm') 
     ! for neutron matter 
     tzmin = 1 
     tzmax = 1
     if ( input_type == 'density' ) THEN
        dens = input 
        kf = (3.d0*pi**2*dens)**(1.d0/3.d0)
     ELSEIF ( input_type == 'kfermi' ) THEN 
        kf = input 
        dens = kf**3/3.d0/pi**2
     end if
  case('snm') 
     tzmin = -1 
     tzmax =  1
     if ( input_type == 'density' ) THEN
        ! for nuclear matter 
        dens = input 
        kf = (3.d0*pi**2*dens/2.d0)**(1.d0/3.d0)
     ELSEIF ( input_type == 'kfermi' ) then 
        kf = input 
        dens = 2.d0*kf**3/3.d0/pi**2
     end if
  end select
  
  lx = (real(below_ef)/dens)**(1.d0/3.d0)
  ly = lx
  lz = lx
  
  if ( iam == 0 ) write(6,*) 'Boxsize', lx 
  if ( iam == 0 ) write(6,*) 'Density', dens
  if ( iam == 0 ) write(6,*) 'K_f', kf
  if ( iam == 0 ) write(6,*) 'Nmax', nmax
  
  if ( ly == 0 .and. lz == 0 ) then
     volume = lx
  elseif ( ly /= 0 .and. lz == 0 ) then
     volume = lx*ly
  else
     volume = lx*ly*lz
  end if
  
  !
  ! For now only neutrons 
  ! 
  
  allocate( kz_mesh(nmax*2+1) ) 
  allocate( ky_mesh(nmax*2+1) ) 
  allocate( kx_mesh(nmax*2+1) ) 
  kz_mesh = 0.d0
  ky_mesh = 0.d0
  kx_mesh = 0.d0
  
  
  nkx = 2*nmax+1
  nky = nkx
  nkz = nkx
  all_orbit%total_orbits = (szmax-szmin+2) * (tzmax-tzmin+2) * (nkx*nky*nkz)/4 
  
  !all_orbit%total_orbits = (nkx*nky*nkz)
  !     Setup all possible orbit information
  CALL allocate_sp_array(all_orbit,all_orbit%total_orbits)
  
  !
  !
  !
  tot_orbs = all_orbit%total_orbits
  if ( iam == 0 ) write(6,*) 'Total number of states', tot_orbs
  
  allocate( erg( nkx*nky*nkz) ) 
  allocate( indx_inv(nkx*nky*nkz) , indx(nkx*nky*nkz) ) 
  erg = 0.D0
  indx_inv = 0
  indx = 0 
  
  
end SUBROUTINE allocate_sp_data


SUBROUTINE setup_sp_data(itwistx, itwisty, itwistz)
  USE PARALLEL
  USE single_particle_orbits
  USE constants
  use one_body_operators
  use wave_functions
  use chiral_potentials
  
  
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: itwistx, itwisty, itwistz  
  INTEGER :: i,j,k, ii, jj, a, nn, b,kk, c,d, ndim, i1,j1, n
  INTEGER :: nx,ny,nz,sz,tzp, nkx, nky, nkz, p, q, r, s, is, it
  REAL*8 ::  sum1, vmom_yukawa_rel, vrspace_yukawa, vrspace_minnesota, kx,ky,kz,rr
  REAL*8 :: theta_x, theta_y, theta_z, sum_sz, sum_tz 
  
  kz_mesh = 0.d0
  ky_mesh = 0.d0
  kx_mesh = 0.d0
  
  
  !
  ! setup k-mesh
  !
  nkz = 1
  nkx = 1 
  nky = 1
    
  theta_x = twist_angle(1,itwistx)
  theta_y = twist_angle(2,itwisty)
  theta_z = twist_angle(3,itwistz)
  
  
  kx_mesh(1) = theta_x/lx 
  ii = 1
  do i = 1, nmax
     do nx = -i,i,2*i
        ii = ii + 1
        kx_mesh(ii) = (2.d0*nx*pi+theta_x)/lx
        
        
     end do
  end do

  do i = 1, size(kx_mesh)
     if ( iam == 0 ) write(6,*) 'K-mesh', kx_mesh(i)
  end do
  nkx = ii 
  
  if ( abs(lz) /= 0.d0 ) then 
     kz_mesh(1) = theta_z/lz  
     ii = 1
     do i = 1, nmax
        do nx = -i,i,2*i
           ii = ii + 1
           kz_mesh(ii) = (2.d0*nx*pi+theta_z )/lz
        end do
     end do
     nkz = ii 
  end if

  if ( abs(ly) /= 0.d0 ) then 
     ky_mesh(1) = theta_y/ly 
     ii = 1
     do i = 1, nmax
        do nx = -i,i,2*i
           ii = ii + 1
           ky_mesh(ii) = (2.d0*nx*pi+theta_y)/ly
           
           !write(6,*) 'k-mesh', ky_mesh(ii) 
        end do
     end do
     nky = ii 
  end if
  
  
  !write(6,*) (szmax-szmin+2)/2
  !stop
  ! total number of states 
  ! XXXXXXXXXXXXXXXXXXX 
  
  erg = 0.d0 
  ii = 0
  do i = 1, nkx
     do j = 1, nky
        do k = 1, nkz
           ii = ii + 1 
           erg(ii) = (kx_mesh(i)**2+ky_mesh(j)**2+kz_mesh(k)**2)
           
        end do
     end do
  end do
  
  call index( ii, erg, indx_inv) 
    
  do i=1, size(indx) 
     ! XXXXXXXXXXXXXXXXXXXXXXXXXXXX
     !indx_inv(i) = i 
     ! XXXXXXXXXXXXXXXXXXXXXXXXXXXX
     
     indx(indx_inv(i))=i
  end do
  
  
  ii = 0
  do i = 1, nkx
     do j = 1, nky
        do k = 1, nkz
           ii = ii + 1 
           jj = indx(ii) 
           
           is=0
           do sz = szmin, szmax, (szmax-szmin+2)/2
              is=is+1
              it=0
              do tzp = tzmin, tzmax, (tzmax-tzmin+2)/2
                 it=it+1
                 kk= ((szmax-szmin+2) * (tzmax-tzmin+2))/4 * (jj-1) +2*(it-1)+is
                 
                 
                 ! periodic boundary condition
                 all_orbit%kx(kk) = kx_mesh(i) 
                 all_orbit%ky(kk) = ky_mesh(j)
                 all_orbit%kz(kk) = kz_mesh(k)
                 
                 all_orbit%nx(kk) = nint( (kx_mesh(i) * lx - theta_x ) /2./pi) 
                 all_orbit%ny(kk) = nint( (ky_mesh(j) * ly - theta_y ) /2./pi)
                 all_orbit%nz(kk) = nint( (kz_mesh(k) * lz - theta_z ) /2./pi)
                 
                 all_orbit%szp(kk) = sz
                 all_orbit%itzp(kk) = tzp
                 all_orbit%e(kk)= (all_orbit%kx(kk)**2 + all_orbit%ky(kk)**2 + all_orbit%kz(kk)**2) * hbarc**2/p_mass/2.
                 
            
                 
              end do
           end do
        end do
     end do
  end do
  DO i=1, all_orbit%total_orbits
     if ( i <= below_ef ) THEN 
        all_orbit%orbit_status(i) = 'hole' 
     else 
        all_orbit%orbit_status(i) = 'particle' 
     end if
  end DO
  
  nkxmax = maxval( all_orbit%nx ) ; nkymax = maxval( all_orbit%ny ); nkzmax = maxval( all_orbit%nz )
  
  allocate( orbitof(-nkxmax:nkxmax, -nkymax:nkymax, -nkzmax:nkzmax, szmin:szmax, tzmin:tzmax) )
  orbitof = 0 
  
  if ( iam == 0 ) write(6,*) 'nkmin, nkmax', nkxmax, nkymax, nkzmax 
  
  if(iam == 0)write(6,*)'Occupied j-orbits' 
  sum_tz = 0 
  sum_sz = 0
  DO i=1,below_ef 
     if(iam == 0)WRITE(6,'((i6,1x),1x,1(g20.10,1x),1x,5(I3,1x),a)')  &
          i, all_orbit%e(i), (all_orbit%nx(i)), (all_orbit%ny(i)), (all_orbit%nz(i)), & 
          (all_orbit%szp(i)), (all_orbit%itzp(i)), all_orbit%orbit_status(i)

     sum_sz = sum_sz + all_orbit%szp(i) 
     sum_tz = sum_tz + all_orbit%itzp(i) 
  ENDDO
  if(iam == 0)write(6,*)'Unoccupied j-orbits' 
  DO i= below_ef+1, below_ef+11 !tot_orbs
     if(iam == 0)WRITE(6,'((i6,1x),1x,1(g20.10,1x),1x,5(I3,1x),a)')  &
          i, all_orbit%e(i), (all_orbit%nx(i)), (all_orbit%ny(i)), (all_orbit%nz(i)), & 
          (all_orbit%szp(i)), (all_orbit%itzp(i)), all_orbit%orbit_status(i)
  ENDDO
  DO i=  1, tot_orbs
     orbitof( (all_orbit%nx(i)), (all_orbit%ny(i)), (all_orbit%nz(i)), &
          (all_orbit%szp(i)), (all_orbit%itzp(i) ) ) = i 
  ENDDO
  
  if ( iam == 0 ) write(6,*) 'Sz and Tz for vacuum state', sum_sz, sum_tz 

END SUBROUTINE setup_sp_data




SUBROUTINE precalc_chp_functions
  USE PARALLEL
  USE single_particle_orbits
  USE constants
  use chiral_potentials
  use chiral_tables 
  use chiral_constants

  implicit none 
  real*8 :: k1(3), k2(3), qtrans(3)
  integer :: p,q, ndim, m1,m2, t1,t2, nx1,ny1,nz1,nx2,ny2,nz2,i1,i2,m3,m4,m5,m6
  integer :: nx3,ny3,nz3,nxmin, nxmax
  
  allocate( sigma_dot_q_tab(-1:1,-1:1,-2*nmax:2*nmax, -2*nmax:2*nmax, -2*nmax:2*nmax) )
  sigma_dot_q_tab = 0.d0 
  allocate( sigma_dot_q_ope_tab(-1:1,-1:1,-2*nmax:2*nmax, -2*nmax:2*nmax, -2*nmax:2*nmax) )
  sigma_dot_q_ope_tab = 0.d0 
  
  nxmin =  1 
  nxmax = -1
  
  do p = 1, all_orbit%total_orbits
     
     m1 = all_orbit%szp(p) 
     t1 = all_orbit%itzp(p) 
     k1(1) = all_orbit%kx(p)
     k1(2) = all_orbit%ky(p)
     k1(3) = all_orbit%kz(p)
     
     nx1 = all_orbit%nx(p)
     ny1 = all_orbit%ny(p)
     nz1 = all_orbit%nz(p)
     
     do q = 1, all_orbit%total_orbits
        

        m2 = all_orbit%szp(q) 
        t2 = all_orbit%itzp(q) 
        k2(1) = all_orbit%kx(q)
        k2(2) = all_orbit%ky(q)
        k2(3) = all_orbit%kz(q)
        
        qtrans= hbarc*( k2-k1 )
        
        nx2 = all_orbit%nx(q)
        ny2 = all_orbit%ny(q)
        nz2 = all_orbit%nz(q)
        
        sigma_dot_q_tab(m1,m2,nx2-nx1, ny2-ny1, nz2-nz1) = chp_sigma_dot_q_mtx(m1,m2,qtrans)
        sigma_dot_q_ope_tab(m1,m2,nx2-nx1, ny2-ny1, nz2-nz1) = chp_sigma1_dot_q_ope(m1,m2,qtrans)
        
        
        nx3 = ny1*nz2-nz1*ny2
        ny3 = nz1*nx2-nx1*nz2
        nz3 = nx1*ny2-ny1*nx2
        
        if ( nx3 > nxmax ) nxmax = nx3 
        if ( nx3 < nxmin ) nxmin = nx3 
!!$        if ( ny3 > nymax ) nymax = ny3 
!!$        if ( ny3 < nymin ) nymin = ny3 
!!$        if ( nz3 > nzmax ) nzmax = nz3 
!!$        if ( nz3 < nzmin ) nzmin = nz3 
!!$        
        
     end do
  end do

  !if ( iam == 0 ) 
  !write(6,*) 'Nxmax, Nxmin', nxmin, nxmax 
  allocate( sigma_dot_qxq_tab(-1:1,-1:1, nxmin:nxmax,   nxmin:nxmax,   nxmin:nxmax) )
  sigma_dot_qxq_tab = 0.d0 
  
  do p = 1, all_orbit%total_orbits
     
     m1 = all_orbit%szp(p) 
     t1 = all_orbit%itzp(p) 
     k1(1) = all_orbit%kx(p)
     k1(2) = all_orbit%ky(p)
     k1(3) = all_orbit%kz(p)
     
     nx1 = all_orbit%nx(p)
     ny1 = all_orbit%ny(p)
     nz1 = all_orbit%nz(p)
     
     do q = 1, all_orbit%total_orbits
        

        m2 = all_orbit%szp(q) 
        t2 = all_orbit%itzp(q) 
        k2(1) = all_orbit%kx(q)
        k2(2) = all_orbit%ky(q)
        k2(3) = all_orbit%kz(q)
        
        !qtrans= hbarc*( k2-k1 )
        qtrans(1) = k1(2)*k2(3)-k1(3)*k2(2) 
        qtrans(2) = k1(3)*k2(1)-k1(1)*k2(3) 
        qtrans(3) = k1(1)*k2(2)-k1(2)*k2(1) 
        
        
        nx2 = all_orbit%nx(q)
        ny2 = all_orbit%ny(q)
        nz2 = all_orbit%nz(q)
        
        nx3 = ny1*nz2-nz1*ny2
        ny3 = nz1*nx2-nx1*nz2
        nz3 = nx1*ny2-ny1*nx2
        
        sigma_dot_qxq_tab(m1,m2,nx3, ny3, nz3) = chp_sigma_dot_q_mtx(m1,m2,qtrans)
        
     end do
  end do
  
  do m1 = -1,1,2 
     do m2 = -1,1,2 

        do m3 = -1,1,2 
           do m4 = -1,1,2 
              
              tau_dot_tau_tab(m1,m2,m3,m4) = chp_tau_dot_tau_mtx(m1,m2,m3,m4)
              do m5 = -1,1,2 
                 do m6 = -1,1,2 
                    
                    tau1_dot_tauXtau_tab(m1,m2,m3,m4,m5,m6) = chp_tau1_dot_tauXtau(m1,m2,m3,m4,m5,m6)
                    
                 end do
              end do
           end do
        end do
     end do
  end do

  delta_tab = 0.d0 
  do m1 = -1,1,2 
     do m2 = -1,1,2 
        if ( m1 == m2 ) delta_tab(m1,m2) = 1.d0
     end do
  end do
  
end SUBROUTINE precalc_chp_functions

SUBROUTINE normal_ordered_hamiltonian
  USE PARALLEL
  USE single_particle_orbits
  USE constants
  use one_body_operators
  use wave_functions
  use chiral_potentials
  use t2_storage
  use ang_mom_functions
  
  IMPLICIT NONE 
  INTEGER :: i,j,k,l,p,q,r, ia,ic,h, nx2,ny2,nz2,tz2,sz2,channel,bra,ket, ndim2, i1
  complex*16 :: sum2, sum3, sum4, norm, sum1, gmat, vint1, vint2 
  !double precision :: xx(100), wxx(100)
  real*8 :: cg1, cg2
  complex*16 :: vsum1, vsum2, vsum3, vsum4,vsum, e0_tmp, e0_3nf
  integer :: m1,m2,m3,m4,spin, ms1, ms2, a,b,c,d, ii,jj
  integer, allocatable :: ijk_confs(:), ijpq_confs(:) 
  integer ::  kk1,kk2,kk3,kk4, nconfs_tot, nconfs, number_mtxel_iam, diff, n
  integer, allocatable :: nconf_low(:), nconf_high(:)
  complex*16, allocatable :: fock_mtx_3nf(:,:), t2_tmp1(:), t2_tmp2(:)
  real*8  :: startwtime , endwtime
  
  allocate( tkin(all_orbit%total_orbits, all_orbit%total_orbits) )
  allocate( fock_mtx(all_orbit%total_orbits, all_orbit%total_orbits) )
  fock_mtx = 0.d0
  tkin = 0.d0 
  e0 = 0.d0 
  do i = 1, all_orbit%total_orbits
     tkin(i,i) = all_orbit%e(i) 
  end do
  
  

  ! 
  ! add external field to kinetic energy
  !
  !tkin = tkin + vext 
  
  fock_mtx = tkin 
  do ia = 1, all_orbit%total_orbits
     do ic = 1, all_orbit%total_orbits
        if ( all_orbit%nx(ia) /= all_orbit%nx(ic) ) cycle 
        if ( all_orbit%ny(ia) /= all_orbit%ny(ic) ) cycle
        if ( all_orbit%nz(ia) /= all_orbit%nz(ic) ) cycle
        if ( all_orbit%itzp(ia) /= all_orbit%itzp(ic) ) cycle
        
        do h = 1, below_ef
           
           nx2 = all_orbit%nx(ia) + all_orbit%nx(h)
           ny2 = all_orbit%ny(ia) + all_orbit%ny(h)
           nz2 = all_orbit%nz(ia) + all_orbit%nz(h)
           tz2 = (all_orbit%itzp(ia) + all_orbit%itzp(h))/2
                      
           if ( ia > below_ef .and. ic > below_ef ) then
              channel = locate_channel(4,tz2, nx2, ny2, nz2) 
              if ( channel == 0 ) cycle 
              
              bra = hp_hphp%ival(h,ia)
              ket = hp_hphp%ival(h,ic)
              if ( bra * ket == 0 ) cycle
              gmat = vnn_hphp(channel)%val(bra,ket) 
           elseif ( ia <= below_ef .and. ic <= below_ef ) then
              channel = locate_channel(1,tz2, nx2, ny2, nz2) 
              if ( channel == 0 ) cycle
              
              bra = hh_hhhh%ival(h,ia)
              ket = hh_hhhh%ival(h,ic)
              if ( bra * ket == 0 ) cycle
              gmat = vnn_hhhh(channel)%val(bra,ket) 
           elseif ( ia > below_ef .and. ic <= below_ef ) then
              channel = locate_channel(2,tz2, nx2, ny2, nz2) 
              if ( channel == 0 ) cycle
              
              bra = hh_hhhp%ival(h,ic)
              ket = hp_hhhp%ival(h,ia)
              if ( bra * ket == 0 ) cycle
              gmat = vnn_hhhp(channel)%val(bra,ket) 
           elseif ( ia <= below_ef .and. ic > below_ef ) then
              channel = locate_channel(2,tz2, nx2, ny2, nz2) 
              if ( channel == 0 ) cycle

              bra = hh_hhhp%ival(h,ia)
              ket = hp_hhhp%ival(h,ic)
              if ( bra * ket == 0 ) cycle
              gmat = vnn_hhhp(channel)%val(bra,ket) 
              
           end if
              
           fock_mtx(ia,ic) = fock_mtx(ia,ic) + gmat
        end DO
     end do
  end do
  

  e0= 0.d0 
  do i = 1, below_ef
     e0 = e0 + tkin(i,i)
  end do
  
  do i = 1, below_ef
     do j = 1, below_ef
        nx2 = all_orbit%nx(i) + all_orbit%nx(j)
        ny2 = all_orbit%ny(i) + all_orbit%ny(j)
        nz2 = all_orbit%nz(i) + all_orbit%nz(j)
        tz2 = (all_orbit%itzp(i) + all_orbit%itzp(j))/2
        
        channel = locate_channel(1,tz2, nx2, ny2, nz2) 
        if ( channel == 0 ) cycle
              
        bra = hh_hhhh%ival(i,j)
        ket = hh_hhhh%ival(i,j)
        if ( bra * ket == 0 ) cycle
        !e0 = e0 + 0.5d0 * ( vmom_minnesota(i,j,i,j) - vmom_minnesota(i,j,j,i) )
        e0 = e0 + 0.5d0 * vnn_hhhh(channel)%val(bra,ket) 
        
     end do
  end do
  if ( iam == 0 ) write(6,*) 'Vacuum expectation value', e0 
  
  

  if ( .not. tnf_switch ) return 
  
  !
  ! Compute 3nf vacuum expectation value 
  !
  nconfs = 0 
  do i = 1, below_ef
     do j = i+1, below_ef
        do k = j+1, below_ef
           
           nconfs= nconfs + 1 
        end do
     end do
  end do
  allocate( ijk_confs(3*nconfs) )
  
  nconfs = 0 
  ijk_confs = 0
  do i = 1, below_ef
     do j = i+1, below_ef
        do k = j+1, below_ef
           
           nconfs= nconfs + 1 
           kk1 = nconfs*3 
           kk2 = nconfs*3 - 1
           kk3 = nconfs*3 - 2
           ijk_confs(kk1) = i
           ijk_confs(kk2) = j
           ijk_confs(kk3) = k
           
        end do
     end do
  end do

  allocate( nconf_low(num_procs), nconf_high(num_procs) )
  nconf_low = 0
  nconf_high = 0

  !
  ! find number of configs on each process
  !
  nconfs_tot = nconfs
  number_mtxel_iam = floor(real(nconfs_tot/num_procs ))
  diff = nconfs - num_procs * number_mtxel_iam 
  
  n = 1 
  do i = 1, num_procs 
     if ( i < num_procs ) then
        nconf_low(i) = n 
        nconf_high(i) = n + number_mtxel_iam - 1
     else
        nconf_low(i) = n 
        nconf_high(i) = n + number_mtxel_iam - 1 + diff
     end if
     n = n + number_mtxel_iam 
  end do
  
  call mpi_barrier(mpi_comm_world,ierror)
  startwtime = MPI_WTIME()
  
  !do i = 1, num_procs
  !   if ( iam == 0 ) write(6,*) 'iam',i,'lower/upper', nconf_low(i), nconf_high(i)
  !end do
  e0_3nf = 0.d0 
  e0_tmp = 0.d0 
  
  !$omp parallel default(shared) private(ii,kk1,kk2,kk3,i,j,k,p)
  !$omp do schedule(dynamic), reduction(+:e0_tmp)
  do ii = nconf_low(iam+1), nconf_high(iam+1) 
     
     kk1 = ii*3 
     kk2 = ii*3 - 1
     kk3 = ii*3 - 2
     i = ijk_confs(kk1)
     j = ijk_confs(kk2)
     k = ijk_confs(kk3) 
     
     e0_tmp = e0_tmp + chiral_3nf_asym(i,j,k,i,j,k)
     
     !e0_tmp = e0_tmp + ( chiral_3nf(i,j,k,i,j,k) - chiral_3nf(i,j,k,j,i,k) &
     !     - chiral_3nf(i,j,k,k,j,i) - chiral_3nf(i,j,k,i,k,j) &
     !     + chiral_3nf(i,j,k,j,k,i) + chiral_3nf(i,j,k,k,i,j) )
  end do
  !$omp end do
  !$omp end parallel
  

  
  
  call mpi_reduce(e0_tmp,e0_3nf,1,mpi_complex16,mpi_sum,master, &
       mpi_comm_world,ierror)
  
  call mpi_bcast(e0_3nf,1,mpi_complex16,master, &
       mpi_comm_world,ierror)
  
  e0 = e0 + e0_3nf 
  if ( iam == 0 ) write(6,*) '3NF vacuum expectation value', e0_3nf/below_ef
  if ( iam == 0 ) write(6,*) 'Vacuum expectation value with 3nf', e0/below_ef
  
  
  
  call mpi_barrier(mpi_comm_world,ierror)
  endwtime = MPI_WTIME()
  if ( iam == 0 ) write(6,*) 'Time for e0', endwtime - startwtime
  
  if ( tnf_approx == 0 ) then
     deallocate( ijk_confs, nconf_low, nconf_high )

     !deallocate( fock_mtx, tkin , orbitof ) 
     return 
  end if
  
!!$  deallocate( ijk_confs, nconf_low, nconf_high )
!!$  deallocate( fock_mtx, tkin , orbitof ) 
!!$  
  


  !
  ! Compute 3nf normal ordered contribution to fock_mtx
  !
  
  nconfs = 0 
  do p = 1, tot_orbs
     do q = 1, tot_orbs
        if ( all_orbit%nx(p) /= all_orbit%nx(q) ) cycle 
        if ( all_orbit%ny(p) /= all_orbit%ny(q) ) cycle
        if ( all_orbit%nz(p) /= all_orbit%nz(q) ) cycle
        if ( all_orbit%itzp(p) /= all_orbit%itzp(q) ) cycle
        
        do i = 1, below_ef
           do j = i, below_ef
                            
              
              nconfs= nconfs + 1 
           end do
        end do
     end do
  end do
  
  allocate( ijpq_confs(4*nconfs) )
  
  nconfs = 0 
  ijpq_confs = 0
  do p = 1, tot_orbs
     do q = 1, tot_orbs
        if ( all_orbit%nx(p) /= all_orbit%nx(q) ) cycle 
        if ( all_orbit%ny(p) /= all_orbit%ny(q) ) cycle
        if ( all_orbit%nz(p) /= all_orbit%nz(q) ) cycle
        if ( all_orbit%itzp(p) /= all_orbit%itzp(q) ) cycle
        
        do i = 1, below_ef
           do j = i, below_ef
                            
              nconfs= nconfs + 1 
              kk1 = nconfs*4 
              kk2 = nconfs*4 - 1
              kk3 = nconfs*4 - 2
              kk4 = nconfs*4 - 3
              ijpq_confs(kk1) = i
              ijpq_confs(kk2) = j
              ijpq_confs(kk3) = p
              ijpq_confs(kk4) = q
              
           end do
        end do
     end do
  end do

  nconf_low = 0
  nconf_high = 0
  
  !
  ! find number of configs on each process
  !
  nconfs_tot = nconfs
  number_mtxel_iam = floor(real(nconfs_tot/num_procs ))
  diff = nconfs - num_procs * number_mtxel_iam 
  
  n = 1 
  do i = 1, num_procs 
     if ( i < num_procs ) then
        nconf_low(i) = n 
        nconf_high(i) = n + number_mtxel_iam - 1
     else
        nconf_low(i) = n 
        nconf_high(i) = n + number_mtxel_iam - 1 + diff
     end if
     n = n + number_mtxel_iam 
  end do
  
  !do i = 1, num_procs
  !   if ( iam == 0 ) write(6,*) 'iam',i,'lower/upper', nconf_low(i), nconf_high(i)
  !end do
  call mpi_barrier(mpi_comm_world,ierror)
  startwtime = MPI_WTIME()
  
  
  allocate( fock_mtx_3nf(tot_orbs, tot_orbs) )
  fock_mtx_3nf = 0.d0 
!  !$omp parallel do reduction (+: fock_mtx_3nf ) private (ii,kk1,kk2,kk3,kk4,i,j,p,q)
  
  !$omp parallel default(shared) private(ii,kk1,kk2,kk3,kk4,i,j,p,q)
  !$omp do schedule(dynamic), reduction(+: fock_mtx_3nf )
  do ii = nconf_low(iam+1), nconf_high(iam+1) 
     
     kk1 = ii*4 
     kk2 = ii*4 - 1
     kk3 = ii*4 - 2
     kk4 = ii*4 - 3
     i = ijpq_confs(kk1)
     j = ijpq_confs(kk2)
     p = ijpq_confs(kk3) 
     q = ijpq_confs(kk4) 
     
     
     fock_mtx_3nf(p,q) = fock_mtx_3nf(p,q) + chiral_3nf_asym(i,j,p,i,j,q) 
     !fock_mtx_3nf(p,q) = fock_mtx_3nf(p,q) + ( chiral_3nf(i,j,p,i,j,q) - chiral_3nf(i,j,p,j,i,q) &
     !     - chiral_3nf(i,j,p,q,j,i) - chiral_3nf(i,j,p,i,q,j) &
     !     + chiral_3nf(i,j,p,j,q,i) + chiral_3nf(i,j,p,q,i,j) )
  end do
  !$omp end do
  !$omp end parallel
  
  call mpi_barrier(mpi_comm_world,ierror)
  endwtime = MPI_WTIME()
  if ( iam == 0 ) write(6,*) 'Time for f3nf', endwtime - startwtime
  
  ndim2 = tot_orbs**2
  allocate( t2_tmp1(ndim2) )
  allocate( t2_tmp2(ndim2) )
  t2_tmp1 = 0.d0
  t2_tmp2 = 0.d0

  i1 = 0
  do p = 1, tot_orbs
     do q = 1, tot_orbs
        i1 = i1 + 1
        t2_tmp1(i1) = fock_mtx_3nf(p,q)
     end do
  end do
  
  call mpi_reduce(t2_tmp1,t2_tmp2,ndim2,mpi_complex16,mpi_sum, &
       master,mpi_comm_world,ierror)
  t2_tmp1 = 0.d0
  t2_tmp1 = t2_tmp2
  call mpi_bcast(t2_tmp1,ndim2,mpi_complex16,master, &
       mpi_comm_world,ierror)

  fock_mtx_3nf = 0.d0 
  i1 = 0
  do p = 1, tot_orbs
     do q = 1, tot_orbs
        i1 = i1 + 1
        fock_mtx_3nf(p,q) = t2_tmp1(i1)
     end do
  end do
  deallocate( t2_tmp1, t2_tmp2 )
  
  e0_tmp = 0.d0 
  do i = 1, below_ef
     e0_tmp = e0_tmp + fock_mtx_3nf(i,i)/3.d0 
  end do
  if ( iam == 0 ) write(6,*) 'e0 3nf', e0_tmp, e0_3nf 
  
  
  fock_mtx = fock_mtx + fock_mtx_3nf 
  
  
  if ( tnf_approx == 1 ) then
     deallocate( fock_mtx_3nf ) 
     deallocate( ijk_confs )
     deallocate( ijpq_confs )
  
     return 
  end if
  
  !  deallocate( fock_mtx, tkin, orbitof ) 
  
  call setup_3nf_2bNO
  deallocate( fock_mtx_3nf ) 
  deallocate( ijk_confs )
  deallocate( ijpq_confs )
  
  
END SUBROUTINE normal_ordered_hamiltonian



SUBROUTINE setup_3nf_2bNO
  USE single_particle_orbits
  USE constants
  use one_body_operators
  use t2_storage
  use configurations
  USE KSPACE 
  use CHIRAL_POTENTIALS

  IMPLICIT NONE
  
  INTEGER :: i,j,k,l,p,a,b,c,d,dim1,dim2, number_channels, nx,ny,nz,sz,tz,  sumx, sumy, sumz  
  INTEGER :: k1,k2,k3,k4,k5, ab_confs, ij_confs, bra,ket, channel
  INTEGER :: nxa, nya, nza, nxb, nyb, nzb,tza,tzb,sza,szb, ichan   
  INTEGER :: numchannels(1:6), channel1, channel2, channel3, channel4, channel5, channel6
  INTEGER :: ket_min, ket_max, bra_min, bra_max, i1, local_ch_id, ndim, bra2, ket2, chunk
  complex*16, allocatable :: t2_tmp1(:), t2_tmp2(:), amat(:,:)
  real*8  :: startwtime , endwtime
  complex*16 :: vsum, v3nf 
  
  startwtime = MPI_WTIME()
  
  
  if ( iam == 0 ) write(6,*) 'Setting up v3nf_2bNO for hhhh block'
  do channel   = 1, channels%number_hhhh_confs 
     nx = channels%hhhh_quantum_numbers(channel*4)
     ny = channels%hhhh_quantum_numbers(channel*4-1)
     nz = channels%hhhh_quantum_numbers(channel*4-2)
     tz = channels%hhhh_quantum_numbers(channel*4-3)
     
     
     dim1 = size( vnn_hhhh(channel)%val, 1)
     allocate( amat( dim1, dim1) ) 
     amat  = 0.d0 
     
     !$omp parallel default(shared) private(bra,i,j,ket,k,l,p,v3nf)
     !$omp do schedule(dynamic), reduction(+:amat)
     do bra = 1, size(  lookup_hhhh_configs(1,channel)%ival, 2) 
        i= lookup_hhhh_configs(1,channel)%ival(1,bra) 
        j = lookup_hhhh_configs(1,channel)%ival(2,bra) 
        
        if ( j > i ) cycle 
        do ket = bra, size(  lookup_hhhh_configs(1,channel)%ival, 2) 
           k = lookup_hhhh_configs(1,channel)%ival(1,ket)
           l = lookup_hhhh_configs(1,channel)%ival(2,ket) 
           if ( l > k ) cycle 
           
           !vsum = 0.d0 
           do p = 1, below_ef
              !v3nf =  ( chiral_3nf(i,j,p,k,l,p) - chiral_3nf(i,j,p,l,k,p) &
              !        - chiral_3nf(i,j,p,p,l,k) - chiral_3nf(i,j,p,k,p,l) &
              !        + chiral_3nf(i,j,p,l,p,k) + chiral_3nf(i,j,p,p,k,l) )
              v3nf =  chiral_3nf_asym(i,j,p,k,l,p) 

              amat(bra,ket) = amat(bra,ket) + v3nf 
              if ( bra /= ket ) amat(ket,bra) = amat(ket,bra) + conjg(v3nf) 
              
           end do
           
     
        end do
     end do
     !$omp end do
     !$omp end parallel
     
     do bra = 1, size(  lookup_hhhh_configs(1,channel)%ival, 2) 
        i= lookup_hhhh_configs(1,channel)%ival(1,bra) 
        j = lookup_hhhh_configs(1,channel)%ival(2,bra) 
        
        if ( j > i ) cycle 
        do ket = 1, size(  lookup_hhhh_configs(1,channel)%ival, 2) 
           k = lookup_hhhh_configs(1,channel)%ival(1,ket)
           l = lookup_hhhh_configs(1,channel)%ival(2,ket) 
           if ( l > k ) cycle 
           
           bra2 = hh_hhhh%ival(j,i)
           amat(bra2,ket) = -amat(bra,ket) 
           ket2 = hh_hhhh%ival(l,k)
           amat(bra,ket2) = -amat(bra,ket) 
           amat(bra2,ket2) = amat(bra,ket) 
        end do
     end do
     
     vnn_hhhh(channel)%val = vnn_hhhh(channel)%val + amat 
     deallocate( amat ) 
  end do
  
  endwtime = MPI_WTIME()
  if ( iam == 0 ) write(6,*) 'Time to setup v3nf_2bNO hhhh block', endwtime - startwtime
  startwtime = MPI_WTIME()
  
  allocate( v3nf_hhpp(channels%number_hhpp_confs) )
  do channel   = 1, channels%number_hhpp_confs 
     dim1 = size(  lookup_hhpp_configs(1,channel)%ival, 2) 
     dim2 = size(  lookup_hhpp_configs(2,channel)%ival, 2) 
     allocate( v3nf_hhpp(channel)%val(dim1, dim2) )
     v3nf_hhpp(channel)%val = 0.d0 
  end do
  
  !
  ! setup hhpp NN interaction
  !
  if ( iam == 0 ) write(6,*) 'Setting up v3nf_2bNO for hhpp block'
  do channel   = 1, channels%number_hhpp_confs 
     nx = channels%hhpp_quantum_numbers(channel*4)
     ny = channels%hhpp_quantum_numbers(channel*4-1)
     nz = channels%hhpp_quantum_numbers(channel*4-2)
     tz = channels%hhpp_quantum_numbers(channel*4-3)
     
     if ( check_my_channel_hhpp(channel) == 0 ) cycle
     local_ch_id = channel
     !
     ! ket side if fully stored on each proc
     !
     bra_min = mapping_hhpp(iam+1,local_ch_id,2)
     bra_max = mapping_hhpp(iam+1,local_ch_id,3)
     !
     ! bra side is distributed 
     !
     ket_min = mapping_hhpp(iam+1,local_ch_id,4)
     ket_max = mapping_hhpp(iam+1,local_ch_id,5)
     
     allocate( amat(size(v3nf_hhpp(channel)%val,1), size(v3nf_hhpp(channel)%val,2) ) ) 
     amat = 0.d0 
     !$omp parallel default(shared) private(bra,i,j,ket,k,l,p,v3nf)
     !$omp do schedule(dynamic), reduction(+:amat)
     do ket = bra_min, bra_max !1, size(  lookup_hhpp_configs(2,channel)%ival, 2) 
        !do ket = 1, size(  lookup_hhpp_configs(2,channel)%ival, 2) 
        
        i = lookup_hhpp_configs(2,channel)%ival(1,ket)
        j = lookup_hhpp_configs(2,channel)%ival(2,ket) 
        do bra = 1, size(  lookup_hhpp_configs(1,channel)%ival, 2) 
           k = lookup_hhpp_configs(1,channel)%ival(1,bra) 
           l = lookup_hhpp_configs(1,channel)%ival(2,bra) 
           
           !vsum = 0.d0 
           do p = 1, below_ef
              !v3nf =  ( chiral_3nf(i,j,p,k,l,p) - chiral_3nf(i,j,p,l,k,p) &
              !        - chiral_3nf(i,j,p,p,l,k) - chiral_3nf(i,j,p,k,p,l) &
              !        + chiral_3nf(i,j,p,l,p,k) + chiral_3nf(i,j,p,p,k,l) )
              v3nf =  chiral_3nf_asym(i,j,p,k,l,p) 
              
              amat(bra,ket) = amat(bra,ket) + v3nf
              !vsum = vsum + v3nf 
           end do
           
           !v3nf_hhpp(channel)%val(bra,ket) = vsum 
           
        end do
     end do
     !$omp end do
     !$omp end parallel
     
     v3nf_hhpp(channel)%val = v3nf_hhpp(channel)%val + amat 
     deallocate( amat ) 
     
  end do
  
  number_channels = size( v3nf_hhpp )
  do channel = 1, number_channels 
     bra = size(v3nf_hhpp(channel)%val,1)
     ket = size(v3nf_hhpp(channel)%val,2)
     
     ndim = bra * ket
     allocate( t2_tmp1(ndim) )
     allocate( t2_tmp2(ndim) )
     t2_tmp1 = 0.d0
     t2_tmp2 = 0.d0
     
     i1 = 0 
     do i =1, bra
        do j = 1, ket
           i1 = i1 + 1
           t2_tmp1(i1) = v3nf_hhpp(channel)%val(i,j)
        end do
     end do
     
     call mpi_reduce(t2_tmp1,t2_tmp2,ndim,mpi_complex16,mpi_sum, &
          master,mpi_comm_world,ierror)
     t2_tmp1 = 0.d0 
     t2_tmp1 = t2_tmp2 
     call mpi_bcast(t2_tmp1,ndim,mpi_complex16,master, &
          mpi_comm_world,ierror)
     
     i1 = 0 
     do i =1, bra
        do j = 1, ket
           i1 = i1 + 1
           v3nf_hhpp(channel)%val(i,j) = t2_tmp1(i1)
        end do
     end do
     
     deallocate( t2_tmp1, t2_tmp2 ) 
  end do

  do channel   = 1, channels%number_hhpp_confs 
     vnn_hhpp(channel)%val = vnn_hhpp(channel)%val + v3nf_hhpp(channel)%val 
  end do
  
  do channel   = 1, channels%number_hhpp_confs 
     deallocate ( v3nf_hhpp(channel)%val )
  end do
  deallocate( v3nf_hhpp ) 
  
  endwtime = MPI_WTIME()
  if ( iam == 0 ) write(6,*) 'Time to setup v3nf_2bNO hhpp block', endwtime - startwtime
  startwtime = MPI_WTIME()
  
  !
  ! setup hhpp NN interaction
  !
  if ( iam == 0 ) write(6,*) 'Setting up v3nf_2bNO for hphp block'
  allocate( v3nf_hphp(channels%number_hphp_confs) )
  do channel   = 1, channels%number_hphp_confs 
     dim1 = size(  lookup_hphp_configs(1,channel)%ival, 2) 
     allocate( v3nf_hphp(channel)%val(dim1, dim1) )
     v3nf_hphp(channel)%val = 0.d0 
  end do
  
  do channel   = 1, channels%number_hphp_confs 
     nx = channels%hphp_quantum_numbers(channel*4)
     ny = channels%hphp_quantum_numbers(channel*4-1)
     nz = channels%hphp_quantum_numbers(channel*4-2)
     tz = channels%hphp_quantum_numbers(channel*4-3)
     
     if ( check_my_channel_hphp(channel) == 0 ) cycle
     local_ch_id = channel
     !
     ! ket side if fully stored on each proc
     !
     bra_min = mapping_hphp(iam+1,local_ch_id,2)
     bra_max = mapping_hphp(iam+1,local_ch_id,3)
     !
     ! bra side is distributed 
     !
     ket_min = mapping_hphp(iam+1,local_ch_id,4)
     ket_max = mapping_hphp(iam+1,local_ch_id,5)
     
     allocate( amat(size(v3nf_hphp(channel)%val,1), size(v3nf_hphp(channel)%val,2) ) ) 
     amat = 0.d0 
     !$omp parallel default(shared) private(bra,i,j,ket,k,l,p,v3nf)
     !$omp do schedule(dynamic), reduction(+:amat)
     !do bra = bra_min, bra_max !1, size(  lookup_hhpp_configs(2,channel)%ival, 2) 
     do bra = 1, size(  lookup_hphp_configs(1,channel)%ival, 2) 
        i= lookup_hphp_configs(1,channel)%ival(1,bra) 
        j = lookup_hphp_configs(1,channel)%ival(2,bra) 
        
        do ket = bra, size(  lookup_hphp_configs(1,channel)%ival, 2) 
           k = lookup_hphp_configs(1,channel)%ival(1,ket)
           l = lookup_hphp_configs(1,channel)%ival(2,ket) 
           
           !vsum = 0.d0 
           do p = 1, below_ef
              !v3nf =  ( chiral_3nf(i,j,p,k,l,p) - chiral_3nf(i,j,p,l,k,p) &
              !        - chiral_3nf(i,j,p,p,l,k) - chiral_3nf(i,j,p,k,p,l) &
              !        + chiral_3nf(i,j,p,l,p,k) + chiral_3nf(i,j,p,p,k,l) )
              v3nf =  chiral_3nf_asym(i,j,p,k,l,p) 
              
              amat(bra,ket) = amat(bra,ket) + v3nf 
              if ( ket /= bra ) amat(ket,bra) = amat(ket,bra) + conjg(v3nf)
              !vsum = vsum + v3nf 
           end do
           !v3nf_hphp(channel)%val(bra,ket) = v3nf_hphp(channel)%val(bra,ket) + vsum 
           
        end do
     end do
     !$omp end do
     !$omp end parallel
     
     v3nf_hphp(channel)%val(bra_min:bra_max,ket_min:ket_max) = & 
          v3nf_hphp(channel)%val(bra_min:bra_max, ket_min:ket_max) + amat(bra_min:bra_max, ket_min:ket_max)
     deallocate( amat ) 
  end do
  
  number_channels = size( v3nf_hphp )
  do channel = 1, number_channels 
     bra = size(v3nf_hphp(channel)%val,1)
     ket = size(v3nf_hphp(channel)%val,2)
     
     ndim = bra * ket
     allocate( t2_tmp1(ndim) )
     allocate( t2_tmp2(ndim) )
     t2_tmp1 = 0.d0
     t2_tmp2 = 0.d0

     i1 = 0 
     do i =1, bra
        do j = 1, ket
           i1 = i1 + 1
           t2_tmp1(i1) = v3nf_hphp(channel)%val(i,j)
        end do
     end do
     
     call mpi_reduce(t2_tmp1,t2_tmp2,ndim,mpi_complex16,mpi_sum, &
          master,mpi_comm_world,ierror)
     t2_tmp1 = 0.d0 
     t2_tmp1 = t2_tmp2 
     call mpi_bcast(t2_tmp1,ndim,mpi_complex16,master, &
          mpi_comm_world,ierror)
     
     i1 = 0 
     do i =1, bra
        do j = 1, ket
           i1 = i1 + 1
           v3nf_hphp(channel)%val(i,j) = t2_tmp1(i1)
        end do
     end do
     
     deallocate( t2_tmp1, t2_tmp2 ) 
  end do
  
  do channel   = 1, channels%number_hphp_confs 
     vnn_hphp(channel)%val = vnn_hphp(channel)%val + v3nf_hphp(channel)%val 
  end do
  
  do channel   = 1, channels%number_hphp_confs 
     deallocate ( v3nf_hphp(channel)%val ) 
  end do
  deallocate ( v3nf_hphp )
  
  endwtime = MPI_WTIME()
  if ( iam == 0 ) write(6,*) 'Time to setup v3nf_2bNO hphp block', endwtime - startwtime
  startwtime = MPI_WTIME()

  
  !
  ! setup pppp NN interaction
  !
  if ( iam == 0 ) write(6,*) 'Setting up v3nf_NO for pppp block'
  !do channel   = 1, channels%number_pppp_confs 
  do channel   = my_channel_low(iam), my_channel_high(iam) !channels%number_pppp_confs 
     nx = channels%hhpp_quantum_numbers(channel*4)
     ny = channels%hhpp_quantum_numbers(channel*4-1)
     nz = channels%hhpp_quantum_numbers(channel*4-2)
     tz = channels%hhpp_quantum_numbers(channel*4-3)
     
     !
     ! check that pppp channel is restricted to hhpp channels. 
     !
     channel2 = channel !locate_channel(3,tz, nx, ny, nz)
     !if ( channel2 == 0 ) cycle 
     !if ( check_my_channel(channel2) == 0 ) cycle
     local_ch_id = channel2
     !if ( channel2 /= channel ) write(6,*) 'error in channel' 
     
     !
     ! ket side if fully stored on each proc
     !
     ket_min = mapping(iam+1,local_ch_id,2)
     ket_max = mapping(iam+1,local_ch_id,3)
     !
     ! bra side is distributed 
     !
     bra_min = mapping(iam+1,local_ch_id,4)
     bra_max = mapping(iam+1,local_ch_id,5)
     
     !allocate( amat(bra_min:bra_max, ket_min:ket_max) )
     
     !chunk = nint( size( lookup_pppp_configs(1,channel)%ival, 2)/4.)
     
     allocate( amat(size(  lookup_pppp_configs(1,channel)%ival, 2), & 
          size(  lookup_pppp_configs(1,channel)%ival, 2) ) )
     
     amat = 0.d0 
!     ! $omp parallel do schedule(dynamic), reduction (+: amat ) private (bra,i,j,ket,k,l,p,v3nf)
     !$omp parallel default(shared) private(bra,i,j,ket,k,l,p,v3nf)
     !$omp do schedule(dynamic)
     
     !, reduction(+: amat )
     !do bra = bra_min, bra_max !1, size(  lookup_pppp_configs(1,channel)%ival, 2) 
     do bra = 1, size(  lookup_pppp_configs(1,channel)%ival, 2) 
        i = lookup_pppp_configs(1,channel)%ival(1,bra) 
        j = lookup_pppp_configs(1,channel)%ival(2,bra) 
        
        if ( j > i ) cycle 
        
        do ket = bra, size(  lookup_pppp_configs(1,channel)%ival, 2) 
           k = lookup_pppp_configs(1,channel)%ival(1,ket)
           l = lookup_pppp_configs(1,channel)%ival(2,ket) 
           
           if ( l > k ) cycle 
           do p = 1, below_ef
              
              !v3nf = 0.d0 
              v3nf =  chiral_3nf_asym(i,j,p,k,l,p) 
              amat(bra,ket) = amat(bra,ket) + v3nf 
              if ( ket /= bra ) amat(ket,bra) = amat(ket,bra) + conjg(v3nf)
              
           end do
                      
        end do
     end do
     !$omp end do
     !$omp end parallel
     
     do bra = 1, size(  lookup_pppp_configs(1,channel)%ival, 2) 
        i = lookup_pppp_configs(1,channel)%ival(1,bra) 
        j = lookup_pppp_configs(1,channel)%ival(2,bra) 
        if ( j > i ) cycle 
        
        do ket = 1, size(  lookup_pppp_configs(1,channel)%ival, 2) 
           k = lookup_pppp_configs(1,channel)%ival(1,ket)
           l = lookup_pppp_configs(1,channel)%ival(2,ket) 
           if ( l > k ) cycle 
                      
           bra2 = pp_pppp%ival(j,i)
           amat(bra2,ket) = -amat(bra,ket) 
           ket2 = pp_pppp%ival(l,k)
           amat(bra,ket2) = -amat(bra,ket) 
           amat(bra2,ket2) = amat(bra,ket) 
        end do
     end do
     
     vnn_pppp(channel2)%val(bra_min:bra_max,:) = vnn_pppp(channel2)%val(bra_min:bra_max,:) + & 
          amat(bra_min:bra_max,:)
     deallocate( amat )
  end do
  
  endwtime = MPI_WTIME()
  if ( iam == 0 ) write(6,*) 'Time to setup v3nf_2bNO pppp block', endwtime - startwtime
end SUBROUTINE setup_3nf_2bNO



! jobvl, jobvr, np, temp1, lda, omega2_eig, vl, ldvl, eigen_vec, ldvr, &
!       work1, lwork, rwork, info
SUBROUTINE lapack_diag(h, cvec, ceig, n )
  

  implicit none
  integer, intent(in) :: n
  complex*16, intent(in) :: h(n,n)
  COMPLEX*16, intent(out) :: cvec(n,n)
  COMPLEX*16 ::  vl(n,n)
  COMPLEX*16, intent(out) :: ceig(n)
  DOUBLE PRECISION :: rwork(2*n)
  COMPLEX*16, DIMENSION(10000) :: work1
  INTEGER :: lda, ldb, ldvl, ldvr, info, lwork, ilo , ihi
  CHARACTER*1 :: jobvl, jobvr
  DOUBLE PRECISION, DIMENSION(n) :: scale, rconde, rcondv
  complex*16 :: norm
  integer :: i, j

  jobvl = 'N' ;  jobvr = 'V';  lda = n
  ldvl = 1;  ldvr = n;  lwork = 10000
  ceig = 0.; cvec = 0.
  
  CALL zgeev( jobvl, jobvr, n, h, lda, ceig, vl, ldvl, cvec, ldvr, &
       work1, lwork, rwork, info )
  
  ! berggren normalization
  do i = 1, n
     norm = sum( cvec(:,i)*cvec(:,i) )
     cvec(:,i) = cvec(:,i)/sqrt(norm)
     !write(6,*) i, norm

  end do
  
  do i = 1, n
     do j = 1, n
!        if ( i /=j .and. abs(sum( (cvec(:,i))*cvec(:,j) ) ) > 1.e-5 ) then
!           write(6,*)  i,j, sum( (cvec(:,i))*cvec(:,j) ) 
!           !stop
!        end if
     end do
  end do
  
  call eig_sort(ceig,cvec, n)

end SUBROUTINE lapack_diag

!
! eigenvalue sort
! sort cmplx vector real(a(1))<real(a(2)) < ... < real(a(n))
!
SUBROUTINE eig_sort(a,b, n)
  IMPLICIT NONE
  
  INTEGER :: i, j, n
  complex*16, DIMENSION(n), INTENT(INOUT) :: a
  complex*16, DIMENSION(n,n), INTENT(INOUT) :: b
  complex*16 :: temp1, temp2
  complex*16, DIMENSION(n) :: temp3
  
  DO i = 1, n
     DO j = 1, n
        IF ( real( a(i) )  < real( a(j) ) ) THEN
           
           temp1 = a(i)
           a(i) = a(j) 
           a(j) = temp1
           
           temp3(:) = b(:,i)
           b(:,i) = b(:,j) 
           b(:,j) = temp3(:)
           


        END IF
     END DO
  END DO

END SUBROUTINE eig_sort


!*****7**-*********-*********-*********-*********-*********-*********-72
      SUBROUTINE index(n,arr,indx)
! integer array arr is indexed in ascending order
      INTEGER :: n,indx(*),M,NSTACK
      real*8 :: arr(*)
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,indxt,ir,itemp,j,jstack,k,l,istack(NSTACK)
      real*8 :: a
      do 11 j=1,n
        indx(j)=j
11    continue
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 13 j=l+1,ir
          indxt=indx(j)
          a=arr(indxt)
          do 12 i=j-1,1,-1
            if(arr(indx(i)).le.a)goto 2
            indx(i+1)=indx(i)
12        continue
          i=0
2         indx(i+1)=indxt
13      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        itemp=indx(k)
        indx(k)=indx(l+1)
        indx(l+1)=itemp
        if(arr(indx(l+1)).gt.arr(indx(ir)))then
          itemp=indx(l+1)
          indx(l+1)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l)).gt.arr(indx(ir)))then
          itemp=indx(l)
          indx(l)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l+1)).gt.arr(indx(l)))then
          itemp=indx(l+1)
          indx(l+1)=indx(l)
          indx(l)=itemp
        endif
        i=l+1
        j=ir
        indxt=indx(l)
        a=arr(indxt)
3       continue
          i=i+1
        if(arr(indx(i)).lt.a)goto 3
4       continue
          j=j-1
        if(arr(indx(j)).gt.a)goto 4
        if(j.lt.i)goto 5
        itemp=indx(i)
        indx(i)=indx(j)
        indx(j)=itemp
        goto 3
5       indx(l)=indx(j)
        indx(j)=indxt
        jstack=jstack+2
        if(jstack.gt.NSTACK)pause 'NSTACK too small in indexx'
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      END subroutine index
!*****7**-*********-*********-*********-*********-*********-*********-72





!
! check one-body basis in momentum/position space
! 
subroutine diagonalize_h_osc_basis
  Use constants
  use single_particle_orbits
  USE wave_functions
  use ang_mom_functions

  
  implicit none
  complex*16, allocatable  :: ham(:,:), tkin(:,:), ceig(:), cvec(:,:) 
  real*8 :: int_sum, vrspace_minnesota
  real*8, allocatable, dimension(:,:) :: osc_funcs
  double precision, allocatable, dimension(:) :: rr, wrr
  integer :: a, na,la,ja, nc,tza, lang, jang, iso, i,j, n_max, nn
  double precision :: factor, yukawa, eff_mass, ekin1, ph, alpha, charge, mu
  double precision :: cx(0:200), xp, zz, dz
  real*8 :: oscl_r, osc_a, osc_c, zlab1, zlab2, zlab, sum_rel, e_ho
  integer :: nosc_max, n, np, lb
  
  hbar_omega = 20.d0
  oscl=hbarc/SQRT(p_mass*hbar_omega)
  oscl_r=oscl*SQRT(2.)

  nosc_max = 50
  ! specify quantum numbers
  lang = 0; jang =1; iso = 1; lb = 0
  nn = 250
  allocate( rr(nn), wrr(nn) )
  call gauss_legendre(0.d0, 20.d0, rr, wrr, nn)
  
  write(6,*) hbar_omega, oscl_r 
  allocate( osc_funcs(0:nosc_max, nn) )
  allocate( cvec(nosc_max+1, nosc_max+1), ceig(nosc_max+1) )
  allocate( tkin(nosc_max+1, nosc_max+1),  ham(nosc_max+1, nosc_max+1) ) 
  tkin = 0.d0 
  ham = 0.d0 
  cvec = 0.d0 
  ceig = 0.d0 
  
  osc_funcs = 0.d0 
  la = lang  
  do na = 0, nosc_max, 1 
     
     ph = 1.d0 
     factor = 0.5D0*((na+la+2)*LOG(2.D0)+fac(na)-dfac(2*na+2*la+1)-0.5D0*LOG(pi))
     factor = EXP(factor)
     sum_rel=0.
     DO i=1,nn
        
        zlab= rr(i)/oscl_r
        
        CALL laguerre_general( na, la+0.5D0, zlab*zlab, cx )
        xp = cx(na)*exp(-zlab*zlab*0.5d0)*(zlab**la)
        osc_funcs(na,i) = xp*factor*(1.d0/oscl_r**(1.5D0))
        sum_rel = sum_rel + osc_funcs(na, i)**2 * rr(i)**2*wrr(i)
        
     end DO
     write(6,*) 'Norm', sum_rel 
     
  end do
  
  !
  ! kinetic energy
  !
  tkin = 0.d0
  do n = 0,  nosc_max
     do np = 0,  nosc_max
        
        e_ho = 0.d0
        IF ( n == np-1) e_ho = dSQRT(np*(np+lb+0.5d0)) 
        IF ( n == np)   e_ho =  2*np+lb+1.5 ! harmonic oscillator diagonal energy
        IF ( n == np+1) e_ho = dSQRT((np+1.d0)*(np+lb+1.5d0))
        tkin(n+1,np+1) = e_ho*hbar_omega/2.d0
                
     end do
  end do
  
  !
  ! potential setup
  !
  ham = 0.d0
  do na = 0, nosc_max, 1
     do nc = 0,  nosc_max, 1
        
        ! r-space
        int_sum = 0.d0
        do i = 1, nn
           
           zlab = dcmplx( rr(i),0.d0 )
                   
           int_sum = int_sum + wrr(i)*rr(i)**2 * osc_funcs(na,i) * & 
                osc_funcs(nc,i)  * vrspace_minnesota(rr(i))
           
        end do
        ham(na+1,nc+1) =  int_sum + tkin(na+1,nc+1)
        
        
     end do
  end do
  
  
  call lapack_diag( ham, cvec, ceig, nosc_max+1 )

  
  do i = 1, 10 ! nosc_max+1
     write(6,*) dble(ceig(i)), aimag(ceig(i))
  end do
  
end subroutine diagonalize_h_osc_basis

