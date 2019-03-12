SUBROUTINE ccd_iter
  USE diis_mod
  USE PARALLEL
  USE CONSTANTS
  USE one_body_operators
  USE t2_storage
  
  IMPLICIT NONE
  INTEGER :: count, itimes, ntimes, channel, bra, ket, i
  complex*16 :: ener1, ener2, dener
  logical :: switch 
  real*8  ::  startwtime , endwtime
  
  !
  ! Setup initial t1 and t2 amplitudes
  !
  
  call ccd_energy_save(ener1 ) 
  mbpt2 = e0 + ener1 
  

  !
  ! read in subspace dimension and at what step diis is called
  !
  !read(5,*);read(5,*) subspace, diis_step
  if ( cc_approx == 'CCDT1' ) then 
     
     call setup_vnn_hppp_block
     !call setup_t3_channel_structures
     call setup_t3full_channel_structures
     call mpi_barrier(mpi_comm_world, ierror)
     !call setup_proc_t3_mappings
     if ( iam == 0 ) write(6,*) 'Setting up t3 and v3' 
     call setup_t3
     call mpi_barrier(mpi_comm_world, ierror)
     if ( iam == 0 ) write(6,*) 'Setting up t3 and v3 done' 
     !call t3_intermediate
     call t3_energy(ener1) 
     
     
  end if
  
  
  subspace = 6
  diis_step = 6 
  
  
  
  
  call diis_setup
  !
  ! Allocate energy denominators
  !
  allocate( d1ik(1:below_ef), d1ac(below_ef+1:tot_orbs) )
  allocate( t1_denom(below_ef+1:tot_orbs, 1:below_ef) )
  allocate( d2ij(1:below_ef), d2ab(below_ef+1:tot_orbs) )
  t1_denom = 0.d0
  D1IK = 0.D0; D1AC = 0.D0
  d2ij = 0.d0; d2ab = 0.d0
  
  
  !call setup_v3nf_channel_structures
  
  !
  ! here is the main iteration loop
  !
  startwtime = MPI_WTIME()
  nstep=0
  count = 0
  ener2 = 0.d0
  ntimes =  100
  dener=4.0
  do itimes=1, ntimes
     if(abs(dener) > 1.0e-6 )then
        
        !if ( iam == 0 ) write(6,*) 't2 eqn'
        
        if ( cc_approx == 'CCDT1' ) call t3_intermediate
        call t2_intermediate(0)
        
        
        
        count = count + 1
        nstep = nstep + 1
        
        !call linear_t2(0.d0) !
        call diis_t1t2(count)
        call ccd_energy_save(ener2)
        if ( cc_approx == 'CCDT1' ) call t3_energy(ener2)
        
        dener=ener2-ener1
        ener1=ener2
        
        if(iam ==0)write(6,'(a,2i5,2e16.8)')'dener',iam,itimes,real(dener) 
     end if
  end do
  
!!$  if ( iam == 0 ) write(6,*) '------------------------'
!!$  if ( iam == 0 ) write(6,*) 'Adding T2V3 contribution to T2'
!!$  
!!$  nstep=0
!!$  count = 0
!!$  ener2 = 0.d0
!!$  ntimes =  80
!!$  dener=4.0
!!$  do itimes=1, ntimes
!!$     if(abs(dener) > 1.0e-6 )then
!!$        
!!$        !if ( iam == 0 ) write(6,*) 't2 eqn'
!!$        
!!$        if ( cc_approx == 'CCDT1' ) call t3_intermediate
!!$        call t2_intermediate(1)
!!$        
!!$        
!!$        count = count + 1
!!$        nstep = nstep + 1
!!$        
!!$        !call linear_t2(0.d0) !
!!$        call diis_t1t2(count)
!!$        call ccd_energy_save(ener2)
!!$        if ( cc_approx == 'CCDT1' ) call t3_energy(ener2)
!!$        
!!$        dener=ener2-ener1
!!$        ener1=ener2
!!$        
!!$        if(iam ==0)write(6,'(a,2i5,2e16.8)')'dener',iam,itimes,real(dener) 
!!$     end if
!!$  end do
!!$  
  
  call diis_take_down
  
  eccsd = e0 + ener2 
  endwtime = MPI_WTIME()
  if ( iam == 0 ) write(6,*) 'Total execution time CCD', endwtime - startwtime
  
  IF ( cc_approx == 'CCD(T)') then 
     if ( iam == 0 ) write(*,*) "preparing to do CCD(T): "
     do i = 1, size( vnn_hhpp ) 
        deallocate( vnn_hhpp(i)%val )
     end do
     do i = 1, size( vnn_hphp ) 
        deallocate( vnn_hphp(i)%val )
     end do
     deallocate( vnn_hhpp ) 
     deallocate( vnn_hphp ) 
     
     do i = 1, size( t2_ccm )
        deallocate( t2_ccm_eqn(i)%val )
     end do
     deallocate( t2_ccm_eqn )
     
     do i = 1, size( t2_ph_ccm) 
        deallocate( t2_ph_ccm(i)%val )
        deallocate( t2_ph_ccm_eqn(i)%val )
     end do
     deallocate( t2_ph_ccm, t2_ph_ccm_eqn )
     
     call setup_vnn_hppp_block
     call setup_t3_channel_structures
     call pert_triples2
     
  end if

  if(iam==0) write(6,*) "Finished with ccd_iter"

  deallocate( d1ik, d1ac )
  deallocate( t1_denom )
  deallocate( d2ij, d2ab )
END SUBROUTINE ccd_iter

!
! Setup two-body channels and calculate and store nucleon-nucleon interaction
!
SUBROUTINE setup_channel_structures
  USE single_particle_orbits
  USE constants
  use one_body_operators
  use t2_storage
  use configurations
  USE KSPACE 
  use CHIRAL_POTENTIALS

  IMPLICIT NONE
  
  INTEGER :: i,j,k,l,a,b,c,d,dim1,dim2, number_channels, nx,ny,nz,sz,tz,  sumx, sumy, sumz  
  INTEGER :: k1,k2,k3,k4,k5, ab_confs, ij_confs, bra,ket, channel
  real*8 :: memory,  denom
  integer, allocatable :: numstates(:,:,:,:,:)
  INTEGER :: nxa, nya, nza, nxb, nyb, nzb,tza,tzb,sza,szb, ichan   
  INTEGER :: numchannels(1:6), channel1, channel2, channel3, channel4, channel5, channel6
  INTEGER, ALLOCATABLE :: numchannels1(:), numchannels2(:), numchannels3(:), & 
       numchannels4(:), numchannels5(:), numchannels6(:)
  INTEGER :: ket_min, ket_max, bra_min, bra_max, i1, local_ch_id, ndim
  complex*16, allocatable :: t2_tmp1(:), t2_tmp2(:) 
  real*8  :: startwtime , endwtime
  

  nx_min = 1000
  nx_max = -1000
  ny_min = 1000
  ny_max = -1000
  nz_min = 1000
  nz_max = -1000
  do i = 1, tot_orbs
     do j = 1, tot_orbs
        
        sumx = all_orbit%nx(i) +  all_orbit%nx(j) 
        sumy = all_orbit%ny(i) +  all_orbit%ny(j) 
        sumz = all_orbit%nz(i) +  all_orbit%nz(j) 
        if ( sumx > nx_max ) nx_max = sumx 
        if ( sumy > ny_max ) ny_max = sumy 
        if ( sumz > nz_max ) nz_max = sumz 
        if ( sumx < nx_min ) nx_min = sumx 
        if ( sumy < ny_min ) ny_min = sumy 
        if ( sumz < nz_min ) nz_min = sumz 
     end do
  end do
  
  if ( iam == 0 ) write(6,*) nx_min, nx_max,  ny_min, ny_max,  nz_min, nz_max
  allocate( locate_t2channel(-1:1,  nx_min:nx_max, ny_min:ny_max, nz_min:nz_max) )
  allocate( locate_channel(1:6,-1:1, nx_min:nx_max, ny_min:ny_max, nz_min:nz_max) )
  allocate( numstates(1:3,-1:1, nx_min:nx_max, ny_min:ny_max, nz_min:nz_max) )
  locate_t2channel =  0 
  locate_channel = 0 
  numstates = 0 
  do a = 1,tot_orbs
     nxa = all_orbit%nx(a)
     nya = all_orbit%ny(a)
     nza = all_orbit%nz(a)
     tza=all_orbit%itzp(a)
     sza=all_orbit%szp(a)
     
     do b = 1, tot_orbs
        nxb = all_orbit%nx(b)
        nyb = all_orbit%ny(b)
        nzb = all_orbit%nz(b)
        tzb = all_orbit%itzp(b)
        szb = all_orbit%szp(b)
        
        tz = (tza + tzb)/2 
        sz = (sza + szb)/2 
        Nx = nxa + nxb 
        ny = nya + nyb 
        nz = nza + nzb
        
        if ( a <= below_ef .and. b <= below_ef ) numstates(1,tz, nx, ny, nz) = numstates(1,tz, nx, ny, nz)  + 1 
        if ( a <= below_ef .and. b >  below_ef ) numstates(2,tz, nx, ny, nz) = numstates(2,tz, nx, ny, nz)  + 1 
        if ( a > below_ef  .and. b >  below_ef ) numstates(3,tz, nx, ny, nz) = numstates(3,tz, nx, ny, nz)  + 1 
        
     end do
  end do
  
  number_channels = 0 
  numchannels = 0 
  !     loop over isospin projection
  DO tz=-1,1 
     !     loop over total momentum 
     DO nx=nx_min, nx_max 
        !     loop over total momentum 
        DO ny=ny_min, ny_max 
           !     loop over total momentum 
           DO nz=nz_min, nz_max 
              
              
              if ( numstates(1,tz, nx, ny, nz) /= 0 ) then
                 numchannels(1) = numchannels(1) + 1 
                 locate_channel(1,tz, nx, ny, nz) = numchannels(1)
                 !write(6,*) numchannels(1), numstates(1,tz, sz, nx, ny, nz)
              end if
              
              if ( numstates(1,tz, nx, ny, nz)* numstates(2,tz, nx, ny, nz) /= 0 ) then
                 numchannels(2) = numchannels(2) + 1 
                 locate_channel(2,tz, nx, ny, nz) = numchannels(2)
              end if
              
              if ( numstates(1,tz, nx, ny, nz)* numstates(3,tz, nx, ny, nz) /= 0 ) then
                 numchannels(3) = numchannels(3) + 1 
                 locate_channel(3,tz, nx, ny, nz) = numchannels(3)
                 !locate_channel(6,tz, nx, ny, nz) = numchannels(3)
              end if
              
              if ( numstates(2,tz, nx, ny, nz) /= 0 ) then 
                 numchannels(4) = numchannels(4) + 1 
                 locate_channel(4,tz, nx, ny, nz) = numchannels(4)
              end if
              
!!$              if ( numstates(2,tz, nx, ny, nz)* numstates(3,tz, nx, ny, nz) /= 0 ) then
!!$                 numchannels(5) = numchannels(5) + 1 
!!$                 locate_channel(5,tz, nx, ny, nz) = numchannels(5)
!!$              end if
              
              ! pppp has the same channels at t2 
              if ( numstates(1,tz, nx, ny, nz)*numstates(3,tz, nx, ny, nz)  /= 0 ) then
                 numchannels(6) = numchannels(6) + 1 
                 locate_channel(6,tz, nx, ny, nz) = numchannels(6)
              end if
              
           end DO
        end DO
     end DO
  end DO


  ALLOCATE( lookup_hhhh_configs(1:1,numchannels(1)))
  ALLOCATE( lookup_hhhp_configs(1:2,numchannels(2)))
  ALLOCATE( lookup_hhpp_configs(1:2,numchannels(3)))
  ALLOCATE( lookup_hphp_configs(1:1,numchannels(4)))
!!$  ALLOCATE( lookup_hppp_configs(1:2,numchannels(5)))
  if(cc_approx .ne. 'mbpt2') then
     ALLOCATE( lookup_pppp_configs(1:1,numchannels(6)))
  end if

  memory = 0.d0 
  !     loop over isospin projection
  DO tz=-1,1 
     !     loop over total momentum 
     DO nx=nx_min, nx_max 
        !     loop over total momentum 
        DO ny=ny_min, ny_max 
           !     loop over total momentum 
           DO nz=nz_min, nz_max 
                 
              channel1 = locate_channel(1, tz, nx, ny, nz) 
              channel2 = locate_channel(2, tz, nx, ny, nz) 
              channel3 = locate_channel(3, tz, nx, ny, nz) 
              channel4 = locate_channel(4, tz, nx, ny, nz) 
              channel5 = locate_channel(5, tz, nx, ny, nz) 
              channel6 = locate_channel(6, tz, nx, ny, nz) 
              if ( channel1 /= 0 ) then 
                 ij_confs = numstates(1,tz, nx, ny, nz) 
                 memory = memory + dble(ij_confs)*2.*4./1.e9 
                 
                 ALLOCATE( lookup_hhhh_configs(1,channel1)%ival(2,ij_confs) )
              end if
              if ( channel2 /= 0 ) then 
                 ij_confs = numstates(1,tz, nx, ny, nz) 
                 ab_confs = numstates(2,tz, nx, ny, nz) 
                 
                 memory = memory + dble(ij_confs+ab_confs)*2.*4./1.e9 
                 ALLOCATE( lookup_hhhp_configs(1,channel2)%ival(2,ij_confs) )
                 ALLOCATE( lookup_hhhp_configs(2,channel2)%ival(2,ab_confs) )
              end if
              if ( channel3 /= 0 ) then 
                 ij_confs = numstates(1,tz, nx, ny, nz) 
                 ab_confs = numstates(3,tz, nx, ny, nz) 
                 
                 memory = memory + dble(ij_confs+ab_confs)*2.*4./1.e9 
                 ALLOCATE( lookup_hhpp_configs(1,channel3)%ival(2,ij_confs) )
                 ALLOCATE( lookup_hhpp_configs(2,channel3)%ival(2,ab_confs) )
              end if
              if ( channel4 /= 0 ) then 
                 ij_confs = numstates(2,tz, nx, ny, nz) 
                 
                 memory = memory + dble(ij_confs)*2.*4./1.e9 
                 ALLOCATE( lookup_hphp_configs(1,channel4)%ival(2,ij_confs) )
              end if
!!$              if ( channel5 /= 0 ) then 
!!$                 ij_confs = numstates(2,tz, nx, ny, nz) 
!!$                 ab_confs = numstates(3,tz, nx, ny, nz) 
!!$                 
!!$                 memory = memory + dble(ij_confs+ab_confs)*2.*4./1.e9 
!!$                 ALLOCATE( lookup_hppp_configs(1,channel5)%ival(2,ij_confs) )
!!$                 ALLOCATE( lookup_hppp_configs(2,channel5)%ival(2,ab_confs) )
!!$              end if
              if ( channel6 /= 0 .and. channel3 /= 0 ) then 
                 ab_confs = numstates(3,tz, nx, ny, nz) 
                 if(cc_approx .ne. 'mbpt2')then
                    memory = memory + dble(ab_confs)*2.*4./1.e9 
                    ALLOCATE( lookup_pppp_configs(1,channel6)%ival(2,ab_confs) )
                 end if
              end if
                            
           end DO
        end DO
     end DO
  end DO
  
  if (iam == 0 ) write(6,*) 'Total memory for lookup configs', memory, 'GByte' 
  
  !
  ! setup hh configurations 
  !
  allocate( hh_hhhh%ival(below_ef,below_ef) )
  allocate( hh_hhhp%ival(below_ef,below_ef) )
  allocate( hh_hhpp%ival(below_ef,below_ef) )
  allocate( hp_hphp%ival(below_ef,below_ef+1:tot_orbs) )
!  allocate( hp_hppp%ival(below_ef,below_ef+1:tot_orbs) )
  allocate( hp_hhhp%ival(below_ef,below_ef+1:tot_orbs) )
  allocate( pp_hhpp%ival(below_ef+1:tot_orbs,below_ef+1:tot_orbs) )
!  allocate( pp_hppp%ival(below_ef+1:tot_orbs,below_ef+1:tot_orbs) )
  if(cc_approx .ne. 'mbpt2')then
     allocate( pp_pppp%ival(below_ef+1:tot_orbs,below_ef+1:tot_orbs) )
     pp_pppp%ival = 0   
  end if

  hh_hhhh%ival = 0; hh_hhhp%ival = 0; hh_hhpp%ival = 0; hp_hphp%ival = 0  
  hp_hhhp%ival = 0; pp_hhpp%ival = 0
 ! hp_hppp%ival = 0; pp_hppp%ival = 0  
 
  
  
  allocate( numchannels1( maxval( locate_channel(1,:,:,:,:) ) ) )
  allocate( numchannels2( maxval( locate_channel(2,:,:,:,:) ) ) )
  allocate( numchannels3( maxval( locate_channel(3,:,:,:,:) ) ) )
  allocate( numchannels4( maxval( locate_channel(4,:,:,:,:) ) ) )
 ! allocate( numchannels5( maxval( locate_channel(5,:,:,:,:) ) ) )
  allocate( numchannels6( maxval( locate_channel(6,:,:,:,:) ) ) )
  numchannels1 = 0
  numchannels2 = 0
  numchannels3 = 0
  numchannels4 = 0
 ! numchannels5 = 0
  numchannels6 = 0
  do a = 1,below_ef
     nxa = all_orbit%nx(a)
     nya = all_orbit%ny(a)
     nza = all_orbit%nz(a)
     tza=all_orbit%itzp(a)
     sza=all_orbit%szp(a)
     
     do b = 1, below_ef
        nxb = all_orbit%nx(b)
        nyb = all_orbit%ny(b)
        nzb = all_orbit%nz(b)
        tzb = all_orbit%itzp(b)
        szb = all_orbit%szp(b)
        
        tz = (tza + tzb)/2 
        sz = (sza + szb)/2 
        Nx = nxa + nxb 
        ny = nya + nyb 
        nz = nza + nzb
        

        if ( locate_channel(1, tz, nx, ny, nz) /= 0 ) then 
           
           channel = locate_channel(1, tz, nx, ny, nz) 
           numchannels1(channel) = numchannels1(channel) + 1 
           ichan = numchannels1(channel)
           
           hh_hhhh%ival(a,b) = ichan
           lookup_hhhh_configs(1,channel)%ival(1,ichan) = a 
           lookup_hhhh_configs(1,channel)%ival(2,ichan) = b 
        end if
        
        if ( locate_channel(2, tz, nx, ny, nz) /= 0 ) then
           channel = locate_channel(2, tz, nx, ny, nz) 
           numchannels2(channel) = numchannels2(channel) + 1 
           ichan = numchannels2(channel)
           
           hh_hhhp%ival(a,b) = ichan
           lookup_hhhp_configs(1,channel)%ival(1,ichan) = a 
           lookup_hhhp_configs(1,channel)%ival(2,ichan) = b 
           
        end if
        
        if ( locate_channel(3, tz, nx, ny, nz) /= 0 ) then
           channel = locate_channel(3, tz, nx, ny, nz) 
           numchannels3(channel) = numchannels3(channel) + 1 
           ichan = numchannels3(channel)
           
           hh_hhpp%ival(a,b) = ichan
           lookup_hhpp_configs(1,channel)%ival(1,ichan) = a 
           lookup_hhpp_configs(1,channel)%ival(2,ichan) = b 
       
        end if
        
     end do
  end do

  !
  ! setup hp configurations 
  !
  numchannels1 = 0
  numchannels2 = 0
  numchannels3 = 0
  numchannels4 = 0
!  numchannels5 = 0
  numchannels6 = 0
  do a = 1,below_ef
     nxa = all_orbit%nx(a)
     nya = all_orbit%ny(a)

     nza = all_orbit%nz(a)
     tza=all_orbit%itzp(a)
     sza=all_orbit%szp(a)
     
     do b = below_ef+1, tot_orbs
        nxb = all_orbit%nx(b)
        nyb = all_orbit%ny(b)
        nzb = all_orbit%nz(b)
        tzb = all_orbit%itzp(b)
        szb = all_orbit%szp(b)
        
        tz = (tza + tzb)/2 
        sz = (sza + szb)/2 
        Nx = nxa + nxb 
        ny = nya + nyb 
        nz = nza + nzb
        

        if ( locate_channel(2, tz, nx, ny, nz) /= 0 ) then
           channel = locate_channel(2, tz, nx, ny, nz) 
           numchannels2(channel) = numchannels2(channel) + 1 
           ichan = numchannels2(channel)

           hp_hhhp%ival(a,b) = ichan
           lookup_hhhp_configs(2,channel)%ival(1,ichan) = a 
           lookup_hhhp_configs(2,channel)%ival(2,ichan) = b 
           
        end if
        
        if ( locate_channel(4, tz, nx, ny, nz) /= 0 ) then
           channel = locate_channel(4, tz, nx, ny, nz) 
           numchannels4(channel) = numchannels4(channel) + 1 
           ichan = numchannels4(channel)
           
           hp_hphp%ival(a,b) = ichan
           lookup_hphp_configs(1,channel)%ival(1,ichan) = a 
           lookup_hphp_configs(1,channel)%ival(2,ichan) = b 
           
        end if
        
!!$        if ( locate_channel(5, tz, nx, ny, nz) /= 0 ) then
!!$           channel = locate_channel(5, tz, nx, ny, nz) 
!!$           numchannels5(channel) = numchannels5(channel) + 1 
!!$           ichan = numchannels5(channel)
!!$           
!!$           hp_hppp%ival(a,b) = ichan
!!$           lookup_hppp_configs(1,channel)%ival(1,ichan) = a 
!!$           lookup_hppp_configs(1,channel)%ival(2,ichan) = b 
!!$       
!!$        end if
        
     end do
  end do


  !
  ! setup pp configurations 
  !
  numchannels1 = 0
  numchannels2 = 0
  numchannels3 = 0
  numchannels4 = 0
!  numchannels5 = 0
  numchannels6 = 0
  do a = below_ef+1, tot_orbs
     nxa = all_orbit%nx(a)
     nya = all_orbit%ny(a)
     nza = all_orbit%nz(a)
     tza=all_orbit%itzp(a)
     sza=all_orbit%szp(a)
     
     do b = below_ef+1, tot_orbs
        nxb = all_orbit%nx(b)
        nyb = all_orbit%ny(b)
        nzb = all_orbit%nz(b)
        tzb = all_orbit%itzp(b)
        szb = all_orbit%szp(b)
        
        tz = (tza + tzb)/2 
        sz = (sza + szb)/2 
        Nx = nxa + nxb 
        ny = nya + nyb 
        nz = nza + nzb
        

        if ( locate_channel(3, tz, nx, ny, nz) /= 0 ) then
           channel = locate_channel(3, tz, nx, ny, nz) 
           numchannels3(channel) = numchannels3(channel) + 1 
           ichan = numchannels3(channel)

           pp_hhpp%ival(a,b) = ichan
           lookup_hhpp_configs(2,channel)%ival(1,ichan) = a 
           lookup_hhpp_configs(2,channel)%ival(2,ichan) = b 
           
        end if
        
!!$        if ( locate_channel(5, tz, nx, ny, nz) /= 0 ) then
!!$           channel = locate_channel(5, tz, nx, ny, nz) 
!!$           numchannels5(channel) = numchannels5(channel) + 1 
!!$           ichan = numchannels5(channel)
!!$
!!$           pp_hppp%ival(a,b) = ichan
!!$           lookup_hppp_configs(2,channel)%ival(1,ichan) = a 
!!$           lookup_hppp_configs(2,channel)%ival(2,ichan) = b 
!!$       
!!$        end if
        if ( locate_channel(6, tz, nx, ny, nz) /= 0 .and. & 
             locate_channel(3, tz, nx, ny, nz) /= 0 ) then
           channel = locate_channel(6, tz, nx, ny, nz) 
           numchannels6(channel) = numchannels6(channel) + 1 
           ichan = numchannels6(channel)
           if(cc_approx .ne. 'mbpt2') then
              pp_pppp%ival(a,b) = ichan
              lookup_pppp_configs(1,channel)%ival(1,ichan) = a 
              lookup_pppp_configs(1,channel)%ival(2,ichan) = b 
           end if
        end if
     end do
  end do
  

  memory = 0.d0 
  number_channels = 0 
  number_channels = maxval(locate_channel(1,:,:,:,:) ) 
  channels%number_hhhh_confs  = number_channels
  ALLOCATE( channels%hhhh_quantum_numbers(4*number_channels) )
  memory = memory + dble(number_channels)*4.*4./1.e9 
  
  number_channels = 0 
  number_channels = maxval(locate_channel(2,:,:,:,:) ) 
  channels%number_hhhp_confs  = number_channels
  ALLOCATE( channels%hhhp_quantum_numbers(4*number_channels) )
  memory = memory + dble(number_channels)*4.*4./1.e9 

  number_channels = 0 
  number_channels = maxval(locate_channel(3,:,:,:,:) ) 
  channels%number_hhpp_confs  = number_channels
  ALLOCATE( channels%hhpp_quantum_numbers(4*number_channels) )
  memory = memory + dble(number_channels)*4.*4./1.e9 

  number_channels = 0 
  number_channels = maxval(locate_channel(4,:,:,:,:) ) 
  channels%number_hphp_confs  = number_channels
  ALLOCATE( channels%hphp_quantum_numbers(4*number_channels) )
  memory = memory + dble(number_channels)*4.*4./1.e9 
  
!!$  number_channels = 0 
!!$  number_channels = maxval(locate_channel(5,:,:,:,:) ) 
!!$  channels%number_hppp_confs  = number_channels
!!$  ALLOCATE( channels%hppp_quantum_numbers(4*number_channels) )
!!$  memory = memory + dble(number_channels)*4.*4./1.e9 
  
  number_channels = 0 
  number_channels = maxval(locate_channel(6,:,:,:,:) ) 
  channels%number_pppp_confs  = number_channels
  ALLOCATE( channels%pppp_quantum_numbers(4*number_channels) )
  memory = memory + dble(number_channels)*4.*4./1.e9 
  
  if (iam == 0 ) write(6,*) 'Total memory for channels', memory, 'GByte' 
  
  
  !     loop over isospin projection
  DO tz=-1,1 
     !     loop over total momentum 
     DO nx=nx_min, nx_max 
        !     loop over total momentum 
        DO ny=ny_min, ny_max 
           !     loop over total momentum 
           DO nz=nz_min, nz_max 
              
              number_channels = locate_channel(1,tz, nx, ny, nz) 
              if ( number_channels == 0 ) cycle 
                    
              
              k1 = number_channels*4 
              k2 = number_channels*4 - 1
              k3 = number_channels*4 - 2
              k4 = number_channels*4 - 3
              channels%hhhh_quantum_numbers(k1) = nx
              channels%hhhh_quantum_numbers(k2) = ny
              channels%hhhh_quantum_numbers(k3) = nz
              channels%hhhh_quantum_numbers(k4) = tz
                               
              
           end DO
        end DO
     end DO
  end DO

  !     loop over isospin projection
  DO tz=-1,1 
     !     loop over total momentum 
     DO nx=nx_min, nx_max 
        !     loop over total momentum 
        DO ny=ny_min, ny_max 
           !     loop over total momentum 
           DO nz=nz_min, nz_max 
                 
              number_channels = locate_channel(2,tz, nx, ny, nz) 
              if ( number_channels == 0 ) cycle 
                    
                    
              k1 = number_channels*4 
              k2 = number_channels*4 - 1
              k3 = number_channels*4 - 2
              k4 = number_channels*4 - 3
              channels%hhhp_quantum_numbers(k1) = nx
              channels%hhhp_quantum_numbers(k2) = ny
              channels%hhhp_quantum_numbers(k3) = nz
              channels%hhhp_quantum_numbers(k4) = tz
              
           end DO
        end DO
     end DO
  end DO

  !     loop over isospin projection
  DO tz=-1,1 
     !     loop over total momentum 
     DO nx=nx_min, nx_max 
        !     loop over total momentum 
        DO ny=ny_min, ny_max 
           !     loop over total momentum 
           DO nz=nz_min, nz_max 
              
              number_channels = locate_channel(3,tz, nx, ny, nz) 
              if ( number_channels == 0 ) cycle 
                    
              
              k1 = number_channels*4 
              k2 = number_channels*4 - 1
              k3 = number_channels*4 - 2
              k4 = number_channels*4 - 3
              channels%hhpp_quantum_numbers(k1) = nx
              channels%hhpp_quantum_numbers(k2) = ny
              channels%hhpp_quantum_numbers(k3) = nz
              channels%hhpp_quantum_numbers(k4) = tz
                               
              
           end DO
        end DO
     end DO
  end DO

  !     loop over isospin projection
  DO tz=-1,1 
     !     loop over total momentum 
     DO nx=nx_min, nx_max 
        !     loop over total momentum 
        DO ny=ny_min, ny_max 
           !     loop over total momentum 
           DO nz=nz_min, nz_max 
                 
              number_channels = locate_channel(4,tz, nx, ny, nz) 
              if ( number_channels == 0 ) cycle 
              
                    
              k1 = number_channels*4 
              k2 = number_channels*4 - 1
              k3 = number_channels*4 - 2
              k4 = number_channels*4 - 3
              channels%hphp_quantum_numbers(k1) = nx
              channels%hphp_quantum_numbers(k2) = ny
              channels%hphp_quantum_numbers(k3) = nz
              channels%hphp_quantum_numbers(k4) = tz
              
           end DO
        end DO
     end DO
  end DO

!!$  !     loop over isospin projection
!!$  DO tz=-1,1 
!!$     !     loop over total momentum 
!!$     DO nx=nx_min, nx_max 
!!$        !     loop over total momentum 
!!$        DO ny=ny_min, ny_max 
!!$           !     loop over total momentum 
!!$           DO nz=nz_min, nz_max 
!!$              
!!$              number_channels = locate_channel(5,tz, nx, ny, nz) 
!!$              if ( number_channels == 0 ) cycle 
!!$              
!!$                    
!!$              k1 = number_channels*4 
!!$              k2 = number_channels*4 - 1
!!$              k3 = number_channels*4 - 2
!!$              k4 = number_channels*4 - 3
!!$              channels%hppp_quantum_numbers(k1) = nx
!!$              channels%hppp_quantum_numbers(k2) = ny
!!$              channels%hppp_quantum_numbers(k3) = nz
!!$              channels%hppp_quantum_numbers(k4) = tz
!!$              
!!$              
!!$              
!!$           end DO
!!$        end DO
!!$     end DO
!!$  end DO

     !     loop over isospin projection
     DO tz=-1,1 
        !     loop over total momentum 
        DO nx=nx_min, nx_max 
           !     loop over total momentum 
           DO ny=ny_min, ny_max 
              !     loop over total momentum 
              DO nz=nz_min, nz_max 
                 
                 
                 number_channels = locate_channel(6,tz, nx, ny, nz) 
                 if ( number_channels == 0 ) cycle 
                 
                 k1 = number_channels*4 
                 k2 = number_channels*4 - 1
                 k3 = number_channels*4 - 2
                 k4 = number_channels*4 - 3
                 channels%pppp_quantum_numbers(k1) = nx
                 channels%pppp_quantum_numbers(k2) = ny
                 channels%pppp_quantum_numbers(k3) = nz
                 channels%pppp_quantum_numbers(k4) = tz
                 
              end DO
           end DO
        end DO
     end DO
     
 

 call setup_proc_mappings_pp

  !
  ! setup chiral interactions 
  ! 
  allocate( vnn_hhhh(channels%number_hhhh_confs) )
  allocate( vnn_hhhp(channels%number_hhhp_confs) )
  allocate( vnn_hhpp(channels%number_hhpp_confs) )
  allocate( vnn_hphp(channels%number_hphp_confs) )
  !allocate( vnn_hppp(channels%number_hppp_confs) )
  !  allocate( vnn_pppp(channels%number_pppp_confs) )

  !
  ! setup hhhh NN interaction
  !

  startwtime = MPI_WTIME()
  if ( iam == 0 ) write(6,*) 'Setting up vnn for hhhh block'
  memory = 0.d0 
 
  if(chiral_delta_flag)then
     do channel   = 1, channels%number_hhhh_confs 
        nx = channels%hhhh_quantum_numbers(channel*4)
        ny = channels%hhhh_quantum_numbers(channel*4-1)
        nz = channels%hhhh_quantum_numbers(channel*4-2)
        tz = channels%hhhh_quantum_numbers(channel*4-3)
        
        dim1 = size(  lookup_hhhh_configs(1,channel)%ival, 2) 
        memory = memory + dble(dim1*dim1)*8./1.e9
        
        allocate( vnn_hhhh(channel)%val(dim1, dim1) )
        vnn_hhhh(channel)%val = 0.d0 
        do bra = 1, size(  lookup_hhhh_configs(1,channel)%ival, 2) 
           i= lookup_hhhh_configs(1,channel)%ival(1,bra) 
           j = lookup_hhhh_configs(1,channel)%ival(2,bra) 
           
           do ket = 1, size(  lookup_hhhh_configs(1,channel)%ival, 2) 
              k = lookup_hhhh_configs(1,channel)%ival(1,ket)
              l = lookup_hhhh_configs(1,channel)%ival(2,ket) 
              
              vnn_hhhh(channel)%val(bra,ket) = ( chiral_NN_with_delta(i,j,k,l) - chiral_NN_with_delta(i,j,l,k) )
              
           end do
        end do
     end do
     
  else
     
     do channel   = 1, channels%number_hhhh_confs 
        nx = channels%hhhh_quantum_numbers(channel*4)
        ny = channels%hhhh_quantum_numbers(channel*4-1)
        nz = channels%hhhh_quantum_numbers(channel*4-2)
        tz = channels%hhhh_quantum_numbers(channel*4-3)
        
        dim1 = size(  lookup_hhhh_configs(1,channel)%ival, 2) 
        memory = memory + dble(dim1*dim1)*8./1.e9
     
        allocate( vnn_hhhh(channel)%val(dim1, dim1) )
        vnn_hhhh(channel)%val = 0.d0 
        do bra = 1, size(  lookup_hhhh_configs(1,channel)%ival, 2) 
           i= lookup_hhhh_configs(1,channel)%ival(1,bra) 
           j = lookup_hhhh_configs(1,channel)%ival(2,bra) 
           
           do ket = 1, size(  lookup_hhhh_configs(1,channel)%ival, 2) 
              k = lookup_hhhh_configs(1,channel)%ival(1,ket)
              l = lookup_hhhh_configs(1,channel)%ival(2,ket) 
              
              vnn_hhhh(channel)%val(bra,ket) = ( chiral_pot(i,j,k,l) - chiral_pot(i,j,l,k) )
              
           end do
        end do
     end do
  end if


  endwtime = MPI_WTIME()
  if ( iam == 0 ) write(6,*) 'Time to setup hhhh block', endwtime - startwtime
  startwtime = MPI_WTIME()
  !
  ! setup hhhp NN interaction
  !
  if ( iam == 0 ) write(6,*) 'Setting up vnn for hhhp block'
  if(chiral_delta_flag) then
     do channel   = 1, channels%number_hhhp_confs 
        nx = channels%hhhp_quantum_numbers(channel*4)
        ny = channels%hhhp_quantum_numbers(channel*4-1)
        nz = channels%hhhp_quantum_numbers(channel*4-2)
        tz = channels%hhhp_quantum_numbers(channel*4-3)
        
        dim1 = size(  lookup_hhhp_configs(1,channel)%ival, 2) 
        dim2 = size(  lookup_hhhp_configs(2,channel)%ival, 2) 
        memory = memory + dble(dim1*dim2)*8./1.e9
        
        allocate( vnn_hhhp(channel)%val(dim1, dim2) )
        vnn_hhhp(channel)%val = 0.d0 
        do bra = 1, size(  lookup_hhhp_configs(1,channel)%ival, 2) 
           i= lookup_hhhp_configs(1,channel)%ival(1,bra) 
           j = lookup_hhhp_configs(1,channel)%ival(2,bra) 
           
           do ket = 1, size(  lookup_hhhp_configs(2,channel)%ival, 2) 
              k = lookup_hhhp_configs(2,channel)%ival(1,ket)
              d = lookup_hhhp_configs(2,channel)%ival(2,ket) 
              
              vnn_hhhp(channel)%val(bra,ket) = ( chiral_NN_with_delta(i,j,k,d) - chiral_NN_with_delta(i,j,d,k) )
              
              
           end do
        end do
     end do
  else  
     do channel   = 1, channels%number_hhhp_confs 
        nx = channels%hhhp_quantum_numbers(channel*4)
        ny = channels%hhhp_quantum_numbers(channel*4-1)
        nz = channels%hhhp_quantum_numbers(channel*4-2)
        tz = channels%hhhp_quantum_numbers(channel*4-3)
        
        dim1 = size(  lookup_hhhp_configs(1,channel)%ival, 2) 
        dim2 = size(  lookup_hhhp_configs(2,channel)%ival, 2) 
        memory = memory + dble(dim1*dim2)*8./1.e9
        
        allocate( vnn_hhhp(channel)%val(dim1, dim2) )
        vnn_hhhp(channel)%val = 0.d0 
        do bra = 1, size(  lookup_hhhp_configs(1,channel)%ival, 2) 
           i= lookup_hhhp_configs(1,channel)%ival(1,bra) 
           j = lookup_hhhp_configs(1,channel)%ival(2,bra) 
           
           do ket = 1, size(  lookup_hhhp_configs(2,channel)%ival, 2) 
              k = lookup_hhhp_configs(2,channel)%ival(1,ket)
              d = lookup_hhhp_configs(2,channel)%ival(2,ket) 
              
              vnn_hhhp(channel)%val(bra,ket) = ( chiral_pot(i,j,k,d) - chiral_pot(i,j,d,k) )
              
           end do
        end do
     end do
  end if


  endwtime = MPI_WTIME()
  if ( iam == 0 ) write(6,*) 'Time to setup hhhp block', endwtime - startwtime
  startwtime = MPI_WTIME()
  !
  ! setup hhpp NN interaction
  !
  if ( iam == 0 ) write(6,*) 'Setting up vnn for hhpp block'
 
  if(chiral_delta_flag) then
     do channel   = 1, channels%number_hhpp_confs 
        dim1 = size(  lookup_hhpp_configs(1,channel)%ival, 2) 
        dim2 = size(  lookup_hhpp_configs(2,channel)%ival, 2) 
        allocate( vnn_hhpp(channel)%val(dim1, dim2) )
        vnn_hhpp(channel)%val = 0.d0 
     end do
  
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
        
        do ket = bra_min, bra_max !1, size(  lookup_hhpp_configs(2,channel)%ival, 2) 
           !do ket = 1, size(  lookup_hhpp_configs(2,channel)%ival, 2) 
           
           c = lookup_hhpp_configs(2,channel)%ival(1,ket)
           d = lookup_hhpp_configs(2,channel)%ival(2,ket) 
           do bra = 1, size(  lookup_hhpp_configs(1,channel)%ival, 2) 
              i= lookup_hhpp_configs(1,channel)%ival(1,bra) 
              j = lookup_hhpp_configs(1,channel)%ival(2,bra) 
              
              vnn_hhpp(channel)%val(bra,ket) = ( chiral_NN_with_delta(i,j,c,d) -  chiral_NN_with_delta(i,j,d,c) )
              
           end do
        end do
     end do
  else
     do channel   = 1, channels%number_hhpp_confs 
        dim1 = size(  lookup_hhpp_configs(1,channel)%ival, 2) 
        dim2 = size(  lookup_hhpp_configs(2,channel)%ival, 2) 
        allocate( vnn_hhpp(channel)%val(dim1, dim2) )
        vnn_hhpp(channel)%val = 0.d0 
     end do
     
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
        
        do ket = bra_min, bra_max !1, size(  lookup_hhpp_configs(2,channel)%ival, 2) 
           !do ket = 1, size(  lookup_hhpp_configs(2,channel)%ival, 2) 
           
           c = lookup_hhpp_configs(2,channel)%ival(1,ket)
           d = lookup_hhpp_configs(2,channel)%ival(2,ket) 
           do bra = 1, size(  lookup_hhpp_configs(1,channel)%ival, 2) 
              i= lookup_hhpp_configs(1,channel)%ival(1,bra) 
              j = lookup_hhpp_configs(1,channel)%ival(2,bra) 
              
              vnn_hhpp(channel)%val(bra,ket) = ( chiral_pot(i,j,c,d) -  chiral_pot(i,j,d,c) )
              
           end do
        end do
     end do
  end if

  number_channels = size( vnn_hhpp )
  do channel = 1, number_channels 
     bra = size(vnn_hhpp(channel)%val,1)
     ket = size(vnn_hhpp(channel)%val,2)
     
     ndim = bra * ket
     allocate( t2_tmp1(ndim) )
     allocate( t2_tmp2(ndim) )
     t2_tmp1 = 0.d0
     t2_tmp2 = 0.d0

     i1 = 0 
     do i =1, bra
        do j = 1, ket
           i1 = i1 + 1
           t2_tmp1(i1) = vnn_hhpp(channel)%val(i,j)
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
           vnn_hhpp(channel)%val(i,j) = t2_tmp1(i1)
        end do
     end do
     
     deallocate( t2_tmp1, t2_tmp2 ) 
  end do

  endwtime = MPI_WTIME()
  if ( iam == 0 ) write(6,*) 'Time to setup hhpp block', endwtime - startwtime
  startwtime = MPI_WTIME()
  !
  ! setup hhpp NN interaction
  !
  if ( iam == 0 ) write(6,*) 'Setting up vnn for hphp block'

   do channel   = 1, channels%number_hphp_confs 

     dim1 = size(  lookup_hphp_configs(1,channel)%ival, 2) 
     memory = memory + dble(dim1*dim1)*8./1.e9
     allocate( vnn_hphp(channel)%val(dim1, dim1) )
     vnn_hphp(channel)%val = 0.d0 
  end do

  if(chiral_delta_flag) then
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
        
        do bra = bra_min, bra_max !1, size(  lookup_hhpp_configs(2,channel)%ival, 2) 
           !do bra = 1, size(  lookup_hphp_configs(1,channel)%ival, 2) 
           i= lookup_hphp_configs(1,channel)%ival(1,bra) 
           b = lookup_hphp_configs(1,channel)%ival(2,bra) 
           
           do ket = 1, size(  lookup_hphp_configs(1,channel)%ival, 2) 
              k = lookup_hphp_configs(1,channel)%ival(1,ket)
              d = lookup_hphp_configs(1,channel)%ival(2,ket) 
              
              vnn_hphp(channel)%val(bra,ket) = ( chiral_NN_with_delta(i,b,k,d) - chiral_NN_with_delta(i,b,d,k) )
           end do
        end do
     end do
  else 
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
        
        do bra = bra_min, bra_max !1, size(  lookup_hhpp_configs(2,channel)%ival, 2) 
           !do bra = 1, size(  lookup_hphp_configs(1,channel)%ival, 2) 
           i= lookup_hphp_configs(1,channel)%ival(1,bra) 
           b = lookup_hphp_configs(1,channel)%ival(2,bra) 
           
           do ket = 1, size(  lookup_hphp_configs(1,channel)%ival, 2) 
              k = lookup_hphp_configs(1,channel)%ival(1,ket)
              d = lookup_hphp_configs(1,channel)%ival(2,ket) 
              
              vnn_hphp(channel)%val(bra,ket) = ( chiral_pot(i,b,k,d) - chiral_pot(i,b,d,k) )
           end do
        end do
     end do
  end if 

  number_channels = size( vnn_hphp )
  do channel = 1, number_channels 
     bra = size(vnn_hphp(channel)%val,1)
     ket = size(vnn_hphp(channel)%val,2)
     
     ndim = bra * ket
     allocate( t2_tmp1(ndim) )
     allocate( t2_tmp2(ndim) )
     t2_tmp1 = 0.d0
     t2_tmp2 = 0.d0

     i1 = 0 
     do i =1, bra
        do j = 1, ket
           i1 = i1 + 1
           t2_tmp1(i1) = vnn_hphp(channel)%val(i,j)
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
           vnn_hphp(channel)%val(i,j) = t2_tmp1(i1)
        end do
     end do
     
     deallocate( t2_tmp1, t2_tmp2 ) 
  end do

  endwtime = MPI_WTIME()
  if ( iam == 0 ) write(6,*) 'Time to setup hphp block', endwtime - startwtime
  startwtime = MPI_WTIME()


  
  !
  ! setup pppp NN interaction
  !
 !!Changed 
  if(cc_approx .ne. 'mbpt2') then  
  if ( iam == 0 ) write(6,*) 'Setting up vnn for pppp block'

  if(chiral_delta_flag) then
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
        
        dim1 = size(  lookup_pppp_configs(1,channel)%ival, 2) 
        memory = memory + dble((bra_max-bra_min+1)*dim1)*16./1.e9
        
        allocate( vnn_pppp(channel2)%val(bra_min:bra_max, ket_min:ket_max) )
        vnn_pppp(channel2)%val = 0.d0 
        do bra = bra_min, bra_max !1, size(  lookup_pppp_configs(1,channel)%ival, 2) 
           a = lookup_pppp_configs(1,channel)%ival(1,bra) 
           b = lookup_pppp_configs(1,channel)%ival(2,bra) 
           
           do ket = ket_min, ket_max !1, size(  lookup_pppp_configs(1,channel)%ival, 2) 
              c = lookup_pppp_configs(1,channel)%ival(1,ket)
              d = lookup_pppp_configs(1,channel)%ival(2,ket) 
              ! XXXXXXXXXXXXXX
              vnn_pppp(channel2)%val(bra,ket) = ( chiral_NN_with_delta(a,b,c,d) - chiral_NN_with_delta(a,b,d,c) )
              
           end do
        end do
     end do
  else 
     
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
        
        dim1 = size(  lookup_pppp_configs(1,channel)%ival, 2) 
        memory = memory + dble((bra_max-bra_min+1)*dim1)*16./1.e9
        
        allocate( vnn_pppp(channel2)%val(bra_min:bra_max, ket_min:ket_max) )
        vnn_pppp(channel2)%val = 0.d0 
        do bra = bra_min, bra_max !1, size(  lookup_pppp_configs(1,channel)%ival, 2) 
           a = lookup_pppp_configs(1,channel)%ival(1,bra) 
           b = lookup_pppp_configs(1,channel)%ival(2,bra) 
           
           do ket = ket_min, ket_max !1, size(  lookup_pppp_configs(1,channel)%ival, 2) 
              c = lookup_pppp_configs(1,channel)%ival(1,ket)
              d = lookup_pppp_configs(1,channel)%ival(2,ket) 
              
              vnn_pppp(channel2)%val(bra,ket) = ( chiral_pot(a,b,c,d) - chiral_pot(a,b,d,c) )
              
              !if ( bra == ket ) vnn_pppp(channel2)%val(bra,ket) = 1.d0 
           
              
           end do
        end do
     end do
  end if
  
  endwtime = MPI_WTIME()
  if ( iam == 0 ) write(6,*) 'Time to setup pppp block', endwtime - startwtime
end IF
  if (iam == 0 ) write(6,*) 'Total memory for vnn storage', memory, 'GByte' 
!!$  
!!$  allocate( interm_pppp(channels%number_hhpp_confs) ) 
!!$  allocate( pp_pppp_rest%ival(below_ef+1:tot_orbs, below_ef+1:tot_orbs))
!!$  pp_pppp_rest%ival = 0
!!$  memory = 0.d0 
!!$  if ( iam == 0 ) write(6,*) 'Setting up intermediate for pppp block'
!!$  do channel   = 1, channels%number_pppp_confs 
!!$     nx = channels%pppp_quantum_numbers(channel*4)
!!$     ny = channels%pppp_quantum_numbers(channel*4-1)
!!$     nz = channels%pppp_quantum_numbers(channel*4-2)
!!$     tz = channels%pppp_quantum_numbers(channel*4-3)
!!$     
!!$     !
!!$     ! check that pppp channel is restricted to hhpp channels. 
!!$     !
!!$     channel2 = locate_channel(3,tz, nx, ny, nz)
!!$     if ( channel2 == 0 ) cycle 
!!$     if ( channel2 /= channel ) then
!!$        write(6,*) 'error', channel, channel2
!!$     end if
!!$          
!!$     dim1 = 0 
!!$     do bra = 1, size(  lookup_pppp_configs(1,channel)%ival, 2) 
!!$        a = lookup_pppp_configs(1,channel)%ival(1,bra) 
!!$        b = lookup_pppp_configs(1,channel)%ival(2,bra) 
!!$        
!!$!        if ( a >= b ) cycle 
!!$        dim1 = dim1 + 1 
!!$        pp_pppp_rest%ival(a,b) = dim1 
!!$     end do
!!$     
!!$     memory = memory + dble(dim1*dim1)*16./1.e9
!!$     allocate( interm_pppp(channel2)%val(dim1, dim1) ) 
!!$     interm_pppp(channel2)%val = 0.d0 
!!$     
!!$     
!!$  end do
!!$  if (iam == 0 ) write(6,*) 'Total memory for pppp intermediate storage', memory, 'GByte' 
!!$  
  
  deallocate( numstates, numchannels1, numchannels2, & 
       numchannels3, numchannels4, numchannels6 )
 
end SUBROUTINE setup_channel_structures

!
!
!
SUBROUTINE setup_vnn_hppp_block
  USE single_particle_orbits
  USE constants
  use one_body_operators
  use t2_storage
  use configurations
  USE KSPACE 
  use CHIRAL_POTENTIALS

  IMPLICIT NONE
  
  INTEGER :: i,j,k,l,a,b,c,d,dim1,dim2, number_channels, nx,ny,nz,sz,tz,  sumx, sumy, sumz  
  INTEGER :: k1,k2,k3,k4,k5, ab_confs, ij_confs, bra,ket, channel
  real*8 :: memory,  denom
  integer, allocatable :: numstates(:,:,:,:,:)
  INTEGER :: nxa, nya, nza, nxb, nyb, nzb,tza,tzb,sza,szb, ichan, p
  INTEGER :: numchannels(1:6), channel1, channel2, channel3, channel4, channel5, channel6
  INTEGER, ALLOCATABLE :: numchannels1(:), numchannels2(:), numchannels3(:), & 
       numchannels4(:), numchannels5(:), numchannels6(:)
  INTEGER :: ket_min, ket_max, bra_min, bra_max, i1, local_ch_id, ndim
  complex*16, allocatable :: t2_tmp1(:), t2_tmp2(:), amat(:,:)
  real*8  :: startwtime , endwtime
  

  nx_min = 1000
  nx_max = -1000
  ny_min = 1000
  ny_max = -1000
  nz_min = 1000
  nz_max = -1000
  do i = 1, tot_orbs
     do j = 1, tot_orbs
        
        sumx = all_orbit%nx(i) +  all_orbit%nx(j) 
        sumy = all_orbit%ny(i) +  all_orbit%ny(j) 
        sumz = all_orbit%nz(i) +  all_orbit%nz(j) 
        if ( sumx > nx_max ) nx_max = sumx 
        if ( sumy > ny_max ) ny_max = sumy 
        if ( sumz > nz_max ) nz_max = sumz 
        if ( sumx < nx_min ) nx_min = sumx 
        if ( sumy < ny_min ) ny_min = sumy 
        if ( sumz < nz_min ) nz_min = sumz 
     end do
  end do
  
  if ( iam == 0 ) write(6,*) nx_min, nx_max,  ny_min, ny_max,  nz_min, nz_max
  allocate( numstates(1:3,-1:1, nx_min:nx_max, ny_min:ny_max, nz_min:nz_max) )
  numstates = 0 
  do a = 1,tot_orbs
     nxa = all_orbit%nx(a)
     nya = all_orbit%ny(a)
     nza = all_orbit%nz(a)
     tza=all_orbit%itzp(a)
     sza=all_orbit%szp(a)
     
     do b = 1, tot_orbs
        nxb = all_orbit%nx(b)
        nyb = all_orbit%ny(b)
        nzb = all_orbit%nz(b)
        tzb = all_orbit%itzp(b)
        szb = all_orbit%szp(b)
        
        tz = (tza + tzb)/2 
        sz = (sza + szb)/2 
        Nx = nxa + nxb 
        ny = nya + nyb 
        nz = nza + nzb
        
        if ( a <= below_ef .and. b <= below_ef ) numstates(1,tz, nx, ny, nz) = numstates(1,tz, nx, ny, nz)  + 1 
        if ( a <= below_ef .and. b >  below_ef ) numstates(2,tz, nx, ny, nz) = numstates(2,tz, nx, ny, nz)  + 1 
        if ( a > below_ef  .and. b >  below_ef ) numstates(3,tz, nx, ny, nz) = numstates(3,tz, nx, ny, nz)  + 1 
        
     end do
  end do
  
  number_channels = 0 
  numchannels = 0 
  !     loop over isospin projection
  DO tz=-1,1 
     !     loop over total momentum 
     DO nx=nx_min, nx_max 
        !     loop over total momentum 
        DO ny=ny_min, ny_max 
           !     loop over total momentum 
           DO nz=nz_min, nz_max 
              
              
              if ( numstates(2,tz, nx, ny, nz)* numstates(3,tz, nx, ny, nz) /= 0 ) then
                 numchannels(5) = numchannels(5) + 1 
                 locate_channel(5,tz, nx, ny, nz) = numchannels(5)
              end if
              
           end DO
        end DO
     end DO
  end DO
  
  ALLOCATE( lookup_hppp_configs(1:2,numchannels(5)))
  
  memory = 0.d0 
  !     loop over isospin projection
  DO tz=-1,1 
     !     loop over total momentum 
     DO nx=nx_min, nx_max 
        !     loop over total momentum 
        DO ny=ny_min, ny_max 
           !     loop over total momentum 
           DO nz=nz_min, nz_max 
              
              channel5 = locate_channel(5, tz, nx, ny, nz) 
              if ( channel5 /= 0 ) then 
                 ij_confs = numstates(2,tz, nx, ny, nz) 
                 ab_confs = numstates(3,tz, nx, ny, nz) 
                 
                 memory = memory + dble(ij_confs+ab_confs)*2.*4./1.e9 
                 ALLOCATE( lookup_hppp_configs(1,channel5)%ival(2,ij_confs) )
                 ALLOCATE( lookup_hppp_configs(2,channel5)%ival(2,ab_confs) )
              end if
              
           end DO
        end DO
     end DO
  end DO
  
  if (iam == 0 ) write(6,*) 'Total memory for lookup configs', memory, 'GByte' 
  
  !
  ! setup hh configurations 
  !
  allocate( hp_hppp%ival(below_ef,below_ef+1:tot_orbs) )
  allocate( pp_hppp%ival(below_ef+1:tot_orbs,below_ef+1:tot_orbs) )
  hp_hppp%ival = 0; pp_hppp%ival = 0  
  allocate( numchannels5( maxval( locate_channel(5,:,:,:,:) ) ) )
  numchannels5 = 0
  
  !
  ! setup hp configurations 
  !
  numchannels5 = 0
  do a = 1,below_ef
     nxa = all_orbit%nx(a)
     nya = all_orbit%ny(a)

     nza = all_orbit%nz(a)
     tza=all_orbit%itzp(a)
     sza=all_orbit%szp(a)
     
     do b = below_ef+1, tot_orbs
        nxb = all_orbit%nx(b)
        nyb = all_orbit%ny(b)
        nzb = all_orbit%nz(b)
        tzb = all_orbit%itzp(b)
        szb = all_orbit%szp(b)
        
        tz = (tza + tzb)/2 
        sz = (sza + szb)/2 
        Nx = nxa + nxb 
        ny = nya + nyb 
        nz = nza + nzb
        
        if ( locate_channel(5, tz, nx, ny, nz) /= 0 ) then
           channel = locate_channel(5, tz, nx, ny, nz) 
           numchannels5(channel) = numchannels5(channel) + 1 
           ichan = numchannels5(channel)
           
           hp_hppp%ival(a,b) = ichan
           lookup_hppp_configs(1,channel)%ival(1,ichan) = a 
           lookup_hppp_configs(1,channel)%ival(2,ichan) = b 
       
        end if
        
     end do
  end do


  !
  ! setup pp configurations 
  !
  numchannels5 = 0
  do a = below_ef+1, tot_orbs
     nxa = all_orbit%nx(a)
     nya = all_orbit%ny(a)
     nza = all_orbit%nz(a)
     tza=all_orbit%itzp(a)
     sza=all_orbit%szp(a)
     
     do b = below_ef+1, tot_orbs
        nxb = all_orbit%nx(b)
        nyb = all_orbit%ny(b)
        nzb = all_orbit%nz(b)
        tzb = all_orbit%itzp(b)
        szb = all_orbit%szp(b)
        
        tz = (tza + tzb)/2 
        sz = (sza + szb)/2 
        Nx = nxa + nxb 
        ny = nya + nyb 
        nz = nza + nzb
        


        
        if ( locate_channel(5, tz, nx, ny, nz) /= 0 ) then
           channel = locate_channel(5, tz, nx, ny, nz) 
           numchannels5(channel) = numchannels5(channel) + 1 
           ichan = numchannels5(channel)

           pp_hppp%ival(a,b) = ichan
           lookup_hppp_configs(2,channel)%ival(1,ichan) = a 
           lookup_hppp_configs(2,channel)%ival(2,ichan) = b 
       
        end if
        
     end do
  end do
  

 
  number_channels = 0 
  number_channels = maxval(locate_channel(5,:,:,:,:) ) 
  channels%number_hppp_confs  = number_channels
  ALLOCATE( channels%hppp_quantum_numbers(4*number_channels) )
  memory = memory + dble(number_channels)*4.*4./1.e9 
  
  
  if (iam == 0 ) write(6,*) 'Total memory for channels', memory, 'GByte' 
  
  




  !     loop over isospin projection
  DO tz=-1,1 
     !     loop over total momentum 
     DO nx=nx_min, nx_max 
        !     loop over total momentum 
        DO ny=ny_min, ny_max 
           !     loop over total momentum 
           DO nz=nz_min, nz_max 
              
              number_channels = locate_channel(5,tz, nx, ny, nz) 
              if ( number_channels == 0 ) cycle 
              
                    
              k1 = number_channels*4 
              k2 = number_channels*4 - 1
              k3 = number_channels*4 - 2
              k4 = number_channels*4 - 3
              channels%hppp_quantum_numbers(k1) = nx
              channels%hppp_quantum_numbers(k2) = ny
              channels%hppp_quantum_numbers(k3) = nz
              channels%hppp_quantum_numbers(k4) = tz
              
              
              
           end DO
        end DO
     end DO
  end DO
  allocate( vnn_hppp(channels%number_hppp_confs) )
  
  call setup_proc_mappings_hppp
  
  startwtime = MPI_WTIME()

  !
  ! setup hppp NN interaction
  !
  
  do channel   = 1, channels%number_hppp_confs 
     
     dim1 = size(  lookup_hppp_configs(1,channel)%ival, 2) 
     dim2 = size(  lookup_hppp_configs(2,channel)%ival, 2) 
     memory = memory + dble(dim1*dim2)*8./1.e9
     
     allocate( vnn_hppp(channel)%val(dim1, dim2) )
     vnn_hppp(channel)%val = 0.d0 
  end do
  
  if ( tnf_approx >= 2 ) then 
     if ( iam == 0 ) write(6,*) 'Setting up vnn (NN+3NF) for hppp block' 
     do channel   = 1, channels%number_hppp_confs 
        nx = channels%hppp_quantum_numbers(channel*4)
        ny = channels%hppp_quantum_numbers(channel*4-1)
        nz = channels%hppp_quantum_numbers(channel*4-2)
        tz = channels%hppp_quantum_numbers(channel*4-3)
        
        if ( check_my_channel_hppp(channel) == 0 ) cycle
        local_ch_id = channel
        !
        ! ket side if fully stored on each proc
        !
        bra_min = mapping_hppp(iam+1,local_ch_id,2)
        bra_max = mapping_hppp(iam+1,local_ch_id,3)
        !
        ! bra side is distributed 
        !
        ket_min = mapping_hppp(iam+1,local_ch_id,4)
        ket_max = mapping_hppp(iam+1,local_ch_id,5)
        
        allocate( amat(size(vnn_hppp(channel)%val,1), size(vnn_hppp(channel)%val,2) ) ) 
        amat = 0.d0 
        if(chiral_delta_flag) then
           !$omp parallel default(shared) private(bra,i,j,ket,k,l,p)
           !$omp do schedule(dynamic), reduction(+:amat)
           do ket = bra_min, bra_max 
              
              i = lookup_hppp_configs(2,channel)%ival(1,ket)
              j = lookup_hppp_configs(2,channel)%ival(2,ket) 
              do bra = 1, size(  lookup_hppp_configs(1,channel)%ival, 2) 
                 k= lookup_hppp_configs(1,channel)%ival(1,bra) 
                 l = lookup_hppp_configs(1,channel)%ival(2,bra) 
                 
                 do p = 1, below_ef
                    
                    amat(bra,ket) = amat(bra,ket) + chiral_3nf_with_delta_asym(i,j,p,k,l,p) 
                    
                 end do
                 
                 amat(bra,ket) = amat(bra,ket) + ( chiral_NN_with_delta(i,j,k,l) - chiral_NN_with_delta(i,j,l,k) )           
                 
                 
              end do
           end do
           !$omp end do
           !$omp end parallel
        ELSE
           !$omp parallel default(shared) private(bra,i,j,ket,k,l,p)
           !$omp do schedule(dynamic), reduction(+:amat)
           do ket = bra_min, bra_max 
              
              i = lookup_hppp_configs(2,channel)%ival(1,ket)
              j = lookup_hppp_configs(2,channel)%ival(2,ket) 
              do bra = 1, size(  lookup_hppp_configs(1,channel)%ival, 2) 
                 k= lookup_hppp_configs(1,channel)%ival(1,bra) 
                 l = lookup_hppp_configs(1,channel)%ival(2,bra) 
                 
                 do p = 1, below_ef
                    amat(bra,ket) = amat(bra,ket) + chiral_3nf_asym(i,j,p,k,l,p) 
                    
                 end do
                 
                 amat(bra,ket) = amat(bra,ket) + ( chiral_pot(i,j,k,l) - chiral_pot(i,j,l,k) )
                 
              end do
           end do
           !$omp end do
           !$omp end parallel
        end if
        vnn_hppp(channel)%val = vnn_hppp(channel)%val + amat 
        deallocate( amat ) 
        
     end do
     
  else 
     if ( iam == 0 ) write(6,*) 'Setting up vnn (NN) for hppp block' 
     do channel   = 1, channels%number_hppp_confs 
        nx = channels%hppp_quantum_numbers(channel*4)
        ny = channels%hppp_quantum_numbers(channel*4-1)
        nz = channels%hppp_quantum_numbers(channel*4-2)
        tz = channels%hppp_quantum_numbers(channel*4-3)
        
        if ( check_my_channel_hppp(channel) == 0 ) cycle
        local_ch_id = channel
        !
        ! ket side if fully stored on each proc
        !
        bra_min = mapping_hppp(iam+1,local_ch_id,2)
        bra_max = mapping_hppp(iam+1,local_ch_id,3)
        !
        ! bra side is distributed 
        !
        ket_min = mapping_hppp(iam+1,local_ch_id,4)
        ket_max = mapping_hppp(iam+1,local_ch_id,5)
        
        allocate( amat(size(vnn_hppp(channel)%val,1), size(vnn_hppp(channel)%val,2) ) ) 
        amat = 0.d0 
        !$omp parallel default(shared) private(bra,i,j,ket,k,l,p)
        !$omp do schedule(dynamic), reduction(+:amat)
        do ket = bra_min, bra_max 
           
           i = lookup_hppp_configs(2,channel)%ival(1,ket)
           j = lookup_hppp_configs(2,channel)%ival(2,ket) 
           do bra = 1, size(  lookup_hppp_configs(1,channel)%ival, 2) 
              k= lookup_hppp_configs(1,channel)%ival(1,bra) 
              l = lookup_hppp_configs(1,channel)%ival(2,bra) 

              if(chiral_delta_flag) THEN
                 amat(bra,ket) = amat(bra,ket) + ( chiral_NN_with_delta(i,j,k,l) - chiral_NN_with_delta(i,j,l,k) )           
              else  
                 amat(bra,ket) = amat(bra,ket) + ( chiral_pot(i,j,k,l) - chiral_pot(i,j,l,k) )           
              end if

           end do
        end do
        !$omp end do
        !$omp end parallel
     
        vnn_hppp(channel)%val = vnn_hppp(channel)%val + amat 
        deallocate( amat ) 
          
     end do
  end if
  
  number_channels = size( vnn_hppp )
  do channel = 1, number_channels 
     bra = size(vnn_hppp(channel)%val,1)
     ket = size(vnn_hppp(channel)%val,2)
     
     ndim = bra * ket
     allocate( t2_tmp1(ndim) )
     allocate( t2_tmp2(ndim) )
     t2_tmp1 = 0.d0
     t2_tmp2 = 0.d0

     i1 = 0 
     do i =1, bra
        do j = 1, ket
           i1 = i1 + 1
           t2_tmp1(i1) = vnn_hppp(channel)%val(i,j)
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
           vnn_hppp(channel)%val(i,j) = t2_tmp1(i1)
        end do
     end do
     
     deallocate( t2_tmp1, t2_tmp2 ) 
  end do
  
  
  endwtime = MPI_WTIME()
  if ( iam == 0 ) write(6,*) 'Time to setup hppp block', endwtime - startwtime
  if (iam == 0 ) write(6,*) 'Total memory for vnn storage', memory, 'GByte' 
  
  deallocate( numstates, numchannels5 )
  
end SUBROUTINE setup_vnn_hppp_block


!
!
!
SUBROUTINE setup_ph_channel_structures
  USE single_particle_orbits
  USE constants
  use one_body_operators
  use t2_storage
  use configurations
  USE KSPACE 
  use CHIRAL_POTENTIALS

  IMPLICIT NONE
  
  INTEGER :: i,j,k,l,a,b,c,d,dim1,dim2, number_channels, nx,ny,nz,sz,tz,  sumx, sumy, sumz  
  INTEGER :: k1,k2,k3,k4,k5, ab_confs, ij_confs, bra,ket, channel
  real*8 :: memory,  denom
  integer, allocatable :: numstates(:,:,:,:,:)
  INTEGER :: nxa, nya, nza, nxb, nyb, nzb,tza,tzb,sza,szb, ichan   
  INTEGER :: numchannels(1:3), channel1, channel2, channel3, channel4, channel5, channel6
  INTEGER, ALLOCATABLE :: numchannels1(:), numchannels2(:), numchannels3(:)
  INTEGER :: ket_min, ket_max, bra_min, bra_max, i1, local_ch_id, ndim
  INTEGER :: nx2, ny2, nz2, tz2, bra2, ket2
  complex*16, allocatable :: t2_tmp1(:), t2_tmp2(:) 
  real*8  :: startwtime , endwtime
  

  nx_min = 1000
  nx_max = -1000
  ny_min = 1000
  ny_max = -1000
  nz_min = 1000
  nz_max = -1000
  do i = 1, tot_orbs
     do j = 1, tot_orbs
        
        sumx = all_orbit%nx(i) -  all_orbit%nx(j) 
        sumy = all_orbit%ny(i) -  all_orbit%ny(j) 
        sumz = all_orbit%nz(i) -  all_orbit%nz(j) 
        if ( sumx > nx_max ) nx_max = sumx 
        if ( sumy > ny_max ) ny_max = sumy 
        if ( sumz > nz_max ) nz_max = sumz 
        if ( sumx < nx_min ) nx_min = sumx 
        if ( sumy < ny_min ) ny_min = sumy 
        if ( sumz < nz_min ) nz_min = sumz 
     end do
  end do
  
  if ( iam == 0 ) write(6,*) nx_min, nx_max,  ny_min, ny_max,  nz_min, nz_max
  allocate( locate_ph_channel(1:3,-1:1, nx_min:nx_max, ny_min:ny_max, nz_min:nz_max) )
  allocate( numstates(1:2,-1:1, nx_min:nx_max, ny_min:ny_max, nz_min:nz_max) )
  locate_ph_channel = 0 
  numstates = 0 
  do a = 1,tot_orbs
     nxa = all_orbit%nx(a)
     nya = all_orbit%ny(a)
     nza = all_orbit%nz(a)
     tza=all_orbit%itzp(a)
     sza=all_orbit%szp(a)
     
     do b = 1, tot_orbs
        nxb = all_orbit%nx(b)
        nyb = all_orbit%ny(b)
        nzb = all_orbit%nz(b)
        tzb = all_orbit%itzp(b)
        szb = all_orbit%szp(b)
        
        tz = (tza - tzb)/2 
        Nx = nxa - nxb 
        ny = nya - nyb 
        nz = nza - nzb
        
        if ( a <= below_ef .and. b >  below_ef ) numstates(1,tz, nx, ny, nz) = numstates(1,tz, nx, ny, nz)  + 1 
        if ( a >  below_ef .and. b <= below_ef ) numstates(2,tz, nx, ny, nz) = numstates(2,tz, nx, ny, nz)  + 1 
        
     end do
  end do
  
  number_channels = 0 
  numchannels = 0 
  !     loop over isospin projection
  DO tz=-1,1 
     !     loop over total momentum 
     DO nx=nx_min, nx_max 
        !     loop over total momentum 
        DO ny=ny_min, ny_max 
           !     loop over total momentum 
           DO nz=nz_min, nz_max 
              
              
              ! phhp channel 
              if ( numstates(2,tz, nx, ny, nz)* numstates(1,tz, nx, ny, nz) /= 0 ) then
                 numchannels(1) = numchannels(1) + 1 
                 locate_ph_channel(1,tz, nx, ny, nz) = numchannels(1)
              end if
              
              ! hphp channel 
              if ( numstates(1,tz, nx, ny, nz) /= 0 ) then
                 numchannels(2) = numchannels(2) + 1 
                 locate_ph_channel(2,tz, nx, ny, nz) = numchannels(2)
              end if
              
              ! hpph channel 
              if ( numstates(1,tz, nx, ny, nz)* numstates(2,tz, nx, ny, nz) /= 0 ) then
                 numchannels(3) = numchannels(3) + 1 
                 locate_ph_channel(3,tz, nx, ny, nz) = numchannels(3)
              end if
              
              
           end DO
        end DO
     end DO
  end DO


  ALLOCATE( lookup_ph_phhp_configs(1:2,numchannels(1)))
  ALLOCATE( lookup_ph_hphp_configs(1:1,numchannels(2)))
  ALLOCATE( lookup_ph_hpph_configs(1:2,numchannels(3)))
  
  memory = 0.d0 
  !     loop over isospin projection
  DO tz=-1,1 
     !     loop over total momentum 
     DO nx=nx_min, nx_max 
        !     loop over total momentum 
        DO ny=ny_min, ny_max 
           !     loop over total momentum 
           DO nz=nz_min, nz_max 
                 
              channel1 = locate_ph_channel(1, tz, nx, ny, nz) 
              channel2 = locate_ph_channel(2, tz, nx, ny, nz) 
              channel3 = locate_ph_channel(3, tz, nx, ny, nz) 
              if ( channel1 /= 0 ) then 
                 ij_confs = numstates(2,tz, nx, ny, nz) 
                 ab_confs = numstates(1,tz, nx, ny, nz) 
                 
                 memory = memory + dble(ij_confs+ab_confs)*2.*4./1.e9 
                 ALLOCATE( lookup_ph_phhp_configs(1,channel1)%ival(2,ij_confs) )
                 ALLOCATE( lookup_ph_phhp_configs(2,channel1)%ival(2,ab_confs) )
              end if
              if ( channel2 /= 0 ) then 
                 ij_confs = numstates(1,tz, nx, ny, nz) 
                 
                 memory = memory + dble(ij_confs)*2.*4./1.e9 
                 ALLOCATE( lookup_ph_hphp_configs(1,channel2)%ival(2,ij_confs) )
              end if
              if ( channel3 /= 0 ) then 
                 ij_confs = numstates(1,tz, nx, ny, nz) 
                 ab_confs = numstates(2,tz, nx, ny, nz) 
                 
                 memory = memory + dble(ij_confs+ab_confs)*2.*4./1.e9 
                 ALLOCATE( lookup_ph_hpph_configs(1,channel3)%ival(2,ij_confs) )
                 ALLOCATE( lookup_ph_hpph_configs(2,channel3)%ival(2,ab_confs) )
              end if
              
           end DO
        end DO
     end DO
  end DO
  
  if (iam == 0 ) write(6,*) 'Total memory for lookup ph configs', memory, 'GByte' 
  
  !
  ! setup hh configurations 
  !
  allocate( ph_ph_phhp%ival(below_ef+1:tot_orbs,below_ef))
  allocate( hp_ph_phhp%ival(below_ef,below_ef+1:tot_orbs))
  allocate( hp_ph_hphp%ival(below_ef,below_ef+1:tot_orbs))
  allocate( hp_ph_hpph%ival(below_ef,below_ef+1:tot_orbs))
  allocate( ph_ph_hpph%ival(below_ef+1:tot_orbs,below_ef))
  
  ph_ph_phhp%ival = 0; hp_ph_phhp%ival = 0
  hp_ph_hphp%ival = 0
  hp_ph_hpph%ival = 0; ph_ph_hpph%ival = 0
  
  
  !
  ! setup hp configurations 
  !
  
  allocate( numchannels1( maxval( locate_ph_channel(1,:,:,:,:) ) ) )
  allocate( numchannels2( maxval( locate_ph_channel(2,:,:,:,:) ) ) )
  allocate( numchannels3( maxval( locate_ph_channel(3,:,:,:,:) ) ) )
  numchannels1 = 0
  numchannels2 = 0
  numchannels3 = 0
  do a = 1, below_ef
     nxa = all_orbit%nx(a)
     nya = all_orbit%ny(a)
     nza = all_orbit%nz(a)
     tza=all_orbit%itzp(a)
     sza=all_orbit%szp(a)
     
     do b = below_ef+1, tot_orbs
        nxb = all_orbit%nx(b)
        nyb = all_orbit%ny(b)
        nzb = all_orbit%nz(b)
        tzb = all_orbit%itzp(b)
        szb = all_orbit%szp(b)
        
        tz = (tza - tzb)/2 
        Nx = nxa - nxb 
        ny = nya - nyb 
        nz = nza - nzb
        
        
        if ( locate_ph_channel(1, tz, nx, ny, nz) /= 0 ) then 
           
           channel = locate_ph_channel(1, tz, nx, ny, nz) 
           numchannels1(channel) = numchannels1(channel) + 1 
           ichan = numchannels1(channel)
           
           hp_ph_phhp%ival(a,b) = ichan
           lookup_ph_phhp_configs(2,channel)%ival(1,ichan) = a 
           lookup_ph_phhp_configs(2,channel)%ival(2,ichan) = b 
        end if
        
        if ( locate_ph_channel(2, tz, nx, ny, nz) /= 0 ) then
           channel = locate_ph_channel(2, tz, nx, ny, nz) 
           numchannels2(channel) = numchannels2(channel) + 1 
           ichan = numchannels2(channel)
           
           hp_ph_hphp%ival(a,b) = ichan
           lookup_ph_hphp_configs(1,channel)%ival(1,ichan) = a 
           lookup_ph_hphp_configs(1,channel)%ival(2,ichan) = b 
           
        end if
        
        if ( locate_ph_channel(3, tz, nx, ny, nz) /= 0 ) then
           channel = locate_ph_channel(3, tz, nx, ny, nz) 
           numchannels3(channel) = numchannels3(channel) + 1 
           ichan = numchannels3(channel)
           
           hp_ph_hpph%ival(a,b) = ichan
           lookup_ph_hpph_configs(1,channel)%ival(1,ichan) = a 
           lookup_ph_hpph_configs(1,channel)%ival(2,ichan) = b 
           
        end if
        
     end do
  end do

  !
  ! setup ph configurations 
  !
  numchannels1 = 0
  numchannels2 = 0
  numchannels3 = 0
  do a = below_ef+1, tot_orbs
     nxa = all_orbit%nx(a)
     nya = all_orbit%ny(a)
     nza = all_orbit%nz(a)
     tza=all_orbit%itzp(a)
     sza=all_orbit%szp(a)
     
     do b = 1, below_ef
        nxb = all_orbit%nx(b)
        nyb = all_orbit%ny(b)
        nzb = all_orbit%nz(b)
        tzb = all_orbit%itzp(b)
        szb = all_orbit%szp(b)
        
        tz = (tza - tzb)/2 
        Nx = nxa - nxb 
        ny = nya - nyb 
        nz = nza - nzb
        
        
        if ( locate_ph_channel(1, tz, nx, ny, nz) /= 0 ) then 
           
           channel = locate_ph_channel(1, tz, nx, ny, nz) 
           numchannels1(channel) = numchannels1(channel) + 1 
           ichan = numchannels1(channel)
           
           ph_ph_phhp%ival(a,b) = ichan
           lookup_ph_phhp_configs(1,channel)%ival(1,ichan) = a 
           lookup_ph_phhp_configs(1,channel)%ival(2,ichan) = b 
        end if
        
        
        if ( locate_ph_channel(3, tz, nx, ny, nz) /= 0 ) then
           channel = locate_ph_channel(3, tz, nx, ny, nz) 
           numchannels3(channel) = numchannels3(channel) + 1 
           ichan = numchannels3(channel)
           
           ph_ph_hpph%ival(a,b) = ichan
           lookup_ph_hpph_configs(2,channel)%ival(1,ichan) = a 
           lookup_ph_hpph_configs(2,channel)%ival(2,ichan) = b 
           
        end if
        
     end do
  end do
  
  
  memory = 0.d0 
  number_channels = 0 
  number_channels = maxval(locate_ph_channel(1,:,:,:,:) ) 
  channels%number_ph_phhp_confs  = number_channels
  ALLOCATE( channels%ph_phhp_quantum_numbers(4*number_channels) )
  memory = memory + dble(number_channels)*4.*4./1.e9 
  
  number_channels = 0 
  number_channels = maxval(locate_ph_channel(2,:,:,:,:) ) 
  channels%number_ph_hphp_confs  = number_channels
  ALLOCATE( channels%ph_hphp_quantum_numbers(4*number_channels) )
  memory = memory + dble(number_channels)*4.*4./1.e9 

  number_channels = 0 
  number_channels = maxval(locate_ph_channel(3,:,:,:,:) ) 
  channels%number_ph_hpph_confs  = number_channels
  ALLOCATE( channels%ph_hpph_quantum_numbers(4*number_channels) )
  memory = memory + dble(number_channels)*4.*4./1.e9 
  
  if (iam == 0 ) write(6,*) 'Total memory for ph channels', memory, 'GByte' 
    
  !     loop over isospin projection
  DO tz=-1,1 
     !     loop over total momentum 
     DO nx=nx_min, nx_max 
        !     loop over total momentum 
        DO ny=ny_min, ny_max 
           !     loop over total momentum 
           DO nz=nz_min, nz_max 
              
              number_channels = locate_ph_channel(1,tz, nx, ny, nz) 
              if ( number_channels == 0 ) cycle 
                    
              
              k1 = number_channels*4 
              k2 = number_channels*4 - 1
              k3 = number_channels*4 - 2
              k4 = number_channels*4 - 3
              channels%ph_phhp_quantum_numbers(k1) = nx
              channels%ph_phhp_quantum_numbers(k2) = ny
              channels%ph_phhp_quantum_numbers(k3) = nz
              channels%ph_phhp_quantum_numbers(k4) = tz
                               
              
           end DO
        end DO
     end DO
  end DO

  !     loop over isospin projection
  DO tz=-1,1 
     !     loop over total momentum 
     DO nx=nx_min, nx_max 
        !     loop over total momentum 
        DO ny=ny_min, ny_max 
           !     loop over total momentum 
           DO nz=nz_min, nz_max 
                 
              number_channels = locate_ph_channel(2,tz, nx, ny, nz) 
              if ( number_channels == 0 ) cycle 
              
                    
              k1 = number_channels*4 
              k2 = number_channels*4 - 1
              k3 = number_channels*4 - 2
              k4 = number_channels*4 - 3
              channels%ph_hphp_quantum_numbers(k1) = nx
              channels%ph_hphp_quantum_numbers(k2) = ny
              channels%ph_hphp_quantum_numbers(k3) = nz
              channels%ph_hphp_quantum_numbers(k4) = tz
              
           end DO
        end DO
     end DO
  end DO

  !     loop over isospin projection
  DO tz=-1,1 
     !     loop over total momentum 
     DO nx=nx_min, nx_max 
        !     loop over total momentum 
        DO ny=ny_min, ny_max 
           !     loop over total momentum 
           DO nz=nz_min, nz_max 
              
              number_channels = locate_ph_channel(3,tz, nx, ny, nz) 
              if ( number_channels == 0 ) cycle 
                    
              
              k1 = number_channels*4 
              k2 = number_channels*4 - 1
              k3 = number_channels*4 - 2
              k4 = number_channels*4 - 3
              channels%ph_hpph_quantum_numbers(k1) = nx
              channels%ph_hpph_quantum_numbers(k2) = ny
              channels%ph_hpph_quantum_numbers(k3) = nz
              channels%ph_hpph_quantum_numbers(k4) = tz
              
              
           end DO
        end DO
     end DO
  end DO
  
  call setup_proc_mappings_ph
  
  !
  ! setup chiral interactions 
  ! 
!!$  allocate( vnn_ph_hpph(channels%number_ph_hpph_confs) )
  
  !
  ! setup hhhh NN interaction
  !

  startwtime = MPI_WTIME()
  if ( iam == 0 ) write(6,*) 'Setting up vnn for ph hpph block'
  memory = 0.d0 
!!$  do channel   = 1, channels%number_ph_hpph_confs 
!!$     nx = channels%ph_hpph_quantum_numbers(channel*4)
!!$     ny = channels%ph_hpph_quantum_numbers(channel*4-1)
!!$     nz = channels%ph_hpph_quantum_numbers(channel*4-2)
!!$     tz = channels%ph_hpph_quantum_numbers(channel*4-3)
!!$     
!!$     dim1 = size(  lookup_ph_hpph_configs(1,channel)%ival, 2) 
!!$     dim2 = size(  lookup_ph_hpph_configs(2,channel)%ival, 2) 
!!$     memory = memory + dble(dim1*dim2)*16./1.e9
!!$     
!!$     allocate( vnn_ph_hpph(channel)%val(dim1, dim2) )
!!$     vnn_ph_hpph(channel)%val = 0.d0 
!!$     do bra = 1, size(  lookup_ph_hpph_configs(1,channel)%ival, 2) 
!!$        i = lookup_ph_hpph_configs(1,channel)%ival(1,bra) 
!!$        a = lookup_ph_hpph_configs(1,channel)%ival(2,bra) 
!!$        
!!$        do ket = 1, size(  lookup_ph_hpph_configs(2,channel)%ival, 2) 
!!$           b = lookup_ph_hpph_configs(2,channel)%ival(1,ket)
!!$           j = lookup_ph_hpph_configs(2,channel)%ival(2,ket) 
!!$           
!!$           nx2 = all_orbit%nx(i) + all_orbit%nx(j)
!!$           ny2 = all_orbit%ny(i) + all_orbit%ny(j)
!!$           nz2 = all_orbit%nz(i) + all_orbit%nz(j)
!!$           tz2 = (all_orbit%itzp(i) + all_orbit%itzp(j))/2
!!$           channel2 = locate_channel(3,tz2, nx2, ny2, nz2) 
!!$           if ( channel2 == 0 ) cycle 
!!$           
!!$           bra2 = hh_hhpp%ival(i,j)
!!$           ket2 = pp_hhpp%ival(a,b)
!!$           if ( bra2*ket2 == 0 ) cycle 
!!$           
!!$           vnn_ph_hpph(channel)%val(bra,ket) = vnn_hhpp(channel2)%val(bra2,ket2)
!!$           
!!$        end do
!!$     end do
!!$  end do
  
  allocate( t2_ph_ccm(channels%number_ph_phhp_confs) )
  allocate( t2_ph_ccm_eqn(channels%number_ph_phhp_confs) )
  do channel   = 1, channels%number_ph_phhp_confs 
     nx = channels%ph_phhp_quantum_numbers(channel*4)
     ny = channels%ph_phhp_quantum_numbers(channel*4-1)
     nz = channels%ph_phhp_quantum_numbers(channel*4-2)
     tz = channels%ph_phhp_quantum_numbers(channel*4-3)
     
     dim1 = size(  lookup_ph_phhp_configs(1,channel)%ival, 2) 
     dim2 = size(  lookup_ph_phhp_configs(2,channel)%ival, 2) 
     memory = memory + dble(dim1*dim2)*16.*2./1.e9
     
     allocate( t2_ph_ccm(channel)%val(dim1, dim2) )
     allocate( t2_ph_ccm_eqn(channel)%val(dim1, dim2) )
     t2_ph_ccm(channel)%val = 0.d0 
     t2_ph_ccm_eqn(channel)%val = 0.d0 
     
  end do
  
!!$  number_channels = size( vnn_hhpp )
!!$  do channel = 1, number_channels 
!!$     bra = size(vnn_hhpp(channel)%val,1)
!!$     ket = size(vnn_hhpp(channel)%val,2)
!!$     
!!$     ndim = bra * ket
!!$     allocate( t2_tmp1(ndim) )
!!$     allocate( t2_tmp2(ndim) )
!!$     t2_tmp1 = 0.d0
!!$     t2_tmp2 = 0.d0
!!$
!!$     i1 = 0 
!!$     do i =1, bra
!!$        do j = 1, ket
!!$           i1 = i1 + 1
!!$           t2_tmp1(i1) = vnn_hhpp(channel)%val(i,j)
!!$        end do
!!$     end do
!!$     
!!$     call mpi_reduce(t2_tmp1,t2_tmp2,ndim,mpi_complex16,mpi_sum, &
!!$          master,mpi_comm_world,ierror)
!!$     t2_tmp1 = 0.d0 
!!$     t2_tmp1 = t2_tmp2 
!!$     call mpi_bcast(t2_tmp1,ndim,mpi_complex16,master, &
!!$          mpi_comm_world,ierror)
!!$     
!!$     i1 = 0 
!!$     do i =1, bra
!!$        do j = 1, ket
!!$           i1 = i1 + 1
!!$           vnn_hhpp(channel)%val(i,j) = t2_tmp1(i1)
!!$        end do
!!$     end do
!!$     
!!$     deallocate( t2_tmp1, t2_tmp2 ) 
!!$  end do
  endwtime = MPI_WTIME()
  if ( iam == 0 ) write(6,*) 'Time to setup ph hpph block', endwtime - startwtime
  
  if (iam == 0 ) write(6,*) 'Total memory for vnn/t2 ph storage', memory, 'GByte' 
  
  deallocate( numstates, numchannels1, numchannels2, numchannels3 )
  
end SUBROUTINE setup_ph_channel_structures



SUBROUTINE setup_t_amplitudes 
  USE single_particle_orbits
  USE constants
  use one_body_operators
  use t2_storage
  use configurations
  USE KSPACE 
  use CHIRAL_POTENTIALS

  IMPLICIT NONE 
  REAL*8 :: memory, denom
  INTEGER :: channel,bra,ket, dim1, dim2, a,b,i,j, nx2, ny2, nz2, tz2, bra2, ket2, channel2

  memory = 0.d0 
  allocate( t2_ccm(channels%number_hhpp_confs ) )
  allocate( t2_ccm_eqn(channels%number_hhpp_confs ) )
  do channel   = 1, channels%number_hhpp_confs 
     
     dim1 = size(  lookup_hhpp_configs(1,channel)%ival, 2) 
     dim2 = size(  lookup_hhpp_configs(2,channel)%ival, 2) 
     memory = memory + dble(dim1*dim2)*16./1.e9
     
     allocate( t2_ccm(channel)%val(dim2, dim1) )
     allocate( t2_ccm_eqn(channel)%val(dim2, dim1) )
     t2_ccm(channel)%val = 0.d0 
     t2_ccm_eqn(channel)%val = 0.d0 
     do bra = 1, size(  lookup_hhpp_configs(1,channel)%ival, 2) 
        i= lookup_hhpp_configs(1,channel)%ival(1,bra) 
        j = lookup_hhpp_configs(1,channel)%ival(2,bra) 
        
        do ket = 1, size(  lookup_hhpp_configs(2,channel)%ival, 2) 
           a = lookup_hhpp_configs(2,channel)%ival(1,ket)
           b = lookup_hhpp_configs(2,channel)%ival(2,ket) 
           
           denom = -fock_mtx(a,a) - fock_mtx(b,b) + fock_mtx(i,i) + fock_mtx(j,j) 
           if ( abs( denom ) == 0.d0 ) cycle 
           
           t2_ccm(channel)%val(ket,bra) = conjg(vnn_hhpp(channel)%val(bra,ket)) / denom 
        end do
     end do
  end do
  
  
  do channel   = 1, channels%number_ph_phhp_confs 
     
     t2_ph_ccm(channel)%val = 0.d0 
     t2_ph_ccm_eqn(channel)%val = 0.d0 

     do bra = 1, size(  lookup_ph_phhp_configs(1,channel)%ival, 2) 
        a = lookup_ph_phhp_configs(1,channel)%ival(1,bra) 
        i = lookup_ph_phhp_configs(1,channel)%ival(2,bra) 
        
        do ket = 1, size(  lookup_ph_phhp_configs(2,channel)%ival, 2) 
           j = lookup_ph_phhp_configs(2,channel)%ival(1,ket)
           b = lookup_ph_phhp_configs(2,channel)%ival(2,ket) 
           
           nx2 = all_orbit%nx(i) + all_orbit%nx(j)
           ny2 = all_orbit%ny(i) + all_orbit%ny(j)
           nz2 = all_orbit%nz(i) + all_orbit%nz(j)
           tz2 = (all_orbit%itzp(i) + all_orbit%itzp(j))/2
           channel2 = locate_channel(3,tz2, nx2, ny2, nz2) 
           if ( channel2 == 0 ) cycle 
           
           bra2 = hh_hhpp%ival(i,j)
           ket2 = pp_hhpp%ival(a,b)
           if ( bra2*ket2 == 0 ) cycle 
           
           t2_ph_ccm(channel)%val(bra,ket) = t2_ccm(channel2)%val(ket2,bra2)
           
        end do
     end do
  end do
  
  if (iam == 0 ) write(6,*) 'Total memory for t2 storage', 2.*memory, 'GByte' 


end SUBROUTINE setup_t_amplitudes



SUBROUTINE diag_2particles
  USE single_particle_orbits
  USE constants
  use one_body_operators
  use t2_storage
  use configurations
  USE KSPACE 
  use chiral_potentials 

  IMPLICIT NONE
  
  INTEGER :: a,b, p,q,r,s,number_channels, nx,ny,nz,sz,tz
  INTEGER :: ij_confs, bra,ket, channel, ndim
  real*8 :: delta, kin_energy, memory
  integer, allocatable :: numstates(:,:,:,:)
  INTEGER :: nxa, nya, nza, nxb, nyb, nzb,tza,tzb,sza,szb, ichan   
  INTEGER :: numchannels(1), channel1
  INTEGER, ALLOCATABLE :: numchannels1(:)
  integer, allocatable :: locate_2part_channel(:,:,:,:)
  complex*16, allocatable :: hf_mtx(:,:), hf_vect(:,:), kinetic_energy(:,:), hf_eigen(:)
  

  allocate( locate_2part_channel(-1:1, nx_min:nx_max, ny_min:ny_max, nz_min:nz_max) )
  allocate( numstates(-1:1, nx_min:nx_max, ny_min:ny_max, nz_min:nz_max) )
  locate_2part_channel =  0 
  numstates = 0 
  do a = 1,tot_orbs
     nxa = all_orbit%nx(a)
     nya = all_orbit%ny(a)
     nza = all_orbit%nz(a)
     tza=all_orbit%itzp(a)
     sza=all_orbit%szp(a)
     
     do b = a+1, tot_orbs
        nxb = all_orbit%nx(b)
        nyb = all_orbit%ny(b)
        nzb = all_orbit%nz(b)
        tzb = all_orbit%itzp(b)
        szb = all_orbit%szp(b)
        
        tz = (tza + tzb)/2 
        sz = (sza + szb)/2 
        Nx = nxa + nxb 
        ny = nya + nyb 
        nz = nza + nzb
        
        numstates(tz, nx, ny, nz) = numstates(tz, nx, ny, nz)  + 1 
                
     end do
  end do
  
  number_channels = 0 
  numchannels = 0 
  !     loop over isospin projection
  DO tz=-1,1 
     !     loop over total momentum 
     DO nx=nx_min, nx_max 
        !     loop over total momentum 
        DO ny=ny_min, ny_max 
           !     loop over total momentum 
           DO nz=nz_min, nz_max 

              if ( numstates(tz, nx, ny, nz) /= 0 ) then
                 numchannels(1) = numchannels(1) + 1 
                 locate_2part_channel(tz, nx, ny, nz) = numchannels(1)
              end if
              
           end DO
        end DO
     end DO
  end DO
  
  ALLOCATE( lookup_pqrs_configs(1:1,numchannels(1)))
  memory = 0.d0 
  !     loop over isospin projection
  DO tz=-1,1 
     !     loop over total momentum 
     DO nx=nx_min, nx_max 
        !     loop over total momentum 
        DO ny=ny_min, ny_max 
           !     loop over total momentum 
           DO nz=nz_min, nz_max 
              
              channel1 = locate_2part_channel( tz, nx, ny, nz) 
              if ( channel1 /= 0 ) then 
                 ij_confs = numstates(tz, nx, ny, nz) 
                 ALLOCATE( lookup_pqrs_configs(1,channel1)%ival(2,ij_confs) )
                 
                 memory = memory + ij_confs*ij_confs*8./1.e9 
              end if
              
           end DO
        end DO
     end DO
  end DO
  
  write(6,*) 'memory', memory, 'Gbyte'
  
  allocate( numchannels1( maxval( locate_2part_channel(:,:,:,:) ) ) )
  numchannels1 = 0
  do a = 1, tot_orbs
     nxa = all_orbit%nx(a)
     nya = all_orbit%ny(a)
     nza = all_orbit%nz(a)
     tza=all_orbit%itzp(a)
     sza=all_orbit%szp(a)
     
     do b = a+1, tot_orbs
        nxb = all_orbit%nx(b)
        nyb = all_orbit%ny(b)
        nzb = all_orbit%nz(b)
        tzb = all_orbit%itzp(b)
        szb = all_orbit%szp(b)
        
        tz = (tza + tzb)/2 
        sz = (sza + szb)/2 
        Nx = nxa + nxb 
        ny = nya + nyb 
        nz = nza + nzb
        
        
        if ( locate_2part_channel(tz, nx, ny, nz) /= 0 ) then 
           
           channel = locate_2part_channel(tz, nx, ny, nz) 
           numchannels1(channel) = numchannels1(channel) + 1 
           ichan = numchannels1(channel)
           
           lookup_pqrs_configs(1,channel)%ival(1,ichan) = a 
           lookup_pqrs_configs(1,channel)%ival(2,ichan) = b 
        end if
        
     end do
  end do
  
  
  tz = 0
  sz = 0   
  nx = 0 
  ny = 0 
  nz = 0 
  channel = locate_2part_channel(tz, nx, ny, nz) 
  
  write(6,*) 'channel', channel,  size(  lookup_pqrs_configs(1,channel)%ival, 2)  
  ndim =  size(  lookup_pqrs_configs(1,channel)%ival, 2)  
  allocate( kinetic_energy(ndim, ndim),  hf_mtx(ndim, ndim), hf_vect(ndim, ndim), hf_eigen(ndim) )
                 
  kinetic_energy = 0.d0; hf_mtx = 0.d0; hf_vect = 0.d0; hf_eigen = 0.d0 
  do bra = 1, size(  lookup_pqrs_configs(1,channel)%ival, 2) 
     p = lookup_pqrs_configs(1,channel)%ival(1,bra) 
     q = lookup_pqrs_configs(1,channel)%ival(2,bra) 
                                   
     do ket = 1, size(  lookup_pqrs_configs(1,channel)%ival, 2) 
        r = lookup_pqrs_configs(1,channel)%ival(1,ket) 
        s = lookup_pqrs_configs(1,channel)%ival(2,ket) 
        if(chiral_delta_flag)THEN
           hf_mtx(bra,ket) = ( chiral_NN_with_delta(p,q,r,s) - chiral_NN_with_delta(p,q,s,r) ) + kin_energy(p,q,r,s)
        else
           hf_mtx(bra,ket) = ( chiral_pot(p,q,r,s) - chiral_pot(p,q,s,r) ) + kin_energy(p,q,r,s)
        end if
        
     end do
  end do
  CALL lapack_diag (hf_mtx, hf_vect, hf_eigen, ndim )
  do p = 1, ndim 
     if ( real(hf_eigen(p)) < -0.001 ) write(66,*) channel, p, real(hf_eigen(p)) 
     if ( real(hf_eigen(p)) < -0.001 ) write(6,*) channel, p, real(hf_eigen(p)) 
  end do
                 
  deallocate( kinetic_energy, hf_mtx, hf_vect, hf_eigen ) 
  
  deallocate( numstates, numchannels1 ) 

end SUBROUTINE diag_2particles
