SUBROUTINE ccd_iter
  USE diis_mod
  USE PARALLEL
  USE CONSTANTS
  USE one_body_operators
  USE t2_storage

  IMPLICIT NONE
  INTEGER :: count, itimes, ntimes
  REAL*8 :: ener1, ener2, dener
  logical :: switch 
  !
  ! Setup initial t1 and t2 amplitudes
  !
  call setup_tchannels 
  call setup_channel_structures
  
  read(5,*);read(5,*) switch 
  if ( switch ) then  
     call diag_2particles
     write(*,*) 'diagonalization done'
     stop
  end if
  
  !
  ! read in subspace dimension and at what step diis is called
  !
  read(5,*);read(5,*) subspace, diis_step
  write(6,*) subspace, diis_step
  
  call ccd_energy_save(ener1 ) 
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
    
  
  !
  ! here is the main iteration loop
  !
  nstep=0
  count = 0
  ener2 = 0.d0
  ntimes =   1000
  dener=4.0
  do itimes=1, ntimes
     if(abs(dener) > 1.0e-8 )then
        
        if ( iam == 0 ) write(6,*) 't2 eqn'
        call t2_intermediate
        
        count = count + 1
        nstep = nstep + 1
        
        !call linear_t2(0.d0) !
        call diis_t1t2(count)
        call ccd_energy_save(ener2)
        
        dener=ener2-ener1
        ener1=ener2
        
        if(iam ==0)write(6,'(a,2i5,2e16.8)')'dener',iam,itimes,real(dener) 
     end if
  end do
  
  
  
  
  
END SUBROUTINE ccd_iter

!
!
!
SUBROUTINE setup_tchannels 
  USE single_particle_orbits
  USE constants
  use one_body_operators
  use t2_storage
  use configurations
  USE KSPACE 
  use CHIRAL_POTENTIALS

  IMPLICIT NONE
  
  TYPE (configuration_descriptor) :: t2_ket_configs 
  TYPE (configuration_descriptor) :: t2_bra_configs 
  INTEGER :: i,j,a,b, number_channels, nx,ny,nz,sz,tz,  sumx, sumy, sumz  
  INTEGER :: k1,k2,k3,k4,k5, ab_confs, ij_confs, bra,ket, channel
  real*8 :: memory, vmom_yukawa, denom
  integer, allocatable :: numstates(:,:,:,:,:,:)
  INTEGER :: nxa, nya, nza, nxb, nyb, nzb,tza,tzb,sza,szb, ichan   
  INTEGER :: numchannels(1:6), channel1, channel2, channel3, channel4, channel5, channel6
  INTEGER, ALLOCATABLE :: numchannels1(:), numchannels2(:), numchannels3(:), & 
       numchannels4(:), numchannels5(:), numchannels6(:)
  
  
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
  allocate( locate_t2channel(-1:1, -1:1, nx_min:nx_max, ny_min:ny_max, nz_min:nz_max) )
  allocate( locate_channel(1:6,-1:1, -1:1, nx_min:nx_max, ny_min:ny_max, nz_min:nz_max) )
  allocate( numstates(1:3,-1:1, -1:1, nx_min:nx_max, ny_min:ny_max, nz_min:nz_max) )
  locate_t2channel =  0 
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
        
        if ( a <= below_ef .and. b <= below_ef ) numstates(1,tz, sz, nx, ny, nz) = numstates(1,tz, sz, nx, ny, nz)  + 1 
        if ( a <= below_ef .and. b >  below_ef ) numstates(2,tz, sz, nx, ny, nz) = numstates(2,tz, sz, nx, ny, nz)  + 1 
        if ( a > below_ef  .and. b >  below_ef ) numstates(3,tz, sz, nx, ny, nz) = numstates(3,tz, sz, nx, ny, nz)  + 1 
        
     end do
  end do
  
  number_channels = 0 
  numchannels = 0 
  !     loop over isospin projection
  DO tz=-1,1 
     !     loop over spin projection
     DO sz=-1,1           
        !     loop over total momentum 
        DO nx=nx_min, nx_max 
           !     loop over total momentum 
           DO ny=ny_min, ny_max 
              !     loop over total momentum 
              DO nz=nz_min, nz_max 

                 
                 if ( numstates(1,tz, sz, nx, ny, nz) /= 0 ) then
                    numchannels(1) = numchannels(1) + 1 
                    locate_channel(1,tz, sz, nx, ny, nz) = numchannels(1)
                    !write(6,*) numchannels(1), numstates(1,tz, sz, nx, ny, nz)
                 end if

                 if ( numstates(1,tz, sz, nx, ny, nz)* numstates(2,tz, sz, nx, ny, nz) /= 0 ) then
                    numchannels(2) = numchannels(2) + 1 
                    locate_channel(2,tz, sz, nx, ny, nz) = numchannels(2)
                 end if 
                 
                 if ( numstates(1,tz, sz, nx, ny, nz)* numstates(3,tz, sz, nx, ny, nz) /= 0 ) then
                    numchannels(3) = numchannels(3) + 1 
                    locate_channel(3,tz, sz, nx, ny, nz) = numchannels(3)
                 end if 
                 
                 if ( numstates(2,tz, sz, nx, ny, nz) /= 0 ) then 
                    numchannels(4) = numchannels(4) + 1 
                    locate_channel(4,tz, sz, nx, ny, nz) = numchannels(4)
                 end if
                 
                 if ( numstates(2,tz, sz, nx, ny, nz)* numstates(3,tz, sz, nx, ny, nz) /= 0 ) then
                    numchannels(5) = numchannels(5) + 1 
                    locate_channel(5,tz, sz, nx, ny, nz) = numchannels(5)
                 end if
                 
                 if ( numstates(3,tz, sz, nx, ny, nz) /= 0 ) then
                    numchannels(6) = numchannels(6) + 1 
                    locate_channel(6,tz, sz, nx, ny, nz) = numchannels(6)
                 end if
                 
              end DO
           end DO
        end DO
     end DO
  end DO


  ALLOCATE( lookup_hhhh_configs(1:1,numchannels(1)))
  ALLOCATE( lookup_hhhp_configs(1:2,numchannels(2)))
  ALLOCATE( lookup_hhpp_configs(1:2,numchannels(3)))
  ALLOCATE( lookup_hphp_configs(1:1,numchannels(4)))
  ALLOCATE( lookup_hppp_configs(1:2,numchannels(5)))
  ALLOCATE( lookup_pppp_configs(1:1,numchannels(6)))
  
  memory = 0.d0 
  !     loop over isospin projection
  DO tz=-1,1 
     !     loop over spin projection
     DO sz=-1,1           
        !     loop over total momentum 
        DO nx=nx_min, nx_max 
           !     loop over total momentum 
           DO ny=ny_min, ny_max 
              !     loop over total momentum 
              DO nz=nz_min, nz_max 
                 
                 channel1 = locate_channel(1, tz, sz, nx, ny, nz) 
                 channel2 = locate_channel(2, tz, sz, nx, ny, nz) 
                 channel3 = locate_channel(3, tz, sz, nx, ny, nz) 
                 channel4 = locate_channel(4, tz, sz, nx, ny, nz) 
                 channel5 = locate_channel(5, tz, sz, nx, ny, nz) 
                 channel6 = locate_channel(6, tz, sz, nx, ny, nz) 
                 if ( channel1 /= 0 ) then 
                    ij_confs = numstates(1,tz, sz, nx, ny, nz) 
                    memory = memory + dble(ij_confs)*2.*4./1.e9 
                    
                    ALLOCATE( lookup_hhhh_configs(1,channel1)%ival(2,ij_confs) )
                 end if
                 if ( channel2 /= 0 ) then 
                    ij_confs = numstates(1,tz, sz, nx, ny, nz) 
                    ab_confs = numstates(2,tz, sz, nx, ny, nz) 
                    
                    memory = memory + dble(ij_confs+ab_confs)*2.*4./1.e9 
                    ALLOCATE( lookup_hhhp_configs(1,channel2)%ival(2,ij_confs) )
                    ALLOCATE( lookup_hhhp_configs(2,channel2)%ival(2,ab_confs) )
                 end if
                 if ( channel3 /= 0 ) then 
                    ij_confs = numstates(1,tz, sz, nx, ny, nz) 
                    ab_confs = numstates(3,tz, sz, nx, ny, nz) 
                    
                    memory = memory + dble(ij_confs+ab_confs)*2.*4./1.e9 
                    ALLOCATE( lookup_hhpp_configs(1,channel3)%ival(2,ij_confs) )
                    ALLOCATE( lookup_hhpp_configs(2,channel3)%ival(2,ab_confs) )
                 end if
                 if ( channel4 /= 0 ) then 
                    ij_confs = numstates(2,tz, sz, nx, ny, nz) 
                    
                    memory = memory + dble(ij_confs)*2.*4./1.e9 
                    ALLOCATE( lookup_hphp_configs(1,channel4)%ival(2,ij_confs) )
                 end if
                 if ( channel5 /= 0 ) then 
                    ij_confs = numstates(2,tz, sz, nx, ny, nz) 
                    ab_confs = numstates(3,tz, sz, nx, ny, nz) 
                    
                    memory = memory + dble(ij_confs+ab_confs)*2.*4./1.e9 
                    ALLOCATE( lookup_hppp_configs(1,channel5)%ival(2,ij_confs) )
                    ALLOCATE( lookup_hppp_configs(2,channel5)%ival(2,ab_confs) )
                 end if
                 if ( channel6 /= 0 ) then 
                    ab_confs = numstates(3,tz, sz, nx, ny, nz) 

                    memory = memory + dble(ab_confs)*2.*4./1.e9 
                    ALLOCATE( lookup_pppp_configs(1,channel6)%ival(2,ab_confs) )
                 end if
               
                 
              end DO
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
  allocate( hp_hppp%ival(below_ef,below_ef+1:tot_orbs) )
  allocate( hp_hhhp%ival(below_ef,below_ef+1:tot_orbs) )
  allocate( pp_hhpp%ival(below_ef+1:tot_orbs,below_ef+1:tot_orbs) )
  allocate( pp_hppp%ival(below_ef+1:tot_orbs,below_ef+1:tot_orbs) )
  allocate( pp_pppp%ival(below_ef+1:tot_orbs,below_ef+1:tot_orbs) )
  

  hh_hhhh%ival = 0; hh_hhhp%ival = 0; hh_hhpp%ival = 0; hp_hphp%ival = 0  
  hp_hppp%ival = 0; hp_hhhp%ival = 0; pp_hhpp%ival = 0; pp_hppp%ival = 0  
  pp_pppp%ival = 0  
  
  
  allocate( numchannels1( maxval( locate_channel(1,:,:,:,:,:) ) ) )
  allocate( numchannels2( maxval( locate_channel(2,:,:,:,:,:) ) ) )
  allocate( numchannels3( maxval( locate_channel(3,:,:,:,:,:) ) ) )
  allocate( numchannels4( maxval( locate_channel(4,:,:,:,:,:) ) ) )
  allocate( numchannels5( maxval( locate_channel(5,:,:,:,:,:) ) ) )
  allocate( numchannels6( maxval( locate_channel(6,:,:,:,:,:) ) ) )
  numchannels1 = 0
  numchannels2 = 0
  numchannels3 = 0
  numchannels4 = 0
  numchannels5 = 0
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
        

        if ( locate_channel(1, tz, sz, nx, ny, nz) /= 0 ) then 
           
           channel = locate_channel(1, tz, sz, nx, ny, nz) 
           numchannels1(channel) = numchannels1(channel) + 1 
           ichan = numchannels1(channel)
           
           hh_hhhh%ival(a,b) = ichan
           lookup_hhhh_configs(1,channel)%ival(1,ichan) = a 
           lookup_hhhh_configs(1,channel)%ival(2,ichan) = b 
        end if
        
        if ( locate_channel(2, tz, sz, nx, ny, nz) /= 0 ) then
           channel = locate_channel(2, tz, sz, nx, ny, nz) 
           numchannels2(channel) = numchannels2(channel) + 1 
           ichan = numchannels2(channel)
           
           hh_hhhp%ival(a,b) = ichan
           lookup_hhhp_configs(1,channel)%ival(1,ichan) = a 
           lookup_hhhp_configs(1,channel)%ival(2,ichan) = b 
           
        end if
        
        if ( locate_channel(3, tz, sz, nx, ny, nz) /= 0 ) then
           channel = locate_channel(3, tz, sz, nx, ny, nz) 
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
  numchannels5 = 0
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
        

        if ( locate_channel(2, tz, sz, nx, ny, nz) /= 0 ) then
           channel = locate_channel(2, tz, sz, nx, ny, nz) 
           numchannels2(channel) = numchannels2(channel) + 1 
           ichan = numchannels2(channel)

           hp_hhhp%ival(a,b) = ichan
           lookup_hhhp_configs(2,channel)%ival(1,ichan) = a 
           lookup_hhhp_configs(2,channel)%ival(2,ichan) = b 
           
        end if
        
        if ( locate_channel(4, tz, sz, nx, ny, nz) /= 0 ) then
           channel = locate_channel(4, tz, sz, nx, ny, nz) 
           numchannels4(channel) = numchannels4(channel) + 1 
           ichan = numchannels4(channel)
           
           hp_hphp%ival(a,b) = ichan
           lookup_hphp_configs(1,channel)%ival(1,ichan) = a 
           lookup_hphp_configs(1,channel)%ival(2,ichan) = b 
           
        end if
        
        if ( locate_channel(5, tz, sz, nx, ny, nz) /= 0 ) then
           channel = locate_channel(5, tz, sz, nx, ny, nz) 
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
  numchannels1 = 0
  numchannels2 = 0
  numchannels3 = 0
  numchannels4 = 0
  numchannels5 = 0
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
        

        if ( locate_channel(3, tz, sz, nx, ny, nz) /= 0 ) then
           channel = locate_channel(3, tz, sz, nx, ny, nz) 
           numchannels3(channel) = numchannels3(channel) + 1 
           ichan = numchannels3(channel)

           pp_hhpp%ival(a,b) = ichan
           lookup_hhpp_configs(2,channel)%ival(1,ichan) = a 
           lookup_hhpp_configs(2,channel)%ival(2,ichan) = b 
           
        end if
        
        if ( locate_channel(5, tz, sz, nx, ny, nz) /= 0 ) then
           channel = locate_channel(5, tz, sz, nx, ny, nz) 
           numchannels5(channel) = numchannels5(channel) + 1 
           ichan = numchannels5(channel)

           pp_hppp%ival(a,b) = ichan
           lookup_hppp_configs(2,channel)%ival(1,ichan) = a 
           lookup_hppp_configs(2,channel)%ival(2,ichan) = b 
       
        end if
        if ( locate_channel(6, tz, sz, nx, ny, nz) /= 0 ) then
           channel = locate_channel(6, tz, sz, nx, ny, nz) 
           numchannels6(channel) = numchannels6(channel) + 1 
           ichan = numchannels6(channel)
           
           pp_pppp%ival(a,b) = ichan
           lookup_pppp_configs(1,channel)%ival(1,ichan) = a 
           lookup_pppp_configs(1,channel)%ival(2,ichan) = b 
       
        end if
     end do
  end do
  

  

  
  !
  ! setup t-amps 
  !

  number_channels = 0
  memory = 0.d0 
  !     loop over isospin projection
  DO tz=-1,1 
     !     loop over spin projection
     DO sz=-1,1           
        !     loop over total momentum 
        DO nx=nx_min, nx_max 
           !     loop over total momentum 
           DO ny=ny_min, ny_max 
              !     loop over total momentum 
              DO nz=nz_min, nz_max 
                 
                 !
                 ! call number of bra_confs
                 !
                 CALL  number_confs(nx,ny,nz,sz,tz,1,t2_ket_configs)
                 IF ( t2_ket_configs%number_confs <= 0 ) CYCLE
                 number_channels = number_channels + 1 
                 locate_t2channel(tz, sz, nx, ny, nz) = number_channels
                 
                 
                 CALL  number_confs(nx,ny,nz,sz,tz,3,t2_bra_configs)
                 IF ( t2_bra_configs%number_confs <= 0 ) CYCLE
                 
                 memory = memory + t2_bra_configs%number_confs*t2_ket_configs%number_confs * 8./10**9
                 
              end DO
           end DO
        end DO
     end DO
  end DO
  
  if ( iam == 0 ) write(6,*) 'total number of t2 channels', number_channels, numchannels(3) 
  channels%number_confs  = number_channels
  ALLOCATE( channels%t2_quantum_numbers(5*number_channels) )
  allocate( t2_ccm(number_channels) )
  allocate( t2_ccm_eqn(number_channels) )

  
  if ( iam == 0 ) write(6,*) 'storage for t2 amplitudes', 2.*memory, 'Gbyte' 
  
  !     loop over isospin projection
  DO tz=-1,1 
     !     loop over spin projection
     DO sz=-1,1           
        !     loop over total momentum 
        DO nx=nx_min, nx_max 
           !     loop over total momentum 
           DO ny=ny_min, ny_max 
              !     loop over total momentum 
              DO nz=nz_min, nz_max 
                 
                 number_channels = locate_t2channel(tz, sz, nx, ny, nz) 
                 if ( number_channels == 0 ) cycle 
                 
                 !
                 ! call number of bra_confs
                 !
                 CALL  number_confs(nx,ny,nz,sz,tz,1,t2_ket_configs)
                 IF ( t2_ket_configs%number_confs <= 0 ) CYCLE
                 
                 k1 = number_channels*5 
                 k2 = number_channels*5 - 1
                 k3 = number_channels*5 - 2
                 k4 = number_channels*5 - 3
                 k5 = number_channels*5 - 4
                 channels%t2_quantum_numbers(k1) = nx
                 channels%t2_quantum_numbers(k2) = ny
                 channels%t2_quantum_numbers(k3) = nz
                 channels%t2_quantum_numbers(k4) = sz
                 channels%t2_quantum_numbers(k5) = tz
                 
                 CALL  number_confs(nx,ny,nz,sz,tz,3,t2_bra_configs)
                 IF ( t2_bra_configs%number_confs <= 0 ) CYCLE
                 
                 ab_confs = t2_bra_configs%number_confs 
                 ij_confs = t2_ket_configs%number_confs 
                 
                 
                 allocate( t2_ccm(number_channels)%val(0:ab_confs, 0:ij_confs) )
                 allocate( t2_ccm_eqn(number_channels)%val(0:ab_confs, 0:ij_confs) )
                 
                 t2_ccm(number_channels)%val = 0.d0
                 t2_ccm_eqn(number_channels)%val = 0.d0
                 
                 
              end DO
           end DO
        end DO
     end DO
  end DO
  
  number_channels = channels%number_confs 
  ALLOCATE( lookup_t2_configs(1:3,number_channels))
  memory = 0.d0 
  do channel   = 1, channels%number_confs 
     nx = channels%t2_quantum_numbers(channel*5)
     ny = channels%t2_quantum_numbers(channel*5-1)
     nz = channels%t2_quantum_numbers(channel*5-2)
     sz = channels%t2_quantum_numbers(channel*5-3)
     tz = channels%t2_quantum_numbers(channel*5-4)
     
     CALL  number_confs(nx,ny,nz,sz,tz,1,t2_ket_configs)
     IF ( t2_ket_configs%number_confs <= 0 ) CYCLE
                 
     CALL  number_confs(nx,ny,nz,sz,tz,3,t2_bra_configs)
     IF ( t2_bra_configs%number_confs <= 0 ) CYCLE
                 
     ab_confs = t2_bra_configs%number_confs 
     ij_confs = t2_ket_configs%number_confs 
                
     ALLOCATE(t2_ket_configs%config_ab(t2_ket_configs%number_confs*2))
     ALLOCATE(t2_bra_configs%config_ab(t2_bra_configs%number_confs*2))
     
     CALL setup_confs(nx,ny,nz,sz,tz,1,t2_ket_configs)
     CALL setup_confs(nx,ny,nz,sz,tz,3,t2_bra_configs)
     
     ALLOCATE( lookup_t2_configs(1,channel)%ival(2,ij_confs) )
     ALLOCATE( lookup_t2_configs(3,channel)%ival(2,ab_confs) )
     
     memory = memory + dble(2.*ab_confs+2.*ij_confs)*4./1.e9 
     DO bra=1,t2_bra_configs%number_confs
        a= t2_bra_configs%config_ab(bra*2-1)
        b= t2_bra_configs%config_ab(bra*2)
        
        lookup_t2_configs(3,channel)%ival(1,bra) = a
        lookup_t2_configs(3,channel)%ival(2,bra) = b
     end DO
     
     DO ket=1,t2_ket_configs%number_confs
        i = t2_ket_configs%config_ab(ket*2-1)
        j = t2_ket_configs%config_ab(ket*2)

        lookup_t2_configs(1,channel)%ival(1,ket) = i
        lookup_t2_configs(1,channel)%ival(2,ket) = j
     end DO
     

     
     DO bra=1,t2_bra_configs%number_confs
        a= t2_bra_configs%config_ab(bra*2-1)
        b= t2_bra_configs%config_ab(bra*2)
        
        DO ket=1,t2_ket_configs%number_confs
           i = t2_ket_configs%config_ab(ket*2-1)
           j = t2_ket_configs%config_ab(ket*2)
           
           denom = -fock_mtx(a,a) - fock_mtx(b,b) + fock_mtx(i,i) + fock_mtx(j,j)  
           if ( abs( denom ) == 0.d0 ) cycle 
           
           t2_ccm(channel)%val(bra,ket) =   chiral_pot(a,b,i,j) /denom
           t2_ccm_eqn(channel)%val(bra,ket) = chiral_pot(a,b,i,j) /denom
           
                       
        end DO
     end DO
  end do
  
  
  do channel   = 1, channels%number_confs 
     nx = channels%t2_quantum_numbers(channel*5)
     ny = channels%t2_quantum_numbers(channel*5-1)
     nz = channels%t2_quantum_numbers(channel*5-2)
     sz = channels%t2_quantum_numbers(channel*5-3)
     tz = channels%t2_quantum_numbers(channel*5-4)
     
     CALL  number_confs(nx,ny,nz,sz,tz,2,t2_ket_configs)
     IF ( t2_ket_configs%number_confs <= 0 ) CYCLE
     
     ALLOCATE(t2_ket_configs%config_ab(t2_ket_configs%number_confs*2))
     ab_confs = t2_ket_configs%number_confs
     CALL setup_confs(nx,ny,nz,sz,tz,2,t2_ket_configs)
     
     ALLOCATE( lookup_t2_configs(2,channel)%ival(2,ab_confs) )
     
     memory = memory + dble(2.*ab_confs)*4./1.e9 
     DO ket=1,t2_ket_configs%number_confs
        i= t2_ket_configs%config_ab(ket*2-1)
        a= t2_ket_configs%config_ab(ket*2)
        
        lookup_t2_configs(2,channel)%ival(1,ket) = i
        lookup_t2_configs(2,channel)%ival(2,ket) = a
     end DO
  end do
  

  if ( iam == 0 ) write(6,*) 'Memory for lookup ab configs', memory, 'GByte'



  
  deallocate( numstates, numchannels1, numchannels2, & 
       numchannels3, numchannels4, numchannels5, numchannels6 )

end SUBROUTINE setup_tchannels



SUBROUTINE setup_channel_structures
  USE single_particle_orbits
  USE constants
  use one_body_operators
  use t2_storage
  use configurations
  use kspace 
  use chiral_potentials

  IMPLICIT NONE
  
  INTEGER :: number_channels, nx,ny,nz,sz,tz
  INTEGER :: k1,k2,k3,k4,k5, dim1, dim2
  INTEGER :: channel, channel2, i,j,k,l,a,b,c,d, bra,ket
  REAL*8 :: memory 

  
end SUBROUTINE setup_channel_structures



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
  real*8 :: vmom_yukawa, delta, kin_energy, memory
  integer, allocatable :: numstates(:,:,:,:,:)
  INTEGER :: nxa, nya, nza, nxb, nyb, nzb,tza,tzb,sza,szb, ichan   
  INTEGER :: numchannels(1), channel1
  INTEGER, ALLOCATABLE :: numchannels1(:)
  integer, allocatable :: locate_2part_channel(:,:,:,:,:)
  complex*16, allocatable :: hf_mtx(:,:), hf_vect(:,:), kinetic_energy(:,:), hf_eigen(:)
  

  allocate( locate_2part_channel(-1:1, -1:1, nx_min:nx_max, ny_min:ny_max, nz_min:nz_max) )
  allocate( numstates(-1:1, -1:1, nx_min:nx_max, ny_min:ny_max, nz_min:nz_max) )
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
        
        numstates(tz, sz, nx, ny, nz) = numstates(tz, sz, nx, ny, nz)  + 1 
                
     end do
  end do
  
  number_channels = 0 
  numchannels = 0 
  !     loop over isospin projection
  DO tz=-1,1 
     !     loop over spin projection
     DO sz=-1,1           
        !     loop over total momentum 
        DO nx=nx_min, nx_max 
           !     loop over total momentum 
           DO ny=ny_min, ny_max 
              !     loop over total momentum 
              DO nz=nz_min, nz_max 

                 if ( numstates(tz, sz, nx, ny, nz) /= 0 ) then
                    numchannels(1) = numchannels(1) + 1 
                    locate_2part_channel(tz, sz, nx, ny, nz) = numchannels(1)
                 end if
                 
              end DO
           end DO
        end DO
     end DO
  end DO
  
  ALLOCATE( lookup_pqrs_configs(1:1,numchannels(1)))
  memory = 0.d0 
  !     loop over isospin projection
  DO tz=-1,1 
     !     loop over spin projection
     DO sz=-1,1           
        !     loop over total momentum 
        DO nx=nx_min, nx_max 
           !     loop over total momentum 
           DO ny=ny_min, ny_max 
              !     loop over total momentum 
              DO nz=nz_min, nz_max 
                 
                 channel1 = locate_2part_channel( tz, sz, nx, ny, nz) 
                 if ( channel1 /= 0 ) then 
                    ij_confs = numstates(tz, sz, nx, ny, nz) 
                    ALLOCATE( lookup_pqrs_configs(1,channel1)%ival(2,ij_confs) )
                    
                    memory = memory + ij_confs*ij_confs*8./1.e9 
                 end if
                 
              end DO
           end DO
        end DO
     end DO
  end DO
  
  write(6,*) 'memory', memory, 'Gbyte'
  
  allocate( numchannels1( maxval( locate_2part_channel(:,:,:,:,:) ) ) )
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
        
        
        if ( locate_2part_channel(tz, sz, nx, ny, nz) /= 0 ) then 
           
           channel = locate_2part_channel(tz, sz, nx, ny, nz) 
           numchannels1(channel) = numchannels1(channel) + 1 
           ichan = numchannels1(channel)
           
           lookup_pqrs_configs(1,channel)%ival(1,ichan) = a 
           lookup_pqrs_configs(1,channel)%ival(2,ichan) = b 
        end if
        
     end do
  end do
  

  tz = 1
  sz = 0   
  nx = 0 
  ny = 0 
  nz = 0 
  channel = locate_2part_channel(tz, sz, nx, ny, nz) 
  
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
        
        hf_mtx(bra,ket) = (chiral_pot(p,q,r,s)) + kin_energy(p,q,r,s)
                       
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
