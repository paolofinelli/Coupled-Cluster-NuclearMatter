!
!
!
SUBROUTINE setup_v3nf_channel_structures
  USE single_particle_orbits
  USE constants
  use one_body_operators
  use t2_storage
  use configurations
  USE KSPACE 
  use CHIRAL_POTENTIALS
  use t3_constants
  
  IMPLICIT NONE
  INTEGER :: i,j,k, sumx, sumy, sumz, a,b,c,nxa,nxb,nxc,nya,nyb,nyc,nza,nzb,nzc,tza,tzb,tzc,sza,szb,szc,tz,sz,nx,ny,nz
  INTEGER :: d, nxd, nyd, nzd, szd, tzd
  INTEGER :: numchannels, ij_confs, ab_confs, channel, sz_min, sz_max, tz_min, tz_max, sumtz, sumsz, ichan
  integer(8), allocatable :: numstates(:,:,:,:,:)
  INTEGER, ALLOCATABLE :: numchannels1(:),amin(:), amax(:), bmin(:), bmax(:), cmin(:), cmax(:) 
  real*8 :: mem, mem_3nf,delta 
  INTEGER :: k1,k2,k3,k4,k5, dim1, bra, ket, ket2 
  INTEGER(8) :: total_work, work_per_proc, curr_work, curr_proc, curr_channel
  integer :: processor_work(num_procs), mtx_dim, mtx_dim_work, temp, local_ch_id
  integer(8) :: row_dim, column_dim, ket_dim, bra_dim
  INTEGER :: bra_min, bra_max, ket_min, ket_max, number_channels, dim2, channel2
  LOGICAL :: has_ch_been_added
  REAL*8, allocatable :: work(:), memory(:)
  integer(8):: ltot_work, tot_work, work_pr_proc

  allocate( work(1:num_procs) )
  allocate( memory(1:num_procs) )
  work = 0.d0
  memory = 0.d0 
  
  nx_min3 = 1000
  nx_max3 = -1000
  ny_min3 = 1000
  ny_max3 = -1000
  nz_min3 = 1000
  nz_max3 = -1000
  tz_min3 = 100 
  tz_max3 = -100
  sz_min3 = 100 
  sz_max3 = -100
  do i = 1, tot_orbs 
     do j = 1, tot_orbs 
        do k = 1, tot_orbs 
           
           sumx = all_orbit%nx(i) +  all_orbit%nx(j) +  all_orbit%nx(k) 
           sumy = all_orbit%ny(i) +  all_orbit%ny(j) +  all_orbit%ny(k) 
           sumz = all_orbit%nz(i) +  all_orbit%nz(j) +  all_orbit%nz(k) 
           
           sumtz = all_orbit%itzp(i) +  all_orbit%itzp(j) +  all_orbit%itzp(k) 
           sumsz = all_orbit%szp(i) +  all_orbit%szp(j) +  all_orbit%szp(k) 
           
           
           if ( sumx > nx_max3 ) nx_max3 = sumx 
           if ( sumy > ny_max3 ) ny_max3 = sumy 
           if ( sumz > nz_max3 ) nz_max3 = sumz 
           if ( sumx < nx_min3 ) nx_min3 = sumx 
           if ( sumy < ny_min3 ) ny_min3 = sumy 
           if ( sumz < nz_min3 ) nz_min3 = sumz 
           if ( sumtz > tz_max3 ) tz_max3 = sumtz 
           if ( sumsz > sz_max3 ) sz_max3 = sumsz 
           if ( sumtz < tz_min3 ) tz_min3 = sumtz 
           if ( sumsz < sz_min3 ) sz_min3 = sumsz 
           
           
           
        end do
     end do
  end do
!!$  
!!$  do a = 1,below_ef
!!$     nxa = all_orbit%nx(a)
!!$     nya = all_orbit%ny(a)
!!$     nza = all_orbit%nz(a)
!!$     tza=all_orbit%itzp(a)
!!$     sza=all_orbit%szp(a)
!!$     
!!$     do b = below_ef+1, tot_orbs
!!$        nxb = all_orbit%nx(b)
!!$        nyb = all_orbit%ny(b)
!!$        nzb = all_orbit%nz(b)
!!$        tzb = all_orbit%itzp(b)
!!$        szb = all_orbit%szp(b)
!!$        
!!$        do c = below_ef+1, tot_orbs
!!$           nxc = all_orbit%nx(c)
!!$           nyc = all_orbit%ny(c)
!!$           nzc = all_orbit%nz(c)
!!$           tzc = all_orbit%itzp(c)
!!$           szc = all_orbit%szp(c)
!!$        
!!$           do d = 1, below_ef
!!$              nxd = all_orbit%nx(d)
!!$              nyd = all_orbit%ny(d)
!!$              nzd = all_orbit%nz(d)
!!$              tzd = all_orbit%itzp(d)
!!$              szd = all_orbit%szp(d)
!!$              
!!$              sumtz = (tza + tzb + tzc - tzd)
!!$              sumsz = (sza + szb + szc - szd)
!!$              sumx = nxa + nxb + nxc - nxd 
!!$              sumy = nya + nyb + nyc - nyd 
!!$              sumz = nza + nzb + nzc - nzd 
!!$              
!!$              if ( sumx > nx_max3 ) nx_max3 = sumx 
!!$              if ( sumy > ny_max3 ) ny_max3 = sumy 
!!$              if ( sumz > nz_max3 ) nz_max3 = sumz 
!!$              if ( sumx < nx_min3 ) nx_min3 = sumx 
!!$              if ( sumy < ny_min3 ) ny_min3 = sumy 
!!$              if ( sumz < nz_min3 ) nz_min3 = sumz 
!!$              if ( sumtz > tz_max3 ) tz_max3 = sumtz 
!!$              if ( sumsz > sz_max3 ) sz_max3 = sumsz 
!!$              if ( sumtz < tz_min3 ) tz_min3 = sumtz 
!!$              if ( sumsz < sz_min3 ) sz_min3 = sumsz 
!!$              
!!$              
!!$           end do
!!$        end do
!!$     end do
!!$  end do
!!$ 
  
  if ( iam == 0 ) write(6,*) nx_max3,ny_max3, nz_max3, tz_min3,tz_max3
  allocate( locate_v3channel(6, tz_min3:tz_max3, nx_min3:nx_max3, ny_min3:ny_max3, nz_min3:nz_max3) )
  locate_v3channel = 0  
  
  allocate( numstates(1:6,tz_min3:tz_max3, nx_min3:nx_max3, ny_min3:ny_max3, nz_min3:nz_max3) )
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
  
        do c = 1, tot_orbs
           nxc = all_orbit%nx(c)
           nyc = all_orbit%ny(c)
           nzc = all_orbit%nz(c)
           tzc = all_orbit%itzp(c)
           szc = all_orbit%szp(c)
           
           tz = (tza + tzb + tzc)
           sz = (sza + szb + szc)
           Nx = nxa + nxb + nxc
           ny = nya + nyb + nyc
           nz = nza + nzb + nzc
             
           ! hhh
           if ( a <= below_ef .and. b <= below_ef .and. c <= below_ef ) & 
                numstates(1,tz, nx, ny, nz) = numstates(1,tz, nx, ny, nz)  + 1 
           ! hhp
           if ( a <= below_ef .and. b <= below_ef .and. c > below_ef  ) & 
                numstates(2,tz, nx, ny, nz) = numstates(2,tz, nx, ny, nz)  + 1 
           ! hpp
           if ( a <= below_ef .and. b > below_ef  .and. c > below_ef  ) & 
                numstates(3,tz, nx, ny, nz) = numstates(3,tz, nx, ny, nz)  + 1 
           ! ppp
           if ( a > below_ef  .and. b >  below_ef .and. c > below_ef  ) & 
                numstates(4,tz, nx, ny, nz) = numstates(4,tz, nx, ny, nz)  + 1 
           
        end do
     end do
  end do
  
!!$  !
!!$  ! hpph-1 
!!$  !
!!$  do a = 1,below_ef
!!$     nxa = all_orbit%nx(a)
!!$     nya = all_orbit%ny(a)
!!$     nza = all_orbit%nz(a)
!!$     tza=all_orbit%itzp(a)
!!$     sza=all_orbit%szp(a)
!!$     
!!$     do b = below_ef+1, tot_orbs
!!$        nxb = all_orbit%nx(b)
!!$        nyb = all_orbit%ny(b)
!!$        nzb = all_orbit%nz(b)
!!$        tzb = all_orbit%itzp(b)
!!$        szb = all_orbit%szp(b)
!!$        
!!$        do c = below_ef+1, tot_orbs
!!$           nxc = all_orbit%nx(c)
!!$           nyc = all_orbit%ny(c)
!!$           nzc = all_orbit%nz(c)
!!$           tzc = all_orbit%itzp(c)
!!$           szc = all_orbit%szp(c)
!!$        
!!$           do d = 1, below_ef
!!$              nxd = all_orbit%nx(d)
!!$              nyd = all_orbit%ny(d)
!!$              nzd = all_orbit%nz(d)
!!$              tzd = all_orbit%itzp(d)
!!$              szd = all_orbit%szp(d)
!!$              
!!$              tz = (tza + tzb + tzc - tzd)
!!$              sz = (sza + szb + szc - szd)
!!$              Nx = nxa + nxb + nxc - nxd 
!!$              ny = nya + nyb + nyc - nyd 
!!$              nz = nza + nzb + nzc - nzd 
!!$              
!!$              numstates(5,tz, nx, ny, nz) = numstates(5,tz, nx, ny, nz)  + 1 
!!$              
!!$           end do
!!$        end do
!!$     end do
!!$  end do
!!$  !
!!$  ! pp 
!!$  ! 
!!$  do a = below_ef+1, tot_orbs
!!$     nxa = all_orbit%nx(a)
!!$     nya = all_orbit%ny(a)
!!$     nza = all_orbit%nz(a)
!!$     tza=all_orbit%itzp(a)
!!$     sza=all_orbit%szp(a)
!!$     
!!$     do b = below_ef+1, tot_orbs
!!$        nxb = all_orbit%nx(b)
!!$        nyb = all_orbit%ny(b)
!!$        nzb = all_orbit%nz(b)
!!$        tzb = all_orbit%itzp(b)
!!$        szb = all_orbit%szp(b)
!!$              
!!$        tz = (tza + tzb)
!!$        sz = (sza + szb)
!!$        Nx = nxa + nxb
!!$        ny = nya + nyb
!!$        nz = nza + nzb
!!$        
!!$        numstates(6,tz, nx, ny, nz) = numstates(6,tz, nx, ny, nz)  + 1 
!!$        
!!$     end do
!!$  end do
!!$  
!!$
!!$  !
!!$  ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!!$  ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!!$  !
!!$  ! setup <hpph-1||pp> channels and procs mapping
!!$  !
!!$  numchannels = 0 
!!$  !     loop over isospin projection
!!$  DO tz=tz_min3,tz_max3,2 
!!$     !     loop over total momentum 
!!$     DO nx=nx_min3, nx_max3 
!!$        !     loop over total momentum 
!!$        DO ny=ny_min3, ny_max3 
!!$           !     loop over total momentum 
!!$           DO nz=nz_min3, nz_max3 
!!$              
!!$              if ( numstates(5,tz, nx, ny, nz) == 0 ) cycle
!!$              if ( numstates(6,tz, nx, ny, nz) == 0 ) cycle
!!$              !
!!$              ! check that hpph-1 channel is restricted to hhpp channels. 
!!$              !
!!$              if ( locate_channel(3,tz/2, nx, ny, nz) == 0 ) cycle 
!!$              numchannels = numchannels + 1 
!!$              locate_v3channel(5,tz, nx, ny, nz) = numchannels
!!$              
!!$           end DO
!!$        end DO
!!$     end DO
!!$  end DO
!!$  
!!$  ltot_work=0
!!$  tot_work = 0
!!$  
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
!!$     
!!$     
!!$     dim1 = numstates(5,tz*2, nx, ny, nz) 
!!$     dim2 = numstates(6,tz*2, nx, ny, nz) 
!!$     
!!$     tot_work = tot_work + dim1*dim2
!!$     ltot_work = ltot_work + int(dim1,8)*int(dim2, 8)
!!$     
!!$  end do
!!$
!!$  
!!$  work_pr_proc = int( tot_work/num_procs )
!!$  work_pr_proc = int( ltot_work/int(num_procs,8) )
!!$  
!!$  if ( iam == 0 ) write(6,*) ltot_work, work_pr_proc, number_channels 
!!$  
!!$  work = 0.d0
!!$  memory = 0.d0
!!$  
!!$  number_channels = channels%number_pppp_confs 
!!$  allocate( mapping_ph_hpphpp(1:num_procs, 1:number_channels, 1:5) )
!!$  mapping_ph_hpphpp = 0 
!!$  
!!$  curr_work = 0
!!$  curr_proc = 0
!!$  local_ch_id = 1
!!$  number_channels = 0
!!$  do channel   = 1, channels%number_pppp_confs 
!!$     nx = channels%pppp_quantum_numbers(channel*4)
!!$     ny = channels%pppp_quantum_numbers(channel*4-1)
!!$     nz = channels%pppp_quantum_numbers(channel*4-2)
!!$     tz = channels%pppp_quantum_numbers(channel*4-3)
!!$     
!!$     
!!$     !
!!$     ! check that pppp channel is restricted to hhpp channels. 
!!$     !
!!$     channel2 = locate_channel(3,tz, nx, ny, nz)
!!$     if ( channel2 == 0 ) cycle 
!!$
!!$     dim1 = numstates(6,tz*2, nx, ny, nz) 
!!$     dim2 = numstates(5,tz*2, nx, ny, nz) 
!!$          
!!$     bra_dim = dim1 
!!$     ket_dim = dim2
!!$     
!!$     number_channels = number_channels + 1 
!!$     has_ch_been_added = .FALSE. 
!!$           
!!$     i = 0
!!$     DO WHILE ( i < ket_dim )
!!$              
!!$        i = i + 1
!!$        temp = curr_work + bra_dim
!!$              
!!$        if ( temp <= work_pr_proc .or. curr_proc == num_procs - 1 ) then
!!$           curr_work = temp
!!$           if (  has_ch_been_added ) then
!!$                    
!!$              mapping_ph_hpphpp(curr_proc+1,local_ch_id,3 ) = i 
!!$              
!!$           else
!!$              
!!$              
!!$              has_ch_been_added = .true. 
!!$              mapping_ph_hpphpp(curr_proc+1,local_ch_id,1 ) = number_channels 
!!$              mapping_ph_hpphpp(curr_proc+1,local_ch_id,2 ) = i
!!$              mapping_ph_hpphpp(curr_proc+1,local_ch_id,3 ) = i
!!$              mapping_ph_hpphpp(curr_proc+1,local_ch_id,4 ) = 1
!!$              mapping_ph_hpphpp(curr_proc+1,local_ch_id,5 ) = bra_dim
!!$              
!!$           end if
!!$
!!$           if ( i == ket_dim ) then
!!$              
!!$              local_ch_id = local_ch_id + 1
!!$              
!!$           end if
!!$           work(curr_proc+1) = curr_work
!!$           
!!$        else 
!!$                 
!!$           has_ch_been_added = .false.
!!$           work(curr_proc+1) = curr_work
!!$           curr_work = 0
!!$           curr_proc = curr_proc + 1
!!$           i = i - 1
!!$           
!!$        end if
!!$     end DO
!!$  end DO
!!$  
!!$  number_channels = channels%number_pppp_confs 
!!$  work = 0.d0
!!$  do i = 1, num_procs
!!$     do local_ch_id = 1, number_channels 
!!$        if ( mapping_ph_hpphpp(i,local_ch_id,1) == 0 ) cycle
!!$        
!!$        row_dim = mapping_ph_hpphpp(i,local_ch_id,5)-mapping_ph_hpphpp(i,local_ch_id,4)+1
!!$        column_dim = mapping_ph_hpphpp(i,local_ch_id,3)-mapping_ph_hpphpp(i,local_ch_id,2)+1
!!$        work(i) = work(i) + real( row_dim*column_dim )
!!$        
!!$     end do
!!$  end do
!!$  
!!$  call mpi_barrier(mpi_comm_world,ierror)
!!$  number_channels = channels%number_pppp_confs 
!!$  allocate( check_my_channel_ph_hpphpp(number_channels) ) 
!!$  allocate( my_ph_hpphpp_channel_low(0:num_procs-1), my_ph_hpphpp_channel_high(0:num_procs-1) )
!!$  
!!$  my_ph_hpphpp_channel_low(iam)  = number_channels + 1
!!$  my_ph_hpphpp_channel_high(iam) = 0 
!!$  check_my_channel_ph_hpphpp = 0 
!!$  do local_ch_id = 1, number_channels 
!!$     if ( mapping_ph_hpphpp(iam+1,local_ch_id,1) == 0 ) cycle
!!$     
!!$     check_my_channel_ph_hpphpp(local_ch_id) = 1 
!!$     if ( local_ch_id > my_ph_hpphpp_channel_high(iam) ) my_ph_hpphpp_channel_high(iam) = local_ch_id
!!$     if ( local_ch_id < my_ph_hpphpp_channel_low(iam) ) my_ph_hpphpp_channel_low(iam) = local_ch_id
!!$     
!!$  end do
!!$  
!!$  call mpi_barrier(mpi_comm_world, ierror)
!!$  allocate( v3nf_ph_hpphpp(my_ph_hpphpp_channel_low(iam):my_ph_hpphpp_channel_high(iam)) )
!!$  call mpi_barrier(mpi_comm_world, ierror)
!!$  
!!$  mem = 0.d0  
!!$  mem_3nf = 0.d0 
!!$  do channel = my_ph_hpphpp_channel_low(iam), my_ph_hpphpp_channel_high(iam)
!!$     nx = channels%pppp_quantum_numbers(channel*4)
!!$     ny = channels%pppp_quantum_numbers(channel*4-1)
!!$     nz = channels%pppp_quantum_numbers(channel*4-2)
!!$     tz = channels%pppp_quantum_numbers(channel*4-3)
!!$     
!!$     local_ch_id = channel 
!!$     !
!!$     ! ket side if fully stored on each proc
!!$     !
!!$     bra_min = mapping_ph_hpphpp(iam+1,local_ch_id,2)
!!$     bra_max = mapping_ph_hpphpp(iam+1,local_ch_id,3)
!!$     !
!!$     ! bra side is distributed 
!!$     !
!!$     ket_min = mapping_ph_hpphpp(iam+1,local_ch_id,4)
!!$     ket_max = mapping_ph_hpphpp(iam+1,local_ch_id,5)
!!$     
!!$     ij_confs = numstates(5,tz*2, nx, ny, nz) 
!!$     ab_confs = numstates(6,tz*2, nx, ny, nz) 
!!$     mem = mem + dble(ab_confs)*3.*4./1.e9 
!!$     mem_3nf = mem_3nf + dble((bra_max-bra_min+1)*ab_confs)*16./1.e9 
!!$     
!!$  end DO
!!$  if ( iam == 0 ) write(6,*) 'Total memory for lookup hpphpp-confs', mem, 'Gbyte'
!!$  if ( iam == 0 ) write(6,*) 'Total memory for v3nf_hpphpp', mem_3nf, 'Gbyte'
!!$  
!!$  call mpi_barrier(mpi_comm_world, ierror )
!!$  
!!$  ALLOCATE( lookup_ph_hpphpp_configs(1:2,my_ph_hpphpp_channel_low(iam):my_ph_hpphpp_channel_high(iam)) )
!!$  mem = 0.d0  
!!$  mem_3nf = 0.d0 
!!$  do channel = my_ph_hpphpp_channel_low(iam), my_ph_hpphpp_channel_high(iam)
!!$     nx = channels%pppp_quantum_numbers(channel*4)
!!$     ny = channels%pppp_quantum_numbers(channel*4-1)
!!$     nz = channels%pppp_quantum_numbers(channel*4-2)
!!$     tz = channels%pppp_quantum_numbers(channel*4-3)
!!$     
!!$     if ( locate_channel(3,tz, nx, ny, nz) == 0 ) cycle 
!!$     
!!$     local_ch_id = channel 
!!$     !
!!$     ! ket side if fully stored on each proc
!!$     !
!!$     bra_min = mapping_ph_hpphpp(iam+1,local_ch_id,2)
!!$     bra_max = mapping_ph_hpphpp(iam+1,local_ch_id,3)
!!$     !
!!$     ! bra side is distributed 
!!$     !
!!$     ket_min = mapping_ph_hpphpp(iam+1,local_ch_id,4)
!!$     ket_max = mapping_ph_hpphpp(iam+1,local_ch_id,5)
!!$     
!!$     ij_confs = numstates(5,tz*2, nx, ny, nz) 
!!$     ab_confs = numstates(6,tz*2, nx, ny, nz) 
!!$     mem = mem + dble(ij_confs+ab_confs)*3.*4./1.e9 
!!$     mem_3nf = mem_3nf + dble(ij_confs*ab_confs)*16./1.e9 
!!$     ALLOCATE( lookup_ph_hpphpp_configs(1,channel)%ival(4,ij_confs))
!!$     ALLOCATE( lookup_ph_hpphpp_configs(2,channel)%ival(2,ab_confs))
!!$     
!!$     ALLOCATE( v3nf_ph_hpphpp(channel)%val(bra_min:bra_max, ket_min:ket_max) )
!!$     v3nf_ph_hpphpp(channel)%val = 0.d0 
!!$     
!!$  end DO
!!$  if ( iam == 0 ) write(6,*) 'Total memory for lookup ph hpphpp-confs', mem, 'Gbyte'
!!$  if ( iam == 0 ) write(6,*) 'Total memory for v3nf_ph_hpphpp', mem_3nf, 'Gbyte'
!!$  
!!$  
!!$  allocate( numchannels1( channels%number_pppp_confs ) ) 
!!$  
!!$  !
!!$  ! hpp configurations
!!$  !
!!$  
!!$  numchannels1  = 0
!!$
!!$  !
!!$  ! hpph-1 
!!$  !
!!$  do a = 1,below_ef
!!$     nxa = all_orbit%nx(a)
!!$     nya = all_orbit%ny(a)
!!$     nza = all_orbit%nz(a)
!!$     tza=all_orbit%itzp(a)
!!$     sza=all_orbit%szp(a)
!!$     
!!$     do b = below_ef+1, tot_orbs
!!$        nxb = all_orbit%nx(b)
!!$        nyb = all_orbit%ny(b)
!!$        nzb = all_orbit%nz(b)
!!$        tzb = all_orbit%itzp(b)
!!$        szb = all_orbit%szp(b)
!!$        
!!$        do c = below_ef+1, tot_orbs
!!$           nxc = all_orbit%nx(c)
!!$           nyc = all_orbit%ny(c)
!!$           nzc = all_orbit%nz(c)
!!$           tzc = all_orbit%itzp(c)
!!$           szc = all_orbit%szp(c)
!!$        
!!$           do d = 1, below_ef
!!$              nxd = all_orbit%nx(d)
!!$              nyd = all_orbit%ny(d)
!!$              nzd = all_orbit%nz(d)
!!$              tzd = all_orbit%itzp(d)
!!$              szd = all_orbit%szp(d)
!!$              
!!$              tz = (tza + tzb + tzc - tzd)
!!$              sz = (sza + szb + szc - szd)
!!$              Nx = nxa + nxb + nxc - nxd 
!!$              ny = nya + nyb + nyc - nyd 
!!$              nz = nza + nzb + nzc - nzd 
!!$              
!!$              if ( numstates(5,tz, nx, ny, nz) == 0 ) cycle
!!$              if ( numstates(6,tz, nx, ny, nz) == 0 ) cycle
!!$              
!!$              channel = locate_channel(3,tz/2, nx, ny, nz)
!!$              if ( channel == 0 ) cycle 
!!$              if ( check_my_channel_ph_hpphpp(channel) == 0 ) cycle 
!!$              
!!$              numchannels1(channel) = numchannels1(channel) + 1
!!$              ichan = numchannels1(channel) 
!!$              
!!$              lookup_ph_hpphpp_configs(1,channel)%ival(1,ichan) = a 
!!$              lookup_ph_hpphpp_configs(1,channel)%ival(2,ichan) = b 
!!$              lookup_ph_hpphpp_configs(1,channel)%ival(3,ichan) = c 
!!$              lookup_ph_hpphpp_configs(1,channel)%ival(4,ichan) = d 
!!$              
!!$           end do
!!$        end do
!!$     end do
!!$  end do
!!$
!!$  !
!!$  ! pp 
!!$  ! 
!!$  numchannels1  = 0
!!$  do a = below_ef+1, tot_orbs
!!$     nxa = all_orbit%nx(a)
!!$     nya = all_orbit%ny(a)
!!$     nza = all_orbit%nz(a)
!!$     tza=all_orbit%itzp(a)
!!$     sza=all_orbit%szp(a)
!!$     
!!$     do b = below_ef+1, tot_orbs
!!$        nxb = all_orbit%nx(b)
!!$        nyb = all_orbit%ny(b)
!!$        nzb = all_orbit%nz(b)
!!$        tzb = all_orbit%itzp(b)
!!$        szb = all_orbit%szp(b)
!!$              
!!$        tz = (tza + tzb)
!!$        sz = (sza + szb)
!!$        Nx = nxa + nxb
!!$        ny = nya + nyb
!!$        nz = nza + nzb
!!$
!!$        if ( numstates(5,tz, nx, ny, nz) == 0 ) cycle
!!$        if ( numstates(6,tz, nx, ny, nz) == 0 ) cycle
!!$              
!!$        channel = locate_channel(3,tz/2, nx, ny, nz)
!!$        if ( channel == 0 ) cycle 
!!$        if ( check_my_channel_ph_hpphpp(channel) == 0 ) cycle 
!!$        
!!$        numchannels1(channel) = numchannels1(channel) + 1
!!$        ichan = numchannels1(channel) 
!!$           
!!$        lookup_ph_hpphpp_configs(2,channel)%ival(1,ichan) = a 
!!$        lookup_ph_hpphpp_configs(2,channel)%ival(2,ichan) = b 
!!$                      
!!$     end do
!!$  end do
!!$  
!!$  do channel = my_ph_hpphpp_channel_low(iam), my_ph_hpphpp_channel_high(iam)
!!$     nx = channels%pppp_quantum_numbers(channel*4)
!!$     ny = channels%pppp_quantum_numbers(channel*4-1)
!!$     nz = channels%pppp_quantum_numbers(channel*4-2)
!!$     tz = channels%pppp_quantum_numbers(channel*4-3)
!!$     
!!$     local_ch_id = channel 
!!$     !
!!$     ! ket side if fully stored on each proc
!!$     !
!!$     bra_min = mapping_ph_hpphpp(iam+1,local_ch_id,2)
!!$     bra_max = mapping_ph_hpphpp(iam+1,local_ch_id,3)
!!$     !
!!$     ! bra side is distributed 
!!$     !
!!$     ket_min = mapping_ph_hpphpp(iam+1,local_ch_id,4)
!!$     ket_max = mapping_ph_hpphpp(iam+1,local_ch_id,5)
!!$     
!!$     !$omp parallel default(shared) private(bra,a,b,c,ket,i,j,k)
!!$     !$omp do schedule(dynamic)
!!$     do bra = bra_min, bra_max 
!!$        a = lookup_ph_hpphpp_configs(1,channel)%ival(1,bra)
!!$        b = lookup_ph_hpphpp_configs(1,channel)%ival(2,bra) 
!!$        c = lookup_ph_hpphpp_configs(1,channel)%ival(3,bra) 
!!$        i = lookup_ph_hpphpp_configs(1,channel)%ival(4,bra) 
!!$        
!!$        do ket = 1, size(  lookup_ph_hpphpp_configs(2,channel)%ival, 2) 
!!$           j = lookup_ph_hpphpp_configs(2,channel)%ival(1,ket) 
!!$           k = lookup_ph_hpphpp_configs(2,channel)%ival(2,ket) 
!!$                      
!!$           v3nf_ph_hpphpp(channel)%val(bra,ket) = chiral_3nf_asym(a,b,c,i,j,k)
!!$        end do
!!$     end do
!!$     !$omp end do
!!$     !$omp end parallel
!!$     
!!$  end do
!!$  
!!$
!!$
!!$  deallocate( numchannels1 )

  !
  ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !
  ! setup <hpp||hpp> channels and procs mapping
  !
  numchannels = 0 
  !     loop over isospin projection
  DO tz=tz_min3,tz_max3,2 
     !     loop over total momentum 
     DO nx=nx_min3, nx_max3 
        !     loop over total momentum 
        DO ny=ny_min3, ny_max3 
           !     loop over total momentum 
           DO nz=nz_min3, nz_max3 
              
              if ( numstates(3,tz, nx, ny, nz) /= 0 ) then
                 numchannels = numchannels + 1 
                 locate_v3channel(1,tz, nx, ny, nz) = numchannels
              end if
              
           end DO
        end DO
     end DO
  end DO
  
  
  mem = 0.d0 
  number_channels = 0 
  number_channels = maxval(locate_v3channel(1,:,:,:,:) ) 
  channels%number_hpphpp_confs  = number_channels
  ALLOCATE( channels%hpphpp_quantum_numbers(4*number_channels) )
  mem = mem + dble(number_channels)*4.*4./1.e9 
  
  if ( iam == 0 ) write(6,*) 'Total memory for t3-channels', mem, 'Gbyte' 
  
  !     loop over isospin projection
  DO tz=tz_min3,tz_max3,2 
     !     loop over total momentum 
     DO nx=nx_min3, nx_max3 
        !     loop over total momentum 
        DO ny=ny_min3, ny_max3 
           !     loop over total momentum 
           DO nz=nz_min3, nz_max3 
                 
              number_channels = locate_v3channel(1,tz,nx,ny,nz) 
              if ( number_channels == 0 ) cycle 
              
              k1 = number_channels*4 
              k2 = number_channels*4 - 1
              k3 = number_channels*4 - 2
              k4 = number_channels*4 - 3
              channels%hpphpp_quantum_numbers(k1) = nx
              channels%hpphpp_quantum_numbers(k2) = ny
              channels%hpphpp_quantum_numbers(k3) = nz
              channels%hpphpp_quantum_numbers(k4) = tz
                            
           end DO
        end DO
     end DO
  end DO

  ! 
  !  setup processor mapping and 
  !  calculate work for hhhppp blocks 
  !
  ltot_work=0
  tot_work = 0
  do channel   = 1, channels%number_hpphpp_confs 
     nx = channels%hpphpp_quantum_numbers(channel*4)
     ny = channels%hpphpp_quantum_numbers(channel*4-1)
     nz = channels%hpphpp_quantum_numbers(channel*4-2)
     tz = channels%hpphpp_quantum_numbers(channel*4-3)

     dim1 = numstates(3,tz, nx, ny, nz) 
     dim2 = numstates(3,tz, nx, ny, nz)/2
     
     !dim1 = size(  lookup_t3_configs(1,channel)%ival,2 )
     !dim2 = size(  lookup_t3_configs(2,channel)%ival,2 )
     
     
     tot_work = tot_work + dim1*dim2
     ltot_work = ltot_work + int(dim1,8)*int(dim2, 8)
     
  end do
  
  number_channels = channels%number_hpphpp_confs 
  
  work = 0.d0
  memory = 0.d0
  
  work_pr_proc = int( ltot_work/int(num_procs,8),8 )
  work = 0.d0
  
  allocate( mapping_hpphpp(1:num_procs, 1:number_channels, 1:5) )
  mapping_hpphpp = 0 
  if ( iam == 0 ) write(6,*) ltot_work, work_pr_proc, number_channels 

  curr_work = 0
  curr_proc = 0
  local_ch_id = 1
  number_channels = 0
  do channel   = 1, channels%number_hpphpp_confs 
     nx = channels%hpphpp_quantum_numbers(channel*4)
     ny = channels%hpphpp_quantum_numbers(channel*4-1)
     nz = channels%hpphpp_quantum_numbers(channel*4-2)
     tz = channels%hpphpp_quantum_numbers(channel*4-3)
     
     dim1 = numstates(3,tz, nx, ny, nz)
     dim2 = numstates(3,tz, nx, ny, nz) 
     
     !dim1 = size(  lookup_t3_configs(1,channel)%ival, 2) 
     !dim2 = size(  lookup_t3_configs(2,channel)%ival, 2) 
     
     bra_dim = dim1/2 
     ket_dim = dim2
     
     number_channels = number_channels + 1 
     has_ch_been_added = .FALSE. 
           
     i = 0
     DO WHILE ( i < ket_dim )
              
        i = i + 1
        temp = curr_work + bra_dim
              
        if ( temp <= work_pr_proc .or. curr_proc == num_procs - 1 ) then
           curr_work = temp
           if (  has_ch_been_added ) then
                    
              mapping_hpphpp(curr_proc+1,local_ch_id,3 ) = i 
              
           else
              
              
              has_ch_been_added = .true. 
              mapping_hpphpp(curr_proc+1,local_ch_id,1 ) = number_channels 
              mapping_hpphpp(curr_proc+1,local_ch_id,2 ) = i
              mapping_hpphpp(curr_proc+1,local_ch_id,3 ) = i
              mapping_hpphpp(curr_proc+1,local_ch_id,4 ) = 1
              mapping_hpphpp(curr_proc+1,local_ch_id,5 ) = bra_dim
              
           end if

           if ( i == ket_dim ) then
              
              local_ch_id = local_ch_id + 1
              
           end if
           work(curr_proc+1) = curr_work
           
        else 
                 
           has_ch_been_added = .false.
           work(curr_proc+1) = curr_work
           curr_work = 0
           curr_proc = curr_proc + 1
           i = i - 1
           
        end if
     end DO
  end DO
  
  number_channels = channels%number_hpphpp_confs 
  work = 0.d0
  do i = 1, num_procs
     do local_ch_id = 1, number_channels 
        if ( mapping_hpphpp(i,local_ch_id,1) == 0 ) cycle
        
        row_dim = mapping_hpphpp(i,local_ch_id,5)-mapping_hpphpp(i,local_ch_id,4)+1
        column_dim = mapping_hpphpp(i,local_ch_id,3)-mapping_hpphpp(i,local_ch_id,2)+1
        work(i) = work(i) + real( row_dim*column_dim )
        
     end do
  end do
  
  call mpi_barrier(mpi_comm_world,ierror)
  number_channels = channels%number_hpphpp_confs 
  
  allocate( my_hpphpp_channel_low(0:num_procs-1), my_hpphpp_channel_high(0:num_procs-1) )
  allocate( check_my_channel_hpphpp(1:number_channels) )
  check_my_channel_hpphpp = 0
  
  my_hpphpp_channel_low(iam)  = number_channels + 1
  my_hpphpp_channel_high(iam) = 0 
  do local_ch_id = 1, number_channels 
     if ( mapping_hpphpp(iam+1,local_ch_id,1) == 0 ) cycle
     check_my_channel_hpphpp(local_ch_id) = 1
      
     if ( local_ch_id > my_hpphpp_channel_high(iam) ) my_hpphpp_channel_high(iam) = local_ch_id
     if ( local_ch_id < my_hpphpp_channel_low(iam) ) my_hpphpp_channel_low(iam) = local_ch_id
     
  end do
  
  call mpi_barrier(mpi_comm_world, ierror)
  allocate( v3nf_hpphpp(my_hpphpp_channel_low(iam):my_hpphpp_channel_high(iam)) )
  call mpi_barrier(mpi_comm_world, ierror)
  
  mem = 0.d0  
  mem_3nf = 0.d0 
  do channel = my_hpphpp_channel_low(iam), my_hpphpp_channel_high(iam)
     nx = channels%hpphpp_quantum_numbers(channel*4)
     ny = channels%hpphpp_quantum_numbers(channel*4-1)
     nz = channels%hpphpp_quantum_numbers(channel*4-2)
     tz = channels%hpphpp_quantum_numbers(channel*4-3)
     
     local_ch_id = channel 
     !
     ! ket side if fully stored on each proc
     !
     bra_min = mapping_hpphpp(iam+1,local_ch_id,2)
     bra_max = mapping_hpphpp(iam+1,local_ch_id,3)
     !
     ! bra side is distributed 
     !
     ket_min = mapping_hpphpp(iam+1,local_ch_id,4)
     ket_max = mapping_hpphpp(iam+1,local_ch_id,5)
     
     ij_confs = numstates(3,tz, nx, ny, nz) 
     ab_confs = numstates(3,tz, nx, ny, nz)/2 
     mem = mem + dble(ab_confs)*3.*4./1.e9 
     mem_3nf = mem_3nf + dble((bra_max-bra_min+1)*ab_confs)*16./1.e9 
          
  end DO
  if ( iam == 0 ) write(6,*) 'Total memory for lookup hpphpp-confs', mem, 'Gbyte'
  if ( iam == 0 ) write(6,*) 'Total memory for v3nf_hpphpp', mem_3nf, 'Gbyte'
  
  call mpi_barrier(mpi_comm_world, ierror )
  
  ALLOCATE( lookup_hpphpp_configs(1:1,my_hpphpp_channel_low(iam):my_hpphpp_channel_high(iam)) )
  mem = 0.d0  
  mem_3nf = 0.d0 
  do channel = my_hpphpp_channel_low(iam), my_hpphpp_channel_high(iam)
     nx = channels%hpphpp_quantum_numbers(channel*4)
     ny = channels%hpphpp_quantum_numbers(channel*4-1)
     nz = channels%hpphpp_quantum_numbers(channel*4-2)
     tz = channels%hpphpp_quantum_numbers(channel*4-3)
     
     local_ch_id = channel 
     !
     ! ket side if fully stored on each proc
     !
     bra_min = mapping_hpphpp(iam+1,local_ch_id,2)
     bra_max = mapping_hpphpp(iam+1,local_ch_id,3)
     !
     ! bra side is distributed 
     !
     ket_min = mapping_hpphpp(iam+1,local_ch_id,4)
     ket_max = mapping_hpphpp(iam+1,local_ch_id,5)
     
     ij_confs = numstates(3,tz, nx, ny, nz) 
     ab_confs = numstates(3,tz, nx, ny, nz) 
     mem = mem + dble(ij_confs+ab_confs)*3.*4./1.e9 
     mem_3nf = mem_3nf + dble(ij_confs*ab_confs)*16./1.e9 
     ALLOCATE( lookup_hpphpp_configs(1,channel)%ival(3,ij_confs))
     !ALLOCATE( lookup_hpphpp_configs(2,channel)%ival(3,ab_confs))
     ALLOCATE( v3nf_hpphpp(channel)%val(bra_min:bra_max, ket_min:ket_max) )
     v3nf_hpphpp(channel)%val = 0.d0 
     
  end DO
  if ( iam == 0 ) write(6,*) 'Total memory for lookup hpphpp-confs', mem, 'Gbyte'
  if ( iam == 0 ) write(6,*) 'Total memory for v3nf_hpphpp', mem_3nf, 'Gbyte'
  
  
  allocate( numchannels1( maxval( locate_v3channel(1,:,:,:,:) ) ) )
  
  !
  ! hpp configurations
  !
  
  numchannels1  = 0
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
        do c = below_ef+1, tot_orbs
           nxc = all_orbit%nx(c)
           nyc = all_orbit%ny(c)
           nzc = all_orbit%nz(c)
           tzc = all_orbit%itzp(c)
           szc = all_orbit%szp(c)
           
           
           tz = (tza + tzb + tzc)
           sz = (sza + szb + szc)
           Nx = nxa + nxb + nxc
           ny = nya + nyb + nyc
           nz = nza + nzb + nzc
           
           
           channel = locate_v3channel(1,tz, nx, ny, nz) 
           if ( channel == 0 ) cycle 
           if ( check_my_channel_hpphpp(channel) == 0 ) cycle 
           
           numchannels1(channel) = numchannels1(channel) + 1
           ichan = numchannels1(channel) 
           
           lookup_hpphpp_configs(1,channel)%ival(1,ichan) = a 
           lookup_hpphpp_configs(1,channel)%ival(2,ichan) = b 
           lookup_hpphpp_configs(1,channel)%ival(3,ichan) = c 
           
        end do
     end do
  end do
  
  
  if ( iam == 0 ) write(6,*) 'ok1'
  call mpi_barrier(mpi_comm_world, ierror) 
  
  do channel = my_hpphpp_channel_low(iam), my_hpphpp_channel_high(iam)
     nx = channels%hpphpp_quantum_numbers(channel*4)
     ny = channels%hpphpp_quantum_numbers(channel*4-1)
     nz = channels%hpphpp_quantum_numbers(channel*4-2)
     tz = channels%hpphpp_quantum_numbers(channel*4-3)
     
     local_ch_id = channel 
     !
     ! ket side if fully stored on each proc
     !
     bra_min = mapping_hpphpp(iam+1,local_ch_id,2)
     bra_max = mapping_hpphpp(iam+1,local_ch_id,3)
     !
     ! bra side is distributed 
     !
     ket_min = mapping_hpphpp(iam+1,local_ch_id,4)
     ket_max = mapping_hpphpp(iam+1,local_ch_id,5)
     
     !$omp parallel default(shared) private(bra,a,b,c,ket,ket2,i,j,k)
     !$omp do schedule(dynamic)
     do bra = bra_min, bra_max 
        a = lookup_hpphpp_configs(1,channel)%ival(1,bra)
        b = lookup_hpphpp_configs(1,channel)%ival(2,bra) 
        c = lookup_hpphpp_configs(1,channel)%ival(3,bra) 
        
        ket2 = 0 
        do ket = 1, size(  lookup_hpphpp_configs(1,channel)%ival, 2) 
           i = lookup_hpphpp_configs(1,channel)%ival(1,ket) 
           j = lookup_hpphpp_configs(1,channel)%ival(2,ket) 
           k = lookup_hpphpp_configs(1,channel)%ival(3,ket) 
           
           if ( j >= k ) cycle 
           ket2 = ket2 + 1 
           v3nf_hpphpp(channel)%val(bra,ket2) = chiral_3nf_asym(a,b,c,i,j,k)
        end do
     end do
     !$omp end do
     !$omp end parallel
     
  end do
  deallocate( numchannels1 ) 

  if ( iam == 0 ) write(6,*) 'ok2'
  call mpi_barrier(mpi_comm_world, ierror) 

  !
  ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !
  ! setup <hhp||hhp> channels and procs mapping
  !
  numchannels = 0 
  !     loop over isospin projection
  DO tz=tz_min3,tz_max3,2 
     !     loop over total momentum 
     DO nx=nx_min3, nx_max3 
        !     loop over total momentum 
        DO ny=ny_min3, ny_max3 
           !     loop over total momentum 
           DO nz=nz_min3, nz_max3 
              
              if ( numstates(2,tz, nx, ny, nz) /= 0 ) then
                 numchannels = numchannels + 1 
                 locate_v3channel(2,tz, nx, ny, nz) = numchannels
              end if
              
           end DO
        end DO
     end DO
  end DO
  
  
  mem = 0.d0 
  number_channels = 0 
  number_channels = maxval(locate_v3channel(2,:,:,:,:) ) 
  channels%number_hhphhp_confs  = number_channels
  ALLOCATE( channels%hhphhp_quantum_numbers(4*number_channels) )
  mem = mem + dble(number_channels)*4.*4./1.e9 
  
  if ( iam == 0 ) write(6,*) 'Total memory for t3-channels', mem, 'Gbyte' 
  
  !     loop over isospin projection
  DO tz=tz_min3,tz_max3,2 
     !     loop over total momentum 
     DO nx=nx_min3, nx_max3 
        !     loop over total momentum 
        DO ny=ny_min3, ny_max3 
           !     loop over total momentum 
           DO nz=nz_min3, nz_max3 
                 
              number_channels = locate_v3channel(2,tz,nx,ny,nz) 
              if ( number_channels == 0 ) cycle 
              
              k1 = number_channels*4 
              k2 = number_channels*4 - 1
              k3 = number_channels*4 - 2
              k4 = number_channels*4 - 3
              channels%hhphhp_quantum_numbers(k1) = nx
              channels%hhphhp_quantum_numbers(k2) = ny
              channels%hhphhp_quantum_numbers(k3) = nz
              channels%hhphhp_quantum_numbers(k4) = tz
                            
           end DO
        end DO
     end DO
  end DO

  ! 
  !  setup processor mapping and 
  !  calculate work for hhhppp blocks 
  !
  ltot_work=0
  tot_work = 0
  do channel   = 1, channels%number_hhphhp_confs 
     nx = channels%hhphhp_quantum_numbers(channel*4)
     ny = channels%hhphhp_quantum_numbers(channel*4-1)
     nz = channels%hhphhp_quantum_numbers(channel*4-2)
     tz = channels%hhphhp_quantum_numbers(channel*4-3)

     dim1 = numstates(2,tz, nx, ny, nz) 
     dim2 = numstates(2,tz, nx, ny, nz) 
     
     tot_work = tot_work + dim1*dim2
     ltot_work = ltot_work + int(dim1,8)*int(dim2, 8)
     
  end do
  
  number_channels = channels%number_hhphhp_confs 
  
  work = 0.d0
  memory = 0.d0
  
  work_pr_proc = int( ltot_work/int(num_procs,8),8 )
  work = 0.d0
  
  allocate( mapping_hhphhp(1:num_procs, 1:number_channels, 1:5) )
  mapping_hhphhp = 0 
  if ( iam == 0 ) write(6,*) ltot_work, work_pr_proc, number_channels 

  curr_work = 0
  curr_proc = 0
  local_ch_id = 1
  number_channels = 0
  do channel   = 1, channels%number_hhphhp_confs 
     nx = channels%hhphhp_quantum_numbers(channel*4)
     ny = channels%hhphhp_quantum_numbers(channel*4-1)
     nz = channels%hhphhp_quantum_numbers(channel*4-2)
     tz = channels%hhphhp_quantum_numbers(channel*4-3)
     
     dim1 = numstates(2,tz, nx, ny, nz) 
     dim2 = numstates(2,tz, nx, ny, nz) 
     
     !dim1 = size(  lookup_t3_configs(1,channel)%ival, 2) 
     !dim2 = size(  lookup_t3_configs(2,channel)%ival, 2) 
     
     bra_dim = dim1 
     ket_dim = dim2
     
     number_channels = number_channels + 1 
     has_ch_been_added = .FALSE. 
           
     i = 0
     DO WHILE ( i < ket_dim )
              
        i = i + 1
        temp = curr_work + bra_dim
              
        if ( temp <= work_pr_proc .or. curr_proc == num_procs - 1 ) then
           curr_work = temp
           if (  has_ch_been_added ) then
                    
              mapping_hhphhp(curr_proc+1,local_ch_id,3 ) = i 
              
           else
              
              
              has_ch_been_added = .true. 
              mapping_hhphhp(curr_proc+1,local_ch_id,1 ) = number_channels 
              mapping_hhphhp(curr_proc+1,local_ch_id,2 ) = i
              mapping_hhphhp(curr_proc+1,local_ch_id,3 ) = i
              mapping_hhphhp(curr_proc+1,local_ch_id,4 ) = 1
              mapping_hhphhp(curr_proc+1,local_ch_id,5 ) = bra_dim
              
           end if

           if ( i == ket_dim ) then
              
              local_ch_id = local_ch_id + 1
              
           end if
           work(curr_proc+1) = curr_work
           
        else 
                 
           has_ch_been_added = .false.
           work(curr_proc+1) = curr_work
           curr_work = 0
           curr_proc = curr_proc + 1
           i = i - 1
           
        end if
     end DO
  end DO
  
  number_channels = channels%number_hhphhp_confs 
  work = 0.d0
  do i = 1, num_procs
     do local_ch_id = 1, number_channels 
        if ( mapping_hhphhp(i,local_ch_id,1) == 0 ) cycle
        
        row_dim = mapping_hhphhp(i,local_ch_id,5)-mapping_hhphhp(i,local_ch_id,4)+1
        column_dim = mapping_hhphhp(i,local_ch_id,3)-mapping_hhphhp(i,local_ch_id,2)+1
        work(i) = work(i) + real( row_dim*column_dim )
        
     end do
  end do
  
  call mpi_barrier(mpi_comm_world,ierror)
  number_channels = channels%number_hhphhp_confs 
  
  allocate( my_hhphhp_channel_low(0:num_procs-1), my_hhphhp_channel_high(0:num_procs-1) )
  allocate( check_my_channel_hhphhp(1:number_channels) )
  check_my_channel_hhphhp = 0
  
  my_hhphhp_channel_low(iam)  = number_channels + 1
  my_hhphhp_channel_high(iam) = 0 
  do local_ch_id = 1, number_channels 
     if ( mapping_hhphhp(iam+1,local_ch_id,1) == 0 ) cycle
     check_my_channel_hhphhp(local_ch_id) = 1
      
     if ( local_ch_id > my_hhphhp_channel_high(iam) ) my_hhphhp_channel_high(iam) = local_ch_id
     if ( local_ch_id < my_hhphhp_channel_low(iam) ) my_hhphhp_channel_low(iam) = local_ch_id
     
  end do
  
  call mpi_barrier(mpi_comm_world, ierror)
  allocate( v3nf_hhphhp(my_hhphhp_channel_low(iam):my_hhphhp_channel_high(iam)) )
  call mpi_barrier(mpi_comm_world, ierror)
  
  

  
  ALLOCATE( lookup_hhphhp_configs(1:1,my_hhphhp_channel_low(iam):my_hhphhp_channel_high(iam)) )
  mem = 0.d0  
  mem_3nf = 0.d0 
  do channel = my_hhphhp_channel_low(iam), my_hhphhp_channel_high(iam)
     nx = channels%hhphhp_quantum_numbers(channel*4)
     ny = channels%hhphhp_quantum_numbers(channel*4-1)
     nz = channels%hhphhp_quantum_numbers(channel*4-2)
     tz = channels%hhphhp_quantum_numbers(channel*4-3)
     
     local_ch_id = channel 
     !
     ! ket side if fully stored on each proc
     !
     bra_min = mapping_hhphhp(iam+1,local_ch_id,2)
     bra_max = mapping_hhphhp(iam+1,local_ch_id,3)
     !
     ! bra side is distributed 
     !
     ket_min = mapping_hhphhp(iam+1,local_ch_id,4)
     ket_max = mapping_hhphhp(iam+1,local_ch_id,5)
     
     ij_confs = numstates(2,tz, nx, ny, nz) 
     ab_confs = numstates(2,tz, nx, ny, nz) 
     mem = mem + dble(ij_confs+ab_confs)*3.*4./1.e9 
     mem_3nf = mem_3nf + dble(ij_confs*ab_confs)*16./1.e9 
     ALLOCATE( lookup_hhphhp_configs(1,channel)%ival(3,ij_confs))
     !ALLOCATE( lookup_hhphhp_configs(2,channel)%ival(3,ab_confs))
     ALLOCATE( v3nf_hhphhp(channel)%val(bra_min:bra_max, ket_min:ket_max) )
     v3nf_hhphhp(channel)%val = 0.d0 
       
  end DO
  if ( iam == 0 ) write(6,*) 'Total memory for lookup hhphhp-confs', mem, 'Gbyte'
  if ( iam == 0 ) write(6,*) 'Total memory for v3nf_hhphhp', mem_3nf, 'Gbyte'
  
  
  allocate( numchannels1( maxval( locate_v3channel(2,:,:,:,:) ) ) )
  !
  ! hhp configurations
  !
  numchannels1 = 0
  do a = 1, below_ef
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
        do c = below_ef+1, tot_orbs
           nxc = all_orbit%nx(c)
           nyc = all_orbit%ny(c)
           nzc = all_orbit%nz(c)
           tzc = all_orbit%itzp(c)
           szc = all_orbit%szp(c)
           
           
           tz = (tza + tzb + tzc)
           sz = (sza + szb + szc)
           Nx = nxa + nxb + nxc
           ny = nya + nyb + nyc
           nz = nza + nzb + nzc
           
           
           channel = locate_v3channel(2,tz, nx, ny, nz) 
           if ( channel == 0 ) cycle 
           if ( check_my_channel_hhphhp(channel) == 0 ) cycle 
           
           numchannels1(channel) = numchannels1(channel)+1
           ichan = numchannels1(channel)
           
           lookup_hhphhp_configs(1,channel)%ival(1,ichan) = a 
           lookup_hhphhp_configs(1,channel)%ival(2,ichan) = b 
           lookup_hhphhp_configs(1,channel)%ival(3,ichan) = c 
           
        end do
     end do
  end do

  deallocate( numchannels1 )
  
  do channel = my_hhphhp_channel_low(iam), my_hhphhp_channel_high(iam)
     nx = channels%hhphhp_quantum_numbers(channel*4)
     ny = channels%hhphhp_quantum_numbers(channel*4-1)
     nz = channels%hhphhp_quantum_numbers(channel*4-2)
     tz = channels%hhphhp_quantum_numbers(channel*4-3)
     
     local_ch_id = channel 
     !
     ! ket side if fully stored on each proc
     !
     bra_min = mapping_hhphhp(iam+1,local_ch_id,2)
     bra_max = mapping_hhphhp(iam+1,local_ch_id,3)
     !
     ! bra side is distributed 
     !
     ket_min = mapping_hhphhp(iam+1,local_ch_id,4)
     ket_max = mapping_hhphhp(iam+1,local_ch_id,5)
     
     !$omp parallel default(shared) private(bra,a,b,c,ket,i,j,k)
     !$omp do schedule(dynamic)
     do bra = bra_min, bra_max 
        a = lookup_hhphhp_configs(1,channel)%ival(1,bra)
        b = lookup_hhphhp_configs(1,channel)%ival(2,bra) 
        c = lookup_hhphhp_configs(1,channel)%ival(3,bra) 
        
        if ( b > a ) cycle 
        do ket = 1, size(  lookup_hhphhp_configs(1,channel)%ival, 2) 
           i = lookup_hhphhp_configs(1,channel)%ival(1,ket) 
           j = lookup_hhphhp_configs(1,channel)%ival(2,ket) 
           k = lookup_hhphhp_configs(1,channel)%ival(3,ket) 
           
           v3nf_hhphhp(channel)%val(bra,ket) = chiral_3nf_asym(a,b,c,i,j,k)
        end do
     end do
     !$omp end do
     !$omp end parallel
     
  end do


  !
  ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !
  ! setup <hhp||ppp> channels and procs mapping
  !
  numchannels = 0 
  !     loop over isospin projection
  DO tz=tz_min3,tz_max3,2 
     !     loop over total momentum 
     DO nx=nx_min3, nx_max3 
        !     loop over total momentum 
        DO ny=ny_min3, ny_max3 
           !     loop over total momentum 
           DO nz=nz_min3, nz_max3 
              
              if ( numstates(2,tz, nx, ny, nz)*numstates(4,tz, nx, ny, nz) /= 0 ) then
                 numchannels = numchannels + 1 
                 locate_v3channel(3,tz, nx, ny, nz) = numchannels
              end if
              
           end DO
        end DO
     end DO
  end DO
  
  
  mem = 0.d0 
  number_channels = 0 
  number_channels = maxval(locate_v3channel(3,:,:,:,:) ) 
  channels%number_hhpppp_confs  = number_channels
  ALLOCATE( channels%hhpppp_quantum_numbers(4*number_channels) )
  mem = mem + dble(number_channels)*4.*4./1.e9 
  
  if ( iam == 0 ) write(6,*) 'Total memory for t3-channels', mem, 'Gbyte' 
  
  !     loop over isospin projection
  DO tz=tz_min3,tz_max3,2 
     !     loop over total momentum 
     DO nx=nx_min3, nx_max3 
        !     loop over total momentum 
        DO ny=ny_min3, ny_max3 
           !     loop over total momentum 
           DO nz=nz_min3, nz_max3 
                 
              number_channels = locate_v3channel(3,tz,nx,ny,nz) 
              if ( number_channels == 0 ) cycle 
              
              k1 = number_channels*4 
              k2 = number_channels*4 - 1
              k3 = number_channels*4 - 2
              k4 = number_channels*4 - 3
              channels%hhpppp_quantum_numbers(k1) = nx
              channels%hhpppp_quantum_numbers(k2) = ny
              channels%hhpppp_quantum_numbers(k3) = nz
              channels%hhpppp_quantum_numbers(k4) = tz
                            
           end DO
        end DO
     end DO
  end DO

  ! 
  !  setup processor mapping and 
  !  calculate work for hhhppp blocks 
  !
  ltot_work=0
  tot_work = 0
  do channel   = 1, channels%number_hhpppp_confs 
     nx = channels%hhpppp_quantum_numbers(channel*4)
     ny = channels%hhpppp_quantum_numbers(channel*4-1)
     nz = channels%hhpppp_quantum_numbers(channel*4-2)
     tz = channels%hhpppp_quantum_numbers(channel*4-3)

     dim1 = numstates(2,tz, nx, ny, nz)/2
     dim2 = numstates(4,tz, nx, ny, nz) 
     
     tot_work = tot_work + dim1*dim2
     ltot_work = ltot_work + int(dim1,8)*int(dim2, 8)
     
  end do
  
  number_channels = channels%number_hhpppp_confs 
  
  work = 0.d0
  memory = 0.d0
  
  work_pr_proc = int( ltot_work/int(num_procs,8),8 )
  work = 0.d0
  
  allocate( mapping_hhpppp(1:num_procs, 1:number_channels, 1:5) )
  mapping_hhpppp = 0 
  if ( iam == 0 ) write(6,*) ltot_work, work_pr_proc, number_channels 

  curr_work = 0
  curr_proc = 0
  local_ch_id = 1
  number_channels = 0
  do channel   = 1, channels%number_hhpppp_confs 
     nx = channels%hhpppp_quantum_numbers(channel*4)
     ny = channels%hhpppp_quantum_numbers(channel*4-1)
     nz = channels%hhpppp_quantum_numbers(channel*4-2)
     tz = channels%hhpppp_quantum_numbers(channel*4-3)
     
     dim1 = numstates(2,tz, nx, ny, nz) 
     dim2 = numstates(4,tz, nx, ny, nz) 
     
     !dim1 = size(  lookup_t3_configs(1,channel)%ival, 2) 
     !dim2 = size(  lookup_t3_configs(2,channel)%ival, 2) 
     
     bra_dim = dim1/2 
     ket_dim = dim2
     
     number_channels = number_channels + 1 
     has_ch_been_added = .FALSE. 
           
     i = 0
     DO WHILE ( i < ket_dim )
              
        i = i + 1
        temp = curr_work + bra_dim
              
        if ( temp <= work_pr_proc .or. curr_proc == num_procs - 1 ) then
           curr_work = temp
           if (  has_ch_been_added ) then
                    
              mapping_hhpppp(curr_proc+1,local_ch_id,3 ) = i 
              
           else
              
              
              has_ch_been_added = .true. 
              mapping_hhpppp(curr_proc+1,local_ch_id,1 ) = number_channels 
              mapping_hhpppp(curr_proc+1,local_ch_id,2 ) = i
              mapping_hhpppp(curr_proc+1,local_ch_id,3 ) = i
              mapping_hhpppp(curr_proc+1,local_ch_id,4 ) = 1
              mapping_hhpppp(curr_proc+1,local_ch_id,5 ) = bra_dim
              
           end if

           if ( i == ket_dim ) then
              
              local_ch_id = local_ch_id + 1
              
           end if
           work(curr_proc+1) = curr_work
           
        else 
                 
           has_ch_been_added = .false.
           work(curr_proc+1) = curr_work
           curr_work = 0
           curr_proc = curr_proc + 1
           i = i - 1
           
        end if
     end DO
  end DO
  
  number_channels = channels%number_hhpppp_confs 
  work = 0.d0
  do i = 1, num_procs
     do local_ch_id = 1, number_channels 
        if ( mapping_hhpppp(i,local_ch_id,1) == 0 ) cycle
        
        row_dim = mapping_hhpppp(i,local_ch_id,5)-mapping_hhpppp(i,local_ch_id,4)+1
        column_dim = mapping_hhpppp(i,local_ch_id,3)-mapping_hhpppp(i,local_ch_id,2)+1
        work(i) = work(i) + real( row_dim*column_dim )
        
     end do
  end do
  
  call mpi_barrier(mpi_comm_world,ierror)
  number_channels = channels%number_hhpppp_confs 
  
  allocate( my_hhpppp_channel_low(0:num_procs-1), my_hhpppp_channel_high(0:num_procs-1) )
  allocate( check_my_channel_hhpppp(1:number_channels) )
  check_my_channel_hhpppp = 0
  
  my_hhpppp_channel_low(iam)  = number_channels + 1
  my_hhpppp_channel_high(iam) = 0 
  do local_ch_id = 1, number_channels 
     if ( mapping_hhpppp(iam+1,local_ch_id,1) == 0 ) cycle
     check_my_channel_hhpppp(local_ch_id) = 1
      
     if ( local_ch_id > my_hhpppp_channel_high(iam) ) my_hhpppp_channel_high(iam) = local_ch_id
     if ( local_ch_id < my_hhpppp_channel_low(iam) ) my_hhpppp_channel_low(iam) = local_ch_id
     
  end do
  
  call mpi_barrier(mpi_comm_world, ierror)
  allocate( v3nf_hhpppp(my_hhpppp_channel_low(iam):my_hhpppp_channel_high(iam)) )
  call mpi_barrier(mpi_comm_world, ierror)
  
  

  
  ALLOCATE( lookup_hhpppp_configs(1:2,my_hhpppp_channel_low(iam):my_hhpppp_channel_high(iam)) )
  mem = 0.d0  
  mem_3nf = 0.d0 
  do channel = my_hhpppp_channel_low(iam), my_hhpppp_channel_high(iam)
     nx = channels%hhpppp_quantum_numbers(channel*4)
     ny = channels%hhpppp_quantum_numbers(channel*4-1)
     nz = channels%hhpppp_quantum_numbers(channel*4-2)
     tz = channels%hhpppp_quantum_numbers(channel*4-3)
     
     local_ch_id = channel 
     !
     ! ket side if fully stored on each proc
     !
     bra_min = mapping_hhpppp(iam+1,local_ch_id,2)
     bra_max = mapping_hhpppp(iam+1,local_ch_id,3)
     !
     ! bra side is distributed 
     !
     ket_min = mapping_hhpppp(iam+1,local_ch_id,4)
     ket_max = mapping_hhpppp(iam+1,local_ch_id,5)
     
     ij_confs = numstates(2,tz, nx, ny, nz)/2 
     ab_confs = numstates(4,tz, nx, ny, nz) 
     mem = mem + dble(ij_confs+ab_confs)*3.*4./1.e9 
     mem_3nf = mem_3nf + dble(ij_confs*ab_confs)*16./1.e9 
     ALLOCATE( lookup_hhpppp_configs(1,channel)%ival(3,ij_confs))
     ALLOCATE( lookup_hhpppp_configs(2,channel)%ival(3,ab_confs))
     ALLOCATE( v3nf_hhpppp(channel)%val(ket_min:ket_max,bra_min:bra_max) )
     v3nf_hhpppp(channel)%val = 0.d0 
     
  end DO
  if ( iam == 0 ) write(6,*) 'Total memory for lookup hhpppp-confs', mem, 'Gbyte'
  if ( iam == 0 ) write(6,*) 'Total memory for v3nf_hhpppp', mem_3nf, 'Gbyte'
  
  
  allocate( numchannels1( maxval( locate_v3channel(3,:,:,:,:) ) ) )
  !
  ! hhp configurations
  !
  numchannels1 = 0
  do a = 1, below_ef
     nxa = all_orbit%nx(a)
     nya = all_orbit%ny(a)
     nza = all_orbit%nz(a)
     tza=all_orbit%itzp(a)
     sza=all_orbit%szp(a)
     
!     do b = 1, below_ef
     do b = a+1, below_ef
        nxb = all_orbit%nx(b)
        nyb = all_orbit%ny(b)
        nzb = all_orbit%nz(b)
        tzb = all_orbit%itzp(b)
        szb = all_orbit%szp(b)
        do c = below_ef+1, tot_orbs
           nxc = all_orbit%nx(c)
           nyc = all_orbit%ny(c)
           nzc = all_orbit%nz(c)
           tzc = all_orbit%itzp(c)
           szc = all_orbit%szp(c)
           
           
           tz = (tza + tzb + tzc)
           sz = (sza + szb + szc)
           Nx = nxa + nxb + nxc
           ny = nya + nyb + nyc
           nz = nza + nzb + nzc
           
           
           channel = locate_v3channel(3,tz, nx, ny, nz) 
           if ( channel == 0 ) cycle 
           if ( check_my_channel_hhpppp(channel) == 0 ) cycle 
           
           numchannels1(channel) = numchannels1(channel)+1
           ichan = numchannels1(channel)
           
           lookup_hhpppp_configs(1,channel)%ival(1,ichan) = a 
           lookup_hhpppp_configs(1,channel)%ival(2,ichan) = b 
           lookup_hhpppp_configs(1,channel)%ival(3,ichan) = c 
           
        end do
     end do
  end do

  !
  ! ppp configurations
  !
  numchannels1 = 0
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
        do c = below_ef+1, tot_orbs
           nxc = all_orbit%nx(c)
           nyc = all_orbit%ny(c)
           nzc = all_orbit%nz(c)
           tzc = all_orbit%itzp(c)
           szc = all_orbit%szp(c)
           
           
           tz = (tza + tzb + tzc)
           sz = (sza + szb + szc)
           Nx = nxa + nxb + nxc
           ny = nya + nyb + nyc
           nz = nza + nzb + nzc
           
           
           channel = locate_v3channel(3,tz, nx, ny, nz) 
           if ( channel == 0 ) cycle 
           if ( check_my_channel_hhpppp(channel) == 0 ) cycle 
           
           numchannels1(channel) = numchannels1(channel)+1
           ichan = numchannels1(channel)
           
           lookup_hhpppp_configs(2,channel)%ival(1,ichan) = a 
           lookup_hhpppp_configs(2,channel)%ival(2,ichan) = b 
           lookup_hhpppp_configs(2,channel)%ival(3,ichan) = c 
           
        end do
     end do
  end do

  deallocate( numchannels1 )
  
  do channel = my_hhpppp_channel_low(iam), my_hhpppp_channel_high(iam)
     nx = channels%hhpppp_quantum_numbers(channel*4)
     ny = channels%hhpppp_quantum_numbers(channel*4-1)
     nz = channels%hhpppp_quantum_numbers(channel*4-2)
     tz = channels%hhpppp_quantum_numbers(channel*4-3)
     
     local_ch_id = channel 
     !
     ! ket side if fully stored on each proc
     !
     bra_min = mapping_hhpppp(iam+1,local_ch_id,2)
     bra_max = mapping_hhpppp(iam+1,local_ch_id,3)
     !
     ! bra side is distributed 
     !
     ket_min = mapping_hhpppp(iam+1,local_ch_id,4)
     ket_max = mapping_hhpppp(iam+1,local_ch_id,5)
     
     !$omp parallel default(shared) private(bra,a,b,c,ket,i,j,k)
     !$omp do schedule(dynamic)
     do ket = bra_min, bra_max 
        a = lookup_hhpppp_configs(2,channel)%ival(1,ket)
        b = lookup_hhpppp_configs(2,channel)%ival(2,ket) 
        c = lookup_hhpppp_configs(2,channel)%ival(3,ket) 
        
        do bra = 1, size(  lookup_hhpppp_configs(1,channel)%ival, 2) 
           i = lookup_hhpppp_configs(1,channel)%ival(1,bra) 
           j = lookup_hhpppp_configs(1,channel)%ival(2,bra) 
           k = lookup_hhpppp_configs(1,channel)%ival(3,bra) 
           
           v3nf_hhpppp(channel)%val(bra,ket) = chiral_3nf_asym(i,j,k,a,b,c)
        end do
     end do
     !$omp end do
     !$omp end parallel
     
  end do
  
  !
  ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !
  ! setup <hhh||hpp> channels and procs mapping
  !
  numchannels = 0 
  !     loop over isospin projection
  DO tz=tz_min3,tz_max3,2 
     !     loop over total momentum 
     DO nx=nx_min3, nx_max3 
        !     loop over total momentum 
        DO ny=ny_min3, ny_max3 
           !     loop over total momentum 
           DO nz=nz_min3, nz_max3 
              
              if ( numstates(1,tz, nx, ny, nz)*numstates(3,tz, nx, ny, nz) /= 0 ) then
                 numchannels = numchannels + 1 
                 locate_v3channel(4,tz, nx, ny, nz) = numchannels
              end if
              
           end DO
        end DO
     end DO
  end DO
  
  
  mem = 0.d0 
  number_channels = 0 
  number_channels = maxval(locate_v3channel(4,:,:,:,:) ) 
  channels%number_hhhhpp_confs  = number_channels
  ALLOCATE( channels%hhhhpp_quantum_numbers(4*number_channels) )
  mem = mem + dble(number_channels)*4.*4./1.e9 
  
  if ( iam == 0 ) write(6,*) 'Total memory for t3-channels', mem, 'Gbyte' 
  
  !     loop over isospin projection
  DO tz=tz_min3,tz_max3,2 
     !     loop over total momentum 
     DO nx=nx_min3, nx_max3 
        !     loop over total momentum 
        DO ny=ny_min3, ny_max3 
           !     loop over total momentum 
           DO nz=nz_min3, nz_max3 
                 
              number_channels = locate_v3channel(4,tz,nx,ny,nz) 
              if ( number_channels == 0 ) cycle 
              
              k1 = number_channels*4 
              k2 = number_channels*4 - 1
              k3 = number_channels*4 - 2
              k4 = number_channels*4 - 3
              channels%hhhhpp_quantum_numbers(k1) = nx
              channels%hhhhpp_quantum_numbers(k2) = ny
              channels%hhhhpp_quantum_numbers(k3) = nz
              channels%hhhhpp_quantum_numbers(k4) = tz
                            
           end DO
        end DO
     end DO
  end DO

  ! 
  !  setup processor mapping and 
  !  calculate work for hhhppp blocks 
  !
  ltot_work=0
  tot_work = 0
  do channel   = 1, channels%number_hhhhpp_confs 
     nx = channels%hhhhpp_quantum_numbers(channel*4)
     ny = channels%hhhhpp_quantum_numbers(channel*4-1)
     nz = channels%hhhhpp_quantum_numbers(channel*4-2)
     tz = channels%hhhhpp_quantum_numbers(channel*4-3)

     dim1 = numstates(1,tz, nx, ny, nz) 
     dim2 = numstates(3,tz, nx, ny, nz) 
     
     tot_work = tot_work + dim1*dim2
     ltot_work = ltot_work + int(dim1,8)*int(dim2, 8)
     
  end do
  
  number_channels = channels%number_hhhhpp_confs 
  
  work = 0.d0
  memory = 0.d0
  
  work_pr_proc = int( ltot_work/int(num_procs,8),8 )
  work = 0.d0
  
  allocate( mapping_hhhhpp(1:num_procs, 1:number_channels, 1:5) )
  mapping_hhhhpp = 0 
  if ( iam == 0 ) write(6,*) ltot_work, work_pr_proc, number_channels 

  curr_work = 0
  curr_proc = 0
  local_ch_id = 1
  number_channels = 0
  do channel   = 1, channels%number_hhhhpp_confs 
     nx = channels%hhhhpp_quantum_numbers(channel*4)
     ny = channels%hhhhpp_quantum_numbers(channel*4-1)
     nz = channels%hhhhpp_quantum_numbers(channel*4-2)
     tz = channels%hhhhpp_quantum_numbers(channel*4-3)
     
     dim1 = numstates(1,tz, nx, ny, nz) 
     dim2 = numstates(3,tz, nx, ny, nz) 
     
     !dim1 = size(  lookup_t3_configs(1,channel)%ival, 2) 
     !dim2 = size(  lookup_t3_configs(2,channel)%ival, 2) 
     
     bra_dim = dim1 
     ket_dim = dim2
     
     number_channels = number_channels + 1 
     has_ch_been_added = .FALSE. 
           
     i = 0
     DO WHILE ( i < ket_dim )
              
        i = i + 1
        temp = curr_work + bra_dim
              
        if ( temp <= work_pr_proc .or. curr_proc == num_procs - 1 ) then
           curr_work = temp
           if (  has_ch_been_added ) then
                    
              mapping_hhhhpp(curr_proc+1,local_ch_id,3 ) = i 
              
           else
              
              
              has_ch_been_added = .true. 
              mapping_hhhhpp(curr_proc+1,local_ch_id,1 ) = number_channels 
              mapping_hhhhpp(curr_proc+1,local_ch_id,2 ) = i
              mapping_hhhhpp(curr_proc+1,local_ch_id,3 ) = i
              mapping_hhhhpp(curr_proc+1,local_ch_id,4 ) = 1
              mapping_hhhhpp(curr_proc+1,local_ch_id,5 ) = bra_dim
              
           end if

           if ( i == ket_dim ) then
              
              local_ch_id = local_ch_id + 1
              
           end if
           work(curr_proc+1) = curr_work
           
        else 
                 
           has_ch_been_added = .false.
           work(curr_proc+1) = curr_work
           curr_work = 0
           curr_proc = curr_proc + 1
           i = i - 1
           
        end if
     end DO
  end DO
  
  number_channels = channels%number_hhhhpp_confs 
  work = 0.d0
  do i = 1, num_procs
     do local_ch_id = 1, number_channels 
        if ( mapping_hhhhpp(i,local_ch_id,1) == 0 ) cycle
        
        row_dim = mapping_hhhhpp(i,local_ch_id,5)-mapping_hhhhpp(i,local_ch_id,4)+1
        column_dim = mapping_hhhhpp(i,local_ch_id,3)-mapping_hhhhpp(i,local_ch_id,2)+1
        work(i) = work(i) + real( row_dim*column_dim )
        
     end do
  end do
  
  call mpi_barrier(mpi_comm_world,ierror)
  number_channels = channels%number_hhhhpp_confs 
  
  allocate( my_hhhhpp_channel_low(0:num_procs-1), my_hhhhpp_channel_high(0:num_procs-1) )
  allocate( check_my_channel_hhhhpp(1:number_channels) )
  check_my_channel_hhhhpp = 0
  
  my_hhhhpp_channel_low(iam)  = number_channels + 1
  my_hhhhpp_channel_high(iam) = 0 
  do local_ch_id = 1, number_channels 
     if ( mapping_hhhhpp(iam+1,local_ch_id,1) == 0 ) cycle
     check_my_channel_hhhhpp(local_ch_id) = 1
      
     if ( local_ch_id > my_hhhhpp_channel_high(iam) ) my_hhhhpp_channel_high(iam) = local_ch_id
     if ( local_ch_id < my_hhhhpp_channel_low(iam) ) my_hhhhpp_channel_low(iam) = local_ch_id
     
  end do
  
  call mpi_barrier(mpi_comm_world, ierror)
  allocate( v3nf_hhhhpp(my_hhhhpp_channel_low(iam):my_hhhhpp_channel_high(iam)) )
  call mpi_barrier(mpi_comm_world, ierror)
  
  

  
  ALLOCATE( lookup_hhhhpp_configs(1:2,my_hhhhpp_channel_low(iam):my_hhhhpp_channel_high(iam)) )
  mem = 0.d0  
  mem_3nf = 0.d0 
  do channel = my_hhhhpp_channel_low(iam), my_hhhhpp_channel_high(iam)
     nx = channels%hhhhpp_quantum_numbers(channel*4)
     ny = channels%hhhhpp_quantum_numbers(channel*4-1)
     nz = channels%hhhhpp_quantum_numbers(channel*4-2)
     tz = channels%hhhhpp_quantum_numbers(channel*4-3)
     
     local_ch_id = channel 
     !
     ! ket side if fully stored on each proc
     !
     bra_min = mapping_hhhhpp(iam+1,local_ch_id,2)
     bra_max = mapping_hhhhpp(iam+1,local_ch_id,3)
     !
     ! bra side is distributed 
     !
     ket_min = mapping_hhhhpp(iam+1,local_ch_id,4)
     ket_max = mapping_hhhhpp(iam+1,local_ch_id,5)
     
     ij_confs = numstates(1,tz, nx, ny, nz) 
     ab_confs = numstates(3,tz, nx, ny, nz) 
     mem = mem + dble(ij_confs+ab_confs)*3.*4./1.e9 
     mem_3nf = mem_3nf + dble(ij_confs*ab_confs)*16./1.e9 
     ALLOCATE( lookup_hhhhpp_configs(1,channel)%ival(3,ij_confs))
     ALLOCATE( lookup_hhhhpp_configs(2,channel)%ival(3,ab_confs))
     ALLOCATE( v3nf_hhhhpp(channel)%val(ket_min:ket_max,bra_min:bra_max) )
     v3nf_hhhhpp(channel)%val = 0.d0 
     
  end DO
  if ( iam == 0 ) write(6,*) 'Total memory for lookup hhhhpp-confs', mem, 'Gbyte'
  if ( iam == 0 ) write(6,*) 'Total memory for v3nf_hhhhpp', mem_3nf, 'Gbyte'
  
  
  allocate( numchannels1( maxval( locate_v3channel(4,:,:,:,:) ) ) )
  !
  ! hhh configurations
  !
  numchannels1 = 0
  do a = 1, below_ef
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
        do c = 1, below_ef
           nxc = all_orbit%nx(c)
           nyc = all_orbit%ny(c)
           nzc = all_orbit%nz(c)
           tzc = all_orbit%itzp(c)
           szc = all_orbit%szp(c)
           
           
           tz = (tza + tzb + tzc)
           sz = (sza + szb + szc)
           Nx = nxa + nxb + nxc
           ny = nya + nyb + nyc
           nz = nza + nzb + nzc
           
           
           channel = locate_v3channel(4,tz, nx, ny, nz) 
           if ( channel == 0 ) cycle 
           if ( check_my_channel_hhhhpp(channel) == 0 ) cycle 
           
           numchannels1(channel) = numchannels1(channel)+1
           ichan = numchannels1(channel)
           
           lookup_hhhhpp_configs(1,channel)%ival(1,ichan) = a 
           lookup_hhhhpp_configs(1,channel)%ival(2,ichan) = b 
           lookup_hhhhpp_configs(1,channel)%ival(3,ichan) = c 
           
        end do
     end do
  end do

  !
  ! hpp configurations
  !
  numchannels1 = 0
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
        do c = below_ef+1, tot_orbs
           nxc = all_orbit%nx(c)
           nyc = all_orbit%ny(c)
           nzc = all_orbit%nz(c)
           tzc = all_orbit%itzp(c)
           szc = all_orbit%szp(c)
           
           
           tz = (tza + tzb + tzc)
           sz = (sza + szb + szc)
           Nx = nxa + nxb + nxc
           ny = nya + nyb + nyc
           nz = nza + nzb + nzc
           
           
           channel = locate_v3channel(4,tz, nx, ny, nz) 
           if ( channel == 0 ) cycle 
           if ( check_my_channel_hhhhpp(channel) == 0 ) cycle 
           
           numchannels1(channel) = numchannels1(channel)+1
           ichan = numchannels1(channel)
           
           lookup_hhhhpp_configs(2,channel)%ival(1,ichan) = a 
           lookup_hhhhpp_configs(2,channel)%ival(2,ichan) = b 
           lookup_hhhhpp_configs(2,channel)%ival(3,ichan) = c 
           
        end do
     end do
  end do

  deallocate( numchannels1 )
  
  do channel = my_hhhhpp_channel_low(iam), my_hhhhpp_channel_high(iam)
     nx = channels%hhhhpp_quantum_numbers(channel*4)
     ny = channels%hhhhpp_quantum_numbers(channel*4-1)
     nz = channels%hhhhpp_quantum_numbers(channel*4-2)
     tz = channels%hhhhpp_quantum_numbers(channel*4-3)
     
     local_ch_id = channel 
     !
     ! ket side if fully stored on each proc
     !
     bra_min = mapping_hhhhpp(iam+1,local_ch_id,2)
     bra_max = mapping_hhhhpp(iam+1,local_ch_id,3)
     !
     ! bra side is distributed 
     !
     ket_min = mapping_hhhhpp(iam+1,local_ch_id,4)
     ket_max = mapping_hhhhpp(iam+1,local_ch_id,5)
     
     !$omp parallel default(shared) private(bra,a,b,c,ket,i,j,k)
     !$omp do schedule(dynamic)
     do ket = bra_min, bra_max 
        a = lookup_hhhhpp_configs(2,channel)%ival(1,ket)
        b = lookup_hhhhpp_configs(2,channel)%ival(2,ket) 
        c = lookup_hhhhpp_configs(2,channel)%ival(3,ket) 
        
        if ( b > c ) cycle 
        do bra = 1, size(  lookup_hhhhpp_configs(1,channel)%ival, 2) 
           i = lookup_hhhhpp_configs(1,channel)%ival(1,bra) 
           j = lookup_hhhhpp_configs(1,channel)%ival(2,bra) 
           k = lookup_hhhhpp_configs(1,channel)%ival(3,bra) 
           
           v3nf_hhhhpp(channel)%val(bra,ket) = chiral_3nf_asym(i,j,k,a,b,c)
        end do
     end do
     !$omp end do
     !$omp end parallel
     
  end do
  
  
  
  deallocate( work, memory )
  deallocate( numstates ) 
  
end SUBROUTINE setup_v3nf_channel_structures

SUBROUTINE compute_v3nf_memory
  USE single_particle_orbits
  USE constants
  use one_body_operators
  use t2_storage
  use configurations
  USE KSPACE 
  use CHIRAL_POTENTIALS
  use t3_constants
  
  IMPLICIT NONE
  INTEGER :: i,j,k, sumx, sumy, sumz, a,b,c,nxa,nxb,nxc,nya,nyb,nyc,nza,nzb,nzc,tza,tzb,tzc,sza,szb,szc,tz,sz,nx,ny,nz
  INTEGER :: d, nxd, nyd, nzd, szd, tzd
  INTEGER :: numchannels, ij_confs, ab_confs, channel, sz_min, sz_max, tz_min, tz_max, sumtz, sumsz, ichan
  integer(8), allocatable :: numstates(:,:,:,:,:)
  INTEGER, ALLOCATABLE :: numchannels1(:),amin(:), amax(:), bmin(:), bmax(:), cmin(:), cmax(:) 
  real*8 :: mem, mem_3nf,delta 
  INTEGER :: k1,k2,k3,k4,k5, dim1, bra, ket, ket2 
  INTEGER(8) :: total_work, work_per_proc, curr_work, curr_proc, curr_channel
  integer :: processor_work(num_procs), mtx_dim, mtx_dim_work, temp, local_ch_id
  integer(8) :: row_dim, column_dim, ket_dim, bra_dim
  INTEGER :: bra_min, bra_max, ket_min, ket_max, number_channels, dim2, channel2
  LOGICAL :: has_ch_been_added
  REAL*8, allocatable :: work(:), memory(:)
  integer(8):: ltot_work, tot_work, work_pr_proc

  allocate( work(1:num_procs) )
  allocate( memory(1:num_procs) )
  work = 0.d0
  memory = 0.d0 
  
  nx_min3 = 1000
  nx_max3 = -1000
  ny_min3 = 1000
  ny_max3 = -1000
  nz_min3 = 1000
  nz_max3 = -1000
  tz_min3 = 100 
  tz_max3 = -100
  sz_min3 = 100 
  sz_max3 = -100
  do i = 1, tot_orbs 
     do j = 1, tot_orbs 
        do k = 1, tot_orbs 
           
           sumx = all_orbit%nx(i) +  all_orbit%nx(j) +  all_orbit%nx(k) 
           sumy = all_orbit%ny(i) +  all_orbit%ny(j) +  all_orbit%ny(k) 
           sumz = all_orbit%nz(i) +  all_orbit%nz(j) +  all_orbit%nz(k) 
           
           sumtz = all_orbit%itzp(i) +  all_orbit%itzp(j) +  all_orbit%itzp(k) 
           sumsz = all_orbit%szp(i) +  all_orbit%szp(j) +  all_orbit%szp(k) 
           
           
           if ( sumx > nx_max3 ) nx_max3 = sumx 
           if ( sumy > ny_max3 ) ny_max3 = sumy 
           if ( sumz > nz_max3 ) nz_max3 = sumz 
           if ( sumx < nx_min3 ) nx_min3 = sumx 
           if ( sumy < ny_min3 ) ny_min3 = sumy 
           if ( sumz < nz_min3 ) nz_min3 = sumz 
           if ( sumtz > tz_max3 ) tz_max3 = sumtz 
           if ( sumsz > sz_max3 ) sz_max3 = sumsz 
           if ( sumtz < tz_min3 ) tz_min3 = sumtz 
           if ( sumsz < sz_min3 ) sz_min3 = sumsz 
           
           
           
        end do
     end do
  end do
  
  if ( iam == 0 ) write(6,*) nx_max3,ny_max3, nz_max3, tz_min3,tz_max3
  allocate( locate_v3channel(6, tz_min3:tz_max3, nx_min3:nx_max3, ny_min3:ny_max3, nz_min3:nz_max3) )
  locate_v3channel = 0  
  
  allocate( numstates(1:6,tz_min3:tz_max3, nx_min3:nx_max3, ny_min3:ny_max3, nz_min3:nz_max3) )
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
  
        do c = 1, tot_orbs
           nxc = all_orbit%nx(c)
           nyc = all_orbit%ny(c)
           nzc = all_orbit%nz(c)
           tzc = all_orbit%itzp(c)
           szc = all_orbit%szp(c)
           
           tz = (tza + tzb + tzc)
           sz = (sza + szb + szc)
           Nx = nxa + nxb + nxc
           ny = nya + nyb + nyc
           nz = nza + nzb + nzc
             
           ! hhh
           if ( a <= below_ef .and. b <= below_ef .and. c <= below_ef ) & 
                numstates(1,tz, nx, ny, nz) = numstates(1,tz, nx, ny, nz)  + 1 
           ! hhp
           if ( a <= below_ef .and. b <= below_ef .and. c > below_ef  ) & 
                numstates(2,tz, nx, ny, nz) = numstates(2,tz, nx, ny, nz)  + 1 
           ! hpp
           if ( a <= below_ef .and. b > below_ef  .and. c > below_ef  ) & 
                numstates(3,tz, nx, ny, nz) = numstates(3,tz, nx, ny, nz)  + 1 
           ! ppp
           if ( a > below_ef  .and. b >  below_ef .and. c > below_ef  ) & 
                numstates(4,tz, nx, ny, nz) = numstates(4,tz, nx, ny, nz)  + 1 
           
        end do
     end do
  end do
  

  !
  ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !
  ! setup <hpp||hpp> channels and procs mapping
  !
  numchannels = 0 
  !     loop over isospin projection
  DO tz=tz_min3,tz_max3,2 
     !     loop over total momentum 
     DO nx=nx_min3, nx_max3 
        !     loop over total momentum 
        DO ny=ny_min3, ny_max3 
           !     loop over total momentum 
           DO nz=nz_min3, nz_max3 
              
              if ( numstates(3,tz, nx, ny, nz) /= 0 ) then
                 numchannels = numchannels + 1 
                 locate_v3channel(1,tz, nx, ny, nz) = numchannels
              end if
              
           end DO
        end DO
     end DO
  end DO
  
  
  mem = 0.d0 
  number_channels = 0 
  number_channels = maxval(locate_v3channel(1,:,:,:,:) ) 
  channels%number_hpphpp_confs  = number_channels
  ALLOCATE( channels%hpphpp_quantum_numbers(4*number_channels) )
  mem = mem + dble(number_channels)*4.*4./1.e9 
  
  if ( iam == 0 ) write(6,*) 'Total memory for t3-channels', mem, 'Gbyte' 
  
  !     loop over isospin projection
  DO tz=tz_min3,tz_max3,2 
     !     loop over total momentum 
     DO nx=nx_min3, nx_max3 
        !     loop over total momentum 
        DO ny=ny_min3, ny_max3 
           !     loop over total momentum 
           DO nz=nz_min3, nz_max3 
                 
              number_channels = locate_v3channel(1,tz,nx,ny,nz) 
              if ( number_channels == 0 ) cycle 
              
              k1 = number_channels*4 
              k2 = number_channels*4 - 1
              k3 = number_channels*4 - 2
              k4 = number_channels*4 - 3
              channels%hpphpp_quantum_numbers(k1) = nx
              channels%hpphpp_quantum_numbers(k2) = ny
              channels%hpphpp_quantum_numbers(k3) = nz
              channels%hpphpp_quantum_numbers(k4) = tz
                            
           end DO
        end DO
     end DO
  end DO

  ! 
  !  setup processor mapping and 
  !  calculate work for hhhppp blocks 
  !
  ltot_work=0
  tot_work = 0
  do channel   = 1, channels%number_hpphpp_confs 
     nx = channels%hpphpp_quantum_numbers(channel*4)
     ny = channels%hpphpp_quantum_numbers(channel*4-1)
     nz = channels%hpphpp_quantum_numbers(channel*4-2)
     tz = channels%hpphpp_quantum_numbers(channel*4-3)

     dim1 = numstates(3,tz, nx, ny, nz) 
     dim2 = numstates(3,tz, nx, ny, nz)/2
     
     !dim1 = size(  lookup_t3_configs(1,channel)%ival,2 )
     !dim2 = size(  lookup_t3_configs(2,channel)%ival,2 )
     
     
     tot_work = tot_work + dim1*dim2
     ltot_work = ltot_work + int(dim1,8)*int(dim2, 8)
     
  end do
  
  number_channels = channels%number_hpphpp_confs 
  
  work = 0.d0
  memory = 0.d0
  
  work_pr_proc = int( ltot_work/int(num_procs,8),8 )
  work = 0.d0
  
  allocate( mapping_hpphpp(1:num_procs, 1:number_channels, 1:5) )
  mapping_hpphpp = 0 
  if ( iam == 0 ) write(6,*) ltot_work, work_pr_proc, number_channels 

  curr_work = 0
  curr_proc = 0
  local_ch_id = 1
  number_channels = 0
  do channel   = 1, channels%number_hpphpp_confs 
     nx = channels%hpphpp_quantum_numbers(channel*4)
     ny = channels%hpphpp_quantum_numbers(channel*4-1)
     nz = channels%hpphpp_quantum_numbers(channel*4-2)
     tz = channels%hpphpp_quantum_numbers(channel*4-3)
     
     dim1 = numstates(3,tz, nx, ny, nz)
     dim2 = numstates(3,tz, nx, ny, nz) 
     
     !dim1 = size(  lookup_t3_configs(1,channel)%ival, 2) 
     !dim2 = size(  lookup_t3_configs(2,channel)%ival, 2) 
     
     bra_dim = dim1/2 
     ket_dim = dim2
     
     number_channels = number_channels + 1 
     has_ch_been_added = .FALSE. 
           
     i = 0
     DO WHILE ( i < ket_dim )
              
        i = i + 1
        temp = curr_work + bra_dim
              
        if ( temp <= work_pr_proc .or. curr_proc == num_procs - 1 ) then
           curr_work = temp
           if (  has_ch_been_added ) then
                    
              mapping_hpphpp(curr_proc+1,local_ch_id,3 ) = i 
              
           else
              
              
              has_ch_been_added = .true. 
              mapping_hpphpp(curr_proc+1,local_ch_id,1 ) = number_channels 
              mapping_hpphpp(curr_proc+1,local_ch_id,2 ) = i
              mapping_hpphpp(curr_proc+1,local_ch_id,3 ) = i
              mapping_hpphpp(curr_proc+1,local_ch_id,4 ) = 1
              mapping_hpphpp(curr_proc+1,local_ch_id,5 ) = bra_dim
              
           end if

           if ( i == ket_dim ) then
              
              local_ch_id = local_ch_id + 1
              
           end if
           work(curr_proc+1) = curr_work
           
        else 
                 
           has_ch_been_added = .false.
           work(curr_proc+1) = curr_work
           curr_work = 0
           curr_proc = curr_proc + 1
           i = i - 1
           
        end if
     end DO
  end DO
  
  number_channels = channels%number_hpphpp_confs 
  work = 0.d0
  do i = 1, num_procs
     do local_ch_id = 1, number_channels 
        if ( mapping_hpphpp(i,local_ch_id,1) == 0 ) cycle
        
        row_dim = mapping_hpphpp(i,local_ch_id,5)-mapping_hpphpp(i,local_ch_id,4)+1
        column_dim = mapping_hpphpp(i,local_ch_id,3)-mapping_hpphpp(i,local_ch_id,2)+1
        work(i) = work(i) + real( row_dim*column_dim )
        
     end do
  end do
  
  call mpi_barrier(mpi_comm_world,ierror)
  number_channels = channels%number_hpphpp_confs 
  
  allocate( my_hpphpp_channel_low(0:num_procs-1), my_hpphpp_channel_high(0:num_procs-1) )
  allocate( check_my_channel_hpphpp(1:number_channels) )
  check_my_channel_hpphpp = 0
  
  my_hpphpp_channel_low(iam)  = number_channels + 1
  my_hpphpp_channel_high(iam) = 0 
  do local_ch_id = 1, number_channels 
     if ( mapping_hpphpp(iam+1,local_ch_id,1) == 0 ) cycle
     check_my_channel_hpphpp(local_ch_id) = 1
      
     if ( local_ch_id > my_hpphpp_channel_high(iam) ) my_hpphpp_channel_high(iam) = local_ch_id
     if ( local_ch_id < my_hpphpp_channel_low(iam) ) my_hpphpp_channel_low(iam) = local_ch_id
     
  end do
  
  call mpi_barrier(mpi_comm_world, ierror)
   
  mem = 0.d0  
  mem_3nf = 0.d0 
  do channel = my_hpphpp_channel_low(iam), my_hpphpp_channel_high(iam)
     nx = channels%hpphpp_quantum_numbers(channel*4)
     ny = channels%hpphpp_quantum_numbers(channel*4-1)
     nz = channels%hpphpp_quantum_numbers(channel*4-2)
     tz = channels%hpphpp_quantum_numbers(channel*4-3)
     
     local_ch_id = channel 
     !
     ! ket side if fully stored on each proc
     !
     bra_min = mapping_hpphpp(iam+1,local_ch_id,2)
     bra_max = mapping_hpphpp(iam+1,local_ch_id,3)
     !
     ! bra side is distributed 
     !
     ket_min = mapping_hpphpp(iam+1,local_ch_id,4)
     ket_max = mapping_hpphpp(iam+1,local_ch_id,5)
     
     ij_confs = numstates(3,tz, nx, ny, nz) 
     ab_confs = numstates(3,tz, nx, ny, nz)/2 
     mem = mem + dble(ab_confs)*3.*4./1.e9 
     mem_3nf = mem_3nf + dble((bra_max-bra_min+1)*ab_confs)*16./1.e9 
          
  end DO
  if ( iam == 0 ) write(6,*) 'Total memory for lookup hpphpp-confs', mem, 'Gbyte'
  if ( iam == 0 ) write(6,*) 'Total memory for v3nf_hpphpp', mem_3nf, 'Gbyte'
  
  call mpi_barrier(mpi_comm_world, ierror )
  mem = 0.d0  
  mem_3nf = 0.d0 
  do channel = my_hpphpp_channel_low(iam), my_hpphpp_channel_high(iam)
     nx = channels%hpphpp_quantum_numbers(channel*4)
     ny = channels%hpphpp_quantum_numbers(channel*4-1)
     nz = channels%hpphpp_quantum_numbers(channel*4-2)
     tz = channels%hpphpp_quantum_numbers(channel*4-3)
     
     local_ch_id = channel 
     !
     ! ket side if fully stored on each proc
     !
     bra_min = mapping_hpphpp(iam+1,local_ch_id,2)
     bra_max = mapping_hpphpp(iam+1,local_ch_id,3)
     !
     ! bra side is distributed 
     !
     ket_min = mapping_hpphpp(iam+1,local_ch_id,4)
     ket_max = mapping_hpphpp(iam+1,local_ch_id,5)
     
     ij_confs = numstates(3,tz, nx, ny, nz) 
     ab_confs = numstates(3,tz, nx, ny, nz) 
     mem = mem + dble(ij_confs+ab_confs)*3.*4./1.e9 
     mem_3nf = mem_3nf + dble((bra_max-bra_min+1)*(ket_max-ket_min+1))*16./1.e9 
     
  end DO
  if ( iam == 0 ) write(6,*) 'Total memory for lookup hpphpp-confs', mem, 'Gbyte'
  if ( iam == 0 ) write(6,*) 'Total memory for v3nf_hpphpp', mem_3nf, 'Gbyte'
  
  
  
  !
  ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !
  ! setup <hhp||hhp> channels and procs mapping
  !
  numchannels = 0 
  !     loop over isospin projection
  DO tz=tz_min3,tz_max3,2 
     !     loop over total momentum 
     DO nx=nx_min3, nx_max3 
        !     loop over total momentum 
        DO ny=ny_min3, ny_max3 
           !     loop over total momentum 
           DO nz=nz_min3, nz_max3 
              
              if ( numstates(2,tz, nx, ny, nz) /= 0 ) then
                 numchannels = numchannels + 1 
                 locate_v3channel(2,tz, nx, ny, nz) = numchannels
              end if
              
           end DO
        end DO
     end DO
  end DO
  
  
  mem = 0.d0 
  number_channels = 0 
  number_channels = maxval(locate_v3channel(2,:,:,:,:) ) 
  channels%number_hhphhp_confs  = number_channels
  ALLOCATE( channels%hhphhp_quantum_numbers(4*number_channels) )
  mem = mem + dble(number_channels)*4.*4./1.e9 
  
  if ( iam == 0 ) write(6,*) 'Total memory for t3-channels', mem, 'Gbyte' 
  
  !     loop over isospin projection
  DO tz=tz_min3,tz_max3,2 
     !     loop over total momentum 
     DO nx=nx_min3, nx_max3 
        !     loop over total momentum 
        DO ny=ny_min3, ny_max3 
           !     loop over total momentum 
           DO nz=nz_min3, nz_max3 
                 
              number_channels = locate_v3channel(2,tz,nx,ny,nz) 
              if ( number_channels == 0 ) cycle 
              
              k1 = number_channels*4 
              k2 = number_channels*4 - 1
              k3 = number_channels*4 - 2
              k4 = number_channels*4 - 3
              channels%hhphhp_quantum_numbers(k1) = nx
              channels%hhphhp_quantum_numbers(k2) = ny
              channels%hhphhp_quantum_numbers(k3) = nz
              channels%hhphhp_quantum_numbers(k4) = tz
                            
           end DO
        end DO
     end DO
  end DO

  ! 
  !  setup processor mapping and 
  !  calculate work for hhhppp blocks 
  !
  ltot_work=0
  tot_work = 0
  do channel   = 1, channels%number_hhphhp_confs 
     nx = channels%hhphhp_quantum_numbers(channel*4)
     ny = channels%hhphhp_quantum_numbers(channel*4-1)
     nz = channels%hhphhp_quantum_numbers(channel*4-2)
     tz = channels%hhphhp_quantum_numbers(channel*4-3)

     dim1 = numstates(2,tz, nx, ny, nz) 
     dim2 = numstates(2,tz, nx, ny, nz) 
     
     tot_work = tot_work + dim1*dim2
     ltot_work = ltot_work + int(dim1,8)*int(dim2, 8)
     
  end do
  
  number_channels = channels%number_hhphhp_confs 
  
  work = 0.d0
  memory = 0.d0
  
  work_pr_proc = int( ltot_work/int(num_procs,8),8 )
  work = 0.d0
  
  allocate( mapping_hhphhp(1:num_procs, 1:number_channels, 1:5) )
  mapping_hhphhp = 0 
  if ( iam == 0 ) write(6,*) ltot_work, work_pr_proc, number_channels 

  curr_work = 0
  curr_proc = 0
  local_ch_id = 1
  number_channels = 0
  do channel   = 1, channels%number_hhphhp_confs 
     nx = channels%hhphhp_quantum_numbers(channel*4)
     ny = channels%hhphhp_quantum_numbers(channel*4-1)
     nz = channels%hhphhp_quantum_numbers(channel*4-2)
     tz = channels%hhphhp_quantum_numbers(channel*4-3)
     
     dim1 = numstates(2,tz, nx, ny, nz) 
     dim2 = numstates(2,tz, nx, ny, nz) 
     
     !dim1 = size(  lookup_t3_configs(1,channel)%ival, 2) 
     !dim2 = size(  lookup_t3_configs(2,channel)%ival, 2) 
     
     bra_dim = dim1 
     ket_dim = dim2
     
     number_channels = number_channels + 1 
     has_ch_been_added = .FALSE. 
           
     i = 0
     DO WHILE ( i < ket_dim )
              
        i = i + 1
        temp = curr_work + bra_dim
              
        if ( temp <= work_pr_proc .or. curr_proc == num_procs - 1 ) then
           curr_work = temp
           if (  has_ch_been_added ) then
                    
              mapping_hhphhp(curr_proc+1,local_ch_id,3 ) = i 
              
           else
              
              
              has_ch_been_added = .true. 
              mapping_hhphhp(curr_proc+1,local_ch_id,1 ) = number_channels 
              mapping_hhphhp(curr_proc+1,local_ch_id,2 ) = i
              mapping_hhphhp(curr_proc+1,local_ch_id,3 ) = i
              mapping_hhphhp(curr_proc+1,local_ch_id,4 ) = 1
              mapping_hhphhp(curr_proc+1,local_ch_id,5 ) = bra_dim
              
           end if

           if ( i == ket_dim ) then
              
              local_ch_id = local_ch_id + 1
              
           end if
           work(curr_proc+1) = curr_work
           
        else 
                 
           has_ch_been_added = .false.
           work(curr_proc+1) = curr_work
           curr_work = 0
           curr_proc = curr_proc + 1
           i = i - 1
           
        end if
     end DO
  end DO
  
  number_channels = channels%number_hhphhp_confs 
  work = 0.d0
  do i = 1, num_procs
     do local_ch_id = 1, number_channels 
        if ( mapping_hhphhp(i,local_ch_id,1) == 0 ) cycle
        
        row_dim = mapping_hhphhp(i,local_ch_id,5)-mapping_hhphhp(i,local_ch_id,4)+1
        column_dim = mapping_hhphhp(i,local_ch_id,3)-mapping_hhphhp(i,local_ch_id,2)+1
        work(i) = work(i) + real( row_dim*column_dim )
        
     end do
  end do
  
  call mpi_barrier(mpi_comm_world,ierror)
  number_channels = channels%number_hhphhp_confs 
  
  allocate( my_hhphhp_channel_low(0:num_procs-1), my_hhphhp_channel_high(0:num_procs-1) )
  allocate( check_my_channel_hhphhp(1:number_channels) )
  check_my_channel_hhphhp = 0
  
  my_hhphhp_channel_low(iam)  = number_channels + 1
  my_hhphhp_channel_high(iam) = 0 
  do local_ch_id = 1, number_channels 
     if ( mapping_hhphhp(iam+1,local_ch_id,1) == 0 ) cycle
     check_my_channel_hhphhp(local_ch_id) = 1
      
     if ( local_ch_id > my_hhphhp_channel_high(iam) ) my_hhphhp_channel_high(iam) = local_ch_id
     if ( local_ch_id < my_hhphhp_channel_low(iam) ) my_hhphhp_channel_low(iam) = local_ch_id
     
  end do
  
  
  

  
  ALLOCATE( lookup_hhphhp_configs(1:1,my_hhphhp_channel_low(iam):my_hhphhp_channel_high(iam)) )
  mem = 0.d0  
  mem_3nf = 0.d0 
  do channel = my_hhphhp_channel_low(iam), my_hhphhp_channel_high(iam)
     nx = channels%hhphhp_quantum_numbers(channel*4)
     ny = channels%hhphhp_quantum_numbers(channel*4-1)
     nz = channels%hhphhp_quantum_numbers(channel*4-2)
     tz = channels%hhphhp_quantum_numbers(channel*4-3)
     
     local_ch_id = channel 
     !
     ! ket side if fully stored on each proc
     !
     bra_min = mapping_hhphhp(iam+1,local_ch_id,2)
     bra_max = mapping_hhphhp(iam+1,local_ch_id,3)
     !
     ! bra side is distributed 
     !
     ket_min = mapping_hhphhp(iam+1,local_ch_id,4)
     ket_max = mapping_hhphhp(iam+1,local_ch_id,5)
     
     ij_confs = numstates(2,tz, nx, ny, nz) 
     ab_confs = numstates(2,tz, nx, ny, nz) 
     mem = mem + dble(ij_confs+ab_confs)*3.*4./1.e9 
     mem_3nf = mem_3nf + dble((bra_max-bra_min+1)*(ket_max-ket_min+1))*16./1.e9 
     
  end DO
  if ( iam == 0 ) write(6,*) 'Total memory for lookup hhphhp-confs', mem, 'Gbyte'
  if ( iam == 0 ) write(6,*) 'Total memory for v3nf_hhphhp', mem_3nf, 'Gbyte'
  
  
  !
  ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !
  ! setup <hhp||ppp> channels and procs mapping
  !
  numchannels = 0 
  !     loop over isospin projection
  DO tz=tz_min3,tz_max3,2 
     !     loop over total momentum 
     DO nx=nx_min3, nx_max3 
        !     loop over total momentum 
        DO ny=ny_min3, ny_max3 
           !     loop over total momentum 
           DO nz=nz_min3, nz_max3 
              
              if ( numstates(2,tz, nx, ny, nz)*numstates(4,tz, nx, ny, nz) /= 0 ) then
                 numchannels = numchannels + 1 
                 locate_v3channel(3,tz, nx, ny, nz) = numchannels
              end if
              
           end DO
        end DO
     end DO
  end DO
  
  
  mem = 0.d0 
  number_channels = 0 
  number_channels = maxval(locate_v3channel(3,:,:,:,:) ) 
  channels%number_hhpppp_confs  = number_channels
  ALLOCATE( channels%hhpppp_quantum_numbers(4*number_channels) )
  mem = mem + dble(number_channels)*4.*4./1.e9 
  
  if ( iam == 0 ) write(6,*) 'Total memory for t3-channels', mem, 'Gbyte' 
  
  !     loop over isospin projection
  DO tz=tz_min3,tz_max3,2 
     !     loop over total momentum 
     DO nx=nx_min3, nx_max3 
        !     loop over total momentum 
        DO ny=ny_min3, ny_max3 
           !     loop over total momentum 
           DO nz=nz_min3, nz_max3 
                 
              number_channels = locate_v3channel(3,tz,nx,ny,nz) 
              if ( number_channels == 0 ) cycle 
              
              k1 = number_channels*4 
              k2 = number_channels*4 - 1
              k3 = number_channels*4 - 2
              k4 = number_channels*4 - 3
              channels%hhpppp_quantum_numbers(k1) = nx
              channels%hhpppp_quantum_numbers(k2) = ny
              channels%hhpppp_quantum_numbers(k3) = nz
              channels%hhpppp_quantum_numbers(k4) = tz
                            
           end DO
        end DO
     end DO
  end DO

  ! 
  !  setup processor mapping and 
  !  calculate work for hhhppp blocks 
  !
  ltot_work=0
  tot_work = 0
  do channel   = 1, channels%number_hhpppp_confs 
     nx = channels%hhpppp_quantum_numbers(channel*4)
     ny = channels%hhpppp_quantum_numbers(channel*4-1)
     nz = channels%hhpppp_quantum_numbers(channel*4-2)
     tz = channels%hhpppp_quantum_numbers(channel*4-3)

     dim1 = numstates(2,tz, nx, ny, nz)/2
     dim2 = numstates(4,tz, nx, ny, nz) 
     
     tot_work = tot_work + dim1*dim2
     ltot_work = ltot_work + int(dim1,8)*int(dim2, 8)
     
  end do
  
  number_channels = channels%number_hhpppp_confs 
  
  work = 0.d0
  memory = 0.d0
  
  work_pr_proc = int( ltot_work/int(num_procs,8),8 )
  work = 0.d0
  
  allocate( mapping_hhpppp(1:num_procs, 1:number_channels, 1:5) )
  mapping_hhpppp = 0 
  if ( iam == 0 ) write(6,*) ltot_work, work_pr_proc, number_channels 

  curr_work = 0
  curr_proc = 0
  local_ch_id = 1
  number_channels = 0
  do channel   = 1, channels%number_hhpppp_confs 
     nx = channels%hhpppp_quantum_numbers(channel*4)
     ny = channels%hhpppp_quantum_numbers(channel*4-1)
     nz = channels%hhpppp_quantum_numbers(channel*4-2)
     tz = channels%hhpppp_quantum_numbers(channel*4-3)
     
     dim1 = numstates(2,tz, nx, ny, nz) 
     dim2 = numstates(4,tz, nx, ny, nz) 
     
     !dim1 = size(  lookup_t3_configs(1,channel)%ival, 2) 
     !dim2 = size(  lookup_t3_configs(2,channel)%ival, 2) 
     
     bra_dim = dim1/2 
     ket_dim = dim2
     
     number_channels = number_channels + 1 
     has_ch_been_added = .FALSE. 
           
     i = 0
     DO WHILE ( i < ket_dim )
              
        i = i + 1
        temp = curr_work + bra_dim
              
        if ( temp <= work_pr_proc .or. curr_proc == num_procs - 1 ) then
           curr_work = temp
           if (  has_ch_been_added ) then
                    
              mapping_hhpppp(curr_proc+1,local_ch_id,3 ) = i 
              
           else
              
              
              has_ch_been_added = .true. 
              mapping_hhpppp(curr_proc+1,local_ch_id,1 ) = number_channels 
              mapping_hhpppp(curr_proc+1,local_ch_id,2 ) = i
              mapping_hhpppp(curr_proc+1,local_ch_id,3 ) = i
              mapping_hhpppp(curr_proc+1,local_ch_id,4 ) = 1
              mapping_hhpppp(curr_proc+1,local_ch_id,5 ) = bra_dim
              
           end if

           if ( i == ket_dim ) then
              
              local_ch_id = local_ch_id + 1
              
           end if
           work(curr_proc+1) = curr_work
           
        else 
                 
           has_ch_been_added = .false.
           work(curr_proc+1) = curr_work
           curr_work = 0
           curr_proc = curr_proc + 1
           i = i - 1
           
        end if
     end DO
  end DO
  
  number_channels = channels%number_hhpppp_confs 
  work = 0.d0
  do i = 1, num_procs
     do local_ch_id = 1, number_channels 
        if ( mapping_hhpppp(i,local_ch_id,1) == 0 ) cycle
        
        row_dim = mapping_hhpppp(i,local_ch_id,5)-mapping_hhpppp(i,local_ch_id,4)+1
        column_dim = mapping_hhpppp(i,local_ch_id,3)-mapping_hhpppp(i,local_ch_id,2)+1
        work(i) = work(i) + real( row_dim*column_dim )
        
     end do
  end do
  
  call mpi_barrier(mpi_comm_world,ierror)
  number_channels = channels%number_hhpppp_confs 
  
  allocate( my_hhpppp_channel_low(0:num_procs-1), my_hhpppp_channel_high(0:num_procs-1) )
  allocate( check_my_channel_hhpppp(1:number_channels) )
  check_my_channel_hhpppp = 0
  
  my_hhpppp_channel_low(iam)  = number_channels + 1
  my_hhpppp_channel_high(iam) = 0 
  do local_ch_id = 1, number_channels 
     if ( mapping_hhpppp(iam+1,local_ch_id,1) == 0 ) cycle
     check_my_channel_hhpppp(local_ch_id) = 1
      
     if ( local_ch_id > my_hhpppp_channel_high(iam) ) my_hhpppp_channel_high(iam) = local_ch_id
     if ( local_ch_id < my_hhpppp_channel_low(iam) ) my_hhpppp_channel_low(iam) = local_ch_id
     
  end do
  
  mem = 0.d0  
  mem_3nf = 0.d0 
  do channel = my_hhpppp_channel_low(iam), my_hhpppp_channel_high(iam)
     nx = channels%hhpppp_quantum_numbers(channel*4)
     ny = channels%hhpppp_quantum_numbers(channel*4-1)
     nz = channels%hhpppp_quantum_numbers(channel*4-2)
     tz = channels%hhpppp_quantum_numbers(channel*4-3)
     
     local_ch_id = channel 
     !
     ! ket side if fully stored on each proc
     !
     bra_min = mapping_hhpppp(iam+1,local_ch_id,2)
     bra_max = mapping_hhpppp(iam+1,local_ch_id,3)
     !
     ! bra side is distributed 
     !
     ket_min = mapping_hhpppp(iam+1,local_ch_id,4)
     ket_max = mapping_hhpppp(iam+1,local_ch_id,5)
     
     ij_confs = numstates(2,tz, nx, ny, nz)/2 
     ab_confs = numstates(4,tz, nx, ny, nz) 
     mem = mem + dble(ij_confs+ab_confs)*3.*4./1.e9 
     mem_3nf = mem_3nf + dble((bra_max-bra_min+1)*(ket_max-ket_min+1))*16./1.e9 
     
  end DO
  if ( iam == 0 ) write(6,*) 'Total memory for lookup hhpppp-confs', mem, 'Gbyte'
  if ( iam == 0 ) write(6,*) 'Total memory for v3nf_hhpppp', mem_3nf, 'Gbyte'
  
  
  !
  ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !
  ! setup <hhh||hpp> channels and procs mapping
  !
  numchannels = 0 
  !     loop over isospin projection
  DO tz=tz_min3,tz_max3,2 
     !     loop over total momentum 
     DO nx=nx_min3, nx_max3 
        !     loop over total momentum 
        DO ny=ny_min3, ny_max3 
           !     loop over total momentum 
           DO nz=nz_min3, nz_max3 
              
              if ( numstates(1,tz, nx, ny, nz)*numstates(3,tz, nx, ny, nz) /= 0 ) then
                 numchannels = numchannels + 1 
                 locate_v3channel(4,tz, nx, ny, nz) = numchannels
              end if
              
           end DO
        end DO
     end DO
  end DO
  
  
  mem = 0.d0 
  number_channels = 0 
  number_channels = maxval(locate_v3channel(4,:,:,:,:) ) 
  channels%number_hhhhpp_confs  = number_channels
  ALLOCATE( channels%hhhhpp_quantum_numbers(4*number_channels) )
  mem = mem + dble(number_channels)*4.*4./1.e9 
  
  if ( iam == 0 ) write(6,*) 'Total memory for t3-channels', mem, 'Gbyte' 
  
  !     loop over isospin projection
  DO tz=tz_min3,tz_max3,2 
     !     loop over total momentum 
     DO nx=nx_min3, nx_max3 
        !     loop over total momentum 
        DO ny=ny_min3, ny_max3 
           !     loop over total momentum 
           DO nz=nz_min3, nz_max3 
                 
              number_channels = locate_v3channel(4,tz,nx,ny,nz) 
              if ( number_channels == 0 ) cycle 
              
              k1 = number_channels*4 
              k2 = number_channels*4 - 1
              k3 = number_channels*4 - 2
              k4 = number_channels*4 - 3
              channels%hhhhpp_quantum_numbers(k1) = nx
              channels%hhhhpp_quantum_numbers(k2) = ny
              channels%hhhhpp_quantum_numbers(k3) = nz
              channels%hhhhpp_quantum_numbers(k4) = tz
                            
           end DO
        end DO
     end DO
  end DO

  ! 
  !  setup processor mapping and 
  !  calculate work for hhhppp blocks 
  !
  ltot_work=0
  tot_work = 0
  do channel   = 1, channels%number_hhhhpp_confs 
     nx = channels%hhhhpp_quantum_numbers(channel*4)
     ny = channels%hhhhpp_quantum_numbers(channel*4-1)
     nz = channels%hhhhpp_quantum_numbers(channel*4-2)
     tz = channels%hhhhpp_quantum_numbers(channel*4-3)

     dim1 = numstates(1,tz, nx, ny, nz) 
     dim2 = numstates(3,tz, nx, ny, nz) 
     
     tot_work = tot_work + dim1*dim2
     ltot_work = ltot_work + int(dim1,8)*int(dim2, 8)
     
  end do
  
  number_channels = channels%number_hhhhpp_confs 
  
  work = 0.d0
  memory = 0.d0
  
  work_pr_proc = int( ltot_work/int(num_procs,8),8 )
  work = 0.d0
  
  allocate( mapping_hhhhpp(1:num_procs, 1:number_channels, 1:5) )
  mapping_hhhhpp = 0 
  if ( iam == 0 ) write(6,*) ltot_work, work_pr_proc, number_channels 

  curr_work = 0
  curr_proc = 0
  local_ch_id = 1
  number_channels = 0
  do channel   = 1, channels%number_hhhhpp_confs 
     nx = channels%hhhhpp_quantum_numbers(channel*4)
     ny = channels%hhhhpp_quantum_numbers(channel*4-1)
     nz = channels%hhhhpp_quantum_numbers(channel*4-2)
     tz = channels%hhhhpp_quantum_numbers(channel*4-3)
     
     dim1 = numstates(1,tz, nx, ny, nz) 
     dim2 = numstates(3,tz, nx, ny, nz) 
     
     !dim1 = size(  lookup_t3_configs(1,channel)%ival, 2) 
     !dim2 = size(  lookup_t3_configs(2,channel)%ival, 2) 
     
     bra_dim = dim1 
     ket_dim = dim2
     
     number_channels = number_channels + 1 
     has_ch_been_added = .FALSE. 
           
     i = 0
     DO WHILE ( i < ket_dim )
              
        i = i + 1
        temp = curr_work + bra_dim
              
        if ( temp <= work_pr_proc .or. curr_proc == num_procs - 1 ) then
           curr_work = temp
           if (  has_ch_been_added ) then
                    
              mapping_hhhhpp(curr_proc+1,local_ch_id,3 ) = i 
              
           else
              
              
              has_ch_been_added = .true. 
              mapping_hhhhpp(curr_proc+1,local_ch_id,1 ) = number_channels 
              mapping_hhhhpp(curr_proc+1,local_ch_id,2 ) = i
              mapping_hhhhpp(curr_proc+1,local_ch_id,3 ) = i
              mapping_hhhhpp(curr_proc+1,local_ch_id,4 ) = 1
              mapping_hhhhpp(curr_proc+1,local_ch_id,5 ) = bra_dim
              
           end if

           if ( i == ket_dim ) then
              
              local_ch_id = local_ch_id + 1
              
           end if
           work(curr_proc+1) = curr_work
           
        else 
                 
           has_ch_been_added = .false.
           work(curr_proc+1) = curr_work
           curr_work = 0
           curr_proc = curr_proc + 1
           i = i - 1
           
        end if
     end DO
  end DO
  
  number_channels = channels%number_hhhhpp_confs 
  work = 0.d0
  do i = 1, num_procs
     do local_ch_id = 1, number_channels 
        if ( mapping_hhhhpp(i,local_ch_id,1) == 0 ) cycle
        
        row_dim = mapping_hhhhpp(i,local_ch_id,5)-mapping_hhhhpp(i,local_ch_id,4)+1
        column_dim = mapping_hhhhpp(i,local_ch_id,3)-mapping_hhhhpp(i,local_ch_id,2)+1
        work(i) = work(i) + real( row_dim*column_dim )
        
     end do
  end do
  
  call mpi_barrier(mpi_comm_world,ierror)
  number_channels = channels%number_hhhhpp_confs 
  
  allocate( my_hhhhpp_channel_low(0:num_procs-1), my_hhhhpp_channel_high(0:num_procs-1) )
  allocate( check_my_channel_hhhhpp(1:number_channels) )
  check_my_channel_hhhhpp = 0
  
  my_hhhhpp_channel_low(iam)  = number_channels + 1
  my_hhhhpp_channel_high(iam) = 0 
  do local_ch_id = 1, number_channels 
     if ( mapping_hhhhpp(iam+1,local_ch_id,1) == 0 ) cycle
     check_my_channel_hhhhpp(local_ch_id) = 1
      
     if ( local_ch_id > my_hhhhpp_channel_high(iam) ) my_hhhhpp_channel_high(iam) = local_ch_id
     if ( local_ch_id < my_hhhhpp_channel_low(iam) ) my_hhhhpp_channel_low(iam) = local_ch_id
     
  end do
  
  mem = 0.d0  
  mem_3nf = 0.d0 
  do channel = my_hhhhpp_channel_low(iam), my_hhhhpp_channel_high(iam)
     nx = channels%hhhhpp_quantum_numbers(channel*4)
     ny = channels%hhhhpp_quantum_numbers(channel*4-1)
     nz = channels%hhhhpp_quantum_numbers(channel*4-2)
     tz = channels%hhhhpp_quantum_numbers(channel*4-3)
     
     local_ch_id = channel 
     !
     ! ket side if fully stored on each proc
     !
     bra_min = mapping_hhhhpp(iam+1,local_ch_id,2)
     bra_max = mapping_hhhhpp(iam+1,local_ch_id,3)
     !
     ! bra side is distributed 
     !
     ket_min = mapping_hhhhpp(iam+1,local_ch_id,4)
     ket_max = mapping_hhhhpp(iam+1,local_ch_id,5)
     
     ij_confs = numstates(1,tz, nx, ny, nz) 
     ab_confs = numstates(3,tz, nx, ny, nz) 
     mem = mem + dble(ij_confs+ab_confs)*3.*4./1.e9 
     mem_3nf = mem_3nf + dble((bra_max-bra_min+1)*(Ket_max-ket_min+1))*16./1.e9 
  
     
  end DO
  if ( iam == 0 ) write(6,*) 'Total memory for lookup hhhhpp-confs', mem, 'Gbyte'
  if ( iam == 0 ) write(6,*) 'Total memory for v3nf_hhhhpp', mem_3nf, 'Gbyte'
  
  
  
  deallocate( work, memory )
  deallocate( numstates ) 
  
end SUBROUTINE compute_v3nf_memory

