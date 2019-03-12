!
!
!
SUBROUTINE setup_t3full_channel_structures
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
  INTEGER :: numchannels, ij_confs, ab_confs, channel, sz_min, sz_max, tz_min, tz_max, sumtz, sumsz, ichan
  integer(8), allocatable :: numstates(:,:,:,:,:)
  INTEGER, ALLOCATABLE :: numchannels1(:),amin(:), amax(:), bmin(:), bmax(:), cmin(:), cmax(:) 
  real*8 :: mem, delta 
  INTEGER :: k1,k2,k3,k4,k5, dim1
  INTEGER(8) :: total_work, work_per_proc, curr_work, curr_proc, curr_channel
  integer :: processor_work(num_procs), mtx_dim, mtx_dim_work, temp, local_ch_id
  integer(8) :: row_dim, column_dim, ket_dim, bra_dim
  INTEGER :: bra_min, bra_max, ket_min, ket_max, number_channels, dim2
  LOGICAL :: has_ch_been_added
  REAL*8, allocatable :: work(:), memory(:)
  integer(8):: ltot_work, tot_work, work_pr_proc
  
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
  do i = 1, tot_orbs !below_ef
     do j = 1, tot_orbs !below_ef
        do k = 1, tot_orbs !below_ef
           
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
  allocate( locate_t3channel(tz_min3:tz_max3, nx_min3:nx_max3, ny_min3:ny_max3, nz_min3:nz_max3) )
  locate_t3channel = 0  
  
  allocate( numstates(1:2,tz_min3:tz_max3, nx_min3:nx_max3, ny_min3:ny_max3, nz_min3:nz_max3) )
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
             
           if ( a <= below_ef .and. b <= below_ef .and. c <= below_ef ) & 
                numstates(1,tz, nx, ny, nz) = numstates(1,tz, nx, ny, nz)  + 1 
           if ( a > below_ef  .and. b >  below_ef .and. c > below_ef )  & 
                numstates(2,tz, nx, ny, nz) = numstates(2,tz, nx, ny, nz)  + 1 
           
        end do
     end do
  end do


  numchannels = 0 
  !     loop over isospin projection
  DO tz=tz_min3,tz_max3,2 
     !     loop over total momentum 
     DO nx=nx_min3, nx_max3 
        !     loop over total momentum 
        DO ny=ny_min3, ny_max3 
           !     loop over total momentum 
           DO nz=nz_min3, nz_max3 
                 
              if ( numstates(1,tz, nx, ny, nz)* numstates(2,tz, nx, ny, nz) /= 0 ) then
                 numchannels = numchannels + 1 
                 locate_t3channel(tz, nx, ny, nz) = numchannels
              end if
              
           end DO
        end DO
     end DO
  end DO
  
  
  mem = 0.d0 
  number_channels = 0 
  number_channels = maxval(locate_t3channel(:,:,:,:) ) 
  channels%number_t3_confs  = number_channels
  ALLOCATE( channels%t3_quantum_numbers(4*number_channels) )
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
                 
              number_channels = locate_t3channel(tz,nx,ny,nz) 
              if ( number_channels == 0 ) cycle 
              
              k1 = number_channels*4 
              k2 = number_channels*4 - 1
              k3 = number_channels*4 - 2
              k4 = number_channels*4 - 3
              channels%t3_quantum_numbers(k1) = nx
              channels%t3_quantum_numbers(k2) = ny
              channels%t3_quantum_numbers(k3) = nz
              channels%t3_quantum_numbers(k4) = tz
                            
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
  do channel   = 1, channels%number_t3_confs 
     nx = channels%t3_quantum_numbers(channel*4)
     ny = channels%t3_quantum_numbers(channel*4-1)
     nz = channels%t3_quantum_numbers(channel*4-2)
     tz = channels%t3_quantum_numbers(channel*4-3)

     dim1 = numstates(1,tz, nx, ny, nz) 
     dim2 = numstates(2,tz, nx, ny, nz) 
     
     !dim1 = size(  lookup_t3_configs(1,channel)%ival,2 )
     !dim2 = size(  lookup_t3_configs(2,channel)%ival,2 )
     
     
     tot_work = tot_work + dim1*dim2
     ltot_work = ltot_work + int(dim1,8)*int(dim2, 8)
     
  end do
  
  number_channels = channels%number_t3_confs 
  
  allocate( work(1:num_procs) )
  allocate( memory(1:num_procs) )
  work = 0.d0
  memory = 0.d0
  
  work_pr_proc = int( ltot_work/int(num_procs,8),8 )
  work = 0.d0
  
  allocate( mapping_t3(1:num_procs, 1:number_channels, 1:5) )
  mapping_t3 = 0 
  if ( iam == 0 ) write(6,*) ltot_work, work_pr_proc, number_channels 

  curr_work = 0
  curr_proc = 0
  local_ch_id = 1
  number_channels = 0
  do channel   = 1, channels%number_t3_confs 
     nx = channels%t3_quantum_numbers(channel*4)
     ny = channels%t3_quantum_numbers(channel*4-1)
     nz = channels%t3_quantum_numbers(channel*4-2)
     tz = channels%t3_quantum_numbers(channel*4-3)
     
     dim1 = numstates(1,tz, nx, ny, nz) 
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
                    
              mapping_t3(curr_proc+1,local_ch_id,3 ) = i 
              
           else
              
              
              has_ch_been_added = .true. 
              mapping_t3(curr_proc+1,local_ch_id,1 ) = number_channels 
              mapping_t3(curr_proc+1,local_ch_id,2 ) = i
              mapping_t3(curr_proc+1,local_ch_id,3 ) = i
              mapping_t3(curr_proc+1,local_ch_id,4 ) = 1
              mapping_t3(curr_proc+1,local_ch_id,5 ) = bra_dim
              
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
  
  number_channels = channels%number_t3_confs 
  work = 0.d0
  do i = 1, num_procs
     do local_ch_id = 1, number_channels 
        if ( mapping_t3(i,local_ch_id,1) == 0 ) cycle
        
        row_dim = mapping_t3(i,local_ch_id,5)-mapping_t3(i,local_ch_id,4)+1
        column_dim = mapping_t3(i,local_ch_id,3)-mapping_t3(i,local_ch_id,2)+1
        work(i) = work(i) + real( row_dim*column_dim )
        
     end do
  end do
  
!!$  write(6,*);write(6,*);write(6,*)
!!$  
  do i = 1, num_procs
     
!     if ( iam ==0 ) write(6,*)
!     if ( iam ==0 ) write(6,'(a12,2x,i2,2x,a10,2x,(E12.6))') 'proc number', & 
!          i-1, '% work', work(i) *100.d0/ real( work_pr_proc ) 
     do local_ch_id = 1, number_channels 
        if ( mapping_t3(i,local_ch_id,1) == 0 ) cycle
        
!        if ( iam ==0 ) write(6,*)
!        if ( iam ==0 ) write(6,'(a20,2x,i6)') 'Channel #:', mapping_t3(i,local_ch_id,1)
!        if ( iam ==0 ) write(6,'(a20,2x,i6)') 'Row index start:', mapping_t3(i,local_ch_id,4)
!        if ( iam ==0 ) write(6,'(a20,2x,i6)') 'Row index end:', mapping_t3(i,local_ch_id,5)
!        if ( iam ==0 ) write(6,'(a20,2x,i6)') 'Column index start:', mapping_t3(i,local_ch_id,2)
!        if ( iam ==0 ) write(6,'(a20,2x,i6)') 'Column index end:', mapping_t3(i,local_ch_id,3)
        
        
     end do
  end do
  
  call mpi_barrier(mpi_comm_world,ierror)
  number_channels = channels%number_t3_confs 

  allocate( my_t3channel_low(0:num_procs-1), my_t3channel_high(0:num_procs-1) )
  allocate( check_my_t3channel(1:number_channels) )
  check_my_t3channel = 0
  
  my_t3channel_low(iam)  = number_channels + 1
  my_t3channel_high(iam) = 0 
  do local_ch_id = 1, number_channels 
     if ( mapping_t3(iam+1,local_ch_id,1) == 0 ) cycle
     check_my_t3channel(local_ch_id) = 1
      
     if ( local_ch_id > my_t3channel_high(iam) ) my_t3channel_high(iam) = local_ch_id
     if ( local_ch_id < my_t3channel_low(iam) ) my_t3channel_low(iam) = local_ch_id
     
  end do
  
  call mpi_barrier(mpi_comm_world, ierror)
  !write(6,*) iam, my_t3channel_low(iam), my_t3channel_high(iam)
  allocate( t3_ccm(my_t3channel_low(iam):my_t3channel_high(iam)) )
  allocate( v3nf_ppphhh(my_t3channel_low(iam):my_t3channel_high(iam)) )
  call mpi_barrier(mpi_comm_world, ierror)
  
  deallocate( work, memory )

  
  ALLOCATE( lookup_t3_configs(1:2,my_t3channel_low(iam):my_t3channel_high(iam)) )
  mem = 0.d0  
  do channel = my_t3channel_low(iam), my_t3channel_high(iam)
     nx = channels%t3_quantum_numbers(channel*4)
     ny = channels%t3_quantum_numbers(channel*4-1)
     nz = channels%t3_quantum_numbers(channel*4-2)
     tz = channels%t3_quantum_numbers(channel*4-3)
     
     local_ch_id = channel 
     !
     ! ket side if fully stored on each proc
     !
     bra_min = mapping_t3(iam+1,local_ch_id,2)
     bra_max = mapping_t3(iam+1,local_ch_id,3)
     !
     ! bra side is distributed 
     !
     ket_min = mapping_t3(iam+1,local_ch_id,4)
     ket_max = mapping_t3(iam+1,local_ch_id,5)
     
     ij_confs = numstates(1,tz, nx, ny, nz) 
     ab_confs = numstates(2,tz, nx, ny, nz) 
     mem = mem + dble(ij_confs+ab_confs)*3.*4./1.e9 
     ALLOCATE( lookup_t3_configs(1,channel)%ival(3,ij_confs))
     ALLOCATE( lookup_t3_configs(2,channel)%ival(3,ab_confs))
     
  end DO
  if ( iam == 0 ) write(6,*) 'Total memory for lookup t3-confs', mem, 'Gbyte'
  deallocate( numstates ) 
  

  allocate( numchannels1( maxval( locate_t3channel(:,:,:,:) ) ) )

  
  !
  ! ppp configurations
  !
  allocate( ppp_hhhppp%ival3(below_ef+1:tot_orbs, below_ef+1:tot_orbs, 0:1) )
  ppp_hhhppp%ival3 = 0 
  
  numchannels1 = 0
  do a = below_ef+1,tot_orbs
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
           
           
           channel = locate_t3channel(tz, nx, ny, nz) 
           if ( channel == 0 ) cycle 
           if ( check_my_t3channel(channel) == 0 ) cycle 
           
           numchannels1(channel) = numchannels1(channel) + 1 
           ichan = numchannels1(channel)
           !ppp_hhhppp%ival3(a,b,(szc+1)/2) = ichan 
!!$           if ( a > amax(channel) ) amax(channel) = a 
!!$           if ( a < amin(channel) ) amin(channel) = a 
!!$           if ( b > bmax(channel) ) bmax(channel) = b 
!!$           if ( b < bmin(channel) ) bmin(channel) = b 
!!$           if ( c > cmax(channel) ) cmax(channel) = c 
!!$           if ( c < cmin(channel) ) cmin(channel) = c 
           
           
           lookup_t3_configs(2,channel)%ival(1,ichan) = a 
           lookup_t3_configs(2,channel)%ival(2,ichan) = b 
           lookup_t3_configs(2,channel)%ival(3,ichan) = c 
           
        end do
     end do
  end do
  
!!$  do channel = 1,  channels%number_t3_confs
!!$     write(6,'(9(I3,1x))') channel, below_ef, tot_orbs, amin(channel), amax(channel), bmin(channel), bmax(channel), cmin(channel), cmax(channel) 
!!$  end do

  allocate( hhh_hhhppp%ival3(below_ef, below_ef, below_ef) )
  hhh_hhhppp%ival3 = 0 
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
           
           channel = locate_t3channel(tz, nx, ny, nz) 
           if ( channel == 0 ) cycle 
           if ( check_my_t3channel(channel) == 0 ) cycle 
           
           numchannels1(channel) = numchannels1(channel) + 1 
           ichan = numchannels1(channel)
           
           hhh_hhhppp%ival3(a,b,c) = ichan
           lookup_t3_configs(1,channel)%ival(1,ichan) = a 
           lookup_t3_configs(1,channel)%ival(2,ichan) = b 
           lookup_t3_configs(1,channel)%ival(3,ichan) = c 
           
        end do
     end do
  end do
  
  deallocate( numchannels1 )

  
end SUBROUTINE setup_t3full_channel_structures



subroutine mbpt2_3nf
  USE PARALLEL
  USE single_particle_orbits
  USE constants
  use one_body_operators
  use wave_functions
  use chiral_potentials
  use configurations
  use t3_constants
  use t2_storage
  
  IMPLICIT NONE
  INTEGER :: d, bra_confs, ket_confs, bra,ket, ii
  complex*16  :: spen, et_4tmp,et_4tmp2, et_4, et_4t, tmp_sum1, tmp_sum2, et_5,  & 
       est_5, sum1,sum2,sum3,sum4,sum5,sum6, sum7,sum8,sum9,tc_term, td_term,lc_term
  integer, allocatable :: nconf_low(:), nconf_high(:)
  integer :: nconfs, nconfs_tot, number_mtxel_iam, diff, n
  real*8  ::  delta_p, startwtime , endwtime
  complex*16 :: diag1 
  INTEGER :: i,j,k, sumx, sumy, sumz, a,b,c,nxa,nxb,nxc,nya,nyb,nyc,nza,nzb,nzc,tza,tzb,tzc,sza,szb,szc,tz,sz,nx,ny,nz
  INTEGER :: numchannels, ij_confs, ab_confs, channel, sz_min, sz_max, tz_min, tz_max, sumtz, sumsz, ichan
  integer(8), allocatable :: numstates(:,:,:,:,:)
  INTEGER, ALLOCATABLE :: numchannels1(:),amin(:), amax(:), bmin(:), bmax(:), cmin(:), cmax(:) 
  real*8 :: mem,dum
  INTEGER :: k1,k2,k3,k4,k5, dim1
  INTEGER(8) :: total_work, work_per_proc, curr_work, curr_proc, curr_channel
  integer :: processor_work(num_procs), mtx_dim, mtx_dim_work, temp, local_ch_id
  integer(8) :: row_dim, column_dim, ket_dim, bra_dim
  INTEGER :: bra_min, bra_max, ket_min, ket_max, number_channels, dim2
  LOGICAL :: has_ch_been_added
  REAL*8, allocatable :: work(:), memory(:)
  integer(8):: ltot_work, tot_work, work_pr_proc
  complex*16, allocatable :: tmp_sum(:), et_4tmp22(:)
  
  
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
  do i = 1, tot_orbs !below_ef
     do j = i+1, tot_orbs !below_ef
        do k = j+1, tot_orbs !below_ef
           
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
  allocate( locate_t3channel(tz_min3:tz_max3, nx_min3:nx_max3, ny_min3:ny_max3, nz_min3:nz_max3) )
  locate_t3channel = 0  
  
  allocate( numstates(1:2,tz_min3:tz_max3, nx_min3:nx_max3, ny_min3:ny_max3, nz_min3:nz_max3) )
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
  
        do c = b+1, tot_orbs
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
             
           if ( a <= below_ef .and. b <= below_ef .and. c <= below_ef ) & 
                numstates(1,tz, nx, ny, nz) = numstates(1,tz, nx, ny, nz)  + 1 
           if ( a > below_ef  .and. b >  below_ef .and. c > below_ef ) & 
                numstates(2,tz, nx, ny, nz) = numstates(2,tz, nx, ny, nz)  + 1 
           
        end do
     end do
  end do


  numchannels = 0 
  !     loop over isospin projection
  DO tz=tz_min3,tz_max3,2 
     !     loop over total momentum 
     DO nx=nx_min3, nx_max3 
        !     loop over total momentum 
        DO ny=ny_min3, ny_max3 
           !     loop over total momentum 
           DO nz=nz_min3, nz_max3 
                 
              if ( numstates(1,tz, nx, ny, nz)* numstates(2,tz, nx, ny, nz) /= 0 ) then
                 numchannels = numchannels + 1 
                 locate_t3channel(tz, nx, ny, nz) = numchannels
              end if
              
           end DO
        end DO
     end DO
  end DO
  
  
  mem = 0.d0 
  number_channels = 0 
  number_channels = maxval(locate_t3channel(:,:,:,:) ) 
  channels%number_t3_confs  = number_channels
  ALLOCATE( channels%t3_quantum_numbers(4*number_channels) )
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
                 
              number_channels = locate_t3channel(tz,nx,ny,nz) 
              if ( number_channels == 0 ) cycle 
              
              k1 = number_channels*4 
              k2 = number_channels*4 - 1
              k3 = number_channels*4 - 2
              k4 = number_channels*4 - 3
              channels%t3_quantum_numbers(k1) = nx
              channels%t3_quantum_numbers(k2) = ny
              channels%t3_quantum_numbers(k3) = nz
              channels%t3_quantum_numbers(k4) = tz
                            
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
  do channel   = 1, channels%number_t3_confs 
     nx = channels%t3_quantum_numbers(channel*4)
     ny = channels%t3_quantum_numbers(channel*4-1)
     nz = channels%t3_quantum_numbers(channel*4-2)
     tz = channels%t3_quantum_numbers(channel*4-3)

     dim1 = numstates(1,tz, nx, ny, nz) 
     dim2 = numstates(2,tz, nx, ny, nz) 
     
     !dim1 = size(  lookup_t3_configs(1,channel)%ival,2 )
     !dim2 = size(  lookup_t3_configs(2,channel)%ival,2 )
     
     
     tot_work = tot_work + dim1*dim2
     ltot_work = ltot_work + int(dim1,8)*int(dim2, 8)
     
  end do
  
  number_channels = channels%number_t3_confs 
  
  allocate( work(1:num_procs) )
  allocate( memory(1:num_procs) )
  work = 0.d0
  memory = 0.d0
  
  work_pr_proc = int( ltot_work/int(num_procs,8),8 )
  work = 0.d0
  
  allocate( mapping_t3(1:num_procs, 1:number_channels, 1:5) )
  mapping_t3 = 0 
  if ( iam == 0 ) write(6,*) ltot_work, work_pr_proc, number_channels 

  curr_work = 0
  curr_proc = 0
  local_ch_id = 1
  number_channels = 0
  do channel   = 1, channels%number_t3_confs 
     nx = channels%t3_quantum_numbers(channel*4)
     ny = channels%t3_quantum_numbers(channel*4-1)
     nz = channels%t3_quantum_numbers(channel*4-2)
     tz = channels%t3_quantum_numbers(channel*4-3)
     
     dim1 = numstates(1,tz, nx, ny, nz) 
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
                    
              mapping_t3(curr_proc+1,local_ch_id,3 ) = i 
              
           else
              
              
              has_ch_been_added = .true. 
              mapping_t3(curr_proc+1,local_ch_id,1 ) = number_channels 
              mapping_t3(curr_proc+1,local_ch_id,2 ) = i
              mapping_t3(curr_proc+1,local_ch_id,3 ) = i
              mapping_t3(curr_proc+1,local_ch_id,4 ) = 1
              mapping_t3(curr_proc+1,local_ch_id,5 ) = bra_dim
              
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
  
  number_channels = channels%number_t3_confs 
  work = 0.d0
  do i = 1, num_procs
     do local_ch_id = 1, number_channels 
        if ( mapping_t3(i,local_ch_id,1) == 0 ) cycle
        
        row_dim = mapping_t3(i,local_ch_id,5)-mapping_t3(i,local_ch_id,4)+1
        column_dim = mapping_t3(i,local_ch_id,3)-mapping_t3(i,local_ch_id,2)+1
        work(i) = work(i) + real( row_dim*column_dim )
        
     end do
  end do
  
!!$  write(6,*);write(6,*);write(6,*)
!!$  
  do i = 1, num_procs
     
!     if ( iam ==0 ) write(6,*)
!     if ( iam ==0 ) write(6,'(a12,2x,i2,2x,a10,2x,(E12.6))') 'proc number', & 
!          i-1, '% work', work(i) *100.d0/ real( work_pr_proc ) 
     do local_ch_id = 1, number_channels 
        if ( mapping_t3(i,local_ch_id,1) == 0 ) cycle
        
!        if ( iam ==0 ) write(6,*)
!        if ( iam ==0 ) write(6,'(a20,2x,i6)') 'Channel #:', mapping_t3(i,local_ch_id,1)
!        if ( iam ==0 ) write(6,'(a20,2x,i6)') 'Row index start:', mapping_t3(i,local_ch_id,4)
!        if ( iam ==0 ) write(6,'(a20,2x,i6)') 'Row index end:', mapping_t3(i,local_ch_id,5)
!        if ( iam ==0 ) write(6,'(a20,2x,i6)') 'Column index start:', mapping_t3(i,local_ch_id,2)
!        if ( iam ==0 ) write(6,'(a20,2x,i6)') 'Column index end:', mapping_t3(i,local_ch_id,3)
        
        
     end do
  end do
  
  call mpi_barrier(mpi_comm_world,ierror)
  number_channels = channels%number_t3_confs 

  allocate( my_t3channel_low(0:num_procs-1), my_t3channel_high(0:num_procs-1) )
  allocate( check_my_t3channel(1:number_channels) )
  check_my_t3channel = 0
  
  my_t3channel_low(iam)  = number_channels + 1
  my_t3channel_high(iam) = 0 
  do local_ch_id = 1, number_channels 
     if ( mapping_t3(iam+1,local_ch_id,1) == 0 ) cycle
     check_my_t3channel(local_ch_id) = 1
      
     if ( local_ch_id > my_t3channel_high(iam) ) my_t3channel_high(iam) = local_ch_id
     if ( local_ch_id < my_t3channel_low(iam) ) my_t3channel_low(iam) = local_ch_id
     
  end do
  
  call mpi_barrier(mpi_comm_world, ierror)
  !write(6,*) iam, my_t3channel_low(iam), my_t3channel_high(iam)
    
  deallocate( work, memory )

  
  ALLOCATE( lookup_t3_configs(1:2,my_t3channel_low(iam):my_t3channel_high(iam)) )
  mem = 0.d0  
  do channel = my_t3channel_low(iam), my_t3channel_high(iam)
     nx = channels%t3_quantum_numbers(channel*4)
     ny = channels%t3_quantum_numbers(channel*4-1)
     nz = channels%t3_quantum_numbers(channel*4-2)
     tz = channels%t3_quantum_numbers(channel*4-3)
     
     local_ch_id = channel 
     !
     ! ket side if fully stored on each proc
     !
     bra_min = mapping_t3(iam+1,local_ch_id,2)
     bra_max = mapping_t3(iam+1,local_ch_id,3)
     !
     ! bra side is distributed 
     !
     ket_min = mapping_t3(iam+1,local_ch_id,4)
     ket_max = mapping_t3(iam+1,local_ch_id,5)
     
     ij_confs = numstates(1,tz, nx, ny, nz) 
     ab_confs = numstates(2,tz, nx, ny, nz) 
     mem = mem + dble(ij_confs+ab_confs)*3.*4./1.e9 
     ALLOCATE( lookup_t3_configs(1,channel)%ival(3,ij_confs))
     ALLOCATE( lookup_t3_configs(2,channel)%ival(3,ab_confs))
     
  end DO
  if ( iam == 0 ) write(6,*) 'Total memory for lookup t3-confs', mem, 'Gbyte'
  deallocate( numstates ) 
  
  


  allocate( numchannels1( maxval( locate_t3channel(:,:,:,:) ) ) )

  
  numchannels1 = 0
  do a = below_ef+1,tot_orbs
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
        do c = b+1, tot_orbs
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
           
           
           channel = locate_t3channel(tz, nx, ny, nz) 
           if ( channel == 0 ) cycle 
           if ( check_my_t3channel(channel) == 0 ) cycle 
           
           numchannels1(channel) = numchannels1(channel) + 1 
           ichan = numchannels1(channel)
           
           
           lookup_t3_configs(2,channel)%ival(1,ichan) = a 
           lookup_t3_configs(2,channel)%ival(2,ichan) = b 
           lookup_t3_configs(2,channel)%ival(3,ichan) = c 
           
        end do
     end do
  end do
  
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
     
     do b = a+1, below_ef
        nxb = all_orbit%nx(b)
        nyb = all_orbit%ny(b)
        nzb = all_orbit%nz(b)
        tzb = all_orbit%itzp(b)
        szb = all_orbit%szp(b)
        
        do c = b+1, below_ef
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
           
           channel = locate_t3channel(tz, nx, ny, nz) 
           if ( channel == 0 ) cycle 
           if ( check_my_t3channel(channel) == 0 ) cycle 
           
           numchannels1(channel) = numchannels1(channel) + 1 
           ichan = numchannels1(channel)
           
           lookup_t3_configs(1,channel)%ival(1,ichan) = a 
           lookup_t3_configs(1,channel)%ival(2,ichan) = b 
           lookup_t3_configs(1,channel)%ival(3,ichan) = c 
           
        end do
     end do
  end do
  
  deallocate( numchannels1 )
  
  allocate( et_4tmp22(0:30), tmp_sum(0:30) ) 
  tmp_sum = 0.d0 
  delta_p = (6.-kf)/30. 
  
  call mpi_barrier(mpi_comm_world,ierror)
  startwtime = MPI_WTIME()
  ket_confs = 0 
  bra_confs = 0 
  et_4 = 0.d0 
  et_5 = 0.d0 
  et_4tmp = 0.d0
  et_4tmp2  = 0.d0
  et_4tmp22  = 0.d0
  do channel   =  my_t3channel_low(iam), my_t3channel_high(iam)
     nx = channels%t3_quantum_numbers(channel*4)
     ny = channels%t3_quantum_numbers(channel*4-1)
     nz = channels%t3_quantum_numbers(channel*4-2)
     tz = channels%t3_quantum_numbers(channel*4-3)
     
     
     local_ch_id = channel
     !
     ! ket side if fully stored on each proc
     !
     bra_min = mapping_t3(iam+1,local_ch_id,2)
     bra_max = mapping_t3(iam+1,local_ch_id,3)
     !
     ! bra side is distributed 
     !
     ket_min = mapping_t3(iam+1,local_ch_id,4)
     ket_max = mapping_t3(iam+1,local_ch_id,5)
     
     if ( iam == 0 ) write(6,*) channel, size(  lookup_t3_configs(2,channel)%ival, 2), size(  lookup_t3_configs(1,channel)%ival, 2) 
     
     !tmp_sum2 =dcmplx(0.d0,0.d0) 
     tmp_sum = 0.d0 
     !$omp parallel default(shared) private(ket,a,b,c,bra,i,j,k,spen, & 
     !$omp td_term, ii)
     !$omp do schedule(dynamic), reduction(+:tmp_sum)
     do ket = bra_min, bra_max 
        a = lookup_t3_configs(2,channel)%ival(1,ket)
        b = lookup_t3_configs(2,channel)%ival(2,ket) 
        c = lookup_t3_configs(2,channel)%ival(3,ket) 
        
        
        do bra = 1, size(  lookup_t3_configs(1,channel)%ival, 2) 
           i = lookup_t3_configs(1,channel)%ival(1,bra) 
           j = lookup_t3_configs(1,channel)%ival(2,bra) 
           k = lookup_t3_configs(1,channel)%ival(3,bra) 
           
           
           spen = ( fock_mtx(i,i) + fock_mtx(j,j)  + fock_mtx(k,k) & 
                - fock_mtx(a,a)-fock_mtx(b,b) - fock_mtx(c,c) )
           
           td_term = 0.d0 
           td_term = chiral_3nf_asym(a,b,c,i,j,k)
           !tmp_sum2 = tmp_sum2 + conjg(td_term)*(td_term) / spen 
           do ii=0, 30      
              if ( sqrt( all_orbit%kx(a)**2+all_orbit%ky(a)**2 +all_orbit%kz(a)**2 ) <= kf+delta_p*ii .and. & 
                   sqrt( all_orbit%kx(b)**2+all_orbit%ky(b)**2 +all_orbit%kz(b)**2 ) <= kf+delta_p*ii .and. & 
                   sqrt( all_orbit%kx(c)**2+all_orbit%ky(c)**2 +all_orbit%kz(c)**2 ) <= kf+delta_p*ii ) then 
                 tmp_sum(ii) = tmp_sum(ii) + conjg(td_term)*(td_term) / spen 
              end if
           end do
           
        end DO
     end DO
     !$omp end do
     !$omp end parallel
     
     et_4tmp22  = et_4tmp22  + tmp_sum
  end do
  
  
!!$  dum = -1.d0
!!$
!!$  startwtime = MPI_WTIME()
!!$  ket_confs = 0 
!!$  bra_confs = 0 
!!$  et_4 = 0.d0 
!!$  et_5 = 0.d0 
!!$  et_4tmp = 0.d0
!!$  et_4tmp2  = 0.d0
!!$  do channel   =  my_t3channel_low(iam), my_t3channel_high(iam)
!!$     nx = channels%t3_quantum_numbers(channel*4)
!!$     ny = channels%t3_quantum_numbers(channel*4-1)
!!$     nz = channels%t3_quantum_numbers(channel*4-2)
!!$     tz = channels%t3_quantum_numbers(channel*4-3)
!!$     
!!$     
!!$     local_ch_id = channel
!!$     !
!!$     ! ket side if fully stored on each proc
!!$     !
!!$     bra_min = mapping_t3(iam+1,local_ch_id,2)
!!$     bra_max = mapping_t3(iam+1,local_ch_id,3)
!!$     !
!!$     ! bra side is distributed 
!!$     !
!!$     ket_min = mapping_t3(iam+1,local_ch_id,4)
!!$     ket_max = mapping_t3(iam+1,local_ch_id,5)
!!$     
!!$     if ( iam == 0 ) write(6,*) channel, size(  lookup_t3_configs(2,channel)%ival, 2), size(  lookup_t3_configs(1,channel)%ival, 2) 
!!$     
!!$     tmp_sum2 =dcmplx(0.d0,0.d0) 
!!$     
!!$     do bra = 1, size(  lookup_t3_configs(1,channel)%ival, 2) 
!!$        i = lookup_t3_configs(1,channel)%ival(1,bra) 
!!$        j = lookup_t3_configs(1,channel)%ival(2,bra) 
!!$        k = lookup_t3_configs(1,channel)%ival(3,bra) 
!!$        
!!$        sum1 = 0.d0 
!!$        !$omp parallel default(shared) private(ket,a,b,c,spen, & 
!!$        !$omp td_term)
!!$        !$omp do schedule(dynamic), reduction(+:sum1)
!!$        do ket = bra_min, bra_max 
!!$           a = lookup_t3_configs(2,channel)%ival(1,ket)
!!$           b = lookup_t3_configs(2,channel)%ival(2,ket) 
!!$           c = lookup_t3_configs(2,channel)%ival(3,ket) 
!!$           
!!$           
!!$           spen = ( fock_mtx(i,i) + fock_mtx(j,j)  + fock_mtx(k,k) & 
!!$                - fock_mtx(a,a)-fock_mtx(b,b) - fock_mtx(c,c) )
!!$           
!!$           td_term = 0.d0 
!!$           td_term = chiral_3nf_asym(a,b,c,i,j,k)
!!$           sum1 = sum1 + conjg(td_term)*(td_term) / spen 
!!$           
!!$        end DO
!!$        !$omp end do
!!$        !$omp end parallel
!!$        
!!$        if ( abs(sum1) > dum ) then 
!!$           dum = abs(sum1) 
!!$           write(62,'(3g,1x,2I,1x,3g,1x,2I,3g,1x,2I,g)') all_orbit%kx(i), all_orbit%ky(i), all_orbit%kz(i), all_orbit%szp(i), all_orbit%itzp(i), & 
!!$                all_orbit%kx(j), all_orbit%ky(j), all_orbit%kz(j), all_orbit%szp(j), all_orbit%itzp(j), & 
!!$                all_orbit%kx(k), all_orbit%ky(k), all_orbit%kz(k), all_orbit%szp(k), all_orbit%itzp(k), & 
!!$                dble(sum1) 
!!$
!!$        end if
!!$        
!!$     end DO
!!$  end do
!!$  
  
  call mpi_allreduce(mpi_in_place,et_4tmp22,size(et_4tmp22),mpi_complex16,mpi_sum, &
       mpi_comm_world,ierror)
  
  call mpi_barrier(mpi_comm_world,ierror)
  endwtime = MPI_WTIME()
  IF ( IAM == 0 ) WRITE(6,*) 'EXECUTION TIME FOR CCD(T)', endwtime-startwtime
  !if ( iam == 0 ) write(6,'(a,1x,2(F10.5,2x),a,1x,4(F10.5,3x))') 'Kf, Dens', kf, dens, 'HF, MBPT2, CCD, CCD(T)', e0, mbpt2, e0 + eccsd, e0 + eccsd + et_4
  if ( iam == 0 ) write(6,'(a,1x,2(F16.5,2x),a,1x,2(F16.5,3x))') 'Kf, Dens', kf, dens, 'MBPT2 3NF', & 
       real(et_5/below_ef), aimag(et_5/below_ef)
  do ii = 0, 30
     if (iam == 0 ) write(6,*) kf+delta_p*ii, dble(et_4tmp22(ii)/below_ef) 
  end do
  
end subroutine mbpt2_3nf


