
SUBROUTINE setup_proc_mappings_pp
  USE single_particle_orbits
  USE constants
  use one_body_operators
  use t2_storage
  use configurations
  USE KSPACE 
  use CHIRAL_POTENTIALS
  use parallel

  IMPLICIT NONE
  INTEGER :: channel, channel2, dim1, nx,ny,nz,sz,tz
  INTEGER :: total_work, work_per_proc, curr_work, curr_proc, curr_channel
  integer :: processor_work(num_procs), mtx_dim, mtx_dim_work, temp, tot_work, work_pr_proc, local_ch_id
  integer :: row_dim, column_dim, ket_dim, bra_dim
  INTEGER :: bra_min, bra_max, ket_min, ket_max, number_channels, i, dim2
  LOGICAL :: has_ch_been_added
  REAL*8, allocatable :: work(:), memory(:)
  integer(8):: ltot_work
  
  allocate( work(1:num_procs) )
  allocate( memory(1:num_procs) )

  ltot_work=0
  tot_work = 0
if(cc_approx .ne. 'mbpt2')then
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
     
     
     dim1 = size(  lookup_pppp_configs(1,channel)%ival, 2) 
          
     tot_work = tot_work + dim1**2
     ltot_work = ltot_work + int(dim1,8)**2 
     
  end do
  
  work_pr_proc = int( tot_work/num_procs )
  work_pr_proc = int( ltot_work/int(num_procs,8) )
  
  work = 0.d0
  memory = 0.d0
  
  number_channels = channels%number_pppp_confs 
  allocate( mapping(1:num_procs, 1:number_channels, 1:5) )
  mapping = 0 
  
  curr_work = 0
  curr_proc = 0
  local_ch_id = 1
  number_channels = 0
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
     
     mtx_dim = size(  lookup_pppp_configs(1,channel)%ival, 2) 
     mtx_dim_work = size(  lookup_pppp_configs(1,channel)%ival, 2) 
     number_channels = number_channels + 1 
     has_ch_been_added = .FALSE. 
     
     i = 0
     DO WHILE ( i < mtx_dim )
              
        i = i + 1
        temp = curr_work + mtx_dim_work 
        
        if ( temp <= work_pr_proc .or. curr_proc == num_procs - 1 ) then
           curr_work = temp
           if (  has_ch_been_added ) then
                    
              mapping(curr_proc+1,local_ch_id,5 ) = i 
              
           else
              
                    
              has_ch_been_added = .true. 
              mapping(curr_proc+1,local_ch_id,1 ) = number_channels 
              mapping(curr_proc+1,local_ch_id,2 ) = 1
              mapping(curr_proc+1,local_ch_id,3 ) = mtx_dim
              mapping(curr_proc+1,local_ch_id,4 ) = i
              mapping(curr_proc+1,local_ch_id,5 ) = i 
                    
           end if

           if ( i == mtx_dim ) then
                    
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
  end do
  work = 0.d0
  memory = 0.d0
  do i = 1, num_procs
     do local_ch_id = 1, number_channels 
        if ( mapping(i,local_ch_id,1) == 0 ) cycle
        
        row_dim = mapping(i,local_ch_id,5)-mapping(i,local_ch_id,4)+1
        column_dim = mapping(i,local_ch_id,3)-mapping(i,local_ch_id,2)+1
        work(i) = work(i) + real( row_dim*column_dim )
        !work(i) = work(i) + real( row_dim*column_dim * ket_dim_channel(local_ch_id) )
        memory(i) = memory(i)  + real( row_dim*column_dim )*8.d0/1.e9
        
     end do
  end do

!!$  do i = 1, num_procs
!!$     
!!$     if ( iam ==0 ) write(6,*)
!!$     if ( iam ==0 ) write(6,'(a12,2x,i2,2x,a10,2x,(E12.6))') 'proc number', & 
!!$          i-1, '% work', work(i) *100.d0/ real( work_pr_proc ) 
!!$     if ( iam ==0 ) write(6,'(a12,2x,i2,2x,a10,2x,(E12.6))') 'proc number', & 
!!$          i-1, 'memory in Gb', memory(i)
!!$     
!!$     do local_ch_id = 1, number_channels 
!!$        if ( mapping(i,local_ch_id,1) == 0 ) cycle
!!$        
!!$        if ( iam ==0 ) write(6,*)
!!$        if ( iam ==0 ) write(6,'(a20,2x,i4)') 'Channel #:', mapping(i,local_ch_id,1)
!!$        if ( iam ==0 ) write(6,'(a20,2x,i4)') 'Column index start:', mapping(i,local_ch_id,2)
!!$        if ( iam ==0 ) write(6,'(a20,2x,i4)') 'Column index end:', mapping(i,local_ch_id,3)
!!$        if ( iam ==0 ) write(6,'(a20,2x,i4)') 'Row index start:', mapping(i,local_ch_id,4)
!!$        if ( iam ==0 ) write(6,'(a20,2x,i4)') 'Row index end:', mapping(i,local_ch_id,5)
!!$        
!!$        
!!$     end do
!!$  end do

  number_channels = channels%number_pppp_confs 
  
  
  allocate( check_my_channel(1:number_channels) )
  allocate( fully_stored_channel(1:number_channels) )
  allocate( my_channel_low(0:num_procs-1), my_channel_high(0:num_procs-1))
  fully_stored_channel = 0
  check_my_channel = 0
  
  my_channel_low(iam)  = number_channels + 1
  my_channel_high(iam) = 0 
  do local_ch_id = 1, number_channels 
     if ( mapping(iam+1,local_ch_id,1) == 0 ) cycle
     check_my_channel(local_ch_id) = 1
     
     if ( local_ch_id > my_channel_high(iam) ) my_channel_high(iam) = local_ch_id
     if ( local_ch_id < my_channel_low(iam) ) my_channel_low(iam) = local_ch_id
     
  end do
  !write(6,*) my_channel_low, my_channel_high
  if(cc_approx .ne. 'mbpt2')then
     allocate( vnn_pppp(my_channel_low(iam):my_channel_high(iam)) )
  end if
end if

  !
  ! calculate work for hhpp block 
  !
  ltot_work=0
  tot_work = 0
  do channel   = 1, channels%number_hhpp_confs 
     nx = channels%hhpp_quantum_numbers(channel*4)
     ny = channels%hhpp_quantum_numbers(channel*4-1)
     nz = channels%hhpp_quantum_numbers(channel*4-2)
     tz = channels%hhpp_quantum_numbers(channel*4-3)
     
     dim1 = size(  lookup_hhpp_configs(1,channel)%ival, 2) 
     dim2 = size(  lookup_hhpp_configs(2,channel)%ival, 2) 
     
     tot_work = tot_work + dim1*dim2
     ltot_work = ltot_work + int(dim1,8)*int(dim2, 8)
     
  end do

  number_channels = channels%number_hhpp_confs 
  work_pr_proc = int( ltot_work/int(num_procs,8) )
  work = 0.d0
  
  allocate( mapping_hhpp(1:num_procs, 1:number_channels, 1:5) )
  mapping_hhpp = 0 
  if ( iam == 0 ) write(6,*) ltot_work, work_pr_proc, number_channels 

  curr_work = 0
  curr_proc = 0
  local_ch_id = 1
  number_channels = 0
  do channel   = 1, channels%number_hhpp_confs 
     nx = channels%hhpp_quantum_numbers(channel*4)
     ny = channels%hhpp_quantum_numbers(channel*4-1)
     nz = channels%hhpp_quantum_numbers(channel*4-2)
     tz = channels%hhpp_quantum_numbers(channel*4-3)
     
     dim1 = size(  lookup_hhpp_configs(1,channel)%ival, 2) 
     dim2 = size(  lookup_hhpp_configs(2,channel)%ival, 2) 
     
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
                    
              mapping_hhpp(curr_proc+1,local_ch_id,3 ) = i 
              
           else
              
              
              has_ch_been_added = .true. 
              mapping_hhpp(curr_proc+1,local_ch_id,1 ) = number_channels 
              mapping_hhpp(curr_proc+1,local_ch_id,2 ) = i
              mapping_hhpp(curr_proc+1,local_ch_id,3 ) = i
              mapping_hhpp(curr_proc+1,local_ch_id,4 ) = 1
              mapping_hhpp(curr_proc+1,local_ch_id,5 ) = bra_dim
              
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
  
  number_channels = channels%number_hhpp_confs 
  work = 0.d0
  do i = 1, num_procs
     do local_ch_id = 1, number_channels 
        if ( mapping_hhpp(i,local_ch_id,1) == 0 ) cycle
        
        row_dim = mapping_hhpp(i,local_ch_id,5)-mapping_hhpp(i,local_ch_id,4)+1
        column_dim = mapping_hhpp(i,local_ch_id,3)-mapping_hhpp(i,local_ch_id,2)+1
        work(i) = work(i) + real( row_dim*column_dim )
        
     end do
  end do
  
!!$  write(6,*);write(6,*);write(6,*)
!!$  
!!$  do i = 1, num_procs
!!$     
!!$     if ( iam ==0 ) write(6,*)
!!$     if ( iam ==0 ) write(6,'(a12,2x,i2,2x,a10,2x,(E12.6))') 'proc number', & 
!!$          i-1, '% work', work(i) *100.d0/ real( work_pr_proc ) 
!!$     do local_ch_id = 1, number_channels 
!!$        if ( mapping_hhpp(i,local_ch_id,1) == 0 ) cycle
!!$        
!!$        if ( iam ==0 ) write(6,*)
!!$        if ( iam ==0 ) write(6,'(a20,2x,i4)') 'Channel #:', mapping_hhpp(i,local_ch_id,1)
!!$        if ( iam ==0 ) write(6,'(a20,2x,i4)') 'Row index start:', mapping_hhpp(i,local_ch_id,4)
!!$        if ( iam ==0 ) write(6,'(a20,2x,i4)') 'Row index end:', mapping_hhpp(i,local_ch_id,5)
!!$        if ( iam ==0 ) write(6,'(a20,2x,i4)') 'Column index start:', mapping_hhpp(i,local_ch_id,2)
!!$        if ( iam ==0 ) write(6,'(a20,2x,i4)') 'Column index end:', mapping_hhpp(i,local_ch_id,3)
!!$        
!!$        
!!$     end do
!!$  end do
  
  call mpi_barrier(mpi_comm_world,ierror)
  number_channels = channels%number_hhpp_confs 

  allocate( my_hhpp_channel_low(0:num_procs-1), my_hhpp_channel_high(0:num_procs-1) )
  allocate( check_my_channel_hhpp(1:number_channels) )
  check_my_channel_hhpp = 0
  
  my_hhpp_channel_low(iam)  = number_channels + 1
  my_hhpp_channel_high(iam) = 0 
  do local_ch_id = 1, number_channels 
     if ( mapping_hhpp(iam+1,local_ch_id,1) == 0 ) cycle
     check_my_channel_hhpp(local_ch_id) = 1
      
     if ( local_ch_id > my_hhpp_channel_high(iam) ) my_hhpp_channel_high(iam) = local_ch_id
     if ( local_ch_id < my_hhpp_channel_low(iam) ) my_hhpp_channel_low(iam) = local_ch_id
     
  end do


  !
  ! calculate work for hphp block 
  !
  ltot_work=0
  tot_work = 0
  do channel   = 1, channels%number_hphp_confs 
     nx = channels%hphp_quantum_numbers(channel*4)
     ny = channels%hphp_quantum_numbers(channel*4-1)
     nz = channels%hphp_quantum_numbers(channel*4-2)
     tz = channels%hphp_quantum_numbers(channel*4-3)
     
     dim1 = size(  lookup_hphp_configs(1,channel)%ival, 2) 
     dim2 = dim1 
     
     tot_work = tot_work + dim1*dim2
     ltot_work = ltot_work + int(dim1,8)*int(dim2, 8)
     
  end do

  number_channels = channels%number_hphp_confs 
  work_pr_proc = int( ltot_work/int(num_procs,8) )
  work = 0.d0

  allocate( mapping_hphp(1:num_procs, 1:number_channels, 1:5) )
  mapping_hphp = 0 
  if ( iam == 0 ) write(6,*) ltot_work, work_pr_proc, number_channels 

  curr_work = 0
  curr_proc = 0
  local_ch_id = 1
  number_channels = 0
  do channel   = 1, channels%number_hphp_confs 
     nx = channels%hphp_quantum_numbers(channel*4)
     ny = channels%hphp_quantum_numbers(channel*4-1)
     nz = channels%hphp_quantum_numbers(channel*4-2)
     tz = channels%hphp_quantum_numbers(channel*4-3)
     
     dim1 = size(  lookup_hphp_configs(1,channel)%ival, 2) 
     dim2 = dim1
     
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
                    
              mapping_hphp(curr_proc+1,local_ch_id,3 ) = i 
              
           else
              
              
              has_ch_been_added = .true. 
              mapping_hphp(curr_proc+1,local_ch_id,1 ) = number_channels 
              mapping_hphp(curr_proc+1,local_ch_id,2 ) = i
              mapping_hphp(curr_proc+1,local_ch_id,3 ) = i
              mapping_hphp(curr_proc+1,local_ch_id,4 ) = 1
              mapping_hphp(curr_proc+1,local_ch_id,5 ) = bra_dim
              
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
  
  number_channels = channels%number_hphp_confs 
  work = 0.d0
  do i = 1, num_procs
     do local_ch_id = 1, number_channels 
        if ( mapping_hphp(i,local_ch_id,1) == 0 ) cycle
        
        row_dim = mapping_hphp(i,local_ch_id,5)-mapping_hphp(i,local_ch_id,4)+1
        column_dim = mapping_hphp(i,local_ch_id,3)-mapping_hphp(i,local_ch_id,2)+1
        work(i) = work(i) + real( row_dim*column_dim )
        
     end do
  end do
  
!!$  write(6,*);write(6,*);write(6,*)
!!$  
!!$  do i = 1, num_procs
!!$     
!!$     if ( iam ==0 ) write(6,*)
!!$     if ( iam ==0 ) write(6,'(a12,2x,i2,2x,a10,2x,(E12.6))') 'proc number', & 
!!$          i-1, '% work', work(i) *100.d0/ real( work_pr_proc ) 
!!$     do local_ch_id = 1, number_channels 
!!$        if ( mapping_hphp(i,local_ch_id,1) == 0 ) cycle
!!$        
!!$        if ( iam ==0 ) write(6,*)
!!$        if ( iam ==0 ) write(6,'(a20,2x,i4)') 'Channel #:', mapping_hphp(i,local_ch_id,1)
!!$        if ( iam ==0 ) write(6,'(a20,2x,i4)') 'Row index start:', mapping_hphp(i,local_ch_id,4)
!!$        if ( iam ==0 ) write(6,'(a20,2x,i4)') 'Row index end:', mapping_hphp(i,local_ch_id,5)
!!$        if ( iam ==0 ) write(6,'(a20,2x,i4)') 'Column index start:', mapping_hphp(i,local_ch_id,2)
!!$        if ( iam ==0 ) write(6,'(a20,2x,i4)') 'Column index end:', mapping_hphp(i,local_ch_id,3)
!!$        
!!$        
!!$     end do
!!$  end do
  
  call mpi_barrier(mpi_comm_world,ierror)
  number_channels = channels%number_hphp_confs 
  allocate( check_my_channel_hphp(1:number_channels) )
  check_my_channel_hphp = 0
  
  do local_ch_id = 1, number_channels 
     if ( mapping_hphp(iam+1,local_ch_id,1) == 0 ) cycle
     check_my_channel_hphp(local_ch_id) = 1
  end do

!!$
!!$  !
!!$  ! calculate work for hhhh block 
!!$  !
!!$  
!!$  ltot_work=0
!!$  tot_work = 0
!!$  do channel   = 1, channels%number_hhhh_confs 
!!$     nx = channels%hhhh_quantum_numbers(channel*4)
!!$     ny = channels%hhhh_quantum_numbers(channel*4-1)
!!$     nz = channels%hhhh_quantum_numbers(channel*4-2)
!!$     tz = channels%hhhh_quantum_numbers(channel*4-3)
!!$     
!!$     dim1 = size(  lookup_hhhh_configs(1,channel)%ival, 2) 
!!$     dim2 = dim1 
!!$     
!!$     tot_work = tot_work + dim1*dim2
!!$     ltot_work = ltot_work + int(dim1,8)*int(dim2, 8)
!!$     
!!$  end do
!!$
!!$  number_channels = channels%number_hhhh_confs 
!!$  work_pr_proc = int( ltot_work/int(num_procs,8) )
!!$  work = 0.d0
!!$
!!$  allocate( mapping_hhhh(1:num_procs, 1:number_channels, 1:5) )
!!$  mapping_hhhh = 0 
!!$  if ( iam == 0 ) write(6,*) ltot_work, work_pr_proc, number_channels 
!!$
!!$  curr_work = 0
!!$  curr_proc = 0
!!$  local_ch_id = 1
!!$  number_channels = 0
!!$  do channel   = 1, channels%number_hhhh_confs 
!!$     nx = channels%hhhh_quantum_numbers(channel*4)
!!$     ny = channels%hhhh_quantum_numbers(channel*4-1)
!!$     nz = channels%hhhh_quantum_numbers(channel*4-2)
!!$     tz = channels%hhhh_quantum_numbers(channel*4-3)
!!$     
!!$     dim1 = size(  lookup_hhhh_configs(1,channel)%ival, 2) 
!!$     dim2 = dim1
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
!!$              mapping_hhhh(curr_proc+1,local_ch_id,3 ) = i 
!!$              
!!$           else
!!$              
!!$              
!!$              has_ch_been_added = .true. 
!!$              mapping_hhhh(curr_proc+1,local_ch_id,1 ) = number_channels 
!!$              mapping_hhhh(curr_proc+1,local_ch_id,2 ) = i
!!$              mapping_hhhh(curr_proc+1,local_ch_id,3 ) = i
!!$              mapping_hhhh(curr_proc+1,local_ch_id,4 ) = 1
!!$              mapping_hhhh(curr_proc+1,local_ch_id,5 ) = bra_dim
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
!!$  number_channels = channels%number_hhhh_confs 
!!$  work = 0.d0
!!$  do i = 1, num_procs
!!$     do local_ch_id = 1, number_channels 
!!$        if ( mapping_hhhh(i,local_ch_id,1) == 0 ) cycle
!!$        
!!$        row_dim = mapping_hhhh(i,local_ch_id,5)-mapping_hhhh(i,local_ch_id,4)+1
!!$        column_dim = mapping_hhhh(i,local_ch_id,3)-mapping_hhhh(i,local_ch_id,2)+1
!!$        work(i) = work(i) + real( row_dim*column_dim )
!!$        
!!$     end do
!!$  end do
!!$  

!!$  
!!$  do i = 1, num_procs
!!$     
!!$     if ( iam ==0 ) write(6,*)
!!$     if ( iam ==0 ) write(6,'(a12,2x,i2,2x,a10,2x,(E12.6))') 'proc number', & 
!!$          i-1, '% work', work(i) *100.d0/ real( work_pr_proc ) 
!!$     do local_ch_id = 1, number_channels 
!!$        if ( mapping_hhhh(i,local_ch_id,1) == 0 ) cycle
!!$        
!!$        if ( iam ==0 ) write(6,*)
!!$        if ( iam ==0 ) write(6,'(a20,2x,i4)') 'Channel #:', mapping_hhhh(i,local_ch_id,1)
!!$        if ( iam ==0 ) write(6,'(a20,2x,i4)') 'Row index start:', mapping_hhhh(i,local_ch_id,4)
!!$        if ( iam ==0 ) write(6,'(a20,2x,i4)') 'Row index end:', mapping_hhhh(i,local_ch_id,5)
!!$        if ( iam ==0 ) write(6,'(a20,2x,i4)') 'Column index start:', mapping_hhhh(i,local_ch_id,2)
!!$        if ( iam ==0 ) write(6,'(a20,2x,i4)') 'Column index end:', mapping_hhhh(i,local_ch_id,3)
!!$        
!!$        
!!$     end do
!!$  end do
!!$  
!!$  call mpi_barrier(mpi_comm_world,ierror)
!!$  number_channels = channels%number_hhhh_confs 
!!$  allocate( check_my_channel_hhhh(1:number_channels) )
!!$  check_my_channel_hhhh = 0
!!$  
!!$  do local_ch_id = 1, number_channels 
!!$     if ( mapping_hhhh(iam+1,local_ch_id,1) == 0 ) cycle
!!$     check_my_channel_hhhh(local_ch_id) = 1
!!$  end do


  deallocate( work, memory )

  
end SUBROUTINE setup_proc_mappings_pp


SUBROUTINE setup_proc_mappings_hppp
  USE single_particle_orbits
  USE constants
  use one_body_operators
  use t2_storage
  use configurations
  USE KSPACE 
  use CHIRAL_POTENTIALS
  use parallel
  
  IMPLICIT NONE
  INTEGER :: channel, channel2, dim1, nx,ny,nz,sz,tz
  INTEGER :: total_work, work_per_proc, curr_work, curr_proc, curr_channel
  integer :: processor_work(num_procs), mtx_dim, mtx_dim_work, temp, tot_work, work_pr_proc, local_ch_id
  integer :: row_dim, column_dim, ket_dim, bra_dim
  INTEGER :: bra_min, bra_max, ket_min, ket_max, number_channels, i, dim2
  LOGICAL :: has_ch_been_added
  REAL*8, allocatable :: work(:), memory(:)
  integer(8):: ltot_work
  
  ltot_work=0
  tot_work = 0
  
  allocate( work(1:num_procs) )
  allocate( memory(1:num_procs) )
  work = 0.d0
  memory = 0.d0
  
  
  !
  ! calculate work for hppp block 
  !
  ltot_work=0
  tot_work = 0
  do channel   = 1, channels%number_hppp_confs 
     nx = channels%hppp_quantum_numbers(channel*4)
     ny = channels%hppp_quantum_numbers(channel*4-1)
     nz = channels%hppp_quantum_numbers(channel*4-2)
     tz = channels%hppp_quantum_numbers(channel*4-3)
     
     dim1 = size(  lookup_hppp_configs(1,channel)%ival, 2) 
     dim2 = size(  lookup_hppp_configs(2,channel)%ival, 2) 
     
     tot_work = tot_work + dim1*dim2
     ltot_work = ltot_work + int(dim1,8)*int(dim2, 8)
     
  end do

  number_channels = channels%number_hppp_confs 
  work_pr_proc = int( ltot_work/int(num_procs,8) )
  work = 0.d0
  
  allocate( mapping_hppp(1:num_procs, 1:number_channels, 1:5) )
  mapping_hppp = 0 
  if ( iam == 0 ) write(6,*) ltot_work, work_pr_proc, number_channels 

  curr_work = 0
  curr_proc = 0
  local_ch_id = 1
  number_channels = 0
  do channel   = 1, channels%number_hppp_confs 
     nx = channels%hppp_quantum_numbers(channel*4)
     ny = channels%hppp_quantum_numbers(channel*4-1)
     nz = channels%hppp_quantum_numbers(channel*4-2)
     tz = channels%hppp_quantum_numbers(channel*4-3)
     
     dim1 = size(  lookup_hppp_configs(1,channel)%ival, 2) 
     dim2 = size(  lookup_hppp_configs(2,channel)%ival, 2) 
     
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
                    
              mapping_hppp(curr_proc+1,local_ch_id,3 ) = i 
              
           else
              
              
              has_ch_been_added = .true. 
              mapping_hppp(curr_proc+1,local_ch_id,1 ) = number_channels 
              mapping_hppp(curr_proc+1,local_ch_id,2 ) = i
              mapping_hppp(curr_proc+1,local_ch_id,3 ) = i
              mapping_hppp(curr_proc+1,local_ch_id,4 ) = 1
              mapping_hppp(curr_proc+1,local_ch_id,5 ) = bra_dim
              
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
  
  number_channels = channels%number_hppp_confs 
  work = 0.d0
  do i = 1, num_procs
     do local_ch_id = 1, number_channels 
        if ( mapping_hppp(i,local_ch_id,1) == 0 ) cycle
        
        row_dim = mapping_hppp(i,local_ch_id,5)-mapping_hppp(i,local_ch_id,4)+1
        column_dim = mapping_hppp(i,local_ch_id,3)-mapping_hppp(i,local_ch_id,2)+1
        work(i) = work(i) + real( row_dim*column_dim )
        
     end do
  end do
  
!!$  write(6,*);write(6,*);write(6,*)
!!$  
!!$  do i = 1, num_procs
!!$     
!!$     if ( iam ==0 ) write(6,*)
!!$     if ( iam ==0 ) write(6,'(a12,2x,i2,2x,a10,2x,(E12.6))') 'proc number', & 
!!$          i-1, '% work', work(i) *100.d0/ real( work_pr_proc ) 
!!$     do local_ch_id = 1, number_channels 
!!$        if ( mapping_hppp(i,local_ch_id,1) == 0 ) cycle
!!$        
!!$        if ( iam ==0 ) write(6,*)
!!$        if ( iam ==0 ) write(6,'(a20,2x,i4)') 'Channel #:', mapping_hppp(i,local_ch_id,1)
!!$        if ( iam ==0 ) write(6,'(a20,2x,i4)') 'Row index start:', mapping_hppp(i,local_ch_id,4)
!!$        if ( iam ==0 ) write(6,'(a20,2x,i4)') 'Row index end:', mapping_hppp(i,local_ch_id,5)
!!$        if ( iam ==0 ) write(6,'(a20,2x,i4)') 'Column index start:', mapping_hppp(i,local_ch_id,2)
!!$        if ( iam ==0 ) write(6,'(a20,2x,i4)') 'Column index end:', mapping_hppp(i,local_ch_id,3)
!!$        
!!$        
!!$     end do
!!$  end do
  
  call mpi_barrier(mpi_comm_world,ierror)
  number_channels = channels%number_hppp_confs 
  allocate( check_my_channel_hppp(1:number_channels) )
  check_my_channel_hppp = 0
  
  do local_ch_id = 1, number_channels 
     if ( mapping_hppp(iam+1,local_ch_id,1) == 0 ) cycle
     check_my_channel_hppp(local_ch_id) = 1
  end do

  
  deallocate( work, memory )
end SUBROUTINE setup_proc_mappings_hppp

SUBROUTINE setup_proc_mappings_ph
  USE single_particle_orbits
  USE constants
  use one_body_operators
  use t2_storage
  use configurations
  USE KSPACE 
  use CHIRAL_POTENTIALS
  use parallel
  
  IMPLICIT NONE
  INTEGER :: channel, channel2, dim1, nx,ny,nz,sz,tz
  INTEGER :: total_work, work_per_proc, curr_work, curr_proc, curr_channel
  integer :: processor_work(num_procs), mtx_dim, mtx_dim_work, temp, tot_work, work_pr_proc, local_ch_id
  integer :: row_dim, column_dim, ket_dim, bra_dim
  INTEGER :: bra_min, bra_max, ket_min, ket_max, number_channels, i, dim2
  LOGICAL :: has_ch_been_added
  REAL*8, allocatable :: work(:), memory(:)
  integer(8):: ltot_work
  
  allocate( work(1:num_procs) )
  allocate( memory(1:num_procs) )
  work = 0.d0
  memory = 0.d0
  
  !
  ! calculate work for ph_hphp block 
  !
  ltot_work=0
  tot_work = 0
  do channel   = 1, channels%number_ph_hphp_confs 
     nx = channels%ph_hphp_quantum_numbers(channel*4)
     ny = channels%ph_hphp_quantum_numbers(channel*4-1)
     nz = channels%ph_hphp_quantum_numbers(channel*4-2)
     tz = channels%ph_hphp_quantum_numbers(channel*4-3)
     
     dim1 = size(  lookup_ph_hphp_configs(1,channel)%ival, 2) 
     dim2 = dim1 
     
     tot_work = tot_work + dim1*dim2
     ltot_work = ltot_work + int(dim1,8)*int(dim2, 8)
     
  end do

  number_channels = channels%number_ph_hphp_confs 
  work_pr_proc = int( ltot_work/int(num_procs,8) )
  work = 0.d0

  allocate( mapping_ph_hphp(1:num_procs, 1:number_channels, 1:5) )
  mapping_ph_hphp = 0 
  if ( iam == 0 ) write(6,*) ltot_work, work_pr_proc, number_channels 

  curr_work = 0
  curr_proc = 0
  local_ch_id = 1
  number_channels = 0
  do channel   = 1, channels%number_ph_hphp_confs 
     nx = channels%ph_hphp_quantum_numbers(channel*4)
     ny = channels%ph_hphp_quantum_numbers(channel*4-1)
     nz = channels%ph_hphp_quantum_numbers(channel*4-2)
     tz = channels%ph_hphp_quantum_numbers(channel*4-3)
     
     dim1 = size(  lookup_ph_hphp_configs(1,channel)%ival, 2) 
     dim2 = dim1
     
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
                    
              mapping_ph_hphp(curr_proc+1,local_ch_id,3 ) = i 
              
           else
              
              
              has_ch_been_added = .true. 
              mapping_ph_hphp(curr_proc+1,local_ch_id,1 ) = number_channels 
              mapping_ph_hphp(curr_proc+1,local_ch_id,2 ) = i
              mapping_ph_hphp(curr_proc+1,local_ch_id,3 ) = i
              mapping_ph_hphp(curr_proc+1,local_ch_id,4 ) = 1
              mapping_ph_hphp(curr_proc+1,local_ch_id,5 ) = bra_dim
              
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
  
  number_channels = channels%number_ph_hphp_confs 
  work = 0.d0
  do i = 1, num_procs
     do local_ch_id = 1, number_channels 
        if ( mapping_ph_hphp(i,local_ch_id,1) == 0 ) cycle
        
        row_dim = mapping_ph_hphp(i,local_ch_id,5)-mapping_ph_hphp(i,local_ch_id,4)+1
        column_dim = mapping_ph_hphp(i,local_ch_id,3)-mapping_ph_hphp(i,local_ch_id,2)+1
        work(i) = work(i) + real( row_dim*column_dim )
        
     end do
  end do
  
  
  !do i = 1, num_procs
     
  !   if ( iam ==0 ) write(6,*)
  !   if ( iam ==0 ) write(6,'(a12,2x,i2,2x,a10,2x,(E12.6))') 'proc number', & 
  !        i-1, '% work', work(i) *100.d0/ real( work_pr_proc ) 
  !   do local_ch_id = 1, number_channels 
  !      if ( mapping_ph_hphp(i,local_ch_id,1) == 0 ) cycle
  !      
  !      if ( iam ==0 ) write(6,*)
  !      if ( iam ==0 ) write(6,'(a20,2x,i4)') 'Channel #:', mapping_ph_hphp(i,local_ch_id,1)
  !      if ( iam ==0 ) write(6,'(a20,2x,i4)') 'Row index start:', mapping_ph_hphp(i,local_ch_id,4)
  !      if ( iam ==0 ) write(6,'(a20,2x,i4)') 'Row index end:', mapping_ph_hphp(i,local_ch_id,5)
  !      if ( iam ==0 ) write(6,'(a20,2x,i4)') 'Column index start:', mapping_ph_hphp(i,local_ch_id,2)
  !      if ( iam ==0 ) write(6,'(a20,2x,i4)') 'Column index end:', mapping_ph_hphp(i,local_ch_id,3)
  !      
  !      
  !   end do
  !end do
  
  call mpi_barrier(mpi_comm_world,ierror)
  number_channels = channels%number_ph_hphp_confs 
  allocate( check_my_channel_ph_hphp(1:number_channels) )
  check_my_channel_ph_hphp = 0
  

  allocate( my_ph_hphp_channel_low(0:num_procs-1), my_ph_hphp_channel_high(0:num_procs-1) )
  my_ph_hphp_channel_low(iam)  = number_channels + 1
  my_ph_hphp_channel_high(iam) = 0 
  do local_ch_id = 1, number_channels 
     if ( mapping_ph_hphp(iam+1,local_ch_id,1) == 0 ) cycle
     check_my_channel_ph_hphp(local_ch_id) = 1
     
     if ( local_ch_id > my_ph_hphp_channel_high(iam) ) my_ph_hphp_channel_high(iam) = local_ch_id
     if ( local_ch_id < my_ph_hphp_channel_low(iam) ) my_ph_hphp_channel_low(iam) = local_ch_id
     
  end do
  
  deallocate( work, memory )
  
end SUBROUTINE setup_proc_mappings_ph




SUBROUTINE setup_proc_t3_mappings
  USE single_particle_orbits
  USE constants
  use one_body_operators
  use t2_storage
  use configurations
  USE KSPACE 
  use CHIRAL_POTENTIALS
  use parallel
  
  IMPLICIT NONE
  INTEGER :: channel, channel2, dim1, nx,ny,nz,sz,tz
  INTEGER :: total_work, work_per_proc, curr_work, curr_proc, curr_channel
  integer :: processor_work(num_procs), mtx_dim, mtx_dim_work, temp, tot_work, work_pr_proc, local_ch_id
  integer :: row_dim, column_dim, ket_dim, bra_dim
  INTEGER :: bra_min, bra_max, ket_min, ket_max, number_channels, i, dim2
  LOGICAL :: has_ch_been_added
  REAL*8, allocatable :: work(:), memory(:)
  integer(8):: ltot_work
  
  !
  ! calculate work for hhpp block 
  !
  ltot_work=0
  tot_work = 0
  do channel   = 1, channels%number_t3_confs 
     nx = channels%t3_quantum_numbers(channel*4)
     ny = channels%t3_quantum_numbers(channel*4-1)
     nz = channels%t3_quantum_numbers(channel*4-2)
     tz = channels%t3_quantum_numbers(channel*4-3)
     
     dim1 = size(  lookup_t3_configs(1,channel)%ival,2 )
     dim2 = size(  lookup_t3_configs(2,channel)%ival,2 )
  
     
     tot_work = tot_work + dim1*dim2
     ltot_work = ltot_work + int(dim1,8)*int(dim2, 8)
     
  end do
  
  number_channels = channels%number_t3_confs 
  
  allocate( work(1:num_procs) )
  allocate( memory(1:num_procs) )
  work = 0.d0
  memory = 0.d0
  
  work_pr_proc = int( ltot_work/int(num_procs,8) )
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
     
     dim1 = size(  lookup_t3_configs(1,channel)%ival, 2) 
     dim2 = size(  lookup_t3_configs(2,channel)%ival, 2) 
     
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
  
  !write(6,*);write(6,*);write(6,*)
  
  do i = 1, num_procs

     if ( i /= 1 ) cycle 
     if ( iam ==0 ) write(6,*)
     if ( iam ==0 ) write(6,'(a12,2x,i2,2x,a10,2x,(E12.6))') 'proc number', & 
          i-1, '% work', work(i) *100.d0/ real( work_pr_proc ) 
     do local_ch_id = 1, number_channels 
        if ( mapping_t3(i,local_ch_id,1) == 0 ) cycle
        
        if ( iam ==0 ) write(6,*)
        if ( iam ==0 ) write(6,'(a20,2x,i4)') 'Channel #:', mapping_t3(i,local_ch_id,1)
        if ( iam ==0 ) write(6,'(a20,2x,i4)') 'Row index start:', mapping_t3(i,local_ch_id,4)
        if ( iam ==0 ) write(6,'(a20,2x,i4)') 'Row index end:', mapping_t3(i,local_ch_id,5)
        if ( iam ==0 ) write(6,'(a20,2x,i4)') 'Column index start:', mapping_t3(i,local_ch_id,2)
        if ( iam ==0 ) write(6,'(a20,2x,i4)') 'Column index end:', mapping_t3(i,local_ch_id,3)
        
        
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

  
end SUBROUTINE setup_proc_t3_mappings
