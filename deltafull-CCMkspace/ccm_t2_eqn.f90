!
! 
!
SUBROUTINE t2_intermediate(switch) 
  USE single_particle_orbits
  USE configurations
  USE constants
  use t2_storage
  use one_body_operators
  use parallel 
  use chiral_potentials 
  
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: switch
  INTEGER :: i,j,k,l,a,b,c,d, nx,ny,nz,sz,tz
  INTEGER :: bra,ket, channel, bra2,ket2, bra3,ket3,dim1, dim2, dim3, number_channels
  COMPLEX*16 ::  sum1, sum2 
  INTEGER :: nxd, nyd, nzd, szd, tzd, nxc, nyc, nzc, szc, tzc, nxk, nyk, nzk, szk, tzk
  COMPLEX*16, ALLOCATABLE :: temp1(:,:), temp2(:,:), temp(:,:), t3_tmp1(:,:), t3_tmp2(:,:),t2_tmp1(:), t2_tmp2(:), amat(:,:)
  real*8  ::  startwtime , endwtime
  INTEGER :: nx2, ny2, nz2, sz2, tz2, channel2, channel3, channel1
  INTEGER :: nx3, ny3, nz3, sz3, tz3, channel4, ndim, i1
  INTEGER :: bra_min, bra_max, ket_min, ket_max, local_ch_id
  !
  ! set t2 amps to zero
  !
  do channel = 1, channels%number_hhpp_confs 
     t2_ccm_eqn(channel)%val = 0.d0
  end do
  
  do channel   = 1, channels%number_ph_phhp_confs 
     t2_ph_ccm_eqn(channel)%val = 0.d0 
  end do
  !
  ! setup ph t-amplitudes 
  !
  do channel   = 1, channels%number_ph_phhp_confs 
     
     t2_ph_ccm(channel)%val = 0.d0 
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

!!$  call mpi_barrier(mpi_comm_world,ierror)
!!$  startwtime = MPI_WTIME()
  


  call setup_I6I7_intermediate



!!$  call mpi_barrier(mpi_comm_world,ierror)
!!$  endwtime = MPI_WTIME()
!!$  if ( iam == 0 ) write(6,*) 'Time for I6/I7', endwtime - startwtime
!!$  startwtime = MPI_WTIME()

  !
  ! setup dynamic denominators d2ab and d2ij
  !
  do i = 1, below_ef
     d2ij(i) = I7(i,i)
     I7(i,i) = 0.d0
  end do
  
  do a = below_ef+1, tot_orbs
     d2ab(a) = I6(a,a)
     I6(a,a) = 0.d0
  end do
!!$  
!!$  do a = below_ef+1, tot_orbs
!!$     do c = below_ef+1, tot_orbs
!!$        if ( abs(I6(a,c)) > 1.e-6 ) write(6,*) a,c,i6(a,c)
!!$     end do
!!$  end do
!!$  
!!$  write(6,*) 
!!$  write(6,*)
!!$  do i = 1, below_ef
!!$     do j = 1, below_ef
!!$        if ( abs(i7(i,j)) > 1.e-6 ) write(6,*) i,j,i7(i,j)
!!$     end do
!!$  end do
!!$  
  
  !do channel = 1, channels%number_hhpp_confs 
  do channel = my_hhpp_channel_low(iam), my_hhpp_channel_high(iam)
     nx = channels%hhpp_quantum_numbers(channel*4)
     ny = channels%hhpp_quantum_numbers(channel*4-1)
     nz = channels%hhpp_quantum_numbers(channel*4-2)
     tz = channels%hhpp_quantum_numbers(channel*4-3)
     
     !if ( check_my_channel_hhpp(channel) == 0 ) cycle
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
     
     do bra = bra_min, bra_max !1, size(  lookup_hhpp_configs(2,channel)%ival, 2) 
        a = lookup_hhpp_configs(2,channel)%ival(1,bra)
        b = lookup_hhpp_configs(2,channel)%ival(2,bra) 

        do ket = 1, size(  lookup_hhpp_configs(1,channel)%ival, 2) 
           i= lookup_hhpp_configs(1,channel)%ival(1,ket) 
           j = lookup_hhpp_configs(1,channel)%ival(2,ket) 
           
           t2_ccm_eqn(channel)%val(bra,ket) = conjg(vnn_hhpp(channel)%val(ket,bra)) 
                      
        end do
     end do
  end do

  !
  ! I6 and I7 contribution to t2
  !
  !do channel   = 1, channels%number_hhpp_confs 
  do channel = my_hhpp_channel_low(iam), my_hhpp_channel_high(iam)
     nx = channels%hhpp_quantum_numbers(channel*4)
     ny = channels%hhpp_quantum_numbers(channel*4-1)
     nz = channels%hhpp_quantum_numbers(channel*4-2)
     tz = channels%hhpp_quantum_numbers(channel*4-3)
     
     !if ( check_my_channel_hhpp(channel) == 0 ) cycle
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
     
     do bra = bra_min, bra_max !1, size(  lookup_hhpp_configs(2,channel)%ival, 2) 
        a = lookup_hhpp_configs(2,channel)%ival(1,bra)
        b = lookup_hhpp_configs(2,channel)%ival(2,bra) 

        do ket = 1, size(  lookup_hhpp_configs(1,channel)%ival, 2) 
           i= lookup_hhpp_configs(1,channel)%ival(1,ket) 
           j = lookup_hhpp_configs(1,channel)%ival(2,ket) 
           
           sum1 = 0.d0
           do szc = -1, 1 ,2 
              nxc = all_orbit%nx(b)
              nyc = all_orbit%ny(b)
              nzc = all_orbit%nz(b)
              tzc = all_orbit%itzp(b)
              c = orbitof(nxc, nyc, nzc, szc, tzc)
              if  ( c <= below_ef ) cycle 
              
              bra2 = pp_hhpp%ival(a,c) 
              if ( bra2 == 0 ) cycle 
              sum1 = sum1 + I6(b,c)*t2_ccm(channel)%val(bra2,ket)
           end do
           
           t2_ccm_eqn(channel)%val(bra,ket) = t2_ccm_eqn(channel)%val(bra,ket) + sum1 
           bra3 = pp_hhpp%ival(b,a) 
           t2_ccm_eqn(channel)%val(bra3,ket) = t2_ccm_eqn(channel)%val(bra3,ket) - sum1 
           
           sum1 = 0.d0
           do szk = -1, 1, 2
              nxk = all_orbit%nx(j)
              nyk = all_orbit%ny(j)
              nzk = all_orbit%nz(j)
              tzk = all_orbit%itzp(j)
              k = orbitof(nxk, nyk, nzk, szk, tzk)
              if  ( k > below_ef ) cycle 
              

              ket2 = hh_hhpp%ival(i,k) 
              if ( ket2 == 0 ) cycle 
              sum1 = sum1 - I7(k,j)*t2_ccm(channel)%val(bra,ket2)
           end do
           
           t2_ccm_eqn(channel)%val(bra,ket) = t2_ccm_eqn(channel)%val(bra,ket) + sum1 
           ket3 = hh_hhpp%ival(j,i)
           t2_ccm_eqn(channel)%val(bra,ket3) = t2_ccm_eqn(channel)%val(bra,ket3) - sum1 
           
           
        end do
     end do
  end do
  
!!$  call mpi_barrier(mpi_comm_world,ierror)
!!$  endwtime = MPI_WTIME()
!!$  if ( iam == 0 ) write(6,*) 'Time for I6/I7->t2', endwtime - startwtime
!!$  startwtime = MPI_WTIME()


  
!!$  number_channels = channels%number_hhhh_confs  
!!$  allocate( I9(number_channels) )
!!$  do channel   = 1, channels%number_hhhh_confs 
!!$     dim1 = size(  lookup_hhhh_configs(1,channel)%ival, 2) 
!!$     allocate( i9(channel)%val(dim1, dim1) ) 
!!$     i9(channel)%val = 0.d0 
!!$  end do
!!$  
!!$
!!$  
!!$  do channel   = 1, channels%number_hhhh_confs 
!!$     nx = channels%hhhh_quantum_numbers(channel*4)
!!$     ny = channels%hhhh_quantum_numbers(channel*4-1)
!!$     nz = channels%hhhh_quantum_numbers(channel*4-2)
!!$     tz = channels%hhhh_quantum_numbers(channel*4-3)
!!$     
!!$     do bra = 1, size(  lookup_hhhh_configs(1,channel)%ival, 2) 
!!$        k= lookup_hhhh_configs(1,channel)%ival(1,bra) 
!!$        l = lookup_hhhh_configs(1,channel)%ival(2,bra) 
!!$        ket3 = hh_hhpp%ival(k,l)
!!$        if ( ket3 == 0 ) cycle 
!!$        
!!$        do ket = 1, size(  lookup_hhhh_configs(1,channel)%ival, 2) 
!!$           i = lookup_hhhh_configs(1,channel)%ival(1,ket)
!!$           j = lookup_hhhh_configs(1,channel)%ival(2,ket) 
!!$           
!!$           !i9(channel)%val(bra,ket) = chiral_pot(k,l,i,j)
!!$           !i9(channel)%val(bra,ket) = vnn_hhhh(channel)%val(bra,ket) 
!!$           
!!$           nx2 = all_orbit%nx(i) + all_orbit%nx(j)
!!$           ny2 = all_orbit%ny(i) + all_orbit%ny(j)
!!$           nz2 = all_orbit%nz(i) + all_orbit%nz(j)
!!$           tz2 = (all_orbit%itzp(i) + all_orbit%itzp(j))/2
!!$           sz2 = (all_orbit%szp(i) + all_orbit%szp(j))/2
!!$           
!!$           channel2 = locate_channel(3,tz2, nx2, ny2, nz2) 
!!$           if ( channel2 == 0 ) cycle 
!!$           
!!$           if ( check_my_channel_hhpp(channel2) == 0 ) cycle
!!$           local_ch_id = channel2
!!$           !
!!$           ! ket side if fully stored on each proc
!!$           !
!!$           bra_min = mapping_hhpp(iam+1,local_ch_id,2)
!!$           bra_max = mapping_hhpp(iam+1,local_ch_id,3)
!!$                      
!!$           ket2 = hh_hhpp%ival(i,j)
!!$           sum1 = 0.d0 
!!$           do bra2 = bra_min, bra_max !1, size(  lookup_hhpp_configs(2,channel2)%ival, 2) 
!!$              c = lookup_hhpp_configs(2,channel2)%ival(1,bra2)
!!$              d = lookup_hhpp_configs(2,channel2)%ival(2,bra2) 
!!$              
!!$              !sum1 = sum1 + 0.5d0 * chiral_pot(k,l,c,d)*t2_ccm(channel2)%val(bra2,ket2)
!!$              sum1 = sum1 + 0.5d0 * vnn_hhpp(channel2)%val(ket3,bra2)*t2_ccm(channel2)%val(bra2,ket2)
!!$           end do
!!$           
!!$           i9(channel)%val(bra,ket) = i9(channel)%val(bra,ket) + sum1 
!!$           
!!$        end do
!!$     end do
!!$  end do
!!$
!!$  number_channels = size( i9 )
!!$  do channel = 1, number_channels 
!!$     bra = size(i9(channel)%val,1)
!!$     ket = size(i9(channel)%val,2)
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
!!$           t2_tmp1(i1) = i9(channel)%val(i,j)
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
!!$     i9(channel)%val = 0.d0
!!$     i1 = 0 
!!$     do i =1, bra
!!$        do j = 1, ket
!!$           i1 = i1 + 1
!!$           i9(channel)%val(i,j) = t2_tmp1(i1)
!!$        end do
!!$     end do
!!$     
!!$     deallocate( t2_tmp1, t2_tmp2 ) 
!!$  end do
!!$  
!!$  do channel   = 1, channels%number_hhhh_confs 
!!$     i9(channel)%val =i9(channel)%val + vnn_hhhh(channel)%val
!!$  end do
!!$  
!!$  !
!!$  ! I8 and I9 contribution to t2-amps
!!$  !
!!$
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
!!$     if ( check_my_channel(channel2) == 0 ) cycle
!!$     local_ch_id = channel2
!!$     
!!$     !
!!$     ! ket side if fully stored on each proc
!!$     !
!!$     ket_min = mapping(iam+1,local_ch_id,2)
!!$     ket_max = mapping(iam+1,local_ch_id,3)
!!$     !
!!$     ! bra side is distributed 
!!$     !
!!$     bra_min = mapping(iam+1,local_ch_id,4)
!!$     bra_max = mapping(iam+1,local_ch_id,5)
!!$     
!!$     dim1 = size( vnn_pppp(channel2)%val,1 ) 
!!$     dim2 = size( t2_ccm(channel2)%val, 2)
!!$     dim3 = size( vnn_pppp(channel2)%val, 2) 
!!$     
!!$     call ZGEMM ( 'n', 'n', dim1, dim2, dim3, DCMPLX(0.5d0,0.D0), vnn_pppp(channel2)%val(bra_min:bra_max,ket_min:ket_max), & 
!!$          dim1,t2_ccm(channel2)%val(ket_min:ket_max,:), dim3, & 
!!$          DCMPLX(1.d0,0.D0), t2_ccm_eqn(channel2)%val(bra_min:bra_max,:), dim1 )
!!$     
!!$  end do
!!$  !stop
!!$  do channel   = 1, channels%number_hhpp_confs 
!!$     nx = channels%hhpp_quantum_numbers(channel*4)
!!$     ny = channels%hhpp_quantum_numbers(channel*4-1)
!!$     nz = channels%hhpp_quantum_numbers(channel*4-2)
!!$     tz = channels%hhpp_quantum_numbers(channel*4-3)
!!$          
!!$     if ( check_my_channel_hhpp(channel) == 0 ) cycle
!!$     local_ch_id = channel
!!$     !
!!$     ! ket side if fully stored on each proc
!!$     !
!!$     bra_min = mapping_hhpp(iam+1,local_ch_id,2)
!!$     bra_max = mapping_hhpp(iam+1,local_ch_id,3)
!!$     !
!!$     ! bra side is distributed 
!!$     !
!!$     ket_min = mapping_hhpp(iam+1,local_ch_id,4)
!!$     ket_max = mapping_hhpp(iam+1,local_ch_id,5)
!!$     
!!$     do bra = bra_min, bra_max !1, size(  lookup_hhpp_configs(2,channel)%ival, 2) 
!!$        a = lookup_hhpp_configs(2,channel)%ival(1,bra)
!!$        b = lookup_hhpp_configs(2,channel)%ival(2,bra) 
!!$        
!!$        do ket = 1, size(  lookup_hhpp_configs(1,channel)%ival, 2) 
!!$           i= lookup_hhpp_configs(1,channel)%ival(1,ket) 
!!$           j = lookup_hhpp_configs(1,channel)%ival(2,ket) 
!!$                             
!!$
!!$           sum2 = 0.d0 
!!$           do ket2 = 1, size(  lookup_hhpp_configs(1,channel)%ival, 2) 
!!$              k = lookup_hhpp_configs(1,channel)%ival(1,ket2) 
!!$              l = lookup_hhpp_configs(1,channel)%ival(2,ket2) 
!!$              !sum2 = sum2 + 0.5d0* chiral_pot(k,l,i,j)*t2_ccm(channel)%val(bra,ket2)
!!$              
!!$              nx3 = all_orbit%nx(l) + all_orbit%nx(k)
!!$              ny3 = all_orbit%ny(l) + all_orbit%ny(k)
!!$              nz3 = all_orbit%nz(l) + all_orbit%nz(k)
!!$              tz3 = (all_orbit%itzp(l) + all_orbit%itzp(k))/2
!!$              channel4 = locate_channel(1,tz3, nx3, ny3, nz3) 
!!$              ket3 = hh_hhhh%ival(i,j)
!!$              bra3 = hh_hhhh%ival(k,l)
!!$              if ( channel4 == 0 ) cycle 
!!$              
!!$              !sum2 = sum2 + 0.5d0* chiral_pot(k,l,i,j)*t2_ccm(channel)%val(bra,ket2)
!!$              sum2 = sum2 + 0.5d0* i9(channel4)%val(bra3,ket3)*t2_ccm(channel)%val(bra,ket2)
!!$           end do
!!$           
!!$           t2_ccm_eqn(channel)%val(bra,ket) = t2_ccm_eqn(channel)%val(bra,ket) + sum2
!!$           
!!$        end do
!!$     end do
!!$  end do
!!$
!!$  endwtime = MPI_WTIME()

  
!!$  call mpi_barrier(mpi_comm_world,ierror)
!!$  startwtime = MPI_WTIME()

  number_channels = channels%number_hhhh_confs  
  allocate( I9(number_channels) )
  do channel   = 1, channels%number_hhhh_confs 
     dim1 = size(  lookup_hhhh_configs(1,channel)%ival, 2) 
     allocate( i9(channel)%val(dim1, dim1) ) 
     i9(channel)%val = 0.d0 
  end do
  
  do channel = my_hhpp_channel_low(iam), my_hhpp_channel_high(iam)
     !do channel   = 1, channels%number_hhpp_confs 
     nx = channels%hhpp_quantum_numbers(channel*4)
     ny = channels%hhpp_quantum_numbers(channel*4-1)
     nz = channels%hhpp_quantum_numbers(channel*4-2)
     tz = channels%hhpp_quantum_numbers(channel*4-3)
     
     channel2 = locate_channel(1,tz, nx, ny, nz) 
     if ( channel2 == 0 ) cycle 
     
     !if ( check_my_channel_hhpp(channel) == 0 ) cycle
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
     
     dim1 = size( vnn_hhpp(channel)%val, 1)
     dim2 = size( t2_ccm(channel)%val,2 ) 
     dim3 = bra_max-bra_min+1 
     
     ! Nov. 10
     call ZGEMM ( 'n', 'n', dim1, dim2, dim3, dcmplx(0.5d0,0.d0), vnn_hhpp(channel)%val(:,bra_min:bra_max), & 
          dim1,t2_ccm(channel)%val(bra_min:bra_max,:), dim3, & 
          dcmplx(1.d0,0.d0), I9(channel2)%val, dim1 )
     
     
  end do
  
!!$  call mpi_barrier(mpi_comm_world,ierror)
!!$  endwtime = MPI_WTIME()
!!$  if ( iam == 0 ) write(6,*) 'Time for I9a', endwtime - startwtime
!!$  startwtime = MPI_WTIME()
  
  number_channels = size( i9 )
  do channel = 1, number_channels 
     bra = size(i9(channel)%val,1)
     ket = size(i9(channel)%val,2)
     
     ndim = bra * ket
     allocate( t2_tmp1(ndim) )
     allocate( t2_tmp2(ndim) )
     t2_tmp1 = 0.d0
     t2_tmp2 = 0.d0

     i1 = 0 
     do i =1, bra
        do j = 1, ket
           i1 = i1 + 1
           t2_tmp1(i1) = i9(channel)%val(i,j)
        end do
     end do
     
     call mpi_allreduce(t2_tmp1,t2_tmp2,ndim,mpi_complex16,mpi_sum, &
          mpi_comm_world,ierror)
     
     !call mpi_reduce(t2_tmp1,t2_tmp2,ndim,mpi_complex16,mpi_sum, &
     !     master,mpi_comm_world,ierror)
     t2_tmp1 = 0.d0 
     t2_tmp1 = t2_tmp2 
     !call mpi_bcast(t2_tmp1,ndim,mpi_complex16,master, &
     !     mpi_comm_world,ierror)
     
     i9(channel)%val = 0.d0
     i1 = 0 
     do i =1, bra
        do j = 1, ket
           i1 = i1 + 1
           i9(channel)%val(i,j) = t2_tmp1(i1)
        end do
     end do
     
     deallocate( t2_tmp1, t2_tmp2 ) 
  end do
!!$
!!$  number_channels = size( i9 )
!!$  ndim = 0 
!!$  do channel = 1, number_channels 
!!$     bra = size(i9(channel)%val,1)
!!$     ket = size(i9(channel)%val,2)
!!$     
!!$     ndim = ndim +  bra * ket
!!$  end do
!!$  
!!$  allocate( t2_tmp1(ndim) )
!!$  allocate( t2_tmp2(ndim) )
!!$  t2_tmp1 = 0.d0
!!$  t2_tmp2 = 0.d0
!!$  
!!$  i1 = 0 
!!$  do channel = 1, number_channels 
!!$     bra = size(i9(channel)%val,1)
!!$     ket = size(i9(channel)%val,2)
!!$     
!!$     do i =1, bra
!!$        do j = 1, ket
!!$           i1 = i1 + 1
!!$           t2_tmp1(i1) = i9(channel)%val(i,j)
!!$        end do
!!$     end do
!!$  end do
!!$  
!!$  call mpi_reduce(t2_tmp1,t2_tmp2,ndim,mpi_complex16,mpi_sum, &
!!$       master,mpi_comm_world,ierror)
!!$  t2_tmp1 = 0.d0 
!!$  t2_tmp1 = t2_tmp2 
!!$  call mpi_bcast(t2_tmp1,ndim,mpi_complex16,master, &
!!$       mpi_comm_world,ierror)
!!$  
!!$  i1 = 0 
!!$  do channel = 1, number_channels 
!!$     i9(channel)%val = 0.d0
!!$     do i =1, bra
!!$        do j = 1, ket
!!$           i1 = i1 + 1
!!$           i9(channel)%val(i,j) = t2_tmp1(i1)
!!$        end do
!!$     end do
!!$  end do
!!$  deallocate( t2_tmp1, t2_tmp2 ) 
  
  do channel   = 1, channels%number_hhhh_confs 
     i9(channel)%val =i9(channel)%val + vnn_hhhh(channel)%val
  end do



!!$  call mpi_barrier(mpi_comm_world,ierror)
!!$  endwtime = MPI_WTIME()
!!$  if ( iam == 0 ) write(6,*) 'Time for I9b', endwtime - startwtime
!!$  startwtime = MPI_WTIME()


  !
  ! I8 and I9 contribution to t2-amps
  !
  
  do channel   = my_channel_low(iam), my_channel_high(iam) !1, channels%number_pppp_confs 
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
     
     dim1 = size( vnn_pppp(channel2)%val,1 ) 
     dim2 = size( t2_ccm(channel2)%val, 2)
     dim3 = size( vnn_pppp(channel2)%val, 2) 
     
     call ZGEMM ( 'n', 'n', dim1, dim2, dim3, DCMPLX(0.5d0,0.D0), vnn_pppp(channel2)%val(bra_min:bra_max,ket_min:ket_max), & 
          dim1,t2_ccm(channel2)%val(ket_min:ket_max,:), dim3, & 
          DCMPLX(1.d0,0.D0), t2_ccm_eqn(channel2)%val(bra_min:bra_max,:), dim1 )
     
  end do

!!$  call mpi_barrier(mpi_comm_world,ierror)
!!$  endwtime = MPI_WTIME()
!!$  if ( iam == 0 ) write(6,*) 'Time for I8->t2', endwtime - startwtime
!!$  startwtime = MPI_WTIME()

  do channel = my_hhpp_channel_low(iam), my_hhpp_channel_high(iam)
     !do channel   = 1, channels%number_hhpp_confs 
     nx = channels%hhpp_quantum_numbers(channel*4)
     ny = channels%hhpp_quantum_numbers(channel*4-1)
     nz = channels%hhpp_quantum_numbers(channel*4-2)
     tz = channels%hhpp_quantum_numbers(channel*4-3)
     
     channel2 = locate_channel(1,tz, nx, ny, nz) 
     if ( channel2 == 0 ) cycle 
     
     !if ( check_my_channel_hhpp(channel) == 0 ) cycle
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
     
     dim1 = bra_max-bra_min + 1 
     dim2 = size( I9(channel2)%val,2 )
     dim3 = size( t2_ccm(channel)%val,2) 
     
     ! Nov. 10
     call ZGEMM ( 'n', 'n', dim1, dim2, dim3, dcmplx(0.5d0,0.d0), t2_ccm(channel)%val(bra_min:bra_max,:), & 
          dim1,I9(channel2)%val, dim3, & 
          dcmplx(1.d0,0.d0), t2_ccm_eqn(channel)%val(bra_min:bra_max,:), dim1 )
     
     
  end do

!!$  call mpi_barrier(mpi_comm_world,ierror)
!!$  endwtime = MPI_WTIME()
!!$  if ( iam == 0 ) write(6,*) 'Time for I9->t2', endwtime - startwtime
!!$  startwtime = MPI_WTIME()



  
  number_channels = channels%number_ph_hphp_confs
  allocate( I4_ph(number_channels) )
  do channel   = 1, channels%number_ph_hphp_confs 
     dim1 = size(  lookup_ph_hphp_configs(1,channel)%ival, 2) 
     allocate( i4_ph(channel)%val(dim1, dim1) ) 
     i4_ph(channel)%val = 0.d0 
  end do
  
  !do channel   = 1, channels%number_ph_hphp_confs 
  do channel   = my_ph_hphp_channel_low(iam), my_ph_hphp_channel_high(iam) 
     nx = channels%ph_hphp_quantum_numbers(channel*4)
     ny = channels%ph_hphp_quantum_numbers(channel*4-1)
     nz = channels%ph_hphp_quantum_numbers(channel*4-2)
     tz = channels%ph_hphp_quantum_numbers(channel*4-3)
     
     !if ( check_my_channel_ph_hphp(channel) == 0 ) cycle
     local_ch_id = channel
     !
     ! ket side if fully stored on each proc
     !
     bra_min = mapping_ph_hphp(iam+1,local_ch_id,2)
     bra_max = mapping_ph_hphp(iam+1,local_ch_id,3)
     !
     ! bra side is distributed 
     !
     ket_min = mapping_ph_hphp(iam+1,local_ch_id,4)
     ket_max = mapping_ph_hphp(iam+1,local_ch_id,5)
     
     channel1 = locate_ph_channel(1,tz, nx, ny, nz)
     channel3 = locate_ph_channel(3,tz, nx, ny, nz)
     
     if ( channel1*channel3 == 0 ) cycle 
     
     dim1 = size(  lookup_ph_hpph_configs(1,channel3)%ival, 2) 
     dim2 = size(  lookup_ph_hpph_configs(2,channel3)%ival, 2) 
     allocate( amat(dim1, dim2) )
     amat = 0.d0 
     do bra = 1, size(  lookup_ph_hpph_configs(1,channel3)%ival, 2) 
        i = lookup_ph_hpph_configs(1,channel3)%ival(1,bra) 
        a = lookup_ph_hpph_configs(1,channel3)%ival(2,bra) 
        
        do ket = 1, size(  lookup_ph_hpph_configs(2,channel3)%ival, 2) 
           b = lookup_ph_hpph_configs(2,channel3)%ival(1,ket)
           j = lookup_ph_hpph_configs(2,channel3)%ival(2,ket) 
           
           nx2 = all_orbit%nx(i) + all_orbit%nx(j)
           ny2 = all_orbit%ny(i) + all_orbit%ny(j)
           nz2 = all_orbit%nz(i) + all_orbit%nz(j)
           tz2 = (all_orbit%itzp(i) + all_orbit%itzp(j))/2
           channel2 = locate_channel(3,tz2, nx2, ny2, nz2) 
           if ( channel2 == 0 ) cycle 
           
           bra2 = hh_hhpp%ival(i,j)
           ket2 = pp_hhpp%ival(a,b)
           if ( bra2*ket2 == 0 ) cycle 
           
           amat(bra,ket) = vnn_hhpp(channel2)%val(bra2,ket2)
           
        end do
     end do
     
     dim1 = size( amat,1 )
     dim2 = bra_max-bra_min + 1 !size( t2_ph_ccm(channel1)%val, 2)
     dim3 = size( amat,2 )
     
     call ZGEMM ( 'n', 'n', dim1, dim2, dim3, dcmplx(0.5d0,0.d0), amat, & 
          dim1,t2_ph_ccm(channel1)%val(:,bra_min:bra_max), dim3, & 
          dcmplx(1.d0,0.d0), I4_ph(channel)%val(:,bra_min:bra_max), dim1 )
     
     !call ZGEMM ( 'n', 'n', dim1, dim2, dim3, dcmplx(0.5d0,0.d0), vnn_ph_hpph(channel3)%val, & 
     !     dim1,t2_ph_ccm(channel1)%val(:,bra_min:bra_max), dim3, & 
     !     dcmplx(1.d0,0.d0), I4_ph(channel)%val(:,bra_min:bra_max), dim1 )
     
     deallocate( amat ) 
     
  end do
  
!!$  call mpi_barrier(mpi_comm_world,ierror)
!!$  endwtime = MPI_WTIME()
!!$  if ( iam == 0 ) write(6,*) 'Time for I4a', endwtime - startwtime
!!$  startwtime = MPI_WTIME()

  
  do channel   = my_ph_hphp_channel_low(iam), my_ph_hphp_channel_high(iam) 
     !do channel   = 1, channels%number_ph_hphp_confs 
     
     !if ( check_my_channel_ph_hphp(channel) == 0 ) cycle
     local_ch_id = channel
     !
     ! ket side if fully stored on each proc
     !
     bra_min = mapping_ph_hphp(iam+1,local_ch_id,2)
     bra_max = mapping_ph_hphp(iam+1,local_ch_id,3)
     !
     ! bra side is distributed 
     !
     ket_min = mapping_ph_hphp(iam+1,local_ch_id,4)
     ket_max = mapping_ph_hphp(iam+1,local_ch_id,5)
     
     do ket = bra_min, bra_max !1, size(  lookup_ph_hphp_configs(1,channel)%ival, 2) 
        j = lookup_ph_hphp_configs(1,channel)%ival(1,ket)
        b = lookup_ph_hphp_configs(1,channel)%ival(2,ket) 
     
        do bra = 1, size(  lookup_ph_hphp_configs(1,channel)%ival, 2) 
           i = lookup_ph_hphp_configs(1,channel)%ival(1,bra) 
           a = lookup_ph_hphp_configs(1,channel)%ival(2,bra) 
        
           
           nx2 = all_orbit%nx(i) + all_orbit%nx(b)
           ny2 = all_orbit%ny(i) + all_orbit%ny(b)
           nz2 = all_orbit%nz(i) + all_orbit%nz(b)
           tz2 = (all_orbit%itzp(i) + all_orbit%itzp(b))/2
           channel2 = locate_channel(4,tz2, nx2, ny2, nz2) 
           if ( channel2 == 0 ) cycle 
           
           ket2 = hp_hphp%ival(i,b)
           bra2 = hp_hphp%ival(j,a)
           if ( bra2*ket2 == 0 ) cycle 
           
           I4_ph(channel)%val(bra,ket) = I4_ph(channel)%val(bra,ket) - vnn_hphp(channel2)%val(ket2,bra2)  

        end do
     end do
  end do

  if ( switch == 1 ) then 
     call mpi_barrier(mpi_comm_world,ierror)
     !if ( iam == 0 ) write(6,*) 'Time for I4b', endwtime - startwtime
     startwtime = MPI_WTIME()
     
     ! add linear terms in T2 and V3nf 
     if ( tnf_approx == 3 ) call t2diags_v3t2_linear
     call mpi_barrier(mpi_comm_world,ierror)
     endwtime = MPI_WTIME()
     if ( iam == 0 ) write(6,*) 'Time for V3T2 diags', endwtime - startwtime

  end if
!!$  startwtime = MPI_WTIME()



!!$
!!$  number_channels = size( i4_ph )
!!$  do channel = 1, number_channels 
!!$     bra = size(i4_ph(channel)%val,1)
!!$     ket = size(i4_ph(channel)%val,2)
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
!!$           t2_tmp1(i1) = i4_ph(channel)%val(i,j)
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
!!$     i4_ph(channel)%val = 0.d0
!!$     i1 = 0 
!!$     do i =1, bra
!!$        do j = 1, ket
!!$           i1 = i1 + 1
!!$           i4_ph(channel)%val(i,j) = t2_tmp1(i1)
!!$        end do
!!$     end do
!!$     
!!$     deallocate( t2_tmp1, t2_tmp2 ) 
!!$  end do
!!$  
!!$  call mpi_barrier(mpi_comm_world,ierror)
!!$  endwtime = MPI_WTIME()
!!$  if ( iam == 0 ) write(6,*) 'Time for I4c', endwtime - startwtime
!!$  startwtime = MPI_WTIME()

!!$   number_channels = channels%number_hphp_confs  
!!$  allocate( I4(number_channels) )
!!$  do channel   = 1, channels%number_hphp_confs 
!!$     dim1 = size(  lookup_hphp_configs(1,channel)%ival, 2) 
!!$     allocate( i4(channel)%val(dim1, dim1) ) 
!!$     i4(channel)%val = 0.d0 
!!$  end do
!!$
!!$  do channel   = 1, channels%number_hhpp_confs 
!!$     nx = channels%hhpp_quantum_numbers(channel*4)
!!$     ny = channels%hhpp_quantum_numbers(channel*4-1)
!!$     nz = channels%hhpp_quantum_numbers(channel*4-2)
!!$     tz = channels%hhpp_quantum_numbers(channel*4-3)
!!$     
!!$     if ( check_my_channel_hhpp(channel) == 0 ) cycle
!!$     local_ch_id = channel
!!$     !
!!$     ! ket side if fully stored on each proc
!!$     !
!!$     bra_min = mapping_hhpp(iam+1,local_ch_id,2)
!!$     bra_max = mapping_hhpp(iam+1,local_ch_id,3)
!!$     !
!!$     ! bra side is distributed 
!!$     !
!!$     ket_min = mapping_hhpp(iam+1,local_ch_id,4)
!!$     ket_max = mapping_hhpp(iam+1,local_ch_id,5)
!!$     
!!$     do bra = bra_min, bra_max !1, size(  lookup_hhpp_configs(2,channel)%ival, 2) 
!!$        d = lookup_hhpp_configs(2,channel)%ival(1,bra)
!!$        b = lookup_hhpp_configs(2,channel)%ival(2,bra) 
!!$        
!!$        do ket = 1, size(  lookup_hhpp_configs(1,channel)%ival, 2) 
!!$           l= lookup_hhpp_configs(1,channel)%ival(1,ket) 
!!$           j = lookup_hhpp_configs(1,channel)%ival(2,ket) 
!!$  
!!$           do k = 1, below_ef
!!$              do szc = -1, 1, 2
!!$                 nxc = all_orbit%nx(k) + all_orbit%nx(b) - all_orbit%nx(j)
!!$                 nyc = all_orbit%ny(k) + all_orbit%ny(b) - all_orbit%ny(j) 
!!$                 nzc = all_orbit%nz(k) + all_orbit%nz(b) - all_orbit%nz(j) 
!!$                 tzc = all_orbit%itzp(k) + all_orbit%itzp(b) - all_orbit%itzp(j)
!!$                 if ( abs(nxc) > nkxmax ) cycle 
!!$                 if ( abs(nyc) > nkymax ) cycle 
!!$                 if ( abs(nzc) > nkzmax ) cycle 
!!$                 if ( abs(tzc) > 1 ) cycle 
!!$                 c = orbitof(nxc, nyc, nzc, szc, tzc)
!!$                 if  ( c <= below_ef ) cycle 
!!$                 
!!$                 nx2 = all_orbit%nx(k) + all_orbit%nx(b)
!!$                 ny2 = all_orbit%ny(k) + all_orbit%ny(b)
!!$                 nz2 = all_orbit%nz(k) + all_orbit%nz(b)
!!$                 tz2 = (all_orbit%itzp(k) + all_orbit%itzp(b))/2
!!$                 channel2 = locate_channel(4,tz2, nx2, ny2, nz2) 
!!$                 if ( channel2 == 0 ) cycle 
!!$                            
!!$                 ket2 = hp_hphp%ival(k,b)
!!$                 bra2 = hp_hphp%ival(j,c)
!!$              
!!$                 nx3 = all_orbit%nx(d) + all_orbit%nx(c)
!!$                 ny3 = all_orbit%ny(d) + all_orbit%ny(c)
!!$                 nz3 = all_orbit%nz(d) + all_orbit%nz(c)
!!$                 tz3 = (all_orbit%itzp(d) + all_orbit%itzp(c))/2
!!$                 channel3 = locate_channel(3,tz3, nx3, ny3, nz3) 
!!$                 if ( channel3 == 0 ) cycle 
!!$                 ket3 = hh_hhpp%ival(k,l)
!!$                 bra3 = pp_hhpp%ival(c,d)
!!$
!!$                 i4(channel2)%val(ket2,bra2) = i4(channel2)%val(ket2,bra2) -  & 
!!$                      0.5d0 * vnn_hhpp(channel3)%val(ket3,bra3)*t2_ccm(channel)%val(bra,ket)
!!$                 
!!$              end do
!!$           end do
!!$        end do
!!$     end do
!!$  end do
!!$
!!$  number_channels = size( i4 )
!!$  do channel = 1, number_channels 
!!$     bra = size(i4(channel)%val,1)
!!$     ket = size(i4(channel)%val,2)
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
!!$           t2_tmp1(i1) = i4(channel)%val(i,j)
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
!!$     i4(channel)%val = 0.d0
!!$     i1 = 0 
!!$     do i =1, bra
!!$        do j = 1, ket
!!$           i1 = i1 + 1
!!$           i4(channel)%val(i,j) = t2_tmp1(i1)
!!$        end do
!!$     end do
!!$     
!!$     deallocate( t2_tmp1, t2_tmp2 ) 
!!$  end do
!!$  
!!$  
!!$  do channel   = 1, channels%number_hphp_confs 
!!$     i4(channel)%val = i4(channel)%val + vnn_hphp(channel)%val
!!$  end do
!!$  
!!$
!!$  do channel   = 1, channels%number_ph_hphp_confs 
!!$     nx = channels%ph_hphp_quantum_numbers(channel*4)
!!$     ny = channels%ph_hphp_quantum_numbers(channel*4-1)
!!$     nz = channels%ph_hphp_quantum_numbers(channel*4-2)
!!$     tz = channels%ph_hphp_quantum_numbers(channel*4-3)
!!$     
!!$     dim1 = size(  lookup_ph_hphp_configs(1,channel)%ival, 2) 
!!$     do bra = 1, size(  lookup_ph_hphp_configs(1,channel)%ival, 2) 
!!$        i = lookup_ph_hphp_configs(1,channel)%ival(1,bra) 
!!$        a = lookup_ph_hphp_configs(1,channel)%ival(2,bra) 
!!$        
!!$        do ket = 1, size(  lookup_ph_hphp_configs(1,channel)%ival, 2) 
!!$           j = lookup_ph_hphp_configs(1,channel)%ival(1,ket)
!!$           b = lookup_ph_hphp_configs(1,channel)%ival(2,ket) 
!!$           
!!$           nx2 = all_orbit%nx(i) + all_orbit%nx(b)
!!$           ny2 = all_orbit%ny(i) + all_orbit%ny(b)
!!$           nz2 = all_orbit%nz(i) + all_orbit%nz(b)
!!$           tz2 = (all_orbit%itzp(i) + all_orbit%itzp(b))/2
!!$           channel2 = locate_channel(4,tz2, nx2, ny2, nz2) 
!!$           if ( channel2 == 0 ) cycle 
!!$           
!!$           ket2 = hp_hphp%ival(i,b)
!!$           bra2 = hp_hphp%ival(j,a)
!!$           if ( bra2*ket2 == 0 ) cycle 
!!$           
!!$           if ( abs( I4(channel2)%val(ket2,bra2) + I4_ph(channel)%val(bra,ket) ) > 1.e-10 ) write(6,*) & 
!!$                'error i4', real(-I4(channel2)%val(ket2,bra2)), real(I4_ph(channel)%val(bra,ket))
!!$        end do
!!$     end do
!!$  end do
!!$  !stop
  
  do channel   = my_ph_hphp_channel_low(iam), my_ph_hphp_channel_high(iam) 
     !do channel   = 1, channels%number_ph_hphp_confs 
     nx = channels%ph_hphp_quantum_numbers(channel*4)
     ny = channels%ph_hphp_quantum_numbers(channel*4-1)
     nz = channels%ph_hphp_quantum_numbers(channel*4-2)
     tz = channels%ph_hphp_quantum_numbers(channel*4-3)
     
     !if ( check_my_channel_ph_hphp(channel) == 0 ) cycle
     local_ch_id = channel
     !
     ! ket side if fully stored on each proc
     !
     bra_min = mapping_ph_hphp(iam+1,local_ch_id,2)
     bra_max = mapping_ph_hphp(iam+1,local_ch_id,3)
     !
     ! bra side is distributed 
     !
     ket_min = mapping_ph_hphp(iam+1,local_ch_id,4)
     ket_max = mapping_ph_hphp(iam+1,local_ch_id,5)
     
     channel2 = locate_ph_channel(1,tz, nx, ny, nz)
     if ( channel2 == 0 ) cycle 
     
     dim1 = size( t2_ph_ccm(channel2)%val,1 )
     dim2 = bra_max-bra_min+1 
     dim3 = size( t2_ph_ccm(channel2)%val, 2) 
     
     call ZGEMM ( 'n', 'n', dim1, dim2, dim3, dcmplx(1.d0,0.d0), t2_ph_ccm(channel2)%val, & 
          dim1,I4_ph(channel)%val(:,bra_min:bra_max), dim3, & 
          dcmplx(1.d0,0.d0), t2_ph_ccm_eqn(channel2)%val(:,bra_min:bra_max), dim1 )
     
  end do

!!$  call mpi_barrier(mpi_comm_world,ierror)
!!$  endwtime = MPI_WTIME()
!!$  if ( iam == 0 ) write(6,*) 'Time for I4d', endwtime - startwtime
!!$  startwtime = MPI_WTIME()

!!$  number_channels = size( t2_ph_ccm_eqn )
!!$  do channel = 1, number_channels 
!!$     bra = size(t2_ph_ccm_eqn(channel)%val,1)
!!$     ket = size(t2_ph_ccm_eqn(channel)%val,2)
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
!!$           t2_tmp1(i1) = t2_ph_ccm_eqn(channel)%val(i,j)
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
!!$     t2_ph_ccm_eqn(channel)%val = 0.d0
!!$     i1 = 0 
!!$     do i =1, bra
!!$        do j = 1, ket
!!$           i1 = i1 + 1
!!$           t2_ph_ccm_eqn(channel)%val(i,j) = t2_tmp1(i1)
!!$        end do
!!$     end do
!!$     
!!$     deallocate( t2_tmp1, t2_tmp2 ) 
!!$  end do
!!$
!!$  call mpi_barrier(mpi_comm_world,ierror)
!!$  endwtime = MPI_WTIME()
!!$  if ( iam == 0 ) write(6,*) 'Time for I4e', endwtime - startwtime
!!$  startwtime = MPI_WTIME()

  !
  ! I4 contribution to t2-amps
  !
  !do channel   = 1, channels%number_hhpp_confs 
!!$  do channel = my_hhpp_channel_low(iam), my_hhpp_channel_high(iam)
!!$     nx = channels%hhpp_quantum_numbers(channel*4)
!!$     ny = channels%hhpp_quantum_numbers(channel*4-1)
!!$     nz = channels%hhpp_quantum_numbers(channel*4-2)
!!$     tz = channels%hhpp_quantum_numbers(channel*4-3)
!!$     
!!$     !if ( check_my_channel_hhpp(channel) == 0 ) cycle
!!$     local_ch_id = channel
!!$     !
!!$     ! ket side if fully stored on each proc
!!$     !
!!$     bra_min = mapping_hhpp(iam+1,local_ch_id,2)
!!$     bra_max = mapping_hhpp(iam+1,local_ch_id,3)
!!$     !
!!$     ! bra side is distributed 
!!$     !
!!$     ket_min = mapping_hhpp(iam+1,local_ch_id,4)
!!$     ket_max = mapping_hhpp(iam+1,local_ch_id,5)
!!$     
!!$     do bra = bra_min, bra_max 
!!$        a = lookup_hhpp_configs(2,channel)%ival(1,bra)
!!$        b = lookup_hhpp_configs(2,channel)%ival(2,bra) 
!!$        
!!$        do ket = 1, size(  lookup_hhpp_configs(1,channel)%ival, 2) 
!!$           i = lookup_hhpp_configs(1,channel)%ival(1,ket) 
!!$           j = lookup_hhpp_configs(1,channel)%ival(2,ket) 
!!$           
!!$           nx2 = all_orbit%nx(a) - all_orbit%nx(i)
!!$           ny2 = all_orbit%ny(a) - all_orbit%ny(i)
!!$           nz2 = all_orbit%nz(a) - all_orbit%nz(i)
!!$           tz2 = (all_orbit%itzp(a) - all_orbit%itzp(i))/2
!!$           channel2 = locate_ph_channel(1,tz2, nx2, ny2, nz2) 
!!$           if ( channel2 == 0 ) cycle 
!!$           
!!$           bra2 = ph_ph_phhp%ival(a,i)
!!$           ket2 = hp_ph_phhp%ival(j,b)
!!$           if ( bra2*ket2 == 0 ) cycle 
!!$           
!!$           sum2 = t2_ph_ccm_eqn(channel2)%val(bra2,ket2)
!!$           t2_ccm_eqn(channel)%val(bra,ket) = t2_ccm_eqn(channel)%val(bra,ket) + sum2 
!!$           bra2 = pp_hhpp%ival(b,a)
!!$           t2_ccm_eqn(channel)%val(bra2,ket) = t2_ccm_eqn(channel)%val(bra2,ket) - sum2 
!!$           ket2 = hh_hhpp%ival(j,i)
!!$           t2_ccm_eqn(channel)%val(bra,ket2) = t2_ccm_eqn(channel)%val(bra,ket2) - sum2 
!!$           t2_ccm_eqn(channel)%val(bra2,ket2) = t2_ccm_eqn(channel)%val(bra2,ket2) + sum2 
!!$           
!!$        end do
!!$     end do
!!$  end do


  

  do channel   = my_ph_hphp_channel_low(iam), my_ph_hphp_channel_high(iam) 
     !do channel   = 1, channels%number_ph_hphp_confs 
     nx = channels%ph_hphp_quantum_numbers(channel*4)
     ny = channels%ph_hphp_quantum_numbers(channel*4-1)
     nz = channels%ph_hphp_quantum_numbers(channel*4-2)
     tz = channels%ph_hphp_quantum_numbers(channel*4-3)
     
     channel2 = locate_ph_channel(1,tz, nx, ny, nz)
     !if ( channel2 /= channel ) write(6,*) 'iam', iam, 'Error',  channel, channel2 
     if ( channel2 == 0 ) cycle 
     local_ch_id = channel
     
     !
     ! ket side if fully stored on each proc
     !
     bra_min = mapping_ph_hphp(iam+1,local_ch_id,2)
     bra_max = mapping_ph_hphp(iam+1,local_ch_id,3)
     !
     ! bra side is distributed 
     !
     ket_min = mapping_ph_hphp(iam+1,local_ch_id,4)
     ket_max = mapping_ph_hphp(iam+1,local_ch_id,5)
     
     do ket = bra_min, bra_max !1, size(  lookup_ph_hphp_configs(1,channel)%ival, 2) 
        j = lookup_ph_phhp_configs(2,channel2)%ival(1,ket)
        b = lookup_ph_phhp_configs(2,channel2)%ival(2,ket) 
        
        do bra = 1, size(  lookup_ph_phhp_configs(1,channel2)%ival, 2) 
           a = lookup_ph_phhp_configs(1,channel2)%ival(1,bra) 
           i = lookup_ph_phhp_configs(1,channel2)%ival(2,bra) 
        
           nx2 = all_orbit%nx(i) + all_orbit%nx(j)
           ny2 = all_orbit%ny(i) + all_orbit%ny(j)
           nz2 = all_orbit%nz(i) + all_orbit%nz(j)
           tz2 = (all_orbit%itzp(i) + all_orbit%itzp(j))/2
           channel1 = locate_channel(3,tz2, nx2, ny2, nz2) 
           if ( channel1 == 0 ) cycle 
           
           ket2 = hh_hhpp%ival(i,j)
           bra2 = pp_hhpp%ival(a,b)
           if ( bra2*ket2 == 0 ) cycle 
           
           sum2 = t2_ph_ccm_eqn(channel2)%val(bra,ket)
           t2_ccm_eqn(channel1)%val(bra2,ket2) = t2_ccm_eqn(channel1)%val(bra2,ket2) + sum2 
           bra3 = pp_hhpp%ival(b,a)
           t2_ccm_eqn(channel1)%val(bra3,ket2) = t2_ccm_eqn(channel1)%val(bra3,ket2) - sum2 
           ket3 = hh_hhpp%ival(j,i)
           t2_ccm_eqn(channel1)%val(bra2,ket3) = t2_ccm_eqn(channel1)%val(bra2,ket3) - sum2 
           t2_ccm_eqn(channel1)%val(bra3,ket3) = t2_ccm_eqn(channel1)%val(bra3,ket3) + sum2 
           
        end do
     end do
  end do
  
  !if ( cc_approx == 'CCDT1' ) CALL t2diags_v2t3_linear
  
  
  

!!$     
!!$
!!$  call mpi_barrier(mpi_comm_world,ierror)
!!$  
  
  

  !endwtime = MPI_WTIME()
  !if ( iam == 0 ) write(6,*) 'Time for I4->t2', endwtime - startwtime
!!$  startwtime = MPI_WTIME()
  
  number_channels = size( t2_ccm )
  do channel = 1, number_channels 
     bra = size(t2_ccm(channel)%val,1)
     ket = size(t2_ccm(channel)%val,2)
     
     ndim = bra * ket
     allocate( t2_tmp1(ndim) )
     allocate( t2_tmp2(ndim) )
     t2_tmp1 = 0.d0
     t2_tmp2 = 0.d0

     i1 = 0 
     do i =1, bra
        do j = 1, ket
           i1 = i1 + 1
           t2_tmp1(i1) = t2_ccm_eqn(channel)%val(i,j)
        end do
     end do
     
     call mpi_allreduce(t2_tmp1,t2_tmp2,ndim,mpi_complex16,mpi_sum, &
          mpi_comm_world,ierror)
     
     !call mpi_reduce(t2_tmp1,t2_tmp2,ndim,mpi_complex16,mpi_sum, &
     !     master,mpi_comm_world,ierror)
     t2_tmp1 = 0.d0 
     t2_tmp1 = t2_tmp2 
     !call mpi_bcast(t2_tmp1,ndim,mpi_complex16,master, &
     !     mpi_comm_world,ierror)
     
     i1 = 0 
     do i =1, bra
        do j = 1, ket
           i1 = i1 + 1
           t2_ccm_eqn(channel)%val(i,j) = t2_tmp1(i1)
        end do
     end do
     
     deallocate( t2_tmp1, t2_tmp2 ) 
  end do

!!$  call mpi_barrier(mpi_comm_world,ierror)
!!$  endwtime = MPI_WTIME()
!!$  if ( iam == 0 ) write(6,*) 'Time for t2 recuce 1', endwtime - startwtime
!!$  startwtime = MPI_WTIME()
!!$  
!!$  
!!$  number_channels = size( t2_ccm )
!!$  do channel = 1, number_channels 
!!$     bra = size(t2_ccm(channel)%val,1)
!!$     ket = size(t2_ccm(channel)%val,2)
!!$     
!!$     do j = 1, ket
!!$        ndim = bra 
!!$        
!!$        allocate( t2_tmp1(ndim) )
!!$        allocate( t2_tmp2(ndim) )
!!$        t2_tmp1 = 0.d0
!!$        t2_tmp2 = 0.d0
!!$
!!$        do i =1, bra
!!$           t2_tmp1(i) = t2_ccm_eqn(channel)%val(i,j)
!!$        end do
!!$        
!!$        call mpi_reduce(t2_tmp1,t2_tmp2,ndim,mpi_complex16,mpi_sum, &
!!$             master,mpi_comm_world,ierror)
!!$        t2_tmp1 = 0.d0 
!!$        t2_tmp1 = t2_tmp2 
!!$        call mpi_bcast(t2_tmp1,ndim,mpi_complex16,master, &
!!$             mpi_comm_world,ierror)
!!$     
!!$        do i =1, bra
!!$           t2_ccm_eqn(channel)%val(i,j) = t2_tmp1(i)
!!$        end do
!!$        deallocate( t2_tmp1, t2_tmp2 ) 
!!$        
!!$        
!!$     end do
!!$  end do
!!$  call mpi_barrier(mpi_comm_world,ierror)
!!$  endwtime = MPI_WTIME()
!!$  if ( iam == 0 ) write(6,*) 'Time for t2 recuce 2', endwtime - startwtime
  
  deallocate( i6, i7 ) 
  
  do i = 1, size(i4_ph) 
     deallocate(i4_ph(i)%val )
  ENDDO
  deallocate( i4_ph ) 

  do i = 1, size(i9)
     deallocate(i9(i)%val)
  end do
  deallocate( i9 )
  
    
!!$  do i = 1, size(i4) 
!!$     deallocate(i4(i)%val )
!!$  ENDDO
!!$  deallocate( i4 ) 
!!$  

end  SUBROUTINE t2_intermediate




!
!
!
SUBROUTINE setup_I6I7_intermediate
  
  USE single_particle_orbits
  USE configurations
  USE constants
  use t2_storage
  use one_body_operators
  use chiral_potentials
  
  IMPLICIT NONE
  
  INTEGER :: i,a,i1, c,d,k,l,bra,bra2,ket, ndim2, nxc, nyc,nzc, tzc,szc
  COMPLEX*16 :: gmat,t2, sum1
  INTEGER :: bra_min, bra_max, ket_min, ket_max, local_ch_id
  COMPLEX*16, allocatable :: t2_tmp1(:), t2_tmp2(:) 
  real*8 :: endwtime, startwtime
  !
  ! Setup I6ac intermediate
  !
  allocate( I6ac(below_ef+1:tot_orbs, below_ef+1:tot_orbs) ) 
  allocate( I6(below_ef+1:tot_orbs, below_ef+1:tot_orbs) ) 
  allocate( I7( 1:below_ef, 1:below_ef) )
  
  
  !
  ! I5b  intermediate
  !
  i6 = 0.d0 
  I7 = 0.d0
  I6ac = 0.d0

!!$
!!$  call mpi_barrier(mpi_comm_world,ierror)
!!$  startwtime = MPI_WTIME()
!!$  
!!$
!!$  
  !
  ! I7(b) intermediate
  !
  do i1 = my_hhpp_channel_low(iam), my_hhpp_channel_high(iam)
     !do i1 = 1, channels%number_hhpp_confs 
     
     !if ( check_my_channel_hhpp(i1) == 0 ) cycle
     local_ch_id = i1
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

     DO ket= bra_min, bra_max !1, size(  lookup_hhpp_configs(2,i1)%ival, 2) 
        c = lookup_hhpp_configs(2,i1)%ival(1,ket) 
        d = lookup_hhpp_configs(2,i1)%ival(2,ket) 
        
        DO bra=1, size(  lookup_hhpp_configs(1,i1)%ival, 2) 
           i = lookup_hhpp_configs(1,i1)%ival(1,bra) 
           l = lookup_hhpp_configs(1,i1)%ival(2,bra) 
           
           do szc = -1, 1 ,2 
              nxc = all_orbit%nx(i)
              nyc = all_orbit%ny(i)
              nzc = all_orbit%nz(i)
              tzc = all_orbit%itzp(i)
              k = orbitof(nxc, nyc, nzc, szc, tzc)
              if  ( k > below_ef ) cycle 
           
              bra2 = hh_hhpp%ival(k,l)
              if ( bra2 == 0 ) cycle 
              gmat = vnn_hhpp(i1)%val(bra2,ket) 
              i7(k,i) =  i7(k,i) + 0.5d0 * gmat*t2_ccm(i1)%val(ket,bra) 
              
           end do
        end DO
     end DO
  end do


  !
  ! I6(abc) intermediate
  ! Here I6(b) is calculated, I6(a+c) calculated in t1-eqn
  !
  do i1 = my_hhpp_channel_low(iam), my_hhpp_channel_high(iam)
     !do i1 = 1, channels%number_hhpp_confs 
     
     !if ( check_my_channel_hhpp(i1) == 0 ) cycle
     local_ch_id = i1
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
     
     DO bra= bra_min, bra_max !1, size( lookup_hhpp_configs(2,i1)%ival, 2) 
        d = lookup_hhpp_configs(2,i1)%ival(1,bra) 
        a = lookup_hhpp_configs(2,i1)%ival(2,bra) 
        
        DO ket=1, size(  lookup_hhpp_configs(1,i1)%ival, 2) 
           k = lookup_hhpp_configs(1,i1)%ival(1,ket) 
           l = lookup_hhpp_configs(1,i1)%ival(2,ket)  
           
           t2 = t2_ccm(i1)%val(bra,ket) 
           do szc = -1, 1 ,2 
              nxc = all_orbit%nx(a)
              nyc = all_orbit%ny(a)
              nzc = all_orbit%nz(a)
              tzc = all_orbit%itzp(a)
              c = orbitof(nxc, nyc, nzc, szc, tzc)
              if  ( c <= below_ef ) cycle 
              
              bra2 = pp_hhpp%ival(d,c)
              if ( bra2 == 0 ) cycle 
              
              gmat = vnn_hhpp(i1)%val(ket,bra2)
              sum1 = - 0.5d0* gmat * t2 
              I6(a,c) = I6(a,c) + sum1
           end do
        end DO
     end DO
  end do

!!$  call mpi_barrier(mpi_comm_world,ierror)
!!$  endwtime = MPI_WTIME()
!!$  if ( iam == 0 ) write(6,*) 'Time for I6/I7 a', endwtime - startwtime
!!$  startwtime = MPI_WTIME()

  
  ndim2 = 0 
  do a = 1, tot_orbs
     do c = 1, tot_orbs
        
        if ( all_orbit%nx(a) /= all_orbit%nx(c) ) cycle 
        if ( all_orbit%ny(a) /= all_orbit%ny(c) ) cycle
        if ( all_orbit%nz(a) /= all_orbit%nz(c) ) cycle
        if ( all_orbit%itzp(a) /= all_orbit%itzp(c) ) cycle
        
        ndim2 = ndim2 + 1 
     END DO
  end do
 
  
  allocate( t2_tmp1(ndim2) )
  allocate( t2_tmp2(ndim2) )
  t2_tmp1 = 0.d0
  t2_tmp2 = 0.d0
  
  i1 = 0
  do a = 1, tot_orbs
     do c = 1, tot_orbs
        if ( all_orbit%nx(a) /= all_orbit%nx(c) ) cycle 
        if ( all_orbit%ny(a) /= all_orbit%ny(c) ) cycle
        if ( all_orbit%nz(a) /= all_orbit%nz(c) ) cycle
        if ( all_orbit%itzp(a) /= all_orbit%itzp(c) ) cycle
        
        i1 = i1 + 1
        if ( a > below_ef .and. c > below_ef ) t2_tmp1(i1) = I6(a,c)
        if ( a <= below_ef .and. c <= below_ef ) t2_tmp1(i1) = I7(a,c)
     end do
  end do
  
  call mpi_reduce(t2_tmp1,t2_tmp2,ndim2,mpi_complex16,mpi_sum, &
       master,mpi_comm_world,ierror)
  t2_tmp1 = 0.d0
  t2_tmp1 = t2_tmp2
  call mpi_bcast(t2_tmp1,ndim2,mpi_complex16,master, &
       mpi_comm_world,ierror)

  i6 = 0.d0
  I7 = 0.d0
  i1 = 0
  do a = 1, tot_orbs
     do c = 1, tot_orbs
        if ( all_orbit%nx(a) /= all_orbit%nx(c) ) cycle 
        if ( all_orbit%ny(a) /= all_orbit%ny(c) ) cycle
        if ( all_orbit%nz(a) /= all_orbit%nz(c) ) cycle
        if ( all_orbit%itzp(a) /= all_orbit%itzp(c) ) cycle
                
        i1 = i1 + 1
        if ( a > below_ef .and. c > below_ef )   I6(a,c) = t2_tmp1(i1)
        if ( a <= below_ef .and. c <= below_ef ) I7(a,c) = t2_tmp1(i1)
     end do
  end do
  deallocate( t2_tmp1, t2_tmp2 )
  
  
  I7 = I7 + fock_mtx(1:below_ef, 1:below_ef) 
  I6ac = fock_mtx(below_ef+1:tot_orbs, below_ef+1:tot_orbs) 
  
  I6 = i6 + I6ac

!!$  call mpi_barrier(mpi_comm_world,ierror)
!!$  endwtime = MPI_WTIME()
!!$  if ( iam == 0 ) write(6,*) 'Time for I6/I7 b', endwtime - startwtime

  deallocate( i6ac )

end SUBROUTINE setup_I6I7_intermediate

subroutine t2diags_v2t3_linear
  USE single_particle_orbits
  USE constants
  use one_body_operators
  use t2_storage
  use configurations
  USE KSPACE 
  use CHIRAL_POTENTIALS
  use t3_constants
  use parallel
  
  IMPLICIT NONE
  INTEGER :: i,j,k,l, sumx, sumy, sumz, a,b,c,nxd,nyd,nzd,tzd,szd,tz,sz,nx,ny,nz
  INTEGER :: local_ch_id, dim1, dim2, bra_min, bra_max, ket_min, ket_max, bra, ket
  INTEGER :: channel, d, bra2, ket2, nx2, ny2, nz2,tz2, channel2, nx3, ny3, nz3,tz3,channel5, bra3,ket3
  COMPLEX*16 :: SUM1 

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
     
     

     do bra = bra_min, bra_max  !1, size(lookup_t3_configs(2,channel)%ival,2) ! 
        a = lookup_t3_configs(2,channel)%ival(1,bra)
        c = lookup_t3_configs(2,channel)%ival(2,bra) 
        d = lookup_t3_configs(2,channel)%ival(3,bra) 
        
        if ( c > d ) cycle 
        
        nx3 = all_orbit%nx(d) + all_orbit%nx(c)
        ny3 = all_orbit%ny(d) + all_orbit%ny(c)
        nz3 = all_orbit%nz(d) + all_orbit%nz(c)
        tz3 = (all_orbit%itzp(d) + all_orbit%itzp(c))/2
        channel5 = locate_channel(5,tz3, nx3, ny3, nz3) 
        if ( channel5 == 0 ) cycle

        ket3 = pp_hppp%ival(c,d)
        if ( ket3 == 0 ) cycle
     
        do ket = 1, size(  lookup_t3_configs(1,channel)%ival, 2) 
           i = lookup_t3_configs(1,channel)%ival(1,ket) 
           j = lookup_t3_configs(1,channel)%ival(2,ket) 
           k = lookup_t3_configs(1,channel)%ival(3,ket) 
           
           nx2 = all_orbit%nx(i) + all_orbit%nx(j)
           ny2 = all_orbit%ny(i) + all_orbit%ny(j)
           nz2 = all_orbit%nz(i) + all_orbit%nz(j)
           tz2 = (all_orbit%itzp(i) + all_orbit%itzp(j))/2
           channel2 = locate_channel(3,tz2, nx2, ny2, nz2) 
           
           if ( channel2 == 0 ) cycle
           ket2 = hh_hhpp%ival(i,j)
           if ( ket2 == 0 ) cycle
           
           nxd = all_orbit%nx(i) + all_orbit%nx(j) - all_orbit%nx(a)
           nyd = all_orbit%ny(i) + all_orbit%ny(j) - all_orbit%ny(a) 
           nzd = all_orbit%nz(i) + all_orbit%nz(j) - all_orbit%nz(a) 
           tzd = all_orbit%itzp(i) + all_orbit%itzp(j) - all_orbit%itzp(a)
           
           if ( abs(nxd) > nkxmax ) cycle
           if ( abs(nyd) > nkymax ) cycle
           if ( abs(nzd) > nkzmax ) cycle
           if ( abs(tzd) > 1 ) cycle
           
           do szd = -1, 1 ,2 
              
              b = orbitof(nxd, nyd, nzd, szd, tzd)
              if ( b <= below_ef ) cycle 
              
              bra2 = pp_hhpp%ival(a,b)
              if ( bra2 == 0 ) cycle
              bra3 = hp_hppp%ival(k,b)
              if ( bra3 == 0 ) cycle
              
              sum1 = -vnn_hppp(channel5)%val(bra3,ket3)*t3_ccm(channel)%val(bra,ket)

              t2_ccm_eqn(channel2)%val(bra2,ket2) = t2_ccm_eqn(channel2)%val(bra2,ket2) + sum1 
              bra2 = pp_hhpp%ival(b,a)
              t2_ccm_eqn(channel2)%val(bra2,ket2) = t2_ccm_eqn(channel2)%val(bra2,ket2) - sum1 
              
           end do
           
        end do
     end do

  end do



  
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
     
     

     do bra = bra_min, bra_max  !1, size(lookup_t3_configs(2,channel)%ival,2) ! 
        a = lookup_t3_configs(2,channel)%ival(1,bra)
        b = lookup_t3_configs(2,channel)%ival(2,bra) 
        c = lookup_t3_configs(2,channel)%ival(3,bra) 
        
        nx2 = all_orbit%nx(a) + all_orbit%nx(b)
        ny2 = all_orbit%ny(a) + all_orbit%ny(b)
        nz2 = all_orbit%nz(a) + all_orbit%nz(b)
        tz2 = (all_orbit%itzp(a) + all_orbit%itzp(b))/2
        channel2 = locate_channel(3,tz2, nx2, ny2, nz2) 
        if ( channel2 == 0 ) cycle
        bra2 = pp_hhpp%ival(a,b)
        if ( bra2 == 0 ) cycle
              
        !if ( c > d ) cycle 
        
        do ket = 1, size(  lookup_t3_configs(1,channel)%ival, 2) 
           i = lookup_t3_configs(1,channel)%ival(1,ket) 
           k = lookup_t3_configs(1,channel)%ival(2,ket) 
           l = lookup_t3_configs(1,channel)%ival(3,ket) 
           
           if ( k > l ) cycle 
           
           nx3 = all_orbit%nx(l) + all_orbit%nx(k)
           ny3 = all_orbit%ny(l) + all_orbit%ny(k)
           nz3 = all_orbit%nz(l) + all_orbit%nz(k)
           tz3 = (all_orbit%itzp(l) + all_orbit%itzp(k))/2
           channel5 = locate_channel(2,tz3, nx3, ny3, nz3) 
           if ( channel5 == 0 ) cycle
           bra3 = hh_hhhp%ival(k,l)
           if ( bra3 == 0 ) cycle
           
           nxd = all_orbit%nx(a) + all_orbit%nx(b) - all_orbit%nx(i)
           nyd = all_orbit%ny(a) + all_orbit%ny(b) - all_orbit%ny(i) 
           nzd = all_orbit%nz(a) + all_orbit%nz(b) - all_orbit%nz(i) 
           tzd = all_orbit%itzp(a) + all_orbit%itzp(b) - all_orbit%itzp(i)
           
           if ( abs(nxd) > nkxmax ) cycle
           if ( abs(nyd) > nkymax ) cycle
           if ( abs(nzd) > nkzmax ) cycle
           if ( abs(tzd) > 1 ) cycle
           
           do szd = -1, 1 ,2 
              
              j = orbitof(nxd, nyd, nzd, szd, tzd)
              if ( j > below_ef ) cycle 
              
              ket2 = hh_hhpp%ival(i,j)
              if ( ket2 == 0 ) cycle
              ket3 = hp_hhhp%ival(j,c)
              if ( ket3 == 0 ) cycle
              
              sum1 = -vnn_hhhp(channel5)%val(bra3,ket3)*t3_ccm(channel)%val(bra,ket)
              
              t2_ccm_eqn(channel2)%val(bra2,ket2) = t2_ccm_eqn(channel2)%val(bra2,ket2) + sum1 
              ket2 = hh_hhpp%ival(j,i)
              t2_ccm_eqn(channel2)%val(bra2,ket2) = t2_ccm_eqn(channel2)%val(bra2,ket2) - sum1 
              
           end do
           
        end do
     end do

  end do
end subroutine t2diags_v2t3_linear



SUBROUTINE t2diags_v3t2_linear2
  USE single_particle_orbits
  USE configurations
  USE constants
  use t2_storage
  use one_body_operators
  use parallel 
  use chiral_potentials 
  
  IMPLICIT NONE
  INTEGER :: i,j,k,l,m,a,b,c,d,e, nx,ny,nz,sz,tz
  INTEGER :: bra,ket, ket4,channel, bra2,ket2, bra3,ket3,dim1, dim2, dim3, number_channels
  COMPLEX*16 ::  sum1, sum2, t2dum 
  INTEGER :: nxd, nyd, nzd, szd, tzd,nxe, nye, nze, sze, tze, nxc, nyc, nzc, szc, tzc, nxk, nyk, nzk, szk, tzk
  COMPLEX*16, ALLOCATABLE :: temp1(:,:), temp2(:,:), temp(:,:), t3_tmp1(:,:), t3_tmp2(:,:),t2_tmp1(:), t2_tmp2(:), amat(:,:)
  real*8  ::  startwtime , endwtime
  INTEGER :: nx2, ny2, nz2, sz2, tz2, channel2, channel3, channel1
  INTEGER :: nx3, ny3, nz3, sz3, tz3, channel4, ndim, i1
  INTEGER :: bra_min, bra_max, ket_min, ket_max, local_ch_id, bra4
  complex*16, allocatable :: x_ae(:,:), x_mi(:,:)
  integer, allocatable :: ijc_conf(:,:,:), ij_conf(:,:), ab_conf(:,:)
  
  allocate( ab_conf(below_ef+1:tot_orbs, below_ef+1:tot_orbs) )
  allocate( ij_conf(1:below_ef, 1:below_ef) )
  ab_conf = 0 
  ij_conf = 0 
  
  call mpi_barrier(mpi_comm_world,ierror)
  startwtime = MPI_WTIME()
  
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
!!$     do bra = bra_min, bra_max 
!!$        k = lookup_ph_hpphpp_configs(1,channel)%ival(1,bra)
!!$        a = lookup_ph_hpphpp_configs(1,channel)%ival(2,bra) 
!!$        b = lookup_ph_hpphpp_configs(1,channel)%ival(3,bra) 
!!$        i = lookup_ph_hpphpp_configs(1,channel)%ival(4,bra) 
!!$        
!!$        nx2 = all_orbit%nx(a) + all_orbit%nx(b)
!!$        ny2 = all_orbit%ny(a) + all_orbit%ny(b)
!!$        nz2 = all_orbit%nz(a) + all_orbit%nz(b)
!!$        tz2 = (all_orbit%itzp(a) + all_orbit%itzp(b))/2
!!$        channel2 = locate_channel(3,tz2, nx2, ny2, nz2) 
!!$        if ( channel2 == 0 ) cycle 
!!$           
!!$        bra2 = pp_hhpp%ival(a,b)
!!$        if ( bra2 == 0 ) cycle 
!!$
!!$        do ket = 1, size(  lookup_ph_hpphpp_configs(2,channel)%ival, 2) 
!!$           c = lookup_ph_hpphpp_configs(2,channel)%ival(1,ket) 
!!$           d = lookup_ph_hpphpp_configs(2,channel)%ival(2,ket) 
!!$           
!!$           bra3 = pp_hhpp%ival(c,d)
!!$           if ( bra3 == 0 ) cycle 
!!$           
!!$           nxd = all_orbit%nx(a) + all_orbit%nx(b) - all_orbit%nx(i)
!!$           nyd = all_orbit%ny(a) + all_orbit%ny(b) - all_orbit%ny(i) 
!!$           nzd = all_orbit%nz(a) + all_orbit%nz(b) - all_orbit%nz(i) 
!!$           tzd = all_orbit%itzp(a) + all_orbit%itzp(b) - all_orbit%itzp(i)
!!$           
!!$           if ( abs(nxd) > nkxmax ) cycle
!!$           if ( abs(nyd) > nkymax ) cycle
!!$           if ( abs(nzd) > nkzmax ) cycle
!!$           if ( abs(tzd) > 1 ) cycle
!!$           
!!$           do szd = -1, 1 ,2 
!!$              
!!$              j = orbitof(nxd, nyd, nzd, szd, tzd)
!!$              if ( j > below_ef ) cycle 
!!$                            
!!$              ket2 = hh_hhpp%ival(i,j)
!!$              if ( ket2 == 0) cycle 
!!$              
!!$              ket3 = hh_hhpp%ival(k,j)
!!$              if ( ket3 == 0) cycle 
!!$              
!!$              
!!$              sum1 = -t2_ccm(channel)%val(bra3,ket3)*v3nf_ph_hpphpp(channel)%val(bra,ket)
!!$              t2_ccm_eqn(channel2)%val(bra2,ket2)  = t2_ccm_eqn(channel2)%val(bra2,ket2) + sum1 
!!$              ket2 = hh_hhpp%ival(j,i)
!!$              t2_ccm_eqn(channel2)%val(bra2,ket2)  = t2_ccm_eqn(channel2)%val(bra2,ket2) - sum1 
!!$              
!!$           end do
!!$        end do
!!$     end do
!!$  end do
!!$  
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
     
     
     do bra = bra_min, bra_max 
        k = lookup_hpphpp_configs(1,channel)%ival(1,bra)
        a = lookup_hpphpp_configs(1,channel)%ival(2,bra) 
        b = lookup_hpphpp_configs(1,channel)%ival(3,bra) 
        
        nx2 = all_orbit%nx(a) + all_orbit%nx(b)
        ny2 = all_orbit%ny(a) + all_orbit%ny(b)
        nz2 = all_orbit%nz(a) + all_orbit%nz(b)
        tz2 = (all_orbit%itzp(a) + all_orbit%itzp(b))/2
        channel2 = locate_channel(3,tz2, nx2, ny2, nz2) 
        if ( channel2 == 0 ) cycle 
           
        bra2 = pp_hhpp%ival(a,b)
        if ( bra2 == 0 ) cycle 
        
        ket4 = 0 
        do ket = 1, size(  lookup_hpphpp_configs(1,channel)%ival, 2) 
           i = lookup_hpphpp_configs(1,channel)%ival(1,ket) 
           c = lookup_hpphpp_configs(1,channel)%ival(2,ket) 
           d = lookup_hpphpp_configs(1,channel)%ival(3,ket) 
           
           if ( c >= d ) cycle 
           ket4 = ket4 + 1
           
           nx3 = all_orbit%nx(c) + all_orbit%nx(d)
           ny3 = all_orbit%ny(c) + all_orbit%ny(d)
           nz3 = all_orbit%nz(c) + all_orbit%nz(d)
           tz3 = (all_orbit%itzp(c) + all_orbit%itzp(d))/2
           channel3 = locate_channel(3,tz3, nx3, ny3, nz3) 
           if ( channel3 == 0 ) cycle 
           
           bra3 = pp_hhpp%ival(c,d)
           if ( bra3 == 0 ) cycle 
           
           nxd = all_orbit%nx(a) + all_orbit%nx(b) - all_orbit%nx(i)
           nyd = all_orbit%ny(a) + all_orbit%ny(b) - all_orbit%ny(i) 
           nzd = all_orbit%nz(a) + all_orbit%nz(b) - all_orbit%nz(i) 
           tzd = all_orbit%itzp(a) + all_orbit%itzp(b) - all_orbit%itzp(i)
           
           if ( abs(nxd) > nkxmax ) cycle
           if ( abs(nyd) > nkymax ) cycle
           if ( abs(nzd) > nkzmax ) cycle
           if ( abs(tzd) > 1 ) cycle
           
           do szd = -1, 1 ,2 
              
              j = orbitof(nxd, nyd, nzd, szd, tzd)
              if ( j > below_ef ) cycle 
                            
              ket2 = hh_hhpp%ival(i,j)
              if ( ket2 == 0) cycle 
              
              ket3 = hh_hhpp%ival(k,j)
              if ( ket3 == 0) cycle 
              
              
              sum1 = -t2_ccm(channel3)%val(bra3,ket3)*v3nf_hpphpp(channel)%val(bra,ket4)
              t2_ccm_eqn(channel2)%val(bra2,ket2)  = t2_ccm_eqn(channel2)%val(bra2,ket2) + sum1 
              ket2 = hh_hhpp%ival(j,i)
              t2_ccm_eqn(channel2)%val(bra2,ket2)  = t2_ccm_eqn(channel2)%val(bra2,ket2) - sum1 
           end do
        end do
                
     end do
  end do
  
  
  call mpi_barrier(mpi_comm_world,ierror)
  endwtime = MPI_WTIME()
  if ( iam == 0 ) write(6,*) 'Time for diag1', endwtime - startwtime
  startwtime = MPI_WTIME()
  
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
     
     do bra = bra_min, bra_max 
        k = lookup_hhphhp_configs(1,channel)%ival(1,bra)
        l = lookup_hhphhp_configs(1,channel)%ival(2,bra) 
        a = lookup_hhphhp_configs(1,channel)%ival(3,bra) 
        
        if ( l > k ) cycle 
        
        nx2 = all_orbit%nx(k) + all_orbit%nx(l)
        ny2 = all_orbit%ny(k) + all_orbit%ny(l)
        nz2 = all_orbit%nz(k) + all_orbit%nz(l)
        tz2 = (all_orbit%itzp(k) + all_orbit%itzp(l))/2
        channel3 = locate_channel(3,tz2, nx2, ny2, nz2) 
        if ( channel3 == 0 ) cycle 
           
        ket3 = hh_hhpp%ival(k,l)
        if ( ket3 == 0 ) cycle 
           
        do ket = 1, size(  lookup_hhphhp_configs(1,channel)%ival, 2) 
           i = lookup_hhphhp_configs(1,channel)%ival(1,ket) 
           j = lookup_hhphhp_configs(1,channel)%ival(2,ket) 
           c = lookup_hhphhp_configs(1,channel)%ival(3,ket) 
           
           
           nx3 = all_orbit%nx(i) + all_orbit%nx(j)
           ny3 = all_orbit%ny(i) + all_orbit%ny(j)
           nz3 = all_orbit%nz(i) + all_orbit%nz(j)
           tz3 = (all_orbit%itzp(i) + all_orbit%itzp(j))/2
           channel2 = locate_channel(3,tz3, nx3, ny3, nz3) 
           if ( channel2 == 0 ) cycle 
           
           ket2 = hh_hhpp%ival(i,j)
           if ( ket2 == 0 ) cycle 
           
           nxd = all_orbit%nx(i) + all_orbit%nx(j) - all_orbit%nx(a)
           nyd = all_orbit%ny(i) + all_orbit%ny(j) - all_orbit%ny(a) 
           nzd = all_orbit%nz(i) + all_orbit%nz(j) - all_orbit%nz(a) 
           tzd = all_orbit%itzp(i) + all_orbit%itzp(j) - all_orbit%itzp(a)
           
           if ( abs(nxd) > nkxmax ) cycle
           if ( abs(nyd) > nkymax ) cycle
           if ( abs(nzd) > nkzmax ) cycle
           if ( abs(tzd) > 1 ) cycle
           
           do szd = -1, 1 ,2 
              
              b = orbitof(nxd, nyd, nzd, szd, tzd)
              if ( b <= below_ef ) cycle 
              
              bra3 = pp_hhpp%ival(c,b)
              if ( bra3 == 0) cycle 
              
              bra2 = pp_hhpp%ival(a,b)
              if ( bra2 == 0) cycle 
              
              
              sum1 = t2_ccm(channel3)%val(bra3,ket3)*v3nf_hhphhp(channel)%val(bra,ket)
              t2_ccm_eqn(channel2)%val(bra2,ket2)  = t2_ccm_eqn(channel2)%val(bra2,ket2) + sum1 
              bra2 = pp_hhpp%ival(b,a)
              t2_ccm_eqn(channel2)%val(bra2,ket2)  = t2_ccm_eqn(channel2)%val(bra2,ket2) - sum1 
           end do
        end do
     end do
  end do
  
  
  call mpi_barrier(mpi_comm_world,ierror)
  endwtime = MPI_WTIME()
  if ( iam == 0 ) write(6,*) 'Time for diag2', endwtime - startwtime
  startwtime = MPI_WTIME()

   
  allocate( x_ae(below_ef+1:tot_orbs, below_ef+1:tot_orbs)) 
  x_ae = 0.d0 
  
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
     

     do ket = bra_min, bra_max 
        c = lookup_hhpppp_configs(2,channel)%ival(1,ket)
        d = lookup_hhpppp_configs(2,channel)%ival(2,ket) 
        e = lookup_hhpppp_configs(2,channel)%ival(3,ket) 
        
        if ( d > c ) cycle 
        nx2 = all_orbit%nx(c) + all_orbit%nx(d)
        ny2 = all_orbit%ny(c) + all_orbit%ny(d)
        nz2 = all_orbit%nz(c) + all_orbit%nz(d)
        tz2 = (all_orbit%itzp(c) + all_orbit%itzp(d))/2
        channel2 = locate_channel(3,tz2, nx2, ny2, nz2) 
        if ( channel2 == 0 ) cycle 
        
        
        bra2 = pp_hhpp%ival(c,d)
        if ( bra2 == 0 ) cycle 
        

        do bra = 1, size(  lookup_hhpppp_configs(1,channel)%ival, 2) 
           k = lookup_hhpppp_configs(1,channel)%ival(1,bra) 
           l = lookup_hhpppp_configs(1,channel)%ival(2,bra) 
           a = lookup_hhpppp_configs(1,channel)%ival(3,bra) 
           
           
           if ( l > k ) cycle 

           nx2 = all_orbit%nx(k) + all_orbit%nx(l)
           ny2 = all_orbit%ny(k) + all_orbit%ny(l)
           nz2 = all_orbit%nz(k) + all_orbit%nz(l)
           tz2 = (all_orbit%itzp(k) + all_orbit%itzp(l))/2
           if ( locate_channel(3,tz2, nx2, ny2, nz2) /= channel2 ) cycle  

           ket2 = hh_hhpp%ival(k,l)
           if ( ket2 == 0 ) cycle 
           
           x_ae( a,e) = x_ae( a,e) + t2_ccm(channel2)%val(bra2,ket2)*v3nf_hhpppp(channel)%val(bra,ket)
           !x_ae( a,e) = x_ae( a,e) + v3nf_hhpppp(channel)%val(bra,ket)
           
        end do
     end do
  end do


  
  call mpi_allreduce(mpi_in_place,x_ae,size(x_ae),mpi_complex16,mpi_sum, &
       mpi_comm_world,ierror)
  
  do channel = my_hhpp_channel_low(iam), my_hhpp_channel_high(iam)
     nx = channels%hhpp_quantum_numbers(channel*4)
     ny = channels%hhpp_quantum_numbers(channel*4-1)
     nz = channels%hhpp_quantum_numbers(channel*4-2)
     tz = channels%hhpp_quantum_numbers(channel*4-3)
     
     !if ( check_my_channel_hhpp(channel) == 0 ) cycle
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
     
     do bra = bra_min, bra_max !1, size(  lookup_hhpp_configs(2,channel)%ival, 2) 
        a = lookup_hhpp_configs(2,channel)%ival(1,bra)
        b = lookup_hhpp_configs(2,channel)%ival(2,bra) 

        do ket = 1, size(  lookup_hhpp_configs(1,channel)%ival, 2) 
           i= lookup_hhpp_configs(1,channel)%ival(1,ket) 
           j = lookup_hhpp_configs(1,channel)%ival(2,ket) 
           
           sum1 = 0.d0
           do szc = -1, 1 ,2 
              nxc = all_orbit%nx(a)
              nyc = all_orbit%ny(a)
              nzc = all_orbit%nz(a)
              tzc = all_orbit%itzp(a)
              c = orbitof(nxc, nyc, nzc, szc, tzc)
              if  ( c <= below_ef ) cycle 
              
              bra2 = pp_hhpp%ival(c,b)
              if ( bra2 == 0 ) cycle 
              sum1 = sum1 + x_ae(a,c)*t2_ccm(channel)%val(bra2,ket)
           end do
           
           t2_ccm_eqn(channel)%val(bra,ket) = t2_ccm_eqn(channel)%val(bra,ket) + sum1 
           bra3 = pp_hhpp%ival(b,a) 
           t2_ccm_eqn(channel)%val(bra3,ket) = t2_ccm_eqn(channel)%val(bra3,ket) - sum1 
           
        end do
     end do
  end do
  
  deallocate( x_ae )

  call mpi_barrier(mpi_comm_world,ierror)
  endwtime = MPI_WTIME()
  if ( iam == 0 ) write(6,*) 'Time for diag3', endwtime - startwtime
  startwtime = MPI_WTIME()

  allocate( x_mi(below_ef, below_ef)) 
  x_mi = 0.d0 
  
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
     
     do ket = bra_min, bra_max 
        i = lookup_hhhhpp_configs(2,channel)%ival(1,ket)
        c = lookup_hhhhpp_configs(2,channel)%ival(2,ket) 
        d = lookup_hhhhpp_configs(2,channel)%ival(3,ket) 
        
        if ( c > d ) cycle 
        nx2 = all_orbit%nx(c) + all_orbit%nx(d)
        ny2 = all_orbit%ny(c) + all_orbit%ny(d)
        nz2 = all_orbit%nz(c) + all_orbit%nz(d)
        tz2 = (all_orbit%itzp(c) + all_orbit%itzp(d))/2
        channel2 = locate_channel(3,tz2, nx2, ny2, nz2) 
        if ( channel2 == 0 ) cycle 
        
        
        bra2 = pp_hhpp%ival(c,d)
        if ( bra2 == 0 ) cycle 
        
        
        do bra = 1, size(  lookup_hhhhpp_configs(1,channel)%ival, 2) 
           k = lookup_hhhhpp_configs(1,channel)%ival(1,bra) 
           l = lookup_hhhhpp_configs(1,channel)%ival(2,bra) 
           m = lookup_hhhhpp_configs(1,channel)%ival(3,bra) 
           
           if ( k > l ) cycle 
     
           nx2 = all_orbit%nx(k) + all_orbit%nx(l)
           ny2 = all_orbit%ny(k) + all_orbit%ny(l)
           nz2 = all_orbit%nz(k) + all_orbit%nz(l)
           tz2 = (all_orbit%itzp(k) + all_orbit%itzp(l))/2
           if ( locate_channel(3,tz2, nx2, ny2, nz2) /= channel2 ) cycle  
           
           ket2 = hh_hhpp%ival(k,l)
           if ( ket2 == 0 ) cycle 
           
           x_mi(m,i) = x_mi(m,i) + t2_ccm(channel2)%val(bra2,ket2)*v3nf_hhhhpp(channel)%val(bra,ket) 
           !x_mi(m,i) = x_mi(m,i) + v3nf_hhhhpp(channel)%val(bra,ket) 
           
        end do
     end do
  end do
  
  call mpi_allreduce(mpi_in_place,x_mi,size(x_mi),mpi_complex16,mpi_sum, &
       mpi_comm_world,ierror)
  
  do channel = my_hhpp_channel_low(iam), my_hhpp_channel_high(iam)
     nx = channels%hhpp_quantum_numbers(channel*4)
     ny = channels%hhpp_quantum_numbers(channel*4-1)
     nz = channels%hhpp_quantum_numbers(channel*4-2)
     tz = channels%hhpp_quantum_numbers(channel*4-3)
     
     !if ( check_my_channel_hhpp(channel) == 0 ) cycle
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
     
     do bra = bra_min, bra_max !1, size(  lookup_hhpp_configs(2,channel)%ival, 2) 
        a = lookup_hhpp_configs(2,channel)%ival(1,bra)
        b = lookup_hhpp_configs(2,channel)%ival(2,bra) 
        
        do ket = 1, size(  lookup_hhpp_configs(1,channel)%ival, 2) 
           i= lookup_hhpp_configs(1,channel)%ival(1,ket) 
           j = lookup_hhpp_configs(1,channel)%ival(2,ket) 
        
           sum1 = 0.d0
           do szk = -1, 1, 2
              nxk = all_orbit%nx(i)
              nyk = all_orbit%ny(i)
              nzk = all_orbit%nz(i)
              tzk = all_orbit%itzp(i)
              m = orbitof(nxk, nyk, nzk, szk, tzk)
              if  ( m > below_ef ) cycle 
              
              
              ket2 = hh_hhpp%ival(m,j)
              if ( ket2 == 0 ) cycle 
              sum1 = sum1 - x_mi(m,i)*t2_ccm(channel)%val(bra,ket2)
           end do
           
           t2_ccm_eqn(channel)%val(bra,ket) = t2_ccm_eqn(channel)%val(bra,ket) + sum1 
           ket3 = hh_hhpp%ival(j,i)
           t2_ccm_eqn(channel)%val(bra,ket3) = t2_ccm_eqn(channel)%val(bra,ket3) - sum1 
           
           
        end do
     end do
  end do

  deallocate( x_mi)
  call mpi_barrier(mpi_comm_world,ierror)
  endwtime = MPI_WTIME()
  if ( iam == 0 ) write(6,*) 'Time for diag4', endwtime - startwtime
  startwtime = MPI_WTIME()


  
  allocate( interm_hhhh(channels%number_hhhh_confs) )
  do channel   = 1, channels%number_hhhh_confs 
     nx = channels%hhhh_quantum_numbers(channel*4)
     ny = channels%hhhh_quantum_numbers(channel*4-1)
     nz = channels%hhhh_quantum_numbers(channel*4-2)
     tz = channels%hhhh_quantum_numbers(channel*4-3)

     dim1 = size(  lookup_hhhh_configs(1,channel)%ival, 2) 
     allocate( interm_hhhh(channel)%val(dim1, dim1) )
     interm_hhhh(channel)%val = 0.d0 
  end do
  
  
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
     
     do ket = bra_min, bra_max 
        j = lookup_hhhhpp_configs(2,channel)%ival(1,ket)
        c = lookup_hhhhpp_configs(2,channel)%ival(2,ket) 
        d = lookup_hhhhpp_configs(2,channel)%ival(3,ket) 
        
        if ( c > d ) cycle 
        nx2 = all_orbit%nx(c) + all_orbit%nx(d)
        ny2 = all_orbit%ny(c) + all_orbit%ny(d)
        nz2 = all_orbit%nz(c) + all_orbit%nz(d)
        tz2 = (all_orbit%itzp(c) + all_orbit%itzp(d))/2
        channel2 = locate_channel(3,tz2, nx2, ny2, nz2) 
        if ( channel2 == 0 ) cycle 
        
        bra2 = pp_hhpp%ival(c,d)
        if ( bra2 == 0 ) cycle 
        
        
        do bra = 1, size(  lookup_hhhhpp_configs(1,channel)%ival, 2) 
           m = lookup_hhhhpp_configs(1,channel)%ival(1,bra) 
           k = lookup_hhhhpp_configs(1,channel)%ival(2,bra) 
           l = lookup_hhhhpp_configs(1,channel)%ival(3,bra) 
           
           
           nx3 = all_orbit%nx(k) + all_orbit%nx(l)
           ny3 = all_orbit%ny(k) + all_orbit%ny(l)
           nz3 = all_orbit%nz(k) + all_orbit%nz(l)
           tz3 = (all_orbit%itzp(k) + all_orbit%itzp(l))/2
           channel3 = locate_channel(1,tz3, nx3, ny3, nz3) 
           if ( channel3 == 0 ) cycle 
                 
           bra3 = hh_hhhh%ival(k,l)
           if ( bra3 == 0) cycle 
           
           nxd = all_orbit%nx(k) + all_orbit%nx(l) - all_orbit%nx(j)
           nyd = all_orbit%ny(k) + all_orbit%ny(l) - all_orbit%ny(j) 
           nzd = all_orbit%nz(k) + all_orbit%nz(l) - all_orbit%nz(j) 
           tzd = all_orbit%itzp(k) + all_orbit%itzp(l) - all_orbit%itzp(j)
                 
           if ( abs(nxd) > nkxmax ) cycle
           if ( abs(nyd) > nkymax ) cycle
           if ( abs(nzd) > nkzmax ) cycle
           if ( abs(tzd) > 1 ) cycle
           
           do szd = -1, 1 ,2 
              
              i = orbitof(nxd, nyd, nzd, szd, tzd)
              if ( i > below_ef ) cycle 
              
              ket3 = hh_hhhh%ival(i,j)
              if ( ket3 == 0) cycle 
              ket2 = hh_hhpp%ival(m,i)
              if ( ket2 == 0 ) cycle 

              
              sum1 =  t2_ccm(channel2)%val(bra2,ket2)*v3nf_hhhhpp(channel)%val(bra,ket) 

              interm_hhhh(channel3)%val(bra3,ket3) = interm_hhhh(channel3)%val(bra3,ket3) + sum1 
              ket3 = hh_hhhh%ival(j,i)
              interm_hhhh(channel3)%val(bra3,ket3) = interm_hhhh(channel3)%val(bra3,ket3) - sum1 
              
           end do
        end do
     end do
  end do
  
  
  do i = 1, size(interm_hhhh)
     call mpi_allreduce(mpi_in_place,interm_hhhh(i)%val,size(interm_hhhh(i)%val),mpi_complex16,mpi_sum, &
          mpi_comm_world,ierror)

  ENDDO
  
  do channel = my_hhpp_channel_low(iam), my_hhpp_channel_high(iam)
     !do channel   = 1, channels%number_hhpp_confs 
     nx = channels%hhpp_quantum_numbers(channel*4)
     ny = channels%hhpp_quantum_numbers(channel*4-1)
     nz = channels%hhpp_quantum_numbers(channel*4-2)
     tz = channels%hhpp_quantum_numbers(channel*4-3)
     
     channel2 = locate_channel(1,tz, nx, ny, nz) 
     if ( channel2 == 0 ) cycle 
     
     !if ( check_my_channel_hhpp(channel) == 0 ) cycle
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
     
     dim1 = bra_max-bra_min + 1 
     dim2 = size( I9(channel2)%val,2 )
     dim3 = size( t2_ccm(channel)%val,2) 
     
     ! Nov. 10
     call ZGEMM ( 'n', 'n', dim1, dim2, dim3, dcmplx(0.5d0,0.d0), t2_ccm(channel)%val(bra_min:bra_max,:), & 
          dim1,interm_hhhh(channel2)%val, dim3, & 
          dcmplx(1.d0,0.d0), t2_ccm_eqn(channel)%val(bra_min:bra_max,:), dim1 )
     
  end do
  
  call mpi_barrier(mpi_comm_world,ierror)
  endwtime = MPI_WTIME()
  if ( iam == 0 ) write(6,*) 'Time for diag6', endwtime - startwtime
  startwtime = MPI_WTIME()
  
  do i = 1, size(interm_hhhh )
     deallocate( interm_hhhh(i)%val ) 
  end do
  deallocate( interm_hhhh )
  
  
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
     interm_pppp(channel2)%val = 0.d0 
  end do

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
     
     !write(6,*) bra_max-bra_min+1, size(v3nf_hhpppp(channel)%val, 2 )
     do ket = bra_min, bra_max 
        e = lookup_hhpppp_configs(2,channel)%ival(1,ket)
        c = lookup_hhpppp_configs(2,channel)%ival(2,ket) 
        d = lookup_hhpppp_configs(2,channel)%ival(3,ket) 
        
        if ( c >= d ) cycle 
        nx2 = all_orbit%nx(c) + all_orbit%nx(d)
        ny2 = all_orbit%ny(c) + all_orbit%ny(d)
        nz2 = all_orbit%nz(c) + all_orbit%nz(d)
        tz2 = (all_orbit%itzp(c) + all_orbit%itzp(d))/2
        channel2 = locate_channel(3,tz2, nx2, ny2, nz2) 
        if ( channel2 == 0 ) cycle 
        
        ket2 = pp_pppp_rest%ival(c,d)
        if ( ket2 == 0 ) cycle 
        
        do bra = 1, size(  lookup_hhpppp_configs(1,channel)%ival, 2) 
           k = lookup_hhpppp_configs(1,channel)%ival(1,bra) 
           l = lookup_hhpppp_configs(1,channel)%ival(2,bra) 
           b = lookup_hhpppp_configs(1,channel)%ival(3,bra) 
  
           if ( k >= l ) cycle 
           nx3 = all_orbit%nx(k) + all_orbit%nx(l)
           ny3 = all_orbit%ny(k) + all_orbit%ny(l)
           nz3 = all_orbit%nz(k) + all_orbit%nz(l)
           tz3 = (all_orbit%itzp(k) + all_orbit%itzp(l))/2
           channel3 = locate_channel(3,tz3, nx3, ny3, nz3) 
           if ( channel3 == 0 ) cycle 
           ket3 = hh_hhpp%ival(k,l)
           if ( ket3 == 0 ) cycle 
           
           nxd = all_orbit%nx(k) + all_orbit%nx(l) - all_orbit%nx(e)
           nyd = all_orbit%ny(k) + all_orbit%ny(l) - all_orbit%ny(e) 
           nzd = all_orbit%nz(k) + all_orbit%nz(l) - all_orbit%nz(e) 
           tzd = all_orbit%itzp(k) + all_orbit%itzp(l) - all_orbit%itzp(e)
                 
           if ( abs(nxd) > nkxmax ) cycle
           if ( abs(nyd) > nkymax ) cycle
           if ( abs(nzd) > nkzmax ) cycle
           if ( abs(tzd) > 1 ) cycle
           
           do szd = -1, 1 ,2 
              
              a = orbitof(nxd, nyd, nzd, szd, tzd)
              if ( a <= below_ef ) cycle 
              if ( a >= b ) cycle 
              
              bra3 = pp_hhpp%ival(e,a) 
              if ( bra3 == 0 ) cycle 
              bra2 = pp_pppp_rest%ival(a,b)
              if ( bra2 == 0 ) cycle 
              
              sum1 = - t2_ccm(channel3)%val(bra3,ket3) * v3nf_hhpppp(channel)%val(bra,ket)
              interm_pppp(channel2)%val(bra2,ket2) = interm_pppp(channel2)%val(bra2,ket2) + sum1 
                            
           end do
        end do
     end do
  end do

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
     
     !write(6,*) bra_max-bra_min+1, size(v3nf_hhpppp(channel)%val, 2 )
     do ket = bra_min, bra_max 
        e = lookup_hhpppp_configs(2,channel)%ival(1,ket)
        c = lookup_hhpppp_configs(2,channel)%ival(2,ket) 
        d = lookup_hhpppp_configs(2,channel)%ival(3,ket) 
        
        if ( c >= d ) cycle 
        nx2 = all_orbit%nx(c) + all_orbit%nx(d)
        ny2 = all_orbit%ny(c) + all_orbit%ny(d)
        nz2 = all_orbit%nz(c) + all_orbit%nz(d)
        tz2 = (all_orbit%itzp(c) + all_orbit%itzp(d))/2
        channel2 = locate_channel(3,tz2, nx2, ny2, nz2) 
        if ( channel2 == 0 ) cycle 
        
        ket2 = pp_pppp_rest%ival(c,d)
        if ( ket2 == 0 ) cycle 
        
        do bra = 1, size(  lookup_hhpppp_configs(1,channel)%ival, 2) 
           k = lookup_hhpppp_configs(1,channel)%ival(1,bra) 
           l = lookup_hhpppp_configs(1,channel)%ival(2,bra) 
           a = lookup_hhpppp_configs(1,channel)%ival(3,bra) 
  
           if ( k >= l ) cycle 
           nx3 = all_orbit%nx(k) + all_orbit%nx(l)
           ny3 = all_orbit%ny(k) + all_orbit%ny(l)
           nz3 = all_orbit%nz(k) + all_orbit%nz(l)
           tz3 = (all_orbit%itzp(k) + all_orbit%itzp(l))/2
           channel3 = locate_channel(3,tz3, nx3, ny3, nz3) 
           if ( channel3 == 0 ) cycle 
           ket3 = hh_hhpp%ival(k,l)
           if ( ket3 == 0 ) cycle 
           
           nxd = all_orbit%nx(k) + all_orbit%nx(l) - all_orbit%nx(e)
           nyd = all_orbit%ny(k) + all_orbit%ny(l) - all_orbit%ny(e) 
           nzd = all_orbit%nz(k) + all_orbit%nz(l) - all_orbit%nz(e) 
           tzd = all_orbit%itzp(k) + all_orbit%itzp(l) - all_orbit%itzp(e)
                 
           if ( abs(nxd) > nkxmax ) cycle
           if ( abs(nyd) > nkymax ) cycle
           if ( abs(nzd) > nkzmax ) cycle
           if ( abs(tzd) > 1 ) cycle
           
           do szd = -1, 1 ,2 
              
              b = orbitof(nxd, nyd, nzd, szd, tzd)
              if ( b <= below_ef ) cycle 
              if ( a >= b ) cycle 
              
              bra3 = pp_hhpp%ival(e,b) 
              if ( bra3 == 0 ) cycle 
              bra2 = pp_pppp_rest%ival(a,b)
              if ( bra2 == 0 ) cycle 
              
              sum1 =  t2_ccm(channel3)%val(bra3,ket3) * v3nf_hhpppp(channel)%val(bra,ket)
              interm_pppp(channel2)%val(bra2,ket2) = interm_pppp(channel2)%val(bra2,ket2) + sum1 
              
           end do
        end do
     end do
  end do

  do i = 1, size( interm_pppp ) 
     call mpi_allreduce(mpi_in_place, interm_pppp(i)%val, size(interm_pppp(i)%val), mpi_complex16, mpi_sum, &
          mpi_comm_world,ierror)
  end do
  
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
!!$     if ( check_my_channel(channel2) == 0 ) cycle
!!$     local_ch_id = channel2
!!$     
!!$     !
!!$     ! ket side if fully stored on each proc
!!$     !
!!$     ket_min = mapping(iam+1,local_ch_id,2)
!!$     ket_max = mapping(iam+1,local_ch_id,3)
!!$     !
!!$     ! bra side is distributed 
!!$     !
!!$     bra_min = mapping(iam+1,local_ch_id,4)
!!$     bra_max = mapping(iam+1,local_ch_id,5)
!!$     
!!$     do bra = bra_min, bra_max !1, size(  lookup_pppp_configs(1,channel)%ival, 2) 
!!$        a = lookup_pppp_configs(1,channel)%ival(1,bra) 
!!$        b = lookup_pppp_configs(1,channel)%ival(2,bra) 
!!$        
!!$        do ket = ket_min, ket_max !1, size(  lookup_pppp_configs(1,channel)%ival, 2) 
!!$           c = lookup_pppp_configs(1,channel)%ival(1,ket)
!!$           d = lookup_pppp_configs(1,channel)%ival(2,ket) 
!!$  
!!$           
!!$           sum1 = 0.d0 
!!$           
!!$           !$omp parallel default(shared) private(k,l,nx2,ny2,nz2,tz2,channel3,ket2,nxd,nyd,nzd,tzd,szd,e,bra2)
!!$           !$omp do schedule(dynamic), reduction(+:sum1)
!!$           do k = 1, below_ef
!!$              do l = k+1, below_ef
!!$                 
!!$                 nx2 = all_orbit%nx(k) + all_orbit%nx(l)
!!$                 ny2 = all_orbit%ny(k) + all_orbit%ny(l)
!!$                 nz2 = all_orbit%nz(k) + all_orbit%nz(l)
!!$                 tz2 = (all_orbit%itzp(k) + all_orbit%itzp(l))/2
!!$                 channel3 = locate_channel(3,tz2, nx2, ny2, nz2) 
!!$                 if ( channel3 == 0 ) cycle 
!!$                 
!!$                 ket2 = hh_hhpp%ival(k,l)
!!$                 if ( ket2 == 0) cycle 
!!$                 
!!$                 
!!$                 nxd = all_orbit%nx(k) + all_orbit%nx(l) - all_orbit%nx(a)
!!$                 nyd = all_orbit%ny(k) + all_orbit%ny(l) - all_orbit%ny(a) 
!!$                 nzd = all_orbit%nz(k) + all_orbit%nz(l) - all_orbit%nz(a) 
!!$                 tzd = all_orbit%itzp(k) + all_orbit%itzp(l) - all_orbit%itzp(a)
!!$                 
!!$                 if ( abs(nxd) > nkxmax ) cycle
!!$                 if ( abs(nyd) > nkymax ) cycle
!!$                 if ( abs(nzd) > nkzmax ) cycle
!!$                 if ( abs(tzd) > 1 ) cycle
!!$                 
!!$                 do szd = -1, 1 ,2 
!!$              
!!$                    e = orbitof(nxd, nyd, nzd, szd, tzd)
!!$                    if ( e <= below_ef ) cycle 
!!$                    
!!$                    bra2 = pp_hhpp%ival(e,a)
!!$                    if ( bra2 == 0) cycle 
!!$                    sum1 = sum1 - t2_ccm(channel3)%val(bra2,ket2) *chiral_3nf_asym(k,l,b,e,c,d)
!!$                    
!!$                 end do
!!$              end do
!!$           end do
!!$           !$omp end do
!!$           !$omp end parallel
!!$           
!!$           sum2 = 0.d0 
!!$           !$omp parallel default(shared) private(k,l,nx2,ny2,nz2,tz2,channel3,ket2,nxd,nyd,nzd,tzd,szd,e,bra2)
!!$           !$omp do schedule(dynamic), reduction(+:sum2)
!!$           do k = 1, below_ef
!!$              do l = k+1, below_ef
!!$                 
!!$                 nx2 = all_orbit%nx(k) + all_orbit%nx(l)
!!$                 ny2 = all_orbit%ny(k) + all_orbit%ny(l)
!!$                 nz2 = all_orbit%nz(k) + all_orbit%nz(l)
!!$                 tz2 = (all_orbit%itzp(k) + all_orbit%itzp(l))/2
!!$                 channel3 = locate_channel(3,tz2, nx2, ny2, nz2) 
!!$                 if ( channel3 == 0 ) cycle 
!!$                 
!!$                 ket2 = hh_hhpp%ival(k,l)
!!$                 if ( ket2 == 0) cycle 
!!$                 
!!$
!!$                 nxd = all_orbit%nx(k) + all_orbit%nx(l) - all_orbit%nx(b)
!!$                 nyd = all_orbit%ny(k) + all_orbit%ny(l) - all_orbit%ny(b) 
!!$                 nzd = all_orbit%nz(k) + all_orbit%nz(l) - all_orbit%nz(b) 
!!$                 tzd = all_orbit%itzp(k) + all_orbit%itzp(l) - all_orbit%itzp(b)
!!$                 
!!$                 if ( abs(nxd) > nkxmax ) cycle
!!$                 if ( abs(nyd) > nkymax ) cycle
!!$                 if ( abs(nzd) > nkzmax ) cycle
!!$                 if ( abs(tzd) > 1 ) cycle
!!$                 
!!$                 do szd = -1, 1 ,2 
!!$              
!!$                    e = orbitof(nxd, nyd, nzd, szd, tzd)
!!$                    if ( e <= below_ef ) cycle 
!!$                    
!!$                    bra2 = pp_hhpp%ival(e,b)
!!$                    if ( bra2 == 0) cycle 
!!$                    sum2 = sum2 - t2_ccm(channel3)%val(bra2,ket2) *chiral_3nf_asym(k,l,a,e,c,d)
!!$                    
!!$                 end do
!!$              end do
!!$           end do
!!$           !$omp end do
!!$           !$omp end parallel
!!$           
!!$           interm_pppp(channel2)%val(bra,ket) =interm_pppp(channel2)%val(bra,ket) + sum1 - sum2 
!!$
!!$        end do
!!$     end do
!!$  end do
!!$  


  

  
  do channel   = my_channel_low(iam), my_channel_high(iam) !1, channels%number_pppp_confs 
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
     
     
     
     allocate( amat( size(  lookup_pppp_configs(1,channel)%ival, 2), & 
          size(  lookup_pppp_configs(1,channel)%ival, 2) ) ) 
     amat = 0.d0 
     
     dim1 = 0 
     do bra = 1, size(  lookup_pppp_configs(1,channel)%ival, 2) 
        a = lookup_pppp_configs(1,channel)%ival(1,bra) 
        b = lookup_pppp_configs(1,channel)%ival(2,bra) 
        
        if ( a >= b ) cycle 
        dim1 = dim1 + 1 
        

        dim2 = 0 
        do ket = 1, size(  lookup_pppp_configs(1,channel)%ival, 2) 
           c = lookup_pppp_configs(1,channel)%ival(1,ket)
           d = lookup_pppp_configs(1,channel)%ival(2,ket) 

           if ( c >= d ) cycle 
           dim2 = dim2 + 1
           sum1 = interm_pppp(channel)%val(dim1, dim2) 

           amat(bra,ket) = amat(bra,ket) + sum1 
           bra2 = pp_pppp%ival(b,a) 
           amat(bra2,ket)  = amat(bra2,ket)  - sum1 
           ket2 = pp_pppp%ival(d,c) 
           amat(bra,ket2)  = amat(bra,ket2)  - sum1 
           amat(bra2,ket2) = amat(bra2,ket2) + sum1 
           
           
        end do
     end do
     
     dim1 = size( vnn_pppp(channel2)%val,1 ) 
     dim2 = size( t2_ccm(channel2)%val, 2)
     dim3 = size( vnn_pppp(channel2)%val, 2) 
     
     call ZGEMM ( 'n', 'n', dim1, dim2, dim3, DCMPLX(0.5d0,0.D0), amat(bra_min:bra_max,ket_min:ket_max), & 
          dim1,t2_ccm(channel2)%val(ket_min:ket_max,:), dim3, & 
          DCMPLX(1.d0,0.D0), t2_ccm_eqn(channel2)%val(bra_min:bra_max,:), dim1 )
     
     deallocate( amat ) 

  end do

  call mpi_barrier(mpi_comm_world,ierror)
  endwtime = MPI_WTIME()
  if ( iam == 0 ) write(6,*) 'Time for diag5', endwtime - startwtime
  startwtime = MPI_WTIME()

  
  allocate( i4_tmp(size(i4_ph)) ) 
  do i = 1, size(I4_ph ) 
     dim1 = size(i4_ph(i)%val, 1)
     dim2 = size(i4_ph(i)%val, 2)
     allocate( i4_tmp(i)%val(dim1, dim2) ) 
     i4_tmp(i)%val = 0.d0 
  end do
  
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
     
     !write(6,*) bra_max-bra_min+1, size(v3nf_hhpppp(channel)%val, 2 )
     do ket = bra_min, bra_max 
        d = lookup_hhpppp_configs(2,channel)%ival(1,ket)
        e = lookup_hhpppp_configs(2,channel)%ival(2,ket) 
        c = lookup_hhpppp_configs(2,channel)%ival(3,ket) 
        
        if ( e >= d ) cycle 
        nx2 = all_orbit%nx(e) + all_orbit%nx(d)
        ny2 = all_orbit%ny(e) + all_orbit%ny(d)
        nz2 = all_orbit%nz(e) + all_orbit%nz(d)
        tz2 = (all_orbit%itzp(e) + all_orbit%itzp(d))/2
        channel2 = locate_channel(3,tz2, nx2, ny2, nz2) 
        if ( channel2 == 0 ) cycle 
        
        bra2 = pp_hhpp%ival(d,e)
        if ( bra2 == 0 ) cycle 
        

        do bra = 1, size(  lookup_hhpppp_configs(1,channel)%ival, 2) 
           l = lookup_hhpppp_configs(1,channel)%ival(1,bra) 
           k = lookup_hhpppp_configs(1,channel)%ival(2,bra) 
           b = lookup_hhpppp_configs(1,channel)%ival(3,bra) 
           
           nx3 = all_orbit%nx(k) - all_orbit%nx(c)
           ny3 = all_orbit%ny(k) - all_orbit%ny(c)
           nz3 = all_orbit%nz(k) - all_orbit%nz(c)
           tz3 = (all_orbit%itzp(k) - all_orbit%itzp(c))/2
           channel3 = locate_ph_channel(2,tz3, nx3, ny3, nz3)
           if ( channel3 == 0 ) cycle 
           bra3 = hp_ph_hphp%ival(k,c) 
           if ( bra3 == 0 ) cycle 
           
           nxd = all_orbit%nx(d) + all_orbit%nx(e) - all_orbit%nx(l)
           nyd = all_orbit%ny(d) + all_orbit%ny(e) - all_orbit%ny(l) 
           nzd = all_orbit%nz(d) + all_orbit%nz(e) - all_orbit%nz(l) 
           tzd = all_orbit%itzp(d) + all_orbit%itzp(e) - all_orbit%itzp(l)
           
           if ( abs(nxd) > nkxmax ) cycle
           if ( abs(nyd) > nkymax ) cycle
           if ( abs(nzd) > nkzmax ) cycle
           if ( abs(tzd) > 1 ) cycle
                 
           do szd = -1, 1 ,2 
              
              i = orbitof(nxd, nyd, nzd, szd, tzd)
              if ( i > below_ef ) cycle 
              
              ket2 = hh_hhpp%ival(l,i)
              if ( ket2 == 0 ) cycle 
              ket3 = hp_ph_hphp%ival(i,b) 
              if ( ket3 == 0 ) cycle 
              
              sum1 = t2_ccm(channel2)%val(bra2,ket2)*v3nf_hhpppp(channel)%val(bra,ket)
              I4_tmp(channel3)%val(bra3,ket3) = I4_tmp(channel3)%val(bra3,ket3) + sum1
           end do
        end do
     end do
  end do
  
  
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
     
     !write(6,*) bra_max-bra_min+1, size(v3nf_hhhhpp(channel)%val, 2 )
     do ket = bra_min, bra_max 
        i = lookup_hhhhpp_configs(2,channel)%ival(1,ket)
        c = lookup_hhhhpp_configs(2,channel)%ival(2,ket) 
        d = lookup_hhhhpp_configs(2,channel)%ival(3,ket) 
        
        

        do bra = 1, size(  lookup_hhhhpp_configs(1,channel)%ival, 2) 
           l = lookup_hhhhpp_configs(1,channel)%ival(1,bra) 
           k = lookup_hhhhpp_configs(1,channel)%ival(2,bra) 
           m = lookup_hhhhpp_configs(1,channel)%ival(3,bra) 
           
           if ( l > m ) cycle 
           nx2 = all_orbit%nx(l) + all_orbit%nx(m)
           ny2 = all_orbit%ny(l) + all_orbit%ny(m)
           nz2 = all_orbit%nz(l) + all_orbit%nz(m)
           tz2 = (all_orbit%itzp(l) + all_orbit%itzp(m))/2
           channel2 = locate_channel(3,tz2, nx2, ny2, nz2) 
           if ( channel2 == 0 ) cycle 
        
        
           ket2 = hh_hhpp%ival(l,m)
           if ( ket2 == 0 ) cycle 
              
           
           nx3 = all_orbit%nx(k) - all_orbit%nx(c)
           ny3 = all_orbit%ny(k) - all_orbit%ny(c)
           nz3 = all_orbit%nz(k) - all_orbit%nz(c)
           tz3 = (all_orbit%itzp(k) - all_orbit%itzp(c))/2
           channel3 = locate_ph_channel(2,tz3, nx3, ny3, nz3)
           if ( channel3 == 0 ) cycle 
           bra3 = hp_ph_hphp%ival(k,c) 
           if ( bra3 == 0 ) cycle 
           
           nxd = all_orbit%nx(l) + all_orbit%nx(m) - all_orbit%nx(d)
           nyd = all_orbit%ny(l) + all_orbit%ny(m) - all_orbit%ny(d) 
           nzd = all_orbit%nz(l) + all_orbit%nz(m) - all_orbit%nz(d) 
           tzd = all_orbit%itzp(l) + all_orbit%itzp(m) - all_orbit%itzp(d)
           
           if ( abs(nxd) > nkxmax ) cycle
           if ( abs(nyd) > nkymax ) cycle
           if ( abs(nzd) > nkzmax ) cycle
           if ( abs(tzd) > 1 ) cycle
                 
           do szd = -1, 1 ,2 
              
              b = orbitof(nxd, nyd, nzd, szd, tzd)
              if ( b <= below_ef ) cycle 
              
              bra2 = pp_hhpp%ival(d,b)
              if ( bra2 == 0 ) cycle 

              ket3 = hp_ph_hphp%ival(i,b) 
              if ( ket3 == 0 ) cycle 
              
              sum1 = -t2_ccm(channel2)%val(bra2,ket2)*v3nf_hhhhpp(channel)%val(bra,ket)
              I4_tmp(channel3)%val(bra3,ket3) = I4_tmp(channel3)%val(bra3,ket3) + sum1
           end do
        end do
     end do
  end do
  

  do i = 1, size( i4_tmp ) 
     
     call mpi_allreduce(mpi_in_place,i4_tmp(i)%val,size(i4_tmp(i)%val),mpi_complex16,mpi_sum, &
          mpi_comm_world,ierror)
  end do
  
  do channel   = my_ph_hphp_channel_low(iam), my_ph_hphp_channel_high(iam) 
     !do channel   = 1, channels%number_ph_hphp_confs 
     nx = channels%ph_hphp_quantum_numbers(channel*4)
     ny = channels%ph_hphp_quantum_numbers(channel*4-1)
     nz = channels%ph_hphp_quantum_numbers(channel*4-2)
     tz = channels%ph_hphp_quantum_numbers(channel*4-3)
     
     !if ( check_my_channel_ph_hphp(channel) == 0 ) cycle
     local_ch_id = channel
     !
     ! ket side if fully stored on each proc
     !
     bra_min = mapping_ph_hphp(iam+1,local_ch_id,2)
     bra_max = mapping_ph_hphp(iam+1,local_ch_id,3)
     !
     ! bra side is distributed 
     !
     ket_min = mapping_ph_hphp(iam+1,local_ch_id,4)
     ket_max = mapping_ph_hphp(iam+1,local_ch_id,5)
     
     channel2 = locate_ph_channel(1,tz, nx, ny, nz)
     if ( channel2 == 0 ) cycle 
     
     dim1 = size( t2_ph_ccm(channel2)%val,1 )
     dim2 = bra_max-bra_min+1 
     dim3 = size( t2_ph_ccm(channel2)%val, 2) 
     
     call ZGEMM ( 'n', 'n', dim1, dim2, dim3, dcmplx(1.d0,0.d0), t2_ph_ccm(channel2)%val, & 
          dim1,I4_tmp(channel)%val(:,bra_min:bra_max), dim3, & 
          dcmplx(1.d0,0.d0), t2_ph_ccm_eqn(channel2)%val(:,bra_min:bra_max), dim1 )
     
  end do
  
  do i = 1, size( i4_tmp )
     deallocate(i4_tmp(i)%val )
  end do
  deallocate( i4_tmp ) 
  
  call mpi_barrier(mpi_comm_world,ierror)
  endwtime = MPI_WTIME()
  if ( iam == 0 ) write(6,*) 'Time for diag7/8', endwtime - startwtime
  
  deallocate( ab_conf, ij_conf )

  
end SUBROUTINE t2diags_v3t2_linear2

SUBROUTINE t2diags_v3t2_linear
  USE single_particle_orbits
  USE configurations
  USE constants
  use t2_storage
  use one_body_operators
  use parallel 
  use chiral_potentials 
  
  IMPLICIT NONE
  INTEGER :: i,j,k,l,m,a,b,c,d,e, nx,ny,nz,sz,tz
  INTEGER :: bra,ket, channel, bra2,ket2, bra3,ket3,dim1, dim2, dim3, number_channels
  COMPLEX*16 ::  sum1, sum2, t2dum 
  INTEGER :: nxd, nyd, nzd, szd, tzd, nxc, nyc, nzc, szc, tzc, nxk, nyk, nzk, szk, tzk
  COMPLEX*16, ALLOCATABLE :: temp1(:,:), temp2(:,:), temp(:,:), t3_tmp1(:,:), t3_tmp2(:,:),t2_tmp1(:), t2_tmp2(:), amat(:,:)
  real*8  ::  startwtime , endwtime
  INTEGER :: nx2, ny2, nz2, sz2, tz2, channel2, channel3, channel1
  INTEGER :: nx3, ny3, nz3, sz3, tz3, channel4, ndim, i1
  INTEGER :: bra_min, bra_max, ket_min, ket_max, local_ch_id
  complex*16, allocatable :: x_ae(:,:), x_mi(:,:)
  
  
  call mpi_barrier(mpi_comm_world,ierror)
  startwtime = MPI_WTIME()
  
  do channel = my_hhpp_channel_low(iam), my_hhpp_channel_high(iam)
     nx = channels%hhpp_quantum_numbers(channel*4)
     ny = channels%hhpp_quantum_numbers(channel*4-1)
     nz = channels%hhpp_quantum_numbers(channel*4-2)
     tz = channels%hhpp_quantum_numbers(channel*4-3)
     
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
     
     do bra = bra_min, bra_max 
        a = lookup_hhpp_configs(2,channel)%ival(1,bra)
        b = lookup_hhpp_configs(2,channel)%ival(2,bra) 

        do ket = 1, size(  lookup_hhpp_configs(1,channel)%ival, 2) 
           i= lookup_hhpp_configs(1,channel)%ival(1,ket) 
           j = lookup_hhpp_configs(1,channel)%ival(2,ket) 
           
           
           sum1 = 0.d0 
           !$omp parallel default(shared) private(k,nx2,ny2,nz2,tz2,channel2,ket2,d,szd, & 
           !$omp nxd, nyd,nzd,tzd,c,bra2 )
           !$omp do schedule(dynamic), reduction(+:sum1)
           do k = 1, below_ef
              
              nx2 = all_orbit%nx(j) + all_orbit%nx(k)
              ny2 = all_orbit%ny(j) + all_orbit%ny(k)
              nz2 = all_orbit%nz(j) + all_orbit%nz(k)
              tz2 = (all_orbit%itzp(j) + all_orbit%itzp(k))/2
              channel2 = locate_channel(3,tz2, nx2, ny2, nz2) 
              if ( channel2 == 0 ) cycle 
              
              ket2 = hh_hhpp%ival(k,j)
              if ( ket2 == 0) cycle 
              do d = below_ef+1, tot_orbs
                 
                 nxd = all_orbit%nx(k) + all_orbit%nx(j) - all_orbit%nx(d)
                 nyd = all_orbit%ny(k) + all_orbit%ny(j) - all_orbit%ny(d) 
                 nzd = all_orbit%nz(k) + all_orbit%nz(j) - all_orbit%nz(d) 
                 tzd = all_orbit%itzp(k) + all_orbit%itzp(j) - all_orbit%itzp(d)
                 
                 if ( abs(nxd) > nkxmax ) cycle
                 if ( abs(nyd) > nkymax ) cycle
                 if ( abs(nzd) > nkzmax ) cycle
                 if ( abs(tzd) > 1 ) cycle
                 
                 do szd = -1, 1 ,2 
                    
                    c = orbitof(nxd, nyd, nzd, szd, tzd)
                    if ( c <= below_ef ) cycle 
                    if ( d > c ) cycle 
                    
                    bra2 = pp_hhpp%ival(c,d)
                    if ( bra2 == 0 ) cycle 
                    if ( abs(  t2_ccm(channel2)%val(bra2,ket2) ) < 1.e-8 ) cycle 
                    
                    sum1 = sum1 + t2_ccm(channel2)%val(bra2,ket2)*chiral_3nf_asym(a,k,b,i,c,d)
                 end do
              end do
           end do
           !$omp end do
           !$omp end parallel

           ket3 = hh_hhpp%ival(j,i)
           t2_ccm_eqn(channel)%val(bra,ket)  = t2_ccm_eqn(channel)%val(bra,ket) + sum1 
           t2_ccm_eqn(channel)%val(bra,ket3) = t2_ccm_eqn(channel)%val(bra,ket3) - sum1 
           
        end do
     end do
  end do
  
  
  call mpi_barrier(mpi_comm_world,ierror)
  endwtime = MPI_WTIME()
  if ( iam == 0 ) write(6,*) 'Time for diag1', endwtime - startwtime
  startwtime = MPI_WTIME()
  
  do channel = my_hhpp_channel_low(iam), my_hhpp_channel_high(iam)
     nx = channels%hhpp_quantum_numbers(channel*4)
     ny = channels%hhpp_quantum_numbers(channel*4-1)
     nz = channels%hhpp_quantum_numbers(channel*4-2)
     tz = channels%hhpp_quantum_numbers(channel*4-3)
     
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
     
     do bra = bra_min, bra_max 
        a = lookup_hhpp_configs(2,channel)%ival(1,bra)
        b = lookup_hhpp_configs(2,channel)%ival(2,bra) 

        do ket = 1, size(  lookup_hhpp_configs(1,channel)%ival, 2) 
           i= lookup_hhpp_configs(1,channel)%ival(1,ket) 
           j = lookup_hhpp_configs(1,channel)%ival(2,ket) 
           
           
           sum1 = 0.d0 
           !$omp parallel default(shared) private(c,nx2,ny2,nz2,tz2,channel2,bra2,k,szd,nxd,nyd,nzd,tzd,l,ket2)
           !$omp do schedule(dynamic), reduction(+:sum1)
           do k = 1, below_ef
              do l = k+1, below_ef
                 
                 nx2 = all_orbit%nx(k) + all_orbit%nx(l)
                 ny2 = all_orbit%ny(k) + all_orbit%ny(l)
                 nz2 = all_orbit%nz(k) + all_orbit%nz(l)
                 tz2 = (all_orbit%itzp(k) + all_orbit%itzp(l))/2
                 channel2 = locate_channel(3,tz2, nx2, ny2, nz2) 
                 if ( channel2 == 0 ) cycle 
                 
                 ket2 = hh_hhpp%ival(k,l)
                 if ( ket2 == 0 ) cycle 
                 
                 nxd = all_orbit%nx(k) + all_orbit%nx(l) - all_orbit%nx(b)
                 nyd = all_orbit%ny(k) + all_orbit%ny(l) - all_orbit%ny(b) 
                 nzd = all_orbit%nz(k) + all_orbit%nz(l) - all_orbit%nz(b) 
                 tzd = all_orbit%itzp(k) + all_orbit%itzp(l) - all_orbit%itzp(b)
                 
                 if ( abs(nxd) > nkxmax ) cycle
                 if ( abs(nyd) > nkymax ) cycle
                 if ( abs(nzd) > nkzmax ) cycle
                 if ( abs(tzd) > 1 ) cycle
                 
                 do szd = -1, 1 ,2 
                    
                    c = orbitof(nxd, nyd, nzd, szd, tzd)
                    if ( c <= below_ef ) cycle 
                    bra2 = pp_hhpp%ival(c,b)
                    if ( bra2 == 0) cycle 
                    if ( abs(  t2_ccm(channel2)%val(bra2,ket2) ) < 1.e-8 ) cycle 

                    sum1 = sum1 - t2_ccm(channel2)%val(bra2,ket2)*chiral_3nf_asym(a,k,l,i,c,j)
                 end do
              end do
           end do
           !$omp end do
           !$omp end parallel
           
           bra3 = pp_hhpp%ival(b,a)
           t2_ccm_eqn(channel)%val(bra,ket)  = t2_ccm_eqn(channel)%val(bra,ket) + sum1 
           t2_ccm_eqn(channel)%val(bra3,ket) = t2_ccm_eqn(channel)%val(bra3,ket) - sum1 
           
        end do
     end do
  end do

  
  call mpi_barrier(mpi_comm_world,ierror)
  endwtime = MPI_WTIME()
  if ( iam == 0 ) write(6,*) 'Time for diag2', endwtime - startwtime
  startwtime = MPI_WTIME()
  
  
  allocate( x_ae(below_ef+1:tot_orbs, below_ef+1:tot_orbs)) 
  x_ae = 0.d0 
  
  !
  ! nonlinear in t2 
  !
  do channel = my_hhpp_channel_low(iam), my_hhpp_channel_high(iam)
     nx = channels%hhpp_quantum_numbers(channel*4)
     ny = channels%hhpp_quantum_numbers(channel*4-1)
     nz = channels%hhpp_quantum_numbers(channel*4-2)
     tz = channels%hhpp_quantum_numbers(channel*4-3)
     
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
     
     do bra = bra_min, bra_max 
        c = lookup_hhpp_configs(2,channel)%ival(1,bra)
        d = lookup_hhpp_configs(2,channel)%ival(2,bra) 
        
        if ( d > c ) cycle 
        
        do ket = 1, size(  lookup_hhpp_configs(1,channel)%ival, 2) 
           k= lookup_hhpp_configs(1,channel)%ival(1,ket) 
           l = lookup_hhpp_configs(1,channel)%ival(2,ket) 
           
           if ( l > k ) cycle 
   
           !$omp parallel default(shared) private(a,szc,nxc,nyc,nzc,tzc,e)
           !$omp do schedule(dynamic), reduction(+:x_ae)
           do a = below_ef+1, tot_orbs
              do szc = -1, 1 ,2 
                 nxc = all_orbit%nx(a)
                 nyc = all_orbit%ny(a)
                 nzc = all_orbit%nz(a)
                 tzc = all_orbit%itzp(a)
                 e = orbitof(nxc, nyc, nzc, szc, tzc)
                 if  ( e <= below_ef ) cycle 
                 
                 !x_ae( a,e) = x_ae( a,e) +  0.25d0 * t2_ccm(channel)%val(bra,ket)*chiral_3nf_asym(k,l,a,c,d,e)
                 x_ae( a,e) = x_ae( a,e) + t2_ccm(channel)%val(bra,ket)*chiral_3nf_asym(k,l,a,c,d,e)
              end do
           end do
           !$omp end do
           !$omp end parallel
           
           
           
        end do
     end do
  end do
  
  call mpi_allreduce(mpi_in_place,x_ae,size(x_ae),mpi_complex16,mpi_sum, &
       mpi_comm_world,ierror)
  
  do channel = my_hhpp_channel_low(iam), my_hhpp_channel_high(iam)
     nx = channels%hhpp_quantum_numbers(channel*4)
     ny = channels%hhpp_quantum_numbers(channel*4-1)
     nz = channels%hhpp_quantum_numbers(channel*4-2)
     tz = channels%hhpp_quantum_numbers(channel*4-3)
     
     !if ( check_my_channel_hhpp(channel) == 0 ) cycle
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
     
     do bra = bra_min, bra_max !1, size(  lookup_hhpp_configs(2,channel)%ival, 2) 
        a = lookup_hhpp_configs(2,channel)%ival(1,bra)
        b = lookup_hhpp_configs(2,channel)%ival(2,bra) 

        do ket = 1, size(  lookup_hhpp_configs(1,channel)%ival, 2) 
           i= lookup_hhpp_configs(1,channel)%ival(1,ket) 
           j = lookup_hhpp_configs(1,channel)%ival(2,ket) 
           
           sum1 = 0.d0
           do szc = -1, 1 ,2 
              nxc = all_orbit%nx(a)
              nyc = all_orbit%ny(a)
              nzc = all_orbit%nz(a)
              tzc = all_orbit%itzp(a)
              c = orbitof(nxc, nyc, nzc, szc, tzc)
              if  ( c <= below_ef ) cycle 
              
              bra2 = pp_hhpp%ival(c,b)
              if ( bra2 == 0 ) cycle 
              sum1 = sum1 + x_ae(a,c)*t2_ccm(channel)%val(bra2,ket)
           end do
           
           t2_ccm_eqn(channel)%val(bra,ket) = t2_ccm_eqn(channel)%val(bra,ket) + sum1 
           bra3 = pp_hhpp%ival(b,a) 
           t2_ccm_eqn(channel)%val(bra3,ket) = t2_ccm_eqn(channel)%val(bra3,ket) - sum1 
           
        end do
     end do
  end do
  
  deallocate( x_ae )

  call mpi_barrier(mpi_comm_world,ierror)
  endwtime = MPI_WTIME()
  if ( iam == 0 ) write(6,*) 'Time for diag3', endwtime - startwtime
  startwtime = MPI_WTIME()
  

  allocate( x_mi(below_ef, below_ef)) 
  x_mi = 0.d0 
  
  do channel = my_hhpp_channel_low(iam), my_hhpp_channel_high(iam)
     nx = channels%hhpp_quantum_numbers(channel*4)
     ny = channels%hhpp_quantum_numbers(channel*4-1)
     nz = channels%hhpp_quantum_numbers(channel*4-2)
     tz = channels%hhpp_quantum_numbers(channel*4-3)
     
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
     
     do bra = bra_min, bra_max 
        c = lookup_hhpp_configs(2,channel)%ival(1,bra)
        d = lookup_hhpp_configs(2,channel)%ival(2,bra) 

        IF ( d > c ) cycle
        do ket = 1, size(  lookup_hhpp_configs(1,channel)%ival, 2) 
           k= lookup_hhpp_configs(1,channel)%ival(1,ket) 
           l = lookup_hhpp_configs(1,channel)%ival(2,ket) 
           
           if ( l > k ) cycle 
           !$omp parallel default(shared) private(m,szc,nxc,nyc,nzc,tzc,i)
           !$omp do schedule(dynamic), reduction(+:x_mi)
           do m = 1, below_ef
              do szc = -1, 1 ,2 
                 nxc = all_orbit%nx(m)
                 nyc = all_orbit%ny(m)
                 nzc = all_orbit%nz(m)
                 tzc = all_orbit%itzp(m)
                 i = orbitof(nxc, nyc, nzc, szc, tzc)
                 if  ( i > below_ef ) cycle 
                 
                 !x_mi(m,i) = x_mi(m,i) +  0.25d0 * t2_ccm(channel)%val(bra,ket)*chiral_3nf_asym(k,l,m,c,d,i)
                 x_mi(m,i) = x_mi(m,i) + t2_ccm(channel)%val(bra,ket)*chiral_3nf_asym(k,l,m,c,d,i)
              end do
           end do
           !$omp end do
           !$omp end parallel
           
           
           
        end do
     end do
  end do

  call mpi_allreduce(mpi_in_place,x_mi,size(x_mi),mpi_complex16,mpi_sum, &
       mpi_comm_world,ierror)
  
  do channel = my_hhpp_channel_low(iam), my_hhpp_channel_high(iam)
     nx = channels%hhpp_quantum_numbers(channel*4)
     ny = channels%hhpp_quantum_numbers(channel*4-1)
     nz = channels%hhpp_quantum_numbers(channel*4-2)
     tz = channels%hhpp_quantum_numbers(channel*4-3)
     
     !if ( check_my_channel_hhpp(channel) == 0 ) cycle
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
     
     do bra = bra_min, bra_max !1, size(  lookup_hhpp_configs(2,channel)%ival, 2) 
        a = lookup_hhpp_configs(2,channel)%ival(1,bra)
        b = lookup_hhpp_configs(2,channel)%ival(2,bra) 
        
        do ket = 1, size(  lookup_hhpp_configs(1,channel)%ival, 2) 
           i= lookup_hhpp_configs(1,channel)%ival(1,ket) 
           j = lookup_hhpp_configs(1,channel)%ival(2,ket) 
        
           sum1 = 0.d0
           do szk = -1, 1, 2
              nxk = all_orbit%nx(i)
              nyk = all_orbit%ny(i)
              nzk = all_orbit%nz(i)
              tzk = all_orbit%itzp(i)
              m = orbitof(nxk, nyk, nzk, szk, tzk)
              if  ( m > below_ef ) cycle 
              
              
              ket2 = hh_hhpp%ival(m,j)
              if ( ket2 == 0 ) cycle 
              sum1 = sum1 - x_mi(m,i)*t2_ccm(channel)%val(bra,ket2)
           end do
           
           t2_ccm_eqn(channel)%val(bra,ket) = t2_ccm_eqn(channel)%val(bra,ket) + sum1 
           ket3 = hh_hhpp%ival(j,i)
           t2_ccm_eqn(channel)%val(bra,ket3) = t2_ccm_eqn(channel)%val(bra,ket3) - sum1 
           
           
        end do
     end do
  end do

  deallocate( x_mi)
  call mpi_barrier(mpi_comm_world,ierror)
  endwtime = MPI_WTIME()
  if ( iam == 0 ) write(6,*) 'Time for diag4', endwtime - startwtime
  startwtime = MPI_WTIME()
  
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

     interm_pppp(channel2)%val = 0.d0 
  end do
  
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
     
     do bra = bra_min, bra_max !1, size(  lookup_pppp_configs(1,channel)%ival, 2) 
        a = lookup_pppp_configs(1,channel)%ival(1,bra) 
        b = lookup_pppp_configs(1,channel)%ival(2,bra) 
        
        do ket = ket_min, ket_max !1, size(  lookup_pppp_configs(1,channel)%ival, 2) 
           c = lookup_pppp_configs(1,channel)%ival(1,ket)
           d = lookup_pppp_configs(1,channel)%ival(2,ket) 
  
           
           sum1 = 0.d0 
           
           !$omp parallel default(shared) private(k,l,nx2,ny2,nz2,tz2,channel3,ket2,nxd,nyd,nzd,tzd,szd,e,bra2)
           !$omp do schedule(dynamic), reduction(+:sum1)
           do k = 1, below_ef
              do l = k+1, below_ef
                 
                 nx2 = all_orbit%nx(k) + all_orbit%nx(l)
                 ny2 = all_orbit%ny(k) + all_orbit%ny(l)
                 nz2 = all_orbit%nz(k) + all_orbit%nz(l)
                 tz2 = (all_orbit%itzp(k) + all_orbit%itzp(l))/2
                 channel3 = locate_channel(3,tz2, nx2, ny2, nz2) 
                 if ( channel3 == 0 ) cycle 
                 
                 ket2 = hh_hhpp%ival(k,l)
                 if ( ket2 == 0) cycle 
                 
                 
                 nxd = all_orbit%nx(k) + all_orbit%nx(l) - all_orbit%nx(a)
                 nyd = all_orbit%ny(k) + all_orbit%ny(l) - all_orbit%ny(a) 
                 nzd = all_orbit%nz(k) + all_orbit%nz(l) - all_orbit%nz(a) 
                 tzd = all_orbit%itzp(k) + all_orbit%itzp(l) - all_orbit%itzp(a)
                 
                 if ( abs(nxd) > nkxmax ) cycle
                 if ( abs(nyd) > nkymax ) cycle
                 if ( abs(nzd) > nkzmax ) cycle
                 if ( abs(tzd) > 1 ) cycle
                 
                 do szd = -1, 1 ,2 
              
                    e = orbitof(nxd, nyd, nzd, szd, tzd)
                    if ( e <= below_ef ) cycle 
                    
                    bra2 = pp_hhpp%ival(e,a)
                    if ( bra2 == 0) cycle 
                    sum1 = sum1 - t2_ccm(channel3)%val(bra2,ket2) *chiral_3nf_asym(k,l,b,e,c,d)
                    
                 end do
              end do
           end do
           !$omp end do
           !$omp end parallel
           
           sum2 = 0.d0 
           !$omp parallel default(shared) private(k,l,nx2,ny2,nz2,tz2,channel3,ket2,nxd,nyd,nzd,tzd,szd,e,bra2)
           !$omp do schedule(dynamic), reduction(+:sum2)
           do k = 1, below_ef
              do l = k+1, below_ef
                 
                 nx2 = all_orbit%nx(k) + all_orbit%nx(l)
                 ny2 = all_orbit%ny(k) + all_orbit%ny(l)
                 nz2 = all_orbit%nz(k) + all_orbit%nz(l)
                 tz2 = (all_orbit%itzp(k) + all_orbit%itzp(l))/2
                 channel3 = locate_channel(3,tz2, nx2, ny2, nz2) 
                 if ( channel3 == 0 ) cycle 
                 
                 ket2 = hh_hhpp%ival(k,l)
                 if ( ket2 == 0) cycle 
                 

                 nxd = all_orbit%nx(k) + all_orbit%nx(l) - all_orbit%nx(b)
                 nyd = all_orbit%ny(k) + all_orbit%ny(l) - all_orbit%ny(b) 
                 nzd = all_orbit%nz(k) + all_orbit%nz(l) - all_orbit%nz(b) 
                 tzd = all_orbit%itzp(k) + all_orbit%itzp(l) - all_orbit%itzp(b)
                 
                 if ( abs(nxd) > nkxmax ) cycle
                 if ( abs(nyd) > nkymax ) cycle
                 if ( abs(nzd) > nkzmax ) cycle
                 if ( abs(tzd) > 1 ) cycle
                 
                 do szd = -1, 1 ,2 
              
                    e = orbitof(nxd, nyd, nzd, szd, tzd)
                    if ( e <= below_ef ) cycle 
                    
                    bra2 = pp_hhpp%ival(e,b)
                    if ( bra2 == 0) cycle 
                    sum2 = sum2 - t2_ccm(channel3)%val(bra2,ket2) *chiral_3nf_asym(k,l,a,e,c,d)
                    
                 end do
              end do
           end do
           !$omp end do
           !$omp end parallel
           
           interm_pppp(channel2)%val(bra,ket) =interm_pppp(channel2)%val(bra,ket) + sum1 - sum2 
           
        end do
     end do
  end do
  


  

  
  do channel   = my_channel_low(iam), my_channel_high(iam) !1, channels%number_pppp_confs 
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
     
     dim1 = size( vnn_pppp(channel2)%val,1 ) 
     dim2 = size( t2_ccm(channel2)%val, 2)
     dim3 = size( vnn_pppp(channel2)%val, 2) 
     
     call ZGEMM ( 'n', 'n', dim1, dim2, dim3, DCMPLX(0.5d0,0.D0), interm_pppp(channel2)%val(bra_min:bra_max,ket_min:ket_max), & 
          dim1,t2_ccm(channel2)%val(ket_min:ket_max,:), dim3, & 
          DCMPLX(1.d0,0.D0), t2_ccm_eqn(channel2)%val(bra_min:bra_max,:), dim1 )
     
  end do

  call mpi_barrier(mpi_comm_world,ierror)
  endwtime = MPI_WTIME()
  if ( iam == 0 ) write(6,*) 'Time for diag5', endwtime - startwtime
  startwtime = MPI_WTIME()

!!$  
  
  allocate( interm_hhhh(channels%number_hhhh_confs) )
  do channel   = 1, channels%number_hhhh_confs 
     nx = channels%hhhh_quantum_numbers(channel*4)
     ny = channels%hhhh_quantum_numbers(channel*4-1)
     nz = channels%hhhh_quantum_numbers(channel*4-2)
     tz = channels%hhhh_quantum_numbers(channel*4-3)

     dim1 = size(  lookup_hhhh_configs(1,channel)%ival, 2) 
     allocate( interm_hhhh(channel)%val(dim1, dim1) )
     interm_hhhh(channel)%val = 0.d0 
  end do
  
  
  
  do channel = my_hhpp_channel_low(iam), my_hhpp_channel_high(iam)
     nx = channels%hhpp_quantum_numbers(channel*4)
     ny = channels%hhpp_quantum_numbers(channel*4-1)
     nz = channels%hhpp_quantum_numbers(channel*4-2)
     tz = channels%hhpp_quantum_numbers(channel*4-3)
     
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
     
     do bra = bra_min, bra_max 
        c = lookup_hhpp_configs(2,channel)%ival(1,bra)
        d = lookup_hhpp_configs(2,channel)%ival(2,bra) 

        if ( d > c ) cycle 
        
        do ket = 1, size(  lookup_hhpp_configs(1,channel)%ival, 2) 
           m= lookup_hhpp_configs(1,channel)%ival(1,ket) 
           i = lookup_hhpp_configs(1,channel)%ival(2,ket) 
           
           
           do k = 1, below_ef
              do l = 1, below_ef
                 
                 nx2 = all_orbit%nx(k) + all_orbit%nx(l)
                 ny2 = all_orbit%ny(k) + all_orbit%ny(l)
                 nz2 = all_orbit%nz(k) + all_orbit%nz(l)
                 tz2 = (all_orbit%itzp(k) + all_orbit%itzp(l))/2
                 channel3 = locate_channel(1,tz2, nx2, ny2, nz2) 
                 if ( channel3 == 0 ) cycle 
                 
                 bra2 = hh_hhhh%ival(k,l)
                 if ( bra2 == 0) cycle 
                 
                 
                 nxd = all_orbit%nx(k) + all_orbit%nx(l) - all_orbit%nx(i)
                 nyd = all_orbit%ny(k) + all_orbit%ny(l) - all_orbit%ny(i) 
                 nzd = all_orbit%nz(k) + all_orbit%nz(l) - all_orbit%nz(i) 
                 tzd = all_orbit%itzp(k) + all_orbit%itzp(l) - all_orbit%itzp(i)
                 
                 if ( abs(nxd) > nkxmax ) cycle
                 if ( abs(nyd) > nkymax ) cycle
                 if ( abs(nzd) > nkzmax ) cycle
                 if ( abs(tzd) > 1 ) cycle
                 
                 do szd = -1, 1 ,2 
              
                    j = orbitof(nxd, nyd, nzd, szd, tzd)
                    if ( j > below_ef ) cycle 
                    
                    ket2 = hh_hhhh%ival(i,j)
                    if ( ket2 == 0) cycle 
                    
                    sum1 = t2_ccm(channel)%val(bra,ket)*chiral_3nf_asym(m,k,l,c,d,j)
                    
                    interm_hhhh(channel3)%val(bra2,ket2) = interm_hhhh(channel3)%val(bra2,ket2) + sum1 
                    ket2 = hh_hhhh%ival(j,i)
                    interm_hhhh(channel3)%val(bra2,ket2) = interm_hhhh(channel3)%val(bra2,ket2) - sum1 
                                        
                 end do
              end do
           end do
        end do
     end do
  end do
  
  do i = 1, size(interm_hhhh)
     call mpi_allreduce(mpi_in_place,interm_hhhh(i)%val,size(interm_hhhh(i)%val),mpi_complex16,mpi_sum, &
          mpi_comm_world,ierror)

  ENDDO

  do channel = my_hhpp_channel_low(iam), my_hhpp_channel_high(iam)
     !do channel   = 1, channels%number_hhpp_confs 
     nx = channels%hhpp_quantum_numbers(channel*4)
     ny = channels%hhpp_quantum_numbers(channel*4-1)
     nz = channels%hhpp_quantum_numbers(channel*4-2)
     tz = channels%hhpp_quantum_numbers(channel*4-3)
     
     channel2 = locate_channel(1,tz, nx, ny, nz) 
     if ( channel2 == 0 ) cycle 
     
     !if ( check_my_channel_hhpp(channel) == 0 ) cycle
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
     
     dim1 = bra_max-bra_min + 1 
     dim2 = size( I9(channel2)%val,2 )
     dim3 = size( t2_ccm(channel)%val,2) 
     
     ! Nov. 10
     call ZGEMM ( 'n', 'n', dim1, dim2, dim3, dcmplx(0.5d0,0.d0), t2_ccm(channel)%val(bra_min:bra_max,:), & 
          dim1,interm_hhhh(channel2)%val, dim3, & 
          dcmplx(1.d0,0.d0), t2_ccm_eqn(channel)%val(bra_min:bra_max,:), dim1 )
     
  end do
  
  call mpi_barrier(mpi_comm_world,ierror)
  endwtime = MPI_WTIME()
  if ( iam == 0 ) write(6,*) 'Time for diag6', endwtime - startwtime
  startwtime = MPI_WTIME()
  
  do i = 1, size(interm_hhhh )
     deallocate( interm_hhhh(i)%val ) 
  end do
  deallocate( interm_hhhh )

  
  



  
  do channel   = my_ph_hphp_channel_low(iam), my_ph_hphp_channel_high(iam) 
     
     local_ch_id = channel
     !
     ! ket side if fully stored on each proc
     !
     bra_min = mapping_ph_hphp(iam+1,local_ch_id,2)
     bra_max = mapping_ph_hphp(iam+1,local_ch_id,3)
     !
     ! bra side is distributed 
     !
     ket_min = mapping_ph_hphp(iam+1,local_ch_id,4)
     ket_max = mapping_ph_hphp(iam+1,local_ch_id,5)
     
     do ket = bra_min, bra_max !1, size(  lookup_ph_hphp_configs(1,channel)%ival, 2) 
        i = lookup_ph_hphp_configs(1,channel)%ival(1,ket)
        b = lookup_ph_hphp_configs(1,channel)%ival(2,ket) 
     
        do bra = 1, size(  lookup_ph_hphp_configs(1,channel)%ival, 2) 
           k = lookup_ph_hphp_configs(1,channel)%ival(1,bra) 
           c = lookup_ph_hphp_configs(1,channel)%ival(2,bra) 
           
           sum1= 0.d0 
           !$omp parallel default(shared) private(l,d,nx2,ny2,nz2,tz2,channel2,ket2,nxd,nyd,nzd,tzd,szd,e,bra2)
           !$omp do schedule(dynamic), reduction(+:sum1)
           do l = 1, below_ef
              nx2 = all_orbit%nx(i) + all_orbit%nx(l)
              ny2 = all_orbit%ny(i) + all_orbit%ny(l)
              nz2 = all_orbit%nz(i) + all_orbit%nz(l)
              tz2 = (all_orbit%itzp(i) + all_orbit%itzp(l))/2
              channel2 = locate_channel(3,tz2, nx2, ny2, nz2) 
              if ( channel2 == 0 ) cycle 
              ket2 = hh_hhpp%ival(l,i)
              if ( ket2 == 0 ) cycle
              
              do d = below_ef+1, tot_orbs
                 
                 
                 nxd = all_orbit%nx(i) + all_orbit%nx(l) - all_orbit%nx(d)
                 nyd = all_orbit%ny(i) + all_orbit%ny(l) - all_orbit%ny(d) 
                 nzd = all_orbit%nz(i) + all_orbit%nz(l) - all_orbit%nz(d) 
                 tzd = all_orbit%itzp(i) + all_orbit%itzp(l) - all_orbit%itzp(d)
                 
                 if ( abs(nxd) > nkxmax ) cycle
                 if ( abs(nyd) > nkymax ) cycle
                 if ( abs(nzd) > nkzmax ) cycle
                 if ( abs(tzd) > 1 ) cycle
                 
                 do szd = -1, 1 ,2 
              
                    e = orbitof(nxd, nyd, nzd, szd, tzd)
                    if ( e <= below_ef ) cycle 
                    bra2 = pp_hhpp%ival(d,e)
                    if ( bra2 == 0 ) cycle 
                    if ( e >= d ) cycle
                    sum1 = sum1 + t2_ccm(channel2)%val(bra2,ket2)*chiral_3nf_asym(l,k,b,d,e,c)
                    
                 end do
              end do
           end do
           !$omp end do
           !$omp end parallel
           
           I4_ph(channel)%val(bra,ket) = I4_ph(channel)%val(bra,ket) + sum1

           sum1= 0.d0 
           !$omp parallel default(shared) private(l,m,nx2,ny2,nz2,tz2,channel2,ket2,nxd,nyd,nzd,tzd,szd,d,bra2)
           !$omp do schedule(dynamic), reduction(+:sum1)
           do l = 1, below_ef
              do m = l+1, below_ef
                 nx2 = all_orbit%nx(m) + all_orbit%nx(l)
                 ny2 = all_orbit%ny(m) + all_orbit%ny(l)
                 nz2 = all_orbit%nz(m) + all_orbit%nz(l)
                 tz2 = (all_orbit%itzp(m) + all_orbit%itzp(l))/2
                 channel2 = locate_channel(3,tz2, nx2, ny2, nz2) 
                 if ( channel2 == 0 ) cycle 
                 ket2 = hh_hhpp%ival(l,m)
                 if ( ket2 == 0 ) cycle
                 
              
                 
                 nxd = all_orbit%nx(m) + all_orbit%nx(l) - all_orbit%nx(b)
                 nyd = all_orbit%ny(m) + all_orbit%ny(l) - all_orbit%ny(b) 
                 nzd = all_orbit%nz(m) + all_orbit%nz(l) - all_orbit%nz(b) 
                 tzd = all_orbit%itzp(m) + all_orbit%itzp(l) - all_orbit%itzp(b)
                 
                 if ( abs(nxd) > nkxmax ) cycle
                 if ( abs(nyd) > nkymax ) cycle
                 if ( abs(nzd) > nkzmax ) cycle
                 if ( abs(tzd) > 1 ) cycle
                 
                 do szd = -1, 1 ,2 
              
                    d = orbitof(nxd, nyd, nzd, szd, tzd)
                    if ( d <= below_ef ) cycle 
                    bra2 = pp_hhpp%ival(d,b)
                    if ( bra2 == 0 ) cycle 
                    
                    sum1 = sum1 - t2_ccm(channel2)%val(bra2,ket2)*chiral_3nf_asym(l,k,m,d,i,c)
                    
                 end do
              end do
           end do
           !$omp end do
           !$omp end parallel
           
           I4_ph(channel)%val(bra,ket) = I4_ph(channel)%val(bra,ket) + sum1
        end do
     end do
  end do
  
  
  call mpi_barrier(mpi_comm_world,ierror)
  endwtime = MPI_WTIME()
  if ( iam == 0 ) write(6,*) 'Time for diag7/8', endwtime - startwtime



!!$  do channel   = my_ph_hphp_channel_low(iam), my_ph_hphp_channel_high(iam) 
!!$     !do channel   = 1, channels%number_ph_hphp_confs 
!!$     nx = channels%ph_hphp_quantum_numbers(channel*4)
!!$     ny = channels%ph_hphp_quantum_numbers(channel*4-1)
!!$     nz = channels%ph_hphp_quantum_numbers(channel*4-2)
!!$     tz = channels%ph_hphp_quantum_numbers(channel*4-3)
!!$     
!!$     !if ( check_my_channel_ph_hphp(channel) == 0 ) cycle
!!$     local_ch_id = channel
!!$     !
!!$     ! ket side if fully stored on each proc
!!$     !
!!$     bra_min = mapping_ph_hphp(iam+1,local_ch_id,2)
!!$     bra_max = mapping_ph_hphp(iam+1,local_ch_id,3)
!!$     !
!!$     ! bra side is distributed 
!!$     !
!!$     ket_min = mapping_ph_hphp(iam+1,local_ch_id,4)
!!$     ket_max = mapping_ph_hphp(iam+1,local_ch_id,5)
!!$     
!!$     channel2 = locate_ph_channel(1,tz, nx, ny, nz)
!!$     if ( channel2 == 0 ) cycle 
!!$     
!!$     dim1 = size( t2_ph_ccm(channel2)%val,1 )
!!$     dim2 = bra_max-bra_min+1 
!!$     dim3 = size( t2_ph_ccm(channel2)%val, 2) 
!!$     
!!$     call ZGEMM ( 'n', 'n', dim1, dim2, dim3, dcmplx(1.d0,0.d0), t2_ph_ccm(channel2)%val, & 
!!$          dim1,I4_ph(channel)%val(:,bra_min:bra_max), dim3, & 
!!$          dcmplx(1.d0,0.d0), t2_ph_ccm_eqn(channel2)%val(:,bra_min:bra_max), dim1 )
!!$     
!!$  end do

  
end SUBROUTINE t2diags_v3t2_linear
