!
! This routine has been checked with Davids, and is correct.
!
SUBROUTINE t2_intermediate
  USE single_particle_orbits
  USE configurations
  USE constants
  use t2_storage
  use one_body_operators
  use parallel 
  use chiral_potentials 
  
  IMPLICIT NONE
  INTEGER :: i,j,k,l,a,b,c,d, nx,ny,nz,sz,tz
  INTEGER :: bra,ket, channel, bra2,ket2, bra3,ket3,dim1, dim2, dim3, number_channels
  COMPLEX*16 ::  sum1, sum2 
  INTEGER :: nxd, nyd, nzd, szd, tzd, nxc, nyc, nzc, szc, tzc, nxk, nyk, nzk, szk, tzk
  COMPLEX*16, ALLOCATABLE :: temp1(:,:), temp2(:,:), temp(:,:), t2_tmp1(:), t2_tmp2(:)
  real*8  ::  startwtime , endwtime
  INTEGER :: nx2, ny2, nz2, sz2, tz2, channel2, channel3
  INTEGER :: nx3, ny3, nz3, sz3, tz3, channel4, ndim, i1
  INTEGER :: bra_min, bra_max, ket_min, ket_max, local_ch_id
  !
  ! set t2 amps to zero
  !
  do channel = 1, channels%number_hhpp_confs 
     t2_ccm_eqn(channel)%val = 0.d0
  end do
    
  
  call setup_I6I7_intermediate
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
  
  
  startwtime = MPI_WTIME()
  
  number_channels = channels%number_hhhh_confs  
  allocate( I9(number_channels) )
  do channel   = 1, channels%number_hhhh_confs 
     dim1 = size(  lookup_hhhh_configs(1,channel)%ival, 2) 
     allocate( i9(channel)%val(dim1, dim1) ) 
     i9(channel)%val = 0.d0 
  end do
  

  
  do channel   = 1, channels%number_hhhh_confs 
     nx = channels%hhhh_quantum_numbers(channel*4)
     ny = channels%hhhh_quantum_numbers(channel*4-1)
     nz = channels%hhhh_quantum_numbers(channel*4-2)
     tz = channels%hhhh_quantum_numbers(channel*4-3)
     
     do bra = 1, size(  lookup_hhhh_configs(1,channel)%ival, 2) 
        k= lookup_hhhh_configs(1,channel)%ival(1,bra) 
        l = lookup_hhhh_configs(1,channel)%ival(2,bra) 
        ket3 = hh_hhpp%ival(k,l)
        if ( ket3 == 0 ) cycle 
        
        do ket = 1, size(  lookup_hhhh_configs(1,channel)%ival, 2) 
           i = lookup_hhhh_configs(1,channel)%ival(1,ket)
           j = lookup_hhhh_configs(1,channel)%ival(2,ket) 
           
           !i9(channel)%val(bra,ket) = chiral_pot(k,l,i,j)
           !i9(channel)%val(bra,ket) = vnn_hhhh(channel)%val(bra,ket) 
           
           nx2 = all_orbit%nx(i) + all_orbit%nx(j)
           ny2 = all_orbit%ny(i) + all_orbit%ny(j)
           nz2 = all_orbit%nz(i) + all_orbit%nz(j)
           tz2 = (all_orbit%itzp(i) + all_orbit%itzp(j))/2
           sz2 = (all_orbit%szp(i) + all_orbit%szp(j))/2
           
           channel2 = locate_channel(3,tz2, nx2, ny2, nz2) 
           if ( channel2 == 0 ) cycle 
           
           if ( check_my_channel_hhpp(channel2) == 0 ) cycle
           local_ch_id = channel2
           !
           ! ket side if fully stored on each proc
           !
           bra_min = mapping_hhpp(iam+1,local_ch_id,2)
           bra_max = mapping_hhpp(iam+1,local_ch_id,3)
                      
           ket2 = hh_hhpp%ival(i,j)
           sum1 = 0.d0 
           do bra2 = bra_min, bra_max !1, size(  lookup_hhpp_configs(2,channel2)%ival, 2) 
              c = lookup_hhpp_configs(2,channel2)%ival(1,bra2)
              d = lookup_hhpp_configs(2,channel2)%ival(2,bra2) 
              
              !sum1 = sum1 + 0.5d0 * chiral_pot(k,l,c,d)*t2_ccm(channel2)%val(bra2,ket2)
              sum1 = sum1 + 0.5d0 * vnn_hhpp(channel2)%val(ket3,bra2)*t2_ccm(channel2)%val(bra2,ket2)
           end do
           
           i9(channel)%val(bra,ket) = i9(channel)%val(bra,ket) + sum1 
           
        end do
     end do
  end do

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
     
     call mpi_reduce(t2_tmp1,t2_tmp2,ndim,mpi_complex16,mpi_sum, &
          master,mpi_comm_world,ierror)
     t2_tmp1 = 0.d0 
     t2_tmp1 = t2_tmp2 
     call mpi_bcast(t2_tmp1,ndim,mpi_complex16,master, &
          mpi_comm_world,ierror)
     
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
  
  do channel   = 1, channels%number_hhhh_confs 
     i9(channel)%val =i9(channel)%val + vnn_hhhh(channel)%val
  end do
  
  !
  ! I8 and I9 contribution to t2-amps
  !
!!$  do channel   = 1, channels%number_hhpp_confs 
!!$     nx = channels%hhpp_quantum_numbers(channel*4)
!!$     ny = channels%hhpp_quantum_numbers(channel*4-1)
!!$     nz = channels%hhpp_quantum_numbers(channel*4-2)
!!$     tz = channels%hhpp_quantum_numbers(channel*4-3)
!!$     
!!$     !channel2 = locate_channel(6,tz, nx, ny, nz) 
!!$     channel2 = locate_channel(3,tz,  nx, ny, nz) 
!!$     if ( channel2 == 0 ) cycle 
!!$     
!!$     do ket = 1, size(  lookup_hhpp_configs(1,channel)%ival, 2) 
!!$        i= lookup_hhpp_configs(1,channel)%ival(1,ket) 
!!$        j = lookup_hhpp_configs(1,channel)%ival(2,ket) 
!!$        
!!$        do bra = 1, size(  lookup_hhpp_configs(2,channel)%ival, 2) 
!!$           a = lookup_hhpp_configs(2,channel)%ival(1,bra)
!!$           b = lookup_hhpp_configs(2,channel)%ival(2,bra) 
!!$           
!!$           bra3 = pp_pppp%ival(a,b) 
!!$           if ( bra3 == 0 ) cycle 
!!$             
!!$           sum1 = 0.d0 
!!$           do bra2 = 1, size(  lookup_hhpp_configs(2,channel)%ival, 2) 
!!$              c = lookup_hhpp_configs(2,channel)%ival(1,bra2)
!!$              d = lookup_hhpp_configs(2,channel)%ival(2,bra2) 
!!$              ket3 = pp_pppp%ival(c,d)
!!$              if ( ket3 == 0 ) cycle 
!!$              !sum1 = sum1 + 0.5d0* chiral_pot(a,b,c,d)*t2_ccm(channel)%val(bra2,ket)
!!$              sum1 = sum1 + 0.5d0* vnn_pppp(channel2)%val(bra3,ket3)*t2_ccm(channel)%val(bra2,ket)
!!$           end do
!!$           t2_ccm_eqn(channel)%val(bra,ket) = t2_ccm_eqn(channel)%val(bra,ket) + sum1 
!!$           
!!$        end do
!!$     end do
!!$  end do

  
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
     
     dim1 = size( vnn_pppp(channel2)%val,1 ) 
     dim2 = size( t2_ccm(channel2)%val, 2)
     dim3 = size( vnn_pppp(channel2)%val, 2) 
     
     call ZGEMM ( 'n', 'n', dim1, dim2, dim3, DCMPLX(0.5d0,0.D0), vnn_pppp(channel2)%val(bra_min:bra_max,ket_min:ket_max), & 
          dim1,t2_ccm(channel2)%val(ket_min:ket_max,:), dim3, & 
          DCMPLX(1.d0,0.D0), t2_ccm_eqn(channel2)%val(bra_min:bra_max,:), dim1 )
     
  end do
  !stop
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
     
     do bra = bra_min, bra_max !1, size(  lookup_hhpp_configs(2,channel)%ival, 2) 
        a = lookup_hhpp_configs(2,channel)%ival(1,bra)
        b = lookup_hhpp_configs(2,channel)%ival(2,bra) 
        
        do ket = 1, size(  lookup_hhpp_configs(1,channel)%ival, 2) 
           i= lookup_hhpp_configs(1,channel)%ival(1,ket) 
           j = lookup_hhpp_configs(1,channel)%ival(2,ket) 
                             

           sum2 = 0.d0 
           do ket2 = 1, size(  lookup_hhpp_configs(1,channel)%ival, 2) 
              k = lookup_hhpp_configs(1,channel)%ival(1,ket2) 
              l = lookup_hhpp_configs(1,channel)%ival(2,ket2) 
              !sum2 = sum2 + 0.5d0* chiral_pot(k,l,i,j)*t2_ccm(channel)%val(bra,ket2)
              
              nx3 = all_orbit%nx(l) + all_orbit%nx(k)
              ny3 = all_orbit%ny(l) + all_orbit%ny(k)
              nz3 = all_orbit%nz(l) + all_orbit%nz(k)
              tz3 = (all_orbit%itzp(l) + all_orbit%itzp(k))/2
              channel4 = locate_channel(1,tz3, nx3, ny3, nz3) 
              ket3 = hh_hhhh%ival(i,j)
              bra3 = hh_hhhh%ival(k,l)
              if ( channel4 == 0 ) cycle 
              
              !sum2 = sum2 + 0.5d0* chiral_pot(k,l,i,j)*t2_ccm(channel)%val(bra,ket2)
              sum2 = sum2 + 0.5d0* i9(channel4)%val(bra3,ket3)*t2_ccm(channel)%val(bra,ket2)
           end do
           
           t2_ccm_eqn(channel)%val(bra,ket) = t2_ccm_eqn(channel)%val(bra,ket) + sum2
           
        end do
     end do
  end do

  endwtime = MPI_WTIME()
  !if ( iam == 0 ) write(6,*) 'Total execution time for I8 v1', endwtime - startwtime


!!$
!!$  startwtime = MPI_WTIME()

!!$  endwtime = MPI_WTIME()
!!$  if ( iam == 0 ) write(6,*) 'Total execution time for I8 v2', endwtime - startwtime
!!$  
 
  
  
  number_channels = channels%number_hphp_confs  
  allocate( I4(number_channels) )
  do channel   = 1, channels%number_hphp_confs 
     dim1 = size(  lookup_hphp_configs(1,channel)%ival, 2) 
     allocate( i4(channel)%val(dim1, dim1) ) 
     i4(channel)%val = 0.d0 
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
     
     do bra = bra_min, bra_max !1, size(  lookup_hhpp_configs(2,channel)%ival, 2) 
        d = lookup_hhpp_configs(2,channel)%ival(1,bra)
        b = lookup_hhpp_configs(2,channel)%ival(2,bra) 
        
        do ket = 1, size(  lookup_hhpp_configs(1,channel)%ival, 2) 
           l= lookup_hhpp_configs(1,channel)%ival(1,ket) 
           j = lookup_hhpp_configs(1,channel)%ival(2,ket) 
  
           do k = 1, below_ef
              do szc = -1, 1, 2
                 nxc = all_orbit%nx(k) + all_orbit%nx(b) - all_orbit%nx(j)
                 nyc = all_orbit%ny(k) + all_orbit%ny(b) - all_orbit%ny(j) 
                 nzc = all_orbit%nz(k) + all_orbit%nz(b) - all_orbit%nz(j) 
                 tzc = all_orbit%itzp(k) + all_orbit%itzp(b) - all_orbit%itzp(j)
                 if ( abs(nxc) > nkxmax ) cycle 
                 if ( abs(nyc) > nkymax ) cycle 
                 if ( abs(nzc) > nkzmax ) cycle 
                 if ( abs(tzc) > 1 ) cycle 
                 c = orbitof(nxc, nyc, nzc, szc, tzc)
                 if  ( c <= below_ef ) cycle 
                 
                 nx2 = all_orbit%nx(k) + all_orbit%nx(b)
                 ny2 = all_orbit%ny(k) + all_orbit%ny(b)
                 nz2 = all_orbit%nz(k) + all_orbit%nz(b)
                 tz2 = (all_orbit%itzp(k) + all_orbit%itzp(b))/2
                 channel2 = locate_channel(4,tz2, nx2, ny2, nz2) 
                 if ( channel2 == 0 ) cycle 
                            
                 ket2 = hp_hphp%ival(k,b)
                 bra2 = hp_hphp%ival(j,c)
              
                 nx3 = all_orbit%nx(d) + all_orbit%nx(c)
                 ny3 = all_orbit%ny(d) + all_orbit%ny(c)
                 nz3 = all_orbit%nz(d) + all_orbit%nz(c)
                 tz3 = (all_orbit%itzp(d) + all_orbit%itzp(c))/2
                 channel3 = locate_channel(3,tz3, nx3, ny3, nz3) 
                 if ( channel3 == 0 ) cycle 
                 ket3 = hh_hhpp%ival(k,l)
                 bra3 = pp_hhpp%ival(c,d)

                 i4(channel2)%val(ket2,bra2) = i4(channel2)%val(ket2,bra2) -  & 
                      0.5d0 * vnn_hhpp(channel3)%val(ket3,bra3)*t2_ccm(channel)%val(bra,ket)
                 
              end do
           end do
        end do
     end do
  end do

  number_channels = size( i4 )
  do channel = 1, number_channels 
     bra = size(i4(channel)%val,1)
     ket = size(i4(channel)%val,2)
     
     ndim = bra * ket
     allocate( t2_tmp1(ndim) )
     allocate( t2_tmp2(ndim) )
     t2_tmp1 = 0.d0
     t2_tmp2 = 0.d0

     i1 = 0 
     do i =1, bra
        do j = 1, ket
           i1 = i1 + 1
           t2_tmp1(i1) = i4(channel)%val(i,j)
        end do
     end do
     
     call mpi_reduce(t2_tmp1,t2_tmp2,ndim,mpi_complex16,mpi_sum, &
          master,mpi_comm_world,ierror)
     t2_tmp1 = 0.d0 
     t2_tmp1 = t2_tmp2 
     call mpi_bcast(t2_tmp1,ndim,mpi_complex16,master, &
          mpi_comm_world,ierror)
     
     i4(channel)%val = 0.d0
     i1 = 0 
     do i =1, bra
        do j = 1, ket
           i1 = i1 + 1
           i4(channel)%val(i,j) = t2_tmp1(i1)
        end do
     end do
     
     deallocate( t2_tmp1, t2_tmp2 ) 
  end do
  
  
  do channel   = 1, channels%number_hphp_confs 
     i4(channel)%val = i4(channel)%val + vnn_hphp(channel)%val
  end do
  
  !
  ! I4 contribution to t2-amps
  !
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
     
     do bra = bra_min, bra_max !1, size(  lookup_hhpp_configs(2,channel)%ival, 2) 
        a = lookup_hhpp_configs(2,channel)%ival(1,bra)
        b = lookup_hhpp_configs(2,channel)%ival(2,bra) 
        
        do ket = 1, size(  lookup_hhpp_configs(1,channel)%ival, 2) 
           i= lookup_hhpp_configs(1,channel)%ival(1,ket) 
           j = lookup_hhpp_configs(1,channel)%ival(2,ket) 
        
           
           sum2 = 0.d0 
           do k = 1, below_ef
              do szc = -1, 1, 2
                 nxc = all_orbit%nx(k) + all_orbit%nx(b) - all_orbit%nx(j)
                 nyc = all_orbit%ny(k) + all_orbit%ny(b) - all_orbit%ny(j) 
                 nzc = all_orbit%nz(k) + all_orbit%nz(b) - all_orbit%nz(j) 
                 tzc = all_orbit%itzp(k) + all_orbit%itzp(b) - all_orbit%itzp(j)
                 if ( abs(nxc) > nkxmax ) cycle 
                 if ( abs(nyc) > nkymax ) cycle 
                 if ( abs(nzc) > nkzmax ) cycle 
                 if ( abs(tzc) > 1 ) cycle 
                 c = orbitof(nxc, nyc, nzc, szc, tzc)
                 if ( c <= below_ef ) cycle 

                 nx2 = all_orbit%nx(i) + all_orbit%nx(k)
                 ny2 = all_orbit%ny(i) + all_orbit%ny(k)
                 nz2 = all_orbit%nz(i) + all_orbit%nz(k)
                 tz2 = (all_orbit%itzp(i) + all_orbit%itzp(k))/2
                 channel2 = locate_channel(3,tz2, nx2, ny2, nz2) 
                 if ( channel2 == 0 ) cycle 
                 ket2 = hh_hhpp%ival(i,k)
                 bra2 = pp_hhpp%ival(a,c)
             
                 nx3 = all_orbit%nx(b) + all_orbit%nx(k)
                 ny3 = all_orbit%ny(b) + all_orbit%ny(k)
                 nz3 = all_orbit%nz(b) + all_orbit%nz(k)
                 tz3 = (all_orbit%itzp(b) + all_orbit%itzp(k))/2
                 channel4 = locate_channel(4,tz3, nx3, ny3, nz3) 
                 ket3 = hp_hphp%ival(j,c)
                 bra3 = hp_hphp%ival(k,b)
                 if ( channel4 == 0 ) cycle 
              
                 !sum2 = sum2 - chiral_pot(k,b,j,c)*t2_ccm(channel2)%val(bra2,ket2)
                 sum2 = sum2 - I4(channel4)%val(bra3,ket3)*t2_ccm(channel2)%val(bra2,ket2)
              end do
           end do
           
           t2_ccm_eqn(channel)%val(bra,ket) = t2_ccm_eqn(channel)%val(bra,ket) + sum2 
           bra2 = pp_hhpp%ival(b,a)
           t2_ccm_eqn(channel)%val(bra2,ket) = t2_ccm_eqn(channel)%val(bra2,ket) - sum2 
           ket2 = hh_hhpp%ival(j,i)
           t2_ccm_eqn(channel)%val(bra,ket2) = t2_ccm_eqn(channel)%val(bra,ket2) - sum2 
           t2_ccm_eqn(channel)%val(bra2,ket2) = t2_ccm_eqn(channel)%val(bra2,ket2) + sum2 
           
        end do
     end do
  end do

  
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
           t2_ccm_eqn(channel)%val(i,j) = t2_tmp1(i1)
        end do
     end do
     
     deallocate( t2_tmp1, t2_tmp2 ) 
  end do

  deallocate( i6, i7 ) 
  do i = 1, size(i4) 
     deallocate(i4(i)%val )
  ENDDO
  deallocate( i4 ) 
  
  do i = 1, size(i9)
     deallocate(i9(i)%val)
  end do
  deallocate( i9 )

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
  
  INTEGER :: i,a,i1, c,d,k,l,bra,bra2,ket, ndim2
  COMPLEX*16 :: gmat,t2, sum1
  INTEGER :: bra_min, bra_max, ket_min, ket_max, local_ch_id
  COMPLEX*16, allocatable :: t2_tmp1(:), t2_tmp2(:) 
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
  !
  ! Add fock matrix to I5,  I7 and I6 
  !
  
  !
  ! I7(b) intermediate
  !
  do i1 = 1, channels%number_hhpp_confs 
     
     if ( check_my_channel_hhpp(i1) == 0 ) cycle
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
           
           
           do k = 1, below_ef 
              if ( all_orbit%nx(k) /= all_orbit%nx(i) ) cycle 
              if ( all_orbit%ny(k) /= all_orbit%ny(i) ) cycle
              if ( all_orbit%nz(k) /= all_orbit%nz(i) ) cycle
              if ( all_orbit%itzp(k) /= all_orbit%itzp(i) ) cycle
              
              !k = i 
              !gmat =  chiral_pot(k,l,c,d)
              bra2 = hh_hhpp%ival(k,l)
              if ( bra2 == 0 ) cycle 
              !write(6,*) bra2, size(vnn_hhpp(i1)%val, 1) 
              gmat = vnn_hhpp(i1)%val(bra2,ket) 
              !gmat = 0.d0 
              i7(k,i) =  i7(k,i) + 0.5d0 * gmat*t2_ccm(i1)%val(ket,bra) 
              
           end do
        end DO
     end DO
  end do


  !
  ! I6(abc) intermediate
  ! Here I6(b) is calculated, I6(a+c) calculated in t1-eqn
  !
  do i1 = 1, channels%number_hhpp_confs 
     
     if ( check_my_channel_hhpp(i1) == 0 ) cycle
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
           do c = below_ef+1, tot_orbs
              if ( all_orbit%nx(a) /= all_orbit%nx(c) ) cycle 
              if ( all_orbit%ny(a) /= all_orbit%ny(c) ) cycle
              if ( all_orbit%nz(a) /= all_orbit%nz(c) ) cycle
              if ( all_orbit%itzp(a) /= all_orbit%itzp(c) ) cycle
              
              
              bra2 = pp_hhpp%ival(d,c)
              if ( bra2 == 0 ) cycle 
              !c = a   
              !gmat = chiral_pot(k,l,d,c)
              gmat = vnn_hhpp(i1)%val(ket,bra2)
              sum1 = - 0.5d0* gmat * t2 
              ! XXXXXXXXXXXXXXXXXXXXXXXX
              I6(a,c) = I6(a,c) + sum1
           end do
        end DO
     end DO
  end do
  
  bra = tot_orbs
  ket = tot_orbs
  
  ndim2 = bra * ket
  allocate( t2_tmp1(ndim2) )
  allocate( t2_tmp2(ndim2) )
  t2_tmp1 = 0.d0
  t2_tmp2 = 0.d0

  i1 = 0
  do a = 1, tot_orbs
     do c = 1, tot_orbs
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
        i1 = i1 + 1
        if ( a > below_ef .and. c > below_ef )   I6(a,c) = t2_tmp1(i1)
        if ( a <= below_ef .and. c <= below_ef ) I7(a,c) = t2_tmp1(i1)
     end do
  end do
  deallocate( t2_tmp1, t2_tmp2 )
  
  
  I7 = I7 + fock_mtx(1:below_ef, 1:below_ef) 
  I6ac = fock_mtx(below_ef+1:tot_orbs, below_ef+1:tot_orbs) 
  
  I6 = i6 + I6ac

  
  deallocate( i6ac )

end SUBROUTINE setup_I6I7_intermediate
