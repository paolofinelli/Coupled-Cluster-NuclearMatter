!
! This routine has been checked with Davids, and is correct.
!
SUBROUTINE t2_intermediate
  
  USE single_particle_orbits
  USE configurations
  USE constants
  use gmat_storage
  use one_body_operators
  use ang_mom_functions

  IMPLICIT NONE
  TYPE (configuration_descriptor) :: bra_configs
  TYPE (configuration_descriptor) :: ket_configs
  
  TYPE (configuration_descriptor) :: t2_bra_configs 
  TYPE (configuration_descriptor) :: t2_ket_configs 
  INTEGER :: t2_channel, number_t2_channels, channel, number_channels, bra3
  INTEGER :: i,j,a,b,c,k,l,d,bra,ket, bra2,ket2, i1, ang_mom, p_parity, isospin_z
  INTEGER :: I3_channel, I8_channel, I9_channel
  REAL*8 ::  fnorm, t2_old, cross, t2_angmom_fact, ang_mom_factor, angmom_fact
  REAL*8 :: sum1, sum2, phase, t2, gmat, phase_ab, phase_ij, const, phase_cd
  INTEGER :: ph_channel, ph_channel1, ph_channel2, jj, ji, ja, jb,ket3 
  INTEGER :: bra_conf, ket_conf, jt, j_min, j_max, iph, jc, jk
  INTEGER :: tza, tzb, tzi, tzj ,tzc, tzk, tz_ac, tz_ik, tz_kb, tz_cj, li, lj, la, lb, lc, lk
  INTEGER :: ipar_ac, ipar_cj, ph_channel3, ipar_ik, ipar_kb, I10_channel, I11_channel, gmat_channel
  INTEGER :: I4_channel, ipar, tz, channel1, channel2, channel3, bra4, ket4, I11b_channel, dim1, dim2, dim3
  INTEGER :: ld1, ld2, ld3
  LOGICAL :: triag
  REAL*4 :: TM1, TM2, tm3
  real*8, allocatable :: amat(:,:), bmat(:,:), t2_tmp1(:), t2_tmp2(:)
  INTEGER :: bra_min, bra_max, ket_min, ket_max, local_ch_id, ndim
  integer, allocatable :: nconf_low(:), nconf_high(:)
  integer :: nconfs, nconfs_tot, number_mtxel_iam, diff, n
  integer(8) :: count, count1
  real*8  ::  startwtime , endwtime
  integer, parameter :: chunk = 100
  
  count = 0 
  !
  ! set t2 amps to zero
  !
  do i1 = 1, channels%number_confs 
     ang_mom   = channels%config_J_ipar_tz(i1*3)
     p_parity  = channels%config_J_ipar_tz(i1*3-1)
     isospin_z = channels%config_J_ipar_tz(i1*3-2)


     t2_channel = locate_gmatchannel(3,isospin_z, p_parity, ang_mom)
     if ( t2_channel == 0  ) cycle
     t2_ccm_eqn%depot(t2_channel)%val = 0.d0
  end do
  
  do channel = 1, ph_channels%number_confs 
     ang_mom   = ph_channels%config_J_ipar_tz(channel*3)
     p_parity  = ph_channels%config_J_ipar_tz(channel*3-1)
     isospin_z = ph_channels%config_J_ipar_tz(channel*3-2)

     if ( locate_crosscoupled_channel(1, isospin_z, p_parity, ang_mom) == 0 ) cycle
     ph_channel1 = locate_crosscoupled_channel(1, isospin_z, p_parity, ang_mom ) 
     t2_cross_ccm_eqn%depot(ph_channel1)%val = 0.d0
     t2_cross_ccm%depot(ph_channel1)%val = 0.d0
  end do
  
  
  !
  ! Setup cross-coupled matrix elements gmat_ph_hhpp%depot(:)%val(:,:)
  ! 
  ! Coupling direction: ingoing -> outgoing 
  ! <ij||ab> (ai->J) and (bj->J) if one changes coupling direction 
  ! to outgoing -> ingoing one gets a phase (-1)**(ji+jj+ja+jb - 2J ) 
  ! 
  do channel   = 1, ph_channels%number_confs 
     ang_mom   = ph_channels%config_J_ipar_tz(channel*3)
     p_parity  = ph_channels%config_J_ipar_tz(channel*3-1)
     isospin_z = ph_channels%config_J_ipar_tz(channel*3-2)
     
     if ( locate_crosscoupled_channel(1, isospin_z, p_parity, ang_mom) == 0 ) cycle
     ph_channel = locate_crosscoupled_channel(1, isospin_z, p_parity, ang_mom ) 
     
     allocate( amat( size(t2_cross_ccm%depot(ph_channel)%val,1 ),  & 
          size(t2_cross_ccm%depot(ph_channel)%val,2 ) ) )
     amat = 0.d0 
     
     !$omp parallel default(shared) private(bra,a,i,ja,ji,ket,j,b,jj,jb,tz,ipar,j_min,j_max, & 
     !$omp fnorm,jt,t2_channel,bra2,ket2,ang_mom_factor,t2_old)
     !$omp do schedule(static), reduction(+:amat) 
     DO bra=1, size(  lookup_ph_configs(3,channel)%ival, 2) 
        a = lookup_ph_configs(3,channel)%ival(1,bra) 
        i = lookup_ph_configs(3,channel)%ival(2,bra) 

        ja = all_orbit%jj(a)
        ji = all_orbit%jj(i)

        DO ket=1, size(  lookup_ph_configs(2,channel)%ival, 2) 
           j = lookup_ph_configs(2,channel)%ival(1,ket) 
           b = lookup_ph_configs(2,channel)%ival(2,ket) 

           jj=all_orbit%jj(j)
           jb=all_orbit%jj(b)
                         
           tz = (all_orbit%itzp(i) + all_orbit%itzp(j))/2
           ipar = ( 1 - (-1)**( all_orbit%ll(i) + all_orbit%ll(j) ))/2
           
           
           j_min=max( ABS((all_orbit%jj(i)-all_orbit%jj(j))/2), ABS((all_orbit%jj(a)-all_orbit%jj(b))/2) )
           j_max=min( (all_orbit%jj(i)+all_orbit%jj(j))/2, (all_orbit%jj(a)+all_orbit%jj(b))/2 ) 
                 
           !
           ! Additional phase since we couple in reverse directon (ia->J) and (jb->J)
           !
           
           
           fnorm=iph((ja+jj)/2+ang_mom) 
           
           IF(j_min > j_max) CYCLE
           !cross=0.
           DO jt=j_min,j_max
              
              IF(( i == j ).AND.(MOD(jt,2)/=0)) CYCLE
              IF(( a == b ).AND.(MOD(jt,2)/=0)) CYCLE
              IF(triag(all_orbit%jj(i),all_orbit%jj(j),2*jt)) CYCLE
              IF(triag(all_orbit%jj(a),all_orbit%jj(b),2*jt)) CYCLE
              t2_channel = locate_gmatchannel(3,tz,ipar,jt)
              if ( t2_channel == 0 ) cycle
              bra2 = pp_hhpp(t2_channel)%ival(a,b)
              ket2 = hh_hhpp(t2_channel)%ival(i,j)
              
              if ( bra2 == 0 ) cycle
              if ( ket2 == 0 ) cycle
           
              ang_mom_factor=sjs(ji,ja,2*ang_mom, jb,jj,2*jt)* &
                   (2.*jt+1.)*fnorm*iph(jt)
              t2_old = t2_ccm%depot(t2_channel)%val(bra2,ket2) 
              
              amat(bra,ket) = amat(bra,ket) + ang_mom_factor*t2_old
           end DO
           
        end do
     end do
     !$omp end do                                                                                                                                      
     !$omp end parallel 

     t2_cross_ccm%depot(ph_channel)%val = amat 
     deallocate ( amat ) 
     
  end DO

  
  !
  ! I4(a) 
  !
  number_channels = size(gmat_hp_hphp%depot)
  allocate( I4abc(1:number_channels) )
  
  do i=1, number_channels
     bra = size(gmat_hp_hphp%depot(i)%val,1)
     ket = size(gmat_hp_hphp%depot(i)%val,2)
     
     !allocate( I4abc(i)%val(0:bra-1,0:ket-1) )
     allocate( I4abc(i)%val(bra,ket) )
     I4abc(i)%val = 0.d0
     
  end do
  
  
  !
  ! allocate and setup matrix for I4(b) intermediate
  !
  number_channels = size(gmat_hp_hpph%depot)
  allocate( temp_mat(1:number_channels) )
  do i = 1, number_channels 
     
     bra2 = size(gmat_hp_hpph%depot(i)%val,1)
     ket2 = size(gmat_hp_hpph%depot(i)%val,2)

     allocate( temp_mat(i)%val(bra2,ket2) )
     temp_mat(i)%val = 0.d0
     
  end do
  

  do channel   = 1, ph_channels%number_confs 
     ang_mom   = ph_channels%config_J_ipar_tz(channel*3)
     p_parity  = ph_channels%config_J_ipar_tz(channel*3-1)
     isospin_z = ph_channels%config_J_ipar_tz(channel*3-2)
     
     channel1 = locate_crosscoupled_channel(4,isospin_z,p_parity,ang_mom)
     if ( channel1 == 0 ) cycle
     
     
     DO ket=1, size(  lookup_ph_configs(3,channel)%ival, 2) 
        b = lookup_ph_configs(3,channel)%ival(1,ket) 
        j = lookup_ph_configs(3,channel)%ival(2,ket) 
           
        phase = iph( (all_orbit%jj(b) + all_orbit%jj(j))/2 )
        temp_mat(channel1)%val(:,ket) = gmat_hp_hpph%depot(channel1)%val(:,ket) * phase 
           
     end do

  end do

  !
  ! I4(b) (matmul) 
  !
  count = 0
  do channel   = 1, ph_channels%number_confs 
     ang_mom   = ph_channels%config_J_ipar_tz(channel*3)
     p_parity  = ph_channels%config_J_ipar_tz(channel*3-1)
     isospin_z = ph_channels%config_J_ipar_tz(channel*3-2)
     
     channel1 = locate_crosscoupled_channel(1,isospin_z,p_parity,ang_mom)
     channel2 = locate_crosscoupled_channel(4,isospin_z,p_parity,ang_mom)
     channel3 = locate_crosscoupled_channel(2,isospin_z,p_parity,ang_mom)

     
     
     
     if ( channel1*channel2*channel3 == 0 ) cycle
     if ( check_my_channel_hp_4a(channel1) == 0 ) cycle
          
               
     local_ch_id = channel1
          
     !
     ! bra side if fully stored on each proc
     !
     bra_min = mapping_hp_4a(iam+1,local_ch_id,4)
     bra_max = mapping_hp_4a(iam+1,local_ch_id,5)
     !
     ! ket side is distributed 
     !
     ket_min = mapping_hp_4a(iam+1,local_ch_id,2)
     ket_max = mapping_hp_4a(iam+1,local_ch_id,3)
     
     allocate( amat(bra_min:bra_max, ket_min:ket_max) )
     
     
     DO ket=ket_min, ket_max !1, size(  lookup_ph_configs(2,channel)%ival, 2) 
        j = lookup_ph_configs(2,channel)%ival(1,ket) 
        b = lookup_ph_configs(2,channel)%ival(2,ket) 
     
        DO bra=1, size(  lookup_ph_configs(3,channel)%ival, 2) 
           d = lookup_ph_configs(3,channel)%ival(1,bra) 
           l = lookup_ph_configs(3,channel)%ival(2,bra) 

           count = count + 1
           amat(bra, ket) = -0.5d0* iph(ang_mom) * t2_cross_ccm%depot(channel1)%val(bra, ket) - & 
                t1_ccm(d,j)*t1_ccm(b,l) * iph( (all_orbit%jj(b) + all_orbit%jj(j))/2 )
        end do
     end DO

     dim1 = size( temp_mat(channel2)%val,1 )
     !dim2 = size( t2_cross_ccm%depot(channel1)%val, 2)
     dim2 = size ( amat, 2 ) 
     dim3 = size( temp_mat(channel2)%val,2 )
     
     const = -0.5d0* iph(ang_mom) 
     

     ! Nov 10
     call DGEMM ( 'n', 'n', dim1, dim2, dim3, 1.d0, temp_mat(channel2)%val, & 
          dim1,amat, dim3, & 
          1.d0, I4abc(channel3)%val(:, ket_min:ket_max), dim1 )
     
     I4abc(channel3)%val(:,ket_min:ket_max) = I4abc(channel3)%val(:,ket_min:ket_max) + &
          gmat_hp_hphp%depot(channel3)%val(:, ket_min:ket_max) 

     deallocate ( amat ) 
  end do

  
  
  
  number_channels = size(gmat_hp_hpph%depot)
  do i = 1, number_channels 

     deallocate( temp_mat(i)%val )
  end do
  deallocate( temp_mat ) 

  !
  ! I4(c)
  ! Needs to be parallelized separately
  ! 

  number_channels = size( I4abc)
  allocate( temp_mat(1:number_channels) )
  do i = 1, number_channels 
     
     bra2 = size(I4abc(i)%val,1)
     ket2 = size(I4abc(i)%val,2)
     
     !allocate( temp_mat(i)%val(0:bra2-1,0:ket2-1) )
     allocate( temp_mat(i)%val(bra2,ket2) )
     temp_mat(i)%val = 0.d0
     
  end do
  
  
  do channel   = 1, ph_channels%number_confs 
     ang_mom   = ph_channels%config_J_ipar_tz(channel*3)
     p_parity  = ph_channels%config_J_ipar_tz(channel*3-1)
     isospin_z = ph_channels%config_J_ipar_tz(channel*3-2)

     
     channel1 = locate_crosscoupled_channel(2,isospin_z,p_parity,ang_mom)
     channel2 = locate_crosscoupled_channel(5,isospin_z,p_parity,ang_mom)
     channel3 = locate_crosscoupled_channel(1,isospin_z,p_parity,ang_mom)

     if ( channel1*channel2*channel3 == 0 ) cycle
     !local_ch_id = channel3
     !if ( check_my_channel_hp_4a(channel3) == 0 ) cycle
     local_ch_id = channel2
     if ( check_my_channel_hp_hppp(channel2) == 0 ) cycle
     
     
     
     !
     ! bra side if fully stored on each proc
     !
     bra_min = mapping_hp_hppp(iam+1,local_ch_id,4)
     bra_max = mapping_hp_hppp(iam+1,local_ch_id,5)
     !
     ! ket side is distributed 
     !
     ket_min = mapping_hp_hppp(iam+1,local_ch_id,2)
     ket_max = mapping_hp_hppp(iam+1,local_ch_id,3)
     
     DO ket=ket_min, ket_max !1, size(  lookup_ph_configs(2,channel)%ival, 2) 
        d = lookup_ph_configs(4,channel)%ival(1,ket) 
        b = lookup_ph_configs(4,channel)%ival(2,ket) 
     
        DO bra=1, size(  lookup_ph_configs(2,channel)%ival, 2) 
           k = lookup_ph_configs(2,channel)%ival(1,bra) 
           c = lookup_ph_configs(2,channel)%ival(2,bra) 

           
           jc = all_orbit%jj(c)
           jk = all_orbit%jj(k)
           
           
           bra2 = hp_hphp_cross(channel1)%ival(k,c)
           if ( bra2 == 0 ) cycle
                      
           do j = 1, below_ef 
                         
              ket2 = hp_hphp_cross(channel1)%ival(j,b)
              if ( ket2 == 0 ) cycle
              
              count = count + 1
              gmat = gmat_hp_hppp%depot(channel2)%val(bra,ket)
              sum1 = gmat * t1_ccm(d,j) 
              temp_mat(channel1)%val(bra2,ket2) = temp_mat(channel1)%val(bra2,ket2) + sum1 
              
           end do
        end DO
     end DO

     
  end do

  
  number_channels = size( I4abc )
  do t2_channel = 1, number_channels 
     bra = size(temp_mat(t2_channel)%val,1)
     ket = size(temp_mat(t2_channel)%val,2)
     
     ndim = bra * ket
     allocate( t2_tmp1(ndim) )
     allocate( t2_tmp2(ndim) )
     t2_tmp1 = 0.d0
     t2_tmp2 = 0.d0

     i1 = 0 
     do j =1, ket 
        do i = 1, bra
           i1 = i1 + 1
           t2_tmp1(i1) = temp_mat(t2_channel)%val(i,j) 
        end do
     end do
     
     call mpi_reduce(t2_tmp1,t2_tmp2,ndim,mpi_real8,mpi_sum, &
          master,mpi_comm_world,ierror)
     call mpi_bcast(t2_tmp2,ndim,mpi_real8,master, &
          mpi_comm_world,ierror)
     
     i1 = 0 
     do j = 1, ket
        do i =1, bra
           i1 = i1 + 1
           temp_mat(t2_channel)%val(i,j) = 0.d0
           temp_mat(t2_channel)%val(i,j) = t2_tmp2(i1)
        end do
     end do

     deallocate( t2_tmp1, t2_tmp2 ) 
  end do
  !
  ! add phase (-1)**( jk + jc - J ) to I4abc intermediate
  !
  do channel   = 1, ph_channels%number_confs 
     ang_mom   = ph_channels%config_J_ipar_tz(channel*3)
     p_parity  = ph_channels%config_J_ipar_tz(channel*3-1)
     isospin_z = ph_channels%config_J_ipar_tz(channel*3-2)
     
     channel1 = locate_crosscoupled_channel(2,isospin_z,p_parity,ang_mom)
     channel2 = locate_crosscoupled_channel(1,isospin_z,p_parity,ang_mom)
     if ( channel1*channel2 == 0 ) cycle
     if ( check_my_channel_hp_4a(channel2) == 0 ) cycle
     local_ch_id = channel2
     !
     ! ket side if fully stored on each proc
     !
     bra_min = mapping_hp_4a(iam+1,local_ch_id,4)
     bra_max = mapping_hp_4a(iam+1,local_ch_id,5)
     !
     ! bra side is distributed 
     !
     ket_min = mapping_hp_4a(iam+1,local_ch_id,2)
     ket_max = mapping_hp_4a(iam+1,local_ch_id,3)
     

     I4abc(channel1)%val(:,ket_min:ket_max) = I4abc(channel1)%val(:,ket_min:ket_max) + & 
          temp_mat(channel1)%val(:,ket_min:ket_max)
     
  end do


  do i = 1, size( temp_mat ) 
     deallocate( temp_mat(i)%val )
  end do
  deallocate( temp_mat ) 
  
  !
  ! I11(b) intermediate
  ! 
  do channel   = 1, ph_channels%number_confs 
     ang_mom   = ph_channels%config_J_ipar_tz(channel*3)
     p_parity  = ph_channels%config_J_ipar_tz(channel*3-1)
     isospin_z = ph_channels%config_J_ipar_tz(channel*3-2)
     
     channel1 = locate_crosscoupled_channel(2,isospin_z,p_parity,ang_mom)
     channel2 = locate_crosscoupled_channel(3,isospin_z,p_parity,ang_mom)
     channel3 = locate_crosscoupled_channel(1,isospin_z,p_parity,ang_mom)
     
     if ( channel1*channel2*channel3 == 0 ) cycle
     local_ch_id = channel3
     if ( check_my_channel_hp_4a(channel3) == 0 ) cycle
     
     
     !
     ! bra side if fully stored on each proc
     !
     bra_min = mapping_hp_4a(iam+1,local_ch_id,4)
     bra_max = mapping_hp_4a(iam+1,local_ch_id,5)
     !
     ! ket side is distributed 
     !
     ket_min = mapping_hp_4a(iam+1,local_ch_id,2)
     ket_max = mapping_hp_4a(iam+1,local_ch_id,3)
     
     
     DO ket=ket_min, ket_max !1, size(  lookup_ph_configs(2,channel)%ival, 2) 
        j = lookup_ph_configs(2,channel)%ival(1,ket) 
        b = lookup_ph_configs(2,channel)%ival(2,ket) 
     
        DO bra=1, size(  lookup_ph_configs(2,channel)%ival, 2) 
           k = lookup_ph_configs(2,channel)%ival(1,bra) 
           c = lookup_ph_configs(2,channel)%ival(2,bra) 

           
           jc = all_orbit%jj(c)
           jk = all_orbit%jj(k)
           
           bra2 = hp_hphh_cross(channel2)%ival(k,c)
           if ( bra2 == 0 ) cycle
                      
           sum1 = 0.d0
           do l = 1, below_ef 
              
              ket2 = hh_hphh_cross(channel2)%ival(j,l)
              if ( ket2 == 0 ) cycle
              
              count = count + 1
              gmat = gmat_hp_hphh%depot(channel2)%val(bra2,ket2)
              sum1 = sum1 - gmat * t1_ccm(b,l)
              
           end do
           
           I4abc(channel1)%val(bra,ket) = I4abc(channel1)%val(bra,ket) + sum1 

        end DO
     end DO
     
  end do
  !call mpi_barrier(mpi_comm_world,ierror)
  
  
  !
  ! add phase (-1)**( jk + jc - J ) to I4abc intermediate
  !
  do channel   = 1, ph_channels%number_confs 
     ang_mom   = ph_channels%config_J_ipar_tz(channel*3)
     p_parity  = ph_channels%config_J_ipar_tz(channel*3-1)
     isospin_z = ph_channels%config_J_ipar_tz(channel*3-2)
     
     channel1 = locate_crosscoupled_channel(2,isospin_z,p_parity,ang_mom)
     channel2 = locate_crosscoupled_channel(1,isospin_z,p_parity,ang_mom)
     if ( channel1*channel2 == 0 ) cycle
     if ( check_my_channel_hp_4a(channel2) == 0 ) cycle

     
     local_ch_id = channel2
     !
     ! ket side if fully stored on each proc
     !
     bra_min = mapping_hp_4a(iam+1,local_ch_id,4)
     bra_max = mapping_hp_4a(iam+1,local_ch_id,5)
     !
     ! bra side is distributed 
     !
     ket_min = mapping_hp_4a(iam+1,local_ch_id,2)
     ket_max = mapping_hp_4a(iam+1,local_ch_id,3)
     

     
     DO ket= 1, size(  lookup_ph_configs(2,channel)%ival, 2) 
        k = lookup_ph_configs(2,channel)%ival(1,ket) 
        c = lookup_ph_configs(2,channel)%ival(2,ket) 

        phase = iph( (all_orbit%jj(k) + all_orbit%jj(c))/2 - ang_mom )
        I4abc(channel1)%val(ket,ket_min:ket_max) = I4abc(channel1)%val(ket,ket_min:ket_max)*phase 
        !I4abc(channel1)%val(ket,:) = I4abc(channel1)%val(ket,:)*phase 
        
     end do

  end do

  
  


  !
  ! I4abc contriubtion to t2-amps 
  !
  do channel   = 1, ph_channels%number_confs 
     ang_mom   = ph_channels%config_J_ipar_tz(channel*3)
     p_parity  = ph_channels%config_J_ipar_tz(channel*3-1)
     isospin_z = ph_channels%config_J_ipar_tz(channel*3-2)

     channel1 = locate_crosscoupled_channel(1,isospin_z,p_parity,ang_mom)
     channel2 = locate_crosscoupled_channel(2,isospin_z,p_parity,ang_mom)
     
     if ( channel1*channel2 == 0 ) cycle
     if ( check_my_channel_hp_4a(channel1) == 0 ) cycle

     local_ch_id = channel1
     !
     ! ket side if fully stored on each proc
     !
     bra_min = 1 !mapping_hp_4a(iam+1,local_ch_id,4)
     bra_max = size( I4abc(channel2)%val, 1 ) !mapping_hp_4a(iam+1,local_ch_id,5)
     !
     ! bra side is distributed 
     !
     ket_min = mapping_hp_4a(iam+1,local_ch_id,2)
     ket_max = mapping_hp_4a(iam+1,local_ch_id,3)
     
     allocate( amat(bra_min:bra_max, ket_min:ket_max) )
     amat(:, ket_min:ket_max) = I4abc(channel2)%val(:, ket_min:ket_max) 
     
     dim1 = size( t2_cross_ccm%depot(channel1)%val,1 )
     dim2 = size( amat,2)
     dim3 = size( t2_cross_ccm%depot(channel1)%val, 2) 



     
     call DGEMM ( 'n', 'n', dim1, dim2, dim3, -1.d0, t2_cross_ccm%depot(channel1)%val, & 
          dim1,amat, dim3, & 
          1.d0, t2_cross_ccm_eqn%depot(channel1)%val(:,ket_min:ket_max), dim1 )
     
     deallocate ( amat ) 
  end do


  do i = 1, size( i4abc(:) )
     deallocate( i4abc(i)%val ) 
  end do
  deallocate( I4abc )



  number_channels = size( t2_cross_ccm_eqn%depot )
  do t2_channel = 1, number_channels 
     bra = size(t2_cross_ccm_eqn%depot(t2_channel)%val,1)
     ket = size(t2_cross_ccm_eqn%depot(t2_channel)%val,2)
     
     ndim = bra * ket
     allocate( t2_tmp1(ndim) )
     allocate( t2_tmp2(ndim) )
     t2_tmp1 = 0.d0
     t2_tmp2 = 0.d0

     i1 = 0 
     do j =1, ket 
        do i = 1, bra
           i1 = i1 + 1
           t2_tmp1(i1) = t2_cross_ccm_eqn%depot(t2_channel)%val(i,j) 
        end do
     end do
     
     call mpi_reduce(t2_tmp1,t2_tmp2,ndim,mpi_real8,mpi_sum, &
          master,mpi_comm_world,ierror)
     call mpi_bcast(t2_tmp2,ndim,mpi_real8,master, &
          mpi_comm_world,ierror)
     
     i1 = 0 
     do j = 1, ket
        do i = 1, bra 
           i1 = i1 + 1
           t2_cross_ccm_eqn%depot(t2_channel)%val(i,j) = 0.d0
           t2_cross_ccm_eqn%depot(t2_channel)%val(i,j) = t2_tmp2(i1)
        end do
     end do

     deallocate( t2_tmp1, t2_tmp2 ) 
  end do
  

  !call mpi_barrier(mpi_comm_world,ierror)
  !endwtime = MPI_WTIME()
  !if ( iam == 0 ) write(6,*) 'Total execution time for I4', endwtime - startwtime
  

  
  !call mpi_barrier(mpi_comm_world,ierror)
  !startwtime = MPI_WTIME()

  
  
  
  
  !
  ! setup dynamic denominators d2ab and d2ij
  !
  do i = 1, below_ef
     d2ij(i) = d1ik(i) 
  end do
  
  do a = below_ef+1, tot_orbs
     d2ab(a) = I6abc(a,a)
     I6abc(a,a) = 0.d0
  end do
  
  !
  ! Outer loop over all hhpp channels
  !
  do i1 = 1, channels%number_confs 
     ang_mom   = channels%config_J_ipar_tz(i1*3)
     p_parity  = channels%config_J_ipar_tz(i1*3-1)
     isospin_z = channels%config_J_ipar_tz(i1*3-2)
     

     t2_channel = locate_gmatchannel(3,isospin_z, p_parity, ang_mom)
     
     if ( t2_channel == 0 ) cycle
     if ( check_my_channel_hhpp(t2_channel) == 0 ) cycle
     
     
     local_ch_id = t2_channel
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
     

     t2_angmom_fact = 2.d0*j_hhpp(t2_channel) + 1.d0 
     t2_ccm_eqn%depot(t2_channel)%val(bra_min:bra_max,:) =t2_ccm_eqn%depot(t2_channel)%val(bra_min:bra_max,:) + &
          transpose( gmat_mtx_hhpp%depot(t2_channel)%val(:,bra_min:bra_max) ) 
     
     DO bra=bra_min, bra_max
        a = lookup_ab_configs(3,i1)%ival(1,bra) 
        b = lookup_ab_configs(3,i1)%ival(2,bra) 
        phase_ab = (-1.d0)**( (all_orbit%jj(a) + all_orbit%jj(b))/2 - j_hhpp(t2_channel) + 1) 
        
        DO ket=1, size(  lookup_ab_configs(1,i1)%ival, 2) 
           i = lookup_ab_configs(1,i1)%ival(1,ket) 
           j = lookup_ab_configs(1,i1)%ival(2,ket)  
           
           phase_ij = (-1.d0)**( (all_orbit%jj(i) + all_orbit%jj(j))/2 - j_hhpp(t2_channel) + 1) 

           !
           ! I6 contribution to t2-amps
           !
           sum1 = 0.d0
           do c = below_ef+1, tot_orbs
              bra2 = pp_hhpp(t2_channel)%ival(c,b) 
              if ( bra2 == 0 ) cycle

              t2 = t2_ccm%depot(t2_channel)%val(bra2,ket) 
              sum1 = sum1 + t2 * I6abc(a,c)  
           end do
           sum2 = 0.d0
           do c = below_ef+1, tot_orbs
              bra2 = pp_hhpp(t2_channel)%ival(c,a) 
              if ( bra2 == 0 ) cycle
              
              t2 = t2_ccm%depot(t2_channel)%val(bra2,ket) 

              sum2 = sum2 + t2 * I6abc(b,c)  

           end do

           
           t2_ccm_eqn%depot(t2_channel)%val(bra,ket)  = t2_ccm_eqn%depot(t2_channel)%val(bra,ket) + sum1 + phase_ab*sum2
           
           !
           ! I7 contribution to t2-amps
           !
           sum1 = 0.d0
           do k = 1, below_ef
              ket2 = hh_hhpp(t2_channel)%ival(k,j)
              if ( ket2 == 0 ) cycle
              t2 = t2_ccm%depot(t2_channel)%val(bra,ket2) 
              sum1 = sum1 + t2 * I7abcde(k,i)
           end do
           sum2 = 0.d0
           do k = 1, below_ef
              ket2 = hh_hhpp(t2_channel)%ival(k,i)
              if ( ket2 == 0 ) cycle
              
              t2 = t2_ccm%depot(t2_channel)%val(bra,ket2) 

              sum2 = sum2 + t2 * I7abcde(k,j)

           end do
           t2_ccm_eqn%depot(t2_channel)%val(bra,ket)  = t2_ccm_eqn%depot(t2_channel)%val(bra,ket) - (sum1 + phase_ij*sum2)
           

        end do
     end DO
  END do


  
  
  
  
  
  !call mpi_barrier(mpi_comm_world,ierror)
  !endwtime = MPI_WTIME()
  !if ( iam == 0 ) write(6,*) 'Total execution time for I8', endwtime - startwtime
  deallocate( I7abcde, I6abc )
  
  
  ! 
  ! I8 contribution to t2-amps 
  !
  do i1 = 1, channels%number_confs 
     ang_mom   = channels%config_J_ipar_tz(i1*3)
     p_parity  = channels%config_J_ipar_tz(i1*3-1)
     isospin_z = channels%config_J_ipar_tz(i1*3-2)
     

     t2_channel = locate_gmatchannel(3,isospin_z, p_parity, ang_mom)
     I8_channel = locate_gmatchannel(6,isospin_z, p_parity, ang_mom)
     
     if ( t2_channel*I8_channel == 0 ) cycle
     if ( check_my_channel(I8_channel) == 0 ) cycle
     
     local_ch_id = I8_channel
     
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
     
     
     allocate( amat(bra_min:bra_max, ket_min:ket_max) )
     allocate( bmat(size(  lookup_ab_configs(3,i1)%ival, 2),  size(  lookup_ab_configs(1,i1)%ival, 2) ) )
     
     bmat = 0.d0
     DO bra=1, size(  lookup_ab_configs(3,i1)%ival, 2) 
        a = lookup_ab_configs(3,i1)%ival(1,bra) 
        b = lookup_ab_configs(3,i1)%ival(2,bra) 
        phase_ab = iph( (all_orbit%jj(b)+all_orbit%jj(a))/2 + ang_mom + 1)
        
        DO ket=1, size(  lookup_ab_configs(1,i1)%ival, 2) 
           i = lookup_ab_configs(1,i1)%ival(1,ket) 
           j = lookup_ab_configs(1,i1)%ival(2,ket)  
           
           phase_ij = iph( (all_orbit%jj(i)+all_orbit%jj(j))/2 + ang_mom + 1)

           bmat(bra,ket) = bmat(bra,ket) + t1_ccm(a,i) * t1_ccm(b,j) & 
                + t1_ccm(b,i) * t1_ccm(a,j) * phase_ab & 
                + t1_ccm(a,j) * t1_ccm(b,i) * phase_ij & 
                + t1_ccm(b,j) * t1_ccm(a,i) * phase_ab * phase_ij 
           
        end DO
     end DO
     
     bmat = 0.5d0 * bmat + t2_ccm%depot(t2_channel)%val  
     
     amat(bra_min:bra_max, ket_min:ket_max) = gmat_mtx_pppp%depot(I8_channel)%val(bra_min:bra_max, ket_min:ket_max)
     
     dim1 = size( amat,1 ) 
     !dim1 = size( gmat_mtx_pppp%depot(I8_channel)%val,1 ) !size( amat, 1)
     dim2 = size( t2_ccm%depot(t2_channel)%val, 2)
     !dim3 =size( gmat_mtx_pppp%depot(I8_channel)%val,2 ) !size( amat, 2)
     dim3 = size( amat, 2) 
     
     
     ! Nov. 10
     call DGEMM ( 'n', 'n', dim1, dim2, dim3, 0.5d0, amat, & 
          dim1,bmat, dim3, & 
          1.d0, t2_ccm_eqn%depot(t2_channel)%val(bra_min:bra_max,:), dim1 )

     deallocate( amat, bmat )
  end do
  
  !
  ! I9(a) intermediate
  !
  number_channels = size(gmat_mtx_hhhh%depot)
  allocate( I9abc(1:number_channels) )

  do i=1, number_channels
     bra = size(gmat_mtx_hhhh%depot(i)%val,1)
     ket = size(gmat_mtx_hhhh%depot(i)%val,2)

     allocate( I9abc(i)%val(bra,ket) )
     !
     ! I9a =<kl||ij>
     !
     I9abc(i)%val = gmat_mtx_hhhh%depot(i)%val 

  end do
  

  !
  ! I9(b) intermediate
  !
  do i1 = 1, channels%number_confs 
     ang_mom   = channels%config_J_ipar_tz(i1*3)
     p_parity  = channels%config_J_ipar_tz(i1*3-1)
     isospin_z = channels%config_J_ipar_tz(i1*3-2)


     t2_channel = locate_gmatchannel(3,isospin_z, p_parity, ang_mom)
     I9_channel = locate_gmatchannel(1,isospin_z, p_parity, ang_mom)

     if ( t2_channel*I9_channel == 0 ) cycle

     
     dim1 = size( gmat_mtx_hhpp%depot(t2_channel)%val, 1)
     dim2 = size( t2_ccm%depot(t2_channel)%val,2 ) 
     dim3 = size( gmat_mtx_hhpp%depot(t2_channel)%val, 2 ) 
     
     
     call DGEMM ( 'n', 'n', dim1, dim2, dim3, 0.5d0, gmat_mtx_hhpp%depot(t2_channel)%val, & 
          dim1,t2_ccm%depot(t2_channel)%val, dim3, & 
          1.d0, I9abc(I9_channel)%val, dim1 )
     
  end do
  

  
  !
  ! I9(c) intermediate
  !
  do i1 = 1, channels%number_confs 
     ang_mom   = channels%config_J_ipar_tz(i1*3)
     p_parity  = channels%config_J_ipar_tz(i1*3-1)
     isospin_z = channels%config_J_ipar_tz(i1*3-2)


     I9_channel = locate_gmatchannel(1,isospin_z, p_parity, ang_mom)
     I3_channel = locate_gmatchannel(2,isospin_z, p_parity, ang_mom)
     
     if ( I9_channel*I3_channel == 0 ) cycle
     
     allocate( amat( size(i9abc(I9_channel)%val,1),size(i9abc(I9_channel)%val,2) ) )
     amat = 0.d0 

     !$omp parallel default(shared) private(bra3,k,l,bra,bra2,ket3,i,j,sum1,ket2,gmat)
     !$omp do schedule(static), reduction(+:amat)
     DO bra3=1, size(  lookup_ab_configs(1,i1)%ival, 2) 
        k = lookup_ab_configs(1,i1)%ival(1,bra3) 
        l = lookup_ab_configs(1,i1)%ival(2,bra3) 
        
        bra = hh_hhhh(I9_channel)%ival(k,l)
        if ( bra == 0 ) cycle
        bra2 = hh_hhhp(I3_channel)%ival(k,l)
        if ( bra2 == 0 ) cycle

        DO ket3=1, size(  lookup_ab_configs(1,i1)%ival, 2) 
           i = lookup_ab_configs(1,i1)%ival(1,ket3) 
           j = lookup_ab_configs(1,i1)%ival(2,ket3)  
           
           ket = hh_hhhh(I9_channel)%ival(i,j)
           if ( ket == 0 ) cycle


           sum1 = 0.d0
           do c = below_ef+1, tot_orbs
              ket2 = hp_hhhp(I3_channel)%ival(i,c)
              if ( ket2 == 0 ) cycle
              gmat = gmat_mtx_hhhp%depot(I3_channel)%val(bra2,ket2) 
                            
              sum1 = 0.5d0*( gmat + I3ab(i3_channel)%val(bra2, ket2) ) * t1_ccm(c,j) 
              amat(bra,ket) = amat(bra,ket) + sum1 
           end do
           
        end DO
     end DO
     !$omp end do                                                                                                                                      
     !$omp end parallel 
   

     I9abc(I9_channel)%val = I9abc(I9_channel)%val + amat 
     amat = 0.d0
     !$omp parallel default(shared) private(bra3,k,l,bra,bra2,ket3,i,j,sum2,ket2,gmat)
     !$omp do schedule(static), reduction(+:amat)
     DO bra3=1, size(  lookup_ab_configs(1,i1)%ival, 2) 
        k = lookup_ab_configs(1,i1)%ival(1,bra3) 
        l = lookup_ab_configs(1,i1)%ival(2,bra3) 
        
        bra = hh_hhhh(I9_channel)%ival(k,l)
        if ( bra == 0 ) cycle
        bra2 = hh_hhhp(I3_channel)%ival(k,l)
        if ( bra2 == 0 ) cycle

        DO ket3=1, size(  lookup_ab_configs(1,i1)%ival, 2) 
           i = lookup_ab_configs(1,i1)%ival(1,ket3) 
           j = lookup_ab_configs(1,i1)%ival(2,ket3)  
           ! 
           ! (ij)
           !          
           ket = hh_hhhh(I9_channel)%ival(i,j)
           if ( ket == 0 ) cycle
           phase_ij = (-1.d0)**( (all_orbit%jj(i) + all_orbit%jj(j))/2 - j_hhhh(I9_channel) + 1) 

           sum2 = 0.d0
           do c = below_ef+1, tot_orbs
              ket2 = hp_hhhp(I3_channel)%ival(j,c)
              if ( ket2 == 0 ) cycle
              gmat = gmat_mtx_hhhp%depot(I3_channel)%val(bra2,ket2) 
              
              sum2 =  0.5d0* phase_ij*(gmat+ I3ab(i3_channel)%val(bra2, ket2)) * t1_ccm(c,i) 
              amat(bra,ket) = amat(bra,ket) + sum2
           end do
           
        end do
     end do
     !$omp end do                                                                                                                                      
     !$omp end parallel 

     I9abc(I9_channel)%val = I9abc(I9_channel)%val + amat 
     deallocate( amat ) 
  end do
  
  
  do i = 1, size( I3ab(:) )
     deallocate( I3ab(i)%val )
  end do
  deallocate( I3ab )

  ! 
  ! I9(abc) contribution to t2-amps
  !
  do i1 = 1, channels%number_confs 
     ang_mom   = channels%config_J_ipar_tz(i1*3)
     p_parity  = channels%config_J_ipar_tz(i1*3-1)
     isospin_z = channels%config_J_ipar_tz(i1*3-2)


     t2_channel = locate_gmatchannel(3,isospin_z, p_parity, ang_mom)
     I9_channel = locate_gmatchannel(1,isospin_z, p_parity, ang_mom)
     
     
     if ( t2_channel*I9_channel == 0 ) cycle
     if ( check_my_channel_hhpp( t2_channel) == 0 ) cycle

     local_ch_id = t2_channel
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
     !dim1 = size( t2_ccm%depot(t2_channel)%val,1) 
     dim2 = size( I9abc(I9_channel)%val,2 )
     dim3 = size( t2_ccm%depot(t2_channel)%val,2) 

     ! Nov. 10
     call DGEMM ( 'n', 'n', dim1, dim2, dim3, 0.5d0, t2_ccm%depot(t2_channel)%val(bra_min:bra_max,:), & 
          dim1,I9abc(I9_channel)%val, dim3, & 
          1.d0, t2_ccm_eqn%depot(t2_channel)%val(bra_min:bra_max,:), dim1 )
     
     
  end do

  !call mpi_barrier(mpi_comm_world,ierror)
  !endwtime = MPI_WTIME()
  !if ( iam == 0 ) write(6,*) 'Total execution time for I9', endwtime - startwtime

  

  !call mpi_barrier(mpi_comm_world,ierror)
  !startwtime = MPI_WTIME()
  !
  ! I10 (af) contribution to t2-amps
  !
  do i1 = 1, channels%number_confs 
     ang_mom   = channels%config_J_ipar_tz(i1*3)
     p_parity  = channels%config_J_ipar_tz(i1*3-1)
     isospin_z = channels%config_J_ipar_tz(i1*3-2)
     

     t2_channel = locate_gmatchannel(3,isospin_z, p_parity, ang_mom)
     if ( t2_channel == 0 ) cycle
     if ( check_my_channel_hhpp(t2_channel) == 0 ) cycle

     I10_channel = locate_gmatchannel(5,isospin_z, p_parity, ang_mom)
     if ( I10_channel == 0 ) cycle

     local_ch_id = t2_channel
     
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
     
     

     t2_angmom_fact = 2.d0*j_hhpp(t2_channel) + 1.d0 
     
     DO bra=bra_min, bra_max
        a = lookup_ab_configs(3,i1)%ival(1,bra) 
        b = lookup_ab_configs(3,i1)%ival(2,bra) 
        
        ket2 = pp_hppp(I10_channel)%ival(a,b)
        phase_ab = 1.d0
        if ( a > b ) then
           phase_ab = iph( (all_orbit%jj(a)+all_orbit%jj(b))/2 + ang_mom + 1)
           ket2 = pp_hppp(I10_channel)%ival(b,a)
        end if
        if ( ket2 == 0 ) cycle

        DO ket=1, size(  lookup_ab_configs(1,i1)%ival, 2) 
           i = lookup_ab_configs(1,i1)%ival(1,ket) 
           j = lookup_ab_configs(1,i1)%ival(2,ket)  
           
           phase_ij = (-1.d0)**( (all_orbit%jj(i) + all_orbit%jj(j))/2 - j_hhpp(t2_channel) + 1) 

           !
           ! I10a contribution to t2-amps
           !
           sum1 = 0.d0
           do c = below_ef+1, tot_orbs

              bra2 = hp_hppp(I10_channel)%ival(i,c)
              if ( bra2 == 0 ) cycle

              gmat = phase_ab * gmat_mtx_hppp%depot(I10_channel)%val(bra2,ket2) ! I10af(I10_channel)%val(bra2,ket2) 
              sum1 = sum1 + gmat * t1_ccm(c,j) 
           end do

           t2_ccm_eqn%depot(t2_channel)%val(bra,ket) = t2_ccm_eqn%depot(t2_channel)%val(bra,ket) + sum1 

           !
           ! P(ij)
           !
           ket3 = hh_hhpp(t2_channel)%ival(j,i)
           t2_ccm_eqn%depot(t2_channel)%val(bra,ket3) = t2_ccm_eqn%depot(t2_channel)%val(bra,ket3) +phase_ij* sum1 

        end do
     end DO
  end do
  
  !  DEALLOCATE( I10af )
  !call mpi_barrier(mpi_comm_world,ierror)
  !endwtime = MPI_WTIME()
  !if ( iam == 0 ) write(6,*) 'Total execution time for I10', endwtime - startwtime

  !call mpi_barrier(mpi_comm_world,ierror)
  !startwtime = MPI_WTIME()
    
  
  !
  ! I11(a) intermediate
  !
  number_channels = size(gmat_mtx_hhhp%depot)
  allocate( I11abdef(1:number_channels) )
  
  do i=1, size( I11abdef ) 
     bra = size(gmat_mtx_hhhp%depot(i)%val,1)
     ket = size(gmat_mtx_hhhp%depot(i)%val,2)

     allocate( I11abdef(i)%val(bra,ket) )
     !
     ! I11a =<ak||ij>
     !
     I11abdef(i)%val = 0.d0 !gmat_mtx_hhhp%depot(i)%val 

  end do
  
  
  !
  ! allocate and setup matrix for I11(d) intermediate
  !
  number_channels = size(gmat_mtx_hhpp%depot)
  allocate( temp_mat(1:number_channels) )
  do channel   = 1, channels%number_confs 
     ang_mom   = channels%config_J_ipar_tz(channel*3)
     p_parity  = channels%config_J_ipar_tz(channel*3-1)
     isospin_z = channels%config_J_ipar_tz(channel*3-2)
     
     channel1 = locate_gmatchannel(3,isospin_z, p_parity, ang_mom)
     if ( channel1 == 0 ) cycle
     
     bra2 = size(gmat_mtx_hhpp%depot(channel1)%val,2)
     ket2 = size(gmat_mtx_hhpp%depot(channel1)%val,1)
     
     allocate( temp_mat(channel1)%val(bra2,ket2) )
     temp_mat(channel1)%val = 0.d0
     

     
     DO bra=1, size(  lookup_ab_configs(3,channel)%ival, 2) 
        c = lookup_ab_configs(3,channel)%ival(1,bra) 
        d = lookup_ab_configs(3,channel)%ival(2,bra) 
        
        DO ket=1, size(  lookup_ab_configs(1,channel)%ival, 2) 
           i = lookup_ab_configs(1,channel)%ival(1,ket) 
           j = lookup_ab_configs(1,channel)%ival(2,ket)  

           
           phase_ij = (-1.d0)**( (all_orbit%jj(i) + all_orbit%jj(j))/2 - ang_mom + 1) 

           temp_mat(channel1)%val(bra,ket) = t1_ccm(c,i)*t1_ccm(d,j) + phase_ij * t1_ccm(c,j)*t1_ccm(d,i)
            
        end do
     end DO
     
     temp_mat(channel1)%val = temp_mat(channel1)%val + t2_ccm%depot(channel1)%val

  end do


  
  
  local_ch_id = 0 
  do i1 = 1, channels%number_confs 
     ang_mom   = channels%config_J_ipar_tz(i1*3)
     p_parity  = channels%config_J_ipar_tz(i1*3-1)
     isospin_z = channels%config_J_ipar_tz(i1*3-2)
     

     I11_channel = locate_gmatchannel(2,isospin_z, p_parity, ang_mom)
     t2_channel = locate_gmatchannel(3,isospin_z, p_parity, ang_mom)
     I8_channel = locate_gmatchannel(5,isospin_z, p_parity, ang_mom)

     if ( t2_channel*I11_channel*I8_channel == 0 ) cycle
     local_ch_id = I8_channel !local_ch_id + 1 
     !if ( check_my_channel_hppp(I8_channel) == 0 ) cycle
     if ( check_my_channel_hppp(local_ch_id) == 0 ) cycle
     
     !local_ch_id = I8_channel
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
     
   
     dim1 = size( temp_mat(t2_channel)%val,2 )
     dim2 = size( gmat_mtx_hppp%depot(I8_channel)%val,1 )
     dim3 = bra_max-bra_min+1 
     
     allocate( amat(dim2, size(lookup_ab_configs(3,i1)%ival, 2) )) 
     amat= 0.d0

     DO i=1, size(  lookup_ab_configs(3,i1)%ival, 2) 
        a = lookup_ab_configs(3,i1)%ival(1,i) 
        b = lookup_ab_configs(3,i1)%ival(2,i) 
        
        
        ket = pp_hppp(I8_channel)%ival(a,b)
        phase_ab = 1.d0
        if ( a > b ) then
           phase_ab = iph( (all_orbit%jj(a) + all_orbit%jj(b))/2 + ang_mom +1 )
           ket = pp_hppp(I8_channel)%ival(b,a)
        end if
        if ( ket == 0 ) cycle
        amat(:,i) = gmat_mtx_hppp%depot(I8_channel)%val(:,ket) * phase_ab 
     end DO
     
     !dim3 = size( t2_ccm%depot(t2_channel)%val,1) 

     ld1 = bra_max-bra_min+1 ! 
     !ld1 = size( temp_mat(t2_channel)%val,1 ) 
     ld2 = size( gmat_mtx_hppp%depot(I8_channel)%val,1 )
     ld3 = size( I11abdef(I11_channel)%val,1 )
     
     ! Nov. 10
     call DGEMM ( 't', 't', dim1, dim2, dim3, 0.5d0, temp_mat(t2_channel)%val(bra_min:bra_max,:), & 
          ld1, amat(:,bra_min:bra_max), ld2, & 
          1.d0, I11abdef(I11_channel)%val, ld3 )
     
     
     deallocate( amat )
  
  end do
  
  do i = 1, size(temp_mat(:) )
     deallocate( temp_mat(i)%val )
  end do
  deallocate( temp_mat )

  !write(6,*) 'iam', iam, 'count for I11', count
            
  number_channels = size( I11abdef )
  do t2_channel = 1, number_channels 
     bra = size(I11abdef(t2_channel)%val,1)
     ket = size(I11abdef(t2_channel)%val,2)
     
     ndim = bra * ket
     allocate( t2_tmp1(ndim) )
     allocate( t2_tmp2(ndim) )
     t2_tmp1 = 0.d0
     t2_tmp2 = 0.d0

     i1 = 0 
     do i =1, bra
        do j = 1, ket
           i1 = i1 + 1
           !t2_tmp1(i1) =  temp_mat(t2_channel)%val(i,j) 
           ! t2_tmp1(i1) =  temp_mat(t2_channel)%val(i,j) + t2_ccm_eqn%depot(t2_channel)%val(i,j)
           t2_tmp1(i1) = I11abdef(t2_channel)%val(i,j)
        end do
     end do
     
     call mpi_reduce(t2_tmp1,t2_tmp2,ndim,mpi_real8,mpi_sum, &
          master,mpi_comm_world,ierror)
     t2_tmp1 = 0.d0 
     t2_tmp1 = t2_tmp2 
     call mpi_bcast(t2_tmp1,ndim,mpi_real8,master, &
          mpi_comm_world,ierror)
     
     i1 = 0 
     do i =1, bra
        do j = 1, ket
           i1 = i1 + 1
           I11abdef(t2_channel)%val(i,j) = 0.d0
           !t2_ccm_eqn%depot(t2_channel)%val(i,j) = t2_ccm_eqn%depot(t2_channel)%val(i,j) + t2_tmp1(i1)
           I11abdef(t2_channel)%val(i,j) = t2_tmp1(i1)
        end do
     end do
     
     deallocate( t2_tmp1, t2_tmp2 ) 
  end do

  !call mpi_barrier(mpi_comm_world,ierror)
  

  !
  ! I11(b) is rewritten and added to I4abc intermediate
  !

  !
  ! I11(d) intermediate
  !
  do i1 = 1, channels%number_confs 
     ang_mom   = channels%config_J_ipar_tz(i1*3)
     p_parity  = channels%config_J_ipar_tz(i1*3-1)
     isospin_z = channels%config_J_ipar_tz(i1*3-2)

     I11_channel = locate_gmatchannel(2,isospin_z, p_parity, ang_mom)
     I8_channel = locate_gmatchannel(4,isospin_z, p_parity, ang_mom)
     if ( I11_channel * I8_channel == 0 ) cycle

     allocate( amat(size(I11abdef(I11_channel)%val,1), size(I11abdef(I11_channel)%val,2) ) )
     amat = 0.d0 
     !$omp parallel default(shared) private(bra,i,j,phase_ij,ket,k,b,bra2,sum1,c,ket2,phase,gmat)
     !$omp do schedule(static), reduction(+:amat)
     DO bra=1, size(  lookup_ab_configs(1,i1)%ival, 2) 
        i = lookup_ab_configs(1,i1)%ival(1,bra) 
        j = lookup_ab_configs(1,i1)%ival(2,bra) 

        phase_ij = (-1.d0)**( (all_orbit%jj(i) + all_orbit%jj(j))/2 - ang_mom + 1) 

        DO ket=1, size(  lookup_ab_configs(2,i1)%ival, 2) 
           k = lookup_ab_configs(2,i1)%ival(1,ket) 
           b = lookup_ab_configs(2,i1)%ival(2,ket) 
           
        
           bra2 = hp_hphp(I8_channel)%ival(k,b)
           if ( bra2 == 0 ) cycle
        
           sum1 = 0.d0
           do c = below_ef+1, tot_orbs
                    
              ket2 = hp_hphp(I8_channel)%ival(j,c)
              if ( ket2 == 0 ) cycle
                    
              phase = (-1.d0)**( (all_orbit%jj(j) + all_orbit%jj(c))/2 - ang_mom + 1) 
              
              gmat = gmat_mtx_hphp%depot(I8_channel)%val(bra2,ket2)
              sum1 = gmat * t1_ccm(c,i) * phase
              amat( bra,ket) = amat(bra,ket) + sum1 
           end do
        end DO
     end DO
     !$omp end do                                                                                                                                      
     !$omp end parallel 
   

     I11abdef(I11_channel)%val = I11abdef(I11_channel)%val + amat 
     
     amat = 0.d0 
     !$omp parallel default(shared) private(bra,i,j,phase_ij,ket,k,b,bra2,sum1,c,ket2,phase,gmat)
     !$omp do schedule(static), reduction(+:amat)
     DO bra=1, size(  lookup_ab_configs(1,i1)%ival, 2) 
        i = lookup_ab_configs(1,i1)%ival(1,bra) 
        j = lookup_ab_configs(1,i1)%ival(2,bra) 

        phase_ij = (-1.d0)**( (all_orbit%jj(i) + all_orbit%jj(j))/2 - ang_mom + 1) 

        DO ket=1, size(  lookup_ab_configs(2,i1)%ival, 2) 
           k = lookup_ab_configs(2,i1)%ival(1,ket) 
           b = lookup_ab_configs(2,i1)%ival(2,ket) 
           
        
           bra2 = hp_hphp(I8_channel)%ival(k,b)
           if ( bra2 == 0 ) cycle
        
           sum1 = 0.d0
           do c = below_ef+1, tot_orbs
                    
              ket2 = hp_hphp(I8_channel)%ival(i,c)
              if ( ket2 == 0 ) cycle
                    
              phase = (-1.d0)**( (all_orbit%jj(i) + all_orbit%jj(c))/2 - ang_mom + 1) 
              
              gmat = gmat_mtx_hphp%depot(I8_channel)%val(bra2,ket2)
              sum1 = gmat * t1_ccm(c,j) * phase
              amat( bra,ket) = amat(bra,ket) + sum1 * phase_ij 
           end do
        end DO
     end DO
     !$omp end do                                                                                                                                      
     !$omp end parallel 
     
     I11abdef(I11_channel)%val = I11abdef(I11_channel)%val + amat 
     deallocate( amat ) 
     

  end do

  !
  ! I11(e) intermediate
  !
  do i1 = 1, channels%number_confs 
     ang_mom   = channels%config_J_ipar_tz(i1*3)
     p_parity  = channels%config_J_ipar_tz(i1*3-1)
     isospin_z = channels%config_J_ipar_tz(i1*3-2)


     I11_channel = locate_gmatchannel(2,isospin_z, p_parity, ang_mom)
     t2_channel = locate_gmatchannel(3,isospin_z, p_parity, ang_mom)
     if ( I11_channel*t2_channel == 0 ) cycle
     
     allocate( amat(size(I11abdef(I11_channel)%val,1), size(I11abdef(I11_channel)%val,2) ) )
     amat = 0.d0 
     !$omp parallel default(shared) private(bra,i,j,ket2,ket,k,b,sum1,c,bra2,t2)
     !$omp do schedule(static), reduction(+:amat)
     DO bra=1, size(  lookup_ab_configs(1,i1)%ival, 2) 
        i = lookup_ab_configs(1,i1)%ival(1,bra) 
        j = lookup_ab_configs(1,i1)%ival(2,bra) 
        
        ket2 = hh_hhpp(t2_channel)%ival(i,j)
        if ( ket2 == 0 ) cycle
        DO ket=1, size(  lookup_ab_configs(2,i1)%ival, 2) 
           k = lookup_ab_configs(2,i1)%ival(1,ket) 
           b = lookup_ab_configs(2,i1)%ival(2,ket) 
    
           sum1 = 0.d0
           do c = below_ef+1, tot_orbs

              bra2 = pp_hhpp(t2_channel)%ival(c,b)  
              if ( bra2 == 0 ) cycle

              t2 = t2_ccm%depot(t2_channel)%val(bra2, ket2)
              sum1 = t2 * I5ab(k,c) 
              amat( bra,ket) = amat(bra,ket) + sum1 
           end do
        end do
     end DO
     !$omp end do                                                                                                                                      
     !$omp end parallel 
     
     I11abdef(I11_channel)%val = I11abdef(I11_channel)%val + amat 
     deallocate( amat ) 
  end do
  deallocate( I5ab )

  !
  ! I11(f) intermediate
  !
  do i1 = 1, channels%number_confs 
     ang_mom   = channels%config_J_ipar_tz(i1*3)
     p_parity  = channels%config_J_ipar_tz(i1*3-1)
     isospin_z = channels%config_J_ipar_tz(i1*3-2)

     
     I11_channel = locate_gmatchannel(2,isospin_z, p_parity, ang_mom)
     I9_channel = locate_gmatchannel(1,isospin_z, p_parity, ang_mom)

     if ( I11_channel*I9_channel == 0 ) cycle
     
     allocate( amat(size(I11abdef(I11_channel)%val,1), size(I11abdef(I11_channel)%val,2) ) )
     amat = 0.d0 
     !$omp parallel default(shared) private(bra,i,j,ket2,ket,k,b,sum1,l,bra2)
     !$omp do schedule(static), reduction(+:amat)
     DO bra=1, size(  lookup_ab_configs(1,i1)%ival, 2) 
        i = lookup_ab_configs(1,i1)%ival(1,bra) 
        j = lookup_ab_configs(1,i1)%ival(2,bra) 
        
        ket2 = hh_hhhh(I9_channel)%ival(i,j)
        if ( ket2 == 0 ) cycle

        DO ket=1, size(  lookup_ab_configs(2,i1)%ival, 2) 
           k = lookup_ab_configs(2,i1)%ival(1,ket) 
           b = lookup_ab_configs(2,i1)%ival(2,ket) 
    
           sum1 = 0.d0
           do l = 1, below_ef

              bra2 = hh_hhhh(I9_channel)%ival(k,l)
              if ( bra2 == 0 ) cycle


              sum1 =  -0.5d0*I9abc(I9_channel)%val(bra2,ket2)*t1_ccm(b,l)
              amat( bra,ket) = amat(bra,ket) + sum1 
              
           end do

        end DO
     end DO
     !$omp end do                                                                                                                                      
     !$omp end parallel 
     
     I11abdef(I11_channel)%val = I11abdef(I11_channel)%val + amat 
     deallocate( amat ) 

  end do

  do i = 1, size( I9abc(:) )
     deallocate( I9abc(i)%val )
  end do
  deallocate( I9abc )

  !write(6,*) 'time10'
  !
  ! I11 (a) contribution to t2-amps
  !
  do i1 = 1, channels%number_confs 
     ang_mom   = channels%config_J_ipar_tz(i1*3)
     p_parity  = channels%config_J_ipar_tz(i1*3-1)
     isospin_z = channels%config_J_ipar_tz(i1*3-2)
     
     
     t2_channel = locate_gmatchannel(3,isospin_z, p_parity, ang_mom)
     
     if ( t2_channel == 0 ) cycle
     I11_channel = locate_gmatchannel(2,isospin_z, p_parity, ang_mom)
     if ( I11_channel == 0 ) cycle
     if ( check_my_channel_hhpp( t2_channel ) == 0 ) cycle

     local_ch_id = t2_channel
     
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
     
     
     
     DO bra=bra_min, bra_max
        a = lookup_ab_configs(3,i1)%ival(1,bra) 
        b = lookup_ab_configs(3,i1)%ival(2,bra) 
        phase_ab = (-1.d0)**( (all_orbit%jj(a) + all_orbit%jj(b))/2 - j_hhpp(t2_channel) + 1) 
              
        DO ket=1, size(  lookup_ab_configs(1,i1)%ival, 2) 
           i = lookup_ab_configs(1,i1)%ival(1,ket) 
           j = lookup_ab_configs(1,i1)%ival(2,ket) 
    
           bra2 = hh_hhhp(I11_channel)%ival(i,j)
           if ( bra2 == 0 ) cycle
           

           !
           ! I11abdef contribution to t2-amps
           !
           sum1 = 0.d0
           do k = 1, below_ef

              ket2 = hp_hhhp(I11_channel)%ival(k,b)
              if ( ket2 == 0 ) cycle
              

              phase = iph( (all_orbit%jj(b)+all_orbit%jj(k))/2 - ang_mom )
              gmat = I11abdef(I11_channel)%val(bra2,ket2) + gmat_mtx_hhhp%depot(I11_channel)%val(bra2,ket2)  
              sum1 = sum1 + gmat * t1_ccm(a,k) 
           end do
           
           
           !
           ! P(ab)
           !
           sum2 = 0.d0
           do k = 1, below_ef

              ket2 = hp_hhhp(I11_channel)%ival(k,a)
              if ( ket2 == 0 ) cycle
              

              phase = iph( (all_orbit%jj(a)+all_orbit%jj(k))/2 - ang_mom )
              gmat = I11abdef(I11_channel)%val(bra2,ket2) + gmat_mtx_hhhp%depot(I11_channel)%val(bra2,ket2)  
              sum2 = sum2 + gmat * t1_ccm(b,k) 
           end do
           
           t2_ccm_eqn%depot(t2_channel)%val(bra,ket) = t2_ccm_eqn%depot(t2_channel)%val(bra,ket) - & 
                sum1 - sum2 * phase_ab 
           
        end do
     end DO
    
  end do

  do i = 1, size( I11abdef(:) ) 
     DEALLOCATE( I11abdef(i)%val )
  end do
  DEALLOCATE( I11abdef ) 

  !call mpi_barrier(mpi_comm_world,ierror)
  !endwtime = MPI_WTIME()
  !if ( iam == 0 ) write(6,*) 'Total execution time for I11', endwtime - startwtime



  !write(6,*) 'time12'
  !
  ! go from cross-coupled matrix elements to pp coupled mtx-elements 
  ! 
  ! Coupling direction: ingoing -> outgoing 
  ! <ij||ab> (ai->J) and (bj->J) if one changes coupling direction 
  ! to outgoing -> ingoing one gets a phase (-1)**(ji+jj+ja+jb - 2J ) 
  ! 
  bra_conf = 3
  ket_conf = 1
  do channel   = 1, channels%number_confs 
     ang_mom   = channels%config_J_ipar_tz(channel*3)
     p_parity  = channels%config_J_ipar_tz(channel*3-1)
     isospin_z = channels%config_J_ipar_tz(channel*3-2)

     t2_channel = locate_gmatchannel(3,isospin_z,p_parity,ang_mom)
     if ( t2_channel == 0 ) cycle
     if ( check_my_channel_hhpp( t2_channel ) == 0 ) cycle
     
     
     local_ch_id = t2_channel
     
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
     
     DO bra=bra_min, bra_max
        a = lookup_ab_configs(3,channel)%ival(1,bra) 
        b = lookup_ab_configs(3,channel)%ival(2,bra) 
        
        ja=all_orbit%jj(a)
        jb=all_orbit%jj(b)
                      
        DO ket=1, size(  lookup_ab_configs(1,channel)%ival, 2) 
           i = lookup_ab_configs(1,channel)%ival(1,ket) 
           j = lookup_ab_configs(1,channel)%ival(2,ket) 

           tz = (all_orbit%itzp(a) - all_orbit%itzp(i))/2
           ipar = (1 - (-1)**(all_orbit%ll(a)-all_orbit%ll(i) ) )/2
           
           ji=all_orbit%jj(i)
           jj=all_orbit%jj(j)
           
           j_min=max( ABS((all_orbit%jj(i)-all_orbit%jj(a))/2), ABS((all_orbit%jj(j)-all_orbit%jj(b))/2) )
           j_max=min( (all_orbit%jj(i)+all_orbit%jj(a))/2, (all_orbit%jj(j)+all_orbit%jj(b))/2 ) 

           phase_ij = (-1.d0)**( (all_orbit%jj(i) + all_orbit%jj(j))/2 - ang_mom + 1) 
           phase_ab = (-1.d0)**( (all_orbit%jj(a) + all_orbit%jj(b))/2 - ang_mom + 1) 
           
           
           !
           ! Additional phase since we couple in reverse directon (ia->J) and (jb->J)
           !
           
           fnorm= iph((jj+ja)/2+ang_mom) 
           
           
           IF(j_min > j_max) CYCLE
           cross=0.
           DO jt=j_min,j_max
              
              IF(triag(all_orbit%jj(a),all_orbit%jj(i),2*jt)) CYCLE
              IF(triag(all_orbit%jj(b),all_orbit%jj(j),2*jt)) CYCLE
              
              if ( locate_crosscoupled_channel(1, tz, ipar, jt) == 0 ) cycle
              ph_channel = locate_crosscoupled_channel(1, tz, ipar, jt ) 
              
              bra2 = ph_phhp_cross(ph_channel)%ival(a,i)
              !if ( bra2 == 0 ) cycle
              
              ket2 = hp_phhp_cross(ph_channel)%ival(j,b)
              !if ( ket2 == 0 ) cycle
              
              ang_mom_factor=sjs(ji,jj,2*ang_mom, jb,ja,2*jt)* &
                   (2.*jt+1.)*fnorm*iph(jt)
              
              
              count = count + 1
              !t2_old = t2_cross_ccm%depot(ph_channel)%val(bra2,ket2)
              t2_old = t2_cross_ccm_eqn%depot(ph_channel)%val(bra2,ket2)
              cross = cross + ang_mom_factor*t2_old
              
           end DO

           t2_ccm_eqn%depot(t2_channel)%val(bra,ket)  = t2_ccm_eqn%depot(t2_channel)%val(bra,ket) + cross
           
           bra3 = pp_hhpp(t2_channel)%ival(b,a) 
           t2_ccm_eqn%depot(t2_channel)%val(bra3,ket)  = t2_ccm_eqn%depot(t2_channel)%val(bra3,ket) + &
                phase_ab * cross
           
           ket3 = hh_hhpp(t2_channel)%ival(j,i)
           t2_ccm_eqn%depot(t2_channel)%val(bra,ket3)  = t2_ccm_eqn%depot(t2_channel)%val(bra,ket3) + &
                phase_ij * cross
           
           bra3 = pp_hhpp(t2_channel)%ival(b,a)
           ket3 = hh_hhpp(t2_channel)%ival(j,i)
           t2_ccm_eqn%depot(t2_channel)%val(bra3,ket3)  = t2_ccm_eqn%depot(t2_channel)%val(bra3,ket3) + &
                phase_ab * phase_ij * cross
                      
        end do
     end do

  end do
  

  
  number_channels = size( t2_ccm%depot )
  do t2_channel = 1, number_channels 
     bra = size(t2_ccm%depot(t2_channel)%val,1)
     ket = size(t2_ccm%depot(t2_channel)%val,2)
     
     ndim = bra * ket
     allocate( t2_tmp1(ndim) )
     allocate( t2_tmp2(ndim) )
     t2_tmp1 = 0.d0
     t2_tmp2 = 0.d0

     i1 = 0 
     do i =1, bra
        do j = 1, ket
           i1 = i1 + 1
           !t2_tmp1(i1) =  temp_mat(t2_channel)%val(i,j) 
           ! t2_tmp1(i1) =  temp_mat(t2_channel)%val(i,j) + t2_ccm_eqn%depot(t2_channel)%val(i,j)
           t2_tmp1(i1) = t2_ccm_eqn%depot(t2_channel)%val(i,j)
        end do
     end do
     
     call mpi_reduce(t2_tmp1,t2_tmp2,ndim,mpi_real8,mpi_sum, &
          master,mpi_comm_world,ierror)
     t2_tmp1 = 0.d0 
     t2_tmp1 = t2_tmp2 
     call mpi_bcast(t2_tmp1,ndim,mpi_real8,master, &
          mpi_comm_world,ierror)
     
     i1 = 0 
     do i =1, bra
        do j = 1, ket
           i1 = i1 + 1
           !t2_ccm_eqn%depot(t2_channel)%val(i,j) = t2_ccm_eqn%depot(t2_channel)%val(i,j) + t2_tmp1(i1)
           t2_ccm_eqn%depot(t2_channel)%val(i,j) = t2_tmp1(i1)
        end do
     end do
     
     deallocate( t2_tmp1, t2_tmp2 ) 
  end do

  


end  SUBROUTINE t2_intermediate
