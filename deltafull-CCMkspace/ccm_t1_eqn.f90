!
! This routine has been checked with Davids, and is correct.
!
SUBROUTINE t1_intermediate
  
  USE single_particle_orbits
  USE configurations
  USE constants
  use t2_storage
  use one_body_operators
  
  IMPLICIT NONE
  
  !real*8  ::  startwtime , endwtime
  !double precision :: ddot
  integer :: channel, nx,ny,nz,sz,tz,bra,ket, l
  integer :: k,a,i,c,d, channel2, bra2,ket2 
  real*8 :: vmom_yukawa, sum1, gmat, t2  
  double precision :: ddot 

  t1_ccm_eqn = 0.d0
  

  !
  ! setup intermediates 
  !
  !startwtime = MPI_WTIME()
  call setup_I3_intermediate
  call setup_I5I6I7_intermediate
  
  
  
  
  !
  ! I2a contribution to t1 ( ok )
  ! 
  do channel2   = 1, channels%number_hhpp_confs 
     nx = channels%hhpp_quantum_numbers(channel2*5)
     ny = channels%hhpp_quantum_numbers(channel2*5-1)
     nz = channels%hhpp_quantum_numbers(channel2*5-2)
     sz = channels%hhpp_quantum_numbers(channel2*5-3)
     tz = channels%hhpp_quantum_numbers(channel2*5-4)
     
     channel = LOCATE_CHANNEL(5, tz, sz, nx, ny, nz)
     if ( channel == 0 ) cycle 
     
     do bra = 1, size(  lookup_hppp_configs(1,channel)%ival, 2) 
        k = lookup_hppp_configs(1,channel)%ival(1,bra) 
        a = lookup_hppp_configs(1,channel)%ival(2,bra) 
        do ket = 1, size(  lookup_hppp_configs(2,channel)%ival, 2) 
           c = lookup_hppp_configs(2,channel)%ival(1,ket)
           d = lookup_hppp_configs(2,channel)%ival(2,ket) 
           bra2 = pp_hhpp%ival(c,d)
           if ( bra2 == 0 ) cycle
           do i = below_ef+1, tot_orbs 
              
              ket2 = hh_hhpp%ival(k,i)
              t2 = t2_ccm(channel)%val(bra2,ket2)
              gmat = vmom_yukawa(k,a,c,d) 
              sum1 = 0.5d0* gmat * t2 
              t1_ccm_eqn(a,i) = t1_ccm_eqn(a,i) + sum1 
              
           end do
        end do
     end do
  end do
  
  
  
  
  !
  ! I3ab contribution to t1(ok )
  ! 
  do channel   = 1, channels%number_hhpp_confs 
     nx = channels%hhpp_quantum_numbers(channel*5)
     ny = channels%hhpp_quantum_numbers(channel*5-1)
     nz = channels%hhpp_quantum_numbers(channel*5-2)
     sz = channels%hhpp_quantum_numbers(channel*5-3)
     tz = channels%hhpp_quantum_numbers(channel*5-4)
     
     channel2 = locate_channel(2,  tz, sz, nx, ny, nz)
     if ( channel2 == 0 ) cycle 

     do bra = 1, size(  lookup_hhpp_configs(1,channel)%ival, 2) 
        k = lookup_hhpp_configs(1,channel)%ival(1,bra) 
        l = lookup_hhpp_configs(1,channel)%ival(2,bra) 
        
        bra2  = hh_hhhp%ival(k,l)
        if ( bra2 == 0 ) cycle 
        do ket = 1, size(  lookup_hhpp_configs(2,channel)%ival, 2) 
           a = lookup_hhpp_configs(2,channel)%ival(1,ket)
           c = lookup_hhpp_configs(2,channel)%ival(2,ket) 
           
           
           do i = 1, below_ef
              ket2 = hp_hhhp%ival(i,c)
              t2 = t2_ccm(channel)%val(ket,bra)
              gmat = I3ab(channel2)%val(bra2, ket2) 
              sum1 = -0.5d0* gmat * t2 
              t1_ccm_eqn(a,i) = t1_ccm_eqn(a,i) + sum1 
              
           end do
        end do
     end do
  end do
  
  
  !deallocate( I3ab )
  !
  ! I5a contribution to t1 (ok)
  !
  do channel   = 1, channels%number_hhpp_confs 
     nx = channels%hhpp_quantum_numbers(channel*5)
     ny = channels%hhpp_quantum_numbers(channel*5-1)
     nz = channels%hhpp_quantum_numbers(channel*5-2)
     sz = channels%hhpp_quantum_numbers(channel*5-3)
     tz = channels%hhpp_quantum_numbers(channel*5-4)
  
     do bra = 1, size(  lookup_hhpp_configs(1,channel)%ival, 2) 
        k = lookup_hhpp_configs(1,channel)%ival(1,bra) 
        i = lookup_hhpp_configs(1,channel)%ival(2,bra) 
        
        do ket = 1, size(  lookup_hhpp_configs(2,channel)%ival, 2) 
           c = lookup_hhpp_configs(2,channel)%ival(1,ket)
           a = lookup_hhpp_configs(2,channel)%ival(2,ket) 
           
           t2 = t2_ccm(channel)%val(ket,bra) 
           sum1 =  t2 * I5ab(k,c) 
           t1_ccm_eqn(a,i) = t1_ccm_eqn(a,i) + sum1
        end do
     end do
  end do

  !deallocate( I5ab )
  
  !
  ! Subtract diagonal parts of I6 and I7
  !
  do a = below_ef+1, tot_orbs
     d1ac(a) = I6ac(a,a) 
     I6ac(a,a) = 0.d0

  end do

  do i = 1, below_ef
     d1ik(i) = I7abcde(i,i)
     I7abcde(i,i) = 0.d0
  end do

  
  
  !
  ! I6a and I7a contribution to t1 (ok)
  !
  do i = 1, below_ef
     do a = below_ef+1, tot_orbs
        
        if ( all_orbit%szp(a) /= all_orbit%szp(i) ) cycle
        if ( all_orbit%itzp(a) /= all_orbit%itzp(i) ) cycle

        t1_ccm_eqn(a,i) = t1_ccm_eqn(a,i) + & 
             DDOT(size(i6ac,2),i6ac(a,:),1,t1_ccm(:,i),1) + & 
             DDOT(size(t1_ccm,2),t1_ccm(a,:),1,i7abcde(:,i),1)
     end do
  end do
  deallocate( I6ac )
 
!!$  !
!!$  ! I4a contribution to t1(ok)
!!$  !
!!$  do i1 = 1, channels%number_confs 
!!$     ang_mom   = channels%config_J_ipar_tz(i1*3)
!!$     p_parity  = channels%config_J_ipar_tz(i1*3-1)
!!$     isospin_z = channels%config_J_ipar_tz(i1*3-2)
!!$       
!!$     
!!$     channel = locate_gmatchannel(4,isospin_z, p_parity, ang_mom)
!!$     if ( channel == 0 ) cycle
!!$     
!!$     angmom_fact = 2.d0*ang_mom + 1.d0 
!!$     !$omp parallel default(shared) private(bra, k, a, ket,i,c,phase,gmat, sum1 ) 
!!$     !$omp do schedule(static), reduction(+:t1_ccm_eqn) 
!!$     DO bra=1, size(  lookup_ab_configs(2,i1)%ival, 2) 
!!$        k = lookup_ab_configs(2,i1)%ival(1,bra) 
!!$        a = lookup_ab_configs(2,i1)%ival(2,bra) 
!!$       
!!$        DO ket=1, size(  lookup_ab_configs(2,i1)%ival, 2) 
!!$           i = lookup_ab_configs(2,i1)%ival(1,ket) 
!!$           c = lookup_ab_configs(2,i1)%ival(2,ket)  
!!$
!!$           
!!$           if ( all_orbit%jj(a) /= all_orbit%jj(i) ) cycle
!!$           if ( all_orbit%itzp(a) /= all_orbit%itzp(i) ) cycle
!!$        
!!$           phase = (-1.d0)**( (all_orbit%jj(c) + all_orbit%jj(i))/2 - ang_mom + 1) 
!!$           gmat = gmat_mtx_hphp%depot(channel)%val(bra,ket) * phase 
!!$           sum1 = angmom_fact * gmat * t1_ccm(c,k) 
!!$           t1_ccm_eqn(a,i) = t1_ccm_eqn(a,i) + sum1/(all_orbit%jj(a) + 1.d0)
!!$           
!!$        end DO
!!$     end DO
!!$     !$omp end do                                                                                                                                      
!!$     !$omp end parallel 
!!$     
!!$  end do
!!$  
!!$  !
!!$  !Add fock matrix to T1 
!!$  ! 
!!$  t1_ccm_eqn(below_ef+1:tot_orbs, 1:below_ef) = t1_ccm_eqn(below_ef+1:tot_orbs, 1:below_ef) + & 
!!$       fock_mtx(below_ef+1:tot_orbs, 1:below_ef)
!!$
!!$  !
!!$  ! Setup energy denominator
!!$  !
!!$  t1_denom = 0.d0
!!$  do i = 1, below_ef
!!$     do a = below_ef+1, tot_orbs
!!$          
!!$        t1_denom(a,i) = t1_denom(a,i) + d1ik(i) - d1ac(a)  
!!$                
!!$     end do
!!$  end do
!!$  
!!$  !endwtime = MPI_WTIME()
!!$  !call mpi_barrier(mpi_comm_world,ierror)
!!$  !if ( iam == 0 ) write(6,*) 'Total execution time T1 intermediates', endwtime - startwtime
!!$  

end SUBROUTINE t1_intermediate




!
! This routine has been checked with Davids, and is correct.
!
SUBROUTINE setup_I3_intermediate
  
  USE single_particle_orbits
  USE configurations
  USE constants
  use t2_storage
  use one_body_operators
  
  IMPLICIT NONE
  INTEGER :: ket_confs, bra_confs
  INTEGER :: i,number_channels 
  INTEGER :: bra, ket, c, d, k, l, i1
  REAL*8 :: gmat, sum1 !dum, denom, sum1,sum2,t2, angmom_fact, phase, interm
  !INTEGER :: local_ch_id, bra_min, bra_max, ket_min, ket_max, ndim2, ndim
  real*8 :: vmom_yukawa

  !
  ! I3ab intermediate
  !
  number_channels = channels%number_hhhp_confs 
  allocate( I3ab(1:number_channels) )
  do i=1, number_channels
     ket_confs = size( lookup_hhhp_configs(2,i)%ival, 2 )
     bra_confs = size( lookup_hhhp_configs(1,i)%ival, 2 )
     
     allocate( I3ab(i)%val(0:bra_confs,0:ket_confs) )
     i3ab(i)%val = 0.d0 
  end do
  
  
  do i1 = 1, channels%number_hhhp_confs 
     
     DO ket= 1, size(  lookup_hhhp_configs(2,i1)%ival, 2) 
        i = lookup_hhhp_configs(2,i1)%ival(1,ket) 
        c = lookup_hhhp_configs(2,i1)%ival(2,ket) 
        
        DO bra=1, size(  lookup_hhhp_configs(1,i1)%ival, 2) 
           k = lookup_hhhp_configs(1,i1)%ival(1,bra) 
           l = lookup_hhhp_configs(1,i1)%ival(2,bra) 
           
           I3ab(i1)%val(bra, ket) = vmom_yukawa(k,l,i,c) 
           sum1 = 0.d0 
           do d = below_ef+1, tot_orbs
              
              gmat = vmom_yukawa(k,l,c,d) 
              sum1 = - gmat * t1_ccm(d,i) 
              
           end do
           I3ab(i1)%val(bra, ket) = I3ab(i1)%val(bra, ket) + sum1  

        end DO
     end DO
  end do
  
  

end SUBROUTINE setup_I3_intermediate
!
!
!
SUBROUTINE setup_I5I6I7_intermediate
  
  USE single_particle_orbits
  USE configurations
  USE constants
  use t2_storage
  use one_body_operators
  
  IMPLICIT NONE

  INTEGER :: i,a,i1, c,d,k,l,bra,ket, channel  
  REAL*8 :: gmat,t2, sum1
  double precision :: ddot
  real*8 :: vmom_yukawa
  integer :: tz, sz, nx,ny,nz

  !
  ! Setup I6ac intermediate
  !
  allocate( I6ac(below_ef+1:tot_orbs, below_ef+1:tot_orbs) ) 
  allocate( I6abc(below_ef+1:tot_orbs, below_ef+1:tot_orbs) ) 
  allocate( I7abcde( 1:below_ef, 1:below_ef) )
  allocate( I5ab(1:below_ef, below_ef+1:tot_orbs) )
  
  !
  ! I5b  intermediate
  !
  i6abc = 0.d0 
  I7abcde = 0.d0
  I6ac = 0.d0
  I5ab = 0.d0
  write(6,*) 'ok1'
  do i1 = 1, channels%number_hhpp_confs 
     
     DO ket= 1, size(  lookup_hhpp_configs(2,i1)%ival, 2) 
        c = lookup_hhpp_configs(2,i1)%ival(1,ket) 
        d = lookup_hhpp_configs(2,i1)%ival(2,ket) 
        
        DO bra=1, size(  lookup_hhpp_configs(1,i1)%ival, 2) 
           k = lookup_hhpp_configs(1,i1)%ival(1,bra) 
           l = lookup_hhpp_configs(1,i1)%ival(2,bra) 
           
           gmat =  vmom_yukawa(c,d,k,l)
           I5ab(k,c) = I5ab(k,c) + gmat * t1_ccm(d,l)
           
        end DO
     end DO
  end do
  
  write(6,*) 'ok2'
  
  do channel   = 1, channels%number_hppp_confs 
     nx = channels%hppp_quantum_numbers(channel*5)
     ny = channels%hppp_quantum_numbers(channel*5-1)
     nz = channels%hppp_quantum_numbers(channel*5-2)
     sz = channels%hppp_quantum_numbers(channel*5-3)
     tz = channels%hppp_quantum_numbers(channel*5-4)
     
     
     do bra = 1, size(  lookup_hppp_configs(1,channel)%ival, 2) 
        k = lookup_hppp_configs(1,channel)%ival(1,bra) 
        d = lookup_hppp_configs(1,channel)%ival(2,bra) 
        do ket = 1, size(  lookup_hppp_configs(2,channel)%ival, 2) 
           a = lookup_hppp_configs(2,channel)%ival(1,ket)
           c = lookup_hppp_configs(2,channel)%ival(2,ket) 
                       
           sum1 = t1_ccm(d,k) * vmom_yukawa(k,a,d,c)
           I6ac(a,c) = I6ac(a,c) + sum1 
        end do
     end do
  end DO
  
  write(6,*) 'ok3'
  
  !
  ! Setup I7abcde intermediate
  !
  i7abcde = 0.d0 
  !
  ! I7(b) intermediate
  !
  do i1 = 1, channels%number_hhpp_confs 
     
     DO ket= 1, size(  lookup_hhpp_configs(2,i1)%ival, 2) 
        c = lookup_hhpp_configs(2,i1)%ival(1,ket) 
        d = lookup_hhpp_configs(2,i1)%ival(2,ket) 
        
        DO bra=1, size(  lookup_hhpp_configs(1,i1)%ival, 2) 
           i = lookup_hhpp_configs(1,i1)%ival(1,bra) 
           l = lookup_hhpp_configs(1,i1)%ival(2,bra) 
           
           
           do k = 1, below_ef 
              gmat =  vmom_yukawa(k,l,c,d)
              i7abcde(k,i) =  i7abcde(k,i) + 0.5d0 * gmat*t2_ccm(i1)%val(ket,bra) 
           end do
        end DO
     end DO
  end do
  
  write(6,*) 'ok4'
  !
  ! I7(c) intermediate
  !
  do channel   = 1, channels%number_hhhp_confs 
     nx = channels%hhhp_quantum_numbers(channel*5)
     ny = channels%hhhp_quantum_numbers(channel*5-1)
     nz = channels%hhhp_quantum_numbers(channel*5-2)
     sz = channels%hhhp_quantum_numbers(channel*5-3)
     tz = channels%hhhp_quantum_numbers(channel*5-4)
     
     
     do bra = 1, size(  lookup_hhhp_configs(1,channel)%ival, 2) 
        k = lookup_hhhp_configs(1,channel)%ival(1,bra) 
        l = lookup_hhhp_configs(1,channel)%ival(2,bra) 
        do ket = 1, size(  lookup_hhhp_configs(2,channel)%ival, 2) 
           i = lookup_hhhp_configs(2,channel)%ival(1,ket)
           c = lookup_hhhp_configs(2,channel)%ival(2,ket) 

           sum1 = vmom_yukawa(k,l,i,c) * t1_ccm(c,l) 
           I7abcde(k,i) = I7abcde(k,i) + sum1 
        end do
     end do
  end DO
  
  write(6,*) 'ok5'
  !
  ! I7(d,e) intermediate
  !
  do k = 1, below_ef
     do i = 1, below_ef
        
        I7abcde(k,i) =  I7abcde(k,i) + & 
             DDOT(size(i5ab,2),i5ab(k,below_ef+1:tot_orbs),1,t1_ccm(below_ef+1:tot_orbs,i),1)
        
     end do
  end do
  
  write(6,*) 'ok6'
  !
  ! Add fock matrix to I5,  I7 and I6 
  !
  I7abcde = I7abcde + fock_mtx(1:below_ef, 1:below_ef) 
  I6ac = I6ac + fock_mtx(below_ef+1:tot_orbs, below_ef+1:tot_orbs) 
  I5ab = I5ab + fock_mtx(1:below_ef, below_ef+1:tot_orbs) 
  
  !
  ! I6(abc) intermediate
  ! Here I6(b) is calculated, I6(a+c) calculated in t1-eqn
  !
  do i1 = 1, channels%number_hhpp_confs 
     
     DO bra=1, size( lookup_hhpp_configs(2,i1)%ival, 2) 
        d = lookup_hhpp_configs(2,i1)%ival(1,bra) 
        a = lookup_hhpp_configs(2,i1)%ival(2,bra) 
        
        DO ket=1, size(  lookup_hhpp_configs(1,i1)%ival, 2) 
           k = lookup_hhpp_configs(1,i1)%ival(1,ket) 
           l = lookup_hhpp_configs(1,i1)%ival(2,ket)  
           
           t2 = t2_ccm(i1)%val(bra,ket) 
           do c = below_ef+1, tot_orbs
              
              gmat = vmom_yukawa(k,l,d,c)
              sum1 = - 0.5d0* gmat * t2 
              I6abc(a,c) = I6abc(a,c) + sum1
              
           end do
        end DO
     end DO
  end do
  
  write(6,*) 'ok7'
  I6abc = i6abc + I6ac
  

end SUBROUTINE setup_I5I6I7_intermediate
