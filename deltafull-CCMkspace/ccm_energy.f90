!
! subroutine calculates the ccsd energy
!
SUBROUTINE ccd_energy_save(ener1)
  USE single_particle_orbits
  USE constants
  use one_body_operators
  use t2_storage
  use configurations
  use chiral_potentials

  IMPLICIT NONE
  
  complex*16, intent(inout) :: ener1 
  INTEGER :: i,j,a,b, bra2, ket2, nx, ny, nz, tz, channel3, tz2, channel2, bra3,ket3
  INTEGER :: bra,ket, channel, local_ch_id, bra_min, bra_max, ket_min, ket_max
  complex*16 ::  sum3, sum2, sum33
  INTEGER :: nx1, ny1, nz1, nx2, ny2, nz2, nx3, ny3, nz3, nx4, ny4, nz4 
  real*8  ::  startwtime , endwtime
  
  sum33  = 0.d0; sum2 = 0.d0 

  sum3 = 0.d0 

  do channel = my_hhpp_channel_low(iam), my_hhpp_channel_high(iam)
     
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
     !!$omp parallel default(shared) private(bra,i,j,ket,a,b)
     !!$omp do schedule(dynamic), reduction(+:sum3)
     DO bra= bra_min, bra_max !1, size( lookup_hhpp_configs(2,channel)%ival, 2) 
        a = lookup_hhpp_configs(2,channel)%ival(1,bra) 
        b = lookup_hhpp_configs(2,channel)%ival(2,bra) 
        
        DO ket=1, size(  lookup_hhpp_configs(1,channel)%ival, 2) 
           i = lookup_hhpp_configs(1,channel)%ival(1,ket) 
           j = lookup_hhpp_configs(1,channel)%ival(2,ket)  
           
           sum3 = sum3 + 0.25d0 * t2_ccm(channel)%val(bra,ket)*vnn_hhpp(channel)%val(ket,bra)
           
        end DO
     end DO
     !!$omp end do
     !!$omp end parallel
  end do
  
  call mpi_allreduce(sum3,sum2,1,mpi_complex16,mpi_sum, &
          mpi_comm_world,ierror)
  
  !call mpi_reduce(sum3,sum2,1,mpi_complex16,mpi_sum, &
  !     master,mpi_comm_world,ierror)
  
  sum3 = sum2 
  !call mpi_bcast(sum3,1,mpi_complex16,master, &
  !     mpi_comm_world,ierror)

  ener1 =  sum3 
  if ( iam == 0 ) write(6,*) 'CCD energy', ener1 + e0 
  
  !stop
  
END SUBROUTINE ccd_energy_save
