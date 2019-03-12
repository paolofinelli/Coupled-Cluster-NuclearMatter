!
! simple linear mixing 
! 
subroutine linear_t1(ener)
  USE parallel
  USE CONSTANTS
  USE one_body_operators
  USE t2_storage
  USE DIIS_MOD
  USE configurations
 
  IMPLICIT NONE 
  REAL*8 :: a1, b1 
  complex*16 :: d2f5d, ener
  INTEGER :: i,j,a,b,number_channels, channel, bra, ket
  
  
  a1 = 0.25d0
  b1 = 0.75d0
  do i=1,below_ef
     do a=below_ef+1,tot_orbs
        
        t1_ccm(a,i)= a1*t1_ccm(a,i) + b1* t1_ccm_eqn(a,i)/ ( t1_denom(a,i) + ener )
                
     end do
  end do
  
end subroutine linear_t1
!
! simple linear mixing 
! 
subroutine linear_t2(ener)
  USE parallel
  USE CONSTANTS
  USE one_body_operators
  USE t2_storage
  USE DIIS_MOD
  USE configurations
 
  IMPLICIT NONE 
  REAL*8 :: a1, b1 
  COMPLEX*16 :: d2f5d, ener
  INTEGER :: i,j,a,b,number_channels, channel, bra, ket
  
  
  a1 = 0.5d0
  b1 = 0.5d0
  
  do channel   = 1, channels%number_hhpp_confs  
     DO bra=1, size(  lookup_hhpp_configs(2,channel)%ival, 2) 
        a = lookup_hhpp_configs(2,channel)%ival(1,bra)
        b = lookup_hhpp_configs(2,channel)%ival(2,bra) 
        
        DO ket=1, size(  lookup_hhpp_configs(1,channel)%ival, 2) 
           i= lookup_hhpp_configs(1,channel)%ival(1,ket) 
           j = lookup_hhpp_configs(1,channel)%ival(2,ket) 
           
           d2f5d= ( d2ij(i) + d2ij(j) - d2ab(a)-d2ab(b) ) + ener 
           
           t2_ccm(channel)%val(bra,ket) = a1*t2_ccm(channel)%val(bra,ket) + & 
                b1*t2_ccm_eqn(channel)%val(bra,ket)/d2f5d
  
        end DO
     end DO
  end do
  
  
end subroutine linear_t2


!
! Diis mixing
!
subroutine diis_t1t2( count)
  USE parallel
  USE CONSTANTS
  USE one_body_operators
  USE t2_storage
  USE DIIS_MOD
  USE configurations
  
  IMPLICIT NONE 
  REAL*8 :: a1, b1 
  complex*16 :: d2f5d
  INTEGER :: i,j,k,a,b,number_channels, channel, bra, ket
  integer :: ii,i1, sub, count
  integer, allocatable,dimension(:) :: ipiv
  complex*16,  allocatable,dimension(:) :: work
  integer :: lwork,info, aa, bb, jj, conf
  double precision :: t1old, t2old
  integer :: nx2, ny2, nz2, tz2, channel2, bra2,ket2
  
  a1 = 0.75d0
  b1 = 0.25d0
  
  
  
  if ( nstep > subspace ) then
        
     nstep = subspace
     do i = 1, subspace -1 
        t2_diis(:,i) = t2_diis(:,i+1)
     end do
  end if
  

  i1 = 0 
  do channel   = 1, channels%number_hhpp_confs  
     DO bra=1, size(  lookup_hhpp_configs(2,channel)%ival, 2) 
        a = lookup_hhpp_configs(2,channel)%ival(1,bra)
        b = lookup_hhpp_configs(2,channel)%ival(2,bra) 
        
        DO ket=1, size(  lookup_hhpp_configs(1,channel)%ival, 2) 
           i= lookup_hhpp_configs(1,channel)%ival(1,ket) 
           j = lookup_hhpp_configs(1,channel)%ival(2,ket) 
           
           d2f5d= ( d2ij(i) + d2ij(j) - d2ab(a)-d2ab(b) )
                 
           t2_ccm(channel)%val(bra,ket) = a1*t2_ccm(channel)%val(bra,ket) + & 
                b1*t2_ccm_eqn(channel)%val(bra,ket)/d2f5d
           
           i1=i1+1
           t2_diis(i1,nstep) = t2_ccm(channel)%val(bra,ket)
        end DO
     end DO
  end do
  
  
  if ( count == diis_step   ) then


     count = 0

     do i=1,nstep-1
        do k=1,n_diis
           t2_diisd(k,i)=t2_diis(k,i+1)-t2_diis(k,i)
        end do
     end do

     allocate(diis_mat(nstep,nstep))
     allocate(sol_vect(nstep))
     
     diis_mat=0.0
     do i=1,nstep-1
        do j=1,nstep-1
           do k=1,n_diis
              diis_mat(i,j)=diis_mat(i,j)+t2_diisd(k,i)*t2_diisd(k,j)
           end do
        end do
     end do
     do i=1,nstep-1
        diis_mat(i,nstep)=-1.0
        diis_mat(nstep,i)=-1.0
     end do
     diis_mat(nstep,nstep)=0.0
     
     if( iam == 0 )then
        do i=1,nstep
           do j=1,nstep
              write(6,'(a,2i5,2e12.4)')'diis_mat',i,j,real(diis_mat(i,j)) 
           end do
        end do
     end if

        

     !                                                                                                                             
     sol_vect=0.0
     sol_vect(nstep)=-1.0
     allocate(ipiv(nstep))
     lwork=nstep
     
     call ZGESV( nstep, 1, diis_mat, nstep, IPIV, sol_vect, nstep, INFO )
     
     if(iam == 0)then
        write(6,*)'csysv info',info
        do i=1,nstep
           write(6,'(a,i5,2e12.4)')'sol_vec',i,real(sol_vect(i)),aimag(sol_vect(i))
	end do
     end if
     
     
     !
     ! use t2_diisd(:,1) to store the p=sum_{i=1,m}c_i t2_diis(:,i)
     !
     do i1=1,n_diis
        t2_diisd(i1,1)=sol_vect(1)*t2_diis(i1,1)
     end do
     do i=2,nstep-1
        do i1=1,n_diis
           t2_diisd(i1,1)=t2_diisd(i1,1)+sol_vect(i)*t2_diis(i1,i)
        end do
     end do
     
     i1 = 0 
     do channel   = 1, channels%number_hhpp_confs  
        DO bra=1, size(  lookup_hhpp_configs(2,channel)%ival, 2) 
           a = lookup_hhpp_configs(2,channel)%ival(1,bra)
           b = lookup_hhpp_configs(2,channel)%ival(2,bra) 
           
           DO ket=1, size(  lookup_hhpp_configs(1,channel)%ival, 2) 
              i= lookup_hhpp_configs(1,channel)%ival(1,ket) 
              j = lookup_hhpp_configs(1,channel)%ival(2,ket) 

              i1=i1+1
              
              d2f5d= ( d2ij(i) + d2ij(j) - d2ab(a)-d2ab(b) )
              
              t2_ccm(channel)%val(bra,ket) = t2_diisd(i1,1)                     
           end DO
        end DO
     end do
     
     deallocate(diis_mat)
     deallocate(sol_vect)
     deallocate(ipiv)
     !deallocate(work)
  end if
  
!!$
!!$  do channel   = 1, channels%number_ph_phhp_confs 
!!$     
!!$     t2_ph_ccm(channel)%val = 0.d0 
!!$     do bra = 1, size(  lookup_ph_phhp_configs(1,channel)%ival, 2) 
!!$        a = lookup_ph_phhp_configs(1,channel)%ival(1,bra) 
!!$        i = lookup_ph_phhp_configs(1,channel)%ival(2,bra) 
!!$        
!!$        do ket = 1, size(  lookup_ph_phhp_configs(2,channel)%ival, 2) 
!!$           j = lookup_ph_phhp_configs(2,channel)%ival(1,ket)
!!$           b = lookup_ph_phhp_configs(2,channel)%ival(2,ket) 
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
!!$           t2_ph_ccm(channel)%val(bra,ket) = t2_ccm(channel2)%val(ket2,bra2)
!!$           
!!$        end do
!!$     end do
!!$  end do
  
     
end subroutine diis_t1t2



!
! performs a diis step 
! 

subroutine diis_setup
  USE parallel
  USE CONSTANTS
  USE one_body_operators
  USE t2_storage
  use diis_mod
  USE configurations
  
  implicit none
  integer :: ii, i,j,a,b, bra,ket, channel
  
  ii=0
  do channel   = 1, channels%number_hhpp_confs  
     DO bra=1, size(  lookup_hhpp_configs(2,channel)%ival, 2) 
        a = lookup_hhpp_configs(2,channel)%ival(1,bra)
        b = lookup_hhpp_configs(2,channel)%ival(2,bra) 
        
        DO ket=1, size(  lookup_hhpp_configs(1,channel)%ival, 2) 
           i= lookup_hhpp_configs(1,channel)%ival(1,ket) 
           j = lookup_hhpp_configs(1,channel)%ival(2,ket) 
           
           ii = ii + 1 

        end DO
     end DO
  end do
  
  n_diis=ii
  nn_diis= subspace
  if(iam == 0) write(6,*)'diis space',n_diis
  allocate(t2_diis(n_diis,nn_diis+1))
  allocate(t2_diisd(n_diis,nn_diis))
   
  if ( iam == 0 ) write(6,*) 'Storage of diis vectors', 2.d0*n_diis*(nn_diis+1.d0)*8.d0/1.e9,'Gb'
    

end subroutine diis_setup


subroutine diis_take_down

  USE parallel
  USE CONSTANTS
  USE one_body_operators
  USE t2_storage
  use diis_mod

  implicit none
  
  deallocate(t2_diis)
  deallocate(t2_diisd)
  
end subroutine diis_take_down
