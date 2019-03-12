subroutine t3_intermediate
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
  INTEGER :: i,j,k, sumx, sumy, sumz, a,b,c,nxa,nxb,nxc,nya,nyb,nyc,nza,nzb,nzc,tza,tzb,tzc,sza,szb,szc,tz,sz,nx,ny,nz
  INTEGER :: numchannels, ij_confs, ab_confs, channel, sz_min, sz_max, tz_min, tz_max, sumtz, sumsz, ichan
  INTEGER :: number_channels, k1,k2,k3,k4,k5,bra2,ket2
  INTEGER :: local_ch_id, dim1, dim2, bra_min, bra_max, ket_min, ket_max, bra, ket
  COMPLEX*16 :: spen, t3_ener
  complex*16 :: sum1,sum2,sum3,sum4,sum5,sum6, sum7,sum8,sum9,tc_term, td_term 
  complex*16, allocatable :: amat(:,:)
  
   
  do channel = my_t3channel_low(iam), my_t3channel_high(iam)
     t3_ccm(channel)%val = 0.D0
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
     
     
!!$     ppp_hhhppp%ival3 = 0 
!!$     do bra = 1, size(lookup_t3_configs(2,channel)%ival,2)
!!$        a = lookup_t3_configs(2,channel)%ival(1,bra)
!!$        b = lookup_t3_configs(2,channel)%ival(2,bra) 
!!$        c = lookup_t3_configs(2,channel)%ival(3,bra) 
!!$        
!!$        ppp_hhhppp%ival3(a,b,(all_orbit%szp(c)+1)/2) = bra
!!$     end do
!!$     
     hhh_hhhppp%ival3  = 0
     do ket = 1, size(  lookup_t3_configs(1,channel)%ival, 2) 
        i = lookup_t3_configs(1,channel)%ival(1,ket) 
        j = lookup_t3_configs(1,channel)%ival(2,ket) 
        k = lookup_t3_configs(1,channel)%ival(3,ket) 
        hhh_hhhppp%ival3(i,j,k) = ket 
     end do
     
     
     !$omp parallel default(shared) private(bra,a,b,c,ket,i,j,k,ket2,spen, & 
     !$omp sum1,sum2,sum3,sum4,sum5,sum6,sum7,sum8,sum9,tc_term)
     !$omp do schedule(dynamic) 
     do bra = bra_min, bra_max  !1, size(lookup_t3_configs(2,channel)%ival,2) ! 
        a = lookup_t3_configs(2,channel)%ival(1,bra)
        b = lookup_t3_configs(2,channel)%ival(2,bra) 
        c = lookup_t3_configs(2,channel)%ival(3,bra) 
        
        do ket = 1, size(  lookup_t3_configs(1,channel)%ival, 2) 
           i = lookup_t3_configs(1,channel)%ival(1,ket) 
           j = lookup_t3_configs(1,channel)%ival(2,ket) 
           k = lookup_t3_configs(1,channel)%ival(3,ket) 
           
           spen = ( fock_mtx(i,i) + fock_mtx(j,j)  + fock_mtx(k,k) & 
                - fock_mtx(a,a)-fock_mtx(b,b) - fock_mtx(c,c) )
           
           !t3_ccm(channel)%val(bra,ket) = t3_ccm(channel)%val(bra,ket) + (v3nf_ppphhh(channel)%val(bra,ket))/spen 
           t3_ccm(channel)%val(bra,ket) = (v3nf_ppphhh(channel)%val(bra,ket))/spen 
           tc_term = 0.d0 
           
           sum1 = 0.d0; sum2 = 0.d0; sum3 = 0.d0
           sum4 = 0.d0; sum5 = 0.d0; sum6 = 0.d0
           sum7 = 0.d0; sum8 = 0.d0; sum9 = 0.d0
           ! 
           ! P(a/bc)P(k/ij) <bc|V|dk><ad|t|ij>
           !  1    
           call three_body_oneparticle_1(a,b,c,i,j,k,sum1)
           ! (ab)
           call three_body_oneparticle_1(b,a,c,i,j,k,sum2)
           ! (ac)
           call three_body_oneparticle_1(c,b,a,i,j,k,sum3)
           ! (ik)
           call three_body_oneparticle_1(a,b,c,k,j,i,sum4)
           ! (kj)
           call three_body_oneparticle_1(a,b,c,i,k,j,sum5)
           ! (ac)(jk)
           call three_body_oneparticle_1(c,b,a,i,k,j,sum6)
           ! (ac)(ik)
           call three_body_oneparticle_1(c,b,a,k,j,i,sum7)
           ! (ab)(jk)
           call three_body_oneparticle_1(b,a,c,i,k,j,sum8)
           ! (ab)(ik)
           call three_body_oneparticle_1(b,a,c,k,j,i,sum9)
           
           tc_term = ( sum1 - sum2 - sum3  - sum4 - sum5 + sum6 + sum7 + sum8 + sum9) 
           
!!$           
!!$           t3_ccm(channel)%val(bra,ket) = t3_ccm(channel)%val(bra,ket) + (sum1-sum2-sum3) /spen 
!!$           ! (ik) 
!!$           ket2 =  hhh_hhhppp%ival3(k,j,i)
!!$           t3_ccm(channel)%val(bra,ket2) = t3_ccm(channel)%val(bra,ket2) - sum1/spen 
!!$           ! (kj) 
!!$           ket2 =  hhh_hhhppp%ival3(i,k,j)
!!$           t3_ccm(channel)%val(bra,ket2) = t3_ccm(channel)%val(bra,ket2) - sum1/spen 
!!$           ! (ac)(kj) 
!!$           ket2 =  hhh_hhhppp%ival3(i,k,j)
!!$           t3_ccm(channel)%val(bra,ket2) = t3_ccm(channel)%val(bra,ket2) + sum3/spen 
!!$           ! (ac)(ik) 
!!$           ket2 =  hhh_hhhppp%ival3(k,j,i)
!!$           t3_ccm(channel)%val(bra,ket2) = t3_ccm(channel)%val(bra,ket2) + sum3/spen 
!!$           ! (ab)(kj) 
!!$           ket2 =  hhh_hhhppp%ival3(i,k,j)
!!$           t3_ccm(channel)%val(bra,ket2) = t3_ccm(channel)%val(bra,ket2) + sum2/spen 
!!$           ! (ab)(ik) 
!!$           ket2 =  hhh_hhhppp%ival3(k,j,i)
!!$           t3_ccm(channel)%val(bra,ket2) = t3_ccm(channel)%val(bra,ket2) + sum2/spen 
!!$           
!!$           
           
           
           sum1 = 0.d0; sum2 = 0.d0; sum3 = 0.d0
           sum4 = 0.d0; sum5 = 0.d0; sum6 = 0.d0
           sum7 = 0.d0; sum8 = 0.d0; sum9 = 0.d0
           ! 
           ! P(c/ab)P(i/jk) <ab|t|il><lc|v|jk>
           !  1    ok 
           call three_body_onehole_1(a,b,c,i,j,k,sum1)
           ! (cb)
           call three_body_onehole_1(a,c,b,i,j,k,sum2)
           ! (ac)
           call three_body_onehole_1(c,b,a,i,j,k,sum3)
           ! (ik)
           call three_body_onehole_1(a,b,c,k,j,i,sum4)
           ! (ij)
           call three_body_onehole_1(a,b,c,j,i,k,sum5)
           ! (ac)(ij)
           call three_body_onehole_1(c,b,a,j,i,k,sum6)
           ! (ac)(ik)
           call three_body_onehole_1(c,b,a,k,j,i,sum7)
           ! (cb)(ij)
           call three_body_onehole_1(a,c,b,j,i,k,sum8)
           ! (cb)(ik)
           call three_body_onehole_1(a,c,b,k,j,i,sum9)
           
           tc_term = tc_term -( sum1 - sum2 - sum3  - sum4 - sum5 + sum6 + sum7 + sum8 + sum9 )
           !t3_ccm(channel)%val(bra,ket) =  ( tc_term+v3nf_ppphhh(channel)%val(bra,ket)) /spen 
           !t3_ccm(channel)%val(bra,ket) = t3_ccm(channel)%val(bra,ket) + tc_term /spen 
           t3_ccm(channel)%val(bra,ket) = t3_ccm(channel)%val(bra,ket) + tc_term /spen 
           

!!$           
!!$           t3_ccm(channel)%val(bra,ket) = t3_ccm(channel)%val(bra,ket) - (sum1-sum2-sum3) /spen 
!!$           ! (ik) 
!!$           ket2 =  hhh_hhhppp%ival3(k,j,i)
!!$           t3_ccm(channel)%val(bra,ket2) = t3_ccm(channel)%val(bra,ket2) + sum1/spen 
!!$           ! (ij) 
!!$           ket2 =  hhh_hhhppp%ival3(j,i,k)
!!$           t3_ccm(channel)%val(bra,ket2) = t3_ccm(channel)%val(bra,ket2) + sum1/spen 
!!$           ! (ac)(ij) 
!!$           ket2 =  hhh_hhhppp%ival3(j,i,k)
!!$           t3_ccm(channel)%val(bra,ket2) = t3_ccm(channel)%val(bra,ket2) - sum3/spen 
!!$           ! (ac)(ik) 
!!$           ket2 =  hhh_hhhppp%ival3(k,j,i)
!!$           t3_ccm(channel)%val(bra,ket2) = t3_ccm(channel)%val(bra,ket2) - sum3/spen 
!!$           ! (cb)(ij) 
!!$           ket2 =  hhh_hhhppp%ival3(j,i,k)
!!$           t3_ccm(channel)%val(bra,ket2) = t3_ccm(channel)%val(bra,ket2) - sum2/spen 
!!$           ! (cb)(ik) 
!!$           ket2 =  hhh_hhhppp%ival3(k,j,i)
!!$           t3_ccm(channel)%val(bra,ket2) = t3_ccm(channel)%val(bra,ket2) - sum2/spen 

        end DO
     end DO
     !$omp end do
     !$omp end parallel
 
  end do


end subroutine t3_intermediate

