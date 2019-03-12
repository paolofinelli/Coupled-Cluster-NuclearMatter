!
!     <abc|V|ide><de|t|jk>
!
subroutine three_body_v3t2pp(a,b,c,i,j,k,sum1)

  USE single_particle_orbits
  USE configurations
  USE constants
  use t2_storage
  use one_body_operators
  use parallel 
  use chiral_potentials 
  
  IMPLICIT NONE
  complex*16, INTENT(OUT) :: sum1 
  INTEGER, intent(in) :: i,j,k,a,b,c
  INTEGER :: e, d, nxd, nyd, nzd,szd,tzd 
  INTEGER :: nx2, ny2, nz2, sz2, tz2, channel2,channel5
  INTEGER :: nx3, ny3, nz3, sz3, tz3, bra2,ket2,bra3,ket3
  
  sum1 = 0.d0 
  nx2 = all_orbit%nx(k) + all_orbit%nx(j)
  ny2 = all_orbit%ny(k) + all_orbit%ny(j)
  nz2 = all_orbit%nz(k) + all_orbit%nz(j)
  tz2 = (all_orbit%itzp(k) + all_orbit%itzp(j))/2
  channel2 = locate_channel(3,tz2, nx2, ny2, nz2) 
  
  if ( channel2 == 0 ) return 
  ket2 = hh_hhpp%ival(j,k)
  if ( ket2 == 0 ) return 
  
  do e = below_ef+1, tot_orbs
     
     nxd = all_orbit%nx(j) + all_orbit%nx(k) - all_orbit%nx(e)
     nyd = all_orbit%ny(j) + all_orbit%ny(k) - all_orbit%ny(e) 
     nzd = all_orbit%nz(j) + all_orbit%nz(k) - all_orbit%nz(e) 
     tzd = all_orbit%itzp(j) + all_orbit%itzp(k) - all_orbit%itzp(e)
     
     if ( abs(nxd) > nkxmax ) CYCLE 
     if ( abs(nyd) > nkymax ) CYCLE 
     if ( abs(nzd) > nkzmax ) CYCLE 
     if ( abs(tzd) > 1 ) CYCLE  
     do szd = -1, 1 ,2 
     
        d = orbitof(nxd, nyd, nzd, szd, tzd)
        if ( d <= below_ef ) cycle 
        
        bra2 = pp_hhpp%ival(d,e)
        if ( bra2 == 0 ) cycle
        
        sum1 = sum1 + 0.5d0 * chiral_3nf_asym(a,b,c,i,d,e)*t2_ccm(channel2)%val(bra2,ket2)
     end do
  end do
  
end subroutine three_body_v3t2pp

!
!     <alm|V|ijk><bc|t|lm>
!
subroutine three_body_v3t2hh(a,b,c,i,j,k,sum1)

  USE single_particle_orbits
  USE configurations
  USE constants
  use t2_storage
  use one_body_operators
  use parallel 
  use chiral_potentials 
  
  IMPLICIT NONE
  complex*16, INTENT(OUT) :: sum1 
  INTEGER, intent(in) :: i,j,k,a,b,c
  INTEGER :: m,l, nxd, nyd, nzd,szd,tzd 
  INTEGER :: nx2, ny2, nz2, sz2, tz2, channel2,channel5
  INTEGER :: nx3, ny3, nz3, sz3, tz3, bra2,ket2,bra3,ket3
  
  sum1 = 0.d0 
  nx2 = all_orbit%nx(b) + all_orbit%nx(c)
  ny2 = all_orbit%ny(b) + all_orbit%ny(c)
  nz2 = all_orbit%nz(b) + all_orbit%nz(c)
  tz2 = (all_orbit%itzp(b) + all_orbit%itzp(c))/2
  channel2 = locate_channel(3,tz2, nx2, ny2, nz2) 
  
  if ( channel2 == 0 ) return 
  bra2 = pp_hhpp%ival(b,c)
  if ( bra2 == 0 ) return 
  
  do m = 1, below_ef
     
     nxd = all_orbit%nx(b) + all_orbit%nx(c) - all_orbit%nx(m)
     nyd = all_orbit%ny(b) + all_orbit%ny(c) - all_orbit%ny(m) 
     nzd = all_orbit%nz(b) + all_orbit%nz(c) - all_orbit%nz(m) 
     tzd = all_orbit%itzp(b) + all_orbit%itzp(c) - all_orbit%itzp(m)
     
     if ( abs(nxd) > nkxmax ) CYCLE 
     if ( abs(nyd) > nkymax ) CYCLE 
     if ( abs(nzd) > nkzmax ) CYCLE 
     if ( abs(tzd) > 1 ) CYCLE  
     do szd = -1, 1 ,2 
     
        l = orbitof(nxd, nyd, nzd, szd, tzd)
        if ( l > below_ef ) cycle 
        
        ket2 = hh_hhpp%ival(l,m)
        if ( ket2 == 0 ) cycle
        
        sum1 = sum1 + 0.5d0 * chiral_3nf_asym(a,l,m,i,j,k)*t2_ccm(channel2)%val(bra2,ket2)
     end do
  end do
  
end subroutine three_body_v3t2hh

!
!     <abl|V|ijd><dc|t|lk>
!
subroutine three_body_v3t2ph(a,b,c,i,j,k,sum1)

  USE single_particle_orbits
  USE configurations
  USE constants
  use t2_storage
  use one_body_operators
  use parallel 
  use chiral_potentials 
  
  IMPLICIT NONE
  complex*16, INTENT(OUT) :: sum1 
  INTEGER, intent(in) :: i,j,k,a,b,c
  INTEGER :: d, l, nxd, nyd, nzd,szd,tzd 
  INTEGER :: nx2, ny2, nz2, sz2, tz2, channel2,channel5
  INTEGER :: nx3, ny3, nz3, sz3, tz3, bra2,ket2,bra3,ket3
  
  sum1 = 0.d0 
  do l = 1, below_ef
  
     nx2 = all_orbit%nx(l) + all_orbit%nx(k)
     ny2 = all_orbit%ny(l) + all_orbit%ny(k)
     nz2 = all_orbit%nz(l) + all_orbit%nz(k)
     tz2 = (all_orbit%itzp(l) + all_orbit%itzp(k))/2
     channel2 = locate_channel(3,tz2, nx2, ny2, nz2) 
     
     if ( channel2 == 0 ) cycle
     ket2 = hh_hhpp%ival(l,k)
     if ( ket2 == 0 ) cycle 
  
  
     nxd = all_orbit%nx(l) + all_orbit%nx(k) - all_orbit%nx(c)
     nyd = all_orbit%ny(l) + all_orbit%ny(k) - all_orbit%ny(c) 
     nzd = all_orbit%nz(l) + all_orbit%nz(k) - all_orbit%nz(c) 
     tzd = all_orbit%itzp(l) + all_orbit%itzp(k) - all_orbit%itzp(c)
     
     if ( abs(nxd) > nkxmax ) CYCLE 
     if ( abs(nyd) > nkymax ) CYCLE 
     if ( abs(nzd) > nkzmax ) CYCLE 
     if ( abs(tzd) > 1 ) CYCLE  
     do szd = -1, 1 ,2 
     
        d = orbitof(nxd, nyd, nzd, szd, tzd)
        if ( d <= below_ef ) cycle 
        
        bra2 = pp_hhpp%ival(d,c)
        if ( bra2 == 0 ) cycle
        if ( abs(t2_ccm(channel2)%val(bra2,ket2) ) < 1.e-8 ) cycle
        
        sum1 = sum1 + chiral_3nf_asym(a,b,l,i,j,d)*t2_ccm(channel2)%val(bra2,ket2)
     end do
  end do
  
end subroutine three_body_v3t2ph

!
!     <cb|V|kd><ad|t|ij>
!
subroutine three_body_oneparticle_1(a,b,c,i,j,k,sum1)

  USE single_particle_orbits
  USE configurations
  USE constants
  use t2_storage
  use one_body_operators
  use parallel 
  use chiral_potentials 
  
  IMPLICIT NONE
  complex*16, INTENT(OUT) :: sum1 
  INTEGER, intent(in) :: i,j,k,a,b,c
  INTEGER :: d, nxd, nyd, nzd,szd,tzd 
  INTEGER :: nx2, ny2, nz2, sz2, tz2, channel2,channel5
  INTEGER :: nx3, ny3, nz3, sz3, tz3, bra2,ket2,bra3,ket3
  
  sum1 = 0.d0 
  nxd = all_orbit%nx(i) + all_orbit%nx(j) - all_orbit%nx(a)
  nyd = all_orbit%ny(i) + all_orbit%ny(j) - all_orbit%ny(a) 
  nzd = all_orbit%nz(i) + all_orbit%nz(j) - all_orbit%nz(a) 
  tzd = all_orbit%itzp(i) + all_orbit%itzp(j) - all_orbit%itzp(a)
  
  if ( abs(nxd) > nkxmax ) RETURN 
  if ( abs(nyd) > nkymax ) RETURN 
  if ( abs(nzd) > nkzmax ) RETURN 
  if ( abs(tzd) > 1 ) RETURN  
  
  do szd = -1, 1 ,2 
     
     d = orbitof(nxd, nyd, nzd, szd, tzd)
     if ( d <= below_ef ) cycle 
     
     nx2 = all_orbit%nx(i) + all_orbit%nx(j)
     ny2 = all_orbit%ny(i) + all_orbit%ny(j)
     nz2 = all_orbit%nz(i) + all_orbit%nz(j)
     tz2 = (all_orbit%itzp(i) + all_orbit%itzp(j))/2
     channel2 = locate_channel(3,tz2, nx2, ny2, nz2) 
     
     if ( channel2 == 0 ) cycle
     ket2 = hh_hhpp%ival(i,j)
     bra2 = pp_hhpp%ival(a,d)
     if ( bra2*ket2 == 0 ) cycle
  
     nx3 = all_orbit%nx(b) + all_orbit%nx(c)
     ny3 = all_orbit%ny(b) + all_orbit%ny(c)
     nz3 = all_orbit%nz(b) + all_orbit%nz(c)
     tz3 = (all_orbit%itzp(b) + all_orbit%itzp(c))/2
     channel5 = locate_channel(5,tz3, nx3, ny3, nz3) 
     bra3 = hp_hppp%ival(k,d)
     ket3 = pp_hppp%ival(c,b)
     if ( channel5 == 0 ) cycle
     if ( bra3*ket3 == 0 ) cycle
     
     sum1 = sum1 + conjg(vnn_hppp(channel5)%val(bra3,ket3))*t2_ccm(channel2)%val(bra2,ket2)
     !sum1 = sum1 + (chiral_pot(k,d,c,b)-chiral_pot(k,d,b,c))*t2_ccm(channel2)%val(bra2,ket2)
  end do
     
end subroutine three_body_oneparticle_1


!
! <lc|V|jk><ab|t|il> = conjg(<jk|V|lc>)<ab|t|il> 
!
subroutine three_body_onehole_1(a,b,c,i,j,k,sum1)
  
  USE single_particle_orbits
  USE configurations
  USE constants
  use t2_storage
  use one_body_operators
  use parallel 
  use chiral_potentials 
  
  IMPLICIT NONE
  complex*16, INTENT(OUT) :: sum1 
  INTEGER, intent(in) :: i,j,k,a,b,c
  INTEGER :: l, nxl, nyl, nzl,szl,tzl 
  INTEGER :: nx2, ny2, nz2, sz2, tz2, channel2,channel5
  INTEGER :: nx3, ny3, nz3, sz3, tz3, bra2,ket2,bra3,ket3
  
  sum1 = 0.d0 
  nxl = all_orbit%nx(k) + all_orbit%nx(j) - all_orbit%nx(c)
  nyl = all_orbit%ny(k) + all_orbit%ny(j) - all_orbit%ny(c) 
  nzl = all_orbit%nz(k) + all_orbit%nz(j) - all_orbit%nz(c) 
  tzl = all_orbit%itzp(k) + all_orbit%itzp(j) - all_orbit%itzp(c)
  if ( abs(nxl) > nkxmax ) RETURN 
  if ( abs(nyl) > nkymax ) RETURN 
  if ( abs(nzl) > nkzmax ) RETURN 
  if ( abs(tzl) > 1 ) RETURN  
  
  do szl = -1, 1, 2 
     l = orbitof(nxl, nyl, nzl, szl, tzl)
     if ( l > below_ef ) cycle
     
     nx2 = all_orbit%nx(i) + all_orbit%nx(l)
     ny2 = all_orbit%ny(i) + all_orbit%ny(l)
     nz2 = all_orbit%nz(i) + all_orbit%nz(l)
     tz2 = (all_orbit%itzp(i) + all_orbit%itzp(l))/2
     channel2 = locate_channel(3,tz2, nx2, ny2, nz2) 
     if ( channel2 == 0 )  cycle
     ket2 = hh_hhpp%ival(i,l)
     bra2 = pp_hhpp%ival(a,b)
     if ( bra2*ket2 == 0 ) cycle
     
     ! <jk|V|lc><ab|t|il>
     nx3 = all_orbit%nx(l) + all_orbit%nx(c)
     ny3 = all_orbit%ny(l) + all_orbit%ny(c)
     nz3 = all_orbit%nz(l) + all_orbit%nz(c)
     tz3 = (all_orbit%itzp(l) + all_orbit%itzp(c))/2
     channel5 = locate_channel(2,tz3, nx3, ny3, nz3) 
     bra3 = hh_hhhp%ival(j,k)
     ket3 = hp_hhhp%ival(l,c)
     if ( channel5 == 0 ) cycle
     if ( bra3*ket3 == 0 ) cycle
     
     sum1 = sum1 + conjg(vnn_hhhp(channel5)%val(bra3,ket3))*t2_ccm(channel2)%val(bra2,ket2)
     !sum1 = sum1 + (chiral_pot(j,k,l,c)-chiral_pot(j,k,c,l)) *t2_ccm(channel2)%val(bra2,ket2)
  end do
  
end subroutine three_body_onehole_1



subroutine pert_triples
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
  INTEGER :: i,j,k,l,a,b,c,d, number_channels, nx,ny,nz,sz,tz,  sumx, sumy, sumz,k1,k2,k3
  INTEGER :: bra_confs, ket_confs, bra,ket, ii, channel
  INTEGER :: sz_min, sz_max, tz_min, tz_max, sumtz, sumsz
  complex*16  :: spen, et_4tmp,et_4tmp2, et_4, et_4t, tmp_sum1, tmp_sum2, et_5,  & 
       est_5, sum1,sum2,sum3,sum4,sum5,sum6, sum7,sum8,sum9,tc_term, td_term,lc_term
  integer, allocatable :: nconf_low(:), nconf_high(:)
  integer :: nconfs, nconfs_tot, number_mtxel_iam, diff, n
  real*8  ::  startwtime , endwtime
  complex*16 :: diag1 
  
  startwtime = MPI_WTIME()
  
  allocate( nconf_low(num_procs), nconf_high(num_procs) )
  
  ket_confs = 0 
  bra_confs = 0 
  et_4 = 0.d0 
  et_5 = 0.d0 
  et_4tmp = 0.d0
  et_4tmp2  = 0.d0
  do channel   = 1, channels%number_t3_confs 
     nx = channels%t3_quantum_numbers(channel*4)
     ny = channels%t3_quantum_numbers(channel*4-1)
     nz = channels%t3_quantum_numbers(channel*4-2)
     tz = channels%t3_quantum_numbers(channel*4-3)
     
     
     
     nconf_low = 0
     nconf_high = 0
     !
     ! find number of configs on each process
     !
     bra_confs = size(  lookup_t3_configs(2,channel)%ival, 2)
     nconfs_tot = bra_confs
     number_mtxel_iam = floor(real(nconfs_tot/num_procs ))
     !if ( iam == 0 ) write(6,*) number_mtxel_iam 
     diff = bra_confs - num_procs * number_mtxel_iam 
     
     !if ( iam == 0 ) write(6,*) diff
     !if ( iam == 0 ) write(6,*) 'Total number of mtx-elements', bra_confs
     !if ( iam == 0 ) write(6,*) 'Total number of mtx-elements', number_mtxel_iam*num_procs + diff
     
     if ( iam == num_procs - 1 ) nconfs = number_mtxel_iam + diff
     if ( iam /= num_procs - 1 ) nconfs = number_mtxel_iam  
           
     n = 1 
     do i = 1, num_procs 
        if ( i < num_procs ) then
           nconf_low(i) = n 
           nconf_high(i) = n + number_mtxel_iam - 1
        else
           nconf_low(i) = n 
           nconf_high(i) = n + number_mtxel_iam - 1 + diff
        end if
        n = n + number_mtxel_iam 
     end do
     
     !do i = 1, num_procs
     !   if ( iam == 0 ) write(6,*) 'iam',i,'lower/upper', nconf_low(i), nconf_high(i)
     !end do
     
     
     if ( iam == 0 ) write(6,*) channel, size(  lookup_t3_configs(1,channel)%ival, 2) , size(  lookup_t3_configs(1,channel)%ival, 1) 
     
     tmp_sum1 =dcmplx(0.d0,0.d0) 
     tmp_sum2 =dcmplx(0.d0,0.d0) 
     !$omp parallel default(shared) private(ket,a,b,c,bra,i,j,k,spen, & 
     !$omp sum1,sum2,sum3,sum4,sum5,sum6,sum7,sum8,sum9,tc_term,td_term)
     !$omp do schedule(dynamic), reduction(+:tmp_sum1,tmp_sum2) 
     do ket = nconf_low(iam+1), nconf_high(iam+1) ! 1, size(  lookup_t3_configs(2,channel)%ival, 2) ! 
        a = lookup_t3_configs(2,channel)%ival(1,ket)
        b = lookup_t3_configs(2,channel)%ival(2,ket) 
        c = lookup_t3_configs(2,channel)%ival(3,ket) 
        
        do bra = 1, size(  lookup_t3_configs(1,channel)%ival, 2) 
           i = lookup_t3_configs(1,channel)%ival(1,bra) 
           j = lookup_t3_configs(1,channel)%ival(2,bra) 
           k = lookup_t3_configs(1,channel)%ival(3,bra) 
           
           
           spen = ( fock_mtx(i,i) + fock_mtx(j,j)  + fock_mtx(k,k) & 
                - fock_mtx(a,a)-fock_mtx(b,b) - fock_mtx(c,c) )
                       
           tc_term = 0.d0
           td_term = 0.d0 
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
           
           tc_term = ( sum1 - sum2 - sum3  - sum4 - sum5 + sum6 + sum7 + sum8 + sum9 )
           
           sum1 = 0.d0; sum2 = 0.d0; sum3 = 0.d0
           sum4 = 0.d0; sum5 = 0.d0; sum6 = 0.d0
           sum7 = 0.d0; sum8 = 0.d0; sum9 = 0.d0
           ! 
           ! P(c/ab)P(i/jk) <bc|V|dk><ad|t|ij>
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
           if ( tnf_approx == 3 ) then 
              td_term = tc_term + chiral_3nf_asym(a,b,c,i,j,k)
              
              sum1 = 0.d0; sum2 = 0.d0; sum3 = 0.d0
              sum4 = 0.d0; sum5 = 0.d0; sum6 = 0.d0
              sum7 = 0.d0; sum8 = 0.d0; sum9 = 0.d0
              !  1    ok 
              call three_body_v3t2pp(a,b,c,i,j,k,sum1)
              !  (ij)    ok 
              call three_body_v3t2pp(a,b,c,j,i,k,sum2)
              !  (ik)    ok 
              call three_body_v3t2pp(a,b,c,k,j,i,sum3)
              td_term = td_term + sum1 - sum2 - sum3
              
              sum1 = 0.d0; sum2 = 0.d0; sum3 = 0.d0
              sum4 = 0.d0; sum5 = 0.d0; sum6 = 0.d0
              sum7 = 0.d0; sum8 = 0.d0; sum9 = 0.d0
              !  1    ok 
              call three_body_v3t2hh(a,b,c,i,j,k,sum1)
              !  (ab)    ok 
              call three_body_v3t2hh(b,a,c,i,j,k,sum2)
              !  (ac)    ok 
              call three_body_v3t2hh(c,b,a,i,j,k,sum3)
              td_term = td_term + sum1 - sum2 - sum3
              
              sum1 = 0.d0; sum2 = 0.d0; sum3 = 0.d0
              sum4 = 0.d0; sum5 = 0.d0; sum6 = 0.d0
              sum7 = 0.d0; sum8 = 0.d0; sum9 = 0.d0
              !  1    ok 
              call three_body_v3t2ph(a,b,c,i,j,k,sum1)
              !  (ik)    ok 
              call three_body_v3t2ph(a,b,c,k,j,i,sum2)
              !  (jk)    ok 
              call three_body_v3t2ph(a,b,c,i,k,i,sum3)
              !  (ac)    ok 
              call three_body_v3t2ph(c,b,a,i,j,k,sum4)
              !  (bc)    ok 
              call three_body_v3t2ph(a,c,b,i,j,k,sum5)
              !  (ik)(ac)    ok 
              call three_body_v3t2ph(c,b,a,k,j,i,sum6)
              !  (ik)(bc)    ok 
              call three_body_v3t2ph(a,c,b,k,j,i,sum7)
              !  (jk)(ac)    ok 
              call three_body_v3t2ph(c,b,a,i,k,i,sum8)
              !  (jk)(bc)    ok 
              call three_body_v3t2ph(a,c,b,i,k,i,sum9)
              
              td_term = td_term + ( sum1 - sum2 - sum3  - sum4 - sum5 + sum6 + sum7 + sum8 + sum9 )
              
           end if
           tmp_sum1 = tmp_sum1 + conjg(tc_term)*(tc_term) / spen 
           tmp_sum2 = tmp_sum2 + conjg(td_term)*(td_term) / spen 
           
        end DO
     end DO
     !$omp end do
     !$omp end parallel
     
     et_4tmp  = et_4tmp  + tmp_sum1
     et_4tmp2  = et_4tmp2  + tmp_sum2
  end do
  
  
  
  call mpi_allreduce(et_4tmp,et_4,1,mpi_complex16,mpi_sum, &
       mpi_comm_world,ierror)
  call mpi_allreduce(et_4tmp2,et_5,1,mpi_complex16,mpi_sum, &
       mpi_comm_world,ierror)
   
  endwtime = MPI_WTIME()
  IF ( IAM == 0 ) WRITE(6,*) 'EXECUTION TIME FOR CCD(T)', endwtime-startwtime
  !if ( iam == 0 ) write(6,'(a,1x,2(F10.5,2x),a,1x,4(F10.5,3x))') 'Kf, Dens', kf, dens, 'HF, MBPT2, CCD, CCD(T)', e0, mbpt2, e0 + eccsd, e0 + eccsd + et_4
  if ( iam == 0 ) write(6,'(a,1x,2(F16.5,2x),a,1x,2(F16.5,3x))') 'Kf, Dens', & 
       kf, dens, 'CCSD(T)', real(et_4/below_ef), real(et_5/below_ef)
  if ( tnf_approx == 3 ) eccsdt = eccsd + et_5
  if ( tnf_approx  < 3 ) eccsdt = eccsd + et_4
  

end subroutine pert_triples


subroutine pert_triples2
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
  INTEGER :: i,j,k,l,a,b,c,d, number_channels, nx,ny,nz,sz,tz,  sumx, sumy, sumz,k1,k2,k3
  INTEGER :: bra_confs, ket_confs, bra,ket, ii, channel
  INTEGER :: sz_min, sz_max, tz_min, tz_max, sumtz, sumsz
  complex*16  :: spen, et_4tmp,et_4tmp2, et_4, et_4t, tmp_sum1, tmp_sum2, et_5,  & 
       est_5, sum1,sum2,sum3,sum4,sum5,sum6, sum7,sum8,sum9,tc_term, td_term,lc_term
  integer, allocatable :: nconf_low(:), nconf_high(:)
  integer :: nconfs, nconfs_tot, number_mtxel_iam, diff, n
  real*8  ::  startwtime , endwtime
  complex*16 :: diag1 
  integer :: bra_min, bra_max, ket_min, ket_max, local_ch_id
  
  
  
  call setup_proc_t3_mappings
  
  call mpi_barrier(mpi_comm_world,ierror)
  startwtime = MPI_WTIME()
  ket_confs = 0 
  bra_confs = 0 
  et_4 = 0.d0 
  et_5 = 0.d0 
  et_4tmp = 0.d0
  et_4tmp2  = 0.d0
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
     
     tmp_sum1 =dcmplx(0.d0,0.d0) 
     tmp_sum2 =dcmplx(0.d0,0.d0) 
     !$omp parallel default(shared) private(ket,a,b,c,bra,i,j,k,spen, & 
     !$omp sum1,sum2,sum3,sum4,sum5,sum6,sum7,sum8,sum9,tc_term,td_term)
     !$omp do schedule(dynamic), reduction(+:tmp_sum1,tmp_sum2) 
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
                       
           tc_term = 0.d0
           td_term = 0.d0 
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
           
           tc_term = ( sum1 - sum2 - sum3  - sum4 - sum5 + sum6 + sum7 + sum8 + sum9 )
           
           sum1 = 0.d0; sum2 = 0.d0; sum3 = 0.d0
           sum4 = 0.d0; sum5 = 0.d0; sum6 = 0.d0
           sum7 = 0.d0; sum8 = 0.d0; sum9 = 0.d0
           ! 
           ! P(c/ab)P(i/jk) <bc|V|dk><ad|t|ij>
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
           if ( tnf_approx == 3 ) then 
              td_term = tc_term + chiral_3nf_asym(a,b,c,i,j,k)
              
!!$              sum1 = 0.d0; sum2 = 0.d0; sum3 = 0.d0
!!$              sum4 = 0.d0; sum5 = 0.d0; sum6 = 0.d0
!!$              sum7 = 0.d0; sum8 = 0.d0; sum9 = 0.d0
!!$              !  1    ok 
!!$              !call three_body_v3t2pp(a,b,c,i,j,k,sum1)
!!$              !  (ij)    ok 
!!$              !call three_body_v3t2pp(a,b,c,j,i,k,sum2)
!!$              !  (ik)    ok 
!!$              !call three_body_v3t2pp(a,b,c,k,j,i,sum3)
!!$              td_term = td_term + sum1 - sum2 - sum3
!!$              
!!$              sum1 = 0.d0; sum2 = 0.d0; sum3 = 0.d0
!!$              sum4 = 0.d0; sum5 = 0.d0; sum6 = 0.d0
!!$              sum7 = 0.d0; sum8 = 0.d0; sum9 = 0.d0
!!$              !  1    ok 
!!$              !call three_body_v3t2hh(a,b,c,i,j,k,sum1)
!!$              !  (ab)    ok 
!!$              !call three_body_v3t2hh(b,a,c,i,j,k,sum2)
!!$              !  (ac)    ok 
!!$              !call three_body_v3t2hh(c,b,a,i,j,k,sum3)
!!$              td_term = td_term + sum1 - sum2 - sum3
!!$              
!!$              sum1 = 0.d0; sum2 = 0.d0; sum3 = 0.d0
!!$              sum4 = 0.d0; sum5 = 0.d0; sum6 = 0.d0
!!$              sum7 = 0.d0; sum8 = 0.d0; sum9 = 0.d0
!!$              !  1    ok 
!!$              call three_body_v3t2ph(a,b,c,i,j,k,sum1)
!!$              !  (ik)    ok 
!!$              call three_body_v3t2ph(a,b,c,k,j,i,sum2)
!!$              !  (jk)    ok 
!!$              call three_body_v3t2ph(a,b,c,i,k,i,sum3)
!!$              !  (ac)    ok 
!!$              call three_body_v3t2ph(c,b,a,i,j,k,sum4)
!!$              !  (bc)    ok 
!!$              call three_body_v3t2ph(a,c,b,i,j,k,sum5)
!!$              !  (ik)(ac)    ok 
!!$              call three_body_v3t2ph(c,b,a,k,j,i,sum6)
!!$              !  (ik)(bc)    ok 
!!$              call three_body_v3t2ph(a,c,b,k,j,i,sum7)
!!$              !  (jk)(ac)    ok 
!!$              call three_body_v3t2ph(c,b,a,i,k,i,sum8)
!!$              !  (jk)(bc)    ok 
!!$              call three_body_v3t2ph(a,c,b,i,k,i,sum9)
!!$              
!!$              td_term = td_term + ( sum1 - sum2 - sum3  - sum4 - sum5 + sum6 + sum7 + sum8 + sum9 )
              
           end if
           tmp_sum1 = tmp_sum1 + conjg(tc_term)*(tc_term) / spen 
           tmp_sum2 = tmp_sum2 + conjg(td_term)*(td_term) / spen 
           
        end DO
     end DO
     !$omp end do
     !$omp end parallel

     et_4tmp  = et_4tmp  + tmp_sum1
     et_4tmp2  = et_4tmp2  + tmp_sum2
  end do
  
  
  
  call mpi_allreduce(et_4tmp,et_4,1,mpi_complex16,mpi_sum, &
       mpi_comm_world,ierror)
  call mpi_allreduce(et_4tmp2,et_5,1,mpi_complex16,mpi_sum, &
       mpi_comm_world,ierror)
  
  call mpi_barrier(mpi_comm_world,ierror)
  endwtime = MPI_WTIME()
  IF ( IAM == 0 ) WRITE(6,*) 'EXECUTION TIME FOR CCD(T)', endwtime-startwtime
  !if ( iam == 0 ) write(6,'(a,1x,2(F10.5,2x),a,1x,4(F10.5,3x))') 'Kf, Dens', kf, dens, 'HF, MBPT2, CCD, CCD(T)', e0, mbpt2, e0 + eccsd, e0 + eccsd + et_4
  if ( iam == 0 ) write(6,'(a,1x,2(F16.5,2x),a,1x,2(F16.5,3x))') 'Kf, & 
       Dens', kf, dens, 'CCSD(T)', real(et_4/below_ef), real(et_5/below_ef)
  if ( tnf_approx == 3 ) eccsdt = eccsd + et_5
  if ( tnf_approx  < 3 ) eccsdt = eccsd + et_4
  

end subroutine pert_triples2

