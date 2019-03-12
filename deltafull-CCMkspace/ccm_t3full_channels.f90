!
!
!
SUBROUTINE setup_t3_channel_structures
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
  INTEGER, ALLOCATABLE :: numchannels1(:)
  real*8 :: memory 
  INTEGER :: number_channels, k1,k2,k3,k4,k5
  
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
  
  
  ALLOCATE( lookup_t3_configs(1:2,numchannels))
  memory = 0.d0 
  !     loop over isospin projection
  DO tz=tz_min3,tz_max3,2 
     !     loop over total momentum 
     DO nx=nx_min3, nx_max3 
        !     loop over total momentum 
        DO ny=ny_min3, ny_max3 
           !     loop over total momentum 
           DO nz=nz_min3, nz_max3 
                 
              channel = locate_t3channel(tz,nx,ny,nz) 
              if ( channel == 0 ) cycle 
              ij_confs = numstates(1,tz, nx, ny, nz) 
              ab_confs = numstates(2,tz, nx, ny, nz) 
              memory = memory + dble(ij_confs+ab_confs)*3.*4./1.e9 
              ALLOCATE( lookup_t3_configs(1,channel)%ival(3,ij_confs) )
              ALLOCATE( lookup_t3_configs(2,channel)%ival(3,ab_confs) )
              
           end DO
        end do
     end DO
  end DO
  
  if ( iam == 0 ) write(6,*) 'Total memory for lookup t3-confs', memory, 'Gbyte'
    
  allocate( numchannels1( maxval( locate_t3channel(:,:,:,:) ) ) )
  
  !
  ! ppp configurations
  !
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
           
           
           if ( locate_t3channel(tz, nx, ny, nz) /= 0 ) then
              channel = locate_t3channel(tz, nx, ny, nz) 
              numchannels1(channel) = numchannels1(channel) + 1 
              ichan = numchannels1(channel)
              
              !ppp_hhhppp%ival3(a,b,c) = ichan
              !hh_hhhp%ival(a,b) = ichan
              lookup_t3_configs(2,channel)%ival(1,ichan) = a 
              lookup_t3_configs(2,channel)%ival(2,ichan) = b 
              lookup_t3_configs(2,channel)%ival(3,ichan) = c 
              
           end if
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
           
           
           if ( locate_t3channel(tz, nx, ny, nz) /= 0 ) then
              channel = locate_t3channel(tz, nx, ny, nz) 
              numchannels1(channel) = numchannels1(channel) + 1 
              ichan = numchannels1(channel)
              
              !hhh_hhhppp%ival3(a,b,c) = ichan
              lookup_t3_configs(1,channel)%ival(1,ichan) = a 
              lookup_t3_configs(1,channel)%ival(2,ichan) = b 
              lookup_t3_configs(1,channel)%ival(3,ichan) = c 
              
           end if
        end do
     end do
  end do
  
  memory = 0.d0 
  number_channels = 0 
  number_channels = maxval(locate_t3channel(:,:,:,:) ) 
  channels%number_t3_confs  = number_channels
  ALLOCATE( channels%t3_quantum_numbers(4*number_channels) )
  memory = memory + dble(number_channels)*5.*4./1.e9 
  

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
  

  
end SUBROUTINE setup_t3_channel_structures



!
!
!
SUBROUTINE setup_t3
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
  INTEGER :: number_channels, k1,k2,k3,k4,k5
  INTEGER :: local_ch_id, dim1, dim2, bra_min, bra_max, ket_min, ket_max, bra, ket, bra2, ket2
  COMPLEX*16 :: spen, v3nf 
  complex*16, allocatable :: amat(:,:)
  real*8 :: mem 
  
  call mpi_barrier(mpi_comm_world, ierror) 
  
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
     
     
     !write(6,*) iam , channel, size(lookup_t3_configs(2,channel)%ival,2), size(lookup_t3_configs(1,channel)%ival,2)
!!$     dim1 = size(lookup_t3_configs(2,channel)%ival,2)
!!$     dim2 = size(lookup_t3_configs(1,channel)%ival,2)
!!$     allocate( t3_ccm(channel)%val(dim1, dim2))
!!$     allocate( v3nf_ppphhh(channel)%val(dim1, dim2))
!!$     t3_ccm(channel)%val = 0.d0
!!$     v3nf_ppphhh(channel)%val = 0.d0 
     
     dim1 = bra_max-bra_min+1 
     dim2 = ket_max-ket_min+1
     allocate( t3_ccm(channel)%val(bra_min:bra_max, ket_min:ket_max))
     allocate( v3nf_ppphhh(channel)%val(bra_min:bra_max, ket_min:ket_max))
     t3_ccm(channel)%val(bra_min:bra_max, ket_min:ket_max) = 0.d0
     v3nf_ppphhh(channel)%val(bra_min:bra_max, ket_min:ket_max) = 0.d0 
     
     mem = mem + dble(dim1*dim2)*16.*2./1.e9  
     
  end do
  call mpi_barrier(mpi_comm_world, ierror) 
  !write(61,*) 'iam', iam, 'Total memory for t3', mem, 'Gbyte' 
  if ( iam == 0 ) write(6,*) 'iam', iam, 'Total memory for v3nf and t3', mem, 'Gbyte' 
  
  call mpi_barrier(mpi_comm_world, ierror) 
  
  
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
     
     
     ppp_hhhppp%ival3 = 0 
     do bra = 1, size(lookup_t3_configs(2,channel)%ival,2)
        a = lookup_t3_configs(2,channel)%ival(1,bra)
        b = lookup_t3_configs(2,channel)%ival(2,bra) 
        c = lookup_t3_configs(2,channel)%ival(3,bra) 
        
        ppp_hhhppp%ival3(a,b,(all_orbit%szp(c)+1)/2) = bra
     end do
     
     hhh_hhhppp%ival3 = 0 
     do ket = 1, size(  lookup_t3_configs(1,channel)%ival, 2) 
        i = lookup_t3_configs(1,channel)%ival(1,ket) 
        j = lookup_t3_configs(1,channel)%ival(2,ket) 
        k = lookup_t3_configs(1,channel)%ival(3,ket) 
        hhh_hhhppp%ival3(i,j,k) = ket
     end do
     
!!$     dim1 = size(lookup_t3_configs(2,channel)%ival,2)
!!$     dim2 = size(lookup_t3_configs(1,channel)%ival,2)
!!$     allocate( amat(dim1,dim2) ) 
!!$     amat = 0.d0 

     !$omp parallel default(shared) private(bra,a,b,c,ket,i,j,k,v3nf,bra2,ket2) 
     !$omp do schedule(dynamic)
     do bra = bra_min, bra_max  ! 1, size(lookup_t3_configs(2,channel)%ival,2) 
        a = lookup_t3_configs(2,channel)%ival(1,bra)
        b = lookup_t3_configs(2,channel)%ival(2,bra) 
        c = lookup_t3_configs(2,channel)%ival(3,bra) 
        
        !if ( a > b ) cycle 
        !if ( b > c ) cycle 
        !if ( a > c ) cycle 
        do ket = 1, size(  lookup_t3_configs(1,channel)%ival, 2) 
           i = lookup_t3_configs(1,channel)%ival(1,ket) 
           j = lookup_t3_configs(1,channel)%ival(2,ket) 
           k = lookup_t3_configs(1,channel)%ival(3,ket) 
           
           if (  i > j ) cycle    
           if (  j > k ) cycle    
           if (  i > k ) cycle    
           v3nf = chiral_3nf_asym(a,b,c,i,j,k) 
           v3nf_ppphhh(channel)%val(bra,ket) = v3nf 
           !amat(bra,ket) = v3nf  
           
           
        end do
     end do
     !$omp end do
     !$omp end parallel
     
     
     do bra =  bra_min, bra_max 
        a = lookup_t3_configs(2,channel)%ival(1,bra)
        b = lookup_t3_configs(2,channel)%ival(2,bra) 
        c = lookup_t3_configs(2,channel)%ival(3,bra) 
        
        do ket = 1, size(  lookup_t3_configs(1,channel)%ival, 2) 
           i = lookup_t3_configs(1,channel)%ival(1,ket) 
           j = lookup_t3_configs(1,channel)%ival(2,ket) 
           k = lookup_t3_configs(1,channel)%ival(3,ket) 
        
           if (  i > j ) cycle    
           if (  j > k ) cycle    
           if (  i > k ) cycle    
                   
           v3nf = v3nf_ppphhh(channel)%val(bra,ket) 
           ! (ij) 
           ket2 = hhh_hhhppp%ival3(j,i,k)
           v3nf_ppphhh(channel)%val(bra,ket2) = v3nf_ppphhh(channel)%val(bra,ket2) - v3nf 
           ! (ik) 
           ket2 = hhh_hhhppp%ival3(k,j,i)
           v3nf_ppphhh(channel)%val(bra,ket2) = v3nf_ppphhh(channel)%val(bra,ket2) - v3nf 
           ! (jk) 
           ket2 = hhh_hhhppp%ival3(i,k,j)
           v3nf_ppphhh(channel)%val(bra,ket2) = v3nf_ppphhh(channel)%val(bra,ket2) - v3nf 
           ! (jki)
           ket2 = hhh_hhhppp%ival3(j,k,i)
           v3nf_ppphhh(channel)%val(bra,ket2) = v3nf_ppphhh(channel)%val(bra,ket2) + v3nf 
           ! (kij)
           ket2 = hhh_hhhppp%ival3(k,i,j)
           v3nf_ppphhh(channel)%val(bra,ket2) = v3nf_ppphhh(channel)%val(bra,ket2) + v3nf 
           
        end do
     end do

!!$     
!!$     do bra =  1, size(lookup_t3_configs(2,channel)%ival,2) !bra_min, bra_max 
!!$        a = lookup_t3_configs(2,channel)%ival(1,bra)
!!$        b = lookup_t3_configs(2,channel)%ival(2,bra) 
!!$        c = lookup_t3_configs(2,channel)%ival(3,bra) 
!!$        
!!$        if ( a > b ) cycle 
!!$        if ( b > c ) cycle 
!!$        if ( a > c ) cycle 
!!$        do ket = 1, size(  lookup_t3_configs(1,channel)%ival, 2) 
!!$           i = lookup_t3_configs(1,channel)%ival(1,ket) 
!!$           j = lookup_t3_configs(1,channel)%ival(2,ket) 
!!$           k = lookup_t3_configs(1,channel)%ival(3,ket) 
!!$           
!!$           if (  i > j ) cycle    
!!$           
!!$           v3nf = amat(bra,ket)
!!$           ! (ab) 
!!$           bra2 = ppp_hhhppp%ival3(b,a,(all_orbit%szp(c)+1)/2)
!!$           amat(bra2,ket) = amat(bra2,ket) - v3nf 
!!$           ! (ij) 
!!$           ket2 = hhh_hhhppp%ival3(j,i,k)
!!$           amat(bra,ket2) = amat(bra,ket2) - v3nf 
!!$           ! (ab)(ij) 
!!$           bra2 = ppp_hhhppp%ival3(b,a,(all_orbit%szp(c)+1)/2)
!!$           ket2 = hhh_hhhppp%ival3(j,i,k)
!!$           amat(bra2,ket2) = amat(bra2,ket2) + v3nf 
!!$           ! (ac) 
!!$           bra2 = ppp_hhhppp%ival3(c,b,(all_orbit%szp(a)+1)/2)
!!$           amat(bra2,ket) = amat(bra2,ket) - v3nf 
!!$           ! (bc) 
!!$           bra2 = ppp_hhhppp%ival3(a,c,(all_orbit%szp(b)+1)/2)
!!$           amat(bra2,ket) = amat(bra2,ket) - v3nf 
!!$           ! (ac)(ij) 
!!$           bra2 = ppp_hhhppp%ival3(c,b,(all_orbit%szp(a)+1)/2)
!!$           ket2 = hhh_hhhppp%ival3(j,i,k)
!!$           amat(bra2,ket2) = amat(bra2,ket2) + v3nf 
!!$           ! (bc)(ij)
!!$           bra2 = ppp_hhhppp%ival3(a,c,(all_orbit%szp(b)+1)/2)
!!$           ket2 = hhh_hhhppp%ival3(j,i,k)
!!$           amat(bra2,ket2) = amat(bra2,ket2) + v3nf 
!!$           ! (bca) 
!!$           bra2 = ppp_hhhppp%ival3(b,c,(all_orbit%szp(a)+1)/2)
!!$           amat(bra2,ket) = amat(bra2,ket) + v3nf 
!!$           ! (cab) 
!!$           bra2 = ppp_hhhppp%ival3(c,a,(all_orbit%szp(b)+1)/2)
!!$           amat(bra2,ket) = amat(bra2,ket) + v3nf 
!!$           ! (bca)(ij) 
!!$           bra2 = ppp_hhhppp%ival3(b,c,(all_orbit%szp(a)+1)/2)
!!$           amat(bra2,ket2) = amat(bra2,ket2) - v3nf 
!!$           ! (cab)(ij) 
!!$           bra2 = ppp_hhhppp%ival3(c,a,(all_orbit%szp(b)+1)/2)
!!$           amat(bra2,ket2) = amat(bra2,ket2) - v3nf 
!!$        end do
!!$     end do
!!$
!!$     v3nf_ppphhh(channel)%val = amat
!!$     deallocate (amat)
     
     
  end do
  
  call mpi_barrier(mpi_comm_world, ierror) 
  if ( iam == 0 ) write(6,*) 'setting up t3'
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
     

     t3_ccm(channel)%val = 0.d0
     do bra = bra_min, bra_max  ! 1, size(  lookup_t3_configs(2,channel)%ival, 2)  !
        a = lookup_t3_configs(2,channel)%ival(1,bra)
        b = lookup_t3_configs(2,channel)%ival(2,bra) 
        c = lookup_t3_configs(2,channel)%ival(3,bra) 
                
        do ket = 1, size(  lookup_t3_configs(1,channel)%ival, 2) 
           i = lookup_t3_configs(1,channel)%ival(1,ket) 
           j = lookup_t3_configs(1,channel)%ival(2,ket) 
           k = lookup_t3_configs(1,channel)%ival(3,ket) 
           
           spen = ( fock_mtx(i,i) + fock_mtx(j,j)  + fock_mtx(k,k) & 
                - fock_mtx(a,a)-fock_mtx(b,b) - fock_mtx(c,c) )
           
           t3_ccm(channel)%val(bra,ket) = v3nf_ppphhh(channel)%val(bra,ket)/spen 
        end do
     end do
     
  end do
  call mpi_barrier(mpi_comm_world, ierror) 
  if ( iam == 0 ) write(6,*) 'setting up t3 done'
  
end SUBROUTINE setup_t3


SUBROUTINE t3_energy(ener3)
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
  COMPLEX*16, INTENT(INOUT) :: ENER3
  INTEGER :: i,j,k, sumx, sumy, sumz, a,b,c,nxa,nxb,nxc,nya,nyb,nyc,nza,nzb,nzc,tza,tzb,tzc,sza,szb,szc,tz,sz,nx,ny,nz
  INTEGER :: numchannels, ij_confs, ab_confs, channel, sz_min, sz_max, tz_min, tz_max, sumtz, sumsz, ichan
  INTEGER :: number_channels, k1,k2,k3,k4,k5
  INTEGER :: local_ch_id, dim1, dim2, bra_min, bra_max, ket_min, ket_max, bra, ket
  COMPLEX*16 :: spen, sum1, t3_ener, tmp_sum1
  integer :: cc1, cc2
  
  sum1 = 0.d0 
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
     
     tmp_sum1 = dcmplx(0.d0)
     !$omp parallel default(shared) private(bra,a,b,c,ket,i,j,k,spen)
     !$omp do schedule(dynamic), reduction(+:tmp_sum1)
     do bra = bra_min, bra_max 
        a = lookup_t3_configs(2,channel)%ival(1,bra)
        b = lookup_t3_configs(2,channel)%ival(2,bra) 
        c = lookup_t3_configs(2,channel)%ival(3,bra) 
     
        if ( a > b ) cycle
        if ( a > c ) cycle
        if ( b > c ) cycle 
        
        do ket = 1, size(  lookup_t3_configs(1,channel)%ival, 2) 
           i = lookup_t3_configs(1,channel)%ival(1,ket) 
           j = lookup_t3_configs(1,channel)%ival(2,ket) 
           k = lookup_t3_configs(1,channel)%ival(3,ket) 

           if ( i > j ) cycle
           if ( i > k ) cycle
           if ( j > k ) cycle 
           
!!$           nxc = all_orbit%nx(i) + all_orbit%nx(j) + all_orbit%nx(k) - all_orbit%nx(a) - all_orbit%nx(b)
!!$           nyc = all_orbit%ny(i) + all_orbit%ny(j) + all_orbit%ny(k) - all_orbit%ny(a) - all_orbit%ny(b)
!!$           nzc = all_orbit%nz(i) + all_orbit%nz(j) + all_orbit%nz(k) - all_orbit%nz(a) - all_orbit%nz(b)
!!$           tzc = all_orbit%itzp(i) + all_orbit%itzp(j) + all_orbit%itzp(k) - all_orbit%itzp(a) - all_orbit%itzp(b)
!!$           if ( abs(nxc) > nkxmax ) cycle 
!!$           if ( abs(nyc) > nkymax ) cycle 
!!$           if ( abs(nzc) > nkzmax ) cycle 
!!$           if ( abs(tzc) > 1 ) cycle 
!!$           cc1 = orbitof(nxc, nyc, nzc, -1, tzc)
!!$           cc2 = orbitof(nxc, nyc, nzc,  1, tzc)
!!$           
!!$           write(6,*) cc1,cc2,c, all_orbit%szp(c)   
!!$           

           spen = ( fock_mtx(i,i) + fock_mtx(j,j)  + fock_mtx(k,k) & 
                - fock_mtx(a,a)-fock_mtx(b,b) - fock_mtx(c,c) )
           
           tmp_sum1 = tmp_sum1 + t3_ccm(channel)%val(bra,ket)*conjg(v3nf_ppphhh(channel)%val(bra,ket)) 
           !sum1 = sum1 + conjg(chiral_3nf_asym(a,b,c,i,j,k))*chiral_3nf_asym(a,b,c,i,j,k)/spen 
           
        end do
     end do
     !$omp end do
     !$omp end parallel
     sum1 = sum1 + tmp_sum1
     
  end do
  
  call mpi_allreduce(sum1,t3_ener,1,mpi_complex16,mpi_sum, &
       mpi_comm_world,ierror)
  
  ener3 = ener3 + t3_ener 
  if ( iam == 0 ) write(6,*) 'CCDT-1 energy', t3_ener 
  
  
end SUBROUTINE t3_energy
