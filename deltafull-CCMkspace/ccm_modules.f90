!             Program block heff-modules.f   
!
!             Author:   Morten Hjorth-Jensen
!             ADDRESS:  Dept. Physics, University Oslo, N-0316 OSLO
!             E-MAIL:   morten.hjorth-jensen@fys.uio.no
!             LANGUAGE: F90/F95
!             LAST UPGRADE : April 2005
!             This program block contains the definiton of
!             all modules used by the various program blocks.
!    This module contains all constants and declarations 
!    of variables read in by the function read_data. These
!    variables are used by many functions.

MODULE parallel
  
  include 'mpif.h'
  INTEGER, PUBLIC   :: ierror,iam,num_procs,master,nc_break
  
  INTEGER :: mystatus(MPI_STATUS_SIZE)

  TYPE, PUBLIC :: topology
    INTEGER :: my_size, my_start, my_stop
    INTEGER, DIMENSION(:), POINTER :: all_starts, all_stops, all_sizes
  END TYPE topology
  TYPE(topology) :: tfull, tmin
  
  INTEGER :: DSIZE =8 ! 16 for complex
CONTAINS
  SUBROUTINE setup_parallel_lookup(this, n)
        TYPE(topology) :: this
        INTEGER, INTENT(IN) :: n
        INTEGER :: overflow, proc


        overflow = MOD(n,num_procs)
        ALLOCATE(this%all_starts(num_procs), this%all_stops(num_procs), &
                this%all_sizes(num_procs))

        !IF (iam == master .and. debug > 1)  WRITE(*,*) 'Num procs: ',  num_procs
        !IF (iam == master.AND. debug > 1) WRITE(*,*) 'Matrix dimension: ', n

        DO proc = 0, num_procs - 1
            IF (proc < overflow) THEN
                this%all_sizes(proc+1) = n/num_procs + 1
                this%all_starts(proc+1) = proc*this%all_sizes(proc+1)+ 1
            ELSE
                this%all_sizes(proc+1) = n/num_procs
                this%all_starts(proc+1) = overflow*(this%all_sizes(proc+1)+ 1)
                this%all_starts(proc+1) = this%all_starts(proc+1) + (proc-overflow)*this%all_sizes(proc+1) +1
            ENDIF
            this%all_stops(proc+1) = this%all_starts(proc+1) + this%all_sizes(proc+1) -1
        ENDDO

        this%my_size = this%all_sizes(iam + 1)
        this%my_start = this%all_starts(iam + 1)
        this%my_stop = this%all_stops(iam + 1)
        !WRITE(*,*) 'iam', iam, my_size, my_start, my_stop

        DO proc = 0, num_procs - 1
        !    IF (iam == master .and. debug > 1) WRITE(*,*) 'Local:', proc, this%all_sizes(proc+1), this%all_starts(proc+1), this%all_stops(proc+1)
        ENDDO

    END SUBROUTINE setup_parallel_lookup

    SUBROUTINE takedown_parallel_lookup(this)
        TYPE(topology) :: this
        IF(ASSOCIATED(this%all_starts)) DEALLOCATE(this%all_starts)
        IF(ASSOCIATED(this%all_stops)) DEALLOCATE(this%all_stops)
        IF(ASSOCIATED(this%all_sizes)) DEALLOCATE(this%all_sizes)
    END SUBROUTINE takedown_parallel_lookup


END MODULE parallel


MODULE constants
  INTEGER,  PARAMETER :: dp = KIND(1.0D0)
  INTEGER, PARAMETER :: dpc = KIND((1.0D0,1.0D0))
  ! min and max isospin projection
  INTEGER, PUBLIC :: itzmin, itzmax
  ! min and max total two-body angular momentum in lab frame
  INTEGER, PUBLIC :: j_lab_min, j_lab_max
  ! min and max total two-body angular momentum in Rel-CoM frame
  INTEGER, PUBLIC :: nkxmax, nkymax, nkzmax
  INTEGER, PUBLIC :: occ_protons, occ_neutrons, ntot
  real*8, public :: com_beta,com_switch, j2_beta, volume, lx, ly, lz
  REAL*8, public :: hcom_value, kin, kf, dens, cutoff_3nf_reg, lambda_delta 
  COMPLEX*16 :: E0,  eccsd, eccsdt, mbpt2, total_energy
  INTEGER, PUBLIC :: mass_nucleus, tot_orbs, below_ef, above_ef
  REAL(DP), PUBLIC ::  oscl, hbar_omega, kinfactor
  REAL(DP) , PARAMETER, PUBLIC :: p_mass = 939.565 !938.926_dp
  REAL(DP), PARAMETER, PUBLIC :: theta_rot = 0.0_dp! 0.125_dp
  REAL(DP), PARAMETER, PUBLIC :: hbarc = 197.326968_dp
  REAL(DP), PARAMETER, PUBLIC :: hb2ip = hbarc*hbarc/p_mass
  REAL(DP), PUBLIC, PARAMETER :: pi = 3.14159265358979_dp  !3.141592741012573_dp
  REAL(DP), PUBLIC, PARAMETER :: pi_2 = 1.570796370506287_dp
  REAL(DP), PUBLIC, PARAMETER :: pi_4 = 0.7853981852531433_dp
  LOGICAL :: switch_density, tnf_switch, chiral_delta_flag
  CHARACTER(LEN=100) :: boundary_conditions, cc_approx
  INTEGER :: tnf_approx, nexp_3nf_nonlocal, delta_chiral_order
  
END MODULE constants


MODULE kspace 

  INTEGER, PUBLIC :: nx_min, nx_max, ny_min, ny_max, nz_min, nz_max 

end MODULE kspace
 
MODULE t3_constants 

  INTEGER, PUBLIC :: nx_min3, nx_max3, ny_min3, ny_max3, nz_min3, nz_max3
  INTEGER, PUBLIC :: sz_min3, sz_max3, tz_min3, tz_max3 
  INTEGER, ALLOCATABLE, public :: LOCATE_T3CHANNEL(:,:,:,:)
  
end MODULE t3_constants

!
!            
!     This module contains the angular momentun functions
!     and transformation coefficients when going from 
!     lab system  <--> cm system
!
MODULE ang_mom_functions
  REAL*8, PRIVATE :: f_mb(50),g_mb(50),w_mb(50)
  INTEGER, PRIVATE :: kh(200)
  REAL*8, PARAMETER, PRIVATE :: pi=3.141592654
  REAL*8, PRIVATE :: q(50,50), cn(0:51,0:51)

CONTAINS
  !
  !     factorials for 3j,6j and 9j symbols             
  !     for moshinsky trans brackets and for            
  !     vector brackets                                 
  !
  SUBROUTINE commons_to_angmom
    IMPLICIT NONE 
    INTEGER :: l, k, i, j
    REAL*8 :: a , sq_pi, fj, tfj, fk
    !    3j, 6j and 9j symbols
    kh=1
    kh(100) =0
    DO l=1,50
       q(l,1)=1.0d0
       q(l,l)=1.0d0
       kh(l+l+100)=0
    ENDDO
    DO l=2,49
       DO k=2,l
          q(l+1,k)=q(l,k-1)+q(l,k)
       ENDDO
    ENDDO
    !    Moshinsky brackets
    f_mb(1)=0.
    g_mb(1)=LOG(0.5D0)
    w_mb(1)=0.
    DO i=2,50
       a=i-1
       f_mb(i)=f_mb(i-1)+LOG(a)
       g_mb(i)=g_mb(i-1)+LOG(a+0.5D0)
       w_mb(i)=LOG(a+a+1.)
    ENDDO
    !    spherical harmonics
    cn=0.
    sq_pi=1./SQRT(2.*pi)
    DO j=0,51
       cn(0,j)=SQRT(0.5*(2.*j+1.))
    ENDDO
    DO j=1,51
       tfj=2.*j
       cn(j,j)=cn(j-1,j-1)*SQRT((tfj+1.)/tfj)
    ENDDO
    DO j=0,51
       fj=FLOAT(j)
       DO k=1,j-1
          fk=FLOAT(k)
          cn(k,j)=cn(k-1,j)*SQRT((fj+fk)*(fj-fk+1.))*0.5/fk
       ENDDO
    ENDDO
    cn=cn*sq_pi

  END SUBROUTINE commons_to_angmom
  !
  !     calculates 3j-symbols            
  !
  REAL*8 FUNCTION tjs(j_a,j_b,j_c,m_a,m_b,m_c)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: j_a,j_b,j_c,m_a,m_b,m_c
    INTEGER :: ja, jb, jc, mb, ma, mc, la, lb, lc, lt, ld, ja2, jb2, &
         jc2, i, k0, k1, k, ip
    REAL*8 :: x, fn, p

    tjs=0.
    ja=(j_a+m_a)/2+1
    ma=(j_a-m_a)/2+1
    jb=(j_b+m_b)/2+1
    mb=(j_b-m_b)/2+1
    jc=(j_c+m_c)/2+1
    mc=(j_c-m_c)/2+1
    la=(j_b+j_c-j_a)/2+1
    lb=(j_c+j_a-j_b)/2+1
    lc=(j_a+j_b-j_c)/2+1
    lt=(j_a+j_b+j_c)/2+1
    ld=MIN(ja,jb,jc,ma,mb,mc,la,lb,lc)
    IF(((m_a+m_b+m_c) <= 0).AND.(ld > 0)) THEN
       ja2=j_a+m_a
       jb2=j_b+m_b
       jc2=j_c+m_c
       i=ja2+jb2+jc2-ja2/2*2-jb2/2*2-jc2/2*2
       IF(i == 0) then 
          fn=q(ja+ma-1,lc)*q(jb+mb-1,lc)/(q(lt,jc+mc-1)*q(lt+1,2) &
               *q(ja+ma-1,ja)*q(jb+mb-1,jb)*q(jc+mc-1,jc))
          k0=MAX(0,lc-ja,lc-mb)+1
          k1=MIN(lc,ma,jb)
          x=0.
          DO k=k0,k1
             x=-x-q(lc,k)*q(lb,ma-k+1)*q(la,jb-k+1)
          ENDDO
          ip=k1+lb+jc
          p=1-2*(ip-ip/2*2)
          tjs=p*x*SQRT(fn)
       ENDIF
    ENDIF

  END FUNCTION tjs
  !
  !     calculates 6j-symbols            
  !
  REAL*8 FUNCTION sjs(j_a,j_b,j_c,l_a,l_b,l_c)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: j_a,j_b,j_c,l_a,l_b,l_c
    INTEGER :: ja,jb,jc,la,lb,lc,i,mt,ma,mb,mc,na,nb,nc,ka,&
         kb,kc,l,l0,l1
    REAL*8 :: x, fs, fss

    sjs=0.0d0
    ja=j_a + 1
    jb=j_b + 1
    jc=j_c + 1
    la=l_a + 1
    lb=l_b + 1
    lc=l_c + 1
    i=kh(ja+jb-jc+99)+kh(jb+jc-ja+99)+kh(jc+ja-jb+99)+kh(ja+lb-lc+99) &
         +kh(lb+lc-ja+99)+kh(lc+ja-lb+99)+kh(la+jb-lc+99)+kh(jb+lc-la+99) &
         +kh(lc+la-jb+99)+kh(la+lb-jc+99)+kh(lb+jc-la+99)+kh(jc+la-lb+99)
    IF(i <= 0) THEN
       mt=(j_a+j_b+j_c)/2 + 2
       ma=(j_a+l_b+l_c)/2+ 2
       mb=(l_a+j_b+l_c)/2+ 2
       mc=(l_a+l_b+j_c)/2+ 2
       na=mt-ja
       nb=mt-jb
       nc=mt-jc
       ka=ma-lc
       kb=mb-lc
       kc=mc-jc
       fss=q(mt,ja+1)*q(ja,nc)/(q(ma,ja+1)*q(ja,ka)*q(mb,la+1)* &
            q(la,kb)*q(mc,la+1)*q(la,kc))
       fs=SQRT(fss)/(l_a + 1.)
       l0=MAX(mt,ma,mb,mc)+1
       l1=MIN(ma+na,mb+nb,mc+nc)
       x=0.
       DO l=l0,l1
          x=-x+q(l-1,mt)*q(na,l-ma)*q(nb,l-mb)*q(nc,l-mc)
       ENDDO
       sjs=-(1+2*(l1/2*2-l1))*fs*x
    ENDIF

  END FUNCTION sjs
  !
  !     calculates ninej-symbols
  !       
  REAL*8 FUNCTION snj (ia,ib,ie,ic,id,if,ig,ih,it)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ia,ib,ie,ic,id,if,ig,ih,it
    INTEGER :: ja,jb,je,jc,jd,jf,jg,jh,jt,i,la,ld,ma,mc,na,nb,le,lf,&
         lg,me,mf,mg,ne,nf,ng,lx,mx,nx,jsi,jsf, js,is,lb, lc, &
         mb, ly, my,ny,l,l0,m0,n0,l1,m1,n1,m,n,ihx
    REAL*8 :: x, fn, fd, ps, fs, u, y, z, ud, p

    snj=0.
    ja=ia+1
    jb=ib+1
    jc=ic+1
    jd=id+1
    je=ie+1
    jf=IF+1
    jg=ig+1
    jh=ih+1
    jt=it+1
    i=kh(ja+jb-je+99)+kh(jb+je-ja+99)+kh(je+ja-jb+99)+kh(jc+jd-jf+99) &
         +kh(jd+jf-jc+99)+kh(jf+jc-jd+99)+kh(jg+jh-jt+99)+kh(jh+jt-jg+99) &
         +kh(jt+jg-jh+99)+kh(ja+jc-jg+99)+kh(jc+jg-ja+99)+kh(jg+ja-jc+99) &
         +kh(jb+jd-jh+99)+kh(jd+jh-jb+99)+kh(jh+jb-jd+99)+kh(je+jf-jt+99) &
         +kh(jf+jt-je+99)+kh(jt+je-jf+99)
    IF(i <= 0) THEN
       la=(ie+IF+it)/2+2
       ld=(ig+ih+it)/2+2
       ma=(ia+ic+ig)/2+2
       mc=(IF+ic+id)/2+2
       na=(ib+id+ih)/2+2
       nb=(ib+ie+ia)/2+2
       le=(ie+IF-it)/2+1
       lf=(IF+it-ie)/2+1
       lg=(it+ie-IF)/2+1
       me=(ia+ic-ig)/2+1
       mf=(ic+ig-ia)/2+1
       mg=(ig+ia-ic)/2+1
       ne=(ib+id-ih)/2+1
       nf=(id+ih-ib)/2+1
       ng=(ih+ib-id)/2+1
       lx=(it+ig-ih)/2+1
       mx=(ic+id-IF)/2+1
       nx=(ib+ie-ia)/2+1
       fn=q(la,jt+1)*q(jt,lg)*q(ma,jc+1)*q(jc,mf)*q(na,jb+1)*q(jb,ne)
       fd=q(ld,jt+1)*q(jt,lx)*q(mc,jc+1)*q(jc,mx)*q(nb,jb+1)*q(jb,nx)
       jsi=MAX(ABS(je-jh),ABS(jg-jf),ABS(ja-jd))+1
       jsf=MIN(je+jh,jg+jf,ja+jd)-1
       ps=-1-2*(jsi/2*2-jsi)
       fs=ps*SQRT(fn/fd)/FLOAT((ig+1)*(ie+1))
       u=0.
       DO js=jsi,jsf,2
          is=js-1
          lb=(ie+ih+is)/2+2
          lc=(ig+IF+is)/2+2
          mb=(ia+id+is)/2+2
          ly=(ie+ih-is)/2+1
          my=(ig+IF-is)/2+1
          ny=(ia-id+is)/2+1
          ud=q(lb,je+1)*q(je,ly)*q(lc,jg+1)*q(jg,my)*q(mb,js+1)*q(js,ny)
          l0=MAX(la,lb,lc,ld)+1
          m0=MAX(ma,mb,mc,lc)+1
          n0=MAX(na,nb,mb,lb)+1
          l1=MIN(le+ld,lf+lb,lg+lc)
          m1=MIN(me+lc,mf+mb,mg+mc)
          n1=MIN(ne+lb,nf+nb,ng+mb)
          x=0.
          DO l=l0,l1
             x=-x-q(l-1,la)*q(le,l-ld)*q(lf,l-lb)*q(lg,l-lc)
          ENDDO
          y=0.
          DO m=m0,m1
             y=-y-q(m-1,ma)*q(me,m-lc)*q(mf,m-mb)*q(mg,m-mc)
          ENDDO
          z=0.
          DO n=n0,n1
             z=-z-q(n-1,na)*q(ne,n-lb)*q(nf,n-nb)*q(ng,n-mb)
          ENDDO
          ihx=l1+m1+n1
          p=1+2*(ihx/2*2-ihx)
          u=u+p*x*y*z/ud
       ENDDO
       snj=u*fs
    ENDIF

  END FUNCTION snj

  !
  !     This routine calculates the moshinsky vector bracket       
  !     Note that D=mass1/mass2                                    
  !     Ref  m.sotona and m.gmitro  comp.phys.comm 3(1972)53       
  !
  REAL*8 FUNCTION gmosh &
       (n,l,nc,lc,n1,l1,n2,l2,lr,d)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n,l,nc,lc,n1,l1,n2,l2,lr
    REAL*8, INTENT(IN) :: d 
    INTEGER :: ip,ixf,ix, iyi, iyf, j1f,j2,k1i,k1f,m1f,iy,m2f,k2, &
         m2,m2i,m1,j1,k2f,k2i,k1
    REAL*8 :: dl, d1l, bb, ba, anorm, y, p, bc, cfac, bm , &
         sm, s, sxy, bxy

    gmosh=0.
    IF(n+n+nc+nc+l+lc-n1-n1-n2-n2-l1-l2 /= 0 ) RETURN
    IF(l+lc-lr < 0 ) RETURN
    IF(l1+l2-lr < 0 ) RETURN
    IF(ABS(l-lc)-lr > 0 ) RETURN
    IF(ABS(l1-l2)-lr > 0 ) RETURN
    DL=LOG(D)
    D1L=LOG(D+1.)
    bb=f_mb(n1+1)+f_mb(n2+1)+f_mb(n+1)-f_mb(nc+1)+ &
         g_mb(n1+l1+1)+g_mb(n2+l2+1) &
         -g_mb(n+l+1)-g_mb(nc+lc+1) 
    ba=w_mb(l1+1)+w_mb(l2+1)+w_mb(lc+1)+w_mb(l+1)+ &
         f_mb(l1+l2-lr+1)+f_mb(l+lc+lr+2) &
         +f_mb(l+lc-lr+1)+f_mb(lc+lr-l+1)+ &
         f_mb(lr+l-lc+1)-f_mb(l1+l2+lr+2) &
         -f_mb(l1+lr-l2+1)-f_mb(l2+lr-l1+1)-DBLE(l)*d1l
    ip=lr+n+n1+n2
    p=1+2*(ip/2*2-ip)
    anorm=p*EXP(0.5D0*(bb+ba))
    y=0.
    j1f=l+1
    DO j1=1,j1f
       j2=l+2-j1
       k1i=ABS(l1-j1+1)+1
       k1f=l1+j1
       DO k1=k1i,k1f,2
          m1f=n1-(j1+k1-l1)/2+2
          IF(m1f-1 < 0 )  CYCLE
          k2i=MAX(ABS(l2-j2+1),ABS(lc-k1+1))+1
          k2f=MIN(l2+j2,lc+k1)
          IF(k2i-k2f > 0 ) CYCLE
          DO k2=k2i,k2f,2
             m2f=n2-(j2+k2-l2)/2+2
             IF(m2f-1 < 0 )  CYCLE
             ip=j2-1+(l1+k1+j1+l2+k2+j2)/2
             p=1+2*(ip/2*2-ip)
             bc=0.5D0*(DBLE(k1+j2-2)*dl-DBLE(k1+k2-2)*d1l) &
                  +f_mb(k1+l1-j1+1)+f_mb(k1+k2-lc-1)+ &
                  f_mb(k2+l2-j2+1)-f_mb(k1+l1+j1)-f_mb(k1+k2+lc)- &
                  f_mb(k2+l2+j2)+w_mb(k1)+w_mb(k2)+f_mb((k1+l1+j1)/2)+ &
                  f_mb((k1+k2+lc)/2)+f_mb((k2+l2+j2)/2)- &
                  f_mb((k1+l1-j1)/2+1)-f_mb((l1+j1-k1)/2+1)- &
                  f_mb((j1+k1-l1)/2)-f_mb((k1+k2-lc)/2)- &
                  f_mb((k2+lc-k1)/2+1)-f_mb((lc+k1-k2)/2+1) &
                  -f_mb((k2+l2-j2)/2+1)-f_mb((l2+j2-k2)/2+1)- &
                  f_mb((j2+k2-l2)/2)
             cfac=p*EXP(bc)
             sxy=0.
             ixf=MIN(k1+k1,k1+k2-lc)-1
             DO ix=1,ixf
                iyi=MAX(1,ix+j1+l2-k1-lr)
                iyf=MIN(l2+l2+1,l1+l2-lr+1,l2+lc+ix-k1-j2+2)
                IF(iyi-iyf > 0 ) CYCLE
                DO iy=iyi,iyf
                   ip=ix+iy
                   p=1+2*(ip/2*2-ip)
                   bxy=f_mb(k1+k1-ix)+f_mb(l2+l2-iy+2)+ &
                        f_mb(k2+lc-k1+ix)+f_mb(l1+lr-l2+iy) &
                        -f_mb(ix)-f_mb(iy)-f_mb(k1+k2-lc-ix)- &
                        f_mb(l1+l2-lr-iy+2)-f_mb(k1-l2+lr-j1+iy-ix+1)- &
                        f_mb(l2-k1+lc-j2+ix-iy+3)
                   sxy=sxy+p*EXP(bxy)
                ENDDO
             ENDDO
             s=cfac*sxy
             sm=0.
             DO m1=1,m1f
                m2i=MAX(1,nc-m1-(k1+k2-lc)/2+3)
                IF(m2i-m2f > 0 ) CYCLE
                DO m2=m2i,m2f
                   ip=m1+m2
                   p=1+2*(ip/2*2-ip)
                   bm=DBLE(m1-1)*DL-DBLE(m1+m2-2)*d1l+g_mb(1) &
                        +g_mb(m1+m2+(k1+k2+lc)/2-2)-g_mb(k1+m1-1)- &
                        g_mb(k2+m2-1)+f_mb(m1+m2+(k1+k2-lc)/2-2)- &
                        f_mb(m1)-f_mb(m2)-f_mb(n1-m1-(j1+k1-l1)/2+3)- &
                        f_mb(n2-m2-(j2+k2-l2)/2+3) &
                        -f_mb(m1+m2-nc+(k1+k2-lc)/2-2)
                   sm=sm+p*EXP(bm)
                ENDDO
             ENDDO
             y=y+s*sm
          ENDDO
       ENDDO
    ENDDO
    gmosh=anorm*y

  END FUNCTION  gmosh


  ! parity function

  REAL*8 FUNCTION parity(m)
     IMPLICIT NONE
     INTEGER, INTENT(IN) :: m
     parity= -1.0
     if(mod(m,2) == 0)parity=1.0
   END FUNCTION parity

  !  Spherical harmonics from Num. Recipes  

  REAL*8 FUNCTION spherical_harmonics(m1,l,x)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: m1, l
    REAL*8, INTENT(IN) ::  x
    REAL*8, DIMENSION(0:51) :: y
    INTEGER :: iphase, m, j
    REAL*8 :: fj, z, fac, div, sum, a, b, c
    spherical_harmonics=0.
    m=IABS(m1)
    IF(m.LT.0) m=-m1
    y(0)=1.
    IF(l.EQ.0) THEN
       sum=y(0)
    ELSE
       a=m-l
       b=l+m+1
       c=m+1
       z=0.5-x*0.5
       DO j=1,l-m+1
          fj=j-1
          y(j)=y(j-1)*(a+fj)*(b+fj)*z
          div=(c+fj)*(fj+1.)
          y(j)=y(j)/div
       ENDDO
       IF(m > 0) then
          fac=(1.-x*x)**m
          fac=SQRT(fac)
       ELSE
          fac=1.
       ENDIF
       sum=0.
       DO j=0,l-m
          sum=sum+y(j)
       ENDDO
       iphase=m
       IF(m1.LT.0) then
          iphase=0
       ENDIF
       sum=sum*fac*((-1)**iphase)
    ENDIF
    spherical_harmonics=cn(m,l)*sum

  END FUNCTION spherical_harmonics

END MODULE ang_mom_functions





MODULE one_body_operators
  real*8, public, allocatable, dimension(:) :: x_mesh, wx_mesh
  real*8, public, allocatable, dimension(:,:) :: tkin, hcom_1p, vxy
  complex*16, public, allocatable, dimension(:,:) :: sp_basis, fock_mtx
  complex*16, public, allocatable, dimension(:,:,:,:) :: v2body
  complex*16, public, allocatable, dimension(:,:) :: t1_ccm, t1_ccm_eqn, t1_denom, u_hf
END MODULE one_body_operators

!
!
!
module diis_mod
  
  complex*16, public, allocatable, dimension(:,:) :: t2_diis
  complex*16, public, allocatable, dimension(:,:) :: t2_diisd
  integer, public  :: n_diis,nn_diis,nstep, diis_step, subspace 
  complex*16, public, allocatable, dimension(:,:) :: diis_mat, fvec_sub, xvec_sub
  complex*16, public, allocatable, dimension(:)   :: sol_vect
end module diis_mod


MODULE single_particle_orbits
  TYPE, PUBLIC :: single_particle_descript
     INTEGER :: total_orbits
     INTEGER, DIMENSION(:), POINTER :: nx, ny, nz, itzp, szp
     CHARACTER (LEN=10), DIMENSION(:), POINTER :: orbit_status, orb_type, model_space
     REAL*8, DIMENSION(:), POINTER :: e, kx, ky, kz
     
  END TYPE single_particle_descript
  TYPE (single_particle_descript), PUBLIC :: all_orbit, neutron_data, &
       proton_data
  integer, ALLOCATABLE, PUBLIC :: orbitof(:,:,:,:,:)
  real*8, allocatable :: kx_mesh(:), ky_mesh(:), kz_mesh(:), erg(:), xx(:), wxx(:)
  integer, allocatable :: indx(:), indx_inv(:)
  integer :: szmin, szmax, tzmin, tzmax, nmax, my_twists
  REAL*8, ALLOCATABLE, PUBLIC :: twist_angle(:,:)
  
CONTAINS
  SUBROUTINE allocate_sp_array(this_array,n)
    TYPE (single_particle_descript), INTENT(INOUT) :: this_array
    INTEGER , INTENT(IN) :: n
    INTEGER :: I
    IF (ASSOCIATED (this_array%kx) ) DEALLOCATE(this_array%kx)
    ALLOCATE(this_array%kx(n))
    IF (ASSOCIATED (this_array%ky) ) DEALLOCATE(this_array%ky)
    ALLOCATE(this_array%ky(n))
    IF (ASSOCIATED (this_array%kz) ) DEALLOCATE(this_array%kz)
    ALLOCATE(this_array%kz(n))
    IF (ASSOCIATED (this_array%nx) ) DEALLOCATE(this_array%nx)
    ALLOCATE(this_array%nx(n))
    IF (ASSOCIATED (this_array%ny) ) DEALLOCATE(this_array%ny)
    ALLOCATE(this_array%ny(n))
    IF (ASSOCIATED (this_array%nz) ) DEALLOCATE(this_array%nz)
    ALLOCATE(this_array%nz(n))
    IF (ASSOCIATED (this_array%itzp) ) DEALLOCATE(this_array%itzp)
    ALLOCATE(this_array%itzp(n))
    IF (ASSOCIATED (this_array%szp) ) DEALLOCATE(this_array%szp)
    ALLOCATE(this_array%szp(n))
    IF (ASSOCIATED (this_array%e) ) DEALLOCATE(this_array%e)
    ALLOCATE(this_array%e(n))
    IF (ASSOCIATED (this_array%orbit_status) ) DEALLOCATE(this_array%orbit_status)
    ALLOCATE(this_array%orbit_status(n))
    IF (ASSOCIATED (this_array%model_space) ) DEALLOCATE(this_array%model_space)
    ALLOCATE(this_array%model_space(n))
    !           blank all characters and zero all other values
    DO i= 1, n
       this_array%model_space(i)= ' '
       this_array%orbit_status(i)= ' '
       this_array%e(i)=0.
       this_array%kx(i)=0
       this_array%ky(i)=0
       this_array%kz(i)=0
       this_array%nx(i)=0
       this_array%ny(i)=0
       this_array%nz(i)=0
       this_array%itzp(i)=0
       this_array%szp(i)=0
    ENDDO

  END SUBROUTINE allocate_sp_array

  SUBROUTINE deallocate_sp_array(this_array)
    TYPE (single_particle_descript), INTENT(INOUT) :: this_array
    DEALLOCATE(this_array%kx) ; DEALLOCATE(this_array%ky)
    DEALLOCATE(this_array%kz); DEALLOCATE(this_array%nx) ; DEALLOCATE(this_array%ny)
    DEALLOCATE(this_array%nz); DEALLOCATE(this_array%szp) ;DEALLOCATE(this_array%itzp)
    DEALLOCATE(this_array%e);
    DEALLOCATE(this_array%orbit_status); DEALLOCATE(this_array%model_space)
  END SUBROUTINE deallocate_sp_array
END MODULE single_particle_orbits

MODULE energy_variables
  CHARACTER(LEN=100), PUBLIC :: basis_type
  REAL*8, PUBLIC :: hbar_omega, oscl
  INTEGER, PUBLIC :: mass_closed_shell_core
  REAL*8 , PUBLIC :: p_mass, hbarc, hb2ip
  PARAMETER (p_mass =938.926D0, hbarc = 197.327D0)
END MODULE energy_variables




!     Modules specific to the g-matrix calculation and effective operators
!     only
!     arrays containing mesh points and harmonic oscillator wave functions
!     In addition, routines for various wave functions are also included

!     Modules specific to the g-matrix calculation and effective operators
!     only
!     arrays containing mesh points and harmonic oscillator wave functions
!     In addition, routines for various wave functions are also included

MODULE wave_functions

  use constants
      
  INTEGER , PUBLIC :: n_k1, n_k, n_k2
  DOUBLE PRECISION, ALLOCATABLE, PUBLIC :: ra(:),wra(:), krel(:), wkrel(:)
  DOUBLE PRECISION, PUBLIC :: k_cutoff, k_max
  DOUBLE PRECISION, ALLOCATABLE, PUBLIC :: rnlr(:,:,:)
CONTAINS



  
  SUBROUTINE gauss_laguerre(n,x,w,alf)

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n
    INTEGER :: i, its, j
    REAL(dp), INTENT(IN) :: alf 
    REAL(dp), DIMENSION(n), INTENT(INOUT) :: x,w
    REAL(dp), PARAMETER :: eps=3.E-13
    INTEGER, PARAMETER :: maxit=50
    REAL(KIND = 8) :: p1,p2,p3,pp,z,z1, ai

    DO i=1, n
       IF  (i == 1) THEN                         
          z=(1.0D0+alf)*(3.0D0+0.92D0*alf)/(1.0D0+2.4D0*n+1.8D0*alf)
       ELSEIF (i == 2) THEN
          z = z+(15.0D0+6.25D0*alf)/(1.0D0+0.9D0*alf+2.5D0*n)
       ELSE                               
          ai=i-2;
          z = z+((1.0D0+2.55D0*ai)/(1.9D0*ai)+1.26D0*ai*  &
               alf/(1.0D0+3.5D0*ai))*(z-x(i-2))/(1.0D0+0.3D0*alf)
       ENDIF
       DO its=1, MAXIT
          p1=1.0D0;
          p2=0.0D0;
          DO j=1, n
             p3=p2                               
             p2=p1
             p1=((2*j-1+alf-z)*p2-(j-1+alf)*p3)/j
          ENDDO
          pp=(n*p1-(n+alf)*p2)/z
          z1=z
          z=z1-p1/pp                          
          IF (ABS(z-z1) <= EPS) GOTO 10
       ENDDO
10     CONTINUE
       IF (its > MAXIT) THEN
          WRITE(6,*) 'Too many iterations in gauss-Laguerre'
          STOP 
       ENDIF
       x(i)=z                             
       w(i) = -EXP(gammln(alf+n)-gammln(n*1.D0))/(pp*n*p2)
    ENDDO
  END SUBROUTINE gauss_laguerre

  DOUBLE PRECISION FUNCTION gammln(xx)
    IMPLICIT NONE
    INTEGER :: j
    REAL (dp), INTENT(IN) :: xx
    REAL (dp)  cof(6),stp,half,one,fpf,x,tmp,ser
    DATA cof,stp/76.18009173d0,-86.50532033d0,24.01409822d0, &
         -1.231739516d0,.120858003d-2,-.536382d-5,2.50662827465d0/
    DATA half,one,fpf/0.5d0,1.0d0,5.5d0/
    x=xx-one
    tmp=x+fpf
    tmp=(x+half)*LOG(tmp)-tmp
    ser=one
    DO j=1,6
       x=x+one
       ser=ser+cof(j)/x
    ENDDO
    gammln=tmp+LOG(stp*ser)

  END FUNCTION gammln

  DOUBLE PRECISION FUNCTION factrl(n)
    IMPLICIT NONE
    REAL(dp), DIMENSION(33) :: a
    INTEGER, INTENT(IN)  :: n
    INTEGER :: j, ntop
    DATA ntop,a(1)/0,1./

    IF (n <= ntop) THEN
       factrl=a(n+1)
    ELSEIF (n <= 32) THEN
       DO j=ntop+1,n
          a(j+1)=j*a(j)
       ENDDO
       ntop=n
       factrl=a(n+1)
    ELSE
       factrl=exp(gammln(n+1.D0))
    ENDIF

  END  FUNCTION factrl

  SUBROUTINE laguerre_general( n, alpha, x, cx )
    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: n
    REAL (dp ) ::  alpha
    REAL ( dp ) :: cx(0:n)
    !complex ( dpc ) :: cx(0:n)
    INTEGER :: i
    REAL*8, INTENT(IN) ::  x
    !complex ( dpc ), INTENT(IN) ::  x

    IF ( alpha <= -1.0D+00 ) THEN
       WRITE ( *, '(a)' ) ' '
       WRITE ( *, '(a)' ) 'LAGUERRE_GENERAL - Fatal error!'
       WRITE ( *, '(a,g14.6)' ) '  The input value of ALPHA is ', alpha
       WRITE ( *, '(a)' ) '  but ALPHA must be greater than -1.'
       STOP
    END IF
    IF ( n < 0 ) THEN
       RETURN
    END IF
    cx(0) = 1.0D+00
    IF ( n == 0 ) THEN
       RETURN
    END IF
    cx(1) = 1.0D+00 + alpha - x
    DO i = 2, n
       cx(i) = ( ( REAL ( 2 * i - 1, kind = 8 ) + alpha - x ) * cx(i-1)   &
            + ( REAL (   - i + 1, kind = 8 ) - alpha     ) * cx(i-2) ) &
            / REAL (     i,     kind = 8 )
    END DO

  END SUBROUTINE laguerre_general

  DOUBLE PRECISION FUNCTION  fac(m)
    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: m
    INTEGER :: i

    fac = 0.0D0
    IF(m == 0) return
    DO i=1,m
       fac=fac+LOG(FLOAT(i))
    ENDDO

  END FUNCTION  fac

  DOUBLE PRECISION FUNCTION  dfac(m)
    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: m
    INTEGER :: i

    IF (MOD(m,2).ne.1) stop 'wrong argument to dfac'
    dfac = 0.0D0
    IF (m == 1)return
    DO i=3,m,2
       dfac=dfac+LOG(FLOAT(i))
    ENDDO
  END FUNCTION  dfac



  !
  !     H.O. functions using Kummers function   
  !
  REAL*8 FUNCTION rnl(n,l,z)
    IMPLICIT NONE
    INTEGER :: lll, nn
    INTEGER, INTENT(IN) :: l, n
    REAL*8 :: y, dl, gamfaa, pi, dfll, gamfab, dfnn
    REAL*8, INTENT(IN) :: z

    rnl=0. ; y=0.5*z*z
    IF(y > 60.0) RETURN
    dl = l
    IF((ABS(z) < 1.0d-6) .AND. (l == 0)) rnl = 1.0d0
    IF( ABS(z) > 1.0d-6) rnl = (z**l) * EXP(-y) * hypkum(n,dl+1.5,z*z)
    pi = 3.1415926535897932
    gamfaa = 0.5 * SQRT(pi)
    IF(l /= 0) THEN
       DO lll = 1, l
          dfll = lll - 1
          gamfaa = gamfaa * (dfll + 1.5)
       ENDDO
    ENDIF
    gamfab = gamfaa
    IF(n /= 0) THEN
       dfll = dl + 0.5
       DO nn = 1, n
          dfnn = nn
          gamfab = gamfab * ((dfnn + dfll) / dfnn)
       ENDDO
    ENDIF
    rnl = rnl * (SQRT(2.0 * gamfab) / gamfaa)

  END FUNCTION rnl
  !
  !     Kummers function, Abramowitz & Stegun   
  !     exp. 13.1.2. a(there) equals (-n)       
  !  
  REAL*8 FUNCTION hypkum(n,b,z)
    IMPLICIT NONE
    INTEGER :: nmax, nf
    INTEGER, INTENT(IN)  :: n
    REAL*8 :: af, bf, zf, term, dfnf, xadd, sum
    REAL*8, INTENT(IN) :: b, z

    IF(n < 0) WRITE (6,*)' error exit in hypkum ',  n,b,z
    hypkum = 1.0
    IF(n == 0) RETURN
    nmax = n ; af = - n ; bf = b ; zf = z ; sum = 1.0 ; term = 1.0
    DO nf = 1, nmax
       dfnf = nf
       xadd = dfnf - 1.0
       term = term * ((af + xadd) / (bf + xadd)) * (zf / dfnf)
       IF(ABS(term) <  1.0d-12) EXIT
       sum = sum + term
    ENDDO
    hypkum = sum

  END FUNCTION hypkum


  !  This function sets up the recursive relation
  !  for the associated Legendre polynomials

  DOUBLE PRECISION FUNCTION legendre_polynomials(l, m, x)
    IMPLICIT NONE
    DOUBLE PRECISION ::  fact,pll,pmm,pmmp1,somx2
    DOUBLE PRECISION, INTENT(IN)  :: x
    INTEGER ::  i,ll
    INTEGER, INTENT(IN) :: l, m

    !  check whether m, l and x are ok

    IF((M < 0).OR.(M > L).OR.(ABS(X) > 1.)) THEN
       WRITE(6,*) 'bad arguments', m, l, x; RETURN
    ENDIF

    !  calculate now pmm as starting point for iterations

    pmm=1.0
    IF (m > 0) THEN
       somx2=SQRT((1.0-x)*(1.0+x))
       fact=1.0;
       DO i=1, m
          pmm = -fact*somx2*pmm
          fact = fact+2.0
       ENDDO
    ENDIF

    !  if l == m we do not need to use recursion relation

    IF (l == m) THEN
       legendre_polynomials=pmm

       !  recursive relation for associated Legendre polynomials

    ELSE
       pmmp1=x*(2*m+1)*pmm

       !  analytical formula for the case l == m+1

       IF (l == (m+1)) THEN
          legendre_polynomials=pmmp1
       ELSE
          DO ll=m+2, l
             pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m)
             pmm=pmmp1
             pmmp1=pll
          ENDDO
          legendre_polynomials= pll
       ENDIF
    ENDIF

  END FUNCTION legendre_polynomials

  !
  !
  !      This routine calculates gauss-legendre mesh points and weights      
  !      input:                                                              
  !      x1   : lower limit of the integration interval                      
  !      x2   : upper limit ---------- "" -------------                      
  !      n    : the desired number of mesh points                            
  !      output :                                                            
  !      x     : gauss-legendre mesh points on the interval (x1,x2)          
  !      w     : the corresponding weights                                   
  !      From  : Numerical recipes
  !      F90 version : M. Hjorth-Jensen
  !
  SUBROUTINE gauss_legendre(x1,x2,x,w,n)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n
    INTEGER :: i, j, m
    DOUBLE PRECISION, INTENT(IN) :: x1, x2
    DOUBLE PRECISION, INTENT(INOUT) :: x, w
    DOUBLE PRECISION :: eps
    DIMENSION :: x(n), w(n)
    PARAMETER (eps=3.D-14)
    DOUBLE PRECISION :: p1,p2,p3,pp,xl,xm,z,z1

    m=(n+1)/2
    xm=0.5d0*(x2+x1)
    xl=0.5d0*(x2-x1)
    DO i=1,m
       z1=0.
       z=COS(3.141592654d0*(i-.25d0)/(n+.5d0))
       DO WHILE ( ABS(z-z1) > EPS)
          p1=1.
          p2=0.
          DO j=1,n
             p3=p2
             p2=p1
             p1=((2.*j-1.)*z*p2-(j-1.)*p3)/j
          ENDDO
          pp=n*(z*p1-p2)/(z*z-1.)
          z1=z
          z=z-p1/pp
       ENDDO
       x(i)=xm-xl*z
       x(n+1-i)=xm+xl*z
       w(i)=2.*xl/((1.-z*z)*pp*pp)
       w(n+1-i)=w(i)
    ENDDO

  END SUBROUTINE gauss_legendre

END MODULE wave_functions


!     module which defines configuration type, general structure
!     either for lab frame case or relcm system
MODULE configurations
  use parallel
  USE single_particle_orbits
  
  TYPE configuration_descriptor        
     INTEGER :: number_hpphpp_confs, number_hhphhp_confs, number_hhpppp_confs, number_hhhhpp_confs
     INTEGER :: number_hhhh_confs, number_hhhp_confs, number_hhpp_confs, & 
          number_hphp_confs, number_hppp_confs, number_pppp_confs
     INTEGER :: number_ph_phhp_confs, number_ph_hphp_confs, number_ph_hpph_confs
     INTEGER :: number_confs, model_space_confs, number_t3_confs
     INTEGER, DIMENSION(:), POINTER :: config_ab, config_abc
     INTEGER, DIMENSION(:), POINTER :: t2_quantum_numbers
     INTEGER, DIMENSION(:), POINTER :: t3_quantum_numbers
     
     INTEGER, DIMENSION(:), POINTER :: hhhh_quantum_numbers
     INTEGER, DIMENSION(:), POINTER :: hhhp_quantum_numbers
     INTEGER, DIMENSION(:), POINTER :: hhpp_quantum_numbers
     INTEGER, DIMENSION(:), POINTER :: hphp_quantum_numbers
     INTEGER, DIMENSION(:), POINTER :: hppp_quantum_numbers
     INTEGER, DIMENSION(:), POINTER :: pppp_quantum_numbers
     
     INTEGER, DIMENSION(:), POINTER :: ph_phhp_quantum_numbers
     INTEGER, DIMENSION(:), POINTER :: ph_hphp_quantum_numbers
     INTEGER, DIMENSION(:), POINTER :: ph_hpph_quantum_numbers

     INTEGER, DIMENSION(:), POINTER :: hpphpp_quantum_numbers
     INTEGER, DIMENSION(:), POINTER :: hhphhp_quantum_numbers
     INTEGER, DIMENSION(:), POINTER :: hhpppp_quantum_numbers
     INTEGER, DIMENSION(:), POINTER :: hhhhpp_quantum_numbers
     
  END TYPE configuration_descriptor
  
  TYPE (configuration_descriptor), PUBLIC :: channels
    
CONTAINS

  
  
  !                       
  !     setting up all configurations for given J, Tz and parity for
  !     the g-matrix model space

  SUBROUTINE number_confs(nx,ny,nz,sz,tz,struct,this)
    USE constants
    USE single_particle_orbits
    IMPLICIT NONE
    TYPE (configuration_descriptor), INTENT(INOUT) :: this
    INTEGER :: a,b,tza,tzb,sza,szb,sz, tz, struct, a_max, b_max, a_min, b_min, nconfs
    INTEGER :: nx, ny, nz, nxa, nya, nza, nxb, nyb, nzb  
    
    a_max = tot_orbs
    b_max = tot_orbs
    
    a_min = 1
    b_min = 1
    
    if ( struct == 1 ) THEN
       a_max = below_ef
       b_max = below_ef
    ELSEIF ( struct == 2 ) THEN
       a_min = 1
       b_min = below_ef + 1
       
       a_max = below_ef
       b_max = tot_orbs
       
    ELSEIF ( struct == 3 ) THEN
       a_min = below_ef + 1
       b_min = below_ef + 1
       
       a_max = tot_orbs
       b_max = tot_orbs
    end if
    
    nconfs=0
    DO a=a_min, a_max
       nxa = all_orbit%nx(a)
       nya = all_orbit%ny(a)
       nza = all_orbit%nz(a)
       tza=all_orbit%itzp(a)
       sza=all_orbit%szp(a)
       
       DO b=b_min, b_max
          
          nxb = all_orbit%nx(b)
          nyb = all_orbit%ny(b)
          nzb = all_orbit%nz(b)
          tzb=all_orbit%itzp(b)
          szb=all_orbit%szp(b)
          
          IF (tzb+tza /= tz*2 ) CYCLE
          IF (szb+sza /= sz*2 ) CYCLE
          if ( nxa + nxb /= Nx ) CYCLE
          if ( nya + nyb /= Ny ) CYCLE
          if ( nza + nzb /= Nz ) CYCLE
          

          nconfs=nconfs+1
       ENDDO
    ENDDO
    this%number_confs=nconfs
    
  END SUBROUTINE number_confs

  !
  !
  !
  SUBROUTINE setup_confs(nx,ny,nz,sz,tz,struct,this)
    USE constants
    USE single_particle_orbits
    IMPLICIT NONE
    TYPE (configuration_descriptor), INTENT(INOUT) :: this
    INTEGER :: a,b,tza,tzb,sza,szb,sz, tz, struct, a_max, b_max, a_min, b_min, nconfs
    INTEGER :: nx, ny, nz, nxa, nya, nza, nxb, nyb, nzb , k1, k2
    
    a_max = tot_orbs
    b_max = tot_orbs
    
    a_min = 1
    b_min = 1
    
    if ( struct == 1 ) THEN
       a_max = below_ef
       b_max = below_ef
    ELSEIF ( struct == 2 ) THEN
       a_min = 1
       b_min = below_ef + 1
       
       a_max = below_ef
       b_max = tot_orbs
       
    ELSEIF ( struct == 3 ) THEN
       a_min = below_ef + 1
       b_min = below_ef + 1
       
       a_max = tot_orbs
       b_max = tot_orbs
    end if
    
    nconfs=0
    DO a=a_min, a_max
       nxa = all_orbit%nx(a)
       nya = all_orbit%ny(a)
       nza = all_orbit%nz(a)
       tza=all_orbit%itzp(a)
       sza=all_orbit%szp(a)
       
       DO b=b_min, b_max
          
          nxb = all_orbit%nx(b)
          nyb = all_orbit%ny(b)
          nzb = all_orbit%nz(b)
          tzb=all_orbit%itzp(b)
          szb=all_orbit%szp(b)
          
          IF (tzb+tza /= tz*2 ) CYCLE
          IF (szb+sza /= sz*2 ) CYCLE
          if ( nxa + nxb /= Nx ) CYCLE
          if ( nya + nyb /= Ny ) CYCLE
          if ( nza + nzb /= Nz ) CYCLE
          
          
          nconfs=nconfs+1
          k2=nconfs*2
          k1=k2-1
          this%config_ab(k1)=b
          this%config_ab(k2)=a
 
       ENDDO
    ENDDO
  end SUBROUTINE setup_confs

END MODULE configurations


!
!
!
MODULE t2_storage
  
  INTEGER, ALLOCATABLE, PUBLIC :: my_t3channel_low(:), my_t3channel_high(:), my_ph_hphp_channel_low(:), my_ph_hphp_channel_high(:)
  INTEGER, ALLOCATABLE, PUBLIC :: my_channel_low(:), my_channel_high(:), my_hhpp_channel_low(:), my_hhpp_channel_high(:)
  INTEGER, ALLOCATABLE, PUBLIC :: mapping(:,:,:),mapping_t3(:,:,:), mapping_hhpp(:,:,:),mapping_hppp(:,:,:),mapping_hphp(:,:,:)
  INTEGER, ALLOCATABLE, PUBLIC :: my_hpphpp_channel_low(:), my_hpphpp_channel_high(:)
  INTEGER, ALLOCATABLE, PUBLIC :: my_ph_hpphpp_channel_low(:), my_ph_hpphpp_channel_high(:)
  INTEGER, ALLOCATABLE, PUBLIC :: my_hhphhp_channel_low(:), my_hhphhp_channel_high(:)
  INTEGER, ALLOCATABLE, PUBLIC :: my_hhpppp_channel_low(:), my_hhpppp_channel_high(:)
  INTEGER, ALLOCATABLE, PUBLIC :: my_hhhhpp_channel_low(:), my_hhhhpp_channel_high(:)
  INTEGER, ALLOCATABLE, PUBLIC :: mapping_hpphpp(:,:,:), mapping_ph_hpphpp(:,:,:),  & 
       mapping_hhphhp(:,:,:), mapping_hhpppp(:,:,:), mapping_hhhhpp(:,:,:)
  INTEGER, ALLOCATABLE, PUBLIC :: mapping_hhhh(:,:,:), check_my_channel_hhhh(:)
  INTEGER, ALLOCATABLE, PUBLIC :: mapping_ph_hphp(:,:,:), check_my_channel_ph_hphp(:)
  INTEGER, ALLOCATABLE, PUBLIC :: check_my_t3channel(:), fully_stored_channel(:)
  INTEGER, ALLOCATABLE, PUBLIC :: check_my_channel(:), check_my_channel_hhpp(:), check_my_channel_hppp(:), check_my_channel_hphp(:)
  INTEGER, ALLOCATABLE, PUBLIC :: check_my_channel_hpphpp(:), check_my_channel_hhphhp(:), & 
       check_my_channel_hhpppp(:), check_my_channel_hhhhpp(:), check_my_channel_ph_hpphpp(:)
  INTEGER, ALLOCATABLE, public :: LOCATE_T2CHANNEL(:,:,:,:)
  INTEGER, ALLOCATABLE, public :: LOCATE_CHANNEL(:,:,:,:,:)
  INTEGER, ALLOCATABLE, public :: LOCATE_ph_CHANNEL(:,:,:,:,:)
  INTEGER, ALLOCATABLE, public :: LOCATE_V3CHANNEL(:,:,:,:,:)
  
  TYPE, PUBLIC :: block_storage
     complex*16, DIMENSION(:), ALLOCATABLE :: val1
     complex*16, DIMENSION(:,:), ALLOCATABLE :: val
     complex*16, DIMENSION(:,:,:), ALLOCATABLE :: val3
  END TYPE block_storage
  
  TYPE, PUBLIC :: integer_storage
     INTEGER, DIMENSION(:,:), ALLOCATABLE :: ival
     INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: ival3
  END TYPE integer_storage
  
  
  
  TYPE (integer_storage), PUBLIC :: hhh_hhhppp
  TYPE (integer_storage), PUBLIC :: ppp_hhhppp
  
  TYPE (integer_storage), ALLOCATABLE, PUBLIC :: pq_conf(:)
  TYPE (integer_storage), PUBLIC :: hh_hhhh
  TYPE (integer_storage), PUBLIC :: hh_hhhp
  TYPE (integer_storage), PUBLIC :: hh_hhpp
  TYPE (integer_storage), PUBLIC :: hp_hphp
  TYPE (integer_storage), PUBLIC :: hp_hppp
  TYPE (integer_storage), PUBLIC :: hp_hhhp
  TYPE (integer_storage), PUBLIC :: pp_hhpp
  TYPE (integer_storage), PUBLIC :: pp_hppp
  TYPE (integer_storage), PUBLIC :: pp_pppp
  TYPE (integer_storage), PUBLIC :: pp_pppp_rest

  TYPE (integer_storage), PUBLIC :: ph_ph_phhp
  TYPE (integer_storage), PUBLIC :: hp_ph_phhp
  TYPE (integer_storage), PUBLIC :: hp_ph_hphp
  TYPE (integer_storage), PUBLIC :: hp_ph_hpph
  TYPE (integer_storage), PUBLIC :: ph_ph_hpph
  
  TYPE (integer_storage), ALLOCATABLE, PUBLIC :: lookup_t2_configs(:,:) 
  TYPE (integer_storage), ALLOCATABLE, PUBLIC :: lookup_t3_configs(:,:) 
  TYPE (integer_storage), ALLOCATABLE, PUBLIC :: lookup_hhhh_configs(:,:) 
  TYPE (integer_storage), ALLOCATABLE, PUBLIC :: lookup_hhhp_configs(:,:) 
  TYPE (integer_storage), ALLOCATABLE, PUBLIC :: lookup_hhpp_configs(:,:) 
  TYPE (integer_storage), ALLOCATABLE, PUBLIC :: lookup_hphp_configs(:,:) 
  TYPE (integer_storage), ALLOCATABLE, PUBLIC :: lookup_hppp_configs(:,:) 
  TYPE (integer_storage), ALLOCATABLE, PUBLIC :: lookup_pppp_configs(:,:) 
  TYPE (integer_storage), ALLOCATABLE, PUBLIC :: lookup_pqrs_configs(:,:) 
  
  TYPE (integer_storage), ALLOCATABLE, PUBLIC :: lookup_ph_phhp_configs(:,:) 
  TYPE (integer_storage), ALLOCATABLE, PUBLIC :: lookup_ph_hphp_configs(:,:) 
  TYPE (integer_storage), ALLOCATABLE, PUBLIC :: lookup_ph_hpph_configs(:,:) 
  
  TYPE (integer_storage), ALLOCATABLE, PUBLIC :: lookup_hpphpp_configs(:,:) 
  TYPE (integer_storage), ALLOCATABLE, PUBLIC :: lookup_hhphhp_configs(:,:) 
  TYPE (integer_storage), ALLOCATABLE, PUBLIC :: lookup_hhpppp_configs(:,:) 
  TYPE (integer_storage), ALLOCATABLE, PUBLIC :: lookup_hhhhpp_configs(:,:) 
  
  TYPE (integer_storage), ALLOCATABLE, PUBLIC :: lookup_ph_hpphpp_configs(:,:) 
  
  ! 
  ! Intermediates for t1 and t2 amps
  ! 
  complex*16, ALLOCATABLE, PUBLIC :: I5ab(:,:), I6abc(:,:), I6ac(:,:), I7abcde(:,:) 
  complex*16, ALLOCATABLE, PUBLIC :: I5(:,:), I6(:,:), I7(:,:) 
  complex*16, ALLOCATABLE, PUBLIC :: D1ik(:), D1ac(:) 
  complex*16, ALLOCATABLE, PUBLIC :: d2ij(:),d2ab(:) 
  TYPE (block_storage), ALLOCATABLE, PUBLIC :: interm_hhphhp(:), v3nf_hpphpp(:), v3nf_ph_hpphpp(:)
  TYPE (block_storage), ALLOCATABLE, PUBLIC :: v3nf_hhphhp(:),v3nf_hhpppp(:),v3nf_hhhhpp(:)
  TYPE (block_storage), ALLOCATABLE, PUBLIC :: vnn_hhhh(:), v3nf_hhhh(:), v3nf_ppphhh(:)
  TYPE (block_storage), ALLOCATABLE, PUBLIC :: vnn_hhhp(:)
  TYPE (block_storage), ALLOCATABLE, PUBLIC :: vnn_hhpp(:), v3nf_hhpp(:)
  TYPE (block_storage), ALLOCATABLE, PUBLIC :: vnn_hphp(:), v3nf_hphp(:)
  TYPE (block_storage), ALLOCATABLE, PUBLIC :: vnn_hppp(:)
  TYPE (block_storage), ALLOCATABLE, PUBLIC :: vnn_pppp(:)
  TYPE (block_storage), ALLOCATABLE, PUBLIC :: interm_pppp(:)
  TYPE (block_storage), ALLOCATABLE, PUBLIC :: interm_hhhh(:)
  TYPE (block_storage), ALLOCATABLE, PUBLIC :: interm_hphp(:)
  
  TYPE (block_storage), ALLOCATABLE, PUBLIC :: vnn_ph_hpph(:)
  TYPE (block_storage), ALLOCATABLE, PUBLIC :: I2(:)
  TYPE (block_storage), ALLOCATABLE, PUBLIC :: I3ab(:)
  TYPE (block_storage), ALLOCATABLE, PUBLIC :: I3(:)
  TYPE (block_storage), ALLOCATABLE, PUBLIC :: I4abc(:)
  TYPE (block_storage), ALLOCATABLE, PUBLIC :: I4(:), I4_tmp(:), I4_ph(:)
  TYPE (block_storage), ALLOCATABLE, PUBLIC :: I11b_hp(:)
  TYPE (block_storage), ALLOCATABLE, PUBLIC :: I8(:)
  TYPE (block_storage), ALLOCATABLE, PUBLIC :: I9abc(:)
  TYPE (block_storage), ALLOCATABLE, PUBLIC :: I9(:)
  TYPE (block_storage), ALLOCATABLE, PUBLIC :: I10(:)
  TYPE (block_storage), ALLOCATABLE, PUBLIC :: I10af(:)
  TYPE (block_storage), ALLOCATABLE, PUBLIC :: I11abdef(:), I11(:)
  TYPE (block_storage), ALLOCATABLE, PUBLIC :: temp_mat(:)
  TYPE (block_storage), ALLOCATABLE, PUBLIC :: temp_mat2(:)

  TYPE (block_storage), ALLOCATABLE, PUBLIC :: t2_ccm(:)
  TYPE (block_storage), ALLOCATABLE, PUBLIC :: t2_ccm_eqn(:)
  
  TYPE (block_storage), ALLOCATABLE, PUBLIC :: t2_ph_ccm(:)
  TYPE (block_storage), ALLOCATABLE, PUBLIC :: t2_ph_ccm_eqn(:)
  
  TYPE (block_storage), ALLOCATABLE, PUBLIC :: t3_ccm(:)

END MODULE t2_storage
