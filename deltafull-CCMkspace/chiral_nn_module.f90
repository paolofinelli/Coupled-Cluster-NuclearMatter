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




module chiral_constants 
  use constants, only : pi  
  
  ! chiral order definition
  INTEGER, PARAMETER :: LO   = 0
  ! all contributions vanish at order 1
  ! due to parity and time-reversal  invariance
  INTEGER, PARAMETER :: NLO  = 2
  INTEGER, PARAMETER :: NNLO = 3
  INTEGER, PARAMETER :: N3LO = 4
  
  real*8, parameter :: LEC_c1 =  -0.9186395287347203d0
  real*8, parameter :: LEC_c2 =   0.d0 
  real*8, parameter :: LEC_c3 =  -3.888687492763241d0
  real*8, parameter :: LEC_c4 =   4.310327160829740d0
  real*8, parameter :: LEC_c1_3NF =  -0.9186395287347203d0
  real*8, parameter :: LEC_c2_3NF =   0.d0 
  real*8, parameter :: LEC_c3_3NF =  -3.888687492763241d0
  real*8, parameter :: LEC_c4_3NF =   4.310327160829740d0
  real*8, parameter :: gA = 1.29d0 
  real*8, parameter :: gA2 = gA*gA 
  real*8, parameter :: gA4 = gA2*gA2 
  real*8, parameter :: fpi = 92.4d0 
  real*8, parameter :: fpi2 = fpi*fpi
  real*8, parameter :: fpi4 = fpi2*fpi2
  real*8, parameter :: twopi = pi*2.d0
  real*8, parameter :: pi2 = pi*pi
  real*8 :: const(100) 
  real*8, parameter :: proton_mass = 938.272d0 
  real*8, parameter :: neutron_mass = 939.5653d0
  real*8 :: mnuc(-1:1), mpi(-1:2)
  REAL*8 :: mnuc2(-1:1)   ! nucleon mass squared
  REAL*8 :: mpi2(-1:2)    ! pion mass squared
  REAL*8 :: mpi4(-1:2)    ! mpi2 squared
  REAL*8 :: mpi5(-1:2)    ! mpi^5
  REAL*8 :: twompi(-1:2)  ! two times pion mass
  REAL*8 :: fourmpi2(-1:2)! four times pion mass squared
  COMPLEX*16 :: sigma_x(2,2), sigma_y(2,2), sigma_z(2,2)
  REAL*8 :: c1,c2,c3,c4
  REAL*8 :: C1_3NF, C2_3NF, C3_3NF, C4_3NF
  REAL*8, PARAMETER :: lambda = 500.D0 
  REAL*8, PARAMETER :: sfr = 700.d0 
  REAL*8, PARAMETER :: sfr2 = sfr*sfr
  CHARACTER(LEN=2), PARAMETER :: chp_itope = 'EM'
  REAL*8 :: sfr_heavyside(-1:2) ! THETA(sfr-twompi)
  REAL*8 :: CS(-1:1), CT(-1:1), cnlo(1:7) 
  REAL*8 :: CE, CD, Econst, Dconst 
  
CONTAINS 
  
  subroutine init_chp_constants 
    implicit none   
    double precision :: c1s0(-1:1), c3s1(-1:1), cnlo_pw(1:7)
    real*8 :: lambdachi
    integer :: i 
    
    mnuc(-1) = proton_mass
    mnuc(0)  = 0.5d0*(proton_mass+neutron_mass)
    mnuc(+1) = neutron_mass
    
    mpi(-1)  = 139.5702d0
    mpi(0)   = 134.9766d0
    mpi(+1)  = 139.5702d0
    mpi(+2)  = ( mpi(-1)+mpi(0)+mpi(+1) ) /3.d0  
    
    
    mnuc2(:)   = mnuc*mnuc 
    mpi2(:)    = mpi*mpi   
    mpi4(:)    = mpi2*mpi2
    mpi5(:)    = mpi4*mpi 
    twompi(:)  = 2.0D0*mpi 
    fourmpi2(:)= 4.0D0*mpi2
     
    c1 = LEC_c1*1.0D-3
    c2 = LEC_c2*1.0D-3
    c3 = LEC_c3*1.0D-3
    c4 = LEC_c4*1.0D-3

    c1_3NF = LEC_c1_3NF*1.0D-3
    c2_3NF = LEC_c2_3NF*1.0D-3
    c3_3NF = LEC_c3_3NF*1.0D-3
    c4_3NF = LEC_c4_3NF*1.0D-3
  
    const = 0.0D0
    !1: 1/(2pi)^3
    const(1) = 1.0D0/(twopi**3)
    !2: gA^2/(4*fpi^2)
    const(2) = gA2/(4.0D0*fpi2)
    !3: 384pi^2*fpi^4
    const(3) = 384.0D0*pi2*fpi4
    !4: 5gA^4-4gA^2-1
    const(4) = 5.0D0*gA4-4.0D0*gA2-1.0D0
    !5: 23gA^4-10gA^2-1
    const(5) = 23.0D0*gA4-10.0D0*gA2 - 1.0D0
    !6: 48gA^4
    const(6) = 48.0D0*gA4
    !7: 3gA^4
    const(7) = 3.0D0*gA4
    !8: 64pi^2fpi^4
    const(8) = 64.0D0*pi2*fpi4
    !9: 3gA^2/16pifpi^4
    const(9) = 3.0D0*gA2/(16.0D0*pi*fpi4)
    !10: ga^2/16
    const(10) = ga2/16.0D0
    !11: 2.0D0*(2c1-c3)
    const(11) = 2.0D0*(2.0D0*c1-c3)
    !12 : const(7)/256pifpi^4
    const(12) = const(7)/(256.0D0*pi*fpi4)
    !13: gA^2/(128pifpi^4)
    const(13) = gA2/(128.0D0*pi*fpi4)
    !14: gA^4/(128pifpi^4)
    const(14) = gA4/(128.0D0*pi*fpi4)
    !15: 3gA^4/(512pifpi^4)
    const(15) = 3.0D0*gA4/(512.0D0*pi*fpi4)
    !16: gA2/(32pifpi^4)
    const(16) = gA2/(32.0D0*pi*fpi4)
    !17: gA2/8
    const(17) = gA2/8.0D0
    !18: gA4/(256pifpi^4)
    const(18) = gA4/(256.0D0*pi*fpi4)
    !19: 3gA4/(32pifpi^4)
    const(19) = 3.0D0*gA4/(32.0D0*pi*fpi4)
    !20: const(16)*(1-gA2)
    const(20) = const(16)*(1.0D0-gA2)
    !21: gA/(8*fpi^2)
    const(21) = gA/(8.0D0*fpi2)
    
    sigma_x = 0.d0
    sigma_y = 0.d0
    sigma_z = 0.d0 
    
    sigma_x(1,2) = 1.d0 
    sigma_x(2,1) = 1.d0 
    
    sigma_z(1,1) = 1.d0 
    sigma_z(2,2) = -1.d0 
    
    sigma_y(1,2) = dcmplx(0.d0,-1.d0) 
    sigma_y(2,1) = dcmplx(0.d0, 1.d0) 
    
    DO i=-1,2
       IF ( (sfr-twompi(i)) <  0.0D0 ) sfr_heavyside(i) = 0.0D0
       IF ( (sfr-twompi(i)) >= 0.0D0 ) sfr_heavyside(i) = 1.0D0
    END DO
    

    !
    ! leading order contacts in PW
    !
    c1s0(-1) = -0.1513660372031080D+00
    c1s0(0)  = -0.1521410882366787D+00
    c1s0(1)  = -0.1517647459006913D+00
    
    c3s1(-1) = -0.1584341766228121D+00
    c3s1(0)  = -0.1584341766228121D+00
    c3s1(1)  = -0.1584341766228121D+00
    
    ! LO-CONTACTS IF THIS IS MAXIMUM ORDER
    ! THEN NO CIB
    !c1s0 =  -0.1521410882366787D+00
    !c3s1 =  -0.1584341766228121D+00

    ! the LO CIBcontacts are input in units of 10^4/GeV^2
    ! this = 10^-2/MeV^2 
    c1s0 = c1s0*0.01d0 
    c3s1 = c3s1*0.01d0 
    
    CS = (c1s0 + 3.d0*c3s1) /16.d0/pi
    CT = (c3s1 - c1s0) /16.d0/pi
    
    
    !
    ! next-to-leading order contacts in PW
    !
    cnlo_pw(1) =  0.2404021944134705D+01
    cnlo_pw(2) =  0.1263390763475578D+01
    cnlo_pw(3) =  0.4170455420556486D+00
    cnlo_pw(4) = -0.7826584999752046D+00
    cnlo_pw(5) =  0.9283846626623043D+00
    cnlo_pw(6) =  0.6181414190474576D+00
    cnlo_pw(7) = -0.6778085114063558D+00
    
    ! NLO
    ! the NLO contacts are input in units of 10^4/GeV^4
    ! this = 10^-8/MeV^4
    DO i=1,7
       CNLO_PW(i) = CNLO_PW(i) * 1.D-08
    END DO
    
    cnlo(1) = (-5.d0*cnlo_pw(7)+6.d0*cnlo_pw(5)-3.d0*cnlo_pw(4)-3.d0*cnlo_pw(3)-cnlo_pw(2)+2.d0*cnlo_pw(1))/(64.d0*pi)
    cnlo(2) = ( 5.d0*cnlo_pw(7)+6.d0*cnlo_pw(5)+3.d0*cnlo_pw(4)+3.d0*cnlo_pw(3)+cnlo_pw(2)+2.d0*cnlo_pw(1))/(16.d0*pi)
    cnlo(3) = -( 2.d0*cnlo_pw(7)-2.d0*dsqrt(2.d0)*cnlo_pw(6)-2.d0*cnlo_pw(5)-3.d0*cnlo_pw(3)+cnlo_pw(2)+2.d0*cnlo_pw(1))/(64.d0*pi)
    cnlo(4) = -(-2.d0*cnlo_pw(7)-2.d0*dsqrt(2.d0)*cnlo_pw(6)-2.d0*cnlo_pw(5)+3.d0*cnlo_pw(3)-cnlo_pw(2)+2.d0*cnlo_pw(1))/(16.d0*pi)
    cnlo(5) = -(-5.d0*cnlo_pw(7)+3.d0*cnlo_pw(4)+2.d0*cnlo_pw(2))/(16.d0*pi)
    cnlo(6) = ( cnlo_pw(7)-6.d0*dsqrt(2.d0)*cnlo_pw(6)-3.d0*cnlo_pw(4)+2.d0*cnlo_pw(2))/(64.d0*pi)
    cnlo(7) = -(cnlo_pw(7)+6.d0*dsqrt(2.d0)*cnlo_pw(6)-3.d0*cnlo_pw(4)+2.d0*cnlo_pw(2))/(16.d0*pi) 
    
    
    !WRITE(6,"(A12,F30.16)") 'C1', cnlo(1)* 1.D8
    !WRITE(6,"(A12,F30.16)") 'C2', cnlo(2)* 1.D8
    !WRITE(6,"(A12,F30.16)") 'C3', cnlo(3)* 1.D8
    !WRITE(6,"(A12,F30.16)") 'C4', cnlo(4)* 1.D8
    !WRITE(6,"(A12,F30.16)") 'C5', cnlo(5)* 1.D8
    !WRITE(6,"(A12,F30.16)") 'C6', cnlo(6)* 1.D8
    !WRITE(6,"(A12,F30.16)") 'C7', cnlo(7)* 1.D8
    
    
    ! NNLO 3NF constants
    !
    !
    !read(5,*);read(5,*) cE, cD 
    

    !cE = -0.398
    !cD = -0.39
    lambdaChi = 700 ! MeV
    Econst = cE/(fpi4*lambdaChi)
    Dconst = cD/(fpi2*lambdaChi)
    

  end subroutine init_chp_constants
 
  subroutine init_cEcD
    
    real*8 :: lambdachi
    lambdaChi = 700 ! MeV
    Econst = cE/(fpi4*lambdaChi)
    Dconst = cD/(fpi2*lambdaChi)


  end subroutine init_cEcD

  ! static one pion exchange, [1] eq 4.5
  ! without isospin structure
  ! q2  : momentum transfer squared
  ! impi: determines which mpi2 to use, 
  FUNCTION chp_one_pion_exchange(q2, impi) RESULT(res)
  
    IMPLICIT NONE 
    REAL*8 , INTENT(IN) :: q2
    INTEGER, INTENT(IN) :: impi
    REAL*8 :: res
    
    res = 0.0D0
    res = -1.0D0 * const(2)/ (q2 + mpi2(impi))
    
        
  END FUNCTION chp_one_pion_exchange
  

  ! NLO function Eq 4.9
  ! q2  : momentum transfer squared
  ! impi: determines which mpi2 to use, 
  FUNCTION chp_NLO_two_pion_exchange_Wc(q2, L, w, impi) RESULT(res)
    
    REAL*8 , INTENT(IN) :: q2
    REAL*8 , INTENT(IN) :: L
    REAL*8 , INTENT(IN) :: w
    INTEGER, INTENT(IN) :: impi
    REAL*8 :: res
    
    res = 0.0D0
    res = -L *(fourmpi2(impi)*const(4) + q2*const(5) + const(6)*mpi4(impi)/(w*w))/const(3) 
!    write(6,*) L,fourmpi2(impi),const(4)
    
  END FUNCTION chp_NLO_two_pion_exchange_Wc
  
  ! NLO function Eq 4.10
  ! q2  : momentum transfer squared
  ! impi: determines which mpi2 to use, 
  FUNCTION chp_NLO_two_pion_exchange_Vs(q2,L,impi) RESULT(res)
    
    REAL*8 , INTENT(IN) :: q2
    REAL*8 , INTENT(IN) :: L
    INTEGER, INTENT(IN) :: impi
    REAL*8 :: res
    
    res = 0.0D0
    res = +const(7)*L*q2/const(8)
    
  END FUNCTION chp_NLO_two_pion_exchange_Vs

  ! NLO Vt function Eq 4.10
  ! q2  : momentum transfer squared
  ! impi: determines which mpi2 to use, 
  FUNCTION chp_NLO_two_pion_exchange_VT(L,impi) RESULT(res)
    
    REAL*8 , INTENT(IN) :: L
    INTEGER, INTENT(IN) :: impi
    REAL*8 :: res
    
    res = 0.0D0 
    res = -const(7)*L/const(8)
    
  END FUNCTION chp_NLO_two_pion_exchange_VT

  ! NLO loop function w [1] Eq 4.12 (DR and SFR)
  ! q2  : momentum transfer squared
  ! impi: determines which mpi2 to use, 
  FUNCTION chp_NLO_two_pion_exchange_loop_w(q2, impi) RESULT(res)
    
    REAL*8 , INTENT(IN) :: q2
    INTEGER, INTENT(IN) :: impi
    REAL*8 :: res
    
    res = 0.0D0
    res = DSQRT(fourmpi2(impi) + q2)
    
  END FUNCTION chp_NLO_two_pion_exchange_loop_w

  ! NLO SFR loop function s [2] Eq 2.16
  ! impi: determines which mpi2 to use, 
  FUNCTION chp_NLO_sfr_two_pion_exchange_loop_s(impi) RESULT(res)
    
    INTEGER, INTENT(IN) :: impi
    REAL*8 :: res
    
    res = 0.0D0
    res = DSQRT(sfr2 - fourmpi2(impi))

  END FUNCTION chp_NLO_sfr_two_pion_exchange_loop_s
  
  ! NLO SFR loop function L [2] Eq 2.16
  ! q   : momentum transfer
  ! q2  : momentum transfer squared
  ! w,s : SFR loop functions
  ! impi: determines which mpi2 to use, 
  FUNCTION chp_NLO_sfr_two_pion_exchange_loop_L(q, q2, w, s, impi) RESULT(res)
    
    REAL*8 , INTENT(INOUT) :: q
    REAL*8 , INTENT(INOUT) :: q2
    REAL*8 , INTENT(IN) :: w
    REAL*8 , INTENT(IN) :: s
    INTEGER, INTENT(IN) :: impi
    REAL*8 :: res, eps 
    
    eps = 1.0D-8 
    res = 0.0D0
    IF (sfr_heavyside(impi) == 0.0D0) return
    
    if ( q == 0.d0 ) then  
       q = q + eps 
       q2 = q*q 
    end if

    res = w * dlog( (sfr2*w*w+q2*s*s+2.0D0*sfr*q*w*s)/( fourmpi2(impi)*(sfr2+q2) ) )/(2.0D0*q)
    !write(6,'(3g)') q, q2, res !(sfr-twompi(i)), sfr, sfr2, q,q2 ! (sfr2*w*w+q2*s*s+2.0D0*sfr*q*w*s), ( fourmpi2(impi)*(sfr2+q2) )
    !write(6,*) 0.5d0*( 2.d0*s/sfr/dsqrt(fourmpi2(impi)) - sfr2*fourmpi2(impi)) * w 
    
  END FUNCTION chp_NLO_sfr_two_pion_exchange_loop_L
  
   ! NNLO loop function wtilde SQUARED [1] Eq 4.20 (DR)
  ! q2  : momentum transfer squared
  ! impi: determines which mpi2 to use, 
  FUNCTION chp_NNLO_two_pion_exchange_loop_wtilde2(q2, impi) RESULT(res)
    
    REAL*8 , INTENT(IN) :: q2
    INTEGER, INTENT(IN) :: impi
    REAL*8 :: res
    
    res = 0.0D0
    res = 2.0D0*mpi2(impi) + q2
    
  END FUNCTION chp_NNLO_two_pion_exchange_loop_wtilde2

  ! NNLO loop function wtilde [2] Eq 2.17 (SFR)
  ! q   : momentum transfer
  ! q2  : momentum transfer squared
  ! impi: determines which mpi to use, 
  FUNCTION chp_NNLO_sfr_two_pion_exchange_loop_A(q, q2, impi) RESULT(res)
    
    REAL*8 , INTENT(INOUT) :: q
    REAL*8 , INTENT(INOUT) :: q2
    INTEGER, INTENT(IN) :: impi
    REAL*8 :: res
    REAL*8 :: eps
    
    res = 0.0D0
    eps = 1.0D-8
    IF (sfr_heavyside(impi) == 0.0D0) return
    
    if ( q == 0.d0 ) then  
       q = q + eps 
       q2 = q*q 
    end if
    
    res = datan( q*(sfr-twompi(impi) )/(q2 + sfr*twompi(impi) ) )/(2.0D0*q)
    
  END FUNCTION chp_NNLO_sfr_two_pion_exchange_loop_A

  ! NNLO function Eq [1] 4.13 and 4.21
  ! q2    : momentum transfer squared
  ! w     : Eq 4.12
  ! A     : NNLO loop function A Eq. 4.19
  ! wt2   : NNLO loop function wtilde Squared(!) (Eq. 4.20)
  ! impi  : determines which mpi2 to use
  ! chp_itope;  EM or KW, choose which iterated 2pe 'model' to use
  ! this is set in the module header
  FUNCTION chp_NNLO_two_pion_exchange_Vc(q2, w, A, wt2, impi) RESULT(res)
    
    REAL*8 , INTENT(IN) :: q2
    REAL*8 , INTENT(IN) :: w
    REAL*8 , INTENT(IN) :: A
    REAL*8 , INTENT(IN) :: wt2
    INTEGER, INTENT(IN) :: impi
    REAL*8 :: res
    
    res = 0.0D0
    res = const(9)*(const(10)*mpi5(impi)/(mnuc(0)*w*w) - &
         (mpi2(impi)*const(11) - q2*c3 - q2*3.0D0*const(10)/mnuc(0) ) *wt2*A)
    
    ! IF twopion exchange in BbS , i.e. EM format
    ! add correction to the NNLO central term
    IF (chp_itope == 'EM') THEN
       res = res - const(12)*(mpi(impi)*w*w+wt2*wt2*A)/mnuc(0)
    END IF

  END FUNCTION chp_NNLO_two_pion_exchange_Vc

  ! NNLO function Eq [1] 4.14 and 4.22
  ! q2    : momentum transfer squared
  ! w     : Eq 4.12
  ! A     : NNLO loop function A Eq. 4.19
  ! wt2   : NNLO loop function wtilde Squared(!) (Eq. 4.20)
  ! impi  : determines which mpi2 to use
  ! chp_itope;  EM or KW, choose which iterated 2pe 'model' to use
  ! this is set in the module header
  FUNCTION chp_NNLO_two_pion_exchange_Wc(q2, w, A, wt2, impi) RESULT(res)
    
    REAL*8 , INTENT(IN) :: q2
    REAL*8 , INTENT(IN) :: w
    REAL*8 , INTENT(IN) :: A
    REAL*8 , INTENT(IN) :: wt2
    INTEGER, INTENT(IN) :: impi
    REAL*8 :: res
    
    res = 0.0D0
    res = const(13) * (3.0D0*gA2*mpi5(impi)/(w*w) - & 
         (fourmpi2(impi) + 2.0D0*q2 - gA2*(fourmpi2(impi)+3.0D0*q2))*wt2*A)/mnuc(0)
             
    ! IF twopion exchange in BbS , i.e. EM format
    ! add correction to the NNLO central term
    IF (chp_itope == 'EM') THEN
       res = res + const(14)*(mpi(impi)*w*w + wt2*wt2*A)/mnuc(0)
    END IF
    
  END FUNCTION chp_NNLO_two_pion_exchange_Wc
  
  ! NNLO function Eq [1] 4.15 and 4.23
  ! w     : Eq 4.12
  ! A     : NNLO loop function A Eq. 4.19
  ! wt2   : NNLO loop function wtilde Squared(!) (Eq. 4.20)
  ! impi  : determines which mpi2 to use
  ! chp_itope;  EM or KW, choose which iterated 2pe 'model' to use
  ! this is set in the module header
  FUNCTION chp_NNLO_two_pion_exchange_VT(w, A, wt2, impi) RESULT(res)
    
    REAL*8 , INTENT(IN) :: w
    REAL*8 , INTENT(IN) :: A
    REAL*8 , INTENT(IN) :: wt2
    INTEGER, INTENT(IN) :: impi
    REAL*8 :: res
    
    res = 0.0D0
    res = 3.0D0*const(15)*wt2*A/mnuc(0)
             
    ! IF twopion exchange in BbS , i.e. EM format
    ! add correction to the NNLO central term
    IF (chp_itope == 'EM') THEN
       res = res + const(15)*(mpi(impi) + w*w*A )/mnuc(0)
    END IF
    
  END FUNCTION chp_NNLO_two_pion_exchange_VT

  ! NNLO function Eq [1] 4.15 and 4.23
  ! q2    : momentum transfer squared
  ! w     : Eq 4.12
  ! A     : NNLO loop function A Eq. 4.19
  ! wt2   : NNLO loop function wtilde Squared(!) (Eq. 4.20)
  ! impi  : determines which mpi2 to use
  ! chp_itope;  EM or KW, choose which iterated 2pe 'model' to use
  ! this is set in the module header
  FUNCTION chp_NNLO_two_pion_exchange_Vs(q2, w, A, wt2, impi) RESULT(res)
    
    REAL*8 , INTENT(IN) :: q2
    REAL*8 , INTENT(IN) :: w
    REAL*8 , INTENT(IN) :: A
    REAL*8 , INTENT(IN) :: wt2
    INTEGER, INTENT(IN) :: impi
    REAL*8 :: res
    
    res = 0.0D0
    res = -3.0D0*q2*const(15)*wt2*A/mnuc(0)

    ! IF twopion exchange in BbS , i.e. EM format
    ! add correction to the NNLO central term
    IF (chp_itope == 'EM') THEN
       res = res - 1.0D0*const(15)*q2*(mpi(impi) + w*w*A)/mnuc(0)
    END IF
    
    !WRITE(*,*)q2
    !WRITE(*,*)const(15)
    !WRITE(*,*)wt2
    !WRITE(*,*)A
    !WRITE(*,*)mnuc(0)
    !WRITE(*,*)const(15)
    !WRITE(*,*)mpi(impi)
    !WRITE(*,*)w*w
    !WRITE(*,*)res
    
  END FUNCTION chp_NNLO_two_pion_exchange_Vs
  
  ! NNLO function Eq [1] 4.16 and 4.24
  ! q2    : momentum transfer squared
  ! w     : Eq 4.12
  ! A     : NNLO loop function A Eq. 4.19
  ! impi  : determines which mpi2 to use
  ! chp_itope;  EM or KW, choose which iterated 2pe 'model' to use
  ! this is set in the module header
  FUNCTION chp_NNLO_two_pion_exchange_WT(q2, w, A, impi) RESULT(res)
    
    REAL*8 , INTENT(IN) :: q2
    REAL*8 , INTENT(IN) :: w
    REAL*8 , INTENT(IN) :: A
    INTEGER, INTENT(IN) :: impi
    REAL*8 :: res
    
    res = 0.0D0
    res = -1.0D0*const(16)*A*( (c4 + 1.0D0/(4.0D0*mnuc(0)))*w*w - &
         const(17)*(10.0D0*mpi2(impi) + 3.0D0*q2)/mnuc(0))
             
    ! IF twopion exchange in BbS , i.e. EM format
    ! add correction to the NNLO central term
    IF (chp_itope == 'EM') THEN
       res = res - const(18)*(mpi(impi) + w*w*A)/mnuc(0)
    END IF
    
  END FUNCTION chp_NNLO_two_pion_exchange_WT 

  ! NNLO function Eq [1] 4.16 and 4.24
  ! q2    : momentum transfer squared
  ! w     : Eq 4.12
  ! A     : NNLO loop function A Eq. 4.19
  ! impi  : determines which mpi2 to use
  ! chp_itope;  EM or KW, choose which iterated 2pe 'model' to use
  ! this is set in the module header
  FUNCTION chp_NNLO_two_pion_exchange_Ws(q2, w, A, impi) RESULT(res)
    
    REAL*8 , INTENT(IN) :: q2
    REAL*8 , INTENT(IN) :: w
    REAL*8 , INTENT(IN) :: A
    INTEGER, INTENT(IN) :: impi
    REAL*8 :: res
    
    res = 0.0D0
    res = q2*const(16)*A*( (c4+1.0D0/(4.0D0*mnuc(0)))*w*w - & 
         const(17)*(10.0D0*mpi2(impi) +3.0D0*q2)/mnuc(0))
             
    ! IF twopion exchange in BbS , i.e. EM format
    ! add correction to the NNLO central term
    IF (chp_itope == 'EM') THEN
       res = res + const(18)*q2*(mpi(impi) + w*w*A)/mnuc(0)
    END IF
    
  END FUNCTION chp_NNLO_two_pion_exchange_Ws

  ! NNLO function Eq [1] 4.17
  ! A     : NNLO loop function A Eq. 4.19
  ! wt2   : NNLO loop function wtilde Squared(!) (Eq. 4.20)
  ! impi  : determines which mpi2 to use
  FUNCTION chp_NNLO_two_pion_exchange_VLS(A, wt2, impi) RESULT(res)
    
    REAL*8 , INTENT(IN) :: A
    REAL*8 , INTENT(IN) :: wt2
    INTEGER, INTENT(IN) :: impi
    REAL*8 :: res
    
    res = 0.0D0
    res = const(19) * wt2*A/mnuc(0)
             
  END FUNCTION chp_NNLO_two_pion_exchange_VLS

  ! NNLO function Eq [1] 4.18
  ! w   : NNLO loop function wtilde Squared(!) (Eq. 4.20)
  ! A     : NNLO loop function A Eq. 4.19
  ! impi  : determines which mpi2 to use
  FUNCTION chp_NNLO_two_pion_exchange_WLS(w, A, impi) RESULT(res)
    
    REAL*8 , INTENT(IN) :: w
    REAL*8 , INTENT(IN) :: A
    INTEGER, INTENT(IN) :: impi
    REAL*8 :: res
    
    res = 0.0D0
    res = const(20)*w*w*A/mnuc(0)
             
  END FUNCTION chp_NNLO_two_pion_exchange_WLS

  FUNCTION chp_sigma_dot_sigma_mtx(ms1,ms2,ms3,ms4) RESULT(res)
    
    IMPLICIT NONE 
    INTEGER, INTENT(IN) :: ms1,ms2,ms3,ms4
    REAL*8 :: res
    complex*16 :: res1
    COMPLEX*16 :: chi1(2), chi2(2), chi3(2), chi4(2)
    INTEGER :: i1 
    
    res = 0.0D0
    chi1 = 0.d0; chi2 = 0.d0; chi3 = 0.d0; chi4 = 0.d0
    i1 = nint(1.5-0.5*ms1) 
    chi1(i1) = 1.d0 
    i1 = nint(1.5-0.5*ms2) 
    chi2(i1) = 1.d0 
    i1 = nint(1.5-0.5*ms3) 
    chi3(i1) = 1.d0 
    i1 = nint(1.5-0.5*ms4) 
    chi4(i1) = 1.d0 
    
 
    res1 = dot_product(chi1, matmul( sigma_x, chi3)) * dot_product(chi2, matmul( sigma_x, chi4)) &
         + dot_product(chi1, matmul( sigma_y, chi3)) * dot_product(chi2, matmul( sigma_y, chi4)) &
         + dot_product(chi1, matmul( sigma_z, chi3)) * dot_product(chi2, matmul( sigma_z, chi4)) 

    !if ( abs( aimag(res1)) > 1.e-10 ) write(6,*)  'yep', res1
    res = res1
    
    
  END FUNCTION chp_sigma_dot_sigma_mtx

  FUNCTION chp_sigma_dot_sigma(ms1,ms2,ms3,ms4) RESULT(res)
    
    IMPLICIT NONE 
    INTEGER, INTENT(IN) :: ms1,ms2,ms3,ms4
    REAL*8 :: res, spin 
    COMPLEX*16 :: res1, jx1, jx2, jy1, jy2, jz1, jz2 
    REAL*8 :: m1,m2,m3,m4
    INTEGER :: i1 
    
    m1 = 0.5d0*dble(ms1); m2 = 0.5d0*dble(ms2)
    m3 = 0.5d0*dble(ms3); m4 = 0.5d0*dble(ms4)
    
    jx1 = 0.d0
    jx2 = 0.d0 
    jy1 = 0.d0
    jy2 = 0.d0 
    jz1 = 0.d0
    jz2 = 0.d0 
    
    spin = 0.5d0 
    res = 0.0D0
    if ( ms1 == ms3 ) THEN
       jx1 = 0.d0
       jy1 = 0.d0 
       jz1 = 2.d0*m1
    end if
    if ( ms1 == ms3 + 2 ) then 
       jx1 = sqrt( (spin+m3+1.d0)*(spin-m3) )
       jy1 = -dcmplx(0.d0,1.d0)*sqrt( (spin+m3+1.d0)*(spin-m3) )
       jz1 = 0.d0
    end if
    if ( ms1 == ms3 - 2 ) then 
       jx1 = sqrt( (spin-m3+1.d0)*(spin+m3) )
       jy1 = dcmplx(0.d0,1.d0)*sqrt( (spin-m3+1.d0)*(spin+m3) )
       jz1 = 0.d0
    end if
    
    if ( ms2 == ms4 ) THEN
       jx2 = 0.d0
       jy2 = 0.d0 
       jz2 = 2.d0*m2 
    end if
    if ( ms2 == ms4 + 2 ) then 
       jx2 = sqrt( (spin+m4+1.d0)*(spin-m4) )
       jy2 = -dcmplx(0.d0,1.d0)*sqrt( (spin+m4+1.d0)*(spin-m4) )
       jz2 = 0.d0
    end if
    if ( ms2 == ms4 - 2 ) then 
       jx2 = sqrt( (spin-m4+1.d0)*(spin+m4) )
       jy2 = dcmplx(0.d0,1.d0)*sqrt( (spin-m4+1.d0)*(spin+m4) )
       jz2 = 0.d0
    end if
            
    res1 = jx1*jx2 + jy1*jy2 + jz1*jz2
    res = res1
    
  END FUNCTION chp_sigma_dot_sigma

  FUNCTION chp_tau_dot_tau_mtx(ms1,ms2,ms3,ms4) RESULT(res)
    
    IMPLICIT NONE 
    INTEGER, INTENT(IN) :: ms1,ms2,ms3,ms4
    REAL*8 :: res, delta
    complex*16 :: res1
    COMPLEX*16 :: chi1(2), chi2(2), chi3(2), chi4(2)
    INTEGER :: i1 
    
    res = 0.0D0
    chi1 = 0.d0; chi2 = 0.d0; chi3 = 0.d0; chi4 = 0.d0
    i1 = nint(1.5-0.5*ms1) 
    chi1(i1) = 1.d0 
    i1 = nint(1.5-0.5*ms2) 
    chi2(i1) = 1.d0 
    i1 = nint(1.5-0.5*ms3) 
    chi3(i1) = 1.d0 
    i1 = nint(1.5-0.5*ms4) 
    chi4(i1) = 1.d0 
    
    
    res1 = dot_product(chi1, matmul( sigma_x, chi3)) * dot_product(chi2, matmul( sigma_x, chi4)) &
         + dot_product(chi1, matmul( sigma_y, chi3)) * dot_product(chi2, matmul( sigma_y, chi4)) &
         + dot_product(chi1, matmul( sigma_z, chi3)) * dot_product(chi2, matmul( sigma_z, chi4)) 
    res = res1
    
    !if ( abs( aimag(res1)) > 1.e-10 ) write(6,*)  'yep', res1
    
  END FUNCTION chp_tau_dot_tau_mtx
  
  
  FUNCTION chp_tau_dot_tau(ms1,ms2,ms3,ms4) RESULT(res)
    
    IMPLICIT NONE 
    INTEGER, INTENT(IN) :: ms1,ms2,ms3,ms4
    REAL*8 :: res, spin 
    COMPLEX*16 :: res1, jx1, jx2, jy1, jy2, jz1, jz2 
    REAL*8 :: m1,m2,m3,m4
    INTEGER :: i1 
    
    m1 = 0.5d0*dble(ms1); m2 = 0.5d0*dble(ms2)
    m3 = 0.5d0*dble(ms3); m4 = 0.5d0*dble(ms4)
    
    jx1 = 0.d0
    jx2 = 0.d0 
    jy1 = 0.d0
    jy2 = 0.d0 
    jz1 = 0.d0
    jz2 = 0.d0 
    
    spin = 0.5d0 
    res = 0.0D0
    if ( ms1 == ms3 ) THEN
       jx1 = 0.d0
       jy1 = 0.d0 
       jz1 = 2.d0*m1
    end if
    if ( ms1 == ms3 + 2 ) then 
       jx1 = sqrt( (spin+m3+1.d0)*(spin-m3) )
       jy1 = -dcmplx(0.d0,1.d0)*sqrt( (spin+m3+1.d0)*(spin-m3) )
       jz1 = 0.d0
    end if
    if ( ms1 == ms3 - 2 ) then 
       jx1 = sqrt( (spin-m3+1.d0)*(spin+m3) )
       jy1 = dcmplx(0.d0,1.d0)*sqrt( (spin-m3+1.d0)*(spin+m3) )
       jz1 = 0.d0
    end if
    
    if ( ms2 == ms4 ) THEN
       jx2 = 0.d0
       jy2 = 0.d0 
       jz2 = 2.d0*m2 
    end if
    if ( ms2 == ms4 + 2 ) then 
       jx2 = sqrt( (spin+m4+1.d0)*(spin-m4) )
       jy2 = -dcmplx(0.d0,1.d0)*sqrt( (spin+m4+1.d0)*(spin-m4) )
       jz2 = 0.d0
    end if
    if ( ms2 == ms4 - 2 ) then 
       jx2 = sqrt( (spin-m4+1.d0)*(spin+m4) )
       jy2 = dcmplx(0.d0,1.d0)*sqrt( (spin-m4+1.d0)*(spin+m4) )
       jz2 = 0.d0
    end if
            
    res1 = jx1*jx2 + jy1*jy2 + jz1*jz2
    res = res1
    
  END FUNCTION chp_tau_dot_tau
  
  FUNCTION chp_sigma_dot_q_mtx(ms1,ms3,q) RESULT(res)

    IMPLICIT NONE 
    INTEGER, INTENT(IN) :: ms1, ms3
    REAL*8, INTENT(IN) :: q(3) 
    COMPLEX*16 :: res
    COMPLEX*16 :: res1, res2 
    COMPLEX*16 :: chi1(2), chi3(2), mat(2,2) 
    INTEGER :: i1 
    
    chi1 = 0.d0; chi3 = 0.d0
    i1 = nint(1.5-0.5*ms1) 
    chi1(i1) = 1.d0 
    i1 = nint(1.5-0.5*ms3) 
    chi3(i1) = 1.d0 
    
    mat = q(1)*sigma_x+q(2)*sigma_y+q(3)*sigma_z 
    res1 = dot_product(chi1, matmul( mat, chi3))
    res = res1 
    
  end FUNCTION chp_sigma_dot_q_mtx

  
  FUNCTION chp_sigma1_dot_q_sigma2_dot_q_mtx(ms1,ms2,ms3,ms4,q) RESULT(res)
    
    IMPLICIT NONE 
    INTEGER, INTENT(IN) :: ms1,ms2,ms3,ms4
    REAL*8, INTENT(IN) :: q(3) 
    COMPLEX*16 :: res
    COMPLEX*16 :: res1, res2 
    COMPLEX*16 :: chi1(2), chi2(2), chi3(2), chi4(2), mat(2,2) 
    INTEGER :: i1 
    
    chi1 = 0.d0; chi2 = 0.d0; chi3 = 0.d0; chi4 = 0.d0
    i1 = nint(1.5-0.5*ms1) 
    chi1(i1) = 1.d0 
    i1 = nint(1.5-0.5*ms2) 
    chi2(i1) = 1.d0 
    i1 = nint(1.5-0.5*ms3) 
    chi3(i1) = 1.d0 
    i1 = nint(1.5-0.5*ms4) 
    chi4(i1) = 1.d0 
    
    mat = q(1)*sigma_x+q(2)*sigma_y+q(3)*sigma_z 
    res1 = dot_product(chi1, matmul( mat, chi3))
    res2 = dot_product(chi2, matmul( mat, chi4))
    
    res = res1 * res2  
    !if ( abs( aimag(res1*res2) ) > 1.e-10 ) write(6,*) res1*res2 
    

  END FUNCTION chp_sigma1_dot_q_sigma2_dot_q_mtx

  FUNCTION chp_sigma1_dot_q_sigma2_dot_q(ms1,ms2,ms3,ms4,q) RESULT(res)

    IMPLICIT NONE 
    INTEGER, INTENT(IN) :: ms1,ms2,ms3,ms4
    REAL*8, INTENT(IN) :: q(3) 
    COMPLEX*16 :: res
    COMPLEX*16 :: res1,res2, jx1, jx2, jy1, jy2, jz1, jz2 
    REAL*8 :: m1,m2,m3,m4, spin 
    INTEGER :: i1 
    
    m1 = 0.5d0*dble(ms1); m2 = 0.5d0*dble(ms2)
    m3 = 0.5d0*dble(ms3); m4 = 0.5d0*dble(ms4)
    
    jx1 = 0.d0
    jx2 = 0.d0 
    jy1 = 0.d0
    jy2 = 0.d0 
    jz1 = 0.d0
    jz2 = 0.d0 
    
    spin = 0.5d0 
    res = 0.0D0
    if ( ms1 == ms3 ) THEN
       jx1 = 0.d0
       jy1 = 0.d0 
       jz1 = 2.d0*m1
    end if
    if ( ms1 == ms3 + 2 ) then 
       jx1 = sqrt( (spin+m3+1.d0)*(spin-m3) )
       jy1 = -dcmplx(0.d0,1.d0)*sqrt( (spin+m3+1.d0)*(spin-m3) )
       jz1 = 0.d0
    end if
    if ( ms1 == ms3 - 2 ) then 
       jx1 = sqrt( (spin-m3+1.d0)*(spin+m3) )
       jy1 = dcmplx(0.d0,1.d0)*sqrt( (spin-m3+1.d0)*(spin+m3) )
       jz1 = 0.d0
    end if
    
    if ( ms2 == ms4 ) THEN
       jx2 = 0.d0
       jy2 = 0.d0 
       jz2 = 2.d0*m2 
    end if
    if ( ms2 == ms4 + 2 ) then 
       jx2 = sqrt( (spin+m4+1.d0)*(spin-m4) )
       jy2 = -dcmplx(0.d0,1.d0)*sqrt( (spin+m4+1.d0)*(spin-m4) )
       jz2 = 0.d0
    end if
    if ( ms2 == ms4 - 2 ) then 
       jx2 = sqrt( (spin-m4+1.d0)*(spin+m4) )
       jy2 = dcmplx(0.d0,1.d0)*sqrt( (spin-m4+1.d0)*(spin+m4) )
       jz2 = 0.d0
    end if
            
    res1 = jx1*q(1) + jy1*q(2) + jz1*q(3)
    res2 = jx2*q(1) + jy2*q(2) + jz2*q(3)
    res = res1*res2 
    
  end FUNCTION chp_sigma1_dot_q_sigma2_dot_q

  FUNCTION chp_spin_dot_qxk_mtx(ms1,ms2,ms3,ms4,qxk) RESULT(res)

    IMPLICIT NONE 
    INTEGER, INTENT(IN) :: ms1,ms2,ms3,ms4
    REAL*8, INTENT(IN) :: qxk(3)
    real*8 :: delta 
    COMPLEX*16 :: res
    COMPLEX*16 :: res1, res2 
    COMPLEX*16 :: chi1(2), chi2(2), chi3(2), chi4(2), mat(2,2) 
    INTEGER :: i1 
    
    chi1 = 0.d0; chi2 = 0.d0; chi3 = 0.d0; chi4 = 0.d0
    i1 = nint(1.5-0.5*ms1) 
    chi1(i1) = 1.d0 
    i1 = nint(1.5-0.5*ms2) 
    chi2(i1) = 1.d0 
    i1 = nint(1.5-0.5*ms3) 
    chi3(i1) = 1.d0 
    i1 = nint(1.5-0.5*ms4) 
    chi4(i1) = 1.d0 
    
    mat =-dcmplx(0.d0,1.d0) * 0.5D0*(qxk(1)*sigma_x+qxk(2)*sigma_y+qxk(3)*sigma_z)
    
    res1 = dot_product(chi1, matmul( mat, chi3))*delta(ms2,ms4)
    res2 = dot_product(chi2, matmul( mat, chi4))*delta(ms1,ms3)
    
    res = (res1+res2)
    
    
  END FUNCTION chp_spin_dot_qxk_mtx
  
  
  FUNCTION chp_spin_dot_qxk(ms1,ms2,ms3,ms4,qxk) RESULT(res)

    IMPLICIT NONE 

    INTEGER, INTENT(IN) :: ms1,ms2,ms3,ms4
    REAL*8, INTENT(IN) :: qxk(3) 
    REAL*8 ::  spin, delta, d1, d2
    COMPLEX*16 :: res, res1,res2, jx1, jx2, jy1, jy2, jz1, jz2 
    REAL*8 :: m1,m2,m3,m4
    INTEGER :: i1 
    
    m1 = 0.5d0*dble(ms1); m2 = 0.5d0*dble(ms2)
    m3 = 0.5d0*dble(ms3); m4 = 0.5d0*dble(ms4)
    
    jx1 = 0.d0
    jx2 = 0.d0 
    jy1 = 0.d0
    jy2 = 0.d0 
    jz1 = 0.d0
    jz2 = 0.d0 
    
    spin = 0.5d0 
    res = 0.0D0
    if ( ms1 == ms3 ) THEN
       jx1 = 0.d0
       jy1 = 0.d0 
       jz1 = 2.d0*m1
    end if
    if ( ms1 == ms3 + 2 ) then 
       jx1 = sqrt( (spin+m3+1.d0)*(spin-m3) )
       jy1 = -dcmplx(0.d0,1.d0)*sqrt( (spin+m3+1.d0)*(spin-m3) )
       jz1 = 0.d0
    end if
    if ( ms1 == ms3 - 2 ) then 
       jx1 = sqrt( (spin-m3+1.d0)*(spin+m3) )
       jy1 = dcmplx(0.d0,1.d0)*sqrt( (spin-m3+1.d0)*(spin+m3) )
       jz1 = 0.d0
    end if
    
    if ( ms2 == ms4 ) THEN
       jx2 = 0.d0
       jy2 = 0.d0 
       jz2 = 2.d0*m2 
    end if
    if ( ms2 == ms4 + 2 ) then 
       jx2 = sqrt( (spin+m4+1.d0)*(spin-m4) )
       jy2 = -dcmplx(0.d0,1.d0)*sqrt( (spin+m4+1.d0)*(spin-m4) )
       jz2 = 0.d0
    end if
    if ( ms2 == ms4 - 2 ) then 
       jx2 = sqrt( (spin-m4+1.d0)*(spin+m4) )
       jy2 = dcmplx(0.d0,1.d0)*sqrt( (spin-m4+1.d0)*(spin+m4) )
       jz2 = 0.d0
    end if
    
    d1 = delta(ms1,ms3)
    d2 = delta(ms2,ms4)
    res1 = -dcmplx(0.d0,0.5d0) * ( (jx1*d2+jx2*d1)*qxk(1) + (jy1*d2+jy2*d1)*qxk(2) + (jz1*d2+jz2*d1)*qxk(3) )
    res = res1 
    
    
  END FUNCTION chp_spin_dot_qxk
  
  
  ! regulator, eq 4.63 
  ! pfinal : final momentum
  ! pinit  : initial momentum
  ! n      : cutoff order 
  ! LAMBDA is accessed from the chp constant chp_lambda
  FUNCTION freg(pfinal, pinit, n) RESULT(res)
    
    REAL*8 , INTENT(IN) :: pfinal(3), pinit(3)
    REAL*8 , INTENT(IN) :: n
    REAL*8 :: res,exponent, p2, pp2 
    
    res = 0.0D0
    p2 = sum( pfinal*pfinal) 
    pp2 = sum( pinit*pinit )
    
    exponent = (p2**n/lambda**(2.d0*n) + &
         pp2**n/lambda**(2.0D0*n) )
    
    res = dexp(-exponent)
    
  END FUNCTION freg
  

  

end module chiral_constants



module chiral_potentials 

  implicit none 
  
  TYPE, PUBLIC :: chp_int_type
     INTEGER          :: val
     CHARACTER(LEN=12):: name
     LOGICAL          :: set
  END TYPE chp_int_type

  TYPE, PUBLIC :: chp_real_type
     REAL*8           :: val
     CHARACTER(LEN=12):: name
     LOGICAL          :: set
  END TYPE chp_real_type

  TYPE(chp_int_type), PRIVATE :: chp_chiral_order ! LO, NLO, NNLO, N3LO
  TYPE(chp_real_type), PRIVATE :: chp_regcut, chp_regcut_nnlo
  
  
contains 
  
  complex*16 function chiral_pot(p,q,r,s)
    USE single_particle_orbits
    USE constants
    use chiral_constants
    use ang_mom_functions, only : tjs 
    
    implicit none 
    INTEGER, INTENT(IN) :: p,q,r,s 
    INTEGER :: m1,m2,m3,m4, spin, iph, t1,t2,t3,t4, Tiso
    INTEGER :: nx1, ny1, nz1, nx2, ny2, nz2, nx3, ny3, nz3, nx4, ny4, nz4, mt1, mt2, mt3, mt4 
    REAL*8 :: k1(3), k2(3), k3(3), k4(3), kmean(3) 
    REAL*8 :: qtrans(3), prel(3), pprel(3), qxk(3)
    REAL*8 :: q2, p2, kmean2, qabs, pp2, cg1, cg2 
    REAL*8 :: delta, nucleon_mass, relativity_factor
    
    COMPLEX*16 :: vlo, vnlo, vnnlo, vn3lo, vdir
    COMPLEX*16 :: cont_lo, cont_nlo, cont_n3lo, v1pe_tiso(0:1)
    ! LIST OF OPERATOR STRUCTURES
    COMPLEX*16 :: Vc
    COMPLEX*16 :: Wc
    COMPLEX*16 :: Vs
    COMPLEX*16 :: Ws
    COMPLEX*16 :: VLS
    COMPLEX*16 :: WLS
    COMPLEX*16 :: VsigL
    COMPLEX*16 :: WsigL
    COMPLEX*16 :: VT
    COMPLEX*16 :: WT
    COMPLEX*16 :: Vsigk
    COMPLEX*16 :: Wsigk
    COMPLEX*16 :: term_CS
    COMPLEX*16 :: term_CT
    COMPLEX*16 :: term_C(1:7)
    REAL*8 :: loop_w
    REAL*8 :: loop_wtilde2
    REAL*8 :: loop_s
    REAL*8 :: loop_L
    REAL*8 :: loop_A
    
    chp_chiral_order%val = NNLO 
    
    Vc      = 0.0D0
    Wc      = 0.0D0
    Vs      = 0.0D0
    Ws      = 0.0D0
    VLS     = 0.0D0
    WLS     = 0.0D0
    VsigL   = 0.0D0
    WsigL   = 0.0D0
    VT      = 0.0D0
    WT      = 0.0D0
    Vsigk   = 0.0D0
    Wsigk   = 0.0D0
  
    term_C = 0.0D0

    vlo = 0.0D0 ; vnlo = 0.0D0 ; vnnlo = 0.0D0 ; vn3lo = 0.0D0
    
    cont_lo = 0.0D0 ; cont_nlo = 0.0D0 ; cont_n3lo = 0.0D0

    nx1 = all_orbit%nx(p)
    ny1 = all_orbit%ny(p)
    nz1 = all_orbit%nz(p)
    nx2 = all_orbit%nx(q)
    ny2 = all_orbit%ny(q)
    nz2 = all_orbit%nz(q)
    nx3 = all_orbit%nx(r)
    ny3 = all_orbit%ny(r)
    nz3 = all_orbit%nz(r)
    nx4 = all_orbit%nx(s)
    ny4 = all_orbit%ny(s)
    nz4 = all_orbit%nz(s)
    ! 
    ! Conservation of linear momentum
    !
    if ( nx1 + nx2 /= nx3 + nx4 ) return 
    if ( ny1 + ny2 /= ny3 + ny4 ) return 
    if ( nz1 + nz2 /= nz3 + nz4 ) return 
    
    k1(1) = all_orbit%kx(p)
    k1(2) = all_orbit%ky(p)
    k1(3) = all_orbit%kz(p)
    k2(1) = all_orbit%kx(q)
    k2(2) = all_orbit%ky(q)
    k2(3) = all_orbit%kz(q)
    k3(1) = all_orbit%kx(r)
    k3(2) = all_orbit%ky(r)
    k3(3) = all_orbit%kz(r)
    k4(1) = all_orbit%kx(s)
    k4(2) = all_orbit%ky(s)
    k4(3) = all_orbit%kz(s)
    
    ! momenta in MeV
    k1 = k1*hbarc
    k2 = k2*hbarc
    k3 = k3*hbarc
    k4 = k4*hbarc

    ! 
    ! conservation of isospin 
    !
    if ( all_orbit%itzp(p) + all_orbit%itzp(q) /= all_orbit%itzp(r) + all_orbit%itzp(s) ) return 
  
    m1 = all_orbit%szp(p) 
    m2 = all_orbit%szp(q) 
    m3 = all_orbit%szp(r) 
    m4 = all_orbit%szp(s) 
  
    t1 = all_orbit%itzp(p) 
    t2 = all_orbit%itzp(q) 
    t3 = all_orbit%itzp(r) 
    t4 = all_orbit%itzp(s) 
  
    ! 
    ! RELATIVE MOMENTA <prel |v| pprel > 
    ! 
    prel  = 0.5d0*(k1-k2)
    pprel = 0.5d0*(k3-k4)
    p2  = sum(prel*prel) 
    pp2 = sum(pprel*pprel)
    !
    ! AVERAGE MOMENTA kav = 1/2(prel+pprel)
    !
    kmean = 0.5d0*( prel + pprel ) 
    kmean2 = sum(kmean*kmean)
    
    nucleon_mass = mnuc((t1+t2)/2)
    relativity_factor = nucleon_mass/ & 
         sqrt( sqrt(( nucleon_mass**2 + p2) * (nucleon_mass**2 + pp2 )) )
    
    !
    ! momentum transfer 
    !
    qtrans = prel - pprel
    q2 = sum(qtrans*qtrans) 
    qabs = dsqrt(q2)
    
    !
    !  cross product between momentum transfer and average momenta q X k 
    !
    qxk(1) = qtrans(2)*kmean(3)-qtrans(3)*kmean(2) 
    qxk(2) = qtrans(3)*kmean(1)-qtrans(1)*kmean(3) 
    qxk(3) = qtrans(1)*kmean(2)-qtrans(2)*kmean(1) 
    
    IF (chp_chiral_order%val == LO ) THEN
       
       ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
       !
       ! LEADING ORDER 
       ! 
       ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
       
       chp_regcut%val = 2.0D0
       
       !
       ! leading order one-pion exchange 
       ! 
       WT = chp_sigma1_dot_q_sigma2_dot_q_mtx(m1,m2,m3,m4,qtrans)*&
            chp_one_pion_exchange(q2, 2)*chp_tau_dot_tau_mtx(t1,t2,t3,t4)
       
       vlo = WT
       
       
       !
       ! leading order contacts 
       !
       term_CS = CS(0)*delta(m1,m3)*delta(m2,m4)*delta(t1,t3)*delta(t2,t4)
       term_CT = CT(0)*chp_sigma_dot_sigma_mtx(m1,m2,m3,m4)*delta(t1,t3)*delta(t2,t4) 
       cont_lo = term_CS + term_CT
       
    END IF
    
    IF (chp_chiral_order%val >= NLO) THEN
       
       ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
       !
       ! NEXT TO LEADING ORDER 
       ! 
       ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
       
       chp_regcut%val = 2.0D0
       
       !
       ! NLO pion-exchange 
       ! 
       ! irreducible, or non-polynomial, NLO two-pion exchanges
       
       loop_w = chp_NLO_two_pion_exchange_loop_w(q2,2)
       loop_s = chp_NLO_sfr_two_pion_exchange_loop_s(2)
       loop_L = chp_NLO_sfr_two_pion_exchange_loop_L(qabs, q2, loop_w, loop_s, 2)
       
       ! LO CIB one-pion exchanges
       ! [1] Eq. 4.77-4.79
       
       !
       ! Note: the LO contribution is modified at NLO
       ! to include CIB effects. Therefore, should you go beyond
       ! LO, make sure to not include the vlo contribution 
       ! in the sum defining vdir.
       !
       
       
       ! 
       ! Below tjs is the 3-j symbol and is given in the module "ang_mom_functions" 
       ! You need to call "commons_to_angmom" to initialize the 3-j, 6-j symboles etc. 
       ! 
       WT = 0.d0 
       IF (t1 + t2 == 0) THEN 
          
          WT = 0.d0 
          DO Tiso = 0, 2, 2
             cg1 = tjs(1,1,Tiso,t1,t2,-(t1+t2))*iph( (t1+t2)/2 )*sqrt(Tiso+1.d0)
             cg2 = tjs(1,1,Tiso,t3,t4,-(t3+t4))*iph( (t3+t4)/2 )*sqrt(Tiso+1.d0)
             
             WT = WT + cg1*cg2*(-chp_one_pion_exchange(q2, 0) + iph(Tiso/2+1)*2.0D0*chp_one_pion_exchange(q2,1)) 
          end DO
          

       END IF
       IF (t1 + t2 /= 0) THEN
          WT = chp_one_pion_exchange(q2, 0) 
       END IF
       vlo = WT*chp_sigma1_dot_q_sigma2_dot_q_mtx(m1,m2,m3,m4,qtrans)
       
       
       !
       ! leading order CIB contacts 
       !
       term_CS = CS((t1+t2)/2)*delta(m1,m3)*delta(m2,m4)*delta(t1,t3)*delta(t2,t4)
       term_CT = CT((t1+t2)/2) * chp_sigma_dot_sigma_mtx(m1,m2,m3,m4) *delta(t1,t3)*delta(t2,t4) 
       cont_lo = term_CS + term_CT
       !cont_lo = 4.*cont_lo
       
       Wc = chp_tau_dot_tau_mtx(t1,t2,t3,t4)*delta(m1,m3)*delta(m2,m4)* & 
            chp_NLO_two_pion_exchange_Wc(q2, loop_L, loop_w, 2)

       Vs = chp_NLO_two_pion_exchange_Vs(q2, loop_L, 2)* &
            chp_sigma_dot_sigma_mtx(m1,m2,m3,m4)*delta(t1,t3)*delta(t2,t4) 

       VT = chp_NLO_two_pion_exchange_VT(loop_L, 2)* & 
            chp_sigma1_dot_q_sigma2_dot_q_mtx(m1,m2,m3,m4,qtrans)*delta(t1,t3)*delta(t2,t4) 
       
       vnlo = ( WT + WC + VS + VT ) 
       
       !
       ! next-to-leading order contacts 
       !
       term_C(1) = cnlo(1)*q2*delta(m1,m3)*delta(m2,m4)*delta(t1,t3)*delta(t2,t4) 

       term_C(2) = cnlo(2)*kmean2*delta(m1,m3)*delta(m2,m4)*delta(t1,t3)*delta(t2,t4)

       term_C(3) = cnlo(3)*q2*chp_sigma_dot_sigma_mtx(m1,m2,m3,m4)*delta(t1,t3)*delta(t2,t4)

       term_C(4) = cnlo(4)*kmean2*chp_sigma_dot_sigma_mtx(m1,m2,m3,m4)*delta(t1,t3)*delta(t2,t4)

       term_C(5) = cnlo(5)*chp_spin_dot_qxk_mtx(m1,m2,m3,m4,qxk)*delta(t1,t3)*delta(t2,t4)
       
       term_C(6) = cnlo(6)*chp_sigma1_dot_q_sigma2_dot_q_mtx(m1,m2,m3,m4,qtrans)*delta(t1,t3)*delta(t2,t4)

       term_C(7) = cnlo(7)*chp_sigma1_dot_q_sigma2_dot_q_mtx(m1,m2,m3,m4,kmean)*delta(t1,t3)*delta(t2,t4) 
       
       cont_nlo = SUM(term_C)

    END IF
    
    IF (chp_chiral_order%val >= NNLO) THEN
       
       ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
       !
       ! NEXT TO NEXT TO LEADING ORDER 
       ! 
       ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
       
       chp_regcut_nnlo%val = 3.0D0
       
       ! NNLO loop functions
       
       loop_wtilde2 = chp_NNLO_two_pion_exchange_loop_wtilde2(q2,2)
       loop_A = chp_NNLO_sfr_two_pion_exchange_loop_A(qabs, q2, 2)
       
       ! chp_itope type is set in the header of this function
       Vc  = chp_NNLO_two_pion_exchange_Vc (q2, loop_w, loop_A, loop_wtilde2, 2)*&
             delta(m1,m3)*delta(m2,m4)*delta(t1,t3)*delta(t2,t4)

       Wc  = chp_NNLO_two_pion_exchange_Wc (q2, loop_w, loop_A, loop_wtilde2, 2)*&
             chp_tau_dot_tau_mtx(t1,t2,t3,t4)*delta(m1,m3)*delta(m2,m4)

       VLS = chp_NNLO_two_pion_exchange_VLS(loop_A, loop_wtilde2, 2)*&
            chp_spin_dot_qxk_mtx(m1,m2,m3,m4,qxk)*delta(t1,t3)*delta(t2,t4)
       
       WLS = chp_NNLO_two_pion_exchange_WLS(loop_w, loop_A, 2)*&
            chp_spin_dot_qxk_mtx(m1,m2,m3,m4,qxk)*chp_tau_dot_tau_mtx(t1,t2,t3,t4)

       VT  = chp_NNLO_two_pion_exchange_VT (loop_w, loop_A, loop_wtilde2, 2)*&
             chp_sigma1_dot_q_sigma2_dot_q_mtx(m1,m2,m3,m4,qtrans)*delta(t1,t3)*delta(t2,t4) 

       WT  = chp_NNLO_two_pion_exchange_WT (q2, loop_w, loop_A, 2)*&
             chp_sigma1_dot_q_sigma2_dot_q_mtx(m1,m2,m3,m4,qtrans)*chp_tau_dot_tau_mtx(t1,t2,t3,t4)

       Vs  = chp_NNLO_two_pion_exchange_Vs (q2, loop_w, loop_A, loop_wtilde2, 2)*&
             chp_sigma_dot_sigma_mtx(m1,m2,m3,m4)*delta(t1,t3)*delta(t2,t4) 

       Ws  = chp_NNLO_two_pion_exchange_Ws (q2, loop_w, loop_A, 2)*&
             chp_sigma_dot_sigma_mtx(m1,m2,m3,m4)*chp_tau_dot_tau_mtx(t1,t2,t3,t4)

       vnnlo = (Vc + Wc + VLS + WLS + VT + WT + Vs + Ws)
       !vnnlo = ( VLS + WLS )
       
    END IF
    
    IF (chp_chiral_order%val >= N3LO) THEN
       
       ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
       !
       ! NEXT TO NEXT TO NEXT TO LEADING ORDER 
       ! 
       ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

       chp_regcut_nnlo%val = 3.0D0
       
       vn3lo = 0.0D0
       cont_n3lo = 0.0D0

    END IF
    
    ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    !
    ! DONE WITH CHIRAL ORDER CONTRIBUTIONS
    !
    ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    
    ! sum up all orders 
    ! regulator and relativity factor 
    IF (chp_chiral_order%val >= NNLO) THEN
       vdir = (vlo + cont_lo + vnlo + cont_nlo + vnnlo + cont_n3lo + vn3lo) * &
           relativity_factor * freg(prel, pprel, chp_regcut_nnlo%val)
    ELSE
       vdir = (vlo + cont_lo + vnlo + cont_nlo + vnnlo + cont_n3lo + vn3lo) * &
            relativity_factor * freg(prel, pprel, chp_regcut%val)
    ENDIF
        
    chiral_pot =  hbarc**3 * (vdir)/volume 
    
        
  end function chiral_pot


  

end module chiral_potentials

