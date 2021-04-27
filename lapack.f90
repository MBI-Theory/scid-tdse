!
!   SCID-TDSE: Simple 1-electron atomic TDSE solver
!   Copyright (C) 2015-2021 Serguei Patchkovskii, Serguei.Patchkovskii@mbi-berlin.de
!
!   This program is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with this program.  If not, see <https://www.gnu.org/licenses/>.
!
module lapack
!
!  Simplistic type-agnostic LAPACK and LINPACK interface
!
  use accuracy
  implicit none

  character(len=clen), save :: rcsid_lapack = "$Id: lapack.f90,v 1.10 2021/04/26 15:44:44 ps Exp ps $"

  interface lapack_gelss
     module procedure lapack_cgelss
     module procedure lapack_zgelss
     module procedure lapack_sgelss
     module procedure lapack_dgelss
!*lq module procedure lapack_quad_dgelss
!*lq module procedure lapack_quad_zgelss
  end interface lapack_gelss

  interface lapack_stev
     module procedure lapack_sstev
     module procedure lapack_dstev
  end interface lapack_stev

  interface lapack_sterf
     module procedure lapack_dsterf
     module procedure lapack_ssterf
  end interface lapack_sterf

  interface lapack_geev
     module procedure lapack_cgeev
     module procedure lapack_zgeev
     module procedure lapack_cgeev2
     module procedure lapack_zgeev2
!*lq module procedure lapack_quad_zgeev2
  end interface lapack_geev

  interface lapack_heev
     module procedure lapack_cheev
     module procedure lapack_zheev
  end interface lapack_heev

  interface lapack_syev
     module procedure lapack_dsyev
     module procedure lapack_ssyev
!*lq module procedure lapack_quad_dsyev
  end interface lapack_syev

  interface lapack_ginverse
     module procedure lapack_ginverse_real
     module procedure lapack_ginverse_double
!*lq module procedure lapack_ginverse_quad
     module procedure lapack_ginverse_complex
     module procedure lapack_ginverse_doublecomplex
  end interface lapack_ginverse

  interface lapack_svd
     module procedure lapack_dgesvd
     module procedure lapack_zgesvd
  end interface lapack_svd  

  interface lapack_gesv
     module procedure lapack_dgesv
     module procedure lapack_cgesv
     module procedure lapack_zgesv
!*lq module procedure lapack_quad_zgesv
!*lq module procedure lapack_quad_dgesv
  end interface lapack_gesv

  interface linpack_determinant
     module procedure linpack_determinant_double
  end interface linpack_determinant

  interface linpack_determinant_trash_input
     module procedure linpack_determinant_double_trash_input
  end interface linpack_determinant_trash_input

  interface lapack_determinant
     module procedure lapack_determinant_double
!*lq module procedure lapack_determinant_quad
  end interface lapack_determinant

  contains

  subroutine lapack_cgelss(a,b)
    complex(srk), intent(inout) :: a(:,:)
    complex(srk), intent(inout) :: b(:,:)

    external cgelss
    real(srk)    :: s    (   min(size(a,dim=1),size(a,dim=2)))
    complex(srk) :: work (50*max(size(a,dim=1),size(a,dim=2)))
    real(srk)    :: rwork( 5*min(size(a,dim=1),size(a,dim=2)))
    integer      :: rank, info           ! Must be of default integer kind
    integer      :: na1, na2, nb1, nb2   ! Must be of default integer kind
    
    na1 = size(a,dim=1) ; na2 = size(a,dim=2)
    nb1 = size(b,dim=1) ; nb2 = size(b,dim=2)
    call cgelss(na1,na2,nb2,a(1:na1,1:na2),na1,b(1:nb1,1:nb2),nb1, &
                s, 100.0*spacing(1.0), rank, work, size(work), rwork, info)

    if (info/=0) then
      write (out,"(' cgelss returned ',i8)") info
      stop 'lapack_cgelss - cgelss failed'
    end if
  end subroutine lapack_cgelss

  subroutine lapack_zgelss(a,b)
    complex(kind=drk), intent(inout) :: a(:,:)
    complex(kind=drk), intent(inout) :: b(:,:)

    external zgelss
    real(drk)         :: s    (   min(size(a,dim=1),size(a,dim=2)))
    complex(kind=drk) :: work (50*max(size(a,dim=1),size(a,dim=2)))
    real(drk)         :: rwork( 5*min(size(a,dim=1),size(a,dim=2)))
    integer           :: rank, info         ! Must be of default integer kind
    integer           :: na1, na2, nb1, nb2 ! Must be of default integer kind
    
    na1 = size(a,dim=1) ; na2 = size(a,dim=2)
    nb1 = size(b,dim=1) ; nb2 = size(b,dim=2)
    call zgelss(na1,na2,nb2,a(1:na1,1:na2),na1,b(1:nb1,1:nb2),nb1, &
                s, 100.0d0*spacing(1.0d0), rank, work, size(work), rwork, info)

    if (info/=0) then
      write (out,"(' cgelss returned ',i8)") info
      stop 'lapack_cgelss - cgelss failed'
    end if
  end subroutine lapack_zgelss

!*lq subroutine lapack_quad_zgelss(a,b)
!*lq  complex(qrk), intent(inout) :: a(:,:)
!*lq  complex(qrk), intent(inout) :: b(:,:)
!*lq  
!*lq  external quad_zgelss
!*lq  real(qrk)    :: s    (   min(size(a,dim=1),size(a,dim=2)))
!*lq  complex(qrk) :: work (50*max(size(a,dim=1),size(a,dim=2)))
!*lq  real(qrk)    :: rwork( 5*min(size(a,dim=1),size(a,dim=2)))
!*lq  real(qrk)    :: eps
!*lq  integer      :: rank, info          ! Must be of default integer kind
!*lq  integer      :: na1, na2, nb1, nb2  ! Must be of default integer kind
!*lq 
!*lq  na1 = size(a,dim=1) ; na2 = size(a,dim=2)
!*lq  nb1 = size(b,dim=1) ; nb2 = size(b,dim=2)
!*lq  eps = 100*spacing(1._qrk)
!*lq  call quad_zgelss(na1,na2,nb2,a(1:na1,1:na2),na1,b(1:nb1,1:nb2),nb1, &
!*lq                   s, eps, rank, work, size(work), rwork, info)
!*lq  
!*lq  if (info/=0) then
!*lq    write (out,"(' quad_cgelss returned ',i8)") info
!*lq    stop 'lapack_quad_cgelss - quad_cgelss failed'
!*lq  end if
!*lq end subroutine lapack_quad_zgelss

  subroutine lapack_sgelss(a,b)
    real(srk), intent(inout) :: a(:,:)
    real(srk), intent(inout) :: b(:,:)

    external sgelss
    real(srk) :: s    (   min(size(a,dim=1),size(a,dim=2)))
    real(srk) :: work (50*max(size(a,dim=1),size(a,dim=2),size(b,dim=2)))
    integer   :: rank, info          ! Must be of default integer kind
    integer   :: na1, na2, nb1, nb2  ! Must be of default integer kind
    
    na1 = size(a,dim=1) ; na2 = size(a,dim=2)
    nb1 = size(b,dim=1) ; nb2 = size(b,dim=2)
    call sgelss(na1,na2,nb2,a(1:na1,1:na2),na1,b(1:nb1,1:nb2),nb1, &
                s, 100.0*spacing(1.0), rank, work, size(work), info)

    if (info/=0) then
      write (out,"(' sgelss returned ',i8)") info
      stop 'lapack_sgelss - sgelss failed'
    end if
  end subroutine lapack_sgelss

  subroutine lapack_dgelss(a,b)
    real(drk), intent(inout) :: a(:,:)
    real(drk), intent(inout) :: b(:,:)

    external dgelss
    real(drk) :: s    (   min(size(a,dim=1),size(a,dim=2)))
    real(drk) :: work (50*max(size(a,dim=1),size(a,dim=2),size(b,dim=2)))
    integer   :: rank, info         ! Must be of default integer kind
    integer   :: na1, na2, nb1, nb2 ! Must be of default integer kind
    
    na1 = size(a,dim=1) ; na2 = size(a,dim=2)
    nb1 = size(b,dim=1) ; nb2 = size(b,dim=2)
    call dgelss(na1,na2,nb2,a(1:na1,1:na2),na1,b(1:nb1,1:nb2),nb1, &
                s, 100.0d0*spacing(1.0d0), rank, work, size(work), info)

    if (info/=0) then
      write (out,"(' dgelss returned ',i8)") info
      stop 'lapack_dgelss - dgelss failed'
    end if
  end subroutine lapack_dgelss

!*lq subroutine lapack_quad_dgelss(a,b)
!*lq  real(qrk), intent(inout) :: a(:,:)
!*lq  real(qrk), intent(inout) :: b(:,:)
  
!*lq  external quad_dgelss
!*lq  real(qrk) :: s    (   min(size(a,dim=1),size(a,dim=2)))
!*lq  real(qrk) :: work (50*max(size(a,dim=1),size(a,dim=2),size(b,dim=2)))
!*lq  integer   :: rank, info         ! Must be of default integer kind
!*lq  integer   :: na1, na2, nb1, nb2 ! Must be of default integer kind
!*lq  
!*lq  na1 = size(a,dim=1) ; na2 = size(a,dim=2)
!*lq  nb1 = size(b,dim=1) ; nb2 = size(b,dim=2)
!*lq  call quad_dgelss(na1,na2,nb2,a(1:na1,1:na2),na1,b(1:nb1,1:nb2),nb1, &
!*lq                   s, 100._qrk*spacing(1.0_qrk), rank, work, size(work), info)

!*lq  if (info/=0) then
!*lq    write (out,"(' quad_dgelss returned ',i8)") info
!*lq    stop 'lapack_quad_dgelss - quad_dgelss failed'
!*lq  end if
!*lq end subroutine lapack_quad_dgelss

  subroutine lapack_sstev(d,e,z)
    real(srk), intent(inout) :: d(:)   ! In:  Diagonal elements of the matrix 
                                       ! Out: Eigenvalues, ascending order
    real(srk), intent(inout) :: e(:)   ! In:  Sub-/super-diagonal elements of the matrix
                                       ! Out: Destroyed
    real(srk), intent(out)   :: z(:,:) ! Out: Eigenvectors

    real(srk) :: work(max(1,2*size(d)-2))
    integer   :: info      ! Must be of default integer kind
    integer   :: nz1, nz2  ! Must be of default integer kind
    
    nz1 = size(z,dim=1) ; nz2 = size(z,dim=2)
    call sstev('V',size(d),d,e,z(1:nz1,1:nz2),nz1,work,info)

    if (info/=0) then
      write (out,"(' sstev returned ',i8)") info
      stop 'lapack_sstev - sstev failed'
    end if
  end subroutine lapack_sstev

  subroutine lapack_dstev(d,e,z)
    real(drk), intent(inout) :: d(:)   ! In:  Diagonal elements of the matrix 
                                              ! Out: Eigenvalues, ascending order
    real(drk), intent(inout) :: e(:)   ! In:  Sub-/super-diagonal elements of the matrix
                                              ! Out: Destroyed
    real(drk), intent(out)   :: z(:,:) ! Out: Eigenvectors

    real(drk) :: work(max(1,2*size(d)-2))
    integer   :: info      ! Must be of default integer kind
    integer   :: nz1, nz2  ! Must be of default integer kind
    
    nz1 = size(z,dim=1) ; nz2 = size(z,dim=2)
    call dstev('V',size(d),d,e,z(1:nz1,1:nz2),nz1,work,info)

    if (info/=0) then
      write (out,"(' dstev returned ',i8)") info
      stop 'lapack_dstev - dstev failed'
    end if
  end subroutine lapack_dstev

  subroutine lapack_cgeev(h,e)
    complex(srk), intent(inout) :: h(:,:)  ! In:  Hermitian matrix to be diagonalized
                                           ! Out: Eigenvectors
    complex(srk), intent(out)   :: e(:)    ! Out: Eigenvalues

    complex(srk) :: work(50*size(h,dim=2))
    real(srk)    :: rwork(3*size(h,dim=2))
    complex(srk) :: vl(1,1)
    complex(srk) :: vr(size(h,dim=2),size(h,dim=2))
    integer      :: info      ! Must be of default integer kind
    integer      :: nh1, nh2  ! Must be of default integer kind
    
    nh1 = size(h,dim=1) ; nh2 = size(h,dim=2)
    call cgeev('N','V',nh2,h(1:nh1,1:nh2),nh1,e(:),vl,1,vr,nh2,work,size(work),rwork,info)
    if (info/=0) then
      write (out,"(' cgeev returned ',i8)") info
      stop 'lapack_cgeev - cgeev failed'
    end if
    h(1:nh2,1:nh2) = vr
  end subroutine lapack_cgeev

  subroutine lapack_zgeev(h,e)
    complex(kind=drk), intent(inout) :: h(:,:)  ! In:  Hermitian matrix to be diagonalized
                                             ! Out: Eigenvectors
    complex(kind=drk), intent(out)   :: e(:)    ! Out: Eigenvalues

    complex(kind=drk) :: work(50*size(h,dim=2))
    real(drk)         :: rwork(3*size(h,dim=2))
    complex(kind=drk) :: vl(1,1)
    complex(kind=drk) :: vr(size(h,dim=2),size(h,dim=2))
    integer           :: info     ! Must be of default integer kind
    integer           :: nh1, nh2 ! Must be of default integer kind
    
    nh1 = size(h,dim=1) ; nh2 = size(h,dim=2)
    call zgeev('N','V',nh2,h(1:nh1,1:nh2),nh1,e(:),vl,1,vr,nh2,work,size(work),rwork,info)
    if (info/=0) then
      write (out,"(' zgeev returned ',i8)") info
      stop 'lapack_zgeev - zgeev failed'
    end if
    h(1:nh2,1:nh2) = vr
  end subroutine lapack_zgeev

  subroutine lapack_cgeev2(h,e)
    complex(kind=srk), intent(inout) :: h(:,:,:) ! In:  h(:,:,1) = hermitian matrix to be diagonalized
                                                 ! Out: h(:,:,1) = Left eigenvectors
                                                 !      h(:,:,2) = Right eigenvectors
    complex(kind=srk), intent(out)   :: e(:)     ! Out: Eigenvalues

    complex(kind=srk) :: work(50*size(h,dim=2))
    real(srk)         :: rwork(3*size(h,dim=2))
    complex(kind=srk) :: vl(size(h,dim=2),size(h,dim=2))
    integer           :: info     ! Must be of default integer kind
    integer           :: nh1, nh2 ! Must be of default integer kind
    
    nh1 = size(h,dim=1) ; nh2 = size(h,dim=2)
    if (nh1<nh2) stop 'lapack_cgeev2 - oops'
    call cgeev('V','V',nh2,h(1:nh1,1:nh2,1),nh1,e(:),vl,nh2,h(1:nh1,1:nh2,2),nh1,work,size(work),rwork,info)
    if (info/=0) then
      write (out,"(' zgeev returned ',i8)") info
      stop 'lapack_cgeev2 - zgeev failed'
    end if
    h(1:nh2,1:nh2,1) = vl
  end subroutine lapack_cgeev2

  subroutine lapack_zgeev2(h,e)
    complex(kind=drk), intent(inout) :: h(:,:,:) ! In:  h(:,:,1) = hermitian matrix to be diagonalized
                                                 ! Out: h(:,:,1) = Left eigenvectors
                                                 !      h(:,:,2) = Right eigenvectors
    complex(kind=drk), intent(out)   :: e(:)     ! Out: Eigenvalues

    complex(kind=drk) :: work(50*size(h,dim=2))
    real(drk)         :: rwork(3*size(h,dim=2))
    complex(kind=drk) :: vl(size(h,dim=2),size(h,dim=2))
    integer           :: info     ! Must be of default integer kind
    integer           :: nh1, nh2 ! Must be of default integer kind
    
    nh1 = size(h,dim=1) ; nh2 = size(h,dim=2)
    if (nh1<nh2) stop 'lapack_zgeev2 - oops'
    call zgeev('V','V',nh2,h(1:nh1,1:nh2,1),nh1,e(:),vl,nh2,h(1:nh1,1:nh2,2),nh1,work,size(work),rwork,info)
    if (info/=0) then
      write (out,"(' zgeev returned ',i8)") info
      stop 'lapack_zgeev2 - zgeev failed'
    end if
    h(1:nh2,1:nh2,1) = vl
  end subroutine lapack_zgeev2

!*lq subroutine lapack_quad_zgeev2(h,e)
!*lq   complex(qrk), intent(inout) :: h(:,:,:) ! In:  h(:,:,1) = hermitian matrix to be diagonalized
!*lq                                           ! Out: h(:,:,1) = Left eigenvectors
!*lq                                           !      h(:,:,2) = Right eigenvectors
!*lq   complex(qrk), intent(out)   :: e(:)     ! Out: Eigenvalues
   
!*lq   complex(qrk) :: work(50*size(h,dim=2))
!*lq   real(qrk)    :: rwork(3*size(h,dim=2))
!*lq   complex(qrk) :: vl(size(h,dim=2),size(h,dim=2))
!*lq   integer      :: info      ! Must be of default integer kind
!*lq   integer      :: nh1, nh2  ! Must be of default integer kind
!*lq   
!*lq   nh1 = size(h,dim=1) ; nh2 = size(h,dim=2)
!*lq   if (nh1<nh2) stop 'lapack%lapack_quad_zgeev2 - oops'
!*lq   call quad_zgeev('V','V',nh2,h(1:nh1,1:nh2,1),nh1,e(:),vl,nh2,h(1:nh1,1:nh2,2),nh1,work,size(work),rwork,info)
!*lq   if (info/=0) then
!*lq     write (out,"(' quad_zgeev returned ',i8)") info
!*lq     stop 'lapack%lapack_quad_zgeev2 - quad_zgeev failed'
!*lq   end if
!*lq   h(1:nh2,1:nh2,1) = vl
!*lq end subroutine lapack_quad_zgeev2

  subroutine lapack_cheev(h,e)
    complex(srk), intent(inout) :: h(:,:)  ! In:  Hermitian matrix to be diagonalized
                                           ! Out: Eigenvectors
    real(srk), intent(out)   :: e(:)       ! Out: Eigenvalues

    complex(srk) :: work(50*size(h,dim=1))
    real(srk)    :: rwork(3*size(h,dim=1))
    integer      :: info                   ! Must be of default integer kind
    integer      :: nh1, nh2               ! Must be of default integer kind
    
    nh1 = size(h,dim=1) ; nh2 = size(h,dim=2)
    call cheev('V','U',nh1,h(1:nh1,1:nh2),nh1,e(:),work,size(work),rwork,info)
    if (info/=0) then
      write (out,"(' cheev returned ',i8)") info
      stop 'lapack_cheev - cheev failed'
    end if
  end subroutine lapack_cheev

  subroutine lapack_zheev(h,e)
    complex(drk), intent(inout) :: h(:,:)  ! In:  Hermitian matrix to be diagonalized
                                           ! Out: Eigenvectors
    real(drk), intent(out)   :: e(:)       ! Out: Eigenvalues

    complex(drk) :: work(50*size(h,dim=1))
    real(drk)    :: rwork(3*size(h,dim=1))
    integer      :: info            ! Must be of default integer kind
    integer      :: nh1, nh2        ! Must be of default integer kind
    
    nh1 = size(h,dim=1) ; nh2 = size(h,dim=2)
    call zheev('V','U',nh1,h(1:nh1,1:nh2),nh1,e(:),work,size(work),rwork,info)
    if (info/=0) then
      write (out,"(' zheev returned ',i8)") info
      stop 'lapack_zheev - zheev failed'
    end if
  end subroutine lapack_zheev

  subroutine lapack_dsyev(h,e)
    real(drk), intent(inout) :: h(:,:)  ! In:  symmetric matrix to be diagonalized
                                               ! Out: Eigenvectors
    real(drk), intent(out)   :: e(:)    ! Out: Eigenvalues

    real(drk) :: work(50*size(h,dim=1))
    integer   :: info      ! Must be of default integer kind
    integer   :: nh1, nh2  ! Must be of default integer kind
    
    nh1 = size(h,dim=1) ; nh2 = size(h,dim=2)
    call dsyev('V','U',nh1,h(1:nh1,1:nh2),nh1,e(:),work,size(work),info)
    if (info/=0) then
      write (out,"(' dsyev returned ',i8)") info
      stop 'lapack_dsyev - dsyev failed'
    end if
  end subroutine lapack_dsyev

!*lq subroutine lapack_quad_dsyev(h,e)
!*lq   real(qrk), intent(inout) :: h(:,:)  ! In:  symmetric matrix to be diagonalized
!*lq                                       ! Out: Eigenvectors
!*lq   real(qrk), intent(out)   :: e(:)    ! Out: Eigenvalues
   
!*lq   real(qrk) :: work(50*size(h,dim=1))
!*lq   integer   :: info       ! Must be of default integer kind
!*lq   integer   :: nh1, nh2   ! Must be of default integer kind
!*lq   
!*lq   nh1 = size(h,dim=1) ; nh2 = size(h,dim=2)
!*lq   call quad_dsyev('V','U',nh1,h(1:nh1,1:nh2),nh1,e(:),work,size(work),info)
!*lq   if (info/=0) then
!*lq     write (out,"(' quad_dsyev returned ',i8)") info
!*lq     stop 'lapack_quad_dsyev - dsyev failed'
!*lq   end if
!*lq end subroutine lapack_quad_dsyev

  subroutine lapack_ssyev(h,e)
    real(srk), intent(inout) :: h(:,:)  ! In:  symmetric matrix to be diagonalized
                                        ! Out: Eigenvectors
    real(srk), intent(out)   :: e(:)    ! Out: Eigenvalues

    real(srk)        :: work(50*size(h,dim=1))
    integer          :: info       ! Must be of default integer kind
    integer          :: nh1, nh2   ! Must be of default integer kind
    
    nh1 = size(h,dim=1) ; nh2 = size(h,dim=2)
    call ssyev('V','U',nh1,h(1:nh1,1:nh2),nh1,e(:),work,size(work),info)
    if (info/=0) then
      write (out,"(' ssyev returned ',i8)") info
      stop 'lapack_ssyev - ssyev failed'
    end if
  end subroutine lapack_ssyev

  subroutine lapack_dgesvd(a,s,u,vth)
    real(drk), intent(inout) :: a  (:,:)  ! In:  Matrix to be decomposed
                                          ! Out: Content destroyed
    real(drk), intent(out)   :: s  (:)    ! Out: Singular values
    real(drk), intent(out)   :: u  (:,:)  ! Out: Left singular vectors
    real(drk), intent(out)   :: vth(:,:)  ! Out: Right singular vectors, transposed & conjugated
                                          ! The overall result is A = U S VTH

    character(len=1)       :: jobu        ! Either 'A' (all) or 'S' (singular), depending on the
    character(len=1)       :: jobvth      ! sizes of u and vth arrays
    real(drk)              :: lwq(1)
    real(drk), allocatable :: work(:)
    integer                :: info, lwork                    ! Must be of default integer kind
    integer                :: m, n, lda, ldu, ldvth, nsing   ! Must be of default integer kind
   
    m     = size(a,dim=1) 
    n     = size(a,dim=2) 
    lda   = m
    nsing = min(m,n)
    ldu   = size(u,dim=1)
    ldvth = size(vth,dim=1)
    !
    if (size(s)<nsing) stop 'lapack%lapack_dgesvd - array s is too small'
    !
    if (size(u,dim=1)<m    ) stop 'lapack%lapack_dgesvd - array u is too small (1)'
    if (size(u,dim=2)<nsing) stop 'lapack%lapack_dgesvd - array u is too small (2)'
    jobu = 'S'
    if (size(u,dim=2)>=m) jobu = 'A'
    !
    if (size(vth,dim=2)<n    ) stop 'lapack%lapack_dgesvd - array vth is too small (2)'
    if (size(vth,dim=1)<nsing) stop 'lapack%lapack_dgesvd - array vth is too small (1)'
    jobvth = 'S'
    if (size(vth,dim=1)>=n) jobvth = 'A'
    !
    call dgesvd(jobu,jobvth,m,n,a,lda,s,u,ldu,vth,ldvth,lwq, -1,   info)
    if (info/=0) then
      write (out,"(' dgesvd returned ',i8,' for workspace query')") info
      stop 'lapack_dgesvd - dgesvd failed'
    end if
    !
    lwork = 1+nint(lwq(1))
    allocate (work(lwork),stat=info)
    if (info/=0) then
      write (out,"(' Error ',i8,' allocating ',i10,'-element array for dgesvd')") info, lwork
      stop 'lapack_dgesvd - dgesvd failed'
    end if
    !
    call dgesvd(jobu,jobvth,m,n,a,lda,s,u,ldu,vth,ldvth,work,lwork,info)
    if (info/=0) then
      write (out,"(' dgesvd returned ',i8)") info
      stop 'lapack_dgesvd - dgesvd failed'
    end if
    deallocate (work)
  end subroutine lapack_dgesvd

  subroutine lapack_zgesvd(a,s,u,vth)
    complex(kind=drk), intent(inout) :: a  (:,:)  ! In:  Matrix to be decomposed
                                                  ! Out: Content destroyed
    real(drk), intent(out)    :: s  (:)    ! Out: Singular values
    complex(kind=drk), intent(out)   :: u  (:,:)  ! Out: Left singular vectors
    complex(kind=drk), intent(out)   :: vth(:,:)  ! Out: Right singular vectors, transposed & conjugated
                                                  ! The overall result is A = U S VTH

    complex(drk) :: work(50*max(size(a,dim=1),size(a,dim=2)))
    real(drk)    :: rwork(5*size(a,dim=1))
    integer      :: info, lwork           ! Must be of default integer kind
    integer      :: m, n, lda, ldu, ldvth ! Must be of default integer kind
   
    m     = size(a,dim=1) ; lda = m
    n     = size(a,dim=2) 
    ldu   = size(u,dim=1)
    ldvth = size(vth,dim=1)
    lwork = size(work)
    if (size(s)<min(m,n))               stop 'lapack%lapack_zgesvd - array s is too small'
    if (ldu<m .or. size(u,dim=2)<m)     stop 'lapack%lapack_zgesvd - array u is too small'
    if (ldvth<n .or. size(vth,dim=2)<n) stop 'lapack%lapack_zgesvd - array vth is too small'
    !
    call zgesvd('A','A',m,n,a,lda,s,u,ldu,vth,ldvth,work,lwork,rwork,info)
    if (info/=0) then
      write (out,"(' zgesvd returned ',i8)") info
      stop 'lapack_zgesvd - zgesvd failed'
    end if
  end subroutine lapack_zgesvd

  subroutine lapack_dsterf(a,b)
    real(drk), intent(inout) :: a(:) ! In: Diagonal elements of the tri-diagonal matrix
                                            ! Out: Eigenvalues, in the ascending order
    real(drk), intent(inout) :: b(:) ! In: Sub-diagonal elements of the tri-diagonal matirx
                                            ! Out: Destroyed
    !
    integer :: na, nb ! Must be of default integer kind
    integer :: info   ! Must be of default integer kind
    !
    na = size(a)
    nb = size(b)
    if (na/=nb+1) then
      write (out,"('lapack_dsterf: inconsistent array sizes: diagonal ',i6,' subdiagonal ',i6)") na, nb
      stop 'lapack_dsterf - bad input'
    end if
    call dsterf(na,a,b,info)
    if (info/=0) then
      write (out,"(' dsterf returned ',i8)") info
      stop 'lapack_dsterf - dsterf failed'
    end if
  end subroutine lapack_dsterf

  subroutine lapack_ssterf(a,b)
    real(srk), intent(inout) :: a(:) ! In: Diagonal elements of the tri-diagonal matrix
                                     ! Out: Eigenvalues, in the ascending order
    real(srk), intent(inout) :: b(:) ! In: Sub-diagonal elements of the tri-diagonal matirx
                                     ! Out: Destroyed
    !
    integer :: na, nb  ! Must be of default integer kind
    integer :: info    ! Must be of default integer kind
    !
    na = size(a)
    nb = size(b)
    if (na/=nb+1) then
      write (out,"('lapack_ssterf: inconsistent array sizes: diagonal ',i6,' subdiagonal ',i6)") na, nb
      stop 'lapack_ssterf - bad input'
    end if
    call ssterf(na,a,b,info)
    if (info/=0) then
      write (out,"(' ssterf returned ',i8)") info
      stop 'lapack_ssterf - ssterf failed'
    end if
  end subroutine lapack_ssterf

  subroutine lapack_ginverse_real(amat,power_,eps_)
    real(srk), intent(inout)        :: amat(:,:) ! In: matrix to invert
                                                 ! Out: generalized inverse of the matrix
    real(srk), intent(in), optional :: power_    ! In: power of the inverse; 1/2 if omitted
                                                 !     (corresponding to A^(-1))
    real(srk), intent(in), optional :: eps_      ! In: Eigenvalue tolerance; negative values
                                                 !     will use the default.
    !
    real(srk) :: eval(size(amat,dim=1))
    real(srk) :: evec(size(amat,dim=1),size(amat,dim=1))
    real(srk) :: eps
    real(srk) :: power
    !
    power = 0.5
    if (present(power_)) power = power_
    !
    evec  = amat
    call lapack_syev(evec(:,:),eval(:))
    eps = -1
    if (present(eps_)) eps = eps_
    if (eps<0) eps = 100.0*spacing(maxval(eval))
    where (abs(eval)>eps) 
      eval = 1.0 / eval**power
    elsewhere  
      eval = 0.0
    end where
    evec = evec * spread(eval,dim=1,ncopies=size(evec,dim=1))
    amat = matmul(evec,transpose(evec))
    !
  end subroutine lapack_ginverse_real

  subroutine lapack_ginverse_complex(amat,power_,eps_)
    complex(srk), intent(inout)     :: amat(:,:) ! In: matrix to invert
                                                 ! Out: generalized inverse of the matrix
    real(srk), intent(in), optional :: power_    ! In: power of the inverse; 1/2 if omitted
                                                 !     (corresponding to A^(-1))
    real(srk), intent(in), optional :: eps_      ! In: Eigenvalue tolerance; negative values
                                                 !     will use the default.
    !
    real(srk)    :: eval(size(amat,dim=1))
    complex(srk) :: evec(size(amat,dim=1),size(amat,dim=1))
    real(srk)    :: eps
    real(srk)    :: power
    !
    power = 0.5
    if (present(power_)) power = power_
    !
    evec  = amat
    call lapack_heev(evec(:,:),eval(:))
    eps = -1
    if (present(eps_)) eps = eps_
    if (eps<0) eps = 100.0*spacing(maxval(eval))
    where (abs(eval)>eps) 
      eval = 1.0 / eval**power
    elsewhere  
      eval = 0.0
    end where
    evec = evec * spread(eval,dim=1,ncopies=size(evec,dim=1))
    amat = matmul(evec,transpose(conjg(evec)))
    !
  end subroutine lapack_ginverse_complex

  subroutine lapack_ginverse_double(amat,power_,eps_)
    real(drk), intent(inout) :: amat(:,:)     ! In: matrix to invert
                                                     ! Out: generalized inverse of the matrix
    real(drk), intent(in), optional :: power_ ! In: power of the inverse; 1/2 if omitted
                                                     !     (corresponding to A^(-1))
    real(drk), intent(in), optional :: eps_      ! In: Eigenvalue tolerance; negative values
                                                 !     will use the default.

    real(drk) :: eval(size(amat,dim=1))
    real(drk) :: evec(size(amat,dim=1),size(amat,dim=1))
    real(drk) :: eps
    real(drk) :: power
    !
    power = 0.5d0
    if (present(power_)) power = power_
    !
    evec = amat
    call lapack_syev(evec(:,:),eval(:))
    eps = -1
    if (present(eps_)) eps = eps_
    if (eps<0) eps = 100.0d0*spacing(maxval(eval))
    where (abs(eval)>eps) 
      eval = 1.0d0 / eval**power
    elsewhere  
      eval = 0.0d0
    end where
    evec = evec * spread(eval,dim=1,ncopies=size(evec,dim=1))
    amat = matmul(evec,transpose(evec))
    !
  end subroutine lapack_ginverse_double

!*lq  subroutine lapack_ginverse_quad(amat,power_,eps_)
!*lq    real(qrk), intent(inout) :: amat(:,:)     ! In: matrix to invert
!*lq                                              ! Out: generalized inverse of the matrix
!*lq    real(qrk), intent(in), optional :: power_ ! In: power of the inverse; 1/2 if omitted
!*lq                                              !     (corresponding to A^(-1))
!*lq    real(qrk), intent(in), optional :: eps_   ! In: Eigenvalue tolerance; negative values
!*lq                                              !     will use the default.
!*lq
!*lq    real(qrk) :: eval(size(amat,dim=1))
!*lq    real(qrk) :: evec(size(amat,dim=1),size(amat,dim=1))
!*lq    real(qrk) :: eps
!*lq    real(qrk) :: power
!*lq    !
!*lq    power = 0.5_qrk
!*lq    if (present(power_)) power = power_
!*lq    !
!*lq    evec = amat
!*lq    call lapack_syev(evec(:,:),eval(:))
!*lq    eps = -1
!*lq    if (present(eps_)) eps = eps_
!*lq    if (eps<0) eps = 100.0_qrk*spacing(maxval(eval))
!*lq    where (abs(eval)>eps) 
!*lq      eval = 1.0_qrk / eval**power
!*lq    elsewhere  
!*lq      eval = 0.0_qrk
!*lq    end where
!*lq    evec = evec * spread(eval,dim=1,ncopies=size(evec,dim=1))
!*lq    amat = matmul(evec,transpose(evec))
!*lq    !
!*lq  end subroutine lapack_ginverse_quad

  subroutine lapack_ginverse_doublecomplex(amat,power_)
    complex(kind=drk), intent(inout)     :: amat(:,:)   ! In: matrix to invert
                                                     ! Out: generalized inverse of the matrix
    real(drk), intent(in), optional :: power_ ! In: power of the inverse; 1/2 if omitted
                                                     !     (corresponding to A^(-1))
    !
    real(drk) :: eval(size(amat,dim=1))
    complex(kind=drk)   :: evec(size(amat,dim=1),size(amat,dim=1))
    real(drk) :: eps
    real(drk) :: power
    !
    power = 0.5d0
    if (present(power_)) power = power_
    !
    evec  = amat
    call lapack_heev(evec(:,:),eval(:))
    eps = 100.0d0*spacing(maxval(eval))
    where (abs(eval)>eps) 
      eval = 1.0d0 / eval**power
    elsewhere  
      eval = 0.0
    end where
    evec = evec * spread(eval,dim=1,ncopies=size(evec,dim=1))
    amat = matmul(evec,transpose(conjg(evec)))
    !
  end subroutine lapack_ginverse_doublecomplex

  subroutine lapack_dgesv(a,b)
    real(drk), intent(inout) :: a(:,:) ! In: Linear system matrix
                                       ! Out: Destroyed
    real(drk), intent(inout) :: b(:,:) ! In: Right-half side
                                       ! Out: Solutions
    !
    integer :: info, n, nrhs, lda, ldb
    integer :: ipiv(size(a,dim=2))
    !
    n    = size(a,dim=2)
    nrhs = size(b,dim=2)
    lda  = size(a,dim=1)
    ldb  = size(b,dim=1)
    if (lda<n .or. ldb<n) stop 'lapack%lapack_dgesv - crazy inputs'
    call dgesv(n,nrhs,a,lda,ipiv,b,ldb,info)
    if (info/=0) then
      write (out,"('lapack%lapack_dgesv failed with code ',i0)") info
      stop 'lapack%lapack_dgesv failed'
    end if
  end subroutine lapack_dgesv

  subroutine lapack_cgesv(a,b)
    complex(srk), intent(inout) :: a(:,:) ! In: Linear system matrix
                                          ! Out: Destroyed
    complex(srk), intent(inout) :: b(:,:) ! In: Right-half side
                                          ! Out: Solutions
    !
    integer :: info, n, nrhs, lda, ldb
    integer :: ipiv(size(a,dim=2))
    !
    n    = size(a,dim=2)
    nrhs = size(b,dim=2)
    lda  = size(a,dim=1)
    ldb  = size(b,dim=1)
    if (lda<n .or. ldb<n) stop 'lapack%lapack_cgesv - crazy inputs'
    call cgesv(n,nrhs,a,lda,ipiv,b,ldb,info)
    if (info/=0) then
      write (out,"('lapack%lapack_cgesv failed with code ',i0)") info
      stop 'lapack%lapack_cgesv failed'
    end if
  end subroutine lapack_cgesv

  subroutine lapack_zgesv(a,b)
    complex(drk), intent(inout) :: a(:,:) ! In: Linear system matrix
                                          ! Out: Destroyed
    complex(drk), intent(inout) :: b(:,:) ! In: Right-half side
                                          ! Out: Solutions
    !
    integer :: info, n, nrhs, lda, ldb
    integer :: ipiv(size(a,dim=2))
    !
    n    = size(a,dim=2)
    nrhs = size(b,dim=2)
    lda  = size(a,dim=1)
    ldb  = size(b,dim=1)
    if (lda<n .or. ldb<n) stop 'lapack%lapack_zgesv - crazy inputs'
    call zgesv(n,nrhs,a,lda,ipiv,b,ldb,info)
    if (info/=0) then
      write (out,"('lapack%lapack_zgesv failed with code ',i0)") info
      stop 'lapack%lapack_zgesv failed'
    end if
  end subroutine lapack_zgesv

!*lq subroutine lapack_quad_dgesv(a,b)
!*lq   real(qrk), intent(inout) :: a(:,:) ! In: Linear system matrix
!*lq                                      ! Out: Destroyed
!*lq   real(qrk), intent(inout) :: b(:,:) ! In: Right-half side
!*lq                                      ! Out: Solutions
!*lq   !
!*lq   integer :: info, n, nrhs, lda, ldb
!*lq   integer :: ipiv(size(a,dim=2))
!*lq   !
!*lq   n    = size(a,dim=2)
!*lq   nrhs = size(b,dim=2)
!*lq   lda  = size(a,dim=1)
!*lq   ldb  = size(b,dim=1)
!*lq   if (lda<n .or. ldb<n) stop 'lapack%lapack_quad_dgesv - crazy inputs'
!*lq   call quad_dgesv(n,nrhs,a,lda,ipiv,b,ldb,info)
!*lq   if (info/=0) then
!*lq     write (out,"('lapack%lapack_quad_dgesv failed with code ',i0)") info
!*lq     stop 'lapack%lapack_quad_dgesv failed'
!*lq   end if
!*lq end subroutine lapack_quad_dgesv

!*lq subroutine lapack_quad_zgesv(a,b)
!*lq   complex(qrk), intent(inout) :: a(:,:) ! In: Linear system matrix
!*lq                                         ! Out: Destroyed
!*lq   complex(qrk), intent(inout) :: b(:,:) ! In: Right-half side
!*lq                                         ! Out: Solutions
!*lq   !
!*lq   integer :: info, n, nrhs, lda, ldb
!*lq   integer :: ipiv(size(a,dim=2))
!*lq   !
!*lq   n    = size(a,dim=2)
!*lq   nrhs = size(b,dim=2)
!*lq   lda  = size(a,dim=1)
!*lq   ldb  = size(b,dim=1)
!*lq   if (lda<n .or. ldb<n) stop 'lapack%lapack_quad_zgesv - crazy inputs'
!*lq   call quad_zgesv(n,nrhs,a,lda,ipiv,b,ldb,info)
!*lq   if (info/=0) then
!*lq     write (out,"('lapack%lapack_quad_zgesv failed with code ',i0)") info
!*lq     stop 'lapack%lapack_quad_zgesv failed'
!*lq   end if
!*lq end subroutine lapack_quad_zgesv

  function lapack_determinant_double(mat) result(det)
    real(drk), intent(in)  :: mat(:,:) ! Matrix to compute determinant for
    real(drk)              :: det
    !
    real(drk), allocatable :: tm(:,:)
    integer                :: order, info ! Must be of default integer kind
    integer, allocatable   :: ipiv(:)     ! Must be of default integer kind
    integer                :: i
    !
    order = size(mat,dim=1)
    if (order<=0) stop 'lapack%lapack_determinant_double - zero-size matrix'
    if (size(mat,dim=2)/=order) stop 'lapack%lapack_determinant_double - matrix not square'
    !
    allocate (tm(order,order),ipiv(order),stat=info)
    if (info/=0) stop 'lapack%lapack_determinant_double - no memory'
    !
    tm = mat
    call dgetrf(order,order,tm,order,ipiv,info)
    if (info<0) stop 'lapack%lapack_determinant_double - dgetrf failed'
    !
    det = 1
    compute_determinant: do i=1,order
      det = det * tm(i,i)
      if (ipiv(i)/=i) det = -det
    end do compute_determinant
    !
    deallocate(tm,ipiv,stat=info)
    if (info/=0) stop 'lapack%lapack_determinant_double - memory deallocation failed'
    !
    ! write (out,"('DEBUG: double lapack det = ',g30.20,' diff from linpack = ',g30.20)") &
    !        det, det - linpack_determinant_double(mat)
  end function lapack_determinant_double
  !
!*lq function lapack_determinant_quad(mat) result(det)
!*lq   real(qrk), intent(in)  :: mat(:,:)    ! Matrix to compute determinant for
!*lq   real(qrk)              :: det
!*lq   !
!*lq   real(qrk), allocatable :: tm(:,:)
!*lq   integer                :: order, info ! Must be of default integer kind
!*lq   integer, allocatable   :: ipiv(:)     ! Must be of default integer kind
!*lq   integer                :: i
!*lq   !
!*lq   order = size(mat,dim=1)
!*lq   if (order<=0) stop 'lapack%lapack_determinant_quad - zero-size matrix'
!*lq   if (size(mat,dim=2)/=order) stop 'lapack%lapack_determinant_quad - matrix not square'
!*lq   !
!*lq   allocate (tm(order,order),ipiv(order),stat=info)
!*lq   if (info/=0) stop 'lapack%lapack_determinant_quad - no memory'
!*lq   !
!*lq   tm = mat
!*lq   call quad_dgetrf(order,order,tm,order,ipiv,info)
!*lq   if (info<0) stop 'lapack%lapack_determinant_quad - dgetrf failed'
!*lq   !
!*lq   det = 1
!*lq   compute_determinant: do i=1,order
!*lq     det = det * tm(i,i)
!*lq     if (ipiv(i)/=i) det = -det
!*lq   end do compute_determinant
!*lq   !
!*lq   deallocate(tm,ipiv,stat=info)
!*lq   if (info/=0) stop 'lapack%lapack_determinant_quad - memory deallocation failed'
!*lq   !
!*lq   ! write (out,"('DEBUG: quad lapack det = ',g30.20,' diff from linpack = ',g30.20)") &
!*lq   !        det, det - linpack_determinant_double(real(mat,kind=drk))
!*lq end function lapack_determinant_quad
  !
  function linpack_determinant_double(mat) result(det)
    real(drk), intent(in)  :: mat(:,:) ! Matrix to compute determinant for
    real(drk)              :: det
    !
    real(drk), allocatable :: tm(:,:)
    integer                :: order, info            ! Must be of default integer kind
    integer                :: ipvt(size(mat,dim=1))  ! Must be of default integer kind
    real(drk)              :: work(size(mat,dim=1))
    real(drk)              :: detx(2)
    external               :: dgefa, dgedi
    !
    order = size(mat,dim=1)
    if (size(mat,dim=2)/=order) then
      write (out,"('Determinant requested for a non-square matrix: ',2i5,'. Bummer.')") &
             order, size(mat,dim=2)
      stop 'lapack%linpack_determinant_double - bad input'
    end if
    !
    allocate (tm(order,order),stat=info)
    if (info/=0) then
      write (out,"('Error ',i5,' allocating order-',i5,' matrix.')") info, order
      stop 'lapack%linpack_determinant_double - no memory'
    end if
    tm = mat
    !
    call dgefa(tm,order,order,ipvt,info)
    !
    call dgedi(tm,order,order,ipvt,detx,work,10)
    !
    ! tm = mat
    ! call lapack_dsyev(tm,work)
    ! write (out,"('Diagonalization gives: ',g40.20/' linpack gives ',g40.20/' Diff = ',g40.20)") &
    !        product(work), detx(1) * 10.0d0**detx(2), product(work) - detx(1) * 10.0d0**detx(2)
    !
    deallocate(tm,stat=info)
    if (info/=0) then
      write (out,"('Error ',i5,' deallocating order-',i5,' matrix.')") info, order
      stop 'lapack%linpack_determinant_double - memory deallocation failed'
    end if
    !
    det = detx(1) * 10.0d0**detx(2)
  end function linpack_determinant_double

  function linpack_determinant_double_trash_input(mat) result(det)
    real(drk), intent(inout)  :: mat(:,:) ! Matrix to compute determinant for
    real(drk)                 :: det
    !
    integer                   :: order, info            ! Must be of default integer kind
    integer                   :: ipvt(size(mat,dim=1))  ! Must be of default integer kind
    real(drk)                 :: work(size(mat,dim=1))
    real(drk)                 :: detx(2)
    external                  :: dgefa, dgedi
    !
    order = size(mat,dim=1)
    if (size(mat,dim=2)/=order) then
      write (out,"('Determinant requested for a non-square matrix: ',2i5,'. Bummer.')") &
             order, size(mat,dim=2)
      stop 'lapack%linpack_determinant_double_trash - bad input'
    end if
    !
    call dgefa(mat,order,order,ipvt,info)
    if (info>0) then
      !  Zero pivot; determinant vanishes
      det = 0
      return
    end if
    !
    call dgedi(mat,order,order,ipvt,detx,work,10)
    !
    det = detx(1) * 10.0d0**detx(2)
  end function linpack_determinant_double_trash_input

end module lapack
