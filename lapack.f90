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
  !  Simplistic type-agnostic LAPACK interface (subset)
  !
  use accuracy
  implicit none

  character(len=clen), save :: rcsid_lapack = "$Id: lapack.f90,v 1.13 2025/07/11 15:08:35 ps Exp $"

  interface lapack_geev
     module procedure lapack_cgeev2
     module procedure lapack_zgeev2
!*lq module procedure lapack_quad_zgeev2
  end interface lapack_geev

  contains

  subroutine lapack_cgeev2(nh,h,e)
    integer(ik), intent(in)     :: nh         ! Dimensions of the arrays; needed to avoid temporaries
    complex(srk), intent(inout) :: h(nh,nh,2) ! In:  h(:,:,1) = General matrix to be diagonalized
                                              ! Out: h(:,:,1) = Left eigenvectors
                                              !      h(:,:,2) = Right eigenvectors
    complex(srk), intent(out)   :: e(nh)      ! Out: Eigenvalues

    complex(srk), allocatable :: work(:), vl(:,:)
    real(srk), allocatable    :: rwork(:)
    !
    integer(lik) :: lnh      ! Lapack-type version of nh
    integer(lik) :: info     ! Error code; must be of default integer kind
    integer(lik) :: lwork    ! Optimal work size; must be of default integer kind
    external     :: cgeev
   
    lnh = nh
    if (int(lnh,kind=ik)/=nh) then
      write (out,"('lapack_cgeev2: Bad integer conversion: nh = ',i0,' -> ',i0,' -> ',i0)") nh, lnh, int(lnh,kind=ik)
      stop 'lapack%lapack_cgeev2 - Bad conversion to LAPACK integer type'
    end if
    allocate (vl(lnh,lnh),work(1),rwork(2*lnh),stat=info)
    if (info/=0) then
      write (out,"('lapack_cgeev2: allocate failed (1) with code ',i0)") info
      stop 'lapack_cgeev2 - allocate failed (1)'
    end if

    call cgeev('V','V',lnh,h(:,:,1),lnh,e,vl,lnh,h(:,:,2),lnh,work,-1,rwork,info) ! Expect gfortran warning
    if (info/=0) then
      write (out,"(' cgeev (1) returned ',i0)") info
      stop 'lapack_cgeev2 - cgeev failed (1)'
    end if

    lwork = int(work(1))
    deallocate (work)
    allocate (work(lwork),stat=info)
    if (info/=0) then
      write (out,"('lapack_cgeev2: allocate failed (2) with code ',i0)") info
      stop 'lapack_cgeev2 - allocate failed (2)'
    end if

    call cgeev('V','V',lnh,h(:,:,1),lnh,e,vl,lnh,h(:,:,2),lnh,work,lwork,rwork,info) ! Expect gfortran warning
    if (info/=0) then
      write (out,"(' cgeev (2) returned ',i0)") info
      stop 'lapack_cgeev2 - cgeev failed (2)'
    end if

    h(:,:,1) = vl

    deallocate (vl,work,rwork)
  end subroutine lapack_cgeev2

  subroutine lapack_zgeev2(nh,h,e)
    integer(ik), intent(in)     :: nh         ! Dimensions of the arrays; needed to avoid temporaries
    complex(drk), intent(inout) :: h(nh,nh,2) ! In:  h(:,:,1) = General matrix to be diagonalized
                                              ! Out: h(:,:,1) = Left eigenvectors
                                              !      h(:,:,2) = Right eigenvectors
    complex(drk), intent(out)   :: e(nh)      ! Out: Eigenvalues

    complex(drk), allocatable :: work(:), vl(:,:)
    real(drk), allocatable    :: rwork(:)
    !
    integer(lik) :: lnh      ! Lapack-type version of nh
    integer(lik) :: info     ! Error code; must be of default integer kind
    integer(lik) :: lwork    ! Optimal work size; must be of default integer kind
    external     :: zgeev
    
    lnh = nh
    if (int(lnh,kind=ik)/=nh) then
      write (out,"('lapack_zgeev2: Bad integer conversion: nh = ',i0,' -> ',i0,' -> ',i0)") nh, lnh, int(lnh,kind=ik)
      stop 'lapack%lapack_cgeev2 - Bad conversion to LAPACK integer type'
    end if
    allocate (vl(lnh,lnh),work(1),rwork(2*lnh),stat=info)
    if (info/=0) then
      write (out,"('lapack_zgeev2: allocate failed (1) with code ',i0)") info
      stop 'lapack_zgeev2 - allocate failed (1)'
    end if

    call zgeev('V','V',lnh,h(:,:,1),lnh,e,vl,lnh,h(:,:,2),lnh,work,-1,rwork,info) ! Expect gfortran warning
    if (info/=0) then
      write (out,"(' zgeev (1) returned ',i0)") info
      stop 'lapack_zgeev2 - zgeev failed (1)'
    end if

    lwork = int(work(1))
    deallocate (work)
    allocate (work(lwork),stat=info)
    if (info/=0) then
      write (out,"('lapack_zgeev2: allocate failed (2) with code ',i0)") info
      stop 'lapack_zgeev2 - allocate failed (2)'
    end if

    call zgeev('V','V',lnh,h(:,:,1),lnh,e,vl,lnh,h(:,:,2),lnh,work,lwork,rwork,info) ! Expect gfortran warning
    if (info/=0) then
      write (out,"(' zgeev (2) returned ',i0)") info
      stop 'lapack_zgeev2 - zgeev failed (2)'
    end if

    h(:,:,1) = vl

    deallocate (vl,work,rwork)
  end subroutine lapack_zgeev2

!*lq  subroutine lapack_quad_zgeev2(nh,h,e)
!*lq    integer(ik), intent(in)     :: nh         ! Dimensions of the arrays; needed to avoid temporaries
!*lq    complex(qrk), intent(inout) :: h(nh,nh,2) ! In:  h(:,:,1) = General matrix to be diagonalized
!*lq                                              ! Out: h(:,:,1) = Left eigenvectors
!*lq                                              !      h(:,:,2) = Right eigenvectors
!*lq    complex(qrk), intent(out)   :: e(nh)      ! Out: Eigenvalues
!*lq
!*lq    complex(qrk), allocatable :: work(:), vl(:,:)
!*lq    real(qrk), allocatable    :: rwork(:)
!*lq    !
!*lq    integer(lik) :: lnh      ! Lapack-type version of nh
!*lq    integer      :: info     ! Error code; must be of default integer kind
!*lq    integer      :: lwork    ! Optimal work size; must be of default integer kind
!*lq    external     :: lapack_quad_zgeev2
!*lq    
!*lq    lnh = nh
!*lq    if (int(lnh,kind=ik)/=nh) then
!*lq      write (out,"('lapack_quad_cgeev2: Bad integer conversion: nh = ',i0,' -> ',i0,' -> ',i0)") nh, lnh, int(lnh,kind=ik)
!*lq      stop 'lapack%lapack_quad_cgeev2 - Bad conversion to LAPACK integer type'
!*lq    end if
!*lq    allocate (vl(lnh,lnh),work(1),rwork(2*lnh),stat=info)
!*lq    if (info/=0) then
!*lq      write (out,"('lapack_quad_zgeev2: allocate failed (1) with code ',i0)") info
!*lq      stop 'lapack_quad_zgeev2 - allocate failed (1)'
!*lq    end if
!*lq
!*lq    call quad_zgeev('V','V',lnh,h(:,:,1),lnh,e,vl,lnh,h(:,:,2),lnh,work,-1,rwork,info) ! Expect gfortran warning
!*lq    if (info/=0) then
!*lq      write (out,"(' quad_zgeev (1) returned ',i0)") info
!*lq      stop 'lapack_quad_zgeev2 - quad_zgeev failed (1)'
!*lq    end if
!*lq
!*lq    lwork = int(work(1))
!*lq    deallocate (work)
!*lq    allocate (work(lwork),stat=info)
!*lq    if (info/=0) then
!*lq      write (out,"('lapack_quad_zgeev2: allocate failed (2) with code ',i0)") info
!*lq      stop 'lapack_quad_zgeev2 - allocate failed (2)'
!*lq    end if
!*lq
!*lq    call quad_zgeev('V','V',lnh,h(:,:,1),lnh,e,vl,lnh,h(:,:,2),lnh,work,lwork,rwork,info) ! Expect gfortran warning
!*lq    if (info/=0) then
!*lq      write (out,"(' quad_zgeev (2) returned ',i0)") info
!*lq      stop 'lapack_quad_zgeev2 - quad_zgeev failed (2)'
!*lq    end if
!*lq
!*lq    h(:,:,1) = vl
!*lq
!*lq    deallocate (vl,work,rwork)
!*lq  end subroutine lapack_quad_zgeev2

end module lapack
