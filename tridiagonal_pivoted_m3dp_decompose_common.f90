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
! !
! !  m3dp_decompose() is a clone of LAPACK _GTTRF(). It computes LU factorization of
! !  a tridiagonal matrix using elimination with partial pivoting and row interchanges.
! !
! subroutine m3dp_decompose_r(m,mf,lpiv,fail)
!   real(rk), intent(in)           :: m (:,:) ! Tridiagonal matrix to factorize;
!                                             ! m(i,1) is main diagonal (i,i)
!                                             ! m(i,2) is sub-diagonal (i+1,i)
!                                             ! m(i,3) is super-diagonal (i,i+1)
!   real(rk), intent(out)          :: mf(:,:) ! LU factorization;
!                                             ! mf(i,1) is the inverse of the diagonal of U [LAPACK's 1/D]
!                                             ! mf(i,2) is the first super-diagonal of U [LAPACK's DU]
!                                             ! mf(i,3) is the second superdiagonal of U [LAPACK's DU2]
!                                             ! mf(i,4) contains L matrix factors [LAPACK's DL]
!   logical, intent(out)           :: lpiv(:) ! lpiv(i) is true if row i was interchanged with row i+1
!   logical, intent(out), optional :: fail    ! Set to .true. if decomposition fails.
!   !
!   real(rk)      :: fact, temp         ! Must match the type of mf()
    !
!   character(len=clen), save :: rcsid_tridiagonal_pivoted_m3dp_decompose_common = "$Id: tridiagonal_pivoted_m3dp_decompose_common.f90,v 1.7 2021/04/26 15:44:44 ps Exp $"
    !
    integer(ik)   :: n                  ! Dimensionality of the linear system
    integer(ik)   :: i
    !
    n = size(m,dim=1)
    if (size(m,dim=2)<3 .or. size(mf,dim=2)<4 .or. size(mf,dim=1)/=n .or. size(lpiv)/=n) then
      write (out,"(/'tridiagonal_pivoted%m3dp_decompose: Bad array dimensions in pivoted solver'/)") 
      write (out,"('   m: ',2i8)") ubound(m)
      write (out,"('  mf: ',2i8)") ubound(mf)
      write (out,"('lpiv: ',2i8)") ubound(lpiv)
      call flush_wrapper(out)
      stop 'tridiagonal_pivoted%m3dp_decompose - bad array dimensions'
    end if
    ! Copy the matrix to the output space, and initialize the pivot and second superdiagonal
    mf(:,1) = m(:,1)
    mf(:,4) = m(:,2)
    mf(:,2) = m(:,3)
    mf(:,3) = 0       ! Zero out the second superdiagonal
    mf(n,2) = 0       ! Zero out the last element of the superdiagonal as well
    lpiv(:) = .false. ! Default to no interchanges
    !
    eliminate_columns: do i=1,n-1
      if( abs(mf(i,1))>=abs(mf(i,4)) ) then
         !
         !  No row interchange required, eliminate sub-diagonal
         !  If the diagonal is zero, the subdiagonal will also vanish.
         !
         if( abs(mf(i,1))/=0 ) then
            fact      = mf(i,4) / mf(i,1)         ! subdiagonal / diagonal
            mf(i,  4) = fact                      ! Remember L factor
            mf(i+1,1) = mf(i+1,1) - fact*mf(i,2)  ! Subtract the superdiagonal from the next row
         end if
      else
         !
         !   Interchange rows I and I+1, eliminate former diagonal
         !
         fact      = mf(i,1) / mf(i,4)             ! diagonal / subdiagonal
         mf(i,1)   = mf(i,4)                       ! move former subdiagonal to the new diagonal
         mf(i,4)   = fact                          ! remember L factor
         temp      = mf(i,2)                       ! Remember former super-diagonal; this will be the diagonal for the next row
         mf(i,2)   = mf(i+1,1)                     ! Diagonal from the next row is the new super-diagonal
         mf(i+1,1) = temp - fact*mf(i+1,1)         ! Subtract new super-diagonal from the new diagonal of the next row
         !  The next two statements are meaningless for (i==n-1) [next-to-last column]
         !  However, since we initialized mf(2,n) to zero, they are also harmless.
         mf(i,3)   = mf(i+1,2)                     ! Superdiagonal of the next row becomes the second superdiagonal
         mf(i+1,2) = -fact*mf(i+1,2)               ! ... and we need to subtract it from the next superdiagonal
         !
         lpiv(i)   = .true.                        ! remember we pivoted here
      end if
    end do eliminate_columns
    !
    !  If any of the elements on the diagonal of U matrix are zero, decomposition has failed
    ! 
    if (any(mf(:,1)==0)) then
      if (present(fail)) then
        fail = .true.
        return
      else
        write (out,"('tridiagonal_pivoted%m3dp_decompose - LU decomposition failed')") 
        stop 'tridiagonal_pivoted%m3dp_decompose - failed decomposition'
      end if
    end if
    !
    !  Take the inverse, since this is what we need in the solver.
    !
    if (present(fail)) fail = .false.
    mf(:,1) = 1._rk/mf(:,1)
! end subroutine m3dp_decompose_r
