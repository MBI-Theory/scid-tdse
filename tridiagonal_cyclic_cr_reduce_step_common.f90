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
  !
  !  A single step of cyclic reduction. Looking at the notes in:
  !  "2016_Dec_22 - Cyclic Reduction Solver Notes.pdf"
  !  is helpful if you are trying to understand what is going on.
  !
! subroutine cr_reduce_step(m,f,fail)
!   real(rk), intent(in)  :: m (:,:) ! Tridiagonal matrix to be reduced
!                                    ! First index:
!                                    !   m(1,i) = Main diagonal; element (i,i)
!                                    !   m(2,i) = Sub-diagonal; element (i+1,i)
!                                    !   m(3,i) = Super-diagonal; element (i,i+1)
!                                    ! Second index: 1 .. nm = nf + nr
!   real(rk), intent(out) :: f(:,:)  ! Reduction data.
!                                    ! First index:
!                                    !   f(1..3,i) = Coefficients for back-substitution.
!                                    !      1 = Weight of the odd resudies in the solutions (aka inverse diagonal)
!                                    !      2 = Weight of the higher-numbered even solution
!                                    !      3 = Weight of the lower-numbered even solution
!                                    !   f(4..5,i) = Coefficients for forward reduction.
!                                    !      4 = Weight of the lower-numbered odd residue
!                                    !      5 = Weight of the higher-numbered odd residue
!                                    !   f(6..8,i) = (scratch) Reduced tridiagonal matrix
!                                    !      Same format as m()
!   logical, intent(out)  :: fail
    !
    integer(ik) :: nm, nf, nr
    !
    nm = size(m,dim=1)
    nf = size(f,dim=1)
    nr = nm - nf
    !
    if (size(m,dim=2)<3 .or. size(f,dim=2)<8) stop 'tridiagonal_cyclic%cr_step - Bad argument dimensions (1)'
    if (nf/=nr .and. nf/=nr+1) stop 'tridiagonal_cyclic%cr_step - Bad argument dimensions (2)'
    !
    !  Begin by calculating reciprocals of the odd elements on the main diagonal.
    !  IF any one of these is zero, the decomposition failed.
    !
    if ( any(abs(m(1:2*nf-1:2,1))<100*tiny(1._rk)) ) then
      fail = .true.
      return
    end if
    fail = .false.
    f(1:nf,1) = 1._rk / m(1:2*nf-1:2,1)
    !
    !  Coefficients for back-substitution. 
    !
    f(1:nr,2)   = -m(1:2*nr-1:2,3) * f(1:nr,1)   ! Weight of the higher-numbered even solution
    f(1:nf-1,3) = -m(2:2*nf-2:2,2) * f(2:nf,1)   ! Weight of the lower-numbered even solution
    !
    !  Coefficients for forward rhs reduction
    !
    f(1:nr,4)   = -m(1:2*nr-1:2,2) * f(1:nr,1)   ! Weight of the lower-numbered odd residue
    f(1:nf-1,5) = -m(2:2*nf-2:2,3) * f(2:nf,1)   ! Weight of the higher-numbered odd residue
    !
    !  Reduced linear-system matrix
    !
    f(1:nr,6)   = m(2:2*nr:2,1) - m(1:2*nr-1:2,2)*m(1:2*nr-1:2,3)*f(1:nr,1)  
    f(1:nf-1,6) = f(1:nf-1,6)   - m(2:2*nf-2:2,3)*m(2:2*nf-2:2,2)*f(2:nf,1) 
    f(1:nf-1,7) =               - m(3:2*nf-1:2,2)*m(2:2*nf-2:2,2)*f(2:nf,1)   ! Sub-diagonal
    f(1:nf-1,8) =               - m(3:2*nf-1:2,3)*m(2:2*nf-2:2,3)*f(2:nf,1)   ! Super-diagonal
! end subroutine cr_reduce_step
