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
! subroutine m3u_decompose_r(m,mf,fail)
!   real(rk), intent(in)           :: m (:,:) ! Tridiagonal matrix, assumed to be diagonally dominant
!   real(rk), intent(out)          :: mf(:,:) ! Factorized tridiagonal matrix
!   logical, intent(out), optional :: fail    ! Set to .true. if decomposition fails; 
!                                             ! if fail is absent, abort on decomposition failure.
!   !
!   character(len=clen), save :: rcsid_tridiagonal_tools_m3u_decompose_common = "$Id: tridiagonal_tools_m3u_decompose_common.f90,v 1.5 2021/04/26 15:44:44 ps Exp $"
!   real(rk)    :: denom
    integer(ik) :: i, sz
    logical     :: l_fail
    !
    sz = size(m,dim=1)
    if (size(m,dim=2)<3 .or. size(mf,dim=2)<3 .or. size(mf,dim=1)/=sz) then
      stop 'tridiagonal_tools%m3u_decompose_common - bad input sizes'
    end if
    if (present(fail)) fail = .false.
    !
    l_fail = .false.
    catch_block: do
      denom   = m(1,1)
      if ( abs(denom)<=100*tiny(abs(denom)) ) then
        l_fail = .true.
        exit catch_block
      end if
      mf(1,1) = 1._rk/denom
      mf(1,2) = 0._rk
      mf(1,3) = m(1,3)*mf(1,1)
      factor_m3d: do i=2,sz-1
        denom   = m(i,1)-m(i-1,2)*mf(i-1,3)
        if ( abs(denom)<=100*tiny(abs(denom)) ) then
          l_fail = .true.
          exit catch_block
        end if
        mf(i,1) = 1._rk/denom
        mf(i,2) = -m(i-1,2)*mf(i,1)
        mf(i,3) = m(i,3)*mf(i,1)
      end do factor_m3d
      if (sz<2) exit catch_block
      denom    = m(sz,1)-m(sz-1,2)*mf(sz-1,3)
      if ( abs(denom)<=100*tiny(abs(denom)) ) then
        l_fail = .true.
        exit catch_block
      end if
      mf(sz,1) = 1._rk/denom
      mf(sz,2) = -m(sz-1,2)*mf(sz,1)
      mf(sz,3) = 0._rk
      exit catch_block
    end do catch_block
    if (present(fail)) then
      fail = l_fail
    else if (l_fail) then
      write (out,"('Fatal error in m3u_decompose_common: denominator ',g34.16e3,' is too small.')") denom
      stop 'tridiagonal_tools%m3u_decompose_common - decomposition failed'
    end if
    !
! end subroutine m3u_decompose_r
