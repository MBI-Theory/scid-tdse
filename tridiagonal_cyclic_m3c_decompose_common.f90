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
! subroutine m3c_decompose(m,mf,fail)
!   real(rk), intent(in)          :: m (:,:) ! Tridiagonal matrix, assumed to be diagonally dominant
!   real(rk), intent(out)         :: mf(:,:) ! Factorized tridiagonal matrix
!   logical, intent(out),optional :: fail    ! Set to .true. if decomposition fails; 
!                                            ! if fail is absent, abort on decomposition failure.
    !
    integer(ik) :: nm        ! Current number of variables
    integer(ik) :: nf        ! Number of variables eliminated in this round
    integer(ik) :: nr        ! Number of variables in the reduced system of equations
    integer(ik) :: ip        ! Postion of the decomposition data in the mf() array
    integer(ik) :: ip0       ! Postion of the the reduced linear system in the mf() array
    integer(ik) :: ipx       ! Largest position in the reduction buffer
    logical     :: lfail
    !
    nm  = size(m,dim=1)
    if (nm<1) stop 'tridiagonal_cyclic%m3c_decompose - must have at least one linear equation!'
    ip  = 1
    ipx = size(mf,dim=1)
    !
    catch_fail: do
      !
      !  First reduction step is special: The source matrix is in the input array.
      !
      nr = nm/2 
      nf = nm - nr
      if (ip+nf>ipx) stop 'tridiagonal_cyclic%m3c_decompose - blown the reduction buffer (1)'
      call cr_reduce_step(m,mf(ip:ip+nf-1,:),lfail)
      if (lfail) exit catch_fail
      !
      ip0 = ip
      ip  = ip + nf  ! Advance decomposition buffer position
      nm  = nr       ! Now our linear system is half the size
      !
      recursive_decomposition: do 
        nr = nm/2 
        nf = nm - nr
        if (ip+nf-1>ipx) stop 'tridiagonal_cyclic%m3c_decompose - blown the reduction buffer (2)'
        call cr_reduce_step(mf(ip0:ip0+nm-1,6:8),mf(ip:ip+nf-1,:),lfail)
        if (lfail) exit catch_fail ! Decomposition failed
        if (nr==0) exit catch_fail ! Decomposition completed successfully
        ip0 = ip
        ip  = ip + nf              ! Advance decomposition buffer position
        nm  = nr                   ! Now our linear system is half the size
      end do recursive_decomposition
      return
    end do catch_fail
    !
    if (lfail) then
      if (.not.present(fail)) stop 'tridiagonal_cyclic%m3c_decompose - decomposition failed'
    end if
    if (present(fail)) fail = lfail
    return
! end subroutine m3c_decompose
