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
!  Hacks: Things which should not work, but do by interfering with compiler's decisions
!
module hacks
!DIR$ NOOPTIMIZE 
!
! Force Intel compiler to suppress optimizations
!
  use accuracy
  implicit none
  private
  public hack_store
  !
  !  This subroutine does nothing at all. It exists to force compiler to store
  !  a value in memory.
  !
  interface hack_store
    module procedure hack_store_integer
    module procedure hack_store_real
    module procedure hack_store_complex
    module procedure hack_store_integer_vector
    module procedure hack_store_real_vector
    module procedure hack_store_complex_vector
  end interface hack_store
  !
  !
  contains
  !
  subroutine hack_store_integer(x)
    integer(ik), intent(inout) :: x
!DIR$ ATTRIBUTES NOINLINE :: hack_store_integer
  end subroutine hack_store_integer
  !
  subroutine hack_store_real(x)
    real(rk), intent(inout) :: x
!DIR$ ATTRIBUTES NOINLINE :: hack_store_real
  end subroutine hack_store_real
  !
  subroutine hack_store_complex(x)
    complex(rk), intent(inout) :: x
!DIR$ ATTRIBUTES NOINLINE :: hack_store_complex
  end subroutine hack_store_complex
  !
  subroutine hack_store_integer_vector(x)
    integer(ik), intent(inout) :: x(:)
!DIR$ ATTRIBUTES NOINLINE :: hack_store_integer_vector
  end subroutine hack_store_integer_vector
  !
  subroutine hack_store_real_vector(x)
    real(rk), intent(inout) :: x(:)
!DIR$ ATTRIBUTES NOINLINE :: hack_store_real_vector
  end subroutine hack_store_real_vector
  !
  subroutine hack_store_complex_vector(x)
    complex(rk), intent(inout) :: x(:)
!DIR$ ATTRIBUTES NOINLINE :: hack_store_complex_vector
  end subroutine hack_store_complex_vector
end module hacks
