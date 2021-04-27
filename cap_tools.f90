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
!  Complex absorbing potential
!  These routines should not be used within a critical path; if their output is
!  needed repeatedly, then cache it.
!
module cap_tools
  use accuracy
  use constants
  use timer
  implicit none
  private
  public cap_name, cap_param
  public cap_initialize
  public cap_evaluate_potential
  public rcsid_cap_tools
  !
  character(len=clen) :: rcsid_cap_tools = "$Id: cap_tools.f90,v 1.6 2021/04/26 15:44:44 ps Exp ps $"
  !
  !  We recognize
  !
  character(len=clen), save :: cap_name     = 'manolopoulos' ! Name of the potential; can be one of:
                                                             ! 'none'         - Do not absorb
                                                             ! 'manolopoulos' - See D.E. Manolopoulos, JCP 117, 9552 (2002)
  real(rk), save            :: cap_param(2) = (/ 0.2_rk, 0.2_rk /)
                                                             ! For cap_name=='manolopoulos':
                                                             ! [1] = kmin, Minimal momentum to be absorbed
                                                             ! [2] = delta, JWKB scaling constant. Will affect the efficiency with 
                                                             !       with kmin momentum is absorbed. Default is 0.2, corresponding 
                                                             !       to 1% reflection probability at k=kmin. Smaller values of delta 
                                                             !       decrease reflection probability.

  !
  !  Parameters for the Manolopoulos CAP
  !
  real(rk), parameter :: ma_c = 2.622057554292119810464840_rk ! The position of the CAP singularity 
  real(rk), parameter :: ma_a = 1._rk - 16._rk*ma_c**(-3)
  real(rk), parameter :: ma_b = (1._rk - 17._rk*ma_c**(-3)) * ma_c**(-2)
  real(rk), save      :: ma_scale                             !
  real(rk), save      :: ma_width 
  !
  contains
  !
  subroutine cap_initialize
    real(rk) :: delta, kmin
    !
    select case (cap_name)
      case default
        write (out,"('cap_tools%cap_setup_potential: potential ',a,' is not recognized')") trim(cap_name)
        stop 'cap_tools%cap_setup_potential - unknown potential'
      case ('none')
        write (out,"(/'Using reflecting outer boundary')") 
      case ('manolopoulos')
        delta    = cap_param(2)
        kmin     = cap_param(1)
        ma_width = ma_c / ( 2._rk * delta * kmin )
        ma_scale = 0.5_rk * kmin**2 / electron_mass
        write (out,"(/'Using absorbing outer boundary from: D.E. Manolopoulos, JCP 117, 9552 (2002)')")
        write (out,"(t5,' delta = ',g24.13)") delta
        write (out,"(t5,'  kmin = ',g24.13)") kmin
        write (out,"(t5,' width = ',g24.13)") ma_width
        write (out,"(t5,' scale = ',g24.13)") ma_scale
    end select 
  end subroutine cap_initialize
  !
  function cap_evaluate_potential(r,rmax) result(v)
    real(rk), intent(in) :: r    ! Radial position where potential is needed
    real(rk), intent(in) :: rmax ! Maximum extent of the simulation region (NOT the position 
                                 ! of the last grid point; absorbing potential may diverge at the edge
                                 ! of the simulation volume).
    complex(rk)          :: v    ! Complex absorbing potential at a grid point
    !
    real(rk) :: xeff, rmin
    !
    if (r<=0._rk .or. r>=rmax) then
      stop 'cap_tools%cap_evaluate_potential - domain error (1)'
    end if
    !
    v = 0
    select case (cap_name)
      case default
        write (out,"('cap_tools%cap_evaluate_potential: potential ',a,' is not recognized')") trim(cap_name)
        stop 'cap_tools%cap_evaluate_potential - unknown potential'
      case ('none')
      case ('manolopoulos')
        rmin = rmax - ma_width
        if (rmin<=0) then
          stop 'cap_tools%cap_evaluate_potential - domain error (2)'
        end if
        if (r>=rmin) then
          xeff = (r-rmin)/ma_width
          v    = -(0._rk,1._rk)*ma_scale*maPotential(xeff)
        end if
    end select 
    ! 
  end function cap_evaluate_potential
  !
  !  The JWKB reflection-free CAP shape, in renormalized coordinates
  !
  function maPotential(x) result(v)
    real(rk), intent(in) :: x ! Dimensionless penetration coordinate
    real(rk)             :: v ! The potential
    real(rk)             :: xt
    !
    xt = ma_c * min(max(0._rk,x),1._rk-spacing(100._rk))
    v  = ma_a * xt - ma_b * xt**3 + 4._rk * (ma_c-xt)**(-2) - 4._rk * (ma_c+xt)**(-2)
  end function maPotential
  !
end module cap_tools
