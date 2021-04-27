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
!  Report versions of all modules in this build
!
 subroutine versions
   use accuracy
   use bicg_tools
   use cap_tools
   use checkpoint_tools
   use composition_analysis
   use constants
   use coulomb_functions
   use cubic_spline
   use math
   use node_tools
   use potential_tools
   use propagator_tools
   use rotation_tools
   use sort_tools
   use spherical_bessel
   use spherical_data
   use spherical_data_initialize
   use spherical_tsurf
   use spherical_tsurf_data
   use test_tools
   use timer
   use tridiagonal_cyclic
   use tridiagonal_pivoted
   use tridiagonal_tools
   use vectorpotential_tools
   use wavefunction_tools
   !
   write (out,"(t5,a)") trim(rcsid_accuracy)
   write (out,"(t5,a)") trim(rcsid_bicg_tools)
   write (out,"(t5,a)") trim(rcsid_cap_tools)
   write (out,"(t5,a)") trim(rcsid_checkpoint_tools)
   write (out,"(t5,a)") trim(rcsid_composition_analysis)
   write (out,"(t5,a)") trim(rcsid_constants)
   write (out,"(t5,a)") trim(rcsid_coulomb_functions)
   write (out,"(t5,a)") trim(rcsid_cubic_spline)
   write (out,"(t5,a)") trim(rcsid_math)
   write (out,"(t5,a)") trim(rcsid_node_tools)
   write (out,"(t5,a)") trim(rcsid_potential_tools)
   write (out,"(t5,a)") trim(rcsid_propagator_tools)
   write (out,"(t5,a)") trim(rcsid_rotation_tools)
   write (out,"(t5,a)") trim(rcsid_sort_tools)
   write (out,"(t5,a)") trim(rcsid_spherical_bessel)
   write (out,"(t5,a)") trim(rcsid_spherical_data)
   write (out,"(t5,a)") trim(rcsid_spherical_data_initialize)
   write (out,"(t5,a)") trim(rcsid_spherical_tsurf)
   write (out,"(t5,a)") trim(rcsid_spherical_tsurf_data)
   write (out,"(t5,a)") trim(rcsid_test_tools)
   write (out,"(t5,a)") trim(rcsid_timer)
   write (out,"(t5,a)") trim(rcsid_tridiagonal_cyclic)
   write (out,"(t5,a)") trim(rcsid_tridiagonal_pivoted)
   write (out,"(t5,a)") trim(rcsid_tridiagonal_tools)
   write (out,"(t5,a)") trim(rcsid_vectorpotential_tools)
   write (out,"(t5,a)") trim(rcsid_wavefunction_tools)
 end subroutine versions
