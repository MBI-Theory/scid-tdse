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
!  Mathematical and physical constants not related to numerical precision
!
module constants
  use accuracy
  implicit none
  public
  !
  character(len=clen), save :: rcsid_constants = "$Id: constants.f90,v 1.7 2021/04/26 15:44:44 ps Exp $"
  !
  !  Mathematical constants
  !
  real(rk), parameter  :: pi      = 3.1415926535897932384626433832795028841972_rk ! pi, of course
  real(xk), parameter  :: pi_xk   = 3.1415926535897932384626433832795028841972_xk ! pi, of course
  !
  !  Powers of imaginary unity
  !
  !  I**N
  complex(rk), parameter ::  ipow(0:3) = (/ (1._rk,0._rk), (0._rk,1._rk), (-1._rk,0._rk), (0._rk,-1._rk) /)
  !  (-I)**N
  complex(rk), parameter :: mipow(0:3) = (/ (1._rk,0._rk), (0._rk,-1._rk), (-1._rk,0._rk), (0._rk,1._rk) /)
  !
  !  Physical constants
  !
  real(rk), parameter :: femtosecond = 41.341373336561364_rk           ! Conversion factor: from femtoseconds to au[t]
  real(rk), parameter :: abohr       = 0.5291772083_rk                 ! Conversion factor: from au[l] (Bohrs) to Angstrom
  real(rk), parameter :: h2ev        = 27.2113845_rk                   ! Conversion factor: from Hartree to electron-volts
  real(rk), parameter :: h2cm        = 219474.6313705_rk               ! Conversion factor: from Hartree to cm^-1 (wavenumbers)
  real(rk), parameter :: k_Boltzmann = 1.0_rk/315774.65_rk             ! Conversion factor: from Kelvin to Hartree 
                                                                       !                    (aka Boltzmann constant)
  real(rk), parameter :: Hartree     = 4.35974417e-18_rk               ! Conversion factor: from au[e] (Hartree) to Joules
  real(rk), parameter :: bar         = 101325._rk                      ! Conversion factor: from bars to Pascal 
                                                                       !                    (aka standard pressure)
  real(rk), parameter :: vlight      = 137.0359991_rk                  ! Speed of light in atomic units; also the
                                                                       ! inverse of the fine structure constant
  real(rk), parameter :: R_gas       = 8.314472_rk                     ! Molar gas constant, in J/mol-K
  real(rk), parameter :: N_Avogadro  = 6.0221415e23_rk                 ! Avogadro number, in particles/mole
  real(rk), parameter :: unified_atomic_mass_unit = 1.660538782e-27_rk ! Unified atomic mass unit, in kg
  real(rk), parameter :: electron_mass_in_amu     = 5.4857990943e-4_rk ! Electron mass in u-amu
  real(rk), parameter :: au2wcm2     = 3.5e16_rk                       ! Conversion factor: from square of peak intensity 
                                                                       ! (linear polarization) to cycle-average intensity in W/cm^2
  real(rk), parameter :: au2tesla = 2.350517382e5_rk                   ! Conversion factor from "zeeman" to SI Tesla
  real(rk), parameter :: gfactor  = 2.0023193043622_rk                 ! Electron g-factor
  !
  !  Units
  !
  real(rk), parameter :: electron_mass   =  1._rk                      ! Mass of the electron in our system of units
  real(rk), parameter :: electron_charge = -1._rk                      ! Charhe of the electron in our system of units
end module constants
