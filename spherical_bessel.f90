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
!  Calculate spherical bessel functions using downward recurrences,
!  as suggested in A&S section 10.5
!
!  Code is this module is not intended to be efficient, but must 
!  be numerically accurate for all inputs and real kinds.
!
 module spherical_bessel
   use accuracy
   implicit none
   private
   public besselj_table
   public rcsid_spherical_bessel
   !
   character(len=clen), save :: rcsid_spherical_bessel = "$Id: spherical_bessel.f90,v 1.7 2023/06/09 14:10:24 ps Exp $"
   !
   ! integer, parameter :: ik = selected_int_kind(15)
   ! integer, parameter :: rk = selected_real_kind(28,50)
   !
   contains
   !
   subroutine besselj_table(lmax,z,jtab)
     integer(ik), intent(in) :: lmax          ! Maximum desired order of the Bessel function
     real(rk), intent(in)    :: z             ! Argument of the Bessel function. Should not be zero!
     real(rk)                :: jtab(-1:lmax) ! Table of the Bessel functions. The index is the
                                              ! order of Bessel function; we'll start at order -1,
                                              ! defined as cos(z)/z
     !
     integer(ik) :: nmax  ! Starting order for the recursion
     integer(ik) :: n     ! Current order
     real(rk)    :: jnm1, jn, jnp1
     real(rk)    :: zm1, scl
     real(rk)    :: big
     !
     nmax = starting_besselj_order(lmax,z)
     !
     !  Evaluate recursions; remember values in the range [-1:lmax]
     !
     zm1  = 1._rk/z
     jn   = 1._rk 
     jnp1 = 0._rk
     big  = sqrt(huge(1._rk))
     downward_recursion: do n=nmax+1,0,-1
       jnm1 = real(2*n+1,kind=rk)*zm1*jn - jnp1
       ! Try rescaling if jnm1 grows too large
       if (abs(jnm1)>=big) then
         jnm1 = jnm1 / big 
         jn   = jn   / big
         jnp1 = jnp1 / big
         if (n<=lmax) jtab(n:lmax) = jtab(n:lmax) / big
       end if
       !
       if (n-1<=lmax) jtab(n-1) = jnm1
       jnp1 = jn
       jn   = jnm1
     end do downward_recursion
     !
     !  Scale the results, using the known value for j0(z) = sin(z)/z
     !
     scl  = sin(z)*zm1/jtab(0)
     jtab = scl*jtab
   end subroutine besselj_table
   !
   function starting_besselj_order(lmax,z) result(nmax)
     integer(ik), intent(in) :: lmax
     real(rk), intent(in)    :: z
     integer(ik)             :: nmax
     !
     real(rk) :: euler_e, eps, err
     !
     !  Go for something better than machine precision; the result should then be good to
     !  machine precision.
     !
     eps = 1e-3_rk * spacing(1._rk)
     !
     !  Start with a guess for N. It has to be bigger than lmax, and bigger than (E/2)*z
     !
     euler_e = exp(1._rk)
     nmax = max(lmax+1,nint(0.5_rk*euler_e*z+1))
     !
     !  Keep increasing N until estimated error is smaller than the target error
     !  We do not need to get the exact minimum starting order, so we'll grow by
     !  5% on each step.
     !
     grow_n: do 
       ! write (*,*) ' nmax = ', nmax
       err = (0.5_rk*euler_e*z/real(nmax,kind=rk))**(nmax-lmax)
       if (err<eps) exit grow_n
       if (nmax>=huge(1_ik)/3_ik) stop 'spherical_bessel%starting_besselj_order - starting order too large'
       nmax = nint(1.05_rk*real(nmax,kind=rk)+1._rk)
     end do grow_n
     ! write (*,*) ' nmax = ', nmax
   end function starting_besselj_order
 end module spherical_bessel
!
!program test
!  use spherical_bessel
!  integer(ik)           :: lmax, i
!  real(rk)              :: z
!  real(rk), allocatable :: jtab(:)
!  !
!  read(*,*) lmax, z
!  allocate (jtab(-1:lmax))
!  call besselj_table(lmax,z,jtab)
!  write (*,"((1x,i5,1x,g45.32e4))") (i, jtab(i), i=-1,lmax)
!end program test
