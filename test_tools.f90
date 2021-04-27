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
!  Accuracy/convergence checks; not directly necessary for time propagation.
!
module test_tools
  use accuracy
  use constants
  use timer
  use spherical_data
  use potential_tools
  use tridiagonal_tools
  use propagator_tools
  use wavefunction_tools
  use node_tools
  use lapack
  use sort_tools
  implicit none
  private
  public derivatives_test, fieldfree_test, imaginary_propagation_test
  public rcsid_test_tools
  !
  character(len=clen), save :: rcsid_test_tools = "$Id: test_tools.f90,v 1.29 2021/04/26 15:44:44 ps Exp ps $"
  !
  contains
  !
  !  Accuracy test/debugging aid. 
  !  Calculate first and second derivatives of a function known analytically,
  !  go gauge accuracy of our grid and correctness of implicit derivative
  !  operators obtained in sd_initialize()
  !
  subroutine derivatives_test(verbose)
    integer(ik), intent(in) :: verbose ! Level of the output
    !
    real(rk) :: func (sd_nradial)             ! Function values at grid points
    real(rk) :: ref_g(sd_nradial)             ! Reference (analytical) gradient at grid points
    real(rk) :: ref_l(sd_nradial)             ! Reference (analytical) laplacian at grid points
    real(rk) :: grad (sd_nradial)             ! Numerical gradient at grid points
    real(rk) :: lap  (sd_nradial)             ! Numerical laplacian at grid points
    real(rk) :: tmp  (sd_nradial)             ! Temporary array for implicit operator evaluation
    real(rk) :: scr  (sd_nradial,m3d_sc_size) ! Scratch for iterative refinement of linear equations solutions
    !
    call TimerStart('Radial derivatives test')
    write (out,"(/'Testing derivatives in the L=0 channel (unnormalized 1S function)'/)")
    func  = test_1S(sd_rtab(:sd_nradial))
    ref_g = grad_1S(sd_rtab(:sd_nradial))
    ref_l =  lap_1S(sd_rtab(:sd_nradial))
    call m3d_multiply(sd_d1n_l0,func,tmp)
    call m3d_solve(sd_m1n_l0,sd_m1nf_l0,sd_m1np_l0,tmp,grad,scr)
    call m3d_multiply(sd_d2n_l0,func,tmp)
    call m3d_solve(sd_m2n_l0,sd_m2nf_l0,sd_m2np_l0,tmp,lap,scr) 
    call compare
    !
    write (out,"(/'Testing derivatives in the L>0 channel (unnormalized 2P function)'/)")
    func  = test_2P(sd_rtab(:sd_nradial))
    ref_g = grad_2P(sd_rtab(:sd_nradial))
    ref_l =  lap_2P(sd_rtab(:sd_nradial))
    call m3d_multiply(sd_d1n_lx,func,tmp)
    call m3d_solve(sd_m1n_lx,sd_m1nf_lx,sd_m1np_lx,tmp,grad,scr)
    call m3d_multiply(sd_d2n_lx,func,tmp)
    call m3d_solve(sd_m2n_lx,sd_m2nf_lx,sd_m2np_lx,tmp,lap,scr) 
    call compare
    call TimerStop('Radial derivatives test')
    !
    contains
    !
    subroutine compare
      integer(ik) :: ipt
      real(rk)    :: abserr_g, abserr_l, relerr_g, relerr_l
      real(rk)    :: max_abserr_g, max_abserr_l, max_relerr_g, max_relerr_l
      integer(ik) :: pos_maxabs_g, pos_maxabs_l, pos_maxrel_g, pos_maxrel_l
      real(rk)    :: eps_grad, eps_lap
      !
      if (verbose>=1) then
        write (out,"((1x,a4,6(1x,a14)))") &
               ' pt ', ' R, Bohr ', ' Func ', ' An.Grad. ', ' An.Lapl. ', ' Err.Grad. ', ' Err.Lapl. ', &
               '----', '---------', '------', '----------', '----------', '-----------', '-----------'
      end if
      !
      !  We would like to avoid divergent error estimates where quantities drop to nearly zero
      !  Therefore, we'll shift the denominator by the spacing of the largest value of the quantity.
      !
      eps_grad = 100._rk*spacing(maxval(abs(ref_g)))
      eps_lap  = 100._rk*spacing(maxval(abs(ref_l)))
      !
      max_abserr_g = 0 ; max_abserr_l = 0 ; pos_maxabs_g = 0 ; pos_maxabs_l = 0 ;
      max_relerr_g = 0 ; max_relerr_l = 0 ; pos_maxrel_g = 0 ; pos_maxrel_l = 0 ;
      scan_points: do ipt=1,sd_nradial
        abserr_g = abs(grad(ipt) - ref_g(ipt))
        abserr_l = abs(lap (ipt) - ref_l(ipt))
        relerr_g = abserr_g / (abs(ref_g(ipt)) + eps_grad)
        relerr_l = abserr_l / (abs(ref_l(ipt)) + eps_lap )
        if (abserr_g>max_abserr_g) then ; max_abserr_g = abserr_g ; pos_maxabs_g = ipt ; end if
        if (abserr_l>max_abserr_l) then ; max_abserr_l = abserr_l ; pos_maxabs_l = ipt ; end if
        if (relerr_g>max_relerr_g) then ; max_relerr_g = relerr_g ; pos_maxrel_g = ipt ; end if
        if (relerr_l>max_relerr_l) then ; max_relerr_l = relerr_l ; pos_maxrel_l = ipt ; end if
        if (verbose>=1) then
          write (out,"(1x,i4,6(1x,g14.7))") &
                 ipt, sd_rtab(ipt), func(ipt), ref_g(ipt), ref_l(ipt), grad(ipt)-ref_g(ipt), lap(ipt)-ref_l(ipt)
        end if
      end do scan_points
      if (verbose>=1) write (out,"()")
      !
      call report_error('absolute','gradient', max_abserr_g,pos_maxabs_g,ref_g,grad)
      call report_error('relative','gradient', max_relerr_g,pos_maxrel_g,ref_g,grad)
      call report_error('absolute','laplacian',max_abserr_l,pos_maxabs_l,ref_l,lap )
      call report_error('relative','laplacian',max_relerr_l,pos_maxrel_l,ref_l,lap )
      call flush_wrapper(out)
    end subroutine compare
    !
    subroutine report_error(code,name,err,pos,ref,val)
      character(len=*), intent(in) :: code   ! Type of the error
      character(len=*), intent(in) :: name   ! Name of the quantity
      real(rk), intent(in)         :: err    ! Value of the error
      integer(ik), intent(in)      :: pos    ! Grid position
      real(rk), intent(in)         :: ref(:) ! Table of the reference analytical values
      real(rk), intent(in)         :: val(:) ! Table of the calculated numerical values
      !
      write (out,"(     ' Max ',a8,' error in the ',a9,' is ',g24.13)") code, name, err
      write (out,"('       it is found at radial grid point ',i6)") pos
      write (out,"('                   radial coordinate is ',g24.13)") sd_rtab(pos)
      write (out,"('        function value at this point is ',g24.13)") func(pos)
      write (out,"(   '  analytical ',a9,' at this point is ',g24.13)") name, ref(pos)
      write (out,"(   '   numerical ',a9,' at this point is ',g24.13)") name, val(pos)
      write (out,"()")
    end subroutine report_error
    !
    !  Note that test functions include the overall r prefactor, to compensate
    !  for the implicit 1/r in the Ansatz
    !
    elemental real(rk) function test_1S(r)
      real(rk), intent(in) :: r ! Radial coordinate
      test_1S = r * exp(-sd_rgrid_zeta*r)
    end function test_1S
    elemental real(rk) function grad_1S(r)
      real(rk), intent(in) :: r ! Radial coordinate
      grad_1S = (1-r*sd_rgrid_zeta) * exp(-sd_rgrid_zeta*r)
    end function grad_1S
    elemental real(rk) function lap_1S(r)
      real(rk), intent(in) :: r ! Radial coordinate
      lap_1S = sd_rgrid_zeta*(-2+r*sd_rgrid_zeta) * exp(-sd_rgrid_zeta*r)
    end function lap_1S
    !
    elemental real(rk) function test_2P(r)
      real(rk), intent(in) :: r ! Radial coordinate
      test_2P = r**2 * exp(-0.5_rk*sd_rgrid_zeta*r)
    end function test_2P
    elemental real(rk) function grad_2P(r)
      real(rk), intent(in) :: r ! Radial coordinate
      grad_2P = 0.5_rk * r * (4-r*sd_rgrid_zeta) * exp(-0.5_rk*sd_rgrid_zeta*r)
    end function grad_2P
    elemental real(rk) function lap_2P(r)
      real(rk), intent(in) :: r ! Radial coordinate
      lap_2P = 0.25_rk * (8-r*sd_rgrid_zeta*(8-r*sd_rgrid_zeta)) * exp(-0.5_rk*sd_rgrid_zeta*r)
    end function lap_2P
  end subroutine derivatives_test
  !
  !  Calculation of the field-free eigenspectrum
  !
  subroutine fieldfree_test(verbose)
    integer(ik), intent(in)            :: verbose        ! Level of the output
    !
    integer(ik) :: lval, iev, nrep
    complex(rk) :: block_evec(sd_nradial,sd_nradial,2) ! Eigenvectors (used to check for bad grids)
    complex(rk) :: block_eval(sd_nradial)              ! Eigenvalues
    complex(rk) :: all_eval  (sd_nradial,0:sd_lmax)    ! All eigenvalues; sorted in increading real part order
    !
    call TimerStart('Field-free solutions')
    !$omp parallel do schedule(guided) default(none) &
    !$omp& private(lval,block_evec,block_eval) &
    !$omp& shared(sd_lmax,all_eval,verbose)
    scan_l_channels: do lval=0,sd_lmax
      call wt_atomic_solutions(verbose,lval,block_eval,block_evec)
      all_eval(:,lval) = block_eval(:)
      call check_effective_grid(verbose,lval,block_evec)
    end do scan_l_channels
    !$omp end parallel do
    !
    write (out,"((1x,a4,1x,a4,1x,2(1x,a24)))") &
               ' L= ', '  i ', '  Re(eps)  ', '  Im(eps)  ', &
               '----', '----', '-----------', '-----------'
    report_l_channels: do lval=0,sd_lmax
       nrep = sd_nradial
       if (verbose<3) nrep = min(sd_nradial,20)
       if (verbose<2) nrep = min(sd_nradial,10)
       if (verbose<1) nrep = min(sd_nradial,2)
       report_eigenvalues: do iev=1,nrep
         write (out,"(1x,i4,1x,i4,1x,2(1x,g24.13e3))") lval, iev, all_eval(iev,lval)
       end do report_eigenvalues
       write (out,"()")
    end do report_l_channels
    write (out,"()")
    call flush_wrapper(out)
    !
    call TimerStop('Field-free solutions')
  end subroutine fieldfree_test
  !
  subroutine check_effective_grid(verbose,lval,evec)
    integer(ik), intent(in) :: verbose     ! Level of the output
    integer(ik), intent(in) :: lval        ! Angular momentum L
    complex(rk), intent(in) :: evec(:,:,:) ! Eigenvectors
    !
    integer(ik) :: is, ipt
    real(rk)    :: rhol_p, rhor_p   ! Left and right "density" at the previous point
    real(rk)    :: rhol_c, rhor_c   ! Ditto, current point
    real(rk)    :: p1, p2, eps
    !
    scan_states: do is=1,sd_nradial
      eps     = spacing(maxval(abs(evec(:,is,2))**2))
      rhol_p  = abs(evec(1,is,1))**2
      rhor_p  = abs(evec(1,is,2))**2
      scan_points: do ipt=2,sd_nradial
        rhol_c = abs(evec(ipt,is,1))**2
        rhor_c = abs(evec(ipt,is,2))**2
        !
        ! Compare the ratio rhol/rhor, trying to avoid division by zero
        !
        if (rhor_c>=eps .and. rhor_p>=eps) then
          p1 = rhol_p*rhor_c
          p2 = rhor_p*rhol_c
          if (abs(p1-p2)>0.5_rk*(p1+p2)) then
            if (verbose>=0) then
              write (out,"(/'Possible grid problem detected: L= ',i0,' state ',i0,' grid point ',i0,' R= ',g24.13)") &
                     lval, is, ipt, sd_rtab(ipt)
              write (out,"('Previous point: <L,L>= ',g24.13,' <R,R>= ',g24.13,' rat= ',g24.13)") rhol_p, rhor_p, rhol_p/rhor_p
              write (out,"(' Current point: <L,L>= ',g24.13,' <R,R>= ',g24.13,' rat= ',g24.13)") rhol_c, rhor_c, rhol_c/rhor_c
              if (verbose<=2) then
                write (out,"(' ... skipping further warnings for this L'/)") 
                exit scan_states
              end if
            end if
          end if
        end if
        !
        rhol_p = rhol_c
        rhor_p = rhor_c
      end do scan_points
    end do scan_states
  end subroutine check_effective_grid
  !
  !  Propagate wavefunction in imaginary time under the action of atomic Hamiltonian
  !
  subroutine imaginary_propagation_test(verbose,psi_l,psi_r,apot,dt,nstep)
    integer(ik), intent(in)     :: verbose ! How much output to produce
    type(sd_wfn), intent(inout) :: psi_l   ! Left wavefunction to propagate
    type(sd_wfn), intent(inout) :: psi_r   ! Right wavefunction to propagate
    real(xk), intent(in)        :: apot    ! Vector-potential, possibly in higher accuracy
    real(xk), intent(in)        :: dt      ! Magnitude of the time step, possibly in higher accuracy
    integer(ik), intent(in)     :: nstep   ! Number of imaginary time step to perform
    !
    integer(ik) :: step
    complex(rk) :: norm(2), normx
    complex(rk) :: energy(2)
    complex(xk) :: cdt
    !
    call TimerStart('Imaginary propagation')
    !
    call wt_normalize(psi_l,psi_r,norm) 
    call nt_merge_borders(psi_l)
    call nt_merge_borders(psi_r)
    call wt_energy(psi_l,psi_r,apot=apot,energy=energy,norm=normx)
    write (out,"('Initial norm= ',g24.14,1x,g24.13,' en= ',g24.13,1x,g24.13)") norm(1), energy(1)
    !
    cdt = -(0._rk,1._rk)*dt
    imaginary_propagate: do step=1,nstep
      call nt_merge_borders(psi_r)
      call pt_fwd_laser_n(psi_r,apot,0.5_xk*cdt)
      call pt_fwd_atomic_n(psi_r,cdt)
      call pt_rev_atomic_n(psi_r,cdt)
      call nt_merge_borders(psi_r)
      call pt_rev_laser_n(psi_r,apot,0.5_xk*cdt)
      !
      call nt_merge_borders(psi_l)
      call pt_fwd_laser_t(psi_l,apot,-0.5_xk*cdt)
      call pt_fwd_atomic_t(psi_l,-cdt)
      call pt_rev_atomic_t(psi_l,-cdt)
      call nt_merge_borders(psi_l)
      call pt_rev_laser_t(psi_l,apot,-0.5_xk*cdt)
      !
      call wt_normalize(psi_l,psi_r,norm) 
      !
      if (sd_adaptive) then
        call wt_update_lrmax(verbose,psi_r)
        call wt_update_lrmax(verbose,psi_l)
      end if
      !
      if (nt_rebalance_needed(max(psi_l%lmax,psi_r%lmax))) then
        call nt_merge_all(psi_r)
        call nt_merge_all(psi_l)
        call nt_rebalance(max(psi_l%lmax,psi_r%lmax))
      end if
      !
      if (mod(step,100)==1 .or. step==nstep) then
        call nt_merge_borders(psi_l)
        call nt_merge_borders(psi_r)
        call wt_energy(psi_l,psi_r,apot=apot,energy=energy,norm=normx)
        write (out,"('Time step ',i0,' norm= ',g34.24e3,1x,g24.13e3,' en= ',g34.23e3,1x,g34.23e3)") step, norm(1), energy(1)
        call flush_wrapper(out)
      end if
    end do imaginary_propagate
    call TimerStop('Imaginary propagation')
  end subroutine imaginary_propagation_test
end module test_tools
