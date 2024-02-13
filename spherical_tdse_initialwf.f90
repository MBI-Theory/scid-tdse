!
!   SCID-TDSE: Simple 1-electron atomic TDSE solver
!   Copyright (C) 2015-2024 Serguei Patchkovskii, Serguei.Patchkovskii@mbi-berlin.de
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
!   Preparation of the initial guess
!
module spherical_tdse_initialwf
  use spherical_tdse_data
  use spherical_tdse_io
  implicit none
  private
  public rcsid_spherical_tdse_initialwf
  public prepare_one_initial_wavefunction
  !
  character(len=clen), save :: rcsid_spherical_tdse_initialwf = "$Id: spherical_tdse_initialwf.f90,v 1.1 2024/02/13 14:22:14 ps Exp $"
  !
  contains
  !
  subroutine random_initial_wfn(tag,wfn)
    character(len=*), intent(in) :: tag
    type(sd_wfn), intent(inout)  :: wfn
    integer(ik)                  :: lval, mval, sval, ipt
    real(rk)                     :: hr(2)
    !
    call random_seed()
    loop_m: do mval=sd_mmin,sd_mmax
      loop_l: do lval=abs(mval),sd_lmax
        loop_s: do sval=1,sd_nspin
          loop_pt: do ipt=1,sd_nradial
            call random_number(hr)
            wfn%wfn(:,sval,lval,mval) = cmplx(hr(1),hr(2),kind=rk)
          end do loop_pt
        end do loop_s
      end do loop_l
    end do loop_m
    write (out,"('Initial ',a,' wavefunction set to random values')") trim(tag)
  end subroutine random_initial_wfn
  !
  subroutine atomic_initial_wfn(wfn_l,wfn_r)
    type(sd_wfn), intent(inout) :: wfn_l, wfn_r
    !
    integer(ik)              :: ipt, lval, mval, ind, nvec, alloc
    complex(rk), allocatable :: evec(:,:,:)   ! Eigenvectors
    complex(rk), allocatable :: eval(:)       ! Eigenvalues
    !
    lval = initial_wfn_index(1)
    mval = initial_wfn_index(2)
    ind  = initial_wfn_index(3)
    nvec = sd_nradial
    if (initial_wfn=='single') then
      ind  = 1
      nvec = 1
    end if
    !
    write (out,"('Choosing atomic solution with L=',i0,' M=',i0,' I=',i0)") lval, mval, ind
    if (lval<0 .or. lval>sd_lmax .or. mval<sd_mmin .or. mval>sd_mmax .or. &
        abs(mval)>lval .or. ind<1 .or. ind>sd_nradial) then
      stop 'spherical_tdse_initialwf%atomic_initial_wfn - bad initial_wfn_index'
    end if
    !
    allocate (evec(sd_nradial,nvec,2),eval(nvec),stat=alloc)
    if (alloc/=0) stop 'spherical_tdse_initialwf%atomic_initial_wfn - allocation failed'
    select case (initial_wfn)
      case default
        stop 'spherical_tdse_initialwf%atomic_initial_wfn - bad initial_wfn'
      case ('atomic')
        call wt_atomic_solutions(verbose,lval,eval,evec)
      case ('single')
        eval(ind) = initial_wfn_energy
        call wt_one_atomic_solution(verbose,lval,eval(ind),evec(:,ind,:))
    end select
    !
    write (out,"('Energy of the atomic solution is ',2(g28.16e3,1x),'Hartree')") eval(ind)
    wfn_l%wfn = 0
    wfn_l%wfn(:,1,lval,mval) = evec(:,ind,1)
    wfn_r%wfn = 0
    wfn_r%wfn(:,1,lval,mval) = evec(:,ind,2)
    !
    if (verbose>=2) then 
      write (out,"(/'Initial radial wavefunction was:'/)")
      write (out,"((1x,a6,5(1x,a24)))") &
             ' I ', '  R,Bohr ', '  Re(psi_l)  ', '  Im(psi_l)  ', '  Re(psi_r)  ', '  Im(psi_r)  ', &
             '---', '---------', '-------------', '-------------', '-------------', '-------------'
      print_radial: do ipt=1,sd_nradial
        write (out,"(1x,i6,5(1x,g24.13e3))") ipt, sd_rtab(ipt), evec(ipt,ind,:)
      end do print_radial
      write (out,"()")
    end if
    !
    deallocate (evec,eval)
  end subroutine atomic_initial_wfn
  !
  subroutine prepare_one_initial_wavefunction(wfn_l,wfn_r)
    type(sd_wfn), intent(inout) :: wfn_l, wfn_r
    !
    real(rk)    :: ram
    complex(rk) :: norm(2)
    integer(ik) :: alloc
    !
    call TimerStart('Initial wavefunction')
    ram = 2*2*rk_bytes()*real(sd_nradial,kind=rk)*sd_nspin*real(sd_lmax,kind=rk)*real(sd_mmax-sd_mmin+1,kind=rk)
    write (out,"('Wavefunction requires ',f12.3,' Mbytes of memory')") ram/(1024._rk**2)
    call flush_wrapper(out)
    allocate (wfn_l%wfn(sd_nradial,sd_nspin,0:sd_lmax,sd_mmin:sd_mmax), &
              wfn_r%wfn(sd_nradial,sd_nspin,0:sd_lmax,sd_mmin:sd_mmax),stat=alloc)
    if (alloc/=0) then
      write (out,"('Allocation of wavefunction array failed, code = ',i0)") alloc
      stop 'spherical_tdse%prepare_one_initial_wavefunction - out of memory'
    end if
    !
    !  Wavefunction initialization happens on the master node. The wavefunction is broadcast
    !  to other nodes if necessary.
    !
    if (nts%this_node==1) then
      select case (initial_wfn)
        case default
          write (out,"('Initial wavefunction choice ',a,' is not recognized')") trim(initial_wfn)
          stop 'spherical_tdse%prepare_one_initial_wavefunction - bad initial_wfn'
        case ('random')
          call random_initial_wfn('left',wfn_l)
          call random_initial_wfn('right',wfn_r)
        case ('unit')
          wfn_l%wfn = 1
          wfn_r%wfn = 1
        case ('atomic','single')
          call atomic_initial_wfn(wfn_l,wfn_r)
        case ('read')
          call fetch_wavefunction(initial_wfn_file,wfn_l,wfn_r)
      end select
    end if
    call nt_broadcast(wfn_l)
    call nt_broadcast(wfn_r)
    !
    !  If we are running adaptive Lmax, we do not yet know the appropriate lmax value here.
    !  Before we can figure it out, we need to normalize
    !
    wfn_l%lmax     = sd_lmax
    wfn_r%lmax     = sd_lmax
    ! sd_tolerance_* arrays are in the right/left order
    wfn_l%sd_tol_l = sd_tolerance_l(2)
    wfn_l%sd_tol_r = sd_tolerance_r(2)
    wfn_r%sd_tol_l = sd_tolerance_l(1)
    wfn_r%sd_tol_r = sd_tolerance_r(1)
    !
    call wt_normalize(wfn_l,wfn_r,norm)
    !
    !  The next subroutine is not parallelized, and runs on the master!
    !  This is inefficient, but is only done once.
    !
    call nt_merge_all(wfn_l)
    call nt_merge_all(wfn_r)
    if (nts%this_node==1) then
      call wt_reset_lrmax(wfn_l)
      call wt_reset_lrmax(wfn_r)
      wfn_l%lmax_top = wfn_l%lmax
      wfn_r%lmax_top = wfn_r%lmax
    end if
    call nt_broadcast(wfn_l)
    call nt_broadcast(wfn_r)
    !
    write (out,"('        Initial wavefunction norm was ',g24.13e3,1x,g24.13e3)") norm(1)
    write (out,"('Right wavefunction Cartesian norm was ',g24.13e3)") real(norm(2),kind=rk)
    write (out,"('                  Left/Right Lmax was ',2i4)") wfn_l%lmax, wfn_r%lmax
    write (out,"('               Left/Right nradial was ',2i8)") wfn_l%nradial, wfn_r%nradial
    call TimerStop('Initial wavefunction')
  end subroutine prepare_one_initial_wavefunction
end module spherical_tdse_initialwf
