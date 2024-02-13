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
!   Wavefunction I/O subroutines
!
module spherical_tdse_io
  use spherical_tdse_data
  implicit none
  private
  public rcsid_spherical_tdse_io
  public visualize_wavefunctions
  public dump_all_wavefunctions
  public fetch_wavefunction
  !
  character(len=clen), save :: rcsid_spherical_tdse_io = "$Id: spherical_tdse_io.f90,v 1.2 2024/02/13 16:10:04 ps Exp $"
  !
  contains
  !
  subroutine visualize_wavefunctions(wfn_l,wfn_r,tsurf,prefix,iseq,timestep,th,ph)
    type(sd_wfn), intent(inout)   :: wfn_l, wfn_r
    type(sts_data), intent(inout) :: tsurf
    character(len=*), intent(in)  :: prefix
    integer(ik), intent(inout)    :: iseq
    integer(ik), intent(in)       :: timestep
    real(xk), intent(in)          :: th, ph    ! Rotation angles for the current localm coordinate system
    !
    integer(ik)          :: lval, mval, sval, iu_temp
    character(len=clen)  :: filename
    real(rk)             :: rm(3,3)
    complex(rk)          :: l_lap(sd_nradial), l_grad(sd_nradial)
    complex(rk)          :: r_lap(sd_nradial), r_grad(sd_nradial)
    complex(rk)          :: tmp(sd_nradial), scr(sd_nradial,m3d_sc_size)
    !
    if (prefix==' ') return
    !
    call TimerStart('Visualize wavefunction')
    !
    if (skip_left_propagation .and. fake_left_propagation) then
      call wt_fake_left(wfn_l,wfn_r)
    end if
    call nt_merge_all(wfn_l)
    call nt_merge_all(wfn_r)
    if (nts%this_node/=1) then
      call TimerStop('Visualize wavefunction')
      return
    end if
    write (filename,"(a,i10.10,'.data')") trim(prefix), iseq
    iseq = iseq + 1
    open(newunit=iu_temp,form='formatted',recl=256,action='write',position='rewind',status='replace',file=trim(filename))
    write (iu_temp,"(6(1x,i0))") sd_nradial, sd_nspin, sd_lmax, sd_mmin, sd_mmax, timestep
    !
    call MathRotationMatrix((/real(ph,kind=rk),real(th,kind=rk),0._rk/),rm)
    write (iu_temp,"(3(1x,g27.18e3))") rm(1,:)
    write (iu_temp,"(3(1x,g27.18e3))") rm(2,:)
    write (iu_temp,"(3(1x,g27.18e3))") rm(3,:)
    !
    write (iu_temp,"(5(1x,g22.13e3))") sd_rtab(1:sd_nradial)
    !
    dump_m_channels: do mval=sd_mmin,sd_mmax
      dump_l_channels: do lval=0,sd_lmax
        dump_s_channels: do sval=1,sd_nspin
          if (lval==0) then
            !$omp parallel sections default(none) &
            !$omp&  shared(sd_d2n_l0,sd_m2n_l0,sd_m2np_l0,sd_m2nf_l0) &
            !$omp&  shared(sd_d1n_l0,sd_m1n_l0,sd_m1np_l0,sd_m1nf_l0) &
            !$omp&  shared(wfn_l,wfn_r,sval,lval,mval,l_lap,r_lap,l_grad,r_grad) &
            !$omp&  private(tmp,scr)
            ! Left laplacian
            !$omp section
            call m3d_multiply(sd_d2n_l0,wfn_l%wfn(:,sval,lval,mval),tmp)
            call m3d_solve(sd_m2n_l0,sd_m2nf_l0,sd_m2np_l0,tmp,l_lap,scr)
            ! Right laplacian
            !$omp section
            call m3d_multiply(sd_d2n_l0,wfn_r%wfn(:,sval,lval,mval),tmp)
            call m3d_solve(sd_m2n_l0,sd_m2nf_l0,sd_m2np_l0,tmp,r_lap,scr)
            ! Left gradient
            !$omp section
            call m3d_multiply(sd_d1n_l0,wfn_l%wfn(:,sval,lval,mval),tmp)
            call m3d_solve(sd_m1n_l0,sd_m1nf_l0,sd_m1np_l0,tmp,l_grad,scr)
            ! Right gradient
            !$omp section
            call m3d_multiply(sd_d1n_l0,wfn_r%wfn(:,sval,lval,mval),tmp)
            call m3d_solve(sd_m1n_l0,sd_m1nf_l0,sd_m1np_l0,tmp,r_grad,scr)
            !$omp end parallel sections
          else ! lval /= 0
            !$omp parallel sections default(none) &
            !$omp&  shared(sd_d2n_lx,sd_m2n_lx,sd_m2np_lx,sd_m2nf_lx) &
            !$omp&  shared(sd_d1n_lx,sd_m1n_lx,sd_m1np_lx,sd_m1nf_lx) &
            !$omp&  shared(wfn_l,wfn_r,sval,lval,mval,l_lap,r_lap,l_grad,r_grad) &
            !$omp&  private(tmp,scr)
            ! Left laplacian
            !$omp section
            call m3d_multiply(sd_d2n_lx,wfn_l%wfn(:,sval,lval,mval),tmp)
            call m3d_solve(sd_m2n_lx,sd_m2nf_lx,sd_m2np_lx,tmp,l_lap,scr)
            ! Right laplacian
            !$omp section
            call m3d_multiply(sd_d2n_lx,wfn_r%wfn(:,sval,lval,mval),tmp)
            call m3d_solve(sd_m2n_lx,sd_m2nf_lx,sd_m2np_lx,tmp,r_lap,scr)
            ! Left gradient
            !$omp section
            call m3d_multiply(sd_d1n_lx,wfn_l%wfn(:,sval,lval,mval),tmp)
            call m3d_solve(sd_m1n_lx,sd_m1nf_lx,sd_m1np_lx,tmp,l_grad,scr)
            ! Right gradient
            !$omp section
            call m3d_multiply(sd_d1n_lx,wfn_r%wfn(:,sval,lval,mval),tmp)
            call m3d_solve(sd_m1n_lx,sd_m1nf_lx,sd_m1np_lx,tmp,r_grad,scr)
            !$omp end parallel sections
          end if
          write (iu_temp,"(5(1x,g22.13e3))") (1._rk/tsurf%wfn_scale) * wfn_l%wfn(1:sd_nradial,sval,lval,mval)
          write (iu_temp,"(5(1x,g22.13e3))") (1._rk/tsurf%wfn_scale) * l_grad
          write (iu_temp,"(5(1x,g22.13e3))") (1._rk/tsurf%wfn_scale) * l_lap
          write (iu_temp,"(5(1x,g22.13e3))") tsurf%wfn_scale * wfn_r%wfn(1:sd_nradial,sval,lval,mval)
          write (iu_temp,"(5(1x,g22.13e3))") tsurf%wfn_scale * r_grad
          write (iu_temp,"(5(1x,g22.13e3))") tsurf%wfn_scale * r_lap
        end do dump_s_channels
      end do dump_l_channels
    end do dump_m_channels
    close(iu_temp)
    call TimerStop('Visualize wavefunction')
  end subroutine visualize_wavefunctions
  !
  subroutine dump_one_wavefunction(prefix,wfn_l,wfn_r)
    character(len=*), intent(in) :: prefix
    type(sd_wfn), intent(inout)  :: wfn_l, wfn_r
    integer(ik)                  :: lval, mval, sval, ir, ios, iu_temp
    character(len=clen)          :: filename
    !
    call TimerStart('Dump wavefunction')
    dump_l_channels: do lval=0,sd_lmax
      dump_m_channels: do mval=max(sd_mmin,-lval),min(lval,sd_mmax)
        if (mval>=0) then
          write (filename,"(a,'-L',i0.3,'-M+',i0.4)") trim(prefix), lval, mval
        else
          write (filename,"(a,'-L',i0.3,'-M-',i0.4)") trim(prefix), lval, -mval
        end if
        open(newunit=iu_temp,form='formatted',recl=256,action='write',position='rewind', &
             status='replace',file=trim(filename),iostat=ios)
        if (ios/=0) then
          write (out,"('WARNING: File ',a,' cannot be opened for writing. Code = ',i0)") trim(filename), ios
          cycle dump_m_channels
        end if
        write (iu_temp,"('#',a2,1x,a24,2(2x,a34,1x,a34))") &
              'S', ' R, Bohr ', ' Re[left wfn] ', ' Im[left wfn] ', ' Re[right wfn] ', ' Im[right wfn] '
        dump_s_channels: do sval=1,sd_nspin
          dump_radial: do ir=1,sd_nradial
            write (iu_temp,"(1x,i2,1x,g24.13,2(2x,g34.21e3,1x,g34.21e3))") &
                   sval, sd_rtab(ir), wfn_l%wfn(ir,sval,lval,mval), wfn_r%wfn(ir,sval,lval,mval)
          end do dump_radial
        end do dump_s_channels
        close(iu_temp)
      end do dump_m_channels
    end do dump_l_channels
    call TimerStop('Dump wavefunction')
  end subroutine dump_one_wavefunction
  !
  subroutine dump_all_wavefunctions(tag,names)
    character(len=*), intent(in) :: tag       ! What kind of dump is this 
    character(len=*), intent(in) :: names(:)  ! Names of the file prefixes, one per ensemble member
    !
    integer(ik) :: ie
    !
    dump_ensemble: do ie=1,ensemble_size
      if (names(ie)==' ') cycle dump_ensemble
      write (out,"(/'Dumping ',a,' wavefunction to disk, prefix = ',a/)") trim(tag), trim(names(ie))
      call flush_wrapper(out)
      call nt_merge_all(wfns_l(ie))
      call nt_merge_all(wfns_r(ie))
      call flush_wrapper(out)
      if (nts%this_node==1) then
        call dump_one_wavefunction(names(ie),wfns_l(ie),wfns_r(ie))
      end if
    end do dump_ensemble
  end subroutine dump_all_wavefunctions
  !
  subroutine fetch_wavefunction(prefix,wfn_l,wfn_r)
    character(len=*), intent(in) :: prefix
    type(sd_wfn), intent(inout)  :: wfn_l, wfn_r
    integer(ik)                  :: lval, mval, sval, ir, iu_temp
    character(len=clen)          :: filename
    character(len=1)             :: buf
    integer(ik)                  :: line, ios
    integer(ik)                  :: stmp
    real(rk)                     :: rtmp(5)
    !
    call TimerStart('Fetch wavefunction')
    fetch_l_channels: do lval=0,sd_lmax
      fetch_m_channels: do mval=max(sd_mmin,-lval),min(lval,sd_mmax)
        if (mval>=0) then
          write (filename,"(a,'-L',i0.3,'-M+',i0.4)") trim(prefix), lval, mval
        else
          write (filename,"(a,'-L',i0.3,'-M-',i0.4)") trim(prefix), lval, -mval
        end if
        open(newunit=iu_temp,form='formatted',recl=256,action='read',position='rewind',status='old',file=trim(filename))
        line = 1
        read (iu_temp,"(a1)",iostat=ios) buf
        if (ios/=0) then
          write (out,"('Error ',i0,' skipping line ',i0,' of file ',a)") ios, line, trim(filename)
          stop 'spherical_tdse%fetch_wavefunctions - skip error'
        end if
        fetch_s_channels: do sval=1,sd_nspin
          fetch_radial: do ir=1,sd_nradial
            line = line + 1
            read (iu_temp,*,iostat=ios) stmp, rtmp
            if (ios/=0) then
              write (out,"('Error ',i0,' reading line ',i0,' of file ',a)") ios, line, trim(filename)
              stop 'spherical_tdse%fetch_wavefunctions - read error'
            end if
            if (stmp/=sval .or. abs(sd_rtab(ir)-rtmp(1))>1e-10_rk*max(1._rk,sd_rtab(ir))) then
              write (out,"('Grid mismatch reading line ',i0,' of file ',a)") line, trim(filename)
              stop 'spherical_tdse%fetch_wavefunctions - data error'
            end if
            wfn_l%wfn(ir,sval,lval,mval) = cmplx(rtmp(2),rtmp(3),kind=rk)
            wfn_r%wfn(ir,sval,lval,mval) = cmplx(rtmp(4),rtmp(5),kind=rk)
          end do fetch_radial
        end do fetch_s_channels
        close(iu_temp)
      end do fetch_m_channels
    end do fetch_l_channels
    call TimerStop('Fetch wavefunction')
  end subroutine fetch_wavefunction
end module spherical_tdse_io
