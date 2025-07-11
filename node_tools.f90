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
!  A limited distributed-memory parallelization of scid_tdse
!
!  The parallelization scheme assumes that the individual nodes (which may be multi-threaded 
!  with OpenMP) have sufficient memory to hold all data (wavefunctions etc) on each node. 
!  Each node is expected to have uniform access to all input and output files.
!  Each node is "responsible" for a contiguous range of L values. We will try to avoid
!  synchronizing the entire wavefunction as much as possible; instead, we'll maintain a
!  "buffer" of the immediate adjacent L values on each node. This buffer is chosen to be
!  sufficient to advance the calculation on each node. The full-synchronizaton step is
!  still unavoidable for checkpoints and final analysis and reporting. 
!
!  This parallelization is not intended to scale to a large number of nodes.
!
!  The underlying transport is MPI; however, no routines outside this module should
!  access the transport protocol directly, or make additional assumptions about its 
!  properties.
!
!  The implementation is extremely naive. We deliberately refrain from using MPI
!  data types, and always communicate using the MPI_BYTE type - all we are interested
!  from the communication linrary is serving as a dumb data conduit. Going this route,
!  we give up hererogeneous execution, but gain the ability to use any Fortran datatype
!  - even those not supported by the MPI implementation.  As a side effect, we do not 
!  use reduction collectives (or any other advanced MPI features) at all.
!
!  All communications are over MPI_COMM_WORLD.
!
!  All nodes are assumed to have the same speed. This assumption is in nt_rebalance();
!  the rest of the code will (or should) accept any partitioning set up there.
!
!  We assime that the MPI implementation provides some means of copying the standard
!  input to all nodes. This has to happen -before- we do MPI initialization.
!
!  Unfortunately, MPI seems to be deliberately designed to lead to ugly Fortran
!  and unreadable code, no matter what you do.
!
module node_tools
  use accuracy
  use constants
  use hacks
  use spherical_data
  use spherical_tsurf_data
  use timer
  !
  !  MPI 3 and higher declares "use mpi" to be a legacy interface. However, with Intel 18 compilers,
  !  it allows the final binary to link MPI libraries statically. Furthermore, using the legacy
  !  interface allows integer types to be used for most operations and constants, reducing the
  !  number of conditionals we are forces to use if MPI is not present/not required.
  !
!*mp use mpi
  !
  implicit none
  private
  public rcsid_node_tools
  public nt_node_output, nt_rebalance_interval, nt_use_multinode, nt_verbose, nt_max_requests
  public nts
  public nt_initialize, nt_broadcast, nt_add, nt_max, nt_merge_all, nt_merge_borders
  public nt_rebalance, nt_rebalance_needed, nt_finalize
  public nt_force_stdin
  !
  character(len=clen), save :: rcsid_node_tools = "$Id: node_tools.f90,v 1.13 2025/07/11 15:08:35 ps Exp $"
  !
  !  All local state must be assembled in the nt_state structure, to potentiallly 
  !  allow multiple simultaneous instances of the solver at some far away future time
  !
  !  Values marked with [GBL] are synchronized across nodes
  !  Values marked with [LOC] are unique on each node
  !
  type nt_state
    logical                  :: mpi_active    ! [GBL] .True. if MPI initialization succeeded (even if we only have 1 node!)
    integer(ik)              :: n_nodes       ! [GBL] Number of nodes. Node IDs are from 1 to n_nodes
    integer(ik)              :: this_node     ! [LOC] ID of the currently-executing node, from 1 to n_nodes
                                              !       Node 1 is the "master" node. 
    integer(ik), allocatable :: lrng(:,:)     ! [GBL] First [(1,:)] and last [(2,:)] L value owned by each node
    integer(ik), allocatable :: owns(:)       ! [GBL] Node ID "owning" the block (L value) and responsible for
                                              !       updating it. Each block must be allocated to a node. The indices
                                              !       run from 0 to sd_lmax
    integer(ik), allocatable :: uses(:)       ! [GBL] Node ID "using" hte block (L value) in the boundary region
                                              !       Value of 0 means the block is not a part of a bondary region
    integer(ik), allocatable :: dist(:)       ! [LOC] Distance between each block and the nearest block "owned"
                                              !       by the currently-executing node.
    integer(ik)              :: cnt_balance   ! [GBL] Number of calls to nt_rebalance since the last rebalancing
    integer(ik)              :: lmax_balance  ! [GBL] lmax value at last rebalance
    !
    !  Asynchronous MPI queue. Use mpq_get_slot() and mpq_complete_all to manage
    !
    integer, allocatable     :: mpi_req(:)    ! [LOC] Table of MPI request IDs, used for non-blocking communications
    integer, allocatable     :: mpi_stat(:,:) ! [LOC] Table of MPI status variables, used for non-blocking communications
    integer(ik)              :: req_head      ! [LOC] Head of the FIFO request queue. This is the first free position
                                              !       in mpi_req/mpi_stat. It grows upwards.
    integer(ik)              :: req_tail      ! [LOC] Tail of the FIFO request queue. req_tail==req_head indicates an
                                              !       empty FIFO. Otherwise, req_tail points to the oldest element of
                                              !       mpi_req/mpi_stat.
    !
    !  Our load-balancing scheme relies on a very simple cost model. We assume that
    !  the cost of propagating each L block is given by:
    !
    !    cost(l) = a + b*l
    !
    !  where the constants are determined by the type of the calculation we do. For
    !  linear polarization (sd_mmin==sd_mmax), we have a=1, b=0. For general 
    !  polarization, a=1, b=2. We then distribute the available L blocks over the
    !  nodes in contiguous chunks, trying to keep the cost on each node the same.
    !
    !  In the future, we may also keep track of each node's relative performance on 
    !  previous balancing cycles, and adjust the workload accordingly. However, this
    !  is an enhancement for the future.
    !
    real(rk)                  :: cost_a
    real(rk)                  :: cost_b
  end type nt_state
  !
  !  The default is to try using multiple nodes if compiled with multi-node support enabled.
  !
!*mp logical, save          :: nt_use_multinode      = .true.  ! Enables/disables distributed memory parallel execution
!*nm logical, save          :: nt_use_multinode      = .false. ! Enables/disables distributed memory parallel execution
  integer(ik), save         :: nt_rebalance_interval = 500_ik  ! Number of time steps to perform before each rebalanancing
  integer(ik), save         :: nt_max_requests       = 128_ik  ! Maximum number of asymchronous MPI requests in flight + 1
                                                               ! This number should be high enough to saturate the communication
                                                               ! channel, but not so high so as to trigger resource exhaustion.
  character(len=clen), save :: nt_node_output        = "node_" ! For all nodes other than the master, the output is 
                                                               ! redirected to the file (nt_node_output//this_node//".out").
  integer(ik), save         :: nt_verbose            = 2       ! Debugging level for internode communications
  type(nt_state), save      :: nts
  !
  integer, save             :: mpi_integersize                 ! Size of integer(ik) in MPI "bytes". Has to be default integer
  integer, save             :: mpi_realsize                    ! Size of real(rk) in MPI "bytes"
  integer, save             :: mpi_xrealsize                   ! Size of real(xk) in MPI "bytes"
  integer, save             :: mpi_complexsize                 ! Size of complex(rk) MPI "bytes"
  !
  interface nt_broadcast
    module procedure nt_broadcast_wavefunction
    module procedure nt_broadcast_tsurf
    module procedure nt_broadcast_integer
    module procedure nt_broadcast_xreal_array
  end interface nt_broadcast
  !
  interface nt_merge_all
    module procedure nt_merge_all_wavefunction
    module procedure nt_merge_all_tsurf
  end interface nt_merge_all
  !
  interface nt_add
    module procedure nt_add_real
    module procedure nt_add_complex
    module procedure nt_add_complex_array
    module procedure nt_add_complex_array2
  end interface nt_add
  !
  interface nt_max
    module procedure nt_max_real
    module procedure nt_max_real_array
  end interface nt_max
  !
  !  There's a bit of an embuggerance here: Some MPI implementations/compilers want 
  !  the external declation below. And some consider it to be an error.
  !
  !  external MPI_IBcast, MPI_Bcast, MPI_ISend, MPI_IRecv, MPI_IGather, MPI_AllGather, MPI_WaitAll
  !
  contains
  !
  !  Externally-visible subroutines
  !
  !  See whether there is a command-line argument present in the form:
  !
  !    --scid-stdin filename
  !
  !  If it is, remap input stream to this file. The argument could appear
  !  at any position in the argument list; all other arguments will be
  !  ignored.
  !
  !  This routine relies on intrinsics which first appeared in Fortran-2003 
  !
  subroutine nt_force_stdin
    intrinsic command_argument_count, get_command_argument
    integer(ik)         :: ncmd, icmd, ios
    character(len=clen) :: cmd
    !
    ncmd = command_argument_count()
    find_stdin_name: do
      scan_arguments: do icmd=1,ncmd-1
        call get_command_argument(icmd,value=cmd)
        if (cmd=='--scid-stdin' .or. cmd=='--SCID-STDIN') then
          call get_command_argument(icmd+1,value=cmd)
          exit find_stdin_name
        end if
      end do scan_arguments
      ! The stdin name wasn't found. Return doing nothing
      return
    end do find_stdin_name
    !
    if (nt_verbose>=1) then
      write (out,"('Attaching standard input to ""',a,'""')") trim(cmd)
    end if
    open (input,form='formatted',action='read',position='rewind',status='old',file=trim(cmd),iostat=ios)
    if (ios/=0) then
      write (out,"('Error ',i0,' opening ""',a,'"" as the standard input.')") ios, trim(cmd)
      stop 'node_tools%nt_force_stdin - no file'
    end if
  end subroutine nt_force_stdin
  !
  subroutine nt_initialize
    integer(ik)           :: alloc, ios
!*mp  integer             :: ierror   ! Required to be of the default integer type
!*mp  integer             :: rank     ! MPI rank (0-based!)
!*mp  integer             :: size     ! MPI communicator size
    character(len=2*clen) :: name_buf
    !
    call TimerStart('Node init')
    !
    !  Default: single node (ie not running node-parallel). 
    !
    nts%n_nodes    = 1
    nts%mpi_active = .false.
    if (nt_use_multinode) then
!*nm  write (out,"('WARNING: MPI parallelization requested, but is not compiled-in.')")
!*nm  write (out,"('WARNING: Continuing on a single node')")
!*mp  mpi_init_block: do
!*mp    call MPI_Init(ierror)
!*mp    if (ierror/=MPI_SUCCESS) then 
!*mp      write (out,"('MPI_Init failed. Error code = ',i0)") ierror
!*mp      stop 'node_tools%nt_initialize - MPI_Init failed'
!*mp    end if
!*mp    call MPI_Comm_Rank(MPI_COMM_WORLD,rank,ierror)
!*mp    if (ierror/=MPI_SUCCESS) then
!*mp      write (out,"('MPI_Comm_Size failed. Error code = ',i0)") ierror
!*mp      stop 'node_tools%nt_initialize - MPI_Comm_Rank failed'
!*mp    end if
!*mp    call MPI_Comm_Size(MPI_COMM_WORLD,size,ierror)
!*mp    if (ierror/=MPI_SUCCESS) then
!*mp      write (out,"('MPI_Comm_Size failed. Error code = ',i0)") ierror
!*mp      stop 'node_tools%nt_initialize - MPI_Comm_Size failed'
!*mp    end if
        !
        !  In principle, the sizes of the primary types should be calculated using
        !  MPI_SizeOf macros/functions. Unfortunately, the set of types and kinds 
        !  for which these macros are available is not guaranteed to cover the types
        !  we need here. Instead, let's use our own routines, which have the advantage
        !  of working ...
        !
!!      !*mp    call MPI_SizeOf(1_ik,mpi_integersize,ierror)
!!      !*mp    call MPI_SizeOf(1.0_rk,mpi_realsize,ierror)
!!      !*mp    call MPI_SizeOf(1.0_xk,mpi_xrealsize,ierror)
!!      !*mp    mpi_complexsize = 2 * mpi_realsize ! Intel's mpi module implementation is crap, and does not include all kinds.
        mpi_integersize = ik_bytes ()
        mpi_realsize    = rk_bytes ()
        mpi_xrealsize   = xk_bytes ()
        mpi_complexsize = 2 * mpi_realsize
!*mp    nts%this_node  = rank + 1  ! MPI ranks are zero-based
!*mp    nts%n_nodes    = size
!*mp    nts%mpi_active = .true.
!*mp    exit mpi_init_block
!*mp  end do mpi_init_block
    end if
    !
    allocate (nts%owns(0:sd_lmax),nts%uses(0:sd_lmax),nts%dist(0:sd_lmax),nts%lrng(2,nts%n_nodes),stat=alloc)
    if (alloc/=0) then
      stop 'node_tools%nt_initialize - out of memory (1)'
    end if
!*mp allocate (nts%mpi_req(nt_max_requests),nts%mpi_stat(MPI_STATUS_SIZE,nt_max_requests),stat=alloc)
    if (alloc/=0) then
      stop 'node_tools%nt_initialize - out of memory (2)'
    end if
    !
    if (nts%n_nodes /= 1) then
      !
      !  Finish MPI process initialization. At the least, we need to redirect the
      !  output of processes other than the master.
      !
      if (nts%this_node/=1) then
        write (name_buf,"(a,i0,'.out')") trim(nt_node_output), nts%this_node
        open (out,form='formatted',action='write',position='rewind',status='replace',file=trim(name_buf),iostat=ios)
        if (ios/=0) then
          ! We can no longer use out - this I/O unit is in an undefined state now!
          write (0,"('Error ',i0,' redirecting output on node ',i0,' to file ',a)") ios, nts%this_node, trim(name_buf)
          stop 'node_tools%nt_initialize - redirection error'
        end if
      end if
      !
      !  Initialize non-blocking request queue
      !
!*mp  nts%mpi_req(:) = MPI_REQUEST_NULL
      nts%req_head   = 1
      nts%req_tail   = 1
      !
      !  Set up the cost model 
      !
      nts%cost_a = 1.0_rk
      if (sd_mmin==sd_mmax) then
        nts%cost_b = 0.0_rk
      else
        nts%cost_b = 2.0_rk
      end if
      !
      !  Bit of reporting
      !
      if (nt_verbose>=0) then
        write (out,"('Executing on node ',i0,' of ',i0)") nts%this_node, nts%n_nodes
        write (out,"('Cost model: a = ',f0.3,' b = ',f0.3)") nts%cost_a, nts%cost_b
      end if
      if (nt_verbose>=1) then
        write (out,"('MPI integer(ik) size is ',i0,' bytes')") mpi_integersize
        write (out,"('   MPI real(rk) size is ',i0,' bytes')") mpi_realsize
        write (out,"('MPI complex(rk) size is ',i0,' bytes')") mpi_complexsize
      end if
    else ! nts%n_nodes == 1
      if (nts%mpi_active) then
        write (out,"(/'WARNING: Running with a single MPI node. MPI capability will not be used.')")
        write (out,"( 'WARNING: If this is not the desired result, please check MPI process startup.'/)")
      end if
      !
      !  Special case: We are running on a single node; regardless of whether MPI was
      !  initialized or not, we won't be using it. All remaining nt_calls become no-ops.
      !
      nts%this_node = 1
      nts%owns(:)   = 1
      nts%uses(:)   = 0
      nts%dist(:)   = 0
      nts%lrng(1,:) = 0
      nts%lrng(2,:) = sd_lmax
    end if
    call nt_rebalance(sd_lmax) ! nt_rebalance() does nothing if running on a single node
    call TimerStop('Node init')
  end subroutine nt_initialize
  !
  subroutine nt_finalize
!*mp  integer  :: ierror   ! Required to be of the default integer type
    if (nts%mpi_active) then
      !
      !  Ignore all errors: there is nothing constructive we could do at this
      !  point if an error occurs.
      !
!*mp  call MPI_Barrier(MPI_COMM_WORLD,ierror)
!*mp  call MPI_Finalize(ierror)
    end if
  end subroutine nt_finalize
  !
  !  A few routines for dealing with the wavefunctions.
  !
  !  Distribute wavefunction from node 1 to all other nodes. The nodes do not yet
  !  know the L and R range of the wavefunction.
  !
  subroutine nt_broadcast_wavefunction(wfn)
    type(sd_wfn), intent(inout) :: wfn
    integer(ik)                 :: buf(3)
    integer                     :: count   ! MPI variables; must be default-integer
!*mp  integer                   :: ierror  ! MPI variables; must be default-integer
    integer(ik)                 :: slot
    integer(ik)                 :: lv, mv
    !
    if (nts%n_nodes==1) return
    !
    call TimerStart('Node broadcast WF')
    !
    !  There is an MPI-native way of sending the structure - but it is much more
    !  verbose, and looks uglier. We won't bother!
    !
    !  The only values which matter are on the process root (this_node==1), but
    !  there is no harm in initializing the array on all nodes - and this turns
    !  off a gfortran warning, too.
    !
    ! if (nts%this_node==1) then
      buf(1) = wfn%lmax
      buf(2) = wfn%nradial
      buf(3) = wfn%lmax_top
    ! end if
    count = size(buf) * mpi_integersize
    slot  = mpq_get_slot()
!*mp call MPI_IBcast(buf,count,MPI_BYTE,0,MPI_COMM_WORLD,nts%mpi_req(slot),ierror)
!*mp if (ierror/=MPI_SUCCESS) then
!*mp   write (out,"('nt_broadcast_wavefunction: MPI_IBcast (1) returned ',i0)") ierror
!*mp   stop 'node_tools%nt_broadcast_wavefunction - MPI_IBcast failed (1)'
!*mp end if
    !
    !  It is not a good idea to broadcast the entire wavefunction array here, for two reasons:
    !  for sd_mmin/=sd_mmax, we will be moving around a lot of garbage. Secondly, for large
    !  calculations, the total size of (wfn) may exceed the range of the default integer type.
    !  Therefore, we'll fire off a bunch of non-blocking broadcasts, each of a fairly small
    !  size. 
    !
    count = sd_nradial * sd_nspin * mpi_complexsize
    !
    !  The wavefunction arrays on the nodes are not initialized yet; transfer the complete wavefunction
    !
    bcast_mval: do mv=sd_mmin,sd_mmax
      bcast_lval: do lv=0,sd_lmax
        if (abs(mv)>lv) cycle bcast_lval
        slot = mpq_get_slot()
        ! Magic "0" below is the MPI rank of the master node (1-1=0)
!*mp    call MPI_IBcast(wfn%wfn(:,:,lv,mv),count,MPI_BYTE,0,MPI_COMM_WORLD,nts%mpi_req(slot),ierror)
!*mp    if (ierror/=MPI_SUCCESS) then
!*mp      write (out,"('nt_broadcast_wavefunction: MPI_IBcast (2) returned ',i0)") ierror
!*mp      stop 'node_tools%nt_broadcast_wavefunction - MPI_IBcast failed (2)'
!*mp    end if
      end do bcast_lval
    end do bcast_mval
    !
    call mpq_complete_all
    !
    call hack_store(buf)  ! Necessary to prevent caching of the values in buf()
    !
    wfn%lmax     = buf(1)
    wfn%nradial  = buf(2)
    wfn%lmax_top = buf(3)
    !
    if (nt_verbose>=3) then
      write (out,"('nt_broadcast_wavefunction: node = ',i0,' lmax = ',i0,' nradial = ',i0)") &
             nts%this_node, wfn%lmax, wfn%nradial
    end if
    call TimerStop('Node broadcast WF')
  end subroutine nt_broadcast_wavefunction
  !
  !  Pull wavefunction on each node from the authoritative node, synchronizing the wavefunction
  !  on all nodes.
  !
  subroutine nt_merge_all_wavefunction(wfn)
    type(sd_wfn), intent(inout) :: wfn
    !
    integer(ik)  :: lv, mv, slot, m_min, m_max
    integer      :: count, owner  ! MPI variables; must be default-integer type
!*mp  integer    :: ierror        ! ditto
    !
    if (nts%n_nodes==1) return
    !
    call TimerStart('Node merge WF')
    !
    !  We assume that all nodes have a consistent idea of the range of L and R values now
    !  (if they don't, we have bigger issues!). We still choose to transfer the full radial
    !  part - easily fixed if necessary ...
    !
    count = sd_nradial * sd_nspin * mpi_complexsize
    merge_lval: do lv=0,wfn%lmax
      m_min = max(-lv,sd_mmin)
      m_max = min( lv,sd_mmax)
      merge_mval: do mv=m_min,m_max
        slot  = mpq_get_slot()
        owner = int(nts%owns(lv)-1)
!*mp    call MPI_IBcast(wfn%wfn(:,:,lv,mv),count,MPI_BYTE,owner,MPI_COMM_WORLD,nts%mpi_req(slot),ierror)
!*mp    if (ierror/=MPI_SUCCESS) then
!*mp      write (out,"('nt_merge_all_wavefunction: MPI_IBcast returned ',i0)") ierror
!*mp      stop 'node_tools%nt_merge_all_wavefunction - MPI_IBcast failed'
!*mp    end if
      end do merge_mval
    end do merge_lval
    call mpq_complete_all
    call TimerStop('Node merge WF')
  end subroutine nt_merge_all_wavefunction
  !
  subroutine nt_merge_borders(wfn)
    type(sd_wfn), intent(inout) :: wfn
    integer(ik)  :: lv, mv, slot
    integer(ik)  :: m_min, m_max
    integer      :: count, owner, user, tag  ! MPI variables. Must be default-integer
!*mp  integer    :: ierror                   ! MPI variables. Must be default-integer
    !
    if (nts%n_nodes==1) return
    !
    call TimerStart('Node merge WF border')
    count = sd_nradial * sd_nspin * mpi_complexsize
    merge_border_lval: do lv=0,wfn%lmax
      ! For L values which do not belong to the border, there is nothing to do.
      if (nts%uses(lv)==0) cycle merge_border_lval
      ! If this node does not own and does not use this L value, there is nothing to do.
      if (all((/nts%uses(lv),nts%owns(lv)/) .ne. nts%this_node)) cycle merge_border_lval
      !
      m_min = max(-lv,sd_mmin)
      m_max = min( lv,sd_mmax)
      merge_border_mval: do mv=m_min,m_max
        slot  = mpq_get_slot()
        owner = int(nts%owns(lv)-1)
        user  = int(nts%uses(lv)-1)
        tag   = 1000*lv + mv
        if (nts%owns(lv)==nts%this_node) then
!         write (out,"(i0,' sends ',i0,', ',i0,' to ',i0)") nts%this_node, lv, mv, nts%uses(lv)
          call flush_wrapper(out)
!*mp      call MPI_ISend(wfn%wfn(:,:,lv,mv),count,MPI_BYTE,user,tag,MPI_COMM_WORLD,nts%mpi_req(slot),ierror)
!*mp      if (ierror/=MPI_SUCCESS) then
!*mp        write (out,"('nt_merge_borders: MPI_ISend returned ',i0)") ierror
!*mp        stop 'node_tools%nt_merge_border - MPI_ISend failed'
!*mp      end if
        else ! nts%uses(lv) == nts%this_node
!         write (out,"(i0,' receives ',i0,', ',i0,' from ',i0)") nts%this_node, lv, mv, nts%owns(lv)
          call flush_wrapper(out)
!*mp      call MPI_IRecv(wfn%wfn(:,:,lv,mv),count,MPI_BYTE,owner,tag,MPI_COMM_WORLD,nts%mpi_req(slot),ierror)
!*mp      if (ierror/=MPI_SUCCESS) then
!*mp        write (out,"('nt_merge_borders: MPI_IRecv returned ',i0)") ierror
!*mp        stop 'node_tools%nt_merge_border - MPI_IRecv failed'
!*mp      end if
        end if
      end do merge_border_mval
    end do merge_border_lval
    call mpq_complete_all
    call TimerStop('Node merge WF border')
  end subroutine nt_merge_borders
  !
  subroutine nt_broadcast_tsurf(ts)
    type(sts_data), intent(inout) :: ts
    !
    integer     :: count  ! MPI variables. Must be default-integer
!*mp  integer   :: ierror ! MPI variables. Must be default-integer
    integer(ik) :: slot
    integer(ik) :: ikd
    !
    if (nts%n_nodes==1) return
    !
    call TimerStart('Node broadcast TS')
    !
    !  We need to broadcast the following quantities to all nodes:
    !
    !    wfn_scale, vphase, pref
    !
    count = mpi_realsize
    slot  = mpq_get_slot()
!*mp call MPI_IBcast(ts%wfn_scale,count,MPI_BYTE,0,MPI_COMM_WORLD,nts%mpi_req(slot),ierror)
!*mp if (ierror/=MPI_SUCCESS) then
!*mp   write (out,"('nt_broadcast_tsurf: MPI_IBcast (1) returned ',i0)") ierror
!*mp   stop 'node_tools%nt_broadcast_tsurf- MPI_IBcast failed (1)'
!*mp end if
    !
    !  If we are accumulating Volkov phase, running phase and the prefactor
    !  of Volkov states needs to be broadcast. We'll use chunks (sts_kgrid_count)
    !  in size
    !
    if (allocated(ts%vphase)) then
      !
      !  We'll broadcast the Volkov phases in chunks of (sts_kgrid_count) in size
      !
      count = sts_kgrid_count * mpi_realsize
      bcast_vphase: do ikd=1,sts_dgrid_count
        slot = mpq_get_slot()
        ! Magic "0" below is the MPI rank of the master node (1-1=0)
!*mp    call MPI_IBcast(ts%vphase(:,ikd),count,MPI_BYTE,0,MPI_COMM_WORLD,nts%mpi_req(slot),ierror)
!*mp    if (ierror/=MPI_SUCCESS) then
!*mp      write (out,"('nt_broadcast_tsurf: MPI_IBcast (2) returned ',i0)") ierror
!*mp      stop 'node_tools%nt_broadcast_tsurf - MPI_IBcast failed (2)'
!*mp    end if
      end do bcast_vphase
    end if
    !
    if (allocated(ts%pref)) then
      !
      !  We'll broadcast the Volkov phases in chunks of (sts_kgrid_count) in size
      !
      count = sts_kgrid_count * mpi_complexsize
      bcast_pref: do ikd=1,sts_dgrid_count
        slot = mpq_get_slot()
        ! Magic "0" below is the MPI rank of the master node (1-1=0)
!*mp    call MPI_IBcast(ts%pref(:,ikd),count,MPI_BYTE,0,MPI_COMM_WORLD,nts%mpi_req(slot),ierror)
!*mp    if (ierror/=MPI_SUCCESS) then
!*mp      write (out,"('nt_broadcast_tsurf: MPI_IBcast (3) returned ',i0)") ierror
!*mp      stop 'node_tools%nt_broadcast_tsurf - MPI_IBcast failed (3)'
!*mp    end if
      end do bcast_pref
    end if
    !
    !  On the master node, that's all we need to do. On all other nodes, we need to zero out
    !  amplitude and fcamp as well.
    !
    if (nts%this_node/=1) then
      if (allocated(ts%amplitude)) ts%amplitude = 0
      if (allocated(ts%fcamp)    ) ts%fcamp     = 0
    end if
    !
    call mpq_complete_all
    !
    if (nt_verbose>=3) then
      write (out,"('nt_broadcast_tsurf: node = ',i0)") nts%this_node
    end if
    call TimerStop('Node broadcast TS')
  end subroutine nt_broadcast_tsurf
  !
  subroutine nt_broadcast_integer(val)
    integer(ik), intent(inout) :: val
    !
    integer     :: count  ! MPI variables. Must be default-integer
!*mp  integer   :: ierror ! MPI variables. Must be default-integer
    !
    if (nts%n_nodes==1) return
    !
    call TimerStart('Node broadcast integer')
    !
    count = mpi_integersize
!*mp call MPI_Bcast(val,count,MPI_BYTE,0,MPI_COMM_WORLD,ierror)
!*mp if (ierror/=MPI_SUCCESS) then
!*mp   write (out,"('nt_broadcast_integer: MPI_Bcast (1) returned ',i0)") ierror
!*mp   stop 'node_tools%nt_broadcast_integer MPI_Bcast failed (1)'
!*mp end if
    call TimerStop('Node broadcast integer')
  end subroutine nt_broadcast_integer
  !
  subroutine nt_broadcast_xreal_array(val)
    real(xk), intent(inout) :: val(:)
    !
    integer     :: count  ! MPI variables. Must be default-integer
!*mp  integer   :: ierror ! MPI variables. Must be default-integer
    !
    if (nts%n_nodes==1) return
    !
    call TimerStart('Node broadcast xr array')
    !
    count = size(val) * mpi_xrealsize
!*mp call MPI_Bcast(val,count,MPI_BYTE,0,MPI_COMM_WORLD,ierror)
!*mp if (ierror/=MPI_SUCCESS) then
!*mp   write (out,"('nt_broadcast_xreal_array: MPI_Bcast (1) returned ',i0)") ierror
!*mp   stop 'node_tools%nt_broadcast_xreal_array MPI_Bcast failed (1)'
!*mp end if
    call TimerStop('Node broadcast xr array')
  end subroutine nt_broadcast_xreal_array
  !
  subroutine nt_merge_all_tsurf(ts)
    type(sts_data), intent(inout) :: ts
    !
    integer                   :: count  ! MPI variables. Must be default-integer
!*mp  integer                 :: ierror ! MPI variables. Must be default-integer
    integer(ik)               :: slot, done
    integer(ik)               :: ikd, lv, mv, alloc, ic
    !
    !  mbuf is the merge buffer; it only needs to be present on the master node. 
    !  The first index  [1:sts_kgrid_count] is the k magnitude
    !  The second index [1:nts%n_nodes]     is the node index
    !  The last index   [1:nt_max_requests] is the communication-request index
    !
    complex(rk), allocatable  :: mbuf   (:,:,:)
    integer(ik)               :: slot_id(2,nt_max_requests)
    !
    if (nts%n_nodes==1) return
    !
    if (.not.allocated(ts%amplitude) .and. .not.allocated(ts%fcamp)) return
    !
    call TimerStart('Node merge TS')
    !
    if (nts%this_node==1) then
      allocate (mbuf(sts_kgrid_count,nts%n_nodes,nt_max_requests),stat=alloc)
    else
      allocate (mbuf(1,1,1),stat=alloc)  ! Just enough to make it safe to pass to a subroutine
    end if
    if (alloc/=0) then
      write (out,"('Error ',i0,' allocating reduction buffer in node_tools%nt_merge_all_tsurf')") alloc
      stop 'node_tools%nt_merge_all_tsurf - allocate failed'
    end if
    !
    !  We need to merge (add up) the values in ts%amplitude and ts%fcamp.
    !  The result is stored on the master node; once the merge is complete,
    !  the values on all other nodes must be zeroed out.
    !
    count = sts_kgrid_count * mpi_complexsize
    if (allocated(ts%amplitude)) then
      collect_amp: do ikd=1,sts_dgrid_count
        slot = mpq_get_slot(done)
        call finish_amp_reduction(done)  ! The slot number of a completed communication is in (done)
        slot_id(1,slot) = ikd            ! Rememeber what this slot is working on - we'll need it later
        ! The magic zero constant is the master process rank
!*mp    call MPI_IGather(ts%amplitude(:,ikd),count,MPI_BYTE,mbuf(:,:,slot),count,MPI_BYTE, &
!*mp                     0,MPI_COMM_WORLD,nts%mpi_req(slot),ierror)
!*mp    if (ierror/=MPI_SUCCESS) then
!*mp      write (out,"('nt_merge_all_tsurf: MPI_IGather (1) failed with code ',i0)") ierror
!*mp      stop 'node_tools%nt_merge_all_tsurf - MPI_IGather (1)'
!*mp    end if
      end do collect_amp
      !
      !  There may be some requests in flight still. Lets harvest them and finish the reduction
      !  We know that at most nt_max_requests could be flying about
      !
      mopup_collect_amp: do ic=1,nt_max_requests
        slot = mpq_get_slot(done)
        call finish_amp_reduction(done)
      end do mopup_collect_amp
    end if
    !
    if (allocated(ts%fcamp)) then
      collect_fcamp_m: do mv=sd_mmin,sd_mmax
        collect_fcamp_l: do lv=abs(mv),sd_lmax
          slot = mpq_get_slot(done)
          call finish_fcamp_reduction(done)
          slot_id(1,slot) = lv
          slot_id(2,slot) = mv
          ! The magic zero constant is the master process rank
!*mp      call MPI_IGather(ts%fcamp(:,lv,mv),count,MPI_BYTE,mbuf(:,:,slot),count,MPI_BYTE, &
!*mp                       0,MPI_COMM_WORLD,nts%mpi_req(slot),ierror)
!*mp      if (ierror/=MPI_SUCCESS) then
!*mp        write (out,"('nt_merge_all_tsurf: MPI_IGather (2) failed with code ',i0)") ierror
!*mp        stop 'node_tools%nt_merge_all_tsurf - MPI_IGather (2)'
!*mp      end if
        end do collect_fcamp_l
      end do collect_fcamp_m
      !
      !  There may be some requests in flight still. Lets harvest them and finish the reduction
      !
      mopup_collect_fcamp: do ic=1,nt_max_requests
        slot = mpq_get_slot(done)
        call finish_fcamp_reduction(done)
      end do mopup_collect_fcamp
    end if
    !
    if (nt_verbose>=3) then
      write (out,"('nt_merge_tsurf: node = ',i0)") nts%this_node
    end if
    !
    if (allocated(mbuf)) deallocate (mbuf)
    !
    call TimerStop('Node merge TS')
    !
    contains
      subroutine finish_amp_reduction(done)
        integer(ik), intent(in) :: done
        integer(ik)             :: ikd2
        !
        if (done==0) return ! Nothing completed, just return
        ikd2 = slot_id(1,done)
        if (nts%this_node==1) then
          ts%amplitude(:,ikd2) = sum(mbuf(:,:,done),dim=2)
        else
          ts%amplitude(:,ikd2) = 0
        end if
      end subroutine finish_amp_reduction
      !
      subroutine finish_fcamp_reduction(done)
        integer(ik), intent(in) :: done
        integer(ik)             :: lv2, mv2
        !
        if (done==0) return ! Nothing completed, just return
        lv2 = slot_id(1,done)
        mv2 = slot_id(2,done)
        if (nts%this_node==1) then
          ts%fcamp(:,lv2,mv2) = sum(mbuf(:,:,done),dim=2)
        else
          ts%fcamp(:,lv2,mv2) = 0
        end if
      end subroutine finish_fcamp_reduction
  end subroutine nt_merge_all_tsurf
  !
  function nt_rebalance_needed (l_max) result (go)
    integer(ik), intent(in)       :: l_max
    logical                       :: go
    !
    nts%cnt_balance = nts%cnt_balance + 1
    !
    if (l_max/=nts%lmax_balance) then
      !
      !  We have static load balancing model; since l_max haven't changed,
      !  there is nothing to do.
      !
      go = nts%cnt_balance>nt_rebalance_interval
    else
      go = .false.
    end if
  end function nt_rebalance_needed
  !
  !  Unconditionally rebalance the distributed-memory workload. Wavefunctions
  !  on all nodes must be in a consistent state before calling nt_rebalance(),
  !  or bad things will happen!
  !
  subroutine nt_rebalance(l_max)
    integer(ik), intent(in)       :: l_max
    !
    real(rk)    :: work_left             ! Total work left to allocate
    real(rk)    :: work_chunk            ! Work we need to allocate to the current node
    real(rk)    :: work_cur_l0           ! Work if the block starts at the current L
    real(rk)    :: work_cur_l1           ! Work if the block starts at L+1
    integer(ik) :: node, l_val, l_top, l_bnd, l_low
    integer(ik) :: n_own                 ! Owner node
!*mp  integer   :: ierror                ! MPI status variable
    !
    if (nts%n_nodes==1) return
    !
    call TimerStart('Node rebalance')
    nts%cnt_balance  = 0
    nts%lmax_balance = l_max
    !
    !  Reset existing data tables
    !
    nts%lrng = -1
    nts%owns =  0
    nts%uses =  0
    nts%dist = -1
    !
    !  Work unit assignment is a little tricky: we need to make sure the calculation 
    !  is as balanced as possible, while making sure that each work unit is assigned 
    !  to a node. We also need to make sure that each node is assigned at least four
    !  work units, so that we do not violate our single-source, single-consumer
    !  boundary model. If l_max is still too low to give every node a sufficiently
    !  large work chunk (which could happen early in a calculation), we'll increase
    !  it to four times the number of nodes, but not more than sd_lmax.
    !
    !  Because work units tend to increase for higher L (cost_b>=0), we'll start at
    !  the top. All L values greater than l_max belong to the top node. These work
    !  units have zero cost at the time of rebalancing.
    !
    if (sd_lmax<4*nts%n_nodes-1) then
      write (out,"(/'FATAL: Parallel execution requires that sd_lmax>=4*n_nodes-1')")
      write (out,"( 'FATAL: Please increase sd_lmax or decrease the number of nodes')")
      write (out,"(t5,'sd_lmax = ',i0)") sd_lmax
      write (out,"(t5,'n_nodes = ',i0)") nts%n_nodes
      stop 'node_tools%nt_rebalance - Too small'
    end if
    !
    l_top      = min(sd_lmax,max(l_max,4*nts%n_nodes-1))
    work_left  = block_cost(0_ik,l_top)
    nts%owns(l_top+1:sd_lmax) = nts%n_nodes
    nts%lrng(2,nts%n_nodes)   = sd_lmax
    node_assign: do node=nts%n_nodes,2,-1
      work_chunk  = work_left / node ! This is how much work we'd like to allocate
      l_val       = l_top - 4        ! Guarantee a certain number of L values in the work chunk
      if (l_val<0) then
        write (out,"(/'Work chunk allocation failed for node ',i0,' (top). l_val = ',i0,' l_top = ',i0)") node,l_val,l_top
        call report_distribution
        stop 'node_tools%nt_rebalance - work allocation failed (1)'
      end if
      work_cur_l1 = block_cost(l_val+1_ik,l_top)
      work_cur_l0 = block_cost(l_val+0_ik,l_top)
      find_l_bottom: do while (abs(work_cur_l0-work_chunk)<=abs(work_cur_l1-work_chunk)) 
        l_val = l_val - 1
        if (l_val<0) then
          write (out,"(/'Work chunk allocation failed for node ',i0,' (bottom). l_val = ',i0,' l_top = ',i0)") node,l_val,l_top
          call report_distribution
          stop 'node_tools%nt_rebalance - work allocation failed (2)'
        end if
        work_cur_l1 = work_cur_l0
        work_cur_l0 = block_cost(l_val,l_top)
      end do find_l_bottom
      !
      !  At this point, l_val + 1 is the start of the best-fit work chunk.
      !  It's work value is in work_cur_l1
      !
      nts%owns(l_val+1:l_top) = node
      nts%lrng(1,node+0) = l_val+1
      nts%lrng(2,node-1) = l_val
      l_top     = l_val
      work_left = work_left - work_cur_l1
    end do node_assign
    !
    !  Node 1 (the master) gets all remaining work units. Make sure there are
    !  enough (at least 4) units here first. Remeber that we start at l=0!
    !
    if (l_top<3) then
      write (out,"(/'Work chunk allocation failed for node 1. l_top = ',i0)") l_top
      call report_distribution
      stop 'node_tools%nt_rebalance - work allocation failed (3)'
    end if
    nts%owns(:l_top) = 1
    nts%lrng(1,1)    = 0
    !
    !  For each L value we also need to know which nodes require it in the boundary
    !  region. Potential "users" can be nodes one higher or one lower than the owner.
    !
    assign_boundary: do l_val=0,sd_lmax
      n_own = nts%owns(l_val)
      if (n_own>1) then
        !
        !  Node (n_own-1) will need access to the two lowest L values
        !
        l_bnd = nts%lrng(1,n_own)
        nts%uses(l_bnd+0) = n_own-1
        nts%uses(l_bnd+1) = n_own-1
      end if
      if (n_own<nts%n_nodes) then
        !
        !  Node (n_own+1) will need access to the two highest L values
        !
        l_bnd = nts%lrng(2,n_own)
        nts%uses(l_bnd-0) = n_own+1
        nts%uses(l_bnd-1) = n_own+1
      end if
    end do assign_boundary
    !
    !  Finally, for this node we need to know the distance from all L values to the
    !  nearest L value it "owns"
    !
    l_low = nts%lrng(1,nts%this_node)
    l_top = nts%lrng(2,nts%this_node)
    assign_distance: do l_val=0,sd_lmax
      if (l_val<l_low) then
        nts%dist(l_val) = l_low - l_val
      else if (l_val>l_top) then
        nts%dist(l_val) = l_val - l_top
      else
        nts%dist(l_val) = 0
      end if
    end do assign_distance
    !
    if (nt_verbose>=2) then
      write (out,"(/'nt_rebalance: Work chunk distribution succeeded.')")
      call report_distribution
    end if
    !
    !  This routine completely rearranges inter-node communication patterns.
    !  It requires a barrier, or bad things will happen.
    !
!*mp call MPI_Barrier(MPI_COMM_WORLD,ierror)
!*mp if (ierror/=MPI_SUCCESS) then
!*mp   write (out,"('MPI_Barrier failed. Error code = ',i0)") ierror
!*mp   stop 'node_tools%nt_rebalance - MPI_Barrier failed'
!*mp end if
    call TimerStop('Node rebalance')
    !
    contains
      !
      !  Total cost of a block of calculations
      !
      real(rk) function block_cost(lmin,lmax)
        integer(ik), intent(in) :: lmin, lmax ! Start and end (inclusing) of an L block
        !
        block_cost = 0.5_rk*(lmax-lmin+1._rk)*(2._rk*nts%cost_a + nts%cost_b*(lmax+lmin))
      end function block_cost
      !
      subroutine report_distribution
        integer(ik) :: lv, nd
        !
        write (out,"(/t5,'Allocation of L blocks per execution node. L_max = ',i0)") l_max
        write (out,"(/t5,'Properties of each L block:'/)")
        write (out,"((4(1x,a5)))") ' L ', 'Owns', 'Uses', 'Dist', &
                                   '---', '----', '----', '----'
        write (out,"((4(1x,i5)))") (lv, nts%owns(lv), nts%uses(lv), nts%dist(lv), lv=0,sd_lmax)
        write (out,"(/t5,'L ranges for each node. L_max = ',i0/)") l_max
        write (out,"((4(1x,a5),1x,a12))") &
             'Node', ' L1 ', ' Ln ', 'Count', 'Est. Cost', &
             '----', '----', '----', '-----', '---------'
        write (out,"((4(1x,i5),1x,f12.2))") &
             (nd, nts%lrng(1,nd), nts%lrng(2,nd), nts%lrng(2,nd)-nts%lrng(1,nd)+1, &
                       block_cost(nts%lrng(1,nd),nts%lrng(2,nd)), nd=1,nts%n_nodes)
        write (out,"()")
      end subroutine report_distribution
  end subroutine nt_rebalance
  !
  subroutine nt_add_real(val)
    real(rk), intent(inout) :: val
    !
!*mp real(rk)  :: buf(nts%n_nodes)
!*mp integer   :: count
!*mp  integer :: ierror
    !
    if (nts%n_nodes==1) return
    !
    call TimerStart('Node add real')
!*mp count = mpi_realsize
!*mp call MPI_AllGather(val,count,MPI_BYTE,buf,count,MPI_BYTE,MPI_COMM_WORLD,ierror)
!*mp if (ierror/=MPI_SUCCESS) then
!*mp   write (out,"('nt_add_real: MPI_AllGather failed with code ',i0)") ierror
!*mp   stop 'node_tools%nt_add_real - MPI_AllGather'
!*mp end if
!*mp val = sum(buf)
    !
    if (nt_verbose>=5) then
      write (out,"('nt_add_real: node = ',i0,' result = ',g48.24e3)") nts%this_node, val
    end if
    call TimerStop('Node add real')
  end subroutine nt_add_real
  !
  subroutine nt_add_complex(val)
    complex(rk), intent(inout) :: val
    !
!*mp complex(rk) :: buf(nts%n_nodes)
!*mp integer     :: count
!*mp integer   :: ierror
    !
    if (nts%n_nodes==1) return
    !
    call TimerStart('Node add complex')
!*mp count = mpi_complexsize
!*mp call MPI_AllGather(val,count,MPI_BYTE,buf,count,MPI_BYTE,MPI_COMM_WORLD,ierror)
!*mp if (ierror/=MPI_SUCCESS) then
!*mp   write (out,"('nt_add_complex: MPI_AllGather failed with code ',i0)") ierror
!*mp   stop 'node_tools%nt_add_complex - MPI_AllGather'
!*mp end if
!*mp val = sum(buf)
    !
    if (nt_verbose>=5) then
      write (out,"('nt_add_complex: node = ',i0,' result = ',2g38.24e3)") nts%this_node, val
    end if
    call TimerStop('Node add complex')
  end subroutine nt_add_complex
  !
  subroutine nt_add_complex_array(val)
    complex(rk), intent(inout) :: val(:)
    !
!*mp complex(rk) :: buf(size(val),nts%n_nodes)
!*mp integer     :: count
!*mp integer   :: ierror
    !
    if (nts%n_nodes==1) return
    !
    call TimerStart('Node add complex array')
!*mp count = size(val) * mpi_complexsize
!*mp call MPI_AllGather(val,count,MPI_BYTE,buf,count,MPI_BYTE,MPI_COMM_WORLD,ierror)
!*mp if (ierror/=MPI_SUCCESS) then
!*mp   write (out,"('nt_add_complex_array: MPI_AllGather failed with code ',i0)") ierror
!*mp   stop 'node_tools%nt_add_complex_array - MPI_AllGather'
!*mp end if
!*mp val = sum(buf,dim=2)
    !
    if (nt_verbose>=5) then
      write (out,"('nt_add_complex_array: node = ',i0,' result:')") nts%this_node
      write (out,"(4g38.24e3)") val
    end if
    call TimerStop('Node add complex array')
  end subroutine nt_add_complex_array
  !
  subroutine nt_add_complex_array2(val)
    complex(rk), intent(inout) :: val(:,:)
    !
!*mp complex(rk) :: buf(size(val,dim=1),size(val,dim=2),nts%n_nodes)
!*mp integer     :: count
!*mp integer     :: ierror
    !
    if (nts%n_nodes==1) return
    !
    call TimerStart('Node add complex array2')
!*mp count = size(val) * mpi_complexsize
!*mp call MPI_AllGather(val,count,MPI_BYTE,buf,count,MPI_BYTE,MPI_COMM_WORLD,ierror)
!*mp if (ierror/=MPI_SUCCESS) then
!*mp   write (out,"('nt_add_complex2: MPI_AllGather failed with code ',i0)") ierror
!*mp   stop 'node_tools%nt_add_complex_array2 - MPI_AllGather'
!*mp end if
!*mp val = sum(buf,dim=3)
    !
    if (nt_verbose>=5) then
      write (out,"('nt_add_complex_array2: node = ',i0,' result:')") nts%this_node
      write (out,"(4g38.24e3)") val
    end if
    call TimerStop('Node add complex array2')
  end subroutine nt_add_complex_array2
  !
  subroutine nt_max_real(val)
    real(rk), intent(inout) :: val
    !
    real(rk)  :: buf(nts%n_nodes)
    integer   :: count
!*mp  integer :: ierror
    !
    if (nts%n_nodes==1) return
    !
    call TimerStart('Node max real')
    count = mpi_realsize
!*mp call MPI_AllGather(val,count,MPI_BYTE,buf,count,MPI_BYTE,MPI_COMM_WORLD,ierror)
!*mp if (ierror/=MPI_SUCCESS) then
!*mp   write (out,"('nt_max_real: MPI_AllGather failed with code ',i0)") ierror
!*mp   stop 'node_tools%nt_max_real - MPI_AllGather'
!*mp end if
    !
    val = maxval(buf,dim=1)
    if (nt_verbose>=5) then
      write (out,"('nt_max_real: node = ',i0,' result = ',g38.24e3)") nts%this_node, val
    end if
    call TimerStop('Node max real')
  end subroutine nt_max_real
  !
  subroutine nt_max_real_array(val)
    real(rk), intent(inout) :: val(:)
    !
    real(rk)  :: buf(size(val),nts%n_nodes)
    integer   :: count
!*mp  integer :: ierror
    !
    if (nts%n_nodes==1) return
    !
    call TimerStart('Node max real array')
    count = size(val) * mpi_realsize
!*mp call MPI_AllGather(val,count,MPI_BYTE,buf,count,MPI_BYTE,MPI_COMM_WORLD,ierror)
!*mp if (ierror/=MPI_SUCCESS) then
!*mp   write (out,"('nt_max_real_array: MPI_AllGather failed with code ',i0)") ierror
!*mp   stop 'node_tools%nt_max_real_array - MPI_AllGather'
!*mp end if
    !
    val = maxval(buf,dim=2)
    if (nt_verbose>=5) then
      write (out,"('nt_max_real_array: node = ',i0,' result:')") nts%this_node
      write (out,"(4g38.24e3)") val
    end if
    call TimerStop('Node max real array')
  end subroutine nt_max_real_array
  !
  !  Internal subroutines
  !
  !
  !  mpq_get_slot() and mpq_complete_all() are used for managing the 
  !  non-blocking communication queue.
  !
  function mpq_get_slot(retired) result (slot)
    integer(ik), intent(out), optional :: retired  ! If present, will contail slot index which was retired
                                                   ! at this call, or zero if no requests were retired
    integer(ik)                        :: slot
!*mp  integer                          :: ierror
    !
    slot = nts%req_head
    nts%req_head = nts%req_head + 1
    if (nts%req_head>nt_max_requests) nts%req_head = 1
    !
    !  If FIFO is full, we'll need to retire the oldest request
    !
    if (present(retired)) retired = 0
    if (nts%req_head == nts%req_tail) then
!*mp  if (nts%mpi_req(nts%req_tail)/=MPI_REQUEST_NULL) then
!*mp    if (present(retired)) then
!*mp      retired = nts%req_tail
!*mp    end if
        !
        !  Note that MPI_Wait will complete immediately and without an error
        !  if the request ID is MPI_REQUEST_NULL (ie the request was never started)
        !
!*mp    call MPI_Wait(nts%mpi_req(nts%req_tail),nts%mpi_stat(:,nts%req_tail),ierror) 
!*mp    if (ierror/=MPI_SUCCESS) then 
!*mp      write (out,"('MPI_Wait returned ',i0)") ierror
!*mp      if (ierror==MPI_ERR_IN_STATUS) then
!*mp        write (out,"('Extended error status:')")
!*mp        write (out,"(t5,'source = ',i0)") nts%mpi_stat(MPI_SOURCE,nts%req_tail)
!*mp        write (out,"(t5,'   tag = ',i0)") nts%mpi_stat(MPI_TAG,nts%req_tail)
!*mp        write (out,"(t5,' error = ',i0)") nts%mpi_stat(MPI_ERROR,nts%req_tail)
!*mp      end if
!*mp      stop 'node_tools%mpq_get_slot - MPI_Wait failed'
!*mp    end if
        !
!*mp    nts%mpi_req(nts%req_tail) = MPI_REQUEST_NULL
!*mp  end if
      nts%req_tail = nts%req_tail + 1
      if (nts%req_tail>nt_max_requests) nts%req_tail = 1
    end if
  end function mpq_get_slot
  !
  subroutine mpq_complete_all
    integer     :: nreq
!*mp  integer   :: ierror
    !
    nreq = int(nt_max_requests)
!*mp call MPI_WaitAll(nreq,nts%mpi_req,nts%mpi_stat,ierror)
!*mp if (ierror/=MPI_SUCCESS) then 
       !
       !  At least some of the requests experienced errors.
       !
!*mp   write (out,"('MPI_WaitAll returned ',i0)") ierror
!*mp   if (ierror==MPI_ERR_IN_STATUS) then
         scan_requests: do while(nts%req_tail/=nts%req_head)
!*mp       if (nts%mpi_stat(MPI_ERROR,nts%req_tail)/=MPI_SUCCESS) then
!*mp         write (out,"('Extended error status:')")
!*mp         write (out,"(t5,'source = ',i0)") nts%mpi_stat(MPI_SOURCE,nts%req_tail)
!*mp         write (out,"(t5,'   tag = ',i0)") nts%mpi_stat(MPI_TAG,nts%req_tail)
!*mp         write (out,"(t5,' error = ',i0)") nts%mpi_stat(MPI_ERROR,nts%req_tail)
!*mp       end if
           nts%req_tail = nts%req_tail + 1
           if (nts%req_tail>nt_max_requests) nts%req_tail = 1
         end do scan_requests
!*mp   end if
!*mp   stop 'node_tools%mpq_complete_all - MPI_WaitAll failed'
!*mp end if
!*mp nts%mpi_req(:) = MPI_REQUEST_NULL
    nts%req_head   = 1
    nts%req_tail   = 1
  end subroutine mpq_complete_all
end module node_tools
