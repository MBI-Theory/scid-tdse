#!/bin/bash
#
# mode can be "single file_name", "cheap", "medium", or "all"
#
mode="${1:-"cheap"}"
test="${2:-"hydrogen_1S_2P0_uniform.inp"}"
#
passed=""
failed=""
dontknow=""
#
#  MacOS requires a wrapper to let program run at full speed
#
wrapper=""
[ "$(uname -s)" = "Darwin" ] && wrapper="caffeinate -i "
#
function run_test () {
  local inp="$1" ;
  local out ;
  local trc ;
  local chk ;

  out="$(echo "${inp}" | sed -e 's/\.inp/.out/')"
  trc="$(echo "${inp}" | sed -e 's/\.inp/.trace/')"
  if [ -r "${out}" ] ; then
    echo "${out} already exists; can't run ${inp}"
  else
    echo "Executing ${inp}"
    export MALLOC_TRACE="$trc"  # Will do nothing unless mtrace() call is activated in spherical_tdse.f90
                                # Use "mtrace" script (part of glibc-utils package) to analyze the traces
    export HUGETLB_MORECORE=thp
    # export OMP_STACKSIZE=32768
    export OMP_STACKSIZE=3072
    # Please note that if there are enough physical cores to support all running threads
    # (both OMP and MPI) setting OMP_WAIT_POLICY=passive and/or I_MPI_WAIT_MODE=on will
    # likely reduce the performance. The difference is quite dramatic on AMD Zen nodes.
    export OMP_WAIT_POLICY=passive
    export OMP_NESTED=false
    export I_MPI_WAIT_MODE=1
    export I_MPI_DEBUG=5
    # export I_MPI_FABRICS=shm:ofi  # May be needed by Intel MPI 2020 and later
    export I_MPI_OFI_PROVIDER=shm
    # ulimit -s 32768
    ulimit -s 3072
    ${wrapper} ../spherical_tdse.x < "${inp}" > "${out}" 2>&1 
    # Intel MPI
    # mpirun -genvall -np 2 -s all $(pwd)/../spherical_tdse.x < "${inp}" > "${out}" 2>&1
    # OpenMPI
    #       -mca mpi_show_mca_params all  <- Using this option causes OpenMPI4 to segfault. Classy.
    # /usr/lib64/mpi/gcc/openmpi/bin/mpirun \
    #       -mca btl ^tcp \
    #       -bind-to none -np 2 $(pwd)/../spherical_tdse.x --scid-stdin "${inp}" > "${out}" 2>&1
  fi
  chk="$(echo "${inp}" | sed -e 's/\.inp/.chk/')"
  if [ -x "$chk" ] ; then
    ./"${chk}" "${inp}" true
    if [ $? -eq 0 ] ; then
      passed="$passed $inp"
    else
      failed="$failed $inp"
    fi
  else
    dontknow="$dontknow $inp"
  fi
  }
#
function summary () {
  [ ! -z "$1"        ] && echo "$1"
  [ ! -z "$passed"   ] && echo "Tests PASSED: $passed"
  [ ! -z "$failed"   ] && echo "Tests FAILED: $failed"
  [ ! -z "$dontknow" ] && echo "Can't CHECK: $dontknow"
  }
#
if [ "$mode" == "single" ] ; then
  run_test "$test"
  summary "Singe-test run"
  exit 0
fi
#
if [ "$mode" == "mpi" ] ; then
  run_test "argon_3P1m_ell_ckpt_mpi.inp"
  run_test "argon_3P1m_ell_rstrt_mpi.inp"
  summary "Minimal MPI tests"
  exit 0
fi
#
echo "Cheap tests (expected runtime < 1 minute each)"
for inp in hydrogen_1S_2P0_uniform.inp hydrogen_1S_2P0_uniform_restart.inp hydrogen_1S_2P0.inp \
           hydrogen_2P0_ion.inp hydrogen_2P0_ion_restart.inp helium_GJG75.inp \
           helium_triplet_spline_linear.inp helium_triplet_vptable_linear.inp \
           argon_hhg_full.inp argon_hhg_fakeleft.inp argon_hhg_fakeleft2.inp \
           hydrogen_TME.inp \
           ; do
  run_test $inp
done
summary ""
[ "$mode" == "cheap" ] && exit 0
#
echo ""
echo "Intermediate tests (expected runtime < 10 minutes each)"
for inp in helium_1S_adiabatic.inp argon_3P1_cooper.inp argon_3P1_offcooper.inp \
           argon_3P1m_circ_l.inp argon_3P1m_circ_r.inp hydrogen_1S_hhg_linear.inp \
           hydrogen_2P0_sfi_tsurf.inp argon_3P1m_ell_ckpt_mpi.inp argon_3P1m_ell_rstrt_mpi.inp \
           hydrogen_1S_hhg_spline.inp argon_circ-hhg_full.inp argon_circ-hhg_fakeleft2.inp \
           hydrogen_static.inp lithium_ensemble.inp lithium_ensemble_restart.inp \
           ; do
  run_test $inp
done
summary
[ "$mode" != "all" ] && exit 0
#
echo ""
echo "Expensive tests (expected runtime > 1 hour each)"
for inp in hydrogen_2P0_sfi.inp hydrogen_1S_hhg_elliptical.inp hydrogen_1S_hhg_circular.inp \
           hydrogen_1S_hhg_ell_refsol.inp \
           ; do
  run_test $inp
done
summary
#
