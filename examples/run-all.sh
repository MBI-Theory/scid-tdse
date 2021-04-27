#!/bin/bash
#
# mode can be "cheap", "medium", or "all"
#
mode="${1:-"cheap"}"
#
passed=""
failed=""
dontknow=""
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
    export OMP_STACKSIZE=500M
    # Please note that if there are enough physical cores to support all running threads
    # (both OMP and MPI) setting OMP_WAIT_POLICY=passive and/or I_MPI_WAIT_MODE=on will
    # likely reduce the performance. The difference is quite dramatic on AMD Zen nodes.
    export OMP_WAIT_POLICY=passive
    export OMP_NESTED=false
    export I_MPI_WAIT_MODE=1
    export I_MPI_DEBUG=5
    # export I_MPI_FABRICS=shm:ofi  # May be needed by Intel MPI 2020 and later
    export I_MPI_OFI_PROVIDER=shm
    ulimit -s 1024000
      ../spherical_tdse.x < "${inp}" > "${out}" 2>&1 
    # Intel MPI
    # mpirun -genvall -np 2 -s all $(pwd)/../spherical_tdse.x < "${inp}" > "${out}" 2>&1
    # OpenMPI
    # /usr/lib64/mpi/gcc/openmpi/bin/mpirun \
    #       -mca btl ^tcp \
    #       -mca mpi_show_mca_params all \
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
echo "Cheap tests (expected runtime < 1 minute each)"
for inp in hydrogen_1S_2P0_uniform.inp hydrogen_1S_2P0_uniform_restart.inp hydrogen_1S_2P0.inp \
           hydrogen_2P0_ion.inp hydrogen_2P0_ion_restart.inp helium_GJG75.inp \
           helium_triplet_spline_linear.inp helium_triplet_vptable_linear.inp \
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
           hydrogen_1S_hhg_spline.inp \
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
