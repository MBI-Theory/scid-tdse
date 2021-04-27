#!/bin/bash
#
#  Common-case output tester for SCID-TDSE. Maximum-accuracy check.
#  WARNING: Only the most basic functionality is compared to the reference output!
#
base="$(echo "$1" | sed -e 's/\.inp//' -e 's/\.out//')"
quiet="${2:-false}"
#
#  The name of the output file to check
#
out="${base}.out"
if [ ! -r "${out}" ] ; then
  echo "ERROR: Output file ${out} does not exist."
  exit 1
fi
#
#  Figure out the accuracy of the calculation - the name of the reference output will depend on it
#
#
#  Setup defaults for the double precision (this is our most common case)
#
extraext="_ref"
fmt_energy="%.10f"     # Format used for energy conversion
fmt_eigenenergy="%.8f" # Format used for eigenstate energy conversion
fmt_population="%.9f"  # Format used for population conversion
fmt_dipole="%.9f"      # Format used for dipole conversion
fmt_pes="%.7f"         # Format used for photoelectron amplitudes
#
if { grep -qs 'FATAL: Please increase sd_lmax or decrease the number of nodes' "${out}" ; } ; then
  echo "ERROR: Test case is too small for a distributed-memory run"
  exit 1
fi
#
decimals="$(awk '/^  *Real kind =/{print $6}' FS='[\t ()]+' "${out}")"
if [ -z "$decimals" ] ; then
  echo "ERROR: Can't determine working precision of the calculation."
  exit 1
fi
#
  if [ "$decimals" -le  8 ] ; then
    # 
    #  Single precision
    #
    extraext="_refsingle"
    fmt_energy="%.2f"
    fmt_eigenenergy="%.2f"
    fmt_population="%.1f"
    fmt_dipole="%.1f"
    fmt_pes="%.1f"
elif [ "$decimals" -le 16 ] ; then
    #
    #  Double precision
    #
    true ;
elif [ "$decimals" -le 34 ] ; then
    #
    #  Quadprecision
    #
    extraext="_refquad"
    fmt_energy="%.30f"
    fmt_eigenenergy="%.28f"
    fmt_population="%.28f"
    fmt_dipole="%.28f"
    fmt_pes="%.28f"
    echo "WARNING: Standard Unix tools do not support quadruple-precision floating point"
    echo "WARNING: Can't check the full dynamic range of quadruple-precision results"
else
    echo "WARNING: Number of decimals in the real type (${decimals}) is unexpected. Trying to continue anyway"
fi
#
#  Continue with the analysis; we may need to try couple different paths
#
nextscript="$(dirname "$0")/check_common.sh"
[ -x "${nextscript}" ] && source "${nextscript}"
nextscript="$(dirname "$0")/scripts/check_common.sh"
[ -x "${nextscript}" ] && source "${nextscript}"
