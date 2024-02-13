#!/bin/bash
#
#  Output tester for SCID-TDSE.
#  This script should not be called directly; invoke it AFTER choosing the test precision
#
if [ -z "${quiet}" -o -z "${out}" -o -z "${extraext}" -o \
     -z "${fmt_energy}" -o -z "${fmt_eigenenergy}" -o -z "${fmt_population}" -o \
     -z "${fmt_dipole}" -o -z "${fmt_pes}" ] ; then
  echo "${0}: Please do not invoke this script directly from the command line!"
  exit 99
fi
#
#  Intermediate scratch files we'll need
#
tmpref="/tmp/scid_check_general.$$.ref"
tmpout="/tmp/scid_check_general.$$.out"
trap "rm -f ${tmpref} ${tmpout}" 0
#
#  Do we have the reference output?
#
ref="${out}${extraext}"
if [ ! -r "${ref}" ] ; then
  echo "WARNING: Precision-specific reference output (${ref}) does not exist. Falling back to double precision reference"
  extraext="_ref"
  ref="${out}${extraext}"
fi
if [ ! -r "${ref}" ] ; then
  echo "ERROR: Reference output (${ref}) not found. Aborting"
  exit 1
fi
#
if { ! $quiet ; } ; then
  echo "Using reference output: ${ref}"
  echo "Using format strings: energy = ${fmt_energy}, population = ${fmt_population}, pes = ${fmt_pes}"
fi
#
#  We'll need some cosmetic editing as well!
#
function fix_zero () {
  sed -e 's/[+-]\(0\.0*\)\( \)/\1\2/g' -e 's/[+-]\(0\.0*\)$/\1/g'
  }
#
# Some files are just too large!
#
function take_100th () {
  awk 'BEGIN{nr=0}/^#/{next}{nr++}((nr%100)==1){print}' "$1"
  }
function take_1000th () {
  awk 'BEGIN{nr=0}/^#/{next}{nr++}((nr%1000)==1){print}' "$1"
  }
#
#  0. Start with the failure count being zero
#
total=0
failed=0
#
#  NaN results
#
total=$((total + 1))
nancount="$(awk 'BEGIN{c=-0}/ nan /{c++}END{print c}' IGNORECASE=1 "$out")"
if { ! $quiet ; } ; then
  echo "Test $total: $nancount lines contain NaN (not-a-number) value"
fi
if [ "${nancount}" != "0" ] ; then
  echo "Encountered $nancount NaN values. Something's wrong with the numerics."
  failed=$((failed + 1))
fi
#
#  Total energy of the initial solution
#
total=$((total + 1))
energyre_ref="$(awk '/^Energy of the atomic solution is/{printf "'${fmt_energy}' ", $7}' "$ref" | fix_zero)"
energyim_ref="$(awk '/^Energy of the atomic solution is/{printf "'${fmt_energy}' ", $8}' "$ref" | fix_zero)"
energyre_out="$(awk '/^Energy of the atomic solution is/{printf "'${fmt_energy}' ", $7}' "$out" | fix_zero)"
energyim_out="$(awk '/^Energy of the atomic solution is/{printf "'${fmt_energy}' ", $8}' "$out" | fix_zero)"
if { ! $quiet ; } ; then
  echo "Test $total: Initial energy ${energyre_out}, ${energyim_out} (expected ${energyre_ref}, ${energyim_ref})"
fi
if [ "${energyre_ref}" != "${energyre_out}" -o "${energyim_ref}" != "${energyim_out}" ] ; then
  echo "Energy(initial). Ref = ${energyre_ref}, ${energyim_ref}. Calc = ${energyre_out}, ${energyim_out}"
  failed=$((failed + 1))
fi
#
#  Right to real-space conversion factor
#
total=$((total + 1))
right2real_ref="$(awk '/^ Right to real-space wavefunction conversion factor: /{printf "'${fmt_energy}' ", $7}' "$ref" | fix_zero)"
right2real_out="$(awk '/^ Right to real-space wavefunction conversion factor: /{printf "'${fmt_energy}' ", $7}' "$out" | fix_zero)"
if { ! $quiet ; } ; then
  echo "Test $total: Right to real-space conversion factor: ${right2real_out} (expected ${right2real_ref})"
fi
if [ "${right2real_ref}" != "${right2real_out}" ] ; then
  echo "Right-to-real. Ref = ${right2real_ref}. Calc = ${right2real_out}"
  failed=$((failed + 1))
fi
#
#  Extract data from the running output (^@). 
#
function fetch_running () {
  awk '/^@ /{v=sprintf(sprintf("%s ",fmt),$(c))}
       /^@w=1 /{v="";}
       /^@w=/  {v=v sprintf(sprintf("%s ",fmt),$(c))}
       END {printf "%s",v}' c="$1" fmt="$2" "$3" | fix_zero
  }
#
#  Total norm at the end of the simulation
#
total=$((total + 1))
normre_ref="$(fetch_running 12 "${fmt_population}" "$ref")"
normim_ref="$(fetch_running 13 "${fmt_population}" "$ref")"
normre_out="$(fetch_running 12 "${fmt_population}" "$out")"
normim_out="$(fetch_running 13 "${fmt_population}" "$out")"
if { ! $quiet ; } ; then
  echo "Test $total: Final total norm ${normre_out},${normim_out} (expected ${normre_ref},${normim_ref})"
fi
if [ "${normre_ref}" != "${normre_out}" -o "${normim_ref}" != "${normim_out}" ] ; then 
  echo "Norm(final). Ref = ${normre_ref},${normim_ref}. Calc = ${normre_out},${normim_out}"
  failed=$((failed + 1))
fi
#
#  Total energy at the end of the simulation
#
total=$((total + 1))
finere_ref="$(fetch_running 15 "${fmt_energy}" "$ref")"
fineim_ref="$(fetch_running 16 "${fmt_energy}" "$ref")"
finere_out="$(fetch_running 15 "${fmt_energy}" "$out")"
fineim_out="$(fetch_running 16 "${fmt_energy}" "$out")"
if { ! $quiet ; } ; then
  echo "Test $total: Final energy ${finere_out},${fineim_out} (expected ${finere_ref},${fineim_ref})"
fi
if [ "${finere_ref}" != "${finere_out}" -o "${fineim_ref}" != "${fineim_out}" ] ; then 
  echo "Energy(final). Ref = ${finere_ref},${fineim_ref}. Calc = ${finere_out},${fineim_out}"
  failed=$((failed + 1))
fi
#
#  Dipole modulus at the end of the simulation
#
total=$((total + 1))
findip_ref="$(fetch_running 18 "${fmt_dipole}" "$ref")"
findip_out="$(fetch_running 18 "${fmt_dipole}" "$out")"
if { ! $quiet ; } ; then
  echo "Test $total: Final dipole modulus ${findip_out} (expected ${findip_ref})"
fi
if [ "${findip_ref}" != "${findip_out}" ] ; then 
  echo "Dipole(final). Ref = ${findip_ref}. Calc = ${findip_out}"
  failed=$((failed + 1))
fi
#
#  Final norm, by angular momentum channel
#
function reduce_angular () {
  awk 'BEGIN {infty=1e10;line=infty}
       /^  *Final norm, by total angular momentum/{ line = NR + 4 }
       (NR>=line)&&(NF!=4) {line=infty}
       (NR>=line) { printf " %d %d '${fmt_population}' '${fmt_population}'\n", $1, $2, $3, $4 }' "$1"
  }
total=$((total + 1))
reduce_angular "$ref" | fix_zero > "${tmpref}"
reduce_angular "$out" | fix_zero > "${tmpout}"
if { ! $quiet ; } ; then
  echo "Test $total: Final channel population. First five lines (out vs ref):"
  diff --side-by-side "${tmpout}" "${tmpref}" | head -5
fi
if { ! diff -q "${tmpref}" "${tmpout}" > /dev/null ; } ; then
  echo "Final populations of angular momentum channels differ"
  diff -c "${tmpref}" "${tmpout}"
  failed=$((failed + 1))
fi
#
#  Total population, bound states
#
total=$((total + 1))
boundpopre_ref="$(awk '/^  *Bound states: /{printf "'${fmt_population}' ",$3}' "$ref" | fix_zero)"
boundpopim_ref="$(awk '/^  *Bound states: /{printf "'${fmt_population}' ",$4}' "$ref" | fix_zero)"
boundpopre_out="$(awk '/^  *Bound states: /{printf "'${fmt_population}' ",$3}' "$out" | fix_zero)"
boundpopim_out="$(awk '/^  *Bound states: /{printf "'${fmt_population}' ",$4}' "$out" | fix_zero)"
if { ! $quiet ; } ; then
  echo "Test $total: Final bound state population ${boundpopre_out},${boundpopim_out} (expected ${boundpopre_ref},${boundpopim_ref})"
fi
if [ "${boundpopre_ref}" != "${boundpopre_out}" -o "${boundpopim_ref}" != "${boundpopim_out}" ] ; then 
  echo "Bound population. Ref = ${boundpopre_ref},${boundpopim_ref}. Calc = ${boundpopre_out},${boundpopim_out}"
  failed=$((failed + 1))
fi
#
#  Total population, continuum states
#
total=$((total + 1))
contpopre_ref="$(awk '/^  *Continuum states: /{printf "'${fmt_population}' ",$3}' "$ref" | fix_zero)"
contpopim_ref="$(awk '/^  *Continuum states: /{printf "'${fmt_population}' ",$4}' "$ref" | fix_zero)"
contpopre_out="$(awk '/^  *Continuum states: /{printf "'${fmt_population}' ",$3}' "$out" | fix_zero)"
contpopim_out="$(awk '/^  *Continuum states: /{printf "'${fmt_population}' ",$4}' "$out" | fix_zero)"
if { ! $quiet ; } ; then
  echo "Test $total: Final continuum state population ${contpopre_out},${contpopim_out} (expected ${contpopre_ref},${contpopim_ref})"
fi
if [ "${contpopre_ref}" != "${contpopre_out}" -o "${contpopim_ref}" != "${contpopim_out}" ] ; then 
  echo "Continuum population. Ref = ${contpopre_ref},${contpopim_ref}. Calc = ${contpopre_out},${contpopim_out}"
  failed=$((failed + 1))
fi
#
#  Bound/continuum population per angular momentum channel
#
function reduce_boundpop () {
  awk 'BEGIN {infty=1e10;line=infty}
       /^  *Final populations, by total angular momentum and angular momentum projection/{ line = NR + 4 }
       (NR>=line)&&(NF!=8) { line=infty }
       (NR>=line) { printf " %d %d '${fmt_population}' '${fmt_population}' '${fmt_population}' '${fmt_population}'\n", $1, $2, $5, $6, $7, $8 }' "$1"
  }
total=$((total + 1))
reduce_boundpop "$ref" | fix_zero > "${tmpref}"
reduce_boundpop "$out" | fix_zero > "${tmpout}"
if { ! $quiet ; } ; then
  echo "Test $total: Per-channel bound and continuum population. First five lines (out vs ref):"
  diff --side-by-side "${tmpout}" "${tmpref}" | head -5
fi
if { ! diff -q "${tmpref}" "${tmpout}" > /dev/null ; } ; then
  echo "Final bound/continuum populations of angular momentum channels differ"
  diff -c "${tmpref}" "${tmpout}"
  failed=$((failed + 1))
fi
#
#  Individual-state population
#
function reduce_statepop () {
  awk 'BEGIN {infty=1e10;line=infty}
       /^  *Large amplitudes of individual field-free states/{ line = NR + 4 }
       (NR>=line)&&(NF!=11) { line=infty }
       (NR>=line) { printf " %d %d %d '${fmt_eigenenergy}' '${fmt_eigenenergy}' '${fmt_population}' '${fmt_population}'\n", $1, $2, $3, $4, $5, $6, $7 }' "$1" | \
  awk '($6==0.0)&&($7==0.0){next}{print}'
  }
total=$((total + 1))
reduce_statepop "$ref" | fix_zero > "${tmpref}"
reduce_statepop "$out" | fix_zero > "${tmpout}"
if { ! $quiet ; } ; then
  echo "Test $total: Final eigenstate population. First five lines (out vs ref):"
  diff --side-by-side "${tmpout}" "${tmpref}" | head -5
fi
if { ! diff -q "${tmpref}" "${tmpout}" > /dev/null ; } ; then
  echo "Final state populations differ. First 25 lines of the output are:"
  diff -c "${tmpref}" "${tmpout}" | head -25
  failed=$((failed + 1))
fi
#
#  Check for spherical-wave photoelectron spectra being calculated
#
function reduce_spherical_waves () {
  awk '/^#/{next} (NF<5){next}
       { printf " %d %d '${fmt_energy}' '${fmt_pes}' '${fmt_pes}'\n", $1, $2, $3, $4, $5 ; }' "$1"
  }
for wave_spec_out in $(awk '/^  *Coulomb-wave phases are written to /{print $6}' "${out}") ; do
  if [ ! -z "${wave_spec_out}" ] ; then
    wave_spec_ref="${wave_spec_out}${extraext}"
    if [ ! -r "${wave_spec_ref}" ] ; then
      echo "WARNING: Reference spectrum ${wave_spec_ref} is not present. Skipping test"
    else
      total=$((total + 1))
      reduce_spherical_waves "${wave_spec_ref}" | fix_zero > "${tmpref}"   
      reduce_spherical_waves "${wave_spec_out}" | fix_zero > "${tmpout}"   
      if { ! $quiet ; } ; then
        echo "Test $total: Spherical-wave PES. First five lines (out vs ref):"
        diff --side-by-side "${tmpout}" "${tmpref}" | head -5
      fi
      if { ! diff -q "${tmpref}" "${tmpout}" > /dev/null ; } ; then
        echo "Spherical PES (${wave_spec_out}) differ. First 25 lines of the output are:"
        diff -c "${tmpref}" "${tmpout}" | head -25
        failed=$((failed + 1))
      fi
    fi
  fi
done
#
#  Check for PES spectra in OpenDX files - Volkov & Coulomb
#
function reduce_opendx_pes () {
  awk 'BEGIN {go=0} /^end/{go=1;next} (!go){next}
       { for (i=1;i<=NF;i++) { printf "'${fmt_pes}' ", $(i) } ; printf "\n" ; }' "$1"
  }
for dx_spec_out in $(awk '/^  *Generating OpenDX file for photoelectron spectrum: /{print $7}' "${out}") ; do
  if [ ! -z "${dx_spec_out}" ] ; then
    dx_spec_ref="${dx_spec_out}${extraext}"
    if [ ! -r "${dx_spec_ref}" ] ; then
      echo "WARNING: Reference spectrum ${dx_spec_ref} is not present. Skipping test"
    else
      total=$((total + 1))
      reduce_opendx_pes "${dx_spec_ref}" | fix_zero > "${tmpref}"   
      reduce_opendx_pes "${dx_spec_out}" | fix_zero > "${tmpout}"   
      if { ! $quiet ; } ; then
        echo "Test $total: OpenDX PES. First five lines (out vs ref):"
        diff --side-by-side "${tmpout}" "${tmpref}" | head -5
      fi
      if { ! diff -q "${tmpref}" "${tmpout}" > /dev/null ; } ; then
        echo "Cartesian PES (${dx_spec_out}) differ. First 25 lines of the output are:"
        diff -c "${tmpref}" "${tmpout}" | head -25
        failed=$((failed + 1))
      fi
    fi
  fi
done
#
#  Check PE spectra as tables - Volkov & Coulomb
#
function reduce_table_pes () {
  awk '/^#/{next}{ printf "'${fmt_pes}' '${fmt_pes}'\n", $7, $8 }' "$1"
  }
for table_spec_out in $(awk '/^  *Photoelectron spectrum written to /{print $5}' "${out}") ; do
  if [ ! -z "${table_spec_out}" ] ; then
    table_spec_ref="${table_spec_out}${extraext}"
    if [ ! -r "${table_spec_ref}" ] ; then
      echo "WARNING: Reference spectrum ${table_spec_ref} is not present. Skipping test"
    else
      total=$((total + 1))
      reduce_table_pes "${table_spec_ref}" | fix_zero > "${tmpref}"   
      reduce_table_pes "${table_spec_out}" | fix_zero > "${tmpout}"   
      if { ! $quiet ; } ; then
        echo "Test $total: PES First five lines (out vs ref):"
        diff --side-by-side "${tmpout}" "${tmpref}" | head -5
      fi
      if { ! diff -q "${tmpref}" "${tmpout}" > /dev/null ; } ; then
        echo "Cartesian PES (${table_spec_out}) differ. First 25 lines of the output are:"
        diff -c "${tmpref}" "${tmpout}" | head -25
        failed=$((failed + 1))
      fi
    fi
  fi
done

#
#  Check for the dipole, dipole velocity, and dipole acceleration in the detailed output
#
function reduce_detail () {
  awk '/^#/{next}
       function pd(v) { printf "'${fmt_dipole}' ", v }
       {printf "%d ", $1 ; for(i=12;i<=29;i+=2){ pd($(i)) ; } ; printf "\n"}' "$1"
  }
function detail_list () {
  awk '/^Saving detailed output to /{print $5}
       /^Saving detailed output for w.f. /{print $8}' "$1"
  }
for detail_out in $(detail_list "${out}") ; do
  if [ ! -z "${detail_out}" ] ; then
    detail_ref="${detail_out}_sum${extraext}"
    if [ ! -r "${detail_ref}" ] ; then
      echo "WARNING: Reference output ${detail_ref} is not present. Skipping test"
    else
      total=$((total + 1))
      take_1000th   "${detail_out}" > "${detail_out}_sum"
      reduce_detail "${detail_out}_sum" | fix_zero > "${tmpref}"
      reduce_detail "${detail_ref}"     | fix_zero > "${tmpout}"
      if { ! $quiet ; } ; then
        echo "Test $total: Detailed output first five lines (out vs ref):"
        diff --side-by-side "${tmpout}" "${tmpref}" | head -5
      fi
      if { ! diff -q "${tmpref}" "${tmpout}" > /dev/null ; } ; then
        echo "Detailed output (${detail_out}) differs. First 25 lines of the output are:"
        diff -c "${tmpref}" "${tmpout}" | head -25
        failed=$((failed + 1))
      fi
    fi
  fi
done
#
#  Ensemble-output check. We verify the overlap, the energy, and the dipole
#
function reduce_ensemble () {
  awk '/^#/{next}
       function pd(v) { printf "'${fmt_dipole}' ",     v }
       function pe(v) { printf "'${fmt_energy}' ",     v }
       function ps(v) { printf "'${fmt_population}' ", v }
       {printf "%d %d %d", $1, $2, $3 ; 
        ps($8) ; ps($9) ; pe($10) ; pe($11) ;
        for(i=14;i<=31;i+=2){ pd($(i)) ; } ; 
        printf "\n"}' "$1"
  }
function ensemble_list () {
  awk '/^Saving ensemble output to /{print $5}' "$1"
  }
for ensemble_out in $(ensemble_list "${out}") ; do
  if [ ! -z "${ensemble_out}" ] ; then
    ensemble_ref="${ensemble_out}_sum${extraext}"
    if [ ! -r "${ensemble_ref}" ] ; then
      echo "WARNING: Reference output ${ensemble_ref} is not present. Skipping test"
    else
      total=$((total + 1))
      take_1000th   "${ensemble_out}" > "${ensemble_out}_sum"
      reduce_ensemble "${ensemble_out}_sum" | fix_zero > "${tmpref}"
      reduce_ensemble "${ensemble_ref}"     | fix_zero > "${tmpout}"
      if { ! $quiet ; } ; then
        echo "Test $total: Detailed output first five lines (out vs ref):"
        diff --side-by-side "${tmpout}" "${tmpref}" | head -5
      fi
      if { ! diff -q "${tmpref}" "${tmpout}" > /dev/null ; } ; then
        echo "Detailed output (${ensemble_out}) differs. First 25 lines of the output are:"
        diff -c "${tmpref}" "${tmpout}" | head -25
        failed=$((failed + 1))
      fi
    fi
  fi
done
#
#  Check transition dipole calculation (separate file)
#
function reduce_dipole () {
  awk '/^#/{next}
       function pd(v) { printf "'${fmt_dipole}' ", v }
       {for(i=1;i<=6;i++) printf "%d ", $(i) ;
        for(i=7;i<=16;i+=2){ pd($(i)) ; } ; printf "\n"}' "$1"
  }
dipole_out="$(awk '/^Saving dipole matrix elements to file /{print $7}' "${out}")"
if [ ! -z "${dipole_out}" ] ; then
  dipole_ref="${dipole_out}${extraext}"
  if [ ! -r "${dipole_ref}" ] ; then
    echo "WARNING: Reference output ${dipole_ref} is not present. Skipping test"
  else
    total=$((total + 1))
    reduce_dipole "${dipole_out}" | fix_zero > "${tmpref}"
    reduce_dipole "${dipole_ref}" | fix_zero > "${tmpout}"
    if { ! $quiet ; } ; then
      echo "Test $total: transition dipole output first five lines (out vs ref):"
      diff --side-by-side "${tmpout}" "${tmpref}" | head -5
    fi
    if { ! diff -q "${tmpref}" "${tmpout}" > /dev/null ; } ; then
      echo "Transition dipole output (${dipole_out}) differs. First 25 lines of the output are:"
      diff -c "${tmpref}" "${tmpout}" | head -25
      failed=$((failed + 1))
    fi
  fi
fi
#
#  Execution times
#
walltime_ref="$(awk '/^ start  *[1-9]/{print $3}' "${ref}")"
walltime_out="$(awk '/^ start  *[1-9]/{print $3}' "${out}")"
#
if [ "$failed" -eq 0 ] ; then
  echo "OK. Runtime ref ${walltime_ref} calc ${walltime_out}"
  exit 0
else
  echo "FAIL ($failed out of $total). Runtime ref ${walltime_ref} calc ${walltime_out}"
  exit 1
fi
