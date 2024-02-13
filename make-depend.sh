#!/bin/bash
#
#  Generate a simple Fortran-90 dependency file for module includes
#
for file in "$@" ; do
  base="$(echo "$file" | sed -e 's/.f90$//')"
  echo -n "${base}.o: ${file} "
  awk 'BEGIN{self["OMP_LIB"]="OMP_LIB";
             self["ISO_C_BINDING"]="ISO_C_BINDING";
             self["IEEE_ARITHMETIC"]="IEEE_ARITHMETIC";
             self["ISO_FORTRAN_ENV"]="ISO_FORTRAN_ENV"}
       /^[ \t]*module[ \t]+/{self[$2]=$2}
       /^[ \t]*use[ \t]+/{if (!($2 in self)){printf "%s.o\n", $2}}
       /^[ \t]*include[ \t]+/{printf "%s\n", substr($2,2,length($2)-2)}' "${file}" | sort -u | tr -s '\12' ' '
  echo "" ; echo ""
done
