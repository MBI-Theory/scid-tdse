#!/bin/bash
f1="$1"
f2="$2"
lval="${3:-"0"}"
mval="${4:-"0"}"
#
{
cat <<eoi
twopi=2*4*atan2(1.,1.)
pi=4*atan2(1.,1.)
h2ev=27.2113845
set terminal x11 0 noenhanced persist
set title "${file}"
set xlabel "Energy, Hartree"
set ylabel "log_{10}(|P(k(E))|)"
set grid
plot "< ./combine_coulomb_waves.awk '${f1}'" u (\$1**2/2):(log10(\$2)) ti "$f1 Total" w l ls 1 lw 1.5, \
     "< awk '(\$1==l)&&(\$2==m)' l=${lval} m=${mval} '${f1}'" u (\$3**2/2):(log10(\$4**2+\$5**2)) ti "$f1 L=${lval} M=${mval}" w l ls 2 lw 1.5, \
     "< ./combine_coulomb_waves.awk '${f2}'" u (\$1**2/2):(log10(\$2)) ti "$f2 Total" w l ls 3 lw 0.5, \
     "< awk '(\$1==l)&&(\$2==m)' l=${lval} m=${mval} '${f2}'" u (\$3**2/2):(log10(\$4**2+\$5**2)) ti "$f2 L=${lval} M=${mval}" w l ls 4 lw 0.5
eoi
} | tee "last.gpl" | gnuplot
