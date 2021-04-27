#!/bin/bash
{
cat <<eoi
set grid
set xlabel "R, Bohr"
set ylabel "log_{10}|{/Symbol Y}|"
set title "$1"
set terminal x11 0 enhanced persist
set key right top
plot "$1" u 2:(log10(sqrt(\$3**2+\$4**2))) title "Left"  with lines ls 1, \
     "$1" u 2:(log10(sqrt(\$5**2+\$6**2))) title "Right" with lines ls 2
#
set terminal x11 1 enhanced persist
set xrange [0:0.5]
set key right bottom
replot
eoi
} | tee "last.gpl" | gnuplot 
