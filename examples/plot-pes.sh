#!/bin/bash
file="$1"
l="${2:-0}"
m="${3:-0}"
wcut="${4:-1e-30}"
#
dump="${TMPDIR:-/tmp}/.plot-pes.tmp.$$"
trap "rm -f ${dump}" 0
#
awk \
'BEGIN { infty = 1e20 ; go = infty ; npeak = 0 ; gauss = 1 ; }
/Large amplitudes of individual field-free states/ { go = NR + 4 }
(NR>=go)&&(NF<6){go=infty;next;}
(NR>=go)&&($6>=wcut)&&((($1==l)&&($2==m))||(l==-1)) {
  npeak ++ ; reen[npeak] = $4 ; imen[npeak] = $5 ; wgt[npeak] = $6 ;
  }
END {
  printf "spec(e)=0 \\\n"
  for (i=1;i<=npeak;i++) {
    e =    reen[i] ;
    g = -2*imen[i] ;
    w =    wgt [i] ;
    if (g<=0) { g = 1e-15 } ;
    if (gauss) {
      printf "   +(%.14g*(2/%.14g)*sqrt(log(2.)/pi)*exp(-4*log(2.)*((e-(%.14g))/%.14g)**2))\\\n", w, g, e, g ;
      }
    else {
      printf "   +(%.14g*((%.14g)/twopi)/((e-(%.14g))**2+((%.14g)**2)/4.))\\\n", w, g, e, g ;
      }
    }
  printf "   +0\n"
  }
' l=${l} m=${m} wcut=${wcut} ${file} > ${dump}
#
{
cat ${dump}
cat <<eoi
twopi=2*4*atan2(1.,1.)
pi=4*atan2(1.,1.)
h2ev=27.2113845
set terminal x11 0 noenhanced persist
set title "${file} L=${l} M=${m}"
set xlabel "Energy, Hartree"
set ylabel "log_{10}(P(E))"
set xrange [-1:4]
set yrange [-18:1]
set samples 20000
set grid
plot log10(spec(x)) noti w l
eoi
} | tee "last.gpl" | gnuplot
