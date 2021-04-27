#!/usr/bin/awk -f
BEGIN {
  nk = 0 ;
  }
/^#/ { next }
(NF<5) { next } 
/---/ { next }
{ k=$3 ; a2 = $4**2 + $5**2 ; 
  if (!khave[k]) { 
    khave[k] = 1 ; 
    klist[++nk] = k ;
    kamp[k] = 0 ;
    }
  kamp[k] += a2
  }
END {
  for (ik=1;ik<=nk;ik++) {
    k  = klist[ik] ;
    a2 = kamp[k] ;
    printf " %22.14g %22.14g\n", k, a2 ;
    }
  }
