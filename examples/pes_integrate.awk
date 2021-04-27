#!/usr/bin/awk -f
#
#  Naive, trapezoidal-rule integration of the photoelectron spectrum.
#
#  We assume that all of r, theta, and phi grids are uniformly-spaced. 
#  The spacing is extracted from the first two points where parameters differ.
#
BEGIN {
  pi = 4.0 * atan2(1.0,1.0) ;
  dk = -1 ; dth = -1 ; dph = -1 ;
  ok = -1 ; oth = -1 ; oph = -1 ;
  }
/^#/ { next ; }
(NF>=8) {
   k  = $4 ;
   th = $5*pi/180.0 ;
   ph = $6*pi/180.0 ;
   a2 = $7**2 + $8**2 ;
   #
   s2 += k**2 * sin(th) * a2 ;
   #
   if ( (dk <0) && (ok >0) && (ok !=k ) ) { dk  = k  - ok  }
   if ( (dth<0) && (oth>0) && (oth!=th) ) { dth = th - oth }
   if ( (dph<0) && (oph>0) && (oph!=ph) ) { dph = ph - oph }
   ok  = k ;
   oth = th ;
   oph = ph ;
   }
END {
   printf "%14.7g\n", s2 * dk * dth * dph ;
   }
