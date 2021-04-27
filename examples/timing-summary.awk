#!/usr/bin/awk -f 
#
#  Simple script to extract timings from the output of run-all.sh
#
function report (n,t) {
  printf " %-40s %10.1f\n", n, t ;
  }
/^Executing /{ name = $2 ; }
/already exists; can't run/{ name = $(NF) ; }
/^OK./   { time = $(NF) ; report(name,time) ; total += time ; }
/^FAIL./ { time = $(NF) ; report(name,time) ; total += time ; }
END { report("Total",total) ; }
