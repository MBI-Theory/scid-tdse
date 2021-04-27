#!/bin/bash
make clean
( cd examples ; ./clean.sh ; cd .. )
#
now="$(date '+%Y_%m_%d')"
#
files="README.txt FILES.txt Makefile vanilla.mak make-depend.sh *.f90 *.f preprocess/ configs/ extras/ doc/ dx/"
sha256sum $(find ${files} \! -type d \! -name '*,v' -print) > MANIFEST
tar jcf distributions/scid-tdse_${now}.tbz --exclude-vcs MANIFEST ${files}
#
files_ex="examples/"
sha256sum $(find ${files} \! -type d \! -name '*,v' -print) > MANIFEST_TESTS
tar jcf distributions/scid-tdse-tests_${now}.tbz --exclude-vcs MANIFEST_TESTS ${files_ex}
