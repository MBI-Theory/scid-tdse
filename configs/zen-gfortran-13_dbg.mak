BUILD_ID :="Debugging gfortran-13 (ZEN/EPYC), built on $(shell hostname) at $(shell date)"
 ACT = sed -e 's/^!\*qd/    /' # Enable quad-math statements
#ACT = sed -e 's/^!\*nq/    /' # Disable quad-math statements
# -Wconversion-extra produces a huge number of irrelevant warnings (and very occasionally relevant ones)
# -Wcompare-reals produces exclusively irrelewant warnings
# -Wunused-dummy-argument produces warnings we really don't care about
F90 = gfortran-13 -I. \
      -std=gnu -Wall -Wcharacter-truncation -Wextra -Wno-unused-dummy-argument -Wno-compare-reals \
         -Wimplicit-interface -Wimplicit-procedure -Wintrinsics-std \
         -Wrealloc-lhs-all \
         -fcheck=all \
      -O0 -fprotect-parens -march=native -mtune=native -fopenmp \
      -ffast-math -fcx-fortran-rules -freciprocal-math \
      -floop-block \
      -fno-realloc-lhs -fbacktrace -g \
      -static-libgfortran -static-libquadmath \
      -cpp -D__BUILD_ID__='$(BUILD_ID)' -ffree-line-length-none
#     -fexternal-blas -fblas-matmul-limit=50 # Causes memory errors
#     -floop-block # Requires Graphite isl, not configured on many distributions
#     -ffast-math -fcx-fortran-rules -mrecip 
F90L = $(F90) 
LAPACK = -lopenblas_openmp # -lquadlapack_gfortran
LIBEXTRA = -B/usr/share/libhugetlbfs -lhugetlbfs -Wl,--hugetlbfs-align
