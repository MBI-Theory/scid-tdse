BUILD_ID :="Optimized gfortran, built on $(shell hostname) at $(shell date)"
ACT = sed -e 's/^!\*qd/    /' # Enable quad-math statements
F90 = gfortran-6 -I. -I${MKLROOT}/include \
      -m64 -O3 -march=knl -mtune=knl -fopenmp \
      -ffast-math -fcx-fortran-rules -mrecip \
      -fexternal-blas -fblas-matmul-limit=50 \
      -fno-realloc-lhs -fbacktrace -g \
      -cpp -D__BUILD_ID__='$(BUILD_ID)' -ffree-line-length-none
F90L = $(F90) 
LAPACK =  -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_gf_lp64.a $(MKLROOT)/lib/intel64/libmkl_core.a $(MKLROOT)/lib/intel64/libmkl_gnu_thread.a -Wl,--end-group -lpthread -lm -ldl
LIBEXTRA = 
