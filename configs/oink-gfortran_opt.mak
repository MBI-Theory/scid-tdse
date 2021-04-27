BUILD_ID :="Optimized gfortran, built on $(shell hostname) at $(shell date)"
ACT = sed -e 's/^!\*qd/    /' # Enable quad-math statements
F90 = gfortran -I. -I${MKLROOT}/include \
      -m64 -mavx -O3 -march=native -mtune=native -fopenmp \
      -floop-block \
      -ffast-math -fcx-fortran-rules -mrecip \
      -fexternal-blas -fblas-matmul-limit=50 \
      -fno-realloc-lhs -fbacktrace -g \
      -cpp -D__BUILD_ID__='$(BUILD_ID)' -ffree-line-length-none
#     -ffast-math -fcx-fortran-rules -mrecip 
F90L = $(F90) 
#LAPACK = -L/home/ps/lib64 -llapack_gfortran -lblas 
LAPACK =  -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_gf_lp64.a $(MKLROOT)/lib/intel64/libmkl_core.a $(MKLROOT)/lib/intel64/libmkl_gnu_thread.a -Wl,--end-group -lpthread -lm -ldl
LIBEXTRA = -B/usr/share/libhugetlbfs -lhugetlbfs -Wl,--hugetlbfs-align
