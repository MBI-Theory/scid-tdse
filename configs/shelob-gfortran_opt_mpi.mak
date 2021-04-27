BUILD_ID :="Optimized gfortran with OpenMPI, built on $(shell hostname) at $(shell date)"
ACT = sed -e 's/^!\*qd/    /' # Enable quad-math statements
ACT2 = -e 's/^!\*mp/    /' # Enable MPI statements
F90 = /usr/lib64/mpi/gcc/openmpi3/bin/mpif90 -I. \
      -m64 -mavx -O3 -fprotect-parens -march=native -mtune=native -fopenmp \
      -ffast-math -fcx-fortran-rules -mrecip \
      -fexternal-blas -fblas-matmul-limit=50 \
      -fno-realloc-lhs -fbacktrace -g \
      -static-libgfortran \
      -cpp -D__BUILD_ID__='$(BUILD_ID)' -ffree-line-length-none
F90L = $(F90) 
LAPACK = /home/ps/lib64/liblapack_gfortran.a -lopenblas_openmp # -lquadlapack_gfortran
LIBEXTRA = -B/usr/share/libhugetlbfs -lhugetlbfs -Wl,--hugetlbfs-align -Wl,-rpath,/usr/lib64/mpi/gcc/openmpi3/lib64/
