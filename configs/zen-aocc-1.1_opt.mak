BUILD_ID :="Optimized AOCC 1.1 (ZEN/EPYC), built on $(shell hostname) at $(shell date)"
ACT = sed -e 's/^!\*qd/    /' # Enable quad-math statements
F90 = gfortran -I. -fopenmp \
     -cpp -D__BUILD_ID__='$(BUILD_ID)' -ffree-line-length-none \
     -O3 -mavx -madx -funroll-loops -ffast-math \
     -fcx-fortran-rules -mrecip -fno-realloc-lhs -fbacktrace -g \
     -fplugin=dragonegg.so -fplugin-arg-dragonegg-llvm-option="-merge-constant -lsr-in-nested-loop"
#
F90L = $(F90) -flto -Wl,-plugin-opt=-merge-constant \
                    -Wl,-plugin-opt=-lsr-in-nested-loop \
                    -Wl,-plugin-opt=-disable-vect-cmp
LAPACK = /home/ps/lib64/liblapack_gfortran.a -lopenblas
LIBEXTRA = -B/usr/share/libhugetlbfs -lhugetlbfs -Wl,--hugetlbfs-align \
           -ljemalloc -lgfortran /home/ps/AOCC/amdlibm-3.2.1/lib/static/libamdlibm.a -latomic
