BUILD_ID :="Optimized ifort, built on $(shell hostname) at $(shell date)"
ACT = sed -e 's/^!\*qd/    /' # Enable quad-math statements
# This config file is tested with ifort 17.0.1.132 Build 20161005
F90 = ifort \
	    -I${MKLROOT}/include -heap-arrays 32 -warn -assume buffered_io \
            -qopenmp -qopenmp-simd -simd -qopt-report=5 \
            -xmic-avx512 -Ofast \
            -debug full -debug extended -traceback \
            -cpp -D__BUILD_ID__='$(BUILD_ID)'
F90L = $(F90) 
LAPACK = # -L/home/ps/lib64 -lquadlapack_intel
LIBEXTRA = -lhugetlbfs -Wl,-z,common-page-size=2097152 -Wl,-z,max-page-size=2097152 \
           -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl
#          -Wl,-rpath,$(MKLROOT)/lib/intel64/ -Wl,-rpath,/opt/intel/composer_xe_2015.1.133/compiler/lib/intel64/
