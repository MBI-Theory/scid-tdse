BUILD_ID :="Optimized Intel oneAPI(ifx), built on $(shell hostname) at $(shell date)"
ACT = sed -e 's/^!\*qd/    /' # Enable quad-math statements
#ACT2 = -e 's/^!\*mp/    /' # Enable MPI statements
# -axMIC-AVX512 is not supported by recent oneapi
F90 = ifx \
            -nogen-interfaces \
            -qopenmp -O3 -ipo -xAVX2 -qmkl=sequential -align all -falign-functions \
            -assume buffered_io -assume nobuffered_stdout -assume protect_parens -heap-arrays 32 \
            -warn -debug full -debug extended -traceback -qopt-report=5 \
            -cpp -D__BUILD_ID__='$(BUILD_ID)'
F90L = $(F90) 
LAPACK = # Lapack is in Intel MKL, 
LIBEXTRA = \
           -static-intel -qopenmp-link=static \
           -lhugetlbfs -Wl,-z,common-page-size=2097152 -Wl,-z,max-page-size=2097152 \
           -lpthread -lm 
#          -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_core.a $(MKLROOT)/lib/intel64/libmkl_sequential.a -Wl,--end-group -lpthread -lm 
