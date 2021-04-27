BUILD_ID :="Optimized ifort, built on $(shell hostname) at $(shell date)"
ACT = sed -e 's/^!\*qd/    /' # Enable quad-math statements
ACT2 = -e 's/^!\*mp/    /' # Enable MPI statements
#
#  -O3 -ipo tends to get very expensive ...
#  -axAVX2,MIC-AVX512 tends to get quite expensive, too!
#
F90 = mpiifort \
            -qopenmp -align all -align array256byte -pad \
            -warn -assume buffered_io \
            -O3 -ipo8 -no-prec-div -xAVX -axAVX2,MIC-AVX512 -complex_limited_range -fp-model fast=1 -ftz -assume protect_parens -heap-arrays 32 \
            -qopt-prefetch=0 -qopt-matmul -qopt-subscript-in-range -qopt-dynamic-align \
            -qopt-mem-layout-trans=3 -qopt-multi-version-aggressive \
            -debug full -debug extended -traceback -qopt-report=5 \
            -cpp -D__BUILD_ID__='$(BUILD_ID)'
#  Conservative: -complex_limited_range -fp-model fast=1 -ftz -assume protect_parens
#  Aggressive:   -complex_limited_range -fp-model fast=2 -ftz
# -fp-model precise -assume protect_parens # is necessary to make ifort fully compliant with the Fortran semantics
# -complex_limited_range -fp-model strict -ftz # produces code similar to gfortran in speed, but still less accurate
# -complex_limited_range -fp-model fast=1 -ftz # produces faster code, but dramatically decreases the accuracy
# -complex_limited_range -fp-model fast=2 -ftz # produces still slightly faster code, but dramatically decreases the accuracy
F90L = $(F90) 
LAPACK = # -L/home/ps/lib64 -lquadlapack_intel
# -static_mpi -static-intel  # -static_mpi is not compatible with using Fortran-08 bindings (use mpi_f08)
# -check_mpi is a debugging option. It requires Intel Trace Collector
LIBEXTRA = \
           -static_mpi -static-intel -qopenmp-link=static \
           -lhugetlbfs -Wl,-z,common-page-size=2097152 -Wl,-z,max-page-size=2097152 \
           -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_core.a $(MKLROOT)/lib/intel64/libmkl_sequential.a -Wl,--end-group -lpthread -lm 
#          -Wl,-rpath,$(MKLROOT)/lib/intel64/ -Wl,-rpath,/opt/intel/compilers_and_libraries_2018.0.128/linux/compiler/lib/intel64/
