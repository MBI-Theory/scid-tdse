.PHONY: goal clear

goal: makefile.dep
	  $(MAKE) spherical_tdse.x
	  # make build_pes.x

MAKEFLAGS = -r -j8

.SUFFIXES: .f90 .o .x .c .dep

#
#  This is the default; a config file may override it.
#
#ACT = sed -e 's/^!\*qd/    /' # Enable quad-math statements
ACT = sed -e 's/^!\*nq/    /' # Disable quad-math statements
#ACT2 = -e 's/^!\*mp/    /' # Enable MPI statements
ACT2 = -e 's/^!\*nm/    /' # Disable MPI statements

#
# System-specific overrides
#
  include vanilla.mak
# include configs/babel-gfortran_opt.mak
# include configs/babel-ifort18_opt.mak
# include configs/zen-gfortran-7_opt.mak
# include configs/zen-gfortran-11_opt.mak
# include configs/zen-gfortran-12_opt.mak
# include configs/zen-gfortran-13_opt.mak
# include configs/zen-gfortran-13_dbg.mak
# include configs/zen-oneapi_opt.mak
# include configs/zen-oneapi-i8_opt.mak
# include configs/zen-oneapi_opt_mpi.mak
# include configs/zen-aocc-1.1_opt.mak      # VERY SLOW CODE. DO NOT USE.
# include configs/oink-gfortran_opt.mak
# include configs/macos_m1-gfortran_opt.mak
# include configs/shelob-ifort_opt.mak
# include configs/shelob-ifort18_opt.mak
# include configs/shelob-ifort18_opt_mpi.mak
# include configs/shelob-ifort_noquad_opt.mak
# include configs/shelob-gfortran_opt.mak
# include configs/shelob-gfortran_opt_mpi.mak
# include configs/shelob-gfortran-7_opt.mak
# include configs/shelob-gfortran-10_opt.mak
# include configs/shelob-gfortran_dbg.mak
# include configs/shelob-gfortran_noquad_opt.mak
# include configs/shelob-gfortran_noquad_dbg.mak
# include configs/shelob-aocc20_opt.mak    # VERY SLOW & BUGGY CODE. DO NOT USE.
# include configs/vulcan_mic-ifort_opt.mak
# include configs/smeagol-ifort_opt.mak
# include configs/smeagol-gfortran_opt.mak
# include configs/sedna-ifort_opt.mak

#
# Finish the set-up
#
LIBS = $(LAPACK) $(LAPACK) $(LIBEXTRA)

#
# Compiling and archiving rules
#
.f90.o:
	$(ACT) $(ACT2) $< >preprocess/$<
	$(F90) -c preprocess/$<

#hacks.o:	hacks.f90 accuracy.o
#	$(F90) -O0 -c hacks.f90

clean:
	-/bin/rm -f *.{o,mod,x,il,a} *__genmod.f90 checkpoint_{field,main}.* makefile.dep *.optrpt *.dbg ./preprocess/*.f90

makefile.dep: $(shell echo *.f90)
	./make-depend.sh $^ > $@

#
# Explicit dependencies
#

LIBSPHERICAL += accuracy.o
LIBSPHERICAL += bicg_tools.o
LIBSPHERICAL += cap_tools.o
LIBSPHERICAL += checkpoint_tools.o
LIBSPHERICAL += composition_analysis.o
LIBSPHERICAL += constants.o
LIBSPHERICAL += coulomb_functions.o
LIBSPHERICAL += cubic_spline.o
LIBSPHERICAL += hacks.o
LIBSPHERICAL += lapack.o
LIBSPHERICAL += math.o
LIBSPHERICAL += node_tools.o
LIBSPHERICAL += potential_tools.o
LIBSPHERICAL += propagator_tools.o
LIBSPHERICAL += rotation_tools.o
LIBSPHERICAL += sort_tools.o
LIBSPHERICAL += spherical_data.o
LIBSPHERICAL += spherical_data_initialize.o
LIBSPHERICAL += spherical_bessel.o
LIBSPHERICAL += spherical_tdse_data.o
LIBSPHERICAL += spherical_tdse_field.o
LIBSPHERICAL += spherical_tdse_initialwf.o
LIBSPHERICAL += spherical_tdse_io.o
LIBSPHERICAL += spherical_tdse_propagate.o
LIBSPHERICAL += spherical_tsurf.o
LIBSPHERICAL += spherical_tsurf_data.o
LIBSPHERICAL += test_tools.o
LIBSPHERICAL += timer.o
LIBSPHERICAL += tridiagonal_cyclic.o
LIBSPHERICAL += tridiagonal_pivoted.o
LIBSPHERICAL += tridiagonal_tools.o
LIBSPHERICAL += vectorpotential_tools.o
LIBSPHERICAL += wavefunction_tools.o

#
# Building the binaries
#
spherical_tdse.x: spherical_tdse.o versions.o $(LIBSPHERICAL)
	$(F90L) -o spherical_tdse.x spherical_tdse.o versions.o $(LIBSPHERICAL) $(LIBS)
	-hugeedit --text --data spherical_tdse.x

build_pes.x: build_pes.o $(LIBSPHERICAL)
	$(F90L) -o build_pes.x build_pes.o $(LIBSPHERICAL) $(LIBS)
#
# Automatically-generated dependencies
#
include makefile.dep
