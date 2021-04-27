exit 1;

Last updated: 2015 July 11
------------

Contents
========

0. Installation pre-requisites
1. Installation
2. Advanced installation
3. Further reading

Spherical-Coordinate Implicit Derivatives TDSE (SCID-TDSE)
=========================================================

This code implements a numerical solution of time-dependent Schroedinger equation 
(TDSE) for a one-electron atom subject to a laser field. TDSE is solved in spherical
coordinates, with the laser field treated in the velocity-gauge dipole approximation.

The code employs a norm-conserving propagation scheme of Harm Geert Muller
[H.G. Muller, "An efficient Propagation Scheme for the Time-Dependent,
Schroedinger Equation in the Velocity Gauge", Laser Physics 9, 138-148 (1999)],
adapted to support arbitrary radial grids and laser fields.

The code is parallelized for OpenMP, and can utilize 2-4 CPUs for runs using
linearly polarized light (axial symmetry). Many more CPUs can be used effectively 
for fully-3D runs; the actual efficiency depends on the size of the angular basis.

The code can be built in single, double, or quadruple precision (if supported by the
Fortran compiler used). Using other REAL kinds should also be possible, but may
require changes to the LAPACK interface (lapack.f90). 

All questions, comments, and requests should be sent to:

Serguei Patchkovskii, Serguei.Patchkovskii@mbi-berlin.de

0. Installation pre-requisites

Building SCID-TDSE requires a Fortran-95 compiler, and a working LAPACK library.
GNU Fortran 4.8.2 and Intel Fortran 15.0.1.133 are known to work. LAPACK 3.5.0 
is known to work. 

Building for parallel execution requires OpenMP (2.0 or better) support, and
LAPACK/BLAS built with re-entrancy support.

WARNING: Many vendor-supplied LAPACK produce incorrect results or hang 
WARNING: when used inside a parallel region.

On system with small default page size (e.g. x86_64), linking with libhugetlbfs
(http://libhugetlbfs.sourceforge.net/) is highly recommended.

1. Installation

By default, all system-specific configuration variables controlling SCID-TDSE 
build process are assembled in the file "vanilla.mak", found in the top-level
installation directory.

A build configuration file contains the following variables:

 BUILD_ID - Character string identifying a specific build; it will be printed
            each time the resulting binary is executed. See subroutine "start"
            in "spherical_tdse.f90" for more details.

 ACT      - Defines a shell command for processing quad-precision statements.
            The command must read a Fortran source file on the standard input,
            and produce an "activated" source on standard output. If quadruple-
            precision support is not needed, the command should uncomment all
            input lines starting with a string "!*nq". If quadruple-precision 
            support is to be enabled, the command should uncomment lines stating 
            with "!*qd". For example, on a Unix-like system:
              ACT = sed -e 's/^!\*qd/    /' # Enable quad-math statements
            or:
              ACT = sed -e 's/^!\*nq/    /' # Disable quad-math statements
            Activated sources are kept in the ./preprocess/ subdirectory.

 F90      - A command line used to invoke Fortran compiler, including all 
            optimization options. Among other things, this command line should
            enable Fortran source pre-processor, and define the variable
            __BUILD_ID__.
            For quadruple-precision builds, Intel Fortran is highly recommended.
            For single- and double-precision builds, GNU Fortran with OpenBLAS 
            is recommended.

 F90L     - A command used to invoke Fortran linker. Normally, this is the same
            command line used to invoke the compiler:
              F90L = $(F90)

 LAPACK   - Link instructions for LAPACK and BLAS libraries. If quad-precision
            support is enabled, this line should also include the special
            quad-precision LAPACK library (see "Prerequisites" above). If OpenMP
            parallelization is enabled, the libraries MUST be re-entrant;  they
            will be called from inside the parallel region.

 LIBEXTRA - Any other libraries needed to build the executable. For example,
            libhugetlbfs can appear here.

Please note that the default build configuration file emphasizes portability
over the performance and features. In particular, it disables quadruple-precision 
support, OpenMP support, as well as all system-specific optimizations. 

You will almost certainly get a faster code by using system-specific optimization 
options.

If desired, choose the desired real and integer kinds by editing "accuracy.f90".
The kind of the integer type used in the code is controlled by the "ik" parameter.
The real kind used in the calculations is controller by the "rk" parameter.
The default is to use (at least) a 4-byte integer and a double-precision
real. 

Finally, say "make" and wait for the dust to settle. If nothing goes wrong, make
will produce an executable file "spherical_tdse.x"

It is highly recommended to run at least some test cases, supplied in the "examples/"
subdirectory. The instructions for running the test cases and a bried description
of the tests are available in "doc/EXAMPLES.txt"

2. Advanced installation

It is possible to build SCID-TDSE to use quadruple-precision arithmetics on
systems which support it. Doing so requires a quadruple-precision versions of 
LAPACK and BLAS. Building these libraries from source is described in the file
"doc/README-quad.txt". Please ignore README-quad.txt unless you wish to build with
quadruple-precision support (which is not the default).

If you are feeling adventurous, a few examples of advanced build configuration 
files can be found in ./configs/ subdirectory, including:

 shelob-gfortran_opt.mak        - GNU Fortran, optimized
 shelob-gfortran_dbg.mak        - GNU Fortran, debug
 shelob-gfortran_noquad_opt.mak - GNU Fortran, optimized, no quad-precision libraries
 shelob-gfortran_noquad_dbg.mak - GNU Fortran, debug, no quad-precision libraries
 shelob-ifort_opt.mak           - Intel Fortran, optimized

All these examples assume a 64-bit Intel Linux system, with AVX instruction
support and libhugetlbfs installed.

Please note that the advanced build configurations will likely require further 
debugging.

2. Further reading

Additional documentation files are found in the "doc/" subdirectory.
