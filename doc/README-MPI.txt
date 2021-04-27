Last updated: 2018 June 22
------------

1. General remarks

SCID-TDSE supports parallel execution on distributed-memory systems using
MPI communication library. In principle, all types of calculations and
input options are supported, with the following caveats:

 a. MPI parallelization should not be used within a single, shared-memory
    node. The OpenMP parallelization in SCID-TDSE is more efficient and
    less resource-intensive on shared-memory systems (including systems
    with moderately non-uniform memory access profiles).

 b. All nodes should have the same number of cores, should have the
    same architecture, and run at the same speed. SCID-TDSE partitions
    the worload using these assumptions, so that the faster nodes will
    end up waiting on the slowest node to complete its workshare.

 c. Each node should have sufficient RAM to complete the calculation
    on its own. All in-memory objects (wavefunctions, photoelectron
    spectra, etc.) are replicated on each node.

 d. The workload is distributed over the angular momentum L. A calculation
    must use a sufficiently large number of angular momenta to permit an
    equitable distribution of work between the nodes. Given (NN) nodes
    with (NC) CPU cores each, it is necessary to include at least (4*NN)
    angular momenta in the calculation. Optimal performance requires at 
    least (2*NC*NN) L values. Calculations using general field polarizations
    require more L values to be included to achieve optimal load 
    partitioning.

    Not all parts of the calculation are parallelized on distributed-memory
    configurations. Most importantly, propagation of Volkov phases (needed 
    for STS_VOLKOV=.TRUE.) are replicated on each node. In situations where
    many photoelectron amplitudes are needed, calculation of these phases
    may limit node scalability.

 e. All nodes must receive a copy of the input file on the standard input.
    This input processing happens _before_ the MPI library in initialized
    on the nodes. All nodes must have read access to all input files.

 f. The output files (spectra, tables, detailed output, visualization,
    checkpoints) are produced only on the master node, with the following
    exceptions:

    i) Each node will produce its own copy of the standard output. The
       name of this output can be controlled by (nt_node_output) input
       variable.
    ii) If caching of atomic field-free solutions is enabled
       (wt_atomic_cache_prefix/=' '), each node will update the cache
       locally.

 g. Checkpoint files can be used interchangeably between single-node
    and multi-node runs on the same architecture/compiler combination.
    MPI-enabled binary can be used for single-node calculations by
    setting nt_use_multinode=.false. in the input file.

 h. MPI communications and additional synchronization overhead have
    non-negligible cost. There is generally no guarantee of speedup 
    over a single-node run; please make sure to test on _your_ inputs.

 i. You need a fast, dedicated MPI network to run SCID-TDSE on multiple 
    nodes. You will almost certainly not get a speedup using standard
    ethernet networks and TCP/IP transports.

2. Notes for specific MPI implementations:

2.1 Intel MPI

The code was tested with Intel MPI 2018.0.128 and Intel Fortran 18.0.0.128
(see config file "shelob-ifort18_opt_mpi.mak").

This MPI implementation defaults to aggressive busy-waiting, both in MPI
and OpenMP liberaries. We obtain better performance by using the less
aggressive settings, e.g.:

    export HUGETLB_MORECORE=thp
    export OMP_STACKSIZE=500M
    export OMP_WAIT_POLICY=passive
    export I_MPI_WAIT_MODE=on
    ulimit -s 1024000
    mpirun -genvall -np N -s all spherical_tdse.x < input.inp > output.out

The "-s all" option instructs mpirun to present a copy of the standard
input on all nodes. It is also possible to encode the environment variables
directly in the mpirun command line, using the "-genv VAR=value" syntax.

Intel MPI will try hard to find a communication channel even if the fast (e.g.
the Infiniband) network fails to initialize. While this behaviour is very
user-friendly, it can also lead to mysterious slowdowns for applications which
do require a fast network. Please consider ether running with I_MPI_DEBUG=5
(which reports the low-level transport used by the library), or forcing the
library to abort if fast network is unavailable (e.g. by using -RDMA switch in
the mpirun command line).

2.2 OpenMPI

The code was tested with OpenMPI 1.10.6 and GNU Fortran 4.8.5 (see config file
"shelob-gfortran_opt_mpi.mak").

OpenMPI does not support redirection of standard input to all ranks simultaneosly.
It is therefore necessary to use a magic command-like switch "--scid-stdin":

   mpirun -mca btl ^tcp -mca mpi_show_mca_params all -bind-to none -np N spherical_tdse.x --scid-stdin input.inp > output.out

OpenMPI by default binds processes it starts to a specific socket. This is almost
certainly *not* the desired behaviour for an OpenMP process like scid-tdse. Using
"-bind-to none" will disable process/core binding altogether. The incantation
"-mca btl ^tcp" will stop OpenMPI from using the TCP transport (which is useless
for scid-tdse).

At least with OpenMPI 1.10.6, the library overhead is significantly worse than
with Intel MPI, so it is harder to achieve a speedup on the same hardware.
