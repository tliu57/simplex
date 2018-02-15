# simplex
simplex.c takes as input an MPS file and optimizes using the 2-phase standard simplex method.
simplex_mpi.c parallelizes the code in simplex.c using MPI.
simplex_hybrid.c parallelizes the code in simplex.c using both MPI and OpenMP.

To run the serial version, type ./runSerial.
To run the MPI version, type ./runMPI.
To run the hybrid version with both MPI and OpenMP, type ./runHybrid.
