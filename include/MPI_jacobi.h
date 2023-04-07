#ifndef MPI_JACOBI_H_INCLUDE
#define MPI_JACOBI_H_INCLUDE

#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

void MPI_evolve(double *, double *, int); // Evolve Jacobi
void MPI_initialize(double *, double *, int); // Initialize matrices
void get_counts(int *, int *, int);
void get_dimensions(int *, int *, int *, int);
void MPI_jacobi(double *, double *, int, int, double *);

#endif
