#ifndef MPI_JACOBI_H_INCLUDE
#define MPI_JACOBI_H_INCLUDE

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include "serial.h"
#include "data.h"

#ifdef MPI
#include <mpi.h>

void get_counts(int *, int *, int);
void get_dimensions(int *, int *, int *, int);
#endif

void evolve(double *, double *, int); // Evolve Jacobi
void initialize(double *, double *, int); // Initialize matrices
void jacobi(double *, double *, int, int, double *);

#endif
