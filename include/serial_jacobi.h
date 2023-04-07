#ifndef SERIAL_JACOBI_H_INCLUDE
#define SERIAL_JACOBI_H_INCLUDE

#include <stdlib.h>
#include <stdio.h>

void serial_evolve(double *, double *, int); // Evolve Jacobi
void serial_initialize(double *, double *, int); // Initialize matrices
void serial_jacobi(double *, double *, int, int, double *);

#endif
