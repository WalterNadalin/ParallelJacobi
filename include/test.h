#ifndef SERIAL_JACOBI_H_INCLUDE
#define SERIAL_JACOBI_H_INCLUDE

#include <stdlib.h>
#include <stdio.h>
#include "data.h"

void serial_evolve(double *, double *, size_t); // Evolve Jacobi
void serial_initialize(double *, double *, size_t); // Initialize matrices
void serial_jacobi(double *, double *, size_t, size_t, double *);
size_t test(const char *, size_t, size_t);

#endif
