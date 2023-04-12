#ifndef DATA_H_INCLUDE
#define DATA_H_INCLUDE
#define COMMANDS 9

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <unistd.h>
#include <sys/time.h>
#include <string.h>

#ifdef MPI
#include <mpi.h>
#endif

void save(double *, size_t, char *); // Save matrix to file
void plot(double *, size_t, char *);
double seconds(void); // Return the elapsed time
size_t get_parameters(size_t *, char **, size_t *, size_t *);
void gather(double*, size_t, double *);
void get_counts(int *, int *, size_t);

#endif
