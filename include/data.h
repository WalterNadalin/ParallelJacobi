#ifndef DATA_H_INCLUDE
#define DATA_H_INCLUDE

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <unistd.h>
#include <sys/time.h>
#include <string.h>
#include "serial.h"
#include "simulation.h"

#ifdef MPI
#include <mpi.h>
#endif

void save(double *, int, char *); // Save matrix to file
void plot(double *, int, char *);
double seconds(void); // Return the elapsed time
int get_parameters(int *, char **, int *, int *);
void gather(double*, int, double *);

#endif
