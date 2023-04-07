#ifndef DATA_H_INCLUDE
#define DATA_H_INCLUDE

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>

#ifdef MPI
#include <mpi.h>
#endif

void save(double *, int, char *); // Save matrix to filei
double seconds(void); // Return the elapsed time
void distributed_save(double*, int, char *);
int get_parameters(int *, char **, int *, int *);
int test(char *, int, int);

#endif
