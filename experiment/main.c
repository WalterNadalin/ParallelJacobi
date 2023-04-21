#include <complex.h>
#include <fftw3-mpi.h> // http://www.fftw.org/doc/MPI-Files-and-Data-Types.html#MPI-Files-and-Data-Types
#include <fftw3.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#define RED "\x1b[31m"
#define GREEN "\x1b[32m"
#define NORMAL "\x1b[m"

int index_f(int i1, int i2, int i3, int n2, int n3) {
  return n3 * (n2 * i1 + i2) + i3;
}

void print_3D(fftw_complex *mat, int nx, int ny, int nz) {
  for (int ix = 0; ix < nx; ix++) {
    for (int iy = 0; iy < ny; iy++) {
      for (int iz = 0; iz < nz; iz++) {
        printf("%f + i%f\t", creal(mat[index_f(ix, iy, iz, ny, nz)]),
               cimag(mat[index_f(ix, iy, iz, ny, nz)]));
      }
      printf("\n");
    }
    printf("\n");
  }
}

int equal(fftw_complex x, fftw_complex y, double eps) {
  if (creal(x) - creal(y) > eps || creal(x) - creal(y) < -eps ||
      cimag(x) - cimag(y) > eps || cimag(x) - cimag(y) < -eps)
    return 0;
  else
    return 1;
}

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);

  int size, rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  const int n1 = 8, n2 = 4, n3 = 6;
  const double L1 = 10., L2 = 10., L3 = 20., rad_conc = 0.6; 
  const int n1_local = n1 / size, local_size_grid = n1_local * n2 * n3,
            n2_local = n2 / size, local_size_slice = n2_local * n1_local * n3,
            global_size_grid = n1 * n2 * n3, n1_local_offset = rank * n1_local;
  int i1, i2, i3, index, index_buf;
  double f1conc, f2conc, f3conc, fac, ss = 0, ss_all, x1, x2, x3;

  if ((n1 % size || n2 % size) && !rank) {
    fprintf(stdout, "First two dimensions must be multiple of the number of "
                    "processes. The program will be aborted...\n\n");
    MPI_Abort(MPI_COMM_WORLD, 1);
    return 1;
  }

  fftw_complex *recv_buf =
      (fftw_complex *)fftw_malloc(local_size_grid * sizeof(fftw_complex));
  fftw_complex *send_buf =
      (fftw_complex *)fftw_malloc(local_size_grid * sizeof(fftw_complex));
  fftw_complex *fft_data =
      (fftw_complex *)fftw_malloc(local_size_grid * sizeof(fftw_complex));
  fftw_complex *check =
      (fftw_complex *)fftw_malloc(local_size_grid * sizeof(fftw_complex));

  MPI_Datatype blocks;
  MPI_Type_vector(n1_local, n2_local * n3, n2 * n3, MPI_C_DOUBLE_COMPLEX, &blocks);
  MPI_Type_commit(&blocks);

  for (i3 = 0; i3 < n3; ++i3) {
    x3 = L3 * ((double)i3) / n3;
    f3conc = exp(-pow((x3 - 0.5 * L3) / rad_conc, 2));

    for (i2 = 0; i2 < n2; ++i2) {
      x2 = L2 * ((double)i2) / n2;;
      f2conc = exp(-pow((x2 - 0.5 * L2) / rad_conc, 2));

      for (i1 = 0; i1 < n1_local; ++i1) {
        x1 = L1 * ((double)(i1 + n1_local_offset)) / n1;
        f1conc = exp(-pow((x1 - 0.5 * L1) / rad_conc, 2));
        index = index_f(i1, i2, i3, n2, n3);
        fft_data[index] = I * 0.0 + 1.f * rank + 10.f * i2 + 100.f * i3 + i1 * 1e-1; // f1conc * f2conc * f3conc + I * 0.0; 
        //index_buf = n3 * ((i2 % n2_local) * (1 - n1_local) + i2 * n1_local +
        //                  i1 * n2_local) +
        //            i3;
        //send_buf[index_buf] = fft_data[index];
        
        ss += fft_data[index];
      }
    }
  }

  // Normalize the concentration
  fac = L1 * L2 * L3 / global_size_grid;
  MPI_Allreduce(&ss, &ss_all, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  ss_all = 1.0 / (ss_all * fac);

  for (int i = 0; i < local_size_grid; ++i) {
   // fft_data[i] = creal(fft_data[i]) * ss_all + 0.0 * I;
  }
  
  int n_2d[] = {n2, n3};

 /* fftw_plan fw_plan_2d =
      fftw_plan_many_dft(2, n_2d, n1_local, fft_data, n_2d, 1, n_2d[0] * n_2d[1],
                         fft_data, n_2d, 1, n_2d[0] * n_2d[1], FFTW_FORWARD, FFTW_ESTIMATE);*/
 // fftw_execute(fw_plan_2d);
  
	for (i3 = 0; i3 < n3; ++i3) {
    for (i2 = 0; i2 < n2; ++i2) {
      for (i1 = 0; i1 < n1_local; ++i1) {
        index = index_f(i1, i2, i3, n2, n3);
        index_buf = n3 * ((i2 % n2_local) * (1 - n1_local) + i2 * n1_local +
                          i1 * n2_local) + i3;
       send_buf[index_buf] = fft_data[index];
      }
    }
  }

  // Reorder the different data blocks to be contigous in memory.
  // The new distribution will allow to use the Alltoall function
  int choose = atoi(argv[1]); // , try = atoi(argv[2]);

  if ((choose > size - 1)) {
    if (!rank)
      fprintf(stdout, "Process chosen not available, process set to 0.\n\n");
    choose = 0;
  }

  if (rank == choose)
    print_3D(fft_data, n1_local, n2, n3);

  // Perform an Alltoall communication
  int *counts_send = (int *)malloc(size * sizeof(int));
  int *displacements_send = (int *)malloc(size * sizeof(int));
  int *counts_recv = (int *)malloc(size * sizeof(int));
  int *displacements_recv = (int *)malloc(size * sizeof(int));
  MPI_Datatype *sendtypes = (MPI_Datatype *)malloc(size * sizeof(MPI_Datatype));
  MPI_Datatype *recvtypes = (MPI_Datatype *)malloc(size * sizeof(MPI_Datatype));

	displacements_send[0] = displacements_recv[0] = 0;
	sendtypes[0] = MPI_C_DOUBLE_COMPLEX;
	recvtypes[0] = MPI_C_DOUBLE_COMPLEX;
	
	for(int i = 0; i < size - 1; i++) {
		counts_send[i] = local_size_slice;
		counts_recv[i] = local_size_slice;
		displacements_send[i + 1] = 0; //displacements_send[i] + n2_local * n3 * sizeof(fftw_complex); 
		displacements_recv[i + 1] = 0; //displacements_recv[i] + local_size_slice * sizeof(fftw_complex);
		sendtypes[i] = MPI_C_DOUBLE_COMPLEX;
		recvtypes[i] = MPI_C_DOUBLE_COMPLEX;
	}
	
	counts_send[size - 1] = local_size_slice;
	counts_recv[size - 1] = local_size_slice;
	
 MPI_Alltoallw(fft_data, counts_send, displacements_send, sendtypes, recv_buf, counts_recv, displacements_recv, recvtypes, MPI_COMM_WORLD);
 
 /* MPI_Alltoall(fft_data, 1, blocks, recv_buf, local_size_slice,
               MPI_C_DOUBLE_COMPLEX, MPI_COMM_WORLD);*/
  MPI_Alltoall(send_buf, local_size_slice, MPI_C_DOUBLE_COMPLEX, check,
               local_size_slice, MPI_C_DOUBLE_COMPLEX, MPI_COMM_WORLD);

  // MPI_C_DOUBLE_COMPLEX
  if (rank == choose) {
     printf(" -----\n\n");
     print_3D(recv_buf, n1, n2_local, n3);
     printf(" .....\n\n");
   }
   
       MPI_Barrier(MPI_COMM_WORLD);
  //  Among i1 dimension
  double eps = 1e-6;
  int flag = 1;

    for (i3 = 0; i3 < n3; i3++) {
            if(!flag) break;
      for (i2 = 0; i2 < n2_local; i2++) {
                if(!flag) break;
        for (i1 = 0; i1 < n1; i1++) {
          index = (i1 * n3 + i3) * n2_local + i2;

          if (!equal(recv_buf[index], check[index], eps)) {
            //printf("%f + i%f\t\n", creal(recv_buf[index]),       cimag(recv_buf[index]));
           	//printf("%f + i%f\t\n", creal(check[index]),   cimag(check[index]));
            printf("%s %d: Results are NOT compatible\n\n", RED, rank);
            flag = 0;

            break;
          }
        }
      }
    }
    
    MPI_Barrier(MPI_COMM_WORLD);

    if (flag)
      printf("%s\n %d: Results are compatible\n\n", GREEN, rank);

    printf("%s", NORMAL);

  /*
     // Perform an Alltoall communication
                 // MPI_Alltoall(recv_buf, local_size_slice,
     MPI_C_DOUBLE_COMPLEX, send_buf, 1, blocks, fft->mpi_comm);

     //if(!rank) print_3D(send_buf, n1_local, n2, n3);
     // Reoder the different data blocks to be consistent with the initial
     // distribution.


     for (i1 = 0; i1 < n1_local; i1++) {
       for (i2 = 0; i2 < n2; i2++) {
         for (i3 = 0; i3 < n3; i3++) {
           index = index_f(i1, i2, i3, n1_local, n2, n3);
           data_rec[index] = send_buf[index];
         }
       }
     }*/
  MPI_Type_free(&blocks);
  fftw_free(recv_buf);
  fftw_free(send_buf);
  fftw_free(fft_data);
  fftw_free(check);
 // fftw_destroy_plan(fw_plan_2d);
  MPI_Finalize();
}
