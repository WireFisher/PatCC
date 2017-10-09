#ifndef P_DELAUNAY_TESTCASE_H
#define P_DELAUNAY_TESTCASE_H

#include "mpi.h"

void gethostname(char *, int);
int get_mpi_rank();
int get_mpi_size();
MPI_Comm get_mpi_comm();
int get_openmp_rank();
int get_openmp_size();

double** get_grid_coord_values(int);
int get_grid_num_points(int);
void get_grid_boundry(int, double*, double*, double*, double*);
int get_grid_size(int);
int get_polar_points(char);


#endif
