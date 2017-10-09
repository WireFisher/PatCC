#ifndef P_DELAUNAY_TESTCASE_H
#define P_DELAUNAY_TESTCASE_H

#include "mpi.h"

void gethostname(char *, int);
int get_mpi_rank();
int get_mpi_size();
MPI_Comm get_mpi_comm();
int get_openmp_rank();
int get_openmp_size();

#endif
