#include "testcase.h"
#include <cstdio>
#include <cstring>
#include "processing_unit_mgt.h"

void gethostname(char *hostname, int len) {
    strncpy(hostname, "default", len);
    return;
}

int get_mpi_rank() {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    return rank;
}

int get_mpi_size() {
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    return size;
}

MPI_Comm get_mpi_comm() {
    return MPI_COMM_WORLD;
}

int get_openmp_rank() {
    return 0;

}

int get_openmp_size() {
    return 1;
}



double** get_grid_coord_values(int grid_id)
{
}
int get_grid_num_points(int grid_id)
{
}
void get_grid_boundry(int grid_id, double* min_lat, double* max_lat, double* min_lon, double* max_lon)
{
}
int get_grid_size(int grid_id)
{
}
int get_polar_points(char polar)
{
}


int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);
    Processing_info *processing_info = new Processing_info();
    processing_info->pick_out_actived_processing_units(1,2, 200.0);
    MPI_Finalize();
    return 0;
}
