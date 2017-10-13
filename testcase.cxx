#include "processing_unit_mgt.h"
#include "delaunay_grid_decomposition_mgt.h"
#include <mpi.h>
#include <cstdio>
#include <cstring>
#include "omp.h"

Grid_info_manager *grid_info_mgr;
Process_thread_manager *process_thread_mgr;

int main(int argc, char* argv[]) {
    grid_info_mgr = new Grid_info_manager();
    process_thread_mgr = new Process_thread_manager();
    int mpi_rank;
    Processing_info *processing_info;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    //printf("omp_get_num_threads: %d\n", omp_get_num_threads());
    //printf("omp_get_max_threads: %d\n", omp_get_max_threads());
    //omp_set_num_threads(2);
    //printf("omp_get_max_threads: %d\n", omp_get_max_threads());
    #pragma omp parallel
    {
        #pragma omp master
        {
            processing_info = new Processing_info();

            processing_info->pick_out_actived_processing_units(1,2, 200.0);

            //if(mpi_rank == 0)
            //    processing_info->print_all_nodes_info();

            //decomposition
            
            //Delaunay
        }
    }
    MPI_Finalize();
    return 0;
}
