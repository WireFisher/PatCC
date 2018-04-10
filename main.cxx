#include "component.h"
#include <cassert>
#include <cstdio>
#include <sys/time.h>

Grid_info_manager *grid_info_mgr;
Process_thread_manager *process_thread_mgr;

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    timeval start, end;
    gettimeofday(&start, NULL);
    process_thread_mgr = new Process_thread_manager();
    grid_info_mgr = new Grid_info_manager();

    Component* comp;
    comp = new Component(0);
    comp->register_grid(new Grid(1));
    comp->generate_delaunay_trianglulation(1);

    delete process_thread_mgr;
    delete grid_info_mgr;
    delete comp;
    gettimeofday(&end, NULL);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if(rank == 0)
        printf("Total time: %ldms\n", ((end.tv_sec - start.tv_sec) * 1000000 + (end.tv_usec - start.tv_usec)) / 1000);
    MPI_Finalize();
}
