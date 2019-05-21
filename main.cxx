#include "component.h"
#include <cassert>
#include <cstdio>
#include <sys/time.h>
#include <unistd.h>
#include <sys/stat.h>

Grid_info_manager *grid_info_mgr;
Process_thread_manager *process_thread_mgr;

char usage[] = "usage: OMP_NUM_THREADS=nt mpiexec -n np ./patcc gridFile\n";

void redirect_stdout()
{
    int rank;
    char log_path[128];

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0 && access("log", F_OK) != 0)
        mkdir("log", 0755);

    MPI_Barrier(MPI_COMM_WORLD);

    if(rank != 0) {
        snprintf(log_path, 128, "log/log.%d", rank);
        freopen(log_path, "w", stdout);
    }

#ifndef TIME_PERF
    setvbuf(stdout, 0, 2, 0);
#endif
}


int main(int argc, char** argv)
{
    if (argc != 2) {
        perror(usage);
        return -1;
    }

    MPI_Init(&argc, &argv);

    redirect_stdout();

    process_thread_mgr = new Process_thread_manager();
    grid_info_mgr = new Grid_info_manager();

    if(!grid_info_mgr->read_grid_from_text(argv[1])) {
        log(LOG_ERROR, "Failed in reading grid file\n");
        return -1;
    }

    Patcc* patcc = new Patcc(0);
    patcc->register_grid(new Grid(1));
    patcc->generate_delaunay_trianglulation(1);

    delete process_thread_mgr;
    delete grid_info_mgr;
    delete patcc;
    MPI_Finalize();
}
