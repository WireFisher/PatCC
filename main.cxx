#include "component.h"
#include <cassert>
#include <cstdio>
#include <sys/time.h>
#include <unistd.h>
#include <sys/stat.h>

Grid_info_manager *grid_info_mgr;
Process_thread_manager *process_thread_mgr;

char usage[] = "usage: OMP_NUM_THREADS=nt mpiexec -n np ./patcc gridFile\n";
int main(int argc, char** argv)
{
    if (argc != 2) {
        perror(usage);
        return -1;
    }

    if (access("log", F_OK) != 0)
        mkdir("log", 0755);

    MPI_Init(&argc, &argv);

    process_thread_mgr = new Process_thread_manager();
    grid_info_mgr = new Grid_info_manager();

    if(!grid_info_mgr->read_grid_from_text(argv[1])) {
        perror("Failed in reading grid file\n");
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
