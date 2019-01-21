#include "component.h"
#include <cassert>
#include <cstdio>
#include <sys/time.h>

Grid_info_manager *grid_info_mgr;
Process_thread_manager *process_thread_mgr;

char usage[] = "usage:  mpiexec -n np ./patcc grid_file\n";
int main(int argc, char** argv)
{
    if (argc != 2) {
        perror(usage);
        return -1;
    }

    MPI_Init(&argc, &argv);

    process_thread_mgr = new Process_thread_manager();
    grid_info_mgr = new Grid_info_manager();

    if(!grid_info_mgr->read_grid_from_text(argv[1])) {
        perror("Failed in reading grid file\n");
        return -1;
    }

    Component* comp;
    comp = new Component(0);
    comp->register_grid(new Grid(1));
    comp->generate_delaunay_trianglulation(1);

    delete process_thread_mgr;
    delete grid_info_mgr;
    delete comp;
    MPI_Finalize();
}
