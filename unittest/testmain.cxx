#include "mpi.h"
#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include "../processing_unit_mgt.h"
#include "../delaunay_grid_decomposition_mgt.h"

#include <omp.h>
#include <sched.h>
#include <sys/syscall.h>

Grid_info_manager *grid_info_mgr;
Process_thread_manager *process_thread_mgr;


int main(int argc, char** argv) {
    int result = 0;
    int rank;
    char *log_path;

    ::testing::InitGoogleMock(&argc, argv);
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    log_path = new char[32];
    snprintf(log_path, 32, "log/log.%d", rank);
    if(rank != 0)
        freopen(log_path, "w", stdout);
#ifndef TIME_PERF
    setvbuf(stdout, 0, 2, 0);
#endif

    result = RUN_ALL_TESTS();
    MPI_Finalize();

    return result;
}
