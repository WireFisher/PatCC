#include "mpi.h"
#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include "../processing_unit_mgt.h"
#include "../delaunay_grid_decomposition_mgt.h"


int main(int argc, char** argv) {
    int result = 0;
    int rank;

    ::testing::InitGoogleMock(&argc, argv);
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if(rank != 0)
        freopen("/dev/null", "w", stdout);
    result = RUN_ALL_TESTS();
    MPI_Finalize();

    return result;
}
