#include "mpi.h"
#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include "processing_unit_mgt.h"
#include "delaunay_grid_decomposition_mgt.h"

Grid_info_manager *grid_info_mgr;
Process_thread_manager_base *process_thread_mgr;

class Mock_Process_thread_manager : public Process_thread_manager_base
{
public:
    MOCK_METHOD2(get_hostname, void(char* hostname, int len));
    MOCK_METHOD0(get_mpi_rank, int());
    MOCK_METHOD0(get_mpi_size, int());
    MOCK_METHOD0(get_mpi_comm, MPI_Comm());
    MOCK_METHOD0(get_openmp_size, int());
};

using ::testing::Return;
using ::testing::NiceMock;
using ::testing::_;
using ::testing::Invoke;


void get_default_hostname(char* hostname, int len)
{
    memcpy(hostname, "default", len);
}

TEST(Test_processing_units_mgt, Basic) {
    NiceMock<Mock_Process_thread_manager> *mock_process_thread_manager = new NiceMock<Mock_Process_thread_manager>;
    process_thread_mgr = mock_process_thread_manager;
    grid_info_mgr = new Grid_info_manager;

    char *hostname;
    int len;
    ON_CALL(*mock_process_thread_manager, get_hostname(_, _))
        .WillByDefault(Invoke(get_default_hostname));

    int mpi_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    ON_CALL(*mock_process_thread_manager, get_mpi_rank())
        .WillByDefault(Return(mpi_rank));

    int mpi_size;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    ON_CALL(*mock_process_thread_manager, get_mpi_size())
        .WillByDefault(Return(mpi_size));

    ON_CALL(*mock_process_thread_manager, get_mpi_comm())
        .WillByDefault(Return(MPI_COMM_WORLD));

    ON_CALL(*mock_process_thread_manager, get_openmp_size())
        .WillByDefault(Return(1));

    Processing_info *processing_info;
    Processing_unit **processing_units;
    int num_total_processing_units;
    unsigned int hostname_checksum;
    int common_id_prev;

    processing_info = new Processing_info();
    //processing_info->pick_out_actived_processing_units(1,1, 200.0);
    processing_units = processing_info->get_processing_units();
    num_total_processing_units = processing_info->get_num_total_processing_units();

    ASSERT_GT(num_total_processing_units, 0);
    ASSERT_LE(num_total_processing_units, mpi_size);
    ASSERT_NE(processing_units, (void*)NULL);

    hostname_checksum = processing_units[0]->hostname_checksum;
    common_id_prev = 0;
    for(int i = 0; i < num_total_processing_units; i ++) {
        ASSERT_EQ(processing_units[i]->hostname_checksum, hostname_checksum);
        ASSERT_GE(processing_units[i]->process_id, 0);
        ASSERT_LT(processing_units[i]->process_id, mpi_size);
        ASSERT_EQ(processing_units[i]->thread_id, 0);
        if(i != num_total_processing_units - 1)
            ASSERT_EQ(processing_units[i]->common_id, processing_units[i+1]->common_id - 1);
    }

    delete grid_info_mgr;
    delete process_thread_mgr;
};

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
