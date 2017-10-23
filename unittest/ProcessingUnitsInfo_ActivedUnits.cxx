#include "mpi.h"
#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include "../processing_unit_mgt.h"
#include "../delaunay_grid_decomposition_mgt.h"

extern Grid_info_manager *grid_info_mgr;
extern Process_thread_manager *process_thread_mgr;

class Mock_Process_thread_manager : public Process_thread_manager//_base
{
public:
    MOCK_METHOD2(get_hostname, void(char* hostname, int len));
    MOCK_METHOD0(get_mpi_rank, int());
    MOCK_METHOD0(get_mpi_size, int());
    MOCK_METHOD0(get_mpi_comm, MPI_Comm());
    MOCK_METHOD0(get_openmp_size, int());
    //MOCK_METHOD7(allgather, int(void*, int, MPI_Datatype, void*, int, MPI_Datatype, MPI_Comm));
};

using ::testing::Return;
using ::testing::NiceMock;
using ::testing::_;
using ::testing::Invoke;
using ::testing::ExitedWithCode;


static int mpi_rank = -1;
static int mpi_size = 0;

static void get_default_hostname(char* hostname, int len)
{
    memcpy(hostname, "default", len);
}

static void get_empty_hostname(char* hostname, int len)
{
    memcpy(hostname, "", len);
}

static void get_different_hostname(char* hostname, int len)
{
    if(mpi_rank < 0)
        ASSERT_TRUE(false);
    else if(mpi_rank < mpi_size * 0.25)
        memcpy(hostname, "test_node_1", len);
    else if(mpi_rank < mpi_size * 0.5)
        memcpy(hostname, "test_node_2", len);
    else if(mpi_rank < mpi_size * 0.75)
        memcpy(hostname, "test_node_3", len);
    else if(mpi_rank < mpi_size)
        memcpy(hostname, "test_node_4", len);
    else
        ASSERT_TRUE(false);
}

/*
int allgather_num_threads(void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm)
{
    int* buf = (int*)recvbuf;
    for(int i = 0; i < mpi_size; i++)
        buf[i] = 4;
}

int allgather_hostname_checksum(void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm)
{
    unsigned int* buf = (unsigned int*)recvbuf;
    for(int i = 0; i < mpi_size; i++)
        buf[i] = 567;
}
    
    ON_CALL(*mock_process_thread_manager, allgather(_, _, MPI_INT, _, _, _, _))
        .WillByDefault(Invoke(allgather_num_threads));

    ON_CALL(*mock_process_thread_manager, allgather(_, _, MPI_UNSIGNED, _, _, _, _))
        .WillByDefault(Invoke(allgather_hostname_checksum));
    */

/*
TEST(ProcessingUnitsInfoDeathTest, ThreadNumZero) {
    NiceMock<Mock_Process_thread_manager> *mock_process_thread_manager = new NiceMock<Mock_Process_thread_manager>;
    process_thread_mgr = mock_process_thread_manager;
    grid_info_mgr = new Grid_info_manager;

    char *hostname;
    int len;
    ON_CALL(*mock_process_thread_manager, get_hostname(_, _))
        .WillByDefault(Invoke(get_different_hostname));

    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    ON_CALL(*mock_process_thread_manager, get_mpi_rank())
        .WillByDefault(Return(mpi_rank));

    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    ON_CALL(*mock_process_thread_manager, get_mpi_size())
        .WillByDefault(Return(mpi_size));

    ON_CALL(*mock_process_thread_manager, get_mpi_comm())
        .WillByDefault(Return(MPI_COMM_WORLD));

    ON_CALL(*mock_process_thread_manager, get_openmp_size())
        .WillByDefault(Return(0));

    Processing_resource *processing_info;
    Processing_unit **processing_units;
    int num_total_processing_units;
    unsigned int hostname_checksum;
    int common_id_prev;

    
    ASSERT_DEATH({processing_info = new Processing_resource(); }, "");
    //processing_info->pick_out_actived_processing_units(1,1, 200.0);

*
    ASSERT_GE(num_total_processing_units, 0);
    ASSERT_LE(num_total_processing_units, mpi_size * 4);
    ASSERT_NE(processing_units, (void*)NULL);

    unsigned int checksums[128];
    int num_checksums = 0;
    bool have_this_checksum;
    common_id_prev = 0;
    for(int i = 0; i < num_total_processing_units; i ++) {
        have_this_checksum = false;
        for(int j = 0; j < num_checksums; j++)
            if(checksums[j] == processing_units[i]->hostname_checksum)
                have_this_checksum = true;
        if(!have_this_checksum)
            checksums[num_checksums++] = processing_units[i]->hostname_checksum;
        ASSERT_LE(num_checksums, 5);
        ASSERT_GE(processing_units[i]->process_id, 0);
        ASSERT_LT(processing_units[i]->process_id, mpi_size);
        ASSERT_GE(processing_units[i]->thread_id, 0);
        ASSERT_LT(processing_units[i]->thread_id, 4);
        if(i != num_total_processing_units - 1)
            ASSERT_EQ(processing_units[i]->common_id, processing_units[i+1]->common_id - 1);
    }
*
    delete grid_info_mgr;
    delete process_thread_mgr;
};
*/

TEST(ProcessingUnitsInfoTest, ActivedUnitsLT) {
    int nums_thread[10] = {16, 32, 10, 1, 40, 11, 17, 19, 7, 4};
    int num_thread;
    int total_num_threads = 0;
    NiceMock<Mock_Process_thread_manager> *mock_process_thread_manager = new NiceMock<Mock_Process_thread_manager>;
    process_thread_mgr = mock_process_thread_manager;
    grid_info_mgr = new Grid_info_manager;

    char *hostname;
    int len;
    ON_CALL(*mock_process_thread_manager, get_hostname(_, _))
        .WillByDefault(Invoke(get_different_hostname));

    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    ON_CALL(*mock_process_thread_manager, get_mpi_rank())
        .WillByDefault(Return(mpi_rank));

    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    ON_CALL(*mock_process_thread_manager, get_mpi_size())
        .WillByDefault(Return(mpi_size));
    
    num_thread = nums_thread[mpi_rank%10];
    for(int i = 0; i < mpi_size; i++)
        total_num_threads += nums_thread[i%10];
    ON_CALL(*mock_process_thread_manager, get_mpi_comm())
        .WillByDefault(Return(MPI_COMM_WORLD));

    ON_CALL(*mock_process_thread_manager, get_openmp_size())
        .WillByDefault(Return(num_thread));

    Processing_resource *processing_info;
    Processing_unit **processing_units;
    int num_total_processing_units;
    unsigned int hostname_checksum;
    int common_id_prev;
    /*
    Workload_info *workload_info;
    
    processing_info = new Processing_resource();
    processing_info->pick_out_actived_processing_units(1, total_num_threads/2, 200.0);
    workload_info = processing_info->search_grid_workload_info(1);
    processing_units = processing_info->get_processing_units();
    num_total_processing_units = processing_info->get_num_total_processing_units();

    ASSERT_NE(workload_info, (void*)NULL);
    ASSERT_GE(num_total_processing_units, 0);
    ASSERT_LE(num_total_processing_units, total_num_threads);
    ASSERT_NE(processing_units, (void*)NULL);

    ASSERT_EQ(workload_info->grid_id, 1);
    ASSERT_EQ(workload_info->size, total_num_threads);
    ASSERT_NE(workload_info->is_actived, (void*)NULL);
    ASSERT_NE(workload_info->workloads, (void*)NULL);
    ASSERT_NE(workload_info->actived_common_id, (void*)NULL);

    unsigned int checksums[128];
    int num_actived[128] = {0};
    int num_total[128] = {0};
    int current_node_index;
    int num_checksums = 0;
    bool have_this_checksum;
    int current_actived_common_id_index = 0;
    common_id_prev = 0;
    for(int i = 0; i < num_total_processing_units; i ++) {
        have_this_checksum = false;
        for(int j = 0; j < num_checksums; j++)
            if(checksums[j] == processing_units[i]->hostname_checksum) {
                have_this_checksum = true;
                current_node_index = j;
                break;
            }
        if(!have_this_checksum) {
            checksums[num_checksums++] = processing_units[i]->hostname_checksum;
            current_node_index = num_checksums - 1;
        }

        if(workload_info->is_actived[i]) {
            num_actived[current_node_index] ++;
            ASSERT_EQ(i, workload_info->actived_common_id[current_actived_common_id_index++]);
        }
        num_total[current_node_index] ++;

        ASSERT_LE(num_checksums, 5);
        ASSERT_GE(processing_units[i]->process_id, 0);
        ASSERT_LT(processing_units[i]->process_id, mpi_size);
        ASSERT_GE(processing_units[i]->thread_id, 0);
        ASSERT_LT(processing_units[i]->thread_id, nums_thread[processing_units[i]->process_id%10]);
        if(i != num_total_processing_units - 1)
            ASSERT_EQ(processing_units[i]->common_id, processing_units[i+1]->common_id - 1);
    }

    int total_num_actived = 0;
    for(int i = 0; i < num_checksums; i ++) {
        ASSERT_LE(num_actived[i], num_total[i]/2 + 1);
        //printf("%d: %d/%d\n", i, num_actived[i], num_total[i]);
        ASSERT_GT(num_actived[i], 0);
        total_num_actived += num_actived[i];
    }
    ASSERT_EQ(workload_info->size_actived, total_num_actived);
    ASSERT_EQ(workload_info->size_actived, total_num_threads/2);
    //printf("Total: %d/%d\n", workload_info->size_actived, total_num_threads/2);
    */
    delete grid_info_mgr;
    delete process_thread_mgr;
};
