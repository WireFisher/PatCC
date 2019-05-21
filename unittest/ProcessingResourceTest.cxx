#include "mpi.h"
#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include "../processing_unit_mgt.h"
#include "../grid_decomposition.h"

extern Grid_info_manager *grid_info_mgr;
extern Process_thread_manager *process_thread_mgr;

class Mock_Process_thread_manager : public Process_thread_manager
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


int mpi_rank = -1;
int mpi_size = 0;

void get_default_hostname(char* hostname, int len)
{
    memcpy(hostname, "default", len);
};

void get_empty_hostname(char* hostname, int len)
{
    memcpy(hostname, "", len);
};

void get_different_hostname(char* hostname, int len)
{
    if(mpi_rank < 0)
        ASSERT_TRUE(false);
    else if(mpi_rank < mpi_size * 0.25)
        memcpy(hostname, "default01", len);
    else if(mpi_rank < mpi_size * 0.5)
        memcpy(hostname, "default02", len);
    else if(mpi_rank < mpi_size * 0.75)
        memcpy(hostname, "default03", len);
    else if(mpi_rank < mpi_size)
        memcpy(hostname, "default04", len);
    else
        ASSERT_TRUE(false);
};


void setup_dependency(void(*func)(char*, int), int num_threads)
{
    NiceMock<Mock_Process_thread_manager> *mock_process_thread_manager = new NiceMock<Mock_Process_thread_manager>;
    process_thread_mgr = mock_process_thread_manager;

    ON_CALL(*mock_process_thread_manager, get_hostname(_, _))
        .WillByDefault(Invoke(*func));

    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    ON_CALL(*mock_process_thread_manager, get_mpi_rank())
        .WillByDefault(Return(mpi_rank));

    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    ON_CALL(*mock_process_thread_manager, get_mpi_size())
        .WillByDefault(Return(mpi_size));

    ON_CALL(*mock_process_thread_manager, get_mpi_comm())
        .WillByDefault(Return(MPI_COMM_WORLD));

    ON_CALL(*mock_process_thread_manager, get_openmp_size())
        .WillByDefault(Return(num_threads));
};

void general_check(Processing_resource *processing_info, int total_num_threads)
{
    ASSERT_GT(processing_info->get_num_total_processing_units(), 0);
    ASSERT_EQ(processing_info->get_num_total_processing_units(), total_num_threads);
    ASSERT_NE(processing_info->get_processing_units(), (void*)NULL);
};

void check_proc_units(Processing_resource *processing_info, int num_hostnames, int max_num_threads_per_proc)
{
    unsigned int checksums[128];
    int num_checksums = 0;
    bool have_this_checksum;
    Processing_unit **processing_units;

    processing_units = processing_info->get_processing_units();

    for(int i = 0; i < processing_info->get_num_total_processing_units(); i ++) {
        have_this_checksum = false;
        for(int j = 0; j < num_checksums; j++)
            if(checksums[j] == processing_units[i]->hostname_checksum)
                have_this_checksum = true;
        if(!have_this_checksum)
            checksums[num_checksums++] = processing_units[i]->hostname_checksum;
        ASSERT_LE(num_checksums, num_hostnames);
        ASSERT_GE(processing_units[i]->process_id, 0);
        ASSERT_LT(processing_units[i]->process_id, mpi_size);
        ASSERT_GE(processing_units[i]->thread_id, 0);
        ASSERT_LT(processing_units[i]->thread_id, max_num_threads_per_proc);
        if(i != processing_info->get_num_total_processing_units() - 1)
            ASSERT_EQ(processing_units[i]->common_id, processing_units[i+1]->common_id - 1);
    }
    ASSERT_EQ(num_hostnames, num_checksums);
};

void check_active_units(Processing_resource *processing_info, int num_active_units, int total_num_threads, bool *active_flag)
{
    Processing_unit **processing_units;
    int num_total_processing_units;
    double ratio;
    
    processing_units = processing_info->get_processing_units();
    num_total_processing_units = processing_info->get_num_total_processing_units();

    ratio = (double)num_active_units / total_num_threads;
    ASSERT_GE(num_total_processing_units, 0);
    ASSERT_LE(num_total_processing_units, total_num_threads);
    ASSERT_NE(processing_units, (void*)NULL);

    unsigned int checksums[128];
    int num_active[128] = {0};
    int num_total[128] = {0};
    int current_node_index;
    int num_checksums = 0;
    bool have_this_checksum;

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

        if(active_flag[i])
            num_active[current_node_index] ++;
        num_total[current_node_index] ++;
    }

    int total_num_active = 0;
    for(int i = 0; i < num_checksums; i ++) {
        ASSERT_LE(num_active[i], num_total[i] * ratio + 1);
        //printf("%d: %d/%d\n", i, num_active[i], num_total[i]);
        ASSERT_GE(num_active[i], 0);
        total_num_active += num_active[i];
    }
    ASSERT_EQ(total_num_active, num_active_units);
    //printf("Total: %d/%d\n", total_num_active, total_num_threads);
};

TEST(ProcessingResourceTest, EmptyHostname) {
    int num_different_hostnames = 1;
    int max_num_threads_per_proc = 4;
    Processing_resource *proc_resrc;

    setup_dependency(get_empty_hostname, 4);
    proc_resrc = new Processing_resource();

    general_check(proc_resrc, mpi_size * 4);
    check_proc_units(proc_resrc, num_different_hostnames, max_num_threads_per_proc);

    delete proc_resrc;
    delete process_thread_mgr;
};

TEST(ProcessingResourceTest, DefaultHostname) {
    int num_different_hostnames = 1;
    int max_num_threads_per_proc = 4;
    Processing_resource *proc_resrc;

    setup_dependency(get_default_hostname, 4);
    proc_resrc = new Processing_resource();

    general_check(proc_resrc, mpi_size * 4);
    check_proc_units(proc_resrc, num_different_hostnames, max_num_threads_per_proc);

    delete proc_resrc;
    delete process_thread_mgr;
};

TEST(ProcessingResourceTest, DifferentHostname) {
    int num_different_hostnames = 4;
    int max_num_threads_per_proc = 4;
    Processing_resource *proc_resrc;

    setup_dependency(get_different_hostname, 4);
    proc_resrc = new Processing_resource();

    general_check(proc_resrc, mpi_size * 4);
    check_proc_units(proc_resrc, num_different_hostnames, max_num_threads_per_proc);

    delete proc_resrc;
    delete process_thread_mgr;
};

TEST(ProcessingResourceTest, ThreadNumOne) {
    int num_different_hostnames = 4;
    int num_thread = 1;
    Processing_resource *proc_resrc;

    setup_dependency(get_different_hostname, num_thread);
    proc_resrc = new Processing_resource();

    general_check(proc_resrc, mpi_size * num_thread);
    check_proc_units(proc_resrc, num_different_hostnames, num_thread);

    delete proc_resrc;
    delete process_thread_mgr;
};

TEST(ProcessingResourceTest, ThreadNumHundred) {
    int num_different_hostnames = 4;
    int num_thread = 131;
    Processing_resource *proc_resrc;

    setup_dependency(get_different_hostname, num_thread);
    proc_resrc = new Processing_resource();

    general_check(proc_resrc, mpi_size * num_thread);
    check_proc_units(proc_resrc, num_different_hostnames, num_thread);

    delete proc_resrc;
    delete process_thread_mgr;
};

TEST(ProcessingResourceTest, ThreadNumThousand) {
    int num_different_hostnames = 4;
    int num_thread = 1313;
    Processing_resource *proc_resrc;

    setup_dependency(get_different_hostname, num_thread);
    proc_resrc = new Processing_resource();

    general_check(proc_resrc, mpi_size * num_thread);
    check_proc_units(proc_resrc, num_different_hostnames, num_thread);

    delete proc_resrc;
    delete process_thread_mgr;
};

TEST(ProcessingResourceTest, ThreadNumMillion) {
    int num_different_hostnames = 4;
    int num_thread = 1313131;
    Processing_resource *proc_resrc;

    setup_dependency(get_different_hostname, num_thread);
    proc_resrc = new Processing_resource();

    general_check(proc_resrc, mpi_size * num_thread);
    check_proc_units(proc_resrc, num_different_hostnames, num_thread);

    delete proc_resrc;
    delete process_thread_mgr;
};

TEST(ProcessingResourceTest, ThreadNumDifferent1) {
    int num_different_hostnames = 4;
    int nums_thread[10] = {16, 32, 10, 40, 4, 11, 17, 1, 7, 19};
    int max_num_threads_per_proc = 40;
    int num_thread;
    int total_num_threads = 0;
    Processing_resource *proc_resrc;

    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    num_thread = nums_thread[mpi_rank%10];
    for(int i = 0; i < mpi_size; i++)
        total_num_threads += nums_thread[i%10];

    setup_dependency(get_different_hostname, num_thread);

    proc_resrc = new Processing_resource();

    general_check(proc_resrc, total_num_threads);
    check_proc_units(proc_resrc, num_different_hostnames, max_num_threads_per_proc);

    delete proc_resrc;
    delete process_thread_mgr;
};

TEST(ProcessingResourceTest, ThreadNumDifferent2) {
    int num_different_hostnames = 4;
    int nums_thread[10] = {9, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    int max_num_threads_per_proc = 9;
    int num_thread;
    int total_num_threads = 0;
    Processing_resource *proc_resrc;

    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    num_thread = nums_thread[mpi_rank%10];
    for(int i = 0; i < mpi_size; i++)
        total_num_threads += nums_thread[i%10];

    setup_dependency(get_different_hostname, num_thread);

    proc_resrc = new Processing_resource();

    general_check(proc_resrc, total_num_threads);
    check_proc_units(proc_resrc, num_different_hostnames, max_num_threads_per_proc);

    delete proc_resrc;
    delete process_thread_mgr;
};

TEST(ProcessingResourceTest, DISABLED_ActivedUnitsLT1) {
    int num_different_hostnames = 4;
    int nums_thread[10] = {1, 1, 4, 1, 1, 1, 1, 1, 1, 1};
    int max_num_threads_per_proc = 9;
    int num_thread;
    int total_num_threads = 0;
    Processing_resource *proc_resrc;
    bool *active_flag;

    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    num_thread = nums_thread[mpi_rank%10];
    for(int i = 0; i < mpi_size; i++)
        total_num_threads += nums_thread[i%10];

    setup_dependency(get_different_hostname, num_thread);

    proc_resrc = new Processing_resource();
    general_check(proc_resrc, total_num_threads);
    check_proc_units(proc_resrc, num_different_hostnames, max_num_threads_per_proc);

    active_flag = new bool[total_num_threads];
    proc_resrc->pick_out_active_processing_units(total_num_threads/2, active_flag);
    check_active_units(proc_resrc, total_num_threads/2, total_num_threads, active_flag);

    delete active_flag;
    delete proc_resrc;
    delete process_thread_mgr;
};

TEST(ProcessingResourceTest, DISABLED_ActivedUnitsLT2) {
    int num_different_hostnames = 4;
    int nums_thread[10] = {16, 32, 10, 40, 4, 11, 17, 1, 7, 19};
    int max_num_threads_per_proc = 40;
    int num_thread;
    int total_num_threads = 0;
    Processing_resource *proc_resrc;
    bool *active_flag;

    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    num_thread = nums_thread[mpi_rank%10];
    for(int i = 0; i < mpi_size; i++)
        total_num_threads += nums_thread[i%10];

    setup_dependency(get_different_hostname, num_thread);

    proc_resrc = new Processing_resource();
    general_check(proc_resrc, total_num_threads);
    check_proc_units(proc_resrc, num_different_hostnames, max_num_threads_per_proc);

    active_flag = new bool[total_num_threads];
    proc_resrc->pick_out_active_processing_units(total_num_threads/2, active_flag);
    check_active_units(proc_resrc, total_num_threads/2, total_num_threads, active_flag);

    delete active_flag;
    delete proc_resrc;
    delete process_thread_mgr;
};

TEST(ProcessingResourceTest, DISABLED_ActivedUnitsEQ) {
    int num_different_hostnames = 4;
    int nums_thread[10] = {16, 32, 10, 40, 4, 11, 17, 1, 7, 19};
    int max_num_threads_per_proc = 40;
    int num_thread;
    int total_num_threads = 0;
    Processing_resource *proc_resrc;
    bool *active_flag;

    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    num_thread = nums_thread[mpi_rank%10];
    for(int i = 0; i < mpi_size; i++)
        total_num_threads += nums_thread[i%10];

    setup_dependency(get_different_hostname, num_thread);

    proc_resrc = new Processing_resource();
    general_check(proc_resrc, total_num_threads);
    check_proc_units(proc_resrc, num_different_hostnames, max_num_threads_per_proc);

    active_flag = new bool[total_num_threads];
    proc_resrc->pick_out_active_processing_units(total_num_threads, active_flag);
    check_active_units(proc_resrc, total_num_threads, total_num_threads, active_flag);

    delete active_flag;
    delete proc_resrc;
    delete process_thread_mgr;
};
