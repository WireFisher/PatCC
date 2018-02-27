#include "mpi.h"
#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include "../processing_unit_mgt.h"
#include "../delaunay_grid_decomposition_mgt.h"

extern Grid_info_manager *grid_info_mgr;
extern Process_thread_manager *process_thread_mgr;

class Mock_Process_thread_manager2 : public Process_thread_manager//_base
{
public:
    MOCK_METHOD2(get_hostname, void(char* hostname, int len));
    MOCK_METHOD0(get_openmp_size, int());
};

class Mock_Grid_info_manager : public Grid_info_manager
{
public:
    MOCK_METHOD1(get_grid_coord_values, double**(int));
    MOCK_METHOD1(get_grid_num_points, int(int));
    MOCK_METHOD5(get_grid_boundry, void(int, double*, double*, double*, double*));
    MOCK_METHOD2(get_polar_points, int(int, char));
    MOCK_METHOD1(is_grid_cyclic, bool(int));
};

using ::testing::Return;
using ::testing::NiceMock;
using ::testing::_;
using ::testing::Invoke;
using ::testing::ExitedWithCode;


static int mpi_rank = -1;
static int mpi_size = 0;

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

static void get_boundry(int grid_id, double* min_lon, double* max_lon, double* min_lat, double* max_lat)
{
    *min_lon =   0.0;
    *max_lon = 360.0;
    *min_lat = -90.0;
    *max_lat =  90.0;
}

TEST(GridDecompositionTest, GlobalManyDifferent) {
    int nums_thread[10] = {16, 32, 10, 1, 40, 11, 17, 19, 7, 4};
    int num_thread;
    int total_num_threads = 0;
    NiceMock<Mock_Process_thread_manager2> *mock_process_thread_manager = new NiceMock<Mock_Process_thread_manager2>;
    NiceMock<Mock_Grid_info_manager> *mock_grid_info_manager = new NiceMock<Mock_Grid_info_manager>;
    process_thread_mgr = mock_process_thread_manager;
    grid_info_mgr = mock_grid_info_manager;

    ON_CALL(*mock_process_thread_manager, get_hostname(_, _))
        .WillByDefault(Invoke(get_different_hostname));

    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    num_thread = nums_thread[mpi_rank%10];
    for(int i = 0; i < mpi_size; i++)
        total_num_threads += nums_thread[i%10];

    ON_CALL(*mock_process_thread_manager, get_openmp_size())
        .WillByDefault(Return(num_thread));

    double *coord_values[2];
    int size = 300;
    int num_points = size * size;
    coord_values[0] = new double[num_points]();
    coord_values[1] = new double[num_points]();
    for(int i = 0; i < size; i++)
        for(int j = 0; j < size; j++) {
            coord_values[0][i * size + j] =   0.0 + 360.0 * j / size;
            coord_values[1][i * size + j] = -90.0 + 180.0 * i / size;
        }
    ON_CALL(*mock_grid_info_manager, get_grid_coord_values(1))
        .WillByDefault(Return(coord_values));

    ON_CALL(*mock_grid_info_manager, get_grid_num_points(1))
        .WillByDefault(Return(num_points));

    ON_CALL(*mock_grid_info_manager, get_grid_boundry(1, _, _, _, _))
        .WillByDefault(Invoke(get_boundry));

    ON_CALL(*mock_grid_info_manager, get_polar_points(1, _))
        .WillByDefault(Return(3));

    ON_CALL(*mock_grid_info_manager, is_grid_cyclic(1))
        .WillByDefault(Return(true));

    Processing_resource *processing_info;
    Delaunay_grid_decomposition *grid_decomp;
    
    processing_info = new Processing_resource();
    grid_decomp = new Delaunay_grid_decomposition(1, processing_info, 100);
    grid_decomp->generate_grid_decomposition();

    FILE *log_file;
    char log_path[64];
    snprintf(log_path, 64, "log/GlobalManyDifferent_log.%d", mpi_rank);
    log_file = fopen(log_path, "w");
    fprintf(log_file, "size: %zu\n", grid_decomp->get_local_leaf_nodes().size());
    for(unsigned int i = 0; i < grid_decomp->get_local_leaf_nodes().size(); i++) {
        fprintf(log_file, "-%u: %d-\n",i ,grid_decomp->get_local_leaf_nodes()[i]->get_num_kernel_points());
        for(int j=0; j < grid_decomp->get_local_leaf_nodes()[i]->get_num_kernel_points(); j++)
            fprintf(log_file, "  (%lf, %lf)\n", grid_decomp->get_local_leaf_nodes()[i]->get_points_coord()[0][j], grid_decomp->get_local_leaf_nodes()[i]->get_points_coord()[1][j]);
    }
    fclose(log_file);

    delete grid_info_mgr;
    delete process_thread_mgr;
};
