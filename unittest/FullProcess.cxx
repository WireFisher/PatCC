#include "mpi.h"
#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include "../component.h"
#include "../processing_unit_mgt.h"
#include "../delaunay_grid_decomposition_mgt.h"

extern Grid_info_manager *grid_info_mgr;
extern Process_thread_manager *process_thread_mgr;

class Mock_Process_thread_manager3 : public Process_thread_manager
{
public:
    MOCK_METHOD2(get_hostname, void(char* hostname, int len));
    MOCK_METHOD0(get_openmp_size, int());
    MOCK_METHOD0(get_mpi_comm, MPI_Comm());
};

class Mock_Grid_info_manager2 : public Grid_info_manager
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
static double *coord_values[2] = {NULL, NULL};
static int num_points = 0;
static double min_lat, max_lat, min_lon, max_lon;


static void get_boundry(int grid_id, double* mi_lon, double* ma_lon, double* mi_lat, double* ma_lat)
{
    *mi_lon = min_lon;
    *ma_lon = max_lon;
    *mi_lat = min_lat;
    *ma_lat = max_lat;
}


static void prepare_grid()
{
    int size = 300;
    num_points = size * size;
    if (coord_values[0] || coord_values[1]) {
        delete coord_values[0];
        delete coord_values[1];
        coord_values[0] = coord_values[1] = NULL;
    }

    coord_values[0] = new double[num_points]();
    coord_values[1] = new double[num_points]();

    for(int i = 0; i < size; i++)
        for(int j = 0; j < size; j++) {
            coord_values[0][i * size + j] =   0.0 + 359.0 * j / size;
            coord_values[1][i * size + j] = -89.0 + 179.0 * i / size;
        }

    min_lon =   0.0;
    max_lon = 360.0;
    min_lat = -90.0;
    max_lat =  90.0;
}


void get_default_hostname2(char* hostname, int len)
{
    memcpy(hostname, "default", len);
};


TEST(FullProcess, Basic) {
    const int num_thread = 1;
    MPI_Comm comm = MPI_COMM_WORLD;

    NiceMock<Mock_Process_thread_manager3> *mock_process_thread_manager = new NiceMock<Mock_Process_thread_manager3>;
    NiceMock<Mock_Grid_info_manager2> *mock_grid_info_manager = new NiceMock<Mock_Grid_info_manager2>;
    process_thread_mgr = mock_process_thread_manager;
    grid_info_mgr = mock_grid_info_manager;

    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    ON_CALL(*mock_process_thread_manager, get_openmp_size())
        .WillByDefault(Return(num_thread));
    ON_CALL(*mock_process_thread_manager, get_mpi_comm())
        .WillByDefault(Return(comm));

    prepare_grid();

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

    ON_CALL(*mock_process_thread_manager, get_hostname(_, _))
        .WillByDefault(Invoke(get_default_hostname2));

    Component* comp;
    comp = new Component(0);
    comp->register_grid(new Grid(1));
    comp->generate_delaunay_trianglulation(1);

    delete comp;

    MPI_Comm split_world;
    if (mpi_size/3 > 1)
        MPI_Comm_split(MPI_COMM_WORLD, mpi_rank%3, mpi_rank/3, &split_world);

    if (mpi_size/3 > 1 && mpi_rank%3 == 0) {
        comm = split_world;
        ON_CALL(*mock_process_thread_manager, get_mpi_comm())
            .WillByDefault(Return(comm));

        comp = new Component(1);
        comp->register_grid(new Grid(1));
        comp->generate_delaunay_trianglulation(1);

        delete comp;

        int new_mpi_size;
        MPI_Comm_size(split_world, &new_mpi_size);

        FILE *fp;
        char fmt[] = "md5sum log/global_triangles_%d log/global_triangles_%d|awk -F\" \" '{print $1}'";
        char cmd[256];
        snprintf(cmd, 256, fmt, mpi_size, new_mpi_size);

        char md5[2][64];
        memset(md5[0], 0, 64);
        memset(md5[1], 0, 64);
        fp = popen(cmd, "r");
        fgets(md5[0], 64, fp);
        fgets(md5[1], 64, fp);
        ASSERT_EQ(strcmp(md5[0], md5[1]), 0);
    }
    delete process_thread_mgr;
    delete grid_info_mgr;
};
