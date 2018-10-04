#include "mpi.h"
#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include "../component.h"
#include "../processing_unit_mgt.h"
#include "../delaunay_grid_decomposition_mgt.h"
#include "../netcdf_utils.h"
#include "../ccpl_utils.h"

extern Grid_info_manager *grid_info_mgr;
extern Process_thread_manager *process_thread_mgr;

double max_point_lat;
class Mock_Process_thread_manager3 : public Process_thread_manager
{
public:
    MOCK_METHOD0(get_mpi_comm, MPI_Comm());
};

class Mock_Grid_info_manager2 : public Grid_info_manager
{
public:
    MOCK_METHOD1(get_grid_coord_values, double**(int));
    MOCK_METHOD1(get_grid_num_points, int(int));
    MOCK_METHOD5(get_grid_boundry, void(int, double*, double*, double*, double*));
    MOCK_METHOD5(set_grid_boundry, void(int, double, double, double, double));
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
static bool is_cyclic = false;


class FullProcess : public ::testing::Test
{
public:
    virtual void SetUp() {
        mock_process_thread_manager = new NiceMock<Mock_Process_thread_manager3>;
        mock_grid_info_manager = new NiceMock<Mock_Grid_info_manager2>;

        process_thread_mgr = mock_process_thread_manager;
        grid_info_mgr = mock_grid_info_manager;
    }
    virtual void TearDown() {
        delete mock_process_thread_manager;
        delete mock_grid_info_manager;
        process_thread_mgr = NULL;
        grid_info_mgr = NULL;
    }

    NiceMock<Mock_Process_thread_manager3> *mock_process_thread_manager;
    NiceMock<Mock_Grid_info_manager2> *mock_grid_info_manager;
};


static void get_boundry(int grid_id, double* mi_lon, double* ma_lon, double* mi_lat, double* ma_lat)
{
    *mi_lon = min_lon;
    *ma_lon = max_lon;
    *mi_lat = min_lat;
    *ma_lat = max_lat;
}


static void set_boundry(int grid_id, double mi_lon, double ma_lon, double mi_lat, double ma_lat)
{
    min_lon = mi_lon;
    max_lon = ma_lon;
    min_lat = mi_lat;
    max_lat = ma_lat;
}


static void prepare_grid()
{
    int size = 300;
    num_points = size * size;
    delete coord_values[0];
    delete coord_values[1];
    coord_values[0] = coord_values[1] = NULL;

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


#ifdef NETCDF
void prepare_three_polar_grid()
{
    int num_dims;
    int *dim_size_ptr;
    int field_size;
    int field_size2;
    void *coord_buf0, *coord_buf1;
    bool squeeze = true;

    read_file_field_as_float("gridfile/three_polars_grid.nc", "nav_lon", &coord_buf0, &num_dims, &dim_size_ptr, &field_size);
    delete dim_size_ptr;
    read_file_field_as_float("gridfile/three_polars_grid.nc", "nav_lat", &coord_buf1, &num_dims, &dim_size_ptr, &field_size2);
    delete dim_size_ptr;
    assert(field_size == field_size2);
    num_points = field_size;
    coord_values[PDLN_LON] = (double*)coord_buf0;
    coord_values[PDLN_LAT] = (double*)coord_buf1;

    for(int i = 0; i < num_points; i++)
        if(coord_values[PDLN_LON][i] < 0.0)
            coord_values[PDLN_LON][i] += 360.0;

    for(int i = 0; i < num_points; i++)
        if(std::abs(coord_values[PDLN_LON][i] - 360.0) < PDLN_FLOAT_EQ_ERROR) {
            coord_values[PDLN_LON][i] = 0.0;
        }

    delete_redundent_points(coord_values[PDLN_LON], coord_values[PDLN_LAT], num_points);
    assert(have_redundent_points(coord_values[PDLN_LON], coord_values[PDLN_LAT], num_points) == false);

    if(squeeze) {
        for(int i = 0; i < num_points/100; i++) {
            coord_values[PDLN_LON][i] = coord_values[PDLN_LON][i*100];
            coord_values[PDLN_LAT][i] = coord_values[PDLN_LAT][i*100];
        }
        num_points /= 100;
    }

    min_lon =   0.0;
    max_lon = 360.0;
    min_lat = -80.0;
    max_lat =  90.0;
}

void prepare_latlon_grid()
{
    int num_dims;
    int *dim_size_ptr;
    int field_size;
    int field_size2;
    void *coord_buf0, *coord_buf1;

    read_file_field_as_double("gridfile/lonlat.nc", "lon", &coord_buf0, &num_dims, &dim_size_ptr, &field_size);
    delete dim_size_ptr;
    read_file_field_as_double("gridfile/lonlat.nc", "lat", &coord_buf1, &num_dims, &dim_size_ptr, &field_size2);
    delete dim_size_ptr;

    num_points = field_size*field_size2;
    coord_values[PDLN_LON] = new double [num_points];
    coord_values[PDLN_LAT] = new double [num_points];

    int count = 0;
    for(int i = 0; i < field_size; i ++)
        for(int j = 0; j < field_size2; j++) {
            coord_values[PDLN_LON][count] = ((double*)coord_buf0)[i];
            coord_values[PDLN_LAT][count++] = ((double*)coord_buf1)[j];
        }

    assert(count == num_points);
    assert(!have_redundent_points(coord_values[PDLN_LON], coord_values[PDLN_LAT], num_points));

    min_lon =   1.0;
    max_lon = 360.0;
    min_lat = -90.0;
    max_lat =  90.0;
}


void prepare_latlon_mutipolars()
{
    int num_dims;
    int *dim_size_ptr;
    int field_size;
    int field_size2;
    void *coord_buf0, *coord_buf1;

    read_file_field_as_double("gridfile/lonlat_90.nc", "lon", &coord_buf0, &num_dims, &dim_size_ptr, &field_size);
    delete dim_size_ptr;
    read_file_field_as_double("gridfile/lonlat_90.nc", "lat", &coord_buf1, &num_dims, &dim_size_ptr, &field_size2);
    delete dim_size_ptr;

    num_points = field_size*field_size2;
    coord_values[PDLN_LON] = new double [num_points];
    coord_values[PDLN_LAT] = new double [num_points];

    int count = 0;
    for(int j = field_size2-1; j >= 0; j--)
        for(int i = 0; i < field_size; i ++) {
            coord_values[PDLN_LON][count] = ((double*)coord_buf0)[i];
            coord_values[PDLN_LAT][count++] = ((double*)coord_buf1)[j];
        }

    assert(count == num_points);
    assert(!have_redundent_points(coord_values[PDLN_LON], coord_values[PDLN_LAT], num_points));

    min_lon =   0.0;
    max_lon = 360.0;
    min_lat = -90.0;
    max_lat =  90.0;
}


void prepare_latlon_singlepolar()
{
    int num_dims;
    int *dim_size_ptr;
    int field_size;
    int field_size2;
    void *coord_buf0, *coord_buf1;

    read_file_field_as_double("gridfile/lonlat_90.nc", "lon", &coord_buf0, &num_dims, &dim_size_ptr, &field_size);
    delete dim_size_ptr;
    read_file_field_as_double("gridfile/lonlat_90.nc", "lat", &coord_buf1, &num_dims, &dim_size_ptr, &field_size2);
    delete dim_size_ptr;

    num_points = field_size*field_size2;
    coord_values[PDLN_LON] = new double [num_points+2];
    coord_values[PDLN_LAT] = new double [num_points+2];

    int count = 0;
    for(int j = field_size2-1; j >= 0; j--)
        for(int i = 0; i < field_size; i ++) {
            if(fabs(((double*)coord_buf1)[j] - 90) > 1e-8 && fabs(((double*)coord_buf1)[j] - -90) > 1e-8) {
                coord_values[PDLN_LON][count] = ((double*)coord_buf0)[i];
                coord_values[PDLN_LAT][count++] = ((double*)coord_buf1)[j];
            }
        }

    coord_values[PDLN_LON][count] = 0.0;
    coord_values[PDLN_LAT][count++] = 90.0;
    coord_values[PDLN_LON][count] = 0.0;
    coord_values[PDLN_LAT][count++] = -90.0;
    num_points = count;
    assert(!have_redundent_points(coord_values[PDLN_LON], coord_values[PDLN_LAT], num_points));

    min_lon =   0.0;
    max_lon = 360.0;
    min_lat = -90.0;
    max_lat =  90.0;
}
#endif


TEST_F(FullProcess, Basic) {
    MPI_Comm comm = MPI_COMM_WORLD;


    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    ON_CALL(*mock_process_thread_manager, get_mpi_comm())
        .WillByDefault(Return(comm));

    prepare_grid();

    ON_CALL(*mock_grid_info_manager, get_grid_coord_values(1))
        .WillByDefault(Return(coord_values));

    ON_CALL(*mock_grid_info_manager, get_grid_num_points(1))
        .WillByDefault(Return(num_points));

    ON_CALL(*mock_grid_info_manager, get_grid_boundry(1, _, _, _, _))
        .WillByDefault(Invoke(get_boundry));

    ON_CALL(*mock_grid_info_manager, is_grid_cyclic(1))
        .WillByDefault(Return(true));

    Component* comp;
    comp = new Component(0);
    comp->register_grid(new Grid(1));
    comp->generate_delaunay_trianglulation(1, true);

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
        comp->generate_delaunay_trianglulation(1, true);

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
        ASSERT_STREQ(md5[0], md5[1]);
    }
};


#ifdef NETCDF
TEST_F(FullProcess, LatLonGrid) {
    MPI_Comm comm = MPI_COMM_WORLD;

    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    ON_CALL(*mock_process_thread_manager, get_mpi_comm())
        .WillByDefault(Return(comm));

    prepare_latlon_grid();

    ON_CALL(*mock_grid_info_manager, get_grid_coord_values(1))
        .WillByDefault(Return(coord_values));

    ON_CALL(*mock_grid_info_manager, get_grid_num_points(1))
        .WillByDefault(Return(num_points));

    ON_CALL(*mock_grid_info_manager, get_grid_boundry(1, _, _, _, _))
        .WillByDefault(Invoke(get_boundry));

    ON_CALL(*mock_grid_info_manager, is_grid_cyclic(1))
        .WillByDefault(Return(true));

    Component* comp;
    comp = new Component(0);
    comp->register_grid(new Grid(1));
    comp->generate_delaunay_trianglulation(1, true);

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
        comp->generate_delaunay_trianglulation(1, true);

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
        ASSERT_STREQ(md5[0], md5[1]);
    }
};


TEST_F(FullProcess, LatLonSinglePolar) {
    MPI_Comm comm = MPI_COMM_WORLD;

    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    ON_CALL(*mock_process_thread_manager, get_mpi_comm())
        .WillByDefault(Return(comm));

    prepare_latlon_singlepolar();

    ON_CALL(*mock_grid_info_manager, get_grid_coord_values(1))
        .WillByDefault(Return(coord_values));

    ON_CALL(*mock_grid_info_manager, get_grid_num_points(1))
        .WillByDefault(Return(num_points));

    ON_CALL(*mock_grid_info_manager, get_grid_boundry(1, _, _, _, _))
        .WillByDefault(Invoke(get_boundry));

    ON_CALL(*mock_grid_info_manager, is_grid_cyclic(1))
        .WillByDefault(Return(true));

    Component* comp;
    comp = new Component(0);
    comp->register_grid(new Grid(1));
    comp->generate_delaunay_trianglulation(1, true);

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
        comp->generate_delaunay_trianglulation(1, true);

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
        ASSERT_STREQ(md5[0], md5[1]);
    }
};



TEST_F(FullProcess, LatLonMutiPolars) {
    MPI_Comm comm = MPI_COMM_WORLD;

    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    ON_CALL(*mock_process_thread_manager, get_mpi_comm())
        .WillByDefault(Return(comm));

    prepare_latlon_mutipolars();

    ON_CALL(*mock_grid_info_manager, get_grid_coord_values(1))
        .WillByDefault(Return(coord_values));

    ON_CALL(*mock_grid_info_manager, get_grid_num_points(1))
        .WillByDefault(Return(num_points));

    ON_CALL(*mock_grid_info_manager, get_grid_boundry(1, _, _, _, _))
        .WillByDefault(Invoke(get_boundry));

    ON_CALL(*mock_grid_info_manager, is_grid_cyclic(1))
        .WillByDefault(Return(true));

    Component* comp;
    comp = new Component(0);
    comp->register_grid(new Grid(1));
    comp->generate_delaunay_trianglulation(1, true);

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
        comp->generate_delaunay_trianglulation(1, true);

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
        ASSERT_STREQ(md5[0], md5[1]);
    }
};


TEST_F(FullProcess, ThreePolar) {
    MPI_Comm comm = MPI_COMM_WORLD;

    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    ON_CALL(*mock_process_thread_manager, get_mpi_comm())
        .WillByDefault(Return(comm));

    prepare_three_polar_grid();

    ON_CALL(*mock_grid_info_manager, get_grid_coord_values(1))
        .WillByDefault(Return(coord_values));

    ON_CALL(*mock_grid_info_manager, get_grid_num_points(1))
        .WillByDefault(Return(num_points));

    ON_CALL(*mock_grid_info_manager, get_grid_boundry(1, _, _, _, _))
        .WillByDefault(Invoke(get_boundry));

    ON_CALL(*mock_grid_info_manager, is_grid_cyclic(1))
        .WillByDefault(Return(true));

    Component* comp;
    comp = new Component(0);
    comp->register_grid(new Grid(1));
    comp->generate_delaunay_trianglulation(1, true);

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
        comp->generate_delaunay_trianglulation(1, true);

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
        ASSERT_STREQ(md5[0], md5[1]);
    }
};


#define CHECK_PARALLEL_CONSISTENCY (true)
const char dim1_grid_path[] = "gridfile/many_types_of_grid/one_dimension/%s";
//"V3_Greenland_pole_x1_T_grid.nc", //md5sum wrong
const char dim1_grid_name[][64] = {
    /*
    */
    "ne30np4-t2.nc",  //assert false | 360point: not assert but false | got wrong fake cyclic triangles: OK | md5sum wrong | 15 pass 5 assert false|md5sum wrong: OK by 1e-10
    "ne60np4_pentagons_100408.nc", //x assert false | 360point: OK | md5sum wrong |15 pass 5 assert false|md5sum wrong: OK by 1e-10
    "gx3v5_Present_DP_x3.nc", //x | md5sum wrong: OK
    "Gamil_360x180_Grid.nc", // x ok
    "licom_eq1x1_degree_Grid.nc", //x ok
    "licom_gr1x1_degree_Grid.nc", //x ok
    "LICOM_P5_Grid.nc",
    "ll1deg_grid.nc",
    "ll2.5deg_grid.nc",
    "R42_Gaussian_Grid.nc",
    "T85_Gaussian_Grid.nc",
    "Gamil_128x60_Grid.nc", // x | deleting outter triangle: ok
    "fv1.9x2.5_050503.nc", // x ok
    "T42_Gaussian_Grid.nc",
    "T42_grid.nc",
    "T62_Gaussian_Grid.nc",
    "ar9v4_100920.nc", // x can't pass check cause extreme triangles: introducing threshold OK | md5sum wrong: virtual p
    "wr50a_090301.nc", //assert length false: wrong support for non-0~360 grid| can't pass check cause extreme triangles|15 pass 5 max|grid preteatment: OK
    "Version_3_of_Greenland_pole_x1_T-grid.nc", //x | md5sum wrong
    "R05_Grid.nc",
    /*
    */
};
const char dim1_global_grid_name[][64] = {
    "Gamil_128x60_Grid.nc",
    "fv1.9x2.5_050503.nc",
    "Gamil_360x180_Grid.nc",
    "R05_Grid.nc",
    "ll1deg_grid.nc",
    "ll2.5deg_grid.nc",
    "R42_Gaussian_Grid.nc",
    "T42_Gaussian_Grid.nc",
    "T42_grid.nc",
    "T62_Gaussian_Grid.nc",
    "T85_Gaussian_Grid.nc",
};


void prepare_dim1_grid(const char grid_name[])
{
    char fullname[128];
    int num_dims;
    int *dim_size_ptr;
    int field_size;
    int field_size2;
    void *coord_buf0, *coord_buf1;
    char lon_unit[32];
    char lat_unit[32];
    int squeeze_ratio = 0;

    if (mpi_rank == 0) {
        if(strncmp(grid_name, "ar9v4_100920.nc", 15) == 0)
            squeeze_ratio = 100;
        if(strncmp(grid_name, "Version_3_of_Greenland_pole_x1_T-grid.nc", 64) == 0)
            squeeze_ratio = 10;

        snprintf(fullname, 128, dim1_grid_path, grid_name);
        read_file_field_as_double(fullname, "grid_center_lon", &coord_buf0, &num_dims, &dim_size_ptr, &field_size, lon_unit);
        delete dim_size_ptr;
        read_file_field_as_double(fullname, "grid_center_lat", &coord_buf1, &num_dims, &dim_size_ptr, &field_size2, lat_unit);
        delete dim_size_ptr;

        ASSERT_EQ(field_size, field_size2);
        num_points = field_size;
        coord_values[PDLN_LON] = (double*)coord_buf0;
        coord_values[PDLN_LAT] = (double*)coord_buf1;

        min_lon = 1e10;
        max_lon = -1e10;
        min_lat = 1e10;
        max_lat = -1e10;

        bool convert = strncmp(lon_unit, "radians", 32) == 0;
        for(int i = 0; i < num_points; i ++) {
            if(convert) {
                coord_values[PDLN_LON][i] = RADIAN_TO_DEGREE(coord_values[PDLN_LON][i]);
                coord_values[PDLN_LAT][i] = RADIAN_TO_DEGREE(coord_values[PDLN_LAT][i]);
            }
            while(coord_values[PDLN_LON][i] >= 360)
                coord_values[PDLN_LON][i] -= 360;
            if(coord_values[PDLN_LON][i] < min_lon) min_lon = coord_values[PDLN_LON][i];
            if(coord_values[PDLN_LON][i] > max_lon) max_lon = coord_values[PDLN_LON][i];
            if(coord_values[PDLN_LAT][i] < min_lat) min_lat = coord_values[PDLN_LAT][i];
            if(coord_values[PDLN_LAT][i] > max_lat) max_lat = coord_values[PDLN_LAT][i];
        }
     
        if(squeeze_ratio > 0) {
            for(int i = 0; i < num_points/squeeze_ratio; i++) {
                coord_values[PDLN_LON][i] = coord_values[PDLN_LON][i*squeeze_ratio];
                coord_values[PDLN_LAT][i] = coord_values[PDLN_LAT][i*squeeze_ratio];
            }
            num_points = num_points/squeeze_ratio;
        }

        //printf("num points: %d\n", num_points);
        //printf("point range: %lf, %.60lf, %lf, %lf\n", min_lon, max_lon, min_lat, max_lat);
        max_point_lat = max_lat;
        max_lon += 0.0001;
        if(max_lon > 360) max_lon = 360;
        max_lat += 0.0001;
        if(max_lat > 90) max_lat = 90;
        assert(!have_redundent_points(coord_values[PDLN_LON], coord_values[PDLN_LAT], num_points));

        for(unsigned i = 0; i < sizeof(dim1_global_grid_name)/64; i++)
            if(strncmp(grid_name, dim1_global_grid_name[i], 64) == 0) {
                min_lon = 0;
                max_lon = 360;
                min_lat = -90;
                max_lat = 90;
                break;
            }

        if(strncmp(grid_name, "V3_Greenland_pole_x1_T_grid.nc", 64) == 0 ||
           strncmp(grid_name, "Version_3_of_Greenland_pole_x1_T-grid.nc", 64) == 0 ||
           strncmp(grid_name, "ar9v4_100920.nc", 64) == 0 ||
           strncmp(grid_name, "gx3v5_Present_DP_x3.nc", 64) == 0) {
                min_lon = 0;
                max_lon = 360;
                max_lat = 90;
        }
        if(strncmp(grid_name, "licom_eq1x1_degree_Grid.nc", 64) == 0 ||
           strncmp(grid_name, "licom_gr1x1_degree_Grid.nc", 64) == 0 ||
           strncmp(grid_name, "LICOM_P5_Grid.nc", 64) == 0) {
                min_lon = 0;
                max_lon = 360;
                max_lat = 90;
        }
        if(strncmp(grid_name, "wr50a_090301.nc", 64) == 0) {
            min_lon = -180;
            max_lon = 180;
            max_lat = 90;
        }

        if(fabs((max_lon - min_lon) - 360) < 0.5)
            is_cyclic = true;
        else
            is_cyclic = false;
    }

    MPI_Bcast(&num_points, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (mpi_rank != 0) {
        coord_values[PDLN_LON] = new double[num_points];
        coord_values[PDLN_LAT] = new double[num_points];
    }
    MPI_Bcast(coord_values[PDLN_LON], num_points, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(coord_values[PDLN_LAT], num_points, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&min_lon, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&max_lon, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&min_lat, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&max_lat, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&is_cyclic, 1, MPI_CHAR, 0, MPI_COMM_WORLD);
    assert(sizeof(bool) == 1);
#ifdef TIME_PERF
    printf("[ - ] Total points: %d\n", num_points);
#endif
};


#include <unistd.h>
TEST_F(FullProcess, ManyTypesOfGrids) {
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Comm comm = MPI_COMM_WORLD;

    printf("pid: %d\n", getpid());
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    MPI_Comm split_world;
    if (mpi_size/3 > 1)
        MPI_Comm_split(MPI_COMM_WORLD, mpi_rank%3, mpi_rank/3, &split_world);

    for(unsigned i = 0; i < sizeof(dim1_grid_name)/64; i++) {
        MPI_Barrier(MPI_COMM_WORLD);
        printf("processing: %s\n", dim1_grid_name[i]);
        prepare_dim1_grid(dim1_grid_name[i]);

        ON_CALL(*mock_process_thread_manager, get_mpi_comm())
            .WillByDefault(Return(comm));

        ON_CALL(*mock_grid_info_manager, get_grid_coord_values(1))
            .WillByDefault(Return(coord_values));

        ON_CALL(*mock_grid_info_manager, get_grid_num_points(1))
            .WillByDefault(Return(num_points));

        ON_CALL(*mock_grid_info_manager, get_grid_boundry(1, _, _, _, _))
            .WillByDefault(Invoke(get_boundry));

        ON_CALL(*mock_grid_info_manager, set_grid_boundry(1, _, _, _, _))
            .WillByDefault(Invoke(set_boundry));

        ON_CALL(*mock_grid_info_manager, is_grid_cyclic(1))
            .WillByDefault(Return(is_cyclic));

        Component* comp;
        comp = new Component(0);
        comp->register_grid(new Grid(1));
        int ret = comp->generate_delaunay_trianglulation(1, true);
        EXPECT_EQ(ret, 0);
        delete comp;

        if (CHECK_PARALLEL_CONSISTENCY && mpi_size/3 > 1 && mpi_rank%3 == 0) {
            printf("spliting world: %s\n", dim1_grid_name[i]);
            MPI_Barrier(split_world);
            ON_CALL(*mock_process_thread_manager, get_mpi_comm())
                .WillByDefault(Return(split_world));

            comp = new Component(1);
            comp->register_grid(new Grid(1));
            int ret = comp->generate_delaunay_trianglulation(1, true);
            EXPECT_EQ(ret, 0);

            delete comp;

            if (mpi_rank == 0) {
                int new_mpi_size;
                MPI_Comm_size(split_world, &new_mpi_size);

                FILE *fp;
                char fmt[] = "md5sum log/global_triangles_* | awk -F\" \" '{print $1}'";
                char cmd[256];
                snprintf(cmd, 256, fmt, mpi_size, new_mpi_size);

                char md5[2][64];
                memset(md5[0], 0, 64);
                memset(md5[1], 0, 64);
                fp = popen(cmd, "r");
                fgets(md5[0], 64, fp);
                fgets(md5[1], 64, fp);
                EXPECT_STREQ(md5[0], md5[1]);

                if(strncmp(md5[0], md5[1], 64) == 0) {
                    char cmd[256];
                    snprintf(cmd, 256, "test -e log/image_global_triangles_15.png && mv log/image_global_triangles_15.png log/image_%s.png", dim1_grid_name[i]);
                    system(cmd);
                    snprintf(cmd, 256, "test -e log/original_input_points.png && mv log/original_input_points.png log/input_%s.png", dim1_grid_name[i]);
                    system(cmd);
                }
            }
            MPI_Barrier(split_world);
        }

    }
};
#endif


const int autogen_grid_size[] = {
                               //100000,
                               //1000000,
                               10000000,
                               //100000,
                               //1000000,
                               //10000000,
                               //100000,
                               //1000000,
                               //10000000,
                               //99540,
                               //999000,
                               //9995082,
                               //99225,
                               //998001,
                             };
const char autogen_grid_name[][64] = { 
                                    //"lonlat_random_global_100000.dat",
                                    //"lonlat_random_global_1000000.dat",
                                    "lonlat_random_global_10000000.dat",
                                    //"lonlat_random_regional_100000.dat",
                                    //"lonlat_random_regional_1000000.dat",
                                    //"lonlat_random_regional_10000000.dat",
                                    //"MonteCarlo_100000.dat",
                                    //"MonteCarlo_1000000.dat",
                                    //"MonteCarlo_10000000.dat",
                                    //"lonlat_uniform_global_100000.dat",
                                    //"lonlat_uniform_global_1000000.dat",
                                    //"lonlat_uniform_global_1000000.dat",
                                    //"lonlat_non-uniform_global_100000.dat",
                                    //"lonlat_non-uniform_global_1000000.dat",
                                  };
const char autogen_grid_path[] = "gridfile/performence_evaluation/%s";
void prepare_autogen_grid(const char grid_name[], int grid_size)
{
    char fullname[128];
    int squeeze_ratio = 0;

    if (mpi_rank == 0) {
        snprintf(fullname, 128, autogen_grid_path, grid_name);
        FILE *fp = fopen(fullname, "r");
        if(!fp) {
            fprintf(stderr, "can not find grid file\n");
            return;
        }

        coord_values[PDLN_LON] = new double[grid_size];
        coord_values[PDLN_LAT] = new double[grid_size];

        if (strstr(grid_name, "lonlat_uniform_global") ||
            strstr(grid_name, "lonlat_non-uniform_global")) {
            for(int i = 0; i < grid_size; i++)
                fscanf(fp, "%lf %lf\n", &coord_values[PDLN_LON][i], &coord_values[PDLN_LAT][i]);
        } else {
            double *x = new double[grid_size];
            double *y = new double[grid_size];
            double *z = new double[grid_size];
            for(int i = 0; i < grid_size; i++)
                fscanf(fp, "%lf %lf %lf\n", &x[i], &y[i], &z[i]);

            for(int i = 0; i < grid_size; i++) {
                coord_values[PDLN_LON][i] = atan2(y[i], x[i]);
                coord_values[PDLN_LAT][i] = asin(z[i]);
            }
            delete[] x;
            delete[] y;
            delete[] z;

            for(int i = 0; i < num_points; i ++) {
                coord_values[PDLN_LON][i] = RADIAN_TO_DEGREE(coord_values[PDLN_LON][i]);
                coord_values[PDLN_LAT][i] = RADIAN_TO_DEGREE(coord_values[PDLN_LAT][i]);
                while(coord_values[PDLN_LON][i] >= 360)
                    coord_values[PDLN_LON][i] -= 360;
                while(coord_values[PDLN_LON][i] < 0)
                    coord_values[PDLN_LON][i] += 360;
            }
        }
        fclose(fp);

        num_points = grid_size;
     
        if(squeeze_ratio > 0) {
            for(int i = 0; i < num_points/squeeze_ratio; i++) {
                coord_values[PDLN_LON][i] = coord_values[PDLN_LON][i*squeeze_ratio];
                coord_values[PDLN_LAT][i] = coord_values[PDLN_LAT][i*squeeze_ratio];
            }
            num_points = num_points/squeeze_ratio;
        }

        //printf("num points: %d\n", num_points);
        assert(!have_redundent_points(coord_values[PDLN_LON], coord_values[PDLN_LAT], num_points));

        min_lon = 0;
        max_lon = 360;
        min_lat = -90;
        max_lat = 90;
        is_cyclic = true;
    }

    MPI_Bcast(&num_points, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (mpi_rank != 0) {
        coord_values[PDLN_LON] = new double[num_points];
        coord_values[PDLN_LAT] = new double[num_points];
    }
    MPI_Bcast(coord_values[PDLN_LON], num_points, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(coord_values[PDLN_LAT], num_points, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&min_lon, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&max_lon, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&min_lat, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&max_lat, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&is_cyclic, 1, MPI_CHAR, 0, MPI_COMM_WORLD);

    assert(sizeof(bool) == 1);
#ifdef TIME_PERF
    printf("[ - ] Total points: %d\n", num_points);
#endif
};


TEST_F(FullProcess, Performance) {
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Comm comm = MPI_COMM_WORLD;

    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    MPI_Comm split_world;
    if (mpi_size/3 > 1)
        MPI_Comm_split(MPI_COMM_WORLD, mpi_rank%3, mpi_rank/3, &split_world);

    for(unsigned i = 0; i < sizeof(autogen_grid_name)/64; i++) {
        printf("processing: %s\n", autogen_grid_name[i]);
        prepare_autogen_grid(autogen_grid_name[i], autogen_grid_size[i]);

        ON_CALL(*mock_process_thread_manager, get_mpi_comm())
            .WillByDefault(Return(comm));

        ON_CALL(*mock_grid_info_manager, get_grid_coord_values(1))
            .WillByDefault(Return(coord_values));

        ON_CALL(*mock_grid_info_manager, get_grid_num_points(1))
            .WillByDefault(Return(num_points));

        ON_CALL(*mock_grid_info_manager, get_grid_boundry(1, _, _, _, _))
            .WillByDefault(Invoke(get_boundry));

        ON_CALL(*mock_grid_info_manager, set_grid_boundry(1, _, _, _, _))
            .WillByDefault(Invoke(set_boundry));

        ON_CALL(*mock_grid_info_manager, is_grid_cyclic(1))
            .WillByDefault(Return(is_cyclic));

        Component* comp;
        comp = new Component(0);
        comp->register_grid(new Grid(1));
        int ret = comp->generate_delaunay_trianglulation(1);
        EXPECT_EQ(ret, 0);
        delete comp;
    }
};
