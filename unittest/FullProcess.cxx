#include "mpi.h"
#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "netcdf.h"
#include "../component.h"
#include "../processing_unit_mgt.h"
#include "../delaunay_grid_decomposition_mgt.h"
#include "../netcdf_utils.h"
#include "../ccpl_utils.h"

#include <cmath>

#define ROUND_VALUE (10000000.0)

extern Grid_info_manager *grid_info_mgr;
extern Process_thread_manager *process_thread_mgr;

FILE* test_log = fopen("all_grid_checksum", "a+");
double max_point_lat;
extern int my_argc;
extern char** my_argv;
class Mock_Process_thread_manager3 : public Process_thread_manager
{
public:
    MOCK_METHOD0(get_mpi_comm, MPI_Comm());
};

class Mock_Grid_info_manager2 : public Grid_info_manager
{
public:
    MOCK_METHOD1(get_grid_coord_values, double**(int));
    MOCK_METHOD1(get_grid_mask, bool*(int));
    MOCK_METHOD1(get_grid_num_points, int(int));
    MOCK_METHOD5(get_grid_boundry, void(int, double*, double*, double*, double*));
    MOCK_METHOD5(set_grid_boundry, void(int, double, double, double, double));
    MOCK_METHOD1(is_grid_cyclic, bool(int));
    MOCK_METHOD4(get_disabled_points_info, void(int, DISABLING_POINTS_METHOD*, int*, void**));
};

using ::testing::Return;
using ::testing::ReturnPointee;
using ::testing::NiceMock;
using ::testing::_;
using ::testing::Invoke;
using ::testing::ExitedWithCode;

static int mpi_rank = -1;
static int mpi_size = 0;
static double *coord_values[2] = {NULL, NULL};
static bool *grid_mask = NULL;
static int num_points = 0;
static double min_lat, max_lat, min_lon, max_lon;
static bool is_cyclic = false;
static MPI_Comm comm = MPI_COMM_WORLD;
static DISABLING_POINTS_METHOD disabling_method = NO_DISABLED_POINTS;
static int disabling_num = 0;
static void* disabling_data = NULL;

static void get_boundry(int grid_id, double* mi_lon, double* ma_lon, double* mi_lat, double* ma_lat);
static void set_boundry(int grid_id, double mi_lon, double ma_lon, double mi_lat, double ma_lat);
static void get_disabled(int grid_id, DISABLING_POINTS_METHOD* method, int* num, void** data);


class FullProcess : public ::testing::Test
{
public:
    virtual void SetUp() {
        mock_process_thread_manager = new NiceMock<Mock_Process_thread_manager3>;
        mock_grid_info_manager = new NiceMock<Mock_Grid_info_manager2>;

        process_thread_mgr = mock_process_thread_manager;
        grid_info_mgr = mock_grid_info_manager;

        ON_CALL(*mock_grid_info_manager, get_grid_coord_values(1))
            .WillByDefault(Return(coord_values));

        ON_CALL(*mock_grid_info_manager, get_grid_boundry(1, _, _, _, _))
            .WillByDefault(Invoke(get_boundry));

        ON_CALL(*mock_grid_info_manager, get_grid_num_points(1))
            .WillByDefault(ReturnPointee(&num_points));

        ON_CALL(*mock_grid_info_manager, is_grid_cyclic(1))
            .WillByDefault(ReturnPointee(&is_cyclic));

        ON_CALL(*mock_process_thread_manager, get_mpi_comm())
            .WillByDefault(ReturnPointee(&comm));

        ON_CALL(*mock_grid_info_manager, set_grid_boundry(1, _, _, _, _))
            .WillByDefault(Invoke(set_boundry));

        ON_CALL(*mock_grid_info_manager, get_disabled_points_info(1, _, _, _))
            .WillByDefault(Invoke(get_disabled));
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

    is_cyclic = true;
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
        if(std::abs(coord_values[PDLN_LON][i] - 360.0) < PDLN_ABS_TOLERANCE) {
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

    is_cyclic = true;
}


void prepare_big_three_polar_grid()
{
    int num_dims;
    int *dim_size_ptr;
    int field_size;
    int field_size2;
    void *coord_buf0, *coord_buf1;
    bool squeeze = false;

    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    if (mpi_rank == 0) {
        read_file_field_as_double("gridfile/many_types_of_grid/big_three_poles.nc", "x_T", &coord_buf0, &num_dims, &dim_size_ptr, &field_size);
        delete dim_size_ptr;
        read_file_field_as_double("gridfile/many_types_of_grid/big_three_poles.nc", "y_T", &coord_buf1, &num_dims, &dim_size_ptr, &field_size2);
        delete dim_size_ptr;
        assert(field_size == field_size2);
        num_points = field_size;
        coord_values[PDLN_LON] = (double*)coord_buf0;
        coord_values[PDLN_LAT] = (double*)coord_buf1;

        for(int i = 0; i < num_points; i++) {
            while(coord_values[PDLN_LON][i] < 0.0)
                coord_values[PDLN_LON][i] += 360.0;
            while(coord_values[PDLN_LON][i] >= 360.0)
                coord_values[PDLN_LON][i] -= 360.0;
        }

        delete_redundent_points(coord_values[PDLN_LON], coord_values[PDLN_LAT], num_points);

        if(squeeze) {
            for(int i = 0; i < num_points/100; i++) {
                coord_values[PDLN_LON][i] = coord_values[PDLN_LON][i*100];
                coord_values[PDLN_LAT][i] = coord_values[PDLN_LAT][i*100];
            }
            num_points /= 100;
        }
    }

    MPI_Bcast(&num_points, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (mpi_rank != 0) {
        coord_values[PDLN_LON] = new double[num_points];
        coord_values[PDLN_LAT] = new double[num_points];
    }

    MPI_Bcast(coord_values[PDLN_LON], num_points, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(coord_values[PDLN_LAT], num_points, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    min_lon =   0.0;
    max_lon = 360.0;
    min_lat = -83.0;
    max_lat =  90.0;

    is_cyclic = true;

    double* circle = new double[6];
    disabling_method = DISABLE_POINTS_BY_RANGE;
    disabling_num = 2;
    circle[0] = 80;
    circle[1] = 65.6;
    circle[2] = 0.75;
    circle[3] = 260;
    circle[4] = 65.6;
    circle[5] = 0.75;
    disabling_data = circle;
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

    is_cyclic = true;
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

    is_cyclic = true;
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

    is_cyclic = true;
}
#endif


TEST_F(FullProcess, OLD_Basic) {
    comm = MPI_COMM_WORLD;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    prepare_grid();

    Patcc* comp;
    comp = new Patcc(0);
    comp->register_grid(new Grid(1));
    comp->generate_delaunay_trianglulation(1, true);

    delete comp;

<<<<<<< HEAD
//    if (mpi_size/3 > 1)
//        MPI_Comm_split(MPI_COMM_WORLD, mpi_rank%3, mpi_rank/3, &comm);
//
//    if (mpi_size/3 > 1 && mpi_rank%3 == 0) {
//        comp = new Component(1);
//        comp->register_grid(new Grid(1));
//        comp->generate_delaunay_trianglulation(1, true);
//
//        delete comp;
//
//        int new_mpi_size;
//        MPI_Comm_size(comm, &new_mpi_size);
//
=======
    if (mpi_size/3 > 1)
        MPI_Comm_split(MPI_COMM_WORLD, mpi_rank%3, mpi_rank/3, &comm);

    if (mpi_size/3 > 1 && mpi_rank%3 == 0) {
        comp = new Patcc(1);
        comp->register_grid(new Grid(1));
        comp->generate_delaunay_trianglulation(1, true);

        delete comp;

        int new_mpi_size;
        MPI_Comm_size(comm, &new_mpi_size);

>>>>>>> mask
        if (mpi_rank == 0) {
            FILE *fp;
            char cmd[] = "md5sum log/global_triangles_* | awk -F\" \" '{print $1}'";

            char md5[2][64];
            memset(md5[0], 0, 64);
//            memset(md5[1], 0, 64);
            fp = popen(cmd, "r");
            fgets(md5[0], 64, fp);
//            fgets(md5[1], 64, fp);
//            EXPECT_STREQ(md5[0], md5[1]);
			fprintf(test_log, "Basic %s\n",md5[0]);        
			fflush(test_log);
			char mv[] = "mv log/global_triangles_* log/check/global_triangles_Basic";
			fp = popen(mv, "r");
        }
//    }
};


#ifdef NETCDF
TEST_F(FullProcess, OLD_LatLonGrid) {
    comm = MPI_COMM_WORLD;

    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);


    prepare_latlon_grid();

    Patcc* comp;
    comp = new Patcc(0);
    comp->register_grid(new Grid(1));
    comp->generate_delaunay_trianglulation(1, true);

    delete comp;

<<<<<<< HEAD
//    if (mpi_size/3 > 1)
//        MPI_Comm_split(MPI_COMM_WORLD, mpi_rank%3, mpi_rank/3, &comm);
//
//    if (mpi_size/3 > 1 && mpi_rank%3 == 0) {
//        comp = new Component(1);
//        comp->register_grid(new Grid(1));
//        comp->generate_delaunay_trianglulation(1, true);
//
//        delete comp;
//
//        int new_mpi_size;
//        MPI_Comm_size(comm, &new_mpi_size);
//
=======
    if (mpi_size/3 > 1)
        MPI_Comm_split(MPI_COMM_WORLD, mpi_rank%3, mpi_rank/3, &comm);

    if (mpi_size/3 > 1 && mpi_rank%3 == 0) {
        comp = new Patcc(1);
        comp->register_grid(new Grid(1));
        comp->generate_delaunay_trianglulation(1, true);

        delete comp;

        int new_mpi_size;
        MPI_Comm_size(comm, &new_mpi_size);

>>>>>>> mask
        if (mpi_rank == 0) {
            FILE *fp;
            char cmd[] = "md5sum log/global_triangles_* | awk -F\" \" '{print $1}'";

            char md5[2][64];
            memset(md5[0], 0, 64);
//            memset(md5[1], 0, 64);
            fp = popen(cmd, "r");
            fgets(md5[0], 64, fp);
//            fgets(md5[1], 64, fp);
//            EXPECT_STREQ(md5[0], md5[1]);
            fprintf(test_log, "LatLonGrid %s\n",md5[0]);
			fflush(test_log);
			char mv[] = "mv log/global_triangles_* log/check/global_triangles_LatLonGrid";
            fp = popen(mv, "r");
        }
//    }
};


TEST_F(FullProcess, OLD_LatLonSinglePolar) {
    comm = MPI_COMM_WORLD;

    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    prepare_latlon_singlepolar();

    Patcc* comp;
    comp = new Patcc(0);
    comp->register_grid(new Grid(1));
    comp->generate_delaunay_trianglulation(1, true);

    delete comp;

<<<<<<< HEAD
//    if (mpi_size/3 > 1)
//        MPI_Comm_split(MPI_COMM_WORLD, mpi_rank%3, mpi_rank/3, &comm);
//
//    if (mpi_size/3 > 1 && mpi_rank%3 == 0) {
//        comp = new Component(1);
//        comp->register_grid(new Grid(1));
//        comp->generate_delaunay_trianglulation(1, true);
//
//        delete comp;
//
//        int new_mpi_size;
//        MPI_Comm_size(comm, &new_mpi_size);
//
=======
    if (mpi_size/3 > 1)
        MPI_Comm_split(MPI_COMM_WORLD, mpi_rank%3, mpi_rank/3, &comm);

    if (mpi_size/3 > 1 && mpi_rank%3 == 0) {
        comp = new Patcc(1);
        comp->register_grid(new Grid(1));
        comp->generate_delaunay_trianglulation(1, true);

        delete comp;

        int new_mpi_size;
        MPI_Comm_size(comm, &new_mpi_size);

>>>>>>> mask
        if (mpi_rank == 0) {
            FILE *fp;
            char cmd[] = "md5sum log/global_triangles_* | awk -F\" \" '{print $1}'";

            char md5[2][64];
            memset(md5[0], 0, 64);
//            memset(md5[1], 0, 64);
            fp = popen(cmd, "r");
            fgets(md5[0], 64, fp);
//            fgets(md5[1], 64, fp);
//            EXPECT_STREQ(md5[0], md5[1]);
			fprintf(test_log,"LatLonSinglePolar %s\n",md5[0]);           
			fflush(test_log);
			char mv[] = "mv log/global_triangles_* log/check/global_triangles_LatLonSinglePolar";
            fp = popen(mv, "r");
        }
//    }
};



TEST_F(FullProcess, OLD_LatLonMutiPolars) {
    comm = MPI_COMM_WORLD;

    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    prepare_latlon_mutipolars();

    Patcc* comp;
    comp = new Patcc(0);
    comp->register_grid(new Grid(1));
    comp->generate_delaunay_trianglulation(1, true);

    delete comp;

<<<<<<< HEAD
//    if (mpi_size/3 > 1)
//        MPI_Comm_split(MPI_COMM_WORLD, mpi_rank%3, mpi_rank/3, &comm);
//
//    if (mpi_size/3 > 1 && mpi_rank%3 == 0) {
//        comp = new Component(1);
//        comp->register_grid(new Grid(1));
//        comp->generate_delaunay_trianglulation(1, true);
//
//        delete comp;
//
//        int new_mpi_size;
//        MPI_Comm_size(comm, &new_mpi_size);
//
=======
    if (mpi_size/3 > 1)
        MPI_Comm_split(MPI_COMM_WORLD, mpi_rank%3, mpi_rank/3, &comm);

    if (mpi_size/3 > 1 && mpi_rank%3 == 0) {
        comp = new Patcc(1);
        comp->register_grid(new Grid(1));
        comp->generate_delaunay_trianglulation(1, true);

        delete comp;

        int new_mpi_size;
        MPI_Comm_size(comm, &new_mpi_size);

>>>>>>> mask
        if (mpi_rank == 0) {
            FILE *fp;
            char cmd[] = "md5sum log/global_triangles_* | awk -F\" \" '{print $1}'";

            char md5[2][64];
            memset(md5[0], 0, 64);
//            memset(md5[1], 0, 64);
            fp = popen(cmd, "r");
            fgets(md5[0], 64, fp);
//            fgets(md5[1], 64, fp);
//            EXPECT_STREQ(md5[0], md5[1]);
			fprintf(test_log,"LatLonMutiPolars %s\n",md5[0]);           
			fflush(test_log);
			char mv[] = "mv log/global_triangles_* log/check/global_triangles_LatLonMutiPolar";
            fp = popen(mv, "r");

		}
//    }
};


TEST_F(FullProcess, OLD_ThreePolar) {
    comm = MPI_COMM_WORLD;

    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    prepare_three_polar_grid();

    Patcc* comp;
    comp = new Patcc(0);
    comp->register_grid(new Grid(1));
    comp->generate_delaunay_trianglulation(1, true);

    delete comp;

<<<<<<< HEAD
//    if (mpi_size/3 > 1)
//        MPI_Comm_split(MPI_COMM_WORLD, mpi_rank%3, mpi_rank/3, &comm);
//
//    if (mpi_size/3 > 1 && mpi_rank%3 == 0) {
//        comp = new Component(1);
//        comp->register_grid(new Grid(1));
//        comp->generate_delaunay_trianglulation(1, true);
//
//        delete comp;
//
//        int new_mpi_size;
//        MPI_Comm_size(comm, &new_mpi_size);
//
=======
    if (mpi_size/3 > 1)
        MPI_Comm_split(MPI_COMM_WORLD, mpi_rank%3, mpi_rank/3, &comm);

    if (mpi_size/3 > 1 && mpi_rank%3 == 0) {
        comp = new Patcc(1);
        comp->register_grid(new Grid(1));
        comp->generate_delaunay_trianglulation(1, true);

        delete comp;

        int new_mpi_size;
        MPI_Comm_size(comm, &new_mpi_size);

>>>>>>> mask
        if (mpi_rank == 0) {
            FILE *fp;
            char cmd[] = "md5sum log/global_triangles_* | awk -F\" \" '{print $1}'";

            char md5[2][64];
            memset(md5[0], 0, 64);
//            memset(md5[1], 0, 64);
            fp = popen(cmd, "r");
            fgets(md5[0], 64, fp);
//            fgets(md5[1], 64, fp);
//            EXPECT_STREQ(md5[0], md5[1]);
			fprintf(test_log,"ThreePolar %s\n",md5[0]);
			fflush(test_log);
			char mv[] = "mv log/global_triangles_* log/check/global_triangles_ThreePolar";
            fp = popen(mv, "r");

		}
//    }
};


TEST_F(FullProcess, OLD_ThreePolarBig) {
    comm = MPI_COMM_WORLD;

    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    prepare_big_three_polar_grid();

<<<<<<< HEAD
    Component* comp;
    comp = new Component(0);
=======
    printf("processing...\n");
    Patcc* comp;
    comp = new Patcc(0);
>>>>>>> mask
    comp->register_grid(new Grid(1));
    comp->generate_delaunay_trianglulation(1, true);

    delete comp;

<<<<<<< HEAD
//    if (mpi_size/3 > 1)
//        MPI_Comm_split(MPI_COMM_WORLD, mpi_rank%3, mpi_rank/3, &comm);
//
//    if (mpi_size/3 > 1 && mpi_rank%3 == 0) {
//        printf("part of processes...\n");
//        comp = new Component(1);
//        comp->register_grid(new Grid(1));
//        comp->generate_delaunay_trianglulation(1, true);
//
//        delete comp;
//
/*        if (mpi_rank == 0) {
=======
    if (mpi_size/3 > 1)
        MPI_Comm_split(MPI_COMM_WORLD, mpi_rank%3, mpi_rank/3, &comm);

    if (mpi_size/3 > 1 && mpi_rank%3 == 0) {
        printf("part of processes...\n");
        comp = new Patcc(1);
        comp->register_grid(new Grid(1));
        comp->generate_delaunay_trianglulation(1, true);

        delete comp;

        if (mpi_rank == 0) {
>>>>>>> mask
            FILE *fp;
            char cmd[] = "md5sum log/global_triangles_* | awk -F\" \" '{print $1}'";

            char md5[2][64];
            memset(md5[0], 0, 64);
//          memset(md5[1], 0, 64);
            fp = popen(cmd, "r");
            fgets(md5[0], 64, fp);
//          fgets(md5[1], 64, fp);
//          EXPECT_STREQ(md5[0], md5[1]);	
			fprintf(test_log,"ThreePolarBig %s\n",md5[0]);
			fflush(test_log);
			char mv[] = "mv log/global_triangles_* log/check/global_triangles_ThreePolarBig";
			fp = popen(mv, "r");

		}
		*/
        MPI_Barrier(MPI_COMM_WORLD);
//    }
};


#define CHECK_PARALLEL_CONSISTENCY (true)
const char dim1_grid_path[] = "gridfile/many_types_of_grid/one_dimension/%s";
const char dim1_grid_name[][64] = {
    "ne30np4-t2.nc",
    "ne60np4_pentagons_100408.nc",
    "gx3v5_Present_DP_x3.nc",
    "Gamil_360x180_Grid.nc",
    "licom_eq1x1_degree_Grid.nc",
    "licom_gr1x1_degree_Grid.nc",
    "LICOM_P5_Grid.nc",
    "ll1deg_grid.nc",
    "ll2.5deg_grid.nc",
    "R42_Gaussian_Grid.nc",
    "T85_Gaussian_Grid.nc",
    "Gamil_128x60_Grid.nc",
    "fv1.9x2.5_050503.nc",
    "T42_Gaussian_Grid.nc",
    "T42_grid.nc",
    "T62_Gaussian_Grid.nc",
    "ar9v4_100920.nc",
    "wr50a_090301.nc",
    "Version_3_of_Greenland_pole_x1_T-grid.nc",
    "R05_Grid.nc",
    "CUBE_grid_2.5.nc",
    "CUBE_grid_1.nc",
    "CUBE_grid_0.3.nc",
    "CUBE_grid_0.1.nc",
};
const char dim1_global_grid_name[][64] = {
    "ne30np4-t2.nc",
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
    "CUBE_grid_2.5.nc",
    "CUBE_grid_1.nc",
    "CUBE_grid_0.3.nc",
    "CUBE_grid_0.1.nc",
};


void save_dim1_grid(const char grid_name[]);
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

    grid_mask = NULL;
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
        if (strncmp(grid_name, "mom_h2d_T_grid@mom.nc", 64) == 0) {
            void* raw_mask;
            int* int_mask;
            read_file_field_as_int(fullname, "grid_imask", &raw_mask, &num_dims, &dim_size_ptr, &field_size2, NULL);
            int_mask = (int*)raw_mask;
            grid_mask = new bool[field_size2];
            for (int i = 0; i < field_size2; i ++)
                grid_mask[i] = int_mask[i];
        }

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
        if (strstr(grid_name, "CUBE_grid_")) {
            for (int i = 0; i < num_points; i++) {
                coord_values[PDLN_LON][i] = round(coord_values[PDLN_LON][i]*ROUND_VALUE)/ROUND_VALUE;
                coord_values[PDLN_LAT][i] = round(coord_values[PDLN_LAT][i]*ROUND_VALUE)/ROUND_VALUE;
                while(coord_values[PDLN_LON][i] >= 360)
                    coord_values[PDLN_LON][i] -= 360;
                while(coord_values[PDLN_LON][i] < 0)
                    coord_values[PDLN_LON][i] += 360;
            }
            delete_redundent_points(coord_values[PDLN_LON], coord_values[PDLN_LAT], num_points);
        } else
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
        if(strncmp(grid_name, "mom_h2d_T_grid@mom.nc", 64) == 0) {
            min_lon = 0;
            max_lon = 360;
            min_lat = -83;
            max_lat = 90;
        }

        if(fabs((max_lon - min_lon) - 360) < 0.5)
            is_cyclic = true;
        else
            is_cyclic = false;

        for (int i = 0; i < num_points; i++) {
            coord_values[PDLN_LON][i] = round(coord_values[PDLN_LON][i]*ROUND_VALUE)/ROUND_VALUE;
            coord_values[PDLN_LAT][i] = round(coord_values[PDLN_LAT][i]*ROUND_VALUE)/ROUND_VALUE;
            while(coord_values[PDLN_LON][i] >= 360)
                coord_values[PDLN_LON][i] -= 360;
            while(coord_values[PDLN_LON][i] < 0)
                coord_values[PDLN_LON][i] += 360;
        }

        assert(!have_redundent_points(coord_values[PDLN_LON], coord_values[PDLN_LAT], num_points));
        //save_dim1_grid(grid_name);
    }

    MPI_Bcast(&num_points, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (mpi_rank != 0) {
        coord_values[PDLN_LON] = new double[num_points];
        coord_values[PDLN_LAT] = new double[num_points];
        if (strncmp(grid_name, "mom_h2d_T_grid@mom.nc", 64) == 0)
            grid_mask = new bool[num_points];
    }

    MPI_Bcast(coord_values[PDLN_LON], num_points, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(coord_values[PDLN_LAT], num_points, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&min_lon, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&max_lon, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&min_lat, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&max_lat, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&is_cyclic, 1, MPI_CHAR, 0, MPI_COMM_WORLD);

    if (strncmp(grid_name, "mom_h2d_T_grid@mom.nc", 64) == 0)
        MPI_Bcast(grid_mask, num_points, MPI_CHAR, 0, MPI_COMM_WORLD);

    //double* circle = new double[3];
    //circle[0] = 270;
    //circle[1] = 45;
    //circle[2] = 10;
    //disabling_num = 1;
    //disabling_data = circle;
    disabling_method = NO_DISABLED_POINTS;

    assert(sizeof(bool) == 1);
#ifdef TIME_PERF
    printf("[ - ] Total points: %d\n", num_points);
#endif
};


const char dim1_grid_save_path[] = "gridfile/for_scvt/%s";
void save_dim1_grid(const char grid_name[])
{
    char fullname[128];

    snprintf(fullname, 128, dim1_grid_save_path, grid_name);

    FILE* fp = fopen(fullname, "w");
    for (int i = 0; i < num_points; i++) {
        double x = cos(DEGREE_TO_RADIAN(coord_values[PDLN_LON][i]))*cos(DEGREE_TO_RADIAN(coord_values[PDLN_LAT][i]));
        double y = sin(DEGREE_TO_RADIAN(coord_values[PDLN_LON][i]))*cos(DEGREE_TO_RADIAN(coord_values[PDLN_LAT][i]));
        double z = sin(DEGREE_TO_RADIAN(coord_values[PDLN_LAT][i]));
        fprintf(fp, "%.20lf %.20lf %.20lf\n", x, y, z);
    }
    fclose(fp);
    printf("File writen\n");
}


#include <unistd.h>
TEST_F(FullProcess, ManyTypesOfGrids) {
    MPI_Barrier(MPI_COMM_WORLD);
    printf("pid: %d\n", getpid());
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    ON_CALL(*mock_grid_info_manager, get_grid_mask(1))
        .WillByDefault(ReturnPointee(&grid_mask));

    MPI_Barrier(MPI_COMM_WORLD);
    printf("processing: %s\n", my_argv[1]);
    comm = MPI_COMM_WORLD;

    prepare_dim1_grid(my_argv[1]);

<<<<<<< HEAD
    Component* comp;
    comp = new Component(0);
    comp->register_grid(new Grid(1));
    int ret = comp->generate_delaunay_trianglulation(1, true);
    EXPECT_EQ(ret, 0);
    delete comp;
=======
    for(unsigned i = 0; i < sizeof(dim1_grid_name)/64; i++) {
        MPI_Barrier(MPI_COMM_WORLD);
        printf("processing: %s\n", dim1_grid_name[i]);
        comm = MPI_COMM_WORLD;
        prepare_dim1_grid(dim1_grid_name[i]);

        Patcc* comp;
        comp = new Patcc(0);
        comp->register_grid(new Grid(1));
        int ret = comp->generate_delaunay_trianglulation(1, true);
        EXPECT_EQ(ret, 0);
        delete comp;

        if (CHECK_PARALLEL_CONSISTENCY && mpi_size/3 > 1 && mpi_rank%3 == 0) {
            printf("spliting world: %s\n", dim1_grid_name[i]);
            comm = split_world;
            MPI_Barrier(split_world);
            comp = new Patcc(1);
            comp->register_grid(new Grid(1));
            int ret = comp->generate_delaunay_trianglulation(1, true);
            EXPECT_EQ(ret, 0);

            delete comp;

            if (mpi_rank == 0) {
                FILE *fp;
                char cmd[] = "md5sum log/global_triangles_* | awk -F\" \" '{print $1}'";

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
>>>>>>> mask

    MPI_Barrier(MPI_COMM_WORLD);
};
#endif


const int autogen_grid_size[] = {
                               10000,
                               100000,
                               1000000,
                               10000000,

                               100000,
                               1000000,
                               10000000,

                               99225,
                               998001,

                               64800,
                               720000,
                               6480000,
                               10800,

                               0,
                               0,
                               0,
                               0,
                             };
const char autogen_grid_name[][64] = { 
                                    "lonlat_random_global_10000.dat",
                                    "lonlat_random_global_100000.dat",
                                    "lonlat_random_global_1000000.dat",
                                    "lonlat_random_global_10000000.dat",

                                    "MonteCarlo_100000.dat",
                                    "MonteCarlo_1000000.dat",
                                    "MonteCarlo_10000000.dat",

                                    "lonlat_non-uniform_global_100000.dat",
                                    "lonlat_non-uniform_global_1000000.dat",

                                    "lonlat_grid_1.dat",
                                    "lonlat_grid_0.3.dat",
                                    "lonlat_grid_0.1.dat",
                                    "lonlat_grid_small.dat",

                                    "CUBE_grid_2.5.nc",
                                    "CUBE_grid_1.nc",
                                    "CUBE_grid_0.3.nc",
                                    "CUBE_grid_0.1.nc",
                                  };
const char autogen_grid_path[] = "gridfile/performence_evaluation/%s";
void prepare_autogen_grid(const char grid_name[], int grid_size)
{
    char fullname[128];
    int squeeze_ratio = 0;

    if (mpi_rank == 0) {
        snprintf(fullname, 128, autogen_grid_path, grid_name);

        if (strstr(grid_name, "CUBE_grid_")) {
            int num_dims;
            int *dim_size_ptr;
            int field_size, field_size2;
            void *coord_buf0, *coord_buf1;
            char lon_unit[32];
            char lat_unit[32];

            snprintf(fullname, 128, autogen_grid_path, grid_name);
            read_file_field_as_double(fullname, "grid_center_lon", &coord_buf0, &num_dims, &dim_size_ptr, &field_size, lon_unit);
            delete dim_size_ptr;
            read_file_field_as_double(fullname, "grid_center_lat", &coord_buf1, &num_dims, &dim_size_ptr, &field_size2, lat_unit);
            delete dim_size_ptr;

            ASSERT_EQ(field_size, field_size2);
            num_points = field_size;
            coord_values[PDLN_LON] = (double*)coord_buf0;
            coord_values[PDLN_LAT] = (double*)coord_buf1;

            for (int i = 0; i < num_points; i++) {
                coord_values[PDLN_LON][i] = round(coord_values[PDLN_LON][i]*ROUND_VALUE)/ROUND_VALUE;
                coord_values[PDLN_LAT][i] = round(coord_values[PDLN_LAT][i]*ROUND_VALUE)/ROUND_VALUE;
                while(coord_values[PDLN_LON][i] >= 360)
                    coord_values[PDLN_LON][i] -= 360;
                while(coord_values[PDLN_LON][i] < 0)
                    coord_values[PDLN_LON][i] += 360;
            }
            delete_redundent_points(coord_values[PDLN_LON], coord_values[PDLN_LAT], num_points);

            /*
            FILE* fp = fopen("CUBE_grid_for_SCVT_0.1.dat", "w");
            for (int i = 0; i < num_points; i++) {
                double x = cos(DEGREE_TO_RADIAN(coord_values[PDLN_LON][i]))*cos(DEGREE_TO_RADIAN(coord_values[PDLN_LAT][i]));
                double y = sin(DEGREE_TO_RADIAN(coord_values[PDLN_LON][i]))*cos(DEGREE_TO_RADIAN(coord_values[PDLN_LAT][i]));
                double z = sin(DEGREE_TO_RADIAN(coord_values[PDLN_LAT][i]));
                fprintf(fp, "%.20lf %.20lf %.20lf\n", x, y, z);
            }
            fclose(fp);
            printf("File writen\n");
            */
        } else {
            FILE *fp = fopen(fullname, "r");
            if(!fp) {
                fprintf(stderr, "can not find grid file\n");
                return;
            }

            coord_values[PDLN_LON] = new double[grid_size];
            coord_values[PDLN_LAT] = new double[grid_size];

            if (strstr(grid_name, "lonlat_uniform_global") ||
                strstr(grid_name, "lonlat_non-uniform_global") ||
                strstr(grid_name, "lonlat_grid_")) {
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

                for(int i = 0; i < grid_size; i ++) {
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
        }
     
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
    comm = MPI_COMM_WORLD;
	
	int grid_size=0;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

	for(unsigned i = 0; i < sizeof(autogen_grid_name)/64; i++) {
		if(strcmp(my_argv[1],autogen_grid_name[i])==0){
			grid_size = autogen_grid_size[i];
			break;
		}
	}
    printf("processing: %s\n", my_argv[1]);
    prepare_autogen_grid(my_argv[1],grid_size);

<<<<<<< HEAD
    Component* comp;
    comp = new Component(0);
    comp->register_grid(new Grid(1));
    int ret = comp->generate_delaunay_trianglulation(1, true);
    EXPECT_EQ(ret, 0);
    delete comp;
=======
        Patcc* comp;
        comp = new Patcc(0);
        comp->register_grid(new Grid(1));
        int ret = comp->generate_delaunay_trianglulation(1);
        EXPECT_EQ(ret, 0);
        delete comp;
    }
>>>>>>> mask
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

static void get_disabled(int grid_id, DISABLING_POINTS_METHOD* method, int* num, void** data)
{
    *method = disabling_method;
    *num = disabling_num;
    *data = disabling_data;
}

const char real_grid_name[][64]={
		"ice_grid_at_ice_comp@cesm@ICE.nc",
		"ocn_grid_at_ocn_comp@cesm@OCN.nc",
		"MITgcm_H2D_grid@mitgcm.nc",
		"gigc_H2D_grid@GIGC.nc",
		"gigc_H2D_grid_temp@GIGC.nc",
		"pom_grid_via_local@pom@POM_TOP_d01.nc",
		"wrf_H2D_grid_no_mask@wrf.nc",
		"wrf_H2D_grid_ocn_mask@wrf.nc",
		"wrf_grid_via_local@wrf@WRF_TOP_d01.nc",
		"wwatch3_H2D_grid@wwatch3.nc",
		"wwatch3_H2D_grid_temp@wwatch3.nc",
		"mpas_H2D_grid@mpas.nc",
		"mpas_H2D_grid_mask@mpas.nc",
		"rof_grid_at_rof_comp@cesm@ROF.nc",
		"lnd_grid_at_lnd_comp@cesm@LND.nc",
		"atm_grid_at_gamil_comp@gamil.nc",
		"atm_grid_at_atm_comp@cesm@ATM.nc",
		"cfs_atm_H2D_grid@cfs_atm.nc",
		"gamil_H2D_grid@gamil.nc",
		"gamil_H2D_grid_temp@gamil.nc",
		"cfs_ocn_H2D_grid@cfs_ocn.nc",
		"bccam_2D_grid@bccam.nc",
		"bccam_2D_grid_temp@bccam.nc",
		"bccavim_lnd_grid@bccavim.nc",
		"bccavim_river_grid@bccavim.nc",
		"bccavim_river_grid_temp@bccavim.nc",
		"cpl_atm_grid@cpl.nc",
		"cpl_ice_grid@cpl.nc",
		"cpl_lnd_grid@cpl.nc",
		"cpl_lnd_grid_temp@cpl.nc",
		"cpl_ocn_grid@cpl.nc",
		"cpl_river_grid@cpl.nc",
		"mom_ocn_grid@mom.nc",
		"mom_ocn_grid_temp@mom.nc",
		"sis_ice_grid@sis.nc",
		"sis_ice_grid_temp@sis.nc",
		"grapes_h2d_ocn_grid@grapes.nc",
		"grapes_h2d_T_grid@grapes.nc",

};

const char real_global_grid_name[][64]={
	    "mpas_H2D_grid@mpas.nc",
		"mpas_H2D_grid_mask@mpas.nc",
        "rof_grid_at_rof_comp@cesm@ROF.nc",
        "lnd_grid_at_lnd_comp@cesm@LND.nc",
		"atm_grid_at_atm_comp@cesm@ATM.nc",
        "cfs_atm_H2D_grid@cfs_atm.nc",
        "gamil_H2D_grid@gamil.nc",
        "gamil_H2D_grid_temp@gamil.nc",
		"cfs_ocn_H2D_grid@cfs_ocn.nc",
		
};

void prepare_real_grid(const char grid_name[]){
	int num_dims;
	int *dim_size_ptr;
	int field_size, field_size2;
	void *coord_buf0, *coord_buf1;
	char lon_unit[32];
	char lat_unit[32];
	char cyclic[32];
	const char* field_lat_name;
	const char* field_lon_name;
	int ncfile_id, variable_lon_id, variable_lat_id;
	nc_type xtype;
	char file_name[64]="coupler_grid/";
	sprintf(file_name,"%s%s",file_name,grid_name);
	if(strcmp(grid_name,"atm_grid_at_atm_comp@cesm@ATM.nc")==0 || 
	   strcmp(grid_name,"ice_grid_at_ice_comp@cesm@ICE.nc")==0 ||
	   strcmp(grid_name,"lnd_grid_at_lnd_comp@cesm@LND.nc")==0 ||
	   strcmp(grid_name,"ocn_grid_at_ocn_comp@cesm@OCN.nc")==0 ||
	   strcmp(grid_name,"rof_grid_at_rof_comp@cesm@ROF.nc")==0){
		field_lat_name="lat";
		field_lon_name="lon";
	}
	else{
		field_lat_name="grid_center_lat";
		field_lon_name="grid_center_lon";
	}
	
	report_nc_error(nc_open(file_name, NC_NOWRITE, &ncfile_id));
	report_nc_error(nc_inq_varid(ncfile_id, field_lon_name, &variable_lon_id));
	report_nc_error(nc_inq_varid(ncfile_id, field_lat_name, &variable_lat_id));
	report_nc_error(nc_inq_var(ncfile_id, variable_lon_id, NULL, &xtype, NULL, NULL, NULL));
	
	report_nc_error(nc_get_att_double(ncfile_id,NC_GLOBAL,"min_lon",&min_lon));
	report_nc_error(nc_get_att_double(ncfile_id,NC_GLOBAL,"min_lat",&min_lat));
	report_nc_error(nc_get_att_double(ncfile_id,NC_GLOBAL,"max_lon",&max_lon));
	report_nc_error(nc_get_att_double(ncfile_id,NC_GLOBAL,"max_lat",&max_lat));
	report_nc_error(nc_get_att_text(ncfile_id,NC_GLOBAL,"cyclic_or_acyclic",cyclic));

	if(strcmp(cyclic,"acyclic")==0)
		is_cyclic = false;
	else
		is_cyclic = true;

	if(xtype == NC_DOUBLE){
		read_file_field_as_double(file_name, field_lon_name, &coord_buf0, &num_dims, &dim_size_ptr, &field_size, lon_unit);
		delete dim_size_ptr;
		read_file_field_as_double(file_name, field_lat_name, &coord_buf1, &num_dims, &dim_size_ptr, &field_size2, lat_unit);
		delete dim_size_ptr;
		ASSERT_EQ(field_size, field_size2);
		num_points = field_size;
		coord_values[PDLN_LON] = (double*)coord_buf0;
		coord_values[PDLN_LAT] = (double*)coord_buf1;
		for (int i = 0; i < num_points; i++) {
		    coord_values[PDLN_LON][i] = round(coord_values[PDLN_LON][i]*ROUND_VALUE)/ROUND_VALUE;
		    coord_values[PDLN_LAT][i] = round(coord_values[PDLN_LAT][i]*ROUND_VALUE)/ROUND_VALUE;
		    while(coord_values[PDLN_LON][i] >= 360)
		        coord_values[PDLN_LON][i] -= 360;
		    while(coord_values[PDLN_LON][i] < 0)
		        coord_values[PDLN_LON][i] += 360;
        if(coord_values[PDLN_LON][i] < min_lon) min_lon = coord_values[PDLN_LON][i];
        if(coord_values[PDLN_LON][i] > max_lon) max_lon = coord_values[PDLN_LON][i];
        if(coord_values[PDLN_LAT][i] < min_lat) min_lat = coord_values[PDLN_LAT][i];
        if(coord_values[PDLN_LAT][i] > max_lat) max_lat = coord_values[PDLN_LAT][i];
		
		}
		delete_redundent_points(coord_values[PDLN_LON], coord_values[PDLN_LAT], num_points);
	}
	else if(xtype == NC_FLOAT){
		read_file_field_as_float(file_name, "grid_center_lon", &coord_buf0, &num_dims, &dim_size_ptr, &field_size);
		delete dim_size_ptr;
		read_file_field_as_float(file_name, "grid_center_lat", &coord_buf1, &num_dims, &dim_size_ptr, &field_size2);
		delete dim_size_ptr;
		ASSERT_EQ(field_size, field_size2);
		num_points = field_size;
		coord_values[PDLN_LON] = (double*)coord_buf0;
		coord_values[PDLN_LAT] = (double*)coord_buf1;

		for (int i = 0; i < num_points; i++) {
		    coord_values[PDLN_LON][i] = (double)round(coord_values[PDLN_LON][i]*ROUND_VALUE)/ROUND_VALUE;
		    coord_values[PDLN_LAT][i] = (double)round(coord_values[PDLN_LAT][i]*ROUND_VALUE)/ROUND_VALUE;
		    while(coord_values[PDLN_LON][i] >= 360)
		        coord_values[PDLN_LON][i] -= 360;
		    while(coord_values[PDLN_LON][i] < 0)
		        coord_values[PDLN_LON][i] += 360;
		
		if(coord_values[PDLN_LON][i] < min_lon) min_lon = coord_values[PDLN_LON][i];
        if(coord_values[PDLN_LON][i] > max_lon) max_lon = coord_values[PDLN_LON][i];
        if(coord_values[PDLN_LAT][i] < min_lat) min_lat = coord_values[PDLN_LAT][i];
		if(coord_values[PDLN_LAT][i] > max_lat) max_lat = coord_values[PDLN_LAT][i];
		
		}
		delete_redundent_points(coord_values[PDLN_LON], coord_values[PDLN_LAT], num_points);
	}
	
	for(unsigned int i = 0; i < sizeof(real_global_grid_name)/64; i ++){
		if(strcmp(grid_name,real_grid_name[i])==0){
			min_lon = 0;
			max_lon = 360;
			min_lat = -90;
			max_lat = 90;
			is_cyclic = true;
			break;
		}
	}	
	if(strcmp(grid_name,"gigc_H2D_grid_temp@GIGC.nc")==0 || 
	   strcmp(grid_name,"datamodel_GIGC_H2D_grid@gamil.nc")==0 ||		
	   strcmp(grid_name,"gigc_H2D_grid@GIGC.nc")==0){
		min_lon = 0.0;
		max_lon = 360.0;
		min_lat = -90;
		max_lat = 90;
	}else if(strcmp(grid_name,"atm_grid_at_gamil_comp@gamil.nc")==0){
		min_lon = 0.0;
		max_lon = 358.0;
		min_lat = -90;
		max_lat = 90;
	}else if(strcmp(grid_name,"wwatch3_H2D_grid_temp@wwatch3.nc")==0 ||
	   		 strcmp(grid_name,"wwatch3_H2D_grid@wwatch3.nc")==0){
	    min_lon = 0.0;
	    max_lon = 360.0;
	    min_lat = -79;
	    max_lat = 79;
		is_cyclic = true;
	}else if(strcmp(grid_name,"wrf_H2D_grid_ocn_mask@wrf.nc")==0 ||
	   		 strcmp(grid_name,"wrf_H2D_grid_ocn_no_mask@wrf.nc")==0 || 
			 strcmp(grid_name,"wrf_H2D_grid_no_mask@wrf.nc")==0){
		min_lon = 0.0;
		max_lon = 360.0;
		min_lat = 0.0;
		max_lat = 90.0;
		is_cyclic = true;
	}else if(strcmp(grid_name,"ocn_grid_at_ocn_comp@cesm@OCN.nc")==0 ||
			 strcmp(grid_name,"ice_grid_at_ice_comp@cesm@ICE.nc")==0){
	    min_lon = 0.0;
	    max_lon = 360.0;
	    min_lat = -80.0;
	    max_lat = 90.0;
		is_cyclic = true;
	}else if(strcmp(grid_name,"MITgcm_H2D_grid@mitgcm.nc")==0){
		min_lon = 0.0;
        max_lon = 360.0;
	    min_lat = 30.0;
		max_lat = 90.0;
	}else if(strcmp(grid_name,"pom_grid_via_local@pom@POM_TOP_d01.nc")==0 || 
			 strcmp(grid_name,"wrf_grid_via_local@wrf@WRF_TOP_d01.nc")==0){
		min_lon = 30.0;
		max_lon = 100.0;
		min_lat = -15.0;
		max_lat = 30.0;
		is_cyclic = false;
	}else if(strcmp(grid_name,"cpl_ocn_grid@cpl.nc")==0 ||
			 strcmp(grid_name,"sis_ice_grid@sis.nc")==0 || 
			 strcmp(grid_name,"mom_ocn_grid@mom.nc")==0 || 
			 strcmp(grid_name,"mom_ocn_grid_temp@mom.nc")==0 || 
			 strcmp(grid_name,"sis_ice_grid_temp@sis.nc")==0 || 
			 strcmp(grid_name,"cpl_ice_grid@cpl.nc")==0){
        min_lon = 0.0;
        max_lon = 360.0;
        min_lat = -82.0;
        max_lat = 90.0;
        is_cyclic = true;
	}else if(strcmp(grid_name,"cfs_ocn_H2D_grid@cfs_ocn.nc")==0){
		min_lon = 0.0;
		max_lon = 360.0;
		min_lat = -81.0;
		max_lat = 90.0;
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
}

TEST_F(FullProcess, RealGridTest) {
    MPI_Barrier(MPI_COMM_WORLD);
    comm = MPI_COMM_WORLD;

    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
	
	printf("processing real grid: %s\n", my_argv[1]);
    prepare_real_grid(my_argv[1]);
    
    Component* comp;
    comp = new Component(0);
    comp->register_grid(new Grid(1));
	int ret = comp->generate_delaunay_trianglulation(1, true);
    EXPECT_EQ(ret, 0);
    delete comp;
/*
	if (mpi_rank == 0) {
		FILE *fp;
		char cmd[] = "md5sum log/global_triangles_* | awk -F\" \" '{print $1}'";
		char md5[2][64];
		memset(md5[0], 0, 64);
		memset(md5[1], 0, 64);
		fp = popen(cmd, "r");
		fgets(md5[0], 64, fp);
        fprintf(test_log,"%s", my_argv[1]);
        fflush(test_log);
        char mv[] = "mv log/global_triangles_* log/check/global_triangle_Real_";
        char tep[100];
        sprintf(tep,"%s%s",mv,my_argv[1]);
        fp = popen(tep, "r");
	}
*/
}
