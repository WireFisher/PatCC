#include "component.h"
#include "common_utils.h"
#include <cstdio>
#include <sys/time.h>
#include "timer.h"
#define DEFAULT_MIN_NUM_POINTS 100
#define MULTIPLICATION_COEFFICIENT 2

long time_proc_mgt = 0;
long time_pretreat = 0;
long time_decomose = 0;
long time_expand = 0;
long time_local_tri = 0;
long time_consisty_check = 0;
long time_total = 0;

Grid* Component::search_grid_by_id(int id)
{
    for(unsigned int i = 0; i < this->grids.size(); i++)
        if(grids[i]->get_grid_id() == id)
            return grids[i];

    return NULL;
}

Grid::~Grid()
{
    delete delaunay_triangulation;
}

int Grid::generate_delaunay_trianglulation(Processing_resource *proc_resource)
{
    int min_num_points_per_chunk, ret;
    bool do_sequentially = false;
    if(delaunay_triangulation)
        return 0;

    min_num_points_per_chunk = DEFAULT_MIN_NUM_POINTS;
    do {
        min_num_points_per_chunk *= MULTIPLICATION_COEFFICIENT;
        if(delaunay_triangulation) {
            delete delaunay_triangulation;
            delaunay_triangulation = NULL;
        }
        delaunay_triangulation = new Delaunay_grid_decomposition(this->grid_id, proc_resource, min_num_points_per_chunk);

        timeval start, end;
        MPI_Barrier(proc_resource->get_mpi_comm());
        gettimeofday(&start, NULL);
        if(delaunay_triangulation->generate_grid_decomposition()) {
            do_sequentially = true;
            break;
        }

        MPI_Barrier(proc_resource->get_mpi_comm());
        gettimeofday(&end, NULL);
#ifdef TIME_PERF
        printf("[ - ] Grid Decomposition: %ld ms\n", ((end.tv_sec - start.tv_sec) * 1000000 + (end.tv_usec - start.tv_usec)) / 1000);
        time_decomose = ((end.tv_sec - start.tv_sec) * 1000000 + (end.tv_usec - start.tv_usec)) / 1000;
#endif
        //delaunay_triangulation->plot_grid_decomposition("log/grid_decomp_info.png");

        gettimeofday(&start, NULL);
        ret = delaunay_triangulation->generate_trianglulation_for_local_decomp();
        MPI_Barrier(proc_resource->get_mpi_comm());
        gettimeofday(&end, NULL);
        int all_ret = 0;
        MPI_Allreduce(&ret, &all_ret, 1, MPI_UNSIGNED, MPI_LOR, proc_resource->get_mpi_comm());
        ret = all_ret;
#ifdef TIME_PERF
        printf("[ - ] All Trianglulation: %ld ms\n", ((end.tv_sec - start.tv_sec) * 1000000 + (end.tv_usec - start.tv_usec)) / 1000);
#endif
        /* Return Values: 0 - success
         *                1 - fail, normal decomp's expanded_boundry exceeded too large (expanding fail)
         *                2 - fail, polar  decomp's expanded_boundry exceeded threshold */
        if(ret == 1)
            return -1;
        if(ret == 2) {
            do_sequentially = true;
            break;
        }
    } while(ret);

    //mpi allgather do_sequentially
    if(do_sequentially) { 
        if(delaunay_triangulation) {
            delete delaunay_triangulation;
            delaunay_triangulation = NULL;
        }
        this->delaunay_triangulation = new Delaunay_grid_decomposition(this->grid_id, NULL, 0);
        assert(false);
    }
    return 0;
}


#ifdef OPENCV
void Grid::plot_triangles_into_file()
{
    this->delaunay_triangulation->plot_local_triangles("log/chunk");
}
#endif


void Grid::merge_all_triangles(bool sort)
{
    delaunay_triangulation->merge_all_triangles(sort);
}


Component::Component(int id): component_id(id)
{
    proc_resource = NULL;
}

/* Debugging staff */
double global_p_lon[4];
double global_p_lat[4];
#define swap(a, b) {tmp = a; a = b; b = tmp;}
void Component::grid_pretreatment(int grid_id)
{
    double min_lon, max_lon, min_lat, max_lat;
    double **coord_values;
    int num_points;
    bool is_cyclic;

    grid_info_mgr->get_grid_boundry(grid_id, &min_lon, &max_lon, &min_lat, &max_lat);
    coord_values = grid_info_mgr->get_grid_coord_values(grid_id);
    num_points = grid_info_mgr->get_grid_num_points(grid_id);
    is_cyclic = grid_info_mgr->is_grid_cyclic(grid_id);

    //int count = 0;
    //for(int i = 0; i < num_points; i++) {//(0.000000, -90.000000), (315.000000, -88.827445), (360.000000, -89.170820)
    //    //if(coord_values[PDLN_LON][i] < 348.5 && coord_values[PDLN_LON][i] > 344.7 && coord_values[PDLN_LAT][i] > 3.2 && coord_values[PDLN_LAT][i] < 4.1) {
    //    if(fabs(coord_values[PDLN_LON][i] - 0) < 1e-5 && fabs(coord_values[PDLN_LAT][i] - -90.000000) < 1e-5) {
    //        global_p_lon[0] = coord_values[PDLN_LON][i];
    //        global_p_lat[0] = coord_values[PDLN_LAT][i];
    //        printf("found1\n");
    //    }
    //    if(fabs(coord_values[pdln_lon][i] - 315.000000) < 1e-5 && fabs(coord_values[pdln_lat][i] - -88.827445) < 1e-5) {
    //        global_p_lon[1] = coord_values[pdln_lon][i];
    //        global_p_lat[1] = coord_values[pdln_lat][i];
    //        printf("found2\n");
    //    }
    //    if(fabs(coord_values[pdln_lon][i] - 360.000000) < 1e-5 && fabs(coord_values[pdln_lat][i] - -89.170820) < 1e-5) {
    //        global_p_lon[2] = coord_values[pdln_lon][i];
    //        global_p_lat[2] = coord_values[pdln_lat][i];
    //        printf("found3\n");
    //    }
    //}
            //printf("coord: %lf, %lf\n", coord_values[pdln_lon][i], coord_values[pdln_lat][i]);
            //global_p_lon[count] = coord_values[pdln_lon][i];
            //global_p_lat[count++] = coord_values[pdln_lat][i];
            //PDASSERT(count <= 4);
    //double tmp;
    //swap(global_p_lon[2], global_p_lon[3]);
    //swap(global_p_lat[2], global_p_lat[3]);
    if(is_cyclic) {
        if(!(min_lon >= 0 && min_lon <=360 && max_lon >=0 && max_lon <= 360)) {
            for(int i = 0; i < num_points; i++) {
                while(coord_values[PDLN_LON][i] >= 360) coord_values[PDLN_LON][i] -= 360;
                while(coord_values[PDLN_LON][i] < 0) coord_values[PDLN_LON][i] += 360;
            }
            min_lon = 0;
            max_lon = 360;
            grid_info_mgr->set_grid_boundry(grid_id, min_lon, max_lon, min_lat, max_lat);
        }
    }
    /* deal with non-monotonic grid */
    else if(min_lon > max_lon) {
        PDASSERT(min_lon >= 0 && min_lon <= 360);
        PDASSERT(max_lon >= 0 && max_lon <= 360);

        double split_line = (min_lon + max_lon) * 0.5;
        for(int i = 0; i < num_points; i++) {
            if(coord_values[PDLN_LON][i] > split_line) coord_values[PDLN_LON][i] -= 360;
        }
        min_lon -= 360;
        grid_info_mgr->set_grid_boundry(grid_id, min_lon, max_lon, min_lat, max_lat);
    }
}


Component::~Component()
{
    delete proc_resource;
    for(unsigned i = 0; i < grids.size(); i ++)
        delete grids[i];
}


int Component::generate_delaunay_trianglulation(int grid_id, bool sort)
{
    timeval start, end;

    Grid *operating_grid = this->search_grid_by_id(grid_id);
    if (operating_grid == NULL)
        return -1;

    MPI_Barrier(process_thread_mgr->get_mpi_comm());
    gettimeofday(&start, NULL);
    if(proc_resource == NULL)
        proc_resource = new Processing_resource();
    MPI_Barrier(proc_resource->get_mpi_comm());
    gettimeofday(&end, NULL);
#ifdef TIME_PERF
    printf("[ - ] Procs Resource MGR: %ld us\n", (end.tv_sec - start.tv_sec) * 1000000 + (end.tv_usec - start.tv_usec));
    time_proc_mgt = (end.tv_sec - start.tv_sec) * 1000000 + (end.tv_usec - start.tv_usec);
#endif

    //proc_resource->print_all_nodes_info();

    gettimeofday(&start, NULL);
    grid_pretreatment(grid_id);
    MPI_Barrier(proc_resource->get_mpi_comm());
    gettimeofday(&end, NULL);
#ifdef TIME_PERF
    printf("[ - ] Grid Pre-treatment: %ld ms\n", ((end.tv_sec - start.tv_sec) * 1000000 + (end.tv_usec - start.tv_usec)) / 1000);
    time_pretreat += ((end.tv_sec - start.tv_sec) * 1000000 + (end.tv_usec - start.tv_usec)) / 1000;
#endif

    gettimeofday(&start, NULL);
    if(operating_grid->generate_delaunay_trianglulation(this->proc_resource))
        return -1;
#ifdef OPENCV
    operating_grid->plot_triangles_into_file();
#endif
    operating_grid->merge_all_triangles(sort);
    MPI_Barrier(proc_resource->get_mpi_comm());
    gettimeofday(&end, NULL);
#ifdef TIME_PERF
    printf("[ - ] Total Time Elapsed: %ld ms\n", ((end.tv_sec - start.tv_sec) * 1000000 + (end.tv_usec - start.tv_usec)) / 1000);
    time_total = ((end.tv_sec - start.tv_sec) * 1000000 + (end.tv_usec - start.tv_usec)) / 1000;
#endif

#ifdef TIME_PERF
    printf("%ld\n", time_proc_mgt);
    printf("%ld\n", time_pretreat);
    printf("%ld\n", time_decomose);
    printf("%ld\n", time_expand);
    printf("%ld\n", time_local_tri);
    printf("%ld\n", time_consisty_check);
    printf("%ld\n", time_total);
#endif
    return 0;
}
