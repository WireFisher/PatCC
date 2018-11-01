#include "component.h"
#include "common_utils.h"
#include <cstdio>
#include <sys/time.h>
#include "timer.h"
#include <omp.h>

#define PDLN_DEFAULT_MIN_NUM_POINTS (100)

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
    delaunay_triangulation = new Delaunay_grid_decomposition(grid_id, proc_resource, PDLN_DEFAULT_MIN_NUM_POINTS);

    timeval start, end;
    MPI_Barrier(proc_resource->get_mpi_comm());
    gettimeofday(&start, NULL);
    if(delaunay_triangulation->generate_grid_decomposition()) {
        return -1;
    }

    gettimeofday(&end, NULL);
#ifdef TIME_PERF
    printf("[ - ] Grid Decomposition: %ld ms\n", (end.tv_sec - start.tv_sec) * 1000000 + (end.tv_usec - start.tv_usec));
    time_decomose += (end.tv_sec - start.tv_sec) * 1000000 + (end.tv_usec - start.tv_usec);
#endif
    //delaunay_triangulation->plot_grid_decomposition("log/grid_decomp_info.png");

    gettimeofday(&start, NULL);
    int ret = delaunay_triangulation->generate_trianglulation_for_local_decomp();

    gettimeofday(&end, NULL);
    int all_ret = 0;
    MPI_Allreduce(&ret, &all_ret, 1, MPI_UNSIGNED, MPI_LOR, proc_resource->get_mpi_comm());
#ifdef TIME_PERF
    printf("[ - ] All Trianglulation: %ld ms\n", (end.tv_sec - start.tv_sec) * 1000000 + (end.tv_usec - start.tv_usec));
#endif
    if(all_ret == 1)
        return -1;

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

    int total_threads = omp_get_max_threads();

    if(is_cyclic) {
        if(!(min_lon >= 0 && min_lon <=360 && max_lon >=0 && max_lon <= 360)) {
            #pragma omp parallel for
            for (int k = 0; k < total_threads; k++) {
                int local_start = k * (num_points / total_threads);
                int local_num   = k==total_threads-1 ? num_points/total_threads+num_points%total_threads : num_points / total_threads;

                for(int i = local_start; i < local_start+local_num; i++) {
                    while(coord_values[PDLN_LON][i] >= 360) coord_values[PDLN_LON][i] -= 360;
                    while(coord_values[PDLN_LON][i] < 0) coord_values[PDLN_LON][i] += 360;
                }
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
        #pragma omp parallel for
        for (int k = 0; k < total_threads; k++) {
            int local_start = k * (num_points / total_threads);
            int local_num   = k==total_threads-1 ? num_points/total_threads+num_points%total_threads : num_points / total_threads;

            for(int i = local_start; i < local_start+local_num; i++)
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

    //PDLN_Timer timer;
    //timer.tick();
    MPI_Barrier(process_thread_mgr->get_mpi_comm());
    gettimeofday(&start, NULL);
    if(proc_resource == NULL)
        proc_resource = new Processing_resource();
    gettimeofday(&end, NULL);
    //double time = timer.tick();
#ifdef TIME_PERF
    printf("[ - ] Procs Resource MGR: %ld us\n", (end.tv_sec - start.tv_sec) * 1000000 + (end.tv_usec - start.tv_usec));
    time_proc_mgt += (end.tv_sec - start.tv_sec) * 1000000 + (end.tv_usec - start.tv_usec);
#endif

    //proc_resource->print_all_nodes_info();

    gettimeofday(&start, NULL);
    grid_pretreatment(grid_id);
    gettimeofday(&end, NULL);
#ifdef TIME_PERF
    printf("[ - ] Gri Pre-treatment: %ld ms\n", (end.tv_sec - start.tv_sec) * 1000000 + (end.tv_usec - start.tv_usec));
    time_pretreat += (end.tv_sec - start.tv_sec) * 1000000 + (end.tv_usec - start.tv_usec);
#endif

    gettimeofday(&start, NULL);
    if(operating_grid->generate_delaunay_trianglulation(this->proc_resource))
        return -1;
#ifdef OPENCV
    operating_grid->plot_triangles_into_file();
#endif
    operating_grid->merge_all_triangles(sort);
    gettimeofday(&end, NULL);
#ifdef TIME_PERF
    printf("[ - ] Total Time Elapsed: %ld ms\n", (end.tv_sec - start.tv_sec) * 1000000 + (end.tv_usec - start.tv_usec));
    time_total = (end.tv_sec - start.tv_sec) * 1000000 + (end.tv_usec - start.tv_usec);
#endif

#ifdef TIME_PERF

    long g_time_proc_mgt, g_time_pretreat, g_time_decomose, g_time_expand, g_time_local_tri, g_time_consisty_check, g_time_total;
    MPI_Comm comm = process_thread_mgr->get_mpi_comm();
    MPI_Reduce(&time_proc_mgt, &g_time_proc_mgt, 1, MPI_LONG, MPI_MAX, 0, comm);
    MPI_Reduce(&time_pretreat, &g_time_pretreat, 1, MPI_LONG, MPI_MAX, 0, comm);
    MPI_Reduce(&time_decomose, &g_time_decomose, 1, MPI_LONG, MPI_MAX, 0, comm);
    MPI_Reduce(&time_expand, &g_time_expand, 1, MPI_LONG, MPI_MAX, 0, comm);
    MPI_Reduce(&time_local_tri, &g_time_local_tri, 1, MPI_LONG, MPI_MAX, 0, comm);
    MPI_Reduce(&time_consisty_check, &g_time_consisty_check, 1, MPI_LONG, MPI_MAX, 0, comm);
    MPI_Reduce(&time_total, &g_time_total, 1, MPI_LONG, MPI_MAX, 0, comm);
    printf("%ld\n", g_time_proc_mgt);
    printf("%ld\n", g_time_pretreat);
    printf("%ld\n", g_time_decomose);
    printf("%ld\n", g_time_expand);
    printf("%ld\n", g_time_local_tri);
    printf("%ld\n", g_time_consisty_check);
    printf("%ld\n", g_time_total);
#endif
    return 0;
}
