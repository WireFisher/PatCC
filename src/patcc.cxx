#include "patcc.h"
#include "common_utils.h"
#include "projection.h"
#include "timer.h"
#include <cstdio>
#include <sys/time.h>
#include <omp.h>
#include <cstring>

#define PDLN_DEFAULT_MIN_NUM_POINTS (100)

long time_proc_mgt = 0;
long time_pretreat = 0;
long time_decomose = 0;
long time_expand = 0;
long time_local_tri = 0;
long time_consisty_check = 0;
long time_total = 0;

Grid* Patcc::search_grid_by_id(int id)
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

int Grid::generate_delaunay_trianglulation(Processing_resource *proc_resource, Grid_info grid_info)
{
    delaunay_triangulation = new Delaunay_grid_decomposition(grid_info, proc_resource, PDLN_DEFAULT_MIN_NUM_POINTS);

    timeval start, end;
    MPI_Barrier(proc_resource->get_mpi_comm());
    log(LOG_INFO, "decomposing grid\n");
    gettimeofday(&start, NULL);

    if(delaunay_triangulation->generate_grid_decomposition()) {
        return -1;
    }

    gettimeofday(&end, NULL);
    time_decomose += (end.tv_sec - start.tv_sec) * 1000000 + (end.tv_usec - start.tv_usec);
    //delaunay_triangulation->plot_grid_decomposition("log/grid_decomp_info.png");

    log(LOG_INFO, "generating trianglulation\n");
    gettimeofday(&start, NULL);
    int ret = delaunay_triangulation->generate_trianglulation_for_local_decomp();

    gettimeofday(&end, NULL);
    int all_ret = 0;
    MPI_Allreduce(&ret, &all_ret, 1, MPI_UNSIGNED, MPI_LOR, proc_resource->get_mpi_comm());
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


Patcc::Patcc(int id): component_id(id)
{
    proc_resource = NULL;
}

/* Debugging staff */
double global_p_lon[4];
double global_p_lat[4];
#define swap(a, b) {tmp = a; a = b; b = tmp;}

#define PAT_GVPOINT_DENSITY  (1)
#define PAT_INSERT_EXPAND_RATIO (0.01)
#define PAT_DUP_USER_INPUT (true)
void Patcc::grid_preprocessing(int grid_id)
{
    double min_lon, max_lon, min_lat, max_lat;
    double **user_coord_values;
    double *coord_values[2];
    int num_points;
    bool is_cyclic;

    grid_info_mgr->get_grid_boundry(grid_id, &min_lon, &max_lon, &min_lat, &max_lat);
    user_coord_values = grid_info_mgr->get_grid_coord_values(grid_id);
    num_points = grid_info_mgr->get_grid_num_points(grid_id);
    is_cyclic = grid_info_mgr->is_grid_cyclic(grid_id);

    if (PAT_DUP_USER_INPUT) {
        coord_values[0] = new double[num_points];
        memcpy(coord_values[0], user_coord_values[0], sizeof(double)*num_points);
        coord_values[1] = new double[num_points];
        memcpy(coord_values[1], user_coord_values[1], sizeof(double)*num_points);
    } else {
        coord_values[0] = user_coord_values[0];
        coord_values[1] = user_coord_values[1];
    }

    int total_threads = omp_get_max_threads();

    DISABLING_POINTS_METHOD mask_method;
    int num;
    void* data;
    grid_info_mgr->get_disabled_points_info(grid_id, &mask_method, &num, &data);

    bool* mask = NULL;
    if (mask_method != NO_DISABLED_POINTS)
        mask = new bool[num_points];

    if (mask_method == DISABLE_POINTS_BY_INDEX) {
        memset(mask, 1, num_points);
        int *disabled_idx = (int*) data;
        for (int i = 0; i < num; i++)
            mask[disabled_idx[i]] = false;
    }

    bool do_normalize = is_cyclic;
    bool do_monotone  = !is_cyclic && min_lon > max_lon;
    bool do_spole_processing = float_eq(min_lat, -90.0);
    bool do_npole_processing = float_eq(max_lat,  90.0);
    bool do_disabled_point_making = mask_method == DISABLE_POINTS_BY_RANGE;

    double split_line = (min_lon + max_lon) * 0.5;

    double min_lat_except_pole_public = 1e10;
    double max_lat_except_pole_public = -1e10;
    double minX_public = 1e10;
    double maxX_public = -1e10;
    double minY_public = 1e10;
    double maxY_public = -1e10;

    double* all_min_lats = new double[total_threads];
    double* all_max_lats = new double[total_threads];
    double* all_minX = new double[total_threads];
    double* all_maxX = new double[total_threads];
    double* all_minY = new double[total_threads];
    double* all_maxY = new double[total_threads];

    for (int i = 0; i < total_threads; i++) {
        all_min_lats[i] = 1e10;
        all_max_lats[i] = -1e10;
        all_minX[i] = 1e10;
        all_maxX[i] = -1e10;
        all_minY[i] = 1e10;
        all_maxY[i] = -1e10;
    }

    #pragma omp parallel for
    for (int k = 0; k < total_threads; k++) {
        int local_start = k * (num_points / total_threads);
        int local_num   = k==total_threads-1 ? num_points/total_threads+num_points%total_threads : num_points / total_threads;

        double min_lat_except_pole = 1e10;
        double max_lat_except_pole = -1e10;
        double minX = 1e10;
        double maxX = -1e10;
        double minY = 1e10;
        double maxY = -1e10;

        for(int i = local_start; i < local_start+local_num; i++) {
            if (do_normalize) {
                while(coord_values[PDLN_LON][i] >= 360) coord_values[PDLN_LON][i] -= 360;
                while(coord_values[PDLN_LON][i] < 0) coord_values[PDLN_LON][i] += 360;
            }

            if (do_monotone) {
                if(coord_values[PDLN_LON][i] > split_line) coord_values[PDLN_LON][i] -= 360;
            }

            if (do_spole_processing) {
                if(float_eq(coord_values[PDLN_LAT][i], -90.0)) {
                    #pragma omp critical
                    shifted_spoles_index.push_back(i);
                } else if(min_lat_except_pole > coord_values[PDLN_LAT][i])
                    min_lat_except_pole = coord_values[PDLN_LAT][i];
            }

            if (do_npole_processing) {
                if(float_eq(coord_values[PDLN_LAT][i], 90.0)) {
                    #pragma omp critical
                    shifted_npoles_index.push_back(i);
                } else if(max_lat_except_pole < coord_values[PDLN_LAT][i])
                    max_lat_except_pole = coord_values[PDLN_LAT][i];
            }

            if (do_disabled_point_making) {
                mask[i] = true;
                for (int j = 0; j < num; j++) {
                    double *disabled_circle = (double*) data;
                    if (point_in_circle(coord_values[0][i], coord_values[1][i], &disabled_circle[j*3])) {
                        mask[i] = false;
                        break;
                    }
                }
            }

            if (coord_values[PDLN_LON][i] < minX) minX = coord_values[PDLN_LON][i];
            if (coord_values[PDLN_LON][i] > maxX) maxX = coord_values[PDLN_LON][i];
            if (coord_values[PDLN_LAT][i] < minY) minY = coord_values[PDLN_LAT][i];
            if (coord_values[PDLN_LAT][i] > maxY) maxY = coord_values[PDLN_LAT][i];
        }

        all_min_lats[k] = min_lat_except_pole;
        all_max_lats[k] = max_lat_except_pole;
        all_minX[k] = minX;
        all_maxX[k] = maxX;
        all_minY[k] = minY;
        all_maxY[k] = maxY;
    }

    for (int i = 0; i < total_threads; i++) {
        if (min_lat_except_pole_public > all_min_lats[i]) min_lat_except_pole_public = all_min_lats[i];
        if (max_lat_except_pole_public < all_max_lats[i]) max_lat_except_pole_public = all_max_lats[i];
        if (minX_public > all_minX[i]) minX_public = all_minX[i];
        if (maxX_public < all_maxX[i]) maxX_public = all_maxX[i];
        if (minY_public > all_minY[i]) minY_public = all_minY[i];
        if (maxY_public < all_maxY[i]) maxY_public = all_maxY[i];
    }

    if(is_cyclic) {
        min_lon = 0;
        max_lon = 360;
    } else if(min_lon > max_lon) {
        /* deal with non-monotonic grid */
        PDASSERT(min_lon >= 0 && min_lon <= 360);
        PDASSERT(max_lon >= 0 && max_lon <= 360);
        min_lon -= 360;
    }

    if(do_spole_processing && shifted_spoles_index.size() != 1) {
        double shifting_lat = (-90.0 + min_lat_except_pole_public) * 0.5;

        for(unsigned i = 0; i < shifted_spoles_index.size(); i++)
            coord_values[PDLN_LAT][shifted_spoles_index[i]] = shifting_lat;
    }

    if(do_npole_processing && shifted_npoles_index.size() != 1) {
        double shifting_lat = (90.0 + max_lat_except_pole_public) * 0.5;

        for(unsigned i = 0; i < shifted_npoles_index.size(); i++)
            coord_values[PDLN_LAT][shifted_npoles_index[i]] = shifting_lat;
    }

    bool do_fence_point_inserting = !float_eq(min_lat, -90) || !float_eq(max_lat, 90) || !is_cyclic;
    bool do_virtual_pole_inserting = (do_spole_processing && shifted_spoles_index.size() != 1) ||
                                     (do_npole_processing && shifted_npoles_index.size() != 1);

    double* extended_coord[2];
    bool*   extended_mask = NULL;

    int num_vpoles, num_current;
    if (!do_fence_point_inserting && !do_virtual_pole_inserting) {
        extended_coord[0] = coord_values[0];
        extended_coord[1] = coord_values[1];
        extended_mask = mask;
        num_vpoles = 0;
        num_current = num_points;
    } else {
        /* do_fence_point_inserting */
        double dx   = maxX_public - minX_public;
        double dy   = maxY_public - minY_public;
        double dMax = std::max(dx, dy);

        double v_minx = minX_public-dMax*PAT_INSERT_EXPAND_RATIO;
        double v_maxx = maxX_public+dMax*PAT_INSERT_EXPAND_RATIO;
        double v_miny = minY_public-dMax*PAT_INSERT_EXPAND_RATIO;
        double v_maxy = maxY_public+dMax*PAT_INSERT_EXPAND_RATIO;

        double r_minx = is_cyclic ? std::max(v_minx, 0.) : v_minx;
        double r_maxx = is_cyclic ? std::min(v_maxx, 360.) : v_maxx;
        double r_miny = std::max(v_miny, -89.9999);
        double r_maxy = std::min(v_maxy,  89.9999);

        /* x * y = num_points, x : y = dx : dy */
        unsigned num_x = (unsigned)sqrt(num_points * dx / dy);
        unsigned num_y = num_x * dy / dx;
        num_x /= PAT_GVPOINT_DENSITY;
        num_y /= PAT_GVPOINT_DENSITY;

        unsigned num_new_points = 0;
        if (!float_eq(min_lat, -90))
            num_new_points += num_x;
        if (!float_eq(max_lat, 90))
            num_new_points += num_x;
        if (!is_cyclic)
            num_new_points += num_y*2;
        if (do_virtual_pole_inserting)
            num_new_points += 2;

        extended_coord[0] = new double[num_points + num_new_points];
        extended_coord[1] = new double[num_points + num_new_points];
        if (mask)
            extended_mask = new bool[num_points + num_new_points];

        /* Firstly, store all original points */
        memcpy(extended_coord[PDLN_LON], coord_values[PDLN_LON], num_points*sizeof(double));
        memcpy(extended_coord[PDLN_LAT], coord_values[PDLN_LAT], num_points*sizeof(double));
        if (mask)
            memcpy(extended_mask, mask, num_points*sizeof(bool));

        num_current = num_points;
        /* Then, virtual poles follow */
        if(do_spole_processing && shifted_spoles_index.size() != 1) {
            extended_coord[PDLN_LON][num_current] = 0;
            extended_coord[PDLN_LAT][num_current] = -90;
            if (mask)
                extended_mask[num_current] = true;
            num_current++;
        }

        if(do_npole_processing && shifted_npoles_index.size() != 1) {
            extended_coord[PDLN_LON][num_current] = 0;
            extended_coord[PDLN_LAT][num_current] = 90;
            if (mask)
                extended_mask[num_current] = true;
            num_current++;
        }
        num_vpoles = num_current - num_points;

        /* At last, fence points follow */
        if(!float_eq(min_lat, -90)) {
            /* use (num_x+1) and (i+1) to avoid generating end points, because
             * end points will overlap with those of other inserting line. */
            for(unsigned i = 0; i < num_x; i++) {
                extended_coord[PDLN_LON][num_current] = r_minx+(r_maxx-r_minx)/(num_x+1)*(i+1);
                extended_coord[PDLN_LAT][num_current++] = r_miny;
            }
            if(min_lat > r_miny) min_lat = r_miny;
        }

        if(!float_eq(max_lat, 90)) {
            for(unsigned i = 0; i < num_x; i++) {
                extended_coord[PDLN_LON][num_current] = r_minx+(r_maxx-r_minx)/(num_x+1)*(i+1);
                extended_coord[PDLN_LAT][num_current++] = r_maxy;
            }
            if(max_lat < r_maxy) max_lat = r_maxy;
        }

        if(!is_cyclic) {
            for(unsigned i = 0; i < num_y; i++) {
                extended_coord[PDLN_LON][num_current] = r_minx;
                extended_coord[PDLN_LAT][num_current++] = r_miny+(r_maxy-r_miny)/(num_y+1)*(i+1);
            }
            if(min_lon > r_minx) min_lon = r_minx;

            for(unsigned i = 0; i < num_y; i++) {
                extended_coord[PDLN_LON][num_current] = r_maxx;
                extended_coord[PDLN_LAT][num_current++] = r_miny+(r_maxy-r_miny)/(num_y+1)*(i+1);
            }
            if(max_lon < r_maxx) max_lon = r_maxx;
        }

        if (mask)
            memset(&extended_mask[num_points], 1, num_current - num_points);

        PDASSERT(num_current - num_points <= num_new_points && num_current - num_points >= num_new_points - 2);
    }

    grid_info.coord_values[PDLN_LON] = extended_coord[PDLN_LON];
    grid_info.coord_values[PDLN_LAT] = extended_coord[PDLN_LAT];
    grid_info.mask = extended_mask;
    grid_info.num_total_points = num_current;
    grid_info.num_vitual_poles = num_vpoles;
    grid_info.num_fence_points = num_current - num_points - num_vpoles;
    grid_info.boundary.min_lon = min_lon;
    grid_info.boundary.max_lon = max_lon;
    grid_info.boundary.min_lat = min_lat;
    grid_info.boundary.max_lat = max_lat;
    grid_info.is_cyclic = is_cyclic;
}


Patcc::~Patcc()
{
    delete proc_resource;
    for(unsigned i = 0; i < grids.size(); i ++)
        delete grids[i];
}


int Patcc::generate_delaunay_trianglulation(int grid_id, bool sort)
{
    timeval start, end;

    Grid *operating_grid = search_grid_by_id(grid_id);
    if (operating_grid == NULL)
        return -1;

    //PDLN_Timer timer;
    //timer.tick();
    log(LOG_INFO, "initializing processing resource manager\n");
    MPI_Barrier(process_thread_mgr->get_mpi_comm());
    gettimeofday(&start, NULL);
    if(proc_resource == NULL)
        proc_resource = new Processing_resource();
    gettimeofday(&end, NULL);
    //double time = timer.tick();
    time_proc_mgt += (end.tv_sec - start.tv_sec) * 1000000 + (end.tv_usec - start.tv_usec);

    //proc_resource->print_all_nodes_info();

    log(LOG_INFO, "preprocessing grid\n");
    gettimeofday(&start, NULL);
    grid_preprocessing(grid_id);
    gettimeofday(&end, NULL);
    time_pretreat += (end.tv_sec - start.tv_sec) * 1000000 + (end.tv_usec - start.tv_usec);

    gettimeofday(&start, NULL);
    if(operating_grid->generate_delaunay_trianglulation(proc_resource, grid_info)) {
        log(LOG_ERROR, "failed\n");
        return -1;
    }

#ifdef OPENCV
    operating_grid->plot_triangles_into_file();
#endif

    log(LOG_INFO, "collecting results\n");
    operating_grid->merge_all_triangles(sort);
    gettimeofday(&end, NULL);
    time_total = (end.tv_sec - start.tv_sec) * 1000000 + (end.tv_usec - start.tv_usec);

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
    log(LOG_INFO, "%ld\n", g_time_proc_mgt);
    log(LOG_INFO, "%ld\n", g_time_pretreat);
    log(LOG_INFO, "%ld\n", g_time_decomose);
    log(LOG_INFO, "%ld\n", g_time_expand);
    log(LOG_INFO, "%ld\n", g_time_local_tri);
    log(LOG_INFO, "%ld\n", g_time_consisty_check);
    log(LOG_INFO, "%ld\n", g_time_total);
#endif
    return 0;
}
