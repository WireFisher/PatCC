/***************************************************************
  *  Copyright (c) 2013, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Haoyu Yang,
  *  supervised by Dr. Li Liu. If you have any problem,
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/

#include "delaunay_grid_decomposition_mgt.h"
#include <cstdio>
#include <cstddef>
#include <cstring>
#include <cassert>
#include <algorithm>
#include <cmath>
#include <vector>
#include <tr1/unordered_map>
#include <sys/time.h>
#include "ccpl_utils.h"
#include "netcdf_utils.h"
#include "opencv_utils.h"

#define PDLN_DEFAULT_EXPANGDING_RATIO (0.2)
#define PDLN_EXPECTED_EXPANDING_LOOP_TIMES (3)

#define PDLN_SPOLAR_MAX_LAT (-19.47)
#define PDLN_NPOLAR_MIN_LAT (19.47)

#define PDLN_MAX_ITER_COUNT (10)

#define PDLN_TOLERABLE_ERROR (0.0001)

#define PDLN_DECOMPOSE_COMMON_MODE (0)
#define PDLN_DECOMPOSE_SPOLAR_MODE (1)
#define PDLN_DECOMPOSE_NPOLAR_MODE (2)

#define PDLN_NODE_TYPE_COMMON PDLN_DECOMPOSE_COMMON_MODE
#define PDLN_NODE_TYPE_SPOLAR PDLN_DECOMPOSE_SPOLAR_MODE
#define PDLN_NODE_TYPE_NPOLAR PDLN_DECOMPOSE_NPOLAR_MODE

#define PDLN_DOUBLE_INVALID_VALUE ((double)0xDEADBEEFDEADBEEF)

#define PDLN_HIGH_BOUNDRY_SHIFTING (1e-6)

#define PDLN_MAX_NUM_PROCESSING_UNITS 512

#define PDLN_UP     0
#define PDLN_LEFT   1
#define PDLN_DOWN   2
#define PDLN_RIGHT  3

bool operator == (Triangle_ID_Only t1, Triangle_ID_Only t2)
{
    if(t2.id[0] != t1.id[0] && t2.id[0] != t1.id[1] && t2.id[0] != t1.id[2])
        return false;
    if(t2.id[1] != t1.id[0] && t2.id[1] != t1.id[1] && t2.id[1] != t1.id[2])
        return false;
    if(t2.id[2] != t1.id[0] && t2.id[2] != t1.id[1] && t2.id[2] != t1.id[2])
        return false;
    return true;
}


bool Boundry::operator== (const Boundry &boundry) const
{
    return min_lat == boundry.min_lat && min_lon == boundry.min_lon && max_lat == boundry.max_lat && max_lon == boundry.max_lon;
}


bool Boundry::operator!= (const Boundry &boundry) const
{
    return !(min_lat == boundry.min_lat && min_lon == boundry.min_lon && max_lat == boundry.max_lat && max_lon == boundry.max_lon);
}


bool Boundry::operator<= (const Boundry &boundry) const
{
    return min_lat >= boundry.min_lat && min_lon >= boundry.min_lon && max_lat <= boundry.max_lat && max_lon <= boundry.max_lon;
}


Boundry& Boundry::operator* (double ratio)
{
    min_lat -= (max_lat - min_lat) * ratio * 0.5;
    max_lat += (max_lat - min_lat) * ratio * 0.5;
    min_lon -= (max_lon - min_lon) * ratio * 0.5;
    max_lon += (max_lon - min_lon) * ratio * 0.5;
    return *this;
}


void Boundry::squeeze(const Boundry *inner, double area_ratio)
{
    double area_outter = (max_lat - min_lat) * (max_lon - min_lon);
    double area_inner  = (inner->max_lat - inner->min_lat) * (inner->max_lon - inner->min_lon);
    double edge_ratio  = sqrt(area_ratio * area_outter / (area_outter + (area_ratio-1) * area_inner));
    printf("edge_ratio %lf\n", edge_ratio);
    assert(edge_ratio >= 1);

#ifdef DEBUG
    assert(*inner <= *this);
    //assert((max_lat - min_lat) * edge_ratio * (inner->min_lat - min_lat) / (max_lat - min_lat - inner->max_lat + inner->min_lat) >= 0);
    //assert((max_lat - min_lat) * edge_ratio * (max_lat - inner->max_lat) / (max_lat - min_lat - inner->max_lat + inner->min_lat) >= 0);
    //assert((max_lon - min_lon) * edge_ratio * (inner->min_lon - min_lon) / (max_lon - min_lon - inner->max_lon + inner->min_lon) >= 0);
    //assert((max_lon - min_lon) * edge_ratio * (max_lon - inner->max_lon) / (max_lon - min_lon - inner->max_lon + inner->min_lon) >= 0);
#endif
    
    printf("old [%lf, %lf], [%lf, %lf]\n", min_lon, max_lon, min_lat, max_lat);
    //double delta = (1 - edge_ratio) / edge_ratio;
    double delta = (edge_ratio - 1) / edge_ratio;
    double lef_rate = (inner->min_lon - min_lon) / (max_lon - min_lon - inner->max_lon + inner->min_lon);
    double rit_rate = (max_lon - inner->max_lon) / (max_lon - min_lon - inner->max_lon + inner->min_lon);
    double bot_rate = (inner->min_lat - min_lat) / (max_lat - min_lat - inner->max_lat + inner->min_lat);
    double top_rate = (max_lat - inner->max_lat) / (max_lat - min_lat - inner->max_lat + inner->min_lat);
    printf("delta: %lf\n", delta);
    if(fabs(max_lon - min_lon - inner->max_lon + inner->min_lon) < FLOAT_ERROR)
        lef_rate = rit_rate = 0;
    if(fabs(max_lat - min_lat - inner->max_lat + inner->min_lat) < FLOAT_ERROR)
        bot_rate = top_rate = 0;

    assert(delta >= 0);
    min_lon += (max_lon - min_lon) * delta * lef_rate;
    max_lon -= (max_lon - min_lon) * delta * rit_rate;
    min_lat += (max_lat - min_lat) * delta * bot_rate;
    max_lat -= (max_lat - min_lat) * delta * top_rate;
    printf("new [%lf, %lf], [%lf, %lf]\n", min_lon, max_lon, min_lat, max_lat);
    printf("ref [%lf, %lf], [%lf, %lf]\n", inner->min_lon, inner->max_lon, inner->min_lat, inner->max_lat);

    max(*inner);
    printf("max [%lf, %lf], [%lf, %lf]\n", min_lon, max_lon, min_lat, max_lat);
}


void Boundry::legalize(const Boundry *outer_boundry, bool is_cyclic)
{
    min_lat = std::max(min_lat, outer_boundry->min_lat);
    max_lat = std::min(max_lat, outer_boundry->max_lat);
    if(!is_cyclic) {
        min_lon = std::max(min_lon, outer_boundry->min_lon);
        max_lon = std::min(max_lon, outer_boundry->max_lon);
    }
}


void Boundry::max(const Boundry b)
{
    double min_lo = b.min_lon;
    double max_lo = b.max_lon;
    double min_la = b.min_lat;
    double max_la = b.max_lat;
    min_lon = std::min(min_lon, min_lo);
    max_lon = std::max(max_lon, max_lo);
    min_lat = std::min(min_lat, min_la);
    max_lat = std::max(max_lat, max_la);
}


void Boundry::max(double min_lo, double max_lo, double min_la, double max_la)
{
    min_lon = std::min(min_lon, min_lo);
    max_lon = std::max(max_lon, max_lo);
    min_lat = std::min(min_lat, min_la);
    max_lat = std::max(max_lat, max_la);
}

Search_tree_node::Search_tree_node(Search_tree_node *p, double *coord_value[2], int *global_index, int num_points, Boundry boundry, int type)
    : parent(p)
    , node_type(type)
    , region_id(-1)
    , kernel_boundry(NULL)
    , expanded_boundry(NULL)
    , rotated_kernel_boundry(NULL)
    , rotated_expanded_boundry(NULL)
    , real_boundry(NULL)
    , len_expanded_points_coord_buf(0)
    , num_kernel_points(num_points)
    , num_expanded_points(0)
    , num_projected_points(0)
    , midline(Midline{-1, -361.0})
    , group_intervals(NULL)
    , triangulation(NULL)
    , virtual_point_local_index(-1)
{
    children[0] = NULL;
    children[1] = NULL;
    children[2] = NULL;
    projected_coord[0] = NULL;
    projected_coord[1] = NULL;

    kernel_boundry    = new Boundry();
    expanded_boundry  = new Boundry();
    *kernel_boundry   = boundry;
    *expanded_boundry = boundry;

    points_coord[0]     = new double[num_points];
    points_coord[1]     = new double[num_points];
    points_global_index = new int[num_points];
    memcpy(points_coord[0], coord_value[0], num_points * sizeof(double));
    memcpy(points_coord[1], coord_value[1], num_points * sizeof(double));
    memcpy(points_global_index, global_index, num_points * sizeof(int));

    if(type == PDLN_NODE_TYPE_COMMON) {
        center[PDLN_LON] = (boundry.min_lon + boundry.max_lon) * 0.5;
        center[PDLN_LAT] = (boundry.min_lat + boundry.max_lat) * 0.5;
        //fix_view_point();
    }
    else if(type == PDLN_NODE_TYPE_SPOLAR) {
        center[PDLN_LON] = 0.0;
        center[PDLN_LAT] = -90.0;
    }
    else if(type == PDLN_NODE_TYPE_NPOLAR) {
        center[PDLN_LON] = 0.0;
        center[PDLN_LAT] = 90.0;
    }

    non_monotonic = kernel_boundry->min_lon > kernel_boundry->max_lon;

    expanding_scale[0] = 2;
    expanding_scale[1] = 2;
    expanding_scale[2] = 2;
    expanding_scale[3] = 2;
    num_neighbors_on_boundry[0] = 0;
    num_neighbors_on_boundry[1] = 0;
    num_neighbors_on_boundry[2] = 0;
    num_neighbors_on_boundry[3] = 0;
    edge_expanding_count[0] = 0;
    edge_expanding_count[1] = 0;
    edge_expanding_count[2] = 0;
    edge_expanding_count[3] = 0;
}


Search_tree_node::~Search_tree_node()
{
    delete kernel_boundry;
    delete expanded_boundry;
    delete[] points_coord[0];
    delete[] points_coord[1];
    delete[] points_global_index;
    delete rotated_kernel_boundry;
    delete rotated_expanded_boundry;
    for(int i = 0; i < 3; i ++)
        delete children[i];
    delete triangulation;
}


//const double fixed_view_points[2][8] = {{45, 135, 225, 315, 45, 135, 225, 315},
//                                       {-45, -45, -45, -45, 45,  45,  45,  45}};
const double fixed_view_points[2][4] = {{45, 135, 225, 315},
                                       {0, 0, 0, 0}};
void Search_tree_node::fix_view_point()
{
    unsigned current_index = 0;
    double current_distence;

    while(center[PDLN_LON] <    0) center[PDLN_LON] += 360;
    while(center[PDLN_LON] >= 360) center[PDLN_LON] -= 360;
    current_distence = (center[PDLN_LON] - fixed_view_points[PDLN_LON][0]) * (center[PDLN_LON] - fixed_view_points[PDLN_LON][0]) +
                       (center[PDLN_LAT] - fixed_view_points[PDLN_LAT][0]) * (center[PDLN_LAT] - fixed_view_points[PDLN_LAT][0]); //TODO: merge code
    for(unsigned i = 1; i < sizeof(fixed_view_points)/sizeof(double)/2; i++) {
        double dist = (center[PDLN_LON] - fixed_view_points[PDLN_LON][i]) * (center[PDLN_LON] - fixed_view_points[PDLN_LON][i]) +
                      (center[PDLN_LAT] - fixed_view_points[PDLN_LAT][i]) * (center[PDLN_LAT] - fixed_view_points[PDLN_LAT][i]);
        if (dist < current_distence) {
            current_distence = dist;
            current_index = i;
        }
    }
    center[PDLN_LON] = fixed_view_points[PDLN_LON][current_index];
    center[PDLN_LAT] = fixed_view_points[PDLN_LAT][current_index];
}


/* assumption: ids are already sorted. */
void Search_tree_node::update_region_ids(int num)
{
    region_ids.clear();
    region_ids.reserve(num);
    for(int i = 0; i < num; i ++)
        region_ids.push_back(i);
}


void Search_tree_node::update_region_ids(vector<int> proc_units_id)
{
    region_ids.clear();
    region_ids = proc_units_id;
}


inline void calculate_circle_center(double x[3], double y[3], double *center_x, double *center_y)
{
    double mid_x[2], mid_y[2];
    double k[2];

    mid_x[0] = (x[0] + x[1]) * 0.5;
    mid_y[0] = (y[0] + y[1]) * 0.5;
    mid_x[1] = (x[0] + x[2]) * 0.5;
    mid_y[1] = (y[0] + y[2]) * 0.5;
    k[0] = - (x[1] - x[0]) / (y[1] - y[0]);
    k[1] = - (x[2] - x[0]) / (y[2] - y[0]);

    *center_x = (mid_y[1] - mid_y[0] - k[1]*mid_x[1] + k[0]*mid_x[0]) / (k[0] - k[1]);
    *center_y = mid_y[0] + k[0]*(mid_y[1] - mid_y[0] - k[1]*mid_x[1] + k[1]*mid_x[0]) / (k[0] - k[1]);
}

extern double global_p_lon[4];
extern double global_p_lat[4];
#define PDLN_INSERT_VIRTUAL_POINT (true)
#define PDLN_REMOVE_UNNECESSARY_TRIANGLES (false)
void Search_tree_node::generate_local_triangulation(bool is_cyclic, int num_inserted)
{
    delete triangulation;

    //char filename[64];
    //snprintf(filename, 64, "log/original_input_points_%d.png", region_id);
    //plot_points_into_file(filename, points_coord[PDLN_LON], points_coord[PDLN_LAT], num_kernel_points + num_expanded_points, PDLN_PLOT_GLOBAL);
    if(rotated_expanded_boundry != NULL) {
            //calculate_stereographic_projection(global_p_lon[0], global_p_lat[0], center[PDLN_LON], center[PDLN_LAT], global_p_lon[0], global_p_lat[0]);
            //calculate_stereographic_projection(global_p_lon[1], global_p_lat[1], center[PDLN_LON], center[PDLN_LAT], global_p_lon[1], global_p_lat[1]);
            //calculate_stereographic_projection(global_p_lon[2], global_p_lat[2], center[PDLN_LON], center[PDLN_LAT], global_p_lon[2], global_p_lat[2]);
            //calculate_stereographic_projection(global_p_lon[3], global_p_lat[3], center[PDLN_LON], center[PDLN_LAT], global_p_lon[3], global_p_lat[3]);
        triangulation = new Delaunay_Voronoi(num_kernel_points + num_expanded_points,
                                             projected_coord[PDLN_LON], projected_coord[PDLN_LAT], points_coord[PDLN_LON], points_coord[PDLN_LAT],
                                             points_global_index, false,
                                             rotated_expanded_boundry->min_lon, rotated_expanded_boundry->max_lon,
                                             rotated_expanded_boundry->min_lat, rotated_expanded_boundry->max_lat, NULL, virtual_point_local_index);

        if(node_type != PDLN_NODE_TYPE_COMMON) {
            //triangulation->plot_original_points_into_file(filename);

            /*
            char filename[64];
            int rank, mpi_size;
            MPI_Comm_rank(process_thread_mgr->get_mpi_comm(), &rank);
            MPI_Comm_size(process_thread_mgr->get_mpi_comm(), &mpi_size);
            snprintf(filename, 64, "log/projected_triangles_%d-%d.png", mpi_size, rank);
            triangulation->plot_projection_into_file(filename);
            */

            triangulation->update_all_points_coord(points_coord[PDLN_LON], points_coord[PDLN_LAT], num_kernel_points + num_expanded_points);
            //triangulation->correct_cyclic_triangles(cyclic_triangles, is_cyclic);
            triangulation->recognize_cyclic_triangles();
            triangulation->relegalize_all_triangles();

            if(PDLN_INSERT_VIRTUAL_POINT && polars_local_index->size() > 1)
                reset_polars();
        }
        else {
            Point circle_center;
            double radius;
            double x[3], y[3];

            calculate_real_boundary();
            //triangulation->plot_projection_into_file(filename);

            if(PDLN_REMOVE_UNNECESSARY_TRIANGLES && real_boundry->min_lat < 0) {
                calculate_stereographic_projection(0, real_boundry->min_lat, this->center[PDLN_LON], this->center[PDLN_LAT], x[0], y[0]);
                calculate_stereographic_projection(90, real_boundry->min_lat, this->center[PDLN_LON], this->center[PDLN_LAT], x[1], y[1]);
                calculate_stereographic_projection(180, real_boundry->min_lat, this->center[PDLN_LON], this->center[PDLN_LAT], x[2], y[2]);

                calculate_circle_center(x, y, &circle_center.x, &circle_center.y);
                radius = sqrt((x[2]-circle_center.x)*(x[2]-circle_center.x)+(y[2]-circle_center.y)*(y[2]-circle_center.y));

                if(radius < 100) {
                    //printf("[%d] + circle_center: (%lf, %lf), point: (%lf, %lf) radius: %lf\n", rank, circle_center.x, circle_center.y, x[2], y[2], radius);
                    triangulation->remove_triangles_in_circle(circle_center, radius);
                    double lon = (real_boundry->max_lon + real_boundry->min_lon + 360.0) * 0.5;
                    double head_lon, head_lat, tail_lon, tail_lat;
                    calculate_stereographic_projection(lon, -center[PDLN_LAT]-20 , center[PDLN_LON], center[PDLN_LAT], head_lon, head_lat);
                    calculate_stereographic_projection(lon, real_boundry->min_lat, center[PDLN_LON], center[PDLN_LAT], tail_lon, tail_lat);
                    //printf("[%d]real: (%lf, %lf, %lf, %lf), (%lf, %lf) -- (%lf, %lf)\n", rank, real_boundry->min_lon, real_boundry->max_lon,
                    //                                                               real_boundry->min_lat, real_boundry->max_lat,
                    //                                                             head_lon, head_lat, tail_lon, tail_lat);
                    triangulation->remove_triangles_on_segment(Point(head_lon, head_lat), Point(tail_lon, tail_lat));
                }
            }

            if(PDLN_REMOVE_UNNECESSARY_TRIANGLES && real_boundry->max_lat > 0) {
                calculate_stereographic_projection(0, real_boundry->max_lat, this->center[PDLN_LON], this->center[PDLN_LAT], x[0], y[0]);
                calculate_stereographic_projection(90, real_boundry->max_lat, this->center[PDLN_LON], this->center[PDLN_LAT], x[1], y[1]);
                calculate_stereographic_projection(180, real_boundry->max_lat, this->center[PDLN_LON], this->center[PDLN_LAT], x[2], y[2]);

                calculate_circle_center(x, y, &circle_center.x, &circle_center.y);
                radius = sqrt((x[2]-circle_center.x)*(x[2]-circle_center.x)+(y[2]-circle_center.y)*(y[2]-circle_center.y));

                if(radius < 100) {
                    triangulation->remove_triangles_in_circle(circle_center, radius);
                    double lon = (real_boundry->max_lon + real_boundry->min_lon + 360.0) * 0.5;
                    double head_lon, head_lat, tail_lon, tail_lat;
                    calculate_stereographic_projection(lon, real_boundry->max_lat+0.1, center[PDLN_LON], center[PDLN_LAT], head_lon, head_lat);
                    calculate_stereographic_projection(lon, -center[PDLN_LAT]+20 , center[PDLN_LON], center[PDLN_LAT], tail_lon, tail_lat);
                    triangulation->remove_triangles_on_segment(Point(head_lon, head_lat), Point(tail_lon, tail_lat));
                }
            }

            /*
            char filename[64];
            int rank, mpi_size;
            MPI_Comm_rank(process_thread_mgr->get_mpi_comm(), &rank);
            MPI_Comm_size(process_thread_mgr->get_mpi_comm(), &mpi_size);
            snprintf(filename, 64, "log/projected_triangles_%d-%d.png", mpi_size, rank);
            triangulation->plot_projection_into_file(filename);
            */

            triangulation->update_all_points_coord(points_coord[PDLN_LON], points_coord[PDLN_LAT], num_kernel_points + num_expanded_points);
            triangulation->remove_triangles_on_or_out_of_boundary(real_boundry->min_lon, real_boundry->max_lon, real_boundry->min_lat, real_boundry->max_lat);
            //triangulation->relegalize_all_triangles();
            triangulation->uncyclic_all_points();
            triangulation->recognize_cyclic_triangles();
        }
        if(num_inserted > 0)
            triangulation->remove_triangles_till(num_inserted);
    }
    else {
        triangulation = new Delaunay_Voronoi(num_kernel_points + num_expanded_points,
                                             points_coord[PDLN_LON], points_coord[PDLN_LAT], NULL, NULL, points_global_index, false,
                                             expanded_boundry->min_lon, expanded_boundry->max_lon,
                                             expanded_boundry->min_lat, expanded_boundry->max_lat, NULL, virtual_point_local_index);
        triangulation->uncyclic_all_points();
        triangulation->recognize_cyclic_triangles();
        if(num_inserted > 0)
            triangulation->remove_triangles_till(num_inserted);
    }
}


void Search_tree_node::reset_polars()
{
#ifdef DEBUG
    assert(node_type != PDLN_NODE_TYPE_COMMON);
#endif
    double reset_lat_value = node_type == PDLN_NODE_TYPE_NPOLAR ? 90 : -90;
    triangulation->update_points_coord_y(reset_lat_value, polars_local_index);
    triangulation->remove_triangles_only_containing_virtual_polar();
}


void Search_tree_node::count_points(double *coord[2], int *idx, int offset, int num, Midline midline,
                                    //double *c_points_coord[4], int *child_points_idx[2],
                                    int c_num_points[2])
{
    c_num_points[0] = 0;
    c_num_points[1] = 0;

    #pragma omp parallel for
    for(int i = offset; i < offset+num; i ++) {
        if(coord[midline.type][i] < midline.value) {
            //child_coord[midline.type][c_num_points[0]] = coord[midline.type][i];
            //child_coord[(midline.type+1)%2][c_num_points[0]] = coord[(midline.type+1)%2][i];
            //child_points_idx[0][c_num_points[0]] = idx[i];
            #pragma omp critical
            c_num_points[0]++;
        }
        else {
            //child_coord[2+midline.type][c_num_points[1]] = coord[midline.type][i];
            //child_coord[2+(midline.type+1)%2][c_num_points[1]] = coord[(midline.type+1)%2][i];
            //child_points_idx[1][c_num_points[1]] = idx[i];
            #pragma omp critical
            c_num_points[1]++;
        }
    }
}


void Search_tree_node::split_local_points(Midline midline, double *c_points_coord[4], int *child_points_idx[2], int c_num_points[2])
{
    c_num_points[0] = 0;
    c_num_points[1] = 0;

    if(non_monotonic && midline.type == PDLN_LON)
        assert(false);
    else {
        #pragma omp parallel for
        for(int i = 0; i < num_kernel_points; i ++) {
            if(points_coord[midline.type][i] < midline.value) {
                c_points_coord[midline.type][c_num_points[0]] = points_coord[midline.type][i];
                c_points_coord[(midline.type+1)%2][c_num_points[0]] = points_coord[(midline.type+1)%2][i];
                #pragma omp critical
                child_points_idx[0][c_num_points[0]++] = points_global_index[i];
            }
            else {
                c_points_coord[2+midline.type][c_num_points[1]] = points_coord[midline.type][i];
                c_points_coord[2+(midline.type+1)%2][c_num_points[1]] = points_coord[(midline.type+1)%2][i];
                #pragma omp critical 
                child_points_idx[1][c_num_points[1]++] = points_global_index[i];
            }
        }
    }
}


void Search_tree_node::decompose_by_processing_units_number(double *workloads, double *c_points_coord[4], int *child_points_idx[2], int c_num_points[2],
                                   Boundry c_boundry[2], vector<int> child_proc_units_id[2], int mode, int *c_intervals[2], int c_num_intervals[2])
{
    double length[2], boundry_values[4], child_total_workload[2];
    Midline midline;
    unsigned i;
    int iteration_count;

    assert(region_ids.size() > 1);
    boundry_values[PDLN_LON] = kernel_boundry->min_lon;
    boundry_values[PDLN_LAT] = kernel_boundry->min_lat;
    boundry_values[PDLN_LON+2] = kernel_boundry->max_lon;
    boundry_values[PDLN_LAT+2] = kernel_boundry->max_lat;

    if(non_monotonic) {
        assert(false);
        boundry_values[PDLN_LON] -= 360.0;
    }

    assert(boundry_values[PDLN_LON] != boundry_values[PDLN_LON+2]);
    assert(boundry_values[PDLN_LAT] < boundry_values[PDLN_LAT+2]);
    length[0] = boundry_values[PDLN_LON+2] - boundry_values[PDLN_LON];
    length[1] = boundry_values[PDLN_LAT+2] - boundry_values[PDLN_LAT];
    assert(length[0] >= 0 || length[1] >= 0);
    assert(length[0] <= (360.0+PDLN_HIGH_BOUNDRY_SHIFTING) && length[0] >= 0.0 && length[1] <= (180.0+PDLN_HIGH_BOUNDRY_SHIFTING) && length[1] >= 0.0);
    
    if(mode == PDLN_DECOMPOSE_SPOLAR_MODE || mode == PDLN_DECOMPOSE_NPOLAR_MODE)
        midline.type = PDLN_LAT;
    else if(length[1] > length[0])
        midline.type = PDLN_LAT;
    else
        midline.type = PDLN_LON;

    child_proc_units_id[0].clear();
    child_proc_units_id[1].clear();

    /* PDLN_DECOMPOSE_COMMON_MODE: 0   1   2   3 | 4   5   6   7
     * PDLN_DECOMPOSE_SPOLAR_MODE: 0 | 1   2   3   4   5   6   7
     * PDLN_DECOMPOSE_NPOLAR_MODE: 0   1   2   3   4   5   6 | 7 */
    if(mode == PDLN_DECOMPOSE_COMMON_MODE) {
        assert(c_intervals);
        assert(c_num_intervals);
        int mid_idx = 0; 
        if(num_groups == 1) {
            mid_idx = region_ids.size()/2;
            c_intervals[0] = c_intervals[1] = NULL;
            c_num_intervals[0] = c_num_intervals[1] = 1;
        }
        else {
            for(i = 0; i < num_groups/2; i++)
                mid_idx += group_intervals[i];
            c_intervals[0] = group_intervals;
            c_intervals[1] = group_intervals+num_groups/2;
            c_num_intervals[0] = num_groups/2;
            c_num_intervals[1] = num_groups - c_num_intervals[0];
        }

        for(i = 0; i < mid_idx; i++)
            child_proc_units_id[0].push_back(region_ids[i]);
        for(; i < region_ids.size(); i++)
            child_proc_units_id[1].push_back(region_ids[i]);
    }
    else if(mode == PDLN_DECOMPOSE_SPOLAR_MODE) {
        child_proc_units_id[0].push_back(region_ids[0]);
        for(i = 1; i < region_ids.size(); i++)
            child_proc_units_id[1].push_back(region_ids[i]);
    }
    else if(mode == PDLN_DECOMPOSE_NPOLAR_MODE) {
        for(i = 0; i < region_ids.size()-1; i++)
            child_proc_units_id[0].push_back(region_ids[i]);
        child_proc_units_id[1].push_back(region_ids[i]);
    }
    else
        assert(false);

#ifdef DEBUG
    assert(child_proc_units_id[0].size() + child_proc_units_id[1].size() == region_ids.size());
#endif

    if(region_ids.size() > 1) {
        for(i = 0, child_total_workload[0] = 0.0; i < child_proc_units_id[0].size(); i++) {
            child_total_workload[0] += workloads[child_proc_units_id[0][i]];
        }
        for(i = 0, child_total_workload[1] = 0.0; i < child_proc_units_id[1].size(); i++)
            child_total_workload[1] += workloads[child_proc_units_id[1][i]];

        midline.value = boundry_values[midline.type] + length[midline.type] * child_total_workload[0] / (child_total_workload[0] + child_total_workload[1]);
        split_local_points(midline, c_points_coord, child_points_idx, c_num_points);

        iteration_count = 1;
        while (c_num_points[0] == 0 || c_num_points[1] == 0 ||
               fabs(c_num_points[0]/c_num_points[1] - child_total_workload[0]/child_total_workload[1]) > PDLN_TOLERABLE_ERROR) {

            if(iteration_count++ > PDLN_MAX_ITER_COUNT)
                break;

            if(c_num_points[0] < child_total_workload[0]) {
#ifdef DEBUG
                assert(c_num_points[1] >= child_total_workload[1] || fabs(child_total_workload[1] - c_num_points[1]) < PDLN_FLOAT_EQ_ERROR);
#endif
                midline.value += (boundry_values[2+midline.type] - midline.value) * (c_num_points[1] - child_total_workload[1]) / c_num_points[1];
            }
            else {
#ifdef DEBUG
                assert(c_num_points[1] <= child_total_workload[1] || fabs(child_total_workload[1] - c_num_points[1]) < PDLN_FLOAT_EQ_ERROR);
#endif
                midline.value -= (midline.value - boundry_values[midline.type]) * (c_num_points[0] - child_total_workload[0]) / c_num_points[0];
            }
            assert(midline.value > boundry_values[midline.type]);
            assert(midline.value < boundry_values[2+midline.type]);
            /* TODO: Search only half of the whole points, but not the whole points */
            split_local_points(midline, c_points_coord, child_points_idx, c_num_points);
        }
    }
    else
        midline.value = boundry_values[2+midline.type];

    if(midline.type == PDLN_LON) {
        c_boundry[0].min_lat = c_boundry[1].min_lat = kernel_boundry->min_lat;
        c_boundry[0].max_lat = c_boundry[1].max_lat = kernel_boundry->max_lat;
        c_boundry[0].min_lon = kernel_boundry->min_lon;
        c_boundry[0].max_lon = c_boundry[1].min_lon = midline.value;
        c_boundry[1].max_lon = kernel_boundry->max_lon;
    }
    else if(midline.type == PDLN_LAT) {
        c_boundry[0].min_lon = c_boundry[1].min_lon = kernel_boundry->min_lon;
        c_boundry[0].max_lon = c_boundry[1].max_lon = kernel_boundry->max_lon;
        c_boundry[0].min_lat = kernel_boundry->min_lat;
        c_boundry[0].max_lat = c_boundry[1].min_lat = midline.value;
        c_boundry[1].max_lat = kernel_boundry->max_lat;
    }
    else
        assert(false);

}


void Search_tree_node::decompose_by_fixed_longitude(double fixed_lon, double *workloads, double *c_points_coord[4], int *child_points_idx[2], int c_num_points[2],
                                                    Boundry c_boundry[2], vector<int> child_proc_units_id[2])
{
    Midline midline;

    assert(region_ids.size() > 0);
    assert(fixed_lon >= 0 && fixed_lon < 360);
    
    midline.type = PDLN_LON;
    midline.value = fixed_lon;

    split_local_points(midline, c_points_coord, child_points_idx, c_num_points);

    split_processing_units_by_points_number(workloads, c_num_points[0], c_num_points[1], region_ids, child_proc_units_id);

    assert(child_proc_units_id[0].size() + child_proc_units_id[1].size() == region_ids.size());

    c_boundry[0] = c_boundry[1] = *kernel_boundry;
    c_boundry[0].max_lon = c_boundry[1].min_lon = fixed_lon;
}


void Search_tree_node::split_processing_units_by_points_number(double *workloads, int left_num_points, int right_num_points,
                                                               vector<int> parent_proc_units_id, vector<int> child_proc_units_id[2])
{
    unsigned int splitting_index;
    double left_workloads, prev_left_workloads;

#ifdef DEBUG
    double total_workloads = 0;

    for(unsigned int i = 0; i < parent_proc_units_id.size(); i++)
        total_workloads += workloads[parent_proc_units_id[i]];
    assert(std::abs(total_workloads - left_num_points - right_num_points) < PDLN_FLOAT_EQ_ERROR);
#endif

    child_proc_units_id[0].clear();
    child_proc_units_id[1].clear();

    for(splitting_index = 0, left_workloads = 0; splitting_index < parent_proc_units_id.size(); splitting_index++) {
        left_workloads += workloads[parent_proc_units_id[splitting_index]];
        if(left_workloads > left_num_points)
            break;
    }
    prev_left_workloads = left_workloads - workloads[parent_proc_units_id[splitting_index]];

    if(std::abs(prev_left_workloads - left_num_points) < std::abs(left_workloads - left_num_points))
        splitting_index--;

    for(unsigned int i = 0; i <= splitting_index; i++)
        child_proc_units_id[0].push_back(parent_proc_units_id[i]);

    for(unsigned int i = splitting_index+1; i < parent_proc_units_id.size(); i++)
        child_proc_units_id[1].push_back(parent_proc_units_id[i]);
}


void Search_tree_node::add_expanded_points(double *lon_value, double *lat_value, int *global_idx, int num_points)
{
    double *coord_value[2];
    coord_value[PDLN_LON] = lon_value;
    coord_value[PDLN_LAT] = lat_value;
    add_expanded_points(coord_value, global_idx, num_points);
}


void Search_tree_node::add_expanded_points(double *coord_value[2], int *global_idx, int num_points)
{
    double *tmp_coord_value[2];
    int *tmp_idx;
    if(num_expanded_points + num_points > len_expanded_points_coord_buf) {
        len_expanded_points_coord_buf = num_expanded_points + num_points * 4 * PDLN_EXPECTED_EXPANDING_LOOP_TIMES;
        tmp_coord_value[0] = new double[num_kernel_points + len_expanded_points_coord_buf];
        tmp_coord_value[1] = new double[num_kernel_points + len_expanded_points_coord_buf];
        tmp_idx = new int[num_kernel_points + len_expanded_points_coord_buf];
        memcpy(tmp_coord_value[0], points_coord[0], sizeof(double) * (num_kernel_points + num_expanded_points));
        memcpy(tmp_coord_value[1], points_coord[1], sizeof(double) * (num_kernel_points + num_expanded_points));
        memcpy(tmp_idx, points_global_index, sizeof(int) * (num_kernel_points + num_expanded_points));
        delete[] points_coord[0];
        delete[] points_coord[1];
        delete[] points_global_index;
        points_coord[0] = tmp_coord_value[0];
        points_coord[1] = tmp_coord_value[1];
        points_global_index = tmp_idx;
        if(projected_coord[0] != NULL) {
            tmp_coord_value[0] = new double[num_kernel_points + len_expanded_points_coord_buf];
            tmp_coord_value[1] = new double[num_kernel_points + len_expanded_points_coord_buf];
            memcpy(tmp_coord_value[0], projected_coord[0], sizeof(double) * (num_kernel_points + num_expanded_points));
            memcpy(tmp_coord_value[1], projected_coord[1], sizeof(double) * (num_kernel_points + num_expanded_points));
            delete[] projected_coord[0];
            delete[] projected_coord[1];
            projected_coord[0] = tmp_coord_value[0];
            projected_coord[1] = tmp_coord_value[1];
        }
    }
    memcpy(points_coord[0] + num_kernel_points + num_expanded_points, coord_value[0], sizeof(double) * num_points);
    memcpy(points_coord[1] + num_kernel_points + num_expanded_points, coord_value[1], sizeof(double) * num_points);
    memcpy(points_global_index + num_kernel_points + num_expanded_points, global_idx, sizeof(int) * num_points);

    //MPI_Comm_rank(process_thread_mgr->get_mpi_comm(), &rank);
    fix_expanded_boundry(num_kernel_points + num_expanded_points, num_points);
    //int rank;
    //printf("[%d] fix: %lf, %lf, %lf, %lf\n", rank, expanded_boundry->min_lon, expanded_boundry->max_lon, expanded_boundry->min_lat, expanded_boundry->max_lat);

    num_expanded_points += num_points;
    assert(num_expanded_points >= num_points);
}


void Search_tree_node::calculate_real_boundary()
{
    Boundry boundry;
    boundry.min_lat = 1e10;
    boundry.max_lat = -1e10;
    boundry.min_lon = 1e10;
    boundry.max_lon = -1e10;
    for(int i = 0; i < num_kernel_points + num_expanded_points; i++) {
        if(points_coord[PDLN_LON][i] < boundry.min_lon) boundry.min_lon = points_coord[PDLN_LON][i];
        if(points_coord[PDLN_LON][i] > boundry.max_lon) boundry.max_lon = points_coord[PDLN_LON][i];
    }
    for(int i = 0; i < num_kernel_points + num_expanded_points; i++) {
        if(points_coord[PDLN_LAT][i] < boundry.min_lat) boundry.min_lat = points_coord[PDLN_LAT][i];
        if(points_coord[PDLN_LAT][i] > boundry.max_lat) boundry.max_lat = points_coord[PDLN_LAT][i];
    }

    //printf("%lf, %lf, %lf, %lf\n", boundry.min_lon, boundry.max_lon, boundry.min_lat, boundry.max_lat);

    if(real_boundry == NULL)
        real_boundry = new Boundry();
    *real_boundry = boundry;
}


void Boundry::move_close(double *coord[2], int offset, int count)
{
    double min_lon = 1e10;
    double max_lon = -1e10;
    double min_lat = 1e10;
    double max_lat = -1e10;

    for(int i = offset; i < offset + count; i++) {
        if (min_lon > coord[PDLN_LON][i]) min_lon = coord[PDLN_LON][i];
        if (max_lon < coord[PDLN_LON][i]) max_lon = coord[PDLN_LON][i];
        if (min_lat > coord[PDLN_LAT][i]) min_lat = coord[PDLN_LAT][i];
        if (max_lat < coord[PDLN_LAT][i]) max_lat = coord[PDLN_LAT][i];
    }
}


void Search_tree_node::fix_expanded_boundry(int index, int count)
{
    double min_lon = 1e10;
    double max_lon = -1e10;
    double min_lat = 1e10;
    double max_lat = -1e10;

    for(int i = index; i < index + count; i++) {
        if (min_lon > points_coord[PDLN_LON][i]) min_lon = points_coord[PDLN_LON][i];
        if (max_lon < points_coord[PDLN_LON][i]) max_lon = points_coord[PDLN_LON][i];
        if (min_lat > points_coord[PDLN_LAT][i]) min_lat = points_coord[PDLN_LAT][i];
        if (max_lat < points_coord[PDLN_LAT][i]) max_lat = points_coord[PDLN_LAT][i];
    }
    if (node_type == PDLN_DECOMPOSE_COMMON_MODE)
        expanded_boundry->max(min_lon, max_lon + PDLN_HIGH_BOUNDRY_SHIFTING, min_lat, max_lat + PDLN_HIGH_BOUNDRY_SHIFTING);
    else if (node_type == PDLN_DECOMPOSE_SPOLAR_MODE)
        expanded_boundry->max_lat = std::max(expanded_boundry->max_lat, max_lat + PDLN_HIGH_BOUNDRY_SHIFTING);
    else if (node_type == PDLN_DECOMPOSE_NPOLAR_MODE)
        expanded_boundry->min_lat = std::min(expanded_boundry->min_lat, min_lat);
}


bool operator == (pair<Search_tree_node*, bool> p1, Search_tree_node* p2)
{
    return p1.first == p2;
}


void Search_tree_node::add_neighbors(vector<Search_tree_node*> ns)
{
    for(unsigned int i = 0; i < ns.size(); i++) {
        if(region_id == ns[i]->region_id)
            continue;

        if(find(neighbors.begin(), neighbors.end(), ns[i]) == neighbors.end())
            neighbors.push_back(pair<Search_tree_node*, bool>(ns[i], false));
    }
}


void Search_tree_node::project_grid()
{
    if(projected_coord[0] == NULL) { /* first time to generate rotated grid */
        projected_coord[0] = new double[num_kernel_points + len_expanded_points_coord_buf];
        projected_coord[1] = new double[num_kernel_points + len_expanded_points_coord_buf];

        #pragma omp parallel for
        for(int i = 0; i < num_kernel_points + num_expanded_points; i++) {
            calculate_stereographic_projection(points_coord[PDLN_LON][i], points_coord[PDLN_LAT][i], center[PDLN_LON], center[PDLN_LAT],
                                               projected_coord[PDLN_LON][i], projected_coord[PDLN_LAT][i]);
        }

        num_projected_points = num_kernel_points + num_expanded_points;
    }
    else {
        #pragma omp parallel for
        for(int i = num_projected_points; i < num_kernel_points + num_expanded_points; i++) {
            calculate_stereographic_projection(points_coord[PDLN_LON][i], points_coord[PDLN_LAT][i], center[PDLN_LON], center[PDLN_LAT],
                                               projected_coord[PDLN_LON][i], projected_coord[PDLN_LAT][i]);
        }
        num_projected_points = num_kernel_points + num_expanded_points;
    }

    /* recalculate expanded boundary */
    double top = -1e20, bot = 1e20, left = 1e20, right = -1e20;
    for(int i = 0; i < num_projected_points; i++) { //TODO: i can be started from non-zero
        if (projected_coord[PDLN_LON][i] < left)  left = projected_coord[PDLN_LON][i];
        if (projected_coord[PDLN_LON][i] > right) right = projected_coord[PDLN_LON][i];
        if (projected_coord[PDLN_LAT][i] < bot) bot = projected_coord[PDLN_LAT][i];
        if (projected_coord[PDLN_LAT][i] > top) top = projected_coord[PDLN_LAT][i];
    }
    if(rotated_expanded_boundry != NULL) {
        assert(rotated_expanded_boundry->min_lon >= left && rotated_expanded_boundry->max_lon <= right && rotated_expanded_boundry->min_lat >= bot && rotated_expanded_boundry->max_lat <= top);
        delete rotated_expanded_boundry;
    }
    rotated_expanded_boundry = new Boundry(left, right, bot, top);
}


Delaunay_grid_decomposition::Delaunay_grid_decomposition(int grid_id, Processing_resource *proc_info, int min_points_per_chunk)
    : original_grid(grid_id)
    , processing_info(proc_info)
    , min_points_per_chunk(min_points_per_chunk)
    , workloads(NULL)
    , average_workload(0)
    , regionID_to_unitID(NULL)
    , all_group_intervals(NULL)
{
    double **coords;
    double *coord_values[2];
    int *global_index;
    Boundry boundry;

    assert(processing_info != NULL);

    is_cyclic  = grid_info_mgr->is_grid_cyclic(grid_id);
    coords     = grid_info_mgr->get_grid_coord_values(grid_id);
    num_points = grid_info_mgr->get_grid_num_points(grid_id);
    grid_info_mgr->get_grid_boundry(grid_id, &boundry.min_lon, &boundry.max_lon, &boundry.min_lat, &boundry.max_lat);
    coord_values[0] = coords[0];
    coord_values[1] = coords[1];

    active_processing_units_flag = new bool[processing_info->get_num_total_processing_units()];
    /* coord_values[0] coord_values[1] will be changed */
    num_inserted = insert_virtual_points(coord_values, &boundry, num_points);
    num_points += num_inserted;

    global_index = new int[num_points];
    for(int i = 0; i < num_points; i++)
        global_index[i] = i;

    if(boundry.max_lon - boundry.min_lon < 360.0)
        boundry.max_lon += PDLN_HIGH_BOUNDRY_SHIFTING;
    boundry.max_lat += PDLN_HIGH_BOUNDRY_SHIFTING;

    assert(boundry.max_lon - boundry.min_lon <= 360.0);
    search_tree_root = new Search_tree_node(NULL, coord_values, global_index, num_points, boundry, PDLN_NODE_TYPE_COMMON);
    search_tree_root->calculate_real_boundary();

    //if(processing_info->get_local_process_id() == 0) {
    //    char filename[64];
    //    snprintf(filename, 64, "log/original_input_points.png");
    //    plot_points_into_file(filename, coord_values[PDLN_LON], coord_values[PDLN_LAT], num_points, PDLN_PLOT_GLOBAL);
    //}
    
    if(num_inserted > 0) {
        delete coord_values[0];
        delete coord_values[1];
    }
    delete global_index;

    c_points_coord[0] = c_points_coord[1] = c_points_coord[2] = c_points_coord[3] = NULL;
    c_points_index[0] = c_points_index[1] = NULL;
}


Delaunay_grid_decomposition::~Delaunay_grid_decomposition()
{
    delete search_tree_root;
    if(workloads != NULL)
        delete[] workloads; 
    delete[] active_processing_units_flag;
    delete[] all_group_intervals;

    for(int i = 0; i < 4; i++)
        delete[] c_points_coord[i];
    delete[] c_points_index[0];
    delete[] c_points_index[1];
}


#define PDLN_GVPOINT_DENSITY  (20)
#define PDLN_INSERT_EXPAND_RATIO (0.01)
int Delaunay_grid_decomposition::insert_virtual_points(double *coord_values[2], Boundry *boundry, int num_points)
{
    if(is_cyclic && fabs(boundry->max_lat - 90) < PDLN_FLOAT_EQ_ERROR && fabs(boundry->min_lat - -90) < PDLN_FLOAT_EQ_ERROR)
        return 0;

    double minX = coord_values[PDLN_LON][0];
    double maxX = coord_values[PDLN_LON][0];
    double minY = coord_values[PDLN_LAT][0];
    double maxY = coord_values[PDLN_LAT][0];
    for(int i = 0; i < num_points; i++)
    {
        if (coord_values[PDLN_LON][i] < minX) minX = coord_values[PDLN_LON][i];
        if (coord_values[PDLN_LON][i] > maxX) maxX = coord_values[PDLN_LON][i];
        if (coord_values[PDLN_LAT][i] < minY) minY = coord_values[PDLN_LAT][i];
        if (coord_values[PDLN_LAT][i] > maxY) maxY = coord_values[PDLN_LAT][i];
    }

    double dx   = maxX - minX;
    double dy   = maxY - minY;
    double dMax = std::max(dx, dy);

    double v_minx = minX-dMax*PDLN_INSERT_EXPAND_RATIO;
    double v_maxx = maxX+dMax*PDLN_INSERT_EXPAND_RATIO;
    double v_miny = minY-dMax*PDLN_INSERT_EXPAND_RATIO;
    double v_maxy = maxY+dMax*PDLN_INSERT_EXPAND_RATIO;

    double r_minx = is_cyclic && v_minx < 0   ?   1 : v_minx;
    double r_maxx = is_cyclic && v_maxx > 360 ? 359 : v_maxx;
    double r_miny = v_miny < -90 ? -89 : v_miny;
    double r_maxy = v_maxy > 90  ?  89 : v_maxy;
    
    /*
     * x * y = num_points
     * x : y = dx : dy
     */
    unsigned num_x = (unsigned)sqrt(num_points * dx / (double)dy);
    unsigned num_y = num_x * dy / dx;
    num_x /= PDLN_GVPOINT_DENSITY;
    num_y /= PDLN_GVPOINT_DENSITY;

    int num_inserted = 0;
    double *inserted_coord[2];
    inserted_coord[0] = new double[num_points + 2*num_x + 2*num_y];
    inserted_coord[1] = new double[num_points + 2*num_x + 2*num_y];
    if(fabs(boundry->max_lat - 90) >= PDLN_FLOAT_EQ_ERROR && v_maxy < 90) {
        for(unsigned i = 1; i < num_x-1; i++) {
            inserted_coord[PDLN_LON][num_inserted] = r_minx+(r_maxx-r_minx)/num_x*i;
            inserted_coord[PDLN_LAT][num_inserted++] = r_maxy;
        }
        if(boundry->max_lat < r_maxy) boundry->max_lat = r_maxy;
    }
    if(fabs(boundry->min_lat - -90) >= PDLN_FLOAT_EQ_ERROR && v_miny > -90) {
        for(unsigned i = 1; i < num_x-1; i++) {
            inserted_coord[PDLN_LON][num_inserted] = r_minx+(r_maxx-r_minx)/num_x*i;
            inserted_coord[PDLN_LAT][num_inserted++] = r_miny;
        }
        if(boundry->min_lat > r_miny) boundry->min_lat = r_miny;
    }
    if(!is_cyclic) {
        for(unsigned i = 1; i < num_y-1; i++) {
            inserted_coord[PDLN_LON][num_inserted] = v_minx;
            inserted_coord[PDLN_LAT][num_inserted++] = v_miny+(v_maxy-v_miny)/num_y*i;
        }
        for(unsigned i = 1; i < num_y-1; i++) {
            inserted_coord[PDLN_LON][num_inserted] = v_maxx;
            inserted_coord[PDLN_LAT][num_inserted++] = v_miny+(v_maxy-v_miny)/num_y*i;
        }
        if(boundry->min_lon > r_minx) boundry->min_lon = r_minx;
        if(boundry->max_lon < r_maxx) boundry->max_lon = r_maxx;
    }
    assert(num_inserted <= 2*num_x + 2*num_y);

    if(num_inserted > 0) {
        /* store inserted points first. Then original points follow */
        memcpy(inserted_coord[PDLN_LON]+num_inserted, coord_values[PDLN_LON], num_points*sizeof(double));
        memcpy(inserted_coord[PDLN_LAT]+num_inserted, coord_values[PDLN_LAT], num_points*sizeof(double));

        coord_values[PDLN_LON] = inserted_coord[PDLN_LON];
        coord_values[PDLN_LAT] = inserted_coord[PDLN_LAT];
    }
    else {
        delete inserted_coord[0];
        delete inserted_coord[1];
    }
    return num_inserted;
}


void Delaunay_grid_decomposition::initialze_workload()
{
    int max_punits   = (num_points + min_points_per_chunk - 1) / min_points_per_chunk;
    int total_punits = processing_info->get_num_total_processing_units();
    num_regions      = std::max(std::min(total_punits, max_punits), 4);
    average_workload = (double)num_points / num_regions;
    std::vector<int> all_regions_id;

    assert(min_points_per_chunk > 0);
    delete[] regionID_to_unitID;
    delete[] workloads;

    /* +2: space for extra regions of polars' processes */
    regionID_to_unitID = new int[num_regions+2];
    workloads          = new double[num_regions+2];

    /* the first and the last ID(0 and num_regions+1) are preserved, for the same purpose */
    for(int i = 1; i < num_regions+1; i++) {
        all_regions_id.push_back(i);
        workloads[i] = average_workload;
    }

    processing_info->pick_out_active_processing_units(num_regions, active_processing_units_flag);

    /* offset for the same purpose above */
    const int offset = 1;
    if(total_punits > 3) {
        for(int i = 0, regionID = 1; i < total_punits; i++)
            if(active_processing_units_flag[i])
                regionID_to_unitID[regionID++] = i;
    }
    else if(total_punits == 1)
        regionID_to_unitID[0+offset] = regionID_to_unitID[1+offset] = regionID_to_unitID[2+offset] = regionID_to_unitID[3+offset] = 0;
    else if(total_punits == 2) {
        regionID_to_unitID[0+offset] = regionID_to_unitID[1+offset] = 0;
        regionID_to_unitID[2+offset] = regionID_to_unitID[3+offset] = 1;
    }
    else if(total_punits == 3) {
        regionID_to_unitID[0+offset] = 0;
        regionID_to_unitID[1+offset] = regionID_to_unitID[2+offset] = 1;
        regionID_to_unitID[3+offset] = 2;
    }

    search_tree_root->update_region_ids(all_regions_id);
    assert(search_tree_root->region_ids.size() > 0);
}


void Delaunay_grid_decomposition::update_workloads(int total_workload, vector<int> &ids)
{
    /* NOTE: This will make the final workload of leaf node not exactly more than min_points_per_chunk */
    if(ids.size() == 1) {
        workloads[ids[0]] = total_workload;
        return;
    }

    double old_total_workload = 0.0;
    for(unsigned int i = 0; i < ids.size(); i++)
        old_total_workload += workloads[ids[i]];

    for(unsigned int i = 0; i < ids.size(); i++)
        workloads[ids[i]] = workloads[ids[i]] * total_workload / old_total_workload;

    for(unsigned int i = 0; i < ids.size();) {
        if(ids.size() < 2)
            break;

        if(workloads[ids[i]] < min_points_per_chunk) {
            double unassigned_workload = workloads[ids[i]];
            workloads[ids[i]] = 0;
            active_processing_units_flag[ids[i]] = false;
            ids.erase(ids.begin() + i);
            for(unsigned int j = 0; j < ids.size(); j++)
                workloads[ids[j]] += unassigned_workload / ids.size();
        }
        else
            i++;
    }

#ifdef DEBUG
    assert(ids.size() > 0);
#endif
}


void Delaunay_grid_decomposition::initialze_buffer()
{
    for(int i = 0; i < 4; i++)
        delete[] c_points_coord[i];
    delete[] c_points_index[0];
    delete[] c_points_index[1];

    c_regions_id[0].reserve(num_regions+2);
    c_regions_id[1].reserve(num_regions+2);
    for(int i = 0; i < 4; i++)
        c_points_coord[i] = new double[search_tree_root->num_kernel_points];
    c_points_index[0] = new int[search_tree_root->num_kernel_points];
    c_points_index[1] = new int[search_tree_root->num_kernel_points];
}


/* "common_node" means non-polar node */
void Delaunay_grid_decomposition::decompose_common_node_recursively(Search_tree_node *node, bool lazy_mode)
{
    int c_num_points[2];
    Boundry c_boundry[2];
    int *c_intervals[2];
    int c_num_intervals[2];

    assert(node->region_ids.size() > 0);
    if(node->region_ids.size() == 1) {
        if(have_local_region_ids(node->region_ids))
            local_leaf_nodes.push_back(node);
        all_leaf_nodes.push_back(node);
        return;
    }

    node->decompose_by_processing_units_number(workloads, c_points_coord, c_points_index, 
                                               c_num_points, c_boundry, c_regions_id, 
                                               PDLN_DECOMPOSE_COMMON_MODE, c_intervals, 
                                               c_num_intervals);

    node->children[0] = alloc_search_tree_node(node, c_points_coord,   c_points_index[0], c_num_points[0], c_boundry[0], c_regions_id[0], PDLN_NODE_TYPE_COMMON);
    node->children[2] = alloc_search_tree_node(node, c_points_coord+2, c_points_index[1], c_num_points[1], c_boundry[1], c_regions_id[1], PDLN_NODE_TYPE_COMMON);

    node->children[0]->set_groups(c_intervals[0], c_num_intervals[0]);
    node->children[2]->set_groups(c_intervals[1], c_num_intervals[1]);
    //int rank;
    //MPI_Comm_rank(processing_info->get_mpi_comm(), &rank);
    //printf("[Rank%d]x[ST-INFO-PRE] p: %p, first: %p, third: %p\n", rank, node, node->children[0], node->children[2]);

    if(!lazy_mode || have_local_region_ids(node->children[0]->region_ids))
        decompose_common_node_recursively(node->children[0], lazy_mode);
    if(!lazy_mode || have_local_region_ids(node->children[2]->region_ids))
        decompose_common_node_recursively(node->children[2], lazy_mode);
}


Search_tree_node* Delaunay_grid_decomposition::alloc_search_tree_node(Search_tree_node* parent, double *coord_values[2], int *index, 
                                                                      int num_points, Boundry boundary, vector<int> &common_ids, int type)
{
    assert(common_ids.size() > 0);
    Search_tree_node *new_node = new Search_tree_node(parent, coord_values, index, num_points, boundary, type);

    /* common_ids can be modified by update_workloads */
    update_workloads(num_points, common_ids);

    new_node->update_region_ids(common_ids);
    if(common_ids.size() == 1)
        new_node->region_id = common_ids[0];
    return new_node;
}


bool Delaunay_grid_decomposition::do_two_regions_overlap(Boundry region1, Boundry region2)
{
    if(region1.max_lat <= region2.min_lat || region1.min_lat >= region2.max_lat)
        return false;
    if(region1.max_lon <= region2.min_lon || region1.min_lon >= region2.max_lon)
        return false;

    return true;
}

#define PDLN_BOUNDRY_TYPE_CLEAR     (0x0FFFFFFF)
#define PDLN_BOUNDRY_TYPE_NON       (0x00000000)
#define PDLN_BOUNDRY_TYPE_U         (0x10000000)
#define PDLN_BOUNDRY_TYPE_D         (0x20000000)
#define PDLN_BOUNDRY_TYPE_L         (0x40000000)
#define PDLN_BOUNDRY_TYPE_R         (0x80000000)
#define PDLN_BOUNDRY_TYPE_LR        (0xC0000000)
#define PDLN_BOUNDRY_TYPE_INVALID   (0xF0000000)
#define set_boundry_type(val, type) ((val & PDLN_BOUNDRY_TYPE_CLEAR) | type)
#define get_boundry_type(val)       (val & ~PDLN_BOUNDRY_TYPE_CLEAR)
/*
 *      ┌───┐
 *      │ 0 │
 *  ┌───┼───┼───┐
 *  │ 1 │   │ 3 │
 *  └───┼───┼───┘
 *      │ 2 │
 *      └───┘
 */
unsigned Delaunay_grid_decomposition::compute_common_boundry(Search_tree_node *tree_node, Search_tree_node *neighbor_node, Point *boundry_head,
                                                             Point *boundry_tail, Point *cyclic_boundry_head, Point *cyclic_boundry_tail)
{
    unsigned boundry_type = 0;
    double coord_value[2][2];

    coord_value[0][PDLN_LAT] = coord_value[0][PDLN_LON] = coord_value[1][PDLN_LAT] = coord_value[1][PDLN_LON] = PDLN_DOUBLE_INVALID_VALUE;
    if(tree_node->kernel_boundry->max_lat == neighbor_node->kernel_boundry->min_lat) { // Case 0
        if(std::max(tree_node->kernel_boundry->min_lon, neighbor_node->kernel_boundry->min_lon) <
           std::min(tree_node->kernel_boundry->max_lon, neighbor_node->kernel_boundry->max_lon)) {
            tree_node->num_neighbors_on_boundry[PDLN_UP]++;
            boundry_type |= PDLN_BOUNDRY_TYPE_U;
            coord_value[0][PDLN_LAT] = coord_value[1][PDLN_LAT] = tree_node->kernel_boundry->max_lat;
            coord_value[0][PDLN_LON] = std::max(tree_node->kernel_boundry->min_lon, neighbor_node->kernel_boundry->min_lon);
            coord_value[1][PDLN_LON] = std::min(tree_node->kernel_boundry->max_lon, neighbor_node->kernel_boundry->max_lon);
        }
    }
    else if(tree_node->kernel_boundry->min_lon == neighbor_node->kernel_boundry->max_lon) { // Case 1
        if(std::max(tree_node->kernel_boundry->min_lat, neighbor_node->kernel_boundry->min_lat) <
           std::min(tree_node->kernel_boundry->max_lat, neighbor_node->kernel_boundry->max_lat)) {
            tree_node->num_neighbors_on_boundry[PDLN_LEFT]++;
            boundry_type |= PDLN_BOUNDRY_TYPE_L;
            coord_value[0][PDLN_LON] = coord_value[1][PDLN_LON] = neighbor_node->kernel_boundry->max_lon;
            coord_value[0][PDLN_LAT] = std::max(tree_node->kernel_boundry->min_lat, neighbor_node->kernel_boundry->min_lat);
            coord_value[1][PDLN_LAT] = std::min(tree_node->kernel_boundry->max_lat, neighbor_node->kernel_boundry->max_lat);
        }
    }
    else if(tree_node->kernel_boundry->min_lat == neighbor_node->kernel_boundry->max_lat)  { // Case 2
        if(std::max(tree_node->kernel_boundry->min_lon, neighbor_node->kernel_boundry->min_lon) < 
           std::min(tree_node->kernel_boundry->max_lon, neighbor_node->kernel_boundry->max_lon)) {
            tree_node->num_neighbors_on_boundry[PDLN_DOWN]++;
            boundry_type |= PDLN_BOUNDRY_TYPE_D;
            coord_value[0][PDLN_LAT] = coord_value[1][PDLN_LAT] = neighbor_node->kernel_boundry->max_lat;
            coord_value[0][PDLN_LON] = std::max(tree_node->kernel_boundry->min_lon, neighbor_node->kernel_boundry->min_lon);
            coord_value[1][PDLN_LON] = std::min(tree_node->kernel_boundry->max_lon, neighbor_node->kernel_boundry->max_lon);
        }
    }
    else if(tree_node->kernel_boundry->max_lon == neighbor_node->kernel_boundry->min_lon) { // Case 3
        if(std::max(tree_node->kernel_boundry->min_lat, neighbor_node->kernel_boundry->min_lat) <
           std::min(tree_node->kernel_boundry->max_lat, neighbor_node->kernel_boundry->max_lat)) {
            tree_node->num_neighbors_on_boundry[PDLN_RIGHT]++;
            boundry_type |= PDLN_BOUNDRY_TYPE_R;
            coord_value[0][PDLN_LON] = coord_value[1][PDLN_LON] = tree_node->kernel_boundry->max_lon;
            coord_value[0][PDLN_LAT] = std::max(tree_node->kernel_boundry->min_lat, neighbor_node->kernel_boundry->min_lat);
            coord_value[1][PDLN_LAT] = std::min(tree_node->kernel_boundry->max_lat, neighbor_node->kernel_boundry->max_lat);
        }
    }
    *boundry_head = Point(coord_value[0][PDLN_LON], coord_value[0][PDLN_LAT]);
    *boundry_tail = Point(coord_value[1][PDLN_LON], coord_value[1][PDLN_LAT]);

    /* cyclic boundry detecting */
    coord_value[0][PDLN_LAT] = coord_value[0][PDLN_LON] = coord_value[1][PDLN_LAT] = coord_value[1][PDLN_LON] = PDLN_DOUBLE_INVALID_VALUE;
    if(fabs(fabs(tree_node->kernel_boundry->min_lon - neighbor_node->kernel_boundry->max_lon) - 360.0) < PDLN_FLOAT_EQ_ERROR ||
       fabs(fabs(tree_node->kernel_boundry->max_lon - neighbor_node->kernel_boundry->min_lon) - 360.0) < PDLN_FLOAT_EQ_ERROR) { // Case 1 or case 3
        if(std::max(tree_node->kernel_boundry->min_lat, neighbor_node->kernel_boundry->min_lat) < 
           std::min(tree_node->kernel_boundry->max_lat, neighbor_node->kernel_boundry->max_lat)) {
            coord_value[0][PDLN_LON] = coord_value[1][PDLN_LON] = std::min(tree_node->kernel_boundry->min_lon, neighbor_node->kernel_boundry->min_lon);
            coord_value[0][PDLN_LAT] = std::max(tree_node->kernel_boundry->min_lat, neighbor_node->kernel_boundry->min_lat);
            coord_value[1][PDLN_LAT] = std::min(tree_node->kernel_boundry->max_lat, neighbor_node->kernel_boundry->max_lat);
        }
    }
    if(std::max(tree_node->kernel_boundry->min_lat, neighbor_node->kernel_boundry->min_lat) < 
       std::min(tree_node->kernel_boundry->max_lat, neighbor_node->kernel_boundry->max_lat)) {
        if(fabs(fabs(tree_node->kernel_boundry->min_lon - neighbor_node->kernel_boundry->max_lon) - 360.0) < PDLN_FLOAT_EQ_ERROR) { // Case 1
            tree_node->num_neighbors_on_boundry[PDLN_LEFT]++;
            boundry_type |= PDLN_BOUNDRY_TYPE_L;
        }
        if(fabs(fabs(tree_node->kernel_boundry->max_lon - neighbor_node->kernel_boundry->min_lon) - 360.0) < PDLN_FLOAT_EQ_ERROR) { // Case 3
            tree_node->num_neighbors_on_boundry[PDLN_RIGHT]++;
            boundry_type |= PDLN_BOUNDRY_TYPE_R;
        }
    }
    *cyclic_boundry_head = Point(coord_value[0][PDLN_LON], coord_value[0][PDLN_LAT]);
    *cyclic_boundry_tail = Point(coord_value[1][PDLN_LON], coord_value[1][PDLN_LAT]);
    return boundry_type;
}


/* non-block */
void Delaunay_grid_decomposition::send_triangles_to_remote(int src_common_id, int dst_common_id, Triangle_Transport *triangles_buf, int num_triangles, int tag)
{
    MPI_Request request;
    if(processing_info->get_local_process_id() == processing_info->get_processing_unit(dst_common_id)->process_id)
        processing_info->send_to_local_thread(triangles_buf, num_triangles, sizeof(Triangle_Transport),
                                              processing_info->get_processing_unit(src_common_id)->thread_id,
                                              processing_info->get_processing_unit(dst_common_id)->thread_id, tag);
    else {
        MPI_Isend(triangles_buf, num_triangles*sizeof(Triangle_Transport), MPI_CHAR,
                  processing_info->get_processing_unit(dst_common_id)->process_id, 
                  tag, processing_info->get_mpi_comm(), &request);
    }
}


/* block */
int Delaunay_grid_decomposition::recv_triangles_from_remote(int src_common_id, int dst_common_id, Triangle_Transport *triangles_buf, int num_max_triangles, int tag)
{
    if(processing_info->get_local_process_id() == processing_info->get_processing_unit(src_common_id)->process_id) {
        assert(false);
        return processing_info->recv_from_local_thread(triangles_buf, num_max_triangles, sizeof(Triangle_Transport),
                                                       processing_info->get_processing_unit(src_common_id)->thread_id,
                                                       processing_info->get_processing_unit(dst_common_id)->thread_id, tag) / sizeof(Triangle_Transport);
    }
    else {
        MPI_Status status;
        int count;
        MPI_Recv(triangles_buf, num_max_triangles*sizeof(Triangle_Transport), MPI_CHAR,
                 processing_info->get_processing_unit(src_common_id)->process_id, 
                 tag, processing_info->get_mpi_comm(), &status);
        MPI_Get_count(&status, MPI_CHAR, &count);
        assert(count%sizeof(Triangle_Transport) == 0);
        return count/sizeof(Triangle_Transport);
    }
}


namespace std   
{
    namespace tr1
    {
        template <>
        struct hash<Point>
        {
            std::size_t operator()(const Point &p) const
            {
                return (hash<double>()(p.x) ^ hash<double>()(p.y) << 1) >> 1;
            }
        };

        template <>  
        struct hash<Triangle_Transport>  
        {  
            std::size_t operator()(const Triangle_Transport &t) const  
            {  
                return hash<Point>()(t.v[0]) ^ hash<Point>()(t.v[1]) ^ hash<Point>()(t.v[2]);
            }
        }; 

        template <>  
        struct hash<Triangle_ID_Only>  
        {  
            std::size_t operator()(const Triangle_ID_Only &t) const  
            {  
                return hash<int>()(t.id[0]+t.id[1]+t.id[2]+t.id[0]*t.id[1]*t.id[2]);
            }
        }; 
    }
} 


#define PDLN_SET_TAG_ITER(tag)          ( tag  &0x000000FF)
#define PDLN_SET_TAG_SRC(tag, id)       ((id<<20&0xFFF00000) | tag)
#define PDLN_SET_TAG_DST(tag, id)       ((id<< 8&0x000FFF00) | tag)
#define PDLN_SET_TAG(src, dst, iter)    (PDLN_SET_TAG_DST(PDLN_SET_TAG_SRC(PDLN_SET_TAG_ITER(iter), src), dst))

void Delaunay_grid_decomposition::send_recv_checksums_with_neighbors(Search_tree_node *leaf_node, unsigned *local_checksums,
                                                                     unsigned *remote_checksums, vector<MPI_Request*> *waiting_list, int iter)
{
    /* calculate local checksum and send to neighbor */
    Point common_boundary_head, common_boundary_tail, cyclic_common_boundary_head, cyclic_common_boundary_tail;

    double threshold = std::min(search_tree_root->kernel_boundry->max_lon - search_tree_root->kernel_boundry->min_lon,
                                search_tree_root->kernel_boundry->max_lat - search_tree_root->kernel_boundry->min_lat) /
                       sqrt(processing_info->get_num_total_processing_units()) / 2.0;

    leaf_node->init_num_neighbors_on_boundry(0);
    for(unsigned i = 0; i < leaf_node->neighbors.size(); i++) {
        //printf("[%2d] neighbor ID: %2lu\n", leaf_node->region_id, leaf_node->neighbors[i].first->region_id);
        if(leaf_node->neighbors[i].second)
            continue;

        /* compute shared boundry of leaf_node and its neighbor */
#ifdef DEBUG
        assert(leaf_node->neighbors[i].first->region_ids.size() == 1);
        assert(iter <= 0xFF);
        assert(leaf_node->region_id <= 0xFFF);
        assert(leaf_node->neighbors[i].first->region_id <= 0xFFF);
#endif
        unsigned boundry_type = compute_common_boundry(leaf_node, leaf_node->neighbors[i].first, &common_boundary_head, &common_boundary_tail,
                                              &cyclic_common_boundary_head, &cyclic_common_boundary_tail);
        
        local_checksums[i] = 0;
        if(common_boundary_head.x != PDLN_DOUBLE_INVALID_VALUE) {
            unsigned checksum = leaf_node->triangulation->calculate_triangles_intersected_checksum(common_boundary_head, common_boundary_tail, threshold);
            local_checksums[i] ^= checksum;
        }

        if(cyclic_common_boundary_head.x != PDLN_DOUBLE_INVALID_VALUE) {
            unsigned checksum = leaf_node->triangulation->calculate_triangles_intersected_checksum(cyclic_common_boundary_head, cyclic_common_boundary_tail, threshold);
            local_checksums[i] ^= checksum;
        }
        local_checksums[i] = set_boundry_type(local_checksums[i], boundry_type);

        if(common_boundary_head.x != PDLN_DOUBLE_INVALID_VALUE || cyclic_common_boundary_head.x != PDLN_DOUBLE_INVALID_VALUE) {
            MPI_Request *req = new MPI_Request;
#ifdef DEBUG
            waiting_list->push_back(req);
#endif
            MPI_Isend(&local_checksums[i], 1, MPI_UNSIGNED, processing_info->get_processing_unit(regionID_to_unitID[leaf_node->neighbors[i].first->region_ids[0]])->process_id,
                      PDLN_SET_TAG(leaf_node->region_id, leaf_node->neighbors[i].first->region_id, iter),
                      processing_info->get_mpi_comm(), req); 

            waiting_list->push_back(new MPI_Request);
            MPI_Irecv(&remote_checksums[i], 1, MPI_UNSIGNED, processing_info->get_processing_unit(regionID_to_unitID[leaf_node->neighbors[i].first->region_ids[0]])->process_id,
                      PDLN_SET_TAG(leaf_node->neighbors[i].first->region_id, leaf_node->region_id, iter),
                      processing_info->get_mpi_comm(), waiting_list->back());
        }
        else
            remote_checksums[i] = 0;
    }
}


bool Delaunay_grid_decomposition::are_checksums_identical(Search_tree_node *leaf_node, unsigned *local_checksums, unsigned *remote_checksums)
{
    if(leaf_node->neighbors.size() == 0) {
        printf("some region has no neighbor, that's weird\n");
        return false;
    }

    bool ok = true;
    for(unsigned i = 0; i < leaf_node->neighbors.size(); i++) {
        if(leaf_node->neighbors[i].second) {
            /* local_checksums[i] stay unchanged if the i-th neighbor has passed the check,
             * so we can just get the boundry type from local_checksums[i] */
            leaf_node->clear_expanding_count(get_boundry_type(local_checksums[i]));
            //printf("[%d] neighbor %d already done\n", leaf_node->region_id, leaf_node->neighbors[i].first->region_id);
            continue;
        }

        if((local_checksums[i] & PDLN_BOUNDRY_TYPE_CLEAR) == (remote_checksums[i] & PDLN_BOUNDRY_TYPE_CLEAR)) {
            //printf("[%d] neighbor %d done, %x vs %x\n", leaf_node->region_id, leaf_node->neighbors[i].first->region_id, local_checksums[i], remote_checksums[i]);
            leaf_node->neighbors[i].second = true;
            leaf_node->reduce_num_neighbors_on_boundry(get_boundry_type(local_checksums[i]));
            leaf_node->clear_expanding_count(get_boundry_type(local_checksums[i]));
        }
        else {
            //printf("[%d] neighbor %d not , %x vs %x\n", leaf_node->region_id, leaf_node->neighbors[i].first->region_id, local_checksums[i], remote_checksums[i]);
            ok = false;
        }
    }

    return ok;
}


int Delaunay_grid_decomposition::assign_polars(bool assign_south_polar, bool assign_north_polar)
{
    double*     c_points_coord[4];
    int*        c_points_index[2];
    Boundry     c_boundry[2];
    int         c_num_points[2];
    vector<int> c_regions_id[2];
    Midline     midline;

    if(!assign_south_polar && !assign_north_polar)
        return 0;
    
    for(int i = 0; i < 4; i++)
        c_points_coord[i] = new double[search_tree_root->num_kernel_points];
    c_points_index[0] = new int[search_tree_root->num_kernel_points];
    c_points_index[1] = new int[search_tree_root->num_kernel_points];

    if(assign_south_polar) {
        //printf("South polar need rotating.\n");
        current_tree_node->decompose_by_processing_units_number(workloads, c_points_coord, c_points_index, c_num_points, c_boundry, c_regions_id, PDLN_DECOMPOSE_SPOLAR_MODE);
        if(c_boundry[0].max_lat > PDLN_SPOLAR_MAX_LAT) {
            midline.type = PDLN_LAT;
            midline.value = PDLN_SPOLAR_MAX_LAT;
            search_tree_root->split_local_points(midline, c_points_coord, c_points_index, c_num_points);

            if(c_num_points[0] < min_points_per_chunk)
                goto fail;

            c_boundry[0].max_lat = PDLN_SPOLAR_MAX_LAT;
            c_boundry[1].min_lat = PDLN_SPOLAR_MAX_LAT;

            /* 0 is preserved for this situation */
            c_regions_id[0].clear();
            c_regions_id[0].push_back(0);
            workloads[1] -= c_num_points[0];
            c_regions_id[1].insert(c_regions_id[1].begin(), 1);
        }
        search_tree_root->children[0] = alloc_search_tree_node(search_tree_root, c_points_coord,   c_points_index[0], c_num_points[0],
                                                               c_boundry[0], c_regions_id[0], PDLN_NODE_TYPE_SPOLAR);
        search_tree_root->children[1] = alloc_search_tree_node(search_tree_root, c_points_coord+2, c_points_index[1], c_num_points[1],
                                                               c_boundry[1], c_regions_id[1], PDLN_NODE_TYPE_COMMON);

        current_tree_node = search_tree_root->children[1];
        
        if(have_local_region_ids(search_tree_root->children[0]->region_ids))
            local_leaf_nodes.push_back(search_tree_root->children[0]);
        all_leaf_nodes.push_back(search_tree_root->children[0]);

        /* multiple polars are shifting only on polar node, points of search_tree_root won't change */
        double shifted_polar_lat = search_tree_root->children[0]->load_polars_info();
        if(shifted_polar_lat != PDLN_DOUBLE_INVALID_VALUE)
            search_tree_root->real_boundry->min_lat = shifted_polar_lat;
    }
    
    if(assign_north_polar) {
        //printf("Nouth polar need rotating.\n");
        current_tree_node->decompose_by_processing_units_number(workloads, c_points_coord, c_points_index, c_num_points, c_boundry, c_regions_id, PDLN_DECOMPOSE_NPOLAR_MODE);
        if(c_boundry[1].min_lat < PDLN_NPOLAR_MIN_LAT) {
            midline.type = PDLN_LAT;
            midline.value = PDLN_NPOLAR_MIN_LAT;
            current_tree_node->split_local_points(midline, c_points_coord, c_points_index, c_num_points);

            if(c_num_points[1] < min_points_per_chunk)
                goto fail;

            c_boundry[0].max_lat = PDLN_NPOLAR_MIN_LAT;
            c_boundry[1].min_lat = PDLN_NPOLAR_MIN_LAT;

            /* num_regions+1 is preserved for this situation */
            c_regions_id[1].clear();
            c_regions_id[1].push_back(num_regions+1);
            workloads[num_regions] -= c_num_points[1];
            c_regions_id[0].push_back(num_regions);
        }
        if(search_tree_root->children[1] != NULL)
            delete search_tree_root->children[1];

        search_tree_root->children[2] = alloc_search_tree_node(search_tree_root, c_points_coord+2, c_points_index[1], c_num_points[1], c_boundry[1],
                                                               c_regions_id[1], PDLN_NODE_TYPE_NPOLAR);

        search_tree_root->children[1] = alloc_search_tree_node(search_tree_root, c_points_coord,   c_points_index[0], c_num_points[0], c_boundry[0],
                                                               c_regions_id[0], PDLN_NODE_TYPE_COMMON);

        current_tree_node = search_tree_root->children[1];

        if(have_local_region_ids(search_tree_root->children[2]->region_ids))
            local_leaf_nodes.push_back(search_tree_root->children[2]);
        all_leaf_nodes.push_back(search_tree_root->children[2]);
        /* multiple polars are shifting only on polar node, points of search_tree_root won't change */
        double shifted_polar_lat = search_tree_root->children[2]->load_polars_info();
        if(shifted_polar_lat != PDLN_DOUBLE_INVALID_VALUE)
            search_tree_root->real_boundry->max_lat = shifted_polar_lat;
    }

    for(int i = 0; i < 4; i++)
        delete[] c_points_coord[i];
    delete[] c_points_index[0];
    delete[] c_points_index[1];
    return 0;

fail:
    for(int i = 0; i < 4; i++)
        delete[] c_points_coord[i];
    delete[] c_points_index[0];
    delete[] c_points_index[1];
    printf("assign polars fault, %d, %d vs %d\n", c_num_points[0], c_num_points[1], min_points_per_chunk);
    return 1;
}


void Delaunay_grid_decomposition::decompose_with_fixed_longitude(double fixed_lon)
{
    double*     c_points_coord[4];
    int*        c_points_index[2];
    Boundry     c_boundry[2];
    int         c_num_points[2];
    vector<int> c_regions_id[2];

    for(int i = 0; i < 4; i++)
        c_points_coord[i] = new double[current_tree_node->num_kernel_points];
    c_points_index[0] = new int[search_tree_root->num_kernel_points];
    c_points_index[1] = new int[search_tree_root->num_kernel_points];

    current_tree_node->decompose_by_fixed_longitude(fixed_lon, workloads, c_points_coord, c_points_index, c_num_points, c_boundry, c_regions_id);

    current_tree_node->children[0] = alloc_search_tree_node(current_tree_node, c_points_coord,   c_points_index[0], c_num_points[0], c_boundry[0],
                                                            c_regions_id[0], PDLN_NODE_TYPE_COMMON);

    current_tree_node->children[2] = alloc_search_tree_node(current_tree_node, c_points_coord+2, c_points_index[1], c_num_points[1], c_boundry[1],
                                                            c_regions_id[1], PDLN_NODE_TYPE_COMMON);

    for(int i = 0; i < 4; i++)
        delete[] c_points_coord[i];
    delete[] c_points_index[0];
    delete[] c_points_index[1];
}


//TODO: get faster
bool Delaunay_grid_decomposition::have_local_region_ids(vector<int> chunk_id)
{
    for(int j = 0; j < processing_info->get_num_local_proc_processing_units(); j++)
        for(unsigned i = 0; i < chunk_id.size(); i++)
            if(regionID_to_unitID[chunk_id[i]] == processing_info->get_local_proc_common_id()[j])
                return true;

    return false;
}


int Delaunay_grid_decomposition::generate_grid_decomposition(bool lazy_mode)
{
    initialze_workload();
    initialze_buffer();
    current_tree_node = search_tree_root;

    double min_lon, max_lon, min_lat, max_lat;
    grid_info_mgr->get_grid_boundry(original_grid, &min_lon, &max_lon, &min_lat, &max_lat);

    if(assign_polars(std::abs(min_lat - -90.0) < PDLN_FLOAT_EQ_ERROR, std::abs(max_lat -  90.0) < PDLN_FLOAT_EQ_ERROR))
        return 1;

    int num_computing_nodes = processing_info->get_num_computing_nodes();
    Processing_unit** units = processing_info->get_processing_units();

    delete all_group_intervals;
    all_group_intervals    = new int[num_computing_nodes]();
    unsigned old_checksum  = units[regionID_to_unitID[current_tree_node->region_ids[0]]]->hostname_checksum;
    int cur_group          = 0;
    all_group_intervals[0] = 1;
    for(int i = 1; i < current_tree_node->region_ids.size(); i++)
        if (old_checksum == units[regionID_to_unitID[current_tree_node->region_ids[i]]]->hostname_checksum)
            all_group_intervals[cur_group]++;
        else {
            all_group_intervals[++cur_group]++;
            old_checksum = units[regionID_to_unitID[current_tree_node->region_ids[i]]]->hostname_checksum;
        }
    assert(cur_group+1 <= num_computing_nodes);

    current_tree_node->set_groups(all_group_intervals, cur_group+1);
    decompose_common_node_recursively(current_tree_node, lazy_mode);
    return 0;
}


double Search_tree_node::load_polars_info()
{
    polars_local_index = new vector<int>();

    if(node_type == PDLN_NODE_TYPE_NPOLAR) {
        double nearest_point_lat = -1e10;
        for(int i = 0; i < num_kernel_points+num_expanded_points; i++) {
            if(std::abs(points_coord[PDLN_LAT][i] - 90.0) < PDLN_FLOAT_EQ_ERROR)
                polars_local_index->push_back(i);
            else if(nearest_point_lat < points_coord[PDLN_LAT][i])
                nearest_point_lat = points_coord[PDLN_LAT][i];
        }

        if(PDLN_INSERT_VIRTUAL_POINT && polars_local_index->size() != 1) {
            /* Note: only polar nodes' coord value is changed, search tree root's coord value is unchanged */
            for(unsigned i = 0; i < polars_local_index->size(); i++)
                points_coord[PDLN_LAT][(*polars_local_index)[i]] = (90.0 + nearest_point_lat) * 0.5;

            double vpoint_lon = 0;
            double vpoint_lat = 90;
            int vpoint_idx = -1;
            virtual_point_local_index = num_kernel_points+num_expanded_points;
            add_expanded_points(&vpoint_lon, &vpoint_lat, &vpoint_idx, 1);
#ifdef DEBUG
            assert(points_coord[PDLN_LON][virtual_point_local_index] == vpoint_lon);
            assert(points_coord[PDLN_LAT][virtual_point_local_index] == vpoint_lat);
            assert(points_global_index[virtual_point_local_index] == vpoint_idx);
#endif
            return (90.0 + nearest_point_lat) * 0.5;
        }
        return PDLN_DOUBLE_INVALID_VALUE;
    }
    else if(node_type == PDLN_NODE_TYPE_SPOLAR) {
        double nearest_point_lat = 1e10;
        for(int i = 0; i < num_kernel_points+num_expanded_points; i++) {
            if(std::abs(points_coord[PDLN_LAT][i] - -90.0) < PDLN_FLOAT_EQ_ERROR)
                polars_local_index->push_back(i);
            else if(nearest_point_lat > points_coord[PDLN_LAT][i])
                nearest_point_lat = points_coord[PDLN_LAT][i];
        }

        if(PDLN_INSERT_VIRTUAL_POINT && polars_local_index->size() != 1) {
            for(unsigned i = 0; i < polars_local_index->size(); i++)
                points_coord[PDLN_LAT][(*polars_local_index)[i]] = (-90.0 + nearest_point_lat) * 0.5;

            double vpoint_lon = 0;
            double vpoint_lat = -90;
            int vpoint_idx = -1;
            virtual_point_local_index = num_kernel_points+num_expanded_points;
            add_expanded_points(&vpoint_lon, &vpoint_lat, &vpoint_idx, 1);
#ifdef DEBUG
            assert(points_coord[PDLN_LON][virtual_point_local_index] == vpoint_lon);
            assert(points_coord[PDLN_LAT][virtual_point_local_index] == vpoint_lat);
            assert(points_global_index[virtual_point_local_index] == vpoint_idx);
#endif
            return (-90.0 + nearest_point_lat) * 0.5;
        }
        return PDLN_DOUBLE_INVALID_VALUE;
    }
    return PDLN_DOUBLE_INVALID_VALUE;
}


Boundry Search_tree_node::expand()
{
    Boundry expanded = *expanded_boundry;

    if(node_type == PDLN_NODE_TYPE_SPOLAR) {
        if(num_neighbors_on_boundry[PDLN_UP] > 0) expanded.max_lat += (expanded.max_lat - expanded.min_lat) * expanding_scale[PDLN_UP]++ / 10.;
        return expanded;
    }
    if(node_type == PDLN_NODE_TYPE_NPOLAR) {
        if(num_neighbors_on_boundry[PDLN_DOWN] > 0) expanded.min_lat -= (expanded.max_lat - expanded.min_lat) * expanding_scale[PDLN_DOWN]++ / 10.;
        return expanded;
    }
    if(node_type == PDLN_NODE_TYPE_COMMON) {
        if(num_neighbors_on_boundry[PDLN_UP] > 0) {
            expanded.max_lat += (expanded.max_lat - expanded.min_lat) * expanding_scale[PDLN_UP]++ / 10.;
            edge_expanding_count[PDLN_UP]++;
        }
        if(num_neighbors_on_boundry[PDLN_LEFT] > 0) { 
            expanded.min_lon -= (expanded.max_lon - expanded.min_lon) * expanding_scale[PDLN_LEFT]++ / 10.;
            edge_expanding_count[PDLN_LEFT]++;
        }
        if(num_neighbors_on_boundry[PDLN_DOWN] > 0) {
            expanded.min_lat -= (expanded.max_lat - expanded.min_lat) * expanding_scale[PDLN_DOWN]++ / 10.;
            edge_expanding_count[PDLN_DOWN]++;
        }
        if(num_neighbors_on_boundry[PDLN_RIGHT] > 0) {
            expanded.max_lon += (expanded.max_lon - expanded.min_lon) * expanding_scale[PDLN_RIGHT]++ / 10.;
            edge_expanding_count[PDLN_RIGHT]++;
        }
    }
    return expanded;
}


bool Search_tree_node::expanding_success(Boundry *prev_expanded, Boundry *max_expanded, bool is_cyclic)
{
    if(expanded_boundry->min_lat < prev_expanded->min_lat || expanded_boundry->min_lat <= max_expanded->min_lat)
        num_neighbors_on_boundry[PDLN_DOWN] = 0;
    else {
        //printf("down not\n");
        //num_neighbors_on_boundry[PDLN_DOWN] = 1;
    }
    if(expanded_boundry->max_lat > prev_expanded->max_lat || expanded_boundry->max_lat >= max_expanded->max_lat)
        num_neighbors_on_boundry[PDLN_UP] = 0;
    else {
        //printf("up not\n");
        //num_neighbors_on_boundry[PDLN_UP] = 1;
    }
    if(node_type != PDLN_NODE_TYPE_COMMON || expanded_boundry->min_lon < prev_expanded->min_lon || (!is_cyclic && expanded_boundry->min_lon <= max_expanded->min_lon))
        num_neighbors_on_boundry[PDLN_LEFT] = 0;
    else {
        //printf("left not\n");
        //num_neighbors_on_boundry[PDLN_LEFT] = 1;
    }
    if(node_type != PDLN_NODE_TYPE_COMMON || expanded_boundry->max_lon > prev_expanded->max_lon || (!is_cyclic && expanded_boundry->max_lon >= max_expanded->max_lon))
        num_neighbors_on_boundry[PDLN_RIGHT] = 0;
    else {
        //printf("right not\n");
        //num_neighbors_on_boundry[PDLN_RIGHT] = 1;
    }
    return num_neighbors_on_boundry[PDLN_DOWN] == 0 && num_neighbors_on_boundry[PDLN_UP] == 0 &&
           num_neighbors_on_boundry[PDLN_LEFT] == 0 && num_neighbors_on_boundry[PDLN_RIGHT] == 0;
}


void Search_tree_node::set_groups(int *intervals, int num)
{
    group_intervals = intervals;
    num_groups = num;
}


#define PDLN_MIN_QUOTA (400.0)
int Delaunay_grid_decomposition::expand_tree_node_boundry(Search_tree_node* tree_node, double expanding_ratio)
{
    double *expanded_coord[2];
    int *expanded_index;
    int num_found;

    expanded_coord[0] = new double[search_tree_root->num_kernel_points]; //FIXME: buf too large
    expanded_coord[1] = new double[search_tree_root->num_kernel_points];
    expanded_index    = new int[search_tree_root->num_kernel_points];

    for (int i = 0; i < 4; i++)
        if (tree_node->edge_expanding_count[i] > 2)
            tree_node->num_neighbors_on_boundry[i+1] = tree_node->num_neighbors_on_boundry[(i+3)%4] = 1;

    int num_edge_to_expand = 0;
    for (int i = 0; i < 4; i++)
        if (tree_node->num_neighbors_on_boundry[i] > 0)
            num_edge_to_expand++;

    assert(num_edge_to_expand > 0);

    double quota = std::max(average_workload * expanding_ratio, PDLN_MIN_QUOTA) / num_edge_to_expand;
    Boundry old_boundry, new_boundry;
    old_boundry = *tree_node->expanded_boundry;
    int fail_count = 0;
    do {
        Boundry last_boundry = *tree_node->expanded_boundry;
        new_boundry = tree_node->expand();
        new_boundry.legalize(search_tree_root->kernel_boundry, is_cyclic);

        //printf("kern boundary: %lf, %lf, %lf, %lf\n", tree_node->kernel_boundry->min_lon, tree_node->kernel_boundry->max_lon, tree_node->kernel_boundry->min_lat, tree_node->kernel_boundry->max_lat);
        //printf("before adjust: %lf, %lf, %lf, %lf\n", new_boundry.min_lon, new_boundry.max_lon, new_boundry.min_lat, new_boundry.max_lat);

        //printf("adjusting ID: %d\n", tree_node->region_id);
        adjust_expanding_boundry(&old_boundry, &new_boundry, quota, expanded_coord, expanded_index, &num_found);
        //printf("after  adjust: %lf, %lf, %lf, %lf\n", new_boundry.min_lon, new_boundry.max_lon, new_boundry.min_lat, new_boundry.max_lat);
        //add_halo_points(tree_node, &old_boundry, &new_boundry);
        //printf("expanded boundry : %lf, %lf, %lf, %lf\n", tree_node->expanded_boundry->min_lon, tree_node->expanded_boundry->max_lon, tree_node->expanded_boundry->min_lat, tree_node->expanded_boundry->max_lat);
        //printf("root real boundry: %lf, %lf, %lf, %lf\n", search_tree_root->real_boundry->min_lon, search_tree_root->real_boundry->max_lon, search_tree_root->real_boundry->min_lat, search_tree_root->real_boundry->max_lat);
        if(last_boundry == *tree_node->expanded_boundry)
            fail_count ++;
        else
            fail_count = 0;
        if(fail_count > 5) {
            printf("expanding failed, max3\n");
            return -1;
        }
        if(processing_info->get_num_total_processing_units() > 4 &&
           new_boundry.max_lon - new_boundry.min_lon > (search_tree_root->kernel_boundry->max_lon - search_tree_root->kernel_boundry->min_lon) * 0.75 &&
           new_boundry.max_lat - new_boundry.min_lat > (search_tree_root->kernel_boundry->max_lat - search_tree_root->kernel_boundry->min_lat) * 0.75) {
            printf("expanded to the max1\n");
            return -1;
        }
        if(new_boundry == *search_tree_root->kernel_boundry || new_boundry.max_lon - new_boundry.min_lon > 360.0) {
            printf("expanded to the max2\n");
            return -1;
        }
        *tree_node->expanded_boundry = new_boundry;
    }while(!tree_node->expanding_success(&old_boundry, search_tree_root->real_boundry, is_cyclic));

    vector<Search_tree_node*> leaf_nodes_found = search_halo_points_from_top(&old_boundry, &new_boundry, expanded_coord, expanded_index, &num_found);
    tree_node->add_expanded_points(expanded_coord, expanded_index, num_found);
    tree_node->add_neighbors(leaf_nodes_found);

    delete[] expanded_coord[0];
    delete[] expanded_coord[1];
    delete[] expanded_index;

    return 0;
}


/*
 *  ┌───┬───────┐
 *  │ 0 │     3 │
 *  │   ├───┬───┤
 *  │   │   │   │
 *  ├───┴───┤   │
 *  │ 1     │ 2 │
 *  └───────┴───┘
 */
void Delaunay_grid_decomposition::halo_to_rectangles(Boundry inner_boundry, Boundry outer_boundry, Boundry sub_rectangle[4])
{
    sub_rectangle[0].min_lon = sub_rectangle[1].min_lon = outer_boundry.min_lon;
    sub_rectangle[1].min_lat = sub_rectangle[2].min_lat = outer_boundry.min_lat;
    sub_rectangle[2].max_lon = sub_rectangle[3].max_lon = outer_boundry.max_lon;
    sub_rectangle[3].max_lat = sub_rectangle[0].max_lat = outer_boundry.max_lat;

    sub_rectangle[0].min_lat = sub_rectangle[1].max_lat = inner_boundry.min_lat;
    sub_rectangle[1].max_lon = sub_rectangle[2].min_lon = inner_boundry.max_lon;
    sub_rectangle[2].max_lat = sub_rectangle[3].min_lat = inner_boundry.max_lat;
    sub_rectangle[3].min_lon = sub_rectangle[0].max_lon = inner_boundry.min_lon;
}


void Delaunay_grid_decomposition::rectangles_to_halo(Boundry sub_rectangle[4], Boundry* outer_boundry)
{
    outer_boundry->min_lon = sub_rectangle[0].min_lon;
    outer_boundry->min_lat = sub_rectangle[1].min_lat;
    outer_boundry->max_lon = sub_rectangle[2].max_lon;
    outer_boundry->max_lat = sub_rectangle[3].max_lat;
}


static inline bool is_in_region(double x, double y, Boundry region)
{
    return x >= region.min_lon && x < region.max_lon && y >= region.min_lat && y < region.max_lat;
}


int Delaunay_grid_decomposition::classify_points(double *coord[2], int *index, int num, Boundry region, int start)
{
    int j = start;
    for (int i = start; i < num; i++)
        if (is_in_region(coord[PDLN_LON][i], coord[PDLN_LAT][i], region)) {
            std::swap(coord[PDLN_LON][i], coord[PDLN_LON][j]);
            std::swap(coord[PDLN_LAT][i], coord[PDLN_LAT][j]);
            std::swap(index[i], index[j++]);
        }
    return j - start;
}


#define PDLN_SAVE_N (0)
#define PDLN_SAVE_L (1)
#define PDLN_SAVE_R (2)
double Delaunay_grid_decomposition::adjust_subrectangle(double l, double r, double *coord[2], int *idx, int offset,
                                                        int num, Boundry *bound, int linetype, int save_option)
{
    int c_num_points[2];
    double length[2], boundry_values[4];
    Midline midline;

    if (bound->min_lon == bound->max_lon || bound->min_lat == bound->max_lat)
        return 0.0;

    assert(l > 0);
    assert(r > 0);

    boundry_values[PDLN_LON] = bound->min_lon;
    boundry_values[PDLN_LAT] = bound->min_lat;
    boundry_values[PDLN_LON+2] = bound->max_lon;
    boundry_values[PDLN_LAT+2] = bound->max_lat;

    assert(boundry_values[PDLN_LON] != boundry_values[PDLN_LON+2]);
    assert(boundry_values[PDLN_LAT] < boundry_values[PDLN_LAT+2]);
    length[0] = boundry_values[PDLN_LON+2] - boundry_values[PDLN_LON];
    length[1] = boundry_values[PDLN_LAT+2] - boundry_values[PDLN_LAT];
    //printf("l: %lf, r: %lf, total: %d, [%lf, %lf],[%lf, %lf]\n", l, r, num, boundry_values[PDLN_LON], boundry_values[PDLN_LON+2], boundry_values[PDLN_LAT], boundry_values[PDLN_LAT+2]);
    midline.type = linetype;
    midline.value = boundry_values[midline.type] + length[midline.type] * l / (l + r);
    Search_tree_node::count_points(coord, idx, offset, num, midline, c_num_points);
    //printf("midline.value: %lf, (%d, %d)\n", midline.value, c_num_points[0], c_num_points[1]);

    int iteration_count = 1;
    while (c_num_points[0] == 0 || c_num_points[1] == 0 ||
           fabs(c_num_points[0]/c_num_points[1] - l/r) > PDLN_TOLERABLE_ERROR) {
        //printf("loop: %d\n", iteration_count);

        if (save_option == PDLN_SAVE_L) {
            if(iteration_count++ > PDLN_MAX_ITER_COUNT && c_num_points[0] > 0)
                break;
        } else {
            if(iteration_count++ > PDLN_MAX_ITER_COUNT && c_num_points[1] > 0)
                break;
        }

        if(c_num_points[0] < l) {
#ifdef DEBUG
            assert(c_num_points[1] >= r || fabs(r - c_num_points[1]) < PDLN_FLOAT_EQ_ERROR);
#endif
            midline.value += (boundry_values[2+midline.type] - midline.value) * (c_num_points[1] - r) / c_num_points[1];
        }
        else {
#ifdef DEBUG
            assert(c_num_points[1] <= r || fabs(r - c_num_points[1]) < PDLN_FLOAT_EQ_ERROR);
#endif
            midline.value -= (midline.value - boundry_values[midline.type]) * (c_num_points[0] - l) / c_num_points[0];
        }
        //printf("midline.value: %lf, (%d, %d)\n", midline.value, c_num_points[0], c_num_points[1]);
        assert(midline.value > boundry_values[midline.type]);
        assert(midline.value < boundry_values[2+midline.type]);
        /* TODO: Search only half of the whole points, but not the whole points */
        Search_tree_node::count_points(coord, idx, offset, num, midline, c_num_points);
    }
    if (save_option == PDLN_SAVE_L) {
        assert(c_num_points[0] > 0);
        if (linetype == PDLN_LON) bound->max_lon = midline.value;
        else bound->max_lat = midline.value;
    } else {
        assert(c_num_points[1] > 0);
        if (linetype == PDLN_LON) bound->min_lon = midline.value;
        else bound->min_lat = midline.value;
    }
    return midline.value;
}


vector<Search_tree_node*> Delaunay_grid_decomposition::adjust_expanding_boundry(const Boundry* inner, Boundry* outer, double quota,
                                                           double *expanded_coord[2], int *expanded_index, int *num_found)
{
    vector<Search_tree_node*> leaf_nodes_found = search_halo_points_from_top(inner, outer, expanded_coord, expanded_index, num_found);
    int squeeze_count = 0;
    Boundry old_outer = *outer;
    Boundry sub_rectangles[4];
    halo_to_rectangles(*inner, *outer, sub_rectangles);
    //for (int i = 0; i < 4; i++)
    //    printf("sub rectangles%d [%lf, %lf], [%lf, %lf]\n", i, sub_rectangles[i].min_lon, sub_rectangles[i].max_lon, sub_rectangles[i].min_lat, sub_rectangles[i].max_lat);
    int sorted = 0;
    for (int i = 0; i < 4; i++) {
        int num = classify_points(expanded_coord, expanded_index, *num_found, sub_rectangles[i], sorted);
        if (num <= 0 && sub_rectangles[i].min_lat != sub_rectangles[i].max_lat && sub_rectangles[i].min_lon != sub_rectangles[i].max_lon)
            assert(false);
        if (num <= 0 || num <= quota)
            continue;
        int linetype = i%2 ? PDLN_LAT : PDLN_LON;
        double l_num = i<2 ? num - quota : quota;
        double r_num = i<2 ? quota : num - quota;
        int savetype = i<2 ? PDLN_SAVE_R : PDLN_SAVE_L;
        adjust_subrectangle(l_num, r_num, expanded_coord, expanded_index, sorted, num, &sub_rectangles[i], linetype, savetype);
        sorted += num;
    }
    //for (int i = 0; i < 4; i++)
    //    printf("sub rectangles%d [%lf, %lf], [%lf, %lf]\n", i, sub_rectangles[i].min_lon, sub_rectangles[i].max_lon, sub_rectangles[i].min_lat, sub_rectangles[i].max_lat);
    rectangles_to_halo(sub_rectangles, outer);
    //printf("fnl [%lf, %lf], [%lf, %lf]\n", outer->min_lon, outer->max_lon, outer->min_lat, outer->max_lat);
    return leaf_nodes_found;
}


void Delaunay_grid_decomposition::add_halo_points(Search_tree_node* dst_tree_node, Boundry* inner_boundary, Boundry* outer_boundary)
{
    double *expanded_points_coord[2];
    int *expanded_index;
    int num_found;

    expanded_points_coord[0] = new double[search_tree_root->num_kernel_points]; //FIXME: buf too large
    expanded_points_coord[1] = new double[search_tree_root->num_kernel_points];
    expanded_index           = new int[search_tree_root->num_kernel_points];

    vector<Search_tree_node*> leaf_nodes_found = search_halo_points_from_top(inner_boundary, outer_boundary, expanded_points_coord, expanded_index, &num_found);
    dst_tree_node->add_expanded_points(expanded_points_coord, expanded_index, num_found);
    dst_tree_node->add_neighbors(leaf_nodes_found);

    delete[] expanded_points_coord[0];
    delete[] expanded_points_coord[1];
    delete[] expanded_index;
}


vector<Search_tree_node*> Delaunay_grid_decomposition::search_halo_points_from_top(const Boundry* inner_boundary, const Boundry* outer_boundary,
                                                                             double *coord_values[2], int *global_idx,
                                                                             int *num_found)
{
    vector<Search_tree_node*> leaf_nodes_found;

    *num_found = 0;

    if(*inner_boundary == *outer_boundary)
        return leaf_nodes_found;

    search_down_for_points_in_halo(search_tree_root, inner_boundary, outer_boundary, leaf_nodes_found, coord_values, global_idx, num_found);

    return leaf_nodes_found;
}


void Delaunay_grid_decomposition::search_halo_points_from_buf(Boundry* inner_boundary, Boundry* outer_boundary,
                                                              double *coord_values[2], int *global_idx, int *num_found)
{
    if(*inner_boundary == *outer_boundary)
        return;

    int num_points = *num_found;
    *num_found = 0;
    Search_tree_node::search_points_in_halo(inner_boundary, outer_boundary, coord_values, global_idx, num_points, coord_values, global_idx, num_found);
}


void Delaunay_grid_decomposition::search_down_for_points_in_halo(Search_tree_node *node, const Boundry *inner_boundary,
                                                                 const Boundry *outer_boundary, vector<Search_tree_node*> &leaf_nodes_found,
                                                                 double *coord_values[2], int *global_idx,
                                                                 int *num_found)
{
    double*     c_points_coord[4];
    int*        c_points_index[2];
    int         c_num_points[2];
    Boundry     c_boundry[2];
    vector<int> c_regions_id[2];
    int*        c_intervals[2];
    int         c_num_intervals[2];

    assert(node->region_ids.size() > 0);

    Boundry region = *outer_boundary;
    if(node->region_ids.size() == 1) {
        if(do_two_regions_overlap(region, *node->kernel_boundry) ||
           do_two_regions_overlap(Boundry(region.min_lon + 360.0, region.max_lon + 360.0, region.min_lat, region.max_lat), *node->kernel_boundry) ||
           do_two_regions_overlap(Boundry(region.min_lon - 360.0, region.max_lon - 360.0, region.min_lat, region.max_lat), *node->kernel_boundry)) {
            leaf_nodes_found.push_back(node);
            node->search_points_in_halo(inner_boundary, outer_boundary, coord_values, global_idx, num_found);
        }
        all_leaf_nodes.push_back(node);
        return;
    }

    if(node->children[0] == NULL && node->children[1] == NULL && node->children[2] == NULL) {
        for(int i = 0; i < 4; i++)
            c_points_coord[i] = new double[node->num_kernel_points];
        c_points_index[0] = new int[node->num_kernel_points];
        c_points_index[1] = new int[node->num_kernel_points];
        
        node->decompose_by_processing_units_number(workloads, c_points_coord, c_points_index, 
                                                   c_num_points, c_boundry, c_regions_id, 
                                                   PDLN_DECOMPOSE_COMMON_MODE, c_intervals,
                                                   c_num_intervals);
        assert(c_regions_id[0].size() > 0);

        node->children[0] = alloc_search_tree_node(node, c_points_coord,   c_points_index[0], c_num_points[0],
                                                   c_boundry[0], c_regions_id[0], PDLN_NODE_TYPE_COMMON);

        node->children[2] = alloc_search_tree_node(node, c_points_coord+2, c_points_index[1], c_num_points[1],
                                                   c_boundry[1], c_regions_id[1], PDLN_NODE_TYPE_COMMON);
        node->children[0]->set_groups(c_intervals[0], c_num_intervals[0]);
        node->children[2]->set_groups(c_intervals[1], c_num_intervals[1]);

        assert(node->children[0]->region_ids.size() > 0);
        assert(node->children[2]->region_ids.size() > 0);

        for(int i = 0; i < 4; i++)
            delete[] c_points_coord[i];
        delete[] c_points_index[0];
        delete[] c_points_index[1];
        //TODO: optimize new delete
    }

    for(int i = 0; i < 3; i ++)
        if(node->children[i] != NULL) {
            if(do_two_regions_overlap(region, *node->children[i]->kernel_boundry) ||
               do_two_regions_overlap(Boundry(region.min_lon + 360.0, region.max_lon + 360.0, region.min_lat, region.max_lat), *node->children[i]->kernel_boundry) ||
               do_two_regions_overlap(Boundry(region.min_lon - 360.0, region.max_lon - 360.0, region.min_lat, region.max_lat), *node->children[i]->kernel_boundry)) {
                search_down_for_points_in_halo(node->children[i], inner_boundary, outer_boundary, leaf_nodes_found, coord_values, global_idx, num_found);
            }
        }
}


void Search_tree_node::search_points_in_halo(const Boundry *inner_boundary, const Boundry *outer_boundary,
                                             double *const coord[2], const int *idx, int num_points,
                                             double *coord_values[2], int *global_idx, int *num_found)
{
    Boundry l_inner = *inner_boundary;
    Boundry l_outer = *outer_boundary;
    Boundry r_inner = *inner_boundary;
    Boundry r_outer = *outer_boundary;
    l_inner.min_lon -= 360.0;
    l_inner.max_lon -= 360.0;
    l_outer.min_lon -= 360.0;
    l_outer.max_lon -= 360.0;
    r_inner.min_lon += 360.0;
    r_inner.max_lon += 360.0;
    r_outer.min_lon += 360.0;
    r_outer.max_lon += 360.0;

    #pragma omp parallel for
    for(int j = 0; j < num_points; j++) {
        if (is_coordinate_in_halo(coord[PDLN_LON][j], coord[PDLN_LAT][j], inner_boundary, outer_boundary)) {
            coord_values[PDLN_LON][*num_found] = coord[PDLN_LON][j];
            coord_values[PDLN_LAT][*num_found] = coord[PDLN_LAT][j];
            global_idx[(*num_found)++] = idx[j]; //FIXME: there may be data race
            continue;
        }
        if (is_coordinate_in_halo(coord[PDLN_LON][j], coord[PDLN_LAT][j], &l_inner, &l_outer)) {
            coord_values[PDLN_LON][*num_found] = coord[PDLN_LON][j] + 360.0;
            coord_values[PDLN_LAT][*num_found] = coord[PDLN_LAT][j];
            global_idx[(*num_found)++] = idx[j];
            continue;
        }
        if (is_coordinate_in_halo(coord[PDLN_LON][j], coord[PDLN_LAT][j], &r_inner, &r_outer)) {
            coord_values[PDLN_LON][*num_found] = coord[PDLN_LON][j] - 360.0;
            coord_values[PDLN_LAT][*num_found] = coord[PDLN_LAT][j];
            global_idx[(*num_found)++] = idx[j];
            continue;
        }
    }
}


void Search_tree_node::search_points_in_halo(const Boundry *inner_boundary, const Boundry *outer_boundary,
                                             double *coord_values[2], int *global_idx,
                                             int *num_found)
{
    if(*kernel_boundry <= *inner_boundary)
        return;
    search_points_in_halo(inner_boundary, outer_boundary, points_coord, points_global_index, num_kernel_points ,coord_values, global_idx, num_found);
}


inline bool Search_tree_node::is_coordinate_in_halo(double x, double y, const Boundry *inner, const Boundry *outer)
{
    return !(x < inner->max_lon && x >= inner->min_lon && y < inner->max_lat && y >= inner->min_lat) &&
           (x < outer->max_lon && x >= outer->min_lon && y < outer->max_lat && y >= outer->min_lat);
}


void Search_tree_node::init_num_neighbors_on_boundry(int n)
{
    num_neighbors_on_boundry[0] = num_neighbors_on_boundry[1] = num_neighbors_on_boundry[2] = num_neighbors_on_boundry[3] = n;
}


void Search_tree_node::reduce_num_neighbors_on_boundry(unsigned type)
{
    switch(type) {
        case PDLN_BOUNDRY_TYPE_U: num_neighbors_on_boundry[PDLN_UP]--; break;
        case PDLN_BOUNDRY_TYPE_D: num_neighbors_on_boundry[PDLN_DOWN]--; break;
        case PDLN_BOUNDRY_TYPE_L: num_neighbors_on_boundry[PDLN_LEFT]--; break;
        case PDLN_BOUNDRY_TYPE_R: num_neighbors_on_boundry[PDLN_RIGHT]--; break;
        case PDLN_BOUNDRY_TYPE_LR:
                                  num_neighbors_on_boundry[PDLN_LEFT]--;
                                  num_neighbors_on_boundry[PDLN_RIGHT]--;
                                  break;
        case PDLN_BOUNDRY_TYPE_NON: break;
        default: assert(false);
    }
}


void Search_tree_node::clear_expanding_count(unsigned type)
{
    switch(type) {
        case PDLN_BOUNDRY_TYPE_U: edge_expanding_count[PDLN_UP] = 0; break;
        case PDLN_BOUNDRY_TYPE_D: edge_expanding_count[PDLN_DOWN] = 0; break;
        case PDLN_BOUNDRY_TYPE_L: edge_expanding_count[PDLN_LEFT] = 0; break;
        case PDLN_BOUNDRY_TYPE_R: edge_expanding_count[PDLN_RIGHT] = 0; break;
        case PDLN_BOUNDRY_TYPE_LR:
                                  edge_expanding_count[PDLN_LEFT] = 0;
                                  edge_expanding_count[PDLN_RIGHT] = 0;
                                  break;
        case PDLN_BOUNDRY_TYPE_NON: break;
        default: assert(false);
    }
}


bool Delaunay_grid_decomposition::is_polar_node(Search_tree_node *node) const
{
    return node != NULL && node->node_type != PDLN_NODE_TYPE_COMMON;
}


#define PDLN_MAX_NUM_NEIGHBORS 128
int Delaunay_grid_decomposition::generate_trianglulation_for_local_decomp()
{
    bool* is_local_leaf_node_finished = new bool[local_leaf_nodes.size()]();
    unsigned** local_leaf_checksums   = new unsigned*[local_leaf_nodes.size()];
    unsigned** remote_leaf_checksums  = new unsigned*[local_leaf_nodes.size()];

#ifdef DEBUG
    for(unsigned i = 0; i < local_leaf_nodes.size(); i++) {
        assert(is_local_leaf_node_finished[i] == false);
    }
#endif

    for(unsigned i = 0; i < local_leaf_nodes.size(); i++) {
        local_leaf_checksums[i] = new unsigned[PDLN_MAX_NUM_NEIGHBORS];
        remote_leaf_checksums[i] = new unsigned[PDLN_MAX_NUM_NEIGHBORS];
    }

    vector<MPI_Request*> *waiting_lists = new vector<MPI_Request*> [local_leaf_nodes.size()];

    int iter = 0;
    bool all_finished = false;
    double expanding_ratio = PDLN_DEFAULT_EXPANGDING_RATIO;
    for(unsigned i = 0; i < local_leaf_nodes.size(); i++)
        local_leaf_nodes[i]->init_num_neighbors_on_boundry(1);

    while(iter < 10) {
        int ret = 0;
        for(unsigned i = 0; i < local_leaf_nodes.size(); i++)
            if(!is_local_leaf_node_finished[i]) {
                ret |= expand_tree_node_boundry(local_leaf_nodes[i], expanding_ratio);

                //if (((local_leaf_nodes[i]->kernel_boundry->min_lat+local_leaf_nodes[i]->kernel_boundry->max_lat)*0.5 < -20 && is_polar_node(search_tree_root->children[0])) ||
                //    ((local_leaf_nodes[i]->kernel_boundry->min_lat+local_leaf_nodes[i]->kernel_boundry->max_lat)*0.5 >  20 && is_polar_node(search_tree_root->children[2])))
                if (is_polar_node(search_tree_root->children[0]) || is_polar_node(search_tree_root->children[2]))
                    local_leaf_nodes[i]->project_grid();
            }

        int all_ret = 0;
        MPI_Allreduce(&ret, &all_ret, 1, MPI_UNSIGNED, MPI_LOR, processing_info->get_mpi_comm());
        if(all_ret) {
            all_finished = false;
            break;
        }

        #pragma omp parallel for
        for(unsigned i = 0; i < local_leaf_nodes.size(); i++) {
            if(!is_local_leaf_node_finished[i]) {

                timeval start, end;
                gettimeofday(&start, NULL);

                local_leaf_nodes[i]->generate_local_triangulation(is_cyclic, num_inserted);

                gettimeofday(&end, NULL);
                int rank;
                MPI_Comm_rank(processing_info->get_mpi_comm(), &rank);
#ifdef TIME_PERF
                printf("[%3d] %dth \"generate local triangulation\": %ldms, number of points: %d\n", rank, iter,
                                                                                    ((end.tv_sec - start.tv_sec) * 1000000 + (end.tv_usec - start.tv_usec)) / 1000,
                                                                                    local_leaf_nodes[i]->num_kernel_points + local_leaf_nodes[i]->num_expanded_points);
#endif
            }
        }

        #pragma omp parallel for
        for(unsigned i = 0; i < local_leaf_nodes.size(); i++) {
            if(!is_local_leaf_node_finished[i]) {
#ifdef DEBUG
                assert(local_leaf_nodes[i]->neighbors.size() < PDLN_MAX_NUM_NEIGHBORS);
#endif
                send_recv_checksums_with_neighbors(local_leaf_nodes[i], local_leaf_checksums[i], remote_leaf_checksums[i], waiting_lists + i, iter);
            }
        }

        for(unsigned i = 0; i < local_leaf_nodes.size(); i++) {
            if(!is_local_leaf_node_finished[i]) {
                for(unsigned j = 0; j < waiting_lists[i].size(); j++) {
                    MPI_Wait(waiting_lists[i][j], MPI_STATUS_IGNORE);
                    delete waiting_lists[i][j];
                }
                waiting_lists[i].clear();
            }
            else
                assert(waiting_lists[i].size() == 0);
        }

        #pragma omp parallel for
        for(unsigned i = 0; i < local_leaf_nodes.size(); i++) {
            if(!is_local_leaf_node_finished[i]) {
                is_local_leaf_node_finished[i] = are_checksums_identical(local_leaf_nodes[i], local_leaf_checksums[i], remote_leaf_checksums[i]);
                //if(iter>=1)
                //    is_local_leaf_node_finished[i] = true;

            }
        }

        all_finished = true;
        for(unsigned i = 0; i < local_leaf_nodes.size(); i++)
            if(!is_local_leaf_node_finished[i])
                all_finished = false;

        if(!all_finished)
            expanding_ratio += 0.1;
        iter++;
    }

    delete [] is_local_leaf_node_finished;

    for(unsigned i = 0; i < local_leaf_nodes.size(); i++) {
        delete [] local_leaf_checksums[i];
        delete [] remote_leaf_checksums[i];
        local_leaf_checksums[i] = NULL;
        remote_leaf_checksums[i] = NULL;
    }

    delete [] local_leaf_checksums;
    delete [] remote_leaf_checksums;

    delete [] waiting_lists;
    
    if(all_finished)
        return 0;
    else
        return 1;
}


int Delaunay_grid_decomposition::generate_trianglulation_for_whole_grid()
{
    assert(false);
    return 0;
}


void Delaunay_grid_decomposition::print_tree_node_info_recursively(Search_tree_node *node)
{
    if(node->region_ids.size() == 1){
        printf("[ID%d]x[ST-Info] LEAF %p\n", local_leaf_nodes[0]->region_id, node);
        return;
    }

    printf("[ID%d]x[ST-Info] %p: %p, %p, %p\n", local_leaf_nodes[0]->region_ids[0], node, node->children[0], node->children[1], node->children[2]);
    for(int i = 0; i < 3; i ++)
        if(node->children[i])
            print_tree_node_info_recursively(node->children[i]);
}


void Delaunay_grid_decomposition::print_whole_search_tree_info()
{
    printf("[ID%d]x[ST-Info] ROOT %p\n", local_leaf_nodes[0]->region_ids[0], search_tree_root);
    print_tree_node_info_recursively(search_tree_root);
}


#ifdef OPENCV
void Delaunay_grid_decomposition::plot_grid_decomposition(const char *filename)
{
    if (processing_info->get_local_process_id() == 0) {
        plot_points_into_file(filename, search_tree_root->points_coord[PDLN_LON], search_tree_root->points_coord[PDLN_LAT], search_tree_root->num_kernel_points, PDLN_PLOT_GLOBAL);
        for(unsigned i = 0; i < all_leaf_nodes.size(); i++) {
            Boundry b = *all_leaf_nodes[i]->kernel_boundry;
            plot_rectangle_into_file(filename, b.min_lon, b.max_lon, b.min_lat, b.max_lat, PDLN_PLOT_COLOR_RED, PDLN_PLOT_FILEMODE_APPEND);
            char number[8];
            snprintf(number, 8, "%d", all_leaf_nodes[i]->region_id);
            plot_text_into_file(filename, number, b.min_lon, b.max_lon, b.min_lat, b.max_lat, PDLN_PLOT_COLOR_RED);
        }
    }
}


void Delaunay_grid_decomposition::plot_local_triangles(const char *perfix)
{
    for(unsigned int i = 0; i < local_leaf_nodes.size(); i++) {
        char filename[64];
        snprintf(filename, 64, "%s_%d.png", perfix, local_leaf_nodes[i]->region_id);
        local_leaf_nodes[i]->triangulation->plot_into_file(filename, local_leaf_nodes[i]->kernel_boundry->min_lon,
                                                                     local_leaf_nodes[i]->kernel_boundry->max_lon,
                                                                     local_leaf_nodes[i]->kernel_boundry->min_lat,
                                                                     local_leaf_nodes[i]->kernel_boundry->max_lat);
    }
}
#endif


void delete_redundent_triangles(Triangle_Transport *&all_triangles, int &num)
{
    std::tr1::unordered_map<Triangle_Transport, std::list<int> > hash_table;
    std::tr1::unordered_map<Triangle_Transport, std::list<int> >::iterator it_hash;

    if(num == 0)
        return;

    Triangle_Transport *tmp_triangles = new Triangle_Transport[num];

    int count = 0;
    for(int i = 0; i < num; i++) {
        it_hash = hash_table.find(all_triangles[i]);
        if(it_hash != hash_table.end()) {
            bool same = false;
            for(std::list<int>::iterator it_list = it_hash->second.begin(); it_list != it_hash->second.end(); it_list ++)
                if(all_triangles[*it_list] == all_triangles[i]) {
                    same = true;
                    break;
                }
            if(same)
                continue;
            else {
                it_hash->second.push_back(i);
                tmp_triangles[count++] = all_triangles[i];
            }
        }
        else {
            hash_table[all_triangles[i]].push_back(i);
            tmp_triangles[count++] = all_triangles[i];
        }
    }

    delete[] all_triangles;
    all_triangles = tmp_triangles;
    num = count;

    return;
}


void Delaunay_grid_decomposition::save_unique_triangles_into_file(Triangle_Transport *&triangles, int num_triangles, bool sort)
{
    int num_different_triangles;
    if (sort) {
        sort_points_in_triangle(triangles, num_triangles);
        sort_triangles(triangles, num_triangles);
        int i, j;
        for(i = 0, j = 1; j < num_triangles; j++) {
            if(triangles[i].v[0].id == triangles[j].v[0].id &&
               triangles[i].v[1].id == triangles[j].v[1].id &&
               triangles[i].v[2].id == triangles[j].v[2].id) {
                continue;
            }
            else
                triangles[++i] = triangles[j];
        }
        num_different_triangles = i + 1;
    } else {
        delete_redundent_triangles(triangles, num_triangles);
        num_different_triangles = num_triangles;
    }
    
    char file_fmt[] = "log/global_triangles_%d";
    char filename[64];
    snprintf(filename, 64, file_fmt, processing_info->get_num_total_processing_units());
    FILE *fp = fopen(filename, "w");
    for(int i = 0; i < num_different_triangles; i++)
        fprintf(fp, "%d, %d, %d, (%lf, %lf), (%lf, %lf), (%lf, %lf)\n", triangles[i].v[0].id, triangles[i].v[1].id, triangles[i].v[2].id,
                                                                        triangles[i].v[0].x, triangles[i].v[0].y,
                                                                        triangles[i].v[1].x, triangles[i].v[1].y,
                                                                        triangles[i].v[2].x, triangles[i].v[2].y);
    fclose(fp);

#ifdef OPENCV
    char file_fmt2[] = "log/image_global_triangles_%d";
    snprintf(filename, 64, file_fmt2, processing_info->get_num_total_processing_units());
    plot_triangles_into_file(filename, triangles, num_different_triangles, true);
#endif
}


#define PDLN_MERGE_TAG_MASK 0x0200
void Delaunay_grid_decomposition::merge_all_triangles(bool sort)
{
    /* Let n be the number of points, if there are b vertices on the convex hull,
     * then any triangulation of the points has at most 2n − 2 − b triangles,
     * plus one exterior face */
    int local_buf_len = 0;
    for(unsigned i = 0; i < local_leaf_nodes.size(); i++)
        local_buf_len += local_leaf_nodes[i]->num_kernel_points * 3; 

    Triangle_Transport* local_triangles = new Triangle_Transport[local_buf_len];
    int num_local_triangles = 0;
    int num_triangles = 0;
    for(unsigned int i = 0; i < local_leaf_nodes.size(); i++) {
        local_leaf_nodes[i]->triangulation->get_triangles_in_region(local_leaf_nodes[i]->kernel_boundry->min_lon, local_leaf_nodes[i]->kernel_boundry->max_lon,
                                                                    local_leaf_nodes[i]->kernel_boundry->min_lat, local_leaf_nodes[i]->kernel_boundry->max_lat,
                                                                    local_triangles + num_local_triangles, &num_triangles, local_buf_len - num_local_triangles); // This function call is un-class-elegant
        num_local_triangles += num_triangles;
    }


    MPI_Barrier(processing_info->get_mpi_comm());

    if(processing_info->get_local_process_id() == 0) {
        int *num_remote_triangles = new int[processing_info->get_num_total_processes()];
        int remote_buf_len = 0;
        MPI_Status status;

        for(int i = 1; i < processing_info->get_num_total_processes(); i++)
            MPI_Recv(&num_remote_triangles[i], 1, MPI_INT, i, PDLN_MERGE_TAG_MASK, processing_info->get_mpi_comm(), &status);
        for(int i = 1; i < processing_info->get_num_total_processes(); i++) {
            //assert(num_remote_triangles[i] > min_points_per_chunk/2);
            remote_buf_len += num_remote_triangles[i];
        }
        Triangle_Transport *remote_triangles = new Triangle_Transport[remote_buf_len + num_local_triangles];

        int count = 0;
        for(int i = 1; i < processing_info->get_num_total_processes(); i++) {
            MPI_Recv(remote_triangles + count, num_remote_triangles[i] * sizeof(Triangle_Transport), MPI_CHAR, i, PDLN_MERGE_TAG_MASK, processing_info->get_mpi_comm(), &status);
            int tmp_count;
            MPI_Get_count(&status, MPI_CHAR, &tmp_count);
            
            /*
            char filename[64];
            snprintf(filename, 64, "log/process_local_triangles%d", i);
            plot_triangles_into_file(filename, remote_triangles+count, tmp_count/sizeof(Triangle_Transport));

            for(int j = 0; j < tmp_count/sizeof(Triangle_Transport); j++)
                printf("%d, %d, %d\n", (remote_triangles + count)[j].v[0].id, (remote_triangles + count)[j].v[1].id, (remote_triangles + count)[j].v[2].id);
            printf("==============\n");
            */
#ifdef DEBUG
            assert(tmp_count % sizeof(Triangle_Transport) == 0);
#endif
            count += tmp_count / sizeof(Triangle_Transport);
        }
        assert(count == remote_buf_len);
        memcpy(remote_triangles + remote_buf_len, local_triangles, num_local_triangles * sizeof(Triangle_Transport));
        save_unique_triangles_into_file(remote_triangles, remote_buf_len + num_local_triangles, sort);
        delete[] remote_triangles;
        delete[] num_remote_triangles;
    }
    else {
        MPI_Send(&num_local_triangles, 1, MPI_INT, 0, PDLN_MERGE_TAG_MASK, processing_info->get_mpi_comm());
        MPI_Send(local_triangles, num_local_triangles * sizeof(Triangle_Transport), MPI_CHAR, 0, PDLN_MERGE_TAG_MASK, processing_info->get_mpi_comm());
    }
}


static double fRand(double fMin, double fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}


void Grid_info_manager::gen_basic_grid()
{
    int size = 300;

    num_points = size * size;
    coord_values[0] = new double[num_points]();
    coord_values[1] = new double[num_points]();
    srand(0);
    for(int i = 0; i < size; i++)
        for(int j = 0; j < size; j++) {
            coord_values[0][i * size + j] =  0.0  + 360.0 * j / size;
            //coord_values[1][i * size + j] = -89.0 + 178.0 * i / size;
            coord_values[1][i * size + j] = -89.0 + 178.0 * i / size;
            //coord_values[0][i * size + j] = fRand(0.0, 359.0);
            //coord_values[1][i * size + j] = fRand(-89.0, 89.0);
        } 

    coord_values[0][0] = 0.0;
    coord_values[1][0] = -90.0;
    coord_values[0][299] = 0.0;
    coord_values[1][299] = 90.0;

    //min_lat = -89.0;
    //max_lat =  89.0;
    //max_lon = 359.0;
}


#ifdef NETCDF
void Grid_info_manager::gen_three_polar_grid()
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
    //printf("num points: %d\n", num_points);
    assert(have_redundent_points(coord_values[PDLN_LON], coord_values[PDLN_LAT], num_points) == false);

    if(squeeze) {
        for(int i = 0; i < num_points/100; i++) {
            coord_values[PDLN_LON][i] = coord_values[PDLN_LON][i*100];
            coord_values[PDLN_LAT][i] = coord_values[PDLN_LAT][i*100];
        }
        num_points /= 100;
        //printf("num points: %d\n", num_points);
    }

    min_lon =   0.0;
    max_lon = 360.0;
    min_lat = -80.0;
    max_lat =  90.0;
}

void Grid_info_manager::gen_latlon_grid()
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
    min_lat = -89.0;
    max_lat =  89.0;
}


void Grid_info_manager::gen_latlon_90_grid()
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
#endif


Grid_info_manager::Grid_info_manager()
{
    coord_values[0] = coord_values[1] = NULL;
    //gen_three_polar_grid();
    //gen_latlon_grid();
    //gen_latlon_90_grid();
}


Grid_info_manager::Grid_info_manager(double *coord[2], int num)
{

}


Grid_info_manager::~Grid_info_manager()
{
    delete coord_values[0];
    delete coord_values[1];
}
double** Grid_info_manager::get_grid_coord_values(int grid_id)
{
    return coord_values;
}
int Grid_info_manager::get_grid_num_points(int grid_id)
{
    return num_points;
}
void Grid_info_manager::get_grid_boundry(int grid_id, double* mi_lon, double* ma_lon, double* mi_lat, double* ma_lat)
{
    *mi_lon = min_lon;
    *ma_lon = max_lon;
    *mi_lat = min_lat;
    *ma_lat = max_lat;
}


void Grid_info_manager::set_grid_boundry(int grid_id, double mi_lon, double ma_lon, double mi_lat, double ma_lat)
{
    min_lon = mi_lon;
    max_lon = ma_lon;
    min_lat = mi_lat;
    max_lat = ma_lat;
}


bool Grid_info_manager::is_grid_cyclic(int grid_id)
{
    return true;
}
