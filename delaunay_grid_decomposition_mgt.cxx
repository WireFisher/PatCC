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

#define DEFAULT_EXPANGDING_RATIO (0.2)
#define PDLN_EXPECTED_EXPANDING_LOOP_TIMES (3)

#define PDLN_SPOLAR_MAX_LAT (0.0)
#define PDLN_NPOLAR_MIN_LAT (0.0)

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

Boundry& Boundry::operator= (Boundry &boundry)
{
    min_lat = boundry.min_lat;
    min_lon = boundry.min_lon;
    max_lat = boundry.max_lat;
    max_lon = boundry.max_lon;
    return *this;
}


bool Boundry::operator== (Boundry &boundry)
{
    return min_lat == boundry.min_lat && min_lon == boundry.min_lon && max_lat == boundry.max_lat && max_lon == boundry.max_lon;
}


bool Boundry::operator!= (Boundry &boundry)
{
    return !(min_lat == boundry.min_lat && min_lon == boundry.min_lon && max_lat == boundry.max_lat && max_lon == boundry.max_lon);
}


bool Boundry::operator<= (Boundry &boundry)
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


void Boundry::legalize(const Boundry *outer_boundry, bool is_cyclic)
{
    min_lat = std::max(min_lat, outer_boundry->min_lat);
    max_lat = std::min(max_lat, outer_boundry->max_lat);
    if(!is_cyclic) {
        min_lon = std::max(min_lon, outer_boundry->min_lon);
        max_lon = std::min(max_lon, outer_boundry->max_lon);
    }
}


void Boundry::max(double min_lo, double max_lo, double min_la, double max_la)
{
    min_lon = std::min(min_lon, min_lo);
    max_lon = std::max(max_lon, max_lo);
    min_lat = std::min(min_lat, min_la);
    max_lat = std::max(max_lat, max_la);
}

Search_tree_node::Search_tree_node(Search_tree_node *parent, double *coord_value[2], int *global_index, int num_points, Boundry boundry, int type) {
    this->parent = parent;
    this->children[0] = this->children[1] = this->children[2] = NULL;

    this->kernel_boundry = new Boundry();
    this->expanded_boundry = new Boundry();
    this->real_boundry = NULL;
    *this->kernel_boundry = *this->expanded_boundry = boundry;
    this->rotated_kernel_boundry = this->rotated_expanded_boundry = NULL;

    this->midline.type = -1;
    this->midline.value = -361.0;
    this->points_coord[0] = new double[num_points];
    this->points_coord[1] = new double[num_points];
    this->points_global_index = new int[num_points];
    memcpy(this->points_coord[0], coord_value[0], num_points * sizeof(double));
    memcpy(this->points_coord[1], coord_value[1], num_points * sizeof(double));
    memcpy(this->points_global_index, global_index, num_points * sizeof(int));
    this->len_expanded_points_coord_buf = 0;
    this->projected_coord[0] = this->projected_coord[1] = NULL;
    this->num_rotated_points = 0;

    this->num_kernel_points = num_points;
    this->num_expanded_points = 0;
    this->node_type = type;
    this->triangulation = NULL;

    if(type == PDLN_NODE_TYPE_COMMON) {
        this->center[PDLN_LON] = (boundry.min_lon + boundry.max_lon) * 0.5;
        this->center[PDLN_LAT] = (boundry.min_lat + boundry.max_lat) * 0.5;
    }
    else if(type == PDLN_NODE_TYPE_SPOLAR) {
        this->center[PDLN_LON] = 0.0;
        this->center[PDLN_LAT] = -90.0;
    }
    else if(type == PDLN_NODE_TYPE_NPOLAR) {
        this->center[PDLN_LON] = 0.0;
        this->center[PDLN_LAT] = 90.0;
    }

    this->non_monotonic = kernel_boundry->min_lon > kernel_boundry->max_lon;
    this->virtual_point_local_index = -1;
}


Search_tree_node::~Search_tree_node()
{
    delete kernel_boundry;
    delete expanded_boundry;
    delete[] points_coord[0];
    delete[] points_coord[1];
    delete[] points_global_index;
    if(rotated_kernel_boundry != NULL)
        delete rotated_kernel_boundry;
    if(rotated_expanded_boundry != NULL)
        delete rotated_expanded_boundry;
    for(int i = 0; i < 3; i ++)
        if(children[i] != NULL)
            delete children[i];

    if(triangulation != NULL)
        delete triangulation;
}


/* assumption: ids are already sorted. */
void Search_tree_node::update_processing_units_id(int num) 
{
    processing_units_id.clear();
    for(int i = 0; i < num; i ++)
        processing_units_id.push_back(i);
}


void Search_tree_node::update_processing_units_id(vector<int> proc_units_id) 
{
    processing_units_id.clear();
    processing_units_id = proc_units_id;
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


#define PDLN_INSERT_VIRTUAL_POINT true
void Search_tree_node::generate_local_triangulation(bool is_cyclic)
{
    if(triangulation != NULL)
        delete triangulation;

    //plot_points_info_file(filename, points_coord[PDLN_LON], points_coord[PDLN_LAT], num_kernel_points + num_expanded_points);

    if(rotated_expanded_boundry != NULL) {
        triangulation = new Delaunay_Voronoi(num_kernel_points + num_expanded_points,
                                             projected_coord[PDLN_LON], projected_coord[PDLN_LAT], points_global_index, false,
                                             rotated_expanded_boundry->min_lon, rotated_expanded_boundry->max_lon,
                                             rotated_expanded_boundry->min_lat, rotated_expanded_boundry->max_lat, NULL, virtual_point_local_index);

        if(node_type != PDLN_NODE_TYPE_COMMON) {
            //triangulation->plot_original_points_into_file(filename);

            double lon, head_lon, head_lat, tail_lon, tail_lat;

            calculate_real_boundary(); // TODO: delete it
            if (node_type == PDLN_NODE_TYPE_NPOLAR) {
                lon = 360 - 1e-8;
                calculate_stereographic_projection(lon, real_boundry->max_lat-0.1, center[PDLN_LON], center[PDLN_LAT], head_lon, head_lat);
                calculate_stereographic_projection(lon, real_boundry->min_lat,     center[PDLN_LON], center[PDLN_LAT], tail_lon, tail_lat);
            }
            else {
                lon = 360 - 1e-8; // TODO: is this OK?
                calculate_stereographic_projection(lon, real_boundry->min_lat+0.1, center[PDLN_LON], center[PDLN_LAT], head_lon, head_lat);
                calculate_stereographic_projection(lon, real_boundry->max_lat,     center[PDLN_LON], center[PDLN_LAT], tail_lon, tail_lat);
            }

            std::vector<Triangle*> cyclic_triangles = triangulation->search_cyclic_triangles_for_rotated_grid(Point(head_lon, head_lat), Point(tail_lon, tail_lat));


            char filename[64];
            int rank;
            MPI_Comm_rank(process_thread_mgr->get_mpi_comm(), &rank);
            snprintf(filename, 64, "log/projected_triangles_%d.png", rank);
            triangulation->plot_projection_into_file(filename);

            //snprintf(filename, 64, "log/fake_cyclic_triangles_%d.png", rank);
            //plot_triangles_info_file(filename, cyclic_triangles);

            triangulation->update_all_points_coord(points_coord[PDLN_LON], points_coord[PDLN_LAT], num_kernel_points + num_expanded_points);
            triangulation->correct_cyclic_triangles(cyclic_triangles, is_cyclic);
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

            if(real_boundry->min_lat < 0) {
                calculate_stereographic_projection(0, real_boundry->min_lat, this->center[PDLN_LON], this->center[PDLN_LAT], x[0], y[0]);
                calculate_stereographic_projection(90, real_boundry->min_lat, this->center[PDLN_LON], this->center[PDLN_LAT], x[1], y[1]);
                calculate_stereographic_projection(180, real_boundry->min_lat, this->center[PDLN_LON], this->center[PDLN_LAT], x[2], y[2]);

                calculate_circle_center(x, y, &circle_center.x, &circle_center.y);

                radius = sqrt((x[2]-circle_center.x)*(x[2]-circle_center.x)+(y[2]-circle_center.y)*(y[2]-circle_center.y));

                if(radius < 100) {
                    //printf("[%d] + circle_center: (%lf, %lf), point: (%lf, %lf) radius: %lf\n", rank, circle_center.x, circle_center.y, x[2], y[2], radius);
                    triangulation->remove_triangles_in_circle(circle_center, radius);

                    double lon, head_lon, head_lat, tail_lon, tail_lat;
                    lon = (real_boundry->max_lon + real_boundry->min_lon + 360.0) * 0.5;
                    calculate_stereographic_projection(lon, -center[PDLN_LAT]-20 , center[PDLN_LON], center[PDLN_LAT], head_lon, head_lat);
                    calculate_stereographic_projection(lon, real_boundry->min_lat, center[PDLN_LON], center[PDLN_LAT], tail_lon, tail_lat);

                    //printf("[%d]real: (%lf, %lf, %lf, %lf), (%lf, %lf) -- (%lf, %lf)\n", rank, real_boundry->min_lon, real_boundry->max_lon,
                    //                                                               real_boundry->min_lat, real_boundry->max_lat,
                    //                                                             head_lon, head_lat, tail_lon, tail_lat);
                    triangulation->remove_triangles_on_segment(Point(head_lon, head_lat), Point(tail_lon, tail_lat));
                }
            }

            if(real_boundry->max_lat > 0) {
                calculate_stereographic_projection(0, real_boundry->max_lat, this->center[PDLN_LON], this->center[PDLN_LAT], x[0], y[0]);
                calculate_stereographic_projection(90, real_boundry->max_lat, this->center[PDLN_LON], this->center[PDLN_LAT], x[1], y[1]);
                calculate_stereographic_projection(180, real_boundry->max_lat, this->center[PDLN_LON], this->center[PDLN_LAT], x[2], y[2]);

                calculate_circle_center(x, y, &circle_center.x, &circle_center.y);

                radius = sqrt((x[2]-circle_center.x)*(x[2]-circle_center.x)+(y[2]-circle_center.y)*(y[2]-circle_center.y));

                if(radius < 100) {
                    triangulation->remove_triangles_in_circle(circle_center, radius);

                    double lon, head_lon, head_lat, tail_lon, tail_lat;
                    lon = (real_boundry->max_lon + real_boundry->min_lon + 360.0) * 0.5;
                    calculate_stereographic_projection(lon, real_boundry->max_lat+0.1, center[PDLN_LON], center[PDLN_LAT], head_lon, head_lat);
                    calculate_stereographic_projection(lon, -center[PDLN_LAT]+20 , center[PDLN_LON], center[PDLN_LAT], tail_lon, tail_lat);

                    triangulation->remove_triangles_on_segment(Point(head_lon, head_lat), Point(tail_lon, tail_lat));
                }
            }

            triangulation->update_all_points_coord(points_coord[PDLN_LON], points_coord[PDLN_LAT], num_kernel_points + num_expanded_points);
            triangulation->remove_triangles_on_or_out_of_boundary(real_boundry->min_lon, real_boundry->max_lon, real_boundry->min_lat, real_boundry->max_lat);
            triangulation->relegalize_all_triangles();
        }
    }
    else {
        triangulation = new Delaunay_Voronoi(num_kernel_points + num_expanded_points,
                                             points_coord[PDLN_LON], points_coord[PDLN_LAT], points_global_index, false,
                                             expanded_boundry->min_lon, expanded_boundry->max_lon,
                                             expanded_boundry->min_lat, expanded_boundry->max_lat, NULL, virtual_point_local_index);

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


void Search_tree_node::split_local_points(Midline midline, double *child_points_coord[4], int *child_points_idx[2], int child_num_points[2])
{
    child_num_points[0] = 0;
    child_num_points[1] = 0;

    if(non_monotonic && midline.type == PDLN_LON)
        assert(false);
    else {
        #pragma omp parallel for
        for(int i = 0; i < num_kernel_points; i ++) {
            if(points_coord[midline.type][i] < midline.value) {
                child_points_coord[midline.type][child_num_points[0]] = points_coord[midline.type][i];
                child_points_coord[(midline.type+1)%2][child_num_points[0]] = points_coord[(midline.type+1)%2][i];
                #pragma omp critical
                child_points_idx[0][child_num_points[0]++] = points_global_index[i];
            }
            else {
                child_points_coord[2+midline.type][child_num_points[1]] = points_coord[midline.type][i];
                child_points_coord[2+(midline.type+1)%2][child_num_points[1]] = points_coord[(midline.type+1)%2][i];
                #pragma omp critical 
                child_points_idx[1][child_num_points[1]++] = points_global_index[i];
            }
        }
    }
}


void Search_tree_node::decompose_by_processing_units_number(double *workloads, double *child_points_coord[4], int *child_points_idx[2], int child_num_points[2],
                                   Boundry child_boundry[2], vector<int> child_proc_units_id[2], int mode)
{
    double length[2], boundry_values[4], child_total_workload[2];
    Midline midline;
    unsigned int i;
    int iteration_count;

    assert(processing_units_id.size() > 1);
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
        for(i = 0; i < processing_units_id.size()/2; i++)
            child_proc_units_id[0].push_back(processing_units_id[i]);
        for(; i < processing_units_id.size(); i++)
            child_proc_units_id[1].push_back(processing_units_id[i]);
    }
    else if(mode == PDLN_DECOMPOSE_SPOLAR_MODE) {
        child_proc_units_id[0].push_back(processing_units_id[0]);
        for(i = 1; i < processing_units_id.size(); i++)
            child_proc_units_id[1].push_back(processing_units_id[i]);
    }
    else if(mode == PDLN_DECOMPOSE_NPOLAR_MODE) {
        for(i = 0; i < processing_units_id.size()-1; i++)
            child_proc_units_id[0].push_back(processing_units_id[i]);
        child_proc_units_id[1].push_back(processing_units_id[i]);
    }
    else
        assert(false);

#ifdef DEBUG
    assert(child_proc_units_id[0].size() + child_proc_units_id[1].size() == processing_units_id.size());
#endif

    if(processing_units_id.size() > 1) {
        for(i = 0, child_total_workload[0] = 0.0; i < child_proc_units_id[0].size(); i++)
            child_total_workload[0] += workloads[child_proc_units_id[0][i]];
        for(i = 0, child_total_workload[1] = 0.0; i < child_proc_units_id[1].size(); i++)
            child_total_workload[1] += workloads[child_proc_units_id[1][i]];

        midline.value = boundry_values[midline.type] + length[midline.type] * child_total_workload[0] / (child_total_workload[0] + child_total_workload[1]);
        split_local_points(midline, child_points_coord, child_points_idx, child_num_points);

        iteration_count = 1;
        while (child_num_points[0] == 0 || child_num_points[1] == 0 ||
               fabs(child_num_points[0]/child_num_points[1] - child_total_workload[0]/child_total_workload[1]) > PDLN_TOLERABLE_ERROR) {

            if(iteration_count++ > PDLN_MAX_ITER_COUNT)
                break;

            if(child_num_points[0] < child_total_workload[0]) {
#ifdef DEBUG
                assert(child_num_points[1] >= child_total_workload[1] || fabs(child_total_workload[1] - child_num_points[1]) < PDLN_FLOAT_EQ_ERROR);
#endif
                midline.value += (boundry_values[2+midline.type] - midline.value) * (child_num_points[1] - child_total_workload[1]) / child_num_points[1];
            }
            else {
#ifdef DEBUG
                assert(child_num_points[1] <= child_total_workload[1] || fabs(child_total_workload[1] - child_num_points[1]) < PDLN_FLOAT_EQ_ERROR);
#endif
                midline.value -= (midline.value - boundry_values[midline.type]) * (child_num_points[0] - child_total_workload[0]) / child_num_points[0];
            }
            assert(midline.value > boundry_values[midline.type]);
            assert(midline.value < boundry_values[2+midline.type]);
            /* TODO: Search only half of the whole points, but not the whole points */
            split_local_points(midline, child_points_coord, child_points_idx, child_num_points);
        }
    }
    else
        midline.value = boundry_values[2+midline.type];

    if(midline.type == PDLN_LON) {
        child_boundry[0].min_lat = child_boundry[1].min_lat = kernel_boundry->min_lat;
        child_boundry[0].max_lat = child_boundry[1].max_lat = kernel_boundry->max_lat;
        child_boundry[0].min_lon = kernel_boundry->min_lon;
        child_boundry[0].max_lon = child_boundry[1].min_lon = midline.value;
        child_boundry[1].max_lon = kernel_boundry->max_lon;
    }
    else if(midline.type == PDLN_LAT) {
        child_boundry[0].min_lon = child_boundry[1].min_lon = kernel_boundry->min_lon;
        child_boundry[0].max_lon = child_boundry[1].max_lon = kernel_boundry->max_lon;
        child_boundry[0].min_lat = kernel_boundry->min_lat;
        child_boundry[0].max_lat = child_boundry[1].min_lat = midline.value;
        child_boundry[1].max_lat = kernel_boundry->max_lat;
    }
    else
        assert(false);

}


void Search_tree_node::decompose_by_fixed_longitude(double fixed_lon, double *workloads, double *child_points_coord[4], int *child_points_idx[2], int child_num_points[2],
                                                    Boundry child_boundry[2], vector<int> child_proc_units_id[2])
{
    Midline midline;

    assert(processing_units_id.size() > 0);
    assert(fixed_lon >= 0 && fixed_lon < 360);
    
    midline.type = PDLN_LON;
    midline.value = fixed_lon;

    split_local_points(midline, child_points_coord, child_points_idx, child_num_points);

    split_processing_units_by_points_number(workloads, child_num_points[0], child_num_points[1], processing_units_id, child_proc_units_id);

    assert(child_proc_units_id[0].size() + child_proc_units_id[1].size() == processing_units_id.size());

    child_boundry[0] = child_boundry[1] = *kernel_boundry;
    child_boundry[0].max_lon = child_boundry[1].min_lon = fixed_lon;
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
        if(processing_units_id[0] == ns[i]->processing_units_id[0])
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

        num_rotated_points = num_kernel_points + num_expanded_points;
    }
    else {
        #pragma omp parallel for
        for(int i = num_rotated_points; i < num_kernel_points + num_expanded_points; i++) {
            calculate_stereographic_projection(points_coord[PDLN_LON][i], points_coord[PDLN_LAT][i], center[PDLN_LON], center[PDLN_LAT],
                                               projected_coord[PDLN_LON][i], projected_coord[PDLN_LAT][i]);
        }
        num_rotated_points = num_kernel_points + num_expanded_points;
    }

    /* recalculate expanded boundary */
    double top = -1e20, bot = 1e20, left = 1e20, right = -1e20;
    for(int i = 0; i < num_rotated_points; i++) { //TODO: i can be started from non-zero
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


Delaunay_grid_decomposition::Delaunay_grid_decomposition(int grid_id, Processing_resource *proc_info, int min_num_points_per_chunk)
{
    double **coord_values;
    int *global_index;
    int num_points;
    Boundry boundry;

    this->original_grid = grid_id;
    this->processing_info = proc_info;
    assert(this->processing_info != NULL);

    this->min_num_points_per_chunk = min_num_points_per_chunk;
    
    coord_values = grid_info_mgr->get_grid_coord_values(grid_id);
    num_points = grid_info_mgr->get_grid_num_points(grid_id);


    this->is_cyclic = grid_info_mgr->is_grid_cyclic(grid_id);
    grid_info_mgr->get_grid_boundry(grid_id, &boundry.min_lon, &boundry.max_lon, &boundry.min_lat, &boundry.max_lat);
    this->processing_info->get_num_total_processing_units();

    this->active_processing_units_flag = new bool[processing_info->get_num_total_processing_units()];

    global_index = new int[num_points];
    for(int i = 0; i < num_points; i++)
        global_index[i] = i;

    if(boundry.max_lon - boundry.min_lon < 360.0)
        boundry.max_lon += PDLN_HIGH_BOUNDRY_SHIFTING;
    boundry.max_lat += PDLN_HIGH_BOUNDRY_SHIFTING;

    assert(boundry.max_lon - boundry.min_lon <= 360.0);
    search_tree_root = new Search_tree_node(NULL, coord_values, global_index, num_points, boundry, PDLN_NODE_TYPE_COMMON);

    //if(processing_info->get_local_process_id() == 0) {
    //    char filename[64];
    //    snprintf(filename, 64, "log/original_input_points.png");
    //    plot_points_info_file(filename, coord_values[PDLN_LON], coord_values[PDLN_LAT], num_points, PDLN_PLOT_GLOBAL);
    //}

    this->workloads = NULL;
    //this->initialze_workload();
    
    delete global_index;
}


Delaunay_grid_decomposition::~Delaunay_grid_decomposition()
{
    delete search_tree_root;
    if(workloads != NULL)
        delete[] workloads; 
    delete[] active_processing_units_flag;
}


void Delaunay_grid_decomposition::initialze_workload()
{
    int max_num_processing_units, num_active_processing_units;
    double average_workload;
    std::vector<int> active_processing_common_id;

    assert(min_num_points_per_chunk > 0);
    max_num_processing_units = (grid_info_mgr->get_grid_num_points(original_grid) + min_num_points_per_chunk - 1) / min_num_points_per_chunk;

    num_active_processing_units = std::min(processing_info->get_num_total_processing_units(), max_num_processing_units);

    average_workload = (double)grid_info_mgr->get_grid_num_points(original_grid) / num_active_processing_units;

    if(workloads != NULL)
        delete[] workloads;
    workloads = new double[processing_info->get_num_total_processing_units()];

    processing_info->pick_out_active_processing_units(num_active_processing_units, active_processing_units_flag);

    for(int i = 0; i < processing_info->get_num_total_processing_units(); i++) {
        if(active_processing_units_flag[i]) {
            workloads[i] = average_workload;
            active_processing_common_id.push_back(i);
        }
        else
            workloads[i] = 0;
    }

    assert(active_processing_common_id.size() == (unsigned int)num_active_processing_units);

    search_tree_root->update_processing_units_id(active_processing_common_id);
}


void Delaunay_grid_decomposition::update_workloads(int total_workload, vector<int> &ids)
{
    double old_total_workload = 0.0;
    double unassigned_workload;

    /* NOTE: This will make the final workload of leaf node not exactly more than min_num_points_per_chunk */
    if(ids.size() == 1) {
        workloads[ids[0]] = total_workload;
        return;
    }

    for(unsigned int i = 0; i < ids.size(); i++)
        old_total_workload += workloads[ids[i]];

    for(unsigned int i = 0; i < ids.size(); i++)
        workloads[ids[i]] = workloads[ids[i]] * total_workload / old_total_workload;

    for(unsigned int i = 0; i < ids.size();) {
        if(ids.size() < 2)
            break;

        if(workloads[ids[i]] < min_num_points_per_chunk) {
            unassigned_workload = workloads[ids[i]];
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


/* "common_node" means non-polar node */
void Delaunay_grid_decomposition::decompose_common_node_recursively(Search_tree_node *node)
{
    double *child_points_coord[4];
    int *child_points_index[2];
    int child_num_points[2];
    Boundry child_boundry[2];
    vector<int> child_proc_id[2];

    assert(node->processing_units_id.size() > 0);
    if(node->processing_units_id.size() == 1) {
        if(have_local_processing_units_id(node->processing_units_id))
            local_leaf_nodes.push_back(node);
        return;
    }
    
    for(int i = 0; i < 4; i++)
        child_points_coord[i] = new double[node->num_kernel_points];
    child_points_index[0] = new int[node->num_kernel_points];
    child_points_index[1] = new int[node->num_kernel_points];
    
    node->decompose_by_processing_units_number(workloads, child_points_coord, child_points_index, child_num_points, child_boundry, child_proc_id, PDLN_DECOMPOSE_COMMON_MODE);

    node->children[0] = alloc_search_tree_node(node, child_points_coord,   child_points_index[0], child_num_points[0], child_boundry[0], child_proc_id[0], PDLN_NODE_TYPE_COMMON);
    node->children[2] = alloc_search_tree_node(node, child_points_coord+2, child_points_index[1], child_num_points[1], child_boundry[1], child_proc_id[1], PDLN_NODE_TYPE_COMMON);

    for(int i = 0; i < 4; i++)
        delete[] child_points_coord[i];
    delete[] child_points_index[0];
    delete[] child_points_index[1];
    //TODO: optimize new delete
    
    //int rank;
    //MPI_Comm_rank(processing_info->get_mpi_comm(), &rank);
    //printf("[Rank%d]x[ST-INFO-PRE] p: %p, first: %p, third: %p\n", rank, node, node->children[0], node->children[2]);

    if(have_local_processing_units_id(node->children[0]->processing_units_id)) 
        decompose_common_node_recursively(node->children[0]);
    if(have_local_processing_units_id(node->children[2]->processing_units_id)) 
        decompose_common_node_recursively(node->children[2]);
}


Search_tree_node* Delaunay_grid_decomposition::alloc_search_tree_node(Search_tree_node* parent, double *coord_values[2], int *index, 
                                                                      int num_points, Boundry boundary, vector<int> &common_ids, int type)
{
    assert(common_ids.size() > 0);
    Search_tree_node *new_node = new Search_tree_node(parent, coord_values, index, num_points, boundary, type);

    /* common_ids can be modified by update_workloads */
    update_workloads(num_points, common_ids);

    new_node->update_processing_units_id(common_ids);
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


void Delaunay_grid_decomposition::compute_common_boundry(Search_tree_node *tree_node1, Search_tree_node *tree_node2,
                                                         Point *boundry_head, Point *boundry_tail, Point *boundry2_head, Point *boundry2_tail)
{
    double coord_value[2][2];

    coord_value[0][PDLN_LAT] = coord_value[0][PDLN_LON] = coord_value[1][PDLN_LAT] = coord_value[1][PDLN_LON] = PDLN_DOUBLE_INVALID_VALUE;
    if(fabs(tree_node1->kernel_boundry->min_lat - tree_node2->kernel_boundry->max_lat) < PDLN_FLOAT_EQ_ERROR) {
        if(std::max(tree_node1->kernel_boundry->min_lon, tree_node2->kernel_boundry->min_lon) < 
           std::min(tree_node1->kernel_boundry->max_lon, tree_node2->kernel_boundry->max_lon)) {
            coord_value[0][PDLN_LAT] = coord_value[1][PDLN_LAT] = tree_node2->kernel_boundry->max_lat;
            coord_value[0][PDLN_LON] = std::max(tree_node1->kernel_boundry->min_lon, tree_node2->kernel_boundry->min_lon);
            coord_value[1][PDLN_LON] = std::min(tree_node1->kernel_boundry->max_lon, tree_node2->kernel_boundry->max_lon);
        }
    }
    else if(fabs(tree_node1->kernel_boundry->max_lat - tree_node2->kernel_boundry->min_lat) < PDLN_FLOAT_EQ_ERROR) {
        if(std::max(tree_node1->kernel_boundry->min_lon, tree_node2->kernel_boundry->min_lon) <
           std::min(tree_node1->kernel_boundry->max_lon, tree_node2->kernel_boundry->max_lon)) {
            coord_value[0][PDLN_LAT] = coord_value[1][PDLN_LAT] = tree_node1->kernel_boundry->max_lat;
            coord_value[0][PDLN_LON] = std::max(tree_node1->kernel_boundry->min_lon, tree_node2->kernel_boundry->min_lon);
            coord_value[1][PDLN_LON] = std::min(tree_node1->kernel_boundry->max_lon, tree_node2->kernel_boundry->max_lon);
        }
    }
    else if(fabs(tree_node1->kernel_boundry->min_lon - tree_node2->kernel_boundry->max_lon) < PDLN_FLOAT_EQ_ERROR) {
        if(std::max(tree_node1->kernel_boundry->min_lat, tree_node2->kernel_boundry->min_lat) <
           std::min(tree_node1->kernel_boundry->max_lat, tree_node2->kernel_boundry->max_lat)) {
            coord_value[0][PDLN_LON] = coord_value[1][PDLN_LON] = tree_node2->kernel_boundry->max_lon;
            coord_value[0][PDLN_LAT] = std::max(tree_node1->kernel_boundry->min_lat, tree_node2->kernel_boundry->min_lat);
            coord_value[1][PDLN_LAT] = std::min(tree_node1->kernel_boundry->max_lat, tree_node2->kernel_boundry->max_lat);
        }
    }
    else if(fabs(tree_node1->kernel_boundry->max_lon - tree_node2->kernel_boundry->min_lon) < PDLN_FLOAT_EQ_ERROR) {
        if(std::max(tree_node1->kernel_boundry->min_lat, tree_node2->kernel_boundry->min_lat) <
           std::min(tree_node1->kernel_boundry->max_lat, tree_node2->kernel_boundry->max_lat)) {
            coord_value[0][PDLN_LON] = coord_value[1][PDLN_LON] = tree_node1->kernel_boundry->max_lon;
            coord_value[0][PDLN_LAT] = std::max(tree_node1->kernel_boundry->min_lat, tree_node2->kernel_boundry->min_lat);
            coord_value[1][PDLN_LAT] = std::min(tree_node1->kernel_boundry->max_lat, tree_node2->kernel_boundry->max_lat);
        }
    }
    *boundry_head = Point(coord_value[0][PDLN_LON], coord_value[0][PDLN_LAT]);
    *boundry_tail = Point(coord_value[1][PDLN_LON], coord_value[1][PDLN_LAT]);

    coord_value[0][PDLN_LAT] = coord_value[0][PDLN_LON] = coord_value[1][PDLN_LAT] = coord_value[1][PDLN_LON] = PDLN_DOUBLE_INVALID_VALUE;
    if(fabs(fabs(tree_node1->kernel_boundry->min_lon - tree_node2->kernel_boundry->max_lon) - 360.0) < PDLN_FLOAT_EQ_ERROR) {
        if(std::max(tree_node1->kernel_boundry->min_lat, tree_node2->kernel_boundry->min_lat) < 
           std::min(tree_node1->kernel_boundry->max_lat, tree_node2->kernel_boundry->max_lat)) {
            coord_value[0][PDLN_LON] = coord_value[1][PDLN_LON] = tree_node1->kernel_boundry->min_lon;
            coord_value[0][PDLN_LAT] = std::max(tree_node1->kernel_boundry->min_lat, tree_node2->kernel_boundry->min_lat);
            coord_value[1][PDLN_LAT] = std::min(tree_node1->kernel_boundry->max_lat, tree_node2->kernel_boundry->max_lat);
        }
    }
    else if(fabs(fabs(tree_node1->kernel_boundry->max_lon - tree_node2->kernel_boundry->min_lon) - 360.0) < PDLN_FLOAT_EQ_ERROR) {
        if(std::max(tree_node1->kernel_boundry->min_lat, tree_node2->kernel_boundry->min_lat) <
           std::min(tree_node1->kernel_boundry->max_lat, tree_node2->kernel_boundry->max_lat)) {
            coord_value[0][PDLN_LON] = coord_value[1][PDLN_LON] = tree_node1->kernel_boundry->max_lon;
            coord_value[0][PDLN_LAT] = std::max(tree_node1->kernel_boundry->min_lat, tree_node2->kernel_boundry->min_lat);
            coord_value[1][PDLN_LAT] = std::min(tree_node1->kernel_boundry->max_lat, tree_node2->kernel_boundry->max_lat);
        }
    }
    *boundry2_head = Point(coord_value[0][PDLN_LON], coord_value[0][PDLN_LAT]);
    *boundry2_tail = Point(coord_value[1][PDLN_LON], coord_value[1][PDLN_LAT]);
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
    unsigned checksum;

    double threshold = std::min(search_tree_root->kernel_boundry->max_lon - search_tree_root->kernel_boundry->min_lon,
                                search_tree_root->kernel_boundry->max_lat - search_tree_root->kernel_boundry->min_lat) /
                       sqrt(processing_info->get_num_total_processing_units()) / 2.0;
    //double threshold = std::min(leaf_node->kernel_boundry->max_lon - leaf_node->kernel_boundry->min_lon,
    //                            leaf_node->kernel_boundry->max_lat - leaf_node->kernel_boundry->min_lat) / 2.0;
    //threshold = 0;

    for(unsigned i = 0; i < leaf_node->neighbors.size(); i++) {
        if(leaf_node->neighbors[i].second)
            continue;

        /* compute shared boundry of leaf_node and its neighbor */
#ifdef DEBUG
        assert(leaf_node->neighbors[i].first->processing_units_id.size() == 1);
        assert(iter <= 0xFF);
        assert(leaf_node->processing_units_id[0] <= 0xFFF);
        assert(leaf_node->neighbors[i].first->processing_units_id[0] <= 0xFFF);
#endif
        compute_common_boundry(leaf_node, leaf_node->neighbors[i].first, &common_boundary_head, &common_boundary_tail,
                               &cyclic_common_boundary_head, &cyclic_common_boundary_tail);
        
        local_checksums[i] = 0;
        if(common_boundary_head.x != PDLN_DOUBLE_INVALID_VALUE) {
            checksum = leaf_node->triangulation->calculate_triangles_intersected_checksum(common_boundary_head, common_boundary_tail, threshold);

            if(checksum == PDLN_CHECKSUM_FALSE)
                local_checksums[i] = PDLN_CHECKSUM_FALSE;
            else
                local_checksums[i] ^= checksum;
        }

        if(cyclic_common_boundary_head.x != PDLN_DOUBLE_INVALID_VALUE) {
            checksum = leaf_node->triangulation->calculate_triangles_intersected_checksum(cyclic_common_boundary_head, cyclic_common_boundary_tail, threshold);

            if(checksum == PDLN_CHECKSUM_FALSE || local_checksums[i] == PDLN_CHECKSUM_FALSE)
                local_checksums[i] = PDLN_CHECKSUM_FALSE;
            else
                local_checksums[i] ^= checksum;
        }

        if(common_boundary_head.x != PDLN_DOUBLE_INVALID_VALUE || cyclic_common_boundary_head.x != PDLN_DOUBLE_INVALID_VALUE) {
            waiting_list->push_back(new MPI_Request);
            MPI_Isend(&local_checksums[i], 1, MPI_UNSIGNED, processing_info->get_processing_unit(leaf_node->neighbors[i].first->processing_units_id[0])->process_id, 
                      PDLN_SET_TAG(leaf_node->processing_units_id[0], leaf_node->neighbors[i].first->processing_units_id[0], iter),
                      processing_info->get_mpi_comm(), waiting_list->back()); 

            waiting_list->push_back(new MPI_Request);
            MPI_Irecv(&remote_checksums[i], 1, MPI_UNSIGNED, processing_info->get_processing_unit(leaf_node->neighbors[i].first->processing_units_id[0])->process_id,
                      PDLN_SET_TAG(leaf_node->neighbors[i].first->processing_units_id[0], leaf_node->processing_units_id[0], iter),
                      processing_info->get_mpi_comm(), waiting_list->back());
        }
        else
            remote_checksums[i] = 0;
    }
}


bool Delaunay_grid_decomposition::are_checksums_identical(Search_tree_node *leaf_node, unsigned *local_checksums, unsigned *remote_checksums)
{
    if(leaf_node->neighbors.size() == 0) {
        //printf("are_checksums_identical fase ending\n");
        return false;
    }

    bool ok = true;
    for(unsigned i = 0; i < leaf_node->neighbors.size(); i++) {
        if(leaf_node->neighbors[i].second) {
            //printf("[%d] neighbor %d already done\n", leaf_node->processing_units_id[0], leaf_node->neighbors[i].first->processing_units_id[0]);
            continue;
        }

        if(local_checksums[i] == remote_checksums[i] && local_checksums[i] != PDLN_CHECKSUM_FALSE) {
            //printf("[%d] neighbor %d done, %d vs %d\n", leaf_node->processing_units_id[0], leaf_node->neighbors[i].first->processing_units_id[0], local_checksums[i], remote_checksums[i]);
            leaf_node->neighbors[i].second = true;
        }
        else {
            //printf("[%d] neighbor %d not , %d vs %d\n", leaf_node->processing_units_id[0], leaf_node->neighbors[i].first->processing_units_id[0], local_checksums[i], remote_checksums[i]);
            ok = false;
        }
    }

    return ok;
}


int Delaunay_grid_decomposition::assign_polars(bool assign_south_polar, bool assign_north_polar)
{
    double *child_points_coord[4];
    int *child_points_index[2];
    Boundry child_boundry[2];
    int child_num_points[2];
    vector<int> child_proc_id[2];
    Midline midline;

    if(!assign_south_polar && !assign_north_polar)
        return 0;
    
    for(int i = 0; i < 4; i++)
        child_points_coord[i] = new double[search_tree_root->num_kernel_points];
    child_points_index[0] = new int[search_tree_root->num_kernel_points];
    child_points_index[1] = new int[search_tree_root->num_kernel_points];

    if(assign_south_polar) {
        //printf("South polar need rotating.\n");
        search_tree_root->decompose_by_processing_units_number(workloads, child_points_coord, child_points_index, child_num_points, child_boundry, child_proc_id, PDLN_DECOMPOSE_SPOLAR_MODE);
        if(child_boundry[0].max_lat > PDLN_SPOLAR_MAX_LAT || search_tree_root->processing_units_id.size() == 1) {
            midline.type = PDLN_LAT;
            midline.value = PDLN_SPOLAR_MAX_LAT;
            search_tree_root->split_local_points(midline, child_points_coord, child_points_index, child_num_points);

            if(child_num_points[0] < min_num_points_per_chunk) {
                for(int i = 0; i < 4; i++)
                    delete[] child_points_coord[i];
                delete[] child_points_index[0];
                delete[] child_points_index[1];
                printf("assign_polars fault, %d vs %d\n", child_num_points[0], min_num_points_per_chunk);
                return 1;
            }
            child_boundry[0].max_lat = PDLN_SPOLAR_MAX_LAT;
            child_boundry[1].min_lat = PDLN_SPOLAR_MAX_LAT;
            workloads[search_tree_root->processing_units_id[0]] -= child_num_points[0];
            child_proc_id[1].insert(child_proc_id[1].begin(), search_tree_root->processing_units_id[0]);
        }
        search_tree_root->children[0] = alloc_search_tree_node(search_tree_root, child_points_coord,   child_points_index[0], child_num_points[0],
                                                               child_boundry[0], child_proc_id[0], PDLN_NODE_TYPE_SPOLAR);

        search_tree_root->children[1] = alloc_search_tree_node(search_tree_root, child_points_coord+2, child_points_index[1], child_num_points[1],
                                                               child_boundry[1], child_proc_id[1], PDLN_NODE_TYPE_COMMON);

        current_tree_node = search_tree_root->children[1];
        
        if(have_local_processing_units_id(search_tree_root->children[0]->processing_units_id))
            local_leaf_nodes.push_back(search_tree_root->children[0]);

        search_tree_root->children[0]->load_polars_info();
    }
    
    if(assign_north_polar) {
        //printf("Nouth polar need rotating.\n");
        current_tree_node->decompose_by_processing_units_number(workloads, child_points_coord, child_points_index, child_num_points, child_boundry, child_proc_id, PDLN_DECOMPOSE_NPOLAR_MODE);
        if(child_boundry[1].min_lat < PDLN_NPOLAR_MIN_LAT || current_tree_node->processing_units_id.size() == 1) {
            midline.type = PDLN_LAT;
            midline.value = PDLN_NPOLAR_MIN_LAT;
            current_tree_node->split_local_points(midline, child_points_coord, child_points_index, child_num_points);

            if(child_num_points[1] < min_num_points_per_chunk) {
                for(int i = 0; i < 4; i++)
                    delete[] child_points_coord[i];
                delete[] child_points_index[0];
                delete[] child_points_index[1];
                printf("assign_polars fault, %d vs %d\n", child_num_points[1], min_num_points_per_chunk);
                return 1;
            }
            child_boundry[0].max_lat = PDLN_NPOLAR_MIN_LAT;
            child_boundry[1].min_lat = PDLN_NPOLAR_MIN_LAT;
            workloads[search_tree_root->processing_units_id.back()] -= child_num_points[1];
            child_proc_id[0].push_back(search_tree_root->processing_units_id.back());
        }
        if(search_tree_root->children[1] != NULL)
            delete search_tree_root->children[1];
        search_tree_root->children[1] = alloc_search_tree_node(search_tree_root, child_points_coord,   child_points_index[0], child_num_points[0], child_boundry[0], 
                                                               child_proc_id[0], PDLN_NODE_TYPE_COMMON);

        search_tree_root->children[2] = alloc_search_tree_node(search_tree_root, child_points_coord+2, child_points_index[1], child_num_points[1], child_boundry[1],
                                                               child_proc_id[1], PDLN_NODE_TYPE_NPOLAR);

        current_tree_node = search_tree_root->children[1];

        if(have_local_processing_units_id(search_tree_root->children[2]->processing_units_id))
            local_leaf_nodes.push_back(search_tree_root->children[2]);

        search_tree_root->children[2]->load_polars_info();
    }

    for(int i = 0; i < 4; i++)
        delete[] child_points_coord[i];
    delete[] child_points_index[0];
    delete[] child_points_index[1];
    return 0;
}


void Delaunay_grid_decomposition::assign_cyclic_grid_for_single_processing_unit()
{
    double *child_points_coord[4];
    int *child_points_index[2];
    Boundry child_boundry[2];
    int child_num_points[2];
    vector<int> child_proc_id[2];
    Midline midline;

    for(int i = 0; i < 4; i++)
        child_points_coord[i] = new double[current_tree_node->num_kernel_points];
    child_points_index[0] = new int[search_tree_root->num_kernel_points];
    child_points_index[1] = new int[search_tree_root->num_kernel_points];

    assert(current_tree_node->processing_units_id.size() == 1);
    midline.type = PDLN_LON;
    midline.value = 180.0;
    child_boundry[0] = child_boundry[1] = *current_tree_node->kernel_boundry;
    child_boundry[0].max_lon = child_boundry[1].min_lon = 180.0;
    child_proc_id[0] = child_proc_id[1] = current_tree_node->processing_units_id;
    current_tree_node->split_local_points(midline, child_points_coord, child_points_index, child_num_points);
    workloads[current_tree_node->processing_units_id[0]] -= child_num_points[0] + child_num_points[1];

    current_tree_node->children[0] = alloc_search_tree_node(current_tree_node, child_points_coord,   child_points_index[0], child_num_points[0], child_boundry[0],
                                                                  child_proc_id[0], PDLN_NODE_TYPE_COMMON);

    current_tree_node->children[2] = alloc_search_tree_node(current_tree_node, child_points_coord+2, child_points_index[1], child_num_points[1], child_boundry[1],
                                                                  child_proc_id[1], PDLN_NODE_TYPE_COMMON);

    assert(have_local_processing_units_id(child_proc_id[0]));
    local_leaf_nodes.push_back(current_tree_node->children[0]);
    local_leaf_nodes.push_back(current_tree_node->children[2]);

    for(int i = 0; i < 4; i++)
        delete[] child_points_coord[i];
    delete[] child_points_index[0];
    delete[] child_points_index[1];
}


void Delaunay_grid_decomposition::decompose_with_fixed_longitude(double fixed_lon)
{
    double *child_points_coord[4];
    int *child_points_index[2];
    Boundry child_boundry[2];
    int child_num_points[2];
    vector<int> child_proc_id[2];

    for(int i = 0; i < 4; i++)
        child_points_coord[i] = new double[current_tree_node->num_kernel_points];
    child_points_index[0] = new int[search_tree_root->num_kernel_points];
    child_points_index[1] = new int[search_tree_root->num_kernel_points];

    current_tree_node->decompose_by_fixed_longitude(fixed_lon, workloads, child_points_coord, child_points_index, child_num_points, child_boundry, child_proc_id);

    current_tree_node->children[0] = alloc_search_tree_node(current_tree_node, child_points_coord,   child_points_index[0], child_num_points[0], child_boundry[0],
                                                            child_proc_id[0], PDLN_NODE_TYPE_COMMON);

    current_tree_node->children[2] = alloc_search_tree_node(current_tree_node, child_points_coord+2, child_points_index[1], child_num_points[1], child_boundry[1],
                                                            child_proc_id[1], PDLN_NODE_TYPE_COMMON);

    for(int i = 0; i < 4; i++)
        delete[] child_points_coord[i];
    delete[] child_points_index[0];
    delete[] child_points_index[1];
}


//TODO: get faster
bool Delaunay_grid_decomposition::have_local_processing_units_id(vector<int> chunk_id)
{
    for(int j = 0; j < processing_info->get_num_local_proc_processing_units(); j++)
        for(unsigned int i = 0; i < chunk_id.size(); i++)
            if(chunk_id[i] == processing_info->get_local_proc_common_id()[j])
                return true;

    return false;
}

int Delaunay_grid_decomposition::generate_grid_decomposition()
{

    initialze_workload();
    current_tree_node = search_tree_root;

    double min_lon, max_lon, min_lat, max_lat;
    bool is_non_monotonic;

    grid_info_mgr->get_grid_boundry(original_grid, &min_lon, &max_lon, &min_lat, &max_lat);
    is_non_monotonic = min_lon > max_lon;

    if(assign_polars(std::abs(min_lat - -90.0) < PDLN_FLOAT_EQ_ERROR, std::abs(max_lat -  90.0) < PDLN_FLOAT_EQ_ERROR))
        return 1;

    if(is_cyclic && current_tree_node->processing_units_id.size() == 1) {
        assign_cyclic_grid_for_single_processing_unit();
        return 0;
    }

    decompose_common_node_recursively(current_tree_node);
    return 0;
}


void Search_tree_node::load_polars_info()
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
        }
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
        }
    }
}


int Delaunay_grid_decomposition::expand_tree_node_boundry(Search_tree_node* tree_node, double expanding_ratio)
{
    Boundry old_boundry = *tree_node->expanded_boundry;
    Boundry new_boundry = *tree_node->expanded_boundry;

    if(tree_node->node_type == PDLN_NODE_TYPE_COMMON)
        new_boundry = new_boundry * expanding_ratio;
    else if(tree_node->node_type == PDLN_NODE_TYPE_SPOLAR)
        new_boundry.max_lat += (new_boundry.max_lat - new_boundry.min_lat) * expanding_ratio * 2;
    else if(tree_node->node_type == PDLN_NODE_TYPE_NPOLAR)
        new_boundry.min_lat -= (new_boundry.max_lat - new_boundry.min_lat) * expanding_ratio * 2;

    new_boundry.legalize(search_tree_root->kernel_boundry, is_cyclic);

    int rank;
    MPI_Comm_rank(processing_info->get_mpi_comm(), &rank);
    //printf("[%d] old boundary: %lf, %lf, %lf, %lf\n", rank, old_boundry.min_lon, old_boundry.max_lon, old_boundry.min_lat, old_boundry.max_lat);
    //printf("[%d] new boundary: %lf, %lf, %lf, %lf\n", rank, new_boundry.min_lon, new_boundry.max_lon, new_boundry.min_lat, new_boundry.max_lat);
    if(processing_info->get_num_total_processing_units() > 4 &&
       new_boundry.max_lon - new_boundry.min_lon > (search_tree_root->kernel_boundry->max_lon - search_tree_root->kernel_boundry->min_lon) * 0.75 &&
       new_boundry.max_lat - new_boundry.min_lat > (search_tree_root->kernel_boundry->max_lat - search_tree_root->kernel_boundry->min_lat) * 0.75) {
        printf("expanged to the max1\n");
        return -1;
    }
    if(new_boundry == *search_tree_root->kernel_boundry || new_boundry.max_lon - new_boundry.min_lon > 360.0) {
        printf("expanged to the max2\n");
        return -1;
    }

    add_halo_points(tree_node, &old_boundry, &new_boundry);

    return 0;
}


void Delaunay_grid_decomposition::add_halo_points(Search_tree_node* dst_tree_node, Boundry* inner_boundary, Boundry* outer_boundary)
{
    double *expanded_points_coord[2];
    int *expanded_index;
    int num_points_found;

    expanded_points_coord[0] = new double[search_tree_root->num_kernel_points]; //FIXME: buf too large
    expanded_points_coord[1] = new double[search_tree_root->num_kernel_points];
    expanded_index = new int[search_tree_root->num_kernel_points];

    vector<Search_tree_node*> leaf_nodes_found = search_points_in_halo(inner_boundary, outer_boundary, expanded_points_coord, expanded_index, &num_points_found);
    dst_tree_node->add_expanded_points(expanded_points_coord, expanded_index, num_points_found);
    dst_tree_node->add_neighbors(leaf_nodes_found);

    delete[] expanded_points_coord[0];
    delete[] expanded_points_coord[1];
    delete[] expanded_index;
}


vector<Search_tree_node*> Delaunay_grid_decomposition::search_points_in_halo(Boundry* inner_boundary, Boundry* outer_boundary,
                                                                             double *coord_values[2], int *global_idx,
                                                                             int *num_points_found)
{
    vector<Search_tree_node*> leaf_nodes_found;

    *num_points_found = 0;

    if(*inner_boundary == *outer_boundary)
        return leaf_nodes_found;

    search_down_for_points_in_halo(search_tree_root, inner_boundary, outer_boundary, leaf_nodes_found, coord_values, global_idx, num_points_found);

    return leaf_nodes_found;
}


void Delaunay_grid_decomposition::search_down_for_points_in_halo(Search_tree_node *node, Boundry *inner_boundary,
                                                                 Boundry *outer_boundary, vector<Search_tree_node*> &leaf_nodes_found,
                                                                 double *coord_values[2], int *global_idx,
                                                                 int *num_points_found)
{
    double *child_points_coord[4];
    int *child_points_index[2];
    int child_num_points[2];
    Boundry child_boundry[2];
    vector<int> child_proc_id[2];

    assert(node->processing_units_id.size() > 0);

    Boundry region = *outer_boundary;
    if(node->processing_units_id.size() == 1) {
        if(do_two_regions_overlap(region, *node->kernel_boundry) ||
           do_two_regions_overlap(Boundry(region.min_lon + 360.0, region.max_lon + 360.0, region.min_lat, region.max_lat), *node->kernel_boundry) ||
           do_two_regions_overlap(Boundry(region.min_lon - 360.0, region.max_lon - 360.0, region.min_lat, region.max_lat), *node->kernel_boundry)) {
            leaf_nodes_found.push_back(node);
            node->search_points_in_halo(inner_boundary, outer_boundary, coord_values, global_idx, num_points_found);
        }
        return;
    }

    if(node->children[0] == NULL && node->children[1] == NULL && node->children[2] == NULL) {
        for(int i = 0; i < 4; i++)
            child_points_coord[i] = new double[node->num_kernel_points];
        child_points_index[0] = new int[node->num_kernel_points];
        child_points_index[1] = new int[node->num_kernel_points];
        
        node->decompose_by_processing_units_number(workloads, child_points_coord, child_points_index, child_num_points, child_boundry, child_proc_id, PDLN_DECOMPOSE_COMMON_MODE);
        assert(child_proc_id[0].size() > 0);

        node->children[0] = alloc_search_tree_node(node, child_points_coord,   child_points_index[0], child_num_points[0],
                                                   child_boundry[0], child_proc_id[0], PDLN_NODE_TYPE_COMMON);

        node->children[2] = alloc_search_tree_node(node, child_points_coord+2, child_points_index[1], child_num_points[1],
                                                   child_boundry[1], child_proc_id[1], PDLN_NODE_TYPE_COMMON);

        assert(node->children[0]->processing_units_id.size() > 0);
        assert(node->children[2]->processing_units_id.size() > 0);

        for(int i = 0; i < 4; i++)
            delete[] child_points_coord[i];
        delete[] child_points_index[0];
        delete[] child_points_index[1];
        //TODO: optimize new delete
    }

    for(int i = 0; i < 3; i ++)
        if(node->children[i] != NULL) {
            if(do_two_regions_overlap(region, *node->children[i]->kernel_boundry) ||
               do_two_regions_overlap(Boundry(region.min_lon + 360.0, region.max_lon + 360.0, region.min_lat, region.max_lat), *node->children[i]->kernel_boundry) ||
               do_two_regions_overlap(Boundry(region.min_lon - 360.0, region.max_lon - 360.0, region.min_lat, region.max_lat), *node->children[i]->kernel_boundry)) {
                search_down_for_points_in_halo(node->children[i], inner_boundary, outer_boundary, leaf_nodes_found, coord_values, global_idx, num_points_found);
            }
        }
}



void Search_tree_node::search_points_in_halo(Boundry *inner_boundary, Boundry *outer_boundary,
                                             double *coord_values[2], int *global_idx,
                                             int *num_points_found)
{
    if(*kernel_boundry <= *inner_boundary)
        return;

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
    for(int j = 0; j < num_kernel_points; j++) {
        if (is_coordinate_in_halo(points_coord[PDLN_LON][j], points_coord[PDLN_LAT][j], inner_boundary, outer_boundary)) {
            coord_values[PDLN_LON][*num_points_found] = points_coord[PDLN_LON][j];
            coord_values[PDLN_LAT][*num_points_found] = points_coord[PDLN_LAT][j];
            global_idx[(*num_points_found)++] = points_global_index[j];
            continue;
        }
        if (is_coordinate_in_halo(points_coord[PDLN_LON][j], points_coord[PDLN_LAT][j], &l_inner, &l_outer)) {
            coord_values[PDLN_LON][*num_points_found] = points_coord[PDLN_LON][j] + 360.0;
            coord_values[PDLN_LAT][*num_points_found] = points_coord[PDLN_LAT][j];
            global_idx[(*num_points_found)++] = points_global_index[j];
            continue;
        }
        if (is_coordinate_in_halo(points_coord[PDLN_LON][j], points_coord[PDLN_LAT][j], &r_inner, &r_outer)) {
            coord_values[PDLN_LON][*num_points_found] = points_coord[PDLN_LON][j] - 360.0;
            coord_values[PDLN_LAT][*num_points_found] = points_coord[PDLN_LAT][j];
            global_idx[(*num_points_found)++] = points_global_index[j];
            continue;
        }
    }
}

inline bool Search_tree_node::is_coordinate_in_halo(double x, double y, Boundry *inner, Boundry *outer)
{
    return !(x < inner->max_lon && x >= inner->min_lon && y < inner->max_lat && y >= inner->min_lat) &&
           (x < outer->max_lon && x >= outer->min_lon && y < outer->max_lat && y >= outer->min_lat);
}

bool Search_tree_node::check_if_all_outer_edge_out_of_kernel_boundry(Boundry *origin_grid_boundry, bool is_cyclic)
{
    double left, right, top, bot;
    if(triangulation == NULL)
        return false;

    top = fabs(kernel_boundry->max_lat - origin_grid_boundry->max_lat) < PDLN_FLOAT_EQ_ERROR ?
          (kernel_boundry->max_lat+kernel_boundry->min_lat)/2.0 : kernel_boundry->max_lat; //Make the value small enough to pass the check. This is not strict. May need FIXME.
    bot = fabs(kernel_boundry->min_lat - origin_grid_boundry->min_lat) < PDLN_FLOAT_EQ_ERROR ?
          (kernel_boundry->max_lat+kernel_boundry->min_lat)/2.0 : kernel_boundry->min_lat;
    if(is_cyclic) {
        left = kernel_boundry->max_lon;
        right = kernel_boundry->min_lon;
    }
    else {
        left = fabs(kernel_boundry->max_lon - origin_grid_boundry->max_lon) < PDLN_FLOAT_EQ_ERROR ?
               (kernel_boundry->max_lon+kernel_boundry->min_lon)/2.0 : kernel_boundry->max_lon;
        right = fabs(kernel_boundry->min_lon - origin_grid_boundry->min_lon) < PDLN_FLOAT_EQ_ERROR ?
                (kernel_boundry->max_lon+kernel_boundry->min_lon)/2.0 : kernel_boundry->min_lon;
    } 
    //int rank;
    //MPI_Comm_rank(processing_info->get_mpi_comm(), &rank);
    //printf("[%d] checking boundary: %lf, %lf, %lf, %lf\n", rank, top, bot, left, right);

    //return triangulation->check_if_all_outer_edge_out_of_region(left, right, bot, top);
    if(triangulation->check_if_all_outer_edge_out_of_region(left, right, bot, top))
        return true;
    else {
        printf("checkfailed for %d\n", processing_units_id[0]);
        return false;
    }
}


bool Delaunay_grid_decomposition::is_polar_node(Search_tree_node *node) const
{
    return node != NULL && node->node_type != PDLN_NODE_TYPE_COMMON;
}


#define PDLN_MAX_NUM_NEIGHBORS 128
int Delaunay_grid_decomposition::generate_trianglulation_for_local_decomp()
{
    bool *is_local_leaf_node_finished;
    unsigned **local_leaf_checksums, **remote_leaf_checksums;

    is_local_leaf_node_finished = new bool[local_leaf_nodes.size()]();
    local_leaf_checksums = new unsigned*[local_leaf_nodes.size()];
    remote_leaf_checksums = new unsigned*[local_leaf_nodes.size()];

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
    double expanding_ratio = DEFAULT_EXPANGDING_RATIO;
    while(iter < 10) {
        int ret = 0;
        int all_ret;
        for(unsigned i = 0; i < local_leaf_nodes.size(); i++)
            if(!is_local_leaf_node_finished[i]) {
                ret |= expand_tree_node_boundry(local_leaf_nodes[i], expanding_ratio);

                if (is_polar_node(search_tree_root->children[0]) || is_polar_node(search_tree_root->children[2])) // TODO: can be more accurate
                    local_leaf_nodes[i]->project_grid();
            }

        MPI_Allreduce(&ret, &all_ret, 1, MPI_UNSIGNED, MPI_LOR, processing_info->get_mpi_comm());
        if(all_ret) {
            all_finished = false;
            break;
        }
        /* TODO: allgather ret */

        #pragma omp parallel for
        for(unsigned i = 0; i < local_leaf_nodes.size(); i++) {
            if(!is_local_leaf_node_finished[i]) {

                timeval start, end;
                gettimeofday(&start, NULL);

                local_leaf_nodes[i]->generate_local_triangulation(is_cyclic);

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
                is_local_leaf_node_finished[i] = are_checksums_identical(local_leaf_nodes[i], local_leaf_checksums[i], remote_leaf_checksums[i]) &&
                                                 local_leaf_nodes[i]->check_if_all_outer_edge_out_of_kernel_boundry(search_tree_root->kernel_boundry, is_cyclic);
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
/* 
 * Return Values: 0 - Success
 *                1 - Failed in expanding
 *                2 x Fail, polar decomp's expanded_boundry exceeded threshold
 *                2 - Fail, expanding loop too many times
 */


int Delaunay_grid_decomposition::generate_trianglulation_for_whole_grid()
{
    assert(false);
    return 0;
}


void Delaunay_grid_decomposition::plot_local_triangles(const char *perfix)
{
    for(unsigned int i = 0; i < local_leaf_nodes.size(); i++) {
        char filename[64];
        snprintf(filename, 64, "%s_%d.png", perfix, local_leaf_nodes[i]->processing_units_id[0]);
        local_leaf_nodes[i]->triangulation->plot_into_file(filename, local_leaf_nodes[i]->kernel_boundry->min_lon,
                                                                     local_leaf_nodes[i]->kernel_boundry->max_lon,
                                                                     local_leaf_nodes[i]->kernel_boundry->min_lat,
                                                                     local_leaf_nodes[i]->kernel_boundry->max_lat);
    }
}


void Delaunay_grid_decomposition::print_tree_node_info_recursively(Search_tree_node *node)
{
    if(node->processing_units_id.size() == 1){
        printf("[ID%d]x[ST-Info] LEAF %p\n", local_leaf_nodes[0]->processing_units_id[0], node);
        return;
    }

    printf("[ID%d]x[ST-Info] %p: %p, %p, %p\n", local_leaf_nodes[0]->processing_units_id[0], node, node->children[0], node->children[1], node->children[2]);
    for(int i = 0; i < 3; i ++)
        if(node->children[i])
            print_tree_node_info_recursively(node->children[i]);
}


void Delaunay_grid_decomposition::print_whole_search_tree_info()
{
    printf("[ID%d]x[ST-Info] ROOT %p\n", local_leaf_nodes[0]->processing_units_id[0], search_tree_root);
    print_tree_node_info_recursively(search_tree_root);
}


void Delaunay_grid_decomposition::save_ordered_triangles_into_file(Triangle_Transport *triangles, int num_triangles)
{
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
    int num_different_triangles = i + 1;
    
    char file_fmt[] = "log/global_triangles_%d";
    char filename[64];
    snprintf(filename, 64, file_fmt, processing_info->get_num_total_processing_units());
    FILE *fp = fopen(filename, "w");
    for(i = 0; i < num_different_triangles; i++)
        fprintf(fp, "%d, %d, %d, (%lf, %lf), (%lf, %lf), (%lf, %lf)\n", triangles[i].v[0].id, triangles[i].v[1].id, triangles[i].v[2].id,
                                                                        triangles[i].v[0].x, triangles[i].v[0].y,
                                                                        triangles[i].v[1].x, triangles[i].v[1].y,
                                                                        triangles[i].v[2].x, triangles[i].v[2].y);
    fclose(fp);

    char file_fmt2[] = "log/image_global_triangles_%d";
    snprintf(filename, 64, file_fmt2, processing_info->get_num_total_processing_units());
    plot_triangles_info_file(filename, triangles, num_different_triangles);
}


#define PDLN_MERGE_TAG_MASK 0x0200
void Delaunay_grid_decomposition::merge_all_triangles()
{
    int local_buf_len = 0;
    int num_triangles, num_local_triangles;
    Triangle_Transport *local_triangles;

    /* Let n be the number of points, if there are b vertices on the convex hull,
     * then any triangulation of the points has at most 2n − 2 − b triangles,
     * plus one exterior face */
    for(unsigned int i = 0; i < local_leaf_nodes.size(); i++)
        local_buf_len += local_leaf_nodes[i]->num_kernel_points * 3; 

    local_triangles = new Triangle_Transport[local_buf_len];
    num_local_triangles = 0;
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
        Triangle_Transport *remote_triangles;
        MPI_Status status;

        for(int i = 1; i < processing_info->get_num_total_processes(); i++)
            MPI_Recv(&num_remote_triangles[i], 1, MPI_INT, i, PDLN_MERGE_TAG_MASK, processing_info->get_mpi_comm(), &status);
        for(int i = 1; i < processing_info->get_num_total_processes(); i++) {
            //assert(num_remote_triangles[i] > min_num_points_per_chunk/2);
            remote_buf_len += num_remote_triangles[i];
        }
        remote_triangles = new Triangle_Transport[remote_buf_len + num_local_triangles];

        int count = 0;
        for(int i = 1; i < processing_info->get_num_total_processes(); i++) {
            MPI_Recv(remote_triangles + count, num_remote_triangles[i] * sizeof(Triangle_Transport), MPI_CHAR, i, PDLN_MERGE_TAG_MASK, processing_info->get_mpi_comm(), &status);
            int tmp_count;
            MPI_Get_count(&status, MPI_CHAR, &tmp_count);
            
            //char filename[64];
            //snprintf(filename, 64, "log/process_local_triangles%d", i);
            //plot_triangles_info_file(filename, remote_triangles+count, tmp_count/sizeof(Triangle_Transport));

            //for(int j = 0; j < tmp_count/sizeof(Triangle_Transport); j++)
            //    printf("%d, %d, %d\n", (remote_triangles + count)[j].v[0].id, (remote_triangles + count)[j].v[1].id, (remote_triangles + count)[j].v[2].id);
            //printf("==============\n");
#ifdef DEBUG
            assert(tmp_count % sizeof(Triangle_Transport) == 0);
#endif
            count += tmp_count / sizeof(Triangle_Transport);
        }
        assert(count == remote_buf_len);
        memcpy(remote_triangles + remote_buf_len, local_triangles, num_local_triangles * sizeof(Triangle_Transport));
        save_ordered_triangles_into_file(remote_triangles, remote_buf_len + num_local_triangles);
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


Grid_info_manager::Grid_info_manager()
{
    gen_three_polar_grid();
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
int Grid_info_manager::get_polar_points(int grid_id, char polar)
{
    if(polar == 'S')
        return 0;
    if(polar == 'N')
        return 0;
    return 0;
}
