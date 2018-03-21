#include "component.h"
#include <cassert>
#include <cstdio>
#include <sys/time.h>
#define DEFAULT_MIN_NUM_POINTS 100
#define MULTIPLICATION_COEFFICIENT 2

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
        gettimeofday(&start, NULL);
        if(delaunay_triangulation->generate_grid_decomposition()) {
            do_sequentially = true;
            break;
        }

        gettimeofday(&end, NULL);
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        printf("[%3d] Grid decomposition: %ldms\n", rank, ((end.tv_sec - start.tv_sec) * 1000000 + (end.tv_usec - start.tv_usec)) / 1000);

        gettimeofday(&start, NULL);
        ret = delaunay_triangulation->generate_trianglulation_for_local_decomp();
        gettimeofday(&end, NULL);
        //printf("[%3d] Trianglulation for local decomp: %ldms\n", rank, ((end.tv_sec - start.tv_sec) * 1000000 + (end.tv_usec - start.tv_usec)) / 1000);
        /* Return Values: 0 - success
         *                1 - fail, normal decomp's expanded_boundry exceeded -/+90 (expanding fail)
         *                2 - fail, polar  decomp's expanded_boundry exceeded threshold */
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
        this->delaunay_triangulation->generate_trianglulation_for_whole_grid();
    }
    return 0;
}

void Grid::plot_triangles_into_file()
{
    this->delaunay_triangulation->plot_local_triangles("log/chunk");
}

void Grid::merge_all_triangles()
{
    delaunay_triangulation->merge_all_triangles();
}


Component::Component(int id): component_id(id)
{
    proc_resource = NULL;
}


void Component::grid_pretreatment(int grid_id)
{
    double min_lon, max_lon, min_lat, max_lat;
    double **coord_values;
    int num_points;

    int num_npolar, num_spolar;
    double max_lat_without_polar, min_lat_without_polar;

    grid_info_mgr->get_grid_boundry(grid_id, &min_lon, &max_lon, &min_lat, &max_lat);
    coord_values = grid_info_mgr->get_grid_coord_values(grid_id);
    num_points = grid_info_mgr->get_grid_num_points(grid_id);

    /* deal with non-monotonic grid */
    if(min_lon > max_lon) {
        assert(min_lon >= 0 && min_lon <= 360);
        assert(max_lon >= 0 && max_lon <= 360);

        double split_line = (min_lon + max_lon) * 0.5;
        for(int i = 0; i < num_points; i++) {
            if(coord_values[PDLN_LON][i] > split_line) coord_values[PDLN_LON][i] -= 360;
        }
        min_lon -= 360;
        grid_info_mgr->set_grid_boundry(grid_id, min_lon, max_lon, min_lat, max_lat);
    }

    /* deal with redundent polar points */
    num_npolar = 0;
    num_spolar = 0;
    max_lat_without_polar = -1e10;
    min_lat_without_polar = 1e10;
    for(int i = 0; i < num_points; i++) {
        if(coord_values[PDLN_LAT][i] == 90.0) {
            polars_index.push_back(std::pair<int, bool>(i, 0));
            num_npolar++;
        }
        else if(coord_values[PDLN_LAT][i] == -90.0) {
            polars_index.push_back(std::pair<int, bool>(i, 1));
            num_spolar++;
        }
        else {
            if(max_lat_without_polar < coord_values[PDLN_LAT][i]) max_lat_without_polar = coord_values[PDLN_LAT][i];
            if(min_lat_without_polar > coord_values[PDLN_LAT][i]) min_lat_without_polar = coord_values[PDLN_LAT][i];
        }
    }

    for(unsigned i = 0; i < polars_index.size(); i++) {
        if(polars_index[i].second == 0)
            coord_values[PDLN_LON][polars_index[i].first] = (90.0 + max_lat_without_polar) * 0.5;
        else
            coord_values[PDLN_LON][polars_index[i].first] = (-90.0 + min_lat_without_polar) * 0.5;
    }
}


Component::~Component()
{
    delete proc_resource;
    for(unsigned i = 0; i < grids.size(); i ++)
        delete grids[i];
}


void Component::generate_delaunay_trianglulation(int grid_id)
{
    Grid *operating_grid;
    operating_grid = this->search_grid_by_id(grid_id);
    assert(operating_grid);

    if(proc_resource == NULL) {
        delete proc_resource;
        proc_resource = new Processing_resource();
    }
    //this->proc_resource->print_all_nodes_info();

    grid_pretreatment(grid_id);

    operating_grid->generate_delaunay_trianglulation(this->proc_resource);
    operating_grid->plot_triangles_into_file();
    //operating_grid->merge_all_triangles();
}

Grid_info_manager *grid_info_mgr;
Process_thread_manager *process_thread_mgr;

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    timeval start, end;
    gettimeofday(&start, NULL);
    process_thread_mgr = new Process_thread_manager();
    grid_info_mgr = new Grid_info_manager();

    Component* comp;
    comp = new Component(0);
    comp->register_grid(new Grid(1));
    comp->generate_delaunay_trianglulation(1);

    delete process_thread_mgr;
    delete grid_info_mgr;
    delete comp;
    gettimeofday(&end, NULL);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if(rank == 0)
        printf("Total time: %ldms\n", ((end.tv_sec - start.tv_sec) * 1000000 + (end.tv_usec - start.tv_usec)) / 1000);
    MPI_Finalize();
}
