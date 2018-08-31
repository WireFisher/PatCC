#ifndef __COMPONENT__
#define __COMPONENT__

#include "processing_unit_mgt.h"
#include "delaunay_grid_decomposition_mgt.h"

class Grid
{
private:
    int grid_id;
    Delaunay_grid_decomposition *delaunay_triangulation;
public:
    Grid(int id):grid_id(id){ delaunay_triangulation = NULL; };
    ~Grid();
    int get_grid_id(){ return grid_id; };
    int generate_delaunay_trianglulation(Processing_resource*);
    bool have_delaunay_trianglulation(){return delaunay_triangulation != NULL; };
    void merge_all_triangles();
#ifdef OPENCV
    void plot_triangles_into_file();
#endif
};

class Component
{
private:
    int component_id;
    vector<Grid*> grids;
    Processing_resource *proc_resource;

    Grid* search_grid_by_id(int);
    void grid_pretreatment(int);
public:
    Component(int id);
    ~Component();
    void register_grid(Grid* grid){this->grids.push_back(grid); };
    int generate_delaunay_trianglulation(int);
};

#endif
