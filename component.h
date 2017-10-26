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
    Grid(int id):grid_id(id){};
    ~Grid();
    int get_grid_id(){ return grid_id; };
    int generate_delaunay_trianglulation(Processing_resource*);
    bool have_delaunay_trianglulation(){return delaunay_triangulation != NULL; };
};

class Component
{
private:
    int component_id;
    vector<Grid*> grids;
    Processing_resource *proc_resource;

    Grid* search_grid_by_id(int);
public:
    Component(int id):component_id(id){};
    ~Component();
    void register_grid(Grid* grid){this->grids.push_back(grid); };
    void generate_delaunay_trianglulation(int);
};

#endif
