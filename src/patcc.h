/***************************************************************
  *  Copyright (c) 2019, Tsinghua University.
  *  This is a source file of PatCC.
  *  This file was initially finished by Dr. Li Liu and
  *  Haoyu Yang. If you have any problem,
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#ifndef __COMPONENT__
#define __COMPONENT__

#include "processing_unit_mgt.h"
#include "grid_decomposition.h"

class Grid
{
private:
    int grid_id;
    Delaunay_grid_decomposition *delaunay_triangulation;
public:
    Grid(int id):grid_id(id){ delaunay_triangulation = NULL; };
    ~Grid();
    int get_grid_id(){ return grid_id; };
    int generate_delaunay_trianglulation(Processing_resource*, Grid_info);
    bool have_delaunay_trianglulation(){return delaunay_triangulation != NULL; };
    void merge_all_triangles(bool);
#ifdef OPENCV
    void plot_triangles_into_file();
#endif
};


class Patcc
{
public:
    Patcc(int id);
    ~Patcc();
    void register_grid(Grid* grid){this->grids.push_back(grid); };
    int generate_delaunay_trianglulation(int, bool=false);

private:
    Grid* search_grid_by_id(int);
    void grid_preprocessing(int);

    int component_id;
    vector<Grid*> grids;
    Processing_resource *proc_resource;
    vector<int> shifted_spoles_index;
    vector<int> shifted_npoles_index;
    Grid_info grid_info;
};

#endif
