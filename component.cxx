#include "component.h"
#define DEFAULT_MIN_NUM_POINTS 100
Grid* Component::search_grid_by_id(int id)
{
    for(int i = 0; i < this->grids.size(); i++)
        if(grids[i]->get_grid_id() == id)
            return grids[i];

    return NULL;
}

int Grid::generate_delaunay_trianglulation(Processing_resource *proc_resource)
{
    int min_num_points_per_chunk, ret;
    bool do_sequentially;
    if(this->delaunay_triangulation)
        return 0;

    min_num_points_per_chunk = DEFAULT_MIN_NUM_POINTS;
    do {
        min_num_points_per_chunk *= 2;
        if(this->delaunay_triangulation)
            delete this->delaunay_triangulation;
        this->delaunay_triangulation = new Delaunay_grid_decomposition(this->grid_id, proc_resource, min_num_points_per_chunk);
        if(this->delaunay_triangulation->generate_grid_decomposition()) {
            do_sequentially = true;
            break;
        }
        ret = this->delaunay_triangulation->generate_trianglulation_for_local_decomp();
        /* Return Values: 0 - success
         *                1 - fail, normal decomp's expanded_boundry exceeded -/+90 
         *                2 - fail, polar  decomp's expanded_boundry exceeded threshold */
        if(ret == 2) {
            do_sequentially = true;
            break;
        }
    } while(ret);

    if(do_sequentially) { 
        if(this->delaunay_triangulation)
            delete this->delaunay_triangulation;
        this->delaunay_triangulation = new Delaunay_grid_decomposition(this->grid_id, NULL, 0);
        this->delaunay_triangulation->generate_trianglulation_for_whole_grid();
    }
}

int Component::generate_delaunay_trianglulation(int grid_id)
{
    Grid *operating_grid;
    operating_grid = this->search_grid_by_id(grid_id);
    if(!operating_grid)
        return -1;

    if(!this->proc_resource)
        this->proc_resource = new Processing_resource();
    //this->proc_resource->print_all_nodes_info();

    if(operating_grid->have_delaunay_trianglulation())
        return 0;

    /* grid pretreatment */

    operating_grid->generate_delaunay_trianglulation(this->proc_resource);
}

Grid_info_manager *grid_info_mgr;
Process_thread_manager *process_thread_mgr;

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    process_thread_mgr = new Process_thread_manager();
    grid_info_mgr = new Grid_info_manager();

    Component* comp;
    comp = new Component(0);
    comp->register_grid(new Grid(1));
    comp->generate_delaunay_trianglulation(1);

    delete process_thread_mgr;
    delete grid_info_mgr;
    MPI_Finalize();
}
