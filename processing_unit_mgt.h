/***************************************************************
  *  Copyright (c) 2013, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Haoyu Yang,
  *  supervised by Dr. Li Liu. If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/

#ifndef PROCESSING_UNIT_MGT_H
#define PROCESSING_UNIT_MGT_H

#include <map>
#include <vector>
using std::vector;
using std::map;

struct Processing_unit;
typedef map <unsigned int, vector <Processing_unit*> > MAP_UINT_VECTOR_T;

struct Workload_info {
    int grid_id;
    int size;
    bool *is_actived;
    double *workloads;
    int *actived_common_id;
    int size_actived;

    Workload_info(int, int);
    ~Workload_info();
    void update_actived_common_id();
};


/*
struct Workload_info {
    int grid_id;
    bool is_active;
    double workload;
};
*/

struct Processing_unit {
    //char *hostname;
    unsigned int hostname_checksum;
    int process_id;
    int thread_id;
    int common_id;
    //vector<struct Workload_info*> grids_workload_info;
    //struct Workload_info* current_grid_workload_info;
    Processing_unit(unsigned int checksum, int proc_id, int thread_id): hostname_checksum(checksum),
                                                                        process_id(proc_id),
                                                                        thread_id(thread_id) {};
};

/*
struct Computing_node {
    char *hostname;
    vector <struct Processing_unit*> processing_units;
};
*/

class Processing_info {
private:
    int component_id;
    map <unsigned int, vector <Processing_unit*> > computing_nodes;
    vector<Workload_info*> grids_workload_info;
    Processing_unit **processing_units; //(All threads_info sorted by common_id)
    //vector <struct Computing_node*> computing_nodes;
    //int *processing_units_id_actived;
    
    int num_total_processing_units;
    //int num_actived_processing_units;
    Processing_unit **local_proc_processing_units; //(length is number of threads in local process)
    int num_local_proc_processing_units;
    
    int identify_processing_units_by_hostname();
public:
    Processing_info();
    ~Processing_info();
    void pick_out_actived_processing_units(int, int);
    Workload_info* search_grid_workload_info(int);
    Workload_info* search_or_add_grid_workload_info(int, int);
    Processing_unit* get_processing_unit(int common_id) { return this->processing_units[common_id]; };
    int get_num_total_processing_units(){ return num_total_processing_units; };
    //int get_num_actived_processing_units(){ return num_actived_processing_units; };
    Processing_unit** get_local_proc_processing_units(){ return local_proc_processing_units; };
};


class Processing_info_mgt {
private:
    map <int, Processing_info*> all_processing_info;
public:
    Processing_info_mgt();
    ~Processing_info_mgt();
    get_processing_info(int component_id) {return all_processing_info[component_id]; };
}
#endif
