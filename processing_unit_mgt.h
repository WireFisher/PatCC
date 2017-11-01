/***************************************************************
  *  Copyright (c) 2013, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Haoyu Yang,
  *  supervised by Dr. Li Liu. If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/

#ifndef PROCESSING_UNIT_MGT_H
#define PROCESSING_UNIT_MGT_H

#include "mpi.h"
#include <map>
#include <vector>
using std::vector;
using std::map;

class Processing_unit;


class Processing_unit {
public:
    friend class Processing_resource;
    //char *hostname;
    unsigned int hostname_checksum;
    int process_id;
    int thread_id;
    int common_id;
    Processing_unit(unsigned int checksum, int proc_id, int thread_id): hostname_checksum(checksum),
                                                                        process_id(proc_id),
                                                                        thread_id(thread_id) {};
};


class Processing_resource {
private:
    int component_id;
    map <unsigned int, vector <Processing_unit*> > computing_nodes;
    Processing_unit **processing_units; //(All threads_info sorted by common_id)
    int num_total_processing_units;
    int *local_proc_common_id;
    int num_local_proc_processing_units;
    
    int identify_processing_units_by_hostname();

public:
    Processing_resource();
    ~Processing_resource();
    void pick_out_active_processing_units(int, bool*);
    Processing_unit* get_processing_unit(int common_id) { return this->processing_units[common_id]; };
    Processing_unit** get_processing_units() { return this->processing_units; };
    int get_num_total_processing_units(){ return num_total_processing_units; };
    int* get_local_proc_common_id(){ return local_proc_common_id; };
    int get_num_local_proc_processing_units(){ return num_local_proc_processing_units; };
    void print_all_nodes_info();
};


class Processing_resource_mgt {
private:
    map <int, Processing_resource*> all_processing_info;
public:
    Processing_resource_mgt();
    ~Processing_resource_mgt();
    Processing_resource* get_processing_info(int component_id) {return all_processing_info[component_id]; };
};


class Process_thread_manager
{
public:
    virtual ~Process_thread_manager();
    virtual void get_hostname(char*, int);
    virtual int get_mpi_rank();
    virtual int get_mpi_size();
    virtual MPI_Comm get_mpi_comm();
    virtual int get_openmp_size();
    virtual int allgather(void*, int, MPI_Datatype, void*, int, MPI_Datatype, MPI_Comm);
};

extern Process_thread_manager *process_thread_mgr;
#endif
