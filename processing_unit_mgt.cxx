/***************************************************************
  *  Copyright (c) 2013, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Haoyu Yang,
  *  supervised by Dr. Li Liu. If you have any problem,
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/

#include "mpi.h"
#include "processing_unit_mgt.h"
#include <unistd.h>
#include <map>
#include "testcase.h"
#include "cstdio"
#include <cstring>

#define MAX_HOSTNAME_LEN 128


Workload_info::Workload_info(int id, int num_processing_units): grid_id(id), size(num_processing_units)
{
    //assert size > 0
    this->is_actived = new bool[this->size]();
    this->workloads = new double[this->size]();
    this->actived_common_id = new int[this->size]();
}


Workload_info::~Workload_info()
{
    delete[] is_actived;
    delete[] workloads;
    delete[] actived_common_id;
}


/* BKDR Hash Function */
unsigned int BKDRHash(char *str, unsigned int n)
{
    unsigned int seed = 131; // 31 131 1313 13131 131313 etc..
    unsigned int hash = 0;
    char *head = str;

    while (*str && str - head < n)
    {
        hash = hash * seed + (*str++);
    }

    return (hash & 0x7FFFFFFF);
}


Processing_info::Processing_info() {
    char hostname[MAX_HOSTNAME_LEN];
    int local_proc_id, num_procs, local_thread_id, num_local_threads;
    int *num_threads_per_process;
    unsigned int local_hostname_checksum, *hostname_checksum_per_process;
    MPI_Comm comm;

    processing_units = NULL;
    processing_units_id_actived = NULL;
    component_id = -1;

    num_total_processing_units = 0;
    num_actived_processing_units = 0;
    local_proc_processing_units = NULL;
    num_local_proc_processing_units = 0;

    gethostname(hostname, MAX_HOSTNAME_LEN);
    local_hostname_checksum = BKDRHash(hostname, MAX_HOSTNAME_LEN);
    local_proc_id = get_mpi_rank();
    num_procs = get_mpi_size();
    comm = get_mpi_comm();
    local_thread_id = get_openmp_rank();
    num_local_threads = get_openmp_size();
    
    //assertion num_procs > 0
    num_threads_per_process = new int[num_procs];
    hostname_checksum_per_process = new unsigned int[num_procs];

    MPI_Allgather(&num_local_threads, 1, MPI_INT, num_threads_per_process, 1, MPI_INT, comm);
    MPI_Allgather(&local_hostname_checksum, 1, MPI_UNSIGNED, hostname_checksum_per_process, 1, MPI_UNSIGNED, comm);
    
    /* assume that "num_threads_per_process" and "hostname_checksum_per_process" are sorted by process id */
    for(int i=0; i < num_procs; i ++)
        //if(computing_nodes->find(hostname_checksum_per_process[i]) == computing_nodes->end())
        for(int j=0; j < num_threads_per_process[i]; j ++)
            computing_nodes[hostname_checksum_per_process[i]].push_back(new Processing_unit(hostname_checksum_per_process[i], i, j));

    this->identify_processing_units_by_hostname();

    delete[] num_threads_per_process;
    delete[] hostname_checksum_per_process;
}


Processing_info::~Processing_info() {
}


int Processing_info::identify_processing_units_by_hostname() {

    MAP_UINT_VECTOR_T::iterator it;
    int i, current_common_id, num_total_processing_units;
    
    current_common_id = 0;
    num_total_processing_units = 0;
    for(it = computing_nodes.begin(); it != computing_nodes.end(); it ++)
        num_total_processing_units += it->second.size();

    this->processing_units = new Processing_unit*[num_total_processing_units];
    for(it = computing_nodes.begin(); it != computing_nodes.end(); it ++)
        for(i = 0; i < it->second.size(); i ++) {
            it->second[i]->common_id = current_common_id;
            this->processing_units[current_common_id++] = it->second[i];
        }
    
    this->num_total_processing_units = current_common_id;
    return current_common_id;
}


Workload_info* Processing_info::search_grid_workload_info(int grid_id) {

    for(int i = 0; i < this->grids_workload_info.size(); i++)
        if(this->grids_workload_info[i]->grid_id == grid_id)
            return this->grids_workload_info[i];

    return NULL;
}


Workload_info* Processing_info::search_or_add_grid_workload_info(int grid_id, int num_points) {

    for(int i = 0; i < this->grids_workload_info.size(); i++)
        if(this->grids_workload_info[i]->grid_id == grid_id)
            return this->grids_workload_info[i];

    return this->grids_workload_info.push_back(new Workload_info((grid_id, num_points));
}


/*
 * Return Values 0: Workload info already exists
 */
void Processing_info::pick_out_actived_processing_units(int grid_id, int num, double average_workload) {
    double ratio;
    int i;
    int num_actived_units, num_units;
    int *num_actived_units_per_node;
    MAP_UINT_VECTOR_T::iterator it;
    Workload_info *current_workload_info;

    if(this->search_grid_workload_info(grid_id))
        return 0;

    ratio = ((double)num/this->num_total_processing_units);
    printf("%lf\n", ratio);
    num_actived_units = 0; 

    num_actived_units_per_node = new int[computing_nodes.size()];
    for(it = computing_nodes.begin(), i = 0; it != computing_nodes.end(); it ++, i++) {
        num_actived_units_per_node[i] = it->second.size() * ratio;
        num_actived_units += num_actived_units_per_node[i];
    }
    for(it = computing_nodes.begin(), i = 0; num_actived_units < num; it ++, i ++) {
        if(it == computing_nodes.end()) {
            it = computing_nodes.begin();
            i = 0;
        }
        
        if(num_actived_units_per_node[i] < it->second.size()) {
            num_actived_units_per_node[i];
            num_actived_units++;
        }
    }
    //assert num_actived_units = num
    current_workload_info = this->search_or_add_grid_workload_info(grid_id, this->num_total_processing_units);
    //this->grids_workload_info.push_back(new Workload_info(grid_id, this->num_total_processing_units));
    for(it = computing_nodes.begin(), i = 0; it != computing_nodes.end(); it ++, i ++)
        for(int j=0; j < num_actived_units_per_node[i]; j ++) {
            current_workload_info->is_actived[it->second[j]->common_id] = true;
            current_workload_info->workloads[it->second[j]->common_id]= average_workload;
        }
    this->grids_workload_info.back()->size_actived = num_actived_units;
    this->grids_workload_info.back()->update_actived_common_id();

    delete[] num_actived_units_per_node;
}


void Workload_info::update_actived_common_id()
{
    int i, j;
    if(actived_common_id)
        delete[] actived_common_id;

    actived_common_id = new int[this->size_actived];
    for(i = 0, j = 0; i < this->size; i++)
        if(this->is_actived[i] == true)
            actived_common_id[j++] = i;
    //assert j == this->size_actived
}
