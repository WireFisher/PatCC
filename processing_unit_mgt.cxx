/***************************************************************
  *  Copyright (c) 2013, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Haoyu Yang,
  *  supervised by Dr. Li Liu. If you have any problem,
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/

#include "processing_unit_mgt.h"
#include <unistd.h>
#include <map>
#include <cstdio>
#include <cstring>
#include <cassert>
#include <omp.h>

#define MAX_HOSTNAME_LEN 32

/*
Workload_info::Workload_info(int id, int num_processing_units): grid_id(id), size(num_processing_units)
{
    assert(size > 0);
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
*/

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


Processing_resource::Processing_resource() {
    char hostname[MAX_HOSTNAME_LEN];
    int local_proc_id, num_procs, num_local_threads;
    int *num_threads_per_process;
    unsigned int local_hostname_checksum, *hostname_checksum_per_process;
    MPI_Comm comm;

    this->processing_units = NULL;
    this->component_id = -1;

    this->num_total_processing_units = 0;
    this->num_local_proc_processing_units = 0;

    process_thread_mgr->get_hostname(hostname, MAX_HOSTNAME_LEN);
    local_hostname_checksum = BKDRHash(hostname, MAX_HOSTNAME_LEN);
    //printf("%s\nhash:%d\n", hostname, local_hostname_checksum);
    local_proc_id = process_thread_mgr->get_mpi_rank();
    num_procs = process_thread_mgr->get_mpi_size();
    comm = process_thread_mgr->get_mpi_comm();
    num_local_threads = process_thread_mgr->get_openmp_size();
    
    assert(num_procs > 0);
    if(num_local_threads <= 0)
        assert(false);

    this->local_proc_common_id= new int[num_local_threads];
    num_threads_per_process = new int[num_procs];
    hostname_checksum_per_process = new unsigned int[num_procs];

    process_thread_mgr->allgather(&num_local_threads, 1, MPI_INT, num_threads_per_process, 1, MPI_INT, comm);
    process_thread_mgr->allgather(&local_hostname_checksum, 1, MPI_UNSIGNED, hostname_checksum_per_process, 1, MPI_UNSIGNED, comm);
    //MPI_Allgather();
    //MPI_Allgather(&local_hostname_checksum, 1, MPI_UNSIGNED, hostname_checksum_per_process, 1, MPI_UNSIGNED, comm);
    
    /* assume that "num_threads_per_process" and "hostname_checksum_per_process" are sorted by process id */
    for(int i=0; i < num_procs; i ++)
        //if(computing_nodes->find(hostname_checksum_per_process[i]) == computing_nodes->end())
        for(int j=0; j < num_threads_per_process[i]; j ++) 
            computing_nodes[hostname_checksum_per_process[i]].push_back(new Processing_unit(hostname_checksum_per_process[i], i, j));

    this->identify_processing_units_by_hostname();

    MAP_UINT_VECTOR_T::iterator it;
    for(it = computing_nodes.begin(); it != computing_nodes.end(); it ++)
        for(int i = 0; i < it->second.size(); i++)
            if(it->second[i]->process_id == local_proc_id)
                this->local_proc_common_id[this->num_local_proc_processing_units++] = it->second[i]->common_id;

    assert(num_local_proc_processing_units == num_local_threads);
    //printf("num_local_proc_processing_units: %d\n", this->num_local_proc_processing_units);
    //for(int i=0; i < this->num_local_proc_processing_units; i++)
        //printf("local_common_id: %d\n", this->local_proc_common_id[i]);

    delete[] num_threads_per_process;
    delete[] hostname_checksum_per_process;
}


Processing_resource::~Processing_resource() {
}


int Processing_resource::identify_processing_units_by_hostname() {

    MAP_UINT_VECTOR_T::iterator it;
    int i, current_common_id, num_total_proc_units;
    
    current_common_id = 0;
    num_total_proc_units = 0;
    for(it = computing_nodes.begin(); it != computing_nodes.end(); it ++)
        num_total_proc_units += it->second.size();

    this->processing_units = new Processing_unit*[num_total_proc_units];
    for(it = computing_nodes.begin(); it != computing_nodes.end(); it ++)
        for(i = 0; i < it->second.size(); i ++) {
            it->second[i]->common_id = current_common_id;
            this->processing_units[current_common_id++] = it->second[i];
        }
    
    this->num_total_processing_units = current_common_id;
    return current_common_id;
}

/*
Workload_info* Processing_resource::search_grid_workload_info(int grid_id) {

    for(int i = 0; i < this->grids_workload_info.size(); i++)
        if(this->grids_workload_info[i]->grid_id == grid_id)
            return this->grids_workload_info[i];

    return NULL;
}


Workload_info* Processing_resource::search_or_add_grid_workload_info(int grid_id, int num_points) {

    for(int i = 0; i < this->grids_workload_info.size(); i++)
        if(this->grids_workload_info[i]->grid_id == grid_id)
            return this->grids_workload_info[i];

    this->grids_workload_info.push_back(new Workload_info(grid_id, num_points));
    return this->grids_workload_info.back();
}
*/

/*
 * 
 */
void Processing_resource::pick_out_actived_processing_units(int num_total_actived_units, bool* is_actived) {
    double ratio;
    int i;
    int num_actived_units, num_zero_actived_nodes_index;
    int *num_actived_units_per_node, *zero_actived_nodes_index;
    MAP_UINT_VECTOR_T::iterator it;

    ratio = ((double)num_total_actived_units/this->num_total_processing_units);
    num_actived_units = 0; 

    num_actived_units_per_node = new int[computing_nodes.size()];
    zero_actived_nodes_index = new int[computing_nodes.size()];
    num_zero_actived_nodes_index = 0;
    for(it = computing_nodes.begin(), i = 0; it != computing_nodes.end(); it ++, i++) {
        num_actived_units_per_node[i] = it->second.size() * ratio;
        if(num_actived_units_per_node[i] == 0) {
            zero_actived_nodes_index[num_zero_actived_nodes_index++] = i;
        }
        num_actived_units += num_actived_units_per_node[i];
    }

    for(i = 0; i < num_zero_actived_nodes_index; i++) {
        if(num_actived_units >= num_total_actived_units)
            break;
        num_actived_units += ++num_actived_units_per_node[zero_actived_nodes_index[i]];
    }

    it = computing_nodes.begin();
    i = 0;
    while(num_actived_units < num_total_actived_units) {
        if(num_actived_units_per_node[i] < it->second.size()) {
            num_actived_units_per_node[i]++;
            num_actived_units++;
        }
        it ++;
        i ++;
        if(it == computing_nodes.end()) {
            it = computing_nodes.begin();
            i = 0;
        }

    }

    assert(num_actived_units = num_total_actived_units);
    
    for(i = 0; i < this->num_total_processing_units; i++)
        is_actived[i] = false;

    for(it = computing_nodes.begin(), i = 0; it != computing_nodes.end(); it ++, i ++)
        for(int j=0; j < num_actived_units_per_node[i]; j ++)
            is_actived[it->second[j]->common_id] = true;

    //current_workload_info->size_actived = num_actived_units;
    //current_workload_info->update_actived_common_id();
    //for(int i = 0; i < num_actived_units; i ++)
    //    printf("updated_actived_common_id: %d\n", current_workload_info->actived_common_id[i]);
    
    delete[] num_actived_units_per_node;
    delete[] zero_actived_nodes_index;
}

/*
void Workload_info::update_actived_common_id()
{
    int i, j;
    if(actived_common_id)
        delete[] actived_common_id;

    actived_common_id = new int[this->size_actived];
    for(i = 0, j = 0; i < this->size; i++)
        if(this->is_actived[i] == true)
            actived_common_id[j++] = i;

    assert(j == this->size_actived);
}

void Workload_info::update_workloads(int total_workload, vector<int> &ids, int min_workload_for_single_unit)
{
    double old_total_workload = 0.0;
    double unassigned_workload = 0.0;

    for(int i = 0; i < ids.size(); i++)
        old_total_workload += this->workloads[ids[i]];

    for(int i = 0; i < ids.size(); i++)
        this->workloads[ids[i]] = this->workloads[ids[i]] * total_workload / old_total_workload;

    for(int i = 0; i < ids.size(); i++)
        if(this->workloads[ids[i]] < min_workload_for_single_unit) {
            unassigned_workload += this->workloads[ids[i]];
            ids.erase(ids.begin() + i);
        }

    assert(!0.0 > 0);
    if(unassigned_workload > 0)
        for(int i = 0; i < ids.size(); i++)
            this->workloads[ids[i]] += unassigned_workload / ids.size();
}
*/


void Processing_resource::print_all_nodes_info()
{
    for(int i = 0; i < this->num_total_processing_units; i++)
        printf("hostname: %u\nproc_id: %d\nthread_id: %d\ncommon_id: %d\n========\n", this->processing_units[i]->hostname_checksum, this->processing_units[i]->process_id, this->processing_units[i]->thread_id, this->processing_units[i]->common_id);
}

void Process_thread_manager::get_hostname(char *hostname, int len) {
    strncpy(hostname, "default", len);
    return;
}

Process_thread_manager::~Process_thread_manager() {
}

int Process_thread_manager::get_mpi_rank() {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    return rank;
}

int Process_thread_manager::get_mpi_size() {
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    return size;
}

MPI_Comm Process_thread_manager::get_mpi_comm() {
    return MPI_COMM_WORLD;
}

int Process_thread_manager::get_openmp_size() {
    return omp_get_max_threads();
}

int Process_thread_manager::allgather(void *sendbuf, int sendcount, MPI_Datatype sendtype, 
                                      void *recvbuf, int recvcount, MPI_Datatype recvtype,
                                      MPI_Comm comm) {
    MPI_Allgather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm);// Haven't been tested
}
