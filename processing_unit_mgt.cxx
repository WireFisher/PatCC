/***************************************************************
  *  Copyright (c) 2013, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Haoyu Yang,
  *  supervised by Dr. Li Liu. If you have any problem,
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/

#include "processing_unit_mgt.h"
#include <unistd.h>
#include <cstdio>
#include <cstring>
#include <omp.h>
#include "common_utils.h"

#include <sched.h>
#include <sys/syscall.h>
#include <sys/sysinfo.h>

#define MAX_HOSTNAME_LEN 32
typedef std::map <unsigned int, vector <Processing_unit*> > MAP_UINT_VECTOR_T;

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
    int *num_threads_per_process;
    unsigned int local_hostname_checksum, *hostname_checksum_per_process;

    this->processing_units = NULL;
    this->component_id = -1;

    this->num_total_processing_units = 0;
    this->num_local_proc_processing_units = 0;

    process_thread_mgr->get_hostname(hostname, MAX_HOSTNAME_LEN);
    local_hostname_checksum = BKDRHash(hostname, MAX_HOSTNAME_LEN);
    local_process_id = process_thread_mgr->get_mpi_rank();
    num_total_processes = process_thread_mgr->get_mpi_size();
    mpi_comm = process_thread_mgr->get_mpi_comm();
    num_local_threads = process_thread_mgr->get_openmp_size();
    
    PDASSERT(num_total_processes > 0);
    if(num_local_threads <= 0)
        PDASSERT(false);

    local_proc_common_id = new int[num_local_threads];
    num_threads_per_process = new int[num_total_processes];
    hostname_checksum_per_process = new unsigned int[num_total_processes];

    process_thread_mgr->allgather(&num_local_threads, 1, MPI_INT, num_threads_per_process, 1, MPI_INT, mpi_comm);
    process_thread_mgr->allgather(&local_hostname_checksum, 1, MPI_UNSIGNED, hostname_checksum_per_process, 1, MPI_UNSIGNED, mpi_comm);
    //MPI_Allgather();
    //MPI_Allgather(&local_hostname_checksum, 1, MPI_UNSIGNED, hostname_checksum_per_process, 1, MPI_UNSIGNED, mpi_comm);
    
    /* assume that "num_threads_per_process" and "hostname_checksum_per_process" are sorted by process id */
    for(int i=0; i < num_total_processes; i ++)
        for(int j=0; j < num_threads_per_process[i]; j ++) 
            computing_nodes[hostname_checksum_per_process[i]].push_back(new Processing_unit(hostname_checksum_per_process[i], i, j));

    identify_processing_units_by_hostname();

    MAP_UINT_VECTOR_T::iterator it;
    for(it = computing_nodes.begin(); it != computing_nodes.end(); it ++)
        for(unsigned int i = 0; i < it->second.size(); i++)
            if(it->second[i]->process_id == local_process_id)
                local_proc_common_id[num_local_proc_processing_units++] = it->second[i]->common_id;

    PDASSERT(num_local_proc_processing_units == num_local_threads);

    delete[] num_threads_per_process;
    delete[] hostname_checksum_per_process;

    set_cpu_affinity();
}


Processing_resource::~Processing_resource() {
    MAP_UINT_VECTOR_T::iterator it;

    delete[] local_proc_common_id;
    for(it = computing_nodes.begin(); it != computing_nodes.end(); it ++)
        it->second.clear();
    computing_nodes.clear();
    delete[] processing_units;
}


void Processing_resource::set_cpu_affinity()
{
    #pragma omp parallel
    {
        cpu_set_t mask;

        int omp_rank = omp_get_thread_num();
        PDASSERT(omp_rank < get_num_local_proc_processing_units());

        int total_cores = get_nprocs();

        CPU_ZERO(&mask);
        CPU_SET(get_local_proc_common_id()[omp_rank]%total_cores, &mask);

        pid_t tid = syscall(__NR_gettid);
        if (sched_setaffinity(tid, sizeof(cpu_set_t), &mask) == -1)
            perror("Error: sched_setaffinity");

        //fprintf(stderr, "Processing Unit %d on CPU %d\n", get_local_proc_common_id()[omp_rank], sched_getcpu());
    }
}

int Processing_resource::identify_processing_units_by_hostname() {

    MAP_UINT_VECTOR_T::iterator it;
    unsigned int i, current_common_id, num_total_proc_units;
    
    current_common_id = 0;
    num_total_proc_units = 0;
    for(it = computing_nodes.begin(); it != computing_nodes.end(); it ++)
        num_total_proc_units += it->second.size();

    processing_units = new Processing_unit*[num_total_proc_units];
    for(it = computing_nodes.begin(); it != computing_nodes.end(); it ++)
        for(i = 0; i < it->second.size(); i ++) {
            it->second[i]->common_id = current_common_id;
            processing_units[current_common_id++] = it->second[i];
        }
    
    num_total_processing_units = current_common_id;
    return current_common_id;
}


void Processing_resource::pick_out_active_processing_units(int num_total_active_units, bool* is_active) {
    double ratio;
    int i;
    int num_active_units;
    int *num_active_units_per_node;
    MAP_UINT_VECTOR_T::iterator it;

    if(num_total_active_units >= num_total_processing_units) {
        for(i = 0; i < num_total_processing_units; i++)
            is_active[i] = true;
        return;
    }

    ratio = ((double)num_total_active_units/num_total_processing_units);
    num_active_units = 0; 

    num_active_units_per_node = new int[computing_nodes.size()];
    for(it = computing_nodes.begin(), i = 0; it != computing_nodes.end(); it ++, i++) {
        num_active_units_per_node[i] = it->second.size() * ratio;
        if(num_active_units_per_node[i] == 0)
            num_active_units_per_node[i] = 1;
        num_active_units += num_active_units_per_node[i];
    }

    for(it = computing_nodes.begin(), i = 0; it != computing_nodes.end(); it++, i++) {
        if(num_active_units == num_total_active_units)
            break;
        if(num_active_units < num_total_active_units) {
            num_active_units_per_node[i]++;
            num_active_units++;
        }
        else {
            num_active_units_per_node[i]--;
            num_active_units--;
        }
#ifdef DEBUG
        PDASSERT(num_active_units_per_node[i] >= 0);
        PDASSERT((unsigned)num_active_units_per_node[i] <= it->second.size());
#endif
    }

    PDASSERT(num_active_units == num_total_active_units);
    
    for(i = 0; i < num_total_processing_units; i++)
        is_active[i] = false;

    for(it = computing_nodes.begin(), i = 0; it != computing_nodes.end(); it ++, i ++)
        for(int j=0; j < num_active_units_per_node[i]; j ++)
            is_active[it->second[j]->common_id] = true;

    delete[] num_active_units_per_node;
}


bool Processing_resource::send_to_local_thread(const void *buf, int count, int size, int src, int dst, int tag)
{
    local_thread_comm.push_back(Thread_comm_packet(buf, count*size, src, dst, tag));
    return true;
}


int Processing_resource::recv_from_local_thread(void *buf, int max_count, int size, int src, int dst, int tag)
{
    for(unsigned int i = 0; i < local_thread_comm.size(); i++)
        if(local_thread_comm[i].src == src && local_thread_comm[i].dest == dst && local_thread_comm[i].tag == tag) {
            memcpy(buf, local_thread_comm[i].buf, std::min(max_count, local_thread_comm[i].len));
            return std::min(max_count, local_thread_comm[i].len);
        }
    return -1;
}


void Processing_resource::print_all_nodes_info()
{
    for(int i = 0; i < num_total_processing_units; i++)
        printf("hostname: %u\nproc_id: %d\nthread_id: %d\ncommon_id: %d\n========\n", processing_units[i]->hostname_checksum, processing_units[i]->process_id, 
                                                                                      processing_units[i]->thread_id, processing_units[i]->common_id);
}

void Process_thread_manager::get_hostname(char *hostname, int len) {
    strncpy(hostname, "default", len);
    return;
}

Process_thread_manager::~Process_thread_manager() {
}

int Process_thread_manager::get_mpi_rank() {
    int rank;
    MPI_Comm_rank(get_mpi_comm(), &rank);
    return rank;
}

int Process_thread_manager::get_mpi_size() {
    int size;
    MPI_Comm_size(get_mpi_comm(), &size);
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
    MPI_Allgather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm);
    return 0;
}
