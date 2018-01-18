#ifndef NETCDF_UTILS_H
#define NETCDF_UTILS_H

#include "netcdf.h"

void report_nc_error(int rcode)
{
    if(rcode != NC_NOERR)
        printf("Netcdf error: %s\n", nc_strerror(rcode));
}

void read_file_field(const char *file_name, const char *field_name, void **data_array_ptr, int *num_dims, int **dim_size_ptr, int *field_size)
{
    int i, ncfile_id, rcode, variable_id, *dim_ids, *dim_size;
    long total_size;
    size_t dim_len;
    nc_type nc_var_type;
    float *data_array;
    double *data_array_l;

    *data_array_ptr = NULL;
    *num_dims = 0;
    *dim_size_ptr = NULL;
    *field_size = -1;
    
    rcode = nc_open(file_name, NC_NOWRITE, &ncfile_id);
    report_nc_error(rcode);
    rcode = nc_inq_varid(ncfile_id, field_name, &variable_id);
    if (rcode != NC_NOERR) {
        rcode = nc_close(ncfile_id);
        report_nc_error(rcode);
        return;
    }
    rcode = nc_inq_varndims(ncfile_id, variable_id, num_dims);
    report_nc_error(rcode);
    dim_ids = new int [*num_dims];
    dim_size = new int [*num_dims];
    rcode = nc_inq_vardimid(ncfile_id, variable_id, dim_ids);
    report_nc_error(rcode);
    for (i = 0; i < *num_dims; i ++) {
        rcode = nc_inq_dimlen(ncfile_id, dim_ids[i], &dim_len);
        report_nc_error(rcode);
        dim_size[i] = dim_len;
    }
    rcode = nc_inq_vartype(ncfile_id, variable_id, &nc_var_type);
    report_nc_error(rcode);

    total_size = 1;
    for (i = 0; i < *num_dims; i ++)
        total_size *= dim_size[i];
    data_array = new float [total_size];
    rcode = nc_get_var(ncfile_id, variable_id, data_array);

    data_array_l = new double [total_size];
    for(i = 0; i < total_size; i ++)
        data_array_l[i] = data_array[i];

    delete data_array;
    *data_array_ptr = data_array_l;
    *dim_size_ptr = dim_size;
    if (total_size > 0)
        *field_size = total_size;

    delete [] dim_ids;

    rcode = nc_close(ncfile_id);
    report_nc_error(rcode);
}

#endif
