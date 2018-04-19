#ifndef NETCDF_UTILS_H
#define NETCDF_UTILS_H

extern void report_nc_error(int rcode);
extern void read_file_field_as_float(const char *file_name, const char *field_name, void **data_array_ptr, int *num_dims, int **dim_size_ptr, int *field_size);
extern void read_file_field_as_double(const char *file_name, const char *field_name, void **data_array_ptr, int *num_dims, int **dim_size_ptr, int *field_size);

#endif
