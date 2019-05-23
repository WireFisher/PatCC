#ifndef NETCDF_UTILS_H
#define NETCDF_UTILS_H

extern void report_nc_error(int rcode);
extern void read_file_field_as_float(const char *, const char *, void **, int *, int **, int *);
extern void read_file_field_as_double(const char *, const char *, void **, int *, int **, int *, char* unit=NULL);
extern void read_file_field_as_int(const char *, const char *, void **, int *, int **, int *, char* unit=NULL);

#endif
