/***************************************************************
  *  Copyright (c) 2019, Tsinghua University.
  *  This is a source file of PatCC.
  *  This file was initially finished by Dr. Li Liu and
  *  Haoyu Yang. If you have any problem,
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#ifndef NETCDF_UTILS_H
#define NETCDF_UTILS_H

extern void report_nc_error(int rcode);
extern void read_file_field_as_float(const char *, const char *, void **, int *, int **, int *);
extern void read_file_field_as_double(const char *, const char *, void **, int *, int **, int *, char* unit=NULL);
extern void read_file_field_as_int(const char *, const char *, void **, int *, int **, int *, char* unit=NULL);

#endif
