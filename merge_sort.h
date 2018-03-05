/***************************************************************
  *  Copyright (c) 2013, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu and Haoyu Yang. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn or
  *  Haoyu Yang via yanghy16@mails.tsinghua.edu.cn
  ***************************************************************/


#ifndef MERGE_SORT_H
#define MERGE_SORT_H

#include "cstdlib"

#define ASSIGN(a, b, size)                              \
  do                                                    \
    {                                                   \
      size_t __size = (size);                           \
      char *__a = (a), *__b = (b);                      \
      do                                                \
        {                                               \
              *__a++ = *__b++;                          \
        } while (--__size > 0);                         \
    } while (0)

typedef int (*_compare_func) (const void *, const void *);

void merge(char *sorted_values, char *sorted_values_buf, size_t size, size_t lo, size_t mi, size_t hi, _compare_func cmp)
{
    size_t length_b = mi - lo,length_c = hi - mi;
 
    for( size_t i = 0; i < length_b; i++) {
        sorted_values_buf[i] = sorted_values[i + lo];
    }
    char *C = sorted_values + mi;

    for ( size_t i = 0, j = 0, k = 0; j < length_b; ){
        //if (k < length_c && C[k] < sorted_values_buf[j]) {
        if (k < length_c && ((*cmp)(C+k, sorted_values_buf+j) < 0)) {
            //sorted_values[(i++) + lo] = C[k++];
            ASSIGN(sorted_values+i+lo, C+k, size);
            i += size;
            k += size;
        }
        else {
            //sorted_values[(i++) + lo] = sorted_values_buf[j++];
            ASSIGN(sorted_values+i+lo, sorted_values_buf+j, size);
            i += size;
            j += size;
        }
    }
}

//typedef int (*_compare_func) (const void *, const void *, void *);

void merge_sort(void *const pbase, size_t total_elems, size_t size, _compare_func cmp)
{
    size_t scale, i;

    char *base_ptr = (char *) pbase;
    char *sorted_values_buf = (char *) malloc(total_elems * size);

    for(scale = 1; scale < total_elems; scale *= 2)
    {   
        for(i = 0; i < total_elems; i += scale * 2)
        {
            if(i + scale * 2 < total_elems)
                merge(base_ptr, sorted_values_buf, size, i*size, (i+scale)*size, (i+scale*2)*size, cmp);
            else if(i + scale < total_elems)
                merge(base_ptr, sorted_values_buf, size, i*size, (i+scale)*size, total_elems*size, cmp);
        }
    }

    free(sorted_values_buf);
}


#endif

