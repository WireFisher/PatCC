#ifndef OPENCV_UTILS_H
#define OPENCV_UTILS_H

void plot_edge_into_file(const char *filename, double *head_coord[2], double *tail_coord[2], int num_edges,
                         double min_x = 0.0, double max_x = 0.0, double min_y = 0.0, double max_y = 0.0);

void plot_points_info_file(const char *filename, double *x, double *y, int num);


#endif
