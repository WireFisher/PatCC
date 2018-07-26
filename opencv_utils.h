#ifndef OPENCV_UTILS_H
#define OPENCV_UTILS_H

#define PDLN_PLOT_COLOR_WHITE   0
#define PDLN_PLOT_COLOR_RED     1
#define PDLN_PLOT_COLOR_DEFAULT PDLN_PLOT_COLOR_WHITE

#define PDLN_PLOT_FILEMODE_NEW     0
#define PDLN_PLOT_FILEMODE_APPEND  1

#define PDLN_PLOT_GLOBAL        0
#define PDLN_PLOT_PROJECTION    1

void plot_edge_into_file(const char *filename, double *head_coord[2], double *tail_coord[2], int num_edges,
                         double min_x = 0.0, double max_x = 0.0, double min_y = 0.0, double max_y = 0.0);
void plot_projected_edge_into_file(const char *filename, double *head_coord[2], double *tail_coord[2], int num_edges, int, int);

void plot_points_into_file(const char *filename, double *x, double *y, int num, int mode = PDLN_PLOT_PROJECTION,
                           double min_x = 0.0, double max_x = 0.0, double min_y = 0.0, double max_y = 0.0);
void plot_rectangle_into_file(const char *filename, double x_min, double x_max, double y_min, double y_lax, int color, int filemode);
void plot_text_into_file(const char*, const char*, double, double, double, double, int);
#endif
