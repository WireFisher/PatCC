#ifndef PDLN_CCPL_UTILS_H
#define PDLN_CCPL_UTILS_H

#define DEGREE_TO_RADIAN(data)    (data*PI/180.0)
#define RADIAN_TO_DEGREE(data)    (data*180.0/PI)

extern void rotate_sphere_coordinate(double, double, double&, double&);
extern void calculate_orthographic_projection(double, double, double &, double &);
extern void calculate_stereographic_projection(double, double, double, double, double &, double &);

#endif
