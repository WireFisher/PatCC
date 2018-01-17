#ifndef PDLN_CCPL_UTILS_H
#define PDLN_CCPL_UTILS_H

#define DEGREE_TO_RADIAN(data)    (data*PI/180)
#define RADIAN_TO_DEGREE(data)    (data*180/PI)

extern void rotate_sphere_coordinate(double, double, double&, double&);

#endif