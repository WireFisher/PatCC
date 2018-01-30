#include "ccpl_utils.h"
#include <cmath>
#include <cassert>
#include <algorithm>

#define PI (3.1415926535898)

void rotate_sphere_coordinate(double lon_original, double lat_original, double &lon_rotated, double &lat_rotated)
{
    double lon_rotated_radian, lat_rotated_radian, temp1_value, temp2_value;


    temp1_value = cos(DEGREE_TO_RADIAN(lon_original))*cos(DEGREE_TO_RADIAN(lat_original));
    temp2_value = sin(DEGREE_TO_RADIAN(lon_original))*cos(DEGREE_TO_RADIAN(lat_original))/sqrt(1-temp1_value*temp1_value);
    if (temp1_value < -1.0)
        temp1_value = -1.0;
    if (temp1_value > 1.0)
        temp1_value = 1.0;
    if (temp2_value < -1.0)
        temp2_value = -1.0;
    if (temp2_value > 1.0)
        temp2_value = 1.0;    
    lon_rotated_radian = asin(temp2_value);
    lat_rotated_radian = -asin(temp1_value);    
    if (cos(lon_rotated_radian)*cos(lat_rotated_radian)*sin(DEGREE_TO_RADIAN(lat_original)) < 0)
        lon_rotated_radian = PI - lon_rotated_radian;
    if (lon_rotated_radian < 0)
        lon_rotated_radian += 2*PI;
    
    lon_rotated = RADIAN_TO_DEGREE(lon_rotated_radian);
    lat_rotated = RADIAN_TO_DEGREE(lat_rotated_radian);
    
    if (lat_original == 90) {
        lat_rotated = 0;
        lon_rotated = 0;
    }
    else if (lat_original == -90) {
        lat_rotated = 0;
        lon_rotated = 180;
    }

    if (temp1_value == 1.0) {
        lat_rotated = -90;
        lon_rotated = 0;        
    }
    if (temp1_value == -1.0) {
        lat_rotated = 90;
        lon_rotated = 0;        
    }

    lon_rotated = (double) ((float) lon_rotated);
    lat_rotated = (double) ((float) lat_rotated);

    assert(lon_rotated >= 0 && lon_rotated <= 360);
    assert(lat_rotated >= -90 && lat_rotated <= 90);

    if (lon_rotated == 360)
        lon_rotated = 0;
}


/* reference: https://www.uwgb.edu/dutchs/structge/sphproj.htm */
void calculate_orthographic_projection(double lon_original, double lat_original, double &X, double &Z)
{
    X = cos(DEGREE_TO_RADIAN(lat_original)) * sin(DEGREE_TO_RADIAN(lon_original)) * 100;
    Z = cos(DEGREE_TO_RADIAN(lat_original)) * cos(DEGREE_TO_RADIAN(lon_original)) * 100;

    /* To make fake-cyclic triangles' edges be in proper order */
    if(lat_original > 0.0)
        X *= -1;
}

/* Original formula:
 *      X = 2 * x / (1 + y)
 *      Z = 2 * z / (1 + y)
 */

static inline void lonlat2xyz(double lon, double lat, double *x, double *y, double *z)
{
    *x = cos(DEGREE_TO_RADIAN(lat)) * sin(DEGREE_TO_RADIAN(lon));
    *y = sin(DEGREE_TO_RADIAN(lat));
    *z = cos(DEGREE_TO_RADIAN(lat)) * cos(DEGREE_TO_RADIAN(lon));
}


static inline void stereo_lonlat2xyz(double lon_p, double lat_p, double lon_t, double lat_t, double *x, double *y, double *z)
{
    double p_x, p_y, p_z;
    double t_x, t_y, t_z;
    double s;

    lonlat2xyz(lon_p, lat_p, &p_x, &p_y, &p_z);
    lonlat2xyz(lon_t, lat_t, &t_x, &t_y, &t_z);

    s = 2.0 / (t_x*(t_x+p_x) + t_y*(t_y+p_y) + t_z*(t_z+p_z));

    *x = s * p_x + (s - 1.0) * t_x;
    *y = s * p_y + (s - 1.0) * t_y;
    *z = s * p_z + (s - 1.0) * t_z;
}


static inline void normalize_vector(double *x, double *y, double *z)
{
    double length;
    length = std::sqrt(*x * *x + *y * *y + *z * *z);
    *x /= length;
    *y /= length;
    *z /= length;
}


static inline void calculate_unit_vectors(double lon_tan, double lat_tan, 
                                          double *v1_x, double *v1_y, double *v1_z,
                                          double *v2_x, double *v2_y, double *v2_z)
{
    double t_x, t_y, t_z;
    double axis_x, axis_y, axis_z;
    double min_dir;

    lonlat2xyz(lon_tan, lat_tan, &t_x, &t_y, &t_z);
    min_dir = std::min(std::abs(t_x),std::min(std::abs(t_y),std::abs(t_z)));

    if(min_dir == std::abs(t_x)){
        axis_x = 1.0;
        axis_y = 0.0;
        axis_z = 0.0;
    } else if (min_dir == std::abs(t_y)){
        axis_x = 0.0;
        axis_y = 1.0;
        axis_z = 0.0;
    } else if (min_dir == std::abs(t_z)){
        axis_x = 0.0;
        axis_y = 0.0;
        axis_z = 1.0;
    }

    *v1_x = t_y * axis_z - t_z * axis_y;
    *v1_y = t_z * axis_x - t_x * axis_z;
    *v1_z = t_x * axis_y - t_y * axis_x;

    normalize_vector(v1_x, v1_y, v1_z);

    *v2_x = t_y * *v1_z - t_z * *v1_y;
    *v2_y = t_z * *v1_x - t_x * *v1_z;
    *v2_z = t_x * *v1_y - t_y * *v1_x;

    normalize_vector(v2_x, v2_y, v2_z);
}


void calculate_stereographic_projection(double lon_original, double lat_original, double lon_tan, double lat_tan, double &X, double &Y)
{
    double q_x, q_y, q_z;
    double e1_x, e1_y, e1_z;
    double e2_x, e2_y, e2_z;

    stereo_lonlat2xyz(lon_original, lat_original, lon_tan, lat_tan, &q_x, &q_y, &q_z);
    calculate_unit_vectors(lon_tan, lat_tan, &e1_x, &e1_y, &e1_z, &e2_x, &e2_y, &e2_z);

    X = q_x * e1_x + q_y * e1_y + q_z * e1_z;
    Y = q_x * e2_x + q_y * e2_y + q_z * e2_z;

    X *= 50.0;
    Y *= 50.0;
    /* To make fake-cyclic triangles' edges be in proper order */
    //if(lat_original > 0.0)
    //    X *= -1;
}
