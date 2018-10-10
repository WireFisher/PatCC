#include "ccpl_utils.h"
#include <cmath>
#include "common_utils.h"
#include <algorithm>

/* reference: https://www.uwgb.edu/dutchs/structge/sphproj.htm */
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

    lonlat2xyz(lon_p, lat_p, &p_x, &p_y, &p_z);
    lonlat2xyz(lon_t, lat_t, &t_x, &t_y, &t_z);

    double s  = 2.0 / (t_x*(t_x+p_x) + t_y*(t_y+p_y) + t_z*(t_z+p_z));

    *x = s * p_x + (s - 1.0) * t_x;
    *y = s * p_y + (s - 1.0) * t_y;
    *z = s * p_z + (s - 1.0) * t_z;
}


static inline void stereo_lonlat2xyz(double lon_p, double lat_p, double t_x, double t_y, double t_z, double *x, double *y, double *z)
{
    double p_x, p_y, p_z;

    lonlat2xyz(lon_p, lat_p, &p_x, &p_y, &p_z);

    double s  = 2.0 / (t_x*(t_x+p_x) + t_y*(t_y+p_y) + t_z*(t_z+p_z));

    *x = s * p_x + (s - 1.0) * t_x;
    *y = s * p_y + (s - 1.0) * t_y;
    *z = s * p_z + (s - 1.0) * t_z;
}


static inline void normalize_vector(double *x, double *y, double *z)
{
    double length = std::sqrt(*x * *x + *y * *y + *z * *z);
    *x /= length;
    *y /= length;
    *z /= length;
}


static inline void calculate_unit_vectors(double lon_tan, double lat_tan, 
                                          double *v1_x, double *v1_y, double *v1_z,
                                          double *v2_x, double *v2_y, double *v2_z)
{
    double t_x, t_y, t_z;

    lonlat2xyz(lon_tan, lat_tan, &t_x, &t_y, &t_z);
    double min_dir = std::min(std::abs(t_x),std::min(std::abs(t_y),std::abs(t_z)));

    double axis_x = 1.0;
    double axis_y = 0.0;
    double axis_z = 0.0;

    *v1_x = t_y * axis_z - t_z * axis_y;
    *v1_y = t_z * axis_x - t_x * axis_z;
    *v1_z = t_x * axis_y - t_y * axis_x;

    normalize_vector(v1_x, v1_y, v1_z);

    *v2_x = t_y * *v1_z - t_z * *v1_y;
    *v2_y = t_z * *v1_x - t_x * *v1_z;
    *v2_z = t_x * *v1_y - t_y * *v1_x;

    normalize_vector(v2_x, v2_y, v2_z);
}


static inline void calculate_unit_vectors(double t_x, double t_y, double t_z,
                                          double *v1_x, double *v1_y, double *v1_z,
                                          double *v2_x, double *v2_y, double *v2_z)
{
    double min_dir = std::min(std::abs(t_x),std::min(std::abs(t_y),std::abs(t_z)));

    double axis_x = 1.0;
    double axis_y = 0.0;
    double axis_z = 0.0;

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
}


void fast_stereographic_projection(double lon_original, double lat_original,
                                   double x_tan, double y_tan, double z_tan,
                                   double e1_x, double e1_y, double e1_z,
                                   double e2_x, double e2_y, double e2_z,
                                   double &X, double &Y)
{
    double q_x, q_y, q_z;

    stereo_lonlat2xyz(lon_original, lat_original, x_tan, y_tan, z_tan, &q_x, &q_y, &q_z);

    X = q_x * e1_x + q_y * e1_y + q_z * e1_z;
    Y = q_x * e2_x + q_y * e2_y + q_z * e2_z;

    X *= 50.0;
    Y *= 50.0;
}
