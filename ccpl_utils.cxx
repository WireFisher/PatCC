#include "ccpl_utils.h"
#include <cmath>
#include "common_utils.h"
#include <algorithm>

/* reference: https://www.uwgb.edu/dutchs/structge/sphproj.htm */
/* Original formula:
 *      X = 2 * x / (1 + y)
 *      Z = 2 * z / (1 + y)
 */
static inline void lonlat2xyz(double lon, double lat, PAT_REAL *x, PAT_REAL *y, PAT_REAL *z)
{
    *x = cosl(DEGREE_TO_RADIAN((PAT_REAL)lat)) * sinl(DEGREE_TO_RADIAN((PAT_REAL)lon));
    *y = sinl(DEGREE_TO_RADIAN((PAT_REAL)lat));
    *z = cosl(DEGREE_TO_RADIAN((PAT_REAL)lat)) * cosl(DEGREE_TO_RADIAN((PAT_REAL)lon));
}


static inline void stereo_lonlat2xyz(double lon_p, double lat_p, double lon_t, double lat_t, PAT_REAL *x, PAT_REAL *y, PAT_REAL *z)
{
    PAT_REAL p_x, p_y, p_z;
    PAT_REAL t_x, t_y, t_z;

    lonlat2xyz(lon_p, lat_p, &p_x, &p_y, &p_z);
    lonlat2xyz(lon_t, lat_t, &t_x, &t_y, &t_z);

    PAT_REAL s  = 2.0 / (t_x*(t_x+p_x) + t_y*(t_y+p_y) + t_z*(t_z+p_z));

    *x = s * p_x + (s - 1.0) * t_x;
    *y = s * p_y + (s - 1.0) * t_y;
    *z = s * p_z + (s - 1.0) * t_z;
}


static inline void fast_stereo_lonlat2xyz(double lon_p, double lat_p, PAT_REAL t_x, PAT_REAL t_y, PAT_REAL t_z, PAT_REAL *x, PAT_REAL *y, PAT_REAL *z)
{
    PAT_REAL p_x, p_y, p_z;

    lonlat2xyz(lon_p, lat_p, &p_x, &p_y, &p_z);

    PAT_REAL s  = 2.0 / (t_x*(t_x+p_x) + t_y*(t_y+p_y) + t_z*(t_z+p_z));

    *x = s * p_x + (s - 1.0) * t_x;
    *y = s * p_y + (s - 1.0) * t_y;
    *z = s * p_z + (s - 1.0) * t_z;
}


static inline void normalize_vector(PAT_REAL *x, PAT_REAL *y, PAT_REAL *z)
{
    PAT_REAL length = std::sqrt(*x * *x + *y * *y + *z * *z);
    *x /= length;
    *y /= length;
    *z /= length;
}


static inline void calculate_unit_vectors(double t_lon, double t_lat, 
                                          PAT_REAL *v1_x, PAT_REAL *v1_y, PAT_REAL *v1_z,
                                          PAT_REAL *v2_x, PAT_REAL *v2_y, PAT_REAL *v2_z)
{
    PAT_REAL t_x, t_y, t_z;

    lonlat2xyz(t_lon, t_lat, &t_x, &t_y, &t_z);

    PAT_REAL axis_x = 1.0;
    PAT_REAL axis_y = 0.0;
    PAT_REAL axis_z = 0.0;

    *v1_x = t_y * axis_z - t_z * axis_y;
    *v1_y = t_z * axis_x - t_x * axis_z;
    *v1_z = t_x * axis_y - t_y * axis_x;

    normalize_vector(v1_x, v1_y, v1_z);

    *v2_x = t_y * *v1_z - t_z * *v1_y;
    *v2_y = t_z * *v1_x - t_x * *v1_z;
    *v2_z = t_x * *v1_y - t_y * *v1_x;

    normalize_vector(v2_x, v2_y, v2_z);
}


const PAT_REAL multip_ratio = 100;
void calculate_stereographic_projection(double p_lon, double p_lat, double t_lon, double t_lat, double &X, double &Y)
{
    PAT_REAL q_x, q_y, q_z;
    PAT_REAL e1_x, e1_y, e1_z;
    PAT_REAL e2_x, e2_y, e2_z;

    stereo_lonlat2xyz(p_lon, p_lat, t_lon, t_lat, &q_x, &q_y, &q_z);
    calculate_unit_vectors(t_lon, t_lat, &e1_x, &e1_y, &e1_z, &e2_x, &e2_y, &e2_z);

    X = (q_x * e1_x + q_y * e1_y + q_z * e1_z) * multip_ratio;
    Y = (q_x * e2_x + q_y * e2_y + q_z * e2_z) * multip_ratio;
}


/* 
 * p: point need to be projected
 * t: tangent point
 * e1 & e2: unit vector on the projection plane
 * X & Y: coordinate values of point p on the projection plane
 */
void fast_stereographic_projection(double p_lon, double p_lat,
                                   PAT_REAL t_x, PAT_REAL t_y, PAT_REAL t_z,
                                   PAT_REAL e1_x, PAT_REAL e1_y, PAT_REAL e1_z,
                                   PAT_REAL e2_x, PAT_REAL e2_y, PAT_REAL e2_z,
                                   double &X, double &Y)
{
    PAT_REAL q_x, q_y, q_z;

    fast_stereo_lonlat2xyz(p_lon, p_lat, t_x, t_y, t_z, &q_x, &q_y, &q_z);

    X = (q_x * e1_x + q_y * e1_y + q_z * e1_z) * multip_ratio;
    Y = (q_x * e2_x + q_y * e2_y + q_z * e2_z) * multip_ratio;
}


bool point_in_circle(double lon_deg, double lat_deg, double circle_data[3])
{
    double lon1 = DEGREE_TO_RADIAN(lon_deg);
    double lat1 = DEGREE_TO_RADIAN(lat_deg);
    double lon2 = DEGREE_TO_RADIAN(circle_data[0]);
    double lat2 = DEGREE_TO_RADIAN(circle_data[1]);
    return acos(sin(lat1) * sin(lat2) + cos(lat1) * cos(lat2) * cos(lon1 - lon2)) <= DEGREE_TO_RADIAN(circle_data[2]);
}
