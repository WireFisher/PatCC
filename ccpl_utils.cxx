#include "ccpl_utils.h"
#include <cmath>
#include <cassert>

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


void calculate_orthographic_projection(double lon_original, double lat_original, double &x, double &y)
{
    x = cos(DEGREE_TO_RADIAN(lat_original)) * sin(DEGREE_TO_RADIAN(lon_original)) * 100;
    y = cos(DEGREE_TO_RADIAN(lat_original)) * cos(DEGREE_TO_RADIAN(lon_original)) * 100;

    /* To make fake-cyclic triangles' edges be in proper order */
    if(lat_original > 0.0)
        x *= -1;
}
