#include "mpi.h"
#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include "../delaunay_voronoi.h"
#include <opencv2/opencv.hpp>

using ::testing::Return;
using ::testing::_;
using ::testing::Invoke;

static double fRand(double fMin, double fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

void draw_line(cv::Mat img, Edge *e, double min_lon, double max_lon, double min_lat, double max_lat, cv::Scalar scalar)
{
    int thickness = 1;
    int line_type = cv::LINE_8;

    cv::line(img, cv::Point(e->head->x * 10, e->head->y * 10), cv::Point(e->tail->x * 10, e->tail->y * 10), scalar, thickness, line_type);
}

void write_to_file(cv::Mat mat, char *filename)
{
    std::vector<int> compression_params;
    compression_params.push_back(cv::IMWRITE_PNG_COMPRESSION);
    compression_params.push_back(9);
    try {
        cv::imwrite(filename, mat, compression_params);
    }
    catch (cv::Exception& ex) {
        fprintf(stderr, "Exception converting image to PNG format: %s\n", ex.what());
    }
    fprintf(stdout, "Saved PNG file with alpha data.\n");
}

static inline void do_triangulation_plot(int num_points, double *lat_values, double *lon_values, 
                                         double min_lon, double max_lon, double min_lat, double max_lat,
                                         bool *redundant_cell_mark, char* img_path)
{
    std::vector<Edge*> edges;

    Delaunay_Voronoi* delau;

    delau = new Delaunay_Voronoi(num_points, lon_values, lat_values, false, min_lon, max_lon, min_lat, max_lat, redundant_cell_mark);
    edges = delau->get_all_delaunay_edge();

    cv::Mat mat = cv::Mat::zeros(max_lat*10, max_lon*10, CV_8UC3);
    for(unsigned int i = 0; i < edges.size(); i++) {
        draw_line(mat, edges[i], min_lon, max_lon, min_lat, max_lat, cv::Scalar(255, 255, 255));
    }

    /*
    edges = delau->get_all_legal_delaunay_edge();
    for(unsigned int i = 0; i < edges.size(); i++) {
        draw_line(mat, edges[i], min_lon, max_lon, min_lat, max_lat, cv::Scalar(255, 0, 0));
    }
    */

    write_to_file(mat, img_path);
};


TEST(DelaunayTriangulationTest, Rectangle) {
    int len_points = 10;
    int num_points = len_points * len_points;
    int num_part;
    double min_lat = 0.0;
    double max_lat = 150.0;
    double min_lon = 0.0;
    double max_lon = 150.0;
    double *lat, *lon, *lat_part, *lon_part;
    bool *redundant_cell_mark;

    lat = new double[num_points];
    lon = new double[num_points];
    redundant_cell_mark = new bool[num_points];

    for(int i = 0, index = 0; i < len_points; i++)
        for(int j = 0; j < len_points; j++) {
        lat[index] = 15 * i;
        lon[index] = 15 * j;
        redundant_cell_mark[index++] = false;
    }

    do_triangulation_plot(num_points, lat, lon, min_lon, max_lon, min_lat, max_lat, redundant_cell_mark, "log/image0-0.png");

    lat_part = new double[num_points];
    lon_part = new double[num_points];

    /*
    num_part = 0;
    for(int i = 0; i < num_points; i++)
        if(lat[i] < min_lat+(max_lat-min_lat)*5.0/8.0 && lon[i] < min_lon+(max_lon-min_lon)*5.0/8.0) {
            lat_part[num_part]   = lat[i];
            lon_part[num_part++] = lon[i];
        }
    do_triangulation_plot(num_part, lat_part, lon_part, min_lon, max_lon, min_lat, max_lat, redundant_cell_mark, "log/image1-1.png");

    num_part = 0;
    for(int i = 0; i < num_points; i++)
        if(lat[i] < min_lat+(max_lat-min_lat)*5.0/8.0 && lon[i] > min_lon+(max_lon-min_lon)*3.0/8.0) {
            lat_part[num_part]   = lat[i];
            lon_part[num_part++] = lon[i];
        }
    do_triangulation_plot(num_part, lat_part, lon_part, min_lon, max_lon, min_lat, max_lat, redundant_cell_mark, "log/image1-2.png");

    num_part = 0;
    for(int i = 0; i < num_points; i++)
        if(lat[i] > min_lat+(max_lat-min_lat)*3.0/8.0 && lon[i] < min_lon+(max_lon-min_lon)*5.0/8.0) {
            lat_part[num_part]   = lat[i];
            lon_part[num_part++] = lon[i];
        }
    do_triangulation_plot(num_part, lat_part, lon_part, min_lon, max_lon, min_lat, max_lat, redundant_cell_mark, "log/image1-3.png");

    num_part = 0;
    for(int i = 0; i < num_points; i++)
        if(lat[i] > min_lat+(max_lat-min_lat)*3.0/8.0 && lon[i] > min_lon+(max_lon-min_lon)*3.0/8.0) {
            lat_part[num_part]   = lat[i];
            lon_part[num_part++] = lon[i];
        }
    do_triangulation_plot(num_part, lat_part, lon_part, min_lon, max_lon, min_lat, max_lat, redundant_cell_mark, "log/image1-4.png");
    */

    delete lat;
    delete lon;
    delete lat_part;
    delete lon_part;
};


TEST(DelaunayTriangulationTest, ObliqueRectangle) {
    int len_points = 10;
    int num_points = len_points * len_points;
    int num_part;
    double min_lat = 0.0;
    double max_lat = 150.0;
    double min_lon = 0.0;
    double max_lon = 315.0;
    double *lat, *lon, *lat_part, *lon_part;
    bool *redundant_cell_mark;

    lat = new double[num_points];
    lon = new double[num_points];
    redundant_cell_mark = new bool[num_points];

    for(int i = 0, index = 0; i < len_points; i++)
        for(int j = 0; j < len_points; j++) {
        lat[index] = 15 * i;
        if(i%2 == 0)
            lon[index] = 30 * j;
        else
            lon[index] = 15 + 30 * j;
        redundant_cell_mark[index++] = false;
    }

    do_triangulation_plot(num_points, lat, lon, min_lon, max_lon, min_lat, max_lat, redundant_cell_mark, "log/image1-0.png");

    lat_part = new double[num_points];
    lon_part = new double[num_points];

    /*
    num_part = 0;
    for(int i = 0; i < num_points; i++)
        if(lat[i] < min_lat+(max_lat-min_lat)*5.0/8.0 && lon[i] < min_lon+(max_lon-min_lon)*5.0/8.0) {
            lat_part[num_part]   = lat[i];
            lon_part[num_part++] = lon[i];
        }
    do_triangulation_plot(num_part, lat_part, lon_part, min_lon, max_lon, min_lat, max_lat, redundant_cell_mark, "log/image1-1.png");

    num_part = 0;
    for(int i = 0; i < num_points; i++)
        if(lat[i] < min_lat+(max_lat-min_lat)*5.0/8.0 && lon[i] > min_lon+(max_lon-min_lon)*3.0/8.0) {
            lat_part[num_part]   = lat[i];
            lon_part[num_part++] = lon[i];
        }
    do_triangulation_plot(num_part, lat_part, lon_part, min_lon, max_lon, min_lat, max_lat, redundant_cell_mark, "log/image1-2.png");

    num_part = 0;
    for(int i = 0; i < num_points; i++)
        if(lat[i] > min_lat+(max_lat-min_lat)*3.0/8.0 && lon[i] < min_lon+(max_lon-min_lon)*5.0/8.0) {
            lat_part[num_part]   = lat[i];
            lon_part[num_part++] = lon[i];
        }
    do_triangulation_plot(num_part, lat_part, lon_part, min_lon, max_lon, min_lat, max_lat, redundant_cell_mark, "log/image1-3.png");

    num_part = 0;
    for(int i = 0; i < num_points; i++)
        if(lat[i] > min_lat+(max_lat-min_lat)*3.0/8.0 && lon[i] > min_lon+(max_lon-min_lon)*3.0/8.0) {
            lat_part[num_part]   = lat[i];
            lon_part[num_part++] = lon[i];
        }
    do_triangulation_plot(num_part, lat_part, lon_part, min_lon, max_lon, min_lat, max_lat, redundant_cell_mark, "log/image1-4.png");
    */

    delete lat;
    delete lon;
    delete lat_part;
    delete lon_part;
};


TEST(DelaunayTriangulationTest, Random) {
    int num_points = 3600;
    int num_part;
    double min_lat = 0.0;
    double max_lat = 350.0;
    double min_lon = 0.0;
    double max_lon = 350.0;
    double *lat, *lon, *lat_part, *lon_part;
    bool *redundant_cell_mark;

    lat = new double[num_points];
    lon = new double[num_points];
    redundant_cell_mark = new bool[num_points];

    for(int i = 0; i < num_points; i++) {
        lat[i] = fRand(min_lat, max_lat);
        lon[i] = fRand(min_lon, max_lon);
        redundant_cell_mark[i] = false;
    }

    do_triangulation_plot(num_points, lat, lon, min_lon, max_lon, min_lat, max_lat, redundant_cell_mark, "log/image2-0.png");

    lat_part = new double[num_points];
    lon_part = new double[num_points];

    num_part = 0;
    for(int i = 0; i < num_points; i++)
        if(lat[i] < min_lat+(max_lat-min_lat)*5.0/8.0 && lon[i] < min_lon+(max_lon-min_lon)*5.0/8.0) {
            lat_part[num_part]   = lat[i];
            lon_part[num_part++] = lon[i];
        }
    do_triangulation_plot(num_part, lat_part, lon_part, min_lon, max_lon, min_lat, max_lat, redundant_cell_mark, "log/image2-1.png");

    num_part = 0;
    for(int i = 0; i < num_points; i++)
        if(lat[i] < min_lat+(max_lat-min_lat)*5.0/8.0 && lon[i] > min_lon+(max_lon-min_lon)*3.0/8.0) {
            lat_part[num_part]   = lat[i];
            lon_part[num_part++] = lon[i];
        }
    do_triangulation_plot(num_part, lat_part, lon_part, min_lon, max_lon, min_lat, max_lat, redundant_cell_mark, "log/image2-2.png");

    num_part = 0;
    for(int i = 0; i < num_points; i++)
        if(lat[i] > min_lat+(max_lat-min_lat)*3.0/8.0 && lon[i] < min_lon+(max_lon-min_lon)*5.0/8.0) {
            lat_part[num_part]   = lat[i];
            lon_part[num_part++] = lon[i];
        }
    do_triangulation_plot(num_part, lat_part, lon_part, min_lon, max_lon, min_lat, max_lat, redundant_cell_mark, "log/image2-3.png");

    num_part = 0;
    for(int i = 0; i < num_points; i++)
        if(lat[i] > min_lat+(max_lat-min_lat)*3.0/8.0 && lon[i] > min_lon+(max_lon-min_lon)*3.0/8.0) {
            lat_part[num_part]   = lat[i];
            lon_part[num_part++] = lon[i];
        }
    do_triangulation_plot(num_part, lat_part, lon_part, min_lon, max_lon, min_lat, max_lat, redundant_cell_mark, "log/image2-4.png");

    delete lat;
    delete lon;
    delete lat_part;
    delete lon_part;
};
