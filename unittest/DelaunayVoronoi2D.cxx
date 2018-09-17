#include "mpi.h"
#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include "../delaunay_voronoi_2D.h"

using ::testing::Return;
using ::testing::_;
using ::testing::Invoke;

static double fRand(double fMin, double fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

#ifdef OPENCV
#include <opencv2/opencv.hpp>
void draw_point(cv::Mat img, double x, double y, cv::Scalar scalar)
{
    cv::circle(img, cv::Point(x * 10, y * 10), 6, scalar, -1, 8);
}

void draw_line(cv::Mat img, Edge *e, double min_x, double max_x, double min_y, double max_y, cv::Scalar scalar)
{
    int thickness = 2;
    //int line_type = cv::LINE_8;

    /*
    if(e->head->x < min_x + (max_x - min_x) * 0.1 &&
       e->head->y < min_y + (max_y - min_y) * 0.1 &&
       e->tail->x < min_x + (max_x - min_x) * 0.1 &&
       e->tail->y < min_y + (max_y - min_y) * 0.1 )
        cv::line(img, cv::Point(e->head->x * 100, e->head->y * 100), cv::Point(e->tail->x * 100, e->tail->y * 100), scalar, thickness, line_type);
        */

    cv::line(img, cv::Point(e->head->x * 10, e->head->y * 10), cv::Point(e->tail->x * 10, e->tail->y * 10), scalar, thickness, 8);
}


void write_to_file(cv::Mat mat,const char *filename)
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
#endif


static inline void do_triangulation_plot(int num_points, double *lat_values, double *lon_values, 
                                         double min_lon, double max_lon, double min_lat, double max_lat,
                                         const char* img_path)
{
    std::vector<Edge*> edges;
    Delaunay_Voronoi* delau;

    int *index = new int[num_points];
    for(int i = 0; i < num_points; i++)
        index[i] = i;

    delau = new Delaunay_Voronoi();
    delau->add_points(lon_values, lat_values, index, num_points);
    delau->triangulate();

    delete index;
    edges = delau->get_all_delaunay_edge();

#ifdef OPENCV
    cv::Mat mat = cv::Mat(max_lat*10, max_lon*10, CV_8UC3, cv::Scalar(255, 255, 255));

    for(int i = 0; i < num_points; i++)
        //draw_point(mat, lon_values[i], lat_values[i], cv::Scalar(0x44, 0xBB, 0xBB));
        draw_point(mat, lon_values[i], lat_values[i], cv::Scalar(0, 0, 0));

    write_to_file(mat, "log/input_pp.png");

    for(unsigned int i = 0; i < edges.size(); i++)
        draw_line(mat, edges[i], min_lon, max_lon, min_lat, max_lat, cv::Scalar(0, 0, 0));

    /*
    edges = delau->get_all_legal_delaunay_edge();
    for(unsigned int i = 0; i < edges.size(); i++) {
        draw_line(mat, edges[i], min_lon, max_lon, min_lat, max_lat, cv::Scalar(255, 0, 0));
    }
    */

    write_to_file(mat, img_path);
#endif
    delete delau;
};


static inline void do_triangulation_plot_small_part(int num_points, double *lat_values, double *lon_values, 
                                                    double min_lon, double max_lon, double min_lat, double max_lat,
                                                    const char* img_path)
{
    std::vector<Edge*> edges;
    Delaunay_Voronoi* delau;

    int *index = new int[num_points];
    for(int i = 0; i < num_points; i++)
        index[i] = i;
    delau = new Delaunay_Voronoi();
    delau->add_points(lon_values, lat_values, index, num_points);
    delau->triangulate();

    delete index;
    edges = delau->get_all_delaunay_edge();

#ifdef OPENCV
    cv::Mat mat = cv::Mat::zeros(max_lat*10, max_lon*10, CV_8UC3);

    for(unsigned int i = 0; i < edges.size(); i++)
        draw_line(mat, edges[i], min_lon, max_lon, min_lat, max_lat, cv::Scalar(255, 255, 255));

    //for(int i = 0; i < num_points; i++)
    //    draw_point(mat, lon_values[i], lat_values[i], cv::Scalar(0x44, 0xBB, 0xBB));
    /*
    edges = delau->get_all_legal_delaunay_edge();
    for(unsigned int i = 0; i < edges.size(); i++) {
        draw_line(mat, edges[i], min_lon, max_lon, min_lat, max_lat, cv::Scalar(255, 0, 0));
    }
    */
    write_to_file(mat, img_path);
#endif

    delete delau;
};


TEST(DelaunayTriangulationTest, Rectangle) {
    int len_points = 10;
    int num_points = len_points * len_points;
    double min_lat = 0.0;
    double max_lat = 150.0;
    double min_lon = 0.0;
    double max_lon = 150.0;
    double *lat, *lon;

    lat = new double[num_points];
    lon = new double[num_points];

    for(int i = 0, index = 0; i < len_points; i++)
        for(int j = 0; j < len_points; j++) {
            lat[index] = 15 * i;
            lon[index++] = 15 * j;
    }

    do_triangulation_plot(num_points, lat, lon, min_lon, max_lon, min_lat, max_lat, "log/image1-0.png");

    delete lat;
    delete lon;
};


TEST(DelaunayTriangulationTest, ObliqueRectangle) {
    int len_points = 10;
    int num_points = len_points * len_points;
    double min_lat = 0.0;
    double max_lat = 150.0;
    double min_lon = 0.0;
    double max_lon = 315.0;
    double *lat, *lon;

    lat = new double[num_points];
    lon = new double[num_points];

    for(int i = 0, index = 0; i < len_points; i++)
        for(int j = 0; j < len_points; j++) {
        lat[index] = 15 * i;
        if(i%2 == 0)
            lon[index] = 30 * j;
        else
            lon[index] = 15 + 30 * j;
        index++;
    }

    do_triangulation_plot(num_points, lat, lon, min_lon, max_lon, min_lat, max_lat, "log/image2-0.png");

    delete lat;
    delete lon;
};


TEST(DelaunayTriangulationTest, RandomSmallSize) {
    int num_points = 25;
    double min_lat = 0.0;
    double max_lat = 160.0;
    double min_lon = 0.0;
    double max_lon = 160.0;
    double *lat, *lon;

    lat = new double[num_points];
    lon = new double[num_points];

    for(int i = 0; i < num_points; i++) {
        lat[i] = fRand(min_lat, max_lat);
        lon[i] = fRand(min_lon, max_lon);
    }

    do_triangulation_plot(num_points, lat, lon, min_lon, max_lon, min_lat, max_lat, "log/image3-0.png");

    delete lat;
    delete lon;
};

TEST(DelaunayTriangulationTest, RandomMiddleSize) {
    int num_points = 3600;
    int num_part;
    double min_lat = 0.0;
    double max_lat = 350.0;
    double min_lon = 0.0;
    double max_lon = 350.0;
    double *lat, *lon, *lat_part, *lon_part;

    lat = new double[num_points];
    lon = new double[num_points];

    for(int i = 0; i < num_points; i++) {
        lat[i] = fRand(min_lat, max_lat);
        lon[i] = fRand(min_lon, max_lon);
    }

    do_triangulation_plot(num_points, lat, lon, min_lon, max_lon, min_lat, max_lat, "log/image4-0.png");

    lat_part = new double[num_points];
    lon_part = new double[num_points];

    num_part = 0;
    for(int i = 0; i < num_points; i++)
        if(lat[i] < min_lat+(max_lat-min_lat)*5.0/8.0 && lon[i] < min_lon+(max_lon-min_lon)*5.0/8.0) {
            lat_part[num_part]   = lat[i];
            lon_part[num_part++] = lon[i];
        }
    do_triangulation_plot(num_part, lat_part, lon_part, min_lon, max_lon, min_lat, max_lat, "log/image4-1.png");

    num_part = 0;
    for(int i = 0; i < num_points; i++)
        if(lat[i] < min_lat+(max_lat-min_lat)*5.0/8.0 && lon[i] > min_lon+(max_lon-min_lon)*3.0/8.0) {
            lat_part[num_part]   = lat[i];
            lon_part[num_part++] = lon[i];
        }
    do_triangulation_plot(num_part, lat_part, lon_part, min_lon, max_lon, min_lat, max_lat, "log/image4-2.png");

    num_part = 0;
    for(int i = 0; i < num_points; i++)
        if(lat[i] > min_lat+(max_lat-min_lat)*3.0/8.0 && lon[i] < min_lon+(max_lon-min_lon)*5.0/8.0) {
            lat_part[num_part]   = lat[i];
            lon_part[num_part++] = lon[i];
        }
    do_triangulation_plot(num_part, lat_part, lon_part, min_lon, max_lon, min_lat, max_lat, "log/image4-3.png");

    num_part = 0;
    for(int i = 0; i < num_points; i++)
        if(lat[i] > min_lat+(max_lat-min_lat)*3.0/8.0 && lon[i] > min_lon+(max_lon-min_lon)*3.0/8.0) {
            lat_part[num_part]   = lat[i];
            lon_part[num_part++] = lon[i];
        }
    do_triangulation_plot(num_part, lat_part, lon_part, min_lon, max_lon, min_lat, max_lat, "log/image4-4.png");

    delete lat;
    delete lon;
    delete lat_part;
    delete lon_part;
};


TEST(DelaunayTriangulationTest, RandomLargeSize) {
    //int num_points = 1000000;
    int num_points = 1600000;
    double min_lat = 0.0;
    double max_lat = 350.0;
    double min_lon = 0.0;
    double max_lon = 350.0;
    double *lat, *lon;

    lat = new double[num_points];
    lon = new double[num_points];

    for(int i = 0; i < num_points; i++) {
        lat[i] = fRand(min_lat, max_lat);
        lon[i] = fRand(min_lon, max_lon);
    }

    do_triangulation_plot_small_part(num_points, lat, lon, min_lon, max_lon, min_lat, max_lat, "log/image5-0.png");

    delete lat;
    delete lon;
};
