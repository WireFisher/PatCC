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

void createAlphaMat(cv::Mat &mat)
{
    CV_Assert(mat.channels() == 4);
    for (int i = 0; i < mat.rows; ++i) {
        for (int j = 0; j < mat.cols; ++j) {
            cv::Vec4b& bgra = mat.at<cv::Vec4b>(i, j);
            bgra[0] = UCHAR_MAX; // Blue
            bgra[1] = cv::saturate_cast<uchar>((float (mat.cols - j)) / ((float)mat.cols) * UCHAR_MAX); // Green
            bgra[2] = cv::saturate_cast<uchar>((float (mat.rows - i)) / ((float)mat.rows) * UCHAR_MAX); // Red
            bgra[3] = cv::saturate_cast<uchar>(0.5 * (bgra[1] + bgra[2])); // Alpha
        }
    }
};

template <typename T>
void draw_line(cv::Mat img, Point start, Point end)
{
    int thickness = 1;
    int line_type = cv::LINE_8;

    cv::line(img, cv::Point(start.lon, start.lat), cv::Point(end.lon, end.lat), cv::Scalar(255, 255, 255), thickness, line_type);
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
/*
TEST(DelaunayTriangulationTest, Basic) {
    std::vector<Point<double> > points; 
    
    for(int i = 0; i < 60; i++)
        for(int j = 0; j < 60; j++)
            points.push_back(Point<double>(0 + i * 20, 0 + j * 20));

    Delaunay<double> triangulation;
    triangulation.triangulate(points);
    std::vector<Edge<double> > edges = triangulation.getEdges();

    cv::Mat mat = cv::Mat::zeros(1200, 1200, CV_8UC3);
    //createAlphaMat(mat);
    for(std::vector<Edge<double> >::iterator e = edges.begin(); e != edges.end(); e++)
        draw_line<double>(mat, *e->p1, *e->p2);

    write_to_file(mat, "log/image1.png");
};
*/
static inline void do_triangulation_plot(int num_points, double *lat_values, double *lon_values, 
                                         double min_lon, double max_lon, double min_lat, double max_lat,
                                         bool *redundant_cell_mark, char* img_path)
{
    /*
    Delaunay<double> triangulation;
    triangulation.triangulate(points);
    std::vector<Edge<double> > edges = triangulation.getEdges();
    */
    double *vertex_lon_values, *vertex_lat_values;
    int num_vertexes;
    std::vector<Pure_Edge> edges;

    Delaunay_Voronoi* delau;

    delau = new Delaunay_Voronoi(num_points, lat_values, lon_values, false, min_lon, max_lon, min_lat, max_lat, redundant_cell_mark, &vertex_lon_values, &vertex_lat_values, &num_vertexes);
    edges = delau->get_all_delaunay_edge();

    cv::Mat mat = cv::Mat::zeros(max_lat*10, max_lon*10, CV_8UC3);
    //createAlphaMat(mat);
    for(unsigned int i = 0; edges.size(); i++)
        draw_line<double>(mat, edges[i].v[0], edges[i].v[1]);

    write_to_file(mat, img_path);
};

TEST(DelaunayTriangulationTest, Random) {
    
    int num_points = 360;
    double min_lat = 10.0;
    double max_lat = 80.0;
    double min_lon = 10.0;
    double max_lon = 80.0;
    double *lat, *lon;
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

    /*
    for(unsigned int i = 0; i < points.size(); i++)
        if(points[i].x < 1800 && points[i].y < 1800)
            points_part.push_back(points[i]);

    do_triangulation_plot(points_part, "log/image2-1.png");

    points_part.clear();
    for(unsigned int i = 0; i < points.size(); i++)
        if(points[i].x > 600 && points[i].y > 600)
            points_part.push_back(points[i]);

    do_triangulation_plot(points_part, "log/image2-2.png");

    points_part.clear();
    for(unsigned int i = 0; i < points.size(); i++)
        if(points[i].x > 600 && points[i].y < 1800)
            points_part.push_back(points[i]);

    do_triangulation_plot(points_part, "log/image2-3.png");

    points_part.clear();
    for(unsigned int i = 0; i < points.size(); i++)
        if(points[i].x < 1800 && points[i].y > 600)
            points_part.push_back(points[i]);

    do_triangulation_plot(points_part, "log/image2-4.png");
    */
};
