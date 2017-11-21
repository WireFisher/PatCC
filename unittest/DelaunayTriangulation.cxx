#include "mpi.h"
#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include "../processing_unit_mgt.h"
#include "../delaunay_grid_decomposition_mgt.h"
#include "../delaunay_triangulation_2D.h"

#include <opencv2/opencv.hpp>

extern Grid_info_manager *grid_info_mgr;
extern Process_thread_manager *process_thread_mgr;

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
void draw_line(cv::Mat img, Point<T> start, Point<T> end)
{
    int thickness = 1;
    int line_type = cv::LINE_8;

    cv::line(img, cv::Point(start.x, start.y), cv::Point(end.x, end.y), cv::Scalar(255, 255, 255), thickness, line_type);
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

TEST(DelaunayTriangulationTest, Basic) {
    std::vector<Point<double> > points; 
    
    for(int i = 0; i < 20; i++)
        for(int j = 0; j < 20; j++)
            points.push_back(Point<double>(0 + i * 10, 0 + j * 10));

    Delaunay<double> triangulation;
    triangulation.triangulate(points);
    std::vector<Edge<double> > edges = triangulation.getEdges();

    cv::Mat mat = cv::Mat::zeros(220, 220, CV_8UC3);
    //createAlphaMat(mat);
    for(std::vector<Edge<double> >::iterator e = edges.begin(); e != edges.end(); e++)
        draw_line<double>(mat, e->p1, e->p2);

    write_to_file(mat, "log/image1.png");
};

TEST(DelaunayTriangulationTest, Random) {
    std::vector<Point<double> > points; 
    
    for(int i = 0; i < 20; i++)
        for(int j = 0; j < 20; j++)
            points.push_back(Point<double>(fRand(0.0, 219.0), fRand(0.0, 219.0)));

    Delaunay<double> triangulation;
    triangulation.triangulate(points);
    std::vector<Edge<double> > edges = triangulation.getEdges();

    cv::Mat mat = cv::Mat::zeros(220, 220, CV_8UC3);
    //createAlphaMat(mat);
    for(std::vector<Edge<double> >::iterator e = edges.begin(); e != edges.end(); e++)
        draw_line<double>(mat, e->p1, e->p2);

    write_to_file(mat, "log/image2.png");
};
