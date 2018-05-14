#include "opencv_utils.h"
#include "delaunay_voronoi_2D.h"
#include <opencv2/opencv.hpp>

void plot_edge_into_file(const char *filename, double *head_coord[2], double *tail_coord[2], int num_edges,
                         double min_x, double max_x, double min_y, double max_y)
{
    std::vector<Edge*> edges;

    cv::Mat mat = cv::Mat::zeros(180*10, 400*10, CV_8UC3); /* rows, columns*/

    cv::line(mat, cv::Point((0.0+20.0) * 10, (-90.0+90.0) * 10), cv::Point((0.0+20.0) * 10, (90.0+90.0) * 10), cv::Scalar(0, 255, 255), 1, cv::LINE_8);
    cv::line(mat, cv::Point((0.0+20.0) * 10, (90.0+90.0) * 10), cv::Point((360.0+20.0) * 10, (90.0+90.0) * 10), cv::Scalar(0, 255, 255), 1, cv::LINE_8);
    cv::line(mat, cv::Point((360.0+20.0) * 10, (90.0+90.0) * 10), cv::Point((360.0+20.0) * 10, (-90.0+90.0) * 10), cv::Scalar(0, 255, 255), 1, cv::LINE_8);
    cv::line(mat, cv::Point((360.0+20.0) * 10, (-90.0+90.0) * 10), cv::Point((0.0+20.0) * 10, (-90.0+90.0) * 10), cv::Scalar(0, 255, 255), 1, cv::LINE_8);

    /*
    edges = delau->get_all_legal_delaunay_edge();
    for(unsigned int i = 0; i < edges.size(); i++) {
        draw_line(mat, edges[i], min_lon, max_lon, min_lat, max_lat, cv::Scalar(255, 0, 0));
    }
    */

    if(min_x != 0.0 || max_x != 0.0 || min_y != 0.0 || max_y != 0.0) {
        cv::line(mat, cv::Point((min_x+20.0) * 10, (min_y+90.0) * 10), cv::Point((min_x+20.0) * 10, (max_y+90.0) * 10), cv::Scalar(0, 0, 255), 2, cv::LINE_8);
        cv::line(mat, cv::Point((min_x+20.0) * 10, (max_y+90.0) * 10), cv::Point((max_x+20.0) * 10, (max_y+90.0) * 10), cv::Scalar(0, 0, 255), 2, cv::LINE_8);
        cv::line(mat, cv::Point((max_x+20.0) * 10, (max_y+90.0) * 10), cv::Point((max_x+20.0) * 10, (min_y+90.0) * 10), cv::Scalar(0, 0, 255), 2, cv::LINE_8);
        cv::line(mat, cv::Point((max_x+20.0) * 10, (min_y+90.0) * 10), cv::Point((min_x+20.0) * 10, (min_y+90.0) * 10), cv::Scalar(0, 0, 255), 2, cv::LINE_8);
    }

    for(int i = 0; i < num_edges; i++)
        cv::line(mat, cv::Point((head_coord[0][i]+20.0) * 10, (head_coord[1][i]+90.0) * 10), cv::Point((tail_coord[0][i]+20.0) * 10, (tail_coord[1][i]+90.0) * 10),
                 cv::Scalar(255, 255, 255), 1, cv::LINE_8);
    //double x1 = 270.010468;
    //double y1 = -0.999945;
    //double x2 = 270.550415;
    //double y2 = -42.530453;
    //cv::line(mat, cv::Point((x1+20.0) * 10, (y1+90.0) * 10), cv::Point((x2+20.0) * 10, (y2+90.0) * 10), cv::Scalar(0, 0, 255), 2, cv::LINE_8);
    //cv::line(mat, cv::Point(215.4 * 10, 90.0 * 10), cv::Point(143.6 * 10, 90 * 10), cv::Scalar(0, 255, 0), 2, cv::LINE_8);

    std::vector<int> compression_params;
    compression_params.push_back(cv::IMWRITE_PNG_COMPRESSION);
    compression_params.push_back(9);
    try {
        cv::imwrite(filename, mat, compression_params);
    }
    catch (cv::Exception& ex) {
        fprintf(stderr, "Exception converting image to PNG format: %s\n", ex.what());
    }
}


void plot_projected_edge_into_file(const char *filename, double *head_coord[2], double *tail_coord[2], int num_edges,
                                   int color, int filemode)
{
    std::vector<Edge*> edges;
    int x_shift = 1500;
    int y_shift = 1500;
    int scale = 10;
    
    cv::Mat mat;
    if(filemode == PDLN_PLOT_FILEMODE_NEW)
        mat = cv::Mat::zeros(300*10, 300*10, CV_8UC3); /* rows, columns*/
    else if(filemode == PDLN_PLOT_FILEMODE_APPEND)
        mat = cv::imread(filename, CV_LOAD_IMAGE_COLOR);

    if(!mat.data) {
        perror("PLOT in APPEND MODE but file doesn't exist.\n");
        return;
    }

    cv::line(mat, cv::Point(-150*scale+x_shift,    0*scale+y_shift), cv::Point(150*scale+x_shift,   0*scale+y_shift), cv::Scalar(0, 255, 255), 1, cv::LINE_8);
    cv::line(mat, cv::Point(   0*scale+x_shift, -150*scale+y_shift), cv::Point(  0*scale+x_shift, 150*scale+y_shift), cv::Scalar(0, 255, 255), 1, cv::LINE_8);

    cv::Scalar cv_color;
    switch(color) {
        case PDLN_PLOT_COLOR_WHITE: cv_color = cv::Scalar(255, 255, 255); break;
        case PDLN_PLOT_COLOR_RED  : cv_color = cv::Scalar(0x80, 0x80, 0xF0); break;
        default : cv_color = cv::Scalar(255, 255, 255);
    }

    for(int i = 0; i < num_edges; i++)
        cv::line(mat, cv::Point(head_coord[0][i] * scale + x_shift, head_coord[1][i] * scale + y_shift),
                      cv::Point(tail_coord[0][i] * scale + x_shift, tail_coord[1][i] * scale + y_shift),
                 cv_color, 1, cv::LINE_8);

    std::vector<int> compression_params;
    compression_params.push_back(cv::IMWRITE_PNG_COMPRESSION);
    compression_params.push_back(9);
    try {
        cv::imwrite(filename, mat, compression_params);
    }
    catch (cv::Exception& ex) {
        fprintf(stderr, "Exception converting image to PNG format: %s\n", ex.what());
    }
}


void plot_points_info_file(const char *filename, double *x, double *y, int num, int mode,
                           double min_x, double max_x, double min_y, double max_y)
{
    int x_shift, y_shift, scale;
    cv::Mat mat;

    if(mode == PDLN_PLOT_GLOBAL) {
        x_shift = 200;
        y_shift = 900;
        scale = 10;
        mat = cv::Mat::zeros(180*10, 400*10, CV_8UC3); /* rows, columns*/

        cv::line(mat, cv::Point(  0.0*scale+x_shift, -90*scale+y_shift), cv::Point(  0.0*scale+x_shift,  90*scale+y_shift), cv::Scalar(0, 255, 255), 1, cv::LINE_8);
        cv::line(mat, cv::Point(  0.0*scale+x_shift,  90*scale+y_shift), cv::Point(360.0*scale+x_shift,  90*scale+y_shift), cv::Scalar(0, 255, 255), 1, cv::LINE_8);
        cv::line(mat, cv::Point(360.0*scale+x_shift,  90*scale+y_shift), cv::Point(360.0*scale+x_shift, -90*scale+y_shift), cv::Scalar(0, 255, 255), 1, cv::LINE_8);
        cv::line(mat, cv::Point(360.0*scale+x_shift, -90*scale+y_shift), cv::Point(  0.0*scale+x_shift, -90*scale+y_shift), cv::Scalar(0, 255, 255), 1, cv::LINE_8);
    }
    else if(mode == PDLN_PLOT_PROJECTION) {
        x_shift = 1500;
        y_shift = 1500;
        scale = 500;
        mat = cv::Mat::zeros(300*10, 300*10, CV_8UC3); /* rows, columns*/

        cv::line(mat, cv::Point(-150*scale+x_shift,    0*scale+y_shift), cv::Point(150*scale+x_shift,   0*scale+y_shift), cv::Scalar(0, 255, 255), 1, cv::LINE_8);
        cv::line(mat, cv::Point(   0*scale+x_shift, -150*scale+y_shift), cv::Point(  0*scale+x_shift, 150*scale+y_shift), cv::Scalar(0, 255, 255), 1, cv::LINE_8);
    }

    for(int i = 0; i < num; i++)
        cv::circle(mat, cv::Point(x[i]*scale+x_shift, y[i]*scale+y_shift), 1, cv::Scalar(255, 255, 255), -1, 8);

    if(min_x != 0.0 || max_x != 0.0 || min_y != 0.0 || max_y != 0.0) {
        cv::line(mat, cv::Point(min_x*scale+x_shift, min_y*scale+y_shift), cv::Point(min_x*scale+x_shift, max_y*scale+y_shift), cv::Scalar(0, 0, 255), 2, cv::LINE_8);
        cv::line(mat, cv::Point(min_x*scale+x_shift, max_y*scale+y_shift), cv::Point(max_x*scale+x_shift, max_y*scale+y_shift), cv::Scalar(0, 0, 255), 2, cv::LINE_8);
        cv::line(mat, cv::Point(max_x*scale+x_shift, max_y*scale+y_shift), cv::Point(max_x*scale+x_shift, min_y*scale+y_shift), cv::Scalar(0, 0, 255), 2, cv::LINE_8);
        cv::line(mat, cv::Point(max_x*scale+x_shift, min_y*scale+y_shift), cv::Point(min_x*scale+x_shift, min_y*scale+y_shift), cv::Scalar(0, 0, 255), 2, cv::LINE_8);
    }

    std::vector<int> compression_params;
    compression_params.push_back(cv::IMWRITE_PNG_COMPRESSION);
    compression_params.push_back(9);
    try {
        cv::imwrite(filename, mat, compression_params);
    }
    catch (cv::Exception& ex) {
        fprintf(stderr, "Exception converting image to PNG format: %s\n", ex.what());
    }
}
