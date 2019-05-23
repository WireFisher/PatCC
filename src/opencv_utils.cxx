/***************************************************************
  *  Copyright (c) 2019, Tsinghua University.
  *  This is a source file of PatCC.
  *  This file was initially finished by Dr. Li Liu and
  *  Haoyu Yang. If you have any problem,
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#include "opencv_utils.h"
#include <opencv2/opencv.hpp>

void plot_edge_into_file(const char *filename, double *head_coord[2], double *tail_coord[2], int num_edges,
                         double min_x, double max_x, double min_y, double max_y)
{
    int x_shift = 200;
    int y_shift = 900;
    int x_scale = 10;
    int y_scale = -10;

    cv::Mat mat = cv::Mat::zeros(180*10, 400*10, CV_8UC3); /* rows, columns*/

    cv::rectangle(mat, cv::Point(0*x_scale+x_shift, -90*y_scale+y_shift), cv::Point(360*x_scale+x_shift, 90*y_scale+y_shift), cv::Scalar(0, 255, 255), 1, 8);

    if(min_x != 0.0 || max_x != 0.0 || min_y != 0.0 || max_y != 0.0)
        cv::rectangle(mat, cv::Point(min_x*x_scale+x_shift, min_y*y_scale+y_shift), cv::Point(max_x*x_scale+x_shift, max_y*y_scale+y_shift), cv::Scalar(0, 0, 255), 2, 8);

    for(int i = 0; i < num_edges; i++)
        cv::line(mat, cv::Point(head_coord[0][i]*x_scale+x_shift, head_coord[1][i]*y_scale+y_shift), cv::Point(tail_coord[0][i]*x_scale+x_shift, tail_coord[1][i]*y_scale+y_shift),
                 cv::Scalar(255, 255, 255), 1, 8);

    std::vector<int> compression_params;
    compression_params.push_back(cv::IMWRITE_PNG_COMPRESSION);
    compression_params.push_back(9);
    try {
        cv::imwrite(filename, mat, compression_params);
    } catch (cv::Exception& ex) {
        fprintf(stderr, "Exception converting image to PNG format: %s\n", ex.what());
    }
}


void plot_projected_edge_into_file(const char *filename, double *head_coord[2], double *tail_coord[2], int num_edges,
                                   int color, int filemode)
{
    int x_shift = 1500;
    int y_shift = 1500;
    int x_scale = 10;
    int y_scale = -10;
    
    cv::Mat mat;
    if(filemode == PDLN_PLOT_FILEMODE_NEW)
        mat = cv::Mat::zeros(300*10, 300*10, CV_8UC3); /* rows, columns*/
    else if(filemode == PDLN_PLOT_FILEMODE_APPEND)
        mat = cv::imread(filename, CV_LOAD_IMAGE_COLOR);

    if(!mat.data) {
        perror("PLOT in APPEND MODE but file doesn't exist.\n");
        return;
    }

    cv::line(mat, cv::Point(-150*x_scale+x_shift,    0*y_scale+y_shift), cv::Point(150*x_scale+x_shift,   0*y_scale+y_shift), cv::Scalar(0, 255, 255), 1, 8);
    cv::line(mat, cv::Point(   0*x_scale+x_shift, -150*y_scale+y_shift), cv::Point(  0*x_scale+x_shift, 150*y_scale+y_shift), cv::Scalar(0, 255, 255), 1, 8);

    cv::Scalar cv_color;
    switch(color) {
        case PDLN_PLOT_COLOR_WHITE: cv_color = cv::Scalar(255, 255, 255); break;
        case PDLN_PLOT_COLOR_RED  : cv_color = cv::Scalar(0x80, 0x80, 0xF0); break;
        default : cv_color = cv::Scalar(255, 255, 255);
    }

    for(int i = 0; i < num_edges; i++)
        cv::line(mat, cv::Point(head_coord[0][i] * x_scale + x_shift, head_coord[1][i] * y_scale + y_shift),
                      cv::Point(tail_coord[0][i] * x_scale + x_shift, tail_coord[1][i] * y_scale + y_shift),
                 cv_color, 1, 8);

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


void plot_points_into_file(const char *filename, double *x, double *y, bool* mask, int num, int mode,
                           double min_x, double max_x, double min_y, double max_y)
{
    int x_shift, y_shift, x_scale, y_scale;
    cv::Mat mat;

    if(mode == PDLN_PLOT_GLOBAL) {
        x_shift = 200;
        y_shift = 900;
        x_scale = 10;
        y_scale = -10;
        mat = cv::Mat(180*10, 400*10, CV_8UC3, cv::Scalar(255, 255, 255)); /* rows, columns*/

        cv::rectangle(mat, cv::Point(  0.0*x_scale+x_shift, -90*y_scale+y_shift), cv::Point(360.0*x_scale+x_shift,  90*y_scale+y_shift), cv::Scalar(0, 255, 255));
    }
    else if(mode == PDLN_PLOT_PROJECTION) {
        x_shift = 1500;
        y_shift = 1500;
        x_scale = 500;
        y_scale = -500;
        mat = cv::Mat(300*10, 300*10, CV_8UC3, cv::Scalar(255, 255, 255)); /* rows, columns*/

        cv::line(mat, cv::Point(-150*x_scale+x_shift,    0*y_scale+y_shift), cv::Point(150*x_scale+x_shift,   0*y_scale+y_shift), cv::Scalar(0, 255, 255), 1, 8);
        cv::line(mat, cv::Point(   0*x_scale+x_shift, -150*y_scale+y_shift), cv::Point(  0*x_scale+x_shift, 150*y_scale+y_shift), cv::Scalar(0, 255, 255), 1, 8);
    }

    for(int i = 0; i < num; i++) {
        if (!mask || mask[i])
            cv::circle(mat, cv::Point(x[i]*x_scale+x_shift, y[i]*y_scale+y_shift), 2, cv::Scalar(0, 0, 0), -1);
    }

    if(min_x != 0.0 || max_x != 0.0 || min_y != 0.0 || max_y != 0.0)
        cv::rectangle(mat, cv::Point(min_x*x_scale+x_shift, min_y*y_scale+y_shift), cv::Point(max_x*x_scale+x_shift, max_y*y_scale+y_shift), cv::Scalar(0, 0, 255), 2);

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


void plot_rectangle_into_file(const char *filename, double x_min, double x_max, double y_min, double y_max, int color, int filemode)
{
    int x_shift = 200;
    int y_shift = 900;
    int x_scale = 10;
    int y_scale = -10;
    
    cv::Mat mat;
    if(filemode == PDLN_PLOT_FILEMODE_NEW)
        mat = cv::Mat::zeros(300*10, 300*10, CV_8UC3); /* rows, columns*/
    else if(filemode == PDLN_PLOT_FILEMODE_APPEND)
        mat = cv::imread(filename, CV_LOAD_IMAGE_COLOR);

    if(!mat.data) {
        perror("PLOT in APPEND MODE but file doesn't exist.\n");
        return;
    }

    cv::Scalar cv_color;
    switch(color) {
        case PDLN_PLOT_COLOR_RED  : cv_color = cv::Scalar(0x80, 0x80, 0xF0); break;
        case PDLN_PLOT_COLOR_WHITE: cv_color = cv::Scalar(255, 255, 255); break;
        default : cv_color = cv::Scalar(255, 255, 255);
    }

    cv::rectangle(mat, cv::Point(x_min*x_scale+x_shift, y_min*y_scale+y_shift), cv::Point(x_max*x_scale+x_shift, y_max*y_scale+y_shift), cv_color, 5);

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

void plot_text_into_file(const char *filename, const char *text, double x_min, double x_max, double y_min, double y_max, int color)
{
    int x_shift = 200;
    int y_shift = 900;
    int x_scale = 10;
    int y_scale = -10;

    x_min = x_min*x_scale+x_shift;
    x_max = x_max*x_scale+x_shift;
    y_min = y_min*y_scale+y_shift;
    y_max = y_max*y_scale+y_shift;

    cv::Mat mat = cv::imread(filename, CV_LOAD_IMAGE_COLOR);
    if(!mat.data) {
        perror("PLOT in APPEND MODE but file doesn't exist.\n");
        return;
    }

    cv::Scalar cv_color;
    switch(color) {
        case PDLN_PLOT_COLOR_RED  : cv_color = cv::Scalar(0x80, 0x80, 0xF0); break;
        case PDLN_PLOT_COLOR_WHITE: cv_color = cv::Scalar(255, 255, 255); break;
        default : cv_color = cv::Scalar(255, 255, 255);
    }

    const int font_face = cv::FONT_HERSHEY_SIMPLEX;
    const int thickness = 3;
    double font_size = 1;
    int baseline;
    cv::Point origin;

    cv::Size text_size = cv::getTextSize(text, font_face, font_size, thickness, &baseline);
    font_size = std::floor(std::min((x_max-x_min)/text_size.width, (y_max-y_min)/text_size.height));
    font_size /= 3.;
    text_size = cv::getTextSize(text, font_face, font_size, thickness, &baseline);
    origin.x = x_min + (x_max-x_min)/2. - text_size.width/2.;
    origin.y = y_min + (y_max-y_min)/2. + text_size.height/2.;

    cv::putText(mat, text, origin, font_face, font_size, cv_color, thickness);

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
