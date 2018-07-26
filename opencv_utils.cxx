#include "opencv_utils.h"
#include "delaunay_voronoi_2D.h"
#include <opencv2/opencv.hpp>

void plot_edge_into_file(const char *filename, double *head_coord[2], double *tail_coord[2], int num_edges,
                         double min_x, double max_x, double min_y, double max_y)
{
    int x_shift = 200;
    int y_shift = 900;
    int scale = 10;
    std::vector<Edge*> edges;

    cv::Mat mat = cv::Mat::zeros(180*10, 400*10, CV_8UC3); /* rows, columns*/

    cv::line(mat, cv::Point(  0.0*scale+x_shift, -90.0*scale+y_shift), cv::Point(  0.0*scale+x_shift,  90.0*scale+y_shift), cv::Scalar(0, 255, 255), 1, cv::LINE_8);
    cv::line(mat, cv::Point(  0.0*scale+x_shift,  90.0*scale+y_shift), cv::Point(360.0*scale+x_shift,  90.0*scale+y_shift), cv::Scalar(0, 255, 255), 1, cv::LINE_8);
    cv::line(mat, cv::Point(360.0*scale+x_shift,  90.0*scale+y_shift), cv::Point(360.0*scale+x_shift, -90.0*scale+y_shift), cv::Scalar(0, 255, 255), 1, cv::LINE_8);
    cv::line(mat, cv::Point(360.0*scale+x_shift, -90.0*scale+y_shift), cv::Point(  0.0*scale+x_shift, -90.0*scale+y_shift), cv::Scalar(0, 255, 255), 1, cv::LINE_8);

    if(min_x != 0.0 || max_x != 0.0 || min_y != 0.0 || max_y != 0.0) {
        cv::line(mat, cv::Point(min_x*scale+x_shift, min_y*scale+y_shift), cv::Point(min_x*scale+x_shift, max_y*scale+y_shift), cv::Scalar(0, 0, 255), 2, cv::LINE_8);
        cv::line(mat, cv::Point(min_x*scale+x_shift, max_y*scale+y_shift), cv::Point(max_x*scale+x_shift, max_y*scale+y_shift), cv::Scalar(0, 0, 255), 2, cv::LINE_8);
        cv::line(mat, cv::Point(max_x*scale+x_shift, max_y*scale+y_shift), cv::Point(max_x*scale+x_shift, min_y*scale+y_shift), cv::Scalar(0, 0, 255), 2, cv::LINE_8);
        cv::line(mat, cv::Point(max_x*scale+x_shift, min_y*scale+y_shift), cv::Point(min_x*scale+x_shift, min_y*scale+y_shift), cv::Scalar(0, 0, 255), 2, cv::LINE_8);
    }

    for(int i = 0; i < num_edges; i++)
        cv::line(mat, cv::Point(head_coord[0][i]*scale+x_shift, head_coord[1][i]*scale+y_shift), cv::Point(tail_coord[0][i]*scale+x_shift, tail_coord[1][i]*scale+y_shift),
                 cv::Scalar(255, 255, 255), 1, cv::LINE_8);

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


void plot_points_into_file(const char *filename, double *x, double *y, int num, int mode,
                           double min_x, double max_x, double min_y, double max_y)
{
    int x_shift, y_shift, scale;
    cv::Mat mat;

    if(mode == PDLN_PLOT_GLOBAL) {
        x_shift = 200;
        y_shift = 900;
        scale = 10;
        mat = cv::Mat(180*10, 400*10, CV_8UC3, cv::Scalar(255, 255, 255)); /* rows, columns*/

        cv::rectangle(mat, cv::Point(  0.0*scale+x_shift, -90*scale+y_shift), cv::Point(360.0*scale+x_shift,  90*scale+y_shift), cv::Scalar(0, 255, 255));
    }
    else if(mode == PDLN_PLOT_PROJECTION) {
        x_shift = 1500;
        y_shift = 1500;
        scale = 500;
        mat = cv::Mat(300*10, 300*10, CV_8UC3, cv::Scalar(255, 255, 255)); /* rows, columns*/

        cv::line(mat, cv::Point(-150*scale+x_shift,    0*scale+y_shift), cv::Point(150*scale+x_shift,   0*scale+y_shift), cv::Scalar(0, 255, 255), 1, cv::LINE_8);
        cv::line(mat, cv::Point(   0*scale+x_shift, -150*scale+y_shift), cv::Point(  0*scale+x_shift, 150*scale+y_shift), cv::Scalar(0, 255, 255), 1, cv::LINE_8);
    }

    for(int i = 0; i < num; i++)
        cv::circle(mat, cv::Point(x[i]*scale+x_shift, y[i]*scale+y_shift), 2, cv::Scalar(0, 0, 0), -1);

    if(min_x != 0.0 || max_x != 0.0 || min_y != 0.0 || max_y != 0.0)
        cv::rectangle(mat, cv::Point(min_x*scale+x_shift, min_y*scale+y_shift), cv::Point(max_x*scale+x_shift, max_y*scale+y_shift), cv::Scalar(0, 0, 255), 2);

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

    cv::Scalar cv_color;
    switch(color) {
        case PDLN_PLOT_COLOR_RED  : cv_color = cv::Scalar(0x80, 0x80, 0xF0); break;
        case PDLN_PLOT_COLOR_WHITE: cv_color = cv::Scalar(255, 255, 255); break;
        default : cv_color = cv::Scalar(255, 255, 255);
    }

    cv::rectangle(mat, cv::Point(x_min*scale+x_shift, y_min*scale+y_shift), cv::Point(x_max*scale+x_shift, y_max*scale+y_shift), cv_color, 5);

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
    int scale = 10;
    x_min = x_min*scale+x_shift;
    x_max = x_max*scale+x_shift;
    y_min = y_min*scale+y_shift;
    y_max = y_max*scale+y_shift;

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
    double font_scale = 1;
    int baseline;
    cv::Point origin;

    cv::Size text_size = cv::getTextSize(text, font_face, font_scale, thickness, &baseline);
    font_scale = std::floor(std::min((x_max-x_min)/text_size.width, (y_max-y_min)/text_size.height));
    font_scale /= 3.;
    text_size = cv::getTextSize(text, font_face, font_scale, thickness, &baseline);
    origin.x = x_min + (x_max-x_min)/2. - text_size.width/2.;
    origin.y = y_min + (y_max-y_min)/2. + text_size.height/2.;

    cv::putText(mat, text, origin, font_face, font_scale, cv_color, thickness);

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
