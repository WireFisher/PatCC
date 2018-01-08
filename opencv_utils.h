#include <opencv2/opencv.hpp>

void plot_edge_into_file(const char *filename, double *head_coord[2], double *tail_coord[2], int num_edges,
                         double min_x = 0.0, double max_x = 0.0, double min_y = 0.0, double max_y = 0.0)
{
    std::vector<Edge*> edges;

    cv::Mat mat = cv::Mat::zeros(180*10, 360*10, CV_8UC3); /* rows, columns*/

    for(int i = 0; i < num_edges; i++)
        cv::line(mat, cv::Point(head_coord[0][i] * 10, (head_coord[1][i]+90.0) * 10), cv::Point(tail_coord[0][i] * 10, (tail_coord[1][i]+90.0) * 10),
                 cv::Scalar(255, 255, 255), 1, cv::LINE_8);

    /*
    edges = delau->get_all_legal_delaunay_edge();
    for(unsigned int i = 0; i < edges.size(); i++) {
        draw_line(mat, edges[i], min_lon, max_lon, min_lat, max_lat, cv::Scalar(255, 0, 0));
    }
    */

    if(min_x != 0.0 || max_x != 0.0 || min_y != 0.0 || max_y != 0.0) {
        cv::line(mat, cv::Point(min_x * 10, (min_y+90.0) * 10), cv::Point(min_x * 10, (max_y+90.0) * 10), cv::Scalar(0, 0, 255), 2, cv::LINE_8);
        cv::line(mat, cv::Point(min_x * 10, (max_y+90.0) * 10), cv::Point(max_x * 10, (max_y+90.0) * 10), cv::Scalar(0, 0, 255), 2, cv::LINE_8);
        cv::line(mat, cv::Point(max_x * 10, (max_y+90.0) * 10), cv::Point(max_x * 10, (min_y+90.0) * 10), cv::Scalar(0, 0, 255), 2, cv::LINE_8);
        cv::line(mat, cv::Point(max_x * 10, (min_y+90.0) * 10), cv::Point(min_x * 10, (min_y+90.0) * 10), cv::Scalar(0, 0, 255), 2, cv::LINE_8);
    }
    //cv::line(mat, cv::Point(0.0 * 10, 180.0 * 10), cv::Point(360.0 * 10, 180 * 10), cv::Scalar(0, 0, 255), 2, cv::LINE_8);
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
    fprintf(stdout, "Saved PNG file: %s\n", filename);
}
