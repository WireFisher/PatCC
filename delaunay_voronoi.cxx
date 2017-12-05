#include "delaunay_voronoi.h"
#include <cassert>
#include <cmath>
#include <cstdio>
#include <sys/time.h>


#define FLOAT_ERROR 1e-10

Delaunay_Voronoi *current_delaunay_voronoi = NULL; //thread safe?

/*
 *            o Center
 *            ^
 *           /
 *          /
 *         /
 *     V1 o -----> o V2
 */
double compute_three_2D_points_cross_product(double center_coord1_value,
                                             double center_coord2_value,
                                             double vertex1_coord1_value,
                                             double vertex1_coord2_value,
                                             double vertex2_coord1_value,
                                             double vertex2_coord2_value)
{
    double vertexes_coord1_value_diff, vertexes_coord2_value_diff;
    double center_vertex_coord1_value_diff, center_vertex_coord2_value_diff;
    
    
    vertexes_coord1_value_diff = vertex2_coord1_value - vertex1_coord1_value;
    vertexes_coord2_value_diff = vertex2_coord2_value - vertex1_coord2_value;
    center_vertex_coord1_value_diff = center_coord1_value - vertex1_coord1_value;
    center_vertex_coord2_value_diff = center_coord2_value - vertex1_coord2_value;
    

    return vertexes_coord1_value_diff*center_vertex_coord2_value_diff - center_vertex_coord1_value_diff*vertexes_coord2_value_diff;
}

inline double det(const Point *pt1, const Point *pt2, const Point *pt3)
{

//    return (pt1.lon-pt3.lon)*(pt2.lat-pt3.lat) - (pt1.lat-pt3.lat)*(pt2.lon-pt3.lon);

//    return compute_three_3D_points_cross_product(pt3->x, pt3->y, pt3->z, pt1->x, pt1->y, pt1->z, pt2->x, pt2->y, pt2->z);
    return compute_three_2D_points_cross_product(pt3->x, pt3->y, pt1->x, pt1->y, pt2->x, pt2->y);
}

Point::Point()
{
    this->current_triangle = NULL;
}
Point::Point(double x, double y, int id)
{
    this->x = x;
    this->y = y;
    this->id = id;
    this->current_triangle = NULL;
}


double Point::calculate_distance(const Point *pt) const
{
    double dx = pt->x - x;
    double dy = pt->y - y;
    return sqrt(dx * dx + dy * dy);
}

Point::~Point()
{
}

Point operator-(Point p1, Point p2)
{
    return Point(p1.x - p2.x, p1.y - p2.y, -1);
}

/**
 * Check point's position relative to an edge<pt1, pt2>
 * Points should be distinct
 * @param  pt1    the head of the edge
 * @param  pt2    the head of the edge 
 * @return    1    left
 *            0    on
 *            -1    right
 */
int Point::position_to_edge(const Point *pt1, const Point *pt2) const
{
    double res1 = det(pt1, pt2, this);

    if (std::fabs(res1) <= FLOAT_ERROR)
        return 0;
    if (res1 > 0)
        return 1;
    return -1;
}

/**
 * Check point's position relative to a triangle
 * This point and points of the triangle should be distinct
 * @return     0    inside
 *            -1    outside
 *             1    lies on the edge <pt1, pt2>
 *             2    lies on the edge <pt2, pt3>
 *             3    lies on the edge <pt3, pt1>
 */
int Point::position_to_triangle(const Triangle *triangle) const
{
    int pos, ret = 0;
    pos = position_to_edge(triangle->v[0], triangle->v[1]);
    if (pos == -1)
        return -1;
    else if (pos == 0)
        ret = 1;
    pos = position_to_edge(triangle->v[1], triangle->v[2]);
    if (pos == -1)
        return -1;
    else if (pos == 0)
        ret = 2;
    pos = position_to_edge(triangle->v[2], triangle->v[0]);
    if (pos == -1)
        return -1;
    else if (pos == 0)
        ret = 3;
    return ret;
}


Edge::Edge(Point *head, Point *tail) 
{
    this->head = head;
    this->tail = tail;
    twin_edge = NULL; 
    prev_edge_in_triangle = NULL; 
    next_edge_in_triangle = NULL;
    triangle = NULL;
}

Edge *Edge::generate_twins_edge()
{
    Edge *twins_edge = current_delaunay_voronoi->allocate_edge(tail, head);
    twins_edge->twin_edge = this;
    this->twin_edge = twins_edge;

    return twins_edge;
}

Edge::~Edge()
{
}

bool Delaunay_Voronoi::is_triangle_legal(const Point *pt, const Edge *edge)
{
    if (!edge->twin_edge)
        return true;

    return !edge->triangle->circum_circle_contains(edge->twin_edge->prev_edge_in_triangle->head);
    /*
    const Point *vi = edge->head;
    const Point *vj = edge->next_edge_in_triangle->head;
    const Point *vk = edge->twin_edge->prev_edge_in_triangle->head;

    Point temp_point1, temp_point2, temp_point3;

    temp_point1 = *vi - *pt;
    temp_point2 = *vj - *pt;
    temp_point3 = *vk - *pt;

    if (det(&temp_point1, &temp_point2, &temp_point3) >= FLOAT_ERROR)
        return false;
    else return true;
    */
}

void Delaunay_Voronoi::legalize_triangles(Point *vr, Edge *edge, vector<Triangle*> *leaf_triangles)
{
    if (is_triangle_legal(vr, edge))
        return;

    //EXECUTION_REPORT(REPORT_ERROR, -1, edge->triangle->is_leaf, "remap software error1 in legalize_triangles\n");
    assert(edge->triangle->is_leaf);
    //EXECUTION_REPORT(REPORT_ERROR, -1, edge->twin_edge->triangle->is_leaf, "remap software error2 in legalize_triangles %lx\n", (long)(edge->twin_edge->triangle));
    assert(edge->twin_edge->triangle->is_leaf);
    leaf_triangles->push_back(edge->twin_edge->triangle);
    edge->twin_edge->triangle->reference_count ++;
    edge->triangle->is_leaf = false;
    edge->twin_edge->triangle->is_leaf = false;

    Point *vk = edge->twin_edge->prev_edge_in_triangle->head;
    Edge *eij = edge;
    Edge *ejr = eij->next_edge_in_triangle;
    Edge *eri = ejr->next_edge_in_triangle;
    Edge *eji = eij->twin_edge;
    Edge *eik = eji->next_edge_in_triangle;
    Edge *ekj = eik->next_edge_in_triangle;
    Edge *erk = current_delaunay_voronoi->allocate_edge(vr, vk);
    Edge *ekr = erk->generate_twins_edge();
    Triangle* tikr = current_delaunay_voronoi->allocate_Triangle(eik,ekr,eri);
    Triangle* tjrk = current_delaunay_voronoi->allocate_Triangle(ejr,erk,ekj);
    leaf_triangles->push_back(tikr);
    leaf_triangles->push_back(tjrk);
    tikr->reference_count ++;
    tjrk->reference_count ++;
    legalize_triangles(vr, eik, leaf_triangles);
    legalize_triangles(vr, ekj, leaf_triangles);
}

Triangle::Triangle()
{
    edge[0] = NULL;
    edge[1] = NULL;
    edge[2] = NULL;
}

Triangle::Triangle(Point *point1, Point *point2, Point *point3)
{
    if (point1->position_to_edge(point2, point3) == 1)
        initialize_triangle_with_edges(current_delaunay_voronoi->allocate_edge(point1,point2),
                                       current_delaunay_voronoi->allocate_edge(point2,point3),
                                       current_delaunay_voronoi->allocate_edge(point3,point1));
    else initialize_triangle_with_edges(current_delaunay_voronoi->allocate_edge(point3,point2),
                                        current_delaunay_voronoi->allocate_edge(point2,point1),
                                        current_delaunay_voronoi->allocate_edge(point1,point3));

    calulate_circum_circle();
}

Triangle::Triangle(Edge *edge1, Edge *edge2, Edge *edge3)
{
    initialize_triangle_with_edges(edge1, edge2, edge3);
    calulate_circum_circle();
}

Triangle::~Triangle()
{
}

void Triangle::calulate_circum_circle()
{
    double ab, cd, ef;

    ab = (v[0]->x * v[0]->x) + (v[0]->y * v[0]->y);
    cd = (v[1]->x * v[1]->x) + (v[1]->y * v[1]->y);
    ef = (v[2]->x * v[2]->x) + (v[2]->y * v[2]->y);

    circum_center[0] = (ab * (v[2]->y - v[1]->y) + cd * (v[0]->y - v[2]->y) + ef * (v[1]->y - v[0]->y)) /
                       (v[0]->x * (v[2]->y - v[1]->y) + v[1]->x * (v[0]->y - v[2]->y) + v[2]->x * (v[1]->y - v[0]->y)) / 2.f;
    circum_center[1] = (ab * (v[2]->x - v[1]->x) + cd * (v[0]->x - v[2]->x) + ef * (v[1]->x - v[0]->x)) /
                       (v[0]->y * (v[2]->x - v[1]->x) + v[1]->y * (v[0]->x - v[2]->x) + v[2]->y * (v[1]->x - v[0]->x)) / 2.f;
    circum_radius = sqrtf(((v[0]->x - circum_center[0]) * (v[0]->x - circum_center[0])) + ((v[0]->y - circum_center[1]) * (v[0]->y - circum_center[1])));
}

void Triangle::initialize_triangle_with_edges(Edge *edge1, Edge *edge2, Edge *edge3)
{
    Point *pt1, *pt2, *pt3;


    is_leaf = true;
    reference_count = 1;
    
    pt1 = edge1->head;
    pt2 = edge2->head;
    pt3 = edge3->head;

    //EXECUTION_REPORT(REPORT_ERROR, -1, fabs(det(pt1, pt2, pt3)) > e && fabs(det(pt2, pt3, pt1)) > e && fabs(det(pt3, pt1, pt2)) > e,
    //                 "points given to construct triangle are on the same line.");
    assert(std::fabs(det(pt1, pt2, pt3)) > FLOAT_ERROR && std::fabs(det(pt2, pt3, pt1)) > FLOAT_ERROR && std::fabs(det(pt3, pt1, pt2)) > FLOAT_ERROR);
    //EXECUTION_REPORT(REPORT_ERROR, -1, edge1->tail==edge2->head && edge2->tail==edge3->head && edge3->tail==edge1->head, "edges given to construct triangle is invalid.");
    assert(edge1->tail==edge2->head && edge2->tail==edge3->head && edge3->tail==edge1->head);
       
    v[0] = pt1;
    if (pt1->position_to_edge(pt2, pt3) == 1) {
        v[1] = pt2;
        v[2] = pt3;
        this->edge[0] = edge1;
        this->edge[1] = edge2;
        this->edge[2] = edge3;
    }
    else {
        v[1] = pt3;
        v[2] = pt2;
        this->edge[0] = edge3->twin_edge;
        this->edge[1] = edge2->twin_edge;
        this->edge[2] = edge1->twin_edge;
        //EXECUTION_REPORT(REPORT_ERROR, -1, edge3->twin_edge != NULL && edge2->twin_edge != NULL && edge1->twin_edge != NULL, "remap software error3 in new Triangle");
        assert(edge3->twin_edge != NULL && edge2->twin_edge != NULL && edge1->twin_edge != NULL);
    }
    
    this->edge[0]->next_edge_in_triangle = this->edge[1];
    this->edge[1]->next_edge_in_triangle = this->edge[2];
    this->edge[2]->next_edge_in_triangle = this->edge[0];
    this->edge[0]->prev_edge_in_triangle = this->edge[2];
    this->edge[1]->prev_edge_in_triangle = this->edge[0];
    this->edge[2]->prev_edge_in_triangle = this->edge[1];

    this->edge[0]->triangle = this;
    this->edge[1]->triangle = this;
    this->edge[2]->triangle = this;
}

bool Triangle::circum_circle_contains(Point *p)
{
    return sqrt(((p->x - circum_center[0]) * (p->x - circum_center[0])) + ((p->y - circum_center[1]) * (p->y - circum_center[1]))) <= circum_radius;
}

int Triangle::find_best_candidate_point()
{
    double max_min_dist = -1, min_dist, dist;
    int best_candidate_id;

    
    if (remained_points_in_triangle.size() == 0)
        return -1;

    for (int i = 0; i < remained_points_in_triangle.size(); i ++) {
        min_dist = remained_points_in_triangle[i]->calculate_distance(v[0]);
        dist = remained_points_in_triangle[i]->calculate_distance(v[1]);
        if (min_dist > dist)
            min_dist = dist;
        dist = remained_points_in_triangle[i]->calculate_distance(v[2]);
        if (min_dist > dist)
            min_dist = dist;
        if (max_min_dist < min_dist) {
            max_min_dist = min_dist;
            best_candidate_id = i;
        }
    }

    return best_candidate_id;
}


bool Triangle::contain_vertex(Point* pt)
{
    for(int i = 0; i < 3; i++)
        if(v[i] == pt)
            return true;
    return false;
}


void Delaunay_Voronoi::distribute_points_into_triangles(vector<Point*> *pnts, vector<Triangle*> *triangles)
{
    bool find_triangle;


    for (int i = 0; i < pnts->size(); i ++) {
        find_triangle = false;
        for (int j = 0; j < triangles->size(); j ++) {
            if (!((*triangles)[j])->is_leaf)
                continue;
            if ((*pnts)[i]->position_to_triangle(((*triangles)[j])) >= 0) {
                (*pnts)[i]->current_triangle = (*triangles)[j];
                (*triangles)[j]->remained_points_in_triangle.push_back((*pnts)[i]);
                find_triangle = true;
                break;
            }
        }
        if (!find_triangle) 
            if (is_global_grid)
                //EXECUTION_REPORT(REPORT_ERROR, -1, false, "CoR may have bugs, please contact liuli-cess@tsinghua.edu.cn");
                assert(false);
            else //EXECUTION_REPORT(REPORT_ERROR, -1, false, "please enlarge the boundary of the regional grid"); 
                assert(false);
    }
}

void Delaunay_Voronoi::triangularization_process(Triangle *triangle)
{
    int best_candidate_point_id;
    Point *best_candidate_point;
    vector<Triangle *> leaf_triangles;
    vector<Triangle *> existing_triangles;


    if (!triangle->is_leaf) {
        triangle->reference_count --;
        return;
    }
        
    if (triangle->remained_points_in_triangle.size() == 0) {
        result_leaf_triangles.push_back(triangle);
        return;
    }

    triangle->is_leaf = false;
    
    best_candidate_point_id = triangle->find_best_candidate_point();
    best_candidate_point = triangle->remained_points_in_triangle[best_candidate_point_id];
    triangle->remained_points_in_triangle.erase(triangle->remained_points_in_triangle.begin()+best_candidate_point_id);

    if (best_candidate_point->position_to_triangle(triangle) == 0) { //inside
        Edge *e_v1_can = current_delaunay_voronoi->allocate_edge(triangle->v[0], best_candidate_point);
        Edge *e_can_v1 = e_v1_can->generate_twins_edge();
        Edge *e_v2_can = current_delaunay_voronoi->allocate_edge(triangle->v[1], best_candidate_point);
        Edge *e_can_v2 = e_v2_can->generate_twins_edge();
        Edge *e_v3_can = current_delaunay_voronoi->allocate_edge(triangle->v[2], best_candidate_point);
        Edge *e_can_v3 = e_v3_can->generate_twins_edge();
        Triangle *t_v1_v2_can = current_delaunay_voronoi->allocate_Triangle(triangle->edge[0], e_v2_can, e_can_v1);
        Triangle *t_v2_v3_can = current_delaunay_voronoi->allocate_Triangle(triangle->edge[1], e_v3_can, e_can_v2);
        Triangle *t_v3_v1_can = current_delaunay_voronoi->allocate_Triangle(triangle->edge[2], e_v1_can, e_can_v3);
        leaf_triangles.push_back(triangle);
        leaf_triangles.push_back(t_v1_v2_can);
        leaf_triangles.push_back(t_v2_v3_can);
        leaf_triangles.push_back(t_v3_v1_can);
        triangle->reference_count ++;
        t_v1_v2_can->reference_count ++;
        t_v2_v3_can->reference_count ++;
        t_v3_v1_can->reference_count ++;
        legalize_triangles(best_candidate_point, triangle->edge[0], &leaf_triangles);
        legalize_triangles(best_candidate_point, triangle->edge[1], &leaf_triangles);
        legalize_triangles(best_candidate_point, triangle->edge[2], &leaf_triangles);
    }
    else {
        Point *vi, *vj, *vk, *vl;
        Edge *eij, *ejk, *eki;
        Edge *eil, *elj, *eji;
        switch (best_candidate_point->position_to_triangle(triangle)) {
            case 1:
                vi = triangle->v[0];
                vj = triangle->v[1];
                vk = triangle->v[2];
                eij = triangle->edge[0];
                break;
            case 2:
                vi = triangle->v[1];
                vj = triangle->v[2];
                vk = triangle->v[0];
                eij = triangle->edge[1];
                break;
            case 3:
                vi = triangle->v[2];
                vj = triangle->v[0];
                vk = triangle->v[1];
                eij = triangle->edge[2];
                break;
            default:
                //EXECUTION_REPORT(REPORT_ERROR, -1, false, "remap software error2 in triangularization_process");
                assert(false);
                break;
        }
        //EXECUTION_REPORT(REPORT_ERROR, -1, best_candidate_point->position_to_edge(vi, vj) == 0, "remap software error3 in triangularization_process");
        assert(best_candidate_point->position_to_edge(vi, vj) == 0);
        if (eij->twin_edge != NULL)
            //EXECUTION_REPORT(REPORT_ERROR, -1, eij->twin_edge->triangle->is_leaf, "remap software error3 in triangularization_process");
            assert(eij->twin_edge->triangle->is_leaf);
        ejk = eij->next_edge_in_triangle;
        eki = ejk->next_edge_in_triangle;
        if (eij->twin_edge != NULL) { 
            eji = eij->twin_edge;
            eij->twin_edge->triangle->is_leaf = false;
            eil = eji->next_edge_in_triangle;
            elj = eil->next_edge_in_triangle;
            vl = elj->head;
        }
        Edge *eri = current_delaunay_voronoi->allocate_edge(best_candidate_point, vi);
        Edge *eir = eri->generate_twins_edge();
        Edge *erk = current_delaunay_voronoi->allocate_edge(best_candidate_point, vk);
        Edge *ekr = erk->generate_twins_edge();
        Edge *erj = current_delaunay_voronoi->allocate_edge(best_candidate_point, vj);
        Triangle* tirk = current_delaunay_voronoi->allocate_Triangle(eir, erk, eki);
        Triangle* tjkr = current_delaunay_voronoi->allocate_Triangle(ejk, ekr, erj);
        leaf_triangles.push_back(triangle);
        leaf_triangles.push_back(tirk);
        leaf_triangles.push_back(tjkr);
        triangle->reference_count ++;
        tirk->reference_count ++;
        tjkr->reference_count ++;
        legalize_triangles(best_candidate_point, ejk, &leaf_triangles);
        legalize_triangles(best_candidate_point, eki, &leaf_triangles); 
        if (eij->twin_edge != NULL) {
            Edge *ejr = erj->generate_twins_edge();
            Edge *erl = current_delaunay_voronoi->allocate_edge(best_candidate_point, vl);
            Edge *elr = erl->generate_twins_edge();
            Triangle* tilr = current_delaunay_voronoi->allocate_Triangle(eil, elr, eri);
            Triangle* tjrl = current_delaunay_voronoi->allocate_Triangle(ejr, erl, elj);
            leaf_triangles.push_back(eij->twin_edge->triangle);
            leaf_triangles.push_back(tilr);
            leaf_triangles.push_back(tjrl);
            tilr->reference_count ++;
            tjrl->reference_count ++;
            legalize_triangles(best_candidate_point, eil, &leaf_triangles);
            legalize_triangles(best_candidate_point, elj, &leaf_triangles);
        }
        else {
            eir->twin_edge = NULL;
        }
    }

    for (int i = 0; i < leaf_triangles.size(); i ++) {
        if (leaf_triangles[i]->is_leaf)
            //EXECUTION_REPORT(REPORT_ERROR, -1, leaf_triangles[i]->remained_points_in_triangle.size() == 0, "remap software error1 in triangularization_process");
            assert(leaf_triangles[i]->remained_points_in_triangle.size() == 0);
    }
    for (int i = 0; i < leaf_triangles.size(); i ++) {
        if (leaf_triangles[i]->is_leaf)
            continue;
        distribute_points_into_triangles(&(leaf_triangles[i]->remained_points_in_triangle), &leaf_triangles);
    }       
    for (int i = 0; i < leaf_triangles.size(); i ++)
        triangularization_process(leaf_triangles[i]);

    triangle->reference_count --;
}

/* This function should be call only once for the root virtual triangle,
 * becase it will alloc new memory for all points and cells. */
Triangle* Delaunay_Voronoi::initialize_super_triangle(int num_points, double *x, double *y, bool *redundant_cell_mark)
{
    double minX, maxX, minY, maxY;
    double dx, dy, deltaMax, midx, midy;
    Edge *e1, *e2, *e3;
    Triangle *super;

    assert(x != NULL);
    assert(y != NULL);

    maxX = minX = x[0];
    maxY = minY = y[0];
    for(int i = 0; i < num_points; i++) 
    {
        if (x[i] < minX) minX = x[i];
        if (x[i] > maxX) maxX = x[i];
        if (y[i] < minY) minY = y[i];
        if (y[i] > maxY) maxY = y[i];
    }
    
    dx = maxX - minX;
    dy = maxY - minY;
    deltaMax = std::max(dx, dy);
    midx = (minX + maxX) / 2.0;
    midy = (minY + maxY) / 2.0;

    virtual_point[0] = new Point(midx, midy + 20 * deltaMax , -1);
    virtual_point[1] = new Point(midx - 20 * deltaMax, midy - deltaMax, -1);
    virtual_point[2] = new Point(midx + 20 * deltaMax, midy - deltaMax, -1);

    //e1 = current_delaunay_voronoi->allocate_edge(virtual_point[1], virtual_point[0]);
    //e2 = current_delaunay_voronoi->allocate_edge(virtual_point[0], virtual_point[2]);
    //e3 = current_delaunay_voronoi->allocate_edge(virtual_point[2], virtual_point[1]);
    e1 = current_delaunay_voronoi->allocate_edge(virtual_point[0], virtual_point[1]);
    e2 = current_delaunay_voronoi->allocate_edge(virtual_point[1], virtual_point[2]);
    e3 = current_delaunay_voronoi->allocate_edge(virtual_point[2], virtual_point[0]);
    super = current_delaunay_voronoi->allocate_Triangle(e1, e2, e3);

    cells = new Cell[num_points];
    for (int i = 0; i < num_points; i ++) {
        Point *point = new Point(x[i], y[i], i);
        cells[i].center = point;
        if (!redundant_cell_mark[i]) {
            point->current_triangle = super;
            super->remained_points_in_triangle.push_back(point);
        }
    }

    return super;
}


void Delaunay_Voronoi::clear_triangle_containing_virtual_point()
{
    for(vector<Triangle*>::iterator t = result_leaf_triangles.begin(); t != result_leaf_triangles.end(); )
        if((*t)->contain_vertex(virtual_point[0]) || (*t)->contain_vertex(virtual_point[1]) || (*t)->contain_vertex(virtual_point[2]))
            result_leaf_triangles.erase(t);
        else
            t++;
}


Delaunay_Voronoi::Delaunay_Voronoi(int num_points, double *lat_values, double *lon_values, bool is_global_grid,
                   double min_lon, double max_lon, double min_lat, double max_lat, bool *redundant_cell_mark)
{
    Triangle *root;
    timeval start, end;
    int i, j;
    bool cyclic = min_lon==0 && max_lon==360;

    gettimeofday(&start, NULL);
    current_delaunay_voronoi = this;

    num_cells = num_points;

    this->is_global_grid = is_global_grid;

    root = initialize_super_triangle(num_points, lat_values, lon_values, redundant_cell_mark);

    triangularization_process(root);

    clear_triangle_containing_virtual_point();

    //generate_Voronoi_diagram();
    //extract_vertex_coordinate_values(num_points, output_vertex_lon_values, output_vertex_lat_values, output_num_vertexes);

    /* Below is for testing */
    gettimeofday(&end, NULL);
    printf("Delaunay time consuming: %ldms\n", ((end.tv_sec - start.tv_sec) * 1000000 + (end.tv_usec - start.tv_usec)) / 1000);   
}

Edge *Delaunay_Voronoi::allocate_edge(Point *head, Point *tail)
{
    Edge *new_edge = new Edge(head, tail);
    edge_pool.push_back(new_edge);

    return new_edge;
}


Triangle *Delaunay_Voronoi::allocate_Triangle(Point *point1, Point *point2, Point *point3)
{
    Triangle *new_triangle = new Triangle(point1, point2, point3);
    triangle_pool.push_back(new_triangle);

    return new_triangle;
}


Triangle *Delaunay_Voronoi::allocate_Triangle(Edge *edge1, Edge *edge2, Edge *edge3)
{
    Triangle *new_triangle = new Triangle(edge1, edge2, edge3);
    triangle_pool.push_back(new_triangle);

    return new_triangle;
}


vector<Edge*> Delaunay_Voronoi::get_all_delaunay_edge()
{
    vector<Edge*> all_edges;

    for(unsigned int i = 0; i < result_leaf_triangles.size(); i ++)
        if(result_leaf_triangles[i]->is_leaf) {
            all_edges.push_back(result_leaf_triangles[i]->edge[0]);
            all_edges.push_back(result_leaf_triangles[i]->edge[1]);
            all_edges.push_back(result_leaf_triangles[i]->edge[2]);
        }

    return all_edges;
}
