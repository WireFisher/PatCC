#include "mpi.h"
#include "delaunay_voronoi_2D.h"
#include "opencv_utils.h"
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <sys/time.h>
#include <tr1/unordered_map>
#include <list>
#include "merge_sort.h"


/*
 *           o Center
 *           ^
 *          /
 *         /
 *     V1 o -----> o V2
 */
double compute_three_2D_points_cross_product(double center_x, double center_y,
                                             double v1_x, double v1_y,
                                             double v2_x, double v2_y)
{
    double delta1_x = v2_x - v1_x;
    double delta1_y = v2_y - v1_y;
    double delta2_x = center_x - v1_x;
    double delta2_y = center_y - v1_y;

    return delta1_x*delta2_y - delta2_x*delta1_y;
}


double det(const Point *pt1, const Point *pt2, const Point *pt3)
{
    return compute_three_2D_points_cross_product(pt3->x, pt3->y, pt1->x, pt1->y, pt2->x, pt2->y);
}


static inline void swap(Point *p1, Point *p2)
{
    Point tmp = *p1;
    *p1 = *p2;
    *p2 = tmp;
}


static int compare_v2(const void* a, const void* b)
{
    Triangle_Transport t1 = *(const Triangle_Transport*)a;
    Triangle_Transport t2 = *(const Triangle_Transport*)b;

    if(t1.v[2].id < t2.v[2].id) return -1;
    if(t1.v[2].id > t2.v[2].id) return  1;
    return 0;
}


static int compare_v1(const void* a, const void* b)
{
    Triangle_Transport t1 = *(const Triangle_Transport*)a;
    Triangle_Transport t2 = *(const Triangle_Transport*)b;

    if(t1.v[1].id < t2.v[1].id) return -1;
    if(t1.v[1].id > t2.v[1].id) return  1;
    return 0;
}


static int compare_v0(const void* a, const void* b)
{
    Triangle_Transport t1 = *(const Triangle_Transport*)a;
    Triangle_Transport t2 = *(const Triangle_Transport*)b;

    if(t1.v[0].id < t2.v[0].id) return -1;
    if(t1.v[0].id > t2.v[0].id) return  1;
    return 0;
}


static int compare_lon(const void* a, const void* b)
{
    Triangle_Transport t1 = *(const Triangle_Transport*)a;
    Triangle_Transport t2 = *(const Triangle_Transport*)b;

    if(t1.v[0].x < t2.v[0].x) return -1;
    if(t1.v[0].x > t2.v[0].x) return  1;
    return 0;
}


static inline void radix_sort(Triangle_Transport *triangles, int num_triangles)
{
    assert(sizeof(Triangle_Transport) > sizeof(void *)/2);
    merge_sort(triangles, num_triangles, sizeof(Triangle_Transport), compare_lon);
    merge_sort(triangles, num_triangles, sizeof(Triangle_Transport), compare_v2);
    merge_sort(triangles, num_triangles, sizeof(Triangle_Transport), compare_v1);
    merge_sort(triangles, num_triangles, sizeof(Triangle_Transport), compare_v0);
}


void sort_points_in_triangle(Triangle_Transport *triangles, int num_triangles)
{
    for(int i = 0; i < num_triangles; i++) {
        if(triangles[i].v[0].id > triangles[i].v[1].id) swap(&triangles[i].v[0], &triangles[i].v[1]);
        if(triangles[i].v[1].id > triangles[i].v[2].id) swap(&triangles[i].v[1], &triangles[i].v[2]);
        if(triangles[i].v[0].id > triangles[i].v[1].id) swap(&triangles[i].v[0], &triangles[i].v[1]);
    }
}


void sort_triangles(Triangle_Transport *triangles, int num_triangles)
{
    radix_sort(triangles, num_triangles);
}


Point::Point()
{
}


Point::Point(double x, double y, int id)
{
    this->x = x;
    this->y = y;
    this->id = id;
}


double Point::calculate_distance(const Point *pt) const
{
    double dx = pt->x - x;
    double dy = pt->y - y;
    return sqrt(dx * dx + dy * dy);
}


double Point::calculate_distance(double pt_x, double pt_y) const
{
    double dx = pt_x - x;
    double dy = pt_y - y;
    return sqrt(dx * dx + dy * dy);
}


Point::~Point()
{
}


Point operator - (Point p1, Point p2)
{
    return Point(p1.x - p2.x, p1.y - p2.y, -1);
}


bool operator == (Point p1, Point p2)
{
    if(p1.id != -1 && p2.id != -1)
        return p1.id == p2.id;

    if(std::abs(p1.x - p2.x) < FLOAT_ERROR && std::abs(p1.y - p2.y) < FLOAT_ERROR )
        return true;
    else if(std::abs(std::abs(p1.x - p2.x) - 360.0) < FLOAT_ERROR && std::abs(p1.y - p2.y) < FLOAT_ERROR )
        return true;
    else
        return false;
}


bool operator != (Point p1, Point p2)
{
    return !(p1 == p2);
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

    if (std::fabs(res1) < FLOAT_ERROR_HI)
        return 0;
    else if (res1 > 0)
        return 1;
    else return -1;
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
#ifdef DDEBUG
    bool on1 = position_to_edge(triangle->v[0], triangle->v[1]) == 0;
    bool on2 = position_to_edge(triangle->v[1], triangle->v[2]) == 0;
    bool on3 = position_to_edge(triangle->v[2], triangle->v[0]) == 0;
    assert(!(on1 && on2));
    assert(!(on2 && on3));
    assert(!(on3 && on1));
#endif
    int ret = 0;
    int pos = position_to_edge(triangle->v[0], triangle->v[1]);
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


inline int Point::is_in_region(double min_x, double max_x, double min_y, double max_y) const
{
    return x > min_x && x < max_x && y > min_y && y < max_y;
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


Edge::~Edge()
{
}


Edge* Delaunay_Voronoi::generate_twins_edge(Edge *e)
{
    Edge *twins_edge = allocate_edge(e->tail, e->head);
    twins_edge->twin_edge = e;
    e->twin_edge = twins_edge;

    return twins_edge;
}


double Delaunay_Voronoi::calculate_angle(const Point *center, const Point *p1, const Point *p2)
{
    double ax, ay, bx, by;
    assert(center != p1 && center != p2) ;

    ax = p1->x - center->x;
    ay = p1->y - center->y;
    bx = p2->x - center->x;
    by = p2->y - center->y;
    return acos((ax*bx+ay*by) / sqrt((ax*ax+ay*ay)*(bx*bx+by*by)));
}


const Point *Delaunay_Voronoi::get_lowest_point_of_four(const Point *p1, const Point *p2, const Point *p3, const Point *p4)
{
    const Point *p = p1;
    if(!x_ref) {
        if(p2->x < p->x || (p2->x == p->x && p2->y < p->y)) p = p2;
        if(p3->x < p->x || (p3->x == p->x && p3->y < p->y)) p = p3;
        if(p4->x < p->x || (p4->x == p->x && p4->y < p->y)) p = p4;
    } else {
        if(x_ref[p2->id] < x_ref[p->id] || (x_ref[p2->id] == x_ref[p->id] && y_ref[p2->id] < y_ref[p->id])) p = p2;
        if(x_ref[p3->id] < x_ref[p->id] || (x_ref[p3->id] == x_ref[p->id] && y_ref[p3->id] < y_ref[p->id])) p = p3;
        if(x_ref[p4->id] < x_ref[p->id] || (x_ref[p4->id] == x_ref[p->id] && y_ref[p4->id] < y_ref[p->id])) p = p4;
    }
    return p;
}


bool Delaunay_Voronoi::is_angle_ambiguous(const Point *pt, const Edge *edge)
{
    if(pt == get_lowest_point_of_four(pt, edge->head, edge->tail, edge->twin_edge->prev_edge_in_triangle->head) ||
       edge->twin_edge->prev_edge_in_triangle->head == get_lowest_point_of_four(pt, edge->head, edge->tail, edge->twin_edge->prev_edge_in_triangle->head)) {
        if(edge->triangle->is_cyclic)
            return true;
        else
            return false;
    } else {
        if(edge->triangle->is_cyclic)
            return false;
        else
            return true;
    }
}


/* Debug staff */
int triangulate_count = 0;
bool Print_Error_info=false;
int on_circle_count = 0;
int not_on_circle_count = 0;
extern double global_p_lon[4];
extern double global_p_lat[4];
#define point_is(pt, a, b) (fabs(pt->x - a) < 1e-5 && fabs(pt->y - b) < 1e-5)
#define edge_is(e, a, b, c, d) ((point_is(e->head, a, b) && point_is(e->tail, c, d)) || (point_is(e->tail, a, b) && point_is(e->head, c, d)))
#define triangle_is(t, a, b, c, d ,e, f) ((point_is(t->v[0], a, b) || point_is(t->v[0], c, d) || point_is(t->v[0], e, f)) && \
                                          (point_is(t->v[1], a, b) || point_is(t->v[1], c, d) || point_is(t->v[1], e, f)) && \
                                          (point_is(t->v[2], a, b) || point_is(t->v[2], c, d) || point_is(t->v[2], e, f)) )
/*
    bool is = false;
    for(int i = 0 ;i < 4; i++) {
        if(point_is(pt, global_p_lon[i], global_p_lat[i]))
            is = true;
    }
    if(is) {
        is = false;
        for(int i = 0 ;i < 4; i++) {
            if(edge_is(edge, global_p_lon[i], global_p_lat[i], global_p_lon[(i+1)%4], global_p_lat[(i+1)%4]))
                is = true;
        }
        if(edge_is(edge, global_p_lon[0], global_p_lat[0], global_p_lon[2], global_p_lat[2]))
            is = true;
        if(edge_is(edge, global_p_lon[1], global_p_lat[1], global_p_lon[3], global_p_lat[3]))
            is = true;
    }
*/

enum REASON {
    SUCCESS,
    IN_CIRCUMCIRCLE,
    AMBIGUOUS
} reason;
bool Delaunay_Voronoi::is_triangle_legal(const Point *pt, const Edge *edge)
{
#ifdef DEBUG
    reason = SUCCESS;
#endif
    if (!edge->twin_edge) {
        return true;
    }

    if(!edge->twin_edge->triangle) {
        return true;
    }

    if(!edge->twin_edge->triangle->is_leaf) {
        return true;
    }

    int ret = edge->triangle->circum_circle_contains(edge->twin_edge->prev_edge_in_triangle->head, tolerance);
    if (ret == -1) {
        return true;
    }

    if (ret == 0) {
        return !is_angle_too_large(pt, edge);
    }

#ifdef DEBUG
    reason = IN_CIRCUMCIRCLE;
#endif
    return false;
}

bool Delaunay_Voronoi::is_angle_too_large(const Point *pt, const Edge *edge)
{
    if(pt == get_lowest_point_of_four(pt, edge->head, edge->tail, edge->twin_edge->prev_edge_in_triangle->head) ||
       edge->twin_edge->prev_edge_in_triangle->head == get_lowest_point_of_four(pt, edge->head, edge->tail, edge->twin_edge->prev_edge_in_triangle->head))
        return false;
    else {
#ifdef DEBUG
        reason = AMBIGUOUS;
#endif
        return true;
    }
}

bool Delaunay_Voronoi::is_triangle_ambiguous(const Point *pt, const Edge *edge)
{
    if (!edge->twin_edge) {
        return false;
    }

    if(!edge->twin_edge->triangle) {
        return false;
    }

    if(!edge->twin_edge->triangle->is_leaf) {
        return false;
    }

    int ret = edge->triangle->circum_circle_contains(edge->twin_edge->prev_edge_in_triangle->head);
    if (ret == 0) {
        return is_angle_ambiguous(pt, edge);
    }

    return false;
}


bool Delaunay_Voronoi::is_triangle_legal(const Triangle *t)
{
    for(int i = 0; i < 3; i++)
        if(!is_triangle_legal(t->edge[i]->prev_edge_in_triangle->head, t->edge[i])) {
            printf("[%d] +illegal triangle: (%lf, %lf), (%lf, %lf), (%lf, %lf)\n", 1, t->v[0]->x, t->v[0]->y, t->v[1]->x, t->v[1]->y, t->v[2]->x, t->v[2]->y);
            Triangle *tt = t->edge[i]->twin_edge->triangle;
            printf("[%d] -illegal triangle: (%lf, %lf), (%lf, %lf), (%lf, %lf)\n", 1, tt->v[0]->x, tt->v[0]->y, tt->v[1]->x, tt->v[1]->y, tt->v[2]->x, tt->v[2]->y);
            printf("[%d] +: %d, -: %d\n", 1, t->is_leaf, tt->is_leaf);
            printf("===============================================================\n");
            return false;
        }
    return true;
}


bool Delaunay_Voronoi::is_all_leaf_triangle_legal()
{
    //Print_Error_info = true;
    //int rank;
    //MPI_Comm_rank(process_thread_mgr->get_mpi_comm(), &rank);

    bool is_legal = true;
    for(unsigned i = 0; i < result_leaf_triangles.size(); i++)
        if(result_leaf_triangles[i]->is_leaf)
            if(!is_triangle_legal(result_leaf_triangles[i])) {
                is_legal = false;
                //printf("illegal reason: %d\n", reason);
                //fprintf(stderr, "[%d] illegal: (%lf, %lf), (%lf, %lf), (%lf, %lf)\n", rank, result_leaf_triangles[i]->v[0]->x, result_leaf_triangles[i]->v[0]->y, result_leaf_triangles[i]->v[1]->x, result_leaf_triangles[i]->v[1]->y, result_leaf_triangles[i]->v[2]->x, result_leaf_triangles[i]->v[2]->y);
            }
    if(is_legal)
        return true;
    else
        return false;
}


void Delaunay_Voronoi::legalize_triangles(Point *vr, Edge *edge, vector<Triangle*> *leaf_triangles)
{
    //int rank;
    //MPI_Comm_rank(process_thread_mgr->get_mpi_comm(), &rank);
    //if(triangulate_count == 261 && rank == 79)
    //    printf("(%lf, %lf)->(%lf, %lf)-(%lf, %lf): %d\n", vr->x, vr->y, edge->head->x, edge->head->y, edge->tail->x, edge->tail->y, is_triangle_legal(vr, edge));
    if (is_triangle_legal(vr, edge))
        return;

    //EXECUTION_REPORT(REPORT_ERROR, -1, edge->triangle->is_leaf, "remap software error1 in legalize_triangles\n");
    assert(edge->triangle->is_leaf);
    //EXECUTION_REPORT(REPORT_ERROR, -1, edge->twin_edge->triangle->is_leaf, "remap software error2 in legalize_triangles %lx\n", (long)(edge->twin_edge->triangle));
    assert(edge->twin_edge->triangle->is_leaf);
    leaf_triangles->push_back(edge->twin_edge->triangle);
    edge->triangle->is_leaf = false;
    edge->twin_edge->triangle->is_leaf = false;

    Point *vk = edge->twin_edge->prev_edge_in_triangle->head;
    Edge *eij = edge;
    Edge *ejr = eij->next_edge_in_triangle;
    Edge *eri = ejr->next_edge_in_triangle;
    Edge *eji = eij->twin_edge;
    Edge *eik = eji->next_edge_in_triangle;
    Edge *ekj = eik->next_edge_in_triangle;
    Edge *erk = allocate_edge(vr, vk);
    Edge *ekr = generate_twins_edge(erk);
    Triangle* tikr = allocate_Triangle(eik,ekr,eri);
    Triangle* tjrk = allocate_Triangle(ejr,erk,ekj);
    leaf_triangles->push_back(tikr);
    leaf_triangles->push_back(tjrk);
    legalize_triangles(vr, eik, leaf_triangles);
    legalize_triangles(vr, ekj, leaf_triangles);
}


void Delaunay_Voronoi::relegalize_triangles(Point *vr, Edge *edge)
{
    if (!is_triangle_ambiguous(vr, edge))
        return;

    if (edge->triangle->is_cyclic != edge->twin_edge->triangle->is_cyclic)
        return;
    //EXECUTION_REPORT(REPORT_ERROR, -1, edge->triangle->is_leaf, "remap software error1 in legalize_triangles\n");
    assert(edge->triangle->is_leaf);
    //EXECUTION_REPORT(REPORT_ERROR, -1, edge->twin_edge->triangle->is_leaf, "remap software error2 in legalize_triangles %lx\n", (long)(edge->twin_edge->triangle));
    assert(edge->twin_edge->triangle->is_leaf);

    Triangle *triangle_origin, *triangle_twin;
    triangle_origin = edge->triangle;
    triangle_twin = edge->twin_edge->triangle;

    Point *vk = edge->twin_edge->prev_edge_in_triangle->head;
    Edge *eij = edge;
    Edge *ejr = eij->next_edge_in_triangle;
    Edge *eri = ejr->next_edge_in_triangle;
    Edge *eji = eij->twin_edge;
    Edge *eik = eji->next_edge_in_triangle;
    Edge *ekj = eik->next_edge_in_triangle;
    Edge *erk = allocate_edge(vr, vk);
    Edge *ekr = generate_twins_edge(erk);

    bool force = edge->triangle->is_cyclic ? true : false;
    triangle_origin->initialize_triangle_with_edges(eik, ekr, eri, force);
    triangle_twin->initialize_triangle_with_edges(ejr, erk, ekj, force);
    if(force)
        triangle_origin->is_cyclic = triangle_twin->is_cyclic = true;

    //relegalize_triangles(vr, eik);
    //relegalize_triangles(vr, ekj);
}


Triangle::Triangle()
{
    edge[0] = NULL;
    edge[1] = NULL;
    edge[2] = NULL;
}


Triangle::Triangle(Point *p1, Point *p2, Point *p3)
{
    initialize_triangle_with_edges(new Edge(p1, p2), new Edge(p2, p3), new Edge(p3, p1));
}


Triangle::Triangle(Edge *edge1, Edge *edge2, Edge *edge3)
{
    initialize_triangle_with_edges(edge1, edge2, edge3);
}


Triangle::~Triangle()
{
}


double Triangle::area()
{
    double a = v[0]->calculate_distance(v[1]);
    double b = v[1]->calculate_distance(v[2]);
    double c = v[2]->calculate_distance(v[0]);
    double p = (a + b + c) * 0.5;
    return sqrt(p*(p-a)*(p-b)*(p-c));
}


void Triangle::check_and_set_twin_edge_relationship(Triangle *another_triangle)
{
    for (int i = 0; i < 3; i ++)
        for (int j = 0; j < 3; j ++)
            if (edge[i]->head == another_triangle->edge[j]->tail && edge[i]->tail == another_triangle->edge[j]->head) {
                edge[i]->twin_edge = another_triangle->edge[j];
                another_triangle->edge[j]->twin_edge = edge[i];
            }
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
    circum_radius = sqrt(((v[0]->x - circum_center[0]) * (v[0]->x - circum_center[0])) + ((v[0]->y - circum_center[1]) * (v[0]->y - circum_center[1])));
}


void Triangle::initialize_triangle_with_edges(Edge *edge1, Edge *edge2, Edge *edge3, bool force)
{
    if(force) {
        v[0] = edge1->head;
        v[1] = edge2->head;
        v[2] = edge3->head;
        edge[0] = edge1;
        edge[1] = edge2;
        edge[2] = edge3;
    } else {
        Point *pt1, *pt2, *pt3;

        is_leaf = true;
        is_cyclic = false;
        
        pt1 = edge1->head;
        pt2 = edge2->head;
        pt3 = edge3->head;

        //EXECUTION_REPORT(REPORT_ERROR, -1, fabs(det(pt1, pt2, pt3)) > e && fabs(det(pt2, pt3, pt1)) > e && fabs(det(pt3, pt1, pt2)) > e,
        //                 "points given to construct triangle are on the same line.");

#ifdef DEBUG
        if(!(std::fabs(det(pt1, pt2, pt3)) > FLOAT_ERROR_HI && std::fabs(det(pt2, pt3, pt1)) > FLOAT_ERROR_HI && std::fabs(det(pt3, pt1, pt2)) > FLOAT_ERROR_HI)) {
            printf("(%.20lf, %.20lf), (%.20lf, %.20lf), (%.20lf, %.20lf)\n", pt1->x, pt1->y, pt2->x, pt2->y, pt3->x, pt3->y);
            printf("std::fabs(det(pt1, pt2, pt3)): %.20lf\nstd::fabs(det(pt2, pt3, pt1)): %.20lf\nstd::fabs(det(pt3, pt1, pt2)): %.20lf\n", std::fabs(det(pt1, pt2, pt3)),
                                                                                                                               std::fabs(det(pt2, pt3, pt1)),
                                                                                                                               std::fabs(det(pt3, pt1, pt2)));
        }
#endif
        /* if there are unmarked redundant points, the assertion may fail */
        assert(std::fabs(det(pt1, pt2, pt3)) > FLOAT_ERROR_HI && std::fabs(det(pt2, pt3, pt1)) > FLOAT_ERROR_HI && std::fabs(det(pt3, pt1, pt2)) > FLOAT_ERROR_HI);
        //EXECUTION_REPORT(REPORT_ERROR, -1, edge1->tail==edge2->head && edge2->tail==edge3->head && edge3->tail==edge1->head, "edges given to construct triangle is invalid.");
#ifdef DEBUG
        if(!(edge1->tail==edge2->head && edge2->tail==edge3->head && edge3->tail==edge1->head)) {
            printf("edge1: %p->%p, edge2: %p->%p, edge3: %p->%p\n", edge1->head, edge1->tail, edge2->head, edge2->tail, edge3->head, edge3->tail);
            printf("1: (%lf, %lf)->(%lf, %lf) 2: (%lf, %lf)->(%lf, %lf) 3:(%lf, %lf)->(%lf, %lf)\n", edge1->head->x, edge1->head->y, edge1->tail->x, edge1->tail->y,
                                                                                                     edge2->head->x, edge2->head->y, edge2->tail->x, edge2->tail->y,
                                                                                                     edge3->head->x, edge3->head->y, edge3->tail->x, edge3->tail->y);
        }
#endif
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
            //int rank;
            //MPI_Comm_rank(process_thread_mgr->get_mpi_comm(), &rank);
            //printf("[%d] assert false\n", rank);
            assert(false);
            v[1] = pt3;
            v[2] = pt2;
            this->edge[0] = edge3->twin_edge;
            this->edge[1] = edge2->twin_edge;
            this->edge[2] = edge1->twin_edge;
            //EXECUTION_REPORT(REPORT_ERROR, -1, edge3->twin_edge != NULL && edge2->twin_edge != NULL && edge1->twin_edge != NULL, "remap software error3 in new Triangle");
            assert(edge3->twin_edge != NULL && edge2->twin_edge != NULL && edge1->twin_edge != NULL);
        }
    }
    remained_points = NULL;
    num_remained_points = 0;
    
    this->edge[0]->next_edge_in_triangle = this->edge[1];
    this->edge[1]->next_edge_in_triangle = this->edge[2];
    this->edge[2]->next_edge_in_triangle = this->edge[0];
    this->edge[0]->prev_edge_in_triangle = this->edge[2];
    this->edge[1]->prev_edge_in_triangle = this->edge[0];
    this->edge[2]->prev_edge_in_triangle = this->edge[1];

    this->edge[0]->triangle = this;
    this->edge[1]->triangle = this;
    this->edge[2]->triangle = this;
    calulate_circum_circle();
}


/*
 * Input : Point to be checked
 * Return:  1    point is in circum circle
 *          0    point is on circum circle
 *         -1    point is out of circum circle
 */
int Triangle::circum_circle_contains(Point *p, double tolerance)
{
    calulate_circum_circle();
    double dist2 = ((p->x - circum_center[0]) * (p->x - circum_center[0])) + ((p->y - circum_center[1]) * (p->y - circum_center[1]));
    if(std::fabs(dist2 - circum_radius*circum_radius) < tolerance &&
       really_on_circum_circle(p, tolerance))
        return 0;
    else if(dist2 < circum_radius*circum_radius)
        return 1;
    else // (dist > circum_radius)
        return -1;
}


bool Triangle::really_on_circum_circle(Point *p, double tolerance)
{
    Point *pt[4];

    for(int i = 0; i < 3; i++)
        pt[i] = v[i];
    pt[3] = p;

    for(int j = 3; j > 0; j--)
        for(int i = 0; i < j; i++)
            if(pt[i]->id > pt[i+1]->id) {
                Point *tmp = pt[i];
                pt[i] = pt[i+1];
                pt[i+1] = tmp;
            }

    double ab = (pt[0]->x * pt[0]->x) + (pt[0]->y * pt[0]->y);
    double cd = (pt[1]->x * pt[1]->x) + (pt[1]->y * pt[1]->y);
    double ef = (pt[2]->x * pt[2]->x) + (pt[2]->y * pt[2]->y);

    double center[2];
    center[0] = (ab * (pt[2]->y - pt[1]->y) + cd * (pt[0]->y - pt[2]->y) + ef * (pt[1]->y - pt[0]->y)) /
                       (pt[0]->x * (pt[2]->y - pt[1]->y) + pt[1]->x * (pt[0]->y - pt[2]->y) + pt[2]->x * (pt[1]->y - pt[0]->y)) / 2.f;
    center[1] = (ab * (pt[2]->x - pt[1]->x) + cd * (pt[0]->x - pt[2]->x) + ef * (pt[1]->x - pt[0]->x)) /
                       (pt[0]->y * (pt[2]->x - pt[1]->x) + pt[1]->y * (pt[0]->x - pt[2]->x) + pt[2]->y * (pt[1]->x - pt[0]->x)) / 2.f;
    double radius2 = ((pt[0]->x - center[0]) * (pt[0]->x - center[0])) + ((pt[0]->y - center[1]) * (pt[0]->y - center[1]));
    double dist2 = ((pt[3]->x - center[0]) * (pt[3]->x - center[0])) + ((pt[3]->y - center[1]) * (pt[3]->y - center[1]));
    return std::fabs(dist2 - radius2) < tolerance;
}


int Triangle::find_best_candidate_point()
{
    double min_dist=1e10, dist;
    int best_candidate_id;

    
    if (num_remained_points == 0)
        return -1;

    for (unsigned i = 0; i < num_remained_points; i ++) {
        dist = remained_points[i]->calculate_distance((v[0]->x+v[1]->x+v[2]->x)/3.0, (v[0]->y+v[1]->y+v[2]->y)/3.0);
        if (min_dist > dist) {
            min_dist = dist;
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


void Delaunay_Voronoi::swap_link_node(int idx1, int idx2)
{
    double tmp_x  = all_points[idx1].x;
    double tmp_y  = all_points[idx1].y;
    double tmp_id = all_points[idx1].id;

    all_points[idx1].x  = all_points[idx2].x;
    all_points[idx1].y  = all_points[idx2].y;
    all_points[idx1].id = all_points[idx2].id;

    all_points[idx2].x  = tmp_x;
    all_points[idx2].y  = tmp_y;
    all_points[idx2].id = tmp_id;
}


void Delaunay_Voronoi::distribute_points_into_triangles(int head, int tail, vector<Triangle*> *triangles)
{
    if (num_pnts < 1)
        return;

    int start = head;
    for (unsigned i = 0; i < triangles->size(); i++) {
        if (!(*triangles)[i]->is_leaf)
            continue;

        int end = tail;
        for (int j = start; j <= end;) {
            if (all_points[j] && all_points[j]->position_to_triangle((*triangles)[i]) >= 0) {
                j = all_points[j].next;
            } else {
                swap(pnts[j], pnts[end]);
                end = all_points[end].prev;
            }
        }
        (*triangles)[i]->set_remained_points(start, end);
        start = end + 1;
    }

    for (unsigned i = start; i < num_pnts; i++)
        assert(!pnts[i]);
    return;
}


void Delaunay_Voronoi::triangularization_process(Triangle *triangle)
{
    int best_candidate_point_id;
    Point *best_candidate_point;
    vector<Triangle *> leaf_triangles;

#ifdef DEBUG
    assert(triangle->is_leaf);
#endif
    if(!triangle->is_leaf)
        return;

    if (triangle->num_remained_points == 0) {
        result_leaf_triangles.push_back(triangle);
        return;
    }

    triangle->is_leaf = false;

    best_candidate_point_id = triangle->find_best_candidate_point();
    best_candidate_point = triangle->remained_points[best_candidate_point_id];

    triangle->remained_points[best_candidate_point_id] = NULL;

    if (best_candidate_point->position_to_triangle(triangle) == 0) { //inside
        Edge *e_v1_can = allocate_edge(triangle->v[0], best_candidate_point);
        Edge *e_can_v1 = generate_twins_edge(e_v1_can);
        Edge *e_v2_can = allocate_edge(triangle->v[1], best_candidate_point);
        Edge *e_can_v2 = generate_twins_edge(e_v2_can);
        Edge *e_v3_can = allocate_edge(triangle->v[2], best_candidate_point);
        Edge *e_can_v3 = generate_twins_edge(e_v3_can);
        Triangle *t_can_v1_v2 = allocate_Triangle(e_can_v1, triangle->edge[0], e_v2_can);
        Triangle *t_can_v2_v3 = allocate_Triangle(e_can_v2, triangle->edge[1], e_v3_can);
        Triangle *t_can_v3_v1 = allocate_Triangle(e_can_v3, triangle->edge[2], e_v1_can);
        leaf_triangles.push_back(triangle);
        leaf_triangles.push_back(t_can_v1_v2);
        leaf_triangles.push_back(t_can_v2_v3);
        leaf_triangles.push_back(t_can_v3_v1);
        
        /* t_can_v1_v2->v[0] is best_candidate_point, actually */
        legalize_triangles(t_can_v1_v2->v[0], t_can_v1_v2->edge[1], &leaf_triangles);
        legalize_triangles(t_can_v2_v3->v[0], t_can_v2_v3->edge[1], &leaf_triangles);
        legalize_triangles(t_can_v3_v1->v[0], t_can_v3_v1->edge[1], &leaf_triangles);
    }
    else { // on the side
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
                //EXECUTION_REPORT(REPORT_ERROR, -1, false, "point, which should be found in triangle, is outside of triangle");
                assert(false);
                break;
        }
        //EXECUTION_REPORT(REPORT_ERROR, -1, best_candidate_point->position_to_edge(vi, vj) == 0, "point, which should be on the edge, is not on the edge");
        assert(best_candidate_point->position_to_edge(vi, vj) == 0);
        if (eij->twin_edge != NULL)
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
        Edge *eir = allocate_edge(vi, best_candidate_point);
        Edge *erk = allocate_edge(best_candidate_point, vk);
        Edge *ekr = generate_twins_edge(erk);
        Edge *erj = allocate_edge(best_candidate_point, vj);
        Triangle* tirk = allocate_Triangle(eir, erk, eki);
        Triangle* tjkr = allocate_Triangle(ejk, ekr, erj);
        leaf_triangles.push_back(triangle);
        leaf_triangles.push_back(tirk);
        leaf_triangles.push_back(tjkr);
        legalize_triangles(best_candidate_point, ejk, &leaf_triangles);
        legalize_triangles(best_candidate_point, eki, &leaf_triangles); 
        if (eij->twin_edge != NULL) {
            Edge *eri = generate_twins_edge(eir);
            Edge *ejr = generate_twins_edge(erj);
            Edge *erl = allocate_edge(best_candidate_point, vl);
            Edge *elr = generate_twins_edge(erl);
            Triangle* tilr = allocate_Triangle(eil, elr, eri);
            Triangle* tjrl = allocate_Triangle(ejr, erl, elj);
            leaf_triangles.push_back(eij->twin_edge->triangle);
            leaf_triangles.push_back(tilr);
            leaf_triangles.push_back(tjrl);
            legalize_triangles(best_candidate_point, eil, &leaf_triangles);
            legalize_triangles(best_candidate_point, elj, &leaf_triangles);
        }
        else {
            //eir->twin_edge = NULL;
        }
    }

#ifdef DEBUG
    for (unsigned i = 0; i < leaf_triangles.size(); i ++) {
        if (leaf_triangles[i]->is_leaf)
            assert(leaf_triangles[i]->num_remained_points == 0);
    }
#endif

    Point **combined_buf = NULL;
    int num_combined = 0;
    merge_buffer(&leaf_triangles, &combined_buf, &num_combined);

    distribute_points_into_triangles(combined_buf, num_combined, &leaf_triangles);

        //char filename[64];
        //snprintf(filename, 64, "log/single_step/step_%d.png", triangulate_count++);
        //plot_current_step_into_file(filename);
        //printf("plot step %d\n", triangulate_count);
        //triangulate_count++;

    for (unsigned i = 0; i < leaf_triangles.size(); i ++)
        if(leaf_triangles[i]->is_leaf)
            triangularization_process(leaf_triangles[i]);

}


void Delaunay_Voronoi::merge_buffer(vector<Triangle*> *ts, Point*** out_buf, int *out_num)
{
    Point** buf2 = *out_buf;
    int     num2 = *out_num;

    vector<Triangle*> non_leaf;
    for (unsigned i = 0; i < ts->size(); i ++)
        if (!(*ts)[i]->is_leaf && (*ts)[i]->remained_points)
            non_leaf.push_back((*ts)[i]);

    printf("=========================\n");
    for (unsigned i = 0; i < non_leaf.size(); i ++)
        printf("%p, %p\n", non_leaf[i]->remained_points, non_leaf[i]->remained_points+non_leaf[i]->num_remained_points);
    printf("+++++++++++++++++++++++++\n");
    bool changed = true;
    while(changed) {
        changed = false;
        for (unsigned i = 0; i < non_leaf.size(); i ++) {
            Point** buf1 = non_leaf[i]->remained_points;
            int     num1 = non_leaf[i]->num_remained_points;
            if (buf2 == NULL) {
                buf2 = buf1;
                num2 = num1;
                changed = true;
            }
            else if (buf1 + num1 == buf2) {
                buf2 = buf1;
                num2 += num1;
                changed = true;
            }
            else if (buf2 + num2 == buf1) {
                num2 += num1;
                changed = true;
            }

            if (changed) {
                non_leaf.erase(non_leaf.begin() + i);
                break;
            }
        }
    }
    if (non_leaf.size() > 0) {
        if ()
    }
    *out_buf = buf2;
    *out_num = num2;
    printf("%p, %p\n", buf2, buf2+num2);
    printf("=========================\n");
}


/* This function should be call only once for the root virtual triangle,
 * becase it will alloc new memory for all points and cells. */
#define PDLN_INSERT_EXTRA_VPOINT (false)
#define PDLN_VPOINT_DENSITY  (20)
vector<Triangle*> Delaunay_Voronoi::generate_initial_triangles(int num_points, double *x, double *y, bool *redundant_cell_mark)
{
    double minX, maxX, minY, maxY;
    double dx, dy, deltaMax;
    vector<Triangle *> virtual_triangles;

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

    double ratio  = 0.1;
    double v_minx = minX-deltaMax*ratio;
    double v_maxx = maxX+deltaMax*ratio;
    double v_miny = minY-deltaMax*ratio;
    double v_maxy = maxY+deltaMax*ratio;

    virtual_point[0] = new Point(minX-deltaMax*ratio, minY-deltaMax*ratio, -1);
    virtual_point[1] = new Point(minX-deltaMax*ratio, maxY+deltaMax*ratio, -1);
    virtual_point[2] = new Point(maxX+deltaMax*ratio, minY-deltaMax*ratio, -1);
    virtual_point[3] = new Point(maxX+deltaMax*ratio, maxY+deltaMax*ratio, -1);

    virtual_triangles.push_back(allocate_Triangle(allocate_edge(virtual_point[0], virtual_point[2]),
                                                  allocate_edge(virtual_point[2], virtual_point[1]),
                                                  allocate_edge(virtual_point[1], virtual_point[0])));
    virtual_triangles.push_back(allocate_Triangle(allocate_edge(virtual_point[3], virtual_point[1]),
                                                  allocate_edge(virtual_point[1], virtual_point[2]),
                                                  allocate_edge(virtual_point[2], virtual_point[3])));
    virtual_triangles[0]->edge[1]->twin_edge = virtual_triangles[1]->edge[1];
    virtual_triangles[1]->edge[1]->twin_edge = virtual_triangles[0]->edge[1];

    all_points = new Point[num_points];

    for (int i = 0; i < num_points; i ++) {
        all_points[i].x    = x[i];
        all_points[i].y    = y[i];
        all_points[i].id   = i;
        all_points[i].next = i + 1;
        all_points[i].prev = i - 1;
    }
    all_points[0].prev = all_points[num_points-1].next = -1;

    distribute_points_into_triangles(0, num_points-1, &virtual_triangles);

    return virtual_triangles;
}


void Delaunay_Voronoi::clear_triangle_containing_virtual_point()
{
    extra_virtual_point.push_back(virtual_point[0]);
    extra_virtual_point.push_back(virtual_point[1]);
    extra_virtual_point.push_back(virtual_point[2]);
    extra_virtual_point.push_back(virtual_point[3]);
    for(vector<Triangle*>::iterator t = result_leaf_triangles.begin(); t != result_leaf_triangles.end(); ) {
        bool contain = false;
        for(unsigned i = 0; i < extra_virtual_point.size(); i++)
            if((*t)->is_leaf && (*t)->contain_vertex(extra_virtual_point[i])) {
                for(unsigned j = 0; j < 3; j++) {
                    (*t)->edge[j]->triangle = NULL;
                    if((*t)->edge[j]->twin_edge)
                        (*t)->edge[j]->twin_edge->twin_edge = NULL;
                }
                (*t)->is_leaf = false;
                result_leaf_triangles.erase(t); //TODO: erase is slow
                contain = true;
                break;
            }
        if(!contain)
            t++;
    }
}


Delaunay_Voronoi::Delaunay_Voronoi(int num_points, double *x_values, double *y_values, const double *x_origin, const double *y_origin,
                                   int *global_idx, bool is_global_grid, double min_lon, double max_lon, 
                                   double min_lat, double max_lat, bool *redundant_cell_mark, int virtual_polar_local_index)
    : tolerance(FLOAT_ERROR)
    , x_ref(x_origin)
    , y_ref(y_origin)
{
    timeval start, end;

    assert((x_ref && y_ref) || (!x_ref && !x_ref));
#ifdef DEBUG
    assert(have_redundent_points(x_values, y_values, num_points) == false);
#endif
    int triangles_count_estimate = 2*num_points;
    triangle_pool.reserve(triangles_count_estimate);
    edge_pool.reserve(triangles_count_estimate*3);

    gettimeofday(&start, NULL);

    num_cells = num_points;

    this->is_global_grid = is_global_grid;

    if(global_idx != NULL) {
        global_index = new int[num_points];
        memcpy(global_index, global_idx, num_points*sizeof(int));
    }

    this->vpolar_local_index = virtual_polar_local_index;
    vector<Triangle*> initial_triangles = generate_initial_triangles(num_points, x_values, y_values, redundant_cell_mark);

    //save_original_points_into_file();

    for(unsigned i = 0; i < initial_triangles.size(); i++)
        triangularization_process(initial_triangles[i]);

    clear_triangle_containing_virtual_point();

    update_virtual_polar_info();
    //generate_Voronoi_diagram();
    //extract_vertex_coordinate_values(num_points, output_vertex_lon_values, output_vertex_lat_values, output_num_vertexes);

    //if(!is_all_leaf_triangle_legal())
    //    printf("warning: illegal\n");

    gettimeofday(&end, NULL);
}


Delaunay_Voronoi::~Delaunay_Voronoi()
{
    for (int i = 0; i < num_cells; i ++)
        delete cells[i].center;
    delete [] cells;
    /*
    for (unsigned i = 0; i < edge_pool.size(); i ++)
        delete edge_pool[i];
        //edge_allocator.deleteElement(edge_pool[i]);
    for (unsigned i = 0; i < triangle_pool.size(); i ++)
        delete triangle_pool[i];
        //triangle_allocator.deleteElement(triangle_pool[i]);
    */
}


void Delaunay_Voronoi::check_and_set_twin_edge_relationship(vector<Triangle*> *triangles)
{
    for (unsigned i = 0; i < triangles->size(); i ++)
        for (unsigned j = i+1; j < triangles->size(); j ++)
            (*triangles)[i]->check_and_set_twin_edge_relationship((*triangles)[j]);
}


bool have_redundent_points(const double *x, const double *y, int num)
{
    std::tr1::unordered_map<double, std::list<int> > hash_table;
    std::tr1::unordered_map<double, std::list<int> >::iterator it_hash;

    if(num == 0)
        return false;

    bool have_redundent = false;
    for(int i = 0; i < num; i++) {
        it_hash = hash_table.find(x[i] * 1000.0 + y[i]);
        if(it_hash != hash_table.end()) {
            bool same = false;
            for(std::list<int>::iterator it_list = it_hash->second.begin(); it_list != it_hash->second.end(); it_list ++)
                if(x[*it_list] == x[i] && y[*it_list] == y[i]) {
                    same = true;
                    break;
                }
            if(same){
                printf("redundent_point: %lf, %lf\n", x[i], y[i]);
                have_redundent = true;
            }
            else {
                it_hash->second.push_back(i);
            }
        }
        else {
            hash_table[x[i] * 1000.0 + y[i]].push_back(i);
        }
    }

    if(have_redundent)
        return true;

    return false;
}


void delete_redundent_points(double *&x, double *&y, int &num)
{
    std::tr1::unordered_map<double, std::list<int> > hash_table;
    std::tr1::unordered_map<double, std::list<int> >::iterator it_hash;

    double *tmp_x, *tmp_y;
    int count = 0;

    if(num == 0)
        return;

    tmp_x = new double[num];
    tmp_y = new double[num];

    for(int i = 0; i < num; i++) {
        it_hash = hash_table.find(x[i] * 1000.0 + y[i]);
        if(it_hash != hash_table.end()) {
            bool same = false;
            for(std::list<int>::iterator it_list = it_hash->second.begin(); it_list != it_hash->second.end(); it_list ++)
                if(x[*it_list] == x[i] && y[*it_list] == y[i]) {
                    same = true;
                    break;
                }
            if(same)
                continue;
            else {
                it_hash->second.push_back(i);
                tmp_x[count] = x[i];
                tmp_y[count++] = y[i];
            }
        }
        else {
            hash_table[x[i] * 1000.0 + y[i]].push_back(i);
            tmp_x[count] = x[i];
            tmp_y[count++] = y[i];
        }
    }

    delete x;
    delete y;
    x = tmp_x;
    y = tmp_y;
    num = count;

    return;
}


Edge *Delaunay_Voronoi::allocate_edge(Point *head, Point *tail)
{
    //Edge *new_edge = new Edge(head, tail);
    Edge *new_edge = edge_allocator.newElement(Edge(head, tail));
    edge_pool.push_back(new_edge);

    return new_edge;
}


Triangle *Delaunay_Voronoi::allocate_Triangle(Edge *edge1, Edge *edge2, Edge *edge3)
{
    //Triangle *new_triangle = new Triangle(edge1, edge2, edge3);
    
    //Triangle *new_triangle = new Triangle();
    Triangle *new_triangle = triangle_allocator.newElement(edge1, edge2, edge3);
    //Triangle *new_triangle = triangle_allocator.newElement(Triangle());
    //new_triangle->initialize_triangle_with_edges(edge1, edge2, edge3);
    triangle_pool.push_back(new_triangle);

    return new_triangle;
}


vector<Edge*> Delaunay_Voronoi::get_all_delaunay_edge()
{
    vector<Edge*> all_edges;

    for(unsigned i = 0; i < result_leaf_triangles.size(); i ++)
        if(result_leaf_triangles[i]->is_leaf) {
            all_edges.push_back(result_leaf_triangles[i]->edge[0]);
            all_edges.push_back(result_leaf_triangles[i]->edge[1]);
            all_edges.push_back(result_leaf_triangles[i]->edge[2]);
        }

    return all_edges;
}


vector<Edge*> Delaunay_Voronoi::get_all_legal_delaunay_edge()
{
    vector<Edge*> all_edges;

    for(unsigned i = 0; i < result_leaf_triangles.size(); i ++)
        if(result_leaf_triangles[i]->is_leaf)
            if(is_triangle_legal(result_leaf_triangles[i])) {
                all_edges.push_back(result_leaf_triangles[i]->edge[0]);
                all_edges.push_back(result_leaf_triangles[i]->edge[1]);
                all_edges.push_back(result_leaf_triangles[i]->edge[2]);
            }

    return all_edges;
}


bool Delaunay_Voronoi::check_if_all_outer_edge_out_of_region(double min_x, double max_x, double min_y, double max_y)
{
    for(unsigned i = 0; i < result_leaf_triangles.size(); i++) {
        if(!result_leaf_triangles[i]->is_leaf)
            continue;

        for(int j = 0; j < 3; j++)
            if(result_leaf_triangles[i]->edge[j]->twin_edge == NULL || result_leaf_triangles[i]->edge[j]->twin_edge->triangle == NULL)
                if((result_leaf_triangles[i]->edge[j]->head->x >= min_x && result_leaf_triangles[i]->edge[j]->head->x <= max_x &&
                    result_leaf_triangles[i]->edge[j]->head->y >= min_y && result_leaf_triangles[i]->edge[j]->head->y <= max_y) ||
                   (result_leaf_triangles[i]->edge[j]->tail->x >= min_x && result_leaf_triangles[i]->edge[j]->tail->x <= max_x &&
                    result_leaf_triangles[i]->edge[j]->tail->y >= min_y && result_leaf_triangles[i]->edge[j]->tail->y <= max_y)) {
                    printf("check_if_all_outer_edge_out_of_region: false (%lf, %lf) (%lf, %lf) (%lf, %lf)\n", result_leaf_triangles[i]->v[0]->x,
                                                                                                              result_leaf_triangles[i]->v[0]->y,
                                                                                                              result_leaf_triangles[i]->v[1]->x,
                                                                                                              result_leaf_triangles[i]->v[1]->y,
                                                                                                              result_leaf_triangles[i]->v[2]->x,
                                                                                                              result_leaf_triangles[i]->v[2]->y);
                    return false;
                }
    }
    //printf("check_if_all_outer_edge_out_of_region: true\n");
    return true;
}


static inline double point_distence_to_line(double px, double py, double x1, double y1, double x2, double y2)
{
    //TODO: may be slow
    if(x1 - x2 == 0)
        return fabs(px - x1);
    else if(y1 - y2 == 0)
        return fabs(py - y1);
    else {
        double A=(y1-y2)/(x1-x2);
        double B=y1-A*x1;
        return fabs((A*px+B-py)/sqrt(A*A+1));
    }
}


#define Point_distence_in_threshold(v, p1, p2, threshold) (point_distence_to_line(v->x, v->y, p1.x, p1.y, p2.x, p2.y) <= threshold)
#define All_distence_in_threshold(triangle, p1, p2, threshold) (        \
    Point_distence_in_threshold(triangle->v[0], p1, p2, threshold) &&   \
    Point_distence_in_threshold(triangle->v[1], p1, p2, threshold) &&   \
    Point_distence_in_threshold(triangle->v[2], p1, p2, threshold))

#define Vertexs_not_on_same_side(triangle, p1, p2) (                                             \
    triangle->v[0]->position_to_edge(&p1, &p2) * triangle->v[1]->position_to_edge(&p1, &p2) <= 0 || \
    triangle->v[1]->position_to_edge(&p1, &p2) * triangle->v[2]->position_to_edge(&p1, &p2) <= 0 )

#define Is_segment_in_triangle(triangle, p1, p2) \
    (p1.position_to_triangle(triangle) >= 0 && p2.position_to_triangle(triangle) >= 0)

#define Is_segment_intersected_with_edge(p_e1, p_e2, p1, p2) (                      \
    (p1.position_to_edge(p_e1, p_e2) * p2.position_to_edge(p_e1, p_e2) < 0   ) && \
    (p1.position_to_edge(p_e1, p_e2)!=0 || p2.position_to_edge(p_e1, p_e2)!=0) && \
    (p_e1->position_to_edge(&p1, &p2) * p_e2->position_to_edge(&p1, &p2) <= 0      ) )

#define Is_segment_intersected_with_one_of_edges(triangle, p1, p2) (            \
    Is_segment_intersected_with_edge(triangle->v[0], triangle->v[1], p1, p2) || \
    Is_segment_intersected_with_edge(triangle->v[1], triangle->v[2], p1, p2) || \
    Is_segment_intersected_with_edge(triangle->v[2], triangle->v[0], p1, p2) )

#define Is_triangle_intersecting_with_segment(triangle, p1, p2, threshold) (   \
    (threshold == 0 || All_distence_in_threshold(triangle, p1, p2, threshold)) && \
    Vertexs_not_on_same_side(triangle, p1, p2) &&       \
    (Is_segment_intersected_with_one_of_edges(triangle, p1, p2) || \
     Is_segment_in_triangle(triangle, p1, p2)) )

#define copy_points_value(dst, src) {   \
    dst->v[0]->x = src->v[0]->x;        \
    dst->v[0]->y = src->v[0]->y;        \
    dst->v[1]->x = src->v[1]->x;        \
    dst->v[1]->y = src->v[1]->y;        \
    dst->v[2]->x = src->v[2]->x;        \
    dst->v[2]->y = src->v[2]->y; }

std::vector<Triangle*> Delaunay_Voronoi::find_triangles_intersecting_with_segment(Point head, Point tail, double threshold)
{
    std::vector<Triangle*> triangles_found;
    Point uncyclic_head = head, uncyclic_tail = tail;
    Triangle *triangle_l, *triangle_r;

    if(uncyclic_head.x == uncyclic_tail.x) {
        while(uncyclic_head.x >= 360) uncyclic_head.x -= 360;
        while(uncyclic_head.x < 0) uncyclic_head.x += 360;
        while(uncyclic_tail.x >= 360) uncyclic_tail.x -= 360;
        while(uncyclic_tail.x < 0) uncyclic_tail.x += 360;
    }

    triangle_l = new Triangle(new Point(0, 0, -1), new Point(1, 0, -1), new Point(0, 1, -1));
    triangle_r = new Triangle(new Point(0, 0, -1), new Point(1, 0, -1), new Point(0, 1, -1));

    for(unsigned i = 0; i < result_leaf_triangles.size(); i++) {
        if(!result_leaf_triangles[i]->is_leaf)
            continue;

        if(result_leaf_triangles[i]->is_cyclic) {
            copy_points_value(triangle_l, result_leaf_triangles[i]);
            copy_points_value(triangle_r, result_leaf_triangles[i]);

            for(int j = 0; j < 3; j++) if(triangle_l->v[j]->x >= 180) triangle_l->v[j]->x -= 360;
            for(int j = 0; j < 3; j++) if(triangle_r->v[j]->x <  180) triangle_r->v[j]->x += 360;

            if(Is_triangle_intersecting_with_segment(triangle_l, uncyclic_head, uncyclic_tail, threshold) ||
               Is_triangle_intersecting_with_segment(triangle_r, uncyclic_head, uncyclic_tail, threshold) )
                triangles_found.push_back(result_leaf_triangles[i]);
        } else {
            if(Is_triangle_intersecting_with_segment(result_leaf_triangles[i], head, tail, threshold))
                triangles_found.push_back(result_leaf_triangles[i]);
        }
    }
    return triangles_found;
}


void Delaunay_Voronoi::remove_triangles_on_segment(Point head, Point tail)
{
    for(unsigned i = 0; i < result_leaf_triangles.size(); i++) {
        if(!result_leaf_triangles[i]->is_leaf)
            continue;

        if(result_leaf_triangles[i]->v[0]->position_to_edge(&head, &tail) * result_leaf_triangles[i]->v[1]->position_to_edge(&head, &tail) > 0 &&
           result_leaf_triangles[i]->v[1]->position_to_edge(&head, &tail) * result_leaf_triangles[i]->v[2]->position_to_edge(&head, &tail) > 0 )
            continue;

        /* two points of segment is in/on triangle */
        if(head.position_to_triangle(result_leaf_triangles[i]) >= 0 && tail.position_to_triangle(result_leaf_triangles[i]) >= 0) {
            remove_leaf_triangle(result_leaf_triangles[i]);
            continue;
        }

        /* segment is intersected with at least one edge of triangle */
        for(int j = 0; j < 3; j++)
            if((head.position_to_edge(result_leaf_triangles[i]->v[j], result_leaf_triangles[i]->v[(j+1)%3]) != 0 ||
               tail.position_to_edge(result_leaf_triangles[i]->v[j], result_leaf_triangles[i]->v[(j+1)%3]) != 0 ) &&
               head.position_to_edge(result_leaf_triangles[i]->v[j], result_leaf_triangles[i]->v[(j+1)%3]) *
               tail.position_to_edge(result_leaf_triangles[i]->v[j], result_leaf_triangles[i]->v[(j+1)%3]) < 0 &&
               result_leaf_triangles[i]->v[j]->position_to_edge(&head, &tail) *
               result_leaf_triangles[i]->v[(j+1)%3]->position_to_edge(&head, &tail) <= 0) {
                remove_leaf_triangle(result_leaf_triangles[i]);
                break;
            }
    }
}


std::vector<Triangle*> Delaunay_Voronoi::search_cyclic_triangles_for_rotated_grid(Point cyclic_boundary_head, Point cyclic_bounary_tail)
{
    return find_triangles_intersecting_with_segment(cyclic_boundary_head, cyclic_bounary_tail, 0);
}


inline double calculate_distence_square(double x1, double y1, double x2, double y2)
{
    return (x1-x2)*(x1-x2)+(y1-y2)*(y1-y2);
}


inline bool Delaunay_Voronoi::is_triangle_in_circle(Triangle* t, Point center, double radius)
{
    return calculate_distence_square(center.x,
                                     center.y,
                                     (t->v[0]->x + t->v[1]->x + t->v[2]->x) * 0.33333333333,
                                     (t->v[0]->y + t->v[1]->y + t->v[2]->y) * 0.33333333333) < radius*radius;
}


void Delaunay_Voronoi::remove_triangles_in_circle(Point center, double radius)
{
    for(unsigned i = 0; i < result_leaf_triangles.size(); i++)
        if(result_leaf_triangles[i]->is_leaf) {
            if (is_triangle_in_circle(result_leaf_triangles[i], center, radius))
                remove_leaf_triangle(result_leaf_triangles[i]);
        }
}


void Delaunay_Voronoi::correct_cyclic_triangles(std::vector<Triangle*> cyclic_triangles, bool is_grid_cyclic)
{
    if (is_grid_cyclic) {
        std::vector<Triangle*> right_triangles;
        std::vector<Point*> vl_pool, vr_pool;
        for (unsigned i = 0; i < cyclic_triangles.size(); i++) {
            Point *vl[3], *vr[3];
            Edge *el[3], *er[3];
            for (unsigned j = 0; j < 3; j++) {
                if (cyclic_triangles[i]->v[j]->x < 180.0) {
                    vl[j] = cyclic_triangles[i]->v[j];

                    /* search in the pool for the same point firstly */
                    vr[j] = NULL;
                    for (unsigned k = 0; k < vr_pool.size(); k++)
                        if (vr_pool[k]->id == vl[j]->id) {
                            vr[j] = vr_pool[k];
                            break;
                        }

                    if(vr[j] == NULL) {
                        vr[j] = new Point(vl[j]->x + 360.0, vl[j]->y, vl[j]->id);
                        vr_pool.push_back(vr[j]);
                    }
                }
                else {
                    vr[j] = cyclic_triangles[i]->v[j];

                    /* search in the pool for the same point firstly */
                    vl[j] = NULL;
                    for (unsigned k = 0; k < vl_pool.size(); k++)
                        if (vl_pool[k]->id == vr[j]->id) {
                            vl[j] = vl_pool[k];
                            break;
                        }

                    if(vl[j] == NULL) {
                        vl[j] = new Point(vr[j]->x - 360.0, vr[j]->y, vr[j]->id);
                        vl_pool.push_back(vl[j]);
                    }
                }
            }

            for (unsigned j = 0; j < 3; j++) {
                if (vl[j]->x >= 0.0 && vl[(j+1)%3]->x >= 0.0)
                    el[j] = cyclic_triangles[i]->edge[j];
                else
                    el[j] = allocate_edge(vl[j], vl[(j+1)%3]);

                if (vr[j]->x < 360.0 && vr[(j+1)%3]->x < 360.0)
                    er[j] = cyclic_triangles[i]->edge[j];
                else
                    er[j] = allocate_edge(vr[j], vr[(j+1)%3]);
                
            }
            
            /* Change cyclic triangles directly into left triangles */
            //int rank;
            //MPI_Comm_rank(process_thread_mgr->get_mpi_comm(), &rank);
            //printf("[%d] (%lf, %lf), (%lf, %lf), (%lf, %lf)\n", rank, el[0]->head->x, el[0]->head->y, el[1]->head->x, el[1]->head->y, el[2]->head->x, el[2]->head->y);
            cyclic_triangles[i]->initialize_triangle_with_edges(el[0], el[1], el[2]);

            /* Alloc new triangles for right triangles */
            Triangle *t = allocate_Triangle(er[0], er[1], er[2]);
            right_triangles.push_back(t);
            result_leaf_triangles.push_back(t);
        }

        check_and_set_twin_edge_relationship(&cyclic_triangles);
        check_and_set_twin_edge_relationship(&right_triangles);

    }
    else {
        for (unsigned i = 0; i < cyclic_triangles.size(); i++)
            remove_leaf_triangle(cyclic_triangles[i]);
    }
}


void Delaunay_Voronoi::recognize_cyclic_triangles()
{
    for(unsigned i = 0; i < result_leaf_triangles.size(); i++)
        if(result_leaf_triangles[i]->is_leaf) {
            if (result_leaf_triangles[i]->v[0]->calculate_distance(result_leaf_triangles[i]->v[1]) > 180 ||
                result_leaf_triangles[i]->v[1]->calculate_distance(result_leaf_triangles[i]->v[2]) > 180 ||
                result_leaf_triangles[i]->v[2]->calculate_distance(result_leaf_triangles[i]->v[0]) > 180 )
                result_leaf_triangles[i]->is_cyclic = true;
        }
}


inline void Delaunay_Voronoi::remove_leaf_triangle(Triangle* t)
{
    for (unsigned j = 0; j < 3; j++)
        if (t->edge[j]->twin_edge != NULL)
            t->edge[j]->twin_edge->twin_edge = NULL;
    t->is_leaf = false;
}

void Delaunay_Voronoi::relegalize_all_triangles()
{
    for(unsigned i = 0; i < result_leaf_triangles.size(); i++) {
        if (!result_leaf_triangles[i]->is_leaf)
            continue;

        //int rank;
        //MPI_Comm_rank(process_thread_mgr->get_mpi_comm(), &rank);
        
        //printf("legalizing %d: (%lf, %lf) (%lf, %lf) (%lf, %lf)\n", i, result_leaf_triangles[i]->v[0]->x, result_leaf_triangles[i]->v[0]->y, result_leaf_triangles[i]->v[1]->x,
                                                    //result_leaf_triangles[i]->v[1]->y, result_leaf_triangles[i]->v[2]->x, result_leaf_triangles[i]->v[2]->y);
        for(unsigned j = 0; j < 3; j++) {
            //printf("[%d] legalizing: (%lf, %lf) vs (%lf, %lf)--(%lf, %lf)\n", rank, result_leaf_triangles[i]->v[j]->x, result_leaf_triangles[i]->v[j]->y,
            //                                                             result_leaf_triangles[i]->edge[(j+1)%3]->head->x, result_leaf_triangles[i]->edge[(j+1)%3]->head->y,
            //                                                             result_leaf_triangles[i]->edge[(j+1)%3]->tail->x, result_leaf_triangles[i]->edge[(j+1)%3]->tail->y);
            //if(result_leaf_triangles[i]->edge[(j+1)%3]->twin_edge && result_leaf_triangles[i]->edge[(j+1)%3]->twin_edge->triangle)
            //    printf("twin: (%lf, %lf)-(%lf, %lf)-(%lf, %lf) is_leaf: %d\n", result_leaf_triangles[i]->edge[(j+1)%3]->twin_edge->triangle->v[0]->x,
            //                                                             result_leaf_triangles[i]->edge[(j+1)%3]->twin_edge->triangle->v[0]->y,
            //                                                             result_leaf_triangles[i]->edge[(j+1)%3]->twin_edge->triangle->v[1]->x,
            //                                                             result_leaf_triangles[i]->edge[(j+1)%3]->twin_edge->triangle->v[1]->y,
            //                                                             result_leaf_triangles[i]->edge[(j+1)%3]->twin_edge->triangle->v[2]->x,
            //                                                             result_leaf_triangles[i]->edge[(j+1)%3]->twin_edge->triangle->v[2]->y,
            //                                                             result_leaf_triangles[i]->edge[(j+1)%3]->twin_edge->triangle->is_leaf);
            relegalize_triangles(result_leaf_triangles[i]->v[j], result_leaf_triangles[i]->edge[(j+1)%3]/*, &result_leaf_triangles*/);
            //legalize_triangles(result_leaf_triangles[i]->v[j], result_leaf_triangles[i]->edge[(j+1)%3], &result_leaf_triangles);
        }
    }
}


void Delaunay_Voronoi::remove_triangles_on_or_out_of_boundary(double min_x, double max_x, double min_y, double max_y)
{
    for(unsigned i = 0; i < result_leaf_triangles.size(); i++)
        if(result_leaf_triangles[i]->is_leaf) {
            if (result_leaf_triangles[i]->v[0]->is_in_region(min_x, max_x, min_y, max_y) ||
                result_leaf_triangles[i]->v[1]->is_in_region(min_x, max_x, min_y, max_y) ||
                result_leaf_triangles[i]->v[2]->is_in_region(min_x, max_x, min_y, max_y))
                continue;

            remove_leaf_triangle(result_leaf_triangles[i]);
        }
}


void Delaunay_Voronoi::get_triangles_intersecting_with_segment(Point head, Point tail, Triangle_Transport *output_triangles, int *num_triangles, int buf_len, double threshold)
{
    int current = 0;

    std::vector<Triangle*> ts = find_triangles_intersecting_with_segment(head, tail, threshold);

    for(unsigned i = 0; i < ts.size(); i++)
        output_triangles[current++] = Triangle_Transport(Point(ts[i]->v[0]->x, ts[i]->v[0]->y, global_index[ts[i]->v[0]->id]),  // TODO: use global_index as id directly
                                                         Point(ts[i]->v[1]->x, ts[i]->v[1]->y, global_index[ts[i]->v[1]->id]), 
                                                         Point(ts[i]->v[2]->x, ts[i]->v[2]->y, global_index[ts[i]->v[2]->id]));

    assert(current < buf_len);
    *num_triangles = current;
}


static inline unsigned hash_triangle_by_id(Triangle_Transport triangle)
{
    return triangle.v[0].id ^ (triangle.v[1].id << 10) ^ (triangle.v[2].id << 20);
}

unsigned Delaunay_Voronoi::calculate_triangles_checksum(Triangle_Transport *triangles, int num_triangles)
{
    if(num_triangles < 1)
        return PDLN_CHECKSUM_FALSE;

    sort_points_in_triangle(triangles, num_triangles);
    //sort_triangles(triangles, num_triangles);

    int checksum = hash_triangle_by_id(triangles[0]);
    for(int i = 1; i < num_triangles; i ++)
        checksum = checksum ^ hash_triangle_by_id(triangles[i]);

    return checksum;
}


unsigned Delaunay_Voronoi::calculate_triangles_intersected_checksum(Point head, Point tail, double threshold)
{
    /* Let n be the number of points, if there are b vertices on the convex hull,
     * then any triangulation of the points has at most 2n − 2 − b triangles,
     * plus one exterior face */
    Triangle_Transport *triangles = new Triangle_Transport[2*num_cells];
    int num_triangles;
    unsigned checksum;

    get_triangles_intersecting_with_segment(head, tail, triangles, &num_triangles, 2*num_cells, threshold);

    checksum = calculate_triangles_checksum(triangles, num_triangles);

    //char filename[64];
    //int rank;
    //MPI_Comm_rank(process_thread_mgr->get_mpi_comm(), &rank);
    //snprintf(filename, 64, "log/boundary_triangles_%d_%x", rank, checksum);
    //plot_triangles_into_file(filename, triangles, num_triangles);


    delete[] triangles;
    return checksum;
}


/* including triangles intersecting with boundary */
void Delaunay_Voronoi::get_triangles_in_region(double min_x, double max_x, double min_y, double max_y, 
                                               Triangle_Transport *output_triangles, int *output_num_triangles, int buf_len)
{
    int current = 0;
    int num_triangles = 0;

    //printf("result_leaf_triangles.size: %lu\n", result_leaf_triangles.size());
    for(unsigned i = 0; i < result_leaf_triangles.size(); i++) {
        if(!result_leaf_triangles[i]->is_leaf)
            continue;

        bool in = true;
        for(int j = 0; j < 3; j++)
            if(!(result_leaf_triangles[i]->v[j]->x < max_x && result_leaf_triangles[i]->v[j]->x > min_x &&
                 result_leaf_triangles[i]->v[j]->y < max_y && result_leaf_triangles[i]->v[j]->y > min_y)) {
                in = false;
                break;
            }

        if(in) {
            output_triangles[current++] = Triangle_Transport(Point(result_leaf_triangles[i]->v[0]->x, result_leaf_triangles[i]->v[0]->y, global_index[result_leaf_triangles[i]->v[0]->id]),
                                                             Point(result_leaf_triangles[i]->v[1]->x, result_leaf_triangles[i]->v[1]->y, global_index[result_leaf_triangles[i]->v[1]->id]),
                                                             Point(result_leaf_triangles[i]->v[2]->x, result_leaf_triangles[i]->v[2]->y, global_index[result_leaf_triangles[i]->v[2]->id]));
            //printf("inner: %d, %d, %d\n", result_leaf_triangles[i]->v[0]->id, result_leaf_triangles[i]->v[1]->id, result_leaf_triangles[i]->v[2]->id);
        }
    }
    assert(current < buf_len);
    //printf("kernel: %d, buf_len: %d\n", current, buf_len);

    get_triangles_intersecting_with_segment(Point(min_x, min_y), Point(max_x, min_y), output_triangles+current, &num_triangles, buf_len-current);
    current += num_triangles;
    get_triangles_intersecting_with_segment(Point(max_x, min_y), Point(max_x, max_y), output_triangles+current, &num_triangles, buf_len-current);
    current += num_triangles;
    get_triangles_intersecting_with_segment(Point(max_x, max_y), Point(min_x, max_y), output_triangles+current, &num_triangles, buf_len-current);
    current += num_triangles;
    get_triangles_intersecting_with_segment(Point(min_x, max_y), Point(min_x, min_y), output_triangles+current, &num_triangles, buf_len-current);
    *output_num_triangles = current + num_triangles;
}


void Delaunay_Voronoi::update_virtual_polar_info()
{
    for(unsigned i = 0; i < result_leaf_triangles.size(); i++) {
        if(!result_leaf_triangles[i]->is_leaf)
            continue;
        if(result_leaf_triangles[i]->contain_vertex(cells[vpolar_local_index].center))
            triangles_containing_vpolar.push_back(result_leaf_triangles[i]);
    }
}


void Delaunay_Voronoi::remove_triangles_only_containing_virtual_polar()
{
    double common_lat = cells[vpolar_local_index].center->y;

    for(unsigned i = 0; i < result_leaf_triangles.size(); i++) {
        if(!result_leaf_triangles[i]->is_leaf)
            continue;
        if(result_leaf_triangles[i]->v[0]->y == common_lat &&
           result_leaf_triangles[i]->v[1]->y == common_lat &&
           result_leaf_triangles[i]->v[2]->y == common_lat)
            remove_leaf_triangle(result_leaf_triangles[i]);
    }
}


void Delaunay_Voronoi::remove_triangles_till(int max_global_index)
{
    for(unsigned i = 0; i < result_leaf_triangles.size(); i++) {
        if(!result_leaf_triangles[i]->is_leaf)
            continue;
        if((global_index[result_leaf_triangles[i]->v[0]->id] < max_global_index && global_index[result_leaf_triangles[i]->v[0]->id] >= 0) ||
           (global_index[result_leaf_triangles[i]->v[1]->id] < max_global_index && global_index[result_leaf_triangles[i]->v[1]->id] >= 0) ||
           (global_index[result_leaf_triangles[i]->v[2]->id] < max_global_index && global_index[result_leaf_triangles[i]->v[2]->id] >= 0))
            remove_leaf_triangle(result_leaf_triangles[i]);
    }
}

#ifdef OPENCV
void Delaunay_Voronoi::plot_into_file(const char *filename, double min_x, double max_x, double min_y, double max_y)
{
    unsigned num_edges;
    double *head_coord[2], *tail_coord[2];

    num_edges = 3 * result_leaf_triangles.size();
    head_coord[0] = new double[num_edges];
    head_coord[1] = new double[num_edges];
    tail_coord[0] = new double[num_edges];
    tail_coord[1] = new double[num_edges];

    num_edges = 0;
    for(unsigned i = 0; i < result_leaf_triangles.size(); i ++)
        if(result_leaf_triangles[i]->is_leaf) {
            for(unsigned j = 0; j < 3; j++) {
                head_coord[0][num_edges] = result_leaf_triangles[i]->edge[j]->head->x;
                head_coord[1][num_edges] = result_leaf_triangles[i]->edge[j]->head->y;
                tail_coord[0][num_edges] = result_leaf_triangles[i]->edge[j]->tail->x;
                tail_coord[1][num_edges++] = result_leaf_triangles[i]->edge[j]->tail->y;
            }
            if(result_leaf_triangles[i]->is_cyclic)
                for(unsigned j = num_edges-1; j > num_edges-4; j--) {
                    if(head_coord[0][j] > 180) head_coord[0][j] -= 360;
                    if(tail_coord[0][j] > 180) tail_coord[0][j] -= 360;
                }
        }

    assert(num_edges%3 == 0);
    assert(num_edges <= 3 * result_leaf_triangles.size());
    plot_edge_into_file(filename, head_coord, tail_coord, num_edges, min_x, max_x, min_y, max_y);

    delete head_coord[0];
    delete head_coord[1];
    delete tail_coord[0];
    delete tail_coord[1];
}


void Delaunay_Voronoi::plot_current_step_into_file(const char *filename)
{
    unsigned num_edges;
    double *head_coord[2], *tail_coord[2];

    num_edges = 3 * triangle_pool.size();
    head_coord[0] = new double[num_edges];
    head_coord[1] = new double[num_edges];
    tail_coord[0] = new double[num_edges];
    tail_coord[1] = new double[num_edges];

    num_edges = 0;
    for(unsigned i = 0; i < triangle_pool.size(); i ++)
        if(triangle_pool[i]->is_leaf)
            for(unsigned j = 0; j < 3; j++) {
                head_coord[0][num_edges] = triangle_pool[i]->edge[j]->head->x;
                head_coord[1][num_edges] = triangle_pool[i]->edge[j]->head->y;
                tail_coord[0][num_edges] = triangle_pool[i]->edge[j]->tail->x;
                tail_coord[1][num_edges++] = triangle_pool[i]->edge[j]->tail->y;
            }

    assert(num_edges%3 == 0);
    assert(num_edges <= 3 * triangle_pool.size());
    plot_projected_edge_into_file(filename, head_coord, tail_coord, num_edges, PDLN_PLOT_COLOR_DEFAULT, PDLN_PLOT_FILEMODE_NEW);

    delete head_coord[0];
    delete head_coord[1];
    delete tail_coord[0];
    delete tail_coord[1];
}

void Delaunay_Voronoi::plot_projection_into_file(const char *filename, double min_x, double max_x, double min_y, double max_y)
{
    unsigned num_edges;
    double *head_coord[2], *tail_coord[2];

    num_edges = 3 * result_leaf_triangles.size();
    head_coord[0] = new double[num_edges];
    head_coord[1] = new double[num_edges];
    tail_coord[0] = new double[num_edges];
    tail_coord[1] = new double[num_edges];

    num_edges = 0;
    for(unsigned i = 0; i < result_leaf_triangles.size(); i ++)
        if(result_leaf_triangles[i]->is_leaf)
            for(unsigned j = 0; j < 3; j++) {
                head_coord[0][num_edges] = result_leaf_triangles[i]->edge[j]->head->x;
                head_coord[1][num_edges] = result_leaf_triangles[i]->edge[j]->head->y;
                tail_coord[0][num_edges] = result_leaf_triangles[i]->edge[j]->tail->x;
                tail_coord[1][num_edges++] = result_leaf_triangles[i]->edge[j]->tail->y;
            }

    assert(num_edges <= 3 * result_leaf_triangles.size());
    plot_projected_edge_into_file(filename, head_coord, tail_coord, num_edges, PDLN_PLOT_COLOR_DEFAULT, PDLN_PLOT_FILEMODE_NEW);

    num_edges = 0;
    for(unsigned i = 0; i < triangles_containing_vpolar.size(); i++)
        for(unsigned j = 0; j < 3; j++) {
            head_coord[0][num_edges] = triangles_containing_vpolar[i]->edge[j]->head->x;
            head_coord[1][num_edges] = triangles_containing_vpolar[i]->edge[j]->head->y;
            tail_coord[0][num_edges] = triangles_containing_vpolar[i]->edge[j]->tail->x;
            tail_coord[1][num_edges++] = triangles_containing_vpolar[i]->edge[j]->tail->y;
        }
    plot_projected_edge_into_file(filename, head_coord, tail_coord, num_edges, PDLN_PLOT_COLOR_RED, PDLN_PLOT_FILEMODE_APPEND);

    /* Debug staff */
    //num_edges = 0;
    //for(unsigned i = 0; i < 3; i++) {
    //    head_coord[0][num_edges] = global_p_lon[i];
    //    head_coord[1][num_edges] = global_p_lat[i];
    //    tail_coord[0][num_edges] = global_p_lon[(i+1)%3];
    //    tail_coord[1][num_edges++] = global_p_lat[(i+1)%3];
    //}
    //plot_projected_edge_into_file(filename, head_coord, tail_coord, num_edges, PDLN_PLOT_COLOR_RED, PDLN_PLOT_FILEMODE_APPEND);

    delete head_coord[0];
    delete head_coord[1];
    delete tail_coord[0];
    delete tail_coord[1];
}


void Delaunay_Voronoi::plot_original_points_into_file(const char *filename, double min_x, double max_x, double min_y, double max_y)
{
    double *coord[2];

    coord[0] = new double[num_cells];
    coord[1] = new double[num_cells];

    int num = 0;
    for(int i = 0; i < num_cells; i ++) {
        coord[0][num] = cells[i].center->x;
        coord[1][num++] = cells[i].center->y;
    }

    assert(num == num_cells);
    plot_points_into_file(filename, coord[0], coord[1], num_cells, min_x, max_x, min_y, max_y);

    delete coord[0];
    delete coord[1];
}


void plot_triangles_into_file(const char *prefix, Triangle_Transport *t, int num, bool plot_cyclic_triangles)
{
    int num_edges;
    double *head_coord[2], *tail_coord[2];
    char filename[128];

    num_edges = 3 * num;
    head_coord[0] = new double[num_edges];
    head_coord[1] = new double[num_edges];
    tail_coord[0] = new double[num_edges];
    tail_coord[1] = new double[num_edges];

    num_edges = 0;
    for(int i = 0; i < num; i ++)
        if (t[i].v[0].calculate_distance(&t[i].v[1]) < 180 && t[i].v[1].calculate_distance(&t[i].v[2]) < 180 && t[i].v[2].calculate_distance(&t[i].v[0]) < 180)
            for(int j = 0; j < 3; j++) {
                head_coord[0][num_edges] = t[i].v[j].x;
                head_coord[1][num_edges] = t[i].v[j].y;
                tail_coord[0][num_edges] = t[i].v[(j+1)%3].x;
                tail_coord[1][num_edges++] = t[i].v[(j+1)%3].y;
            }
        else if(plot_cyclic_triangles) {
            for(int j = 0; j < 3; j++)
                if(t[i].v[j].x > 180) t[i].v[j].x -= 360;
            for(int j = 0; j < 3; j++) {
                head_coord[0][num_edges] = t[i].v[j].x;
                head_coord[1][num_edges] = t[i].v[j].y;
                tail_coord[0][num_edges] = t[i].v[(j+1)%3].x;
                tail_coord[1][num_edges++] = t[i].v[(j+1)%3].y;
            }
        }

    assert(num_edges%3 == 0);
    snprintf(filename, 128, "%s.png", prefix);
    plot_edge_into_file(filename, head_coord, tail_coord, num_edges);

    delete head_coord[0];
    delete head_coord[1];
    delete tail_coord[0];
    delete tail_coord[1];

}


void plot_triangles_into_file(const char *prefix, std::vector<Triangle*> t)
{
    unsigned num = t.size();
    int num_edges;
    double *head_coord[2], *tail_coord[2];
    char filename[128];

    num_edges = 3 * num;
    head_coord[0] = new double[num_edges];
    head_coord[1] = new double[num_edges];
    tail_coord[0] = new double[num_edges];
    tail_coord[1] = new double[num_edges];

    num_edges = 0;
    for(unsigned i = 0; i < num; i ++)
        for(int j = 0; j < 3; j++) {
            head_coord[0][num_edges] = t[i]->v[j]->x;
            head_coord[1][num_edges] = t[i]->v[j]->y;
            tail_coord[0][num_edges] = t[i]->v[(j+1)%3]->x;
            tail_coord[1][num_edges++] = t[i]->v[(j+1)%3]->y;
        }

    assert(num_edges%3 == 0);
    snprintf(filename, 128, "%s.png", prefix);
    //plot_edge_into_file(filename, head_coord, tail_coord, num_edges);
    plot_projected_edge_into_file(filename, head_coord, tail_coord, num_edges, PDLN_PLOT_COLOR_WHITE, PDLN_PLOT_FILEMODE_NEW);

    delete head_coord[0];
    delete head_coord[1];
    delete tail_coord[0];
    delete tail_coord[1];

}
#endif


static int compare_cell_x(const void* a, const void* b)
{
    Cell t1 = *(const Cell*)a;
    Cell t2 = *(const Cell*)b;

    if(t1.center->x < t2.center->x) return -1;
    if(t1.center->x > t2.center->x) return  1;
    return 0;
}


static int compare_cell_y(const void* a, const void* b)
{
    Cell t1 = *(const Cell*)a;
    Cell t2 = *(const Cell*)b;

    if(t1.center->y < t2.center->y) return -1;
    if(t1.center->y > t2.center->y) return  1;
    return 0;
}


void sort_cells(Cell *cells, int num_cells)
{
    merge_sort(cells, num_cells, sizeof(Cell), compare_cell_x);
    merge_sort(cells, num_cells, sizeof(Cell), compare_cell_y);
}


void Delaunay_Voronoi::save_original_points_into_file()
{
    Cell *tmp_cells;

#ifdef CHECK_PARALLEL_CONSISTENCY
    tmp_cells = new Cell[num_cells];
    memcpy(tmp_cells, cells, num_cells*sizeof(Cell));
    sort_cells(tmp_cells, num_cells);
#else
    tmp_cells = cells;
#endif

    FILE *fp;
    fp = fopen("log/original_points.txt", "w");
    for(int i = 0; i < num_cells; i++)
        if(tmp_cells[i].center->x < 5 && tmp_cells[i].center->x > -5 && tmp_cells[i].center->y < 5 && tmp_cells[i].center->y > -5)
            fprintf(fp, "%.20lf, %.20lf\n", tmp_cells[i].center->x, tmp_cells[i].center->y);
    fclose(fp);

#ifdef CHECK_PARALLEL_CONSISTENCY
    delete [] tmp_cells;
#endif
}


void Delaunay_Voronoi::update_all_points_coord(double *x_values, double *y_values, int num)
{
    assert(num == num_cells);
    for(int i = 0; i < num; i++) {
        cells[i].center->x = x_values[i];
        cells[i].center->y = y_values[i];
    }

    for(unsigned i = 0; i < result_leaf_triangles.size(); i++)
        if(result_leaf_triangles[i]->is_leaf)
            result_leaf_triangles[i]->calulate_circum_circle();
}


void Delaunay_Voronoi::update_points_coord_y(double reset_lat_value, vector<int> *polars_local_index)
{
    for(unsigned i = 0; i < polars_local_index->size(); i++)
        cells[(*polars_local_index)[i]].center->y = reset_lat_value;
    // WARNING: circum_circle will not been recalculated
}


void Delaunay_Voronoi::uncyclic_all_points()
{
    for(unsigned i = 0; i < num_cells; i++) {
        while(cells[i].center->x >= 360)
            cells[i].center->x -= 360;
        while(cells[i].center->x < 0)
            cells[i].center->x += 360;
    }
}

Triangle_Transport::Triangle_Transport(Point p0, Point p1, Point p2)
{
    v[0] = p0;
    v[1] = p1;
    v[2] = p2;
}


bool operator == (Triangle_Transport t1, Triangle_Transport t2)
{
#ifdef DEBUG
    assert(t1.v[0] != t1.v[1] && t1.v[1] != t1.v[2] && t1.v[2] != t1.v[0]);
    assert(t2.v[0] != t2.v[1] && t2.v[1] != t2.v[2] && t2.v[2] != t2.v[0]);
#endif
    if(t2.v[0] != t1.v[0] && t2.v[0] != t1.v[1] && t2.v[0] != t1.v[2])
        return false;
    if(t2.v[1] != t1.v[0] && t2.v[1] != t1.v[1] && t2.v[1] != t1.v[2])
        return false;
    if(t2.v[2] != t1.v[0] && t2.v[2] != t1.v[1] && t2.v[2] != t1.v[2])
        return false;
    return true;
}


void save_triangles_info_file(const char *prefix, Triangle_Transport *t, int num)
{
    char filename[128];
    FILE* fp;

    snprintf(filename, 128, "%s.txt", prefix);
    fp = fopen(filename, "w");
    for(int i = 0; i < num; i++)
        fprintf(fp, "[%d] (%.15lf, %.15lf), (%.15lf, %.15lf), (%.15lf, %.15lf)\n", i, t[i].v[0].x, t[i].v[0].y, t[i].v[1].x, t[i].v[1].y, t[i].v[2].x, t[i].v[2].y);
    fclose(fp);
}
