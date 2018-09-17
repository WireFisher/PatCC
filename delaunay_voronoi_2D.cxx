#include "mpi.h"
#include "delaunay_voronoi_2D.h"
#include "opencv_utils.h"
#include "pd_assert.h"
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
    Triangle_pack t1 = *(const Triangle_pack*)a;
    Triangle_pack t2 = *(const Triangle_pack*)b;

    if(t1.v[2].id < t2.v[2].id) return -1;
    if(t1.v[2].id > t2.v[2].id) return  1;
    return 0;
}


static int compare_v1(const void* a, const void* b)
{
    Triangle_pack t1 = *(const Triangle_pack*)a;
    Triangle_pack t2 = *(const Triangle_pack*)b;

    if(t1.v[1].id < t2.v[1].id) return -1;
    if(t1.v[1].id > t2.v[1].id) return  1;
    return 0;
}


static int compare_v0(const void* a, const void* b)
{
    Triangle_pack t1 = *(const Triangle_pack*)a;
    Triangle_pack t2 = *(const Triangle_pack*)b;

    if(t1.v[0].id < t2.v[0].id) return -1;
    if(t1.v[0].id > t2.v[0].id) return  1;
    return 0;
}


static int compare_lon(const void* a, const void* b)
{
    Triangle_pack t1 = *(const Triangle_pack*)a;
    Triangle_pack t2 = *(const Triangle_pack*)b;

    if(t1.v[0].x < t2.v[0].x) return -1;
    if(t1.v[0].x > t2.v[0].x) return  1;
    return 0;
}


static inline void radix_sort(Triangle_pack *triangles, int num_triangles)
{
    PDASSERT(sizeof(Triangle_pack) > sizeof(void *)/2);
    merge_sort(triangles, num_triangles, sizeof(Triangle_pack), compare_lon);
    merge_sort(triangles, num_triangles, sizeof(Triangle_pack), compare_v2);
    merge_sort(triangles, num_triangles, sizeof(Triangle_pack), compare_v1);
    merge_sort(triangles, num_triangles, sizeof(Triangle_pack), compare_v0);
}


void sort_points_in_triangle(Triangle_pack *triangles, int num_triangles)
{
    for(int i = 0; i < num_triangles; i++) {
        if(triangles[i].v[0].id > triangles[i].v[1].id) swap(&triangles[i].v[0], &triangles[i].v[1]);
        if(triangles[i].v[1].id > triangles[i].v[2].id) swap(&triangles[i].v[1], &triangles[i].v[2]);
        if(triangles[i].v[0].id > triangles[i].v[1].id) swap(&triangles[i].v[0], &triangles[i].v[1]);
    }
}


inline void sort_points_in_triangle(Triangle_pack& triangle)
{
    if(triangle.v[0].id > triangle.v[1].id) std::swap(triangle.v[0], triangle.v[1]);
    if(triangle.v[1].id > triangle.v[2].id) std::swap(triangle.v[1], triangle.v[2]);
    if(triangle.v[0].id > triangle.v[1].id) std::swap(triangle.v[0], triangle.v[1]);
}


void sort_triangles(Triangle_pack *triangles, int num_triangles)
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
    PDASSERT(!(on1 && on2));
    PDASSERT(!(on2 && on3));
    PDASSERT(!(on3 && on1));
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


int Point::position_to_triangle(const Triangle_pack *triangle) const
{
#ifdef DDEBUG
    bool on1 = position_to_edge(&triangle->v[0], &triangle->v[1]) == 0;
    bool on2 = position_to_edge(&triangle->v[1], &triangle->v[2]) == 0;
    bool on3 = position_to_edge(&triangle->v[2], &triangle->v[0]) == 0;
    PDASSERT(!(on1 && on2));
    PDASSERT(!(on2 && on3));
    PDASSERT(!(on3 && on1));
#endif
    int ret = 0;
    int pos = position_to_edge(&triangle->v[0], &triangle->v[1]);
    if (pos == -1)
        return -1;
    else if (pos == 0)
        ret = 1;
    pos = position_to_edge(&triangle->v[1], &triangle->v[2]);
    if (pos == -1)
        return -1;
    else if (pos == 0)
        ret = 2;
    pos = position_to_edge(&triangle->v[2], &triangle->v[0]);
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
    PDASSERT(center != p1 && center != p2) ;

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


void Delaunay_Voronoi::legalize_triangles(Point *vr, Edge *edge, unsigned *sp)
{
    if (is_triangle_legal(vr, edge))
        return;

    //EXECUTION_REPORT(REPORT_ERROR, -1, edge->triangle->is_leaf, "remap software error1 in legalize_triangles\n");
    PDASSERT(edge->triangle->is_leaf);
    //EXECUTION_REPORT(REPORT_ERROR, -1, edge->twin_edge->triangle->is_leaf, "remap software error2 in legalize_triangles %lx\n", (long)(edge->twin_edge->triangle));
    PDASSERT(edge->twin_edge->triangle->is_leaf);
    push(sp, edge->twin_edge->triangle);
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
    push(sp, tikr);
    push(sp, tjrk);
    legalize_triangles(vr, eik, sp);
    legalize_triangles(vr, ekj, sp);
}


void Delaunay_Voronoi::relegalize_triangles(Point *vr, Edge *edge)
{
    if (!is_triangle_ambiguous(vr, edge))
        return;

    if (edge->triangle->is_cyclic != edge->twin_edge->triangle->is_cyclic)
        return;
    //EXECUTION_REPORT(REPORT_ERROR, -1, edge->triangle->is_leaf, "remap software error1 in legalize_triangles\n");
    PDASSERT(edge->triangle->is_leaf);
    //EXECUTION_REPORT(REPORT_ERROR, -1, edge->twin_edge->triangle->is_leaf, "remap software error2 in legalize_triangles %lx\n", (long)(edge->twin_edge->triangle));
    PDASSERT(edge->twin_edge->triangle->is_leaf);

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
                       (v[0]->x * (v[2]->y - v[1]->y) + v[1]->x * (v[0]->y - v[2]->y) + v[2]->x * (v[1]->y - v[0]->y)) * 0.5;
    circum_center[1] = (ab * (v[2]->x - v[1]->x) + cd * (v[0]->x - v[2]->x) + ef * (v[1]->x - v[0]->x)) /
                       (v[0]->y * (v[2]->x - v[1]->x) + v[1]->y * (v[0]->x - v[2]->x) + v[2]->y * (v[1]->x - v[0]->x)) * 0.5;
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
        /* if there are unmarked redundant points, the PDASSERTion may fail */
        PDASSERT(std::fabs(det(pt1, pt2, pt3)) > FLOAT_ERROR_HI && std::fabs(det(pt2, pt3, pt1)) > FLOAT_ERROR_HI && std::fabs(det(pt3, pt1, pt2)) > FLOAT_ERROR_HI);
        //EXECUTION_REPORT(REPORT_ERROR, -1, edge1->tail==edge2->head && edge2->tail==edge3->head && edge3->tail==edge1->head, "edges given to construct triangle is invalid.");
#ifdef DEBUG
        if(!(edge1->tail==edge2->head && edge2->tail==edge3->head && edge3->tail==edge1->head)) {
            printf("edge1: %p->%p, edge2: %p->%p, edge3: %p->%p\n", edge1->head, edge1->tail, edge2->head, edge2->tail, edge3->head, edge3->tail);
            printf("1: (%lf, %lf)->(%lf, %lf) 2: (%lf, %lf)->(%lf, %lf) 3:(%lf, %lf)->(%lf, %lf)\n", edge1->head->x, edge1->head->y, edge1->tail->x, edge1->tail->y,
                                                                                                     edge2->head->x, edge2->head->y, edge2->tail->x, edge2->tail->y,
                                                                                                     edge3->head->x, edge3->head->y, edge3->tail->x, edge3->tail->y);
        }
#endif
        PDASSERT(edge1->tail==edge2->head && edge2->tail==edge3->head && edge3->tail==edge1->head);
           
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
            //printf("[%d] PDASSERT false\n", rank);
            PDASSERT(false);
            v[1] = pt3;
            v[2] = pt2;
            this->edge[0] = edge3->twin_edge;
            this->edge[1] = edge2->twin_edge;
            this->edge[2] = edge1->twin_edge;
            //EXECUTION_REPORT(REPORT_ERROR, -1, edge3->twin_edge != NULL && edge2->twin_edge != NULL && edge1->twin_edge != NULL, "remap software error3 in new Triangle");
            PDASSERT(edge3->twin_edge != NULL && edge2->twin_edge != NULL && edge1->twin_edge != NULL);
        }
    }
    remained_points_head = -1;
    remained_points_tail = -1;
    
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
                       (pt[0]->x * (pt[2]->y - pt[1]->y) + pt[1]->x * (pt[0]->y - pt[2]->y) + pt[2]->x * (pt[1]->y - pt[0]->y)) * 0.5;
    center[1] = (ab * (pt[2]->x - pt[1]->x) + cd * (pt[0]->x - pt[2]->x) + ef * (pt[1]->x - pt[0]->x)) /
                       (pt[0]->y * (pt[2]->x - pt[1]->x) + pt[1]->y * (pt[0]->x - pt[2]->x) + pt[2]->y * (pt[1]->x - pt[0]->x)) * 0.5;
    double radius2 = ((pt[0]->x - center[0]) * (pt[0]->x - center[0])) + ((pt[0]->y - center[1]) * (pt[0]->y - center[1]));
    double dist2 = ((pt[3]->x - center[0]) * (pt[3]->x - center[0])) + ((pt[3]->y - center[1]) * (pt[3]->y - center[1]));
    return std::fabs(dist2 - radius2) < tolerance;
}


int Triangle::find_best_candidate_point(Point* buf) const
{
    double min_dist=1e10, dist;
    PDASSERT(remained_points_head != -1 || remained_points_tail != -1);

    double center_x = (v[0]->x+v[1]->x+v[2]->x) * 0.3333333333333333333333333333;
    double center_y = (v[0]->y+v[1]->y+v[2]->y) * 0.3333333333333333333333333333;
    int best_candidate_id = remained_points_head;
    for (int i = remained_points_head; i > -1; i = buf[i].next) {
        dist = buf[i].calculate_distance(center_x, center_y);
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


void Delaunay_Voronoi::swap_points(int idx1, int idx2)
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


void Delaunay_Voronoi::distribute_points_into_triangles(int head, int tail, unsigned base, unsigned top)
{
    int start = head;
    if (tail == -1)
        return;
    PDASSERT(head != -1);

    for (unsigned i = base+1; i <= top; i++) {
        if (!triangle_stack[i]->is_leaf)
            continue;

        int end = tail;
        int j = start;
        for (; j != end && j > -1;) {
            if (all_points[j].position_to_triangle(triangle_stack[i]) >= 0) {
                j = all_points[j].next;
            } else {
                swap_points(j, end);
                end = all_points[end].prev;
            }
        }
        if (j > -1 && all_points[end].position_to_triangle(triangle_stack[i]) < 0) // Case "j == end"
            end = all_points[end].prev;
        triangle_stack[i]->set_remained_points(start, end);
        if (end > -1) {
            start = all_points[end].next;
            all_points[end].next = -1;
            if (start > -1)
                all_points[start].prev = -1;
            else
                break;
        }
    }
    PDASSERT(start == -1);
}


void Triangle::set_remained_points(int head, int tail)
{
    if (tail > -1 && head > -1) {
        remained_points_head = head;
        remained_points_tail = tail;
    } else {
        remained_points_head = -1;
        remained_points_tail = -1;
    }
}


Point* Triangle::pop_tail(Point* buf)
{
    PDASSERT(remained_points_tail != -1);
    Point* old_tail = &buf[remained_points_tail];

    remained_points_tail = buf[remained_points_tail].prev;
    if (remained_points_tail > -1)
        buf[remained_points_tail].next = -1;
    else
        remained_points_head = -1;

    return old_tail;
}


void Delaunay_Voronoi::push(unsigned *stack_top, Triangle* t)
{
    *stack_top = *stack_top + 1;
    if (*stack_top >= stack_size) {
        stack_size *= 2;
        Triangle** tmp = new Triangle*[stack_size];
        memcpy(tmp, triangle_stack, sizeof(Triangle*) * *stack_top);
        delete[] triangle_stack;
        triangle_stack = tmp;
    }

    triangle_stack[*stack_top] = t;
}


void Delaunay_Voronoi::triangularization_process(Triangle *triangle, unsigned stack_base)
{
    unsigned stack_top = stack_base;

#ifdef DEBUG
    PDASSERT(triangle->is_leaf);
#endif
    if (triangle->remained_points_tail == -1) {
        result_leaf_triangles.push_back(triangle);
        return;
    }

    triangle->is_leaf = false;

    int candidate_id = triangle->find_best_candidate_point(all_points);
    swap_points(candidate_id, triangle->remained_points_tail);
    Point* dividing_point = triangle->pop_tail(all_points);

    if (dividing_point->position_to_triangle(triangle) == 0) { // inside
        Edge *e_v1_can = allocate_edge(triangle->v[0], dividing_point);
        Edge *e_can_v1 = generate_twins_edge(e_v1_can);
        Edge *e_v2_can = allocate_edge(triangle->v[1], dividing_point);
        Edge *e_can_v2 = generate_twins_edge(e_v2_can);
        Edge *e_v3_can = allocate_edge(triangle->v[2], dividing_point);
        Edge *e_can_v3 = generate_twins_edge(e_v3_can);
        Triangle *t_can_v1_v2 = allocate_Triangle(e_can_v1, triangle->edge[0], e_v2_can);
        Triangle *t_can_v2_v3 = allocate_Triangle(e_can_v2, triangle->edge[1], e_v3_can);
        Triangle *t_can_v3_v1 = allocate_Triangle(e_can_v3, triangle->edge[2], e_v1_can);
        push(&stack_top, triangle);
        push(&stack_top, t_can_v1_v2);
        push(&stack_top, t_can_v2_v3);
        push(&stack_top, t_can_v3_v1);
        
        /* Actually, t_can_v1_v2->v[0] is dividing_point */
        legalize_triangles(t_can_v1_v2->v[0], t_can_v1_v2->edge[1], &stack_top);
        legalize_triangles(t_can_v2_v3->v[0], t_can_v2_v3->edge[1], &stack_top);
        legalize_triangles(t_can_v3_v1->v[0], t_can_v3_v1->edge[1], &stack_top);
    }
    else { // on the side
        Point *vi, *vj, *vk, *vl;
        Edge *eij, *ejk, *eki;
        Edge *eil, *elj, *eji;
        switch (dividing_point->position_to_triangle(triangle)) {
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
                PDASSERT(false);
                break;
        }
        //EXECUTION_REPORT(REPORT_ERROR, -1, dividing_point->position_to_edge(vi, vj) == 0, "point, which should be on the edge, is not on the edge");
        PDASSERT(dividing_point->position_to_edge(vi, vj) == 0);
        if (eij->twin_edge != NULL)
            PDASSERT(eij->twin_edge->triangle->is_leaf);
        ejk = eij->next_edge_in_triangle;
        eki = ejk->next_edge_in_triangle;
        if (eij->twin_edge != NULL) { 
            eji = eij->twin_edge;
            eij->twin_edge->triangle->is_leaf = false;
            eil = eji->next_edge_in_triangle;
            elj = eil->next_edge_in_triangle;
            vl = elj->head;
        }
        Edge *eir = allocate_edge(vi, dividing_point);
        Edge *erk = allocate_edge(dividing_point, vk);
        Edge *ekr = generate_twins_edge(erk);
        Edge *erj = allocate_edge(dividing_point, vj);
        Triangle* tirk = allocate_Triangle(eir, erk, eki);
        Triangle* tjkr = allocate_Triangle(ejk, ekr, erj);
        push(&stack_top, triangle);
        push(&stack_top, tirk);
        push(&stack_top, tjkr);
        legalize_triangles(dividing_point, ejk, &stack_top);
        legalize_triangles(dividing_point, eki, &stack_top); 
        if (eij->twin_edge != NULL) {
            Edge *eri = generate_twins_edge(eir);
            Edge *ejr = generate_twins_edge(erj);
            Edge *erl = allocate_edge(dividing_point, vl);
            Edge *elr = generate_twins_edge(erl);
            Triangle* tilr = allocate_Triangle(eil, elr, eri);
            Triangle* tjrl = allocate_Triangle(ejr, erl, elj);
            push(&stack_top, eij->twin_edge->triangle);
            push(&stack_top, tilr);
            push(&stack_top, tjrl);
            legalize_triangles(dividing_point, eil, &stack_top);
            legalize_triangles(dividing_point, elj, &stack_top);
        }
        else {
            //eir->twin_edge = NULL;
        }
    }

#ifdef DEBUG
    for (unsigned i = stack_base+1; i <= stack_top; i ++)
        if (triangle_stack[i]->is_leaf)
            PDASSERT(triangle_stack[i]->remained_points_head == -1 && triangle_stack[i]->remained_points_tail == -1);
#endif

    int list_head, list_tail;
    link_remained_list(stack_base, stack_top, &list_head, &list_tail);

    distribute_points_into_triangles(list_head, list_tail, stack_base, stack_top);

    //char filename[64];
    //snprintf(filename, 64, "log/single_step/step_%d.png", triangulate_count);
    //plot_current_step_into_file(filename);
    //printf("plot step %d\n", triangulate_count);
    //triangulate_count++;

    for (unsigned i = stack_base+1; i <= stack_top; i ++)
        if(triangle_stack[i]->is_leaf) {
            triangularization_process(triangle_stack[i], stack_top);
        }

}


static int compare_node_index(const void* a, const void* b)
{
    if(*(const int*)a < *(const int*)b) return -1;
    if(*(const int*)a > *(const int*)b) return 1;
    return 0;
}


void Delaunay_Voronoi::link_remained_list(unsigned base, unsigned top, int* head, int* tail)
{
    const unsigned max_leaf_triangles = 64;
    unsigned count;
    unsigned i;
    int head_tail[max_leaf_triangles*2]; // [head1, tail1, head2, tail2, ...]

    for (i = base+1, count = 0; i <= top; i ++)
        if (!triangle_stack[i]->is_leaf && triangle_stack[i]->remained_points_tail > -1) {
            head_tail[count * 2] = triangle_stack[i]->remained_points_head;
            head_tail[count * 2 + 1] = triangle_stack[i]->remained_points_tail;
            count++;
        }
    PDASSERT(count <= max_leaf_triangles);

#ifdef DEBUG
    bool* map = new bool[num_points]();
    for (i = base+1; i <= top; i ++)
        if (!triangle_stack[i]->is_leaf && triangle_stack[i]->remained_points_tail > -1) {
            for (int j = triangle_stack[i]->remained_points_head; j > -1; j = all_points[j].next) {
                PDASSERT(map[j] == false);
                map[j] = true;
            }
        }
    delete[] map;
#endif

    unsigned point_count = 0;
    if (count > 0) {
        merge_sort(head_tail, count, sizeof(int)*2, compare_node_index);
        for (unsigned i = 0; i < count; i++)
            for (int j = head_tail[i*2]; j > -1; j = all_points[j].next)
                point_count++;

        for (unsigned i = 0; i < count - 1; i++) {
            int cur_tail_id = head_tail[i*2+1];
            int nxt_head_id = head_tail[i*2+2];
            all_points[cur_tail_id].next = nxt_head_id;
            all_points[nxt_head_id].prev = cur_tail_id;
        }
        *head = head_tail[0];
        *tail = head_tail[count*2-1];
    } else {
        *head = -1;
        *tail = -1;
    }
#ifdef DEBUG
    unsigned point_count2 = 0;
    for (int i = *head; i > -1; i = all_points[i].next)
        point_count2++;
    PDASSERT(point_count == point_count2);
#endif
}


#define PDLN_INSERT_EXTRA_VPOINT (false)
#define PDLN_VPOINT_DENSITY  (20)
unsigned Delaunay_Voronoi::generate_initial_triangles(int num_points, double *x, double *y)
{
    double minX, maxX, minY, maxY;
    double dx, dy, deltaMax;
    unsigned stack_base = 0;
    unsigned stack_top  = 0;

    PDASSERT(x != NULL);
    PDASSERT(y != NULL);

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
    virtual_point[0] = new Point(minX-deltaMax*ratio, minY-deltaMax*ratio, -1);
    virtual_point[1] = new Point(minX-deltaMax*ratio, maxY+deltaMax*ratio, -1);
    virtual_point[2] = new Point(maxX+deltaMax*ratio, minY-deltaMax*ratio, -1);
    virtual_point[3] = new Point(maxX+deltaMax*ratio, maxY+deltaMax*ratio, -1);

    push(&stack_top, allocate_Triangle(allocate_edge(virtual_point[0], virtual_point[2]),
                                                  allocate_edge(virtual_point[2], virtual_point[1]),
                                                  allocate_edge(virtual_point[1], virtual_point[0])));
    push(&stack_top, allocate_Triangle(allocate_edge(virtual_point[3], virtual_point[1]),
                                                  allocate_edge(virtual_point[1], virtual_point[2]),
                                                  allocate_edge(virtual_point[2], virtual_point[3])));
    triangle_stack[stack_top-1]->edge[1]->twin_edge = triangle_stack[stack_top]->edge[1];
    triangle_stack[stack_top]->edge[1]->twin_edge = triangle_stack[stack_top-1]->edge[1];

    all_points = new Point[num_points];

    for (int i = 0; i < num_points; i ++) {
        all_points[i].x    = x[i];
        all_points[i].y    = y[i];
        all_points[i].id   = i;
        all_points[i].next = i + 1;
        all_points[i].prev = i - 1;
    }
    all_points[0].prev = all_points[num_points-1].next = -1;

    distribute_points_into_triangles(0, num_points-1, stack_base, stack_top);

    return stack_top;
}


void Delaunay_Voronoi::map_buffer_index_to_point_index()
{
    point_idx_to_buf_idx = new int[num_points]();
    for (int i = 0; i < num_points; i++) {
#ifdef DEBUG
        PDASSERT(point_idx_to_buf_idx[all_points[i].id] == 0);
#endif
        point_idx_to_buf_idx[all_points[i].id] = i;
    }

#ifdef DEBUG
    for (int i = 0; i < num_points; i++) {
        PDASSERT(x_store[i] == all_points[point_idx_to_buf_idx[i]].x);
        PDASSERT(y_store[i] == all_points[point_idx_to_buf_idx[i]].y);
    }
#endif
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
                                   double min_lat, double max_lat, int virtual_polar_local_index)
    : num_points(num_points)
    , tolerance(FLOAT_ERROR)
    , x_ref(x_origin)
    , y_ref(y_origin)
{
    timeval start, end;

    PDASSERT((x_ref && y_ref) || (!x_ref && !x_ref));
#ifdef DEBUG
    PDASSERT(have_redundent_points(x_values, y_values, num_points) == false);
    x_store = x_values;
    y_store = y_values;

    for(int i = 0; i < num_points; i++)
        PDASSERT(x_values[i] >= min_lon && x_values[i] <= max_lon && y_values[i] >= min_lat && y_values[i] <= max_lat);
#endif
    int triangles_count_estimate = 2*(num_points+4);
    triangle_pool.reserve(triangles_count_estimate);
    edge_pool.reserve(triangles_count_estimate*3);
    result_leaf_triangles.reserve(triangles_count_estimate);
    result_triangles_pack.reserve(triangles_count_estimate);

    stack_size = triangles_count_estimate * 2;
    triangle_stack = new Triangle*[stack_size];

    gettimeofday(&start, NULL);

    this->is_global_grid = is_global_grid;

    if(global_idx != NULL) {
        global_index = new int[num_points];
        memcpy(global_index, global_idx, num_points*sizeof(int));
    }

    this->vpolar_local_index = virtual_polar_local_index;
    unsigned stack_top = generate_initial_triangles(num_points, x_values, y_values);

    //save_original_points_into_file();

    for(unsigned i = 1; i <= stack_top; i++)
        if(triangle_stack[i]->is_leaf)
            triangularization_process(triangle_stack[i], stack_top);

    map_buffer_index_to_point_index();
    clear_triangle_containing_virtual_point();

    update_virtual_polar_info();

    gettimeofday(&end, NULL);
}


Delaunay_Voronoi::~Delaunay_Voronoi()
{
    delete[] all_points;
    delete[] global_index;
    for (int i = 0; i < 4; i++)
        delete virtual_point[i];
    delete[] point_idx_to_buf_idx;
    delete[] triangle_stack;
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

        for(unsigned j = 0; j < 3; j++) {
            relegalize_triangles(result_leaf_triangles[i]->v[j], result_leaf_triangles[i]->edge[(j+1)%3]);
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


void Delaunay_Voronoi::get_triangles_intersecting_with_segment(Point head, Point tail, Triangle_pack *output_triangles, int *num_triangles, int buf_len, double threshold)
{
    int current = 0;

    std::vector<Triangle*> ts = find_triangles_intersecting_with_segment(head, tail, threshold);

    for(unsigned i = 0; i < ts.size(); i++)
        output_triangles[current++] = Triangle_pack(Point(ts[i]->v[0]->x, ts[i]->v[0]->y, global_index[ts[i]->v[0]->id]),  // TODO: use global_index as id directly
                                                         Point(ts[i]->v[1]->x, ts[i]->v[1]->y, global_index[ts[i]->v[1]->id]), 
                                                         Point(ts[i]->v[2]->x, ts[i]->v[2]->y, global_index[ts[i]->v[2]->id]));

    PDASSERT(current < buf_len);
    *num_triangles = current;
}


static inline unsigned hash_triangle_by_id(Triangle_pack triangle)
{
    return triangle.v[0].id ^ (triangle.v[1].id << 10) ^ (triangle.v[2].id << 20);
}


inline double ppoint_distence_to_line(double px, double py, double x1, double y1, double x2, double y2)
{
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

inline bool point_distence_in_threshold(Point v, Point p1, Point p2, double threshold) {
    return ppoint_distence_to_line(v.x, v.y, p1.x, p1.y, p2.x, p2.y) <= threshold;
}

inline bool all_distence_in_threshold(Triangle_pack* triangle, Point p1, Point p2, double threshold) {
    return point_distence_in_threshold(triangle->v[0], p1, p2, threshold) &&
           point_distence_in_threshold(triangle->v[1], p1, p2, threshold) &&
           point_distence_in_threshold(triangle->v[2], p1, p2, threshold);
}

inline bool vertexs_not_on_same_side(Triangle_pack* triangle, Point p1, Point p2) {
    return triangle->v[0].position_to_edge(&p1, &p2) * triangle->v[1].position_to_edge(&p1, &p2) <= 0 ||
           triangle->v[1].position_to_edge(&p1, &p2) * triangle->v[2].position_to_edge(&p1, &p2) <= 0;
}

inline bool is_segment_in_triangle(Triangle_pack* triangle, Point p1, Point p2) {
    return p1.position_to_triangle(triangle) >= 0 && p2.position_to_triangle(triangle) >= 0;
}

inline bool is_segment_intersected_with_edge(Point p_e1, Point p_e2, Point p1, Point p2) {
    return (p1.position_to_edge(&p_e1, &p_e2) * p2.position_to_edge(&p_e1, &p_e2) < 0) &&
           (p1.position_to_edge(&p_e1, &p_e2)!=0 || p2.position_to_edge(&p_e1, &p_e2)!=0 ) &&
           (p_e1.position_to_edge(&p1, &p2) * p_e2.position_to_edge(&p1, &p2) <= 0);
}

inline bool is_segment_intersected_with_one_of_edges(Triangle_pack* triangle, Point p1, Point p2) {
    return is_segment_intersected_with_edge(triangle->v[0], triangle->v[1], p1, p2) ||
           is_segment_intersected_with_edge(triangle->v[1], triangle->v[2], p1, p2) ||
           is_segment_intersected_with_edge(triangle->v[2], triangle->v[0], p1, p2);
}

inline bool is_triangle_intersecting_with_segment(Triangle_pack* triangle, Point p1, Point p2, double threshold) {
    return (threshold == 0 || all_distence_in_threshold(triangle, p1, p2, threshold)) &&
           vertexs_not_on_same_side(triangle, p1, p2) &&
           (is_segment_intersected_with_one_of_edges(triangle, p1, p2) || is_segment_in_triangle(triangle, p1, p2));
}

unsigned Delaunay_Voronoi::calculate_checksum(Point head, Point tail, double threshold)
{
    PDASSERT(result_triangles_pack.size() > 0);
    Point uncyclic_head = head, uncyclic_tail = tail;

    if(uncyclic_head.x == uncyclic_tail.x) {
        while(uncyclic_head.x >= 360) uncyclic_head.x -= 360;
        while(uncyclic_head.x < 0) uncyclic_head.x += 360;
        while(uncyclic_tail.x >= 360) uncyclic_tail.x -= 360;
        while(uncyclic_tail.x < 0) uncyclic_tail.x += 360;
    }

    unsigned checksum = 0;
    for(unsigned i = 0; i < result_triangles_pack.size();) {
        if (result_triangles_pack[i].is_cyclic) {
            if(is_triangle_intersecting_with_segment(&result_triangles_pack[i+1], uncyclic_head, uncyclic_tail, threshold) ||
               is_triangle_intersecting_with_segment(&result_triangles_pack[i+2], uncyclic_head, uncyclic_tail, threshold) )
                checksum ^= hash_triangle_by_id(result_triangles_pack[i]);
            i += 3;
        } else {
            if (is_triangle_intersecting_with_segment(&result_triangles_pack[i], head, tail, threshold))
                checksum ^= hash_triangle_by_id(result_triangles_pack[i]);
            i++;
        }
    }

    return checksum;
}


/* including triangles intersecting with boundary */
void Delaunay_Voronoi::get_triangles_in_region(double min_x, double max_x, double min_y, double max_y, 
                                               Triangle_pack *output_triangles, int *output_num_triangles, int buf_len)
{
    int current = 0;
    int num_triangles = 0;

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
            output_triangles[current++] = Triangle_pack(Point(result_leaf_triangles[i]->v[0]->x, result_leaf_triangles[i]->v[0]->y, global_index[result_leaf_triangles[i]->v[0]->id]),
                                                             Point(result_leaf_triangles[i]->v[1]->x, result_leaf_triangles[i]->v[1]->y, global_index[result_leaf_triangles[i]->v[1]->id]),
                                                             Point(result_leaf_triangles[i]->v[2]->x, result_leaf_triangles[i]->v[2]->y, global_index[result_leaf_triangles[i]->v[2]->id]));
        }
    }
    PDASSERT(current < buf_len);

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
        if(result_leaf_triangles[i]->contain_vertex(&all_points[point_idx_to_buf_idx[vpolar_local_index]]))
            triangles_containing_vpolar.push_back(result_leaf_triangles[i]);
    }
}


void Delaunay_Voronoi::make_final_triangle_pack()
{
    for(unsigned i = 0; i < result_leaf_triangles.size(); i++) {
        if(!result_leaf_triangles[i]->is_leaf)
            continue;

        Point**   v = result_leaf_triangles[i]->v;
        bool cyclic = result_leaf_triangles[i]->is_cyclic;
        result_triangles_pack.push_back(Triangle_pack(Point(v[0]->x, v[0]->y, global_index[v[0]->id]),
                                                      Point(v[1]->x, v[1]->y, global_index[v[1]->id]),
                                                      Point(v[2]->x, v[2]->y, global_index[v[2]->id]),
                                                      cyclic));
        sort_points_in_triangle(result_triangles_pack.back());
        if (cyclic) {
            double tmp_x[3];
            Triangle_pack t_cur = result_triangles_pack.back();
            for(int j = 0; j < 3; j++) {
                if(t_cur.v[j].x >= 180)
                    tmp_x[j] = t_cur.v[j].x - 360;
                else
                    tmp_x[j] = t_cur.v[j].x;
            }
            result_triangles_pack.push_back(Triangle_pack(Point(tmp_x[0], t_cur.v[0].y, t_cur.v[0].id),
                                                          Point(tmp_x[1], t_cur.v[1].y, t_cur.v[1].id),
                                                          Point(tmp_x[2], t_cur.v[2].y, t_cur.v[2].id),
                                                          cyclic));
            for(int j = 0; j < 3; j++) {
                if(t_cur.v[j].x < 180)
                    tmp_x[j] = t_cur.v[j].x + 360;
                else
                    tmp_x[j] = t_cur.v[j].x;
            }
            result_triangles_pack.push_back(Triangle_pack(Point(tmp_x[0], t_cur.v[0].y, t_cur.v[0].id),
                                                          Point(tmp_x[1], t_cur.v[1].y, t_cur.v[1].id),
                                                          Point(tmp_x[2], t_cur.v[2].y, t_cur.v[2].id),
                                                          cyclic));
        }
    }
}


void Delaunay_Voronoi::remove_triangles_only_containing_virtual_polar()
{
    double common_lat = all_points[point_idx_to_buf_idx[vpolar_local_index]].y;

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

    PDASSERT(num_edges%3 == 0);
    PDASSERT(num_edges <= 3 * result_leaf_triangles.size());
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

    PDASSERT(num_edges%3 == 0);
    PDASSERT(num_edges <= 3 * triangle_pool.size());
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

    PDASSERT(num_edges <= 3 * result_leaf_triangles.size());
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

    coord[0] = new double[num_points];
    coord[1] = new double[num_points];

    int num = 0;
    for(int i = 0; i < num_points; i ++) {
        coord[0][num] = all_points[i].x;
        coord[1][num++] = all_points[i].y;
    }

    PDASSERT(num == num_points);
    plot_points_into_file(filename, coord[0], coord[1], num_points, min_x, max_x, min_y, max_y);

    delete coord[0];
    delete coord[1];
}


void plot_triangles_into_file(const char *prefix, Triangle_pack *t, int num, bool plot_cyclic_triangles)
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

    PDASSERT(num_edges%3 == 0);
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

    PDASSERT(num_edges%3 == 0);
    snprintf(filename, 128, "%s.png", prefix);
    //plot_edge_into_file(filename, head_coord, tail_coord, num_edges);
    plot_projected_edge_into_file(filename, head_coord, tail_coord, num_edges, PDLN_PLOT_COLOR_WHITE, PDLN_PLOT_FILEMODE_NEW);

    delete head_coord[0];
    delete head_coord[1];
    delete tail_coord[0];
    delete tail_coord[1];

}
#endif


void Delaunay_Voronoi::save_original_points_into_file()
{
    FILE *fp;
    fp = fopen("log/original_points.txt", "w");
    for(int i = 0; i < num_points; i++)
        if(all_points[i].x < 5 && all_points[i].x > -5 && all_points[i].y < 5 && all_points[i].y > -5)
            fprintf(fp, "%.20lf, %.20lf\n", all_points[i].x, all_points[i].y);
    fclose(fp);
}


void Delaunay_Voronoi::update_all_points_coord(double *x_values, double *y_values, int num)
{
    PDASSERT(num == num_points);
    for(int i = 0; i < num; i++) {
        all_points[point_idx_to_buf_idx[i]].x = x_values[i];
        all_points[point_idx_to_buf_idx[i]].y = y_values[i];
    }

    for(unsigned i = 0; i < result_leaf_triangles.size(); i++)
        if(result_leaf_triangles[i]->is_leaf)
            result_leaf_triangles[i]->calulate_circum_circle();
}


/* NOTE: circum_circle will not been recalculated */
void Delaunay_Voronoi::update_points_coord_y(double reset_lat_value, vector<int> *polars_local_index)
{
    for(unsigned i = 0; i < polars_local_index->size(); i++)
        all_points[point_idx_to_buf_idx[(*polars_local_index)[i]]].y = reset_lat_value;
}


void Delaunay_Voronoi::uncyclic_all_points()
{
    for(int i = 0; i < num_points; i++) {
        while(all_points[i].x >= 360)
            all_points[i].x -= 360;
        while(all_points[i].x < 0)
            all_points[i].x += 360;
    }
}

Triangle_pack::Triangle_pack(Point p0, Point p1, Point p2, bool cyclic)
{
    v[0] = p0;
    v[1] = p1;
    v[2] = p2;
    is_cyclic = cyclic;
}


bool operator == (Triangle_pack t1, Triangle_pack t2)
{
#ifdef DEBUG
    PDASSERT(t1.v[0] != t1.v[1] && t1.v[1] != t1.v[2] && t1.v[2] != t1.v[0]);
    PDASSERT(t2.v[0] != t2.v[1] && t2.v[1] != t2.v[2] && t2.v[2] != t2.v[0]);
#endif
    if(t2.v[0] != t1.v[0] && t2.v[0] != t1.v[1] && t2.v[0] != t1.v[2])
        return false;
    if(t2.v[1] != t1.v[0] && t2.v[1] != t1.v[1] && t2.v[1] != t1.v[2])
        return false;
    if(t2.v[2] != t1.v[0] && t2.v[2] != t1.v[1] && t2.v[2] != t1.v[2])
        return false;
    return true;
}


void save_triangles_info_file(const char *prefix, Triangle_pack *t, int num)
{
    char filename[128];
    FILE* fp;

    snprintf(filename, 128, "%s.txt", prefix);
    fp = fopen(filename, "w");
    for(int i = 0; i < num; i++)
        fprintf(fp, "[%d] (%.15lf, %.15lf), (%.15lf, %.15lf), (%.15lf, %.15lf)\n", i, t[i].v[0].x, t[i].v[0].y, t[i].v[1].x, t[i].v[1].y, t[i].v[2].x, t[i].v[2].y);
    fclose(fp);
}
