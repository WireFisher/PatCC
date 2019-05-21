#include "mpi.h"
#include "delaunay_voronoi_2D.h"
#include "opencv_utils.h"
#include "common_utils.h"
#include <cmath>
#include <cstdio>
#include <cstring>
#include <sys/time.h>
#include <tr1/unordered_map>
#include <list>
#include "merge_sort.h"
#include <utility>
#include "coordinate_hash.h"

#define PAT_NUM_LOCAL_VPOINTS (4)
#define PAT_CYCLIC_EDGE_THRESHOLD (180)

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
    Triangle_inline t1 = *(const Triangle_inline*)a;
    Triangle_inline t2 = *(const Triangle_inline*)b;

    if(t1.v[2].id < t2.v[2].id) return -1;
    if(t1.v[2].id > t2.v[2].id) return  1;
    return 0;
}


static int compare_v1(const void* a, const void* b)
{
    Triangle_inline t1 = *(const Triangle_inline*)a;
    Triangle_inline t2 = *(const Triangle_inline*)b;

    if(t1.v[1].id < t2.v[1].id) return -1;
    if(t1.v[1].id > t2.v[1].id) return  1;
    return 0;
}


static int compare_v0(const void* a, const void* b)
{
    Triangle_inline t1 = *(const Triangle_inline*)a;
    Triangle_inline t2 = *(const Triangle_inline*)b;

    if(t1.v[0].id < t2.v[0].id) return -1;
    if(t1.v[0].id > t2.v[0].id) return  1;
    return 0;
}


static inline void radix_sort(Triangle_inline *triangles, int num_triangles)
{
    PDASSERT(sizeof(Triangle_inline) > sizeof(void *)/2);
    merge_sort(triangles, num_triangles, sizeof(Triangle_inline), compare_v2);
    merge_sort(triangles, num_triangles, sizeof(Triangle_inline), compare_v1);
    merge_sort(triangles, num_triangles, sizeof(Triangle_inline), compare_v0);
}


void sort_points_in_triangle(Triangle_inline *triangles, int num_triangles)
{
    for(int i = 0; i < num_triangles; i++) {
        if(triangles[i].v[0].id > triangles[i].v[1].id) swap(&triangles[i].v[0], &triangles[i].v[1]);
        if(triangles[i].v[1].id > triangles[i].v[2].id) swap(&triangles[i].v[1], &triangles[i].v[2]);
        if(triangles[i].v[0].id > triangles[i].v[1].id) swap(&triangles[i].v[0], &triangles[i].v[1]);
    }
}


inline void sort_points_in_triangle(Triangle_inline& triangle)
{
    if(triangle.v[0].id > triangle.v[1].id) std::swap(triangle.v[0], triangle.v[1]);
    if(triangle.v[1].id > triangle.v[2].id) std::swap(triangle.v[1], triangle.v[2]);
    if(triangle.v[0].id > triangle.v[1].id) std::swap(triangle.v[0], triangle.v[1]);
}


void sort_triangles(Triangle_inline *triangles, int num_triangles)
{
    radix_sort(triangles, num_triangles);
}


Point::Point()
{
}


Point::Point(double x, double y)
    : x(x)
    , y(y)
    , mask(true)
{
}


Point::Point(double x, double y, int id, int fd, int bk)
    : x(x)
    , y(y)
    , id(id)
    , mask(true)
    , next(fd)
    , prev(bk)
{
}


Point::Point(double x, double y, int id, bool msk, int fd, int bk)
    : x(x)
    , y(y)
    , id(id)
    , mask(msk)
    , next(fd)
    , prev(bk)
{
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


static inline double calculate_distance(double x0, double y0, double x1, double y1)
{
    double dx = x0 - x1;
    double dy = y0 - y1;
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

    if(float_eq(p1.x, p2.x) && float_eq(p1.y, p2.y))
        return true;
    else if(float_eq(fabs(p1.x - p2.x), 360.0) && float_eq(p1.y, p2.y))
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
    //double delta1_x = pt2->x - pt1->x;
    //double delta1_y = pt2->y - pt1->y;
    //double delta2_x = x - pt1->x;
    //double delta2_y = y - pt1->y;

    //double left  = delta1_x*delta2_y;
    //double right = delta2_x*delta1_y;

    //if (relative_eq_hi(left, right))
    double res1 = det(pt1, pt2, this);
    if (float_eq_hi(res1, 0))
        return 0;
    //else if (left > right)
    else if (res1 > 0)
        return 1;
    else
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
int Point::position_to_triangle(const Point* v0, const Point* v1, const Point* v2) const
{
#ifdef DDEBUG
    bool on1 = position_to_edge(v0, v1) == 0;
    bool on2 = position_to_edge(v1, v2) == 0;
    bool on3 = position_to_edge(v2, v0) == 0;
    PDASSERT(!(on1 && on2));
    PDASSERT(!(on2 && on3));
    PDASSERT(!(on3 && on1));
#endif
    int ret = 0;
    int pos = position_to_edge(v0, v1);
    if (pos == -1)
        return -1;
    else if (pos == 0)
        ret = 1;
    pos = position_to_edge(v1, v2);
    if (pos == -1)
        return -1;
    else if (pos == 0)
        ret = 2;
    pos = position_to_edge(v2, v0);
    if (pos == -1)
        return -1;
    else if (pos == 0)
        ret = 3;
    return ret;
}


int Point::position_to_triangle(const Triangle_inline *triangle) const
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

Edge::Edge()
{
}


Edge::~Edge()
{
}


Edge* Delaunay_Voronoi::make_twins_edge(Edge *e)
{
    Edge *twins_edge = allocate_edge(e->tail, e->head);
    twins_edge->twin_edge = e;
    e->twin_edge = twins_edge;

    return twins_edge;
}


int Delaunay_Voronoi::get_lowest_point_of_four(int idx1, int idx2, int idx3, int idx4)
{
    int idx = idx1;

    Point* p  = &all_points[idx1];
    Point* p2 = &all_points[idx2];
    Point* p3 = &all_points[idx3];
    Point* p4 = &all_points[idx4];

    if(!x_ref) {
        if(p2->x < p->x || (p2->x == p->x && p2->y < p->y)) { p = p2; idx = idx2;}
        if(p3->x < p->x || (p3->x == p->x && p3->y < p->y)) { p = p3; idx = idx3;}
        if(p4->x < p->x || (p4->x == p->x && p4->y < p->y)) { p = p4; idx = idx4;}
    } else {
        if(x_ref[p2->id] < x_ref[p->id] || (x_ref[p2->id] == x_ref[p->id] && y_ref[p2->id] < y_ref[p->id])) { p = p2; idx = idx2;}
        if(x_ref[p3->id] < x_ref[p->id] || (x_ref[p3->id] == x_ref[p->id] && y_ref[p3->id] < y_ref[p->id])) { p = p3; idx = idx3;}
        if(x_ref[p4->id] < x_ref[p->id] || (x_ref[p4->id] == x_ref[p->id] && y_ref[p4->id] < y_ref[p->id])) { p = p4; idx = idx4;}
        //if ((x_ref[p1->id] > 345 || x_ref[p1->id] < 0.5) && y_ref[p1->id] > 77.6 && y_ref[p1->id] < 78.2 &&
        //    (x_ref[p2->id] > 345 || x_ref[p2->id] < 0.5) && y_ref[p2->id] > 77.6 && y_ref[p2->id] < 78.2 &&
        //    (x_ref[p3->id] > 345 || x_ref[p3->id] < 0.5) && y_ref[p3->id] > 77.6 && y_ref[p3->id] < 78.2 &&
        //    (x_ref[p4->id] > 345 || x_ref[p4->id] < 0.5) && y_ref[p4->id] > 77.6 && y_ref[p4->id] < 78.2)
        //    printf("lowest: %p, [%lf, %lf] in [%.15lf, %.15lf] [%.15lf, %.15lf] [%.15lf, %.15lf] [%.15lf, %.15lf]\n", p, x_ref[p->id], y_ref[p->id],
        //                                                                                         x_ref[p1->id], y_ref[p1->id],
        //                                                                                         x_ref[p2->id], y_ref[p2->id],
        //                                                                                         x_ref[p3->id], y_ref[p3->id],
        //                                                                                         x_ref[p4->id], y_ref[p4->id]);
    }
    return idx;
}




/* Debug staff */
int triangulate_count = 0;
bool Print_Error_info=false;
int on_circle_count = 0;
int not_on_circle_count = 0;
extern double global_p_lon[4];
extern double global_p_lat[4];
/*
#define point_is(pt, a, b) (fabs(pt->x - a) < 1e-5 && fabs(pt->y - b) < 1e-5)
#define edge_is(e, a, b, c, d) ((point_is(e->head, a, b) && point_is(e->tail, c, d)) || (point_is(e->tail, a, b) && point_is(e->head, c, d)))
#define triangle_is(t, a, b, c, d ,e, f) ((point_is(t->v[0], a, b) || point_is(t->v[0], c, d) || point_is(t->v[0], e, f)) && \
                                          (point_is(t->v[1], a, b) || point_is(t->v[1], c, d) || point_is(t->v[1], e, f)) && \
                                          (point_is(t->v[2], a, b) || point_is(t->v[2], c, d) || point_is(t->v[2], e, f)) )
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


bool Delaunay_Voronoi::is_edge_legal(int p_idx, const Edge *edge)
{
#ifdef DEBUG
    reason = SUCCESS;
#endif
    PDASSERT(edge->triangle->is_leaf);
    if (!edge->twin_edge) {
        return true;
    }

    if(!edge->twin_edge->triangle) {
        return true;
    }

    if(!edge->twin_edge->triangle->is_leaf) {
        return true;
    }

    int ret = circum_circle_contains_reliably(edge, head(edge->twin_edge->prev_edge_in_triangle), tolerance);

    if (ret == -1) {
        return true;
    }

    if (ret == 0) {
        return check_uniqueness(p_idx, edge);
    }

#ifdef DEBUG
    reason = IN_CIRCUMCIRCLE;
#endif
    return false;
}

bool Delaunay_Voronoi::check_uniqueness(int p_idx, const Edge *edge)
{
    int lowest = get_lowest_point_of_four(p_idx, edge->head, edge->tail, edge->twin_edge->prev_edge_in_triangle->head);
    bool is_lowest = (p_idx == lowest || edge->twin_edge->prev_edge_in_triangle->head == lowest);

    //if(edge->triangle->is_cyclic && edge->twin_edge->triangle->is_cyclic)
    //    is_lowest = !is_lowest;

#ifdef DEBUG
    if(!is_lowest)
        reason = AMBIGUOUS;
#endif

    return is_lowest;
}


bool Delaunay_Voronoi::is_angle_ambiguous(int p_idx, const Edge *edge)
{
    int lowest = get_lowest_point_of_four(p_idx, edge->head, edge->tail, edge->twin_edge->prev_edge_in_triangle->head);
    if(p_idx == lowest || edge->twin_edge->prev_edge_in_triangle->head == lowest) {
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

bool Delaunay_Voronoi::is_triangle_ambiguous(int p_idx, Edge *edge)
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

    if(edge->twin_edge->triangle->is_virtual) {
        return false;
    }

    Point* p[4];
    bool modified[4] = {0};
    p[0] = &all_points[p_idx];
    p[1] = head(edge);
    p[2] = tail(edge);
    p[3] = head(edge->twin_edge->prev_edge_in_triangle);

    PDASSERT(!modified[0] && !modified[1] && !modified[2] && !modified[3]);

    int ret;
    if (edge->triangle->is_cyclic && edge->twin_edge->triangle->is_cyclic) {
        for (int i = 0; i < 4; i++)
            if (p[i]->x > 180) {
                p[i]->x -= 360;
                modified[i] = true;
            }
        ret = circum_circle_contains_reliably(edge, head(edge->twin_edge->prev_edge_in_triangle), tolerance);
        for (int i = 0; i < 4; i++)
            if (modified[i])
                p[i]->x += 360;
    } else {
        ret = circum_circle_contains_reliably(edge, head(edge->twin_edge->prev_edge_in_triangle), tolerance);
    }

    if (ret == 0) {
        return is_angle_ambiguous(p_idx, edge);
    }

    return false;
}


bool Delaunay_Voronoi::is_triangle_legal(const Triangle *t)
{
    for(int i = 0; i < 3; i++)
        if(!is_edge_legal(t->edge[i]->prev_edge_in_triangle->head, t->edge[i])) {
            //printf("[%d] +illegal triangle: (%lf, %lf), (%lf, %lf), (%lf, %lf)\n", 1, vertex(t, 0)->x, vertex(t, 0)->y, vertex(t, 1)->x, vertex(t, 1)->y, vertex(t, 2)->x, vertex(t, 2)->y);
            Triangle *tt = t->edge[i]->twin_edge->triangle;
            //printf("[%d] -illegal triangle: (%lf, %lf), (%lf, %lf), (%lf, %lf)\n", 1, vertex(tt, 0)->x, vertex(tt, 0)->y, vertex(tt, 1)->x, vertex(tt, 1)->y, vertex(tt, 2)->x, vertex(tt, 2)->y);
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
    for(unsigned i = 0; i < all_leaf_triangles.size(); i++)
        if(all_leaf_triangles[i]->is_leaf && !all_leaf_triangles[i]->is_virtual)
            if(!is_triangle_legal(all_leaf_triangles[i])) {
                is_legal = false;
                //printf("illegal reason: %d\n", reason);
                //fprintf(stderr, "[%d] illegal: (%lf, %lf), (%lf, %lf), (%lf, %lf)\n", rank, vertex(all_leaf_triangles[i], 0)->x, vertex(all_leaf_triangles[i], 0)->y, vertex(all_leaf_triangles[i], 1)->x, vertex(all_leaf_triangles[i], 1)->y, vertex(all_leaf_triangles[i], 2)->x, vertex(all_leaf_triangles[i], 2)->y);
            }
    if(is_legal)
        return true;
    else
        return false;
}


bool Delaunay_Voronoi::is_delaunay_legal(const Triangle *t)
{
    for(int i = 0; i < 3; i++)
        if(!is_delaunay_legal(head(t->edge[i]->prev_edge_in_triangle), t->edge[i]))
            return false;
    return true;
}


bool Delaunay_Voronoi::is_delaunay_legal(const Point *pt, const Edge *edge)
{
    PDASSERT(edge->triangle->is_leaf);
    if (!edge->twin_edge) {
        return true;
    }

    if(!edge->twin_edge->triangle) {
        return true;
    }

    if(!edge->twin_edge->triangle->is_leaf) {
        return true;
    }

    int ret = circum_circle_contains_reliably(edge, head(edge->twin_edge->prev_edge_in_triangle), tolerance);
    if (ret == 1)
        return false;
    
    return true;
}

/*
 *
 *    Vi ________ Vk
 *      /\      /
 *     /  \ 2  /
 *    / 1  \  /
 *   /______\/
 *  Vr      Vj
 *
 *  1: triangle
 *  2: twin triangle
 */
void Delaunay_Voronoi::legalize_triangles(int vr_idx, Edge *edge, unsigned stack_base, unsigned *stack_top)
{
    if (is_edge_legal(vr_idx, edge))
        return;

    PDASSERT(edge->triangle->is_leaf);
    PDASSERT(edge->twin_edge->triangle->is_leaf);
    push(stack_top, edge->twin_edge->triangle);

    edge->triangle->is_leaf = false;
    edge->twin_edge->triangle->is_leaf = false;

    int vk_idx = edge->twin_edge->prev_edge_in_triangle->head;
    Edge *eij = edge;
    Edge *ejr = eij->next_edge_in_triangle;
    Edge *eri = ejr->next_edge_in_triangle;
    Edge *eji = eij->twin_edge;
    Edge *eik = eji->next_edge_in_triangle;
    Edge *ekj = eik->next_edge_in_triangle;
    Edge *erk = allocate_edge(vr_idx, vk_idx);
    Edge *ekr = make_twins_edge(erk);
    Triangle* tikr = allocate_triangle(eik,ekr,eri);
    Triangle* tjrk = allocate_triangle(ejr,erk,ekj);
    push(stack_top, tikr);
    push(stack_top, tjrk);

    legalize_triangles(vr_idx, eik, stack_base, stack_top);
    legalize_triangles(vr_idx, ekj, stack_base, stack_top);
}


void Delaunay_Voronoi::relegalize_triangles(int vr_idx, Edge *edge)
{
    if (!is_triangle_ambiguous(vr_idx, edge))
        return;

    if (edge->triangle->is_cyclic != edge->twin_edge->triangle->is_cyclic)
        return;
    PDASSERT(edge->triangle->is_leaf);
    PDASSERT(edge->twin_edge->triangle->is_leaf);

    Triangle *triangle_origin, *triangle_twin;
    triangle_origin = edge->triangle;
    triangle_twin = edge->twin_edge->triangle;

    int vk_idx = edge->twin_edge->prev_edge_in_triangle->head;
    Edge *eij = edge;
    Edge *ejr = eij->next_edge_in_triangle;
    Edge *eri = ejr->next_edge_in_triangle;
    Edge *eji = eij->twin_edge;
    Edge *eik = eji->next_edge_in_triangle;
    Edge *ekj = eik->next_edge_in_triangle;
    Edge *erk = allocate_edge(vr_idx, vk_idx);
    Edge *ekr = make_twins_edge(erk);

    bool force = edge->triangle->is_cyclic ? true : false;
    initialize_triangle_with_edges(triangle_origin, eik, ekr, eri, force);
    initialize_triangle_with_edges(triangle_twin, ejr, erk, ekj, force);
    if(force)
        triangle_origin->is_cyclic = triangle_twin->is_cyclic = true;

}


Triangle::Triangle()
    : stack_ref_count(0)
{
    edge[0] = NULL;
    edge[1] = NULL;
    edge[2] = NULL;
}


Triangle::~Triangle()
{
}


void Triangle::calulate_circum_circle(const Point* v0, const Point* v1, const Point* v2)
{
    PAT_REAL ab, cd, ef;

    ab = ((PAT_REAL)v0->x * (PAT_REAL)v0->x) + ((PAT_REAL)v0->y * (PAT_REAL)v0->y);
    cd = ((PAT_REAL)v1->x * (PAT_REAL)v1->x) + ((PAT_REAL)v1->y * (PAT_REAL)v1->y);
    ef = ((PAT_REAL)v2->x * (PAT_REAL)v2->x) + ((PAT_REAL)v2->y * (PAT_REAL)v2->y);

    circum_center[0] = (ab * ((PAT_REAL)v2->y - (PAT_REAL)v1->y) + cd * ((PAT_REAL)v0->y - (PAT_REAL)v2->y) + ef * ((PAT_REAL)v1->y - (PAT_REAL)v0->y)) /
                       ((PAT_REAL)v0->x * ((PAT_REAL)v2->y - (PAT_REAL)v1->y) + (PAT_REAL)v1->x * ((PAT_REAL)v0->y - (PAT_REAL)v2->y) + (PAT_REAL)v2->x * ((PAT_REAL)v1->y - (PAT_REAL)v0->y)) * 0.5;
    circum_center[1] = (ab * ((PAT_REAL)v2->x - (PAT_REAL)v1->x) + cd * ((PAT_REAL)v0->x - (PAT_REAL)v2->x) + ef * ((PAT_REAL)v1->x - (PAT_REAL)v0->x)) /
                       ((PAT_REAL)v0->y * ((PAT_REAL)v2->x - (PAT_REAL)v1->x) + (PAT_REAL)v1->y * ((PAT_REAL)v0->x - (PAT_REAL)v2->x) + (PAT_REAL)v2->y * ((PAT_REAL)v1->x - (PAT_REAL)v0->x)) * 0.5;
    circum_radius2 = (((PAT_REAL)v0->x - circum_center[0]) * ((PAT_REAL)v0->x - circum_center[0])) + (((PAT_REAL)v0->y - circum_center[1]) * ((PAT_REAL)v0->y - circum_center[1]));
}


void Delaunay_Voronoi::initialize_triangle_with_edges(Triangle* t, Edge *edge1, Edge *edge2, Edge *edge3, bool force)
{
    Point *pt1, *pt2, *pt3;

    t->is_leaf = true;
    t->is_cyclic = false;
    t->is_virtual = false;
    if(force) {
        t->v[0] = edge1->head;
        t->v[1] = edge2->head;
        t->v[2] = edge3->head;
        t->edge[0] = edge1;
        t->edge[1] = edge2;
        t->edge[2] = edge3;
        pt1 = vertex(t, 0);
        pt2 = vertex(t, 1);
        pt3 = vertex(t, 2);
    } else {
        pt1 = head(edge1);
        pt2 = head(edge2);
        pt3 = head(edge3);

#ifdef DEBUG
        if(float_eq_hi(det(pt1, pt2, pt3), 0) || float_eq_hi(det(pt2, pt3, pt1), 0) || float_eq_hi(det(pt3, pt1, pt2), 0)) {
            fprintf(stderr, "(%.20lf, %.20lf), (%.20lf, %.20lf), (%.20lf, %.20lf)\n", pt1->x, pt1->y, pt2->x, pt2->y, pt3->x, pt3->y);
            //fprintf(stderr, "(%.20lf, %.20lf, %d), (%.20lf, %.20lf, %d), (%.20lf, %.20lf, %d)\n", x_ref[pt1->id], y_ref[pt1->id], global_index[pt1->id], x_ref[pt2->id], y_ref[pt2->id], global_index[pt2->id], x_ref[pt3->id], y_ref[pt3->id], global_index[pt3->id]);
            fprintf(stderr, "std::fabs(det(pt1, pt2, pt3)): %.20lf\nstd::fabs(det(pt2, pt3, pt1)): %.20lf\nstd::fabs(det(pt3, pt1, pt2)): %.20lf\n", std::fabs(det(pt1, pt2, pt3)),
                                                                                                                               std::fabs(det(pt2, pt3, pt1)),
                                                                                                                               std::fabs(det(pt3, pt1, pt2)));
            PDASSERT(false);
        }
#endif
        /* if there are unmarked redundant points, the PDASSERTion may fail */
        PDASSERT(!float_eq_hi(det(pt1, pt2, pt3), 0) && !float_eq_hi(det(pt2, pt3, pt1), 0) && !float_eq_hi(det(pt3, pt1, pt2), 0));
#ifdef DEBUG
        if(!(edge1->tail==edge2->head && edge2->tail==edge3->head && edge3->tail==edge1->head)) {
            printf("1: (%lf, %lf)->(%lf, %lf) 2: (%lf, %lf)->(%lf, %lf) 3:(%lf, %lf)->(%lf, %lf)\n", head(edge1)->x, head(edge1)->y, tail(edge1)->x, tail(edge1)->y,
                                                                                                     head(edge2)->x, head(edge2)->y, tail(edge2)->x, tail(edge2)->y,
                                                                                                     head(edge3)->x, head(edge3)->y, tail(edge3)->x, tail(edge3)->y);
        }
#endif
        PDASSERT(edge1->tail==edge2->head && edge2->tail==edge3->head && edge3->tail==edge1->head);
           
        t->v[0] = edge1->head;
        if (pt1->position_to_edge(pt2, pt3) == 1) {
            t->v[1] = edge2->head;
            t->v[2] = edge3->head;
            t->edge[0] = edge1;
            t->edge[1] = edge2;
            t->edge[2] = edge3;
        }
        else {
            fprintf(stderr, "not counterclockwise (%.20lf, %.20lf), (%.20lf, %.20lf), (%.20lf, %.20lf)\n", x_ref[pt1->id], y_ref[pt1->id], x_ref[pt2->id], y_ref[pt2->id], x_ref[pt3->id], y_ref[pt3->id]);
            fprintf(stderr, "not counterclockwise (%.20lf, %.20lf), (%.20lf, %.20lf), (%.20lf, %.20lf)\n", pt1->x, pt1->y, pt2->x, pt2->y, pt3->x, pt3->y);
            PDASSERT(false);
            t->v[1] = edge3->head;
            t->v[2] = edge2->head;
            if (edge1->twin_edge == NULL)
                edge1->twin_edge = make_twins_edge(edge1);
            if (edge2->twin_edge == NULL)
                edge2->twin_edge = make_twins_edge(edge2);
            if (edge3->twin_edge == NULL)
                edge3->twin_edge = make_twins_edge(edge3);
            t->edge[0] = edge3->twin_edge;
            t->edge[1] = edge2->twin_edge;
            t->edge[2] = edge1->twin_edge;
            assert(edge3->twin_edge != NULL && edge2->twin_edge != NULL && edge1->twin_edge != NULL);
        }

        /* recognizing cyclic triangle */
        //if (polar_mode) {
        //    int id[3];
        //    for (int j = 0; j < 3; j++)
        //        id[j] = vertex(t, j)->id;

        //    if (id[0] != -1 && id[1] != -1 && id[2] != -1)
        //        for (int j = 0; j < 3; j++) {
        //            if (calculate_distance(x_ref[id[j]], y_ref[id[j]], x_ref[id[(j+1)%3]], y_ref[id[(j+1)%3]]) > 180) {
        //                t->is_cyclic = true;
        //                break;
        //            }
        //        }
        //}
    }

    t->remained_points_head = -1;
    t->remained_points_tail = -1;
    
    t->edge[0]->next_edge_in_triangle = t->edge[1];
    t->edge[1]->next_edge_in_triangle = t->edge[2];
    t->edge[2]->next_edge_in_triangle = t->edge[0];
    t->edge[0]->prev_edge_in_triangle = t->edge[2];
    t->edge[1]->prev_edge_in_triangle = t->edge[0];
    t->edge[2]->prev_edge_in_triangle = t->edge[1];

    t->edge[0]->triangle = t;
    t->edge[1]->triangle = t;
    t->edge[2]->triangle = t;

    t->edge[0]->ref_inc();
    t->edge[1]->ref_inc();
    t->edge[2]->ref_inc();
    t->calulate_circum_circle(pt1, pt2, pt3);

}


void Delaunay_Voronoi::initialize_edge(Edge* e, int head, int tail)
{
    e->head = head;
    e->tail = tail;
    e->twin_edge = NULL;
    e->next_edge_in_triangle = NULL;
    e->prev_edge_in_triangle = NULL;
    e->ref_count = 0;
    e->triangle = NULL;
}


/*
 * Input : Point to be checked
 * Return:  1    point is in circum circle
 *          0    point is on circum circle
 *         -1    point is out of circum circle
 */
int Delaunay_Voronoi::circum_circle_contains_reliably(const Edge *edge, Point *p, double tolerance)
{
    int ret1 = edge->triangle->circum_circle_position_to(p, tolerance);
    int ret2 = edge->twin_edge->triangle->circum_circle_position_to(&all_points[edge->prev_edge_in_triangle->head], tolerance);
    int ret3 = edge->triangle->circum_radius2 > edge->twin_edge->triangle->circum_radius2 ? ret1 : ret2;
    //int ret3 = ret1;

    if (ret1 > 0 && ret2 > 0) {
        return 0;
    } else if (ret3 == 0 || ret3 == 2) {
        return 1;
    } else {
        return -1;
    }
}


/*
 * Input : Point to be checked
 * Return:  2  on/in
 *          1  on/out
 *          0  in
 *         -1  out
 */
int Triangle::circum_circle_position_to(Point *p, double tolerance)
{
    PAT_REAL dist2 = (((PAT_REAL)p->x - circum_center[0]) * ((PAT_REAL)p->x - circum_center[0])) + (((PAT_REAL)p->y - circum_center[1]) * ((PAT_REAL)p->y - circum_center[1]));

    int ret = dist2 > circum_radius2 ? -1 : 0;

    if(relative_eq_int(dist2, circum_radius2, tolerance))
        ret += 2;

    return ret;
}


int Triangle::find_best_candidate_point(Point* buf)
{
    if (remained_points_tail == -1)
        return -1;

    double min_dist=1e10, dist;
    double center_x = (buf[v[0]].x+buf[v[1]].x+buf[v[2]].x) * 0.3333333333333333333333333333;
    double center_y = (buf[v[0]].y+buf[v[1]].y+buf[v[2]].y) * 0.3333333333333333333333333333;
    int best_candidate_id = -1;
    bool no_more_mask = true;
    for (int i = remained_points_head; i > -1; i = buf[i].next) {
        if (buf[i].mask)
            no_more_mask = false;
        else
            continue;

        dist = buf[i].calculate_distance(center_x, center_y);
        if (min_dist > dist) {
            min_dist = dist;
            best_candidate_id = i;
        }
    }

    if (no_more_mask) {
        remained_points_head = remained_points_tail = -1;
        return -1;
    }

    return best_candidate_id;
}


bool Triangle::contain_vertex(int pt)
{
    for(int i = 0; i < 3; i++)
        if(v[i] == pt)
            return true;
    return false;
}


void Delaunay_Voronoi::swap_points(int idx1, int idx2)
{
    double tmp_x    = all_points[idx1].x;
    double tmp_y    = all_points[idx1].y;
    int    tmp_id   = all_points[idx1].id;
    bool   tmp_mask = all_points[idx1].mask;

    all_points[idx1].x    = all_points[idx2].x;
    all_points[idx1].y    = all_points[idx2].y;
    all_points[idx1].id   = all_points[idx2].id;
    all_points[idx1].mask = all_points[idx2].mask;

    all_points[idx2].x    = tmp_x;
    all_points[idx2].y    = tmp_y;
    all_points[idx2].id   = tmp_id;
    all_points[idx2].mask = tmp_mask;
}


/*
 *   head:  head of linked list of points.
 *   tail:  tail of linked list of points.
 *   base:  base of the triangle stack.
 *   top:   top of the triangle stack.
 */
void Delaunay_Voronoi::distribute_points_into_triangles(int head, int tail, unsigned base, unsigned top)
{
    int start = head;
    int end, j;

    if (tail == -1)
        return;

    PDASSERT(head != -1);

    for (unsigned i = base+1; i <= top; i++) {
        if (!triangle_stack[i] || !triangle_stack[i]->is_leaf)
            continue;

        end = tail;

        /* put points in current triangle together */
        for (j = start; j != end && j > -1;) {
            int* v_idx = triangle_stack[i]->v;
            if (all_points[j].position_to_triangle(&all_points[v_idx[0]], &all_points[v_idx[1]], &all_points[v_idx[2]]) >= 0) { // in triangles or on the edge
                j = all_points[j].next;
            } else {
                swap_points(j, end);
                end = all_points[end].prev;
            }
        }

        int* v_idx = triangle_stack[i]->v;
        if (j > -1 && all_points[end].position_to_triangle(&all_points[v_idx[0]], &all_points[v_idx[1]], &all_points[v_idx[2]]) < 0) // Case "j == end"
            end = all_points[end].prev;

        /* set triangle remained points info */
        triangle_stack[i]->set_remained_points(start, end);

        /* go on */
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


int Triangle::pop_tail(Point* buf)
{
    PDASSERT(remained_points_tail != -1);
    int old_tail = remained_points_tail;

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
    t->stack_ref_count++;
}


inline Triangle_inline Delaunay_Voronoi::pack_triangle(Triangle* t)
{
    if (fast_mode) {
        Triangle_inline a = Triangle_inline(Point(vertex(t, 0)->x, vertex(t, 0)->y, global_index[vertex(t, 0)->id]),
                                            Point(vertex(t, 1)->x, vertex(t, 1)->y, global_index[vertex(t, 1)->id]),
                                            Point(vertex(t, 2)->x, vertex(t, 2)->y, global_index[vertex(t, 2)->id]),
                                            false);
        a.check_cyclic();
        return a;
    } else if (polar_mode) {
        Triangle_inline a = Triangle_inline(Point(x_ref[vertex(t, 0)->id], y_ref[vertex(t, 0)->id], global_index[vertex(t, 0)->id]),
                                            Point(x_ref[vertex(t, 1)->id], y_ref[vertex(t, 1)->id], global_index[vertex(t, 1)->id]),
                                            Point(x_ref[vertex(t, 2)->id], y_ref[vertex(t, 2)->id], global_index[vertex(t, 2)->id]),
                                            false);
        a.check_cyclic();
        return a;
    } else if (x_ref) {
        return Triangle_inline(Point(x_ref[vertex(t, 0)->id], y_ref[vertex(t, 0)->id], global_index[vertex(t, 0)->id]),
                               Point(x_ref[vertex(t, 1)->id], y_ref[vertex(t, 1)->id], global_index[vertex(t, 1)->id]),
                               Point(x_ref[vertex(t, 2)->id], y_ref[vertex(t, 2)->id], global_index[vertex(t, 2)->id]),
                               false);
    } else {
        return Triangle_inline(Point(vertex(t, 0)->x, vertex(t, 0)->y, global_index[vertex(t, 0)->id]),
                               Point(vertex(t, 1)->x, vertex(t, 1)->y, global_index[vertex(t, 1)->id]),
                               Point(vertex(t, 2)->x, vertex(t, 2)->y, global_index[vertex(t, 2)->id]),
                               false);
    }
}


static inline double point_distence_to_line(double px, double py, double x1, double y1, double x2, double y2)
{
    //TODO: may be slow
    if(x1 == x2)
        return fabs(px - x1);
    else if(y1 == y2)
        return fabs(py - y1);
    else {
        double A=(y1-y2)/(x1-x2);
        double B=y1-A*x1;
        return fabs((A*px+B-py)/sqrt(A*A+1));
    }
}

static inline bool point_distence_in_threshold(Point v, Point p1, Point p2, double threshold) {
    return point_distence_to_line(v.x, v.y, p1.x, p1.y, p2.x, p2.y) <= threshold;
}

static inline bool all_distence_in_threshold(Triangle_inline* triangle, Point p1, Point p2, double threshold) {
    return point_distence_in_threshold(triangle->v[0], p1, p2, threshold) &&
           point_distence_in_threshold(triangle->v[1], p1, p2, threshold) &&
           point_distence_in_threshold(triangle->v[2], p1, p2, threshold);
}

static inline bool is_segment_in_triangle(Triangle_inline* triangle, Point p1, Point p2) {
    return p1.position_to_triangle(triangle) >= 0 && p2.position_to_triangle(triangle) >= 0;
}

static inline bool is_segment_intersected_with_edge(Point p_e1, Point p_e2, Point p1, Point p2) {
    return (p1.position_to_edge(&p_e1, &p_e2) * p2.position_to_edge(&p_e1, &p_e2) < 0) &&
           (p_e1.position_to_edge(&p1, &p2) * p_e2.position_to_edge(&p1, &p2) <= 0);
}

static inline bool is_segment_intersected_with_any_edges(Triangle_inline* triangle, Point p1, Point p2) {
    return is_segment_intersected_with_edge(triangle->v[0], triangle->v[1], p1, p2) ||
           is_segment_intersected_with_edge(triangle->v[1], triangle->v[2], p1, p2) ||
           is_segment_intersected_with_edge(triangle->v[2], triangle->v[0], p1, p2);
}

inline void get_bounding_box(Triangle_inline* t, double& x_min, double& x_max, double& y_min, double& y_max) {
    Point* v = t->v;
    x_min = v[0].x;
    x_max = v[0].x;
    y_min = v[0].y;
    y_max = v[0].y;
    for (int i = 1; i < 3; i++) {
        if (v[i].x < x_min) x_min = v[i].x;
        if (v[i].x > x_max) x_max = v[i].x;
        if (v[i].y < y_min) y_min = v[i].y;
        if (v[i].y > y_max) y_max = v[i].y;
    }
}

inline bool is_triangle_intersecting_with_segment(Triangle_inline* triangle, Point p1, Point p2, double threshold) {
    PDASSERT(p1.x == p2.x || p1.y == p2.y);
    double x_min, x_max, y_min, y_max;
    get_bounding_box(triangle, x_min, x_max, y_min, y_max);

    double seg_min, seg_max;
    if (p1.x == p2.x) {
        if (p1.y > p2.y) {
            seg_min = p2.y;
            seg_max = p1.y;
        } else {
            seg_min = p1.y;
            seg_max = p2.y;
        }

        if(!(x_min <= p1.x && x_max >= p1.x && !(y_max < seg_min || y_min > seg_max)))
            return false;
    } else {
        if (p1.x > p2.x) {
            seg_min = p2.x;
            seg_max = p1.x;
        } else {
            seg_min = p1.x;
            seg_max = p2.x;
        }

        if(!(y_min <= p1.y && y_max >= p1.y && !(x_max < seg_min || x_min > seg_max)))
            return false;
    }
    return //(threshold == 0 || all_distence_in_threshold(triangle, p1, p2, threshold)) &&
           (is_segment_intersected_with_any_edges(triangle, p1, p2) || is_segment_in_triangle(triangle, p1, p2));
}


unsigned Delaunay_Voronoi::triangulating_process(Triangle *triangle, unsigned stack_base)
{
    unsigned stack_top = stack_base;

#ifdef DEBUG
    PDASSERT(triangle->is_leaf);
#endif
    int candidate_id = triangle->find_best_candidate_point(&all_points[0]);

    if (candidate_id == -1) {
        if (polar_mode) {
            int id[3];
            for (int j = 0; j < 3; j++)
                id[j] = vertex(triangle, j)->id;

            if (id[0] != -1 && id[1] != -1 && id[2] != -1)
                for (int j = 0; j < 3; j++) {
                    if (calculate_distance(x_ref[id[j]], y_ref[id[j]], x_ref[id[(j+1)%3]], y_ref[id[(j+1)%3]]) > PAT_CYCLIC_EDGE_THRESHOLD) {
                        triangle->is_cyclic = true;
                        break;
                    }
                }

            if (triangle->is_cyclic) {
                for(int j = 0; j < 3; j++)
                    relegalize_triangles(triangle->v[j], triangle->edge[(j+1)%3]);
            }
        }

        return stack_top;
    }

    triangle->is_leaf = false;
    swap_points(candidate_id, triangle->remained_points_tail);
    int dividing_idx = triangle->pop_tail(&all_points[0]);
    Point* dividing_point = &all_points[dividing_idx];

    int* v_idx = triangle->v;
    if (dividing_point->position_to_triangle(&all_points[v_idx[0]], &all_points[v_idx[1]], &all_points[v_idx[2]]) == 0) { // inside
        Edge *e_v1_can = allocate_edge(triangle->v[0], dividing_idx);
        Edge *e_can_v1 = make_twins_edge(e_v1_can);
        Edge *e_v2_can = allocate_edge(triangle->v[1], dividing_idx);
        Edge *e_can_v2 = make_twins_edge(e_v2_can);
        Edge *e_v3_can = allocate_edge(triangle->v[2], dividing_idx);
        Edge *e_can_v3 = make_twins_edge(e_v3_can);
        Triangle *t_can_v1_v2 = allocate_triangle(e_can_v1, triangle->edge[0], e_v2_can);
        Triangle *t_can_v2_v3 = allocate_triangle(e_can_v2, triangle->edge[1], e_v3_can);
        Triangle *t_can_v3_v1 = allocate_triangle(e_can_v3, triangle->edge[2], e_v1_can);
        push(&stack_top, triangle);
        push(&stack_top, t_can_v1_v2);
        push(&stack_top, t_can_v2_v3);
        push(&stack_top, t_can_v3_v1);
        
        /* Actually, vertex(t_can_v1_v2, 0) is dividing_point */
        legalize_triangles(t_can_v1_v2->v[0], t_can_v1_v2->edge[1], stack_base, &stack_top);
        legalize_triangles(t_can_v2_v3->v[0], t_can_v2_v3->edge[1], stack_base, &stack_top);
        legalize_triangles(t_can_v3_v1->v[0], t_can_v3_v1->edge[1], stack_base, &stack_top);
    } else { // on the side
        int idx_i, idx_j, idx_k;
        Edge *eij;
        switch (dividing_point->position_to_triangle(&all_points[v_idx[0]], &all_points[v_idx[1]], &all_points[v_idx[2]])) {
            case 1:
                idx_i = triangle->v[0];
                idx_j = triangle->v[1];
                idx_k = triangle->v[2];
                eij = triangle->edge[0];
                break;
            case 2:
                idx_i = triangle->v[1];
                idx_j = triangle->v[2];
                idx_k = triangle->v[0];
                eij = triangle->edge[1];
                break;
            case 3:
                idx_i = triangle->v[2];
                idx_j = triangle->v[0];
                idx_k = triangle->v[1];
                eij = triangle->edge[2];
                break;
            default:
                printf("point, which should be found in triangle, is outside of triangle\n");
                PDASSERT(false);
                break;
        }
        Point* vi = &all_points[idx_i];
        Point* vj = &all_points[idx_j];
        PDASSERT(dividing_point->position_to_edge(vi, vj) == 0);
        PDASSERT(eij->twin_edge == NULL || eij->twin_edge->triangle->is_leaf);

        Edge* ejk = eij->next_edge_in_triangle;
        Edge* eki = ejk->next_edge_in_triangle;
        Edge *eil, *elj, *eji;

        Edge *eir = allocate_edge(idx_i, dividing_idx);
        Edge *erk = allocate_edge(dividing_idx, idx_k);
        Edge *ekr = make_twins_edge(erk);
        Edge *erj = allocate_edge(dividing_idx, idx_j);
        Triangle* tirk = allocate_triangle(eir, erk, eki);
        Triangle* tjkr = allocate_triangle(ejk, ekr, erj);
        push(&stack_top, triangle);
        push(&stack_top, tirk);
        push(&stack_top, tjkr);
        legalize_triangles(dividing_idx, ejk, stack_base, &stack_top);
        legalize_triangles(dividing_idx, eki, stack_base, &stack_top); 
        if (eij->twin_edge != NULL) {
            eji = eij->twin_edge;
            eij->twin_edge->triangle->is_leaf = false;
            eil = eji->next_edge_in_triangle;
            elj = eil->next_edge_in_triangle;
            int idx_l = elj->head;
            Edge *eri = make_twins_edge(eir);
            Edge *ejr = make_twins_edge(erj);
            Edge *erl = allocate_edge(dividing_idx, idx_l);
            Edge *elr = make_twins_edge(erl);
            Triangle* tilr = allocate_triangle(eil, elr, eri);
            Triangle* tjrl = allocate_triangle(ejr, erl, elj);
            push(&stack_top, eij->twin_edge->triangle);
            push(&stack_top, tilr);
            push(&stack_top, tjrl);
            legalize_triangles(dividing_idx, eil, stack_base, &stack_top);
            legalize_triangles(dividing_idx, elj, stack_base, &stack_top);
        }
    }

#ifdef DEBUG
    for (unsigned i = stack_base+1; i <= stack_top; i ++)
        if (triangle_stack[i] && triangle_stack[i]->is_leaf)
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

    //printf("base: %d, top: %d\n", stack_base, stack_top);
    for (unsigned i = stack_base+1; i <= stack_top; i ++)
        if(triangle_stack[i] && triangle_stack[i]->is_leaf)
            triangulating_process(triangle_stack[i], stack_top);

    Triangle* killed;
    for (unsigned i = stack_base+1; i <= stack_top; i ++) {
        triangle_stack[i]->stack_ref_count--;
        if(triangle_stack[i]->stack_ref_count <= 0 && !triangle_stack[i]->is_leaf) {
            killed = triangle_stack[i];
            //delete killed;
            clean_triangle(killed);
            triangle_allocator.deleteElement(killed);
        }
    }

    return stack_top;
}

void Delaunay_Voronoi::clean_triangle(Triangle* t)
{
    for (int i = 0; i < 3; i++)
        ref_dec(t->edge[i]);
}


static int compare_node_index(const void* a, const void* b)
{
    if(*(const int*)a < *(const int*)b) return -1;
    if(*(const int*)a > *(const int*)b) return 1;
    return 0;
}


void Delaunay_Voronoi::link_remained_list(unsigned base, unsigned top, int* head, int* tail)
{
    unsigned count;
    unsigned i;
    int* head_tail = (int*) alloca(sizeof(int) * (top - base) * 2); // [head1, tail1, head2, tail2, ...]

    for (i = base+1, count = 0; i <= top; i ++)
        if (triangle_stack[i] && !triangle_stack[i]->is_leaf && triangle_stack[i]->remained_points_tail > -1) {
            head_tail[count * 2] = triangle_stack[i]->remained_points_head;
            head_tail[count * 2 + 1] = triangle_stack[i]->remained_points_tail;
            count++;
        }

#ifdef DEBUG
    bool* map = new bool[num_points]();
    for (i = base+1; i <= top; i ++)
        if (triangle_stack[i] && !triangle_stack[i]->is_leaf && triangle_stack[i]->remained_points_tail > -1) {
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


void Delaunay_Voronoi::map_buffer_index_to_point_index()
{
    delete[] point_idx_to_buf_idx;
    point_idx_to_buf_idx = new int[num_points - PAT_NUM_LOCAL_VPOINTS]();
    for (int i = PAT_NUM_LOCAL_VPOINTS; i < num_points; i++) {
#ifdef DEBUG
        PDASSERT(point_idx_to_buf_idx[all_points[i].id] == 0);
#endif
        point_idx_to_buf_idx[all_points[i].id] = i;
    }

#ifdef DEBUG
    //for (int i = 0; i < num_points; i++) {
    //    PDASSERT(x_store[i] == all_points[point_idx_to_buf_idx[i]].x);
    //    PDASSERT(y_store[i] == all_points[point_idx_to_buf_idx[i]].y);
    //}
#endif
}


void Delaunay_Voronoi::mark_special_triangles()
{
#ifdef DEBUG
    triangles_containing_vpolar.clear();
#endif
    for (unsigned i = 0; i < all_leaf_triangles.size(); i++) {
        if(!all_leaf_triangles[i]->is_leaf)
            continue;

        for(unsigned j = 0; j < PAT_NUM_LOCAL_VPOINTS; j++)
            if(all_leaf_triangles[i]->contain_vertex(j)) {
                all_leaf_triangles[i]->is_virtual = true;
                break;
            }

#ifdef DEBUG
        if(!all_leaf_triangles[i]->is_virtual &&
           all_leaf_triangles[i]->contain_vertex(point_idx_to_buf_idx[vpolar_local_index]))
            triangles_containing_vpolar.push_back(all_leaf_triangles[i]);
#endif
    }
}


Delaunay_Voronoi::Delaunay_Voronoi()
    : triangle_stack(NULL)
    , stack_size(0)
    , polar_mode(false)
    , fast_mode(false)
    , tolerance(1e-9)
    , num_points(0)
    , vpolar_local_index(-1)
    , x_ref(NULL)
    , y_ref(NULL)
    , global_index(NULL)
    , point_idx_to_buf_idx(NULL)
    , have_bound(false)
{
    avoiding_line_head[0].x = avoiding_line_head[0].y = 0;
    avoiding_line_head[1].x = avoiding_line_head[1].y = 0;
    avoiding_line_tail[0].x = avoiding_line_tail[0].y = 0;
    avoiding_line_tail[1].x = avoiding_line_tail[1].y = 0;
    avoiding_line_head[0].id = avoiding_line_head[1].id = 0;

    avoiding_circle_center[0].x = avoiding_circle_center[0].y = 0;
    avoiding_circle_center[1].x = avoiding_circle_center[1].y = 0;
    avoiding_circle_center[0].id = avoiding_circle_center[1].id = 0;
    avoiding_circle_radius[0] = avoiding_circle_radius[1] = 0.0;
}


Delaunay_Voronoi::~Delaunay_Voronoi()
{
    delete[] global_index;
    delete[] point_idx_to_buf_idx;
    delete[] triangle_stack;
    delete[] x_ref;
    delete[] y_ref;
    for (unsigned i = 0; i < all_leaf_triangles.size(); i ++) {
        clean_triangle(all_leaf_triangles[i]);
        triangle_allocator.deleteElement(all_leaf_triangles[i]);
    }
}


inline bool Delaunay_Voronoi::point_in_triangle(double x, double y, Triangle* t)
{
    Point p(x, y);
    int* v_idx = t->v;
    int ret = p.position_to_triangle(&all_points[v_idx[0]], &all_points[v_idx[1]], &all_points[v_idx[2]]);
    return ret != -1;
}


struct Bound {
    double min_x;
    double max_x;
    double min_y;
    double max_y;
};


inline bool Delaunay_Voronoi::point_in_bound(double x, double y, Bound* b)
{
    return x >= b->min_x && x <= b->max_x && y >= b->min_y && y <= b->max_y;
}


Bound* Delaunay_Voronoi::make_bounding_box()
{
    unsigned num_triangles = all_leaf_triangles.size();
    Bound* bound = new Bound[num_triangles];

    for (unsigned i = 0; i < num_triangles; i++){
        Triangle* lf = all_leaf_triangles[i];
        double x[3], y[3];
        for (int j = 0; j < 3; j++) {
            x[j] = vertex(lf, j)->x;
            y[j] = vertex(lf, j)->y;
        }
        if (x[0] > x[1]) std::swap(x[0], x[1]);
        if (x[1] > x[2]) std::swap(x[1], x[2]);
        if (x[0] > x[1]) std::swap(x[0], x[1]);
        if (y[0] > y[1]) std::swap(y[0], y[1]);
        if (y[1] > y[2]) std::swap(y[1], y[2]);
        if (y[0] > y[1]) std::swap(y[0], y[1]);

        bound[i].min_x = x[0];
        bound[i].max_x = x[2];
        bound[i].min_y = y[0];
        bound[i].max_y = y[2];
    }

    return bound;
}


static inline unsigned hash(double x, double y, double block_size, double min_x, double len_x, double min_y, double len_y)
{
    unsigned x_idx = (x-min_x) / len_x * block_size;
    unsigned y_idx = (y-min_y) / len_y * block_size;
    PDASSERT(x_idx < block_size && y_idx < block_size);
    return x_idx+y_idx*block_size;
}


static inline void hash(const Bound* bound, double block_size, double min_x, double len_x, double min_y, double len_y,
                   unsigned* x_idx_bgn, unsigned* x_idx_end, unsigned* y_idx_bgn, unsigned* y_idx_end)
{
    *x_idx_bgn = (bound->min_x-min_x) / len_x * block_size; // FIXME: optimise
    *y_idx_bgn = (bound->min_y-min_y) / len_y * block_size;
    *x_idx_end = (bound->max_x-min_x) / len_x * block_size;
    *y_idx_end = (bound->max_y-min_y) / len_y * block_size;
    if (*x_idx_end == block_size) (*x_idx_end)--;
    if (*y_idx_end == block_size) (*y_idx_end)--;
}


#define PAT_TRIANGLES_PER_BLOCK (50)
#define PAT_MAX_BLOCK (10000)
void Delaunay_Voronoi::distribute_initial_points(const double* x, const double* y, int num, int** output_nexts)
{
    double min_x = all_points[0].x;
    double max_y = all_points[0].y;
    double max_x = all_points[2].x;
    double min_y = all_points[2].y;

    double len_x = max_x - min_x;
    double len_y = max_y - min_y;

    unsigned num_triangles = all_leaf_triangles.size();
    int* nexts = new int[num];

    memset(nexts, -1, num*sizeof(int));

    unsigned block_size = std::sqrt(std::min(num_triangles / PAT_TRIANGLES_PER_BLOCK, (unsigned)PAT_MAX_BLOCK));

    if (block_size*block_size > 2) {
        Bound* bound = make_bounding_box();

        for (unsigned i = 0; i < num_triangles; i++) {
            PDASSERT(bound[i].min_x >= min_x);
            PDASSERT(bound[i].max_x <= max_x);
            PDASSERT(bound[i].min_y >= min_y);
            PDASSERT(bound[i].max_y <= max_y);
        }

        unsigned* num_including_points = new unsigned[block_size*block_size]();
        unsigned* num_including_triangles = new unsigned[block_size*block_size]();

        /* first scan: counting */
        for (int i = 0; i < num; i++) {
            unsigned idx = hash(x[i], y[i], block_size, min_x, len_x, min_y, len_y);
            num_including_points[idx]++;
        }

        for (unsigned i = 0; i < num_triangles; i++) {
            unsigned x_idx_bgn, y_idx_bgn, x_idx_end, y_idx_end;
            hash(&bound[i], block_size, min_x, len_x, min_y, len_y, &x_idx_bgn, &x_idx_end, &y_idx_bgn, &y_idx_end);
            for (unsigned j = y_idx_bgn; j <= y_idx_end; j++)
                for (unsigned k = x_idx_bgn; k <= x_idx_end; k++) {
                    unsigned idx = k + j * block_size;
                    PDASSERT(j < block_size);
                    PDASSERT(k < block_size);
                    PDASSERT(idx < block_size*block_size);
                    num_including_triangles[idx]++;
                }
        }

        /* allocing memory */
        unsigned** including_points = new unsigned*[block_size*block_size];
        unsigned** including_triangles = new unsigned*[block_size*block_size];
        for (unsigned i = 0; i < block_size*block_size; i++) {
            if (num_including_points[i] > 0) {
                including_points[i] = new unsigned[num_including_points[i]];
                including_triangles[i] = new unsigned[num_including_triangles[i]];
            } else {
                including_points[i] = NULL;
                including_triangles[i] = NULL;
            }
        }

        /* second scan: distributing points and triangles into mesh */
        memset(num_including_points, 0, sizeof(unsigned)*block_size*block_size);
        memset(num_including_triangles, 0, sizeof(unsigned)*block_size*block_size);

        for (int i = 0; i < num; i++) {
            unsigned idx = hash(x[i], y[i], block_size, min_x, len_x, min_y, len_y);
            including_points[idx][num_including_points[idx]] = i;
            num_including_points[idx]++;
        }

        for (unsigned i = 0; i < num_triangles; i++) {
            unsigned x_idx_bgn, y_idx_bgn, x_idx_end, y_idx_end;
            hash(&bound[i], block_size, min_x, len_x, min_y, len_y, &x_idx_bgn, &x_idx_end, &y_idx_bgn, &y_idx_end);
            for (unsigned j = y_idx_bgn; j <= y_idx_end; j++)
                for (unsigned k = x_idx_bgn; k <= x_idx_end; k++) {
                    unsigned idx = k + j * block_size;
                    if (including_triangles[idx]) {
                        including_triangles[idx][num_including_triangles[idx]] = i;
                        num_including_triangles[idx]++;
                    }
                }
        }

        /* distributing points into triangles in the same mesh */
        for (unsigned m = 0; m < block_size*block_size; m++) {
            if (including_points[m] == NULL)
                continue;

            unsigned* p_idxs = including_points[m];
            unsigned* t_idxs = including_triangles[m];
            unsigned  p_num = num_including_points[m];
            unsigned  t_num = num_including_triangles[m];

            for (unsigned i = 0; i < p_num; i++) {
                unsigned p_idx = p_idxs[i];
                for (unsigned j = 0; j < t_num; j++) {
                    unsigned t_idx = t_idxs[j];
                    if (!all_leaf_triangles[t_idx]->is_leaf)
                        continue;

                    if (!point_in_bound(x[p_idx], y[p_idx], &bound[t_idx]))
                        continue;

                    if (point_in_triangle(x[p_idx], y[p_idx], all_leaf_triangles[t_idx])) {
                        if (all_leaf_triangles[t_idx]->remained_points_head == -1)
                            all_leaf_triangles[t_idx]->remained_points_head = all_leaf_triangles[t_idx]->remained_points_tail = p_idx;
                        else
                            all_leaf_triangles[t_idx]->remained_points_tail = nexts[all_leaf_triangles[t_idx]->remained_points_tail] = p_idx;
                        break;
                    }
                }
            }
        }

        /* freeing memory */
        for (unsigned i = 0; i < block_size*block_size; i++) {
            delete[] including_points[i];
            delete[] including_triangles[i];
        }
        delete[] including_points;
        delete[] including_triangles;

        delete[] num_including_points;
        delete[] num_including_triangles;
        delete[] bound;
    } else {
        for (int i = 0; i < num; i++) {
            for (unsigned j = 0; j < num_triangles; j++) {
               if (!all_leaf_triangles[j]->is_leaf)
                   continue;

               if (point_in_triangle(x[i], y[i], all_leaf_triangles[j])) {
                   if (all_leaf_triangles[j]->remained_points_head == -1)
                       all_leaf_triangles[j]->remained_points_head = all_leaf_triangles[j]->remained_points_tail = i;
                   else
                       all_leaf_triangles[j]->remained_points_tail = nexts[all_leaf_triangles[j]->remained_points_tail] = i;
                   break;
               }
            }
        }
    }

    *output_nexts = nexts;
}


void Delaunay_Voronoi::enlarge_super_rectangle(const double* x, const double* y, int num)
{
    double min_x = 1e10;
    double max_x = -1e10;
    double min_y = 1e10;
    double max_y = -1e10;

    for (int i = 0; i < num; i++) {
        if (x[i] < min_x) min_x = x[i];
        if (x[i] > max_x) max_x = x[i];
        if (y[i] < min_y) min_y = y[i];
        if (y[i] > max_y) max_y = y[i];
    }

    const double ratio  = 0.1;
    double dx = max_x - min_x;
    double dy = max_y - min_y;
    double delta = std::max(dx, dy) * ratio;

    min_x -= delta;
    max_x += delta;
    min_y -= delta;
    max_y += delta;

    all_points[0].x = std::min(min_x, all_points[0].x);
    all_points[0].y = std::max(max_y, all_points[0].y);
    all_points[1].x = std::min(min_x, all_points[1].x);
    all_points[1].y = std::min(min_y, all_points[1].y);
    all_points[2].x = std::max(max_x, all_points[2].x);
    all_points[2].y = std::min(min_y, all_points[2].y);
    all_points[3].x = std::max(max_x, all_points[3].x);
    all_points[3].y = std::max(max_y, all_points[3].y);
}


void Delaunay_Voronoi::initialize(int num)
{
    int triangles_count_estimate = 2*(num+PAT_NUM_LOCAL_VPOINTS);
    all_leaf_triangles.reserve(triangles_count_estimate);

    stack_size = triangles_count_estimate * 2;
    triangle_stack = new Triangle*[stack_size];

    max_points = num*3/2;
    all_points = (Point *)::operator new(max_points * sizeof(Point));

    new(&all_points[0]) Point(-0.1,  0.1, -1, -1, -1);
    new(&all_points[1]) Point(-0.1, -0.1, -1, -1, -1);
    new(&all_points[2]) Point( 0.1, -0.1, -1, -1, -1);
    new(&all_points[3]) Point( 0.1,  0.1, -1, -1, -1);

    num_points = PAT_NUM_LOCAL_VPOINTS;

    all_leaf_triangles.push_back(allocate_triangle(allocate_edge(0, 1), allocate_edge(1, 2), allocate_edge(2, 0)));
    all_leaf_triangles.push_back(allocate_triangle(allocate_edge(0, 2), allocate_edge(2, 3), allocate_edge(3, 0)));

    Edge* e1 = all_leaf_triangles[0]->edge[2];
    Edge* e2 = all_leaf_triangles[1]->edge[0];
    e1->twin_edge = e2;
    e2->twin_edge = e1;
}


void Delaunay_Voronoi::extend_points_buffer(int introduced_points)
{
    int full_size = num_points + introduced_points;
    int realloc_size = full_size + introduced_points * PDLN_EXPECTED_EXPANDING_TIMES;

    void* tmp_buf = ::operator new(realloc_size * sizeof(Point));
    memcpy(tmp_buf, all_points, num_points * sizeof(Point));

    ::operator delete(all_points);
    all_points = (Point *)tmp_buf;
    max_points = realloc_size;
}


void Delaunay_Voronoi::add_points(const double* x, const double* y, const bool* mask, int num)
{
#ifdef DEBUG
    PDASSERT(have_redundent_points(x, y, num) == false);
#endif
    PDASSERT(x != NULL);
    PDASSERT(y != NULL);

    if(stack_size == 0)
        initialize(num);
    else if (num_points + num > max_points)
        extend_points_buffer(num);

    enlarge_super_rectangle(x, y, num);

    int *nexts;

    distribute_initial_points(x, y, num, &nexts);

    dirty_triangles_count = 0;
    int buf_idx_cur = num_points;
    int local_idx_start = num_points - PAT_NUM_LOCAL_VPOINTS;
    for (unsigned i = 0; i < all_leaf_triangles.size(); i++) {
        int head = buf_idx_cur;
        for (int p = all_leaf_triangles[i]->remained_points_head; p != -1; p = nexts[p]) {
            new(&all_points[buf_idx_cur]) Point(x[p], y[p], local_idx_start+p, mask ? mask[p] : true, buf_idx_cur+1, buf_idx_cur-1);
            buf_idx_cur++;
        }

        if (head != buf_idx_cur) {
            all_points[head].prev = -1;
            all_points[buf_idx_cur-1].next = -1;
            all_leaf_triangles[i]->set_remained_points(head, buf_idx_cur-1);
            push(&dirty_triangles_count, all_leaf_triangles[i]);
        }
    }
    num_points += num;

#ifdef DEBUG
    x_store = x;
    y_store = y;

    bool* check = new bool[num_points]();
    for (int i = 4; i < num_points; i++) {
        PDASSERT(check[all_points[i].id] == 0);
        check[all_points[i].id] = 1;
    }
#endif

    delete[] nexts;
}


void Delaunay_Voronoi::map_global_index(const int* global_idx)
{
    global_index = global_idx;
}


void Delaunay_Voronoi::set_virtual_polar_index(int idx)
{
    if (idx == -1)
        return;

    vpolar_local_index = idx + PAT_NUM_LOCAL_VPOINTS;
}


void Delaunay_Voronoi::set_origin_coord(const double *x_origin, const double *y_origin, int num)
{
    x_ref = x_origin;
    y_ref = y_origin;

    if (num_points == 0)
        num_points = num;
    else
        PDASSERT(num + PAT_NUM_LOCAL_VPOINTS == num_points);
}


#define PAT_COORDHASH_FACTOR (2)
bool Delaunay_Voronoi::try_fast_triangulate(double min_x, double max_x, double min_y, double max_y)
{
    double area = (max_x - min_x) * (max_y - min_y);
    double density = std::sqrt(num_points / area);

    int lon_expected_points = density * (max_x - min_x);
    int lat_expected_points = density * (max_y - min_y);

    Coord_Hash lon_table(lon_expected_points * PAT_COORDHASH_FACTOR);
    Coord_Hash lat_table(lat_expected_points * PAT_COORDHASH_FACTOR);

    lon_table.set_hashing_params(min_x, max_x);
    lat_table.set_hashing_params(min_y, max_y);

    bool have_vpolar = vpolar_local_index != -1;
    int vpolar_index = -1;
    if (have_vpolar)
        vpolar_index = vpolar_local_index - PAT_NUM_LOCAL_VPOINTS;

    for (int i = 0; i < num_points; i++) {
        if (vpolar_index != -1 && vpolar_index == i)
            continue;

        lon_table.put(x_ref[i]);
        lat_table.put(y_ref[i]);
        if (lon_table.get_num_unique_values() * lat_table.get_num_unique_values() > num_points)
            return false;
    }

    int num_lon = lon_table.get_num_unique_values();
    int num_lat = lat_table.get_num_unique_values();

    if (!have_vpolar && num_lon * num_lat != num_points)
        return false;
    if (have_vpolar && num_lon * num_lat != num_points-1)
        return false;

    fast_mode = true;

    lon_table.make_sorted_index();
    lat_table.make_sorted_index();

    if (!have_vpolar)
        all_points = (Point *)::operator new(num_lon * num_lat * sizeof(Point));
    else
        all_points = (Point *)::operator new((num_lon * num_lat + 1) * sizeof(Point));

    for (int i = 0; i < num_points; i++) {
        if (vpolar_index != -1 && vpolar_index == i)
            continue;

        double lon = x_ref[i];
        double lat = y_ref[i];
        int buf_idx = lon_table.get_index(lon) + lat_table.get_index(lat) * num_lon;
        new(&all_points[buf_idx]) Point(lon, lat, buf_idx, -1, -1);
    }

    /* put vpolar at last space */
    if (have_vpolar)
        new(&all_points[num_lon * num_lat]) Point(x_ref[vpolar_index], y_ref[vpolar_index], num_lon * num_lat, -1, -1);

    fast_triangulate(num_lon, num_lat, have_vpolar);

    return true; 
}


void Delaunay_Voronoi::fast_triangulate(int num_lon, int num_lat, bool have_vpolar)
{
    if (have_vpolar)
        PDASSERT(num_lon * num_lat + 1 == num_points);
    else
        PDASSERT(num_lon * num_lat == num_points);

    for (int j = 0; j < num_lat - 1; j++) {
        for (int i = 0; i < num_lon - 1; i++) {
            all_leaf_triangles.push_back(allocate_triangle(j*num_lon+i,     j*num_lon+i+1, (j+1)*num_lon+i));
            all_leaf_triangles.push_back(allocate_triangle((j+1)*num_lon+i, j*num_lon+i+1, (j+1)*num_lon+i+1));
        }
        int i = num_lon - 1;
        all_leaf_triangles.push_back(allocate_triangle(j*num_lon+i,     j*num_lon, (j+1)*num_lon+i, true));
        all_leaf_triangles.push_back(allocate_triangle((j+1)*num_lon+i, j*num_lon, (j+1)*num_lon,   true));
        all_leaf_triangles[all_leaf_triangles.size()-1]->is_cyclic = true;
        all_leaf_triangles[all_leaf_triangles.size()-2]->is_cyclic = true;
    }

    if (have_vpolar) {
        int vpole_idx = num_points-1;
        if (float_eq(all_points[vpole_idx].y, -90)) {
            /* south pole */
            for (int i = 0; i < num_lon - 1; i++)
                all_leaf_triangles.push_back(allocate_triangle(vpole_idx, i+1, i));
            all_leaf_triangles.push_back(allocate_triangle(vpole_idx, 0, num_lon-1, true));
        } else if (float_eq(all_points[vpole_idx].y, 90)) {
            /* north pole */
            for (int i = 0; i < num_lon - 1; i++)
                all_leaf_triangles.push_back(allocate_triangle(vpole_idx, (num_lat-1)*num_lon+i, (num_lat-1)*num_lon+i+1));
            all_leaf_triangles.push_back(allocate_triangle(vpole_idx, (num_lat-1)*num_lon+num_lon-1, (num_lat-1)*num_lon+0, true));
        } else {
            PDASSERT(false);
        }
    }
}


void Delaunay_Voronoi::triangulate()
{
    PDASSERT(num_points > 0);

    //save_original_points_into_file();

    for(unsigned i = 1; i <= dirty_triangles_count; i++)
        if(triangle_stack[i] && triangle_stack[i]->is_leaf)
            triangulating_process(triangle_stack[i], dirty_triangles_count);

    map_buffer_index_to_point_index();

    make_final_triangle();

    mark_special_triangles();
#ifdef DEBUG
    //validate_result();
#endif
}


void Delaunay_Voronoi::make_final_triangle()
{
    all_leaf_triangles.clear();
    triangle_allocator.get_all_leaf_triangle(all_leaf_triangles);
}


void Delaunay_Voronoi::validate_result()
{
    for (unsigned i = 0; i < all_leaf_triangles.size(); i ++)
        if (all_leaf_triangles[i]->is_virtual) {
            Triangle* lt = all_leaf_triangles[i];
            all_leaf_triangles[i]->calulate_circum_circle(vertex(lt, 0), vertex(lt, 1), vertex(lt, 2));
        }

    bool valid = true;
    for (unsigned i = 0; i < all_leaf_triangles.size(); i ++) {
        if (!all_leaf_triangles[i]->is_leaf || all_leaf_triangles[i]->is_virtual)
            continue;

        if (!is_delaunay_legal(all_leaf_triangles[i])) {
            valid = false;
            break;
        }
    }
    PDASSERT(valid);

    bool* have_triangle = new bool[num_points - PAT_NUM_LOCAL_VPOINTS]();
    for (unsigned i = 0; i < all_leaf_triangles.size(); i ++) {
        if (!all_leaf_triangles[i]->is_leaf || all_leaf_triangles[i]->is_virtual)
            continue;

        have_triangle[vertex(all_leaf_triangles[i], 0)->id] = true;
        have_triangle[vertex(all_leaf_triangles[i], 1)->id] = true;
        have_triangle[vertex(all_leaf_triangles[i], 2)->id] = true;
    }

    for (int i = 0; i < num_points - PAT_NUM_LOCAL_VPOINTS; i ++)
        PDASSERT(have_triangle[i]);
}


void Delaunay_Voronoi::set_checksum_bound(double min_x, double max_x, double min_y, double max_y, double threshold)
{
    have_bound = true;
    bound_vertexes[0] = Point(min_x, min_y);
    bound_vertexes[1] = Point(max_x, min_y);
    bound_vertexes[2] = Point(max_x, max_y);
    bound_vertexes[3] = Point(min_x, max_y);
    checking_threshold = threshold;
}


static inline unsigned long long hash_two_double(double a, double b)
{
    unsigned long long *aa = reinterpret_cast<unsigned long long*>(&a);
    unsigned long long *bb = reinterpret_cast<unsigned long long*>(&b);
    return (*aa > 1) ^ *bb;
}


bool have_redundent_points(const double *x, const double *y, int num)
{
    std::tr1::unordered_map<unsigned long long, std::list<int> > hash_table;
    std::tr1::unordered_map<unsigned long long, std::list<int> >::iterator it_hash;

    if(num == 0)
        return false;

    bool have_redundent = false;
    for(int i = 0; i < num; i++) {
        it_hash = hash_table.find(hash_two_double(x[i], y[i]));
        if(it_hash != hash_table.end()) {
            bool same = false;
            for(std::list<int>::iterator it_list = it_hash->second.begin(); it_list != it_hash->second.end(); it_list ++)
                if(x[*it_list] == x[i] && y[*it_list] == y[i]) {
                    same = true;
                    break;
                }
            if(same){
                printf("redundent_point: %.20lf, %.20lf\n", x[i], y[i]);
                have_redundent = true;
            }
            else {
                it_hash->second.push_back(i);
            }
        }
        else {
            hash_table[hash_two_double(x[i], y[i])].push_back(i);
        }
    }

    if(have_redundent)
        return true;

    return false;
}


void report_redundent_points(const double *x, const double *y, const int *index, int num)
{
    std::tr1::unordered_map<unsigned long long, std::list<int> > hash_table;
    std::tr1::unordered_map<unsigned long long, std::list<int> >::iterator it_hash;

    if (num == 0)
        return;

    bool have_redundent = false;
    for (int i = 0; i < num; i++) {
        it_hash = hash_table.find(hash_two_double(x[i], y[i]));
        if (it_hash != hash_table.end()) {
            bool same = false;
            for (std::list<int>::iterator it_list = it_hash->second.begin(); it_list != it_hash->second.end(); it_list ++)
                if (x[*it_list] == x[i] && y[*it_list] == y[i]) {
                    same = true;
                    have_redundent = true;
                    printf("Point %d same as point %d: %.20lf, %.20lf\n", index[*it_list], index[i], x[i], y[i]);
                    break;
                }
            if (!same)
                it_hash->second.push_back(i);
        } else {
            hash_table[hash_two_double(x[i], y[i])].push_back(i);
        }
    }
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

Edge* Delaunay_Voronoi::allocate_edge(int head, int tail)
{
    //Edge *new_edge = new Edge();
    Edge *new_edge = edge_allocator.newElement();
    initialize_edge(new_edge, head, tail);

    return new_edge;
}


Triangle* Delaunay_Voronoi::allocate_triangle(Edge *edge1, Edge *edge2, Edge *edge3, bool force)
{
    //Triangle *new_triangle = new Triangle();
    Triangle *new_triangle = triangle_allocator.newElement();
    initialize_triangle_with_edges(new_triangle, edge1, edge2, edge3, force);

    return new_triangle;
}


Triangle* Delaunay_Voronoi::allocate_triangle(int idx1, int idx2, int idx3, bool force)
{
    Edge* e1 = allocate_edge(idx1, idx2);
    Edge* e2 = allocate_edge(idx2, idx3);
    Edge* e3 = allocate_edge(idx3, idx1);
    return allocate_triangle(e1, e2, e3, force);
}


vector<Edge*> Delaunay_Voronoi::get_all_delaunay_edge()
{
    vector<Edge*> all_edges;

    for(unsigned i = 0; i < all_leaf_triangles.size(); i ++)
        if(all_leaf_triangles[i]->is_leaf && !all_leaf_triangles[i]->is_virtual) {
            all_edges.push_back(all_leaf_triangles[i]->edge[0]);
            all_edges.push_back(all_leaf_triangles[i]->edge[1]);
            all_edges.push_back(all_leaf_triangles[i]->edge[2]);
        }

    return all_edges;
}


vector<Edge*> Delaunay_Voronoi::get_all_legal_delaunay_edge()
{
    vector<Edge*> all_edges;

    for(unsigned i = 0; i < all_leaf_triangles.size(); i ++)
        if(all_leaf_triangles[i]->is_leaf && !all_leaf_triangles[i]->is_virtual)
            if(is_triangle_legal(all_leaf_triangles[i])) {
                all_edges.push_back(all_leaf_triangles[i]->edge[0]);
                all_edges.push_back(all_leaf_triangles[i]->edge[1]);
                all_edges.push_back(all_leaf_triangles[i]->edge[2]);
            }

    return all_edges;
}


void Delaunay_Voronoi::remove_triangles_on_segment(Point head, Point tail)
{
    for(unsigned i = 0; i < all_leaf_triangles.size(); i++) {
        if(!all_leaf_triangles[i]->is_leaf)
            continue;

        if(vertex(all_leaf_triangles[i], 0)->position_to_edge(&head, &tail) * vertex(all_leaf_triangles[i], 1)->position_to_edge(&head, &tail) > 0 &&
           vertex(all_leaf_triangles[i], 1)->position_to_edge(&head, &tail) * vertex(all_leaf_triangles[i], 2)->position_to_edge(&head, &tail) > 0 )
            continue;

        /* two points of segment is in/on triangle */
        int* v_idx = all_leaf_triangles[i]->v;
        if(head.position_to_triangle(&all_points[v_idx[0]], &all_points[v_idx[1]], &all_points[v_idx[2]]) >= 0 &&
           tail.position_to_triangle(&all_points[v_idx[0]], &all_points[v_idx[1]], &all_points[v_idx[2]]) >= 0) {
            remove_leaf_triangle(all_leaf_triangles[i]);
            continue;
        }

        /* segment is intersected with at least one edge of triangle */
        for(int j = 0; j < 3; j++)
            if ((head.position_to_edge(vertex(all_leaf_triangles[i], j), vertex(all_leaf_triangles[i], (j+1)%3)) != 0 ||
                 tail.position_to_edge(vertex(all_leaf_triangles[i], j), vertex(all_leaf_triangles[i], (j+1)%3)) != 0) &&
                 head.position_to_edge(vertex(all_leaf_triangles[i], j), vertex(all_leaf_triangles[i], (j+1)%3)) *
                 tail.position_to_edge(vertex(all_leaf_triangles[i], j), vertex(all_leaf_triangles[i], (j+1)%3)) < 0 &&
                 vertex(all_leaf_triangles[i], j)->position_to_edge(&head, &tail) *
                 vertex(all_leaf_triangles[i], (j+1)%3)->position_to_edge(&head, &tail) <= 0) {
                remove_leaf_triangle(all_leaf_triangles[i]);
                break;
            }
    }
}


inline double calculate_distence_square(double x1, double y1, double x2, double y2)
{
    return (x1-x2)*(x1-x2)+(y1-y2)*(y1-y2);
}


inline bool Delaunay_Voronoi::is_triangle_in_circle(Triangle* t, Point center, double radius)
{
    return calculate_distence_square(center.x,
                                     center.y,
                                     (vertex(t, 0)->x + vertex(t, 1)->x + vertex(t, 2)->x) * 0.33333333333,
                                     (vertex(t, 0)->y + vertex(t, 1)->y + vertex(t, 2)->y) * 0.33333333333) < radius*radius;
}


void Delaunay_Voronoi::remove_triangles_in_circle(Point center, double radius)
{
    for(unsigned i = 0; i < all_leaf_triangles.size(); i++)
        if(all_leaf_triangles[i]->is_leaf) {
            if (is_triangle_in_circle(all_leaf_triangles[i], center, radius))
                remove_leaf_triangle(all_leaf_triangles[i]);
        }
}


void Delaunay_Voronoi::mark_cyclic_triangles()
{
    PDASSERT(x_ref && y_ref);

    for(unsigned i = 0; i < all_leaf_triangles.size(); i++)
        if(all_leaf_triangles[i]->is_leaf && !all_leaf_triangles[i]->is_virtual) {
            int id[3];
            for (int j = 0; j < 3; j++)
                id[j] = vertex(all_leaf_triangles[i], j)->id;

            for (int j = 0; j < 3; j++)
                if (calculate_distance(x_ref[id[j]], y_ref[id[j]], x_ref[id[(j+1)%3]], y_ref[id[(j+1)%3]]) > PAT_CYCLIC_EDGE_THRESHOLD) {
                    all_leaf_triangles[i]->is_cyclic = true;
                    break;
                }
        }
}


inline void Delaunay_Voronoi::remove_leaf_triangle(Triangle* t)
{
    t->is_virtual = true;
}

void Delaunay_Voronoi::relegalize_all_triangles()
{
    for(unsigned i = 0; i < all_leaf_triangles.size(); i++) {
        if (!all_leaf_triangles[i]->is_leaf || all_leaf_triangles[i]->is_virtual)
            continue;

        for(unsigned j = 0; j < 3; j++) {
            relegalize_triangles(all_leaf_triangles[i]->v[j], all_leaf_triangles[i]->edge[(j+1)%3]);
        }
    }
}


void Delaunay_Voronoi::remove_triangles_on_or_out_of_boundary(double min_x, double max_x, double min_y, double max_y)
{
    for(unsigned i = 0; i < all_leaf_triangles.size(); i++)
        if(all_leaf_triangles[i]->is_leaf) {
            if (vertex(all_leaf_triangles[i], 0)->is_in_region(min_x, max_x, min_y, max_y) ||
                vertex(all_leaf_triangles[i], 1)->is_in_region(min_x, max_x, min_y, max_y) ||
                vertex(all_leaf_triangles[i], 2)->is_in_region(min_x, max_x, min_y, max_y))
                continue;

            remove_leaf_triangle(all_leaf_triangles[i]);
        }
}


static inline unsigned long hash_triangle_by_id(Triangle_inline triangle)
{
    return ((unsigned long)triangle.v[0].id * (unsigned long)triangle.v[1].id) * (unsigned long)triangle.v[2].id;
}


int Delaunay_Voronoi::bound_direction(const Point* a, const Point* b)
{
    if (a->x == b->x) {
        if (a->x == bound_vertexes[0].x)
            return PDLN_LEFT;
        else if (a->x == bound_vertexes[2].x)
            return PDLN_RIGHT;
        else
            return -1;
    } else if (a->y == b->y) {
        if (a->y == bound_vertexes[0].y)
            return PDLN_DOWN;
        else if (a->y == bound_vertexes[2].y)
            return PDLN_UP;
        else
            return -1;
    } else
        return -1;
}


unsigned long Delaunay_Voronoi::cal_checksum(Point head, Point tail, double threshold)
{
    checksum_storage.clear(); //FIXME
    if (float_eq(head.x, tail.x) && float_eq(head.y, tail.y))
        return 0;

    /* searching for previously calculated checksums */
    for (unsigned i = 0; i < checksum_storage.size(); i++)
        if (head == checksum_storage[i].first.first && tail == checksum_storage[i].first.second)
            return checksum_storage[i].second;

    /* calculating checksum */
    int dir = bound_direction(&head, &tail);
    PDASSERT(dir > -1);

    unsigned long checksum = 0;

    for(unsigned i = 0; i < bound_triangles[dir].size();) {
        if (bound_triangles[dir][i].is_cyclic) {
            if(is_triangle_intersecting_with_segment(&bound_triangles[dir][i+1], head, tail, threshold) ||
               is_triangle_intersecting_with_segment(&bound_triangles[dir][i+2], head, tail, threshold)) {
                checksum ^= hash_triangle_by_id(bound_triangles[dir][i]);
            }
            i += 3;
        } else {
            if (is_triangle_intersecting_with_segment(&bound_triangles[dir][i], head, tail, threshold)) {
                checksum ^= hash_triangle_by_id(bound_triangles[dir][i]);
            }
            i++;
        }
    }

    /* highest 4 bits are reserved */
    checksum &= 0x0FFFFFFFFFFFFFFF;


    /* storing checksum */
    checksum_storage.push_back(std::make_pair(std::make_pair(head, tail), checksum));

    /* Debug: plot triangles on common edge */
    /*
    //printf("head tail [%lf, %lf], [%lf, %lf] %u\n", head.x, head.y, tail.x, tail.y, checksum);
    unsigned size_plot = 0;
    Triangle_inline* plot_triangles = new Triangle_inline[bound_triangles[dir].size()];
    for(unsigned i = 0; i < bound_triangles[dir].size();) {
        if (bound_triangles[dir][i].is_cyclic) {
            if(is_triangle_intersecting_with_segment(&bound_triangles[dir][i+1], head, tail, threshold) ||
               is_triangle_intersecting_with_segment(&bound_triangles[dir][i+2], head, tail, threshold) )
                plot_triangles[size_plot++] = bound_triangles[dir][i+1];
            i += 3;
        } else {
            if (is_triangle_intersecting_with_segment(&bound_triangles[dir][i], head, tail, threshold))
                plot_triangles[size_plot++] = bound_triangles[dir][i];
            i++;
        }
    }

    char filename[64];
    int rank;
    MPI_Comm_rank(process_thread_mgr->get_mpi_comm(), &rank);
    snprintf(filename, 64, "log/image_boundary_triangles_%03d_%016lx", rank, checksum);
    plot_triangles_into_file(filename, plot_triangles, size_plot);

    snprintf(filename, 64, "log/boundary_triangles_%03d_%016lx", rank, checksum);
    FILE* fp = fopen(filename, "w");
    for (int i = 0; i < size_plot; i++)
        fprintf(fp, "%d, %d, %d\n", plot_triangles[i].v[0].id, 
                                    plot_triangles[i].v[1].id,
                                    plot_triangles[i].v[2].id);
    fclose(fp);

    delete[] plot_triangles;
    */

    return checksum;
}


/* including triangles intersecting with boundary */
void Delaunay_Voronoi::get_triangles_in_region(double min_x, double max_x, double min_y, double max_y, 
                                               Triangle_inline *output_triangles, int *output_num_triangles, int buf_len)
{
    int current = 0;

    if(x_ref) {
        for(unsigned i = 0; i < all_leaf_triangles.size(); i++) {
            if(!all_leaf_triangles[i]->is_leaf || all_leaf_triangles[i]->is_virtual)
                continue;

            bool in = true;
            int* v_idx = all_leaf_triangles[i]->v;
            for(int j = 0; j < 3; j++)
                if(!(x_ref[all_points[v_idx[j]].id] < max_x && x_ref[all_points[v_idx[j]].id] > min_x &&
                     y_ref[all_points[v_idx[j]].id] < max_y && y_ref[all_points[v_idx[j]].id] > min_y)) {
                    in = false;
                    break;
                }

            if(in) {
                output_triangles[current++] = Triangle_inline(Point(x_ref[all_points[v_idx[0]].id], y_ref[all_points[v_idx[0]].id], global_index[all_points[v_idx[0]].id]),
                                                              Point(x_ref[all_points[v_idx[1]].id], y_ref[all_points[v_idx[1]].id], global_index[all_points[v_idx[1]].id]),
                                                              Point(x_ref[all_points[v_idx[2]].id], y_ref[all_points[v_idx[2]].id], global_index[all_points[v_idx[2]].id]));
            }
        }
    } else {
        for(unsigned i = 0; i < all_leaf_triangles.size(); i++) {
            if(!all_leaf_triangles[i]->is_leaf || all_leaf_triangles[i]->is_virtual)
                continue;

            bool in = true;
            int* v_idx = all_leaf_triangles[i]->v;
            for(int j = 0; j < 3; j++)
                if(!(all_points[v_idx[j]].x < max_x && all_points[v_idx[j]].x > min_x && all_points[v_idx[j]].y < max_y && all_points[v_idx[j]].y > min_y)) {
                    in = false;
                    break;
                }

            if(in) {
                output_triangles[current++] = Triangle_inline(Point(all_points[v_idx[0]].x, all_points[v_idx[0]].y, global_index[all_points[v_idx[0]].id]),
                                                              Point(all_points[v_idx[1]].x, all_points[v_idx[1]].y, global_index[all_points[v_idx[1]].id]),
                                                              Point(all_points[v_idx[2]].x, all_points[v_idx[2]].y, global_index[all_points[v_idx[2]].id]));
            }
        }
    }

    for (unsigned j = 0; j < 4; j++)
        for (unsigned i = 0; i < bound_triangles[j].size();) {
            output_triangles[current++] = bound_triangles[j][i];
            if (bound_triangles[j][i].is_cyclic)
                i += 3;
            else
                i += 1;
        }

    PDASSERT(current < buf_len);
    *output_num_triangles = current;
}


void Delaunay_Voronoi::set_polar_mode(bool mode)
{
    polar_mode = mode;
}


void Delaunay_Voronoi::make_bounding_triangle_pack()
{
    bound_triangles[PDLN_DOWN].clear();
    bound_triangles[PDLN_RIGHT].clear();
    bound_triangles[PDLN_UP].clear();
    bound_triangles[PDLN_LEFT].clear();

    bool u, d, l ,r;
    if (polar_mode) {
        for(unsigned i = 0; i < all_leaf_triangles.size(); i++) {
            if(!all_leaf_triangles[i]->is_leaf || all_leaf_triangles[i]->is_virtual)
                continue;

            Triangle_inline tp = pack_triangle(all_leaf_triangles[i]);
            sort_points_in_triangle(tp);
            d = is_triangle_intersecting_with_segment(&tp, bound_vertexes[0], bound_vertexes[1], checking_threshold);
            r = is_triangle_intersecting_with_segment(&tp, bound_vertexes[1], bound_vertexes[2], checking_threshold);
            u = is_triangle_intersecting_with_segment(&tp, bound_vertexes[2], bound_vertexes[3], checking_threshold);
            l = is_triangle_intersecting_with_segment(&tp, bound_vertexes[3], bound_vertexes[0], checking_threshold);
            if ((d || r || u || l ) && is_triangle_valid(all_leaf_triangles[i])) {
                if (d) add_to_bound_triangles(tp, PDLN_DOWN);
                if (r) add_to_bound_triangles(tp, PDLN_RIGHT);
                if (u) add_to_bound_triangles(tp, PDLN_UP);
                if (l) add_to_bound_triangles(tp, PDLN_LEFT);
            }
        }
    } else {
        for(unsigned i = 0; i < all_leaf_triangles.size(); i++) {
            if(!all_leaf_triangles[i]->is_leaf || all_leaf_triangles[i]->is_virtual)
                continue;

            Triangle_inline tp = pack_triangle(all_leaf_triangles[i]);
            sort_points_in_triangle(tp);
            d = is_triangle_intersecting_with_segment(&tp, bound_vertexes[0], bound_vertexes[1], checking_threshold);
            r = is_triangle_intersecting_with_segment(&tp, bound_vertexes[1], bound_vertexes[2], checking_threshold);
            u = is_triangle_intersecting_with_segment(&tp, bound_vertexes[2], bound_vertexes[3], checking_threshold);
            l = is_triangle_intersecting_with_segment(&tp, bound_vertexes[3], bound_vertexes[0], checking_threshold);
            if ((d || r || u || l ) && is_triangle_valid(all_leaf_triangles[i])) {
                if (d) bound_triangles[PDLN_DOWN].push_back(tp);
                if (r) bound_triangles[PDLN_RIGHT].push_back(tp);
                if (u) bound_triangles[PDLN_UP].push_back(tp);
                if (l) bound_triangles[PDLN_LEFT].push_back(tp);
            }
        }
    }
}


inline bool Delaunay_Voronoi::is_triangle_valid(Triangle* tri)
{
    if (avoiding_line_head[0].id && is_triangle_on_line(tri, &avoiding_line_head[0], &avoiding_line_tail[0]))
        return false;
    if (avoiding_line_head[1].id && is_triangle_on_line(tri, &avoiding_line_head[1], &avoiding_line_tail[1]))
        return false;
    if (avoiding_circle_center[0].id && is_triangle_in_circle(tri, avoiding_circle_center[0], avoiding_circle_radius[0]))
        return false;
    if (avoiding_circle_center[1].id && is_triangle_in_circle(tri, avoiding_circle_center[1], avoiding_circle_radius[1]))
        return false;
    return true;
}


inline bool Delaunay_Voronoi::is_triangle_on_line(Triangle* tri, Point* head, Point* tail)
{
    if(vertex(tri, 0)->position_to_edge(head, tail) * vertex(tri, 1)->position_to_edge(head, tail) > 0 &&
       vertex(tri, 1)->position_to_edge(head, tail) * vertex(tri, 2)->position_to_edge(head, tail) > 0 )
        return false;

    /* two points of segment is in/on triangle */
    int* v_idx = tri->v;
    if(head->position_to_triangle(&all_points[v_idx[0]], &all_points[v_idx[1]], &all_points[v_idx[2]]) >= 0 &&
       tail->position_to_triangle(&all_points[v_idx[0]], &all_points[v_idx[1]], &all_points[v_idx[2]]) >= 0) {
        return true;
    }

    /* segment is intersected with at least one edge of triangle */
    for(int j = 0; j < 3; j++)
        if ((head->position_to_edge(vertex(tri, j), vertex(tri, (j+1)%3)) != 0 ||
             tail->position_to_edge(vertex(tri, j), vertex(tri, (j+1)%3)) != 0) &&
             head->position_to_edge(vertex(tri, j), vertex(tri, (j+1)%3)) *
             tail->position_to_edge(vertex(tri, j), vertex(tri, (j+1)%3)) < 0 &&
             vertex(tri, j)->position_to_edge(head, tail) *
             vertex(tri, (j+1)%3)->position_to_edge(head, tail) <= 0) {
            return true;
        }

    return false;
}


void Delaunay_Voronoi::set_avoiding_line(unsigned id, Point head, Point tail)
{
    avoiding_line_head[id] = head;
    avoiding_line_tail[id] = tail;
    avoiding_line_head[id].id = 1; // just a marking flag
}


void Delaunay_Voronoi::set_avoiding_circle(unsigned id, Point center, double radius)
{
    avoiding_circle_center[id] = center;
    avoiding_circle_center[id].id = 1; // just a marking flag
    avoiding_circle_radius[id] = radius;
}


inline void Delaunay_Voronoi::add_to_bound_triangles(Triangle_inline& t, unsigned type)
{
    if (t.is_cyclic) {
        bound_triangles[type].push_back(t);

        double tmp_x[3];
        for(int j = 0; j < 3; j++)
            tmp_x[j] = t.v[j].x >= 180 ? t.v[j].x - 360 : t.v[j].x;
        bound_triangles[type].push_back(Triangle_inline(Point(tmp_x[0], t.v[0].y, t.v[0].id), Point(tmp_x[1], t.v[1].y, t.v[1].id),
                                                      Point(tmp_x[2], t.v[2].y, t.v[2].id), true));

        for(int j = 0; j < 3; j++)
            tmp_x[j] = t.v[j].x < 180 ? t.v[j].x + 360 : t.v[j].x;
        bound_triangles[type].push_back(Triangle_inline(Point(tmp_x[0], t.v[0].y, t.v[0].id), Point(tmp_x[1], t.v[1].y, t.v[1].id),
                                                      Point(tmp_x[2], t.v[2].y, t.v[2].id), true));
    } else
        bound_triangles[type].push_back(t);
}


void Delaunay_Voronoi::remove_triangles_only_containing_virtual_polar()
{
    double common_lat = all_points[point_idx_to_buf_idx[vpolar_local_index]].y;

    for(unsigned i = 0; i < all_leaf_triangles.size(); i++) {
        if(!all_leaf_triangles[i]->is_leaf)
            continue;
        if(vertex(all_leaf_triangles[i], 0)->y == common_lat &&
           vertex(all_leaf_triangles[i], 1)->y == common_lat &&
           vertex(all_leaf_triangles[i], 2)->y == common_lat)
            remove_leaf_triangle(all_leaf_triangles[i]);
    }
}


void Delaunay_Voronoi::remove_triangles_containing_vertexs(int idx_start, int num)
{
    int idx_end = idx_start + num;
    for(unsigned i = 0; i < all_leaf_triangles.size(); i++) {
        if(!all_leaf_triangles[i]->is_leaf)
            continue;
        if((global_index[vertex(all_leaf_triangles[i], 0)->id] < idx_end && global_index[vertex(all_leaf_triangles[i], 0)->id] >= idx_start) ||
           (global_index[vertex(all_leaf_triangles[i], 1)->id] < idx_end && global_index[vertex(all_leaf_triangles[i], 1)->id] >= idx_start) ||
           (global_index[vertex(all_leaf_triangles[i], 2)->id] < idx_end && global_index[vertex(all_leaf_triangles[i], 2)->id] >= idx_start))
            remove_leaf_triangle(all_leaf_triangles[i]);
    }
}

#ifdef OPENCV
void Delaunay_Voronoi::plot_into_file(const char *filename, double min_x, double max_x, double min_y, double max_y)
{
    unsigned num_edges;
    double *head_coord[2], *tail_coord[2];

    num_edges = 3 * all_leaf_triangles.size();
    head_coord[0] = new double[num_edges];
    head_coord[1] = new double[num_edges];
    tail_coord[0] = new double[num_edges];
    tail_coord[1] = new double[num_edges];

    num_edges = 0;
    if (x_ref) {
        for(unsigned i = 0; i < all_leaf_triangles.size(); i ++) {
            if(all_leaf_triangles[i]->is_leaf && !all_leaf_triangles[i]->is_virtual) {
                for(unsigned j = 0; j < 3; j++) {
                    head_coord[0][num_edges]   = x_ref[head(all_leaf_triangles[i]->edge[j])->id];
                    head_coord[1][num_edges]   = y_ref[head(all_leaf_triangles[i]->edge[j])->id];
                    tail_coord[0][num_edges]   = x_ref[tail(all_leaf_triangles[i]->edge[j])->id];
                    tail_coord[1][num_edges++] = y_ref[tail(all_leaf_triangles[i]->edge[j])->id];
                }
                if(all_leaf_triangles[i]->is_cyclic)
                    for(unsigned j = num_edges-1; j > num_edges-4; j--) {
                        if(head_coord[0][j] > 180) head_coord[0][j] -= 360;
                        if(tail_coord[0][j] > 180) tail_coord[0][j] -= 360;
                    }
            }
        }
    } else {
        for(unsigned i = 0; i < all_leaf_triangles.size(); i ++) {
            if(all_leaf_triangles[i]->is_leaf && !all_leaf_triangles[i]->is_virtual) {
                for(unsigned j = 0; j < 3; j++) {
                    head_coord[0][num_edges]   = head(all_leaf_triangles[i]->edge[j])->x;
                    head_coord[1][num_edges]   = head(all_leaf_triangles[i]->edge[j])->y;
                    tail_coord[0][num_edges]   = tail(all_leaf_triangles[i]->edge[j])->x;
                    tail_coord[1][num_edges++] = tail(all_leaf_triangles[i]->edge[j])->y;
                }
                if(all_leaf_triangles[i]->is_cyclic)
                    for(unsigned j = num_edges-1; j > num_edges-4; j--) {
                        if(head_coord[0][j] > 180) head_coord[0][j] -= 360;
                        if(tail_coord[0][j] > 180) tail_coord[0][j] -= 360;
                    }
            }
        }
    }

    PDASSERT(num_edges%3 == 0);
    PDASSERT(num_edges <= 3 * all_leaf_triangles.size());
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

    num_edges = 3 * all_leaf_triangles.size();
    head_coord[0] = new double[num_edges];
    head_coord[1] = new double[num_edges];
    tail_coord[0] = new double[num_edges];
    tail_coord[1] = new double[num_edges];

    num_edges = 0;
    for(unsigned i = 0; i < all_leaf_triangles.size(); i ++)
        if(all_leaf_triangles[i] && all_leaf_triangles[i]->is_leaf)
            for(unsigned j = 0; j < 3; j++) {
                head_coord[0][num_edges] = head(all_leaf_triangles[i]->edge[j])->x;
                head_coord[1][num_edges] = head(all_leaf_triangles[i]->edge[j])->y;
                tail_coord[0][num_edges] = tail(all_leaf_triangles[i]->edge[j])->x;
                tail_coord[1][num_edges] = tail(all_leaf_triangles[i]->edge[j])->y;
                num_edges++;
            }

    PDASSERT(num_edges%3 == 0);
    PDASSERT(num_edges <= 3 * all_leaf_triangles.size());
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

    num_edges = 3 * all_leaf_triangles.size();
    head_coord[0] = new double[num_edges];
    head_coord[1] = new double[num_edges];
    tail_coord[0] = new double[num_edges];
    tail_coord[1] = new double[num_edges];

    num_edges = 0;
    for(unsigned i = 0; i < all_leaf_triangles.size(); i ++)
        if(all_leaf_triangles[i]->is_leaf && !all_leaf_triangles[i]->is_virtual)
            for(unsigned j = 0; j < 3; j++) {
                head_coord[0][num_edges] = head(all_leaf_triangles[i]->edge[j])->x;
                head_coord[1][num_edges] = head(all_leaf_triangles[i]->edge[j])->y;
                tail_coord[0][num_edges] = tail(all_leaf_triangles[i]->edge[j])->x;
                tail_coord[1][num_edges] = tail(all_leaf_triangles[i]->edge[j])->y;
                num_edges++;
            }

    PDASSERT(num_edges <= 3 * all_leaf_triangles.size());
    plot_projected_edge_into_file(filename, head_coord, tail_coord, num_edges, PDLN_PLOT_COLOR_DEFAULT, PDLN_PLOT_FILEMODE_NEW);

#ifdef DEBUG
    num_edges = 0;
    for(unsigned i = 0; i < triangles_containing_vpolar.size(); i++)
        for(unsigned j = 0; j < 3; j++) {
            head_coord[0][num_edges] = head(triangles_containing_vpolar[i]->edge[j])->x;
            head_coord[1][num_edges] = head(triangles_containing_vpolar[i]->edge[j])->y;
            tail_coord[0][num_edges] = tail(triangles_containing_vpolar[i]->edge[j])->x;
            tail_coord[1][num_edges] = tail(triangles_containing_vpolar[i]->edge[j])->y;
            num_edges++;
        }
    plot_projected_edge_into_file(filename, head_coord, tail_coord, num_edges, PDLN_PLOT_COLOR_RED, PDLN_PLOT_FILEMODE_APPEND);
#endif

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

    coord[0] = new double[num_points - PAT_NUM_LOCAL_VPOINTS];
    coord[1] = new double[num_points - PAT_NUM_LOCAL_VPOINTS];

    int num = 0;
    for(int i = PAT_NUM_LOCAL_VPOINTS; i < num_points; i ++) {
        coord[0][num] = all_points[i].x;
        coord[1][num++] = all_points[i].y;
    }

    plot_points_into_file(filename, coord[0], coord[1], NULL, num_points - PAT_NUM_LOCAL_VPOINTS, min_x, max_x, min_y, max_y);

    delete coord[0];
    delete coord[1];
}


void plot_triangles_into_file(const char *prefix, Triangle_inline *t, int num, bool plot_cyclic_triangles)
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
    for(int i = 0; i < num; i ++) {
        if (t[i].v[0].x > 360 || t[i].v[1].x > 360 || t[i].v[2].x > 360) {
            for(int j = 0; j < 3; j++)
                t[i].v[j].x -= 360;
        } else if (plot_cyclic_triangles &&
                   (t[i].v[0].calculate_distance(&t[i].v[1]) >= PAT_CYCLIC_EDGE_THRESHOLD ||
                   t[i].v[1].calculate_distance(&t[i].v[2]) >= PAT_CYCLIC_EDGE_THRESHOLD ||
                   t[i].v[2].calculate_distance(&t[i].v[0]) >= PAT_CYCLIC_EDGE_THRESHOLD )) {
            for(int j = 0; j < 3; j++)
                if(t[i].v[j].x > 180) t[i].v[j].x -= 360;
        }

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


void plot_triangles_into_file(const char *prefix, std::vector<Triangle*> t, Point* whole_points)
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
            head_coord[0][num_edges] = whole_points[t[i]->v[j]].x;
            head_coord[1][num_edges] = whole_points[t[i]->v[j]].y;
            tail_coord[0][num_edges] = whole_points[t[i]->v[(j+1)%3]].x;
            tail_coord[1][num_edges] = whole_points[t[i]->v[(j+1)%3]].y;
            num_edges++;
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
    for(int i = PAT_NUM_LOCAL_VPOINTS; i < num_points; i++)
        fprintf(fp, "%.20lf, %.20lf\n", all_points[i].x, all_points[i].y);
    fclose(fp);
}


void Delaunay_Voronoi::update_all_points_coord(double *x_values, double *y_values, int num)
{
    PDASSERT(num == num_points - PAT_NUM_LOCAL_VPOINTS);
    for(int i = 0; i < num; i++) {
        all_points[point_idx_to_buf_idx[i]].x = x_values[i];
        all_points[point_idx_to_buf_idx[i]].y = y_values[i];
    }

    for(unsigned i = 0; i < all_leaf_triangles.size(); i++) {
        if(!all_leaf_triangles[i]->is_leaf)
            continue;

        Triangle* lt = all_leaf_triangles[i];
        lt->calulate_circum_circle(vertex(lt, 0), vertex(lt, 1), vertex(lt, 2));
    }
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

Triangle_inline::Triangle_inline(Point p0, Point p1, Point p2, bool cyclic)
{
    v[0] = p0;
    v[1] = p1;
    v[2] = p2;
    is_cyclic = cyclic;
}


void Triangle_inline::check_cyclic()
{
    if (v[0].calculate_distance(&v[1]) > PAT_CYCLIC_EDGE_THRESHOLD ||
        v[1].calculate_distance(&v[2]) > PAT_CYCLIC_EDGE_THRESHOLD ||
        v[2].calculate_distance(&v[0]) > PAT_CYCLIC_EDGE_THRESHOLD )
        is_cyclic = true;
}

bool operator == (Triangle_inline t1, Triangle_inline t2)
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


void save_triangles_info_file(const char *prefix, Triangle_inline *t, int num)
{
    char filename[128];
    FILE* fp;

    snprintf(filename, 128, "%s.txt", prefix);
    fp = fopen(filename, "w");
    for(int i = 0; i < num; i++)
        fprintf(fp, "[%d] (%.15lf, %.15lf), (%.15lf, %.15lf), (%.15lf, %.15lf)\n", i, t[i].v[0].x, t[i].v[0].y, t[i].v[1].x, t[i].v[1].y, t[i].v[2].x, t[i].v[2].y);
    fclose(fp);
}
