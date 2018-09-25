#include "memory_pool.h"
#include <cstring>

#define PDLN_MEMORY_POOL_SIZE (0x10000)
Triangle_pool::Triangle_pool()
    : pagesize(PDLN_MEMORY_POOL_SIZE)
    , cur_page(NULL)
    , top_chunk(NULL)
    , end_chunk(NULL)
    , bins(NULL)
{
}


Triangle_pool::~Triangle_pool()
{
    char* cur = cur_page;
    while(cur) {
        char* next = (char*)((Bin*)cur)->next;
        operator delete(cur);
        cur = next;
    }
}


size_t padding(char* p, size_t align)
{
    size_t result = (size_t)p;
    return ((align - result) % align);
}


Triangle* Triangle_pool::get_page_body(char* page)
{
    char* body = page + sizeof(Bin);
    body += padding(body, sizeof(Triangle));
    return (Triangle*)body;
}


Triangle* Triangle_pool::get_page_end(char* page)
{
    return (Triangle*)(page + pagesize - sizeof(Triangle) + 1);
}


void Triangle_pool::allocNewPage()
{
    char* new_page = (char*)operator new(pagesize);
    memset(new_page, 0, pagesize);
    ((Bin*)new_page)->next = (Bin*)cur_page;
    cur_page = new_page;

    char* page_body = new_page + sizeof(Bin);
    page_body += padding(page_body, sizeof(Triangle));
    top_chunk = (Triangle*)page_body;
    end_chunk = (Triangle*)(new_page + pagesize - sizeof(Triangle) + 1);
}


Triangle* Triangle_pool::alloc()
{
    if (bins) {
        Triangle* c = (Triangle*)bins;
        bins = bins->next;
        return c;
    } else {
        if (top_chunk >= end_chunk)
            allocNewPage();
        return top_chunk++;
    }
}


void Triangle_pool::free(Triangle* chunk)
{
    Bin* bin = (Bin*)chunk;
    bin->next = bins;
    bins = bin;
}


Triangle* Triangle_pool::newElement()
{
    Triangle* result = alloc();
    new (result) Triangle();
    return result;
}

void Triangle_pool::deleteElement(Triangle* c)
{
    if(c) {
        c->~Triangle();
        free(c);
    }
}


void Triangle_pool::get_all_leaf_triangle(std::vector<Triangle*>& all/*, std::vector<Triangle*>& boundary*/)
{
    char* curpage = cur_page;
    while (curpage) {
        Triangle* begin = get_page_body(curpage);
        Triangle* end   = get_page_end(curpage);
        for (;begin < end; begin++)
            if (begin->is_leaf) {
                all.push_back(begin);

                //Triangle_pack tp = pack_triangle(triangle);
                //if (is_triangle_intersecting_with_segment(&tp, bound_vertexes[0], bound_vertexes[1], checking_threshold) ||
                //    is_triangle_intersecting_with_segment(&tp, bound_vertexes[1], bound_vertexes[2], checking_threshold) ||
                //    is_triangle_intersecting_with_segment(&tp, bound_vertexes[2], bound_vertexes[3], checking_threshold) ||
                //    is_triangle_intersecting_with_segment(&tp, bound_vertexes[3], bound_vertexes[0], checking_threshold) )
                //    all_leaf_triangles_on_boundary.push_back(triangle);

            }
        curpage = (char*)((Bin*)curpage)->next;
    }
}


Edge_pool::Edge_pool()
    : pagesize(PDLN_MEMORY_POOL_SIZE)
    , cur_page(NULL)
    , top_chunk(NULL)
    , end_chunk(NULL)
    , bins(NULL)
{
}


Edge_pool::~Edge_pool()
{
    char* cur = cur_page;
    while(cur) {
        char* next = (char*)((Bin*)cur)->next;
        operator delete(cur);
        cur = next;
    }
}


void Edge_pool::allocNewPage()
{
    char* new_page = (char*)operator new(pagesize);
    memset(new_page, 0, pagesize);
    ((Bin*)new_page)->next = (Bin*)cur_page;
    cur_page = new_page;

    char* page_body = new_page + sizeof(Bin);
    page_body += padding(page_body, sizeof(Edge));
    top_chunk = (Edge*)page_body;
    end_chunk = (Edge*)(new_page + pagesize - sizeof(Edge) + 1);
}


Edge* Edge_pool::alloc()
{
    if (bins) {
        Edge* c = (Edge*)bins;
        bins = bins->next;
        return c;
    } else {
        if (top_chunk >= end_chunk)
            allocNewPage();
        return top_chunk++;
    }
}


void Edge_pool::free(Edge* chunk)
{
    Bin* bin = (Bin*)chunk;
    bin->next = bins;
    bins = bin;
}

Edge* Edge_pool::newElement(Point* p1, Point* p2)
{
    Edge* result = alloc();
    new (result) Edge(p1, p2);
    return result;
}

void Edge_pool::deleteElement(Edge* c)
{
    if(c) {
        c->~Edge();
        free(c);
    }
}
