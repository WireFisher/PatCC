#ifndef PDLN_MEMORY_POOL_H
#define PDLN_MEMORY_POOL_H

#include <cstddef>
#include <vector>
#include "triangle.h"

/* freed chunk */
struct Bin {
    Bin* next;
};


class Triangle_pool {
    public:
        Triangle_pool();
        ~Triangle_pool();

        Triangle* alloc();
        void free(Triangle*);
        Triangle* newElement();
        void deleteElement(Triangle*);

        void get_all_leaf_triangle(std::vector<Triangle*>&);
    private:

        Triangle* get_page_body(char*);
        Triangle* get_page_end(char*);
        void allocNewPage();

        size_t    pagesize;
        char*     cur_page;     // current operating page
        Triangle* top_chunk;    // first unallocated chunk
        Triangle* end_chunk;
        Bin*      bins;         // freed chunk list
};


class Edge_pool {
    public:
        Edge_pool();
        ~Edge_pool();

        Edge* alloc();
        void  free(Edge*);
        Edge* newElement();
        void  deleteElement(Edge*);

    private:

        void allocNewPage();

        size_t pagesize;
        char*  cur_page;     // current operating page
        Edge*  top_chunk;    // first unallocated chunk
        Edge*  end_chunk;
        Bin*   bins;         // freed chunk list
};
#endif
