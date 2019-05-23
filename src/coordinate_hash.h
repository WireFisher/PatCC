/***************************************************************
  *  Copyright (c) 2019, Tsinghua University.
  *  This is a source file of PatCC.
  *  This file was initially finished by Dr. Li Liu and
  *  Haoyu Yang. If you have any problem,
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#ifndef COORDINATE_HASH_H
#define COORDINATE_HASH_H

struct Entry_Node {
    double   value;
    Entry_Node* next;
    unsigned index;
};

class Coord_Hash {
    public:
        Coord_Hash(int);

        void set_hashing_params(double, double);

        void put(double);

        int get_index(double);

        int get_num_unique_values() { return unique_values; };

        void make_sorted_index();

    private:
        int hash(double val);

        void insert_to_list(int, double);

        int    capacity;
        int    unique_values;
        double range_min;
        double range_len;
        Entry_Node* table;
        double latest_val;
        int    latest_key;
        int    latest_idx;
};
#endif
