/***************************************************************
  *  Copyright (c) 2013, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Haoyu Yang,
  *  supervised by Dr. Li Liu. If you have any problem,
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/

#include "delaunay_triangulation_2D.h"
#include <iostream>
#include <cmath>
#include <cassert>
#include <algorithm>

//typedef std::vector<Vector2<T> > std::vector<Vector2<T> >;
//typedef std::vector<Triangle<T> > std::vector<Triangle<T> >;
//typedef std::vector<Edge<T> > std::vector<Edge<T> >;

template <typename T>
Vector2<T>::Vector2()


template <typename T>
Vector2<T>::Vector2(T _x, T _y) 


template <typename T>
Vector2<T>::Vector2(const Vector2 &v)


template <typename T>
void Vector2<T>::set(const Vector2 &v)


template <typename T>
T Vector2<T>::dist2(const Vector2 &v)


template <typename T>
float Vector2<T>::dist(const Vector2 &v)



	



template <class T>
bool Triangle<T>::circumCircleContains(const Vector2<T> &v)






template <class T>
const std::vector<Triangle<T> >& Delaunay<T>::triangulate(std::vector<Vector2<T> > &vertices)


