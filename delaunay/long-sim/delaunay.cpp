#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>

#include <cstdio>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef CGAL::Delaunay_triangulation_2<K>  Triangulation;
typedef Triangulation::Edge_iterator  Edge_iterator;
typedef Triangulation::Point          Point;

typedef CGAL::Triangulation_data_structure_2
    <CGAL::Triangulation_vertex_base_2<K> >::Face_handle  Face_handle;

int main() {
    Triangulation T;
    Edge_iterator iter;

    T.insert(Point(2.0, 1.0));
    T.insert(Point(0.0, 0.0));
    T.insert(Point(-1.0, -2.0));
    T.insert(Point(-2.0, -1.0));
    T.insert(Point(1.0, 2.0));

    for(Edge_iterator iter = T.edges_begin(); iter != T.edges_end(); ++iter) {
        Face_handle curr = iter->first;
        int n = iter->second;

        std::printf(
            "P1: (%.3f, %.3f) P2: (%.3f, %.3f)\n",
            curr->vertex(curr->cw(n))->point().x(),
            curr->vertex(curr->cw(n))->point().y(),
            curr->vertex(curr->ccw(n))->point().x(),
            curr->vertex(curr->ccw(n))->point().y());
    } // end for(Edge_iterator iter = T.edges_Begin(); iter != T.edges_end(); ++iter) 

    return 0;
} // end int main() 
