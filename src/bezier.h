
#ifndef BEZIER_H
#define BEZIER_H

#include <iostream>
#include <vector>

#include "util.h"



class Surface
{
    public:
        Surface();
        Surface(int n, int m);
        ~Surface() {};

        // Returns the order of the surface (_n, _m)
        std::pair<int, int> getOrder() const;

        // Returns the vector of control points
        std::vector<Point3D> const & getControlPoints() const;
        // Returns the number of control points
        uint getNumControlPoints();
        
        // Return reference, mutable for surface[{i,j}] = Point3D(x, y, z)
        Point3D& controlPoint(uint const i, uint const j);
        Point3D& operator[](std::initializer_list<uint> const indices);

        // Return const reference, for cases the surface should be const, but we still want it to be readable
        const Point3D& controlPoint(uint const i, uint const j) const;
        const Point3D& operator[](std::initializer_list<uint> const indices) const;

        void removeControlPoint(int i, int j);
        Point3D generatePoint(float u, float v) const;
    
    private:
        int _n;
        int _m;
        
        // Should be of size (_n+1)(_m+1)
        std::vector<Point3D> _control_points;
};


#endif // BEZIER_H
