
#ifndef BEZIER_H
#define BEZIER_H

#include <iostream>
#include <vector>



// typedef struct Point3D
// {
//     float x;
//     float y;
//     float z;
// } Point3D;

class Point3D
{
    public:
        Point3D();
        Point3D(float x, float y, float z);
        ~Point3D() {};

        float xCoord() const;
        float yCoord() const;
        float zCoord() const;

        Point3D& operator=(Point3D const & other);
        Point3D& operator+=(Point3D const & other);
        Point3D& operator-=(Point3D const & other);
        Point3D& operator*=(float scalar);
        Point3D operator+(Point3D const & rhs);
        Point3D operator-(Point3D const & rhs);
        
    private:
        float _x;
        float _y;
        float _z;
};


class Surface
{
    public:
        Surface();
        Surface(int n, int m);
        ~Surface() {};

        std::vector<Point3D> const & getControlPoints() const;
        
        void insertControlPoint(Point3D point);
        void removeControlPoint(int i, int j);
        Point3D generatePoint(float u, float v);
    
    private:
        int _n;
        int _m;
        std::vector<Point3D> _control_points; // Should be of size (_n+1)(_m+1)

};


#endif // BEZIER_H
