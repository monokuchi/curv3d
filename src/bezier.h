
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

        float xCoord();
        float yCoord();
        float zCoord();

        Point3D& operator=(const Point3D& other);
        Point3D& operator+=(const Point3D& other);
        Point3D& operator-=(const Point3D& other);
        Point3D& operator*=(float scalar);
        Point3D operator+(const Point3D& rhs);
        Point3D operator-(const Point3D& rhs);
        
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

        void insertControlPoint(Point3D point);
        void removeControlPoint(int i, int j);
        Point3D generatePoint(float u, float v);
    
    private:
        int _n;
        int _m;
        std::vector<Point3D> _control_points; // Should be of size (_n+1)(_m+1)

};


#endif // BEZIER_H
