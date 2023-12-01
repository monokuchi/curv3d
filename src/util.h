
#ifndef UTIL_H
#define UTIL_H

#include <iostream>



class Point2D
{
    public:
        Point2D();
        Point2D(float x, float y);
        ~Point2D() {};

        float xCoord() const;
        float yCoord() const;

        Point2D& operator=(Point2D const & other);
        Point2D& operator+=(Point2D const & other);
        Point2D& operator-=(Point2D const & other);
        Point2D& operator*=(float scalar);
        Point2D operator+(Point2D const & rhs);
        Point2D operator-(Point2D const & rhs);
        Point2D operator*(float scalar);

        float mag();
        void printPoint2D();
        
    private:
        float _x;
        float _y;
};


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
        Point3D operator*(float scalar);

        float mag();
        void printPoint3D();
        
    private:
        float _x;
        float _y;
        float _z;
};


Point3D crossProduct(Point3D& vec_1, Point3D& vec_2);
float dotProduct(Point3D& vec_1, Point3D& vec_2);


#endif // UTIL_H
