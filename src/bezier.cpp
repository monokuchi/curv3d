
#include <cmath>
#include <fstream>

#include "bezier.h"



float factorial(int n)
{
    return std::tgamma(n+1);
}

float bernstein(int n, int i, float u)
{
    return (factorial(n) / (factorial(i)*factorial(n-i))) * std::pow(u, i) * std::pow((1-u), n-i);
}


Surface::Surface()
{
    _n = 0;
    _m = 0;

    _control_points.resize((_n+1)*(_m+1));
}

Surface::Surface(int n, int m)
{
    _n = n;
    _m = m;

    _control_points.resize((_n+1)*(_m+1));
}

std::pair<int, int> Surface::getOrder() const
{
    return std::make_pair(_n, _m);
}

const std::vector<Point3D>& Surface::getControlPoints() const
{
    return _control_points;
}

uint Surface::getNumControlPoints()
{
    return (_n+1)*(_m+1);
}

Point3D& Surface::controlPoint(uint const i, uint const j)
{
    // surface.controlPoint(i, j) = Point3D(x, y, z);
    if (i*(_n+1) + j > (_n+1)*(_m+1))
    {
        std::cout << "Exceeded total number of allowed control points!" << std::endl;
        throw std::exception();
    }
    else
    {
        return _control_points[i*(_n+1) + j];
    }
}

Point3D& Surface::operator[](std::initializer_list<uint> const indices)
{
    // surface[{i, j}] = Point3D(x, y, z);
    if (indices.size() == 2)
    {
        return controlPoint(*indices.begin(), *(indices.end()-1));
    }
    else
    {
        std::cout << "Number of indices is not 2!" << std::endl;
        throw std::exception();
    }
}

const Point3D& Surface::controlPoint(uint const i, uint const j) const
{
    // surface.controlPoint(i, j) = Point3D(x, y, z);
    if (i*(_n+1) + j > (_n+1)*(_m+1))
    {
        std::cout << "Exceeded total number of allowed control points!" << std::endl;
        throw std::exception();
    }
    else
    {
        return _control_points[i*(_n+1) + j];
    }
}

const Point3D& Surface::operator[](std::initializer_list<uint> const indices) const
{
    // surface[{i, j}] = Point3D(x, y, z);
    if (indices.size() == 2)
    {
        return controlPoint(*indices.begin(), *(indices.end()-1));
    }
    else
    {
        std::cout << "Number of indices is not 2!" << std::endl;
        throw std::exception();
    }
}

void Surface::removeControlPoint(int i, int j)
{
    size_t pos = (i*(_n+1) + j) + 1;
    if (pos <= _control_points.size())
    {
        _control_points.erase(_control_points.begin() + pos);
    }
    else
    {
        std::cout << "Invalid index for _control_points!" << std::endl;
    }
}

Point3D Surface::generatePoint(float u, float v) const
{
    Point3D point;

    for (size_t i=0; i<=_n; i++)
    {
        for (size_t j=0; j<=_m; j++)
        {
            // Bernstein Basis Polynomials
            float b_i = bernstein(_n, i, u);
            float b_j = bernstein(_m, j, v);

            // std::cout << "b_i: " << b_i << std::endl;
            // std::cout << "b_j: " << b_j << std::endl;

            // Calculate the scaled control point
            Point3D scaled_control_point = _control_points[i*(_n+1) + j];
            scaled_control_point *= (b_i * b_j);

            // Summate
            point += scaled_control_point;
        }
    }

    return point;
}


// int main()
// {
//     /* PLANE */
//     // Surface surface(1, 1);
//     // surface[{0, 0}] = Point3D(0, 0, 0);
//     // surface[{1, 0}] = Point3D(1, 0, 0);
//     // surface[{0, 1}] = Point3D(0, 1, 0);
//     // surface[{1, 1}] = Point3D(1, 1, 0);

//     // surface.controlPoint(0, 0) = Point3D(0, 0, 0);
//     // surface.controlPoint(1, 0) = Point3D(1, 0, 0);
//     // surface.controlPoint(0, 1) = Point3D(0, 1, 0);
//     // surface.controlPoint(1, 1) = Point3D(1, 1, 0);


//     /* DIP */
//     Surface surface(1, 1);
//     surface[{0, 0}] = Point3D(0, 0, .5);
//     surface[{1, 0}] = Point3D(1, .3, .5);
//     surface[{0, 1}] = Point3D(.3, 1, .5);
//     surface[{1, 1}] = Point3D(.7, .6, .2);


//     /* UMBRELLA */
//     // Surface surface(2, 2);
//     // surface[{0, 0}] = Point3D(0, 0, 0);
//     // surface[{1, 0}] = Point3D(.5, 0, .5);
//     // surface[{2, 0}] = Point3D(1, 0, 0);
//     // surface[{0, 1}] = Point3D(0, .5, .5);
//     // surface[{1, 1}] = Point3D(.5, .5, 1);
//     // surface[{2, 1}] = Point3D(1, .5, .5);
//     // surface[{0, 2}] = Point3D(0, 1, 0);
//     // surface[{1, 2}] = Point3D(.5, 1, .5);
//     // surface[{2, 2}] = Point3D(1, 1, 0);




    
//     // Generate the points
//     std::vector<Point3D> points;
//     for (size_t i=0; i<=100; i++)
//     {
//         for (size_t j=0; j<=100; j++)
//         {
//             float u = i / 100.0;
//             float v = j / 100.0;

//             Point3D point = surface.generatePoint(u, v);
//             point.printPoint3D();
//             points.push_back(point);
//         }
//     }

//     // Get the control points
//     std::vector<Point3D> control_points = surface.getControlPoints();

//     // Store points into a text file
//     std::fstream file;
//     file.open("../scripts/points.txt", std::fstream::in | std::fstream::out | std::fstream::trunc);
//     if (file.is_open())
//     {
//         for (auto& point : points)
//         {
//             file << point.xCoord() << ", " << point.yCoord() << ", " << point.zCoord() << std::endl;
//         }

//         file << "ControlPoints" << std::endl;

//         for (auto& control_point : control_points)
//         {
//             file << control_point.xCoord() << ", " << control_point.yCoord() << ", " << control_point.zCoord() << std::endl;
//         }
//     }
//     else
//     {
//         std::cout << "Can't open file!" << std::endl;
//     }
//     file.close();
// }
