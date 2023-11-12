
#include <bezier.h>
#include <cmath>
#include <fstream>



float factorial(int n)
{
    return std::tgamma(n+1);
}

float bernstein(int n, int i, float u)
{
    return (factorial(n) / (factorial(i)*factorial(n-i))) * std::pow(u, i) * std::pow((1-u), n-i);
}


void printPoint3D(Point3D point)
{
    std::cout << "x: " << point.xCoord() << "  y: " << point.yCoord() << "  z: " << point.zCoord() << std::endl;
}

Point3D::Point3D()
{
    _x = 0.0;
    _y = 0.0;
    _z = 0.0;
}

Point3D::Point3D(float x, float y, float z)
{
    _x = x;
    _y = y;
    _z = z;
}

float Point3D::xCoord() const
{
    return _x;
}

float Point3D::yCoord() const
{
    return _y;
}

float Point3D::zCoord() const
{
    return _z;
}

Point3D& Point3D::operator=(Point3D const & other)
{
    if (this == &other)
    {
        return *this;
    }
    
    this->_x = other._x;
    this->_y = other._y;
    this->_z = other._z;

    return *this;
}

Point3D& Point3D::operator+=(Point3D const & other)
{
    this->_x += other._x;
    this->_y += other._y;
    this->_z += other._z;

    return *this;
}

Point3D& Point3D::operator-=(Point3D const & other)
{
    this->_x += other._x;
    this->_y += other._y;
    this->_z += other._z;

    return *this;
}

Point3D& Point3D::operator*=(float scalar)
{
    this->_x *= scalar;
    this->_y *= scalar;
    this->_z *= scalar;

    return *this;
}

Point3D Point3D::operator+(Point3D const & rhs)
{
    Point3D point = *this;

    point._x += rhs._x;
    point._y += rhs._y;
    point._z += rhs._z;

    return point;
}

Point3D Point3D::operator-(Point3D const & rhs)
{
    Point3D point = *this;

    point._x -= rhs._x;
    point._y -= rhs._y;
    point._z -= rhs._z;

    return point;
}


Surface::Surface()
{
    _n = 0;
    _m = 0;

    _control_points.reserve((_n+1)*(_m+1));
}

Surface::Surface(int n, int m)
{
    _n = n;
    _m = m;

    _control_points.reserve((_n+1)*(_m+1));
}

std::vector<Point3D> const & Surface::getControlPoints() const
{
    return _control_points;
}

void Surface::insertControlPoint(Point3D point)
{
    if (_control_points.size() < (_n+1)*(_m+1))
    {
        _control_points.push_back(point);
    }
    else
    {
        std::cout << "Exceeded total number of allowed control points!" << std::endl;
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

Point3D Surface::generatePoint(float u, float v)
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


int main()
{
    /* PLANE */
    // Surface surface(1, 1);
    // surface.insertControlPoint(Point3D(0, 0, 0));
    // surface.insertControlPoint(Point3D(1, 0, 0));
    // surface.insertControlPoint(Point3D(0, 1, 0));
    // surface.insertControlPoint(Point3D(1, 1, 0));


    /* CUBE */
    // Surface surface(2, 2);
    // surface.insertControlPoint(Point3D(0, 0, 0));
    // surface.insertControlPoint(Point3D(1, 0, 0));
    // surface.insertControlPoint(Point3D(1, 1, 0));
    // surface.insertControlPoint(Point3D(0, 1, 0));
    // surface.insertControlPoint(Point3D(0, 0, 1));
    // surface.insertControlPoint(Point3D(1, 0, 1));
    // surface.insertControlPoint(Point3D(1, 1, 1));
    // surface.insertControlPoint(Point3D(0, 1, 1));


    /* DIP */
    // Surface surface(1, 1);
    // surface.insertControlPoint(Point3D(0, 0, .5));
    // surface.insertControlPoint(Point3D(1, .3, .5));
    // surface.insertControlPoint(Point3D(.3, 1, .5));
    // surface.insertControlPoint(Point3D(.7, .6, .2));


    /* Umbrella */
    Surface surface(2, 2);
    surface.insertControlPoint(Point3D(0, 0, 0));
    surface.insertControlPoint(Point3D(.5, 0, .5));
    surface.insertControlPoint(Point3D(1, 0, 0));
    surface.insertControlPoint(Point3D(0, .5, .5));
    surface.insertControlPoint(Point3D(0, 1, 0));
    surface.insertControlPoint(Point3D(.5, 1, .5));
    surface.insertControlPoint(Point3D(1, 1, 0));
    surface.insertControlPoint(Point3D(1, .5, .5));
    surface.insertControlPoint(Point3D(.5, .5, 1));



    
    // Generate the points
    std::vector<Point3D> points;
    for (int i=0; i<=100; i++)
    {
        for (int j=0; j<=100; j++)
        {
            float u = i / 100.0;
            float v = j / 100.0;

            Point3D point = surface.generatePoint(u, v);
            printPoint3D(point);
            points.push_back(point);
        }
    }

    // Get the control points
    std::vector<Point3D> control_points = surface.getControlPoints();

    // Store points into a text file
    std::fstream file;
    file.open("../scripts/points.txt", std::fstream::in | std::fstream::out | std::fstream::trunc);
    if (file.is_open())
    {
        for (auto& point : points)
        {
            file << point.xCoord() << ", " << point.yCoord() << ", " << point.zCoord() << std::endl;
        }

        file << "ControlPoints" << std::endl;

        for (auto& control_point : control_points)
        {
            file << control_point.xCoord() << ", " << control_point.yCoord() << ", " << control_point.zCoord() << std::endl;
        }
    }
    else
    {
        std::cout << "Can't open file!" << std::endl;
    }
    file.close();
}
