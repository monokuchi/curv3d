
#include <bezier.h>
#include <cmath>



float factorial(int n)
{
    return std::tgamma(n+1);
}

float bernstein(int n, int i, float u)
{
    return (factorial(n) / (factorial(i)*factorial(n-i))) * std::pow(u, i) * std::pow((1-u), n-i);
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

float Point3D::xCoord()
{
    return _x;
}

float Point3D::yCoord()
{
    return _y;
}

float Point3D::zCoord()
{
    return _z;
}

Point3D& Point3D::operator=(const Point3D& other)
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

Point3D& Point3D::operator+=(const Point3D& other)
{
    this->_x += other._x;
    this->_y += other._y;
    this->_z += other._z;

    return *this;
}

Point3D& Point3D::operator-=(const Point3D& other)
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

Point3D Point3D::operator+(const Point3D& rhs)
{
    Point3D point = *this;

    point._x += rhs._x;
    point._y += rhs._y;
    point._z += rhs._z;

    return point;
}

Point3D Point3D::operator-(const Point3D& rhs)
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
    size_t pos = (j*_n + i) + 1;
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
    Point3D point = {0.0, 0.0, 0.0};

    for (size_t i=0; i<u; i++)
    {
        for (size_t j=0; j<v; j++)
        {
            // Bernstein Polynomials
            float b_i = bernstein(_n, i, u);
            float b_j = bernstein(_m, j, v);

            std::cout << "i: " << b_i << std::endl;
            std::cout << "j: " << b_j << std::endl;

            // Calculate the scaled control point
            Point3D scaled_control_point = _control_points[j*_n + i];
            scaled_control_point *= (b_i * b_j);

            // Summate
            point += scaled_control_point;
        }
    }

    return point;
}

void printPoint3D(Point3D point)
{
    std::cout << "x: " << point.xCoord() << std::endl;
    std::cout << "y: " << point.yCoord() << std::endl;
    std::cout << "z: " << point.zCoord() << std::endl;
}


int main()
{
    Surface surface(2, 2);
    surface.insertControlPoint(Point3D(0, 0, 0));
    surface.insertControlPoint(Point3D(1, 0, 0));
    surface.insertControlPoint(Point3D(0, 1, 0));
    surface.insertControlPoint(Point3D(1, 1, 0));
    
    printPoint3D(surface.generatePoint(1, 2));
}
