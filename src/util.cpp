
#include <util.h>
#include <cmath>



Point2D::Point2D()
{
    _x = 0.0f;
    _y = 0.0f;
}

Point2D::Point2D(float x, float y)
{
    _x = x;
    _y = y;
}

float Point2D::xCoord() const
{
    return _x;
}

float Point2D::yCoord() const
{
    return _y;
}

Point2D& Point2D::operator=(const Point2D& other)
{
    if (this == &other)
    {
        return *this;
    }
    
    this->_x = other._x;
    this->_y = other._y;

    return *this;
}

Point2D& Point2D::operator+=(const Point2D& other)
{
    this->_x += other._x;
    this->_y += other._y;

    return *this;
}

Point2D& Point2D::operator-=(const Point2D& other)
{
    this->_x -= other._x;
    this->_y -= other._y;

    return *this;
}

Point2D& Point2D::operator*=(float scalar)
{
    this->_x *= scalar;
    this->_y *= scalar;

    return *this;
}

Point2D Point2D::operator+(const Point2D& rhs)
{
    Point2D point = *this;

    point._x += rhs._x;
    point._y += rhs._y;

    return point;
}

Point2D Point2D::operator-(const Point2D& rhs)
{
    Point2D point = *this;

    point._x -= rhs._x;
    point._y -= rhs._y;

    return point;
}

Point2D Point2D::operator*(float scalar)
{
    Point2D point = *this;

    point._x *= scalar;
    point._y *= scalar;

    return point;
}

float Point2D::mag()
{
    return std::sqrt(std::pow(_x, 2.0f) + std::pow(_y, 2.0f));
}

void Point2D::printPoint2D()
{
    std::cout << "x: " << this->xCoord() << "  y: " << this->yCoord() << std::endl;
}


Point3D::Point3D()
{
    _x = 0.0f;
    _y = 0.0f;
    _z = 0.0f;
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
    this->_x -= other._x;
    this->_y -= other._y;
    this->_z -= other._z;

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

Point3D Point3D::operator*(float scalar)
{
    Point3D point = *this;

    point._x *= scalar;
    point._y *= scalar;
    point._z *= scalar;

    return point;
}

float Point3D::mag()
{
    return std::sqrt(std::pow(_x, 2.0f) + std::pow(_y, 2.0f) + std::pow(_z, 2.0f));
}

void Point3D::printPoint3D()
{
    std::cout << "x: " << this->xCoord() << "  y: " << this->yCoord() << "  z: " << this->zCoord() << std::endl;
}


Point3D crossProduct(Point3D& vec_1, Point3D& vec_2)
{
    return Point3D(vec_1.yCoord()*vec_2.zCoord() - vec_1.zCoord()*vec_2.yCoord(), vec_1.zCoord()*vec_2.xCoord() - vec_1.xCoord()*vec_2.zCoord(), vec_1.xCoord()*vec_2.yCoord() - vec_1.yCoord()*vec_2.xCoord());
}

float dotProduct(Point3D& vec_1, Point3D& vec_2)
{
    return vec_1.xCoord()*vec_2.xCoord() + vec_1.yCoord()*vec_2.yCoord() + vec_1.zCoord()*vec_2.zCoord();
}

float dotProduct(Point2D& vec_1, Point2D& vec_2)
{
    return vec_1.xCoord()*vec_2.xCoord() + vec_1.yCoord()*vec_2.yCoord();
}
