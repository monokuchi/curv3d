
#ifndef RASTERIZE_H
#define RASTERIZE_H

#include <iostream>
#include <vector>

#include <bezier.h>
#include <eigen3/Eigen/Dense>



typedef struct RGB
{
    RGB(float r, float g, float b);

    // Color (RGB)
    float r;
    float g;
    float b;
} RGB;


typedef struct Pixel
{
    Pixel(int x, int y, RGB rgb);

    // Position (x, y)
    Point2D position;

    // Color (RGB)
    RGB color;
} Pixel;


typedef struct BoundingBox
{
    // These two points define the rectangular bounding box
    Point2D upper_left_corner;
    Point2D lower_right_corner;

    // Returns true if point is in bounding box
    bool inBoundingBox(Point2D const & point) const;
} BoundingBox;


typedef struct ZBuffer
{
    ZBuffer(int width, int height);

    // Should be of size (width * height)
    std::vector<float> z_buffer;
} ZBuffer;


typedef struct Image
{
    std::vector<Pixel> pixels;
} Image;


class Camera
{
    public:
        Camera();
        Camera(Point3D focal_point, Point3D orientation_vector, float orientation_angle, float size, uint res_x, uint res_y);
        ~Camera() {};

        
        // Sets focal point
        void position(float new_x, float new_y, float new_z);
        // Moves focal point based off a change in x, y, z
        void move(float dx, float dy, float dz);

        // Sets the rotation of the camera (both orientation vector and angle), uses mathematics convention for spherical coordinates
        void rotation(float theta, float phi, float orientation_angle);
        // Sets rotation of the camera based off a change in angles, uses mathematics convention for spherical coordinates
        void rotate(float d_theta, float d_phi, float d_orientation_angle);

        // Changes focal length (orientation vector)
        void focalLength(float length);
        // Sets new orientation vector going in (x, y, z) direction (can be any magnitude)
        void direction(float x, float y, float z);

        void calculateIntersections(Point3D pixel_to_camera, Surface& surface);
        Point2D calculateGradient(Surface& surface, float u, float v);
        uint gradientDescent(Point2D point, Surface& surface, float gamma, float u, float v);
        Image projectSurface(Surface& surface);

    private:
        // Focal point of the camera in the world space
        Point3D _focal_point;
        // Orientation vector with respect to the camera space
        Point3D _orientation_vector;
        // Orientation angle (0 degrees is upright, 180 degrees is upside down)
        float _orientation_angle;

        // Physical size of the camera (width)
        float _size;
        // Number of horizontal pixels on the screen
        uint _resolution_x;
        // Number of vertical pixels on the screen
        uint _resolution_y;

        // Returns a set of normalized basis vectors of the screen space
        std::pair<Point3D, Point3D> getScreenBasis();
};


// class Screen
// {
//     public:
//         Screen();
//         Screen(int width, int height);
//         ~Screen() {};

//     private:
//         std::vector<Pixel> _pixels;
//         int _width;
//         int _height;
// };


class World
{
    public:
        World();
        ~World() {};

        std::vector<Surface> getSurfaces();
        void insertSurface(Surface surface);

    private:
        std::vector<Surface> _surfaces;
};


class Transform
{
    public:
        Transform();
        ~Transform() {};

    private:
        Eigen::MatrixXf _objectTransform;
        Eigen::MatrixXf _worldTransform;
        Eigen::MatrixXf _cameraTransform;
        Eigen::MatrixXf _screenTransform;
};



#endif // RASTERIZE_H