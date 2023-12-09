
#include "Scene.hpp"
#include "Utils.hpp"

#include <fstream>

using namespace curv3d;

void surface_points(const Surface& surf) {
    std::fstream file;
    file.open("scripts/points.txt", std::fstream::in | std::fstream::out | std::fstream::trunc);
    ERROR(!file.is_open(), "Can't open file");

    uint divisions = 100;
    for (uint i = 0; i < divisions; ++i) {
        for (uint j = 0; j < divisions; ++j) {
            SPoint uv(i / (real)divisions, j / (real)divisions);
            LPoint point = surf.at(uv);
            file << point.x << ", " << point.y << ", " << point.z << std::endl;
        }
    }

    file << "ControlPoints" << std::endl;

    for (uint iPoint = 0; iPoint < surf.num_points; ++iPoint) {
        GPoint point = surf.control(iPoint) + surf.origin;
        file << point.x << ", " << point.y << ", " << point.z << std::endl;
    }
    file.close();
}

int main() {
    TwoIndex resolution(640, 480);
    curv3d::Camera camera(GPoint(-5, 0, 0), resolution, 1, 1.2, Angles(0, 0, 0));
    ImageBuffer buffer(resolution);
    Scene scene(camera, &buffer);

    // Make object.
    Object obj(GPoint(-0.5, -0.5, -0.5));
    for (uint i = 0; i < 3; ++i) {
        for (uint j = 0; j < 3; ++j) {
            if (i == 1 and j == 1) {
                obj.add_point(Vector(0.5 * i, 0.5 * j, 1));
            } else {
                obj.add_point(Vector(0.5 * i, 0.5 * j, ((i + j) % 2) ? 0.5 : 0));
            }
        }
    }
    obj.add_surface(Order(2, 2), {{0, 1, 2}, {3, 4, 5}, {6, 7, 8}});
    scene.add_object(obj);

    surface_points(obj.surface(0));

    std::vector<uint> iters = scene.slow_rasterize();
    // Can use for analyzing convergences.

    return 0;
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
// }
