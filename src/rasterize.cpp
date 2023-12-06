
#include "rasterize.h"



// Screen rasterize(World world)
// {
//     Screen screen;

//     for (auto surface : world.getSurfaces())
//     {

//     }

//     return screen;
// }






Camera::Camera()
{
    Point3D _focal_point(0.0f, 0.0f, 0.0f);
    Point3D _orientation_vector(0.0f, 0.0f, 0.0f);
    float _orientation_angle = 0.0f;

    float _size = 0.0f;
    uint _resolution_x = 0;
    uint _resolution_y = 0;
}

Camera::Camera(Point3D focal_point, Point3D orientation_vector, float orientation_angle, float size, uint res_x, uint res_y)
{
    Point3D _focal_point = focal_point;
    Point3D _orientation_vector = orientation_vector;
    float _orientation_angle = orientation_angle;

    float _size = size;
    uint _resolution_x = res_x;
    uint _resolution_y = res_y;
}

void Camera::position(float new_x, float new_y, float new_z)
{
    _focal_point = Point3D(new_x, new_y, new_z);
}

void Camera::move(float dx, float dy, float dz)
{
    position(_focal_point.xCoord()+dx, _focal_point.yCoord()+dy, _focal_point.zCoord()+dz);
}

void Camera::rotation(float theta, float phi, float orientation_angle)
{
    _orientation_vector = Point3D(std::cos(theta)*std::sin(phi), std::sin(theta)*std::sin(phi), std::cos(phi)) * _orientation_vector.mag();
    _orientation_angle = orientation_angle;
}

void Camera::rotate(float d_theta, float d_phi, float d_orientation_angle)
{
    float x = _orientation_vector.xCoord();
    float y = _orientation_vector.yCoord();
    float z = _orientation_vector.zCoord();

    float theta = std::atan2(y, x);
    float phi = acos(z / _orientation_vector.mag());

    rotation(theta+d_theta, phi+d_phi, _orientation_angle+d_orientation_angle);
}

void Camera::focalLength(float length)
{
    _orientation_vector *= (length / _orientation_vector.mag());
}

void Camera::direction(float x, float y, float z)
{
    _orientation_vector = Point3D(x, y, z);
}

void Camera::resolution(uint res_x, uint res_y)
{
    _resolution_x = res_x;
    _resolution_y = res_y;
}

void Camera::size(float width)
{
    _size = width;
}

std::pair<Point3D, Point3D> Camera::getScreenBasis()
{
    // Basis vector of screen space with respect to camera space
    Point3D u(_orientation_vector.yCoord(), -_orientation_vector.xCoord(), 0);
    // If we want to normalize u to width / _resolution_x we need to multiply u by width / (mag(u) * _resolution_x)
    u *= (_size / (u.mag() * _resolution_x));

    Point3D f = _orientation_vector;
    f *= 1 / _orientation_vector.mag();

    // v = u x f
    Point3D v = crossProduct(u, f);

    Point3D rotated_u = (u*std::cos(_orientation_angle)) + (v*std::sin(_orientation_angle));
    // rotated_v = rotated_u x f
    Point3D rotated_v = crossProduct(rotated_u, f);

    return std::make_pair(rotated_u, rotated_v);
}

void Camera::calculateIntersections(Surface& surface)
{
    std::pair<Point3D, Point3D> basis_vectors = getScreenBasis();
    // std::cout << "Screen Basis Vectors: " << std::endl;
    // basis_vectors.first.printPoint3D();
    // basis_vectors.second.printPoint3D();


    // Implement bounding box later
    std::pair<double, double> target_u_v(.5, .2);
    for (size_t i=0; i<_resolution_x; i++)
    {
        for (size_t j=0; j<_resolution_y; j++)
        {
            // screen_point is the vector from the center of the screen to Pixel(i, j)
            // Point3D screen_point = basis_vectors.first*(i-_resolution_x/2) + basis_vectors.second*(j-_resolution_y/2);
            Point3D screen_point = basis_vectors.first*(i) + basis_vectors.second*(j);
            // pixel is the vector from focal_point to Pixel(i, j)
            Point3D pixel = _orientation_vector + screen_point;
            // pixel.printPoint3D();
            

            // min |p x (f(u,v) - foc)|^2, subject to u,v in [0,1]^2 
            // Assume order of the bezier surface, f(u, v), is fixed at 2

            double gamma = .001;
            
            Point3D first = pixel - _focal_point;
            Point3D second(pixel.xCoord()+.5, pixel.yCoord(), pixel.zCoord());
            Point3D third = crossProduct(first, second);
            // epsilon = |(pixel(i,j) - foc) x (pixel(i+.5, j))|^2
            double epsilon = dotProduct(third, third);
            
            uint iteration_count = gradientDescent(surface, pixel, target_u_v, gamma, epsilon);

            std::cout << "Gradient Descent completed in " << iteration_count << " steps!" << std::endl;
            std::cout << "u: " << target_u_v.first << " v: " << target_u_v.second << std::endl;
        }
    }
}

Point2D Camera::calculateGradient(Surface& surface, Point3D& pixel, std::pair<double, double> u_v_pair)
{
    /* 
        Calculates gradient of |p x (f(u, v) - foc)|^2
        p: Pixel(i,j) - foc
        f(u, v): 3-D point on bezier surface
        foc: focal point of camera
    */

    Eigen::RowVectorXf r(surface.getNumControlPoints());
    Eigen::MatrixXf q(surface.getNumControlPoints(), 2);

    double u = u_v_pair.first;
    double v = u_v_pair.second;
    // Since we fixed order to be (n, m) = (2, 2) we can manually fill in q
    q << -2*(1-u)*std::pow(1-v, 2.0f), -2*(1-v)*std::pow(1-u, 2.0f),
    -4*(v-std::pow(v, 2.0f))*(1-u), 2*std::pow(1-u, 2.0f)*(1-(2*v)),
    -2*std::pow(v, 2.0f)*(1-u), 2*std::pow(1-u, 2.0f)*v,

    2*std::pow(1-v, 2.0f)*(1-(2*u)), -4*(u-std::pow(u, 2.0f))*(1-v),
    4*(v-std::pow(v, 2.0f))*(1-(2*u)), 4*(u-std::pow(u, 2.0f))*(1-(2*v)),
    2*std::pow(v, 2.0f)*(1-(2*u)), 4*(u-std::pow(u, 2.0f))*v,

    2*std::pow(1-v, 2.0f)*u, -2*std::pow(u, 2.0f)*(1-v),
    4*(v-std::pow(v, 2.0f))*u, 2*std::pow(u, 2.0f)*(1-(2*v)),
    2*std::pow(v, 2.0f)*u, 2*std::pow(u, 2.0f)*v;


    std::cout << "The matrix q is of size " << q.rows() << "x" << q.cols() << std::endl;
    std::cout << q << std::endl;


    Point3D p = pixel - _focal_point;
    Point3D f = surface.generatePoint(u, v);


    std::pair<int, int> order = surface.getOrder();
    std::vector<float> vector_of_r;
    vector_of_r.reserve(surface.getNumControlPoints());
    for (size_t i=0; i<=order.first; i++) // Fix order to be (n, m) = (2, 2) so we have 9 control points
    {
        for (size_t j=0; j<=order.second; j++)
        {
            Point3D k = surface.controlPoint(i, j);
            Point3D h = crossProduct(p, k);
            
            // c = p x (f(u, v) - foc) = (p x f(u, v)) - (p x foc)
            Point3D c = crossProduct(p, f) - crossProduct(p, _focal_point);

            vector_of_r.push_back(dotProduct(c, h));
        }
    }

    // Fill in r vector
    for (size_t n=0; n<vector_of_r.size(); n++)
    {
        r(n) = vector_of_r[n];
    }

    std::cout << "The vector r is of size " << r.size() << std::endl;
    std::cout << r << std::endl;

    Eigen::Vector2f gradient = 2.0*(r*q);


    return Point2D(gradient(0), gradient(1));
}

uint Camera::gradientDescent(Surface& surface, Point3D& pixel, std::pair<double, double>& target_u_v, float gamma, double epsilon)
{
    Point2D curr_u_v = Point2D(target_u_v.first, target_u_v.second);
    Point2D prev_u_v = curr_u_v;

    Point2D diff_u_v;
    double diff;
    uint num_iterations = 0;

    do
    {
        curr_u_v.printPoint2D();
        curr_u_v -= (calculateGradient(surface, pixel, std::make_pair(curr_u_v.xCoord(), curr_u_v.yCoord())) * gamma);
        diff_u_v = curr_u_v-prev_u_v;
        diff = dotProduct(diff_u_v, diff_u_v);
        std::cout << "diff: " << diff << std::endl;
        prev_u_v = curr_u_v;
        num_iterations++;
    } while (diff > epsilon);
    
    target_u_v.first = curr_u_v.xCoord();
    target_u_v.second = curr_u_v.yCoord();

    return num_iterations;
}

BoundingBox Camera::calculateBoundingBox(Surface& surface)
{
    /* Finds the 3-D bounding box for a particular surface */


    return BoundingBox();
}









float Camera::calcualteC(Surface& surface, Point3D& pixel, float u, float v)
{
    // Calculates |c|^2
    Point3D p = pixel - _focal_point;
    Point3D f = surface.generatePoint(u, v);
    Point3D c = crossProduct(p, f) - crossProduct(p, _focal_point);
    return dotProduct(c, c);
}


std::pair<float, float> Camera::validateGradient(Surface& surface, Point3D& pixel, float u, float v, float h)
{
    Point2D gradient = calculateGradient(surface, pixel, std::make_pair(u, v));


    float du = (calcualteC(surface, pixel, u+h, v) - calcualteC(surface, pixel, u-h, v)) / (2*h);
    float dv = (calcualteC(surface, pixel, u, v+h) - calcualteC(surface, pixel, u, v-h)) / (2*h);

    // float du = (calcualteC(surface, pixel, u+h, v) - calcualteC(surface, pixel, u, v)) / h;
    // float dv = (calcualteC(surface, pixel, u, v+h) - calcualteC(surface, pixel, u, v)) / h;

    return std::make_pair(std::abs(gradient.xCoord() - du), std::abs(gradient.yCoord() - dv));
}










int main()
{
    // World world;
    // world.insertSurface(Surface());

    // rasterize(world);

    Camera camera;
    camera.resolution(2, 2);
    camera.size(10);
    camera.direction(-1.634, .52652, 6.427);
    // camera.direction(.5, -.5, .5);
    camera.focalLength(10.424);
    

    Surface surface(2, 2);
    surface[{0, 0}] = Point3D(0, 0, 0);
    surface[{1, 0}] = Point3D(.5, 0, .5);
    surface[{2, 0}] = Point3D(1, 0, 0);
    surface[{0, 1}] = Point3D(0, .5, .5);
    surface[{1, 1}] = Point3D(.5, .5, 1);
    surface[{2, 1}] = Point3D(1, .5, .5);
    surface[{0, 2}] = Point3D(0, 1, 0);
    surface[{1, 2}] = Point3D(.5, 1, .5);
    surface[{2, 2}] = Point3D(1, 1, 0);








    // VALIDATES THE GRADIENT
    // uint num_of_iterations = 50000;
    // float mean_error_u = 0.0;
    // float mean_error_v = 0.0;
    // Point3D pixel(1, 2, 3);
    // for (size_t i=0; i<num_of_iterations; i++)
    // {
    //     float u = static_cast<float> (rand()) / static_cast<float> (RAND_MAX);
    //     float v = static_cast<float> (rand()) / static_cast<float> (RAND_MAX);
    //     // std::cout << u << " " << v << std::endl;
    //     std::pair<float, float> error_pair = camera.validateGradient(surface, pixel, u, v, .001);
    //     mean_error_u += error_pair.first;
    //     mean_error_v += error_pair.second;
    // }

    // std::cout << "Mean Error u: " << mean_error_u/num_of_iterations << std::endl;
    // std::cout << "Mean Error v: " << mean_error_v/num_of_iterations << std::endl;







    // RUNS GRADIENT DESCENT
    // Point3D pixel(1, 2, 3);
    // std::pair<double, double> u_v_pair(1.4, .2);
    // double gamma = .01;
    // double epsilon = .001;

    // uint iteration_count = camera.gradientDescent(surface, pixel, u_v_pair, gamma, epsilon);

    // std::cout << "Gradient Descent completed in " << iteration_count << " steps!" << std::endl;
    // std::cout << "u: " << u_v_pair.first << " v: " << u_v_pair.second << std::endl;







    camera.calculateIntersections(surface);


}
