
#include <rasterize.h>



float calculateIntersection(Camera camera, Surface surface, Pixel pixel)
{
    return 2;
}

float newtonsMethod()
{
    return 2;
}

void calculateBoundingBox()
{

}

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

std::pair<Point3D, Point3D> Camera::getScreenBasis()
{
    // float f_x = _orientation_vector.xCoord();
    // float f_y = _orientation_vector.yCoord();
    // float f_z = _orientation_vector.zCoord();


    // Basis vector of screen space with respect to camera space
    Point3D u(_orientation_vector.yCoord(), -_orientation_vector.xCoord(), 0);

    // If we want to normalize u to width / _resolution_x we need to multiply u by width / (mag(u) * _resolution_x)
    float normalizing_factor = _size / (u.mag() * _resolution_x);
    u *= normalizing_factor;


    Point3D f = _orientation_vector;
    f *= 1 / _orientation_vector.mag();


    // float u_x = u.xCoord();
    // float u_y = u.yCoord();
    // float u_z = u.zCoord();

    // v = u x f
    // Point3D v(u_y*f_z - u_z*f_y, u_z*f_x - f_z*u_x, u_x*f_y - u_y*f_x);
    Point3D v = crossProduct(u, f);



    Point3D rotated_u = (u*std::cos(_orientation_angle)) + (v*std::sin(_orientation_angle));

    // float rotated_u_x = rotated_u.xCoord();
    // float rotated_u_y = rotated_u.yCoord();
    // float rotated_u_z = rotated_u.zCoord();
    
    // rotated_v = rotated_u x f
    // Point3D rotated_v(rotated_u_y*f_z - rotated_u_z*f_y, rotated_u_z*f_x - f_z*rotated_u_x, rotated_u_x*f_y - rotated_u_y*f_x);
    Point3D rotated_v = crossProduct(rotated_u, f);

    return std::make_pair(rotated_u, rotated_v);
}

void Camera::calculateIntersections(Point3D pixel_to_camera, Surface& surface)
{
    std::pair<Point3D, Point3D> basis_vectors = getScreenBasis();


    // Implement bounding box later
    for (size_t i=0; i<_resolution_x; i++)
    {
        for (size_t j=0; j<_resolution_y; j++)
        {
            Point3D screen_point = basis_vectors.first*(i-_resolution_x/2) + basis_vectors.second*(j-_resolution_y/2);

            Point3D p = _orientation_vector + screen_point;
            

            // min |p x (f(u,v) - foc)|^2, subject to u,v in [0,1]^2 
            // Assume order of the bezier surface, f(u, v), is fixed at 2
            


        }
    }
}

Point2D Camera::calculateGradient(Surface& surface, float u, float v)
{
    /* 
        Calculates gradient of |p x (f(u, v) - foc)|^2
        p: orientation vector
        f(u, v): 3-D point on bezier surface
        foc: focal point of camera
    */

    Eigen::RowVectorXf r(surface.getNumControlPoints());
    Eigen::MatrixXf q(surface.getNumControlPoints(), 2);


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


    Point3D f = surface.generatePoint(u, v);


    std::pair<int, int> order = surface.getOrder();
    std::vector<float> vector_of_r;
    vector_of_r.reserve(surface.getNumControlPoints());
    for (size_t i=0; i<=order.first; i++) // Fix order to be (n, m) = (2, 2) so we have 9 control points
    {
        for (size_t j=0; j<=order.second; j++)
        {
            Point3D k = surface.controlPoint(i, j);
            Point3D h = crossProduct(_orientation_vector, k);
            
            // c = p x (f(u, v) - foc) = (p x f(u, v)) - (p x foc)
            Point3D c = crossProduct(_orientation_vector, f) - crossProduct(_orientation_vector, _focal_point);

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

    Eigen::Vector2f gradient = 2*(r*q);


    return Point2D(gradient(0), gradient(1));
}

uint Camera::gradientDescent(Point2D& point, Surface& surface, float gamma, float u, float v)
{
    Point2D prev_point = point;

    double diff;
    uint num_iterations = 0;

    do
    {
        point.printPoint2D();
        point -= calculateGradient(surface, u, v) * gamma;
        diff = dotProduct(point, prev_point);
        std::cout << "diff: " << diff << std::endl;
        num_iterations++;
    } while (diff > std::numeric_limits<double>::epsilon());
    

    return num_iterations;
}




float Camera::calcualteC(Surface& surface, float u, float v)
{
    // Calculates |c|^2
    Point3D f = surface.generatePoint(u, v);
    Point3D c = crossProduct(_orientation_vector, f) - crossProduct(_orientation_vector, _focal_point);
    return dotProduct(c, c);
}


std::pair<float, float> Camera::validateGradient(Surface& surface, float u, float v, float h)
{
    Point2D gradient = calculateGradient(surface, u, v);


    float du = (calcualteC(surface, u+h, v) - calcualteC(surface, u-h, v)) / (2*h);
    float dv = (calcualteC(surface, u, v+h) - calcualteC(surface, u, v-h)) / (2*h);

    // float du = (calcualteC(surface, u+h, v) - calcualteC(surface, u, v)) / h;
    // float dv = (calcualteC(surface, u, v+h) - calcualteC(surface, u, v)) / h;

    return std::make_pair(std::abs(gradient.xCoord() - du), std::abs(gradient.yCoord() - dv));
}










int main()
{
    // World world;
    // world.insertSurface(Surface());

    // rasterize(world);

    Camera camera;
    camera.direction(1.634, -.52652, 6.427);
    // camera.direction(.5, .5, .5);
    camera.focalLength(19);
    

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










    // RUNS GRADIENT DESCENT
    Point2D gradient_point(1, 1);
    uint iteration_count = camera.gradientDescent(gradient_point, surface, .001, .3, .5);

    std::cout << "Gradient Descent completed in " << iteration_count << " steps!" << std::endl;
    gradient_point.printPoint2D();






    // VALIDATES THE GRADIENT
    // uint num_of_iterations = 50000;
    // float mean_error_u = 0.0;
    // float mean_error_v = 0.0;
    // for (size_t i=0; i<num_of_iterations; i++)
    // {
    //     float u = static_cast<float> (rand()) / static_cast<float> (RAND_MAX);
    //     float v = static_cast<float> (rand()) / static_cast<float> (RAND_MAX);
    //     // std::cout << u << " " << v << std::endl;
    //     std::pair<float, float> error_pair = camera.validateGradient(surface, u, v, .001);
    //     mean_error_u += error_pair.first;
    //     mean_error_v += error_pair.second;
    // }

    // std::cout << "Mean Error u: " << mean_error_u/num_of_iterations << std::endl;
    // std::cout << "Mean Error v: " << mean_error_v/num_of_iterations << std::endl;
}
