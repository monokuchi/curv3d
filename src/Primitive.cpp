#include "Primitive.hpp"
#include "Object.hpp"

namespace curv3d {

static const uint CHOOSE_CACHE_SIZE = (MAX_ORDER + 1) * (MAX_ORDER + 2) / 2 - 1;

// Generates the array for caching the choose functions (n k).
// where (n k), from n = 1 to MAX_ORDER, k = 0 to n, stored consecutively.
constexpr std::array<real, CHOOSE_CACHE_SIZE> make_choose_cache() {
    std::array<real, CHOOSE_CACHE_SIZE> cache{};
    uint start = 0;
    for (uint iOrder = 1; iOrder <= MAX_ORDER; ++iOrder) {
        real* arr = &cache[start];
        arr[0] = 1, arr[iOrder] = 1;
        for (uint i = 1, j = iOrder - 1; i <= j; ++i, --j) {
            real val = arr[i - 1] * (iOrder + 1 - i) / i;
            arr[i] = val, arr[j] = val;
        }
        start += iOrder + 1;
    }
    return cache;
}

static constexpr std::array<real, CHOOSE_CACHE_SIZE> CHOOSE_CACHE = make_choose_cache();
template <uint N> real choose_n_cache(const uint k) {
    return CHOOSE_CACHE[N * (N + 1) / 2 - 1 + k];
};

template <size_t... Ints> real choose(std::index_sequence<Ints...>, const uint n, const uint k) {
    std::function<real(const uint)> choose[] = {&choose_n_cache<Ints>...}; // Force instantiation.
    return choose[n](k);
};

real choose(const uint n, const uint k) {
    return choose(std::make_index_sequence<MAX_ORDER + 1>{}, n, k);
};

/////////////////////
// Surface Methods //
/////////////////////

HD real Surface::bern_coeff(const real u, const uint i, const uint n) {
    return pow(u, i) * pow(1 - u, n - i);
}

HD Surface::Surface(const Order& order, real* const points, uint* const indices,
                    const GPoint& origin)
    : pPoints(points), pIndices(indices), order(order), num_points((order[0] + 1) * (order[1] + 1)),
      origin(origin) {}

HD LPoint Surface::control(const uint i) {
    ERROR_IF(i, >, num_points, "Exceeded total number of control points");
    return LPoint::from_pointer(pPoints + 3 * pIndices[i]);
}

HD LPoint Surface::control(const uint i, const uint j) { return control(i * (order[0] + 1) + j); };

HD LPoint Surface::operator[](const uint i) { return control(i); }

HD const LPoint Surface::control(const uint i) const {
    ERROR_IF(i, >, num_points, "Exceeded total number of control points");
    return LPoint::from_pointer(pPoints + 3 * pIndices[i]);
}
HD const LPoint Surface::control(const uint i, const uint j) const {
    return control(i * (order[0] + 1) + j);
};
HD const LPoint Surface::operator[](const uint i) const { return control(i); }

HD LPoint Surface::at(const SPoint& uv) const {
    LPoint point = origin;
    for (uint i = 0; i <= order[0]; ++i) {
        LPoint temp(0);
        for (uint j = 0; j <= order[1]; ++j) {
            real bernsteinM = choose(order[1], j) * bern_coeff(uv[1], j, order[1]);
            temp += bernsteinM * control(i, j);
        }
        real bernsteinN = choose(order[0], i) * bern_coeff(uv[0], i, order[0]);
        point += bernsteinN * temp;
    }
    return point;
};

BoundingBox Surface::basic_bounding_box(const Camera& camera) const {
    BoundingBox bb;
    for (uint iPoint = 0; iPoint < num_points; ++iPoint) {
        SPoint proj = camera.project(control(iPoint) + origin);
        for (uint i = 0; i < 2; ++i) {
            if (proj[i] < bb.min[i]) bb.min[i] = proj[i];
            if (proj[i] > bb.max[i]) bb.max[i] = proj[i];
        }
    }
    return bb;
}

BoundingBox Surface::extrema_bounding_box(const Camera& camera) const {
    BoundingBox bb;
    // PointMatrix extrema(num_points);

    // Get new points.
    //
    // for (uint iPoint = 0; iPoint < {points}; ++iPoint) {
    //     SPoint proj = camera.project(extrema[i] + origin);
    //     for (uint i = 0; i < 2; ++i) {
    //         if (proj[i] < bb.min[i]) bb.min[i] = proj[i];
    //         if (proj[i] > bb.max[i]) bb.max[i] = proj[i];
    //     }
    // }
    return bb;

    // Eigen::RowVectorXf k(surface.getNumControlPoints());
    // Eigen::MatrixXf q(surface.getNumControlPoints(), 2);

    // double u = u_v_pair.first;
    // double v = u_v_pair.second;

    // // Since we fixed order to be (n, m) = (2, 2) we can manually fill in q
    // q << -2 * (1 - u) * std::pow(1 - v, 2.0f), -2 * (1 - v) * std::pow(1 - u, 2.0f),
    //     -4 * (v - std::pow(v, 2.0f)) * (1 - u), 2 * std::pow(1 - u, 2.0f) * (1 - (2 * v)),
    //     -2 * std::pow(v, 2.0f) * (1 - u), 2 * std::pow(1 - u, 2.0f) * v,

    //     2 * std::pow(1 - v, 2.0f) * (1 - (2 * u)), -4 * (u - std::pow(u, 2.0f)) * (1 - v),
    //     4 * (v - std::pow(v, 2.0f)) * (1 - (2 * u)), 4 * (u - std::pow(u, 2.0f)) * (1 - (2 * v)),
    //     2 * std::pow(v, 2.0f) * (1 - (2 * u)), 4 * (u - std::pow(u, 2.0f)) * v,

    //     2 * std::pow(1 - v, 2.0f) * u, -2 * std::pow(u, 2.0f) * (1 - v),
    //     4 * (v - std::pow(v, 2.0f)) * u, 2 * std::pow(u, 2.0f) * (1 - (2 * v)),
    //     2 * std::pow(v, 2.0f) * u, 2 * std::pow(u, 2.0f) * v;

    // std::vector<Point3D> control_points = surface.getControlPoints();

    // for (size_t i = 0; i < control_points.size(); i++) {
    //     // k(i) = control_points[i];
    //     control_points[i].printPoint3D();
    // }
}

real Surface::alignment(const SPoint& uv, const GPoint& focal, const LPoint& pixel) const {
    Vector v = pixel.cross(at(uv) + origin - focal);
    return v.dot(v);
}

SPoint Surface::alignment_gradient(const SPoint& uv, const GPoint& focal,
                                   const LPoint& pixel) const {
    /*
        Calculates gradient of |p x (f(u, v) - foc)|^2
        p: Pixel(i,j) - foc
        f(u, v): 3-D point on bezier surface
        foc: focal point of camera
    */
    ERROR_IF(order, !=, Order(2, 2), "Cannot compute gradient for non 2nd order");

    Eigen::RowVectorXf r(num_points);
    Eigen::MatrixXf q(num_points, 2);

    real u = uv[0], v = uv[1];

    // Since we fixed order to be (n, m) = (2, 2) we can manually fill in q
    q << -2 * (1 - u) * std::pow(1 - v, 2.0f), -2 * (1 - v) * std::pow(1 - u, 2.0f),
        -4 * (v - std::pow(v, 2.0f)) * (1 - u), 2 * std::pow(1 - u, 2.0f) * (1 - (2 * v)),
        -2 * std::pow(v, 2.0f) * (1 - u), 2 * std::pow(1 - u, 2.0f) * v,

        2 * std::pow(1 - v, 2.0f) * (1 - (2 * u)), -4 * (u - std::pow(u, 2.0f)) * (1 - v),
        4 * (v - std::pow(v, 2.0f)) * (1 - (2 * u)), 4 * (u - std::pow(u, 2.0f)) * (1 - (2 * v)),
        2 * std::pow(v, 2.0f) * (1 - (2 * u)), 4 * (u - std::pow(u, 2.0f)) * v,

        2 * std::pow(1 - v, 2.0f) * u, -2 * std::pow(u, 2.0f) * (1 - v),
        4 * (v - std::pow(v, 2.0f)) * u, 2 * std::pow(u, 2.0f) * (1 - (2 * v)),
        2 * std::pow(v, 2.0f) * u, 2 * std::pow(u, 2.0f) * v;

    for (uint iPoint = 0; iPoint < num_points; ++iPoint) {
        Vector c = pixel.cross(at(uv) + origin - focal);
        r(iPoint) = c.dot(pixel.cross(control(iPoint)));
    }

    Eigen::Vector2f gradient = 2.0 * (r * q);
    return SPoint(gradient(0), gradient(1));
}

SPoint Surface::validate_gradient(const GPoint& focal, const LPoint& pixel) const {
    uint divisions = 100;
    real h = 0.0001;

    SPoint error(0, 0);
    for (uint i = 0; i < divisions; ++i) {
        for (uint j = 0; j < divisions; ++j) {
            SPoint uv(i / (real)divisions, j / (real)divisions);

            SPoint grad = alignment_gradient(uv, focal, pixel);
            SPoint approx(alignment(uv + SPoint(h, 0), focal, pixel) -
                              alignment(uv - SPoint(h, 0), focal, pixel),
                          alignment(uv + SPoint(0, h), focal, pixel) -
                              alignment(uv - SPoint(0, h), focal, pixel));
            approx /= (2 * h);
            error += (grad - approx);
        }
    }
    error /= divisions;
    return error;
}

}; // namespace curv3d
