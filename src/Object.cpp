#include <cmath>

#include "Object.hpp"

namespace curv3d {

////////////////////
// Object Methods //
////////////////////

Object::Object(const GPoint& position, const Angles& rotation, const real scale)
    : mPosition(position), mRotation(rotation), mScale(Vector(scale)), mBasis(Basis::Identity()) {
    transform_points(mScale);
}

Object::Object(const GPoint& position, const Angles& rotation, const Vector& scale)
    : mPosition(position), mRotation(rotation), mScale(scale), mBasis(Basis::Identity()) {
    transform_points(mScale);
}

void Object::transform_points(const Vector& old_scale) {
    Vector cos(std::cos(mRotation[0]), std::cos(mRotation[1]), std::cos(mRotation[2]));
    Vector sin(std::sin(mRotation[0]), std::sin(mRotation[1]), std::sin(mRotation[2]));

    real sin01 = sin[0] * sin[1];

    Basis transform = mBasis;
    for (uint i = 0; i < 3; ++i) {
        Vector::from_pointer(transform.data() + 3 * i) /= old_scale[i];
    }
    transform.transposeInPlace(); // Gets the scaled inverse of the basis.

    // Gets the new basis.
    // Basis X
    mBasis(0, 0) = cos[0] * cos[1];
    mBasis(1, 0) = sin01;
    mBasis(2, 0) = -sin[1];
    // Basis Y
    mBasis(0, 1) = cos[0] * sin[1] * sin[2] - sin[0] * cos[2];
    mBasis(1, 1) = cos[0] * cos[2] + sin01 * sin[2];
    mBasis(2, 1) = cos[1] * sin[2];
    // Basis Z
    mBasis(0, 2) = sin[0] * sin[2] + cos[0] * sin[1] * cos[2];
    mBasis(1, 2) = sin01 * cos[2] - cos[0] * sin[2];
    mBasis(2, 2) = cos[1] * cos[2];

    for (uint i = 0; i < 3; ++i) { // Scale new basis.
        Vector::from_pointer(mBasis.data() + 3 * i) *= mScale[i];
    }
    transform *= mBasis; // Final transform.

    mPoints = transform * mPoints; // Transforms points.
}

void Object::scale_points(const Vector& old_scale) {
    Basis transform = Basis::Zero();
    for (uint i = 0; i < 3; ++i) {
        transform(i, i) *= mScale[i] / old_scale[i];
    }
    mPoints = transform * mPoints;
}

// Getters
GPoint Object::position() const { return mPosition; }
Angles Object::rotation() const { return mRotation; }
Vector Object::scale() const { return mScale; }

template <uint N> Vector Object::basis() const {
    return Vector(mBasis(0, N), mBasis(1, N), mBasis(2, N));
}
template Vector Object::basis<0>() const;
template Vector Object::basis<1>() const;
template Vector Object::basis<2>() const;

uint Object::num_points() const { return static_cast<uint>(mPoints.cols()); }
uint Object::num_surfaces() const { return static_cast<uint>(mOrders.size()); }

void Object::position(const GPoint& position) { mPosition = position; }

void Object::rotation(const Angles& rotation) {
    mRotation = rotation;
    transform_points(mScale);
}

void Object::scale(const Vector& scale) {
    Vector old_scale = mScale;
    mScale = scale;
    scale_points(old_scale);
}
void Object::scale(const real scale) { this->scale(Vector(scale)); }

template <> void Object::scale_axis<0>(const real scale) {
    this->scale(Vector(scale, mScale[1], mScale[2]));
}
template <> void Object::scale_axis<1>(const real scale) {
    this->scale(Vector(mScale[0], scale, mScale[2]));
}
template <> void Object::scale_axis<2>(const real scale) {
    this->scale(Vector(mScale[0], mScale[1], scale));
}

void Object::transform(const Angles& rotation, const Vector& scale) {
    Vector old_scale = mScale;
    mRotation = rotation, mScale = scale;
    transform_points(old_scale);
}
void Object::transform(const Angles& rotation, const real scale) {
    transform(rotation, Vector(scale));
}

void Object::move(const Vector& delta) { mPosition += delta; }
void Object::rotate(const Vector& delta) {
    mRotation += delta;
    transform_points(mRotation);
}

LPoint Object::point(const uint i) {
    ERROR_IF(i, >, num_points(), "Point index out of range");
    return LPoint::from_pointer(mPoints.data() + 3 * i);
}

const LPoint Object::point(const uint i) const {
    ERROR_IF(i, >, num_points(), "Point index out of range");
    return LPoint::from_pointer(const_cast<real*>(mPoints.data() + 3 * i));
}

void Object::add_point(const Vector& v) {
    uint index = static_cast<uint>(mPoints.cols());
    mPoints.conservativeResize(3, mPoints.cols() + 1);

    Eigen::Vector3d new_v(v.x, v.y, v.z);
    new_v = mBasis * new_v;

    mPoints(0, index) = new_v(0);
    mPoints(1, index) = new_v(1);
    mPoints(2, index) = new_v(2);
}

void Object::add_surface(const Order& order, const ControlList& list) {
    auto list_x = list.size(), list_y = list[0].size();
    ERROR_IF(list_x, !=, order[0] + 1, "Order and control list mismatch");
    ERROR_IF(list_y, !=, order[1] + 1, "Order and control list mismatch");

    if (num_surfaces()) {
        Order prev_order = mOrders.back();
        uint prev_num_points = (prev_order[0] + 1) * (prev_order[1] + 1);
        mOffsets.push_back(mOffsets.back() + prev_num_points);
    } else {
        mOffsets.push_back(0);
    }

    mOrders.push_back(order);
    for (auto it_x = list.begin(); it_x != list.end(); ++it_x) {
        for (auto it_y = it_x->begin(); it_y != it_x->end(); ++it_y) {
            mIndices.push_back(*it_y);
        }
    }
}

Surface Object::surface(const uint i) {
    ERROR_IF(i, >=, num_surfaces(), "Surface index out of bounds");
    return Surface(mOrders[i], mPoints.data(), &mIndices[mOffsets[i]], mPosition);
}

const Surface Object::surface(const uint i) const {
    ERROR_IF(i, >=, num_surfaces(), "Surface index out of bounds");
    return Surface(mOrders[i], const_cast<real*>(mPoints.data()),
                   const_cast<uint*>(&mIndices[mOffsets[i]]), mPosition);
}

////////////////////
// Camera Methods //
////////////////////

Camera::Camera(const GPoint& position, const TwoIndex& resolution, const real width,
               const real focal_length, const Angles& rotation)
    : Object(position, rotation,
             Vector(focal_length, -width / resolution[0], width / resolution[0])),
      resolution(resolution) {
    add_point(
        LPoint(1, -static_cast<real>(resolution[0]) / 2, -static_cast<real>(resolution[1]) / 2));
};

real Camera::focal_length() const { return basis<0>().norm(); }
Vector Camera::direction() const { return basis<0>(); }
real Camera::width() const { return basis<1>().norm() * resolution[0]; }

ScreenBasis Camera::screen_basis() const { return std::make_pair(basis<1>(), basis<2>()); }
LPoint Camera::screen_origin() const { return point(0); }

void Camera::focal_length(const real length) { scale_axis<0>(length); }
void Camera::width(const real width) {
    real pix_size = width / resolution[0];
    scale(Vector(mScale[0], pix_size, pix_size));
}

SPoint Camera::project(const Vector& vec) const {
    Vector w = vec - mPosition;
    real f2 = mScale[0] * mScale[0];
    real correction = f2 / (f2 + w.dot(basis<0>()));
    SPoint proj = SPoint(w.dot(basis<1>()), w.dot(basis<2>())) / pow(mScale[1], 2);
    return correction * proj + SPoint(resolution[0] / 2, resolution[1] / 2);
}

LPoint Camera::pixel(const uint i, const uint j) const {
    return screen_origin() + basis<1>() * (i + 0.5) + basis<2>() * (j + 0.5);
}

}; // namespace curv3d
