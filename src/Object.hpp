#ifndef OBJECT_HPP
#define OBJECT_HPP

#include <iostream>
#include <vector>
#include <Eigen/Dense>

#include "Primitive.hpp"
#include "Utils.hpp"

namespace curv3d {

using Basis = Eigen::Matrix<real, 3, 3, Eigen::ColMajor>; /**< Represents a 3D Basis. */
/**
 * Represents a 3xN matrix of points.
 */
using PointMatrix = Eigen::Matrix<real, 3, Eigen::Dynamic, Eigen::ColMajor>;
using ScreenBasis = std::pair<Vector, Vector>; /**< Represents the basis of a camera screen. */

/**
 * Generic class for objects, which are composed of Surfaces.
 */
class Object {
  protected:
    GPoint mPosition; /** The global position of this Object. */
    Angles mRotation; /** The 'spherical' angles of the object. */
    Vector mScale;    /**< The scaling of the object. */
    Basis mBasis;     /**< The basis of the object, defined by the rotation and scale. */

    std::vector<Order> mOrders;       /**< Stores the orders of the Surfaces of this Object. */
    PointMatrix mPoints;              /**< Stores the points this object owns. */
    std::vector<uint> mIndices;       /**< Stores the indices of the control points for Surfaces. */
    std::vector<uint> mOffsets = {0}; /**< Stores the index offsets. */

    /**
     * Updates the points based on the object parameters.
     */
    void transform_points(const Vector& old_scale);

    /**
     * Updates the basis based only on the object scale.
     */
    void scale_points(const Vector& old_scale);

  public:
    /**
     * Constructs an generic object.
     * \param position the position of the object. Default is (0,0,0).
     * \param rotation the rotation of the object. Default is (0,0,0).
     * \param scale the scale of the object. Default is (1,1,1).
     */
    Object(const GPoint& position = GPoint(0.0f), const Angles& rotation = Angles(0.0f),
           const real scale = 1.0f);
    Object(const GPoint& position, const Angles& rotation, const Vector& scale);

    HD GPoint position() const; /**< Gets the position of this Object. */
    HD Angles rotation() const; /**< Gets the rotation of this Object. */
    HD Vector scale() const;    /**< Gets the scale of this Object. */

    HD template <uint N> Vector basis() const;
    uint num_points() const;   /**< Gets the number of points. */
    uint num_surfaces() const; /**< Gets the number of surfaces. */

    HD void position(const GPoint& position); /**< Sets the position of this Object. */
    void rotation(const Angles& rotation);    /**< Sets the rotation of this Object. */
    void scale(const Vector& scale);          /**< Sets the scale of this Object. */
    void scale(const real scale);             /**< Sets the uniform scale of this Object. */
    template <uint N>
    void scale_axis(const real scale); /**< Sets the scale of a given axis of this Object. */

    /**
     * Transforms this object with a given rotation and scale.
     * \param rotation the new rotation of this object.
     * \param scale the new scale of this object.
     */
    void transform(const Angles& rotation, const Vector& scale);
    void transform(const Angles& rotation, const real scale);

    HD void move(const Vector& delta); /**< Moves this object with offset delta. */
    void rotate(const Vector& delta);  /**< Rotates this object with offset delta. */

    LPoint point(const uint i); /**< Gets the i-th point. */
    const LPoint point(const uint i) const;
    void add_point(const Vector& v); /**< Adds a point to the PointMatrix. */

    /**
     * Adds a surface to this Object, given an order and list of indices of the
     * control points.
     * \param order the order of the surface.
     * \param list the list of the indices of the control points.
     */
    void add_surface(const Order& order, const ControlList& list);

    Surface surface(const uint i); /**< Returns a pointer to the i-th surface. */
    const Surface surface(const uint i) const;
};

class Camera : public Object {
  public:
    const TwoIndex resolution; /**< The screen resolution (pixels) of the camera. */

    Camera(const GPoint& position, const TwoIndex& resolution, const real width,
           const real focal_length, const Angles& rotation = Angles(0.0f));

    real focal_length() const; /**< Gets the focal_length. */
    Vector direction() const;  /**< Gets the direction of the camera. */
    real width() const;        /**< Gets the camera width's in physical space. */

    LPoint screen_origin() const;     /**< Gets the origin of the screen. */
    ScreenBasis screen_basis() const; /**< Gets the basis of the screen. */

    void focal_length(const real length); /**< Sets the focal length. */
    void width(const real width);         /**< Sets the camera width. */

    /**
     * Projects a vector onto the screen of the camera.
     */
    SPoint project(const Vector& vec) const;

    /**
     * Gets the location of the center of the (i,j)-th pixel.
     */
    LPoint pixel(const uint i, const uint j) const;
};

};     // namespace curv3d
#endif // OBJECT_HPP
