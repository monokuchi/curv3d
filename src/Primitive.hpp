#ifndef PRIMITIVE_HPP
#define PRIMITIVE_HPP

#include <array>

#include "Utils.hpp"

#define GRID 9
#define MAX_ORDER 5

namespace curv3d {

class Camera; // Forward declaration.

/**
 * Class for Surfaces based on order.
 */
class Surface {
  private:
    real* const pPoints;  /**< Pointer to entire control point array. */
    uint* const pIndices; /**< Pointer to this Surface's indices. */

    HD real static bern_coeff(const real u, const uint i, const uint n);
    //    static constexpr std::array<std::array<real, GRID + 1>, GRID + 1> mQ = make_Q();

  public:
    const Order order;     /**< Stores the order. */
    const uint num_points; /**< Stores the number of control points. */
    const GPoint origin;   /**< The origin of the surface. */

    /**
     * Constructor for a Surface.
     * \param order the order of the surface.
     * \param the pointer to the points array.
     * \param the pointer to the start of the indices. */
    HD Surface(const Order& order, real* const points, uint* const indices, const GPoint& origin);

    HD LPoint control(const uint i);               /**< Returns the i-th control point. */
    HD LPoint control(const uint i, const uint j); /**< Returns the (i,j)-th control point. */
    HD LPoint operator[](const uint i);

    HD const LPoint control(const uint i) const;
    HD const LPoint control(const uint i, const uint j) const;
    HD const LPoint operator[](const uint i) const;

    /**
     * Evaluates the Bezier surface equation B(u,v) and returns the point
     * on the surface at (u,v).
     * \param uv the coordinate (u,v).
     */
    HD LPoint at(const SPoint& uv) const;

    /**
     * Gets the bounding box of the Surface on the Camera using the control points.
     * \param camera the camera.
     */
    BoundingBox basic_bounding_box(const Camera& camera) const;

    /**
     * Gets the bounding box of the Surface on the Camera using the local extrema.
     * \param camera the camera.
     */
    BoundingBox extrema_bounding_box(const Camera& camera) const;

    /**
     * Gets the alignment between a point (u,v) and a camera pixel.
     * \param uv the (u,v) pair.
     * \param focal the focal vector (position of the camera).
     * \param pixel the local pixel vector.
     */
    real alignment(const SPoint& uv, const GPoint& focal, const LPoint& pixel) const;

    /**
     * Gets the gradient of the alignment between a point (u,v) and a camera pixel.
     * \param uv the (u,v) pair.
     * \param focal the focal vector (position of the camera).
     * \param pixel the local pixel vector.
     */
    SPoint alignment_gradient(const SPoint& uv, const GPoint& focal, const LPoint& pixel) const;

    SPoint validate_gradient(const GPoint& focal, const LPoint& pixel) const;
};

};     // namespace curv3d
#endif // PRIMITIVE_HPP
