#ifndef SCENE_HPP
#define SCENE_HPP

#include "Object.hpp"

namespace curv3d {

/**
 * Struct for the Pixel data.
 */
struct Pixel {
    uint8_t r, g, b;

    Pixel(const uint8_t r, const uint8_t g, const uint8_t b) : r(r), g(g), b(b){};
    Pixel(volatile Pixel& pix) : r(pix.r), g(pix.g), b(pix.b){};
    void operator=(const Pixel& pix) volatile;

    template <uint N> uint8_t& operator()();
    template <uint N> const uint8_t& operator()() const;

    static volatile Pixel& from_pointer(volatile uint8_t* const p);
};

struct ImageBuffer {
    TwoIndex resolution;
    volatile uint8_t* pointer;

    ImageBuffer(const TwoIndex& resolution);
    ~ImageBuffer();

    volatile Pixel& operator()(const uint i, const uint j);
};

/**
 * Class that contains the information for one scene or world.
 */
class Scene {
  private:
    std::vector<Object> mObjects;
    Camera mCamera;
    ImageBuffer* mBuffer;

  public:
    Scene(const Camera& camera, ImageBuffer* const buffer) : mCamera(camera), mBuffer(buffer){};

    uint num_objects() const; /**< Gets the number of objects. */

    Camera& camera(); /**< Returns a reference to the Camera.*/
    const Camera& camera() const;

    Object& object(const uint i); /**< Returns a reference to the i-th object. */
    const Object& object(const uint i) const;

    void add_object(const Object& object); /**< Adds an object to the Scene.*/

    /**
     * Rasterizes the image slowly using CPU and looping through each bounding box.
     */
    std::vector<uint> slow_rasterize();

    /**
     * Rasterizes the image by pre-sorting Surfaces into Tiles.
     */
    void tile_rasterize();
};

};     // namespace curv3d
#endif // SCENE_HPP
