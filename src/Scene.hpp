#ifndef SCENE_HPP
#define SCENE_HPP

#include "Object.hpp"

namespace curv3d {

/**
 * Struct for the Pixel data.
 */
struct Pixel {
    uint8_t r, g, b; /**< Stores the R, G, B values. */

    /**
     * Constructor for a Pixel.
     * \param r the red value.
     * \param g the green value.
     * \param b the blue value.
     */
    Pixel(const uint8_t r, const uint8_t g, const uint8_t b) : r(r), g(g), b(b){};
    Pixel(volatile Pixel& pix) : r(pix.r), g(pix.g), b(pix.b){};
    void operator=(const Pixel& pix) volatile;

    uint8_t& operator[](const uint i); /**< Returns the i-th (RGB) color. */
    const uint8_t& operator[](const uint i) const;

    /**
     * Recasts a pointer into a Pixel at that location.
     * \param p the pointer at that location.
     */
    static volatile Pixel& from_pointer(volatile uint8_t* const p);
};

struct ImageBuffer {
    TwoIndex resolution;       /**< The resolution for this buffer.*/
    volatile uint8_t* pointer; /**< The pointer to the actual Pixel data. */

    /**
     * Constructor for an ImageBuffer.
     * \param resolution the desired resolution.
     */
    ImageBuffer(const TwoIndex& resolution);
    ~ImageBuffer();

    /**
     * Returns a reference to the  (i,j)th pixel of the ImageBuffer.
     */
    volatile Pixel& operator()(const uint i, const uint j);
};

/**
 * Class that contains the information for one scene or world.
 */
class Scene {
  private:
    std::vector<Object> mObjects; /**< Stores the objects. */
    Camera mCamera;               /**< The associated camera. */
    ImageBuffer* mBuffer;         /**< Pointer to the ImageBuffer. */

  public:
    /**
     * Constructor for a Scene.
     * \param camera the camera for the scene.
     * \param buffer the pointer to the ImageBuffer.
     */
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
