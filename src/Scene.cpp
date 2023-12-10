#include "Scene.hpp"

namespace curv3d {

///////////////////
// Pixel Methods //
///////////////////

uint8_t& Pixel::operator[](const uint i) {
    ERROR_IF(i, >, 2, "Color number must be [0:2]");
    return *(&r + i);
}
const uint8_t& Pixel::operator[](const uint i) const {
    ERROR_IF(i, >, 2, "Color number must be [0:2]");
    return *(&r + i);
}

void Pixel::operator=(const Pixel& pix) volatile { r = pix.r, g = pix.g, b = pix.b; }

volatile Pixel& Pixel::from_pointer(volatile uint8_t* const p) {
    return *reinterpret_cast<volatile Pixel*>(p);
}

/////////////////////////
// ImageBuffer Methods //
/////////////////////////
ImageBuffer::ImageBuffer(const TwoIndex& resolution)
    : resolution(resolution), pointer(new volatile uint8_t[3 * resolution[0] * resolution[1]]) {
    memset(const_cast<uint8_t*>(pointer), 0u, 3 * resolution[0] * resolution[1] * sizeof(uint8_t));
};

ImageBuffer::~ImageBuffer() { delete[] pointer; }

volatile Pixel& ImageBuffer::operator()(const uint i, const uint j) {
    return Pixel::from_pointer(pointer + 3 * resolution[0] * i + j);
}

///////////////////
// Scene Methods //
///////////////////

uint Scene::num_objects() const { return static_cast<uint>(mObjects.size()); }

Camera& Scene::camera() { return mCamera; }
const Camera& Scene::camera() const { return mCamera; }

Object& Scene::object(const uint i) {
    ERROR_IF(i, >=, num_objects(), "Object index out of bounds");
    return mObjects[i];
}
const Object& Scene::object(const uint i) const {
    ERROR_IF(i, >=, num_objects(), "Object index out of bounds");
    return mObjects[i];
}

void Scene::add_object(const Object& object) { mObjects.push_back(object); }

std::pair<TwoIndex, TwoIndex> bounding_box_indices(const Camera& camera, const BoundingBox& bb) {
    TwoIndex min(std::max(0, static_cast<int>(bb.min[0])),
                 std::max(0, static_cast<int>(bb.min[1])));
    TwoIndex max(std::min(static_cast<int>(camera.resolution[0]), static_cast<int>(bb.max[0])),
                 std::min(static_cast<int>(camera.resolution[1]), static_cast<int>(bb.max[1])));
    return std::make_pair(min, max);
}

uint gradient_descent(SPoint& uv, const Surface& surface, const GPoint& focal, const LPoint& pixel,
                      const real gamma, const real eps) {
    uint max_iterations = 100;
    SPoint grad;
    for (uint iIter = 0; iIter < max_iterations; ++iIter) {
        grad = surface.alignment_gradient(uv, focal, pixel);
        uv -= gamma * grad;
        if (grad.dot(grad) <= eps) return iIter + 1;
    }
    return 0;
}

std::vector<uint> Scene::slow_rasterize() {
    std::vector<uint> iters;

    // No iterators for later parallelism.
    for (uint iObject = 0; iObject < num_objects(); ++iObject) {
        Object obj = object(iObject);
        for (uint iSurface = 0; iSurface < obj.num_surfaces(); ++iSurface) {
            Surface surf = obj.surface(iSurface);

            auto [min, max] = bounding_box_indices(mCamera, surf.basic_bounding_box(mCamera));

            PRINT(min << " " << max);

            for (uint iX = min[0]; iX <= max[0]; ++iX) {
                for (uint iY = min[1]; iY <= max[1]; ++iY) {
                    LPoint pixel = mCamera.pixel(iX, iY);

                    // PRINT("Error " << surf.validate_gradient(mCamera.position(), pixel));

                    // min |p x (f(u,v) - foc)|^2, subject to u,v in [0,1]^2
                    // Assume order of the bezier surface, f(u, v), is fixed at 2
                    // epsilon = |(pixel(i,j) - foc) x (pixel(i+.5, j))|^2
                    real gamma = 0.1;
                    Vector eps_v = pixel.cross(pixel + mCamera.screen_basis().first / 2);
                    real eps = eps_v.dot(eps_v);

                    SPoint uv(0.5, 0.5); // Initial guess.
                    uint iter = gradient_descent(uv, surf, mCamera.position(), pixel, gamma, eps);

                    iters.push_back(iter);
                    // PRINT("uv: " << uv);
                    if (iter) (*mBuffer)(iX, iY) = Pixel(255, 255, 255);
                }
            }
        }
    }
    return iters;
}
}; // namespace curv3d
