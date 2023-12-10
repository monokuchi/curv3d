#ifndef UTILS_HPP
#define UTILS_HPP

#include <chrono>
#include <cstdint>
#include <limits>
#include <vector>

#include <CudaTools/Array.h>
#include <CudaTools/Core.h>
#include <CudaTools/Macros.h>
#include <CudaTools/Types.h>

namespace CT = CudaTools;
using real = CT::Types::real64; /**< Stores real values. */
using uint = uint32_t;          /**< Stores an unsigned integer. */
using Shape = CT::Shape;

template <typename T> using Array = CT::Array<T>;

#define TIME_START(name) auto begin_##name = std::chrono::steady_clock::now()

#define TIME_END(name)                                                                             \
    auto end_##name = std::chrono::steady_clock::now();                                            \
    auto time_ms_##name =                                                                          \
        std::chrono::duration_cast<std::chrono::milliseconds>(end_##name - begin_##name).count();  \
    auto time_mus_##name =                                                                         \
        std::chrono::duration_cast<std::chrono::microseconds>(end_##name - begin_##name).count();  \
    if (time_ms_##name == 0) {                                                                     \
        printf("[%s] Time Elapsed: %ld[µs]\n", #name, time_mus_##name);                            \
    } else {                                                                                       \
        printf("[%s] Time Elapsed: %ld[ms]\n", #name, time_ms_##name);                             \
    }

#define TIME(call, name)                                                                           \
    TIME_START(name);                                                                              \
    call;                                                                                          \
    TIME_END(name);

#define ERROR_IF(a, op, b, msg)                                                                    \
    if (a op b) {                                                                                  \
        std::ostringstream os_a;                                                                   \
        std::ostringstream os_b;                                                                   \
        os_a << a;                                                                                 \
        os_b << b;                                                                                 \
        printf("\033[1;31m[curv3d]\033[0m %s:%d\n | %s: (" #a ")%s " #op " (" #b ")%s.\n",         \
               __FILE__, __LINE__, msg, os_a.str().c_str(), os_b.str().c_str());                   \
        throw std::exception();                                                                    \
    }

#define ERROR(a, msg)                                                                              \
    if (a) {                                                                                       \
        printf("\033[1;31m[curv3d]\033[0m %s:%d\n | %s: " #a ".\n", __FILE__, __LINE__, msg);      \
        throw std::exception();                                                                    \
    }

#define PRINT(items) std::cout << items << std::endl;

namespace curv3d {

real constexpr pow(real base, uint exp) { return (!exp) ? 1 : base * pow(base, exp - 1); }

// For the reinterpret_cast to work correctly, both TwoVector and ThreeVector
// need to be a 'standard-layout' class.

/**
 * A class for a mathematical vector in \f$ \mathbb{R}^2 \f$, supporting the
 * standard operations.
 */
template <typename T> struct TwoVector {
  public:
    template <typename U> friend std::ostream& operator<<(std::ostream&, const TwoVector<U>&);
    T x, y;

    HD TwoVector() = default;
    /**
     * Constructs a TwoVector with provided coordinates.
     * \param x the x coordinate.
     * \param y the y coordinate.
     */
    HD TwoVector(const T x, const T y) : x(x), y(y){};
    /**
     * Constructs a TwoVector with all same provided value.
     * \param x the value in all components.
     */
    HD TwoVector(const T x) : x(x), y(x){};

    // Due to memory layout, dereferencing needs these offsets.

    /**
     * Recasts a pointer into a TwoVector at that pointer.
     * \param p the pointer to the data.
     */
    HD static TwoVector& from_pointer(T* const p) { return *reinterpret_cast<TwoVector<T>*>(p); };

    HD T& operator[](const uint i) {
        ERROR_IF(i, >, 1, "Component number must be [0:1]");
        return *(&x + i);
    }
    HD const T& operator[](const uint i) const {
        ERROR_IF(i, >, 1, "Component number must be [0:1]");
        return *(&x + i);
    }

    HD TwoVector operator+(const TwoVector& v) const { return TwoVector(x + v.x, y + v.y); };
    HD TwoVector operator-(const TwoVector& v) const { return TwoVector(x - v.x, y - v.y); };
    HD TwoVector operator*(const T c) const { return TwoVector(x * c, y * c); };
    HD TwoVector operator/(const T c) const { return TwoVector(x / c, y / c); };
    HD TwoVector operator-() const { return TwoVector(-x, -y); }
    HD TwoVector& operator+=(const TwoVector& v) {
        x += v.x, y += v.y;
        return *this;
    };
    HD TwoVector& operator-=(const TwoVector& v) {
        x -= v.x, y -= v.y;
        return *this;
    };
    HD TwoVector& operator*=(const T c) {
        x *= c, y *= c;
        return *this;
    };
    HD TwoVector& operator/=(const T c) {
        x /= c, y /= c;
        return *this;
    };

    HD bool operator==(const TwoVector& v) const { return (x == v.x) and (y == v.y); };
    HD bool operator!=(const TwoVector& v) const { return (x != v.x) or (y != v.y); };

    HD real norm() const { return std::sqrt(static_cast<real>(dot(*this))); };
    HD T dot(const TwoVector& v) const { return x * v.x + y * v.y; };
    HD TwoVector unit() const { return *this / norm(); };
};

// Right operator overloads.
template <typename T> HD TwoVector<T> operator*(const T c, const TwoVector<T>& v) { return v * c; };
template <typename T> HD TwoVector<T> operator/(const T c, const TwoVector<T>& v) { return v / c; };
template <typename T> std::ostream& operator<<(std::ostream& out, const TwoVector<T>& v) {
    return out << "(" << v.x << ", " << v.y << ")";
};

/**
 * A class for a mathematical vector in \f$ \mathbb{R}^3 \f$, supporting the
 * standard operations.
 */
template <typename T> struct ThreeVector {
  public:
    template <typename U> friend std::ostream& operator<<(std::ostream&, const ThreeVector<U>&);
    T x, y, z;

    HD ThreeVector() = default;
    /**
     * Constructs a ThreeVector with provided coordinates.
     * \param x the x coordinate.
     * \param y the y coordinate.
     * \param z the z coordinate.
     */
    HD ThreeVector(const T x, const T y, const T z) : x(x), y(y), z(z){};
    /**
     * Constructs a ThreeVector with all same provided value.
     * \param x the value in all components.
     */
    HD ThreeVector(const T x) : x(x), y(x), z(x){};

    /**
     * Recasts a pointer into a ThreeVector at that pointer.
     * \param p the pointer to the data.
     */
    HD static ThreeVector& from_pointer(T* const p) {
        return *reinterpret_cast<ThreeVector<T>*>(p);
    };

    HD T& operator[](const uint i) {
        ERROR_IF(i, >, 2, "Component number must be [0:2]");
        return *(&x + i);
    };
    HD const T& operator[](const uint i) const {
        ERROR_IF(i, >, 2, "Component number must be [0:2]");
        return *(&x + i);
    };

    HD ThreeVector operator+(const ThreeVector& v) const {
        return ThreeVector(x + v.x, y + v.y, z + v.z);
    };
    HD ThreeVector operator-(const ThreeVector& v) const {
        return ThreeVector(x - v.x, y - v.y, z - v.z);
    };
    HD ThreeVector operator*(const T c) const { return ThreeVector(x * c, y * c, z * c); };
    HD ThreeVector operator/(const T c) const { return ThreeVector(x / c, y / c, z / c); };
    HD ThreeVector operator-() const { return ThreeVector(-x, -y, -z); }
    HD ThreeVector& operator+=(const ThreeVector& v) {
        x += v.x, y += v.y, z += v.z;
        return *this;
    };
    HD ThreeVector& operator-=(const ThreeVector& v) {
        x -= v.x, y -= v.y, z -= v.z;
        return *this;
    };
    HD ThreeVector& operator*=(const T c) {
        x *= c, y *= c, z *= c;
        return *this;
    };
    HD ThreeVector& operator/=(const T c) {
        x /= c, y /= c, z /= c;
        return *this;
    };

    HD bool operator==(const ThreeVector& v) const {
        return (x == v.x) and (y == v.y) and (z == v.z);
    };
    HD bool operator!=(const ThreeVector& v) const {
        return (x != v.x) and (y != v.y) and (z != v.z);
    };

    HD real norm() const { return std::sqrt(static_cast<real>(dot(*this))); };
    HD T dot(const ThreeVector& v) const { return x * v.x + y * v.y + z * v.z; };
    HD ThreeVector cross(const ThreeVector& v) const {
        return ThreeVector(y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x);
    };
};

// Right operator overloads.
template <typename T> HD ThreeVector<T> operator*(const T c, const ThreeVector<T>& v) {
    return v * c;
};

template <typename T> HD ThreeVector<T> operator/(const T c, const ThreeVector<T>& v) {
    return v / c;
};

template <typename T> std::ostream& operator<<(std::ostream& out, const ThreeVector<T>& v) {
    return out << "(" << v.x << ", " << v.y << ", " << v.z << ")";
};

}; // namespace curv3d

// Overloading hash operation.
namespace std {
template <typename T> struct hash<curv3d::TwoVector<T>> {
    size_t operator()(const curv3d::TwoVector<T>& v) const {
        // Taken from boost::hash_combine, where we just hash each one and combine
        // the hashes.
        size_t val = hash<T>()(v.x);
        val ^= hash<T>()(v.y) + 0x517cc1b727220a95 + (val << 6) + (val >> 2);
        return val;
    }
};
template <typename T> struct hash<curv3d::ThreeVector<T>> {
    size_t operator()(const curv3d::ThreeVector<T>& v) const {
        // Taken from boost::hash_combine, where we just hash each one and combine
        // the hashes.
        size_t val = hash<T>()(v.x);
        val ^= hash<T>()(v.y) + 0x517cc1b727220a95 + (val << 6) + (val >> 2);
        val ^= hash<T>()(v.z) + 0x517cc1b727220a95 + (val << 6) + (val >> 2);
        return val;
    }
};
}; // namespace std

// Type aliases.
namespace curv3d {
using SPoint = TwoVector<real>;   /**< Represents a coordinat (u,v) on a surface primitive. */
using SVector = TwoVector<real>;  /**< Represents a two vector for surface gradients. */
using Order = TwoVector<uint>;    /**< Represents the order of a surface primitive. */
using TwoIndex = TwoVector<uint>; /**< Represents a two-tuple of indices. */

using GPoint = ThreeVector<real>; /**< Represents a coordinate in global space. */
using LPoint = ThreeVector<real>; /**< Represents a coordinate in local space. */
using Vector = ThreeVector<real>; /**< Represents a three-dimensional vector. */
using Angles = ThreeVector<real>; /**< Represents three angles for a basis. */

using ControlList = std::vector<std::vector<uint>>; /**< Represents the indices of control points. */

struct BoundingBox { /**< Represents the two corners of a bounding box. */
    SPoint min = SPoint(std::numeric_limits<real>::max());
    SPoint max = SPoint(std::numeric_limits<real>::min());
};

}; // namespace curv3d

#endif // UTILS_HPP
