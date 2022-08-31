// Rename some type and variables for better
// naming convention and other stuff

#pragma once

#include <Eigen/Eigen>
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <vector>

namespace Shade {
  // fast random generator using xorshift128
  template <typename T> class DefaultRandEngine {
    std::default_random_engine e;
    std::uniform_real_distribution<T> u;

  public:
    DefaultRandEngine() : u(0.0, 1.0) {}

    void setState(int a, int b) { e.seed(a ^ b); }

    T operator()() { return u(e); }
  };

  using RandEngine = DefaultRandEngine<double>;

  constexpr double PI = 3.14159265358979;

  using Eigen::Affine3d;
  using Eigen::AngleAxisd;
  using Eigen::Matrix3d;
  using Eigen::Matrix4d;
  using Eigen::Vector3d;
  using Eigen::Vector4d;
  using Eigen::Quaterniond;

  using std::max;
  using std::min;

  template <typename T>
  inline Eigen::Vector3<T> boxMin(const Eigen::Vector3<T> &a,
                                  const Eigen::Vector3<T> &b) {
    return Eigen::Vector3<T>(min(a.x(), b.x()), min(a.y(), b.y()),
                             min(a.z(), b.z()));
  }

  template <typename T>
  inline Eigen::Vector3<T> boxMax(const Eigen::Vector3<T> &a,
                                  const Eigen::Vector3<T> &b) {
    return Eigen::Vector3<T>(max(a.x(), b.x()), max(a.y(), b.y()),
                             max(a.z(), b.z()));
  }

  template <typename T>
  inline Eigen::Vector3<T> invert(const Eigen::Vector3<T> &vec) {
    return Eigen::Vector3<T>(1.0 / vec.x(), 1.0 / vec.y(), 1.0 / vec.z());
  }

  template <typename T>
  inline Eigen::Vector3<T> operator*(const Eigen::Vector3<T> &a,
                                     const Eigen::Vector3<T> &b) {
    return Eigen::Vector3<T>(a.x() * b.x(), a.y() * b.y(), a.z() * b.z());
  }

  template <typename T>
  inline Eigen::Vector3<T> operator+(const Eigen::Vector3<T> &a, T val) {
    return Eigen::Vector3<T>(a.x() + val, a.y() + val, a.z() + val);
  }

  template <typename T>
  inline Eigen::Vector3<T> operator-(const Eigen::Vector3<T> &a, T val) {
    return Eigen::Vector3<T>(a.x() - val, a.y() - val, a.z() - val);
  }

  template <typename T>
  inline Eigen::Vector3<T> operator%(const Eigen::Vector3<T> &a,
                                     const Eigen::Vector3<T> &b) {
    return a.cross(b);
  }

  template <typename T>
  inline bool operator<=(const Eigen::Vector3<T> &a,
                         const Eigen::Vector3<T> &b) {
    return a.x() <= b.x() && a.y() <= b.y() && a.z() <= b.z();
  }

  template <typename T>
  inline bool operator>=(const Eigen::Vector3<T> &a,
                         const Eigen::Vector3<T> &b) {
    return a.x() >= b.x() && a.y() >= b.y() && a.z() >= b.z();
  }

  template <typename T> inline double norm2(const Eigen::Vector3<T> &a) {
    return a.x() * a.x() + a.y() * a.y() + a.z() * a.z();
  }

  template <typename T>
  inline Eigen::Vector3<T> clamp(const Eigen::Vector3<T> &a) {
    auto x = a.x() > 1.0 ? 1.0 : a.x();
    auto y = a.y() > 1.0 ? 1.0 : a.y();
    auto z = a.z() > 1.0 ? 1.0 : a.z();
    return Eigen::Vector3<T>(x, y, z);
  }

  template <typename T> class RandGroup {
    template <typename U> using Array = std::vector<U>;

    Array<DefaultRandEngine<T>> randomGenerators;

  public:
    RandGroup() {}

    void resize(int n) { randomGenerators.resize(n); }

    DefaultRandEngine<T> &operator[](int i) { return randomGenerators[i]; }
  };

  inline Vector3d unitSphereRandom(RandEngine &gen) {
    Vector3d ret;
    do {
      ret = Vector3d(gen(), gen(), gen());
    } while (ret.norm() > 1);
    return ret;
  }
} // namespace Shade