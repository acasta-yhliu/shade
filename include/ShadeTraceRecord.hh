#pragma once

#include <ShadeMaterial.hh>
#include <ShadeMath.hh>
#include <vector>

namespace Shade {
  struct Ray {
    Vector3d o, d, id;

    Ray() {}

    Ray(const Vector3d &org, const Vector3d &dir)
        : o(org), d(dir.normalized()) {
      id = invert(d);
    }
  };

  struct HPoint {
    Vector3d f, pos, norm, flux;
    bool valid;
    double r2;
    int n;
  };

  struct Hitrec {
    double t;
    Vector3d norm, uv;
    Material *mat;

    void setHit(double t, const Vector3d &n, const Vector3d &uv,
                Material *mat) {
      this->uv = uv;
      this->t = t;
      this->norm = n;
      this->mat = mat;
    }
  };
} // namespace Shade