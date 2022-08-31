#pragma once

#include <ShadeMath.hh>
#include <fstream>
#include <string>
#include <vector>

namespace Shade {
  template <typename T> using Array = std::vector<T>;

  struct BRDF {
    double spec, diff, refr, rhoD, rhoS, phongS, refN;

    BRDF() {}

    constexpr BRDF(double spec, double diff, double refr, double rhoD,
                   double rhoS, double phongS, double refN)
        : spec(spec), diff(diff), refr(refr), rhoD(rhoD), rhoS(rhoS),
          phongS(phongS), refN(refN) {}

    enum DefaultBRDF {
      DIFFUSE,
      MIRROR,
      GLASS,
      LIGHT,
      MARBLE,
      FLOOR,
      WALL,
      DESK,
      STANFORD_MODEL,
      WATER,
      TEAPOT,
      TRUEGLASS,
      SCATTER,
      METAL,
    };
  };

  constexpr BRDF PRE_BRDFS[] = {
      BRDF(0, 1, 0, 0.7, 0, 0, 0),          // DIFFUSE
      BRDF(1, 0, 0, 0, 0, 0, 0),            // MIRROR
      BRDF(0, 0, 1, 0, 0, 0, 1.65),         // GLASS
      BRDF(0, 1, 0, 0, 0, 0, 0),            // LIGHT
      BRDF(0.1, 0.9, 0, 1, 0, 50, 0),       // MARBLE
      BRDF(0.05, 0.95, 0, 0.9, 0.1, 50, 0), // FLOOR
      BRDF(0, 1, 0, 1, 0, 0, 0),            // WALL
      BRDF(0, 1, 0, 1, 0, 0, 0),            // DESK
      BRDF(0, 1, 0, 0.9, 0.1, 10, 1),       // STANFORD_MODEL
      BRDF(0.05, 0, 0.95, 0, 0, 0, 1.3),          // WATER
      BRDF(0, 0, 1, 0, 0, 0, 1.5),          // TEAPOT
      BRDF(0.0, 0.1, 0.9, 0, 0.7, 0, 1.65), // TRUEGLASS
      BRDF(0, 0.2, 0, 0.7, 0, 0, 0),        // SCATTER
      BRDF(0.75, 0.25, 0, 1, 0, 100, 0),       // METAL
  };

  struct Material {
    BRDF::DefaultBRDF brdf;
    Vector3d color;
    bool reqUV = false;

    Material() {}

    Material(BRDF::DefaultBRDF brdf, const Vector3d &color)
        : brdf(brdf), color(color) {}

    virtual Vector3d query(const Vector3d &point) const { return color; }
  };

  struct Texture : public Material {
    Array<Vector3d> data;
    int w, h;

    Texture(BRDF::DefaultBRDF brdf, const char *filename) {
      reqUV = true;
      this->brdf = brdf;
      std::ifstream f(filename);
      std::string buf;
      int colorBit;
      f >> buf; // read P3
      f >> w >> h;
      printf("> Got size %d %d\n", w, h);
      f >> colorBit;
      for (auto i = 0; i < w * h; ++i) {
        double r, g, b;
        f >> r >> g >> b;
        r /= 255.0, g /= 255.0, b /= 255.0;
        data.push_back(Vector3d(r, g, b));
      }
      printf("> Texture '%s' loaded\n", filename);
      printf("> Total size: %lld\n", data.size());
    }

    Vector3d query(const Vector3d &point) const override {
      int batchId = point.z();
      auto x = int(point.x() * w);
      auto y = int(point.y() * w);

      x = x < 0 ? 0 : x >= w ? w - 1 : x;
      y = y < 0 ? 0 : y >= w ? w - 1 : y;
      auto clr = data[x + y * w + batchId * w * w];
      return clr;
    }
  };
} // namespace Shade