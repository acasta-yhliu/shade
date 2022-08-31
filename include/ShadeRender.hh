#pragma once

#include <ShadeAccelerate.hh>
#include <ShadeConfig.hh>
#include <ShadeObject.hh>

#include <fstream>
#include <iostream>
#include <map>
#include <string>

namespace Shade {
  class Scene {
    Array<Object *> objects;
    Array<Object *> lights;
    Array<Material *> materials;
    BVHNode *root;

    Vector3d lorg;

  public:
    Scene() : root(nullptr) {}

    BVHNode *buildScene() {
      if (root)
        delete root;
      root = new BVHNode(objects.data(), objects.size());
      return root;
    }

    void addObject(Object *obj) { objects.push_back(obj); }

    void addLight(Object *light) {
      objects.push_back(light);
      lights.push_back(light);
    }

    size_t objectSize() const { return objects.size(); }
    size_t lightSize() const { return lights.size(); }

    void configLight(const Vector3d &l) { lorg = l; }

    void generateRay(Ray &pr, Vector3d &f, RandEngine &gen) {
      lights[0]->emit(pr, f, gen);
      f *= (4 * PI);
    }

    ~Scene() {
      for (auto i : objects)
        delete i;
      for (auto i : materials)
        delete i;
      if (root)
        delete root;
    }
  };

  class Camera {
    int w, h;
    HashGrid grid;

    Vector3d *c; // canvas
    BVHNode *root;

    Array<RandEngine> pixelRand;
    Array<RandEngine> photonRand;

    void trace(const Ray &r, int dpt, HPoint *m, const Vector3d &fl,
               const Vector3d &adj, RandEngine &gen) {
      dpt++;
      Hitrec h;
      if (!root->intersect(r, h, 1e-4) || (dpt >= 10))
        return;
      const Material *mat = h.mat;
      double t = h.t;
      Vector3d n = h.norm;
      Vector3d x = r.o + r.d * t, f = mat->query(h.uv);
      Vector3d nl = n.dot(r.d) < 0 ? n : n * -1;
      double p = f.x() > f.y() && f.x() > f.z() ? f.x()
                 : f.y() > f.z()                ? f.y()
                                                : f.z();
#define GETBRDF(prop) PRE_BRDFS[mat->brdf].prop
      auto spec = GETBRDF(spec), diff = GETBRDF(diff), refr = GETBRDF(refr);
      auto s = spec + diff + refr;
      auto act = gen(); //* s;

      if (diff > 0 && act <= diff && act > 0) {
        double r1 = 2. * PI * gen(), r2 = gen();
        Vector3d w = nl, u = ((fabs(w.x()) > .1 ? Vector3d(0, 1, 0)
                                                : Vector3d(1, 0, 0)) %
                              w)
                                 .normalized();
        Vector3d v = w % u;

        if (m) {
          m->f = f * adj;
          m->pos = x;
          m->norm = n;
        } else {
          auto a = gen();
          if (a <= GETBRDF(rhoS)) {
            double r2s = pow(r2, 1.0 / (GETBRDF(phongS) + 1));
            auto d = (u * cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1 - r2))
                         .normalized();
            trace(Ray(x, d), dpt, m, (f * fl), adj, gen);
          } else {
            a -= GETBRDF(rhoS);
            auto &hp = grid[x];
            for (auto hitpoint : hp) {
              Vector3d v = hitpoint->pos - x;
              if ((hitpoint->norm.dot(n) > 1e-3) &&
                  (v.dot(v) <= hitpoint->r2)) {
                double g =
                    (hitpoint->n * ALPHA + ALPHA) / (hitpoint->n * ALPHA + 1.0);
                hitpoint->r2 = hitpoint->r2 * g;
                hitpoint->n++;
                hitpoint->flux =
                    (hitpoint->flux + (hitpoint->f * fl) * (1. / PI)) * g;
              }
            }
            if (a <= GETBRDF(rhoD)) {
              double r2s = sqrt(r2);
              auto d =
                  (u * cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1 - r2))
                      .normalized();
              trace(Ray(x, d), dpt, m, (f * fl) /* (1. / GETBRDF(rhoD))*/, adj,
                    gen);
            }
          }
        }
      }
      act -= diff;

      if (spec > 0 && act <= spec && act > 0) {
        trace(Ray(x, r.d - n * 2.0 * n.dot(r.d)), dpt, m, f * fl, f * adj, gen);
      }
      act -= spec;

      if (refr > 0 && act <= refr && act > 0) {
        Ray lr(x, r.d - n * 2.0 * n.dot(r.d));
        bool into = (n.dot(nl) > 0.0);
        double nc = 1.0, nt = GETBRDF(refN), nnt = into ? nc / nt : nt / nc,
               ddn = r.d.dot(nl), cos2t;

        // total internal reflection
        if ((cos2t = 1 - nnt * nnt * (1 - ddn * ddn)) < 0)
          return trace(lr, dpt, m, fl, adj, gen);

        Vector3d td =
            (r.d * nnt - n * ((into ? 1 : -1) * (ddn * nnt + sqrt(cos2t))))
                .normalized();
        double a = nt - nc, b = nt + nc, R0 = a * a / (b * b),
               c = 1 - (into ? -ddn : td.dot(n));
        double Re = R0 + (1 - R0) * c * c * c * c * c, P = Re;
        Ray rr(x, td);
        Vector3d fa = f * adj;
        if (m) {
          trace(lr, dpt, m, fl, fa * Re, gen);
          trace(rr, dpt, m, fl, fa * (1.0 - Re), gen);
        } else {
          (gen() < P) ? trace(lr, dpt, m, fl, fa, gen)
                      : trace(rr, dpt, m, fl, fa, gen);
        }
      }

      if (mat->brdf == BRDF::SCATTER) {
        trace(Ray(r.o, unitSphereRandom(gen)), dpt, m, fl, f, gen);
      }
    }

    int gamma(double x) { return int(pow(1 - exp(-x), 1 / 2.2) * 255 + .5); }

    Vector3d lowerLeft, horizon, vertical, pos, _u, _v, _w;
    double lensR;

    Ray getRay(double s, double t, RandEngine &rand) {
      auto rd = lensR * unitSphereRandom(rand);
      auto offset = _u * rd.x() + _v * rd.y();
      return Ray(
          pos + offset,
          (lowerLeft + s * horizon + t * vertical - pos - offset).normalized());
    }

    void saveBackup(int nr) {
      char name[50];
      sprintf(name, "progressive/backup%d.ppm", nr);
      FILE *f = fopen(name, "w");
      fprintf(f, "P3\n%d %d\n255\n", w, h);
      for (int i = 0; i < w * h; i++) {
        uint8_t r = gamma(c[i].x() / nr);
        uint8_t g = gamma(c[i].y() / nr);
        uint8_t b = gamma(c[i].z() / nr);
        fprintf(f, "%d %d %d ", r, g, b);
      }
      fclose(f);
    }

  public:
    Camera(int w, int h) : w(w), h(h), grid(w, h) {
      c = new Vector3d[w * h];
      pixelRand.resize(w * h);
      for (auto y = 0; y < h; ++y)
        for (auto x = 0; x < w; ++x)
          pixelRand[x + y * w].setState(rand(), rand());
    }

    void configCamera(const Vector3d &pos, const Vector3d &lookAt,
                      const Vector3d &up, double fov, double aperture = 0,
                      double focusDisk = 10.0) {
      lensR = aperture / 2;
      auto ratio = double(w) / double(h);
      this->pos = pos;
      auto theta = fov * PI / 180;
      auto halfH = tan(theta / 2);
      auto halfW = ratio * halfH;
      _w = (pos - lookAt).normalized();
      _u = up.cross(_w).normalized();
      _v = _w.cross(_u).normalized();
      lowerLeft = pos - halfW * focusDisk * _u - halfH * focusDisk * _v -
                  focusDisk * _w;
      horizon = 2 * halfW * _u * focusDisk;
      vertical = 2 * halfH * _v * focusDisk;
    }

    void renderScene(Scene &scene, int nround, int samps) {
      photonRand.resize(samps);
      for (auto i = 0; i < samps; ++i)
        photonRand[i].setState(rand(), rand());
      root = scene.buildScene();

      for (int nr = 0; nr < nround; ++nr) {
        printf("> doing round %d\n", nr);
        // trace eye rays and store measurement points
        for (int y = 0; y < h; y++) {
          fprintf(stderr, "\rHitPointPass %5.2f%%", 100.0 * y / (h - 1));
#pragma omp parallel for schedule(guided)
          for (int x = 0; x < w; x++) {
            int pixel_index = x + y * w;
            RandEngine &gen = pixelRand[pixel_index];
            auto u = double(x + gen()) / double(w);
            auto v = double(y + gen()) / double(h);
            trace(getRay(u, 1 - v, gen), 0, &grid[pixel_index],
                  Vector3d(0, 0, 0), Vector3d(1, 1, 1), gen);
          }
        }
        fprintf(stderr, "\n");
        fprintf(stderr, "> building hash grid\n");

        // build the hash table over the measurement points
        grid.buildHashGrid();
        for (int i = 0; i < samps; i++) {
          double p = 100. * (i + 1) / samps;
          fprintf(stderr, "\rPhotonPass %5.2f%%", p);
          int m = 1000 * i;
          Ray r;
          Vector3d f;
#pragma omp parallel for
          for (int j = 0; j < 1000; j++) {
            scene.generateRay(r, f, photonRand[i]);
            trace(r, 0, nullptr, f, Vector3d(1, 1, 1), photonRand[i]);
          }
        }
        fprintf(stderr, "\rfinish\n");

        // density estimation
        for (int i = 0; i < w * h; ++i) {
          c[i] =
              c[i] + clamp(Vector3d(grid[i].flux * (1.0 / (PI * grid[i].r2 * samps * 1000.0))));
        }

        if (nr % 100 == 0 && nr > 0) {
          saveBackup(nr);
        }
      }

      for (int i = 0; i < w * h; ++i) {
        c[i] /= nround;
      }
    }

    void savePPM(const char *filename) {
      FILE *f = fopen(filename, "w");
      fprintf(f, "P3\n%d %d\n255\n", w, h);
      for (int i = 0; i < w * h; i++) {
        uint8_t r = gamma(c[i].x());
        uint8_t g = gamma(c[i].y());
        uint8_t b = gamma(c[i].z());
        fprintf(f, "%d %d %d ", r, g, b);
      }
      fclose(f);
    }
  };
} // namespace Shade