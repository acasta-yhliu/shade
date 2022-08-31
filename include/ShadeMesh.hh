#pragma once

#include <ShadeAccelerate.hh>
#include <ShadeObject.hh>
#include <iostream>

namespace Shade {

  class BezierSurface : public Object {
    int **binom;
    Vector3d P(double u, double v) const {
      Vector3d res(0, 0, 0);
      for (int i = 0; i <= n; ++i)
        for (int j = 0; j <= m; ++j)
          res = res + p[i][j] * B(n, i, u) * B(m, j, v);
      return res;
    }

    double B(int n, int k, double u) const {
      return binom[n][k] * pow(u, k) * pow(1 - u, n - k);
    }

    double dB(int n, int k, double u) const {
      return binom[n][k] * (k * pow(u, k - 1) * pow(1 - u, n - k) -
                            (n - k) * pow(u, k) * pow(1 - u, n - k - 1));
    }

    Vector3d F(const Vector3d &x, const Ray &ray) const {
      return ray.o + ray.d * x.x() - P(x.y(), x.z());
    }

    Vector3d d(const Vector3d &x, const Ray &ray) const {
      double jac[3][3], jac_inv[3][3];
      jac[0][0] = ray.d.x();
      jac[1][0] = ray.d.y();
      jac[2][0] = ray.d.z();
      Vector3d du(0, 0, 0), dv(0, 0, 0);
      for (int i = 0; i <= n; ++i)
        for (int j = 0; j <= m; ++j) {
          du = du - p[i][j] * dB(n, i, x.y()) * B(m, j, x.z());
          dv = dv - p[i][j] * B(n, i, x.y()) * dB(m, j, x.z());
        }
      jac[0][1] = du.x();
      jac[1][1] = du.y();
      jac[2][1] = du.z();
      jac[0][2] = dv.x();
      jac[1][2] = dv.y();
      jac[2][2] = dv.z();

      jac_inv[0][0] = jac[1][1] * jac[2][2] - jac[1][2] * jac[2][1];
      jac_inv[0][1] = jac[0][2] * jac[2][1] - jac[0][1] * jac[2][2];
      jac_inv[0][2] = jac[0][1] * jac[1][2] - jac[0][2] * jac[1][1];
      jac_inv[1][0] = jac[1][2] * jac[2][0] - jac[1][0] * jac[2][2];
      jac_inv[1][1] = jac[0][0] * jac[2][2] - jac[0][2] * jac[2][0];
      jac_inv[1][2] = jac[1][0] * jac[0][2] - jac[0][0] * jac[1][2];
      jac_inv[2][0] = jac[1][0] * jac[2][1] - jac[1][1] * jac[2][0];
      jac_inv[2][1] = jac[0][1] * jac[2][0] - jac[0][0] * jac[2][1];
      jac_inv[2][2] = jac[0][0] * jac[1][1] - jac[1][0] * jac[0][1];
      double d = jac[0][0] * (jac[1][1] * jac[2][2] - jac[1][2] * jac[2][1]) -
                 jac[1][0] * (jac[0][1] * jac[2][2] - jac[0][2] * jac[2][1]) +
                 jac[2][0] * (jac[0][1] * jac[1][2] - jac[0][2] * jac[1][1]);
      for (int i = 0; i <= 2; ++i)
        for (int j = 0; j <= 2; ++j)
          jac_inv[i][j] /= d;

      Vector3d f = F(x, ray);
      return Vector3d(
          jac_inv[0][0] * f.x() + jac_inv[0][1] * f.y() + jac_inv[0][2] * f.z(),
          jac_inv[1][0] * f.x() + jac_inv[1][1] * f.y() + jac_inv[1][2] * f.z(),
          jac_inv[2][0] * f.x() + jac_inv[2][1] * f.y() +
              jac_inv[2][2] * f.z());
    }

    int n, m;
    Vector3d **p, mMin, mMax, mCenter;

    mutable RandEngine random;

  public:
    BezierSurface(int n, int m, Vector3d *point, const Vector3d &pos,
                  const Vector3d &scale, Material *mat)
        : Object(mat), n(n), m(m) {
      random.setState(rand(), rand());
      this->p = new Vector3d *[n + 1];
      mMin = Vector3d(1e100, 1e100, 1e100);
      mMax = -mMin;
      mCenter = Vector3d(0, 0, 0);
      for (int i = 0; i <= n; ++i) {
        this->p[i] = new Vector3d[m + 1];
        for (int j = 0; j <= m; ++j) {
          this->p[i][j] = point[i + j * n] * scale + pos;
          mMin = boxMin(mMin, this->p[i][j]);
          mMax = boxMax(mMax, this->p[i][j]);
          mCenter = mCenter + p[i][j];
        }
      }
      mCenter = mCenter / ((n + 1) * (m + 1));
      int nm = max(n, m);
      binom = new int *[nm + 1];
      for (int i = 0; i <= nm; ++i)
        binom[i] = new int[nm + 1];
      binom[0][0] = 1;
      for (int i = 1; i <= n || i <= m; ++i) {
        binom[i][0] = 1;
        for (int j = 1; j <= i; ++j)
          binom[i][j] = binom[i - 1][j] + binom[i - 1][j - 1];
      }
    }

    AABB boudingBox() const override { return AABB(mMin, mMax); }

    bool intersect(const Ray &ray, Hitrec &h, double tmin) const override {
      double mint = 1e100;
      double resu, resv;
      for (int i = 1; i <= n; ++i)
        for (int j = 1; j <= m; ++j) {
          Vector3d min = Vector3d(1e100, 1e100, 1e100);
          Vector3d max = min * -1;

          for (int _i = i - 1; _i <= i; ++_i)
            for (int _j = j - 1; _j <= j; ++_j) {
              Vector3d p = this->p[_i][_j];
              min = boxMin(min, p);
              max = boxMax(max, p);
            }

          double t = -1e100;
          if (!(ray.o >= min && ray.o <= max)) { // outside
            if (fabs(ray.d.x()) > 0)
              t = std::max(t, std::min((min.x() - ray.o.x()) / ray.d.x(),
                                       (max.x() - ray.o.x()) / ray.d.x()));
            if (fabs(ray.d.y()) > 0)
              t = std::max(t, std::min((min.y() - ray.o.y()) / ray.d.y(),
                                       (max.y() - ray.o.y()) / ray.d.y()));
            if (fabs(ray.d.z()) > 0)
              t = std::max(t, std::min((min.z() - ray.o.z()) / ray.d.z(),
                                       (max.z() - ray.o.z()) / ray.d.z()));
            if (t < 0)
              continue;
            Vector3d pp = ray.o + ray.d * t;
            if (!(pp >= min && pp <= max))
              continue;
          } else
            t = 0;

          for (int __ = 0; __ < 10; ++__) {
            double u = random(), v = random();
            auto x = Vector3d(t, u, v);
            double lambda = 1;
            double last = norm2(F(x, ray));
            for (int _ = 0; _ < 10; ++_) {
              x = x - d(x, ray);
              double cost = norm2(F(x, ray));
              if (!(norm2(x) <= 1e3)) {
                x.y() = x.z() = 1e100;
                break;
              }
              if (last - cost < 1e-8)
                break;
              last = cost;
            }
            t = x.x();
            u = x.y();
            v = x.z();
            if (0 <= u && u <= 1 && 0 <= v && v <= 1 &&
                norm2(F(x, ray)) < 1e-5) {
              mint = std::min(mint, t);
              resu = u;
              resv = v;
              break;
            }
          }
        }
      if (mint < 1e100) {
        Vector3d du(0, 0, 0), dv(0, 0, 0);
        for (int i = 0; i <= n; ++i)
          for (int j = 0; j <= m; ++j) {
            du = du - p[i][j] * dB(n, i, resu) * B(m, j, resv);
            dv = dv - p[i][j] * B(n, i, resu) * dB(m, j, resv);
          }
        h.setHit(mint, du.cross(dv), Vector3d(), mat);
        return true;
      } else
        return false;
    }
  };

  class Mesh : public Object {
    Array<Object *> faces;
    BVHNode *root;
    int nfaces;

  public:
    Mesh(const char *filename, const Vector3d &pos, const Vector3d &scale,
         const Quaterniond &rot, Material *mat, bool obj = true)
        : Object(mat) {
      if (obj) {
        auto f = fopen(filename, "r");
        char leading;
        double t1, t2, t3;
        int i1, i2, i3;
        std::vector<Vector3d> vertices;
        int id = 0;
        while (fscanf(f, "%c", &leading) != EOF) {
          // printf("scanned %c", leading);
          if (leading == 'v') {
            fscanf(f, "%lf %lf %lf", &t1, &t2, &t3);
            fgetc(f);
            vertices.push_back(scale * (rot * Vector3d(t1, t2, t3)) + pos);
          } else if (leading == 'f') {
            fscanf(f, "%d %d %d", &i1, &i2, &i3);
            fgetc(f);
            faces.push_back(new Triangle(vertices[i1 - 1], vertices[i2 - 1],
                                         vertices[i3 - 1], mat, id++));
          } else
            break;
        }
        nfaces = faces.size();
        printf("> Import '%s' :: %d Faces\n", filename, nfaces);
        root = new BVHNode(faces.data(), nfaces);
      } else {
        std::ifstream f(filename);
        f >> nfaces;
        Array<Vector3d> vertices;
        for (int _ = 0; _ < nfaces; ++_) {
          vertices.clear();
          int n, m;
          double a, b, c;
          f >> n >> m;
          for (int i = 0; i <= n; ++i) {
            for (int j = 0; j <= m; ++j) {
              f >> a >> b >> c;
              vertices.push_back(scale * (rot * Vector3d(a, b, c)) + pos);
            }
          }
          faces.push_back(
              new BezierSurface(n, m, vertices.data(), pos, scale, mat));
        }
        root = new BVHNode(faces.data(), nfaces);
        printf("> Import '%s' :: %d Faces\n", filename, nfaces);
      }
    }

    bool intersect(const Ray &r, Hitrec &h, double tmin) const override {
      return root->intersect(r, h, tmin);
    }

    AABB boudingBox() const override { return root->boudingBox(); }

    virtual ~Mesh() {
      for (auto i : faces)
        delete i;
      delete root;
    }
  };

  class InterplotMesh : public Object {
    Array<Object *> faces;
    BVHNode *root;
    int nfaces;

  public:
    InterplotMesh(const char *filename, const Vector3d &pos,
                  const Vector3d &scale, const Quaterniond &rot, Material *mat)
        : Object(mat) {
      struct VerticeNorm {
        Vector3d vertice;
        Vector3d norm;
        int n;

        VerticeNorm(const Vector3d &v) : vertice(v), norm(0, 0, 0), n(0) {}

        void addNorm(const Vector3d &nm) {
          n++;
          norm += nm;
        }

        Vector3d getNorm() const { return norm / double(n); }
      };

      struct TriIndex {
        int i1, i2, i3, id;
        TriIndex(int i1, int i2, int i3, int id)
            : i1(i1), i2(i2), i3(i3), id(id) {}
      };
      auto f = fopen(filename, "r");
      char leading;
      double t1, t2, t3;
      int i1, i2, i3;
      std::vector<VerticeNorm> vertices;
      std::vector<TriIndex> tris;
      int id = 0;
      while (fscanf(f, "%c", &leading) != EOF) {
        // printf("scanned %c", leading);
        if (leading == 'v') {
          fscanf(f, "%lf %lf %lf", &t1, &t2, &t3);
          fgetc(f);
          vertices.push_back(
              VerticeNorm(scale * (rot * Vector3d(t1, t2, t3)) + pos));
        } else if (leading == 'f') {
          fscanf(f, "%d %d %d", &i1, &i2, &i3);
          fgetc(f);
          auto a = vertices[i1 - 1].vertice, b = vertices[i2 - 1].vertice,
               c = vertices[i3 - 1].vertice;
          Vector3d v0v1 = b - a;
          Vector3d v0v2 = c - a;
          auto norm = (v0v1 % v0v2).normalized();
          vertices[i1 - 1].addNorm(norm);
          vertices[i2 - 1].addNorm(norm);
          vertices[i3 - 1].addNorm(norm);
          tris.push_back(TriIndex(i1, i2, i3, id++));
        } else
          break;
      }
      for (auto t : tris) {
        auto v1 = vertices[t.i1 - 1];
        auto v2 = vertices[t.i2 - 1];
        auto v3 = vertices[t.i3 - 1];
        faces.push_back(new InterplotTriangle(
            v1.vertice, v2.vertice, v3.vertice, v1.getNorm(), v2.getNorm(),
            v3.getNorm(), mat, t.id));
      }
      nfaces = faces.size();
      printf("> Import '%s' :: %d Faces\n", filename, nfaces);
      root = new BVHNode(faces.data(), nfaces);
    }

    bool intersect(const Ray &r, Hitrec &h, double tmin) const override {
      return root->intersect(r, h, tmin);
    }

    AABB boudingBox() const override { return root->boudingBox(); }

    virtual ~InterplotMesh() {
      for (auto i : faces)
        delete i;
      delete root;
    }
  };
} // namespace Shade