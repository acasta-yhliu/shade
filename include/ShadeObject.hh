#pragma once

#include <ShadeMaterial.hh>
#include <ShadeMath.hh>
#include <ShadeTraceRecord.hh>

namespace Shade {
  struct AABB {
    Vector3d min, max; // axis aligned bounding box

    AABB() {}

    AABB(const Vector3d &min, const Vector3d &max) : min(min), max(max) {}

    void fit(const Vector3d &p) {
      min = boxMin(min, p);
      max = boxMax(max, p);
    }

    void reset() {
      min = Vector3d(1e20, 1e20, 1e20);
      max = Vector3d(-1e20, -1e20, -1e20);
    }

    AABB operator+(const AABB &b) {
      AABB ret;
      ret.min = boxMin(min, b.min);
      ret.max = boxMax(max, b.max);
      return ret;
    }

    bool intersect(const Ray &r, double tMin) const {
      Vector3d dirfrac = r.id;
      double t1 = (min.x() - r.o.x()) * dirfrac.x();
      double t2 = (max.x() - r.o.x()) * dirfrac.x();
      double t3 = (min.y() - r.o.y()) * dirfrac.y();
      double t4 = (max.y() - r.o.y()) * dirfrac.y();
      double t5 = (min.z() - r.o.z()) * dirfrac.z();
      double t6 = (max.z() - r.o.z()) * dirfrac.z();
      double tmin = std::max(std::max(std::min(t1, t2), std::min(t3, t4)),
                             std::min(t5, t6));
      double tmax = std::min(std::min(std::max(t1, t2), std::max(t3, t4)),
                             std::max(t5, t6));
      if (tmax < 0)
        return false;
      if (tmin > tmax)
        return false;
      return true;
    }
  };

  class Object {
    friend class Light;

  protected:
    Material *mat;

  public:
    Object(Material *mat) : mat(mat) {}

    virtual bool intersect(const Ray &ray, Hitrec &h, double tmin) const = 0;

    virtual AABB boudingBox() const = 0;

    virtual void emit(Ray &ray, Vector3d &f, RandEngine &rand) const {}

    virtual ~Object() {}
  };

  class Sphere : public Object {
    double rad;
    Vector3d p;

    Vector3d getUV(const Vector3d &pos) const {
      double u = 1 - (atan2(pos.z(), pos.x()) + PI) / (2 * PI);
      double v = pos.y() * 0.5 + 0.5;
      return Vector3d(u, v, 0);
    }

  public:
    Sphere(double rad, const Vector3d &p, Material *mat)
        : Object(mat), rad(rad), p(p) {}

    bool intersect(const Ray &r, Hitrec &h, double tmin = 1e-4) const override {
      Vector3d op = p - r.o;
      double t, b = op.dot(r.d), det = b * b - op.dot(op) + rad * rad;
      if (det < 0) {
        return false;
      } else {
        det = sqrt(det);
      }
      if ((t = b - det) > tmin) {
        auto pos = r.o + r.d * t;
        auto norm = (pos - p).normalized();
        auto uv = mat->reqUV ? getUV(norm) : Vector3d();
        h.setHit(t, norm.dot(r.d) < 0 ? norm : -norm, uv, mat);
        return true;
      } else if ((t = b + det) > tmin) {
        auto pos = r.o + r.d * t;
        auto norm = (pos - p).normalized();
        auto uv = mat->reqUV ? getUV(norm) : Vector3d();
        h.setHit(t, norm.dot(r.d) < 0 ? norm : -norm, uv, mat);
        return true;
      }
      return false;
    }

    AABB boudingBox() const override { return AABB(p - rad, p + rad); }

    void emit(Ray &ray, Vector3d &f, RandEngine &rand) const override {
      // random choice a point on the surface
      double pp = 2 * PI * rand(), t = 2.0 * acos(sqrt(1.0 - rand()));
      double st = sin(t);
      auto d = Vector3d(cos(pp) * st, cos(t), sin(pp) * st);
      f = mat->query(d * rad);
      ray = Ray(p + d * (rad + 0.1), d + Vector3d(0, -rad / 5, 0));
    }
  };

  class Triangle : public Object {
    Vector3d vertices[3];
    Vector3d norm;
    int batchId = 0;

    Vector3d getUV(const Vector3d &pos) const {
      auto toP = (pos - vertices[1]);
      auto l = toP.norm();
      auto x = (vertices[0] - vertices[1]);
      auto y = (vertices[2] - vertices[1]);
      auto xl = x.norm(), yl = y.norm();
      auto u = toP.dot(x) / (xl * xl);
      auto v = toP.dot(y) / (yl * yl);
      return Vector3d(u, v, batchId);
    }

  public:
    Triangle(const Vector3d &a, const Vector3d &b, const Vector3d &c,
             Material *mat, int batchId = 0)
        : Object(mat), vertices{a, b, c}, batchId(batchId) {
      Vector3d v0v1 = b - a;
      Vector3d v0v2 = c - a;
      norm = (v0v1 % v0v2).normalized();
    }

    bool intersect(const Ray &r, Hitrec &h, double tmin) const override {
      auto orig = r.o, dir = r.d;
      auto v0 = vertices[0], v1 = vertices[1], v2 = vertices[2];
      Vector3d e1 = v1 - v0, e2 = v2 - v0;
      auto p = dir % e2;
      auto det = e1.dot(p), abs_det = fabs(det);
      Vector3d t = det > 0 ? orig - v0 : v0 - orig;
      if (abs_det < 1e-6)
        return false;
      auto udet = t.dot(p);
      if (udet < 0.0f || udet > abs_det)
        return false;
      auto q = t % e1;
      auto vdet = dir.dot(q);
      if (vdet < 0.0f || udet + vdet > abs_det)
        return false;
      auto time = fabs(e2.dot(q));
      time /= abs_det;
      if (time < tmin)
        return false;
      auto uv = mat->reqUV ? getUV(r.o + r.d * time) : Vector3d();
      h.setHit(time, norm, uv, mat);
      return true;
    }

    AABB boudingBox() const override {
      auto a = vertices[0], b = vertices[1], c = vertices[2];
      return AABB(boxMin(boxMin(a, b), c), boxMax(boxMax(a, b), c));
    }

    void emit(Ray &ray, Vector3d &f, RandEngine &rand) const override {
      auto a = rand(), b = rand() * (1 - a), c = 1 - a - b;
      auto p = vertices[0] * a + vertices[1] * b + vertices[2] * c;
      f = mat->query(p);
      ray = Ray(p, norm);
    }
  };

  class InterplotTriangle : public Object {
    Vector3d vertices[3];
    Vector3d norm[3];
    int batchId = 0;

    Vector3d getUV(const Vector3d &pos) const {
      auto toP = (pos - vertices[1]);
      auto l = toP.norm();
      auto x = (vertices[0] - vertices[1]);
      auto y = (vertices[2] - vertices[1]);
      auto xl = x.norm(), yl = y.norm();
      auto u = toP.dot(x) / (xl * xl);
      auto v = toP.dot(y) / (yl * yl);
      return Vector3d(u, v, batchId);
    }

    Vector3d interplotNorm(const Vector3d &p) const {
      auto a = vertices[0], b = vertices[1], c = vertices[2];
      auto alpha = (-(p.x() - b.x()) * (c.y() - b.y()) +
                    (p.y() - b.y()) * (c.x() - b.x())) /
                   (-(a.x() - b.x()) * (c.y() - b.y()) +
                    (a.y() - b.y()) * (c.x() - b.x()));
      auto beta = (-(p.x() - c.x()) * (a.y() - c.y()) +
                   (p.y() - c.y()) * (a.x() - c.x())) /
                  (-(b.x() - c.x()) * (a.y() - c.y()) +
                   (b.y() - c.y()) * (a.x() - c.x()));
      auto gamma = 1 - alpha - beta;
      return alpha * norm[0] + beta * norm[1] + gamma * norm[2];
    }

  public:
    InterplotTriangle(const Vector3d &a, const Vector3d &b, const Vector3d &c,
                      const Vector3d &n1, const Vector3d &n2,
                      const Vector3d &n3, Material *mat, int batchId = 0)
        : Object(mat), vertices{a, b, c}, norm{n1, n2, n3}, batchId(batchId) {}

    bool intersect(const Ray &r, Hitrec &h, double tmin) const override {
      auto orig = r.o, dir = r.d;
      auto v0 = vertices[0], v1 = vertices[1], v2 = vertices[2];
      Vector3d e1 = v1 - v0, e2 = v2 - v0;
      auto p = dir % e2;
      auto det = e1.dot(p), abs_det = fabs(det);
      Vector3d t = det > 0 ? orig - v0 : v0 - orig;
      if (abs_det < 1e-6)
        return false;
      auto udet = t.dot(p);
      if (udet < 0.0f || udet > abs_det)
        return false;
      auto q = t % e1;
      auto vdet = dir.dot(q);
      if (vdet < 0.0f || udet + vdet > abs_det)
        return false;
      auto time = fabs(e2.dot(q));
      time /= abs_det;
      if (time < tmin)
        return false;
      auto pos = r.o + r.d * time;
      auto uv = mat->reqUV ? getUV(pos) : Vector3d();
      h.setHit(time, interplotNorm(pos), uv, mat);
      return true;
    }

    AABB boudingBox() const override {
      auto a = vertices[0], b = vertices[1], c = vertices[2];
      return AABB(boxMin(boxMin(a, b), c), boxMax(boxMax(a, b), c));
    }

    void emit(Ray &ray, Vector3d &f, RandEngine &rand) const override {
      auto a = rand(), b = rand() * (1 - a), c = 1 - a - b;
      auto p = vertices[0] * a + vertices[1] * b + vertices[2] * c;
      f = mat->query(p);
      ray = Ray(p, interplotNorm(p));
    }
  };

  class Transform : public Object {
    Object *baseObject;
    Matrix4d transObj;
    Matrix4d trans;

    static Vector3d transPoint(const Matrix4d &mat, const Vector3d &p) {
      auto transed = (mat * Vector4d(p.x(), p.y(), p.z(), 1));
      return Vector3d(transed.x(), transed.y(), transed.z());
    }

    static Vector3d transDir(const Matrix4d &mat, const Vector3d &d) {
      auto transed = (mat * Vector4d(d.x(), d.y(), d.z(), 0));
      return Vector3d(transed.x(), transed.y(), transed.z());
    }

  public:
    Transform(const Matrix4d &trans, const Matrix4d &invTrans, Object *baseObj)
        : Object(nullptr), transObj(trans), trans(invTrans),
          baseObject(baseObj) {}

    bool intersect(const Ray &ray, Hitrec &h, double tmin) const override {
      auto trOrg = transPoint(trans, ray.o);
      auto trDir = transDir(trans, ray.d).normalized();
      bool inter = baseObject->intersect(Ray(trOrg, trDir), h, tmin);
      if (inter)
        h.setHit(h.t, transDir(trans.transpose(), h.norm).normalized(), h.uv,
                 h.mat);
      return inter;
    }

    AABB boudingBox() const override {
      auto bbox = baseObject->boudingBox();
      return AABB(transPoint(transObj, bbox.min),
                  transPoint(transObj, bbox.max));
    }
  };

  class ConstantMedium : public Object {
    double density;
    Object *baseObject;
    mutable RandEngine engine;

  public:
    ConstantMedium(Object *bObject, double d, Material *mat)
        : Object(mat), density(d), baseObject(bObject) {
      engine.setState(rand(), rand());
    }

    AABB boudingBox() const override { return baseObject->boudingBox(); }

    bool intersect(const Ray &ray, Hitrec &h, double tmin) const override {
      Hitrec rec1, rec2;
      if (baseObject->intersect(ray, rec1, tmin)) {
        if (baseObject->intersect(ray, rec2, rec1.t + 0.0001)) {
          double t1 = max(rec1.t, tmin);
          double t2 = rec2.t;
          if (t1 >= t2)
            return false;
          t1 = max(t1, 0.0);

          double distInside = (t2 - t1) * ray.d.norm();
          double pass = exp(-distInside * density);
          if (engine() > pass) {
            h.setHit(t1, rec1.norm, Vector3d(), mat);
            return true;
          } else {
            return false;
          }
        }
      }
      return false;
    }
  };
} // namespace Shade