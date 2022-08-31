#pragma once

#include <ShadeObject.hh>
#include <ShadeTraceRecord.hh>
#include <vector>

namespace Shade {
  template <typename T> using Array = std::vector<T>;

  class HashGrid {
    Array<HPoint> hitpoints;
    Array<Array<HPoint *>> hashGrid;
    AABB hpbbox;

    uint32_t numHash, numPhoton;
    double hashS;

    int w, h;

    uint32_t hash(int ix, int iy, int iz) {
      return (uint32_t)((ix * 73856093) ^ (iy * 19349663) ^ (iz * 83492791)) %
             numHash;
    }

  public:
    HashGrid() {}

    void resize(int w, int h) {
      this->w = w, this->h = h;
      hitpoints.resize(w * h);
      hashGrid.resize(w * h);
    }

    HashGrid(int w, int h) : w(w), h(h) {
      hitpoints.resize(w * h);
      hashGrid.resize(w * h);
    }

    void buildHashGrid() {
      hpbbox.reset();
      for (auto &hp : hitpoints)
        hpbbox.fit(hp.pos);
      Vector3d ssize = hpbbox.max - hpbbox.min;
      double irad =
          ((ssize.x() + ssize.y() + ssize.z()) / 3.0) / ((w + h) / 2.0) * 2.0;
      hpbbox.reset();
      int vphoton = 0;
      for (auto &hp : hitpoints) {
        hp.r2 = irad * irad;
        hp.n = 0;
        hp.flux = Vector3d();
        vphoton++;
        hpbbox.fit(hp.pos - irad);
        hpbbox.fit(hp.pos + irad);
      }

      // make each grid cell two times larger than the initial radius
      hashS = 1.0 / (irad * 2.0);
      numHash = vphoton;

      // build the hash table
      for (uint32_t i = 0; i < numHash; i++)
        hashGrid[i].clear();
      for (auto &hp : hitpoints) {
        Vector3d BMin = ((hp.pos - irad) - hpbbox.min) * hashS;
        Vector3d BMax = ((hp.pos + irad) - hpbbox.min) * hashS;
        for (int iz = abs(int(BMin.z())); iz <= abs(int(BMax.z())); iz++) {
          for (int iy = abs(int(BMin.y())); iy <= abs(int(BMax.y())); iy++) {
            for (int ix = abs(int(BMin.x())); ix <= abs(int(BMax.x())); ix++) {
              int hv = hash(ix, iy, iz);
              hashGrid[hv].push_back(&hp);
            }
          }
        }
      }
    }

    Array<HPoint *> &operator[](const Vector3d &x) {
      auto hh = (x - hpbbox.min) * hashS;
      int ix = abs(int(hh.x())), iy = abs(int(hh.y())), iz = abs(int(hh.z()));
      return hashGrid[hash(ix, iy, iz)];
    }

    HPoint &operator[](int i) { return hitpoints[i]; }
  };

  class BVHNode : public Object {
    Object *left, *right;
    AABB box;

    static bool sortX(const Object *a, const Object *b) {
      return a->boudingBox().min.x() < b->boudingBox().min.x();
    }

    static bool sortY(const Object *a, const Object *b) {
      return a->boudingBox().min.y() < b->boudingBox().min.y();
    }

    static bool sortZ(const Object *a, const Object *b) {
      return a->boudingBox().min.z() < b->boudingBox().min.z();
    }

  public:
    BVHNode(Object **objs, int n, int axis = 0) : Object(nullptr) {
      if (n == 1) {
        left = right = objs[0];
        box = objs[0]->boudingBox();
        return;
      } else {
        switch (axis) {
        case 0:
          std::sort(objs, objs + n, sortX);
          break;
        case 1:
          std::sort(objs, objs + n, sortY);
          break;
        case 2:
          std::sort(objs, objs + n, sortZ);
          break;
        }
        left = new BVHNode(objs, n / 2, (axis + 1) % 3);
        right = new BVHNode(objs + n / 2, n - n / 2, (axis + 1) % 3);
        box = left->boudingBox() + right->boudingBox();
      }
    }

    bool intersect(const Ray &ray, Hitrec &h,
                   double tmin = 1e-4) const override {
      if (box.intersect(ray, tmin)) {
        if (left == right)
          return left->intersect(ray, h, tmin);
        Hitrec hLeft, hRight;
        bool isLeft = left->intersect(ray, hLeft, tmin);
        bool isRight = right->intersect(ray, hRight, tmin);
        if (isLeft && isRight)
          h = hLeft.t < hRight.t ? hLeft : hRight;
        else if (isLeft)
          h = hLeft;
        else if (isRight)
          h = hRight;
        return isLeft || isRight;
      }
      return false;
    }

    AABB boudingBox() const override { return box; }

    virtual ~BVHNode() {
      if (left != right) {
        delete left;
        delete right;
      }
    }
  };
} // namespace Shade