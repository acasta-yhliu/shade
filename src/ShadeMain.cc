#include <Shade>
#include <fstream>
#include <string>
using namespace Shade;

struct ShadeConfig {
  int w, h, samps, nr;
  std::string filename;
  Vector3d pos, lookAt, up;
  double fov, aperture, focus;
};

BRDF::DefaultBRDF getMaterial(const std::string &str) {
#define RETMAT(X)                                                              \
  if (str == #X)                                                               \
    return BRDF::X;
  RETMAT(DIFFUSE)
  RETMAT(MIRROR)
  RETMAT(GLASS)
  RETMAT(LIGHT)
  RETMAT(MARBLE)
  RETMAT(FLOOR)
  RETMAT(WALL)
  RETMAT(DESK)
  RETMAT(STANFORD_MODEL)
  RETMAT(WATER)
  RETMAT(TEAPOT)
  RETMAT(TRUEGLASS)
  RETMAT(SCATTER)
  RETMAT(METAL)
  return BRDF::DIFFUSE;
}

Object *getObject(std::ifstream &f, Material *&mat,
                  std::map<std::string, Material *> &matMap) {
  std::string objName, matName;
  f >> objName;
  if (objName == "<<SPHERE>>") {
    double rad, x, y, z;
    f >> rad >> x >> y >> z >> matName;
    mat = matMap[matName];
    return new Sphere(rad, Vector3d(x, y, z), mat);
  } else if (objName == "<<POLY>>") {
    std::string filename;
    double x, y, z, s1, s2, s3, r1, r2, r3;
    f >> filename >> x >> y >> z >> s1 >> s2 >> s3 >> r1 >> r2 >> r3 >> matName;
    auto rot = AngleAxisd(r3, Vector3d::UnitZ()) *
               AngleAxisd(r2, Vector3d::UnitY()) *
               AngleAxisd(r1, Vector3d::UnitX());
    mat = matMap[matName];
    return new Mesh(filename.c_str(), Vector3d(x, y, z), Vector3d(s1, s2, s3),
                    rot, mat);
  } else if (objName == "<<INTERPOLY>>") {
    std::string filename;
    double x, y, z, s1, s2, s3, r1, r2, r3;
    f >> filename >> x >> y >> z >> s1 >> s2 >> s3 >> r1 >> r2 >> r3 >> matName;
    auto rot = AngleAxisd(r3, Vector3d::UnitZ()) *
               AngleAxisd(r2, Vector3d::UnitY()) *
               AngleAxisd(r1, Vector3d::UnitX());
    mat = matMap[matName];
    return new InterplotMesh(filename.c_str(), Vector3d(x, y, z),
                             Vector3d(s1, s2, s3), rot, mat);
  } else if (objName == "<<BEZIER>>") {
    std::string filename;
    double x, y, z, s1, s2, s3, r1, r2, r3;
    f >> filename >> x >> y >> z >> s1 >> s2 >> s3 >> r1 >> r2 >> r3 >> matName;
    auto rot = AngleAxisd(r3, Vector3d::UnitZ()) *
               AngleAxisd(r2, Vector3d::UnitY()) *
               AngleAxisd(r1, Vector3d::UnitX());
    mat = matMap[matName];
    return new Mesh(filename.c_str(), Vector3d(x, y, z), Vector3d(s1, s2, s3),
                    rot, mat, false);
  } else if (objName == "<<TRI>>") {
    Vector3d v0, v1, v2;
    f >> v0.x() >> v0.y() >> v0.z();
    f >> v1.x() >> v1.y() >> v1.z();
    f >> v2.x() >> v2.y() >> v2.z();
    f >> matName;
    mat = matMap[matName];
    return new Triangle(v0, v1, v2, mat);
  } else if (objName == "<<TRANS>>") {
    double t1, t2, t3, s1, s2, s3, r1, r2, r3;
    f >> t1 >> t2 >> t3 >> s1 >> s2 >> s3 >> r1 >> r2 >> r3;
    Material *baseMat;
    auto baseObj = getObject(f, baseMat, matMap);
    mat = baseMat;
    Affine3d affine;
    affine.setIdentity();
    affine.translate(Vector3d(t1, t2, t3));
    affine.scale(Vector3d(s1, s2, s3));
    affine.rotate(AngleAxisd(r3, Vector3d::UnitZ()) *
                  AngleAxisd(r2, Vector3d::UnitY()) *
                  AngleAxisd(r1, Vector3d::UnitX()));
    std::cout << affine.matrix() << std::endl
              << affine.inverse().matrix() << std::endl;
    auto obj =
        new Transform(affine.matrix(), affine.inverse().matrix(), baseObj);
    return obj;
  } else if (objName == "<<VOLUMN>>") {
    double d;
    f >> d >> matName; // it actually override the base object
    auto baseObj = getObject(f, mat, matMap);
    mat = matMap[matName];
    auto obj = new ConstantMedium(baseObj, d, mat);
    return obj;
  } else {
    fprintf(stdout, "> Unknown object '%s'\n", objName.c_str());
    exit(EXIT_FAILURE);
  }
  return nullptr;
}

void parseScene(const char *filename, Scene &scene, ShadeConfig &config) {
  std::ifstream f(filename);
  std::string configWord;
  std::map<std::string, Material *> matMap;
  while (true) {
    f >> configWord;
    if (configWord == "<<BEGIN>>") {
      while (true) {
        f >> configWord;
        if (configWord == "<<END>>")
          break;
        std::cout << "> Loading :: " << configWord << std::endl;
        if (configWord == "<<CONFIG>>") {
          // folled by size, nround, nsamp, name
          f >> configWord;
          f >> config.w;
          f >> config.h;
          f >> configWord;
          f >> config.nr;
          f >> configWord;
          f >> config.samps;
          f >> configWord;
          f >> config.filename;
        } else if (configWord == "<<MATERIAL>>") {
          int nmat;
          std::string matName;
          std::string matType;
          std::string texType;
          double c1, c2, c3;
          f >> nmat;
          for (auto i = 0; i < nmat; ++i) {
            f >> matName;
            f >> matType;
            f >> texType;
            BRDF::DefaultBRDF mat = getMaterial(matType);
            if (texType == "<<CONST>>") {
              f >> c1 >> c2 >> c3;
              matMap[matName] = new Material(mat, Vector3d(c1, c2, c3));
              printf(">>> Material '%s' :: '%s'(%d) (%lf, %lf, %lf)\n",
                     matName.c_str(), matType.c_str(), mat, c1, c2, c3);
            } else if (texType == "<<MAP>>") {
              std::string fname;
              f >> fname;
              matMap[matName] = new Texture(mat, fname.c_str());
              printf(">>> Material '%s' :: '%s'\n", matName.c_str(),
                     fname.c_str());
            }
          }
        } else if (configWord == "<<SCENE>>") {
          int nobj;
          f >> nobj;
          Material *mat = nullptr;
          for (int i = 0; i < nobj; ++i) {
            auto obj = getObject(f, mat, matMap);
            if (mat->brdf == BRDF::LIGHT)
              scene.addLight(obj);
            else
              scene.addObject(obj);
          }
        } else if (configWord == "<<CAMERA>>") {
          f >> configWord;
          f >> config.pos.x() >> config.pos.y() >> config.pos.z();
          f >> configWord;
          f >> config.lookAt.x() >> config.lookAt.y() >> config.lookAt.z();
          f >> configWord;
          f >> config.up.x() >> config.up.y() >> config.up.z();
          f >> configWord;
          f >> config.fov >> config.aperture >> config.focus;
        }
      }
    } else if (configWord == "<<ENDCONFIG>>")
      break;
  }
}

int main(int argc, char *argv[]) {
  if (argc != 2) {
    fprintf(stderr, "Usage: Shade {Configuration}\n");
    exit(EXIT_FAILURE);
  }

  const char *configFile = argv[1];
  ShadeConfig config;
  Scene scene;
  parseScene(configFile, scene, config);

  printf("> Finish Parsing, Display Details\n");
  printf("> Canvas :: %dx%d\n", config.w, config.h);
  printf("> Render :: %d rounds, %d samps\n", config.nr, config.samps * 1000);

  Camera cam(config.w, config.h);
  cam.configCamera(config.pos, config.lookAt, config.up, config.fov,
                   config.aperture, config.focus);
  cam.renderScene(scene, config.nr, config.samps);
  cam.savePPM(config.filename.c_str());
}