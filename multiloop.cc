#include <fstream>
#include <iostream>
#include <optional>
#include <string>
#include <sstream>

#include <geometry.hh>
#include <harmonic.h>

#define ANSI_DECLARATORS
#define REAL double
#define VOID void
extern "C" {
#include <triangle.h>
}

const double EPSILON = 1.0e-5;
const size_t LEVELS = 9;
const size_t RESOLUTION = 50;

/*
File format:
- n # number of loops, then for each loop i:
  - n_i # number of curves in loop i, then for each curve j:
    - n_j # number of control points in curve j, then for each point k:
      - xk yk # assumed to be in [-100,-100]x[100x100]
- resolution
- m # number of lines
- <what to show> # e.g. 3 0 s 1 l 2 h (s-lines of side 0, lambdas of corner 1, h-lines of side 2)
 */

using namespace Geometry;

using Segment = Point2DVector;
using Curve = Point2DVector;
using Loop = std::vector<Curve>;

struct Setup {
  double resolution;
  size_t number_of_lines;
  enum class LineType { LAMBDA, S, H };
  using Line = std::pair<size_t, LineType>;
  std::vector<Line> lines;
  std::vector<Loop> loops;
};

void bernstein(size_t n, double u, DoubleVector &coeff) {
  coeff.clear(); coeff.reserve(n + 1);
  coeff.push_back(1.0);
  double u1 = 1.0 - u;
  for (size_t j = 1; j <= n; ++j) {
    double saved = 0.0;
    for (size_t k = 0; k < j; ++k) {
      double  tmp = coeff[k];
      coeff[k] = saved + tmp * u1;
      saved = tmp * u;
    }
    coeff.push_back(saved);
  }
}

Point2D bezierEval(const Point2DVector &cp, double u) {
  DoubleVector coeff;
  size_t d = cp.size() - 1;
  bernstein(d, u, coeff);
  Point2D p(0, 0);
  for (size_t j = 0; j <= d; ++j)
    p += cp[j] * coeff[j];
  return p;
}

Setup readSetup(const std::string &filename) {
  Setup result;
  std::ifstream f(filename);
  size_t n;
  f >> n;
  for (size_t i = 0; i < n; ++i) {
    Loop loop;
    size_t ni;
    f >> ni;
    for (size_t j = 0; j < ni; ++j) {
      Curve curve;
      size_t nj;
      f >> nj;
      for (size_t k = 0; k < nj; ++k) {
        double xk, yk;
        f >> xk >> yk;
        curve.emplace_back(xk, yk);
      }
      loop.push_back(curve);
    }
    result.loops.push_back(loop);
  }
  f >> result.resolution >> result.number_of_lines;
  f >> n;
  for (size_t i = 0; i < n; ++i) {
    size_t j;
    char c;
    auto t = Setup::LineType::LAMBDA;
    f >> j >> c;
    switch (c) {
    case 's': t = Setup::LineType::S; break;
    case 'h': t = Setup::LineType::H; break;
    case 'l': break;
    default:
      std::cerr << "Unknown line type: " << c << std::endl;
    }
    result.lines.emplace_back(j, t);
  }
  return result;
}

std::vector<HarmonicMap *> initializeMaps(const Setup &setup) {
  const double min[] = { -100, -100 }, max[] = { 100, 100 };
  std::vector<HarmonicMap *> result;
  for (const auto &l : setup.loops) {
    for (size_t i = 0; i < l.size(); ++i) {
      auto map = harmonic_create(min, max, LEVELS);
      for (const auto &loop : setup.loops) {
        for (size_t j = 0; j < loop.size(); ++j) {
          const auto &curve = loop[j];
          size_t n = curve.size();
          std::vector<double> points;
          for (size_t k = 0; k < n; ++k) {
            const auto &p = curve[k];
            double v = 0.0;
            if (&l == &loop) {
              if (i == j)
                v = (double)k / (n - 1);
              else if ((i + 1) % l.size() == j)
                v = (double)(n - 1 - k) / (n - 1);
            }
            points.push_back(p[0]);
            points.push_back(p[1]);
            points.push_back(v);
          }
          harmonic_add_curve(map, &points[0], n);
        }
      }
      harmonic_solve(map, EPSILON, false);
      result.push_back(map);
    }
  }
  return result;
}

Point2DVector discretizeCurves(const Loop &loop) {
  Point2DVector result;

  for (const auto &curve : loop) {
    if (curve.size() == 2)
      result.push_back(curve[0]);
    else {
      for (size_t i = 0; i < RESOLUTION; ++i) {
        double u = (double)i / RESOLUTION;
        result.push_back(bezierEval(curve, u));
      }
    }
  }

  return result;
}

TriMesh initializeMesh(const Setup &setup) {
  DoubleVector points;
  std::vector<int> segments;

  for (const auto &loop : setup.loops) {
    size_t base = points.size();
    auto polyline = discretizeCurves(loop);
    size_t n = polyline.size();
    for (const auto &p : polyline) {
      points.push_back(p[0]);
      points.push_back(p[1]);
    }
    for (size_t i = 0; i < n; ++i) {
      segments.push_back(base + i);
      segments.push_back(base + i + 1);
    }
    segments.back() = base;
  }

  // Setup output data structure
  struct triangulateio in, out;
  in.pointlist = &points[0];
  in.numberofpoints = points.size() / 2;
  in.numberofpointattributes = 0;
  in.pointmarkerlist = nullptr;
  in.segmentlist = &segments[0];
  in.numberofsegments = segments.size() / 2;
  in.segmentmarkerlist = nullptr;
  in.numberofholes = 0;
  in.numberofregions = 0;

  // Setup output data structure
  out.pointlist = nullptr;
  out.pointattributelist = nullptr;
  out.pointmarkerlist = nullptr;
  out.trianglelist = nullptr;
  out.triangleattributelist = nullptr;
  out.segmentlist = nullptr;
  out.segmentmarkerlist = nullptr;

  // Call the library function [with maximum triangle area = resolution]
  std::ostringstream cmd;
  cmd << "pqa" << std::fixed << setup.resolution << "DBPzQ";
  triangulate(const_cast<char *>(cmd.str().c_str()), &in, &out, (struct triangulateio *)nullptr);

  // Process the result
  TriMesh mesh;
  mesh.resizePoints(out.numberofpoints);
  for (int i = 0; i < out.numberofpoints; ++i) {
    mesh[i][0] = out.pointlist[2*i];
    mesh[i][1] = out.pointlist[2*i+1];
    mesh[i][2] = 0.0;
  }
  for (int i = 0; i < out.numberoftriangles; ++i)
    mesh.addTriangle(out.trianglelist[3*i+0],
                     out.trianglelist[3*i+1],
                     out.trianglelist[3*i+2]);

  return mesh;
}

size_t previousLine(const Setup &setup, size_t i) {
  size_t n = 0;
  for (const auto &l : setup.loops) {
    if (i - n < l.size()) {
      if (i == n)
        return n + l.size() - 1;
      return i - 1;
    }
    n += l.size();
  }
  return -1; // should not come here
}

std::vector<Segment> generateSegments(const Setup &setup,
                                      const std::vector<HarmonicMap *> &maps,
                                      const TriMesh &mesh) {
  std::vector<Segment> result;
  double density = 1.0 / setup.number_of_lines;
  for (const auto &line : setup.lines) {
    size_t i = line.first, im = previousLine(setup, i);
    auto eval = [&](const Point3D &p) {
      double point[] = { p[0], p[1] }, value, value2, denom, result;
      switch (line.second) {
      case Setup::LineType::LAMBDA:
        harmonic_eval(maps[i], point, &value);
        result = value;
        break;
      case Setup::LineType::S:
        harmonic_eval(maps[i], point, &value);
        harmonic_eval(maps[im], point, &value2);
        denom = value + value2;
        result = denom < EPSILON ? 0.0 : value / denom;
        break;
      case Setup::LineType::H:
        harmonic_eval(maps[i], point, &value);
        harmonic_eval(maps[im], point, &value2);
        result = 1.0 - value - value2;
        break;
      default:;
      }
      return result;
    };
    auto slice = [&](size_t i, size_t j) {
      double x = eval(mesh[i]), y = eval(mesh[j]);
      if (y > x) {
        std::swap(x, y);
        std::swap(i, j);
      }
      int q1 = x / density;
      int q2 = y / density;
      if (q1 - q2 == 1) {
        double alpha = (density * q1 - y) / (x - y);
        auto p = mesh[i] * alpha + mesh[j] * (1 - alpha);
        return std::optional<Point2D>(std::in_place, p[0], p[1]);
      } else
        return std::optional<Point2D>();
    };
    for (const auto &tri : mesh.triangles()) {
      auto ab = slice(tri[0], tri[1]);
      auto ac = slice(tri[0], tri[2]);
      auto bc = slice(tri[1], tri[2]);
      if (ab.has_value() + ac.has_value() + bc.has_value() == 2) {
        Segment segment;
        if (ab)
          segment.push_back(ab.value());
        if (ac)
          segment.push_back(ac.value());
        if (bc)
          segment.push_back(bc.value());
        result.push_back(segment);
      }
    }
  }
  return result;
}

void writeSegments(const Setup &setup, const std::vector<Segment> &segments,
                   const std::string &filename) {
  auto scale = [](const Point2D &p) { return p + Vector2D(100, 100); };
  std::ofstream f(filename);
  f << "%!PS-Adobe-2.0\n%%BoundingBox: 0 0 200 200\n%%EndComments\n";
  f << "1 setlinewidth" << std::endl;
  for (const auto &loop : setup.loops) {
    auto polyline = discretizeCurves(loop);
    f << "newpath" << std::endl;
    auto p = scale(polyline.front());
    f << p[0] << ' ' << p[1] << " moveto" << std::endl;
    for (size_t i = 1; i < polyline.size(); ++i) {
      auto p = scale(polyline[i]);
      f << p[0] << ' ' << p[1] << " lineto" << std::endl;
    }
    f << "closepath stroke" << std::endl;
  }
  f << "0.5 setlinewidth" << std::endl;
  for (const auto &s : segments) {
    auto p1 = scale(s[0]), p2 = scale(s[1]);
    f << "newpath " << p1[0] << ' ' << p1[1] << " moveto" << std::endl;
    f << p2[0] << ' ' << p2[1] << " lineto stroke" << std::endl;
  }
  f << "showpage" << std::endl;
}

int main(int argc, char **argv) {
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " <input.mlp>" << std::endl;
    return 1;
  }

  auto setup = readSetup(argv[1]);
  auto maps = initializeMaps(setup);
  auto mesh = initializeMesh(setup);

  auto segments = generateSegments(setup, maps, mesh);
  writeSegments(setup, segments, "output.eps");

  for (auto m : maps)
    harmonic_free(m);
}
