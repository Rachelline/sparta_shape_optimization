/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   svg_shape: implementation. See svg_shape.h.
------------------------------------------------------------------------- */

#include "svg_shape.h"
#include "bezier_geom.h"

#include <cctype>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

namespace {

// Outline of the symmetric body as a flat [x0,y0,x1,y1,...] point list.
std::vector<double> outline(const double alpha[4], double chord, int nseg)
{
  int npt = 2 * nseg + 1;              // total endpoints for the closed body
  std::vector<double> pts(2 * npt);
  std::vector<double> norms(2 * (2 * nseg));
  BezierGeom::symmetric_body_to_lines(alpha, chord, nseg, pts.data(),
                                      norms.data());
  return pts;
}

std::string capitalize(std::string s)
{
  if (!s.empty()) s[0] = (char) std::toupper((unsigned char) s[0]);
  return s;
}

}  // namespace

int write_shapes_svg(const std::string &path,
                     const double alpha_init[4], double value_init,
                     const double alpha_final[4], double value_final,
                     double chord, int nseg,
                     const std::string &objective_name)
{
  std::vector<double> o0 = outline(alpha_init, chord, nseg);
  std::vector<double> o1 = outline(alpha_final, chord, nseg);

  // bounding box over both outlines + the interior control points
  double xmin = 0.0, xmax = chord, ymin = 0.0, ymax = 0.0;
  auto grow = [&](double x, double y) {
    if (x < xmin) xmin = x; if (x > xmax) xmax = x;
    if (y < ymin) ymin = y; if (y > ymax) ymax = y;
  };
  for (size_t i = 0; i + 1 < o0.size(); i += 2) grow(o0[i], o0[i + 1]);
  for (size_t i = 0; i + 1 < o1.size(); i += 2) grow(o1[i], o1[i + 1]);
  grow(alpha_init[0], alpha_init[1]);  grow(alpha_init[2], alpha_init[3]);
  grow(alpha_final[0], alpha_final[1]); grow(alpha_final[2], alpha_final[3]);

  double wspan = xmax - xmin, hspan = ymax - ymin;
  if (wspan <= 0) wspan = 1.0;
  if (hspan <= 0) hspan = 1.0;
  double margin = 0.08 * (wspan > hspan ? wspan : hspan);
  xmin -= margin; xmax += margin; ymin -= margin; ymax += margin;
  wspan = xmax - xmin; hspan = ymax - ymin;

  const double W = 900.0, plotH = 520.0, headerH = 70.0;
  const double H = plotH + headerH;
  double s = (W / wspan < plotH / hspan) ? (W / wspan) : (plotH / hspan);
  double ox = 0.5 * (W - s * wspan);          // center horizontally
  double oy = 0.5 * (plotH - s * hspan);      // center vertically in plot area

  auto sx = [&](double x) { return ox + (x - xmin) * s; };
  auto sy = [&](double y) { return headerH + oy + (ymax - y) * s; };  // flip y

  auto polypoints = [&](const std::vector<double> &o) {
    std::ostringstream ss;
    for (size_t i = 0; i + 1 < o.size(); i += 2)
      ss << sx(o[i]) << "," << sy(o[i + 1]) << " ";
    return ss.str();
  };

  std::ofstream f(path.c_str());
  if (!f) return 1;

  f << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" << W
    << "\" height=\"" << H << "\" viewBox=\"0 0 " << W << " " << H
    << "\" font-family=\"-apple-system,Helvetica,Arial,sans-serif\">\n";
  f << "  <rect width=\"" << W << "\" height=\"" << H
    << "\" fill=\"#0f1117\"/>\n";

  // axis line at y = 0 (symmetry axis)
  f << "  <line x1=\"" << sx(xmin) << "\" y1=\"" << sy(0.0) << "\" x2=\""
    << sx(xmax) << "\" y2=\"" << sy(0.0)
    << "\" stroke=\"#3a3f4b\" stroke-width=\"1\" stroke-dasharray=\"4 4\"/>\n";

  // initial shape: faint dashed grey
  f << "  <polygon points=\"" << polypoints(o0)
    << "\" fill=\"none\" stroke=\"#8a8f9a\" stroke-width=\"1.6\" "
       "stroke-dasharray=\"6 5\" opacity=\"0.9\"/>\n";

  // optimized shape: bold colored fill + outline
  f << "  <polygon points=\"" << polypoints(o1)
    << "\" fill=\"#4fc3f7\" fill-opacity=\"0.18\" stroke=\"#4fc3f7\" "
       "stroke-width=\"2.6\"/>\n";

  // control-point markers (upper P1, P2)
  auto marker = [&](double x, double y, const char *color) {
    f << "  <circle cx=\"" << sx(x) << "\" cy=\"" << sy(y)
      << "\" r=\"4\" fill=\"" << color << "\"/>\n";
  };
  marker(alpha_init[0], alpha_init[1], "#8a8f9a");
  marker(alpha_init[2], alpha_init[3], "#8a8f9a");
  marker(alpha_final[0], alpha_final[1], "#ffb74d");
  marker(alpha_final[2], alpha_final[3], "#ffb74d");

  // header + legend
  char buf[256];
  std::string title = capitalize(objective_name) + " shape optimization";
  f << "  <text x=\"20\" y=\"30\" fill=\"#e6e9ef\" font-size=\"20\" "
       "font-weight=\"700\">" << title << "</text>\n";
  snprintf(buf, sizeof(buf),
           "initial %s = %.4e      optimized %s = %.4e",
           objective_name.c_str(), value_init, objective_name.c_str(),
           value_final);
  f << "  <text x=\"20\" y=\"54\" fill=\"#aab1bd\" font-size=\"14\">" << buf
    << "</text>\n";

  // legend swatches (top-right)
  f << "  <line x1=\"" << (W - 250) << "\" y1=\"28\" x2=\"" << (W - 215)
    << "\" y2=\"28\" stroke=\"#8a8f9a\" stroke-width=\"1.6\" "
       "stroke-dasharray=\"6 5\"/>\n";
  f << "  <text x=\"" << (W - 208) << "\" y=\"33\" fill=\"#aab1bd\" "
       "font-size=\"13\">initial</text>\n";
  f << "  <line x1=\"" << (W - 250) << "\" y1=\"50\" x2=\"" << (W - 215)
    << "\" y2=\"50\" stroke=\"#4fc3f7\" stroke-width=\"2.6\"/>\n";
  f << "  <text x=\"" << (W - 208) << "\" y=\"55\" fill=\"#aab1bd\" "
       "font-size=\"13\">optimized</text>\n";

  f << "</svg>\n";
  f.close();
  return f.good() ? 0 : 1;
}
