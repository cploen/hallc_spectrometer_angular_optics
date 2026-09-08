#ifndef HALLC_ANGULAR_SPECTROMETER_CONFIG_H
#define HALLC_ANGULAR_SPECTROMETER_CONFIG_H

// Campaign-selected centered-sieve profiles. See docs/HMS_SHMS_REVIEW.md.
// No ROOT dependency; the shared .def table is also read by Python.
#include <algorithm>
#include <cmath>
#include <fstream>
#include <vector>
#include <sstream>
#include <string>

namespace hallc {
struct Spectrometer {
  std::string name;
  std::string prefix;
  std::string cer;
  int nx=0, ny=0;
  double sieveDistance=0, dx=0, dy=0, deltaMin=0, deltaMax=0;
  double xMP=0, yMP=0, hbLinear=0, hbQuadratic=0;
  bool shms() const { return name == "SHMS"; }
  std::string lowerName() const { return shms() ? "shms" : "hms"; }
  double xs(int i) const { return (i-(nx-1)/2.0)*dx; }
  double ys(int i) const {
    // Keep the existing HMS operation order for regression parity.
    return shms() ? (i-5)*dy : (i-4)*0.6*2.54;
  }
  double hb(double delta) const { return hbLinear*delta+hbQuadratic*delta*delta; }
  bool acceptsDelta(double delta) const { return delta>deltaMin && delta<deltaMax; }
  double xMis(double angle) const {
    if (shms()) return xMP;
    const double a=std::min(std::abs(angle),50.0);
    return 0.1*(2.37-0.086*a+0.0012*a*a);
  }
  double yMis(double angle) const {
    if (shms()) return yMP;
    const double a=std::min(std::abs(angle),40.0);
    return 0.1*(0.52-0.012*a+0.002*a*a);
  }
  double yPlotMax() const { return shms() ? 9.0 : 7.0; }
  double xPlotMax() const { return shms() ? 14.0 : 12.5; }

  std::string branch(const std::string& suffix) const {
    return prefix + "." + suffix;
  }
  std::string cherenkovBranch() const {
    return branch(cer+".npeSum");
  }
};

inline Spectrometer profileForName(const std::string& name) {
#define HALLC_PROFILE(ARM,PREFIX,CER,NX,NY,L,DX,DY,DMIN,DMAX,XMP,YMP,HBL,HBQ) \
  if (name == #ARM) return {#ARM,#PREFIX,#CER,NX,NY,L,DX,DY,DMIN,DMAX,XMP,YMP,HBL,HBQ};
#include "spectrometer_profiles.def"
#undef HALLC_PROFILE
  return {};
}

struct Truth {
  double xptar, yptar, xtar, ytar;
};
inline Truth targetTruth(const Spectrometer& spec, double angleDeg,
                         double zfoil, double xs, double ys, double delta,
                         double reactx, double reacty, double xbpm) {
  const double angle=angleDeg*std::acos(-1.0)/180.0;
  const double c=std::cos(angle), s=std::sin(angle);
  const double xmis=spec.xMis(angleDeg), ymis=spec.yMis(angleDeg);
  const double L=spec.sieveDistance;
  Truth t;
  if (spec.shms()) {
    // Holly SHMS_optics/shms_optics.cpp, c7d70e8, physical-event equations.
    const double xv=-reacty-xmis, yv=-zfoil*s+reactx*c-ymis;
    const double zv=zfoil*c+reactx*s;
    t.xptar=(xs-xv)/(L-zv);
    t.yptar=(ys-spec.hb(delta)-yv)/(L-zv);
    t.xtar=xv-t.xptar*zv;
    t.ytar=yv-t.yptar*zv;
  } else {
    // Preserve the pre-switch HMS equations, including beam conventions.
    const double xbeam=-xbpm;
    const double yc=zfoil*s+xbeam*c-ymis;
    t.yptar=(ys-yc)/(L-zfoil*c);
    t.ytar=zfoil*(s-t.yptar*c)+xbeam*(c+t.yptar*s)-ymis;
    t.xptar=xs/(L-zfoil*c);
    t.xtar=-reacty-xmis-t.xptar*zfoil*c;
  }
  return t;
}

// Target positions and slice boundaries are run properties for either arm.
// In particular, -10,0,+10 cm is supported equally for HMS and SHMS.
struct RunMetadata {
  std::string opticsId;
  int sieveFlag=0;
  double angle=0;
  std::vector<double> foils, edges;
};
inline bool readRunMetadata(int run, RunMetadata& info, std::string& error,
                            const std::string& filename="DATfiles/list_of_optics_run.dat") {
  std::ifstream input(filename);
  std::string line;
  auto trim=[](std::string x) {
    const auto b=x.find_first_not_of(" \t\r\n");
    return b==std::string::npos ? std::string() : x.substr(b,x.find_last_not_of(" \t\r\n")-b+1);
  };
  while (std::getline(input,line)) {
    std::istringstream row(line); std::string field;
    std::vector<std::string> fields;
    while (std::getline(row,field,',')) fields.push_back(trim(field));
    if (fields.size()<6 || fields[0]!=std::to_string(run)) continue;
    info.opticsId=fields[1];
    int nfoils=0, nedges=0;
    if (!(std::istringstream(fields[2])>>info.angle) ||
        !(std::istringstream(fields[3])>>nfoils) ||
        !(std::istringstream(fields[4])>>info.sieveFlag) ||
        !(std::istringstream(fields[5])>>nedges)) break;
    info.foils.clear(); info.edges.clear();
    for (auto* values : {&info.foils,&info.edges}) {
      if (!std::getline(input,line)) break;
      std::istringstream items(line);
      while (std::getline(items,field,',')) {
        if (trim(field).empty()) continue;
        double value;
        if (!(std::istringstream(field)>>value) || !std::isfinite(value)) {
          error="Invalid foil/delta value in "+filename; return false;
        }
        values->push_back(value);
      }
    }
    if (nfoils>0 && int(info.foils.size())==nfoils && int(info.edges.size())==nedges && nedges>=2 &&
        std::adjacent_find(info.edges.begin(),info.edges.end(),[](double a,double b){return a>=b;})==info.edges.end())
      return true;
    break;
  }
  error="Missing or invalid optics metadata for run "+std::to_string(run)+" in "+filename;
  return false;
}

// Inspect whole path components, never the substring "HMS" inside "SHMS".
// Accept a campaign itself or one of its output paths. Reject conflicting
// campaign components and unnamed campaigns instead of defaulting to HMS.
inline bool spectrometerFromCampaign(const std::string& path,
                                    Spectrometer& result, std::string& error) {
  result = {};
  error.clear();
  std::istringstream stream(path);
  std::string part, name;
  while (std::getline(stream, part, '/')) {
    std::string candidate;
    if (part == "SHMS" || part.compare(0, 5, "SHMS_") == 0)
      candidate = "SHMS";
    else if (part == "HMS" || part.compare(0, 4, "HMS_") == 0)
      candidate = "HMS";
    if (candidate.empty()) continue;
    if (!name.empty() && candidate != name) {
      error = "Conflicting HMS/SHMS campaign components: " + path;
      return false;
    }
    name = candidate;
  }
  if (name.empty()) {
    error = "Campaign must be named HMS_<campaign> or SHMS_<campaign>: " + path;
    return false;
  }
  result = profileForName(name);
  return true;
}
}  // namespace hallc
#endif
