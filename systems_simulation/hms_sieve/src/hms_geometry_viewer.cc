// Geometry/convention display only. All positions come from build_geometry.py
// via scene.tsv, in laboratory (X,Y,Z), centimetres. No HMS constants or Euler
// rotations are maintained here. Markers and plane outlines are not matter.
//
// API references: Geant4 G4VisManager.hh (RegisterRunDurationUserVisAction),
// G4VisCommandsSceneAdd.cc (/vis/scene/add/userAction), and the visualization
// chapter of the Geant4 Application Developers Guide. Native build/runtime
// baseline viewer ran on iFarm (see IFARM_VALIDATION.md). This layered revision
// still needs its own native smoke test; Python tests do not validate rendering.
#include "annotation_layout.hh"
#include "G4Box.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4Geantino.hh"
#include "G4LogicalVolume.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4Point3D.hh"
#include "G4Polyline.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4Text.hh"
#include "G4ThreeVector.hh"
#include "G4UImanager.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VUserPhysicsList.hh"
#include "G4VUserVisAction.hh"
#include "G4VVisManager.hh"
#include "G4VisAttributes.hh"
#include "G4VisExecutive.hh"
#include "G4VisExtent.hh"
#if HMS_WITH_INTERACTIVE
#include "G4UIExecutive.hh"
#endif

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
struct Primitive {
  std::string kind, name;
  G4Point3D a, b;
  G4Colour colour;
};

struct Scene {
  std::vector<Primitive> objects;
  double xmin = 0, xmax = 0, ymin = 0, ymax = 0, zmin = 0, zmax = 0;

  void include(const G4Point3D& p) {
    xmin = std::min(xmin, p.x()); xmax = std::max(xmax, p.x());
    ymin = std::min(ymin, p.y()); ymax = std::max(ymax, p.y());
    zmin = std::min(zmin, p.z()); zmax = std::max(zmax, p.z());
  }

  G4VisExtent extent() const {
    return G4VisExtent(xmin, xmax, ymin, ymax, zmin, zmax);
  }

  double worldHalfSize() const {
    // Visualization envelope derived from scene bounds, not an equipment size.
    return 1.1 * std::max({std::abs(xmin), std::abs(xmax), std::abs(ymin),
                         std::abs(ymax), std::abs(zmin), std::abs(zmax)}) + cm;
  }
};

Scene readScene(const std::string& path) {
  std::ifstream input(path);
  if (!input) throw std::runtime_error("Cannot read generated scene: " + path);
  Scene scene;
  std::string line;
  std::size_t number = 0;
  while (std::getline(input, line)) {
    ++number;
    std::istringstream row(line);
    Primitive item;
    if (!(row >> item.kind) || item.kind[0] == '#') continue;
    double ax = 0, ay = 0, az = 0, bx = 0, by = 0, bz = 0;
    double red = 0, green = 0, blue = 0;
    bool valid = static_cast<bool>(row >> item.name >> ax >> ay >> az);
    if (item.kind == "line") valid = valid && static_cast<bool>(row >> bx >> by >> bz);
    else if (item.kind == "point") { bx = ax; by = ay; bz = az; }
    else valid = false;
    valid = valid && static_cast<bool>(row >> red >> green >> blue);
    std::string extra;
    if (!valid || (row >> extra))
      throw std::runtime_error("Malformed scene row " + std::to_string(number));
    for (double coordinate : {ax, ay, az, bx, by, bz, red, green, blue})
      if (!std::isfinite(coordinate)) throw std::runtime_error("Nonfinite scene coordinate");
    for (double component : {red, green, blue})
      if (component < 0 || component > 1) throw std::runtime_error("Colour outside [0,1]");
    item.a = G4Point3D(ax * cm, ay * cm, az * cm);
    item.b = G4Point3D(bx * cm, by * cm, bz * cm);
    item.colour = G4Colour(red, green, blue);
    scene.include(item.a); scene.include(item.b);
    scene.objects.push_back(item);
  }
  if (scene.objects.empty()) throw std::runtime_error("Generated scene is empty");
  return scene;
}

class VacuumWorld final : public G4VUserDetectorConstruction {
 public:
  explicit VacuumWorld(double halfSize) : halfSize_(halfSize) {}
  G4VPhysicalVolume* Construct() override {
    auto* box = new G4Box("CoordinateDisplayWorldSolid", halfSize_, halfSize_, halfSize_);
    auto* vacuum = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
    auto* logical = new G4LogicalVolume(box, vacuum, "CoordinateDisplayVacuum");
    logical->SetVisAttributes(G4VisAttributes::GetInvisible());
    return new G4PVPlacement(nullptr, G4ThreeVector(), logical,
                             "CoordinateDisplayWorld", nullptr, false, 0);
  }
 private:
  double halfSize_;
};

// Only the initialization required by the G4 kernel. No EM, hadronic, decay,
// detector, field, beam, or scattering model is defined, and no events are run.
class DisplayOnlyPhysics final : public G4VUserPhysicsList {
 public:
  void ConstructParticle() override { G4Geantino::GeantinoDefinition(); }
  void ConstructProcess() override { AddTransportation(); }
  void SetCuts() override {}
};

// Fixed screen offsets below are typography, not physical geometry parameters.
void label3D(G4VVisManager* vis, const std::string& label, const G4Point3D& position,
             const G4Colour& colour, double dx = 6, double dy = 0,
             G4Text::Layout layout = G4Text::left) {
  G4Text text(label, position);
  text.SetScreenSize(13.0);
  text.SetLayout(layout);
  text.SetOffset(dx, dy);
  G4VisAttributes attributes(colour);
  text.SetVisAttributes(attributes);
  vis->Draw(text);
}

void label2D(G4VVisManager* vis, const std::string& label, double x, double y,
             const G4Colour& colour = G4Colour(0.15, 0.15, 0.15), double size = 13) {
  G4Text text(label, G4Point3D(x, y, 0));
  text.SetScreenSize(size);
  text.SetLayout(G4Text::left);
  G4VisAttributes attributes(colour);
  text.SetVisAttributes(attributes);
  vis->Draw2D(text);  // Normalized screen coordinates, independent of camera.
}

void drawPrimitive(G4VVisManager* vis, const Primitive& item) {
  G4VisAttributes attributes(item.colour);
  attributes.SetLineWidth(2.0);
  if (item.kind == "line") {
    G4Polyline line;
    line.push_back(item.a); line.push_back(item.b);
    line.SetVisAttributes(attributes);
    vis->Draw(line);
  } else {
    G4Circle marker(item.a);
    marker.SetScreenSize(6.0);
    marker.SetFillStyle(G4Circle::filled);
    marker.SetVisAttributes(attributes);
    vis->Draw(marker);
  }
}

bool isAxis(const std::string& name) {
  return name == "laboratory_beam_axis" || name == "HMS_central_axis" ||
         hms_display::starts(name, "LAB_") ||
         hms_display::starts(name, "HMS_TRANSPORT_") ||
         hms_display::starts(name, "SIEVE_LOCAL_");
}

void drawAxisArrow(G4VVisManager* vis, const Primitive& item) {
  const G4ThreeVector start(item.a.x(), item.a.y(), item.a.z());
  const G4ThreeVector end(item.b.x(), item.b.y(), item.b.z());
  const G4ThreeVector delta = end - start;
  if (delta.mag2() == 0) return;
  const G4ThreeVector direction = delta.unit();
  G4ThreeVector side = direction.cross(G4ThreeVector(0, 1, 0));
  if (side.mag2() < 1.e-12) side = direction.cross(G4ThreeVector(1, 0, 0));
  side = side.unit();

  const bool beam = item.name == "laboratory_beam_axis";
  const bool hmsAxis = item.name == "HMS_central_axis";
  const double fraction = beam ? .55 : hmsAxis ? .75 : .90;
  const double headLength = (beam || hmsAxis ? 5.0 : 1.8) * cm;
  const double halfWidth = (beam || hmsAxis ? 2.0 : .7) * cm;
  const G4ThreeVector tip = start + fraction * delta;
  const G4ThreeVector base = tip - headLength * direction;
  const G4ThreeVector left = base + halfWidth * side;
  const G4ThreeVector right = base - halfWidth * side;

  G4Polyline arrow;
  arrow.push_back(G4Point3D(tip.x(), tip.y(), tip.z()));
  arrow.push_back(G4Point3D(left.x(), left.y(), left.z()));
  arrow.push_back(G4Point3D(tip.x(), tip.y(), tip.z()));
  arrow.push_back(G4Point3D(right.x(), right.y(), right.z()));
  G4VisAttributes attributes(item.colour);
  attributes.SetLineWidth(2.0);
  arrow.SetVisAttributes(attributes);
  vis->Draw(arrow);
}

class CoordinateDrawing final : public G4VUserVisAction {
 public:
  CoordinateDrawing(const Scene& scene, hms_display::Layer layer)
      : scene_(scene), layer_(layer) {}
  void Draw() override {
    using namespace hms_display;
    auto* vis = G4VVisManager::GetConcreteInstance();
    if (!vis) return;
    if (layer_ == Layer::code || layer_ == Layer::branches) {
      drawNames(vis);
      return;
    }
    if (layer_ == Layer::labels) {
      drawPlaneNames(vis);
      return;
    }
    for (const auto& item : scene_.objects) {
      if (primitiveLayer(item.name) != layer_) continue;
      drawPrimitive(vis, item);
      if (isAxis(item.name)) drawAxisArrow(vis, item);
      if (layer_ == Layer::lab || layer_ == Layer::hms || layer_ == Layer::sieve) {
        const std::string frame = layer_ == Layer::lab ? "LAB " :
                                  layer_ == Layer::hms ? "HMS " : "SIEVE ";
        std::string axis(1, item.name.back());
        if (layer_ == Layer::lab) axis[0] = static_cast<char>(axis[0] - 'a' + 'A');
        const double offset = layer_ == Layer::lab ? 20 : layer_ == Layer::hms ? 0 : -18;
        label3D(vis, frame + axis, item.b, item.colour, 5, offset);
      }
      if (layer_ == Layer::rays && item.kind == "line") {
        const int n = item.name[0] - 'A';
        G4Point3D middle((item.a.x()+item.b.x())/2, (item.a.y()+item.b.y())/2,
                         (item.a.z()+item.b.z())/2);
        label3D(vis, item.name.substr(0,1), middle, item.colour, 6, 22-22*n);
      }
    }
    if (layer_ == Layer::base) {
      label2D(vis, "HMS sieve-slit study | schematic geometry", -.96, -.88);
      label2D(vis, "Nominal centers and projection planes", -.96, -.95, G4Colour(.4,.4,.4), 12);
    }
    if (layer_ == Layer::rays) {
      // Connect existing sieve intersection markers. These are residual
      // annotations, not new rays, displaced points, or revised constructions.
      const auto& a = object("A_endpoint_sieve_intersection");
      for (const auto* name : {"B_fit_target_sieve_intersection", "C_HCANA_coordinate_step_sieve_intersection"}) {
        const auto& other = object(name);
        G4Polyline connector;
        connector.push_back(a.a); connector.push_back(other.a);
        G4VisAttributes attributes(G4Colour(.5,.5,.5));
        connector.SetVisAttributes(attributes);
        vis->Draw(connector);
      }
      label2D(vis, "A: endpoint ray", .10, -.58, object("A_endpoint").colour);
      label2D(vis, "B: fit target construction", .10, -.66, object("B_fit_target").colour);
      label2D(vis, "C: HCANA coordinate step (analytic inputs)", .10, -.74,
              object("C_HCANA_coordinate_step").colour);
      label2D(vis, "Dots: sieve intersections; gray: residual connectors", .10, -.82, G4Colour(.4,.4,.4), 12);
    }
  }
 private:
  const Primitive& object(const std::string& name) const {
    for (const auto& item : scene_.objects) if (item.name == name) return item;
    throw std::runtime_error("Annotation requires scene primitive: " + name);
  }

  void drawPlaneNames(G4VVisManager* vis) const {
    const auto& lab = object("lab_Z0_reference_plane_edge2");
    const auto& hms = object("HMS_target_z0_reference_plane_edge2");
    const auto& sieve = object("HCANA_sieve_zL_projection_plane_edge2");
    // Short labels, with different screen offsets at the nearby target planes.
    label3D(vis, "Lab Z=0", lab.b, lab.colour, 6, 28);
    label3D(vis, "HMS z=0 [T]", hms.b, hms.colour, 6, 0);
    label3D(vis, "Projection plane [S]", sieve.b, sieve.colour, -6, 12, G4Text::right);
    const auto& beam = object("laboratory_beam_axis");
    const auto& axis = object("HMS_central_axis");
    label3D(vis, "Beam", beam.b, beam.colour, -6, 12, G4Text::right);
    label3D(vis, "HMS axis", axis.b, axis.colour, -6, -18, G4Text::right);
  }

  void drawNames(G4VVisManager* vis) const {
    using namespace hms_display;
    const bool code = layer_ == Layer::code;
    const double x = code ? -.96 : .08;
    label2D(vis, code ? "Code names / concepts" : "ROOT aliases (= same stored quantity)", x, .92,
            G4Colour(.1,.1,.1), 15);
    label2D(vis, "Naming key", x, .84, G4Colour(.4,.4,.4), 12);
    // Keep the two naming columns compact and above the main horizontal axes.
    double y = .73;
    for (const auto& q : quantities) {
      label2D(vis, code ? codeLabel(q) : aliasLabel(q), x, y);
      y -= .055;
    }
    if (code) {
      label2D(vis, "Constructed quantities (not saved aliases)", x, .10, G4Colour(.5,0,.4), 14);
      label2D(vis, "[T] xtarT, ytarT, xptarT, yptarT: construction B", x, .04);
      label2D(vis, "[S] xsT, ysT: assigned nominal hole center", x, -.02);
      label2D(vis, "ztarT: nominal foil Z; distinct from reactz / ztar", x, -.08);
      label2D(vis, "Construction B: fit-side HMS target coordinates", x, -.14, G4Colour(.4,.4,.4), 12);
    } else {
      label2D(vis, "[T] target plane/direction; [S] sieve projection", x, -.30, G4Colour(.4,.4,.4), 12);
      label2D(vis, "[V] lab reaction coordinates; [R] raster/BPM", x, -.38, G4Colour(.4,.4,.4), 12);
      label2D(vis, "V and R are different inputs; provider mapping TBD", x, -.46, G4Colour(.4,.4,.4), 12);
    }
    // [T]/[S] are anchored to the SAME planes for both naming layers. No
    // branch name is attached to an A/B/C endpoint as if it were a saved event.
    // V and R deliberately have no invented stored-vertex/BPM-plane marker.
  }

  const Scene& scene_;
  hms_display::Layer layer_;
};

void command(const std::string& value) {
  const auto result = G4UImanager::GetUIpointer()->ApplyCommand(value);
  if (result != 0) throw std::runtime_error("Geant4 command failed (" +
                          std::to_string(result) + "): " + value);
}
}  // namespace

int main(int argc, char** argv) {
  if (argc < 3 || argc > 4 || (argc == 4 && std::string(argv[3]) != "--interactive")) {
    std::cerr << "Usage: hms_geometry_viewer scene.tsv view.mac [--interactive]\n"
                 "       hms_geometry_viewer scene.tsv export_vrml.mac\n";
    return 2;
  }
  try {
    const auto scene = readScene(argv[1]);
#if HMS_WITH_INTERACTIVE
    std::unique_ptr<G4UIExecutive> ui;
    if (argc == 4) {
      int uiArgc = 1;  // Application arguments are not G4 UI arguments.
      ui = std::make_unique<G4UIExecutive>(uiArgc, argv);
    }
#else
    if (argc == 4) throw std::runtime_error("Rebuild with HMS_WITH_INTERACTIVE=ON");
#endif
    auto run = std::make_unique<G4RunManager>();
    run->SetUserInitialization(new VacuumWorld(scene.worldHalfSize()));
    run->SetUserInitialization(new DisplayOnlyPhysics());
    run->Initialize();
    std::vector<std::unique_ptr<CoordinateDrawing>> drawings;  // Outlive vis manager.
    for (const auto& layer : hms_display::layers)
      drawings.push_back(std::make_unique<CoordinateDrawing>(scene, layer.layer));
    auto vis = std::make_unique<G4VisExecutive>("warnings");
    vis->Initialize();
    for (std::size_t i = 0; i < hms_display::layers.size(); ++i)
      vis->RegisterRunDurationUserVisAction(hms_display::layers[i].name, drawings[i].get(), scene.extent());
    std::cout << "Drawing " << scene.objects.size() << " coordinate primitives from " << argv[1]
              << "\nLAB X/Y/Z, input cm; code-convention benchmark, not surveyed solids.\n"
              << "Layers: HMS_base, HMS_labels, HMS_lab_frame, HMS_transport_frame,\n"
                 "        HMS_sieve_frame, HMS_code_names, HMS_branch_names, HMS_rays.\n";
    command(std::string("/control/execute ") + argv[2]);
#if HMS_WITH_INTERACTIVE
    if (ui) ui->SessionStart();
#endif
  } catch (const std::exception& error) {
    std::cerr << "hms_geometry_viewer: " << error.what() << '\n';
    return 1;
  }
  return 0;
}
