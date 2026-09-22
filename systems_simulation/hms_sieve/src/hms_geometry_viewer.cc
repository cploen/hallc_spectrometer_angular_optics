// Geometry/convention display only. All positions come from build_geometry.py
// via scene.tsv, in laboratory (X,Y,Z), centimetres. No HMS constants or Euler
// rotations are maintained here. Markers and plane outlines are not matter.
//
// API references: Geant4 G4VisManager.hh (RegisterRunDurationUserVisAction),
// G4VisCommandsSceneAdd.cc (/vis/scene/add/userAction), and the visualization
// chapter of the Geant4 Application Developers Guide. Native build/runtime
// verification requires a Geant4 installation; Python closure tests do not
// establish that the Geant4 viewer has run.
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

class CoordinateDrawing final : public G4VUserVisAction {
 public:
  explicit CoordinateDrawing(const Scene& scene) : scene_(scene) {}
  void Draw() override {
    auto* vis = G4VVisManager::GetConcreteInstance();
    if (!vis) return;
    for (const auto& item : scene_.objects) {
      G4VisAttributes attributes(item.colour);
      attributes.SetLineWidth(2.0);  // Screen styling, not a physical dimension.
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
      // Keep the 81-hole grid legible: label its central marker and all named
      // axes, origins, rays, and one edge of each reference-plane outline.
      const bool isHole = item.name.find("nominal_hole_") == 0;
      const bool isEdge = item.name.find("_edge") != std::string::npos;
      if ((isHole && item.name != "nominal_hole_4_4") ||
          (isEdge && item.name.find("_edge0") == std::string::npos)) continue;
      std::string label = item.name;
      if (isEdge) label.erase(label.size() - 6);
      std::replace(label.begin(), label.end(), '_', ' ');
      G4Text text(label, item.b);
      text.SetScreenSize(12.0);
      text.SetLayout(G4Text::left);
      text.SetVisAttributes(attributes);
      vis->Draw(text);
    }
  }
 private:
  const Scene& scene_;
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
    CoordinateDrawing drawing(scene);  // Alive until after the vis manager.
    auto vis = std::make_unique<G4VisExecutive>("warnings");
    vis->Initialize();
    vis->RegisterRunDurationUserVisAction("HMS_coordinate_scene", &drawing, scene.extent());
    std::cout << "Drawing " << scene.objects.size() << " coordinate primitives from " << argv[1]
              << "\nLAB X/Y/Z, input cm; code-convention benchmark, not surveyed solids.\n";
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
