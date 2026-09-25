#include "G4Box.hh"
#include "G4Electron.hh"
#include "G4Element.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4Event.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4ParticleGun.hh"
#include "G4PhysicalConstants.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4Transform3D.hh"
#include "G4Tubs.hh"
#include "G4VModularPhysicsList.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4VUserActionInitialization.hh"
#include "G4VSensitiveDetector.hh"
#include "G4UserSteppingAction.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4SubtractionSolid.hh"
#include "G4Track.hh"
#include "G4RotationMatrix.hh"
#include "G4UImanager.hh"
#include "G4VProcess.hh"
#include "G4ios.hh"
#include "Randomize.hh"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>

namespace {
using Values = std::map<std::string, G4double>;

Values ReadValues(const std::string& path) {
  std::ifstream input(path);
  if (!input) throw std::runtime_error("cannot open Geant4 input: " + path);
  Values values;
  std::string key;
  G4double value = 0.;
  while (input >> key >> value) values[key] = value;
  return values;
}

G4double V(const Values& values, const std::string& key) {
  const auto found = values.find(key);
  if (found == values.end()) throw std::runtime_error("missing input key: " + key);
  return found->second;
}

G4ThreeVector Axis(const Values& v, const std::string& name) {
  return {V(v, name + "_x"), V(v, name + "_y"), V(v, name + "_z")};
}

G4RotationMatrix FrameRotation(const Values& v) {
  G4RotationMatrix rotation;
  rotation.set(Axis(v, "hms_x"), Axis(v, "hms_y"), Axis(v, "hms_z"));
  return rotation;
}

G4Material* MakeHeavyMetal() {
  auto* nist = G4NistManager::Instance();
  auto* material = new G4Material("HD17_Tungsten_Alloy", 17.0 * g / cm3, 3);
  material->AddElement(nist->FindOrBuildElement("W"), 0.90);
  material->AddElement(nist->FindOrBuildElement("Ni"), 0.06);
  material->AddElement(nist->FindOrBuildElement("Cu"), 0.04);
  return material;
}

G4Material* MakeCarbon(G4double density) {
  auto* nist = G4NistManager::Instance();
  auto* material = new G4Material("Carbon05_ProvisionalDensity", density * g / cm3, 1);
  material->AddElement(nist->FindOrBuildElement("C"), 1.0);
  return material;
}

class EntranceSensitive : public G4VSensitiveDetector {
 public:
  EntranceSensitive(const G4String& name, Values values)
      : G4VSensitiveDetector(name), values_(std::move(values)), out_("q1_entrance_hits.tsv") {
    out_ << "event\ttrack\tparent\tPDG\tcreator_process\tkinetic_GeV\tmomentum_GeV_c\tdelta_from_6p667\tin_delta_window\tx_lab_cm\ty_lab_cm\tz_lab_cm\tdir_x\tdir_y\tdir_z\n";
    out_ << std::setprecision(12);
  }

  G4bool ProcessHits(G4Step* step, G4TouchableHistory*) override {
    auto* track = step->GetTrack();
    const auto* event = G4RunManager::GetRunManager()->GetCurrentEvent();
    const auto* creator = track->GetCreatorProcess();
    const auto p = track->GetMomentum().mag() / GeV;
    const auto position = step->GetPreStepPoint()->GetPosition();
    const auto direction = track->GetMomentumDirection();
    out_ << event->GetEventID() << '\t' << track->GetTrackID() << '\t'
         << track->GetParentID() << '\t' << track->GetDefinition()->GetPDGEncoding() << '\t'
         << (creator ? creator->GetProcessName() : "primary") << '\t'
         << track->GetKineticEnergy() / GeV << '\t' << p << '\t'
         << (p / V(values_, "hms_p0_GeV_c") - 1.) << '\t'
         << (std::abs(p / V(values_, "hms_p0_GeV_c") - 1.) <= V(values_, "hms_delta_half_width")) << '\t'
         << position.x() / cm << '\t'
         << position.y() / cm << '\t' << position.z() / cm << '\t'
         << direction.x() << '\t' << direction.y() << '\t' << direction.z() << '\n';
    return true;
  }

 private:
  Values values_;
  std::ofstream out_;
};

class Detector : public G4VUserDetectorConstruction {
 public:
  explicit Detector(Values values) : v_(std::move(values)) {}

  G4VPhysicalVolume* Construct() override {
    auto* nist = G4NistManager::Instance();
    auto* vacuum = nist->FindOrBuildMaterial("G4_Galactic");
    auto* worldSolid = new G4Box("World", 400. * cm, 400. * cm, 400. * cm);
    auto* worldLogic = new G4LogicalVolume(worldSolid, vacuum, "World_LV");
    auto* world = new G4PVPlacement(nullptr, {}, worldLogic, "World", nullptr, false, 0, true);

    auto* carbon = MakeCarbon(V(v_, "target_density_g_cm3"));
    auto* targetSolid = new G4Box("Carbon05_ProvisionalPatch",
        V(v_, "target_half_x_cm") * cm, V(v_, "target_half_y_cm") * cm,
        V(v_, "target_thickness_cm") * cm / 2.);
    auto* targetLogic = new G4LogicalVolume(targetSolid, carbon, "Carbon05_Target_LV");
    new G4PVPlacement(nullptr, {0., 0., 0.}, targetLogic, "Carbon05_Target", worldLogic, false, 0, true);

    auto* alloy = MakeHeavyMetal();
    const auto halfX = V(v_, "sieve_half_x_cm") * cm;
    const auto halfY = V(v_, "sieve_half_y_cm") * cm;
    const auto halfZ = V(v_, "sieve_half_z_cm") * cm;
    G4VSolid* sieveSolid = new G4Box("HMS_Sieve_Blank", halfX, halfY, halfZ);
    auto identity = G4RotationMatrix();
    const auto pitchX = V(v_, "pitch_x_cm") * cm;
    const auto pitchY = V(v_, "pitch_y_cm") * cm;
    const auto radius = V(v_, "hole_diameter_cm") * cm / 2.;
    const auto centerRadius = V(v_, "center_hole_diameter_cm") * cm / 2.;
    const auto topOffset = V(v_, "grid_top_offset_cm") * cm;
    for (G4int i = 0; i < 9; ++i) {
      for (G4int j = 0; j < 9; ++j) {
        if (v_.count("blocked_" + std::to_string(i) + "_" + std::to_string(j))) continue;
        const auto holeRadius = (i == 4 && j == 4) ? centerRadius : radius;
        auto* bore = new G4Tubs("HMS_Sieve_Bore", 0., holeRadius, halfZ + 0.01 * cm, 0., twopi);
        const auto x = (i - 4) * pitchX;
        const auto y = (j - 4) * pitchY + topOffset;
        sieveSolid = new G4SubtractionSolid("HMS_Sieve_With_Bore_" + std::to_string(i) + "_" + std::to_string(j),
            sieveSolid, bore, G4Transform3D(identity, {x, y, 0.}));
      }
    }
    auto* sieveLogic = new G4LogicalVolume(sieveSolid, alloy, "HMS_Sieve_Plate_LV");
    auto* rotation = new G4RotationMatrix(FrameRotation(v_));
    const G4ThreeVector sieveCenter(V(v_, "sieve_center_x_lab_cm") * cm,
                                    V(v_, "sieve_center_y_lab_cm") * cm,
                                    V(v_, "sieve_center_z_lab_cm") * cm);
    new G4PVPlacement(rotation, sieveCenter, sieveLogic, "HMS_Sieve_Plate", worldLogic, false, 0, true);

    const auto apertureRadius = V(v_, "q1_aperture_radius_cm") * cm;
    auto* apertureSolid = new G4Tubs("Q1_Entrance_Sensitive_Aperture", 0., apertureRadius,
                                     0.005 * cm, 0., twopi);
    auto* apertureLogic = new G4LogicalVolume(apertureSolid, vacuum, "Q1_Entrance_Aperture_LV");
    new G4PVPlacement(rotation,
        {V(v_, "q1_score_x_lab_cm") * cm, V(v_, "q1_score_y_lab_cm") * cm, V(v_, "q1_score_z_lab_cm") * cm},
        apertureLogic, "Q1_Entrance_Sensitive_Aperture", worldLogic, false, 0, true);
    return world;
  }

  void ConstructSDandField() override {
    auto* detector = new EntranceSensitive("Q1EntranceSensitive", v_);
    G4SDManager::GetSDMpointer()->AddNewDetector(detector);
    SetSensitiveDetector("Q1_Entrance_Aperture_LV", detector);
  }

 private:
  Values v_;
};

class Beam : public G4VUserPrimaryGeneratorAction {
 public:
  explicit Beam(Values values) : v_(std::move(values)), gun_(1) {
    gun_.SetParticleDefinition(G4Electron::Definition());
    gun_.SetParticleEnergy(V(v_, "beam_energy_GeV") * GeV);
    gun_.SetParticleMomentumDirection({0., 0., 1.});
  }

  void GeneratePrimaries(G4Event* event) override {
    const auto x = V(v_, "beam_x_cm") + (2. * G4UniformRand() - 1.) * V(v_, "raster_x_half_cm");
    const auto y = V(v_, "beam_y_cm") + (2. * G4UniformRand() - 1.) * V(v_, "raster_y_half_cm");
    gun_.SetParticlePosition({x * cm, y * cm, V(v_, "source_z_cm") * cm});
    gun_.GeneratePrimaryVertex(event);
  }

 private:
  Values v_;
  G4ParticleGun gun_;
};

class Steps : public G4UserSteppingAction {
 public:
  Steps() : out_("material_steps.tsv") {
    out_ << "event\ttrack\tparent\tPDG\tvolume\tprocess\tkinetic_GeV\tmomentum_GeV_c\tx_lab_cm\ty_lab_cm\tz_lab_cm\n";
    out_ << std::setprecision(12);
  }

  void UserSteppingAction(const G4Step* step) override {
    const auto* pre = step->GetPreStepPoint();
    const auto* volume = pre->GetPhysicalVolume();
    if (!volume || (volume->GetName() != "HMS_Sieve_Plate" && volume->GetName() != "Carbon05_Target")) return;
    const auto* post = step->GetPostStepPoint();
    const auto* process = post->GetProcessDefinedStep();
    const auto* track = step->GetTrack();
    const auto* event = G4RunManager::GetRunManager()->GetCurrentEvent();
    const auto pos = post->GetPosition();
    out_ << event->GetEventID() << '\t' << track->GetTrackID() << '\t' << track->GetParentID() << '\t'
         << track->GetDefinition()->GetPDGEncoding() << '\t' << volume->GetName() << '\t'
         << (process ? process->GetProcessName() : "none") << '\t'
         << track->GetKineticEnergy() / GeV << '\t' << track->GetMomentum().mag() / GeV << '\t'
         << pos.x() / cm << '\t' << pos.y() / cm << '\t' << pos.z() / cm << '\n';
  }

 private:
  std::ofstream out_;
};

class Actions : public G4VUserActionInitialization {
 public:
  explicit Actions(Values values) : v_(std::move(values)) {}
  void Build() const override {
    SetUserAction(new Beam(v_));
    SetUserAction(new Steps());
  }

 private:
  Values v_;
};

class Physics : public G4VModularPhysicsList {
 public:
  Physics() { RegisterPhysics(new G4EmStandardPhysics_option4); }
  void SetCuts() override { SetCutsWithDefault(); }
};
}  // namespace

int main(int argc, char** argv) {
  if (argc < 2 || argc > 3) {
    G4cerr << "Usage: hms_sieve_scattering geant4_input.tsv [run.mac]\n";
    return 2;
  }
  try {
    const auto values = ReadValues(argv[1]);
    auto* runManager = new G4RunManager;
    runManager->SetUserInitialization(new Detector(values));
    auto* physics = new Physics;
    physics->SetDefaultCutValue(1.0 * mm);
    runManager->SetUserInitialization(physics);
    runManager->SetUserInitialization(new Actions(values));
    runManager->Initialize();
    auto* ui = G4UImanager::GetUIpointer();
    if (argc == 3) {
      const auto command = std::string("/control/execute ") + argv[2];
      if (ui->ApplyCommand(command) != 0) throw std::runtime_error("Geant4 macro failed");
    } else {
      const auto command = std::string("/run/beamOn ") + std::to_string(static_cast<long>(V(values, "events")));
      if (ui->ApplyCommand(command) != 0) throw std::runtime_error("Geant4 run failed");
    }
    delete runManager;
  } catch (const std::exception& error) {
    G4cerr << "hms_sieve_scattering: " << error.what() << G4endl;
    return 1;
  }
  return 0;
}
