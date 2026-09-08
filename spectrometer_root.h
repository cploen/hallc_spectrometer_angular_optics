#ifndef HALLC_ANGULAR_SPECTROMETER_ROOT_H
#define HALLC_ANGULAR_SPECTROMETER_ROOT_H
#include "spectrometer_config.h"
#include <TTree.h>
#include <TString.h>
#include <TNamed.h>
#include <TDirectory.h>
#include <iostream>
#include <initializer_list>
namespace hallc {
inline bool loadSpectrometer(const TString& path, Spectrometer& spec) {
  std::string error;
  if (!spectrometerFromCampaign(path.Data(),spec,error)) {
    std::cerr << "ERROR: " << error << std::endl; return false;
  }
  std::cout << "Spectrometer: " << spec.name << "; centered sieve; delta=("
            << spec.deltaMin << "," << spec.deltaMax << ") percent" << std::endl;
  return true;
}
inline bool requireBranches(TTree* tree, std::initializer_list<std::string> names) {
  bool ok=true;
  for (const auto& name : names) {
    if (!tree->GetBranch(name.c_str()) && !tree->GetLeaf(name.c_str())) {
      std::cerr << "ERROR: missing selected-spectrometer branch " << name << std::endl;
      ok=false;
    }
  }
  return ok;
}
inline bool centeredSieveOnly(const Spectrometer& spec, int sieveFlag) {
  if (spec.shms() && sieveFlag!=1) {
    std::cerr << "ERROR: this SHMS version assumes centered sieve (SieveFlag=1); got "
              << sieveFlag << ". Shifted/no-sieve support is not implemented." << std::endl;
    return false;
  }
  return true;
}
inline bool loadRunMetadata(int run, RunMetadata& info, const char* filename="DATfiles/list_of_optics_run.dat") {
  std::string error;
  if (!readRunMetadata(run,info,error,filename)) {
    std::cerr << "ERROR: " << error << std::endl; return false;
  }
  return true;
}
inline void writeProfile(const Spectrometer& spec) {
  TNamed("hallc_spectrometer",spec.name.c_str()).Write();
  TNamed("hallc_geometry_profile",(spec.name+"_centered_v1").c_str()).Write();
  TNamed("hallc_geometry_assumption",spec.shms() ? "Centered sieve; SHMS offsets from 2019 survey note" : "Existing HMS geometry").Write();
}
}
#endif
