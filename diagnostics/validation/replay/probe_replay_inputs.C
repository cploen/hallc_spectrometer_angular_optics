// HMS/SHMS preparation: read-only branch and cut survey, not an optics fit.
// Samples the FIRST maxEvents entries; this is not a representative full-run
// efficiency measurement. Prints TSV records to stdout; creates no ROOT files.
// See docs/HMS_SHMS_REVIEW.md for interpretation and the pending run settings.
#include "../../../spectrometer_config.h"
#include <TFile.h>
#include <TTree.h>
#include <TLeaf.h>
#include <TBranch.h>
#include <TKey.h>
#include <TString.h>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <vector>

namespace hallc_probe {
struct Column {
  std::string name;
  bool required;
  TLeaf* leaf;
  double value = std::numeric_limits<double>::quiet_NaN();
  Long64_t invalid = 0;
  std::vector<double> samples;
  Column(const std::string& n, bool r, TLeaf* l) : name(n), required(r), leaf(l) {}
};
double percentile(const std::vector<double>& sorted, double fraction) {
  const double position = fraction * (sorted.size() - 1);
  const auto lo = static_cast<std::size_t>(position);
  const auto hi = std::min(lo + 1, sorted.size() - 1);
  return sorted[lo] + (sorted[hi] - sorted[lo]) * (position - lo);
}
}

int probe_replay_inputs(TString campaign, TString inputRoot, Long64_t maxEvents=50000) {
  using namespace hallc_probe;
  hallc::Spectrometer spec;
  std::string error;
  if (!hallc::spectrometerFromCampaign(campaign.Data(), spec, error)) {
    std::cerr << "ERROR\t" << error << '\n';
    return 2;
  }
  if (maxEvents <= 0) {
    std::cerr << "ERROR\tmaxEvents must be positive\n";
    return 2;
  }
  std::cout << std::setprecision(10)
            << "FILE\t" << inputRoot << "\nCAMPAIGN\t" << campaign
            << "\nSPECTROMETER\t" << spec.name << "\n";
  std::unique_ptr<TFile> file(TFile::Open(inputRoot, "READ"));
  if (!file || file->IsZombie()) {
    std::cerr << "ERROR\tCannot open replay\t" << inputRoot << '\n';
    return 2;
  }
  TIter next(file->GetListOfKeys());
  while (auto* key = static_cast<TKey*>(next()))
    std::cout << "OBJECT\t" << key->GetName() << '\t' << key->GetClassName() << '\n';
  auto* tree = dynamic_cast<TTree*>(file->Get("T"));
  if (!tree) {
    std::cerr << "ERROR\tNo event tree T\n";
    return 2;
  }
  tree->SetBranchStatus("*", 0);
  std::vector<Column> columns;
  auto add = [&](const std::string& name, bool required) {
    auto* leaf = tree->GetLeaf(name.c_str());
    columns.emplace_back(name, required, leaf);
    std::cout << "BRANCH\t" << name << '\t' << (leaf ? "present" : "missing")
              << '\t' << (required ? "required" : "optional") << '\n';
    if (leaf) {
      tree->SetBranchStatus(leaf->GetBranch()->GetName(), 1);
      if (leaf->GetLeafCount())
        tree->SetBranchStatus(leaf->GetLeafCount()->GetBranch()->GetName(), 1);
    }
  };
  // First three indices are used for the comparisons below.
  add(spec.cherenkovBranch(), true);
  add(spec.branch("cal.etottracknorm"), true);
  add(spec.branch("gtr.dp"), true);
  for (const auto* suffix : {"gtr.x", "gtr.y", "gtr.th", "gtr.ph",
                            "dc.x_fp", "dc.xp_fp", "dc.y_fp", "dc.yp_fp",
                            "react.x", "react.y", "react.z",
                            "extcor.xsieve", "extcor.ysieve",
                            "rb.raster.fr_xbpm_tar", "rb.raster.fr_ybpm_tar"})
    add(spec.branch(suffix), true);
  add(spec.branch("cal.etracknorm"), false);  // Different observable in older configs.
  if (spec.name == "SHMS") add(spec.branch("hgcer.npeSum"), false);
  const Long64_t requested = std::min(maxEvents, tree->GetEntries());
  Long64_t read = 0, finitePID = 0, ridge = 0, currentPID = 0, strictPID = 0;
  Long64_t currentWindow = 0, strictNarrow = 0, strictWide = 0, readErrors = 0;
  Long64_t currentNominalSHMS = 0, strictNominalSHMS = 0;
  for (Long64_t entry=0; entry<requested; ++entry) {
    if (tree->GetEntry(entry) <= 0) { ++readErrors; continue; }
    ++read;
    for (auto& column : columns) {
      column.value = std::numeric_limits<double>::quiet_NaN();
      if (!column.leaf) continue;
      if (column.leaf->GetNdata() == 1) column.value = column.leaf->GetValue();
      if (std::isfinite(column.value)) column.samples.push_back(column.value);
      else ++column.invalid;
    }
    const double npe = columns[0].value, cal = columns[1].value, dp = columns[2].value;
    if (std::isfinite(npe) && std::isfinite(dp) && npe>2 && dp>-10 && dp<10) ++ridge;
    if (!(std::isfinite(npe) && std::isfinite(cal) && std::isfinite(dp))) continue;
    ++finitePID;
    if (npe>6 && cal>0.65) {
      ++currentPID;
      if (dp>-10 && dp<10) ++currentWindow;
      if (dp>-10 && dp<22) ++currentNominalSHMS;
    }
    if (npe>6 && cal>0.8) {
      ++strictPID;
      if (dp>-10 && dp<10) ++strictNarrow;
      if (dp>-10 && dp<22) ++strictNominalSHMS;
      if (dp>-15 && dp<24) ++strictWide;
    }
  }
  std::cout << "ENTRIES\ttotal\t" << tree->GetEntries()
            << "\nENTRIES\tfirst_requested\t" << requested
            << "\nENTRIES\tread\t" << read
            << "\nENTRIES\tread_errors\t" << readErrors
            << "\nENTRIES\tfinite_npe_cal_dp\t" << finitePID << '\n';
  std::cout << "CUT\tridge_npe>2_-10<dp<10\t" << ridge
            << "\nCUT\tnpe>6_cal>0.65_no_dp_window\t" << currentPID
            << "\nCUT\tnpe>6_cal>0.65_-10<dp<10\t" << currentWindow
            << "\nCUT\tnpe>6_cal>0.65_-10<dp<22\t" << currentNominalSHMS
            << "\nCUT\tnpe>6_cal>0.8_no_dp_window\t" << strictPID
            << "\nCUT\tnpe>6_cal>0.8_-10<dp<10\t" << strictNarrow
            << "\nCUT\tnpe>6_cal>0.8_-10<dp<22\t" << strictNominalSHMS
            << "\nCUT\tnpe>6_cal>0.8_-15<dp<24\t" << strictWide << '\n';
  std::cout << "RANGE_HEADER\tbranch\tfinite\tinvalid_or_nonscalar\tmin\tp05\tmedian\tp95\tmax\n";
  bool complete = requested>0 && readErrors==0;
  for (auto& column : columns) {
    if (column.required && (!column.leaf || column.samples.empty() || column.invalid>0))
      complete = false;
    if (!column.leaf) continue;
    std::cout << "RANGE\t" << column.name << '\t' << column.samples.size()
              << '\t' << column.invalid;
    if (!column.samples.empty()) {
      std::sort(column.samples.begin(), column.samples.end());
      std::cout << '\t' << column.samples.front()
                << '\t' << percentile(column.samples, .05)
                << '\t' << percentile(column.samples, .50)
                << '\t' << percentile(column.samples, .95)
                << '\t' << column.samples.back();
    }
    std::cout << '\n';
  }
  std::cout << "STATUS\t" << (complete ? "SCHEMA_SAMPLE_OK" : "NEEDS_REVIEW") << '\n';
  std::cout << "NOTE\tSchema/sample check only; no foil labels, geometry, cuts, or matrix validated.\n";
  return complete ? 0 : 3;
}
