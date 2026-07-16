/*
Diagnostic for the retained SVD fit-sample distribution.

Run from the repository top level with only the campaign directory:

  root -l -b -q 'analyze_fit_sample_balance.C("HMS_6p667GeV")'

Optional arguments:

  analyze_fit_sample_balance(
    campaign,
    fileId,
    maxPerYBin,
    maxPerFoil,
    maxPerRunFoil,
    nfitMax
  )

Example with explicit limits:

  root -l -b -q 'analyze_fit_sample_balance.C("HMS_6p667GeV",-1,1000,15000,30000,200000)'

The macro reads the campaign rungroup TSV and unified numeric-ID TFit
trees automatically. It writes diagnostic TSV files under:

  <campaign>/06b_svd_fit/diagnostics/
*/

#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TSystem.h>

#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

using namespace std;

struct SettingInfo {
  TString rungroup;
  Int_t opticsId;
};

using YKey = tuple<Int_t, Int_t, Int_t>;          // foil, ndel, yscol
using XKey = tuple<Int_t, Int_t, Int_t, Int_t>;   // foil, ndel, yscol, xscol

void analyze_fit_sample_balance(
  TString campaign,
  Int_t fileId = -1,
  Int_t maxPerYBin = 1000,
  Int_t maxPerFoil = 15000,
  Int_t maxPerRunFoil = 30000,
  Long64_t nfitMax = 200000
) {
  TString campaignDir = campaign;
  gSystem->ExpandPathName(campaignDir);

  if (gSystem->AccessPathName(campaignDir)) {
    cout << "ERROR: campaign directory does not exist: "
         << campaignDir << endl;
    return;
  }

  TString configDir = Form("%s/config", campaignDir.Data());
  TString inputTreeDir =
    Form("%s/06a_fit_ntuple/root", campaignDir.Data());
  TString outputDir =
    Form("%s/06b_svd_fit/diagnostics", campaignDir.Data());

  TString findCommand = Form(
    "find '%s' -maxdepth 1 -type f "
    "-name 'rungroups_*_inputs.tsv' -print",
    configDir.Data()
  );

  TString rungroupsTsv = gSystem->GetFromPipe(findCommand);
  rungroupsTsv = rungroupsTsv.Strip(TString::kBoth);

  if (rungroupsTsv.IsNull()) {
    cout << "ERROR: no rungroups_*_inputs.tsv found in "
         << configDir << endl;
    return;
  }

  if (rungroupsTsv.Contains("\n")) {
    cout << "ERROR: multiple rungroups_*_inputs.tsv files found in "
         << configDir << ":" << endl;
    cout << rungroupsTsv << endl;
    return;
  }

  if (gSystem->AccessPathName(inputTreeDir)) {
    cout << "ERROR: fit-tree directory does not exist: "
         << inputTreeDir << endl;
    return;
  }

  cout << "Campaign       : " << campaignDir << endl;
  cout << "Rungroups TSV  : " << rungroupsTsv << endl;
  cout << "Fit-tree dir   : " << inputTreeDir << endl;
  cout << "Diagnostic dir : " << outputDir << endl;

  vector<SettingInfo> settings;

  ifstream in(rungroupsTsv.Data());
  if (!in.is_open()) {
    cout << "ERROR: cannot open " << rungroupsTsv << endl;
    return;
  }

  string line;
  while (getline(in, line)) {
    if (line.empty()) continue;
    if (line.rfind("rungroup\t", 0) == 0) continue;

    string rungroupField;
    string opticsIdField;
    stringstream row(line);

    if (!getline(row, rungroupField, '\t') ||
        !getline(row, opticsIdField, '\t')) {
      cout << "ERROR: malformed TSV row: " << line << endl;
      return;
    }

    SettingInfo setting;
    setting.rungroup = rungroupField.c_str();
    setting.opticsId = TString(opticsIdField.c_str()).Atoi();

    if (setting.rungroup.IsNull() || setting.opticsId <= 0) {
      cout << "ERROR: invalid TSV row: " << line << endl;
      return;
    }

    settings.push_back(setting);
  }
  in.close();

  if (settings.empty()) {
    cout << "ERROR: no settings read from " << rungroupsTsv << endl;
    return;
  }

  gSystem->mkdir(outputDir, kTRUE);

  TString settingOut =
    Form("%s/fit_sample_setting_summary.tsv", outputDir.Data());
  TString ycellOut =
    Form("%s/fit_sample_ycell_counts.tsv", outputDir.Data());
  TString xcellOut =
    Form("%s/fit_sample_xcell_counts.tsv", outputDir.Data());

  ofstream settingFile(settingOut.Data());
  ofstream ycellFile(ycellOut.Data());
  ofstream xcellFile(xcellOut.Data());

  if (!settingFile || !ycellFile || !xcellFile) {
    cout << "ERROR: cannot create diagnostic output files in "
         << outputDir << endl;
    return;
  }

  settingFile
    << "setting_index\trungroup\toptics_id\ttree_entries"
    << "\tavailable_ycells\tavailable_xcells"
    << "\tretained_events\tglobal_before\tglobal_after"
    << "\tglobal_cap_reached\n";

  ycellFile
    << "setting_index\trungroup\toptics_id"
    << "\tfoil\tndel\tyscol"
    << "\tavailable\tretained"
    << "\tavailable_xcols\tretained_xcols"
    << "\tdominant_x_available_fraction"
    << "\tdominant_x_retained_fraction\n";

  xcellFile
    << "setting_index\trungroup\toptics_id"
    << "\tfoil\tndel\tyscol\txscol"
    << "\tavailable\tretained\n";

  Long64_t globalRetained = 0;

  cout << "Current-cap simulation" << endl;
  cout << "  Max per foil/ndel/yscol = " << maxPerYBin << endl;
  cout << "  Max per run/foil        = " << maxPerFoil << endl;
  cout << "  Secondary run/foil max  = " << maxPerRunFoil << endl;
  cout << "  Global NFIT_MAX         = " << nfitMax << endl;

  for (size_t iSetting = 0; iSetting < settings.size(); ++iSetting) {
    const SettingInfo &setting = settings[iSetting];

    TString inputFile = Form(
      "%s/Optics_%d_%d_fit_tree_gmm.root",
      inputTreeDir.Data(),
      setting.opticsId,
      fileId
    );

    TFile *f = TFile::Open(inputFile, "READ");
    if (!f || f->IsZombie()) {
      cout << "ERROR: cannot open " << inputFile << endl;
      if (f) {
        f->Close();
        delete f;
      }
      continue;
    }

    TTree *tree = dynamic_cast<TTree *>(f->Get("TFit"));
    if (!tree) {
      cout << "ERROR: no TFit tree in " << inputFile << endl;
      f->Close();
      delete f;
      continue;
    }

    Int_t foil = -1;
    Int_t ndel = -1;
    Int_t yscol = -1;
    Int_t xscol = -1;

    tree->SetBranchAddress("foil", &foil);
    tree->SetBranchAddress("ndel", &ndel);
    tree->SetBranchAddress("yscol", &yscol);
    tree->SetBranchAddress("xscol", &xscol);

    map<YKey, Long64_t> availableY;
    map<YKey, Long64_t> retainedY;
    map<XKey, Long64_t> availableX;
    map<XKey, Long64_t> retainedX;
    map<Int_t, Long64_t> retainedFoil;

    const Long64_t globalBefore = globalRetained;
    const Long64_t nentries = tree->GetEntries();

    for (Long64_t i = 0; i < nentries; ++i) {
      tree->GetEntry(i);

      YKey ykey = make_tuple(foil, ndel, yscol);
      XKey xkey = make_tuple(foil, ndel, yscol, xscol);

      availableY[ykey]++;
      availableX[xkey]++;

      const bool keep =
        globalRetained < nfitMax &&
        retainedY[ykey] < maxPerYBin &&
        retainedFoil[foil] < maxPerFoil &&
        retainedFoil[foil] < maxPerRunFoil;

      if (keep) {
        retainedY[ykey]++;
        retainedX[xkey]++;
        retainedFoil[foil]++;
        globalRetained++;
      }
    }

    Long64_t settingRetained = globalRetained - globalBefore;
    const bool globalCapReached =
      globalBefore < nfitMax && globalRetained >= nfitMax;

    settingFile
      << iSetting << "\t"
      << setting.rungroup << "\t"
      << setting.opticsId << "\t"
      << nentries << "\t"
      << availableY.size() << "\t"
      << availableX.size() << "\t"
      << settingRetained << "\t"
      << globalBefore << "\t"
      << globalRetained << "\t"
      << (globalCapReached ? 1 : 0)
      << "\n";

    for (const auto &item : availableY) {
      const YKey &ykey = item.first;

      Int_t keyFoil;
      Int_t keyNdel;
      Int_t keyYscol;
      tie(keyFoil, keyNdel, keyYscol) = ykey;

      const Long64_t nAvailable = item.second;
      const Long64_t nRetained = retainedY[ykey];

      set<Int_t> availableXcols;
      set<Int_t> retainedXcols;
      Long64_t dominantAvailable = 0;
      Long64_t dominantRetained = 0;

      for (const auto &xitem : availableX) {
        Int_t xf, xd, xy, xx;
        tie(xf, xd, xy, xx) = xitem.first;

        if (xf != keyFoil || xd != keyNdel || xy != keyYscol) continue;

        availableXcols.insert(xx);
        if (xitem.second > dominantAvailable) {
          dominantAvailable = xitem.second;
        }

        Long64_t xr = retainedX[xitem.first];
        if (xr > 0) retainedXcols.insert(xx);
        if (xr > dominantRetained) dominantRetained = xr;
      }

      const double dominantAvailableFraction =
        nAvailable > 0
          ? static_cast<double>(dominantAvailable) / nAvailable
          : 0.0;

      const double dominantRetainedFraction =
        nRetained > 0
          ? static_cast<double>(dominantRetained) / nRetained
          : 0.0;

      ycellFile
        << iSetting << "\t"
        << setting.rungroup << "\t"
        << setting.opticsId << "\t"
        << keyFoil << "\t"
        << keyNdel << "\t"
        << keyYscol << "\t"
        << nAvailable << "\t"
        << nRetained << "\t"
        << availableXcols.size() << "\t"
        << retainedXcols.size() << "\t"
        << dominantAvailableFraction << "\t"
        << dominantRetainedFraction
        << "\n";
    }

    for (const auto &item : availableX) {
      Int_t keyFoil;
      Int_t keyNdel;
      Int_t keyYscol;
      Int_t keyXscol;
      tie(keyFoil, keyNdel, keyYscol, keyXscol) = item.first;

      xcellFile
        << iSetting << "\t"
        << setting.rungroup << "\t"
        << setting.opticsId << "\t"
        << keyFoil << "\t"
        << keyNdel << "\t"
        << keyYscol << "\t"
        << keyXscol << "\t"
        << item.second << "\t"
        << retainedX[item.first]
        << "\n";
    }

    cout
      << "Setting " << iSetting
      << " " << setting.rungroup
      << " entries=" << nentries
      << " retained=" << settingRetained
      << " global=" << globalRetained
      << endl;

    f->Close();
    delete f;
  }

  settingFile.close();
  ycellFile.close();
  xcellFile.close();

  cout << endl;
  cout << "Final retained total: " << globalRetained << endl;
  cout << "Wrote:" << endl;
  cout << "  " << settingOut << endl;
  cout << "  " << ycellOut << endl;
  cout << "  " << xcellOut << endl;
}
