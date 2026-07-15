#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>

#include "TFile.h"
#include "TTree.h"
#include "TCutG.h"
#include "TString.h"
#include "TROOT.h"
#include "TSystem.h"

double QuantileLocal(std::vector<double> v, double q)
{
  if (v.empty()) return -999.0;

  std::sort(v.begin(), v.end());

  const double pos = q * (v.size() - 1);
  const int i = static_cast<int>(pos);
  const double frac = pos - i;

  if (i + 1 < static_cast<int>(v.size())) {
    return v[i] * (1.0 - frac) + v[i + 1] * frac;
  }

  return v[i];
}


void run_xfp_xpfp_dynamicSplit_batch(
    Int_t metadataRun = 611704,
    TString rungroupTag = "rg04_theta12p395_foilpm3",
    TString campaignDir = "HMS_6p117GeV",
    TString inputRootFile = "",
    Int_t foilIndex = 0,
    Int_t maxBands = 9,
    Double_t thetaStepDeg = 1.0,
    Double_t minPeakSep = 0.12,
    Double_t minPeakFrac = 0.06,
    Double_t minPromFrac = 0.18,
    Int_t smoothPasses = 2,
    Double_t weakFrac = 0.003,
    Double_t weakProm = 0.025,
    Double_t weakSepFactor = 0.30,
    Bool_t writeAutoCuts = true,
    Double_t qlo = 0.05,
    Double_t qhi = 0.95,
    Int_t minEventsPerZone = 1000,
    Double_t minCerNpe = 6.0,
    Double_t minCalEtot = 0.65,
    Long64_t maxEvents = -1)
{
  gROOT->SetBatch(kTRUE);

  campaignDir = campaignDir.Strip(TString::kBoth);
  rungroupTag = rungroupTag.Strip(TString::kBoth);
  inputRootFile = inputRootFile.Strip(TString::kBoth);

  if (campaignDir.EndsWith("/")) campaignDir.Chop();

  if (campaignDir.Length() == 0) {
    std::cerr << "ERROR: campaignDir must be supplied." << std::endl;
    return;
  }

  if (rungroupTag.Length() == 0) {
    std::cerr << "ERROR: rungroupTag must be supplied." << std::endl;
    return;
  }

  if (inputRootFile.Length() == 0) {
    std::cerr << "ERROR: inputRootFile must be supplied explicitly." << std::endl;
    return;
  }

  const TString workerMacro = "assign_xfp_xpfp_angleScanBands_split.C";

  const TString stepDir  = campaignDir + "/02b_angle_scan_x";
  const TString cutsDir  = stepDir + "/cuts";
  const TString plotsDir = stepDir + "/plots";
  const TString rootDir  = stepDir + "/root";
  const TString tsvDir   = stepDir + "/tsv";
  const TString logsDir  = stepDir + "/logs";

  gSystem->mkdir(stepDir.Data(),  kTRUE);
  gSystem->mkdir(cutsDir.Data(),  kTRUE);
  gSystem->mkdir(plotsDir.Data(), kTRUE);
  gSystem->mkdir(rootDir.Data(),  kTRUE);
  gSystem->mkdir(tsvDir.Data(),   kTRUE);
  gSystem->mkdir(logsDir.Data(),  kTRUE);

  const TString ycutFile =
      Form("%s/01_ytar_cuts/cuts/ytar_ridge_cut_%s.root",
           campaignDir.Data(),
           rungroupTag.Data());

  const TString ycutName =
      Form("delta_vs_ytar_cut_foil%d", foilIndex);

  const TString splitTsv =
      Form("%s/xfp_dynamic_splits_%s_foil%d.tsv",
           tsvDir.Data(),
           rungroupTag.Data(),
           foilIndex);

  std::cout << "\n============================================================\n";
  std::cout << "Dynamic XFP/XPFP split batch\n";
  std::cout << "metadata run = " << metadataRun << "\n";
  std::cout << "rungroup     = " << rungroupTag << "\n";
  std::cout << "campaign     = " << campaignDir << "\n";
  std::cout << "foil         = " << foilIndex << "\n";
  std::cout << "split rule   = midpoint of xpfp q" << qlo
            << " and q" << qhi << "\n";
  std::cout << "ROOT file    = " << inputRootFile << "\n";
  std::cout << "ytar cut     = " << ycutFile << " :: " << ycutName << "\n";
  std::cout << "split TSV    = " << splitTsv << "\n";
  std::cout << "============================================================\n\n";

  TFile fin(inputRootFile, "READ");
  if (fin.IsZombie()) {
    std::cerr << "ERROR: cannot open " << inputRootFile << std::endl;
    return;
  }

  TTree* T = nullptr;
  fin.GetObject("T", T);
  if (!T) fin.GetObject("Tout", T);

  if (!T) {
    std::cerr << "ERROR: could not find tree T or Tout in "
              << inputRootFile << std::endl;
    return;
  }

  TFile fcut(ycutFile, "READ");
  if (fcut.IsZombie()) {
    std::cerr << "ERROR: cannot open " << ycutFile << std::endl;
    return;
  }

  TCutG* ycut = nullptr;
  fcut.GetObject(ycutName, ycut);

  if (!ycut) {
    std::cerr << "ERROR: cannot find " << ycutName
              << " in " << ycutFile << std::endl;
    return;
  }

  const double edges[] = {-10.0, -8.0, -5.0, 0.0, 5.0, 10.0};
  const int nEdges = sizeof(edges) / sizeof(edges[0]);
  const int nSlices = nEdges - 1;

  double xpfp = 0.0;
  double delta = 0.0;
  double ytar = 0.0;
  double cer = 0.0;
  double cal = 0.0;

  T->SetBranchStatus("*", 0);
  T->SetBranchStatus("H.dc.xp_fp", 1);
  T->SetBranchStatus("H.gtr.dp", 1);
  T->SetBranchStatus("H.gtr.y", 1);
  T->SetBranchStatus("H.cer.npeSum", 1);
  T->SetBranchStatus("H.cal.etottracknorm", 1);

  T->SetBranchAddress("H.dc.xp_fp", &xpfp);
  T->SetBranchAddress("H.gtr.dp", &delta);
  T->SetBranchAddress("H.gtr.y", &ytar);
  T->SetBranchAddress("H.cer.npeSum", &cer);
  T->SetBranchAddress("H.cal.etottracknorm", &cal);

  std::vector<std::vector<double> > xpfpBySlice(nSlices);

  Long64_t nentries = T->GetEntries();
  if (maxEvents > 0 && maxEvents < nentries) nentries = maxEvents;

  for (Long64_t i = 0; i < nentries; ++i) {
    T->GetEntry(i);

    if (!(cer > minCerNpe)) continue;
    if (!(cal > minCalEtot)) continue;
    if (!ycut->IsInside(ytar, delta)) continue;

    for (int is = 0; is < nSlices; ++is) {
      if (delta >= edges[is] && delta < edges[is + 1]) {
        xpfpBySlice[is].push_back(xpfp);
        break;
      }
    }
  }

  std::ofstream ofs(splitTsv.Data());
  if (!ofs.is_open()) {
    std::cerr << "ERROR: cannot create " << splitTsv << std::endl;
    return;
  }

  ofs << "rungroup\tmetadata_run\tfoil\tndel\tdelta_low\tdelta_high\t"
         "N\tqlo_value\tmedian\tqhi_value\tsplit\tN_low\tN_high\t"
         "ran_low\tran_high\n";

  gROOT->LoadMacro(workerMacro);

  for (int is = 0; is < nSlices; ++is) {
    const double dmin = edges[is];
    const double dmax = edges[is + 1];
    const int N = static_cast<int>(xpfpBySlice[is].size());

    std::cout << "\n============================================================\n";
    std::cout << "Dynamic split slice ndel=" << is
              << " delta=[" << dmin << "," << dmax << ")%"
              << " N=" << N << "\n";

    if (N < minEventsPerZone) {
      std::cout << "Skipping slice: too few total events.\n";

      ofs << rungroupTag << "\t"
          << metadataRun << "\t"
          << foilIndex << "\t"
          << is << "\t"
          << dmin << "\t"
          << dmax << "\t"
          << N << "\t"
          << -999.0 << "\t"
          << -999.0 << "\t"
          << -999.0 << "\t"
          << -999.0 << "\t"
          << 0 << "\t"
          << 0 << "\t"
          << 0 << "\t"
          << 0 << "\n";

      continue;
    }

    const double pLo  = QuantileLocal(xpfpBySlice[is], qlo);
    const double pMed = QuantileLocal(xpfpBySlice[is], 0.50);
    const double pHi  = QuantileLocal(xpfpBySlice[is], qhi);

    const double split = 0.5 * (pLo + pHi);

    int nLow = 0;
    int nHigh = 0;

    for (double v : xpfpBySlice[is]) {
      if (v < split) ++nLow;
      else ++nHigh;
    }

    const bool runLow = (nLow >= minEventsPerZone);
    const bool runHigh = (nHigh >= minEventsPerZone);

    std::cout << "xpfp q" << qlo << " = " << pLo << "\n";
    std::cout << "xpfp median = " << pMed << "\n";
    std::cout << "xpfp q" << qhi << " = " << pHi << "\n";
    std::cout << "dynamic split = " << split << "\n";
    std::cout << "N low/high = " << nLow << " / " << nHigh << "\n";
    std::cout << "============================================================\n";

    ofs << rungroupTag << "\t"
        << metadataRun << "\t"
        << foilIndex << "\t"
        << is << "\t"
        << dmin << "\t"
        << dmax << "\t"
        << N << "\t"
        << pLo << "\t"
        << pMed << "\t"
        << pHi << "\t"
        << split << "\t"
        << nLow << "\t"
        << nHigh << "\t"
        << (runLow ? 1 : 0) << "\t"
        << (runHigh ? 1 : 0) << "\n";

    if (runLow) {
      TString cmdLow = Form(
          "assign_xfp_xpfp_angleScanBands_split("
          "%d,%g,%g,\"%s\","
          "%d,%g,%g,%g,%g,%d,%g,%g,%g,"
          "true,%d,%lld,-1,%d,false,-999.0,-999.0,%g,%s,"
          "\"%s\",\"%s\")",
          metadataRun,
          dmin,
          dmax,
          rungroupTag.Data(),
          maxBands,
          thetaStepDeg,
          minPeakSep,
          minPeakFrac,
          minPromFrac,
          smoothPasses,
          weakFrac,
          weakProm,
          weakSepFactor,
          foilIndex,
          maxEvents,
          is,
          split,
          writeAutoCuts ? "true" : "false",
          campaignDir.Data(),
          inputRootFile.Data());

      std::cout << "\n--- LOW dynamic zone: xpfp < "
                << split << " ---\n";
      std::cout << cmdLow << "\n";

      gROOT->ProcessLine(cmdLow);
    } else {
      std::cout << "Skipping LOW zone: too few events after split.\n";
    }

    if (runHigh) {
      TString cmdHigh = Form(
          "assign_xfp_xpfp_angleScanBands_split("
          "%d,%g,%g,\"%s\","
          "%d,%g,%g,%g,%g,%d,%g,%g,%g,"
          "true,%d,%lld,-1,%d,false,-999.0,%g,999.0,%s,"
          "\"%s\",\"%s\")",
          metadataRun,
          dmin,
          dmax,
          rungroupTag.Data(),
          maxBands,
          thetaStepDeg,
          minPeakSep,
          minPeakFrac,
          minPromFrac,
          smoothPasses,
          weakFrac,
          weakProm,
          weakSepFactor,
          foilIndex,
          maxEvents,
          is,
          split,
          writeAutoCuts ? "true" : "false",
          campaignDir.Data(),
          inputRootFile.Data());

      std::cout << "\n--- HIGH dynamic zone: xpfp >= "
                << split << " ---\n";
      std::cout << cmdHigh << "\n";

      gROOT->ProcessLine(cmdHigh);
    } else {
      std::cout << "Skipping HIGH zone: too few events after split.\n";
    }
  }

  ofs.close();

  std::cout << "\nDone dynamic split batch for "
            << rungroupTag
            << ", foil " << foilIndex << ".\n";
  std::cout << "Wrote split summary: " << splitTsv << "\n";
}
