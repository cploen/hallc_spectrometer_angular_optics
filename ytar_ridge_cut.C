// ytar_ridge_cut.C //
// Multi-foil extension of central-width-guarded ridge-envelope selection.
//
// Method:
//   1. Apply base production cuts.
//   2. In the central/high-density delta region, find dominant ytar ridges.
//   3. Sort ridges by ytar.
//      foil 0 = most negative ytar ridge.
//      foil 1 = next ridge.
//      etc.
//   4. For each detected ridge independently:
//        - track ridge center vs delta
//        - measure left/right width from local density falloff
//        - calibrate max width from that ridge's central/high-density region
//        - cap local widths using that ridge-specific central width
//        - build one TCutG per ridge
//   5. Overlay expert/human cuts if available for benchmarking.
//
// Important:
//   Data-driven parameters are computed separately per ridge.
//   The hand cuts are overlay/benchmark only. They do not control auto cuts.
//
// Output:
//   cuts/ytar_delta_<run>_<tag>_multifoil_cut.root
//     delta_vs_ytar_cut_foil0
//     delta_vs_ytar_cut_foil1
//     ...
//
// Usage:
//   .x ytar_ridge_cut.C(1544,"runtime_tag","-1")

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TCutG.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TLatex.h"
#include "TKey.h"

using namespace std;

struct EventLite {
  double ytar;
  double delta;
};

struct PeakInfo {
  int index = -1;
  int bin = -1;
  double ytar = 0.0;
  double height = 0.0;
  double fracOfMax = 0.0;
  double searchHalfWidth = 0.75;
};

struct CentralWidthCalibration {
  double refLeft = 0.0;
  double refRight = 0.0;
  double maxLeft = 0.0;
  double maxRight = 0.0;
  int nUsed = 0;
};

struct RidgeCutResult {
  int foilIndex = -1;
  PeakInfo peak;
  CentralWidthCalibration calib;

  vector<double> deltaCenters;
  vector<int> rowNVec;
  vector<double> ridgeCenterRaw;
  vector<double> ridgeCenter;
  vector<double> peakValVec;
  vector<double> widthLeftRaw;
  vector<double> widthRightRaw;
  vector<double> widthLeftFinal;
  vector<double> widthRightFinal;
  vector<double> leftBoundary;
  vector<double> rightBoundary;

  Long64_t nAuto = 0;
  Long64_t nHand = 0;
  Long64_t nBoth = 0;

  TCutG *cut = nullptr;
  TGraph *gLeft = nullptr;
  TGraph *gRight = nullptr;
  TGraph *gCenter = nullptr;
  TGraph *gCenterRaw = nullptr;
  TGraph *gWidthLeft = nullptr;
  TGraph *gWidthRight = nullptr;
};

bool HasBranch(TTree *T, const char *bname)
{
  return T && T->GetBranch(bname);
}

TString BuildInputRootPath(Int_t nrun, TString inputFileID)
{
  vector<TString> candidates;

  candidates.push_back(Form("ROOTfiles/OPTICS/newfit_6p667_20260526_1226_no_offsets/rootfiles/nps_hms_optics_%d_1_%s.root", nrun, inputFileID.Data()));
  for (auto &p : candidates) {
    if (!gSystem->AccessPathName(p)) {
      cout << "Found input ROOT file: " << p << endl;
      return p;
    }
  }

  cerr << "ERROR: Could not auto-find input ROOT file. Tried:" << endl;
  for (auto &p : candidates) cerr << "  " << p << endl;

  return candidates[0];
}

double MeanCutX(TCutG *cut)
{
  if (!cut || cut->GetN() <= 0) return 0.0;

  double sum = 0.0;
  for (int i = 0; i < cut->GetN(); i++) {
    double x, y;
    cut->GetPoint(i, x, y);
    sum += x;
  }

  return sum / cut->GetN();
}

vector<TCutG*> LoadExpertCuts(Int_t nrun, TString handFileID)
{
  vector<TCutG*> cuts;

  TString handCutFile = Form("01_ytar_cuts/cuts/ytar_delta_%d_%s_cut.root", nrun, handFileID.Data());

  if (gSystem->AccessPathName(handCutFile)) {
    cerr << "WARNING: expert cut file not found: " << handCutFile << endl;
    return cuts;
  }

  TFile *fcut = TFile::Open(handCutFile, "READ");
  if (!fcut || fcut->IsZombie()) {
    cerr << "WARNING: could not open expert cut file: " << handCutFile << endl;
    return cuts;
  }

  // First try standard names.
  for (int i = 0; i < 20; i++) {
    TString cname = Form("delta_vs_ytar_cut_foil%d", i);
    TCutG *c = (TCutG*)fcut->Get(cname);
    if (c) {
      TCutG *clone = (TCutG*)c->Clone(Form("expert_%s", cname.Data()));
      clone->SetLineColor(kRed);
      clone->SetLineWidth(3);
      cuts.push_back(clone);
    }
  }

  // If standard names failed, collect all TCutG objects.
  if (cuts.empty()) {
    TIter next(fcut->GetListOfKeys());
    TKey *key = nullptr;

    while ((key = (TKey*)next())) {
      TObject *obj = key->ReadObj();
      if (!obj) continue;

      if (obj->InheritsFrom("TCutG")) {
        TCutG *clone = (TCutG*)obj->Clone(Form("expert_%s", obj->GetName()));
        clone->SetLineColor(kRed);
        clone->SetLineWidth(3);
        cuts.push_back(clone);
      }
    }
  }

  sort(cuts.begin(), cuts.end(),
       [](TCutG *a, TCutG *b) {
         return MeanCutX(a) < MeanCutX(b);
       });

  for (size_t i = 0; i < cuts.size(); i++) {
    cuts[i]->SetLineColor(kRed);
    cuts[i]->SetLineWidth(3);
  }

  fcut->Close();

  cout << "Loaded expert cuts: " << cuts.size() << endl;
  return cuts;
}

void SmoothVectorInPlace(vector<double> &v, int passes=2)
{
  if (v.size() < 5) return;

  for (int p = 0; p < passes; p++) {
    vector<double> old = v;

    for (size_t i = 1; i + 1 < v.size(); i++) {
      v[i] = 0.25 * old[i - 1] + 0.50 * old[i] + 0.25 * old[i + 1];
    }
  }
}

double Median(vector<double> v)
{
  if (v.empty()) return 0.0;
  sort(v.begin(), v.end());
  return v[v.size() / 2];
}

double MaxValue(const vector<double> &v)
{
  if (v.empty()) return 0.0;

  double m = v[0];
  for (auto x : v) {
    if (x > m) m = x;
  }

  return m;
}

CentralWidthCalibration GetCentralWidthCalibration(
  const vector<double> &deltaCenters,
  const vector<double> &peakValVec,
  const vector<double> &widthLeftRaw,
  const vector<double> &widthRightRaw,
  double centralDeltaAbsMax,
  double minPeakFracOfGlobal,
  double guardScale
)
{
  CentralWidthCalibration out;

  double globalPeak = MaxValue(peakValVec);

  vector<double> leftVals;
  vector<double> rightVals;

  for (size_t i = 0; i < deltaCenters.size(); i++) {
    if (std::abs(deltaCenters[i]) > centralDeltaAbsMax) continue;
    if (globalPeak > 0.0 && peakValVec[i] < minPeakFracOfGlobal * globalPeak) continue;

    if (widthLeftRaw[i] > 0.0 && widthRightRaw[i] > 0.0) {
      leftVals.push_back(widthLeftRaw[i]);
      rightVals.push_back(widthRightRaw[i]);
    }
  }

  // Fallback: if central/high-density criterion is too strict,
  // use |delta| < centralDeltaAbsMax.
  if (leftVals.size() < 3) {
    leftVals.clear();
    rightVals.clear();

    for (size_t i = 0; i < deltaCenters.size(); i++) {
      if (std::abs(deltaCenters[i]) > centralDeltaAbsMax) continue;

      if (widthLeftRaw[i] > 0.0 && widthRightRaw[i] > 0.0) {
        leftVals.push_back(widthLeftRaw[i]);
        rightVals.push_back(widthRightRaw[i]);
      }
    }
  }

  out.refLeft = Median(leftVals);
  out.refRight = Median(rightVals);
  out.maxLeft = guardScale * out.refLeft;
  out.maxRight = guardScale * out.refRight;
  out.nUsed = leftVals.size();

  return out;
}

double FindWidthAtFractionSide(TH1D *h,
                               int peakBin,
                               int direction,
                               double peakVal,
                               double fracLevel,
                               double fallbackWidth)
{
  if (!h || peakBin < 1 || peakBin > h->GetNbinsX()) return fallbackWidth;
  if (direction != -1 && direction != +1) return fallbackWidth;
  if (peakVal <= 0.0) return fallbackWidth;

  int nb = h->GetNbinsX();
  double peakY = h->GetBinCenter(peakBin);
  double target = fracLevel * peakVal;

  for (int step = 1; step < nb; step++) {
    int b = peakBin + direction * step;
    if (b <= 1 || b >= nb) break;

    double val = h->GetBinContent(b);

    if (val <= target) {
      double y = h->GetBinCenter(b);
      return std::abs(y - peakY);
    }
  }

  return fallbackWidth;
}

vector<PeakInfo> FindMajorPeaks(TH1D *hSmooth,
                                double minPeakFracOfMax,
                                double minPeakSeparation)
{
  vector<PeakInfo> rawPeaks;
  vector<PeakInfo> finalPeaks;

  if (!hSmooth) return finalPeaks;

  int nb = hSmooth->GetNbinsX();
  double maxVal = hSmooth->GetMaximum();

  if (maxVal <= 0.0) return finalPeaks;

  for (int b = 2; b <= nb - 1; b++) {
    double y = hSmooth->GetBinCenter(b);
    double v = hSmooth->GetBinContent(b);
    double vl = hSmooth->GetBinContent(b - 1);
    double vr = hSmooth->GetBinContent(b + 1);

    if (v < minPeakFracOfMax * maxVal) continue;

    if (v >= vl && v >= vr) {
      PeakInfo p;
      p.bin = b;
      p.ytar = y;
      p.height = v;
      p.fracOfMax = v / maxVal;
      rawPeaks.push_back(p);
    }
  }

  // Keep strongest peaks first, suppress close duplicates.
  sort(rawPeaks.begin(), rawPeaks.end(),
       [](const PeakInfo &a, const PeakInfo &b) {
         return a.height > b.height;
       });

  for (auto &p : rawPeaks) {
    bool tooClose = false;

    for (auto &q : finalPeaks) {
      if (std::abs(p.ytar - q.ytar) < minPeakSeparation) {
        tooClose = true;
        break;
      }
    }

    if (!tooClose) finalPeaks.push_back(p);
  }

  // Final index convention: increasing ytar.
  sort(finalPeaks.begin(), finalPeaks.end(),
       [](const PeakInfo &a, const PeakInfo &b) {
         return a.ytar < b.ytar;
       });

  for (size_t i = 0; i < finalPeaks.size(); i++) {
    finalPeaks[i].index = (int)i;
  }

  // Data-driven local search half-width per ridge.
  // For multiple ridges, use a fraction of nearest ridge separation.
  // For one ridge, use a broad but finite window.
  for (size_t i = 0; i < finalPeaks.size(); i++) {
    double nearestSep = 999.0;

    if (finalPeaks.size() == 1) {
      finalPeaks[i].searchHalfWidth = 0.75;
      continue;
    }

    for (size_t j = 0; j < finalPeaks.size(); j++) {
      if (i == j) continue;
      double sep = std::abs(finalPeaks[i].ytar - finalPeaks[j].ytar);
      if (sep < nearestSep) nearestSep = sep;
    }

    double hw = 0.45 * nearestSep;

    if (hw < 0.18) hw = 0.18;
    if (hw > 0.75) hw = 0.75;

    finalPeaks[i].searchHalfWidth = hw;
  }

  return finalPeaks;
}

int FindPeakBinNearExpected(TH1D *h,
                            double expectedY,
                            double searchHalfWidth)
{
  if (!h) return -1;

  int bestBin = -1;
  double bestVal = -1.0;

  for (int b = 1; b <= h->GetNbinsX(); b++) {
    double y = h->GetBinCenter(b);

    if (std::abs(y - expectedY) > searchHalfWidth) continue;

    double v = h->GetBinContent(b);

    if (v > bestVal) {
      bestVal = v;
      bestBin = b;
    }
  }

  if (bestBin < 0) bestBin = h->GetMaximumBin();

  return bestBin;
}

int FindDeltaIndexNear(const vector<double> &deltaCenters,
                       double delta,
                       double tol=1.0e-6)
{
  for (size_t i = 0; i < deltaCenters.size(); i++) {
    if (std::abs(deltaCenters[i] - delta) < tol) return (int)i;
  }

  return -1;
}

int FindPeakBinNearExpectedGuarded(TH1D *h,
                                   double expectedY,
                                   double searchHalfWidth,
                                   double maxCenterStep,
                                   const vector<double> &excludedCenters,
                                   double minExcludedGap)
{
  if (!h) return -1;

  const int nb = h->GetNbinsX();
  if (nb <= 0) return -1;

  if (maxCenterStep <= 0.0) maxCenterStep = searchHalfWidth;

  const double allowedHalfWidth = std::min(searchHalfWidth, maxCenterStep);
  const double hMax = h->GetMaximum();

  auto isExcluded = [&](double y) {
    for (auto yc : excludedCenters) {
      if (std::abs(y - yc) < minExcludedGap) return true;
    }
    return false;
  };

  // First pass: local maxima only.
  // Second pass: any populated bin, still respecting identity guards.
  for (int pass = 0; pass < 2; pass++) {
    int bestBin = -1;
    double bestScore = -1.0e99;

    for (int b = 2; b <= nb - 1; b++) {
      const double y = h->GetBinCenter(b);
      const double dist = std::abs(y - expectedY);

      if (dist > allowedHalfWidth) continue;
      if (isExcluded(y)) continue;

      const double v  = h->GetBinContent(b);
      const double vl = h->GetBinContent(b - 1);
      const double vr = h->GetBinContent(b + 1);

      if (v <= 0.0) continue;

      if (pass == 0 && !(v >= vl && v >= vr)) continue;

      const double densityScore = (hMax > 0.0) ? v / hMax : v;
      const double distanceScore = (allowedHalfWidth > 0.0) ? dist / allowedHalfWidth : 0.0;

      // Prefer strong local peaks, but penalize large jumps.
      // This prevents a weak/missing foil from stealing a neighboring dominant ridge.
      const double score = densityScore - 0.35 * distanceScore;

      if (score > bestScore) {
        bestScore = score;
        bestBin = b;
      }
    }

    if (bestBin >= 1) return bestBin;
  }

  // No guarded candidate found.  Deliberately skip this delta row rather
  // than falling back to the global maximum and risking a foil identity swap.
  return -1;
}

RidgeCutResult BuildCutForRidge(const vector<EventLite> &events,
                                const PeakInfo &peak,
                                int nrun,
                                TString tag,
                                const vector<TCutG*> &expertCuts,
                                double deltaMin,
                                double deltaMax,
                                double deltaStep,
                                double deltaHalfWindow,
                                double widthFracLevel,
                                int smoothPasses1D,
                                int smoothPassesEnvelope,
                                double centralDeltaAbsMax,
                                double minPeakFracOfGlobal,
                                double guardScale,
                                double weakRidgeFracThreshold,
                                double strongRidgeMinWidthScale,
                                double weakRidgeMinWidthScale,
                                double weakRidgeGuardBoost,
                                double rawWideningRetain,
                                int minRowEntries,
                                const vector<RidgeCutResult> *previousResults,
                                double minInterFoilTrackGap,
                                double maxCenterStepFrac)
{
  RidgeCutResult result;
  result.foilIndex = peak.index;
  result.peak = peak;

  bool havePreviousCenter = true;
  double previousCenter = peak.ytar;

  for (double d = deltaMin + 0.5 * deltaStep;
       d <= deltaMax - 0.5 * deltaStep;
       d += deltaStep) {

    double histMin = previousCenter - peak.searchHalfWidth;
    double histMax = previousCenter + peak.searchHalfWidth;

    TH1D *hY = new TH1D(
      Form("hY_foil%d_%g", peak.index, d),
      Form("foil %d delta %.2f;ytar;counts", peak.index, d),
      180, histMin, histMax
    );

    int nrow = 0;

    for (const auto &ev : events) {
      if (std::abs(ev.delta - d) > deltaHalfWindow) continue;
      if (std::abs(ev.ytar - previousCenter) > peak.searchHalfWidth) continue;

      hY->Fill(ev.ytar);
      nrow++;
    }

    if (nrow < minRowEntries) {
      delete hY;
      continue;
    }

    for (int s = 0; s < smoothPasses1D; s++) {
      hY->Smooth(1);
    }

    vector<double> excludedCenters;

    if (previousResults) {
      for (const auto &prev : *previousResults) {
        int j = FindDeltaIndexNear(prev.deltaCenters, d);
        if (j >= 0 && j < (int)prev.ridgeCenterRaw.size()) {
          excludedCenters.push_back(prev.ridgeCenterRaw[j]);
        }
      }
    }

    double maxCenterStep = maxCenterStepFrac * peak.searchHalfWidth;
    if (maxCenterStep < 0.18) maxCenterStep = 0.18;

    int peakBin = FindPeakBinNearExpectedGuarded(
      hY,
      previousCenter,
      peak.searchHalfWidth,
      maxCenterStep,
      excludedCenters,
      minInterFoilTrackGap
    );

    if (peakBin < 1) {
      delete hY;
      continue;
    }

    double peakY = hY->GetBinCenter(peakBin);
    double peakVal = hY->GetBinContent(peakBin);

    if (peakVal <= 0.0) {
      delete hY;
      continue;
    }

    double fallbackWidth = 0.25 * peak.searchHalfWidth;

    double wL = FindWidthAtFractionSide(hY, peakBin, -1, peakVal,
                                        widthFracLevel, fallbackWidth);

    double wR = FindWidthAtFractionSide(hY, peakBin, +1, peakVal,
                                        widthFracLevel, fallbackWidth);

    result.deltaCenters.push_back(d);
    result.rowNVec.push_back(nrow);
    result.ridgeCenterRaw.push_back(peakY);
    result.peakValVec.push_back(peakVal);
    result.widthLeftRaw.push_back(wL);
    result.widthRightRaw.push_back(wR);

    previousCenter = peakY;
    havePreviousCenter = true;

    delete hY;
  }

  if (result.deltaCenters.size() < 4) {
    cerr << "WARNING: too few ridge points for foil "
         << peak.index << endl;
    return result;
  }

  result.ridgeCenter = result.ridgeCenterRaw;
  result.widthLeftFinal = result.widthLeftRaw;
  result.widthRightFinal = result.widthRightRaw;

  SmoothVectorInPlace(result.ridgeCenter, smoothPassesEnvelope);
  SmoothVectorInPlace(result.widthLeftFinal, smoothPassesEnvelope);
  SmoothVectorInPlace(result.widthRightFinal, smoothPassesEnvelope);

  result.calib = GetCentralWidthCalibration(
    result.deltaCenters,
    result.peakValVec,
    result.widthLeftRaw,
    result.widthRightRaw,
    centralDeltaAbsMax,
    minPeakFracOfGlobal,
    guardScale
  );

  // Conservative envelope guard:
  // If a ridge is statistically weaker, local density falloff can make the
  // measured width artificially shrink.  Do not let the envelope collapse.
  bool weakRidge = (peak.fracOfMax < weakRidgeFracThreshold);
  double minWidthScale = weakRidge ? weakRidgeMinWidthScale : strongRidgeMinWidthScale;

  if (weakRidge) {
    result.calib.maxLeft  *= weakRidgeGuardBoost;
    result.calib.maxRight *= weakRidgeGuardBoost;
  }

  cout << "Foil " << peak.index << " central width calibration:" << endl;
  cout << "  peak ytar = " << peak.ytar << endl;
  cout << "  searchHalfWidth = " << peak.searchHalfWidth << endl;
  cout << "  nUsed    = " << result.calib.nUsed << endl;
  cout << "  refLeft  = " << result.calib.refLeft << endl;
  cout << "  refRight = " << result.calib.refRight << endl;
  cout << "  maxLeft  = " << result.calib.maxLeft << endl;
  cout << "  maxRight = " << result.calib.maxRight << endl;

  if (result.calib.nUsed < 3 ||
      result.calib.refLeft <= 0.0 ||
      result.calib.refRight <= 0.0) {
    cerr << "WARNING: central width calibration failed for foil "
         << peak.index << endl;
    return result;
  }

  double minLeft  = minWidthScale * result.calib.refLeft;
  double minRight = minWidthScale * result.calib.refRight;

  for (size_t i = 0; i < result.widthLeftFinal.size(); i++) {
    // Preserve real local widening even after smoothing, but still cap wild tails.
    result.widthLeftFinal[i] = std::max(result.widthLeftFinal[i], rawWideningRetain * result.widthLeftRaw[i]);
    result.widthLeftFinal[i] = std::max(result.widthLeftFinal[i], minLeft);
    if (result.widthLeftFinal[i] > result.calib.maxLeft) result.widthLeftFinal[i] = result.calib.maxLeft;
  }

  for (size_t i = 0; i < result.widthRightFinal.size(); i++) {
    result.widthRightFinal[i] = std::max(result.widthRightFinal[i], rawWideningRetain * result.widthRightRaw[i]);
    result.widthRightFinal[i] = std::max(result.widthRightFinal[i], minRight);
    if (result.widthRightFinal[i] > result.calib.maxRight) result.widthRightFinal[i] = result.calib.maxRight;
  }

  for (size_t i = 0; i < result.deltaCenters.size(); i++) {
    double left = result.ridgeCenter[i] - result.widthLeftFinal[i];
    double right = result.ridgeCenter[i] + result.widthRightFinal[i];

    result.leftBoundary.push_back(left);
    result.rightBoundary.push_back(right);
  }

  // Cosmetic/stability smoothing after final width guards.
  // This removes one-bin boundary jitter without changing the broad envelope philosophy.
  SmoothVectorInPlace(result.leftBoundary, 2);
  SmoothVectorInPlace(result.rightBoundary, 2);

  int npts = 2 * result.deltaCenters.size() + 1;
  result.cut = new TCutG(Form("delta_vs_ytar_cut_foil%d", peak.index), npts);
  result.cut->SetTitle(Form("auto cut foil %d;ytar;delta", peak.index));

  int ip = 0;

  for (size_t i = 0; i < result.deltaCenters.size(); i++) {
    result.cut->SetPoint(ip++, result.leftBoundary[i], result.deltaCenters[i]);
  }

  for (int i = (int)result.deltaCenters.size() - 1; i >= 0; i--) {
    result.cut->SetPoint(ip++, result.rightBoundary[i], result.deltaCenters[i]);
  }

  result.cut->SetPoint(ip++, result.leftBoundary[0], result.deltaCenters[0]);

  int color = kBlue + peak.index;
  if (peak.index == 0) color = kBlue;
  if (peak.index == 1) color = kGreen+2;
  if (peak.index == 2) color = kMagenta+2;

  result.cut->SetLineColor(color);
  result.cut->SetLineWidth(3);

  // Counts against matching expert cut by index, if present.
  TCutG *expert = nullptr;
  if (peak.index >= 0 && peak.index < (int)expertCuts.size()) {
    expert = expertCuts[peak.index];
  }

  for (const auto &ev : events) {
    bool inAuto = result.cut->IsInside(ev.ytar, ev.delta);
    bool inHand = expert ? expert->IsInside(ev.ytar, ev.delta) : false;

    if (inAuto) result.nAuto++;
    if (inHand) result.nHand++;
    if (inAuto && inHand) result.nBoth++;
  }

  result.gLeft = new TGraph();
  result.gRight = new TGraph();
  result.gCenter = new TGraph();
  result.gCenterRaw = new TGraph();
  result.gWidthLeft = new TGraph();
  result.gWidthRight = new TGraph();

  for (size_t i = 0; i < result.deltaCenters.size(); i++) {
    result.gLeft->SetPoint(i, result.deltaCenters[i], result.leftBoundary[i]);
    result.gRight->SetPoint(i, result.deltaCenters[i], result.rightBoundary[i]);
    result.gCenter->SetPoint(i, result.deltaCenters[i], result.ridgeCenter[i]);
    result.gCenterRaw->SetPoint(i, result.deltaCenters[i], result.ridgeCenterRaw[i]);
    result.gWidthLeft->SetPoint(i, result.deltaCenters[i], result.widthLeftFinal[i]);
    result.gWidthRight->SetPoint(i, result.deltaCenters[i], result.widthRightFinal[i]);
  }

  return result;
}

void ytar_ridge_cut(Int_t nrun=1544,
                                             const char *tag="",
                                             const char *inputFileID="-1",
                                             const char *handFileID="-1",
                                             const char *inputRootOverride="",
                                             const char *campaignDir="HMS_6p117GeV")
{
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1,0);

  TString outTag = tag;
  TString inputID = inputFileID;
  TString handID = handFileID;
  TString inputRootOverrideStr = inputRootOverride;
  TString campaign = campaignDir;
  if (campaign.EndsWith("/")) campaign.Chop();
  TString step01 = campaign + "/01_ytar_cuts";

  // ----------------------------
  // Main settings
  // ----------------------------
  const double deltaMin = -10.0;
  const double deltaMax =  10.0;

  const double deltaStep = 0.25;
  const double deltaHalfWindow = 0.35;

  const double centralDeltaAbsMax = 1.5;

  const double plotYtarMin = -15.0;
  const double plotYtarMax =  15.0;

  const double peakFindYtarMin = -15.0;
  const double peakFindYtarMax =  15.0;
  const int peakFindNYBins = 600;

  const int centralSmoothPasses = 3;

  // Hall-C-specific dominant-ridge criteria.
  const double minPeakFracOfMax = 0.20;
  const double minPeakSeparation = 0.20;

  // Ridge-envelope settings.
  const double widthFracLevel = 0.25;
  const int smoothPasses1D = 2;
  const int smoothPassesEnvelope = 8;

  const double minPeakFracOfGlobal = 0.75;
  const double guardScale = 1.75;

  // Width guardrail: ridge finder should be a generous acceptance envelope.
  // GMM cleanup happens later, so avoid clipping real ridge events here.
  const double weakRidgeFracThreshold = 0.60;
  const double strongRidgeMinWidthScale = 1.10;
  const double weakRidgeMinWidthScale = 1.35;
  const double weakRidgeGuardBoost = 1.25;
  const double rawWideningRetain = 0.85;

  // Multi-foil identity protection.
  // The old algorithm tracked each ridge independently and could let a weak
  // foil jump onto a neighboring stronger ridge.  These guards preserve ridge
  // identity across delta rows.
  const double minInterFoilTrackGap = 0.20;
  const double maxCenterStepFrac = 0.60;

  const int minRowEntries = 40;

  TString inputroot = (inputRootOverrideStr.Length() > 0) ? inputRootOverrideStr : BuildInputRootPath(nrun, inputID);

  TString outCutFile = Form("%s/cuts/ytar_ridge_cut_%s.root", step01.Data(), outTag.Data());

  TString outPdf = Form("%s/plots/ytar_ridge_cut_%s.pdf", step01.Data(), outTag.Data());

  TString outCsv = Form("%s/tsv/ytar_ridge_cut_%s.tsv", step01.Data(), outTag.Data());

  gSystem->mkdir(step01 + "/cuts", kTRUE);
  gSystem->mkdir(step01 + "/plots", kTRUE);
  gSystem->mkdir(step01 + "/tsv", kTRUE);

  TFile *fin = TFile::Open(inputroot, "READ");
  if (!fin || fin->IsZombie()) {
    cerr << "ERROR: could not open input ROOT file: " << inputroot << endl;
    return;
  }

  TTree *T = (TTree*)fin->Get("T");
  if (!T) {
    cerr << "ERROR: could not find tree T" << endl;
    return;
  }

  Double_t sumnpe = 0.0;
  Double_t delta = 0.0;
  Double_t ytar = 0.0;

  if (!HasBranch(T, "H.cer.npeSum") ||
      !HasBranch(T, "H.gtr.dp") ||
      !HasBranch(T, "H.gtr.y")) {
    cerr << "ERROR: missing required branch." << endl;
    cerr << "Need H.cer.npeSum, H.gtr.dp, H.gtr.y" << endl;
    return;
  }

  T->SetBranchAddress("H.cer.npeSum", &sumnpe);
  T->SetBranchAddress("H.gtr.dp", &delta);
  T->SetBranchAddress("H.gtr.y", &ytar);

  vector<EventLite> events;

  TH2D *hBase = new TH2D(
    "hBase",
    Form("Run %d: base-cut events;ytar;delta", nrun),
    300, plotYtarMin, plotYtarMax,
    200, deltaMin, deltaMax
  );

  TH1D *hYtarCentral = new TH1D(
    "hYtarCentral",
    Form("Run %d: central-delta ytar density;ytar;counts", nrun),
    peakFindNYBins, peakFindYtarMin, peakFindYtarMax
  );

  Long64_t nentries = T->GetEntries();
  Long64_t nBase = 0;
  Long64_t nCentral = 0;

  for (Long64_t i = 0; i < nentries; i++) {
    T->GetEntry(i);

    if (!(sumnpe > 2.0)) continue;
    if (!(delta > deltaMin && delta < deltaMax)) continue;
    if (!std::isfinite(ytar)) continue;
    if (!std::isfinite(delta)) continue;

    events.push_back({ytar, delta});
    hBase->Fill(ytar, delta);
    nBase++;

    if (std::abs(delta) < centralDeltaAbsMax) {
      hYtarCentral->Fill(ytar);
      nCentral++;
    }
  }

  cout << "Base-cut events loaded: " << nBase << endl;
  cout << "Central events loaded: " << nCentral << endl;

  if (events.size() < 100 || nCentral < 50) {
    cerr << "ERROR: too few events." << endl;
    return;
  }

  TH1D *hYtarSmooth = (TH1D*)hYtarCentral->Clone("hYtarCentralSmooth");
  hYtarSmooth->SetTitle(Form("Run %d: smoothed central ytar density;ytar;smoothed counts", nrun));

  for (int s = 0; s < centralSmoothPasses; s++) {
    hYtarSmooth->Smooth(1);
  }

  vector<PeakInfo> peaks = FindMajorPeaks(
    hYtarSmooth,
    minPeakFracOfMax,
    minPeakSeparation
  );

  cout << "Detected major foil ridges = " << peaks.size() << endl;

  for (auto &p : peaks) {
    cout << "foil " << p.index
         << " peak_ytar=" << p.ytar
         << " height=" << p.height
         << " frac=" << p.fracOfMax
         << " searchHalfWidth=" << p.searchHalfWidth
         << endl;
  }

  if (peaks.empty()) {
    cerr << "ERROR: no foil ridges detected." << endl;
    return;
  }

  if (peaks.size() > 3) {
    cerr << "WARNING: detected more than 3 dominant ridges. Review diagnostic plot." << endl;
  }

  vector<TCutG*> expertCuts = LoadExpertCuts(nrun, handID);

  vector<RidgeCutResult> results;

  for (auto &p : peaks) {
    RidgeCutResult r = BuildCutForRidge(
      events,
      p,
      nrun,
      outTag,
      expertCuts,
      deltaMin,
      deltaMax,
      deltaStep,
      deltaHalfWindow,
      widthFracLevel,
      smoothPasses1D,
      smoothPassesEnvelope,
      centralDeltaAbsMax,
      minPeakFracOfGlobal,
      guardScale,
      weakRidgeFracThreshold,
      strongRidgeMinWidthScale,
      weakRidgeMinWidthScale,
      weakRidgeGuardBoost,
      rawWideningRetain,
      minRowEntries,
      &results,
      minInterFoilTrackGap,
      maxCenterStepFrac
    );

    if (r.cut) results.push_back(r);
  }

  if (results.empty()) {
    cerr << "ERROR: no auto cuts were built." << endl;
    return;
  }

  // ----------------------------
  // CSV output
  // ----------------------------
  ofstream csv(outCsv.Data());

  csv << "run,tag,inputFileID,foil_index,peak_ytar,peak_height,peak_frac_of_max,"
      << "searchHalfWidth,delta_center,nrow,"
      << "ridge_center_raw,peak_density,width_left_raw,width_right_raw,"
      << "ridge_center_smooth,width_left_final,width_right_final,"
      << "left_ytar,right_ytar,width_final,"
      << "refWidthLeft,refWidthRight,maxWidthLeft,maxWidthRight,centralWidthRows,"
      << "nAuto,nHand,nBoth\n";

  for (auto &r : results) {
    for (size_t i = 0; i < r.deltaCenters.size(); i++) {
      csv << nrun << ","
          << outTag << ","
          << inputID << ","
          << r.foilIndex << ","
          << r.peak.ytar << ","
          << r.peak.height << ","
          << r.peak.fracOfMax << ","
          << r.peak.searchHalfWidth << ","
          << r.deltaCenters[i] << ","
          << r.rowNVec[i] << ","
          << r.ridgeCenterRaw[i] << ","
          << r.peakValVec[i] << ","
          << r.widthLeftRaw[i] << ","
          << r.widthRightRaw[i] << ","
          << r.ridgeCenter[i] << ","
          << r.widthLeftFinal[i] << ","
          << r.widthRightFinal[i] << ","
          << r.leftBoundary[i] << ","
          << r.rightBoundary[i] << ","
          << r.rightBoundary[i] - r.leftBoundary[i] << ","
          << r.calib.refLeft << ","
          << r.calib.refRight << ","
          << r.calib.maxLeft << ","
          << r.calib.maxRight << ","
          << r.calib.nUsed << ","
          << r.nAuto << ","
          << r.nHand << ","
          << r.nBoth
          << "\n";
    }
  }

  csv.close();

  // ----------------------------
  // PDF diagnostics
  // ----------------------------
  TCanvas *c = new TCanvas("c", "multifoil central guard cut", 1100, 850);
  c->Print(outPdf + "[");

  // Page 1: central peak detection.
  c->Clear();

  hYtarCentral->SetLineColor(kGray+2);
  hYtarCentral->SetLineWidth(1);
  hYtarCentral->Draw("hist");

  hYtarSmooth->SetLineColor(kBlue);
  hYtarSmooth->SetLineWidth(3);
  hYtarSmooth->Draw("hist same");

  double ymax = hYtarCentral->GetMaximum();

  TLatex latex;
  latex.SetTextSize(0.035);

  vector<TLine*> peakLines;

  for (auto &p : peaks) {
    TLine *line = new TLine(p.ytar, 0.0, p.ytar, ymax);
    line->SetLineColor(kRed);
    line->SetLineWidth(3);
    line->SetLineStyle(2);
    line->Draw("same");
    peakLines.push_back(line);

    latex.DrawLatex(p.ytar + 0.05, 0.82 * ymax,
                    Form("foil %d", p.index));
  }

  TLegend *legPeak = new TLegend(0.56, 0.72, 0.88, 0.90);
  legPeak->SetTextSize(0.025);
  legPeak->SetFillColor(kWhite);
  legPeak->AddEntry(hYtarCentral, "central ytar density", "l");
  legPeak->AddEntry(hYtarSmooth, "smoothed density", "l");
  if (!peakLines.empty()) legPeak->AddEntry(peakLines[0], "detected ridges", "l");
  legPeak->Draw();

  c->Print(outPdf);

  // Page 2: all auto cuts + .
  c->Clear();
  gPad->SetLogz(1);

  hBase->SetTitle(Form("Run %d: auto multifoil cuts with ;ytar;delta", nrun));
  hBase->Draw("colz");

  for (auto &r : results) {
    if (r.cut) r.cut->Draw("L same");
  }

  for (auto *hc : expertCuts) {
    if (hc) hc->Draw("L same");
  }

  TLegend *legAll = new TLegend(0.55, 0.72, 0.88, 0.90);
  legAll->SetTextSize(0.025);
  legAll->SetFillColor(kWhite);

  if (!results.empty() && results[0].cut) legAll->AddEntry(results[0].cut, "auto cuts", "l");
  if (!expertCuts.empty()) legAll->AddEntry(expertCuts[0], "expert cuts", "l");
  legAll->Draw();

  c->Print(outPdf);
  gPad->SetLogz(0);

  // Prevent neighboring ridge cuts from overlapping after width padding/smoothing.
  // Boundary convention: x-axis is ytar. Results are sorted by foilIndex/ytar.
  const double minInterFoilGap = 0.04;
  sort(results.begin(), results.end(), [](const RidgeCutResult &a, const RidgeCutResult &b) { return a.peak.ytar < b.peak.ytar; });

  for (size_t ir = 0; ir + 1 < results.size(); ir++) {
    RidgeCutResult &leftR = results[ir];
    RidgeCutResult &rightR = results[ir+1];

    size_t n = std::min(leftR.deltaCenters.size(), rightR.deltaCenters.size());
    for (size_t i = 0; i < n; i++) {
      double mid = 0.5 * (leftR.ridgeCenter[i] + rightR.ridgeCenter[i]);
      double leftMax = mid - 0.5 * minInterFoilGap;
      double rightMin = mid + 0.5 * minInterFoilGap;

      if (leftR.rightBoundary[i] > leftMax) leftR.rightBoundary[i] = leftMax;
      if (rightR.leftBoundary[i] < rightMin) rightR.leftBoundary[i] = rightMin;
    }
  }

  // Rebuild TCutG objects after overlap protection.
  for (auto &r : results) {
    if (!r.cut) continue;
    delete r.cut;
    int npts = 2 * r.deltaCenters.size() + 1;
    r.cut = new TCutG(Form("delta_vs_ytar_cut_foil%d", r.foilIndex), npts);
    r.cut->SetTitle(Form("auto cut foil %d;ytar;delta", r.foilIndex));
    int ip = 0;
    for (size_t i = 0; i < r.deltaCenters.size(); i++) r.cut->SetPoint(ip++, r.leftBoundary[i], r.deltaCenters[i]);
    for (int i = (int)r.deltaCenters.size() - 1; i >= 0; i--) r.cut->SetPoint(ip++, r.rightBoundary[i], r.deltaCenters[i]);
    r.cut->SetPoint(ip++, r.leftBoundary[0], r.deltaCenters[0]);

    int color = kBlue + r.foilIndex;
    if (r.foilIndex == 0) color = kBlue;
    if (r.foilIndex == 1) color = kGreen+2;
    if (r.foilIndex == 2) color = kMagenta+2;
    r.cut->SetLineColor(color);
    r.cut->SetLineWidth(3);
  }


  // Recompute auto/hand/both counts after overlap protection, then rewrite TSV
  // so the table reflects the final boundaries.
  for (auto &r : results) {
    r.nAuto = 0;
    r.nHand = 0;
    r.nBoth = 0;

    TCutG *expert = nullptr;
    if (r.foilIndex >= 0 && r.foilIndex < (int)expertCuts.size()) {
      expert = expertCuts[r.foilIndex];
    }

    for (const auto &ev : events) {
      bool inAuto = r.cut ? r.cut->IsInside(ev.ytar, ev.delta) : false;
      bool inHand = expert ? expert->IsInside(ev.ytar, ev.delta) : false;

      if (inAuto) r.nAuto++;
      if (inHand) r.nHand++;
      if (inAuto && inHand) r.nBoth++;
    }
  }

  ofstream csvFinal(outCsv.Data());

  csvFinal << "run,tag,inputFileID,foil_index,peak_ytar,peak_height,peak_frac_of_max,"
           << "searchHalfWidth,delta_center,nrow,"
           << "ridge_center_raw,peak_density,width_left_raw,width_right_raw,"
           << "ridge_center_smooth,width_left_final,width_right_final,"
           << "left_ytar,right_ytar,width_final,"
           << "refWidthLeft,refWidthRight,maxWidthLeft,maxWidthRight,centralWidthRows,"
           << "nAuto,nHand,nBoth\n";

  for (auto &r : results) {
    for (size_t i = 0; i < r.deltaCenters.size(); i++) {
      csvFinal << nrun << ","
               << outTag << ","
               << inputID << ","
               << r.foilIndex << ","
               << r.peak.ytar << ","
               << r.peak.height << ","
               << r.peak.fracOfMax << ","
               << r.peak.searchHalfWidth << ","
               << r.deltaCenters[i] << ","
               << r.rowNVec[i] << ","
               << r.ridgeCenterRaw[i] << ","
               << r.peakValVec[i] << ","
               << r.widthLeftRaw[i] << ","
               << r.widthRightRaw[i] << ","
               << r.ridgeCenter[i] << ","
               << r.widthLeftFinal[i] << ","
               << r.widthRightFinal[i] << ","
               << r.leftBoundary[i] << ","
               << r.rightBoundary[i] << ","
               << r.rightBoundary[i] - r.leftBoundary[i] << ","
               << r.calib.refLeft << ","
               << r.calib.refRight << ","
               << r.calib.maxLeft << ","
               << r.calib.maxRight << ","
               << r.calib.nUsed << ","
               << r.nAuto << ","
               << r.nHand << ","
               << r.nBoth
               << "\n";
    }
  }

  csvFinal.close();

  cout << "Rewrote final TSV after overlap protection: " << outCsv << endl;

  // One detail page per foil.
  for (auto &r : results) {
    c->Clear();
    gPad->SetLogz(1);

    hBase->SetTitle(Form("Run %d: foil %d auto cut;ytar;delta", nrun, r.foilIndex));
    hBase->Draw("colz");

    if (r.cut) r.cut->Draw("L same");

    if (r.foilIndex >= 0 && r.foilIndex < (int)expertCuts.size()) {
      expertCuts[r.foilIndex]->Draw("L same");
    }

    TLatex tx;
    tx.SetNDC();
    tx.SetTextSize(0.026);
    tx.DrawLatex(0.14, 0.93,
                 Form("foil %d | peak ytar=%.4f | auto=%lld hand=%lld both=%lld",
                      r.foilIndex, r.peak.ytar, r.nAuto, r.nHand, r.nBoth));
    tx.DrawLatex(0.14, 0.89,
                 Form("refL=%.4f refR=%.4f maxL=%.4f maxR=%.4f nCentralRows=%d",
                      r.calib.refLeft, r.calib.refRight,
                      r.calib.maxLeft, r.calib.maxRight,
                      r.calib.nUsed));

    c->Print(outPdf);
    gPad->SetLogz(0);

    c->Clear();

    r.gLeft->SetTitle(Form("Run %d foil %d: boundaries and ridge center;delta;ytar",
                           nrun, r.foilIndex));

    r.gLeft->SetLineColor(kBlue);
    r.gLeft->SetLineWidth(2);

    r.gRight->SetLineColor(kBlue);
    r.gRight->SetLineWidth(2);
    r.gRight->SetLineStyle(2);

    r.gCenter->SetLineColor(kGreen+2);
    r.gCenter->SetLineWidth(2);

    r.gCenterRaw->SetLineColor(kGray+2);
    r.gCenterRaw->SetLineWidth(1);
    r.gCenterRaw->SetLineStyle(3);

    r.gLeft->Draw("AL");
    r.gLeft->GetXaxis()->SetLimits(-10.0, 10.0);
    r.gLeft->GetYaxis()->SetRangeUser(-3.0, 3.0);
    r.gRight->Draw("L same");
    r.gCenter->Draw("L same");
    r.gCenterRaw->Draw("L same");

    TLegend *leg = new TLegend(0.56, 0.74, 0.88, 0.90);
    leg->SetTextSize(0.025);
    leg->SetFillColor(kWhite);
    leg->AddEntry(r.gLeft, "left boundary", "l");
    leg->AddEntry(r.gRight, "right boundary", "l");
    leg->AddEntry(r.gCenter, "smoothed ridge center", "l");
    leg->AddEntry(r.gCenterRaw, "raw ridge center", "l");
    leg->Draw();

    c->Print(outPdf);
  }

  c->Print(outPdf + "]");

  // ----------------------------
  // ROOT output
  // ----------------------------
  TFile *fout = new TFile(outCutFile, "RECREATE");

  hBase->Write();
  hYtarCentral->Write();
  hYtarSmooth->Write();

  for (auto &r : results) {
    if (r.cut) r.cut->Write();

    if (r.gLeft) r.gLeft->Write(Form("left_boundary_vs_delta_foil%d", r.foilIndex));
    if (r.gRight) r.gRight->Write(Form("right_boundary_vs_delta_foil%d", r.foilIndex));
    if (r.gCenter) r.gCenter->Write(Form("ridge_center_vs_delta_foil%d", r.foilIndex));
    if (r.gCenterRaw) r.gCenterRaw->Write(Form("raw_ridge_center_vs_delta_foil%d", r.foilIndex));
    if (r.gWidthLeft) r.gWidthLeft->Write(Form("left_width_vs_delta_foil%d", r.foilIndex));
    if (r.gWidthRight) r.gWidthRight->Write(Form("right_width_vs_delta_foil%d", r.foilIndex));
  }

  for (size_t i = 0; i < expertCuts.size(); i++) {
    if (expertCuts[i]) expertCuts[i]->Write(Form("expert_cut_sorted_foil%zu", i));
  }

  fout->Close();

  cout << "Wrote cut ROOT file: " << outCutFile << endl;
  cout << "Wrote PDF: " << outPdf << endl;
  cout << "Wrote CSV: " << outCsv << endl;

  fin->Close();
}
