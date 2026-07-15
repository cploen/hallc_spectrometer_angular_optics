// STABLE-THETA PATCH: chooses theta using local-neighborhood stability-adjusted score.\n// DELETE-SLICE PATCH: removes all old yscol cuts for the same nfoil/ndel before rewriting.
// FIXED-THETA PATCH: fixedThetaDeg final argument; fixedThetaDeg >= -900 scans only that angle.
// BATCH PATCH: safe ndel naming + all-delta wrapper. Scoring unchanged.
// assign_yfp_ypfp_angleScanBands_strongWeak_with_pages_compat.C
// COMPAT PATCHED VERSION: writes hand-compatible AUTO BAND cuts.
// CLEAN-GAP SCORE PATCHED VERSION: angle score prefers clean valleys over raw peak count.
// Original algorithm preserved; interface/output layer patched.
//
// Diagnostic: choose the projection angle that maximizes resolved 1D peak structure.
// This directly tests: "rotate until the yfp/ypfp islands separate best."
//
// Coordinates:
//   x = ypfp
//   y = yfp
//   xz = (x - mean_x)/sigma_x
//   yz = (y - mean_y)/sigma_y
//
// For each theta in [0,180):
//   q = xz*cos(theta) + yz*sin(theta)        // candidate separation coordinate
//   p = -xz*sin(theta) + yz*cos(theta)       // perpendicular coordinate
//   histogram q, smooth, find peaks, score
//
//
// Useful test:
//   hcana -l -q 'assign_yfp_ypfp_angleScanBands.C(1544,-8,-5,"auto_ycut",9,1.0,0.18,0.05,0.15,2,0.005,0.04,0.45)'

#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TCutG.h>
#include <TKey.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TMarker.h>
#include <TGraph.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TMath.h>
#include <TObjString.h>
#include <TObjArray.h>
#include <TROOT.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <limits>

using std::cout;
using std::endl;

// ============================================================================
// COMPATIBILITY CONFIG: edit paths here, no other paths hiding in the middle of the macro.
// ============================================================================

static TString YFP_REPO_DIR = ".";
static TString YFP_DAT_DIR = "DATfiles";

static TString YFP_CAMPAIGN_DIR = "HMS_6p117";
static TString YFP_STEP_DIR     = YFP_CAMPAIGN_DIR + "/02a_angle_scan_y";
static TString YFP_CUTS_DIR     = YFP_STEP_DIR + "/cuts";
static TString YFP_PLOTS_DIR    = YFP_STEP_DIR + "/plots";
static TString YFP_ROOT_DIR     = YFP_STEP_DIR + "/root";
static TString YFP_TSV_DIR      = YFP_STEP_DIR + "/tsv";

static TString YFP_OUTPUT_TAG = "";
static TString YFP_INPUT_ROOT_OVERRIDE = "";

const TString YTAR_TAG = "ML_dev";

// Edit this when replay/rootfile location changes.
const TString YFP_TREE_ROOT_DIR =
  "/volatile/hallc/nps/cploen/ROOTfiles/OPTICS/angular_sandbox/newfit_6p667_20260526_1226_no_offsets/rootfiles/";

// qband-to-yscol mapping knobs.
// If labels come out reversed, flip this to kTRUE.
const Bool_t YFP_REVERSE_QBAND_TO_YSCOL = kFALSE;

// If labels are shifted, adjust this.
const Int_t YFP_YSCOL_OFFSET = 0;

// ============================================================================


struct OpticsRunInfo {
  int run = -1;
  TString opticsID = "";
  double centAngle = 0.0;
  int numFoil = 0;
  int sieveFlag = 0;
  int ndelcut = 0;
  std::vector<double> zfoil;
  std::vector<double> delcut;
};

struct QPeak {
  double q = 0.0;
  double height = 0.0;
  double prominence = 0.0;
  int bin = 0;
};

struct AngleResult {
  double thetaDeg = 0.0;
  double thetaRad = 0.0;
  int nPeaks = 0;
  double score = 0.0;
  double totalPeakHeight = 0.0;
  double meanValleyDrop = 0.0;

  // Clean-gap score diagnostics.
  // Prefer projections where adjacent peaks are separated by real low-density
  // valleys rather than fake extra local maxima.
  int nCleanGaps = 0;
  double meanGapQuality = 0.0;
  double minGapQuality = 0.0;
  double meanGapDrop = 0.0;
  double meanGapEmpty = 0.0;
  double meanGapBelowFrac = 0.0;
  double fakePeakPenalty = 0.0;

  std::vector<QPeak> peaks;
  std::vector<double> bounds;
};

static std::vector<TString> SplitCSV_angleScan(const TString& line) {
  std::vector<TString> out;
  TString s(line);
  TObjArray* arr = s.Tokenize(",");
  for (int i = 0; i < arr->GetEntries(); ++i) {
    TString tok = ((TObjString*)arr->At(i))->GetString();
    tok = tok.Strip(TString::kBoth);
    out.push_back(tok);
  }
  delete arr;
  return out;
}

static bool ReadOpticsRunInfo_angleScan(int nrun, OpticsRunInfo& info,
                                         const char* metaFile="DATfiles/list_of_optics_run.dat") {
  std::ifstream fin(metaFile);
  if (!fin.is_open()) {
    cout << "ERROR: cannot open metadata file: " << metaFile << endl;
    return false;
  }

  std::string raw;
  while (std::getline(fin, raw)) {
    TString line(raw.c_str());
    line = line.Strip(TString::kBoth);

    if (line.Length() == 0 || line.BeginsWith("#")) continue;

    auto tok = SplitCSV_angleScan(line);
    if (tok.size() < 6) continue;
    if (tok[0].Atoi() != nrun) continue;

    info.run       = tok[0].Atoi();
    info.opticsID  = tok[1];
    info.centAngle = tok[2].Atof();
    info.numFoil   = tok[3].Atoi();
    info.sieveFlag = tok[4].Atoi();
    info.ndelcut   = tok[5].Atoi();

    // Your file has ymis as token 6 on the header line.
    // We do not need it here.

    // Next non-empty line: z foil positions
    std::string zline_raw;
    if (!std::getline(fin, zline_raw)) {
      cout << "ERROR: missing zfoil line for run " << nrun << endl;
      return false;
    }

    TString zline(zline_raw.c_str());
    zline = zline.Strip(TString::kBoth);
    auto ztok = SplitCSV_angleScan(zline);

    for (int i = 0; i < info.numFoil && i < (int)ztok.size(); ++i) {
      info.zfoil.push_back(ztok[i].Atof());
    }

    // Next non-empty line: delta cut edges
    std::string dline_raw;
    if (!std::getline(fin, dline_raw)) {
      cout << "ERROR: missing delta-cut line for run " << nrun << endl;
      return false;
    }

    TString dline(dline_raw.c_str());
    dline = dline.Strip(TString::kBoth);
    auto dtok = SplitCSV_angleScan(dline);

    for (int i = 0; i < (int)dtok.size(); ++i) {
      info.delcut.push_back(dtok[i].Atof());
    }

    cout << "Parsed run " << nrun
         << " OpticsID=" << info.opticsID
         << " NumFoil=" << info.numFoil
         << " ndelcut=" << info.ndelcut
         << " zfoil entries=" << info.zfoil.size()
         << " delta edges=" << info.delcut.size()
         << endl;

    return true;
  }

  cout << "ERROR: run " << nrun << " not found in " << metaFile << endl;
  return false;
}

static TCutG* GetFirstCutG_angleScan(TFile* f) {
  if (!f || f->IsZombie()) return nullptr;
  TIter next(f->GetListOfKeys());
  TKey* key = nullptr;
  while ((key = (TKey*)next())) {
    TObject* obj = key->ReadObj();
    if (obj && obj->InheritsFrom(TCutG::Class())) return (TCutG*)obj;
  }
  return nullptr;
}

static TCutG* LoadYtarCut_angleScan(int nrun, const TString& ytarTag, int foilIndex=0) {
  TString ytarDir = YFP_CAMPAIGN_DIR + "/01_ytar_cuts/cuts";

  std::vector<TString> files = {
    Form("%s/ytar_ridge_cut_%s.root", ytarDir.Data(), ytarTag.Data()),
    Form("%s/ytar_ridge_cut_%s_run%d.root", ytarDir.Data(), ytarTag.Data(), nrun),
    Form("%s/ytar_ridge_cut_%s_run%d.root", YFP_CUTS_DIR.Data(), ytarTag.Data(), nrun)
  };

  std::vector<TString> names = {
    Form("delta_vs_ytar_cut_foil%d", foilIndex),
    Form("ytar_delta_cut_foil%d", foilIndex),
    Form("foil%d", foilIndex),
    Form("cut_foil%d", foilIndex),
    "delta_vs_ytar_cut",
    "ytar_delta_cut"
  };

  for (auto& fname : files) {
    TFile* f = TFile::Open(fname, "READ");
    if (!f || f->IsZombie()) continue;

    for (auto& name : names) {
      TCutG* c = (TCutG*)f->Get(name);
      if (c) {
        cout << "Loaded ytar cut: " << fname << " :: " << name << endl;
        return c;
      }
    }

    TCutG* first = GetFirstCutG_angleScan(f);
    if (first) {
      cout << "Loaded first available TCutG from " << fname << ": " << first->GetName() << endl;
      return first;
    }
  }

  cout << "WARNING: could not load ytar cut for tag=" << ytarTag
       << " foil=" << foilIndex << endl;
  cout << "Tried:" << endl;
  for (auto& fname : files) cout << "  " << fname << endl;

  return nullptr;
}

static double LocalValleyMin_angleScan(TH1D* h, int b1, int b2) {
  if (b2 < b1) std::swap(b1,b2);
  double valley = std::numeric_limits<double>::max();
  for (int b=b1; b<=b2; ++b)
    valley = std::min(valley, h->GetBinContent(b));
  return valley;
}

static std::vector<QPeak> FindPeaks1D_angleScan(TH1D* hSmooth,
                                                 int maxPeaks,
                                                 double minPeakSep,
                                                 double minPeakFraction,
                                                 double minProminenceFraction) {
  struct CandPeak {
    double q = 0.0;
    double height = 0.0;
    double prominence = 0.0;
    int bin = 0;
  };

  std::vector<CandPeak> localMaxima;
  const double maxContent = hSmooth->GetMaximum();
  const double minHeight = minPeakFraction * maxContent;

  // First collect all local maxima above height threshold.
  for (int ib=2; ib<hSmooth->GetNbinsX(); ++ib) {
    const double ym = hSmooth->GetBinContent(ib-1);
    const double y0 = hSmooth->GetBinContent(ib);
    const double yp = hSmooth->GetBinContent(ib+1);
    if (y0 > ym && y0 >= yp && y0 > minHeight) {
      localMaxima.push_back({hSmooth->GetBinCenter(ib), y0, 0.0, ib});
    }
  }

  // Estimate each maximum's valley prominence relative to neighboring maxima.
  // This rejects shoulder/tail ripples that are not separated by a real dip.
  std::sort(localMaxima.begin(), localMaxima.end(),
            [](const CandPeak& a, const CandPeak& b){ return a.q < b.q; });

  std::vector<QPeak> cands;
  for (size_t i=0; i<localMaxima.size(); ++i) {
    const double h0 = localMaxima[i].height;

    double leftProm = 1.0;
    double rightProm = 1.0;

    if (i > 0) {
      const double valley = LocalValleyMin_angleScan(hSmooth, localMaxima[i-1].bin, localMaxima[i].bin);
      leftProm = (h0 > 0.0) ? (h0 - valley)/h0 : 0.0;
    }

    if (i+1 < localMaxima.size()) {
      const double valley = LocalValleyMin_angleScan(hSmooth, localMaxima[i].bin, localMaxima[i+1].bin);
      rightProm = (h0 > 0.0) ? (h0 - valley)/h0 : 0.0;
    }

    // For interior peaks, require separation on both sides. For edge peaks, one side is enough.
    double prom = 0.0;
    if (i == 0 && localMaxima.size() > 1) prom = rightProm;
    else if (i+1 == localMaxima.size() && localMaxima.size() > 1) prom = leftProm;
    else if (localMaxima.size() == 1) prom = 1.0;
    else prom = std::min(leftProm, rightProm);

    if (prom >= minProminenceFraction) {
      cands.push_back({localMaxima[i].q, h0, prom, localMaxima[i].bin});
    }
  }

  // Keep strongest/prominent peaks first, then suppress duplicates too close in q.
  std::sort(cands.begin(), cands.end(),
            [](const QPeak& a, const QPeak& b){
              if (a.prominence != b.prominence) return a.prominence > b.prominence;
              return a.height > b.height;
            });

  std::vector<QPeak> peaks;
  for (const auto& cand : cands) {
    bool tooClose=false;
    for (const auto& pk : peaks) {
      if (std::fabs(cand.q - pk.q) < minPeakSep) { tooClose=true; break; }
    }
    if (tooClose) continue;
    peaks.push_back(cand);
    if ((int)peaks.size() >= maxPeaks) break;
  }

  std::sort(peaks.begin(), peaks.end(),
            [](const QPeak& a, const QPeak& b){ return a.q < b.q; });

  return peaks;
}

static std::vector<double> BoundariesFromPeaks_angleScan(const std::vector<QPeak>& peaks) {
  std::vector<double> bounds;
  for (size_t i=0; i+1<peaks.size(); ++i)
    bounds.push_back(0.5*(peaks[i].q + peaks[i+1].q));
  return bounds;
}

static double MeanValleyDrop_angleScan(TH1D* hSmooth, const std::vector<QPeak>& peaks) {
  if (peaks.size() < 2) return 0.0;

  double sumDrop = 0.0;
  int n = 0;

  for (size_t i=0; i+1<peaks.size(); ++i) {
    int b1 = peaks[i].bin;
    int b2 = peaks[i+1].bin;
    if (b2 < b1) std::swap(b1,b2);

    double valley = std::numeric_limits<double>::max();
    for (int b=b1; b<=b2; ++b)
      valley = std::min(valley, hSmooth->GetBinContent(b));

    const double lowPeak = std::min(peaks[i].height, peaks[i+1].height);
    if (lowPeak > 0 && std::isfinite(valley)) {
      sumDrop += (lowPeak - valley) / lowPeak; // 0=no valley, 1=deep valley
      n++;
    }
  }

  return (n > 0) ? sumDrop / (double)n : 0.0;
}


// ============================================================================
// Clean-gap scoring helpers
//
// Philosophy:
//   The correct projection is not necessarily the one with the most local
//   maxima. It is the one with the largest number of physically plausible
//   peak-to-peak gaps that are genuinely low-density.
//
// These cuts are intentionally permissive because this stage creates candidate
// row/column bands. GMM cleanup later removes bridge/tail/noisy events.
// ============================================================================

struct CleanGapSummary {
  int nPairs = 0;
  int nCleanGaps = 0;

  double meanQuality = 0.0;
  double minQuality = 0.0;

  double meanDrop = 0.0;
  double meanEmpty = 0.0;
  double meanBelowFrac = 0.0;

  double penalty = 0.0;
};

static double Clamp01_angleScan(double x)
{
  if (x < 0.0) return 0.0;
  if (x > 1.0) return 1.0;
  return x;
}

static CleanGapSummary EvaluateCleanGaps_angleScan(TH1D* hSmooth,
                                                    const std::vector<QPeak>& peaks,
                                                    double minCleanDrop,
                                                    double minCleanEmpty,
                                                    double minPeakSep)
{
  CleanGapSummary out;

  if (!hSmooth || peaks.size() < 2) return out;

  bool firstQuality = true;

  double sumQuality = 0.0;
  double sumDrop = 0.0;
  double sumEmpty = 0.0;
  double sumBelowFrac = 0.0;

  for (size_t i = 0; i + 1 < peaks.size(); ++i) {
    int b1 = peaks[i].bin;
    int b2 = peaks[i+1].bin;
    if (b2 < b1) std::swap(b1, b2);

    if (b2 <= b1) continue;

    const double qSep = std::fabs(peaks[i+1].q - peaks[i].q);
    const double lowPeak = std::min(peaks[i].height, peaks[i+1].height);

    if (lowPeak <= 0.0) continue;

    out.nPairs++;

    double valleyMin = std::numeric_limits<double>::max();
    int valleyBin = b1;

    for (int b = b1; b <= b2; ++b) {
      double v = hSmooth->GetBinContent(b);
      if (v < valleyMin) {
        valleyMin = v;
        valleyBin = b;
      }
    }

    const int nBetween = std::max(1, b2 - b1 + 1);

    // Local valley occupancy near the deepest point.
    // This prevents one empty bin from winning if the whole gap is actually full.
    const int halfWin = std::max(1, std::min(4, nBetween / 8));
    int vlo = std::max(b1, valleyBin - halfWin);
    int vhi = std::min(b2, valleyBin + halfWin);

    double valleyLocalSum = 0.0;
    int valleyLocalN = 0;

    for (int b = vlo; b <= vhi; ++b) {
      valleyLocalSum += hSmooth->GetBinContent(b);
      valleyLocalN++;
    }

    double valleyLocalMean = (valleyLocalN > 0) ? valleyLocalSum / (double)valleyLocalN : valleyMin;

    // Fractional depth of deepest valley relative to the smaller adjacent peak.
    double drop = Clamp01_angleScan((lowPeak - valleyMin) / lowPeak);

    // Local emptiness around the valley bottom.
    double empty = Clamp01_angleScan((lowPeak - valleyLocalMean) / lowPeak);

    // Fraction of the full inter-peak region below a modest low-density threshold.
    // This rewards a visibly open lane but does not require it to be completely empty.
    const double lowDensityThreshold = 0.35 * lowPeak;
    int nBelow = 0;

    for (int b = b1; b <= b2; ++b) {
      if (hSmooth->GetBinContent(b) < lowDensityThreshold) nBelow++;
    }

    double belowFrac = (double)nBelow / (double)nBetween;

    // Separation factor suppresses fake peak splitting where peaks are too close.
    double sepFactor = 1.0;
    if (minPeakSep > 0.0) {
      sepFactor = Clamp01_angleScan((qSep - minPeakSep) / (0.75 * minPeakSep));
    }

    // Gap quality: depth + local emptiness + width of low-density lane.
    double quality =
        (0.45 * drop
       + 0.35 * empty
       + 0.20 * belowFrac) * sepFactor;

    bool clean =
      (qSep >= minPeakSep) &&
      (drop >= minCleanDrop) &&
      (empty >= minCleanEmpty);

    if (clean) out.nCleanGaps++;

    // Penalties are intentionally moderate. We want to demote fake splitting,
    // not throw away imperfect xfp/xpfp cases that GMM can clean.
    if (qSep < minPeakSep && minPeakSep > 0.0) {
      out.penalty += 50000.0 * (1.0 - qSep / minPeakSep);
    }

    if (drop < 0.15) {
      out.penalty += 20000.0 * (0.15 - drop) / 0.15;
    }

    sumQuality += quality;
    sumDrop += drop;
    sumEmpty += empty;
    sumBelowFrac += belowFrac;

    if (firstQuality) {
      out.minQuality = quality;
      firstQuality = false;
    } else {
      out.minQuality = std::min(out.minQuality, quality);
    }
  }

  if (out.nPairs > 0) {
    out.meanQuality = sumQuality / (double)out.nPairs;
    out.meanDrop = sumDrop / (double)out.nPairs;
    out.meanEmpty = sumEmpty / (double)out.nPairs;
    out.meanBelowFrac = sumBelowFrac / (double)out.nPairs;
  }

  return out;
}

// ============================================================================

static AngleResult ScoreAngle_angleScan(const std::vector<double>& xzvals,
                                         const std::vector<double>& yzvals,
                                         double thetaDeg,
                                         int maxBands,
                                         double minPeakSep,
                                         double minPeakFraction,
                                         double minProminenceFraction,
                                         int smoothPasses,
                                         int nBinsQ = 240) {
  const Long64_t N = (Long64_t)xzvals.size();
  const double th = thetaDeg * TMath::Pi()/180.0;
  const double c = std::cos(th);
  const double s = std::sin(th);

  std::vector<double> qvals;
  qvals.reserve(N);
  for (Long64_t i=0; i<N; ++i) qvals.push_back(xzvals[i]*c + yzvals[i]*s);

  double qmin=*std::min_element(qvals.begin(), qvals.end());
  double qmax=*std::max_element(qvals.begin(), qvals.end());
  double qpad=0.05*(qmax-qmin);
  qmin-=qpad; qmax+=qpad;

  TH1D* hQ = new TH1D(Form("hQ_scan_tmp_%g", thetaDeg), "temporary q scan", nBinsQ, qmin, qmax);
  hQ->SetDirectory(nullptr);
  for (double q : qvals) hQ->Fill(q);

  TH1D* hS = (TH1D*)hQ->Clone(Form("hQ_scan_tmp_smooth_%g", thetaDeg));
  hS->SetDirectory(nullptr);
  for (int i=0; i<smoothPasses; ++i) hS->Smooth(1);

  AngleResult res;
  res.thetaDeg = thetaDeg;
  res.thetaRad = th;
  res.peaks = FindPeaks1D_angleScan(hS, maxBands, minPeakSep, minPeakFraction, minProminenceFraction);
  res.bounds = BoundariesFromPeaks_angleScan(res.peaks);
  res.nPeaks = (int)res.peaks.size();
  for (const auto& pk : res.peaks) res.totalPeakHeight += pk.height;
  res.meanValleyDrop = MeanValleyDrop_angleScan(hS, res.peaks);

  // Clean-gap-dominated score.
  //
  // We no longer reward raw peak count as the primary objective because some
  // bad rotations create extra local maxima by piling structures oddly.
  //
  // Instead, reward adjacent peak pairs that have real low-density valleys
  // between them. This still tolerates imperfect valleys because GMM cleanup
  // is expected downstream, especially for xfp/xpfp.
  double meanProm = 0.0;
  for (const auto& pk : res.peaks) meanProm += pk.prominence;
  if (!res.peaks.empty()) meanProm /= (double)res.peaks.size();

  const double minCleanGapDrop  = 0.30;  // permissive: GMM cleans later
  const double minCleanGapEmpty = 0.25;  // permissive: allows bridge/tail events

  CleanGapSummary gapSummary =
    EvaluateCleanGaps_angleScan(hS,
                                 res.peaks,
                                 minCleanGapDrop,
                                 minCleanGapEmpty,
                                 minPeakSep);

  res.nCleanGaps = gapSummary.nCleanGaps;
  res.meanGapQuality = gapSummary.meanQuality;
  res.minGapQuality = gapSummary.minQuality;
  res.meanGapDrop = gapSummary.meanDrop;
  res.meanGapEmpty = gapSummary.meanEmpty;
  res.meanGapBelowFrac = gapSummary.meanBelowFrac;
  res.fakePeakPenalty = gapSummary.penalty;

  // VALLEY-FIRST SCORE PATCH, ported from the evolved XFP/XPFP logic.
  //
  // Prefer angles with clean low-density valleys between bands, not merely
  // large raw peak count. This makes the selected theta less sensitive to
  // fake local maxima and shoulder splitting.
  const double valleyMetric = std::max(res.meanValleyDrop, res.meanGapDrop);

  double lowValleyPenalty = 0.0;
  if (valleyMetric < 0.55) {
    lowValleyPenalty += 3000000.0 * (0.55 - valleyMetric) / 0.55;
  }
  if (res.meanGapQuality < 0.30) {
    lowValleyPenalty += 1500000.0 * (0.30 - res.meanGapQuality) / 0.30;
  }

  res.score = 2500000.0 * valleyMetric
            + 2000000.0 * res.meanGapQuality
            + 1200000.0 * res.meanGapEmpty
            + 800000.0  * res.minGapQuality
            + 250000.0  * res.nCleanGaps
            + 50000.0   * meanProm
            + 10000.0   * res.nPeaks
            + 0.05      * res.totalPeakHeight
            - res.fakePeakPenalty
            - lowValleyPenalty;

  delete hQ;
  delete hS;
  return res;
}


// ============================================================================
// COMPATIBILITY HELPERS
// ============================================================================

static TString BuildYfpTreeRootPath_angleScan(int nrun, int FileID)
{
  if (YFP_INPUT_ROOT_OVERRIDE.Length() > 0) return YFP_INPUT_ROOT_OVERRIDE;

  return Form("%s/nps_hms_optics_%d_1_%d.root",
              YFP_TREE_ROOT_DIR.Data(), nrun, FileID);
}

static TString YfpOutputTag_angleScan(const OpticsRunInfo& info)
{
  return (YFP_OUTPUT_TAG.Length() > 0) ? YFP_OUTPUT_TAG : info.opticsID;
}

static TString BuildYfpAutoBandCutFile_angleScan(const OpticsRunInfo& info)
{
  TString tag = YfpOutputTag_angleScan(info);
  return Form("%s/YpFpYFp_%s_auto_band_cut.root",
              YFP_CUTS_DIR.Data(), tag.Data());
}

static TString BuildYfpLegacyHandCutFile_angleScan(const OpticsRunInfo& info)
{
  TString tag = YfpOutputTag_angleScan(info);
  return Form("%s/YpFpYFp_%s_cut.root",
              YFP_CUTS_DIR.Data(), tag.Data());
}

static TString BuildYfpMergedBandCutFile_angleScan(const OpticsRunInfo& info)
{
  TString tag = YfpOutputTag_angleScan(info);
  return Form("%s/YpFpYFp_%s_merged_band_cut.root",
              YFP_CUTS_DIR.Data(), tag.Data());
}

static TString BuildYfpDiagnosticBase_angleScan(const OpticsRunInfo& info,
                                                 double deltaMin,
                                                 double deltaMax,
                                                 int foilIndex,
                                                 const TString& outDir)
{
  TString tag = YfpOutputTag_angleScan(info);
  TString outbase = Form("%s/yfp_ypfp_angleScan_%s_foil%d_delta_%g_to_%g",
                         outDir.Data(), tag.Data(), foilIndex, deltaMin, deltaMax);
  outbase.ReplaceAll("-", "m");
  outbase.ReplaceAll(".", "p");
  return outbase;
}

static Int_t FindDeltaIndex_angleScan(const OpticsRunInfo& info,
                                       double deltaMin,
                                       double deltaMax,
                                       Int_t ndelOverride=-999)
{
  if (ndelOverride >= 0) return ndelOverride;

  const double tol = 1.0e-4;

  for (int nd = 0; nd + 1 < (int)info.delcut.size(); ++nd) {
    if (std::fabs(info.delcut[nd]   - deltaMin) < tol &&
        std::fabs(info.delcut[nd+1] - deltaMax) < tol) {
      return nd;
    }
  }

  double targetCenter = 0.5 * (deltaMin + deltaMax);
  double bestDist = 1.0e99;
  int bestNd = -1;

  for (int nd = 0; nd + 1 < (int)info.delcut.size(); ++nd) {
    double center = 0.5 * (info.delcut[nd] + info.delcut[nd+1]);
    double dist = std::fabs(center - targetCenter);
    if (dist < bestDist) {
      bestDist = dist;
      bestNd = nd;
    }
  }

  if (bestNd >= 0) {
    cout << "WARNING: exact delta slice not found in metadata." << endl;
    cout << "         Requested [" << deltaMin << ", " << deltaMax << ")" << endl;
    cout << "         Using closest ndel = " << bestNd
         << " from metadata slice ["
         << info.delcut[bestNd] << ", " << info.delcut[bestNd+1] << ")" << endl;
  } else {
    cout << "ERROR: could not determine ndel index." << endl;
  }

  return bestNd;
}

static Int_t YscolFromQBand_angleScan(Int_t qband, Int_t nBands)
{
  Int_t yscol = YFP_REVERSE_QBAND_TO_YSCOL ? (nBands - 1 - qband) : qband;
  yscol += YFP_YSCOL_OFFSET;
  return yscol;
}

static TCutG* MakeQBandCut_angleScan(const TString& cutName,
                                      double qlo,
                                      double qhi,
                                      double pmin,
                                      double pmax,
                                      double mx,
                                      double my,
                                      double sx,
                                      double sy,
                                      double c,
                                      double s)
{
  auto rawX = [&](double q, double p) {
    double xz = q*c - p*s;
    return mx + sx*xz;
  };

  auto rawY = [&](double q, double p) {
    double yz = q*s + p*c;
    return my + sy*yz;
  };

  TCutG* cut = new TCutG(cutName, 5);
  cut->SetTitle(cutName + ";ypfp;yfp");

  cut->SetPoint(0, rawX(qlo, pmin), rawY(qlo, pmin));
  cut->SetPoint(1, rawX(qlo, pmax), rawY(qlo, pmax));
  cut->SetPoint(2, rawX(qhi, pmax), rawY(qhi, pmax));
  cut->SetPoint(3, rawX(qhi, pmin), rawY(qhi, pmin));
  cut->SetPoint(4, rawX(qlo, pmin), rawY(qlo, pmin));

  cut->SetLineColor(kBlue+1);
  cut->SetLineWidth(3);

  return cut;
}

// ============================================================================

// ============================================================================
// Batch patch helper: safe ndel finder.
//
// Prefer DATfile edge matching over stale manual ndelOverride.
// This prevents accidentally writing every delta slice as _ndel_1.
// ============================================================================
static int FindDeltaIndex_angleScan_v2(const OpticsRunInfo& info,
                                       double deltaMin,
                                       double deltaMax,
                                       int ndelOverride)
{
  const double tol = 1.0e-3;
  int autoIndex = -1;

  if (info.delcut.size() >= 2) {
    for (int i = 0; i + 1 < (int)info.delcut.size(); ++i) {
      const double lo = info.delcut[i];
      const double hi = info.delcut[i+1];

      if (std::fabs(deltaMin - lo) < tol &&
          std::fabs(deltaMax - hi) < tol) {
        autoIndex = i;
        break;
      }
    }
  }

  if (autoIndex >= 0) {
    if (ndelOverride >= 0 && ndelOverride != autoIndex) {
      cout << "WARNING: manual ndelOverride=" << ndelOverride
           << " conflicts with auto DATfile ndel=" << autoIndex
           << ". Using auto ndel to avoid cut-name overwrite."
           << endl;
    }
    return autoIndex;
  }

  if (ndelOverride >= 0) {
    cout << "WARNING: could not auto-match delta slice ["
         << deltaMin << ", " << deltaMax
         << ") to DATfile edges. Using manual ndelOverride="
         << ndelOverride << endl;
    return ndelOverride;
  }

  cout << "ERROR: could not determine ndel index for delta slice ["
       << deltaMin << ", " << deltaMax << ")" << endl;

  cout << "DATfile delta edges:";
  for (size_t i = 0; i < info.delcut.size(); ++i) {
    cout << " " << info.delcut[i];
  }
  cout << endl;

  return -1;
}


void run_yfp_ypfp_angleScan_oneDelta(Int_t nrun=1544,
                                                Double_t deltaMin=-10.0,
                                                Double_t deltaMax=-8.0,
                                                TString ytarTag="auto_ycut",
                                                Int_t maxBands=9,
                                                Double_t thetaStepDeg=1.0,
                                                Double_t minPeakSep=0.18,
                                                Double_t minPeakFraction=0.05,
                                                Double_t minProminenceFraction=0.15,
                                                Int_t smoothPasses=2,
                                                Double_t minWeakPeakFraction=0.005,
                                                Double_t minWeakProminenceFraction=0.08,
                                                Double_t weakMinSepFactor=1.0,
                                                Bool_t useYtarCut=true,
                                                Int_t foilIndex=0,
                                                Long64_t maxEvents=-1,
                                                Int_t FileID=-1,
                                                Int_t ndelOverride=-999,
                     Double_t fixedThetaDeg=-999.0,
                     const char *campaignDir="HMS_6p117",
                     const char *inputRootOverride="") {
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kBird);

  YFP_CAMPAIGN_DIR = campaignDir;
  if (YFP_CAMPAIGN_DIR.EndsWith("/")) YFP_CAMPAIGN_DIR.Chop();

  YFP_STEP_DIR  = YFP_CAMPAIGN_DIR + "/02a_angle_scan_y";
  YFP_CUTS_DIR  = YFP_STEP_DIR + "/cuts";
  YFP_PLOTS_DIR = YFP_STEP_DIR + "/plots";
  YFP_ROOT_DIR  = YFP_STEP_DIR + "/root";
  YFP_TSV_DIR   = YFP_STEP_DIR + "/tsv";

  YFP_OUTPUT_TAG = ytarTag;
  YFP_INPUT_ROOT_OVERRIDE = inputRootOverride;

  OpticsRunInfo info;
  TString metaFile = YFP_DAT_DIR + "/list_of_optics_run.dat";
  if (!ReadOpticsRunInfo_angleScan(nrun, info, metaFile.Data())) return;

  Int_t ndelIndex = FindDeltaIndex_angleScan_v2(info, deltaMin, deltaMax, ndelOverride);
  if (ndelIndex < 0) return;

  TString inroot = BuildYfpTreeRootPath_angleScan(nrun, FileID);

  TString outCutRootAutoBand = BuildYfpAutoBandCutFile_angleScan(info);
  TString outCutRootHand     = BuildYfpLegacyHandCutFile_angleScan(info);
  TString outCutRootMerged   = BuildYfpMergedBandCutFile_angleScan(info);

  TString outbasePdf  = BuildYfpDiagnosticBase_angleScan(info, deltaMin, deltaMax, foilIndex, YFP_PLOTS_DIR);
  TString outbaseRoot = BuildYfpDiagnosticBase_angleScan(info, deltaMin, deltaMax, foilIndex, YFP_ROOT_DIR);
  TString outpdf  = outbasePdf + ".pdf";
  TString outroot = outbaseRoot + ".root";

  gSystem->mkdir(YFP_STEP_DIR.Data(), kTRUE);
  gSystem->mkdir(YFP_CUTS_DIR.Data(), kTRUE);
  gSystem->mkdir(YFP_PLOTS_DIR.Data(), kTRUE);
  gSystem->mkdir(YFP_ROOT_DIR.Data(), kTRUE);
  gSystem->mkdir(YFP_TSV_DIR.Data(), kTRUE);

  cout << "\n=== Compatible YFP/YPFP band-cut interface ===" << endl;
  cout << "Input tree ROOT file       : " << inroot << endl;
  cout << "AUTO BAND cut output       : " << outCutRootAutoBand << endl;
  cout << "Legacy HAND cut expected   : " << outCutRootHand << endl;
  cout << "Suggested MERGED BAND file : " << outCutRootMerged << endl;
  cout << "Diagnostic PDF             : " << outpdf << endl;
  cout << "Diagnostic ROOT            : " << outroot << endl;
  cout << "OpticsID                   : " << info.opticsID << endl;
  cout << "FileID                     : " << FileID << endl;
  cout << "foilIndex / nfoil          : " << foilIndex << endl;
  cout << "ndel index                 : " << ndelIndex << endl;
  cout << "=============================================\n" << endl;
  TFile* fin = TFile::Open(inroot, "READ");
  if (!fin || fin->IsZombie()) {
    cout << "ERROR: cannot open input ROOT file: " << inroot << endl;
    return;
  }

  TTree* T = (TTree*)fin->Get("T");
  if (!T) {
    cout << "ERROR: tree T not found in " << inroot << endl;
    return;
  }

  TCutG* ytarCut = nullptr;
  if (useYtarCut) ytarCut = LoadYtarCut_angleScan(nrun, ytarTag, foilIndex);

  Double_t sumnpe=0, etracknorm=0;
  Double_t ytar=0, delta=0, yfp=0, ypfp=0;

  T->SetBranchStatus("*",0);
  T->SetBranchStatus("H.cer.npeSum",1);
  T->SetBranchStatus("H.cal.etottracknorm",1);
  T->SetBranchStatus("H.gtr.y",1);
  T->SetBranchStatus("H.gtr.dp",1);
  T->SetBranchStatus("H.dc.y_fp",1);
  T->SetBranchStatus("H.dc.yp_fp",1);

  T->SetBranchAddress("H.cer.npeSum", &sumnpe);
  T->SetBranchAddress("H.cal.etottracknorm", &etracknorm);
  T->SetBranchAddress("H.gtr.y", &ytar);
  T->SetBranchAddress("H.gtr.dp", &delta);
  T->SetBranchAddress("H.dc.y_fp", &yfp);
  T->SetBranchAddress("H.dc.yp_fp", &ypfp);

  // Plot convention: x = ypfp, y = yfp.
  std::vector<double> xvals; // ypfp
  std::vector<double> yvals; // yfp

  Long64_t nentries = T->GetEntries();
  if (maxEvents > 0 && maxEvents < nentries) nentries = maxEvents;

  Long64_t nPassBasic=0, nPassYtar=0, nPassDelta=0;
  for (Long64_t i = 0; i < nentries; ++i) {
    T->GetEntry(i);

    if (!(sumnpe > 6.0 && etracknorm > 0.65)) continue;
    nPassBasic++;

    if (useYtarCut && ytarCut) {
      if (!ytarCut->IsInside(ytar, delta)) continue; // x=ytar, y=delta
    }
    nPassYtar++;

    if (!(delta >= deltaMin && delta < deltaMax)) continue;
    nPassDelta++;

    xvals.push_back(ypfp);
    yvals.push_back(yfp);
  }

  const Long64_t N = xvals.size();

  cout << "\n=== Angle-scan yfp/ypfp band diagnostic ===" << endl;
  cout << "Run: " << nrun << endl;
  cout << "Input: " << inroot << endl;
  cout << "Metadata opticsID: " << info.opticsID << ", NumFoil: " << info.numFoil << endl;
  cout << "Delta slice: [" << deltaMin << ", " << deltaMax << ") %" << endl;
  cout << "Max allowed bands: " << maxBands << endl;
  cout << "Theta step: " << thetaStepDeg << " deg" << endl;
  cout << "Min projected peak separation: " << minPeakSep << endl;
  cout << "Min projected peak height fraction: " << minPeakFraction << endl;
  cout << "Min projected peak prominence fraction: " << minProminenceFraction << endl;
  cout << "Smoothing passes: " << smoothPasses << endl;
  cout << "Weak peak fraction after angle choice: " << minWeakPeakFraction << endl;
  cout << "Weak peak prominence after angle choice: " << minWeakProminenceFraction << endl;
  cout << "Weak peak min separation factor after angle choice: " << weakMinSepFactor << endl;
  cout << "Events after PID/basic: " << nPassBasic << endl;
  cout << "Events after ytar cut:  " << nPassYtar << endl;
  cout << "Events in delta slice: " << nPassDelta << endl;
  cout << "Events used:           " << N << endl;

  if (N < 100) {
    cout << "ERROR: too few events." << endl;
    return;
  }

  double mx=0, my=0;
  for (Long64_t i=0; i<N; ++i) { mx += xvals[i]; my += yvals[i]; }
  mx /= (double)N;
  my /= (double)N;

  double raw_sxx=0, raw_syy=0, raw_sxy=0;
  for (Long64_t i=0; i<N; ++i) {
    const double dx = xvals[i] - mx;
    const double dy = yvals[i] - my;
    raw_sxx += dx*dx;
    raw_syy += dy*dy;
    raw_sxy += dx*dy;
  }
  raw_sxx /= (double)(N-1);
  raw_syy /= (double)(N-1);
  raw_sxy /= (double)(N-1);

  const double sx = std::sqrt(raw_sxx);
  const double sy = std::sqrt(raw_syy);
  if (sx <= 0 || sy <= 0) {
    cout << "ERROR: zero RMS in one coordinate." << endl;
    return;
  }

  std::vector<double> xzvals, yzvals;
  xzvals.reserve(N); yzvals.reserve(N);
  for (Long64_t i=0; i<N; ++i) {
    xzvals.push_back((xvals[i] - mx)/sx);
    yzvals.push_back((yvals[i] - my)/sy);
  }

  // Scan all projection angles.
  std::vector<AngleResult> scan;
  const double thetaStartScan = (fixedThetaDeg > -900.0 ? fixedThetaDeg : 0.0);
  const double thetaEndScan   = (fixedThetaDeg > -900.0 ? fixedThetaDeg : 180.0);
  const double thetaStepUse   = (fixedThetaDeg > -900.0 ? 1.0 : thetaStepDeg);

  if (fixedThetaDeg > -900.0) {
    cout << "FIXED-THETA MODE: using theta = " << fixedThetaDeg
         << " deg; skipping theta scan." << endl;
  }

  for (double th = thetaStartScan;
       th <= thetaEndScan + 1.0e-9;
       th += thetaStepUse) {
    scan.push_back(ScoreAngle_angleScan(xzvals, yzvals, th,
                                         maxBands, minPeakSep,
                                         minPeakFraction, minProminenceFraction,
                                         smoothPasses));
  }

  if (scan.empty()) {
    cout << "ERROR: no scan results." << endl;
    return;
  }

  // -------------------------------------------------------------------------
  // Stability-aware theta selection.
  //
  // Raw score alone can choose an isolated angle where one cluster splits or a
  // valley-drop diagnostic has a sharp local artifact. Prefer angles sitting
  // on a local plateau: nearby angles should have comparable score, similar
  // clean-gap count, similar peak count, and no sharp valley-drop instability.
  //
  // Fixed-theta mode still works: scan has one entry, so stableScore == rawScore.
  // -------------------------------------------------------------------------
  std::vector<double> stableScores(scan.size(), 0.0);

  const int neighborWindow = 2;          // +/- two scan points, normally +/-2 degrees
  const double scoreRelTol = 0.12;       // neighbor raw score within 12% counts as plateau-like
  const double valleyJumpTol = 0.18;     // large meanValleyDrop jumps are suspicious
  const double plateauReward = 75000.0;  // modest: enough to break bad local spikes
  const double sameCleanGapReward = 90000.0;
  const double samePeakReward = 45000.0;
  const double peakJumpPenalty = 120000.0;
  const double valleyJumpPenalty = 100000.0;
  const double isolatedSpikePenalty = 250000.0;

  for (size_t i = 0; i < scan.size(); ++i) {
    double stableScore = scan[i].score;

    int nNeighbors = 0;
    int nScoreNeighbors = 0;
    int nSameCleanGap = 0;
    int nSamePeak = 0;

    double neighborScoreSum = 0.0;
    double valleyJumpSum = 0.0;
    int peakJumpSum = 0;

    int jlo = std::max<int>(0, (int)i - neighborWindow);
    int jhi = std::min<int>((int)scan.size() - 1, (int)i + neighborWindow);

    for (int j = jlo; j <= jhi; ++j) {
      if (j == (int)i) continue;

      nNeighbors++;
      neighborScoreSum += scan[j].score;

      const double denom = std::max(1.0, std::fabs(scan[i].score));
      const double relDiff = std::fabs(scan[j].score - scan[i].score) / denom;

      if (relDiff <= scoreRelTol) nScoreNeighbors++;
      if (scan[j].nCleanGaps == scan[i].nCleanGaps) nSameCleanGap++;
      if (scan[j].nPeaks == scan[i].nPeaks) nSamePeak++;

      peakJumpSum += std::abs(scan[j].nPeaks - scan[i].nPeaks);
      valleyJumpSum += std::fabs(scan[j].meanValleyDrop - scan[i].meanValleyDrop);
    }

    double meanNeighborScore = (nNeighbors > 0) ? neighborScoreSum / (double)nNeighbors : scan[i].score;
    double meanValleyJump = (nNeighbors > 0) ? valleyJumpSum / (double)nNeighbors : 0.0;

    stableScore += plateauReward * nScoreNeighbors;
    stableScore += sameCleanGapReward * nSameCleanGap;
    stableScore += samePeakReward * nSamePeak;

    stableScore -= peakJumpPenalty * peakJumpSum;

    if (meanValleyJump > valleyJumpTol) {
      stableScore -= valleyJumpPenalty * (meanValleyJump - valleyJumpTol) / valleyJumpTol;
    }

    if (nNeighbors > 0 &&
        scan[i].score > meanNeighborScore * (1.0 + scoreRelTol) &&
        nSameCleanGap < std::max(1, nNeighbors / 2)) {
      stableScore -= isolatedSpikePenalty;
    }

    stableScores[i] = stableScore;
  }

  // VALLEY-FIRST THETA SELECTION PATCH, ported from XFP/XPFP.
  //
  // Prefer angles with genuinely clean valleys. This prevents noisy fake
  // peak-count optima from beating the physically correct rotation.
  std::vector<int> valleyQualified;

  for (size_t i = 0; i < scan.size(); ++i) {
    const double valleyMetric = std::max(scan[i].meanValleyDrop, scan[i].meanGapDrop);

    const bool enoughPeaks = (scan[i].nPeaks >= 4);
    const bool enoughCleanGaps = (scan[i].nCleanGaps >= 4);
    const bool goodValleys = (valleyMetric >= 0.65);
    const bool decentGapQuality = (scan[i].meanGapQuality >= 0.35);

    if (enoughPeaks && enoughCleanGaps && goodValleys && decentGapQuality) {
      valleyQualified.push_back((int)i);
    }
  }

  int bestIdx = 0;

  auto betterTheta = [&](int ia, int ib) {
    const double va = std::max(scan[ia].meanValleyDrop, scan[ia].meanGapDrop);
    const double vb = std::max(scan[ib].meanValleyDrop, scan[ib].meanGapDrop);

    if (stableScores[ia] != stableScores[ib]) return stableScores[ia] > stableScores[ib];
    if (va != vb) return va > vb;
    if (scan[ia].nCleanGaps != scan[ib].nCleanGaps) return scan[ia].nCleanGaps > scan[ib].nCleanGaps;
    if (scan[ia].meanGapQuality != scan[ib].meanGapQuality) return scan[ia].meanGapQuality > scan[ib].meanGapQuality;
    if (scan[ia].nPeaks != scan[ib].nPeaks) return scan[ia].nPeaks > scan[ib].nPeaks;
    return scan[ia].score > scan[ib].score;
  };

  if (!valleyQualified.empty()) {
    bestIdx = valleyQualified[0];
    for (int idx : valleyQualified) {
      if (betterTheta(idx, bestIdx)) bestIdx = idx;
    }
    cout << "VALLEY-FIRST MODE: " << valleyQualified.size()
         << " theta candidates passed valley-quality cuts." << endl;
  } else {
    cout << "VALLEY-FIRST MODE: no theta passed valley-quality cuts; falling back to stable score." << endl;

    for (size_t i = 1; i < scan.size(); ++i) {
      if (betterTheta((int)i, bestIdx)) bestIdx = (int)i;
    }
  }

  // VALLEY-DOMINANT OVERRIDE WITH +/-1 DEGREE STABILITY.
  {
    const double minOverrideValleyDrop = 0.75;
    const double minOverrideGapQuality = 0.45;
    const double requiredDropImprovement = 0.10;

    const double minOneDegValleyFloor = 0.68;
    const double maxOneDegValleySpread = 0.18;

    int valleyBestIdx = -1;
    double valleyBestScore = -1.0e99;

    for (size_t i = 1; i + 1 < scan.size(); ++i) {
      if (scan[i].nPeaks < 4) continue;
      if (scan[i].nCleanGaps < 3) continue;
      if (scan[i].meanValleyDrop < minOverrideValleyDrop) continue;
      if (scan[i].meanGapQuality < minOverrideGapQuality) continue;

      const double vL = scan[i-1].meanValleyDrop;
      const double vC = scan[i].meanValleyDrop;
      const double vR = scan[i+1].meanValleyDrop;

      const double localFloor = std::min(vL, std::min(vC, vR));
      const double localCeil  = std::max(vL, std::max(vC, vR));
      const double localMean  = (vL + vC + vR) / 3.0;
      const double localSpread = localCeil - localFloor;

      if (localFloor < minOneDegValleyFloor) continue;
      if (localSpread > maxOneDegValleySpread) continue;

      const double vscore =
          1600.0 * localFloor
        +  700.0 * localMean
        +  250.0 * scan[i].meanGapQuality
        +   25.0 * scan[i].nCleanGaps
        +    3.0 * scan[i].nPeaks
        +    1.0e-6 * stableScores[i];

      if (vscore > valleyBestScore) {
        valleyBestScore = vscore;
        valleyBestIdx = (int)i;
      }
    }

    if (valleyBestIdx >= 0 &&
        scan[valleyBestIdx].meanValleyDrop >
          scan[bestIdx].meanValleyDrop + requiredDropImprovement) {

      cout << "VALLEY-DOMINANT OVERRIDE +/-1deg stable: theta "
           << scan[bestIdx].thetaDeg
           << " -> " << scan[valleyBestIdx].thetaDeg
           << " because meanValleyDrop improves from "
           << scan[bestIdx].meanValleyDrop
           << " to " << scan[valleyBestIdx].meanValleyDrop
           << endl;

      bestIdx = valleyBestIdx;
    }
  }

  // VALLEY PLATEAU CENTER SNAP.
  {
    double valleyMax = -1.0e99;

    for (size_t i = 1; i + 1 < scan.size(); ++i) {
      if (scan[i].nPeaks < 4) continue;
      if (scan[i].nCleanGaps < 3) continue;
      if (scan[i].meanGapQuality < 0.40) continue;
      if (scan[i].meanValleyDrop > valleyMax) valleyMax = scan[i].meanValleyDrop;
    }

    if (valleyMax > 0.75) {
      const double plateauCut = std::max(0.75, valleyMax - 0.08);

      int snapIdx = -1;
      double snapScoreBest = -1.0e99;

      for (size_t i = 1; i + 1 < scan.size(); ++i) {
        if (scan[i].nPeaks < 4) continue;
        if (scan[i].nCleanGaps < 3) continue;
        if (scan[i].meanGapQuality < 0.40) continue;
        if (scan[i].meanValleyDrop < plateauCut) continue;

        const double vL = scan[i-1].meanValleyDrop;
        const double vC = scan[i].meanValleyDrop;
        const double vR = scan[i+1].meanValleyDrop;

        double localFloor = vL;
        if (vC < localFloor) localFloor = vC;
        if (vR < localFloor) localFloor = vR;

        const double localMean = (vL + vC + vR) / 3.0;
        const double asym = (vL > vR) ? (vL - vR) : (vR - vL);

        const double snapScore =
            2500.0 * localFloor
          + 1000.0 * localMean
          +  250.0 * scan[i].meanGapQuality
          +   40.0 * scan[i].nCleanGaps
          +    4.0 * scan[i].nPeaks
          - 1200.0 * asym
          +    1.0e-6 * stableScores[i];

        if (snapScore > snapScoreBest) {
          snapScoreBest = snapScore;
          snapIdx = (int)i;
        }
      }

      if (snapIdx >= 0 && snapIdx != bestIdx) {
        cout << "VALLEY PLATEAU CENTER SNAP: theta "
             << scan[bestIdx].thetaDeg
             << " -> " << scan[snapIdx].thetaDeg
             << "  valleyMax=" << valleyMax
             << "  chosen meanValleyDrop=" << scan[snapIdx].meanValleyDrop
             << endl;

        bestIdx = snapIdx;
      }
    }
  }

  AngleResult best = scan[bestIdx];

  cout << "\n=== Stability-aware theta selection ===" << endl;
  cout << "Chosen theta: " << best.thetaDeg << " deg" << endl;
  cout << "Raw score at chosen theta: " << best.score << endl;
  cout << "Stable score at chosen theta: " << stableScores[bestIdx] << endl;

  int printLo = std::max<int>(0, bestIdx - 4);
  int printHi = std::min<int>((int)scan.size() - 1, bestIdx + 4);
  cout << "Local theta neighborhood:" << endl;
  cout << "  theta  rawScore  stableScore  nPeaks  cleanGaps  meanValleyDrop" << endl;
  for (int j = printLo; j <= printHi; ++j) {
    cout << "  " << scan[j].thetaDeg
         << "  " << scan[j].score
         << "  " << stableScores[j]
         << "  " << scan[j].nPeaks
         << "  " << scan[j].nCleanGaps
         << "  " << scan[j].meanValleyDrop;
    if (j == bestIdx) cout << "  <-- chosen";
    cout << endl;
  }

  // Save selected YFP/YPFP theta for later fixed-theta reuse.
  // UPSERT: same tag/run/foil/ndel is replaced.
  {
    TString thetaTsv = YFP_TSV_DIR + "/yfp_ypfp_selected_thetas.tsv";
    const TString header =
      "tag\trun\topticsID\tfoil\tndel\tdeltaMin\tdeltaMax\t"
      "thetaDeg\trawScore\tstableScore\tnPeaks\tnCleanGaps\t"
      "meanValleyDrop\tmeanGapQuality\tNevents";

    TString newLine = Form(
      "%s\t%d\t%s\t%d\t%d\t%g\t%g\t%g\t%g\t%g\t%d\t%d\t%g\t%g\t%lld",
      YfpOutputTag_angleScan(info).Data(),
      nrun,
      info.opticsID.Data(),
      foilIndex,
      ndelIndex,
      deltaMin,
      deltaMax,
      best.thetaDeg,
      best.score,
      stableScores[bestIdx],
      best.nPeaks,
      best.nCleanGaps,
      best.meanValleyDrop,
      best.meanGapQuality,
      (Long64_t)N
    );

    TString key = Form(
      "%s\t%d\t%d\t%d",
      YfpOutputTag_angleScan(info).Data(),
      nrun,
      foilIndex,
      ndelIndex
    );

    std::vector<TString> kept;
    bool replaced = false;

    std::ifstream ifs(thetaTsv.Data());
    if (ifs.good()) {
      std::string raw;
      while (std::getline(ifs, raw)) {
        TString line(raw.c_str());
        line = line.Strip(TString::kBoth);

        if (line.Length() == 0 || line.BeginsWith("tag\t")) continue;

        TObjArray *arr = line.Tokenize("\t");
        bool same = false;

        if (arr && arr->GetEntries() >= 5) {
          TString oldKey = Form(
            "%s\t%d\t%d\t%d",
            ((TObjString*)arr->At(0))->GetString().Data(),
            ((TObjString*)arr->At(1))->GetString().Atoi(),
            ((TObjString*)arr->At(3))->GetString().Atoi(),
            ((TObjString*)arr->At(4))->GetString().Atoi()
          );
          same = (oldKey == key);
        }

        if (arr) delete arr;

        if (same) {
          if (!replaced) {
            kept.push_back(newLine);
            replaced = true;
          }
        } else {
          kept.push_back(line);
        }
      }
      ifs.close();
    }

    if (!replaced) kept.push_back(newLine);

    std::ofstream ofs(thetaTsv.Data());
    ofs << header << "\n";
    for (auto &line : kept) ofs << line << "\n";
    ofs.close();

    cout << "Upserted selected Y theta TSV row: " << thetaTsv << endl;
    cout << "  key: " << key << endl;
    cout << "  thetaDeg now = " << best.thetaDeg << endl;
  }

  cout << "\n=== Best angle-scan result ===" << endl;
  cout << "Best theta: " << best.thetaDeg << " deg in standardized xz-yz space" << endl;
  cout << "Accepted peaks: " << best.nPeaks << endl;
  cout << "Score: " << best.score << endl;
  cout << "Mean valley drop: " << best.meanValleyDrop << endl;
  cout << "Clean gaps: " << best.nCleanGaps << " out of "
       << std::max(0, best.nPeaks - 1) << " adjacent peak gaps" << endl;
  cout << "Mean gap quality: " << best.meanGapQuality << endl;
  cout << "Min gap quality: " << best.minGapQuality << endl;
  cout << "Mean gap drop: " << best.meanGapDrop << endl;
  cout << "Mean gap empty: " << best.meanGapEmpty << endl;
  cout << "Fake peak penalty: " << best.fakePeakPenalty << endl;
  cout << "Total peak height: " << best.totalPeakHeight << endl;
  for (size_t i=0; i<best.peaks.size(); ++i) {
    cout << "  band " << i << ": q peak=" << best.peaks[i].q
         << ", height=" << best.peaks[i].height
         << ", prominence=" << best.peaks[i].prominence << endl;
  }
  cout << "Boundaries:";
  for (double b : best.bounds) cout << " " << b;
  cout << endl;

  // Compute final q,p coordinates using best theta.
  const double c = std::cos(best.thetaRad);
  const double s = std::sin(best.thetaRad);
  std::vector<double> qvals, pvals;
  qvals.reserve(N); pvals.reserve(N);
  for (Long64_t i=0; i<N; ++i) {
    const double q =  xzvals[i]*c + yzvals[i]*s;
    const double p = -xzvals[i]*s + yzvals[i]*c;
    qvals.push_back(q);
    pvals.push_back(p);
  }

  double xmin=*std::min_element(xvals.begin(), xvals.end());
  double xmax=*std::max_element(xvals.begin(), xvals.end());
  double ymin=*std::min_element(yvals.begin(), yvals.end());
  double ymax=*std::max_element(yvals.begin(), yvals.end());
  double xpad=0.10*(xmax-xmin), ypad=0.10*(ymax-ymin);
  xmin-=xpad; xmax+=xpad; ymin-=ypad; ymax+=ypad;

  double qmin=*std::min_element(qvals.begin(), qvals.end());
  double qmax=*std::max_element(qvals.begin(), qvals.end());
  double pmin=*std::min_element(pvals.begin(), pvals.end());
  double pmax=*std::max_element(pvals.begin(), pvals.end());
  double qpad=0.10*(qmax-qmin), ppad=0.10*(pmax-pmin);
  qmin-=qpad; qmax+=qpad; pmin-=ppad; pmax+=ppad;

  TH2D* hOrig = new TH2D("hOrig", Form("Run %d, %.1f < #delta < %.1f;ypfp;yfp", nrun, deltaMin, deltaMax),
                         220, xmin, xmax, 220, ymin, ymax);
  TH2D* hQP = new TH2D("hQP", Form("Best angle-scan coords, Run %d;q = best separation coordinate;p = perpendicular coordinate", nrun),
                       220, qmin, qmax, 220, pmin, pmax);
  TH1D* hQ = new TH1D("hQ", Form("Run %d;q = best separation coordinate;counts", nrun), 240, qmin, qmax);
  TH1D* hP = new TH1D("hP", Form("Run %d;p = perpendicular coordinate;counts", nrun), 220, pmin, pmax);

  for (Long64_t i=0; i<N; ++i) {
    hOrig->Fill(xvals[i], yvals[i]);
    hQP->Fill(qvals[i], pvals[i]);
    hQ->Fill(qvals[i]);
    hP->Fill(pvals[i]);
  }

  TH1D* hQSmooth = (TH1D*)hQ->Clone("hQSmooth");
  for (int i=0; i<smoothPasses; ++i) hQSmooth->Smooth(1);

  // -------------------------------------------------------------------------
  // Final band-center selection:
  //   * The angle scan above uses ONLY strong/prominent peaks. This protects
  //     the rotation from tiny edge ripples.
  //   * After theta is fixed, we allow weak-but-separated peaks to become
  //     their own diagnostic bands. These are tagged as weak so a later stage
  //     can keep, reject, or down-weight them.
  // -------------------------------------------------------------------------
  std::vector<QPeak> bandPeaks = best.peaks;
  std::vector<int> bandIsWeak(bandPeaks.size(), 0);

  struct LocalQMax {
    double q = 0.0;
    double height = 0.0;
    double prominence = 0.0;
    int bin = 0;
  };

  std::vector<LocalQMax> localMaximaFinal;
  const double finalMaxContent = hQSmooth->GetMaximum();
  const double weakMinHeight = minWeakPeakFraction * finalMaxContent;

  for (int ib=2; ib<hQSmooth->GetNbinsX(); ++ib) {
    const double ym = hQSmooth->GetBinContent(ib-1);
    const double y0 = hQSmooth->GetBinContent(ib);
    const double yp = hQSmooth->GetBinContent(ib+1);
    if (y0 > ym && y0 >= yp && y0 > weakMinHeight) {
      localMaximaFinal.push_back({hQSmooth->GetBinCenter(ib), y0, 0.0, ib});
    }
  }

  std::sort(localMaximaFinal.begin(), localMaximaFinal.end(),
            [](const LocalQMax& a, const LocalQMax& b){ return a.q < b.q; });

  std::vector<QPeak> weakCandidates;
  for (size_t i=0; i<localMaximaFinal.size(); ++i) {
    const double h0 = localMaximaFinal[i].height;
    double leftProm = 1.0;
    double rightProm = 1.0;

    if (i > 0) {
      const double valley = LocalValleyMin_angleScan(hQSmooth, localMaximaFinal[i-1].bin, localMaximaFinal[i].bin);
      leftProm = (h0 > 0.0) ? (h0 - valley)/h0 : 0.0;
    }

    if (i+1 < localMaximaFinal.size()) {
      const double valley = LocalValleyMin_angleScan(hQSmooth, localMaximaFinal[i].bin, localMaximaFinal[i+1].bin);
      rightProm = (h0 > 0.0) ? (h0 - valley)/h0 : 0.0;
    }

    double prom = 0.0;
    if (i == 0 && localMaximaFinal.size() > 1) prom = rightProm;
    else if (i+1 == localMaximaFinal.size() && localMaximaFinal.size() > 1) prom = leftProm;
    else if (localMaximaFinal.size() == 1) prom = 1.0;
    else prom = std::min(leftProm, rightProm);

    if (prom < minWeakProminenceFraction) continue;

    QPeak pk{localMaximaFinal[i].q, h0, prom, localMaximaFinal[i].bin};

    // Skip if this weak candidate is already represented by a strong peak.
    bool alreadyStrong = false;
    for (const auto& strong : best.peaks) {
      if (std::fabs(pk.q - strong.q) < weakMinSepFactor * minPeakSep) {
        alreadyStrong = true;
        break;
      }
    }
    if (!alreadyStrong) weakCandidates.push_back(pk);
  }

  // Add weak candidates by prominence/height, without allowing tiny duplicate ripples.
  std::sort(weakCandidates.begin(), weakCandidates.end(),
            [](const QPeak& a, const QPeak& b){
              if (a.prominence != b.prominence) return a.prominence > b.prominence;
              return a.height > b.height;
            });

  int nWeakAdded = 0;
  for (const auto& cand : weakCandidates) {
    if ((int)bandPeaks.size() >= maxBands) break;

    bool tooClose = false;
    for (const auto& pk : bandPeaks) {
      if (std::fabs(cand.q - pk.q) < weakMinSepFactor * minPeakSep) {
        tooClose = true;
        break;
      }
    }
    if (tooClose) continue;

    bandPeaks.push_back(cand);
    bandIsWeak.push_back(1);
    nWeakAdded++;
  }

  // Sort final peak list and carry weak/strong labels along with it.
  std::vector<size_t> order(bandPeaks.size());
  for (size_t i=0; i<order.size(); ++i) order[i]=i;
  std::sort(order.begin(), order.end(),
            [&](size_t a, size_t b){ return bandPeaks[a].q < bandPeaks[b].q; });

  std::vector<QPeak> sortedBandPeaks;
  std::vector<int> sortedBandIsWeak;
  for (size_t idx : order) {
    sortedBandPeaks.push_back(bandPeaks[idx]);
    sortedBandIsWeak.push_back(bandIsWeak[idx]);
  }
  bandPeaks.swap(sortedBandPeaks);
  bandIsWeak.swap(sortedBandIsWeak);

  std::vector<double> bandBounds = BoundariesFromPeaks_angleScan(bandPeaks);

  cout << "\n=== Final strong+weak q-band centers ===" << endl;
  cout << "Strong peaks from angle scan: " << best.peaks.size() << endl;
  cout << "Weak candidates considered after fixed angle: " << weakCandidates.size() << endl;
  cout << "Weak peaks added after fixed angle: " << nWeakAdded << endl;
  cout << "Final diagnostic bands: " << bandPeaks.size() << endl;
  for (size_t i=0; i<bandPeaks.size(); ++i) {
    cout << "  qband " << i
         << " [" << (bandIsWeak[i] ? "weak" : "strong") << "]"
         << ": peak q=" << bandPeaks[i].q
         << ", height=" << bandPeaks[i].height
         << ", prominence=" << bandPeaks[i].prominence << endl;
  }

  std::vector<int> bandIndex(N, -1);
  std::vector<int> bandCounts(bandPeaks.size(), 0);
  for (Long64_t i=0; i<N; ++i) {
    int b=0;
    while (b < (int)bandBounds.size() && qvals[i] > bandBounds[b]) b++;
    if (b >= 0 && b < (int)bandPeaks.size()) {
      bandIndex[i]=b;
      bandCounts[b]++;
    }
  }

  cout << "\n=== Final q-band assignment ===" << endl;
  for (size_t i=0; i<bandPeaks.size(); ++i) {
    cout << "  qband " << i << " [" << (bandIsWeak[i] ? "weak" : "strong") << "]: peak q=" << bandPeaks[i].q
         << ", assigned events=" << bandCounts[i] << endl;
  }

  // Diagnostic output paths were built earlier by BuildYfpDiagnosticBase_angleScan(...).

  std::vector<int> colors = {kRed+1, kOrange+7, kSpring+5, kGreen+2, kAzure+7, kBlue+1, kMagenta+1, kViolet+7, kCyan+2, kGray+2};
  std::vector<TGraph*> graphs(bandPeaks.size(), nullptr);
  for (size_t b=0; b<bandPeaks.size(); ++b) graphs[b] = new TGraph();

  for (Long64_t i=0; i<N; ++i) {
    int b = bandIndex[i];
    if (b < 0 || b >= (int)graphs.size()) continue;
    int pnt = graphs[b]->GetN();
    graphs[b]->SetPoint(pnt, xvals[i], yvals[i]);
  }

  std::vector<double> bandLow(bandPeaks.size(), qmin), bandHigh(bandPeaks.size(), qmax);
  for (size_t b=0; b<bandPeaks.size(); ++b) {
    bandLow[b]  = (b==0) ? qmin : bandBounds[b-1];
    bandHigh[b] = (b+1==bandPeaks.size()) ? qmax : bandBounds[b];
  }

  auto drawConstQLine = [&](double q0, int color, int style, int width) {
    // In standardized space, xz = q*c - p*s, yz = q*s + p*c.
    // Map back to raw x,y.
    const double pA = pmin;
    const double pB = pmax;
    const double xzA = q0*c - pA*s;
    const double yzA = q0*s + pA*c;
    const double xzB = q0*c - pB*s;
    const double yzB = q0*s + pB*c;
    TLine* line = new TLine(mx + sx*xzA, my + sy*yzA,
                            mx + sx*xzB, my + sy*yzB);
    line->SetLineColor(color);
    line->SetLineStyle(style);
    line->SetLineWidth(width);
    line->Draw("same");
  };

  TCanvas* c1 = new TCanvas("c_angleScanBand_assign", "Angle-scan q-band assignment", 1200, 950);
  c1->Divide(2,2);

  c1->cd(1);
  gPad->SetLogz();
  hOrig->Draw("COLZ");
  for (const auto& pk : bandPeaks) drawConstQLine(pk.q, kRed+1, 1, 2);
  for (double b : bandBounds) drawConstQLine(b, kRed+1, 2, 2);
  TLatex lat;
  lat.SetNDC(); lat.SetTextSize(0.030);
  lat.DrawLatex(0.12,0.92,Form("Best #theta = %.2f deg z-space; bands=%zu, weak=%d", best.thetaDeg, bandPeaks.size(), nWeakAdded));

  c1->cd(2);
  gPad->SetLogz();
  hQP->Draw("COLZ");
  for (const auto& pk : bandPeaks) {
    TLine* lpk = new TLine(pk.q, pmin, pk.q, pmax);
    lpk->SetLineColor(kRed+1); lpk->SetLineWidth(2); lpk->Draw("same");
  }
  for (double b : bandBounds) {
    TLine* lb = new TLine(b, pmin, b, pmax);
    lb->SetLineColor(kRed+1); lb->SetLineStyle(2); lb->SetLineWidth(2); lb->Draw("same");
  }

  c1->cd(3);
  hQ->Draw("HIST");
  hQSmooth->SetLineColor(kRed+1);
  hQSmooth->SetLineWidth(2);
  hQSmooth->Draw("HIST SAME");
  for (const auto& pk : bandPeaks) {
    TLine* lpk = new TLine(pk.q, 0, pk.q, hQ->GetMaximum()*1.05);
    lpk->SetLineColor(kRed+1); lpk->SetLineWidth(2); lpk->Draw("same");
  }
  for (double b : bandBounds) {
    TLine* lb = new TLine(b, 0, b, hQ->GetMaximum()*1.05);
    lb->SetLineColor(kRed+1); lb->SetLineStyle(2); lb->Draw("same");
  }

  c1->cd(4);
  TH2D* hFrame = new TH2D("hFrameAngleBands", Form("angle-scan q-band assignments, Run %d;ypfp;yfp", nrun), 10, xmin, xmax, 10, ymin, ymax);
  hFrame->SetMinimum(0); hFrame->SetMaximum(1);
  hFrame->Draw("AXIS");
  for (size_t b=0; b<graphs.size(); ++b) {
    graphs[b]->SetMarkerStyle(20);
    graphs[b]->SetMarkerSize(0.18);
    graphs[b]->SetMarkerColor(colors[b % colors.size()]);
    graphs[b]->Draw("P SAME");
  }

  // Score-vs-angle page
  TGraph* gScore = new TGraph();
  TGraph* gStableScore = new TGraph();
  TGraph* gNPeaks = new TGraph();
  TGraph* gValley = new TGraph();
  for (size_t i=0; i<scan.size(); ++i) {
    gScore->SetPoint(i, scan[i].thetaDeg, scan[i].score);
    gStableScore->SetPoint(i, scan[i].thetaDeg, stableScores[i]);
    gNPeaks->SetPoint(i, scan[i].thetaDeg, scan[i].nPeaks);
    gValley->SetPoint(i, scan[i].thetaDeg, scan[i].meanValleyDrop);
  }

  c1->SaveAs(outpdf + "[");
  c1->SaveAs(outpdf);

  TCanvas* cScore = new TCanvas("c_angle_scan_score", "Angle scan score", 1200, 900);
  cScore->Divide(1,3);

  cScore->cd(1);
  gScore->SetTitle(Form("Angle scan score, Run %d, %.1f<#delta<%.1f;theta (deg);score", nrun, deltaMin, deltaMax));
  gScore->SetLineWidth(2);
  gScore->Draw("AL");
  gStableScore->SetLineColor(kRed+1);
  gStableScore->SetLineWidth(2);
  gStableScore->Draw("L SAME");
  TLine* bestScoreLine = new TLine(best.thetaDeg, gPad->GetUymin(), best.thetaDeg, gPad->GetUymax());
  bestScoreLine->SetLineColor(kBlue+1); bestScoreLine->SetLineWidth(2); bestScoreLine->Draw("same");
  TLatex stableLat;
  stableLat.SetNDC();
  stableLat.SetTextSize(0.030);
  stableLat.DrawLatex(0.14,0.84,"black/raw score; red/stability-adjusted score");

  cScore->cd(2);
  gNPeaks->SetTitle("Accepted q peaks vs angle;theta (deg);accepted peaks");
  gNPeaks->SetLineWidth(2);
  gNPeaks->Draw("AL");
  TLine* bestPeakLine = new TLine(best.thetaDeg, gPad->GetUymin(), best.thetaDeg, gPad->GetUymax());
  bestPeakLine->SetLineColor(kRed+1); bestPeakLine->SetLineWidth(2); bestPeakLine->Draw("same");

  cScore->cd(3);
  gValley->SetTitle("Mean valley drop vs angle;theta (deg);mean valley drop");
  gValley->SetLineWidth(2);
  gValley->Draw("AL");
  TLine* bestValleyLine = new TLine(best.thetaDeg, gPad->GetUymin(), best.thetaDeg, gPad->GetUymax());
  bestValleyLine->SetLineColor(kRed+1); bestValleyLine->SetLineWidth(2); bestValleyLine->Draw("same");

  cScore->SaveAs(outpdf);

  std::vector<TH2D*> hBandOrig(bandPeaks.size(), nullptr);
  std::vector<TH2D*> hBandQP(bandPeaks.size(), nullptr);
  std::vector<TH1D*> hBandQ(bandPeaks.size(), nullptr);

  for (size_t b=0; b<bandPeaks.size(); ++b) {
    hBandOrig[b] = new TH2D(Form("hBandOrig_%zu", b), Form("Band %zu only in original coords, Run %d;ypfp;yfp", b, nrun), 220, xmin, xmax, 220, ymin, ymax);
    hBandQP[b] = new TH2D(Form("hBandQP_%zu", b), Form("Band %zu only in q-p coords, Run %d;q;p", b, nrun), 220, qmin, qmax, 220, pmin, pmax);
    hBandQ[b] = new TH1D(Form("hBandQ_%zu", b), Form("Band %zu q distribution, Run %d;q;counts", b, nrun), 240, qmin, qmax);
  }

  for (Long64_t i=0; i<N; ++i) {
    int b = bandIndex[i];
    if (b < 0 || b >= (int)bandPeaks.size()) continue;
    hBandOrig[b]->Fill(xvals[i], yvals[i]);
    hBandQP[b]->Fill(qvals[i], pvals[i]);
    hBandQ[b]->Fill(qvals[i]);
  }

  for (size_t b=0; b<bandPeaks.size(); ++b) {
    TCanvas* cb = new TCanvas(Form("c_angle_band_%zu", b), Form("Angle-scan band %zu diagnostics", b), 1200, 950);
    cb->Divide(2,2);

    cb->cd(1);
    gPad->SetLogz();
    hOrig->Draw("COLZ");
    TGraph* gall = new TGraph();
    TGraph* gsel = new TGraph();
    for (Long64_t i=0; i<N; ++i) {
      int pnt = gall->GetN();
      gall->SetPoint(pnt, xvals[i], yvals[i]);
      if (bandIndex[i] == (int)b) {
        int qn = gsel->GetN();
        gsel->SetPoint(qn, xvals[i], yvals[i]);
      }
    }
    gall->SetMarkerStyle(20); gall->SetMarkerSize(0.10); gall->SetMarkerColor(kGray+1); gall->Draw("P SAME");
    gsel->SetMarkerStyle(20); gsel->SetMarkerSize(0.22); gsel->SetMarkerColor(colors[b % colors.size()]); gsel->Draw("P SAME");
    for (const auto& pk : bandPeaks) drawConstQLine(pk.q, kRed+1, 2, 1);
    drawConstQLine(bandPeaks[b].q, colors[b % colors.size()], 1, 3);
    drawConstQLine(bandLow[b], colors[b % colors.size()], 2, 3);
    drawConstQLine(bandHigh[b], colors[b % colors.size()], 2, 3);
    lat.DrawLatex(0.12,0.92,Form("Angle band %zu [%s] highlighted in original coords", b, bandIsWeak[b] ? "weak" : "strong"));
    lat.DrawLatex(0.12,0.88,Form("peak q=%.3f, range [%.3f, %.3f], N=%d", bandPeaks[b].q, bandLow[b], bandHigh[b], bandCounts[b]));

    cb->cd(2);
    gPad->SetLogz();
    hQP->Draw("COLZ");
    TGraph* gallqp = new TGraph();
    TGraph* gselqp = new TGraph();
    for (Long64_t i=0; i<N; ++i) {
      int pnt = gallqp->GetN();
      gallqp->SetPoint(pnt, qvals[i], pvals[i]);
      if (bandIndex[i] == (int)b) {
        int qn = gselqp->GetN();
        gselqp->SetPoint(qn, qvals[i], pvals[i]);
      }
    }
    gallqp->SetMarkerStyle(20); gallqp->SetMarkerSize(0.10); gallqp->SetMarkerColor(kGray+1); gallqp->Draw("P SAME");
    gselqp->SetMarkerStyle(20); gselqp->SetMarkerSize(0.22); gselqp->SetMarkerColor(colors[b % colors.size()]); gselqp->Draw("P SAME");
    TLine* qlo = new TLine(bandLow[b], pmin, bandLow[b], pmax);
    TLine* qhi = new TLine(bandHigh[b], pmin, bandHigh[b], pmax);
    qlo->SetLineColor(colors[b % colors.size()]); qhi->SetLineColor(colors[b % colors.size()]);
    qlo->SetLineStyle(2); qhi->SetLineStyle(2); qlo->SetLineWidth(3); qhi->SetLineWidth(3);
    qlo->Draw("same"); qhi->Draw("same");

    cb->cd(3);
    hQ->Draw("HIST");
    hQSmooth->SetLineColor(kRed+1); hQSmooth->SetLineWidth(2); hQSmooth->Draw("HIST SAME");
    double ymaxQ = hQ->GetMaximum()*1.08;
    for (const auto& pk : bandPeaks) {
      TLine* lpk = new TLine(pk.q, 0, pk.q, ymaxQ);
      lpk->SetLineColor(kRed+1); lpk->SetLineWidth(pk.q == bandPeaks[b].q ? 3 : 1); lpk->SetLineStyle(pk.q == bandPeaks[b].q ? 1 : 2); lpk->Draw("same");
    }
    TLine* bq1 = new TLine(bandLow[b], 0, bandLow[b], ymaxQ);
    TLine* bq2 = new TLine(bandHigh[b], 0, bandHigh[b], ymaxQ);
    bq1->SetLineColor(colors[b % colors.size()]); bq2->SetLineColor(colors[b % colors.size()]);
    bq1->SetLineWidth(3); bq2->SetLineWidth(3); bq1->SetLineStyle(2); bq2->SetLineStyle(2);
    bq1->Draw("same"); bq2->Draw("same");
    hBandQ[b]->SetLineColor(colors[b % colors.size()]); hBandQ[b]->SetLineWidth(2); hBandQ[b]->Draw("HIST SAME");

    cb->cd(4);
    gPad->SetLogz();
    hBandOrig[b]->Draw("COLZ");
    for (const auto& pk : bandPeaks) drawConstQLine(pk.q, kRed+1, 2, 1);
    drawConstQLine(bandPeaks[b].q, colors[b % colors.size()], 1, 3);
    TLatex lat2; lat2.SetNDC(); lat2.SetTextSize(0.035);
    lat2.DrawLatex(0.12,0.92,Form("Angle band %zu [%s] only", b, bandIsWeak[b] ? "weak" : "strong"));
    lat2.DrawLatex(0.12,0.87,Form("N=%d events", bandCounts[b]));
    lat2.DrawLatex(0.12,0.82,Form("q peak=%.3f", bandPeaks[b].q));
    lat2.DrawLatex(0.12,0.77,Form("q range [%.3f, %.3f]", bandLow[b], bandHigh[b]));

    cb->SaveAs(outpdf);
  }

  c1->SaveAs(outpdf + "]");


  // ----------------------------
  // Hand-compatible AUTO BAND TCutG output
  // ----------------------------
  TFile fAutoCuts(outCutRootAutoBand, "UPDATE");
  if (fAutoCuts.IsZombie()) {
    cout << "ERROR: could not open auto band cut output file: "
         << outCutRootAutoBand << endl;
  } else {
    // Delete all existing cuts for this exact foil/delta slice before writing.
    //
    // This matters when a rerun produces fewer bands than the previous attempt:
    // kOverwrite only replaces same-name objects. It does not remove old
    // yscol objects that are no longer produced.
    std::vector<TString> oldSliceCutNames;
    TString sliceTag = Form("_nfoil_%d_ndel_%d", foilIndex, ndelIndex);

    TIter keyIter(fAutoCuts.GetListOfKeys());
    TKey* oldKey = nullptr;
    while ((oldKey = (TKey*)keyIter())) {
      TString oldName = oldKey->GetName();

      if (!oldName.BeginsWith("hYpFpYFp_cut_yscol_")) continue;
      if (!oldName.Contains(sliceTag)) continue;

      if (std::find(oldSliceCutNames.begin(),
                    oldSliceCutNames.end(),
                    oldName) == oldSliceCutNames.end()) {
        oldSliceCutNames.push_back(oldName);
      }
    }

    if (!oldSliceCutNames.empty()) {
      cout << "Deleting old AUTO BAND cuts for this slice before rewrite:" << endl;
      for (const auto& oldName : oldSliceCutNames) {
        cout << "  deleting " << oldName << ";*" << endl;
        fAutoCuts.Delete(oldName + ";*");
      }
    } else {
      cout << "No pre-existing AUTO BAND cuts found for this foil/delta slice." << endl;
    }

    cout << "\nWriting hand-compatible AUTO BAND YFP/YPFP TCutG cuts:" << endl;

    for (size_t b = 0; b < bandPeaks.size(); ++b) {
      Int_t yscol = YscolFromQBand_angleScan((Int_t)b, (Int_t)bandPeaks.size());

      if (yscol < 0 || yscol > 8) {
        cout << "WARNING: skipping qband " << b
             << " because mapped yscol=" << yscol
             << " is outside [0,8]" << endl;
        continue;
      }

      TString cname = Form("hYpFpYFp_cut_yscol_%d_nfoil_%d_ndel_%d",
                           yscol, foilIndex, ndelIndex);

      TCutG* cband = MakeQBandCut_angleScan(cname,
                                             bandLow[b],
                                             bandHigh[b],
                                             pmin,
                                             pmax,
                                             mx,
                                             my,
                                             sx,
                                             sy,
                                             c,
                                             s);

      fAutoCuts.cd();
      cband->Write(cname, TObject::kOverwrite);

      cout << "  wrote " << cname
           << "  from qband=" << b
           << "  qrange=[" << bandLow[b] << ", " << bandHigh[b] << "]"
           << "  N=" << bandCounts[b]
           << "  weak=" << bandIsWeak[b]
           << endl;

      delete cband;
    }

    fAutoCuts.Close();
    cout << "Wrote hand-compatible AUTO BAND cut file: "
         << outCutRootAutoBand << endl;
  }

  TFile fout(outroot, "RECREATE");
  hOrig->Write();
  hQP->Write();
  hQ->Write();
  hP->Write();
  hQSmooth->Write();
  gScore->Write("g_angle_score");
  gStableScore->Write("g_angle_stable_score");
  gNPeaks->Write("g_angle_npeaks");
  gValley->Write("g_angle_valley");
  for (size_t b=0; b<graphs.size(); ++b) graphs[b]->Write(Form("g_qband_%zu", b));
  for (size_t b=0; b<bandPeaks.size(); ++b) {
    hBandOrig[b]->Write();
    hBandQP[b]->Write();
    hBandQ[b]->Write();
  }
  fout.WriteObject(c1, "c_angleScanBand_assign_overview");
  fout.WriteObject(cScore, "c_angle_scan_score");
  fout.Close();

  cout << "Wrote: " << outpdf << endl;
  cout << "Wrote: " << outroot << endl;
}


// ============================================================================
// New wrapper main function added by patch_yfp_cleanGaps_batch_ndel_ONLY.py
//
// Default mode:
//   process one delta slice: deltaMin to deltaMax
//
// Batch mode:
//   set runAllDeltaSlices=true to process:
//   -10,-8,-5,0,5,10
//
// ndelOverride should usually be left at -999.
// ============================================================================
void assign_yfp_ypfp_angleScanBands(Int_t nrun=1544,
                     Double_t deltaMin=-8.0,
                     Double_t deltaMax=-5.0,
                     const char *ytarTag="ML_dev",
                     Int_t maxBands=9,
                     Double_t thetaStepDeg=1.0,
                     Double_t minPeakSep=0.18,
                     Double_t minPeakFrac=0.08,
                     Double_t minPromFrac=0.25,
                     Int_t smoothPasses=2,
                     Double_t weakMinPeakFrac=0.005,
                     Double_t weakMinPromFrac=0.08,
                     Double_t weakMinSepFactor=1.0,
                     Bool_t useYtarCut=true,
                     Int_t foilIndex=0,
                     Long64_t maxEvents=-1,
                     Int_t FileID=-1,
                     Int_t ndelOverride=-999,
                     Bool_t runAllDeltaSlices=false,
                     Double_t fixedThetaDeg=-999.0,
                     const char *campaignDir="HMS_6p117",
                     const char *inputRootOverride="")
{
  if (!runAllDeltaSlices) {
    run_yfp_ypfp_angleScan_oneDelta(nrun,
                  deltaMin,
                  deltaMax,
                  ytarTag,
                  maxBands,
                  thetaStepDeg,
                  minPeakSep,
                  minPeakFrac,
                  minPromFrac,
                  smoothPasses,
                  weakMinPeakFrac,
                  weakMinPromFrac,
                  weakMinSepFactor,
                  useYtarCut,
                  foilIndex,
                  maxEvents,
                  FileID,
                  ndelOverride,
                  fixedThetaDeg,
                  campaignDir,
                  inputRootOverride);
    return;
  }

  cout << endl;
  cout << "============================================================" << endl;
  cout << "BATCH MODE: YFP/YPFP angle scan over standard delta slices" << endl;
  cout << "Run: " << nrun << "  foilIndex: " << foilIndex << endl;
  cout << "Delta edges: -10,-8,-5,0,5,10" << endl;
  cout << "ytarTag: " << ytarTag << endl;
  cout << "============================================================" << endl;

  const int nEdges = 6;
  const double edges[nEdges] = {-10.0, -8.0, -5.0, 0.0, 5.0, 10.0};

  for (int i = 0; i < nEdges - 1; ++i) {
    cout << endl;
    cout << "------------------------------------------------------------" << endl;
    cout << "Batch delta slice " << i
         << ": [" << edges[i] << ", " << edges[i+1] << ")"
         << endl;
    cout << "------------------------------------------------------------" << endl;

    run_yfp_ypfp_angleScan_oneDelta(nrun,
                  edges[i],
                  edges[i+1],
                  ytarTag,
                  maxBands,
                  thetaStepDeg,
                  minPeakSep,
                  minPeakFrac,
                  minPromFrac,
                  smoothPasses,
                  weakMinPeakFrac,
                  weakMinPromFrac,
                  weakMinSepFactor,
                  useYtarCut,
                  foilIndex,
                  maxEvents,
                  FileID,
                  i,
                  fixedThetaDeg,
                  campaignDir,
                  inputRootOverride);
  }

  cout << endl;
  cout << "BATCH MODE DONE for run " << nrun
       << ", foilIndex " << foilIndex << endl;
}


