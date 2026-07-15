// relabel_xfp_xpfp_autoBands_to_xscol_coloredDensity_batch.C
//
// XFP/XPFP companion to relabel_yfp_ypfp_autoBands_to_yscol_coloredDensity_batch.C.
//
// Purpose:
//   Post-process XFP/XPFP auto-band cuts whose qband/candidate labels are only
//   provisional.  For each cut, project selected events to HMS sieve space,
//   assign the physical xscol from the dominant sieve-row population, and write
//   relabeled component cuts.  Using the mean xsieve alone is unsafe for
//   XFP/XPFP because low/high zone candidates can have tails/overlap that
//   bias the mean into the wrong adjacent row.
//
// Important XFP/XPFP difference from the YFP/YPFP relabeler:
//   The current XFP auto-band step may write separate low/high xpfp-zone
//   candidate cuts, e.g.
//     hXpFpXFp_cut_candidate_02_zone_low_nfoil_0_ndel_3
//     hXpFpXFp_cut_candidate_01_zone_high_nfoil_0_ndel_3
//
//   If two candidate cuts map to the same physical xscol, this macro keeps BOTH.
//   ROOT TCutG cannot safely store a disconnected union as one polygon, so the
//   output uses component names:
//     hXpFpXFp_cut_xscol_4_nfoil_0_ndel_3_part00_zone_low_src_02
//     hXpFpXFp_cut_xscol_4_nfoil_0_ndel_3_part01_zone_high_src_01
//
//   Downstream XFP candidate-tree/GMM code should treat all cuts with the same
//   xscol/nfoil/ndel prefix as an OR/additive collection.
//
// Diagnostics:
//   - low-zone colored-density page
//   - high-zone colored-density page
//   - all-component overlay page
//   - final combined-by-xscol pages showing the summed low+high population
//
// Usage examples:
//   // If the auto-band step used low/high raw xpfp gates, pass the split as
//   // the final argument so relabel counts match the auto-band diagnostics.
//   hcana -b -l -q 'relabel_xfp_xpfp_autoBands_to_xscol_coloredDensity_batch.C(1544,0,3,"ML_dev",-1,"",true,-1,6.0,0.65,25,true,false,"all",false,<xpfpSplit>)'
//
//   // all delta slices for one foil, preserving already relabeled foils if the
//   // final relabeled file exists:
//   hcana -b -l -q 'relabel_xfp_xpfp_autoBands_to_xscol_coloredDensity_batch.C(1544,0,0,"ML_dev",-1,"",true,-1,6.0,0.65,25,true,true,"all")'
//
// Optional override file format:
//   # run foil ndel zone oldIndex newXscol action comment
//   1544 0 3 low  2 4 force good_low_component
//   1544 0 3 high 1 4 force good_high_component_same_xscol
//   1544 0 3 low  5 -1 reject bad_fragment
//
// Backward-compatible override format is also accepted:
//   # run foil ndel oldIndex newXscol action comment
//   1544 0 3 2 4 force old_format_applies_to_any_zone

#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TCutG.h>
#include <TKey.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TText.h>
#include <TLatex.h>
#include <TBox.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TROOT.h>
#include <TColor.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TMath.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <set>
#include <cmath>
#include <algorithm>
#include <limits>

using std::cout;
using std::endl;

// ============================================================================
// Paths. Keep these matched to the current XFP/XPFP angle-scan macro.
// ============================================================================
const TString XRELABEL_REPO_DIR  = "/u/group/nps/cploen/ML_HMS_OPTICS";
const TString XRELABEL_DAT_DIR   = XRELABEL_REPO_DIR + "/DATfiles";
const TString XRELABEL_CUTS_DIR  = XRELABEL_REPO_DIR + "/cuts";
const TString XRELABEL_PLOTS_DIR = XRELABEL_REPO_DIR + "/plots";

const TString XRELABEL_TREE_ROOT_DIR =
  "/volatile/hallc/nps/cploen/ROOTfiles/OPTICS/angular_sandbox/newfit_6p667_20260526_1226_no_offsets/rootfiles/";

// Optional internal overrides used by batch mode so each delta slice can build
// on the relabeled output from the previous slice. Leave empty in normal use.
static TString XRELABEL_INPUT_CUT_OVERRIDE = "";
static TString XRELABEL_OUTPUT_CUT_OVERRIDE = "";

// Campaign-aware overrides. These are set by the shell-generated chain macro.
static TString XRELABEL_RUNGROUP_TAG = "";
static TString XRELABEL_CAMPAIGN_DIR = "";
static TString XRELABEL_INPUT_TREE_OVERRIDE = "";

// ============================================================================

struct XRelabelRunInfo {
  int run = -1;
  TString opticsID = "";
  double centAngle = 0.0;
  int numFoil = 0;
  int sieveFlag = 0;
  int ndelcut = 0;
  std::vector<double> zfoil;
  std::vector<double> delcut;
};

struct XBandCutInfo {
  // oldIndex means candidate index for zone candidates, or old provisional
  // xscol for older global auto-band cuts.
  int oldIndex = -1;
  int oldXscol = -1;
  int candidateIndex = -1;
  int newXscol = -1;
  int forcedXscol = -999;
  int partIndex = -1;

  bool reject = false;
  bool hasOverride = false;
  bool write = false;

  TString action = "auto";
  TString zone = "global";
  TString sourceKind = "unknown"; // candidate or xscol
  TString oldName = "";
  TString newName = "";

  TCutG* cut = nullptr;
  TH2D* hYsXs = nullptr;
  TH2D* hXpX = nullptr;

  Long64_t n = 0;
  double meanXs = 0.0;
  double meanYs = 0.0;
  double rmsXs = 0.0;
  double rmsYs = 0.0;
  double distToGuide = 0.0;

  // Dominant physical sieve-row assignment diagnostics.  This is the primary
  // assignment metric for XFP/XPFP; meanXs is kept only as a diagnostic.
  std::vector<Long64_t> xscolCounts;
  int modalXscol = -1;
  Long64_t modalCount = 0;
  double modalFrac = 0.0;
};

static std::vector<TString> SplitCSV_xrelabel(const TString& line)
{
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

static bool ReadOpticsRunInfo_xrelabel(int nrun, XRelabelRunInfo& info,
                                       const char* metaFile)
{
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

    auto tok = SplitCSV_xrelabel(line);
    if (tok.size() < 6) continue;
    if (tok[0].Atoi() != nrun) continue;

    info.run       = tok[0].Atoi();
    info.opticsID  = tok[1];
    info.centAngle = tok[2].Atof();
    info.numFoil   = tok[3].Atoi();
    info.sieveFlag = tok[4].Atoi();
    info.ndelcut   = tok[5].Atoi();

    std::string zline_raw;
    if (!std::getline(fin, zline_raw)) return false;
    auto ztok = SplitCSV_xrelabel(TString(zline_raw.c_str()).Strip(TString::kBoth));
    for (int i = 0; i < info.numFoil && i < (int)ztok.size(); ++i)
      info.zfoil.push_back(ztok[i].Atof());

    std::string dline_raw;
    if (!std::getline(fin, dline_raw)) return false;
    auto dtok = SplitCSV_xrelabel(TString(dline_raw.c_str()).Strip(TString::kBoth));
    for (int i = 0; i < (int)dtok.size(); ++i)
      info.delcut.push_back(dtok[i].Atof());

    cout << "Parsed run " << nrun
         << " OpticsID=" << info.opticsID
         << " NumFoil=" << info.numFoil
         << " ndelcut=" << info.ndelcut
         << " delta edges=" << info.delcut.size()
         << endl;
    return true;
  }

  cout << "ERROR: run " << nrun << " not found in " << metaFile << endl;
  return false;
}

static TString TreeRootPath_xrelabel(int nrun, int FileID)
{
  if (XRELABEL_INPUT_TREE_OVERRIDE.Length() > 0)
    return XRELABEL_INPUT_TREE_OVERRIDE;

  return Form("%s/nps_hms_optics_%d_1_%d.root",
              XRELABEL_TREE_ROOT_DIR.Data(), nrun, FileID);
}

static TString AutoBandCutFile_xrelabel(const XRelabelRunInfo& info, int FileID)
{
  if (XRELABEL_CAMPAIGN_DIR.Length() > 0 &&
      XRELABEL_RUNGROUP_TAG.Length() > 0) {
    return Form("%s/02b_angle_scan_x/cuts/XpFpXFp_%s_auto_band_cut.root",
                XRELABEL_CAMPAIGN_DIR.Data(),
                XRELABEL_RUNGROUP_TAG.Data());
  }

  return Form("%s/XpFpXFp_%s_%d_auto_band_cut.root",
              XRELABEL_CUTS_DIR.Data(), info.opticsID.Data(), FileID);
}

static TString RelabeledCutFile_xrelabel(const XRelabelRunInfo& info, int FileID)
{
  if (XRELABEL_CAMPAIGN_DIR.Length() > 0 &&
      XRELABEL_RUNGROUP_TAG.Length() > 0) {
    return Form("%s/03b_relabel_x/cuts/XpFpXFp_%s_auto_band_cut_xscolRelabeled.root",
                XRELABEL_CAMPAIGN_DIR.Data(),
                XRELABEL_RUNGROUP_TAG.Data());
  }

  return Form("%s/XpFpXFp_%s_%d_auto_band_cut_xscolRelabeled.root",
              XRELABEL_CUTS_DIR.Data(), info.opticsID.Data(), FileID);
}

static TString DiagnosticPdf_xrelabel(int nrun, int foilIndex, int ndelIndex)
{
  if (XRELABEL_CAMPAIGN_DIR.Length() > 0 &&
      XRELABEL_RUNGROUP_TAG.Length() > 0) {
    return Form("%s/03b_relabel_x/plots/xfp_xpfp_xscolRelabel_%s_foil%d_ndel%d.pdf",
                XRELABEL_CAMPAIGN_DIR.Data(),
                XRELABEL_RUNGROUP_TAG.Data(),
                foilIndex,
                ndelIndex);
  }

  return Form("%s/xfp_xpfp_autoBand_xscolRelabel_coloredDensity_run%d_foil%d_ndel%d.pdf",
              XRELABEL_PLOTS_DIR.Data(), nrun, foilIndex, ndelIndex);
}

static TString DiagnosticRoot_xrelabel(int nrun, int foilIndex, int ndelIndex)
{
  if (XRELABEL_CAMPAIGN_DIR.Length() > 0 &&
      XRELABEL_RUNGROUP_TAG.Length() > 0) {
    return Form("%s/03b_relabel_x/root/xfp_xpfp_xscolRelabel_%s_foil%d_ndel%d.root",
                XRELABEL_CAMPAIGN_DIR.Data(),
                XRELABEL_RUNGROUP_TAG.Data(),
                foilIndex,
                ndelIndex);
  }

  return Form("%s/xfp_xpfp_autoBand_xscolRelabel_coloredDensity_run%d_foil%d_ndel%d.root",
              XRELABEL_PLOTS_DIR.Data(), nrun, foilIndex, ndelIndex);
}

static double XsGuide_xrelabel(int xscol)
{
  return (xscol - 4) * 2.54; // cm, from set_xpfp_xfp_cuts.C / plot_xfp_cuts.C
}

static double YsGuide_xrelabel(int yscol)
{
  return (yscol - 4) * 0.6 * 2.54; // cm, shown faintly for orientation only
}

static int NearestXscol_xrelabel(double xs)
{
  int best = 0;
  double bestDist = 1.0e99;
  for (int i = 0; i < 9; ++i) {
    double d = std::fabs(xs - XsGuide_xrelabel(i));
    if (d < bestDist) {
      bestDist = d;
      best = i;
    }
  }
  return best;
}

static int DominantXscolFromCounts_xrelabel(const std::vector<Long64_t>& counts,
                                            Long64_t& bestCount,
                                            double& bestFrac)
{
  int best = -1;
  bestCount = 0;
  Long64_t total = 0;

  for (int i = 0; i < (int)counts.size(); ++i) {
    total += counts[i];
    if (counts[i] > bestCount) {
      bestCount = counts[i];
      best = i;
    }
  }

  bestFrac = (total > 0) ? ((double)bestCount / (double)total) : 0.0;
  return best;
}

static void DrawXsGuides_xrelabel(double ySieveMin=-7.0, double ySieveMax=7.0)
{
  for (int nxs = 0; nxs < 9; ++nxs) {
    double pos = XsGuide_xrelabel(nxs);
    TLine* line = new TLine(ySieveMin, pos, ySieveMax, pos);
    line->SetLineColor(kRed+1);
    line->SetLineWidth(1);
    line->Draw("same");

    TText* txt = new TText(ySieveMin - 0.45, pos, Form("%d", nxs));
    txt->SetTextColor(kRed+1);
    txt->SetTextSize(0.035);
    txt->SetTextAlign(32);
    txt->Draw("same");
  }
}

static void DrawYsOrientationGuides_xrelabel(double xSieveMin=-12.5, double xSieveMax=12.5)
{
  for (int nys = 0; nys < 9; ++nys) {
    double pos = YsGuide_xrelabel(nys);
    TLine* line = new TLine(pos, xSieveMin, pos, xSieveMax);
    line->SetLineColor(kGray+1);
    line->SetLineStyle(3);
    line->SetLineWidth(1);
    line->Draw("same");
  }
}

static void DrawFrame_xrelabel(const TString& title)
{
  TH2D* frame = new TH2D("frame_xrelabel_tmp", title + ";ysieve (cm);xsieve (cm)",
                         10, -7.0, 7.0, 10, -12.5, 12.5);
  frame->SetDirectory(nullptr);
  frame->SetStats(0);
  frame->Draw("AXIS");
}

static bool ZoneAllowed_xrelabel(const TString& zone, const TString& zoneModeIn)
{
  TString mode = zoneModeIn;
  mode.ToLower();
  TString z = zone;
  z.ToLower();

  if (mode.Length() == 0 || mode == "all" || mode == "*") return true;
  if (mode == "split" || mode == "lowhigh" || mode == "low_high")
    return (z == "low" || z == "high");
  return z == mode;
}

static int ExtractCandidateCut_xrelabel(const TString& name,
                                        int foilIndex,
                                        int ndelIndex,
                                        TString& zoneOut)
{
  int cand=-1, nf=-1, nd=-1;
  char zoneBuf[64];
  zoneBuf[0] = '\0';

  int ok = sscanf(name.Data(),
                  "hXpFpXFp_cut_candidate_%d_zone_%63[^_]_nfoil_%d_ndel_%d",
                  &cand, zoneBuf, &nf, &nd);
  if (ok == 4 && nf == foilIndex && nd == ndelIndex) {
    zoneOut = zoneBuf;
    return cand;
  }
  zoneOut = "";
  return -1;
}


static bool ExtractRelabeledComponent_xrelabel(const TString& name,
                                               int foilIndex,
                                               int ndelIndex,
                                               int& xscolOut,
                                               int& partOut,
                                               TString& zoneOut,
                                               int& srcOut)
{
  int xs=-1, nf=-1, nd=-1, part=-1, src=-1;
  char zoneBuf[64];
  zoneBuf[0] = '\0';

  int ok = sscanf(name.Data(),
                  "hXpFpXFp_cut_xscol_%d_nfoil_%d_ndel_%d_part%d_zone_%63[^_]_src_%d",
                  &xs, &nf, &nd, &part, zoneBuf, &src);
  if (ok == 6 && nf == foilIndex && nd == ndelIndex) {
    xscolOut = xs;
    partOut = part;
    zoneOut = zoneBuf;
    srcOut = src;
    return true;
  }

  xscolOut = -1;
  partOut = -1;
  srcOut = -1;
  zoneOut = "";
  return false;
}

static int ExtractOldXscol_xrelabel(const TString& name, int foilIndex, int ndelIndex)
{
  int xs=-1, nf=-1, nd=-1;
  int ok = sscanf(name.Data(), "hXpFpXFp_cut_xscol_%d_nfoil_%d_ndel_%d", &xs, &nf, &nd);
  if (ok == 3 && nf == foilIndex && nd == ndelIndex) return xs;
  return -1;
}

static TCutG* LoadYtarCut_xrelabel(int nrun, const TString& ytarTag, int foilIndex)
{
  TString fname;

  if (XRELABEL_CAMPAIGN_DIR.Length() > 0 &&
      XRELABEL_RUNGROUP_TAG.Length() > 0) {
    fname = Form("%s/01_ytar_cuts/cuts/ytar_ridge_cut_%s.root",
                 XRELABEL_CAMPAIGN_DIR.Data(),
                 XRELABEL_RUNGROUP_TAG.Data());
  } else {
    fname = Form("%s/ytar_ridge_cut_%s_run%d.root",
                 XRELABEL_CUTS_DIR.Data(),
                 ytarTag.Data(),
                 nrun);
  }
  TFile* f = TFile::Open(fname, "READ");
  if (!f || f->IsZombie()) {
    cout << "WARNING: could not open ytar cut file: " << fname << endl;
    return nullptr;
  }

  std::vector<TString> names = {
    Form("delta_vs_ytar_cut_foil%d", foilIndex),
    Form("ytar_delta_cut_foil%d", foilIndex),
    Form("foil%d", foilIndex),
    Form("cut_foil%d", foilIndex),
    "delta_vs_ytar_cut",
    "ytar_delta_cut"
  };

  for (const auto& nm : names) {
    TCutG* c = (TCutG*)f->Get(nm);
    if (c) {
      TCutG* cc = (TCutG*)c->Clone(Form("ytar_for_xrelabel_foil%d", foilIndex));
      cout << "Loaded ytar cut: " << fname << " :: " << nm << endl;
      return cc;
    }
  }

  cout << "WARNING: no usable ytar cut found in " << fname << endl;
  return nullptr;
}

struct XOverrideEntry {
  int newXscol = -999;
  TString action = "";
  TString comment = "";
};

static TString OverrideKey_xrelabel(int run, int foil, int ndel, const TString& zone, int oldIndex)
{
  return Form("%d_%d_%d_%s_%d", run, foil, ndel, zone.Data(), oldIndex);
}

static std::map<TString, XOverrideEntry> ReadOverrides_xrelabel(const TString& fname)
{
  std::map<TString, XOverrideEntry> out;
  if (fname.Length() == 0) return out;

  std::ifstream fin(fname.Data());
  if (!fin.is_open()) {
    cout << "No override file opened: " << fname << endl;
    return out;
  }

  std::string raw;
  while (std::getline(fin, raw)) {
    TString line(raw.c_str());
    line = line.Strip(TString::kBoth);
    if (line.Length() == 0 || line.BeginsWith("#")) continue;

    // Preferred format:
    //   run foil ndel zone oldIndex newXscol action comment
    std::istringstream iss(raw);
    int run=-1, foil=-1, ndel=-1, oldIndex=-1, newXs=-999;
    std::string zoneStd, actionStd, commentStd;

    bool parsed = false;
    if (iss >> run >> foil >> ndel >> zoneStd >> oldIndex >> newXs >> actionStd) {
      parsed = true;
    } else {
      // Backward-compatible format:
      //   run foil ndel oldIndex newXscol action comment
      iss.clear();
      iss.str(raw);
      if (iss >> run >> foil >> ndel >> oldIndex >> newXs >> actionStd) {
        zoneStd = "*";
        parsed = true;
      }
    }
    if (!parsed) continue;

    std::getline(iss, commentStd);

    TString zone = zoneStd.c_str();
    XOverrideEntry e;
    e.newXscol = newXs;
    e.action = actionStd.c_str();
    e.comment = commentStd.c_str();
    out[OverrideKey_xrelabel(run, foil, ndel, zone, oldIndex)] = e;
  }

  cout << "Loaded override entries: " << out.size() << " from " << fname << endl;
  return out;
}

static bool LookupOverride_xrelabel(const std::map<TString, XOverrideEntry>& overrides,
                                    int run,
                                    int foil,
                                    int ndel,
                                    const TString& zone,
                                    int oldIndex,
                                    XOverrideEntry& out)
{
  TString keyExact = OverrideKey_xrelabel(run, foil, ndel, zone, oldIndex);
  auto it = overrides.find(keyExact);
  if (it != overrides.end()) { out = it->second; return true; }

  TString keyStar = OverrideKey_xrelabel(run, foil, ndel, "*", oldIndex);
  it = overrides.find(keyStar);
  if (it != overrides.end()) { out = it->second; return true; }

  TString keyAny = OverrideKey_xrelabel(run, foil, ndel, "any", oldIndex);
  it = overrides.find(keyAny);
  if (it != overrides.end()) { out = it->second; return true; }

  return false;
}

static int ShadeColor_xrelabel(int baseColor, double frac)
{
  frac = std::max(0.0, std::min(1.0, frac));
  TColor* bc = gROOT->GetColor(baseColor);
  Float_t r=0.0, g=0.0, b=0.0;
  if (bc) bc->GetRGB(r,g,b);
  else { r = 0.2; g = 0.2; b = 0.2; }

  // Blend white -> base color. frac=0 very pale, frac=1 base color.
  Float_t rr = (Float_t)(1.0 - (1.0 - r) * frac);
  Float_t gg = (Float_t)(1.0 - (1.0 - g) * frac);
  Float_t bb = (Float_t)(1.0 - (1.0 - b) * frac);
  return TColor::GetColor(rr, gg, bb);
}

static void DrawDensityBoxes_xrelabel(TH2D* h, int baseColor, bool logScale=true)
{
  if (!h) return;
  double hmax = h->GetMaximum();
  if (hmax <= 0) return;

  for (int ix = 1; ix <= h->GetNbinsX(); ++ix) {
    for (int iy = 1; iy <= h->GetNbinsY(); ++iy) {
      double v = h->GetBinContent(ix, iy);
      if (v <= 0) continue;

      double frac = 0.0;
      if (logScale) frac = std::log(1.0 + v) / std::log(1.0 + hmax);
      else          frac = v / hmax;

      frac = 0.18 + 0.82 * frac; // keep sparse bins visible
      int col = ShadeColor_xrelabel(baseColor, frac);

      double x1 = h->GetXaxis()->GetBinLowEdge(ix);
      double x2 = h->GetXaxis()->GetBinUpEdge(ix);
      double y1 = h->GetYaxis()->GetBinLowEdge(iy);
      double y2 = h->GetYaxis()->GetBinUpEdge(iy);

      TBox* box = new TBox(x1, y1, x2, y2);
      box->SetFillColor(col);
      box->SetLineColor(col);
      box->Draw("same");
    }
  }
}

static void AddBandToCombined_xrelabel(TH2D* hDest, TH2D* hSrc)
{
  if (!hDest || !hSrc) return;
  hDest->Add(hSrc);
}

void relabel_xfp_xpfp_autoBands_to_xscol_coloredDensity_batch_oneDelta(Int_t nrun=1544,
                                                        Int_t foilIndex=0,
                                                        Int_t ndelIndex=0,
                                                        const char* ytarTag="ML_dev",
                                                        Int_t FileID=-1,
                                                        const char* overrideFile="",
                                                        Bool_t useYtarCut=true,
                                                        Long64_t maxEvents=-1,
                                                        Double_t minCerNpe=6.0,
                                                        Double_t minCalEtot=0.65,
                                                        Int_t minEvents=25,
                                                        Bool_t writeRelabeledCuts=true,
                                                        const char* zoneMode="all",
                                                        Bool_t writeSingletonAliases=false,
                                                        Double_t xpfpZoneSplit=-999.0)
{
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);

  if (XRELABEL_CAMPAIGN_DIR.Length() > 0) {
    gSystem->mkdir(Form("%s/03b_relabel_x/cuts", XRELABEL_CAMPAIGN_DIR.Data()), kTRUE);
    gSystem->mkdir(Form("%s/03b_relabel_x/plots", XRELABEL_CAMPAIGN_DIR.Data()), kTRUE);
    gSystem->mkdir(Form("%s/03b_relabel_x/root", XRELABEL_CAMPAIGN_DIR.Data()), kTRUE);
  } else {
    gSystem->mkdir(XRELABEL_CUTS_DIR.Data(), kTRUE);
    gSystem->mkdir(XRELABEL_PLOTS_DIR.Data(), kTRUE);
  }

  XRelabelRunInfo info;
  TString metaFile = (XRELABEL_CAMPAIGN_DIR.Length() > 0)
                   ? "DATfiles/list_of_optics_run.dat"
                   : XRELABEL_DAT_DIR + "/list_of_optics_run.dat";
  if (!ReadOpticsRunInfo_xrelabel(nrun, info, metaFile.Data())) return;

  if (foilIndex < 0 || foilIndex >= info.numFoil) {
    cout << "ERROR: foilIndex out of range: " << foilIndex << endl;
    return;
  }
  if (ndelIndex < 0 || ndelIndex + 1 >= (int)info.delcut.size()) {
    cout << "ERROR: ndelIndex out of range: " << ndelIndex << endl;
    return;
  }

  const double deltaMin = info.delcut[ndelIndex];
  const double deltaMax = info.delcut[ndelIndex+1];

  TString inputTree = TreeRootPath_xrelabel(nrun, FileID);
  TString inputCuts = AutoBandCutFile_xrelabel(info, FileID);
  TString outputCuts = RelabeledCutFile_xrelabel(info, FileID);

  if (XRELABEL_INPUT_CUT_OVERRIDE.Length() > 0) inputCuts = XRELABEL_INPUT_CUT_OVERRIDE;
  if (XRELABEL_OUTPUT_CUT_OVERRIDE.Length() > 0) outputCuts = XRELABEL_OUTPUT_CUT_OVERRIDE;

  TString outpdf = DiagnosticPdf_xrelabel(nrun, foilIndex, ndelIndex);
  TString outroot = DiagnosticRoot_xrelabel(nrun, foilIndex, ndelIndex);

  TString zoneModeStr = zoneMode;

  const bool useXpfpZoneGate = (xpfpZoneSplit > -998.0);
  auto PassesXpfpZoneGate_xrelabel = [&](const TString& zone, double xpfpVal) {
    if (!useXpfpZoneGate) return true;
    if (zone == "low")  return (xpfpVal <  xpfpZoneSplit);
    if (zone == "high") return (xpfpVal >= xpfpZoneSplit);
    return true;
  };

  cout << "\n=== XFP/XPFP auto-band xscol relabel + colored-density diagnostic ===" << endl;
  cout << "Input tree       : " << inputTree << endl;
  cout << "Input cuts       : " << inputCuts << endl;
  cout << "Output cuts      : " << outputCuts << endl;
  cout << "Diagnostic PDF   : " << outpdf << endl;
  cout << "Diagnostic ROOT  : " << outroot << endl;
  cout << "Run/foil/ndel    : " << nrun << " / " << foilIndex << " / " << ndelIndex << endl;
  cout << "Delta slice      : [" << deltaMin << ", " << deltaMax << ")" << endl;
  cout << "Zone mode        : " << zoneModeStr << endl;
  if (useXpfpZoneGate) {
    cout << "XPFp zone split  : " << xpfpZoneSplit
         << "  (low: xpfp < split, high: xpfp >= split)" << endl;
  } else {
    cout << "XPFp zone split  : not applied" << endl;
    cout << "WARNING: if low/high candidate cuts are being relabeled, their raw xpfp" << endl;
    cout << "         gate is NOT stored inside the TCutG. Without xpfpZoneSplit," << endl;
    cout << "         relabel counts can exceed the auto-band diagnostic counts." << endl;
  }
  cout << "Duplicate xscol  : kept as multiple component cuts and combined in plots" << endl;
  cout << "===================================================================\n" << endl;

  TFile* fCuts = TFile::Open(inputCuts, "READ");
  if (!fCuts || fCuts->IsZombie()) {
    cout << "ERROR: cannot open input cut file: " << inputCuts << endl;
    return;
  }

  std::vector<XBandCutInfo> bands;
  std::set<TString> seenCutNames;
  TIter next(fCuts->GetListOfKeys());
  TKey* key = nullptr;
  while ((key = (TKey*)next())) {
    TString name = key->GetName();
    if (seenCutNames.count(name)) continue;
    seenCutNames.insert(name);

    TString zone = "";
    int oldIndex = -1;
    int cand = -1;
    int oldXs = -1;
    TString sourceKind = "";

    if (name.BeginsWith("hXpFpXFp_cut_candidate_")) {
      cand = ExtractCandidateCut_xrelabel(name, foilIndex, ndelIndex, zone);
      if (cand < 0) continue;
      if (!ZoneAllowed_xrelabel(zone, zoneModeStr)) continue;
      oldIndex = cand;
      sourceKind = "candidate";
    } else if (name.BeginsWith("hXpFpXFp_cut_xscol_")) {
      int parsedXs=-1, parsedPart=-1, parsedSrc=-1;
      TString parsedZone="";
      if (ExtractRelabeledComponent_xrelabel(name, foilIndex, ndelIndex,
                                             parsedXs, parsedPart, parsedZone, parsedSrc)) {
        oldXs = parsedXs;
        zone = parsedZone;
        if (!ZoneAllowed_xrelabel(zone, zoneModeStr)) continue;
        oldIndex = parsedSrc;
        sourceKind = "xscol_component";
      } else {
        oldXs = ExtractOldXscol_xrelabel(name, foilIndex, ndelIndex);
        if (oldXs < 0) continue;
        zone = "global";
        if (!ZoneAllowed_xrelabel(zone, zoneModeStr)) continue;
        oldIndex = oldXs;
        sourceKind = "xscol";
      }
    } else {
      continue;
    }

    TObject* obj = key->ReadObj();
    if (!obj || !obj->InheritsFrom(TCutG::Class())) continue;

    XBandCutInfo b;
    b.oldIndex = oldIndex;
    b.oldXscol = oldXs;
    b.candidateIndex = cand;
    b.zone = zone;
    b.sourceKind = sourceKind;
    b.oldName = name;
    b.cut = (TCutG*)obj->Clone(Form("old_xfp_auto_cut_%s_%d_zone_%s_foil%d_ndel%d",
                                    sourceKind.Data(), oldIndex, zone.Data(), foilIndex, ndelIndex));
    b.hYsXs = new TH2D(Form("hYsXs_%s_%d_zone_%s", sourceKind.Data(), oldIndex, zone.Data()),
                       Form("%s old %d zone %s;ysieve (cm);xsieve (cm)",
                            sourceKind.Data(), oldIndex, zone.Data()),
                       100, -7.0, 7.0, 100, -12.5, 12.5);
    b.hYsXs->SetDirectory(nullptr);
    b.hXpX = new TH2D(Form("hXpFpXFp_%s_%d_zone_%s", sourceKind.Data(), oldIndex, zone.Data()),
                      Form("%s old %d zone %s;xpfp;xfp", sourceKind.Data(), oldIndex, zone.Data()),
                      120, -0.08, 0.08, 140, -70.0, 70.0);
    b.hXpX->SetDirectory(nullptr);
    bands.push_back(b);
  }

  std::sort(bands.begin(), bands.end(), [](const XBandCutInfo& a, const XBandCutInfo& b){
    if (a.zone != b.zone) return a.zone < b.zone;
    if (a.sourceKind != b.sourceKind) return a.sourceKind < b.sourceKind;
    return a.oldIndex < b.oldIndex;
  });

  if (bands.empty()) {
    cout << "ERROR: no XFP/XPFP auto-band cuts found for foil " << foilIndex
         << " ndel " << ndelIndex << " zoneMode=" << zoneModeStr
         << " in " << inputCuts << endl;
    return;
  }

  cout << "Found auto-band cuts for this slice:" << endl;
  for (const auto& b : bands) {
    cout << "  " << b.oldName
         << "  source=" << b.sourceKind
         << "  zone=" << b.zone
         << "  oldIndex=" << b.oldIndex
         << endl;
  }

  auto overrides = ReadOverrides_xrelabel(overrideFile);

  TFile* fin = TFile::Open(inputTree, "READ");
  if (!fin || fin->IsZombie()) {
    cout << "ERROR: cannot open input tree file: " << inputTree << endl;
    return;
  }

  TTree* T = (TTree*)fin->Get("T");
  if (!T) {
    cout << "ERROR: tree T not found in " << inputTree << endl;
    return;
  }

  Double_t sumnpe=0, etracknorm=0;
  Double_t ytar=0, delta=0;
  Double_t xfp=0, xpfp=0;
  Double_t ysieve=0, xsieve=0;

  T->SetBranchStatus("*",0);
  T->SetBranchStatus("H.cer.npeSum",1);
  T->SetBranchStatus("H.cal.etottracknorm",1);
  T->SetBranchStatus("H.gtr.y",1);
  T->SetBranchStatus("H.gtr.dp",1);
  T->SetBranchStatus("H.dc.x_fp",1);
  T->SetBranchStatus("H.dc.xp_fp",1);
  T->SetBranchStatus("H.extcor.ysieve",1);
  T->SetBranchStatus("H.extcor.xsieve",1);

  T->SetBranchAddress("H.cer.npeSum", &sumnpe);
  T->SetBranchAddress("H.cal.etottracknorm", &etracknorm);
  T->SetBranchAddress("H.gtr.y", &ytar);
  T->SetBranchAddress("H.gtr.dp", &delta);
  T->SetBranchAddress("H.dc.x_fp", &xfp);
  T->SetBranchAddress("H.dc.xp_fp", &xpfp);
  T->SetBranchAddress("H.extcor.ysieve", &ysieve);
  T->SetBranchAddress("H.extcor.xsieve", &xsieve);

  TCutG* ytarCut = nullptr;
  if (useYtarCut) ytarCut = LoadYtarCut_xrelabel(nrun, ytarTag, foilIndex);

  std::vector<double> sumXs(bands.size(), 0.0), sumYs(bands.size(), 0.0);
  std::vector<double> sumXs2(bands.size(), 0.0), sumYs2(bands.size(), 0.0);

  // Assign XFP/XPFP components by dominant physical sieve row, not by mean.
  // The mean can be pulled across an adjacent row by tails/overlap, which is
  // exactly the failure mode in low/high zone combinations.
  std::vector< std::vector<Long64_t> > xscolCountsByBand(
    bands.size(), std::vector<Long64_t>(9, 0));

  Long64_t nentries = T->GetEntries();
  if (maxEvents > 0 && maxEvents < nentries) nentries = maxEvents;

  Long64_t nPassBasic = 0, nPassYtar = 0, nPassDelta = 0;
  for (Long64_t i = 0; i < nentries; ++i) {
    T->GetEntry(i);

    if (!(sumnpe > minCerNpe && etracknorm > minCalEtot)) continue;
    nPassBasic++;

    if (useYtarCut && ytarCut) {
      if (!ytarCut->IsInside(ytar, delta)) continue; // x=ytar, y=delta
    }
    nPassYtar++;

    if (!(delta >= deltaMin && delta < deltaMax)) continue;
    nPassDelta++;

    for (size_t ib = 0; ib < bands.size(); ++ib) {
      if (!bands[ib].cut) continue;

      // Critical XFP/XPFP zone reconstruction:
      // The angle-scan macro used the raw xpfp low/high gate BEFORE making
      // each candidate TCutG, but that gate is not encoded in the TCutG itself.
      // Reapply it here so relabel N matches the auto-band diagnostic N.
      if (!PassesXpfpZoneGate_xrelabel(bands[ib].zone, xpfp)) continue;

      if (!bands[ib].cut->IsInside(xpfp, xfp)) continue; // x=xpfp, y=xfp

      bands[ib].hYsXs->Fill(ysieve, xsieve);
      bands[ib].hXpX->Fill(xpfp, xfp);
      bands[ib].n++;
      sumXs[ib] += xsieve;
      sumYs[ib] += ysieve;
      sumXs2[ib] += xsieve * xsieve;
      sumYs2[ib] += ysieve * ysieve;

      int nearestPhysicalXscol = NearestXscol_xrelabel(xsieve);
      if (nearestPhysicalXscol >= 0 && nearestPhysicalXscol < 9) {
        xscolCountsByBand[ib][nearestPhysicalXscol]++;
      }
    }
  }

  cout << "\nEvent counts:" << endl;
  cout << "  pass PID/basic : " << nPassBasic << endl;
  cout << "  pass ytar      : " << nPassYtar << endl;
  cout << "  pass delta     : " << nPassDelta << endl;

  for (size_t ib = 0; ib < bands.size(); ++ib) {
    auto& b = bands[ib];

    if (b.n > 0) {
      b.meanXs = sumXs[ib] / (double)b.n;
      b.meanYs = sumYs[ib] / (double)b.n;
      b.rmsXs = std::sqrt(std::max(0.0, sumXs2[ib] / (double)b.n - b.meanXs*b.meanXs));
      b.rmsYs = std::sqrt(std::max(0.0, sumYs2[ib] / (double)b.n - b.meanYs*b.meanYs));

      b.xscolCounts = xscolCountsByBand[ib];
      b.modalXscol = DominantXscolFromCounts_xrelabel(b.xscolCounts,
                                                      b.modalCount,
                                                      b.modalFrac);

      // Primary assignment: dominant physical xscol from event-by-event nearest
      // sieve-row votes.  MeanXs is diagnostic only.
      b.newXscol = b.modalXscol;
      b.distToGuide = (b.newXscol >= 0) ?
        std::fabs(b.meanXs - XsGuide_xrelabel(b.newXscol)) : 1.0e99;
    } else {
      b.newXscol = -1;
      b.distToGuide = 1.0e99;
      b.xscolCounts.assign(9, 0);
      b.modalXscol = -1;
      b.modalCount = 0;
      b.modalFrac = 0.0;
    }

    XOverrideEntry oe;
    if (LookupOverride_xrelabel(overrides, nrun, foilIndex, ndelIndex, b.zone, b.oldIndex, oe)) {
      b.hasOverride = true;
      b.forcedXscol = oe.newXscol;
      b.action = oe.action;
      if (b.action == "reject" || b.forcedXscol < 0) {
        b.reject = true;
        b.write = false;
      } else {
        b.newXscol = b.forcedXscol;
        b.write = true;
      }
    } else {
      b.write = (b.n >= minEvents && b.newXscol >= 0 && b.newXscol <= 8);
    }
  }

  // Preserve duplicate xscol assignments: assign part indices instead of
  // dropping all but the largest component as the YFP/YPFP macro did.
  std::map<int,int> nPartsByXscol;
  for (const auto& b : bands) {
    if (b.write) nPartsByXscol[b.newXscol]++;
  }

  std::map<int,int> nextPartByXscol;
  for (auto& b : bands) {
    if (!b.write) continue;
    b.partIndex = nextPartByXscol[b.newXscol]++;
    b.newName = Form("hXpFpXFp_cut_xscol_%d_nfoil_%d_ndel_%d_part%02d_zone_%s_src_%02d",
                     b.newXscol, foilIndex, ndelIndex, b.partIndex, b.zone.Data(), b.oldIndex);
  }

  cout << "\nAssignment summary:" << endl;
  cout << "  source/zone/old  N  meanXs  rmsXs  modeXscol  modeFrac  newXscol  dist  part  write" << endl;
  for (const auto& b : bands) {
    cout << "  " << b.sourceKind << "/" << b.zone << "/" << b.oldIndex
         << "  N=" << b.n
         << "  meanXs=" << b.meanXs
         << "  rmsXs=" << b.rmsXs
         << "  meanYs=" << b.meanYs
         << "  modeXscol=" << b.modalXscol
         << "  modeFrac=" << b.modalFrac
         << "  newXscol=" << b.newXscol
         << "  dist=" << b.distToGuide
         << "  action=" << b.action
         << "  part=" << b.partIndex
         << "  write=" << b.write
         << "  oldName=" << b.oldName
         << endl;

    if (b.write && b.modalFrac > 0.0 && b.modalFrac < 0.65) {
      cout << "    WARNING: mixed xscol population in " << b.oldName
           << "; dominant fraction=" << b.modalFrac
           << ". Review this component page before trusting it." << endl;
    }
  }

  cout << "\nDuplicate/additive xscol groups:" << endl;
  for (const auto& kv : nPartsByXscol) {
    if (kv.second > 1) {
      cout << "  xscol " << kv.first << " has " << kv.second
           << " component cuts; downstream should OR/add them." << endl;
    }
  }

  std::vector<int> baseColors = {kRed+1, kBlue+1, kGreen+2, kMagenta+1, kOrange+7,
                                 kCyan+2, kViolet+7, kSpring+5, kPink+7, kGray+2,
                                 kAzure+7, kTeal+3, kOrange+1, kViolet+3,
                                 kRed-4, kBlue-4, kGreen-3, kMagenta-4};

  std::vector<TH2D*> hCombinedXscol(9, nullptr);
  std::vector<Long64_t> nCombinedXscol(9, 0);
  for (int xs = 0; xs < 9; ++xs) {
    hCombinedXscol[xs] = new TH2D(Form("hCombinedXscol_%d", xs),
                                  Form("combined xscol %d;ysieve (cm);xsieve (cm)", xs),
                                  100, -7.0, 7.0, 100, -12.5, 12.5);
    hCombinedXscol[xs]->SetDirectory(nullptr);
  }
  for (const auto& b : bands) {
    if (!b.write || b.newXscol < 0 || b.newXscol > 8) continue;
    AddBandToCombined_xrelabel(hCombinedXscol[b.newXscol], b.hYsXs);
    nCombinedXscol[b.newXscol] += b.n;
  }

  // Diagnostic plots.
  TCanvas* c = new TCanvas("c_xrelabel_colored_density", "xscol relabel colored density", 1150, 900);
  c->SaveAs(outpdf + "[");

  auto drawOverlayForZone = [&](const TString& zoneToDraw, const TString& titleTag) {
    c->Clear();
    DrawFrame_xrelabel(Form("Run %d foil %d ndel %d: %s XFP/XPFP candidates projected to sieve",
                            nrun, foilIndex, ndelIndex, titleTag.Data()));
    for (size_t ib = 0; ib < bands.size(); ++ib) {
      if (zoneToDraw != "all" && bands[ib].zone != zoneToDraw) continue;
      DrawDensityBoxes_xrelabel(bands[ib].hYsXs, baseColors[ib % baseColors.size()], true);
    }
    DrawYsOrientationGuides_xrelabel();
    DrawXsGuides_xrelabel();

    TLatex lat; lat.SetNDC(); lat.SetTextSize(0.026);
    double ytxt = 0.92;
    lat.DrawLatex(0.12, ytxt, Form("#delta [%.1f, %.1f), colored boxes = per-cut density; red horizontal lines = physical xscol", deltaMin, deltaMax));
    ytxt -= 0.035;
    for (size_t ib = 0; ib < bands.size() && ytxt > 0.12; ++ib) {
      if (zoneToDraw != "all" && bands[ib].zone != zoneToDraw) continue;
      lat.SetTextColor(baseColors[ib % baseColors.size()]);
      lat.DrawLatex(0.12, ytxt,
                    Form("%s %s %d #rightarrow xscol %d part %d, N=%lld%s",
                         bands[ib].sourceKind.Data(), bands[ib].zone.Data(), bands[ib].oldIndex,
                         bands[ib].newXscol, bands[ib].partIndex, bands[ib].n,
                         bands[ib].hasOverride ? " override" : ""));
      ytxt -= 0.030;
    }
    lat.SetTextColor(kBlack);
    c->SaveAs(outpdf);
  };

  drawOverlayForZone("low",  "LOW-zone verification");
  drawOverlayForZone("high", "HIGH-zone verification");
  drawOverlayForZone("all",  "ALL components before combination");

  c->Clear();
  DrawFrame_xrelabel(Form("Run %d foil %d ndel %d: FINAL combined-by-xscol density", nrun, foilIndex, ndelIndex));
  for (int xs = 0; xs < 9; ++xs) {
    if (nCombinedXscol[xs] <= 0) continue;
    DrawDensityBoxes_xrelabel(hCombinedXscol[xs], baseColors[xs % baseColors.size()], true);
  }
  DrawYsOrientationGuides_xrelabel();
  DrawXsGuides_xrelabel();

  TLatex latC; latC.SetNDC(); latC.SetTextSize(0.027);
  double ytxtC = 0.92;
  latC.DrawLatex(0.12, ytxtC, "Final combination: components with the same xscol are added/overlaid together");
  ytxtC -= 0.035;
  for (int xs = 0; xs < 9 && ytxtC > 0.12; ++xs) {
    if (nCombinedXscol[xs] <= 0) continue;
    latC.SetTextColor(baseColors[xs % baseColors.size()]);
    latC.DrawLatex(0.12, ytxtC,
                   Form("xscol %d: combined N=%lld, components=%d",
                        xs, nCombinedXscol[xs], nPartsByXscol[xs]));
    ytxtC -= 0.030;
  }
  latC.SetTextColor(kBlack);
  c->SaveAs(outpdf);

  for (int xs = 0; xs < 9; ++xs) {
    if (nCombinedXscol[xs] <= 0) continue;
    c->Clear();
    c->Divide(2,2);

    c->cd(1);
    DrawFrame_xrelabel(Form("combined xscol %d in sieve space", xs));
    DrawDensityBoxes_xrelabel(hCombinedXscol[xs], baseColors[xs % baseColors.size()], true);
    DrawYsOrientationGuides_xrelabel();
    DrawXsGuides_xrelabel();

    c->cd(2);
    gPad->SetLogz();
    hCombinedXscol[xs]->Draw("COLZ");
    DrawXsGuides_xrelabel();

    c->cd(3);
    DrawFrame_xrelabel(Form("components contributing to xscol %d", xs));

    // Per-xscol contributor colors: local component 0,1,2,... get distinct
    // colors on this page.  Do not use the global band index here; otherwise
    // the visual meaning can change from page to page and colors may look
    // accidental.
    std::map<int,int> componentColorForBandIndex;
    int localComponentColorIndex = 0;
    for (size_t ib = 0; ib < bands.size(); ++ib) {
      if (!bands[ib].write || bands[ib].newXscol != xs) continue;
      int col = baseColors[localComponentColorIndex % baseColors.size()];
      componentColorForBandIndex[(int)ib] = col;
      DrawDensityBoxes_xrelabel(bands[ib].hYsXs, col, true);
      localComponentColorIndex++;
    }
    DrawYsOrientationGuides_xrelabel();
    DrawXsGuides_xrelabel();

    c->cd(4);
    TLatex t; t.SetNDC(); t.SetTextSize(0.036);
    t.DrawLatex(0.08,0.88,Form("combined xscol_%d", xs));
    t.DrawLatex(0.08,0.80,Form("combined N = %lld", nCombinedXscol[xs]));
    t.DrawLatex(0.08,0.72,Form("component count = %d", nPartsByXscol[xs]));
    t.DrawLatex(0.08,0.64,"lower-left colors = separate contributing cuts");
    t.DrawLatex(0.08,0.58,"components:");
    double yline = 0.51;
    for (size_t ib = 0; ib < bands.size() && yline > 0.08; ++ib) {
      if (!bands[ib].write || bands[ib].newXscol != xs) continue;
      auto colIt = componentColorForBandIndex.find((int)ib);
      if (colIt != componentColorForBandIndex.end()) t.SetTextColor(colIt->second);
      else t.SetTextColor(kBlack);
      t.DrawLatex(0.10,yline,
                  Form("part %02d: %s %s old %d, N=%lld",
                       bands[ib].partIndex, bands[ib].sourceKind.Data(), bands[ib].zone.Data(),
                       bands[ib].oldIndex, bands[ib].n));
      yline -= 0.055;
    }
    t.SetTextColor(kBlack);
    c->SaveAs(outpdf);
  }

  for (size_t ib = 0; ib < bands.size(); ++ib) {
    c->Clear();
    c->Divide(2,2);

    c->cd(1);
    gPad->SetLogz();
    bands[ib].hXpX->Draw("COLZ");
    if (bands[ib].cut) {
      bands[ib].cut->SetLineColor(baseColors[ib % baseColors.size()]);
      bands[ib].cut->SetLineWidth(3);
      bands[ib].cut->Draw("same");
    }

    c->cd(2);
    DrawFrame_xrelabel(Form("%s %s %d projected to sieve", bands[ib].sourceKind.Data(), bands[ib].zone.Data(), bands[ib].oldIndex));
    DrawDensityBoxes_xrelabel(bands[ib].hYsXs, baseColors[ib % baseColors.size()], true);
    DrawYsOrientationGuides_xrelabel();
    DrawXsGuides_xrelabel();

    c->cd(3);
    gPad->SetLogz();
    bands[ib].hYsXs->Draw("COLZ");
    DrawXsGuides_xrelabel();

    c->cd(4);
    TLatex t; t.SetNDC(); t.SetTextSize(0.035);
    t.DrawLatex(0.08,0.88,Form("source: %s", bands[ib].sourceKind.Data()));
    t.DrawLatex(0.08,0.80,Form("zone/old: %s / %d", bands[ib].zone.Data(), bands[ib].oldIndex));
    t.DrawLatex(0.08,0.72,Form("assigned: xscol_%d", bands[ib].newXscol));
    t.DrawLatex(0.08,0.64,Form("part index = %d", bands[ib].partIndex));
    t.DrawLatex(0.08,0.56,Form("N = %lld", bands[ib].n));
    t.DrawLatex(0.08,0.48,Form("dominant xscol = %d  frac = %.2f", bands[ib].modalXscol, bands[ib].modalFrac));
    t.DrawLatex(0.08,0.40,Form("mean xsieve = %.3f cm, rms = %.3f", bands[ib].meanXs, bands[ib].rmsXs));
    t.DrawLatex(0.08,0.32,Form("mean dist to assigned guide = %.3f cm", bands[ib].distToGuide));
    t.DrawLatex(0.08,0.24,Form("action/write = %s / %s", bands[ib].action.Data(), bands[ib].write ? "yes" : "no"));
    t.DrawLatex(0.08,0.16,bands[ib].oldName);
    if (bands[ib].write) t.DrawLatex(0.08,0.09,bands[ib].newName);

    c->SaveAs(outpdf);
  }

  c->SaveAs(outpdf + "]");

  TFile foutDiag(outroot, "RECREATE");
  for (auto& b : bands) {
    if (b.hYsXs) b.hYsXs->Write();
    if (b.hXpX) b.hXpX->Write();
    if (b.cut) b.cut->Write(Form("source_%s", b.oldName.Data()));
  }
  for (int xs = 0; xs < 9; ++xs) {
    if (hCombinedXscol[xs]) hCombinedXscol[xs]->Write();
  }
  foutDiag.Close();

  if (writeRelabeledCuts) {
    TFile* fInCuts = TFile::Open(inputCuts, "READ");
    TFile fOutCuts(outputCuts, "RECREATE");
    if (!fInCuts || fInCuts->IsZombie() || fOutCuts.IsZombie()) {
      cout << "ERROR: could not open cut files for relabeled output." << endl;
    } else {
      TString sliceTag = Form("_nfoil_%d_ndel_%d", foilIndex, ndelIndex);

      TIter copyIter(fInCuts->GetListOfKeys());
      TKey* copyKey = nullptr;
      while ((copyKey = (TKey*)copyIter())) {
        TObject* obj = copyKey->ReadObj();
        if (!obj) continue;
        TString nm = obj->GetName();

        // Replace only this foil/delta slice.  Remove both candidate inputs and
        // any previous relabeled xscol component outputs for this same slice.
        if (nm.Contains(sliceTag) && nm.BeginsWith("hXpFpXFp_cut_candidate_")) continue;
        if (nm.Contains(sliceTag) && nm.BeginsWith("hXpFpXFp_cut_xscol_")) continue;

        fOutCuts.cd();
        obj->Write(nm, TObject::kOverwrite);
      }

      cout << "\nWriting relabeled XFP/XPFP xscol component cuts:" << endl;
      for (const auto& b : bands) {
        if (!b.write || !b.cut) {
          cout << "  skip " << b.oldName << " -> xscol " << b.newXscol << endl;
          continue;
        }

        TCutG* cnew = (TCutG*)b.cut->Clone(b.newName);
        cnew->SetName(b.newName);
        cnew->SetTitle(b.newName + ";xpfp;xfp");
        cnew->SetLineColor(kBlue+1);
        cnew->SetLineWidth(3);
        fOutCuts.cd();
        cnew->Write(b.newName, TObject::kOverwrite);
        cout << "  " << b.oldName << "  ->  " << b.newName
             << "  N=" << b.n
             << "  modeXscol=" << b.modalXscol
             << "  modeFrac=" << b.modalFrac
             << endl;
        delete cnew;

        if (writeSingletonAliases && nPartsByXscol[b.newXscol] == 1) {
          TString aliasName = Form("hXpFpXFp_cut_xscol_%d_nfoil_%d_ndel_%d",
                                   b.newXscol, foilIndex, ndelIndex);
          TCutG* calias = (TCutG*)b.cut->Clone(aliasName);
          calias->SetName(aliasName);
          calias->SetTitle(aliasName + ";xpfp;xfp");
          calias->SetLineColor(kBlue+1);
          calias->SetLineWidth(3);
          fOutCuts.cd();
          calias->Write(aliasName, TObject::kOverwrite);
          cout << "    wrote singleton compatibility alias: " << aliasName << endl;
          delete calias;
        }
      }
      fOutCuts.Close();
    }
  }

  cout << "\nWrote diagnostic PDF : " << outpdf << endl;
  cout << "Wrote diagnostic ROOT: " << outroot << endl;
  if (writeRelabeledCuts) cout << "Wrote relabeled cuts : " << outputCuts << endl;
}

// ============================================================================
// BATCH WRAPPER
//
// Default behavior: process one run/foil/ndel.
// Set runAllDeltaSlices=true to process all DATfile delta slices for this foil.
// Batch mode is cumulative and preserves previously relabeled foils by starting
// from the existing final relabeled cut file if it exists.
// ============================================================================
void relabel_xfp_xpfp_autoBands_to_xscol_coloredDensity_batch(Int_t nrun=1544,
                                                        Int_t foilIndex=0,
                                                        Int_t ndelIndex=0,
                                                        const char* ytarTag="ML_dev",
                                                        Int_t FileID=-1,
                                                        const char* overrideFile="",
                                                        Bool_t useYtarCut=true,
                                                        Long64_t maxEvents=-1,
                                                        Double_t minCerNpe=6.0,
                                                        Double_t minCalEtot=0.65,
                                                        Int_t minEvents=25,
                                                        Bool_t writeRelabeledCuts=true,
                                                        Bool_t runAllDeltaSlices=false,
                                                        const char* zoneMode="all",
                                                        Bool_t writeSingletonAliases=false,
                                                        Double_t xpfpZoneSplit=-999.0)
{
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);

  if (!runAllDeltaSlices) {
    relabel_xfp_xpfp_autoBands_to_xscol_coloredDensity_batch_oneDelta(nrun,
                                                                foilIndex,
                                                                ndelIndex,
                                                                ytarTag,
                                                                FileID,
                                                                overrideFile,
                                                                useYtarCut,
                                                                maxEvents,
                                                                minCerNpe,
                                                                minCalEtot,
                                                                minEvents,
                                                                writeRelabeledCuts,
                                                                zoneMode,
                                                                writeSingletonAliases,
                                                                xpfpZoneSplit);
    return;
  }

  XRelabelRunInfo info;
  TString metaFile = XRELABEL_DAT_DIR + "/list_of_optics_run.dat";
  if (!ReadOpticsRunInfo_xrelabel(nrun, info, metaFile.Data())) return;

  if (foilIndex < 0 || foilIndex >= info.numFoil) {
    cout << "ERROR: foilIndex out of range: " << foilIndex << endl;
    return;
  }

  gSystem->mkdir(XRELABEL_CUTS_DIR.Data(), kTRUE);
  gSystem->mkdir(XRELABEL_PLOTS_DIR.Data(), kTRUE);

  TString originalInputCuts = AutoBandCutFile_xrelabel(info, FileID);
  TString finalOutputCuts   = RelabeledCutFile_xrelabel(info, FileID);
  TString tmpA = Form("%s/.tmp_xrelabel_%s_%d_foil%d_A.root",
                      XRELABEL_CUTS_DIR.Data(), info.opticsID.Data(), FileID, foilIndex);
  TString tmpB = Form("%s/.tmp_xrelabel_%s_%d_foil%d_B.root",
                      XRELABEL_CUTS_DIR.Data(), info.opticsID.Data(), FileID, foilIndex);

  cout << endl;
  cout << "============================================================" << endl;
  cout << "BATCH MODE: relabel XFP/XPFP auto-band cuts over all delta slices" << endl;
  cout << "Run: " << nrun << "  foilIndex: " << foilIndex << endl;
  cout << "Input cuts : " << originalInputCuts << endl;
  cout << "Final cuts : " << finalOutputCuts << endl;
  cout << "Zone mode  : " << zoneMode << endl;
  if (xpfpZoneSplit > -998.0)
    cout << "XPFp split : " << xpfpZoneSplit << "  (low/high raw gate reapplied)" << endl;
  else
    cout << "XPFp split : not applied; low/high cut counts may be inflated" << endl;
  cout << "Duplicate xscol components are preserved and combined in diagnostics." << endl;
  cout << "============================================================" << endl;

  TString currentInput = originalInputCuts;

  if (writeRelabeledCuts && !gSystem->AccessPathName(finalOutputCuts.Data())) {
    currentInput = finalOutputCuts;
    cout << "Continuing from existing relabeled cuts to preserve prior foils: "
         << currentInput << endl;
  }

  const int nSlices = (int)info.delcut.size() - 1;

  for (int nd = 0; nd < nSlices; ++nd) {
    TString currentOutput;
    if (nd == nSlices - 1) currentOutput = finalOutputCuts;
    else currentOutput = (nd % 2 == 0) ? tmpA : tmpB;

    cout << endl;
    cout << "------------------------------------------------------------" << endl;
    cout << "Batch relabel delta slice ndel=" << nd
         << " [" << info.delcut[nd] << ", " << info.delcut[nd+1] << ")" << endl;
    cout << "Read cuts : " << currentInput << endl;
    cout << "Write cuts: " << currentOutput << endl;
    cout << "------------------------------------------------------------" << endl;

    XRELABEL_INPUT_CUT_OVERRIDE = currentInput;
    XRELABEL_OUTPUT_CUT_OVERRIDE = currentOutput;

    relabel_xfp_xpfp_autoBands_to_xscol_coloredDensity_batch_oneDelta(nrun,
                                                                foilIndex,
                                                                nd,
                                                                ytarTag,
                                                                FileID,
                                                                overrideFile,
                                                                useYtarCut,
                                                                maxEvents,
                                                                minCerNpe,
                                                                minCalEtot,
                                                                minEvents,
                                                                writeRelabeledCuts,
                                                                zoneMode,
                                                                writeSingletonAliases,
                                                                xpfpZoneSplit);

    if (writeRelabeledCuts) currentInput = currentOutput;
  }

  XRELABEL_INPUT_CUT_OVERRIDE = "";
  XRELABEL_OUTPUT_CUT_OVERRIDE = "";

  if (writeRelabeledCuts) {
    if (tmpA != finalOutputCuts) gSystem->Unlink(tmpA.Data());
    if (tmpB != finalOutputCuts) gSystem->Unlink(tmpB.Data());
  }

  cout << endl;
  cout << "BATCH MODE DONE for run " << nrun
       << ", foilIndex " << foilIndex << endl;
  if (writeRelabeledCuts) cout << "Final relabeled cuts: " << finalOutputCuts << endl;
}

