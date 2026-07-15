// relabel_yfp_ypfp_autoBands_to_yscol_coloredDensity_batch.C
//
// Post-process current YFP/YPFP AUTO BAND cuts whose embedded yscol labels
// may only be q-band order labels. For each existing hand-compatible cut name,
// project selected events to the HMS sieve plane, assign the nearest physical
// yscol from the standard sieve spacing, and write a relabeled cut file.
//
// Also makes a diagnostic PDF where each cluster projection is drawn with a
// monochromatic density scale in its own base color.
//
// Usage from PROJECT_DIR:
//   hcana -b -l -q 'relabel_yfp_ypfp_autoBands_to_yscol_coloredDensity_batch.C(
//       metadataRun,
//       "rungroupTag",
//       "campaignDir",
//       "inputReplayRoot",
//       foilIndex,
//       ndelIndex,
//       "",
//       true,
//       -1,
//       6.0,
//       0.65,
//       25,
//       true,
//       true)'
//
// Optional override file format:
//   # run foil ndel oldYscol newYscol action comment
//   1540 1 0 0 2 force good_by_eye
//   1540 1 0 5 -1 reject bad_split

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
// Project/campaign path convention
//
// Run from PROJECT_DIR:
//
//   PROJECT_DIR/DATfiles/list_of_optics_run.dat
//   PROJECT_DIR/<campaign>/01_ytar_cuts/cuts
//   PROJECT_DIR/<campaign>/02a_angle_scan_y/cuts   input provisional cuts
//   PROJECT_DIR/<campaign>/03a_relabel_y/cuts
//   PROJECT_DIR/<campaign>/03a_relabel_y/plots
//   PROJECT_DIR/<campaign>/03a_relabel_y/root
//
// Replay input ROOT is supplied explicitly.
// ============================================================================

// Optional internal overrides used by batch mode so each delta slice can
// build on the relabeled output from the previous slice.
static TString RELABEL_INPUT_CUT_OVERRIDE = "";
static TString RELABEL_OUTPUT_CUT_OVERRIDE = "";

struct RelabelRunInfo {
  int run = -1;
  TString opticsID = "";
  double centAngle = 0.0;
  int numFoil = 0;
  int sieveFlag = 0;
  int ndelcut = 0;
  std::vector<double> zfoil;
  std::vector<double> delcut;
};

struct BandCutInfo {
  int oldYscol = -1;
  int newYscol = -1;
  int forcedYscol = -999;
  int part = -1;
  bool reject = false;
  bool hasOverride = false;
  TString action = "auto";
  TString oldName = "";
  TString newName = "";
  TCutG* cut = nullptr;
  TH2D* hYsXs = nullptr;
  TH2D* hYpY = nullptr;
  Long64_t n = 0;
  double meanYs = 0.0;
  double meanXs = 0.0;
  double rmsYs = 0.0;
  double distToGuide = 0.0;
  bool write = false;
};

static std::vector<TString> SplitCSV_relabel(const TString& line)
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

static bool ReadOpticsRunInfo_relabel(int nrun, RelabelRunInfo& info,
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

    auto tok = SplitCSV_relabel(line);
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
    auto ztok = SplitCSV_relabel(TString(zline_raw.c_str()).Strip(TString::kBoth));
    for (int i = 0; i < info.numFoil && i < (int)ztok.size(); ++i)
      info.zfoil.push_back(ztok[i].Atof());

    std::string dline_raw;
    if (!std::getline(fin, dline_raw)) return false;
    auto dtok = SplitCSV_relabel(TString(dline_raw.c_str()).Strip(TString::kBoth));
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

static TString AutoBandCutFile_relabel(const TString& campaignDir,
                                        const TString& rungroupTag)
{
  return Form("%s/02a_angle_scan_y/cuts/YpFpYFp_%s_auto_band_cut.root",
              campaignDir.Data(), rungroupTag.Data());
}

static TString RelabeledCutFile_relabel(const TString& campaignDir,
                                        const TString& rungroupTag)
{
  return Form("%s/03a_relabel_y/cuts/"
              "YpFpYFp_%s_auto_band_cut_yscolRelabeled.root",
              campaignDir.Data(), rungroupTag.Data());
}

static TString EdgeToken_relabel(double value)
{
  TString token;
  const double av = std::fabs(value);

  if (std::fabs(av - TMath::Nint(av)) < 1.0e-6) {
    token = Form("%d", (int)TMath::Nint(av));
  } else {
    token = Form("%.6g", av);
    token.ReplaceAll(".", "p");
  }

  if (value < -1.0e-6) token = "m" + token;
  return token;
}

static TString DeltaTag_relabel(const RelabelRunInfo& info, int ndelIndex)
{
  return Form("delta_%s_to_%s",
              EdgeToken_relabel(info.delcut[ndelIndex]).Data(),
              EdgeToken_relabel(info.delcut[ndelIndex + 1]).Data());
}

static TString DiagnosticPdf_relabel(const TString& campaignDir,
                                     const TString& rungroupTag,
                                     const RelabelRunInfo& info,
                                     int foilIndex,
                                     int ndelIndex)
{
  return Form("%s/03a_relabel_y/plots/"
              "yfp_ypfp_yscolRelabel_%s_foil%d_%s.pdf",
              campaignDir.Data(),
              rungroupTag.Data(),
              foilIndex,
              DeltaTag_relabel(info, ndelIndex).Data());
}

static TString DiagnosticRoot_relabel(const TString& campaignDir,
                                      const TString& rungroupTag,
                                      const RelabelRunInfo& info,
                                      int foilIndex,
                                      int ndelIndex)
{
  return Form("%s/03a_relabel_y/root/"
              "yfp_ypfp_yscolRelabel_%s_foil%d_%s.root",
              campaignDir.Data(),
              rungroupTag.Data(),
              foilIndex,
              DeltaTag_relabel(info, ndelIndex).Data());
}

static double YsGuide_relabel(int yscol)
{
  return (yscol - 4) * 0.6 * 2.54; // cm, from plot_yfp_cuts.C
}

static int NearestYscol_relabel(double ys)
{
  int best = 0;
  double bestDist = 1.0e99;
  for (int i = 0; i < 9; ++i) {
    double d = std::fabs(ys - YsGuide_relabel(i));
    if (d < bestDist) {
      bestDist = d;
      best = i;
    }
  }
  return best;
}

static void DrawYsGuides_relabel(double yMin=-12.5, double yMax=12.5)
{
  for (int nys = 0; nys < 9; ++nys) {
    double pos = YsGuide_relabel(nys);
    TLine* line = new TLine(pos, yMin, pos, yMax);
    line->SetLineColor(kRed+1);
    line->SetLineWidth(1);
    line->Draw("same");

    TText* txt = new TText(pos, yMin - 1.0, Form("%d", nys));
    txt->SetTextColor(kRed+1);
    txt->SetTextSize(0.035);
    txt->SetTextAlign(22);
    txt->Draw("same");
  }
}

static int ExtractOldYscol_relabel(const TString& name, int foilIndex, int ndelIndex)
{
  int ys=-1, nf=-1, nd=-1;
  int ok = sscanf(name.Data(), "hYpFpYFp_cut_yscol_%d_nfoil_%d_ndel_%d", &ys, &nf, &nd);
  if (ok == 3 && nf == foilIndex && nd == ndelIndex) return ys;
  return -1;
}

static TCutG* LoadYtarCut_relabel(const TString& campaignDir,
                                   const TString& rungroupTag,
                                   int foilIndex)
{
  TString fname = Form("%s/01_ytar_cuts/cuts/ytar_ridge_cut_%s.root",
                       campaignDir.Data(),
                       rungroupTag.Data());

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
      TCutG* cc = (TCutG*)c->Clone(Form("ytar_for_relabel_foil%d", foilIndex));
      cout << "Loaded ytar cut: " << fname << " :: " << nm << endl;
      return cc;
    }
  }

  cout << "WARNING: no usable ytar cut found in " << fname << endl;
  return nullptr;
}

struct OverrideEntry {
  int newYscol = -999;
  TString action = "";
  TString comment = "";
};

static std::map<TString, OverrideEntry> ReadOverrides_relabel(const TString& fname)
{
  std::map<TString, OverrideEntry> out;
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

    std::istringstream iss(raw);
    int run=-1, foil=-1, ndel=-1, oldYs=-1, newYs=-999;
    TString action="";
    std::string actionStd, commentStd;
    if (!(iss >> run >> foil >> ndel >> oldYs >> newYs >> actionStd)) continue;
    std::getline(iss, commentStd);

    action = actionStd.c_str();
    TString key = Form("%d_%d_%d_%d", run, foil, ndel, oldYs);
    OverrideEntry e;
    e.newYscol = newYs;
    e.action = action;
    e.comment = commentStd.c_str();
    out[key] = e;
  }

  cout << "Loaded override entries: " << out.size() << " from " << fname << endl;
  return out;
}

static int ShadeColor_relabel(int baseColor, double frac)
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

static void DrawDensityBoxes_relabel(TH2D* h, int baseColor, bool logScale=true)
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

      // Keep sparse bins visible but pale.
      frac = 0.18 + 0.82 * frac;
      int col = ShadeColor_relabel(baseColor, frac);

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

static void DrawFrame_relabel(const TString& title)
{
  TH2D* frame = new TH2D("frame_relabel_tmp", title + ";ysieve (cm);xsieve (cm)",
                         10, -7.0, 7.0, 10, -12.5, 12.5);
  frame->SetDirectory(nullptr);
  frame->SetStats(0);
  frame->Draw("AXIS");
}

void relabel_yfp_ypfp_autoBands_to_yscol_coloredDensity_batch_oneDelta(
                                                        Int_t nrun=1540,
                                                        const char* rungroupTag="",
                                                        const char* campaignDir="",
                                                        const char* inputRootFile="",
                                                        Int_t foilIndex=0,
                                                        Int_t ndelIndex=0,
                                                        const char* overrideFile="",
                                                        Bool_t useYtarCut=true,
                                                        Long64_t maxEvents=-1,
                                                        Double_t minCerNpe=6.0,
                                                        Double_t minCalEtot=0.65,
                                                        Int_t minEvents=25,
                                                        Bool_t writeRelabeledCuts=true)
{
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);

  TString tag = rungroupTag;
  tag = tag.Strip(TString::kBoth);

  TString campaign = campaignDir;
  campaign = campaign.Strip(TString::kBoth);
  if (campaign.EndsWith("/")) campaign.Chop();

  TString inputTree = inputRootFile;
  inputTree = inputTree.Strip(TString::kBoth);

  if (tag.Length() == 0) {
    cout << "ERROR: rungroupTag must be supplied explicitly." << endl;
    return;
  }

  if (campaign.Length() == 0) {
    cout << "ERROR: campaignDir must be supplied explicitly." << endl;
    return;
  }

  if (inputTree.Length() == 0) {
    cout << "ERROR: inputRootFile must be supplied explicitly." << endl;
    return;
  }

  TString cutsDir = campaign + "/03a_relabel_y/cuts";
  TString plotsDir = campaign + "/03a_relabel_y/plots";
  TString rootDir = campaign + "/03a_relabel_y/root";

  gSystem->mkdir(cutsDir.Data(), kTRUE);
  gSystem->mkdir(plotsDir.Data(), kTRUE);
  gSystem->mkdir(rootDir.Data(), kTRUE);

  RelabelRunInfo info;
  TString metaFile = "DATfiles/list_of_optics_run.dat";
  if (!ReadOpticsRunInfo_relabel(nrun, info, metaFile.Data())) return;

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

  TString inputCuts = AutoBandCutFile_relabel(campaign, tag);
  TString outputCuts = RelabeledCutFile_relabel(campaign, tag);

  if (RELABEL_INPUT_CUT_OVERRIDE.Length() > 0)
    inputCuts = RELABEL_INPUT_CUT_OVERRIDE;

  if (RELABEL_OUTPUT_CUT_OVERRIDE.Length() > 0)
    outputCuts = RELABEL_OUTPUT_CUT_OVERRIDE;

  TString outpdf = DiagnosticPdf_relabel(
      campaign, tag, info, foilIndex, ndelIndex);

  TString outroot = DiagnosticRoot_relabel(
      campaign, tag, info, foilIndex, ndelIndex);

  cout << "\n=== YFP/YPFP auto-band yscol relabel + colored-density diagnostic ===" << endl;
  cout << "Input tree       : " << inputTree << endl;
  cout << "Input cuts       : " << inputCuts << endl;
  cout << "Output cuts      : " << outputCuts << endl;
  cout << "Diagnostic PDF   : " << outpdf << endl;
  cout << "Diagnostic ROOT  : " << outroot << endl;
  cout << "Rungroup         : " << tag << endl;
  cout << "Metadata run key : " << nrun << endl;
  cout << "Foil/ndel        : " << foilIndex << " / " << ndelIndex << endl;
  cout << "Delta slice      : [" << deltaMin << ", " << deltaMax << ")" << endl;
  cout << "===================================================================\n" << endl;

  TFile* fCuts = TFile::Open(inputCuts, "READ");
  if (!fCuts || fCuts->IsZombie()) {
    cout << "ERROR: cannot open input cut file: " << inputCuts << endl;
    return;
  }

  std::vector<BandCutInfo> bands;
  std::set<TString> seenCutNames;
  TIter next(fCuts->GetListOfKeys());
  TKey* key = nullptr;
  while ((key = (TKey*)next())) {
    TString name = key->GetName();
    if (!name.BeginsWith("hYpFpYFp_cut_yscol_")) continue;
    if (seenCutNames.count(name)) continue;
    seenCutNames.insert(name);

    int oldYs = ExtractOldYscol_relabel(name, foilIndex, ndelIndex);
    if (oldYs < 0) continue;

    TObject* obj = key->ReadObj();
    if (!obj || !obj->InheritsFrom(TCutG::Class())) continue;

    BandCutInfo b;
    b.oldYscol = oldYs;
    b.oldName = name;
    b.cut = (TCutG*)obj->Clone(Form("old_auto_cut_ys%d_foil%d_ndel%d", oldYs, foilIndex, ndelIndex));
    b.hYsXs = new TH2D(Form("hYsXs_oldYs%d", oldYs),
                       Form("old label yscol %d;ysieve (cm);xsieve (cm)", oldYs),
                       100, -7.0, 7.0, 100, -12.5, 12.5);
    b.hYsXs->SetDirectory(nullptr);
    b.hYpY = new TH2D(Form("hYpFpYFp_oldYs%d", oldYs),
                      Form("old label yscol %d;ypfp;yfp", oldYs),
                      100, -0.035, 0.035, 140, -35.0, 35.0);
    b.hYpY->SetDirectory(nullptr);
    bands.push_back(b);
  }

  std::sort(bands.begin(), bands.end(), [](const BandCutInfo& a, const BandCutInfo& b){
    return a.oldYscol < b.oldYscol;
  });

  if (bands.empty()) {
    cout << "ERROR: no YFP/YPFP auto-band cuts found for foil " << foilIndex
         << " ndel " << ndelIndex << " in " << inputCuts << endl;
    return;
  }

  cout << "Found auto-band cuts for this slice:" << endl;
  for (const auto& b : bands) cout << "  " << b.oldName << endl;

  auto overrides = ReadOverrides_relabel(overrideFile);

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
  Double_t yfp=0, ypfp=0;
  Double_t ysieve=0, xsieve=0;

  T->SetBranchStatus("*",0);
  T->SetBranchStatus("H.cer.npeSum",1);
  T->SetBranchStatus("H.cal.etottracknorm",1);
  T->SetBranchStatus("H.gtr.y",1);
  T->SetBranchStatus("H.gtr.dp",1);
  T->SetBranchStatus("H.dc.y_fp",1);
  T->SetBranchStatus("H.dc.yp_fp",1);
  T->SetBranchStatus("H.extcor.ysieve",1);
  T->SetBranchStatus("H.extcor.xsieve",1);

  T->SetBranchAddress("H.cer.npeSum", &sumnpe);
  T->SetBranchAddress("H.cal.etottracknorm", &etracknorm);
  T->SetBranchAddress("H.gtr.y", &ytar);
  T->SetBranchAddress("H.gtr.dp", &delta);
  T->SetBranchAddress("H.dc.y_fp", &yfp);
  T->SetBranchAddress("H.dc.yp_fp", &ypfp);
  T->SetBranchAddress("H.extcor.ysieve", &ysieve);
  T->SetBranchAddress("H.extcor.xsieve", &xsieve);

  TCutG* ytarCut = nullptr;
  if (useYtarCut)
    ytarCut = LoadYtarCut_relabel(campaign, tag, foilIndex);

  std::vector<double> sumYs(bands.size(), 0.0), sumXs(bands.size(), 0.0);
  std::vector<double> sumYs2(bands.size(), 0.0), sumXs2(bands.size(), 0.0);

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
      if (!bands[ib].cut->IsInside(ypfp, yfp)) continue;

      bands[ib].hYsXs->Fill(ysieve, xsieve);
      bands[ib].hYpY->Fill(ypfp, yfp);
      bands[ib].n++;
      sumYs[ib] += ysieve;
      sumXs[ib] += xsieve;
      sumYs2[ib] += ysieve * ysieve;
      sumXs2[ib] += xsieve * xsieve;
    }
  }

  cout << "\nEvent counts:" << endl;
  cout << "  pass PID/basic : " << nPassBasic << endl;
  cout << "  pass ytar      : " << nPassYtar << endl;
  cout << "  pass delta     : " << nPassDelta << endl;

  std::map<int,int> bestBandForYs;
  for (size_t ib = 0; ib < bands.size(); ++ib) {
    auto& b = bands[ib];

    if (b.n > 0) {
      b.meanYs = sumYs[ib] / (double)b.n;
      b.meanXs = sumXs[ib] / (double)b.n;
      b.rmsYs = std::sqrt(std::max(0.0, sumYs2[ib] / (double)b.n - b.meanYs*b.meanYs));
      b.newYscol = NearestYscol_relabel(b.meanYs);
      b.distToGuide = std::fabs(b.meanYs - YsGuide_relabel(b.newYscol));
    } else {
      b.newYscol = -1;
      b.distToGuide = 1.0e99;
    }

    TString okey = Form("%d_%d_%d_%d", nrun, foilIndex, ndelIndex, b.oldYscol);
    if (overrides.count(okey)) {
      b.hasOverride = true;
      b.forcedYscol = overrides[okey].newYscol;
      b.action = overrides[okey].action;
      if (b.action == "reject" || b.forcedYscol < 0) {
        b.reject = true;
        b.write = false;
      } else {
        b.newYscol = b.forcedYscol;
        b.write = true;
      }
    } else {
      b.write = (b.n >= minEvents && b.newYscol >= 0 && b.newYscol <= 8);
    }

  }

  // Preserve every qualifying provisional band as a separate component.
  // Multiple components may legitimately map to the same physical yscol.
  std::map<int, int> nextPartForYscol;
  std::map<int, int> nComponentsForYscol;

  for (auto& b : bands) {
    if (!b.write) continue;
    b.part = nextPartForYscol[b.newYscol]++;
    nComponentsForYscol[b.newYscol]++;
  }

  for (auto& b : bands) {
    if (!b.write) {
      b.newName = "";
      continue;
    }

    b.newName = Form(
        "hYpFpYFp_cut_yscol_%d_nfoil_%d_ndel_%d_part%02d_src_%02d",
        b.newYscol,
        foilIndex,
        ndelIndex,
        b.part,
        b.oldYscol);
  }

  cout << "\nAssignment summary:" << endl;
  cout << "  oldName  N  meanYs  rmsYs  auto/forced newYscol  dist  write" << endl;
  for (const auto& b : bands) {
    cout << "  " << b.oldName
         << "  N=" << b.n
         << "  meanYs=" << b.meanYs
         << "  rmsYs=" << b.rmsYs
         << "  newYscol=" << b.newYscol
         << "  dist=" << b.distToGuide
         << "  action=" << b.action
         << "  part=" << b.part
         << "  write=" << b.write
         << endl;
  }

  // Diagnostic plots.
  std::vector<int> baseColors = {kRed+1, kBlue+1, kGreen+2, kMagenta+1, kOrange+7,
                                 kCyan+2, kViolet+7, kSpring+5, kPink+7};

  TCanvas* c = new TCanvas("c_relabel_colored_density", "yscol relabel colored density", 1100, 900);
  c->SaveAs(outpdf + "[");

  c->Clear();
  DrawFrame_relabel(
      Form("%s foil %d %s: colored-density ysieve/xsieve overlay",
           tag.Data(),
           foilIndex,
           DeltaTag_relabel(info, ndelIndex).Data()));
  for (size_t ib = 0; ib < bands.size(); ++ib) {
    DrawDensityBoxes_relabel(bands[ib].hYsXs, baseColors[ib % baseColors.size()], true);
  }
  DrawYsGuides_relabel();

  TLatex lat; lat.SetNDC(); lat.SetTextSize(0.027);
  double ytxt = 0.92;
  lat.DrawLatex(0.12, ytxt, Form("#delta [%.1f, %.1f), colored boxes = per-cut density; red lines = physical yscol", deltaMin, deltaMax));
  ytxt -= 0.035;
  for (size_t ib = 0; ib < bands.size() && ib < 10; ++ib) {
    lat.SetTextColor(baseColors[ib % baseColors.size()]);
    lat.DrawLatex(0.12, ytxt,
                  Form("old %d #rightarrow yscol %d, N=%lld%s",
                       bands[ib].oldYscol, bands[ib].newYscol, bands[ib].n,
                       bands[ib].hasOverride ? " override" : ""));
    ytxt -= 0.030;
  }
  lat.SetTextColor(kBlack);
  c->SaveAs(outpdf);

  for (size_t ib = 0; ib < bands.size(); ++ib) {
    c->Clear();
    c->Divide(2,2);

    c->cd(1);
    gPad->SetLogz();
    bands[ib].hYpY->Draw("COLZ");
    if (bands[ib].cut) {
      bands[ib].cut->SetLineColor(baseColors[ib % baseColors.size()]);
      bands[ib].cut->SetLineWidth(3);
      bands[ib].cut->Draw("same");
    }

    c->cd(2);
    DrawFrame_relabel(Form("old yscol %d projected to sieve", bands[ib].oldYscol));
    DrawDensityBoxes_relabel(bands[ib].hYsXs, baseColors[ib % baseColors.size()], true);
    DrawYsGuides_relabel();

    c->cd(3);
    bands[ib].hYsXs->Draw("COLZ");
    DrawYsGuides_relabel();

    c->cd(4);
    TLatex t; t.SetNDC(); t.SetTextSize(0.040);
    t.DrawLatex(0.08,0.88,Form("old cut: yscol_%d", bands[ib].oldYscol));
    t.DrawLatex(0.08,0.80,Form("assigned: yscol_%d", bands[ib].newYscol));
    t.DrawLatex(0.08,0.72,Form("N = %lld", bands[ib].n));
    t.DrawLatex(0.08,0.64,Form("mean ysieve = %.3f cm", bands[ib].meanYs));
    t.DrawLatex(0.08,0.56,Form("rms ysieve = %.3f cm", bands[ib].rmsYs));
    t.DrawLatex(0.08,0.48,Form("dist to guide = %.3f cm", bands[ib].distToGuide));
    t.DrawLatex(0.08,0.40,Form("action = %s", bands[ib].action.Data()));
    t.DrawLatex(0.08,0.32,Form("part = %d", bands[ib].part));
    t.DrawLatex(0.08,0.25,Form("write = %s", bands[ib].write ? "yes" : "no"));
    t.DrawLatex(0.08,0.18,bands[ib].oldName);
    if (bands[ib].write) t.DrawLatex(0.08,0.11,bands[ib].newName);

    c->SaveAs(outpdf);
  }

  c->SaveAs(outpdf + "]");

  TFile foutDiag(outroot, "RECREATE");
  for (auto& b : bands) {
    if (b.hYsXs) b.hYsXs->Write();
    if (b.hYpY) b.hYpY->Write();
    if (b.cut) b.cut->Write(Form("source_%s", b.oldName.Data()));
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

        // Replace only this foil/delta slice.
        if (nm.BeginsWith("hYpFpYFp_cut_yscol_") && nm.Contains(sliceTag)) continue;

        fOutCuts.cd();
        obj->Write(nm, TObject::kOverwrite);
      }

      cout << "\nWriting relabeled cuts:" << endl;

      std::map<int, std::vector<const BandCutInfo*> > writtenByYscol;

      for (const auto& b : bands) {
        if (!b.write || !b.cut) {
          cout << "  skip " << b.oldName << " -> yscol " << b.newYscol << endl;
          continue;
        }

        TCutG* cnew = (TCutG*)b.cut->Clone(b.newName);
        cnew->SetName(b.newName);
        cnew->SetTitle(b.newName + ";ypfp;yfp");
        cnew->SetLineColor(kBlue+1);
        cnew->SetLineWidth(3);

        fOutCuts.cd();
        cnew->Write(b.newName, TObject::kOverwrite);

        writtenByYscol[b.newYscol].push_back(&b);

        cout << "  " << b.oldName << "  ->  " << b.newName
             << "  N=" << b.n << endl;

        delete cnew;
      }

      // Compatibility alias only when exactly one component exists.
      // For multiple components, downstream code must OR the component cuts.
      for (const auto& item : writtenByYscol) {
        int yscol = item.first;
        const auto& comps = item.second;

        if (comps.size() != 1) {
          cout << "  yscol " << yscol
               << " has " << comps.size()
               << " components; no singleton alias written." << endl;
          continue;
        }

        const BandCutInfo* b = comps.front();
        TString alias = Form(
            "hYpFpYFp_cut_yscol_%d_nfoil_%d_ndel_%d",
            yscol,
            foilIndex,
            ndelIndex);

        TCutG* calias = (TCutG*)b->cut->Clone(alias);
        calias->SetName(alias);
        calias->SetTitle(alias + ";ypfp;yfp");
        calias->SetLineColor(kBlue+1);
        calias->SetLineWidth(3);

        fOutCuts.cd();
        calias->Write(alias, TObject::kOverwrite);

        cout << "  singleton alias: " << alias << endl;
        delete calias;
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
// Default behavior is unchanged: one run/foil/ndel is processed.
// Set runAllDeltaSlices=true to process all DATfile delta slices for this foil.
// Batch mode preserves cumulative relabeling by feeding each slice's output cut
// file into the next slice and writing the final combined result at the end.
// ============================================================================
void relabel_yfp_ypfp_autoBands_to_yscol_coloredDensity_batch(
                                                        Int_t nrun=1540,
                                                        const char* rungroupTag="",
                                                        const char* campaignDir="",
                                                        const char* inputRootFile="",
                                                        Int_t foilIndex=0,
                                                        Int_t ndelIndex=0,
                                                        const char* overrideFile="",
                                                        Bool_t useYtarCut=true,
                                                        Long64_t maxEvents=-1,
                                                        Double_t minCerNpe=6.0,
                                                        Double_t minCalEtot=0.65,
                                                        Int_t minEvents=25,
                                                        Bool_t writeRelabeledCuts=true,
                                                        Bool_t runAllDeltaSlices=false)
{
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);

  if (!runAllDeltaSlices) {
    relabel_yfp_ypfp_autoBands_to_yscol_coloredDensity_batch_oneDelta(
                                                                nrun,
                                                                rungroupTag,
                                                                campaignDir,
                                                                inputRootFile,
                                                                foilIndex,
                                                                ndelIndex,
                                                                overrideFile,
                                                                useYtarCut,
                                                                maxEvents,
                                                                minCerNpe,
                                                                minCalEtot,
                                                                minEvents,
                                                                writeRelabeledCuts);
    return;
  }

  TString tag = rungroupTag;
  tag = tag.Strip(TString::kBoth);

  TString campaign = campaignDir;
  campaign = campaign.Strip(TString::kBoth);
  if (campaign.EndsWith("/")) campaign.Chop();

  TString inputTree = inputRootFile;
  inputTree = inputTree.Strip(TString::kBoth);

  if (tag.Length() == 0) {
    cout << "ERROR: rungroupTag must be supplied explicitly." << endl;
    return;
  }

  if (campaign.Length() == 0) {
    cout << "ERROR: campaignDir must be supplied explicitly." << endl;
    return;
  }

  if (inputTree.Length() == 0) {
    cout << "ERROR: inputRootFile must be supplied explicitly." << endl;
    return;
  }

  RelabelRunInfo info;
  TString metaFile = "DATfiles/list_of_optics_run.dat";
  if (!ReadOpticsRunInfo_relabel(nrun, info, metaFile.Data())) return;

  if (foilIndex < 0 || foilIndex >= info.numFoil) {
    cout << "ERROR: foilIndex out of range: " << foilIndex << endl;
    return;
  }

  TString cutsDir = campaign + "/03a_relabel_y/cuts";
  TString plotsDir = campaign + "/03a_relabel_y/plots";
  TString rootDir = campaign + "/03a_relabel_y/root";

  gSystem->mkdir(cutsDir.Data(), kTRUE);
  gSystem->mkdir(plotsDir.Data(), kTRUE);
  gSystem->mkdir(rootDir.Data(), kTRUE);

  TString originalInputCuts = AutoBandCutFile_relabel(campaign, tag);
  TString finalOutputCuts   = RelabeledCutFile_relabel(campaign, tag);

  TString tmpA = Form("%s/.tmp_relabel_%s_foil%d_A.root",
                      cutsDir.Data(), tag.Data(), foilIndex);

  TString tmpB = Form("%s/.tmp_relabel_%s_foil%d_B.root",
                      cutsDir.Data(), tag.Data(), foilIndex);

  cout << endl;
  cout << "============================================================" << endl;
  cout << "BATCH MODE: relabel YFP/YPFP auto-band cuts over all delta slices" << endl;
  cout << "Rungroup: " << tag
       << "  metadata run: " << nrun
       << "  foilIndex: " << foilIndex << endl;
  cout << "Input cuts : " << originalInputCuts << endl;
  cout << "Final cuts : " << finalOutputCuts << endl;
  cout << "Plots are suppressed with gROOT->SetBatch(kTRUE)." << endl;
  cout << "============================================================" << endl;

  TString currentInput = originalInputCuts;

  // IMPORTANT:
  // Batch mode is cumulative over delta slices for one foil. When processing
  // multiple foils one-at-a-time, preserve any previously relabeled foil cuts
  // by starting from the existing final relabeled file if it exists.
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

    RELABEL_INPUT_CUT_OVERRIDE = currentInput;
    RELABEL_OUTPUT_CUT_OVERRIDE = currentOutput;

    relabel_yfp_ypfp_autoBands_to_yscol_coloredDensity_batch_oneDelta(
                                                                nrun,
                                                                tag.Data(),
                                                                campaign.Data(),
                                                                inputTree.Data(),
                                                                foilIndex,
                                                                nd,
                                                                overrideFile,
                                                                useYtarCut,
                                                                maxEvents,
                                                                minCerNpe,
                                                                minCalEtot,
                                                                minEvents,
                                                                writeRelabeledCuts);

    if (writeRelabeledCuts) currentInput = currentOutput;
  }

  RELABEL_INPUT_CUT_OVERRIDE = "";
  RELABEL_OUTPUT_CUT_OVERRIDE = "";

  // Clean up intermediate cumulative files. Leave the final output in place.
  if (writeRelabeledCuts) {
    if (tmpA != finalOutputCuts) gSystem->Unlink(tmpA.Data());
    if (tmpB != finalOutputCuts) gSystem->Unlink(tmpB.Data());
  }

  cout << endl;
  cout << "BATCH MODE DONE for run " << nrun
       << ", foilIndex " << foilIndex << endl;
  if (writeRelabeledCuts) cout << "Final relabeled cuts: " << finalOutputCuts << endl;
}

