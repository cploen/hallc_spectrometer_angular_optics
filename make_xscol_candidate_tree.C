// make_xscol_candidate_tree.C
// v3: preserves low/high xpfp zone gates for relabeled component cuts.
//
// Bridge macro for XFP/XPFP -> (xsieve,ysieve) GMM cleanup.
//
// Purpose:
//   Apply the established foil/delta ytar cuts and final relabeled XFP/XPFP
//   xscol cuts to a replay ROOT tree, then write a flat candidate tree with
//   one row per accepted candidate event/xscol assignment.
//
// This mirrors make_yscol_candidate_tree.C, but handles the X-side relabeler
// convention where one physical xscol may be represented by multiple component
// TCutGs:
//   hXpFpXFp_cut_xscol_4_nfoil_0_ndel_3_part00_zone_low_src_02
//   hXpFpXFp_cut_xscol_4_nfoil_0_ndel_3_part01_zone_high_src_01
//
// Those components are treated as an OR/additive collection for the same
// xscol/foil/delta, but each component also keeps its low/high xpfp-zone gate.
// This matters because the raw xpfp split is not stored inside the TCutG;
// applying only the polygon can double-count the complementary zone.
// Singleton compatibility aliases are used only when no component cuts were
// found for that xscol/foil/delta, so the same cut is not double-counted.
//
// Usage from PROJECT_DIR:
//   hcana -b -l -q 'make_xscol_candidate_tree.C(
//       metadataRun,
//       "rungroupTag",
//       "campaignDir",
//       "inputReplayRoot")'
//
// Example:
//   hcana -b -l -q 'make_xscol_candidate_tree.C(
//       611704,
//       "rg04_theta12p395_foilpm3",
//       "HMS_6p117GeV",
//       "/path/to/nps_hms_optics_6p117_rg04_theta12p395_foilpm3.root")'
//
// Remaining arguments optionally override cut/output paths and baseline cuts.
//
// Output:
//   TTree "TXCand" with event-level candidate rows.
//   TH2D overview_before_foil*_ndel* histograms of xsieve vs ysieve by slice.
//
// Notes:
//   - User-facing xscol is only the final assigned xscol label 0..8.
//   - No old/provisional xscol/qband branch is kept.
//   - If X cuts overlap, the same event can appear more than once with
//     different xscol. nXMatches is stored for diagnosing that case.
//   - nXComponentMatches counts how many component polygons matched for the
//     accepted xscol. For normal non-overlapping low/high components this is 1.

#include <TSystem.h>
#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TH2D.h>
#include <TCutG.h>
#include <TROOT.h>
#include <TMath.h>
#include <TStyle.h>
#include <TKey.h>
#include <TObjArray.h>
#include <TObjString.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <set>
#include <algorithm>

using namespace std;

namespace MakeXCand {

static vector<TString> SplitCSV(const TString& line) {
  vector<TString> out;
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

struct RunInfo {
  Int_t run = -1;
  TString opticsID = "";
  Double_t centAngle = 0.0;
  Int_t numFoil = 0;
  Int_t sieveFlag = 0;
  Int_t ndelcut = 0;
  Double_t ymis = 0.0;
  vector<Double_t> zfoil;
  vector<Double_t> delcut;
};

static bool ReadOpticsRunInfo(Int_t nrun, RunInfo& info, const TString& metaFile) {
  ifstream fin(metaFile.Data());
  if (!fin.is_open()) {
    cout << "ERROR: cannot open metadata file: " << metaFile << endl;
    return false;
  }

  string raw;
  while (getline(fin, raw)) {
    TString line(raw.c_str());
    line = line.Strip(TString::kBoth);
    if (line.Length() == 0 || line.BeginsWith("#")) continue;

    vector<TString> tok = SplitCSV(line);
    if (tok.size() < 6) continue;
    if (tok[0].Atoi() != nrun) continue;

    info.run       = tok[0].Atoi();
    info.opticsID  = tok[1];
    info.centAngle = tok[2].Atof();
    info.numFoil   = tok[3].Atoi();
    info.sieveFlag = tok[4].Atoi();
    info.ndelcut   = tok[5].Atoi();
    info.ymis      = (tok.size() > 6) ? tok[6].Atof() : 0.0;

    string zline_raw;
    if (!getline(fin, zline_raw)) return false;
    vector<TString> ztok = SplitCSV(TString(zline_raw.c_str()).Strip(TString::kBoth));
    for (Int_t i = 0; i < info.numFoil && i < (Int_t)ztok.size(); ++i) {
      info.zfoil.push_back(ztok[i].Atof());
    }

    string dline_raw;
    if (!getline(fin, dline_raw)) return false;
    vector<TString> dtok = SplitCSV(TString(dline_raw.c_str()).Strip(TString::kBoth));
    for (Int_t i = 0; i < (Int_t)dtok.size(); ++i) {
      info.delcut.push_back(dtok[i].Atof());
    }

    cout << "Parsed run " << info.run
         << " OpticsID=" << info.opticsID
         << " CentAngleDeg=" << info.centAngle
         << " NumFoil=" << info.numFoil
         << " ndelcut=" << info.ndelcut
         << " delta edges=" << info.delcut.size()
         << endl;
    return true;
  }

  cout << "ERROR: run " << nrun << " not found in " << metaFile << endl;
  return false;
}

static bool ParseXComponentCutName(const TString& name,
                                   Int_t& xscol,
                                   Int_t& foil,
                                   Int_t& ndel,
                                   Int_t& part,
                                   TString& zone,
                                   Int_t& src) {
  Int_t xs = -1, nf = -1, nd = -1, p = -1, s = -1;
  char zoneBuf[64];
  zoneBuf[0] = '\0';

  Int_t ok = sscanf(name.Data(),
                    "hXpFpXFp_cut_xscol_%d_nfoil_%d_ndel_%d_part%d_zone_%63[^_]_src_%d",
                    &xs, &nf, &nd, &p, zoneBuf, &s);
  if (ok == 6) {
    xscol = xs;
    foil = nf;
    ndel = nd;
    part = p;
    zone = zoneBuf;
    src = s;
    return true;
  }
  return false;
}

static bool ParseXSingletonCutName(const TString& name,
                                   Int_t& xscol,
                                   Int_t& foil,
                                   Int_t& ndel) {
  Int_t xs = -1, nf = -1, nd = -1;
  Int_t ok = sscanf(name.Data(), "hXpFpXFp_cut_xscol_%d_nfoil_%d_ndel_%d", &xs, &nf, &nd);
  if (ok == 3) {
    // Exclude component names: sscanf will still match the prefix.
    if (name.Contains("_part")) return false;
    xscol = xs;
    foil = nf;
    ndel = nd;
    return true;
  }
  return false;
}


static TString EdgeToken(Double_t v) {
  TString s;
  Double_t av = TMath::Abs(v);
  if (TMath::Abs(av - TMath::Nint(av)) < 1.0e-6) {
    s = Form("%d", (Int_t)TMath::Nint(av));
  } else {
    s = Form("%.6g", av);
    s.ReplaceAll(".", "p");
  }
  if (v < -1.0e-6) s = "m" + s;
  return s;
}

static TString DeltaTagFromEdges(Double_t dlo, Double_t dhi) {
  return EdgeToken(dlo) + "_to_" + EdgeToken(dhi);
}

static bool XpfpTokenToDouble(TString tok, Double_t& val) {
  tok = tok.Strip(TString::kBoth);
  if (tok.Length() == 0) return false;
  if (tok.BeginsWith("m")) tok.Replace(0, 1, "-");
  tok.ReplaceAll("p", ".");
  val = tok.Atof();
  return true;
}

static bool FindXpfpSplitFromPlotName(const TString& campaignDir,
                                      const TString& rungroupTag,
                                      Int_t foil,
                                      Int_t ndel,
                                      const RunInfo& info,
                                      Double_t& split) {
  if (ndel < 0 || ndel + 1 >= (Int_t)info.delcut.size()) return false;

  TString dtag = DeltaTagFromEdges(info.delcut[ndel], info.delcut[ndel+1]);

  TString pattern = Form(
      "%s/02b_angle_scan_x/plots/"
      "xfp_xpfp_angleScan_%s_foil%d_delta_%s_xpfp_*_to_999.pdf",
      campaignDir.Data(),
      rungroupTag.Data(),
      foil,
      dtag.Data());

  TString cmd = Form("ls -t %s 2>/dev/null | head -1", pattern.Data());
  TString path = gSystem->GetFromPipe(cmd);
  path = path.Strip(TString::kBoth);

  if (path.Length() == 0) return false;

  TString base = gSystem->BaseName(path);
  TString marker = Form("_delta_%s_xpfp_", dtag.Data());

  Ssiz_t p0 = base.Index(marker);
  if (p0 == kNPOS) return false;

  p0 += marker.Length();

  Ssiz_t p1 = base.Index("_to_999.pdf", p0);
  if (p1 == kNPOS || p1 <= p0) return false;

  TString tok = base(p0, p1 - p0);
  return XpfpTokenToDouble(tok, split);
}

static bool BindBranch(TTree* T, const char* name, void* addr, bool required = true) {
  if (!T->GetBranch(name)) {
    if (required) cout << "ERROR: required branch missing: " << name << endl;
    else          cout << "WARNING: optional branch missing: " << name << endl;
    return false;
  }
  T->SetBranchAddress(name, addr);
  return true;
}

struct XCutComponent {
  TCutG* cut = nullptr;
  Int_t part = -1;
  TString zone = "";
  Int_t src = -1;
  TString name = "";
  Bool_t useZoneGate = kFALSE;
  Double_t xpfpSplit = -999.0;
};

} // namespace MakeXCand

void make_xscol_candidate_tree(Int_t nrun = 1544,
                               TString rungroupTag = "",
                               TString campaignDir = "",
                               TString inputRootFile = "",
                               TString ytarCutOverride = "",
                               TString xCutOverride = "",
                               TString outputRootOverride = "",
                               TString outputTsvOverride = "",
                               Double_t minCerNpe = 6.0,
                               Double_t minCalEtot = 0.65,
                               Double_t deltaMinGlobal = -10.0,
                               Double_t deltaMaxGlobal = 10.0,
                               Bool_t requireYtarCut = kTRUE) {
  using namespace MakeXCand;

  gROOT->SetBatch(kTRUE);
  gStyle->SetPalette(1,0);
  gStyle->SetOptStat(1000011);

  // ------------------------------------------------------------------
  // Read run metadata from DATfiles/list_of_optics_run.dat.
  // ------------------------------------------------------------------
  TString OpticsFile = "DATfiles/list_of_optics_run.dat";
  RunInfo info;
  if (!ReadOpticsRunInfo(nrun, info, OpticsFile)) return;

  if (info.numFoil <= 0 || info.ndelcut <= 0 || (Int_t)info.delcut.size() < 2) {
    cout << "ERROR: invalid NumFoil/ndelcut/delta edge metadata." << endl;
    return;
  }

  // Important: in this workflow DATfiles/list_of_optics_run.dat stores the
  // number of delta boundary values in ndelcut for some runs, not the number
  // of intervals.  The X-side relabel macro uses delcut.size()-1 as the
  // actual number of delta slices, so mirror that here.
  const Int_t nSlices = (Int_t)info.delcut.size() - 1;
  if (nSlices <= 0) {
    cout << "ERROR: fewer than two delta boundaries were parsed." << endl;
    return;
  }
  if (info.ndelcut != nSlices) {
    cout << "NOTE: metadata ndelcut=" << info.ndelcut
         << " but parsed delta intervals=" << nSlices
         << "; using parsed intervals to match relabel workflow." << endl;
  }
  if ((Int_t)info.zfoil.size() < info.numFoil) {
    cout << "ERROR: expected " << info.numFoil << " foil z positions, found " << info.zfoil.size() << endl;
    return;
  }

  for (Int_t nf = 0; nf < info.numFoil; nf++) {
    cout << "  foil " << nf << " ztar=" << info.zfoil[nf] << endl;
  }
  for (Int_t nd = 0; nd < nSlices; nd++) {
    cout << "  ndel " << nd << " delta=[" << info.delcut[nd] << "," << info.delcut[nd+1] << ")" << endl;
  }

  // ------------------------------------------------------------------
  // Campaign-generic file names.
  //
  // Expected invocation location:
  //   PROJECT_DIR = current working directory
  //
  // Layout:
  //   PROJECT_DIR/DATfiles/list_of_optics_run.dat
  //   PROJECT_DIR/<campaign>/01_ytar_cuts
  //   PROJECT_DIR/<campaign>/02b_angle_scan_x
  //   PROJECT_DIR/<campaign>/03b_relabel_x
  //   PROJECT_DIR/<campaign>/04b_candidate_trees_x
  // ------------------------------------------------------------------
  TString tag = rungroupTag;
  if (tag.Length() == 0) tag = info.opticsID;

  TString campaign = campaignDir;
  campaign = campaign.Strip(TString::kBoth);

  if (campaign.Length() == 0) {
    cout << "ERROR: campaignDir must be supplied explicitly." << endl;
    return;
  }

  if (campaign.EndsWith("/")) campaign.Chop();

  if (inputRootFile.Length() == 0) {
    cout << "ERROR: inputRootFile must be supplied explicitly." << endl;
    return;
  }

  TString inputroot = inputRootFile;

  TString ytarCutFile = (ytarCutOverride.Length() > 0)
                      ? ytarCutOverride
                      : Form("%s/01_ytar_cuts/cuts/ytar_ridge_cut_%s.root",
                             campaign.Data(),
                             tag.Data());

  TString xCutFile = (xCutOverride.Length() > 0)
                   ? xCutOverride
                   : Form("%s/03b_relabel_x/cuts/"
                          "XpFpXFp_%s_auto_band_cut_xscolRelabeled.root",
                          campaign.Data(),
                          tag.Data());

  TString outputroot = (outputRootOverride.Length() > 0)
                     ? outputRootOverride
                     : Form("%s/04b_candidate_trees_x/root/"
                            "XscolCandidates_%s.root",
                            campaign.Data(),
                            tag.Data());

  TString outputtsv = (outputTsvOverride.Length() > 0)
                    ? outputTsvOverride
                    : Form("%s/04b_candidate_trees_x/tsv/"
                           "XscolCandidates_%s_summary.tsv",
                           campaign.Data(),
                           tag.Data());

  cout << "Rungroup tag      = " << tag << endl;
  cout << "Campaign dir      = " << campaign << endl;
  cout << "Metadata run key  = " << nrun << endl;
  cout << "Input replay ROOT = " << inputroot << endl;
  cout << "Ytar/delta cuts   = " << ytarCutFile << endl;
  cout << "XFP/XPFP cuts     = " << xCutFile << endl;
  cout << "Output ROOT       = " << outputroot << endl;
  cout << "Output TSV        = " << outputtsv << endl;
  cout << "Baseline cuts     = H.cer.npeSum > " << minCerNpe
       << ", H.cal.etottracknorm > " << minCalEtot
       << ", delta in (" << deltaMinGlobal
       << "," << deltaMaxGlobal << ")" << endl;

  // Optional exact-cut vetoes applied before candidate-tree construction.
  set<string> preGmmVetoCuts;
  set<string> seenPreGmmVetoCuts;
  TString vetoFile = Form("%s/config/pre_gmm_veto.tsv", campaign.Data());

  if (!gSystem->AccessPathName(vetoFile)) {
    ifstream vetoIn(vetoFile.Data());
    string line;

    while (getline(vetoIn, line)) {
      size_t first = line.find_first_not_of(" \t\r");
      if (first == string::npos || line[first] == '#') continue;

      string vetoRungroup, vetoAxis, vetoCutName;
      istringstream iss(line);

      if (!(iss >> vetoRungroup >> vetoAxis >> vetoCutName)) {
        cout << "WARNING: malformed pre-GMM veto line: " << line << endl;
        continue;
      }

      if (vetoRungroup == tag.Data() && vetoAxis == "x") {
        preGmmVetoCuts.insert(vetoCutName);
      }
    }
  }

  cout << "Loaded " << preGmmVetoCuts.size()
       << " X pre-GMM veto cut(s) for " << tag << endl;

  // ------------------------------------------------------------------
  // Recover the low/high raw xpfp split used by the angle-scan stage.
  // The relabeled component TCutGs do NOT contain this 1D zone gate, so
  // downstream event selection must apply it using the zone encoded in the
  // component name.  This makes tree counts match the relabel colored-density
  // diagnostics instead of accepting the complementary xpfp zone too.
  // ------------------------------------------------------------------
  vector<vector<Double_t> > xpfpZoneSplit;
  xpfpZoneSplit.resize(info.numFoil);
  for (Int_t nf = 0; nf < info.numFoil; nf++) {
    xpfpZoneSplit[nf].resize(nSlices, -999.0);
    for (Int_t nd = 0; nd < nSlices; nd++) {
      Double_t split = -999.0;
      if (FindXpfpSplitFromPlotName(campaign, tag, nf, nd, info, split)) {
        xpfpZoneSplit[nf][nd] = split;
        cout << "XPFp zone split foil " << nf << " ndel " << nd
             << " = " << split << endl;
      } else {
        cout << "WARNING: no XPFp zone split PDF found for foil " << nf
             << " ndel " << nd
             << "; low/high component cuts will be applied without the raw zone gate." << endl;
      }
    }
  }

  // ------------------------------------------------------------------
  // Load foil ytar/delta cuts.
  // ------------------------------------------------------------------
  vector<TCutG*> ytar_delta_cut(info.numFoil, nullptr);

  TFile* fYtarDeltaCut = nullptr;
  if (requireYtarCut) {
    fYtarDeltaCut = TFile::Open(ytarCutFile, "READ");
    if (!fYtarDeltaCut || fYtarDeltaCut->IsZombie()) {
      cout << "ERROR: cannot open ytar/delta cut file " << ytarCutFile << endl;
      return;
    }

    for (Int_t nf = 0; nf < info.numFoil; nf++) {
      TCutG* c = (TCutG*)fYtarDeltaCut->Get(Form("delta_vs_ytar_cut_foil%d", nf));
      if (c) {
        ytar_delta_cut[nf] = c;
        cout << "Loaded delta_vs_ytar_cut_foil" << nf << " npts=" << c->GetN() << endl;
      } else {
        cout << "WARNING: missing delta_vs_ytar_cut_foil" << nf << endl;
      }
    }
  }

  // ------------------------------------------------------------------
  // Load XFP/XPFP xscol cuts: [foil][delta][xscol][component].
  // Component cuts are preferred. Singleton aliases are used only as fallback.
  // ------------------------------------------------------------------
  vector<vector<vector<vector<XCutComponent> > > > xpfp_xfp_cut;
  xpfp_xfp_cut.resize(info.numFoil);
  for (Int_t nf = 0; nf < info.numFoil; nf++) {
    xpfp_xfp_cut[nf].resize(nSlices);
    for (Int_t nd = 0; nd < nSlices; nd++) {
      xpfp_xfp_cut[nf][nd].resize(9);
    }
  }

  vector<vector<vector<vector<XCutComponent> > > > singletonFallback;
  singletonFallback.resize(info.numFoil);
  for (Int_t nf = 0; nf < info.numFoil; nf++) {
    singletonFallback[nf].resize(nSlices);
    for (Int_t nd = 0; nd < nSlices; nd++) {
      singletonFallback[nf][nd].resize(9);
    }
  }

  TFile* fXCut = TFile::Open(xCutFile, "READ");
  if (!fXCut || fXCut->IsZombie()) {
    cout << "ERROR: cannot open XFP/XPFP cut file " << xCutFile << endl;
    return;
  }

  Int_t nLoadedComponents = 0;
  Int_t nLoadedSingletons = 0;
  set<TString> seenNames;
  TIter nextKey(fXCut->GetListOfKeys());
  TKey* key = nullptr;
  while ((key = (TKey*)nextKey())) {
    TString name = key->GetName();
    if (seenNames.count(name)) continue;
    seenNames.insert(name);

    if (preGmmVetoCuts.count(name.Data())) {
      cout << "VETO X: " << name << endl;
      seenPreGmmVetoCuts.insert(name.Data());
      continue;
    }

    TObject* obj = key->ReadObj();
    if (!obj || !obj->InheritsFrom(TCutG::Class())) continue;

    Int_t xs = -1, nf = -1, nd = -1, part = -1, src = -1;
    TString zone = "";

    if (ParseXComponentCutName(name, xs, nf, nd, part, zone, src)) {
      if (nf < 0 || nf >= info.numFoil || nd < 0 || nd >= nSlices || xs < 0 || xs > 8) continue;
      XCutComponent comp;
      comp.cut = (TCutG*)obj->Clone(Form("xcand_component_x%d_f%d_d%d_p%d", xs, nf, nd, part));
      comp.part = part;
      comp.zone = zone;
      comp.src = src;
      comp.name = name;
      comp.useZoneGate = ((zone == "low" || zone == "high") && xpfpZoneSplit[nf][nd] > -998.0);
      comp.xpfpSplit = xpfpZoneSplit[nf][nd];
      xpfp_xfp_cut[nf][nd][xs].push_back(comp);
      nLoadedComponents++;
      continue;
    }

    if (ParseXSingletonCutName(name, xs, nf, nd)) {
      if (nf < 0 || nf >= info.numFoil || nd < 0 || nd >= nSlices || xs < 0 || xs > 8) continue;
      XCutComponent comp;
      comp.cut = (TCutG*)obj->Clone(Form("xcand_singleton_x%d_f%d_d%d", xs, nf, nd));
      comp.part = -1;
      comp.zone = "global";
      comp.src = xs;
      comp.name = name;
      comp.useZoneGate = kFALSE;
      comp.xpfpSplit = -999.0;
      singletonFallback[nf][nd][xs].push_back(comp);
      nLoadedSingletons++;
      continue;
    }
  }

  Int_t nFallbackUsed = 0;
  for (Int_t nf = 0; nf < info.numFoil; nf++) {
    for (Int_t nd = 0; nd < nSlices; nd++) {
      for (Int_t xs = 0; xs < 9; xs++) {
        if (!xpfp_xfp_cut[nf][nd][xs].empty()) continue;
        if (singletonFallback[nf][nd][xs].empty()) continue;
        xpfp_xfp_cut[nf][nd][xs] = singletonFallback[nf][nd][xs];
        nFallbackUsed += (Int_t)singletonFallback[nf][nd][xs].size();
      }
    }
  }

  cout << "Loaded " << nLoadedComponents << " XFP/XPFP xscol component cuts." << endl;
  cout << "Loaded " << nLoadedSingletons << " singleton/fallback xscol cuts; used " << nFallbackUsed << "." << endl;

  for (const string& cutName : preGmmVetoCuts) {
    if (!seenPreGmmVetoCuts.count(cutName)) {
      cout << "WARNING: configured X pre-GMM veto cut not found: "
           << cutName << endl;
    }
  }

  cout << "X cuts by foil/ndel/xscol:" << endl;
  for (Int_t nf = 0; nf < info.numFoil; nf++) {
    for (Int_t nd = 0; nd < nSlices; nd++) {
      cout << "  foil " << nf << " ndel " << nd << " :";
      for (Int_t xs = 0; xs < 9; xs++) {
        cout << " x" << xs << "=" << xpfp_xfp_cut[nf][nd][xs].size();
      }
      cout << endl;
    }
  }

  // ------------------------------------------------------------------
  // Open replay tree.
  // ------------------------------------------------------------------
  TFile* fin = TFile::Open(inputroot, "READ");
  if (!fin || fin->IsZombie()) {
    cout << "ERROR: cannot open input replay ROOT " << inputroot << endl;
    return;
  }

  TTree* T = (TTree*)fin->Get("T");
  if (!T) {
    cout << "ERROR: cannot find tree T in " << inputroot << endl;
    return;
  }

  // Input branches. Same branch names as make_yscol_candidate_tree.C.
  Double_t sumnpe = 0, etracknorm = 0;
  Double_t ytar = 0, xtar = 0, reactx = 0, reacty = 0, reactz = 0;
  Double_t delta = 0, yptar = 0, xptar = 0;
  Double_t yfp = 0, ypfp = 0, xfp = 0, xpfp = 0;
  Double_t ysieve = 0, xsieve = 0;
  Double_t xbpm_tar = 0, ybpm_tar = 0, frx = 0, fry = 0;

  bool okRequired = true;
  okRequired &= BindBranch(T, "H.cer.npeSum", &sumnpe, true);
  okRequired &= BindBranch(T, "H.cal.etottracknorm", &etracknorm, true);
  okRequired &= BindBranch(T, "H.gtr.y", &ytar, true);
  okRequired &= BindBranch(T, "H.gtr.x", &xtar, true);
  BindBranch(T, "H.react.x", &reactx, false);
  BindBranch(T, "H.react.y", &reacty, false);
  BindBranch(T, "H.react.z", &reactz, false);
  okRequired &= BindBranch(T, "H.gtr.dp", &delta, true);
  okRequired &= BindBranch(T, "H.gtr.ph", &yptar, true);
  okRequired &= BindBranch(T, "H.gtr.th", &xptar, true);
  okRequired &= BindBranch(T, "H.dc.y_fp", &yfp, true);
  okRequired &= BindBranch(T, "H.dc.yp_fp", &ypfp, true);
  okRequired &= BindBranch(T, "H.dc.x_fp", &xfp, true);
  okRequired &= BindBranch(T, "H.dc.xp_fp", &xpfp, true);
  okRequired &= BindBranch(T, "H.extcor.ysieve", &ysieve, true);
  okRequired &= BindBranch(T, "H.extcor.xsieve", &xsieve, true);
  BindBranch(T, "H.rb.raster.fr_xbpm_tar", &xbpm_tar, false);
  BindBranch(T, "H.rb.raster.fr_ybpm_tar", &ybpm_tar, false);
  BindBranch(T, "H.rb.raster.fr_xa", &frx, false);
  BindBranch(T, "H.rb.raster.fr_ya", &fry, false);

  if (!okRequired) {
    cout << "ERROR: stopping because one or more required replay branches are missing." << endl;
    return;
  }

  // ------------------------------------------------------------------
  // Output tree.
  // ------------------------------------------------------------------
  gSystem->mkdir(gSystem->DirName(outputroot), kTRUE);
  gSystem->mkdir(gSystem->DirName(outputtsv), kTRUE);

  TFile* fout = TFile::Open(outputroot, "RECREATE");
  if (!fout || fout->IsZombie()) {
    cout << "ERROR: cannot create output ROOT " << outputroot << endl;
    return;
  }

  TTree* out = new TTree("TXCand", "XFP/XPFP xscol candidate event tree for GMM cleanup");

  Int_t run_out = 0, foil_out = -1, ndel_out = -1, xscol_out = -1;
  Int_t nXMatches_out = 0;
  Int_t nXComponentMatches_out = 0;
  Long64_t entry_out = -1;
  Double_t zfoil_out = 0;
  Double_t delta_low_out = 0, delta_high_out = 0;
  Int_t base_pass_out = 0, ytar_cut_pass_out = 0, xcut_pass_out = 0;

  out->Branch("run", &run_out, "run/I");
  out->Branch("entry", &entry_out, "entry/L");
  out->Branch("foil", &foil_out, "foil/I");
  out->Branch("ndel", &ndel_out, "ndel/I");
  out->Branch("xscol", &xscol_out, "xscol/I");
  out->Branch("nXMatches", &nXMatches_out, "nXMatches/I");
  out->Branch("nXComponentMatches", &nXComponentMatches_out, "nXComponentMatches/I");
  out->Branch("zfoil", &zfoil_out, "zfoil/D");
  out->Branch("delta_low", &delta_low_out, "delta_low/D");
  out->Branch("delta_high", &delta_high_out, "delta_high/D");
  out->Branch("base_pass", &base_pass_out, "base_pass/I");
  out->Branch("ytar_cut_pass", &ytar_cut_pass_out, "ytar_cut_pass/I");
  out->Branch("xcut_pass", &xcut_pass_out, "xcut_pass/I");

  out->Branch("sumnpe", &sumnpe, "sumnpe/D");
  out->Branch("etracknorm", &etracknorm, "etracknorm/D");
  out->Branch("delta", &delta, "delta/D");
  out->Branch("ytar", &ytar, "ytar/D");
  out->Branch("xtar", &xtar, "xtar/D");
  out->Branch("xptar", &xptar, "xptar/D");
  out->Branch("yptar", &yptar, "yptar/D");
  out->Branch("xfp", &xfp, "xfp/D");
  out->Branch("xpfp", &xpfp, "xpfp/D");
  out->Branch("yfp", &yfp, "yfp/D");
  out->Branch("ypfp", &ypfp, "ypfp/D");
  out->Branch("xsieve", &xsieve, "xsieve/D");
  out->Branch("ysieve", &ysieve, "ysieve/D");
  out->Branch("reactx", &reactx, "reactx/D");
  out->Branch("reacty", &reacty, "reacty/D");
  out->Branch("reactz", &reactz, "reactz/D");
  out->Branch("xbpm_tar", &xbpm_tar, "xbpm_tar/D");
  out->Branch("ybpm_tar", &ybpm_tar, "ybpm_tar/D");
  out->Branch("frx", &frx, "frx/D");
  out->Branch("fry", &fry, "fry/D");

  // Overview histograms: one for each foil/delta, all xscols together.
  vector<vector<TH2D*> > hOverview;
  hOverview.resize(info.numFoil);
  for (Int_t nf = 0; nf < info.numFoil; nf++) {
    hOverview[nf].resize(nSlices, nullptr);
    for (Int_t nd = 0; nd < nSlices; nd++) {
      hOverview[nf][nd] = new TH2D(
          Form("overview_before_foil%d_ndel%d", nf, nd),
          Form("%s foil %d delta %.6g to %.6g percent;xsieve;ysieve",
               tag.Data(),
               nf,
               info.delcut[nd],
               info.delcut[nd+1]),
          200, -15, 15, 200, -10, 10);
    }
  }

  vector<vector<vector<Long64_t> > > counts;
  counts.resize(info.numFoil);
  for (Int_t nf = 0; nf < info.numFoil; nf++) {
    counts[nf].resize(nSlices);
    for (Int_t nd = 0; nd < nSlices; nd++) counts[nf][nd].resize(9, 0);
  }

  Long64_t nentries = T->GetEntries();
  Long64_t nBase = 0, nFoilDelta = 0, nFilled = 0, nOverlapEvents = 0, nComponentOverlapEvents = 0;
  cout << "Start candidate loop over " << nentries << " entries." << endl;

  for (Long64_t i = 0; i < nentries; i++) {
    T->GetEntry(i);
    if (i % 50000 == 0) cout << "  entry " << i << endl;

    base_pass_out = (etracknorm > minCalEtot && sumnpe > minCerNpe &&
                     delta > deltaMinGlobal && delta < deltaMaxGlobal) ? 1 : 0;
    if (!base_pass_out) continue;
    nBase++;

    Int_t nf_found = -1;
    if (requireYtarCut) {
      for (Int_t nf = 0; nf < info.numFoil; nf++) {
        if (ytar_delta_cut[nf] && ytar_delta_cut[nf]->IsInside(ytar, delta)) {
          nf_found = nf;
        }
      }
    } else {
      // Fallback mode is only for debugging. In production, keep requireYtarCut=true.
      nf_found = (info.numFoil == 1) ? 0 : -1;
    }
    if (nf_found < 0) continue;

    Int_t nd_found = -1;
    for (Int_t nd = 0; nd < nSlices; nd++) {
      if (delta >= info.delcut[nd] && delta < info.delcut[nd+1]) {
        nd_found = nd;
        break;
      }
    }
    if (nd_found < 0) continue;
    nFoilDelta++;

    vector<Int_t> xMatches;
    vector<Int_t> xComponentMatchCounts;
    for (Int_t xs = 0; xs < 9; xs++) {
      Int_t nCompMatchThisXs = 0;
      const vector<XCutComponent>& comps = xpfp_xfp_cut[nf_found][nd_found][xs];
      for (size_t ic = 0; ic < comps.size(); ic++) {
        const XCutComponent& comp = comps[ic];
        TCutG* cx = comp.cut;
        if (!cx || !cx->IsInside(xpfp, xfp)) continue;

        if (comp.useZoneGate) {
          if (comp.zone == "low"  && !(xpfp <  comp.xpfpSplit)) continue;
          if (comp.zone == "high" && !(xpfp >= comp.xpfpSplit)) continue;
        }

        nCompMatchThisXs++;
      }
      if (nCompMatchThisXs > 0) {
        xMatches.push_back(xs);
        xComponentMatchCounts.push_back(nCompMatchThisXs);
        if (nCompMatchThisXs > 1) nComponentOverlapEvents++;
      }
    }

    if (xMatches.size() == 0) continue;
    if (xMatches.size() > 1) nOverlapEvents++;

    nXMatches_out = (Int_t)xMatches.size();

    for (size_t im = 0; im < xMatches.size(); im++) {
      run_out = nrun;
      entry_out = i;
      foil_out = nf_found;
      ndel_out = nd_found;
      xscol_out = xMatches[im];
      nXComponentMatches_out = xComponentMatchCounts[im];
      zfoil_out = info.zfoil[nf_found];
      delta_low_out = info.delcut[nd_found];
      delta_high_out = info.delcut[nd_found+1];
      ytar_cut_pass_out = 1;
      xcut_pass_out = 1;

      out->Fill();
      hOverview[nf_found][nd_found]->Fill(xsieve, ysieve);
      counts[nf_found][nd_found][xscol_out]++;
      nFilled++;
    }
  }

  cout << "Done." << endl;
  cout << "  baseline-pass events          = " << nBase << endl;
  cout << "  foil+delta matched events     = " << nFoilDelta << endl;
  cout << "  candidate rows written        = " << nFilled << endl;
  cout << "  events with overlapping X cols = " << nOverlapEvents << endl;
  cout << "  events with overlapping X components inside same xscol = " << nComponentOverlapEvents << endl;

  cout << "Counts by foil/ndel/xscol:" << endl;
  for (Int_t nf = 0; nf < info.numFoil; nf++) {
    for (Int_t nd = 0; nd < nSlices; nd++) {
      cout << "  foil " << nf << " ndel " << nd << " :";
      for (Int_t xs = 0; xs < 9; xs++) cout << " x" << xs << "=" << counts[nf][nd][xs];
      cout << endl;
    }
  }

  fout->cd();
  out->Write();
  for (Int_t nf = 0; nf < info.numFoil; nf++) {
    for (Int_t nd = 0; nd < nSlices; nd++) {
      hOverview[nf][nd]->Write();
    }
  }

  // Save a compact TSV summary in the campaign TSV directory.
  ofstream ofs(outputtsv.Data());

  if (!ofs.is_open()) {
    cout << "ERROR: cannot create output TSV " << outputtsv << endl;
    fout->Close();
    return;
  }

  ofs << "rungroup\trun\tfoil\tndel\txscol\tN\t"
         "delta_low\tdelta_high\tzfoil\n";

  for (Int_t nf = 0; nf < info.numFoil; nf++) {
    for (Int_t nd = 0; nd < nSlices; nd++) {
      for (Int_t xs = 0; xs < 9; xs++) {
        ofs << tag << "\t"
            << nrun << "\t"
            << nf << "\t"
            << nd << "\t"
            << xs << "\t"
            << counts[nf][nd][xs] << "\t"
            << info.delcut[nd] << "\t"
            << info.delcut[nd+1] << "\t"
            << info.zfoil[nf] << "\n";
      }
    }
  }

  ofs.close();

  fout->Close();
  cout << "Wrote " << outputroot << endl;
  cout << "Wrote " << outputtsv << endl;
}

