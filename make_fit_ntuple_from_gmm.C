#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TSystem.h>
#include <TMath.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TStyle.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <map>
#include <string>
#include <algorithm>
#include <cmath>

using namespace std;

struct OpticsInfoGMM {
  Int_t run = -1;
  TString opticsID = "";
  Double_t centAngleDeg = 0.0;
  Int_t numFoil = 0;
  Int_t ndelcut = 0;       // number of delta edges in DB file
  vector<Double_t> zfoil;
  vector<Double_t> delcut; // delta edges
  Bool_t ok = kFALSE;
};

struct ColMaskInfo {
  Int_t foil = -1;
  Int_t ndel = -1;
  Int_t col = -1;          // yscol or xscol
  Double_t score = -1e99;  // use best score if duplicate entry appears
};

static string trim_gmm(const string &s);

struct VetoRuleGMM {
  Int_t run = -1;
  Int_t foil = -1;
  Int_t ndel = -1;
  TString axis = "";
  Int_t col = -1;
  TString reason = "";
};

static TString lower_gmm(TString x) {
  x.ToLower();
  return x;
}

static vector<VetoRuleGMM> read_veto_rules_gmm(TString vetoFile) {
  vector<VetoRuleGMM> rules;

  if (vetoFile.Length() == 0) return rules;
  if (gSystem->AccessPathName(vetoFile.Data())) {
    cout << "No veto file found: " << vetoFile << " ; using no vetoes." << endl;
    return rules;
  }

  ifstream fin(vetoFile.Data());
  string line;
  while (getline(fin,line)) {
    string t = trim_gmm(line);
    if (t.size()==0) continue;
    if (t[0]=='#') continue;

    stringstream ss(t);
    VetoRuleGMM r;
    string axis, reason;
    ss >> r.run >> r.foil >> r.ndel >> axis >> r.col;
    getline(ss, reason);

    if (!ss && axis.size()==0) continue;

    r.axis = lower_gmm(axis.c_str());
    r.reason = trim_gmm(reason).c_str();

    if (!(r.axis=="yscol" || r.axis=="y" ||
          r.axis=="xscol" || r.axis=="x" ||
          r.axis=="both")) {
      cout << "WARNING: bad veto axis in line: " << t << endl;
      continue;
    }

    rules.push_back(r);
  }

  cout << "Loaded " << rules.size() << " veto rules from " << vetoFile << endl;
  for (auto &r : rules) {
    cout << "  veto run=" << r.run
         << " foil=" << r.foil
         << " ndel=" << r.ndel
         << " axis=" << r.axis
         << " col=" << r.col
         << " reason=" << r.reason << endl;
  }

  return rules;
}

static Bool_t veto_match_int_gmm(Int_t ruleVal, Int_t value) {
  return (ruleVal < 0 || ruleVal == value);
}

static Bool_t is_vetoed_gmm(const vector<VetoRuleGMM> &rules,
                            Int_t run, Int_t foil, Int_t ndel,
                            Int_t yscol, Int_t xscol) {
  for (auto &r : rules) {
    if (!veto_match_int_gmm(r.run, run)) continue;
    if (!veto_match_int_gmm(r.foil, foil)) continue;
    if (!veto_match_int_gmm(r.ndel, ndel)) continue;

    TString ax = lower_gmm(r.axis);

    if (ax=="yscol" || ax=="y") {
      if (veto_match_int_gmm(r.col, yscol)) return kTRUE;
    } else if (ax=="xscol" || ax=="x") {
      if (veto_match_int_gmm(r.col, xscol)) return kTRUE;
    } else if (ax=="both") {
      if (veto_match_int_gmm(r.col, yscol) ||
          veto_match_int_gmm(r.col, xscol)) return kTRUE;
    }
  }
  return kFALSE;
}


static string trim_gmm(const string &s) {
  const char *ws = " \t\r\n";
  size_t b = s.find_first_not_of(ws);
  if (b == string::npos) return "";
  size_t e = s.find_last_not_of(ws);
  return s.substr(b, e-b+1);
}

static vector<string> split_csv_gmm(const string &line) {
  vector<string> out;
  string item;
  stringstream ss(line);
  while (getline(ss, item, ',')) out.push_back(trim_gmm(item));
  return out;
}

static Bool_t read_optics_info_gmm(Int_t nrun, const TString &opticsFile, OpticsInfoGMM &info) {
  ifstream fin(opticsFile.Data());
  if (!fin.is_open()) {
    cout << "ERROR: cannot open optics DB: " << opticsFile << endl;
    return kFALSE;
  }

  string line;
  while (getline(fin, line)) {
    string t = trim_gmm(line);
    if (t.size() == 0) continue;
    vector<string> fields = split_csv_gmm(t);
    if (fields.size() < 6) continue;
    if (atoi(fields[0].c_str()) != nrun) continue;

    info.run = nrun;
    info.opticsID = fields[1].c_str();
    info.centAngleDeg = atof(fields[2].c_str());
    info.numFoil = atoi(fields[3].c_str());
    info.ndelcut = atoi(fields[5].c_str());

    string foilLine, delLine;
    if (!getline(fin, foilLine) || !getline(fin, delLine)) {
      cout << "ERROR: optics DB ended while reading foil/delta lines for run " << nrun << endl;
      return kFALSE;
    }

    vector<string> fz = split_csv_gmm(foilLine);
    vector<string> dz = split_csv_gmm(delLine);
    info.zfoil.clear();
    info.delcut.clear();
    for (auto &x : fz) if (x.size()) info.zfoil.push_back(atof(x.c_str()));
    for (auto &x : dz) if (x.size()) info.delcut.push_back(atof(x.c_str()));

    if ((Int_t)info.zfoil.size() != info.numFoil) {
      cout << "WARNING: run " << nrun << " DB NumFoil=" << info.numFoil
           << " but read " << info.zfoil.size() << " foil z positions." << endl;
    }
    if (info.delcut.size() < 2) {
      cout << "ERROR: run " << nrun << " has fewer than 2 delta edges." << endl;
      return kFALSE;
    }

    info.ok = kTRUE;
    return kTRUE;
  }

  cout << "ERROR: did not find run " << nrun << " in " << opticsFile << endl;
  return kFALSE;
}

static Long64_t read_gmm_mask_one_file(const TString &fname,
                                       Bool_t isY,
                                       map<Long64_t, ColMaskInfo> &mask,
                                       Bool_t verbose=kFALSE) {
  if (gSystem->AccessPathName(fname.Data())) return 0;

  TFile *f = TFile::Open(fname.Data(), "READ");
  if (!f || f->IsZombie()) {
    cout << "WARNING: could not open GMM file: " << fname << endl;
    if (f) f->Close();
    return 0;
  }

  TTree *t = (TTree*)f->Get("GMMClean");
  if (!t) {
    cout << "WARNING: no GMMClean tree in " << fname << endl;
    f->Close();
    return 0;
  }

  Int_t entry=0, foil=-1, ndel=-1, col=-1, gmm_keep=0;
  Double_t gmm_score=-1e99;
  t->SetBranchAddress("entry", &entry);
  t->SetBranchAddress("foil", &foil);
  t->SetBranchAddress("ndel", &ndel);
  if (isY) t->SetBranchAddress("yscol", &col);
  else     t->SetBranchAddress("xscol", &col);
  t->SetBranchAddress("gmm_keep", &gmm_keep);
  if (t->GetBranch("gmm_score")) t->SetBranchAddress("gmm_score", &gmm_score);

  Long64_t nkeep = 0;
  Long64_t nent = t->GetEntries();
  for (Long64_t i=0; i<nent; i++) {
    t->GetEntry(i);
    if (gmm_keep != 1) continue;
    ColMaskInfo rec;
    rec.foil = foil;
    rec.ndel = ndel;
    rec.col = col;
    rec.score = gmm_score;

    Long64_t key = (Long64_t)entry;
    auto it = mask.find(key);
    if (it == mask.end() || rec.score > it->second.score) mask[key] = rec;
    nkeep++;
  }

  if (verbose) cout << "Read " << nkeep << " kept entries from " << fname << endl;
  f->Close();
  return nkeep;
}


static TString delta_tag_gmm(Int_t ndel) {
  switch (ndel) {
    case 0: return "m10_to_m8";
    case 1: return "m8_to_m5";
    case 2: return "m5_to_0";
    case 3: return "0_to_5";
    case 4: return "5_to_10";
    default: return "UNKNOWN";
  }
}

static Long64_t read_gmm_masks_for_run(Int_t run,
                                        const TString &rungroup,
                                        Int_t nFoils,
                                        Int_t nDeltaIntervals,
                                        const TString &yBaseDir,
                                        const TString &xBaseDir,
                                        const TString &xTag,
                                        map<Long64_t, ColMaskInfo> &ymask,
                                        map<Long64_t, ColMaskInfo> &xmask,
                                        Bool_t verbose=kFALSE) {
  ymask.clear();
  xmask.clear();
  Long64_t ny=0, nx=0;

  for (Int_t nf=0; nf<nFoils; nf++) {
    for (Int_t nd=0; nd<nDeltaIntervals; nd++) {
      TString dtag = delta_tag_gmm(nd);
      if (dtag == "UNKNOWN") continue;

      TString yfile = Form("%s/gmm_clean_%s_foil%d_delta_%s_allY.root",
                           yBaseDir.Data(), rungroup.Data(), nf, dtag.Data());
      TString xfile = Form("%s/gmm_clean_%s_foil%d_delta_%s_allX.root",
                           xBaseDir.Data(), rungroup.Data(), nf, dtag.Data());
      ny += read_gmm_mask_one_file(yfile, kTRUE,  ymask, verbose);
      nx += read_gmm_mask_one_file(xfile, kFALSE, xmask, verbose);
    }
  }

  cout << "GMM mask summary run " << run
       << ": unique kept Y entries=" << ymask.size()
       << " unique kept X entries=" << xmask.size()
       << " (raw kept rows Y=" << ny << ", X=" << nx << ")" << endl;
  return TMath::Min((Long64_t)ymask.size(), (Long64_t)xmask.size());
}

void plot_fit_ntuple_qa_from_existing(
                              Int_t nrun=611701,
                              Int_t FileID=-1,
                              TString xTag="6p117_noVeto",
                              TString outputDir="HMS_6p117GeV/06a_fit_ntuple");

void make_fit_ntuple_from_gmm(
                              Int_t nrun=611701,
                              Int_t FileID=-1,
                              TString xTag="6p117_noVeto",
                              TString rootDir="/volatile/hallc/nps/cploen/ROOTfiles/OPTICS/angular_sandbox/6p117_from_6p667AutoMatrix_no_offsets_20260630/rootfiles/hadd_rungroups_6p117",
                              TString yBaseDir="HMS_6p117GeV/05a_gmm_cleanup_y/root",
                              TString xBaseDir="HMS_6p117GeV/05b_gmm_cleanup_x/root",
                              Double_t cerCut=6.0,
                              Double_t calCut=0.65,
                              Bool_t savePlots=kTRUE,
                              TString vetoFile="",
                              TString outputDir="HMS_6p117GeV/06a_fit_ntuple",
                              TString inputRootOverride="",
                              TString rungroup="") {
  gStyle->SetOptStat(0);

  OpticsInfoGMM info;
  if (!read_optics_info_gmm(nrun, "DATfiles/list_of_optics_run.dat", info)) return;

  const Int_t nDeltaIntervals = (Int_t)info.delcut.size() - 1;
  cout << "Parsed run " << nrun
       << " OpticsID=" << info.opticsID
       << " angle=" << info.centAngleDeg
       << " NumFoil=" << info.numFoil
       << " delta intervals=" << nDeltaIntervals << endl;

  if (rungroup.Length() == 0) {
    cout << "ERROR: rungroup is required for GMM input filenames." << endl;
    return;
  }

  map<Long64_t, ColMaskInfo> ymask, xmask;
  read_gmm_masks_for_run(nrun, rungroup,
                         info.numFoil, nDeltaIntervals,
                         yBaseDir, xBaseDir, xTag,
                         ymask, xmask, kFALSE);

  vector<VetoRuleGMM> vetoRules = read_veto_rules_gmm(vetoFile);

  if (ymask.empty() || xmask.empty()) {
    cout << "ERROR: empty Y or X GMM mask for run " << nrun << "." << endl;
    cout << "Check the rungroup and Y/X GMM input directories." << endl;
    return;
  }

  TString inputroot;
  if (inputRootOverride.Length() > 0) {
    inputroot = inputRootOverride;
  } else {
    inputroot = Form("%s/nps_hms_optics_6p117_%s.root",
                     rootDir.Data(), info.opticsID.Data());
  }

  cout << "Input ROOT: " << inputroot << endl;

  gSystem->mkdir(outputDir, kTRUE);
  gSystem->mkdir(Form("%s/root", outputDir.Data()), kTRUE);
  gSystem->mkdir(Form("%s/plots", outputDir.Data()), kTRUE);

  TString outputroot = Form("%s/root/Optics_%d_%d_fit_tree_gmm.root",
                            outputDir.Data(), nrun, FileID);
  cout << "Input ROOT : " << inputroot << endl;
  cout << "Output TFit: " << outputroot << endl;

  TFile *fin = TFile::Open(inputroot.Data(), "READ");
  if (!fin || fin->IsZombie()) {
    cout << "ERROR: cannot open input ROOT file: " << inputroot << endl;
    if (fin) fin->Close();
    return;
  }
  TTree *T = (TTree*)fin->Get("T");
  if (!T) {
    cout << "ERROR: cannot find tree T in " << inputroot << endl;
    fin->Close();
    return;
  }

  Double_t sumnpe=0, etracknorm=0;
  Double_t ytar=0, xtar=0, reactx=0, reacty=0, reactz=0;
  Double_t delta=0, yptar=0, xptar=0;
  Double_t yfp=0, ypfp=0, xfp=0, xpfp=0;
  Double_t ysieve=0, xsieve=0, xbpm_tar=0, ybpm_tar=0;

  T->SetBranchAddress("H.cer.npeSum", &sumnpe);
  T->SetBranchAddress("H.cal.etottracknorm", &etracknorm);
  T->SetBranchAddress("H.gtr.y", &ytar);
  T->SetBranchAddress("H.gtr.x", &xtar);
  T->SetBranchAddress("H.react.x", &reactx);
  T->SetBranchAddress("H.react.y", &reacty);
  T->SetBranchAddress("H.react.z", &reactz);
  T->SetBranchAddress("H.gtr.dp", &delta);
  T->SetBranchAddress("H.gtr.ph", &yptar);
  T->SetBranchAddress("H.gtr.th", &xptar);
  T->SetBranchAddress("H.dc.y_fp", &yfp);
  T->SetBranchAddress("H.dc.yp_fp", &ypfp);
  T->SetBranchAddress("H.dc.x_fp", &xfp);
  T->SetBranchAddress("H.dc.xp_fp", &xpfp);
  T->SetBranchAddress("H.extcor.ysieve", &ysieve);
  T->SetBranchAddress("H.extcor.xsieve", &xsieve);
  T->SetBranchAddress("H.rb.raster.fr_xbpm_tar", &xbpm_tar);
  T->SetBranchAddress("H.rb.raster.fr_ybpm_tar", &ybpm_tar);

  TFile fout(outputroot.Data(), "RECREATE");
  TTree *otree = new TTree("TFit", "FitTree from GMM-cleaned yscol/xscol intersection");

  Double_t xptarT=0, ytarT=0, yptarT=0, ysieveT=0, xsieveT=0, ztarT=0, ztar=0, xtarT=0;
  Int_t foilT=-1, ndelT=-1, yscolT=-1, xscolT=-1;
  Long64_t entryT=-1;

  otree->Branch("entry", &entryT, "entry/L");
  otree->Branch("foil", &foilT, "foil/I");
  otree->Branch("ndel", &ndelT, "ndel/I");
  otree->Branch("yscol", &yscolT, "yscol/I");
  otree->Branch("xscol", &xscolT, "xscol/I");

  otree->Branch("ys", &ysieve);
  otree->Branch("ysT", &ysieveT);
  otree->Branch("xs", &xsieve);
  otree->Branch("xsT", &xsieveT);
  otree->Branch("ztarT", &ztarT);
  otree->Branch("ztar", &ztar);
  otree->Branch("xtar", &xtar);
  otree->Branch("xtarT", &xtarT);
  otree->Branch("xptar", &xptar);
  otree->Branch("yptar", &yptar);
  otree->Branch("ytar", &ytar);
  otree->Branch("xptarT", &xptarT);
  otree->Branch("yptarT", &yptarT);
  otree->Branch("ytarT", &ytarT);
  otree->Branch("delta", &delta);
  otree->Branch("xpfp", &xpfp);
  otree->Branch("ypfp", &ypfp);
  otree->Branch("xfp", &xfp);
  otree->Branch("yfp", &yfp);
  otree->Branch("reactxcalc", &reactx);
  otree->Branch("reactycalc", &reacty);

  vector<Double_t> ys_cent;
  vector<Double_t> xs_cent;
  for (Int_t i=0; i<9; i++) {
    ys_cent.push_back((i-4)*0.6*2.54);
    xs_cent.push_back((i-4)*2.54);
  }

  Double_t centAngle = info.centAngleDeg*TMath::Pi()/180.0;
  Double_t y_mis;
  Double_t x_mis;
  if (TMath::Abs(info.centAngleDeg)<40) y_mis = 0.1*(0.52-0.012*TMath::Abs(info.centAngleDeg)+0.002*TMath::Abs(info.centAngleDeg)*TMath::Abs(info.centAngleDeg));
  else y_mis = 0.1*(0.52-0.012*40. + 0.002*40.*40.);
  if (TMath::Abs(info.centAngleDeg)<50) x_mis = 0.1*(2.37-0.086*TMath::Abs(info.centAngleDeg)+0.0012*TMath::Abs(info.centAngleDeg)*TMath::Abs(info.centAngleDeg));
  else x_mis = 0.1*(2.37-0.086*50.+0.0012*50.*50.);

  const Double_t zdis_sieve = 168.0;

  TH2F *hTrueSieve = new TH2F("hTrueSieve", Form("run %d accepted GMM targets; ysieveT (cm); xsieveT (cm)", nrun), 80, -7, 7, 100, -13, 13);
  TH2F *hRecoSieve = new TH2F("hRecoSieve", Form("run %d accepted reconstructed sieve; ysieve (cm); xsieve (cm)", nrun), 80, -7, 7, 100, -13, 13);
  TH1F *hDelta = new TH1F("hDelta", Form("run %d accepted delta; delta; counts", nrun), 80, -10, 10);
  TH1F *hFoil = new TH1F("hFoil", Form("run %d accepted foil; foil; counts", nrun), 5, -0.5, 4.5);

  Long64_t nentries = T->GetEntries();
  Long64_t nBoth=0, nMismatch=0, nPIDFail=0, nVeto=0, nFilled=0;

  for (Long64_t i=0; i<nentries; i++) {
    auto iy = ymask.find(i);
    if (iy == ymask.end()) continue;
    auto ix = xmask.find(i);
    if (ix == xmask.end()) continue;
    nBoth++;

    if (iy->second.foil != ix->second.foil || iy->second.ndel != ix->second.ndel) {
      nMismatch++;
      continue;
    }
    if (iy->second.col < 0 || iy->second.col >= 9 || ix->second.col < 0 || ix->second.col >= 9) continue;
    if (iy->second.foil < 0 || iy->second.foil >= (Int_t)info.zfoil.size()) continue;

    T->GetEntry(i);
    if (!(sumnpe > cerCut && etracknorm > calCut && delta > -10.0 && delta < 10.0)) {
      nPIDFail++;
      continue;
    }

    entryT = i;
    foilT = iy->second.foil;
    ndelT = iy->second.ndel;
    yscolT = iy->second.col;
    xscolT = ix->second.col;

    if (is_vetoed_gmm(vetoRules, nrun, foilT, ndelT, yscolT, xscolT)) {
      nVeto++;
      continue;
    }

    Double_t zf = info.zfoil[foilT];
    Double_t xbeam = -xbpm_tar;
    Double_t ybeam =  ybpm_tar; // retained for possible future diagnostics
    (void)ybeam;

    Double_t ytar_cent = zf*TMath::Sin(centAngle) + xbeam*TMath::Cos(centAngle) - y_mis;
    yptarT = (ys_cent[yscolT] - ytar_cent) / (zdis_sieve - zf*TMath::Cos(centAngle));
    ytarT  = zf*(TMath::Sin(centAngle) - yptarT*TMath::Cos(centAngle))
           + xbeam*(TMath::Cos(centAngle) + yptarT*TMath::Sin(centAngle)) - y_mis;
    xptarT = xs_cent[xscolT] / (zdis_sieve - zf*TMath::Cos(centAngle));
    xtarT  = -reacty - x_mis - xptarT*zf*TMath::Cos(centAngle);

    ysieveT = ys_cent[yscolT];
    xsieveT = xs_cent[xscolT];
    ztarT = zf;
    ztar = reactz;

    otree->Fill();
    hTrueSieve->Fill(ysieveT, xsieveT);
    hRecoSieve->Fill(ysieve, xsieve);
    hDelta->Fill(delta);
    hFoil->Fill(foilT);
    nFilled++;
  }

  cout << "Run " << nrun << " fit ntuple summary" << endl;
  cout << "  entries in both Y and X masks : " << nBoth << endl;
  cout << "  foil/ndel mismatches          : " << nMismatch << endl;
  cout << "  PID fails after mask          : " << nPIDFail << endl;
  cout << "  events vetoed by HIP rules    : " << nVeto << endl;
  cout << "  TFit rows written             : " << nFilled << endl;

  otree->Write();
  hTrueSieve->Write();
  hRecoSieve->Write();
  hDelta->Write();
  hFoil->Write();

  fout.Close();
  fin->Close();

  if (savePlots) {
    plot_fit_ntuple_qa_from_existing(nrun, FileID, xTag, outputDir);
  }
}

void plot_fit_ntuple_qa_from_existing(
                              Int_t nrun,
                              Int_t FileID,
                              TString xTag,
                              TString outputDir) {
  gStyle->SetOptStat(0);

  OpticsInfoGMM info;
  if (!read_optics_info_gmm(
        nrun, "DATfiles/list_of_optics_run.dat", info)) {
    return;
  }

  TString inputroot =
    Form("%s/root/Optics_%d_%d_fit_tree_gmm.root",
         outputDir.Data(), nrun, FileID);

  TFile *f = TFile::Open(inputroot.Data(), "READ");
  if (!f || f->IsZombie()) {
    cout << "ERROR: cannot open existing fit-ntuple ROOT file: "
         << inputroot << endl;
    if (f) {
      f->Close();
      delete f;
    }
    return;
  }

  TH2F *hTrueSource =
    dynamic_cast<TH2F*>(f->Get("hTrueSieve"));
  TH2F *hRecoSource =
    dynamic_cast<TH2F*>(f->Get("hRecoSieve"));
  TH1F *hDeltaSource =
    dynamic_cast<TH1F*>(f->Get("hDelta"));
  TH1F *hFoilSource =
    dynamic_cast<TH1F*>(f->Get("hFoil"));

  if (!hTrueSource || !hRecoSource ||
      !hDeltaSource || !hFoilSource) {
    cout << "ERROR: one or more QA histograms are missing from "
         << inputroot << endl;
    f->ls();
    f->Close();
    delete f;
    return;
  }

  TH2F *hTrueSieve =
    dynamic_cast<TH2F*>(
      hTrueSource->Clone(Form("hTrueSieve_plot_%d", nrun)));
  TH2F *hRecoSieve =
    dynamic_cast<TH2F*>(
      hRecoSource->Clone(Form("hRecoSieve_plot_%d", nrun)));
  TH1F *hDelta =
    dynamic_cast<TH1F*>(
      hDeltaSource->Clone(Form("hDelta_plot_%d", nrun)));
  TH1F *hFoil =
    dynamic_cast<TH1F*>(
      hFoilSource->Clone(Form("hFoil_plot_%d", nrun)));

  hTrueSieve->SetDirectory(nullptr);
  hRecoSieve->SetDirectory(nullptr);
  hDelta->SetDirectory(nullptr);
  hFoil->SetDirectory(nullptr);

  f->Close();
  delete f;

  TString plotDir = Form("%s/plots", outputDir.Data());
  gSystem->mkdir(plotDir.Data(), kTRUE);

  TString pdf =
    Form("%s/gmm_fit_ntuple_%s_%s.pdf",
         plotDir.Data(), info.opticsID.Data(), xTag.Data());

  TCanvas *c =
    new TCanvas(Form("c_gmm_fit_ntuple_%d", nrun),
                Form("GMM fit ntuple QA: %s",
                     info.opticsID.Data()),
                1200, 900);

  c->Divide(2,2);

  c->cd(1);
  hTrueSieve->Draw("COLZ");

  c->cd(2);
  hRecoSieve->Draw("COLZ");

  c->cd(3);
  hDelta->Draw();

  c->cd(4);
  hFoil->Draw();

  c->Modified();
  c->Update();
  c->Print(pdf.Data());

  if (gSystem->AccessPathName(pdf.Data())) {
    cout << "ERROR: PDF was not created: " << pdf << endl;
  } else {
    cout << "Saved fit-ntuple QA PDF: " << pdf << endl;
  }

  delete c;
  delete hTrueSieve;
  delete hRecoSieve;
  delete hDelta;
  delete hFoil;
}


void plot_fit_ntuple_qa_all_existing(
                              Int_t FileID=-1,
                              TString xTag="6p117_noVeto",
                              TString outputDir="HMS_6p117GeV/06a_fit_ntuple") {
  vector<Int_t> runs = {
    611701, 611702, 611703, 611704,
    611705, 611706, 611707, 611708
  };

  for (auto run : runs) {
    cout << "============================================================"
         << endl;
    cout << "Plotting existing fit ntuple for metadata run "
         << run << endl;
    cout << "============================================================"
         << endl;

    plot_fit_ntuple_qa_from_existing(
      run, FileID, xTag, outputDir
    );
  }
}


void make_fit_ntuples_from_gmm_all(
                                   TString xTag="6p117_noVeto",
                                   Int_t FileID=-1,
                                   TString rootDir="/volatile/hallc/nps/cploen/ROOTfiles/OPTICS/angular_sandbox/6p117_from_6p667AutoMatrix_no_offsets_20260630/rootfiles/hadd_rungroups_6p117",
                                   TString yBaseDir="HMS_6p117GeV/05a_gmm_cleanup_y/root",
                                   TString xBaseDir="HMS_6p117GeV/05b_gmm_cleanup_x/root",
                                   TString outputDir="HMS_6p117GeV/06a_fit_ntuple",
                                   TString vetoFile="") {
  vector<Int_t> runs = {
    611701, 611702, 611703, 611704,
    611705, 611706, 611707, 611708
  };

  for (auto run : runs) {
    cout << "============================================================" << endl;
    cout << "Building GMM TFit ntuple for metadata run " << run << endl;
    cout << "============================================================" << endl;

    make_fit_ntuple_from_gmm(
      run, FileID, xTag,
      rootDir, yBaseDir, xBaseDir,
      6.0, 0.65, kTRUE,
      vetoFile, outputDir
    );
  }
}
