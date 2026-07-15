// make_yscol_candidate_tree.C
//
// Bridge macro for YFP/YPFP -> (xsieve,ysieve) GMM cleanup.
//
// Purpose:
//   Apply the established foil/delta ytar cuts and YFP/YPFP yscol cuts
//   to a replay ROOT tree, then write a flat candidate tree containing
//   one row per accepted candidate event/yscol assignment.
//
// This intentionally does NOT require XFP/XPFP row cuts yet.
// It is the staging tree for Python/scikit-learn GMM cleanup in
// projected sieve space: (xsieve, ysieve).
//
// Usage:
//   make_yscol_candidate_tree(
//     metadataRun,
//     rungroupTag,
//     campaignDir,
//     inputRootFile
//   )
//
// The remaining arguments optionally override the derived cut/output paths.
//
// Output:
//   TTree "TYCand" with event-level candidate rows.
//   TH2D overview_before_foil*_ndel* histograms of xsieve vs ysieve by slice.
//
// Notes:
//   - User-facing yscol is only the final assigned yscol label 0..8.
//   - No old_yscol branch is kept.
//   - If Y cuts overlap, the same event can appear more than once with
//     different yscol. nYMatches is stored for diagnosing that case.

#include <TSystem.h>
#include <TString.h>
#include <TFile.h>
#include <TKey.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TCutG.h>
#include <TROOT.h>
#include <TMath.h>
#include <TStyle.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <set>
#include <sstream>
#include <string>

using namespace std;

void make_yscol_candidate_tree(Int_t nrun=1814,
                               TString rungroupTag="",
                               TString campaignDir="",
                               TString inputRootFile="",
                               TString ytarCutOverride="",
                               TString yCutOverride="",
                               TString outputRootOverride="",
                               TString outputTsvOverride="") {
  Bool_t CutYtarFlag = kTRUE;
  Bool_t CutYpFpYFpFlag = kTRUE;

  gStyle->SetPalette(1,0);
  gStyle->SetOptStat(1000011);

  // ------------------------------------------------------------------
  // Read run metadata from DATfiles/list_of_optics_run.dat
  // Mirrors make_fit_ntuple.C so this stays compatible with the optics workflow.
  // ------------------------------------------------------------------
  TString OpticsFile = Form("DATfiles/list_of_optics_run.dat");
  ifstream file_optics(OpticsFile.Data());
  TString OpticsID="";
  Int_t RunNum=0;
  Double_t CentAngleDeg=0.0;
  Int_t SieveFlag=1;
  Double_t ymis=0.0;
  Int_t NumFoil=0;
  Int_t ndelcut=0;
  TString temp;

  vector<Double_t> ztar_foil;
  vector<Double_t> delcut;

  if (!file_optics.is_open()) {
    cout << "ERROR: cannot open " << OpticsFile << endl;
    return;
  }

  cout << "Open run-list file = " << OpticsFile << endl;
  while (RunNum != nrun && !file_optics.eof()) {
    temp.ReadToDelim(file_optics, ',');
    if (temp.Atoi() == nrun) {
      RunNum = temp.Atoi();
    } else {
      temp.ReadLine(file_optics);
    }
  }

  if (RunNum != nrun) {
    cout << "ERROR: run " << nrun << " not found in " << OpticsFile << endl;
    return;
  }

  temp.ReadToDelim(file_optics, ','); OpticsID = temp;
  temp.ReadToDelim(file_optics, ','); CentAngleDeg = temp.Atof();
  temp.ReadToDelim(file_optics, ','); NumFoil = temp.Atoi();
  temp.ReadToDelim(file_optics, ','); SieveFlag = temp.Atoi();
  temp.ReadToDelim(file_optics, ','); ndelcut = temp.Atoi();
  temp.ReadToDelim(file_optics);      ymis = temp.Atof();

  for (Int_t nf=0; nf<NumFoil-1; nf++) {
    temp.ReadToDelim(file_optics, ',');
    ztar_foil.push_back(temp.Atof());
  }
  temp.ReadToDelim(file_optics);
  ztar_foil.push_back(temp.Atof());

  // ndelcut in list_of_optics_run.dat is the number of delta boundaries.
  // Therefore the number of physical delta slices is ndelcut - 1.
  for (Int_t nd=0; nd<ndelcut-1; nd++) {
    temp.ReadToDelim(file_optics, ',');
    delcut.push_back(temp.Atof());
  }
  temp.ReadToDelim(file_optics);
  delcut.push_back(temp.Atof());

  const Int_t nSlices = static_cast<Int_t>(delcut.size()) - 1;

  cout << "Run " << RunNum
       << " OpticsID=" << OpticsID
       << " CentAngleDeg=" << CentAngleDeg
       << " NumFoil=" << NumFoil
       << " nSlices=" << nSlices << endl;

  if (NumFoil <= 0 || nSlices <= 0) {
    cout << "ERROR: NumFoil or nSlices invalid." << endl;
    return;
  }

  for (Int_t nf=0; nf<NumFoil; nf++) cout << "  foil " << nf << " ztar=" << ztar_foil[nf] << endl;
  for (Int_t nd=0; nd<nSlices; nd++) cout << "  ndel " << nd << " delta=[" << delcut[nd] << "," << delcut[nd+1] << ")" << endl;

  // ------------------------------------------------------------------
  // Campaign-generic file names.
  // nrun remains only as the legacy metadata lookup key.
  // ------------------------------------------------------------------
  TString tag = rungroupTag;
  if (tag.Length() == 0) tag = OpticsID;

  TString campaign = campaignDir;
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

  TString ytarCutFile = (ytarCutOverride.Length()>0)
                      ? ytarCutOverride
                      : Form("%s/01_ytar_cuts/cuts/ytar_ridge_cut_%s.root",
                             campaign.Data(), tag.Data());

  TString yCutFile = (yCutOverride.Length()>0)
                   ? yCutOverride
                   : Form("%s/03a_relabel_y/cuts/YpFpYFp_%s_auto_band_cut_yscolRelabeled.root",
                          campaign.Data(), tag.Data());

  TString outputroot = (outputRootOverride.Length()>0)
                     ? outputRootOverride
                     : Form("%s/04a_candidate_trees_y/root/YscolCandidates_%s.root",
                            campaign.Data(), tag.Data());

  TString outputtsv = (outputTsvOverride.Length()>0)
                    ? outputTsvOverride
                    : Form("%s/04a_candidate_trees_y/tsv/YscolCandidates_%s_summary.tsv",
                           campaign.Data(), tag.Data());

  cout << "Rungroup tag      = " << tag << endl;
  cout << "Campaign dir      = " << campaign << endl;
  cout << "Metadata run key  = " << nrun << endl;
  cout << "Input replay ROOT = " << inputroot << endl;
  cout << "Ytar/delta cuts   = " << ytarCutFile << endl;
  cout << "YFP/YPFP cuts     = " << yCutFile << endl;
  cout << "Output ROOT       = " << outputroot << endl;
  cout << "Output TSV        = " << outputtsv << endl;

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

      if (vetoRungroup == tag.Data() && vetoAxis == "y") {
        preGmmVetoCuts.insert(vetoCutName);
      }
    }
  }

  cout << "Loaded " << preGmmVetoCuts.size()
       << " Y pre-GMM veto cut(s) for " << tag << endl;

  // ------------------------------------------------------------------
  // Load foil ytar/delta cuts.
  // ------------------------------------------------------------------
  vector<TCutG*> ytar_delta_cut;
  ytar_delta_cut.resize(NumFoil, nullptr);

  TFile *fYtarDeltaCut = nullptr;
  if (CutYtarFlag) {
    fYtarDeltaCut = TFile::Open(ytarCutFile, "READ");
    if (!fYtarDeltaCut || fYtarDeltaCut->IsZombie()) {
      cout << "ERROR: cannot open ytar/delta cut file " << ytarCutFile << endl;
      return;
    }

    for (Int_t nf=0; nf<NumFoil; nf++) {
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
  // Load YFP/YPFP yscol cuts: [foil][delta][yscol][component].
  //
  // Component cuts are preferred. A singleton compatibility alias is used
  // only when no component cuts exist for that physical yscol.
  // ------------------------------------------------------------------
  vector<vector<vector<vector<TCutG*> > > > ypfp_yfp_cut;
  ypfp_yfp_cut.resize(NumFoil);

  for (Int_t nf=0; nf<NumFoil; nf++) {
    ypfp_yfp_cut[nf].resize(nSlices);
    for (Int_t nd=0; nd<nSlices; nd++) {
      ypfp_yfp_cut[nf][nd].resize(9);
    }
  }

  TFile *fYCut = nullptr;
  if (CutYpFpYFpFlag) {
    fYCut = TFile::Open(yCutFile, "READ");

    if (!fYCut || fYCut->IsZombie()) {
      cout << "ERROR: cannot open YFP/YPFP cut file "
           << yCutFile << endl;
      return;
    }

    Int_t nLoadedComponents = 0;
    Int_t nLoadedSingletons = 0;

    // Load component cuts first.
    TIter nextKey(fYCut->GetListOfKeys());
    TKey* key = nullptr;

    while ((key = (TKey*)nextKey())) {
      TString name = key->GetName();

      Int_t ys=-1, nf=-1, nd=-1, part=-1, src=-1;

      Int_t ok = sscanf(
          name.Data(),
          "hYpFpYFp_cut_yscol_%d_nfoil_%d_ndel_%d_part%d_src_%d",
          &ys, &nf, &nd, &part, &src);

      if (ok != 5) continue;
      if (ys < 0 || ys > 8) continue;
      if (nf < 0 || nf >= NumFoil) continue;
      if (nd < 0 || nd >= nSlices) continue;

      if (preGmmVetoCuts.count(name.Data())) {
        cout << "VETO Y: " << name << endl;
        seenPreGmmVetoCuts.insert(name.Data());
        continue;
      }

      TObject* obj = key->ReadObj();
      if (!obj || !obj->InheritsFrom(TCutG::Class())) continue;

      TCutG* c = (TCutG*)obj->Clone(
          Form("ycand_component_y%d_f%d_d%d_p%d",
               ys, nf, nd, part));

      ypfp_yfp_cut[nf][nd][ys].push_back(c);
      nLoadedComponents++;
    }

    // Use singleton aliases only where no components were found.
    for (Int_t nf=0; nf<NumFoil; nf++) {
      for (Int_t nd=0; nd<nSlices; nd++) {
        for (Int_t ny=0; ny<9; ny++) {
          if (!ypfp_yfp_cut[nf][nd][ny].empty()) continue;

          TString cname = Form(
              "hYpFpYFp_cut_yscol_%d_nfoil_%d_ndel_%d",
              ny, nf, nd);

          if (preGmmVetoCuts.count(cname.Data())) {
            cout << "VETO Y: " << cname << endl;
            seenPreGmmVetoCuts.insert(cname.Data());
            continue;
          }

          TCutG* c = (TCutG*)fYCut->Get(cname);

          if (c) {
            ypfp_yfp_cut[nf][nd][ny].push_back(
                (TCutG*)c->Clone(
                    Form("ycand_singleton_y%d_f%d_d%d",
                         ny, nf, nd)));
            nLoadedSingletons++;
          }
        }
      }
    }

    cout << "Loaded " << nLoadedComponents
         << " YFP/YPFP yscol component cuts and "
         << nLoadedSingletons
         << " singleton fallback cuts." << endl;

    for (const string& cutName : preGmmVetoCuts) {
      if (!seenPreGmmVetoCuts.count(cutName)) {
        cout << "WARNING: configured Y pre-GMM veto cut not found: "
             << cutName << endl;
      }
    }
  }

  // ------------------------------------------------------------------
  // Open replay tree.
  // ------------------------------------------------------------------
  TFile *fin = TFile::Open(inputroot, "READ");
  if (!fin || fin->IsZombie()) {
    cout << "ERROR: cannot open input replay ROOT " << inputroot << endl;
    return;
  }

  TTree *T = (TTree*)fin->Get("T");
  if (!T) {
    cout << "ERROR: cannot find tree T in " << inputroot << endl;
    return;
  }

  // Input branches. Same branch names as make_fit_ntuple.C.
  Double_t sumnpe=0, etracknorm=0;
  Double_t ytar=0, xtar=0, reactx=0, reacty=0, reactz=0;
  Double_t delta=0, yptar=0, xptar=0;
  Double_t yfp=0, ypfp=0, xfp=0, xpfp=0;
  Double_t ysieve=0, xsieve=0;
  Double_t xbpm_tar=0, ybpm_tar=0, frx=0, fry=0;

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
  T->SetBranchAddress("H.rb.raster.fr_xa", &frx);
  T->SetBranchAddress("H.rb.raster.fr_ya", &fry);

  // ------------------------------------------------------------------
  // Output tree.
  // ------------------------------------------------------------------
  gSystem->mkdir(gSystem->DirName(outputroot), kTRUE);
  gSystem->mkdir(gSystem->DirName(outputtsv), kTRUE);
  TFile *fout = TFile::Open(outputroot, "RECREATE");
  if (!fout || fout->IsZombie()) {
    cout << "ERROR: cannot create output ROOT " << outputroot << endl;
    return;
  }

  TTree *out = new TTree("TYCand", "YFP/YPFP yscol candidate event tree for GMM cleanup");

  Int_t run_out=0, foil_out=-1, ndel_out=-1, yscol_out=-1;
  Int_t nYMatches_out=0;
  Long64_t entry_out=-1;
  Double_t zfoil_out=0;
  Double_t delta_low_out=0, delta_high_out=0;
  Int_t base_pass_out=0, ytar_cut_pass_out=0, ycut_pass_out=0;

  out->Branch("run", &run_out, "run/I");
  out->Branch("entry", &entry_out, "entry/L");
  out->Branch("foil", &foil_out, "foil/I");
  out->Branch("ndel", &ndel_out, "ndel/I");
  out->Branch("yscol", &yscol_out, "yscol/I");
  out->Branch("nYMatches", &nYMatches_out, "nYMatches/I");
  out->Branch("zfoil", &zfoil_out, "zfoil/D");
  out->Branch("delta_low", &delta_low_out, "delta_low/D");
  out->Branch("delta_high", &delta_high_out, "delta_high/D");
  out->Branch("base_pass", &base_pass_out, "base_pass/I");
  out->Branch("ytar_cut_pass", &ytar_cut_pass_out, "ytar_cut_pass/I");
  out->Branch("ycut_pass", &ycut_pass_out, "ycut_pass/I");

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

  // Overview histograms: one for each foil/delta, all yscols together.
  vector<vector<TH2D*> > hOverview;
  hOverview.resize(NumFoil);
  for (Int_t nf=0; nf<NumFoil; nf++) {
    hOverview[nf].resize(nSlices, nullptr);
    for (Int_t nd=0; nd<nSlices; nd++) {
      hOverview[nf][nd] = new TH2D(Form("overview_before_foil%d_ndel%d", nf, nd),
                                   Form("Run %d foil %d ndel %d;xsieve;ysieve", nrun, nf, nd),
                                   200, -15, 15, 200, -10, 10);
    }
  }

  vector<vector<vector<Long64_t> > > counts;
  counts.resize(NumFoil);
  for (Int_t nf=0; nf<NumFoil; nf++) {
    counts[nf].resize(nSlices);
    for (Int_t nd=0; nd<nSlices; nd++) counts[nf][nd].resize(9, 0);
  }

  Long64_t nentries = T->GetEntries();
  Long64_t nBase=0, nFoilDelta=0, nFilled=0, nOverlapEvents=0;
  cout << "Start candidate loop over " << nentries << " entries." << endl;

  for (Long64_t i=0; i<nentries; i++) {
    T->GetEntry(i);
    if (i%50000==0) cout << "  entry " << i << endl;

    // Match manual fit_ntuple baseline for now.
    base_pass_out = (etracknorm > 0.65 && sumnpe > 6.0 && delta > -10.0 && delta < 10.0) ? 1 : 0;
    if (!base_pass_out) continue;
    nBase++;

    Int_t nf_found=-1;
    for (Int_t nf=0; nf<NumFoil; nf++) {
      if (ytar_delta_cut[nf] && ytar_delta_cut[nf]->IsInside(ytar, delta)) {
        nf_found = nf;
      }
    }
    if (nf_found < 0) continue;

    Int_t nd_found=-1;
    for (Int_t nd=0; nd<nSlices; nd++) {
      if (delta >= delcut[nd] && delta < delcut[nd+1]) {
        nd_found = nd;
        break;
      }
    }
    if (nd_found < 0) continue;
    nFoilDelta++;

    vector<Int_t> yMatches;

    for (Int_t ny=0; ny<9; ny++) {
      bool matchedThisYscol = false;

      const vector<TCutG*>& components =
          ypfp_yfp_cut[nf_found][nd_found][ny];

      for (size_t ic=0; ic<components.size(); ic++) {
        TCutG* cy = components[ic];
        if (cy && cy->IsInside(ypfp, yfp)) {
          matchedThisYscol = true;
          break;
        }
      }

      if (matchedThisYscol) {
        yMatches.push_back(ny);
      }
    }

    if (yMatches.size() == 0) continue;
    if (yMatches.size() > 1) nOverlapEvents++;

    nYMatches_out = (Int_t)yMatches.size();

    for (size_t im=0; im<yMatches.size(); im++) {
      run_out = nrun;
      entry_out = i;
      foil_out = nf_found;
      ndel_out = nd_found;
      yscol_out = yMatches[im];
      zfoil_out = ztar_foil[nf_found];
      delta_low_out = delcut[nd_found];
      delta_high_out = delcut[nd_found+1];
      ytar_cut_pass_out = 1;
      ycut_pass_out = 1;

      out->Fill();
      hOverview[nf_found][nd_found]->Fill(xsieve, ysieve);
      counts[nf_found][nd_found][yscol_out]++;
      nFilled++;
    }
  }

  cout << "Done." << endl;
  cout << "  baseline-pass events        = " << nBase << endl;
  cout << "  foil+delta matched events   = " << nFoilDelta << endl;
  cout << "  candidate rows written      = " << nFilled << endl;
  cout << "  events with overlapping Y cuts = " << nOverlapEvents << endl;

  cout << "Counts by foil/ndel/yscol:" << endl;
  for (Int_t nf=0; nf<NumFoil; nf++) {
    for (Int_t nd=0; nd<nSlices; nd++) {
      cout << "  foil " << nf << " ndel " << nd << " :";
      for (Int_t ny=0; ny<9; ny++) cout << " y" << ny << "=" << counts[nf][nd][ny];
      cout << endl;
    }
  }

  fout->cd();
  out->Write();
  for (Int_t nf=0; nf<NumFoil; nf++) {
    for (Int_t nd=0; nd<nSlices; nd++) {
      hOverview[nf][nd]->Write();
    }
  }

  // Save a compact TSV summary in the campaign TSV directory.
  ofstream ofs(outputtsv.Data());
  ofs << "rungroup\trun\tfoil\tndel\tyscol\tN\tdelta_low\tdelta_high\tzfoil\n";
  for (Int_t nf=0; nf<NumFoil; nf++) {
    for (Int_t nd=0; nd<nSlices; nd++) {
      for (Int_t ny=0; ny<9; ny++) {
        ofs << tag << "\t" << nrun << "\t" << nf << "\t" << nd << "\t" << ny << "\t"
            << counts[nf][nd][ny] << "\t" << delcut[nd] << "\t" << delcut[nd+1]
            << "\t" << ztar_foil[nf] << "\n";
      }
    }
  }
  ofs.close();

  fout->Close();
  cout << "Wrote " << outputroot << endl;
  cout << "Wrote " << outputtsv << endl;
}

