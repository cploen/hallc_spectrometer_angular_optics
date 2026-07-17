#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TString.h>
#include <TAxis.h>

#include <limits.h>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <set>
#include <string>
#include <cmath>
using namespace std;

TString RealPathOrOriginal(const TString& path) {
  char resolved[PATH_MAX];
  if (realpath(path.Data(), resolved)) return TString(resolved);
  return path;
}

struct RunFoils {
  Int_t run;
  TString opticsID;
  vector<double> foils;
};

struct DeltaSlice {
  double low;
  double high;
  TString label;
  int color;
};

vector<double> ParseFoilList(const string& s) {
  vector<double> foils;
  string item;
  stringstream ss(s);
  while (getline(ss, item, ',')) {
    if (item.size() == 0) continue;
    foils.push_back(atof(item.c_str()));
  }
  return foils;
}

vector<RunFoils> ReadRunList(const TString& runListFile) {
  vector<RunFoils> runs;
  ifstream in(runListFile.Data());

  if (!in.is_open()) {
    cout << "ERROR: cannot open run list: " << runListFile << endl;
    return runs;
  }

  string line;
  while (getline(in, line)) {
    if (line.size() == 0) continue;
    if (line[0] == '#') continue;

    stringstream ss(line);
    Int_t run;
    string opticsID;
    string foilString;

    ss >> run >> opticsID >> foilString;
    if (ss.fail()) {
      cout << "WARNING: could not parse line: " << line << endl;
      continue;
    }

    RunFoils rf;
    rf.run = run;
    rf.opticsID = opticsID.c_str();
    rf.foils = ParseFoilList(foilString);
    runs.push_back(rf);
  }

  return runs;
}

void GetFitWindow(double znom, double& fitLow, double& fitHigh) {
  if (fabs(znom) < 0.01) {
    fitLow = -2.0;
    fitHigh = 2.0;
  } else if (znom < -6.0) {
    fitLow = -10.0;
    fitHigh = -6.0;
  } else if (znom > 6.0) {
    fitLow = 6.0;
    fitHigh = 10.0;
  } else if (znom < 0.0) {
    fitLow = -5.0;
    fitHigh = -1.0;
  } else {
    fitLow = 1.0;
    fitHigh = 5.0;
  }
}

bool RunHasFoil(const RunFoils& rf, double znom) {
  for (double z : rf.foils) {
    if (fabs(z - znom) < 0.01) return true;
  }
  return false;
}

void fit_ztar_resolution_by_delta_slice(
  TString runListFile = "scripts/optics/ztar_stability/ztar_fit_runlist_6p667_combined.txt",
  TString rootDir = "ROOTfiles/OPTICS",
  TString label = "newfit_6p667_ztar_delta_resolution",
  bool includeHighPositiveDelta = false
) {
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  vector<RunFoils> runs = ReadRunList(runListFile);
  if (runs.empty()) {
    cout << "ERROR: no runs loaded." << endl;
    return;
  }

  vector<DeltaSlice> slices = {
    {-10.0, -8.0, "-10 < #delta < -8", kRed + 1},
    { -8.0, -5.0,  "-8 < #delta < -5", kBlue + 1},
    { -5.0,  0.0,   "-5 < #delta < 0", kGreen + 2},
    {  0.0,  5.0,    "0 < #delta < 5", kMagenta + 1},
    {  5.0, 10.0,    "5 < #delta < 10", kOrange + 7}
  };
  set<double> foilSet;
  for (auto& rf : runs) {
    for (double z : rf.foils) foilSet.insert(z);
  }

  TString outPDF = Form("plots/ztar/ztar_delta_slice_resolution_%s", label.Data());
  TString outPDF_individual = Form("plots/ztar/ztar_delta_slice_resolution_by_run_%s", label.Data());
  TString outTSV = Form("ztar_delta_slice_resolution_%s.tsv", label.Data());

  gSystem->mkdir("plots", kTRUE);
  gSystem->mkdir("plots/ztar", kTRUE);

  ofstream out(outTSV.Data());
  out << "nominal_foil_z_cm\tdelta_low\tdelta_high"
      << "\tfit_mean_cm\tfit_mean_err_cm"
      << "\tfit_sigma_cm\tfit_sigma_err_cm"
      << "\tmean_minus_nominal_cm"
      << "\tchi2\tndf\tchi2_ndf"
      << "\tentries\tcontributing_runs\n";

  bool firstPage = true;
  int pageIndex = 0;

  for (double znom : foilSet) {
    pageIndex++;

    double fitLow, fitHigh;
    GetFitWindow(znom, fitLow, fitHigh);

    TCanvas* c = new TCanvas(
      Form("c_foil_%g", znom),
      Form("Delta-slice ztar fits, foil %.1f cm", znom),
      1300, 850
    );

    c->SetRightMargin(0.34);
    c->SetLeftMargin(0.10);
    c->SetBottomMargin(0.12);

    TLegend* leg = new TLegend(0.68, 0.72, 0.96, 0.90);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.030);

    TLatex text;
    text.SetNDC();
    double fitTextSize = (slices.size() > 5) ? 0.023 : 0.027;
    double fitTextStep = (slices.size() > 5) ? 0.105 : 0.125;
    double fitLineStep = (slices.size() > 5) ? 0.030 : 0.035;
    double fitSigmaStep = (slices.size() > 5) ? 0.060 : 0.070;
    text.SetTextSize(fitTextSize);
    text.SetTextAlign(13);

    vector<TH1D*> hists;
    vector<TF1*> fits;

    double globalMax = 0.0;

    for (size_t is = 0; is < slices.size(); is++) {
      auto& sl = slices[is];

      TString hname = Form("hz_foil_%g_slice_%zu", znom, is);
      hname.ReplaceAll("-", "m");
      hname.ReplaceAll(".", "p");

      TH1D* hsum = new TH1D(
        hname,
        Form("Combined runs, foil z_{nom}=%.1f cm; H.react.z [cm]; Counts", znom),
        600, -20.0, 20.0
      );
      hsum->Sumw2();

      vector<int> usedRuns;

      for (auto& rf : runs) {
        if (!RunHasFoil(rf, znom)) continue;

        TString replayFile = Form("%s/nps_hms_optics_%s_1_-1.root",
                                  rootDir.Data(), rf.opticsID.Data());

        TFile* f = TFile::Open(replayFile, "READ");
        if (!f || f->IsZombie()) {
          cout << "WARNING: cannot open " << replayFile << endl;
          continue;
        }

        TTree* T = (TTree*)f->Get("T");
        if (!T) {
          cout << "WARNING: cannot find tree T in " << replayFile << endl;
          f->Close();
          continue;
        }

        TString htmpName = Form("htmp_run%d_foil_%g_slice_%zu", rf.run, znom, is);
        htmpName.ReplaceAll("-", "m");
        htmpName.ReplaceAll(".", "p");

        TH1D* htmp = new TH1D(htmpName, htmpName, 600, -20.0, 20.0);

        TString cut = Form("H.cer.npeSum>2 && H.cal.etottracknorm>0.65 && H.gtr.dp>%g && H.gtr.dp<%g",
                           sl.low, sl.high);

        TString drawCmd = Form("H.react.z >> %s", htmpName.Data());
        T->Draw(drawCmd, cut, "goff");

        if (htmp->GetEntries() > 0) {
          hsum->Add(htmp);
          usedRuns.push_back(rf.run);
        }

        delete htmp;
        f->Close();
      }

      hsum->SetLineColor(sl.color);
      hsum->SetLineWidth(2);

      double thisMax = hsum->GetMaximum();
      if (thisMax > globalMax) globalMax = thisMax;

      hists.push_back(hsum);

      if (hsum->GetEntries() <= 0) {
        fits.push_back(nullptr);
        continue;
      }

      Int_t b1 = hsum->GetXaxis()->FindBin(fitLow);
      Int_t b2 = hsum->GetXaxis()->FindBin(fitHigh);

      Int_t bmax = b1;
      double ymax = 0.0;

      for (Int_t b = b1; b <= b2; b++) {
        double y = hsum->GetBinContent(b);
        if (y > ymax) {
          ymax = y;
          bmax = b;
        }
      }

      double peakGuess = hsum->GetXaxis()->GetBinCenter(bmax);

      double coreHalfWidth = 0.75;
      double coreLow = peakGuess - coreHalfWidth;
      double coreHigh = peakGuess + coreHalfWidth;

      TF1* g = new TF1(
        Form("g_foil_%g_slice_%zu", znom, is),
        "gaus",
        coreLow,
        coreHigh
      );

      g->SetParameters(ymax, peakGuess, 0.5);
      hsum->Fit(g, "RQ0");

      g->SetLineColor(sl.color);
      g->SetLineWidth(3);
      fits.push_back(g);

      TString runList = "";
      for (size_t i = 0; i < usedRuns.size(); i++) {
        if (i > 0) runList += ",";
        runList += Form("%d", usedRuns[i]);
      }

      double mean = g->GetParameter(1);
      double meanErr = g->GetParError(1);
      double sigma = fabs(g->GetParameter(2));
      double sigmaErr = g->GetParError(2);
      double chi2 = g->GetChisquare();
      double ndf = g->GetNDF();
      double chi2ndf = (ndf > 0.0) ? chi2 / ndf : -1.0;

      out << znom << "\t"
          << sl.low << "\t" << sl.high << "\t"
          << mean << "\t" << meanErr << "\t"
          << sigma << "\t" << sigmaErr << "\t"
          << mean - znom << "\t"
          << chi2 << "\t" << ndf << "\t" << chi2ndf << "\t"
          << hsum->GetEntries() << "\t"
          << runList << "\n";
    }

    bool firstDraw = true;

    for (size_t is = 0; is < slices.size(); is++) {
      TH1D* h = hists[is];
      if (!h) continue;

      h->GetXaxis()->SetRangeUser(fitLow - 1.0, fitHigh + 1.0);
      h->SetMaximum(globalMax * 1.25);
      h->SetMinimum(0);

      if (firstDraw) {
        h->Draw("hist");
        firstDraw = false;
      } else {
        h->Draw("hist same");
      }

      if (fits[is]) fits[is]->Draw("same");

      leg->AddEntry(h, slices[is].label, "l");
    }

    leg->Draw();

    TLatex header;
    header.SetNDC();
    header.SetTextAlign(13);
    TString includedRuns = "";
    for (auto& rf_tmp : runs) {
      if (!RunHasFoil(rf_tmp, znom)) continue;
      if (includedRuns.Length() > 0) includedRuns += ", ";
      includedRuns += Form("%d", rf_tmp.run);
    }

    header.SetTextSize(0.026);
    header.DrawLatex(0.12, 0.875, Form("Combined like foil runs: %s", includedRuns.Data()));
    header.DrawLatex(0.12, 0.842, Form("cuts: H.cer.npeSum > 2, H.cal.etottracknorm > 0.65"));

    double yText = 0.69;

    for (size_t is = 0; is < slices.size(); is++) {
      if (!fits[is]) continue;

      TF1* g = fits[is];

      text.SetTextColor(slices[is].color);
      text.DrawLatex(
        0.68,
        yText,
        Form("%s", slices[is].label.Data())
      );

      text.DrawLatex(
        0.68,
        yText - fitLineStep,
        Form("#mu = %.3f #pm %.3f cm", g->GetParameter(1), g->GetParError(1))
      );

      text.DrawLatex(
        0.68,
        yText - fitSigmaStep,
        Form("#sigma = %.3f #pm %.3f cm", fabs(g->GetParameter(2)), g->GetParError(2))
      );

      yText -= fitTextStep;
    }

    text.SetTextColor(kBlack);

    TString end = ".pdf";
    if (firstPage) {
      end = ".pdf(";
      firstPage = false;
    }
    if (pageIndex == (int)foilSet.size()) {
      end = ".pdf)";
    }

    c->Print(outPDF + end);
    delete c;
  }


  // ------------------------------------------------------------------
  // Individual run x foil diagnostic pages.
  // These are not combined. They are for checking bad data by run.
  // ------------------------------------------------------------------
  bool firstIndividualPage = true;
  int individualPageIndex = 0;

  int nIndividualPages = 0;
  for (auto& rf : runs) nIndividualPages += rf.foils.size();

  for (auto& rf : runs) {
    TString replayFile = Form("%s/nps_hms_optics_%s_1_-1.root",
                              rootDir.Data(), rf.opticsID.Data());

    TFile* f = TFile::Open(replayFile, "READ");
    if (!f || f->IsZombie()) {
      cout << "WARNING: cannot open individual file " << replayFile << endl;
      continue;
    }

    TTree* T = (TTree*)f->Get("T");
    if (!T) {
      cout << "WARNING: cannot find tree T in " << replayFile << endl;
      f->Close();
      continue;
    }

    for (double znom : rf.foils) {
      individualPageIndex++;

      double fitLow, fitHigh;
      GetFitWindow(znom, fitLow, fitHigh);

      TCanvas* c = new TCanvas(
        Form("c_indiv_run%d_foil_%g", rf.run, znom),
        Form("Run %d, foil %.1f cm", rf.run, znom),
        1300, 850
      );

      c->SetRightMargin(0.34);
      c->SetLeftMargin(0.10);
      c->SetBottomMargin(0.12);

      TLegend* leg = new TLegend(0.68, 0.72, 0.96, 0.90);
      leg->SetBorderSize(0);
      leg->SetFillStyle(0);
      leg->SetTextSize(0.030);

      vector<TH1D*> hists;
      vector<TF1*> fits;
      double globalMax = 0.0;

      for (size_t is = 0; is < slices.size(); is++) {
        auto& sl = slices[is];

        TString hname = Form("hz_indiv_run%d_foil_%g_slice_%zu", rf.run, znom, is);
        hname.ReplaceAll("-", "m");
        hname.ReplaceAll(".", "p");

        TH1D* h = new TH1D(
          hname,
          Form("Run %d, foil z_{nom}=%.1f cm; H.react.z [cm]; Counts", rf.run, znom),
          600, -20.0, 20.0
        );

        TString cut = Form("H.cer.npeSum>2 && H.cal.etottracknorm>0.65 && H.gtr.dp>%g && H.gtr.dp<%g",
                           sl.low, sl.high);

        TString drawCmd = Form("H.react.z >> %s", hname.Data());
        T->Draw(drawCmd, cut, "goff");

        h->SetLineColor(sl.color);
        h->SetLineWidth(2);

        if (h->GetMaximum() > globalMax) globalMax = h->GetMaximum();

        hists.push_back(h);

        if (h->GetEntries() <= 0) {
          fits.push_back(nullptr);
          continue;
        }

        Int_t b1 = h->GetXaxis()->FindBin(fitLow);
        Int_t b2 = h->GetXaxis()->FindBin(fitHigh);

        Int_t bmax = b1;
        double ymax = 0.0;

        for (Int_t b = b1; b <= b2; b++) {
          double y = h->GetBinContent(b);
          if (y > ymax) {
            ymax = y;
            bmax = b;
          }
        }

        double peakGuess = h->GetXaxis()->GetBinCenter(bmax);
        double coreHalfWidth = 0.75;
        double coreLow = peakGuess - coreHalfWidth;
        double coreHigh = peakGuess + coreHalfWidth;

        TF1* g = new TF1(
          Form("g_indiv_run%d_foil_%g_slice_%zu", rf.run, znom, is),
          "gaus",
          coreLow,
          coreHigh
        );

        g->SetParameters(ymax, peakGuess, 0.5);
        h->Fit(g, "RQ0");

        g->SetLineColor(sl.color);
        g->SetLineWidth(3);
        fits.push_back(g);
      }

      bool firstDraw = true;

      for (size_t is = 0; is < slices.size(); is++) {
        TH1D* h = hists[is];
        if (!h) continue;

        h->GetXaxis()->SetRangeUser(fitLow - 1.0, fitHigh + 1.0);
        h->SetMaximum(globalMax * 1.25);
        h->SetMinimum(0);

        if (firstDraw) {
          h->Draw("hist");
          firstDraw = false;
        } else {
          h->Draw("hist same");
        }

        if (fits[is]) fits[is]->Draw("same");
        leg->AddEntry(h, slices[is].label, "l");
      }

      leg->Draw();

      TLatex header;
      header.SetNDC();
      header.SetTextAlign(13);
      // Canvas title already gives run and foil; do not repeat it inside the plot.
      header.SetTextSize(0.030);
      header.DrawLatex(0.12, 0.885, Form("Individual run diagnostic;"));
      header.DrawLatex(0.12, 0.852, Form("cuts: H.cer.npeSum > 2, H.cal.etottracknorm > 0.65"));

      TLatex text;
      text.SetNDC();
      double fitTextSize = (slices.size() > 5) ? 0.023 : 0.027;
      double fitTextStep = (slices.size() > 5) ? 0.105 : 0.125;
      double fitLineStep = (slices.size() > 5) ? 0.030 : 0.035;
      double fitSigmaStep = (slices.size() > 5) ? 0.060 : 0.070;
      text.SetTextSize(fitTextSize);
      text.SetTextAlign(13);

      double yText = 0.69;

      for (size_t is = 0; is < slices.size(); is++) {
        if (!fits[is]) continue;

        TF1* g = fits[is];

        text.SetTextColor(slices[is].color);
        text.DrawLatex(0.68, yText, Form("%s", slices[is].label.Data()));
        text.DrawLatex(
          0.68,
          yText - fitLineStep,
          Form("#mu = %.3f #pm %.3f cm", g->GetParameter(1), g->GetParError(1))
        );
        text.DrawLatex(
          0.68,
          yText - fitSigmaStep,
          Form("#sigma = %.3f #pm %.3f cm", fabs(g->GetParameter(2)), g->GetParError(2))
        );

        yText -= fitTextStep;
      }

      text.SetTextColor(kBlack);

      TString end = ".pdf";
      if (firstIndividualPage) {
        end = ".pdf(";
        firstIndividualPage = false;
      }
      if (individualPageIndex == nIndividualPages) {
        end = ".pdf)";
      }

      c->Print(outPDF_individual + end);
      delete c;
    }

    f->Close();
  }

  out.close();

  cout << endl;
  cout << "Wrote TSV: " << outTSV << endl;
  cout << "Wrote combined-by-foil PDF: " << outPDF << ".pdf" << endl;
  cout << "Wrote individual-run PDF: " << outPDF_individual << ".pdf" << endl;
}
