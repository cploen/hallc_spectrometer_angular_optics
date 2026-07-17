#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TLine.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TString.h>
#include <limits.h>

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
using namespace std;

TString RealPathOrOriginal(const TString& path) {
  char resolved[PATH_MAX];

  if (realpath(path.Data(), resolved)) {
    return TString(resolved);
  }

  // If file doesn't exist or realpath fails, keep the original path.
  return path;
}

struct RunFoils {
  Int_t run;
  TString opticsID;
  vector<double> foils;
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

void DrawFitText(double x, double y, double znom, TF1* g) {
  if (!g) return;

  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.030);
  latex.SetTextAlign(13);

  latex.DrawLatex(
    x, y,
    Form("z_{nom}=%.1f cm: #mu=%.3f #pm %.3f cm, #sigma=%.3f cm, #chi^{2}/ndf=%.2f",
         znom,
         g->GetParameter(1),
         g->GetParError(1),
         fabs(g->GetParameter(2)),
         (g->GetNDF() > 0 ? g->GetChisquare() / g->GetNDF() : -1.0))
  );
}

void fit_ztar_peaks_from_replay(
  TString runListFile = "REPLACE_WITH_ZTAR_RUNLIST.txt",
  TString rootDir = "REPLACE_WITH_REPLAY_ROOTFILE_DIRECTORY/",
  TString label = "UNLABELED_ZTAR_TEST_DO_NOT_USE_FOR_RESULTS",
  Double_t deltaLow = -999.0,
  Double_t deltaHigh = 999.0
) {
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  if (label == "UNLABELED_ZTAR_TEST_DO_NOT_USE_FOR_RESULTS") {
    cout << "WARNING: using placeholder label: " << label << endl;
    cout << "         Pass a descriptive label, e.g. newfit_6p667_20260526." << endl;
  }
  if (runListFile == "REPLACE_WITH_ZTAR_RUNLIST.txt") {
    cout << "ERROR: runListFile is still a placeholder. Pass an explicit runlist." << endl;
    return;
  }
  if (rootDir == "REPLACE_WITH_REPLAY_ROOTFILE_DIRECTORY/") {
    cout << "ERROR: rootDir is still a placeholder. Pass an explicit replay ROOT directory." << endl;
    return;
  }
  if (deltaLow <= -900.0 || deltaHigh >= 900.0) {
    cout << "ERROR: delta range is still a placeholder. Pass explicit deltaLow and deltaHigh." << endl;
    return;
  }

  vector<RunFoils> runs = ReadRunList(runListFile);
  if (runs.size() == 0) {
    cout << "ERROR: no runs loaded from " << runListFile << endl;
    return;
  }

  TString tag = Form("%s_delta_%g_to_%g", label.Data(), deltaLow, deltaHigh);
  tag.ReplaceAll("-", "m");
  tag.ReplaceAll(".", "p");

  TString outTSV = Form("ztar_peak_fits_%s.tsv", tag.Data());
  TString outPDF = Form("plots/ztar/ztar_peak_fits_%s", tag.Data());

  gSystem->mkdir("plots", kTRUE);
  gSystem->mkdir("plots/ztar", kTRUE);

  ofstream out(outTSV.Data());
  out << "run\topticsID\tnominal_foil_z_cm"
      << "\tfit_mean_cm\tfit_mean_err_cm"
      << "\tfit_sigma_cm\tfit_sigma_err_cm"
      << "\tmean_minus_nominal_cm"
      << "\tchi2\tndf\tchi2_ndf"
      << "\tdelta_low\tdelta_high\tfit_low\tfit_high\tentries\treplay_file\n";

  bool firstPage = true;

  for (size_t irun = 0; irun < runs.size(); irun++) {
    RunFoils& rf = runs[irun];

    TString replayFile = Form("%s/nps_hms_optics_%s_1_-1.root",
                              rootDir.Data(), rf.opticsID.Data());
    TString replayFileReal = RealPathOrOriginal(replayFile);

    TFile* f = new TFile(replayFile);
    if (!f || f->IsZombie()) {
      cout << "ERROR: cannot open " << replayFile << endl;
      continue;
    }

    TTree* T = (TTree*)f->Get("T");
    if (!T) {
      cout << "ERROR: cannot find tree T in " << replayFile << endl;
      f->Close();
      continue;
    }

    TString hname = Form("hz_reactz_run%d", rf.run);
    TH1D* hz = new TH1D(
      hname,
      Form("Run %d: H.react.z, %.1f < #delta < %.1f; H.react.z [cm]; Counts",
           rf.run, deltaLow, deltaHigh),
      600, -20.0, 20.0
    );

    TString cut = Form("H.cer.npeSum>2 && H.gtr.dp>%g && H.gtr.dp<%g",
                       deltaLow, deltaHigh);

    TString drawCmd = Form("H.react.z >> %s", hname.Data());
    T->Draw(drawCmd, cut, "goff");

    if (hz->GetEntries() == 0) {
      cout << "WARNING: no entries for run " << rf.run
           << " with cut " << cut << endl;
      f->Close();
      continue;
    }

    TCanvas* c = new TCanvas(Form("c_run%d", rf.run),
                             Form("Run %d ztar peak fits", rf.run),
                             1000, 750);

    hz->SetLineWidth(2);
    hz->Draw("hist");

    vector<TF1*> fits;

    for (size_t ifoil = 0; ifoil < rf.foils.size(); ifoil++) {
      double znom = rf.foils[ifoil];

      double fitLow, fitHigh;
      GetFitWindow(znom, fitLow, fitHigh);


      Int_t b1 = hz->GetXaxis()->FindBin(fitLow);
      Int_t b2 = hz->GetXaxis()->FindBin(fitHigh);

      Int_t bmax = b1;
      double ymax = 0.0;

      for (Int_t b = b1; b <= b2; b++) {
        double y = hz->GetBinContent(b);
        if (y > ymax) {
          ymax = y;
          bmax = b;
        }
      }

      double peakGuess = hz->GetXaxis()->GetBinCenter(bmax);
      
      double coreHalfWidth = 0.75; // cm; try 0.8, 1.0, or 1.2
      double coreLow  = peakGuess - coreHalfWidth;
      double coreHigh = peakGuess + coreHalfWidth;
      
      TF1* g = new TF1(Form("g_run%d_foil%zu", rf.run, ifoil),
                       "gaus", coreLow, coreHigh);
      
      g->SetParameters(ymax, peakGuess, 0.5);
      hz->Fit(g, "RQ+");
      
      g->SetLineWidth(3);
      g->Draw("same");
      fits.push_back(g);

      double mean = g->GetParameter(1);
      double mean_err = g->GetParError(1);

      double sigma = fabs(g->GetParameter(2));
      double sigma_err = g->GetParError(2);

      double chi2 = g->GetChisquare();
      double ndf = g->GetNDF();
      double chi2_ndf = (ndf > 0.0) ? chi2 / ndf : -1.0;

      out << rf.run << "\t"
          << rf.opticsID << "\t"
          << znom << "\t"
          << mean << "\t"
          << mean_err << "\t"
          << sigma << "\t"
          << sigma_err << "\t"
          << mean - znom << "\t"
          << chi2 << "\t"
          << ndf << "\t"
          << chi2_ndf << "\t"
          << deltaLow << "\t"
          << deltaHigh << "\t"
          << coreLow << "\t"
          << coreHigh << "\t"
          << hz->GetEntries() << "\t"
          << replayFileReal << "\n";

      cout << "Run " << rf.run
           << " znom " << znom
           << " mean " << mean
           << " mean_err " << mean_err
           << " sigma " << sigma
           << " sigma_err " << sigma_err
           << " chi2/ndf " << chi2_ndf
           << " residual " << mean - znom
           << endl;

      TLine* l1 = new TLine(coreLow, 0, coreLow, ymax);
      TLine* l2 = new TLine(coreHigh, 0, coreHigh, ymax);
      l1->SetLineStyle(2);
      l2->SetLineStyle(2);
      l1->Draw("same");
      l2->Draw("same");
    }
     
     TLatex header;
     header.SetNDC();
     header.SetTextAlign(13);
     
     header.SetTextSize(0.018);
     header.DrawLatex(0.08, 0.91, Form("Source: %s", replayFileReal.Data()));
     
     header.SetTextSize(0.032);
     header.DrawLatex(0.13, 0.86, Form("Label: %s", label.Data()));
     header.DrawLatex(0.13, 0.82, Form("Cut: H.cer.npeSum>2, %.1f < H.gtr.dp < %.1f",
                                       deltaLow, deltaHigh));
     
    double yText = 0.75;
    for (size_t i = 0; i < rf.foils.size(); i++) {
      DrawFitText(0.13, yText, rf.foils[i], fits[i]);
      yText -= 0.045;
    }

    TString end = ".pdf";
    if (firstPage) {
      end = ".pdf(";
      firstPage = false;
    }
    if (irun == runs.size() - 1) {
      end = ".pdf)";
    }

    c->Print(outPDF + end);

    f->Close();
  }

  out.close();

  cout << endl;
  cout << "Wrote TSV: " << outTSV << endl;
  cout << "Wrote PDF: " << outPDF << endl;
}
