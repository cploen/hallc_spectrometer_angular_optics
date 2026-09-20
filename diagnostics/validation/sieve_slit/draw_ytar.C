// Render persisted one-dimensional ROOT histograms; no event selection here.
#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TSystem.h>
#include <stdexcept>
#include <string>

void draw_ytar(const char* input, const char* stem, const char* context) {
  TFile file(input, "UPDATE");
  if (file.IsZombie()) throw std::runtime_error("Cannot open histogram file");
  gStyle->SetOptStat(1110);
  gStyle->SetStatX(.89); gStyle->SetStatY(.86);
  gStyle->SetStatW(.23); gStyle->SetStatH(.16);
  gStyle->SetTitleFontSize(.035);
  const char* names[] = {"dense_half", "shoulders", "out_of_core"};
  const int colors[] = {kBlue+1, kRed+1, kMagenta+1};
  for (int i=0; i<3; ++i) {
    auto* hist = dynamic_cast<TH1D*>(file.Get(names[i]));
    if (!hist) throw std::runtime_error("Missing TH1D");
    std::string canvasName = std::string("canvas_") + names[i];
    TCanvas canvas(canvasName.c_str(), hist->GetTitle(), 1200, 750);
    canvas.SetLeftMargin(.12); canvas.SetBottomMargin(.12); canvas.SetTopMargin(.15);
    hist->SetLineColor(colors[i]); hist->SetLineWidth(2);
    hist->GetXaxis()->SetTitle("y_{tar} (cm)");
    hist->GetYaxis()->SetTitle("Events / bin");
    hist->GetYaxis()->SetTitleOffset(1.35);
    hist->Draw("HIST");
    TLatex label; label.SetNDC(); label.SetTextAlign(22); label.SetTextSize(.026);
    label.DrawLatex(.5, .91, context);
    canvas.Update();
    canvas.SaveAs((std::string(stem)+"_"+names[i]+".png").c_str());
    canvas.Write(canvasName.c_str(), TObject::kOverwrite);
  }
  TNamed version("ROOT_version", gROOT->GetVersion());
  version.Write("ROOT_version", TObject::kOverwrite);
  file.Close();
}
