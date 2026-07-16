#include <TSystem.h>
#include <TString.h>
#include "TFile.h"
#include "TTree.h"
#include <TNtuple.h>
#include "TCanvas.h"
#include <iostream>
#include <fstream>
#include "TMath.h"
#include "TH1F.h"
#include <TH2.h>
#include <TCutG.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TROOT.h>
#include <TMath.h>
#include <TLegend.h>
#include <TPaveLabel.h>
#include <TProfile.h>
#include <TPolyLine.h>
#include <TObjArray.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include<math.h>
#include <iostream>
#include <fstream>
using namespace std;

void plot_yptar_residuals(
  Int_t nrun,
  const char *inputTreeFile,
  const char *opticsMetadataFile,
  const char *outputDir,
  const char *tag
){
gStyle->SetPalette(1,0);
 gStyle->SetOptStat(1000011);
 gStyle->SetOptFit(11);
 gStyle->SetTitleOffset(1.,"Y");
 gStyle->SetTitleOffset(.7,"X");
 gStyle->SetLabelSize(0.04,"XY");
 gStyle->SetTitleSize(0.06,"XY");
 gStyle->SetPadLeftMargin(0.12);
  const Int_t Number=3;
  Double_t Red[Number] = { 1.0,0.0,0.0};
  Double_t Blue[Number] = { 1.0,0.0,1.0};
  Double_t Green[Number] = { 0.0,1.0,0.0};
 Double_t Len[Number] = { 0.0,.5,1.0};
 Int_t nb=50;
 TColor::CreateGradientColorTable(Number,Len,Red,Green,Blue,nb);
  //
 TObjArray HList(0);
 //  
 TString  tnrun=Form("%d",nrun);
  //  Get info for that optics run
 TString OpticsFile = opticsMetadataFile;
   ifstream file_optics(OpticsFile.Data());
 TString opticsline;
  TString OpticsID="";
  Double_t ymis =0.0;
  Int_t RunNum=0.;
  Double_t CentAngle=0.;
  Int_t SieveFlag=1;
  Int_t NumFoil=0;
  TString temp;
 //
  vector <Double_t> ztar_foil;
  Int_t ndelcut;
  vector<Double_t > delcut;
  if (file_optics.is_open()) {
    //
    cout << " Open file = " << OpticsFile << endl;
    while (RunNum != nrun && file_optics.good()) {
      temp.ReadToDelim(file_optics,',');
      //cout << temp << endl;
      if (temp.Atoi() == nrun) {
	RunNum = temp.Atoi();
      } else {
	temp.ReadLine(file_optics);
      }
    }
    if (RunNum==nrun) {
      temp.ReadToDelim(file_optics,',');
      OpticsID = temp;
      temp.ReadToDelim(file_optics,',');
      CentAngle = temp.Atof();
      temp.ReadToDelim(file_optics,',');
      NumFoil = temp.Atoi();
      temp.ReadToDelim(file_optics,',');
      SieveFlag = temp.Atoi();
      temp.ReadToDelim(file_optics,',');
      ndelcut = temp.Atoi();
      temp.ReadToDelim(file_optics);
      ymis = temp.Atof();
      for (Int_t nf=0;nf<NumFoil-1;nf++) {
        temp.ReadToDelim(file_optics,',');
	ztar_foil.push_back(temp.Atof());
      }
        temp.ReadToDelim(file_optics);
	ztar_foil.push_back(temp.Atof());
      for (Int_t nd=0;nd<ndelcut;nd++) {
        temp.ReadToDelim(file_optics,',');
	delcut.push_back(temp.Atof());
      }
        temp.ReadToDelim(file_optics);
	delcut.push_back(temp.Atof());
    }
  } else {
    cout << "ERROR: cannot open optics metadata: " << OpticsFile << endl;
    return;
  }
  if (RunNum != nrun) {
    cout << "ERROR: run " << nrun << " not found in optics metadata: "
         << OpticsFile << endl;
    return;
  }
  cout << RunNum << " " << OpticsID << " " << CentAngle << " " << NumFoil << " " << SieveFlag << endl;
  if (NumFoil==0) return;
 //
   TString inputroot = inputTreeFile;
   TString outdir = outputDir;
   TString outtag = tag;
   gSystem->ExpandPathName(inputroot);
   gSystem->ExpandPathName(outdir);
   gSystem->mkdir(outdir, kTRUE);

   TString outputpdf  = Form("%s/%s_yptar_residuals", outdir.Data(), outtag.Data());
   TString outputpdf1 = Form("%s/%s_yptar_residuals_by_sieve", outdir.Data(), outtag.Data());
   TString outputhist = Form("%s/%s_yptar_residuals.root", outdir.Data(), outtag.Data());

    TFile *hroot = TFile::Open(inputroot, "READ");
   if (!hroot || hroot->IsZombie()) {
     cout << "ERROR: cannot open fit tree file: " << inputroot << endl;
     return;
   }
   TTree *FitTree = (TTree*) hroot->Get("TFit");
   if (!FitTree) {
     cout << "ERROR: no TFit tree in: " << inputroot << endl;
     hroot->Close();
     return;
   }
   Double_t  ys,xtar,ztar,xptar,yptar,ytar,delta,xptarT,yptarT,ytarT,ztarT;
   Double_t xfp,xpfp,yfp,ypfp,ysieveT,ysieve;
   FitTree->SetBranchAddress("ys",&ysieve);
   FitTree->SetBranchAddress("ysT",&ysieveT);
   FitTree->SetBranchAddress("xtar",&xtar);
   FitTree->SetBranchAddress("ztar",&ztar);
   FitTree->SetBranchAddress("xptar",&xptar);
   FitTree->SetBranchAddress("yptar",&yptar);
   FitTree->SetBranchAddress("ytar",&ytar);
   FitTree->SetBranchAddress("xptarT",&xptarT);
   FitTree->SetBranchAddress("yptarT",&yptarT);
   FitTree->SetBranchAddress("ytarT",&ytarT);
   FitTree->SetBranchAddress("ztarT",&ztarT);
   FitTree->SetBranchAddress("delta",&delta);
   FitTree->SetBranchAddress("xpfp",&xpfp);
   FitTree->SetBranchAddress("ypfp",&ypfp);
   FitTree->SetBranchAddress("xfp",&xfp);
   FitTree->SetBranchAddress("yfp",&yfp);
 //
  vector <Double_t> ys_cent;
  vector <Double_t> xs_cent;
  for (Int_t nys=0;nys<9;nys++) {
    Double_t ypos=(nys-4)*0.6*2.54;
    ys_cent.push_back(ypos);
    Double_t xpos=(nys-4)*2.54;
    xs_cent.push_back(xpos);
  }
 //
 //
  vector<vector<vector<TH1F*> > > hYpDiff;
  vector<vector<TH2F*> > hYpDiff_YpTrue;
  vector<vector<TGraphErrors*> > gYpDiff_YpTrue;
  vector<TH1F*> hYtar;
  vector<TH1F*> hZtar;
  vector<vector<vector<Double_t > > > hYpDiff_mean;
  vector<vector<vector<Double_t > > > hYpDiff_sigma;
  vector<vector<vector<Double_t > > > hYpDiff_cent;
  vector<vector<vector<Double_t > > > hYpDiff_centerr;
  
  hYpDiff_YpTrue.resize(NumFoil);
  gYpDiff_YpTrue.resize(NumFoil);
  hYpDiff.resize(NumFoil);
  hYpDiff_mean.resize(NumFoil);
  hYpDiff_sigma.resize(NumFoil);
  hYpDiff_cent.resize(NumFoil);
  hYpDiff_centerr.resize(NumFoil);
  hYtar.resize(NumFoil);
  hZtar.resize(NumFoil);
   for (Int_t nf=0;nf<NumFoil;nf++) {
  hYpDiff_YpTrue[nf].resize(ndelcut);
  gYpDiff_YpTrue[nf].resize(ndelcut);
  hYpDiff[nf].resize(ndelcut);
  hYpDiff_mean[nf].resize(ndelcut);
  hYpDiff_sigma[nf].resize(ndelcut);
  hYpDiff_cent[nf].resize(ndelcut);
  hYpDiff_centerr[nf].resize(ndelcut);
   for (Int_t nd=0;nd<ndelcut;nd++) {
      hYpDiff[nf][nd].resize(9);
      hYpDiff_mean[nf][nd].resize(9);
      hYpDiff_sigma[nf][nd].resize(9);
      hYpDiff_cent[nf][nd].resize(9);
      hYpDiff_centerr[nf][nd].resize(9);
   }}
   //  
    for (Int_t nf=0;nf<NumFoil;nf++) {
      hYtar[nf] = new TH1F(Form("hYtar_%d",nf),Form("Run %s Ztar %3.2f ; Ytar",tnrun.Data(),ztar_foil[nf]),100,-5.,5.);
        HList.Add(hYtar[nf]);
     hZtar[nf] = new TH1F(Form("hZtar_%d",nf),Form("Run %s Ztar %3.2f ; Ztar",tnrun.Data(),ztar_foil[nf]),100,-15.,15.);
        HList.Add(hZtar[nf]);
    for (Int_t nd=0;nd<ndelcut;nd++) {
      Double_t DelCent=(delcut[nd+1]+delcut[nd])/2;
      hYpDiff_YpTrue[nf][nd]  = new TH2F(Form("hYpDiff_YpTrue_%d_DelCut_%d",nf,nd),Form("Run %s Ztar %3.2f DelCut %3.1f; YPtar_true (rad); Yptar_true - Yptar (mr) ",tnrun.Data(),ztar_foil[nf],DelCent),90,-.045,.045,80,-20.,20.);
        HList.Add(hYpDiff_YpTrue[nf][nd]);
    for (Int_t ns=0;ns<9;ns++) {
      hYpDiff[nf][nd][ns]  = new TH1F(Form("hYpDiff_%d_DelCut_%d_Ys_%d",nf,nd,ns),Form("Run %s Ztar %3.2f Ys = %d DelCut %3.1f; Yptar_true - Yptar (mr)",tnrun.Data(),ztar_foil[nf],ns,DelCent),80,-20.,20.);
        HList.Add(hYpDiff[nf][nd][ns]);
    }  }}
    //
    
 //
Long64_t nentries = FitTree->GetEntries();
 cout << " start loop " << nentries << endl;
	for (int i = 0; i < nentries; i++) {
      		FitTree->GetEntry(i);
    for (Int_t nf=0;nf<NumFoil;nf++) {
      if (abs(ztarT-ztar_foil[nf])<2) { 
       for (Int_t nd=0;nd<ndelcut;nd++) {
	 if ( delta >=delcut[nd] && delta <delcut[nd+1]) {
	hYtar[nf]->Fill(ytar);
	hZtar[nf]->Fill(ztar);
	   hYpDiff_YpTrue[nf][nd]->Fill(yptarT,(yptar-yptarT)*1000);
           for (Int_t ns=0;ns<9;ns++) {
	     if ( abs(ysieveT-ys_cent[ns])<.5) hYpDiff[nf][nd][ns]->Fill((yptar-yptarT)*1000);
	   }
	 }
       }
    }
    }		
    //	
	}
	//
	TCanvas* can2d[NumFoil];
	for  (Int_t nc=0;nc<NumFoil;nc++) {
	  can2d[nc] = new TCanvas(Form("Can2d_%d",nc),Form("Foil %d",nc), 700,700);
	  can2d[nc]->Divide(2,5);
	for  (Int_t nd=0;nd<ndelcut;nd++) {
	  can2d[nc]->cd(nd+1);
	  hYpDiff_YpTrue[nc][nd]->Draw("colz");
	}
	  TString end = ".pdf";
	  if (nc==0) end=".pdf(";
	  if (nc==NumFoil-1) end=".pdf)";
	  can2d[nc]->Print(outputpdf+end);

	}  //
	//
	//	TCanvas* candel[NumFoil][ndelcut];
	Double_t YP_offset=0.; // Offset to get delta = 0 ztar=0 to about zero.
	for  (Int_t nf=0;nf<NumFoil;nf++) {
	  Int_t colind = 0;
	for  (Int_t nd=0;nd<ndelcut;nd++) {
	  //	  candel[nf][nd] = new TCanvas(Form("Candel_%d_%d",nf,nd),Form("Candel_%d_%d",nf,nd), 700,700);
	  //	  candel[nf][nd]->Divide(2,5);
           for (Int_t ns=0;ns<9;ns++) {
	     //   candel[nf][nd]->cd(ns+1);
	     //  hYpDiff[nf][nd][ns]->Draw();
	     hYpDiff_mean[nf][nd][ns]= hYpDiff[nf][nd][ns]->GetMean()-YP_offset;
	     hYpDiff_sigma[nf][nd][ns]= hYpDiff[nf][nd][ns]->GetRMS();
	     hYpDiff_cent[nf][nd][ns]= ys_cent[ns];
	     hYpDiff_centerr[nf][nd][ns]= 0.0001;
	   }
	   gYpDiff_YpTrue[nf][nd] = new TGraphErrors(9,&hYpDiff_cent[nf][nd][0],&hYpDiff_mean[nf][nd][0],&hYpDiff_centerr[nf][nd][0],&hYpDiff_sigma[nf][nd][0]);
	   colind++;
	   if (colind==5) colind++;
	   if (colind==8) colind=1;
	   gYpDiff_YpTrue[nf][nd]->SetMarkerColor(colind);
	   gYpDiff_YpTrue[nf][nd]->SetMarkerStyle(nd+20);
	}}	  
 gStyle->SetPadRightMargin(0.2);
        
 //
	TCanvas* candel[NumFoil];
        TMultiGraph *mgr[NumFoil];
        TLegend *leg[NumFoil];
	for  (Int_t nf=0;nf<NumFoil;nf++) {
		  candel[nf] = new TCanvas(Form("Candel_%d",nf),Form("Candel_%d",nf), 700,700);
	  	  candel[nf]->Divide(1,1);
		  leg[nf] = new TLegend(.79,.65,.99,.95,"");
		 mgr[nf]=new TMultiGraph();  
	for  (Int_t nd=0;nd<ndelcut;nd++) {
	  mgr[nf]->Add(gYpDiff_YpTrue[nf][nd]);
      Double_t DelCent=(delcut[nd+1]+delcut[nd])/2;
      leg[nf]->AddEntry(gYpDiff_YpTrue[nf][nd],Form("Delta = %3.1f",DelCent),"p");
	    }
	candel[nf]->cd(1);
	mgr[nf]->SetTitle(Form("Ztar = %4.1f HMS Angle = %4.2f; Y_sieve (cm); YPtar -YP_true (mr)",ztar_foil[nf],CentAngle));
	mgr[nf]->SetMinimum(-20);
	mgr[nf]->SetMaximum(+20);
	mgr[nf]->Draw("AP");
	leg[nf]->Draw();
	  TString end = ".pdf";
	  if (nf==0) end=".pdf(";
	  if (nf==NumFoil-1) end=".pdf)";
	  candel[nf]->Print(outputpdf1+end);
	}
	//
 TFile hsimc(outputhist,"recreate");
	HList.Write();

//
}
