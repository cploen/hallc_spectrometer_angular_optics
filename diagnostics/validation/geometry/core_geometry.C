// Core membership is fixed upstream. No additional sieve-coordinate cut here.
#include "../../../spectrometer_config.h"
#include <TFile.h>
#include <TTree.h>
#include <TLeaf.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TGraph.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TPaveStats.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TError.h>
#include <TMatrixDSymEigen.h>
#include <array>
#include <fstream>
#include <iomanip>
#include <map>
#include <set>
#include <memory>
#include <stdexcept>

namespace core_geometry_detail {
struct Event {
  double reacty, residual, predicted, coreScore, ytar, ytarT;
  std::array<double,4> fp;
  int ndel;
};
double value(TTree* tree, const std::string& name) {
  auto leaf=tree->GetLeaf(name.c_str());
  if (!leaf) throw std::runtime_error("Missing branch: "+name);
  return leaf->GetValue();
}
void require(TTree* tree, const std::vector<std::string>& names) {
  for (const auto& name:names) value(tree,name);
}
// Mahalanobis distance in measured FP coordinates, independently in each delta slice.
// Singular directions are omitted. No angle, vertex or sieve position enters the distance.
std::vector<Event> tighter(const std::vector<Event>& events, double fraction) {
  std::map<int,std::vector<size_t>> slices;
  for(size_t i=0;i<events.size();++i) slices[events[i].ndel].push_back(i);
  std::vector<Event> result;
  for(const auto& slice:slices) {
    const auto& ids=slice.second;
    if(ids.size()<6) continue;
    std::array<double,4> mean{}, scale{};
    for(auto i:ids) for(int j=0;j<4;++j) mean[j]+=events[i].fp[j]/ids.size();
    for(auto i:ids) for(int j=0;j<4;++j) scale[j]+=std::pow(events[i].fp[j]-mean[j],2)/ids.size();
    for(auto& s:scale) s=s>0?std::sqrt(s):1.;
    TMatrixDSym covariance(4); covariance.Zero();
    for(auto i:ids) for(int j=0;j<4;++j) for(int k=0;k<4;++k)
      covariance(j,k)+=(events[i].fp[j]-mean[j])*(events[i].fp[k]-mean[k])/(scale[j]*scale[k]*ids.size());
    TMatrixDSymEigen eigen(covariance);
    auto vals=eigen.GetEigenValues(); auto vecs=eigen.GetEigenVectors();
    std::vector<std::pair<double,size_t>> distance;
    for(auto i:ids) {
      double d=0;
      for(int k=0;k<4;++k) if(vals[k]>1e-10) {
        double projection=0;
        for(int j=0;j<4;++j) projection+=vecs(j,k)*(events[i].fp[j]-mean[j])/scale[j];
        d+=projection*projection/vals[k];
      }
      distance.push_back({d,i});
    }
    std::sort(distance.begin(),distance.end());
    for(size_t j=0;j<size_t(std::floor(fraction*ids.size()));++j) result.push_back(events[distance[j].second]);
  }
  return result;
}
std::vector<Event> densestHalf(const std::vector<Event>& events) {
  std::map<int,std::vector<Event>> slices;
  for(const auto& e:events) slices[e.ndel].push_back(e);
  std::vector<Event> result;
  for(auto& item:slices) {
    auto& slice=item.second;
    std::sort(slice.begin(),slice.end(),[](const Event& a,const Event& b){return a.coreScore>b.coreScore;});
    const double threshold=slice[(slice.size()-1)/2].coreScore;
    size_t count=0;
    for(const auto& e:slice) if(e.coreScore>=threshold) {result.push_back(e);++count;}
    std::cout<<"Densest half: delta slice "<<item.first<<", core_score >= "<<threshold
             <<", retained "<<count<<" / "<<slice.size()<<" (ties included)\n";
  }
  return result;
}
void plotYtar(const std::vector<Event>& full,const std::vector<Event>& dense,
              const std::string& name,const std::string& output,TFile& file) {
  double low=full.front().ytar,high=low;
  for(const auto& e:full) {low=std::min({low,e.ytar,e.ytarT});high=std::max({high,e.ytar,e.ytarT});}
  double pad=std::max(.01,.05*(high-low));low-=pad;high+=pad;
  TH1D all((name+"_ytar_core").c_str(),"Full core;ytar (cm);Events",60,low,high);
  TH1D inner((name+"_ytar_dense_half").c_str(),"Densest half;ytar (cm);Events",60,low,high);
  TH1D allTarget((name+"_ytarT_core").c_str(),"",60,low,high);
  TH1D innerTarget((name+"_ytarT_dense_half").c_str(),"",60,low,high);
  for(const auto& e:full) {all.Fill(e.ytar);allTarget.Fill(e.ytarT);}
  for(const auto& e:dense) {inner.Fill(e.ytar);innerTarget.Fill(e.ytarT);}
  double maximum=1.45*std::max({all.GetMaximum(),inner.GetMaximum(),allTarget.GetMaximum(),innerTarget.GetMaximum()});
  gStyle->SetOptStat(1110);gStyle->SetOptFit(0);
  TCanvas canvas((name+"_ytar_canvas").c_str(),"",1200,600);canvas.Divide(2,1);
  TLegend legends[2]={TLegend(.17,.75,.43,.89),TLegend(.17,.75,.43,.89)};
  int panel=0;
  for(auto h:{&all,&inner}) {
    canvas.cd(++panel);gPad->SetLeftMargin(.14);gPad->SetBottomMargin(.14);
    h->SetMaximum(maximum);h->SetMinimum(0);h->SetLineColor(kBlue+2);h->SetLineWidth(2);
    h->Draw("HIST");gPad->Update();
    auto target=panel==1?&allTarget:&innerTarget;
    target->SetLineColor(kRed+1);target->SetLineWidth(2);target->SetStats(false);
    target->Draw("HIST SAME");
    auto& legend=legends[panel-1];legend.SetBorderSize(0);legend.SetTextSize(.035);
    legend.AddEntry(h,"ytar","l");legend.AddEntry(target,"ytarT","l");legend.Draw();
    auto stats=dynamic_cast<TPaveStats*>(h->GetListOfFunctions()->FindObject("stats"));
    if(stats){stats->SetX1NDC(.56);stats->SetX2NDC(.89);stats->SetY1NDC(.72);stats->SetY2NDC(.89);}
    gPad->Modified();gPad->Update();
  }
  canvas.SaveAs((output+"/"+name+"_ytar.pdf").c_str());
  canvas.SaveAs((output+"/"+name+"_ytar.png").c_str());
  file.cd();all.Write();inner.Write();allTarget.Write();innerTarget.Write();canvas.Write();
  std::ofstream summary(output+"/"+name+"_ytar.tsv");
  summary<<"selection\tn\tmean_cm\tstddev_cm\n"<<std::setprecision(10);
  summary<<"core\t"<<all.GetEntries()<<'\t'<<all.GetMean()<<'\t'<<all.GetStdDev()<<'\n';
  summary<<"dense_half\t"<<inner.GetEntries()<<'\t'<<inner.GetMean()<<'\t'<<inner.GetStdDev()<<'\n';
  summary<<"core_ytarT\t"<<allTarget.GetEntries()<<'\t'<<allTarget.GetMean()<<'\t'<<allTarget.GetStdDev()<<'\n';
  summary<<"dense_half_ytarT\t"<<innerTarget.GetEntries()<<'\t'<<innerTarget.GetMean()<<'\t'<<innerTarget.GetStdDev()<<'\n';
  gStyle->SetOptStat(10);gStyle->SetOptFit(111);
}
void plot(const std::vector<Event>& events, const std::string& name,
          const std::string& title, const std::string& output, TFile& file,
          std::ofstream& table, int xscol, int yscol, const std::string& subset, int ndel) {
  std::vector<double> x,y,predicted;
  double sx=0,sy=0,sxx=0,syy=0;
  for(const auto& e:events) {
    x.push_back(e.reacty); y.push_back(e.residual); predicted.push_back(e.predicted);
    sx+=x.back(); sy+=y.back(); sxx+=x.back()*x.back(); syy+=y.back()*y.back();
  }
  const double n=x.size();
  double variance=sxx/n-std::pow(sx/n,2);
  if(variance<=1e-16) {std::cout<<name<<": insufficient reacty range\n";return;}
  auto xr=std::minmax_element(x.begin(),x.end()), yr=std::minmax_element(y.begin(),y.end());
  double xmin=*xr.first,xmax=*xr.second,ymin=*yr.first,ymax=*yr.second;
  const double xpad=std::max(.001,.05*(xmax-xmin)), ypad=std::max(.01,.05*(ymax-ymin));
  TGraph graph(x.size(),x.data(),y.data());
  TF1 line((name+"_line").c_str(),"pol1",xmin,xmax);
  auto fit=graph.Fit(&line,"QSN"); // Unbinned ordinary least squares; no histogram-bin weighting.
  if(int(fit)!=0) throw std::runtime_error("Line fit failed: "+name);
  const double error=line.GetParError(1), slope=line.GetParameter(1);
  TH2D hist(name.c_str(),(title+";reacty (cm);1000 (xptar - xptarT) (mrad)").c_str(),
       70,xmin-xpad,xmax+xpad,80,ymin-ypad,ymax+ypad);
  for(size_t i=0;i<x.size();++i) hist.Fill(x[i],y[i]);
  TCanvas canvas((name+"_canvas").c_str(),"",1000,750);
  canvas.SetLeftMargin(.15);canvas.SetRightMargin(.14);canvas.SetBottomMargin(.13);
  canvas.SetTopMargin(.18);
  line.SetParNames("Intercept (mrad)","Slope (mrad/cm)");
  line.SetLineColor(kRed+1);line.SetLineWidth(2);
  hist.GetListOfFunctions()->Add(line.Clone());
  hist.Draw("COLZ");line.Draw("SAME");
  canvas.Update();
  auto stats=dynamic_cast<TPaveStats*>(hist.GetListOfFunctions()->FindObject("stats"));
  if(!stats) throw std::runtime_error("ROOT fit statistics box was not created");
  stats->SetX1NDC(.48);stats->SetX2NDC(.85);stats->SetY1NDC(.60);stats->SetY2NDC(.815);
  stats->SetTextFont(42);stats->SetTextSize(.024);

  double expected=0;
  for(auto p:predicted) expected+=p/n;
  TLatex subtitle;subtitle.SetNDC();subtitle.SetTextFont(42);subtitle.SetTextSize(.029);
  subtitle.DrawLatex(.15,.88,Form("Predicted slope from geometry: %.2f mrad/cm",expected));
  canvas.Modified();canvas.Update();
  canvas.SaveAs((output+"/"+name+".pdf").c_str());
  canvas.SaveAs((output+"/"+name+".png").c_str());
  file.cd(); hist.Write();line.Write();canvas.Write();
  std::cout<<name<<": N="<<n<<", reacty RMS="<<std::sqrt(variance)<<" cm, range=["<<xmin<<", "<<xmax
    <<"] cm, residual RMS="<<std::sqrt(std::max(0.,syy/n-std::pow(sy/n,2)))<<" mrad\n"
    <<"  slope="<<slope<<" +/- "<<error<<" mrad/cm; predicted="<<expected
    <<" mrad/cm; slope/error="<<(error>0?slope/error:0)
    <<"; (slope-predicted)/error="<<(error>0?(slope-expected)/error:0)<<"\n";
  table<<xscol<<'\t'<<yscol<<'\t'<<subset<<'\t'<<ndel<<'\t'<<n<<'\t'
    <<std::sqrt(variance)<<'\t'<<xmin<<'\t'<<xmax<<'\t'<<line.GetParameter(0)<<'\t'<<slope<<'\t'<<error<<'\t'
    <<std::sqrt(std::max(0.,syy/n-std::pow(sy/n,2)))<<'\t'<<expected<<'\n';
}
}

void core_geometry(const char* corePath,const char* replayPath,const char* output,const char* arm,
                   int opticsId,double angle,double zfoil,const char* holes,int minimum=30,
                   double fraction=.5,double foilWidth=2.,bool denseHalf=false,const char* afterburnerPath="") {
  using namespace core_geometry_detail;
  try {
    gROOT->SetBatch(true);gStyle->SetOptStat(10);gStyle->SetOptFit(111);
    gErrorIgnoreLevel=kWarning;std::cout<<std::unitbuf;
    gStyle->SetTitleFontSize(.035);
    auto spec=hallc::profileForName(arm);
    if(spec.name.empty()) throw std::runtime_error("Unknown spectrometer");
    std::map<Long64_t,std::pair<double,double>> afterburner;
    if(std::string(afterburnerPath).size()) {
      std::ifstream input(afterburnerPath);Long64_t entry;double xp,y;
      if(!input) throw std::runtime_error("Cannot open afterburner predictions");
      while(input>>entry>>xp>>y) afterburner[entry]={xp,y};
      if(afterburner.empty()) throw std::runtime_error("Empty afterburner predictions");
    }
    TFile core(corePath,"READ");
    std::unique_ptr<TFile> replay;
    if(!gSystem->AccessPathName(replayPath)) replay.reset(TFile::Open(replayPath,"READ"));
    if(core.IsZombie()||(replay&&replay->IsZombie())) throw std::runtime_error("Cannot open ROOT inputs");
    auto tree=dynamic_cast<TTree*>(core.Get("CoreSample"));
    auto source=replay?dynamic_cast<TTree*>(replay->Get("T")):tree;
    if(!tree||!source) throw std::runtime_error("Missing CoreSample or T tree");
    std::map<std::string,std::string> saved={{"react.x","reactx"},{"react.y","reacty"},{"react.z","reactz"},
      {"gtr.th","xptar"},{"gtr.dp","delta"},{"dc.x_fp","xfp"},{"dc.xp_fp","xpfp"},
      {"dc.y_fp","yfp"},{"dc.yp_fp","ypfp"},{"cal.etottracknorm","etracknorm"}};
    auto branch=[&](const std::string& suffix) {return replay?spec.branch(suffix):saved.at(suffix);};
    const std::string cerBranch=replay?spec.cherenkovBranch():"sumnpe";
    require(tree,{"entry","run","core_keep","zfoil","xscol","yscol","ndel","delta","delta_low","delta_high","xfp","xpfp","yfp","ypfp"});
    if(denseHalf) require(tree,{"core_score","ytar","xbpm_tar"});
    std::vector<std::string> suffix={"react.x","react.y","react.z","gtr.th","gtr.dp",
       "dc.x_fp","dc.xp_fp","dc.y_fp","dc.yp_fp","cal.etottracknorm"};
    if(replay) source->SetBranchStatus("*",0);
    for(auto& s:suffix) {s=branch(s);require(source,{s});source->SetBranchStatus(s.c_str(),1);}
    require(source,{cerBranch});source->SetBranchStatus(cerBranch.c_str(),1);
    std::vector<std::pair<int,int>> order;
    std::map<std::pair<int,int>,std::vector<Event>> data;
    std::istringstream tokens(holes);std::string token;
    while(std::getline(tokens,token,';')) {
      auto comma=token.find(',');int x=std::stoi(token.substr(0,comma)),y=std::stoi(token.substr(comma+1));
      if(x<0||x>=spec.nx||y<0||y>=spec.ny) throw std::runtime_error("Invalid hole");
      if(data.emplace(std::make_pair(x,y),std::vector<Event>{}).second) order.push_back({x,y});
    }
    std::set<Long64_t> seen;
    double maxClosure=0;size_t rejected=0;
    const double c=std::cos(angle*std::acos(-1.)/180.);
    for(Long64_t i=0;i<tree->GetEntries();++i) {
      tree->GetEntry(i);
      if(value(tree,"core_keep")!=1||std::abs(value(tree,"zfoil")-zfoil)>1e-8) continue;
      std::pair<int,int> hole{int(value(tree,"xscol")),int(value(tree,"yscol"))};
      if(!data.count(hole)) continue;
      auto entry=Long64_t(value(tree,"entry"));
      if(value(tree,"run")!=opticsId||entry<0||(replay&&entry>=source->GetEntries())||!seen.insert(entry).second)
        throw std::runtime_error("Invalid or duplicate core entry / optics ID");
      if(replay) source->GetEntry(entry);
      const std::array<std::string,4> fpNames={"xfp","xpfp","yfp","ypfp"};
      const std::array<std::string,4> fpSuffix={"dc.x_fp","dc.xp_fp","dc.y_fp","dc.yp_fp"};
      Event e{};e.ndel=int(value(tree,"ndel"));
      if(denseHalf) {
        e.coreScore=value(tree,"core_score");
        e.ytar=value(tree,"ytar");
        if(!std::isfinite(e.ytar)) throw std::runtime_error("Nonfinite ytar");
        if(!std::isfinite(e.coreScore)) throw std::runtime_error("Nonfinite core_score");
      }
      for(int j=0;j<4;++j) {
        e.fp[j]=value(source,branch(fpSuffix[j]));
        if(!std::isfinite(e.fp[j])||!std::isfinite(value(tree,fpNames[j]))||std::abs(e.fp[j]-value(tree,fpNames[j]))>1e-7)
          throw std::runtime_error("Core/replay FP mismatch; entry provenance failed");
      }
      double delta=value(source,branch("gtr.dp"));
      if(!std::isfinite(delta)||!std::isfinite(value(tree,"delta"))||std::abs(delta-value(tree,"delta"))>1e-7)
        throw std::runtime_error("Core/replay delta mismatch");
      if(!(delta>spec.deltaMin&&delta<spec.deltaMax&&delta>=value(tree,"delta_low")&&delta<value(tree,"delta_high")&&
           value(source,cerBranch)>2&&value(source,branch("cal.etottracknorm"))>.65&&
           std::abs(value(source,branch("react.z"))-zfoil)<foilWidth)) {++rejected;continue;}
      double reactx=value(source,branch("react.x")),xptar=value(source,branch("gtr.th"));
      if(!afterburner.empty()) {
        auto found=afterburner.find(entry);
        if(found==afterburner.end()) throw std::runtime_error("Missing afterburner event");
        xptar=found->second.first;e.ytar=found->second.second;
      }
      e.reacty=value(source,branch("react.y"));
      if(!std::isfinite(reactx)||!std::isfinite(e.reacty)||!std::isfinite(xptar))
        throw std::runtime_error("Nonfinite replay geometry");
      const double xsT=spec.xs(hole.first);
      // Use the saved event beam coordinate for the ytarT overlay.
      const auto truth=hallc::targetTruth(spec,angle,zfoil,xsT,spec.ys(hole.second),delta,reactx,e.reacty,denseHalf?value(tree,"xbpm_tar"):0.);
      e.ytarT=truth.ytar;
      if(denseHalf&&!std::isfinite(e.ytarT)) throw std::runtime_error("Nonfinite ytarT");
      e.residual=1000*(xptar-truth.xptar);
      e.predicted=spec.shms()?0.:1000/(spec.sieveDistance-zfoil*c);
      double expectedClosure=spec.shms()?0.:-e.reacty-spec.xMis(angle);
      maxClosure=std::max(maxClosure,std::abs(truth.xtar+spec.sieveDistance*truth.xptar-xsT-expectedClosure));
      data[hole].push_back(e);
    }
    gSystem->mkdir(output,true);
    TFile result((std::string(output)+"/geometry.root").c_str(),"RECREATE");
    std::ofstream table(std::string(output)+"/slopes.tsv");table<<std::setprecision(12);
    table<<"xscol\tyscol\tsubset\tndel\tn\treacty_rms_cm\treacty_min_cm\treacty_max_cm\tintercept_mrad\tslope_mrad_per_cm\tslope_error\tresidual_rms_mrad\tpredicted_slope\n";
    std::cout<<std::setprecision(7)
      <<(replay?"Core/replay identity verified in focal-plane coordinates and delta.":
                 "Using saved replay values from CoreSample; source replay comparison unavailable.")
      <<" Quality-cut exclusions="<<rejected
      <<"\nMaximum target-identity algebra discrepancy="<<maxClosure<<" cm\n"
      <<"Predictions assume reconstructed angles follow the reference ray. Core selection can bias slopes.\n"
      <<"Errors are ordinary least-squares statistical errors; no selection systematic is included.\n";
    if(spec.shms()) std::cout<<"SHMS targets already include reacty and xmis: predicted residual slope is zero.\n";
    int plots=0;
    for(auto hole:order) {
      auto& events=data[hole];
      std::string base="x"+std::to_string(hole.first)+"_y"+std::to_string(hole.second);
      if(events.size()<size_t(minimum)) {std::cout<<base<<": skipped, N="<<events.size()<<" < "<<minimum<<"\n";continue;}
      for(int subset=0;subset<2;++subset) {
        auto selected=subset?(denseHalf?densestHalf(events):tighter(events,fraction)):events;
        if(subset&&denseHalf) plotYtar(events,selected,base,output,result);
        std::string label=subset?(denseHalf?"dense_half":"fp_inner"):"core";
        // Pooled slopes and delta-slice slopes expose kinematic mixing.
        std::map<int,std::vector<Event>> slices;
        slices[-1]=selected;
        if(!subset) for(const auto& e:selected) slices[e.ndel].push_back(e);
        for(const auto& slice:slices) {
          if(slice.second.size()<size_t(minimum)) continue;
          std::string name=base+"_"+label+"_"+(slice.first<0?"all_delta":"delta"+std::to_string(slice.first));
          std::string group=gSystem->BaseName(output);group=group.substr(0,group.find('_'));
          std::ostringstream heading;
          if(!afterburner.empty()) heading<<"GMM afterburner: ";
          heading<<arm<<" "<<group<<", "<<angle<<" deg, foil "<<zfoil<<" cm, hole ("
            <<hole.first<<", "<<hole.second<<"), "<<(subset?(denseHalf?"densest half":"inner FP core"):"core");
          if(slice.first>=0) heading<<", delta slice "<<slice.first;
          std::string title=heading.str();
          plot(slice.second,name,title,output,result,table,hole.first,hole.second,label,slice.first);
          ++plots;
        }
      }
    }
    if(!plots) throw std::runtime_error("No sufficiently populated holes to plot");
    result.Close();
  } catch(const std::exception& error) {
    std::cerr<<"ERROR: "<<error.what()<<std::endl;gSystem->Exit(1);
  }
}
