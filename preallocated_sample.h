#ifndef HALLC_PREALLOCATED_SAMPLE_H
#define HALLC_PREALLOCATED_SAMPLE_H
#include <map>
#include <set>
#include <tuple>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <cmath>
#include <TTree.h>
#include <TFile.h>

namespace hallc {
struct AllocatedLeaf {
  int optics, foil, ndel, x, y;
  double z, low, high;
  std::string group;
};
class PreallocatedSample {
 public:
  using ID=std::pair<std::string,Long64_t>;
  std::map<ID,AllocatedLeaf> expected;
  std::set<ID> used;
  std::map<std::tuple<std::string,int,int,int,int>,long long> counts;
  std::map<std::tuple<std::string,int,int,int,int>,AllocatedLeaf> countMetadata;
  std::ofstream admitted, rejected;
  std::string allocationHash;
  bool enabled=false;
  void load(const std::string& path, int maxEvents, const std::string& output) {
    enabled=!path.empty(); if (!enabled) return;
    std::ifstream in(path); std::string line;
    if (!std::getline(in,line) || line!="# core_preallocated_v1") throw std::runtime_error("Require a preallocated build manifest");
    if (!std::getline(in,line) || line.rfind("# allocation_sha256=",0)!=0 || line.size()!=84) throw std::runtime_error("Missing allocation hash");
    allocationHash=line.substr(20);
    const auto slash=path.find_last_of('/');
    std::ifstream build(path.substr(0,slash+1)+"build.json");
    std::string buildText((std::istreambuf_iterator<char>(build)),{});
    if (buildText.find(allocationHash)==std::string::npos || buildText.find("core_preallocated_v1")==std::string::npos)
      throw std::runtime_error("Require matching build.json beside membership manifest");
    std::getline(in,line);
    if (line!="rungroup\tentry\toptics_id\tfoil\tndel\tzfoil\tdelta_low\tdelta_high\txscol\tyscol") throw std::runtime_error("Invalid allocation columns");
    while(std::getline(in,line)) {
      if(line.empty()) continue;
      std::istringstream row(line); AllocatedLeaf leaf; Long64_t entry;
      if (!(row>>leaf.group>>entry>>leaf.optics>>leaf.foil>>leaf.ndel>>leaf.z>>leaf.low>>leaf.high>>leaf.x>>leaf.y)) throw std::runtime_error("Malformed membership row");
      if (!expected.emplace(ID(leaf.group,entry),leaf).second) throw std::runtime_error("Duplicate intended event ID");
    }
    if (expected.empty() || maxEvents<=0 || expected.size()>static_cast<size_t>(maxEvents)) throw std::runtime_error("Preallocated resource limit: zero/oversized sample; no truncation");
    admitted.open(output+"/solver_used_ids.tsv"); rejected.open(output+"/solver_rejected_ids.tsv");
    if (!admitted || !rejected) throw std::runtime_error("Cannot write solver membership audit");
    admitted<<"rungroup\tentry\tfoil\tndel\tzfoil\tdelta_low\tdelta_high\txscol\tyscol\n";
    rejected<<"rungroup\tentry\treason\n";
  }
  void reject(const ID& id,const std::string& why) { rejected<<id.first<<'\t'<<id.second<<'\t'<<why<<'\n'; rejected.flush(); }
  bool accept(const std::string& group,Long64_t entry,int run,int foil,int nd,int x,int y,double z,double delta,double ys,bool valid) {
    ID id(group,entry); auto it=expected.find(id);
    if(it==expected.end()){reject(id,"unexpected_id");return false;}
    const auto& e=it->second;
    if(!valid || run!=e.optics || foil!=e.foil || nd!=e.ndel || x!=e.x || y!=e.y || std::abs(z-e.z)>1e-8 || !(delta>=e.low && delta<e.high)) {reject(id,"invalid_values_or_leaf");return false;}
    if(!used.insert(id).second){reject(id,"duplicate_id");throw std::runtime_error("Duplicate supplied event");}
    admitted<<group<<'\t'<<entry<<'\t'<<foil<<'\t'<<nd<<'\t'<<e.z<<'\t'<<e.low<<'\t'<<e.high<<'\t'<<x<<'\t'<<y<<'\n';
    const auto key=std::make_tuple(group,foil,nd,x,y);
    counts[key]++; countMetadata[key]=e;
    return true;
  }
  void finish(const std::string& output) {
    if(!enabled)return;
    bool mismatch=false;
    for(const auto& kv:expected) if(!used.count(kv.first)){reject(kv.first,"not_admitted");mismatch=true;}
    admitted.flush();
    std::ofstream out(output+"/solver_used_counts.tsv"); out<<"rungroup\tfoil\tndel\txscol\tyscol\tzfoil\tdelta_low\tdelta_high\tselected\n";
    for(const auto& kv:counts) {const auto& e=countMetadata.at(kv.first);out<<std::get<0>(kv.first)<<'\t'<<std::get<1>(kv.first)<<'\t'<<std::get<2>(kv.first)<<'\t'<<std::get<3>(kv.first)<<'\t'<<std::get<4>(kv.first)<<'\t'<<e.z<<'\t'<<e.low<<'\t'<<e.high<<'\t'<<kv.second<<'\n';}
    if(mismatch) throw std::runtime_error("Failed preallocated handoff: intended/admitted membership mismatch");
  }
  void preflight(const std::string& dir,int fileID,const std::map<std::string,int>& settings,int nx,int ny) {
    std::set<ID> supplied;
    const char* doubles[]={"ys","ysT","xtar","xtarT","xptar","yptar","ytar","xptarT","yptarT","ytarT","ztarT","delta","xpfp","ypfp","xfp","yfp"};
    for(const auto& kv:settings) {
      bool needed=false; for(const auto& row:expected) if(row.first.first==kv.first) needed=true;
      std::string path=dir+"/Optics_"+std::to_string(kv.second)+"_"+std::to_string(fileID)+"_fit_tree_gmm.root";
      TFile file(path.c_str(),"READ");
      if(file.IsZombie()){if(needed)throw std::runtime_error("Missing required TFit input: "+path);else continue;}
      auto* tree=dynamic_cast<TTree*>(file.Get("TFit")); if(!tree)throw std::runtime_error("Missing TFit");
      for(const char* key:doubles) if(!tree->GetLeaf(key))throw std::runtime_error(std::string("Missing TFit field: ")+key);
      for(const char* key:{"entry","foil","ndel","xscol","yscol","sample","core_keep"}) if(!tree->GetLeaf(key))throw std::runtime_error(std::string("Missing preallocated field: ")+key);
      for(const char* key:{"foil","ndel","xscol","yscol","sample","core_keep"})
        if(std::string(tree->GetLeaf(key)->GetTypeName())!="Int_t") throw std::runtime_error("Preallocated labels and flags must be Int_t");
      if(std::string(tree->GetLeaf("entry")->GetTypeName())!="Long64_t")throw std::runtime_error("entry must be Long64_t");
      for(Long64_t i=0;i<tree->GetEntries();++i){
        tree->GetEntry(i); ID id(kv.first,tree->GetLeaf("entry")->GetValueLong64());
        auto it=expected.find(id);
        if(it==expected.end()){reject(id,"unexpected_id");throw std::runtime_error("Unexpected supplied ID");}
        if(!supplied.insert(id).second){reject(id,"duplicate_id");throw std::runtime_error("Duplicate supplied ID");}
        bool finite=true;for(const char* key:doubles)finite &= std::isfinite(tree->GetLeaf(key)->GetValue());
        auto val=[&](const char* k){return tree->GetLeaf(k)->GetValue();};
        const auto& e=it->second;
        bool valid=finite && val("sample")==1 && val("core_keep")==1 && val("xscol")==e.x && val("yscol")==e.y && val("foil")==e.foil && val("ndel")==e.ndel && e.optics==kv.second && e.x>=0 && e.x<nx && e.y>=0 && e.y<ny && std::abs(val("ztarT")-e.z)<1e-8 && val("delta")>=e.low && val("delta")<e.high && val("delta")>-15 && val("delta")<30;
        if(!valid){reject(id,"invalid_TFit");throw std::runtime_error("Invalid preallocated TFit event");}
      }
    }
    if(supplied.size()!=expected.size()) {
      for(const auto& e:expected)if(!supplied.count(e.first))reject(e.first,"missing_input");
      throw std::runtime_error("Missing supplied IDs or settings");
    }
  }
};
}
#endif
