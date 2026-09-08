#include "../ytar_ridge_cut.C"
#include <stdexcept>
void require(bool ok) { if (!ok) throw std::runtime_error("ridge overlap regression failed"); }

RidgeCutResult fixtureRidge(int id, double center, vector<double> grid) {
  RidgeCutResult r;
  r.foilIndex=id; r.peak.ytar=center; r.deltaCenters=grid;
  for (double d : grid) {
    double c=center+0.02*d;
    r.rowNVec.push_back(100);
    r.ridgeCenter.push_back(c); r.ridgeCenterRaw.push_back(c);
    r.peakValVec.push_back(20);
    r.widthLeftRaw.push_back(2); r.widthRightRaw.push_back(2);
    r.widthLeftFinal.push_back(2); r.widthRightFinal.push_back(2);
    r.leftBoundary.push_back(c-2); r.rightBoundary.push_back(c+2);
  }
  return r;
}

double atDelta(const RidgeCutResult& r, const vector<double>& v, double d) {
  auto it=lower_bound(r.deltaCenters.begin(),r.deltaCenters.end(),d);
  size_t j=it-r.deltaCenters.begin();
  if (*it==d) return v[j];
  double f=(d-r.deltaCenters[j-1])/(r.deltaCenters[j]-r.deltaCenters[j-1]);
  return v[j-1]+f*(v[j]-v[j-1]);
}

void checkSeparation(const vector<RidgeCutResult>& rows) {
  for (size_t a=0;a<rows.size();++a) for(size_t b=a+1;b<rows.size();++b) {
    const auto& l=rows[a]; const auto& r=rows[b];
    double lo=max(l.deltaCenters.front(),r.deltaCenters.front());
    double hi=min(l.deltaCenters.back(),r.deltaCenters.back());
    for (double d=lo;d<=hi;d+=0.03125)
      require(atDelta(r,r.leftBoundary,d)-atDelta(l,l.rightBoundary,d)>=0.04-1e-10);
  }
}

void test_ridge_overlap() {
  // Reversed foil IDs model SHMS: geometric ordering must govern the guard.
  vector<RidgeCutResult> rows={fixtureRidge(0,2,{-2,-1,1,3,4}),
    fixtureRidge(2,-2,{-3,-1,0,2,4}),fixtureRidge(1,0,{-1,0,1})};
  require(ProtectRidgeEnvelopes(rows));
  require(rows[0].foilIndex==2 && rows[2].foilIndex==0);
  require(rows[0].deltaCenters.front()==-3 && rows[2].deltaCenters.back()==4);
  require(find(rows[0].rowNVec.begin(),rows[0].rowNVec.end(),-1)!=rows[0].rowNVec.end());
  checkSeparation(rows); // Includes outer pair beyond center foil's support.
  for(const auto& r:rows) for(size_t i=0;i<r.deltaCenters.size();++i)
    require(abs(r.widthLeftFinal[i]-(r.ridgeCenter[i]-r.leftBoundary[i]))<1e-10);
  vector<RidgeCutResult> same={fixtureRidge(0,-1,{-1,0,1}),fixtureRidge(1,1,{-1,0,1})};
  require(ProtectRidgeEnvelopes(same)); checkSeparation(same);
  vector<RidgeCutResult> disjoint={fixtureRidge(0,-1,{-3,-2}),fixtureRidge(1,1,{2,3})};
  require(ProtectRidgeEnvelopes(disjoint));
  require(abs(disjoint[0].rightBoundary[0]-0.94)<1e-10);
  vector<RidgeCutResult> crossing={fixtureRidge(0,0,{-1,0,1}),fixtureRidge(1,0.01,{-1,0,1})};
  require(!ProtectRidgeEnvelopes(crossing));
  cout<<"PASS: unequal grids, interpolation, support limits, geometric ordering, crossing rejection"<<endl;
}
