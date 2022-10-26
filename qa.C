using namespace ROOT;
using namespace ROOT::Math;
using namespace ROOT::RDF;
using fourVector = LorentzVector<PtEtaPhiE4D<double>>;

int Mproton(vector<fourVector> recMom, RVec<int> recSimIndex, RVec<int> simPdg, RVec<int> simMotherId)
{
  int mult=0;
  int nTracks=recMom.size();
  for (int i=0; i<nTracks; i++)
  {
    if (recMom.at(i).pt() < .2) continue;
    int simTrackId=recSimIndex.at(i);
    if (simTrackId < 0) continue;
    if (simPdg.at(simTrackId) != 2212) continue;
    if (simMotherId.at(simTrackId) > -1) continue;
    mult++;
  }
  return mult;
}

RVec<int> matchedPdg(RVec<int> recSimIndex, RVec<int> simPdg)
{
  RVec<int> matchedPdg_;
  for (auto& ind:recSimIndex)
  {
    if(ind>-1)
      matchedPdg_.push_back(simPdg.at(ind));
    else
      matchedPdg_.push_back(0);
  }
  return  matchedPdg_;
}

RVec<float> p(vector<fourVector> p)
{
  vector <float> p_;
  for (auto& mom:p) 
    p_.push_back(mom.P());
  return p_;
}

RVec<float> m2(RVec <fourVector> p, RVec<double> beta)
{
  vector<float> m2_;
  int n=p.size();
  for (int i=0;i<n;i++)
  {
    if (beta.at(i)>0)
      m2_.push_back(p.at(i).P()*p.at(i).P()*(1./beta.at(i)/beta.at(i)-1.));
    else
      m2_.push_back(-999);
  }
  return m2_;
};

int vsize(RVec<float> v) {return v.size();};

void qa(const char* in="*tree.root", const char* out="out.root") {
  RDataFrame d("t", in);
  cout << "Number of Events: " << *(d.Count()) << endl;
  vector <RResultPtr<::TH1D >> hists;
  vector <RResultPtr<::TH2D >> hists2d;
  TFile fOut(out,"recreate");
  auto dd=d
    .Filter("vtxChi2>0.")
//    .Define("M", vsize, {"trChi2"})
//    .Define("Mproton",Mproton,{"trMom","trSimIndex","simPdg", "simMotherId"});
    .Define("p", p, {"trMom"})
    .Define("matchedPdg", matchedPdg, {"trSimIndex","simPdg"})
    .Define("isProton", "matchedPdg==2212")
    .Define("isPion", "abs(matchedPdg)==211")
    .Define("isKaon", "abs(matchedPdg)==321")
    .Define("qp", "trCharge*p")
    .Define("m2_tof400",m2,{"trMom","trBetaTof400"})
    .Define("m2_tof700",m2,{"trMom","trBetaTof700"})
    .Define("isTof400", "trBetaTof400>0")
    .Define("isTof700", "trBetaTof700>0")
    .Define("vtxChi2ndf", "vtxChi2/vtxNdf")
//    .Define("deltaN",deltaN,{"m2_tof400","trBetaTof400"})
  ;
//    dd.Display("")->Print();
  dd.Foreach([](uint evtId){if (evtId % 100 == 0) cout << "\r" << evtId;}, {"evtId"}); // progress display
//  hists.push_back(dd.Histo1D({"hMproton","primary proton multiplicity;M_{tracks}",200,0,200}, "Mproton")); 
  hists2d.push_back(dd.Histo2D({"h2_vtx_XY","Reconstructed vertex XY;x (cm);y (cm)",500,-1,1,500,-1,1}, "vtxX", "vtxY")); 
  hists.push_back(dd.Histo1D({"h1_vtx_X","Reconstructed vertex X;x (cm)",500,-1,1}, "vtxX")); 
  hists.push_back(dd.Histo1D({"h1_vtx_Y","Reconstructed vertex Y;y (cm)",500,-1,1}, "vtxY")); 
  hists.push_back(dd.Histo1D({"h1_vtx_Z","Reconstructed vertex Z;z (cm)",500,-1,1}, "vtxZ")); 
  hists2d.push_back(dd.Histo2D({"h2_simVtx_XY","Simulated vertex XY;x (cm);y (cm)",500,-1,1,500,-1,1}, "simVtxX", "simVtxY")); 
  hists.push_back(dd.Histo1D({"h1_simVtx_X","Simulated vertex X;x (cm)",500,-1,1}, "simVtxX")); 
  hists.push_back(dd.Histo1D({"h1_simVtx_Y","Simulated vertex Y;y (cm)",500,-1,1}, "simVtxY")); 
  hists.push_back(dd.Histo1D({"h1_simVtx_Z","Simulated vertex Z;z (cm)",500,-1,1}, "simVtxZ")); 
  hists2d.push_back(dd.Histo2D({"h2_simVtx_X_vtx_X","Simulated vs reconstructed vertex X;x_{sim} (cm);x_{reco} (cm)",500,-1,1,500,-1,1}, "simVtxX", "vtxX")); 
  hists2d.push_back(dd.Histo2D({"h2_simVtx_Y_vtx_Y","Simulated vs reconstructed vertex Y;y_{sim} (cm);y_{reco} (cm)",500,-1,1,500,-1,1}, "simVtxY", "vtxY")); 
  hists2d.push_back(dd.Histo2D({"h2_simVtx_Z_vtx_Z","Simulated vs reconstructed vertex Z;z_{sim} (cm);z_{reco} (cm)",500,-1,1,500,-1,1}, "simVtxZ", "vtxZ")); 
  hists.push_back(dd.Histo1D({"h1_vtxChi2ndf","Reconstructed vertex Chi2;#chi^{2}",500,0,10}, "vtxChi2ndf")); 
//  hists2d.push_back(dd.Histo2D({"hMB","multiplicity vs impact parameter;b (fm);M_{tracks}",160,0,16,200,0,200}, "b", "M"));
  hists2d.push_back(dd.Histo2D({"h2_qp_m2_tof400","m^{2} vs pq (TOF400);p*q (GeV/c);m^{2} (GeV^{2}/c^{4})",1000,-10,10,500,-1,5}, "qp", "m2_tof400", "isTof400"));
  hists2d.push_back(dd.Histo2D({"h2_qp_m2_tof700","m^{2} vs pq (TOF700);p*q (GeV/c);m^{2} (GeV^{2}/c^{4})",1000,-10,10,500,-1,5}, "qp", "m2_tof700", "isTof700"));
  hists2d.push_back(dd.Histo2D({"h2_qp_m2_tof400_isTof700","m^{2} vs pq (TOF700+TOF400);p*q (GeV/c);m^{2} (GeV^{2}/c^{4})",1000,-10,10,500,-1,5}, "qp", "m2_tof400", "isTof700"));
  hists2d.push_back(dd.Histo2D({"h2_qp_m2_tof700_isTof400","m^{2} vs pq (TOF700+TOF400);p*q (GeV/c);m^{2} (GeV^{2}/c^{4})",1000,-10,10,500,-1,5}, "qp", "m2_tof700", "isTof400"));
  hists2d.push_back(dd.Histo2D({"h2_qp_m2_tof400_proton","Proton m^{2} vs pq (TOF400);p*q (GeV/c);m^{2} (GeV^{2}/c^{4})",1000,-10,10,500,-1,5}, "qp", "m2_tof400", "isProton"));
  hists2d.push_back(dd.Histo2D({"h2_qp_m2_tof400_pion","Pion m^{2} vs pq (TOF400);p*q (GeV/c);m^{2} (GeV^{2}/c^{4})",1000,-10,10,500,-1,5}, "qp", "m2_tof400", "isPion"));
  hists2d.push_back(dd.Histo2D({"h2_qp_m2_tof400_kaon","Kaon m^{2} vs pq (TOF400);p*q (GeV/c);m^{2} (GeV^{2}/c^{4})",1000,-10,10,500,-1,5}, "qp", "m2_tof400", "isKaon"));
  hists2d.push_back(dd.Histo2D({"h2_qp_m2_tof700_proton","Proton m^{2} vs pq (TOF700);p*q (GeV/c);m^{2} (GeV^{2}/c^{4})",1000,-10,10,500,-1,5}, "qp", "m2_tof700", "isProton"));
  hists2d.push_back(dd.Histo2D({"h2_qp_m2_tof700_pion","Pion m^{2} vs pq (TOF700);p*q (GeV/c);m^{2} (GeV^{2}/c^{4})",1000,-10,10,500,-1,5}, "qp", "m2_tof700", "isPion"));
  hists2d.push_back(dd.Histo2D({"h2_qp_m2_tof700_kaon","Kaon m^{2} vs pq (TOF700);p*q (GeV/c);m^{2} (GeV^{2}/c^{4})",1000,-10,10,500,-1,5}, "qp", "m2_tof700", "isKaon"));

  for (auto& hist:hists)
    hist->Write();
  for (auto& hist:hists2d)
    hist->Write();
  fOut.Close();
//  hists.back()->DrawClone();
//  hists2d.back()->DrawClone("colz");
}
