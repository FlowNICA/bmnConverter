using namespace ROOT;
using namespace ROOT::Math;
using namespace ROOT::RDF;
using fourVector=LorentzVector<PtEtaPhiE4D<double>>;

BmnFieldMap* magField{nullptr};

TChain* makeChain(string& filename, const char* treename) {
  cout << "Adding files to chain:" << endl;
  TChain *chain = new TChain(treename);
  if (filename.rfind(".root") < filename.size())
    chain->Add(filename.data());
  else {
    TFileCollection fc("fc", "", filename.c_str());
    chain->AddFileInfoList((TCollection*)fc.GetList());
  }
  chain->ls();
  return chain;
}

bool isGoodSimTrack(const CbmMCTrack &track)
{
  if (track.GetStartZ() > 900) return false;
  return true;
}

XYZVector ExtrapolateStraightLine(const FairTrackParam *par, float z)
{
  float dz = z - par->GetZ();
  float x = par->GetX() + par->GetTx() * dz;
  float y = par->GetY() + par->GetTy() * dz;
  return {x, y, z};
}

double getMass(int pdg)
{
  if(pdg<1000000) 
    return TDatabasePDG::Instance()->GetParticle(pdg)->Mass();
  else 
    return 0.931*(pdg/10%1000);
}

int getCharge(int pdg)
{
  if(pdg<1000000)
    return TDatabasePDG::Instance()->GetParticle(pdg)->Charge();
  else 
    return pdg/10000%1000;
}

vector<fourVector> simMomentum(const RVec<CbmMCTrack> tracks)
{
  vector<fourVector> momenta;
  for (auto& track:tracks) {
    if (!isGoodSimTrack(track))
      continue;
    TVector3 mom;
    track.GetMomentum(mom);
    double mass=getMass(track.GetPdgCode());
    momenta.push_back({mom.Pt(),mom.Eta(),mom.Phi(),sqrt(mass*mass+mom.Mag2())});
  }
  return momenta;
}
 
vector<XYZTVector> simPosStart(const RVec<CbmMCTrack> tracks)
{
  vector<XYZTVector> pos;
  for (auto& track:tracks) {
    if (!isGoodSimTrack(track))
      continue;
    pos.push_back({track.GetStartX(),track.GetStartY(),track.GetStartZ(),track.GetStartT()});
  }
  return pos;
}

RVec<int> simMotherId(const RVec<CbmMCTrack> tracks)
{
  vector<int> mothId;
  for (auto& track:tracks) {
    if (!isGoodSimTrack(track))
      continue;
    mothId.push_back(track.GetMotherId());
  }
  return mothId;
}

RVec<int> simPdg(const RVec<CbmMCTrack> tracks)
{
  vector<int> pdg;
  for (auto& track:tracks) {
    if (!isGoodSimTrack(track))
      continue;
    pdg.push_back(track.GetPdgCode());
  }
  return pdg;
}

RVec<short> simCharge (const RVec<int> pdg)
{
  vector<short> ch;
  for (auto &p:pdg)
    ch.push_back(getCharge(p));
  return ch;
}

vector<float> trackP(const RVec<BmnGlobalTrack> tracks)
{
  vector<float> momenta;
  for (auto track:tracks)
    momenta.push_back(1./track.GetParamFirst()->GetQp());
  return momenta;
}

vector<fourVector> trackMomentum(const RVec<BmnGlobalTrack> tracks)
{
  vector<fourVector> momenta;
  for (auto track:tracks) {
    auto *par = track.GetParamFirst();   
    TVector3 mom;
    par->Momentum(mom);
    momenta.push_back({mom.Pt(),mom.Eta(),mom.Phi(),0});
  }   
  return momenta;
}

RVec<short> recCharge(const RVec<BmnGlobalTrack> tracks)
{
  vector<short> charge;
  for (auto track:tracks) {
    int q = track.GetParamFirst()->GetQp() > 0 ? 1 : -1;
    charge.push_back(q);
  }
  return charge;
}

vector<XYZVector> recDca(const RVec<BmnGlobalTrack> tracks, const CbmVertex vtx)
{
  vector<XYZVector> dca;
  for (auto track:tracks) {
    auto par = track.GetParamFirst();
    dca.push_back({par->GetX()-vtx.GetX(),par->GetY()-vtx.GetY(),par->GetZ()-vtx.GetZ()});
  }
  return dca;
}

vector< vector<float> > covMatrix(RVec<BmnGlobalTrack> global_tracks, RVec<CbmStsTrack> tracks)
{
  vector<vector<float>> covariance_matrix;
  for (auto& global_track : global_tracks) {
    auto idx = global_track.GetGemTrackIndex();
    auto track = tracks.at(idx);
    auto* par = track.GetParamFirst();
    covariance_matrix.emplace_back();
    for( int i=0; i<5; ++i ){
      for( int j=0; j<=i; ++j ){
        covariance_matrix.back().push_back( par->GetCovariance(i, j) );
      }
    }
    // Lower triangle of the symmetric covariance matrix
    // C[x, y, tx, ty, Qp]
    // { c_00, c1[0..1], c2[0..2], ... c4[0..4] }
  }
  return covariance_matrix;
}

float determinant3x3( const std::array<std::array<float, 3>, 3>& matrix ){
  auto x_0 = matrix[0][0] * ( matrix[1][1]*matrix[2][2] - matrix[1][2]*matrix[2][1]  );
  auto x_1 = matrix[0][1] * ( matrix[1][0]*matrix[2][2] - matrix[1][2]*matrix[2][0]  );
  auto x_2 = matrix[0][2] * ( matrix[1][0]*matrix[2][1] - matrix[1][1]*matrix[2][0]  );

  return x_0 - x_1 + x_2;
}

std::array<float, 3> cramerFieldSolver3x3( std::array<float, 3> field, std::array<float, 3> coordinate ){
  // Solving the system of equation to extract parameters of quadratic extrapolation of the magnetic field
  // Ax = B
  // xi = detAi / detA
  std::array<std::array<float, 3>, 3> A;
  A[0] = {1.0f, 1.0f, 1.0f };
  A[1] = { coordinate[0], coordinate[1], coordinate[2] };
  A[2] = { coordinate[0]*coordinate[0], coordinate[1]*coordinate[1], coordinate[2]*coordinate[2] };

  auto A0 = A;
  A0[0] = field;
  auto A1 = A;
  A1[1] = field;
  auto A2 = A;
  A2[2] = field;

  auto detA = determinant3x3( A );
  auto detA0 = determinant3x3( A0 );
  auto detA1 = determinant3x3( A1 );
  auto detA2 = determinant3x3( A2 );

  auto p0 = detA0 / detA;
  auto p1 = detA1 / detA;
  auto p2 = detA2 / detA;

  return {p0, p1, p2};
}

vector< vector<float> > magneticField(RVec<BmnGlobalTrack> global_tracks, RVec<CbmStsTrack> tracks, RVec<CbmStsHit> sts_hits)
{
  vector<vector<float>> magnetic_field;
  for (auto& global_track : global_tracks) {
    auto idx = global_track.GetGemTrackIndex();
    auto track = tracks.at(idx);
    std::array<float, 3> hit_z;
    std::array<float, 3> hit_bx;
    std::array<float, 3> hit_by;
    std::array<float, 3> hit_bz;

    for( int i=0; i<3; ++i ){
      // It seems size of the hitmap cannot be less than 4, but just to be safe
      if( i > track.GetStsHits()->GetSize() )
        magnetic_field.push_back( std::vector<float>(10, 0.0f) );

      auto sts_idx = track.GetStsHits()->At(i);
      auto x = sts_hits.at(sts_idx).GetX();
      auto y = sts_hits.at(sts_idx).GetY();
      auto z = sts_hits.at(sts_idx).GetZ();

      hit_z.at(i) = z;
      hit_bx.at(i) = magField->GetBx( x, y, z ); // kGs
      hit_by.at(i) = magField->GetBy( x, y, z ); // kGs
      hit_bz.at(i) = magField->GetBz( x, y, z ); // kGs
    }

    auto parameters_bx = cramerFieldSolver3x3( hit_bx, hit_z );
    auto parameters_by = cramerFieldSolver3x3( hit_by, hit_z );
    auto parameters_bz = cramerFieldSolver3x3( hit_bz, hit_z );

    magnetic_field.emplace_back();
    for( const auto& c : parameters_bx )
      magnetic_field.back().push_back( c );
    for( const auto& c : parameters_by )
      magnetic_field.back().push_back( c );
    for( const auto& c : parameters_bz )
      magnetic_field.back().push_back( c );
    magnetic_field.back().push_back( 0.0 ); // z0
  }
  return magnetic_field;
}

vector<vector<float>> stsTrackParameters(RVec<BmnGlobalTrack> global_tracks, RVec<CbmStsTrack> tracks)
{
  vector<vector<float>> parameters;
  for (auto& global_track : global_tracks) {
    auto idx = global_track.GetGemTrackIndex();
    auto track = tracks.at(idx);
    auto* par = track.GetParamFirst();
    parameters.emplace_back();
    parameters.back().push_back( par->GetX() );
    parameters.back().push_back( par->GetY() );
    parameters.back().push_back( par->GetZ() );
    parameters.back().push_back( par->GetTx() );
    parameters.back().push_back( par->GetTy() );
    parameters.back().push_back( par->GetQp() );
  }
  return parameters;
}

vector<fourVector> stsTrackMomentum(RVec<CbmStsTrack> tracks)
{
  vector<fourVector> momenta;
  for (auto& track : tracks) {
    auto *par = track.GetParamFirst();
    TVector3 mom;
    par->Momentum(mom);
    momenta.push_back({mom.Pt(),mom.Eta(),mom.Phi(),0});
  }
  return momenta;
}

RVec<int> recSimIndex(const RVec<BmnGlobalTrack> recTracks, const RVec<CbmMCTrack> simTracks)
{
  vector<int> newIndex;
  int shift=0;
  int nSimTracks = simTracks.size();
  for (int i=0;i<nSimTracks;i++) {
    if (!isGoodSimTrack(simTracks.at(i)))
    {
      shift++;
      newIndex.push_back(-1);
    }
    else
      newIndex.push_back(i-shift);
  }
  vector<int> simIndex;
  for (auto track:recTracks) {
    int oldIndex=track.GetRefIndex();
    if (oldIndex<0 || oldIndex>=nSimTracks)
      simIndex.push_back(-1);
    else
      simIndex.push_back(newIndex.at(oldIndex));
  }
  return simIndex;
}

vector<XYZVector> recPosLast(const RVec<BmnGlobalTrack> tracks)
{
  vector<XYZVector> pos;
  for (auto track:tracks){
    auto par=track.GetParamLast();
    pos.push_back({par->GetX(), par->GetY(), par->GetZ()});
  }
  return pos;
}

vector<XYZVector> recPos450(const RVec<BmnGlobalTrack> tracks)
{
  vector<XYZVector> pos;
  for (auto track:tracks) 
    pos.push_back(ExtrapolateStraightLine(track.GetParamLast(), 450));
  return pos;
}

vector<XYZVector> tofHitPosition(const TClonesArray hits)
{
  vector<XYZVector> pos;
  for (const auto& hitObj:hits){
    auto hit=(BmnTofHit*)hitObj;
    pos.push_back({hit->GetX(),hit->GetY(),hit->GetZ()});
  }
  return pos;
}

vector<XYZVector> tofRes(const RVec<BmnGlobalTrack> tracks, const TClonesArray hits)
{
  vector<XYZVector> res;
  auto testHit=(BmnTofHit*)hits.At(0);
  if (testHit)
  {
    bool tof400 = testHit->GetZ()<600 ? true : false;
    for (auto track:tracks){
      int hitIndex = tof400 ? track.GetTof1HitIndex() : track.GetTof2HitIndex();
      if (hitIndex<0 || !hits.At(hitIndex)) 
      {
        res.push_back({-999,-999,-999});
        continue;
      }
      auto hit=(BmnTofHit*)hits.At(hitIndex);
      auto par=track.GetParamLast();
      TVector3 pos;
      par->Position(pos);
      auto posAtHitZ_ = ExtrapolateStraightLine(par, hit->GetZ());
      TVector3 posAtHitZ(posAtHitZ_.x(), posAtHitZ_.y(), posAtHitZ_.z());
      TVector3 hitPos;
      hit->Position(hitPos);
      auto posDiff = pos - posAtHitZ;
      res.push_back({posDiff.X(), posDiff.Y(), posDiff.Z()});
    }
  }
  return res;
}

RVec<int> moduleId (const vector<XYZVector> modulePos)
{
  vector <int> moduleIds;
  for (int i=0;i<modulePos.size();i++)
    moduleIds.push_back(i+1);
  return moduleIds;
}

vector<XYZVector> modulePos (const char *geoFile, const char *detectorTag)
{
  bool verbose=false;
  map <int,XYZVector> modulePosMap;
  printf("Reading %s geometry from geometry file\n", detectorTag);
  TGeoManager* geoMan = TGeoManager::Import(geoFile, "FAIRGeom");
  if( !geoMan )
    throw runtime_error(Form("ERROR: No TGeoManager in file %s", geoFile));
  TGeoNode* caveNode = geoMan->GetTopNode();
  if( !caveNode )
    throw runtime_error(Form("ERROR: No cave node in file %s", geoFile));
  TGeoNode* detectorNode = nullptr;
  TString nodeName;
  
  bool nodeFound=false;
  for (int i = 0; i < caveNode->GetNdaughters(); i++) {
    detectorNode = caveNode->GetDaughter(i);
    nodeName = detectorNode->GetName();
    nodeName.ToLower();
    if (nodeName.Contains(detectorTag))
    {	
      nodeFound=true;
      break;
    }
  }
  if( !nodeFound )
    throw runtime_error(Form("ERROR: No detector node %s in cave", detectorTag));
  detectorNode = detectorNode->GetDaughter(0);

  auto geoMatrix = detectorNode->GetMatrix();
  auto geoBox = (TGeoBBox*) detectorNode->GetVolume()->GetShape();
  TVector3 frontFaceLocal(0, 0, -geoBox->GetDZ());
  TVector3 frontFaceGlobal;
  geoMatrix->LocalToMaster(&frontFaceLocal[0], &frontFaceGlobal[0]);

  nodeName=detectorNode->GetName();
  if (nodeName.Contains("box"))
    detectorNode = detectorNode->GetDaughter(detectorNode->GetNdaughters()-1);
  printf("%s node name: %s\n", detectorTag, detectorNode->GetName());

  int nModules = detectorNode->GetNdaughters();
  for (int i = 0; i < nModules; ++i) {
    auto* daughter = detectorNode->GetDaughter(i);
    auto geoMatrix = daughter->GetMatrix();
    TVector3 translation(geoMatrix->GetTranslation());

    int modId = daughter->GetNumber();
    double x  = translation.X();
    double y  = translation.Y();
    translation.SetZ(frontFaceGlobal.Z());
    double z  = translation.Z();
    modulePosMap.insert({modId, {x,y,z}});
  }

  geoMan->GetListOfVolumes()->Delete();
  geoMan->GetListOfShapes()->Delete();
  delete geoMan;
  nModules=modulePosMap.rbegin()->first;
  vector <XYZVector> modulePosVector(nModules,{0.,0.,0.});
  for(auto &modulePos:modulePosMap)
    modulePosVector.at(modulePos.first-1)=modulePos.second;
  if (verbose)
  {
    printf("%d module positions:\n", nModules);
    for(int i=0;i<nModules;i++)
      printf("%d: (%f, %f, %f)\n", i, modulePosVector.at(i).x(), modulePosVector.at(i).y(), modulePosVector.at(i).z());
  }
  return modulePosVector;
}

RVec<float> fhcalModE(BmnFHCalEvent event)
{
  vector<float> fhcalModEnergy_;
  for (int i = 0; i < 54; i++)
    fhcalModEnergy_.push_back(event.GetModule(i+1)->GetEnergy());
  return fhcalModEnergy_;
}

RVec<float> getEloss(const RVec<int> moduleId, const RVec<double> eLossDigis, const vector<XYZVector> modulePos)
{
  int nModules=modulePos.size();
  vector<float> eLossModules(nModules,0);
  int nDigis=eLossDigis.size();
  for(int i=0;i<nDigis;i++)
    if(moduleId.at(i) <= nModules)
      eLossModules.at(moduleId.at(i)-1)=eLossDigis.at(i);
  return eLossModules;
}

RVec<float> fdEloss (const RVec<BmnFDPoint> points)
{
  vector<float> eLoss;
  for (auto p:points)
    if (getCharge(p.GetPdgId()) > 0)
      eLoss.push_back(p.GetEnergyLoss());
  return eLoss;
}

RVec<bool> hasHitFhcal (RVec<CbmMCTrack> particles)
{
  vector<bool> hasHit;
  for(auto &part:particles)
    hasHit.push_back(part.GetNPoints(kFHCAL)>0);
  return hasHit;
}

RVec<short> modNhits (RVec<short> digiModIds, RVec<short> pointModIds)
{
  vector<short> nHits(digiModIds.size(),0);
  for (auto &pointModId:pointModIds)
    for (short i=0;i<digiModIds.size();i++)
      if(pointModId==digiModIds.at(i))
        nHits.at(i)++;
  return nHits;
}

RVec<float> mcPointEloss(const TClonesArray points)
{
  vector<float>el; 
  for (auto p:points)
    el.push_back(((FairMCPoint*)p)->GetEnergyLoss());
  return el;
}

void convertBmn (string inReco="data/run8/rec.root", string inSim="data/run8/sim.root", const char *inGeo="data/run8/full_geometry.root", const char *out="data/run8/tree.root")
{
  TChain *chainRec=makeChain(inReco, "bmndata");
  TChain *chainSim=makeChain(inSim, "bmndata");
//  cout << "Reco: " << chainRec->GetEntries() << " events\n";
//  cout << "Sim: " << chainSim->GetEntries() << " events\n";
  chainRec->AddFriend(chainSim);
  ROOT::RDataFrame d(*chainRec);

  // read first Run Header if present
  int nEvents = chainRec->GetEntries();

  BmnFieldPar* fieldPar{nullptr};
  chainSim->GetFile()->GetObject( "BmnFieldPar", fieldPar );
  magField = new BmnNewFieldMap(fieldPar);
  magField->Init();
//  DstRunHeader* run_header = (DstRunHeader*) chainRec->GetCurrentFile()->Get("DstRunHeader");
//  if (run_header) {
//    cout << "\n|||||||||||||||| RUN SUMMARY |||||||||||||||" << endl;
//    cout << "||\t\t\t\t\t  ||" << endl;
//    cout << "||   Period:        " << run_header->GetPeriodNumber() << "\t\t\t  ||" << endl;
//    cout << "||   Number:        " << run_header->GetRunNumber() << "\t\t\t  ||" << endl;
//    cout << "||   Start Time:    " << run_header->GetStartTime().AsString("s") << "\t  ||" << endl;
//    cout << "||   End Time:      " << run_header->GetFinishTime().AsString("s") << "\t  ||" << endl;
//    cout << "||   Beam:          A = " << run_header->GetBeamA() << ", Z = " << run_header->GetBeamA() << "\t  ||" << endl;
//    cout << "||   Beam energy:   " << run_header->GetBeamEnergy() << " GeV\t\t  ||" << endl;
//    cout << "||   Target:        A = " << run_header->GetTargetA() << ", Z = " << run_header->GetTargetZ() << "\t  ||" << endl;
//    cout << "||   Field voltage: " << setprecision(4) << run_header->GetMagneticField() << " mV\t\t  ||" << endl;
//    cout << "||\t\t\t\t\t  ||" << endl;
//    cout << "||||||||||||||||||||||||||||||||||||||||||||\n" << endl;
//  }
  auto scwallModPos=modulePos(inGeo,"scwall");
  auto hodoModPos=modulePos(inGeo,"hodo");
  auto fhcalModPos=modulePos(inGeo,"fhcal");
  
  auto dd=d
//    .Range(0,1000)
    .Define("psiRP","MCEventHeader.fRotZ")
    .Define("evtId","DstEventHeader.fEventId")
    .Define("b","MCEventHeader.fB")
    .Define("simVtxX","MCEventHeader.fX")
    .Define("simVtxY","MCEventHeader.fY")
    .Define("simVtxZ","MCEventHeader.fZ")
    .Define("simMom",simMomentum,{"MCTrack"})
    .Define("simPdg", simPdg, {"MCTrack"})
    .Define("simMotherId", simMotherId, {"MCTrack"})
//    .Define("simHasHitFhcal", hasHitFhcal, {"MCTrack"})
//    .Define("simPosStart", simPosStart, {"MCTrack"})
//    .Define("simCharge", simCharge, {"simPdg"})
    .Define("vtxX","PrimaryVertex.fX")
    .Define("vtxY","PrimaryVertex.fY")
    .Define("vtxZ","PrimaryVertex.fZ")
    .Define("vtxNtracks","PrimaryVertex.fNTracks")
    .Define("vtxChi2","PrimaryVertex.fChi2")
    .Define("vtxNdf","PrimaryVertex.fNDF")
    .Define("trMom",trackMomentum,{"BmnGlobalTrack"})
    .Define("trNhits","BmnGlobalTrack.fNhits")
    .Define("trNdf","BmnGlobalTrack.fNDF")
    .Define("trChi2","BmnGlobalTrack.fChi2")
    .Define("trP",trackP,{"BmnGlobalTrack"})
    .Define("trChi2vtx","BmnGlobalTrack.fChi2InVertex")
    .Define("trLength","BmnGlobalTrack.fLength")
    .Define("trCharge",recCharge,{"BmnGlobalTrack"})
    .Define("trDca",recDca,{"BmnGlobalTrack","PrimaryVertex."})
    .Define("trTof400hit","BmnGlobalTrack.fTof1Hit")
    .Define("trTof700hit","BmnGlobalTrack.fTof2Hit")
    .Define("trBetaTof400","BmnGlobalTrack.fBeta400")
    .Define("trBetaTof700","BmnGlobalTrack.fBeta700")
    .Define("trPosLast",recPosLast,{"BmnGlobalTrack"})
    .Define("trPos450",recPos450,{"BmnGlobalTrack"})
    .Define("trSimIndex",recSimIndex,{"BmnGlobalTrack","MCTrack"})
    .Define("stsTrackCovMatrix", covMatrix, { "BmnGlobalTrack", "StsTrack" })
    .Define("stsTrackMagField", magneticField, { "BmnGlobalTrack", "StsTrack", "StsHit" })
    .Define("stsTrackParameters", stsTrackParameters, { "BmnGlobalTrack", "StsTrack" })
    .Define("stsTrackMomentum", stsTrackMomentum, { "StsTrack" })
    .Define("tof400hitPos",tofHitPosition,{"BmnTof400Hit"})
    .Define("tof400hitT","BmnTof400Hit.fTimeStamp")
    .Define("tof400hitL","BmnTof400Hit.fLength")
    .Define("tof400hitResX","BmnTof400Hit.fResX")
    .Define("tof400hitResY","BmnTof400Hit.fResY")
    .Define("tof400hitRefIndex","BmnTof400Hit.fRefIndex")
    .Define("tof400hitResCalc",tofRes,{"BmnGlobalTrack","BmnTof400Hit"})
    .Define("tof700hitPos",tofHitPosition,{"BmnTof700Hit"})
    .Define("tof700hitT","BmnTof700Hit.fTimeStamp")
    .Define("tof700hitL","BmnTof700Hit.fLength")
    .Define("tof700hitResX","BmnTof700Hit.fResX")
    .Define("tof700hitResY","BmnTof700Hit.fResY")
    .Define("tof700hitRefIndex","BmnTof700Hit.fRefIndex")
    .Define("tof700hitResCalc",tofRes,{"BmnGlobalTrack","BmnTof700Hit"})
    .Define("scwallModPos",[scwallModPos](){return scwallModPos;})
    .Define("scwallModId",moduleId, {"scwallModPos"})
    .Define("scwallModQGeant", getEloss, {"ScWallDigit.fCellID", "ScWallDigit.fELoss", "scwallModPos"})
    .Define("scwallModQ", getEloss, {"ScWallDigit.fCellID", "ScWallDigit.fELossDigi", "scwallModPos"})
    .Define("hodoModPos",[hodoModPos](){return hodoModPos;})
    .Define("hodoModId",moduleId, {"hodoModPos"})
    .Define("hodoModQGeant", getEloss, {"HodoDigit.fStripID", "HodoDigit.fELoss", "hodoModPos"})
    .Define("hodoModQ", getEloss, {"HodoDigit.fStripID", "HodoDigit.fELossDigi", "hodoModPos"})
    .Define("fhcalModPos",[fhcalModPos](){return fhcalModPos;})
    .Define("fhcalModId",moduleId, {"fhcalModPos"})
    .Define("fhcalModE",fhcalModE,{"FHCalEvent"})
    .Define("fdQ","Sum(FDPoint.fCharge*FDPoint.fCharge)")
    .Define("fdLight","Sum(FDPoint.fLightYield)")
    .Define("fdEloss", fdEloss, {"FDPoint"})
    .Define("bdModId", "BdDigit.fMod")
    .Define("bdModAmp", "BdDigit.fAmp")
    .Define("bdPointEloss", mcPointEloss, {"BdPoint"})
    .Define("bdPointModId", "BdPoint.nCopy")
    .Define("bdPointPdg", "BdPoint.fPdgId")
    .Define("bdPointIsPrimary", "BdPoint.fIsPrimary")
//    .Define("simdModId", "SiMDDigit.fMod")
//    .Define("simdModAmp", "SiMDDigit.fAmp")
//    .Define("simdPointEloss", mcPointEloss, {"SiMDPoint"})
//    .Define("simdPointModId", "SiMDPoint.nCopy")
//    .Define("simdPointPdg", "SiMDPoint.fPdgId")
//    .Define("simdPointIsPrimary", "SiMDPoint.fIsPrimary")
//    .Define("bdModNhits", modNhits, {"bdModId", "BdPoint.nCopy"})
//    .Define("simdModNhits", modNhits, {"simdModId", "SiMDPoint.nCopy"})
//    .Define("bdM","Sum(BdDigit.fMod>=0)")
//    .Define("simdM","Sum(SiMDDigit.fMod>=0)")
  ;
  dd.Foreach([](uint evtId){if (evtId % 100 == 0) cout << "\r" << evtId;}, {"evtId"}); // progress display 
  cout << endl;

  vector<string> definedNames;
  vector<string> toExclude={/*"scwallModPos","fhcalModPos","hodoModPos"*/};
  for (auto& definedName:dd.GetDefinedColumnNames())
  {
    bool exclude=false;
    for (auto &nameToExclude:toExclude)
      if (definedName==nameToExclude)
        exclude=true;
    if (!exclude)
      definedNames.push_back(definedName);
  }
  dd.Snapshot("t",out, definedNames);
}
