using namespace ROOT;
using namespace ROOT::Math;
using namespace ROOT::RDF;
using fourVector=LorentzVector<PtEtaPhiE4D<double>>;

Bool_t GetBmnGeom(const char *fileName)
{
   // ---> Get TGeoManager and top node (TOP)

   TFile *geoFile = new TFile(fileName, "READ");
   if (!geoFile->IsOpen()) {
      cout << "-E- Could not open file!" << endl;
      return kFALSE;
   }
   geoFile->FindObjectAny("SIGEMS"); //"GEMS_geom SIGEMS"); //AZ
   return kTRUE;                     // AZ

   // ---> Get TOP node from file
   TList *     keyList = geoFile->GetListOfKeys();
   TIter       next(keyList);
   TKey *      key = NULL;
   TGeoVolume *top = NULL;
   while ((key = (TKey *)next())) {
      TString className(key->GetClassName());
      if (className.BeginsWith("TGeoVolume")) {
         top = dynamic_cast<TGeoVolume *>(key->ReadObj());
         std::cout << "Found volume " << top->GetName() << endl;
         break;
      }
   }
   if (!top) {
      cout << "-E- Could not find volume object in file" << endl;
      return kFALSE;
   }
   // cout << gGeoManager << endl;
   gGeoManager->GetListOfVolumes()->ls();
   gGeoManager->SetTopVolume(top);

   return kTRUE;
}

int GetNofModules(TGeoNode * station) {
    int nModules = 0; //station->GetNdaughters();
    // --- Modules
    for (int iModule = 0; iModule < station->GetNdaughters(); iModule++) {
        TGeoNode* module = station->GetDaughter(iModule);
        if (TString(module->GetName()).Contains("module")) nModules++;
    }
    return nModules;
}

void ApplyAlignment()
{
  const Int_t ip_r7[6][2]={0,1, 2,3, 4,5, 6,7, 8,9, 10,11};
  // Apply alignment

  TGeoNode* st_node = NULL;
  // Get STS node
  TGeoNode* sts = NULL;
  gGeoManager->CdTop();
  TGeoNode* cave = gGeoManager->GetCurrentNode();
  for (int iNode = 0; iNode < cave->GetNdaughters(); iNode++) {
    TGeoNode* node = cave->GetDaughter(iNode);
    TString name = node->GetName();
    if ( name.Contains("GEMS_0") ) {
      sts = node;
      gGeoManager->CdDown(iNode);
      break;
    }
  }
  if ( ! sts ) {
    cout << "-E- CbmBmnStsDigitize::ApplyAlignment: Cannot find top GEM node"
	 << endl;
    //return kFALSE;
  }

  // Loop over stations in STS
  int statNr = 0;
  int nStat0 = sts->GetNdaughters();
  int nMod0sts=0;
  
  cout<<"STATIONS : "<<nStat0<<endl;
  const double xAlign[12] = {-40.8258, 40.8258, 40.845, -40.845, 40.8608, -40.8608, -40.8578, 40.8578, 40.845, -40.845, -40.8298, 44.8606}; //for 4 GeV
 //double xAlign[] = {40.9414-0.1156,-40.7102-0.1156,  -40.5867-0.2583,41.1033-0.2583,   -39.4213-1.4395,42.3003-1.4395, //after hZKA00Zd.root
//			      43.2731-2.4153,-38.4425-2.4153,  -37.5865-3.2585,44.1035-3.2585, 44.8606-4.03085, -36.7989-4.03085}; //for 4 GeV
                  //9.476, 3.269,-2.7926,-8.796,                     9.476,3.2695,-2.791,-8.7056}; //for 3.5 GeV
   //{9.476,  -8.556,   -8.5656,  9.476,    -2.5526,  3.4695,    -2.551,   3.469};
  //{0, 0, -0.2, 0,       27.94,-34.38,      41.26,-40.44,      41.59,-40.12}; //run5
  //0.24, 0.32,     0.32, 0.50,      0.83, 0.64};//-1.5, -1.70, -1.5,   27.94,-34.38,      -40.44,41.26,      41.59,-40.12};
  const  double yAlign[12] = {22.2753, 22.4706, 22.4851, 22.6689, 22.7851, 22.5916, //after hZKA00Zd.root
			     22.3308, 23.1519, 22.4786, 22.4484, 22.5971, 23.2184};
    //{0,  0, 0.20, 0 ,     0.19,-0.625,      0.31,0.13,      0.276,0.369};//run5
  //0.19, -0.625,     0.13, 0.31,    0.267, 0.369 };//0,    0.20,   0,    0.19,-0.425,      -2.27,-2.09,      -2.133,-2.031};
 // double zAlign2[] = {0, 0.9, 1.6, 0.9,      1.9, 3.9,      2.2,1.7 ,     2.7, 1.7};
/*
 // double zAlign[] = {0, 31.8, 63.8, 95.8, 127.8,159.8,191.8};
  double zAlign3[] = {2.47, 32.85,65.35,96.67,130.635,161.47,194.35};// 2.47};
  double zAlign23[] = {0, 0.0, 0.0, 0.0,-1.03,1.3,0.10,-0.10,0.10,-0.25,
//  double zAlign3[] = {2.47, 32.85,65.25,96.65,130.95,161.45,194.25};// 2.47};
//  double zAlign23[] = {0, 0.0, 0.0, 0.0,-1.2,1.2,0.10,-0.10,0.10,-0.10,
  //3.1, 1.35, 3.59, 1.84, 2.88, 3.37, 2.06, 1.57};   
  0.63,-0.9,0.41,-1.12,   -0.63,0.9,-0.41,   1.12 };//1.12, -0.41, 0.9, -0.63, -1.12, 0.41,-0.9,0.63};         // 0.63,-0.9, 0.41,-1.12,    -0.63,0.9,-0.41,1.12 }; // 
//arZ0[MQP]={ 3.1, 1.35, 3.59, 1.84, 2.88, 3.37, 2.06, 1.57};
 // old
 
*/
   //42.28635/*43.78635*/, 67.9045/*69.4045*/,115.7825/*117.2825*/,138.1185/*139.6185*/163.8475 /*165.3475*/, 186.5085/*188.0085*/
   const double zAlign3[6] ={42.28635,     67.9045,   115.7825, 138.1185,         163.8475,      186.5085};// 2.47}; // for 4 GeV
   //{0,0,0,0,0,0};
    //{43.78635,     69.4045,   117.2825, 139.6185,         165.3475,      188.0085};// 2.47}; // for 4 GeV

              // double zAlign3[] = {2.47, 32.87,65.23,96.67,130.5,161.27,193.83};// 2.47}; // for MC D. Bar
 //double zAlign3[] = {2.47, 32.15,64.65,96.65,130.35,161.25,193.45};// 2.47}; // for MC D. Bar exp value
 //double zAlign3[] = {2.47, 31.87,64.75,96.65,130.5,161.3,193.6};// 2.47}; // for MC D. Bar exp value 2
 //128.65; 131.15
 // double zAlign3[] = {2.47, 32.85,65.35,96.65,130.95,161.45,194.25};// 2.47}; //for 3.5 GeV
  //  double zAlign3[] = {2.47, 32.85,65.25,96.65,130.95,161.45,194.25};// 2.47};
  double zAlign23[] = { 0.14595, -0.14595, -0.1258, 0.1258, -0.2005, 0.2005, // Z after hZKA00Zc.root
			     0.0945, -0.0945, -0.1125, 0.1125, 0.0865, -0.0865 }; //for 4 GeV
 
 // double zAlign23[] = {0, 0.0, 0.0, 0.0,    -1.25,1.25,    0.,0.,0.,0., //-1.45,1.05
 //   0.,0.,0.,0,   0,0.,0,   0 }; // for MC D. Bar



 double driftcorr= 0;//-0.2;
int stn=0;
  //for (int iNode = 1; iNode < sts->GetNdaughters(); iNode++) {
  for (int iNode = 0; iNode < nStat0; iNode++) {
    // Go to station node
    gGeoManager->CdDown(iNode);
    TGeoNode* stationNode = gGeoManager->GetCurrentNode();
    TString statName = stationNode->GetName();
    if ( ! statName.Contains("station") ) {
      gGeoManager->CdUp();
      continue;
    }
    //AZ int statNr = stationNode->GetNumber();
    ++stn; //AZ

    if(stn<4){
gGeoManager->CdUp();
 continue;}

 ++ statNr;

    TGeoHMatrix *matr = gGeoManager->GetCurrentMatrix();
    double* statTrans = matr->GetTranslation();
    //double* statRot = matr->GetRotationMatrix();
    TGeoRotation *r2;
   TGeoTranslation *t2 = new TGeoTranslation(statTrans[0],statTrans[1],statTrans[2]);
  //if( statNr==2 || statNr==3 || statNr==5)r2 = new TGeoRotation("rot",0,180,180);
   // if(statNr==1 || statNr==4 || statNr==6) r2 = new TGeoRotation("rot",0,180,180);
     r2 = new TGeoRotation("rot",0,0,0);

      TGeoCombiTrans *cc1 = new TGeoCombiTrans(*t2,*r2);
    //  statTrans[0] = xAlign[statNr];
    //  statTrans[1] = yAlign[statNr];

    //  cout<<"ST TRANS Z: "<<statTrans[2]<<endl;
  //if(statNr-1 >1)  statTrans[2] = zAlign3[statNr-1]+driftcorr;//zAlign[statNr];
  cout<<"statNr : "<<statNr <<", "<<statTrans[2]  <<" zal: "<<zAlign3[statNr-1]<<endl;
   statTrans[2] = zAlign3[statNr-1];//zAlign[statNr];
   //if(statNr-1 >1)  statTrans[2]=statTrans[2]+driftcorr;
    //cout<<"ST TRANS Z: "<<statTrans[2]<<endl;
    //if(iNode>0) statTrans[2]=statTrans[2]-2.47;
    //matr->SetTranslation(statTrans);
    //matr->SetMatrix(statRot);

    TGeoHMatrix *matr0 = new TGeoHMatrix(*cc1);
    matr0->RegisterYourself();
    
    //  int nModules = stationNode->GetNdaughters();
    int nModules = GetNofModules(stationNode);
    //  cout<<"nModules: "<<nModules<<endl;
    
    //sts->GetVolume()->ReplaceNode(stationNode,0,gGeoManager->GetCurrentMatrix()); //AZ
    //sts->GetVolume()->RemoveNode(stationNode); //AZ
    sts->GetVolume()->AddNode((TGeoVolumeAssembly*)stationNode->GetVolume(),0,matr0); //AZ

    //AZ- hard-coded st_node=sts->GetVolume()->GetNode(iNode+6);
    st_node = (TGeoNode*) sts->GetVolume()->GetNodes()->Last(); //AZ

   // double  statZ = statTrans[2];
   //cout <<"sta: " << statNr << " " << gGeoManager->GetCurrentMatrix()->GetTranslation()[2] << " " << sts->GetNdaughters() << endl;

    //gGeoManager->CdUp();               // to sts

    //-----------------------module translate-------------
    int moduleNr = 0, copy_no = 0;
    //cout<<"nMODULES: "<< nModules<<endl;
    /*if(iNode==0){
      nMod0sts=nModules;
      nModules=nModules*2;
      
      }*/
    // ---> Large sensors
    for (int iStatD = 0; iStatD < nModules; iStatD++) {
      gGeoManager->CdDown(iStatD);

      TGeoNode* module = gGeoManager->GetCurrentNode();
      if ( ! TString(module->GetName()).Contains("module") ) {
	gGeoManager->CdUp();
	continue;
      }
cout << iStatD << " " << module->GetName() << " stn: "<<statNr << endl; 
      
 if(/*iNode>0*/ true){
      if (TString(module->GetName()).Contains("Senso")) {
	//if(iNode>0 && iNode<=3) moduleNr=0;
	//else{
	  if(iStatD==0) moduleNr=0;
	  if(iStatD>0) moduleNr=1;
//	}
	//moduleNr=iStatD;
	//fprintf(parFile, "%4d %4d\n", moduleNr, 1);
	//SaveSensor(geoMan, parFile, phiStat, module);

	int ipNr = 0;
  //cout<<" iNode: "<< iNode<<endl;
	//if (iNode>0 && iNode<=3) ipNr = iNode;
	//else ipNr = iNode*2 - 4 + moduleNr;
//int order=-1;
//if( (statNr==1 || statNr==4 || statNr==6) && moduleNr==0) order=1; 
  

	TGeoHMatrix *matrMod = gGeoManager->GetCurrentMatrix();
	double* modTrans = matrMod->GetTranslation();
  //cout<<" modtransX: "<<modTrans[0]<<" xali: "<<xAlign[ip_r7[statNr-1][moduleNr]]<< "sub: "<<xAlign[ip_r7[statNr-1][moduleNr]]-modTrans[0]  <<endl;
	 modTrans[0] =xAlign[ip_r7[statNr-1][moduleNr]];
	 modTrans[1] = yAlign[ip_r7[statNr-1][moduleNr]] - modTrans[1];
   modTrans[2] = zAlign23[ip_r7[statNr-1][moduleNr]];//zAlign2[ipNr];
  //cout<<"ST TRANS Z: "<<zAlign23[ipNr]<<endl;
	matrMod->SetTranslation(modTrans);
	TGeoHMatrix *matr0Mod = new TGeoHMatrix(*matrMod);
	matr0Mod->RegisterYourself();
  ipNr++;
//cout<<" ip N: "<< ipNr<<endl;
	//sts->GetVolume()->ReplaceNode(stationNode,0,gGeoManager->GetCurrentMatrix()); //AZ
	//sts->GetVolume()->RemoveNode(stationNode); //AZ
	//  stationNode->GetVolume()->AddNode((TGeoVolumeAssembly*)module->GetVolume(),0,matr0Mod); //AZ
	
	//cout<<" 1 st name add: "<< stationNode->GetName()<<" mod name add: "<<module->GetName()<< " "<< module<< " i: "<<iStatD<<" cols: "<<st_node->GetVolume()->GetNdaughters()<<endl;
	//double* sensTrans = matrMod->GetTranslation();
	//cout<<"trans: "<<sensTrans[0]<<" "<<sensTrans[1]<<" "<<sensTrans[2]<< " Nr mod:  "<<moduleNr<<endl;
	//stationNode->GetVolume()->AddNode((TGeoVolumeAssembly*)module->GetVolume(),copy_no,matr0Mod);
	st_node->GetVolume()->AddNode((TGeoVolumeAssembly*)module->GetVolume(),copy_no,matr0Mod);
	//cout<<" 1 st name add: "<< stationNode->GetName()<<" mod name add: "<<module->GetName() << " "<< module<< " i: "<<iStatD<<" cols: "<<st_node->GetVolume()->GetNdaughters()<<endl;
	//double  modZ = modTrans[2];
	//cout <<"mod: " << ipNr << " VEC: " << modTrans[0] << " "<<modTrans[1] << " "<<modTrans[2] << endl;
	copy_no++;
	//delete matr0Mod;
	
      }

      gGeoManager->CdUp();  // back to module

    }} // for (int iStatD = 0; iStatD < nModules;
    //----------------------end module translate----------
    //delete matr0;
    gGeoManager->CdUp();               // to sts
  }                                    // station loop

  int snr=1;
  vector < TGeoNode* > removNodes;
  // Remove extra nodes
  for (int iNode = 0; iNode < nStat0; iNode++) {
    // Go to station node
    //gGeoManager->CdDown(1);
      gGeoManager->CdDown(iNode);
    TGeoNode* stationNode = gGeoManager->GetCurrentNode();
     cout<<" st name del: "<< stationNode->GetName()<< " iNode: "<<iNode<<endl;
    if(iNode>2) removNodes.push_back(stationNode);//sts->GetVolume()->RemoveNode(stationNode); //AZ
    gGeoManager->CdUp();               // to sts
  } 
  for(int o=0; o<removNodes.size(); o++)
  sts->GetVolume()->RemoveNode(removNodes[o]);
removNodes.clear();



  for (int iNode = 0; iNode <nStat0; iNode++) { 

    // Go to station node
    //gGeoManager->CdDown(iNode);
    gGeoManager->CdDown(iNode);
    TGeoNode* stationNode = gGeoManager->GetCurrentNode();
   if(iNode>2){
   int  nMod = 2;
    for (int iStatD = 0; iStatD < nMod; iStatD++) {
      gGeoManager->CdDown(0);//stationNode->GetNdaughters()-1);
      
      TGeoNode* module = gGeoManager->GetCurrentNode();
      //  cout<<" 2 st name del: "<< stationNode->GetName()<<" mod name del: "<<module->GetName() << " i: "<<iStatD<<endl;
      stationNode->GetVolume()->RemoveNode(module); //AZ
      
      gGeoManager->CdUp();               // to sts
    } 
   }
    gGeoManager->CdUp();
  } 
  
  for (int iNode = 0; iNode < sts->GetNdaughters(); iNode++) {

    // Go to station node
    //gGeoManager->CdDown(iNode);
    gGeoManager->CdDown(iNode);
    TGeoNode* stationNode = gGeoManager->GetCurrentNode();
    cout<<" Check-in STATION: "<< stationNode->GetName()<<" zpos : "<<gGeoManager->GetCurrentMatrix()->GetTranslation()[2]<<endl;
    for (int iModule = 0; iModule < stationNode->GetNdaughters() ; iModule++) {
      gGeoManager->CdDown(iModule);
      TGeoNode* moduleNode = gGeoManager->GetCurrentNode();
          cout<<" Check-in st name : "<< stationNode->GetName()<<" mod name : "<<moduleNode->GetName()<<endl;
      double* sensTrans = gGeoManager->GetCurrentMatrix()->GetTranslation();
        cout<<"trans mod: "<<sensTrans[0]<<" "<<sensTrans[1]<<" "<<sensTrans[2]<< endl;
      //  stationNode->GetVolume()->RemoveNode(moduleNode); //AZ
      gGeoManager->CdUp();
    }
    
    gGeoManager->CdUp();               // to sts
  }
 
  
//exit(0);
//gGeoManager->SetVisLevel(500);
//gGeoManager->GetTopVolume()->SetTransparency(0);
//gGeoManager->GetTopVolume()->Draw("ogl");

//int c;
//cin>>c;

}
template <typename T>
int vsize(const RVec<T> vec)
{
  return vec.size();
}

int arrSize(TClonesArray arr)
{
  return arr.GetSize();
}

vector<fourVector> trMomentum(const RVec<CbmStsTrack> tracks)
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

RVec<CbmStsTrack> extrapolateToVertex(const RVec<CbmStsTrack> tracks, const CbmVertex vtx)
{
  vector<CbmStsTrack> extrapolated;
  for (auto track:tracks) {
    CbmKFTrack kfTrack = CbmKFTrack(track,0);
    kfTrack.Extrapolate(vtx.GetZ());
    FairTrackParam parFirst; // Could realy be parLast. See CbmKFTrack() constructor
    kfTrack.GetTrackParam(parFirst);
    track.SetParamFirst(parFirst);
    extrapolated.push_back(track);
  }
  return extrapolated;
}

RVec<int> trCharge(const RVec<CbmStsTrack> tracks)
{
  vector<int> charge;
  for (auto track:tracks) {
    int q = track.GetParamFirst()->GetQp() > 0 ? 1 : -1;
    charge.push_back(q);
  }
  return charge;
}

vector<XYZVector> trDca(const RVec<CbmStsTrack> tracks, const CbmVertex vtx)
{
  vector<ROOT::Math::XYZVector> dca;
  for (auto track:tracks) {
    auto par = track.GetParamFirst();
    dca.push_back({par->GetX()-vtx.GetX(),par->GetY()-vtx.GetY(),par->GetZ()-vtx.GetZ()});
  }
  return dca;
}

RVec<float> trChi2 (const RVec<CbmStsTrack> tracks)
{
  vector<float> chi2;
  for (auto &t:tracks)
    chi2.push_back(t.GetChi2());
  return chi2;
}

RVec<float> trNdf (const RVec<CbmStsTrack> tracks)
{
  vector<float> ndf;
  for (auto &t:tracks)
    ndf.push_back(t.GetNDF());
  return ndf;
}

RVec<float> trNhits (const RVec<CbmStsTrack> tracks)
{
  vector<float> nhits;
  for (auto &t:tracks)
    nhits.push_back(t.GetNStsHits());
  return nhits;
}

BmnEventHeader *header(TClonesArray a){return (BmnEventHeader*)a.At(0);}

void convertBmn_run7 (const char *inTracks="data/run7/plotnikov/3819.root", 
                      const char *inTrigger="data/run7/trigger/3819.root", 
                      const char *out="data/run7/tree/3819.tree.root")
{
  BmnFieldMap* magField = new BmnNewFieldMap("field_sp41v5_ascii_Extrap.dat");
  magField->SetScale(1251.9/900.);
  magField->Init();
  magField->Print();
  if (!GetBmnGeom("SIGEMS_r7_v9.root"))
    return;
  ApplyAlignment();
  auto *run = new FairRunAna();
  run->SetField(magField);
  run->SetSource(new BmnFileSource(inTracks));
  run->SetOutputFile(Form("%s_tmp", out));

  FairTask *pKF = new CbmKF();
  run->AddTask(pKF);
  run->Init();

  TChain chain("BMN_DIGIT");
  chain.Add(inTracks);
  TChain chTrigger("BMN_DIGIT");
  chTrigger.Add(inTrigger);
  chain.AddFriend(&chTrigger);
  ROOT::RDataFrame d(chain);
  // read first Run Header if present
  int nEvents=chain.GetEntries();
  auto dd=d
//    .Range(0,1000)
    .Define("evtId", [](TClonesArray a){return ((BmnEventHeader*)a.At(0))->GetEventId();}, {"EventHeaderBmn"})
    .Define("runId", [](TClonesArray a){return ((BmnEventHeader*)a.At(0))->GetRunId();}, {"EventHeaderBmn"})
    .Define("vtxX","PrimaryVertex.fX")
    .Define("vtxY","PrimaryVertex.fY")
    .Define("vtxZ","PrimaryVertex.fZ")
    .Define("vtxNtracks","PrimaryVertex.fNTracks")
    .Define("vtxChi2","PrimaryVertex.fChi2")
    .Define("vtxNdf","PrimaryVertex.fNDF")
    .Define("trNhits",trNhits,{"StsTrack"})
    .Define("trNdf","StsTrack.fNDF")
    .Define("trChi2",trChi2,{"StsTrack"})
    .Define("trCharge",trCharge,{"StsTrack"})
    .Define("trPidHypo","StsTrack.fPidHypo")
    .Define("StsTrackExtrap",extrapolateToVertex,{"StsTrack", "PrimaryVertex."})
    .Define("trMom",trMomentum,{"StsTrackExtrap"})
    .Define("trDca",trDca,{"StsTrackExtrap","PrimaryVertex."})
//    .Define("trExtrap",trMomentum,{"StsTrackExtrap"})
//    .Define("trExtrapDca",trDca,{"StsTrackExtrap","PrimaryVertex."})
//    .Define("trDxCsc","StsTrackCSC.fDxCsc")
//    .Define("trDyCsc","StsTrackCSC.fDyCsc")
  ;
  dd.Foreach([nEvents](uint evtId){if (evtId % 1000 == 0) cout << "\r" << evtId << "/" << nEvents;}, {"evtId"}); // progress display 
  cout << endl;
  
  vector<string> definedNames;
  vector<string> toExclude={"StsTrackExtrap"};
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
