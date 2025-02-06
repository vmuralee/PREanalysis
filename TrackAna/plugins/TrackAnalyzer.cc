// -*- C++ -*-
//
// Package:    TrackAnalyzer/TrackAnalyzer
// Class:      TrackAnalyzer
//
/**\class TrackAnalyzer TrackAnalyzer.cc TrackAnalyzer/TrackAnalyzer/plugins/TrackAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Vinaya Krishna Nair
//         Created:  Fri, 26 Jul 2024 11:16:33 GMT
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/SiStripDetId/interface/SiStripDetId.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit1D.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2D.h"
#include "DataFormats/SiStripCluster/interface/SiStripClusterTools.h"
#include "DataFormats/TrackerRecHit2D/interface/TrackerSingleRecHit.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1.h"
//
// class declaration
//
//ROOT inclusion
#include "TROOT.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TTree.h"
#include "TMath.h"
#include "TList.h"
#include "TString.h"

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

using reco::TrackCollection;
//using reco::PFJetCollection;

class TrackAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit TrackAnalyzer(const edm::ParameterSet&);
  ~TrackAnalyzer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  //void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  //void endJob() override;
  void SetVarToZero();
  // Example function to calculate deposited charge
  double calculateDepositedCharge(const PSimHit& hit);

  // ----------member data ---------------------------
  bool isMC_;
  edm::EDGetTokenT<reco::TrackCollection> tracksToken_; 
  edm::EDGetTokenT<reco::PFJetCollection> jetToken_;
  edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> m_geomToken_;
  edm::EDGetTokenT<std::vector<PSimHit>> g4SimHitsHighTofTokenTOB_;
  edm::EDGetTokenT<std::vector<PSimHit>> g4SimHitsLowTofTokenTOB_;

  edm::EDGetTokenT<std::vector<PSimHit>> g4SimHitsHighTofTokenTIB_;
  edm::EDGetTokenT<std::vector<PSimHit>> g4SimHitsLowTofTokenTIB_;

  edm::EDGetTokenT<std::vector<PSimHit>> g4SimHitsHighTofTokenTEC_;
  edm::EDGetTokenT<std::vector<PSimHit>> g4SimHitsLowTofTokenTEC_;

  edm::EDGetTokenT<std::vector<PSimHit>> g4SimHitsHighTofTokenTID_;
  edm::EDGetTokenT<std::vector<PSimHit>> g4SimHitsLowTofTokenTID_;

  edm::EDGetTokenT<std::vector<SimTrack>> g4SimHitsTrackToken_;

  TTree* simTreeTOBLoTof;
  TTree* simTreeTOBHiTof;

  TTree* simTreeTIBLoTof;
  TTree* simTreeTIBHiTof;

  TTree* simTreeTECLoTof;
  TTree* simTreeTECHiTof;

  TTree* simTreeTIDLoTof;
  TTree* simTreeTIDHiTof;
  TTree* trackTree;

  edm::Service<TFileService> fs;

  edm::EventNumber_t eventN;
  const static int nMax = 40000000;
  const static int nClus = 4000;

  uint16_t nTracks;
  int nJets;

  int runN;
  int lumi;

  
  float clus_chargePerCM;
  int clus_charge;
  int clus_width;
  int clus_avgcharge;
  int clus_detId;
  uint16_t clus_firstStrip;
  uint16_t clus_endStrip;
  uint8_t subdet;

  float clus_barycenter;

  //simhit TOB
  int simhitTOB_hightof_charge;
  int simhitTOB_lowtof_charge;
  int simhitTOB_charge;
  float simhitTOB_hightof_chargePerCM;
  float simhitTOB_lowtof_chargePerCM;

  int simhitTOB_width;
  int simhitTOB_avgcharge;
  int simhitTOB_hightof_detId;
  int simhitTOB_lowtof_detId;
  float simhitTOB_hightof_pT;
  float simhitTOB_lowtof_pT;
  float simhitTOB_pT;

  float simhitTOB_x;
  float simhitTOB_y;
  float simhitTOB_z;
  
  uint16_t simhitTOB_firstStrip;
  uint16_t simhitTOB_endStrip;
  float simhitTOB_barycenter;


  //simhit TIB
  int simhitTIB_hightof_charge;
  int simhitTIB_lowtof_charge;
  int simhitTIB_charge;
  float simhitTIB_hightof_chargePerCM;
  float simhitTIB_lowtof_chargePerCM;

  int simhitTIB_width;
  int simhitTIB_avgcharge;
  int simhitTIB_hightof_detId;
  int simhitTIB_lowtof_detId;
  float simhitTIB_hightof_pT;
  float simhitTIB_lowtof_pT;
  float simhitTIB_pT;

  float simhitTIB_x;
  float simhitTIB_y;
  float simhitTIB_z;

  //simhit TEC
  int simhitTEC_hightof_charge;
  int simhitTEC_lowtof_charge;
  int simhitTEC_charge;
  float simhitTEC_hightof_chargePerCM;
  float simhitTEC_lowtof_chargePerCM;

  int simhitTEC_width;
  int simhitTEC_avgcharge;
  int simhitTEC_hightof_detId;
  int simhitTEC_lowtof_detId;
  float simhitTEC_hightof_pT;
  float simhitTEC_lowtof_pT;
  float simhitTEC_pT;

  float simhitTEC_x;
  float simhitTEC_y;
  float simhitTEC_z;
  
  uint16_t simhitTEC_firstStrip;
  uint16_t simhitTEC_endStrip;
  float simhitTEC_barycenter;
  
  //simhit TID
  int simhitTID_hightof_charge;
  int simhitTID_lowtof_charge;
  int simhitTID_charge;
  float simhitTID_hightof_chargePerCM;
  float simhitTID_lowtof_chargePerCM;

  int simhitTID_width;
  int simhitTID_avgcharge;
  int simhitTID_hightof_detId;
  int simhitTID_lowtof_detId;
  float simhitTID_hightof_pT;
  float simhitTID_lowtof_pT;
  float simhitTID_pT;

  float simhitTID_x;
  float simhitTID_y;
  float simhitTID_z;

  float hit_pt;
  float hit_eta;
  float hit_phi;

  float trkPt[nMax];
  float trkEta[nMax];
  float trkPhi[nMax];
  float trkDxy1[nMax];
  float trkDxyError1[nMax];
  float trkDz1[nMax];
  float trkDzError1[nMax];
  unsigned char trkAlgo[nMax];
  int  trkNHit[nMax];
  int  trkNdof[nMax];
  uint8_t trkNlayer[nMax];

  float trkChi2[nMax];
  float trkNChi2[nMax];
  float trkPtError[nMax];
  
  float jetPt[nMax];
  float jetEta[nMax];
  float jetPhi[nMax];
  float jetMass[nMax];

  bool fromTrack;  
};


TrackAnalyzer::TrackAnalyzer(const edm::ParameterSet& iConfig){

  isMC_ = iConfig.getUntrackedParameter<bool>("IsMC", false);
  tracksToken_ = consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("tracks"));
  jetToken_   = consumes<reco::PFJetCollection>(iConfig.getParameter<edm::InputTag>("jets"));
  if(isMC_){
    g4SimHitsLowTofTokenTOB_  = consumes<std::vector<PSimHit>>(iConfig.getParameter<edm::InputTag>("g4SimHitsLowTofTagTOB"));
    g4SimHitsHighTofTokenTOB_ = consumes<std::vector<PSimHit>>(iConfig.getParameter<edm::InputTag>("g4SimHitsHighTofTagTOB"));

    g4SimHitsLowTofTokenTIB_  = consumes<std::vector<PSimHit>>(iConfig.getParameter<edm::InputTag>("g4SimHitsLowTofTagTIB"));
    g4SimHitsHighTofTokenTIB_ = consumes<std::vector<PSimHit>>(iConfig.getParameter<edm::InputTag>("g4SimHitsHighTofTagTIB"));

    g4SimHitsLowTofTokenTEC_  = consumes<std::vector<PSimHit>>(iConfig.getParameter<edm::InputTag>("g4SimHitsLowTofTagTEC"));
    g4SimHitsHighTofTokenTEC_ = consumes<std::vector<PSimHit>>(iConfig.getParameter<edm::InputTag>("g4SimHitsHighTofTagTEC"));

    g4SimHitsLowTofTokenTID_  = consumes<std::vector<PSimHit>>(iConfig.getParameter<edm::InputTag>("g4SimHitsLowTofTagTID"));
    g4SimHitsHighTofTokenTID_ = consumes<std::vector<PSimHit>>(iConfig.getParameter<edm::InputTag>("g4SimHitsHighTofTagTID"));

    g4SimHitsTrackToken_ = consumes<std::vector<SimTrack>>(iConfig.getParameter<edm::InputTag>("g4SimHitsTrackTag"));   
  }
  
  m_geomToken_ = esConsumes(); 
  usesResource("TFileService");

  trackTree = fs->make<TTree>("trackTree","trackTree");
  simTreeTOBLoTof = fs->make<TTree>("simTreeTOBLoTof","simTreeTOBLoTof");
  simTreeTOBHiTof = fs->make<TTree>("simTreeTOBHiTof","simTreeTOBHiTof");

  simTreeTIBLoTof = fs->make<TTree>("simTreeTIBLoTof","simTreeTIBLoTof");
  simTreeTIBHiTof = fs->make<TTree>("simTreeTIBHiTof","simTreeTIBHiTof");

  simTreeTECLoTof = fs->make<TTree>("simTreeTECLoTof","simTreeTECLoTof");
  simTreeTECHiTof = fs->make<TTree>("simTreeTECHiTof","simTreeTECHiTof");

  simTreeTIDLoTof = fs->make<TTree>("simTreeTIDLoTof","simTreeTIDLoTof");
  simTreeTIDHiTof = fs->make<TTree>("simTreeTIDHiTof","simTreeTIDHiTof");

  trackTree->Branch("event", &eventN, "event/i");
  trackTree->Branch("run",   &runN, "run/I");
  trackTree->Branch("lumi",  &lumi, "lumi/I");

  trackTree->Branch("nTracks",  &nTracks, "nTracks/I");
  trackTree->Branch("trkPt",  trkPt, "trkPt[nTracks]/F");
  trackTree->Branch("trkEta",  trkEta, "trkEta[nTracks]/F");
  trackTree->Branch("trkPhi",  trkPhi, "trkPhi[nTracks]/F");
  trackTree->Branch("trkDxy1",  trkDxy1, "trkDxy1[nTracks]/F");
  trackTree->Branch("trkDz1",  trkDz1, "trkDz1[nTracks]/F");
  
  trackTree->Branch("trkDxyError1",  trkDxyError1, "trkDxyError1[nTracks]/F");
  trackTree->Branch("trkDzError1",  trkDzError1, "trkDzError1[nTracks]/F");
  
  trackTree->Branch("trkChi2",  trkChi2, "trkChi2[nTracks]/F");
  trackTree->Branch("trkNChi2",  trkNChi2, "trkNChi2[nTracks]/F");
  trackTree->Branch("trkPtError",  trkPtError, "trkPtError[nTracks]/F");

  trackTree->Branch("trkAlgo",trkAlgo,"trkAlgo[nTracks]/s");
  trackTree->Branch("trkNHit",trkNHit,"trkNHit[nTracks]/I");
  trackTree->Branch("trkNdof",trkNdof,"trkNdof[nTracks]/I");
  trackTree->Branch("trkNlayer",trkNlayer,"trkNlayer[nTracks]/I");
  
  trackTree->Branch("nJets",  &nJets, "nJets/I");
  trackTree->Branch("jetPt",  jetPt, "jetPt[nJets]/F");
  trackTree->Branch("jetEta",  jetEta, "jetEta[nJets]/F");
  trackTree->Branch("jetPhi",  jetPhi, "jetPhi[nJets]/F");
  trackTree->Branch("jetMass",  jetMass, "jetMass[nJets]/F");

  // cluster variable
  trackTree->Branch("hit_eta",&hit_eta,"hit_eta/F");
  trackTree->Branch("hit_phi",&hit_phi,"hit_phi/F");
  trackTree->Branch("hit_pt",&hit_pt,"hit_pt/F");

  trackTree->Branch("clus_chargePerCM", &clus_chargePerCM, "clus_chargePerCM/F");
  trackTree->Branch("clus_charge",&clus_charge,"clus_charge/I");
  trackTree->Branch("clus_width",&clus_width,"clus_width/I");
  trackTree->Branch("clus_avgcharge",&clus_avgcharge,"clus_avgcharge/I");
  trackTree->Branch("clus_detId",&clus_detId,"clus_detId/I");
  trackTree->Branch("clus_firstStrip",&clus_firstStrip,"clus_firstStrip/I");
  trackTree->Branch("clus_endStrip",&clus_endStrip,"clus_endStrip/I");
  trackTree->Branch("clus_barycenter",&clus_barycenter,"clus_barycenter/F");

  trackTree->Branch("subdet",&subdet,"subdet/I");
  //sim TOB cluster High Tof
  simTreeTOBHiTof->Branch("simhitTOB_hightof_charge",&simhitTOB_hightof_charge,"simhitTOB_hightof_charge/I");
  simTreeTOBHiTof->Branch("simhitTOB_hightof_pT",&simhitTOB_hightof_pT,"simhitTOB_hightof_pT/F");
  simTreeTOBHiTof->Branch("simhitTOB_hightof_detId",&simhitTOB_hightof_detId,"simhitTOB_hightof_detId/I");
  simTreeTOBHiTof->Branch("simhitTOB_hightof_chargePerCM",&simhitTOB_hightof_chargePerCM,"simhitTOB_hightof_chargePerCM/F");
  //sim TOB cluster Low Tof
  simTreeTOBLoTof->Branch("simhitTOB_lowtof_charge",&simhitTOB_lowtof_charge,"simhitTOB_lowtof_charge/I");
  simTreeTOBLoTof->Branch("simhitTOB_lowtof_pT",&simhitTOB_lowtof_pT,"simhitTOB_lowtof_pT/F");
  simTreeTOBLoTof->Branch("simhitTOB_lowtof_detId",&simhitTOB_lowtof_detId,"simhitTOB_lowtof_detId/I");
  simTreeTOBLoTof->Branch("simhitTOB_lowtof_chargePerCM",&simhitTOB_lowtof_chargePerCM,"simhitTOB_lowtof_chargePerCM/F");

  //sim TIB cluster High Tof
  simTreeTIBHiTof->Branch("simhitTIB_hightof_charge",&simhitTIB_hightof_charge,"simhitTIB_hightof_charge/I");
  simTreeTIBHiTof->Branch("simhitTIB_hightof_pT",&simhitTIB_hightof_pT,"simhitTIB_hightof_pT/F");
  simTreeTIBHiTof->Branch("simhitTIB_hightof_detId",&simhitTIB_hightof_detId,"simhitTIB_hightof_detId/I");
  simTreeTIBHiTof->Branch("simhitTIB_hightof_chargePerCM",&simhitTIB_hightof_chargePerCM,"simhitTIB_hightof_chargePerCM/F");
  //sim TIB cluster Low Tof
  simTreeTIBLoTof->Branch("simhitTIB_lowtof_charge",&simhitTIB_lowtof_charge,"simhitTIB_lowtof_charge/I");
  simTreeTIBLoTof->Branch("simhitTIB_lowtof_pT",&simhitTIB_lowtof_pT,"simhitTIB_lowtof_pT/F");
  simTreeTIBLoTof->Branch("simhitTIB_lowtof_detId",&simhitTIB_lowtof_detId,"simhitTIB_lowtof_detId/I");
  simTreeTIBLoTof->Branch("simhitTIB_lowtof_chargePerCM",&simhitTIB_lowtof_chargePerCM,"simhitTIB_lowtof_chargePerCM/F");

  //sim TIB cluster High Tof
  simTreeTECHiTof->Branch("simhitTEC_hightof_charge",&simhitTEC_hightof_charge,"simhitTEC_hightof_charge/I");
  simTreeTECHiTof->Branch("simhitTEC_hightof_pT",&simhitTEC_hightof_pT,"simhitTEC_hightof_pT/F");
  simTreeTECHiTof->Branch("simhitTEC_hightof_detId",&simhitTEC_hightof_detId,"simhitTEC_hightof_detId/I");
  simTreeTECHiTof->Branch("simhitTEC_hightof_chargePerCM",&simhitTEC_hightof_chargePerCM,"simhitTEC_hightof_chargePerCM/F");
  //sim TIB cluster Low Tof
  simTreeTECLoTof->Branch("simhitTEC_lowtof_charge",&simhitTEC_lowtof_charge,"simhitTEC_lowtof_charge/I");
  simTreeTECLoTof->Branch("simhitTEC_lowtof_pT",&simhitTEC_lowtof_pT,"simhitTEC_lowtof_pT/F");
  simTreeTECLoTof->Branch("simhitTEC_lowtof_detId",&simhitTEC_lowtof_detId,"simhitTEC_lowtof_detId/I");
  simTreeTECLoTof->Branch("simhitTEC_lowtof_chargePerCM",&simhitTEC_lowtof_chargePerCM,"simhitTEC_lowtof_chargePerCM/F");

  //sim TIB cluster High Tof
  simTreeTIDHiTof->Branch("simhitTID_hightof_charge",&simhitTID_hightof_charge,"simhitTID_hightof_charge/I");
  simTreeTIDHiTof->Branch("simhitTID_hightof_pT",&simhitTID_hightof_pT,"simhitTID_hightof_pT/F");
  simTreeTIDHiTof->Branch("simhitTID_hightof_detId",&simhitTID_hightof_detId,"simhitTID_hightof_detId/I");
  simTreeTIDHiTof->Branch("simhitTID_hightof_chargePerCM",&simhitTID_hightof_chargePerCM,"simhitTID_hightof_chargePerCM/F");
  //sim TIB cluster Low Tof
  simTreeTIDLoTof->Branch("simhitTID_lowtof_charge",&simhitTID_lowtof_charge,"simhitTID_lowtof_charge/I");
  simTreeTIDLoTof->Branch("simhitTID_lowtof_pT",&simhitTID_lowtof_pT,"simhitTID_lowtof_pT/F");
  simTreeTIDLoTof->Branch("simhitTID_lowtof_detId",&simhitTID_lowtof_detId,"simhitTID_lowtof_detId/I");
  simTreeTIDLoTof->Branch("simhitTID_lowtof_chargePerCM",&simhitTID_lowtof_chargePerCM,"simhitTID_lowtof_chargePerCM/F");


  // trackTree->Branch("simhit_charge",&simhit_charge,"simhit_charge/I");
  // trackTree->Branch("simhit_pT",&simhit_pT,"simhit_pT/F");
  // trackTree->Branch("simhit_firstStrip",&simhit_firstStrip,"simhit_firstStrip/I");
  // trackTree->Branch("simhit_endStrip",&simhit_endStrip,"simhit_endStrip/I");
  // trackTree->Branch("simhit_barycenter",&simhit_barycenter,"simhit_barycenter/F");
  
  simTreeTOBLoTof->Branch("simhitTOB_x",&simhitTOB_x,"simhitTOB_x/F");
  simTreeTOBLoTof->Branch("simhitTOB_y",&simhitTOB_y,"simhitTOB_y/F");
  simTreeTOBLoTof->Branch("simhitTOB_z",&simhitTOB_z,"simhitTOB_z/F");

  simTreeTIBLoTof->Branch("simhitTIB_x",&simhitTIB_x,"simhitTIB_x/F");
  simTreeTIBLoTof->Branch("simhitTIB_y",&simhitTIB_y,"simhitTIB_y/F");
  simTreeTIBLoTof->Branch("simhitTIB_z",&simhitTIB_z,"simhitTIB_z/F");

  simTreeTECLoTof->Branch("simhitTEC_x",&simhitTEC_x,"simhitTEC_x/F");
  simTreeTECLoTof->Branch("simhitTEC_y",&simhitTEC_y,"simhitTEC_y/F");
  simTreeTECLoTof->Branch("simhitTEC_z",&simhitTEC_z,"simhitTEC_z/F");

  simTreeTIDLoTof->Branch("simhitTID_x",&simhitTID_x,"simhitTID_x/F");
  simTreeTIDLoTof->Branch("simhitTID_y",&simhitTID_y,"simhitTID_y/F");
  simTreeTIDLoTof->Branch("simhitTID_z",&simhitTID_z,"simhitTID_z/F");

  simTreeTOBLoTof->Branch("fromTrack",&fromTrack,"fromTrack/O");
}

TrackAnalyzer::~TrackAnalyzer() {

}

// ------------ method called for each event  ------------
void TrackAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  //SetVarToZero();

  runN = iEvent.id().run();
  eventN = iEvent.id().event();
  lumi = iEvent.id().luminosityBlock();
  const auto& tracksHandle = iEvent.getHandle(tracksToken_);
  const auto& jetHandle = iEvent.getHandle(jetToken_);

  
  const TrackerGeometry* trackerGeometry = &iSetup.getData(m_geomToken_);
  if (!tracksHandle.isValid()) {
    edm::LogError("TrackAnalyzer") << "No valid track collection found";
    return;
  }
  if (!jetHandle.isValid()) {
    edm::LogError("TrackAnalyzer") << "No valid jet collection found";
    return;
  }

  if(isMC_){
    edm::Handle<std::vector<PSimHit>> g4SimHitsTOB_htof = iEvent.getHandle(g4SimHitsHighTofTokenTOB_);
    edm::Handle<std::vector<PSimHit>> g4SimHitsTOB_ltof = iEvent.getHandle(g4SimHitsLowTofTokenTOB_);

    edm::Handle<std::vector<PSimHit>> g4SimHitsTIB_htof = iEvent.getHandle(g4SimHitsHighTofTokenTIB_);
    edm::Handle<std::vector<PSimHit>> g4SimHitsTIB_ltof = iEvent.getHandle(g4SimHitsLowTofTokenTIB_);

    edm::Handle<std::vector<PSimHit>> g4SimHitsTEC_htof = iEvent.getHandle(g4SimHitsHighTofTokenTEC_);
    edm::Handle<std::vector<PSimHit>> g4SimHitsTEC_ltof = iEvent.getHandle(g4SimHitsLowTofTokenTEC_);

    edm::Handle<std::vector<PSimHit>> g4SimHitsTID_htof = iEvent.getHandle(g4SimHitsHighTofTokenTID_);
    edm::Handle<std::vector<PSimHit>> g4SimHitsTID_ltof = iEvent.getHandle(g4SimHitsLowTofTokenTID_);

    
    edm::Handle<std::vector<SimTrack>> g4SimHitsTrack = iEvent.getHandle(g4SimHitsTrackToken_);

    
    if (!g4SimHitsTOB_ltof.isValid()){
      edm::LogError("TrackAnalyzer") << "No valid g4SimHits collection found";
      return;
    }
  

    fromTrack = false;
    for (const auto& simhit : *g4SimHitsTOB_htof){
      auto momentum = simhit.momentumAtEntry(); // momentumEntry is a LorentzVector
      
      float pT = std::sqrt(momentum.x()*momentum.x() + momentum.y()*momentum.y()); 
      simhitTOB_hightof_pT =  pT;
      simhitTOB_hightof_detId = simhit.detUnitId();
      simhitTOB_hightof_charge = calculateDepositedCharge(simhit);
      simhitTOB_hightof_chargePerCM = simhitTOB_hightof_charge*siStripClusterTools::sensorThicknessInverse(simhitTOB_hightof_detId);
    }

    fromTrack = false;
    for (const auto& simhit : *g4SimHitsTOB_ltof){
      auto momentum = simhit.momentumAtEntry(); // momentumEntry is a LorentzVector
      
      float pT = std::sqrt(momentum.x()*momentum.x() + momentum.y()*momentum.y()); 
      simhitTOB_lowtof_pT = pT;
      simhitTOB_lowtof_detId = simhit.detUnitId();
      simhitTOB_lowtof_charge = calculateDepositedCharge(simhit);
      simhitTOB_lowtof_chargePerCM = simhitTOB_lowtof_charge*siStripClusterTools::sensorThicknessInverse(simhitTOB_lowtof_detId);

      
      simhitTOB_x = simhit.localDirection().x();
      simhitTOB_y = simhit.localDirection().y();
      simhitTOB_z = simhit.localDirection().z();
      if(simhit.trackId() > 0)
    	fromTrack = true;
      simTreeTOBLoTof->Fill();
    }

    if (!g4SimHitsTIB_ltof.isValid()){
      edm::LogError("TrackAnalyzer") << "No valid g4SimHits collection found";
      return;
    }
  

    fromTrack = false;
    for (const auto& simhit : *g4SimHitsTIB_htof){
      auto momentum = simhit.momentumAtEntry(); // momentumEntry is a LorentzVector
      
      float pT = std::sqrt(momentum.x()*momentum.x() + momentum.y()*momentum.y()); 
      simhitTIB_hightof_pT =  pT;
      simhitTIB_hightof_detId = simhit.detUnitId();
      simhitTIB_hightof_charge = calculateDepositedCharge(simhit);
      simhitTIB_hightof_chargePerCM = simhitTIB_hightof_charge*siStripClusterTools::sensorThicknessInverse(simhitTIB_hightof_detId);
    }

    fromTrack = false;
    for (const auto& simhit : *g4SimHitsTIB_ltof){
      auto momentum = simhit.momentumAtEntry(); // momentumEntry is a LorentzVector
      
      float pT = std::sqrt(momentum.x()*momentum.x() + momentum.y()*momentum.y()); 
      simhitTIB_lowtof_pT = pT;
      simhitTIB_lowtof_detId = simhit.detUnitId();
      simhitTIB_lowtof_charge = calculateDepositedCharge(simhit);
      simhitTIB_lowtof_chargePerCM = simhitTIB_lowtof_charge*siStripClusterTools::sensorThicknessInverse(simhitTIB_lowtof_detId);

      simhitTIB_x = simhit.localDirection().x();
      simhitTIB_y = simhit.localDirection().y();
      simhitTIB_z = simhit.localDirection().z();
      if(simhit.trackId() > 0)
    	fromTrack = true;
      simTreeTIBLoTof->Fill();
    }

    if (!g4SimHitsTEC_ltof.isValid()){
      edm::LogError("TrackAnalyzer") << "No valid g4SimHits collection found";
      return;
    }
  

    fromTrack = false;
    for (const auto& simhit : *g4SimHitsTEC_htof){
      auto momentum = simhit.momentumAtEntry(); // momentumEntry is a LorentzVector
      
      float pT = std::sqrt(momentum.x()*momentum.x() + momentum.y()*momentum.y()); 
      simhitTEC_hightof_pT =  pT;
      simhitTEC_hightof_detId = simhit.detUnitId();
      simhitTEC_hightof_charge = calculateDepositedCharge(simhit);
      simhitTEC_hightof_chargePerCM = simhitTEC_hightof_charge*siStripClusterTools::sensorThicknessInverse(simhitTEC_hightof_detId);
    }

    fromTrack = false;
    for (const auto& simhit : *g4SimHitsTEC_ltof){
      auto momentum = simhit.momentumAtEntry(); // momentumEntry is a LorentzVector
      
      float pT = std::sqrt(momentum.x()*momentum.x() + momentum.y()*momentum.y()); 
      simhitTEC_lowtof_pT = pT;
      simhitTEC_lowtof_detId = simhit.detUnitId();
      simhitTEC_lowtof_charge = calculateDepositedCharge(simhit);
      simhitTEC_lowtof_chargePerCM = simhitTEC_lowtof_charge*siStripClusterTools::sensorThicknessInverse(simhitTEC_lowtof_detId);

      simhitTEC_x = simhit.localDirection().x();
      simhitTEC_y = simhit.localDirection().y();
      simhitTEC_z = simhit.localDirection().z();
      if(simhit.trackId() > 0)
    	fromTrack = true;
      simTreeTECLoTof->Fill();
    }

    if (!g4SimHitsTID_ltof.isValid()){
      edm::LogError("TrackAnalyzer") << "No valid g4SimHits collection found";
      return;
    }
  

    fromTrack = false;
    for (const auto& simhit : *g4SimHitsTID_htof){
      auto momentum = simhit.momentumAtEntry(); // momentumEntry is a LorentzVector
      
      float pT = std::sqrt(momentum.x()*momentum.x() + momentum.y()*momentum.y()); 
      simhitTID_hightof_pT =  pT;
      simhitTID_hightof_detId = simhit.detUnitId();
      simhitTID_hightof_charge = calculateDepositedCharge(simhit);
      simhitTID_hightof_chargePerCM = simhitTID_hightof_charge*siStripClusterTools::sensorThicknessInverse(simhitTID_hightof_detId);
    }

    fromTrack = false;
    for (const auto& simhit : *g4SimHitsTID_ltof){
      auto momentum = simhit.momentumAtEntry(); // momentumEntry is a LorentzVector
      
      float pT = std::sqrt(momentum.x()*momentum.x() + momentum.y()*momentum.y()); 
      simhitTID_lowtof_pT = pT;
      simhitTID_lowtof_detId = simhit.detUnitId();
      simhitTID_lowtof_charge = calculateDepositedCharge(simhit);
      simhitTID_lowtof_chargePerCM = simhitTID_lowtof_charge*siStripClusterTools::sensorThicknessInverse(simhitTID_lowtof_detId);

      simhitTID_x = simhit.localDirection().x();
      simhitTID_y = simhit.localDirection().y();
      simhitTID_z = simhit.localDirection().z();
      if(simhit.trackId() > 0)
    	fromTrack = true;
      simTreeTIDLoTof->Fill();
    }

    
    
  }

  // Retrieve the actual product from the handle
  const reco::TrackCollection& tracks = *tracksHandle;
  const reco::PFJetCollection& jets = *jetHandle;
  

  for (const auto& track : tracks) {
    
    if (nTracks >= nMax) {
      continue;
    }
    trkPt[nTracks]  = track.pt();
    trkEta[nTracks] = track.eta();
    trkPhi[nTracks] = track.phi();
    trkDxy1[nTracks] = track.dxy();
    trkDz1[nTracks] = track.dz();
    trkDxyError1[nTracks]= track.dxyError();
    trkDzError1[nTracks] = track.dzError();
    trkPtError[nTracks] = track.ptError();
    trkChi2[nTracks] = track.chi2();
    trkNChi2[nTracks] = track.normalizedChi2();
    trkNdof[nTracks] = track.ndof();
    trkAlgo[nTracks] = track.algo();
    trkNHit[nTracks] = track.numberOfValidHits();
    trkNlayer[nTracks] = track.hitPattern().trackerLayersWithMeasurement();
    nTracks++;

    
    for (auto const &hit : track.recHits()){
      if (!hit->isValid())
      	continue;

      DetId detid = hit->geographicalId();
      int subDet = detid.subdetId();
      Point3DBase<float, GlobalTag> modulepos = trackerGeometry->idToDet(detid)->position();
     
      bool hitInStrip = (subDet == SiStripDetId::TIB) || (subDet == SiStripDetId::TID) || (subDet == SiStripDetId::TOB) || (subDet == SiStripDetId::TEC);
      if(subDet == SiStripDetId::TIB)
	subdet = 1;
      else if (subDet == SiStripDetId::TOB)
	subdet = 2;
      else if (subDet == SiStripDetId::TID)
	subdet = 3;
      else if (subDet == SiStripDetId::TEC)
	subdet = 4;
      else
	subdet = 0;
      if (hitInStrip){

  	clus_detId  = detid;
  	hit_eta = modulepos.eta();
  	hit_phi = modulepos.phi();
  	hit_pt = track.pt();
  	const std::type_info &type = typeid(*hit);
  	if (type == typeid(SiStripRecHit1D)){
	  
  	  const SiStripRecHit1D *striphit = dynamic_cast<const SiStripRecHit1D *>(hit);
  	  if(striphit != nullptr){
  	    SiStripRecHit1D::ClusterRef stripclust(striphit->cluster());
  	    clus_chargePerCM = siStripClusterTools::chargePerCM(detid,(*stripclust));
  	    clus_charge = stripclust->charge();
  	    clus_width  = stripclust->size();
  	    clus_avgcharge = (clus_charge+clus_width/2)/clus_width;
  	    clus_barycenter = stripclust->barycenter();
  	    clus_firstStrip = stripclust->firstStrip();
  	    clus_endStrip  = stripclust->endStrip();
  	    trackTree->Fill();
	   
  	  }
	  
  	}
  	else if (type == typeid(SiStripRecHit2D)){
 	  
  	  const SiStripRecHit2D *striphit2D = dynamic_cast<const SiStripRecHit2D *>(hit);
  	  if(striphit2D != nullptr){    

  	    SiStripRecHit2D::ClusterRef stripclust2D(striphit2D->cluster());
  	    clus_chargePerCM = siStripClusterTools::chargePerCM(detid,(*stripclust2D));
  	    clus_charge = stripclust2D->charge();   
  	    clus_width  = stripclust2D->size();
  	    clus_avgcharge = (clus_charge+clus_width/2)/clus_width;
  	    clus_barycenter = stripclust2D->barycenter();
  	    clus_firstStrip = stripclust2D->firstStrip();
  	    clus_endStrip  = stripclust2D->endStrip();
  	    trackTree->Fill();
	    
  	  }
  	}
      }
    }
    
  }



  // for (const auto& jet : jets) {
  //   if (nJets >= nMax) {
  //     continue;
  //   }
  //   jetPt[nJets] = jet.pt();
  //   jetEta[nJets] = jet.eta();
  //   jetPhi[nJets] = jet.phi();
  //   jetMass[nJets] = jet.mass();
  //   nJets++;
  // }
  // trackTree->Fill();

  
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  // if the SetupData is always needed
  auto setup = iSetup.getData(setupToken_);
  // if need the ESHandle to check if the SetupData was there or not
  auto pSetup = iSetup.getHandle(setupToken_);
#endif
}


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void TrackAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  // edm::ParameterSetDescription desc;
  // desc.addUntracked<edm::InputTag>("tracks", edm::InputTag("generalTracks"));
  // descriptions.addWithDefaultLabel(desc);
  
}
double TrackAnalyzer::calculateDepositedCharge(const PSimHit& hit) {
  const double eV_per_pair = 3.6; // Energy needed to create an e-h pair in silicon (eV)
  const double GeV_to_eV = 1e9;  // Conversion factor from GeV to eV

  // Get energy loss in GeV
  double energyLoss_GeV = hit.energyLoss();

  // Convert to eV
  double energyLoss_eV = energyLoss_GeV * GeV_to_eV;

  // Calculate charge in number of e-h pairs
  double charge = energyLoss_eV / eV_per_pair; // Charge in number of e-h pairs

  return charge/(270); // Average no. of electrons per ADC https://cds.cern.ch/record/1379833/files/10.1016_j.phpro.2012.02.423.pdf
}
void TrackAnalyzer::SetVarToZero()
{
  nTracks = 0;
  nJets = 0;
  simhitTOB_charge =0;
  simhitTOB_lowtof_charge =0;
  simhitTOB_hightof_charge =0;
  simhitTOB_lowtof_pT =0;
  simhitTOB_hightof_pT =0;

  simhitTIB_charge =0;
  simhitTIB_lowtof_charge =0;
  simhitTIB_hightof_charge =0;
  simhitTIB_lowtof_pT =0;
  simhitTIB_hightof_pT =0;

  simhitTEC_charge =0;
  simhitTEC_lowtof_charge =0;
  simhitTEC_hightof_charge =0;
  simhitTEC_lowtof_pT =0;
  simhitTEC_hightof_pT =0;

  simhitTID_charge =0;
  simhitTID_lowtof_charge =0;
  simhitTID_hightof_charge =0;
  simhitTID_lowtof_pT =0;
  simhitTID_hightof_pT =0;

  clus_charge = 0;
  clus_width = 0;
  clus_avgcharge = 0;
  clus_barycenter =0;
  clus_firstStrip =0;
  clus_endStrip = 0;
  clus_detId = 0;
  hit_eta = 0;
  hit_phi = 0;
  for (int i = 0; i < nMax; ++i) {
    trkPt[i] = 0;
    trkEta[i] = 0;
    trkPhi[i] = 0;
    trkDxy1[i] = 0;
    trkDz1[i] = 0;
    trkDxyError1[i] = 0;
    trkDzError1[i] = 0;
    trkPtError[i] = 0;
    trkChi2[i] = 0;
    trkNChi2[i] = 0;
    trkNdof[i]=0;
  }
  for (int i=0; i < nMax; ++i){
    jetPt[i] = 0;
    jetEta[i] = 0;
    jetPhi[i] = 0;
    jetMass[i] = 0;
    
  }
}
//define this as a plug-in
DEFINE_FWK_MODULE(TrackAnalyzer);
