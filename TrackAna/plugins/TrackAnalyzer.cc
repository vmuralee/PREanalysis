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
  edm::EDGetTokenT<std::vector<PSimHit>> g4SimHitsHighTofToken_;
  edm::EDGetTokenT<std::vector<PSimHit>> g4SimHitsLowTofToken_;
  edm::EDGetTokenT<std::vector<SimTrack>> g4SimHitsTrackToken_;

  
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
  float clus_barycenter;

  int simhit_hightof_charge;
  int simhit_lowtof_charge;
  int simhit_charge;
  float simhit_chargePerCM;
  int simhit_width;
  int simhit_avgcharge;
  int simhit_hightof_detId;
  int simhit_lowtof_detId;
  float simhit_hightof_pT;
  float simhit_lowtof_pT;
  float simhit_pT;
  uint16_t simhit_firstStrip;
  uint16_t simhit_endStrip;
  float simhit_barycenter;

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

  
};


TrackAnalyzer::TrackAnalyzer(const edm::ParameterSet& iConfig){

  isMC_ = iConfig.getUntrackedParameter<bool>("IsMC", false);
  tracksToken_ = consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("tracks"));
  jetToken_   = consumes<reco::PFJetCollection>(iConfig.getParameter<edm::InputTag>("jets"));
  if(isMC_){
    g4SimHitsLowTofToken_ = consumes<std::vector<PSimHit>>(iConfig.getParameter<edm::InputTag>("g4SimHitsLowTofTag"));
    g4SimHitsHighTofToken_ = consumes<std::vector<PSimHit>>(iConfig.getParameter<edm::InputTag>("g4SimHitsHighTofTag"));
    g4SimHitsTrackToken_ = consumes<std::vector<SimTrack>>(iConfig.getParameter<edm::InputTag>("g4SimHitsTrackTag"));   
  }
  
  m_geomToken_ = esConsumes(); 
  usesResource("TFileService");

  trackTree = fs->make<TTree>("trackTree","trackTree");
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

  trackTree->Branch("simhit_hightof_charge",&simhit_hightof_charge,"simhit_hightof_charge/I");
  trackTree->Branch("simhit_lowtof_charge",&simhit_lowtof_charge,"simhit_lowtof_charge/I");
  trackTree->Branch("simhit_charge",&simhit_charge,"simhit_charge/I");
  trackTree->Branch("simhit_chargePerCM",&simhit_chargePerCM,"simhit_chargePerCM/F");

  trackTree->Branch("simhit_pT",&simhit_pT,"simhit_pT/F");
  trackTree->Branch("simhit_hightof_pT",&simhit_hightof_pT,"simhit_hightof_pT/F");
  trackTree->Branch("simhit_lowtof_pT",&simhit_lowtof_pT,"simhit_lowtof_pT/F");
  trackTree->Branch("simhit_hightof_detId",&simhit_hightof_detId,"simhit_hightof_detId/I");
  trackTree->Branch("simhit_lowtof_detId",&simhit_lowtof_detId,"simhit_lowtof_detId/I");
  trackTree->Branch("simhit_firstStrip",&simhit_firstStrip,"simhit_firstStrip/I");
  trackTree->Branch("simhit_endStrip",&simhit_endStrip,"simhit_endStrip/I");
  trackTree->Branch("simhit_barycenter",&simhit_barycenter,"simhit_barycenter/F");
  
  
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
    edm::Handle<std::vector<PSimHit>> g4SimHits_htof = iEvent.getHandle(g4SimHitsHighTofToken_);
    edm::Handle<std::vector<PSimHit>> g4SimHits_ltof = iEvent.getHandle(g4SimHitsLowTofToken_);
    edm::Handle<std::vector<SimTrack>> g4SimHitsTrack = iEvent.getHandle(g4SimHitsTrackToken_);

    
    if (!g4SimHits_ltof.isValid()){
      edm::LogError("TrackAnalyzer") << "No valid g4SimHits collection found";
      return;
    }
  
    for (const auto& simhit : *g4SimHits_htof){
      auto momentum = simhit.momentumAtEntry(); // momentumEntry is a LorentzVector

      float pT = std::sqrt(momentum.x()*momentum.x() + momentum.y()*momentum.y()); 
      simhit_hightof_pT =  pT;
      simhit_hightof_detId = simhit.detUnitId();
      simhit_hightof_charge = calculateDepositedCharge(simhit);
      simhit_charge = simhit_hightof_charge;
      simhit_chargePerCM = simhit_hightof_charge*siStripClusterTools::sensorThicknessInverse(simhit_hightof_detId);

      trackTree->Fill();
    }
    for (const auto& simhit : *g4SimHits_ltof){
      auto momentum = simhit.momentumAtEntry(); // momentumEntry is a LorentzVector

      float pT = std::sqrt(momentum.x()*momentum.x() + momentum.y()*momentum.y()); 
      simhit_lowtof_pT = pT;
      simhit_lowtof_detId = simhit.detUnitId();
      simhit_lowtof_charge = calculateDepositedCharge(simhit);
      simhit_charge = simhit_lowtof_charge;
      simhit_chargePerCM = simhit_lowtof_charge*siStripClusterTools::sensorThicknessInverse(simhit_lowtof_detId);

      trackTree->Fill();
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
      //bool hitInStrip = (subDet == SiStripDetId::TOB) ;
      if (!hitInStrip)continue;

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



  for (const auto& jet : jets) {
    if (nJets >= nMax) {
      continue;
    }
    jetPt[nJets] = jet.pt();
    jetEta[nJets] = jet.eta();
    jetPhi[nJets] = jet.phi();
    jetMass[nJets] = jet.mass();
    nJets++;
  }
  trackTree->Fill();

  
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
  simhit_charge =0;
  simhit_lowtof_charge =0;
  simhit_hightof_charge =0;
  simhit_lowtof_pT =0;
  simhit_hightof_pT =0;
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
