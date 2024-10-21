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
#include "DataFormats/TrackerRecHit2D/interface/TrackerSingleRecHit.h"

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
  // ----------member data ---------------------------
  edm::EDGetTokenT<reco::TrackCollection> tracksToken_;  //used to select what tracks to read from configuration file
  edm::EDGetTokenT<reco::PFJetCollection> jetToken_;


  //edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
  
  TTree* trackTree;

  edm::Service<TFileService> fs;

  edm::EventNumber_t eventN;
  const static int nMax = 40000000;

  uint16_t nTracks;
  int nJets;

  int runN;
  int lumi;

  int charge_lowpt;
  int charge_highpt;
  int avgcharge_lowpt;
  int avgcharge_highpt;
  int width_lowpt;
  int width_highpt;
  
  
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
  tracksToken_ = consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("tracks"));
  jetToken_   = consumes<reco::PFJetCollection>(iConfig.getParameter<edm::InputTag>("jets"));
    
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
  trackTree->Branch("charge_lowpt",&charge_lowpt,"charge_lowpt/I");
  trackTree->Branch("charge_highpt",&charge_highpt,"charge_highpt/I");
  trackTree->Branch("width_lowpt",&width_lowpt,"width_lowpt/I");
  trackTree->Branch("width_highpt",&width_highpt,"width_highpt/I");
  trackTree->Branch("avgcharge_lowpt",&avgcharge_lowpt,"avgcharge_lowpt/I");
  trackTree->Branch("avgcharge_highpt",&avgcharge_highpt,"avgcharge_highpt/I");
}

TrackAnalyzer::~TrackAnalyzer() {

}

// ------------ method called for each event  ------------
void TrackAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  SetVarToZero();

  runN = iEvent.id().run();
  eventN = iEvent.id().event();
  lumi = iEvent.id().luminosityBlock();
  const auto& tracksHandle = iEvent.getHandle(tracksToken_);
  const auto& jetHandle = iEvent.getHandle(jetToken_);
  
  if (!tracksHandle.isValid()) {
    edm::LogError("TrackAnalyzer") << "No valid track collection found";
    return;
  }
  if (!jetHandle.isValid()) {
    edm::LogError("TrackAnalyzer") << "No valid jet collection found";
    return;
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
      
      bool hitInStrip = (subDet == SiStripDetId::TIB) || (subDet == SiStripDetId::TID) || (subDet == SiStripDetId::TOB) || (subDet == SiStripDetId::TEC);
      
      if (hitInStrip){
      	const std::type_info &type = typeid(*hit);
      	if (type == typeid(SiStripRecHit1D)){
	  
      	  const SiStripRecHit1D *striphit = dynamic_cast<const SiStripRecHit1D *>(hit);
	  if(striphit != nullptr){    
	    SiStripRecHit1D::ClusterRef stripclust(striphit->cluster());
	    if (track.pt() < 0.72 && track.chi2()/track.ndof() < 3){
	      charge_lowpt = stripclust->charge();
	      width_lowpt  = stripclust->size();
	      avgcharge_lowpt = (charge_lowpt+width_lowpt/2)/width_lowpt;
	      
	    }
	    if (track.pt() > 0.72 && track.chi2()/track.ndof() < 3){
	      charge_highpt = stripclust->charge();
	      width_highpt  = stripclust->size();
	      avgcharge_highpt = (charge_highpt+width_highpt/2)/width_highpt;
	    }
      	  }
	}
	if (type == typeid(SiStripRecHit2D)){
	  
      	  const SiStripRecHit2D *striphit2D = dynamic_cast<const SiStripRecHit2D *>(hit);
	  if(striphit2D != nullptr){    
	    SiStripRecHit2D::ClusterRef stripclust2D(striphit2D->cluster());
	    if (track.pt() < 0.72 && track.chi2()/track.ndof() < 3){
	      charge_lowpt = stripclust2D->charge();   
	      width_lowpt  = stripclust2D->size();
	      avgcharge_lowpt = (charge_lowpt+width_lowpt/2)/width_lowpt;
	    }
	    if (track.pt() > 0.72 && track.chi2()/track.ndof() < 3){
	      charge_highpt = stripclust2D->charge();
	      width_highpt  = stripclust2D->size();
	      avgcharge_highpt = (charge_highpt+width_highpt/2)/width_highpt;
	    }
      	  }
      	}
      }
      //++charge;
      trackTree->Fill();
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
  //trackTree->Fill();

  
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
void TrackAnalyzer::SetVarToZero()
{
  nTracks = 0;
  nJets = 0;
  charge_lowpt = 0;
  charge_highpt = 0;
  width_lowpt = 0;
  width_highpt = 0;
  avgcharge_lowpt = 0;
  avgcharge_highpt = 0;
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
