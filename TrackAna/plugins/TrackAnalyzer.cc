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
  //edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
  
  TTree* trackTree;
  edm::Service<TFileService> fs;

  edm::EventNumber_t eventN;
  const static int nMax = 40000000;

  uint16_t nTracks;
  
  int runN;
  int lumi;

  float trkPt[nMax];
  float trkEta[nMax];
  float trkPhi[nMax];
  float trkDxy1[nMax];
  float trkDxyError1[nMax];
  float trkDz1[nMax];
  float trkDzError1[nMax];
  unsigned char trkAlgo[nMax];
  uint8_t trkNHit[nMax];
  uint8_t trkNdof[nMax];
  uint8_t trkNlayer[nMax];

  float trkChi2[nMax];
  float trkPtError[nMax];
  
};


TrackAnalyzer::TrackAnalyzer(const edm::ParameterSet& iConfig){
//: tracksToken_(consumes<TrackCollection>(iConfig.getUntrackedParameter<edm::InputTag>("tracks"))) {
  
  tracksToken_ = consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("tracks"));
  //beamSpotToken_ = consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot"));
    
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
  trackTree->Branch("trkPtError",  trkPtError, "trkPtError[nTracks]/F");

  trackTree->Branch("trkAlgo",trkAlgo,"trkAlgo[nTracks]/s");
  trackTree->Branch("trkNHit",trkNHit,"trkNHit[nTracks]/I");
  trackTree->Branch("trkNdof",trkNdof,"trkNdof[nTracks]/I");
  trackTree->Branch("trkNlayer",trkNlayer,"trkNlayer[nTracks]/I");
}

TrackAnalyzer::~TrackAnalyzer() {

}

// ------------ method called for each event  ------------
void TrackAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  SetVarToZero();

  runN = iEvent.id().run();
  eventN = iEvent.id().event();
  //lumi = iEvent.id().lumi();
  const auto& tracksHandle = iEvent.getHandle(tracksToken_);

  
  // edm::Handle<reco::TrackCollection> tracks;
  // iEvent.getByToken(tracksToken_, tracks);

  // edm::Handle<reco::BeamSpot> beamSpot;
  // iEvent.getByToken(beamSpotToken_, beamSpot);

  if (!tracksHandle.isValid()) {
    edm::LogError("TrackAnalyzer") << "No valid track collection found";
    return;
  }
  
  // Retrieve the actual product from the handle
  const reco::TrackCollection& tracks = *tracksHandle;

  // if (!beamSpot.isValid()) {
  //   edm::LogError("TrackAnalyzer") << "No valid beam spot found";
  //   return;
  // }
 
  
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
    trkChi2[nTracks] = track.normalizedChi2();

    trkAlgo[nTracks] = track.algo();
    trkNHit[nTracks] = track.numberOfValidHits();
    trkNlayer[nTracks] = track.hitPattern().trackerLayersWithMeasurement();
    nTracks++;
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
void TrackAnalyzer::SetVarToZero()
{
  nTracks = 0;
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
    
  }
}
//define this as a plug-in
DEFINE_FWK_MODULE(TrackAnalyzer);
