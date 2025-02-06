// -*- C++ -*-
//
// Package:    PREanalysis/Analyzer
// Class:      RawAnalyzer
//
/**\class RawAnalyzer RawAnalyzer.cc PREanalysis/Analyzer/plugins/RawAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:   Vinaya Krishnan Muraleedharan Nair Bindhu
//         Created:  Tue, 30 Jul 2024 08:35:07 GMT
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

// user include files
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"
#include "DataFormats/Common/interface/DetSet.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/Common/interface/DetSetVectorNew.h"
#include "DataFormats/SiStripCluster/interface/SiStripApproximateCluster.h"
#include "DataFormats/SiStripCluster/interface/SiStripApproximateClusterCollection.h"
#include "DataFormats/SiStripCluster/interface/SiStripCluster.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/SiStripCluster/interface/SiStripClusterTools.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/SiStripDetId/interface/SiStripDetId.h"

#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/StripGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/CommonTopologies/interface/StripTopology.h"

//ROOT inclusion
#include "TROOT.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TTree.h"
#include "TMath.h"
#include "TList.h"
#include "TString.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

using reco::TrackCollection;

class RawAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit RawAnalyzer(const edm::ParameterSet&);
  ~RawAnalyzer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  
  edm::InputTag inputTagClusters;
  
  // Event Data
  edm::EDGetTokenT<edmNew::DetSetVector<SiStripCluster>> clusterToken;
  
  // Event Setup Data
  edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> tkGeomToken_;

  TTree* offlineClusterTree;
  edm::Service<TFileService> fs;

  edm::EventNumber_t eventN;
  int runN;
  int lumi;

  // for offlineClusterTree
  uint32_t    detId;
  uint16_t    firstStrip;
  uint16_t    endStrip;
  float       barycenter;
  uint16_t    size;
  int         charge;
  float       chargePerCM;
  int         subdet;
  const static int nMax = 800000;
  float       hitX[nMax];
  float       hitY[nMax];
  uint16_t    channel[nMax];
  uint16_t    adc[nMax];

  
  // ----------member data ---------------------------
  //edm::EDGetTokenT<TrackCollection> tracksToken_;  //used to select what tracks to read from configuration file
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  edm::ESGetToken<SetupData, SetupRecord> setupToken_;
#endif
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
RawAnalyzer::RawAnalyzer(const edm::ParameterSet& iConfig){
//  : tracksToken_(consumes<TrackCollection>(iConfig.getUntrackedParameter<edm::InputTag>("tracks"))) {

  //tracksToken_ = consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("tracks"));
  inputTagClusters  = iConfig.getParameter<edm::InputTag>("siStripClustersTag");
  clusterToken  = consumes<edmNew::DetSetVector<SiStripCluster>>(inputTagClusters);

  tkGeomToken_ = esConsumes();

  usesResource("TFileService");

  offlineClusterTree = fs->make<TTree>("offlineClusterTree", "offlineClusterTree");
  offlineClusterTree->Branch("event", &eventN, "event/i");
  offlineClusterTree->Branch("run",   &runN, "run/I");
  offlineClusterTree->Branch("lumi",  &lumi, "lumi/I");

  offlineClusterTree->Branch("detId", &detId, "detId/i");
  offlineClusterTree->Branch("firstStrip", &firstStrip, "firstStrip/s");
  offlineClusterTree->Branch("endStrip", &endStrip, "endStrip/s");
  offlineClusterTree->Branch("barycenter", &barycenter, "barycenter/F");
  offlineClusterTree->Branch("size", &size, "size/s");
  offlineClusterTree->Branch("charge", &charge, "charge/I");
  offlineClusterTree->Branch("chargePerCM", &chargePerCM, "chargePerCM/F");
  offlineClusterTree->Branch("subdet", &subdet, "subdet/I");

  offlineClusterTree->Branch("x", hitX, "x[size]/F");
  offlineClusterTree->Branch("y", hitY, "y[size]/F");
  offlineClusterTree->Branch("channel", channel, "channel[size]/s");
  offlineClusterTree->Branch("adc", adc, "adc[size]/s");

  

// #ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
//   setupDataToken_ = esConsumes<SetupData, SetupRecord>();
// #endif
  //now do what ever initialization is needed
}

RawAnalyzer::~RawAnalyzer() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty
}

//
// member functions
//

// ------------ method called for each event  ------------
void RawAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  edm::Handle<edmNew::DetSetVector<SiStripCluster>> clusterCollection = iEvent.getHandle(clusterToken);
  const auto& tkGeom = &iSetup.getData(tkGeomToken_);
  const auto tkDets = tkGeom->dets();

  for (const auto& detSiStripClusters : *clusterCollection) {
    eventN = iEvent.id().event();
    runN   = (int) iEvent.id().run();
    lumi   = (int) iEvent.id().luminosityBlock();
    detId = detSiStripClusters.detId();
    DetId detector_id =  detSiStripClusters.detId();
    int subDet = detector_id.subdetId();
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
    //if (subDet == SiStripDetId::TOB){
      for (const auto& stripCluster : detSiStripClusters) {
      
	firstStrip  = stripCluster.firstStrip();
	endStrip    = stripCluster.endStrip();
	barycenter  = stripCluster.barycenter();
	size        = stripCluster.size();
	charge      = stripCluster.charge();
	chargePerCM = siStripClusterTools::chargePerCM(detId,stripCluster);
	const auto& _detId = detId; // for the capture clause in the lambda function
	auto det = std::find_if(tkDets.begin(), tkDets.end(), [_detId](auto& elem) -> bool {
	    return (elem->geographicalId().rawId() == _detId);
	  });
	const StripTopology& p = dynamic_cast<const StripGeomDetUnit*>(*det)->specificTopology();
	for (int strip = firstStrip; strip < endStrip+1; ++strip)
	  {
	    GlobalPoint gp = (tkGeom->idToDet(detId))->surface().toGlobal(p.localPosition((float) strip));

	    hitX   [strip - firstStrip] = gp.x();
	    hitY   [strip - firstStrip] = gp.y();
	    channel[strip - firstStrip] = strip;
	    adc    [strip - firstStrip] = stripCluster[strip - firstStrip];
	  }

	offlineClusterTree->Fill();
      }
      
      //}
      //}
  }
  // for (const auto& track : iEvent.get(tracksToken_)) {
  //   // do something with track parameters, e.g, plot the charge.
  //   // int charge = track.charge();
  // }
  
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  // if the SetupData is always needed
  auto setup = iSetup.getData(setupToken_);
  // if need the ESHandle to check if the SetupData was there or not
  auto pSetup = iSetup.getHandle(setupToken_);
#endif
}


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void RawAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {


  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("siStripClustersTag", edm::InputTag("siStripClusters"));
  descriptions.add("RawAnalyzer", desc);
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly wh// at you do use, even if it is no parameters
  // edm::ParameterSetDescription desc;
  // desc.setUnknown();
  // descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //edm::ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks", edm::InputTag("ctfWithMaterialTracks"));
  //descriptions.addWithDefaultLabel(desc);
}

#include "FWCore/Framework/interface/MakerMacros.h"
//define this as a plug-in
DEFINE_FWK_MODULE(RawAnalyzer);
