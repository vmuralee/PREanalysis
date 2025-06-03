// -*- C++ -*-
//
// Package:    PREanalysis/StripClusterAnalysis
// Class:      StripClusterAnalysis
//
/**\class StripClusterAnalysis StripClusterAnalysis.cc PREanalysis/StripClusterAnalysis/plugins/StripClusterAnalysis.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Vinaya Krishnan Muraleedharan Nair Bindhu
//         Created:  Mon, 16 Dec 2024 11:10:32 GMT
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

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"

#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/CommonTopologies/interface/StripTopology.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/StripGeomDetUnit.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/SiStripCluster/interface/SiStripCluster.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/Common/interface/DetSet.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/Common/interface/DetSetVectorNew.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit1D.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2D.h"
#include "DataFormats/SiStripCluster/interface/SiStripClusterTools.h"
#include "DataFormats/TrackerRecHit2D/interface/TrackerSingleRecHit.h"

#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/Track/interface/SimTrack.h"

#include "SimTracker/TrackerHitAssociation/interface/TrackerHitAssociator.h"

#include "FWCore/ServiceRegistry/interface/Service.h"

#include "TH1.h"

//ROOT inclusion
#include "TROOT.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TTree.h"
#include "TMath.h"
#include "TList.h"
#include "TString.h"
#include "TVector3.h"
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.



class StripClusterAnalysis : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit StripClusterAnalysis(const edm::ParameterSet&);
  ~StripClusterAnalysis() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;
  GlobalPoint HitPosition(const PSimHit& hit, const DetId detId,const  TrackerGeometry* tkGeom);
  // ----------member data ---------------------------

  // simhits collection
  edm::EDGetTokenT<std::vector<PSimHit>> g4SimHitsCollectionToken_;  // TOB

  edm::EDGetTokenT<reco::TrackCollection> tracksToken_;
  //simTracks
  edm::EDGetTokenT<std::vector<SimTrack>> g4SimHitsTrackToken_;  
  //siStrip cluster
  edm::EDGetTokenT<edmNew::DetSetVector<SiStripCluster>> clusterToken_;

  // Event Setup Data
  edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> tkGeomToken_;
  
  TrackerHitAssociator::Config trackerHitAssociatorConfig_;

  TTree* clusterTree;
  TTree* recohitTree;

  edm::Service<TFileService> fs;
  
  edm::EventNumber_t eventN;
  int runN;
  int lumi;


  std::string subdet_;

  const static int nMax = 800000;
  //float simhit_x[nMax],simhit_y[nMax],simhit_z[nMax];

  uint32_t simhit_detId;
  //,simhit_detId[nMax],simhit_detId[nMax],simhit_detId[nMax];

  float striphit_x,striphit_y,striphit_z;
  float striphit_match_x,striphit_match_y,striphit_match_z;
  float rechit_x,rechit_y,rechit_z;
  float rechit_match_x,rechit_match_y,rechit_match_z;

  int min_tothits;
  float nhitpTrack;
  float simtrk_pt;
  uint32_t    detId;

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  edm::ESGetToken<SetupData, SetupRecord> setupToken_;
#endif
};

StripClusterAnalysis::StripClusterAnalysis(const edm::ParameterSet& iConfig):trackerHitAssociatorConfig_(iConfig, consumesCollector()){

  subdet_ = iConfig.getParameter<std::string>("subdet");
  g4SimHitsCollectionToken_  = consumes<std::vector<PSimHit>>(iConfig.getParameter<edm::InputTag>("g4SimHitsCollection"));
  
  g4SimHitsTrackToken_ = consumes<std::vector<SimTrack>>(iConfig.getParameter<edm::InputTag>("g4SimHitsTrackTag"));   

  clusterToken_  = consumes<edmNew::DetSetVector<SiStripCluster>>(iConfig.getParameter<edm::InputTag>("siStripClustersTag"));
  tracksToken_ = consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("tracks"));
  
  tkGeomToken_ = esConsumes();
  usesResource("TFileService");

  clusterTree = fs->make<TTree>("clusterTree","clusterTree");
  clusterTree->Branch("event", &eventN, "event/i");
  clusterTree->Branch("run",   &runN, "run/I");
  clusterTree->Branch("lumi",  &lumi, "lumi/I");

  clusterTree->Branch("detId", &detId, "detId/i");
 
  clusterTree->Branch("simhit_detId",&simhit_detId,"simhit_detId/I");

  clusterTree->Branch("striphit_match_x",&striphit_match_x,"striphit_match_x/F");
  clusterTree->Branch("striphit_match_y",&striphit_match_y,"striphit_match_y/F");
  clusterTree->Branch("striphit_match_z",&striphit_match_z,"striphit_match_z/F");

 
  clusterTree->Branch("min_tothits",&min_tothits,"min_tothits/I");
  clusterTree->Branch("nhitpTrack",&nhitpTrack,"nhitpTrack/F");
  clusterTree->Branch("simtrk_pt",&simtrk_pt,"simtrk_pt/F");
  
  clusterTree->Branch("striphit_x",&striphit_x,"striphit_x/F");
  clusterTree->Branch("striphit_y",&striphit_y,"striphit_y/F");
  clusterTree->Branch("striphit_z",&striphit_z,"striphit_z/F");

  recohitTree = fs->make<TTree>("recohitTree","recohitTree");
  recohitTree->Branch("event", &eventN, "event/i");
  recohitTree->Branch("run",   &runN, "run/I");
  recohitTree->Branch("lumi",  &lumi, "lumi/I");

  recohitTree->Branch("simhit_detId",&simhit_detId,"simhit_detId/I");
  recohitTree->Branch("rechit_match_x",&rechit_match_x,"rechit_match_x/F");
  recohitTree->Branch("rechit_match_y",&rechit_match_y,"rechit_match_y/F");
  recohitTree->Branch("rechit_match_z",&rechit_match_z,"rechit_match_z/F");



}
StripClusterAnalysis::~StripClusterAnalysis() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty
}

//
// member functions
//

// ------------ method called for each event  ------------
void StripClusterAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  
  const auto& tkGeom = &iSetup.getData(tkGeomToken_);
  const auto tkDets = tkGeom->dets();

  edm::Handle<std::vector<PSimHit>> g4SimHitsCollection = iEvent.getHandle(g4SimHitsCollectionToken_);

  edm::Handle<std::vector<SimTrack>> g4SimHitsTrack = iEvent.getHandle(g4SimHitsTrackToken_);

  edm::Handle<edmNew::DetSetVector<SiStripCluster>> clusterCollection = iEvent.getHandle(clusterToken_);

  const auto& tracksHandle = iEvent.getHandle(tracksToken_);
  TrackerHitAssociator hitAssociator(iEvent, trackerHitAssociatorConfig_);

  
  if (!g4SimHitsCollection.isValid()){
    edm::LogError("TrackAnalyzer") << "No valid g4SimHits collection found";
    return;  
  }
  
  std::map<int,SimTrack> trackMap;
  for (const auto& track : *g4SimHitsTrack){
    trackMap[track.trackId()] = track;

  }


  for (const auto& simhit : *g4SimHitsCollection){
    eventN = iEvent.id().event();
    runN   = (int) iEvent.id().run();
    lumi  = (int) iEvent.id().luminosityBlock();

    int trackId = simhit.trackId();
    auto trackIt = trackMap.find(trackId);

    if(trackIt == trackMap.end())continue;

    for (const auto& detSiStripClusters : *clusterCollection){
    
      detId = detSiStripClusters.id();
      const auto& _detId = detId;
      auto det = std::find_if(tkDets.begin(), tkDets.end(), [_detId](auto& elem)->bool{
          return (elem->geographicalId().rawId() == _detId);
        });
        
      const StripTopology& p = dynamic_cast<const StripGeomDetUnit*>(*det)->specificTopology();
      for (const auto& cluster : detSiStripClusters){
	uint16_t barycenter = cluster.barycenter();
	LocalPoint localPos = p.localPosition((float) barycenter);
	GlobalPoint gp_clus = (tkGeom->idToDet(detId))->surface().toGlobal(localPos);
	GlobalPoint gp_hit = (tkGeom->idToDet(simhit.detUnitId()))->surface().toGlobal(simhit.localPosition());;
      

	simhit_detId = simhit.detUnitId();
	striphit_x = gp_clus.x();  // hit position at barycenter
	striphit_y = gp_clus.y();  // hit position at barycenter
	striphit_z = gp_clus.z();  // hit position at barycenter
      
	if(subdet_ == "TIB" || subdet_ == "TOB"){
	  if(std::fabs(simhit.localPosition().x() - localPos.x()) < 0.01){
	    auto hit_gp = HitPosition(simhit,detId,tkGeom);
	    striphit_match_x = hit_gp.x();
	    striphit_match_y = hit_gp.y();
	    striphit_match_z = hit_gp.z();
	    break;
	  }
	}
	else{
	  if (std::fabs(gp_hit.z() - gp_clus.z()) <= 0.0001){
	    auto hit_gp = HitPosition(simhit,detId,tkGeom);
	    striphit_match_x = hit_gp.x();
	    striphit_match_y = hit_gp.y();
	    striphit_match_z = hit_gp.z();
	    break;
	  }
	}
	clusterTree->Fill();
      }
    }

    
  }

  for(const auto& track : *tracksHandle){
    eventN = iEvent.id().event();
    runN   = (int) iEvent.id().run();
    lumi  = (int) iEvent.id().luminosityBlock();

    for (auto const &hit : track.recHits()){
      if (!hit->isValid())
        continue;
      
      DetId detid = hit->geographicalId();
      int subDet = detid.subdetId();

      bool hitInStrip = false;
      if (subdet_ == "TIB")
	hitInStrip = (subDet == SiStripDetId::TIB);
      else if (subdet_ == "TID")
	hitInStrip = (subDet == SiStripDetId::TID);
      else if (subdet_ == "TOB")
	hitInStrip = (subDet == SiStripDetId::TOB); 
      else
	hitInStrip = (subDet == SiStripDetId::TEC);

      if (hitInStrip){

	const std::type_info &type = typeid(*hit);
        if (type == typeid(SiStripRecHit1D)){

          const SiStripRecHit1D *striphit = dynamic_cast<const SiStripRecHit1D *>(hit);
	  
	  if(striphit != nullptr){
	    std::vector<PSimHit> simHitsAssociated = hitAssociator.associateHit(*striphit);
            for (auto const& simHit : simHitsAssociated){
	      auto hit_gp = HitPosition(simHit,detid,tkGeom);
	      simhit_detId = simHit.detUnitId();
	      rechit_match_x = hit_gp.x();
	      rechit_match_y = hit_gp.y();
	      rechit_match_z = hit_gp.z();
	      recohitTree->Fill();
	    }
	  }
	}

	else if (type == typeid(SiStripRecHit2D)){
	  const SiStripRecHit2D *striphit = dynamic_cast<const SiStripRecHit2D *>(hit);
	  
	  if(striphit != nullptr){
	    std::vector<PSimHit> simHitsAssociated = hitAssociator.associateHit(*striphit);
            for (auto const& simHit : simHitsAssociated){
	      auto hit_gp = HitPosition(simHit,detid,tkGeom);
	      simhit_detId = simHit.detUnitId();
	      rechit_match_x = hit_gp.x();
	      rechit_match_y = hit_gp.y();
	      rechit_match_z = hit_gp.z();
	      recohitTree->Fill();
	    }
	  }
	}
      }
    }
  }
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  // if the SetupData is always needed
  auto setup = iSetup.getData(setupToken_);
  // if need the ESHandle to check if the SetupData was there or not
  auto pSetup = iSetup.getHandle(setupToken_);
#endif
}

// ------------ method called once each job just before starting event loop  ------------
void StripClusterAnalysis::beginJob() {
  // please remove this method if not needed
}

// ------------ method called once each job just after ending the event loop  ------------
void StripClusterAnalysis::endJob() {
  // please remove this method if not needed
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void StripClusterAnalysis::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  // edm::ParameterSetDescription desc;
  // desc.setUnknown();
  // descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //edm::ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks", edm::InputTag("ctfWithMaterialTracks"));
  //descriptions.addWithDefaultLabel(desc);
}

GlobalPoint StripClusterAnalysis::HitPosition(const PSimHit& hit,const DetId detId, const  TrackerGeometry* tkGeom){


  GlobalPoint globalPosition;
  //const DetId detId = hit.detUnitId();
  const auto* geomDet = tkGeom->idToDet(detId);
  if (geomDet){
    const auto localPosition = hit.localPosition();
    globalPosition = geomDet->surface().toGlobal(localPosition);

  }
  return globalPosition;

}
//define this as a plug-in
DEFINE_FWK_MODULE(StripClusterAnalysis);
