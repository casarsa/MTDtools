#include <memory>
#include <iostream>
#include <vector>
#include <unordered_map>
#include <tuple>


#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"

#include "DataFormats/ForwardDetId/interface/MTDDetId.h"
#include "DataFormats/ForwardDetId/interface/BTLDetId.h"
#include "DataFormats/FTLDigi/interface/FTLDigiCollections.h"
#include "DataFormats/FTLRecHit/interface/FTLRecHitCollections.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "CLHEP/Units/GlobalPhysicalConstants.h"

#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TEfficiency.h"


struct MTDinfo {

  float sim_energy;
  float sim_time;

  uint32_t digi_charge;
  uint32_t digi_time1;
  uint32_t digi_time2;

  float ureco_charge;
  float ureco_time;

  float reco_energy;
  float reco_time;

};



bool orderByDetIdThenTime(const std::tuple<PSimHit,uint32_t,float> &a, 
			  const std::tuple<PSimHit,uint32_t,float> &b) {

  unsigned int detId_a(std::get<1>(a)), detId_b(std::get<1>(b));
      
  if(detId_a<detId_b) return true;
  if(detId_a>detId_b) return false;
      
  double time_a(std::get<2>(a)), time_b(std::get<2>(b));
  if(time_a<time_b) return true;
  
  return false;

}


class MTDAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {

public:
  explicit MTDAnalyzer(const edm::ParameterSet&);
  ~MTDAnalyzer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  // ----------member data ---------------------------

  const float btlMinEnergy_;


  //edm::EDGetTokenT<reco::GenParticleCollection> tok_genPart; 
  //edm::EDGetTokenT<edm::SimTrackContainer> tok_simTrack; 


  // --- Tracking Particles
  //edm::EDGetTokenT<TrackingParticleCollection> tok_trkPart; 

  // --- MTD SIM hits
  edm::EDGetTokenT<edm::PSimHitContainer> tok_BTL_sim; 
  edm::EDGetTokenT<edm::PSimHitContainer> tok_ETL_sim; 

  // --- MTD DIGI hits
  edm::EDGetTokenT<BTLDigiCollection> tok_BTL_digi; 
  edm::EDGetTokenT<ETLDigiCollection> tok_ETL_digi; 
  
  // --- MTD uncalibrated RECO hits 
  edm::EDGetTokenT<FTLUncalibratedRecHitCollection> tok_BTL_ureco; 
  edm::EDGetTokenT<FTLUncalibratedRecHitCollection> tok_ETL_ureco; 
  
  // --- MTD RECO hits 
  edm::EDGetTokenT<FTLRecHitCollection> tok_BTL_reco; 
  edm::EDGetTokenT<FTLRecHitCollection> tok_ETL_reco; 

  
  std::unordered_map<uint32_t, std::set<int> > n_btl_simHits;
  std::unordered_map<uint32_t, std::set<int> > n_etl_simHits;
  std::unordered_map<uint32_t, MTDinfo> btl_hits;
  std::unordered_map<uint32_t, MTDinfo> etl_hits;
  

  // --- Histograms declaration

  TH1F *h_n_sim_trk;
  TH1F *h_n_sim_cell;
  TH2F *h_occupancy_sim;
  TH1F *h_t_sim;
  TH1F *h_e_sim;
  TEfficiency* e_e_sim_digi;
  TEfficiency* e_e_sim_ureco;
  TEfficiency* e_e_sim_reco;
  TEfficiency* e_phi_sim_digi;
  TEfficiency* e_phi_sim_ureco;
  TEfficiency* e_phi_sim_reco;
  TEfficiency* e_eta_sim_digi;
  TEfficiency* e_eta_sim_ureco;
  TEfficiency* e_eta_sim_reco;
  
  TH1F *h_phi_sim;
  TH1F *h_eta_sim;

  TProfile *p_e_eta_sim;
  TProfile *p_t_eta_sim;
  TProfile *p_e_phi_sim;
  TProfile *p_t_phi_sim;


  TH1F *h_n_digi;
  TH2F *h_occupancy_digi;
  TH1F *h_t1_digi;
  TH1F *h_t2_digi;
  TH1F *h_e_digi;
  TH1F *h_phi_digi;
  TH1F *h_eta_digi;

  TProfile *p_e_eta_digi;
  TProfile *p_t_eta_digi;
  TProfile *p_e_phi_digi;
  TProfile *p_t_phi_digi;



  TH1F *h_n_ureco;
  TH2F *h_occupancy_ureco;
  TH1F *h_t_ureco;
  TH1F *h_t_ureco_uncorr;
  TH1F *h_e_ureco;

  TH1F *h_n_reco;
  TH2F *h_occupancy_reco;
  TH1F *h_t_reco;
  TH1F *h_t_reco_uncorr;
  TH1F *h_e_reco;

  TH1F *h_t_res;
  TH1F *h_t_res_uncorr;
  TH1F *h_e_res;

  TH2F *h_t_amp_ureco;
  TProfile *p_t_amp_ureco;

  TH2F *h_t_reco_sim;
  TH2F *h_e_reco_sim;

};


MTDAnalyzer::MTDAnalyzer(const edm::ParameterSet& iConfig) :
  btlMinEnergy_( iConfig.getParameter<double>("BTLMinimumEnergy") ) {

  tok_BTL_sim = consumes<edm::PSimHitContainer>(edm::InputTag("g4SimHits","FastTimerHitsBarrel"));
  tok_ETL_sim = consumes<edm::PSimHitContainer>(edm::InputTag("g4SimHits","FastTimerHitsEndcap"));

  tok_BTL_digi = consumes<BTLDigiCollection>(edm::InputTag("mix","FTLBarrel"));
  tok_ETL_digi = consumes<ETLDigiCollection>(edm::InputTag("mix","FTLEndcap"));

  tok_BTL_ureco = consumes<FTLUncalibratedRecHitCollection>(edm::InputTag("mtdUncalibratedRecHits","FTLBarrel"));
  tok_ETL_ureco = consumes<FTLUncalibratedRecHitCollection>(edm::InputTag("mtdUncalibratedRecHits","FTLEndcap"));

  tok_BTL_reco = consumes<FTLRecHitCollection>(edm::InputTag("mtdRecHits","FTLBarrel"));
  tok_ETL_reco = consumes<FTLRecHitCollection>(edm::InputTag("mtdRecHits","FTLEndcap"));

  edm::Service<TFileService> fs;


  // --- Histograms definition

  h_n_sim_trk   = fs->make<TH1F>("h_n_sim_trk", "Number of tracks per BTL cell;N_{trk}", 10, 0., 10.);
  h_n_sim_cell = fs->make<TH1F>("h_n_sim_cell", "Number of BTL cells with SIM hits;N_{BTL cells}", 100, 0., 100.);
  h_occupancy_sim = fs->make<TH2F>("h_occupancy_sim", "SIM hits occupancy;cell #phi;cell #eta", 
				   145, 0., 2305., 86, -43., 43.);
  h_t_sim   = fs->make<TH1F>("h_t_sim", "SIM hits ToA;ToA [ns]", 500, 0., 50.);
  h_e_sim   = fs->make<TH1F>("h_e_sim", "SIM hits energy;E [MeV]", 200, 0., 20.);

  e_e_sim_digi = fs->make<TEfficiency>("e_e_sim_digi", "DIGI efficiency;E_{SIM} [MeV]", 100, 0., 20.);
  e_e_sim_ureco= fs->make<TEfficiency>("e_e_sim_ureco", "URECO efficiency;E_{SIM} [MeV]", 100, 0., 20.);
  e_e_sim_reco = fs->make<TEfficiency>("e_e_sim_reco", "RECO efficiency;E_{SIM} [MeV]", 100, 0., 20.);
  e_phi_sim_digi = fs->make<TEfficiency>("e_phi_sim_digi", "DIGI efficiency;cell #phi", 145, 0., 2305.);
  e_phi_sim_ureco= fs->make<TEfficiency>("e_phi_sim_ureco", "URECO efficiency;cell #phi", 145, 0., 2305.);
  e_phi_sim_reco = fs->make<TEfficiency>("e_phi_sim_reco", "RECO efficiency;cell #phi", 145, 0., 2305.);
  e_eta_sim_digi = fs->make<TEfficiency>("e_eta_sim_digi", "DIGI efficiency;cell #eta", 43, 0., 43.);
  e_eta_sim_ureco= fs->make<TEfficiency>("e_eta_sim_ureco", "URECO efficiency;cell #eta", 43, 0., 43.);
  e_eta_sim_reco = fs->make<TEfficiency>("e_eta_sim_reco", "RECO efficiency;cell #eta", 43, 0., 43.);

  h_phi_sim = fs->make<TH1F>("h_phi_sim", "SIM hits #phi;#phi index", 145, 0., 2305.);
  h_eta_sim = fs->make<TH1F>("h_eta_sim", "SIM hits #eta;#eta index", 86, -43., 43.);

  p_e_eta_sim = fs->make<TProfile>("p_e_eta_sim", "SIM energy vs |#eta|;cell |#eta|;E_{SIM} [MeV]", 
				    43, 0., 43.);
  p_t_eta_sim = fs->make<TProfile>("p_t_eta_sim", "SIM time vs |#eta|;cell |#eta|;T_{SIM} [ns]", 
				    43, 0., 43.);
  p_e_phi_sim = fs->make<TProfile>("p_e_phi_sim", "SIM energy vs #phi;cell #phi;E_{SIM} [MeV]", 
				    145, 0., 2305.);
  p_t_phi_sim = fs->make<TProfile>("p_t_phi_sim", "SIM time vs #phi;cell #phi;T_{SIM} [ns]", 
				    145, 0., 2305.);


  h_n_digi  = fs->make<TH1F>("h_n_digi", "Number of DIGI hits;N_{DIGI hits}", 100, 0., 100.);
  h_occupancy_digi = fs->make<TH2F>("h_occupancy_digi", "DIGI hits occupancy;cell #phi;cell #eta", 
				    145, 0., 2305., 86, -43., 43.);
  h_t1_digi = fs->make<TH1F>("h_t1_digi", "DIGI hits ToA (1);ToA [TDC counts]", 1024, 0., 1024.);
  h_t2_digi = fs->make<TH1F>("h_t2_digi", "DIGI hits ToA (2);ToA [TDC counts]", 1024, 0., 1024.);
  h_e_digi  = fs->make<TH1F>("h_e_digi", "DIGI hits energy;amplitude [ADC counts]", 1024, 0., 1024.);

  h_phi_digi = fs->make<TH1F>("h_phi_digi", "DIGI hits #phi;#phi index", 145, 0., 2305.);
  h_eta_digi = fs->make<TH1F>("h_eta_digi", "DIGI hits #eta;#eta index", 86, -43., 43.);


  p_e_eta_digi = fs->make<TProfile>("p_e_eta_digi", "DIGI charge vs |#eta|;cell |#eta|;ADC counts", 
				    43, 0., 43.);
  p_t_eta_digi = fs->make<TProfile>("p_t_eta_digi", "DIGI time vs |#eta|;cell |#eta|;TDC counts", 
				    43, 0., 43.);
  p_e_phi_digi = fs->make<TProfile>("p_e_phi_digi", "DIGI charge vs #phi;cell #phi;ADC counts", 
				    145, 0., 2305.);
  p_t_phi_digi = fs->make<TProfile>("p_t_phi_digi", "DIGI time vs #phi;cell #phi;TDC counts", 
				    145, 0., 2305.);

  h_n_ureco  = fs->make<TH1F>("h_n_ureco", "Number of URECO hits;N_{URECO hits}", 100, 0., 100.);
  h_occupancy_ureco = fs->make<TH2F>("h_occupancy_ureco", "URECO hits occupancy;cell #phi;cell #eta", 
				     145, 0., 2305., 86, -43., 43.);
  h_t_ureco  = fs->make<TH1F>("h_t_ureco", "URECO hits ToA;ToA [ns]", 250, 0., 25.);
  h_t_ureco_uncorr = fs->make<TH1F>("h_t_ureco_uncorr", "URECO hits ToA;ToA [ns]", 250, 0., 25.);
  h_e_ureco  = fs->make<TH1F>("h_e_ureco", "URECO hits energy;Q [pC]", 300, 0., 600.);

  h_n_reco  = fs->make<TH1F>("h_n_reco", "Number of RECO hits;N_{RECO hits}", 100, 0., 100.);
  h_occupancy_reco = fs->make<TH2F>("h_occupancy_reco", "RECO hits occupancy;cell #phi;cell #eta", 
				    145, 0., 2305., 86, -43., 43.);
  h_t_reco  = fs->make<TH1F>("h_t_reco", "RECO hits ToA;ToA [ns]", 250, 0., 25.);
  h_t_reco_uncorr = fs->make<TH1F>("h_t_reco_uncorr", "RECO hits ToA;ToA [ns]", 250, 0., 25.);
  h_e_reco  = fs->make<TH1F>("h_e_reco", "RECO hits energy;E [MeV]", 200, 0., 20.);

  h_t_res  = fs->make<TH1F>("h_t_res", "ToA resolution;ToA [ns]", 700, -2., 5.);
  h_t_res_uncorr = fs->make<TH1F>("h_t_res_uncorr", "ToA resolution;ToA [ns]", 700, -2., 5.);
  h_e_res  = fs->make<TH1F>("h_e_res", "Energy resolution;E [MeV]", 200, -1., 1.);

  h_t_reco_sim = fs->make<TH2F>("h_t_reco_sim", "ToA reco vs sim;SIM ToA [ns];RECO ToA [ns]", 
				100, -1., 25., 100, 0., 25.);
  h_e_reco_sim = fs->make<TH2F>("h_e_reco_sim", "E reco vs sim;SIM E [MeV];RECO E [MeV]", 
				100, 0., 20., 100, 0., 20.);

  h_t_amp_ureco = fs->make<TH2F>("h_t_amp_ureco", "time vs amplitude;amplitude [pC];time [ns]", 
				100, 0., 600., 400, 0., 20.);

  p_t_amp_ureco = fs->make<TProfile>("p_t_amp_ureco", "time vs amplitude;amplitude [pC];time [ns]", 
				     100, 0., 600.);


}


MTDAnalyzer::~MTDAnalyzer() {}


//
// member functions
//

void
MTDAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  using namespace std;


  edm::Handle<edm::PSimHitContainer>  h_BTL_sim;
  iEvent.getByToken( tok_BTL_sim, h_BTL_sim );
  edm::Handle<edm::PSimHitContainer>  h_ETL_sim;
  iEvent.getByToken( tok_ETL_sim, h_ETL_sim );

  edm::Handle<BTLDigiCollection>   h_BTL_digi;
  iEvent.getByToken( tok_BTL_digi, h_BTL_digi );
  edm::Handle<ETLDigiCollection>   h_ETL_digi;
  iEvent.getByToken( tok_ETL_digi, h_ETL_digi );

  edm::Handle<FTLUncalibratedRecHitCollection> h_BTL_ureco;
  iEvent.getByToken( tok_BTL_ureco, h_BTL_ureco );
  edm::Handle<FTLUncalibratedRecHitCollection> h_ETL_ureco;
  iEvent.getByToken( tok_ETL_ureco, h_ETL_ureco );

  edm::Handle<FTLRecHitCollection> h_BTL_reco;
  iEvent.getByToken( tok_BTL_reco, h_BTL_reco );
  edm::Handle<FTLRecHitCollection> h_ETL_reco;
  iEvent.getByToken( tok_ETL_reco, h_ETL_reco );


  std::set<uint32_t> unique_simHit;



  // ==============================================================================
  //  MTD SIM hits
  // ==============================================================================

  // --- BTL
  if ( h_BTL_sim->size() > 0 ) {
  
    std::vector< std::tuple<PSimHit,uint32_t,float> > hitRefs;
    hitRefs.reserve(h_BTL_sim->size());

    // Sort the SimHits per detector id and time 
    for (auto const& simHit: *h_BTL_sim) {

      DetId id = simHit.detUnitId();
      if (id.rawId() != 0)
	hitRefs.emplace_back( simHit, id.rawId(), simHit.tof() );

    } // simHit loop 
    std::sort(hitRefs.begin(),hitRefs.end(), orderByDetIdThenTime);
    
    // Accumulate the SimHits in the same detector cell
    for (auto const& hitRef: hitRefs) {
    
      const PSimHit &hit = std::get<0>(hitRef);
      DetId id = hit.detUnitId();

      unique_simHit.insert(id.rawId());

      auto simHitIt = btl_hits.emplace(id.rawId(),MTDinfo()).first;

      (simHitIt->second).sim_energy += 1000.*hit.energyLoss();

      n_btl_simHits[id.rawId()].insert(hit.trackId());

      // Get the time of the first SimHit in the cell
      if( (simHitIt->second).sim_time==0 )
	(simHitIt->second).sim_time = hit.tof();

    } // hitRef loop

  } // if ( h_BTL_sim->size() > 0 )



  // ==============================================================================
  //  MTD DIGI hits
  // ==============================================================================

  unsigned int n_digi_btl = 0;

  if ( h_BTL_digi->size() > 0 ) {

    for (const auto& dataFrame: *h_BTL_digi) {

      DetId id =  dataFrame.id();

      // --- loop over the dataFrame samples
      for (int isample = 0; isample<dataFrame.size(); ++isample){

	const auto& sample = dataFrame.sample(isample);
	
	if ( sample.data()!=0 && sample.toa()!=0 ) {

	  // on-time sample
	  if ( isample == 2 ) {

	    btl_hits[id.rawId()].digi_charge =  sample.data();
	    btl_hits[id.rawId()].digi_time1  =  sample.toa();
	    btl_hits[id.rawId()].digi_time2  =  sample.toa2();

	    n_digi_btl++;

	  }
	  //else {
	  //  std::cout << "*** WARNING: Sample " << isample << " is not empty!" << std::endl;
	  //  std::cout << "             amplitude = " << sample.data()
	  //	      << "  time 1 = " << sample.toa() << "  time 2 = " <<  sample.toa2()
	  //	      << std::endl;
	  //
	  //}

	} 

      } // isaple loop

    } // dataFrame loop

  } // if ( h_BTL_digi->size() > 0 )



  // ==============================================================================
  //  MTD Uncalibrated RECO hits
  // ==============================================================================

  unsigned int n_ureco_btl = 0;

  if ( h_BTL_ureco->size() > 0 ) {

    for (const auto& urecHit: *h_BTL_ureco) {

      DetId id = urecHit.id();

      btl_hits[id.rawId()].ureco_charge = urecHit.amplitude();
      btl_hits[id.rawId()].ureco_time   = urecHit.time();

      if ( urecHit.amplitude() > 0. )
	n_ureco_btl++;
      

    } // recHit loop

  } // if ( h_BTL_reco->size() > 0 )



  // ==============================================================================
  //  MTD RECO hits
  // ==============================================================================

  unsigned int n_reco_btl = 0;

  if ( h_BTL_reco->size() > 0 ) {

    for (const auto& recHit: *h_BTL_reco) {

      DetId id = recHit.id();

      btl_hits[id.rawId()].reco_energy = recHit.energy();
      btl_hits[id.rawId()].reco_time   = recHit.time();

      if ( recHit.energy() > 0. )
	n_reco_btl++;


    } // recHit loop

  } // if ( h_BTL_reco->size() > 0 )



  // ==============================================================================
  //  Fill the histograms
  // ==============================================================================

  for (auto const& hit: n_btl_simHits) {
    h_n_sim_trk->Fill((hit.second).size());
  }

  h_n_sim_cell->Fill(unique_simHit.size());
  h_n_digi->Fill(n_digi_btl);
  h_n_ureco->Fill(n_ureco_btl);
  h_n_reco->Fill(n_reco_btl);


  for (auto const& hit: btl_hits) {

    BTLDetId detId(hit.first); 

    //if ( (hit.second).ureco_charge < 10. ) continue;

    h_t_amp_ureco->Fill((hit.second).ureco_charge,(hit.second).ureco_time);
    p_t_amp_ureco->Fill((hit.second).ureco_charge,(hit.second).ureco_time);
    

    //if ( (hit.second).reco_energy < btlMinEnergy_ ) continue;
    
    //if ( (hit.second).sim_energy < 1. ) continue;

    e_e_sim_digi->Fill((hit.second).digi_charge > 0., (hit.second).sim_energy);
    e_e_sim_ureco->Fill((hit.second).ureco_charge > 0., (hit.second).sim_energy);
    e_e_sim_reco->Fill((hit.second).reco_energy > 0., (hit.second).sim_energy);

    e_phi_sim_digi->Fill((hit.second).digi_charge > 0.,   detId.iphi(BTLDetId::CrysLayout::barzflat));
    e_phi_sim_ureco->Fill((hit.second).ureco_charge > 0., detId.iphi(BTLDetId::CrysLayout::barzflat));
    e_phi_sim_reco->Fill((hit.second).reco_energy > 0.,   detId.iphi(BTLDetId::CrysLayout::barzflat));

    e_eta_sim_digi->Fill((hit.second).digi_charge > 0.,   detId.ietaAbs(BTLDetId::CrysLayout::barzflat));
    e_eta_sim_ureco->Fill((hit.second).ureco_charge > 0., detId.ietaAbs(BTLDetId::CrysLayout::barzflat));
    e_eta_sim_reco->Fill((hit.second).reco_energy > 0.,   detId.ietaAbs(BTLDetId::CrysLayout::barzflat));

    if ( (hit.second).sim_energy == 0. ) continue;
    //if ( (hit.second).sim_energy < 2. ) continue;

    h_e_sim->Fill((hit.second).sim_energy);
    h_t_sim->Fill((hit.second).sim_time);

    h_phi_sim->Fill(detId.iphi(BTLDetId::CrysLayout::barzflat));
    h_eta_sim->Fill(detId.ieta(BTLDetId::CrysLayout::barzflat));
    


    h_occupancy_sim->Fill(detId.iphi(BTLDetId::CrysLayout::barzflat),detId.ieta(BTLDetId::CrysLayout::barzflat));

    p_e_eta_sim->Fill(detId.ietaAbs(BTLDetId::CrysLayout::barzflat),(hit.second).sim_energy);
    p_t_eta_sim->Fill(detId.ietaAbs(BTLDetId::CrysLayout::barzflat),(hit.second).sim_time);
    p_e_phi_sim->Fill(detId.iphi(BTLDetId::CrysLayout::barzflat),(hit.second).sim_energy);
    p_t_phi_sim->Fill(detId.iphi(BTLDetId::CrysLayout::barzflat),(hit.second).sim_time);


    if ( (hit.second).digi_charge == 0 ) continue;

    h_occupancy_digi->Fill(detId.iphi(BTLDetId::CrysLayout::barzflat),detId.ieta(BTLDetId::CrysLayout::barzflat));

    h_e_digi ->Fill((hit.second).digi_charge);
    h_t1_digi->Fill((hit.second).digi_time1);
    h_t2_digi->Fill((hit.second).digi_time2);
    
    h_phi_digi->Fill(detId.iphi(BTLDetId::CrysLayout::barzflat));
    h_eta_digi->Fill(detId.ieta(BTLDetId::CrysLayout::barzflat));

    p_e_eta_digi->Fill(detId.ietaAbs(BTLDetId::CrysLayout::barzflat),(hit.second).digi_charge);
    p_t_eta_digi->Fill(detId.ietaAbs(BTLDetId::CrysLayout::barzflat),(hit.second).digi_time1);
    p_e_phi_digi->Fill(detId.iphi(BTLDetId::CrysLayout::barzflat),(hit.second).digi_charge);
    p_t_phi_digi->Fill(detId.iphi(BTLDetId::CrysLayout::barzflat),(hit.second).digi_time1);


    if ( (hit.second).ureco_charge == 0. ) continue;

    h_occupancy_ureco->Fill(detId.iphi(BTLDetId::CrysLayout::barzflat),detId.ieta(BTLDetId::CrysLayout::barzflat));

    h_e_ureco->Fill((hit.second).ureco_charge);
    h_t_ureco->Fill((hit.second).ureco_time);



    if ( (hit.second).reco_energy == 0. ) continue;

    h_occupancy_reco->Fill(detId.iphi(BTLDetId::CrysLayout::barzflat),detId.ieta(BTLDetId::CrysLayout::barzflat));

    h_e_reco->Fill((hit.second).reco_energy);
    h_t_reco->Fill((hit.second).reco_time);
    
    h_e_res->Fill((hit.second).reco_energy-(hit.second).sim_energy);
    h_t_res->Fill((hit.second).reco_time-(hit.second).sim_time);
    
    h_t_reco_sim->Fill((hit.second).sim_time,(hit.second).reco_time);
    h_e_reco_sim->Fill((hit.second).sim_energy,(hit.second).reco_energy);



    // --- Reverse the global correction

    float ureco_charge = (hit.second).ureco_charge;

    const float p0 = 24.8997;
    const float p1 = -0.911385;
    const float p2 =  3.3744717;

    float time_corr = p0*pow(ureco_charge,p1) + p2;

    float ureco_time_uncorr = (hit.second).ureco_time + time_corr;

    float reco_time_uncorr = ureco_time_uncorr;


    h_t_ureco_uncorr->Fill(ureco_time_uncorr);
    h_t_reco_uncorr->Fill(reco_time_uncorr);
    h_t_res_uncorr->Fill(reco_time_uncorr-(hit.second).sim_time);


  } // hit loop


  n_btl_simHits.clear();
  n_etl_simHits.clear();
  btl_hits.clear();
  etl_hits.clear();

}


// ------------ method called once each job just before starting event loop  ------------
void 
MTDAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MTDAnalyzer::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MTDAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MTDAnalyzer);
