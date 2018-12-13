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

#include "Geometry/Records/interface/MTDDigiGeometryRecord.h"
#include "Geometry/MTDGeometryBuilder/interface/MTDGeometry.h"
#include "Geometry/MTDGeometryBuilder/interface/ProxyMTDTopology.h"
#include "Geometry/MTDGeometryBuilder/interface/RectangularMTDTopology.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "CLHEP/Units/GlobalPhysicalConstants.h"

#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"



struct MTDinfo {

  float sim_energy;
  float sim_time;
  float sim_x;
  float sim_y;
  float sim_z;

  uint32_t digi_charge[2];
  uint32_t digi_time1[2];
  uint32_t digi_time2[2];

  float ureco_charge[2];
  float ureco_time[2];

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

  const MTDGeometry* geom_;

  const float btlIntegrationWindow_;
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
  std::unordered_map<uint32_t, std::set<int> > n_etl_simHits[2];
  std::unordered_map<uint32_t, MTDinfo> btl_hits;
  std::unordered_map<uint32_t, MTDinfo> etl_hits[2];
  

  ///////////////////////////////////////////////////////////////////////////////////////////////
  //
  //  Histograms declaration
  //
  ///////////////////////////////////////////////////////////////////////////////////////////////

  // --- BTL -------------------------------------------------------

  // SIM

  TH1F *hb_n_sim_trk;
  TH1F *hb_n_sim_cell;
  TH1F *hb_t_sim;
  TH1F *hb_e_sim;
  TH1F *hb_xloc_sim;
  TH1F *hb_yloc_sim;
  TH1F *hb_zloc_sim;
  TH2F *hb_occupancy_sim;
  TH1F *hb_phi_sim;
  TH1F *hb_z_sim;
  TH1F *hb_eta_sim;
 
  TH2F *hb_t_e_sim;
  TH2F *hb_e_eta_sim;
  TH2F *hb_t_eta_sim;
  TH2F *hb_e_phi_sim;
  TH2F *hb_t_phi_sim;

  TProfile *pb_t_e_sim;
  TProfile *pb_e_eta_sim;
  TProfile *pb_t_eta_sim;
  TProfile *pb_e_phi_sim;
  TProfile *pb_t_phi_sim;


  // DIGI

  TH1F *hb_n_digi[2];
  TH1F *hb_t1_digi[2];
  TH1F *hb_t2_digi[2];
  TH1F *hb_e_digi[2];
  TH2F *hb_occupancy_digi[2];
  TH1F *hb_phi_digi[2];
  TH1F *hb_eta_digi[2];

  TH2F *hb_t1_e_digi[2];
  TH2F *hb_t2_e_digi[2];
  TH2F *hb_e_eta_digi[2];
  TH2F *hb_t1_eta_digi[2];
  TH2F *hb_t2_eta_digi[2];
  TH2F *hb_e_phi_digi[2];
  TH2F *hb_t1_phi_digi[2];
  TH2F *hb_t2_phi_digi[2];

  TProfile *pb_t1_e_digi[2];
  TProfile *pb_t2_e_digi[2];
  TProfile *pb_e_eta_digi[2];
  TProfile *pb_t1_eta_digi[2];
  TProfile *pb_t2_eta_digi[2];
  TProfile *pb_e_phi_digi[2];
  TProfile *pb_t1_phi_digi[2];
  TProfile *pb_t2_phi_digi[2];


  // Uncalibrated RECO

  TH1F *hb_n_ureco[2];
  TH1F *hb_t_ureco[2];
  TH1F *hb_t_ureco_uncorr[2];
  TH1F *hb_e_ureco[2];
  TH2F *hb_occupancy_ureco[2];

  TH2F *hb_t_amp_ureco[2];
  TProfile *pb_t_amp_ureco[2];


  // RECO

  TH1F *hb_n_reco;
  TH2F *hb_occupancy_reco;
  TH1F *hb_t_reco;
  TH1F *hb_t_reco_uncorr;
  TH1F *hb_e_reco;

  TH1F *hb_t_res;
  TH1F *hb_t_res_uncorr;
  TH1F *hb_e_res;

  TH2F *hb_t_reco_sim;
  TH2F *hb_e_reco_sim;


  // --- ETL -------------------------------------------------------

  // SIM

  TH1F *he_n_sim_trk[2];
  TH1F *he_n_sim_cell[2];

  TH1F *he_t_sim[2];
  TH1F *he_e_sim[2];

  TH1F *he_xloc_sim[2];
  TH1F *he_yloc_sim[2];
  TH1F *he_zloc_sim[2];

  TH2F *he_occupancy_sim[2];
  TH1F *he_x_sim[2];
  TH1F *he_y_sim[2];
  TH1F *he_z_sim[2];
  TH1F *he_phi_sim[2];
  TH1F *he_eta_sim[2];

  TH2F *he_t_e_sim[2];
  TH2F *he_e_eta_sim[2];
  TH2F *he_t_eta_sim[2];
  TH2F *he_e_phi_sim[2];
  TH2F *he_t_phi_sim[2];

  TProfile *pe_t_e_sim[2];
  TProfile *pe_e_eta_sim[2];
  TProfile *pe_t_eta_sim[2];
  TProfile *pe_e_phi_sim[2];
  TProfile *pe_t_phi_sim[2];


  // DIGI

  TH1F *he_n_digi[2];
  TH1F *he_t_digi[2];
  TH1F *he_e_digi[2];


  // Uncalibrated RECO

  TH1F *he_n_ureco[2];

  // RECO

  TH1F *he_n_reco[2];


};


MTDAnalyzer::MTDAnalyzer(const edm::ParameterSet& iConfig) :
  geom_(nullptr),
  btlIntegrationWindow_( iConfig.getParameter<double>("BTLIntegrationWindow") ),
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


  ///////////////////////////////////////////////////////////////////////////////////////////////
  //
  //  Histograms definition
  //
  ///////////////////////////////////////////////////////////////////////////////////////////////

  TFileDirectory btl = fs->mkdir( "BTL" );
  TFileDirectory etl = fs->mkdir( "ETL" );

  // ==============================================================================
  //  BTL
  // ==============================================================================

  // --- SIM

  hb_n_sim_trk  = btl.make<TH1F>("h_n_sim_trk", "Number of tracks per BTL cell;N_{trk}", 10, 0., 10.);
  hb_n_sim_cell = btl.make<TH1F>("h_n_sim_cell", "Number of BTL cells with SIM hits;N_{BTL cells}", 250, 0., 5000.);

  hb_t_sim = btl.make<TH1F>("h_t_sim", "BTL SIM hits ToA;ToA_{SIM} [ns]", 250, 0., 25.);
  hb_e_sim = btl.make<TH1F>("h_e_sim", "BTL SIM hits energy;E_{SIM} [MeV]", 200, 0., 20.);
  hb_xloc_sim = btl.make<TH1F>("h_xloc_sim", "BTL SIM local x;x_{SIM} [mm]", 290, -1.45, 1.45);
  hb_yloc_sim = btl.make<TH1F>("h_yloc_sim", "BTL SIM local y;y_{SIM} [mm]", 600, -30., 30.);
  hb_zloc_sim = btl.make<TH1F>("h_zloc_sim", "BTL SIM local z;z_{SIM} [mm]", 400, -2., 2.);

  hb_occupancy_sim = btl.make<TH2F>("h_occupancy_sim", "BTL SIM hits occupancy;z_{SIM} [cm];#phi_{SIM} [rad]",
				    520, -260., 260., 315, -3.15, 3.15 );
  hb_phi_sim = btl.make<TH1F>("h_phi_sim", "BTL SIM hits #phi;#phi_{SIM} [rad]", 315, -3.15, 3.15);
  hb_z_sim   = btl.make<TH1F>("h_z_sim", "BTL SIM hits z;z_{SIM} [cm]", 520, -260., 260.);
  hb_eta_sim = btl.make<TH1F>("h_eta_sim", "BTL SIM hits #eta;#eta_{SIM}", 200, -1.6, 1.6);
  hb_t_e_sim   = btl.make<TH2F>("h_t_e_sim", "BTL SIM time vs energy;E_{SIM} [MeV];T_{SIM} [ns]",
				100, 0., 20., 100, 0., 25.);
  hb_e_eta_sim = btl.make<TH2F>("h_e_eta_sim", "BTL SIM energy vs |#eta|;|#eta_{SIM}|;E_{SIM} [MeV]",
				100, 0., 1.6, 100, 0., 20.);
  hb_t_eta_sim = btl.make<TH2F>("h_t_eta_sim", "BTL SIM time vs |#eta|;|#eta_{SIM}|;T_{SIM} [ns]",
				100, 0., 1.6, 100, 0., 25.);
  hb_e_phi_sim = btl.make<TH2F>("h_e_phi_sim", "BTL SIM energy vs #phi;#phi_{SIM} [rad];E_{SIM} [MeV]",
				100, -3.15, 3.15, 100, 0., 20.);
  hb_t_phi_sim = btl.make<TH2F>("h_t_phi_sim", "BTL SIM time vs #phi;#phi_{SIM} [rad];T_{SIM} [ns]",
				100, -3.15, 3.15, 100, 0., 25.);

  pb_t_e_sim   = btl.make<TProfile>("p_t_e_sim", "BTL SIM time vs energy;E_{SIM} [MeV];T_{SIM} [ns]",
				    100, 0., 20.);
  pb_e_eta_sim = btl.make<TProfile>("p_e_eta_sim", "BTL SIM energy vs |#eta|;|#eta_{SIM}|;E_{SIM} [MeV]",
				    100, 0., 1.6);
  pb_t_eta_sim = btl.make<TProfile>("p_t_eta_sim", "BTL SIM time vs |#eta|;|#eta_{SIM}|;T_{SIM} [ns]",
				    100, 0., 1.6);
  pb_e_phi_sim = btl.make<TProfile>("p_e_phi_sim", "BTL SIM energy vs #phi;#phi_{SIM} [rad];E_{SIM} [MeV]",
				    100, -3.15, 3.15);
  pb_t_phi_sim = btl.make<TProfile>("p_t_phi_sim", "BTL SIM time vs #phi;#phi_{SIM} [rad];T_{SIM} [ns]",
				    100, -3.15, 3.15);


  // --- DIGI

  hb_n_digi[0]  = btl.make<TH1F>("h_n_digi_0", "Number of BTL DIGI hits (L);N_{DIGI hits}", 100, 0., 100.);
  hb_n_digi[1]  = btl.make<TH1F>("h_n_digi_1", "Number of BTL DIGI hits (R);N_{DIGI hits}", 100, 0., 100.);
  hb_t1_digi[0] = btl.make<TH1F>("h_t1_digi_0", "BTL DIGI hits ToA1 (L);ToA [TDC counts]", 1024, 0., 1024.);
  hb_t1_digi[1] = btl.make<TH1F>("h_t1_digi_1", "BTL DIGI hits ToA1 (R);ToA [TDC counts]", 1024, 0., 1024.);
  hb_t2_digi[0] = btl.make<TH1F>("h_t2_digi_0", "BTL DIGI hits ToA2 (L);ToA [TDC counts]", 1024, 0., 1024.);
  hb_t2_digi[1] = btl.make<TH1F>("h_t2_digi_1", "BTL DIGI hits ToA2 (R);ToA [TDC counts]", 1024, 0., 1024.);
  hb_e_digi[0]  = btl.make<TH1F>("h_e_digi_0", "BTL DIGI hits energy (L);amplitude [ADC counts]", 1024, 0., 1024.);
  hb_e_digi[1]  = btl.make<TH1F>("h_e_digi_1", "BTL DIGI hits energy (R);amplitude [ADC counts]", 1024, 0., 1024.);

  hb_occupancy_digi[0] = btl.make<TH2F>("h_occupancy_digi_0", "BTL DIGI hits occupancy (L);cell #phi;cell #eta",
				       145, 0., 2305., 86, -43., 43.);
  hb_occupancy_digi[1] = btl.make<TH2F>("h_occupancy_digi_1", "BTL DIGI hits occupancy (R);cell #phi;cell #eta",
				       145, 0., 2305., 86, -43., 43.);
  hb_phi_digi[0] = btl.make<TH1F>("h_phi_digi_0", "BTL DIGI hits #phi (L);#phi index", 145, 0., 2305.);
  hb_phi_digi[1] = btl.make<TH1F>("h_phi_digi_1", "BTL DIGI hits #phi (R);#phi index", 145, 0., 2305.);
  hb_eta_digi[0] = btl.make<TH1F>("h_eta_digi_0", "BTL DIGI hits #eta (L);#eta index", 86, -43., 43.);
  hb_eta_digi[1] = btl.make<TH1F>("h_eta_digi_1", "BTL DIGI hits #eta (R);#eta index", 86, -43., 43.);


  hb_t1_e_digi[0]   = btl.make<TH2F>("h_t1_e_digi_0", "BTL DIGI time1 vs charge (L);ADC counts;TDC counts",
				    128, 0., 1024., 128, 0., 1024.);
  hb_t1_e_digi[1]   = btl.make<TH2F>("h_t1_e_digi_1", "BTL DIGI time1 vs charge (R);ADC counts;TDC counts",
				    128, 0., 1024., 128, 0., 1024.);
  hb_t2_e_digi[0]   = btl.make<TH2F>("h_t2_e_digi_0", "BTL DIGI time2 vs charge (L);ADC counts;TDC counts",
				    128, 0., 1024., 128, 0., 1024.);
  hb_t2_e_digi[1]   = btl.make<TH2F>("h_t2_e_digi_1", "BTL DIGI time2 vs charge (R);ADC counts;TDC counts",
				    128, 0., 1024., 128, 0., 1024.);
  hb_e_eta_digi[0]  = btl.make<TH2F>("h_e_eta_digi_0", "BTL DIGI charge vs |#eta| (L);cell |#eta|;ADC counts",
				    43, 0., 43., 128, 0., 1024.);
  hb_e_eta_digi[1]  = btl.make<TH2F>("h_e_eta_digi_1", "BTL DIGI charge vs |#eta| (R);cell |#eta|;ADC counts",
				    43, 0., 43., 128, 0., 1024.);
  hb_t1_eta_digi[0] = btl.make<TH2F>("h_t1_eta_digi_0", "BTL DIGI time1 vs |#eta| (L);cell |#eta|;TDC counts",
				    43, 0., 43., 128, 0., 1024.);
  hb_t1_eta_digi[1] = btl.make<TH2F>("h_t1_eta_digi_1", "BTL DIGI time1 vs |#eta| (R);cell |#eta|;TDC counts",
				    43, 0., 43., 128, 0., 1024.);
  hb_t2_eta_digi[0] = btl.make<TH2F>("h_t2_eta_digi_0", "BTL DIGI time2 vs |#eta| (L);cell |#eta|;TDC counts",
				    43, 0., 43., 128, 0., 1024.);
  hb_t2_eta_digi[1] = btl.make<TH2F>("h_t2_eta_digi_1", "BTL DIGI time2 vs |#eta| (R);cell |#eta|;TDC counts",
				    43, 0., 43., 128, 0., 1024.);
  hb_e_phi_digi[0]  = btl.make<TH2F>("h_e_phi_digi_0", "BTL DIGI charge vs #phi (L);cell #phi;ADC counts",
				    145, 0., 2305., 128, 0., 1024.);
  hb_e_phi_digi[1]  = btl.make<TH2F>("h_e_phi_digi_1", "BTL DIGI charge vs #phi (R);cell #phi;ADC counts",
				    145, 0., 2305., 128, 0., 1024.);
  hb_t1_phi_digi[0] = btl.make<TH2F>("h_t1_phi_digi_0", "BTL DIGI time1 vs #phi (L);cell #phi;TDC counts",
				    145, 0., 2305., 128, 0., 1024.);
  hb_t1_phi_digi[1] = btl.make<TH2F>("h_t1_phi_digi_1", "BTL DIGI time1 vs #phi (R);cell #phi;TDC counts",
				    145, 0., 2305., 128, 0., 1024.);
  hb_t2_phi_digi[0] = btl.make<TH2F>("h_t2_phi_digi_0", "BTL DIGI time2 vs #phi (L);cell #phi;TDC counts",
				    145, 0., 2305., 128, 0., 1024.);
  hb_t2_phi_digi[1] = btl.make<TH2F>("h_t2_phi_digi_1", "BTL DIGI time2 vs #phi (R);cell #phi;TDC counts",
				    145, 0., 2305., 128, 0., 1024.);

  pb_t1_e_digi[0]   = btl.make<TProfile>("p_t1_e_digi_0", "BTL DIGI time1 vs charge (L);ADC counts;TDC counts",
					128, 0., 1024.);
  pb_t1_e_digi[1]   = btl.make<TProfile>("p_t1_e_digi_1", "BTL DIGI time1 vs charge (R);ADC counts;TDC counts",
					128, 0., 1024.);
  pb_t2_e_digi[0]   = btl.make<TProfile>("p_t2_e_digi_0", "BTL DIGI time2 vs charge (L);ADC counts;TDC counts",
					128, 0., 1024.);
  pb_t2_e_digi[1]   = btl.make<TProfile>("p_t2_e_digi_1", "BTL DIGI time2 vs charge (R);ADC counts;TDC counts",
					128, 0., 1024.);
  pb_e_eta_digi[0]  = btl.make<TProfile>("p_e_eta_digi_0", "BTL DIGI charge vs |#eta| (L);cell |#eta|;ADC counts",
					43, 0., 43.);
  pb_e_eta_digi[1]  = btl.make<TProfile>("p_e_eta_digi_1", "BTL DIGI charge vs |#eta| (R);cell |#eta|;ADC counts",
					43, 0., 43.);
  pb_t1_eta_digi[0] = btl.make<TProfile>("p_t1_eta_digi_0", "BTL DIGI time1 vs |#eta| (L);cell |#eta|;TDC counts",
					43, 0., 43.);
  pb_t1_eta_digi[1] = btl.make<TProfile>("p_t1_eta_digi_1", "BTL DIGI time1 vs |#eta| (R);cell |#eta|;TDC counts",
					43, 0., 43.);
  pb_t2_eta_digi[0] = btl.make<TProfile>("p_t2_eta_digi_0", "BTL DIGI time2 vs |#eta| (L);cell |#eta|;TDC counts",
					43, 0., 43.);
  pb_t2_eta_digi[1] = btl.make<TProfile>("p_t2_eta_digi_1", "BTL DIGI time2 vs |#eta| (R);cell |#eta|;TDC counts",
					43, 0., 43.);
  pb_e_phi_digi[0]  = btl.make<TProfile>("p_e_phi_digi_0", "BTL DIGI charge vs #phi (L);cell #phi;ADC counts",
					145, 0., 2305.);
  pb_e_phi_digi[1]  = btl.make<TProfile>("p_e_phi_digi_1", "BTL DIGI charge vs #phi (R);cell #phi;ADC counts",
					145, 0., 2305.);
  pb_t1_phi_digi[0] = btl.make<TProfile>("p_t1_phi_digi_0", "BTL DIGI time1 vs #phi (L);cell #phi;TDC counts",
					145, 0., 2305.);
  pb_t1_phi_digi[1] = btl.make<TProfile>("p_t1_phi_digi_1", "BTL DIGI time1 vs #phi (R);cell #phi;TDC counts",
					145, 0., 2305.);
  pb_t2_phi_digi[0] = btl.make<TProfile>("p_t2_phi_digi_0", "BTL DIGI time2 vs #phi (L);cell #phi;TDC counts",
					145, 0., 2305.);
  pb_t2_phi_digi[1] = btl.make<TProfile>("p_t2_phi_digi_1", "BTL DIGI time2 vs #phi (R);cell #phi;TDC counts",
					145, 0., 2305.);


  // --- Uncalibrated RECO

  hb_n_ureco[0]  = btl.make<TH1F>("h_n_ureco_0", "Number of BTL URECO hits (L);N_{URECO hits}", 100, 0., 100.);
  hb_n_ureco[1]  = btl.make<TH1F>("h_n_ureco_1", "Number of BTL URECO hits (R);N_{URECO hits}", 100, 0., 100.);
  hb_occupancy_ureco[0] = btl.make<TH2F>("h_occupancy_ureco_0", "BTL URECO hits occupancy (L);cell #phi;cell #eta",
					145, 0., 2305., 86, -43., 43.);
  hb_occupancy_ureco[1] = btl.make<TH2F>("h_occupancy_ureco_1", "BTL URECO hits occupancy (R);cell #phi;cell #eta",
					145, 0., 2305., 86, -43., 43.);
  hb_t_ureco[0] = btl.make<TH1F>("h_t_ureco_0", "BTL URECO hits ToA (L);ToA [ns]", 250, 0., 25.);
  hb_t_ureco[1] = btl.make<TH1F>("h_t_ureco_1", "BTL URECO hits ToA (R);ToA [ns]", 250, 0., 25.);
  hb_t_ureco_uncorr[0] = btl.make<TH1F>("h_t_ureco_uncorr_0", "BTL URECO hits ToA (L);ToA [ns]", 250, 0., 25.);
  hb_t_ureco_uncorr[1] = btl.make<TH1F>("h_t_ureco_uncorr_1", "BTL URECO hits ToA (R);ToA [ns]", 250, 0., 25.);
  hb_e_ureco[0] = btl.make<TH1F>("h_e_ureco_0", "BTL URECO hits energy (L);Q [pC]", 300, 0., 600.);
  hb_e_ureco[1] = btl.make<TH1F>("h_e_ureco_1", "BTL URECO hits energy (R);Q [pC]", 300, 0., 600.);

  hb_t_amp_ureco[0] = btl.make<TH2F>("h_t_amp_ureco_0", "time vs amplitude (L);amplitude [pC];time [ns]",
				    100, 0., 600., 400, 0., 20.);
  hb_t_amp_ureco[1] = btl.make<TH2F>("h_t_amp_ureco_1", "time vs amplitude (R);amplitude [pC];time [ns]",
				    100, 0., 600., 400, 0., 20.);
  pb_t_amp_ureco[0] = btl.make<TProfile>("p_t_amp_ureco_0", "time vs amplitude (L);amplitude [pC];time [ns]",
					100, 0., 600.);
  pb_t_amp_ureco[1] = btl.make<TProfile>("p_t_amp_ureco_1", "time vs amplitude (R);amplitude [pC];time [ns]",
					100, 0., 600.);


  // --- RECO

  hb_n_reco  = btl.make<TH1F>("h_n_reco", "Number of BTL RECO hits;N_{RECO hits}", 100, 0., 100.);
  hb_occupancy_reco = btl.make<TH2F>("h_occupancy_reco", "BTL RECO hits occupancy;cell #phi;cell #eta",
				    145, 0., 2305., 86, -43., 43.);
  hb_t_reco  = btl.make<TH1F>("h_t_reco", "BTL RECO hits ToA;ToA [ns]", 250, 0., 25.);
  hb_t_reco_uncorr = btl.make<TH1F>("h_t_reco_uncorr", "BTL RECO hits ToA;ToA [ns]", 250, 0., 25.);
  hb_e_reco  = btl.make<TH1F>("h_e_reco", "BTL RECO hits energy;E [MeV]", 200, 0., 20.);

  hb_t_res  = btl.make<TH1F>("h_t_res", "ToA resolution;ToA [ns]", 700, -2., 5.);
  hb_t_res_uncorr = btl.make<TH1F>("h_t_res_uncorr", "ToA resolution;ToA [ns]", 700, -2., 5.);
  hb_e_res  = btl.make<TH1F>("h_e_res", "Energy resolution;E [MeV]", 200, -1., 1.);

  hb_t_reco_sim = btl.make<TH2F>("h_t_reco_sim", "ToA reco vs sim;SIM ToA [ns];BTL RECO ToA [ns]",
				100, -1., 25., 100, 0., 25.);
  hb_e_reco_sim = btl.make<TH2F>("h_e_reco_sim", "E reco vs sim;SIM E [MeV];BTL RECO E [MeV]",
				100, 0., 20., 100, 0., 20.);


  // ==============================================================================
  //  ETL
  // ==============================================================================

  // --- SIM

  he_n_sim_trk[0]   = etl.make<TH1F>("h_n_sim_trk_0", "Number of tracks per ETL cell (-Z);N_{trk}", 10, 0., 10.);
  he_n_sim_trk[1]   = etl.make<TH1F>("h_n_sim_trk_1", "Number of tracks per ETL cell (+Z);N_{trk}", 10, 0., 10.);
  he_n_sim_cell[0] = etl.make<TH1F>("h_n_sim_cell_0", "Number of ETL cells with SIM hits (-Z);N_{ETL cells}", 500, 0., 1000.);
  he_n_sim_cell[1] = etl.make<TH1F>("h_n_sim_cell_1", "Number of ETL cells with SIM hits (+Z);N_{ETL cells}", 500, 0., 1000.);
  he_t_sim[0]   = etl.make<TH1F>("h_t_sim_0", "ETL SIM hits ToA (-Z);ToA_{SIM} [ns]", 250, 0., 25.);
  he_t_sim[1]   = etl.make<TH1F>("h_t_sim_1", "ETL SIM hits ToA (+Z);ToA_{SIM} [ns]", 250, 0., 25.);
  he_e_sim[0]   = etl.make<TH1F>("h_e_sim_0", "ETL SIM hits energy (-Z);E_{SIM} [MIP]", 200, 0., 1.);
  he_e_sim[1]   = etl.make<TH1F>("h_e_sim_1", "ETL SIM hits energy (+Z);E_{SIM} [MIP]", 200, 0., 1.);

  he_xloc_sim[0] = etl.make<TH1F>("h_xloc_sim_0", "ETL SIM local x (-Z);x_{SIM} [mm]", 100, -25., 25.);
  he_xloc_sim[1] = etl.make<TH1F>("h_xloc_sim_1", "ETL SIM local x (+Z);x_{SIM} [mm]", 100, -25., 25.);
  he_yloc_sim[0] = etl.make<TH1F>("h_yloc_sim_0", "ETL SIM local y (-Z);y_{SIM} [mm]", 200, -50., 50.);
  he_yloc_sim[1] = etl.make<TH1F>("h_yloc_sim_1", "ETL SIM local y (+Z);y_{SIM} [mm]", 200, -50., 50.);
  he_zloc_sim[0] = etl.make<TH1F>("h_zloc_sim_0", "ETL SIM local z (-Z);z_{SIM} [mm]", 80, -0.2, 0.2);
  he_zloc_sim[1] = etl.make<TH1F>("h_zloc_sim_1", "ETL SIM local z (+Z);z_{SIM} [mm]", 80, -0.2, 0.2);

  he_occupancy_sim[0] = etl.make<TH2F>("h_occupancy_sim_0", "ETL SIM hits occupancy (-Z);x_{SIM} [cm];y_{SIM} [cm]",
				       135, -135., 135.,  135, -135., 135.);
  he_occupancy_sim[1] = etl.make<TH2F>("h_occupancy_sim_1", "ETL SIM hits occupancy (+Z);x_{SIM} [cm];y_{SIM} [cm]",
				       135, -135., 135.,  135, -135., 135.);
  he_x_sim[0] = etl.make<TH1F>("h_x_sim_0", "ETL SIM hits x (-Z);x_{SIM} [cm]", 135, -135., 135.);
  he_x_sim[1] = etl.make<TH1F>("h_x_sim_1", "ETL SIM hits x (+Z);x_{SIM} [cm]", 135, -135., 135.);
  he_y_sim[0] = etl.make<TH1F>("h_y_sim_0", "ETL SIM hits y (-Z);y_{SIM} [cm]", 135, -135., 135.);
  he_y_sim[1] = etl.make<TH1F>("h_y_sim_1", "ETL SIM hits y (+Z);y_{SIM} [cm]", 135, -135., 135.);
  he_z_sim[0] = etl.make<TH1F>("h_z_sim_0", "ETL SIM hits z (-Z);z_{SIM} [cm]", 100, -304.5, -303.);
  he_z_sim[1] = etl.make<TH1F>("h_z_sim_1", "ETL SIM hits z (+Z);z_{SIM} [cm]", 100,  303., 304.5);
  he_phi_sim[0] = etl.make<TH1F>("h_phi_sim_0", "ETL SIM hits #phi (-Z);#phi_{SIM} [rad]", 315, -3.15, 3.15);
  he_phi_sim[1] = etl.make<TH1F>("h_phi_sim_1", "ETL SIM hits #phi (+Z);#phi_{SIM} [rad]", 315, -3.15, 3.15);
  he_eta_sim[0] = etl.make<TH1F>("h_eta_sim_0", "ETL SIM hits #eta (-Z);#eta_{SIM}", 200, -3.05, -1.55);
  he_eta_sim[1] = etl.make<TH1F>("h_eta_sim_1", "ETL SIM hits #eta (+Z);#eta_{SIM}", 200,  1.55, 3.05);

  he_t_e_sim[0]   = etl.make<TH2F>("h_t_e_sim_0", "ETL SIM time vs energy (-Z);E_{SIM} [MIP];T_{SIM} [ns]",
				   100, 0., 2., 100, 0., 25.);
  he_t_e_sim[1]   = etl.make<TH2F>("h_t_e_sim_1", "ETL SIM time vs energy (+Z);E_{SIM} [MIP];T_{SIM} [ns]",
				   100, 0., 2., 100, 0., 25.);
  he_e_eta_sim[0] = etl.make<TH2F>("h_e_eta_sim_0", "ETL SIM energy vs #eta (-Z);#eta_{SIM};E_{SIM} [MIP]",
				   100, -3.05, -1.55, 100, 0., 2.);
  he_e_eta_sim[1] = etl.make<TH2F>("h_e_eta_sim_1", "ETL SIM energy vs #eta (+Z);#eta_{SIM};E_{SIM} [MIP]",
				   100, 1.55, 3.05, 100, 0., 2.);
  he_t_eta_sim[0] = etl.make<TH2F>("h_t_eta_sim_0", "ETL SIM time vs #eta (-Z);#eta_{SIM};T_{SIM} [ns]",
				   100, -3.05, -1.55, 100, 0., 25.);
  he_t_eta_sim[1] = etl.make<TH2F>("h_t_eta_sim_1", "ETL SIM time vs #eta (+Z);#eta_{SIM};T_{SIM} [ns]",
				   100, 1.55, 3.05, 100, 0., 25.);
  he_e_phi_sim[0] = etl.make<TH2F>("h_e_phi_sim_0", "ETL SIM energy vs #phi (-Z);#phi_{SIM} [rad];E_{SIM} [MIP]",
				   100, -3.15, 3.15, 100, 0., 2.);
  he_e_phi_sim[1] = etl.make<TH2F>("h_e_phi_sim_1", "ETL SIM energy vs #phi (+Z);#phi_{SIM} [rad];E_{SIM} [MIP]",
				   100, -3.15, 3.15, 100, 0., 2.);
  he_t_phi_sim[0] = etl.make<TH2F>("h_t_phi_sim_0", "ETL SIM time vs #phi (-Z);#phi_{SIM} [rad];T_{SIM} [ns]",
				   100, -3.15, 3.15, 100, 0., 25.);
  he_t_phi_sim[1] = etl.make<TH2F>("h_t_phi_sim_1", "ETL SIM time vs #phi (+Z);#phi_{SIM} [rad];T_{SIM} [ns]",
				   100, -3.15, 3.15, 100, 0., 25.);
  pe_t_e_sim[0]   = etl.make<TProfile>("p_t_e_sim_0", "ETL SIM time vs energy (-Z);E_{SIM} [MIP];T_{SIM} [ns]",
				       100, 0., 2.);
  pe_t_e_sim[1]   = etl.make<TProfile>("p_t_e_sim_1", "ETL SIM time vs energy (+Z);E_{SIM} [MIP];T_{SIM} [ns]",
				       100, 0., 2.);
  pe_e_eta_sim[0] = etl.make<TProfile>("p_e_eta_sim_0", "ETL SIM energy vs #eta (-Z);#eta_{SIM};E_{SIM} [MIP]",
				       100, -3.05, -1.55);
  pe_e_eta_sim[1] = etl.make<TProfile>("p_e_eta_sim_1", "ETL SIM energy vs #eta (+Z);#eta_{SIM};E_{SIM} [MIP]",
				       100, 1.55, 3.05);
  pe_t_eta_sim[0] = etl.make<TProfile>("p_t_eta_sim_0", "ETL SIM time vs #eta (-Z);#eta_{SIM};T_{SIM} [ns]",
				       100, -3.05, -1.55);
  pe_t_eta_sim[1] = etl.make<TProfile>("p_t_eta_sim_0", "ETL SIM time vs #eta (+Z);#eta_{SIM};T_{SIM} [ns]",
				       100, 1.55, 3.05);
  pe_e_phi_sim[0] = etl.make<TProfile>("p_e_phi_sim_0", "ETL SIM energy vs #phi (-Z);#phi_{SIM} [rad];E_{SIM} [MIP]",
				       100, -3.15, 3.15);
  pe_e_phi_sim[1] = etl.make<TProfile>("p_e_phi_sim_1", "ETL SIM energy vs #phi (+Z);#phi_{SIM} [rad];E_{SIM} [MIP]",
				       100, -3.15, 3.15);
  pe_t_phi_sim[0] = etl.make<TProfile>("p_t_phi_sim_0", "ETL SIM time vs #phi (-Z);#phi_{SIM} [rad];T_{SIM} [ns]",
				       100, -3.15, 3.15);
  pe_t_phi_sim[1] = etl.make<TProfile>("p_t_phi_sim_1", "ETL SIM time vs #phi (+Z);#phi_{SIM} [rad];T_{SIM} [ns]",
				       100, -3.15, 3.15);



  // --- DIGI

  he_n_digi[0]  = etl.make<TH1F>("h_n_digi_0", "Number of ETL DIGI hits (-Z);N_{DIGI hits}", 100, 0., 100.);
  he_n_digi[1]  = etl.make<TH1F>("h_n_digi_1", "Number of ETL DIGI hits (+Z);N_{DIGI hits}", 100, 0., 100.);

  he_t_digi[0]  = etl.make<TH1F>("h_t_digi_0", "ETL DIGI hits ToA (-Z);ToA [TDC counts]", 1024, 0., 4048.);
  he_t_digi[1]  = etl.make<TH1F>("h_t_digi_1", "ETL DIGI hits ToA (+Z);ToA [TDC counts]", 1024, 0., 4048.);
  he_e_digi[0]  = etl.make<TH1F>("h_e_digi_0", "ETL DIGI hits energy (-Z);amplitude [ADC counts]", 256, 0., 256.);
  he_e_digi[1]  = etl.make<TH1F>("h_e_digi_1", "ETL DIGI hits energy (+Z);amplitude [ADC counts]", 256, 0., 256.);


  // --- Uncalibrated RECO

  he_n_ureco[0]  = etl.make<TH1F>("h_n_ureco_0", "Number of ETL URECO hits (-Z);N_{URECO hits}", 100, 0., 100.);
  he_n_ureco[1]  = etl.make<TH1F>("h_n_ureco_1", "Number of ETL URECO hits (+Z);N_{URECO hits}", 100, 0., 100.);


  // --- RECO

  he_n_reco[0]  = etl.make<TH1F>("h_n_reco_0", "Number of ETL RECO hits (-Z);N_{RECO hits}", 100, 0., 100.);
  he_n_reco[1]  = etl.make<TH1F>("h_n_reco_1", "Number of ETL RECO hits (+Z);N_{RECO hits}", 100, 0., 100.);




}


MTDAnalyzer::~MTDAnalyzer() {}


//
// member functions
//

void
MTDAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  using namespace std;

  edm::ESHandle<MTDGeometry> geom;
  if( geom_ == nullptr ) {
    iSetup.get<MTDDigiGeometryRecord>().get(geom);
    geom_ = geom.product();
  }

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


  ///////////////////////////////////////////////////////////////////////////////////////////////
  //
  //  Get MTD hits
  //
  ///////////////////////////////////////////////////////////////////////////////////////////////

  // ==============================================================================
  //  SIM hits
  // ==============================================================================

  // --- BTL

  std::set<uint32_t> unique_simHit;

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

      if ( hit.tof() < btlIntegrationWindow_ ) // This is to emulate the time integration
	                                       // window in the readout electronics.
	(simHitIt->second).sim_energy += 1000.*hit.energyLoss();

      n_btl_simHits[id.rawId()].insert(hit.trackId());

      // Get the time of the first SimHit in the cell
      if( (simHitIt->second).sim_time==0 ) {

	//auto hit_pos = hit.localPosition();
	auto hit_pos = hit.entryPoint();

	(simHitIt->second).sim_x = hit_pos.x();
	(simHitIt->second).sim_y = hit_pos.y();
	(simHitIt->second).sim_z = hit_pos.z();

	(simHitIt->second).sim_time = hit.tof();

      }


    } // hitRef loop

  } // if ( h_BTL_sim->size() > 0 )


  // --- ETL

  std::set<uint32_t> unique_etl_simHit[2];

  if ( h_ETL_sim->size() > 0 ) {

    std::vector< std::tuple<PSimHit,uint32_t,float> > hitRefs;
    hitRefs.reserve(h_ETL_sim->size());

    // Sort the SimHits per detector id and time
    for (auto const& simHit: *h_ETL_sim) {

      DetId id = simHit.detUnitId();
      if (id.rawId() != 0)
	hitRefs.emplace_back( simHit, id.rawId(), simHit.tof() );

    } // simHit loop
    std::sort(hitRefs.begin(),hitRefs.end(), orderByDetIdThenTime);

    // Accumulate the SimHits in the same detector cell
    for (auto const& hitRef: hitRefs) {

      const PSimHit &hit = std::get<0>(hitRef);
      ETLDetId id = hit.detUnitId();

      int idet = (id.zside()+1)/2;

      unique_etl_simHit[idet].insert(id.rawId());

      auto simHitIt = etl_hits[idet].emplace(id.rawId(),MTDinfo()).first;

      (simHitIt->second).sim_energy += 1000.*hit.energyLoss();

      n_etl_simHits[idet][id.rawId()].insert(hit.trackId());

      // Get the time of the first SimHit in the cell
      if( (simHitIt->second).sim_time==0 ) {

	(simHitIt->second).sim_time = hit.tof();

	//auto hit_pos = hit.localPosition();
	auto hit_pos = hit.entryPoint();

	(simHitIt->second).sim_x = hit_pos.x();
	(simHitIt->second).sim_y = hit_pos.y();
	(simHitIt->second).sim_z = hit_pos.z();

      }

    } // hitRef loop

  } // if ( h_ETL_sim->size() > 0 )


  // ==============================================================================
  //  DIGI hits
  // ==============================================================================

  // --- BTL

  unsigned int n_digi_btl[2] = {0,0};

  if ( h_BTL_digi->size() > 0 ) {

    for (const auto& dataFrame: *h_BTL_digi) {

      DetId id =  dataFrame.id();
      
      const auto& sample_L = dataFrame.sample(0);
      const auto& sample_R = dataFrame.sample(1);

      btl_hits[id.rawId()].digi_charge[0] = sample_L.data();
      btl_hits[id.rawId()].digi_charge[1] = sample_R.data();
      btl_hits[id.rawId()].digi_time1[0]  = sample_L.toa();
      btl_hits[id.rawId()].digi_time1[1]  = sample_R.toa();
      btl_hits[id.rawId()].digi_time2[0]  = sample_L.toa2();
      btl_hits[id.rawId()].digi_time2[1]  = sample_R.toa2();

      if ( sample_L.data() > 0 )
	n_digi_btl[0]++;

      if ( sample_R.data() > 0 )
	n_digi_btl[1]++;

    } // dataFrame loop

  } // if ( h_BTL_digi->size() > 0 )


  // --- ETL

  unsigned int n_digi_etl[2] = {0,0};

  if ( h_ETL_digi->size() > 0 ) {

    for (const auto& dataFrame: *h_ETL_digi) {

      ETLDetId id =  dataFrame.id();
      int idet = (id.zside()+1)/2;

      // --- loop over the dataFrame samples
      for (int isample = 0; isample<dataFrame.size(); ++isample){

	const auto& sample = dataFrame.sample(isample);

	if ( sample.data()!=0 && sample.toa()!=0 ) {

	  // on-time sample
	  if ( isample == 2 ) {

	    etl_hits[idet][id.rawId()].digi_charge[0] =  sample.data();
	    etl_hits[idet][id.rawId()].digi_time1[0]  =  sample.toa();

	    n_digi_etl[idet]++;

	  }

	}

      } // isaple loop

    } // dataFrame loop

  } // if ( h_ETL_digi->size() > 0 )


  // ==============================================================================
  //  Uncalibrated RECO hits
  // ==============================================================================

  // --- BTL

  unsigned int n_ureco_btl[2] = {0,0};

  if ( h_BTL_ureco->size() > 0 ) {

    for (const auto& urecHit: *h_BTL_ureco) {

      DetId id = urecHit.id();

      btl_hits[id.rawId()].ureco_charge[0] = urecHit.amplitude().first;
      btl_hits[id.rawId()].ureco_charge[1] = urecHit.amplitude().second;
      btl_hits[id.rawId()].ureco_time[0]   = urecHit.time().first;
      btl_hits[id.rawId()].ureco_time[1]   = urecHit.time().second;

      if ( urecHit.amplitude().first > 0. )
	n_ureco_btl[0]++;

      if ( urecHit.amplitude().second > 0. )
	n_ureco_btl[1]++;

    } // recHit loop

  } // if ( h_BTL_reco->size() > 0 )


  // --- ETL

  unsigned int n_ureco_etl[2] = {0,0};

  if ( h_ETL_ureco->size() > 0 ) {

    for (const auto& urecHit: *h_ETL_ureco) {

      ETLDetId id = urecHit.id();
      int idet = (id.zside()+1)/2;

      etl_hits[idet][id.rawId()].ureco_charge[0] = urecHit.amplitude().first;
      etl_hits[idet][id.rawId()].ureco_time[0]   = urecHit.time().first;

      if ( urecHit.amplitude().first > 0. )
	n_ureco_etl[idet]++;
      

    } // recHit loop

  } // if ( h_ETL_ureco->size() > 0 )


  // ==============================================================================
  //  RECO hits
  // ==============================================================================

  // --- BTL

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


  // --- ETL

  unsigned int n_reco_etl[2] = {0,0};

  if ( h_ETL_reco->size() > 0 ) {

    for (const auto& recHit: *h_ETL_reco) {

      ETLDetId id = recHit.id();
      int idet = (id.zside()+1)/2;

      etl_hits[idet][id.rawId()].reco_energy = recHit.energy();
      etl_hits[idet][id.rawId()].reco_time   = recHit.time();

      if ( recHit.energy() > 0. )
	n_reco_etl[idet]++;


    } // recHit loop

  } // if ( h_ETL_reco->size() > 0 )


  ///////////////////////////////////////////////////////////////////////////////////////////////
  //
  //  Histograms filling
  //
  ///////////////////////////////////////////////////////////////////////////////////////////////


  // ==============================================================================
  //  BTL
  // ==============================================================================

  for (auto const& hit: n_btl_simHits) {
    hb_n_sim_trk->Fill((hit.second).size());
  }
  hb_n_sim_cell->Fill(unique_simHit.size());
  for (int iside=0; iside<2; ++iside){
    hb_n_digi[iside]->Fill(n_digi_btl[iside]);
    hb_n_ureco[iside]->Fill(n_ureco_btl[iside]);
  }
  hb_n_reco->Fill(n_reco_btl);

  for (auto const& hit: btl_hits) {

    BTLDetId detId(hit.first); 

    if ( (hit.second).reco_energy < btlMinEnergy_ ) continue;
    

    // --- SIM

    hb_e_sim->Fill((hit.second).sim_energy);
    hb_t_sim->Fill((hit.second).sim_time);

    hb_xloc_sim->Fill((hit.second).sim_x);
    hb_yloc_sim->Fill((hit.second).sim_y);
    hb_zloc_sim->Fill((hit.second).sim_z);


    // Get the SIM hit global position
    DetId geoId = BTLDetId(detId.mtdSide(),detId.mtdRR(),detId.module()+14*(detId.modType()-1),0,1);
    const MTDGeomDet* thedet = geom_->idToDet(geoId);
    const ProxyMTDTopology& topoproxy = static_cast<const ProxyMTDTopology&>(thedet->topology());
    const RectangularMTDTopology& topo = static_cast<const RectangularMTDTopology&>(topoproxy.specificTopology()); 
    
    Local3DPoint simscaled(0.1*(hit.second).sim_x,0.1*(hit.second).sim_y,0.1*(hit.second).sim_z);
    simscaled = topo.pixelToModuleLocalPoint(simscaled,detId.row(topo.nrows()),detId.column(topo.nrows()));
    const auto& global_pos = thedet->toGlobal(simscaled);

    hb_occupancy_sim->Fill(global_pos.z(),global_pos.phi());
    hb_phi_sim->Fill(global_pos.phi());
    hb_eta_sim->Fill(global_pos.eta());
    hb_z_sim->Fill(global_pos.z());
  
    hb_t_e_sim->Fill((hit.second).sim_energy,(hit.second).sim_time);
    hb_e_eta_sim->Fill(fabs(global_pos.eta()),(hit.second).sim_energy);
    hb_t_eta_sim->Fill(fabs(global_pos.eta()),(hit.second).sim_time);
    hb_e_phi_sim->Fill(global_pos.phi(),(hit.second).sim_energy);
    hb_t_phi_sim->Fill(global_pos.phi(),(hit.second).sim_time);

    pb_t_e_sim->Fill((hit.second).sim_energy,(hit.second).sim_time);
    pb_e_eta_sim->Fill(fabs(global_pos.eta()),(hit.second).sim_energy);
    pb_t_eta_sim->Fill(fabs(global_pos.eta()),(hit.second).sim_time);
    pb_e_phi_sim->Fill(global_pos.phi(),(hit.second).sim_energy);
    pb_t_phi_sim->Fill(global_pos.phi(),(hit.second).sim_time);




    int hit_iphi = detId.iphi(BTLDetId::CrysLayout::barzflat);
    int hit_ieta = detId.ieta(BTLDetId::CrysLayout::barzflat);



    // Time-walk correction:
    const float p0 =  2.21103;
    const float p1 = -0.933552;
    const float p2 =  0.;
    float time_corr[2];

    for (int iside=0; iside<2; ++iside){

      // --- DIGI

      if ( (hit.second).digi_charge[iside] == 0 ) continue;

      hb_e_digi[iside] ->Fill((hit.second).digi_charge[iside]);
      hb_t1_digi[iside]->Fill((hit.second).digi_time1[iside]);
      hb_t2_digi[iside]->Fill((hit.second).digi_time2[iside]);

      hb_occupancy_digi[iside]->Fill(hit_iphi,hit_ieta);
      hb_phi_digi[iside]->Fill(hit_iphi);
      hb_eta_digi[iside]->Fill(hit_ieta);

      hb_t1_e_digi[iside]->Fill((hit.second).digi_charge[iside], (hit.second).digi_time1[iside]);
      hb_t2_e_digi[iside]->Fill((hit.second).digi_charge[iside], (hit.second).digi_time2[iside]);
      hb_e_eta_digi[iside]->Fill(fabs(hit_ieta),(hit.second).digi_charge[iside]);
      hb_t1_eta_digi[iside]->Fill(fabs(hit_ieta),(hit.second).digi_time1[iside]);
      hb_t2_eta_digi[iside]->Fill(fabs(hit_ieta),(hit.second).digi_time2[iside]);
      hb_e_phi_digi[iside]->Fill(hit_iphi,(hit.second).digi_charge[iside]);
      hb_t1_phi_digi[iside]->Fill(hit_iphi,(hit.second).digi_time1[iside]);
      hb_t2_phi_digi[iside]->Fill(hit_iphi,(hit.second).digi_time2[iside]);

      pb_t1_e_digi[iside]->Fill((hit.second).digi_charge[iside], (hit.second).digi_time1[iside]);
      pb_t2_e_digi[iside]->Fill((hit.second).digi_charge[iside], (hit.second).digi_time2[iside]);
      pb_e_eta_digi[iside]->Fill(fabs(hit_ieta),(hit.second).digi_charge[iside]);
      pb_t1_eta_digi[iside]->Fill(fabs(hit_ieta),(hit.second).digi_time1[iside]);
      pb_t2_eta_digi[iside]->Fill(fabs(hit_ieta),(hit.second).digi_time2[iside]);
      pb_e_phi_digi[iside]->Fill(hit_iphi,(hit.second).digi_charge[iside]);
      pb_t1_phi_digi[iside]->Fill(hit_iphi,(hit.second).digi_time1[iside]);
      pb_t2_phi_digi[iside]->Fill(hit_iphi,(hit.second).digi_time2[iside]);


      // --- Uncalibrated RECO

      if ( (hit.second).ureco_charge[iside] == 0. ) continue;

      hb_e_ureco[iside]->Fill((hit.second).ureco_charge[iside]);
      hb_t_ureco[iside]->Fill((hit.second).ureco_time[iside]);

      hb_occupancy_ureco[iside]->Fill(hit_iphi,hit_ieta);

      hb_t_amp_ureco[iside]->Fill((hit.second).ureco_charge[iside],(hit.second).ureco_time[iside]);
      pb_t_amp_ureco[iside]->Fill((hit.second).ureco_charge[iside],(hit.second).ureco_time[iside]);
    

      // Reverse the time-walk correction

      time_corr[iside] = p0*pow((hit.second).ureco_charge[iside],p1) + p2;

      float ureco_time_uncorr = (hit.second).ureco_time[iside] + time_corr[iside];

      hb_t_ureco_uncorr[iside]->Fill(ureco_time_uncorr);


    } // for iside


    // --- RECO

    if ( (hit.second).reco_energy == 0. ) continue;

    hb_occupancy_reco->Fill(hit_iphi,hit_ieta);

    hb_e_reco->Fill((hit.second).reco_energy);
    hb_t_reco->Fill((hit.second).reco_time);
    
    hb_e_res->Fill((hit.second).reco_energy-(hit.second).sim_energy);
    hb_t_res->Fill((hit.second).reco_time-(hit.second).sim_time);
    
    hb_t_reco_sim->Fill((hit.second).sim_time,(hit.second).reco_time);
    hb_e_reco_sim->Fill((hit.second).sim_energy,(hit.second).reco_energy);

    
    float reco_time_uncorr = (hit.second).reco_time + 0.5*(time_corr[0]+time_corr[1]); 

    hb_t_reco_uncorr->Fill(reco_time_uncorr);
    hb_t_res_uncorr->Fill(reco_time_uncorr-(hit.second).sim_time);


 } // BTL hit loop


  // ==============================================================================
  //  ETL
  // ==============================================================================

  for (int idet=0; idet<2; ++idet){

    for (auto const& hit: n_etl_simHits[idet]) {
      he_n_sim_trk[idet]->Fill((hit.second).size());
    }

    he_n_sim_cell[idet]->Fill(unique_etl_simHit[idet].size());
    he_n_digi[idet]->Fill(n_digi_etl[idet]);
    he_n_ureco[idet]->Fill(n_ureco_etl[idet]);
    he_n_reco[idet]->Fill(n_reco_etl[idet]);


    for (auto const& hit: etl_hits[idet]) {

      ETLDetId detId(hit.first);

      // --- SIM

      if ( (hit.second).sim_energy == 0. ) continue;

      he_e_sim[idet]->Fill((hit.second).sim_energy);
      he_t_sim[idet]->Fill((hit.second).sim_time);
      
      he_xloc_sim[idet]->Fill((hit.second).sim_x);
      he_yloc_sim[idet]->Fill((hit.second).sim_y);
      he_zloc_sim[idet]->Fill((hit.second).sim_z);


      // Get the SIM hit global position
      DetId geoId = ETLDetId(detId.mtdSide(),detId.mtdRR(),detId.module(),0);
      const MTDGeomDet* thedet = geom_->idToDet(geoId);
      Local3DPoint simscaled(0.1*(hit.second).sim_x,0.1*(hit.second).sim_y,0.1*(hit.second).sim_z);
      const auto& global_pos = thedet->toGlobal(simscaled);

      he_occupancy_sim[idet]->Fill(global_pos.x(),global_pos.y());
      he_x_sim[idet]->Fill(global_pos.x());
      he_y_sim[idet]->Fill(global_pos.y());
      he_z_sim[idet]->Fill(global_pos.z());
      he_phi_sim[idet]->Fill(global_pos.phi());
      he_eta_sim[idet]->Fill(global_pos.eta());
  
      he_t_e_sim[idet]->Fill((hit.second).sim_energy,(hit.second).sim_time);
      he_e_eta_sim[idet]->Fill(global_pos.eta(),(hit.second).sim_energy);
      he_t_eta_sim[idet]->Fill(global_pos.eta(),(hit.second).sim_time);
      he_e_phi_sim[idet]->Fill(global_pos.phi(),(hit.second).sim_energy);
      he_t_phi_sim[idet]->Fill(global_pos.phi(),(hit.second).sim_time);

      pe_t_e_sim[idet]->Fill((hit.second).sim_energy,(hit.second).sim_time);
      pe_e_eta_sim[idet]->Fill(global_pos.eta(),(hit.second).sim_energy);
      pe_t_eta_sim[idet]->Fill(global_pos.eta(),(hit.second).sim_time);
      pe_e_phi_sim[idet]->Fill(global_pos.phi(),(hit.second).sim_energy);
      pe_t_phi_sim[idet]->Fill(global_pos.phi(),(hit.second).sim_time);



      // --- DIGI

      if ( (hit.second).digi_charge[0] == 0 ) continue;

      he_e_digi[idet]->Fill((hit.second).digi_charge[0]);
      he_t_digi[idet]->Fill((hit.second).digi_time1[0]);



    } // ETL hit loop

  } // idet loop

  // ---------------------------------------------------------------

  n_btl_simHits.clear();
  n_etl_simHits[0].clear();
  n_etl_simHits[1].clear();
  btl_hits.clear();
  etl_hits[0].clear();
  etl_hits[1].clear();

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
