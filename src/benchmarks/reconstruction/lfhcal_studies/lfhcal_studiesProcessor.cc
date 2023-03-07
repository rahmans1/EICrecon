#include "lfhcal_studiesProcessor.h"
#include "algorithms/tracking/JugTrack/TrackingResultTrajectory.hpp"
#include "edm4eic/vector_utils.h"
#include <Acts/EventData/MultiTrajectoryHelpers.hpp>
#include <Acts/Surfaces/DiscSurface.hpp>
#include <Acts/Surfaces/RadialBounds.hpp>
#include <extensions/spdlog/SpdlogExtensions.h>
#include <services/rootfile/RootFile_service.h>
#include <spdlog/spdlog.h>

#include <extensions/spdlog/SpdlogExtensions.h>
#include <extensions/spdlog/SpdlogMixin.h>
#include <services/log/Log_service.h>
#include <spdlog/fmt/ostr.h>

#include "DD4hep/DetElement.h"
#include "DD4hep/Detector.h"
#include "DD4hep/Objects.h"
#include "DDG4/Geant4Data.h"
#include "DDRec/CellIDPositionConverter.h"
#include "DDRec/Surface.h"
#include "DDRec/SurfaceManager.h"
#include <JANA/JApplication.h>
#include <JANA/JEvent.h>

#include "TCanvas.h"
#include "TChain.h"
#include "TVector3.h"

#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Plugins/DD4hep/ConvertDD4hepDetector.hpp"

#include "clusterizer_MA.h"
// #include <extensions/spdlog/SpdlogMixin.h>

// The following just makes this a JANA plugin
extern "C" {
void InitPlugin(JApplication* app) {
  InitJANAPlugin(app);
  app->Add(new lfhcal_studiesProcessor());
}
}


//-------------------------------------------
// InitWithGlobalRootLock
//-------------------------------------------
void lfhcal_studiesProcessor::InitWithGlobalRootLock() {
  std::string plugin_name = ("lfhcal_studies");

  // InitLogger(plugin_name);
  // Get JANA application
  auto app          = GetApplication();
  auto acts_service = GetApplication()->GetService<ACTSGeo_service>();

  std::string log_level_str = "info";
  m_log                     = app->GetService<Log_service>()->logger(plugin_name);
  app->SetDefaultParameter(plugin_name + ":LogLevel", log_level_str,
                           "LogLevel: trace, debug, info, warn, err, critical, off");
  m_log->set_level(eicrecon::ParseLogLevel(log_level_str));

  m_geo_provider = acts_service->actsGeoProvider();
  m_propagation_algo.init(acts_service->actsGeoProvider(), m_log);

  // Ask service locator a file to write histograms to
  auto root_file_service = app->GetService<RootFile_service>();

  // Get TDirectory for histograms root file
  auto globalRootLock = app->GetService<JGlobalRootLock>();
  globalRootLock->acquire_write_lock();
  auto file = root_file_service->GetHistFile();
  globalRootLock->release_lock();

  // Create a directory for this plugin. And subdirectories for series of histograms
  m_dir_main = file->mkdir(plugin_name.c_str());

  hMCEnergyVsEta = new TH2D("hMCEnergyVsEta", "; E (GeV); #eta", 1500, 0., 150., 400, 1., 5.);
  hMCEnergyVsEta->SetDirectory(m_dir_main);

  // Sum cell clusters rec
  hClusterEcalib_E_eta = new TH3D("hClusterEcalib_E_eta", "; E_{MC} (GeV); E_{rec,rec hit}/E_{MC}; #eta", 1500, 0., 150.0, 200, 0., 2.0, 20, 1, 5);
  hClusterEcalib_E_eta->SetDirectory(m_dir_main);
  hClusterNCells_E_eta = new TH3D("hClusterNCells_E_eta", "; E_{MC} (GeV); N_{cells}; #eta",  1500, 0., 150.0, 500, -0.5, 499.5, 20, 1, 5);
  hClusterNCells_E_eta->SetDirectory(m_dir_main);
  // Sum cell clusters sim
  hClusterESimcalib_E_eta = new TH3D("hClusterESimcalib_E_eta", "; E_{MC} (GeV); E_{rec,sim hit}/E_{MC}; #eta" , 1500, 0., 150.0, 200, 0., 2.0, 20, 1, 5);
  hClusterESimcalib_E_eta->SetDirectory(m_dir_main);
  hClusterSimNCells_E_eta = new TH3D("hClusterSimNCells_E_eta", "; E_{MC} (GeV); N_{cells, sim}; #eta", 1500, 0., 150.0, 500, -0.5, 499.5, 20, 1, 5);
  hClusterSimNCells_E_eta->SetDirectory(m_dir_main);
  // rec cluster
  hRecClusterEcalib_E_eta = new TH3D("hRecClusterEcalib_E_eta", "; E_{MC} (GeV); E_{rec,rec clus}/E_{MC}; #eta", 1500, 0., 150.0, 200, 0., 2.0, 20, 1, 5);
  hRecClusterEcalib_E_eta->SetDirectory(m_dir_main);
  // rec cluster highest 
  hRecClusterEcalib_Ehigh_eta = new TH3D("hRecClusterEcalib_Ehigh_eta", "; E_{MC} (GeV); E_{rec,rec clus high.}/E_{MC}; #eta", 
                                         1500, 0., 150.0, 200, 0., 2.0, 20, 1, 5);
  hRecClusterEcalib_Ehigh_eta->SetDirectory(m_dir_main);
  hRecClusterNCells_Ehigh_eta = new TH3D("hRecClusterNCells_Ehigh_eta", "; E_{MC} (GeV); N_{cells, rec cl., high.}; #eta", 1500, 0., 150.0, 500, -0.5, 499.5, 20, 1, 5);
  hRecClusterNCells_Ehigh_eta->SetDirectory(m_dir_main);
  hRecNClusters_E_eta = new TH3D("hRecNClusters_E_eta", "; E_{MC} (GeV); N_{rec cl.}; #eta",  1500, 0., 150.0, 10, -0.5, 9.5, 20, 1, 5);
  hRecNClusters_E_eta->SetDirectory(m_dir_main);
  // rec cluster framework
  hRecFClusterEcalib_E_eta = new TH3D("hRecFClusterEcalib_E_eta", "; E_{MC} (GeV); E_{rec,fram clus}/E_{MC}; #eta", 1500, 0., 150.0, 200, 0., 2.0, 20, 1, 5);
  hRecFClusterEcalib_E_eta->SetDirectory(m_dir_main);
  // rec cluster framework highest 
  hRecFClusterEcalib_Ehigh_eta = new TH3D("hRecFClusterEcalib_Ehigh_eta", "; E_{MC} (GeV); E_{rec,fram clus high.}/E_{MC}; #eta", 
                                          1500, 0., 150.0, 200, 0., 2.0, 20, 1, 5);
  hRecFClusterEcalib_Ehigh_eta->SetDirectory(m_dir_main);
  hRecFClusterNCells_Ehigh_eta = new TH3D("hRecFClusterNCells_Ehigh_eta", "; E_{MC} (GeV); N_{cells, rec f. cl., high.}; #eta", 1500, 0., 150.0, 500, -0.5, 499.5, 20, 1, 5);
  hRecFClusterNCells_Ehigh_eta->SetDirectory(m_dir_main);
  hRecFNClusters_E_eta = new TH3D("hRecFNClusters_E_eta", "; E_{MC} (GeV); N_{rec f. cl.}; #eta",  1500, 0., 150.0, 10, -0.5, 9.5, 20, 1, 5);
  hRecFNClusters_E_eta->SetDirectory(m_dir_main);

  hClusterEcalib_E_phi = new TH3D("hClusterEcalib_E_phi", "; E_{MC} (GeV); E_{rec,rec hit}/E_{MC}; #varphi (rad)", 
                                  1500, 0., 150.0, 200, 0., 2.0, 360 , 0, 2*TMath::Pi());
  hClusterEcalib_E_phi->SetDirectory(m_dir_main);
  hClusterESimcalib_E_phi = new TH3D("hClusterESimcalib_E_phi", "; E_{MC} (GeV); E_{rec,sim hit}/E_{MC}; #varphi (rad)" , 1500, 0., 150.0, 200, 0., 2.0, 360 , 0, 2*TMath::Pi());
  hClusterESimcalib_E_phi->SetDirectory(m_dir_main);


  hCellESim_layerZ = new TH2D("hCellESim_layerZ", "; #cell ID Z; E_{rec,sim hit} (GeV)" , 70, -0.5, 69.5, 5000, 0, 1);
  hCellESim_layerZ->SetDirectory(m_dir_main);
  hCellESim_layerX = new TH2D("hCellESim_layerX", "; #cell ID X; E_{rec,sim hit} (GeV)" , 240, -0.5, 239.5, 5000, 0, 1);
  hCellESim_layerX->SetDirectory(m_dir_main);
  hCellESim_layerY = new TH2D("hCellESim_layerY", "; #cell ID Y; E_{rec,sim hit} (GeV)" , 240, -0.5, 239.5, 5000, 0, 1);
  hCellESim_layerY->SetDirectory(m_dir_main);
  hCellTSim_layerZ = new TH2D("hCellTSim_layerZ", "; #cell ID Z; t_{rec,sim hit} (GeV)" , 70, -0.5, 69.5, 5000, 0, 1000);
  hCellTSim_layerZ->SetDirectory(m_dir_main);
  
  hSamplingFractionEta = new TH2D("hSamplingFractionEta", "; #eta; f", 400, 1., 5., 500, 0., 0.2);
  hSamplingFractionEta->SetDirectory(m_dir_main);  
  
  hPosCaloModulesXY = new TH2D("hPosCaloModulesXY", "; module ID X; module ID Y", 54, 0., 54., 54, 0., 54.);
  hPosCaloModulesXY->SetDirectory(m_dir_main);

  hCaloCellIDs = new TH3D("hCaloCellIDs", "; id Z; id x; id y", 7, -0.5, 6.5, 54*2, 0., 54.*2, 54*2, 0., 54.*2);
  hCaloCellIDs->Sumw2();
  hCaloCellIDs->SetDirectory(m_dir_main);
  hCaloCellIDs_xy = new TH2D("hCaloCellIDs_xy", "; id x; id y", 54*2, 0., 54.*2, 54*2, 0., 54.*2);
  hCaloCellIDs_xy->SetDirectory(m_dir_main);
  hCaloCellIDsHCluster = new TH3D("hCaloCellIDsHighestCluster", "; id Z; id x; id y", 7, -0.5, 6.5, 54*2, 0., 54.*2, 54*2, 0., 54.*2);
  hCaloCellIDsHCluster->Sumw2();
  hCaloCellIDsHCluster->SetDirectory(m_dir_main);
  
  hPosCaloHitsXY = new TH2D("hPosCaloHitsXY", "; X (cm); Y (cm)", 400, -400., 400., 400, -400., 400.);
  hPosCaloHitsZX = new TH2D("hPosCaloHitsZX", "; Z (cm); X (cm)", 200, 300., 500., 400, -400., 400.);
  hPosCaloHitsZY = new TH2D("hPosCaloHitsZY", "; Z (cm); Y (cm)", 200, 300., 500., 400, -400., 400.);
  hPosCaloHitsXY->SetDirectory(m_dir_main);
  hPosCaloHitsZX->SetDirectory(m_dir_main);
  hPosCaloHitsZY->SetDirectory(m_dir_main);

  hPosCaloSimHitsXY = new TH2D("hPosCaloSimHitsXY", "; X (cm); Y (cm)", 400, -400., 400., 400, -400., 400.);
  hPosCaloSimHitsZX = new TH2D("hPosCaloSimHitsZX", "; Z (cm); X (cm)", 200, 300., 500., 400, -400., 400.);
  hPosCaloSimHitsZY = new TH2D("hPosCaloSimHitsZY", "; Z (cm); Y (cm)", 200, 300., 500., 400, -400., 400.);
  hPosCaloSimHitsXY->SetDirectory(m_dir_main);
  hPosCaloSimHitsZX->SetDirectory(m_dir_main);
  hPosCaloSimHitsZY->SetDirectory(m_dir_main);

  
  hCaloCellIDs_xy8M = new TH2D("hCaloCellIDs_xy8M", "; id x; id x", 54*2, 0., 54.*2, 54*2, 0., 54.*2);
  hCaloCellIDs_xy4M = new TH2D("hCaloCellIDs_xy4M", "; id x; id x", 54*2, 0., 54.*2, 54*2, 0., 54.*2);
  hPosCaloHitsXY4M = new TH2D("hPosCaloHitsXY4M", "; X (cm); Y (cm)", 400, -400., 400., 400, -400., 400.);
  hCaloCellIDs_xy8M->SetDirectory(m_dir_main);    
  hCaloCellIDs_xy4M->SetDirectory(m_dir_main);  
  hPosCaloHitsXY4M->SetDirectory(m_dir_main);
  
  lFHCal_towers_cellE = new float[maxNTowers];
  lFHCal_towers_cellT = new float[maxNTowers];
  lFHCal_towers_cellIDx = new short[maxNTowers];
  lFHCal_towers_cellIDy = new short[maxNTowers];
  lFHCal_towers_cellIDz = new short[maxNTowers];
  lFHCal_towers_clusterIDA = new short[maxNTowers];
  lFHCal_towers_clusterIDB = new short[maxNTowers];
  lFHCal_towers_cellTrueID = new int[maxNTowers];
  
  event_tree = new TTree("event_tree", "event_tree");
  event_tree->SetDirectory(m_dir_main);
  
    // towers LFHCALO
  event_tree->Branch("tower_LFHCAL_N", &lFHCal_towers_N, "tower_LFHCAL_N/I");
  event_tree->Branch("tower_LFHCAL_E", lFHCal_towers_cellE, "tower_LFHCAL_E[tower_LFHCAL_N]/F");
  event_tree->Branch("tower_LFHCAL_T", lFHCal_towers_cellT, "tower_LFHCAL_T[tower_LFHCAL_N]/F");
  event_tree->Branch("tower_LFHCAL_ix", lFHCal_towers_cellIDx, "tower_LFHCAL_ix[tower_LFHCAL_N]/S");
  event_tree->Branch("tower_LFHCAL_iy", lFHCal_towers_cellIDy, "tower_LFHCAL_iy[tower_LFHCAL_N]/S");
  event_tree->Branch("tower_LFHCAL_iz", lFHCal_towers_cellIDz, "tower_LFHCAL_iz[tower_LFHCAL_N]/S");
  event_tree->Branch("tower_LFHCAL_clusIDA", lFHCal_towers_clusterIDA, "tower_LFHCAL_clusIDA[tower_LFHCAL_N]/S");
  event_tree->Branch("tower_LFHCAL_clusIDB", lFHCal_towers_clusterIDB, "tower_LFHCAL_clusIDB[tower_LFHCAL_N]/S");
  event_tree->Branch("tower_LFHCAL_trueID", lFHCal_towers_cellTrueID, "tower_LFHCAL_trueID[tower_LFHCAL_N]/I");

  
  std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
  dd4hep::Detector& detector = dd4hep::Detector::getInstance();
  std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
  // detector.fromCompact("../epic/install/share/epic/epic_gfhcal_only.xml");
  std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
  dd4hep::rec::CellIDPositionConverter cellid_converter(detector);
  std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;

  std::cout << "--------------------------\nID specification:\n";
  try {
  
    m_decoder         = detector.readout("LFHCALHits").idSpec().decoder();
    std::cout <<"1st: "<< m_decoder << std::endl;
    auto module_index_x = m_decoder->index("moduleIDx");
    auto module_index_y = m_decoder->index("moduleIDy");
    iLx = m_decoder->index("towerx");
    iLy = m_decoder->index("towery");
    iLz = m_decoder->index("layerz");
    auto rlayer_index_z = m_decoder->index("rlayerz");
    iPassive = m_decoder->index("passive");
    std::cout << "full list: " << " " << m_decoder->fieldDescription() << std::endl;
  } catch (...) {
    isLFHCal = false;
    std::cout <<"1st: "<< m_decoder << std::endl;
    try {
      m_decoder         = detector.readout("GFHCALHits").idSpec().decoder();
      std::cout <<"2nd: "  << m_decoder << std::endl;
      auto module_index_x = m_decoder->index("modulex");
      auto module_index_y = m_decoder->index("moduley");
      iPassive = m_decoder->index("passive");
      iLx = m_decoder->index("layerx");
      iLy = m_decoder->index("layery");
      iLz = m_decoder->index("layerz");
      std::cout << "full list: " << " " << m_decoder->fieldDescription() << std::endl;
      
      nameSimHits         = "GFHCALHits";
      nameRecHits         = "GFHCALRecHits";
      nameClusters        = "GFHCALClusters";
      nameProtoClusters   = "GFHCALIslandProtoClusters";

    } catch (...) {
        std::cout <<"2nd: "  << m_decoder << std::endl;
        m_log->error("readoutClass not in the output");
        throw std::runtime_error("readoutClass not in the output.");
    }
  }
  
}

//-------------------------------------------
// ProcessSequential
//-------------------------------------------
void lfhcal_studiesProcessor::ProcessSequential(const std::shared_ptr<const JEvent>& event) {
  using namespace std;

  // cout << "lfhcal_studiesProcessor::ProcessSequential" << endl;
  double mceta = 0;
  double mcphi = 0;
  double mcp   = 0;
  double mcenergy = 0;
  for (auto mcparticle : mcParticles()) {
    if (mcparticle->getGeneratorStatus() != 1)
      continue;
    auto& mom = mcparticle->getMomentum();
    // get particle energy
    mcenergy = mcparticle->getEnergy();
    //determine mceta from momentum
    mceta = -log(tan(atan2(sqrt(mom.x * mom.x + mom.y * mom.y), mom.z) / 2.));
    // determine mcphi from momentum
    mcphi = atan2(mom.y, mom.x);
    // determine mc momentum
    mcp = sqrt(mom.x * mom.x + mom.y * mom.y + mom.z * mom.z);

//     std::cout << "MC particle: " << mom.x << " " << mom.y << " " << mom.z << "\ttotmom: " <<
//     mcp << "\t phi: "<< mcphi << "\t eta:" <<  mceta << std::endl; 
    
    hMCEnergyVsEta->Fill(mcp,mceta);
  }
  
  auto simHits = event -> Get<edm4hep::SimCalorimeterHit>(nameSimHits.data());


  std::vector<towersStrct> input_tower_sim;
  int nCaloHitsSim = 0;
  float sumActiveCaloEnergy = 0;
  float sumPassiveCaloEnergy = 0;
  
  
  // process sim hits
  for (auto caloHit : simHits) {
    float x         = caloHit->getPosition().x / 10.;
    float y         = caloHit->getPosition().y / 10.;
    float z         = caloHit->getPosition().z / 10.;
    uint64_t cellID = caloHit->getCellID();
    float energy    = caloHit->getEnergy();
    double time = std::numeric_limits<double>::max();
    for (const auto& c : caloHit->getContributions()) {
        if (c.getTime() <= time) {
            time = c.getTime();
        }
    }
    
    auto detector_module_x  = m_decoder->get(cellID, 1);
    auto detector_module_y  = m_decoder->get(cellID, 2);
    auto detector_passive = m_decoder->get(cellID, iPassive);
    auto detector_layer_x = m_decoder->get(cellID, iLx);
    auto detector_layer_y = m_decoder->get(cellID, iLy);
    int detector_layer_rz = -1;
    if (isLFHCal)
      detector_layer_rz = m_decoder->get(cellID, 7);
    auto detector_layer_z = m_decoder->get(cellID, iLz);
    if(detector_passive == 0) {
      sumActiveCaloEnergy += energy;
    } else {
      sumPassiveCaloEnergy += energy;
    }
    
    if (detector_passive > 0) continue;    
    // calc cell IDs
    int cellIDx = -1;
    int cellIDy = -1;
    int cellIDz = -1;
    if (isLFHCal){
      cellIDx = 54*2 - detector_module_x * 2 + detector_layer_x;
      cellIDy = 54*2 - detector_module_y * 2 + detector_layer_y;
      cellIDz = detector_layer_rz*10+detector_layer_z;
    } else {
      cellIDx = 54*2 - detector_module_x * 4 + detector_layer_x;
      cellIDy = 54*2 - detector_module_y * 2 + detector_layer_y;
      cellIDz = detector_layer_z;      
    }
    nCaloHitsSim++;
    
    hPosCaloSimHitsXY->Fill(x, y);
    hPosCaloSimHitsZX->Fill(z, x);
    hPosCaloSimHitsZY->Fill(z, y);

    hCellESim_layerZ->Fill(cellIDz, energy);
    hCellESim_layerX->Fill(cellIDx, energy);
    hCellESim_layerY->Fill(cellIDy, energy);
    hCellTSim_layerZ->Fill(cellIDz, time);
    
    //loop over input_tower_sim and find if there is already a tower with the same cellID
    bool found = false;
    for (auto& tower : input_tower_sim) {
      if (tower.cellID == cellID) {
        tower.energy += energy;
        found = true;
        break;
      }
    }
    if (!found) {
      towersStrct tempstructT;
      tempstructT.energy        = energy;
      tempstructT.time          = time; 
      tempstructT.posx          = x; 
      tempstructT.posy          = y; 
      tempstructT.posz          = z; 
      tempstructT.cellID        = cellID;
      tempstructT.cellIDx       = cellIDx;
      tempstructT.cellIDy       = cellIDy;
      tempstructT.cellIDz       = cellIDz;
      tempstructT.tower_trueID  = 0; //TODO how to get trueID?
      input_tower_sim.push_back(tempstructT);
    }
  }

  auto recHits = event -> Get<edm4eic::CalorimeterHit>(nameRecHits.data());  
  int nCaloHitsRec = 0;
  std::vector<towersStrct> input_tower_rec;
  std::vector<towersStrct> input_tower_recSav;
  // process rec hits
  for (auto caloHit : recHits) {
    float x         = caloHit->getPosition().x / 10.;
    float y         = caloHit->getPosition().y / 10.;
    float z         = caloHit->getPosition().z / 10.;
    uint64_t cellID = caloHit->getCellID();
    float energy    = caloHit->getEnergy();
    float time      = caloHit->getTime();

    auto detector_module_x  = m_decoder->get(cellID, 1);
    auto detector_module_y  = m_decoder->get(cellID, 2);
    auto detector_passive = m_decoder->get(cellID, iPassive);
    auto detector_layer_x = m_decoder->get(cellID, iLx);
    auto detector_layer_y = m_decoder->get(cellID, iLy);
    int detector_layer_rz = -1;
    int detector_module_t  = 0;
    if (isLFHCal){
      detector_layer_rz   = m_decoder->get(cellID, 7);
      detector_module_t   = m_decoder->get(cellID, 3);
    }
    auto detector_layer_z = m_decoder->get(cellID, iLz);

    if (detector_passive > 0) continue;
    
    // calc cell IDs
    int cellIDx = -1;
    int cellIDy = -1;
    int cellIDz = -1;
    if (isLFHCal){
      cellIDx = 54*2 - detector_module_x * 2 + detector_layer_x;
      cellIDy = 54*2 - detector_module_y * 2 + detector_layer_y;
      cellIDz = detector_layer_rz;
    } else {
      cellIDx = 54*2 - detector_module_x * 4 + detector_layer_x;
      cellIDy = 54*2 - detector_module_y * 2 + detector_layer_y;
      cellIDz = detector_layer_z;      
    }

    
    
    hCaloCellIDs->Fill(cellIDz,cellIDx, cellIDy, energy);
    hCaloCellIDs_xy->Fill(cellIDx, cellIDy);
     
    if (detector_module_t != 0){
      hCaloCellIDs_xy4M->Fill(cellIDx, cellIDy);
      hPosCaloHitsXY4M->Fill(x, y);
    } else {
      hCaloCellIDs_xy8M->Fill(cellIDx, cellIDy);
    }
         
    hPosCaloHitsXY->Fill(x, y);
    hPosCaloHitsZX->Fill(z, x);
    hPosCaloHitsZY->Fill(z, y);

    hPosCaloModulesXY->Fill(detector_module_x, detector_module_y);
    nCaloHitsRec++;
    
    //loop over input_tower_rec and find if there is already a tower with the same cellID
    bool found = false;
    for (auto& tower : input_tower_rec) {
      if (tower.cellID == cellID) {
        tower.energy += energy;
        found = true;
        break;
      }
    }
    if (!found) {
      towersStrct tempstructT;
      tempstructT.energy        = energy; 
      tempstructT.time          = time; 
      tempstructT.posx          = x; 
      tempstructT.posy          = y; 
      tempstructT.posz          = z; 
      tempstructT.cellID        = cellID;
      tempstructT.cellIDx       = cellIDx;
      tempstructT.cellIDy       = cellIDy;
      if (isLFHCal)
        tempstructT.cellIDz       = detector_layer_rz;
      else 
        tempstructT.cellIDz       = cellIDz;
      tempstructT.tower_trueID  = 0; //TODO how to get trueID?
      input_tower_rec.push_back(tempstructT);
      input_tower_recSav.push_back(tempstructT);
    }
  }
//   cout << "nCaloHits sim " << nCaloHitsSim << "\t rec " << nCaloHitsRec << endl;
  if (nCaloHitsRec > 0) nEventsWithCaloHits++;
  
  hSamplingFractionEta->Fill(mceta, sumActiveCaloEnergy / (sumActiveCaloEnergy+sumPassiveCaloEnergy));  
  std::sort(input_tower_rec.begin(), input_tower_rec.end(), &acompare); 
  std::sort(input_tower_recSav.begin(), input_tower_recSav.end(), &acompare); 
  std::sort(input_tower_sim.begin(), input_tower_sim.end(), &acompare); 
  
  // print towers rec hits
  double tot_energyRecHit = 0;
  for (auto& tower : input_tower_rec) {
    tower.energy = tower.energy; // calibrate
    tot_energyRecHit += tower.energy;
  }

  double samplingFractionFe = 0.037;
  double samplingFractionW  = 0.019;
  int minCellIDzDiffSamp    = 5;
  if (!isLFHCal){
    samplingFractionFe = 0.031;
    samplingFractionW  = 0.016;
    minCellIDzDiffSamp = 1; 
  }
  // print towers sim hits
  double tot_energySimHit = 0;
  for (auto& tower : input_tower_sim) {
    if (tower.cellIDz < minCellIDzDiffSamp)
      tower.energy = tower.energy/samplingFractionW; // calibrate
    else 
      tower.energy = tower.energy/samplingFractionFe; // calibrate
    tot_energySimHit += tower.energy;
  }
//   std::cout << "Mc E: " << mcenergy << "\t eta: " << mceta << "\t sim E rec: " << tot_energySimHit << "\t rec E rec: " <<  tot_energyRecHit << std::endl;
  
  hClusterNCells_E_eta->Fill(mcenergy, nCaloHitsRec, mceta); 
  hClusterSimNCells_E_eta->Fill(mcenergy, nCaloHitsSim, mceta); 
  
  hClusterEcalib_E_eta->Fill(mcenergy, tot_energyRecHit/mcenergy, mceta); 
  hClusterESimcalib_E_eta->Fill(mcenergy, tot_energySimHit/mcenergy, mceta);   
  hClusterEcalib_E_phi->Fill(mcenergy, tot_energyRecHit/mcenergy, mcphi); 
  hClusterESimcalib_E_phi->Fill(mcenergy, tot_energySimHit/mcenergy, mcphi);   
  
  
  // MA clusterization
  int removedCells  = 0;
  float minAggE     = 0.001;
  float seedE       = 0.100;
  if (!isLFHCal){
    minAggE     = 0.0005;
    seedE       = 0.20;
  }
  
  
  if (input_tower_rec.size()> 0){
    // clean up rec array for clusterization
    while (input_tower_rec.at(input_tower_rec.size()-1).energy < minAggE ){
      input_tower_rec.pop_back();
      removedCells++;
    }
//     std::cout << "removed " << removedCells << " with E < "  << minAggE << "GeV" << std::endl;
    
    
    int nclusters = 0;
    // vector of clusters
    std::vector<clustersStrct> clusters_calo;
    // vector of towers within the currently found cluster
    std::vector<towersStrct> cluster_towers;
    while (!input_tower_rec.empty() ) {
      cluster_towers.clear();
      clustersStrct tempstructC;
      // always start with highest energetic tower
      if(input_tower_rec.at(0).energy > seedE){
//         std::cout<< "seed: " << input_tower_rec.at(0).energy << "\t" << input_tower_rec.at(0).cellIDx <<  "\t" << input_tower_rec.at(0).cellIDy<<  "\t" << input_tower_rec.at(0).cellIDz<< std::endl;
        tempstructC = findMACluster(seedE, minAggE, input_tower_rec, cluster_towers);

        // determine remaining cluster properties from its towers
        float* showershape_eta_phi = CalculateM02andWeightedPosition(cluster_towers, tempstructC.cluster_E, 4.5);
        tempstructC.cluster_M02 = showershape_eta_phi[0];
        tempstructC.cluster_M20 = showershape_eta_phi[1];
        tempstructC.cluster_Eta = showershape_eta_phi[2];
        tempstructC.cluster_Phi = showershape_eta_phi[3];
        tempstructC.cluster_X = showershape_eta_phi[4];
        tempstructC.cluster_Y = showershape_eta_phi[5];
        tempstructC.cluster_Z = showershape_eta_phi[6];
        tempstructC.cluster_towers = cluster_towers;
//         std::cout <<  "---------> \t " << nclusters << "\tcluster with E = " << tempstructC.cluster_E << "\tEta: " << tempstructC.cluster_Eta<< "\tPhi: " << tempstructC.cluster_Phi
//                                 << "\tX: " << tempstructC.cluster_X<< "\tY: " << tempstructC.cluster_Y<< "\tZ: " << tempstructC.cluster_Z<< "\tntowers: " << tempstructC.cluster_NTowers 
//                                 << "\ttrueID: " << tempstructC.cluster_trueID << std::endl;      
        clusters_calo.push_back(tempstructC);
        
        nclusters++;
      } else {
//         std::cout<< "remaining: "<< (int)input_tower_rec.size() << " largest:" << input_tower_rec.at(0).energy << "\t" << input_tower_rec.at(0).cellIDx <<  "\t" << input_tower_rec.at(0).cellIDy<<  "\t" << input_tower_rec.at(0).cellIDz<< std::endl;
//         for (int ait = 0; ait < (int)input_tower_rec.size(); ait++){
//           std::cout<< input_tower_rec.at(ait).energy << "\t" << input_tower_rec.at(ait).cellIDx <<  "\t" << input_tower_rec.at(ait).cellIDy <<   "\t" << input_tower_rec.at(0).cellIDz << std::endl;
  //         h_clusterizer_nonagg_towers[caloEnum][clusterizerEnum]->Fill(input_tower_rec.size(),input_tower_rec.at(ait).tower_E);
//         }
        input_tower_rec.clear();      
      }
    }
      
    std::sort(clusters_calo.begin(), clusters_calo.end(), &acompareCl);    
//     std::cout << "-----> found " << clusters_calo.size() << " clusters" << std::endl; 
    hRecNClusters_E_eta->Fill(mcenergy, clusters_calo.size(), mceta);  
    int iCl = 0;
    for (auto& cluster : clusters_calo) {
      hRecClusterEcalib_E_eta->Fill(mcenergy, cluster.cluster_E/mcenergy, mceta);
      for (int iCell = 0; iCell < (int)cluster.cluster_towers.size(); iCell++){
        int pSav = 0;
        while(cluster.cluster_towers.at(iCell).cellID !=  input_tower_recSav.at(pSav).cellID && pSav < (int)input_tower_recSav.size() ) pSav++;
        if (cluster.cluster_towers.at(iCell).cellID == input_tower_recSav.at(pSav).cellID)
          input_tower_recSav.at(pSav).tower_clusterIDA = iCl;
      }
      
      if (iCl == 0){
        hRecClusterEcalib_Ehigh_eta->Fill(mcenergy, cluster.cluster_E/mcenergy, mceta);  
        hRecClusterNCells_Ehigh_eta->Fill(mcenergy, cluster.cluster_NTowers, mceta);  
      }
      iCl++;
//       std::cout <<"MA cluster " << iCl << ":\t" << cluster.cluster_E << "\t"<< cluster.cluster_NTowers <<std::endl;
    }
    
    clusters_calo.clear();
  } else {
    hRecNClusters_E_eta->Fill(mcenergy, 0., mceta);  
  }
    
  int iClF          = 0;
  float highestEFr  = 0;
  int iClFHigh      = 0;

  auto lfhcalClustersF = event -> Get<edm4eic::Cluster>(nameClusters.data());  
  for (auto& cluster : lfhcalClustersF) {
    if (cluster->getEnergy() > highestEFr){
      iClFHigh    = iClF;
      highestEFr  = cluster->getEnergy();
    }
    hRecFClusterEcalib_E_eta->Fill(mcenergy, cluster->getEnergy()/mcenergy, mceta);        
    
//     std::cout << "Island cluster "<< iClF << ":\t" << cluster->getEnergy()  << "\t"<< cluster->getNhits()  << std::endl;
//     for (auto& protocluster : lfhcalProtoClustersF()) {
//       if (! )
    for (auto& hit: cluster->getHits()){
//       for (int iCell = 0;  iCell < (int)cluster->getHits().size(); iCell++){
        int pSav = 0;
        while(hit.getCellID() !=  input_tower_recSav.at(pSav).cellID && pSav < (int)input_tower_recSav.size() ) pSav++;
        if (hit.getCellID() == input_tower_recSav.at(pSav).cellID)
          input_tower_recSav.at(pSav).tower_clusterIDB = iClF;
      }
//     }
    iClF++;
  }
  hRecFNClusters_E_eta->Fill(mcenergy, iClF, mceta);  

  iClF          = 0;
  for (auto& cluster : lfhcalClustersF) {
    if (iClF == iClFHigh){
      hRecFClusterEcalib_Ehigh_eta->Fill(mcenergy, cluster->getEnergy()/mcenergy, mceta);  
      hRecFClusterNCells_Ehigh_eta->Fill(mcenergy, cluster->getNhits(), mceta);  
    }
    iClF++;
  }
  
  lFHCal_towers_N = (int)input_tower_recSav.size();
  for (int iCell = 0; iCell < (int)input_tower_recSav.size(); iCell++){
//     std::cout << input_tower_recSav.at(iCell).cellIDx << "\t" << input_tower_recSav.at(iCell).cellIDy << "\t" << input_tower_recSav.at(iCell).cellIDz << "\t" << input_tower_recSav.at(iCell).energy << "\t" << input_tower_recSav.at(iCell).tower_clusterIDA << "\t" << input_tower_recSav.at(iCell).tower_clusterIDB << std::endl;
    
    lFHCal_towers_cellE[iCell]      = (float)input_tower_recSav.at(iCell).energy;
    lFHCal_towers_cellT[iCell]      = (float)input_tower_recSav.at(iCell).time;
    lFHCal_towers_cellIDx[iCell]    = (short)input_tower_recSav.at(iCell).cellIDx;
    lFHCal_towers_cellIDy[iCell]    = (short)input_tower_recSav.at(iCell).cellIDy;
    lFHCal_towers_cellIDz[iCell]    = (short)input_tower_recSav.at(iCell).cellIDz;
    lFHCal_towers_clusterIDA[iCell] = (short)input_tower_recSav.at(iCell).tower_clusterIDA;
    lFHCal_towers_clusterIDB[iCell] = (short)input_tower_recSav.at(iCell).tower_clusterIDB;
    lFHCal_towers_cellTrueID[iCell] = (int)input_tower_recSav.at(iCell).tower_trueID;
  }
  
  event_tree->Fill();
  
  lFHCal_towers_N = 0;
  for (Int_t itow = 0; itow < maxNTowers; itow++){
    lFHCal_towers_cellE[itow]       = 0;
    lFHCal_towers_cellT[itow]       = 0;
    lFHCal_towers_cellIDx[itow]     = 0;
    lFHCal_towers_cellIDy[itow]     = 0;
    lFHCal_towers_cellIDz[itow]     = 0;
    lFHCal_towers_clusterIDA[itow]  = 0;
    lFHCal_towers_clusterIDB[itow]  = 0;
    lFHCal_towers_cellTrueID[itow]  = 0;
  }
  
}

//-------------------------------------------
// FinishWithGlobalRootLock
//-------------------------------------------
void lfhcal_studiesProcessor::FinishWithGlobalRootLock() {
  std::cout << "------> " << nEventsWithCaloHits << " with calo info present"<< std::endl;
  // Do any final calculations here.
    
}
