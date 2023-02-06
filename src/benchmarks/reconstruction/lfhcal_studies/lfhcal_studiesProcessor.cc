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

#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Plugins/DD4hep/ConvertDD4hepDetector.hpp"

// #include <extensions/spdlog/SpdlogMixin.h>

// The following just makes this a JANA plugin
extern "C" {
void InitPlugin(JApplication* app) {
  InitJANAPlugin(app);
  app->Add(new lfhcal_studiesProcessor());
}
}

struct towersStrct{
  towersStrct(): energy(0), cellIDx(-1), cellIDy(-1), cellIDz(-1), tower_trueID(-10000) {}
  float energy;
  int cellIDx;
  int cellIDy;
  int cellIDz;
  int tower_trueID;
} ;
bool acompare(towersStrct lhs, towersStrct rhs) { return lhs.energy > rhs.energy; }

//-------------------------------------------
// InitWithGlobalRootLock
//-------------------------------------------
void lfhcal_studiesProcessor::InitWithGlobalRootLock() {
  std::string plugin_name = ("lfhcal_studies");

  // InitLogger(plugin_name);
  // Get JANA application
  auto app          = GetApplication();
  auto acts_service = GetApplication()->GetService<ACTSGeo_service>();

  std::string log_level_str = "debug";
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
  hClusterEcalib_E_eta = new TH3D("hClusterEcalib_E_eta", "; E_{MC} (GeV); E_{rec,rec hit}/E_{MC}; #eta", 1500, 0., 150.0, 200, 0., 2.0, 20, 1, 5);
  hClusterEcalib_E_eta->SetDirectory(m_dir_main);
  hClusterESimcalib_E_eta = new TH3D("hClusterESimcalib_E_eta", "; E_{MC} (GeV); E_{rec,sim hit}/E_{MC}; #eta" , 1500, 0., 150.0, 200, 0., 2.0, 20, 1, 5);
  hClusterESimcalib_E_eta->SetDirectory(m_dir_main);
  hClusterEcalib_E_phi = new TH3D("hClusterEcalib_E_phi", "; E_{MC} (GeV); E_{rec,rec hit}/E_{MC}; #varphi (rad)", 1500, 0., 150.0, 200, 0., 2.0, 360 , 0, 2*TMath::Pi());
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

  hCaloCellIDs = new TH3D("hCaloCellIDs", "; id Z; id x; id x", 70, 0, 70, 54*2, 0., 54.*2, 54*2, 0., 54.*2);
  hCaloCellIDs->SetDirectory(m_dir_main);
  hCaloCellIDs_xy = new TH2D("hCaloCellIDs_xy", "; id x; id x", 54*2, 0., 54.*2, 54*2, 0., 54.*2);
  hCaloCellIDs_xy->SetDirectory(m_dir_main);
  
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
  
  
  std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
  dd4hep::Detector& detector = dd4hep::Detector::getInstance();
  std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
  // detector.fromCompact("../epic/install/share/epic/epic_gfhcal_only.xml");
  std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
  dd4hep::rec::CellIDPositionConverter cellid_converter(detector);
  std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;

  std::cout << "--------------------------\nID specification:\n";
  m_decoder         = detector.readout("LFHCALHits").idSpec().decoder();
  auto module_index_x = m_decoder->index("moduleIDx");
  std::cout << "modulex index is " << module_index_x << std::endl;
  auto module_index_y = m_decoder->index("moduleIDy");
  std::cout << "moduley index is " << module_index_y << std::endl;
  auto layer_index_x = m_decoder->index("towerx");
  std::cout << "layerx index is " << layer_index_x << std::endl;
  auto layer_index_y = m_decoder->index("towery");
  std::cout << "layery index is " << layer_index_y << std::endl;
  auto layer_index_z = m_decoder->index("layerz");
  std::cout << "layerz index is " << layer_index_z << std::endl;
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
  std::vector<towersStrct> input_tower_sim;
  int nCaloHitsSim = 0;
  float sumActiveCaloEnergy = 0;
  float sumPassiveCaloEnergy = 0;
  
  // process sim hits
  for (auto caloHit : lfhcalSimHits()) {
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
    auto detector_layer_x = m_decoder->get(cellID, 4);
    auto detector_layer_y = m_decoder->get(cellID, 5);
    auto detector_layer_z = m_decoder->get(cellID, 6);
    auto detector_passive = m_decoder->get(cellID, 7);
    if(detector_passive == 0) {
      sumActiveCaloEnergy += energy;
    } else {
      sumPassiveCaloEnergy += energy;
    }
    
    if (detector_passive == 1) continue;    
    // calc cell IDs
    int cellIDx = 54*2 - detector_module_x * 2 + detector_layer_x;
    int cellIDy = 54*2 - detector_module_y * 2 + detector_layer_y;
    int cellIDz = detector_layer_z;
    nCaloHitsSim++;
    
    hPosCaloSimHitsXY->Fill(x, y);
    hPosCaloSimHitsZX->Fill(z, x);
    hPosCaloSimHitsZY->Fill(z, y);

    hCellESim_layerZ->Fill(detector_layer_z, energy);
    hCellESim_layerX->Fill(cellIDx, energy);
    hCellESim_layerY->Fill(cellIDy, energy);
    hCellTSim_layerZ->Fill(detector_layer_z, time);
    
    //loop over input_tower_sim and find if there is already a tower with the same cellID
    bool found = false;
    for (auto& tower : input_tower_sim) {
      if ((tower.cellIDx == cellIDx) && (tower.cellIDy == cellIDy) && (tower.cellIDz == cellIDz)) {
        tower.energy += energy;
        found = true;
        break;
      }
    }
    if (!found) {
      towersStrct tempstructT;
      tempstructT.energy       = energy; 
      tempstructT.cellIDx    = cellIDx;
      tempstructT.cellIDy    = cellIDy;
      tempstructT.cellIDz      = cellIDz;
      tempstructT.tower_trueID  = 0; //TODO how to get trueID?
      input_tower_sim.push_back(tempstructT);
    }
    
  }

  int nCaloHitsRec = 0;
  std::vector<towersStrct> input_tower_rec;
  // process rec hits
  for (auto caloHit : lfhcalRecHits()) {
    float x         = caloHit->getPosition().x / 10.;
    float y         = caloHit->getPosition().y / 10.;
    float z         = caloHit->getPosition().z / 10.;
    uint64_t cellID = caloHit->getCellID();
    float energy    = caloHit->getEnergy();
//     cout << "Calo hit: " << x << " " << y << " " << z << "\tcellID: " << cellID
//          << "\tenergy: " << energy << endl;

    auto detector_module_x  = m_decoder->get(cellID, 1);
    auto detector_module_y  = m_decoder->get(cellID, 2);
    auto detector_module_t  = m_decoder->get(cellID, 3);
    auto detector_layer_x = m_decoder->get(cellID, 4);
    auto detector_layer_y = m_decoder->get(cellID, 5);
    auto detector_layer_z = m_decoder->get(cellID, 6);
    auto detector_passive = m_decoder->get(cellID, 7);
    if (detector_passive == 1) continue;
    
    // calc cell IDs
    int cellIDx = 54*2 - detector_module_x * 2 + detector_layer_x;
    int cellIDy = 54*2 - detector_module_y * 2 + detector_layer_y;
    int cellIDz = detector_layer_z;
    hCaloCellIDs->Fill(cellIDz,cellIDx, cellIDy);
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
      if ((tower.cellIDx == cellIDx) && (tower.cellIDy == cellIDy) && (tower.cellIDz == cellIDz)) {
        tower.energy += energy;
        found = true;
        break;
      }
    }
    if (!found) {
      towersStrct tempstructT;
      tempstructT.energy       = energy; 
      tempstructT.cellIDx    = cellIDx;
      tempstructT.cellIDy    = cellIDy;
      tempstructT.cellIDz      = cellIDz;
      tempstructT.tower_trueID  = 0; //TODO how to get trueID?
      input_tower_rec.push_back(tempstructT);
    }
    
  }
//   cout << "nCaloHits sim " << nCaloHitsSim << "\t rec " << nCaloHitsRec << endl;
  if (nCaloHitsRec > 0) nEventsWithCaloHits++;
  
  hSamplingFractionEta->Fill(mceta, sumActiveCaloEnergy / (sumActiveCaloEnergy+sumPassiveCaloEnergy));  
  std::sort(input_tower_rec.begin(), input_tower_rec.end(), &acompare); 
  std::sort(input_tower_sim.begin(), input_tower_sim.end(), &acompare); 
  
  // print towers rec hits
  double tot_energyRecHit = 0;
  for (auto& tower : input_tower_rec) {
    tower.energy = tower.energy; // calibrate
    tot_energyRecHit += tower.energy;
  }

  double samplingFractionFe = 0.037;
  double samplingFractionW  = 0.019;
  // print towers sim hits
  double tot_energySimHit = 0;
  for (auto& tower : input_tower_sim) {
    if (tower.cellIDz < 5)
      tower.energy = tower.energy/samplingFractionW; // calibrate
    else 
      tower.energy = tower.energy/samplingFractionFe; // calibrate
    tot_energySimHit += tower.energy;
  }
//   std::cout << "Mc E: " << mcenergy << "\t eta: " << mceta << "\t sim E rec: " << tot_energySimHit << "\t rec E rec: " <<  tot_energyRecHit << std::endl;
  
  hClusterEcalib_E_eta->Fill(mcenergy, tot_energyRecHit/mcenergy, mceta); 
  hClusterESimcalib_E_eta->Fill(mcenergy, tot_energySimHit/mcenergy, mceta);   
  hClusterEcalib_E_phi->Fill(mcenergy, tot_energyRecHit/mcenergy, mcphi); 
  hClusterESimcalib_E_phi->Fill(mcenergy, tot_energySimHit/mcenergy, mcphi);   
}

//-------------------------------------------
// FinishWithGlobalRootLock
//-------------------------------------------
void lfhcal_studiesProcessor::FinishWithGlobalRootLock() {
  std::cout << "------> " << nEventsWithCaloHits << " with calo info present"<< std::endl;
  // Do any final calculations here.
}
