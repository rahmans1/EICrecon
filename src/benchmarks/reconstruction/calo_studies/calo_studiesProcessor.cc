#include "calo_studiesProcessor.h"
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
  app->Add(new calo_studiesProcessor());
}
}

//-------------------------------------------
// InitWithGlobalRootLock
//-------------------------------------------
void calo_studiesProcessor::InitWithGlobalRootLock() {
  std::string plugin_name = ("calo_studies");

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

  // Create histograms here. e.g.
  hPosTestX = new TH1D("hPosTestX", "BEMC hit energy (raw)", 100, -100, 100);
  hPosTestX->SetDirectory(m_dir_main);

  nHitsTrackVsEtaVsP =
      new TH3D("nHitsTrackVsEtaVsP", "", 100, -4., 4., 20, -0.5, 19.5, 150, 0., 15.);
  nHitsEventVsEtaVsP =
      new TH3D("nHitsEventVsEtaVsP", "", 100, -4., 4., 20, -0.5, 19.5, 150, 0., 15.);
  nHitsTrackVsEtaVsP->SetDirectory(m_dir_main);
  nHitsEventVsEtaVsP->SetDirectory(m_dir_main);

  hPosCaloModulesXY = new TH2D("hPosCaloModulesXY", "", 54, 0., 54., 54, 0., 54.);
  hPosCaloModulesXY->SetDirectory(m_dir_main);

  hPosCaloHitsXY = new TH2D("hPosCaloHitsXY", "", 400, -400., 400., 400, -400., 400.);
  hPosCaloHitsZX = new TH2D("hPosCaloHitsZX", "", 200, 300., 500., 400, -400., 400.);
  hPosCaloHitsZY = new TH2D("hPosCaloHitsZY", "", 200, 300., 500., 400, -400., 400.);
  hPosCaloHitsXY->SetDirectory(m_dir_main);
  hPosCaloHitsZX->SetDirectory(m_dir_main);
  hPosCaloHitsZY->SetDirectory(m_dir_main);

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
void calo_studiesProcessor::ProcessSequential(const std::shared_ptr<const JEvent>& event) {
  using namespace std;
  nEventsFilled++;
  // cout << "calo_studiesProcessor::ProcessSequential" << endl;
  double mceta = 0;
  double mcphi = 0;
  double mcp   = 0;
  for (auto mcparticle : mcParticles()) {
    if (mcparticle->getGeneratorStatus() != 1)
      continue;
    auto& mom = mcparticle->getMomentum();
    // cout << "MC particle: " << mom.x << " " << mom.y << " " << mom.z << "\ttotmom: " <<
    // sqrt(mom.x*mom.x+mom.y*mom.y+mom.z*mom.z) << endl; determine mceta from momentum
    mceta = -log(tan(atan2(sqrt(mom.x * mom.x + mom.y * mom.y), mom.z) / 2.));
    // determine mcphi from momentum
    mcphi = atan2(mom.y, mom.x);
    // determine mc momentum
    mcp = sqrt(mom.x * mom.x + mom.y * mom.y + mom.z * mom.z);
  }
  int nCaloHits = 0;
  for (auto caloHit : lfhcalHits()) {
    float x         = caloHit->getPosition().x / 10.;
    float y         = caloHit->getPosition().y / 10.;
    float z         = caloHit->getPosition().z / 10.;
    uint64_t cellID = caloHit->getCellID();
    float energy    = caloHit->getEnergy();
    cout << "Calo hit: " << x << " " << y << " " << z << "\tcellID: " << cellID
         << "\tenergy: " << energy << endl;

    auto detector_module_x  = m_decoder->get(cellID, 1);
    auto detector_module_y  = m_decoder->get(cellID, 2);
    auto detector_layer_x = m_decoder->get(cellID, 3);
    auto detector_layer_y = m_decoder->get(cellID, 4);
    auto detector_layer_z = m_decoder->get(cellID, 5);
    cout << "\tmodulex: " << detector_module_x << "\tmoduley: " << detector_module_y << "\tlayerx: " << detector_layer_x
         << "\tlayery: " << detector_layer_y << "\tlayerz: " << detector_layer_z << endl;
    // if ((detector_layer !=4 ) && (detector_layer !=5 )) {
    //   continue;
    // }
    hPosCaloHitsXY->Fill(x, y);
    hPosCaloHitsZX->Fill(z, x);
    hPosCaloHitsZY->Fill(z, y);

    hPosCaloModulesXY->Fill(detector_module_x, detector_module_y);
    nCaloHits++;
  }
  cout << "nCaloHits: " << nCaloHits << endl;
}

//-------------------------------------------
// FinishWithGlobalRootLock
//-------------------------------------------
void calo_studiesProcessor::FinishWithGlobalRootLock() {

  // Do any final calculations here.
  // nHitsTrackVsEtaVsP->Scale(1.0/nTracksFilled);
  // nHitsEventVsEtaVsP->Scale(1.0/nEventsFilled);
}
