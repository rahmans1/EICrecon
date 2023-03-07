
//
// Template for this file generated with eicmkplugin.py
//

#include <JANA/JEventProcessorSequentialRoot.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TFile.h>
#include <TTree.h>

// Include appropirate class headers. e.g.
#include <edm4hep/SimCalorimeterHit.h>
// #include <detectors/BEMC/BEMCRawCalorimeterHit.h>

#include "Acts/Utilities/Logger.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include <Acts/Surfaces/DiscSurface.hpp>
#include "Acts/Definitions/Common.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Plugins/DD4hep/DD4hepDetectorElement.hpp"
#include <Acts/Material/IMaterialDecorator.hpp>
// #include <extensions/spdlog/SpdlogMixin.h>
#include <edm4hep//MCParticle.h>
#include <edm4hep/SimTrackerHit.h>
#include <edm4eic/RawTrackerHit.h>
#include <edm4eic/TrackerHit.h>
#include <edm4eic/CalorimeterHit.h>
#include <edm4eic/ReconstructedParticle.h>
#include <edm4eic/Cluster.h>
#include <edm4eic/ProtoCluster.h>

#include <services/log/Log_service.h>
#include <extensions/spdlog/SpdlogMixin.h>
#include <algorithms/tracking/ActsGeometryProvider.h>
#include <services/geometry/acts/ACTSGeo_service.h>
#include <algorithms/tracking/TrackPropagation.h>
#include <spdlog/spdlog.h>

class lfhcal_studiesProcessor: public JEventProcessorSequentialRoot {
private:

    // Data objects we will need from JANA e.g.
    PrefetchT<edm4hep::MCParticle>  mcParticles   = {this, "MCParticles" };
//     PrefetchT<edm4eic::CalorimeterHit> lfhcalRecHits = {this, "LFHCALRecHits"};
//     PrefetchT<edm4hep::SimCalorimeterHit> lfhcalSimHits = {this, "LFHCALHits"};
//     PrefetchT<edm4eic::Cluster> lfhcalClustersF = {this, "LFHCALClusters"};
//     PrefetchT<edm4eic::ProtoCluster> lfhcalProtoClustersF = {this, "LFHCALIslandProtoClusters"};    
    // Declare histogram and tree pointers here. e.g.

public:
    lfhcal_studiesProcessor() { SetTypeName(NAME_OF_THIS); }
    
    void InitWithGlobalRootLock() override;
    void ProcessSequential(const std::shared_ptr<const JEvent>& event) override;
    void FinishWithGlobalRootLock() override;
    TDirectory *m_dir_main;
    // Declare histogram and tree pointers here. e.g.
    TH2D* hMCEnergyVsEta;
    TH3D* hClusterEcalib_E_eta;
    TH3D* hClusterESimcalib_E_eta;
    TH3D* hRecClusterEcalib_E_eta;
    TH3D* hRecClusterEcalib_Ehigh_eta;
    TH3D* hRecFClusterEcalib_E_eta;
    TH3D* hRecFClusterEcalib_Ehigh_eta;
    TH3D* hClusterNCells_E_eta;
    TH3D* hClusterSimNCells_E_eta;
    TH3D* hRecClusterNCells_Ehigh_eta;
    TH3D* hRecNClusters_E_eta;
    TH3D* hRecFClusterNCells_Ehigh_eta;
    TH3D* hRecFNClusters_E_eta;
    TH3D* hClusterEcalib_E_phi;
    TH3D* hClusterESimcalib_E_phi;
    TH2D* hCellESim_layerZ;
    TH2D* hCellESim_layerX;
    TH2D* hCellESim_layerY;
    TH2D* hCellTSim_layerZ;
    TH2D* hSamplingFractionEta;
//     TH2D* hPosCaloModulesXY;
    TH2D* hPosCaloHitsXY;
    TH2D* hPosCaloHitsZX;
    TH2D* hPosCaloHitsZY;
    TH2D* hPosCaloSimHitsXY;
    TH2D* hPosCaloSimHitsZX;
    TH2D* hPosCaloSimHitsZY;
//     TH2D* hCaloCellIDs_xy;
//     TH3D* hCaloCellIDs;
//     TH2D* hCaloCellIDs_xy8M;
//     TH2D* hCaloCellIDs_xy4M;
//     TH2D* hPosCaloHitsXY4M;

    TTree* event_tree;
    const int maxNTowers = 65000;
    int lFHCal_towers_N;
    bool enableTree = false;
    short* lFHCal_towers_cellIDx;
    short* lFHCal_towers_cellIDy;
    short* lFHCal_towers_cellIDz;
    short* lFHCal_towers_clusterIDA;
    short* lFHCal_towers_clusterIDB;
    float* lFHCal_towers_cellE;
    float* lFHCal_towers_cellT;
    int* lFHCal_towers_cellTrueID;
    
    int nEventsWithCaloHits = 0;
    bool isLFHCal = true;
    std::shared_ptr<spdlog::logger> m_log;
    // Acts::GeometryContext geoContext;

//     /// Tracking propagation algorithm
//     eicrecon::TrackPropagation m_propagation_algo;

//     std::shared_ptr<const ActsGeometryProvider> m_geo_provider;
//     std::shared_ptr<Acts::DiscSurface> m_mRICH_center_surface;
//     std::shared_ptr<Acts::DiscSurface> m_dRICH_center_surface;
    dd4hep::BitFieldCoder* m_decoder;
    std::string nameSimHits         = "LFHCALHits";
    std::string nameRecHits         = "LFHCALRecHits";
    std::string nameClusters        = "LFHCALClusters";
    std::string nameProtoClusters   = "LFHCALIslandProtoClusters";
    short iPassive;
    short iLx;
    short iLy;
    short iLz;
    
};
