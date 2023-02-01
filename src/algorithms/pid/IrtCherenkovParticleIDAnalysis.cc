// Copyright 2022, Christopher Dilks
// Subject to the terms in the LICENSE file found in the top-level directory.

#include "IrtCherenkovParticleIDAnalysis.h"


// RadiatorAnalysis constructor: defines histograms for a radiator
//---------------------------------------------------------------------------
eicrecon::RadiatorAnalysis::RadiatorAnalysis(TString rad_name) : m_rad_name(rad_name) {

  // distributions
  m_npe_dist = new TH1D(
      "npe_dist_"+m_rad_name,
      "Overall NPE for "+m_rad_name+";NPE",
      npe_bins, 0, npe_max
      );
  m_theta_dist = new TH1D(
      "theta_dist_"+m_rad_name,
      "Estimated Cherenkov Angle for "+m_rad_name+";#theta [mrad]",
      theta_bins, 0, theta_max
      );
  m_photon_theta_vs_phi = new TH2D(
      "photon_theta_vs_phi_"+m_rad_name,
      "Estimated Photon #theta vs #phi for "+m_rad_name+";#phi [rad];#theta [mrad]",
      phi_bins, -TMath::Pi(), TMath::Pi(),
      theta_bins, 0, theta_max
      );

  // truth
  m_mc_wavelength = new TH1D(
      "mc_wavelength_"+m_rad_name,
      "MC Photon Wavelength for "+m_rad_name+";#lambda [nm]",
      n_bins, 0, 1000
      );
  m_mc_rindex = new TH1D(
      "mc_rindex_"+m_rad_name,
      "MC Refractive Index for "+m_rad_name+";n",
      10*n_bins, 0.99, 1.03
      );

  // PID
  m_highest_weight_dist = new TH1D(
      "highest_weight_dist_"+m_rad_name,
      "Highest PDG Weight for "+m_rad_name+";PDG",
      pdg_bins(), 0, pdg_bins()
      );

  // momentum scans
  m_npe_vs_p = new TH2D(
      "npe_vs_p_"+m_rad_name,
      "Overall NPE vs. Particle Momentum for "+m_rad_name+";p [GeV];NPE",
      momentum_bins, 0, momentum_max,
      npe_bins, 0, npe_max
      );
  m_theta_vs_p = new TH2D(
      "theta_vs_p_"+m_rad_name,
      "Estimated Cherenkov Angle vs. Particle Momentum for "+m_rad_name+";p [GeV];#theta [mrad]",
      momentum_bins, 0, momentum_max,
      theta_bins, 0, theta_max
      );
  m_highest_weight_vs_p = new TH2D(
      "highest_weight_vs_p_"+m_rad_name,
      "Highest PDG Weight vs. Particle Momentum for "+m_rad_name+";p [GeV]",
      momentum_bins, 0, momentum_max,
      pdg_bins(), 0, pdg_bins()
      );
}


// AlgorithmInit
//---------------------------------------------------------------------------
void eicrecon::IrtCherenkovParticleIDAnalysis::AlgorithmInit(std::shared_ptr<spdlog::logger>& logger) {
  m_log = logger;

  // initialize histograms for each radiator
  for(auto& [id,rad_name] : Tools::GetRadiatorIDs())
    m_radiator_histos.insert({id,std::make_shared<RadiatorAnalysis>(TString(rad_name))});

  // initialize common histograms
  m_nphot_vs_p = new TH2D("nphot_vs_p", "N_{photons} vs. Thrown Momentum;p [GeV];N_{photons}",
      RadiatorAnalysis::momentum_bins, 0, RadiatorAnalysis::momentum_max,
      RadiatorAnalysis::nphot_max,     0, RadiatorAnalysis::nphot_max
      );
  m_nphot_vs_p__transient = new TH1D("nphot_vs_p__transient", "",
      m_nphot_vs_p->GetNbinsX(),
      m_nphot_vs_p->GetXaxis()->GetXmin(),
      m_nphot_vs_p->GetXaxis()->GetXmax()
      );
}


// AlgorithmProcess
//---------------------------------------------------------------------------
void eicrecon::IrtCherenkovParticleIDAnalysis::AlgorithmProcess(
    std::vector<const edm4hep::SimTrackerHit*>       sim_hits,
    std::vector<const edm4eic::CherenkovParticleID*> cherenkov_pids
    )
{
  m_log->trace("{:=^70}"," call IrtCherenkovParticleIDAnalysis::AlgorithmProcess ");

  // get the number of photons vs. momentum for this event
  // - get the true charged particle momentum, and fill 1D histogram `m_nphot_vs_p__transient`
  //   for each (pre-digitized) sensor hit
  // - number of entries in each momentum bin will be the number of photons for this event
  for(const auto& hit : sim_hits) {
    float momentum = 0;
    auto photon = hit->getMCParticle();
    if(photon.parents_size()>0) {
      auto charged_particle = photon.getParents(0);
      momentum = edm4hep::utils::p(charged_particle);
    }
    m_nphot_vs_p__transient->Fill(momentum);
  }
  // - use `m_nphot_vs_p__transient` results to fill 2D hist `m_nphot_vs_p` for this event
  for(int b=1; b<=m_nphot_vs_p__transient->GetNbinsX(); b++) {
    auto nphot = m_nphot_vs_p__transient->GetBinContent(b);
    if(nphot>0) {
      auto momentum = m_nphot_vs_p__transient->GetBinCenter(b);
      m_nphot_vs_p->Fill(momentum,nphot);
    }
  }
  // - clear `m_nphot_vs_p__transient` to be ready for the next event
  m_nphot_vs_p__transient->Reset();


  // loop over `CherenkovParticleID` objects
  for(const auto& pid : cherenkov_pids) {

    // get the histograms for this radiator
    std::shared_ptr<RadiatorAnalysis> radiator_histos;
    TString rad_name;
    try {
      radiator_histos = m_radiator_histos.at(pid->getRadiator());
      rad_name = radiator_histos->GetRadiatorName();
    }
    catch(const std::out_of_range& e) {
      m_log->error("Invalid radiator number {}", pid->getRadiator());
      continue;
    }
    m_log->trace("-> {} Radiator (ID={}):", rad_name, pid->getRadiator());

    // estimate the charged particle energy using the momentum of the first TrackPoint at this radiator's entrance
    auto charged_particle = pid->getChargedParticle();
    if(!charged_particle.isAvailable())   { m_log->warn("Charged particle not available in this radiator");      continue; }
    if(charged_particle.points_size()==0) { m_log->warn("Charged particle has no TrackPoints in this radiator"); continue; }
    auto charged_particle_momentum = edm4hep::utils::magnitude( charged_particle.getPoints(0).momentum );
    m_log->trace("  Charged Particle p = {} GeV at radiator entrance", charged_particle_momentum);
    m_log->trace("  If it is a pion, E = {} GeV", std::hypot(charged_particle_momentum, Tools::GetPDGMass(211)));

    // trace logging for IRT results
    m_log->trace("  Cherenkov Angle Estimate:");
    m_log->trace("    {:>16}:  {:<10}",     "NPE",          pid->getNpe());
    m_log->trace("    {:>16}:  {:<10.3}",   "<theta>",      pid->getTheta());
    m_log->trace("    {:>16}:  {:<10.3}",   "<rindex>",     pid->getRindex());
    m_log->trace("    {:>16}:  {:<10.3}",   "<wavelength>", pid->getWavelength());
    m_log->trace("  Mass Hypotheses:");
    m_log->trace("    {:>6}  {:>10}  {:>10}", "PDG", "Weight", "NPE");
    for(const auto& hyp : pid->getHypotheses())
      m_log->trace("    {:>6}  {:>10.8}  {:>10.8}", hyp.pdg, hyp.weight, hyp.npe);

    // Cherenkov angle estimate and NPE
    radiator_histos->m_npe_dist->Fill(pid->getNpe());
    radiator_histos->m_npe_vs_p->Fill(charged_particle_momentum,pid->getNpe());
    radiator_histos->m_theta_dist->Fill(pid->getTheta()*1e3); // [rad] -> [mrad]
    radiator_histos->m_theta_vs_p->Fill(charged_particle_momentum,pid->getTheta()*1e3); // [rad] -> [mrad]
    for(const auto& [theta,phi] : pid->getThetaPhiPhotons())
      radiator_histos->m_photon_theta_vs_phi->Fill(phi,theta*1e3); // [rad] -> [mrad]

    // fill MC dists
    radiator_histos->m_mc_wavelength->Fill(pid->getWavelength());
    radiator_histos->m_mc_rindex->Fill(pid->getRindex());

    // find the PDG hypothesis with the highest weight
    float max_weight     = -1000;
    int   pdg_max_weight = 0;
    for(const auto& hyp : pid->getHypotheses()) {
      if(hyp.weight > max_weight) {
        max_weight     = hyp.weight;
        pdg_max_weight = hyp.pdg;
      }
    }
    std::string pdg_max_weight_str = "UNKNOWN";
    if(pdg_max_weight!=0 && !std::isnan(pdg_max_weight))
      pdg_max_weight_str = std::to_string(pdg_max_weight);
    m_log->trace(" Highest weight is {} for PDG {} (string='{}')", max_weight, pdg_max_weight, pdg_max_weight_str);
    radiator_histos->m_highest_weight_dist->Fill(pdg_max_weight_str.c_str(), 1);
    radiator_histos->m_highest_weight_vs_p->Fill(charged_particle_momentum, pdg_max_weight_str.c_str(), 1);

  }
}


// AlgorithmFinish
//---------------------------------------------------------------------------
void eicrecon::IrtCherenkovParticleIDAnalysis::AlgorithmFinish() {
  // delete transient histograms, so they don't get written
  delete m_nphot_vs_p__transient;
}
