// Created by Dmitry Romanov
// Subject to the terms in the LICENSE file found in the top-level directory.
//

#pragma once

#include <TH1F.h>
#include <TH3F.h>
#include <TH2F.h>
#include <JANA/JApplication.h>

class TrackingOccupancyAnalysis {

public:
    void init(JApplication *app, TDirectory *plugin_tdir);

    void process(const std::shared_ptr<const JEvent> &event);

private:

    /// This is edm4hep::SimTrackerHits names of different detector readouts
    std::vector<std::string> m_data_names = {
            "SiBarrelHits",         // Barrel Tracker
            "MPGDBarrelHits",       // MPGD
            "VertexBarrelHits",     // Vertex
            "TrackerEndcapHits",    // End Cap tracker
            "TOFEndcapHits",        // End Cap TOF
            "TOFBarrelHits",        // Barrel TOF
    };

    /// Hits count histogram for each hits readout name
    std::vector<std::shared_ptr<TH1F>> m_hits_count_hists;

    /// Hits occupancy histogram for each hits readout name
    std::vector<std::shared_ptr<TH2F>> m_hits_occup_hists;

    /// Total occupancy of all m_data_names
    TH2F * m_total_occup_th2;                 /// MC Particles px,py
};
