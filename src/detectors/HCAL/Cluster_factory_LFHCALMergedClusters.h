// Copyright 2022, Thomas Britton (mod Peter Steinberg)
// Subject to the terms in the LICENSE file found in the top-level directory.
//

#pragma once

#include <random>

#include <JANA/JFactoryT.h>
#include <services/geometry/dd4hep/JDD4hep_service.h>
#include <algorithms/calorimetry/CalorimeterClusterMerger.h>
#include <services/log/Log_service.h>
#include <extensions/spdlog/SpdlogExtensions.h>



// Dummy factory for JFactoryGeneratorT
class Association_factory_LFHCALMergedClustersAssociations : public JFactoryT<edm4eic::MCRecoClusterParticleAssociation> {

public:
    //------------------------------------------
    // Constructor
    Association_factory_LFHCALMergedClustersAssociations(){
        SetTag("LFHCALMergedClustersAssociations");
    }
};



class Cluster_factory_LFHCALMergedClusters : public JFactoryT<edm4eic::Cluster>, CalorimeterClusterMerger {

public:
    //------------------------------------------
    // Constructor
    Cluster_factory_LFHCALMergedClusters(){
        SetTag("LFHCALMergedClusters");
        m_log = japp->GetService<Log_service>()->logger(GetTag());
    }

    //------------------------------------------
    // Init
    void Init() override{
        auto app = GetApplication();
        //-------- Configuration Parameters ------------
        m_input_tag="LFHCALClusters";
        m_inputAssociations_tag="LFHCALClustersAssociations";

        app->SetDefaultParameter("HCAL:LFHCALMergedClusters:input_tag",      m_input_tag, "Name of input collection to use");
        app->SetDefaultParameter("HCAL:LFHCALMergedClusters:inputAssociations_tag",      m_inputAssociations_tag, "Name of input associations collection to use");

        AlgorithmInit(m_log);
    }

    //------------------------------------------
    // ChangeRun
    void ChangeRun(const std::shared_ptr<const JEvent> &event) override{
        AlgorithmChangeRun();
    }

    //------------------------------------------
    // Process
    void Process(const std::shared_ptr<const JEvent> &event) override{

        
        // Prefill inputs
        m_inputClusters=event->Get<edm4eic::Cluster>(m_input_tag);
        m_inputAssociations=event->Get<edm4eic::MCRecoClusterParticleAssociation>(m_inputAssociations_tag); 

        // Call Process for generic algorithm
        AlgorithmProcess();

        //outputs
        // Hand owner of algorithm objects over to JANA
        Set(m_outputClusters);
        event->Insert(m_outputAssociations, "LFHCALMergedClustersAssociations");
        m_outputClusters.clear(); // not really needed, but better to not leave dangling pointers around
        m_outputAssociations.clear();
    }

private:
    // Name of input data type (collection)
    std::string              m_input_tag;
    std::string              m_inputAssociations_tag;
};

