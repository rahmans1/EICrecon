
#ifndef CalorimeterHit_factory_GFHCALRecHits_h_
#define CalorimeterHit_factory_GFHCALRecHits_h_

#include <JANA/JFactoryT.h>

#include <algorithms/calorimetry/CalorimeterHitReco.h>
#include <services/log/Log_service.h>
#include <extensions/spdlog/SpdlogExtensions.h>

class CalorimeterHit_factory_GFHCALRecHits : public JFactoryT<edm4eic::CalorimeterHit>, CalorimeterHitReco {

public:
    //------------------------------------------
    // Constructor
    CalorimeterHit_factory_GFHCALRecHits(){
        SetTag("GFHCALRecHits");
        m_log = japp->GetService<Log_service>()->logger(GetTag());
    }

    //------------------------------------------
    // Init
    void Init() override{
        auto app = GetApplication();

        m_input_tag = "GFHCALRawHits";

        // digitization settings, must be consistent with digi class
        m_capADC=131072;//2^17
        m_dyRangeADC=1 * dd4hep::GeV;//{this, "dynamicRangeADC", 100. * dd4hep::MeV};
        m_pedMeanADC=20;//{this, "pedestalMean", 400};
        m_pedSigmaADC=0.8;//{this, "pedestalSigma", 3.2};
        m_resolutionTDC=10 * dd4hep::picosecond;//{this, "resolutionTDC", 10 * ps};

        // zero suppression values
        m_thresholdFactor=1.0;//{this, "thresholdFactor", 0.0};
        m_thresholdValue=3.0;//{this, "thresholdValue", 0.0};

        // energy correction with sampling fraction
        m_sampFrac=0.033;//{this, "samplingFraction", 1.0};
        m_sampFracLayer[0]=0.016;
        for (int i = 1; i < 13; i++) m_sampFracLayer[i]=0.031;
        // geometry service to get ids, ignored if no names provided
        m_geoSvcName="geoServiceName";
        m_readout="GFHCALHits";       
        m_layerField="";              
        m_sectorField="";             

        m_localDetElement="";         
        u_localDetFields={};          

//        app->SetDefaultParameter("HCAL:tag",              m_input_tag);
        app->SetDefaultParameter("HCAL:GFHCALRecHits:capacityADC",      m_capADC);
        app->SetDefaultParameter("HCAL:GFHCALRecHits:dynamicRangeADC",  m_dyRangeADC);
        app->SetDefaultParameter("HCAL:GFHCALRecHits:pedestalMean",     m_pedMeanADC);
        app->SetDefaultParameter("HCAL:GFHCALRecHits:pedestalSigma",    m_pedSigmaADC);
        app->SetDefaultParameter("HCAL:GFHCALRecHits:resolutionTDC",    m_resolutionTDC);
        app->SetDefaultParameter("HCAL:GFHCALRecHits:thresholdFactor",  m_thresholdFactor);
        app->SetDefaultParameter("HCAL:GFHCALRecHits:thresholdValue",   m_thresholdValue);
        app->SetDefaultParameter("HCAL:GFHCALRecHits:samplingFraction", m_sampFrac);
        app->SetDefaultParameter("HCAL:GFHCALRecHits:geoServiceName",   m_geoSvcName);
        app->SetDefaultParameter("HCAL:GFHCALRecHits:readout",          m_readout);
        app->SetDefaultParameter("HCAL:GFHCALRecHits:layerField",       m_layerField);
        app->SetDefaultParameter("HCAL:GFHCALRecHits:sectorField",      m_sectorField);
        app->SetDefaultParameter("HCAL:GFHCALRecHits:localDetElement",  m_localDetElement);
        app->SetDefaultParameter("HCAL:GFHCALRecHits:localDetFields",   u_localDetFields);
        m_geoSvc = app->template GetService<JDD4hep_service>(); // TODO: implement named geometry service?

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
        rawhits = event->Get<edm4hep::RawCalorimeterHit>(m_input_tag);

        // Call Process for generic algorithm
        AlgorithmProcess();

        // Hand owner of algorithm objects over to JANA
        Set(hits);
        hits.clear(); // not really needed, but better to not leave dangling pointers around
    }

};

#endif // CalorimeterHit_factory_GFHCALRecHits_h_
