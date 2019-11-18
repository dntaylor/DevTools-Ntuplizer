#include "DevTools/Ntuplizer/interface/MuonTree.h"

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/EcalDetId/interface/ESDetId.h"

MuonTree::MuonTree(const edm::ParameterSet &iConfig) :
    muonsToken_(consumes<edm::View<reco::Muon>>(iConfig.getParameter<edm::InputTag>("muonSrc"))),
    genParticlesToken_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genSrc"))),
    verticesToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
    simInfo_(consumes<edm::ValueMap<reco::MuonSimInfo>>(iConfig.getParameter<edm::InputTag>("muonSimInfo")))
{
    // Declare use of TFileService
    usesResource("TFileService");

    edm::Service<TFileService> FS;

    // create tree_
    tree_ = FS->make<TTree>("MuonTree", "MuonTree");

    // now build the branches
    tree_->Branch("num_vertices", &num_vertices_);

    // add branches
    tree_->Branch("muon_pt",     &muon_pt_);
    tree_->Branch("muon_p",      &muon_p_);
    tree_->Branch("muon_eta",    &muon_eta_);
    tree_->Branch("muon_phi",    &muon_phi_);
    tree_->Branch("muon_energy", &muon_energy_);
    tree_->Branch("muon_mass",   &muon_mass_);
    tree_->Branch("muon_charge", &muon_charge_);

    tree_->Branch("muon_isPFMuon",         &muon_isPFMuon_);
    tree_->Branch("muon_isTrackerMuon",    &muon_isTrackerMuon_);
    tree_->Branch("muon_isGlobalMuon",     &muon_isGlobalMuon_);
    tree_->Branch("muon_isStandAloneMuon", &muon_isStandAloneMuon_);
    tree_->Branch("muon_isCaloMuon",       &muon_isCaloMuon_);
    tree_->Branch("muon_CutBasedIdLoose",         &muon_CutBasedIdLoose_);
    tree_->Branch("muon_CutBasedIdMedium",        &muon_CutBasedIdMedium_);
    tree_->Branch("muon_CutBasedIdMediumPrompt",  &muon_CutBasedIdMediumPrompt_);
    tree_->Branch("muon_CutBasedIdTight",         &muon_CutBasedIdTight_);
    tree_->Branch("muon_PFIsoLoose",              &muon_PFIsoLoose_);
    tree_->Branch("muon_PFIsoTight",              &muon_PFIsoTight_);
    tree_->Branch("muon_SoftCutBasedId",          &muon_SoftCutBasedId_);

    tree_->Branch("muon_isolationR03_nTracks",   &muon_isolationR03_nTracks_);
    tree_->Branch("muon_isolationR03_sumPt",     &muon_isolationR03_sumPt_);

    tree_->Branch("muon_innerTrack_pt",          &muon_innerTrack_pt_);
    tree_->Branch("muon_innerTrack_p",           &muon_innerTrack_p_);
    tree_->Branch("muon_innerTrack_eta",         &muon_innerTrack_eta_);
    tree_->Branch("muon_innerTrack_phi",         &muon_innerTrack_phi_);
    tree_->Branch("muon_innerTrack_qoverp",      &muon_innerTrack_qoverp_);
    tree_->Branch("muon_innerTrack_qoverpError", &muon_innerTrack_qoverpError_);
    tree_->Branch("muon_innerTrack_validFraction",                            &muon_innerTrack_validFraction_);
    tree_->Branch("muon_innerTrack_highPurity",                               &muon_innerTrack_highPurity_);
    tree_->Branch("muon_innerTrack_hitPattern_trackerLayersWithMeasurement",  &muon_innerTrack_hitPattern_trackerLayersWithMeasurement_);
    tree_->Branch("muon_innerTrack_hitPattern_pixelLayersWithMeasurement",    &muon_innerTrack_hitPattern_pixelLayersWithMeasurement_);

    tree_->Branch("muon_caloCompatibility", &muon_caloCompatibility_);
    tree_->Branch("muon_calEnergy_ecal_time", &muon_calEnergy_ecal_time_);
    tree_->Branch("muon_calEnergy_em",        &muon_calEnergy_em_);
    tree_->Branch("muon_calEnergy_emMax",     &muon_calEnergy_emMax_);
    tree_->Branch("muon_calEnergy_emS25",     &muon_calEnergy_emS25_);
    tree_->Branch("muon_calEnergy_emS9",      &muon_calEnergy_emS9_);
    tree_->Branch("muon_calEnergy_had",       &muon_calEnergy_had_);
    tree_->Branch("muon_calEnergy_hadMax",    &muon_calEnergy_hadMax_);
    tree_->Branch("muon_calEnergy_hadS9",     &muon_calEnergy_hadS9_);
    tree_->Branch("muon_calEnergy_hcal_time", &muon_calEnergy_hcal_time_);
    tree_->Branch("muon_calEnergy_ho",        &muon_calEnergy_ho_);
    tree_->Branch("muon_calEnergy_hoS9",      &muon_calEnergy_hoS9_);
    tree_->Branch("muon_calEnergy_tower",     &muon_calEnergy_tower_);
    tree_->Branch("muon_calEnergy_towerS9",   &muon_calEnergy_towerS9_);
    tree_->Branch("muon_calEnergy_ecal_ieta", &muon_calEnergy_ecal_ieta_);
    tree_->Branch("muon_calEnergy_ecal_iphi", &muon_calEnergy_ecal_iphi_);
    tree_->Branch("muon_calEnergy_hcal_ieta", &muon_calEnergy_hcal_ieta_);
    tree_->Branch("muon_calEnergy_hcal_iphi", &muon_calEnergy_hcal_iphi_);
    tree_->Branch("muon_calEnergy_hcal_depth",&muon_calEnergy_hcal_depth_);
    tree_->Branch("muon_calEnergy_crossedHadRecHits_ieta",   &muon_calEnergy_crossedHadRecHits_ieta_);
    tree_->Branch("muon_calEnergy_crossedHadRecHits_iphi",   &muon_calEnergy_crossedHadRecHits_iphi_);
    tree_->Branch("muon_calEnergy_crossedHadRecHits_depth",  &muon_calEnergy_crossedHadRecHits_depth_);
    tree_->Branch("muon_calEnergy_crossedHadRecHits_energy", &muon_calEnergy_crossedHadRecHits_energy_);
    tree_->Branch("muon_calEnergy_crossedHadRecHits_time",   &muon_calEnergy_crossedHadRecHits_time_);
    tree_->Branch("muon_calEnergy_crossedHadRecHits_chi2",   &muon_calEnergy_crossedHadRecHits_chi2_);

    tree_->Branch("muon_gen_matches_muon",   &muon_gen_matches_muon_);
    tree_->Branch("muon_gen_matches_pion",   &muon_gen_matches_pion_);
    tree_->Branch("muon_gen_deltaR",         &muon_gen_deltaR_);
    tree_->Branch("muon_gen_pt",             &muon_gen_pt_);
    tree_->Branch("muon_gen_eta",            &muon_gen_eta_);
    tree_->Branch("muon_gen_phi",            &muon_gen_phi_);
    tree_->Branch("muon_gen_mass",           &muon_gen_mass_);

    tree_->Branch("muon_gen_sim_primaryClass",  &muon_gen_sim_primaryClass_);
    tree_->Branch("muon_gen_sim_extendedClass", &muon_gen_sim_extendedClass_);
    tree_->Branch("muon_gen_sim_flavour",       &muon_gen_sim_flavour_);
    tree_->Branch("muon_gen_sim_pdgId",         &muon_gen_sim_pdgId_);
    tree_->Branch("muon_gen_sim_pt",            &muon_gen_sim_pt_);
    tree_->Branch("muon_gen_sim_eta",           &muon_gen_sim_eta_);
    tree_->Branch("muon_gen_sim_phi",           &muon_gen_sim_phi_);
    tree_->Branch("muon_gen_sim_mass",          &muon_gen_sim_mass_);
    tree_->Branch("muon_gen_sim_tpAssoQuality", &muon_gen_sim_tpAssoQuality_);
}

MuonTree::~MuonTree() { }

void MuonTree::beginJob() { }

void MuonTree::endJob() { }

void MuonTree::analyze(const edm::Event &iEvent, const edm::EventSetup &iSetup) {
    edm::Handle<edm::View<reco::Muon>> muons;
    iEvent.getByToken(muonsToken_, muons);

    edm::Handle<reco::GenParticleCollection> genParticles;
    iEvent.getByToken(genParticlesToken_, genParticles);

    edm::Handle<reco::VertexCollection> vertices;
    iEvent.getByToken(verticesToken_, vertices);

    edm::Handle<edm::ValueMap<reco::MuonSimInfo>> simInfo;
    bool simInfoIsAvailable = iEvent.getByToken(simInfo_, simInfo);

    num_vertices_ = vertices->size();

    unsigned int idx = 0;
    for (const auto m: *muons) {
        muon_pt_ = m.pt();
        muon_p_ = m.p();
        muon_eta_ = m.eta();
        muon_phi_ = m.phi();
        muon_energy_ = m.energy();
        muon_mass_ = m.mass();
        muon_charge_ = m.charge();

        muon_isPFMuon_ = m.isPFMuon();
        muon_isTrackerMuon_ = m.isTrackerMuon();
        muon_isGlobalMuon_ = m.isGlobalMuon();
        muon_isStandAloneMuon_ = m.isStandAloneMuon();
        muon_isCaloMuon_ = m.isCaloMuon();
        muon_CutBasedIdLoose_ = m.passed(reco::Muon::CutBasedIdLoose);
        muon_CutBasedIdMedium_ = m.passed(reco::Muon::CutBasedIdMedium);
        muon_CutBasedIdMediumPrompt_ = m.passed(reco::Muon::CutBasedIdMediumPrompt);
        muon_CutBasedIdTight_ = m.passed(reco::Muon::CutBasedIdTight);
        muon_PFIsoLoose_ = m.passed(reco::Muon::PFIsoLoose);
        muon_PFIsoTight_ = m.passed(reco::Muon::PFIsoTight);
        muon_SoftCutBasedId_ = m.passed(reco::Muon::SoftCutBasedId);

        muon_isolationR03_nTracks_ = m.isolationR03().nTracks;
        muon_isolationR03_sumPt_ = m.isolationR03().sumPt;

        if (m.innerTrack().isNonnull()) {
            muon_innerTrack_pt_ = m.innerTrack()->pt();
            muon_innerTrack_p_ = m.innerTrack()->p();
            muon_innerTrack_eta_ = m.innerTrack()->eta();
            muon_innerTrack_phi_ = m.innerTrack()->phi();
            muon_innerTrack_qoverp_ = m.innerTrack()->qoverp();
            muon_innerTrack_qoverpError_ = m.innerTrack()->qoverpError();
            muon_innerTrack_validFraction_ = m.innerTrack()->validFraction();
            muon_innerTrack_highPurity_ = m.innerTrack()->quality(reco::TrackBase::highPurity);
            muon_innerTrack_hitPattern_trackerLayersWithMeasurement_ = m.innerTrack()->hitPattern().trackerLayersWithMeasurement();
            muon_innerTrack_hitPattern_pixelLayersWithMeasurement_ = m.innerTrack()->hitPattern().pixelLayersWithMeasurement();
        } else {
            muon_innerTrack_pt_ = 0;
            muon_innerTrack_p_ = 0;
            muon_innerTrack_eta_ = 0;
            muon_innerTrack_phi_ = 0;
            muon_innerTrack_qoverp_ = 0;
            muon_innerTrack_qoverpError_ = 0;
            muon_innerTrack_validFraction_ = 0;
            muon_innerTrack_highPurity_ = 0;
            muon_innerTrack_hitPattern_trackerLayersWithMeasurement_ = 0;
            muon_innerTrack_hitPattern_pixelLayersWithMeasurement_ = 0;
        }

        muon_caloCompatibility_ = m.caloCompatibility();
        muon_calEnergy_ecal_time_ = m.calEnergy().ecal_time;
        muon_calEnergy_em_ = m.calEnergy().em;
        muon_calEnergy_emMax_ = m.calEnergy().emMax;
        muon_calEnergy_emS25_ = m.calEnergy().emS25;
        muon_calEnergy_emS9_ = m.calEnergy().emS9;
        muon_calEnergy_had_ = m.calEnergy().had;
        muon_calEnergy_hadMax_ = m.calEnergy().hadMax;
        muon_calEnergy_hadS9_ = m.calEnergy().hadS9;
        muon_calEnergy_hcal_time_ = m.calEnergy().hcal_time;
        muon_calEnergy_ho_ = m.calEnergy().ho;
        muon_calEnergy_hoS9_ = m.calEnergy().hoS9;
        muon_calEnergy_tower_ = m.calEnergy().tower;
        muon_calEnergy_towerS9_ = m.calEnergy().towerS9;
        if (m.calEnergy().ecal_id.subdetId() == EBDetId::Subdet) {
            muon_calEnergy_ecal_ieta_ = ((EBDetId)m.calEnergy().ecal_id).ieta();
            muon_calEnergy_ecal_iphi_ = ((EBDetId)m.calEnergy().ecal_id).iphi();
        } else if (m.calEnergy().ecal_id.subdetId() == EEDetId::Subdet) {
            muon_calEnergy_ecal_ieta_ = ((EEDetId)m.calEnergy().ecal_id).ix();
            muon_calEnergy_ecal_iphi_ = ((EEDetId)m.calEnergy().ecal_id).iy();
        } else if (m.calEnergy().ecal_id.subdetId() == ESDetId::Subdet) {
            muon_calEnergy_ecal_ieta_ = ((ESDetId)m.calEnergy().ecal_id).six();
            muon_calEnergy_ecal_iphi_ = ((ESDetId)m.calEnergy().ecal_id).siy();
        } else {
            muon_calEnergy_ecal_ieta_ = 0;
            muon_calEnergy_ecal_iphi_ = 0;
        }
        muon_calEnergy_hcal_ieta_ =  ((HcalDetId)m.calEnergy().hcal_id).ieta();
        muon_calEnergy_hcal_iphi_ =  ((HcalDetId)m.calEnergy().hcal_id).iphi();
        muon_calEnergy_hcal_depth_ = ((HcalDetId)m.calEnergy().hcal_id).depth();
        muon_calEnergy_crossedHadRecHits_ieta_.clear();
        muon_calEnergy_crossedHadRecHits_iphi_.clear();
        muon_calEnergy_crossedHadRecHits_depth_.clear();
        muon_calEnergy_crossedHadRecHits_energy_.clear();
        muon_calEnergy_crossedHadRecHits_time_.clear();
        muon_calEnergy_crossedHadRecHits_chi2_.clear();
        for (auto it: m.calEnergy().crossedHadRecHits) {
            muon_calEnergy_crossedHadRecHits_ieta_.push_back(it.detId.ieta());
            muon_calEnergy_crossedHadRecHits_iphi_.push_back(it.detId.iphi());
            muon_calEnergy_crossedHadRecHits_depth_.push_back(it.detId.depth());
            muon_calEnergy_crossedHadRecHits_energy_.push_back(it.energy);
            muon_calEnergy_crossedHadRecHits_time_.push_back(it.time);
            muon_calEnergy_crossedHadRecHits_chi2_.push_back(it.chi2);
        }

        // For now do deltaR matching, next iteration do sim hit matching
        int isGenMuon = 0;
        int isGenPion = 0;
        float gen_DR = 9.;
        float gen_pt = 0.;
        float gen_eta = 0.;
        float gen_phi = 0.;
        float gen_mass = 0.;
        for (const auto g: *genParticles) {
            if (g.status()!=1) continue;
            if (deltaR(m,g)<gen_DR) {
                gen_DR = deltaR(m,g);
                gen_pt = g.pt();
                gen_eta = g.eta();
                gen_phi = g.phi();
                gen_mass = g.mass();
                isGenMuon = abs(g.pdgId())==13;
                isGenPion = abs(g.pdgId())==211;
            }
        }
        muon_gen_matches_muon_ = isGenMuon;
        muon_gen_matches_pion_ = isGenPion;
        muon_gen_deltaR_ = gen_DR;
        muon_gen_pt_ = gen_pt;
        muon_gen_eta_ = gen_eta;
        muon_gen_phi_ = gen_phi;
        muon_gen_mass_ = gen_mass;

        // Sim hit matching
        edm::RefToBase<reco::Muon> muonRef = muons->refAt(idx);
        reco::CandidateBaseRef muonBaseRef(muonRef);
        if (simInfoIsAvailable) {
          const auto& msi = (*simInfo)[muonBaseRef];
          muon_gen_sim_primaryClass_ = msi.primaryClass;
          muon_gen_sim_extendedClass_ = msi.extendedClass;
          muon_gen_sim_flavour_ = msi.flavour;
          muon_gen_sim_pdgId_ = msi.pdgId;
          muon_gen_sim_pt_ = msi.p4.pt();
          muon_gen_sim_eta_ = msi.p4.eta();
          muon_gen_sim_phi_ = msi.p4.phi();
          muon_gen_sim_mass_ = msi.p4.mass();
          muon_gen_sim_tpAssoQuality_ = msi.tpAssoQuality;
        } else {
          muon_gen_sim_primaryClass_ = 0;
          muon_gen_sim_extendedClass_ = 0;
          muon_gen_sim_flavour_ = 0;
          muon_gen_sim_pdgId_ = 0;
          muon_gen_sim_pt_ = 0;
          muon_gen_sim_eta_ = 0;
          muon_gen_sim_phi_ = 0;
          muon_gen_sim_mass_ = 0;
          muon_gen_sim_tpAssoQuality_ = 0;
        }
        
        if (m.p()>1.0) tree_->Fill();

        idx++;
    }

}
