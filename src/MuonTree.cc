#include "DevTools/Ntuplizer/interface/MuonTree.h"

MuonTree::MuonTree(const edm::ParameterSet &iConfig) :
    muonsToken_(consumes<edm::View<reco::Muon>>(iConfig.getParameter<edm::InputTag>("muonSrc"))),
    genParticlesToken_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genSrc"))),
    simInfo_(consumes<edm::ValueMap<reco::MuonSimInfo>>(iConfig.getParameter<edm::InputTag>("muonSimInfo")))
{
    // Declare use of TFileService
    usesResource("TFileService");

    edm::Service<TFileService> FS;

    // create tree_
    tree_ = FS->make<TTree>("MuonTree", "MuonTree");

    // now build the branches

    // add branches
    tree_->Branch("muon_pt",     &muon_pt_);
    tree_->Branch("muon_p",      &muon_p_);
    tree_->Branch("muon_eta",    &muon_eta_);
    tree_->Branch("muon_phi",    &muon_phi_);
    tree_->Branch("muon_energy", &muon_energy_);
    tree_->Branch("muon_mass",   &muon_mass_);

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

    edm::Handle<edm::ValueMap<reco::MuonSimInfo>> simInfo;
    bool simInfoIsAvailable = iEvent.getByToken(simInfo_, simInfo);

    bool keep = (muons->size()>0);

    muon_pt_.clear();
    muon_p_.clear();
    muon_eta_.clear();
    muon_phi_.clear();
    muon_energy_.clear();
    muon_mass_.clear();

    muon_isPFMuon_.clear();
    muon_isTrackerMuon_.clear();
    muon_isGlobalMuon_.clear();
    muon_isStandAloneMuon_.clear();
    muon_isCaloMuon_.clear();
    muon_CutBasedIdLoose_.clear();
    muon_CutBasedIdMedium_.clear();
    muon_CutBasedIdMediumPrompt_.clear();
    muon_CutBasedIdTight_.clear();
    muon_PFIsoLoose_.clear();
    muon_PFIsoTight_.clear();
    muon_SoftCutBasedId_.clear();

    muon_innerTrack_validFraction_.clear();
    muon_innerTrack_highPurity_.clear();
    muon_innerTrack_hitPattern_trackerLayersWithMeasurement_.clear();
    muon_innerTrack_hitPattern_pixelLayersWithMeasurement_.clear();

    muon_caloCompatibility_.clear();
    muon_calEnergy_ecal_time_.clear();
    muon_calEnergy_em_.clear();
    muon_calEnergy_emMax_.clear();
    muon_calEnergy_emS25_.clear();
    muon_calEnergy_emS9_.clear();
    muon_calEnergy_had_.clear();
    muon_calEnergy_hadMax_.clear();
    muon_calEnergy_hadS9_.clear();
    muon_calEnergy_hcal_time_.clear();
    muon_calEnergy_ho_.clear();
    muon_calEnergy_hoS9_.clear();
    muon_calEnergy_tower_.clear();
    muon_calEnergy_towerS9_.clear();
    muon_calEnergy_crossedHadRecHits_ieta_.clear();
    muon_calEnergy_crossedHadRecHits_iphi_.clear();
    muon_calEnergy_crossedHadRecHits_depth_.clear();
    muon_calEnergy_crossedHadRecHits_energy_.clear();
    muon_calEnergy_crossedHadRecHits_time_.clear();
    muon_calEnergy_crossedHadRecHits_chi2_.clear();

    muon_gen_matches_muon_.clear();
    muon_gen_matches_pion_.clear();
    muon_gen_deltaR_.clear();
    muon_gen_pt_.clear();
    muon_gen_eta_.clear();
    muon_gen_phi_.clear();
    muon_gen_mass_.clear();

    muon_gen_sim_primaryClass_.clear();
    muon_gen_sim_extendedClass_.clear();
    muon_gen_sim_flavour_.clear();
    muon_gen_sim_pdgId_.clear();
    muon_gen_sim_pt_.clear();
    muon_gen_sim_eta_.clear();
    muon_gen_sim_phi_.clear();
    muon_gen_sim_mass_.clear();
    muon_gen_sim_tpAssoQuality_.clear();

    unsigned int idx = 0;
    for (const auto m: *muons) {
        muon_pt_.push_back(m.pt());
        muon_p_.push_back(m.p());
        muon_eta_.push_back(m.eta());
        muon_phi_.push_back(m.phi());
        muon_energy_.push_back(m.energy());
        muon_mass_.push_back(m.mass());

        muon_isPFMuon_.push_back(m.isPFMuon());
        muon_isTrackerMuon_.push_back(m.isTrackerMuon());
        muon_isGlobalMuon_.push_back(m.isGlobalMuon());
        muon_isStandAloneMuon_.push_back(m.isStandAloneMuon());
        muon_isCaloMuon_.push_back(m.isCaloMuon());
        muon_CutBasedIdLoose_.push_back(m.passed(reco::Muon::CutBasedIdLoose));
        muon_CutBasedIdMedium_.push_back(m.passed(reco::Muon::CutBasedIdMedium));
        muon_CutBasedIdMediumPrompt_.push_back(m.passed(reco::Muon::CutBasedIdMediumPrompt));
        muon_CutBasedIdTight_.push_back(m.passed(reco::Muon::CutBasedIdTight));
        muon_PFIsoLoose_.push_back(m.passed(reco::Muon::PFIsoLoose));
        muon_PFIsoTight_.push_back(m.passed(reco::Muon::PFIsoTight));
        muon_SoftCutBasedId_.push_back(m.passed(reco::Muon::SoftCutBasedId));

        if (m.innerTrack().isNonnull()) {
            muon_innerTrack_validFraction_.push_back(m.innerTrack()->validFraction());
            muon_innerTrack_highPurity_.push_back(m.innerTrack()->quality(reco::TrackBase::highPurity));
            muon_innerTrack_hitPattern_trackerLayersWithMeasurement_.push_back(m.innerTrack()->hitPattern().trackerLayersWithMeasurement());
            muon_innerTrack_hitPattern_pixelLayersWithMeasurement_.push_back(m.innerTrack()->hitPattern().pixelLayersWithMeasurement());
        } else {
            muon_innerTrack_validFraction_.push_back(-1);
            muon_innerTrack_highPurity_.push_back(-1);
            muon_innerTrack_hitPattern_trackerLayersWithMeasurement_.push_back(-1);
            muon_innerTrack_hitPattern_pixelLayersWithMeasurement_.push_back(-1);
        }

        muon_caloCompatibility_.push_back(m.caloCompatibility());
        muon_calEnergy_ecal_time_.push_back(m.calEnergy().ecal_time);
        muon_calEnergy_em_.push_back(m.calEnergy().em);
        muon_calEnergy_emMax_.push_back(m.calEnergy().emMax);
        muon_calEnergy_emS25_.push_back(m.calEnergy().emS25);
        muon_calEnergy_emS9_.push_back(m.calEnergy().emS9);
        muon_calEnergy_had_.push_back(m.calEnergy().had);
        muon_calEnergy_hadMax_.push_back(m.calEnergy().hadMax);
        muon_calEnergy_hadS9_.push_back(m.calEnergy().hadS9);
        muon_calEnergy_hcal_time_.push_back(m.calEnergy().hcal_time);
        muon_calEnergy_ho_.push_back(m.calEnergy().ho);
        muon_calEnergy_hoS9_.push_back(m.calEnergy().hoS9);
        muon_calEnergy_tower_.push_back(m.calEnergy().tower);
        muon_calEnergy_towerS9_.push_back(m.calEnergy().towerS9);
        std::vector<int> v_ieta;
        std::vector<int> v_iphi;
        std::vector<int> v_depth;
        std::vector<float> v_energy;
        std::vector<float> v_time;
        std::vector<float> v_chi2;
        for (auto it: m.calEnergy().crossedHadRecHits) {
            v_ieta.push_back(it.detId.ieta());
            v_iphi.push_back(it.detId.iphi());
            v_depth.push_back(it.detId.depth());
            v_energy.push_back(it.energy);
            v_time.push_back(it.time);
            v_chi2.push_back(it.chi2);
        }
        muon_calEnergy_crossedHadRecHits_ieta_.push_back(v_ieta);
        muon_calEnergy_crossedHadRecHits_iphi_.push_back(v_iphi);
        muon_calEnergy_crossedHadRecHits_depth_.push_back(v_depth);
        muon_calEnergy_crossedHadRecHits_energy_.push_back(v_energy);
        muon_calEnergy_crossedHadRecHits_time_.push_back(v_time);
        muon_calEnergy_crossedHadRecHits_chi2_.push_back(v_chi2);

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
        muon_gen_matches_muon_.push_back(isGenMuon);
        muon_gen_matches_pion_.push_back(isGenPion);
        muon_gen_deltaR_.push_back(gen_DR);
        muon_gen_pt_.push_back(gen_pt);
        muon_gen_eta_.push_back(gen_eta);
        muon_gen_phi_.push_back(gen_phi);
        muon_gen_mass_.push_back(gen_mass);

        // Sim hit matching
        edm::RefToBase<reco::Muon> muonRef = muons->refAt(idx);
        reco::CandidateBaseRef muonBaseRef(muonRef);
        if (simInfoIsAvailable) {
          const auto& msi = (*simInfo)[muonBaseRef];
          muon_gen_sim_primaryClass_.push_back(msi.primaryClass);
          muon_gen_sim_extendedClass_.push_back(msi.extendedClass);
          muon_gen_sim_flavour_.push_back(msi.flavour);
          muon_gen_sim_pdgId_.push_back(msi.pdgId);
          muon_gen_sim_pt_.push_back(msi.p4.pt());
          muon_gen_sim_eta_.push_back(msi.p4.eta());
          muon_gen_sim_phi_.push_back(msi.p4.phi());
          muon_gen_sim_mass_.push_back(msi.p4.mass());
          muon_gen_sim_tpAssoQuality_.push_back(msi.tpAssoQuality);
        } else {
          muon_gen_sim_primaryClass_.push_back(0);
          muon_gen_sim_extendedClass_.push_back(0);
          muon_gen_sim_flavour_.push_back(0);
          muon_gen_sim_pdgId_.push_back(0);
          muon_gen_sim_pt_.push_back(0);
          muon_gen_sim_eta_.push_back(0);
          muon_gen_sim_phi_.push_back(0);
          muon_gen_sim_mass_.push_back(0);
          muon_gen_sim_tpAssoQuality_.push_back(0);
        }

        idx++;
    }

    if (keep)
        tree_->Fill();
}
