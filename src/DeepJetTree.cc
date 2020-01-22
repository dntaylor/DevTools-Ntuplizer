#include "DevTools/Ntuplizer/interface/DeepJetTree.h"

DeepJetTree::DeepJetTree(const edm::ParameterSet &iConfig) :
    jetsToken_(consumes<edm::View<pat::Jet>>(iConfig.getParameter<edm::InputTag>("jetSrc"))),
    genParticlesToken_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genSrc"))),
    verticesToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices")))
{

    // Declare use of TFileService
    usesResource("TFileService");

    edm::Service<TFileService> FS;

    // create tree_
    tree_ = FS->make<TTree>("DeepJetTree", "DeepJetTree");

    // now build the branches
    tree_->Branch("num_vertices",                    &num_vertices_);

    // jet global branches
    // kinematics
    tree_->Branch("jet_pt",                          &jet_pt_);
    tree_->Branch("jet_eta",                         &jet_eta_);
    tree_->Branch("jet_phi",                         &jet_phi_);
    tree_->Branch("jet_energy",                      &jet_energy_);
    tree_->Branch("jet_mass",                        &jet_mass_);
    tree_->Branch("jet_jetCharge",                   &jet_jetCharge_);

    // energy sums/fractions/multiplicity
    tree_->Branch("jet_chargedMultiplicity",         &jet_chargedMultiplicity_);
    tree_->Branch("jet_neutralMultiplicity",         &jet_neutralMultiplicity_);

    tree_->Branch("jet_chargedHadronMultiplicity",   &jet_chargedHadronMultiplicity_);
    tree_->Branch("jet_neutralHadronMultiplicity",   &jet_neutralHadronMultiplicity_);
    tree_->Branch("jet_muonMultiplicity",            &jet_muonMultiplicity_);
    tree_->Branch("jet_electronMultiplicity",        &jet_electronMultiplicity_);
    tree_->Branch("jet_photonMultiplicity",          &jet_photonMultiplicity_);

    tree_->Branch("jet_chargedEmEnergy",             &jet_chargedEmEnergy_);
    tree_->Branch("jet_neutralEmEnergy",             &jet_neutralEmEnergy_);
    tree_->Branch("jet_chargedHadronEnergy",         &jet_chargedHadronEnergy_);
    tree_->Branch("jet_neutralHadronEnergy",         &jet_neutralHadronEnergy_);
    tree_->Branch("jet_muonEnergy",                  &jet_muonEnergy_);
    tree_->Branch("jet_electronEnergy",              &jet_electronEnergy_);
    tree_->Branch("jet_photonEnergy",                &jet_photonEnergy_);

    tree_->Branch("jet_chargedEmEnergyFraction",     &jet_chargedEmEnergyFraction_);
    tree_->Branch("jet_neutralEmEnergyFraction",     &jet_neutralEmEnergyFraction_);
    tree_->Branch("jet_chargedHadronEnergyFraction", &jet_chargedHadronEnergyFraction_);
    tree_->Branch("jet_neutralHadronEnergyFraction", &jet_neutralHadronEnergyFraction_);
    tree_->Branch("jet_muonEnergyFraction",          &jet_muonEnergyFraction_);
    tree_->Branch("jet_electronEnergyFraction",      &jet_electronEnergyFraction_);
    tree_->Branch("jet_photonEnergyFraction",        &jet_photonEnergyFraction_);

    tree_->Branch("jet_n60",                         &jet_n60_);
    tree_->Branch("jet_n90",                         &jet_n90_);

    // gen truth
    tree_->Branch("jet_hadronFlavour",               &jet_hadronFlavour_);
    tree_->Branch("jet_partonFlavour",               &jet_partonFlavour_);

    // btag scores AK4
    tree_->Branch("jet_pfJetBProbabilityBJetTags",                       &jet_pfJetBProbabilityBJetTags_);
    tree_->Branch("jet_pfJetProbabilityBJetTags",                        &jet_pfJetProbabilityBJetTags_);
    tree_->Branch("jet_pfTrackCountingHighEffBJetTags",                  &jet_pfTrackCountingHighEffBJetTags_);
    tree_->Branch("jet_pfSimpleSecondaryVertexHighEffBJetTags",          &jet_pfSimpleSecondaryVertexHighEffBJetTags_);
    tree_->Branch("jet_pfSimpleInclusiveSecondaryVertexHighEffBJetTags", &jet_pfSimpleInclusiveSecondaryVertexHighEffBJetTags_);
    tree_->Branch("jet_pfCombinedSecondaryVertexV2BJetTags",             &jet_pfCombinedSecondaryVertexV2BJetTags_);
    tree_->Branch("jet_pfCombinedInclusiveSecondaryVertexV2BJetTags",    &jet_pfCombinedInclusiveSecondaryVertexV2BJetTags_);
    tree_->Branch("jet_softPFMuonBJetTags",                              &jet_softPFMuonBJetTags_);
    tree_->Branch("jet_softPFElectronBJetTags",                          &jet_softPFElectronBJetTags_);
    tree_->Branch("jet_pfCombinedMVAV2BJetTags",                         &jet_pfCombinedMVAV2BJetTags_);
    tree_->Branch("jet_pfCombinedCvsLJetTags",                           &jet_pfCombinedCvsLJetTags_);
    tree_->Branch("jet_pfCombinedCvsBJetTags",                           &jet_pfCombinedCvsBJetTags_);
    tree_->Branch("jet_pfDeepCSVJetTags_probb",                          &jet_pfDeepCSVJetTags_probb_);
    tree_->Branch("jet_pfDeepCSVJetTags_probc",                          &jet_pfDeepCSVJetTags_probc_);
    tree_->Branch("jet_pfDeepCSVJetTags_probudsg",                       &jet_pfDeepCSVJetTags_probudsg_);
    tree_->Branch("jet_pfDeepCSVJetTags_probbb",                         &jet_pfDeepCSVJetTags_probbb_);

    // jet daughters
    tree_->Branch("jet_daughter_pt",                           &jet_daughter_pt_);
    tree_->Branch("jet_daughter_eta",                          &jet_daughter_eta_);
    tree_->Branch("jet_daughter_phi",                          &jet_daughter_phi_);
    tree_->Branch("jet_daughter_energy",                       &jet_daughter_energy_);
    tree_->Branch("jet_daughter_mass",                         &jet_daughter_mass_);
    tree_->Branch("jet_daughter_charge",                       &jet_daughter_charge_);

    tree_->Branch("jet_daughter_etaAtVtx",                     &jet_daughter_etaAtVtx_);
    tree_->Branch("jet_daughter_phiAtVtx",                     &jet_daughter_phiAtVtx_);
    tree_->Branch("jet_daughter_vertexChi2",                   &jet_daughter_vertexChi2_);
    tree_->Branch("jet_daughter_vertexNdof",                   &jet_daughter_vertexNdof_);
    tree_->Branch("jet_daughter_vx",                           &jet_daughter_vx_);
    tree_->Branch("jet_daughter_vy",                           &jet_daughter_vy_);
    tree_->Branch("jet_daughter_vz",                           &jet_daughter_vz_);
    tree_->Branch("jet_daughter_dxy",                          &jet_daughter_dxy_);
    tree_->Branch("jet_daughter_dxyError",                     &jet_daughter_dxyError_);
    tree_->Branch("jet_daughter_dz",                           &jet_daughter_dz_);
    tree_->Branch("jet_daughter_dzError",                      &jet_daughter_dzError_);
    tree_->Branch("jet_daughter_dtime",                        &jet_daughter_dtime_);
    tree_->Branch("jet_daughter_time",                         &jet_daughter_time_);
    tree_->Branch("jet_daughter_timeError",                    &jet_daughter_timeError_);

    tree_->Branch("jet_daughter_pdgId",                        &jet_daughter_pdgId_);

    tree_->Branch("jet_daughter_pixelLayersWithMeasurement",   &jet_daughter_pixelLayersWithMeasurement_);
    tree_->Branch("jet_daughter_stripLayersWithMeasurement",   &jet_daughter_stripLayersWithMeasurement_);
    tree_->Branch("jet_daughter_trackerLayersWithMeasurement", &jet_daughter_trackerLayersWithMeasurement_);
    tree_->Branch("jet_daughter_trackHighPurity",              &jet_daughter_trackHighPurity_);

    tree_->Branch("jet_daughter_caloFraction",                 &jet_daughter_caloFraction_);
    tree_->Branch("jet_daughter_hcalFraction",                 &jet_daughter_hcalFraction_);

    tree_->Branch("jet_daughter_puppiWeight",                  &jet_daughter_puppiWeight_);
    tree_->Branch("jet_daughter_puppiWeightNoLep",             &jet_daughter_puppiWeightNoLep_);

    tree_->Branch("jet_daughter_isIsolatedChargedHadron",      &jet_daughter_isIsolatedChargedHadron_);
    tree_->Branch("jet_daughter_isMuon",                       &jet_daughter_isMuon_);
    tree_->Branch("jet_daughter_isElectron",                   &jet_daughter_isElectron_);
    tree_->Branch("jet_daughter_isPhoton",                     &jet_daughter_isPhoton_);
    tree_->Branch("jet_daughter_isStandAloneMuon",             &jet_daughter_isStandAloneMuon_);
    tree_->Branch("jet_daughter_isTrackerMuon",                &jet_daughter_isTrackerMuon_);
    tree_->Branch("jet_daughter_isGlobalMuon",                 &jet_daughter_isGlobalMuon_);
    tree_->Branch("jet_daughter_isGoodEgamma",                 &jet_daughter_isGoodEgamma_);
    tree_->Branch("jet_daughter_isConvertedPhoton",            &jet_daughter_isConvertedPhoton_);
}

DeepJetTree::~DeepJetTree() { }

void DeepJetTree::beginJob() { }

void DeepJetTree::endJob() { }

void DeepJetTree::analyze(const edm::Event &iEvent, const edm::EventSetup &iSetup) {
    edm::Handle<edm::View<pat::Jet>> jets;
    iEvent.getByToken(jetsToken_, jets);

    edm::Handle<reco::GenParticleCollection> genParticles;
    iEvent.getByToken(genParticlesToken_, genParticles);

    edm::Handle<reco::VertexCollection> vertices;
    iEvent.getByToken(verticesToken_, vertices);

    num_vertices_ = vertices->size();

    unsigned int idx = 0;
    for (const auto j: *jets) {
        jet_pt_ = j.pt();
        jet_eta_ = j.eta();
        jet_phi_ = j.phi();
        jet_energy_ = j.energy();
        jet_mass_ = j.mass();
        jet_jetCharge_ = j.jetCharge();

        jet_numberOfDaughters_ = j.numberOfDaughters();
        jet_chargedMultiplicity_ = j.chargedMultiplicity();
        jet_neutralMultiplicity_ = j.neutralMultiplicity();

        jet_chargedHadronMultiplicity_ = j.chargedHadronMultiplicity();
        jet_neutralHadronMultiplicity_ = j.neutralHadronMultiplicity();
        jet_muonMultiplicity_ = j.muonMultiplicity();
        jet_electronMultiplicity_ = j.electronMultiplicity();
        jet_photonMultiplicity_ = j.photonMultiplicity();

        jet_chargedEmEnergy_ = j.chargedEmEnergy();
        jet_neutralEmEnergy_ = j.neutralEmEnergy();
        jet_chargedHadronEnergy_ = j.chargedHadronEnergy();
        jet_neutralHadronEnergy_ = j.neutralHadronEnergy();
        jet_muonEnergy_ = j.muonEnergy();
        jet_electronEnergy_ = j.electronEnergy();
        jet_photonEnergy_ = j.photonEnergy();

        jet_chargedEmEnergyFraction_ = j.chargedEmEnergyFraction();
        jet_neutralEmEnergyFraction_ = j.neutralEmEnergyFraction();
        jet_chargedHadronEnergyFraction_ = j.chargedHadronEnergyFraction();
        jet_neutralHadronEnergyFraction_ = j.neutralHadronEnergyFraction();
        jet_muonEnergyFraction_ = j.muonEnergyFraction();
        jet_electronEnergyFraction_ = j.electronEnergyFraction();
        jet_photonEnergyFraction_ = j.photonEnergyFraction();

        jet_n60_ = j.n60();
        jet_n90_ = j.n90();

        jet_hadronFlavour_ = j.hadronFlavour();
        jet_partonFlavour_ = j.partonFlavour();

        jet_pfJetBProbabilityBJetTags_ = j.bDiscriminator("pfJetBProbabilityBJetTags");
        jet_pfJetProbabilityBJetTags_ = j.bDiscriminator("pfJetProbabilityBJetTags");
        jet_pfTrackCountingHighEffBJetTags_ = j.bDiscriminator("pfTrackCountingHighEffBJetTags");
        jet_pfSimpleSecondaryVertexHighEffBJetTags_ = j.bDiscriminator("pfSimpleSecondaryVertexHighEffBJetTags");
        jet_pfSimpleInclusiveSecondaryVertexHighEffBJetTags_ = j.bDiscriminator("pfSimpleInclusiveSecondaryVertexHighEffBJetTags");
        jet_pfCombinedSecondaryVertexV2BJetTags_ = j.bDiscriminator("pfCombinedSecondaryVertexV2BJetTags");
        jet_pfCombinedInclusiveSecondaryVertexV2BJetTags_ = j.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
        jet_softPFMuonBJetTags_ = j.bDiscriminator("softPFMuonBJetTags");
        jet_softPFElectronBJetTags_ = j.bDiscriminator("softPFElectronBJetTags");
        jet_pfCombinedMVAV2BJetTags_ = j.bDiscriminator("pfCombinedMVAV2BJetTags");
        jet_pfCombinedCvsLJetTags_ = j.bDiscriminator("pfCombinedCvsLJetTags");
        jet_pfCombinedCvsBJetTags_ = j.bDiscriminator("pfCombinedCvsBJetTags");
        jet_pfDeepCSVJetTags_probb_ = j.bDiscriminator("pfDeepCSVJetTags:probb");
        jet_pfDeepCSVJetTags_probc_ = j.bDiscriminator("pfDeepCSVJetTags:probc");
        jet_pfDeepCSVJetTags_probudsg_ = j.bDiscriminator("pfDeepCSVJetTags:probudsg");
        jet_pfDeepCSVJetTags_probbb_ = j.bDiscriminator("pfDeepCSVJetTags:probbb");

        // jet daughters
        jet_daughter_pt_.clear();
        jet_daughter_eta_.clear();
        jet_daughter_phi_.clear();
        jet_daughter_energy_.clear();
        jet_daughter_mass_.clear();
        jet_daughter_charge_.clear();

        jet_daughter_etaAtVtx_.clear();
        jet_daughter_phiAtVtx_.clear();
        jet_daughter_vertexChi2_.clear();
        jet_daughter_vertexNdof_.clear();
        jet_daughter_vx_.clear();
        jet_daughter_vy_.clear();
        jet_daughter_vz_.clear();
        jet_daughter_dxy_.clear();
        jet_daughter_dxyError_.clear();
        jet_daughter_dz_.clear();
        jet_daughter_dzError_.clear();
        jet_daughter_dtime_.clear();
        jet_daughter_time_.clear();
        jet_daughter_timeError_.clear();

        jet_daughter_pdgId_.clear();

        jet_daughter_pixelLayersWithMeasurement_.clear();
        jet_daughter_stripLayersWithMeasurement_.clear();
        jet_daughter_trackerLayersWithMeasurement_.clear();
        jet_daughter_trackHighPurity_.clear();

        jet_daughter_caloFraction_.clear();
        jet_daughter_hcalFraction_.clear();

        jet_daughter_puppiWeight_.clear();
        jet_daughter_puppiWeightNoLep_.clear();

        jet_daughter_isIsolatedChargedHadron_.clear();
        jet_daughter_isMuon_.clear();
        jet_daughter_isElectron_.clear();
        jet_daughter_isPhoton_.clear();
        jet_daughter_isStandAloneMuon_.clear();
        jet_daughter_isTrackerMuon_.clear();
        jet_daughter_isGlobalMuon_.clear();
        jet_daughter_isGoodEgamma_.clear();
        jet_daughter_isConvertedPhoton_.clear();

        for (size_t d=0; d<j.numberOfDaughters(); d++) {
            const pat::PackedCandidate * dau = (pat::PackedCandidate*)j.daughter(d);
            jet_daughter_pt_.push_back(dau->pt());
            jet_daughter_eta_.push_back(dau->eta());
            jet_daughter_phi_.push_back(dau->phi());
            jet_daughter_energy_.push_back(dau->energy());
            jet_daughter_mass_.push_back(dau->mass());
            jet_daughter_charge_.push_back(dau->charge());

            if (dau->hasTrackDetails()) {
                jet_daughter_etaAtVtx_.push_back(dau->etaAtVtx());
                jet_daughter_phiAtVtx_.push_back(dau->phiAtVtx());
                jet_daughter_dxy_.push_back(dau->dxy());
                jet_daughter_dxyError_.push_back(dau->dxyError());
                jet_daughter_dz_.push_back(dau->dz());
                jet_daughter_dzError_.push_back(dau->dzError());
                jet_daughter_vertexChi2_.push_back(dau->vertexChi2());
                jet_daughter_vertexNdof_.push_back(dau->vertexNdof());
                jet_daughter_vx_.push_back(dau->vx());
                jet_daughter_vy_.push_back(dau->vy());
                jet_daughter_vz_.push_back(dau->vz());

                jet_daughter_pixelLayersWithMeasurement_.push_back(dau->pixelLayersWithMeasurement());
                jet_daughter_stripLayersWithMeasurement_.push_back(dau->stripLayersWithMeasurement());
                jet_daughter_trackerLayersWithMeasurement_.push_back(dau->trackerLayersWithMeasurement());
                jet_daughter_trackHighPurity_.push_back(dau->trackHighPurity());
            } else {
                jet_daughter_etaAtVtx_.push_back(0);
                jet_daughter_phiAtVtx_.push_back(0);
                jet_daughter_dxy_.push_back(0);
                jet_daughter_dxyError_.push_back(0);
                jet_daughter_dz_.push_back(0);
                jet_daughter_dzError_.push_back(0);
                jet_daughter_vertexChi2_.push_back(0);
                jet_daughter_vertexNdof_.push_back(0);
                jet_daughter_vx_.push_back(0);
                jet_daughter_vy_.push_back(0);
                jet_daughter_vz_.push_back(0);

                jet_daughter_pixelLayersWithMeasurement_.push_back(0);
                jet_daughter_stripLayersWithMeasurement_.push_back(0);
                jet_daughter_trackerLayersWithMeasurement_.push_back(0);
                jet_daughter_trackHighPurity_.push_back(0);
            }

            jet_daughter_dtime_.push_back(dau->dtime());
            jet_daughter_time_.push_back(dau->time());
            jet_daughter_timeError_.push_back(dau->timeError());

            jet_daughter_pdgId_.push_back(dau->pdgId());

            jet_daughter_caloFraction_.push_back(dau->caloFraction());
            jet_daughter_hcalFraction_.push_back(dau->hcalFraction());

            jet_daughter_puppiWeight_.push_back(dau->puppiWeight());
            jet_daughter_puppiWeightNoLep_.push_back(dau->puppiWeightNoLep());

            jet_daughter_isIsolatedChargedHadron_.push_back(dau->isIsolatedChargedHadron());
            jet_daughter_isMuon_.push_back(dau->isMuon());
            jet_daughter_isElectron_.push_back(dau->isElectron());
            jet_daughter_isPhoton_.push_back(dau->isPhoton());
            jet_daughter_isStandAloneMuon_.push_back(dau->isStandAloneMuon());
            jet_daughter_isTrackerMuon_.push_back(dau->isTrackerMuon());
            jet_daughter_isGlobalMuon_.push_back(dau->isGlobalMuon());
            jet_daughter_isGoodEgamma_.push_back(dau->isGoodEgamma());
            jet_daughter_isConvertedPhoton_.push_back(dau->isConvertedPhoton());
        }

        if (j.pt()>20.0) tree_->Fill();

        idx++;
    }

}
