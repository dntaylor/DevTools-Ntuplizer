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
    tree_->Branch("num_vertices", &num_vertices_);

    // jet global branches
    tree_->Branch("jet_pt",     &jet_pt_);
    tree_->Branch("jet_eta",    &jet_eta_);
    tree_->Branch("jet_phi",    &jet_phi_);
    tree_->Branch("jet_energy", &jet_energy_);
    tree_->Branch("jet_mass",   &jet_mass_);

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

        
        if (j.pt()>20.0) tree_->Fill();

        idx++;
    }

}
