#ifndef DeepJetTree_h
#define DeepJetTree_h

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TTree.h"

#include "FWCore/Utilities/interface/InputTag.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

class DeepJetTree : public edm::one::EDAnalyzer<edm::one::SharedResources> {
  public:
    explicit DeepJetTree(const edm::ParameterSet&);
    ~DeepJetTree();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  private:
    virtual void beginJob() override;
    virtual void analyze(edm::Event const& iEvent, edm::EventSetup const&) override;
    virtual void endJob() override;

    // trees
    TTree *tree_;

    // tokens
    edm::EDGetTokenT<edm::View<pat::Jet>> jetsToken_;
    edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;
    edm::EDGetTokenT<reco::VertexCollection> verticesToken_;

    // configuration
    edm::ParameterSet muonAssociation_;

    // branches
    
    // event
    int num_vertices_;

    // kinematics
    float jet_pt_;
    float jet_eta_;
    float jet_phi_;
    float jet_energy_;
    float jet_mass_;
    float jet_jetCharge_;

    // energy sums/fractions/multiplicity
    float jet_numberOfDaughters_;
    float jet_chargedMultiplicity_;
    float jet_neutralMultiplicity_;

    float jet_chargedHadronMultiplicity_;
    float jet_neutralHadronMultiplicity_;
    float jet_muonMultiplicity_;
    float jet_electronMultiplicity_;
    float jet_photonMultiplicity_;

    float jet_chargedEmEnergy_;
    float jet_neutralEmEnergy_;
    float jet_chargedHadronEnergy_;
    float jet_neutralHadronEnergy_;
    float jet_muonEnergy_;
    float jet_electronEnergy_;
    float jet_photonEnergy_;

    float jet_chargedEmEnergyFraction_;
    float jet_neutralEmEnergyFraction_;
    float jet_chargedHadronEnergyFraction_;
    float jet_neutralHadronEnergyFraction_;
    float jet_muonEnergyFraction_;
    float jet_electronEnergyFraction_;
    float jet_photonEnergyFraction_;

    float jet_n60_;
    float jet_n90_;

    // gen truth
    float jet_hadronFlavour_;
    float jet_partonFlavour_;

    // btag scores AK4
    float jet_pfJetBProbabilityBJetTags_;
    float jet_pfJetProbabilityBJetTags_;
    float jet_pfTrackCountingHighEffBJetTags_;
    float jet_pfSimpleSecondaryVertexHighEffBJetTags_;
    float jet_pfSimpleInclusiveSecondaryVertexHighEffBJetTags_;
    float jet_pfCombinedSecondaryVertexV2BJetTags_;
    float jet_pfCombinedInclusiveSecondaryVertexV2BJetTags_;
    float jet_softPFMuonBJetTags_;
    float jet_softPFElectronBJetTags_;
    float jet_pfCombinedMVAV2BJetTags_;
    float jet_pfCombinedCvsLJetTags_;
    float jet_pfCombinedCvsBJetTags_;
    float jet_pfDeepCSVJetTags_probb_;
    float jet_pfDeepCSVJetTags_probc_;
    float jet_pfDeepCSVJetTags_probudsg_;
    float jet_pfDeepCSVJetTags_probbb_;

    // jet daughters
    std::vector<float> jet_daughter_pt_;
    std::vector<float> jet_daughter_eta_;
    std::vector<float> jet_daughter_phi_;
    std::vector<float> jet_daughter_energy_;
    std::vector<float> jet_daughter_mass_;
    std::vector<float> jet_daughter_charge_;

    std::vector<float> jet_daughter_etaAtVtx_;
    std::vector<float> jet_daughter_phiAtVtx_;
    std::vector<float> jet_daughter_vertexChi2_;
    std::vector<float> jet_daughter_vertexNdof_;
    std::vector<float> jet_daughter_vx_;
    std::vector<float> jet_daughter_vy_;
    std::vector<float> jet_daughter_vz_;
    std::vector<float> jet_daughter_dxy_;
    std::vector<float> jet_daughter_dxyError_;
    std::vector<float> jet_daughter_dz_;
    std::vector<float> jet_daughter_dzError_;
    std::vector<float> jet_daughter_dtime_;
    std::vector<float> jet_daughter_time_;
    std::vector<float> jet_daughter_timeError_;

    std::vector<float> jet_daughter_pdgId_;

    std::vector<float> jet_daughter_pixelLayersWithMeasurement_;
    std::vector<float> jet_daughter_stripLayersWithMeasurement_;
    std::vector<float> jet_daughter_trackerLayersWithMeasurement_;
    std::vector<float> jet_daughter_trackHighPurity_;

    std::vector<float> jet_daughter_caloFraction_;
    std::vector<float> jet_daughter_hcalFraction_;

    std::vector<float> jet_daughter_puppiWeight_;
    std::vector<float> jet_daughter_puppiWeightNoLep_;

    std::vector<float> jet_daughter_isIsolatedChargedHadron_;
    std::vector<float> jet_daughter_isMuon_;
    std::vector<float> jet_daughter_isElectron_;
    std::vector<float> jet_daughter_isPhoton_;
    std::vector<float> jet_daughter_isStandAloneMuon_;
    std::vector<float> jet_daughter_isTrackerMuon_;
    std::vector<float> jet_daughter_isGlobalMuon_;
    std::vector<float> jet_daughter_isGoodEgamma_;
    std::vector<float> jet_daughter_isConvertedPhoton_;

    // TODO: gen part truth

};

void DeepJetTree::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


DEFINE_FWK_MODULE(DeepJetTree);

#endif
