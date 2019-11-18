#ifndef MuonTree_h
#define MuonTree_h

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

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/MuonReco/interface/MuonSimInfo.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"


class MuonTree : public edm::one::EDAnalyzer<edm::one::SharedResources> {
  public:
    explicit MuonTree(const edm::ParameterSet&);
    ~MuonTree();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  private:
    virtual void beginJob() override;
    virtual void analyze(edm::Event const& iEvent, edm::EventSetup const&) override;
    virtual void endJob() override;

    // trees
    TTree *tree_;

    // tokens
    edm::EDGetTokenT<edm::View<reco::Muon>> muonsToken_;
    edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;
    edm::EDGetTokenT<reco::VertexCollection> verticesToken_;
    edm::EDGetTokenT<edm::ValueMap<reco::MuonSimInfo>> simInfo_;

    // branches
    
    // event
    int num_vertices_;

    // kinematics
    float muon_pt_;
    float muon_p_;
    float muon_eta_;
    float muon_phi_;
    float muon_energy_;
    float muon_mass_;
    int muon_charge_;

    // id
    int muon_isPFMuon_;
    int muon_isTrackerMuon_;
    int muon_isGlobalMuon_;
    int muon_isStandAloneMuon_;
    int muon_isCaloMuon_;
    int muon_CutBasedIdLoose_;
    int muon_CutBasedIdMedium_;
    int muon_CutBasedIdMediumPrompt_;
    int muon_CutBasedIdTight_;
    int muon_PFIsoLoose_;
    int muon_PFIsoTight_;
    int muon_SoftCutBasedId_;

    // isolation
    int muon_isolationR03_nTracks_;
    float muon_isolationR03_sumPt_;

    // inner track
    float muon_innerTrack_pt_;
    float muon_innerTrack_p_;
    float muon_innerTrack_eta_;
    float muon_innerTrack_phi_;
    float muon_innerTrack_qoverp_;
    float muon_innerTrack_qoverpError_;
    float muon_innerTrack_validFraction_;
    int muon_innerTrack_highPurity_;
    int muon_innerTrack_hitPattern_trackerLayersWithMeasurement_;
    int muon_innerTrack_hitPattern_pixelLayersWithMeasurement_;

    // calo compatibility
    float muon_caloCompatibility_;
    float muon_calEnergy_ecal_time_;
    float muon_calEnergy_em_;
    float muon_calEnergy_emMax_;
    float muon_calEnergy_emS25_;
    float muon_calEnergy_emS9_;
    float muon_calEnergy_had_;
    float muon_calEnergy_hadMax_;
    float muon_calEnergy_hadS9_;
    float muon_calEnergy_hcal_time_;
    float muon_calEnergy_ho_;
    float muon_calEnergy_hoS9_;
    float muon_calEnergy_tower_;
    float muon_calEnergy_towerS9_;
    int muon_calEnergy_ecal_ieta_;
    int muon_calEnergy_ecal_iphi_;
    int muon_calEnergy_hcal_ieta_;
    int muon_calEnergy_hcal_iphi_;
    int muon_calEnergy_hcal_depth_;
    std::vector<int> muon_calEnergy_crossedHadRecHits_ieta_;
    std::vector<int> muon_calEnergy_crossedHadRecHits_iphi_;
    std::vector<int> muon_calEnergy_crossedHadRecHits_depth_;
    std::vector<float> muon_calEnergy_crossedHadRecHits_energy_;
    std::vector<float> muon_calEnergy_crossedHadRecHits_time_;
    std::vector<float> muon_calEnergy_crossedHadRecHits_chi2_;

    // gen match
    int muon_gen_matches_muon_;
    int muon_gen_matches_pion_;
    float muon_gen_deltaR_;
    float muon_gen_pt_;
    float muon_gen_eta_;
    float muon_gen_phi_;
    float muon_gen_mass_;

    // gen match via sim hit
    int muon_gen_sim_primaryClass_;
    int muon_gen_sim_extendedClass_;
    int muon_gen_sim_flavour_;
    int muon_gen_sim_pdgId_;
    float muon_gen_sim_pt_;
    float muon_gen_sim_eta_;
    float muon_gen_sim_phi_;
    float muon_gen_sim_mass_;
    float muon_gen_sim_tpAssoQuality_;
};

void MuonTree::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


DEFINE_FWK_MODULE(MuonTree);

#endif
