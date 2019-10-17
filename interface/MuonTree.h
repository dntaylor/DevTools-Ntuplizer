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
    edm::EDGetTokenT<edm::ValueMap<reco::MuonSimInfo>> simInfo_;

    // branches
    // kinematics
    std::vector<float> muon_pt_;
    std::vector<float> muon_p_;
    std::vector<float> muon_eta_;
    std::vector<float> muon_phi_;
    std::vector<float> muon_energy_;
    std::vector<float> muon_mass_;

    // id
    std::vector<int> muon_isPFMuon_;
    std::vector<int> muon_isTrackerMuon_;
    std::vector<int> muon_isGlobalMuon_;
    std::vector<int> muon_isStandAloneMuon_;
    std::vector<int> muon_isCaloMuon_;
    std::vector<int> muon_CutBasedIdLoose_;
    std::vector<int> muon_CutBasedIdMedium_;
    std::vector<int> muon_CutBasedIdMediumPrompt_;
    std::vector<int> muon_CutBasedIdTight_;
    std::vector<int> muon_PFIsoLoose_;
    std::vector<int> muon_PFIsoTight_;
    std::vector<int> muon_SoftCutBasedId_;

    // inner track
    std::vector<float> muon_innerTrack_pt_;
    std::vector<float> muon_innerTrack_p_;
    std::vector<float> muon_innerTrack_eta_;
    std::vector<float> muon_innerTrack_phi_;
    std::vector<float> muon_innerTrack_qoverp_;
    std::vector<float> muon_innerTrack_qoverpError_;
    std::vector<float> muon_innerTrack_validFraction_;
    std::vector<int> muon_innerTrack_highPurity_;
    std::vector<int> muon_innerTrack_hitPattern_trackerLayersWithMeasurement_;
    std::vector<int> muon_innerTrack_hitPattern_pixelLayersWithMeasurement_;

    // calo compatibility
    std::vector<float> muon_caloCompatibility_;
    std::vector<float> muon_calEnergy_ecal_time_;
    std::vector<float> muon_calEnergy_em_;
    std::vector<float> muon_calEnergy_emMax_;
    std::vector<float> muon_calEnergy_emS25_;
    std::vector<float> muon_calEnergy_emS9_;
    std::vector<float> muon_calEnergy_had_;
    std::vector<float> muon_calEnergy_hadMax_;
    std::vector<float> muon_calEnergy_hadS9_;
    std::vector<float> muon_calEnergy_hcal_time_;
    std::vector<float> muon_calEnergy_ho_;
    std::vector<float> muon_calEnergy_hoS9_;
    std::vector<float> muon_calEnergy_tower_;
    std::vector<float> muon_calEnergy_towerS9_;
    std::vector<std::vector<int> > muon_calEnergy_crossedHadRecHits_ieta_;
    std::vector<std::vector<int> > muon_calEnergy_crossedHadRecHits_iphi_;
    std::vector<std::vector<int> > muon_calEnergy_crossedHadRecHits_depth_;
    std::vector<std::vector<float> > muon_calEnergy_crossedHadRecHits_energy_;
    std::vector<std::vector<float> > muon_calEnergy_crossedHadRecHits_time_;
    std::vector<std::vector<float> > muon_calEnergy_crossedHadRecHits_chi2_;

    // gen match
    std::vector<int> muon_gen_matches_muon_;
    std::vector<int> muon_gen_matches_pion_;
    std::vector<float> muon_gen_deltaR_;
    std::vector<float> muon_gen_pt_;
    std::vector<float> muon_gen_eta_;
    std::vector<float> muon_gen_phi_;
    std::vector<float> muon_gen_mass_;

    // gen match via sim hit
    std::vector<int> muon_gen_sim_primaryClass_;
    std::vector<int> muon_gen_sim_extendedClass_;
    std::vector<int> muon_gen_sim_flavour_;
    std::vector<int> muon_gen_sim_pdgId_;
    std::vector<float> muon_gen_sim_pt_;
    std::vector<float> muon_gen_sim_eta_;
    std::vector<float> muon_gen_sim_phi_;
    std::vector<float> muon_gen_sim_mass_;
    std::vector<float> muon_gen_sim_tpAssoQuality_;
};

void MuonTree::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


DEFINE_FWK_MODULE(MuonTree);

#endif
