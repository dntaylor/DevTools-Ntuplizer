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


};

void DeepJetTree::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


DEFINE_FWK_MODULE(DeepJetTree);

#endif
