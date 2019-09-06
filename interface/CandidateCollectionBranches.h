#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/Event.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "CommonTools/Utils/interface/StringObjectFunction.h"
#include "DataFormats/Candidate/interface/Candidate.h"

#include "TTree.h"

template<typename T>
class CandidateCollectionFunction {
  public:
    CandidateCollectionFunction(TTree * tree, std::string functionName, std::string functionString, int maxCount);
    void evaluate(const reco::CandidateView& candidates);

  private:
    StringObjectFunction<reco::Candidate, true> function_;
    std::string functionString_;
    std::string functionName_;
    TBranch * vectorBranch_;
    std::vector<T> values_; 
    int maxCount_;
};

typedef CandidateCollectionFunction<int> CandidateCollectionIntFunction;
typedef CandidateCollectionFunction<float> CandidateCollectionFloatFunction;

template<typename T>
class JetConstituentCollectionFunction : CandidateCollectionFunction<T> {
  public:
    JetConstituentCollectionFunction(TTree * tree, std::string functionName, std::string functionString, int maxCount, int constituentMaxCount);
    void evaluate(const reco::CandidateView& candidates);

  private:
    std::vector<std::vector<T> > values_; 
    int constituentMaxCount_;
};

typedef JetConstituentCollectionFunction<int> JetConstituentCollectionIntFunction;
typedef JetConstituentCollectionFunction<float> JetConstituentCollectionFloatFunction;

class CandidateCollectionBranches {
  public:
    CandidateCollectionBranches(TTree * tree, std::string collectionName,  const edm::ParameterSet& iConfig, edm::ConsumesCollector cc);
    void fill(const edm::Event& iEvent);
    std::string getName() { return collectionName_; }
    int getCount() { return collectionCount_; }
    int keep() { return minCount_ > 0 ? collectionCount_ >= minCount_ : false; }

  private:
    edm::EDGetTokenT<reco::CandidateView> collectionToken_;
    edm::ParameterSet branches_;
    std::vector<std::unique_ptr<CandidateCollectionFloatFunction> > floatFunctions_;
    std::vector<std::unique_ptr<CandidateCollectionIntFunction> > intFunctions_;
    TBranch * collectionCountBranch_;
    std::string collectionName_;
    int collectionCount_;
    int minCount_;
    int maxCount_;
};

class JetCandidateCollectionBranches : CandidateCollectionBranches {
  public:
    JetCandidateCollectionBranches(TTree * tree, std::string collectionName, std::string constituentCollectionName,  const edm::ParameterSet& iConfig, edm::ConsumesCollector cc);
    std::string getConstituentName() { return constituentCollectionName_; }

  private:
    edm::ParameterSet constituentBranches_;
    TBranch * constituentCollectionCountBranch_;
    std::string constituentCollectionName_;
    std::vector<int> constituentCollectionCount_;
    std::vector<std::unique_ptr<JetConstituentCollectionFloatFunction> > constituentFloatFunctions_;
    std::vector<std::unique_ptr<JetConstituentCollectionIntFunction> > constituentIntFunctions_;
    int constituentMinCount_;
    int constituentMaxCount_;
};
