#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include <vector>
#include <memory>
#include <map>
#include <string>
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "helper.h"
#include <limits>
#include <algorithm>
#include "KinVtxFitter.h"

template <typename Lepton>
class SingleLeptonBuilder : public edm::global::EDProducer<>
{

  // perhaps we need better structure here (begin run etc)
public:
  typedef std::vector<Lepton> LeptonCollection;
  typedef std::vector<reco::TransientTrack> TransientTrackCollection;

  explicit SingleLeptonBuilder(const edm::ParameterSet &cfg) : l1_selection_{cfg.getParameter<std::string>("lep1Selection")},
                                                               l2_selection_{cfg.getParameter<std::string>("lep2Selection")},
                                                               pre_vtx_selection_{cfg.getParameter<std::string>("preVtxSelection")},
                                                               post_vtx_selection_{cfg.getParameter<std::string>("postVtxSelection")},
                                                               src_{consumes<LeptonCollection>(cfg.getParameter<edm::InputTag>("src"))},
                                                               ttracks_src_{consumes<TransientTrackCollection>(cfg.getParameter<edm::InputTag>("transientTracksSrc"))}
  {
    produces<pat::CompositeCandidateCollection>("SelectedSingleLeptons");
    produces<std::vector<KinVtxFitter>>("SelectedSingleLeptonKinVtxs");
  }

  ~SingleLeptonBuilder() override {}

  void produce(edm::StreamID, edm::Event &, const edm::EventSetup &) const override;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions) {}

private:
  const StringCutObjectSelector<Lepton> l1_selection_;                        // cut on leading lepton
  const StringCutObjectSelector<Lepton> l2_selection_;                        // cut on sub-leading lepton
  const StringCutObjectSelector<pat::CompositeCandidate> pre_vtx_selection_;  // cut on the di-lepton before the SV fit
  const StringCutObjectSelector<pat::CompositeCandidate> post_vtx_selection_; // cut on the di-lepton after the SV fit
  const edm::EDGetTokenT<LeptonCollection> src_;
  const edm::EDGetTokenT<TransientTrackCollection> ttracks_src_;
};

template <typename Lepton>
void SingleLeptonBuilder<Lepton>::produce(edm::StreamID, edm::Event &evt, edm::EventSetup const &) const
{

  //input
  edm::Handle<LeptonCollection> leptons;
  evt.getByToken(src_, leptons);

  edm::Handle<TransientTrackCollection> ttracks;
  evt.getByToken(ttracks_src_, ttracks);

  // output
  std::unique_ptr<pat::CompositeCandidateCollection> ret_value(new pat::CompositeCandidateCollection());
  std::unique_ptr<std::vector<KinVtxFitter>> kinVtx_out(new std::vector<KinVtxFitter>);

  for (size_t l1_idx = 0; l1_idx < leptons->size(); ++l1_idx)
  {
    edm::Ptr<Lepton> l1_ptr(leptons, l1_idx);
    if (!l1_selection_(*l1_ptr))
      continue;

    for (size_t l2_idx = l1_idx + 1; l2_idx < leptons->size(); ++l2_idx)
    {
      edm::Ptr<Lepton> l2_ptr(leptons, l2_idx);
      if (!l2_selection_(*l2_ptr))
        continue;

      pat::CompositeCandidate lepton_pair;
      lepton_pair.setP4(l1_ptr->p4() + l2_ptr->p4());
      lepton_pair.setCharge(l1_ptr->charge() + l2_ptr->charge());
      lepton_pair.addUserFloat("lep_deltaR", reco::deltaR(*l1_ptr, *l2_ptr));
      int nlowpt = 0;
      if (l1_ptr->hasUserInt("isPF") && l2_ptr->hasUserInt("isPF"))
        nlowpt = 2 - l1_ptr->userInt("isPF") - l2_ptr->userInt("isPF");

      // Put the lepton passing the corresponding selection
      lepton_pair.addUserInt("l1_idx", l1_idx);
      lepton_pair.addUserInt("l2_idx", l2_idx);
      // Use UserCands as they should not use memory but keep the Ptr itself
      lepton_pair.addUserCand("l1", l1_ptr);
      lepton_pair.addUserCand("l2", l2_ptr);
      lepton_pair.addUserInt("nlowpt", nlowpt);
      if (!pre_vtx_selection_(lepton_pair))
        continue; // before making the SV, cut on the info we have

      KinVtxFitter fitter(
          {ttracks->at(l1_idx), ttracks->at(l2_idx)},
          {l1_ptr->mass(), l2_ptr->mass()},
          {LEP_SIGMA, LEP_SIGMA} //some small sigma for the particle mass
      );
      lepton_pair.addUserFloat("sv_chi2", fitter.chi2());
      lepton_pair.addUserFloat("sv_ndof", fitter.dof()); // float??
      lepton_pair.addUserFloat("sv_prob", fitter.prob());
      lepton_pair.addUserFloat("fitted_mass", fitter.success() ? fitter.fitted_candidate().mass() : -1);
      lepton_pair.addUserFloat("fitted_massErr", fitter.success() ? sqrt(fitter.fitted_candidate().kinematicParametersError().matrix()(6, 6)) : -1);
      // if needed, add here more stuff

      // cut on the SV info
      if (!post_vtx_selection_(lepton_pair))
        continue;
      ret_value->push_back(lepton_pair);
      kinVtx_out->push_back(fitter);
    }
  }

  evt.put(std::move(ret_value), "SelectedSingleLeptons");
  evt.put(std::move(kinVtx_out), "SelectedSingleLeptonKinVtxs");
}

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
typedef SingleLeptonBuilder<pat::Muon> SingleMuonBuilder;
// typedef SingleLeptonBuilder<pat::Electron> DiElectronBuilder;

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(SingleMuonBuilder);
// DEFINE_FWK_MODULE(SingleElectronBuilder);
