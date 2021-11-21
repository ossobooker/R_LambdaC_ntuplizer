///////////////////////// Code to produce K* candidates ////////////////////////

#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include <vector>
#include <memory>
#include <map>
#include <string>
#include <cmath>
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "helper.h"
#include <limits>
#include <algorithm>
#include "KinVtxFitter.h"

class LambdaCBuilder : public edm::global::EDProducer<>
{

public:
  typedef std::vector<reco::TransientTrack> TransientTrackCollection;

  explicit LambdaCBuilder(const edm::ParameterSet &cfg) : trk1_selection_{cfg.getParameter<std::string>("trk1Selection")},
                                                          trk2_selection_{cfg.getParameter<std::string>("trk2Selection")},
                                                          trk3_selection_{cfg.getParameter<std::string>("trk3Selection")},
                                                          pre_vtx_selection_{cfg.getParameter<std::string>("preVtxSelection")},
                                                          post_vtx_selection_{cfg.getParameter<std::string>("postVtxSelection")},
                                                          pfcands_{consumes<pat::CompositeCandidateCollection>(cfg.getParameter<edm::InputTag>("pfcands"))},
                                                          ttracks_{consumes<TransientTrackCollection>(cfg.getParameter<edm::InputTag>("transientTracks"))}
  {

    //output
    produces<pat::CompositeCandidateCollection>();
  }

  ~LambdaCBuilder() override {}

  void produce(edm::StreamID, edm::Event &, const edm::EventSetup &) const override;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions) {}

private:
  const StringCutObjectSelector<pat::CompositeCandidate> trk1_selection_;     // cuts on leading cand
  const StringCutObjectSelector<pat::CompositeCandidate> trk2_selection_;     // sub-leading cand
  const StringCutObjectSelector<pat::CompositeCandidate> trk3_selection_;     // sub-leading cand
  const StringCutObjectSelector<pat::CompositeCandidate> pre_vtx_selection_;  // cut on the di-lepton before the SV fit
  const StringCutObjectSelector<pat::CompositeCandidate> post_vtx_selection_; // cut on the di-lepton after the SV fit
  const edm::EDGetTokenT<pat::CompositeCandidateCollection> pfcands_;         //input PF cands this is sorted in pT in previous step
  const edm::EDGetTokenT<TransientTrackCollection> ttracks_;                  //input TTracks of PF cands
};

void LambdaCBuilder::produce(edm::StreamID, edm::Event &evt, edm::EventSetup const &) const
{

  //inputs
  edm::Handle<pat::CompositeCandidateCollection> pfcands;
  evt.getByToken(pfcands_, pfcands);

  edm::Handle<TransientTrackCollection> ttracks;
  evt.getByToken(ttracks_, ttracks);

  // output
  std::unique_ptr<pat::CompositeCandidateCollection> lambda_c_out(new pat::CompositeCandidateCollection());

  // main loop
  for (size_t trk1_idx = 0; trk1_idx < pfcands->size(); ++trk1_idx)
  {

    edm::Ptr<pat::CompositeCandidate> trk1_ptr(pfcands, trk1_idx);
    if (!trk1_selection_(*trk1_ptr))
      continue;

    for (size_t trk2_idx = trk1_idx + 1; trk2_idx < pfcands->size(); ++trk2_idx)
    {
      edm::Ptr<pat::CompositeCandidate> trk2_ptr(pfcands, trk2_idx);
      if (!trk2_selection_(*trk2_ptr))
        continue;

      for (size_t trk3_idx = trk2_idx + 1; trk3_idx < pfcands->size(); ++trk3_idx)
      {
        edm::Ptr<pat::CompositeCandidate> trk3_ptr(pfcands, trk3_idx);
        if (!trk3_selection_(*trk3_ptr))
          continue;

        if (std::abs(trk1_ptr->charge() + trk2_ptr->charge() + trk3_ptr->charge()) != 1)
          continue;

        // create a Lambda_c candidate; add first quantities that can be used for pre fit selection
        pat::CompositeCandidate lambda_c_cand;
        auto trk1_p4 = trk1_ptr->polarP4();
        auto trk2_p4 = trk2_ptr->polarP4();
        auto trk3_p4 = trk3_ptr->polarP4();
        trk1_p4.SetM(PROTON_MASS);
        trk2_p4.SetM(K_MASS);
        trk3_p4.SetM(PI_MASS);

        //adding stuff for pre fit selection
        lambda_c_cand.setP4(trk1_p4 + trk2_p4 + trk3_p4);
        lambda_c_cand.setCharge(trk1_ptr->charge() + trk2_ptr->charge() + trk3_ptr->charge());
        lambda_c_cand.addUserFloat("trk_deltaR12", reco::deltaR(*trk1_ptr, *trk2_ptr));
        lambda_c_cand.addUserFloat("trk_deltaR13", reco::deltaR(*trk1_ptr, *trk3_ptr));
        lambda_c_cand.addUserFloat("trk_deltaR23", reco::deltaR(*trk2_ptr, *trk3_ptr));

        // save indices
        lambda_c_cand.addUserInt("trk1_idx", trk1_idx);
        lambda_c_cand.addUserInt("trk2_idx", trk2_idx);
        lambda_c_cand.addUserInt("trk3_idx", trk3_idx);

        // save cands
        lambda_c_cand.addUserCand("trk1", trk1_ptr);
        lambda_c_cand.addUserCand("trk2", trk2_ptr);
        lambda_c_cand.addUserCand("trk3", trk3_ptr);

        //second mass hypothesis
        trk1_p4.SetM(PROTON_MASS);
        trk2_p4.SetM(PI_MASS);
        trk3_p4.SetM(K_MASS);
        lambda_c_cand.addUserFloat("1barMass", (trk1_p4 + trk2_p4 + trk3_p4).M());

        //THIRD mass hypothesis
        trk1_p4.SetM(PI_MASS);
        trk2_p4.SetM(K_MASS);
        trk3_p4.SetM(PROTON_MASS);
        lambda_c_cand.addUserFloat("2barMass", (trk1_p4 + trk2_p4 + trk3_p4).M());

        //FORTH mass hypothesis
        trk1_p4.SetM(K_MASS);
        trk2_p4.SetM(PROTON_MASS);
        trk3_p4.SetM(PI_MASS);
        lambda_c_cand.addUserFloat("3barMass", (trk1_p4 + trk2_p4 + trk3_p4).M());

        // selection before fit
        if (!pre_vtx_selection_(lambda_c_cand))
          continue;

        KinVtxFitter fitter(
            {ttracks->at(trk1_idx), ttracks->at(trk2_idx), ttracks->at(trk3_idx)},
            {PROTON_MASS, K_MASS, PI_MASS},
            {PROTON_SIGMA, K_SIGMA, K_SIGMA} //PROTON, K and PI sigma equals...
        );
        if (!fitter.success())
          continue;

        // save quantities after fit
        lambda_c_cand.addUserFloat("sv_chi2", fitter.chi2());
        lambda_c_cand.addUserFloat("sv_ndof", fitter.dof());
        lambda_c_cand.addUserFloat("sv_prob", fitter.prob());
        lambda_c_cand.addUserFloat("fitted_mass", fitter.fitted_candidate().mass());
        lambda_c_cand.addUserFloat("fitted_pt", fitter.fitted_candidate().globalMomentum().perp());
        lambda_c_cand.addUserFloat("fitted_eta", fitter.fitted_candidate().globalMomentum().eta());
        lambda_c_cand.addUserFloat("fitted_phi", fitter.fitted_candidate().globalMomentum().phi());

        // second mass hypothesis
        auto fitted_trk1 = fitter.daughter_p4(0);
        auto fitted_trk2 = fitter.daughter_p4(1);
        auto fitted_trk3 = fitter.daughter_p4(2);
        fitted_trk1.SetM(PROTON_MASS);
        fitted_trk2.SetM(PI_MASS);
        fitted_trk3.SetM(K_MASS);
        lambda_c_cand.addUserFloat("fitted_1barMass", (fitted_trk1 + fitted_trk2 + fitted_trk3).M());

        // THIRD mass hypothesis
        fitted_trk1 = fitter.daughter_p4(0);
        fitted_trk2 = fitter.daughter_p4(1);
        fitted_trk3 = fitter.daughter_p4(2);
        fitted_trk1.SetM(PI_MASS);
        fitted_trk2.SetM(K_MASS);
        fitted_trk3.SetM(PROTON_MASS);
        lambda_c_cand.addUserFloat("fitted_2barMass", (fitted_trk1 + fitted_trk2 + fitted_trk3).M());

        // FOURTH mass hypothesis
        fitted_trk1 = fitter.daughter_p4(0);
        fitted_trk2 = fitter.daughter_p4(1);
        fitted_trk3 = fitter.daughter_p4(2);
        fitted_trk1.SetM(K_MASS);
        fitted_trk2.SetM(PROTON_MASS);
        fitted_trk3.SetM(PI_MASS);
        lambda_c_cand.addUserFloat("fitted_3barMass", (fitted_trk1 + fitted_trk2 + fitted_trk3).M());

        // after fit selection
        if (!post_vtx_selection_(lambda_c_cand))
          continue;
        lambda_c_out->emplace_back(lambda_c_cand);
      }
    }
  }

  evt.put(std::move(lambda_c_out));
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(LambdaCBuilder);
