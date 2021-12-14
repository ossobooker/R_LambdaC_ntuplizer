////////////////////// Code to produce K*LL candidates /////////////////////////

#include <vector>
#include <memory>
#include <map>
#include <string>
#include <limits>
#include <algorithm>

#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

#include "helper.h"
#include "KinVtxFitter.h"

class LambdaBToLambdaCMuBuilder : public edm::global::EDProducer<>
{

  // perhaps we need better structure here (begin run etc)
public:
  typedef std::vector<reco::TransientTrack> TransientTrackCollection;

  explicit LambdaBToLambdaCMuBuilder(const edm::ParameterSet &cfg) : // selections
                                                                     muonSelection_{cfg.getParameter<std::string>("muonSelection")},
                                                                     pre_vtx_selection_{cfg.getParameter<std::string>("preVtxSelection")},
                                                                     post_vtx_selection_{cfg.getParameter<std::string>("postVtxSelection")},

                                                                     //inputs
                                                                     muons_{consumes<pat::MuonCollection>(cfg.getParameter<edm::InputTag>("muons"))},
                                                                     lambdas_c_{consumes<pat::CompositeCandidateCollection>(cfg.getParameter<edm::InputTag>("lambdas_c"))},
                                                                     muon_ttracks_{consumes<TransientTrackCollection>(cfg.getParameter<edm::InputTag>("muonTransientTracks"))},
                                                                     lambdas_c_ttracks_{consumes<TransientTrackCollection>(cfg.getParameter<edm::InputTag>("lambdas_cTransientTracks"))},
                                                                     isotracksToken_{consumes<pat::PackedCandidateCollection>(cfg.getParameter<edm::InputTag>("tracks"))},
                                                                     isolostTracksToken_{consumes<pat::PackedCandidateCollection>(cfg.getParameter<edm::InputTag>("lostTracks"))},
                                                                     isotrk_selection_{cfg.getParameter<std::string>("isoTracksSelection")},
                                                                     beamspot_{consumes<reco::BeamSpot>(cfg.getParameter<edm::InputTag>("beamSpot"))}
  {
    //output
    produces<pat::CompositeCandidateCollection>();
  }

  ~LambdaBToLambdaCMuBuilder() override {}

  void produce(edm::StreamID, edm::Event &, const edm::EventSetup &) const override;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions) {}

private:
  // selections
  const StringCutObjectSelector<pat::Muon> muonSelection_;
  const StringCutObjectSelector<pat::CompositeCandidate> pre_vtx_selection_;
  const StringCutObjectSelector<pat::CompositeCandidate> post_vtx_selection_;

  // inputs
  const edm::EDGetTokenT<pat::MuonCollection> muons_;
  const edm::EDGetTokenT<pat::CompositeCandidateCollection> lambdas_c_;
  const edm::EDGetTokenT<TransientTrackCollection> muon_ttracks_;
  const edm::EDGetTokenT<TransientTrackCollection> lambdas_c_ttracks_;
  const edm::EDGetTokenT<pat::PackedCandidateCollection> isotracksToken_;
  const edm::EDGetTokenT<pat::PackedCandidateCollection> isolostTracksToken_;
  const StringCutObjectSelector<pat::PackedCandidate> isotrk_selection_;
  const edm::EDGetTokenT<reco::BeamSpot> beamspot_;
};

void LambdaBToLambdaCMuBuilder::produce(edm::StreamID, edm::Event &evt, edm::EventSetup const &) const
{

  //input
  edm::Handle<pat::MuonCollection> muons;
  evt.getByToken(muons_, muons);

  edm::Handle<TransientTrackCollection> muon_ttracks;
  evt.getByToken(muon_ttracks_, muon_ttracks);

  edm::Handle<pat::CompositeCandidateCollection> lambdas_c;
  evt.getByToken(lambdas_c_, lambdas_c);

  edm::Handle<TransientTrackCollection> lambdas_c_ttracks;
  evt.getByToken(lambdas_c_ttracks_, lambdas_c_ttracks);

  edm::Handle<reco::BeamSpot> beamspot;
  evt.getByToken(beamspot_, beamspot);

  //for isolation
  edm::Handle<pat::PackedCandidateCollection> iso_tracks;
  evt.getByToken(isotracksToken_, iso_tracks);

  edm::Handle<pat::PackedCandidateCollection> iso_lostTracks;
  evt.getByToken(isolostTracksToken_, iso_lostTracks);

  unsigned int nTracks = iso_tracks->size();
  unsigned int totalTracks = nTracks + iso_lostTracks->size();

  // output
  std::unique_ptr<pat::CompositeCandidateCollection> ret_val(new pat::CompositeCandidateCollection());

  for (size_t lambda_c_idx = 0; lambda_c_idx < lambdas_c->size(); ++lambda_c_idx)
  {
    // both Lambda_C already passed cuts; no need for more preselection
    edm::Ptr<pat::CompositeCandidate> lambda_c_ptr(lambdas_c, lambda_c_idx);
    edm::Ptr<reco::Candidate> trk1_ptr = lambda_c_ptr->userCand("trk1");
    edm::Ptr<reco::Candidate> trk2_ptr = lambda_c_ptr->userCand("trk2");
    edm::Ptr<reco::Candidate> trk3_ptr = lambda_c_ptr->userCand("trk3");
    int trk1_idx = lambda_c_ptr->userInt("trk1_idx");
    int trk2_idx = lambda_c_ptr->userInt("trk2_idx");
    int trk3_idx = lambda_c_ptr->userInt("trk3_idx");

    for (size_t muon_idx = 0; muon_idx < muons->size(); ++muon_idx)
    {
      edm::Ptr<pat::Muon> muon_ptr(muons, muon_idx);
      if (!muonSelection_(*muon_ptr))
        continue;
      // edm::Ptr<reco::Candidate> l1_ptr = muon_ptr->userCand("l1");
      // edm::Ptr<reco::Candidate> l2_ptr = muon_ptr->userCand("l2");
      // int l1_idx = muon_ptr->userInt("l1_idx");
      // int l2_idx = muon_ptr->userInt("l2_idx");

      // filter on charge
      if (std::abs(muon_ptr->charge() + lambda_c_ptr->charge()) != 0)
        continue;

      // Lambda_B candidate
      pat::CompositeCandidate cand;
      cand.setP4(muon_ptr->p4() + lambda_c_ptr->p4());
      cand.setCharge(0); //Lambda_B has 0 charge

      // save daughters - unfitted
      cand.addUserCand("trk1", trk1_ptr);
      cand.addUserCand("trk2", trk2_ptr);
      cand.addUserCand("trk3", trk3_ptr);
      cand.addUserCand("lambda_c", lambda_c_ptr);
      cand.addUserCand("muon", muon_ptr);

      // save indices
      cand.addUserInt("muon_idx", muon_idx);
      cand.addUserInt("trk1_idx", trk1_idx);
      cand.addUserInt("trk2_idx", trk2_idx);
      cand.addUserInt("trk3_idx", trk3_idx);
      cand.addUserInt("lambda_c_idx", lambda_c_idx);

      auto dr_info = min_max_dr({muon_ptr, trk1_ptr, trk2_ptr, trk3_ptr});
      cand.addUserFloat("min_dr", dr_info.first);
      cand.addUserFloat("max_dr", dr_info.second);

      // check if pass pre vertex cut
      if (!pre_vtx_selection_(cand))
        continue;

      KinVtxFitter fitter(
          {lambdas_c_ttracks->at(trk1_idx), lambdas_c_ttracks->at(trk2_idx),
           lambdas_c_ttracks->at(trk3_idx), muon_ttracks->at(muon_idx)},
          {trk1_ptr->mass(), trk2_ptr->mass(), trk3_ptr->mass(), muon_ptr->mass()},
          {PROTON_SIGMA, K_SIGMA, PI_SIGMA, LEP_SIGMA} //K_SIGMA == PI_SIGMA == PROTON_SIGMA
      );

      if (!fitter.success())
        continue;

      // Lambda_B position
      cand.setVertex(
          reco::Candidate::Point(
              fitter.fitted_vtx().x(),
              fitter.fitted_vtx().y(),
              fitter.fitted_vtx().z()));

      // vertex vars
      cand.addUserFloat("sv_chi2", fitter.chi2());
      cand.addUserFloat("sv_ndof", fitter.dof());
      cand.addUserFloat("sv_prob", fitter.prob());

      // refitted kinematic vars
      cand.addUserFloat("fitted_lambda_c_mass", (fitter.daughter_p4(0) + fitter.daughter_p4(1) + fitter.daughter_p4(2)).mass());
      cand.addUserFloat("fitted_lambda_c_pt", (fitter.daughter_p4(0) + fitter.daughter_p4(1) + fitter.daughter_p4(2)).pt());
      cand.addUserFloat("fitted_lambda_c_eta", (fitter.daughter_p4(0) + fitter.daughter_p4(1) + fitter.daughter_p4(2)).eta());
      cand.addUserFloat("fitted_lambda_c_phi", (fitter.daughter_p4(0) + fitter.daughter_p4(1) + fitter.daughter_p4(2)).phi());

      auto fit_p4 = fitter.fitted_p4();
      cand.addUserFloat("fitted_pt", fit_p4.pt());
      cand.addUserFloat("fitted_eta", fit_p4.eta());
      cand.addUserFloat("fitted_phi", fit_p4.phi());
      cand.addUserFloat("fitted_mass", fit_p4.mass());
      cand.addUserFloat("fitted_massErr", sqrt(fitter.fitted_candidate().kinematicParametersError().matrix()(6, 6)));

      // refitted daughters (leptons/tracks)
      std::vector<std::string> dnames{"trk1", "trk2", "trk3", "muon"};

      for (size_t idaughter = 0; idaughter < dnames.size(); idaughter++)
      {
        cand.addUserFloat("fitted_" + dnames[idaughter] + "_pt", fitter.daughter_p4(idaughter).pt());
        cand.addUserFloat("fitted_" + dnames[idaughter] + "_eta", fitter.daughter_p4(idaughter).eta());
        cand.addUserFloat("fitted_" + dnames[idaughter] + "_phi", fitter.daughter_p4(idaughter).phi());
      }

      // other vars
      cand.addUserFloat(
          "cos_theta_2D",
          cos_theta_2D(fitter, *beamspot, cand.p4()));
      cand.addUserFloat(
          "fitted_cos_theta_2D",
          cos_theta_2D(fitter, *beamspot, fit_p4));

      auto lxy = l_xy(fitter, *beamspot);
      cand.addUserFloat("l_xy", lxy.value());
      cand.addUserFloat("l_xy_unc", lxy.error());

      // post fit selection
      if (!post_vtx_selection_(cand))
        continue;

      //compute isolation
      float muon_iso03 = 0;
      float muon_iso04 = 0;
      float tk1_iso03 = 0;
      float tk1_iso04 = 0;
      float tk2_iso03 = 0;
      float tk2_iso04 = 0;
      float tk3_iso03 = 0;
      float tk3_iso04 = 0;
      float lambda_b_iso03 = 0;
      float lambda_b_iso04 = 0;

      for (unsigned int iTrk = 0; iTrk < totalTracks; ++iTrk)
      {
        const pat::PackedCandidate &trk = (iTrk < nTracks) ? (*iso_tracks)[iTrk] : (*iso_lostTracks)[iTrk - nTracks];
        // define selections for iso tracks (pT, eta, ...)
        if (!isotrk_selection_(trk))
          continue;
        // check if the track is the kaon, the proton or the pion
        if (trk1_ptr == edm::Ptr<reco::Candidate>(iso_tracks, iTrk))
          continue;
        if (trk2_ptr == edm::Ptr<reco::Candidate>(iso_tracks, iTrk))
          continue;
        if (trk3_ptr == edm::Ptr<reco::Candidate>(iso_tracks, iTrk))
          continue;
        // check if the track is one of the two leptons
        if (track_to_lepton_match(muon_ptr, iso_tracks.id(), iTrk))
          continue;

        // add to final particle iso if dR < cone
        float dr_to_muon = deltaR(cand.userFloat("fitted_muon_eta"), cand.userFloat("fitted_muon_phi"), trk.eta(), trk.phi());
        float dr_to_tk1 = deltaR(cand.userFloat("fitted_trk1_eta"), cand.userFloat("fitted_trk1_phi"), trk.eta(), trk.phi());
        float dr_to_tk2 = deltaR(cand.userFloat("fitted_trk2_eta"), cand.userFloat("fitted_trk2_phi"), trk.eta(), trk.phi());
        float dr_to_tk3 = deltaR(cand.userFloat("fitted_trk3_eta"), cand.userFloat("fitted_trk3_phi"), trk.eta(), trk.phi());
        float dr_to_lambda_b = deltaR(cand.userFloat("fitted_eta"), cand.userFloat("fitted_phi"), trk.eta(), trk.phi());

        // muon iso
        if (dr_to_muon < 0.4)
        {
          muon_iso04 += trk.pt();
          if (dr_to_muon < 0.3)
            muon_iso03 += trk.pt();
        }

        // tk1 iso
        if (dr_to_tk1 < 0.4)
        {
          tk1_iso04 += trk.pt();
          if (dr_to_tk1 < 0.3)
            tk1_iso03 += trk.pt();
        }

        // tk1 iso
        if (dr_to_tk2 < 0.4)
        {
          tk2_iso04 += trk.pt();
          if (dr_to_tk2 < 0.3)
            tk2_iso03 += trk.pt();
        }
        // tk1 iso
        if (dr_to_tk3 < 0.4)
        {
          tk3_iso04 += trk.pt();
          if (dr_to_tk3 < 0.3)
            tk3_iso03 += trk.pt();
        }

        // lambda_b iso
        if (dr_to_lambda_b < 0.4)
        {
          lambda_b_iso04 += trk.pt();
          if (dr_to_lambda_b < 0.3)
            lambda_b_iso03 += trk.pt();
        }
      }
      cand.addUserFloat("muon_iso03", muon_iso03);
      cand.addUserFloat("muon_iso04", muon_iso04);

      cand.addUserFloat("tk1_iso03", tk1_iso03);
      cand.addUserFloat("tk1_iso04", tk1_iso04);

      cand.addUserFloat("tk2_iso03", tk2_iso03);
      cand.addUserFloat("tk2_iso04", tk2_iso04);

      cand.addUserFloat("tk3_iso03", tk2_iso03);
      cand.addUserFloat("tk3_iso04", tk2_iso04);

      cand.addUserFloat("lambda_b_iso03", lambda_b_iso03);
      cand.addUserFloat("lambda_b_iso04", lambda_b_iso04);

      ret_val->push_back(cand);

    } // for(size_t muon_idx = 0; muon_idx < muons->size(); ++muon_idx) {

  } // for(size_t k_idx = 0; k_idx < lambdas_c->size(); ++k_idx)

  evt.put(std::move(ret_val));
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(LambdaBToLambdaCMuBuilder);
