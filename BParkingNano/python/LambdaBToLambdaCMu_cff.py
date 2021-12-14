import FWCore.ParameterSet.Config as cms
from PhysicsTools.BParkingNano.common_cff import *

########## inputs preparation ################

# muonsForLambdaCMu = cms.EDProducer(
#     'SingleMuonBuilder',
#     src=cms.InputTag('muonTrgSelector', 'SelectedMuons'),
#     transientTracksSrc=cms.InputTag('muonTrgSelector',
#                                     'SelectedTransientMuons'),
#     lep1Selection=cms.string('pt > 1.5'),
#     lep2Selection=cms.string(''),
#     preVtxSelection=cms.string(
#         'abs(userCand("l1").vz - userCand("l2").vz) <= 1. && mass() < 5 '
#         '&& mass() > 0 && charge() == 0 && userFloat("lep_deltaR") > 0.03'),
#     postVtxSelection=cms.string(
#         'userFloat("sv_chi2") < 998 && userFloat("sv_prob") > 1.e-5'),
# )

# Lambda_C
LambdaCToProtonKPi = cms.EDProducer(
    'LambdaCBuilder',
    pfcands=cms.InputTag('tracksBPark', 'SelectedTracks'),
    transientTracks=cms.InputTag('tracksBPark', 'SelectedTransientTracks'),
    trk1Selection=cms.string('pt > 0.5 && abs(eta)<2.4'),  #need optimization   
    trk2Selection=cms.string('pt > 0.5 && abs(eta)<2.4'),  #need optimization
    trk3Selection=cms.string('pt > 0.5 && abs(eta)<2.4'),  #need optimization
    preVtxSelection=cms.string(
        '(abs(userCand("trk1").vz - userCand("trk2").vz) < 1.0 && abs(userCand("trk1").vz - userCand("trk3").vz) < 1.0 && abs(userCand("trk2").vz - userCand("trk3").vz) < 1.0)'
        '&& pt() > 2.0 '
        '&& (  (mass() < 2.4 && mass() > 2.1) '
        '|| (userFloat("1barMass") < 2.4 && userFloat("1barMass") > 2.1) '
        '|| (userFloat("2barMass") < 2.4 && userFloat("2barMass") > 2.1) '
        '|| (userFloat("3barMass") < 2.4 && userFloat("3barMass") > 2.1)  )'),
    postVtxSelection=cms.string(
        'userFloat("sv_prob") > 1.e-5'
        ' && (  (userFloat("fitted_mass")<2.4 && userFloat("fitted_mass")>2.1)'
        ' || (userFloat("fitted_1barMass")<2.4 && userFloat("fitted_1barMass")>2.1)'
        ' || (userFloat("fitted_2barMass")<2.4 && userFloat("fitted_2barMass")>2.1)'
        ' || (userFloat("fitted_3barMass")<2.4 && userFloat("fitted_3barMass")>2.1)  )'
    ))

# Lambda_B
LambdaBToLambdaCMu = cms.EDProducer(
    'LambdaBToLambdaCMuBuilder',
    muons=cms.InputTag('muonTrgSelector', 'SelectedMuons'),
    muonTransientTracks=cms.InputTag('muonTrgSelector',
                                     'SelectedTransientMuons'),
    muonSelection=cms.string('pt > 1.5'),
    lambdas_c=cms.InputTag('LambdaCToProtonKPi'),
    lambdas_cTransientTracks=cms.InputTag('tracksBPark',
                                          'SelectedTransientTracks'),
    tracks=cms.InputTag("packedPFCandidates"),
    lostTracks=cms.InputTag("lostTracks"),
    isoTracksSelection=cms.string('pt > 0.7 && abs(eta)<2.5'),
    beamSpot=cms.InputTag("offlineBeamSpot"),
    preVtxSelection=cms.string(
        'pt > 3. && userFloat("min_dr") > 0.03'
        '&& ( (mass < 7. && mass > 4.) '
        '|| (userFloat("1barMass")<7. && userFloat("1barMass")>4.) '
        '|| (userFloat("2barMass")<7. && userFloat("2barMass")>4.) '
        '|| (userFloat("3barMass")<7. && userFloat("3barMass")>4.) )'),
    postVtxSelection=cms.string(
        'userFloat("sv_prob") > 0.001 '
        '&& userFloat("fitted_cos_theta_2D") >= 0'
        '&& ( (userFloat("fitted_mass") > 4 && userFloat("fitted_mass") < 7.)'
        '|| (userFloat("fitted_1barMass") > 4 && userFloat("fitted_1barMass") < 7.) '
        '|| (userFloat("fitted_2barMass") > 4 && userFloat("fitted_2barMass") < 7.)  '
        '|| (userFloat("fitted_3barMass") > 4 && userFloat("fitted_3barMass") < 7.) )'
    ))

################################### Tables #####################################

LambdaCToProtonKPiTable = cms.EDProducer(
    'SimpleCompositeCandidateFlatTableProducer',
    src=cms.InputTag("LambdaCToProtonKPi"),
    cut=cms.string(""),
    name=cms.string("LambdaC"),
    doc=cms.string("LambdaC Variables"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables=cms.PSet(CandVars,
                       nofit_1barMass=ufloat('1barMass'),
                       nofit_2barMass=ufloat('2barMass'),
                       nofit_3barMass=ufloat('3barMass'),
                       fitted_mass=ufloat('fitted_mass'),
                       fitted_1barMass=ufloat('fitted_1barMass'),
                       fitted_2barMass=ufloat('fitted_2barMass'),
                       fitted_3barMass=ufloat('fitted_3barMass'),
                       fitted_pt=ufloat('fitted_pt'),
                       fitted_eta=ufloat('fitted_eta'),
                       fitted_phi=ufloat('fitted_phi'),
                       svprob=ufloat('sv_prob'),
                       trk_deltaR12=ufloat('trk_deltaR12'),
                       trk_deltaR13=ufloat('trk_deltaR13'),
                       trk_deltaR23=ufloat('trk_deltaR23'),
                       trk1_idx=uint('trk1_idx'),
                       trk2_idx=uint('trk2_idx'),
                       trk3_idx=uint('trk3_idx')))

LambdaBToLambdaCMuTable = cms.EDProducer(
    'SimpleCompositeCandidateFlatTableProducer',
    src=cms.InputTag("LambdaBToLambdaCMu"),
    cut=cms.string(""),
    name=cms.string("LambdaBToLambdaCMu"),
    doc=cms.string("LambdaBToLambdaCMu Variables"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables=cms.PSet(
        # pre-fit quantities
        CandVars,
        muon_idx=uint('muon_idx'),
        trk1_idx=uint('trk1_idx'),
        trk2_idx=uint('trk2_idx'),
        trk3_idx=uint('trk3_idx'),
        lambda_c_idx=uint('lambda_c_idx'),
        min_dr=ufloat('min_dr'),
        max_dr=ufloat('max_dr'),
        # fit and vtx info
        chi2=ufloat('sv_chi2'),
        svprob=ufloat('sv_prob'),
        l_xy=ufloat('l_xy'),
        l_xy_unc=ufloat('l_xy_unc'),
        fit_lambda_c_mass=ufloat('fitted_lambda_c_mass'),
        fit_lambda_c_pt=ufloat('fitted_lambda_c_pt'),
        fit_lambda_c_eta=ufloat('fitted_lambda_c_eta'),
        fit_lambda_c_phi=ufloat('fitted_lambda_c_phi'),
        # Cos(theta)
        cos2D=ufloat('cos_theta_2D'),
        fit_cos2D=ufloat('fitted_cos_theta_2D'),
        # post-fit momentum
        fit_mass=ufloat('fitted_mass'),
        fit_massErr=ufloat('fitted_massErr'),
        fit_pt=ufloat('fitted_pt'),
        fit_eta=ufloat('fitted_eta'),
        fit_phi=ufloat('fitted_phi'),

        # 2dn mass hypothesis
        nofit_1barMass=ufloat('1barMass'),
        fit_1barMass=ufloat('fitted_1barMass'),
        fit_1barLambda_c_mass=ufloat('1barMasslambda_c_fullfit'),

        # 3rd mass hypothesis
        nofit_2barMass=ufloat('2barMass'),
        fit_2barMass=ufloat('fitted_2barMass'),
        fit_2barLambda_c_mass=ufloat('2barMasslambda_c_fullfit'),

        # 4th mass hypothesis
        nofit_3barMass=ufloat('3barMass'),
        fit_3barMass=ufloat('fitted_3barMass'),
        fit_3barLambda_c_mass=ufloat('3barMasslambda_c_fullfit'),

        # post-fit tracks/leptons
        #muon
        fit_muon_pt=ufloat('fitted_muon_pt'),
        fit_muon_eta=ufloat('fitted_muon_eta'),
        fit_muon_phi=ufloat('fitted_muon_phi'),

        #trk1
        fit_trk1_pt=ufloat('fitted_trk1_pt'),
        fit_trk1_eta=ufloat('fitted_trk1_eta'),
        fit_trk1_phi=ufloat('fitted_trk1_phi'),

        #trk2
        fit_trk2_pt=ufloat('fitted_trk2_pt'),
        fit_trk2_eta=ufloat('fitted_trk2_eta'),
        fit_trk2_phi=ufloat('fitted_trk2_phi'),

        #trk3
        fit_trk3_pt=ufloat('fitted_trk3_pt'),
        fit_trk3_eta=ufloat('fitted_trk3_eta'),
        fit_trk3_phi=ufloat('fitted_trk3_phi'),

        # isolation
        muon_iso03=ufloat('muon_iso03'),
        muon_iso04=ufloat('muon_iso04'),
        tk1_iso03=ufloat('tk1_iso03'),
        tk1_iso04=ufloat('tk1_iso04'),
        tk2_iso03=ufloat('tk2_iso03'),
        tk2_iso04=ufloat('tk2_iso04'),
        tk3_iso03=ufloat('tk3_iso03'),
        tk3_iso04=ufloat('tk3_iso04'),
        lambda_b_iso03=ufloat('lambda_b_iso03'),
        lambda_b_iso04=ufloat('lambda_b_iso04'),
    ))

CountLambdaBToLambdaCMu = cms.EDFilter("PATCandViewCountFilter",
                                       minNumber=cms.uint32(1),
                                       maxNumber=cms.uint32(999999),
                                       src=cms.InputTag("LambdaBToLambdaCMu"))

########################### Sequencies  ############################

LambdaCToProtonKPiSequence = cms.Sequence(LambdaCToProtonKPi)

# LambdaBToLambdaCMuSequence = cms.Sequence(
#     (muonsForLambdaCMu * LambdaBToLambdaCMu))

LambdaBToLambdaCMuSequence = cms.Sequence(LambdaBToLambdaCMu)
