// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
///
/// \file taskSingleElectron.cxx
/// \brief task for electrons from heavy-flavour hadron decays
/// \author Jonghan Park (Jeonbuk National University), Seul I Jeong (Pusan National University)

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/Expressions.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/SliceCache.h>
#include <Framework/runDataProcessing.h>

#include <TPDGCode.h>

#include <map>

using namespace o2;
using namespace o2::constants::math;
using namespace o2::constants::physics;
using namespace o2::framework;
using namespace o2::framework::expressions;

enum SourceType {
  NotElec = 0,      // not electron
  DirectCharm = 1,  // electrons from prompt charm hadrons
  DirectBeauty = 2, // electrons from primary beauty hadrons
  BeautyCharm = 3,  // electrons from non-prompt charm hadrons
  DirectGamma = 4,  // electrons from direct photon
  GammaPi0 = 5,
  GammaEta = 6,
  GammaOmega = 7,
  GammaPhi = 8,
  GammaEtaPrime = 9,
  GammaRho0 = 10,
  GammaK0s = 11,
  GammaK0l = 12,
  GammaKe3 = 13,
  GammaLambda0 = 14,
  GammaSigma = 15,
  Pi0 = 16,
  Eta = 17,
  Omega = 18,
  Phi = 19,
  EtaPrime = 20,
  Rho0 = 21,
  K0s = 22,
  K0l = 23,
  Ke3 = 24,
  Lambda0 = 25,
  Sigma = 26,
  Else = 27
};

struct HfTaskSingleElectron {

  // Produces

  // Configurable
  Configurable<int> nContribMin{"nContribMin", 2, "min number of contributors"};
  Configurable<float> posZMax{"posZMax", 10., "max posZ cut"};
  Configurable<float> ptTrackMax{"ptTrackMax", 10., "max pt cut"};
  Configurable<float> ptTrackMin{"ptTrackMin", 0.5, "min pt cut"};
  Configurable<float> etaTrackMax{"etaTrackMax", 0.8, "eta cut"};
  Configurable<int> nCrossedRowTpcMin{"nCrossedRowTpcMin", 70, "min # of TPC n cluster crossed rows"};
  Configurable<float> nClsFoundOverFindableTpcMin{"nClsFoundOverFindableTpcMin", 0.8, "min # of TPC found/findable clusters"};
  Configurable<float> chi2PerNClTpcMax{"chi2PerNClTpcMax", 4., "max # of tpc chi2 per clusters"};
  Configurable<int> clsIbItsMin{"clsIbItsMin", 3, "min # of its clusters in IB"};
  Configurable<float> chi2PerNClItsMax{"chi2PerNClItsMax", 6., "min # of its chi2 per clusters"};
  Configurable<float> dcaxyMax{"dcaxyMax", 1., "max of track dca in xy"};
  Configurable<float> dcazMax{"dcazMax", 2., "max of track dca in z"};
  Configurable<float> nSigmaTofMax{"nSigmaTofMax", 3., "max of tof nsigma"};
  Configurable<float> nSigmaTpcMin{"nSigmaTpcMin", -1., "min of tpc nsigma"};
  Configurable<float> nSigmaTpcMax{"nSigmaTpcMax", 3., "max of tpc nsigma"};

  Configurable<int> nBinsP{"nBinsP", 1500, "number of bins of particle momentum"};
  Configurable<int> nBinsPt{"nBinsPt", 100, "N bins in pT histo"};

  Configurable<int> nSigmaTpcHadronMax{"nSigmaTpcHadronMax", -3, "max of tpc hadron nsigma"};
  Configurable<int> nSigmaTpcHadronMin{"nSigmaTpcHadronMin", -5, "min of tpc hadron nsigma"};

  // SliceCache
  SliceCache cache;

  // using declarations
  using MyCollisions = soa::Join<aod::Collisions, aod::EvSels>;
  using TracksEl = soa::Join<aod::Tracks, aod::TrackSelection, aod::TrackSelectionExtension, aod::TracksExtra, aod::TracksDCA, aod::pidTOFFullEl, aod::pidTPCFullEl, aod::pidTOFbeta>;
  using McTracksEl = soa::Join<aod::Tracks, aod::TrackExtra, aod::TracksDCA, aod::pidTOFFullEl, aod::pidTPCFullEl, aod::McTrackLabels, aod::pidTOFbeta>;

  // Filter
  Filter collZFilter = nabs(aod::collision::posZ) < posZMax;

  // Partition

  // ConfigurableAxis
  ConfigurableAxis axisPtEl{"axisPtEl", {VARIABLE_WIDTH, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.75f, 2.0f, 2.25f, 2.5f, 2.75f, 3.f, 3.5f, 4.0f, 5.0f, 6.0f, 8.0f, 10.0f}, "electron pt bins"};

  // Histogram registry
  HistogramRegistry histos{"histos"};

  void init(InitContext const&)
  {
    // AxisSpec
    const AxisSpec axisEvt{4, 0., 4., "nEvents"};
    const AxisSpec axisNCont{100, 0., 100., "nCont"};
    const AxisSpec axisPosZ{600, -30., 30., "Z_{pos}"};
    const AxisSpec axisEta{30, -1.5, +1.5, "#eta"};
    const AxisSpec axisP{nBinsP, 0., 15., "p_{T}"};
    const AxisSpec axisPt{nBinsPt, 0., 15., "p_{T}"};
    const AxisSpec axisNsig{800, -20., 20.};
    const AxisSpec axisTrackIp{4000, -0.2, 0.2, "dca"};
    const AxisSpec axisBeta{600, 0., 1.2, "#beta"};                  // reused from work/32 tof_quality_study
    const AxisSpec axisBetaP{300, 0., 12., "p (GeV/#it{c})"};        // coarser p axis for beta 2D (species split)

    // create histograms
    histos.add("hNEvents", "Number of events", kTH1D, {{1, 0., 1.}});
    histos.add("hVtxZ", "hVtxZ; cm; entries", kTH1D, {axisPosZ});
    histos.add("hEtaTrack", "hEtaTrack; #eta; entries", kTH1D, {axisEta});
    histos.add("hPtTrack", "#it{p}_{T} distribution of selected tracks; #it{p}_{T} (GeV/#it{c}); entries", kTH1D, {axisPt});

    // QA plots for trigger track selection
    histos.add("hNClsTpcTrack", "hNClsTpcTrack", kTH1D, {{200, 0, 200}});
    histos.add("hNClsFoundFindableTpcTrack", "", kTH1D, {{10, 0, 1}});
    histos.add("hChi2TpcTrack", "", kTH1D, {{100, 0, 10}});
    histos.add("hIbClsItsTrack", "", kTH1D, {{10, 0, 10}});
    histos.add("hChi2ItsTrack", "", kTH1D, {{50, 0, 50}});
    histos.add("hDcaXYTrack", "", kTH1D, {{600, -3, 3}});
    histos.add("hDcaZTrack", "", kTH1D, {{600, -3, 3}});

    // ===== [pre-cut QA] BEGIN: track-quality distributions BEFORE trackSel (all charged tracks), for data/MC comparison =====
    histos.add("hEtaTrack_noCut", "#eta before track sel; #eta; entries", kTH1D, {axisEta});
    histos.add("hPtTrack_noCut", "#it{p}_{T} before track sel; #it{p}_{T} (GeV/#it{c}); entries", kTH1D, {axisPt});
    histos.add("hNClsTpcTrack_noCut", "TPC crossed rows before track sel; N_{crossed rows}^{TPC}; entries", kTH1D, {{200, 0, 200}});
    histos.add("hNClsFoundFindableTpcTrack_noCut", "TPC crossed-rows/findable before track sel; N_{crossed rows}/N_{findable}; entries", kTH1D, {{110, 0, 1.1}});
    histos.add("hChi2TpcTrack_noCut", "TPC #chi^{2}/cluster before track sel; #chi^{2}/N_{cls}^{TPC}; entries", kTH1D, {{100, 0, 10}});
    histos.add("hIbClsItsTrack_noCut", "ITS IB clusters before track sel; N_{cls}^{ITS-IB}; entries", kTH1D, {{10, 0, 10}});
    histos.add("hChi2ItsTrack_noCut", "ITS #chi^{2}/cluster before track sel; #chi^{2}/N_{cls}^{ITS}; entries", kTH1D, {{50, 0, 50}});
    histos.add("hDcaXYTrack_noCut", "DCA_{xy} before track sel; DCA_{xy} (cm); entries", kTH1D, {{600, -3, 3}});
    histos.add("hDcaZTrack_noCut", "DCA_{z} before track sel; DCA_{z} (cm); entries", kTH1D, {{600, -3, 3}});
    histos.add("hTpcNSigPt_noCut", "TPC n#sigma_{e} before track sel; #it{p}_{T} (GeV/#it{c}); n#sigma_{e}^{TPC}", kTH2D, {{axisPtEl}, {axisNsig}});
    histos.add("hTofNSigPt_noCut", "TOF n#sigma_{e} before track sel; #it{p}_{T} (GeV/#it{c}); n#sigma_{e}^{TOF}", kTH2D, {{axisPtEl}, {axisNsig}});
    // ===== [pre-cut QA] END =====

    // ===== [pre-cut QA 2D + species] BEGIN: correlations and MC-truth species, before trackSel =====
    // 2D correlations (filled in both data and MC)
    histos.add("hTpcNSigEta_noCut", "n#sigma_{e}^{TPC} vs #eta before track sel; #eta; n#sigma_{e}^{TPC}", kTH2D, {axisEta, axisNsig});                               // eta asymmetry
    histos.add("hNClsCrossedRowsEta_noCut", "TPC crossed rows vs #eta before track sel; #eta; N_{crossed rows}^{TPC}", kTH2D, {axisEta, {200, 0, 200}});             // cluster diffusion (AN#1783 <Ncl> vs eta)
    histos.add("hDcaXYIbClsIts_noCut", "DCA_{xy} vs ITS-IB clusters before track sel; N_{cls}^{ITS-IB}; DCA_{xy} (cm)", kTH2D, {{10, 0, 10}, axisTrackIp});          // d0 resolution vs IB requirement
    histos.add("hTpcNSigNClsCrossedRows_noCut", "n#sigma_{e}^{TPC} vs crossed rows vs p_{T} before track sel; #it{p}_{T} (GeV/#it{c}); N_{crossed rows}^{TPC}; n#sigma_{e}^{TPC}", kTH3F, {axisPtEl, {80, 0, 160}, axisNsig}); // e/pi separation vs crossed rows (momentum-sliceable)
    // MC-truth species (filled only in processMc) — reused from work/15 MCInformation_forEID
    histos.add("hTpcNSigPt_noCut_Ele", "n#sigma_{e}^{TPC} truth e before track sel; #it{p}_{T} (GeV/#it{c}); n#sigma_{e}^{TPC}", kTH2D, {axisPtEl, axisNsig});
    histos.add("hTpcNSigPt_noCut_Pi", "n#sigma_{e}^{TPC} truth #pi before track sel; #it{p}_{T} (GeV/#it{c}); n#sigma_{e}^{TPC}", kTH2D, {axisPtEl, axisNsig});
    histos.add("hTpcNSigPt_noCut_K", "n#sigma_{e}^{TPC} truth K before track sel; #it{p}_{T} (GeV/#it{c}); n#sigma_{e}^{TPC}", kTH2D, {axisPtEl, axisNsig});
    histos.add("hTpcNSigPt_noCut_Pro", "n#sigma_{e}^{TPC} truth p before track sel; #it{p}_{T} (GeV/#it{c}); n#sigma_{e}^{TPC}", kTH2D, {axisPtEl, axisNsig});
    histos.add("hTpcNSigPt_noCut_Other", "n#sigma_{e}^{TPC} truth other before track sel; #it{p}_{T} (GeV/#it{c}); n#sigma_{e}^{TPC}", kTH2D, {axisPtEl, axisNsig});
    histos.add("hDcaZBeauty_noCut", "DCA_{z} of beauty-decay electrons (truth) before track sel; DCA_{z} (cm); entries", kTH1D, {{600, -3, 3}}); // DCA_z signal-loss check
    // TOF beta vs p (reused from work/32 tof_quality_study) + TOF nSigma / beta truth-species
    histos.add("hTofBetaP_noCut", "TOF #beta vs p before track sel; p (GeV/#it{c}); #beta", kTH2D, {axisBetaP, axisBeta});
    histos.add("hTofNSigPt_noCut_Ele", "TOF n#sigma_{e} truth e before track sel; #it{p}_{T} (GeV/#it{c}); n#sigma_{e}^{TOF}", kTH2D, {axisPtEl, axisNsig});
    histos.add("hTofNSigPt_noCut_Pi", "TOF n#sigma_{e} truth #pi before track sel; #it{p}_{T} (GeV/#it{c}); n#sigma_{e}^{TOF}", kTH2D, {axisPtEl, axisNsig});
    histos.add("hTofNSigPt_noCut_K", "TOF n#sigma_{e} truth K before track sel; #it{p}_{T} (GeV/#it{c}); n#sigma_{e}^{TOF}", kTH2D, {axisPtEl, axisNsig});
    histos.add("hTofNSigPt_noCut_Pro", "TOF n#sigma_{e} truth p before track sel; #it{p}_{T} (GeV/#it{c}); n#sigma_{e}^{TOF}", kTH2D, {axisPtEl, axisNsig});
    histos.add("hTofNSigPt_noCut_Other", "TOF n#sigma_{e} truth other before track sel; #it{p}_{T} (GeV/#it{c}); n#sigma_{e}^{TOF}", kTH2D, {axisPtEl, axisNsig});
    histos.add("hTofBetaP_noCut_Ele", "TOF #beta truth e before track sel; p (GeV/#it{c}); #beta", kTH2D, {axisBetaP, axisBeta});
    histos.add("hTofBetaP_noCut_Pi", "TOF #beta truth #pi before track sel; p (GeV/#it{c}); #beta", kTH2D, {axisBetaP, axisBeta});
    histos.add("hTofBetaP_noCut_K", "TOF #beta truth K before track sel; p (GeV/#it{c}); #beta", kTH2D, {axisBetaP, axisBeta});
    histos.add("hTofBetaP_noCut_Pro", "TOF #beta truth p before track sel; p (GeV/#it{c}); #beta", kTH2D, {axisBetaP, axisBeta});
    histos.add("hTofBetaP_noCut_Other", "TOF #beta truth other before track sel; p (GeV/#it{c}); #beta", kTH2D, {axisBetaP, axisBeta});
    // ===== [pre-cut QA 2D + species] END =====

    // pid
    histos.add("hTofNSigPt", "", kTH2D, {{axisPtEl}, {axisNsig}});
    histos.add("hTofNSigPtQA", "", kTH2D, {{axisPtEl}, {axisNsig}});
    histos.add("hTpcNSigP", "", kTH2D, {{axisP}, {axisNsig}});
    histos.add("hTpcNSigPt", "", kTH2D, {{axisPtEl}, {axisNsig}});
    histos.add("hTpcNSigPAfterTofCut", "", kTH2D, {{axisP}, {axisNsig}});
    histos.add("hTpcNSigPtAfterTofCut", "", kTH2D, {{axisPtEl}, {axisNsig}});
    histos.add("hTpcNSigPtQA", "", kTH2D, {{axisPtEl}, {axisNsig}});

    // track impact parameter
    histos.add("hDcaTrack", "", kTH2D, {{axisPtEl}, {axisTrackIp}});
    histos.add("hDcaBeauty", "", kTH2D, {{axisPtEl}, {axisTrackIp}});
    histos.add("hDcaCharm", "", kTH2D, {{axisPtEl}, {axisTrackIp}});
    histos.add("hDcaDalitz", "", kTH2D, {{axisPtEl}, {axisTrackIp}});
    histos.add("hDcaConv", "", kTH2D, {{axisPtEl}, {axisTrackIp}});
    histos.add("hDcaHadron", "", kTH2D, {{axisPtEl}, {axisTrackIp}});

    // QA plots for MC
    histos.add("hPdgC", "", kTH1D, {{10001, -0.5, 10000.5}});
    histos.add("hPdgB", "", kTH1D, {{10001, -0.5, 10000.5}});
    histos.add("hPdgDa", "", kTH1D, {{10001, -0.5, 10000.5}});
    histos.add("hPdgCo", "", kTH1D, {{10001, -0.5, 10000.5}});
  }

  template <typename TrackType>
  bool trackSel(const TrackType& track)
  {
    if ((track.pt() > ptTrackMax) || (track.pt() < ptTrackMin)) {
      return false;
    }
    if (std::abs(track.eta()) > etaTrackMax) {
      return false;
    }

    if (track.tpcNClsCrossedRows() < nCrossedRowTpcMin) {
      return false;
    }

    if (track.tpcCrossedRowsOverFindableCls() < nClsFoundOverFindableTpcMin) {
      return false;
    }

    if (track.tpcChi2NCl() > chi2PerNClTpcMax) {
      return false;
    }

    if (!(track.itsNClsInnerBarrel() == clsIbItsMin)) {
      return false;
    }

    if (track.itsChi2NCl() > chi2PerNClItsMax) {
      return false;
    }

    if (std::abs(track.dcaXY()) > dcaxyMax) {
      return false;
    }

    if (std::abs(track.dcaZ()) > dcazMax) {
      return false;
    }

    return true;
  }

  template <typename TrackType>
  int getElecSource(const TrackType& track, double& mpt, int& mpdg)
  {
    auto mcpart = track.mcParticle();
    if (std::abs(mcpart.pdgCode()) != kElectron) {
      return NotElec;
    }

    int motherPdg = -999;
    int grmotherPdg = -999;
    int ggrmotherPdg = -999; // mother, grand mother, grand grand mother pdg
    int motherPt = -999.;
    int grmotherPt = -999;
    int ggrmotherPt = -999.; // mother, grand mother, grand grand mother pt

    auto partMother = mcpart.template mothers_as<aod::McParticles>(); // first mother particle of electron
    auto partMotherCopy = partMother;                                 // copy of the first mother
    auto mctrack = partMother;                                        // will change all the time

    motherPt = partMother.front().pt();                 // first mother pt
    motherPdg = std::abs(partMother.front().pdgCode()); // first mother pdg
    mpt = motherPt;                                     // copy of first mother pt
    mpdg = motherPdg;                                   // copy of first mother pdg

    // check if electron from charm hadrons
    if ((static_cast<int>(motherPdg / 100.) % 10) == kCharm || (static_cast<int>(motherPdg / 1000.) % 10) == kCharm) {

      // iterate until B hadron is found as an ancestor
      while (partMother.size()) {
        mctrack = partMother.front().template mothers_as<aod::McParticles>();
        if (mctrack.size()) {
          auto const& grmothersIdsVec = mctrack.front().mothersIds();

          if (grmothersIdsVec.empty()) {
            return DirectCharm;
          }
          grmotherPt = mctrack.front().pt();
          grmotherPdg = std::abs(mctrack.front().pdgCode());
          if ((static_cast<int>(grmotherPdg / 100.) % 10) == kBottom || (static_cast<int>(grmotherPdg / 1000.) % 10) == kBottom) {
            mpt = grmotherPt;
            mpdg = grmotherPdg;
            return BeautyCharm;
          }
        }
        partMother = mctrack;
      }
    } else if ((static_cast<int>(motherPdg / 100.) % 10) == kBottom || (static_cast<int>(motherPdg / 1000.) % 10) == kBottom) { // check if electron from beauty hadrons
      return DirectBeauty;
    } else if (motherPdg == kGamma) { // check if electron from photon conversion
      mctrack = partMother.front().template mothers_as<aod::McParticles>();
      if (mctrack.size()) {
        auto const& grmothersIdsVec = mctrack.front().mothersIds();
        if (grmothersIdsVec.empty()) {
          return DirectGamma;
        }
        grmotherPdg = std::abs(mctrack.front().pdgCode());
        mpdg = grmotherPdg;
        mpt = mctrack.front().pt();

        partMother = mctrack;
        mctrack = partMother.front().template mothers_as<aod::McParticles>();
        if (mctrack.size()) {
          auto const& ggrmothersIdsVec = mctrack.front().mothersIds();
          if (ggrmothersIdsVec.empty()) {
            if (grmotherPdg == kPi0) {
              return GammaPi0;
            }
            if (grmotherPdg == Pdg::kEta) {
              return GammaEta;
            }
            if (grmotherPdg == Pdg::kOmega) {
              return GammaOmega;
            }
            if (grmotherPdg == Pdg::kPhi) {
              return GammaPhi;
            }
            if (grmotherPdg == Pdg::kEtaPrime) {
              return GammaEtaPrime;
            }
            if (grmotherPdg == kRho770_0) {
              return GammaRho0;
            }
            return Else;
          }
          ggrmotherPdg = mctrack.front().pdgCode();
          ggrmotherPt = mctrack.front().pt();
          mpdg = ggrmotherPdg;
          mpt = ggrmotherPt;
          if (grmotherPdg == kPi0) {
            if (ggrmotherPdg == kK0Short) {
              return GammaK0s;
            }
            if (ggrmotherPdg == kK0Long) {
              return GammaK0l;
            }
            if (ggrmotherPdg == kKPlus) {
              return GammaKe3;
            }
            if (ggrmotherPdg == kLambda0) {
              return GammaLambda0;
            }
            if (ggrmotherPdg == kSigmaPlus) {
              return GammaSigma;
            }
            mpdg = grmotherPdg;
            mpt = grmotherPt;
            return GammaPi0;
          }
          if (grmotherPdg == Pdg::kEta) {
            mpdg = grmotherPdg;
            mpt = grmotherPt;
            return GammaEta;
          }
          if (grmotherPdg == Pdg::kOmega) {
            mpdg = grmotherPdg;
            mpt = grmotherPt;
            return GammaOmega;
          }
          if (grmotherPdg == Pdg::kPhi) {
            mpdg = grmotherPdg;
            mpt = grmotherPt;
            return GammaPhi;
          }
          if (grmotherPdg == Pdg::kEtaPrime) {
            mpdg = grmotherPdg;
            mpt = grmotherPt;
            return GammaEtaPrime;
          }
          if (grmotherPdg == kRho770_0) {
            mpdg = grmotherPdg;
            mpt = grmotherPt;
            return GammaRho0;
          }
          return Else;
        }
      }
    } else { // check if electron from Dalitz decays
      mctrack = partMother.front().template mothers_as<aod::McParticles>();
      if (mctrack.size()) {
        auto const& grmothersIdsVec = mctrack.front().mothersIds();
        if (grmothersIdsVec.empty()) {
          static const std::map<int, SourceType> pdgToSource = {
            {kPi0, Pi0},
            {Pdg::kEta, Eta},
            {Pdg::kOmega, Omega},
            {Pdg::kPhi, Phi},
            {Pdg::kEtaPrime, EtaPrime},
            {kRho770_0, Rho0},
            {kKPlus, Ke3},
            {kK0Long, K0l}};

          auto it = pdgToSource.find(motherPdg);
          if (it != pdgToSource.end()) {
            return it->second;
          }
          return Else;
        }
        if (motherPdg == kPi0) {
          grmotherPt = mctrack.front().pt();
          grmotherPdg = mctrack.front().pdgCode();
          mpt = grmotherPt;
          mpdg = grmotherPdg;
          if (grmotherPdg == kK0Short) {
            return K0s;
          }
          if (grmotherPdg == kK0Long) {
            return K0l;
          }
          if (grmotherPdg == kKPlus) {
            return Ke3;
          }
          if (grmotherPdg == kLambda0) {
            return Lambda0;
          }
          if (grmotherPdg == kSigmaPlus) {
            return Sigma;
          }
          mpt = motherPt;
          mpdg = motherPdg;
          return Pi0;
        }
        if (motherPdg == Pdg::kEta) {
          return Eta;
        }
        if (motherPdg == Pdg::kOmega) {
          return Omega;
        }
        if (motherPdg == Pdg::kPhi) {
          return Phi;
        }
        if (motherPdg == Pdg::kEtaPrime) {
          return EtaPrime;
        }
        if (motherPdg == kRho770_0) {
          return Rho0;
        }
        if (motherPdg == kKPlus) {
          return Ke3;
        }
        if (motherPdg == kK0Long) {
          return K0l;
        }
        return Else;
      }
    }

    return Else;
  }

  void processData(soa::Filtered<MyCollisions>::iterator const& collision,
                   TracksEl const& tracks)
  {
    float const flagAnalysedEvt = 0.5;

    if (!collision.sel8()) {
      return;
    }

    if (collision.numContrib() < nContribMin) {
      return;
    }

    histos.fill(HIST("hVtxZ"), collision.posZ());
    histos.fill(HIST("hNEvents"), flagAnalysedEvt);

    for (const auto& track : tracks) {

      // ===== [pre-cut QA] BEGIN: fill track-quality distributions BEFORE trackSel (all charged tracks) =====
      histos.fill(HIST("hEtaTrack_noCut"), track.eta());
      histos.fill(HIST("hPtTrack_noCut"), track.pt());
      histos.fill(HIST("hNClsTpcTrack_noCut"), track.tpcNClsCrossedRows());
      histos.fill(HIST("hNClsFoundFindableTpcTrack_noCut"), track.tpcCrossedRowsOverFindableCls());
      histos.fill(HIST("hChi2TpcTrack_noCut"), track.tpcChi2NCl());
      histos.fill(HIST("hIbClsItsTrack_noCut"), track.itsNClsInnerBarrel());
      histos.fill(HIST("hChi2ItsTrack_noCut"), track.itsChi2NCl());
      histos.fill(HIST("hDcaXYTrack_noCut"), track.dcaXY());
      histos.fill(HIST("hDcaZTrack_noCut"), track.dcaZ());
      histos.fill(HIST("hTpcNSigPt_noCut"), track.pt(), track.tpcNSigmaEl());
      histos.fill(HIST("hTofNSigPt_noCut"), track.pt(), track.tofNSigmaEl());
      histos.fill(HIST("hTpcNSigEta_noCut"), track.eta(), track.tpcNSigmaEl());
      histos.fill(HIST("hNClsCrossedRowsEta_noCut"), track.eta(), track.tpcNClsCrossedRows());
      histos.fill(HIST("hDcaXYIbClsIts_noCut"), track.itsNClsInnerBarrel(), track.dcaXY());
      histos.fill(HIST("hTpcNSigNClsCrossedRows_noCut"), track.pt(), track.tpcNClsCrossedRows(), track.tpcNSigmaEl());
      histos.fill(HIST("hTofBetaP_noCut"), track.p(), track.beta());
      // ===== [pre-cut QA] END =====

      if (!trackSel(track)) {
        continue;
      }

      if (!(track.passedITSRefit() && track.passedTPCRefit())) {
        continue;
      }

      histos.fill(HIST("hEtaTrack"), track.eta());
      histos.fill(HIST("hPtTrack"), track.pt());

      histos.fill(HIST("hNClsTpcTrack"), track.tpcNClsCrossedRows());
      histos.fill(HIST("hNClsFoundFindableTpcTrack"), track.tpcCrossedRowsOverFindableCls());
      histos.fill(HIST("hChi2TpcTrack"), track.tpcChi2NCl());
      histos.fill(HIST("hIbClsItsTrack"), track.itsNClsInnerBarrel());
      histos.fill(HIST("hChi2ItsTrack"), track.itsChi2NCl());
      histos.fill(HIST("hDcaXYTrack"), track.dcaXY());
      histos.fill(HIST("hDcaZTrack"), track.dcaZ());

      histos.fill(HIST("hTofNSigPt"), track.pt(), track.tofNSigmaEl());
      histos.fill(HIST("hTpcNSigP"), track.p(), track.tpcNSigmaEl());
      histos.fill(HIST("hTpcNSigPt"), track.pt(), track.tpcNSigmaEl());

      if (std::abs(track.tofNSigmaEl()) > nSigmaTofMax) {
        continue;
      }
      histos.fill(HIST("hTofNSigPtQA"), track.pt(), track.tofNSigmaEl());
      histos.fill(HIST("hTpcNSigPAfterTofCut"), track.p(), track.tpcNSigmaEl());
      histos.fill(HIST("hTpcNSigPtAfterTofCut"), track.pt(), track.tpcNSigmaEl());

      if (track.tpcNSigmaEl() < nSigmaTpcMin || track.tpcNSigmaEl() > nSigmaTpcMax) {
        continue;
      }

      if (track.tpcNSigmaEl() < nSigmaTpcHadronMax && track.tpcNSigmaEl() > nSigmaTpcHadronMin) {

        histos.fill(HIST("hDcaHadron"), track.pt(), track.dcaXY());
      }

      histos.fill(HIST("hTpcNSigPtQA"), track.pt(), track.tpcNSigmaEl());

      histos.fill(HIST("hDcaTrack"), track.pt(), track.dcaXY());
    }
  }
  PROCESS_SWITCH(HfTaskSingleElectron, processData, "For real data", true);

  void processMc(soa::Filtered<MyCollisions>::iterator const& collision,
                 McTracksEl const& tracks,
                 aod::McParticles const&)
  {
    float const flagAnalysedEvt = 0.5;

    if (!collision.sel8()) {
      return;
    }

    if (collision.numContrib() < nContribMin) {
      return;
    }

    histos.fill(HIST("hVtxZ"), collision.posZ());
    histos.fill(HIST("hNEvents"), flagAnalysedEvt);

    for (const auto& track : tracks) {

      // ===== [pre-cut QA] BEGIN: fill track-quality distributions BEFORE trackSel (all charged tracks) =====
      histos.fill(HIST("hEtaTrack_noCut"), track.eta());
      histos.fill(HIST("hPtTrack_noCut"), track.pt());
      histos.fill(HIST("hNClsTpcTrack_noCut"), track.tpcNClsCrossedRows());
      histos.fill(HIST("hNClsFoundFindableTpcTrack_noCut"), track.tpcCrossedRowsOverFindableCls());
      histos.fill(HIST("hChi2TpcTrack_noCut"), track.tpcChi2NCl());
      histos.fill(HIST("hIbClsItsTrack_noCut"), track.itsNClsInnerBarrel());
      histos.fill(HIST("hChi2ItsTrack_noCut"), track.itsChi2NCl());
      histos.fill(HIST("hDcaXYTrack_noCut"), track.dcaXY());
      histos.fill(HIST("hDcaZTrack_noCut"), track.dcaZ());
      histos.fill(HIST("hTpcNSigPt_noCut"), track.pt(), track.tpcNSigmaEl());
      histos.fill(HIST("hTofNSigPt_noCut"), track.pt(), track.tofNSigmaEl());
      histos.fill(HIST("hTpcNSigEta_noCut"), track.eta(), track.tpcNSigmaEl());
      histos.fill(HIST("hNClsCrossedRowsEta_noCut"), track.eta(), track.tpcNClsCrossedRows());
      histos.fill(HIST("hDcaXYIbClsIts_noCut"), track.itsNClsInnerBarrel(), track.dcaXY());
      histos.fill(HIST("hTpcNSigNClsCrossedRows_noCut"), track.pt(), track.tpcNClsCrossedRows(), track.tpcNSigmaEl());
      histos.fill(HIST("hTofBetaP_noCut"), track.p(), track.beta());
      // MC-truth species before trackSel (TPC nSig from work/15 MCInformation_forEID; TOF nSig/beta from work/32 tof_quality_study)
      if (track.has_mcParticle()) {
        auto const mcParticleNoCut = track.mcParticle();
        int const absPdgNoCut = std::abs(mcParticleNoCut.pdgCode());
        if (absPdgNoCut == kElectron) {
          histos.fill(HIST("hTpcNSigPt_noCut_Ele"), track.pt(), track.tpcNSigmaEl());
          histos.fill(HIST("hTofNSigPt_noCut_Ele"), track.pt(), track.tofNSigmaEl());
          histos.fill(HIST("hTofBetaP_noCut_Ele"), track.p(), track.beta());
        } else if (absPdgNoCut == kPiPlus) {
          histos.fill(HIST("hTpcNSigPt_noCut_Pi"), track.pt(), track.tpcNSigmaEl());
          histos.fill(HIST("hTofNSigPt_noCut_Pi"), track.pt(), track.tofNSigmaEl());
          histos.fill(HIST("hTofBetaP_noCut_Pi"), track.p(), track.beta());
        } else if (absPdgNoCut == kKPlus) {
          histos.fill(HIST("hTpcNSigPt_noCut_K"), track.pt(), track.tpcNSigmaEl());
          histos.fill(HIST("hTofNSigPt_noCut_K"), track.pt(), track.tofNSigmaEl());
          histos.fill(HIST("hTofBetaP_noCut_K"), track.p(), track.beta());
        } else if (absPdgNoCut == kProton) {
          histos.fill(HIST("hTpcNSigPt_noCut_Pro"), track.pt(), track.tpcNSigmaEl());
          histos.fill(HIST("hTofNSigPt_noCut_Pro"), track.pt(), track.tofNSigmaEl());
          histos.fill(HIST("hTofBetaP_noCut_Pro"), track.p(), track.beta());
        } else {
          histos.fill(HIST("hTpcNSigPt_noCut_Other"), track.pt(), track.tpcNSigmaEl());
          histos.fill(HIST("hTofNSigPt_noCut_Other"), track.pt(), track.tofNSigmaEl());
          histos.fill(HIST("hTofBetaP_noCut_Other"), track.p(), track.beta());
        }
        int mpdgNoCut{};
        double mptNoCut{};
        int const srcNoCut = getElecSource(track, mptNoCut, mpdgNoCut);
        if (srcNoCut == DirectBeauty || srcNoCut == BeautyCharm) {
          histos.fill(HIST("hDcaZBeauty_noCut"), track.dcaZ());
        }
      }
      // ===== [pre-cut QA] END =====

      if (!trackSel(track)) {
        continue;
      }

      histos.fill(HIST("hEtaTrack"), track.eta());
      histos.fill(HIST("hPtTrack"), track.pt());

      histos.fill(HIST("hNClsTpcTrack"), track.tpcNClsCrossedRows());
      histos.fill(HIST("hNClsFoundFindableTpcTrack"), track.tpcCrossedRowsOverFindableCls());
      histos.fill(HIST("hChi2TpcTrack"), track.tpcChi2NCl());
      histos.fill(HIST("hIbClsItsTrack"), track.itsNClsInnerBarrel());
      histos.fill(HIST("hDcaXYTrack"), track.dcaXY());
      histos.fill(HIST("hDcaZTrack"), track.dcaZ());

      histos.fill(HIST("hTofNSigPt"), track.pt(), track.tofNSigmaEl());
      histos.fill(HIST("hTpcNSigPt"), track.pt(), track.tpcNSigmaEl());

      int mpdg{};   // electron source pdg code
      double mpt{}; // electron source pt
      int const source = getElecSource(track, mpt, mpdg);

      if (source == DirectBeauty || source == BeautyCharm) {
        histos.fill(HIST("hPdgB"), mpdg);
        histos.fill(HIST("hDcaBeauty"), track.pt(), track.dcaXY());
      }

      if (source == DirectCharm) {
        histos.fill(HIST("hPdgC"), mpdg);
        histos.fill(HIST("hDcaCharm"), track.pt(), track.dcaXY());
      }

      if (source >= GammaPi0 && source <= GammaSigma) {
        histos.fill(HIST("hPdgCo"), mpdg);
        histos.fill(HIST("hDcaConv"), track.pt(), track.dcaXY());
      }

      if (source >= Pi0 && source <= Sigma) {
        histos.fill(HIST("hPdgDa"), mpdg);
        histos.fill(HIST("hDcaDalitz"), track.pt(), track.dcaXY());
      }

      if (track.tpcNSigmaEl() < nSigmaTpcHadronMax && track.tpcNSigmaEl() > nSigmaTpcHadronMin)
        histos.fill(HIST("hDcaHadron"), track.pt(), track.dcaXY());

      if (std::abs(track.tofNSigmaEl()) > nSigmaTofMax) {
        continue;
      }
      histos.fill(HIST("hTofNSigPtQA"), track.pt(), track.tofNSigmaEl());
      histos.fill(HIST("hTpcNSigPtAfterTofCut"), track.pt(), track.tpcNSigmaEl());

      if (track.tpcNSigmaEl() < nSigmaTpcMin || track.tpcNSigmaEl() > nSigmaTpcMax) {
        continue;
      }
      histos.fill(HIST("hTpcNSigPtQA"), track.pt(), track.tpcNSigmaEl());

      histos.fill(HIST("hDcaTrack"), track.pt(), track.dcaXY());
    }
  }
  PROCESS_SWITCH(HfTaskSingleElectron, processMc, "For real data", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfTaskSingleElectron>(cfgc)};
}
