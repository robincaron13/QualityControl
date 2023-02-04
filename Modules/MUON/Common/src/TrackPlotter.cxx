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

#include "MUONCommon/TrackPlotter.h"

#include <DataFormatsITSMFT/ROFRecord.h>
#include <DataFormatsMFT/TrackMFT.h>
#include <DataFormatsMCH/TrackMCH.h>
#include <ReconstructionDataFormats/TrackMCHMID.h>
#include <ReconstructionDataFormats/GlobalFwdTrack.h>
#include <Framework/DataRefUtils.h>
#include <Framework/InputRecord.h>
#include <TH2F.h>
#include <gsl/span>
#include <set>

namespace
{

static void setXAxisLabels(TProfile* h)
{
  TAxis* axis = h->GetXaxis();
  for (int i = 1; i <= 10; i++) {
    auto label = fmt::format("CH{}", i);
    axis->SetBinLabel(i, label.c_str());
  }
}

} // namespace

using namespace o2::dataformats;

namespace o2::quality_control_modules::muon
{

TrackPlotter::TrackPlotter(int maxTracksPerTF, GID::Source source, std::string path) : mSrc(source), mPath(path)
{
  createTrackHistos(maxTracksPerTF);
  createTrackPairHistos();
}

void TrackPlotter::createTrackHistos(int maxTracksPerTF)
{
  mNofTracksPerTF[0] = createHisto<TH1F>(TString::Format("%sPositive/TracksPerTF", mPath.c_str()), "Number of tracks per TimeFrame (+);Number of tracks per TF", maxTracksPerTF, 0, maxTracksPerTF, true, "logy");
  mNofTracksPerTF[1] = createHisto<TH1F>(TString::Format("%sNegative/TracksPerTF", mPath.c_str()), "Number of tracks per TimeFrame (-);Number of tracks per TF", maxTracksPerTF, 0, maxTracksPerTF, true, "logy");
  mNofTracksPerTF[2] = createHisto<TH1F>(TString::Format("%sTracksPerTF", mPath.c_str()), "Number of tracks per TimeFrame;Number of tracks per TF", maxTracksPerTF, 0, maxTracksPerTF, true, "logy");

  mTrackDCA[0] = createHisto<TH1F>(TString::Format("%sPositive/TrackDCA", mPath.c_str()), "Track DCA (+);DCA (cm)", 500, 0, 500);
  mTrackDCA[1] = createHisto<TH1F>(TString::Format("%sNegative/TrackDCA", mPath.c_str()), "Track DCA (-);DCA (cm)", 500, 0, 500);
  mTrackDCA[2] = createHisto<TH1F>(TString::Format("%sTrackDCA", mPath.c_str()), "Track DCA;DCA (cm)", 500, 0, 500);

  mTrackPDCA[0] = createHisto<TH1F>(TString::Format("%sPositive/TrackPDCA", mPath.c_str()), "Track p#timesDCA (+);p#timesDCA (GeVcm/c)", 5000, 0, 5000);
  mTrackPDCA[1] = createHisto<TH1F>(TString::Format("%sNegative/TrackPDCA", mPath.c_str()), "Track p#timesDCA (-);p#timesDCA (GeVcm/c)", 5000, 0, 5000);
  mTrackPDCA[2] = createHisto<TH1F>(TString::Format("%sTrackPDCA", mPath.c_str()), "Track p#timesDCA;p#timesDCA (GeVcm/c)", 5000, 0, 5000);

  mTrackPt[0] = createHisto<TH1F>(TString::Format("%sPositive/TrackPt", mPath.c_str()), "Track p_{T} (+);p_{T} (GeV/c)", 300, 0, 30, false, "logy");
  mTrackPt[1] = createHisto<TH1F>(TString::Format("%sNegative/TrackPt", mPath.c_str()), "Track p_{T} (-);p_{T} (GeV/c)", 300, 0, 30, false, "logy");
  mTrackPt[2] = createHisto<TH1F>(TString::Format("%sTrackPt", mPath.c_str()), "Track p_{T};p_{T} (GeV/c)", 300, 0, 30, false, "logy");

  mTrackEta[0] = createHisto<TH1F>(TString::Format("%sPositive/TrackEta", mPath.c_str()), "Track #eta (+);#eta", 200, -4.5, -2);
  mTrackEta[1] = createHisto<TH1F>(TString::Format("%sNegative/TrackEta", mPath.c_str()), "Track #eta (-);#eta", 200, -4.5, -2);
  mTrackEta[2] = createHisto<TH1F>(TString::Format("%sTrackEta", mPath.c_str()), "Track #eta;#eta", 200, -4.5, -2);

  mTrackPhi[0] = createHisto<TH1F>(TString::Format("%sPositive/TrackPhi", mPath.c_str()), "Track #phi (+);#phi (deg)", 360, 0, 360);
  mTrackPhi[1] = createHisto<TH1F>(TString::Format("%sNegative/TrackPhi", mPath.c_str()), "Track #phi (-);#phi (deg)", 360, 0, 360);
  mTrackPhi[2] = createHisto<TH1F>(TString::Format("%sTrackPhi", mPath.c_str()), "Track #phi;#phi (deg)", 360, 0, 360);

  mTrackRAbs[0] = createHisto<TH1F>(TString::Format("%sPositive/TrackRAbs", mPath.c_str()), "Track R_{abs} (+);R_{abs} (cm)", 1000, 0, 100);
  mTrackRAbs[1] = createHisto<TH1F>(TString::Format("%sNegative/TrackRAbs", mPath.c_str()), "Track R_{abs} (-);R_{abs} (cm)", 1000, 0, 100);
  mTrackRAbs[2] = createHisto<TH1F>(TString::Format("%sTrackRAbs", mPath.c_str()), "Track R_{abs};R_{abs} (cm)", 1000, 0, 100);

  mTrackChi2OverNDF[0] = createHisto<TH1F>(TString::Format("%sPositive/TrackChi2OverNDF", mPath.c_str()), "Track #chi^{2}/ndf (+);#chi^{2}/ndf", 500, 0, 50);
  mTrackChi2OverNDF[1] = createHisto<TH1F>(TString::Format("%sNegative/TrackChi2OverNDF", mPath.c_str()), "Track #chi^{2}/ndf (-);#chi^{2}/ndf", 500, 0, 50);
  mTrackChi2OverNDF[2] = createHisto<TH1F>(TString::Format("%sTrackChi2OverNDF", mPath.c_str()), "Track #chi^{2}/ndf;#chi^{2}/ndf", 500, 0, 50);

  mTrackBC = createHisto<TH1F>(TString::Format("%sTrackBC", mPath.c_str()), "Track BC;BC", o2::constants::lhc::LHCMaxBunches, 0, o2::constants::lhc::LHCMaxBunches, true);
  mTrackBCWidth = createHisto<TH1F>(TString::Format("%sTrackBCWidth", mPath.c_str()), "Track BCWidth;BC Width", 400, 0, 400);

  mMatchScoreMFTMCH = createHisto<TH1F>(TString::Format("%sMatchScoreMFTMCH", mPath.c_str()), "Match Score MFT-MCH;score", 1000, 0, 100, true);
  mMatchChi2MFTMCH = createHisto<TH1F>(TString::Format("%sMatchChi2MFTMCH", mPath.c_str()), "Match #chi^{2} MFT-MCH;#chi^{2}", 1000, 0, 100, true);
  mMatchNMFTCandidates = createHisto<TH1F>(TString::Format("%sMatchNMFTCandidates", mPath.c_str()), "MFT Candidates;candidates", 1000, 0, 1000, true);

    std::string titleStr;
    float range = 800;
    if (mSrc == GID::MCHMID) {
      titleStr = "(MID-MCH)";
      range = 20;
    }
    if (mSrc == GID::MFTMCHMID) {
      titleStr = "(MFT-MID)";
    }
    if (mSrc == GID::MFTMCH) {
      titleStr = "(MFT-MCH)";
    }
    mTrackDT = std::make_unique<TH1F>(TString::Format("%sTrackDT", mPath.c_str()),
                                      TString::Format("Track Time Correlation %s;bc", titleStr.c_str()),
                                      1600, -range, range);
    histograms().emplace_back(HistInfo{ mTrackDT.get(), "hist", "logy" });

  mTrackPosAtMFT[0] = std::make_unique<TH2F>(TString::Format("%sTrackPosAtMFT_MFTtrack", mPath.c_str()), "Track position at MFT exit (MFT track);X (cm);Y (cm)", 80, -200, 200, 80, -200, 200);
  histograms().emplace_back(HistInfo{ mTrackPosAtMFT[0].get(), "col", "" });
  mTrackPosAtMFT[1] = std::make_unique<TH2F>(TString::Format("%sTrackPosAtMFT_MCHtrack", mPath.c_str()), "Track position at MFT exit (MCH track);X (cm);Y (cm)", 80, -200, 200, 80, -200, 200);
  histograms().emplace_back(HistInfo{ mTrackPosAtMFT[1].get(), "col", "" });

  mTrackDxMFT = createHisto<TH1F>(TString::Format("%sTrackDxMFT", mPath.c_str()), "Track Dx;cm", 200, -100, 100);
  mTrackDyMFT = createHisto<TH1F>(TString::Format("%sTrackDyMFT", mPath.c_str()), "Track Dy;cm", 200, -100, 100);
  mTrackSxMFT = createHisto<TH1F>(TString::Format("%sTrackSxMFT", mPath.c_str()), "Track Sx;", 400, -2, 2);
  mTrackSyMFT = createHisto<TH1F>(TString::Format("%sTrackSyMFT", mPath.c_str()), "Track Sy;", 400, -2, 2);
}

void TrackPlotter::createTrackPairHistos()
{
  mMinv = createHisto<TH1F>(TString::Format("%sMinv", mPath.c_str()), "#mu^{+}#mu^{-} invariant mass;M_{#mu^{+}#mu^{-}} (GeV/c^{2})", 800, 0, 16);
  mMinvBgd = createHisto<TH1F>(TString::Format("%sMinvBgd", mPath.c_str()), "#mu^{+}#mu^{-} inv. mass background;M_{#mu^{+}#mu^{-}} (GeV/c^{2})", 800, 0, 16);
  mDimuonDT = createHisto<TH1F>(TString::Format("%sDimuonDT", mPath.c_str()), "#mu^{+}#mu^{-} time diff.;ns", 4000, -2000, 2000);
}

void TrackPlotter::fillTrackPairHistos(gsl::span<const MuonTrack> tracks)
{
  if (tracks.size() > 1) {
    for (auto i = 0; i < tracks.size(); i++) {
      auto ti = tracks[i].getMuonMomentumAtVertex();
      if (!tracks[i].canBeMuon()) {
        continue;
      }
      for (auto j = i + 1; j < tracks.size(); j++) {
        if (!tracks[j].canBeMuon()) {
          continue;
        }
        if (tracks[i].getSign() == tracks[j].getSign()) {
          continue;
        }
        auto dt = (tracks[i].getIR() - tracks[j].getIR()).bc2ns();
        auto dtMUS = tracks[i].getTime().getTimeStamp() - tracks[j].getTime().getTimeStamp();
        mDimuonDT->Fill(dtMUS * 1000);
        // tracks are considered to be correlated if they are closer than 1 us in time
        if (std::abs(dtMUS) < 0.1) {
          auto tj = tracks[j].getMuonMomentumAtVertex();
          auto p = ti + tj;
          mMinv->Fill(p.M());
        }
        // the shape of the combinatorial background is derived by combining tracks
        // than belong to different orbits (more than 90 us apart)
        if (std::abs(dtMUS) > 90) {
          auto tj = tracks[j].getMuonMomentumAtVertex();
          auto p = ti + tj;
          mMinvBgd->Fill(p.M());
        }
      }
    }
  }
}

bool TrackPlotter::fillTrackHistos(const MuonTrack& track)
{
  int q = (track.getSign() < 0) ? 1 : 0;
  if (!track.canBeMuon()) {
    return false;
  }

  // only consider

  mTrackBC->Fill(track.getIR().bc);

  switch (mSrc) {
    case GID::MCHMID: {
      auto dtBC = track.getIRMID().toLong() - track.getIRMCH().toLong();
      auto dtMUS = track.getTimeMID().getTimeStamp() - track.getTimeMCH().getTimeStamp();
      mTrackDT->Fill(dtBC);
      //mTrackDT->Fill(dtMUS * 1000);
      break;
    }
    case GID::MFTMCH: {
      auto dtBC = track.getIRMFT().toLong() - track.getIRMCH().toLong();
      auto dtMUS = track.getTimeMFT().getTimeStamp() - track.getTimeMCH().getTimeStamp();
      mTrackDT->Fill(dtBC);
      //mTrackDT->Fill(dtMUS * 1000);
      mMatchScoreMFTMCH->Fill(track.getMatchInfoFwd().getMFTMCHMatchingScore());
      mMatchChi2MFTMCH->Fill(track.getMatchInfoFwd().getMFTMCHMatchingChi2());
      mMatchNMFTCandidates->Fill(track.getMatchInfoFwd().getNMFTCandidates());
      break;
    }
    case GID::MFTMCHMID: {
      auto dtBC1 = track.getIRMFT().toLong() - track.getIR().toLong();
      auto dtMUS1 = track.getTimeMFT().getTimeStamp() - track.getTimeMCH().getTimeStamp();
      mTrackDT->Fill(dtBC1);
      //mTrackDT->Fill(dtMUS1 * 1000);
      mMatchScoreMFTMCH->Fill(track.getMatchInfoFwd().getMFTMCHMatchingScore());
      mMatchChi2MFTMCH->Fill(track.getMatchInfoFwd().getMFTMCHMatchingChi2());
      mMatchNMFTCandidates->Fill(track.getMatchInfoFwd().getNMFTCandidates());
      break;
    }
    default:
      break;
  }

  double dca = track.getDCA();
  mTrackDCA[q]->Fill(dca);
  mTrackDCA[2]->Fill(dca);

  double pdca = track.getPDCAMCH();
  mTrackPDCA[q]->Fill(pdca);
  mTrackPDCA[2]->Fill(pdca);

  auto muon = track.getMuonMomentumAtVertex();

  mTrackEta[q]->Fill(muon.eta());
  mTrackEta[2]->Fill(muon.eta());
  mTrackPhi[q]->Fill(muon.phi() * TMath::RadToDeg() + 180);
  mTrackPhi[2]->Fill(muon.phi() * TMath::RadToDeg() + 180);
  mTrackPt[q]->Fill(muon.pt());
  mTrackPt[2]->Fill(muon.pt());

  auto rAbs = track.getRAbs();
  mTrackRAbs[q]->Fill(rAbs);
  mTrackRAbs[2]->Fill(rAbs);

  float zMFT = sLastMFTPlaneZ;
  o2::mch::TrackParam trackParamAtMFT;

  track.extrapToZMFT(trackParamAtMFT, zMFT);
  double xMFT = trackParamAtMFT.getNonBendingCoor();
  double yMFT = trackParamAtMFT.getBendingCoor();
  double sxMFT = TMath::ATan(trackParamAtMFT.px() / TMath::Abs(trackParamAtMFT.pz()));
  double syMFT = TMath::ATan(trackParamAtMFT.py() / TMath::Abs(trackParamAtMFT.pz()));
  mTrackPosAtMFT[0]->Fill(trackParamAtMFT.getNonBendingCoor(), trackParamAtMFT.getBendingCoor());

  track.extrapToZMCH(trackParamAtMFT, zMFT);
  double xMCH = trackParamAtMFT.getNonBendingCoor();
  double yMCH = trackParamAtMFT.getBendingCoor();
  double sxMCH = TMath::ATan(trackParamAtMFT.px() / TMath::Abs(trackParamAtMFT.pz()));
  double syMCH = TMath::ATan(trackParamAtMFT.py() / TMath::Abs(trackParamAtMFT.pz()));
  mTrackPosAtMFT[1]->Fill(trackParamAtMFT.getNonBendingCoor(), trackParamAtMFT.getBendingCoor());

  mTrackDxMFT->Fill(xMCH - xMFT);
  mTrackDyMFT->Fill(yMCH - yMFT);
  mTrackSxMFT->Fill(sxMCH - sxMFT);
  mTrackSyMFT->Fill(syMCH - syMFT);

  return true;
}

void TrackPlotter::fillHistograms(const o2::globaltracking::RecoContainer& recoCont)
{
  std::vector<MuonTrack> muonTracks;
  if (mSrc == GID::MCH) {
    auto tracksMCH = recoCont.getMCHTracks();
    for (auto& t : tracksMCH) {
      muonTracks.emplace_back(&t, recoCont, mFirstTForbit);
    }
  }
  if (mSrc == GID::MFTMCH || mSrc == GID::MFTMCHMID) {
    auto tracksFwd = recoCont.getGlobalFwdTracks();
    for (auto& t : tracksFwd) {
      MuonTrack mt(&t, recoCont, mFirstTForbit);
      // skip tracks without MID if full matching is requested
      if (mSrc == GID::MFTMCHMID && !mt.hasMID()) {
        continue;
      }
      muonTracks.emplace_back(&t, recoCont, mFirstTForbit);
    }
  }
  if (mSrc == GID::MCHMID) {
    auto tracksMCHMID = recoCont.getMCHMIDMatches();
    for (auto& t : tracksMCHMID) {
      muonTracks.emplace_back(&t, recoCont, mFirstTForbit);
    }
  }

  int nPos{ 0 };
  int nNeg{ 0 };

  for (auto& t : muonTracks) {
    if (t.getSign() < 0) {
      nNeg += 1;
    } else {
      nPos += 1;
    }
  }

  mNofTracksPerTF[0]->Fill(nPos);
  mNofTracksPerTF[1]->Fill(nNeg);
  mNofTracksPerTF[2]->Fill(muonTracks.size());

  decltype(muonTracks.size()) nok{ 0 };

  for (const auto& mt : muonTracks) {
    bool ok = fillTrackHistos(mt);
    if (ok) {
      ++nok;
    }
  }

  fillTrackPairHistos(muonTracks);
}

} // namespace o2::quality_control_modules::muon
