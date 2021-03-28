#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StRefMultCorr/StRefMultCorr.h"
#include "StRefMultCorr/CentralityMaker.h"
#include "StThreeVectorF.hh"
#include "StMessMgr.h"

#include "StRoot/StEventPlaneMaker/StEventPlaneMaker.h"
#include "StRoot/StEventPlaneMaker/StEventPlaneCons.h"
#include "StRoot/StEventPlaneMaker/StEventPlaneUtility.h"
#include "StRoot/StEventPlaneMaker/StEventPlaneCut.h"
#include "StRoot/StEventPlaneMaker/StEventPlaneHistoManager.h"
#include "StRoot/StEventPlaneMaker/StEventPlaneProManager.h"
#include "StRoot/StEventPlaneMaker/StZdcEpManager.h"
#include "StRoot/StEventPlaneMaker/StTpcEpManager.h"

#include <algorithm>

#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <TVector3.h>
#include <TMath.h>

ClassImp(StEventPlaneMaker)

//-----------------------------------------------------------------------------
StEventPlaneMaker::StEventPlaneMaker(const char* name, StPicoDstMaker *picoMaker, const string jobId, const int Mode, const int energy) : StMaker(name)
{
  mPicoDstMaker = picoMaker;
  mPicoDst = NULL;
  mRefMultCorr = NULL;
  mMode = Mode;
  mEnergy = energy;

  if(mMode == 0) // Get Gain Correction Parameters for ZDC-SMD
  {
    mOutPut_GainCorr = Form("./file_%s_GainCorr_%s.root",recoEP::mBeamEnergy[energy].c_str(),jobId.c_str());
  }
  if(mMode == 1) // Get Re-Center Correction Parameters for ZDC-SMD & TPC
  {
    mOutPut_ReCenterPar = Form("./file_%s_ReCenterParameter_%s.root",recoEP::mBeamEnergy[energy].c_str(),jobId.c_str());
  }
}

//----------------------------------------------------------------------------- 
StEventPlaneMaker::~StEventPlaneMaker()
{ /*  */ }

//----------------------------------------------------------------------------- 
int StEventPlaneMaker::Init() 
{
  mEventPlaneCut = new StEventPlaneCut(mEnergy);
  mEventPlaneHistoManager = new StEventPlaneHistoManager();
  mEventPlaneUtility = new StEventPlaneUtility(mEnergy);
  mEventPlaneUtility->initRunIndex(); // initialize std::map for run index
  mEventPlaneProManager = new StEventPlaneProManager();
  mZdcEpManager = new StZdcEpManager(mEnergy); // initialize ZDC EP Manager
  mTpcEpManager = new StTpcEpManager(mEnergy); // initialize TPC EP Manager

  if(!mRefMultCorr)
  {
    if(!mEventPlaneCut->isBES()) mRefMultCorr = CentralityMaker::instance()->getgRefMultCorr_Run14_AuAu200_VpdMB5_P16id(); // 200GeV_2014
    if(mEventPlaneCut->isBES()) mRefMultCorr = CentralityMaker::instance()->getRefMultCorr(); // BESII
  }

  if(mMode == 0)
  { // fill Gain Correction Factors for ZDC-SMD
    mFile_GainCorr = new TFile(mOutPut_GainCorr.c_str(),"RECREATE");
    mEventPlaneHistoManager->initZdcGainCorr();
  }
  if(mMode == 1)
  { // fill ReCenter Correction Parameters for ZDC-SMD & TPC
    mFile_ReCenterPar = new TFile(mOutPut_ReCenterPar.c_str(),"RECREATE");

    mEventPlaneProManager->initZdcReCenter(); // ZDC-SMD
    mEventPlaneHistoManager->initZdcRawEP();
    mZdcEpManager->readGainCorr();

    mEventPlaneProManager->initTpcReCenter(); // TPC
    mEventPlaneHistoManager->initTpcRawEP();
  }

  return kStOK;
}

//----------------------------------------------------------------------------- 
int StEventPlaneMaker::Finish() 
{
  if(mMode == 0)
  {
    if(mOutPut_GainCorr != "")
    {
      mFile_GainCorr->cd();
      mEventPlaneHistoManager->writeZdcGainCorr();
      mFile_GainCorr->Close();
    }
  }
  if(mMode == 1)
  {
    if(mOutPut_ReCenterPar != "")
    {
      mFile_ReCenterPar->cd();

      mEventPlaneProManager->writeZdcReCenter(); // ZDC-SMD
      mEventPlaneHistoManager->writeZdcRawEP();

      mEventPlaneProManager->writeTpcReCenter(); // TPC
      mEventPlaneHistoManager->writeTpcRawEP();

      mFile_ReCenterPar->Close();
    }
  }

  return kStOK;
}

//----------------------------------------------------------------------------- 
void StEventPlaneMaker::Clear(Option_t *opt) {
}

//----------------------------------------------------------------------------- 
int StEventPlaneMaker::Make() 
{
  if(!mPicoDstMaker) 
  {
    LOG_ERROR << " No PicoDstMaker! Skip! " << endm;
    return kStErr;
  }

  mPicoDst = mPicoDstMaker->picoDst();
  if(!mPicoDst) 
  {
    LOG_ERROR << " No PicoDst! Skip! " << endm;
    return kStErr;
  }

  mPicoEvent = (StPicoEvent*)mPicoDst->event();
  if(!mPicoEvent)
  {
    LOG_ERROR << "Error opening picoDst Event, skip!" << endm;
    return kStErr;
  }

  // MinBias trigger
  if( mEventPlaneCut->isMinBias(mPicoEvent) )
  {
    // Event Information
    const int runId = mPicoEvent->runId();
    const int eventId = mPicoEvent->eventId();
    const int refMult = mPicoEvent->refMult();
    const int grefMult = mPicoEvent->grefMult();
    const float vx    = mPicoEvent->primaryVertex().x(); // x works for both TVector3 and StThreeVectorF
    const float vy    = mPicoEvent->primaryVertex().y();
    const float vz    = mPicoEvent->primaryVertex().z();
    const float vzVpd = mPicoEvent->vzVpd();
    const float zdcX  = mPicoEvent->ZDCx();
    const unsigned int numOfBTofHits = mPicoDst->numberOfBTofHits(); // get number of tof hits
    // const unsigned short numOfBTofHits = mPicoEvent->btofTrayMultiplicity();
    const unsigned short numOfBTofMatch = mPicoEvent->nBTOFMatch(); // get number of tof match points
    const unsigned int nTracks = mPicoDst->numberOfTracks(); // get number of tracks

    // StRefMultCorr Cut & centrality
    if(!mRefMultCorr)
    {
      LOG_WARN << " No mRefMultCorr! Skip! " << endl;
      return kStErr;
    }

    mRefMultCorr->init(runId);
    if(!mEventPlaneCut->isBES()) mRefMultCorr->initEvent(grefMult,vz,zdcX); // 200GeV_2014
    if(mEventPlaneCut->isBES()) mRefMultCorr->initEvent(refMult,vz,zdcX); // BES-II might need Luminosity corrections

    // if(mRefMultCorr->isBadRun(runId))
    // {
    //   LOG_ERROR << "Bad Run from StRefMultCorr! Skip!" << endm;
    //   return kStErr;
    // }

    // vz sign
    int vzSign = 0; // 0 for -vz || 1 for vz
    vz > 0.0 ? vzSign = 1 : vzSign = 0;

    const int cent9 = mRefMultCorr->getCentralityBin9(); // get Centrality9
    const double reweight = mRefMultCorr->getWeight(); // get weight
    const int runIndex = mEventPlaneUtility->findRunIndex(runId); // find run index for a specific run
    // const int triggerBin = mEventPlaneCut->getTriggerBin(mPicoEvent);
    // cout << "runId = " << runId << ", runIndex = " << runIndex << endl;
    if(runIndex < 0)
    {
      LOG_ERROR << "Could not find this run Index from StEventPlaneUtility! Skip!" << endm;
      return kStErr;
    }

    bool isPileUpEventStEventPlaneCut = mEventPlaneCut->isPileUpEvent(grefMult,numOfBTofMatch,numOfBTofHits); // 200GeV
    if(mEventPlaneCut->isBES()) isPileUpEventStEventPlaneCut = mEventPlaneCut->isPileUpEvent(refMult,numOfBTofMatch,numOfBTofHits); // 54 GeV | always return false for 27 GeV
    bool isPileUpEventStRefMultCorr = !mRefMultCorr->passnTofMatchRefmultCut(1.0*refMult, 1.0*numOfBTofMatch); // 27 GeV | always return !true for other energies
    bool isPileUpEvent = isPileUpEventStEventPlaneCut || isPileUpEventStRefMultCorr;
    // cout << "isPileUpEvent = " << isPileUpEvent << ", isPileUpEventStEventPlaneCut = " << isPileUpEventStEventPlaneCut << ", isPileUpEventStRefMultCorr = " << isPileUpEventStRefMultCorr << endl;

    if(mEventPlaneCut->passEventCut(mPicoDst) && !isPileUpEvent && cent9 > -0.5)
    { // apply Event Cuts for anlaysis 
      mZdcEpManager->initZdcEp(cent9,runIndex,vzSign); // ZDC-SMD EP
      mTpcEpManager->initTpcEp(cent9,runIndex,vzSign); // TPC EP
      if(mMode == 0)
      { // fill Gain Correction Factors for BBC & ZDC
	for(int i_slat = 0; i_slat < 8; ++i_slat) // read in raw ADC value from ZDC-SMD
	{
	  mZdcEpManager->setZdcSmd(0,0,i_slat,mPicoEvent->ZdcSmdEastVertical(i_slat));
	  mZdcEpManager->setZdcSmd(0,1,i_slat,mPicoEvent->ZdcSmdEastHorizontal(i_slat));
	  mZdcEpManager->setZdcSmd(1,0,i_slat,mPicoEvent->ZdcSmdWestVertical(i_slat));
	  mZdcEpManager->setZdcSmd(1,1,i_slat,mPicoEvent->ZdcSmdWestHorizontal(i_slat));
	}
	for(int i_eastwest = 0; i_eastwest < 2; ++i_eastwest) // fill ZDC Gain Correction Histograms
	{
	  for(int i_verthori = 0; i_verthori < 2; ++i_verthori)
	  {
	    for(int i_slat = 0; i_slat < 8; ++i_slat)
	    {
	      mEventPlaneHistoManager->fillZdcGainCorr(i_eastwest,i_verthori,i_slat,runIndex,mZdcEpManager->getZdcSmd(i_eastwest,i_verthori,i_slat));
	      // cout << "i_eastwest = " << i_eastwest << ", i_verthori = " << i_verthori << ", i_slat = " << i_slat << ", zdc = " << mZdcEpManager->getZdcSmd(i_eastwest,i_verthori,i_slat) << endl;
	    }
	  }
	}
      }

      if(mMode == 1)
      { // fill ReCenter Correction for ZDC-SMD & TPC
	// ZDC-SMD: 
	// apply gain correction 
	// fill recenter correction parameter
	// fill raw ZDC-SMD EP
	for(int i_slat = 0; i_slat < 8; ++i_slat) // read in gain correction factors for ZDC-SMD
	{
	  mZdcEpManager->setZdcSmdGainCorr(0,0,i_slat,mPicoEvent->ZdcSmdEastVertical(i_slat));
	  mZdcEpManager->setZdcSmdGainCorr(0,1,i_slat,mPicoEvent->ZdcSmdEastHorizontal(i_slat));
	  mZdcEpManager->setZdcSmdGainCorr(1,0,i_slat,mPicoEvent->ZdcSmdWestVertical(i_slat));
	  mZdcEpManager->setZdcSmdGainCorr(1,1,i_slat,mPicoEvent->ZdcSmdWestHorizontal(i_slat));
	}

	TVector2 vZdcQ1East = mZdcEpManager->getQEast(mMode);
	TVector2 vZdcQ1West = mZdcEpManager->getQWest(mMode);
	TVector2 vZdcQ1Full = vZdcQ1West-vZdcQ1East;
	if( !(vZdcQ1East.Mod() < 1e-10 || vZdcQ1West.Mod() < 1e-10 || vZdcQ1Full.Mod() < 1e-10) )
	{
	  mEventPlaneProManager->fillZdcReCenterEast(vZdcQ1East,cent9,runIndex,vzSign);
	  mEventPlaneProManager->fillZdcReCenterWest(vZdcQ1West,cent9,runIndex,vzSign);
	  mEventPlaneHistoManager->fillZdcRawSubEP(vZdcQ1East,vZdcQ1West,cent9,runIndex);
	  mEventPlaneHistoManager->fillZdcRawFullEP(vZdcQ1Full,cent9,runIndex);
	}

	// TPC: 
	// fill recenter correction parameter
	// fill raw TPC EP
	for(unsigned int i_track = 0; i_track < nTracks; ++i_track)
	{ // calculate number of tracks used in EP reconstruction
	  StPicoTrack *picoTrack = (StPicoTrack*)mPicoDst->track(i_track); // get picoTrack
	  if(mEventPlaneCut->passTrackEP(picoTrack,mPicoEvent))
	  { // track cut for EP reconstruction
	    TVector3 primMom; // temp fix for StThreeVectorF & TVector3
	    const double primPx    = picoTrack->pMom().x(); // x works for both TVector3 and StThreeVectorF
	    const double primPy    = picoTrack->pMom().y();
	    const double primPz    = picoTrack->pMom().z();
	    primMom.SetXYZ(primPx,primPy,primPz);
	    const double primPt = primMom.Perp(); // track pT
	    if(mTpcEpManager->passTrackEtaEast(picoTrack)) // East sub EP 
	    {
	      mTpcEpManager->addTrackEastRaw(picoTrack);
	    }
	    if(mTpcEpManager->passTrackEtaWest(picoTrack)) // West sub EP 
	    {
	      mTpcEpManager->addTrackWestRaw(picoTrack);
	    }
	    if(mTpcEpManager->passTrackEtaFull(picoTrack)) // Full EP 
	    {
	      mTpcEpManager->addTrackFullRaw(picoTrack);
	    }
	  }
	}
	TVector2 vTpcQ2East = mTpcEpManager->getQVectorRaw(0); // raw QVec East
	TVector2 vTpcQ2West = mTpcEpManager->getQVectorRaw(1); // raw QVec West
	TVector2 vTpcQ2Full = mTpcEpManager->getQVectorRaw(2); // raw QVec Full

	if(mTpcEpManager->passTrackEtaNumRawCut())
	{ // fill recenter correction for east/west sub EP
	  for(unsigned int i_track = 0; i_track < nTracks; ++i_track)
	  {
	    StPicoTrack *picoTrack = (StPicoTrack*)mPicoDst->track(i_track); // get picoTrack
	    if(mEventPlaneCut->passTrackEP(picoTrack,mPicoEvent))
	    { // track cut for EP reconstruction
	      TVector3 primMom; // temp fix for StThreeVectorF & TVector3
	      const double primPx    = picoTrack->pMom().x(); // x works for both TVector3 and StThreeVectorF
	      const double primPy    = picoTrack->pMom().y();
	      const double primPz    = picoTrack->pMom().z();
	      primMom.SetXYZ(primPx,primPy,primPz);
	      const double primPt = primMom.Perp(); // track pT
	      if(mTpcEpManager->passTrackEtaEast(picoTrack)) // East sub EP 
	      {
		TVector2 vTpcq2East = mTpcEpManager->calq2Vector(picoTrack);
		mEventPlaneProManager->fillTpcReCenterEast(vTpcq2East,cent9,runIndex,vzSign,primPt); // fill recenter
	      }
	      if(mTpcEpManager->passTrackEtaWest(picoTrack)) // West sub EP 
	      {
		TVector2 vTpcq2West = mTpcEpManager->calq2Vector(picoTrack);
		mEventPlaneProManager->fillTpcReCenterWest(vTpcq2West,cent9,runIndex,vzSign,primPt); // fill recenter
	      }
	    }
	  }
	  if( !(vTpcQ2East.Mod() < 1e-10 || vTpcQ2West.Mod() < 1e-10) )
	  {
	    mEventPlaneHistoManager->fillTpcRawSubEP(vTpcQ2East,vTpcQ2West,cent9,runIndex);
	  }
	}
	if(mTpcEpManager->passTrackFullNumRawCut())
	{ // fill recenter correction for full EP
	  for(unsigned int i_track = 0; i_track < nTracks; ++i_track)
	  {
	    StPicoTrack *picoTrack = (StPicoTrack*)mPicoDst->track(i_track); // get picoTrack
	    if(mEventPlaneCut->passTrackEP(picoTrack,mPicoEvent))
	    { // track cut for EP reconstruction
	      TVector3 primMom; // temp fix for StThreeVectorF & TVector3
	      const double primPx    = picoTrack->pMom().x(); // x works for both TVector3 and StThreeVectorF
	      const double primPy    = picoTrack->pMom().y();
	      const double primPz    = picoTrack->pMom().z();
	      primMom.SetXYZ(primPx,primPy,primPz);
	      const double primPt = primMom.Perp(); // track pT
	      if(mTpcEpManager->passTrackEtaFull(picoTrack)) // Full EP 
	      {
		TVector2 vTpcq2Full = mTpcEpManager->calq2Vector(picoTrack);
		mEventPlaneProManager->fillTpcReCenterFull(vTpcq2Full,cent9,runIndex,vzSign,primPt); // fill recenter
	      }
	    }
	  }
	  if( !(vTpcQ2Full.Mod() < 1e-10) )
	  {
	    mEventPlaneHistoManager->fillTpcRawFullEP(vTpcQ2Full,cent9,runIndex);
	  }
	}
      }

      mZdcEpManager->clearZdcEp();
      mTpcEpManager->clearTpcEp();
    }
  }

  return kStOK;
}

