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

  if(mMode == 0) // fill Gain Correction Parameters for ZDC-SMD
  {
    mOutPut_GainCorr = Form("./file_%s_GainCorr_%s.root",recoEP::mBeamEnergy[energy].c_str(),jobId.c_str());
  }
  if(mMode == 1) // fill Re-Center Correction Parameters for ZDC-SMD East/West & TPC East/West/Full
  {
    mOutPut_ReCenterPar = Form("./file_%s_ReCenterParameter_%s.root",recoEP::mBeamEnergy[energy].c_str(),jobId.c_str());
  }
  if(mMode == 2) // fill Shift Correction Parameters for ZDC-SMD East/West & TPC East/West/Full
  {
    mOutPut_ShiftPar = Form("./file_%s_ShiftParameter_%s.root",recoEP::mBeamEnergy[energy].c_str(),jobId.c_str());
  }
  if(mMode == 3) // fill Shift Correction Parameters for ZDC-SMD Full EP
  {
    mOutPut_ShiftParFull = Form("./file_%s_ShiftParameterFull_%s.root",recoEP::mBeamEnergy[energy].c_str(),jobId.c_str());
  }
  if(mMode == 4) // fill Event Plane Resolution for ZDC-SMD East/West & TPC East/West
  {
    mOutPut_Resolution = Form("./file_%s_Resolution_%s.root",recoEP::mBeamEnergy[energy].c_str(),jobId.c_str());
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
  { // fill ReCenter Correction Parameters for ZDC-SMD East/West & TPC East/West/Full
    mFile_ReCenterPar = new TFile(mOutPut_ReCenterPar.c_str(),"RECREATE");

    mEventPlaneProManager->initZdcReCenter(); // ZDC-SMD
    mEventPlaneHistoManager->initZdcRawEP();
    mZdcEpManager->readGainCorr(); // read in Gain Correction Parameters

    mEventPlaneProManager->initTpcReCenter(); // TPC
    mEventPlaneHistoManager->initTpcRawEP();
  }
  if(mMode == 2)
  { // fill Shift Correction Parameters for ZDC-SMD East/West & TPC East/West/Full
    mFile_ShiftPar = new TFile(mOutPut_ShiftPar.c_str(),"RECREATE");

    mEventPlaneProManager->initZdcShift(); // ZDC-SMD
    mEventPlaneHistoManager->initZdcReCenterEP();
    mZdcEpManager->readGainCorr(); // read in Gain Correction Parameters
    mZdcEpManager->readReCenterCorr(); // read in ReCenter Correction Parameters

    mEventPlaneProManager->initTpcShift(); // TPC
    mEventPlaneHistoManager->initTpcReCenterEP();
    mTpcEpManager->readReCenterCorr(); // read in ReCenter Correction Parameters
  }
  if(mMode == 3)
  { // fill Shift Correction Parameters for ZDC-SMD Full EP
    mFile_ShiftParFull = new TFile(mOutPut_ShiftParFull.c_str(),"RECREATE");

    mEventPlaneProManager->initZdcShiftFull(); // ZDC-SMD Full EP
    mZdcEpManager->readGainCorr(); // read in Gain Correction Parameters
    mZdcEpManager->readReCenterCorr(); // read in ReCenter Correction Parameters
    mZdcEpManager->readShiftCorr(); // read in Shift Correction Parameters
  }
  if(mMode == 4)
  { // fill EP Resolution for ZDC-SMD Sub & TPC Sub/Ran
    mFile_Resolution = new TFile(mOutPut_Resolution.c_str(),"RECREATE");

    mEventPlaneProManager->initZdcResolution(); // ZDC-SMD
    mEventPlaneHistoManager->initZdcShiftEP();
    mZdcEpManager->readGainCorr(); // read in Gain Correction Parameters
    mZdcEpManager->readReCenterCorr(); // read in ReCenter Correction Parameters
    mZdcEpManager->readShiftCorr(); // read in Shift Correction Parameters
    mZdcEpManager->readShiftCorrFull(); // read in Full Shift Correction Parameters

    mEventPlaneProManager->initTpcResolution(); // TPC
    mEventPlaneHistoManager->initTpcShiftEP();
    mTpcEpManager->readReCenterCorr(); // read in ReCenter Correction Parameters
    mTpcEpManager->readShiftCorr(); // read in ReCenter Correction Parameters
    mUsedTrackCounter = 0; // for Random EP
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
  if(mMode == 2)
  {
    if(mOutPut_ShiftPar != "")
    {
      mFile_ShiftPar->cd();

      mEventPlaneProManager->writeZdcShift(); // ZDC-SMD
      mEventPlaneHistoManager->writeZdcReCenterEP();

      mEventPlaneProManager->writeTpcShift(); // TPC
      mEventPlaneHistoManager->writeTpcReCenterEP();

      mFile_ShiftPar->Close();
    }
  }
  if(mMode == 3)
  {
    if(mOutPut_ShiftParFull != "")
    {
      mFile_ShiftParFull->cd();

      mEventPlaneProManager->writeZdcShiftFull(); // ZDC-SMD 

      mFile_ShiftParFull->Close();
    }
  }
  if(mMode == 4)
  {
    if(mOutPut_Resolution != "")
    {
      mFile_Resolution->cd();

      mEventPlaneProManager->writeZdcResolution(); // ZDC-SMD
      mEventPlaneHistoManager->writeZdcShiftEP();

      mEventPlaneProManager->writeTpcResolution(); // TPC
      mEventPlaneHistoManager->writeTpcShiftEP();

      mFile_Resolution->Close();
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
      { // fill ReCenter Correction for ZDC-SMD East/West & TPC East/West/Full
	// ZDC-SMD: 
	// apply gain correction 
	// fill recenter correction parameter & fill raw ZDC-SMD EP
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
	// fill recenter correction parameter & fill raw TPC EP
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

      if(mMode == 2)
      { // fill Shift Correction for ZDC-SMD East/West & TPC East/West/Full
	// ZDC-SMD: 
	// apply gain correction & apply recenter correction to East/West
	// fill shift correction paramter East/West & fill recenter ZDC-SMD EP
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
	  mEventPlaneProManager->fillZdcShiftEast(vZdcQ1East,cent9,runIndex,vzSign);
	  mEventPlaneProManager->fillZdcShiftWest(vZdcQ1West,cent9,runIndex,vzSign);
	  mEventPlaneHistoManager->fillZdcReCenterSubEP(vZdcQ1East,vZdcQ1West,cent9,runIndex);
	  mEventPlaneHistoManager->fillZdcReCenterFullEP(vZdcQ1Full,cent9,runIndex);
	}

	// TPC: 
	// apply recenter correction to East/West/Full
	// fill shift correction parameter East/West/Full & fill recenter TPC EP
	for(unsigned int i_track = 0; i_track < nTracks; ++i_track)
	{ // calculate QVector after recenter correction
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
	      mTpcEpManager->addTrackEast(picoTrack);
	    }
	    if(mTpcEpManager->passTrackEtaWest(picoTrack)) // West sub EP 
	    {
	      mTpcEpManager->addTrackWest(picoTrack);
	    }
	    if(mTpcEpManager->passTrackEtaFull(picoTrack)) // Full EP 
	    {
	      mTpcEpManager->addTrackFull(picoTrack);
	    }
	  }
	}
	TVector2 vTpcQ2East = mTpcEpManager->getQVector(0); // receter QVec East
	TVector2 vTpcQ2West = mTpcEpManager->getQVector(1); // receter QVec West
	TVector2 vTpcQ2Full = mTpcEpManager->getQVector(2); // receter QVec Full

	if(mTpcEpManager->passTrackEtaNumCut())
	{ // fill shift correction for east/west sub EP
	  if( !(vTpcQ2East.Mod() < 1e-10 || vTpcQ2West.Mod() < 1e-10) )
	  {
	    mEventPlaneProManager->fillTpcShiftEast(vTpcQ2East,cent9,runIndex,vzSign);
	    mEventPlaneProManager->fillTpcShiftWest(vTpcQ2West,cent9,runIndex,vzSign);
	    mEventPlaneHistoManager->fillTpcReCenterSubEP(vTpcQ2East,vTpcQ2West,cent9,runIndex);
	  }
	}
	if(mTpcEpManager->passTrackFullNumCut())
	{ // fill shift correction for full EP
	  if( !(vTpcQ2Full.Mod() < 1e-10) )
	  {
	    mEventPlaneProManager->fillTpcShiftFull(vTpcQ2Full,cent9,runIndex,vzSign);
	    mEventPlaneHistoManager->fillTpcReCenterFullEP(vTpcQ2Full,cent9,runIndex);
	  }
	}
      }
      if(mMode == 3)
      { // fill Shift Correction for ZDC-SMD Full EP
	// ZDC-SMD: 
	// apply gain correction & apply recenter correction East/West & apply shift correction East/West
	// fill shift correction paramter Full
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
	  mEventPlaneProManager->fillZdcShiftFull(vZdcQ1Full,cent9,runIndex,vzSign);
	}
      }
      if(mMode == 4)
      { // fill EP Resolution for ZDC-SMD Sub & TPC Sub/Ran
	// ZDC-SMD: 
	// apply gain correction & apply recenter correction & apply Shift to East/West/Full
	// fill EP Resolution for ZDC-SMD Sub & fill shift ZDC-SMD EP
	for(int i_slat = 0; i_slat < 8; ++i_slat) // read in gain correction factors for ZDC-SMD
	{
	  mZdcEpManager->setZdcSmdGainCorr(0,0,i_slat,mPicoEvent->ZdcSmdEastVertical(i_slat));
	  mZdcEpManager->setZdcSmdGainCorr(0,1,i_slat,mPicoEvent->ZdcSmdEastHorizontal(i_slat));
	  mZdcEpManager->setZdcSmdGainCorr(1,0,i_slat,mPicoEvent->ZdcSmdWestVertical(i_slat));
	  mZdcEpManager->setZdcSmdGainCorr(1,1,i_slat,mPicoEvent->ZdcSmdWestHorizontal(i_slat));
	}

	TVector2 vZdcQ1East = mZdcEpManager->getQEast(mMode);
	TVector2 vZdcQ1West = mZdcEpManager->getQWest(mMode);
	TVector2 vZdcQ1Diff = vZdcQ1West-vZdcQ1East;
	TVector2 vZdcQ1Full = mZdcEpManager->getQFull(vZdcQ1East,vZdcQ1West);
	if( !(vZdcQ1East.Mod() < 1e-10 || vZdcQ1West.Mod() < 1e-10 || vZdcQ1Full.Mod() < 1e-10) )
	{
	  mEventPlaneProManager->fillZdcResSub(vZdcQ1East,vZdcQ1West,cent9,runIndex);
	  mEventPlaneHistoManager->fillZdcShiftSubEP(vZdcQ1East,vZdcQ1West,cent9,runIndex);
	  mEventPlaneHistoManager->fillZdcShiftFullEP(vZdcQ1Diff,vZdcQ1Full,cent9,runIndex);
	}

	// TPC: 
	// apply recenter correction & apply shift correction to East/West/Full
	// fill EP Resolution for Sub/Ran & fill Shift TPC EP
	for(unsigned int i_track = 0; i_track < nTracks; ++i_track)
	{ // calculate QVector after recenter correction
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
	      mTpcEpManager->addTrackEast(picoTrack);
	    }
	    if(mTpcEpManager->passTrackEtaWest(picoTrack)) // West sub EP 
	    {
	      mTpcEpManager->addTrackWest(picoTrack);
	    }
	    if(mTpcEpManager->passTrackEtaFull(picoTrack)) // Full EP 
	    {
	      mTpcEpManager->addTrackFull(picoTrack);
	      mUsedTrackCounter++;
	    }
	  }
	}

	// Random EP
	int ranTrack[mUsedTrackCounter];
	double ranCounter = (double)mUsedTrackCounter/2.0 - 1.0;
	for(int i_track = 0; i_track < mUsedTrackCounter; ++i_track)
	{
	  ranTrack[i_track] = i_track;
	}
	std::srand(time(0));
	std::random_shuffle(ranTrack,ranTrack+mUsedTrackCounter);
	mUsedTrackCounter = 0;
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
	      if((double)ranTrack[mUsedTrackCounter] > ranCounter) // RanA
	      {
		mTpcEpManager->addTrackRanA(picoTrack);
	      }
	      else // Sub Event B
	      {
		mTpcEpManager->addTrackRanB(picoTrack);
	      }
	      mUsedTrackCounter++;
	    }
	  }
	}
	mUsedTrackCounter = 0;

	TVector2 vTpcQ2East = mTpcEpManager->getQVector(0); // receter QVec East
	TVector2 vTpcQ2West = mTpcEpManager->getQVector(1); // receter QVec West
	TVector2 vTpcQ2Full = mTpcEpManager->getQVector(2); // receter QVec Full
	// TVector2 vTpcQ2RanA = mTpcEpManager->getQVector(3); // receter QVec RanA
	// TVector2 vTpcQ2RanB = mTpcEpManager->getQVector(4); // receter QVec RanB 

	if(mTpcEpManager->passTrackEtaNumCut())
	{ // fill shift correction for east/west sub EP
	  if( !(vTpcQ2East.Mod() < 1e-10 || vTpcQ2West.Mod() < 1e-10) )
	  {
	    double TpcPsiEast = mTpcEpManager->calShiftAngle2East();
	    double TpcPsiWest = mTpcEpManager->calShiftAngle2East();
	    mEventPlaneProManager->fillTpcResSub(TpcPsiEast,TpcPsiWest,cent9,runIndex);
	    mEventPlaneHistoManager->fillTpcShiftSubEP(TpcPsiEast,TpcPsiWest,cent9,runIndex);
	  }
	}
	if(mTpcEpManager->passTrackFullNumCut())
	{ // fill shift correction for full EP
	  if( !(vTpcQ2Full.Mod() < 1e-10) )
	  {
	    double TpcPsiRanA = mTpcEpManager->calShiftAngle2RanA();
	    double TpcPsiRanB = mTpcEpManager->calShiftAngle2RanB();
	    double TpcPsiFull = mTpcEpManager->calShiftAngle2Full();
	    mEventPlaneProManager->fillTpcResRan(TpcPsiRanA,TpcPsiRanB,cent9,runIndex);
	    mEventPlaneHistoManager->fillTpcShiftRanEP(TpcPsiRanA,TpcPsiRanB,cent9,runIndex);
	    mEventPlaneHistoManager->fillTpcShiftFullEP(TpcPsiFull,cent9,runIndex);
	  }
	}
      }

      mZdcEpManager->clearZdcEp();
      mTpcEpManager->clearTpcEp();
    }
  }

  return kStOK;
}

