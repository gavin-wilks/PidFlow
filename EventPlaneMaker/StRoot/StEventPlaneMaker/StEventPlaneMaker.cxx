#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StRefMultCorr/StRefMultCorr.h"
#include "StRefMultCorr/CentralityMaker.h"
#include "StEpdUtil/StEpdEpFinder.h"
#include "StEpdUtil/StEpdEpInfo.h"
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
//#include <filesystem>

#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <TVector3.h>
#include <TMath.h>

ClassImp(StEventPlaneMaker)

//-----------------------------------------------------------------------------
StEventPlaneMaker::StEventPlaneMaker(const char* name, StPicoDstMaker *picoMaker, const string jobId, const char* inputDir, const int Mode, const int inputEpdMode, const int EpdMode, const int energy, float mipThresh = 0.3, float maxTile = 2.0) : StMaker(name)
{
  mPicoDstMaker = picoMaker;
  mPicoDst = NULL;
  mRefMultCorr = NULL;
  mEpFinder = NULL;
  mMode = Mode;
  mEpdMode = EpdMode;  
  mInputEpdMode = inputEpdMode;
  mEnergy = energy;
  mMipThresh = mipThresh;
  mMaxTile = maxTile;

  mOutPut_EpdResults = Form("file_%s_EpdResults_%d_%s.root",recoEP::mBeamEnergy[energy].c_str(),mEpdMode,jobId.c_str());  
  mInPut_EpdCorrections = Form("%sfile_%s_EpdCorrections_%d.root",inputDir,recoEP::mBeamEnergy[energy].c_str(),mInputEpdMode);
  mOutPut_EpdCorrections = Form("file_%s_EpdCorrections_%d_%s.root",recoEP::mBeamEnergy[energy].c_str(),mEpdMode,jobId.c_str());

  if(mMode == 0) // fill Re-Center Correction Parameters for TPC East/West/Full
  {
    mOutPut_ReCenterPar = Form("./file_%s_ReCenterParameter_%s.root",recoEP::mBeamEnergy[energy].c_str(),jobId.c_str());
  }
  if(mMode == 1) // fill Shift Correction Parameters for TPC East/West/Full
  {
    mOutPut_ShiftPar = Form("./file_%s_ShiftParameter_%s.root",recoEP::mBeamEnergy[energy].c_str(),jobId.c_str());
  }
  if(mMode == 2) // fill Event Plane Resolution for TPC East/West
  {
    mOutPut_Resolution = Form("./file_%s_Resolution_%s.root",recoEP::mBeamEnergy[energy].c_str(),jobId.c_str());
  }
  if(mMode == 3) // fill Charged Flow
  {
    mOutPut_ChargedFlow = Form("./file_%s_ChargedFlow_%s.root",recoEP::mBeamEnergy[energy].c_str(),jobId.c_str());
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
  mTpcEpManager = new StTpcEpManager(mEnergy); // initialize TPC EP Manager
  
  mEpdHits = new TClonesArray("StPicoEpdHit");

  mPicoDstChain = mPicoDstMaker->chain(); 
  cout << "Chain counts = " << mPicoDstMaker->chain()->GetEntries() << endl;
  
  unsigned int found;
  mPicoDstChain->SetBranchStatus("EpdHit*",1,&found);   // note you need the asterisk
  cout << "EpdHit Branch returned found= " << found << endl;
  mPicoDstChain->SetBranchAddress("EpdHit",&mEpdHits);

  mEpFinder = new StEpdEpFinder(10,mOutPut_EpdCorrections.c_str(),mInPut_EpdCorrections.c_str()); // Set the nEventType to 1, since we do not have separation for centrality yet 
  mEpFinder->SetnMipThreshold(mMipThresh);  // recommended by EPD group
  mEpFinder->SetMaxTileWeight(mMaxTile);    // recommended by EPD group
  mEpFinder->SetEpdHitFormat(2);           // 2=pico
   
  if(mEpdMode > 1)
  {
    mFile_EpdResults = new TFile(mOutPut_EpdResults.c_str(),"RECREATE");
    mEventPlaneHistoManager->initEpdEpResults(); 
    //if(mEpdMode == 2) mEventPlaneProManager->initEpdEtaWeighting();
  }

  //return kStOK; 
  if(!mRefMultCorr)
  {
    //if(!mEventPlaneCut->isBES()) mRefMultCorr = CentralityMaker::instance()->getgRefMultCorr_Run14_AuAu200_VpdMB5_P16id(); // 200GeV_2014
    if(mEventPlaneCut->isBES()) mRefMultCorr = CentralityMaker::instance()->getRefMultCorr(); // BESII
  }

  if(mMode == 0)
  { // fill ReCenter Correction Parameters for ZDC-SMD East/West & TPC East/West/Full
    mFile_ReCenterPar = new TFile(mOutPut_ReCenterPar.c_str(),"RECREATE");

    mEventPlaneProManager->initTpcReCenter(); // TPC
    mEventPlaneHistoManager->initTpcRawEP();
  }
  if(mMode == 1)
  { // fill Shift Correction Parameters for ZDC-SMD East/West & TPC East/West/Full
    mFile_ShiftPar = new TFile(mOutPut_ShiftPar.c_str(),"RECREATE");

    mEventPlaneProManager->initTpcShift(); // TPC
    mEventPlaneHistoManager->initTpcReCenterEP();
    mTpcEpManager->readReCenterCorr(); // read in ReCenter Correction Parameters
    //std::cout << "Passed initialization" << std::endl;
  }
  if(mMode == 2)
  { // fill EP Resolution for ZDC-SMD Sub & TPC Sub/Ran
    mFile_Resolution = new TFile(mOutPut_Resolution.c_str(),"RECREATE");

    mEventPlaneProManager->initTpcResolution(); // TPC
    mEventPlaneHistoManager->initTpcShiftEP();
    mTpcEpManager->readReCenterCorr(); // read in ReCenter Correction Parameters
    mTpcEpManager->readShiftCorr(); // read in ReCenter Correction Parameters
    mUsedTrackCounter = 0; // for Random EP
  }
  if(mMode == 3)
  { // fill Charged Flow
    mFile_ChargedFlow = new TFile(mOutPut_ChargedFlow.c_str(),"RECREATE");

    mEventPlaneProManager->initChargedFlow(); // Charged Flow

    // TPC EP
    mEventPlaneHistoManager->initTpcShiftEP();
    mTpcEpManager->readReCenterCorr(); // read in ReCenter Correction Parameters
    mTpcEpManager->readShiftCorr(); // read in Shift Correction Parameters
    mTpcEpManager->readResolution(); // read in EP Resolution
  }
  
  return kStOK;
}

//----------------------------------------------------------------------------- 
int StEventPlaneMaker::Finish() 
{
  mEpFinder->Finish();
 
  if(mEpdMode > 1)
  { 
    mFile_EpdResults->cd();
    mEventPlaneHistoManager->writeEpdEpResults(); 
    //if (mEpdMode == 2) mEventPlaneHistoManager->writeEpdEtaWeighting(); 
    mFile_EpdResults->Close();
  }
  
  if(mMode == 0)
  {
    if(mOutPut_ReCenterPar != "")
    {
      mFile_ReCenterPar->cd();

      mEventPlaneProManager->writeTpcReCenter(); // TPC
      mEventPlaneHistoManager->writeTpcRawEP();

      mFile_ReCenterPar->Close();
    }
  }
  if(mMode == 1)
  {
    if(mOutPut_ShiftPar != "")
    {
      mFile_ShiftPar->cd();

      mEventPlaneProManager->writeTpcShift(); // TPC
      mEventPlaneHistoManager->writeTpcReCenterEP();

      mFile_ShiftPar->Close();
    }
  }
  if(mMode == 2)
  {
    if(mOutPut_Resolution != "")
    {
      mFile_Resolution->cd();

      mEventPlaneProManager->writeTpcResolution(); // TPC
      mEventPlaneHistoManager->writeTpcShiftEP();

      mFile_Resolution->Close();
    }
  }
  if(mMode == 3)
  {
    if(mOutPut_ChargedFlow != "")
    {
      mFile_ChargedFlow->cd();

      mEventPlaneProManager->writeChargedFlow(); // TPC

      mFile_ChargedFlow->Close();
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
    const TVector3 pv = mPicoEvent->primaryVertex();
    const float vx    = pv.x(); // x works for both TVector3 and StThreeVectorF
    const float vy    = pv.y();
    const float vz    = pv.z();
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
    //if(!mEventPlaneCut->isBES()) mRefMultCorr->initEvent(grefMult,vz,zdcX); // 200GeV_2014
    if(mEventPlaneCut->isBES()) mRefMultCorr->initEvent(refMult,vz,zdcX); // BES-II might need Luminosity corrections
    
    //if(mRefMultCorr->isBadRun(runId))
    //{
    //   LOG_ERROR << "Bad Run from StRefMultCorr! Skip!" << endm;
    //   return kStErr;
    //}
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

    bool isPileUpEventStEventPlaneCut; //= mEventPlaneCut->isPileUpEvent(grefMult,numOfBTofMatch,numOfBTofHits); // 200GeV
    if(mEventPlaneCut->isBES()) isPileUpEventStEventPlaneCut = mEventPlaneCut->isPileUpEvent(refMult,numOfBTofMatch,numOfBTofHits); // 54 GeV | always return false for 27 GeV
    //bool isPileUpEventStRefMultCorr = !mRefMultCorr->passnTofMatchRefmultCut(1.0*refMult, 1.0*numOfBTofMatch); // 27 GeV | always return !true for other energies
    bool isPileUpEvent = isPileUpEventStEventPlaneCut;// || isPileUpEventStRefMultCorr;
    // cout << "isPileUpEvent = " << isPileUpEvent << ", isPileUpEventStEventPlaneCut = " << isPileUpEventStEventPlaneCut << ", isPileUpEventStRefMultCorr = " << isPileUpEventStRefMultCorr << endl;

    //if(mEventPlaneCut->passEventCut(mPicoDst) && !isPileUpEvent && cent9 > -0.5)
    if(mEventPlaneCut->passEventCut(mPicoDst) && !isPileUpEvent && cent9 > -0.5)
    { // apply Event Cuts for anlaysis 
      mTpcEpManager->initTpcEp(cent9,runIndex,vzSign); // TPC EP
      
      
      StEpdEpInfo result = mEpFinder->Results(mEpdHits,pv,cent9); // Centrality bin 0 
      if(mEpdMode > 1) mEventPlaneHistoManager->fillEpdEpResults(result,cent9);

        
      // calculate TPC QVector
      if(mMode == 0)
      { // raw QVector
	for(unsigned int i_track = 0; i_track < nTracks; ++i_track)
	{ // calculate number of tracks used in raw EP reconstruction
	  StPicoTrack *picoTrack = (StPicoTrack*)mPicoDst->track(i_track); // get picoTrack
	  if(mEventPlaneCut->passTrackEp(picoTrack,mPicoEvent))
	  { // track cut for EP reconstruction
	    if(mTpcEpManager->passTrackEpEast(picoTrack)) // East sub EP 
	    {
	      mTpcEpManager->addTrackEastRaw(picoTrack);
	    }
	    if(mTpcEpManager->passTrackEpWest(picoTrack)) // West sub EP 
	    {
	      mTpcEpManager->addTrackWestRaw(picoTrack);
	    }
	    if(mTpcEpManager->passTrackEpFull(picoTrack)) // Full EP 
	    {
	      mTpcEpManager->addTrackFullRaw(picoTrack);
	    }
	  }
	}
      }
      if(mMode > 0)
      { // re-centered QVector
	for(unsigned int i_track = 0; i_track < nTracks; ++i_track)
	{ // calculate QVector after recenter correction
	  StPicoTrack *picoTrack = (StPicoTrack*)mPicoDst->track(i_track); // get picoTrack
	  if(mEventPlaneCut->passTrackEp(picoTrack,mPicoEvent))
	  { // track cut for EP reconstruction
            //std::cout << "passTrackEp" << std::endl;
	    if(mTpcEpManager->passTrackEpEast(picoTrack)) // East sub EP 
	    {
              //std::cout << "passTrackEast" << std::endl;
	      mTpcEpManager->addTrackEast(picoTrack);
              //std::cout << "addTrackEast" << std::endl;
	    }
	    if(mTpcEpManager->passTrackEpWest(picoTrack)) // West sub EP 
	    {
              //std::cout << "passTrackWest" << std::endl;
	      mTpcEpManager->addTrackWest(picoTrack);
              //std::cout << "addTrackWest" << std::endl;
	    }
	    if(mTpcEpManager->passTrackEpFull(picoTrack)) // Full EP 
	    {
	      mTpcEpManager->addTrackFull(picoTrack);
              //std::cout << "addTrackFull" << std::endl;
	      if(mMode == 2) mUsedTrackCounter++; // used in random EP
	    }
	  }
	}
        //std::cout << "Added all the tracks" << std::endl;
      }

      if(mMode == 0)
      { // fill ReCenter Correction for ZDC-SMD East/West & TPC East/West/Full

	// TPC: 
	// fill recenter correction parameter & fill raw TPC EP
	TVector2 vTpcQ2East = mTpcEpManager->getQVectorRaw(0); // raw QVec East
	TVector2 vTpcQ2West = mTpcEpManager->getQVectorRaw(1); // raw QVec West
	TVector2 vTpcQ2Full = mTpcEpManager->getQVectorRaw(2); // raw QVec Full

	if(mTpcEpManager->passTrackEtaNumRawCut())
	{ // fill recenter correction for east/west sub EP
	  for(unsigned int i_track = 0; i_track < nTracks; ++i_track)
	  {
	    StPicoTrack *picoTrack = (StPicoTrack*)mPicoDst->track(i_track); // get picoTrack
	    if(mEventPlaneCut->passTrackEp(picoTrack,mPicoEvent))
	    { // track cut for EP reconstruction
	      TVector3 primMom; // temp fix for StThreeVectorF & TVector3
	      const double primPx    = picoTrack->pMom().x(); // x works for both TVector3 and StThreeVectorF
	      const double primPy    = picoTrack->pMom().y();
	      const double primPz    = picoTrack->pMom().z();
	      primMom.SetXYZ(primPx,primPy,primPz);
	      const double primPt = primMom.Perp(); // track pT
	      if(mTpcEpManager->passTrackEpEast(picoTrack)) // East sub EP 
	      {
		TVector2 vTpcq2East = mTpcEpManager->calq2Vector(picoTrack);
		mEventPlaneProManager->fillTpcReCenterEast(vTpcq2East,cent9,runIndex,vzSign,primPt); // fill recenter
	      }
	      if(mTpcEpManager->passTrackEpWest(picoTrack)) // West sub EP 
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
	    if(mEventPlaneCut->passTrackEp(picoTrack,mPicoEvent))
	    { // track cut for EP reconstruction
	      TVector3 primMom; // temp fix for StThreeVectorF & TVector3
	      const double primPx    = picoTrack->pMom().x(); // x works for both TVector3 and StThreeVectorF
	      const double primPy    = picoTrack->pMom().y();
	      const double primPz    = picoTrack->pMom().z();
	      primMom.SetXYZ(primPx,primPy,primPz);
	      const double primPt = primMom.Perp(); // track pT
	      if(mTpcEpManager->passTrackEpFull(picoTrack)) // Full EP 
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

      if(mMode == 1)
      { // fill Shift Correction for ZDC-SMD East/West & TPC East/West/Full

	// TPC: 
	// apply recenter correction to East/West/Full
	// fill shift correction parameter East/West/Full & fill recenter TPC EP
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
      if(mMode == 2)
      { // fill EP Resolution TPC Sub/Ran
	// TPC: 
	// apply recenter correction & apply shift correction to East/West/Full
	// fill EP Resolution for Sub/Ran & fill Shift TPC EP
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
	  if(mEventPlaneCut->passTrackEp(picoTrack,mPicoEvent))
	  { // track cut for EP reconstruction
	    TVector3 primMom; // temp fix for StThreeVectorF & TVector3
	    const double primPx    = picoTrack->pMom().x(); // x works for both TVector3 and StThreeVectorF
	    const double primPy    = picoTrack->pMom().y();
	    const double primPz    = picoTrack->pMom().z();
	    primMom.SetXYZ(primPx,primPy,primPz);
	    const double primPt = primMom.Perp(); // track pT
	    if(mTpcEpManager->passTrackEpFull(picoTrack)) // Full EP 
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
	    double tpcPsi2East = mTpcEpManager->calShiftAngle2East();
	    double tpcPsi2West = mTpcEpManager->calShiftAngle2West();
	    mEventPlaneProManager->fillTpcResSub(tpcPsi2East,tpcPsi2West,cent9,runIndex);
	    mEventPlaneHistoManager->fillTpcShiftSubEP(tpcPsi2East,tpcPsi2West,cent9,runIndex);
	  }
	}
	if(mTpcEpManager->passTrackFullNumCut())
	{ // fill shift correction for full EP
	  if( !(vTpcQ2Full.Mod() < 1e-10) )
	  {
	    double tpcPsi2RanA = mTpcEpManager->calShiftAngle2RanA();
	    double tpcPsi2RanB = mTpcEpManager->calShiftAngle2RanB();
	    double tpcPsi2Full = mTpcEpManager->calShiftAngle2Full();
	    mEventPlaneProManager->fillTpcResRan(tpcPsi2RanA,tpcPsi2RanB,cent9,runIndex);
	    mEventPlaneHistoManager->fillTpcShiftRanEP(tpcPsi2RanA,tpcPsi2RanB,cent9,runIndex);
	    mEventPlaneHistoManager->fillTpcShiftFullEP(tpcPsi2Full,cent9,runIndex);
	  }
	}
      }
      if(mMode == 3)
      { // fill Charged v1Pp & v2Ep & v2Pp & v3Ep
	// TPC: 
	// apply recenter correction & apply shift correction to East/West/Full
	// fill v2Ep & v3Ep (later)
	TVector2 vTpcQ2East = mTpcEpManager->getQVector(0); // receter QVec East
	TVector2 vTpcQ2West = mTpcEpManager->getQVector(1); // receter QVec West
	TVector2 vTpcQ2Full = mTpcEpManager->getQVector(2); // receter QVec Full

	if(mTpcEpManager->passTrackEtaNumCut())
	{
	  if( !(vTpcQ2East.Mod() < 1e-10 || vTpcQ2West.Mod() < 1e-10) )
	  {
	    const double tpcPsi2East = mTpcEpManager->calShiftAngle2East();
	    const double tpcPsi2West = mTpcEpManager->calShiftAngle2West();
	    const double tpcRes2Sub = mTpcEpManager->getRes2Sub(cent9);
	    for(unsigned int i_track = 0; i_track < nTracks; ++i_track)
	    { // calculate QVector after recenter correction
	      StPicoTrack *picoTrack = (StPicoTrack*)mPicoDst->track(i_track); // get picoTrack
	      TVector3 primMom; // temp fix for StThreeVectorF & TVector3
	      const double primPx    = picoTrack->pMom().x(); // x works for both TVector3 and StThreeVectorF
	      const double primPy    = picoTrack->pMom().y();
	      const double primPz    = picoTrack->pMom().z();
	      primMom.SetXYZ(primPx,primPy,primPz);
	      const double primPt = primMom.Perp(); // track pT
	      const double primEta = primMom.PseudoRapidity(); // track eta
	      const double primPhi = primMom.Phi(); // track eta
	      if(mEventPlaneCut->passTrackChargedFlowEast(picoTrack,mPicoEvent))
	      {
		const double v2Ep = TMath::Cos(2.0*(primPhi-tpcPsi2West));
		mEventPlaneProManager->fillChargedV2Ep(primPt, v2Ep, tpcRes2Sub, cent9, runIndex, reweight);
	      }
	      if(mEventPlaneCut->passTrackChargedFlowWest(picoTrack,mPicoEvent))
	      {
		const double v2Ep = TMath::Cos(2.0*(primPhi-tpcPsi2East));
		mEventPlaneProManager->fillChargedV2Ep(primPt, v2Ep, tpcRes2Sub, cent9, runIndex, reweight);
	      }
	    }
	  }
	}
      }

      mTpcEpManager->clearTpcEp();
    }
  }

  return kStOK;
}

