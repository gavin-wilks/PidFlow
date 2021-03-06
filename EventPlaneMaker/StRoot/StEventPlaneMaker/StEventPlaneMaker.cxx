#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoEvent/StPicoEpdHit.h"
#include "StRefMultCorr/StRefMultCorr.h"
#include "StRefMultCorr/CentralityMaker.h"
#include "StEpdUtil/StEpdEpFinder.h"
#include "StEpdUtil/StEpdEpInfo.h"
#include "StEpdUtil/StEpdGeom.h"
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
#include "StRoot/StEventPlaneMaker/StEpdEpManager.h"

#include <algorithm>
//#include <filesystem>
#include <fstream>
#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <TVector3.h>
#include <TMath.h>

ClassImp(StEventPlaneMaker)

//-----------------------------------------------------------------------------
StEventPlaneMaker::StEventPlaneMaker(const char* name, StPicoDstMaker *picoMaker, const string jobId, const int Mode, const int EpdMode, const int energy, float mipThresh = 0.3, float maxTile = 2.0) : StMaker(name)
{
  mPicoDstMaker = picoMaker;
  mPicoDst = NULL;
  mRefMultCorr = NULL;
  mEpFinder = NULL;
  mMode = Mode;
  mEpdMode = EpdMode;  
  mEnergy = energy;
  mMipThresh = mipThresh;
  mMaxTile = maxTile;
    
  if(mEpdMode == 0)
  {
    mInPut_EpdCorrections = "";
    mOutPut_EpdCorrections = Form("file_%s_EpdCorrections_%d_%s.root",recoEP::mBeamEnergy[energy].c_str(),mEpdMode,jobId.c_str());
  }
  if(mEpdMode == 1)
  {
    mInPut_EpdCorrections = Form("StRoot/StEventPlaneUtility/EpdCorrections/file_%s_EpdCorrections_0.root",recoEP::mBeamEnergy[energy].c_str());
    mOutPut_EpdCorrections = Form("file_%s_EpdCorrections_%d_%s.root",recoEP::mBeamEnergy[energy].c_str(),mEpdMode,jobId.c_str());
  }  
  if(mEpdMode == 2)
  {
    mInPut_EpdCorrections = Form("StRoot/StEventPlaneUtility/EpdCorrections/file_%s_EpdCorrections_1.root",recoEP::mBeamEnergy[energy].c_str());
    mOutPut_EpdCorrections = Form("file_%s_EpdCorrections_%d_%s.root",recoEP::mBeamEnergy[energy].c_str(),mEpdMode,jobId.c_str());
    mOutPut_EpdResults = Form("file_%s_EpdResEta_%d_%s.root",recoEP::mBeamEnergy[energy].c_str(),mEpdMode,jobId.c_str());  
    //cout << "successfully instantised the epd" << endl;
  }
  if(mEpdMode == 3)
  {
    mInPut_EpdCorrections = Form("StRoot/StEventPlaneUtility/EpdCorrections/file_%s_EpdCorrections_2.root",recoEP::mBeamEnergy[energy].c_str());
    mOutPut_EpdCorrections = Form("file_%s_EpdCorrections_%d_%s.root",recoEP::mBeamEnergy[energy].c_str(),mEpdMode,jobId.c_str());
    mOutPut_EpdResults = Form("file_%s_EpdFlowEta_%d_%s.root",recoEP::mBeamEnergy[energy].c_str(),mEpdMode,jobId.c_str());  
  }
  if(mEpdMode == 4)
  {
    mInPut_EpdCorrections = Form("StRoot/StEventPlaneUtility/EpdCorrections/file_%s_EpdCorrections_3.root",recoEP::mBeamEnergy[energy].c_str());
    mOutPut_EpdCorrections = Form("file_%s_EpdCorrections_%d_%s.root",recoEP::mBeamEnergy[energy].c_str(),mEpdMode,jobId.c_str());
    mOutPut_EpdResults = Form("file_%s_EpdFlow_%d_%s.root",recoEP::mBeamEnergy[energy].c_str(),mEpdMode,jobId.c_str());  
  }

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
  mEpdEpManager = new StEpdEpManager(mEnergy); // initialize EPD EP Manager  

  //mEpdHits = new TClonesArray("StPicoEpdHit");

  //mPicoDstChain = mPicoDstMaker->chain(); 
  //cout << "Chain counts = " << mPicoDstMaker->chain()->GetEntries() << endl;
  
  //unsigned int found;
  //mPicoDstChain->SetBranchStatus("EpdHit*",1,&found);   // note you need the asterisk
  //cout << "EpdHit Branch returned found= " << found << endl;
  //mPicoDstChain->SetBranchAddress("EpdHit",&mEpdHits);

  mEpFinder = new StEpdEpFinder(10,mOutPut_EpdCorrections.c_str(),mInPut_EpdCorrections.c_str()); // Set the nEventType to 1, since we do not have separation for centrality yet 
  mEpFinder->SetnMipThreshold(mMipThresh);  // recommended by EPD group
  mEpFinder->SetMaxTileWeight(mMaxTile);    // recommended by EPD group
  mEpFinder->SetEpdHitFormat(2);           // 2=pico
   
  mEpdGeom = new StEpdGeom;

  double mEtaLin[3][9];
  double mEtaCub[3][9];
                     
  if(mEpdMode == 3 || mEpdMode == 5)
  {
    std::ifstream etaweights(Form("StRoot/StEventPlaneUtility/EtaWeights/file_%s_EtaWeights.txt",recoEP::mBeamEnergy[mEnergy].c_str()));
 
    for(int order = 1; order <= recoEP::mEpdEpOrder; order++)
    {
      for(int i_cent = 0; i_cent < 9; i_cent++)
      {
        etaweights >> mEtaLin[order-1][i_cent]; etaweights >> mEtaCub[order-1][i_cent];
      }
    } 
  } 
  if(mEpdMode == 2)
  {
    mFile_EpdResults = new TFile(mOutPut_EpdResults.c_str(),"RECREATE");
 
    mEventPlaneProManager->initEpdRes();
    mEventPlaneProManager->initEpdFlowEta();
    mEventPlaneHistoManager->initEpdEp();
    mEventPlaneHistoManager->initEpdQ();
 
    //cout << "epd is initialized" << endl;
  }
  if(mEpdMode == 3)
  {
    TH2D wt1("Order1etaWeight","Order1etaWeight",100,-5.5,5.5,10,0,10);
    for (int ix=1; ix<101; ix++){
      for (int iy=1; iy<11; iy++){
        double eta = wt1.GetXaxis()->GetBinCenter(ix);
        wt1.SetBinContent(ix,iy,mEtaLin[0][iy-1]*eta+mEtaCub[0][iy-1]*pow(eta,3));
      }
    }
    mEpFinder->SetEtaWeights(1,wt1);    
    
    /*TH2D wt2("Order2etaWeight","Order2etaWeight",100,-5.5,5.5,10,0,10);
    for (int ix=1; ix<101; ix++){
      for (int iy=1; iy<11; iy++){
        double eta = wt2.GetXaxis()->GetBinCenter(ix);
        wt2.SetBinContent(ix,iy,mEtaLin[1][iy-1]*eta+mEtaCub[1][iy-1]*pow(eta,3));
      }
    }
    mEpFinder->SetEtaWeights(2,wt2); 
    
    TH2D wt3("Order3etaWeight","Order3etaWeight",100,-5.5,5.5,10,0,10);
    for (int ix=1; ix<101; ix++){
      for (int iy=1; iy<11; iy++){
        double eta = wt3.GetXaxis()->GetBinCenter(ix);
        wt3.SetBinContent(ix,iy,mEtaLin[2][iy-1]*eta+mEtaCub[2][iy-1]*pow(eta,3));
      }
    }
    mEpFinder->SetEtaWeights(3,wt3);   */
  }
  if(mEpdMode == 4)
  {
    mFile_EpdResults = new TFile(mOutPut_EpdResults.c_str(),"RECREATE");

    mEventPlaneProManager->initEpdRes();
    mEventPlaneProManager->initEpdFlow();
    mEventPlaneProManager->initEpdFlowEta();
    mEventPlaneHistoManager->initEpdEp();
    mEventPlaneHistoManager->initEpdQ();
  }

  if(!mRefMultCorr)
  {
    if(mEventPlaneCut->isBES()) mRefMultCorr = CentralityMaker::instance()->getRefMultCorr(); // BESII
  }

  if(mMode == 0)
  { // fill ReCenter Correction Parameters for TPC East/West/Full
    mFile_ReCenterPar = new TFile(mOutPut_ReCenterPar.c_str(),"RECREATE");

    mEventPlaneProManager->initTpcReCenter(); // TPC
    mEventPlaneHistoManager->initTpcRawEP();
  }
  if(mMode == 1)
  { // fill Shift Correction Parameters for TPC East/West/Full
    mFile_ShiftPar = new TFile(mOutPut_ShiftPar.c_str(),"RECREATE");

    mEventPlaneProManager->initTpcShift(); // TPC
    mEventPlaneHistoManager->initTpcReCenterEP();
    mTpcEpManager->readReCenterCorr(); // read in ReCenter Correction Parameters
  }
  if(mMode == 2)
  { // fill EP Resolution for TPC Sub/Ran
    mFile_Resolution = new TFile(mOutPut_Resolution.c_str(),"RECREATE");

    mEventPlaneProManager->initTpcResolution(); // TPC
    mEventPlaneHistoManager->initTpcShiftEP();
    mEventPlaneProManager->initTpcFlowEta();
    mTpcEpManager->readReCenterCorr(); // read in ReCenter Correction Parameters
    mTpcEpManager->readShiftCorr(); // read in ReCenter Correction Parameters
    mUsedTrackCounter = 0; // for Random EP
    //cout << "initialized the TPC" << endl;
  }
  if(mMode == 3)
  { // fill Charged Flow
    mFile_ChargedFlow = new TFile(mOutPut_ChargedFlow.c_str(),"RECREATE");

    mEventPlaneProManager->initChargedFlow(); // Charged Flow

    // TPC EP
    mEventPlaneHistoManager->initTpcShiftEP();
    //std::cout << "Initialization done" << std::endl;
    mTpcEpManager->readReCenterCorr(); // read in ReCenter Correction Parameters
    //std::cout << "Read in recenter correction" << std::endl;
    mTpcEpManager->readShiftCorr(); // read in Shift Correction Parameters
    //std::cout << "Read in shift correction" << std::endl;
    mTpcEpManager->readResolution(); // read in EP Resolution
    //std::cout << "Read in resolution" << std::endl;
  }
  
  return kStOK;
}

//----------------------------------------------------------------------------- 
int StEventPlaneMaker::Finish() 
{
  mEpFinder->Finish();
 
  if(mEpdMode == 2)
  { 
    mFile_EpdResults->cd();
    mEventPlaneProManager->writeEpdRes();
    mEventPlaneProManager->writeEpdFlowEta();
    mFile_EpdResults->Close();
  }
  if(mEpdMode == 4)
  { 
    mFile_EpdResults->cd();
    mEventPlaneProManager->writeEpdRes();
    mEventPlaneProManager->writeEpdFlow();
    mEventPlaneProManager->writeEpdFlowEta();
    mEventPlaneHistoManager->writeEpdEp(); 
    mEventPlaneHistoManager->writeEpdQ();
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
      mEventPlaneProManager->writeTpcFlowEta();

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
    //cout << "We are in the make function" << endl;
    //if(mEventPlaneCut->passEventCut(mPicoDst) && !isPileUpEvent && cent9 > -0.5)
    if(mEventPlaneCut->passEventCut(mPicoDst) && !isPileUpEvent && cent9 > -0.5)
    { // apply Event Cuts for anlaysis 
      //cout << "Event passed the cut" << endl;
      StEpdEpInfo result;
      if (mEpdMode < 5 && mEpdMode >= 0) result = mEpFinder->Results(mPicoDst->picoArray(StPicoArrays::EpdHit),pv,cent9);
      
      mEpdEpManager->initEpdEp(cent9,runIndex);
      
      //mEpdMode == 0 sets phi-weighting and mEpdMode == 1 sets the psi shift corrections
      if(mEpdMode == 2)
      {// Fill TProfiles for EP resolution and numerator for flow calculations.
       // These two TProfile plots will be used to set eta weights in an external macro.
        //cout << "We are in mEpdMode == 2" << endl;

        mEventPlaneProManager->fillEpdRes(result,cent9,runIndex); 
        
        if(mPicoDst->numberOfEpdHits() > 0)
        {
          cout << "Number of Epd Hits = " << mPicoDst->numberOfEpdHits() << endl;
          for(unsigned int iEpd = 0; iEpd < mPicoDst->numberOfEpdHits(); iEpd++) 
          {
 
            //cout << "Before grabbing the epdHit" << endl;

            StPicoEpdHit *EpdHit = mPicoDst->epdHit(iEpd);
            if (! EpdHit) continue;

            //cout << "We have an EPD hit" << endl;          

            TVector3 pos = mEpdGeom->RandomPointOnTile(EpdHit->id());
            TVector3 StraightLine = pos - pv;
            double phi = StraightLine.Phi();
            while(phi < 0.0) phi += 2.0*TMath::Pi(); 
            double eta = StraightLine.Eta();
            if (!(fabs(eta) > 0) || (fabs(eta) > 1000)) continue;
      
            //cout << "We have found the position of the hit" << endl;
  
            int ew = (EpdHit->id() < 0)? 0 : 1;  //is EPD east or west

            double nMip = EpdHit->nMIP(); // Weight of the track
            double nMipEff = nMip; // Put a max on nMip in calculation to prevent slow tracks from overcontributing 
            if(nMipEff > 3.0) nMipEff = 3.0;
            else if(nMipEff < 0.3) continue;

            //cout << "Determined the weight of the track" << endl;

            double v1, v2, v3; 
            if (1 == ew){  //west
              v1 = TMath::Cos(1.0*(phi-result.EastPhiWeightedAndShiftedPsi(1)));
              v2 = TMath::Cos(2.0*(phi-result.EastPhiWeightedAndShiftedPsi(2)));
              v3 = TMath::Cos(3.0*(phi-result.EastPhiWeightedAndShiftedPsi(3)));
 
              mEventPlaneProManager->fillEpdFlowEta(eta,v1,cent9,1,nMipEff);
              mEventPlaneProManager->fillEpdFlowEta(eta,v2,cent9,2,nMipEff); 
              mEventPlaneProManager->fillEpdFlowEta(eta,v3,cent9,3,nMipEff);

              //cout << "Filling west flow" << endl;
            }
            else if(0 == ew){ //east
              v1 = TMath::Cos(1.0*(phi-result.WestPhiWeightedAndShiftedPsi(1)));
              v2 = TMath::Cos(2.0*(phi-result.WestPhiWeightedAndShiftedPsi(2)));
              v3 = TMath::Cos(3.0*(phi-result.WestPhiWeightedAndShiftedPsi(3)));
 
              mEventPlaneProManager->fillEpdFlowEta(eta,v1,cent9,1,nMipEff);
              mEventPlaneProManager->fillEpdFlowEta(eta,v2,cent9,2,nMipEff); 
              mEventPlaneProManager->fillEpdFlowEta(eta,v3,cent9,3,nMipEff);
              //cout << "Filling east flow" << endl;
            } 
          }
          //cout << "Outside of the loop" << endl;
          mEventPlaneHistoManager->fillEpdEp(result, cent9);
          //cout << "Filled EPD EP" << endl; 
          mEventPlaneHistoManager->fillEpdQ(result, cent9);
          //cout << "Filled EPD Q" << endl;
          //cout << "Successfully ran make for EPD" << endl;
        }
      }
      // mEpdMode == 3 recalculates the shift parameters using the proper eta weights (calculated from the flow output of mEpdMode == 2)
      // One must use the original mEpdMode == 0 corrections file for this step
      if(mEpdMode == 4)
      {
        mEventPlaneProManager->fillEpdRes(result,cent9,runIndex); 

        for(unsigned int iEpd = 0; iEpd < mPicoDst->numberOfEpdHits(); iEpd++) 
        {
          StPicoEpdHit *EpdHit = mPicoDst->epdHit(iEpd);
          if (! EpdHit) continue;

          TVector3 pos = mEpdGeom->RandomPointOnTile(EpdHit->id());
          TVector3 StraightLine = pos - pv;
          double phi = StraightLine.Phi();
          while(phi < 0.0) phi += 2.0*TMath::Pi(); 
          double eta = StraightLine.Eta();
          if (!(fabs(eta) > 0) || (fabs(eta) > 1000)) continue;
          int ew = (EpdHit->id() < 0)? 0 : 1;  //is EPD east or west

          double nMip = EpdHit->nMIP(); // Weight of the track
          double nMipEff = nMip; // Put a max on nMip in calculation to prevent slow tracks from overcontributing
          if(nMipEff > 3.0) nMipEff = 3.0;
          else if(nMipEff < 0.3) continue;
        
          double v1, v2, v3; 
          if (1 == ew){  //west
            v1 = TMath::Cos(1.0*(phi-result.EastPhiWeightedAndShiftedPsi(1)));
            v2 = TMath::Cos(2.0*(phi-result.EastPhiWeightedAndShiftedPsi(2)));
            v3 = TMath::Cos(3.0*(phi-result.EastPhiWeightedAndShiftedPsi(3)));

            mEventPlaneProManager->fillEpdFlowEta(eta,v1,cent9,1,nMipEff);
            mEventPlaneProManager->fillEpdFlowEta(eta,v2,cent9,2,nMipEff); 
            mEventPlaneProManager->fillEpdFlowEta(eta,v3,cent9,3,nMipEff);

            mEventPlaneProManager->fillEpdFlow(v1,cent9,1,runIndex);
            mEventPlaneProManager->fillEpdFlow(v2,cent9,2,runIndex); 
            mEventPlaneProManager->fillEpdFlow(v3,cent9,3,runIndex);
          }
          else if(0 == ew){ //east
            v1 = TMath::Cos(1.0*(phi-result.WestPhiWeightedAndShiftedPsi(1)));
            v2 = TMath::Cos(2.0*(phi-result.WestPhiWeightedAndShiftedPsi(2)));
            v3 = TMath::Cos(3.0*(phi-result.WestPhiWeightedAndShiftedPsi(3)));
   
            mEventPlaneProManager->fillEpdFlowEta(eta,v1,cent9,1,nMipEff);
            mEventPlaneProManager->fillEpdFlowEta(eta,v2,cent9,2,nMipEff); 
            mEventPlaneProManager->fillEpdFlowEta(eta,v3,cent9,3,nMipEff);

            mEventPlaneProManager->fillEpdFlow(v1,cent9,1,runIndex);
            mEventPlaneProManager->fillEpdFlow(v2,cent9,2,runIndex); 
            mEventPlaneProManager->fillEpdFlow(v3,cent9,3,runIndex);
          } 
        }
        mEventPlaneHistoManager->fillEpdEp(result, cent9); 
        mEventPlaneHistoManager->fillEpdQ(result, cent9);
      } 
   
      mTpcEpManager->initTpcEp(cent9,runIndex,vzSign); // TPC EP
      
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
	    if(mTpcEpManager->passTrackEpEast(picoTrack)) // East sub EP 
	    {
	      mTpcEpManager->addTrackEast(picoTrack);
	    }
	    if(mTpcEpManager->passTrackEpWest(picoTrack)) // West sub EP 
	    {
	      mTpcEpManager->addTrackWest(picoTrack);
	    }
	    if(mTpcEpManager->passTrackEpFull(picoTrack)) // Full EP 
	    {
	      mTpcEpManager->addTrackFull(picoTrack);
	      if(mMode == 2) mUsedTrackCounter++; // used in random EP
	    }
	  }
	}
        //std::cout << "Added all the tracks" << std::endl;
      }

      if(mMode == 0)
      { // fill ReCenter Correction for TPC East/West/Full

	// TPC: 
	// fill recenter correction parameter & fill raw TPC EP
        TVector2 vTpcQ1East = mTpcEpManager->getQVectorRaw(0,1); // raw QVec East
        TVector2 vTpcQ1West = mTpcEpManager->getQVectorRaw(1,1); // raw QVec West
        TVector2 vTpcQ1Full = mTpcEpManager->getQVectorRaw(2,1); // raw QVec Full
  
	TVector2 vTpcQ2East = mTpcEpManager->getQVectorRaw(0,2); // raw QVec East
	TVector2 vTpcQ2West = mTpcEpManager->getQVectorRaw(1,2); // raw QVec West
	TVector2 vTpcQ2Full = mTpcEpManager->getQVectorRaw(2,2); // raw QVec Full
  
        TVector2 vTpcQ3East = mTpcEpManager->getQVectorRaw(0,3); // raw QVec East
	TVector2 vTpcQ3West = mTpcEpManager->getQVectorRaw(1,3); // raw QVec West
	TVector2 vTpcQ3Full = mTpcEpManager->getQVectorRaw(2,3); // raw QVec Full
   
        if(mTpcEpManager->passTrackEtaNumRawCut(1))
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
		TVector2 vTpcq1East = mTpcEpManager->calq1Vector(picoTrack);
		mEventPlaneProManager->fillTpcReCenterEast(1,vTpcq1East,cent9,runIndex,vzSign,primPt); // fill recenter
	      }
	      if(mTpcEpManager->passTrackEpWest(picoTrack)) // West sub EP 
	      {
		TVector2 vTpcq1West = mTpcEpManager->calq1Vector(picoTrack);
		mEventPlaneProManager->fillTpcReCenterWest(1,vTpcq1West,cent9,runIndex,vzSign,primPt); // fill recenter
	      }
	    }
	  }
	  if( !(vTpcQ1East.Mod() < 1e-10 || vTpcQ1West.Mod() < 1e-10) )
	  {
	    mEventPlaneHistoManager->fillTpcRawSubEP(1,vTpcQ1East,vTpcQ1West,cent9,runIndex);
	  }
	}
	if(mTpcEpManager->passTrackFullNumRawCut(1))
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
		TVector2 vTpcq1Full = mTpcEpManager->calq1Vector(picoTrack);
		mEventPlaneProManager->fillTpcReCenterFull(1,vTpcq1Full,cent9,runIndex,vzSign,primPt); // fill recenter
	      }
	    }
	  }
	  if( !(vTpcQ1Full.Mod() < 1e-10) )
	  {
	    mEventPlaneHistoManager->fillTpcRawFullEP(1,vTpcQ1Full,cent9,runIndex);
	  }
	}

	if(mTpcEpManager->passTrackEtaNumRawCut(2))
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
		mEventPlaneProManager->fillTpcReCenterEast(2,vTpcq2East,cent9,runIndex,vzSign,primPt); // fill recenter
	      }
	      if(mTpcEpManager->passTrackEpWest(picoTrack)) // West sub EP 
	      {
		TVector2 vTpcq2West = mTpcEpManager->calq2Vector(picoTrack);
		mEventPlaneProManager->fillTpcReCenterWest(2,vTpcq2West,cent9,runIndex,vzSign,primPt); // fill recenter
	      }
	    }
	  }
	  if( !(vTpcQ2East.Mod() < 1e-10 || vTpcQ2West.Mod() < 1e-10) )
	  {
	    mEventPlaneHistoManager->fillTpcRawSubEP(2,vTpcQ2East,vTpcQ2West,cent9,runIndex);
	  }
	}
	if(mTpcEpManager->passTrackFullNumRawCut(2))
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
		mEventPlaneProManager->fillTpcReCenterFull(2,vTpcq2Full,cent9,runIndex,vzSign,primPt); // fill recenter
	      }
	    }
	  }
	  if( !(vTpcQ2Full.Mod() < 1e-10) )
	  {
	    mEventPlaneHistoManager->fillTpcRawFullEP(2,vTpcQ2Full,cent9,runIndex);
	  }
	}

        if(mTpcEpManager->passTrackEtaNumRawCut(3))
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
		TVector2 vTpcq3East = mTpcEpManager->calq3Vector(picoTrack);
		mEventPlaneProManager->fillTpcReCenterEast(3,vTpcq3East,cent9,runIndex,vzSign,primPt); // fill recenter
	      }
	      if(mTpcEpManager->passTrackEpWest(picoTrack)) // West sub EP 
	      {
		TVector2 vTpcq3West = mTpcEpManager->calq3Vector(picoTrack);
		mEventPlaneProManager->fillTpcReCenterWest(3,vTpcq3West,cent9,runIndex,vzSign,primPt); // fill recenter
	      }
	    }
	  }
	  if( !(vTpcQ3East.Mod() < 1e-10 || vTpcQ3West.Mod() < 1e-10) )
	  {
	    mEventPlaneHistoManager->fillTpcRawSubEP(3,vTpcQ3East,vTpcQ3West,cent9,runIndex);
	  }
	}
	if(mTpcEpManager->passTrackFullNumRawCut(3))
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
		TVector2 vTpcq3Full = mTpcEpManager->calq3Vector(picoTrack);
		mEventPlaneProManager->fillTpcReCenterFull(3,vTpcq3Full,cent9,runIndex,vzSign,primPt); // fill recenter
	      }
	    }
	  }
	  if( !(vTpcQ3Full.Mod() < 1e-10) )
	  {
	    mEventPlaneHistoManager->fillTpcRawFullEP(3,vTpcQ3Full,cent9,runIndex);
	  }
	}
      }

      if(mMode == 1)
      { // fill Shift Correction for TPC East/West/Full

	// TPC: 
	// apply recenter correction to East/West/Full
	// fill shift correction parameter East/West/Full & fill recenter TPC EP
	TVector2 vTpcQ1East = mTpcEpManager->getQVector(0,1); // receter QVec East
	TVector2 vTpcQ1West = mTpcEpManager->getQVector(1,1); // receter QVec West
	TVector2 vTpcQ1Full = mTpcEpManager->getQVector(2,1); // receter QVec Full

        TVector2 vTpcQ2East = mTpcEpManager->getQVector(0,2); // receter QVec East
	TVector2 vTpcQ2West = mTpcEpManager->getQVector(1,2); // receter QVec West
	TVector2 vTpcQ2Full = mTpcEpManager->getQVector(2,2); // receter QVec Full

        TVector2 vTpcQ3East = mTpcEpManager->getQVector(0,3); // receter QVec East
	TVector2 vTpcQ3West = mTpcEpManager->getQVector(1,3); // receter QVec West
	TVector2 vTpcQ3Full = mTpcEpManager->getQVector(2,3); // receter QVec Full

        if(mTpcEpManager->passTrackEtaNumCut(1))
	{ // fill shift correction for east/west sub EP
	  if( !(vTpcQ1East.Mod() < 1e-10 || vTpcQ1West.Mod() < 1e-10) )
	  {
	    mEventPlaneProManager->fillTpcShiftEast(1,vTpcQ1East,cent9,runIndex,vzSign);
	    mEventPlaneProManager->fillTpcShiftWest(1,vTpcQ1West,cent9,runIndex,vzSign);
	    mEventPlaneHistoManager->fillTpcReCenterSubEP(1,vTpcQ1East,vTpcQ1West,cent9,runIndex);
	  }
	}
	if(mTpcEpManager->passTrackFullNumCut(1))
	{ // fill shift correction for full EP
	  if( !(vTpcQ1Full.Mod() < 1e-10) )
	  {
	    mEventPlaneProManager->fillTpcShiftFull(1,vTpcQ1Full,cent9,runIndex,vzSign);
	    mEventPlaneHistoManager->fillTpcReCenterFullEP(1,vTpcQ1Full,cent9,runIndex);
	  }
	}

	if(mTpcEpManager->passTrackEtaNumCut(2))
	{ // fill shift correction for east/west sub EP
	  if( !(vTpcQ2East.Mod() < 1e-10 || vTpcQ2West.Mod() < 1e-10) )
	  {
	    mEventPlaneProManager->fillTpcShiftEast(2,vTpcQ2East,cent9,runIndex,vzSign);
	    mEventPlaneProManager->fillTpcShiftWest(2,vTpcQ2West,cent9,runIndex,vzSign);
	    mEventPlaneHistoManager->fillTpcReCenterSubEP(2,vTpcQ2East,vTpcQ2West,cent9,runIndex);
	  }
	}
	if(mTpcEpManager->passTrackFullNumCut(2))
	{ // fill shift correction for full EP
	  if( !(vTpcQ2Full.Mod() < 1e-10) )
	  {
	    mEventPlaneProManager->fillTpcShiftFull(2,vTpcQ2Full,cent9,runIndex,vzSign);
	    mEventPlaneHistoManager->fillTpcReCenterFullEP(2,vTpcQ2Full,cent9,runIndex);
	  }
	}
  
        if(mTpcEpManager->passTrackEtaNumCut(3))
	{ // fill shift correction for east/west sub EP
	  if( !(vTpcQ3East.Mod() < 1e-10 || vTpcQ3West.Mod() < 1e-10) )
	  {
	    mEventPlaneProManager->fillTpcShiftEast(3,vTpcQ3East,cent9,runIndex,vzSign);
	    mEventPlaneProManager->fillTpcShiftWest(3,vTpcQ3West,cent9,runIndex,vzSign);
	    mEventPlaneHistoManager->fillTpcReCenterSubEP(3,vTpcQ3East,vTpcQ3West,cent9,runIndex);
	  }
	}
	if(mTpcEpManager->passTrackFullNumCut(3))
	{ // fill shift correction for full EP
	  if( !(vTpcQ3Full.Mod() < 1e-10) )
	  {
	    mEventPlaneProManager->fillTpcShiftFull(3,vTpcQ3Full,cent9,runIndex,vzSign);
	    mEventPlaneHistoManager->fillTpcReCenterFullEP(3,vTpcQ3Full,cent9,runIndex);
	  }
	}
      }

      if(mMode == 2)
      { // fill EP Resolution TPC Sub/Ran
	// TPC: 
	// apply recenter correction & apply shift correction to East/West/Full
	// fill EP Resolution for Sub/Ran & fill Shift TPC EP
	// Random EP
	//cout << "TPC mode = 2" << endl;
	int ranTrack[mUsedTrackCounter];
	double ranCounter = (double)mUsedTrackCounter/2.0 - 1.0;
	for(int i_track = 0; i_track < mUsedTrackCounter; ++i_track)
	{
	  ranTrack[i_track] = i_track;
	}
        //cout << "Filled ranTrack" << endl;
	std::srand(time(0));
	std::random_shuffle(ranTrack,ranTrack+mUsedTrackCounter);
	mUsedTrackCounter = 0;
	for(unsigned int i_track = 0; i_track < nTracks; ++i_track)
	{
          //cout << "In the track loop" << endl;
	  StPicoTrack *picoTrack = (StPicoTrack*)mPicoDst->track(i_track); // get picoTrack
	  if(mEventPlaneCut->passTrackEp(picoTrack,mPicoEvent))
	  { // track cut for EP reconstruction
            //cout << "passTrackEp" << endl;
	    TVector3 primMom; // temp fix for StThreeVectorF & TVector3
	    const double primPx    = picoTrack->pMom().x(); // x works for both TVector3 and StThreeVectorF
	    const double primPy    = picoTrack->pMom().y();
	    const double primPz    = picoTrack->pMom().z();
	    primMom.SetXYZ(primPx,primPy,primPz);
	    const double primPt = primMom.Perp(); // track pT
	    if(mTpcEpManager->passTrackEpFull(picoTrack)) // Full EP 
	    {
              //cout << "passTrackEpFull" << endl;
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

        TVector2 vTpcQ1East = mTpcEpManager->getQVector(0,1); // receter QVec East
	TVector2 vTpcQ1West = mTpcEpManager->getQVector(1,1); // receter QVec West
	TVector2 vTpcQ1Full = mTpcEpManager->getQVector(2,1); // receter QVec Full

	TVector2 vTpcQ2East = mTpcEpManager->getQVector(0,2); // receter QVec East
	TVector2 vTpcQ2West = mTpcEpManager->getQVector(1,2); // receter QVec West
	TVector2 vTpcQ2Full = mTpcEpManager->getQVector(2,2); // receter QVec Full
        
        TVector2 vTpcQ3East = mTpcEpManager->getQVector(0,3); // receter QVec East
	TVector2 vTpcQ3West = mTpcEpManager->getQVector(1,3); // receter QVec West
	TVector2 vTpcQ3Full = mTpcEpManager->getQVector(2,3); // receter QVec Full

        //cout << "Grabbed Q vectors for each order" << endl;

        if(mTpcEpManager->passTrackEtaNumCut(1))
	{ // fill shift correction for east/west sub EP
	  if( !(vTpcQ1East.Mod() < 1e-10 || vTpcQ1West.Mod() < 1e-10) )
	  {
	    double tpcPsi1East = mTpcEpManager->calShiftAngle1East();
	    double tpcPsi1West = mTpcEpManager->calShiftAngle1West();
            //cout << "ORDER 1:  EastPsi = " << tpcPsi1East << "   WestPsi = " << tpcPsi1West << endl;
	    mEventPlaneProManager->fillTpcResSub(1,tpcPsi1East,tpcPsi1West,cent9,runIndex);
	    mEventPlaneHistoManager->fillTpcShiftSubEP(1,tpcPsi1East,tpcPsi1West,cent9,runIndex);
	    
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
	        const double v1Ep = TMath::Cos(1.0*(primPhi-tpcPsi1West));
	        mEventPlaneProManager->fillTpcFlowEta(primEta, v1Ep, cent9, 1, reweight);
	      }
	      if(mEventPlaneCut->passTrackChargedFlowWest(picoTrack,mPicoEvent))
	      {
	        const double v1Ep = TMath::Cos(1.0*(primPhi-tpcPsi1East));
	        mEventPlaneProManager->fillTpcFlowEta(primEta, v1Ep, cent9, 1, reweight);
	      } 
            }
          }
	}
	if(mTpcEpManager->passTrackFullNumCut(1))
	{ // fill shift correction for full EP
	  if( !(vTpcQ1Full.Mod() < 1e-10) )
	  {
	    double tpcPsi1RanA = mTpcEpManager->calShiftAngle1RanA();
	    double tpcPsi1RanB = mTpcEpManager->calShiftAngle1RanB();
	    double tpcPsi1Full = mTpcEpManager->calShiftAngle1Full();
	    mEventPlaneProManager->fillTpcResRan(1,tpcPsi1RanA,tpcPsi1RanB,cent9,runIndex);
	    mEventPlaneHistoManager->fillTpcShiftRanEP(1,tpcPsi1RanA,tpcPsi1RanB,cent9,runIndex);
	    mEventPlaneHistoManager->fillTpcShiftFullEP(1,tpcPsi1Full,cent9,runIndex);
	  }
	} 

	if(mTpcEpManager->passTrackEtaNumCut(2))
	{ // fill shift correction for east/west sub EP
	  if( !(vTpcQ2East.Mod() < 1e-10 || vTpcQ2West.Mod() < 1e-10) )
	  {
	    double tpcPsi2East = mTpcEpManager->calShiftAngle2East();
	    double tpcPsi2West = mTpcEpManager->calShiftAngle2West();

            //cout << "ORDER 2:  EastPsi = " << tpcPsi2East << "   WestPsi = " << tpcPsi2West << endl;
	    mEventPlaneProManager->fillTpcResSub(2,tpcPsi2East,tpcPsi2West,cent9,runIndex);
	    mEventPlaneHistoManager->fillTpcShiftSubEP(2,tpcPsi2East,tpcPsi2West,cent9,runIndex);

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
                mEventPlaneProManager->fillTpcFlowEta(primEta, v2Ep, cent9, 2, reweight); 
	      }
	      if(mEventPlaneCut->passTrackChargedFlowWest(picoTrack,mPicoEvent))
	      {
	        const double v2Ep = TMath::Cos(2.0*(primPhi-tpcPsi2East));
                mEventPlaneProManager->fillTpcFlowEta(primEta, v2Ep, cent9, 2, reweight);
	      }
            }
	  }
	}
	if(mTpcEpManager->passTrackFullNumCut(2))
	{ // fill shift correction for full EP
	  if( !(vTpcQ2Full.Mod() < 1e-10) )
	  {
	    double tpcPsi2RanA = mTpcEpManager->calShiftAngle2RanA();
	    double tpcPsi2RanB = mTpcEpManager->calShiftAngle2RanB();
	    double tpcPsi2Full = mTpcEpManager->calShiftAngle2Full();
	    mEventPlaneProManager->fillTpcResRan(2,tpcPsi2RanA,tpcPsi2RanB,cent9,runIndex);
	    mEventPlaneHistoManager->fillTpcShiftRanEP(2,tpcPsi2RanA,tpcPsi2RanB,cent9,runIndex);
	    mEventPlaneHistoManager->fillTpcShiftFullEP(2,tpcPsi2Full,cent9,runIndex);
	  }
	}

        if(mTpcEpManager->passTrackEtaNumCut(3))
	{ // fill shift correction for east/west sub EP
	  if( !(vTpcQ3East.Mod() < 1e-10 || vTpcQ3West.Mod() < 1e-10) )
	  {
	    double tpcPsi3East = mTpcEpManager->calShiftAngle3East();
	    double tpcPsi3West = mTpcEpManager->calShiftAngle3West();

            //cout << "ORDER 3:  EastPsi = " << tpcPsi3East << "   WestPsi = " << tpcPsi3West << endl;
	    mEventPlaneProManager->fillTpcResSub(3,tpcPsi3East,tpcPsi3West,cent9,runIndex);
	    mEventPlaneHistoManager->fillTpcShiftSubEP(3,tpcPsi3East,tpcPsi3West,cent9,runIndex);

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
	        const double v3Ep = TMath::Cos(3.0*(primPhi-tpcPsi3West));
                mEventPlaneProManager->fillTpcFlowEta(primEta, v3Ep, cent9, 3, reweight);
	      }
	      if(mEventPlaneCut->passTrackChargedFlowWest(picoTrack,mPicoEvent))
	      {
	        const double v3Ep = TMath::Cos(3.0*(primPhi-tpcPsi3East));
                mEventPlaneProManager->fillTpcFlowEta(primEta, v3Ep, cent9, 3, reweight);
	      }
            }
	  }
	}
	if(mTpcEpManager->passTrackFullNumCut(3))
	{ // fill shift correction for full EP
	  if( !(vTpcQ3Full.Mod() < 1e-10) )
	  {
	    double tpcPsi3RanA = mTpcEpManager->calShiftAngle3RanA();
	    double tpcPsi3RanB = mTpcEpManager->calShiftAngle3RanB();
	    double tpcPsi3Full = mTpcEpManager->calShiftAngle3Full();
	    mEventPlaneProManager->fillTpcResRan(3,tpcPsi3RanA,tpcPsi3RanB,cent9,runIndex);
	    mEventPlaneHistoManager->fillTpcShiftRanEP(3,tpcPsi3RanA,tpcPsi3RanB,cent9,runIndex);
	    mEventPlaneHistoManager->fillTpcShiftFullEP(3,tpcPsi3Full,cent9,runIndex);
	  }
	}
      }
      if(mMode == 3)
      { // fill Charged v1Pp & v2Ep & v2Pp & v3Ep
	// TPC: 
	// apply recenter correction & apply shift correction to East/West/Full
	// fill v2Ep & v3Ep (later)
	TVector2 vTpcQ1East = mTpcEpManager->getQVector(0,1); // receter QVec East
	TVector2 vTpcQ1West = mTpcEpManager->getQVector(1,1); // receter QVec West
	TVector2 vTpcQ1Full = mTpcEpManager->getQVector(2,1); // receter QVec Full
	
	TVector2 vTpcQ2East = mTpcEpManager->getQVector(0,2); // receter QVec East
	TVector2 vTpcQ2West = mTpcEpManager->getQVector(1,2); // receter QVec West
	TVector2 vTpcQ2Full = mTpcEpManager->getQVector(2,2); // receter QVec Full

        TVector2 vTpcQ3East = mTpcEpManager->getQVector(0,3); // receter QVec East
	TVector2 vTpcQ3West = mTpcEpManager->getQVector(1,3); // receter QVec West
	TVector2 vTpcQ3Full = mTpcEpManager->getQVector(2,3); // receter QVec Full
        //std::cout << "Grabbed the Q Values" << std::endl;        
 
        if(mTpcEpManager->passTrackEtaNumCut(1))
	{
	  if( !(vTpcQ1East.Mod() < 1e-10 || vTpcQ1West.Mod() < 1e-10) )
	  {
	    const double tpcPsi1East = mTpcEpManager->calShiftAngle1East();
	    const double tpcPsi1West = mTpcEpManager->calShiftAngle1West();
	    const double tpcRes1Sub = mTpcEpManager->getRes1Sub(cent9);
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
		const double v1Ep = TMath::Cos(1.0*(primPhi-tpcPsi1West));
		mEventPlaneProManager->fillChargedV1Ep(primPt, primEta, v1Ep, tpcRes1Sub, cent9, runIndex, reweight);
	      }
	      if(mEventPlaneCut->passTrackChargedFlowWest(picoTrack,mPicoEvent))
	      {
		const double v1Ep = TMath::Cos(1.0*(primPhi-tpcPsi1East));
		mEventPlaneProManager->fillChargedV1Ep(primPt, primEta, v1Ep, tpcRes1Sub, cent9, runIndex, reweight);
	      }
	    }
	  }
	}

	if(mTpcEpManager->passTrackEtaNumCut(2))
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
		mEventPlaneProManager->fillChargedV2Ep(primPt, primEta, v2Ep, tpcRes2Sub, cent9, runIndex, reweight);
	      }
	      if(mEventPlaneCut->passTrackChargedFlowWest(picoTrack,mPicoEvent))
	      {
		const double v2Ep = TMath::Cos(2.0*(primPhi-tpcPsi2East));
		mEventPlaneProManager->fillChargedV2Ep(primPt, primEta, v2Ep, tpcRes2Sub, cent9, runIndex, reweight);
	      }
	    }
	  }
	}

        if(mTpcEpManager->passTrackEtaNumCut(3))
	{
	  if( !(vTpcQ3East.Mod() < 1e-10 || vTpcQ3West.Mod() < 1e-10) )
	  {
	    const double tpcPsi3East = mTpcEpManager->calShiftAngle3East();
	    const double tpcPsi3West = mTpcEpManager->calShiftAngle3West();
	    const double tpcRes3Sub = mTpcEpManager->getRes3Sub(cent9);
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
		const double v3Ep = TMath::Cos(3.0*(primPhi-tpcPsi3West));
		mEventPlaneProManager->fillChargedV3Ep(primPt, primEta, v3Ep, tpcRes3Sub, cent9, runIndex, reweight);
	      }
	      if(mEventPlaneCut->passTrackChargedFlowWest(picoTrack,mPicoEvent))
	      {
		const double v3Ep = TMath::Cos(3.0*(primPhi-tpcPsi3East));
		mEventPlaneProManager->fillChargedV3Ep(primPt, primEta, v3Ep, tpcRes3Sub, cent9, runIndex, reweight);
	      }
	    }
	  }
	}
        //std::cout << "Finished with 3rd order calculations" << std::endl;
      }

      mTpcEpManager->clearTpcEp();
    }
  }

  return kStOK;
}

