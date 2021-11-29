#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoEvent/StPicoBTofPidTraits.h"
// #include "StRefMultCorr/StRefMultCorr.h"
// #include "StRefMultCorr/CentralityMaker.h"
#include "StMessMgr.h"

#include "StRoot/StEventPlaneMaker/StEventPlaneCut.h"
#include "StRoot/StEventPlaneMaker/StEventPlaneCons.h"

#include "TVector3.h"

ClassImp(StEventPlaneCut)

//---------------------------------------------------------------------------------

StEventPlaneCut::StEventPlaneCut(int energy)
{
  mEnergy = energy;
}

//---------------------------------------------------------------------------------

StEventPlaneCut::~StEventPlaneCut()
{
  /* */
}

//---------------------------------------------------------------------------------

// Event Cuts
bool StEventPlaneCut::isMinBias(StPicoEvent *picoEvent)
{
  // std::cout << "year: " << picoEvent->year() << std::endl;
  // std::cout << "day: " << picoEvent->day() << std::endl;
  // std::cout << "triggerIds: " << picoEvent->triggerIds()[0] << std::endl;  
  if(mEnergy == 0 && recoEP::mBeamYear[mEnergy] == picoEvent->year() && !( picoEvent->isTrigger(650001) || picoEvent->isTrigger(650011) || picoEvent->isTrigger(650021) || picoEvent->isTrigger(650031) || picoEvent->isTrigger(650041) || picoEvent->isTrigger(650051) )) return false; // 14p5GeV_2019
  if(mEnergy == 1 && recoEP::mBeamYear[mEnergy] == picoEvent->year() && !( picoEvent->isTrigger(640001) || picoEvent->isTrigger(640011) || picoEvent->isTrigger(640021) || picoEvent->isTrigger(640031) || picoEvent->isTrigger(640041) || picoEvent->isTrigger(640051) )) return false; // 19p6GeV_2019
  return true;
}

bool StEventPlaneCut::isBES()
{
  if(mEnergy == 0|| mEnergy == 1) return true; // is BES for 14.5 GeV, 19.6 GeV

  return false; // not BES
}

bool StEventPlaneCut::isPileUpEvent(int refMult, int numOfBTofMatch, int numOfBTofHits)
{
  if(mEnergy == 0)
  {        
    if(numOfBTofMatch > recoEP::mMatchedToFMin[mEnergy]) return kFALSE;         
  }

  if(mEnergy == 1)
  {
    if(numOfBTofMatch <= recoEP::mMatchedToFMin[mEnergy]) return kTRUE; 
        
    // 5th order polynomial coefficients for lower cut, l
    double p0l = -13.11    ;   
    double p1l = 0.8207    ;   
    double p2l = -4.241e-3 ;
    double p3l = 2.81e-5   ;   
    double p4l = -6.434e-8 ;
    double p5l = 4.833e-11 ;
    // 5th order polynomial coefficients for higher cut, h
    double p0h = 10.07     ;   
    double p1h = 1.417     ;   
    double p2h = 1.979e-4  ;
    double p3h = -4.87e-6  ;
    double p4h = 1.109e-8  ;
    double p5h = -1.111e-11;

    double refLow  = p0l + p1l*numOfBTofMatch + p2l*pow(numOfBTofMatch,2) + p3l*pow(numOfBTofMatch,3) + p4l*pow(numOfBTofMatch,4) + p5l*pow(numOfBTofMatch,5);
    double refHigh = p0h + p1h*numOfBTofMatch + p2h*pow(numOfBTofMatch,2) + p3h*pow(numOfBTofMatch,3) + p4h*pow(numOfBTofMatch,4) + p5h*pow(numOfBTofMatch,5); 
   
    if(refMult > refLow && refMult < refHigh) return kFALSE;
  }     

  return kTRUE;
}

bool StEventPlaneCut::passEventCut(StPicoDst *picoDst)
{
  StPicoEvent *picoEvent = picoDst->event();
  if(!picoEvent)
  {
    return kFALSE;
  }

  // const int runId   = picoEvent->runId();
  // const int refMult = picoEvent->refMult();
  // cout << "runId = " << runId << ", refMult = " << refMult << endl;

  // event vertex cut
  const float vx    = picoEvent->primaryVertex().x(); // x works for both TVector3 and StThreeVectorF
  const float vy    = picoEvent->primaryVertex().y();
  const float vz    = picoEvent->primaryVertex().z();
  // const float zdcX  = picoEvent->ZDCx();
  const float vzVpd = picoEvent->vzVpd();
  // vz cut
  if(fabs(vz) > recoEP::mVzMaxMap[mEnergy])
  {
    return kFALSE;
  }
  // vr cut
  if(sqrt(vx*vx+vy*vy) > recoEP::mVrMax[mEnergy] || sqrt(vx*vx+vy*vy) <= recoEP::mVrMin[mEnergy])
  {
    return kFALSE;
  }
  // vz-vzVpd cut only for 200 GeV
  if(fabs(vz-vzVpd) > recoEP::mVzVpdDiffMax[mEnergy])
  {
    return kFALSE;
  }

  return kTRUE;
}

//---------------------------------------------------------------------------------

//---------------------------------------------------------------------------------
bool StEventPlaneCut::passTrackBasic(StPicoTrack *picoTrack)
{
  if(!picoTrack) return kFALSE;

  // nHitsFit cut
  if(picoTrack->nHitsFit() < recoEP::mHitsFitTPCMin)
  {
    return kFALSE;
  }

  // nHitsRatio cut
  if(picoTrack->nHitsMax() <= recoEP::mHitsMaxTPCMin)
  {
    return kFALSE;
  }
  if((double)picoTrack->nHitsFit()/(double)picoTrack->nHitsMax() < recoEP::mHitsRatioTPCMin)
  {
    return kFALSE;
  }

  // eta cut
  // float eta = picoTrack->pMom().pseudoRapidity();
  // float eta = picoTrack->pMom().PseudoRapidity();
  // TVector3 primMom; // temp fix for StThreeVectorF & TVector3
  // const double primPx    = picoTrack->pMom().x(); // x works for both TVector3 and StThreeVectorF
  // const double primPy    = picoTrack->pMom().y();
  // const double primPz    = picoTrack->pMom().z();
  // primMom.SetXYZ(primPx,primPy,primPz);

  // const double eta = primMom.PseudoRapidity();
  // if(fabs(eta) > recoEP::mEtaMax)
  // {
  //   return kFALSE;
  // }

  // if(primMom.Pt() < recoEP::mGlobPtMin) // minimum pT cuts
  // {
  //   return kFALSE;
  // }

  return kTRUE;
}

bool StEventPlaneCut::passTrackEp(StPicoTrack *picoTrack, StPicoEvent* picoEvent)
{
  if(!passTrackBasic(picoTrack)) return kFALSE;

  const double vx = picoEvent->primaryVertex().x(); // x works for both TVector3 and StThreeVectorF
  const double vy = picoEvent->primaryVertex().y();
  const double vz = picoEvent->primaryVertex().z();

  // dca cut for event plane reconstruction: 200GeV = 3.0, BES = 1.0
  if(picoTrack->gDCA(vx,vy,vz) > recoEP::mDcaEPMax[mEnergy])
  {
    return kFALSE;
  }

  TVector3 primMom; // temp fix for StThreeVectorF & TVector3
  const double primPx    = picoTrack->pMom().x(); // x works for both TVector3 and StThreeVectorF
  const double primPy    = picoTrack->pMom().y();
  const double primPz    = picoTrack->pMom().z();
  primMom.SetXYZ(primPx,primPy,primPz);

  // eta cut [-1.0,1.0]
  const double eta = primMom.PseudoRapidity();
  if(fabs(eta) > recoEP::mEtaMax)
  {
    return kFALSE;
  }

  // pt cut 0.2 - 2.0 GeV/c
  const double pt = primMom.Perp();
  const double p  = primMom.Mag();
  if(!(pt > recoEP::mPrimPtMin[mEnergy] && pt < recoEP::mPrimPtMax && p < recoEP::mPrimMomMax))
  {
    return kFALSE;
  }

  return kTRUE;
}

bool StEventPlaneCut::passTrackChargedFlowEast(StPicoTrack *picoTrack, StPicoEvent* picoEvent)
{
  if(!passTrackBasic(picoTrack)) return kFALSE;

  const double vx = picoEvent->primaryVertex().x(); // x works for both TVector3 and StThreeVectorF
  const double vy = picoEvent->primaryVertex().y();
  const double vz = picoEvent->primaryVertex().z();

  // dca cut for charged particle flow: 1.0
  if(picoTrack->gDCA(vx,vy,vz) > recoEP::mDcaTrMax)
  {
    return kFALSE;
  }

  TVector3 primMom; // temp fix for StThreeVectorF & TVector3
  const double primPx    = picoTrack->pMom().x(); // x works for both TVector3 and StThreeVectorF
  const double primPy    = picoTrack->pMom().y();
  const double primPz    = picoTrack->pMom().z();
  primMom.SetXYZ(primPx,primPy,primPz);

  // eta cut for charged particle flow in East EP [-1.0,0.0)
  const double eta = primMom.PseudoRapidity();
  if(!(eta >= -1.0*recoEP::mEtaMax && eta < 0.0))
  {
    return kFALSE;
  }

  // pt and momentum cut: pt > 0.2 GeV/c
  const double pt = primMom.Perp();
  const double p  = primMom.Mag();
  if(!(pt > recoEP::mPrimPtMin[mEnergy] && p < recoEP::mPrimMomMax))
  {
    return kFALSE;
  }

  return kTRUE;
}

bool StEventPlaneCut::passTrackChargedFlowWest(StPicoTrack *picoTrack, StPicoEvent* picoEvent)
{
  if(!passTrackBasic(picoTrack)) return kFALSE;

  const double vx = picoEvent->primaryVertex().x(); // x works for both TVector3 and StThreeVectorF
  const double vy = picoEvent->primaryVertex().y();
  const double vz = picoEvent->primaryVertex().z();

  // dca cut for charged particle flow: 1.0
  if(picoTrack->gDCA(vx,vy,vz) > recoEP::mDcaTrMax)
  {
    return kFALSE;
  }

  TVector3 primMom; // temp fix for StThreeVectorF & TVector3
  const double primPx    = picoTrack->pMom().x(); // x works for both TVector3 and StThreeVectorF
  const double primPy    = picoTrack->pMom().y();
  const double primPz    = picoTrack->pMom().z();
  primMom.SetXYZ(primPx,primPy,primPz);

  // eta cut for charged particle flow in West TPC [0.0,1.0]
  const double eta = primMom.PseudoRapidity();
  if(!(eta >= 0.0 && eta <= recoEP::mEtaMax))
  {
    return kFALSE;
  }

  // pt and momentum cut: pt > 0.2 GeV/c
  const double pt = primMom.Perp();
  const double p  = primMom.Mag();
  if(!(pt > recoEP::mPrimPtMin[mEnergy] && p < recoEP::mPrimMomMax))
  {
    return kFALSE;
  }

  return kTRUE;
}

bool StEventPlaneCut::passTrackChargedFlowFull(StPicoTrack *picoTrack, StPicoEvent* picoEvent)
{
  if(!passTrackBasic(picoTrack)) return kFALSE;

  const double vx = picoEvent->primaryVertex().x(); // x works for both TVector3 and StThreeVectorF
  const double vy = picoEvent->primaryVertex().y();
  const double vz = picoEvent->primaryVertex().z();

  // dca cut for charged particle flow: 1.0
  if(picoTrack->gDCA(vx,vy,vz) > recoEP::mDcaTrMax)
  {
    return kFALSE;
  }

  TVector3 primMom; // temp fix for StThreeVectorF & TVector3
  const double primPx    = picoTrack->pMom().x(); // x works for both TVector3 and StThreeVectorF
  const double primPy    = picoTrack->pMom().y();
  const double primPz    = picoTrack->pMom().z();
  primMom.SetXYZ(primPx,primPy,primPz);

  // eta cut for charged particle flow in Full TPC [-1.0,1.0]
  const double eta = primMom.PseudoRapidity();
  if(fabs(eta) > recoEP::mEtaMax)
  {
    return kFALSE;
  }

  // pt and momentum cut: pt > 0.2 GeV/c
  const double pt = primMom.Perp();
  const double p  = primMom.Mag();
  if(!(pt > recoEP::mPrimPtMin[mEnergy] && p < recoEP::mPrimMomMax))
  {
    return kFALSE;
  }

  return kTRUE;
}
//---------------------------------------------------------------------------------
double StEventPlaneCut::getBeta(StPicoTrack *picoTrack, StPicoDst *picoDst)
{
  double beta = -999.9;
  // StPicoTrack *picoTrack = (StPicoTrack*)picoDst->track(i_track); // return ith track
  int tofIndex = picoTrack->bTofPidTraitsIndex(); // return ToF PID traits
  if(tofIndex >= 0)
  {
    StPicoBTofPidTraits *tofTrack = picoDst->btofPidTraits(tofIndex);
    beta = tofTrack->btofBeta();
  }

  return beta;
}

double StEventPlaneCut::getPrimaryMass2(StPicoTrack *picoTrack, StPicoDst *picoDst)
{
  double mass2 = -999.9;
  // StPicoTrack *picoTrack = (StPicoTrack*)picoDst->track(i_track); // return ith track
  int tofIndex = picoTrack->bTofPidTraitsIndex(); // return ToF PID traits
  if(tofIndex >= 0)
  {
    StPicoBTofPidTraits *tofTrack = picoDst->btofPidTraits(tofIndex);
    const double beta = tofTrack->btofBeta();
    // const double Momentum = picoTrack->pMom().mag(); // primary momentum for 54GeV_2017
    // const double Momentum = picoTrack->pMom().Mag(); // primary momentum for 27GeV_2018
    TVector3 primMom; // temp fix for StThreeVectorF & TVector3
    const double primPx    = picoTrack->pMom().x(); // x works for both TVector3 and StThreeVectorF
    const double primPy    = picoTrack->pMom().y();
    const double primPz    = picoTrack->pMom().z();
    primMom.SetXYZ(primPx,primPy,primPz);
    const double Momentum = primMom.Mag(); // primary momentum

    if(tofTrack->btofMatchFlag() > 0 && tofTrack->btof() != 0 && beta != 0)
    {
      mass2 = Momentum*Momentum*(1.0/(beta*beta) - 1.0);
    }
  }

  return mass2;
}

double StEventPlaneCut::getGlobalMass2(StPicoTrack *picoTrack, StPicoDst *picoDst)
{
  double mass2 = -999.9;
  // StPicoTrack *picoTrack = (StPicoTrack*)picoDst->track(i_track); // return ith track
  int tofIndex = picoTrack->bTofPidTraitsIndex(); // return ToF PID traits
  if(tofIndex >= 0)
  {
    StPicoBTofPidTraits *tofTrack = picoDst->btofPidTraits(tofIndex);
    const double beta = tofTrack->btofBeta();
    // const double Momentum = picoTrack->gMom().mag(); // global momentum for 54GeV_2017
    // const double Momentum = picoTrack->gMom().Mag(); // global momentum for 27GeV_2018
    TVector3 globMom; // temp fix for StThreeVectorF & TVector3
    const double globPx     = picoTrack->gMom().x(); // x works for both TVector3 and StThreeVectorF
    const double globPy     = picoTrack->gMom().y();
    const double globPz     = picoTrack->gMom().z();
    globMom.SetXYZ(globPx,globPy,globPz);
    const double Momentum = globMom.Mag(); // global momentum

    if(tofTrack->btofMatchFlag() > 0 && tofTrack->btof() != 0 && beta != 0)
    {
      mass2 = Momentum*Momentum*(1.0/(beta*beta) - 1.0);
    }
  }

  return mass2;
}
//---------------------------------------------------------------------------------
