#ifndef StRunQACons_h
#define StRunQACons_h

#include <string>
#include "TString.h"
// #include "StarClassLibrary/SystemOfUnits.h"

namespace runQA
{
  //--------------------------------------------------
  // used in RunQA
  int const NumBeamEnergy = 2;
  std::string const mBeamEnergy[NumBeamEnergy] = {"14p5GeV_2019","19p6GeV_2019"};
  float const mEnergyValue[NumBeamEnergy] = {14.5,19.6};
  int const mBeamYear[NumBeamEnergy] = {2019,2019};

  // event cut
  float const mVzMaxMap[NumBeamEnergy] = {70.0,70.0}; // 0: 19p6GeV_2019 | 1: 54GeV_2017  
  float const mVrMax[NumBeamEnergy] = {2.0,2.0};
  float const mVrMin[NumBeamEnergy] = {0.0,0.0};
  float const mVzVpdDiffMax[NumBeamEnergy] = {100.0,100.0};
  int const mMatchedToFMin[NumBeamEnergy] = {0,0}; 

  // track cut
  float const mSigScaleMap[NumBeamEnergy] = {1.0,1.0};
  float const mDcaTrQAMax = 3.0; // use primary tracks run-by-run QA 
  int const mHitsDedxMin = 5;
  int const mHitsFitTPCMin = 15;
  int const mHitsMaxTPCMin = 0;
  float const mHitsRatioTPCMin = 0.52;
  float const mEtaMax = 1.0;
  float const mPrimPtMin[NumBeamEnergy] = {0.15,0.15}; // for event plane reconstruction and for pion, kaon, proton: 0.2 for BES
  float const mPrimPtMax = 2.0;
  float const mPrimPtWeight = 2.0;
  float const mPrimMomMax = 10.0; // also use for gMom
  float const mGlobPtMin = 0.1; // for phi, Lambda, K0s

  int const mNumOfRunIndex = 2000;
}

#endif
