#ifndef StEventPlaneCons_h
#define StEventPlaneCons_h

#include <string>
#include "TString.h"
// #include "StarClassLibrary/SystemOfUnits.h"

namespace recoEP
{
  //--------------------------------------------------
  // used in Event Plane Reconstruction
  const int NumBeamEnergy = 2;
  const std::string mBeamEnergy[NumBeamEnergy] = {"14p5GeV_2019","19p6GeV_2019"};
  const double mEnergyValue[NumBeamEnergy] = {14.5,19.6};
  const int mBeamYear[NumBeamEnergy] = {2019,2019};

  // event cut
  const double mVzMaxMap[NumBeamEnergy] = {70.0,70.0}; // 0: 200GeV_2014 | 1: 54GeV_2017 | 2: 27GeV_2018 
  const double mVrMax[NumBeamEnergy] = {2.0,2.0};
  const double mVrMin[NumBeamEnergy] = {0.0,0.0};
  const double mVzVpdDiffMax[NumBeamEnergy] = {10.0,10.0}; // 3.0
  const int mMatchedToFMin[NumBeamEnergy] = {2,2}; // 2

  // track cut
  const double mSigScaleMap[NumBeamEnergy] = {1.0,1.0};
  const double mDcaEPMax[NumBeamEnergy] = {3.0,3.0}; // for event plane reconstruction: 1.0 for BES
  const double mDcaTrMax = 1.0; // for pion, kaon, proton mDcaTrMax = 1.0 for flow
  const int mHitsDedxMin = 5;
  const int mHitsFitTPCMin = 15;
  const int mHitsMaxTPCMin = 0;
  const double mHitsRatioTPCMin = 0.52;
  const double mEtaMax = 1.0;
  const double mEtaGap = 0.05;
  const double mPrimPtMin[NumBeamEnergy] = {0.15,0.15}; // for event plane reconstruction and for pion, kaon, proton: 0.2 for BES
  const double mPrimPtMax = 2.0;
  const double mPrimPtWeight = 2.0;
  const double mPrimMomMax = 10.0; // also use for gMom
  const double mGlobPtMin = 0.1; // for phi, Lambda, K0s
  const int mTrackMin = 2;
  const int mTrackMinFull = 4;

  const int mNumOfRunIndex = 2000;

  // ZDC-SMD Event Plane
  const std::string mEastWest[2] = {"East","West"};
  const std::string mVertHori[2] = {"Vertical","Horizontal"};

  const std::string mVStr[2] = {"pos","neg"};
  const std::string mOrder = "2nd";
  const int mNumShiftOrder = 20;

  // EPD Event Plane
  const int mEpdEpOrder = 3;
}

#endif
