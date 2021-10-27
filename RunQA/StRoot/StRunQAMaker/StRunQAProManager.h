#ifndef StRunQAProManager_h
#define StRunQAProManager_h

#include <TVector3.h>
#include <TString.h>
#include "StMessMgr.h"

class TProfile;

class StRunQAProManager
{
  public:
    StRunQAProManager();
    virtual ~StRunQAProManager();

    // Run-by-Run QA
    void initRunQA();
    void fillRunQA_Event(int triggerBin, int runIdenx, float refMult, float grefMult, float zdcX, float vx, float vy, float vz, int cutSelection);
    void fillRunQA_Track(int triggerBin, int runIdenx, float gDca, int nHitsFit, TVector3 pMom, TVector3 gMom, float dEdx, float beta, float mass2, int nHitsMax, int nHitsDEdx, int cutSelection);
    void fillRunQA_Track_PID(int triggerBin, int runIdenx, float nSigmaPi, float nSigmaK, float nSigmaP, float nSigmaE, float mass2, int cutSelection);
    void writeRunQA();

  private:
    // Run-by-Run QA | x axis is RunIndex
    TProfile *p_mRefMult[2][10]; // 0: before cuts | 1: after cuts
    TProfile *p_mGRefMult[2][10]; // 0-8 for different triggerID | 9 for all triggers
    TProfile *p_mZdcX[2][10];
    TProfile *p_mVz[2][10];
    TProfile *p_mVr[2][10];
    TProfile *p_mGDca[2][10];
    TProfile *p_mNHitsFit[2][10];
    TProfile *p_mPrimPt[2][10];
    TProfile *p_mPrimEta[2][10];
    TProfile *p_mPrimPhi[2][10];
    TProfile *p_mGlobPt[2][10];
    TProfile *p_mGlobEta[2][10];
    TProfile *p_mGlobPhi[2][10];
    TProfile *p_mDEdx[2][10];
    TProfile *p_mBetaInv[2][10];
    TProfile *p_mPrimaryMass2[2][10];
    TProfile *p_mNHitsMax[2][10];
    TProfile *p_mNHitsDEdx[2][10];
    TProfile *p_mNHitsRatio[2][10];
    TProfile *p_mNSigmaPi[2][10];
    TProfile *p_mNSigmaK[2][10];
    TProfile *p_mNSigmaP[2][10];
    TProfile *p_mNSigmaE[2][10];    
    TProfile *p_mMass2Pi[2][10];
    TProfile *p_mMass2K[2][10];
    TProfile *p_mMass2P[2][10];
    TProfile *p_mMass2E[2][10];
 

    std::string mCutsQA[2] = {"Before","After"};

    ClassDef(StRunQAProManager,1)
};

#endif
