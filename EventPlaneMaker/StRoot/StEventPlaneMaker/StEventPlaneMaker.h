#ifndef StEventPlaneMaker_h
#define StEventPlaneMaker_h

#include "StMaker.h"
// #include "TString.h"
#include "TClonesArray.h"
#include "TChain.h"
#include <string>

class StPicoDst;
class StPicoDstMaker;
class StPicoEvent;
class StRefMultCorr;
class StEpdEpFinder;
class StEpdEpInfo;

class StEventPlaneUtility;
class StEventPlaneCut;
class StEventPlaneHistoManager;
class StEventPlaneProManager;
class StZdcEpManager;
class StTpcEpManager;


class StEventPlaneMaker : public StMaker {
  public:
    StEventPlaneMaker(const char *name, StPicoDstMaker *picoMaker, const string jobId, const char* inputDir, const int Mode, const int inputEpdMode, const int EpdMode, const int Energy, const float mipThresh, const float maxTile);
    virtual ~StEventPlaneMaker();
    
    virtual int Init();
    virtual int Make();
    virtual void Clear(Option_t *opt="");
    virtual int Finish();
    
  private:
    StPicoDstMaker *mPicoDstMaker;
    StPicoDst      *mPicoDst;
    StPicoEvent    *mPicoEvent;
    StRefMultCorr  *mRefMultCorr;
    TChain         *mPicoDstChain;

    StEventPlaneUtility *mEventPlaneUtility;
    StEventPlaneCut  *mEventPlaneCut;
    StEventPlaneHistoManager *mEventPlaneHistoManager;
    StEventPlaneProManager *mEventPlaneProManager;
    StZdcEpManager *mZdcEpManager;
    //StBbcEpManager *mBbcEpManager;
    //StEpdEpManager *mEpdEpManager;
    StTpcEpManager *mTpcEpManager;
    StEpdEpFinder *mEpFinder;
    
    TClonesArray *mEpdHits;    
    TClonesArray *mTracks;
    TClonesArray *mEventClonesArray; 

    int mMode; 
    int mEnergy;
    int mEpdMode; // 0 = Phi-weighting corrections | 1 = Psi-shift corrections | 2 = Run with correct phi-weighting and psi-shift corrections
    int mInputEpdMode;
    float mMipThresh;
    float mMaxTile; 

    string mOutPut_GainCorr;
    string mOutPut_ReCenterPar;
    string mOutPut_ShiftPar;
    string mOutPut_ShiftParFull;
    string mOutPut_Resolution;
    string mOutPut_ChargedFlow;
    string mOutPut_EpdResults;
    string mInPut_EpdCorrections;
    string mOutPut_EpdCorrections; // output will be autogenerated by StEpdEpFinder

    TFile *mFile_GainCorr;
    TFile *mFile_ReCenterPar;
    TFile *mFile_ShiftPar;
    TFile *mFile_ShiftParFull;
    TFile *mFile_Resolution;
    TFile *mFile_ChargedFlow;
    TFile *mFile_EpdResults;
    

    int mUsedTrackCounter;

    ClassDef(StEventPlaneMaker, 1)
};

#endif
