#ifndef StTpcEpManager_h
#define StTpcEpManager_h

#include "TObject.h"
#include "TVector2.h"

class StPicoDst;
class StPicoTrack;
class TProfile2D;
class TFile;
class TNtuple;

class StTpcEpManager : public TObject
{
  public:
    StTpcEpManager(int energy);
    virtual ~StTpcEpManager();
    void clearTpcEp();
    void initTpcEp(int Cent9, int RunIndex, int VzSign);

    // ReCenter Correction
    bool passTrackEpEast(StPicoTrack*);
    bool passTrackEpWest(StPicoTrack*);
    bool passTrackEpFull(StPicoTrack*);
    
    TVector2 calq1Vector(StPicoTrack*);
    TVector2 calq2Vector(StPicoTrack*);
    TVector2 calq3Vector(StPicoTrack*);

    double getWeight(StPicoTrack*);

    void addTrackEastRaw(StPicoTrack* picoTrack); // raw EP
    void addTrackWestRaw(StPicoTrack* picoTrack);
    void addTrackFullRaw(StPicoTrack* picoTrack);

    void readReCenterCorr();
    TVector2 getReCenterParEast(int order);
    TVector2 getReCenterParWest(int order);
    TVector2 getReCenterParFull(int order);

    void addTrackEast(StPicoTrack* picoTrack); // re-centered EP
    void addTrackWest(StPicoTrack* picoTrack);
    void addTrackFull(StPicoTrack* picoTrack);
    void addTrackRanA(StPicoTrack* picoTrack); // random sub A
    void addTrackRanB(StPicoTrack* picoTrack); // random sub B
    void Randomization();

    void print(TVector2);

    // Shift Correction
    bool passTrackEtaNumRawCut(int order);
    bool passTrackFullNumRawCut(int order);
    bool passTrackEtaNumCut(int order);
    bool passTrackFullNumCut(int order);

    void readShiftCorr();
    double AngleShift(double PsiRaw, int order);

    double calShiftAngle1East();
    double calShiftAngle1West();
    double calShiftAngle1RanA();
    double calShiftAngle1RanB();
    double calShiftAngle1Full();
    double calShiftAngle1Full(StPicoTrack *picoTrack); // subtract self-correlation
    
    double calShiftAngle2East();
    double calShiftAngle2West();
    double calShiftAngle2RanA();
    double calShiftAngle2RanB();
    double calShiftAngle2Full();
    double calShiftAngle2Full(StPicoTrack *picoTrack); // subtract self-correlation
    
    double calShiftAngle3East();
    double calShiftAngle3West();
    double calShiftAngle3RanA();
    double calShiftAngle3RanB();
    double calShiftAngle3Full();
    double calShiftAngle3Full(StPicoTrack *picoTrack); // subtract self-correlation



    void readResolution();
    double getRes1Sub(int Cent9);
    double getRes1SubErr(int Cent9);
    double getRes1Full(int Cent9);
    double getRes1FullErr(int Cent9);
 
    double getRes2Sub(int Cent9);
    double getRes2SubErr(int Cent9);
    double getRes2Full(int Cent9);
    double getRes2FullErr(int Cent9);

    double getRes3Sub(int Cent9);
    double getRes3SubErr(int Cent9);
    double getRes3Full(int Cent9);
    double getRes3FullErr(int Cent9);


    TVector2 getQVector(int nEP, int order); // east/west/full/subA/subB
    TVector2 getQVectorRaw(int nEP, int order);
    int getNumTrack(int nEP, int order);

  private:
    //Event Plane method
    TVector2 mQ1VecEastRaw, mQ1VecWestRaw, mQ1VecFullRaw;
    int mQ1CounterRawEast, mQ1CounterRawWest, mQ1CounterRawFull;   
    int mQ1CounterRawFull_East, mQ1CounterRawFull_West;

    TVector2 mQ2VecEastRaw, mQ2VecWestRaw, mQ2VecFullRaw;
    int mQ2CounterRawEast, mQ2CounterRawWest, mQ2CounterRawFull;   
    int mQ2CounterRawFull_East, mQ2CounterRawFull_West; 

    TVector2 mQ3VecEastRaw, mQ3VecWestRaw, mQ3VecFullRaw;    
    int mQ3CounterRawEast, mQ3CounterRawWest, mQ3CounterRawFull; 
    int mQ3CounterRawFull_East, mQ3CounterRawFull_West;

    TVector2 mQ1VecEast, mQ1VecWest, mQ1VecFull, mQ1VecRanA, mQ1VecRanB;
    int mQ1CounterEast, mQ1CounterWest, mQ1CounterFull;
    int mQ1CounterFull_East, mQ1CounterFull_West;
    int mQ1CounterRanA, mQ1CounterRanB;

    TVector2 mQ2VecEast, mQ2VecWest, mQ2VecFull, mQ2VecRanA, mQ2VecRanB;
    int mQ2CounterEast, mQ2CounterWest, mQ2CounterFull;
    int mQ2CounterFull_East, mQ2CounterFull_West;
    int mQ2CounterRanA, mQ2CounterRanB;

    TVector2 mQ3VecEast, mQ3VecWest, mQ3VecFull, mQ3VecRanA, mQ3VecRanB;
    int mQ3CounterEast, mQ3CounterWest, mQ3CounterFull;
    int mQ3CounterFull_East, mQ3CounterFull_West;
    int mQ3CounterRanA, mQ3CounterRanB;


    int mEnergy;
    int mCent9;
    int mRunIndex;
    int mVzSign;

    // EP resolution
    double mTpcSubRes1Val[9];
    double mTpcSubRes1Err[9];
    double mTpcFullRes1Val[9];
    double mTpcFullRes1Err[9];
    
    double mTpcSubRes2Val[9];
    double mTpcSubRes2Err[9];
    double mTpcFullRes2Val[9];
    double mTpcFullRes2Err[9];

    double mTpcSubRes3Val[9];
    double mTpcSubRes3Err[9];
    double mTpcFullRes3Val[9];
    double mTpcFullRes3Err[9];


    // TPC ReCenter Correction | x axis is RunIndex, y axis is Centrality
    TProfile2D *p_mTpcq1xEast[2]; // 0 = vertex pos/neg
    TProfile2D *p_mTpcq1yEast[2];
    TProfile2D *p_mTpcq1xWest[2];
    TProfile2D *p_mTpcq1yWest[2];
    TProfile2D *p_mTpcq1xFull[2];
    TProfile2D *p_mTpcq1yFull[2];

    TProfile2D *p_mTpcq2xEast[2]; // 0 = vertex pos/neg
    TProfile2D *p_mTpcq2yEast[2];
    TProfile2D *p_mTpcq2xWest[2];
    TProfile2D *p_mTpcq2yWest[2];
    TProfile2D *p_mTpcq2xFull[2];
    TProfile2D *p_mTpcq2yFull[2];

    TProfile2D *p_mTpcq3xEast[2]; // 0 = vertex pos/neg
    TProfile2D *p_mTpcq3yEast[2];
    TProfile2D *p_mTpcq3xWest[2];
    TProfile2D *p_mTpcq3yWest[2];
    TProfile2D *p_mTpcq3xFull[2];
    TProfile2D *p_mTpcq3yFull[2];

    // Shift Correction | x axis is RunIndex, y axis is Centrality
    TProfile2D *p_mTpcQ1EastCos[2][20]; // 0 = vertex pos/neg, 1 = ShiftOrder
    TProfile2D *p_mTpcQ1EastSin[2][20];
    TProfile2D *p_mTpcQ1WestCos[2][20]; 
    TProfile2D *p_mTpcQ1WestSin[2][20];
    TProfile2D *p_mTpcQ1FullCos[2][20];
    TProfile2D *p_mTpcQ1FullSin[2][20];

    TProfile2D *p_mTpcQ2EastCos[2][20]; // 0 = vertex pos/neg, 1 = ShiftOrder
    TProfile2D *p_mTpcQ2EastSin[2][20];
    TProfile2D *p_mTpcQ2WestCos[2][20]; 
    TProfile2D *p_mTpcQ2WestSin[2][20];
    TProfile2D *p_mTpcQ2FullCos[2][20];
    TProfile2D *p_mTpcQ2FullSin[2][20];

    TProfile2D *p_mTpcQ3EastCos[2][20]; // 0 = vertex pos/neg, 1 = ShiftOrder
    TProfile2D *p_mTpcQ3EastSin[2][20];
    TProfile2D *p_mTpcQ3WestCos[2][20]; 
    TProfile2D *p_mTpcQ3WestSin[2][20];
    TProfile2D *p_mTpcQ3FullCos[2][20];
    TProfile2D *p_mTpcQ3FullSin[2][20];

    TFile *mInPutFile_ReCenter;
    TFile *mInPutFile_Shift;
    TFile *mInPutFile_Res;

  ClassDef(StTpcEpManager,1)
};

#endif
