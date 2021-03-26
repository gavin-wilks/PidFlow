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
    bool passTrackEtaEast(StPicoTrack*);
    bool passTrackEtaWest(StPicoTrack*);
    bool passTrackEtaFull(StPicoTrack*);

    TVector2 calq2Vector(StPicoTrack*);
    double getWeight(StPicoTrack*);

    void addTrackEastRaw(StPicoTrack* picoTrack); // raw EP
    void addTrackWestRaw(StPicoTrack* picoTrack);
    void addTrackFullRaw(StPicoTrack* picoTrack);

    void initReCenterCorr();
    TVector2 getReCenterParEast();
    TVector2 getReCenterParWest();
    TVector2 getReCenterParFull();

    void addTrackEast(StPicoTrack* picoTrack); // re-centered EP
    void addTrackWest(StPicoTrack* picoTrack);
    void addTrackFull(StPicoTrack* picoTrack);
    void addTrackSubA(StPicoTrack* picoTrack); // random sub A
    void addTrackSubB(StPicoTrack* picoTrack); // random sub B
    void Randomization();

    void print(TVector2);

    // Shift Correction
    bool passTrackEtaNumRawCut();
    bool passTrackFullNumRawCut();
    bool passTrackEtaNumCut();
    bool passTrackFullNumCut();

    void initShiftCorr();
    double AngleShift(double Psi2Raw);

    double calShiftAngle2East();
    double calShiftAngle2West();
    double calShiftAngle2SubA();
    double calShiftAngle2SubB();
    double calShiftAngle2Full();
    double calShiftAngle2Full(StPicoTrack *picoTrack); // subtract self-correlation

    void initResolutionCorr();
    double getRes2Sub(int Cent9);
    double getRes2Full(int Cent9);

    TVector2 getQVector(int nEP); // east/west/full/subA/subB
    TVector2 getQVectorRaw(int nEP);
    int getNumTrack(int nEP);

  private:
    //Event Plane method
    TVector2 mQ2VecEastRaw, mQ2VecWestRaw, mQ2VecFullRaw;
    int mQ2CounterRawEast, mQ2CounterRawWest, mQ2CounterRawFull;

    TVector2 mQ2VecEast, mQ2VecWest, mQ2VecFull, mQ2VecSubA, mQ2VecSubB;
    int mQ2CounterEast, mQ2CounterWest, mQ2CounterFull;
    int mQ2CounterFull_East, mQ2CounterFull_West;
    int mQ2CounterSubA, mQ2CounterSubB;

    int mEnergy;
    int mCent9;
    int mRunIndex;
    int mVzSign;

    // TPC ReCenter Correction | x axis is RunIndex, y axis is Centrality
    TProfile2D *p_mTpcq2xEast[2]; // 0 = vertex pos/neg
    TProfile2D *p_mTpcq2yEast[2];
    TProfile2D *p_mTpcq2xWest[2];
    TProfile2D *p_mTpcq2yWest[2];
    TProfile2D *p_mTpcq2xFull[2];
    TProfile2D *p_mTpcq2yFull[2];

    // Shift Correction | x axis is RunIndex, y axis is Centrality
    TProfile2D *p_mTpcQ2EastCos[2][20]; // 0 = vertex pos/neg, 1 = ShiftOrder
    TProfile2D *p_mTpcQ2EastSin[2][20];
    TProfile2D *p_mTpcQ2WestCos[2][20]; 
    TProfile2D *p_mTpcQ2WestSin[2][20];
    TProfile2D *p_mTpcQ2FullCos[2][20];
    TProfile2D *p_mTpcQ2FullSin[2][20];

    TFile *mInPutFile_ReCenter;
    TFile *mInPutFile_Shift;
    TFile *mInPutFile_Res;

  ClassDef(StTpcEpManager,1)
};

#endif
