#ifndef StEventPlaneProManager_h
#define StEventPlaneProManager_h

#include <TString.h>
#include <TVector2.h>
#include "StMessMgr.h"

class TProfile;
class TProfile2D;

class StEventPlaneProManager
{
  public:
    StEventPlaneProManager();
    virtual ~StEventPlaneProManager();

    //--------------ZDC EP---------------
    void initZdcReCenter(); // Re-Center
    void fillZdcReCenterEast(TVector2 qVector, int Cent9, int RunIndex, int VzSign); // VzSign = vertex pos/neg
    void fillZdcReCenterWest(TVector2 qVector, int Cent9, int RunIndex, int VzSign);
    void writeZdcReCenter();

    void initZdcShift(); // Shift
    void fillZdcShiftEast(TVector2 qVector, int Cent9, int RunIndex, int VzSign); // VzSign = vertex pos/neg
    void fillZdcShiftWest(TVector2 qVector, int Cent9, int RunIndex, int VzSign);
    void writeZdcShift();

    void initZdcShiftFull();
    void fillZdcShiftFull(TVector2 qVector, int Cent9, int RunIndex, int VzSign);
    void writeZdcShiftFull();

    void initZdcResolution();
    void fillZdcResSub(TVector2 QEast, TVector2 QWest, int Cent9, int RunIndex);
    void writeZdcResolution();
    //--------------ZDC EP---------------
    
    //--------------TPC EP---------------
    void initTpcReCenter(); // Re-Center
    void fillTpcReCenterEast(TVector2 qVector, int Cent9, int RunIndex, int VzSign, double pt); // VzSign = vertex pos/neg
    void fillTpcReCenterWest(TVector2 qVector, int Cent9, int RunIndex, int VzSign, double pt);
    void fillTpcReCenterFull(TVector2 qVector, int Cent9, int RunIndex, int VzSign, double pt);
    void writeTpcReCenter();

    void initTpcShift(); // Re-Center
    void fillTpcShiftEast(TVector2 qVector, int Cent9, int RunIndex, int VzSign); // VzSign = vertex pos/neg
    void fillTpcShiftWest(TVector2 qVector, int Cent9, int RunIndex, int VzSign);
    void fillTpcShiftFull(TVector2 qVector, int Cent9, int RunIndex, int VzSign);
    void writeTpcShift();

    void initTpcResolution();
    void fillTpcResSub(double PsiEast, double PsiWest, int Cent9, int RunIndex);
    void fillTpcResRan(double PsiRanA, double PsiRanB, int Cent9, int RunIndex);
    void writeTpcResolution();
    //--------------TPC EP---------------

  private:
    //--------------ZDC EP---------------
    // ZDC-SMD ReCenter Correction | x axis is RunIndex, y axis is Centrality
    TProfile2D *p_mZdcQ1EastVertical[2]; // 0 = vertex pos/neg
    TProfile2D *p_mZdcQ1EastHorizontal[2];
    TProfile2D *p_mZdcQ1WestVertical[2];
    TProfile2D *p_mZdcQ1WestHorizontal[2];

    // ZDC-SMD Shift Correction | x axis is RunIndex, y axis is Centrality
    TProfile2D *p_mZdcQ1EastCos[2][20]; // 0 = vertex pos/neg | 1 = shift correction harmonics
    TProfile2D *p_mZdcQ1EastSin[2][20];
    TProfile2D *p_mZdcQ1WestCos[2][20];
    TProfile2D *p_mZdcQ1WestSin[2][20];
    TProfile2D *p_mZdcQ1FullCos[2][20];
    TProfile2D *p_mZdcQ1FullSin[2][20];

    // EP Resolution
    TProfile2D *p_mZdcSubRes1QA; // 1st Res vs runIndex & cent9
    TProfile2D *p_mZdcSubRes2QA; // 2nd Res vs runIndex & cent9
    TProfile *p_mZdcSubRes1; // 1st Res vs cent9
    TProfile *p_mZdcSubRes2; // 2nd Res vs cent9
    //--------------ZDC EP---------------

    //--------------TPC EP---------------
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

    // EP Resolution
    TProfile2D *p_mTpcSubRes2QA; // 2nd Res vs runIndex & cent9
    TProfile2D *p_mTpcRanRes2QA; // 2nd Res vs runIndex & cent9
    TProfile *p_mTpcSubRes2; // 2nd Res vs cent9
    TProfile *p_mTpcRanRes2; // 2nd Res vs cent9
    //--------------TPC EP---------------

    ClassDef(StEventPlaneProManager,1)
};

#endif
