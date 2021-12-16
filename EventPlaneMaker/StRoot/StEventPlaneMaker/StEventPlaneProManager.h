#ifndef StEventPlaneProManager_h
#define StEventPlaneProManager_h

#include "StEventPlaneCons.h"
#include "StEpdUtil/StEpdEpInfo.h"
#include <TString.h>
#include <TVector2.h>
#include "StMessMgr.h"

class TProfile;
class TProfile2D;
class StEpdEpInfo;

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
    void fillTpcReCenterEast(int order, TVector2 qVector, int Cent9, int RunIndex, int VzSign, double pt); // VzSign = vertex pos/neg
    void fillTpcReCenterWest(int order, TVector2 qVector, int Cent9, int RunIndex, int VzSign, double pt);
    void fillTpcReCenterFull(int order, TVector2 qVector, int Cent9, int RunIndex, int VzSign, double pt);
    void writeTpcReCenter();

    void initTpcShift(); // Re-Center
    void fillTpcShiftEast(int order, TVector2 qVector, int Cent9, int RunIndex, int VzSign); // VzSign = vertex pos/neg
    void fillTpcShiftWest(int order, TVector2 qVector, int Cent9, int RunIndex, int VzSign);
    void fillTpcShiftFull(int order, TVector2 qVector, int Cent9, int RunIndex, int VzSign);
    void writeTpcShift();

    void initTpcResolution();
    void fillTpcResSub(int order, double PsiEast, double PsiWest, int Cent9, int RunIndex);
    void fillTpcResRan(int order, double PsiRanA, double PsiRanB, int Cent9, int RunIndex);
    void writeTpcResolution();
    //--------------TPC EP---------------

    //--------------EPD EP---------------
    void initEpdRes();
    void fillEpdRes(StEpdEpInfo result, int Cent9, int runIndex);
    void writeEpdRes();

    void initEpdFlowEta();
    void fillEpdFlowEta(double eta, double v, int Cent9, int order, double weight);
    void writeEpdFlowEta();

    void initEpdFlow();
    void fillEpdFlow(double v, int Cent9, int order, int runIndex);
    void writeEpdFlow();
   
    //--------------EPD EP---------------

    //--------------Charged Flow---------------
    void initChargedFlow();
    void fillChargedV1Pp(double pt, double eta, double v1, double res1, int Cent9, int RunIndex, double reweight);
    void fillChargedV1Ep(double pt, double eta, double v1, double res1, int Cent9, int RunIndex, double reweight);
    void fillChargedV1Epd(double pt, double v1, double res1, int Cent9, int RunIndex);
    void fillChargedV2Pp(double pt, double v2, double res2, int Cent9, int RunIndex, double reweight);
    void fillChargedV2Ep(double pt, double eta, double v2, double res2, int Cent9, int RunIndex, double reweight);
    void fillChargedV2Epd(double pt, double v2, double res2, int Cent9, int RunIndex);
    void fillChargedV3Ep(double pt, double eta, double v3, double res3, int Cent9, int RunIndex, double reweight);
    void fillChargedV3Epd(double pt, double v3, double res3, int Cent9, int RunIndex);
    void writeChargedFlow();
    //--------------Charged Flow---------------

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

    // EP Resolution
    TProfile2D *p_mTpcSubRes1QA; // 1st Res vs runIndex & cent9
    TProfile2D *p_mTpcRanRes1QA; // 1st Res vs runIndex & cent9
    TProfile *p_mTpcSubRes1; // 1st Res vs cent9
    TProfile *p_mTpcRanRes1; // 1st Res vs cent9
    
    TProfile2D *p_mTpcSubRes2QA; // 2nd Res vs runIndex & cent9
    TProfile2D *p_mTpcRanRes2QA; // 2nd Res vs runIndex & cent9
    TProfile *p_mTpcSubRes2; // 2nd Res vs cent9
    TProfile *p_mTpcRanRes2; // 2nd Res vs cent9
    
    TProfile2D *p_mTpcSubRes3QA; // 3rd Res vs runIndex & cent9
    TProfile2D *p_mTpcRanRes3QA; // 3rd Res vs runIndex & cent9
    TProfile *p_mTpcSubRes3; // 3rd Res vs cent9
    TProfile *p_mTpcRanRes3; // 3rd Res vs cent9

    //--------------TPC EP---------------

    //--------------EPD EP---------------
    TProfile2D *p_mEpdSubResQA[recoEP::mEpdEpOrder]; // 1st Res vs runIndex & cent9
    TProfile *p_mEpdSubRes[recoEP::mEpdEpOrder]; // 1st Res vs cent9

    TProfile *p_mEpdAveCos[recoEP::mEpdEpOrder][3];

    TProfile *p_mEpdFlowEta[recoEP::mEpdEpOrder][10];
    TProfile *p_mEpdFlowEtaWeights[recoEP::mEpdEpOrder][10];
    TProfile *p_mEpdFlow[recoEP::mEpdEpOrder];
    TProfile *p_mEpdFlowQA[recoEP::mEpdEpOrder][10];
    //--------------EPD EP---------------

    //--------------Charged Flow---------------
    // charged particle flow for different centrality bins: 0-8 cent9, 9 minBias
    TProfile *p_mChargedV1PpQA[10]; // <v1Pp> vs. runIndex | pt [0.2, 2.0]
    TProfile *p_mChargedV1EpQA[10];
    TProfile *p_mChargedV1EpdQA[10];
    TProfile *p_mChargedV2EpQA[10]; // <v2Ep> vs. runIndex | pt [0.2, 5.2]
    TProfile *p_mChargedV2PpQA[10]; // <v2Pp> vs. runIndex | pt [0.2, 5.2]
    TProfile *p_mChargedV2EpdQA[10];
    TProfile *p_mChargedV3EpQA[10]; // <v3Ep> vs. runIndex | pt [0.2, 5.2]
    TProfile *p_mChargedV3EpdQA[10];
    TProfile *p_mChargedV1Pp[10]; // v1Pp vs. eta
    TProfile *p_mChargedV1Ep[10];
    TProfile *p_mChargedV1Epd[10];
    TProfile *p_mChargedV2Ep[10]; // v2Ep vs. pt
    TProfile *p_mChargedV2Pp[10]; // v2Pp vs. pt
    TProfile *p_mChargedV2Epd[10];
    TProfile *p_mChargedV3Ep[10]; // v3Ep vs. pt
    TProfile *p_mChargedV3Epd[10];
    TProfile *p_mChargedV1EpEta[10];
    TProfile *p_mChargedV2EpEta[10];
    TProfile *p_mChargedV3EpEta[10];


    //--------------Charged Flow---------------

    ClassDef(StEventPlaneProManager,1)
};

#endif
