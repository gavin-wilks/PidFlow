#ifndef StEventPlaneHistoManager_h
#define StEventPlaneHistoManager_h

#include "StMessMgr.h"
#include "TVector2.h"

class TH1F;
class TH2F;

class StEventPlaneHistoManager
{
  public:
    StEventPlaneHistoManager();
    virtual ~StEventPlaneHistoManager();

    //--------------ZDC EP---------------
    void initZdcGainCorr();
    void fillZdcGainCorr(int i_eastwest, int i_verthori, int i_slat, int runIndex, float zdcsmd);
    void writeZdcGainCorr();

    void initZdcRawEP(); // raw ZDC-SMD EP
    void fillZdcRawSubEP(TVector2 QEast, TVector2 QWest, int Cent9, int runIndex);
    void fillZdcRawFullEP(TVector2 QFull, int Cent9, int runIndex);
    void writeZdcRawEP();

    void initZdcReCenterEP(); // recenter ZDC-SMD EP
    void fillZdcReCenterSubEP(TVector2 QEast, TVector2 QWest, int Cent9, int runIndex);
    void fillZdcReCenterFullEP(TVector2 QFull, int Cent9, int runIndex);
    void writeZdcReCenterEP();

    void initZdcShiftEP(); // recenter ZDC-SMD EP
    void fillZdcShiftSubEP(TVector2 QEast, TVector2 QWest, int Cent9, int runIndex);
    void fillZdcShiftFullEP(TVector2 QDiff, TVector2 QFull, int Cent9, int runIndex);
    void writeZdcShiftEP();
    //--------------ZDC EP---------------
    
    //--------------TPC EP---------------
    void initTpcRawEP(); // raw TPC EP
    void fillTpcRawSubEP(TVector2 QEast, TVector2 QWest, int Cent9, int runIndex);
    void fillTpcRawFullEP(TVector2 QFull, int Cent9, int runIndex);
    void writeTpcRawEP();

    void initTpcReCenterEP(); // recenter TPC EP
    void fillTpcReCenterSubEP(TVector2 QEast, TVector2 QWest, int Cent9, int runIndex);
    void fillTpcReCenterFullEP(TVector2 QFull, int Cent9, int runIndex);
    void writeTpcReCenterEP();

    void initTpcShiftEP(); // recenter TPC EP
    void fillTpcShiftSubEP(double PsiEast, double PsiWest, int Cent9, int runIndex);
    void fillTpcShiftRanEP(double PsiRanA, double PsiRanB, int Cent9, int runIndex);
    void fillTpcShiftFullEP(double PsiFull, int Cent9, int runIndex);
    void writeTpcShiftEP();
    //--------------TPC EP---------------

  private:
    //--------------ZDC EP---------------
    TH2F *h_mZdcGainCorr[2][2][8]; // 0: east/west | 1: vertical(x)/horizontal(y) | 2: 7 slats(x)/8 slats(y); | x-axis: runIndex | y-axis: ADC

    // x-axis: runIndex | y-axis: EP
    TH2F *h_mZdcRawEpEast[9]; // raw EP
    TH2F *h_mZdcRawEpWest[9];
    TH2F *h_mZdcRawEpFull[9]; // Qwest-QEast

    TH2F *h_mZdcReCenterEpEast[9]; // recenter EP
    TH2F *h_mZdcReCenterEpWest[9];
    TH2F *h_mZdcReCenterEpFull[9]; // Qwest-QEast

    TH2F *h_mZdcShiftEpEast[9]; // shift EP
    TH2F *h_mZdcShiftEpWest[9];
    TH2F *h_mZdcShiftEpDiff[9]; // Qwest-QEast
    TH2F *h_mZdcShiftEpFull[9]; // Qwest-QEast & shift
    //--------------ZDC EP---------------

    //--------------TPC EP---------------
    // x-axis: runIndex | y-axis: EP
    TH2F *h_mTpcRawEpEast[9]; // raw EP
    TH2F *h_mTpcRawEpWest[9];
    TH2F *h_mTpcRawEpFull[9];

    TH2F *h_mTpcReCenterEpEast[9]; // recenter EP
    TH2F *h_mTpcReCenterEpWest[9];
    TH2F *h_mTpcReCenterEpFull[9];

    TH2F *h_mTpcShiftEpEast[9]; // shift EP
    TH2F *h_mTpcShiftEpWest[9];
    TH2F *h_mTpcShiftEpRanA[9];
    TH2F *h_mTpcShiftEpRanB[9];
    TH2F *h_mTpcShiftEpFull[9];
    //--------------TPC EP---------------

  ClassDef(StEventPlaneHistoManager,1)
};
#endif
