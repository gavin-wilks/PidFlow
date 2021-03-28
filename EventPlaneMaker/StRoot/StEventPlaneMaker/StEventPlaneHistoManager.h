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
    //--------------ZDC EP---------------
    
    //--------------TPC EP---------------
    void initTpcRawEP(); // raw TPC EP
    void fillTpcRawSubEP(TVector2 QEast, TVector2 QWest, int Cent9, int runIndex);
    void fillTpcRawFullEP(TVector2 QFull, int Cent9, int runIndex);
    void writeTpcRawEP();
    //--------------TPC EP---------------

  private:
    //--------------ZDC EP---------------
    TH2F *h_mZdcGainCorr[2][2][8]; // 0: east/west | 1: vertical(x)/horizontal(y) | 2: 7 slats(x)/8 slats(y); | x-axis: runIndex | y-axis: ADC

    // x-axis: runIndex | y-axis: EP
    TH2F *h_mZdcRawEpEast[9]; // raw EP
    TH2F *h_mZdcRawEpWest[9];
    TH2F *h_mZdcRawEpFull[9]; // Qwest-QEast
    //--------------ZDC EP---------------

    //--------------TPC EP---------------
    // x-axis: runIndex | y-axis: EP
    TH2F *h_mTpcRawEpEast[9]; // raw EP
    TH2F *h_mTpcRawEpWest[9];
    TH2F *h_mTpcRawEpFull[9];
    //--------------TPC EP---------------

  ClassDef(StEventPlaneHistoManager,1)
};
#endif
