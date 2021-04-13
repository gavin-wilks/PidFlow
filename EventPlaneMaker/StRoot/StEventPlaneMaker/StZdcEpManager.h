#ifndef StZdcEpManager_h
#define StZdcEpManager_h

#include "TObject.h"
#include "TVector2.h"
#include "TString.h"

class TProfile2D;
class TProfile;
class TFile;

class StZdcEpManager : public TObject
{
  public:
    StZdcEpManager(int energy);
    virtual ~StZdcEpManager();
    void clearZdcEp();
    void initZdcEp(int Cent9, int RunIndex, int VzSign);

    void setZdcSmd(int eastwest,int verthori,int strip,const float zdcsmd);
    float getZdcSmd(int eastwest,int verthori,int strip);

    void readGainCorr();
    void setZdcSmdGainCorr(int eastwest,int verthori,int strip,const float zdcsmd);
    float getZdcSmdGainCorr(int eastwest,int verthori,int strip);
    float getPosition(int eastwest,int verthori,int strip, int mode);
    TVector2 getQEast(int mode);
    TVector2 getQWest(int mode);

    void readReCenterCorr();
    void setZdcSmdCenter();

    void readShiftCorr();
    TVector2 ApplyZdcSmdShiftCorrEast(TVector2 qVector);
    TVector2 ApplyZdcSmdShiftCorrWest(TVector2 qVector);
    float AngleShift(float Psi_shifted);

    void readShiftCorrFull();
    TVector2 ApplyZdcSmdShiftCorrFull(TVector2 qVector);
    TVector2 getQFull(TVector2 QEast, TVector2 QWest);

    void readResolution();
    double getRes1Full(int Cent9);
    double getRes1FullErr(int Cent9);
    double getRes2Full(int Cent9);
    double getRes2FullErr(int Cent9);

  private:

    int mEnergy;
    int mCent9;
    int mRunIndex;
    int mVzSign;
    float mZdcSmd[2][2][8]; // 0: east/west | 1: vertical(x)/horizontal(y) | 2: 7 slats(x)/8 slats(y);
    float mGainCorrFactor[2][2][8];

    // ReCenter Correction | x axis is RunIndex, y axis is Centrality
    TProfile2D *p_mZdcQ1EastVertical[2]; // 0 = vertex pos/neg
    TProfile2D *p_mZdcQ1EastHorizontal[2];
    TProfile2D *p_mZdcQ1WestVertical[2];
    TProfile2D *p_mZdcQ1WestHorizontal[2];
    float mCenterEastVertical, mCenterEastHorizontal, mCenterWestVertical, mCenterWestHorizontal;

    // Shift Correction for East/West/Full
    TProfile2D *p_mZdcQ1EastCos[2][20]; // 0 = vertex pos/neg | 1 = shift correction harmonics
    TProfile2D *p_mZdcQ1EastSin[2][20];
    TProfile2D *p_mZdcQ1WestCos[2][20];
    TProfile2D *p_mZdcQ1WestSin[2][20];
    TProfile2D *p_mZdcQ1FullCos[2][20];
    TProfile2D *p_mZdcQ1FullSin[2][20];

    // EP resolution
    double mZdcFullRes1Val[9];
    double mZdcFullRes1Err[9];
    double mZdcFullRes2Val[9];
    double mZdcFullRes2Err[9];

    TFile *mFile_GainCorrPar;
    TFile *mFile_ReCenterPar;
    TFile *mFile_ShiftPar;
    TFile *mFile_ShiftParFull;
    TFile *mFile_Resolution;

  ClassDef(StZdcEpManager,1)
};

#endif
