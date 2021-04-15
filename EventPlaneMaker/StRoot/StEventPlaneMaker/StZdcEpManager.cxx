// #include "StPicoDstMaker/StPicoDstMaker.h"
// #include "StPicoEvent/StPicoDst.h"
// #include "StPicoEvent/StPicoEvent.h"
// #include "StPicoEvent/StPicoTrack.h"
// #include "StPicoEvent/StPicoBTofPidTraits.h"
#include "StMessMgr.h"

#include "StRoot/StEventPlaneMaker/StZdcEpManager.h"
#include "StRoot/StEventPlaneMaker/StEventPlaneCons.h"

#include "TMath.h"
#include "TString.h"
#include "TFile.h"
#include "TProfile2D.h"
#include "TProfile.h"
#include "TRandom3.h"
#include "TF1.h"

double Resolution_ZdcFull(double *x_val, double *par)
{
  double y;
  double chi = x_val[0];
  double arg = chi*chi/4.0;
  double norm = TMath::Sqrt(TMath::Pi()/2.0)/2.0;

  y = norm*chi*TMath::Exp(-1.0*arg)*(TMath::BesselI0(arg)+TMath::BesselI1(arg));

  return y;
}

ClassImp(StZdcEpManager)

//---------------------------------------------------------------------------------

StZdcEpManager::StZdcEpManager(int energy)
{
  mEnergy = energy;
  clearZdcEp();
  for(int i_eastwest = 0; i_eastwest < 2; ++i_eastwest)
  {
    for(int i_verthori = 0; i_verthori < 2; ++i_verthori)
    {
      for(int i_slat = 0; i_slat < 8; ++i_slat)
      {
	mGainCorrFactor[i_eastwest][i_verthori][i_slat] = -999.9;
      }
    }
  }
  for(int i_cent = 0; i_cent < 9; ++i_cent) 
  {
    mZdcFullRes1Val[i_cent] = 0.0;
    mZdcFullRes1Err[i_cent] = 0.0;
    mZdcFullRes2Val[i_cent] = 0.0;
    mZdcFullRes2Err[i_cent] = 0.0;
  }
}

StZdcEpManager::~StZdcEpManager()
{
  /* */
}

void StZdcEpManager::clearZdcEp()
{
  mCent9 = -1;
  mRunIndex = -1;
  mVzSign = -1;
  for(int i_eastwest = 0; i_eastwest < 2; ++i_eastwest)
  {
    for(int i_verthori = 0; i_verthori < 2; ++i_verthori)
    {
      for(int i_slat = 0; i_slat < 8; ++i_slat)
      {
	mZdcSmd[i_eastwest][i_verthori][i_slat] = 0.0;
      }
    }
  }
  mCenterEastVertical   = -999.9;
  mCenterEastHorizontal = -999.9;
  mCenterWestVertical   = -999.9;
  mCenterWestHorizontal = -999.9;
}

void StZdcEpManager::initZdcEp(int Cent9, int RunIndex, int VzSign)
{
  mCent9 = Cent9;
  mRunIndex = RunIndex;
  mVzSign = VzSign;
}

//---------------------------------------------------------------------------------

void StZdcEpManager::setZdcSmd(int eastwest, int verthori, int slat, const float zdcsmd) 
{
  mZdcSmd[eastwest][verthori][slat] = (zdcsmd > 0.) ? zdcsmd : 0.;
}

float StZdcEpManager::getZdcSmd(int eastwest, int verthori, int slat)
{
  return mZdcSmd[eastwest][verthori][slat];
}

//---------------------------------------------------------------------------------

void StZdcEpManager::readGainCorr()
{
  string InPutFile = Form("StRoot/StEventPlaneUtility/GainCorrPar/file_%s_GainCorrFac.root",recoEP::mBeamEnergy[mEnergy].c_str());
  mFile_GainCorrPar = TFile::Open(InPutFile.c_str());
  for(int i_eastwest = 0; i_eastwest < 2; ++i_eastwest)
  {
    for(int i_verthori = 0; i_verthori < 2; ++i_verthori)
    {
      for(int i_slat = 0; i_slat < 8; ++i_slat)
      {
	string HistName = Form("h_mZdcGainCorrFactor%s%s_%d",recoEP::mEastWest[i_eastwest].c_str(),recoEP::mVertHori[i_verthori].c_str(),i_slat);
	TH1F *h_ZdcGainCorrFac = (TH1F*)mFile_GainCorrPar->Get(HistName.c_str());
	mGainCorrFactor[i_eastwest][i_verthori][i_slat] = h_ZdcGainCorrFac->GetBinContent(1);
	// cout << "i_eastwest = " << i_eastwest << ", i_verthori = " << i_verthori << ", i_slat = " << i_slat << ", mGainCorrFactor = " << mGainCorrFactor[i_eastwest][i_verthori][i_slat] << endl;
      }
    }
  }
}

void StZdcEpManager::setZdcSmdGainCorr(int eastwest, int verthori, int slat, const float zdcsmd)
{
  mZdcSmd[eastwest][verthori][slat] = (zdcsmd > 0.) ? zdcsmd/mGainCorrFactor[eastwest][verthori][slat] : 0.;
  // cout << "input zdc = " << zdcsmd << ", mGainCorrFactor = " << mGainCorrFactor[eastwest][verthori][slat] << ", GainCorred = " << mZdcSmd[eastwest][verthori][slat] << endl;
}

float StZdcEpManager::getZdcSmdGainCorr(int eastwest, int verthori, int slat)
{
  return mZdcSmd[eastwest][verthori][slat];
}

float StZdcEpManager::getPosition(int eastwest, int verthori, int slat, int mode)
{
  //get position of each slat

  float zdcsmd_vert[7] = {0.5,2,3.5,5,6.5,8,9.5};
  float zdcsmd_hori[8] = {1.25,3.25,5.25,7.25,9.25,11.25,13.25,15.25};

  if(mode > 1) // with beam center corrected
  {
    setZdcSmdCenter();
    if(mCenterEastVertical < -100.0 || mCenterEastHorizontal < -100.0 || mCenterWestVertical < -100.0 || mCenterWestHorizontal < -100.0) 
    {
      cout << "Forgot Re-Center!!!!" << endl;
      return 0;
    }
    if(eastwest == 0 && verthori == 0) return zdcsmd_vert[slat]-mCenterEastVertical;
    if(eastwest == 1 && verthori == 0) return mCenterWestVertical-zdcsmd_vert[slat];
    if(eastwest == 0 && verthori == 1) return zdcsmd_hori[slat]/sqrt(2.)-mCenterEastHorizontal;
    if(eastwest == 1 && verthori == 1) return zdcsmd_hori[slat]/sqrt(2.)-mCenterWestHorizontal;
  }
  else // raw beam center returned
  {
    if(eastwest == 0 && verthori == 0) return zdcsmd_vert[slat];
    if(eastwest == 1 && verthori == 0) return -zdcsmd_vert[slat];
    if(eastwest == 0 && verthori == 1) return zdcsmd_hori[slat]/sqrt(2.);
    if(eastwest == 1 && verthori == 1) return zdcsmd_hori[slat]/sqrt(2.);
  }

  return 0;
}

TVector2 StZdcEpManager::getQEast(int mode) 
{
  TVector2 qVector(0.0,0.0);
  float qXsum = 0.; float qYsum = 0.;
  float qXwgt = 0.; float qYwgt = 0.;

  for(int i_vert = 0; i_vert < 7; ++i_vert) // vertical
  {
    qXsum += getPosition(0,0,i_vert,mode)*getZdcSmdGainCorr(0,0,i_vert);
    qXwgt += getZdcSmdGainCorr(0,0,i_vert);
  }
  for(int i_hori = 0; i_hori < 8; ++i_hori) // horizontal
  {
    qYsum += getPosition(0,1,i_hori,mode)*getZdcSmdGainCorr(0,1,i_hori);
    qYwgt += getZdcSmdGainCorr(0,1,i_hori);
  }

  if(qXwgt > 0.0 && qYwgt > 0.0) qVector.Set(qXsum/qXwgt,qYsum/qYwgt);
  if(mode > 2)  qVector = ApplyZdcSmdShiftCorrEast(qVector);

  return qVector;
}

TVector2 StZdcEpManager::getQWest(int mode)
{
  TVector2 qVector(0.0,0.0);
  float qXsum = 0.; float qYsum = 0.;
  float qXwgt = 0.; float qYwgt = 0.;

  for(int i_vert = 0; i_vert < 7; ++i_vert) // vertical
  {
    qXsum += getPosition(1,0,i_vert,mode)*getZdcSmdGainCorr(1,0,i_vert);
    qXwgt += getZdcSmdGainCorr(1,0,i_vert);
  }
  for(int i_hori = 0; i_hori < 8; ++i_hori) // horizontal
  {
    qYsum += getPosition(1,1,i_hori,mode)*getZdcSmdGainCorr(1,1,i_hori);
    qYwgt += getZdcSmdGainCorr(1,1,i_hori);
  }

  if(qXwgt > 0.0 && qYwgt > 0.0) qVector.Set(qXsum/qXwgt,qYsum/qYwgt);
  if(mode > 2) qVector = ApplyZdcSmdShiftCorrWest(qVector);

  return qVector;
}

//---------------------------------------------------------------------------------

void StZdcEpManager::readReCenterCorr()
{
  string InPutFile = Form("StRoot/StEventPlaneUtility/ReCenterParameter/file_%s_ReCenterParameter.root",recoEP::mBeamEnergy[mEnergy].c_str());
  mFile_ReCenterPar = TFile::Open(InPutFile.c_str());

  string ProName;
  for(int i_vz = 0; i_vz < 2; ++i_vz)
  {
    ProName = Form("p_mZdcQ1EastVertical_%s",recoEP::mVStr[i_vz].c_str());
    p_mZdcQ1EastVertical[i_vz] = (TProfile2D*)mFile_ReCenterPar->Get(ProName.c_str());
    ProName = Form("p_mZdcQ1EastHorizontal_%s",recoEP::mVStr[i_vz].c_str());
    p_mZdcQ1EastHorizontal[i_vz] = (TProfile2D*)mFile_ReCenterPar->Get(ProName.c_str());

    ProName = Form("p_mZdcQ1WestVertical_%s",recoEP::mVStr[i_vz].c_str());
    p_mZdcQ1WestVertical[i_vz] = (TProfile2D*)mFile_ReCenterPar->Get(ProName.c_str());
    ProName = Form("p_mZdcQ1WestHorizontal_%s",recoEP::mVStr[i_vz].c_str());
    p_mZdcQ1WestHorizontal[i_vz] = (TProfile2D*)mFile_ReCenterPar->Get(ProName.c_str());
  }
}

void StZdcEpManager::setZdcSmdCenter()
{
  int binEastVertical = p_mZdcQ1EastVertical[mVzSign]->FindBin((double)mRunIndex,(double)mCent9);
  mCenterEastVertical = p_mZdcQ1EastVertical[mVzSign]->GetBinContent(binEastVertical);

  int binEastHorizontal = p_mZdcQ1EastHorizontal[mVzSign]->FindBin((double)mRunIndex,(double)mCent9);
  mCenterEastHorizontal = p_mZdcQ1EastHorizontal[mVzSign]->GetBinContent(binEastHorizontal);

  int binWestVertical = p_mZdcQ1WestVertical[mVzSign]->FindBin((double)mRunIndex,(double)mCent9);
  mCenterWestVertical = -1.0*p_mZdcQ1WestVertical[mVzSign]->GetBinContent(binWestVertical);

  int binWestHorizontal = p_mZdcQ1WestHorizontal[mVzSign]->FindBin((double)mRunIndex,(double)mCent9);
  mCenterWestHorizontal = p_mZdcQ1WestHorizontal[mVzSign]->GetBinContent(binWestHorizontal);
  // cout << "mRunIndex = " << mRunIndex << ", mCent9 = " << mCent9 << ", mVzSign = " << mVzSign << ", mCenterEastVertical = " << mCenterEastVertical << ", mCenterEastHorizontal = " << mCenterEastHorizontal << ", mCenterWestVertical = " << mCenterWestVertical << ", mCenterWestHorizontal = " << mCenterWestHorizontal << endl;
}

//---------------------------------------------------------------------------------

void StZdcEpManager::readShiftCorr()
{
  string InPutFile = Form("StRoot/StEventPlaneUtility/ShiftParameter/file_%s_ShiftParameter.root",recoEP::mBeamEnergy[mEnergy].c_str());
  mFile_ShiftPar = TFile::Open(InPutFile.c_str());

  string ProName;
  for(int i_vz = 0; i_vz < 2; ++i_vz)
  {
    for(int i_shift = 0; i_shift < recoEP::mNumShiftOrder; ++i_shift) // Shift Order
    {
      ProName = Form("p_mZdcQ1EastCos_%s_%d",recoEP::mVStr[i_vz].c_str(),i_shift);
      p_mZdcQ1EastCos[i_vz][i_shift] = (TProfile2D*)mFile_ShiftPar->Get(ProName.c_str());
      ProName = Form("p_mZdcQ1EastSin_%s_%d",recoEP::mVStr[i_vz].c_str(),i_shift);
      p_mZdcQ1EastSin[i_vz][i_shift] = (TProfile2D*)mFile_ShiftPar->Get(ProName.c_str());

      ProName = Form("p_mZdcQ1WestCos_%s_%d",recoEP::mVStr[i_vz].c_str(),i_shift);
      p_mZdcQ1WestCos[i_vz][i_shift] = (TProfile2D*)mFile_ShiftPar->Get(ProName.c_str());
      ProName = Form("p_mZdcQ1WestSin_%s_%d",recoEP::mVStr[i_vz].c_str(),i_shift);
      p_mZdcQ1WestSin[i_vz][i_shift] = (TProfile2D*)mFile_ShiftPar->Get(ProName.c_str());
    }
  }
}

TVector2 StZdcEpManager::ApplyZdcSmdShiftCorrEast(TVector2 qVector)
{
  TVector2 qVecShift(0.0,0.0);
  float Psi_ReCenter = TMath::ATan2(qVector.Y(),qVector.X());
  float delta_Psi = 0.0;
  float Psi_Shift;

  for(Int_t i_shift = 0; i_shift < recoEP::mNumShiftOrder; ++i_shift) // Shift Order loop
  {
    int bin_Cos = p_mZdcQ1EastCos[mVzSign][i_shift]->FindBin((double)mRunIndex,(double)mCent9);
    float mean_Cos = p_mZdcQ1EastCos[mVzSign][i_shift]->GetBinContent(bin_Cos);

    int bin_Sin = p_mZdcQ1EastSin[mVzSign][i_shift]->FindBin((double)mRunIndex,(double)mCent9);
    float mean_Sin = p_mZdcQ1EastSin[mVzSign][i_shift]->GetBinContent(bin_Sin);

    delta_Psi += (2.0/((float)i_shift+1.0))*(-1.0*mean_Sin*TMath::Cos(((float)i_shift+1.0)*Psi_ReCenter)+mean_Cos*TMath::Sin(((float)i_shift+1.0)*Psi_ReCenter));
  }

  float Psi_Shift_raw = Psi_ReCenter + delta_Psi;
  Psi_Shift = AngleShift(Psi_Shift_raw);

  qVecShift.Set(TMath::Cos(Psi_Shift),TMath::Sin(Psi_Shift));

  return qVecShift;
}

TVector2 StZdcEpManager::ApplyZdcSmdShiftCorrWest(TVector2 qVector)
{
  TVector2 qVecShift(0.0,0.0);
  float Psi_ReCenter = TMath::ATan2(qVector.Y(),qVector.X());
  float delta_Psi = 0.0;
  float Psi_Shift;

  for(Int_t i_shift = 0; i_shift < recoEP::mNumShiftOrder; ++i_shift) // Shift Order loop
  {
    int bin_Cos = p_mZdcQ1WestCos[mVzSign][i_shift]->FindBin((double)mRunIndex,(double)mCent9);
    float mean_Cos = p_mZdcQ1WestCos[mVzSign][i_shift]->GetBinContent(bin_Cos);

    int bin_Sin = p_mZdcQ1WestSin[mVzSign][i_shift]->FindBin((double)mRunIndex,(double)mCent9);
    float mean_Sin = p_mZdcQ1WestSin[mVzSign][i_shift]->GetBinContent(bin_Sin);

    delta_Psi += (2.0/((float)i_shift+1.0))*(-1.0*mean_Sin*TMath::Cos(((float)i_shift+1.0)*Psi_ReCenter)+mean_Cos*TMath::Sin(((float)i_shift+1.0)*Psi_ReCenter));
  }

  float Psi_Shift_raw = Psi_ReCenter + delta_Psi;
  Psi_Shift = AngleShift(Psi_Shift_raw);

  qVecShift.Set(TMath::Cos(Psi_Shift),TMath::Sin(Psi_Shift));

  return qVecShift;
}

float StZdcEpManager::AngleShift(float Psi_raw)
{
  float Psi_Corr = Psi_raw;
  if(Psi_raw > 1.0*TMath::Pi())
  {
    Psi_Corr = Psi_raw - TMath::TwoPi();
  }
  if(Psi_raw < -1.0*TMath::Pi())
  {
    Psi_Corr = Psi_raw + TMath::TwoPi();
  }

  return Psi_Corr;
}

//---------------------------------------------------------------------------------


void StZdcEpManager::readShiftCorrFull()
{
  string InPutFile = Form("StRoot/StEventPlaneUtility/ShiftParameterFull/file_%s_ShiftParameterFull.root",recoEP::mBeamEnergy[mEnergy].c_str());
  mFile_ShiftPar = TFile::Open(InPutFile.c_str());

  string ProName;
  for(int i_vz = 0; i_vz < 2; ++i_vz)
  {
    for(int i_shift = 0; i_shift < recoEP::mNumShiftOrder; ++i_shift) // Shift Order
    {
      ProName = Form("p_mZdcQ1FullCos_%s_%d",recoEP::mVStr[i_vz].c_str(),i_shift);
      p_mZdcQ1FullCos[i_vz][i_shift] = (TProfile2D*)mFile_ShiftPar->Get(ProName.c_str());
      ProName = Form("p_mZdcQ1FullSin_%s_%d",recoEP::mVStr[i_vz].c_str(),i_shift);
      p_mZdcQ1FullSin[i_vz][i_shift] = (TProfile2D*)mFile_ShiftPar->Get(ProName.c_str());
    }
  }
}

TVector2 StZdcEpManager::ApplyZdcSmdShiftCorrFull(TVector2 qVector)
{
  TVector2 qVecShift(0.0,0.0);
  float Psi_ReCenter = TMath::ATan2(qVector.Y(),qVector.X());
  float delta_Psi = 0.0;
  float Psi_Shift;

  for(Int_t i_shift = 0; i_shift < recoEP::mNumShiftOrder; ++i_shift) // Shift Order loop
  {
    int bin_Cos = p_mZdcQ1FullCos[mVzSign][i_shift]->FindBin((double)mRunIndex,(double)mCent9);
    float mean_Cos = p_mZdcQ1FullCos[mVzSign][i_shift]->GetBinContent(bin_Cos);

    int bin_Sin = p_mZdcQ1FullSin[mVzSign][i_shift]->FindBin((double)mRunIndex,(double)mCent9);
    float mean_Sin = p_mZdcQ1FullSin[mVzSign][i_shift]->GetBinContent(bin_Sin);

    delta_Psi += (2.0/((float)i_shift+1.0))*(-1.0*mean_Sin*TMath::Cos(((float)i_shift+1.0)*Psi_ReCenter)+mean_Cos*TMath::Sin(((float)i_shift+1.0)*Psi_ReCenter));
  }

  float Psi_Shift_raw = Psi_ReCenter + delta_Psi;
  Psi_Shift = AngleShift(Psi_Shift_raw);

  qVecShift.Set(TMath::Cos(Psi_Shift),TMath::Sin(Psi_Shift));

  return qVecShift;
}

TVector2 StZdcEpManager::getQFull(TVector2 QEast, TVector2 QWest)
{
  TVector2 qVector = QWest-QEast;
  TVector2 qVecShift = ApplyZdcSmdShiftCorrFull(qVector);

  return qVecShift;
}

//---------------------------------------------------------------------------------

void StZdcEpManager::readResolution()
{
  string InPutFile = Form("StRoot/StEventPlaneUtility/Resolution/file_%s_Resolution.root",recoEP::mBeamEnergy[mEnergy].c_str());
  mFile_Resolution = TFile::Open(InPutFile.c_str());

  TF1 *f_res = new TF1("f_res",Resolution_ZdcFull,0,10,0);

  // calculate 1st full event plane resolution
  TProfile *p_mZdcSubRes1 = (TProfile*)mFile_Resolution->Get("p_mZdcSubRes1");
  for(int i_cent = 0; i_cent < 9; ++i_cent)
  {
    const double resRaw = p_mZdcSubRes1->GetBinContent(p_mZdcSubRes1->FindBin(i_cent));
    const double errRaw = p_mZdcSubRes1->GetBinError(p_mZdcSubRes1->FindBin(i_cent));
    if(resRaw > 0)
    {
      const double resSub = TMath::Sqrt(resRaw);
      const double errSub = errRaw/(2.0*TMath::Sqrt(resRaw));

      const double chiSub = f_res->GetX(resSub);
      const double errChiSub = errSub/f_res->Derivative(chiSub);
      const double chiFull = chiSub*TMath::Sqrt(2.0);
      mZdcFullRes1Val[i_cent] = f_res->Eval(chiFull);
      mZdcFullRes1Err[i_cent] = f_res->Derivative(chiFull)*errChiSub*TMath::Sqrt(2.0);
    }
    // cout << "i_cent = " << i_cent << ", resRaw = " << resRaw << ", resFull = " << mZdcFullRes1Val[i_cent] << " +/- " << mZdcFullRes1Err[i_cent] << endl;
  }

  // calculate 2nd full event plane resolution
  TProfile *p_mZdcSubRes2 = (TProfile*)mFile_Resolution->Get("p_mZdcSubRes2");
  for(int i_cent = 0; i_cent < 9; ++i_cent)
  {
    const double resRaw = p_mZdcSubRes2->GetBinContent(p_mZdcSubRes2->FindBin(i_cent));
    const double errRaw = p_mZdcSubRes2->GetBinError(p_mZdcSubRes2->FindBin(i_cent));
    if(resRaw > 0)
    {
      const double resSub = TMath::Sqrt(resRaw);
      const double errSub = errRaw/(2.0*TMath::Sqrt(resRaw));

      const double chiSub = f_res->GetX(resSub);
      const double errChiSub = errSub/f_res->Derivative(chiSub);
      const double chiFull = chiSub*TMath::Sqrt(2.0);
      mZdcFullRes2Val[i_cent] = f_res->Eval(chiFull);
      mZdcFullRes2Err[i_cent] = f_res->Derivative(chiFull)*errChiSub*TMath::Sqrt(2.0);
    }
    // cout << "i_cent = " << i_cent << ", resRaw = " << resRaw << ", resFull = " << mZdcFullRes2Val[i_cent] << " +/- " << mZdcFullRes2Err[i_cent] << endl;
  }
}

double StZdcEpManager::getRes1Full(int Cent9)
{
  return mZdcFullRes1Val[Cent9];
}

double StZdcEpManager::getRes1FullErr(int Cent9)
{
  return mZdcFullRes1Err[Cent9];
}

double StZdcEpManager::getRes2Full(int Cent9)
{
  return mZdcFullRes2Val[Cent9];
}

double StZdcEpManager::getRes2FullErr(int Cent9)
{
  return mZdcFullRes2Err[Cent9];
}
//---------------------------------------------------------------------------------
