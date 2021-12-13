#include "StMessMgr.h"
#include "TMath.h"
#include "TString.h"
#include "TFile.h"
#include "TProfile2D.h"
#include "TProfile.h"
#include "TNtuple.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TVector3.h"

#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"

#include "StRoot/StEventPlaneMaker/StEpdEpManager.h"
#include "StRoot/StEventPlaneMaker/StEventPlaneCons.h"

double Resolution_EpdFull(double *x_val, double *par)
{
  double y;
  double chi = x_val[0];
  double arg = chi*chi/4.0;
  double norm = TMath::Sqrt(TMath::Pi()/2.0)/2.0;

  y = norm*chi*TMath::Exp(-1.0*arg)*(TMath::BesselI0(arg)+TMath::BesselI1(arg));

  return y;
}

ClassImp(StEpdEpManager)

//---------------------------------------------------------------------------------
StEpdEpManager::StEpdEpManager(Int_t energy)
{
  mEnergy = energy;
  clearEpdEp();
}

StEpdEpManager::~StEpdEpManager()
{
  /* */
}

void StEpdEpManager::clearEpdEp()
{
  mCent9 = -1;
  mRunIndex = -1;
}

void StEpdEpManager::initEpdEp(int Cent9, int RunIndex)
{
  mCent9 = Cent9;
  mRunIndex = RunIndex;

  for(int i_cent = 0; i_cent < 9; ++i_cent)
  {
    for(int order = 1; order <= recoEP::mEpdEpOrder; ++order)
    {  
      mEpdSubResVal[i_cent][order-1]  = 0.0;
      mEpdSubResErr[i_cent][order-1]  = 0.0;
    }
  }
}
//---------------------------------------------------------------------------------
void StEpdEpManager::readResolution()
{
  string InPutFile = Form("StRoot/StEventPlaneUtility/EpdResolution/file_%s_EpdRes.root",recoEP::mBeamEnergy[mEnergy].c_str());
  mInPutFile_Res = TFile::Open(InPutFile.c_str());
  
  for(int order = 1; order <= recoEP::mEpdEpOrder; ++order)
  {
    p_mEpdSubRes[order-1] = (TProfile*)mInPutFile_Res->Get(Form("p_mEpdSubRes%d",order));
    for(int i_cent = 0; i_cent < 9; ++i_cent)
    {
      const double resRaw = p_mEpdSubRes[order-1]->GetBinContent(p_mEpdSubRes[order-1]->FindBin(i_cent));
      const double errRaw = p_mEpdSubRes[order-1]->GetBinError(p_mEpdSubRes[order-1]->FindBin(i_cent));
      if(resRaw > 0)
      {
        mEpdSubResVal[i_cent][order-1] = TMath::Sqrt(resRaw);
        mEpdSubResErr[i_cent][order-1] = errRaw/(2.0*TMath::Sqrt(resRaw));
      }
    }
  }
}

double StEpdEpManager::getResSub(int Cent9, int order)
{
  return mEpdSubResVal[Cent9][order-1];
}

double StEpdEpManager::getResSubErr(int Cent9, int order)
{
  return mEpdSubResErr[Cent9][order-1];
}
