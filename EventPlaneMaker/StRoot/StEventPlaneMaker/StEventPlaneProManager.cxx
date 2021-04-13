#include <string>

#include <TProfile.h>
#include <TProfile2D.h>
#include <TMath.h>
#include <TString.h>

#include "StRoot/StEventPlaneMaker/StEventPlaneProManager.h"
#include "StRoot/StEventPlaneMaker/StEventPlaneCons.h"

ClassImp(StEventPlaneProManager)

//---------------------------------------------------------------------------------

StEventPlaneProManager::StEventPlaneProManager()
{
}

StEventPlaneProManager::~StEventPlaneProManager()
{
  /* */
}

//---------------------------------------------------------------------------------

//---------------------------------------------------------------------------------
// ZDC-SMD ReCenter Correction
void StEventPlaneProManager::initZdcReCenter()
{
  string ProName;

  for(int i_vz = 0; i_vz < 2; ++i_vz)
  {
    ProName = Form("p_mZdcQ1EastVertical_%s",recoEP::mVStr[i_vz].c_str());
    p_mZdcQ1EastVertical[i_vz] = new TProfile2D(ProName.c_str(),ProName.c_str(),recoEP::mNumOfRunIndex,-0.5,(float)recoEP::mNumOfRunIndex-0.5,9,-0.5,8.5);
    ProName = Form("p_mZdcQ1EastHorizontal_%s",recoEP::mVStr[i_vz].c_str());
    p_mZdcQ1EastHorizontal[i_vz] = new TProfile2D(ProName.c_str(),ProName.c_str(),recoEP::mNumOfRunIndex,-0.5,(float)recoEP::mNumOfRunIndex-0.5,9,-0.5,8.5);

    ProName = Form("p_mZdcQ1WestVertical_%s",recoEP::mVStr[i_vz].c_str());
    p_mZdcQ1WestVertical[i_vz] = new TProfile2D(ProName.c_str(),ProName.c_str(),recoEP::mNumOfRunIndex,-0.5,(float)recoEP::mNumOfRunIndex-0.5,9,-0.5,8.5);
    ProName = Form("p_mZdcQ1WestHorizontal_%s",recoEP::mVStr[i_vz].c_str());
    p_mZdcQ1WestHorizontal[i_vz] = new TProfile2D(ProName.c_str(),ProName.c_str(),recoEP::mNumOfRunIndex,-0.5,(float)recoEP::mNumOfRunIndex-0.5,9,-0.5,8.5);
  }
}

void StEventPlaneProManager::fillZdcReCenterEast(TVector2 qVector, int Cent9, int RunIndex, int VzSign)
{
  const double qx = qVector.X();
  const double qy = qVector.Y();

  // Event Plane method
  p_mZdcQ1EastVertical[VzSign]->Fill((double)RunIndex,(double)Cent9,(double)qx);
  p_mZdcQ1EastHorizontal[VzSign]->Fill((double)RunIndex,(double)Cent9,(double)qy);
}

void StEventPlaneProManager::fillZdcReCenterWest(TVector2 qVector, int Cent9, int RunIndex, int VzSign)
{
  const double qx = qVector.X();
  const double qy = qVector.Y();

  p_mZdcQ1WestVertical[VzSign]->Fill((double)RunIndex,(double)Cent9,(double)qx);
  p_mZdcQ1WestHorizontal[VzSign]->Fill((double)RunIndex,(double)Cent9,(double)qy);
}

void StEventPlaneProManager::writeZdcReCenter()
{
  for(int i_vz = 0; i_vz < 2; ++i_vz)
  {
    p_mZdcQ1EastVertical[i_vz]->Write();
    p_mZdcQ1EastHorizontal[i_vz]->Write();
    p_mZdcQ1WestVertical[i_vz]->Write();
    p_mZdcQ1WestHorizontal[i_vz]->Write();
  }
}

// ZDC-SMD Shift Correction
void StEventPlaneProManager::initZdcShift()
{
  string ProName;
  for(int i_vz = 0; i_vz < 2; ++i_vz)
  {
    for(int i_shift = 0; i_shift < recoEP::mNumShiftOrder; ++i_shift) // Shift Order
    {
      ProName = Form("p_mZdcQ1EastCos_%s_%d",recoEP::mVStr[i_vz].c_str(),i_shift);
      p_mZdcQ1EastCos[i_vz][i_shift] = new TProfile2D(ProName.c_str(),ProName.c_str(),recoEP::mNumOfRunIndex,-0.5,(float)recoEP::mNumOfRunIndex-0.5,9,-0.5,8.5);
      ProName = Form("p_mZdcQ1EastSin_%s_%d",recoEP::mVStr[i_vz].c_str(),i_shift);
      p_mZdcQ1EastSin[i_vz][i_shift] = new TProfile2D(ProName.c_str(),ProName.c_str(),recoEP::mNumOfRunIndex,-0.5,(float)recoEP::mNumOfRunIndex-0.5,9,-0.5,8.5);

      ProName = Form("p_mZdcQ1WestCos_%s_%d",recoEP::mVStr[i_vz].c_str(),i_shift);
      p_mZdcQ1WestCos[i_vz][i_shift] = new TProfile2D(ProName.c_str(),ProName.c_str(),recoEP::mNumOfRunIndex,-0.5,(float)recoEP::mNumOfRunIndex-0.5,9,-0.5,8.5);
      ProName = Form("p_mZdcQ1WestSin_%s_%d",recoEP::mVStr[i_vz].c_str(),i_shift);
      p_mZdcQ1WestSin[i_vz][i_shift] = new TProfile2D(ProName.c_str(),ProName.c_str(),recoEP::mNumOfRunIndex,-0.5,(float)recoEP::mNumOfRunIndex-0.5,9,-0.5,8.5);
    }
  }
}

void StEventPlaneProManager::fillZdcShiftEast(TVector2 qVector, int Cent9, int RunIndex, int VzSign)
{
  double Psi = TMath::ATan2(qVector.Y(),qVector.X());
  for(int i_shift = 0; i_shift < recoEP::mNumShiftOrder; ++i_shift)
  {
    p_mZdcQ1EastCos[VzSign][i_shift]->Fill((double)RunIndex,(double)Cent9,TMath::Cos((i_shift+1)*Psi));
    p_mZdcQ1EastSin[VzSign][i_shift]->Fill((double)RunIndex,(double)Cent9,TMath::Sin((i_shift+1)*Psi));
  }
}

void StEventPlaneProManager::fillZdcShiftWest(TVector2 qVector, int Cent9, int RunIndex, int VzSign)
{
  float Psi = TMath::ATan2(qVector.Y(),qVector.X());
  for(int i_shift = 0; i_shift < recoEP::mNumShiftOrder; ++i_shift)
  {
    p_mZdcQ1WestCos[VzSign][i_shift]->Fill((double)RunIndex,(double)Cent9,TMath::Cos((i_shift+1)*Psi));
    p_mZdcQ1WestSin[VzSign][i_shift]->Fill((double)RunIndex,(double)Cent9,TMath::Sin((i_shift+1)*Psi));
  }
}

void StEventPlaneProManager::writeZdcShift()
{
  for(int i_vz = 0; i_vz < 2; ++i_vz)
  {
    for(int i_shift = 0; i_shift < recoEP::mNumShiftOrder; ++i_shift) // Shift Order
    {
      p_mZdcQ1EastCos[i_vz][i_shift]->Write();
      p_mZdcQ1EastSin[i_vz][i_shift]->Write();
      p_mZdcQ1WestCos[i_vz][i_shift]->Write();
      p_mZdcQ1WestSin[i_vz][i_shift]->Write();
    }
  }
}

// ZDC-SMD Full Shift Correction
void StEventPlaneProManager::initZdcShiftFull()
{
  string ProName;
  for(int i_vz = 0; i_vz < 2; ++i_vz)
  {
    for(int i_shift = 0; i_shift < recoEP::mNumShiftOrder; ++i_shift) // Shift Order
    {
      ProName = Form("p_mZdcQ1FullCos_%s_%d",recoEP::mVStr[i_vz].c_str(),i_shift);
      p_mZdcQ1FullCos[i_vz][i_shift] = new TProfile2D(ProName.c_str(),ProName.c_str(),recoEP::mNumOfRunIndex,-0.5,(float)recoEP::mNumOfRunIndex-0.5,9,-0.5,8.5);
      ProName = Form("p_mZdcQ1FullSin_%s_%d",recoEP::mVStr[i_vz].c_str(),i_shift);
      p_mZdcQ1FullSin[i_vz][i_shift] = new TProfile2D(ProName.c_str(),ProName.c_str(),recoEP::mNumOfRunIndex,-0.5,(float)recoEP::mNumOfRunIndex-0.5,9,-0.5,8.5);
    }
  }
}

void StEventPlaneProManager::fillZdcShiftFull(TVector2 qVector, int Cent9, int RunIndex, int VzSign)
{
  double Psi = TMath::ATan2(qVector.Y(),qVector.X());
  for(int i_shift = 0; i_shift < recoEP::mNumShiftOrder; ++i_shift)
  {
    p_mZdcQ1FullCos[VzSign][i_shift]->Fill((double)RunIndex,(double)Cent9,TMath::Cos((i_shift+1)*Psi));
    p_mZdcQ1FullSin[VzSign][i_shift]->Fill((double)RunIndex,(double)Cent9,TMath::Sin((i_shift+1)*Psi));
  }
}

void StEventPlaneProManager::writeZdcShiftFull()
{
  for(int i_vz = 0; i_vz < 2; ++i_vz)
  {
    for(int i_shift = 0; i_shift < recoEP::mNumShiftOrder; ++i_shift) // Shift Order
    {
      p_mZdcQ1FullCos[i_vz][i_shift]->Write();
      p_mZdcQ1FullSin[i_vz][i_shift]->Write();
    }
  }
}

// ZDC-SMD Sub EP Resolution
void StEventPlaneProManager::initZdcResolution()
{
  p_mZdcSubRes1QA = new TProfile2D("p_mZdcSubRes1QA","p_mZdcSubRes1QA",recoEP::mNumOfRunIndex,-0.5,(float)recoEP::mNumOfRunIndex-0.5,9,-0.5,8.5);
  p_mZdcSubRes2QA = new TProfile2D("p_mZdcSubRes2QA","p_mZdcSubRes2QA",recoEP::mNumOfRunIndex,-0.5,(float)recoEP::mNumOfRunIndex-0.5,9,-0.5,8.5);
  p_mZdcSubRes1 = new TProfile("p_mZdcSubRes1","p_mZdcRes1",9,-0.5,8.5);
  p_mZdcSubRes2 = new TProfile("p_mZdcSubRes2","p_mZdcRes2",9,-0.5,8.5);
}

void StEventPlaneProManager::fillZdcResSub(TVector2 QEast, TVector2 QWest, int Cent9, int RunIndex)
{
  double PsiEast = TMath::ATan2(QEast.Y(),QEast.X());
  double PsiWest = TMath::ATan2(QWest.Y(),QWest.X());
  double res1 = TMath::Cos(PsiWest-PsiEast+TMath::Pi());
  double res2 = TMath::Cos(2.0*(PsiWest-PsiEast+TMath::Pi()));
  p_mZdcSubRes1QA->Fill((double)RunIndex,(double)Cent9,res1);
  p_mZdcSubRes2QA->Fill((double)RunIndex,(double)Cent9,res2);
  p_mZdcSubRes1->Fill((double)Cent9,res1);
  p_mZdcSubRes2->Fill((double)Cent9,res2);
}

void StEventPlaneProManager::writeZdcResolution()
{
  p_mZdcSubRes1QA->Write();
  p_mZdcSubRes2QA->Write();
  p_mZdcSubRes1->Write();
  p_mZdcSubRes2->Write();
}
//---------------------------------------------------------------------------------

//---------------------------------------------------------------------------------
// TPC ReCenter Correction
void StEventPlaneProManager::initTpcReCenter()
{
  string ProName;
  for(int i_vz = 0; i_vz < 2; ++i_vz)
  {
    ProName = Form("p_mTpcq2xEast_%s",recoEP::mVStr[i_vz].c_str());
    p_mTpcq2xEast[i_vz] = new TProfile2D(ProName.c_str(),ProName.c_str(),recoEP::mNumOfRunIndex,-0.5,(float)recoEP::mNumOfRunIndex-0.5,9,-0.5,8.5); // neg eta
    ProName = Form("p_mTpcq2yEast_%s",recoEP::mVStr[i_vz].c_str());
    p_mTpcq2yEast[i_vz] = new TProfile2D(ProName.c_str(),ProName.c_str(),recoEP::mNumOfRunIndex,-0.5,(float)recoEP::mNumOfRunIndex-0.5,9,-0.5,8.5);
    ProName = Form("p_mTpcq2xWest_%s",recoEP::mVStr[i_vz].c_str());
    p_mTpcq2xWest[i_vz] = new TProfile2D(ProName.c_str(),ProName.c_str(),recoEP::mNumOfRunIndex,-0.5,(float)recoEP::mNumOfRunIndex-0.5,9,-0.5,8.5); // pos eta
    ProName = Form("p_mTpcq2yWest_%s",recoEP::mVStr[i_vz].c_str());
    p_mTpcq2yWest[i_vz] = new TProfile2D(ProName.c_str(),ProName.c_str(),recoEP::mNumOfRunIndex,-0.5,(float)recoEP::mNumOfRunIndex-0.5,9,-0.5,8.5);

    ProName = Form("p_mTpcq2xFull_%s",recoEP::mVStr[i_vz].c_str());
    p_mTpcq2xFull[i_vz] = new TProfile2D(ProName.c_str(),ProName.c_str(),recoEP::mNumOfRunIndex,-0.5,(float)recoEP::mNumOfRunIndex-0.5,9,-0.5,8.5);
    ProName = Form("p_mTpcq2yFull_%s",recoEP::mVStr[i_vz].c_str());
    p_mTpcq2yFull[i_vz] = new TProfile2D(ProName.c_str(),ProName.c_str(),recoEP::mNumOfRunIndex,-0.5,(float)recoEP::mNumOfRunIndex-0.5,9,-0.5,8.5);
  }
}

void StEventPlaneProManager::fillTpcReCenterEast(TVector2 q2Vector, int Cent9, int RunIndex, int VzSign, double pt)
{
  const double q2x = q2Vector.X();
  const double q2y = q2Vector.Y();

  double weight;
  if(pt <= recoEP::mPrimPtWeight)
  {
    weight = pt;
  }
  if(pt > recoEP::mPrimPtWeight)
  {
    weight = recoEP::mPrimPtWeight;
  }

  p_mTpcq2xEast[VzSign]->Fill((double)RunIndex,(double)Cent9,(double)q2x,(double)weight);
  p_mTpcq2yEast[VzSign]->Fill((double)RunIndex,(double)Cent9,(double)q2y,(double)weight);
}

void StEventPlaneProManager::fillTpcReCenterWest(TVector2 q2Vector, int Cent9, int RunIndex, int VzSign, double pt)
{
  const double q2x = q2Vector.X();
  const double q2y = q2Vector.Y();

  double weight;
  if(pt <= recoEP::mPrimPtWeight)
  {
    weight = pt;
  }
  if(pt > recoEP::mPrimPtWeight)
  {
    weight = recoEP::mPrimPtWeight;
  }

  p_mTpcq2xWest[VzSign]->Fill((double)RunIndex,(double)Cent9,(double)q2x,(double)weight);
  p_mTpcq2yWest[VzSign]->Fill((double)RunIndex,(double)Cent9,(double)q2y,(double)weight);
}

void StEventPlaneProManager::fillTpcReCenterFull(TVector2 q2Vector, int Cent9, int RunIndex, int VzSign, double pt)
{
  const double q2x = q2Vector.X();
  const double q2y = q2Vector.Y();

  double weight;
  if(pt <= recoEP::mPrimPtWeight)
  {
    weight = pt;
  }
  if(pt > recoEP::mPrimPtWeight)
  {
    weight = recoEP::mPrimPtWeight;
  }

  p_mTpcq2xFull[VzSign]->Fill((double)RunIndex,(double)Cent9,(double)q2x,(double)weight);
  p_mTpcq2yFull[VzSign]->Fill((double)RunIndex,(double)Cent9,(double)q2y,(double)weight);
}

void StEventPlaneProManager::writeTpcReCenter()
{
  for(int i_vz = 0; i_vz < 2; ++i_vz) // vertex pos/neg
  {
    p_mTpcq2xEast[i_vz]->Write();
    p_mTpcq2yEast[i_vz]->Write();
    p_mTpcq2xWest[i_vz]->Write();
    p_mTpcq2yWest[i_vz]->Write();

    p_mTpcq2xFull[i_vz]->Write();
    p_mTpcq2yFull[i_vz]->Write();
  }
}

// TPC Shift Correction
void StEventPlaneProManager::initTpcShift()
{
  string ProName;
  for(int i_vz = 0; i_vz < 2; ++i_vz)
  {
    for(int i_shift = 0; i_shift < recoEP::mNumShiftOrder; ++i_shift) // Shift Order
    {
      ProName = Form("p_mTpcQ2EastCos_%s_%d",recoEP::mVStr[i_vz].c_str(),i_shift);
      p_mTpcQ2EastCos[i_vz][i_shift] = new TProfile2D(ProName.c_str(),ProName.c_str(),recoEP::mNumOfRunIndex,-0.5,(float)recoEP::mNumOfRunIndex-0.5,9,-0.5,8.5);
      ProName = Form("p_mTpcQ2EastSin_%s_%d",recoEP::mVStr[i_vz].c_str(),i_shift);
      p_mTpcQ2EastSin[i_vz][i_shift] = new TProfile2D(ProName.c_str(),ProName.c_str(),recoEP::mNumOfRunIndex,-0.5,(float)recoEP::mNumOfRunIndex-0.5,9,-0.5,8.5);

      ProName = Form("p_mTpcQ2WestCos_%s_%d",recoEP::mVStr[i_vz].c_str(),i_shift);
      p_mTpcQ2WestCos[i_vz][i_shift] = new TProfile2D(ProName.c_str(),ProName.c_str(),recoEP::mNumOfRunIndex,-0.5,(float)recoEP::mNumOfRunIndex-0.5,9,-0.5,8.5);
      ProName = Form("p_mTpcQ2WestSin_%s_%d",recoEP::mVStr[i_vz].c_str(),i_shift);
      p_mTpcQ2WestSin[i_vz][i_shift] = new TProfile2D(ProName.c_str(),ProName.c_str(),recoEP::mNumOfRunIndex,-0.5,(float)recoEP::mNumOfRunIndex-0.5,9,-0.5,8.5);

      ProName = Form("p_mTpcQ2FullCos_%s_%d",recoEP::mVStr[i_vz].c_str(),i_shift);
      p_mTpcQ2FullCos[i_vz][i_shift] = new TProfile2D(ProName.c_str(),ProName.c_str(),recoEP::mNumOfRunIndex,-0.5,(float)recoEP::mNumOfRunIndex-0.5,9,-0.5,8.5);
      ProName = Form("p_mTpcQ2FullSin_%s_%d",recoEP::mVStr[i_vz].c_str(),i_shift);
      p_mTpcQ2FullSin[i_vz][i_shift] = new TProfile2D(ProName.c_str(),ProName.c_str(),recoEP::mNumOfRunIndex,-0.5,(float)recoEP::mNumOfRunIndex-0.5,9,-0.5,8.5);
    }
  }
}

void StEventPlaneProManager::fillTpcShiftEast(TVector2 qVector, int Cent9, int RunIndex, int VzSign)
{
  double Psi2 = TMath::ATan2(qVector.Y(),qVector.X())/2.0;
  for(int i_shift = 0; i_shift < recoEP::mNumShiftOrder; ++i_shift)
  {
    p_mTpcQ2EastCos[VzSign][i_shift]->Fill((double)RunIndex,(double)Cent9,TMath::Cos(2.0*(i_shift+1)*Psi2));
    p_mTpcQ2EastSin[VzSign][i_shift]->Fill((double)RunIndex,(double)Cent9,TMath::Sin(2.0*(i_shift+1)*Psi2));
  }
}

void StEventPlaneProManager::fillTpcShiftWest(TVector2 qVector, int Cent9, int RunIndex, int VzSign)
{
  float Psi2 = TMath::ATan2(qVector.Y(),qVector.X())/2.0;
  for(int i_shift = 0; i_shift < recoEP::mNumShiftOrder; ++i_shift)
  {
    p_mTpcQ2WestCos[VzSign][i_shift]->Fill((double)RunIndex,(double)Cent9,TMath::Cos(2.0*(i_shift+1)*Psi2));
    p_mTpcQ2WestSin[VzSign][i_shift]->Fill((double)RunIndex,(double)Cent9,TMath::Sin(2.0*(i_shift+1)*Psi2));
  }
}

void StEventPlaneProManager::fillTpcShiftFull(TVector2 qVector, int Cent9, int RunIndex, int VzSign)
{
  float Psi2 = TMath::ATan2(qVector.Y(),qVector.X())/2.0;
  for(int i_shift = 0; i_shift < recoEP::mNumShiftOrder; ++i_shift)
  {
    p_mTpcQ2FullCos[VzSign][i_shift]->Fill((double)RunIndex,(double)Cent9,TMath::Cos(2.0*(i_shift+1)*Psi2));
    p_mTpcQ2FullSin[VzSign][i_shift]->Fill((double)RunIndex,(double)Cent9,TMath::Sin(2.0*(i_shift+1)*Psi2));
  }
}

void StEventPlaneProManager::writeTpcShift()
{
  for(int i_vz = 0; i_vz < 2; ++i_vz)
  {
    for(int i_shift = 0; i_shift < recoEP::mNumShiftOrder; ++i_shift) // Shift Order
    {
      p_mTpcQ2EastCos[i_vz][i_shift]->Write();
      p_mTpcQ2EastSin[i_vz][i_shift]->Write();
      p_mTpcQ2WestCos[i_vz][i_shift]->Write();
      p_mTpcQ2WestSin[i_vz][i_shift]->Write();
      p_mTpcQ2FullCos[i_vz][i_shift]->Write();
      p_mTpcQ2FullSin[i_vz][i_shift]->Write();
    }
  }
}

// TPC Sub/Ran EP Resolution
void StEventPlaneProManager::initTpcResolution()
{
  p_mTpcSubRes2QA = new TProfile2D("p_mTpcSubRes2QA","p_mTpcSubRes2QA",recoEP::mNumOfRunIndex,-0.5,(float)recoEP::mNumOfRunIndex-0.5,9,-0.5,8.5);
  p_mTpcRanRes2QA = new TProfile2D("p_mTpcRanRes2QA","p_mTpcRanRes2QA",recoEP::mNumOfRunIndex,-0.5,(float)recoEP::mNumOfRunIndex-0.5,9,-0.5,8.5);
  p_mTpcSubRes2 = new TProfile("p_mTpcSubRes2","p_mTpcSubRes2",9,-0.5,8.5);
  p_mTpcRanRes2 = new TProfile("p_mTpcRanRes2","p_mTpcRanRes2",9,-0.5,8.5);
}

void StEventPlaneProManager::fillTpcResSub(double PsiEast, double PsiWest, int Cent9, int RunIndex)
{
  double res2 = TMath::Cos(2.0*(PsiWest-PsiEast));
  p_mTpcSubRes2QA->Fill((double)RunIndex,(double)Cent9,res2);
  p_mTpcSubRes2->Fill((double)Cent9,res2);
}

void StEventPlaneProManager::fillTpcResRan(double PsiRanA, double PsiRanB, int Cent9, int RunIndex)
{
  double res2 = TMath::Cos(2.0*(PsiRanB-PsiRanA));
  p_mTpcRanRes2QA->Fill((double)RunIndex,(double)Cent9,res2);
  p_mTpcRanRes2->Fill((double)Cent9,res2);
}

void StEventPlaneProManager::writeTpcResolution()
{
  p_mTpcSubRes2QA->Write();
  p_mTpcRanRes2QA->Write();
  p_mTpcSubRes2->Write();
  p_mTpcRanRes2->Write();
}
//---------------------------------------------------------------------------------

//---------------------------------------------------------------------------------
// Charged Particle Flow
void StEventPlaneProManager::initChargedFlow()
{
  string ProName;
  for(int i_cent = 0; i_cent < 10; ++i_cent)
  {
    ProName = Form("p_mChargedV1PpQA_Cent%d",i_cent);
    p_mChargedV1PpQA[i_cent] = new TProfile(ProName.c_str(),ProName.c_str(),recoEP::mNumOfRunIndex,-0.5,(float)recoEP::mNumOfRunIndex-0.5);
    ProName = Form("p_mChargedV2EpQA_Cent%d",i_cent);
    p_mChargedV2EpQA[i_cent] = new TProfile(ProName.c_str(),ProName.c_str(),recoEP::mNumOfRunIndex,-0.5,(float)recoEP::mNumOfRunIndex-0.5);
    ProName = Form("p_mChargedV2PpQA_Cent%d",i_cent);
    p_mChargedV2PpQA[i_cent] = new TProfile(ProName.c_str(),ProName.c_str(),recoEP::mNumOfRunIndex,-0.5,(float)recoEP::mNumOfRunIndex-0.5);
    ProName = Form("p_mChargedV3EpQA_Cent%d",i_cent);
    p_mChargedV3EpQA[i_cent] = new TProfile(ProName.c_str(),ProName.c_str(),recoEP::mNumOfRunIndex,-0.5,(float)recoEP::mNumOfRunIndex-0.5);
    ProName = Form("p_mChargedV1Pp_Cent%d",i_cent);
    p_mChargedV1Pp[i_cent] = new TProfile(ProName.c_str(),ProName.c_str(),10,-1.0,1.0);
    ProName = Form("p_mChargedV2Ep_Cent%d",i_cent);
    p_mChargedV2Ep[i_cent] = new TProfile(ProName.c_str(),ProName.c_str(),25,0.2,4.2);
    ProName = Form("p_mChargedV2Pp_Cent%d",i_cent);
    p_mChargedV2Pp[i_cent] = new TProfile(ProName.c_str(),ProName.c_str(),25,0.2,4.2);
    ProName = Form("p_mChargedV3Ep_Cent%d",i_cent);
    p_mChargedV3Ep[i_cent] = new TProfile(ProName.c_str(),ProName.c_str(),25,0.2,4.2);
  }
}

void StEventPlaneProManager::fillChargedV1Pp(double pt, double eta, double v1, double res1, int Cent9, int RunIndex, double reweight)
{
  if(res1 > 0 && pt >= 0.2 && pt <= 2.0)
  {
    p_mChargedV1PpQA[Cent9]->Fill(RunIndex,v1/res1,reweight);
    p_mChargedV1PpQA[9]->Fill(RunIndex,v1/res1,reweight);

    p_mChargedV1Pp[Cent9]->Fill(eta,v1/res1,reweight);
    p_mChargedV1Pp[9]->Fill(eta,v1/res1,reweight);
  }
}

void StEventPlaneProManager::fillChargedV2Pp(double pt, double v2, double res2, int Cent9, int RunIndex, double reweight)
{
  if(res2 > 0.0)
  {
    if(pt >= 0.2 && pt <= 5.2)
    {
      p_mChargedV2PpQA[Cent9]->Fill(RunIndex,v2/res2,reweight);
      p_mChargedV2PpQA[9]->Fill(RunIndex,v2/res2,reweight);
    }

    p_mChargedV2Pp[Cent9]->Fill(pt,v2/res2,reweight);
    p_mChargedV2Pp[9]->Fill(pt,v2/res2,reweight);
  }
}

void StEventPlaneProManager::fillChargedV2Ep(double pt, double v2, double res2, int Cent9, int RunIndex, double reweight)
{
  if(res2 > 0.0)
  {
    if(pt >= 0.2 && pt <= 5.2)
    {
      p_mChargedV2EpQA[Cent9]->Fill(RunIndex,v2/res2,reweight);
      p_mChargedV2EpQA[9]->Fill(RunIndex,v2/res2,reweight);
    }

    p_mChargedV2Ep[Cent9]->Fill(pt,v2/res2,reweight);
    p_mChargedV2Ep[9]->Fill(pt,v2/res2,reweight);
  }
}

void StEventPlaneProManager::fillChargedV3Ep(double pt, double v3, double res3, int Cent9, int RunIndex, double reweight)
{
  if(res3 > 0.0)
  {
    if(pt >= 0.2 && pt <= 5.2)
    {
      p_mChargedV3EpQA[Cent9]->Fill(RunIndex,v3/res3,reweight);
      p_mChargedV3EpQA[9]->Fill(RunIndex,v3/res3,reweight);
    }

    p_mChargedV3Ep[Cent9]->Fill(pt,v3/res3,reweight);
    p_mChargedV3Ep[9]->Fill(pt,v3/res3,reweight);
  }
}

void StEventPlaneProManager::writeChargedFlow()
{
  for(int i_cent = 0; i_cent < 10; ++i_cent)
  {
    p_mChargedV1PpQA[i_cent]->Write();
    p_mChargedV2EpQA[i_cent]->Write();
    p_mChargedV2PpQA[i_cent]->Write();
    p_mChargedV3EpQA[i_cent]->Write();
    p_mChargedV1Pp[i_cent]->Write();
    p_mChargedV2Ep[i_cent]->Write();
    p_mChargedV2Pp[i_cent]->Write();
    p_mChargedV3Ep[i_cent]->Write();
  }
}

//---------------------------------------------------------------------------------
