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
    ProName = Form("p_mTpcq1xEast_%s",recoEP::mVStr[i_vz].c_str());
    p_mTpcq1xEast[i_vz] = new TProfile2D(ProName.c_str(),ProName.c_str(),recoEP::mNumOfRunIndex,-0.5,(float)recoEP::mNumOfRunIndex-0.5,9,-0.5,8.5); // neg eta
    ProName = Form("p_mTpcq1yEast_%s",recoEP::mVStr[i_vz].c_str());
    p_mTpcq1yEast[i_vz] = new TProfile2D(ProName.c_str(),ProName.c_str(),recoEP::mNumOfRunIndex,-0.5,(float)recoEP::mNumOfRunIndex-0.5,9,-0.5,8.5);
    ProName = Form("p_mTpcq1xWest_%s",recoEP::mVStr[i_vz].c_str());
    p_mTpcq1xWest[i_vz] = new TProfile2D(ProName.c_str(),ProName.c_str(),recoEP::mNumOfRunIndex,-0.5,(float)recoEP::mNumOfRunIndex-0.5,9,-0.5,8.5); // pos eta
    ProName = Form("p_mTpcq1yWest_%s",recoEP::mVStr[i_vz].c_str());
    p_mTpcq1yWest[i_vz] = new TProfile2D(ProName.c_str(),ProName.c_str(),recoEP::mNumOfRunIndex,-0.5,(float)recoEP::mNumOfRunIndex-0.5,9,-0.5,8.5);

    ProName = Form("p_mTpcq1xFull_%s",recoEP::mVStr[i_vz].c_str());
    p_mTpcq1xFull[i_vz] = new TProfile2D(ProName.c_str(),ProName.c_str(),recoEP::mNumOfRunIndex,-0.5,(float)recoEP::mNumOfRunIndex-0.5,9,-0.5,8.5);
    ProName = Form("p_mTpcq1yFull_%s",recoEP::mVStr[i_vz].c_str());
    p_mTpcq1yFull[i_vz] = new TProfile2D(ProName.c_str(),ProName.c_str(),recoEP::mNumOfRunIndex,-0.5,(float)recoEP::mNumOfRunIndex-0.5,9,-0.5,8.5);

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
  
    ProName = Form("p_mTpcq3xEast_%s",recoEP::mVStr[i_vz].c_str());
    p_mTpcq3xEast[i_vz] = new TProfile2D(ProName.c_str(),ProName.c_str(),recoEP::mNumOfRunIndex,-0.5,(float)recoEP::mNumOfRunIndex-0.5,9,-0.5,8.5); // neg eta
    ProName = Form("p_mTpcq3yEast_%s",recoEP::mVStr[i_vz].c_str());
    p_mTpcq3yEast[i_vz] = new TProfile2D(ProName.c_str(),ProName.c_str(),recoEP::mNumOfRunIndex,-0.5,(float)recoEP::mNumOfRunIndex-0.5,9,-0.5,8.5);
    ProName = Form("p_mTpcq3xWest_%s",recoEP::mVStr[i_vz].c_str());
    p_mTpcq3xWest[i_vz] = new TProfile2D(ProName.c_str(),ProName.c_str(),recoEP::mNumOfRunIndex,-0.5,(float)recoEP::mNumOfRunIndex-0.5,9,-0.5,8.5); // pos eta
    ProName = Form("p_mTpcq3yWest_%s",recoEP::mVStr[i_vz].c_str());
    p_mTpcq3yWest[i_vz] = new TProfile2D(ProName.c_str(),ProName.c_str(),recoEP::mNumOfRunIndex,-0.5,(float)recoEP::mNumOfRunIndex-0.5,9,-0.5,8.5);

    ProName = Form("p_mTpcq3xFull_%s",recoEP::mVStr[i_vz].c_str());
    p_mTpcq3xFull[i_vz] = new TProfile2D(ProName.c_str(),ProName.c_str(),recoEP::mNumOfRunIndex,-0.5,(float)recoEP::mNumOfRunIndex-0.5,9,-0.5,8.5);
    ProName = Form("p_mTpcq3yFull_%s",recoEP::mVStr[i_vz].c_str());
    p_mTpcq3yFull[i_vz] = new TProfile2D(ProName.c_str(),ProName.c_str(),recoEP::mNumOfRunIndex,-0.5,(float)recoEP::mNumOfRunIndex-0.5,9,-0.5,8.5);

  }
}

void StEventPlaneProManager::fillTpcReCenterEast(int order, TVector2 qVector, int Cent9, int RunIndex, int VzSign, double pt)
{
  const double qx = qVector.X();
  const double qy = qVector.Y();

  double weight;
  if(pt <= recoEP::mPrimPtWeight)
  {
    weight = pt;
  }
  if(pt > recoEP::mPrimPtWeight)
  {
    weight = recoEP::mPrimPtWeight;
  }
  
  if(order == 1)
  {
    p_mTpcq1xEast[VzSign]->Fill((double)RunIndex,(double)Cent9,(double)qx,(double)weight);
    p_mTpcq1yEast[VzSign]->Fill((double)RunIndex,(double)Cent9,(double)qy,(double)weight);
  }
  if(order == 2)
  {
    p_mTpcq2xEast[VzSign]->Fill((double)RunIndex,(double)Cent9,(double)qx,(double)weight);
    p_mTpcq2yEast[VzSign]->Fill((double)RunIndex,(double)Cent9,(double)qy,(double)weight);
  }
  if(order == 3)
  {
    p_mTpcq3xEast[VzSign]->Fill((double)RunIndex,(double)Cent9,(double)qx,(double)weight);
    p_mTpcq3yEast[VzSign]->Fill((double)RunIndex,(double)Cent9,(double)qy,(double)weight);
  }
}

void StEventPlaneProManager::fillTpcReCenterWest(int order, TVector2 qVector, int Cent9, int RunIndex, int VzSign, double pt)
{
  const double qx = qVector.X();
  const double qy = qVector.Y();

  double weight;
  if(pt <= recoEP::mPrimPtWeight)
  {
    weight = pt;
  }
  if(pt > recoEP::mPrimPtWeight)
  {
    weight = recoEP::mPrimPtWeight;
  }
  if(order == 1)
  {
    p_mTpcq1xWest[VzSign]->Fill((double)RunIndex,(double)Cent9,(double)qx,(double)weight);
    p_mTpcq1yWest[VzSign]->Fill((double)RunIndex,(double)Cent9,(double)qy,(double)weight);
  }
  if(order == 2)
  {
    p_mTpcq2xWest[VzSign]->Fill((double)RunIndex,(double)Cent9,(double)qx,(double)weight);
    p_mTpcq2yWest[VzSign]->Fill((double)RunIndex,(double)Cent9,(double)qy,(double)weight);
  }
  if(order == 3)
  {
    p_mTpcq3xWest[VzSign]->Fill((double)RunIndex,(double)Cent9,(double)qx,(double)weight);
    p_mTpcq3yWest[VzSign]->Fill((double)RunIndex,(double)Cent9,(double)qy,(double)weight);
  }
}

void StEventPlaneProManager::fillTpcReCenterFull(int order, TVector2 qVector, int Cent9, int RunIndex, int VzSign, double pt)
{
  const double qx = qVector.X();
  const double qy = qVector.Y();

  double weight;
  if(pt <= recoEP::mPrimPtWeight)
  {
    weight = pt;
  }
  if(pt > recoEP::mPrimPtWeight)
  {
    weight = recoEP::mPrimPtWeight;
  }
  if(order == 1)
  {
    p_mTpcq1xFull[VzSign]->Fill((double)RunIndex,(double)Cent9,(double)qx,(double)weight);
    p_mTpcq1yFull[VzSign]->Fill((double)RunIndex,(double)Cent9,(double)qy,(double)weight);
  }
  if(order == 2)
  {
    p_mTpcq2xFull[VzSign]->Fill((double)RunIndex,(double)Cent9,(double)qx,(double)weight);
    p_mTpcq2yFull[VzSign]->Fill((double)RunIndex,(double)Cent9,(double)qy,(double)weight);
  }
  if(order == 3)
  {
    p_mTpcq3xFull[VzSign]->Fill((double)RunIndex,(double)Cent9,(double)qx,(double)weight);
    p_mTpcq3yFull[VzSign]->Fill((double)RunIndex,(double)Cent9,(double)qy,(double)weight);
  }
}

void StEventPlaneProManager::writeTpcReCenter()
{
  for(int i_vz = 0; i_vz < 2; ++i_vz) // vertex pos/neg
  {
    p_mTpcq1xEast[i_vz]->Write();
    p_mTpcq1yEast[i_vz]->Write();
    p_mTpcq1xWest[i_vz]->Write();
    p_mTpcq1yWest[i_vz]->Write();

    p_mTpcq1xFull[i_vz]->Write();
    p_mTpcq1yFull[i_vz]->Write();
    
    p_mTpcq2xEast[i_vz]->Write();
    p_mTpcq2yEast[i_vz]->Write();
    p_mTpcq2xWest[i_vz]->Write();
    p_mTpcq2yWest[i_vz]->Write();

    p_mTpcq2xFull[i_vz]->Write();
    p_mTpcq2yFull[i_vz]->Write();
    
    p_mTpcq3xEast[i_vz]->Write();
    p_mTpcq3yEast[i_vz]->Write();
    p_mTpcq3xWest[i_vz]->Write();
    p_mTpcq3yWest[i_vz]->Write();

    p_mTpcq3xFull[i_vz]->Write();
    p_mTpcq3yFull[i_vz]->Write();
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
      ProName = Form("p_mTpcQ1EastCos_%s_%d",recoEP::mVStr[i_vz].c_str(),i_shift);
      p_mTpcQ1EastCos[i_vz][i_shift] = new TProfile2D(ProName.c_str(),ProName.c_str(),recoEP::mNumOfRunIndex,-0.5,(float)recoEP::mNumOfRunIndex-0.5,9,-0.5,8.5);
      ProName = Form("p_mTpcQ1EastSin_%s_%d",recoEP::mVStr[i_vz].c_str(),i_shift);
      p_mTpcQ1EastSin[i_vz][i_shift] = new TProfile2D(ProName.c_str(),ProName.c_str(),recoEP::mNumOfRunIndex,-0.5,(float)recoEP::mNumOfRunIndex-0.5,9,-0.5,8.5);

      ProName = Form("p_mTpcQ1WestCos_%s_%d",recoEP::mVStr[i_vz].c_str(),i_shift);
      p_mTpcQ1WestCos[i_vz][i_shift] = new TProfile2D(ProName.c_str(),ProName.c_str(),recoEP::mNumOfRunIndex,-0.5,(float)recoEP::mNumOfRunIndex-0.5,9,-0.5,8.5);
      ProName = Form("p_mTpcQ1WestSin_%s_%d",recoEP::mVStr[i_vz].c_str(),i_shift);
      p_mTpcQ1WestSin[i_vz][i_shift] = new TProfile2D(ProName.c_str(),ProName.c_str(),recoEP::mNumOfRunIndex,-0.5,(float)recoEP::mNumOfRunIndex-0.5,9,-0.5,8.5);

      ProName = Form("p_mTpcQ1FullCos_%s_%d",recoEP::mVStr[i_vz].c_str(),i_shift);
      p_mTpcQ1FullCos[i_vz][i_shift] = new TProfile2D(ProName.c_str(),ProName.c_str(),recoEP::mNumOfRunIndex,-0.5,(float)recoEP::mNumOfRunIndex-0.5,9,-0.5,8.5);
      ProName = Form("p_mTpcQ1FullSin_%s_%d",recoEP::mVStr[i_vz].c_str(),i_shift);
      p_mTpcQ1FullSin[i_vz][i_shift] = new TProfile2D(ProName.c_str(),ProName.c_str(),recoEP::mNumOfRunIndex,-0.5,(float)recoEP::mNumOfRunIndex-0.5,9,-0.5,8.5);
    

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
    
   
      ProName = Form("p_mTpcQ3EastCos_%s_%d",recoEP::mVStr[i_vz].c_str(),i_shift);
      p_mTpcQ3EastCos[i_vz][i_shift] = new TProfile2D(ProName.c_str(),ProName.c_str(),recoEP::mNumOfRunIndex,-0.5,(float)recoEP::mNumOfRunIndex-0.5,9,-0.5,8.5);
      ProName = Form("p_mTpcQ3EastSin_%s_%d",recoEP::mVStr[i_vz].c_str(),i_shift);
      p_mTpcQ3EastSin[i_vz][i_shift] = new TProfile2D(ProName.c_str(),ProName.c_str(),recoEP::mNumOfRunIndex,-0.5,(float)recoEP::mNumOfRunIndex-0.5,9,-0.5,8.5);

      ProName = Form("p_mTpcQ3WestCos_%s_%d",recoEP::mVStr[i_vz].c_str(),i_shift);
      p_mTpcQ3WestCos[i_vz][i_shift] = new TProfile2D(ProName.c_str(),ProName.c_str(),recoEP::mNumOfRunIndex,-0.5,(float)recoEP::mNumOfRunIndex-0.5,9,-0.5,8.5);
      ProName = Form("p_mTpcQ3WestSin_%s_%d",recoEP::mVStr[i_vz].c_str(),i_shift);
      p_mTpcQ3WestSin[i_vz][i_shift] = new TProfile2D(ProName.c_str(),ProName.c_str(),recoEP::mNumOfRunIndex,-0.5,(float)recoEP::mNumOfRunIndex-0.5,9,-0.5,8.5);

      ProName = Form("p_mTpcQ3FullCos_%s_%d",recoEP::mVStr[i_vz].c_str(),i_shift);
      p_mTpcQ3FullCos[i_vz][i_shift] = new TProfile2D(ProName.c_str(),ProName.c_str(),recoEP::mNumOfRunIndex,-0.5,(float)recoEP::mNumOfRunIndex-0.5,9,-0.5,8.5);
      ProName = Form("p_mTpcQ3FullSin_%s_%d",recoEP::mVStr[i_vz].c_str(),i_shift);
      p_mTpcQ3FullSin[i_vz][i_shift] = new TProfile2D(ProName.c_str(),ProName.c_str(),recoEP::mNumOfRunIndex,-0.5,(float)recoEP::mNumOfRunIndex-0.5,9,-0.5,8.5);
    }
  }
}

void StEventPlaneProManager::fillTpcShiftEast(int order, TVector2 qVector, int Cent9, int RunIndex, int VzSign)
{
  double Psi = TMath::ATan2(qVector.Y(),qVector.X())/double(order);
  if(order == 1)
  {
    for(int i_shift = 0; i_shift < recoEP::mNumShiftOrder; ++i_shift)
    {
      p_mTpcQ1EastCos[VzSign][i_shift]->Fill((double)RunIndex,(double)Cent9,TMath::Cos(double(order)*(i_shift+1)*Psi));
      p_mTpcQ1EastSin[VzSign][i_shift]->Fill((double)RunIndex,(double)Cent9,TMath::Sin(double(order)*(i_shift+1)*Psi));
    }
  }
  if(order == 2)
  {
    for(int i_shift = 0; i_shift < recoEP::mNumShiftOrder; ++i_shift)
    {
      p_mTpcQ2EastCos[VzSign][i_shift]->Fill((double)RunIndex,(double)Cent9,TMath::Cos(double(order)*(i_shift+1)*Psi));
      p_mTpcQ2EastSin[VzSign][i_shift]->Fill((double)RunIndex,(double)Cent9,TMath::Sin(double(order)*(i_shift+1)*Psi));
    }
  }
  if(order == 3)
  {
    for(int i_shift = 0; i_shift < recoEP::mNumShiftOrder; ++i_shift)
    {
      p_mTpcQ3EastCos[VzSign][i_shift]->Fill((double)RunIndex,(double)Cent9,TMath::Cos(double(order)*(i_shift+1)*Psi));
      p_mTpcQ3EastSin[VzSign][i_shift]->Fill((double)RunIndex,(double)Cent9,TMath::Sin(double(order)*(i_shift+1)*Psi));
    }
  }
}

void StEventPlaneProManager::fillTpcShiftWest(int order, TVector2 qVector, int Cent9, int RunIndex, int VzSign)
{
  float Psi = TMath::ATan2(qVector.Y(),qVector.X())/double(order);
  if(order == 1)
  {
    for(int i_shift = 0; i_shift < recoEP::mNumShiftOrder; ++i_shift)
    {
      p_mTpcQ1WestCos[VzSign][i_shift]->Fill((double)RunIndex,(double)Cent9,TMath::Cos(double(order)*(i_shift+1)*Psi));
      p_mTpcQ1WestSin[VzSign][i_shift]->Fill((double)RunIndex,(double)Cent9,TMath::Sin(double(order)*(i_shift+1)*Psi));
    }
  }
  if(order == 2)
  {
    for(int i_shift = 0; i_shift < recoEP::mNumShiftOrder; ++i_shift)
    {
      p_mTpcQ2WestCos[VzSign][i_shift]->Fill((double)RunIndex,(double)Cent9,TMath::Cos(double(order)*(i_shift+1)*Psi));
      p_mTpcQ2WestSin[VzSign][i_shift]->Fill((double)RunIndex,(double)Cent9,TMath::Sin(double(order)*(i_shift+1)*Psi));
    }
  }
  if(order == 3)
  {
    for(int i_shift = 0; i_shift < recoEP::mNumShiftOrder; ++i_shift)
    {
      p_mTpcQ3WestCos[VzSign][i_shift]->Fill((double)RunIndex,(double)Cent9,TMath::Cos(double(order)*(i_shift+1)*Psi));
      p_mTpcQ3WestSin[VzSign][i_shift]->Fill((double)RunIndex,(double)Cent9,TMath::Sin(double(order)*(i_shift+1)*Psi));
    }
  }
}

void StEventPlaneProManager::fillTpcShiftFull(int order, TVector2 qVector, int Cent9, int RunIndex, int VzSign)
{
  float Psi = TMath::ATan2(qVector.Y(),qVector.X())/double(order);
  if(order == 1)
  {
    for(int i_shift = 0; i_shift < recoEP::mNumShiftOrder; ++i_shift)
    {
      p_mTpcQ1FullCos[VzSign][i_shift]->Fill((double)RunIndex,(double)Cent9,TMath::Cos(double(order)*(i_shift+1)*Psi));
      p_mTpcQ1FullSin[VzSign][i_shift]->Fill((double)RunIndex,(double)Cent9,TMath::Sin(double(order)*(i_shift+1)*Psi));
    }
  }
  if(order == 2)
  {
    for(int i_shift = 0; i_shift < recoEP::mNumShiftOrder; ++i_shift)
    {
      p_mTpcQ2FullCos[VzSign][i_shift]->Fill((double)RunIndex,(double)Cent9,TMath::Cos(double(order)*(i_shift+1)*Psi));
      p_mTpcQ2FullSin[VzSign][i_shift]->Fill((double)RunIndex,(double)Cent9,TMath::Sin(double(order)*(i_shift+1)*Psi));
    }
  }
  if(order == 3)
  {
    for(int i_shift = 0; i_shift < recoEP::mNumShiftOrder; ++i_shift)
    {
      p_mTpcQ3FullCos[VzSign][i_shift]->Fill((double)RunIndex,(double)Cent9,TMath::Cos(double(order)*(i_shift+1)*Psi));
      p_mTpcQ3FullSin[VzSign][i_shift]->Fill((double)RunIndex,(double)Cent9,TMath::Sin(double(order)*(i_shift+1)*Psi));
    }
  }
}

void StEventPlaneProManager::writeTpcShift()
{
  for(int i_vz = 0; i_vz < 2; ++i_vz)
  {
    for(int i_shift = 0; i_shift < recoEP::mNumShiftOrder; ++i_shift) // Shift Order
    {
      p_mTpcQ1EastCos[i_vz][i_shift]->Write();
      p_mTpcQ1EastSin[i_vz][i_shift]->Write();
      p_mTpcQ1WestCos[i_vz][i_shift]->Write();
      p_mTpcQ1WestSin[i_vz][i_shift]->Write();
      p_mTpcQ1FullCos[i_vz][i_shift]->Write();
      p_mTpcQ1FullSin[i_vz][i_shift]->Write();

      p_mTpcQ2EastCos[i_vz][i_shift]->Write();
      p_mTpcQ2EastSin[i_vz][i_shift]->Write();
      p_mTpcQ2WestCos[i_vz][i_shift]->Write();
      p_mTpcQ2WestSin[i_vz][i_shift]->Write();
      p_mTpcQ2FullCos[i_vz][i_shift]->Write();
      p_mTpcQ2FullSin[i_vz][i_shift]->Write();
     
      p_mTpcQ3EastCos[i_vz][i_shift]->Write();
      p_mTpcQ3EastSin[i_vz][i_shift]->Write();
      p_mTpcQ3WestCos[i_vz][i_shift]->Write();
      p_mTpcQ3WestSin[i_vz][i_shift]->Write();
      p_mTpcQ3FullCos[i_vz][i_shift]->Write();
      p_mTpcQ3FullSin[i_vz][i_shift]->Write();
    }
  }
}

// TPC Sub/Ran EP Resolution
void StEventPlaneProManager::initTpcResolution()
{
  p_mTpcSubRes1QA = new TProfile2D("p_mTpcSubRes1QA","p_mTpcSubRes1QA",recoEP::mNumOfRunIndex,-0.5,(float)recoEP::mNumOfRunIndex-0.5,9,-0.5,8.5);
  p_mTpcRanRes1QA = new TProfile2D("p_mTpcRanRes1QA","p_mTpcRanRes1QA",recoEP::mNumOfRunIndex,-0.5,(float)recoEP::mNumOfRunIndex-0.5,9,-0.5,8.5);
  p_mTpcSubRes1 = new TProfile("p_mTpcSubRes1","p_mTpcSubRes1",9,-0.5,8.5);
  p_mTpcRanRes1 = new TProfile("p_mTpcRanRes1","p_mTpcRanRes1",9,-0.5,8.5);

  p_mTpcSubRes2QA = new TProfile2D("p_mTpcSubRes2QA","p_mTpcSubRes2QA",recoEP::mNumOfRunIndex,-0.5,(float)recoEP::mNumOfRunIndex-0.5,9,-0.5,8.5);
  p_mTpcRanRes2QA = new TProfile2D("p_mTpcRanRes2QA","p_mTpcRanRes2QA",recoEP::mNumOfRunIndex,-0.5,(float)recoEP::mNumOfRunIndex-0.5,9,-0.5,8.5);
  p_mTpcSubRes2 = new TProfile("p_mTpcSubRes2","p_mTpcSubRes2",9,-0.5,8.5);
  p_mTpcRanRes2 = new TProfile("p_mTpcRanRes2","p_mTpcRanRes2",9,-0.5,8.5);

  p_mTpcSubRes3QA = new TProfile2D("p_mTpcSubRes3QA","p_mTpcSubRes3QA",recoEP::mNumOfRunIndex,-0.5,(float)recoEP::mNumOfRunIndex-0.5,9,-0.5,8.5);
  p_mTpcRanRes3QA = new TProfile2D("p_mTpcRanRes3QA","p_mTpcRanRes3QA",recoEP::mNumOfRunIndex,-0.5,(float)recoEP::mNumOfRunIndex-0.5,9,-0.5,8.5);
  p_mTpcSubRes3 = new TProfile("p_mTpcSubRes3","p_mTpcSubRes3",9,-0.5,8.5);
  p_mTpcRanRes3 = new TProfile("p_mTpcRanRes3","p_mTpcRanRes3",9,-0.5,8.5);

}

void StEventPlaneProManager::fillTpcResSub(int order, double PsiEast, double PsiWest, int Cent9, int RunIndex)
{
  double res = TMath::Cos(double(order)*(PsiWest-PsiEast));
  if(order == 1)
  {
    p_mTpcSubRes1QA->Fill((double)RunIndex,(double)Cent9,res);
    p_mTpcSubRes1->Fill((double)Cent9,res);
  }
  if(order == 2)
  {
    p_mTpcSubRes2QA->Fill((double)RunIndex,(double)Cent9,res);
    p_mTpcSubRes2->Fill((double)Cent9,res);
  }
  if(order == 3)
  {
    p_mTpcSubRes3QA->Fill((double)RunIndex,(double)Cent9,res);
    p_mTpcSubRes3->Fill((double)Cent9,res);
  }
}

void StEventPlaneProManager::fillTpcResRan(int order, double PsiRanA, double PsiRanB, int Cent9, int RunIndex)
{
  double res = TMath::Cos(double(order)*(PsiRanB-PsiRanA));
  if(order == 1)
  {
    p_mTpcRanRes1QA->Fill((double)RunIndex,(double)Cent9,res);
    p_mTpcRanRes1->Fill((double)Cent9,res);
  }
  if(order == 2)
  {
    p_mTpcRanRes2QA->Fill((double)RunIndex,(double)Cent9,res);
    p_mTpcRanRes2->Fill((double)Cent9,res);
  }
  if(order == 3)
  {
    p_mTpcRanRes3QA->Fill((double)RunIndex,(double)Cent9,res);
    p_mTpcRanRes3->Fill((double)Cent9,res);
  }
}

void StEventPlaneProManager::writeTpcResolution()
{
  p_mTpcSubRes1QA->Write();
  p_mTpcRanRes1QA->Write();
  p_mTpcSubRes1->Write();
  p_mTpcRanRes1->Write();

  p_mTpcSubRes2QA->Write();
  p_mTpcRanRes2QA->Write();
  p_mTpcSubRes2->Write();
  p_mTpcRanRes2->Write();

  p_mTpcSubRes3QA->Write();
  p_mTpcRanRes3QA->Write();
  p_mTpcSubRes3->Write();
  p_mTpcRanRes3->Write();
}

void StEventPlaneProManager::initTpcFlowEta()
{
  for (int order=1; order<=3; order++)
  { 
    for (int icent=0; icent < 9; ++icent)
    { 
      p_mTpcFlowEta[order-1][icent] = new TProfile(Form("p_mTpcFlowEta%d_Cent%d",order,icent),Form("p_mTpcFlowEta%d_Cent%d",order,icent),30,-5.5,5.5);
    }
  }
}

void StEventPlaneProManager::fillTpcFlowEta(double eta, double v, int Cent9, int order, double weight)
{
  p_mTpcFlowEta[order-1][Cent9]->Fill(eta,v,weight);
}

void StEventPlaneProManager::writeTpcFlowEta()
{
  for (int order=1; order<=3; order++)
  { 
    for (int icent=0; icent < 9; ++icent)
    {
      p_mTpcFlowEta[order-1][icent]->Write();
    }
  }
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
    ProName = Form("p_mChargedV1EpQA_Cent%d",i_cent);
    p_mChargedV1EpQA[i_cent] = new TProfile(ProName.c_str(),ProName.c_str(),recoEP::mNumOfRunIndex,-0.5,(float)recoEP::mNumOfRunIndex-0.5);
    ProName = Form("p_mChargedV1EpdQA_Cent%d",i_cent);
    p_mChargedV1EpdQA[i_cent] = new TProfile(ProName.c_str(),ProName.c_str(),recoEP::mNumOfRunIndex,-0.5,(float)recoEP::mNumOfRunIndex-0.5);
    ProName = Form("p_mChargedV2EpQA_Cent%d",i_cent);
    p_mChargedV2EpQA[i_cent] = new TProfile(ProName.c_str(),ProName.c_str(),recoEP::mNumOfRunIndex,-0.5,(float)recoEP::mNumOfRunIndex-0.5);
    ProName = Form("p_mChargedV2PpQA_Cent%d",i_cent);
    p_mChargedV2PpQA[i_cent] = new TProfile(ProName.c_str(),ProName.c_str(),recoEP::mNumOfRunIndex,-0.5,(float)recoEP::mNumOfRunIndex-0.5);
    ProName = Form("p_mChargedV2EpdQA_Cent%d",i_cent);
    p_mChargedV2EpdQA[i_cent] = new TProfile(ProName.c_str(),ProName.c_str(),recoEP::mNumOfRunIndex,-0.5,(float)recoEP::mNumOfRunIndex-0.5);
    ProName = Form("p_mChargedV3EpQA_Cent%d",i_cent);
    p_mChargedV3EpQA[i_cent] = new TProfile(ProName.c_str(),ProName.c_str(),recoEP::mNumOfRunIndex,-0.5,(float)recoEP::mNumOfRunIndex-0.5);
    ProName = Form("p_mChargedV3EpdQA_Cent%d",i_cent);
    p_mChargedV3EpdQA[i_cent] = new TProfile(ProName.c_str(),ProName.c_str(),recoEP::mNumOfRunIndex,-0.5,(float)recoEP::mNumOfRunIndex-0.5);
    ProName = Form("p_mChargedV1Pp_Cent%d",i_cent);
    p_mChargedV1Pp[i_cent] = new TProfile(ProName.c_str(),ProName.c_str(),10,-1.0,1.0);
    ProName = Form("p_mChargedV1Ep_Cent%d",i_cent);
    p_mChargedV1Ep[i_cent] = new TProfile(ProName.c_str(),ProName.c_str(),25,0.2,4.2);
    ProName = Form("p_mChargedV1Epd_Cent%d",i_cent);
    p_mChargedV1Epd[i_cent] = new TProfile(ProName.c_str(),ProName.c_str(),recoEP::mNumOfRunIndex,-0.5,(float)recoEP::mNumOfRunIndex-0.5);
    ProName = Form("p_mChargedV2Ep_Cent%d",i_cent);
    p_mChargedV2Ep[i_cent] = new TProfile(ProName.c_str(),ProName.c_str(),25,0.2,4.2);
    ProName = Form("p_mChargedV2Pp_Cent%d",i_cent);
    p_mChargedV2Pp[i_cent] = new TProfile(ProName.c_str(),ProName.c_str(),25,0.2,4.2);
    ProName = Form("p_mChargedV2Epd_Cent%d",i_cent);
    p_mChargedV2Epd[i_cent] = new TProfile(ProName.c_str(),ProName.c_str(),recoEP::mNumOfRunIndex,-0.5,(float)recoEP::mNumOfRunIndex-0.5);
    ProName = Form("p_mChargedV3Ep_Cent%d",i_cent);
    p_mChargedV3Ep[i_cent] = new TProfile(ProName.c_str(),ProName.c_str(),25,0.2,4.2);
    ProName = Form("p_mChargedV3Epd_Cent%d",i_cent);
    p_mChargedV3Epd[i_cent] = new TProfile(ProName.c_str(),ProName.c_str(),recoEP::mNumOfRunIndex,-0.5,(float)recoEP::mNumOfRunIndex-0.5);
    
    ProName = Form("p_mChargedV1EpEta_Cent%d",i_cent);
    p_mChargedV1EpEta[i_cent] = new TProfile(ProName.c_str(),ProName.c_str(),30,-5.5,5.5);
    ProName = Form("p_mChargedV3EpEta_Cent%d",i_cent);
    p_mChargedV3EpEta[i_cent] = new TProfile(ProName.c_str(),ProName.c_str(),30,-5.5,5.5);
    ProName = Form("p_mChargedV3EpEta_Cent%d",i_cent);
    p_mChargedV3EpEta[i_cent] = new TProfile(ProName.c_str(),ProName.c_str(),30,-5.5,5.5);
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

void StEventPlaneProManager::fillChargedV1Ep(double pt, double eta, double v1, double res1, int Cent9, int RunIndex, double reweight)
{
  if(res1 > 0.0)
  {
    if(pt >= 0.2 && pt <= 5.2)
    {
      p_mChargedV1EpQA[Cent9]->Fill(RunIndex,v1/res1,reweight);
      p_mChargedV1EpQA[9]->Fill(RunIndex,v1/res1,reweight);
    }

    p_mChargedV1Ep[Cent9]->Fill(pt,v1/res1,reweight);
    p_mChargedV1Ep[9]->Fill(pt,v1/res1,reweight);

    p_mChargedV1EpEta[Cent9]->Fill(eta,v1/res1,reweight);
    p_mChargedV1EpEta[9]->Fill(eta,v1/res1,reweight);
  }
}

void StEventPlaneProManager::fillChargedV2Ep(double pt, double eta, double v2, double res2, int Cent9, int RunIndex, double reweight)
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
    
    p_mChargedV2EpEta[Cent9]->Fill(eta,v2/res2,reweight);
    p_mChargedV2EpEta[9]->Fill(eta,v2/res2,reweight);
  }
}

void StEventPlaneProManager::fillChargedV3Ep(double pt, double eta, double v3, double res3, int Cent9, int RunIndex, double reweight)
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
    
    p_mChargedV3EpEta[Cent9]->Fill(eta,v3/res3,reweight);
    p_mChargedV3EpEta[9]->Fill(eta,v3/res3,reweight);
  }
}

void StEventPlaneProManager::fillChargedV1Epd(double pt, double v1, double res1, int Cent9, int RunIndex)
{
  if(res1 > 0.0)
  {
    if(pt >= 0.2 && pt <= 5.2)
    {
      p_mChargedV2EpdQA[Cent9]->Fill(RunIndex,v1/res1);
      p_mChargedV2EpdQA[9]->Fill(RunIndex,v1/res1);
    }

    p_mChargedV1Epd[Cent9]->Fill(pt,v1/res1);
    p_mChargedV1Epd[9]->Fill(pt,v1/res1);
  }
}

void StEventPlaneProManager::fillChargedV2Epd(double pt, double v2, double res2, int Cent9, int RunIndex)
{
  if(res2 > 0.0)
  {
    if(pt >= 0.2 && pt <= 5.2)
    {
      p_mChargedV2EpdQA[Cent9]->Fill(RunIndex,v2/res2);
      p_mChargedV2EpdQA[9]->Fill(RunIndex,v2/res2);
    }

    p_mChargedV2Epd[Cent9]->Fill(pt,v2/res2);
    p_mChargedV2Epd[9]->Fill(pt,v2/res2);
  }
}

void StEventPlaneProManager::fillChargedV3Epd(double pt, double v3, double res3, int Cent9, int RunIndex)
{
  if(res3 > 0.0)
  {
    if(pt >= 0.2 && pt <= 5.2)
    {
      p_mChargedV3EpdQA[Cent9]->Fill(RunIndex,v3/res3);
      p_mChargedV3EpdQA[9]->Fill(RunIndex,v3/res3);
    }

    p_mChargedV3Epd[Cent9]->Fill(pt,v3/res3);
    p_mChargedV3Epd[9]->Fill(pt,v3/res3);
  }
}


void StEventPlaneProManager::writeChargedFlow()
{
  for(int i_cent = 0; i_cent < 10; ++i_cent)
  {
    //std::cout << "we got to the charged flow writing" << std::endl;
    p_mChargedV1PpQA[i_cent]->Write();
    p_mChargedV1EpQA[i_cent]->Write();
    p_mChargedV1EpdQA[i_cent]->Write();
    p_mChargedV2EpQA[i_cent]->Write();
    p_mChargedV2PpQA[i_cent]->Write();
    p_mChargedV2EpdQA[i_cent]->Write();
    p_mChargedV3EpQA[i_cent]->Write();
    p_mChargedV3EpdQA[i_cent]->Write();
    p_mChargedV1Pp[i_cent]->Write();
    p_mChargedV1Ep[i_cent]->Write();
    p_mChargedV1Epd[i_cent]->Write();
    p_mChargedV2Ep[i_cent]->Write();
    p_mChargedV2Pp[i_cent]->Write();
    p_mChargedV2Epd[i_cent]->Write();
    p_mChargedV3Ep[i_cent]->Write();
    p_mChargedV3Epd[i_cent]->Write();
    p_mChargedV1EpEta[i_cent]->Write();
    p_mChargedV2EpEta[i_cent]->Write();
    p_mChargedV3EpEta[i_cent]->Write();
  }
}

//---------------------------------------------------------------------------------

void StEventPlaneProManager::initEpdRes()
{
  std::string corrName[3] = {"Raw", "Pw", "PwS"};
  std::string corrTitle[3] = {"Raw", "#phi-weighted", "#phi-weighted & Shifted"}; 

  for (int order=1; order<=recoEP::mEpdEpOrder; order++)
  {
    for (int  icor=0; icor<3; icor++)
    {
      p_mEpdAveCos[order-1][icor] = new TProfile(Form("p_mEpdCos%d%s",order,corrName[icor].c_str()),  Form("EPD %s #LT cos(%d#Psi_{%d,W}-%d#Psi_{%d,E}) #GT",corrTitle[icor].c_str(),order,order,order,order)  ,9,-0.5,8.5);
    }
 
    p_mEpdSubResQA[order-1] = new TProfile2D(Form("p_mEpdSubRes%dQA",order),Form("p_mEpdSubRes%dQA",order),recoEP::mNumOfRunIndex,-0.5,(float)recoEP::mNumOfRunIndex-0.5,9,-0.5,8.5);
    p_mEpdSubRes[order-1] = new TProfile(Form("p_mEpdSubRes%d",order),Form("p_mEpdSubRes%d",order),9,-0.5,8.5);
  }
}

void StEventPlaneProManager::fillEpdRes(StEpdEpInfo result, int Cent9, int RunIndex)
{
  for (int order=1; order<=recoEP::mEpdEpOrder; order++)
  { 
    p_mEpdAveCos[order-1][0]->Fill(double(Cent9),cos((double)order*(result.EastRawPsi(order)-result.WestRawPsi(order))));
    p_mEpdAveCos[order-1][1]->Fill(double(Cent9),cos((double)order*(result.EastPhiWeightedPsi(order)-result.WestPhiWeightedPsi(order))));
    p_mEpdAveCos[order-1][2]->Fill(double(Cent9),cos((double)order*(result.EastPhiWeightedAndShiftedPsi(order)-result.WestPhiWeightedAndShiftedPsi(order)))); 

    p_mEpdSubResQA[order-1]->Fill((double)RunIndex,(double)Cent9,cos((double)order*(result.EastPhiWeightedAndShiftedPsi(order)-result.WestPhiWeightedAndShiftedPsi(order))));
    p_mEpdSubRes[order-1]->Fill((double)Cent9,cos((double)order*(result.EastPhiWeightedAndShiftedPsi(order)-result.WestPhiWeightedAndShiftedPsi(order))));
  }
}

void StEventPlaneProManager::writeEpdRes()
{
  for (int order=1; order<=recoEP::mEpdEpOrder; order++)
  {
    for (int  icor=0; icor<3; icor++)
    {
      p_mEpdAveCos[order-1][icor]->Write();
    }
    p_mEpdSubResQA[order-1]->Write();
    p_mEpdSubRes[order-1]->Write();
  }
}

void StEventPlaneProManager::initEpdFlowEta()
{
  for (int order=1; order<=recoEP::mEpdEpOrder; order++)
  { 
    for (int icent=0; icent < 9; ++icent)
    { 
      p_mEpdFlowEta[order-1][icent] = new TProfile(Form("p_mEpdFlowEta%d_Cent%d",order,icent),Form("p_mEpdFlowEta%d_Cent%d",order,icent),30,-5.5,5.5);
      p_mEpdFlowEtaWeights[order-1][icent] = new TProfile(Form("p_mEpdFlowEtaWeights%d_Cent%d",order,icent),Form("p_mEpdFlowEtaWeights%d_Cent%d",order,icent),30,-5.5,5.5);
    }
  }
}

void StEventPlaneProManager::fillEpdFlowEta(double eta, double v, int Cent9, int order, double weight)
{
  p_mEpdFlowEta[order-1][Cent9]->Fill(eta,v,weight);
  p_mEpdFlowEtaWeights[order-1][Cent9]->Fill(eta,weight);
}

void StEventPlaneProManager::writeEpdFlowEta()
{
  for (int order=1; order<=recoEP::mEpdEpOrder; order++)
  { 
    for (int icent=0; icent < 9; ++icent)
    {
      p_mEpdFlowEta[order-1][icent]->Write();
      p_mEpdFlowEtaWeights[order-1][icent]->Write();
    }
  }
}

void StEventPlaneProManager::initEpdFlow()
{
  for (int order=1; order<=recoEP::mEpdEpOrder; order++)
  { 
    p_mEpdFlow[order-1] = new TProfile(Form("p_mEpdFlow%d",order),Form("p_mEpdFlow%d",order),9,-0.5,8.5);
    for (int icent=0; icent<10; ++icent)
    {
      p_mEpdFlowQA[order-1][icent] = new TProfile(Form("p_mEpdFlowQA%d",order),Form("p_mEpdFlowQA%d",order),recoEP::mNumOfRunIndex,-0.5,(float)recoEP::mNumOfRunIndex-0.5);   
    }
  }
}

void StEventPlaneProManager::fillEpdFlow(double v, int Cent9, int order, int runIndex)
{
  p_mEpdFlow[order-1]->Fill(Cent9,v);
  p_mEpdFlowQA[order-1][Cent9]->Fill(runIndex,v);
  p_mEpdFlowQA[order-1][9]->Fill(runIndex,v);
}

void StEventPlaneProManager::writeEpdFlow()
{
  for (int order=1; order<=recoEP::mEpdEpOrder; order++)
  { 
    p_mEpdFlow[order-1]->Write();
    for (int icent=0; icent<10; ++icent)
    {
      p_mEpdFlowQA[order-1][icent]->Write();
    }
  }
}
