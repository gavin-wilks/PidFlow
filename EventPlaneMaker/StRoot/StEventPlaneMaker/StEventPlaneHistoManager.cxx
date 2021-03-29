#include <TH2F.h>
#include <TH1F.h>
#include <TMath.h>
#include <TString.h>

#include "StRoot/StEventPlaneMaker/StEventPlaneHistoManager.h"
#include "StRoot/StEventPlaneMaker/StEventPlaneCons.h"

ClassImp(StEventPlaneHistoManager)

//-------------------------------------------------------------------------------------------

StEventPlaneHistoManager::StEventPlaneHistoManager()
{
}

//-------------------------------------------------------------------------------------------

StEventPlaneHistoManager::~StEventPlaneHistoManager()
{
  /* */
}

//-------------------------------------------------------------------------------------------
// ZDC EP
void StEventPlaneHistoManager::initZdcGainCorr()
{
  for(int i_eastwest = 0; i_eastwest < 2; ++i_eastwest)
  {
    for(int i_verthori = 0; i_verthori < 2; ++i_verthori)
    {
      for(int i_slat = 0; i_slat < 8; ++i_slat)
      {
	string HistName = Form("h_mZdcGainCorr%s%s_%d",recoEP::mEastWest[i_eastwest].c_str(),recoEP::mVertHori[i_verthori].c_str(),i_slat);
	h_mZdcGainCorr[i_eastwest][i_verthori][i_slat] = new TH2F(HistName.c_str(),HistName.c_str(),recoEP::mNumOfRunIndex,-0.5,(double)recoEP::mNumOfRunIndex-0.5,5000,-4.5,4995.5);
      }
    }
  }
}

void StEventPlaneHistoManager::fillZdcGainCorr(int i_eastwest, int i_verthori, int i_slat, int runIndex, float zdcsmd)
{
  h_mZdcGainCorr[i_eastwest][i_verthori][i_slat]->Fill((double)runIndex,zdcsmd);
}

void StEventPlaneHistoManager::writeZdcGainCorr()
{
  for(int i_eastwest = 0; i_eastwest < 2; ++i_eastwest)
  {
    for(int i_verthori = 0; i_verthori < 2; ++i_verthori)
    {
      for(int i_slat = 0; i_slat < 8; ++i_slat)
      {
	h_mZdcGainCorr[i_eastwest][i_verthori][i_slat]->Write();
      }
    }
  }
}

// raw ZDC-SMD EP
void StEventPlaneHistoManager::initZdcRawEP()
{
  for(int i_cent = 0; i_cent < 9; ++i_cent)
  {
    string HistName = Form("h_mZdcRawEpEast_%d",i_cent);
    h_mZdcRawEpEast[i_cent] = new TH2F(HistName.c_str(),HistName.c_str(),recoEP::mNumOfRunIndex,-0.5,(double)recoEP::mNumOfRunIndex-0.5,360,-1.0*TMath::Pi(),TMath::Pi());
    HistName = Form("h_mZdcRawEpWest_%d",i_cent);
    h_mZdcRawEpWest[i_cent] = new TH2F(HistName.c_str(),HistName.c_str(),recoEP::mNumOfRunIndex,-0.5,(double)recoEP::mNumOfRunIndex-0.5,360,-1.0*TMath::Pi(),TMath::Pi());
    HistName = Form("h_mZdcRawEpFull_%d",i_cent);
    h_mZdcRawEpFull[i_cent] = new TH2F(HistName.c_str(),HistName.c_str(),recoEP::mNumOfRunIndex,-0.5,(double)recoEP::mNumOfRunIndex-0.5,360,-1.0*TMath::Pi(),TMath::Pi());
  }
}

void StEventPlaneHistoManager::fillZdcRawSubEP(TVector2 QEast, TVector2 QWest, int Cent9, int runIndex)
{
  double PsiEast = TMath::ATan2(QEast.Y(),QEast.X()); h_mZdcRawEpEast[Cent9]->Fill(runIndex,PsiEast);
  double PsiWest = TMath::ATan2(QWest.Y(),QWest.X()); h_mZdcRawEpWest[Cent9]->Fill(runIndex,PsiWest);
}

void StEventPlaneHistoManager::fillZdcRawFullEP(TVector2 QFull, int Cent9, int runIndex)
{
  double PsiFull = TMath::ATan2(QFull.Y(),QFull.X()); h_mZdcRawEpFull[Cent9]->Fill(runIndex,PsiFull);
}

void StEventPlaneHistoManager::writeZdcRawEP()
{
  for(int i_cent = 0; i_cent < 9; ++i_cent)
  {
    h_mZdcRawEpEast[i_cent]->Write();
    h_mZdcRawEpWest[i_cent]->Write();
    h_mZdcRawEpFull[i_cent]->Write();
  }
}

// recenter ZDC-SMD EP
void StEventPlaneHistoManager::initZdcReCenterEP()
{
  for(int i_cent = 0; i_cent < 9; ++i_cent)
  {
    string HistName = Form("h_mZdcReCenterEpEast_%d",i_cent);
    h_mZdcReCenterEpEast[i_cent] = new TH2F(HistName.c_str(),HistName.c_str(),recoEP::mNumOfRunIndex,-0.5,(double)recoEP::mNumOfRunIndex-0.5,360,-1.0*TMath::Pi(),TMath::Pi());
    HistName = Form("h_mZdcReCenterEpWest_%d",i_cent);
    h_mZdcReCenterEpWest[i_cent] = new TH2F(HistName.c_str(),HistName.c_str(),recoEP::mNumOfRunIndex,-0.5,(double)recoEP::mNumOfRunIndex-0.5,360,-1.0*TMath::Pi(),TMath::Pi());
    HistName = Form("h_mZdcReCenterEpFull_%d",i_cent);
    h_mZdcReCenterEpFull[i_cent] = new TH2F(HistName.c_str(),HistName.c_str(),recoEP::mNumOfRunIndex,-0.5,(double)recoEP::mNumOfRunIndex-0.5,360,-1.0*TMath::Pi(),TMath::Pi());
  }
}

void StEventPlaneHistoManager::fillZdcReCenterSubEP(TVector2 QEast, TVector2 QWest, int Cent9, int runIndex)
{
  double PsiEast = TMath::ATan2(QEast.Y(),QEast.X()); h_mZdcReCenterEpEast[Cent9]->Fill(runIndex,PsiEast);
  double PsiWest = TMath::ATan2(QWest.Y(),QWest.X()); h_mZdcReCenterEpWest[Cent9]->Fill(runIndex,PsiWest);
}

void StEventPlaneHistoManager::fillZdcReCenterFullEP(TVector2 QFull, int Cent9, int runIndex)
{
  double PsiFull = TMath::ATan2(QFull.Y(),QFull.X()); h_mZdcReCenterEpFull[Cent9]->Fill(runIndex,PsiFull);
}

void StEventPlaneHistoManager::writeZdcReCenterEP()
{
  for(int i_cent = 0; i_cent < 9; ++i_cent)
  {
    h_mZdcReCenterEpEast[i_cent]->Write();
    h_mZdcReCenterEpWest[i_cent]->Write();
    h_mZdcReCenterEpFull[i_cent]->Write();
  }
}
//-------------------------------------------------------------------------------------------

//-------------------------------------------------------------------------------------------
// TPC EP
// raw TPC EP
void StEventPlaneHistoManager::initTpcRawEP()
{
  for(int i_cent = 0; i_cent < 9; ++i_cent)
  {
    string HistName = Form("h_mTpcRawEpEast_%d",i_cent);
    h_mTpcRawEpEast[i_cent] = new TH2F(HistName.c_str(),HistName.c_str(),recoEP::mNumOfRunIndex,-0.5,(double)recoEP::mNumOfRunIndex-0.5,360,-1.0*TMath::Pi(),TMath::Pi());
    HistName = Form("h_mTpcRawEpWest_%d",i_cent);
    h_mTpcRawEpWest[i_cent] = new TH2F(HistName.c_str(),HistName.c_str(),recoEP::mNumOfRunIndex,-0.5,(double)recoEP::mNumOfRunIndex-0.5,360,-1.0*TMath::Pi(),TMath::Pi());
    HistName = Form("h_mTpcRawEpFull_%d",i_cent);
    h_mTpcRawEpFull[i_cent] = new TH2F(HistName.c_str(),HistName.c_str(),recoEP::mNumOfRunIndex,-0.5,(double)recoEP::mNumOfRunIndex-0.5,360,-1.0*TMath::Pi(),TMath::Pi());
  }
}

void StEventPlaneHistoManager::fillTpcRawSubEP(TVector2 QEast, TVector2 QWest, int Cent9, int runIndex)
{
  double PsiEast = 0.5*TMath::ATan2(QEast.Y(),QEast.X()); h_mTpcRawEpEast[Cent9]->Fill(runIndex,PsiEast);
  double PsiWest = 0.5*TMath::ATan2(QWest.Y(),QWest.X()); h_mTpcRawEpWest[Cent9]->Fill(runIndex,PsiWest);
}

void StEventPlaneHistoManager::fillTpcRawFullEP(TVector2 QFull, int Cent9, int runIndex)
{
  double PsiFull = 0.5*TMath::ATan2(QFull.Y(),QFull.X()); h_mTpcRawEpFull[Cent9]->Fill(runIndex,PsiFull);
}

void StEventPlaneHistoManager::writeTpcRawEP()
{
  for(int i_cent = 0; i_cent < 9; ++i_cent)
  {
    h_mTpcRawEpEast[i_cent]->Write();
    h_mTpcRawEpWest[i_cent]->Write();
    h_mTpcRawEpFull[i_cent]->Write();
  }
}

// recenter TPC EP
void StEventPlaneHistoManager::initTpcReCenterEP()
{
  for(int i_cent = 0; i_cent < 9; ++i_cent)
  {
    string HistName = Form("h_mTpcReCenterEpEast_%d",i_cent);
    h_mTpcReCenterEpEast[i_cent] = new TH2F(HistName.c_str(),HistName.c_str(),recoEP::mNumOfRunIndex,-0.5,(double)recoEP::mNumOfRunIndex-0.5,360,-1.0*TMath::Pi(),TMath::Pi());
    HistName = Form("h_mTpcReCenterEpWest_%d",i_cent);
    h_mTpcReCenterEpWest[i_cent] = new TH2F(HistName.c_str(),HistName.c_str(),recoEP::mNumOfRunIndex,-0.5,(double)recoEP::mNumOfRunIndex-0.5,360,-1.0*TMath::Pi(),TMath::Pi());
    HistName = Form("h_mTpcReCenterEpFull_%d",i_cent);
    h_mTpcReCenterEpFull[i_cent] = new TH2F(HistName.c_str(),HistName.c_str(),recoEP::mNumOfRunIndex,-0.5,(double)recoEP::mNumOfRunIndex-0.5,360,-1.0*TMath::Pi(),TMath::Pi());
  }
}

void StEventPlaneHistoManager::fillTpcReCenterSubEP(TVector2 QEast, TVector2 QWest, int Cent9, int runIndex)
{
  double PsiEast = 0.5*TMath::ATan2(QEast.Y(),QEast.X()); h_mTpcReCenterEpEast[Cent9]->Fill(runIndex,PsiEast);
  double PsiWest = 0.5*TMath::ATan2(QWest.Y(),QWest.X()); h_mTpcReCenterEpWest[Cent9]->Fill(runIndex,PsiWest);
}

void StEventPlaneHistoManager::fillTpcReCenterFullEP(TVector2 QFull, int Cent9, int runIndex)
{
  double PsiFull = 0.5*TMath::ATan2(QFull.Y(),QFull.X()); h_mTpcReCenterEpFull[Cent9]->Fill(runIndex,PsiFull);
}

void StEventPlaneHistoManager::writeTpcReCenterEP()
{
  for(int i_cent = 0; i_cent < 9; ++i_cent)
  {
    h_mTpcReCenterEpEast[i_cent]->Write();
    h_mTpcReCenterEpWest[i_cent]->Write();
    h_mTpcReCenterEpFull[i_cent]->Write();
  }
}
//-------------------------------------------------------------------------------------------
