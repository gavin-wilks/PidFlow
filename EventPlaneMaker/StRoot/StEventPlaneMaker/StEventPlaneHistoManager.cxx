#include <TH2F.h>
#include <TH1F.h>
#include <TProfile.h>
#include <TMath.h>
#include <TString.h>

#include "StRoot/StEventPlaneMaker/StEventPlaneHistoManager.h"
#include "StRoot/StEventPlaneMaker/StEventPlaneCons.h"
#include "StEpdUtil/StEpdEpInfo.h"

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

// shift ZDC-SMD EP
void StEventPlaneHistoManager::initZdcShiftEP()
{
  for(int i_cent = 0; i_cent < 9; ++i_cent)
  {
    string HistName = Form("h_mZdcShiftEpEast_%d",i_cent);
    h_mZdcShiftEpEast[i_cent] = new TH2F(HistName.c_str(),HistName.c_str(),recoEP::mNumOfRunIndex,-0.5,(double)recoEP::mNumOfRunIndex-0.5,360,-1.0*TMath::Pi(),TMath::Pi());
    HistName = Form("h_mZdcShiftEpWest_%d",i_cent);
    h_mZdcShiftEpWest[i_cent] = new TH2F(HistName.c_str(),HistName.c_str(),recoEP::mNumOfRunIndex,-0.5,(double)recoEP::mNumOfRunIndex-0.5,360,-1.0*TMath::Pi(),TMath::Pi());
    HistName = Form("h_mZdcShiftEpDiff_%d",i_cent);
    h_mZdcShiftEpDiff[i_cent] = new TH2F(HistName.c_str(),HistName.c_str(),recoEP::mNumOfRunIndex,-0.5,(double)recoEP::mNumOfRunIndex-0.5,360,-1.0*TMath::Pi(),TMath::Pi());
    HistName = Form("h_mZdcShiftEpFull_%d",i_cent);
    h_mZdcShiftEpFull[i_cent] = new TH2F(HistName.c_str(),HistName.c_str(),recoEP::mNumOfRunIndex,-0.5,(double)recoEP::mNumOfRunIndex-0.5,360,-1.0*TMath::Pi(),TMath::Pi());
  }
}

void StEventPlaneHistoManager::fillZdcShiftSubEP(TVector2 QEast, TVector2 QWest, int Cent9, int runIndex)
{
  double PsiEast = TMath::ATan2(QEast.Y(),QEast.X()); h_mZdcShiftEpEast[Cent9]->Fill(runIndex,PsiEast);
  double PsiWest = TMath::ATan2(QWest.Y(),QWest.X()); h_mZdcShiftEpWest[Cent9]->Fill(runIndex,PsiWest);
}

void StEventPlaneHistoManager::fillZdcShiftFullEP(TVector2 QDiff, TVector2 QFull, int Cent9, int runIndex)
{
  double PsiDiff = TMath::ATan2(QDiff.Y(),QDiff.X()); h_mZdcShiftEpDiff[Cent9]->Fill(runIndex,PsiDiff);
  double PsiFull = TMath::ATan2(QFull.Y(),QFull.X()); h_mZdcShiftEpFull[Cent9]->Fill(runIndex,PsiFull);
}

void StEventPlaneHistoManager::writeZdcShiftEP()
{
  for(int i_cent = 0; i_cent < 9; ++i_cent)
  {
    h_mZdcShiftEpEast[i_cent]->Write();
    h_mZdcShiftEpWest[i_cent]->Write();
    h_mZdcShiftEpDiff[i_cent]->Write();
    h_mZdcShiftEpFull[i_cent]->Write();
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

// shift TPC EP
void StEventPlaneHistoManager::initTpcShiftEP()
{
  for(int i_cent = 0; i_cent < 9; ++i_cent)
  {
    string HistName = Form("h_mTpcShiftEpEast_%d",i_cent);
    h_mTpcShiftEpEast[i_cent] = new TH2F(HistName.c_str(),HistName.c_str(),recoEP::mNumOfRunIndex,-0.5,(double)recoEP::mNumOfRunIndex-0.5,360,-1.0*TMath::Pi(),TMath::Pi());
    HistName = Form("h_mTpcShiftEpWest_%d",i_cent);
    h_mTpcShiftEpWest[i_cent] = new TH2F(HistName.c_str(),HistName.c_str(),recoEP::mNumOfRunIndex,-0.5,(double)recoEP::mNumOfRunIndex-0.5,360,-1.0*TMath::Pi(),TMath::Pi());
    HistName = Form("h_mTpcShiftEpRanA_%d",i_cent);
    h_mTpcShiftEpRanA[i_cent] = new TH2F(HistName.c_str(),HistName.c_str(),recoEP::mNumOfRunIndex,-0.5,(double)recoEP::mNumOfRunIndex-0.5,360,-1.0*TMath::Pi(),TMath::Pi());
    HistName = Form("h_mTpcShiftEpRanB_%d",i_cent);
    h_mTpcShiftEpRanB[i_cent] = new TH2F(HistName.c_str(),HistName.c_str(),recoEP::mNumOfRunIndex,-0.5,(double)recoEP::mNumOfRunIndex-0.5,360,-1.0*TMath::Pi(),TMath::Pi());
    HistName = Form("h_mTpcShiftEpFull_%d",i_cent);
    h_mTpcShiftEpFull[i_cent] = new TH2F(HistName.c_str(),HistName.c_str(),recoEP::mNumOfRunIndex,-0.5,(double)recoEP::mNumOfRunIndex-0.5,360,-1.0*TMath::Pi(),TMath::Pi());
  }
}

void StEventPlaneHistoManager::fillTpcShiftSubEP(double PsiEast, double PsiWest, int Cent9, int runIndex)
{
  h_mTpcShiftEpEast[Cent9]->Fill(runIndex,PsiEast);
  h_mTpcShiftEpWest[Cent9]->Fill(runIndex,PsiWest);
}

void StEventPlaneHistoManager::fillTpcShiftRanEP(double PsiRanA, double PsiRanB, int Cent9, int runIndex)
{
  h_mTpcShiftEpRanA[Cent9]->Fill(runIndex,PsiRanA);
  h_mTpcShiftEpRanB[Cent9]->Fill(runIndex,PsiRanB);
}

void StEventPlaneHistoManager::fillTpcShiftFullEP(double PsiFull, int Cent9, int runIndex)
{
  h_mTpcShiftEpFull[Cent9]->Fill(runIndex,PsiFull);
}

void StEventPlaneHistoManager::writeTpcShiftEP()
{
  for(int i_cent = 0; i_cent < 9; ++i_cent)
  {
    h_mTpcShiftEpEast[i_cent]->Write();
    h_mTpcShiftEpWest[i_cent]->Write();
    h_mTpcShiftEpRanA[i_cent]->Write();
    h_mTpcShiftEpRanB[i_cent]->Write();
    h_mTpcShiftEpFull[i_cent]->Write();
  }
}

void StEventPlaneHistoManager::initEpdEpResults()
{
  // histograms and profiles of event planes
  // EPD "correction level" (last index):  0=raw; 1=shifted
  for (int order=1; order<=recoEP::mEpdEpOrder; order++)
  {
    double PhiMax = 2.0*TMath::Pi()/((double)order);
    h_mEpdEwPsi[order-1][0]  = new TH2F(Form("EpdEwPsi%dRaw",order),  Form("EPD Raw #Psi_{%d,W} vs #Psi_{%d,E}",order,order),  100,0.0,PhiMax,100,0.0,PhiMax);
    h_mEpdEwPsi[order-1][1]  = new TH2F(Form("EpdEwPsi%dPw",order),Form("EPD #phi-weighted #Psi_{%d,W} vs #Psi_{%d,E}",order,order),100,0.0,PhiMax,100,0.0,PhiMax);
    h_mEpdEwPsi[order-1][2]  = new TH2F(Form("EpdEwPsi%dPwS",order),Form("EPD #phi-weighted & Shifted #Psi_{%d,W} vs #Psi_{%d,E}",order,order),100,0.0,PhiMax,100,0.0,PhiMax);    

    //h_mEpdEwPsi_midCentral[order-1][0]  = new TH2F(Form("EpdEwPsi%drawMidCent",order),  Form("MidCentral EPD raw #Psi_{%d,W} vs #Psi_{%d,E}",order,order),  100,0.0,PhiMax,100,0.0,PhiMax);
    //h_mEpdEwPsi_midCentral[order-1][1]  = new TH2F(Form("EpdEwPsi%dshiftMidCent",order),Form("MidCentral EPD shift #Psi_{%d,W} vs #Psi_{%d,E}",order,order),100,0.0,PhiMax,100,0.0,PhiMax);
    for (int icor=0; icor<3; icor++){
      h_mEpdEwPsi[order-1][icor]->GetXaxis()->SetTitle(Form("#Psi_{%d,E}",order));                 h_mEpdEwPsi[order-1][icor]->GetYaxis()->SetTitle(Form("#Psi_{%d,W}",order));
      //h_mEpdEwPsi_midCentral[order-1][icor]->GetXaxis()->SetTitle(Form("#Psi_{%d,E}",order));      h_mEpdEwPsi_midCentral[order-1][icor]->GetYaxis()->SetTitle(Form("#Psi_{%d,W}",order));
    } 
    //
    h_mEpdAveCos[order-1][0] = new TProfile(Form("EpdCos%dRaw",order),  Form("EPD Raw #LT cos(%d#Psi_{%d,W}-%d#Psi_{%d,E}) #GT",order,order,order,order)  ,9,-0.5,8.5);
    h_mEpdAveCos[order-1][1] = new TProfile(Form("EpdCos%dPw",order),Form("EPD #phi-weighted #LT cos(%d#Psi_{%d,W}-%d#Psi_{%d,E}) #GT",order,order,order,order),9,-0.5,8.5);
    h_mEpdAveCos[order-1][2] = new TProfile(Form("EpdCos%dPwS",order),Form("EPD #phi-weighted & Shifted #LT cos(%d#Psi_{%d,W}-%d#Psi_{%d,E}) #GT",order,order,order,order),9,-0.5,8.5);

  }
  
  h_mEpdAveCosD12[0] = new TProfile("EpdD12CosRaw","EPD D12 Raw #LT cos(2#Psi_{1^{st}}-2#Psi_{2^{nd}}) #GT",9,-0.5,8.5); 
  h_mEpdAveCosD12[1] = new TProfile("EpdD12CosPw" ,"EPD D12 #phi-weighted #LT cos(2#Psi_{1^{st}}-2#Psi_{2^{nd}}) #GT",9,-0.5,8.5);
  h_mEpdAveCosD12[2] = new TProfile("EpdD12CosPwS","EPD D12 #phi-weighted & Shifted #LT cos(2#Psi_{1^{st}}-2#Psi_{2^{nd}}) #GT",9,-0.5,8.5);

  h_mEpdFullPsi21[0] = new TH2F("EpdFullPsi21Raw","EPD Raw #Psi_{2} vs #Psi_{1}",100,0.0,2.0*TMath::Pi(),100,0.0,2.0*TMath::Pi());
  h_mEpdFullPsi21[1] = new TH2F("EpdFullPsi21Pw" ,"EPD #phi-weighted #Psi_{2} vs #Psi_{1}",100,0.0,2.0*TMath::Pi(),100,0.0,2.0*TMath::Pi());
  h_mEpdFullPsi21[2] = new TH2F("EpdFullPsi21PwS","EPD #phi-weighted & Shifted #Psi_{2} vs #Psi_{1}",100,0.0,2.0*TMath::Pi(),100,0.0,2.0*TMath::Pi()); 

  h_mEpdFullPsi31[0] = new TH2F("EpdFullPsi31Raw","EPD Raw #Psi_{3} vs #Psi_{1}",100,0.0,2.0*TMath::Pi(),100,0.0,2.0*TMath::Pi());
  h_mEpdFullPsi31[1] = new TH2F("EpdFullPsi31Pw" ,"EPD #phi-weighted #Psi_{3} vs #Psi_{1}",100,0.0,2.0*TMath::Pi(),100,0.0,2.0*TMath::Pi());
  h_mEpdFullPsi31[2] = new TH2F("EpdFullPsi31PwS","EPD #phi-weighted & Shifted #Psi_{3} vs #Psi_{1}",100,0.0,2.0*TMath::Pi(),100,0.0,2.0*TMath::Pi()); 

  h_mEpdFullPsi32[0] = new TH2F("EpdFullPsi32Raw","EPD Raw #Psi_{3} vs #Psi_{2}",100,0.0,2.0*TMath::Pi(),100,0.0,2.0*TMath::Pi());
  h_mEpdFullPsi32[1] = new TH2F("EpdFullPsi32Pw" ,"EPD #phi-weighted #Psi_{3} vs #Psi_{2}",100,0.0,2.0*TMath::Pi(),100,0.0,2.0*TMath::Pi());
  h_mEpdFullPsi32[2] = new TH2F("EpdFullPsi32PwS","EPD #phi-weighted & Shifted #Psi_{3} vs #Psi_{2}",100,0.0,2.0*TMath::Pi(),100,0.0,2.0*TMath::Pi()); 
}
 
void StEventPlaneHistoManager::fillEpdEpResults(StEpdEpInfo result, int CentId)
{
  for (int order=1; order<=recoEP::mEpdEpOrder; order++){
    h_mEpdAveCos[order-1][0]->Fill(double(CentId),cos((double)order*(result.EastRawPsi(order)-result.WestRawPsi(order))));
    h_mEpdAveCos[order-1][1]->Fill(double(CentId),cos((double)order*(result.EastPhiWeightedPsi(order)-result.WestPhiWeightedPsi(order))));
    h_mEpdAveCos[order-1][2]->Fill(double(CentId),cos((double)order*(result.EastPhiWeightedAndShiftedPsi(order)-result.WestPhiWeightedAndShiftedPsi(order))));
    h_mEpdEwPsi[order-1][0]->Fill(result.EastRawPsi(order),result.WestRawPsi(order));
    h_mEpdEwPsi[order-1][1]->Fill(result.EastPhiWeightedPsi(order),result.WestPhiWeightedPsi(order));
    h_mEpdEwPsi[order-1][2]->Fill(result.EastPhiWeightedAndShiftedPsi(order),result.WestPhiWeightedAndShiftedPsi(order));
    //if ((CentId<=6)&&(CentId>=4)){
    //  h_mEpdEwPsi_midCentral[order-1][0]->Fill(result.EastRawPsi(order),result.WestRawPsi(order));
    //  h_mEpdEwPsi_midCentral[order-1][1]->Fill(result.EastPhiWeightedAndShiftedPsi(order),result.WestPhiWeightedAndShiftedPsi(order));
    //}
  }
  h_mEpdAveCosD12[0]->Fill(double(CentId),cos((double)2.0*(result.FullRawPsi(1)-result.FullRawPsi(2))));
  h_mEpdAveCosD12[1]->Fill(double(CentId),cos((double)2.0*(result.FullPhiWeightedPsi(1)-result.FullPhiWeightedPsi(2))));
  h_mEpdAveCosD12[2]->Fill(double(CentId),cos((double)2.0*(result.FullPhiWeightedAndShiftedPsi(1)-result.FullPhiWeightedAndShiftedPsi(2))));

  h_mEpdFullPsi21[0]->Fill(result.FullRawPsi(1),result.FullRawPsi(2));
  h_mEpdFullPsi21[1]->Fill(result.FullPhiWeightedPsi(1),result.FullPhiWeightedPsi(2));
  h_mEpdFullPsi21[2]->Fill(result.FullPhiWeightedAndShiftedPsi(1),result.FullPhiWeightedAndShiftedPsi(2));

  h_mEpdFullPsi31[0]->Fill(result.FullRawPsi(1),result.FullRawPsi(3));
  h_mEpdFullPsi31[1]->Fill(result.FullPhiWeightedPsi(1),result.FullPhiWeightedPsi(3));
  h_mEpdFullPsi31[2]->Fill(result.FullPhiWeightedAndShiftedPsi(1),result.FullPhiWeightedAndShiftedPsi(3));

  h_mEpdFullPsi32[0]->Fill(result.FullRawPsi(2),result.FullRawPsi(3));
  h_mEpdFullPsi32[1]->Fill(result.FullPhiWeightedPsi(2),result.FullPhiWeightedPsi(3));
  h_mEpdFullPsi32[2]->Fill(result.FullPhiWeightedAndShiftedPsi(2),result.FullPhiWeightedAndShiftedPsi(3));
}

void StEventPlaneHistoManager::writeEpdEpResults()
{
  for (int icor=0; icor<3; icor++)
  {
    for (int order=1; order<=recoEP::mEpdEpOrder; order++)
    { 
      h_mEpdEwPsi[order-1][icor]->Write();
      //h_mEpdEwPsi_midCentral[order-1][icor]->Write();
      h_mEpdAveCos[order-1][icor]->Write();
    }
    h_mEpdAveCosD12[icor]->Write();
    h_mEpdFullPsi21[icor]->Write();
    h_mEpdFullPsi31[icor]->Write();
    h_mEpdFullPsi32[icor]->Write();
  } 
}

//-------------------------------------------------------------------------------------------
