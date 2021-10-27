#include "StRoot/StRunQAMaker/StRunQAProManager.h"
#include "StRoot/StRunQAMaker/StRunQACons.h"
#include <TProfile.h>
#include <TMath.h>
#include <TString.h>

#include <string>

ClassImp(StRunQAProManager)

//---------------------------------------------------------------------------------

StRunQAProManager::StRunQAProManager()
{
}

//---------------------------------------------------------------------------------

StRunQAProManager::~StRunQAProManager()
{
  /* */
}

//---------------------------------------------------------------------------------

void StRunQAProManager::initRunQA()
{
  for(int i_cut = 0; i_cut < 2; ++i_cut)
  {
    for(int i_trig = 0; i_trig < 10; ++i_trig)
    {
      std::string ProName = Form("p_mRefMult%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      p_mRefMult[i_cut][i_trig] = new TProfile(ProName.c_str(),ProName.c_str(),runQA::mNumOfRunIndex,-0.5,(float)runQA::mNumOfRunIndex-0.5); 

      ProName = Form("p_mGRefMult%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      p_mGRefMult[i_cut][i_trig] = new TProfile(ProName.c_str(),ProName.c_str(),runQA::mNumOfRunIndex,-0.5,(float)runQA::mNumOfRunIndex-0.5); 

      ProName = Form("p_mZdcX%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      p_mZdcX[i_cut][i_trig] = new TProfile(ProName.c_str(),ProName.c_str(),runQA::mNumOfRunIndex,-0.5,(float)runQA::mNumOfRunIndex-0.5);

      ProName = Form("p_mVz%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      p_mVz[i_cut][i_trig] = new TProfile(ProName.c_str(),ProName.c_str(),runQA::mNumOfRunIndex,-0.5,(float)runQA::mNumOfRunIndex-0.5);

      ProName = Form("p_mVr%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      p_mVr[i_cut][i_trig] = new TProfile(ProName.c_str(),ProName.c_str(),runQA::mNumOfRunIndex,-0.5,(float)runQA::mNumOfRunIndex-0.5);

      ProName = Form("p_mGDca%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      p_mGDca[i_cut][i_trig] = new TProfile(ProName.c_str(),ProName.c_str(),runQA::mNumOfRunIndex,-0.5,(float)runQA::mNumOfRunIndex-0.5);

      ProName = Form("p_mNHitsFit%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      p_mNHitsFit[i_cut][i_trig] = new TProfile(ProName.c_str(),ProName.c_str(),runQA::mNumOfRunIndex,-0.5,(float)runQA::mNumOfRunIndex-0.5);

      ProName = Form("p_mPrimPt%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      p_mPrimPt[i_cut][i_trig] = new TProfile(ProName.c_str(),ProName.c_str(),runQA::mNumOfRunIndex,-0.5,(float)runQA::mNumOfRunIndex-0.5);

      ProName = Form("p_mPrimEta%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      p_mPrimEta[i_cut][i_trig] = new TProfile(ProName.c_str(),ProName.c_str(),runQA::mNumOfRunIndex,-0.5,(float)runQA::mNumOfRunIndex-0.5);

      ProName = Form("p_mPrimPhi%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      p_mPrimPhi[i_cut][i_trig] = new TProfile(ProName.c_str(),ProName.c_str(),runQA::mNumOfRunIndex,-0.5,(float)runQA::mNumOfRunIndex-0.5);

      ProName = Form("p_mGlobPt%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      p_mGlobPt[i_cut][i_trig] = new TProfile(ProName.c_str(),ProName.c_str(),runQA::mNumOfRunIndex,-0.5,(float)runQA::mNumOfRunIndex-0.5);

      ProName = Form("p_mGlobEta%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      p_mGlobEta[i_cut][i_trig] = new TProfile(ProName.c_str(),ProName.c_str(),runQA::mNumOfRunIndex,-0.5,(float)runQA::mNumOfRunIndex-0.5);

      ProName = Form("p_mGlobPhi%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      p_mGlobPhi[i_cut][i_trig] = new TProfile(ProName.c_str(),ProName.c_str(),runQA::mNumOfRunIndex,-0.5,(float)runQA::mNumOfRunIndex-0.5);
  
      ProName = Form("p_mDEdx%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      p_mDEdx[i_cut][i_trig] = new TProfile(ProName.c_str(),ProName.c_str(),runQA::mNumOfRunIndex,-0.5,(float)runQA::mNumOfRunIndex-0.5);

      ProName = Form("p_mBetaInv%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      p_mBetaInv[i_cut][i_trig] = new TProfile(ProName.c_str(),ProName.c_str(),runQA::mNumOfRunIndex,-0.5,(float)runQA::mNumOfRunIndex-0.5);

      ProName = Form("p_mPrimaryMass2%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      p_mPrimaryMass2[i_cut][i_trig] = new TProfile(ProName.c_str(),ProName.c_str(),runQA::mNumOfRunIndex,-0.5,(float)runQA::mNumOfRunIndex-0.5);

      ProName = Form("p_mNHitsMax%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      p_mNHitsMax[i_cut][i_trig] = new TProfile(ProName.c_str(),ProName.c_str(),runQA::mNumOfRunIndex,-0.5,(float)runQA::mNumOfRunIndex-0.5);

      ProName = Form("p_mNHitsDEdx%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      p_mNHitsDEdx[i_cut][i_trig] = new TProfile(ProName.c_str(),ProName.c_str(),runQA::mNumOfRunIndex,-0.5,(float)runQA::mNumOfRunIndex-0.5);

      ProName = Form("p_mNHitsRatio%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      p_mNHitsRatio[i_cut][i_trig] = new TProfile(ProName.c_str(),ProName.c_str(),runQA::mNumOfRunIndex,-0.5,(float)runQA::mNumOfRunIndex-0.5);

      ProName = Form("p_mNSigmaPi%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      p_mNSigmaPi[i_cut][i_trig] = new TProfile(ProName.c_str(),ProName.c_str(),runQA::mNumOfRunIndex,-0.5,(float)runQA::mNumOfRunIndex-0.5);

      ProName = Form("p_mNSigmaK%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      p_mNSigmaK[i_cut][i_trig] = new TProfile(ProName.c_str(),ProName.c_str(),runQA::mNumOfRunIndex,-0.5,(float)runQA::mNumOfRunIndex-0.5);
      
      ProName = Form("p_mNSigmaP%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      p_mNSigmaP[i_cut][i_trig] = new TProfile(ProName.c_str(),ProName.c_str(),runQA::mNumOfRunIndex,-0.5,(float)runQA::mNumOfRunIndex-0.5);
 
      ProName = Form("p_mNSigmaE%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      p_mNSigmaE[i_cut][i_trig] = new TProfile(ProName.c_str(),ProName.c_str(),runQA::mNumOfRunIndex,-0.5,(float)runQA::mNumOfRunIndex-0.5);

      ProName = Form("p_mMass2Pi%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      p_mMass2Pi[i_cut][i_trig] = new TProfile(ProName.c_str(),ProName.c_str(),runQA::mNumOfRunIndex,-0.5,(float)runQA::mNumOfRunIndex-0.5);

      ProName = Form("p_mMass2K%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      p_mMass2K[i_cut][i_trig] = new TProfile(ProName.c_str(),ProName.c_str(),runQA::mNumOfRunIndex,-0.5,(float)runQA::mNumOfRunIndex-0.5);

      ProName = Form("p_mMass2P%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      p_mMass2P[i_cut][i_trig] = new TProfile(ProName.c_str(),ProName.c_str(),runQA::mNumOfRunIndex,-0.5,(float)runQA::mNumOfRunIndex-0.5);
 
      ProName = Form("p_mMass2E%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      p_mMass2E[i_cut][i_trig] = new TProfile(ProName.c_str(),ProName.c_str(),runQA::mNumOfRunIndex,-0.5,(float)runQA::mNumOfRunIndex-0.5);

    }
  }
}

void StRunQAProManager::fillRunQA_Event(int triggerBin, int runIdenx, float refMult, float grefMult, float zdcX, float vx, float vy, float vz, int cutSelection)
{
  // for a specific triggerBin
  p_mRefMult[cutSelection][triggerBin]->Fill(runIdenx, refMult);
  p_mGRefMult[cutSelection][triggerBin]->Fill(runIdenx, grefMult);
  p_mZdcX[cutSelection][triggerBin]->Fill(runIdenx, zdcX);
  p_mVz[cutSelection][triggerBin]->Fill(runIdenx, vz);
  p_mVr[cutSelection][triggerBin]->Fill(runIdenx, TMath::Sqrt(vx*vx+vy*vy));

  // for all triggers
  p_mRefMult[cutSelection][9]->Fill(runIdenx, refMult);
  p_mGRefMult[cutSelection][9]->Fill(runIdenx, grefMult);
  p_mZdcX[cutSelection][9]->Fill(runIdenx, zdcX);
  p_mVz[cutSelection][9]->Fill(runIdenx, vz);
  p_mVr[cutSelection][9]->Fill(runIdenx, TMath::Sqrt(vx*vx+vy*vy));
}

void StRunQAProManager::fillRunQA_Track(int triggerBin, int runIdenx, float gDca, int nHitsFit, TVector3 pMom, TVector3 gMom, float dEdx, float beta, float mass2, int nHitsMax, int nHitsDEdx, int cutSelection)
{
  // for a specific triggerBin
  p_mGDca[cutSelection][triggerBin]->Fill(runIdenx, gDca);
  p_mNHitsFit[cutSelection][triggerBin]->Fill(runIdenx, nHitsFit);
  p_mPrimPt[cutSelection][triggerBin]->Fill(runIdenx, pMom.Pt());
  p_mPrimEta[cutSelection][triggerBin]->Fill(runIdenx, pMom.Eta());
  p_mPrimPhi[cutSelection][triggerBin]->Fill(runIdenx, pMom.Phi());
  p_mGlobPt[cutSelection][triggerBin]->Fill(runIdenx, gMom.Pt());
  p_mGlobEta[cutSelection][triggerBin]->Fill(runIdenx, gMom.Eta());
  p_mGlobPhi[cutSelection][triggerBin]->Fill(runIdenx, gMom.Phi());
  p_mDEdx[cutSelection][triggerBin]->Fill(runIdenx, dEdx);
  if ( beta >= 0.0 ) p_mBetaInv[cutSelection][triggerBin]->Fill(runIdenx, 1.0/beta);
  if ( mass2 != -999.9 ) p_mPrimaryMass2[cutSelection][triggerBin]->Fill(runIdenx, mass2);
  p_mNHitsMax[cutSelection][triggerBin]->Fill(runIdenx, nHitsMax);
  p_mNHitsDEdx[cutSelection][triggerBin]->Fill(runIdenx, nHitsDEdx);
  p_mNHitsRatio[cutSelection][triggerBin]->Fill(runIdenx, (float)nHitsFit/(float)nHitsMax);

  // for all triggers
  p_mGDca[cutSelection][9]->Fill(runIdenx, gDca);
  p_mNHitsFit[cutSelection][9]->Fill(runIdenx, nHitsFit);
  p_mPrimPt[cutSelection][9]->Fill(runIdenx, pMom.Pt());
  p_mPrimEta[cutSelection][9]->Fill(runIdenx, pMom.Eta());
  p_mPrimPhi[cutSelection][9]->Fill(runIdenx, pMom.Phi());
  p_mGlobPt[cutSelection][9]->Fill(runIdenx, gMom.Pt());
  p_mGlobEta[cutSelection][9]->Fill(runIdenx, gMom.Eta());
  p_mGlobPhi[cutSelection][9]->Fill(runIdenx, gMom.Phi());
  p_mDEdx[cutSelection][9]->Fill(runIdenx, dEdx);
  if ( beta >= 0.0 )p_mBetaInv[cutSelection][9]->Fill(runIdenx, 1.0/beta);
  if ( mass2 != -999.9 ) p_mPrimaryMass2[cutSelection][9]->Fill(runIdenx, mass2);
  p_mNHitsMax[cutSelection][9]->Fill(runIdenx, nHitsMax);
  p_mNHitsDEdx[cutSelection][9]->Fill(runIdenx, nHitsDEdx);
  p_mNHitsRatio[cutSelection][9]->Fill(runIdenx, (float)nHitsFit/(float)nHitsMax);
}

void StRunQAProManager::fillRunQA_Track_PID(int triggerBin, int runIdenx, float nSigmaPi, float nSigmaK, float nSigmaP, float nSigmaE, float mass2, int cutSelection)
{
  p_mNSigmaPi[cutSelection][triggerBin]->Fill(runIdenx, nSigmaPi);
  p_mNSigmaK[cutSelection][triggerBin]->Fill(runIdenx, nSigmaK);
  p_mNSigmaP[cutSelection][triggerBin]->Fill(runIdenx, nSigmaP);
  p_mNSigmaE[cutSelection][triggerBin]->Fill(runIdenx, nSigmaE);
  if(nSigmaPi > -2.5 && nSigmaPi < 2.5) p_mMass2Pi[cutSelection][triggerBin]->Fill(runIdenx, mass2);
  if(nSigmaK > -2.5 && nSigmaK < 2.5) p_mMass2K[cutSelection][triggerBin]->Fill(runIdenx, mass2);
  if(nSigmaP > -2.5 && nSigmaP < 2.5) p_mMass2P[cutSelection][triggerBin]->Fill(runIdenx, mass2);
  if(nSigmaE > -2.5 && nSigmaE < 2.5) p_mMass2E[cutSelection][triggerBin]->Fill(runIdenx, mass2);

  p_mNSigmaPi[cutSelection][9]->Fill(runIdenx, nSigmaPi);
  p_mNSigmaK[cutSelection][9]->Fill(runIdenx, nSigmaK);
  p_mNSigmaP[cutSelection][9]->Fill(runIdenx, nSigmaP);
  p_mNSigmaE[cutSelection][9]->Fill(runIdenx, nSigmaE);
  if(nSigmaPi > -2.5 && nSigmaPi < 2.5) p_mMass2Pi[cutSelection][9]->Fill(runIdenx, mass2);
  if(nSigmaK > -2.5 && nSigmaK < 2.5) p_mMass2K[cutSelection][9]->Fill(runIdenx, mass2);
  if(nSigmaP > -2.5 && nSigmaP < 2.5) p_mMass2P[cutSelection][9]->Fill(runIdenx, mass2);
  if(nSigmaE > -2.5 && nSigmaE < 2.5) p_mMass2E[cutSelection][9]->Fill(runIdenx, mass2);
}

void StRunQAProManager::writeRunQA()
{
  for(int i_cut = 0; i_cut < 2; ++i_cut)
  {
    for(int i_trig = 0; i_trig < 10; ++i_trig)
    {
      p_mRefMult[i_cut][i_trig]->Write();
      p_mGRefMult[i_cut][i_trig]->Write();
      p_mZdcX[i_cut][i_trig]->Write();
      p_mVz[i_cut][i_trig]->Write();
      p_mVr[i_cut][i_trig]->Write();
      p_mGDca[i_cut][i_trig]->Write();
      p_mNHitsFit[i_cut][i_trig]->Write();
      p_mPrimPt[i_cut][i_trig]->Write();
      p_mPrimEta[i_cut][i_trig]->Write();
      p_mPrimPhi[i_cut][i_trig]->Write();
      p_mGlobPt[i_cut][i_trig]->Write();
      p_mGlobEta[i_cut][i_trig]->Write();
      p_mGlobPhi[i_cut][i_trig]->Write();
      p_mDEdx[i_cut][i_trig]->Write();
      p_mBetaInv[i_cut][i_trig]->Write();
      p_mPrimaryMass2[i_cut][i_trig]->Write();
      p_mNHitsMax[i_cut][i_trig]->Write();
      p_mNHitsDEdx[i_cut][i_trig]->Write();
      p_mNHitsRatio[i_cut][i_trig]->Write();
      p_mNSigmaPi[i_cut][i_trig]->Write();
      p_mNSigmaK[i_cut][i_trig]->Write();
      p_mNSigmaP[i_cut][i_trig]->Write();
      p_mNSigmaE[i_cut][i_trig]->Write();
      p_mMass2Pi[i_cut][i_trig]->Write();   
      p_mMass2K[i_cut][i_trig]->Write(); 
      p_mMass2P[i_cut][i_trig]->Write(); 
      p_mMass2E[i_cut][i_trig]->Write();    
    }
  }
}
