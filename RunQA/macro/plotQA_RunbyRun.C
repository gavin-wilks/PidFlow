#include <string>
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TLegend.h"
#include "../StRoot/StRunQAMaker/StRunQACons.h"

using namespace std;

static const string mCutsQA[2] = {"Before","After"};

//static const string mParticle[4] = {"Pi","K","P","E"};
//static const string mCutsQA[2] = {"Before","After"};

void plotQA_RunbyRun(int energy = 0, string JobId = "028AC6C12A1F110244470E4CCDFC40FF")
{
  string inputfile = Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/AuAu%s/SpinAlignment/file_%s_RunQA_%s.root",runQA::mBeamEnergy[energy].c_str(),runQA::mBeamEnergy[energy].c_str(),JobId.c_str());
 
  TFile *File_InPut = TFile::Open(inputfile.c_str());

  TProfile *p_mRefMult[2][10]; // 0: before cuts | 1: after cuts
  TProfile *p_mGRefMult[2][10]; // 0-8 for different triggerID | 9 for all triggers
  TProfile *p_mZdcX[2][10];
  TProfile *p_mVz[2][10];
  TProfile *p_mVr[2][10];

  TProfile *p_mGDca[2][10];
  TProfile *p_mNHitsFit[2][10];
  TProfile *p_mPrimPt[2][10];
  TProfile *p_mPrimEta[2][10];
  TProfile *p_mPrimPhi[2][10];
  TProfile *p_mGlobPt[2][10];
  TProfile *p_mGlobEta[2][10];
  TProfile *p_mGlobPhi[2][10];
  TProfile *p_mDEdx[2][10];
  TProfile *p_mBetaInv[2][10];   
  TProfile *p_mPrimaryMass2[2][10]; 
  TProfile *p_mNHitsMax[2][10]; 
  TProfile *p_mNHitsDEdx[2][10]; 
  TProfile *p_mNHitsRatio[2][10];  
  TProfile *p_mNSigmaPi[2][10];
  TProfile *p_mNSigmaK[2][10];
  TProfile *p_mNSigmaP[2][10];
  TProfile *p_mNSigmaE[2][10];
  TProfile *p_mMass2Pi[2][10];
  TProfile *p_mMass2K[2][10];
  TProfile *p_mMass2P[2][10];
  TProfile *p_mMass2E[2][10];
   
  for(int i_cut = 0; i_cut < 2; ++i_cut)
  {
    for(int i_trig = 0; i_trig < 10; ++i_trig)
    {
      std::string ProName = Form("p_mRefMult%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      p_mRefMult[i_cut][i_trig] = (TProfile*)File_InPut->Get(ProName.c_str());

      ProName = Form("p_mGRefMult%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      p_mGRefMult[i_cut][i_trig] = (TProfile*)File_InPut->Get(ProName.c_str());

      ProName = Form("p_mZdcX%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      p_mZdcX[i_cut][i_trig] = (TProfile*)File_InPut->Get(ProName.c_str());

      ProName = Form("p_mVz%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      p_mVz[i_cut][i_trig] = (TProfile*)File_InPut->Get(ProName.c_str());

      ProName = Form("p_mVr%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      p_mVr[i_cut][i_trig] = (TProfile*)File_InPut->Get(ProName.c_str());

      ProName = Form("p_mGDca%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      p_mGDca[i_cut][i_trig] = (TProfile*)File_InPut->Get(ProName.c_str());

      ProName = Form("p_mNHitsFit%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      p_mNHitsFit[i_cut][i_trig] = (TProfile*)File_InPut->Get(ProName.c_str());

      ProName = Form("p_mPrimPt%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      p_mPrimPt[i_cut][i_trig] = (TProfile*)File_InPut->Get(ProName.c_str());

      ProName = Form("p_mPrimEta%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      p_mPrimEta[i_cut][i_trig] = (TProfile*)File_InPut->Get(ProName.c_str());

      ProName = Form("p_mPrimPhi%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      p_mPrimPhi[i_cut][i_trig] = (TProfile*)File_InPut->Get(ProName.c_str());

      ProName = Form("p_mGlobPt%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      p_mGlobPt[i_cut][i_trig] = (TProfile*)File_InPut->Get(ProName.c_str());

      ProName = Form("p_mGlobEta%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      p_mGlobEta[i_cut][i_trig] = (TProfile*)File_InPut->Get(ProName.c_str());

      ProName = Form("p_mGlobPhi%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      p_mGlobPhi[i_cut][i_trig] = (TProfile*)File_InPut->Get(ProName.c_str());
     
      ProName = Form("p_mDEdx%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      p_mDEdx[i_cut][i_trig] = (TProfile*)File_InPut->Get(ProName.c_str());

      ProName = Form("p_mBetaInv%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      p_mBetaInv[i_cut][i_trig] = (TProfile*)File_InPut->Get(ProName.c_str());

      ProName = Form("p_mPrimaryMass2%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      p_mPrimaryMass2[i_cut][i_trig] = (TProfile*)File_InPut->Get(ProName.c_str());

      ProName = Form("p_mNHitsMax%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      p_mNHitsMax[i_cut][i_trig] = (TProfile*)File_InPut->Get(ProName.c_str());

      ProName = Form("p_mNHitsDEdx%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      p_mNHitsDEdx[i_cut][i_trig] = (TProfile*)File_InPut->Get(ProName.c_str());
  
      ProName = Form("p_mNHitsRatio%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      p_mNHitsRatio[i_cut][i_trig] = (TProfile*)File_InPut->Get(ProName.c_str());    

      ProName = Form("p_mNSigmaPi%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      p_mNSigmaPi[i_cut][i_trig] = (TProfile*)File_InPut->Get(ProName.c_str());

      ProName = Form("p_mNSigmaK%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      p_mNSigmaK[i_cut][i_trig] = (TProfile*)File_InPut->Get(ProName.c_str());

      ProName = Form("p_mNSigmaP%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      p_mNSigmaP[i_cut][i_trig] = (TProfile*)File_InPut->Get(ProName.c_str());

      ProName = Form("p_mNSigmaE%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      p_mNSigmaE[i_cut][i_trig] = (TProfile*)File_InPut->Get(ProName.c_str());
      
      ProName = Form("p_mMass2Pi%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      p_mMass2Pi[i_cut][i_trig] = (TProfile*)File_InPut->Get(ProName.c_str());

      ProName = Form("p_mMass2K%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      p_mMass2K[i_cut][i_trig] = (TProfile*)File_InPut->Get(ProName.c_str());

      ProName = Form("p_mMass2P%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      p_mMass2P[i_cut][i_trig] = (TProfile*)File_InPut->Get(ProName.c_str());

      ProName = Form("p_mMass2E%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      p_mMass2E[i_cut][i_trig] = (TProfile*)File_InPut->Get(ProName.c_str());


    }
  }

  const int numOfTriggers = 6;
  //const int triggerID[numOfTriggers] = {650001,650011,650021,650031,650041,650051};
  const int triggerID[numOfTriggers] = {640001,640011,640021,640031,640041,640051};
  const int MarkerColor[numOfTriggers] = {7,6,4,2,16,46};

  TLegend *leg = new TLegend(0.7,0.6,0.9,0.9);

  //---------------------
  TCanvas *c_RunQA_RefMult = new TCanvas("c_RunQA_RefMult","c_RunQA_RefMult",10,10,800,400);
  c_RunQA_RefMult->cd()->SetLeftMargin(0.1);
  c_RunQA_RefMult->cd()->SetRightMargin(0.1);
  c_RunQA_RefMult->cd()->SetBottomMargin(0.1);
  c_RunQA_RefMult->cd()->SetGrid(0,0);
  c_RunQA_RefMult->cd()->SetTicks(1,1);

  p_mRefMult[1][9]->SetTitle("refMult vs. runIndex");
  p_mRefMult[1][9]->SetStats(0);
  p_mRefMult[1][9]->SetMarkerColor(1);
  p_mRefMult[1][9]->SetMarkerStyle(20);
  p_mRefMult[1][9]->SetMarkerSize(1.0);
  p_mRefMult[1][9]->SetLineColor(1);
  p_mRefMult[1][9]->GetXaxis()->SetTitle("runIndex");
  p_mRefMult[1][9]->GetYaxis()->SetTitle("<refMult>");
  //p_mRefMult[1][9]->GetYaxis()->SetRangeUser(0,400);
  p_mRefMult[1][9]->Draw("pE");
  leg->AddEntry(p_mRefMult[1][9],"All Triggers","P");

  for(int i_trig = 0; i_trig < numOfTriggers; ++i_trig)
  {
    p_mRefMult[1][i_trig]->SetMarkerColor(MarkerColor[i_trig]);
    p_mRefMult[1][i_trig]->SetMarkerStyle(24);
    p_mRefMult[1][i_trig]->SetMarkerSize(0.8);
    p_mRefMult[1][i_trig]->SetLineColor(MarkerColor[i_trig]);
    if(p_mRefMult[1][i_trig]->GetEntries() > 0) p_mRefMult[1][i_trig]->Draw("pE same");
    leg->AddEntry(p_mRefMult[1][i_trig],Form("%d",triggerID[i_trig]),"P");
  }

  leg->Draw("same");
  c_RunQA_RefMult->SaveAs(Form("./figures/%s/c_RunQA_RefMult.pdf",runQA::mBeamEnergy[energy].c_str()));
  //---------------------
  TCanvas *c_RunQA_gRefMult = new TCanvas("c_RunQA_gRefMult","c_RunQA_gRefMult",10,10,800,400);
  c_RunQA_gRefMult->cd()->SetLeftMargin(0.1);
  c_RunQA_gRefMult->cd()->SetRightMargin(0.1);
  c_RunQA_gRefMult->cd()->SetBottomMargin(0.1);
  c_RunQA_gRefMult->cd()->SetGrid(0,0);
  c_RunQA_gRefMult->cd()->SetTicks(1,1);

  p_mGRefMult[1][9]->SetTitle("grefMult vs. runIndex");
  p_mGRefMult[1][9]->SetStats(0);
  p_mGRefMult[1][9]->SetMarkerColor(1);
  p_mGRefMult[1][9]->SetMarkerStyle(20);
  p_mGRefMult[1][9]->SetMarkerSize(1.0);
  p_mGRefMult[1][9]->SetLineColor(1);
  p_mGRefMult[1][9]->GetXaxis()->SetTitle("runIndex");
  p_mGRefMult[1][9]->GetYaxis()->SetTitle("<grefMult>");
  //p_mGRefMult[1][9]->GetYaxis()->SetRangeUser(0,400);
  p_mGRefMult[1][9]->Draw("pE");

  for(int i_trig = 0; i_trig < numOfTriggers; ++i_trig)
  {
    p_mGRefMult[1][i_trig]->SetMarkerColor(MarkerColor[i_trig]);
    p_mGRefMult[1][i_trig]->SetMarkerStyle(24);
    p_mGRefMult[1][i_trig]->SetMarkerSize(0.8);
    p_mGRefMult[1][i_trig]->SetLineColor(MarkerColor[i_trig]);
    if(p_mGRefMult[1][i_trig]->GetEntries() > 0) p_mGRefMult[1][i_trig]->Draw("pE same");
  }

  leg->Draw("same");
  c_RunQA_gRefMult->SaveAs(Form("./figures/%s/c_RunQA_gRefMult.pdf",runQA::mBeamEnergy[energy].c_str()));
  //---------------------
  TCanvas *c_RunQA_ZdcX = new TCanvas("c_RunQA_ZdcX","c_RunQA_ZdcX",10,10,800,400);
  c_RunQA_ZdcX->cd()->SetLeftMargin(0.1);
  c_RunQA_ZdcX->cd()->SetRightMargin(0.1);
  c_RunQA_ZdcX->cd()->SetBottomMargin(0.1);
  c_RunQA_ZdcX->cd()->SetGrid(0,0);
  c_RunQA_ZdcX->cd()->SetTicks(1,1);

  p_mZdcX[1][9]->SetTitle("ZdcX vs. runIndex");
  p_mZdcX[1][9]->SetStats(0);
  p_mZdcX[1][9]->SetMarkerColor(1);
  p_mZdcX[1][9]->SetMarkerStyle(20);
  p_mZdcX[1][9]->SetMarkerSize(1.0);
  p_mZdcX[1][9]->SetLineColor(1);
  p_mZdcX[1][9]->GetXaxis()->SetTitle("runIndex");
  p_mZdcX[1][9]->GetYaxis()->SetTitle("<ZdcX>");
  // p_mZdcX[1][9]->GetYaxis()->SetRangeUser(0,400);
  p_mZdcX[1][9]->Draw("pE");

  for(int i_trig = 0; i_trig < numOfTriggers; ++i_trig)
  {
    p_mZdcX[1][i_trig]->SetMarkerColor(MarkerColor[i_trig]);
    p_mZdcX[1][i_trig]->SetMarkerStyle(24);
    p_mZdcX[1][i_trig]->SetMarkerSize(0.8);
    p_mZdcX[1][i_trig]->SetLineColor(MarkerColor[i_trig]);
    if(p_mZdcX[1][i_trig]->GetEntries() > 0) p_mZdcX[1][i_trig]->Draw("pE same");
  } 

  leg->Draw("same");
  c_RunQA_ZdcX->SaveAs(Form("./figures/%s/c_RunQA_ZdcX.pdf",runQA::mBeamEnergy[energy].c_str()));
 
  float y_ll[10]={0.0};
  float y_ul[10]={1000.0,1000.0,1000.0,1000.0,1000.0,1000.0};
  float x_ll[10]={0.0,50.0,150.0,150.0,550.0,900.0,0.0,0.0,0.0,0.0};
  float x_ul[10]={100.0,200.0,250.0,550.0,950.0,1200.0,0.0,0.0,0.0,0.0};

  TProfile *p_mZdcX_wLimits[10];
  for(int i_trig = 0; i_trig < 10; ++i_trig)
  {
    c_RunQA_ZdcX->cd()->Clear();
    p_mZdcX_wLimits[i_trig] = (TProfile*)p_mZdcX[1][i_trig]->Clone(Form("h_ZdcX_wLimits_trig%d", i_trig));
    p_mZdcX_wLimits[i_trig]->GetXaxis()->SetRangeUser(x_ll[i_trig],x_ul[i_trig]);
    p_mZdcX_wLimits[i_trig]->GetYaxis()->SetRangeUser(y_ll[i_trig],y_ul[i_trig]);
    p_mZdcX_wLimits[i_trig]->Draw("pE");
    if(energy==1&&i_trig<6) c_RunQA_ZdcX->SaveAs(Form("./figures/%s/c_RunQA_ZdcX_Trigger%d.pdf",runQA::mBeamEnergy[energy].c_str(),i_trig));
  } 

  //---------------------
  TCanvas *c_RunQA_Vz = new TCanvas("c_RunQA_Vz","c_RunQA_Vz",10,10,800,400);
  c_RunQA_Vz->cd()->SetLeftMargin(0.1);
  c_RunQA_Vz->cd()->SetRightMargin(0.1);
  c_RunQA_Vz->cd()->SetBottomMargin(0.1);
  c_RunQA_Vz->cd()->SetGrid(0,0);
  c_RunQA_Vz->cd()->SetTicks(1,1);

  p_mVz[1][9]->SetTitle("Vz vs. runIndex");
  p_mVz[1][9]->SetStats(0);
  p_mVz[1][9]->SetMarkerColor(1);
  p_mVz[1][9]->SetMarkerStyle(20);
  p_mVz[1][9]->SetMarkerSize(1.0);
  p_mVz[1][9]->SetLineColor(1);
  p_mVz[1][9]->GetXaxis()->SetTitle("runIndex");
  p_mVz[1][9]->GetYaxis()->SetTitle("<Vz>");
  //p_mVz[1][9]->GetYaxis()->SetRangeUser(-3.0,3.0);
  p_mVz[1][9]->Draw("pE");

  for(int i_trig = 0; i_trig < numOfTriggers; ++i_trig)
  {
    p_mVz[1][i_trig]->SetMarkerColor(MarkerColor[i_trig]);
    p_mVz[1][i_trig]->SetMarkerStyle(24);
    p_mVz[1][i_trig]->SetMarkerSize(0.8);
    p_mVz[1][i_trig]->SetLineColor(MarkerColor[i_trig]);
    if(p_mVz[1][i_trig]->GetEntries() > 0) p_mVz[1][i_trig]->Draw("pE same");
  }

  leg->Draw("same");
  c_RunQA_Vz->SaveAs(Form("./figures/%s/c_RunQA_Vz.pdf",runQA::mBeamEnergy[energy].c_str()));
  //---------------------
  TCanvas *c_RunQA_Vr = new TCanvas("c_RunQA_Vr","c_RunQA_Vr",10,10,800,400);
  c_RunQA_Vr->cd()->SetLeftMargin(0.1);
  c_RunQA_Vr->cd()->SetRightMargin(0.1);
  c_RunQA_Vr->cd()->SetBottomMargin(0.1);
  c_RunQA_Vr->cd()->SetGrid(0,0);
  c_RunQA_Vr->cd()->SetTicks(1,1);

  p_mVr[1][9]->SetTitle("Vr vs. runIndex");
  p_mVr[1][9]->SetStats(0);
  p_mVr[1][9]->SetMarkerColor(1);
  p_mVr[1][9]->SetMarkerStyle(20);
  p_mVr[1][9]->SetMarkerSize(1.0);
  p_mVr[1][9]->SetLineColor(1);
  p_mVr[1][9]->GetXaxis()->SetTitle("runIndex");
  p_mVr[1][9]->GetYaxis()->SetTitle("<Vr>");
  //p_mVr[1][9]->GetYaxis()->SetRangeUser(0.1,0.5);
  p_mVr[1][9]->Draw("pE");

  for(int i_trig = 0; i_trig < numOfTriggers; ++i_trig)
  {
    p_mVr[1][i_trig]->SetMarkerColor(MarkerColor[i_trig]);
    p_mVr[1][i_trig]->SetMarkerStyle(24);
    p_mVr[1][i_trig]->SetMarkerSize(0.8);
    p_mVr[1][i_trig]->SetLineColor(MarkerColor[i_trig]);
    if(p_mVr[1][i_trig]->GetEntries() > 0) p_mVr[1][i_trig]->Draw("pE same");
  }

  leg->Draw("same");
  c_RunQA_Vr->SaveAs(Form("./figures/%s/c_RunQA_Vr.pdf",runQA::mBeamEnergy[energy].c_str()));
  //---------------------
  TCanvas *c_RunQA_gDca = new TCanvas("c_RunQA_gDca","c_RunQA_gDca",10,10,800,400);
  c_RunQA_gDca->cd()->SetLeftMargin(0.1);
  c_RunQA_gDca->cd()->SetRightMargin(0.1);
  c_RunQA_gDca->cd()->SetBottomMargin(0.1);
  c_RunQA_gDca->cd()->SetGrid(0,0);
  c_RunQA_gDca->cd()->SetTicks(1,1);

  p_mGDca[1][9]->SetTitle("gDca vs. runIndex");
  p_mGDca[1][9]->SetStats(0);
  p_mGDca[1][9]->SetMarkerColor(1);
  p_mGDca[1][9]->SetMarkerStyle(20);
  p_mGDca[1][9]->SetMarkerSize(1.0);
  p_mGDca[1][9]->SetLineColor(1);
  p_mGDca[1][9]->GetXaxis()->SetTitle("runIndex");
  p_mGDca[1][9]->GetYaxis()->SetTitle("<gDca>");
  //p_mGDca[1][9]->GetYaxis()->SetRangeUser(0.2,0.6);
  p_mGDca[1][9]->Draw("pE");

  for(int i_trig = 0; i_trig < numOfTriggers; ++i_trig)
  {
    p_mGDca[1][i_trig]->SetMarkerColor(MarkerColor[i_trig]);
    p_mGDca[1][i_trig]->SetMarkerStyle(24);
    p_mGDca[1][i_trig]->SetMarkerSize(0.8);
    p_mGDca[1][i_trig]->SetLineColor(MarkerColor[i_trig]);
    if(p_mGDca[1][i_trig]->GetEntries() > 0) p_mGDca[1][i_trig]->Draw("pE same");
  }

  leg->Draw("same");
  c_RunQA_gDca->SaveAs(Form("./figures/%s/c_RunQA_gDca.pdf",runQA::mBeamEnergy[energy].c_str()));
  //---------------------
  TCanvas *c_RunQA_nHitsFit = new TCanvas("c_RunQA_nHitsFit","c_RunQA_nHitsFit",10,10,800,400);
  c_RunQA_nHitsFit->cd()->SetLeftMargin(0.1);
  c_RunQA_nHitsFit->cd()->SetRightMargin(0.1);
  c_RunQA_nHitsFit->cd()->SetBottomMargin(0.1);
  c_RunQA_nHitsFit->cd()->SetGrid(0,0);
  c_RunQA_nHitsFit->cd()->SetTicks(1,1);

  p_mNHitsFit[1][9]->SetTitle("nHitsFit vs. runIndex");
  p_mNHitsFit[1][9]->SetStats(0);
  p_mNHitsFit[1][9]->SetMarkerColor(1);
  p_mNHitsFit[1][9]->SetMarkerStyle(20);
  p_mNHitsFit[1][9]->SetMarkerSize(1.0);
  p_mNHitsFit[1][9]->SetLineColor(1);
  p_mNHitsFit[1][9]->GetXaxis()->SetTitle("runIndex");
  p_mNHitsFit[1][9]->GetYaxis()->SetTitle("<nHitsFit>");
  //p_mNHitsFit[1][9]->GetYaxis()->SetRangeUser(25,40);
  p_mNHitsFit[1][9]->Draw("pE");

  for(int i_trig = 0; i_trig < numOfTriggers; ++i_trig)
  {
    p_mNHitsFit[1][i_trig]->SetMarkerColor(MarkerColor[i_trig]);
    p_mNHitsFit[1][i_trig]->SetMarkerStyle(24);
    p_mNHitsFit[1][i_trig]->SetMarkerSize(0.8);
    p_mNHitsFit[1][i_trig]->SetLineColor(MarkerColor[i_trig]);
    if(p_mNHitsFit[1][i_trig]->GetEntries() > 0) p_mNHitsFit[1][i_trig]->Draw("pE same");
  }

  leg->Draw("same");
  c_RunQA_nHitsFit->SaveAs(Form("./figures/%s/c_RunQA_nHitsFit.pdf",runQA::mBeamEnergy[energy].c_str()));
  //---------------------
  TCanvas *c_RunQA_primPt = new TCanvas("c_RunQA_primPt","c_RunQA_primPt",10,10,800,400);
  c_RunQA_primPt->cd()->SetLeftMargin(0.1);
  c_RunQA_primPt->cd()->SetRightMargin(0.1);
  c_RunQA_primPt->cd()->SetBottomMargin(0.1);
  c_RunQA_primPt->cd()->SetGrid(0,0);
  c_RunQA_primPt->cd()->SetTicks(1,1);

  p_mPrimPt[1][9]->SetTitle("primPt vs. runIndex");
  p_mPrimPt[1][9]->SetStats(0);
  p_mPrimPt[1][9]->SetMarkerColor(1);
  p_mPrimPt[1][9]->SetMarkerStyle(20);
  p_mPrimPt[1][9]->SetMarkerSize(1.0);
  p_mPrimPt[1][9]->SetLineColor(1);
  p_mPrimPt[1][9]->GetXaxis()->SetTitle("runIndex");
  p_mPrimPt[1][9]->GetYaxis()->SetTitle("<p_{T}^{prim}>");
  //p_mPrimPt[1][9]->GetYaxis()->SetRangeUser(0.4,0.8);
  p_mPrimPt[1][9]->Draw("pE");

  for(int i_trig = 0; i_trig < numOfTriggers; ++i_trig)
  {
    p_mPrimPt[1][i_trig]->SetMarkerColor(MarkerColor[i_trig]);
    p_mPrimPt[1][i_trig]->SetMarkerStyle(24);
    p_mPrimPt[1][i_trig]->SetMarkerSize(0.8);
    p_mPrimPt[1][i_trig]->SetLineColor(MarkerColor[i_trig]);
    if(p_mPrimPt[1][i_trig]->GetEntries() > 0) p_mPrimPt[1][i_trig]->Draw("pE same");
  }

  leg->Draw("same");
  c_RunQA_primPt->SaveAs(Form("./figures/%s/c_RunQA_primPt.pdf",runQA::mBeamEnergy[energy].c_str()));
  //---------------------
  TCanvas *c_RunQA_primEta = new TCanvas("c_RunQA_primEta","c_RunQA_primEta",10,10,800,400);
  c_RunQA_primEta->cd()->SetLeftMargin(0.1);
  c_RunQA_primEta->cd()->SetRightMargin(0.1);
  c_RunQA_primEta->cd()->SetBottomMargin(0.1);
  c_RunQA_primEta->cd()->SetGrid(0,0);
  c_RunQA_primEta->cd()->SetTicks(1,1);

  p_mPrimEta[1][9]->SetTitle("primEta vs. runIndex");
  p_mPrimEta[1][9]->SetStats(0);
  p_mPrimEta[1][9]->SetMarkerColor(1);
  p_mPrimEta[1][9]->SetMarkerStyle(20);
  p_mPrimEta[1][9]->SetMarkerSize(1.0);
  p_mPrimEta[1][9]->SetLineColor(1);
  p_mPrimEta[1][9]->GetXaxis()->SetTitle("runIndex");
  p_mPrimEta[1][9]->GetYaxis()->SetTitle("<#eta^{prim}>");
  //p_mPrimEta[1][9]->GetYaxis()->SetRangeUser(-0.05,0.10);
  p_mPrimEta[1][9]->Draw("pE");

  for(int i_trig = 0; i_trig < numOfTriggers; ++i_trig)
  {
    p_mPrimEta[1][i_trig]->SetMarkerColor(MarkerColor[i_trig]);
    p_mPrimEta[1][i_trig]->SetMarkerStyle(24);
    p_mPrimEta[1][i_trig]->SetMarkerSize(0.8);
    p_mPrimEta[1][i_trig]->SetLineColor(MarkerColor[i_trig]);
    if(p_mPrimEta[1][i_trig]->GetEntries() > 0) p_mPrimEta[1][i_trig]->Draw("pE same");
  }

  leg->Draw("same");
  c_RunQA_primEta->SaveAs(Form("./figures/%s/c_RunQA_primEta.pdf",runQA::mBeamEnergy[energy].c_str()));
  //---------------------
  TCanvas *c_RunQA_primPhi = new TCanvas("c_RunQA_primPhi","c_RunQA_primPhi",10,10,800,400);
  c_RunQA_primPhi->cd()->SetLeftMargin(0.1);
  c_RunQA_primPhi->cd()->SetRightMargin(0.1);
  c_RunQA_primPhi->cd()->SetBottomMargin(0.1);
  c_RunQA_primPhi->cd()->SetGrid(0,0);
  c_RunQA_primPhi->cd()->SetTicks(1,1);

  p_mPrimPhi[1][9]->SetTitle("primPhi vs. runIndex");
  p_mPrimPhi[1][9]->SetStats(0);
  p_mPrimPhi[1][9]->SetMarkerColor(1);
  p_mPrimPhi[1][9]->SetMarkerStyle(20);
  p_mPrimPhi[1][9]->SetMarkerSize(1.0);
  p_mPrimPhi[1][9]->SetLineColor(1);
  p_mPrimPhi[1][9]->GetXaxis()->SetTitle("runIndex");
  p_mPrimPhi[1][9]->GetYaxis()->SetTitle("<#phi^{prim}>");
  //p_mPrimPhi[1][9]->GetYaxis()->SetRangeUser(-0.1,0.50);
  p_mPrimPhi[1][9]->Draw("pE");

  for(int i_trig = 0; i_trig < numOfTriggers; ++i_trig)
  {
    p_mPrimPhi[1][i_trig]->SetMarkerColor(MarkerColor[i_trig]);
    p_mPrimPhi[1][i_trig]->SetMarkerStyle(24);
    p_mPrimPhi[1][i_trig]->SetMarkerSize(0.8);
    p_mPrimPhi[1][i_trig]->SetLineColor(MarkerColor[i_trig]);
    if(p_mPrimPhi[1][i_trig]->GetEntries() > 0) p_mPrimPhi[1][i_trig]->Draw("pE same");
  }

  leg->Draw("same");
  c_RunQA_primPhi->SaveAs(Form("./figures/%s/c_RunQA_primPhi.pdf",runQA::mBeamEnergy[energy].c_str()));
  //---------------------
  TCanvas *c_RunQA_globPt = new TCanvas("c_RunQA_globPt","c_RunQA_globPt",10,10,800,400);
  c_RunQA_globPt->cd()->SetLeftMargin(0.1);
  c_RunQA_globPt->cd()->SetRightMargin(0.1);
  c_RunQA_globPt->cd()->SetBottomMargin(0.1);
  c_RunQA_globPt->cd()->SetGrid(0,0);
  c_RunQA_globPt->cd()->SetTicks(1,1);

  p_mGlobPt[1][9]->SetTitle("globPt vs. runIndex");
  p_mGlobPt[1][9]->SetStats(0);
  p_mGlobPt[1][9]->SetMarkerColor(1);
  p_mGlobPt[1][9]->SetMarkerStyle(20);
  p_mGlobPt[1][9]->SetMarkerSize(1.0);
  p_mGlobPt[1][9]->SetLineColor(1);
  p_mGlobPt[1][9]->GetXaxis()->SetTitle("runIndex");
  p_mGlobPt[1][9]->GetYaxis()->SetTitle("<p_{T}^{glob}>");
  p_mGlobPt[1][9]->GetYaxis()->SetRangeUser(0.0,1.0);
  p_mGlobPt[1][9]->Draw("pE");

  for(int i_trig = 0; i_trig < numOfTriggers; ++i_trig)
  {
    p_mGlobPt[1][i_trig]->SetMarkerColor(MarkerColor[i_trig]);
    p_mGlobPt[1][i_trig]->SetMarkerStyle(24);
    p_mGlobPt[1][i_trig]->SetMarkerSize(0.8);
    p_mGlobPt[1][i_trig]->SetLineColor(MarkerColor[i_trig]);
    if(p_mGlobPt[1][i_trig]->GetEntries() > 0) p_mGlobPt[1][i_trig]->Draw("pE same");
  }

  leg->Draw("same");
  c_RunQA_globPt->SaveAs(Form("./figures/%s/c_RunQA_globPt.pdf",runQA::mBeamEnergy[energy].c_str()));
  //---------------------
  TCanvas *c_RunQA_globEta = new TCanvas("c_RunQA_globEta","c_RunQA_globEta",10,10,800,400);
  c_RunQA_globEta->cd()->SetLeftMargin(0.1);
  c_RunQA_globEta->cd()->SetRightMargin(0.1);
  c_RunQA_globEta->cd()->SetBottomMargin(0.1);
  c_RunQA_globEta->cd()->SetGrid(0,0);
  c_RunQA_globEta->cd()->SetTicks(1,1);

  p_mGlobEta[1][9]->SetTitle("globEta vs. runIndex");
  p_mGlobEta[1][9]->SetStats(0);
  p_mGlobEta[1][9]->SetMarkerColor(1);
  p_mGlobEta[1][9]->SetMarkerStyle(20);
  p_mGlobEta[1][9]->SetMarkerSize(1.0);
  p_mGlobEta[1][9]->SetLineColor(1);
  p_mGlobEta[1][9]->GetXaxis()->SetTitle("runIndex");
  p_mGlobEta[1][9]->GetYaxis()->SetTitle("<#eta^{glob}>");
  //p_mGlobEta[1][9]->GetYaxis()->SetRangeUser(-0.05,0.10);
  p_mGlobEta[1][9]->Draw("pE");

  for(int i_trig = 0; i_trig < numOfTriggers; ++i_trig)
  {
    p_mGlobEta[1][i_trig]->SetMarkerColor(MarkerColor[i_trig]);
    p_mGlobEta[1][i_trig]->SetMarkerStyle(24);
    p_mGlobEta[1][i_trig]->SetMarkerSize(0.8);
    p_mGlobEta[1][i_trig]->SetLineColor(MarkerColor[i_trig]);
    if(p_mGlobEta[1][i_trig]->GetEntries() > 0) p_mGlobEta[1][i_trig]->Draw("pE same");
  }

  leg->Draw("same");
  c_RunQA_globEta->SaveAs(Form("./figures/%s/c_RunQA_globEta.pdf",runQA::mBeamEnergy[energy].c_str()));
  //---------------------
  TCanvas *c_RunQA_globPhi = new TCanvas("c_RunQA_globPhi","c_RunQA_globPhi",10,10,800,400);
  c_RunQA_globPhi->cd()->SetLeftMargin(0.1);
  c_RunQA_globPhi->cd()->SetRightMargin(0.1);
  c_RunQA_globPhi->cd()->SetBottomMargin(0.1);
  c_RunQA_globPhi->cd()->SetGrid(0,0);
  c_RunQA_globPhi->cd()->SetTicks(1,1);

  p_mGlobPhi[1][9]->SetTitle("globPhi vs. runIndex");
  p_mGlobPhi[1][9]->SetStats(0);
  p_mGlobPhi[1][9]->SetMarkerColor(1);
  p_mGlobPhi[1][9]->SetMarkerStyle(20);
  p_mGlobPhi[1][9]->SetMarkerSize(1.0);
  p_mGlobPhi[1][9]->SetLineColor(1);
  p_mGlobPhi[1][9]->GetXaxis()->SetTitle("runIndex");
  p_mGlobPhi[1][9]->GetYaxis()->SetTitle("<#phi^{glob}>");
  //p_mGlobPhi[1][9]->GetYaxis()->SetRangeUser(-0.1,0.50);
  p_mGlobPhi[1][9]->Draw("pE");

  for(int i_trig = 0; i_trig < numOfTriggers; ++i_trig)
  {
    p_mGlobPhi[1][i_trig]->SetMarkerColor(MarkerColor[i_trig]);
    p_mGlobPhi[1][i_trig]->SetMarkerStyle(24);
    p_mGlobPhi[1][i_trig]->SetMarkerSize(0.8);
    p_mGlobPhi[1][i_trig]->SetLineColor(MarkerColor[i_trig]);
    if(p_mGlobPhi[1][i_trig]->GetEntries() > 0) p_mGlobPhi[1][i_trig]->Draw("pE same");
  }

  leg->Draw("same");
  c_RunQA_globPhi->SaveAs(Form("./figures/%s/c_RunQA_globPhi.pdf",runQA::mBeamEnergy[energy].c_str()));
  //---------------------
  TCanvas *c_RunQA_dEdx = new TCanvas("c_RunQA_dEdx","c_RunQA_dEdx",10,10,800,400);
  c_RunQA_dEdx->cd()->SetLeftMargin(0.1);
  c_RunQA_dEdx->cd()->SetRightMargin(0.1);
  c_RunQA_dEdx->cd()->SetBottomMargin(0.1);
  c_RunQA_dEdx->cd()->SetGrid(0,0);
  c_RunQA_dEdx->cd()->SetTicks(1,1);

  p_mDEdx[1][9]->SetTitle("dE/dx vs. runIndex");
  p_mDEdx[1][9]->SetStats(0);
  p_mDEdx[1][9]->SetMarkerColor(1);
  p_mDEdx[1][9]->SetMarkerStyle(20);
  p_mDEdx[1][9]->SetMarkerSize(1.0);
  p_mDEdx[1][9]->SetLineColor(1);
  p_mDEdx[1][9]->GetXaxis()->SetTitle("runIndex");
  p_mDEdx[1][9]->GetYaxis()->SetTitle("<dE/dx>");
  //p_mGlobPhi[1][9]->GetYaxis()->SetRangeUser(-0.1,0.50);
  p_mDEdx[1][9]->Draw("pE");

  for(int i_trig = 0; i_trig < numOfTriggers; ++i_trig)
  {
    p_mDEdx[1][i_trig]->SetMarkerColor(MarkerColor[i_trig]);
    p_mDEdx[1][i_trig]->SetMarkerStyle(24);
    p_mDEdx[1][i_trig]->SetMarkerSize(0.8);
    p_mDEdx[1][i_trig]->SetLineColor(MarkerColor[i_trig]);
    if(p_mDEdx[1][i_trig]->GetEntries() > 0) p_mDEdx[1][i_trig]->Draw("pE same");
  }

  leg->Draw("same");
  c_RunQA_dEdx->SaveAs(Form("./figures/%s/c_RunQA_dEdx.pdf",runQA::mBeamEnergy[energy].c_str()));
  //---------------------
  TCanvas *c_RunQA_BetaInv = new TCanvas("c_RunQA_BetaInv","c_RunQA_BetaInv",10,10,800,400);
  c_RunQA_BetaInv->cd()->SetLeftMargin(0.1);
  c_RunQA_BetaInv->cd()->SetRightMargin(0.1);
  c_RunQA_BetaInv->cd()->SetBottomMargin(0.1);
  c_RunQA_BetaInv->cd()->SetGrid(0,0);
  c_RunQA_BetaInv->cd()->SetTicks(1,1);

  p_mBetaInv[1][9]->SetTitle("1/#beta vs. runIndex");
  p_mBetaInv[1][9]->SetStats(0);
  p_mBetaInv[1][9]->SetMarkerColor(1);
  p_mBetaInv[1][9]->SetMarkerStyle(20);
  p_mBetaInv[1][9]->SetMarkerSize(1.0);
  p_mBetaInv[1][9]->SetLineColor(1);
  p_mBetaInv[1][9]->GetXaxis()->SetTitle("runIndex");
  p_mBetaInv[1][9]->GetYaxis()->SetTitle("<1/#beta>");
  //p_mGlobPhi[1][9]->GetYaxis()->SetRangeUser(-0.1,0.50);
  p_mBetaInv[1][9]->Draw("pE");

  for(int i_trig = 0; i_trig < numOfTriggers; ++i_trig)
  {
    p_mBetaInv[1][i_trig]->SetMarkerColor(MarkerColor[i_trig]);
    p_mBetaInv[1][i_trig]->SetMarkerStyle(24);
    p_mBetaInv[1][i_trig]->SetMarkerSize(0.8);
    p_mBetaInv[1][i_trig]->SetLineColor(MarkerColor[i_trig]);
    if(p_mBetaInv[1][i_trig]->GetEntries() > 0) p_mBetaInv[1][i_trig]->Draw("pE same");
  }

  leg->Draw("same");
  c_RunQA_BetaInv->SaveAs(Form("./figures/%s/c_RunQA_BetaInv.pdf",runQA::mBeamEnergy[energy].c_str()));
  
  //---------------------
  TCanvas *c_RunQA_PrimaryMass2 = new TCanvas("c_RunQA_PrimaryMass2","c_RunQA_PrimaryMass2",10,10,800,400);
  c_RunQA_PrimaryMass2->cd()->SetLeftMargin(0.1);
  c_RunQA_PrimaryMass2->cd()->SetRightMargin(0.1);
  c_RunQA_PrimaryMass2->cd()->SetBottomMargin(0.1);
  c_RunQA_PrimaryMass2->cd()->SetGrid(0,0);
  c_RunQA_PrimaryMass2->cd()->SetTicks(1,1);

  p_mPrimaryMass2[1][9]->SetTitle("M^{2} vs. runIndex");
  p_mPrimaryMass2[1][9]->SetStats(0);
  p_mPrimaryMass2[1][9]->SetMarkerColor(1);
  p_mPrimaryMass2[1][9]->SetMarkerStyle(20);
  p_mPrimaryMass2[1][9]->SetMarkerSize(1.0);
  p_mPrimaryMass2[1][9]->SetLineColor(1);
  p_mPrimaryMass2[1][9]->GetXaxis()->SetTitle("runIndex");
  p_mPrimaryMass2[1][9]->GetYaxis()->SetTitle("<M^{2}>");
  //pPrimaryMass2_mGlobPhi[1][9]->GetYaxis()->SetRangeUser(-0.1,0.50);
  p_mPrimaryMass2[1][9]->Draw("pE");

  for(int i_trig = 0; i_trig < numOfTriggers; ++i_trig)
  {
    p_mPrimaryMass2[1][i_trig]->SetMarkerColor(MarkerColor[i_trig]);
    p_mPrimaryMass2[1][i_trig]->SetMarkerStyle(24);
    p_mPrimaryMass2[1][i_trig]->SetMarkerSize(0.8);
    p_mPrimaryMass2[1][i_trig]->SetLineColor(MarkerColor[i_trig]);
    if(p_mPrimaryMass2[1][i_trig]->GetEntries() > 0) p_mPrimaryMass2[1][i_trig]->Draw("pE same");
  }

  leg->Draw("same");
  c_RunQA_PrimaryMass2->SaveAs(Form("./figures/%s/c_RunQA_PrimaryMass2.pdf",runQA::mBeamEnergy[energy].c_str()));
  //---------------------
  TCanvas *c_RunQA_nHitsMax = new TCanvas("c_RunQA_nHitsMax","c_RunQA_nHitsMax",10,10,800,400);
  c_RunQA_nHitsMax->cd()->SetLeftMargin(0.1);
  c_RunQA_nHitsMax->cd()->SetRightMargin(0.1);
  c_RunQA_nHitsMax->cd()->SetBottomMargin(0.1);
  c_RunQA_nHitsMax->cd()->SetGrid(0,0);
  c_RunQA_nHitsMax->cd()->SetTicks(1,1);

  p_mNHitsMax[1][9]->SetTitle("nHitsMax vs. runIndex");
  p_mNHitsMax[1][9]->SetStats(0);
  p_mNHitsMax[1][9]->SetMarkerColor(1);
  p_mNHitsMax[1][9]->SetMarkerStyle(20);
  p_mNHitsMax[1][9]->SetMarkerSize(1.0);
  p_mNHitsMax[1][9]->SetLineColor(1);
  p_mNHitsMax[1][9]->GetXaxis()->SetTitle("runIndex");
  p_mNHitsMax[1][9]->GetYaxis()->SetTitle("<nHitsMax>");
  //pNHitsMax_mGlobPhi[1][9]->GetYaxis()->SetRangeUser(-0.1,0.50);
  p_mNHitsMax[1][9]->Draw("pE");

  for(int i_trig = 0; i_trig < numOfTriggers; ++i_trig)
  {
    p_mNHitsMax[1][i_trig]->SetMarkerColor(MarkerColor[i_trig]);
    p_mNHitsMax[1][i_trig]->SetMarkerStyle(24);
    p_mNHitsMax[1][i_trig]->SetMarkerSize(0.8);
    p_mNHitsMax[1][i_trig]->SetLineColor(MarkerColor[i_trig]);
    if(p_mNHitsMax[1][i_trig]->GetEntries() > 0) p_mNHitsMax[1][i_trig]->Draw("pE same");
  }

  leg->Draw("same");
  c_RunQA_nHitsMax->SaveAs(Form("./figures/%s/c_RunQA_nHitsMax.pdf",runQA::mBeamEnergy[energy].c_str()));
  //---------------------
  TCanvas *c_RunQA_nHitsDEdx = new TCanvas("c_RunQA_nHitsDEdx","c_RunQA_dEdx",10,10,800,400);
  c_RunQA_nHitsDEdx->cd()->SetLeftMargin(0.1);
  c_RunQA_nHitsDEdx->cd()->SetRightMargin(0.1);
  c_RunQA_nHitsDEdx->cd()->SetBottomMargin(0.1);
  c_RunQA_nHitsDEdx->cd()->SetGrid(0,0);
  c_RunQA_nHitsDEdx->cd()->SetTicks(1,1);

  p_mNHitsDEdx[1][9]->SetTitle("nHits dE/dx vs. runIndex");
  p_mNHitsDEdx[1][9]->SetStats(0);
  p_mNHitsDEdx[1][9]->SetMarkerColor(1);
  p_mNHitsDEdx[1][9]->SetMarkerStyle(20);
  p_mNHitsDEdx[1][9]->SetMarkerSize(1.0);
  p_mNHitsDEdx[1][9]->SetLineColor(1);
  p_mNHitsDEdx[1][9]->GetXaxis()->SetTitle("runIndex");
  p_mNHitsDEdx[1][9]->GetYaxis()->SetTitle("<nHits dE/dx>");
  //pNHitsDEdx_mGlobPhi[1][9]->GetYaxis()->SetRangeUser(-0.1,0.50);
  p_mNHitsDEdx[1][9]->Draw("pE");

  for(int i_trig = 0; i_trig < numOfTriggers; ++i_trig)
  {
    p_mNHitsDEdx[1][i_trig]->SetMarkerColor(MarkerColor[i_trig]);
    p_mNHitsDEdx[1][i_trig]->SetMarkerStyle(24);
    p_mNHitsDEdx[1][i_trig]->SetMarkerSize(0.8);
    p_mNHitsDEdx[1][i_trig]->SetLineColor(MarkerColor[i_trig]);
    if(p_mNHitsDEdx[1][i_trig]->GetEntries() > 0) p_mNHitsDEdx[1][i_trig]->Draw("pE same");
  }

  leg->Draw("same");
  c_RunQA_nHitsDEdx->SaveAs(Form("./figures/%s/c_RunQA_nHitsDEdx.pdf",runQA::mBeamEnergy[energy].c_str()));
  //---------------------
  TCanvas *c_RunQA_nHitsRatio = new TCanvas("c_RunQA_nHitsRatio","c_RunQA_dEdx",10,10,800,400);
  c_RunQA_nHitsRatio->cd()->SetLeftMargin(0.1);
  c_RunQA_nHitsRatio->cd()->SetRightMargin(0.1);
  c_RunQA_nHitsRatio->cd()->SetBottomMargin(0.1);
  c_RunQA_nHitsRatio->cd()->SetGrid(0,0);
  c_RunQA_nHitsRatio->cd()->SetTicks(1,1);

  p_mNHitsRatio[1][9]->SetTitle("nHitsRatio vs. runIndex");
  p_mNHitsRatio[1][9]->SetStats(0);
  p_mNHitsRatio[1][9]->SetMarkerColor(1);
  p_mNHitsRatio[1][9]->SetMarkerStyle(20);
  p_mNHitsRatio[1][9]->SetMarkerSize(1.0);
  p_mNHitsRatio[1][9]->SetLineColor(1);
  p_mNHitsRatio[1][9]->GetXaxis()->SetTitle("runIndex");
  p_mNHitsRatio[1][9]->GetYaxis()->SetTitle("<nHitsFit/nHitsMax>");
  //pNHitsRatio_mGlobPhi[1][9]->GetYaxis()->SetRangeUser(-0.1,0.50);
  p_mNHitsRatio[1][9]->Draw("pE");

  for(int i_trig = 0; i_trig < numOfTriggers; ++i_trig)
  {
    p_mNHitsRatio[1][i_trig]->SetMarkerColor(MarkerColor[i_trig]);
    p_mNHitsRatio[1][i_trig]->SetMarkerStyle(24);
    p_mNHitsRatio[1][i_trig]->SetMarkerSize(0.8);
    p_mNHitsRatio[1][i_trig]->SetLineColor(MarkerColor[i_trig]);
    if(p_mNHitsRatio[1][i_trig]->GetEntries() > 0) p_mNHitsRatio[1][i_trig]->Draw("pE same");
  }

  leg->Draw("same");
  c_RunQA_nHitsRatio->SaveAs(Form("./figures/%s/c_RunQA_nHitsRatio.pdf",runQA::mBeamEnergy[energy].c_str()));
  //---------------------
  TCanvas *c_RunQA_nSigmaPi = new TCanvas("c_RunQA_nSigmaPi","c_RunQA_dEdx",10,10,800,400);
  c_RunQA_nSigmaPi->cd()->SetLeftMargin(0.1);
  c_RunQA_nSigmaPi->cd()->SetRightMargin(0.1);
  c_RunQA_nSigmaPi->cd()->SetBottomMargin(0.1);
  c_RunQA_nSigmaPi->cd()->SetGrid(0,0);
  c_RunQA_nSigmaPi->cd()->SetTicks(1,1);

  p_mNSigmaPi[1][9]->SetTitle("nSigmaPi vs. runIndex");
  p_mNSigmaPi[1][9]->SetStats(0);
  p_mNSigmaPi[1][9]->SetMarkerColor(1);
  p_mNSigmaPi[1][9]->SetMarkerStyle(20);
  p_mNSigmaPi[1][9]->SetMarkerSize(1.0);
  p_mNSigmaPi[1][9]->SetLineColor(1);
  p_mNSigmaPi[1][9]->GetXaxis()->SetTitle("runIndex");
  p_mNSigmaPi[1][9]->GetYaxis()->SetTitle("<nSigmaPi>");
  //pNSigmaPi_mGlobPhi[1][9]->GetYaxis()->SetRangeUser(-0.1,0.50);
  p_mNSigmaPi[1][9]->Draw("pE");

  for(int i_trig = 0; i_trig < numOfTriggers; ++i_trig)
  {
    p_mNSigmaPi[1][i_trig]->SetMarkerColor(MarkerColor[i_trig]);
    p_mNSigmaPi[1][i_trig]->SetMarkerStyle(24);
    p_mNSigmaPi[1][i_trig]->SetMarkerSize(0.8);
    p_mNSigmaPi[1][i_trig]->SetLineColor(MarkerColor[i_trig]);
    if(p_mNSigmaPi[1][i_trig]->GetEntries() > 0) p_mNSigmaPi[1][i_trig]->Draw("pE same");
  }

  leg->Draw("same");
  c_RunQA_nSigmaPi->SaveAs(Form("./figures/%s/c_RunQA_nSigmaPi.pdf",runQA::mBeamEnergy[energy].c_str()));
  //---------------------
  TCanvas *c_RunQA_nSigmaK = new TCanvas("c_RunQA_nSigmaK","c_RunQA_dEdx",10,10,800,400);
  c_RunQA_nSigmaK->cd()->SetLeftMargin(0.1);
  c_RunQA_nSigmaK->cd()->SetRightMargin(0.1);
  c_RunQA_nSigmaK->cd()->SetBottomMargin(0.1);
  c_RunQA_nSigmaK->cd()->SetGrid(0,0);
  c_RunQA_nSigmaK->cd()->SetTicks(1,1);

  p_mNSigmaK[1][9]->SetTitle("nSigmaK vs. runIndex");
  p_mNSigmaK[1][9]->SetStats(0);
  p_mNSigmaK[1][9]->SetMarkerColor(1);
  p_mNSigmaK[1][9]->SetMarkerStyle(20);
  p_mNSigmaK[1][9]->SetMarkerSize(1.0);
  p_mNSigmaK[1][9]->SetLineColor(1);
  p_mNSigmaK[1][9]->GetXaxis()->SetTitle("runIndex");
  p_mNSigmaK[1][9]->GetYaxis()->SetTitle("<nSigmaK>");
  //pNSigmaK_mGlobPhi[1][9]->GetYaxis()->SetRangeUser(-0.1,0.50);
  p_mNSigmaK[1][9]->Draw("pE");

  for(int i_trig = 0; i_trig < numOfTriggers; ++i_trig)
  {
    p_mNSigmaK[1][i_trig]->SetMarkerColor(MarkerColor[i_trig]);
    p_mNSigmaK[1][i_trig]->SetMarkerStyle(24);
    p_mNSigmaK[1][i_trig]->SetMarkerSize(0.8);
    p_mNSigmaK[1][i_trig]->SetLineColor(MarkerColor[i_trig]);
    if(p_mNSigmaK[1][i_trig]->GetEntries() > 0) p_mNSigmaK[1][i_trig]->Draw("pE same");
  }

  leg->Draw("same");
  c_RunQA_nSigmaK->SaveAs(Form("./figures/%s/c_RunQA_nSigmaK.pdf",runQA::mBeamEnergy[energy].c_str()));
  //---------------------
  TCanvas *c_RunQA_nSigmaP = new TCanvas("c_RunQA_nSigmaP","c_RunQA_dEdx",10,10,800,400);
  c_RunQA_nSigmaP->cd()->SetLeftMargin(0.1);
  c_RunQA_nSigmaP->cd()->SetRightMargin(0.1);
  c_RunQA_nSigmaP->cd()->SetBottomMargin(0.1);
  c_RunQA_nSigmaP->cd()->SetGrid(0,0);
  c_RunQA_nSigmaP->cd()->SetTicks(1,1);

  p_mNSigmaP[1][9]->SetTitle("nSigmaP vs. runIndex");
  p_mNSigmaP[1][9]->SetStats(0);
  p_mNSigmaP[1][9]->SetMarkerColor(1);
  p_mNSigmaP[1][9]->SetMarkerStyle(20);
  p_mNSigmaP[1][9]->SetMarkerSize(1.0);
  p_mNSigmaP[1][9]->SetLineColor(1);
  p_mNSigmaP[1][9]->GetXaxis()->SetTitle("runIndex");
  p_mNSigmaP[1][9]->GetYaxis()->SetTitle("<nSigmaP>");
  //pNSigmaP_mGlobPhi[1][9]->GetYaxis()->SetRangeUser(-0.1,0.50);
  p_mNSigmaP[1][9]->Draw("pE");

  for(int i_trig = 0; i_trig < numOfTriggers; ++i_trig)
  {
    p_mNSigmaP[1][i_trig]->SetMarkerColor(MarkerColor[i_trig]);
    p_mNSigmaP[1][i_trig]->SetMarkerStyle(24);
    p_mNSigmaP[1][i_trig]->SetMarkerSize(0.8);
    p_mNSigmaP[1][i_trig]->SetLineColor(MarkerColor[i_trig]);
    if(p_mNSigmaP[1][i_trig]->GetEntries() > 0) p_mNSigmaP[1][i_trig]->Draw("pE same");
  }

  leg->Draw("same");
  c_RunQA_nSigmaP->SaveAs(Form("./figures/%s/c_RunQA_nSigmaP.pdf",runQA::mBeamEnergy[energy].c_str()));
  //---------------------
  TCanvas *c_RunQA_nSigmaE = new TCanvas("c_RunQA_nSigmaE","c_RunQA_dEdx",10,10,800,400);
  c_RunQA_nSigmaE->cd()->SetLeftMargin(0.1);
  c_RunQA_nSigmaE->cd()->SetRightMargin(0.1);
  c_RunQA_nSigmaE->cd()->SetBottomMargin(0.1);
  c_RunQA_nSigmaE->cd()->SetGrid(0,0);
  c_RunQA_nSigmaE->cd()->SetTicks(1,1);

  p_mNSigmaE[1][9]->SetTitle("nSigmaE vs. runIndex");
  p_mNSigmaE[1][9]->SetStats(0);
  p_mNSigmaE[1][9]->SetMarkerColor(1);
  p_mNSigmaE[1][9]->SetMarkerStyle(20);
  p_mNSigmaE[1][9]->SetMarkerSize(1.0);
  p_mNSigmaE[1][9]->SetLineColor(1);
  p_mNSigmaE[1][9]->GetXaxis()->SetTitle("runIndex");
  p_mNSigmaE[1][9]->GetYaxis()->SetTitle("<nSigmaE>");
  //pNSigmaE_mGlobPhi[1][9]->GetYaxis()->SetRangeUser(-0.1,0.50);
  p_mNSigmaE[1][9]->Draw("pE");

  for(int i_trig = 0; i_trig < numOfTriggers; ++i_trig)
  {
    p_mNSigmaE[1][i_trig]->SetMarkerColor(MarkerColor[i_trig]);
    p_mNSigmaE[1][i_trig]->SetMarkerStyle(24);
    p_mNSigmaE[1][i_trig]->SetMarkerSize(0.8);
    p_mNSigmaE[1][i_trig]->SetLineColor(MarkerColor[i_trig]);
    if(p_mNSigmaE[1][i_trig]->GetEntries() > 0) p_mNSigmaE[1][i_trig]->Draw("pE same");
  }

  leg->Draw("same");
  c_RunQA_nSigmaE->SaveAs(Form("./figures/%s/c_RunQA_nSigmaE.pdf",runQA::mBeamEnergy[energy].c_str()));
  //---------------------
  TCanvas *c_RunQA_Mass2Pi = new TCanvas("c_RunQA_Mass2Pi","c_RunQA_dEdx",10,10,800,400);
  c_RunQA_Mass2Pi->cd()->SetLeftMargin(0.1);
  c_RunQA_Mass2Pi->cd()->SetRightMargin(0.1);
  c_RunQA_Mass2Pi->cd()->SetBottomMargin(0.1);
  c_RunQA_Mass2Pi->cd()->SetGrid(0,0);
  c_RunQA_Mass2Pi->cd()->SetTicks(1,1);

  p_mMass2Pi[1][9]->SetTitle("Mass2Pi vs. runIndex");
  p_mMass2Pi[1][9]->SetStats(0);
  p_mMass2Pi[1][9]->SetMarkerColor(1);
  p_mMass2Pi[1][9]->SetMarkerStyle(20);
  p_mMass2Pi[1][9]->SetMarkerSize(1.0);
  p_mMass2Pi[1][9]->SetLineColor(1);
  p_mMass2Pi[1][9]->GetXaxis()->SetTitle("runIndex");
  p_mMass2Pi[1][9]->GetYaxis()->SetTitle("<Mass2Pi>");
  //pMass2Pi_mGlobPhi[1][9]->GetYaxis()->SetRangeUser(-0.1,0.50);
  p_mMass2Pi[1][9]->Draw("pE");

  for(int i_trig = 0; i_trig < numOfTriggers; ++i_trig)
  {
    p_mMass2Pi[1][i_trig]->SetMarkerColor(MarkerColor[i_trig]);
    p_mMass2Pi[1][i_trig]->SetMarkerStyle(24);
    p_mMass2Pi[1][i_trig]->SetMarkerSize(0.8);
    p_mMass2Pi[1][i_trig]->SetLineColor(MarkerColor[i_trig]);
    if(p_mMass2Pi[1][i_trig]->GetEntries() > 0) p_mMass2Pi[1][i_trig]->Draw("pE same");
  }

  leg->Draw("same");
  c_RunQA_Mass2Pi->SaveAs(Form("./figures/%s/c_RunQA_Mass2Pi.pdf",runQA::mBeamEnergy[energy].c_str()));
  //---------------------
  TCanvas *c_RunQA_Mass2K = new TCanvas("c_RunQA_Mass2K","c_RunQA_dEdx",10,10,800,400);
  c_RunQA_Mass2K->cd()->SetLeftMargin(0.1);
  c_RunQA_Mass2K->cd()->SetRightMargin(0.1);
  c_RunQA_Mass2K->cd()->SetBottomMargin(0.1);
  c_RunQA_Mass2K->cd()->SetGrid(0,0);
  c_RunQA_Mass2K->cd()->SetTicks(1,1);

  p_mMass2K[1][9]->SetTitle("Mass2K vs. runIndex");
  p_mMass2K[1][9]->SetStats(0);
  p_mMass2K[1][9]->SetMarkerColor(1);
  p_mMass2K[1][9]->SetMarkerStyle(20);
  p_mMass2K[1][9]->SetMarkerSize(1.0);
  p_mMass2K[1][9]->SetLineColor(1);
  p_mMass2K[1][9]->GetXaxis()->SetTitle("runIndex");
  p_mMass2K[1][9]->GetYaxis()->SetTitle("<Mass2K>");
  //pMass2K_mGlobPhi[1][9]->GetYaxis()->SetRangeUser(-0.1,0.50);
  p_mMass2K[1][9]->Draw("pE");

  for(int i_trig = 0; i_trig < numOfTriggers; ++i_trig)
  {
    p_mMass2K[1][i_trig]->SetMarkerColor(MarkerColor[i_trig]);
    p_mMass2K[1][i_trig]->SetMarkerStyle(24);
    p_mMass2K[1][i_trig]->SetMarkerSize(0.8);
    p_mMass2K[1][i_trig]->SetLineColor(MarkerColor[i_trig]);
    if(p_mMass2K[1][i_trig]->GetEntries() > 0) p_mMass2K[1][i_trig]->Draw("pE same");
  }

  leg->Draw("same");
  c_RunQA_Mass2K->SaveAs(Form("./figures/%s/c_RunQA_Mass2K.pdf",runQA::mBeamEnergy[energy].c_str()));
  //---------------------
  TCanvas *c_RunQA_Mass2P = new TCanvas("c_RunQA_Mass2P","c_RunQA_dEdx",10,10,800,400);
  c_RunQA_Mass2P->cd()->SetLeftMargin(0.1);
  c_RunQA_Mass2P->cd()->SetRightMargin(0.1);
  c_RunQA_Mass2P->cd()->SetBottomMargin(0.1);
  c_RunQA_Mass2P->cd()->SetGrid(0,0);
  c_RunQA_Mass2P->cd()->SetTicks(1,1);

  p_mMass2P[1][9]->SetTitle("Mass2P vs. runIndex");
  p_mMass2P[1][9]->SetStats(0);
  p_mMass2P[1][9]->SetMarkerColor(1);
  p_mMass2P[1][9]->SetMarkerStyle(20);
  p_mMass2P[1][9]->SetMarkerSize(1.0);
  p_mMass2P[1][9]->SetLineColor(1);
  p_mMass2P[1][9]->GetXaxis()->SetTitle("runIndex");
  p_mMass2P[1][9]->GetYaxis()->SetTitle("<Mass2P>");
  //pMass2P_mGlobPhi[1][9]->GetYaxis()->SetRangeUser(-0.1,0.50);
  p_mMass2P[1][9]->Draw("pE");

  for(int i_trig = 0; i_trig < numOfTriggers; ++i_trig)
  {
    p_mMass2P[1][i_trig]->SetMarkerColor(MarkerColor[i_trig]);
    p_mMass2P[1][i_trig]->SetMarkerStyle(24);
    p_mMass2P[1][i_trig]->SetMarkerSize(0.8);
    p_mMass2P[1][i_trig]->SetLineColor(MarkerColor[i_trig]);
    if(p_mMass2P[1][i_trig]->GetEntries() > 0) p_mMass2P[1][i_trig]->Draw("pE same");
  }

  leg->Draw("same");
  c_RunQA_Mass2P->SaveAs(Form("./figures/%s/c_RunQA_Mass2P.pdf",runQA::mBeamEnergy[energy].c_str()));
  //---------------------
  TCanvas *c_RunQA_Mass2E = new TCanvas("c_RunQA_Mass2E","c_RunQA_dEdx",10,10,800,400);
  c_RunQA_Mass2E->cd()->SetLeftMargin(0.1);
  c_RunQA_Mass2E->cd()->SetRightMargin(0.1);
  c_RunQA_Mass2E->cd()->SetBottomMargin(0.1);
  c_RunQA_Mass2E->cd()->SetGrid(0,0);
  c_RunQA_Mass2E->cd()->SetTicks(1,1);

  p_mMass2E[1][9]->SetTitle("Mass2E vs. runIndex");
  p_mMass2E[1][9]->SetStats(0);
  p_mMass2E[1][9]->SetMarkerColor(1);
  p_mMass2E[1][9]->SetMarkerStyle(20);
  p_mMass2E[1][9]->SetMarkerSize(1.0);
  p_mMass2E[1][9]->SetLineColor(1);
  p_mMass2E[1][9]->GetXaxis()->SetTitle("runIndex");
  p_mMass2E[1][9]->GetYaxis()->SetTitle("<Mass2E>");
  //pMass2E_mGlobPhi[1][9]->GetYaxis()->SetRangeUser(-0.1,0.50);
  p_mMass2E[1][9]->Draw("pE");

  for(int i_trig = 0; i_trig < numOfTriggers; ++i_trig)
  {
    p_mMass2E[1][i_trig]->SetMarkerColor(MarkerColor[i_trig]);
    p_mMass2E[1][i_trig]->SetMarkerStyle(24);
    p_mMass2E[1][i_trig]->SetMarkerSize(0.8);
    p_mMass2E[1][i_trig]->SetLineColor(MarkerColor[i_trig]);
    if(p_mMass2E[1][i_trig]->GetEntries() > 0) p_mMass2E[1][i_trig]->Draw("pE same");
  }

  leg->Draw("same");
  c_RunQA_Mass2E->SaveAs(Form("./figures/%s/c_RunQA_Mass2E.pdf",runQA::mBeamEnergy[energy].c_str()));  
}
