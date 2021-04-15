#include "../StRoot/StEventPlaneMaker/StEventPlaneCons.h"
#include <string>
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TLegend.h"

using namespace std;

void plotChargedFlow(int beamEnergy = 0, int cent9 = 9)
{
  string inputfile = Form("../StRoot/StEventPlaneUtility/ChargedFlow/file_%s_ChargedFlow.root",recoEP::mBeamEnergy[beamEnergy].c_str());
  TFile *File_InPut = TFile::Open(inputfile.c_str());
  TProfile *p_mChargedV1PpQA[10]; // <v1Pp> vs. runIndex | pt [0.2, 2.0]
  TProfile *p_mChargedV2EpQA[10]; // <v2Ep> vs. runIndex | pt [0.2, 5.2]
  TProfile *p_mChargedV2PpQA[10]; // <v2Pp> vs. runIndex | pt [0.2, 5.2]
  TProfile *p_mChargedV3EpQA[10]; // <v3Ep> vs. runIndex | pt [0.2, 5.2]
  TProfile *p_mChargedV1Pp[10]; // v1Pp vs. eta
  TProfile *p_mChargedV2Ep[10]; // v2Ep vs. pt
  TProfile *p_mChargedV2Pp[10]; // v2Pp vs. pt
  TProfile *p_mChargedV3Ep[10]; // v3Ep vs. pt
  for(int i_cent = 0; i_cent < 10; ++i_cent)
  {
    string ProName;
    ProName = Form("p_mChargedV1PpQA_Cent%d",i_cent);
    p_mChargedV1PpQA[i_cent] = (TProfile*)File_InPut->Get(ProName.c_str()); 
    ProName = Form("p_mChargedV2EpQA_Cent%d",i_cent);
    p_mChargedV2EpQA[i_cent] = (TProfile*)File_InPut->Get(ProName.c_str());
    ProName = Form("p_mChargedV2PpQA_Cent%d",i_cent);
    p_mChargedV2PpQA[i_cent] = (TProfile*)File_InPut->Get(ProName.c_str());
    ProName = Form("p_mChargedV3EpQA_Cent%d",i_cent);
    p_mChargedV3EpQA[i_cent] = (TProfile*)File_InPut->Get(ProName.c_str());
    ProName = Form("p_mChargedV1Pp_Cent%d",i_cent);
    p_mChargedV1Pp[i_cent] = (TProfile*)File_InPut->Get(ProName.c_str());
    ProName = Form("p_mChargedV2Ep_Cent%d",i_cent);
    p_mChargedV2Ep[i_cent] = (TProfile*)File_InPut->Get(ProName.c_str());
    ProName = Form("p_mChargedV2Pp_Cent%d",i_cent);
    p_mChargedV2Pp[i_cent] = (TProfile*)File_InPut->Get(ProName.c_str());
    ProName = Form("p_mChargedV3Ep_Cent%d",i_cent);
    p_mChargedV3Ep[i_cent] = (TProfile*)File_InPut->Get(ProName.c_str());
  }

  TCanvas *c_v1 = new TCanvas("c_v1","c_v1",10,10,800,800);
  c_v1->SetLeftMargin(0.15);
  c_v1->SetBottomMargin(0.15);
  c_v1->SetGrid(0,0);
  c_v1->SetTicks(1,1);
  c_v1->cd();

  TH1F *h_play = new TH1F("h_play","h_play",100,-2.0,8.0);
  for(Int_t i_bin = 0; i_bin < 100; i_bin++)
  {
    h_play->SetBinContent(i_bin+1,-10.0);
    h_play->SetBinError(i_bin+1,1.0);
  }
  h_play->SetTitle("");
  h_play->SetStats(0);
  h_play->GetXaxis()->SetTitle("#eta");
  h_play->GetYaxis()->SetTitle("v_{1} (%)");
  h_play->GetXaxis()->CenterTitle();
  h_play->GetYaxis()->CenterTitle();
  h_play->GetXaxis()->SetTitleSize(0.06);
  h_play->GetYaxis()->SetTitleSize(0.06);
  h_play->GetXaxis()->SetRangeUser(-1.1,1.1);
  h_play->GetYaxis()->SetRangeUser(-0.4,0.4);
  h_play->GetXaxis()->SetLabelSize(0.04);
  h_play->GetYaxis()->SetLabelSize(0.04);
  h_play->SetNdivisions(505,"X");
  h_play->SetNdivisions(505,"Y");
  h_play->DrawCopy("pE");
  p_mChargedV1Pp[cent9]->Scale(100.0);
  p_mChargedV1Pp[cent9]->SetMarkerStyle(24);
  p_mChargedV1Pp[cent9]->SetMarkerColor(kGray+2);
  p_mChargedV1Pp[cent9]->SetMarkerSize(1.5);
  p_mChargedV1Pp[cent9]->SetLineColor(kGray+2);
  p_mChargedV1Pp[cent9]->Draw("pE1X0 same");
  string FigureName = Form("./figures/c_mChargedV1Pp_%s.eps",recoEP::mBeamEnergy[beamEnergy].c_str());
  c_v1->SaveAs(FigureName.c_str());

  TCanvas *c_v2 = new TCanvas("c_v2","c_v2",10,10,800,800);
  c_v2->SetLeftMargin(0.15);
  c_v2->SetBottomMargin(0.15);
  c_v2->SetGrid(0,0);
  c_v2->SetTicks(1,1);
  c_v2->cd();

  h_play->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h_play->GetYaxis()->SetTitle("v_{2}");
  h_play->GetXaxis()->SetRangeUser(0.0,5.0);
  h_play->GetYaxis()->SetRangeUser(-0.05,0.35);
  h_play->DrawCopy("pE");

  p_mChargedV2Ep[cent9]->SetMarkerStyle(24);
  p_mChargedV2Ep[cent9]->SetMarkerColor(kAzure+2);
  p_mChargedV2Ep[cent9]->SetMarkerSize(1.5);
  p_mChargedV2Ep[cent9]->SetLineColor(kAzure+2);
  p_mChargedV2Ep[cent9]->Draw("pE1X0 same");

  p_mChargedV2Pp[cent9]->SetMarkerStyle(24);
  p_mChargedV2Pp[cent9]->SetMarkerColor(kGray+2);
  p_mChargedV2Pp[cent9]->SetMarkerSize(1.5);
  p_mChargedV2Pp[cent9]->SetLineColor(kGray+2);
  p_mChargedV2Pp[cent9]->Draw("pE1X0 same");

  TLegend *leg = new TLegend(0.20,0.70,0.45,0.85);
  leg->SetFillColor(10);
  leg->SetBorderSize(0);
  leg->AddEntry(p_mChargedV2Ep[cent9],"v_{2}^{EP,TPC}","p");
  leg->AddEntry(p_mChargedV2Pp[cent9],"v_{2}^{PP,ZDC-SMD}","p");
  leg->Draw("same");
  FigureName = Form("./figures/c_mChargedV2EpPp_%s.eps",recoEP::mBeamEnergy[beamEnergy].c_str());
  c_v2->SaveAs(FigureName.c_str());

  TCanvas *c_v3 = new TCanvas("c_v3","c_v3",10,10,800,800);
  c_v3->SetLeftMargin(0.15);
  c_v3->SetBottomMargin(0.15);
  c_v3->SetGrid(0,0);
  c_v3->SetTicks(1,1);
  c_v3->cd();

  h_play->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h_play->GetYaxis()->SetTitle("v_{3}");
  h_play->GetXaxis()->SetRangeUser(0.0,5.0);
  h_play->GetYaxis()->SetRangeUser(-0.05,0.35);
  h_play->DrawCopy("pE");

  p_mChargedV3Ep[cent9]->SetMarkerStyle(24);
  p_mChargedV3Ep[cent9]->SetMarkerColor(kAzure+2);
  p_mChargedV3Ep[cent9]->SetMarkerSize(1.5);
  p_mChargedV3Ep[cent9]->SetLineColor(kAzure+2);
  p_mChargedV3Ep[cent9]->Draw("pE1X0 same");
  FigureName = Form("./figures/c_mChargedV3Ep_%s.eps",recoEP::mBeamEnergy[beamEnergy].c_str());
  c_v3->SaveAs(FigureName.c_str());

  TCanvas *c_QA = new TCanvas("c_QA","c_QA",10,10,1200,1200);
  c_QA->Divide(1,4);
  for(int i_pad = 0; i_pad < 4; ++i_pad)
  {
    c_QA->cd(i_pad+1)->SetTopMargin(0.05);
    c_QA->cd(i_pad+1)->SetLeftMargin(0.15);
    c_QA->cd(i_pad+1)->SetBottomMargin(0.15);
    c_QA->cd(i_pad+1)->SetGrid(0,0);
    c_QA->cd(i_pad+1)->SetTicks(1,1);
  }

  TH1F *h_QA = new TH1F("h_QA","h_QA",recoEP::mNumOfRunIndex,-0.5,(float)recoEP::mNumOfRunIndex-0.5);
  for(Int_t i_bin = 0; i_bin < recoEP::mNumOfRunIndex; i_bin++)
  {
    h_QA->SetBinContent(i_bin+1,-10.0);
    h_QA->SetBinError(i_bin+1,1.0);
  }
  h_QA->SetTitle("");
  h_QA->SetStats(0);
  h_QA->GetXaxis()->SetTitle("runIndex");
  h_QA->GetXaxis()->CenterTitle();
  h_QA->GetYaxis()->CenterTitle();
  h_QA->GetXaxis()->SetTitleSize(0.06);
  h_QA->GetYaxis()->SetTitleSize(0.06);
  h_QA->GetXaxis()->SetLabelSize(0.04);
  h_QA->GetYaxis()->SetLabelSize(0.04);
  h_QA->SetNdivisions(505,"X");
  h_QA->SetNdivisions(505,"Y");

  c_QA->cd(1); // <v1Pp> vs. runIndex
  h_QA->GetYaxis()->SetTitle("<v_{1}>");
  h_QA->GetYaxis()->SetRangeUser(-0.01,0.01);
  h_QA->DrawCopy("pE");
  // p_mChargedV1PpQA[cent9]->Scale(100.0);
  p_mChargedV1PpQA[cent9]->SetMarkerStyle(24);
  p_mChargedV1PpQA[cent9]->SetMarkerColor(kGray+2);
  p_mChargedV1PpQA[cent9]->SetMarkerSize(0.5);
  p_mChargedV1PpQA[cent9]->SetLineColor(kGray+2);
  p_mChargedV1PpQA[cent9]->Draw("pE1X0 same");

  c_QA->cd(2); // <v2Pp> vs. runIndex
  h_QA->GetYaxis()->SetTitle("<v_{2}^{PP,ZDC-SMD}>");
  h_QA->GetYaxis()->SetRangeUser(-0.01,0.10);
  h_QA->DrawCopy("pE");
  // p_mChargedV1PpQA[cent9]->Scale(100.0);
  p_mChargedV2PpQA[cent9]->SetMarkerStyle(24);
  p_mChargedV2PpQA[cent9]->SetMarkerColor(kGray+2);
  p_mChargedV2PpQA[cent9]->SetMarkerSize(0.5);
  p_mChargedV2PpQA[cent9]->SetLineColor(kGray+2);
  p_mChargedV2PpQA[cent9]->Draw("pE1X0 same");

  c_QA->cd(3); // <v2Ep> vs. runIndex
  h_QA->GetYaxis()->SetTitle("<v_{2}^{EP,ZDC-SMD}>");
  h_QA->GetYaxis()->SetRangeUser(-0.01,0.10);
  h_QA->DrawCopy("pE");
  // p_mChargedV1PpQA[cent9]->Scale(100.0);
  p_mChargedV2EpQA[cent9]->SetMarkerStyle(24);
  p_mChargedV2EpQA[cent9]->SetMarkerColor(kGray+2);
  p_mChargedV2EpQA[cent9]->SetMarkerSize(0.5);
  p_mChargedV2EpQA[cent9]->SetLineColor(kGray+2);
  p_mChargedV2EpQA[cent9]->Draw("pE1X0 same");

  c_QA->cd(4); // <v3Ep> vs. runIndex
  h_QA->GetYaxis()->SetTitle("<v_{3}^{EP,ZDC-SMD}>");
  h_QA->GetYaxis()->SetRangeUser(-0.01,0.10);
  h_QA->DrawCopy("pE");
  // p_mChargedV1PpQA[cent9]->Scale(100.0);
  p_mChargedV3EpQA[cent9]->SetMarkerStyle(24);
  p_mChargedV3EpQA[cent9]->SetMarkerColor(kGray+2);
  p_mChargedV3EpQA[cent9]->SetMarkerSize(0.5);
  p_mChargedV3EpQA[cent9]->SetLineColor(kGray+2);
  p_mChargedV3EpQA[cent9]->Draw("pE1X0 same");
  FigureName = Form("./figures/c_mChargedFlowQA_%s.eps",recoEP::mBeamEnergy[beamEnergy].c_str());
  c_QA->SaveAs(FigureName.c_str());
}
