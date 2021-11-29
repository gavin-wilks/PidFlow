#include "../StRoot/StEventPlaneMaker/StEventPlaneCons.h"
#include <string>
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TLegend.h"

using namespace std;

void plotEventPlane(int beamEnergy = 0, int cent9 = 3)
{
  //--------------EPD EP---------------
  // x-axis: runIndex | y-axis: EP
  TH2F *h_mZdcRawEpEast[9]; // raw EP
  TH2F *h_mZdcRawEpWest[9];
  TH2F *h_mZdcRawEpFull[9]; // Qwest-QEast

  TH2F *h_mZdcReCenterEpEast[9]; // recenter EP
  TH2F *h_mZdcReCenterEpWest[9];
  TH2F *h_mZdcReCenterEpFull[9]; // Qwest-QEast

  TH2F *h_mZdcShiftEpEast[9]; // shift EP
  TH2F *h_mZdcShiftEpWest[9];
  TH2F *h_mZdcShiftEpDiff[9]; // Qwest-QEast
  TH2F *h_mZdcShiftEpFull[9]; // Qwest-QEast & shift
  //--------------ZDC EP---------------

  //--------------TPC EP---------------
  // x-axis: runIndex | y-axis: EP
  TH2F *h_mTpcRawEpEast[9]; // raw EP
  TH2F *h_mTpcRawEpWest[9];
  TH2F *h_mTpcRawEpFull[9];

  TH2F *h_mTpcReCenterEpEast[9]; // recenter EP
  TH2F *h_mTpcReCenterEpWest[9];
  TH2F *h_mTpcReCenterEpFull[9];

  TH2F *h_mTpcShiftEpEast[9]; // shift EP
  TH2F *h_mTpcShiftEpWest[9];
  TH2F *h_mTpcShiftEpRanA[9];
  TH2F *h_mTpcShiftEpRanB[9];
  TH2F *h_mTpcShiftEpFull[9];
  //--------------TPC EP---------------

  string inputRaw = Form("../StRoot/StEventPlaneUtility/ReCenterParameter/file_%s_ReCenterParameter.root",recoEP::mBeamEnergy[beamEnergy].c_str());
  TFile *File_InPutRaw = TFile::Open(inputRaw.c_str());
  for(int i_cent = 0; i_cent < 9; ++i_cent)
  {
    string HistName = Form("h_mZdcRawEpEast_%d",i_cent);
    h_mZdcRawEpEast[i_cent] = (TH2F*)File_InPutRaw->Get(HistName.c_str());
    HistName = Form("h_mZdcRawEpWest_%d",i_cent);
    h_mZdcRawEpWest[i_cent] = (TH2F*)File_InPutRaw->Get(HistName.c_str());
    HistName = Form("h_mZdcRawEpFull_%d",i_cent);
    h_mZdcRawEpFull[i_cent] = (TH2F*)File_InPutRaw->Get(HistName.c_str());

    HistName = Form("h_mTpcRawEpEast_%d",i_cent);
    h_mTpcRawEpEast[i_cent] = (TH2F*)File_InPutRaw->Get(HistName.c_str());
    HistName = Form("h_mTpcRawEpWest_%d",i_cent);
    h_mTpcRawEpWest[i_cent] = (TH2F*)File_InPutRaw->Get(HistName.c_str());
    HistName = Form("h_mTpcRawEpFull_%d",i_cent);
    h_mTpcRawEpFull[i_cent] = (TH2F*)File_InPutRaw->Get(HistName.c_str());
  }

  string inputReCenter = Form("../StRoot/StEventPlaneUtility/ShiftParameter/file_%s_ShiftParameter.root",recoEP::mBeamEnergy[beamEnergy].c_str());
  TFile *File_InPutReCenter = TFile::Open(inputReCenter.c_str());
  for(int i_cent = 0; i_cent < 9; ++i_cent)
  {
    string HistName = Form("h_mZdcReCenterEpEast_%d",i_cent);
    h_mZdcReCenterEpEast[i_cent] = (TH2F*)File_InPutReCenter->Get(HistName.c_str());
    HistName = Form("h_mZdcReCenterEpWest_%d",i_cent);
    h_mZdcReCenterEpWest[i_cent] = (TH2F*)File_InPutReCenter->Get(HistName.c_str());
    HistName = Form("h_mZdcReCenterEpFull_%d",i_cent);
    h_mZdcReCenterEpFull[i_cent] = (TH2F*)File_InPutReCenter->Get(HistName.c_str());

    HistName = Form("h_mTpcReCenterEpEast_%d",i_cent);
    h_mTpcReCenterEpEast[i_cent] = (TH2F*)File_InPutReCenter->Get(HistName.c_str());
    HistName = Form("h_mTpcReCenterEpWest_%d",i_cent);
    h_mTpcReCenterEpWest[i_cent] = (TH2F*)File_InPutReCenter->Get(HistName.c_str());
    HistName = Form("h_mTpcReCenterEpFull_%d",i_cent);
    h_mTpcReCenterEpFull[i_cent] = (TH2F*)File_InPutReCenter->Get(HistName.c_str());
  }

  string inputShift = Form("../StRoot/StEventPlaneUtility/Resolution/file_%s_Resolution.root",recoEP::mBeamEnergy[beamEnergy].c_str());
  TFile *File_InPutShift = TFile::Open(inputShift.c_str());
  for(int i_cent = 0; i_cent < 9; ++i_cent)
  {
    string HistName = Form("h_mZdcShiftEpEast_%d",i_cent);
    h_mZdcShiftEpEast[i_cent] = (TH2F*)File_InPutShift->Get(HistName.c_str());
    HistName = Form("h_mZdcShiftEpWest_%d",i_cent);
    h_mZdcShiftEpWest[i_cent] = (TH2F*)File_InPutShift->Get(HistName.c_str());
    HistName = Form("h_mZdcShiftEpFull_%d",i_cent);
    h_mZdcShiftEpFull[i_cent] = (TH2F*)File_InPutShift->Get(HistName.c_str());

    HistName = Form("h_mTpcShiftEpEast_%d",i_cent);
    h_mTpcShiftEpEast[i_cent] = (TH2F*)File_InPutShift->Get(HistName.c_str());
    HistName = Form("h_mTpcShiftEpWest_%d",i_cent);
    h_mTpcShiftEpWest[i_cent] = (TH2F*)File_InPutShift->Get(HistName.c_str());
    HistName = Form("h_mTpcShiftEpFull_%d",i_cent);
    h_mTpcShiftEpFull[i_cent] = (TH2F*)File_InPutShift->Get(HistName.c_str());
  }

  TCanvas *c_Zdc = new TCanvas("c_Zdc","c_Zdc",10,10,900,300);
  c_Zdc->Divide(3,1);
  for(int i_pad = 0; i_pad < 3; ++i_pad)
  {
    c_Zdc->cd(i_pad+1)->SetLeftMargin(0.15);
    c_Zdc->cd(i_pad+1)->SetBottomMargin(0.15);
    c_Zdc->cd(i_pad+1)->SetGrid(0,0);
    c_Zdc->cd(i_pad+1)->SetTicks(1,1);
  }

  c_Zdc->cd(1); // East ZDC-SMD
  TH1F *h_mZdcReCenterEpEast1D = (TH1F*)h_mZdcReCenterEpEast[cent9]->ProjectionY()->Clone("h_mZdcReCenterEpEast1D");
  h_mZdcReCenterEpEast1D->SetStats(0);
  h_mZdcReCenterEpEast1D->SetLineColor(4);
  h_mZdcReCenterEpEast1D->GetXaxis()->SetTitle("#Psi_{1}^{ZDC-SMD}");
  h_mZdcReCenterEpEast1D->SetNdivisions(505,"X");
  h_mZdcReCenterEpEast1D->RebinX(4);
  h_mZdcReCenterEpEast1D->GetYaxis()->SetTitle("# Events");
  h_mZdcReCenterEpEast1D->GetYaxis()->SetTitleOffset(1.6);
  h_mZdcReCenterEpEast1D->SetNdivisions(505,"Y");
  h_mZdcReCenterEpEast1D->GetYaxis()->SetRangeUser(0.0,h_mZdcReCenterEpEast1D->GetMaximum()*1.2);
  h_mZdcReCenterEpEast1D->Draw("hE");

  TH1F *h_mZdcShiftEpEast1D = (TH1F*)h_mZdcShiftEpEast[cent9]->ProjectionY()->Clone("h_mZdcShiftEpEast1D");
  h_mZdcShiftEpEast1D->SetLineColor(2);
  h_mZdcShiftEpEast1D->RebinX(4);
  h_mZdcShiftEpEast1D->Draw("hE same");

  TH1F *h_mZdcRawEpEast1D = (TH1F*)h_mZdcRawEpEast[cent9]->ProjectionY()->Clone("h_mZdcRawEpEast1D");
  h_mZdcRawEpEast1D->RebinX(4);
  h_mZdcRawEpEast1D->Draw("hE same");

  TLegend *leg = new TLegend(0.2,0.2,0.5,0.5);
  leg->SetFillColor(10);
  leg->SetBorderSize(0);
  leg->AddEntry(h_mZdcRawEpEast1D,"Raw EP","l");
  leg->AddEntry(h_mZdcReCenterEpEast1D,"ReCenter EP","l");
  leg->AddEntry(h_mZdcShiftEpEast1D,"Shift EP","l");
  leg->Draw("same");

  c_Zdc->cd(2); // West ZDC-SMD
  TH1F *h_mZdcReCenterEpWest1D = (TH1F*)h_mZdcReCenterEpWest[cent9]->ProjectionY()->Clone("h_mZdcReCenterEpWest1D");
  h_mZdcReCenterEpWest1D->SetStats(0);
  h_mZdcReCenterEpWest1D->SetLineColor(4);
  h_mZdcReCenterEpWest1D->GetXaxis()->SetTitle("#Psi_{1}^{ZDC-SMD}");
  h_mZdcReCenterEpWest1D->SetNdivisions(505,"X");
  h_mZdcReCenterEpWest1D->RebinX(4);
  h_mZdcReCenterEpWest1D->GetYaxis()->SetTitle("# Events");
  h_mZdcReCenterEpWest1D->GetYaxis()->SetTitleOffset(1.6);
  h_mZdcReCenterEpWest1D->SetNdivisions(505,"Y");
  h_mZdcReCenterEpWest1D->GetYaxis()->SetRangeUser(0.0,h_mZdcReCenterEpWest1D->GetMaximum()*1.2);
  h_mZdcReCenterEpWest1D->Draw("hE");

  TH1F *h_mZdcShiftEpWest1D = (TH1F*)h_mZdcShiftEpWest[cent9]->ProjectionY()->Clone("h_mZdcShiftEpWest1D");
  h_mZdcShiftEpWest1D->SetLineColor(2);
  h_mZdcShiftEpWest1D->RebinX(4);
  h_mZdcShiftEpWest1D->Draw("hE same");

  TH1F *h_mZdcRawEpWest1D = (TH1F*)h_mZdcRawEpWest[cent9]->ProjectionY()->Clone("h_mZdcRawEpWest1D");
  h_mZdcRawEpWest1D->RebinX(4);
  h_mZdcRawEpWest1D->Draw("hE same");

  c_Zdc->cd(3); // Full ZDC-SMD
  TH1F *h_mZdcReCenterEpFull1D = (TH1F*)h_mZdcReCenterEpFull[cent9]->ProjectionY()->Clone("h_mZdcReCenterEpFull1D");
  h_mZdcReCenterEpFull1D->SetStats(0);
  h_mZdcReCenterEpFull1D->SetLineColor(4);
  h_mZdcReCenterEpFull1D->GetXaxis()->SetTitle("#Psi_{1}^{ZDC-SMD}");
  h_mZdcReCenterEpFull1D->SetNdivisions(505,"X");
  h_mZdcReCenterEpFull1D->RebinX(4);
  h_mZdcReCenterEpFull1D->GetYaxis()->SetTitle("# Events");
  h_mZdcReCenterEpFull1D->GetYaxis()->SetTitleOffset(1.6);
  h_mZdcReCenterEpFull1D->SetNdivisions(505,"Y");
  h_mZdcReCenterEpFull1D->GetYaxis()->SetRangeUser(0.0,h_mZdcReCenterEpFull1D->GetMaximum()*1.2);
  h_mZdcReCenterEpFull1D->Draw("hE");

  TH1F *h_mZdcShiftEpFull1D = (TH1F*)h_mZdcShiftEpFull[cent9]->ProjectionY()->Clone("h_mZdcShiftEpFull1D");
  h_mZdcShiftEpFull1D->SetLineColor(2);
  h_mZdcShiftEpFull1D->RebinX(4);
  h_mZdcShiftEpFull1D->Draw("hE same");

  TH1F *h_mZdcRawEpFull1D = (TH1F*)h_mZdcRawEpFull[cent9]->ProjectionY()->Clone("h_mZdcRawEpFull1D");
  h_mZdcRawEpFull1D->RebinX(4);
  h_mZdcRawEpFull1D->Draw("hE same");

  string FigureName = Form("./figures/c_mZdcEventPlane_%s.eps",recoEP::mBeamEnergy[beamEnergy].c_str());
  c_Zdc->SaveAs(FigureName.c_str());

  TCanvas *c_Tpc = new TCanvas("c_Tpc","c_Tpc",10,10,900,300);
  c_Tpc->Divide(3,1);
  for(int i_pad = 0; i_pad < 3; ++i_pad)
  {
    c_Tpc->cd(i_pad+1)->SetLeftMargin(0.15);
    c_Tpc->cd(i_pad+1)->SetBottomMargin(0.15);
    c_Tpc->cd(i_pad+1)->SetGrid(0,0);
    c_Tpc->cd(i_pad+1)->SetTicks(1,1);
  }

  c_Tpc->cd(1); // East TPC
  TH1F *h_mTpcReCenterEpEast1D = (TH1F*)h_mTpcReCenterEpEast[cent9]->ProjectionY()->Clone("h_mTpcReCenterEpEast1D");
  h_mTpcReCenterEpEast1D->SetStats(0);
  h_mTpcReCenterEpEast1D->SetLineColor(4);
  h_mTpcReCenterEpEast1D->GetXaxis()->SetTitle("#Psi_{2}^{TPC}");
  h_mTpcReCenterEpEast1D->SetNdivisions(505,"X");
  h_mTpcReCenterEpEast1D->RebinX(4);
  h_mTpcReCenterEpEast1D->GetXaxis()->SetRangeUser(-0.5*TMath::Pi(),0.5*TMath::Pi());
  h_mTpcReCenterEpEast1D->GetYaxis()->SetTitle("# Events");
  h_mTpcReCenterEpEast1D->GetYaxis()->SetTitleOffset(1.6);
  h_mTpcReCenterEpEast1D->SetNdivisions(505,"Y");
  h_mTpcReCenterEpEast1D->GetYaxis()->SetRangeUser(0.0,h_mTpcReCenterEpEast1D->GetMaximum()*1.2);
  h_mTpcReCenterEpEast1D->Draw("hE");

  TH1F *h_mTpcShiftEpEast1D = (TH1F*)h_mTpcShiftEpEast[cent9]->ProjectionY()->Clone("h_mTpcShiftEpEast1D");
  h_mTpcShiftEpEast1D->SetLineColor(2);
  h_mTpcShiftEpEast1D->RebinX(4);
  h_mTpcShiftEpEast1D->Draw("hE same");

  TH1F *h_mTpcRawEpEast1D = (TH1F*)h_mTpcRawEpEast[cent9]->ProjectionY()->Clone("h_mTpcRawEpEast1D");
  h_mTpcRawEpEast1D->RebinX(4);
  h_mTpcRawEpEast1D->Draw("hE same");

  leg->Draw("same");

  c_Tpc->cd(2); // West TPC
  TH1F *h_mTpcReCenterEpWest1D = (TH1F*)h_mTpcReCenterEpWest[cent9]->ProjectionY()->Clone("h_mTpcReCenterEpWest1D");
  h_mTpcReCenterEpWest1D->SetStats(0);
  h_mTpcReCenterEpWest1D->SetLineColor(4);
  h_mTpcReCenterEpWest1D->GetXaxis()->SetTitle("#Psi_{2}^{TPC}");
  h_mTpcReCenterEpWest1D->SetNdivisions(505,"X");
  h_mTpcReCenterEpWest1D->RebinX(4);
  h_mTpcReCenterEpWest1D->GetXaxis()->SetRangeUser(-0.5*TMath::Pi(),0.5*TMath::Pi());
  h_mTpcReCenterEpWest1D->GetYaxis()->SetTitle("# Events");
  h_mTpcReCenterEpWest1D->GetYaxis()->SetTitleOffset(1.6);
  h_mTpcReCenterEpWest1D->SetNdivisions(505,"Y");
  h_mTpcReCenterEpWest1D->GetYaxis()->SetRangeUser(0.0,h_mTpcReCenterEpWest1D->GetMaximum()*1.2);
  h_mTpcReCenterEpWest1D->Draw("hE");

  TH1F *h_mTpcShiftEpWest1D = (TH1F*)h_mTpcShiftEpWest[cent9]->ProjectionY()->Clone("h_mTpcShiftEpWest1D");
  h_mTpcShiftEpWest1D->SetLineColor(2);
  h_mTpcShiftEpWest1D->RebinX(4);
  h_mTpcShiftEpWest1D->Draw("hE same");

  TH1F *h_mTpcRawEpWest1D = (TH1F*)h_mTpcRawEpWest[cent9]->ProjectionY()->Clone("h_mTpcRawEpWest1D");
  h_mTpcRawEpWest1D->RebinX(4);
  h_mTpcRawEpWest1D->Draw("hE same");

  c_Tpc->cd(3); // Full TPC
  TH1F *h_mTpcReCenterEpFull1D = (TH1F*)h_mTpcReCenterEpFull[cent9]->ProjectionY()->Clone("h_mTpcReCenterEpFull1D");
  h_mTpcReCenterEpFull1D->SetStats(0);
  h_mTpcReCenterEpFull1D->SetLineColor(4);
  h_mTpcReCenterEpFull1D->GetXaxis()->SetTitle("#Psi_{2}^{TPC}");
  h_mTpcReCenterEpFull1D->SetNdivisions(505,"X");
  h_mTpcReCenterEpFull1D->RebinX(4);
  h_mTpcReCenterEpFull1D->GetXaxis()->SetRangeUser(-0.5*TMath::Pi(),0.5*TMath::Pi());
  h_mTpcReCenterEpFull1D->GetYaxis()->SetTitle("# Events");
  h_mTpcReCenterEpFull1D->GetYaxis()->SetTitleOffset(1.6);
  h_mTpcReCenterEpFull1D->SetNdivisions(505,"Y");
  h_mTpcReCenterEpFull1D->GetYaxis()->SetRangeUser(0.0,h_mTpcReCenterEpFull1D->GetMaximum()*1.2);
  h_mTpcReCenterEpFull1D->Draw("hE");

  TH1F *h_mTpcShiftEpFull1D = (TH1F*)h_mTpcShiftEpFull[cent9]->ProjectionY()->Clone("h_mTpcShiftEpFull1D");
  h_mTpcShiftEpFull1D->SetLineColor(2);
  h_mTpcShiftEpFull1D->RebinX(4);
  h_mTpcShiftEpFull1D->Draw("hE same");

  TH1F *h_mTpcRawEpFull1D = (TH1F*)h_mTpcRawEpFull[cent9]->ProjectionY()->Clone("h_mTpcRawEpFull1D");
  h_mTpcRawEpFull1D->RebinX(4);
  h_mTpcRawEpFull1D->Draw("hE same");

  string FigureName = Form("./figures/c_mTpcEventPlane_%s.eps",recoEP::mBeamEnergy[beamEnergy].c_str());
  c_Tpc->SaveAs(FigureName.c_str());
}
