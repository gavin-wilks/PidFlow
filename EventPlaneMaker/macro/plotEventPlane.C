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

void plotEventPlane(int beamEnergy = 1)
{

  //--------------TPC EP---------------
  // x-axis: runIndex | y-axis: EP
  TH2F *h_mTpcRawEpEast[3][9]; // raw EP
  TH2F *h_mTpcRawEpWest[3][9];
  TH2F *h_mTpcRawEpFull[3][9];

  TH2F *h_mTpcReCenterEpEast[3][9]; // recenter EP
  TH2F *h_mTpcReCenterEpWest[3][9];
  TH2F *h_mTpcReCenterEpFull[3][9];

  TH2F *h_mTpcShiftEpEast[3][9]; // shift EP
  TH2F *h_mTpcShiftEpWest[3][9];
  TH2F *h_mTpcShiftEpRanA[3][9];
  TH2F *h_mTpcShiftEpRanB[3][9];
  TH2F *h_mTpcShiftEpFull[3][9];
  //--------------TPC EP---------------

  string inputRaw = Form("../StRoot/StEventPlaneUtility/ReCenterParameter/file_%s_ReCenterParameter.root",recoEP::mBeamEnergy[beamEnergy].c_str());
  TFile *File_InPutRaw = TFile::Open(inputRaw.c_str());
  
  string HistName;

  for(int order = 1; order <= 3; ++order)
  {
    for(int i_cent = 0; i_cent < 9; ++i_cent)
    {
      HistName = Form("h_mTpcRawEpEast_%d_%d",order,i_cent);
      h_mTpcRawEpEast[order-1][i_cent] = (TH2F*)File_InPutRaw->Get(HistName.c_str());
      HistName = Form("h_mTpcRawEpWest_%d_%d",order,i_cent);
      h_mTpcRawEpWest[order-1][i_cent] = (TH2F*)File_InPutRaw->Get(HistName.c_str());
      HistName = Form("h_mTpcRawEpFull_%d_%d",order,i_cent);
      h_mTpcRawEpFull[order-1][i_cent] = (TH2F*)File_InPutRaw->Get(HistName.c_str());
    }
  }

  string inputReCenter = Form("../StRoot/StEventPlaneUtility/ShiftParameter/file_%s_ShiftParameter.root",recoEP::mBeamEnergy[beamEnergy].c_str());
  TFile *File_InPutReCenter = TFile::Open(inputReCenter.c_str());
  for(int order = 1; order <= 3; ++order)
  {
    for(int i_cent = 0; i_cent < 9; ++i_cent)
    {
      HistName = Form("h_mTpcReCenterEpEast_%d_%d",order,i_cent);
      h_mTpcReCenterEpEast[order-1][i_cent] = (TH2F*)File_InPutReCenter->Get(HistName.c_str());
      HistName = Form("h_mTpcReCenterEpWest_%d_%d",order,i_cent);
      h_mTpcReCenterEpWest[order-1][i_cent] = (TH2F*)File_InPutReCenter->Get(HistName.c_str());
      HistName = Form("h_mTpcReCenterEpFull_%d_%d",order,i_cent);
      h_mTpcReCenterEpFull[order-1][i_cent] = (TH2F*)File_InPutReCenter->Get(HistName.c_str());
    }
  }

  string inputShift = Form("../StRoot/StEventPlaneUtility/Resolution/file_%s_Resolution.root",recoEP::mBeamEnergy[beamEnergy].c_str());
  TFile *File_InPutShift = TFile::Open(inputShift.c_str());
  
  for(int order = 1; order <= 3; ++order)
  {
    for(int i_cent = 0; i_cent < 9; ++i_cent)
    {
      HistName = Form("h_mTpcShiftEpEast_%d_%d",order,i_cent);
      h_mTpcShiftEpEast[order-1][i_cent] = (TH2F*)File_InPutShift->Get(HistName.c_str());
      HistName = Form("h_mTpcShiftEpWest_%d_%d",order,i_cent);
      h_mTpcShiftEpWest[order-1][i_cent] = (TH2F*)File_InPutShift->Get(HistName.c_str());
      HistName = Form("h_mTpcShiftEpFull_%d_%d",order,i_cent);
      h_mTpcShiftEpFull[order-1][i_cent] = (TH2F*)File_InPutShift->Get(HistName.c_str());
    }
  }

  TCanvas *c_Tpc = new TCanvas("c_Tpc","c_Tpc",10,10,900,900);
  c_Tpc->Divide(3,3);
  for(int i_pad = 0; i_pad < 9; ++i_pad)
  {
    c_Tpc->cd(i_pad+1)->SetLeftMargin(0.15);
    c_Tpc->cd(i_pad+1)->SetBottomMargin(0.15);
    c_Tpc->cd(i_pad+1)->SetGrid(0,0);
    c_Tpc->cd(i_pad+1)->SetTicks(1,1);
  }
 
  TH1F *h_mTpcReCenterEpEast1D[3][9];
  TH1F *h_mTpcReCenterEpWest1D[3][9];
  TH1F *h_mTpcShiftEpEast1D[3][9];
  TH1F *h_mTpcShiftEpWest1D[3][9];
  TH1F *h_mTpcRawEpEast1D[3][9];
  TH1F *h_mTpcRawEpWest1D[3][9];
  TH1F *h_mTpcReCenterEpFull1D[3][9];
  TH1F *h_mTpcShiftEpFull1D[3][9];
  TH1F *h_mTpcRawEpFull1D[3][9];
  
  TLegend *leg[3][9];

  gStyle->SetOptStat("ou");

  for(int i_cent = 0; i_cent < 9; ++i_cent)
  { 
    for(int order = 1; order <= 3; ++order)
    {
      c_Tpc->cd(1+3*(order-1)); // East TPC
      h_mTpcReCenterEpEast1D[order-1][i_cent] = (TH1F*)h_mTpcReCenterEpEast[order-1][i_cent]->ProjectionY()->Clone(Form("h_mTpcReCenterEpEast1D_cent%d_order%d",i_cent,order));
      //h_mTpcReCenterEpEast1D[order-1][i_cent]->SetStats(0);
      h_mTpcReCenterEpEast1D[order-1][i_cent]->SetLineColor(4);
      h_mTpcReCenterEpEast1D[order-1][i_cent]->GetXaxis()->SetTitle(Form("#Psi_{%d}^{TPC} Cent%d",order,i_cent));
      h_mTpcReCenterEpEast1D[order-1][i_cent]->SetNdivisions(505,"X");
      h_mTpcReCenterEpEast1D[order-1][i_cent]->RebinX(4);
      h_mTpcReCenterEpEast1D[order-1][i_cent]->GetXaxis()->SetRangeUser(-1.0/double(order)*TMath::Pi(),1.0/double(order)*TMath::Pi());
      h_mTpcReCenterEpEast1D[order-1][i_cent]->GetYaxis()->SetTitle("# Events");
      h_mTpcReCenterEpEast1D[order-1][i_cent]->GetYaxis()->SetTitleOffset(1.6);
      h_mTpcReCenterEpEast1D[order-1][i_cent]->SetNdivisions(505,"Y");
      h_mTpcReCenterEpEast1D[order-1][i_cent]->GetYaxis()->SetRangeUser(0.0,h_mTpcReCenterEpEast1D[order-1][i_cent]->GetMaximum()*1.2);
      h_mTpcReCenterEpEast1D[order-1][i_cent]->Draw("hE");

      h_mTpcShiftEpEast1D[order-1][i_cent] = (TH1F*)h_mTpcShiftEpEast[order-1][i_cent]->ProjectionY()->Clone(Form("h_mTpcShiftEpEast1D_cent%d_order%d",i_cent,order));
      h_mTpcShiftEpEast1D[order-1][i_cent]->SetLineColor(2);
      h_mTpcShiftEpEast1D[order-1][i_cent]->RebinX(4);
      h_mTpcShiftEpEast1D[order-1][i_cent]->Draw("hE same");

      h_mTpcRawEpEast1D[order-1][i_cent] = (TH1F*)h_mTpcRawEpEast[order-1][i_cent]->ProjectionY()->Clone(Form("h_mTpcRawEpEast1D_cent%d_order%d",i_cent,order));
      h_mTpcRawEpEast1D[order-1][i_cent]->RebinX(4);
      h_mTpcRawEpEast1D[order-1][i_cent]->Draw("hE same");

      leg[0][i_cent] = new TLegend(0.2,0.2,0.5,0.5);
      leg[0][i_cent]->SetFillColor(10);
      leg[0][i_cent]->SetBorderSize(0);
      leg[0][i_cent]->AddEntry(h_mTpcRawEpEast1D[order-1][i_cent],"E Raw EP","l");
      leg[0][i_cent]->AddEntry(h_mTpcReCenterEpEast1D[order-1][i_cent],"E ReCenter EP","l");
      leg[0][i_cent]->AddEntry(h_mTpcShiftEpEast1D[order-1][i_cent],"E Shift EP","l"); 
      leg[0][i_cent]->Draw("same");

      cout << "EAST    Order: " << order << "    Centrality: " << i_cent << endl;


      c_Tpc->cd(2+3*(order-1)); // West TPC
      h_mTpcReCenterEpWest1D[order-1][i_cent] = (TH1F*)h_mTpcReCenterEpWest[order-1][i_cent]->ProjectionY()->Clone(Form("h_mTpcReCenterEpWest1D_cent%d_order%d",i_cent,order));
      //h_mTpcReCenterEpWest1D[order-1][i_cent]->SetStats(0);
      h_mTpcReCenterEpWest1D[order-1][i_cent]->SetLineColor(4);
      h_mTpcReCenterEpWest1D[order-1][i_cent]->GetXaxis()->SetTitle(Form("#Psi_{%d}^{TPC} Cent%d",order,i_cent));
      h_mTpcReCenterEpWest1D[order-1][i_cent]->SetNdivisions(505,"X");
      h_mTpcReCenterEpWest1D[order-1][i_cent]->RebinX(4);
      h_mTpcReCenterEpWest1D[order-1][i_cent]->GetXaxis()->SetRangeUser(-1.0/double(order)*TMath::Pi(),1.0/double(order)*TMath::Pi());
      h_mTpcReCenterEpWest1D[order-1][i_cent]->GetYaxis()->SetTitle("# Events");
      h_mTpcReCenterEpWest1D[order-1][i_cent]->GetYaxis()->SetTitleOffset(1.6);
      h_mTpcReCenterEpWest1D[order-1][i_cent]->SetNdivisions(505,"Y");
      h_mTpcReCenterEpWest1D[order-1][i_cent]->GetYaxis()->SetRangeUser(0.0,h_mTpcReCenterEpWest1D[order-1][i_cent]->GetMaximum()*1.2);
      h_mTpcReCenterEpWest1D[order-1][i_cent]->Draw("hE");

      h_mTpcShiftEpWest1D[order-1][i_cent] = (TH1F*)h_mTpcShiftEpWest[order-1][i_cent]->ProjectionY()->Clone(Form("h_mTpcShiftEpWest1D_cent%d_order%d",i_cent,order));
      h_mTpcShiftEpWest1D[order-1][i_cent]->SetLineColor(2);
      h_mTpcShiftEpWest1D[order-1][i_cent]->RebinX(4);
      h_mTpcShiftEpWest1D[order-1][i_cent]->Draw("hE same");

      h_mTpcRawEpWest1D[order-1][i_cent] = (TH1F*)h_mTpcRawEpWest[order-1][i_cent]->ProjectionY()->Clone(Form("h_mTpcRawEpWest1D_cent%d_order%d",i_cent,order));
      h_mTpcRawEpWest1D[order-1][i_cent]->RebinX(4);
      h_mTpcRawEpWest1D[order-1][i_cent]->Draw("hE same");

      leg[1][i_cent] = new TLegend(0.2,0.2,0.5,0.5);
      leg[1][i_cent]->SetFillColor(10);
      leg[1][i_cent]->SetBorderSize(0);
      leg[1][i_cent]->AddEntry(h_mTpcRawEpWest1D[order-1][i_cent],"W Raw EP","l");
      leg[1][i_cent]->AddEntry(h_mTpcReCenterEpWest1D[order-1][i_cent],"W ReCenter EP","l");
      leg[1][i_cent]->AddEntry(h_mTpcShiftEpWest1D[order-1][i_cent],"W Shift EP","l"); 
      leg[1][i_cent]->Draw("same");

      cout << "WEST   Order: " << order << "    Centrality: " << i_cent << endl;

 
      c_Tpc->cd(3+3*(order-1)); // Full TPC
      h_mTpcReCenterEpFull1D[order-1][i_cent] = (TH1F*)h_mTpcReCenterEpFull[order-1][i_cent]->ProjectionY()->Clone(Form("h_mTpcReCenterEpFull1D_cent%d_order%d",i_cent,order));
      //h_mTpcReCenterEpFull1D[order-1][i_cent]->SetStats(0);
      h_mTpcReCenterEpFull1D[order-1][i_cent]->SetLineColor(4);
      h_mTpcReCenterEpFull1D[order-1][i_cent]->GetXaxis()->SetTitle(Form("#Psi_{%d}^{TPC} Cent%d",order,i_cent));
      h_mTpcReCenterEpFull1D[order-1][i_cent]->SetNdivisions(505,"X");
      h_mTpcReCenterEpFull1D[order-1][i_cent]->RebinX(4);
      h_mTpcReCenterEpFull1D[order-1][i_cent]->GetXaxis()->SetRangeUser(-1.0/double(order)*TMath::Pi(),1.0/double(order)*TMath::Pi());
      h_mTpcReCenterEpFull1D[order-1][i_cent]->GetYaxis()->SetTitle("# Events");
      h_mTpcReCenterEpFull1D[order-1][i_cent]->GetYaxis()->SetTitleOffset(1.6);
      h_mTpcReCenterEpFull1D[order-1][i_cent]->SetNdivisions(505,"Y");
      h_mTpcReCenterEpFull1D[order-1][i_cent]->GetYaxis()->SetRangeUser(0.0,h_mTpcReCenterEpFull1D[order-1][i_cent]->GetMaximum()*1.2);
      h_mTpcReCenterEpFull1D[order-1][i_cent]->Draw("hE");

      h_mTpcShiftEpFull1D[order-1][i_cent] = (TH1F*)h_mTpcShiftEpFull[order-1][i_cent]->ProjectionY()->Clone(Form("h_mTpcShiftEpFull1D_cent%d_order%d",i_cent,order));
      h_mTpcShiftEpFull1D[order-1][i_cent]->SetLineColor(2);
      h_mTpcShiftEpFull1D[order-1][i_cent]->RebinX(4);
      h_mTpcShiftEpFull1D[order-1][i_cent]->Draw("hE same");
      h_mTpcRawEpFull1D[order-1][i_cent] = (TH1F*)h_mTpcRawEpFull[order-1][i_cent]->ProjectionY()->Clone(Form("h_mTpcRawEpFull1D_cent%d_order%d",i_cent,order));
      h_mTpcRawEpFull1D[order-1][i_cent]->RebinX(4);
      h_mTpcRawEpFull1D[order-1][i_cent]->Draw("hE same");

      leg[2][i_cent] = new TLegend(0.2,0.2,0.5,0.5);
      leg[2][i_cent]->SetFillColor(10);
      leg[2][i_cent]->SetBorderSize(0);
      leg[2][i_cent]->AddEntry(h_mTpcRawEpFull1D[order-1][i_cent],"F Raw EP","l");
      leg[2][i_cent]->AddEntry(h_mTpcReCenterEpFull1D[order-1][i_cent],"F ReCenter EP","l");
      leg[2][i_cent]->AddEntry(h_mTpcShiftEpFull1D[order-1][i_cent],"F Shift EP","l"); 
      leg[2][i_cent]->Draw("same");

    }
    string FigureName = Form("./figures/c_mTpcEventPlane_%s_Cent%d.pdf",recoEP::mBeamEnergy[beamEnergy].c_str(),i_cent);
    c_Tpc->SaveAs(FigureName.c_str());
  }
}
