#include "../StRoot/StEventPlaneMaker/StEventPlaneCons.h"
#include <string>
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"

using namespace std;

double ResolutionFull(double *x_val, double *par)
{
  double y;
  double chi = x_val[0];
  double arg = chi*chi/4.0;
  double norm = TMath::Sqrt(TMath::Pi()/2.0)/2.0;

  y = norm*chi*TMath::Exp(-1.0*arg)*(TMath::BesselI0(arg)+TMath::BesselI1(arg));

  return y;
}

void plotResolution(int beamEnergy = 1)
{
  string inputfile = Form("../StRoot/StEventPlaneUtility/Resolution/file_%s_Resolution.root",recoEP::mBeamEnergy[beamEnergy].c_str());
  TFile *File_InPut = TFile::Open(inputfile.c_str());
  TF1 *f_res = new TF1("f_res",ResolutionFull,0,10,0);
  double Centrality_start[9] = {0.7,0.6,0.5,0.4,0.3,0.2,0.1,0.05, 0.0};
  double Centrality_stop[9]  = {0.8,0.7,0.6,0.5,0.4,0.3,0.2, 0.1,0.05};

  // ZDC-SMD event plane resolution
  //double mZdcFullRes1Val[9];
  //double mZdcFullRes1Err[9];
  //double mZdcFullRes2Val[9];
  //double mZdcFullRes2Err[9];
  //TGraphAsymmErrors *g_mZdcFullRes1 = new TGraphAsymmErrors();
  //TGraphAsymmErrors *g_mZdcFullRes2 = new TGraphAsymmErrors();

  //// calculate 1st full event plane resolution
  //cout << "ZDC-SMD Event Plane Resolution:" << endl;
  //cout << "1st Full Event Plane:" << endl;
  //TProfile *p_mZdcSubRes1 = (TProfile*)File_InPut->Get("p_mZdcSubRes1");
  //for(int i_cent = 0; i_cent < 9; ++i_cent)
  //{
  //  const double resRaw = p_mZdcSubRes1->GetBinContent(p_mZdcSubRes1->FindBin(i_cent));
  //  const double errRaw = p_mZdcSubRes1->GetBinError(p_mZdcSubRes1->FindBin(i_cent));
  //  if(resRaw > 0)
  //  {
  //    const double resSub = TMath::Sqrt(resRaw); // take square root of <cos(psie-psiw)> for sub event plane resolution
  //    const double errSub = errRaw/(2.0*TMath::Sqrt(resRaw)); // calculate the error for this operation

  //    const double chiSub = f_res->GetX(resSub); // use function to retrieve the chi value (x) for this sub event plane resolution (y)
  //    const double errChiSub = errSub/f_res->Derivative(chiSub); // calculate error for this value
  //    const double chiFull = chiSub*TMath::Sqrt(2.0); // the full event plane chi = sqrt(2)*chi_sub
  //    mZdcFullRes1Val[i_cent] = f_res->Eval(chiFull); // evaluate the function for this chi_full value to retrieve the full psi resolution
  //    mZdcFullRes1Err[i_cent] = f_res->Derivative(chiFull)*errChiSub*TMath::Sqrt(2.0); // calculate the error
  //  }
  //  cout << "i_cent = " << i_cent << ", resRaw = " << resRaw << ", resFull = " << mZdcFullRes1Val[i_cent] << " +/  - " << mZdcFullRes1Err[i_cent] << endl;
  //  g_mZdcFullRes1->SetPoint(i_cent,50.0*(Centrality_start[i_cent]+Centrality_stop[i_cent]),mZdcFullRes1Val[i_cent]*100.0);
  //  g_mZdcFullRes1->SetPointError(i_cent,0.0,0.0,mZdcFullRes1Err[i_cent]*100.0,mZdcFullRes1Err[i_cent]*100.0);
  //}

  //// calculate 2nd full event plane resolution
  //cout << "2nd Full Event Plane:" << endl;
  //TProfile *p_mZdcSubRes2 = (TProfile*)File_InPut->Get("p_mZdcSubRes2");
  //for(int i_cent = 0; i_cent < 9; ++i_cent)
  //{
  //  const double resRaw = p_mZdcSubRes2->GetBinContent(p_mZdcSubRes2->FindBin(i_cent));
  //  const double errRaw = p_mZdcSubRes2->GetBinError(p_mZdcSubRes2->FindBin(i_cent));
  //  if(resRaw > 0)
  //  {
  //    const double resSub = TMath::Sqrt(resRaw);
  //    const double errSub = errRaw/(2.0*TMath::Sqrt(resRaw));

  //    const double chiSub = f_res->GetX(resSub);
  //    const double errChiSub = errSub/f_res->Derivative(chiSub);
  //    const double chiFull = chiSub*TMath::Sqrt(2.0);
  //    mZdcFullRes2Val[i_cent] = f_res->Eval(chiFull);
  //    mZdcFullRes2Err[i_cent] = f_res->Derivative(chiFull)*errChiSub*TMath::Sqrt(2.0);
  //  }
  //  cout << "i_cent = " << i_cent << ", resRaw = " << resRaw << ", resFull = " << mZdcFullRes2Val[i_cent] << " +/  - " << mZdcFullRes2Err[i_cent] << endl;
  //  g_mZdcFullRes2->SetPoint(i_cent,50.0*(Centrality_start[i_cent]+Centrality_stop[i_cent]),mZdcFullRes2Val[i_cent]*100.0);
  //  g_mZdcFullRes2->SetPointError(i_cent,0.0,0.0,mZdcFullRes2Err[i_cent]*100.0,mZdcFullRes2Err[i_cent]*100.0);
  //}

  // TPC event plane resolution
  double mTpcSubRes2Val[9];
  double mTpcSubRes2Err[9];
  double mTpcFullRes2Val[9];
  double mTpcFullRes2Err[9];
  TGraphAsymmErrors *g_mTpcSubRes2 = new TGraphAsymmErrors();
  TGraphAsymmErrors *g_mTpcFullRes2 = new TGraphAsymmErrors();

  // calculate sub event plane resolution
  cout << "TPC Event Plane Resolution:" << endl;
  cout << "Sub Event Plane:" << endl;
  TProfile *p_mTpcSubRes2 = (TProfile*)File_InPut->Get("p_mTpcSubRes2");
  for(int i_cent = 0; i_cent < 9; ++i_cent)
  {
    const double resRaw = p_mTpcSubRes2->GetBinContent(p_mTpcSubRes2->FindBin(i_cent));
    const double errRaw = p_mTpcSubRes2->GetBinError(p_mTpcSubRes2->FindBin(i_cent));
    if(resRaw > 0)
    {
      mTpcSubRes2Val[i_cent] = TMath::Sqrt(resRaw);
      mTpcSubRes2Err[i_cent] = errRaw/(2.0*TMath::Sqrt(resRaw));
    }
    cout << "i_cent = " << i_cent << ", resRaw = " << resRaw << ", resSub = " << mTpcSubRes2Val[i_cent] << " +/- " << mTpcSubRes2Err[i_cent] << endl;
    g_mTpcSubRes2->SetPoint(i_cent,50.0*(Centrality_start[i_cent]+Centrality_stop[i_cent]),mTpcSubRes2Val[i_cent]*100.0);
    g_mTpcSubRes2->SetPointError(i_cent,0.0,0.0,mTpcSubRes2Err[i_cent]*100.0,mTpcSubRes2Err[i_cent]*100.0);
  }

  // calculate full event plane resolution
  cout << "Full Event Plane:" << endl;
  TProfile *p_mTpcRanRes2 = (TProfile*)File_InPut->Get("p_mTpcRanRes2");
  for(int i_cent = 0; i_cent < 9; ++i_cent)
  {
    const double resRaw = p_mTpcRanRes2->GetBinContent(p_mTpcRanRes2->FindBin(i_cent));
    const double errRaw = p_mTpcRanRes2->GetBinError(p_mTpcRanRes2->FindBin(i_cent));
    if(resRaw > 0)
    {
      const double resSub = TMath::Sqrt(resRaw);
      const double errSub = errRaw/(2.0*TMath::Sqrt(resRaw));

      const double chiSub = f_res->GetX(resSub);
      const double errChiSub = errSub/f_res->Derivative(chiSub);
      const double chiFull = chiSub*TMath::Sqrt(2.0);
      mTpcFullRes2Val[i_cent] = f_res->Eval(chiFull);
      mTpcFullRes2Err[i_cent] = f_res->Derivative(chiFull)*errChiSub*TMath::Sqrt(2.0);
    }
    cout << "i_cent = " << i_cent << ", resRaw = " << resRaw << ", resFull = " << mTpcFullRes2Val[i_cent] << " +/- " << mTpcFullRes2Err[i_cent] << endl;
    g_mTpcFullRes2->SetPoint(i_cent,50.0*(Centrality_start[i_cent]+Centrality_stop[i_cent]),mTpcFullRes2Val[i_cent]*100.0);
    g_mTpcFullRes2->SetPointError(i_cent,0.0,0.0,mTpcFullRes2Err[i_cent]*100.0,mTpcFullRes2Err[i_cent]*100.0);
  }

  double mTpcSubRes3Val[9];
  double mTpcSubRes3Err[9];
  double mTpcFullRes3Val[9];
  double mTpcFullRes3Err[9];
  TGraphAsymmErrors *g_mTpcSubRes3 = new TGraphAsymmErrors();
  TGraphAsymmErrors *g_mTpcFullRes3 = new TGraphAsymmErrors();

  // calculate sub event plane resolution
  cout << "TPC Event Plane Resolution:" << endl;
  cout << "Sub Event Plane:" << endl;
  TProfile *p_mTpcSubRes3 = (TProfile*)File_InPut->Get("p_mTpcSubRes3");
  for(int i_cent = 0; i_cent < 9; ++i_cent)
  {
    const double resRaw = p_mTpcSubRes3->GetBinContent(p_mTpcSubRes3->FindBin(i_cent));
    const double errRaw = p_mTpcSubRes3->GetBinError(p_mTpcSubRes3->FindBin(i_cent));
    if(resRaw > 0)
    {
      mTpcSubRes3Val[i_cent] = TMath::Sqrt(resRaw);
      mTpcSubRes3Err[i_cent] = errRaw/(2.0*TMath::Sqrt(resRaw));
    }
    cout << "i_cent = " << i_cent << ", resRaw = " << resRaw << ", resSub = " << mTpcSubRes3Val[i_cent] << " +/- " << mTpcSubRes3Err[i_cent] << endl;
    g_mTpcSubRes3->SetPoint(i_cent,50.0*(Centrality_start[i_cent]+Centrality_stop[i_cent]),mTpcSubRes3Val[i_cent]*100.0);
    g_mTpcSubRes3->SetPointError(i_cent,0.0,0.0,mTpcSubRes3Err[i_cent]*100.0,mTpcSubRes3Err[i_cent]*100.0);
  }

  // calculate full event plane resolution
  cout << "Full Event Plane:" << endl;
  TProfile *p_mTpcRanRes3 = (TProfile*)File_InPut->Get("p_mTpcRanRes3");
  for(int i_cent = 0; i_cent < 9; ++i_cent)
  {
    const double resRaw = p_mTpcRanRes3->GetBinContent(p_mTpcRanRes3->FindBin(i_cent));
    const double errRaw = p_mTpcRanRes3->GetBinError(p_mTpcRanRes3->FindBin(i_cent));
    if(resRaw > 0)
    {
      const double resSub = TMath::Sqrt(resRaw);
      const double errSub = errRaw/(2.0*TMath::Sqrt(resRaw));

      const double chiSub = f_res->GetX(resSub);
      const double errChiSub = errSub/f_res->Derivative(chiSub);
      const double chiFull = chiSub*TMath::Sqrt(2.0);
      mTpcFullRes3Val[i_cent] = f_res->Eval(chiFull);
      mTpcFullRes3Err[i_cent] = f_res->Derivative(chiFull)*errChiSub*TMath::Sqrt(2.0);
    }
    cout << "i_cent = " << i_cent << ", resRaw = " << resRaw << ", resFull = " << mTpcFullRes3Val[i_cent] << " +/- " << mTpcFullRes3Err[i_cent] << endl;
    g_mTpcFullRes3->SetPoint(i_cent,50.0*(Centrality_start[i_cent]+Centrality_stop[i_cent]),mTpcFullRes3Val[i_cent]*100.0);
    g_mTpcFullRes3->SetPointError(i_cent,0.0,0.0,mTpcFullRes3Err[i_cent]*100.0,mTpcFullRes3Err[i_cent]*100.0);
  }

  TCanvas *c_play = new TCanvas("c_play","c_play",10,10,800,800);
  c_play->SetLeftMargin(0.15);
  c_play->SetBottomMargin(0.15);
  c_play->SetGrid(0,0);
  c_play->SetTicks(1,1);
  c_play->cd();

  TH1F *h_play = new TH1F("h_play","h_play",100,0,100);
  for(Int_t i_bin = 0; i_bin < 100; i_bin++)
  {
    h_play->SetBinContent(i_bin+1,-10.0);
    h_play->SetBinError(i_bin+1,1.0);
  }
  h_play->SetTitle("");
  h_play->SetStats(0);
  h_play->GetXaxis()->SetTitle("centrality (%)");
  h_play->GetYaxis()->SetTitle("Resolution (%)");
  h_play->GetXaxis()->CenterTitle();
  h_play->GetYaxis()->CenterTitle();
  h_play->GetXaxis()->SetTitleSize(0.06);
  h_play->GetYaxis()->SetTitleSize(0.06);
  h_play->GetXaxis()->SetRangeUser(0,80);
  h_play->GetYaxis()->SetRangeUser(-1.0,95.0);
  h_play->GetXaxis()->SetLabelSize(0.04);
  h_play->GetYaxis()->SetLabelSize(0.04);
  h_play->SetNdivisions(505,"X");
  h_play->SetNdivisions(505,"Y");
  h_play->Draw("pE");

  //g_mZdcFullRes1->SetMarkerStyle(21);
  //g_mZdcFullRes1->SetMarkerColor(kGray+2);
  //g_mZdcFullRes1->SetMarkerSize(1.5);
  //g_mZdcFullRes1->Draw("pE Same");

  //g_mZdcFullRes2->SetMarkerStyle(25);
  //g_mZdcFullRes2->SetMarkerColor(kGray+2);
  //g_mZdcFullRes2->SetMarkerSize(1.5);
  //g_mZdcFullRes2->Draw("pE Same");

  g_mTpcSubRes2->SetMarkerStyle(24);
  g_mTpcSubRes2->SetMarkerColor(kAzure+2);
  g_mTpcSubRes2->SetMarkerSize(1.5);
  g_mTpcSubRes2->Draw("pE Same");

  g_mTpcFullRes2->SetMarkerStyle(20);
  g_mTpcFullRes2->SetMarkerColor(kAzure+2);
  g_mTpcFullRes2->SetMarkerSize(1.5);
  g_mTpcFullRes2->Draw("pE Same");

  g_mTpcSubRes3->SetMarkerStyle(25);
  g_mTpcSubRes3->SetMarkerColor(kGray+2);
  g_mTpcSubRes3->SetMarkerSize(1.5);
  g_mTpcSubRes3->Draw("pE Same");

  g_mTpcFullRes3->SetMarkerStyle(21);
  g_mTpcFullRes3->SetMarkerColor(kGray+2);
  g_mTpcFullRes3->SetMarkerSize(1.5);
  g_mTpcFullRes3->Draw("pE Same");

  TLegend *leg = new TLegend(0.60,0.70,0.85,0.85);
  leg->SetFillColor(10);
  leg->SetBorderSize(0);
  leg->AddEntry(g_mTpcFullRes2,"2^{nd} Full TPC EP","p");
  leg->AddEntry(g_mTpcSubRes2,"2^{nd} #eta_{sub} TPC EP","p");
  leg->AddEntry(g_mTpcFullRes3,"3^{rd} Full TPC EP","p");
  leg->AddEntry(g_mTpcSubRes3,"3^{rd} #eta_{sub} TPC EP","p");//leg->AddEntry(g_mZdcFullRes1,"1^{st} ZDC-SMD Full EP","p");
  //leg->AddEntry(g_mZdcFullRes2,"2^{nd} ZDC-SMD Full EP","p");
  leg->Draw("same");

  string FigureName = Form("./figures/c_mEpResolution_%s.pdf",recoEP::mBeamEnergy[beamEnergy].c_str());
  c_play->SaveAs(FigureName.c_str());
}
