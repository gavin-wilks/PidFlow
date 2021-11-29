#include "../StRoot/StEventPlaneMaker/StEventPlaneCons.h"
#include <string>
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TGraphAsymmErrors.h"
#include "TProfile.h"
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

void plotEpdEpResolution(int beamEnergy = 1, int EpdMode = 2, string JobId = "1")
{
  string corrName[3] = {"Raw", "Pw", "PwS"}; // Raw, phi-weighted, phi-weighted and shifted
 
  TProfile *p_mEpdSubRes[recoEP::mEpdEpOrder][3]; 
  TProfile *p_mD12[3];
  TProfile *p_mD12R1[3];

  TGraphAsymmErrors *g_mEpdFullRes[recoEP::mEpdEpOrder][3];
 
  const int nCent = 9;
  double Centrality_start[9] = {0.7,0.6,0.5,0.4,0.3,0.2,0.1,0.05, 0.0};
  double Centrality_stop[9]  = {0.8,0.7,0.6,0.5,0.4,0.3,0.2, 0.1,0.05};  

  double mEpdFullResVal[recoEP::mEpdEpOrder][3][nCent];
  double mEpdFullResErr[recoEP::mEpdEpOrder][3][nCent];
 
  int color[4] = {1,2,4,3};
  int mstyle[4] = {21,25,24,20};
  string inputfile = Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/AuAu%s/OutPut/PidFlow/file_%s_EpdResults_%d_%s.root",recoEP::mBeamEnergy[beamEnergy],recoEP::mBeamEnergy[beamEnergy],EpdMode,JobId.c_str());
  TFile *File_InPut = TFile::Open(inputfile.c_str());

  TF1 *f_res = new TF1("f_res",ResolutionFull,0,10,0);

  for(int i_corr = 0; i_corr < 3; ++i_corr)
  {
    for(int order = 1; order <= recoEP::mEpdEpOrder; order++)
    {
      string HistName = Form("EpdCos%d%s",order,corrName[i_corr].c_str());
      p_mEpdSubRes[order-1][i_corr] = (TProfile*) File_InPut->Get(HistName.c_str());
      g_mEpdFullRes[order-1][i_corr] = new TGraphAsymmErrors();
      for(int i_cent = 0; i_cent < nCent; ++i_cent)
      {
        const double resRaw = p_mEpdSubRes[order-1][i_corr]->GetBinContent(p_mEpdSubRes[order-1][i_corr]->FindBin(i_cent));
        const double errRaw = p_mEpdSubRes[order-1][i_corr]->GetBinError(p_mEpdSubRes[order-1][i_corr]->FindBin(i_cent));
        if(resRaw > 0)
        {
          const double resSub = TMath::Sqrt(resRaw);
          const double errSub = errRaw/(2.0*TMath::Sqrt(resRaw));

          const double chiSub = f_res->GetX(resSub);
          const double errChiSub = errSub/f_res->Derivative(chiSub);
          const double chiFull = chiSub*TMath::Sqrt(2.0);
          mEpdFullResVal[order-1][i_corr][i_cent] = f_res->Eval(chiFull);
          mEpdFullResErr[order-1][i_corr][i_cent] = f_res->Derivative(chiFull)*errChiSub*TMath::Sqrt(2.0);
        }
        cout << "i_cent = " << i_cent << ", resRaw = " << resRaw << ", resFull = " << mEpdFullResVal[order-1][i_corr][i_cent] << " +/  - " << mEpdFullResErr[order-1][i_corr][i_cent] << endl;
        cout << "order-1: " << order-1 << "  | i_corr: " << i_corr << "  | i_cent: " << i_cent << endl; 
        g_mEpdFullRes[order-1][i_corr]->SetPoint(i_cent,50.0*(Centrality_start[i_cent]+Centrality_stop[i_cent]),mEpdFullResVal[order-1][i_corr][i_cent]*100.0);
        g_mEpdFullRes[order-1][i_corr]->SetPointError(i_cent,0.0,0.0,mEpdFullResErr[order-1][i_corr][i_cent]*100.0,mEpdFullResErr[order-1][i_corr][i_cent]*100.0);
      }
    }
  }
 
  string outputname = Form("./figures/EpdEpResolutions_%s_%d_%s.pdf",recoEP::mBeamEnergy[beamEnergy].c_str(),EpdMode,JobId.c_str()); 
  string outputstart = Form("%s[",outputname.c_str());
  string outputstop = Form("%s]",outputname.c_str()); 
 
  TCanvas *c_play = new TCanvas("c_play","c_play",10,10,1200,400);
  c_play->Divide(3,1);
  for (int i_pad = 0; i_pad < 3; i_pad++)
  {
    c_play->cd(i_pad+1)->SetLeftMargin(0.15);
    c_play->cd(i_pad+1)->SetBottomMargin(0.15);
    c_play->cd(i_pad+1)->SetGrid(0,0);
    c_play->cd(i_pad+1)->SetTicks(1,1);
  }

  c_play->Print(outputstart.c_str());  
  
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

  for(int i_corr = 0; i_corr < 3; i_corr++) 
  {
    TLegend *leg = new TLegend(0.60,0.70,0.85,0.85);
    leg->SetFillColor(10);
    leg->SetBorderSize(0);
    c_play->cd(1+i_corr);
 
    //h_play->SetTitle(corrName[i_corr].c_str());
    h_play->Draw("pE");

    for(int order = 1; order <= recoEP::mEpdEpOrder; order++)
    { 
      g_mEpdFullRes[order-1][i_corr]->SetMarkerColor(color[order-1]);
      g_mEpdFullRes[order-1][i_corr]->SetLineColor(color[order-1]);
      g_mEpdFullRes[order-1][i_corr]->SetMarkerSize(1.5);
      g_mEpdFullRes[order-1][i_corr]->SetMarkerStyle(mstyle[order-1]);
      //if (order == 1) {
      //  g_mEpdFullRes[order-1][i_corr]->SetTitle("EPD EP Resolutions");
      //  g_mEpdFullRes[order-1][i_corr]->Draw("pE");
      //}
      //else {
      g_mEpdFullRes[order-1][i_corr]->Draw("pE same");
      //} 
      leg->AddEntry(g_mEpdFullRes[order-1][i_corr],Form("R_{%d}",order),"pE");
    }
    leg->Draw("same");
  }

  c_play->Update();
  c_play->Print(outputname.c_str());

  for(int order = 1; order <= recoEP::mEpdEpOrder; order++) 
  {
    TLegend *leg = new TLegend(0.60,0.70,0.85,0.85);
    leg->SetFillColor(10);
    leg->SetBorderSize(0);
    c_play->cd(order)->Clear();
    c_play->cd(order);

    //h_play->SetTitle(Form("EPD EP%d",order));
    h_play->Draw("pE");  

    for(int i_corr = 0; i_corr < 3; i_corr++)
    {
      g_mEpdFullRes[order-1][i_corr]->SetMarkerColor(color[i_corr]);
      g_mEpdFullRes[order-1][i_corr]->SetLineColor(color[i_corr]);
      g_mEpdFullRes[order-1][i_corr]->SetMarkerSize(1.5);
      g_mEpdFullRes[order-1][i_corr]->SetMarkerStyle(mstyle[i_corr]);
      //if (i_corr == 0) {
      //  g_mEpdFullRes[order-1][i_corr]->SetTitle(Form("EPD EP%d Resolutions",order));
      //  g_mEpdFullRes[order-1][i_corr]->Draw("pE");
      //}
      //else {
      g_mEpdFullRes[order-1][i_corr]->Draw("pE same");
      //} 
      leg->AddEntry(g_mEpdFullRes[order-1][i_corr],corrName[i_corr].c_str(),"pE");
    }
    leg->Draw("same");
  }

  c_play->Update();
  c_play->Print(outputname.c_str());
  
  c_play->Print(outputstop.c_str());
}
