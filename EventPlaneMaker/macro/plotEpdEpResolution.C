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

void plotEpdEpResolution(int beamEnergy = 1)
{

  TProfile *p_mResolutionRaw[3];
  TProfile *p_mResolutionShift[3];
  int color[3] = {1,2,4};

  //string inputfile = Form("../StRoot/StEventPlaneUtility/Resolution/file_%s_Resolution.root",recoEP::mBeamEnergy[beamEnergy].c_str());
  string inputfile = Form("../file_%s_EpdResults_0_1.root", recoEP::mBeamEnergy[beamEnergy]);
  TFile *File_InPut = TFile::Open(inputfile.c_str());
 
  for(int o = 1; o <= 3; o++){
    string HistName = Form("EpdCos%draw",o);
    p_mResolutionRaw[o-1] = (TProfile*)File_InPut->Get(HistName.c_str());
    HistName = Form("EpdCos%dshift",o);
    p_mResolutionShift[o-1] = (TProfile*)File_InPut->Get(HistName.c_str());
  }

  TCanvas *c_play = new TCanvas("c_play","c_play",10,10,800,800);
  c_play->SetLeftMargin(0.15);
  c_play->SetBottomMargin(0.15);
  c_play->SetGrid(0,0);
  c_play->SetTicks(1,1);
  c_play->cd();

  TLegend *leg = new TLegend(0.60,0.70,0.85,0.85);
  leg->SetFillColor(10);
  leg->SetBorderSize(0);
  for(int o = 1; o <= 3; o++) {
    p_mResolutionRaw[o-1]->SetMarkerColor(color[o-1]);
    p_mResolutionRaw[o-1]->SetLineColor(color[o-1]);
    if (o == 1) {
      p_mResolutionRaw[o-1]->SetTitle("EP Resolutions");
      p_mResolutionRaw[o-1]->GetYaxis()->SetRangeUser(0,0.1);
      //p_mResolution[o-1]->GetXaxis()->SetTitle("");
      p_mResolutionRaw[o-1]->Draw("pE");
    }
    else {
      p_mResolutionRaw[o-1]->Draw("pE same");
    } 
    p_mResolutionRaw[o-1]->SetStats(0);
    leg->AddEntry(p_mResolutionRaw[o-1],Form("EP%d",o),"pE");
  }   
  leg->Draw("same");

  string FigureName = Form("./figures/EpdEpResolutionRaw_%s.pdf",recoEP::mBeamEnergy[beamEnergy].c_str());
  c_play->SaveAs(FigureName.c_str());

  c_play->cd()->Clear();

  TLegend *leg = new TLegend(0.60,0.70,0.85,0.85);
  leg->SetFillColor(10);
  leg->SetBorderSize(0);
  for(int o = 1; o <= 3; o++) {
    p_mResolutionShift[o-1]->SetMarkerColor(color[o-1]);
    p_mResolutionShift[o-1]->SetLineColor(color[o-1]);
    if (o == 1) {
      p_mResolutionShift[o-1]->SetTitle("EP Resolutions");
      p_mResolutionShift[o-1]->GetYaxis()->SetRangeUser(0,0.1);
      //p_mResolution[o-1]->GetXaxis()->SetTitle("");
      p_mResolutionShift[o-1]->Draw("pE");
    }
    else {
      p_mResolutionShift[o-1]->Draw("pE same");
    } 
    p_mResolutionShift[o-1]->SetStats(0);
    leg->AddEntry(p_mResolutionShift[o-1],Form("EP%d",o),"pE");
  }   
  leg->Draw("same");

  string FigureName = Form("./figures/EpdEpResolutionShift_%s.pdf",recoEP::mBeamEnergy[beamEnergy].c_str());
  c_play->SaveAs(FigureName.c_str()); 

}
