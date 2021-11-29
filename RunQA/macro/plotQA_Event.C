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

double TofMatchLow(double *x, double *par)
{
  // 5th order polynomial coefficients for lower cut, l
  double p0l = -13.11    ;
  double p1l = 0.8207    ;
  double p2l = -4.241e-3 ;
  double p3l = 2.81e-5   ;
  double p4l = -6.434e-8 ;
  double p5l = 4.833e-11 ; 
  return p0l + p1l*x[0] + p2l*pow(x[0],2) + p3l*pow(x[0],3) + p4l*pow(x[0],4) + p5l*pow(x[0],5); 
}

double TofMatchHigh(double *x, double *par)
{
  // 5th order polynomial coefficients for higher cut, h
  double p0h = 10.07     ;
  double p1h = 1.417     ;   
  double p2h = 1.979e-4  ;
  double p3h = -4.87e-6  ;
  double p4h = 1.109e-8  ;
  double p5h = -1.111e-11;
  return p0h + p1h*x[0] + p2h*pow(x[0],2) + p3h*pow(x[0],3) + p4h*pow(x[0],4) + p5h*pow(x[0],5);   
}          
           
static const string mCutsQA[2] = {"Before","After"};
           
void plotQA_Event(int energy = 0, string JobId = "028AC6C12A1F110244470E4CCDFC40FF")
{          
  string inputfile = Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/AuAu%s/SpinAlignment/file_%s_RunQA_%s.root",runQA::mBeamEnergy[energy].c_str(),runQA::mBeamEnergy[energy].c_str(),JobId.c_str());
  

  TFile *File_InPut = TFile::Open(inputfile.c_str());
  TH1F *h_mRefMult[2][10]; // 0: before cuts | 1: after cuts
  TH1F *h_mGRefMult[2][10]; // 0-8 for different triggerID | 9 for all triggers
  TH2F *h_mRefMultGRefMult[2][10];
  TH1F *h_mCentrality9[2][10];
  TH2F *h_mTofMatchRefMult[2][10];
  TH2F *h_mTofHitsRefMult[2][10];
  TH2F *h_mTofMatchGRefMult[2][10];
  TH1F *h_mTofHits[2][10];
  TH1F *h_mTofMatch[2][10];
  TH2F *h_mTofHitsGRefMult[2][10];
  TH2F *h_mVzVzVpd[2][10];
  TH1F *h_mDiffVzVzVpd[2][10];
  TH1F *h_mVertexZ[2][10];
  TH2F *h_mVertexXY[2][10];  


  for(int i_cut = 0; i_cut < 2; ++i_cut)
  {
    for(int i_trig = 0; i_trig < 10; ++i_trig)
    {
      std::string HistName = Form("h_mRefMult%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      h_mRefMult[i_cut][i_trig] = (TH1F*)File_InPut->Get(HistName.c_str());
      h_mRefMult[i_cut][i_trig]->SetLineColor(i_cut+1);
      h_mRefMult[i_cut][i_trig]->GetXaxis()->SetTitle("refMult");

      HistName = Form("h_mGRefMult%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      h_mGRefMult[i_cut][i_trig] = (TH1F*)File_InPut->Get(HistName.c_str());;
      h_mGRefMult[i_cut][i_trig]->SetLineColor(i_cut+1);
      h_mGRefMult[i_cut][i_trig]->GetXaxis()->SetTitle("gRefMult");

      HistName = Form("h_mRefMultGRefMult%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      h_mRefMultGRefMult[i_cut][i_trig] = (TH2F*)File_InPut->Get(HistName.c_str());
      h_mRefMultGRefMult[i_cut][i_trig]->GetXaxis()->SetTitle("refMult");
      h_mRefMultGRefMult[i_cut][i_trig]->GetYaxis()->SetTitle("gRefMult");

      /*HistName = Form("h_mCentrality9%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      h_mCentrality9[i_cut][i_trig] = (TH1F*)File_InPut->Get(HistName.c_str());
      h_mCentrality9[i_cut][i_trig]->SetLineColor(i_cut+1);
      h_mCentrality9[i_cut][i_trig]->GetXaxis()->SetTitle("centrality");
      */
      HistName = Form("h_mTofMatchRefMult%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      h_mTofMatchRefMult[i_cut][i_trig] = (TH2F*)File_InPut->Get(HistName.c_str());
      h_mTofMatchRefMult[i_cut][i_trig]->GetXaxis()->SetTitle("tofMatch");
      h_mTofMatchRefMult[i_cut][i_trig]->GetYaxis()->SetTitle("refMult");

      HistName = Form("h_mTofHitsRefMult%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      h_mTofHitsRefMult[i_cut][i_trig] = (TH2F*)File_InPut->Get(HistName.c_str());
      h_mTofHitsRefMult[i_cut][i_trig]->GetXaxis()->SetTitle("tofHits");
      h_mTofHitsRefMult[i_cut][i_trig]->GetYaxis()->SetTitle("refMult");

      h_mTofHits[i_cut][i_trig] = (TH1F*) h_mTofHitsRefMult[i_cut][i_trig]->ProjectionX();
      h_mTofMatch[i_cut][i_trig] = (TH1F*) h_mTofMatchRefMult[i_cut][i_trig]->ProjectionX();
    
      HistName = Form("h_mTofMatchGRefMult%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      h_mTofMatchGRefMult[i_cut][i_trig] = (TH2F*)File_InPut->Get(HistName.c_str());;
      h_mTofMatchGRefMult[i_cut][i_trig]->GetXaxis()->SetTitle("tofMatch");
      h_mTofMatchGRefMult[i_cut][i_trig]->GetYaxis()->SetTitle("gRefMult");

      HistName = Form("h_mTofHitsGRefMult%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      h_mTofHitsGRefMult[i_cut][i_trig] = (TH2F*)File_InPut->Get(HistName.c_str());
      h_mTofHitsGRefMult[i_cut][i_trig]->GetXaxis()->SetTitle("tofHits");
      h_mTofHitsGRefMult[i_cut][i_trig]->GetYaxis()->SetTitle("gRefMult");

      std::string HistName = Form("h_mVertexXY%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      h_mVertexXY[i_cut][i_trig] = (TH2F*)File_InPut->Get(HistName.c_str());
      h_mVertexXY[i_cut][i_trig]->GetXaxis()->SetTitle("Vx");
      h_mVertexXY[i_cut][i_trig]->GetYaxis()->SetTitle("Vy");

      HistName = Form("h_mVertexZ%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      h_mVertexZ[i_cut][i_trig] = (TH1F*)File_InPut->Get(HistName.c_str());
      h_mVertexZ[i_cut][i_trig]->SetLineColor(i_cut+1);
      h_mVertexZ[i_cut][i_trig]->GetXaxis()->SetTitle("Vz");

      HistName = Form("h_mVzVzVpd%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      h_mVzVzVpd[i_cut][i_trig] = (TH2F*)File_InPut->Get(HistName.c_str());
      h_mVzVzVpd[i_cut][i_trig]->GetXaxis()->SetTitle("Vz");
      h_mVzVzVpd[i_cut][i_trig]->GetYaxis()->SetTitle("VzVpd");

      HistName = Form("h_mDiffVzVzVpd%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      h_mDiffVzVzVpd[i_cut][i_trig] = (TH1F*)File_InPut->Get(HistName.c_str());
      h_mDiffVzVzVpd[i_cut][i_trig]->SetLineColor(i_cut+1);
      h_mDiffVzVzVpd[i_cut][i_trig]->GetXaxis()->SetTitle("Vz-VzVpd");
    }
  }

  TCanvas *c_EventQA[2];
  for(int i_cut = 0; i_cut < 2; ++i_cut)
  {
    string CanName = Form("c_EventQA_%s",mCutsQA[i_cut].c_str());
    c_EventQA[i_cut] = new TCanvas(CanName.c_str(),CanName.c_str(),10,10,2000,800);
    c_EventQA[i_cut]->Divide(5,2);
    for(int i_pad = 0; i_pad < 10; ++i_pad)
    {
      c_EventQA[i_cut]->cd(i_pad+1);
      c_EventQA[i_cut]->cd(i_pad+1)->SetLeftMargin(0.1);
      c_EventQA[i_cut]->cd(i_pad+1)->SetRightMargin(0.1);
      c_EventQA[i_cut]->cd(i_pad+1)->SetBottomMargin(0.1);
      c_EventQA[i_cut]->cd(i_pad+1)->SetGrid(0,0);
      c_EventQA[i_cut]->cd(i_pad+1)->SetTicks(1,1);
    }

    c_EventQA[i_cut]->cd(1);
    c_EventQA[i_cut]->cd(1)->SetLogy();
    if(energy == 0)h_mRefMult[i_cut][9]->Draw("hE");
    if(energy != 0)h_mRefMult[i_cut][9]->Draw("hE");
    /////////////////////////////////////
    TF1 *TofLow = new TF1("TofLow",TofMatchLow,0,550,0);
    TF1 *TofHigh = new TF1("TofHigh",TofMatchHigh,0,550,0);
    TofLow->SetLineWidth(0.5);
    TofHigh->SetLineWidth(0.5);

    c_EventQA[i_cut]->cd(2);
    c_EventQA[i_cut]->cd(2)->SetLogz();
    
    h_mTofMatchRefMult[i_cut][9]->GetXaxis()->SetRangeUser(0.0,700.0);
    h_mTofMatchRefMult[i_cut][9]->GetYaxis()->SetRangeUser(0.0,700.0);
    h_mTofMatchRefMult[i_cut][9]->Draw("colz");
    TofLow->Draw("same");
    TofHigh->Draw("same");
    ////////////////////////////
    c_EventQA[i_cut]->cd(3);
    c_EventQA[i_cut]->cd(3)->SetLogz();
     
    h_mTofHitsRefMult[i_cut][9]->GetYaxis()->SetRangeUser(0.0,700.0);
    h_mTofHitsRefMult[i_cut][9]->GetXaxis()->SetRangeUser(0.0,1500.0);
    h_mTofHitsRefMult[i_cut][9]->Draw("colz");
    //////////////////////
    c_EventQA[i_cut]->cd(4);
    c_EventQA[i_cut]->cd(4)->SetLogy();
    h_mTofMatch[i_cut][9]->Draw("hE");
   
    c_EventQA[i_cut]->cd(5);
    c_EventQA[i_cut]->cd(5)->SetLogy();
    h_mTofHits[i_cut][9]->Draw("hE");

    /*c_EventQA[i_cut]->cd(4);
    c_EventQA[i_cut]->cd(4)->SetLogy();
    h_mCentrality9[i_cut][9]->Draw("hE");
    */

    c_EventQA[i_cut]->cd(7);
    c_EventQA[i_cut]->cd(7)->SetLogz();
    h_mVertexXY[i_cut][9]->Draw("colz");

    c_EventQA[i_cut]->cd(8);
    // c_EventQA[i_cut]->cd(6)->SetLogy();
    h_mVertexZ[i_cut][9]->Draw();

    c_EventQA[i_cut]->cd(9);
    c_EventQA[i_cut]->cd(9)->SetLogz();
    h_mVzVzVpd[i_cut][9]->Draw("colz");

    c_EventQA[i_cut]->cd(10);
    c_EventQA[i_cut]->cd(10)->SetLogy();
    // c_EventQA[i_cut]->cd(8)->SetLogy();
    h_mDiffVzVzVpd[i_cut][9]->Draw();

    string FigName = Form("./figures/%s/c_EventQA_%s_%s_%s.pdf",runQA::mBeamEnergy[energy].c_str(),mCutsQA[i_cut].c_str(),runQA::mBeamEnergy[energy].c_str(),JobId.c_str());
    c_EventQA[i_cut]->SaveAs(FigName.c_str());
  }
}
