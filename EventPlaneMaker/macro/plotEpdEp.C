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

void plotEpdEp(int beamEnergy = 1, int EpdMode = 2, string JobId = "1")
{
  string corrName[3] = {"Raw", "Pw", "PwS"}; // Raw, phi-weighted, phi-weighted and shifted
  string psiName[2] = {"East", "West"};
  cout << "did we get here" << endl;
  TH2F *h_mEwPsi[recoEP::mEpdEpOrder][3];   // load from file
  TH1F *h_mEPsi[recoEP::mEpdEpOrder][3];
  TH1F *h_mWPsi[recoEP::mEpdEpOrder][3];
  TH1F *h_mFullPsi[recoEP::mEpdEpOrder][3]; // load from file
  TH2F *h_mQyQx[2][recoEP::mEpdEpOrder][2]; // load from file
  TH1F *h_mQx[2][recoEP::mEpdEpOrder][2];
  TH1F *h_mQy[2][recoEP::mEpdEpOrder][2]; 

  int color[4] = {1,2,4,3};

  //string inputfile = Form("../StRoot/StEventPlaneUtility/Resolution/file_%s_Resolution.root",recoEP::mBeamEnergy[beamEnergy].c_str());
  string inputfile = Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/AuAu%s/OutPut/PidFlow/file_%s_EpdResults_%d_%s.root",recoEP::mBeamEnergy[beamEnergy],recoEP::mBeamEnergy[beamEnergy],EpdMode,JobId.c_str());
  TFile *File_InPut = TFile::Open(inputfile.c_str());
 cout << "did we get here" << endl;
  for(int c = 0; c < 3; c++) 
  { // loop over correction
    for(int o = 1; o <= 3; o++)
    { // loop over order
      string HistName = Form("EpdEwPsi%d%s",o,corrName[c].c_str());
      h_mEwPsi[o-1][c] = (TH2F*)File_InPut->Get(HistName.c_str());

      HistName = Form("EpdWPsi%d%s",o,corrName[c].c_str());
      string HistTitle = Form("EPD %s #Psi_{w,%d}",corrName[c].c_str(),o);
      h_mWPsi[o-1][c] = (TH1F*)h_mEwPsi[o-1][c]->ProjectionX(HistName.c_str()); 
      h_mWPsi[o-1][c]->SetTitle(HistTitle.c_str());

      HistName = Form("EpdEPsi%d%s",o,corrName[c].c_str());
      HistTitle = Form("EPD %s #Psi_{e,%d}",corrName[c].c_str(),o);      
      h_mEPsi[o-1][c] = (TH1F*)h_mEwPsi[o-1][c]->ProjectionY(HistName.c_str());   
      h_mEPsi[o-1][c]->SetTitle(HistTitle.c_str());
      
      HistName = Form("EpdFullPsi%d%s",o,corrName[c].c_str());
      h_mFullPsi[o-1][c] = (TH1F*)File_InPut->Get(HistName.c_str());
      
      if(c > 1) continue; // there is no shifting correction for Q vector

      for(int p = 0; p < 2; p++)
      { // loop over psi east, west
        HistName = Form("EpdQyQx%s%d%s",psiName[p].c_str(),o,corrName[c].c_str());
        h_mQyQx[p][o-1][c] = (TH2F*)File_InPut->Get(HistName.c_str()); 
 
        HistName = Form("EpdQx%s%d%s",psiName[p].c_str(),o,corrName[c].c_str());
        HistTitle = Form("EPD %s %s Q_{x,%d}",corrName[c].c_str(),psiName[p].c_str(),o);
        h_mQx[p][o-1][c] = (TH1F*)h_mQyQx[p][o-1][c]->ProjectionX(HistName.c_str());
        h_mQx[p][o-1][c]->SetTitle(HistTitle.c_str());        

        HistName = Form("EpdQy%s%d%s",psiName[p].c_str(),o,corrName[c].c_str());
        HistTitle = Form("EPD %s %s Q_{y,%d}",corrName[c].c_str(),psiName[p].c_str(),o);
        h_mQy[p][o-1][c] = (TH1F*)h_mQyQx[p][o-1][c]->ProjectionY(HistName.c_str());
        h_mQy[p][o-1][c]->SetTitle(HistTitle.c_str());
      }
    }
  }
cout << "did we get here" << endl;
  string outputname = Form("./figures/EpdEpResults_%s_%d_%s.pdf",recoEP::mBeamEnergy[beamEnergy].c_str(),EpdMode,JobId.c_str()); 
  string outputstart = Form("%s[",outputname.c_str());
  string outputstop = Form("%s]",outputname.c_str());

  TCanvas *c_play = new TCanvas("c_play","c_play",10,10,1600,1350);
  c_play->Divide(3,3);
  for(int i = 0; i < 9; i++)
  {
    c_play->cd(i+1)->SetLeftMargin(0.15);
    //c_play->cd(i+1)->SetRightMargin(0.10);
    c_play->cd(i+1)->SetBottomMargin(0.15);
    c_play->cd(i+1)->SetGrid(0,0);
    c_play->cd(i+1)->SetTicks(1,1);
  }
  TCanvas *c_playQ = new TCanvas("c_playQ","c_playQ",10,10,1600,900);
  c_playQ->Divide(3,2);
  for(int i = 0; i < 6; i++)
  {
    c_playQ->cd(i+1)->SetLeftMargin(0.15);
    //c_playQ->cd(i+1)->SetRightMargin(0.10);
    c_playQ->cd(i+1)->SetBottomMargin(0.15);
    c_playQ->cd(i+1)->SetGrid(0,0);
    c_playQ->cd(i+1)->SetTicks(1,1);
  }

  c_play->Print(outputstart.c_str());

  for(int c = 0; c < 3; c++)
  {
    for(int o = 1; o <= 3; o++) 
    {
      c_play->cd(3*c+o);
      h_mFullPsi[o-1][c]->Draw("hE"); 
    }
  }

  c_play->Update();
  c_play->Print(outputname.c_str());

  for(int c = 0; c < 3; c++)
  {
    for(int o = 1; o <= 3; o++) 
    {
      c_play->cd(3*c+o);
      c_play->cd(3*c+o)->SetLogz();
      h_mEwPsi[o-1][c]->Draw("colz"); 
    }
  }

  c_play->Update();
  c_play->Print(outputname.c_str());

  for(int c = 0; c < 3; c++)
  {
    for(int o = 1; o <= 3; o++) 
    {
      c_play->cd(3*c+o);
      h_mEPsi[o-1][c]->Draw("hE"); 
    }
  }
  
  c_play->Update();
  c_play->Print(outputname.c_str());

  for(int c = 0; c < 3; c++)
  {
    for(int o = 1; o <= 3; o++) 
    {
      c_play->cd(3*c+o);
      h_mWPsi[o-1][c]->Draw("hE"); 
    }
  } 

  c_play->Update();
  c_play->Print(outputname.c_str());
 
  for(int c = 0; c < 2; c++)
  {
    for(int o = 1; o <= 3; o++) 
    {
      c_playQ->cd(3*c+o);
      c_playQ->cd(3*c+o)->SetLogz();
      h_mQyQx[0][o-1][c]->Draw("colz"); 
    }
  }
  
  c_playQ->Update();
  c_playQ->Print(outputname.c_str());
 
  for(int c = 0; c < 2; c++)
  {
    for(int o = 1; o <= 3; o++) 
    {
      c_playQ->cd(3*c+o);
      c_playQ->cd(3*c+o)->SetLogz();
      h_mQyQx[1][o-1][c]->Draw("colz"); 
    }
  }
 
  c_playQ->Update();
  c_playQ->Print(outputname.c_str());
 
  for(int c = 0; c < 2; c++)
  {
    for(int o = 1; o <= 3; o++) 
    {
      c_playQ->cd(3*c+o);
      c_playQ->cd(3*c+o)->SetLogy();
      h_mQx[0][o-1][c]->Draw("hE"); 
    }
  }
  
  c_playQ->Update();
  c_playQ->Print(outputname.c_str());

  for(int c = 0; c < 2; c++)
  {
    for(int o = 1; o <= 3; o++) 
    {
      c_playQ->cd(3*c+o);
      c_playQ->cd(3*c+o)->SetLogy();
      h_mQx[1][o-1][c]->Draw("hE"); 
    }
  }
 
  c_playQ->Update();
  c_playQ->Print(outputname.c_str());
 
  for(int c = 0; c < 2; c++)
  {
    for(int o = 1; o <= 3; o++) 
    {
      c_playQ->cd(3*c+o);
      c_playQ->cd(3*c+o)->SetLogy();
      h_mQy[0][o-1][c]->Draw("hE"); 
    }
  }
  
  c_playQ->Update();
  c_playQ->Print(outputname.c_str());

  for(int c = 0; c < 2; c++)
  {
    for(int o = 1; o <= 3; o++) 
    {
      c_playQ->cd(3*c+o);
      c_playQ->cd(3*c+o)->SetLogy();
      h_mQy[1][o-1][c]->Draw("hE"); 
    }
  }
  
  c_playQ->Update();
  c_playQ->Print(outputname.c_str());
 
  c_play->Print(outputstop.c_str());
}
