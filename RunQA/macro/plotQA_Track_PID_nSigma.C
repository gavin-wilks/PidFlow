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

static const string mParticle[4] = {"Pi","K","P","E"};

void plotQA_Track_PID_nSigma(int energy = 0, string JobId = "028AC6C12A1F110244470E4CCDFC40FF")
{
  string inputfile = Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/AuAu%s/SpinAlignment/file_%s_RunQA_%s.root",runQA::mBeamEnergy[energy].c_str(),runQA::mBeamEnergy[energy].c_str(),JobId.c_str()); 
 
  TFile *File_InPut = TFile::Open(inputfile.c_str());
  TH2F *h_mNSigmaRig[4][2];
  TH2F *h_mNSigmaEta[4][2];
  TH2F *h_mNSigmaPhi[4][2];

  for(int i_par = 0; i_par < 4; ++i_par)
  {
    for(int i_cut = 0; i_cut < 2; ++i_cut)
    {
      std::string HistName = Form("h_mNSigma%sRig%s_trigger9",mParticle[i_par].c_str(),mCutsQA[i_cut].c_str());
      h_mNSigmaRig[i_par][i_cut] = (TH2F*)File_InPut->Get(HistName.c_str());

      HistName = Form("h_mNSigma%sEta%s_trigger9",mParticle[i_par].c_str(),mCutsQA[i_cut].c_str());
      h_mNSigmaEta[i_par][i_cut] = (TH2F*)File_InPut->Get(HistName.c_str());

      HistName = Form("h_mNSigma%sPhi%s_trigger9",mParticle[i_par].c_str(),mCutsQA[i_cut].c_str());
      h_mNSigmaPhi[i_par][i_cut] = (TH2F*)File_InPut->Get(HistName.c_str());
    }
  }

  //std::string outputname = Form("./figures/%s/c_TrackQA_PID_nSigma_%s_%s.pdf",runQA::mBeamEnergy[energy].c_str(),runQA::mBeamEnergy[energy].c_str(),JobId.c_str());

  TCanvas *c_TrackQA_PID = new TCanvas("c_TrackQA_PID","c_TrackQA_PID",10,10,1200,800);
  c_TrackQA_PID->Divide(3,2);
  for(int i_pad = 0; i_pad < 6; ++i_pad)
  {
    c_TrackQA_PID->cd(i_pad+1);
    c_TrackQA_PID->cd(i_pad+1)->SetLeftMargin(0.1);
    c_TrackQA_PID->cd(i_pad+1)->SetRightMargin(0.1);
    c_TrackQA_PID->cd(i_pad+1)->SetBottomMargin(0.1);
    c_TrackQA_PID->cd(i_pad+1)->SetGrid(0,0);
    c_TrackQA_PID->cd(i_pad+1)->SetTicks(1,1);
    c_TrackQA_PID->cd(i_pad+1)->SetLogz();
  }

  //std::string output_start = Form("%s[",outputname.c_str());
  //c_TrackQA_PID->Print(output_start.c_str()); // open pdf file

  for(int i_par = 0; i_par < 4; ++i_par)
  {
    for(int i_cut = 0; i_cut < 2; ++i_cut)
    {
      c_TrackQA_PID->cd(i_cut*3+1);
      h_mNSigmaRig[i_par][i_cut]->Draw("colz");

      c_TrackQA_PID->cd(i_cut*3+2);
      h_mNSigmaEta[i_par][i_cut]->Draw("colz");

      c_TrackQA_PID->cd(i_cut*3+3);
      h_mNSigmaPhi[i_par][i_cut]->Draw("colz");
    }
    //c_TrackQA_PID->Update();
    std::string outputname = Form("./figures/%s/c_TrackQA_PID_nSigma%s_%s_%s.pdf",runQA::mBeamEnergy[energy].c_str(),mParticle[i_par].c_str(),runQA::mBeamEnergy[energy].c_str(),JobId.c_str());
    c_TrackQA_PID->Print(outputname.c_str());
  }
  
  //std::string output_stop = Form("%s]",outputname.c_str());
  //c_TrackQA_PID->Print(output_stop.c_str());

}
