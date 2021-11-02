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

void plotQA_Track_PID(int energy = 0, string JobId = "028AC6C12A1F110244470E4CCDFC40FF")
{
  string inputfile = Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/AuAu%s/SpinAlignment/file_%s_RunQA_%s.root",runQA::mBeamEnergy[energy].c_str(),runQA::mBeamEnergy[energy].c_str(),JobId.c_str()); 
 
  TFile *File_InPut = TFile::Open(inputfile.c_str());
  TH2F *h_mDEdxMom[2];
  TH2F *h_mBetaMom[2];
  TH2F *h_mMass2Mom[2];
  TH1F *h_mMass2[2];
 
  for(int i_cut = 0; i_cut < 2; ++i_cut)
  {
    std::string HistName = Form("h_mDEdxMom%s_trigger9",mCutsQA[i_cut].c_str());
    h_mDEdxMom[i_cut] = (TH2F*)File_InPut->Get(HistName.c_str());

    HistName = Form("h_mMass2Mom%s_trigger9",mCutsQA[i_cut].c_str());
    h_mMass2Mom[i_cut] = (TH2F*)File_InPut->Get(HistName.c_str());

    HistName = Form("h_mBetaMom%s_trigger9",mCutsQA[i_cut].c_str());
    h_mBetaMom[i_cut] = (TH2F*)File_InPut->Get(HistName.c_str());
   
    h_mMass2[i_cut] = (TH1F*)h_mMass2Mom[i_cut]->ProjectionY();
  }

  TCanvas *c_TrackQA_PID = new TCanvas("c_TrackQA_PID","c_TrackQA_PID",10,10,1600,800);
  c_TrackQA_PID->Divide(4,2);
  for(int i_pad = 0; i_pad < 8; ++i_pad)
  {
    c_TrackQA_PID->cd(i_pad+1);
    c_TrackQA_PID->cd(i_pad+1)->SetLeftMargin(0.1);
    c_TrackQA_PID->cd(i_pad+1)->SetRightMargin(0.1);
    c_TrackQA_PID->cd(i_pad+1)->SetBottomMargin(0.1);
    c_TrackQA_PID->cd(i_pad+1)->SetGrid(0,0);
    c_TrackQA_PID->cd(i_pad+1)->SetTicks(1,1);
    c_TrackQA_PID->cd(i_pad+1)->SetLogz();
  }

  for(int i_cut = 0; i_cut < 2; ++i_cut)
  {
    c_TrackQA_PID->cd(i_cut*4+1);
    h_mDEdxMom[i_cut]->Draw("colz");

    c_TrackQA_PID->cd(i_cut*4+2);
    h_mMass2Mom[i_cut]->Draw("colz");

    c_TrackQA_PID->cd(i_cut*4+3);
    h_mBetaMom[i_cut]->Draw("colz");

    c_TrackQA_PID->cd(i_cut*4+4);
    h_mMass2[i_cut]->Draw("hE");
  }

  string FigName = Form("./figures/%s/c_TrackQA_PID_%s_%s.pdf",runQA::mBeamEnergy[energy].c_str(),runQA::mBeamEnergy[energy].c_str(),JobId.c_str());
  c_TrackQA_PID->SaveAs(FigName.c_str());

}
