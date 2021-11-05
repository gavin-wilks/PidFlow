#include "../StRoot/StRunQAMaker/StRunQACons.h"


char     hname[200];
TH2F*    Multd[200];
TH1D*    Multd1[200];

int tofmu = 18;
int tofmult[19]={ // TOF multiplicity bins to compare with matching TPC value
  1   , 5  , 10 , 20 , 30 , 40 , 50 , 60 , 70 , 80 ,
  100 , 140, 180, 220, 260, 300, 350, 400, 700};

//----------------------------------------------------------------------------------
void TOF_PileUpEventCorr(int energy = 0, string JobId = "028AC6C12A1F110244470E4CCDFC40FF"){
  
  TFile  weightFile(Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/AuAu%s/SpinAlignment/file_%s_RunQA_%s.root",runQA::mBeamEnergy[energy].c_str(),runQA::mBeamEnergy[energy].c_str(),JobId.c_str()), "READ");

  for(int i=0; i<1; i++){
      sprintf(hname,"h_mTofMatchRefMultAfter_trigger9"); // load in RunQA file with TPC Mult vs TPC Match
      TH2F *htmp = (TH2F*)weightFile.Get(hname);
      Multd[i] = (TH2F*)htmp->Clone(hname);
      Multd[i]->GetXaxis()->SetRangeUser(-0.5, 999.5);
      Multd[i]->GetYaxis()->SetRangeUser(-0.5, 999.5);
      Multd[i]->SetDirectory(0);

      Multd1[i] = Multd[i]->ProjectionX(); // project onto x-axis for TOF Match
  }
  //--------------------------------------------------------------------------------
  //--------------------------------------------------------------------------------
  int count[100];
  double XX[100][1000],PR[100][1000],ER[100][1000];
  double XXX[1000];
  //--------------------------------------------------------------------------------
  for(int i=0; i<tofmu; i++){ // loop of 18 TOF multiplicity bins defined at top of file

      count[i] = 0; // initialize count for specific TPC RefMult bin
      int xBinsp = Multd[0]->GetNbinsX(); // grab number of x and y bins for 2D histogram
      int yBinsp = Multd[0]->GetNbinsY();

      for(int y = 1 ; y< yBinsp+1 ; y += 1){ // loop over y bins (TPC RefMult)
          //----------------------------------------------------------------------
          //----------------------------------------------------------------------
          float Yaxis = Multd[0]->GetYaxis()->GetBinCenter(y); // grab y bin center (TPC RefMult)
          //----------------------------------------------------------------------
          float Bincont=0.0,Binerro=0.0,Binwigh=0.0; // initialize the bin error, weight and content
          for(int x = 1 ; x< xBinsp+1 ; x+= 1){ // loop over x bins (TOF Match)
              //------------------------------------------------------------------
              float Xaxis = Multd[0]->GetXaxis()->GetBinCenter(x); // get x bin center (TOF Match)
              //------------------------------------------------------------------
              if(fabs(Xaxis) >= tofmult[i] && fabs(Xaxis) < tofmult[i+1]){ // if the absolute value of the TOF Match bin is within TOF Match bin values defined at top of file 
                  //--------------------------------------------------------------
                  if(Multd[0]->GetBinContent(x,y) < 1.0) continue; // if the bin content is less than one, continue

                  Bincont  +=  Multd[0]->GetBinContent(x,y); // add the bin content to a running total for this y bin (TPC RefMult)
                  Binerro  +=  pow(Multd[0]->GetBinError(x,y),2); // add the bin error^2 to running total for this x bin for adding in quadrature later
                  Binwigh  +=  1.0; // add a count for the number of bins that meet the current if statement
                  //---------------------------------------------------------------
              }
              //-------------------------------------------------------------------
          }
          //-----------------------------------------------------------------------
          if(Binwigh < 1.0 || Bincont < 1.0) continue; // if the total weight or content are less than one, continue
          //-----------------------------------------------------------------------
          PR[i][count[i]]  =  Bincont/Binwigh; // for specific tofmult & y (TPC RefMult bin), set the mean count of the bins
          ER[i][count[i]]  =  sqrt(Binerro)/Binwigh; // calculate error of tofmult & y bin 
          XX[i][count[i]]  =  Yaxis; // set the bin center for this y (TPC RefMult) bin that has content
          count[i] ++; // add a count to the number of y bins that contain bin content
          //-----------------------------------------------------------------------
          //if(Binwigh > 0.0){
          //}
          //-----------------------------------------------------------------------
      }
      //---------------------------------------------------------------------------
      Multd1[0]->GetXaxis()->SetRangeUser( tofmult[i] , tofmult[i+1] ); // set the range of this histogram to the tofmult bin range
      XXX[i]    =  double (Multd1[0]->GetMean()); // get the mean of this histogram to plot the cut values at later (sets location of point on TPC RefMult vs TOFMatch cut)

      cout<<XXX[i]<<"\n"; // spit out the mean
  }
  //--------------------------------------------------------------------------------
  double Sig[200],Ma[200]; // initialize mean and sigma
  //================================================================
  // Move on to plotting
  // Initialize many canvases
  TCanvas * c1 = new TCanvas("Mult1", "Mult1", 1500, 500);
  c1->Divide(3,2);

  TCanvas * c2 = new TCanvas("Mult2", "Mult2", 1500, 500);
  c2->Divide(3,2);

  TCanvas * c3 = new TCanvas("Mult3", "Mult3", 1500, 500);
  c3->Divide(3,2);

  TCanvas * c4 = new TCanvas("Mult4", "Mult4", 1500, 500);
  c4->Divide(3,2);

  TCanvas * c5 = new TCanvas("Mult5", "Mult5", 1500, 500);
  c5->Divide(3,2);


  double CutMax[1000],CutMin[1000];

  string outputname = Form("./figures/%s_%s/PileUpEventCuts_%s_%s.pdf",runQA::mBeamEnergy[energy].c_str(),JobId.c_str(),runQA::mBeamEnergy[energy].c_str(),JobId.c_str()); 
  string outputstart = Form("%s[",outputname.c_str());
  string outputstop = Form("%s]",outputname.c_str());
  c1->Print(outputstart.c_str());  

  for(int xx=0; xx<tofmu; xx++){ // loop over tofbins
    sprintf(hname,"TOFMatch=[%d-%d]",tofmult[xx],tofmult[xx+1]); // set histogram title to be the tofmatch range for this bin

    if(xx < 6){ // way of creating sets of 6 histograms 2x3
      c1->cd(xx+1);
      c1->cd(xx+1)->SetLogy();
    }else if(xx < 6+6){
      c2->cd(xx-5);
      c2->cd(xx-5)->SetLogy();
    }else if(xx < 6+6+6){
      c3->cd(xx-6-5);
      c3->cd(xx-6-5)->SetLogy();
    }else if(xx < 6+6+6+6){
      c4->cd(xx-6-6-5);
      c4->cd(xx-6-6-5)->SetLogy();
    }else if(xx < 6+6+6+6+6){
      c5->cd(xx-6-6-6-5);
      c5->cd(xx-6-6-6-5)->SetLogy();
    }


    TGraphErrors *Rhist = new TGraphErrors(count[xx],XX[xx],PR[xx],0,ER[xx]); // initialize a TGraph with errors
    Rhist->SetTitle(hname);

    Rhist->SetMinimum( 0.9000);
    Rhist->SetMaximum( 19.5e+5);
    Rhist->SetLineColor(kRed);
    Rhist->SetLineWidth(1);
    Rhist->SetMarkerStyle(20);
    Rhist->SetMarkerSize(1);
    Rhist->SetMarkerColor(kRed);
    Rhist->GetXaxis()->SetTitle(" ");
    Rhist->GetXaxis()->SetLabelFont(22);
    Rhist->GetXaxis()->SetLabelSize(0.04);
    Rhist->GetXaxis()->SetTitleSize(0.07);
    Rhist->GetXaxis()->SetTitleOffset(0.67);
    //Rhist->GetXaxis()->SetLimits(-2.8,2.8);
    Rhist->GetXaxis()->SetTitleFont(22);
    Rhist->GetYaxis()->SetTitle(" ");
    Rhist->GetYaxis()->SetLabelFont(22);
    Rhist->GetYaxis()->SetLabelSize(0.04);
    Rhist->GetYaxis()->SetTitleSize(0.07);
    Rhist->GetYaxis()->SetTitleOffset(0.47);
    Rhist->GetYaxis()->SetTitleFont(22);
    Rhist->Draw("AP");

    int NMM1 =0,NMM2 =0; // setting domain for guassian fit based on tofmatch bin
    if(xx <3){
      NMM1 = 2  ;
      NMM2 = 100;
    }else if(xx<6){
      NMM1 = 4  ;
      NMM2 = 100;
    }else if(xx<12){
      NMM1 = 20 ;
      NMM2 = 200;
    }else if(xx<18){
      NMM1 = 100;
      NMM2 = 450;
    }else{
      NMM1 = 100;
      NMM2 = 450;
    }



    TF1* fx = new TF1("fx", "gaus" ,  NMM1, NMM2); // initialize gaussian fit
    fx->SetLineColor(kBlue); 
    fx->SetLineWidth(4);
    Rhist->Fit("fx","R");
  
    Sig[xx] = fx ->GetParameter(2); // grab sigma
    Ma[xx]  = fx ->GetParameter(1); // grab mean

     CutMax[xx] =  4.*(fx ->GetParameter(2)) + (fx ->GetParameter(1)); // cut on mean + 4*sigma for TPC RefMult
     CutMin[xx] = -4.*(fx ->GetParameter(2)) + (fx ->GetParameter(1)); // cut on mean - 4*sigma for TPC RefMult
  
    if(xx == 5){
      c1->Update();
      c1->Print(outputname.c_str());
    }else if(xx == 11){
      c2->Update();
      c2->Print(outputname.c_str());
    }else if(xx == 17){ 
      c3->Update();
      c3->Print(outputname.c_str());
    }else if(xx == 23){
      c4->Update();
      c4->Print(outputname.c_str());
    }else if(xx == 29){ 
      c5->Update();
      c5->Print(outputname.c_str());
    }
  }
  
  // !!!!!!!!!!!! I do not understand manually setting the min, max, and x value for the tofmu+1 bin !!!!!!!!!!!!!!!!!!!!!!!11   
  CutMax[tofmu] =  504.067; // set the final bin max cut 
  CutMin[tofmu] =  350.089; // set the final bin min cut
  XXX[tofmu]    =  520; // set the final bin mean value
 
//  CutMax[tofmu+1] =  725.458;
//  CutMin[tofmu+1] =  500.032;
//  XXX[tofmu+1]    =  3800;
  
  TCanvas * cx8 = new TCanvas("v2", "v2", 1500, 500); // initialize new canvas for final plot of all tofmatch cuts
  cx8->Divide(1);

  TGraphErrors *Rhistxx = new TGraphErrors(tofmu+1,XXX,CutMax,0,0); // graph tofmu+1 points at x mean location of individual bins and maximum cut value (gaussian mean + 4*sigma)
  Rhistxx->SetMinimum( 0.0);
  Rhistxx->SetMaximum( 9.8e+2);
  Rhistxx->SetLineColor(kRed);
  Rhistxx->SetLineWidth(1);
  Rhistxx->SetMarkerStyle(20);
  Rhistxx->SetMarkerSize(1);
  Rhistxx->SetMarkerColor(kRed);
  Rhistxx->GetXaxis()->SetTitle(" ");
  Rhistxx->GetXaxis()->SetLabelFont(22);
  Rhistxx->GetXaxis()->SetLabelSize(0.04);
  Rhistxx->GetXaxis()->SetTitleSize(0.07);
  Rhistxx->GetXaxis()->SetTitleOffset(0.67);
  //Rhistxx->GetXaxis()->SetLimits(-2.8,2.8);
  Rhistxx->GetXaxis()->SetTitleFont(22);
  Rhistxx->GetYaxis()->SetTitle(" ");
  Rhistxx->GetYaxis()->SetLabelFont(22);
  Rhistxx->GetYaxis()->SetLabelSize(0.04);
  Rhistxx->GetYaxis()->SetTitleSize(0.07);
  Rhistxx->GetYaxis()->SetTitleOffset(0.47);
  Rhistxx->GetYaxis()->SetTitleFont(22);
  Rhistxx->Draw("AP");

  TF1* fx1 = new TF1("fx1","[0] + [1]*pow(x,1) + [2]*pow(x,2) + [3]*pow(x,3) + [4]*pow(x,4) + [5]*pow(x,5) ",  0, 550); // fit with 5th degree polynomial
  fx1->SetLineColor(kBlue);
  fx1->SetLineWidth(4);
  Rhistxx->Fit("fx1","M S E R","0"); // M = TMinuit | S = fit returned as TFitResultPtr | E = Error estimation using Minos technique | R = use range specified in the function range

  cx8->Update();
  cx8->Print(outputname.c_str());

  TCanvas * cx9 = new TCanvas("v9", "v9", 1500, 500); // for minimum cuts
  cx9->Divide(1);
  cx9->cd(1);
  TGraphErrors *Rhistxy = new TGraphErrors(tofmu+1,XXX,CutMin,0,0); // graph minimum cuts
  Rhistxy->SetMinimum( 0.0);
  Rhistxy->SetMaximum( 6.8e+2);
  Rhistxy->SetLineColor(kRed);
  Rhistxy->SetLineWidth(1);
  Rhistxy->SetMarkerStyle(20);
  Rhistxy->SetMarkerSize(1);
  Rhistxy->SetMarkerColor(kRed);
  Rhistxy->GetXaxis()->SetTitle(" ");
  Rhistxy->GetXaxis()->SetLabelFont(22);
  Rhistxy->GetXaxis()->SetLabelSize(0.04);
  Rhistxy->GetXaxis()->SetTitleSize(0.07);
  Rhistxy->GetXaxis()->SetTitleOffset(0.67);
  //Rhistxy->GetXaxis()->SetLimits(-2.8,2.8);
  Rhistxy->GetXaxis()->SetTitleFont(22);
  Rhistxy->GetYaxis()->SetTitle(" ");
  Rhistxy->GetYaxis()->SetLabelFont(22);
  Rhistxy->GetYaxis()->SetLabelSize(0.04);
  Rhistxy->GetYaxis()->SetTitleSize(0.07);
  Rhistxy->GetYaxis()->SetTitleOffset(0.47);
  Rhistxy->GetYaxis()->SetTitleFont(22);
  Rhistxy->Draw("AP");

  TF1* fx2 = new TF1("fx2","[0] + [1]*pow(x,1) + [2]*pow(x,2) + [3]*pow(x,3) + [4]*pow(x,4) + [5]*pow(x,5) ",  0, 550);
  fx2->SetLineColor(kBlue) ;
  fx2->SetLineWidth(4);
  Rhistxy->Fit("fx2","M S E R","0");
  
  cx9->Update();
  cx9->Print(outputname.c_str());

  for(int xx=0; xx<tofmu; xx++){ // spit out min and max cuts
    cout<<xx<<"  "<<XXX[xx]<<" "<<CutMax[xx]<<" "<<CutMin[xx]<<"\n";
  }
  
  TCanvas * cx2d = new TCanvas("2d", "2d", 1500, 500); 
  cx2d->Divide(1);
  cx2d->cd(1);
  //Multd[0]->SetMinimum( 0.0);
  //Multd[0]->SetMaximum( 6.8e+2);
  //Multd[0]->SetLineColor(kRed);
  //Multd[0]->SetLineWidth(1);
  //Multd[0]->SetMarkerStyle(20);
  //Multd[0]->SetMarkerSize(1);
  //Multd[0]->SetMarkerColor(kRed);
  Multd[0]->GetXaxis()->SetTitle(" ");
  Multd[0]->GetXaxis()->SetLabelFont(22);
  Multd[0]->GetXaxis()->SetLabelSize(0.04);
  Multd[0]->GetXaxis()->SetTitleSize(0.07);
  Multd[0]->GetXaxis()->SetTitleOffset(0.67);
  Multd[0]->GetXaxis()->SetTitleFont(22);
  Multd[0]->GetYaxis()->SetTitle(" ");
  Multd[0]->GetYaxis()->SetLabelFont(22);
  Multd[0]->GetYaxis()->SetLabelSize(0.04);
  Multd[0]->GetYaxis()->SetTitleSize(0.07);
  Multd[0]->GetYaxis()->SetTitleOffset(0.47);
  Multd[0]->GetYaxis()->SetTitleFont(22);
  Multd[0]->Draw("colz");

  fx1->Draw("same");
  fx2->Draw("same");

  cx2d->Update();
  cx2d->Print(outputname.c_str());
  
  c1->Print(outputstop.c_str());

//return 0;
}
//***********************
//=======================
