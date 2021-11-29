


char     hname[200];
TH2D*    Multd[200];
TH1D*    Multd1[200];

 int tofmu = 18;
 int tofmult[19]={
  1   , 5  , 10 , 20 , 30 , 40 , 50 , 60 , 70 , 80 ,
  100 , 140, 180, 220, 260, 300, 350, 400, 700};

//----------------------------------------------------------------------------------
double TOF_Mult(){
  

  TFile  weightFile("./Au_19_QA_R21.root", "READ");

  for(int i=0; i<1; i++){
      sprintf(hname,"Mylt_Mult1A_4");
      TH2D *htmp = (TH2D*)weightFile.Get(hname);
      Multd[i] = (TH2D*)htmp->Clone(hname);
      Multd[i]->SetDirectory(0);

      Multd1[i]= Multd[i]->ProjectionY();
  }
  //--------------------------------------------------------------------------------
  //--------------------------------------------------------------------------------
  int count[100];
  double XX[100][1000],PR[100][1000],ER[100][1000];
  double XXX[1000];
  //--------------------------------------------------------------------------------
  for(int i=0; i<tofmu; i++){

      count[i] = 0;
      int xBinsp = Multd[0]->GetNbinsX();
      int yBinsp = Multd[0]->GetNbinsY();

      for(int x = 1 ; x< xBinsp+1 ; x += 1){
          //----------------------------------------------------------------------
          //----------------------------------------------------------------------
          float Xaxis = Multd[0]->GetXaxis()->GetBinCenter(x);
          //----------------------------------------------------------------------
          float Bincont=0.0,Binerro=0.0,Binwigh=0.0;
          for(int y = 1 ; y< yBinsp+1 ; y+= 1){
              //------------------------------------------------------------------
              float Yaxis = Multd[0]->GetYaxis()->GetBinCenter(y);
              //------------------------------------------------------------------
              if(fabs(Yaxis) >= tofmult[i] && fabs(Yaxis) < tofmult[i+1]){
                  //--------------------------------------------------------------
                  if(Multd[0]->GetBinContent(x,y) < 1.0) continue;

                  Bincont  +=  Multd[0]->GetBinContent(x,y);
                  Binerro  +=  pow(Multd[0]->GetBinError(x,y),2);
                  Binwigh  +=  1.0;
                  //---------------------------------------------------------------
              }
              //-------------------------------------------------------------------
          }
          //-----------------------------------------------------------------------
          if(Binwigh < 1.0 || Bincont < 1.0) continue;
          //-----------------------------------------------------------------------
          PR[i][count[i]]  =  Bincont/Binwigh;
          ER[i][count[i]]  =  sqrt(Binerro)/Binwigh;
          XX[i][count[i]]  =  Xaxis;
          count[i] ++;
          //-----------------------------------------------------------------------
          //if(Binwigh > 0.0){
          //}
          //-----------------------------------------------------------------------
      }
      //---------------------------------------------------------------------------
      Multd1[0]->GetXaxis()->SetRangeUser( tofmult[i] , tofmult[i+1] );
      XXX[i]    =  double (Multd1[0]->GetMean());

      cout<<XXX[i]<<"\n";
  }
  //--------------------------------------------------------------------------------
  double Sig[200],Ma[200];
  //================================================================
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

  for(int xx=0; xx<tofmu; xx++){
    sprintf(hname,"tofmuch=[%d-%d]",tofmult[xx],tofmult[xx+1]);

    if(xx < 6){
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
    

    TGraphErrors *Rhist = new TGraphErrors(count[xx],XX[xx],PR[xx],0,ER[xx]);
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

    int NMM1 =0,NMM2 =0;
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



    TF1* fx = new TF1("fx", "gaus" ,  NMM1, NMM2);
    fx->SetLineColor(kBlue) ;
    fx->SetLineWidth(4);
    Rhist->Fit("fx","R");

    Sig[xx] = fx ->GetParameter(2);
    Ma[xx]  = fx ->GetParameter(1);

     CutMax[xx] =  4.*(fx ->GetParameter(2)) + (fx ->GetParameter(1));
     CutMin[xx] = -4.*(fx ->GetParameter(2)) + (fx ->GetParameter(1));

  }
  

  CutMax[tofmu] =  504.067;
  CutMin[tofmu] =  350.089;
  XXX[tofmu]    =  520;

//  CutMax[tofmu+1] =  725.458;
//  CutMin[tofmu+1] =  500.032;
//  XXX[tofmu+1]    =  3800;
  
  TCanvas * cx8 = new TCanvas("v2", "v2", 1500, 500);
  cx8->Divide(1);

  TGraphErrors *Rhistxx = new TGraphErrors(tofmu+1,XXX,CutMax,0,0);
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

  TF1* fx1 = new TF1("fx1","[0] + [1]*pow(x,1) + [2]*pow(x,2) + [3]*pow(x,3) + [4]*pow(x,4) + [5]*pow(x,5) ",  0, 550);
  fx1->SetLineColor(kBlue) ;
  fx1->SetLineWidth(4);
  Rhistxx->Fit("fx1","M S E R","0");




  TCanvas * cx9 = new TCanvas("v9", "v9", 1500, 500);
  cx9->Divide(1);
  cx9->cd(1);
  TGraphErrors *Rhistxy = new TGraphErrors(tofmu+1,XXX,CutMin,0,0);
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


  for(int xx=0; xx<tofmu; xx++){
    cout<<xx<<"  "<<XXX[xx]<<" "<<CutMax[xx]<<" "<<CutMin[xx]<<"\n";
  }
  


return 0;
}
//***********************
//=======================
