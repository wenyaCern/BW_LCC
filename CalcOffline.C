#include "TFile.h"
#include "TList.h"
#include "TLine.h"
#include "TH1D.h"
#include "TProfile.h"
#include "TGraphErrors.h"
#include "TDirectoryFile.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TF1.h"

 Double_t StraightFit(Double_t *x, Double_t *par) 
 {
   double fit=par[0]*x[0]+par[1];
   return fit; 
}

void CalcOffline(int cent) {
  TFile* inputFile = TFile::Open(Form("out_Cent%i_formal.root",cent));

  // CMW
  TProfile* pv2posAch[3] = {NULL};
  TProfile* pv2negAch[3] = {NULL};
  pv2posAch[0] = (TProfile*)inputFile->Get("pv2posAch_AftLCC");
  pv2posAch[1] = (TProfile*)inputFile->Get("pv2negAch_AftLCC");
  pv2posAch[2] = (TProfile*)inputFile->Get("pv2posAch_AftCME");
  pv2negAch[0] = (TProfile*)inputFile->Get("pv2negAch_AftLCC");
  pv2negAch[1] = (TProfile*)inputFile->Get("pv2posAch_AftCME");
  pv2negAch[2] = (TProfile*)inputFile->Get("pv2negAch_AftCMW");

  pv2posAch[0]->SetMarkerStyle(20);
  pv2posAch[0]->SetMarkerSize(0.8);
  pv2posAch[0]->SetMarkerColor(kRed);
  pv2posAch[0]->SetLineColor(kRed);
  pv2negAch[0]->SetMarkerStyle(20);
  pv2negAch[0]->SetMarkerSize(0.8);
  pv2negAch[0]->SetMarkerColor(kBlue);
  pv2negAch[0]->SetLineColor(kBlue);

  TLegend* lgResults = new TLegend(0.32,0.71,0.68,0.91);
  lgResults->SetFillStyle(0);
  lgResults->SetBorderSize(0);
  lgResults->SetTextSize(0.04);
  lgResults->AddEntry(pv2posAch[0],"h^{+}","pl"); 
  lgResults->AddEntry(pv2negAch[0],"h^{-}","pl"); 

  TCanvas* cLinear = new TCanvas();
  cLinear->cd();
  pv2posAch[0]->GetXaxis()->SetRangeUser(-0.1, 0.1);
  pv2posAch[0]->GetXaxis()->SetTitle("A_{ch}");
  pv2posAch[0]->GetYaxis()->SetTitle("v_{2}");
  gStyle->SetOptStat(0);
  pv2posAch[0]->Draw("p");
  pv2negAch[0]->Draw("same p");
  lgResults->Draw();
  
  TGraphErrors* g_dv2_ach[3]={NULL};
  TGraphErrors* g_dv2_ach_norm[3]={NULL};
  TF1* f1 = new TF1("f1", "pol1", -0.1, 0.1);

    for (int s=0; s<1; ++s){
      double x[91]={0}, y[91]={0}, yErr[91]={0}, yNorm[91]={0}, yNormErr[91]={0};
      for (int i = 0; i < 91; ++i) {
        x[i]=pv2posAch[s]->GetXaxis()->GetBinCenter(i+1);
        int xBin=pv2posAch[s]->FindBin(x[i]);

        double v2Pos = pv2posAch[s]->GetBinContent(i+1);
        double v2Neg = pv2negAch[s]->GetBinContent(i+1);
        if (v2Pos==0 || v2Neg==0) continue;
        double v2PosErr = pow(pv2posAch[s]->GetBinError(i+1), 2);
        double v2NegErr = pow(pv2negAch[s]->GetBinError(i+1), 2);

        y[i]=v2Neg - v2Pos;
        yErr[i]=sqrt(pow(v2PosErr,2)+pow(v2NegErr,2));


        yNorm[i]=2*(v2Neg - v2Pos)/(v2Neg + v2Pos);
        yNormErr[i]=sqrt( 4*(v2Neg*v2Neg + v2Neg*v2Neg)/pow(v2Neg + v2Pos, 4) *v2PosErr + 
                             4*(v2Neg*v2Neg + v2Neg*v2Neg)/pow(v2Neg + v2Pos, 4) *v2NegErr );
        cout<<x[i]<<"     "<<yNorm[i]<<endl;
      }
      g_dv2_ach[s] = new TGraphErrors(91,x,y,0,yErr);
      if (s==0) g_dv2_ach[s]->SetNameTitle("g_dv2_ach_AftLCC");
      if (s==1) g_dv2_ach[s]->SetNameTitle("g_dv2_ach_AftCME");
      if (s==2) g_dv2_ach[s]->SetNameTitle("g_dv2_ach_AftCMW");
      f1 = new TF1("f1", "pol1", -0.1, 0.1);
      // g_dv2_ach[s]->Fit("f1","R+");

      g_dv2_ach_norm[s] = new TGraphErrors(91,x,yNorm,0,yNormErr);
      if (s==0) g_dv2_ach_norm[s]->SetNameTitle("g_dv2_ach_AftLCC_Norm");
      if (s==1) g_dv2_ach_norm[s]->SetNameTitle("g_dv2_ach_AftCME_Norm");
      if (s==2) g_dv2_ach_norm[s]->SetNameTitle("g_dv2_ach_AftCMW_Norm");
      f1 = new TF1("f1", "pol1", -0.1, 0.1);
      g_dv2_ach_norm[s]->Fit("f1","R+");      
    } // end loop s

    TCanvas* cdv2ach = new TCanvas();
    cdv2ach->cd();
    g_dv2_ach_norm[0]->SetMarkerStyle(28);
    g_dv2_ach_norm[0]->GetXaxis()->SetRangeUser(-0.1, 0.1);
    g_dv2_ach_norm[0]->GetXaxis()->SetTitle("A_{ch}");
    g_dv2_ach_norm[0]->GetYaxis()->SetTitle("#Delta v_{2}");
    g_dv2_ach_norm[0]->Draw("ap");

    //CME
    const Double_t CenBins[] = {0.,5.,10.,20.,30.,40.,50.,60.,70.,80.,90.};
    TProfile * fDnnPro[2][3]={NULL};
    TH1F     * fDnnHist[2][3]={NULL};
    fDnnPro[0][0] = (TProfile*)inputFile->Get("fDnnProSS_AftLCC");
    fDnnPro[1][0] = (TProfile*)inputFile->Get("fDnnProOS_AftLCC");
    if (!fDnnPro[0][0] || !fDnnPro[1][0]){
      cout<<"this root file doesn't contain TProfile hists for DD11 and G112"<<endl;
      return;
    }
    fDnnHist[0][0] = new TH1F("fDnnHistSS_AftLCC","fDnnHistSS_AftLCC",10,CenBins); fDnnHist[0][0]->Sumw2();
    fDnnHist[1][0] = new TH1F("fDnnHistOS_AftLCC","fDnnHistOS_AftLCC",10,CenBins); fDnnHist[1][0]->Sumw2();

    fDnnPro[0][1]= (TProfile*)inputFile->Get("fDnnProSS_AftCME");
    fDnnPro[1][1]= (TProfile*)inputFile->Get("fDnnProOS_AftCME");
    fDnnHist[0][1] = new TH1F("fDnnHistSS_AftCME","fDnnHistSS_AftCME",10,CenBins); fDnnHist[0][1]->Sumw2();
    fDnnHist[1][1] = new TH1F("fDnnHistOS_AftCME","fDnnHistOS_AftCME",10,CenBins); fDnnHist[1][1]->Sumw2();

    fDnnPro[0][2] = (TProfile*)inputFile->Get("fDnnProSS_AftCMW");
    fDnnPro[1][2] = (TProfile*)inputFile->Get("fDnnProOS_AftCMW");
    fDnnHist[0][2] = new TH1F("fDnnHistSS_AftCMW","fDnnHistSS_AftCMW",10,CenBins); fDnnHist[0][2]->Sumw2();
    fDnnHist[1][2] = new TH1F("fDnnHistOS_AftCMW","fDnnHistOS_AftCMW",10,CenBins); fDnnHist[1][2]->Sumw2();

    TProfile * fCMEPro[2][3]={NULL};
    TH1F     * fCMEHist[2][3]={NULL};
  
    fCMEPro[0][0]  = (TProfile*)inputFile->Get("fCMEProSS_AftLCC");
    fCMEPro[1][0]  = (TProfile*)inputFile->Get("fCMEProOS_AftLCC");
    fCMEHist[0][0] = new TH1F("fCMEHistDdeltaSS_AftLCC","fCMEHistDdeltaSS_AftLCC",10,CenBins); fCMEHist[0][0]->Sumw2();
    fCMEHist[1][0] = new TH1F("fCMEHistDdeltaOS_AftLCC","fCMEHistDdeltaOS_AftLCC",10,CenBins); fCMEHist[1][0]->Sumw2();  
   
    fCMEPro[0][1]  = (TProfile*)inputFile->Get("fCMEProSS_AftCME");
    fCMEPro[1][1]  = (TProfile*)inputFile->Get("fCMEProOS_AftCME");
    fCMEHist[0][1] = new TH1F("fCMEHistDdeltaSS_AftCME","fCMEHistDdeltaSS_AftCME",10,CenBins); fCMEHist[0][1]->Sumw2();
    fCMEHist[1][1] = new TH1F("fCMEHistDdeltaOS_AftCME","fCMEHistDdeltaOS_AftCME",10,CenBins); fCMEHist[1][1]->Sumw2();  

    fCMEPro[0][2]  = (TProfile*)inputFile->Get("fCMEProSS_AftCMW");
    fCMEPro[1][2]  = (TProfile*)inputFile->Get("fCMEProOS_AftCMW");
    fCMEHist[0][2] = new TH1F("fCMEHistDdeltaSS_AftCMW","fCMEHistDdeltaSS_AftCMW",10,CenBins); fCMEHist[0][2]->Sumw2();
    fCMEHist[1][2] = new TH1F("fCMEHistDdeltaOS_AftCMW","fCMEHistDdeltaOS_AftCMW",10,CenBins); fCMEHist[1][2]->Sumw2();  

    TH1F* DeltaG112[3]={NULL};
    TH1F* DeltaD11[3]={NULL};
    for(Int_t s=0; s<1; ++s){
      for(Int_t k=0;k<2;k++) {
        for(Int_t c=1;c<=fCMEPro[k][s]->GetNbinsX();c++) { //centrality loop
          Double_t stats[6]={0.};
          fCMEPro[k][s]->GetXaxis()->SetRange(c,c);
          fCMEPro[k][s]->GetStats(stats);
          Double_t SumWeig    = stats[0];
          Double_t SumWeigSq  = stats[1];
          Double_t SumTwo     = stats[4];
          Double_t SumTwoSq   = stats[5];
          if(!SumWeig) continue;
          Double_t Corr   = SumTwo/SumWeig;
          Double_t SqCorr = SumTwoSq/SumWeig;
          Double_t Weig   = SumWeig;
          Double_t SqWeig = SumWeigSq;
          Double_t spread=0., termA=0., termB=0.;
          if(SqCorr-pow(Corr,2.)>=0.) { spread = pow(SqCorr-pow(Corr,2.),0.5); }
          if(TMath::Abs(Weig)>0.) { termA = (pow(SqWeig,0.5)/Weig); }
          if(1.-pow(termA,2.)>0.) { termB = 1./pow(1.-pow(termA,2.),0.5); }
          Double_t CorrErr = termA*spread*termB; // final error (unbiased estimator for standard deviation)

          fCMEHist[k][s]->SetBinContent(c,Corr);
          fCMEHist[k][s]->SetBinError(c,CorrErr);
        } // end loop c : Cent
        fCMEPro[k][s]->GetXaxis()->SetRange(1,-1);
    } // end loop k : SO & SS 
      
    //----- Get (OS - SS) for \Gamma correlator:
    if (s==0) DeltaG112[s] = (TH1F*)(fCMEHist[1][s]->Clone("DeltaG112_AftLCC"));
    if (s==1) DeltaG112[s] = (TH1F*)(fCMEHist[1][s]->Clone("DeltaG112_AftCME"));
    if (s==2) DeltaG112[s] = (TH1F*)(fCMEHist[1][s]->Clone("DeltaG112_AftCMW"));
    DeltaG112[s]->Add(fCMEHist[0][s],-1);
    for (Int_t c=1;c<=DeltaG112[s]->GetNbinsX();c++ ){
      if (DeltaG112[s]->GetBinContent(c)!=0) cout<<"DG11_cent"<<c-1<<"   :   "<<DeltaG112[s]->GetBinContent(c)<<" #pm "<<DeltaG112[s]->GetBinError(c)<<endl;
    };

    // finalise Dnn corr
    for(Int_t k=0;k<2;k++) {
      for(Int_t c=1;c<=fDnnPro[k][s]->GetNbinsX();c++) { //centrality loop
      
        Double_t stats[6]={0.};
      
        fDnnPro[k][s]->GetXaxis()->SetRange(c,c);
        fDnnPro[k][s]->GetStats(stats);
      
        Double_t SumWeig    = stats[0];
        Double_t SumWeigSq  = stats[1];
        Double_t SumTwo     = stats[4];
        Double_t SumTwoSq   = stats[5];
      
        if(!SumWeig) continue;
      
        Double_t Corr   = SumTwo/SumWeig;
        Double_t SqCorr = SumTwoSq/SumWeig;
        Double_t Weig   = SumWeig;
        Double_t SqWeig = SumWeigSq;
        Double_t spread=0., termA=0., termB=0.;
        if(SqCorr-pow(Corr,2.)>=0.) { spread = pow(SqCorr-pow(Corr,2.),0.5); }
        if(TMath::Abs(Weig)>0.) { termA = (pow(SqWeig,0.5)/Weig); }
        if(1.-pow(termA,2.)>0.) { termB = 1./pow(1.-pow(termA,2.),0.5); }
        Double_t CorrErr = termA*spread*termB; // final error (unbiased estimator for standard deviation)

        fDnnHist[k][s]->SetBinContent(c,Corr);
        fDnnHist[k][s]->SetBinError(c,CorrErr);

        } // end loop c
        fDnnPro[k][s]->GetXaxis()->SetRange(1,-1);
      } // end loop k
  
      //------ OS-SS for \Delta_NN ----------
      if (s==0) DeltaD11[s] = (TH1F*)(fDnnHist[1][s]->Clone("DeltaD11_AftLCC"));
      if (s==1) DeltaD11[s] = (TH1F*)(fDnnHist[1][s]->Clone("DeltaD11_AftCME"));
      if (s==2) DeltaD11[s] = (TH1F*)(fDnnHist[1][s]->Clone("DeltaD11_AftCMW"));
      DeltaD11[s]->Add(fDnnHist[0][s],-1);
      for (Int_t c=1;c<=DeltaD11[s]->GetNbinsX();c++ ){
        if (DeltaD11[s]->GetBinContent(c)!=0) cout<<"DD11_cent"<<c-1<<"   :   "<<DeltaD11[s]->GetBinContent(c)<<" #pm "<< DeltaD11[s]->GetBinError(c)<<endl;
      };
   } // end loop s
}


