double RadSys(double rads=0,double phis=0,double Rx=0);
double BosPhi(double phis=0,double Rx=0);
bool calcTypBF=false;
bool calcTypCME=true;
bool calcTypCMW=false;
bool calcOnline=false;
TRandom3 *gRandom = new TRandom3(0);
gRandom->SetSeed(gRandom->Rndm()*100);
const Double_t CenBins[] = {0.,5.,10.,20.,30.,40.,50.,60.,70.,80.,90.};

void simpleBW(int nEvt=1E4, Long64_t seed=1)
{ 
  //---------------------------
  // Paras preparation
  //---------------------------
  const double fMass = 0.1395; // Pion+/-
  double fracLCC[8] = {0.706,0.624,0.58,0.56, 0.54, 0.48, 0.47,0.464}; // LCC fraction for only use LCC to explain Exp data.
  // double fracLCC[8] = {0.600,0.594,0.54,0.52, 0.49, 0.46, 0.45,0.444}; // LCC fraction for explainning exp data together with CMW/CME signal.

  // BW para.
  double rho0Val[9] = {1.26163, 1.2666, 1.25422, 1.22602, 1.19576, 1.14848, 1.08652,0.9936,0.90};
  double rho0Err[9] = {0.010439,0.01011,0.01396,0.01366,0.01115,0.01017,0.009513,0.0114,0.0};
  // double rho2Val[9] = {0.0201526, 0.0319807, 0.0449495, 0.0588247, 0.068093, 0.0698467, 0.0651333,0.0561,0.0}; 
  // Lower rhon2 (old) --> cannot fit v2(Exp) exactly. But we still use the lower rhon2 with adding the CME/CMW signal.
  double rho2Val[9] = {0.0538, 0.063, 0.11, 0.135, 0.15, 0.145, 0.121,0.115,0.0}; // Higher rhon2 to reach up to the v2 in Exp
  double rho2Err[9] = {0.00087638,0.0014616,0.00294638,0.00364378,0.0032272, 0.0027959,0.00233696,0.00245, 0.0};
  double rho4Val[9] = {0.0220175, 0.0203032, 0.0187454, 0.0193533, 0.0200887, 0.0201813, 0.0195, 0.0195, 0.0};
  double rho4Err[9] = {0.0015084,0.00169984,0.00231554,0.00274814,0.00302617,0.00404366, 0.00404,0.004, 0.0};
  double RxVal[9]   = {9.56096, 9.33685, 9.05065, 8.72321, 8.44751, 8.2290, 8.07123, 7.862, 7.000}; // Ry=10.0fm fixed (default value) BW fit.
  double RxErr[9]   = {0.01183, 0.01850, 0.03768, 0.04764, 0.04423, 0.0429, 0.04033, 0.0616, 0.000};
  //rho3 set by hand because BW does not fit v3 simulataneously:
  Double_t rho3Val[9] = {0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05}; // Used same value for all cent. 
  Double_t rho3Err = 0.002; 

  // Multiplicity determination
  Double_t KinCutRatio[9] = {0.43507719, 0.43209614, 0.431945, 0.428783, 0.4336, 0.428906, 0.427373, 0.423507, 0.423523}; // Ratio = multi(after pt,eta cut) / multi(before pt,eta cut)
  Double_t mu5TeV[9] = {2205+74, 1791+61, 1290+50, 851+50, 583.3+20, 359.5+10, 222.9, 106.4+10, 49.07}; // Read raw multiplicity value from ALICE; Also supported by PRL 106, 032301 (2011) 
  // Double_t mu5TeV[9] = {2125, 1780, 1290, 841, 563.3, 350.5, 211, 96.4, 49.07}; // Tuned multiplicity for add CMW/CME signal
  Double_t gMultPre[9];
  for (int i=0; i<11; ++i) gMultPre[i] = mu5TeV[i] / KinCutRatio[i];   

  // NBD distribution of multiplicity
  Double_t NBDLow[9] = {1600, 1300, 900, 500, 400, 300, 0, 0, 0};
  Double_t NBDHigh[9] = {5000, 5000, 3000, 2000, 1200, 1000, 800, 0, 0}; // Low/High cut should be appropriate
  Double_t NStdDev[9] = {300, 300, 180, 124, 87, 63, 47, 28, 17}; // Read roughly from ALICE Pb-Pb 5TeV
  TF1* fNBD[9];
  for (int iCent = 0; iCent < 9; ++iCent){
      double N_Mean  = gMultPre[iCent]/(1.0+fracLCC[iCent]); //  Keep multiplicity to be around mu5TeV[iCent] in different fLCC
      double N_Sigma = TMath::Power(NStdDev[iCent], 2);
      double pNBD = N_Mean / N_Sigma;
      double rNBD = TMath::Power(N_Mean,2) / (N_Sigma - N_Mean);
      fNBD[iCent] = new TF1(Form("NBD%i", iCent),"TMath::Binomial([0]+x-1, x)*TMath::Power(1-[1], x)*TMath::Power([1], [0])", NBDLow[iCent], NBDHigh[iCent]);
      fNBD[iCent]->SetParameter(0, rNBD);
      fNBD[iCent]->SetParameter(1, pNBD);
  };

  // Boltz distribution for thermal
  TF1 *fBoltz1[9];
  double TempVal[9] = {91.3442,86.9569,84.7781,87.3677,91.6305,95.1377,98.1375,108.2,110}; // MeV
  double TempErr[9] = {3.5084, 3.49799, 4.8779, 4.7532, 3.8270, 3.3291, 2.8827, 3.2, 5.0};
  for(int i=0; i<9; i++) {
    TempVal[i] += 20;  // Tuned for v2 and CME observables
    // Note: TempVal is used together with the rhon2. Choos lower rhon2, we apply TempVal without "+=20";
    //           If choose higher rhon2, the TempVal should be added by "+=20".
    TempVal[i] = TempVal[i]/1000.0; // GeV
    TempErr[i] = TempErr[i]/1000.0;
    fBoltz1[i] = new TF1("Boltz1","x*TMath::Sqrt(x*x-[0]*[0])*TMath::Exp(-x/[1])",0.15,10);
    fBoltz1[i]->SetParameter(0,fMass);
    fBoltz1[i]->SetParameter(1,TempVal[i]);
    fBoltz1[i]->SetLineColor(2);
  }

  //----------------
  // QA hists
  //----------------
  TH2D* hXYEmisPoint = new TH2D("hXYEmisPoint","",100,-2,2,100,-2,2); // See charge distribution on XY
  TProfile2D* pChargeXYEmisPoint[10] = {NULL};
  pChargeXYEmisPoint[0] = new TProfile2D("pChargeXYEmisPoint0","",100,-2,2,100,-2,2); // raw
  pChargeXYEmisPoint[1] = new TProfile2D("pChargeXYEmisPoint1","",100,-2,2,100,-2,2); // cme1
  pChargeXYEmisPoint[2] = new TProfile2D("pChargeXYEmisPoint2","",100,-2,2,100,-2,2); // cme2
  pChargeXYEmisPoint[3] = new TProfile2D("pChargeXYEmisPoint3","",100,-2,2,100,-2,2); // cmw1
  pChargeXYEmisPoint[4] = new TProfile2D("pChargeXYEmisPoint4","",100,-2,2,100,-2,2); // cmw2
  pChargeXYEmisPoint[5] = new TProfile2D("pChargeXYEmisPoint5","",100,-2,2,100,-2,2); // cme1+cmw1
  pChargeXYEmisPoint[6] = new TProfile2D("pChargeXYEmisPoint6","",100,-2,2,100,-2,2); // cme1+cmw2
  pChargeXYEmisPoint[7] = new TProfile2D("pChargeXYEmisPoint7","",100,-2,2,100,-2,2); // cme2+cmw1
  pChargeXYEmisPoint[8] = new TProfile2D("pChargeXYEmisPoint8","",100,-2,2,100,-2,2); // cme2+cmw2
   
  TH1D     *fMultEvt[3]; // Multiplicity distribution
  fMultEvt[0]= new TH1D("fMult_AftLCC","fMult_AftLCC",1000,0, 5000);
  fMultEvt[1]= new TH1D("fMult_AftCME","fMult_AftCME",1000,0, 5000);
  fMultEvt[2]= new TH1D("fMult_AftCMW","fMult_AftCMW",1000,0, 5000);
  TH1D     *fLccCheck = new TH1D("fLccCheck","fLccCheck", 57,0, 1);

  TH1D* hAch = new TH1D("ach","",200,-1,1);
  TH1D* hPt    = new TH1D("pt","",100,0,7);
  TH1D* hEta  = new TH1D("eta","",100,-7,7);
  TH1D* hPhi  = new TH1D("phi","",100,-7,7);

  //-----------------------------------
  // Physics v2, pT spectrum 
  //-----------------------------------
  TProfile* v2ptPos[3];
  TProfile* v2ptNeg[3];
  TH1D*    hPtspecPos[3];
  TH1D*    hPtspecNeg[3];

  v2ptPos[0] = new TProfile("v2ptPos_AftLCC","",20,0,2);
  v2ptPos[1] = new TProfile("v2ptPos_AftCME","",20,0,2);
  v2ptPos[2] = new TProfile("v2ptPos_AftCMW","",20,0,2); 

  v2ptNeg[0] = new TProfile("v2ptNeg_AftLCC","",20,0,2);
  v2ptNeg[1] = new TProfile("v2ptNeg_AftCME","",20,0,2);
  v2ptNeg[2] = new TProfile("v2ptNeg_AftCMW","",20,0,2); 

  TProfile*  hV2IntPro = new TProfile("V2IntPro","V2IntPro",10,CenBins,"s"); hV2IntPro->Sumw2();

  double     fPtMin = 0.2;
  double     fPtMax = 5.0; 
  int            fNbinsPt =  50;
  double     fPtWd =  (double)(fPtMax - fPtMin)/fNbinsPt;
  hPtspecPos[0] = new TH1D("hPtspecPosPion_AftLCC", "", fNbinsPt, fPtMin, fPtMax);
  hPtspecPos[1] = new TH1D("hPtspecPosPion_AftCME", "", fNbinsPt, fPtMin, fPtMax);
  hPtspecPos[2] = new TH1D("hPtspecPosPion_AftCMW", "", fNbinsPt, fPtMin, fPtMax);
  hPtspecNeg[0] = new TH1D("hPtspecNegPion_AftLCC", "", fNbinsPt, fPtMin, fPtMax);
  hPtspecNeg[1] = new TH1D("hPtspecNegPion_AftCME", "", fNbinsPt, fPtMin, fPtMax);
  hPtspecNeg[2] = new TH1D("hPtspecNegPion_AftCMW", "", fNbinsPt, fPtMin, fPtMax);

  //----------------------------------------------
  // Physics CME, 3 particle correlator
  //----------------------------------------------
  TProfile * fCMEPro[2][3];
  TH1F     * fCMEHist[2][3];

  fCMEPro[0][0]  = new TProfile("fCMEProSS_AftLCC","fCMEProSS_AftLCC",10,CenBins,"s"); fCMEPro[0][0]->Sumw2();
  fCMEPro[1][0]  = new TProfile("fCMEProOS_AftLCC","fCMEProOS_AftLCC",10,CenBins,"s"); fCMEPro[1][0]->Sumw2();
  fCMEHist[0][0] = new TH1F("fCMEHistDdeltaSS_AftLCC","fCMEHistDdeltaSS_AftLCC",10,CenBins); fCMEHist[0][0]->Sumw2();
  fCMEHist[1][0] = new TH1F("fCMEHistDdeltaOS_AftLCC","fCMEHistDdeltaOS_AftLCC",10,CenBins); fCMEHist[1][0]->Sumw2();  

  fCMEPro[0][1]  = new TProfile("fCMEProSS_AftCME","fCMEProSS_AftCME",10,CenBins,"s"); fCMEPro[0][1]->Sumw2();
  fCMEPro[1][1]  = new TProfile("fCMEProOS_AftCME","fCMEProOS_AftCME",10,CenBins,"s"); fCMEPro[1][1]->Sumw2();
  fCMEHist[0][1] = new TH1F("fCMEHistDdeltaSS_AftCME","fCMEHistDdeltaSS_AftCME",10,CenBins); fCMEHist[0][1]->Sumw2();
  fCMEHist[1][1] = new TH1F("fCMEHistDdeltaOS_AftCME","fCMEHistDdeltaOS_AftCME",10,CenBins); fCMEHist[1][1]->Sumw2();  

  fCMEPro[0][2]  = new TProfile("fCMEProSS_AftCMW","fCMEProSS_AftCMW",10,CenBins,"s"); fCMEPro[0][2]->Sumw2();
  fCMEPro[1][2]  = new TProfile("fCMEProOS_AftCMW","fCMEProOS_AftCMW",10,CenBins,"s"); fCMEPro[1][2]->Sumw2();
  fCMEHist[0][2] = new TH1F("fCMEHistDdeltaSS_AftCMW","fCMEHistDdeltaSS_AftCMW",10,CenBins); fCMEHist[0][2]->Sumw2();
  fCMEHist[1][2] = new TH1F("fCMEHistDdeltaOS_AftCMW","fCMEHistDdeltaOS_AftCMW",10,CenBins); fCMEHist[1][2]->Sumw2();  

  //----------------------------------------------
  // Physics CME, 2 particle correlator
  //----------------------------------------------
  TProfile * fDnnPro[2][3];
  TH1F     * fDnnHist[2][3];

  fDnnPro[0][0] = new TProfile("fDnnProSS_AftLCC","fDnnProSS_AftLCC",10,CenBins,"s"); fDnnPro[0][0]->Sumw2();
  fDnnPro[1][0] = new TProfile("fDnnProOS_AftLCC","fDnnProOS_AftLCC",10,CenBins,"s"); fDnnPro[1][0]->Sumw2();
  fDnnHist[0][0] = new TH1F("fDnnHistSS_AftLCC","fDnnHistSS_AftLCC",10,CenBins); fDnnHist[0][0]->Sumw2();
  fDnnHist[1][0] = new TH1F("fDnnHistOS_AftLCC","fDnnHistOS_AftLCC",10,CenBins); fDnnHist[1][0]->Sumw2();

  fDnnPro[0][1]= new TProfile("fDnnProSS_AftCME","fDnnProSS_AftCME",10,CenBins,"s"); fDnnPro[0][1]->Sumw2();
  fDnnPro[1][1]= new TProfile("fDnnProOS_AftCME","fDnnProOS_AftCME",10,CenBins,"s"); fDnnPro[1][1]->Sumw2();
  fDnnHist[0][1] = new TH1F("fDnnHistSS_AftCME","fDnnHistSS_AftCME",10,CenBins); fDnnHist[0][1]->Sumw2();
  fDnnHist[1][1] = new TH1F("fDnnHistOS_AftCME","fDnnHistOS_AftCME",10,CenBins); fDnnHist[1][1]->Sumw2();

  fDnnPro[0][2] = new TProfile("fDnnProSS_AftCMW","fDnnProSS_AftCMW",10,CenBins,"s"); fDnnPro[0][2]->Sumw2();
  fDnnPro[1][2] = new TProfile("fDnnProOS_AftCMW","fDnnProOS_AftCMW",10,CenBins,"s"); fDnnPro[1][2]->Sumw2();
  fDnnHist[0][2] = new TH1F("fDnnHistSS_AftCMW","fDnnHistSS_AftCMW",10,CenBins); fDnnHist[0][2]->Sumw2();
  fDnnHist[1][2] = new TH1F("fDnnHistOS_AftCMW","fDnnHistOS_AftCMW",10,CenBins); fDnnHist[1][2]->Sumw2();

  //--------------------
  // Physics CMW
  //--------------------
  TProfile* pv2posAch[3] = {NULL};
  TProfile* pv2negAch[3] = {NULL};
  pv2posAch[0] = new TProfile("pv2posAch_AftLCC","",91,-1,1);
  pv2posAch[1] = new TProfile("pv2posAch_AftCME","",91,-1,1);
  pv2posAch[2] = new TProfile("pv2posAch_AftCMW","",91,-1,1);
  pv2negAch[0] = new TProfile("pv2negAch_AftLCC","",91,-1,1);
  pv2negAch[1] = new TProfile("pv2negAch_AftCME","",91,-1,1);
  pv2negAch[2] = new TProfile("pv2negAch_AftCMW","",91,-1,1);

  TProfile* pAch_nSwitchCMW = new TProfile("ach_nSwitchCMW","",91,-1,1); // Switch time depends on Ach  

  //----------------------------------
  // Physics Balance function
  //----------------------------------
  double EtaMax =0.8;
  double dEtaMin[4] = {0,0,0,0};
  double dEtaMax[4] = {2*EtaMax, 2*EtaMax, 2*EtaMax, 2*EtaMax};
  int Nbins[4] = {31,31,31,31};

  THnSparseD* cPP[3] = {NULL}; 
  THnSparseD* fPP[3] = {NULL};  // dim1: PP; dim2: PN; dim3: NP; dim4: NN;
  cPP[0] = new THnSparseD("cValBF_AftLCC","cValBF_AftLCC", 4, Nbins, dEtaMin, dEtaMax); cPP[0]->Sumw2();
  fPP[0]  = new THnSparseD("fValBF_AftLCC","fValBF_AftLCC", 4, Nbins, dEtaMin, dEtaMax);   fPP[0]->Sumw2();
  cPP[1] = new THnSparseD("cValBF_AftCME","cValBF_AftCME", 4, Nbins, dEtaMin, dEtaMax); cPP[1]->Sumw2();
  fPP[1]  = new THnSparseD("fValBF_AftCME","fValBF_AftCME", 4, Nbins, dEtaMin, dEtaMax);   fPP[1]->Sumw2();
  cPP[2] = new THnSparseD("cValBF_AftCMW","cValBF_AftCMW", 4, Nbins, dEtaMin, dEtaMax);cPP[2]->Sumw2();
  fPP[2]  = new THnSparseD("fValBF_AftCMW","fValBF_AftCMW", 4, Nbins, dEtaMin, dEtaMax);  fPP[2]->Sumw2();

  TH1D* BF0[3];
  BF0[0] = new TH1D("BalaFunc0_AftLCC","", 31, 0, 2*EtaMax);
  BF0[1] = new TH1D("BalaFunc0_AftCME","", 31, 0, 2*EtaMax);
  BF0[2] = new TH1D("BalaFunc0_AftCMW","", 31, 0, 2*EtaMax);

  TH1D* hTotEvtMQ = new TH1D("hTotEvtMQ", "", 10, 0, 10);
  hTotEvtMQ->GetXaxis()->SetBinLabel(1,"hTotEvtMQPos_AftLCC");
  hTotEvtMQ->GetXaxis()->SetBinLabel(2,"hTotEvtMQNeg_AftLCC");
  hTotEvtMQ->GetXaxis()->SetBinLabel(3,"hTotEvtMQPos_AftCME");
  hTotEvtMQ->GetXaxis()->SetBinLabel(4,"hTotEvtMQNeg_AftCME");
  hTotEvtMQ->GetXaxis()->SetBinLabel(5,"hTotEvtMQPos_AftCMW");
  hTotEvtMQ->GetXaxis()->SetBinLabel(6,"hTotEvtMQNeg_AftCMW");

  //----------------
  // Loop Evt
  //----------------
  const int NTRKMAX = 5000;
  double TotEvtMQP[3] ={0};
  double TotEvtMQN[3] ={0};

  // Start loop event
  for (int iEvt = 0; iEvt < nEvt; ++iEvt) {
    if (iEvt%10000==0) cout<<"#loop event "<<iEvt<<endl;
    int         centBin = -1;
    double  centrality = 30.5; 
    if(centrality <= 5)                          centBin = 0;
    else if(centrality > 5 &&  centrality <= 10) centBin = 1;    
    else if(centrality > 10 && centrality <= 20) centBin = 2;
    else if(centrality > 20 && centrality <= 30) centBin = 3;
    else if(centrality > 30 && centrality <= 40) centBin = 4;
    else if(centrality > 40 && centrality <= 50) centBin = 5;
    else if(centrality > 50 && centrality <= 60) centBin = 6;
    else if(centrality > 60 && centrality <= 70) centBin = 7;
    
    double Rx = RxVal[centBin]/10.;
    double rho0 = rho0Val[centBin];
    double rho2 = rho2Val[centBin];
    double rho3 = rho3Val[centBin];
    double rho4 = rho4Val[centBin];

    int        nPos=0, nNeg=0, nTrk=0, nPosLowPt=0, nNegLowPt=0;
    double sumCos2phi[6]={0}, sumCosphi[6]={0}, sumSinphi[6]={0}, sumCos2phiLowPt[6]={0};
    double pt[NTRKMAX]={0}, phi[NTRKMAX]={0}, x[NTRKMAX]={0}, y[NTRKMAX]={0}, eta[NTRKMAX]={0};
    int        charge[NTRKMAX]={0}, isTrkLCC[NTRKMAX]={0};

    Int_t nSourceNum = (Int_t)fNBD[centBin]->GetRandom();
    // Loop spatial point
    for (int iSpatialPoint=0; iSpatialPoint<nSourceNum; ++iSpatialPoint) {
      double rads = sqrt(gRandom->Rndm());
      double phis = 2*TMath::Pi()*gRandom->Rndm();
      while (RadSys(rads,phis,Rx)>1.) {
        rads = sqrt(gRandom->Rndm());
        phis = 2*TMath::Pi()*gRandom->Rndm();
      }
      double xThisPoint = rads*cos(phis);
      double yThisPoint = rads*sin(phis);
      hXYEmisPoint->Fill(xThisPoint,yThisPoint);
      
      double csthetas = 2.*(gRandom->Rndm()-0.5);
      double thetas = acos(csthetas);
      double etas = -log(tan(thetas/2.));
      double phib = BosPhi(phis,Rx);
      // double rhob = pow(RadSys(rads,phis,Rx),1)*(rho0+rho2*cos(2.*phib)+rho3*cos(3.*phib)+rho4*cos(4.*phib));
      double rhob = pow(RadSys(rads,phis,Rx),1)*(rho0+rho2*cos(2.*phib)); // Here I do not add rhon3 and rhon4 
      TLorentzVector bvec;
      bvec.SetPxPyPzE(sinh(rhob)*cos(phib),sinh(rhob)*sin(phib),cosh(rhob)*sinh(etas),cosh(rhob)*cosh(etas));
      TVector3 lp1b;
      lp1b = bvec.BoostVector();

      int nParticleThisPoint = -1;
      double fBkgforLCC = fracLCC[centBin];
      do {
        if (gRandom->Rndm() < fBkgforLCC) nParticleThisPoint=2;
        else nParticleThisPoint=1;  
      }while(nParticleThisPoint<0);

      // loop track at each point
      for (int iParticle=0; iParticle<nParticleThisPoint; ++iParticle) {

        double E = fBoltz1[centBin]->GetRandom();
        double p1 = sqrt(E*E-fMass*fMass);
        double cstheta = 2.*(gRandom->Rndm()-0.5);
        double theta = acos(cstheta);
        double pt1 = p1*sin(theta);
        double pz1 = p1*cstheta;
        double phi1 = 2.*TMath::Pi()*gRandom->Rndm();

        TLorentzVector lp1;
        lp1.SetPxPyPzE(pt1*cos(phi1),pt1*sin(phi1),pz1,E);
        lp1.Boost(lp1b);

        double ptThisTrk   = lp1.Pt();
        double etaThisTrk = lp1.PseudoRapidity();
        double phiThisTrk = lp1.Phi();

        int chargeThisTrk=0;
        if (nParticleThisPoint==2 && iParticle==0) chargeThisTrk=1; 
        else if (nParticleThisPoint==2 && iParticle==1) chargeThisTrk=-1;
        else if (nParticleThisPoint==1) chargeThisTrk = (gRandom->Rndm()>0.5) ? 1 : -1; 
        if (chargeThisTrk==0) continue;
        
        if (ptThisTrk<fPtMin || ptThisTrk>fPtMax) continue;
        if (fabs(etaThisTrk)>0.8) continue;

        hPt  ->Fill(ptThisTrk);
        hEta->Fill(etaThisTrk);
        hPhi->Fill(phiThisTrk);
        hV2IntPro->Fill(centrality, cos(2*phiThisTrk));

        if (chargeThisTrk>0) {
          nPos++;
          sumCos2phi[0]+=cos(2*phiThisTrk);
          sumCosphi[0]+=cos(phiThisTrk);
          sumSinphi[0]+=sin(phiThisTrk);
          v2ptPos[0]->Fill(ptThisTrk, cos(2*phiThisTrk));
          hPtspecPos[0]->Fill(ptThisTrk, 1./(ptThisTrk*fPtWd));
        } else {
          nNeg++;
          sumCos2phi[1]+=cos(2*phiThisTrk);
          sumCosphi[1]+=cos(phiThisTrk);
          sumSinphi[1]+=sin(phiThisTrk);
          v2ptNeg[0]->Fill(ptThisTrk, cos(2*phiThisTrk));
          hPtspecNeg[0]->Fill(ptThisTrk, 1./(ptThisTrk*fPtWd));
        }

        if (ptThisTrk<2.0){ // Pt cut for CMW ---> Same with ALICE 5 TeV results 
              if (chargeThisTrk>0) {
                nPosLowPt++;
                sumCos2phiLowPt[0]+=cos(2*phiThisTrk);
              } else {
                nNegLowPt++;
                sumCos2phiLowPt[1]+=cos(2*phiThisTrk);
              }          
        }

        pt[nTrk]   = ptThisTrk; 
        phi[nTrk] = phiThisTrk;
        eta[nTrk] = etaThisTrk;
        charge[nTrk] = chargeThisTrk;
        x[nTrk] = xThisPoint;
        y[nTrk] = yThisPoint;
        pChargeXYEmisPoint[0]->Fill(xThisPoint,yThisPoint,chargeThisTrk);
        if (nParticleThisPoint==2) isTrkLCC[nTrk] = 1;
        nTrk++;
      } // Track
    } // Spatial point

    double totNum = nSourceNum*fracLCC[centBin]*2+nSourceNum*(1-fracLCC[centBin]);
    // cout<<totNum<<"   "<<nTrk<<"  ratio : "<<nTrk/totNum<<endl;
    fLccCheck->Fill(nTrk/totNum);
    fMultEvt[0]->Fill(nPosLowPt+nNegLowPt);

    double ach = -999.;
    ach = (double)(nPosLowPt-nNegLowPt)/(nPosLowPt+nNegLowPt);
    hAch->Fill(ach);
    if (calcTypCMW){
        //===================Fill CMW Profile Histograms===================
        pv2posAch[0]->Fill(ach,sumCos2phiLowPt[0]/nPosLowPt);
        pv2negAch[0]->Fill(ach,sumCos2phiLowPt[1]/nNegLowPt);
        //=================== CMW filling done ==========================
    }

    double CMEDSS = -999, CMEDOS = -999., CMEGSS = -999., CMEGOS = -999.;
    if (calcTypCME){
        //===================Fill CME Profile Histo=======================
        //Delta11:
        CMEDSS = (sumCosphi[0]*sumCosphi[0]+sumSinphi[0]*sumSinphi[0]-nPos)/(nPos*nPos-nPos);
        fDnnPro[0][0]->Fill(centrality,CMEDSS,nPos*nPos-nPos);
        CMEDOS = (sumCosphi[0]*sumCosphi[1]+sumSinphi[0]*sumSinphi[1])/(nPos*nNeg);
        fDnnPro[1][0]->Fill(centrality,CMEDOS,nPos*nNeg);

        //Gamma112:
        CMEGSS = (sumCosphi[0]*sumCosphi[0]-sumSinphi[0]*sumSinphi[0]-sumCos2phi[0])/(nPos*nPos-nPos);
        fCMEPro[0][0]->Fill(centrality,CMEGSS,nPos*nPos-nPos); //Same-Sign
        CMEGOS = (sumCosphi[0]*sumCosphi[1]-sumSinphi[0]*sumSinphi[1])/(nPos*nNeg);
        fCMEPro[1][0]->Fill(centrality,CMEGOS,nPos*nNeg);     //Oppo-Sign

        //=================== CME filling done ==========================
    }

    double mPos = 0, mNeg = 0;
    if (calcTypBF){
        //===================Calc balance function loop trk again!=======================
        mPos = 0, mNeg = 0;
        for (int iTrk = 0; iTrk < nTrk; ++iTrk){

          double pt1   = pt[iTrk];
          double eta1 = eta[iTrk];
          double phi1 = phi[iTrk];
          double ch1 = charge[iTrk];
    
          if (ch1==0) break;
          if (pt1>1.5 || pt1<0.3 ||fabs(eta1)>0.8) continue; // Pt cut same with ALICE 2/5TeV

          for (int jTrk=0; jTrk < nTrk; ++jTrk){
    
            if (jTrk == iTrk) continue;
    
            double pt2 = pt[jTrk];
            double eta2 = eta[jTrk];
            double phi2 = phi[jTrk];
            double ch2 = charge[jTrk];
    
            if (ch2==0) continue;
            if (pt2>1.5 || pt2<0.3 ||fabs(eta2)>0.8) continue;
    
            double DEta = fabs(eta1 - eta2);
            double DPt = pt1 - pt2;
            double DPhi = phi1 - phi2;
    
            // Balance function
            double cpp[4] = {-999., -999., -999., -999.};
            double fpp[4] = {-999., -999., -999., -999.};
            if (ch1  > 0){
              if (ch2 > 0) cpp[0] = DEta; // PP
              if (ch2 < 0) cpp[1] = DEta; // PN
    
              double rndmCharge3 = gRandom->Rndm();
              if (rndmCharge3>0.5) fpp[0]=DEta; // PP
              else fpp[1]=DEta; // PN
            }
            else if (ch1  < 0){
              if (ch2 > 0) cpp[2] = DEta; // NP
              if (ch2 < 0) cpp[3] = DEta; // NN
    
              double rndmCharge3 = gRandom->Rndm();
              if (rndmCharge3>0.5) fpp[2]=DEta; // NP
              else fpp[3]=DEta; // NN
            }
    
            cPP[0]->Fill(cpp);
            fPP[0]->Fill(fpp);
          } // Loop Inner jTck Done
      
          if (ch1 > 0) mPos ++; 
          else mNeg ++;
        }
      
        hTotEvtMQ->Fill(0.5, mPos);
        hTotEvtMQ->Fill(1.5, mNeg);
        TotEvtMQP[0] += mPos;
        TotEvtMQN[0] += mNeg;
        //===================Calc balance function done!=======================      
    }

    //----------------
    // cme
    //----------------
    double fCMEVal[10] = {0.00152293, 0.00229137, 0.00224717, 0.00379572, 0.00430225, 0.00584469, 0.00971766, 0.00905047};
    int nSwitchTotalCME = (int)(nPos+nNeg)*fCMEVal[centBin];
    // cout<<nPos+nNeg<<" "<<nSwitchTotalCME<<endl;
    nPos=0; nNeg=0; nPosLowPt=0; nNegLowPt=0;
    int nSwitchCME = 0;
    int IsSwitchedCME[NTRKMAX]={0};
    int switchTypeCME = (gRandom->Rndm()>0.5) ? 1 : -1; // upper - lower + and upper + lower -
    for (int iTrk = 0; iTrk < NTRKMAX; ++iTrk) {
      // break; // close add CME
      if (charge[iTrk]==0) break;
      // loop for switching charge
      for (int jTrk = iTrk+1; jTrk < NTRKMAX; ++jTrk) {
        // no switch if switching time already reaches the maximum
        if (nSwitchCME>=nSwitchTotalCME) break;
        if (charge[jTrk]==0) continue;
        // no switch if iTrk or jTrk is from LCC
        if (isTrkLCC[iTrk]==1) break;
        if (isTrkLCC[jTrk]==1) continue;
        // no switch if iTrk or jTrk is already switched
        if (IsSwitchedCME[iTrk]==1) break;
        if (IsSwitchedCME[jTrk]==1) continue;
        // no switch if iTrk already matches the desired configuration
        if (switchTypeCME>0 && y[iTrk]>0 && charge[iTrk]<0) break;
        if (switchTypeCME>0 && y[iTrk]<0 && charge[iTrk]>0) break;
        if (switchTypeCME<0 && y[iTrk]>0 && charge[iTrk]>0) break;
        if (switchTypeCME<0 && y[iTrk]<0 && charge[iTrk]<0) break;
        // switch charge
        if (switchTypeCME>0) { // switchTypeCME
          if ((y[iTrk]>0 && charge[iTrk]>0 && y[jTrk]<0 && charge[jTrk]<0) ||
              (y[iTrk]<0 && charge[iTrk]<0 && y[jTrk]>0 && charge[jTrk]>0)) {
            charge[iTrk] = -charge[iTrk];
            charge[jTrk] = -charge[jTrk];
            IsSwitchedCME[iTrk]=1; IsSwitchedCME[jTrk]=1;
            nSwitchCME++;
            break;
          }
        } else {
          if ((y[iTrk]<0 && charge[iTrk]>0 && y[jTrk]>0 && charge[jTrk]<0) ||
              (y[iTrk]>0 && charge[iTrk]<0 && y[jTrk]<0 && charge[jTrk]>0)) {
            charge[iTrk] = -charge[iTrk];
            charge[jTrk] = -charge[jTrk];
            IsSwitchedCME[iTrk]=1; IsSwitchedCME[jTrk]=1;
            nSwitchCME++;
            break;
          }
        } 
      } // jTrk

      if (switchTypeCME>0) {
        pChargeXYEmisPoint[1]->Fill(x[iTrk],y[iTrk],charge[iTrk]); // CME 1 only
        pChargeXYEmisPoint[5]->Fill(x[iTrk],y[iTrk],charge[iTrk]); // CME 1 + CMW 1
        pChargeXYEmisPoint[6]->Fill(x[iTrk],y[iTrk],charge[iTrk]); // CME 1 + CMW 2
      }
      else {
        pChargeXYEmisPoint[2]->Fill(x[iTrk],y[iTrk],charge[iTrk]); // CME 2 only
        pChargeXYEmisPoint[7]->Fill(x[iTrk],y[iTrk],charge[iTrk]); // CME 2 + CMW 1
        pChargeXYEmisPoint[8]->Fill(x[iTrk],y[iTrk],charge[iTrk]); // CME 2 + CMW 2
      }

      if (charge[iTrk]>0) {
        nPos++;
        sumCos2phi[2]+=cos(2*phi[iTrk]);
        sumCosphi[2]+=cos(phi[iTrk]);
        sumSinphi[2]+=sin(phi[iTrk]);
        v2ptPos[1]->Fill(pt[iTrk], cos(2*phi[iTrk]));
        hPtspecPos[1]->Fill(pt[iTrk], 1./(pt[iTrk]*fPtWd));
      } else {
        nNeg++;
        sumCos2phi[3]+=cos(2*phi[iTrk]);
        sumCosphi[3]+=cos(phi[iTrk]);
        sumSinphi[3]+=sin(phi[iTrk]);
        v2ptNeg[1]->Fill(pt[iTrk], cos(2*phi[iTrk]));
        hPtspecNeg[1]->Fill(pt[iTrk], 1./(pt[iTrk]*fPtWd));
      }

      if (pt[iTrk]<2.0){
            if (charge[iTrk]>0) {
                nPosLowPt++;
                sumCos2phiLowPt[2]+=cos(2*phi[iTrk]);
            } else {
                nNegLowPt++;
                sumCos2phiLowPt[3]+=cos(2*phi[iTrk]);
            }          
      }

    } // iTrk

    fMultEvt[1]->Fill(nPosLowPt+nNegLowPt);

    if (calcTypCMW){    
        //===================Fill CMW Profile Histograms===================
        pv2posAch[1]->Fill(ach,sumCos2phiLowPt[2]/nPosLowPt);
        pv2negAch[1]->Fill(ach,sumCos2phiLowPt[3]/nNegLowPt);
        //=================== CMW filling done ==========================
    }

    if (calcTypCME){
        //===================Fill CME Profile Histo=======================
        //Delta11:
        CMEDSS = (sumCosphi[2]*sumCosphi[2]+sumSinphi[2]*sumSinphi[2]-nPos)/(nPos*nPos-nPos);
        fDnnPro[0][1]->Fill((double)(fCMEBin+0.5),CMEDSS,nPos*nPos-nPos);
        CMEDOS = (sumCosphi[2]*sumCosphi[3]+sumSinphi[2]*sumSinphi[3])/(nPos*nNeg);
        fDnnPro[1][1]->Fill((double)(fCMEBin+0.5),CMEDOS,nPos*nNeg);
    
        //Gamma112:
        CMEGSS = (sumCosphi[2]*sumCosphi[2]-sumSinphi[2]*sumSinphi[2]-sumCos2phi[2])/(nPos*nPos-nPos);
        fCMEPro[0][1]->Fill((double)(fCMEBin+0.5),CMEGSS,nPos*nPos-nPos); //Same-Sign
        CMEGOS = (sumCosphi[2]*sumCosphi[3]-sumSinphi[2]*sumSinphi[3])/(nPos*nNeg);
        fCMEPro[1][1]->Fill((double)(fCMEBin+0.5),CMEGOS,nPos*nNeg);     //Oppo-Sign

        //=================== CME filling done ==========================
    }

    if (calcTypBF){
        //===================Calc balance function loop trk again!=======================
        mPos = 0, mNeg = 0;
        for (int iTrk = 0; iTrk < nTrk; ++iTrk){
    
          double pt1   = pt[iTrk];
          double eta1 = eta[iTrk];
          double phi1 = phi[iTrk];
          double ch1 = charge[iTrk];
    
          if (ch1==0) break;
          if (pt1>1.5 || pt1<0.3 ||fabs(eta1)>0.8) continue;
    
          for (int jTrk=0; jTrk < nTrk; ++jTrk){
    
            if (jTrk == iTrk) continue;
    
            double pt2 = pt[jTrk];
            double eta2 = eta[jTrk];
            double phi2 = phi[jTrk];
            double ch2 = charge[jTrk];
    
            if (ch2==0) continue;
            if (pt2>1.5 || pt2<0.3 ||fabs(eta2)>0.8) continue;
    
            double DEta = fabs(eta1 - eta2);
            double DPt = pt1 - pt2;
            double DPhi = phi1 - phi2;
    
            // Balance function
            double cpp[4] = {-999.,-999.,-999.,-999.};
            double fpp[4] = {-999.,-999.,-999.,-999.};
            if (ch1  > 0){
              if (ch2 > 0) cpp[0] = DEta; // PP
              if (ch2 < 0) cpp[1] = DEta; //PN
    
              double rndmCharge3 = gRandom->Rndm();
              if (rndmCharge3>0.5) fpp[0]=DEta; // PP
              else fpp[1]=DEta; // PN
            }
            else if (ch1  < 0){
              if (ch2 > 0) cpp[2] = DEta; // NP
              if (ch2 < 0) cpp[3] = DEta; // NN
    
              double rndmCharge3 = gRandom->Rndm();
              if (rndmCharge3>0.5) fpp[2]=DEta; // NP
              else fpp[3]=DEta; // NN
            }
    
            cPP[1]->Fill(cpp);
            fPP[1]->Fill(fpp);
    
          } // Loop Inner Track Done
          
          if (ch1 > 0) mPos ++; 
          else mNeg ++;
        }
  
        hTotEvtMQ->Fill(2.5, mPos);
        hTotEvtMQ->Fill(3.5, mNeg);
        TotEvtMQP[1] += mPos;
        TotEvtMQN[1] += mNeg;
        //===================Calc balance function done!=======================
    }

    //----------------
    // cmw
    //----------------
    double fCMWVal[10] = {0.00884383, 0.0120555, 0.0141554, 0.00919913, 0.0237261, 0.0304849, 0.0331163, 0.0413417};
    int nSwitchTotalCMW = (int)fabs((nPos-nNeg)*fCMWVal[centBin]);
    // if (nSwitchTotalCMW>0) cout<<nPos<<" "<<nNeg<<" "<<nSwitchTotalCMW<<endl;
    pAch_nSwitchCMW->Fill(ach,nSwitchTotalCMW);
    nPos=0; nNeg=0;
    nPosLowPt=0; nNegLowPt=0;
    int nSwitchCMW = 0;
    int IsSwitchedCMW[NTRKMAX]={0};
    for (int iTrk = 0; iTrk < NTRKMAX; ++iTrk) {
      // break;
      if (charge[iTrk]==0) break;
      // loop for switching charge
      for (int jTrk = iTrk+1; jTrk < NTRKMAX; ++jTrk) {
        // no switch if switching time already reaches the maximum
        if (nSwitchCMW>=nSwitchTotalCMW) break;
        if (charge[jTrk]==0) break;
        // no switch if iTrk or jTrk is from LCC
        if (isTrkLCC[iTrk]==1) break;
        if (isTrkLCC[jTrk]==1) continue;
        // no switch if iTrk or jTrk is already switched
        if (IsSwitchedCME[iTrk]==1) break;
        if (IsSwitchedCMW[iTrk]==1) break;
        if (IsSwitchedCME[jTrk]==1) continue;
        if (IsSwitchedCMW[jTrk]==1) continue;
        if (ach>0) {
          if ((y[iTrk]>0 && charge[iTrk]>0 && y[jTrk]>0 && charge[jTrk]<0 && y[jTrk]>y[iTrk]) ||
              (y[iTrk]>0 && charge[iTrk]<0 && y[jTrk]>0 && charge[jTrk]>0 && y[jTrk]<y[iTrk]) ||
              (y[iTrk]<0 && charge[iTrk]>0 && y[jTrk]<0 && charge[jTrk]<0 && y[jTrk]<y[iTrk]) ||
              (y[iTrk]<0 && charge[iTrk]<0 && y[jTrk]<0 && charge[jTrk]>0 && y[jTrk]>y[iTrk])
             ) {
            charge[iTrk] = -charge[iTrk];
            charge[jTrk] = -charge[jTrk];
            IsSwitchedCMW[iTrk]=1; IsSwitchedCMW[jTrk]=1;
            nSwitchCMW++;
            break;
          }
        } else {
          if ((y[iTrk]>0 && charge[iTrk]<0 && y[jTrk]>0 && charge[jTrk]>0 && y[jTrk]>y[iTrk]) ||
              (y[iTrk]>0 && charge[iTrk]>0 && y[jTrk]>0 && charge[jTrk]<0 && y[jTrk]<y[iTrk]) ||
              (y[iTrk]<0 && charge[iTrk]<0 && y[jTrk]<0 && charge[jTrk]>0 && y[jTrk]<y[iTrk]) ||
              (y[iTrk]<0 && charge[iTrk]>0 && y[jTrk]<0 && charge[jTrk]<0 && y[jTrk]>y[iTrk])
             ) {
            charge[iTrk] = -charge[iTrk];
            charge[jTrk] = -charge[jTrk];
            IsSwitchedCMW[iTrk]=1; IsSwitchedCMW[jTrk]=1;
            nSwitchCMW++;
            break;
          }
        } // ach
      } // jTrk

      if (ach>0) {
        pChargeXYEmisPoint[3]->Fill(x[iTrk],y[iTrk],charge[iTrk]); // CMW 1 only
        pChargeXYEmisPoint[5]->Fill(x[iTrk],y[iTrk],charge[iTrk]); // CME 1 + CMW 1
        pChargeXYEmisPoint[7]->Fill(x[iTrk],y[iTrk],charge[iTrk]); // CME 2 + CMW 1
      } else {
        pChargeXYEmisPoint[4]->Fill(x[iTrk],y[iTrk],charge[iTrk]); // CMW 2 only
        pChargeXYEmisPoint[6]->Fill(x[iTrk],y[iTrk],charge[iTrk]); // CME 1 + CMW 2
        pChargeXYEmisPoint[8]->Fill(x[iTrk],y[iTrk],charge[iTrk]); // CME 2 + CMW 2
      }

      if (charge[iTrk]>0) {
        nPos++;
        sumCos2phi[4]+=cos(2*phi[iTrk]);
        sumCosphi[4]+=cos(phi[iTrk]);
        sumSinphi[4]+=sin(phi[iTrk]);
        v2ptPos[2]->Fill(pt[iTrk], cos(2*phi[iTrk]));
        hPtspecPos[2]->Fill(pt[iTrk], 1./(pt[iTrk]*fPtWd));
      } else {
        nNeg++;
        sumCos2phi[5]+=cos(2*phi[iTrk]);
        sumCosphi[5]+=cos(phi[iTrk]);
        sumSinphi[5]+=sin(phi[iTrk]);
        v2ptNeg[2]->Fill(pt[iTrk], cos(2*phi[iTrk]));
        hPtspecNeg[2]->Fill(pt[iTrk], 1./(pt[iTrk]*fPtWd));
      }

      if (pt[iTrk]<2.0){
            if (charge[iTrk]>0) {
                nPosLowPt++;
                sumCos2phiLowPt[4]+=cos(2*phi[iTrk]);
            } else {
                nNegLowPt++;
                sumCos2phiLowPt[5]+=cos(2*phi[iTrk]);
            }          
      }

    } // iTrk

    fMultEvt[2]->Fill(nPosLowPt+nNegLowPt);

    if (calcTypCMW){
        //===================Fill CMW Profile Histograms===================
        pv2posAch[2]->Fill(ach,sumCos2phiLowPt[4]/nPosLowPt);
        pv2negAch[2]->Fill(ach,sumCos2phiLowPt[5]/nNegLowPt);
        //=================== CMW filling done ==========================

    }

    if (calcTypCME){
        //===================Fill CME Profile Histo=======================
        //Delta11:
        CMEDSS = (sumCosphi[4]*sumCosphi[4]+sumSinphi[4]*sumSinphi[4]-nPos)/(nPos*nPos-nPos);
        fDnnPro[0][2]->Fill(centrality,CMEDSS,nPos*nPos-nPos);
        CMEDOS = (sumCosphi[4]*sumCosphi[5]+sumSinphi[4]*sumSinphi[5])/(nPos*nNeg);
        fDnnPro[1][2]->Fill(centrality,CMEDOS,nPos*nNeg);
    
        //Gamma112:
        CMEGSS = (sumCosphi[4]*sumCosphi[4]-sumSinphi[4]*sumSinphi[4]-sumCos2phi[4])/(nPos*nPos-nPos);
        fCMEPro[0][2]->Fill(centrality,CMEGSS,nPos*nPos-nPos); //Same-Sign
        CMEGOS = (sumCosphi[4]*sumCosphi[5]-sumSinphi[4]*sumSinphi[5])/(nPos*nNeg);
        fCMEPro[1][2]->Fill(centrality,CMEGOS,nPos*nNeg);     //Oppo-Sign
        //=================== CME filling done ==========================
    }

    if (calcTypBF){
        //===================Calc balance function loop trk again!=======================
        mPos = 0, mNeg = 0;
        for (int iTrk = 0; iTrk < NTRKMAX; ++iTrk){
    
          double pt1   = pt[iTrk];
          double eta1 = eta[iTrk];
          double phi1 = phi[iTrk];
          double ch1 = charge[iTrk];
    
          if (ch1==0) break;
          if (pt1>1.5 || pt1<0.3 ||fabs(eta1)>0.8) continue;
    
          for (int jTrk=0; jTrk < NTRKMAX; ++jTrk){
    
            if (jTrk == iTrk) continue;
    
            double pt2 = pt[jTrk];
            double eta2 = eta[jTrk];
            double phi2 = phi[jTrk];
            double ch2 = charge[jTrk];
    
            if (ch2==0) continue;
            if (pt2>1.5 || pt2<0.3 ||fabs(eta2)>0.8) continue;
    
            double DEta = fabs(eta1 - eta2);
            double DPt = pt1 - pt2;
            double DPhi = phi1 - phi2;
    
            // Balance function
            double cpp[4] = {-999.,-999.,-999.,-999.};
            double fpp[4] = {-999.,-999.,-999.,-999.};
            if (ch1  > 0){
              if (ch2 > 0) cpp[0] = DEta; // PP
              if (ch2 < 0) cpp[1] = DEta; //PN
    
              double rndmCharge3 = gRandom->Rndm();
              if (rndmCharge3>0.5) fpp[0]=DEta; // PP
              else fpp[1]=DEta; // PN
            }
            else if (ch1  < 0){
              if (ch2 > 0) cpp[2] = DEta; // NP
              if (ch2 < 0) cpp[3] = DEta; // NN
    
              double rndmCharge3 = gRandom->Rndm();
              if (rndmCharge3>0.5) fpp[2]=DEta; // NP
              else fpp[3]=DEta; // NN
            }
    
            cPP[2]->Fill(cpp);
            fPP[2]->Fill(fpp);
    
          } // Loop Inner Track Done
          
          if (ch1 > 0) mPos ++; 
          else mNeg ++;
        }

        hTotEvtMQ->Fill(4.5, mPos);
        hTotEvtMQ->Fill(5.5, mNeg);
        TotEvtMQP[2] += mPos;
        TotEvtMQN[2] += mNeg;
        //===================Calc balance function done!=======================
    }

  } // evt

  //--------------------
  // Ach-dv2 slope
  //--------------------
  TGraphErrors* g_dv2_ach[3]={NULL};
  TGraphErrors* g_dv2_ach_norm[3]={NULL};
  TF1* f1 = new TF1("f1", "pol1", -0.1, 0.1);

  if (calcTypCMW && calcOnline){
    for (int s=0; s<3; ++s){
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
  }

  //--------------------------
  // finalise CME corr
  //--------------------------
  TH1F* DeltaG112[3]={NULL};
  TH1F* DeltaD11[3]={NULL};
  if (calcTypCME && calcOnline){
      for(Int_t s=0; s<3; ++s){

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
            if (s!=1){
              fCMEHist[k][s]->SetBinContent(c,Corr);
              fCMEHist[k][s]->SetBinError(c,CorrErr);
            }
            else {
              double xfCME = fCMEVal[c-1];
              int binfCME = fCMEHist[k][s]->GetXaxis()->FindBin(xfCME);
              fCMEHist[k][s]->SetBinContent(binfCME, Corr);
              fCMEHist[k][s]->SetBinError(binfCME, CorrErr);              
            }
          } // end loop c : Cent
          fCMEPro[k][s]->GetXaxis()->SetRange(1,-1);
        } // end loop k : SO & SS 
      
        //------ OS-SS for \Gamma_NN ----------
        if (s==0) DeltaG112[s] = (TH1F*)(fCMEHist[1][s]->Clone("DeltaG112_AftLCC"));
        if (s==1) DeltaG112[s] = (TH1F*)(fCMEHist[1][s]->Clone("DeltaG112_AftCME"));
        if (s==2) DeltaG112[s] = (TH1F*)(fCMEHist[1][s]->Clone("DeltaG112_AftCMW"));
        DeltaG112[s]->Add(fCMEHist[0][s],-1);

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
            if (s!=1){
              fDnnHist[k][s]->SetBinContent(c,Corr);
              fDnnHist[k][s]->SetBinError(c,CorrErr);
            } 
            else {
              double xfCME = fCMEVal[c-1];
              int binfCME = fDnnHist[k][s]->GetXaxis()->FindBin(xfCME);
              fDnnHist[k][s]->SetBinContent(binfCME, Corr);
              fDnnHist[k][s]->SetBinError(binfCME, CorrErr);              
            }
          } // end loop c
          fDnnPro[k][s]->GetXaxis()->SetRange(1,-1);
        } // end loop k
  
        //------ OS-SS for \Delta_NN ----------
        if (s==0) DeltaD11[s] = (TH1F*)(fDnnHist[1][s]->Clone("DeltaD11_AftLCC"));
        if (s==1) DeltaD11[s] = (TH1F*)(fDnnHist[1][s]->Clone("DeltaD11_AftCME"));
        if (s==2) DeltaD11[s] = (TH1F*)(fDnnHist[1][s]->Clone("DeltaD11_AftCMW"));
        DeltaD11[s]->Add(fDnnHist[0][s],-1);

      } // end loop s
  }

  //--------------------------
  // Balance func.
  //--------------------------
  if (calcTypBF && calcOnline){
      for (int s=0; s<3; ++s){
        TH1D* fPP_temp0 = (TH1D*)fPP[s]->Projection(0);
        TH1D* fPN_temp0 = (TH1D*)fPP[s]->Projection(1);
        TH1D* fNP_temp0 = (TH1D*)fPP[s]->Projection(2);
        TH1D* fNN_temp0 = (TH1D*)fPP[s]->Projection(3);

        fPN_temp0->Scale(1./fPN_temp0->Integral(1,-1));
        fNP_temp0->Scale(1./fNP_temp0->Integral(1,-1));
        fPP_temp0->Scale(1./fPP_temp0->Integral(1,-1));
        fNN_temp0->Scale(1./fNN_temp0->Integral(1,-1));
      
        TH1D* cPP_temp0 = (TH1D*)cPP[s]->Projection(0);
        TH1D* cPN_temp0 = (TH1D*)cPP[s]->Projection(1);
        TH1D* cNP_temp0 = (TH1D*)cPP[s]->Projection(2);
        TH1D* cNN_temp0 = (TH1D*)cPP[s]->Projection(3);

        cPP_temp0->Scale(1./TotEvtMQP[s]);
        cPN_temp0->Scale(1./TotEvtMQN[s]);
        cNP_temp0->Scale(1./TotEvtMQP[s]);
        cNN_temp0->Scale(1./TotEvtMQN[s]);

        cPN_temp0->Divide(fPN_temp0);
        cNP_temp0->Divide(fNP_temp0);
        cPP_temp0->Divide(fPP_temp0);
        cNN_temp0->Divide(fNN_temp0);
      
        for(Int_t iEta=1;iEta<=cPN_temp0->GetNbinsX();iEta++) { 
          double cpn_ieta = cPN_temp0->GetBinContent(iEta);
          double cpp_ieta = cPP_temp0->GetBinContent(iEta);
          double cnp_ieta = cNP_temp0->GetBinContent(iEta);
          double cnn_ieta = cNN_temp0->GetBinContent(iEta);
          double bf_ieta = (cpn_ieta + cnp_ieta - cpp_ieta - cnn_ieta)/2.;
          // double bf_ieta =( (cpn_ieta - cpp_ieta)/TotEvtMQP + (cnp_ieta - cnn_ieta)/TotEvtMQN ) / 2.;
      
          double cpnErr_ieta = cPN_temp0->GetBinError(iEta);
          double cppErr_ieta = cPP_temp0->GetBinError(iEta);
          double cnpErr_ieta = cNP_temp0->GetBinError(iEta);
          double cnnErr_ieta = cNN_temp0->GetBinError(iEta);
          double bfErr_ieta = sqrt(pow(cpnErr_ieta,2) + pow(cppErr_ieta,2)
                                            + pow(cnpErr_ieta,2) + pow(cnnErr_ieta,2))/2;
          // double bfErr_ieta = sqrt(pow(cpnErr_ieta,2)/pow(TotEvtMQP,2) + pow(cppErr_ieta,2)/pow(TotEvtMQP,2)
          //                                   + pow(cnpErr_ieta,2)/pow(TotEvtMQN,2) + pow(cnnErr_ieta,2)/pow(TotEvtMQN,2))/2;
          BF0[s]->SetBinContent(iEta, bf_ieta);
          BF0[s]->SetBinError(iEta, bfErr_ieta);
        }
        double yVal[31];
        for (int i=0; i<31; ++i){yVal[i]=BF0[s]->GetBinContent(i+1);}
        
        double dE[31];
        for (int i=0; i<31; ++i){dE[i]=BF0[s]->GetBinCenter(i+1);}
        double em;
        double sumbf;
        for (int i=0; i<31; ++i){em += fabs(dE[i])*yVal[i]; sumbf+=yVal[i];}
        cout<<"step "<<s+1<<     "mean val of |Deta| : "<<em/sumbf<<"   "<<BF0[s]->GetMean(1)<<"   "<<BF0[s]->GetMeanError(1)<<endl;
        fPP_temp0->Delete();
        fPN_temp0->Delete();
        fNP_temp0->Delete();
        fNN_temp0->Delete();

        cPP_temp0->Delete();
        cPN_temp0->Delete();
        cNP_temp0->Delete();
        cNN_temp0->Delete();
      }
  } // End loop s

  //----------------
  // done
  //----------------
  TFile* f = new TFile("../output/out_comp_Cent0_BF.root", "RECREATE");
  f->cd();
    hXYEmisPoint->Scale(1./nEvt);
    hXYEmisPoint->Write();
    for (int i = 0; i < 9; ++i) {
      pChargeXYEmisPoint[i]->Write();
    }
    hAch->Write();
    hPt ->Write();
    hEta->Write();
    hPhi->Write();
    hV2IntPro->Write();
    for (int i=0; i<1; ++i){
      fMultEvt[i]->Write();
      v2ptPos[i]->Write();
      v2ptNeg[i]->Write();
      hPtspecPos[i]->Write();
      hPtspecNeg[i]->Write();
    }
  
    if (calcTypCMW){  
      pAch_nSwitchCMW->Write();
      for (int i=0; i<1; ++i){
        pv2posAch[i]->Write();
        pv2negAch[i]->Write();
        if (calcOnline) g_dv2_ach[i]->Write();
        if (calcOnline) g_dv2_ach_norm[i]->Write();
      }
    }
  
    if (calcTypCME){   
      for (int i=0; i<1; ++i){
        fDnnPro[0][i]->Write();
        fDnnPro[1][i]->Write();
        fCMEPro[0][i]->Write();
        fCMEPro[1][i]->Write();
        if (calcOnline) DeltaG112[i]->Write();
        if (calcOnline) DeltaD11[i]->Write();   
      }
    }
    if (calcTypBF){
      for (int i=0; i<1; ++i){
        BF0[i]->Write();
        fPP[i]->Write();
        cPP[i]->Write(); 
        hTotEvtMQ->Write();
      }
    }
  
    f->Close(); // End write into .root file
}

double RadSys(double rad, double phi, double Rx) {
  return sqrt(pow(rad*cos(phi)/Rx,2.)+pow(rad*sin(phi),2.));
}

double BosPhi(double phi, double Rx) {
  double phis = TMath::ATan(Rx*Rx*TMath::Tan(phi));
  if (phi>TMath::Pi()/2. && phi<TMath::Pi()*3./2.)     phis += TMath::Pi();
  else if (phi<2*TMath::Pi() && phi>TMath::Pi()*3./2.) phis += 2*TMath::Pi();
  return phis;
}

