/*****************************************
 *                                       *
 *           Giorgia Rauco               *
 *           July 2014, 7th              *
 *                                       *
 *              VBF Hbb                  * 
 *                                       *
 *  Validation of bjet regression in     * 
 *  Z(->dielectron, dimuon) + 2bjets     * 
 *                                       *
 *****************************************/

#include <TROOT.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TGaxis.h>
#include "TLegend.h"
#include "TH2.h"
#include "TStyle.h"
#include "TColor.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TF1.h"
#include "TAxis.h"
#include "TString.h"
#include "TLatex"

#include <math.h>
#include <stdio.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <stdlib.h>

void Z2bjets(){

  // pt(bb system)/pt(Z) histograms declaration (Unregressed and Regressed)

  int a;
  const unsigned int N = 11; //11 intervals of pt(Z) 
  TH1F *hptRatio_bbmm[N]; 
  TH1F *hptRatio_bbee[N]; 
  TH1F *hptRatio_bbmcmm[N]; 
  TH1F *hptRatio_bbmcee[N]; 
  TH1F *hptRegRatio_bbmm[N];
  TH1F *hptRegRatio_bbee[N]; 
  TH1F *hptRegRatio_bbmcmm[N]; 
  TH1F *hptRegRatio_bbmcee[N]; 
  TH1F *hptRatio_bbmctmm[N];
  TH1F *hptRatio_bbmctee[N] ;
  TH1F *hptRegRatio_bbmctmm[N];
  TH1F *hptRegRatio_bbmctee[N];
 
  char namebbrmm[100], namebbree[100], namebbrmcmm[100], namebbrmcee[100],  namebbrrmm[100], namebbrree[100], namebbrrmcmm[100], namebbrrmcee[100],  namebbrmctmm[100], namebbrmctee[100],  namebbrrmctmm[100], namebbrrmctee[100];

  for(a=0; a<N; a++){
    //declare unregressed pTratio histograms
    sprintf(namebbrmm, "ptRatio_bbmm%d", a);
    hptRatio_bbmm[a] = new TH1F(namebbrmm,namebbrmm,40,0,2);
    sprintf(namebbree, "ptRatio_bbee%d", a);
    hptRatio_bbee[a]  = new TH1F(namebbree,namebbree,40,0,2);
    sprintf(namebbrmcmm, "ptRatio_bbmcmm%d", a);
    hptRatio_bbmcmm[a] = new TH1F(namebbrmcmm,namebbrmcmm,40,0,2);
    sprintf(namebbrmcee, "ptRatio_bbmcee%d", a);
    hptRatio_bbmcee[a]  = new TH1F(namebbrmcee,namebbrmcee,40,0,2);
    printf(namebbrmctmm, "ptRatio_bbmctmm%d", a);
    hptRatio_bbmctmm[a] = new TH1F(namebbrmctmm,namebbrmctmm,40,0,2);
    printf(namebbrmctee, "ptRatio_bbmctee%d", a);
    hptRatio_bbmctee[a]  = new TH1F(namebbrmctee,namebbrmctee,40,0,2);
    printf(namebbrrmctmm, "ptRegRatio_bbmctmm%d", a);

    //declare regressed pTratio histograms
    sprintf(namebbrrmm, "ptRegRatio_bbmm%d", a);
    hptRegRatio_bbmm[a] = new TH1F(namebbrrmm,namebbrrmm,40,0,2);
    sprintf(namebbrree, "ptRegRatio_bbee%d", a);
    hptRegRatio_bbee[a] = new TH1F(namebbrree,namebbrree,40,0,2);
    sprintf(namebbrrmcmm, "ptRegRatio_bbmcmm%d", a);
    hptRegRatio_bbmcmm[a] = new TH1F(namebbrrmcmm,namebbrrmcmm,40,0,2);
    sprintf(namebbrrmcee, "ptRegRatio_bbmcee%d", a);
    hptRegRatio_bbmcee[a] = new TH1F(namebbrrmcee,namebbrrmcee,40,0,2);
    hptRegRatio_bbmctmm[a] = new TH1F(namebbrrmctmm,namebbrrmctmm,40,0,2);
    printf(namebbrrmctee, "ptRegRatio_bbmctee%d", a);
    hptRegRatio_bbmctee[a] = new TH1F(namebbrrmctee,namebbrrmctee,40,0,2);
    
  }

  //declare pt(bb system)/pt(Z), not in different pt(Z) intervals
  hptbalance_ee = new TH1F("ptbalance_ee", "ptbalance_ee",20,0,2);
  hptbalance_mm = new TH1F("ptbalance_mm", "ptbalance_mm",20,0,2);
  hptbalance_mcee = new TH1F("ptbalance_mcee", "ptbalance_mcee",20,0,2);
  hptbalance_mcmm = new TH1F("ptbalance_mcmm", "ptbalance_mcmm",20,0,2);
  hptbalance_mctee = new TH1F("ptbalance_mctee", "ptbalance_mctee",20,0,2);
  hptbalance_mctmm = new TH1F("ptbalance_mctmm", "ptbalance_mctmm",20,0,2);
  hptbalance_ree = new TH1F("ptbalance_ree", "ptbalance_ree",20,0,2);
  hptbalance_rmm = new TH1F("ptbalance_rmm", "ptbalance_rmm",20,0,2);
  hptbalance_rmcee = new TH1F("ptbalance_rmcee", "ptbalance_rmcee",20,0,2);
  hptbalance_rmcmm = new TH1F("ptbalance_rmcmm", "ptbalance_rmcmm",20,0,2);
  hptbalance_rmctee = new TH1F("ptbalance_rmctee", "ptbalance_rmctee",20,0,2);
  hptbalance_rmctmm = new TH1F("ptbalance_rmctmm", "ptbalance_rmctmm",20,0,2);

  //declare bb system invariant mass, not in different pt(Z) intervals
  hmbb_ee = new TH1F("mbb_ee", "mbb_ee", 30,0,300);
  hmbb_mm = new TH1F("mbb_mm", "mbb_mm", 30,0,300);
  hmbb_mcee = new TH1F("mbb_mcee", "mbb_mcee", 30,0,300);
  hmbb_mcmm = new TH1F("mbb_mcmm", "mbb_mcmm", 30,0,300);
  hmbb_mctee = new TH1F("mbb_mctee", "mbb_mctee", 30,0,300);
  hmbb_mctmm = new TH1F("mbb_mctmm", "mbb_mctmm", 30,0,300);
  hmbb_ree = new TH1F("mbb_ree", "mbb_ree", 30,0,300);
  hmbb_rmm = new TH1F("mbb_rmm", "mbb_rmm", 30,0,300);
  hmbb_rmcee = new TH1F("mbb_rmcee", "mbb_rmcee", 30,0,300);
  hmbb_rmcmm = new TH1F("mbb_rmcmm", "mbb_rmcmm", 30,0,300);
  hmbb_rmctee = new TH1F("mbb_rmctee", "mbb_rmctee", 30,0,300);
  hmbb_rmctmm = new TH1F("mbb_rmctmm", "mbb_rmctmm", 30,0,300);

  //open root files
  TFile *DoubleMu = TFile::Open("DoubleMu.root");
  DoubleMu->cd("Hbb");
  TTree * Tree_DoubleMu= (TTree*)gDirectory->Get("events");
  TFile *DoubleElectron = TFile::Open("DoubleElectron.root");
  DoubleElectron->cd("Hbb");
  TTree * Tree_DoubleElectron= (TTree*)gDirectory->Get("events");
  TFile *DYLLJets = TFile::Open("DYLLJets.root");
  DYLLJets->cd("Hbb");
  TTree * Tree_DYLLJets= (TTree*)gDirectory->Get("events");
  TFile *TTJets = TFile::Open("TTJets.root");
  TTJets->cd("Hbb");
  TTree * Tree_TTJets= (TTree*)gDirectory->Get("events");

  //take useful variables from branches
  vector<bool> *triggerResult_DoubleMu, *triggerResult_DoubleElectron, *triggerResult_DYLLJets, *triggerResult_TTJets;
  Tree_DoubleMu->SetBranchAddress("triggerResult", &triggerResult_DoubleMu);
  Tree_DoubleElectron->SetBranchAddress("triggerResult", &triggerResult_DoubleElectron);
  Tree_DYLLJets->SetBranchAddress("triggerResult", &triggerResult_DYLLJets);
  Tree_TTJets->SetBranchAddress("triggerResult", &triggerResult_TTJets);

  vector<float> *jetPhi_mm,*jetPhi_ee,*jetPhi_mc, *jetPhi_mct;
  Tree_DoubleMu->SetBranchAddress("jetPhi", &jetPhi_mm);
  Tree_DoubleElectron->SetBranchAddress("jetPhi", &jetPhi_ee);
  Tree_DYLLJets->SetBranchAddress("jetPhi", &jetPhi_mc);
  Tree_TTJets->SetBranchAddress("jetPhi", &jetPhi_mct);
  vector<float> *jetEta_mm,*jetEta_ee,*jetEta_mc, *jetEta_mct;
  Tree_DoubleMu->SetBranchAddress("jetEta", &jetEta_mm);
  Tree_DoubleElectron->SetBranchAddress("jetEta", &jetEta_ee);
  Tree_DYLLJets->SetBranchAddress("jetEta", &jetEta_mc);
  Tree_TTJets->SetBranchAddress("jetEta", &jetEta_mct);
  vector<float> *jetPt_mm, *jetPt_ee,*jetPt_mc, *jetPt_mct;
  Tree_DoubleMu->SetBranchAddress("jetPt", &jetPt_mm);
  Tree_DoubleElectron->SetBranchAddress("jetPt", &jetPt_ee);
  Tree_DYLLJets->SetBranchAddress("jetPt", &jetPt_mc);
  Tree_TTJets->SetBranchAddress("jetPt", &jetPt_mct);
  vector<float> *jetE_mm,*jetE_ee,*jetE_mc, *jetE_mct;
  Tree_DoubleMu->SetBranchAddress("jetEnergy", &jetE_mm);
  Tree_DoubleElectron->SetBranchAddress("jetEnergy", &jetE_ee);
  Tree_DYLLJets->SetBranchAddress("jetEnergy", &jetE_mc);
  Tree_TTJets->SetBranchAddress("jetEnergy", &jetE_mct);

  vector<float> *lepPhi_mm,*lepPhi_ee,*lepPhi_mc, *lepPhi_mct;
  Tree_DoubleMu->SetBranchAddress("lepPhi", &lepPhi_mm);
  Tree_DoubleElectron->SetBranchAddress("lepPhi", &lepPhi_ee);
  Tree_DYLLJets->SetBranchAddress("lepPhi", &lepPhi_mc);
  Tree_TTJets->SetBranchAddress("lepPhi", &lepPhi_mct);
  vector<float> *lepEta_mm,*lepEta_ee,*lepEta_mc, *lepEta_mct;
  Tree_DoubleMu->SetBranchAddress("lepEta", &lepEta_mm);
  Tree_DoubleElectron->SetBranchAddress("lepEta", &lepEta_ee);
  Tree_DYLLJets->SetBranchAddress("lepEta", &lepEta_mc);
  Tree_TTJets->SetBranchAddress("lepEta", &lepEta_mct);
  vector<float> *lepPt_mm, *lepPt_ee,*lepPt_mc, *lepPt_mct;
  Tree_DoubleMu->SetBranchAddress("lepPt", &lepPt_mm);
  Tree_DoubleElectron->SetBranchAddress("lepPt", &lepPt_ee);
  Tree_DYLLJets->SetBranchAddress("lepPt", &lepPt_mc);
  Tree_TTJets->SetBranchAddress("lepPt", &lepPt_mct);
  vector<float> *lepE_mm,*lepE_ee,*lepE_mc, *lepE_mct;
  Tree_DoubleMu->SetBranchAddress("lepE", &lepE_mm);
  Tree_DoubleElectron->SetBranchAddress("lepE", &lepE_ee);
  Tree_DYLLJets->SetBranchAddress("lepE", &lepE_mc);
  Tree_TTJets->SetBranchAddress("lepE", &lepE_mct);

  vector<bool> *jetBtagM_ee, *jetBtagM_mm, *jetBtagM_mc, *jetBtagM_mct;
  Tree_DoubleMu->SetBranchAddress("jetBtagM", &jetBtagM_mm);
  Tree_DoubleElectron->SetBranchAddress("jetBtagM", &jetBtagM_ee);
  Tree_DYLLJets->SetBranchAddress("jetBtagM", &jetBtagM_mc);
  Tree_TTJets->SetBranchAddress("jetBtagM", &jetBtagM_mct);

  int nBJets_ee, nBJets_mm, nBJets_mc, nBJets_mct;
  Tree_DoubleMu->SetBranchAddress("nBJets", &nBJets_mm);
  Tree_DoubleElectron->SetBranchAddress("nBJets", &nBJets_ee);
  Tree_DYLLJets->SetBranchAddress("nBJets", &nBJets_mc); 
  Tree_TTJets->SetBranchAddress("nBJets", &nBJets_mct);

  int nJets_ee=0, nJets_mm=0, nJets_mc=0, nJets_mct;
  Tree_DoubleMu->SetBranchAddress("nJets", &nJets_mm);
  Tree_DoubleElectron->SetBranchAddress("nJets", &nJets_ee);
  Tree_DYLLJets->SetBranchAddress("nJets", &nJets_mc);  
  Tree_TTJets->SetBranchAddress("nJets", &nJets_mct);  

  vector<bool> *jetpuidm_mm,*jetpuidm_ee,*jetpuidm_mc, *jetpuidm_mct;
  Tree_DoubleMu->SetBranchAddress("jetPuIdM", &jetpuidm_mm);
  Tree_DoubleElectron->SetBranchAddress("jetPuIdM", &jetpuidm_ee);
  Tree_DYLLJets->SetBranchAddress("jetPuIdM", &jetpuidm_mc);
  Tree_TTJets->SetBranchAddress("jetPuIdM", &jetpuidm_mct);

  vector<float> *jetregpt_mm, *jetregpt_ee, *jetregpt_mc, *jetregpt_mct;
  Tree_DoubleMu->SetBranchAddress("jetRegPt", &jetregpt_mm);
  Tree_DoubleElectron->SetBranchAddress("jetRegPt", &jetregpt_ee);
  Tree_DYLLJets->SetBranchAddress("jetRegPt", &jetregpt_mc);
  Tree_TTJets->SetBranchAddress("jetRegPt", &jetregpt_mct);
  
  float puWt;
  Tree_DYLLJets->SetBranchAddress("puWt", &puWt);
  
  //luminosity coefficient to rescale montecarlo
  float dataLumi=19.6;
  float mcDYLumi=8.539;
  float mcTTLumi=27.772;

  //get total number of entries for each dataset
  int n_DoubleMu = Tree_DoubleMu->GetEntries();
  int n_DoubleElectron = Tree_DoubleElectron->GetEntries();
  int n_DYLLJets = Tree_DYLLJets->GetEntries();
  int n_TTJets = Tree_TTJets->GetEntries();
  
  //trigger result boolean
  bool dyismu=false, dyise=false, dyise=false, dataismu=false, dataise=false, ttise=false, ttismu=false;

  //initialize counters
  int countmll_ee=0, countmll_mm=0;
  int countphi_ee=0, countphi_mm=0;
  int countmll_mcee=0, countphi_mcee=0;
  int countmll_mcmm=0, countphi_mcmm=0;
  int countmll_mctee=0, countphi_mctee=0;
  int countmll_mctmm=0, countphi_mctmm=0;
  int countbtag_mm=0, countbtag_ee=0;
  int countbtag_mcee=0, countbtag_mcmm=0;
  int countbtag_mctee=0, countbtag_mctmm=0;
  int countsoftjet_mm=0, countsoftjet_ee=0;
  int countsoftjet_mcmm=0, countsoftjet_mcee=0;
  int countsoftjet_mctmm=0, countsoftjet_mctee=0;
  int countmcmm_tot=0, countmcee_tot=0;
  int countmctmm_tot=0, countmctee_tot=0;
  int countmm_tot=0, countee_tot=0;
  int count1_mm=0, count1_ee=0, count1_mcmm=0, count1_mcee=0, count1_mctmm=0, count1_mctee=0;
  int countbtag2_mm=0, countbtag2_ee=0, countbtag2_mcmm=0, countbtag2_mcee=0, countbtag2_mctmm=0, countbtag2_mctee=0;
  int countmcmm_tot=0, countmcee_tot=0, countmctmm_tot=0, countmctee_tot=0;

  float ptRatiobins[12] ={30,40,50,60,70,80,90,110,130,150,250,500};

  //DIMUON CYCLE
  for (int ev=0; ev<n_DoubleMu; ++ev) {
    Tree_DoubleMu->GetEntry(ev);
    dataismu = (*triggerResult_DoubleMu)[1];
    if(dataismu){
      
      TLorentzVector lepton0_mm;                      
      TLorentzVector lepton1_mm;
      lepton0_mm.SetPtEtaPhiE((*lepPt_mm)[0], (*lepEta_mm)[0],  (*lepPhi_mm)[0], (*lepE_mm)[0]);
      lepton1_mm.SetPtEtaPhiE((*lepPt_mm)[1], (*lepEta_mm)[1] , (*lepPhi_mm)[1], (*lepE_mm)[1]);
      
      TLorentzVector Z_mm;
      Z_mm = (lepton0_mm+lepton1_mm);
      
      TLorentzVector bjet0_mm;
      bjet0_mm.SetPtEtaPhiE( (*jetPt_mm)[0], (*jetEta_mm)[0],  (*jetPhi_mm)[0], (*jetE_mm)[0]);
      
      bool minv_mm = fabs(Z_mm.M()-91)<20;
      bool btag_mm = (*jetBtagM_mm)[0] && (*jetpuidm_mm)[0];
      TLorentzVector bjet1_mm;
      TLorentzVector bjet2_mm;
      TLorentzVector bb_mm;
      TLorentzVector bb_mmReg;
      bool softjet_bbmm =0, softjet_mm=0;
      float pTbb_mm, deltaphi_bbmm;
      bool phi_bbmm=0, btag_bbmm=0;
      double pTveto_mm=0;


      if(nJets_mm>1){  //requiring at least 2 jets                
	bjet1_mm.SetPtEtaPhiE((*jetPt_mm)[1], (*jetEta_mm)[1],  (*jetPhi_mm)[1], (*jetE_mm)[1]);
	bb_mm = (bjet0_mm + bjet1_mm);
	bb_mmReg =(bjet0_mm*(*jetregpt_mm)[0] + bjet1_mm*(*jetregpt_mm)[1]); 
	pTbb_mm=bb_mm.Pt();
	deltaphi_bbmm=  Z_mm.DeltaPhi(bb_mm);
	phi_bbmm = fabs(deltaphi_bbmm)>2.5;
	btag_bbmm = (*jetBtagM_mm)[1] && (*jetpuidm_mm)[1];
      		
	if(nJets_mm>2){ //if there are >=3 jets, then put a veto on the additional (third) jet
	  bjet2_mm.SetPtEtaPhiE((*jetPt_mm)[2], (*jetEta_mm)[2],  (*jetPhi_mm)[2], (*jetE_mm)[2]);
	  if((0.2*Z_mm.Pt())>10){
	    pTveto_mm = 0.2*Z_mm.Pt();
	  }
	  else{
	    pTveto_mm =10;
	  }
	  softjet_bbmm = bjet2_mm.Pt()>pTveto_mm && ((*jetpuidm_mm)[2]);
	}

	//set counters
	count1_mm++;
	if(minv_mm==true){
	  countmll_mm++;
	  if(btag_mm==true && bjet0_mm.Pt()>10){
	    countbtag_mm++;
	    if(btag_bbmm==true && bjet1_mm.Pt()>10){
	      countbtag2_mm++;
	      if(phi_bbmm==true){
		countphi_mm++;
		if(softjet_bbmm==0){
		  countsoftjet_mm++;
		}
	      }
	    }
	  }
	}
	
	//if all the conditions are satisfied, then fill the histograms
	if(minv_mm==1 && phi_bbmm==1 && btag_mm==1 && btag_bbmm==1 && softjet_bbmm==0 && bjet0_mm.Pt()>10 && bjet1_mm.Pt()>10 ){
	  countmm_tot++;	
	  hptbalance_mm->Fill((bb_mm.Pt())/Z_mm.Pt());
	  hptbalance_rmm->Fill((bb_mmReg.Pt())/Z_mm.Pt());
	  hmbb_mm->Fill(bb_mm.M());
	  hmbb_rmm->Fill(bb_mmReg.M());
	  for(int i=0;i<11;i++){
	    if(Z_mm.Pt()>ptRatiobins[i] && Z_mm.Pt()<ptRatiobins[i+1]){
	      hptRatio_bbmm[i]->Fill(bb_mm.Pt()/Z_mm.Pt());
	      hptRegRatio_bbmm[i]->Fill((bb_mmReg.Pt())/Z_mm.Pt());
	     
	    }
	  }
	}
      }
    }
  } 

  //ELECTRON CYCLE
  for (int ev=0; ev<n_DoubleElectron; ++ev) {
    Tree_DoubleElectron->GetEntry(ev);
    dataise= (*triggerResult_DoubleElectron)[0];
    if(dataise) {
     
      TLorentzVector lepton0_ee;                      
      TLorentzVector lepton1_ee;
      lepton0_ee.SetPtEtaPhiE((*lepPt_ee)[0], (*lepEta_ee)[0],  (*lepPhi_ee)[0], (*lepE_ee)[0]);
      lepton1_ee.SetPtEtaPhiE((*lepPt_ee)[1], (*lepEta_ee)[1] , (*lepPhi_ee)[1], (*lepE_ee)[1]);
      
      TLorentzVector Z_ee;
      Z_ee = (lepton0_ee+lepton1_ee);
     
      TLorentzVector bjet0_ee;
      bjet0_ee.SetPtEtaPhiE( (*jetPt_ee)[0], (*jetEta_ee)[0],  (*jetPhi_ee)[0], (*jetE_ee)[0]);

      bool minv_ee = fabs(Z_ee.M()-91)<20;

      bool btag_ee = (*jetBtagM_ee)[0] && (*jetpuidm_ee)[0]; 
      TLorentzVector bjet1_ee;
      TLorentzVector bb_ee;
      TLorentzVector bb_eeReg;
      TLorentzVector bjet2_ee;
      bool softjet_ee = 0, softjet_bbee=0;
      float pTbb_ee, deltaphi_bbee;
      bool phi_bbee=0, btag_bbee=0;
      double pTveto_ee=0;
      
      if(nJets_ee>1){    //requiring at least 2 jets
	bjet1_ee.SetPtEtaPhiE((*jetPt_ee)[1], (*jetEta_ee)[1],  (*jetPhi_ee)[1], (*jetE_ee)[1]);
	bb_ee = (bjet0_ee + bjet1_ee);	
	bb_eeReg =(bjet0_ee*(*jetregpt_ee)[0] + bjet1_ee*(*jetregpt_ee)[1]); 
	pTbb_ee=bb_ee.Pt();
	deltaphi_bbee=  Z_ee.DeltaPhi(bb_ee);
	phi_bbee = fabs(deltaphi_bbee)>2.5;
	btag_bbee = (*jetBtagM_ee)[1] && (*jetpuidm_ee)[1];

	if(nJets_ee>2){  //if there are >=3 jets, then put a veto on the additional (third) jet
	  bjet2_ee.SetPtEtaPhiE((*jetPt_ee)[2], (*jetEta_ee)[2],  (*jetPhi_ee)[2], (*jetE_ee)[2]);
	  if((0.2*Z_ee.Pt())>10){
	    pTveto_ee = 0.2*Z_ee.Pt();
	  }
	  else{
	    pTveto_ee = 10;
	  }
	  softjet_bbee = bjet2_ee.Pt()>pTveto_ee && ((*jetpuidm_ee)[2]);
	}
	
	//set counters
	count1_ee++;
	if(minv_ee==true){
	  countmll_ee++;
	  if(btag_ee==true && bjet0_ee.Pt()>10){
	    countbtag_ee++;
	    if(btag_bbee==true && bjet1_ee.Pt()>10){
	      countbtag2_ee++;
	      if(phi_bbee==true){
		countphi_ee++;
		if(softjet_bbee==0){
		  countsoftjet_ee++;
		}
	      }
	    }
	  }
	}
	//if all the conditions are satisfied, then fill the histograms
	if(minv_ee==1 && phi_bbee==1 && btag_ee==1 && btag_bbee==1 && softjet_bbee==0  && bjet1_ee.Pt()>10 && bjet0_ee.Pt()>10){
	  countee_tot++;
	
	  hptbalance_ee->Fill((bb_ee.Pt())/Z_ee.Pt());
	  hptbalance_ree->Fill((bb_eeReg.Pt())/Z_ee.Pt());
	  hmbb_ee->Fill(bb_ee.M());
	  hmbb_ree->Fill(bb_eeReg.M());
	  for(int i=0;i<11;i++){
	    if(Z_ee.Pt()>ptRatiobins[i] && Z_ee.Pt()<ptRatiobins[i+1]){
	      hptRatio_bbee[i]->Fill(bb_ee.Pt()/Z_ee.Pt());
	      hptRegRatio_bbee[i]->Fill((bb_eeReg.Pt())/Z_ee.Pt());

	    }
	  }
	}
      }
    }
  }
  
  //MONTE-CARLO DRELL YAN CYCLE
  for (int ev=0; ev<n_DYLLJets; ++ev) {
    Tree_DYLLJets->GetEntry(ev);
    
    dyise = (*triggerResult_DYLLJets)[0] ;
    dyismu = (*triggerResult_DYLLJets)[1] ;
    
    TLorentzVector lepton0_mc;                      
    TLorentzVector lepton1_mc;
    lepton0_mc.SetPtEtaPhiE((*lepPt_mc)[0], (*lepEta_mc)[0],  (*lepPhi_mc)[0],  (*lepE_mc)[0] );
    lepton1_mc.SetPtEtaPhiE((*lepPt_mc)[1], (*lepEta_mc)[1] , (*lepPhi_mc)[1], (*lepE_mc)[1]);
    
    TLorentzVector Z_mc;
    Z_mc = (lepton0_mc+lepton1_mc);
    
    TLorentzVector bjet0_mc;
    bjet0_mc.SetPtEtaPhiE( (*jetPt_mc)[0], (*jetEta_mc)[0],  (*jetPhi_mc)[0], (*jetE_mc)[0]); 
        
    // correct MC for pu weigth
    float coeff = dataLumi/mcDYLumi*puWt[1];
    
    bool minv_mc = fabs(Z_mc.M()-91)<10;
    bool btag_mc = (*jetBtagM_mc)[0] && (*jetpuidm_mc)[0];
    TLorentzVector bjet1_mc;
    TLorentzVector bjet2_mc;
    TLorentzVector bb_mc;
    TLorentzVector bb_mcReg;
    bool softjet_mc = 0;
    bool softjet_bbmc=0;
    float pTbb_mc, deltaphi_bbmc;
    bool phi_bbmc=0, btag_bbmc=0;
    double pTveto_mc=0;
    
    if(nJets_mc>1){  //requiring at least 3 jets
      bjet1_mc.SetPtEtaPhiE((*jetPt_mc)[1], (*jetEta_mc)[1],  (*jetPhi_mc)[1], (*jetE_mc)[1]);
      bb_mc = (bjet0_mc + bjet1_mc);
      bb_mcReg =(bjet0_mc*(*jetregpt_mc)[0] + bjet1_mc*(*jetregpt_mc)[1]); 
      pTbb_mc=bb_mc.Pt();
      deltaphi_bbmc= Z_mc.DeltaPhi(bb_mc);
      phi_bbmc=fabs(deltaphi_bbmc)>2.5;
      btag_bbmc = (*jetBtagM_mc)[1] && (*jetpuidm_mc)[1];
      if(nJets_mc>2){ //if there are >=3 jets, then put a veto on the additional (third) jet
	bjet2_mc.SetPtEtaPhiE((*jetPt_mc)[2], (*jetEta_mc)[2],  (*jetPhi_mc)[2], (*jetE_mc)[2]);
	if((0.2*Z_mc.Pt())>10){
	  pTveto_mc = 0.2*Z_mc.Pt();
	}
	else{
	  pTveto_mc = 10;
	}
	softjet_bbmc = bjet2_mc.Pt()>pTveto_mc && ((*jetpuidm_mc)[2]);
      }


      if(dyise){  //dielectron trigger and set counters
	count1_mcee++;
	if(minv_mc==true){
	  countmll_mcee++;
	  if(btag_mc==true && bjet0_mc.Pt()>10){
	    countbtag_mcee++;
	    if(btag_bbmc==true && bjet1_mc.Pt()>10){
	      countbtag2_mcee++;
	      if(phi_bbmc==true){
		countphi_mcee++;
		if(softjet_bbmc==0){
		  countsoftjet_mcee++;
		}
	      }
	    }
	  }
	}
      }
      
      if(dyismu){ //dimuon trigger and set counters
	count1_mcmm++;
	if(minv_mc==true){
	  countmll_mcmm++;
	  if(btag_mc==true && bjet0_mc.Pt()>10){
	    countbtag_mcmm++;
	    if(btag_bbmc==true && bjet1_mc.Pt()>10){
	      countbtag2_mcmm++;
	      if(phi_bbmc==true){
		countphi_mcmm++;
		if(softjet_bbmc==0){
		  countsoftjet_mcmm++;
		}
	      }
	    }
	  }
	}
      }

      //if all the condition are satisfied, fill histograms
      if(minv_mc==1 && phi_bbmc==1 && btag_bbmc==1 && btag_mc==1 && softjet_bbmc==0 && bjet0_mc.Pt()>10 && bjet1_mc.Pt()>10){
	if(dyismu){
	  countmcmm_tot++;
	  hptbalance_mcmm->Fill((bb_mc.Pt())/Z_mc.Pt(),coeff);
	  hptbalance_rmcmm->Fill((bb_mcReg.Pt())/Z_mc.Pt(),coeff);
	  hmbb_mcmm->Fill(bb_mc.M(),coeff);
	  hmbb_rmcmm->Fill(bb_mcReg.M(),coeff);
	  for(int i=0;i<11;i++){
	    if(Z_mc.Pt()>ptRatiobins[i] && Z_mc.Pt()<ptRatiobins[i+1]){
	      hptRatio_bbmcmm[i]->Fill(bb_mc.Pt()/Z_mc.Pt(),coeff);
	      hptRegRatio_bbmcmm[i]->Fill((bb_mcReg.Pt())/Z_mc.Pt(),coeff);	   
	    }
	  }
	}
	if(dyise){	  
	  countmcee_tot++;
	  hptbalance_mcee->Fill((bb_mc.Pt())/Z_mc.Pt(),coeff);
	  hptbalance_rmcee->Fill((bb_mcReg.Pt())/Z_mc.Pt(),coeff);
	  hmbb_mcee->Fill(bb_mc.M(),coeff);
	  hmbb_rmcee->Fill(bb_mcReg.M(),coeff);
	  for(int i=0;i<11;i++){
	    if(Z_mc.Pt()>ptRatiobins[i] && Z_mc.Pt()<ptRatiobins[i+1]){
	      hptRatio_bbmcee[i]->Fill(bb_mc.Pt()/Z_mc.Pt(),coeff);
	      hptRegRatio_bbmcee[i]->Fill((bb_mcReg.Pt())/Z_mc.Pt(),coeff);	     
	    }
	  }
	}	
      }
    }
  }
  

  //MONTE-CARLO TTbar CYCLE
  for (int ev=0; ev<n_TTJets; ++ev) {
    Tree_TTJets->GetEntry(ev);
    
    ttise = (*triggerResult_TTJets)[0] ;
    ttismu = (*triggerResult_TTJets)[1] ;
    
    TLorentzVector lepton0_mct;                      
    TLorentzVector lepton1_mct;
    lepton0_mct.SetPtEtaPhiE((*lepPt_mct)[0], (*lepEta_mct)[0],  (*lepPhi_mct)[0],  (*lepE_mct)[0] );
    lepton1_mct.SetPtEtaPhiE((*lepPt_mct)[1], (*lepEta_mct)[1] , (*lepPhi_mct)[1], (*lepE_mct)[1]);
    
    TLorentzVector Z_mct;
    Z_mct = (lepton0_mct+lepton1_mct);
    
    TLorentzVector bjet0_mct;
    bjet0_mct.SetPtEtaPhiE( (*jetPt_mct)[0], (*jetEta_mct)[0],  (*jetPhi_mct)[0], (*jetE_mct)[0]); 
        
    // correct MC for pu weigth
    float coefft = dataLumi/mcTTLumi*puWt[1];
    
    bool minv_mct = fabs(Z_mct.M()-91)<20;
    bool btag_mct = (*jetBtagM_mct)[0] && (*jetpuidm_mct)[0];
    TLorentzVector bjet1_mct;
    TLorentzVector bjet2_mct;
    TLorentzVector bb_mct;
    TLorentzVector bb_mctReg;
    bool softjet_mct = 0;
    bool softjet_bbmct=0;
    float pTbb_mct, deltaphi_bbmct;
    bool phi_bbmct=0, btag_bbmct=0;
    double pTveto_mct=0;


    if(nJets_mct>1){  //requiring at least 2 jets
      bjet1_mct.SetPtEtaPhiE((*jetPt_mct)[1], (*jetEta_mct)[1],  (*jetPhi_mct)[1], (*jetE_mct)[1]);
      bb_mct = (bjet0_mct + bjet1_mct);
      bb_mctReg =(bjet0_mct*(*jetregpt_mct)[0] + bjet1_mct*(*jetregpt_mct)[1]); 
      pTbb_mct=bb_mct.Pt();
      deltaphi_bbmct= Z_mct.DeltaPhi(bb_mct);
      phi_bbmct=fabs(deltaphi_bbmct)>2.5;
      btag_bbmct = (*jetBtagM_mct)[1] && (*jetpuidm_mct)[1];
      if(nJets_mct>2){ //if there are >=3 jets, then put a veto on the additional (third) jet
	bjet2_mct.SetPtEtaPhiE((*jetPt_mct)[2], (*jetEta_mct)[2],  (*jetPhi_mct)[2], (*jetE_mct)[2]);
	if((0.2*Z_mct.Pt())>10){
	  pTveto_mct = 0.2*Z_mct.Pt();
	}
	else{
	  pTveto_mct = 10;
	}
	softjet_mct =  bjet2_mct.Pt()>pTveto_mct && ((*jetpuidm_mc)[2]);
      }
      if(ttise){  //dielectron trigger and set counters
	count1_mctee++;
	if(minv_mct==true){
	  countmll_mctee++;
	  if(btag_mct==true && bjet0_mct.Pt()>10){
	    countbtag_mctee++;
	    if(btag_bbmct==true && bjet1_mct.Pt()>10){
	      countbtag2_mctee++;
	      if(phi_bbmct==true){
		countphi_mctee++;
		if(softjet_mct==0){
		  countsoftjet_mctee++;
		}
	      }
	    }
	  }
	}
      }
      
      if(ttismu){  //dimuon trigger and set counters
	count1_mctmm++;
	if(minv_mct==true){
	  countmll_mctmm++;
	  if(btag_mct==true && bjet0_mct.Pt()>10){
	    countbtag_mctmm++;
	    if(btag_bbmct==true && bjet1_mct.Pt()>10){
	      countbtag2_mctmm++;
	      if(phi_bbmct==true){
		countphi_mctmm++;
		if(softjet_mct==0){
		    countsoftjet_mctmm++;
		}
	      }
	    }
	  }
	}
      }
      //if all the conditions are satisfied, fill histograms
      if(minv_mct==1 && phi_bbmct==1 && btag_bbmct==1 && btag_mct==1 && softjet_mct==0 && bjet0_mct.Pt()>10 && bjet1_mct.Pt()>10){
	if(ttismu){
	  countmctmm_tot++;
	  hptbalance_mctmm->Fill((bb_mct.Pt())/Z_mct.Pt(),coefft);
	  hptbalance_rmctmm->Fill((bb_mctReg.Pt())/Z_mct.Pt(),coefft);
	  hmbb_mctmm->Fill(bb_mct.M(),coefft);
	  hmbb_rmctmm->Fill(bb_mctReg.M(),coefft);
	  for(int i=0;i<11;i++){
	    if(Z_mct.Pt()>ptRatiobins[i] && Z_mct.Pt()<ptRatiobins[i+1]){
	      hptRatio_bbmctmm[i]->Fill(bb_mct.Pt()/Z_mct.Pt(),coefft);
	      hptRegRatio_bbmctmm[i]->Fill((bb_mctReg.Pt())/Z_mct.Pt(),coefft);
	    
	    }
	  }
	}
	if(ttise){	  
	  countmctee_tot++;
	  hptbalance_mctee->Fill((bb_mct.Pt())/Z_mct.Pt(),coefft);
	  hptbalance_rmctee->Fill((bb_mctReg.Pt())/Z_mct.Pt(),coefft);
	  hmbb_mctee->Fill(bb_mct.M(),coefft);
	  hmbb_rmctee->Fill(bb_mctReg.M(),coefft);
	  for(int i=0;i<11;i++){
	    if(Z_mct.Pt()>ptRatiobins[i] && Z_mct.Pt()<ptRatiobins[i+1]){
	      hptRatio_bbmctee[i]->Fill(bb_mct.Pt()/Z_mct.Pt(),coefft);
	      hptRegRatio_bbmctee[i]->Fill((bb_mctReg.Pt())/Z_mct.Pt(),coefft);
	    
	    }
	  }
	}	
      }
    }
  }

  //cout counters for dimuon and dielectron

  cout << "*************** MUONS ********************" << endl;
  cout << "Trig+1bjet >>>> \t" << count1_mm << "\t" << count1_mcmm << "\t" << count1_mctmm << endl; 
  cout << "Minv >>>> \t" << countmll_mm << "\t" << countmll_mcmm << "\t" << countmll_mctmm << endl; 
  cout << "Btag + pt>10 >>>> \t" << countbtag_mm << "\t" << countbtag_mcmm << "\t" << countbtag_mctmm << endl;
  cout << "Btag2 + pt>10 >>>> \t" << countbtag2_mm << "\t" << countbtag2_mcmm << "\t" << countbtag2_mctmm << endl;
  cout << "DeltaPhi >>>> \t" << countphi_mm << "\t" << countphi_mcmm << "\t" << countphi_mctmm << endl;
  cout << "Additional jet >>>> \t" << countsoftjet_mm << "\t" << countsoftjet_mcmm << "\t" << countsoftjet_mctmm << endl;
  cout << "All selections >>>> \t" << countmm_tot << "\t" << countmcmm_tot << "\t" << countmctmm_tot << endl;

  cout << "*************** ELECTRONS ********************" << endl;
  cout << "Trig+1bjet >>>> \t" << count1_ee << "\t" << count1_mcee << "\t" << count1_mctee << endl; 
  cout << "Minv >>>> \t" << countmll_ee << "\t" << countmll_mcee << "\t" << countmll_mctee << endl; 
  cout << "Btag + pt>10 >>>> \t" << countbtag_ee << "\t" << countbtag_mcee << "\t" << countbtag_mctee << endl;
  cout << "Btag2 + pt>10 >>>> \t" << countbtag2_ee << "\t" << countbtag2_mcee << "\t" << countbtag2_mctee << endl;
  cout << "DeltaPhi >>>> \t" << countphi_ee << "\t" << countphi_mcee << "\t" << countphi_mctee << endl;
  cout << "Additional jet >>>> \t" << countsoftjet_ee << "\t" << countsoftjet_mcee << "\t" << countsoftjet_mctee << endl;
  cout << "All selections >>>> \t" << countee_tot << "\t" << countmcee_tot << "\t" << countmctee_tot << endl;

  //scale MC DY histograms with respect to Data-TTbar
  hptbalance_mcee->Scale((hptbalance_ee->Integral()-hptbalance_mctee->Integral())/hptbalance_mcee->Integral());
  hptbalance_mcmm->Scale((hptbalance_mm->Integral()-hptbalance_mctmm->Integral())/hptbalance_mcmm->Integral());

  hptbalance_rmcee->Scale((hptbalance_ree->Integral()-hptbalance_rmctee->Integral())/hptbalance_rmcee->Integral());
  hptbalance_rmcmm->Scale((hptbalance_rmm->Integral()-hptbalance_rmctmm->Integral())/hptbalance_rmcmm->Integral());

  hmbb_mcee->Scale((hmbb_ee->Integral()-hmbb_mctee->Integral())/hmbb_mcee->Integral());
  hmbb_mcmm->Scale((hmbb_mm->Integral()-hmbb_mctmm->Integral())/hmbb_mcmm->Integral());

  hmbb_rmcee->Scale((hmbb_ree->Integral()-hmbb_rmctee->Integral())/hmbb_rmcee->Integral());
  hmbb_rmcmm->Scale((hmbb_rmm->Integral()-hmbb_rmctmm->Integral())/hmbb_rmcmm->Integral());

  for(int c=0; c<11; c++){
    hptRatio_bbmcmm[c]->Scale((hptRatio_bbmm[c]->Integral()-hptRatio_bbmctmm[c]->Integral())/hptRatio_bbmcmm[c]->Integral());
    hptRatio_bbmcee[c]->Scale((hptRatio_bbee[c]->Integral()-hptRatio_bbmctee[c]->Integral())/hptRatio_bbmcee[c]->Integral());
    hptRegRatio_bbmcmm[c]->Scale((hptRegRatio_bbmm[c]->Integral()-hptRegRatio_bbmctmm[c]->Integral())/hptRegRatio_bbmcmm[c]->Integral());
    hptRegRatio_bbmcee[c]->Scale((hptRegRatio_bbee[c]->Integral()-hptRegRatio_bbmctee[c]->Integral())/hptRegRatio_bbmcee[c]->Integral());
  }

  //get RMS and MEAN of pTratio response in each pt(Z) interval

  float RMS_bbree[11], RMS_bbrmm[11], RMS_bbrmcee[11], RMS_bbrmcmm[11], RMS_bbrree[11], RMS_bbrrmm[11], RMS_bbrrmcee[11], RMS_bbrrmcmm[11], RMS_bbrmctee[11], RMS_bbrmctmm[11],  RMS_bbrrmctee[11], RMS_bbrrmctmm[11];

  float RMSE_bbree[11], RMSE_bbrmm[11], RMSE_bbrmcee[11], RMSE_bbrmcmm[11], RMSE_bbrree[11], RMSE_bbrrmm[11], RMSE_bbrrmcee[11], RMSE_bbrrmcmm[11],  RMSE_bbrmctee[11], RMSE_bbrmctmm[11], RMSE_bbrrmctee[11], RMSE_bbrrmctmm[11];

  float MEAN_bbree[11], MEAN_bbrmm[11], MEAN_bbrmcee[11], MEAN_bbrmcmm[11], MEAN_bbrree[11], MEAN_bbrrmm[11], MEAN_bbrrmcee[11], MEAN_bbrrmcmm[11],  MEAN_bbrmctee[11], MEAN_bbrmctmm[11], MEAN_bbrrmctee[11], MEAN_bbrrmctmm[11];

  float MEANE_bbree[11], MEANE_bbrmm[11], MEANE_bbrmcee[11], MEANE_bbrmcmm[11], MEANE_bbrree[11], MEANE_bbrrmm[11], MEANE_bbrrmcee[11], MEANE_bbrrmcmm[11],  MEANE_bbrmctee[11], MEANE_bbrmctmm[11], MEANE_bbrrmctee[11], MEANE_bbrrmctmm[11];
  
  float RMS_MCee[11] ,  RMS_MCree[11], RMS_MCmm[11], RMS_MCrmm[11];
  float RMSE_MCee[11] ,   RMSE_MCree[11], RMSE_MCmm[11], RMSE_MCrmm[11];

  float MEAN_MCee[11] ,  MEAN_MCree[11], MEAN_MCmm[11], MEAN_MCrmm[11];
  float MEANE_MCee[11] ,   MEANE_MCree[11],MEANE_MCmm[11], MEANE_MCrmm[11];

  for(int d=0; d<11; d++){
 
    //RMS FOR PTRATIO

    RMS_bbrmm[d]=hptRatio_bbmm[d]->GetRMS(1);
    RMS_bbree[d]=hptRatio_bbee[d]->GetRMS(1);
    RMS_bbrmcmm[d]=hptRatio_bbmcmm[d]->GetRMS(1);
    RMS_bbrmcee[d]=hptRatio_bbmcee[d]->GetRMS(1);
    RMS_bbrmctmm[d]=hptRatio_bbmctmm[d]->GetRMS(1);
    RMS_bbrmctee[d]=hptRatio_bbmctee[d]->GetRMS(1); 
    RMS_bbrrmm[d]=hptRegRatio_bbmm[d]->GetRMS(1);
    RMS_bbrree[d]=hptRegRatio_bbee[d]->GetRMS(1);
    RMS_bbrrmcmm[d]=hptRegRatio_bbmcmm[d]->GetRMS(1);
    RMS_bbrrmcee[d]=hptRegRatio_bbmcee[d]->GetRMS(1);
    RMS_bbrrmctmm[d]=hptRegRatio_bbmctmm[d]->GetRMS(1);
    RMS_bbrrmctee[d]=hptRegRatio_bbmctee[d]->GetRMS(1);
    RMSE_bbrmm[d]=hptRatio_bbmm[d]->GetRMSError(1);
    RMSE_bbree[d]=hptRatio_bbee[d]->GetRMSError(1);
    RMSE_bbrmcmm[d]=hptRatio_bbmcmm[d]->GetRMSError(1);
    RMSE_bbrmcee[d]=hptRatio_bbmcee[d]->GetRMSError(1);
    RMSE_bbrmctmm[d]=hptRatio_bbmctmm[d]->GetRMSError(1);
    RMSE_bbrmctee[d]=hptRatio_bbmctee[d]->GetRMSError(1);
    RMSE_bbrrmm[d]=hptRegRatio_bbmm[d]->GetRMSError(1);
    RMSE_bbrree[d]=hptRegRatio_bbee[d]->GetRMSError(1);
    RMSE_bbrrmcmm[d]=hptRegRatio_bbmcmm[d]->GetRMSError(1);
    RMSE_bbrrmcee[d]=hptRegRatio_bbmcee[d]->GetRMSError(1);
    RMSE_bbrrmctmm[d]=hptRegRatio_bbmctmm[d]->GetRMSError(1);
    RMSE_bbrrmctee[d]=hptRegRatio_bbmctee[d]->GetRMSError(1);  

    //MEAN FOR PTRATIO
    MEAN_bbrmm[d]=hptRatio_bbmm[d]->GetMean(1); 
    MEAN_bbree[d]=hptRatio_bbee[d]->GetMean(1);
    MEAN_bbrmcmm[d]=hptRatio_bbmcmm[d]->GetMean(1);
    MEAN_bbrmcee[d]=hptRatio_bbmcee[d]->GetMean(1);
    MEAN_bbrmctmm[d]=hptRatio_bbmctmm[d]->GetMean(1);
    MEAN_bbrmctee[d]=hptRatio_bbmctee[d]->GetMean(1);
    MEAN_bbrrmm[d]=hptRegRatio_bbmm[d]->GetMean(1); 
    MEAN_bbrree[d]=hptRegRatio_bbee[d]->GetMean(1);
    MEAN_bbrrmcmm[d]=hptRegRatio_bbmcmm[d]->GetMean(1);
    MEAN_bbrrmcee[d]=hptRegRatio_bbmcee[d]->GetMean(1);
    MEAN_bbrrmctmm[d]=hptRegRatio_bbmctmm[d]->GetMean(1);
    MEAN_bbrrmctee[d]=hptRegRatio_bbmctee[d]->GetMean(1);
    MEANE_bbrmm[d]=hptRatio_bbmm[d]->GetMeanError(1); 
    MEANE_bbree[d]=hptRatio_bbee[d]->GetMeanError(1);
    MEANE_bbrmcmm[d]=hptRatio_bbmcmm[d]->GetMeanError(1);
    MEANE_bbrmcee[d]=hptRatio_bbmcee[d]->GetMeanError(1);
    MEANE_bbrmctmm[d]=hptRatio_bbmctmm[d]->GetMeanError(1);
    MEANE_bbrmctee[d]=hptRatio_bbmctee[d]->GetMeanError(1);
    MEANE_bbrrmm[d]=hptRegRatio_bbmm[d]->GetMeanError(1); 
    MEANE_bbrree[d]=hptRegRatio_bbee[d]->GetMeanError(1);
    MEANE_bbrrmcmm[d]=hptRegRatio_bbmcmm[d]->GetMeanError(1);
    MEANE_bbrrmcee[d]=hptRegRatio_bbmcee[d]->GetMeanError(1);
    MEANE_bbrrmctmm[d]=hptRegRatio_bbmctmm[d]->GetMeanError(1);
    MEANE_bbrrmctee[d]=hptRegRatio_bbmctee[d]->GetMeanError(1);

    //in order to plot the RMS(mean) vs pt(Z) it is important to add the to MC components

    /*hptRatio_bbmctee[d]->Add(hptRatio_bbmcee[d],(hptRatio_bbee[d]->Integral()-hptRatio_bbmctee[d]->Integral()));
    RMS_MCee[d]=   hptRatio_bbmctee[d]->GetRMS(1);
    RMSE_MCee[d]=   hptRatio_bbmctee[d]->GetRMSError(1);
    MEAN_MCee[d]=   hptRatio_bbmctee[d]->GetMean(1);
    MEANE_MCee[d]=   hptRatio_bbmctee[d]->GetMeanError(1); 
    hptRegRatio_bbmctee[d]->Add(hptRegRatio_bbmcee[d],(hptRegRatio_bbee[d]->Integral()-hptRegRatio_bbmctee[d]->Integral()));
    RMS_MCree[d]=   hptRegRatio_bbmctee[d]->GetRMS(1);
    RMSE_MCree[d]=   hptRegRatio_bbmctee[d]->GetRMSError(1);
    MEAN_MCree[d]=   hptRegRatio_bbmctee[d]->GetMean(1);
    MEANE_MCree[d]=   hptRegRatio_bbmctee[d]->GetMeanError(1);
    hptRatio_bbmctmm[d]->Add(hptRatio_bbmcmm[d],(hptRatio_bbmm[d]->Integral()-hptRatio_bbmctmm[d]->Integral()));
    RMS_MCmm[d]=   hptRatio_bbmctmm[d]->GetRMS(1);
    RMSE_MCmm[d]=   hptRatio_bbmctmm[d]->GetRMSError(1);
    MEAN_MCmm[d]=   hptRatio_bbmctmm[d]->GetMean(1);
    MEANE_MCmm[d]=   hptRatio_bbmctmm[d]->GetMeanError(1);
    hptRegRatio_bbmctmm[d]->Add(hptRegRatio_bbmcmm[d],(hptRegRatio_bbmm[d]->Integral()-hptRegRatio_bbmctmm[d]->Integral()));
    RMS_MCrmm[d]=   hptRegRatio_bbmctmm[d]->GetRMS(1);
    RMSE_MCrmm[d]=   hptRegRatio_bbmctmm[d]->GetRMSError(1);
    MEAN_MCrmm[d]=   hptRegRatio_bbmctmm[d]->GetMean(1);
    MEANE_MCrmm[d]=   hptRegRatio_bbmctmm[d]->GetMeanError(1);*/

  }

  gStyle->SetOptStat(0);

  //plot pTratio canvas
  char name_rmm[100], name_ree[100], name_rrmm[100], name_rree[100];
  char name_bbrmm[100], name_bbree[100], name_bbrrmm[100], name_bbrree[100];
 
  for(int b=0; b<11; b++){
    sprintf(name_bbrmm, "ptRatio_bbmm%d", b);
    //histo_ptRatio(name_bbrmm, hptRatio_bbmm[b], hptRatio_bbmcmm[b], hptRatio_bbmctmm[b]);
    sprintf(name_bbree, "ptRatio_bbee%d", b);
    //histo_ptRatio(name_bbree,  hptRatio_bbee[b], hptRatio_bbmcee[b], hptRatio_bbmctee[b]);
    sprintf(name_bbrrmm, "ptRatio_bbrrmm%d", b);
    //histo_ptRatio(name_bbrrmm, hptRegRatio_bbmm[b], hptRegRatio_bbmcmm[b], hptRegRatio_bbmctmm[b]);
    sprintf(name_bbrree, "ptRatio_bbrree%d", b);
    //histo_ptRatio(name_bbrree, hptRegRatio_bbee[b], hptRegRatio_bbmcee[b], hptRegRatio_bbmctee[b]);
  }
  
  //plot RMS and MEAN canvas

  //RMS("RMSbbee ", "rms_bbee", "rms_bbmcee", RMS_bbree, RMS_MCee, RMSE_bbree, RMSE_MCee, ptRatiobins, true,false,false);
  //RMS("RMSbbmm ", "rms_bbmm", "rms_bbmcmm", RMS_bbrmm, RMS_MCmm, RMSE_bbrmm, RMSE_MCmm, ptRatiobins, true,false,false);
  //RMS("RMSrbbee ", "rms_rbbee", "rms_rbbmcee", RMS_bbree, RMS_MCree, RMSE_bbree, RMSE_MCree, ptRatiobins, true,false,false);
  //RMS("RMSrbbmm ", "rms_rbbmm", "rms_rbbmcmm", RMS_bbrmm, RMS_MCrmm, RMSE_bbrmm, RMSE_MCrmm, ptRatiobins, true,false,false);

  // RMS("MEANbbee ", "mean_bbee", "mean_bbmcee", MEAN_bbree, MEAN_MCee, MEANE_bbree, MEANE_MCee, ptRatiobins, false,false,false);
  // RMS("MEANbbmm ", "mean_bbmm", "mean_bbmcmm", MEAN_bbrmm, MEAN_MCmm, MEANE_bbrmm, MEANE_MCmm, ptRatiobins, false,false,false);
  // RMS("MEANrbbee ", "mean_rbbee", "mean_rbbmcee", MEAN_bbree, MEAN_MCree, MEANE_bbree, MEANE_MCree, ptRatiobins, false,false,false);
  //RMS("MEANrbbmm ", "mean_rbbmm", "mean_rbbmcmm", MEAN_bbrmm, MEAN_MCrmm, MEANE_bbrmm, MEANE_MCrmm, ptRatiobins, false,false,false);

  //plot pt(bb system)/pt(Z) (only data) comparing before and after regression
  
  //ptbalanceFIT("ptbalance_eeFIT", hptbalance_ee, hptbalance_ree);
  //ptbalanceFIT("ptbalance_mmFIT", hptbalance_mm, hptbalance_rmm);

  // plot bb system invarian mass (only data) comparing before and after regression

  //ptbalance("Mbb_mm", hmbb_mm, hmbb_rmm);
  //ptbalance("Mbb_ee", hmbb_ee, hmbb_ree);

  //plot ptratio and mbb  comparing montecarlo and data

  /* ptbalance_sep("mbb_ee", hmbb_ee, hmbb_mctee, hmbb_mcee, true, false);
   ptbalance_sep("mbb_ree", hmbb_ree, hmbb_rmctee, hmbb_rmcee, true, true);
   ptbalance_sep("mbb_mm", hmbb_mm, hmbb_mctmm, hmbb_mcmm, true, false);
   ptbalance_sep("mbb_rmm", hmbb_rmm, hmbb_rmctmm, hmbb_rmcmm, true, true);
   ptbalance_sep("ptbalance_ee", hptbalance_ee, hptbalance_mctee, hptbalance_mcee, false, false);
   ptbalance_sep("ptbalance_mm", hptbalance_mm, hptbalance_mctmm, hptbalance_mcmm, false, false);
   ptbalance_sep("ptbalance_ree", hptbalance_ree, hptbalance_rmctee, hptbalance_rmcee, false, true);
   ptbalance_sep("ptbalance_rmm", hptbalance_rmm, hptbalance_rmctmm, hptbalance_rmcmm, false, true);*/

}

//function to plot mbb and ptratio with data, dy and ttbar
void ptbalance_sep(TString HISTO, TH1F *h_data,  TH1F *h_mct, TH1F *h_mc, bool ismbb, bool isreg){
  TCanvas *canvas = new TCanvas("canvas_"+HISTO, "canvas_"+HISTO, 900, 600);
  TPad *pad1 = new TPad("pad1"+HISTO, "pad1"+HISTO, 0., 0.3, 1.,1.);
  pad1->SetBottomMargin(0.);
  pad1->Draw();
  pad1->cd();
  h_mct->SetLineColor(kBlack);
  h_mct->SetFillColor(kGreen-6);
  h_mc->SetLineColor(kBlack);
  h_mc->SetFillColor(kRed-10);
  THStack *hs = new THStack("hs","");
  hs->Add(h_mct);
  hs->Add(h_mc);
  hs->Draw(" histo");
  canvas->Update();
  hs->SetMinimum(-7);
  if(ismbb==true) hs->SetMaximum(300);
  if(ismbb==false) hs->SetMaximum(400);
  hs->GetYaxis()->SetTitle("Events");
  hs->SetTitle(""); 
  hs->GetYaxis()->SetLabelSize(0.04);
  hs->GetYaxis()->SetTitleSize(0.06);
  hs->GetYaxis()->SetTitleOffset(0.7);
  h_data->SetLineColor(kBlack);
  h_data->SetMarkerStyle(20);
  h_data->SetMarkerSize(1.0);
  h_data->DrawCopy("same ep");
  leg = new TLegend(0.63,0.65,0.78,0.85); 
  leg->SetTextSize(0.03);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->AddEntry(h_data,"Data","lep");
  leg->AddEntry(h_mc,"MC DY","f");
  leg->AddEntry(h_mct,"MC TT","f");
  leg->Draw();
  gPad->SetTickx();
  gPad->SetTicky();
  canvas->cd();
  TPad *pad2 = new TPad("pad2","pad2",0,0.,1,0.3);
  pad2->SetBottomMargin(0.4);
  pad2->SetTopMargin(0.);
  pad2->SetGridy();
  pad2->Draw();
  pad2->cd();
  h_data->Sumw2();
  h_data->SetStats(0);
  h_mc->Add(h_mct);
  h_data->Divide((h_mc));
  h_data->SetMarkerStyle(20);
  h_data->GetYaxis()->SetTitle("DATA/MC");
  if(isreg==true && ismbb==false)h_data->GetXaxis()->SetTitle("Regressed p_{T}^{b#bar{b}}/p_{T}^{Z}");
  if(isreg==false && ismbb==false)h_data->GetXaxis()->SetTitle("Unregressed p_{T}^{b#bar{b}}/p_{T}^{Z}");
  if(isreg==true && ismbb==true)h_data->GetXaxis()->SetTitle("Regressed m_{b#barb} [GeV]");
  if(isreg==false && ismbb==true)h_data->GetXaxis()->SetTitle("Unregressed m_{b#barb} [GeV]");
  h_data->GetYaxis()->SetLabelSize(0.1);
  h_data->GetYaxis()->SetTitleSize(0.1);
  h_data->GetXaxis()->SetLabelSize(0.1);
  h_data->GetXaxis()->SetTitleSize(0.1);
  h_data->GetYaxis()->SetTitleOffset(0.4);
  h_data->GetXaxis()->SetTitleOffset(1);
  h_data->GetYaxis()->SetNdivisions(505,-1);
  h_data->SetMinimum(-0.15);
  h_data->SetMaximum(2.15);
  h_data->Draw("ep");
  h_data->SetTitle("");
  gPad->SetTickx();
  gPad->SetTicky(); 
  canvas->cd();
  canvas->Print(HISTO+".pdf");
  }

//function to fit with a gaussian in the interval [0.7,1.3] the ptratio (or with a 3-iterative gaussian fit), putting together unregressed and regressed data
void ptbalanceFIT(TString HISTO, TH1F *h_bef, TH1F *h_aft){
  TCanvas *canvas = new TCanvas("canvas_"+HISTO, "canvas_"+HISTO, 900, 900);
  TF1* Gaussian_bef = new TF1("Gaussian_bef","gaus",0.7,1.3);
  Gaussian_bef->SetLineStyle(1);
  Gaussian_bef->SetLineWidth(3);
  Gaussian_bef->SetLineColor(kGray+2);
  TF1* Gaussian_aft = new TF1("Gaussian_aft","gaus",0.7,1.3);
  Gaussian_aft->SetLineStyle(1);
  Gaussian_aft->SetLineWidth(3);
  Gaussian_aft->SetLineColor(kRed);
  h_aft->Fit("Gaussian_aft", "", "", 0.7, 1.3);
  h_bef->Fit("Gaussian_bef", "", "", 0.7, 1.3);
  /* double sigma_aft  = Gaussian_aft->GetParameter(2);
  double sigmaE_aft = Gaussian_aft->GetParError(2);
  double sigma_bef = Gaussian_bef->GetParameter(2);
  double sigmaE_bef = Gaussian_bef->GetParError(2);
  double  mu_aft = Gaussian_aft->GetParameter(1);
  double muE_aft = Gaussian_aft->GetParError(1);
  double mu_bef = Gaussian_bef->GetParameter(1);
  double muE_bef = Gaussian_bef->GetParError(1);
  h_aft->Fit("Gaussian_aft", "", "",mu_aft - 1.5*sigma_aft, mu_aft + 1.5*sigma_aft);
  h_bef->Fit("Gaussian_bef", "", "", mu_bef - 1.5*sigma_bef, mu_bef + 1.5*sigma_bef);
  sigma_aft  = Gaussian_aft->GetParameter(2);
  sigmaE_aft = Gaussian_aft->GetParError(2);
  sigma_bef = Gaussian_bef->GetParameter(2);
  sigmaE_bef = Gaussian_bef->GetParError(2);
  mu_aft = Gaussian_aft->GetParameter(1);
  muE_aft = Gaussian_aft->GetParError(1);
  mu_bef = Gaussian_bef->GetParameter(1);
  muE_bef = Gaussian_bef->GetParError(1);
  h_aft->Fit("Gaussian_aft", "", "",mu_aft - 1.5*sigma_aft, mu_aft + 1.5*sigma_aft);
  h_bef->Fit("Gaussian_bef", "", "", mu_bef - 1.5*sigma_bef, mu_bef + 1.5*sigma_bef);
  cout << "***********************************ECCOCI******************************************"<< endl;
  cout << Gaussian_bef->GetParameter(0) << "+/-" << Gaussian_bef->GetParError(0) << "\t *-*-" << Gaussian_aft->GetParameter(0) << "+/-" << Gaussian_aft->GetParError(0) << endl;
  cout << sigma_bef << " +/- " << sigmaE_bef << "\t -*-*-"<< sigma_aft << " +/- "  << sigmaE_aft << endl;
  cout << mu_bef  << " +/- " << muE_bef << "\t -*-*-"<< mu_aft << " +/- " << muE_aft << endl;
  cout << sigma_bef/mu_bef << "//////////" << sigma_aft/mu_aft << endl;
  cout << "*****************************************************************************"<< endl;*/
  h_bef->GetYaxis()->SetTitle("Events");
  h_bef->GetXaxis()->SetTitle("p_{T}^{b#bar{b}}/p_{T}^{Z}");
  h_bef->SetTitle("");
  h_bef->SetLineColor(kGray);
  h_bef->SetFillColor(kGray);
  h_bef->Draw("histo");
  h_bef->GetXaxis()->SetTitleOffset(1.2);
  h_bef->GetYaxis()->SetTitleOffset(1.2);
  h_aft->SetLineColor(kRed);
  h_aft->Draw("same hist");
  Gaussian_bef->Draw("same");
  Gaussian_aft->Draw("same");
  leg = new TLegend(0.65,0.65,0.80,0.85); 
  leg->SetTextSize(0.03);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->AddEntry(h_aft,"Regressed","l");
  leg->AddEntry(h_bef,"Baseline","f");
  leg->Draw();
  gPad->SetTickx();
  gPad->SetTicky();
  canvas->cd();
  canvas->Print(HISTO+".pdf");
}

// function to plot unregressed and regressed ptratio (not fitting it)
void ptbalance(TString HISTO, TH1F *h_bef, TH1F *h_aft){
  TCanvas *canvas = new TCanvas("canvas_"+HISTO, "canvas_"+HISTO, 900, 900);
  h_bef->GetYaxis()->SetTitle("Events");
  h_bef->GetXaxis()->SetTitle("m_{b#bar{b}} [GeV]");
  h_bef->SetTitle("");
  h_bef->SetLineColor(kGray);
  h_bef->SetFillColor(kGray);
  h_bef->Draw("histo");
  h_bef->GetXaxis()->SetTitleOffset(1.2);
  h_bef->GetYaxis()->SetTitleOffset(1.2);
  h_aft->SetLineColor(kRed);
  h_aft->Draw("same histo");
  leg = new TLegend(0.65,0.65,0.80,0.85); 
  leg->SetTextSize(0.03);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->AddEntry(h_aft,"Regressed","l");
  leg->AddEntry(h_bef,"Baseline","f");
  leg->Draw();
  gPad->SetTickx();
  gPad->SetTicky();
  canvas->cd();
  canvas->Print(HISTO+".pdf");
}

//function to plot ptratio in each pt(Z) interval
void histo_ptRatio(TString HISTO, TH1F *h_data, TH1F *h_mc, TH1F *h_mct){
  TCanvas *canvas = new TCanvas("canvas_"+HISTO, "canvas_"+HISTO, 900, 600);
  h_mct->SetTitle("");
  h_mct->SetLineColor(kBlack);
  h_mct->SetFillColor(kGreen-6);
  h_mc->SetLineColor(kBlack);
  h_mc->SetFillColor(kRed-10);
  THStack *hs1 = new THStack("hs1","");
  hs1->Add(h_mct);
  hs1->Add(h_mc);
  hs1->Draw("");
  canvas->Update();
  hs1->GetYaxis()->SetTitle("Events");
  hs1->GetXaxis()->SetTitle("Response");
  h_data->SetMarkerSize(0.8);
  h_data->SetMarkerStyle(20);
  h_data->Draw("same ep");
  leg = new TLegend(0.73,0.65,0.88,0.85); 
  leg->SetTextSize(0.05);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->AddEntry(h_data,"Data","lep");
  leg->AddEntry(h_mc,"DY MC","f");
  leg->AddEntry(h_mct, "TT MC", "f");
  leg->Draw();
  gPad->SetTickx();
  gPad->SetTicky();
  canvas->cd();
  canvas->Print(HISTO+".pdf");
}

//function to plot RMS and MEAN of ptratio distribution
void RMS(TString HISTO, TString HISTO_DATA, TString HISTO_MC, float *v_data, float *v_mc,float *ve_data, float *ve_mc,float *bins,bool isrms, bool islog, bool ismbb) {
  TH1F *h_data = new TH1F(HISTO_DATA,HISTO_DATA,11, bins);
  TH1F *h_mc = new TH1F(HISTO_MC,HISTO_MC,11, bins);
  for(int i=0; i<11; i++){
    h_data->SetBinContent(i+1,v_data[i]);
    h_data->SetBinError(i+1,ve_data[i]);
    h_mc->SetBinContent(i+1,v_mc[i]);
    h_mc->SetBinError(i+1,ve_mc[i]);
  }
  TCanvas *canvas = new TCanvas("canvas_"+HISTO, "canvas_"+HISTO, 900, 700);
  TPad *pad1 = new TPad("pad1"+HISTO, "pad1"+HISTO, 0., 0.3, 1.,1.);
  pad1->SetBottomMargin(0.);
  pad1->Draw();
  pad1->cd();
  h_data->SetTitle("");
  if(ismbb==false){
    if(isrms==true) h_data->GetYaxis()->SetTitle("Unregressed RMS(p_{T}^{bb}/p_{T}^{Z})");
    if(isrms==false) h_data->GetYaxis()->SetTitle("Unregressed Mean(p_{T}^{bb}/p_{T}^{Z})");
  }  
  if(ismbb==true){
    if(isrms==true) h_data->GetYaxis()->SetTitle("Unregressed RMS(M_{bb}) [GeV]");
    if(isrms==false) h_data->GetYaxis()->SetTitle("Unregressed Mean(M_{bb}) [GeV]");
  }
  h_data->SetMinimum(0.83);
  h_data->SetMaximum(1.13);
  h_data->GetYaxis()->SetLabelSize(0.04);
  h_data->GetYaxis()->SetTitleSize(0.06);
  h_data->GetYaxis()->SetTitleOffset(0.7);
  h_data->SetMarkerColor(kBlue);
  h_data->SetMarkerSize(0.8);
  h_data->SetMarkerStyle(20);
  h_data->DrawCopy("ep");
  h_mc->SetMarkerColor(kMagenta);
  h_mc->SetMarkerSize(0.8);
  h_mc->SetMarkerStyle(21);
  h_mc->DrawCopy("same ep");
 
  if(ismbb==false){
    leg = new TLegend(0.73,0.65,0.88,0.85); 
    leg->SetTextSize(0.05);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->AddEntry(h_data,"Data","lep");
    leg->AddEntry(h_mc," MC","lep");
    leg->Draw();
  }
  if(ismbb==true){
    leg = new TLegend(0.91,0.60,0.99,0.80);
    leg->SetTextSize(0.04);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->AddEntry(h_data,"Data","lep");
    leg->AddEntry(h_mc," MC","lep");
    leg->Draw();
  }
  gPad->SetLogx();
  gPad->SetTickx();
  gPad->SetTicky();
  canvas->cd();
  TPad *pad2 = new TPad("pad2","pad2",0,0.,1,0.3);
  pad2->SetBottomMargin(0.4);
  pad2->SetTopMargin(0.);
  pad2->SetGridy();
  pad2->Draw();
  pad2->cd();
  h_data->SetStats(0);
  h_data->Divide(h_mc);
  h_data->SetMarkerStyle(20);
  h_data->GetYaxis()->SetTitle("DATA/MC");
  h_data->GetXaxis()->SetTitle("p_{T}^{Z} [GeV]");
  h_data->GetYaxis()->SetLabelSize(0.1);
  h_data->GetYaxis()->SetTitleSize(0.1);
  h_data->GetXaxis()->SetLabelSize(0.1);
  h_data->GetXaxis()->SetTitleSize(0.12);
  h_data->GetYaxis()->SetTitleOffset(0.4);
  h_data->GetXaxis()->SetTitleOffset(1);
  h_data->GetYaxis()->SetNdivisions(505,-1);
  h_data->GetXaxis()->SetNoExponent(1);
  if(isrms==false){
    h_data->SetMinimum(0.77);
    h_data->SetMaximum(1.23);
    }  
     if(isrms==true){
    h_data->SetMinimum(0.18);
    h_data->SetMaximum(1.82);
    }
  h_data->Draw("ep");
  h_data->SetTitle("");
  gPad->SetLogx();
  gPad->SetTickx();
  gPad->SetTicky();
  canvas->cd();
  canvas->Print(HISTO+".pdf");
}

