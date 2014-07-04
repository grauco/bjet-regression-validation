/*****************************************
 *                                       *
 *           Giorgia Rauco               *
 *           July 2014, 7th              *
 *                                       *
 *              VBF Hbb                  * 
 *                                       *
 *  Validation of bjet regression in     * 
 *  Z(->dielectron, dimuon) + 1bjet      * 
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
#include "THStack"

#include <math.h>
#include <stdio.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <stdlib.h>

void Z1bjet(){

  //histograms declaration

  TH1F* hjetbtag_mm = new TH1F("jetbtag_mm", "jetbtag_mm", 80, 0., 1);
  TH1F* hjetbtag_ee = new TH1F("jetbtag_ee", "jetbtag_ee", 80, 0., 1);
  TH1F* hjetbtag_mcmm = new TH1F("jetbtag_mcmm", "jetbtag_mcmm", 80, 0., 1);
  TH1F* hjetbtag_mcee = new TH1F("jetbtag_mcee", "jetbtag_mcee", 80, 0., 1);  
  TH1F* hjetbtag_mctmm = new TH1F("jetbtag_mctmm", "jetbtag_mctmm", 80, 0., 1);
  TH1F* hjetbtag_mctee = new TH1F("jetbtag_mctee", "jetbtag_mctee", 80, 0., 1); 

  TH1F* hmll_DoubleMu = new TH1F("DoubleMu", "mll_DoubleMu", 100, 60, 160);
  TH1F* hmll_DoubleElectron = new TH1F("DoubleElectron", "mll_DoubleElectron", 100, 60, 160); 
  TH1F* hmll_DYLLJetsee = new TH1F("DYLLJetsDoubleElectron", "mll_DYLLJetsee", 100, 60, 160);
  TH1F* hmll_DYLLJetsmm = new TH1F("DYLLJetsDoubleMuon", "mll_DYLLJetsee", 100, 60, 160);
  TH1F* hmll_TTJetsee = new TH1F("TTJetsDoubleElectron", "mll_TTJetsee", 100, 60, 160);
  TH1F* hmll_TTJetsmm = new TH1F("TTJetsDoubleMuon", "mll_TTJetsee", 100, 60, 160);

  //validation of input variables

  TH1F* hjetpt_ee = new TH1F("jetpt_ee", "jetpt_ee", 80, 0, 400);
  TH1F* hjetpt_mm = new TH1F("jetpt_mm", "jetpt_mm", 80, 0, 400);
  TH1F* hjetpt_mcmm = new TH1F("jetpt_mcmm", "jetpt_mcmm", 80,0, 400);
  TH1F* hjetpt_mcee = new TH1F("jetpt_mcee", "jetpt_mcee", 80,0, 400);
  TH1F* hjetpt_mctmm = new TH1F("jetpt_mctmm", "jetpt_mctmm", 80,0, 400);
  TH1F* hjetpt_mctee = new TH1F("jetpt_mctee", "jetpt_mctee", 80,0, 400);

  TH1F* hjeteta_ee = new TH1F("jeteta_ee", "jeteta_ee", 80, -4, 4);
  TH1F* hjeteta_mm = new TH1F("jeteta_mm", "jeteta_mm", 80, -4, 4);
  TH1F* hjeteta_mcmm = new TH1F("jeteta_mcmm", "jeteta_mcmm", 80,-4, 4);
  TH1F* hjeteta_mcee = new TH1F("jeteta_mcee", "jeteta_mcee", 80,-4, 4);
  TH1F* hjeteta_mctmm = new TH1F("jeteta_mctmm", "jeteta_mctmm", 80,-4, 4);
  TH1F* hjeteta_mctee = new TH1F("jeteta_mctee", "jeteta_mctee", 80,-4, 4);

  TH1F* hjetnhf_ee = new TH1F("jetnhf_ee", "jetnhf_ee", 100, 0, 1);
  TH1F* hjetnhf_mm = new TH1F("jetnhf_mm", "jetnhf_mm", 100, 0, 1);
  TH1F* hjetnhf_mcee = new TH1F("jetnhf_mcee", "jetnhf_mcee", 100, 0, 1);
  TH1F* hjetnhf_mcmm = new TH1F("jetnhf_mcmm", "jetnhf_mcmm", 100, 0, 1);
  TH1F* hjetnhf_mctee = new TH1F("jetnhf_mctee", "jetnhf_mctee", 100, 0, 1);
  TH1F* hjetnhf_mctmm = new TH1F("jetnhf_mctmm", "jetnhf_mctmm", 100, 0, 1);
  
  TH1F* hjetphf_ee = new TH1F("jetphf_ee", "jetphf_ee", 100, 0, 1);
  TH1F* hjetphf_mm = new TH1F("jetphf_mm", "jetphf_mm", 100, 0, 1);
  TH1F* hjetphf_mcee = new TH1F("jetphf_mcee", "jetphf_mcee", 100, 0, 1);
  TH1F* hjetphf_mcmm = new TH1F("jetphf_mcmm", "jetphf_mcmm", 100, 0, 1);
  TH1F* hjetphf_mctee = new TH1F("jetphf_mctee", "jetphf_mctee", 100, 0, 1);
  TH1F* hjetphf_mctmm = new TH1F("jetphf_mctmm", "jetphf_mctmm", 100, 0, 1);

  TH1F* hjetmetphi_ee = new TH1F("jetmetphi_ee", "jetmetphi_ee", 63, 0, 3.15);
  TH1F* hjetmetphi_mm = new TH1F("jetmetphi_mm", "jetmetphi_mm", 63, 0, 3.15);
  TH1F* hjetmetphi_mcee = new TH1F("jetmetphi_mcee", "jetmetphi_mcee", 63, 0, 3.15);
  TH1F* hjetmetphi_mcmm = new TH1F("jetmetphi_mcmm", "jetmetphi_mcmm", 63, 0, 3.15);
  TH1F* hjetmetphi_mctee = new TH1F("jetmetphi_mctee", "jetmetphi_mctee", 63, 0, 3.15);
  TH1F* hjetmetphi_mctmm = new TH1F("jetmetphi_mctmm", "jetmetphi_mctmm", 63, 0, 3.15);
 
  TH1F* hmet_ee = new TH1F("met_ee", "met_ee", 100, 0, 100);
  TH1F* hmet_mm = new TH1F("met_mm", "met_mm", 100, 0, 100);
  TH1F* hmet_mcee = new TH1F("met_mcee", "met_mcee", 100, 0, 100);
  TH1F* hmet_mcmm = new TH1F("met_mcmm", "met_mcmm", 100, 0, 100);
  TH1F* hmet_mctee = new TH1F("met_mctee", "met_mctee", 100, 0, 100);
  TH1F* hmet_mctmm = new TH1F("met_mctmm", "met_mctmm", 100, 0, 100);

  TH1F* hrho_ee = new TH1F("rho_ee", "rho_ee", 70, 0, 35);
  TH1F* hrho_mm = new TH1F("rho_mm", "rho_mm", 70, 0, 35);
  TH1F* hrho_mcee = new TH1F("rho_mcee", "rho_mcee", 70, 0, 35);
  TH1F* hrho_mcmm = new TH1F("rho_mcmm", "rho_mcmm", 70, 0, 35);
  TH1F* hrho_mctee = new TH1F("rho_mctee", "rho_mctee", 70, 0, 35);
  TH1F* hrho_mctmm = new TH1F("rho_mctmm", "rho_mctmm", 70, 0, 35);

  TH1F *hjetvtx3del_ee = new TH1F("jetvtx3del_ee", "jetvtx3del_ee", 75, 0, 1.5);
  TH1F *hjetvtx3del_mm = new TH1F("jetvtx3del_mm", "jetvtx3del_mm", 75, 0, 1.5);
  TH1F *hjetvtx3del_mcee = new TH1F("jetvtx3del_mcee", "jetvtx3del_mcee", 75, 0, 1.5);
  TH1F *hjetvtx3del_mcmm = new TH1F("jetvtx3del_mcmm", "jetvtx3del_mcmm", 75, 0, 1.5);
  TH1F *hjetvtx3del_mctee = new TH1F("jetvtx3del_mctee", "jetvtx3del_mctee", 75, 0, 1.5);
  TH1F *hjetvtx3del_mctmm = new TH1F("jetvtx3del_mctmm", "jetvtx3del_mctmm", 75, 0, 1.5);

  TH1F *hjetmass_ee = new TH1F("jetmass_ee", "jetmass_ee", 100, 0, 100);
  TH1F *hjetmass_mm = new TH1F("jetmass_mm", "jetmass_mm", 100, 0, 100);
  TH1F *hjetmass_mcee = new TH1F("jetmass_mcee", "jetmass_mcee", 100, 0, 100);
  TH1F *hjetmass_mcmm = new TH1F("jetmass_mcmm", "jetmass_mcmm", 100, 0, 100);
  TH1F *hjetmass_mctee = new TH1F("jetmass_mctee", "jetmass_mctee", 100, 0, 100);
  TH1F *hjetmass_mctmm = new TH1F("jetmass_mctmm", "jetmass_mctmm", 100, 0, 100);

  TH1F *hjetvtxmass_ee = new TH1F("jetvtxmass_ee", "jetvtxmass_ee", 50, 0, 5);
  TH1F *hjetvtxmass_mm = new TH1F("jetvtxmasss_mm", "jetvtxmass_mm", 50, 0, 5);
  TH1F *hjetvtxmass_mcee = new TH1F("jetvtxmass_mcee", "jetvtxmass_mcee", 50, 0, 5);
  TH1F *hjetvtxmass_mcmm = new TH1F("jetvtxmass_mcmm", "jetvtxmass_mcmm", 50, 0, 5);
  TH1F *hjetvtxmass_mctee = new TH1F("jetvtxmass_mctee", "jetvtxmass_mctee", 50, 0, 5);
  TH1F *hjetvtxmass_mctmm = new TH1F("jetvtxmass_mctmm", "jetvtxmass_mctmm", 50, 0, 5);

  TH1F *hjetsoftleppt_ee = new TH1F("jetsoftleppt_ee", "jetsoftleppt_ee", 70, 0,70);
  TH1F *hjetsoftleppt_mm = new TH1F("jetsoftleppt_mm", "jetsoftleppt_mm", 70, 0, 70);
  TH1F *hjetsoftleppt_mcee = new TH1F("jetsoftleppt_mcee", "jetsoftleppt_mcee", 70, 0, 70);
  TH1F *hjetsoftleppt_mcmm = new TH1F("jetsoftleppt_mcmm", "jetsoftleppt_mcmm", 70, 0, 70);
  TH1F *hjetsoftleppt_mctee = new TH1F("jetsoftleppt_mctee", "jetsoftleppt_mctee", 70, 0, 70);
  TH1F *hjetsoftleppt_mctmm = new TH1F("jetsoftleppt_mctmm", "jetsoftleppt_mctmm", 70, 0, 70);

  TH1F *hjetsoftlepptrel_ee = new TH1F("jetsoftlepptrel_ee", "jetsoftlepptrel_ee", 40, 0, 20);
  TH1F *hjetsoftlepptrel_mm = new TH1F("jetsoftlepptrel_mm", "jetsoftlepptrel_mm", 40, 0, 20);
  TH1F *hjetsoftlepptrel_mcee = new TH1F("jetsoftlepptrel_mcee", "jetsoftlepptrel_mcee", 40, 0, 20);
  TH1F *hjetsoftlepptrel_mcmm = new TH1F("jetsoftlepptrel_mcmm", "jetsoftlepptrel_mcmm", 40, 0, 20);
  TH1F *hjetsoftlepptrel_mctee = new TH1F("jetsoftlepptrel_mctee", "jetsoftlepptrel_mctee", 40, 0, 20);
  TH1F *hjetsoftlepptrel_mctmm = new TH1F("jetsoftlepptrel_mctmm", "jetsoftlepptrel_mctmm", 40, 0, 20);

  TH1F *hjetleadtrkpt_ee = new TH1F("jetleadtrkpt_ee", "jetleadtrkpt_ee", 50, 0, 100);
  TH1F *hjetleadtrkpt_mm = new TH1F("jetleadtrkpt_mm", "jetleadtrkpt_mm", 50, 0, 100);
  TH1F *hjetleadtrkpt_mcee = new TH1F("jetleadtrkpt_mcee", "jetleadtrkpt_mcee", 50, 0, 100);
  TH1F *hjetleadtrkpt_mcmm = new TH1F("jetleadtrkpt_mcmm", "jetleadtrkpt_mcmm", 50, 0, 100);
  TH1F *hjetleadtrkpt_mctee = new TH1F("jetleadtrkpt_mctee", "jetleadtrkpt_mctee", 50, 0, 100);
  TH1F *hjetleadtrkpt_mctmm = new TH1F("jetleadtrkpt_mctmm", "jetleadtrkpt_mctmm", 50, 0, 100);

  TH1F *hjetpart_ee = new TH1F("jetpart_ee", "jetpart_ee", 70, 0, 70);
  TH1F *hjetpart_mm = new TH1F("jetpart_mm", "jetpart_mm", 70, 0, 70);
  TH1F *hjetpart_mcee = new TH1F("jetpart_mcee", "jetpart_mcee", 70, 0, 70);
  TH1F *hjetpart_mcmm = new TH1F("jetpart_mcmm", "jetpart_mcmm", 70, 0, 70);
  TH1F *hjetpart_mctee = new TH1F("jetpart_mctee", "jetpart_mctee", 70, 0, 70);
  TH1F *hjetpart_mctmm = new TH1F("jetpart_mctmm", "jetpart_mctmm", 70, 0, 70);


  // pt(jet)/pt(Z) histograms declaration (Unregressed and Regressed)
  int a;
  const unsigned int N = 11; //11 intervals of pt(Z) range
  TH1F *hptRatio_mm[N] ;
  TH1F *hptRatio_ee[N] ;
  TH1F *hptRatio_mcmm[N];
  TH1F *hptRatio_mcee[N] ;
  TH1F *hptRatio_mctmm[N];
  TH1F *hptRatio_mctee[N] ;
  TH1F *hptRegRatio_mm[N];
  TH1F *hptRegRatio_ee[N];
  TH1F *hptRegRatio_mcmm[N];
  TH1F *hptRegRatio_mcee[N];
  TH1F *hptRegRatio_mctmm[N];
  TH1F *hptRegRatio_mctee[N];

  char namermm[100], nameree[100], namermcmm[100], namermcee[100],  namermctmm[100], namermctee[100], namerrmm[100], namerree[100], namerrmcmm[100], namerrmcee[100],  namerrmctmm[100], namerrmctee[100];

  for(a=0; a<N; a++){
    //declare unregressed pTratio histograms
    sprintf(namermm, "ptRatio_mm%d", a);
    hptRatio_mm[a] = new TH1F(namermm,namermm,40,0,2);
    printf(nameree, "ptRatio_ee%d", a);
    hptRatio_ee[a]  = new TH1F(nameree,nameree,40,0,2);
    printf(namermcmm, "ptRatio_mcmm%d", a);
    hptRatio_mcmm[a] = new TH1F(namermcmm,namermcmm,40,0,2);
    printf(namermcee, "ptRatio_mcee%d", a);
    hptRatio_mcee[a]  = new TH1F(namermcee,namermcee,40,0,2);
    printf(namermctmm, "ptRatio_mctmm%d", a);
    hptRatio_mctmm[a] = new TH1F(namermctmm,namermctmm,40,0,2);
    printf(namermctee, "ptRatio_mctee%d", a);
    hptRatio_mctee[a]  = new TH1F(namermctee,namermctee,40,0,2);

    //declare regressed pTratio histograms
    printf(namerrmm, "ptRegRatio_mm%d", a);
    hptRegRatio_mm[a] = new TH1F(namerrmm,namerrmm,40,0,2);
    printf(namerree, "ptRegRatio_ee%d", a);
    hptRegRatio_ee[a] = new TH1F(namerree,namerree,40,0,2);
    printf(namerrmcmm, "ptRegRatio_mcmm%d", a);
    hptRegRatio_mcmm[a] = new TH1F(namerrmcmm,namerrmcmm,40,0,2);
    printf(namerrmcee, "ptRegRatio_mcee%d", a);
    hptRegRatio_mcee[a] = new TH1F(namerrmcee,namerrmcee,40,0,2);
    printf(namerrmctmm, "ptRegRatio_mctmm%d", a);
    hptRegRatio_mctmm[a] = new TH1F(namerrmctmm,namerrmctmm,40,0,2);
    printf(namerrmctee, "ptRegRatio_mctee%d", a);
    hptRegRatio_mctee[a] = new TH1F(namerrmctee,namerrmctee,40,0,2);
  }
  
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
  vector<float> *jetbtag_mm,*jetbtag_ee,*jetbtag_mc, *jetbtag_mct;
  Tree_DoubleMu->SetBranchAddress("jetBtag", &jetbtag_mm);
  Tree_DoubleElectron->SetBranchAddress("jetBtag", &jetbtag_ee);
  Tree_DYLLJets->SetBranchAddress("jetBtag", &jetbtag_mc);
  Tree_TTJets->SetBranchAddress("jetBtag", &jetbtag_mct);
  
  float mll_DoubleMu=0., mll_DoubleElectron=0., mll_DYLLJets=0, mll_TTJets;
  Tree_DoubleMu->SetBranchAddress("mll", &mll_DoubleMu);
  Tree_DoubleElectron->SetBranchAddress("mll", &mll_DoubleElectron);
  Tree_DYLLJets->SetBranchAddress("mll", &mll_DYLLJets);
  Tree_TTJets->SetBranchAddress("mll", &mll_TTJets);
  
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

  vector<float> *jetnhf_mm,*jetnhf_ee,*jetnhf_mc, *jetnhf_mct;
  Tree_DoubleMu->SetBranchAddress("jetNhf", &jetnhf_mm);
  Tree_DoubleElectron->SetBranchAddress("jetNhf", &jetnhf_ee);
  Tree_DYLLJets->SetBranchAddress("jetNhf", &jetnhf_mc);
  Tree_TTJets->SetBranchAddress("jetNhf", &jetnhf_mct);

  vector<float> *jetphf_mm,*jetphf_ee,*jetphf_mc, *jetphf_mct;
  Tree_DoubleMu->SetBranchAddress("jetPhf", &jetphf_mm);
  Tree_DoubleElectron->SetBranchAddress("jetPhf", &jetphf_ee);
  Tree_DYLLJets->SetBranchAddress("jetPhf", &jetphf_mc);
  Tree_TTJets->SetBranchAddress("jetPhf", &jetphf_mct);

  vector<float> *jetmetphi_mm,*jetmetphi_ee,*jetmetphi_mc, *jetmetphi_mct;
  Tree_DoubleMu->SetBranchAddress("jetMetPhi", &jetmetphi_mm);
  Tree_DoubleElectron->SetBranchAddress("jetMetPhi", &jetmetphi_ee);
  Tree_DYLLJets->SetBranchAddress("jetMetPhi", &jetmetphi_mc);
  Tree_TTJets->SetBranchAddress("jetMetPhi", &jetmetphi_mct);

  float met_ee, met_mm, met_mc, met_mct;
  Tree_DoubleMu->SetBranchAddress("met", &met_mm);
  Tree_DoubleElectron->SetBranchAddress("met", &met_ee);
  Tree_DYLLJets->SetBranchAddress("met", &met_mc); 
  Tree_TTJets->SetBranchAddress("met", &met_mct); 

  float rho_ee, rho_mm, rho_mc, rho_mct;
  Tree_DoubleMu->SetBranchAddress("rho", &rho_mm);
  Tree_DoubleElectron->SetBranchAddress("rho", &rho_ee);
  Tree_DYLLJets->SetBranchAddress("rho", &rho_mc); 
  Tree_TTJets->SetBranchAddress("rho", &rho_mct);

  vector<float> *jetvtx3del_mm,*jetvtx3del_ee,*jetvtx3del_mc, *jetvtx3del_mct;
  Tree_DoubleMu->SetBranchAddress("jetVtx3deL", &jetvtx3del_mm);
  Tree_DoubleElectron->SetBranchAddress("jetVtx3deL", &jetvtx3del_ee);
  Tree_DYLLJets->SetBranchAddress("jetVtx3deL", &jetvtx3del_mc);
  Tree_TTJets->SetBranchAddress("jetVtx3deL", &jetvtx3del_mct);

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

  vector<float> *jetmass_mm, *jetmass_ee, *jetmass_mc, *jetmass_mct;
  Tree_DoubleMu->SetBranchAddress("jetMass", &jetmass_mm);
  Tree_DoubleElectron->SetBranchAddress("jetMass", &jetmass_ee);
  Tree_DYLLJets->SetBranchAddress("jetMass", &jetmass_mc);
  Tree_TTJets->SetBranchAddress("jetMass", &jetmass_mct);
  
  vector<float> *jetvtxmass_mm, *jetvtxmass_ee, *jetvtxmass_mc, *jetvtxmass_mct;
  Tree_DoubleMu->SetBranchAddress("jetVtxMass", &jetvtxmass_mm);
  Tree_DoubleElectron->SetBranchAddress("jetVtxMass", &jetvtxmass_ee);
  Tree_DYLLJets->SetBranchAddress("jetVtxMass", &jetvtxmass_mc);
  Tree_TTJets->SetBranchAddress("jetVtxMass", &jetvtxmass_mct);

  vector<float> *jetsoftleppt_mm, *jetsoftleppt_ee, *jetsoftleppt_mc, *jetsoftleppt_mct;
  Tree_DoubleMu->SetBranchAddress("jetSoftLepPt", &jetsoftleppt_mm);
  Tree_DoubleElectron->SetBranchAddress("jetSoftLepPt", &jetsoftleppt_ee);
  Tree_DYLLJets->SetBranchAddress("jetSoftLepPt", &jetsoftleppt_mc);
  Tree_TTJets->SetBranchAddress("jetSoftLepPt", &jetsoftleppt_mct);

  vector<float> *jetsoftlepptrel_mm, *jetsoftlepptrel_ee, *jetsoftlepptrel_mc, *jetsoftlepptrel_mct;
  Tree_DoubleMu->SetBranchAddress("jetSoftLepPtRel", &jetsoftlepptrel_mm);
  Tree_DoubleElectron->SetBranchAddress("jetSoftLepPtRel", &jetsoftlepptrel_ee);
  Tree_DYLLJets->SetBranchAddress("jetSoftLepPtRel", &jetsoftlepptrel_mc);
  Tree_TTJets->SetBranchAddress("jetSoftLepPtRel", &jetsoftlepptrel_mct);

  vector<float> *jetleadtrkpt_mm, *jetleadtrkpt_ee, *jetleadtrkpt_mc, *jetleadtrkpt_mct;
  Tree_DoubleMu->SetBranchAddress("jetLeadTrkPt", &jetleadtrkpt_mm);
  Tree_DoubleElectron->SetBranchAddress("jetLeadTrkPt", &jetleadtrkpt_ee);
  Tree_DYLLJets->SetBranchAddress("jetLeadTrkPt", &jetleadtrkpt_mc);
  Tree_TTJets->SetBranchAddress("jetLeadTrkPt", &jetleadtrkpt_mct);

  vector<float> *jetpart_mm, *jetpart_ee, *jetpart_mc, *jetpart_mct;
  Tree_DoubleMu->SetBranchAddress("jetPart", &jetpart_mm);
  Tree_DoubleElectron->SetBranchAddress("jetPart", &jetpart_ee);
  Tree_DYLLJets->SetBranchAddress("jetPart", &jetpart_mc);
  Tree_TTJets->SetBranchAddress("jetPart", &jetpart_mct);
  
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
  bool dyismu=false, dyise=false, ttise=false, ttismu=false, dataismu=false, dataise=false;

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
      
      float deltaphi_mm = Z_mm.DeltaPhi(bjet0_mm);
      bool minv_mm = fabs(Z_mm.M()-91)<20;
      bool phi_mm = fabs(deltaphi_mm)>2.5;
      bool btag_mm = (*jetBtagM_mm)[0] && (*jetpuidm_mm)[0];
      TLorentzVector bjet1_mm;
      bool softjet_mm = 0;
      double pTveto_mm=0;
    
      //if there are >= 2 jets, then put a veto on the additional jet
      if(nJets_mm>1){                 
	bjet1_mm.SetPtEtaPhiE((*jetPt_mm)[1], (*jetEta_mm)[1],  (*jetPhi_mm)[1], (*jetE_mm)[1]);
	if((0.2*Z_mm.Pt())>10){
	  pTveto_mm = 0.2*Z_mm.Pt();
	}
	else{
	  pTveto_mm = 10;
	}
	softjet_mm =  bjet1_mm.Pt()>pTveto_mm && ((*jetpuidm_mm)[1]);
      }
      //set counters
      count1_mm++;
      if(minv_mm==true){
	countmll_mm++;
	if(btag_mm==true && bjet0_mm.Pt()>10){
	  countbtag_mm++;
	  if(phi_mm==true){
	    countphi_mm++;
	    if(softjet_mm==false){
	      countsoftjet_mm++;
	    }
	  }
	}
      }

      //if all the conditions are satisfied, fill the histograms
      if(minv_mm==1 && phi_mm==1 && btag_mm==1 && softjet_mm==0 && bjet0_mm.Pt()>10){
	countmm_tot++;
	hmll_DoubleMu->Fill(mll_DoubleMu);
	hjetbtag_mm->Fill((*jetbtag_mm)[0]);
	hjetpt_mm->Fill(bjet0_mm.Pt());
	hjeteta_mm->Fill(bjet0_mm.Eta());
	hjetnhf_mm->Fill((*jetnhf_mm)[0]);
	hjetphf_mm->Fill((*jetphf_mm)[0]);	
	hjetmetphi_mm->Fill((*jetmetphi_mm)[0]);
	hmet_mm->Fill(met_mm);
	hrho_mm->Fill(rho_mm);
	hjetvtx3del_mm->Fill((*jetvtx3del_mm)[0]);
	hjetmass_mm->Fill((*jetmass_mm)[0]);
	hjetvtxmass_mm->Fill((*jetvtxmass_mm)[0]);
	hjetsoftleppt_mm->Fill((*jetsoftleppt_mm)[0]);
	hjetsoftlepptrel_mm->Fill((*jetsoftlepptrel_mm)[0]);
	hjetleadtrkpt_mm->Fill((*jetleadtrkpt_mm)[0]);
	hjetpart_mm->Fill((*jetpart_mm)[0]);
	
	
	for(int i=0;i<11;i++){     //for each pt(Z) interval
	  if(Z_mm.Pt()>ptRatiobins[i] && Z_mm.Pt()<ptRatiobins[i+1]){ //if pt(Z) is in the selected interval
	    hptRatio_mm[i]->Fill(bjet0_mm.Pt()/Z_mm.Pt()); //fill unregressed pTratio
	    hptRegRatio_mm[i]->Fill((bjet0_mm.Pt()*(*jetregpt_mm)[0])/Z_mm.Pt());  //fill regressed pTratio
	  }
	}
      }  
    }
  } 

  //DIELECTRON CYCLE
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
     
      float deltaphi_ee = Z_ee.DeltaPhi(bjet0_ee);

      bool minv_ee = fabs(Z_ee.M()-91)<20;
      bool phi_ee = fabs(deltaphi_ee)>2.5;
      bool btag_ee = (*jetBtagM_ee)[0] && (*jetpuidm_ee)[0]; 
      TLorentzVector bjet1_ee;
      bool softjet_ee = 0;
      double pTveto_ee=0;
	
      //if there are >=2 jets, set the veto on the additional jet
      if(nJets_ee>1){                   
	bjet1_ee.SetPtEtaPhiE((*jetPt_ee)[1], (*jetEta_ee)[1],  (*jetPhi_ee)[1], (*jetE_ee)[1]);
	if((0.2*Z_ee.Pt())>10){
	  pTveto_ee = 0.2*Z_ee.Pt();
	}
	else{
	  pTveto_ee = 10;
	}
	softjet_ee =  bjet1_ee.Pt()>pTveto_ee && ((*jetpuidm_ee)[1]);
      }	

      //activate counters
      count1_ee++;
      if(minv_ee==true){
	countmll_ee++;
	if(btag_ee==true && bjet0_ee.Pt()>10){
	  countbtag_ee++;
	  if(phi_ee==true){
	    countphi_ee++;
	    if(softjet_ee==false){
	      countsoftjet_ee++;
	    }
	  }
	}
      }

      //if all the conditions are satisfied, fill histograms
      if(minv_ee==1 && btag_ee==1 && phi_ee==1 && softjet_ee==0 && bjet0_ee.Pt()>10){
	countee_tot++;
	hmll_DoubleElectron->Fill(mll_DoubleElectron);
	hjetbtag_ee->Fill((*jetbtag_ee)[0]);
	hjetpt_ee->Fill(bjet0_ee.Pt());
	hjeteta_ee->Fill(bjet0_ee.Eta());
	hjetnhf_ee->Fill((*jetnhf_ee)[0]);
	hjetphf_ee->Fill((*jetphf_ee)[0]);	
	hjetmetphi_ee->Fill((*jetmetphi_ee)[0]);
	hmet_ee->Fill(met_ee);
	hrho_ee->Fill(rho_ee);	
	hjetvtx3del_ee->Fill((*jetvtx3del_ee)[0]);
	hjetmass_ee->Fill((*jetmass_ee)[0]);
	hjetvtxmass_ee->Fill((*jetvtxmass_ee)[0]);
	hjetsoftleppt_ee->Fill((*jetsoftleppt_ee)[0]);
	hjetsoftlepptrel_ee->Fill((*jetsoftlepptrel_ee)[0]);
	hjetleadtrkpt_ee->Fill((*jetleadtrkpt_ee)[0]);
	hjetpart_ee->Fill((*jetpart_ee)[0]);
	
	for(int i=0;i<11;i++){ //for each pt(Z) interval
	  if(Z_ee.Pt()>ptRatiobins[i] && Z_ee.Pt()<ptRatiobins[i+1]){ //if pt(Z) is in the pt range
	    hptRatio_ee[i]->Fill(bjet0_ee.Pt()/Z_ee.Pt());  //fill unregressed pTratio histogram
	    hptRegRatio_ee[i]->Fill((bjet0_ee.Pt()*(*jetregpt_ee)[0])/Z_ee.Pt());  //fill regressed pTratio histogram
	  }
	}
      } 
    }
  }  
  
  //MONTE-CARLO DRELL-YAN CYCLE
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
    
    float deltaphi_mc = Z_mc.DeltaPhi(bjet0_mc);
        
    // correct MC for pu weigth
    float coeff = dataLumi/mcDYLumi*puWt[1];

    bool minv_mc = fabs(Z_mc.M()-91)<20;
    bool phi_mc = fabs(deltaphi_mc)>2.5;    
    bool btag_mc = (*jetBtagM_mc)[0] && (*jetpuidm_mc)[0];
    TLorentzVector bjet1_mc;
    bool softjet_mc = 0;
    double pTveto_mc =0;
    
    //if there are >=2 jets, set the veto on the additional jet
    if(nJets_mc>1){              
      bjet1_mc.SetPtEtaPhiE((*jetPt_mc)[1], (*jetEta_mc)[1],  (*jetPhi_mc)[1], (*jetE_mc)[1]);
      if((0.2*Z_mc.Pt())>10){
	pTveto_mc = 0.2*Z_mc.Pt();
      }
      else{
	pTveto_mc = 10;
      }
      softjet_mc =  bjet1_mc.Pt()>pTveto_mc && ((*jetpuidm_mc)[1]);
    }

    //activate counter for dielectron
    if(dyise){
      count1_mcee++;
      if(minv_mc==true){
	countmll_mcee++;
	if(btag_mc==true && bjet0_mc.Pt()>10){
	  countbtag_mcee++;
	  if(phi_mc==true){
	    countphi_mcee++;
	    if(softjet_mc==false){
	      countsoftjet_mcee++;
	    }
	  }
	}
      }
    }
    
    //activate counters for dimuon
    if(dyismu){
      count1_mcmm++;
      if(minv_mc==true){
	countmll_mcmm++;
	if(btag_mc==true && bjet0_mc.Pt()>10){
	  countbtag_mcmm++;
	  if(phi_mc==true){
	    countphi_mcmm++;
	    if(softjet_mc==false){
	      countsoftjet_mcmm++;
	    }
	  }
	}
      }
    }

    //if all the conditions are satisfied
    if(minv_mc==1 && phi_mc==1  && btag_mc==1 && softjet_mc==0 && bjet0_mc.Pt()>10){
      if(dyise){//dielectron trigger
	countmcee_tot++;
	hmll_DYLLJetsee->Fill(mll_DYLLJets,coeff);
	hjetbtag_mcee->Fill((*jetbtag_mc)[0],coeff);
      	hjetpt_mcee->Fill(bjet0_mc.Pt(), coeff);
	hjeteta_mcee->Fill(bjet0_mc.Eta(),coeff);
	hjetnhf_mcee->Fill((*jetnhf_mc)[0],coeff);
	hjetphf_mcee->Fill((*jetphf_mc)[0],coeff);	
	hjetmetphi_mcee->Fill((*jetmetphi_mc)[0],coeff);
	hmet_mcee->Fill(met_mc,coeff);
	hrho_mcee->Fill(rho_mc,coeff);	
	hjetvtx3del_mcee->Fill((*jetvtx3del_mc)[0],coeff);
	hjetmass_mcee->Fill((*jetmass_mc)[0],coeff);
	hjetvtxmass_mcee->Fill((*jetvtxmass_mc)[0],coeff);
	hjetsoftleppt_mcee->Fill((*jetsoftleppt_mc)[0],coeff);
	hjetsoftlepptrel_mcee->Fill((*jetsoftlepptrel_mc)[0],coeff);
	hjetleadtrkpt_mcee->Fill((*jetleadtrkpt_mc)[0],coeff);
	hjetpart_mcee->Fill((*jetpart_mc)[0],coeff);
	
	for(int i=0;i<11;i++){
	  if(Z_mc.Pt()>ptRatiobins[i] && Z_mc.Pt()<ptRatiobins[i+1]){
	    hptRatio_mcee[i]->Fill(bjet0_mc.Pt()/Z_mc.Pt(),coeff);
	    hptRegRatio_mcee[i]->Fill((bjet0_mc.Pt()*(*jetregpt_mc)[0])/Z_mc.Pt(),coeff);	    
	  }	
	}
      } 
      if(dyismu){ //dimuon trigger
	countmcmm_tot++;
	hmll_DYLLJetsmm->Fill(mll_DYLLJets, coeff);
	hjetbtag_mcmm->Fill((*jetbtag_mc)[0],coeff);
	hjetpt_mcmm->Fill(bjet0_mc.Pt(),coeff);
	hjeteta_mcmm->Fill(bjet0_mc.Eta(),coeff);
	hjetnhf_mcmm->Fill((*jetnhf_mc)[0],coeff);
	hjetphf_mcmm->Fill((*jetphf_mc)[0],coeff);	
	hjetmetphi_mcmm->Fill((*jetmetphi_mc)[0],coeff);
	hmet_mcmm->Fill(met_mc,coeff);
	hrho_mcmm->Fill(rho_mc,coeff);
	hjetvtx3del_mcmm->Fill((*jetvtx3del_mc)[0],coeff);
	hjetmass_mcmm->Fill((*jetmass_mc)[0],coeff);
	hjetvtxmass_mcmm->Fill((*jetvtxmass_mc)[0],coeff);
	hjetsoftleppt_mcmm->Fill((*jetsoftleppt_mc)[0],coeff);
	hjetsoftlepptrel_mcmm->Fill((*jetsoftlepptrel_mc)[0],coeff);
	hjetleadtrkpt_mcmm->Fill((*jetleadtrkpt_mc)[0],coeff);
	hjetpart_mcmm->Fill((*jetpart_mc)[0],coeff);

	for(int i=0;i<11;i++){
	  if(Z_mc.Pt()>ptRatiobins[i] && Z_mc.Pt()<ptRatiobins[i+1]){
	    hptRatio_mcmm[i]->Fill(bjet0_mc.Pt()/Z_mc.Pt(),coeff);
	    hptRegRatio_mcmm[i]->Fill((bjet0_mc.Pt()*(*jetregpt_mc)[0])/Z_mc.Pt(),coeff);
	  }	
	}
      }    
    }   
  }
  
 
  //MONTE-CARLo TTbar CYCLE 
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
    
    float deltaphi_mct = Z_mct.DeltaPhi(bjet0_mct);
        
    // correct MC for pu weigth
    float coefft = dataLumi/mcTTLumi*puWt[1];

    bool minv_mct = fabs(Z_mct.M()-91)<20;
    bool phi_mct = fabs(deltaphi_mct)>2.5;    
    bool btag_mct = (*jetBtagM_mct)[0] && (*jetpuidm_mct)[0];
    TLorentzVector bjet1_mct;
    bool softjet_mct = 0;
    double pTveto_mct=0;
    
    //if there are >=2 jets, set the veto on the additional jet
    if(nJets_mct>1){             
      bjet1_mct.SetPtEtaPhiE((*jetPt_mct)[1], (*jetEta_mct)[1],  (*jetPhi_mct)[1], (*jetE_mct)[1]);
      if((0.2*Z_mct.Pt())>10){
	pTveto_mct = 0.2*Z_mct.Pt();
      }
      else{
	pTveto_mct = 10;
      }
      softjet_mct =  bjet1_mct.Pt()>pTveto_mct && ((*jetpuidm_mct)[1]);      
    }

    //activate counter for dielectron
    if(ttise){
      count1_mctee++;
      if(minv_mct==true){
	countmll_mctee++;
	if(btag_mct==true && bjet0_mct.Pt()>10){
	  countbtag_mctee++;
	  if(phi_mct==true){
	    countphi_mctee++;
	    if(softjet_mct==false){
	      countsoftjet_mctee++;
	    }
	  }
	}
      }
    }

    //activate counter for dimuon
    if(ttismu){
      count1_mctmm++;
      if(minv_mct==true){
	countmll_mctmm++;
	if(btag_mct==true && bjet0_mct.Pt()>10){
	  countbtag_mctmm++;
	  if(phi_mct==true){
	    countphi_mctmm++;
	    if(softjet_mct==false){
	      countsoftjet_mctmm++;
	    }
	  }
	}
      }
    }
    
    //if all the conditins are satisfied
    if(minv_mct==1 && phi_mct==1  && btag_mct==1 && softjet_mct==0  && bjet0_mct.Pt()>10){
      if(ttise){  //dielectron trigger
	countmctee_tot++;
	hmll_TTJetsee->Fill(mll_TTJets,coefft);
	hjetbtag_mctee->Fill((*jetbtag_mct)[0],coefft);
      	hjetpt_mctee->Fill(bjet0_mct.Pt(), coefft);
	hjeteta_mctee->Fill(bjet0_mct.Eta(),coefft);
	hjetnhf_mctee->Fill((*jetnhf_mct)[0],coefft);
	hjetphf_mctee->Fill((*jetphf_mct)[0],coefft);	
	hjetmetphi_mctee->Fill((*jetmetphi_mct)[0],coefft);
	hmet_mctee->Fill(met_mct,coefft);
	hrho_mctee->Fill(rho_mct,coefft);	
	hjetvtx3del_mctee->Fill((*jetvtx3del_mct)[0],coefft);
	hjetmass_mctee->Fill((*jetmass_mct)[0],coefft);
	hjetvtxmass_mctee->Fill((*jetvtxmass_mct)[0],coefft);
	hjetsoftleppt_mctee->Fill((*jetsoftleppt_mct)[0],coefft);
	hjetsoftlepptrel_mctee->Fill((*jetsoftlepptrel_mct)[0],coefft);
	hjetleadtrkpt_mctee->Fill((*jetleadtrkpt_mct)[0],coefft);
	hjetpart_mctee->Fill((*jetpart_mct)[0],coefft);
	
	for(int i=0;i<11;i++){
	  if(Z_mct.Pt()>ptRatiobins[i] && Z_mct.Pt()<ptRatiobins[i+1]){
	    hptRatio_mctee[i]->Fill(bjet0_mct.Pt()/Z_mct.Pt(),coefft);
	    hptRegRatio_mctee[i]->Fill((bjet0_mct.Pt()*(*jetregpt_mct)[0])/Z_mct.Pt(),coefft);	    
	  }	
	}
      } 
      
      if(ttismu){ //dimuon trigger
	countmctmm_tot++;
	hmll_TTJetsmm->Fill(mll_TTJets, coefft);
	hjetbtag_mctmm->Fill((*jetbtag_mct)[0],coefft);
	hjetpt_mctmm->Fill(bjet0_mct.Pt(),coefft);
	hjeteta_mctmm->Fill(bjet0_mct.Eta(),coefft);
	hjetnhf_mctmm->Fill((*jetnhf_mct)[0],coefft);
	hjetphf_mctmm->Fill((*jetphf_mct)[0],coefft);	
	hjetmetphi_mctmm->Fill((*jetmetphi_mct)[0],coefft);
	hmet_mctmm->Fill(met_mct,coefft);
	hrho_mctmm->Fill(rho_mct,coefft);
	hjetvtx3del_mctmm->Fill((*jetvtx3del_mct)[0],coefft);
	hjetmass_mctmm->Fill((*jetmass_mct)[0],coefft);
	hjetvtxmass_mctmm->Fill((*jetvtxmass_mct)[0],coefft);
	hjetsoftleppt_mctmm->Fill((*jetsoftleppt_mct)[0],coefft);
	hjetsoftlepptrel_mctmm->Fill((*jetsoftlepptrel_mct)[0],coefft);
	hjetleadtrkpt_mctmm->Fill((*jetleadtrkpt_mct)[0],coefft);
	hjetpart_mctmm->Fill((*jetpart_mct)[0],coefft);

	for(int i=0;i<11;i++){
	  if(Z_mct.Pt()>ptRatiobins[i] && Z_mct.Pt()<ptRatiobins[i+1]){
	    hptRatio_mctmm[i]->Fill(bjet0_mct.Pt()/Z_mct.Pt(),coefft);
	    hptRegRatio_mctmm[i]->Fill((bjet0_mct.Pt()*(*jetregpt_mct)[0])/Z_mct.Pt(),coefft);
	  }	
	}
      }    
    }   
  }

  //cout counter for dimuon and dielectron
  
  cout << "*************** MUONS ********************" << endl;
  cout << "Trig+1bjet >>>> \t" << count1_mm << "\t" << count1_mcmm << "\t" << count1_mctmm << endl; 
  cout << "Minv >>>> \t" << countmll_mm << "\t" << countmll_mcmm << "\t" << countmll_mctmm << endl; 
  cout << "Btag + pt>10 >>>> \t" << countbtag_mm << "\t" << countbtag_mcmm << "\t" << countbtag_mctmm << endl;
  cout << "DeltaPhi >>>> \t" << countphi_mm << "\t" << countphi_mcmm << "\t" << countphi_mctmm << endl;  
  cout << "Additional Jet >>>> \t" << countsoftjet_mm << "\t" << countsoftjet_mcmm << "\t" << countsoftjet_mctmm << endl; 
  cout << "Total >>>> \t" << countmm_tot << "\t" << countmcmm_tot << "\t" << countmctmm_tot << endl;


  cout << "*************** ELECTRONS ********************" << endl;
  cout << "Trig+1bjet >>>> \t" << count1_ee << "\t" << count1_mcee << "\t" << count1_mctee << endl; 
  cout << "Minv >>>> \t" << countmll_ee << "\t" << countmll_mcee << "\t" << countmll_mctee << endl; 
  cout << "Btag + pt>10 >>>> \t" << countbtag_ee << "\t" << countbtag_mcee << "\t" << countbtag_mctee << endl;
  cout << "DeltaPhi >>>> \t" << countphi_ee << "\t" << countphi_mcee << "\t" << countphi_mctee << endl;  
  cout << "Additional Jet >>>> \t" << countsoftjet_ee << "\t" << countsoftjet_mcee << "\t" << countsoftjet_mctee << endl; 
  cout << "Total >>>> \t" << countee_tot << "\t" << countmcee_tot << "\t" << countmctee_tot << endl;
 
  //scale MC DY histograms with respect to Data-TTbar

  hjetbtag_mcee->Scale((hjetbtag_ee->Integral()-hjetbtag_mctee->Integral())/hjetbtag_mcee->Integral());
  hjetbtag_mcmm->Scale((hjetbtag_mm->Integral()-hjetbtag_mctmm->Integral())/hjetbtag_mcmm->Integral());
 
  hjetpt_mcee->Scale((hjetpt_ee->Integral()-hjetpt_mctee->Integral())/hjetpt_mcee->Integral());
  hjetpt_mcmm->Scale((hjetpt_mm->Integral()-hjetpt_mctmm->Integral())/hjetpt_mcmm->Integral());
  
  hjeteta_mcee->Scale((hjeteta_ee->Integral()-hjeteta_mctee->Integral())/hjeteta_mcee->Integral());
  hjeteta_mcmm->Scale((hjeteta_mm->Integral()-hjeteta_mctmm->Integral())/hjeteta_mcmm->Integral());
  
  hjetnhf_mcee->Scale((hjetnhf_ee->Integral()-hjetnhf_mctee->Integral())/hjetnhf_mcee->Integral());
  hjetnhf_mcmm->Scale((hjetnhf_mm->Integral()-hjetnhf_mctmm->Integral())/hjetnhf_mcmm->Integral());
  
  hjetphf_mcee->Scale((hjetphf_ee->Integral()-hjetphf_mctee->Integral())/hjetphf_mcee->Integral());
  hjetphf_mcmm->Scale((hjetphf_mm->Integral()-hjetphf_mctmm->Integral())/hjetphf_mcmm->Integral());

  hjetmetphi_mcee->Scale((hjetmetphi_ee->Integral()-hjetmetphi_mctee->Integral())/hjetmetphi_mcee->Integral());
  hjetmetphi_mcmm->Scale((hjetmetphi_mm->Integral()-hjetmetphi_mctmm->Integral())/hjetmetphi_mcmm->Integral());
  
  hmet_mcee->Scale((hmet_ee->Integral()-hmet_mctee->Integral())/hmet_mcee->Integral());
  hmet_mcmm->Scale((hmet_mm->Integral()-hmet_mctmm->Integral())/hmet_mcmm->Integral());
  
  hrho_mcee->Scale((hrho_ee->Integral()-hrho_mctee->Integral())/hrho_mcee->Integral());
  hrho_mcmm->Scale((hrho_mm->Integral()-hrho_mctmm->Integral())/hrho_mcmm->Integral());

  hjetvtx3del_mcee->Scale((hjetvtx3del_ee->Integral()-hjetvtx3del_mctee->Integral())/hjetvtx3del_mcee->Integral());
  hjetvtx3del_mcmm->Scale((hjetvtx3del_mm->Integral()-hjetvtx3del_mctmm->Integral())/hjetvtx3del_mcmm->Integral());

  hjetmass_mcee->Scale((hjetmass_ee->Integral()-hjetmass_mctee->Integral())/hjetmass_mcee->Integral());
  hjetmass_mcmm->Scale((hjetmass_mm->Integral()-hjetmass_mctmm->Integral())/hjetmass_mcmm->Integral());

  hjetvtxmass_mcee->Scale((hjetvtxmass_ee->Integral()-hjetvtxmass_mctee->Integral())/hjetvtxmass_mcee->Integral());
  hjetvtxmass_mcmm->Scale((hjetvtxmass_mm->Integral()-hjetvtxmass_mctmm->Integral())/hjetvtxmass_mcmm->Integral());

  hjetsoftleppt_mcee->Scale((hjetsoftleppt_ee->Integral()-hjetsoftleppt_mctee->Integral())/hjetsoftleppt_mcee->Integral());
  hjetsoftleppt_mcmm->Scale((hjetsoftleppt_mm->Integral()-hjetsoftleppt_mctmm->Integral())/hjetsoftleppt_mcmm->Integral());

  hjetsoftlepptrel_mcee->Scale((hjetsoftlepptrel_ee->Integral()-hjetsoftlepptrel_mctee->Integral())/hjetsoftlepptrel_mcee->Integral());
  hjetsoftlepptrel_mcmm->Scale((hjetsoftlepptrel_mm->Integral()-hjetsoftlepptrel_mctmm->Integral())/hjetsoftlepptrel_mcmm->Integral());
  
  hjetleadtrkpt_mcee->Scale((hjetleadtrkpt_ee->Integral()-hjetleadtrkpt_mctee->Integral())/hjetleadtrkpt_mcee->Integral());
  hjetleadtrkpt_mcmm->Scale(hjetleadtrkpt_mm->Integral()/hjetleadtrkpt_mcmm->Integral());

  hjetpart_mcee->Scale((hjetpart_ee->Integral()-hjetpart_mctee->Integral())/hjetpart_mcee->Integral());
  hjetpart_mcmm->Scale((hjetpart_mm->Integral()-hjetpart_mctmm->Integral())/hjetpart_mcmm->Integral());

  for(int c=0; c<11; c++){
    hptRatio_mcmm[c]->Scale((hptRatio_mm[c]->Integral()-hptRatio_mctmm[c]->Integral())/hptRatio_mcmm[c]->Integral());
    hptRatio_mcee[c]->Scale((hptRatio_ee[c]->Integral()-hptRatio_mctee[c]->Integral())/hptRatio_mcee[c]->Integral());
    hptRegRatio_mcmm[c]->Scale((hptRegRatio_mm[c]->Integral()-hptRegRatio_mctmm[c]->Integral())/hptRegRatio_mcmm[c]->Integral());
    hptRegRatio_mcee[c]->Scale((hptRegRatio_ee[c]->Integral()-hptRegRatio_mctee[c]->Integral())/hptRegRatio_mcee[c]->Integral());
  }
 

  //get RMS and MEAN  of pTratio response in each pt(Z) interval
  
  float RMS_ree[11], RMS_rmm[11], RMS_rmcee[11], RMS_rmcmm[11],  RMS_rmctee[11], RMS_rmctmm[11], RMS_rree[11], RMS_rrmm[11], RMS_rrmcee[11], RMS_rrmcmm[11],  RMS_rrmctee[11], RMS_rrmctmm[11];
  float RMSE_ree[11], RMSE_rmm[11], RMSE_rmcee[11], RMSE_rmcmm[11],RMSE_rmctee[11], RMSE_rmctmm[11],, RMSE_rree[11], RMSE_rrmm[11], RMSE_rrmcee[11], RMSE_rrmcmm[11], RMSE_rrmctee[11], RMSE_rrmctmm[11];
  
  float MEAN_ree[11], MEAN_rmm[11], MEAN_rmcee[11], MEAN_rmcmm[11], MEAN_rmctee[11], MEAN_rmctmm[11], MEAN_rree[11], MEAN_rrmm[11], MEAN_rrmcee[11], MEAN_rrmcmm[11], MEAN_rrmctee[11], MEAN_rrmctmm[11];
  float MEANE_ree[11], MEANE_rmm[11], MEANE_rmcee[11], MEANE_rmcmm[11], MEANE_rmctee[11], MEANE_rmctmm[11], MEANE_rree[11], MEANE_rrmm[11], MEANE_rrmcee[11], MEANE_rrmcmm[11],  MEANE_rrmctee[11], MEANE_rrmctmm[11];
  
  float RMS_MCee[11] ,  RMS_MCree[11], RMS_MCmm[11], RMS_MCrmm[11];
  float RMSE_MCee[11] ,   RMSE_MCree[11], RMSE_MCmm[11], RMSE_MCrmm[11];
  float MEAN_MCee[11] ,  MEAN_MCree[11], MEAN_MCmm[11], MEAN_MCrmm[11];
  float MEANE_MCee[11] ,   MEANE_MCree[11],MEANE_MCmm[11], MEANE_MCrmm[11];
  
  for(int d=0; d<11; d++){
    RMS_rmm[d]=hptRatio_mm[d]->GetRMS(1);
    RMS_ree[d]=hptRatio_ee[d]->GetRMS(1);
    RMS_rmcmm[d]=hptRatio_mcmm[d]->GetRMS(1);
    RMS_rmcee[d]=hptRatio_mcee[d]->GetRMS(1);
    RMS_rrmm[d]=hptRegRatio_mm[d]->GetRMS(1);
    RMS_rree[d]=hptRegRatio_ee[d]->GetRMS(1);
    RMS_rrmcmm[d]=hptRegRatio_mcmm[d]->GetRMS(1);
    RMS_rrmcee[d]=hptRegRatio_mcee[d]->GetRMS(1);
    RMSE_rmm[d]=hptRatio_mm[d]->GetRMSError(1);
    RMSE_ree[d]=hptRatio_ee[d]->GetRMSError(1);
    RMSE_rmcmm[d]=hptRatio_mcmm[d]->GetRMSError(1);
    RMSE_rmcee[d]=hptRatio_mcee[d]->GetRMSError(1);
    RMSE_rrmm[d]=hptRegRatio_mm[d]->GetRMSError(1);
    RMSE_rree[d]=hptRegRatio_ee[d]->GetRMSError(1);
    RMSE_rrmcmm[d]=hptRegRatio_mcmm[d]->GetRMSError(1);
    RMSE_rrmcee[d]=hptRegRatio_mcee[d]->GetRMSError(1);
    
    MEAN_rmm[d]=hptRatio_mm[d]->GetMean(1); 
    MEAN_ree[d]=hptRatio_ee[d]->GetMean(1);
    MEAN_rmcmm[d]=hptRatio_mcmm[d]->GetMean(1);
    MEAN_rmcee[d]=hptRatio_mcee[d]->GetMean(1);
    MEAN_rrmm[d]=hptRegRatio_mm[d]->GetMean(1); 
    MEAN_rree[d]=hptRegRatio_ee[d]->GetMean(1);
    MEAN_rrmcmm[d]=hptRegRatio_mcmm[d]->GetMean(1);
    MEAN_rrmcee[d]=hptRegRatio_mcee[d]->GetMean(1);
    MEANE_rmm[d]=hptRatio_mm[d]->GetMeanError(1); 
    MEANE_ree[d]=hptRatio_ee[d]->GetMeanError(1);
    MEANE_rmcmm[d]=hptRatio_mcmm[d]->GetMeanError(1);
    MEANE_rmcee[d]=hptRatio_mcee[d]->GetMeanError(1);
    MEANE_rrmm[d]=hptRegRatio_mm[d]->GetMeanError(1); 
    MEANE_rree[d]=hptRegRatio_ee[d]->GetMeanError(1);
    MEANE_rrmcmm[d]=hptRegRatio_mcmm[d]->GetMeanError(1);
    MEANE_rrmcee[d]=hptRegRatio_mcee[d]->GetMeanError(1);
    
    //in order to plot the RMS(mean) vs pt(Z) it is important to add the to MC components

    /*  hptRatio_mctee[d]->Add(hptRatio_mcee[d],(hptRatio_ee[d]->Integral()-hptRatio_mctee[d]->Integral()));
	RMS_MCee[d]=   hptRatio_mctee[d]->GetRMS(1);
	RMSE_MCee[d]=   hptRatio_mctee[d]->GetRMSError(1);
	MEAN_MCee[d]=   hptRatio_mctee[d]->GetMean(1);
	MEANE_MCee[d]=   hptRatio_mctee[d]->GetMeanError(1); 
	hptRegRatio_mctee[d]->Add(hptRegRatio_mcee[d],(hptRegRatio_ee[d]->Integral()-hptRegRatio_mctee[d]->Integral()));
	RMS_MCree[d]=   hptRegRatio_mctee[d]->GetRMS(1);
	RMSE_MCree[d]=   hptRegRatio_mctee[d]->GetRMSError(1);
	MEAN_MCree[d]=   hptRegRatio_mctee[d]->GetMean(1);
	MEANE_MCree[d]=   hptRegRatio_mctee[d]->GetMeanError(1);
	hptRatio_mctmm[d]->Add(hptRatio_mcmm[d],(hptRatio_mm[d]->Integral()-hptRatio_mctmm[d]->Integral()));
	RMS_MCmm[d]=   hptRatio_mctmm[d]->GetRMS(1);
	RMSE_MCmm[d]=   hptRatio_mctmm[d]->GetRMSError(1);
	MEAN_MCmm[d]=   hptRatio_mctmm[d]->GetMean(1);
	MEANE_MCmm[d]=   hptRatio_mctmm[d]->GetMeanError(1);
	hptRegRatio_mctmm[d]->Add(hptRegRatio_mcmm[d],(hptRegRatio_mm[d]->Integral()-hptRegRatio_mctmm[d]->Integral()));
	RMS_MCrmm[d]=   hptRegRatio_mctmm[d]->GetRMS(1);
	RMSE_MCrmm[d]=   hptRegRatio_mctmm[d]->GetRMSError(1);
	MEAN_MCrmm[d]=   hptRegRatio_mctmm[d]->GetMean(1);
	MEANE_MCrmm[d]=   hptRegRatio_mctmm[d]->GetMeanError(1);*/
    
  }
  
  gStyle->SetOptStat(0);
  
  //plot input variables canvas

  /* histo_inputvariables("jetpt_ee", hjetpt_ee, hjetpt_mcee, hjetpt_mctee, "p_{T} [GeV]", false, true);
     histo_inputvariables("jetpt_mm",hjetpt_mm, hjetpt_mcmm, hjetpt_mctmm, "p_{T} [GeV]", true, true);
     histo_inputvariables("jeteta_ee", hjeteta_ee, hjeteta_mcee, hjeteta_mctee, "#eta", false, false);
     histo_inputvariables("jeteta_mm",hjeteta_mm, hjeteta_mcmm, hjeteta_mctmm, "#eta", true, false);
     histo_inputvariables("jetnhf_ee",hjetnhf_ee, hjetnhf_mcee, hjetnhf_mctee,"Neutral Hadron Fraction", false, true);
     histo_inputvariables("jetnhf_mm",hjetnhf_mm, hjetnhf_mcmm, hjetnhf_mctmm,"Neutral Hadron Fraction", true, true);
     histo_inputvariables("jetphf_ee",hjetphf_ee, hjetphf_mcee, hjetphf_mctee,"Photon Fraction", false, false);
     histo_inputvariables("jetphf_mm",hjetphf_mm, hjetphf_mcmm, hjetphf_mctmm,"Photon Fraction", true, false);
     histo_inputvariables("jetmetphi_ee",hjetmetphi_ee, hjetmetphi_mcee, hjetmetphi_mctee, "#Delta#phi(MET,jet)", false, false);
     histo_inputvariables("jetmetphi_mm", hjetmetphi_mm, hjetmetphi_mcmm, hjetmetphi_mctmm,"#Delta#phi(MET,jet)", true, false);
     histo_inputvariables("met_ee",hmet_ee, hmet_mcee, hmet_mctee,"MET [GeV]", false, false);
     histo_inputvariables("met_mm",hmet_mm, hmet_mcmm, hmet_mctmm,"MET [GeV]", true, false);
     histo_inputvariables("rho_ee",hrho_ee, hrho_mcee, hrho_mctee,"#rho", false, false);
     histo_inputvariables("rho_mm",hrho_mm, hrho_mcmm, hrho_mctmm, "#rho", true, false);
     histo_inputvariables("jetvtx3del_ee", hjetvtx3del_ee, hjetvtx3del_mcee, hjetvtx3del_mctee, "Vtx 3deL", false,true);
     histo_inputvariables("jetvtx3del_mm",hjetvtx3del_mm, hjetvtx3del_mcmm, hjetvtx3del_mctmm,"Vtx 3deL", true, true);
     histo_inputvariables("jetmass_ee", hjetmass_ee, hjetmass_mcee, hjetmass_mctee,"Mass [GeV]", false,true);
     histo_inputvariables("jetmass_mm", hjetmass_mm, hjetmass_mcmm, hjetmass_mctmm,"Mass [GeV]", true,true);
     histo_inputvariables("jetvtxmass_ee", hjetvtxmass_ee, hjetvtxmass_mcee, hjetvtxmass_mctee,"Vtx Mass [GeV]", false,true);
     histo_inputvariables("jetvtxmass_mm", hjetvtxmass_mm, hjetvtxmass_mcmm, hjetvtxmass_mctmm, "Vtx Mass [GeV]", true,true);
     histo_inputvariables("jetsoftleppt_ee", hjetsoftleppt_ee, hjetsoftleppt_mcee, hjetsoftleppt_mctee,"SL p_{T} [GeV]", false, true);
     histo_inputvariables("jetsoftleppt_mm", hjetsoftleppt_mm, hjetsoftleppt_mcmm, hjetsoftleppt_mctmm, "SL p_{T} [GeV]", true,true);
     histo_inputvariables("jetsoftlepptrel_ee", hjetsoftlepptrel_ee, hjetsoftlepptrel_mcee, hjetsoftlepptrel_mctee, "SL p_{T,rel} [GeV]", false,true);
     histo_inputvariables("jetsoftlepptrel_mm", hjetsoftlepptrel_mm, hjetsoftlepptrel_mcmm, hjetsoftlepptrel_mctmm,"SL p_{T,rel} [GeV]", true,true);
     histo_inputvariables("jetleadtrkpt_ee", hjetleadtrkpt_ee, hjetleadtrkpt_mcee, hjetleadtrkpt_mctee,"p_{T} leading track [GeV]", false,false);
     histo_inputvariables("jetleadtrkpt_mm", hjetleadtrkpt_mm, hjetleadtrkpt_mcmm, hjetleadtrkpt_mctmm,"p_{T} leading track [GeV]", true,false);
     histo_inputvariables("jetpart_ee", hjetpart_ee, hjetpart_mcee, hjetpart_mctee,"N_{const}", false,false);
     histo_inputvariables("jetpart_mm", hjetpart_mm, hjetpart_mcmm, hjetpart_mctmm, "N_{const}", true,false);*/
  
  
  //plot ptratio canvas

  char name_rmm[100], name_ree[100], name_rrmm[100], name_rree[100];
  char name_rmmfit[100], name_reefit[100], name_rrmmfit[100], name_rreefit[100];
  
  for(int b=0; b<11; b++){
    //ptratio canvas

    sprintf(name_rmm, "ptRatio_mm%d", b);
    //histo_ptRatio(name_rmm, hptRatio_mm[b], hptRatio_mcmm[b], hptRatio_mctmm[b]);
    sprintf(name_ree, "ptRatio_ee%d", b);
    // histo_ptRatio(name_ree,  hptRatio_ee[b], hptRatio_mcee[b], hptRatio_mctee[b]);
    sprintf(name_rrmm, "ptRatio_rrmm%d", b);
    // histo_ptRatio(name_rrmm, hptRegRatio_mm[b], hptRegRatio_mcmm[b], hptRegRatio_mctmm[b]);
    sprintf(name_rree, "ptRatio_rree%d", b);
    //histo_ptRatio(name_rree, hptRegRatio_ee[b], hptRegRatio_mcee[b], hptRegRatio_mctee[b]);
    
    //fitted ptratio canvas

    sprintf(name_rmmfit, "ptRatio_mmfit%d", b);
    // histo_ptRatioFIT(name_rmmfit, hptRatio_mm[b], hptRatio_mcmm[b], hptRatio_mctmm[b]);
    sprintf(name_reefit, "ptRatio_eefit%d", b);
    // histo_ptRatioFIT(name_reefit,  hptRatio_ee[b], hptRatio_mcee[b], hptRatio_mctee[b]);
    sprintf(name_rrmmfit, "ptRatio_rrmmfit%d", b);
    // histo_ptRatioFIT(name_rrmmfit, hptRegRatio_mm[b], hptRegRatio_mcmm[b], hptRegRatio_mctmm[b]);
    sprintf(name_rreefit, "ptRatio_rreefit%d", b);
    //histo_ptRatioFIT(name_rreefit, hptRegRatio_ee[b], hptRegRatio_mcee[b], hptRegRatio_mctee[b]);
    
  }
  
  //plot RMS and MEAN canvas
  
  /* RMS("RMSmm", "rms_mm", "rms_mcmm",RMS_rmm, RMS_MCmm, RMSE_rmm, RMSE_MCmm, ptRatiobins, true,false, true) ;
     RMS("RMSee", "rms_ee", "rms_mcee",RMS_ree, RMS_MCee, RMSE_ree, RMSE_MCee, ptRatiobins, true,false, false);
     RMS("RMSrmm", "rms_rmm", "rms_rmcmm",RMS_rrmm, RMS_MCrmm, RMSE_rrmm, RMSE_MCrmm, ptRatiobins, true,false, true);
     RMS("RMSree", "rms_ree", "rms_rmcee",RMS_rree, RMS_MCree, RMSE_rree, RMSE_MCree, ptRatiobins, true,false , false);
     RMS("MEANmm", "mean_mm", "mean_mcmm", MEAN_rmm, MEAN_MCmm, MEANE_rmm, MEANE_MCmm, ptRatiobins, false,false,true);
     RMS("MEANee", "mean_ee", "mean_mcee", MEAN_ree, MEAN_MCee, MEANE_ree, MEANE_MCee, ptRatiobins, false,false,false);
     RMS("MEANrmm", "mean_rmm", "mean_rmcmm", MEAN_rrmm, MEAN_MCrmm, MEANE_rrmm, MEANE_MCrmm, ptRatiobins, false,false,true);
     RMS("MEANree", "mean_ree", "mean_rmcee", MEAN_rree, MEAN_MCree, MEANE_rree, MEANE_MCree, ptRatiobins, false,false,false);*/
  
  
  //plot mu e sigma fit canvas
  
  /****FIT CON CONSTRAINT GAUSSIAN******/

  //initialize sigma and mean values, obtained from a 3-iteration fitting procedure

  //Unregressed

  /*float sigma_ee[11]={0.29736,0.274222,0.256071,0.216321,0.218728,0.19779,0.171276,0.168842,0.153795,0.125787,0.0882091};
  float sigma_mcee[11]={0.272615,0.260555,0.239276,0.212059,0.182509,0.17925,0.16528,0.151015,0.14336,0.128601,0.0970573};
  float mu_ee[11]={0.884901,0.87535,0.88293,0.908215,0.90865,0.907541,0.956426,0.927406,0.932846,0.92904,0.901273};
  float mu_mcee[11]={0.882241,0.883599,0.93457,0.949847,0.927692,0.941246,0.934652,0.947158,0.933988,0.938529,0.997161};
  float sigmaE_ee[11]={0.00797647,0.00926866,0.00936242,0.0109899,0.0125523,0.0108234,0.0109574,0.0105579,0.0154998,0.0116246,0.0153667};
  float sigmaE_mcee[11]={0.00777331,0.00936694,0.0144346,0.0123008,0.0119157,0.0115876,0.0112599,0.00880707,0.00790019,0.0128607,0.012286};
  float muE_ee[11]={0.00666491, 0.00772431,0.00824543,0.00886438,0.0102443,0.00989532,0.00846675,0.00886028,0.0115753,0.00797888,0.0122471};
  float muE_mcee[11]={0.00632917,0.00746413,0.00936321,0.00922784,0.00861763,0.00959471,0.00804013,0.00803751,0.00879619,0.00763909,0.013236};

  float sigma_mm[11] = {0.293007, 0.295533,0.254358,0.224627,0.191867,0.18497,0.157033,0.15018,0.185007,0.128487,0.145502};
  float sigma_mcmm[11] = {0.298319,0.261788,0.239426,0.213531,0.215616, 0.177791,0.159256,0.163094,0.162755,0.129154,0.141777};
  float mu_mm[11]={0.88288,0.886656,0.904412,0.92714,0.931066,0.921377,0.93759,0.926672,0.914714,0.935972,0.941325} ;
  float mu_mcmm[11]= {0.903324,0.895319,0.897387,0.917813, 0.922454,0.936056,0.965166,0.940281,0.961347,0.942799,-0.961571 };
  float sigmaE_mm[11] = {0.00882296,0.011846, 0.00960378,0.012197,0.0102434,0.0126199,0.0103156,0.0109096,0.0170314, 0.00615303,0.0142586 };
  float sigmaE_mcmm[11]={0.00691078,0.00902151,0.0106514, 0.0103264,0.0112357, 0.0103453,0.0105199,0.0109532,0.00126835,0.0104876,0.026171};
  float muE_mm[11]={0.00647094,0.00855227,0.00787848,0.00889729, 0.00889336,0.00953539,0.00719224,0.00831924,0.0126732,0.00611884,0.0160031	};
  float muE_mcmm[11]={0.00615469, 0.00698799,0.00796392,0.00817547,0.00928446,0.00892856,0.0073572,0.00882684,0.0114695,0.0078852,0.0199524};
  
  FIT("SIGMA_ee", "sigma_ee", "sigma_mcee", sigma_ee, sigma_mcee, sigmaE_ee, sigmaE_mcee, ptRatiobins, true, false, false);
  FIT("MU_ee", "mu_ee", "mu_mcee", mu_ee, mu_mcee, muE_ee, muE_mcee, ptRatiobins, false, false, false);
  FIT("SIGMA_mm", "sigma_mm", "sigma_mcmm", sigma_mm, sigma_mcmm, sigmaE_mm, sigmaE_mcmm, ptRatiobins, true, false, true);
  FIT("MU_mm", "mu_mm", "mu_mcmm", mu_mm, mu_mcmm, muE_mm, muE_mcmm, ptRatiobins, false, false, true);

  //Regressed 

  float sigma_ree[11] = {0.156183, 0.233329,0.243974,0.21602 ,0.167794,0.17407, 0.167525,0.145103,0.136287,0.0954597,0.100845};
  float sigma_rmcee[11] = {0.152238,0.205075,0.245213,0.192292,0.164877, 0.163367,0.157227,0.122102,0.138148 , 0.127753,0.0991002};
  float mu_ree[11]={1.09437,0.929335 ,0.926245,0.950408,0.955053,0.9513,0.968895, 0.951784 ,0.948884 ,0.945479,0.943855} ;
  float mu_rmcee[11]= {1.09271,0.969991,0.968061,0.980063 ,0.955774,0.971021,0.957778, 0.968649,0.931667,0.947795 ,0.977895};
  float sigmaE_ree[11] = {0.00356422,0.0138309,0.0114703,0.0127471,0.0100034,0.0124613,0.0100916,0.0127218,0.0130619,0.00685046,0.0208502};
  float sigmaE_rmcee[11]={ 0.0033977,0.00648185,0.0118072,0.0063654,0.0116, 0.0082735,0.0109592,0.0109577,0.00825011, 0.00842912,0.0190861};
  float muE_ree[11]={0.00356545,0.00868194,0.00861663,0.00927936, 0.00808937,0.00954522,0.00788624, 0.00890743,0.0101587, 0.0055768,0.014377 };
  float muE_rmcee[11]={0.0034686,0.00579472, 0.00875276,0.00692643, 0.00778102,0.00788153,0.00768663,0.00818119,0.00895392,0.00684163,0.0129105};
  
  float sigma_rmm[11] = {0.16679,0.300328,0.245721,0.212596,0.177496,0.17175,0.149731 ,0.146295,0.146016,0.125597,0.130412};
  float sigma_rmcmm[11] = {0.157132,0.183311,0.224912,0.224227,0.157315,0.161165 ,0.127247,0.107455,0.141682,0.121981,0.12762};
  float mu_rmm[11]={1.08804,0.889683,0.947964,0.970115,0.970115 ,0.964401,0.966129,0.945581,0.942239,0.950776,0.952592} ;
  float mu_rmcmm[11]= {1.09623,0.959229,0.944575,0.977389,0.96154,0.983017,0.985509,0.96909,0.957543,0.945045,0.947456};
  float sigmaE_rmm[11] = {0.00481662,0.0318698,0.0107335,0.0100521,0.00980127,0.0121806,0.00785038,0.00777864,0.0118868,0.008301,0.0146816};
  float sigmaE_rmcmm[11]={0.00399414,0.00649037,0.00853642,0.00855095,0.0106261,0.00967026,0.00496237, 0.00499754,0.0127637,0.00745672,0.0150479};
  float muE_rmm[11]={0.00367115,0.0248354,0.00804365,0.00789985,0.00784971,0.00917975,0.00634604, 0.00721388,0.0103019,0.00666075,0.0154179};
  float muE_rmcmm[11]={0.00327412,0.00507957, 0.00695091,0.00752621,0.0079169,0.00831693,0.00483483,0.00597001,0.0104036,0.00630685,0.0151346};
  
  FIT("MU_rmm", "mu_rmm", "mu_rmcmm", mu_rmm, mu_rmcmm, muE_rmm, muE_rmcmm, ptRatiobins, false, false, true);
  FIT("SIGMA_rmm", "sigma_rmm", "sigma_rmcmm", sigma_rmm, sigma_rmcmm, sigmaE_rmm, sigmaE_rmcmm, ptRatiobins, true, false, false);
  FIT("MU_ree", "mu_ree", "mu_rmcee", mu_ree, mu_rmcee, muE_ree, muE_rmcee, ptRatiobins, false, false, true);
  FIT("SIGMA_ree", "sigma_ree", "sigma_rmcee", sigma_ree, sigma_rmcee, sigmaE_ree, sigmaE_rmcee, ptRatiobins, true, false, false);*/
  
  

}

//function to plot input variables histograms

void histo_inputvariables(TString HISTO, TH1F *h_data, TH1F *h_mc, TH1F *h_mct, TString yAxis, bool ismu, bool logscale){
  TCanvas *can = new TCanvas("can"+HISTO, "can"+HISTO, 900, 600);
  TPad *pad1 = new TPad("pad1"+HISTO, "pad1"+HISTO, 0., 0.3, 1.,1.);
  pad1->SetBottomMargin(0.);
  pad1->Draw();
  pad1->cd();
  h_mct->SetTitle("");
  h_mct->SetLineColor(kBlack);
  h_mct->SetFillColor(kGreen-2);
  h_mc->SetFillColor(kRed-10);
  h_mc->SetLineColor(kBlack);
  THStack *hs = new THStack("hs","");
  hs->Add(h_mct);
  hs->Add(h_mc);
  hs->Draw("");
  can->Update();
  hs->SetMinimum(0.5);
  hs->GetYaxis()->SetTitle("Events");
  hs->GetYaxis()->SetLabelSize(0.04);
  hs->GetYaxis()->SetTitleSize(0.06);
  hs->GetYaxis()->SetTitleOffset(0.6);
  h_data->SetMarkerStyle(20);
  h_data->SetMarkerSize(0.8);
  h_data->DrawCopy("same ep");
  leg = new TLegend(0.91,0.60,0.99,0.80); 
  if(ismu) leg->SetHeader("Z(#mu#mu)+jets");
  if(!ismu) leg->SetHeader("Z(ee)+jets");
  leg->SetTextSize(0.03);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->AddEntry(h_data,"Data","lep");
  leg->AddEntry(h_mc,"DY MC","f");
  leg->AddEntry(h_mct, "t#bar{t} MC", "f");
  leg->Draw();
  if(logscale) gPad->SetLogy();
  gPad->SetTickx();
  gPad->SetTicky();
  can->cd();
  TPad *pad2 = new TPad("pad2","pad2",0,0.,1,0.3);
  pad2->SetBottomMargin(0.4);
  pad2->SetTopMargin(0.);
  pad2->SetGridy();
  pad2->Draw();
  pad2->cd();
  h_data->Sumw2();
  h_mc->Sumw2();
  h_mct->Sumw2();
  h_data->SetStats(0);
  h_mc->Add(h_mct);
  h_data->Divide((h_mc));
  h_data->SetMarkerStyle(20);
  h_data->GetYaxis()->SetTitle("Data/MC");
  h_data->GetXaxis()->SetTitle(""+yAxis);
  h_data->GetYaxis()->SetLabelSize(0.1);
  h_data->GetYaxis()->SetTitleSize(0.1);
  h_data->GetXaxis()->SetLabelSize(0.1);
  h_data->GetXaxis()->SetTitleSize(0.1);
  h_data->GetYaxis()->SetTitleOffset(0.4);
  h_data->GetXaxis()->SetTitleOffset(1);
  h_data->GetYaxis()->SetNdivisions(505,-1);
  h_data->SetMinimum(0.35);
  h_data->SetMaximum(1.65);
  h_data->Draw("ep");
  h_data->SetTitle("");
  gPad->SetTickx();
  gPad->SetTicky(); 
  can->cd();
  can->Print(HISTO+".pdf");
}


//function to plot pTratio histograms

void histo_ptRatio(TString HISTO, TH1F *h_data, TH1F *h_mc, TH1F *h_mct){
  TCanvas *canvas = new TCanvas("canvas_"+HISTO, "canvas_"+HISTO, 900, 900);
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
  leg = new TLegend(0.63,0.65,0.78,0.85); 
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

//function to plot pTratio fitted histograms (3-interation Gaussian)

void histo_ptRatioFIT(TString HISTO, TH1F *h_data, TH1F *h_mc, TH1F *h_mct){
  TCanvas *canvas = new TCanvas("canvas_"+HISTO, "canvas_"+HISTO, 900, 900);
  TF1* Gaussian_data = new TF1("Gaussian_data","gaus",0.,2.);
  Gaussian_data->SetLineStyle(1);
  Gaussian_data->SetLineWidth(3);
  Gaussian_data->SetLineColor(1);
  TF1* Gaussian_mc = new TF1("Gaussian_mc","gaus",0.,2.);
  Gaussian_mc->SetLineStyle(10);
  Gaussian_mc->SetLineWidth(3);
  Gaussian_mc->SetLineColor(2);
  h_data->Fit("Gaussian_data", "", "", 0.7, 1.3);
  h_mc->Add(h_mct);
  h_mc->Fit("Gaussian_mc", "", "", 0.7, 1.3);
  double sigma_data = Gaussian_data->GetParameter(2);
  double sigmaE_data = Gaussian_data->GetParError(2);
  double sigma_mc = Gaussian_mc->GetParameter(2);
  double sigmaE_mc = Gaussian_mc->GetParError(2);
  double  mu_data = Gaussian_data->GetParameter(1);
  double muE_data = Gaussian_data->GetParError(1);
  double mu_mc = Gaussian_mc->GetParameter(1);
  double muE_mc = Gaussian_mc->GetParError(1);
  h_data->Fit("Gaussian_data", "", "",mu_data - 1.5*sigma_data, mu_data + 1.5*sigma_data);
  h_mc->Fit("Gaussian_mc", "", "", mu_mc - 1.5*sigma_mc, mu_mc + 1.5*sigma_mc);
  sigma_data = Gaussian_data->GetParameter(2);
  sigmaE_data = Gaussian_data->GetParError(2);
  sigma_mc = Gaussian_mc->GetParameter(2);
  sigmaE_mc = Gaussian_mc->GetParError(2);
  mu_data = Gaussian_data->GetParameter(1);
  muE_data = Gaussian_data->GetParError(1);
  mu_mc = Gaussian_mc->GetParameter(1);
  muE_mc = Gaussian_mc->GetParError(1);
  h_data->Fit("Gaussian_data", "", "",mu_data - 1.5*sigma_data, mu_data + 1.5*sigma_data);
  h_mc->Fit("Gaussian_mc", "", "", mu_mc - 1.5*sigma_mc, mu_mc + 1.5*sigma_mc);
  sigma_data = Gaussian_data->GetParameter(2);
  sigmaE_data = Gaussian_data->GetParError(2);
  sigma_mc = Gaussian_mc->GetParameter(2);
  sigmaE_mc = Gaussian_mc->GetParError(2);
  mu_data = Gaussian_data->GetParameter(1);
  muE_data = Gaussian_data->GetParError(1);
  mu_mc = Gaussian_mc->GetParameter(1);
  muE_mc = Gaussian_mc->GetParError(1);
  cout << "***********************************ECCOCI******************************************"<< endl;
  cout << sigma_data << " +/- " << sigmaE_data << "\t -*-*-"<< sigma_mc << " +/- "  << sigmaE_mc << endl;
  cout << mu_data  << " +/- " << muE_data << "\t -*-*-"<< mu_mc << " +/- " << muE_mc << endl;
  cout << "*****************************************************************************"<< endl;
  h_mc->GetYaxis()->SetTitle("Events");
  h_mc->GetXaxis()->SetTitle("Response");
  h_mc->SetTitle("");
  h_mc->SetLineColor(kBlack);
  h_mc->SetFillColor(kRed-10);
  h_mc->Draw("histo");
  h_data->SetMarkerSize(1.2);
  h_data->SetMarkerStyle(20);
  h_data->Draw("same ep");
  Gaussian_data->Draw("same");
  Gaussian_mc->Draw("same");
  leg = new TLegend(0.63,0.65,0.78,0.85); 
  leg->SetTextSize(0.03);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->AddEntry(h_data,"Data","lep");
  leg->AddEntry(Gaussian_data, "Fit to Data", "l");
  leg->AddEntry(h_mc,"DY MC","f");
  leg->AddEntry(Gaussian_mc, "Fit to MC", "l");
  leg->Draw();
  gPad->SetTickx();
  gPad->SetTicky();
  canvas->cd();
  canvas->Print(HISTO+".pdf");
  }

//function to plot RMS and Mean distribution with respect to the pt(Z)

void RMS(TString HISTO, TString HISTO_DATA, TString HISTO_MC, float *v_data, float *v_mc, float *ve_data, float *ve_mc, float *bins,bool isrms, bool islog, bool ismu) {
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
  if(isrms==true) h_data->GetYaxis()->SetTitle("Regressed RMS(p_{T}^{b}/p_{T}^{Z})");
  if(isrms==false) h_data->GetYaxis()->SetTitle("Regressed Mean(p_{T}^{b}/p_{T}^{Z})");
  h_data->SetTitle("");
  if (isrms==true) h_data->SetMinimum(0.11);
  if (isrms==false) h_data->SetMinimum(0.87);
  if (isrms==false) h_data->SetMaximum(1.17);
  h_data->GetYaxis()->SetLabelSize(0.04);
  h_data->GetYaxis()->SetTitleSize(0.06);
  h_data->GetYaxis()->SetTitleOffset(0.7);
  h_data->SetMarkerStyle(20);
  h_data->SetMarkerColor(kBlue);
  h_data->SetMarkerSize(0.8);
  h_data->DrawCopy("ep");
  h_mc->SetMarkerStyle(21);
  h_mc->SetMarkerColor(kMagenta);
  h_mc->SetMarkerSize(0.8);
  h_mc->DrawCopy("same ep");
  leg = new TLegend(0.73,0.65,0.88,0.85); 
  leg->SetTextSize(0.05);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->AddEntry(h_data,"Data","lep");
  leg->AddEntry(h_mc,"MC","lep");
  leg->Draw();
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
  if(ismu==true){
    h_data->SetMinimum(0.75);
    h_data->SetMaximum(1.25);
  }
  if(ismu==false){
    h_data->SetMinimum(0.75);
    h_data->SetMaximum(1.25);
  }
  h_data->Draw("ep");
  h_data->SetTitle("");
  gPad->SetLogx();
  gPad->SetTickx();
  gPad->SetTicky();
  canvas->cd();
  canvas->Print(HISTO+".pdf");
}

//function to plot sigma e mu of the ptRatio distribution from the 3-interation fit procedure

void FIT(TString HISTO, TString HISTO_DATA, TString HISTO_MC, float *v_data, float *v_mc,float *ve_data, float *ve_mc,  float *bins,bool isrms, bool islog, bool ismu) {
  TH1F *h_data = new TH1F(HISTO_DATA,HISTO_DATA,11, bins);
  TH1F *h_mc = new TH1F(HISTO_MC,HISTO_MC,11, bins);
  for(int i=0; i<11; i++){
    h_data->SetBinContent(i+1,v_data[i]);
    h_mc->SetBinContent(i+1,v_mc[i]);
    h_data->SetBinError(i+1, ve_data[i]);
    h_mc->SetBinError(i+1, ve_mc[i]);
  }  
   TCanvas *canvas = new TCanvas("canvas_"+HISTO, "canvas_"+HISTO, 900, 700);
  TPad *pad1 = new TPad("pad1"+HISTO, "pad1"+HISTO, 0., 0.3, 1.,1.);
  pad1->SetBottomMargin(0.);
  pad1->Draw();
  pad1->cd();
  if(isrms==true) h_data->GetYaxis()->SetTitle("Unregressed #sigma(p_{T}^{b}/p_{T}^{Z})");
  if(isrms==false) h_data->GetYaxis()->SetTitle("Unregressed #mu(p_{T}^{b}/p_{T}^{Z})");
  h_data->SetTitle("");
  h_data->GetYaxis()->SetLabelSize(0.04);
  h_data->GetYaxis()->SetTitleSize(0.06);
  h_data->GetYaxis()->SetTitleOffset(0.7);
  h_data->SetMarkerStyle(20);
  h_data->SetMarkerColor(kBlue);
  h_data->SetMarkerSize(0.8);
  h_data->DrawCopy("ep");
  h_mc->SetMarkerStyle(21);
  h_mc->SetMarkerColor(kMagenta);
  h_mc->SetMarkerSize(0.8);
  h_mc->DrawCopy("same ep");
  leg = new TLegend(0.73,0.65,0.88,0.85); 
  leg->SetTextSize(0.05);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->AddEntry(h_data,"Data","lep");
  leg->AddEntry(h_mc,"MC","lep");
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
  h_data->SetMinimum(0.77);
  h_data->SetMaximum(1.23);
  h_data->Draw("ep");
  h_data->SetTitle("");
  gPad->SetTickx();
  gPad->SetTicky();
  canvas->cd();
  canvas->Print(HISTO+".pdf");
}

