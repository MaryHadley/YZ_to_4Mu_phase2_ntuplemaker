// h o w   t o   r u n
// root -l
// .L loop.C++
// run("mc_ZUpsi.root") 

#include <vector>
#include <iostream>
#include <string>
using namespace std;
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TLorentzVector.h"

#include <math.h>
#include <TMath.h>

TFile *root_file;
#include "tree.C"
#include "treeMC.C"

tree *TREE;
treeMC *TREEMC;

//test function
//You can add functions in this space outside of the definition of run!
int myAdd(int a, int b){
  int c;
  c = a + b;
  return c;

}

void run(string file){//, string file2){

  // l o a d   t h e   t r e e s 
  root_file = new TFile(file.c_str(),"READ");
  TREE   = new tree((TTree*)root_file->Get("tree"));
  TREEMC = new treeMC((TTree*)root_file->Get("treemc"));
  
  //Announce what root_file is, thanks to Maxi for showing me how to do this
  std::cout << "////////////////////////////////////////" << std::endl;
  std::cout << "Processing file:  " << file.c_str() << std::endl;
  std::cout << "////////////////////////////////////////" << std::endl;
  

  // h i s t o g r a m s
  TH1F *h_reco_Z_mass_noNewCuts    = new TH1F("h_reco_Z_mass_noNewCuts",    "h_reco_Z_mass_noNewCuts", 20, 66., 116.);  h_reco_Z_mass_noNewCuts   ->SetXTitle("m_{#mu#mu} [GeV]"); //might want to change the binning, this is currently 20 bins to cover a range of 5
  h_reco_Z_mass_noNewCuts->Sumw2();
  
  TH1F *h_reco_Upsi_mass_noNewCuts = new TH1F("h_reco_Upsi_mass_noNewCuts", "h_reco_Upsi_mass_noNewCuts", 20, 8., 12.); h_reco_Upsi_mass_noNewCuts->SetXTitle("m_{#mu#mu} [GeV]"); //20 bins here to cover a range of 4 
  h_reco_Upsi_mass_noNewCuts->Sumw2();
  
  TH1F *h_big4MuVtxProb_before_big4MuVtx_Prob_Cut = new TH1F("h_big4MuVtxProb_before_big4MuVtx_Prob_Cut", "h_big4MuVtxProb_before_big4MuVtx_Prob_Cut",200, 0, 1); h_big4MuVtxProb_before_big4MuVtx_Prob_Cut->SetXTitle("big4MuVtxProb_before_big4MuVtx_Prob_Cut");
  h_big4MuVtxProb_before_big4MuVtx_Prob_Cut->Sumw2();
  
//  TH1F *h_dimuon_from_Z_Prob_before_Cut = new TH1F("h_dimuon_from_Z_Prob_before_Cut", "h_dimuon_from_Z_Prob_before_Cut", 200, 0, 1); h_dimuon_from_Z_Prob_before_Cut->SetXTitle("h_dimuon_from_Z_Prob_before_Cut");
  
//  TH1F *h_dimuon_from_upsi_before_Cut  = new TH1F("h_dimuon_from_upsi_before_Cut ",  "h_dimuon_from_upsi_before_Cut ", 200, 0, 1);  h_dimuon_from_upsi_before_Cut ->SetXTitle("h_dimuon_from_upsi_before_Cut ");
  
  TH1F *h_ambig_quad = new TH1F("h_ambi_quad",    "h_ambi_quad", 5, -0.5, 4.5);  h_ambig_quad  ->SetXTitle("Sum of pair_12_34_56, pair_12_34_56, pair_13_24_56, pair_14_23_56");
  h_ambig_quad->Sumw2();
  
  TH1F *h_cutflow_allQuadCuts = new TH1F("h_cutflow_allQuadCuts", "h_cutflow_allQuadCounts", 5, -0.5, 4.5); h_cutflow_allQuadCuts->SetXTitle("Cuts involving overall quad");
  h_cutflow_allQuadCuts->Sumw2();
  
  TH1F *h_cutflow_Z_first_upsi_phase1_second_pair_12_34_56 = new TH1F("h_cutflow_Z_first_upsi_phase1_second_pair_12_34_56","h_cutflow_Z_first_upsi_phase1_second_pair_12_34_56", 15, -0.5, 14.5); h_cutflow_Z_first_upsi_phase1_second_pair_12_34_56->SetXTitle("Cutflow for the Z_first_upsi_phase1_second_pair_12_34_56 case");
  h_cutflow_Z_first_upsi_phase1_second_pair_12_34_56->Sumw2();
  
  TH1F *h_pfIso_lep1 = new TH1F("h_pfIso_lep1", "h_pfIso_lep1", 60, 0, 3); h_pfIso_lep1->SetXTitle("PF Isolation for lep1");
  h_pfIso_lep1->Sumw2();
  
  TH1F *h_pfIso_lep2 = new TH1F("h_pfIso_lep2", "h_pfIso_lep2", 60, 0, 3); h_pfIso_lep2->SetXTitle("PF Isolation for lep2");
  h_pfIso_lep2->Sumw2();
  
  TH1F *h_pfIso_lep3 = new TH1F("h_pfIso_lep3", "h_pfIso_lep3", 60, 0, 3); h_pfIso_lep3->SetXTitle("PF Isolation for lep3");
  h_pfIso_lep3->Sumw2();
  
  TH1F *h_pfIso_lep4 = new TH1F("h_pfIso_lep4", "h_pfIso_lep4", 60, 0, 3); h_pfIso_lep4->SetXTitle("PF Isolation for lep4");
  h_pfIso_lep4 ->Sumw2();
  
  //Ignoring MC for the moment 
  TH1F *h_truth_Z_mass    = new TH1F("h_truth_Z_mass",    "h_truth_Z_mass", 20, 66., 116.);  h_truth_Z_mass->SetMarkerSize(0); //If I change the binning above, would also want to change it here so the truth and recovered plots have same scale 
  h_truth_Z_mass->Sumw2();
  
  TH1F *h_truth_Upsi_mass = new TH1F("h_truth_Upsi_mass", "h_truth_Upsi_mass", 20, 8., 12.); h_truth_Upsi_mass->SetMarkerSize(0); //same comment as above for the upsi truth and recovered mass plots 
  h_truth_Upsi_mass->Sumw2(); 

  // v a r i a b l e s
  double muon_mass = 105.6583 / 1000.; //get mass in GeV
  
//  std::cout << "myAdd(3,5) "<< myAdd(3,5) << std::endl;
  
  //boolean flags
  
  bool doMCTruthMatching = false; //code working for !doMCTruthMatching and doMCTruthMatching :)
  
  
  if (!doMCTruthMatching){
     std::cout << "NOT performing MC truth matching" << std::endl;
     std::cout << "////////////////////////////////" << std::endl;

  }
  
  if (doMCTruthMatching){
     std::cout << "Performing MC truth matching! Make sure you are running on MC!" << std::endl; 
     std::cout << "/////////////////////////////////////////////////////////////" << std::endl; 
  }
  
  //counters
  int pair_12_34_56_count = 0;
  
  int pair_13_24_56_count = 0;
  
  int pair_14_23_56_count = 0; 
  
  int pair_AMBIGUOUS_muQuad_count = 0;
  
  int big4MuVtx_Prob_Cut_fail_count = 0;
  
  int Z_first_upsi_phase1_second_pair_12_34_56_count = 0;
  
  int upsi_phase1_first_Z_second_pair_12_34_56_count = 0;
  
  int Z_first_upsi_phase1_second_pair_13_24_56_count = 0;
  
  int upsi_phase1_first_Z_second_pair_13_24_56_count = 0;
  
  int Z_first_upsi_phase1_second_pair_14_23_56_count = 0;
  
  int upsi_phase1_first_Z_second_pair_14_23_56_count = 0;
  
  int GotHereCount_Z_first_upsi_phase1_second_pair_12_34_56 = 0;
  
  int pfIso_Fail_Count = 0;
  
  int FailureCount = 0;
  
  int QuickCheckCount = 0;
  
  int fillCount =0; 
  
  int gotToEndCount = 0;
  
  //int poodleCount2
  int poodleCount = 0; 
  
  //counters that only have meaning in the doMCTruthMatching case
  int mcSanityCheckCount = 0;
  
  int matchedZCount = 0;
  
  int matchedCount = 0; //matchedCount indicates that we have matched all 4 muons appropriately, 2 to the Z, 2 to an UpsiN, and we have determined what that N is
  
  
  
  
  //Cuts 
  double big4MuVtx_Prob_Cut = 0.01; 
  
  double Z_mass_low  = 66.;
  
  double Z_mass_high = 116.;
  
  double upsi_mass_low_phase1 = 8;
  
  double upsi_mass_high_phase1 = 12;
  
  double upsi_mass_low_phase2 = 8.5;
  
  double upsi_mass_high_phase2 = 11; 
  
  double lead_mu_from_Z_pT_Cut = 30;
  
  double sublead_mu_from_Z_pT_Cut = 15;
  
  double mu_from_upsi_pT_Cut = 4; //We could lower this and gain back a lot of muons, is it worth it in terms of signal to background trade off? TO DISCUSS
  
  double mu_from_Z_eta_Cut = 2.4;
  
  double mu_from_upsi_eta_Cut = 2.4;
  
  double mu_from_Z_3DIPSig_Cut = 4;
  
 // double mu_from_upsi_RAPIDITY_Cut = 2.4; // keep commented in for the moment so things compile
  
  double upsi_RAPIDITY_Cut = 2.4; //we are cutting on the rapidity of the upsi, not on the rapidity of its daughter muons
  
  
  //Trying to put back at least the mu_mu_from_upsi_Prob_Cut to see if it helps us see the upsi 1. In the general case, we don't store the info to associate a given dimuon to a given quad correctly, 
  // but we use it only after we have thrown out all the quads that have ambiguous pairings, which is what leads to there being more dimuon(1,2) vertices than overall quads,
  //and so we have ensured that we have the same number of dimuon(1,2) vertices as we do quads, aka that the length of dimuon(1,2) vectors is equal to the length of the quad quantity vectors
  //so we can use the ->at(i) and know we are getting the dimuon(1,2) that matches that quad without ambiguity 
  
 // double mu_mu_from_Z_Prob_Cut = 0.05; 
  
  double mu_mu_from_upsi_Prob_Cut = 0.05; // 0.1; //suggested by S.L. in the command referenced in the google document "Investigations_into_recovering_upsi_signal_in_data"
  
  double mu_mu_from_Z_Prob_Cut = 0.05; // 0.1;
  
  double pfIso_Cut_Mu_from_Z = 0.35;
  
  double pfIso_Cut_Mu_from_Upsi = 1.65;
  
  double deltaRCut = 0.01;
  
 
  // n e w  s k i m m e d   r o o t   f i l e
  double mass1_quickAndDirty = -99;
  double mass2_quickAndDirty  = -99;
  double Z_mass = -99;
  double upsi_mass = -99;
  double Z_pT = -99;
  double Z_eta = -99;
  double  Z_RAPIDITY = -99;
  double Z_phi = -99;
  double upsi_pT = -99;
  double upsi_eta = -99;
  double upsi_RAPIDITY = -99; 
  double upsi_phi = -99;
  double lead_pT_mu_from_Z_pT = -99;
  double lead_pT_mu_from_Z_eta = -99;
  double lead_pT_mu_from_Z_RAPIDITY = -99;
  double lead_pT_mu_from_Z_phi = -99;
  double sublead_pT_mu_from_Z_pT = -99;
  double sublead_pT_mu_from_Z_eta = -99;
  double sublead_pT_mu_from_Z_RAPIDITY = -99;
  double sublead_pT_mu_from_Z_phi = -99;
  double lead_pT_mu_from_upsi_pT = -99;
  double lead_pT_mu_from_upsi_eta = -99;
  double lead_pT_mu_from_upsi_RAPIDITY = -99;
  double lead_pT_mu_from_upsi_phi = -99;
  double sublead_pT_mu_from_upsi_pT = -99;
  double sublead_pT_mu_from_upsi_eta = -99;
  double sublead_pT_mu_from_upsi_RAPIDITY = -99; 
  double sublead_pT_mu_from_upsi_phi = -99;
  
  //this variable will only be meaningful when doing MC Truth Matching
  int upsi_type = -1; //upsi_type = 1 corresponds to the upsi(1S), upsi_type = 2 corresponds to the upsi(2S), upsi_type = 3 corresponds to the upsi(3S)
  
 
   
  
  TFile *ntuple = new TFile("ntuple_skimmed_inputFileIs_Run2016_Total_16Dec2021_pfIso0p35_forZmu_and1p65_forUpsimu.root", "RECREATE");
  TTree *aux;
  aux = new TTree("tree", "tree");
  aux->Branch("mass1_quickAndDirty", &mass1_quickAndDirty);
  aux->Branch("mass2_quickAndDirty", &mass2_quickAndDirty);
  
  //NOTE: variables are recovered aka what we get from applying all cuts and looking at the reco values (as opposed to MC truth) unless otherwise noted 
  
  //Z and upsi masses
  aux->Branch("Z_mass", &Z_mass);
  aux->Branch("upsi_mass", &upsi_mass);
  
  //Z and upsi pT, eta,phi, RAPIDITY
  aux->Branch("Z_pT", &Z_pT);
  aux->Branch("Z_eta", &Z_eta);
  aux->Branch("Z_RAPIDITY", &Z_RAPIDITY);
  aux->Branch("Z_phi", &Z_phi);
  aux->Branch("upsi_pT", &upsi_pT);
  aux->Branch("upsi_eta", &upsi_eta);
  aux->Branch("upsi_RAPIDITY", &upsi_RAPIDITY);
  aux->Branch("upsi_phi", &upsi_phi);
  
  //muons from Z branches
  aux->Branch("lead_pT_mu_from_Z_pT", &lead_pT_mu_from_Z_pT);
  aux->Branch("lead_pT_mu_from_Z_eta", &lead_pT_mu_from_Z_eta);
  aux->Branch("lead_pT_mu_from_Z_RAPIDITY", &lead_pT_mu_from_Z_RAPIDITY);
  aux->Branch("lead_pT_mu_from_Z_phi", &lead_pT_mu_from_Z_phi);
  aux->Branch("sublead_pT_mu_from_Z_pT", &sublead_pT_mu_from_Z_pT);
  aux->Branch("sublead_pT_mu_from_Z_eta", &sublead_pT_mu_from_Z_eta);
  aux->Branch("sublead_pT_mu_from_Z_RAPIDITY", &sublead_pT_mu_from_Z_RAPIDITY);
  aux->Branch("sublead_pT_mu_from_Z_phi", &sublead_pT_mu_from_Z_phi);
  
  //muons from upsi branches, in the data scenario where we can't distinguish one species of upsi from another (i.e. in the non MC truth matching scenario)
  aux->Branch("lead_pT_mu_from_upsi_pT", &lead_pT_mu_from_upsi_pT);
  aux->Branch("lead_pT_mu_from_upsi_eta", &lead_pT_mu_from_upsi_eta);
  aux->Branch("lead_pT_mu_from_upsi_RAPIDITY", &lead_pT_mu_from_upsi_RAPIDITY);
  aux->Branch("lead_pT_mu_from_upsi_phi", &lead_pT_mu_from_upsi_phi);
  aux->Branch("sublead_pT_mu_from_upsi_pT", &sublead_pT_mu_from_upsi_pT);
  aux->Branch("sublead_pT_mu_from_upsi_eta", &sublead_pT_mu_from_upsi_eta);
  aux->Branch("sublead_pT_mu_from_upsi_RAPIDITY", &sublead_pT_mu_from_upsi_RAPIDITY); 
  aux->Branch("sublead_pT_mu_from_upsi_phi", &sublead_pT_mu_from_upsi_phi);
  
  aux->Branch("upsi_type", &upsi_type);

///////////////////////////
//////    D A T A    //////
///////////////////////////
  int eventCounter = 0;
  int entries = (TREE->fChain)->GetEntries(); //might want to change to GetEntriesFast by that can wait
  for(int iEntry=0; iEntry<entries; iEntry++) {
    (TREE->fChain)->GetEntry(iEntry);
     eventCounter += 1;
    
    double temp_comparison_pt_upsilon = 0;
    double temp_comparison_pt_z = 0;
    
    //I think this is where I want in my temp_<someVar> vectors!! CHECK ME
    std::vector<double> temp_Z_mass;
    std::vector<double> temp_upsi_mass;
    
    std::vector<double> temp_Z_pT;
    std::vector<double> temp_Z_eta;
    std::vector<double>temp_Z_RAPIDITY;
    std::vector<double>temp_Z_phi;
    std::vector<double> temp_upsi_pT;
    std::vector<double> temp_upsi_eta;
    std::vector<double> temp_upsi_RAPIDITY;
    std::vector<double> temp_upsi_phi;
    
    std::vector<double> temp_lead_pT_mu_from_Z_pT;
    std::vector<double> temp_lead_pT_mu_from_Z_eta;
    std::vector<double> temp_lead_pT_mu_from_Z_RAPIDITY;
    std::vector<double> temp_lead_pT_mu_from_Z_phi;
    std::vector<double> temp_sublead_pT_mu_from_Z_pT;
    std::vector<double> temp_sublead_pT_mu_from_Z_eta;
    std::vector<double> temp_sublead_pT_mu_from_Z_RAPIDITY;
    std::vector<double> temp_sublead_pT_mu_from_Z_phi;
    
    std::vector<double> temp_lead_pT_mu_from_upsi_pT;
    std::vector<double> temp_lead_pT_mu_from_upsi_eta;
    std::vector<double> temp_lead_pT_mu_from_upsi_RAPIDITY;
    std::vector<double> temp_lead_pT_mu_from_upsi_phi;
    std::vector<double> temp_sublead_pT_mu_from_upsi_pT;
    std::vector<double> temp_sublead_pT_mu_from_upsi_eta;
    std::vector<double> temp_sublead_pT_mu_from_upsi_RAPIDITY;
    std::vector<double> temp_sublead_pT_mu_from_upsi_phi;
    
    std::vector<int> temp_upsi_type;
    
    temp_Z_mass.clear();
    
    
//    std::cout <<"temp_Z_mass.size()" <<  temp_Z_mass.size() << std::endl; 
    temp_upsi_mass.clear();
    
    temp_Z_pT.clear();
    temp_Z_eta.clear();
    temp_Z_RAPIDITY.clear();
    temp_Z_phi.clear();
    
    temp_upsi_pT.clear();
    temp_upsi_eta.clear();
    temp_upsi_RAPIDITY.clear();
    temp_upsi_phi.clear();
    
    temp_lead_pT_mu_from_Z_pT.clear();
    temp_lead_pT_mu_from_Z_eta.clear();
    temp_lead_pT_mu_from_Z_RAPIDITY.clear();
    temp_lead_pT_mu_from_Z_phi.clear();
    
    temp_sublead_pT_mu_from_Z_pT.clear();
    temp_sublead_pT_mu_from_Z_eta.clear();
    temp_sublead_pT_mu_from_Z_RAPIDITY.clear();
    temp_sublead_pT_mu_from_Z_phi.clear();
    
    temp_lead_pT_mu_from_upsi_pT.clear();
    temp_lead_pT_mu_from_upsi_eta.clear();
    temp_lead_pT_mu_from_upsi_RAPIDITY.clear();
    temp_lead_pT_mu_from_upsi_phi.clear();
    
    temp_sublead_pT_mu_from_upsi_pT.clear();
    temp_sublead_pT_mu_from_upsi_eta.clear();
    temp_sublead_pT_mu_from_upsi_RAPIDITY.clear();
    temp_sublead_pT_mu_from_upsi_phi.clear();
    
    temp_upsi_type.clear();
    
    mass1_quickAndDirty = 0.; mass2_quickAndDirty = 0.;
    
 //   survivor_Z_first_upsi_phase1_second_pair_12_34_56 = false; //I think I don't need this 

    for (int i=0; i<(int)TREE->lepton1_pt->size(); i++) { //note to self: check that lepton1 size will always = lepton 2 size = lepton 3 size = lepton 4 size, I think so but need to double check  //this should be true by definition 
      // lepton1_pt > lepton2_pt > lepton3_pt > lepton4_pt, this is an artifact of how things were done in phase 1 code 
      TLorentzVector lepton1, lepton2, lepton3, lepton4;
      lepton1.SetPtEtaPhiM ( TREE->lepton1_pt->at(i), TREE->lepton1_eta->at(i), TREE->lepton1_phi->at(i), muon_mass);
      lepton2.SetPtEtaPhiM ( TREE->lepton2_pt->at(i), TREE->lepton2_eta->at(i), TREE->lepton2_phi->at(i), muon_mass);
      lepton3.SetPtEtaPhiM ( TREE->lepton3_pt->at(i), TREE->lepton3_eta->at(i), TREE->lepton3_phi->at(i), muon_mass);
      lepton4.SetPtEtaPhiM ( TREE->lepton4_pt->at(i), TREE->lepton4_eta->at(i), TREE->lepton4_phi->at(i), muon_mass);
      
      unsigned int runNumThisQuad = TREE->run_number->at(i);
      unsigned int evNumThisQuad  = TREE->event_number->at(i);
      unsigned int  LSThisQuad     = TREE->lumi_section->at(i);
      //std::cout << "runNumThisQuad  " << runNumThisQuad << std::endl;
      
      h_reco_Upsi_mass_noNewCuts->Fill( (lepton3+lepton4).M() );
      h_reco_Z_mass_noNewCuts->Fill( (lepton1+lepton2).M() );
      
      h_cutflow_allQuadCuts->Fill(1); //here are all the candidate quads
      //Cuts involving the overall quad //
      /////////////////////////////////////
      
      //deal with pairing ambiguous muon quads, eliminate those quads from our consideration 
   //    if ( (TREE->pair_12_34_56->at(i) == 1 && TREE->pair_13_24_56->at(i) == 1) || (TREE->pair_12_34_56->at(i) == 1 && TREE->pair_14_23_56->at(i) == 1)
//            || (TREE->pair_13_24_56->at(i) == 1 && TREE->pair_14_23_56->at(i) == 1)
//            || (TREE->pair_12_34_56->at(i) == 1 && TREE->pair_13_24_56->at(i) == 1 && TREE->pair_14_23_56->at(i) == 1) )
     
     int theSum;
     
     theSum = TREE->pair_12_34_56->at(i) + TREE->pair_13_24_56->at(i) + TREE->pair_14_23_56->at(i);
//     std::cout << "theSum: " << theSum << std::endl;
     h_ambig_quad->Fill(theSum);

     if ( TREE->pair_12_34_56->at(i) + TREE->pair_13_24_56->at(i) + TREE->pair_14_23_56->at(i) > 1) { //cleaner way suggested by S.L., equivalent to what I tried above but shorter!
             std::cout << "FOUND PAIRING AMBIGUOUS QUAD OF MUONS, WILL THROW IT AWAY" << std::endl;
             pair_AMBIGUOUS_muQuad_count += 1;
             continue;
        }
      h_cutflow_allQuadCuts->Fill(2); // here are the quads that survive the ambiguous pair cut
      
      h_big4MuVtxProb_before_big4MuVtx_Prob_Cut->Fill(TREE->big4MuVtx->at(i)); //fill it  before we cut on it
        
      if (TREE->big4MuVtx->at(i) < big4MuVtx_Prob_Cut){ //KEEP HIGH P VALUES! THINK P FOR PROBABILITY! CONFIRMED that TMath::Prob returns a p value, so low indicates stat significance. CONFIRMED WE WANT TO THROW AWAY LOW P VALUES, SEE NOTES IN JULY_2021_LAB_NOTEBOOK!
//         std::cout << "FAILED big4MuVtx_Prob_Cut! Throwing away this quad!" << std::endl; 
//         std::cout << TREE->big4MuVtx->at(i) << std::endl; 
         big4MuVtx_Prob_Cut_fail_count +=1;
         continue;
         }  
      
      h_cutflow_allQuadCuts->Fill(3); // here are the quads that survive the prob cut
      
      
 //     std:: cout << "Checking what TMath::Prob gives, let's try TMath::Prob(3.84, 1)   " << TMath::Prob(3.84, 1) << std::endl; //https://en.wikipedia.org/wiki/Chi-square_distribution //confirmed that this gives out what we think it should, aka this returns .05
//       std:: cout << "Checking what TMath::Prob gives, let's try TMath::Prob(3.32, 9)   " << TMath::Prob(3.32, 9) << std::endl;
      
      //Put in iso003 cuts here, that's the last cut involving the quad //this part taken from https://github.com/cms-ljmet/FWLJMET/blob/c319f38c1e34cf9f0277bd00231e6b75c889523b/LJMet/plugins/MultiLepEventSelector.cc#L542-L548 thank you Sinan for pointing me to this!
      
      //lepton1 
      
      double chIso_lep1 = TREE->lepton1_iso03hadron->at(i);
  //    std::cout << "chIso_lep1:  " << chIso_lep1 << std::endl;
 
      double nhIso_lep1 = TREE->lepton1_iso03neutralHadron->at(i);
 //     std::cout << "nhIso_lep1  " << nhIso_lep1 << std::endl; 

      double gIso_lep1 = TREE->lepton1_iso03photon->at(i); //g for gamma aka the photon
//      std::cout << "gIso_lep1  " << gIso_lep1 << std::endl;
      
      double puIso_lep1 = TREE->lepton1_iso03PU->at(i);
//      std::cout << "puIso_lep1  " << puIso_lep1 << std::endl;
      
      double pT_lep1 = TREE->lepton1_pt->at(i);
      
      double pfIso_lep1 = (chIso_lep1 + std::max(0.,nhIso_lep1 + gIso_lep1 - 0.5*puIso_lep1))/pT_lep1;
      
//      std::cout << "pfIso_lep1 " << pfIso_lep1 << std::endl; 
      h_pfIso_lep1->Fill(pfIso_lep1);
     
      //lepton2
      
      double chIso_lep2 = TREE->lepton2_iso03hadron->at(i);
//      std::cout << "chIso_lep2:  " << chIso_lep2 << std::endl;
 
      double nhIso_lep2 = TREE->lepton2_iso03neutralHadron->at(i);
//      std::cout << "nhIso_lep2  " << nhIso_lep2 << std::endl; 

      double gIso_lep2 = TREE->lepton2_iso03photon->at(i); //g for gamma aka the photon
 //     std::cout << "gIso_lep2  " << gIso_lep2 << std::endl;
      
      double puIso_lep2 = TREE->lepton2_iso03PU->at(i);
 //     std::cout << "puIso_lep2  " << puIso_lep2 << std::endl;
      
      double pT_lep2 = TREE->lepton2_pt->at(i);
      
      double pfIso_lep2 = (chIso_lep2 + std::max(0.,nhIso_lep2 + gIso_lep2 - 0.5*puIso_lep2))/pT_lep2;
      
      h_pfIso_lep2->Fill(pfIso_lep2);
      
      //lepton3
      
      double chIso_lep3 = TREE->lepton3_iso03hadron->at(i);
//      std::cout << "chIso_lep3:  " << chIso_lep3 << std::endl;
 
      double nhIso_lep3 = TREE->lepton3_iso03neutralHadron->at(i);
  //    std::cout << "nhIso_lep3  " << nhIso_lep3 << std::endl; 

      double gIso_lep3 = TREE->lepton3_iso03photon->at(i); //g for gamma aka the photon
 //     std::cout << "gIso_lep3  " << gIso_lep3 << std::endl;
      
      double puIso_lep3 = TREE->lepton3_iso03PU->at(i);
 //     std::cout << "puIso_lep3  " << puIso_lep3 << std::endl;
      
      double pT_lep3 = TREE->lepton3_pt->at(i);
      
      double pfIso_lep3 = (chIso_lep3 + std::max(0.,nhIso_lep3 + gIso_lep3 - 0.5*puIso_lep3))/pT_lep3;
      
      h_pfIso_lep3->Fill(pfIso_lep3);
      
      //lepton4
      double chIso_lep4 = TREE->lepton4_iso03hadron->at(i);
  //    std::cout << "chIso_lep4:  " << chIso_lep4 << std::endl;
 
      double nhIso_lep4 = TREE->lepton4_iso03neutralHadron->at(i);
 //     std::cout << "nhIso_lep4  " << nhIso_lep4 << std::endl; 

      double gIso_lep4 = TREE->lepton4_iso03photon->at(i); //g for gamma aka the photon
 //     std::cout << "gIso_lep4  " << gIso_lep4 << std::endl;
      
      double puIso_lep4 = TREE->lepton4_iso03PU->at(i);
 //     std::cout << "puIso_lep4  " << puIso_lep4 << std::endl;
      
      double pT_lep4 = TREE->lepton4_pt->at(i);
      
      double pfIso_lep4 = (chIso_lep4 + std::max(0.,nhIso_lep4 + gIso_lep4 - 0.5*puIso_lep4))/pT_lep4;
      
      h_pfIso_lep4->Fill(pfIso_lep4);
      
      //Quick and dirty look at how applying the pfIso cuts to only the lead 2 pT muons changes what we see, will implement this cut properly after the quick look
  //     if (pfIso_lep1 > pfIso_Cut || pfIso_lep2 > pfIso_Cut){ // || pfIso_lep3 > pfIso_Cut || pfIso_lep4 > pfIso_Cut){
//          std::cout << "FAILED pfIso Cut!" << std::endl;
//          pfIso_Fail_Count += 1;
//          continue;
//       }
      
   //   h_cutflow_allQuadCuts->Fill(4); //here are the quads that survive pfIso cut
   
      //  std::cout << "TREE->quadHasHowManyTrigMatches->at(i)  " << TREE->quadHasHowManyTrigMatches->at(i) << std::endl;
     
      if (TREE->quadHasHowManyTrigMatches->at(i) < 2) {
        continue;
      }
      
       h_cutflow_allQuadCuts->Fill(4); //here are the quads that survive the trigger matching requirement 
     
       //////////////////////////////////////////////
      //end cuts involving overall quad /////////////
      //////////////////////////////////////////////
    
  
     //start pair specific cuts ////
     ///////////////////////////// 
     
     //Note to self: see here: http://www.hep.shef.ac.uk/edaw/PHY206/Site/2012_course_files/phy206rlec7.pdf, equation 6, for rapidity definition (on page 5)
      
      //these are mutually exclusive, one quad can only be one of these 3 things 
      if (TREE->pair_12_34_56->at(i) ==1){
         
         bool Z_first_upsi_phase1_second_pair_12_34_56 = false;
         bool upsi_phase1_first_Z_second_pair_12_34_56 = false; 
        // std::cout << "TREE->pair_12_34_56->at(i) ==1" << std::endl;
         pair_12_34_56_count += 1;
        
         
//         std::cout << (lepton1 + lepton2).M()<< std::endl; 
         
         if ( (lepton1 + lepton2).M() > Z_mass_low && (lepton1 + lepton2).M() < Z_mass_high && (lepton3+lepton4).M()  > upsi_mass_low_phase1 && (lepton3+lepton4).M() < upsi_mass_high_phase1){
            if ( (TREE->lepton1_charge->at(i) + TREE->lepton2_charge->at(i) == 0) && (TREE->lepton3_charge->at(i) + TREE->lepton4_charge->at(i) == 0) ) {
              Z_first_upsi_phase1_second_pair_12_34_56 = true;
              Z_first_upsi_phase1_second_pair_12_34_56_count +=1;
              std::cout << "Z_first_upsi_phase1_second_pair_12_34_56 = true!" <<std::endl; 
              h_cutflow_Z_first_upsi_phase1_second_pair_12_34_56->Fill(1);
           }
        }
         
         
         if ( (lepton1 + lepton2).M() > upsi_mass_low_phase1 && (lepton1 + lepton2).M() < upsi_mass_high_phase1 && (lepton3+lepton4).M()  > Z_mass_low && (lepton3+lepton4).M() < Z_mass_high ){
            if ( (TREE->lepton1_charge->at(i) + TREE->lepton2_charge->at(i) == 0) && (TREE->lepton3_charge->at(i) + TREE->lepton4_charge->at(i) == 0) ) {
              upsi_phase1_first_Z_second_pair_12_34_56 = true;
              upsi_phase1_first_Z_second_pair_12_34_56_count +=1;
              std::cout << "upsi_phase1_first_Z_second_pair_12_34_56 is true!" << std::endl; 
            }
         }
 
 
           
         if (Z_first_upsi_phase1_second_pair_12_34_56) {
            
            //Start Z cuts 
            if (lepton1.Pt() < lead_mu_from_Z_pT_Cut  || lepton2.Pt() < sublead_mu_from_Z_pT_Cut){
                std::cout << "FAILED Z mu Pt Cuts" << std::endl;
                FailureCount += 1;
                continue;
             }
            h_cutflow_Z_first_upsi_phase1_second_pair_12_34_56->Fill(2);
            
            if (fabs(lepton1.Eta()) > mu_from_Z_eta_Cut || fabs(lepton2.Eta()) > mu_from_Z_eta_Cut){
                std::cout << "FAILED Z mu eta Cuts!" << std::endl; 
                FailureCount += 1;
                continue;
            }
            h_cutflow_Z_first_upsi_phase1_second_pair_12_34_56->Fill(3);
            //std::cout << "TREE->lepton1_isTightMuon->at(i): " << TREE->lepton1_isTightMuon->at(i) << std::endl;
         
           if (TREE->lepton1_isTightMuon->at(i) + TREE->lepton2_isTightMuon->at(i) != 2){  //both of them need to be tight, which has a value of 1, 1 +1 =2 
               std::cout << "AT LEAST ONE OF THE MUS FROM A Z WAS NOT TIGHT, FAILED THE Z->MU MU BOTH MU MOST BE TIGHT CUT" << std::endl;
               FailureCount += 1; 
               continue;
           
           } 
           h_cutflow_Z_first_upsi_phase1_second_pair_12_34_56->Fill(4);
           
           if (fabs(TREE->lepton1_impactParameterSignificance->at(i)) > mu_from_Z_3DIPSig_Cut || fabs(TREE->lepton2_impactParameterSignificance->at(i)) > mu_from_Z_3DIPSig_Cut){
               std::cout << "FAILED mu froom Z IP sig cut!" << std::endl;
               FailureCount += 1; 
               continue; 
            }
            h_cutflow_Z_first_upsi_phase1_second_pair_12_34_56->Fill(5);
            
           if (pfIso_lep1 > pfIso_Cut_Mu_from_Z || pfIso_lep2 > pfIso_Cut_Mu_from_Z){
              std::cout << "FAILED particle flow iso cut" << std::endl;
              pfIso_Fail_Count += 1;
              continue;
           }
           
           h_cutflow_Z_first_upsi_phase1_second_pair_12_34_56->Fill(6);
           
           if (TREE->dimuon1vtx->at(i) < mu_mu_from_Z_Prob_Cut){
               std::cout << "FAILED mu_mu_from_Z_Prob_Cut" << std::endl;
               continue; 
           }
           
           h_cutflow_Z_first_upsi_phase1_second_pair_12_34_56->Fill(7);
           
                      
           //End Z cuts
           
           //Start upsi cuts
           
           if (lepton3.Pt() < mu_from_upsi_pT_Cut || lepton4.Pt() < mu_from_upsi_pT_Cut){
               std::cout << "FAILED  upsi mu Pt Cuts" << std::endl;
               FailureCount += 1; 
               continue; 
           }
           h_cutflow_Z_first_upsi_phase1_second_pair_12_34_56->Fill(8);
           
           if (fabs(lepton3.Eta()) > mu_from_upsi_eta_Cut || fabs(lepton4.Eta()) > mu_from_upsi_eta_Cut){
               std::cout << "FAILED upsi  mu eta cuts!" << std::endl; 
               FailureCount +=1;
               continue; 
           }
           
           h_cutflow_Z_first_upsi_phase1_second_pair_12_34_56->Fill(9);
           
           if (  (lepton3 + lepton4).M()    < upsi_mass_low_phase2 || (lepton3 + lepton4).M() > upsi_mass_high_phase2 ){
               std::cout << "FAILED the tighter phase2 upsi mass cuts!" << std::endl;
               FailureCount +=1; 
               continue;  
           }
          
           h_cutflow_Z_first_upsi_phase1_second_pair_12_34_56->Fill(10);
          
           if (TREE->lepton3_isSoftMuon->at(i) + TREE->lepton4_isSoftMuon->at(i) !=2){
               std::cout << "FAILED mu from upsi must be soft cut" << std::endl;
               FailureCount += 1; 
               continue; 
           }
           
           h_cutflow_Z_first_upsi_phase1_second_pair_12_34_56->Fill(11);
          // if ( fabs(lepton3.Rapidity()) > mu_from_upsi_RAPIDITY_Cut || fabs(lepton4.Rapidity()) > mu_from_upsi_RAPIDITY_Cut ){
               if (   fabs((lepton3 + lepton4).Rapidity()) > upsi_RAPIDITY_Cut ){
               std::cout << "FAILED upsi RAPIDITY cut!" << std::endl; 
               FailureCount +=1;
               continue; 
           }
           h_cutflow_Z_first_upsi_phase1_second_pair_12_34_56->Fill(12);
        
           if (TREE->dimuon2vtx->at(i) < mu_mu_from_upsi_Prob_Cut){
               std::cout << "FAILED mu_mu_from_upsi_Prob_Cut" << std::endl;
               continue; 
           }
           
           h_cutflow_Z_first_upsi_phase1_second_pair_12_34_56->Fill(13);
           
           if (pfIso_lep3 > pfIso_Cut_Mu_from_Upsi || pfIso_lep4 > pfIso_Cut_Mu_from_Upsi){
             continue;
           
           }
           
           h_cutflow_Z_first_upsi_phase1_second_pair_12_34_56->Fill(14);
           
           //If we get here, we have a survivor 
           
 //          survivor_Z_first_upsi_phase1_second_pair_12_34_56 = true;
          
         //  Z_mass = (lepton1 + lepton2).M();
         //  upsi_mass = (lepton3 + lepton4).M();
           
           //Here I would have to turn this block into if !doMCTruthMatching , do the stuff below //flagPoodle
           if (!doMCTruthMatching){
          //   std::cout << "KANGAROO" << std::endl; 
           //Z and upsi masses
              temp_Z_mass.push_back((lepton1 + lepton2).M());
              temp_upsi_mass.push_back((lepton3 + lepton4).M());
           
           //Z pT, eta, RAPIDITY, phi
             temp_Z_pT.push_back((lepton1+lepton2).Pt());
             temp_Z_eta.push_back((lepton1+lepton2).Eta());
             temp_Z_RAPIDITY.push_back((lepton1+lepton2).Rapidity());
             temp_Z_phi.push_back((lepton1+lepton2).Phi());
 //          std::cout << "temp_Z_phi.at(0): " << temp_Z_phi.at(0) << std::endl;
          
           //Upsi pT, eta, RAPIDITY, phi
             temp_upsi_pT.push_back((lepton3+lepton4).Pt());
             temp_upsi_eta.push_back((lepton3+lepton4).Eta());
             temp_upsi_RAPIDITY.push_back((lepton3+lepton4).Rapidity());
             temp_upsi_phi.push_back((lepton3+lepton4).Phi());
           
           //Muons from Z, pT, eta, RAPIDITY,phi
            temp_lead_pT_mu_from_Z_pT.push_back(lepton1.Pt());
            temp_lead_pT_mu_from_Z_eta.push_back(lepton1.Eta());
            temp_lead_pT_mu_from_Z_RAPIDITY.push_back(lepton1.Rapidity());
            temp_lead_pT_mu_from_Z_phi.push_back(lepton1.Phi());
            temp_sublead_pT_mu_from_Z_pT.push_back(lepton2.Pt());
            temp_sublead_pT_mu_from_Z_eta.push_back(lepton2.Eta());
            temp_sublead_pT_mu_from_Z_RAPIDITY.push_back(lepton2.Rapidity());
            temp_sublead_pT_mu_from_Z_phi.push_back(lepton2.Phi());
           
           //Muons from upsi, pT, eta, Rapidity, phi
            temp_lead_pT_mu_from_upsi_pT.push_back(lepton3.Pt());
            temp_lead_pT_mu_from_upsi_eta.push_back(lepton3.Eta());
            temp_lead_pT_mu_from_upsi_RAPIDITY.push_back(lepton3.Rapidity());
            temp_lead_pT_mu_from_upsi_phi.push_back(lepton3.Phi());
           
            temp_sublead_pT_mu_from_upsi_pT.push_back(lepton4.Pt());
            temp_sublead_pT_mu_from_upsi_eta.push_back(lepton4.Eta());
            temp_sublead_pT_mu_from_upsi_RAPIDITY.push_back(lepton4.Rapidity());
            temp_sublead_pT_mu_from_upsi_phi.push_back(lepton4.Phi());            
           }
            GotHereCount_Z_first_upsi_phase1_second_pair_12_34_56 += 1;
           //then I would have to write a new if doMCTruthMatching block, do the truth matching // flagPoodle
    
           if (doMCTruthMatching){
            // std::cout << "Poodles! Doing MC Truth Matching" << std::endl;
             //int found1Index = -1;
             
             
             int entriesMC = (TREEMC->fChain)->GetEntries();
             for(int iEntry=0; iEntry<entriesMC; iEntry++) {
                //std::cout << "Koala bear"<< std::endl; 
                (TREEMC->fChain)->GetEntry(iEntry);
                mcSanityCheckCount++;
                
                if (TREEMC->mc_event_number->at(0) == evNumThisQuad && TREEMC->mc_run_number->at(0) == runNumThisQuad && TREEMC->mc_lumi_section->at(0) == LSThisQuad ){
                  //std::cout << "mc_event_number  " << TREEMC->mc_event_number->at(0) << std::endl;
                  //std::cout << "evNumThisQuad  " << evNumThisQuad << std::endl;
                  //std::cout << "dog" << std::endl;
                  bool matched = false;
                  
                  bool found1 = false; //1 short for lep1
                  bool found2 = false; //2 short for lep2
                  int found1Index = -1;
                  for (int i=0; i<(int)TREEMC->truth_Zmuon_pt->size(); i++){
                    TLorentzVector Z_mu_truth;
                    Z_mu_truth.SetPtEtaPhiM(TREEMC->truth_Zmuon_pt->at(i), TREEMC->truth_Zmuon_eta->at(i), TREEMC->truth_Zmuon_phi->at(i), muon_mass);
                    if (Z_mu_truth.DeltaR(lepton1) < deltaRCut){
                    //  std::cout << "EUREKA" << std::endl;
                      found1Index = i;
                      found1 = true;
                    }
                    
                    if (Z_mu_truth.DeltaR(lepton2) < deltaRCut && i != found1Index){
                     // std::cout << "EUREKA 2" << std::endl;
                      found2 = true;
                    }
                  }
                  
                  if (found1 && found2){
                   // std::cout << "Matched Z" << std::endl;
                    matchedZCount++;
                  }
                  
                  bool found3_in_Upsi1 = false; //3 short for lep3
                  bool found4_in_Upsi1 = false; //4 short for lep4
                  int found3_in_Upsi1_Index = -1;
                  
                  for (int j = 0; j < (int)TREEMC->truth_Upsimuon_pt->size(); j++){
                    //std::cout << "Poodle" << std::endl;
                    TLorentzVector Upsi1_mu_truth;
                    Upsi1_mu_truth.SetPtEtaPhiM(TREEMC->truth_Upsimuon_pt->at(j), TREEMC->truth_Upsimuon_eta->at(j), TREEMC->truth_Upsimuon_phi->at(j), muon_mass);
                    if (Upsi1_mu_truth.DeltaR(lepton3) < deltaRCut){
                     // std::cout << "EUREKA 3" << std::endl;
                      found3_in_Upsi1 = true;
                      found3_in_Upsi1_Index = j;
                    }
                    
                    if (Upsi1_mu_truth.DeltaR(lepton4) < deltaRCut && j != found3_in_Upsi1_Index){
                     // std::cout << "EUREKA 4" << std::endl;
                      found4_in_Upsi1 = true;
                    }
                    
                  }
                 
                 // std::cout << "TREEMC->truth_Upsi2muon_pt->size() " << TREEMC->truth_Upsi2muon_pt->size() << std::endl; 
                 
                  bool found3_in_Upsi2 = false; //3 short for lep3
                  bool found4_in_Upsi2 = false; //4 short for lep4
                  int found3_in_Upsi2_Index =  -1;
                  
                  for (int k = 0; k < (int)TREEMC->truth_Upsi2muon_pt->size(); k++){
                    //std::cout << "Great Ape" << std::endl;
                    TLorentzVector Upsi2_mu_truth;
                    Upsi2_mu_truth.SetPtEtaPhiM(TREEMC->truth_Upsi2muon_pt->at(k), TREEMC->truth_Upsi2muon_eta->at(k), TREEMC->truth_Upsi2muon_phi->at(k), muon_mass);
                    if (Upsi2_mu_truth.DeltaR(lepton3) < deltaRCut){
                     // std::cout << "found3_in_Upsi2" << std::endl;
                      found3_in_Upsi2 = true;
                      found3_in_Upsi2_Index = k;
                    }
                    
                    if (Upsi2_mu_truth.DeltaR(lepton4) < deltaRCut && k != found3_in_Upsi2_Index){
                     // std::cout << "found4_in_Upsi2" << std::endl;
                      found4_in_Upsi2 = true;
                    }
                  }
                  
                  bool found3_in_Upsi3 = false; //3 short for lep3
                  bool found4_in_Upsi3 = false; //4 short for lep4
                  int found3_in_Upsi3_Index = -1;
                  
                 // std::cout << "TREEMC->truth_Upsi3muon_pt->size() " << TREEMC->truth_Upsi3muon_pt->size() << std::endl;
                 for (int l = 0; l < (int)TREEMC->truth_Upsi3muon_pt->size(); l++){
                   //std::cout << "Ron Swanson" << std::endl; 
                   TLorentzVector Upsi3_mu_truth;
                   Upsi3_mu_truth.SetPtEtaPhiM(TREEMC->truth_Upsi3muon_pt->at(l), TREEMC->truth_Upsi3muon_eta->at(l), TREEMC->truth_Upsi3muon_phi->at(l), muon_mass);
                   if (Upsi3_mu_truth.DeltaR(lepton3) < deltaRCut){
                    // std::cout << "found3_in_Upsi3 " << std::endl;
                     found3_in_Upsi3 = true;
                     found3_in_Upsi3_Index = l;
                   }
                   
                   if (Upsi3_mu_truth.DeltaR(lepton4) < deltaRCut && l != found3_in_Upsi3_Index){
                    // std::cout << "found4_in_Upsi3" << std::endl;
                     found4_in_Upsi3 = true;
                   }
                 } 
                 
                 //Note: I really shouldn't need the protection of the !found3_in_Upsi2, etc, since if the Upsi1 vector quantities have size not equal to 0
                 //the other two upsi type vectors should have size 0, based on the way the MC was made (they are Z + one upsi_type samples) 
                 //trying the protection out though because it doesn't cost much
                 //Thank you Bjorn for valuable insight into this point
                 if (found1 && found2 && found3_in_Upsi1 && found4_in_Upsi1 && !found3_in_Upsi2 && !found4_in_Upsi2 && !found3_in_Upsi3 && !found4_in_Upsi3){
                   upsi_type = 1;
                   //std::cout << "upsi_type: " << upsi_type << std::endl;
                   matched = true;
                   matchedCount++;
                 }
                 
                 if (found1 && found2 && found3_in_Upsi2 && found4_in_Upsi2 && !found3_in_Upsi1 && !found4_in_Upsi1 && !found3_in_Upsi3 && !found4_in_Upsi3){
                   upsi_type = 2;
                   //std::cout << "upsi_type: " << upsi_type << std::endl; 
                   matched = true;
                   matchedCount++;
                 }
                 
                 if  (found1 && found2 && found3_in_Upsi3 && found4_in_Upsi3 && !found3_in_Upsi1 && !found4_in_Upsi1 && !found3_in_Upsi2 && !found4_in_Upsi2){
                   upsi_type = 3;
                  // std::cout << "upsi_type: " << upsi_type << std::endl; 
                   matched = true;
                   matchedCount++;
                 }
                 
                 if (matched){
                   //Z and upsi masses
                   temp_Z_mass.push_back((lepton1 + lepton2).M());
                   temp_upsi_mass.push_back((lepton3 + lepton4).M());
           
                  //Z pT, eta, RAPIDITY, phi
                  temp_Z_pT.push_back((lepton1+lepton2).Pt());
                  temp_Z_eta.push_back((lepton1+lepton2).Eta());
                  temp_Z_RAPIDITY.push_back((lepton1+lepton2).Rapidity());
                  temp_Z_phi.push_back((lepton1+lepton2).Phi());
 //          std::cout << "temp_Z_phi.at(0): " << temp_Z_phi.at(0) << std::endl;
          
                 //Upsi pT, eta, RAPIDITY, phi
                  temp_upsi_pT.push_back((lepton3+lepton4).Pt());
                  temp_upsi_eta.push_back((lepton3+lepton4).Eta());
                  temp_upsi_RAPIDITY.push_back((lepton3+lepton4).Rapidity());
                  temp_upsi_phi.push_back((lepton3+lepton4).Phi());
           
                //Muons from Z, pT, eta, RAPIDITY,phi
                 temp_lead_pT_mu_from_Z_pT.push_back(lepton1.Pt());
                 temp_lead_pT_mu_from_Z_eta.push_back(lepton1.Eta());
                 temp_lead_pT_mu_from_Z_RAPIDITY.push_back(lepton1.Rapidity());
                 temp_lead_pT_mu_from_Z_phi.push_back(lepton1.Phi());
                 temp_sublead_pT_mu_from_Z_pT.push_back(lepton2.Pt());
                 temp_sublead_pT_mu_from_Z_eta.push_back(lepton2.Eta());
                 temp_sublead_pT_mu_from_Z_RAPIDITY.push_back(lepton2.Rapidity());
                 temp_sublead_pT_mu_from_Z_phi.push_back(lepton2.Phi());
           
                //Muons from upsi, pT, eta, Rapidity, phi
                 temp_lead_pT_mu_from_upsi_pT.push_back(lepton3.Pt());
                 temp_lead_pT_mu_from_upsi_eta.push_back(lepton3.Eta());
                 temp_lead_pT_mu_from_upsi_RAPIDITY.push_back(lepton3.Rapidity());
                 temp_lead_pT_mu_from_upsi_phi.push_back(lepton3.Phi());
           
                 temp_sublead_pT_mu_from_upsi_pT.push_back(lepton4.Pt());
                 temp_sublead_pT_mu_from_upsi_eta.push_back(lepton4.Eta());
                 temp_sublead_pT_mu_from_upsi_RAPIDITY.push_back(lepton4.Rapidity());
                 temp_sublead_pT_mu_from_upsi_phi.push_back(lepton4.Phi()); 
                 
                 temp_upsi_type.push_back(upsi_type);
                
                 }
                 
                }
                
             }
             
           } //Keep me, I am important and separate from the new stuff you're trying to add

            
      //      
           
         
         }
         
        
      
         
         if  (upsi_phase1_first_Z_second_pair_12_34_56) { 
             
            //start Z cuts  
            if (lepton3.Pt() < lead_mu_from_Z_pT_Cut || lepton4.Pt() <sublead_mu_from_Z_pT_Cut){
                std::cout << "FAILED Z mu Pt Cuts!" << std::endl;
                continue; 
            }
            
            if (fabs(lepton3.Eta()) > mu_from_Z_eta_Cut || fabs(lepton4.Eta()) > mu_from_Z_eta_Cut){
                std::cout << "FAILED Z mu eta Cuts!" << std::endl; 
                continue;
            }
            
            if (TREE->lepton3_isTightMuon->at(i) + TREE->lepton4_isTightMuon->at(i) != 2){  //both of them need to be tight, which has a value of 1, 1 +1 =2 
               std::cout << "AT LEAST ONE OF THE MUS FROM A Z WAS NOT TIGHT, FAILED THE Z->MU MU BOTH MU MOST BE TIGHT CUT" << std::endl;
               continue;
           
           } 
           
            if (fabs(TREE->lepton3_impactParameterSignificance->at(i)) > mu_from_Z_3DIPSig_Cut || fabs(TREE->lepton4_impactParameterSignificance->at(i)) > mu_from_Z_3DIPSig_Cut){
               std::cout << "FAILED mu froom Z IP sig cut!" << std::endl;
               continue; 
            }
            
            if (pfIso_lep3 > pfIso_Cut_Mu_from_Z || pfIso_lep4 > pfIso_Cut_Mu_from_Z){
              std::cout << "FAILED particle flow iso cut" << std::endl;
              pfIso_Fail_Count += 1;
              continue;
           }
            
            
            if (TREE->dimuon2vtx->at(i) < mu_mu_from_Z_Prob_Cut){
               std::cout << "FAILED mu_mu_from_Z_Prob_Cut" << std::endl;
               continue; 
            }
            
            //End Z cuts
            
            //begin upsi cuts 
            if (lepton1.Pt() < mu_from_upsi_pT_Cut || lepton2.Pt() < mu_from_upsi_pT_Cut){
                std::cout << "FAILED upsi mu Pt cuts!" << std::endl;
                continue;
            }
            
            if (fabs(lepton1.Eta()) > mu_from_upsi_eta_Cut || fabs(lepton2.Eta()) > mu_from_upsi_eta_Cut){
                std::cout << "FAILED upsi mu eta cuts!" << std::endl; 
                continue;
            }
            
            if ( (lepton1 + lepton2).M() < upsi_mass_low_phase2 || (lepton1 + lepton2).M() > upsi_mass_high_phase2 ){
                std::cout << "FAILED the tighter phase2 upsi mass cuts!" << std::endl;
                continue; 
            }
            
            if (TREE->lepton1_isSoftMuon->at(i) + TREE->lepton2_isSoftMuon->at(i) !=2){
               std::cout << "FAILED mu from upsi must be soft cut" << std::endl;
               continue; 
           }
            
            if //( fabs(lepton1.Rapidity()) > mu_from_upsi_RAPIDITY_Cut || fabs(lepton2.Rapidity()) > mu_from_upsi_RAPIDITY_Cut ){
                
               ( fabs(  (lepton1 + lepton2).Rapidity() ) > upsi_RAPIDITY_Cut ) {
                
                std::cout << "FAILED  upsi RAPIDITY cut!" << std::endl;
                continue;           
            }
            
            if (TREE->dimuon1vtx->at(i) < mu_mu_from_upsi_Prob_Cut){
                std::cout << "FAILED mu_mu_from_upsi_Prob_Cut" << std::endl;
                continue; 
            }
            
            if (pfIso_lep1 > pfIso_Cut_Mu_from_Upsi || pfIso_lep2 > pfIso_Cut_Mu_from_Upsi){
              continue; 
            }
        
             //If we get here, we have a survivor
           // Z_mass_ = (lepton3 + lepton4).M();
          // upsi_mass = (lepton1 + lepton2).M();
           
           if (!doMCTruthMatching){
           
             temp_Z_mass.push_back((lepton3 + lepton4).M());
             temp_upsi_mass.push_back((lepton1 + lepton2).M());
           
           //pT, eta, RAPIDITY, phi
             temp_Z_pT.push_back((lepton3+lepton4).Pt());
             temp_Z_eta.push_back((lepton3+lepton4).Eta());
             temp_Z_RAPIDITY.push_back((lepton3+lepton4).Rapidity());
             temp_Z_phi.push_back((lepton3+lepton4).Phi());
           
             temp_upsi_pT.push_back((lepton1+lepton2).Pt());
             temp_upsi_eta.push_back((lepton1+lepton2).Eta());
             temp_upsi_RAPIDITY.push_back((lepton1+lepton2).Rapidity());
             temp_upsi_phi.push_back((lepton1+lepton2).Phi());
           
             temp_lead_pT_mu_from_Z_pT.push_back(lepton3.Pt());
             temp_lead_pT_mu_from_Z_eta.push_back(lepton3.Eta());
             temp_lead_pT_mu_from_Z_RAPIDITY.push_back(lepton3.Rapidity());
             temp_lead_pT_mu_from_Z_phi.push_back(lepton3.Phi());
           
             temp_sublead_pT_mu_from_Z_pT.push_back(lepton4.Pt());
             temp_sublead_pT_mu_from_Z_eta.push_back(lepton4.Eta());
             temp_sublead_pT_mu_from_Z_RAPIDITY.push_back(lepton4.Rapidity());
             temp_sublead_pT_mu_from_Z_phi.push_back(lepton4.Phi());
           
             temp_lead_pT_mu_from_upsi_pT.push_back(lepton1.Pt());
             temp_lead_pT_mu_from_upsi_eta.push_back(lepton1.Eta());
             temp_lead_pT_mu_from_upsi_RAPIDITY.push_back(lepton1.Rapidity());
             temp_lead_pT_mu_from_upsi_phi.push_back(lepton1.Phi());
           
             temp_sublead_pT_mu_from_upsi_pT.push_back(lepton2.Pt());
             temp_sublead_pT_mu_from_upsi_eta.push_back(lepton2.Eta());
             temp_sublead_pT_mu_from_upsi_RAPIDITY.push_back(lepton2.Rapidity());
             temp_sublead_pT_mu_from_upsi_phi.push_back(lepton2.Phi());
             }
           
           //flag Begin MC Truth Matching section here
           if (doMCTruthMatching){
             std::cout << "Poodles! Doing MC Truth Matching!" << std::endl; 
             
             int entriesMC = (TREEMC->fChain)->GetEntries();
             for(int iEntry=0; iEntry<entriesMC; iEntry++) {
                //std::cout << "Koala bear"<< std::endl; 
                (TREEMC->fChain)->GetEntry(iEntry);
                if (TREEMC->mc_event_number->at(0) == evNumThisQuad && TREEMC->mc_run_number->at(0) == runNumThisQuad && TREEMC->mc_lumi_section->at(0) == LSThisQuad ){
                  std::cout << "Poodles again!" << std::endl; 
                  
                  bool matched = false; 
                  
                  bool found3 = false; //short for lep3
                  bool found4 = false; //short for lep4
                  int found3Index = -1; 
                  
                  for (int i=0; i<(int)TREEMC->truth_Zmuon_pt->size(); i++){
                    TLorentzVector Z_mu_truth;
                    Z_mu_truth.SetPtEtaPhiM(TREEMC->truth_Zmuon_pt->at(i), TREEMC->truth_Zmuon_eta->at(i), TREEMC->truth_Zmuon_phi->at(i), muon_mass);
                    if (Z_mu_truth.DeltaR(lepton3) < deltaRCut){
                    //  std::cout << "EUREKA" << std::endl;
                      found3Index = i;
                      found3 = true;
                    }
                    
                    if (Z_mu_truth.DeltaR(lepton4) < deltaRCut && i != found3Index){
                     // std::cout << "EUREKA 2" << std::endl;
                      found4 = true;
                    }
                  }
                  
                  if (found3 && found4){
                   // std::cout << "Matched Z" << std::endl;
                    matchedZCount++;
                  }
                  
                  bool found1_in_Upsi1 = false; //1 short for lep1
                  bool found2_in_Upsi1 = false; //2 short for lep2
                  int found1_in_Upsi1_Index = -1;
                  
                  for (int j = 0; j < (int)TREEMC->truth_Upsimuon_pt->size(); j++){
                    //std::cout << "Poodle" << std::endl;
                    TLorentzVector Upsi1_mu_truth;
                    Upsi1_mu_truth.SetPtEtaPhiM(TREEMC->truth_Upsimuon_pt->at(j), TREEMC->truth_Upsimuon_eta->at(j), TREEMC->truth_Upsimuon_phi->at(j), muon_mass);
                    if (Upsi1_mu_truth.DeltaR(lepton1) < deltaRCut){
                     // std::cout << "EUREKA 3" << std::endl;
                      found1_in_Upsi1 = true;
                      found1_in_Upsi1_Index = j;
                    }
                    
                    if (Upsi1_mu_truth.DeltaR(lepton2) < deltaRCut && j != found1_in_Upsi1_Index){
                     // std::cout << "EUREKA 4" << std::endl;
                      found2_in_Upsi1 = true;
                    }
                    
                  }
                  
                  bool found1_in_Upsi2 = false;
                  bool found2_in_Upsi2 = false;
                  int  found1_in_Upsi2_Index = -1; 
                  
                  for (int k = 0; k < (int)TREEMC->truth_Upsi2muon_pt->size(); k++){
                    //std::cout << "Great Ape" << std::endl;
                    TLorentzVector Upsi2_mu_truth;
                    Upsi2_mu_truth.SetPtEtaPhiM(TREEMC->truth_Upsi2muon_pt->at(k), TREEMC->truth_Upsi2muon_eta->at(k), TREEMC->truth_Upsi2muon_phi->at(k), muon_mass);
                    if (Upsi2_mu_truth.DeltaR(lepton1) < deltaRCut){
                     // std::cout << "found1_in_Upsi2" << std::endl;
                      found1_in_Upsi2 = true;
                      found1_in_Upsi2_Index = k;
                    }
                    
                    if (Upsi2_mu_truth.DeltaR(lepton2) < deltaRCut && k != found1_in_Upsi2_Index){
                     // std::cout << "found2_in_Upsi2" << std::endl;
                      found2_in_Upsi2 = true;
                    }
                  }
                  
                  bool found1_in_Upsi3 = false;
                  bool found2_in_Upsi3 = false;
                  int found1_in_Upsi3_Index = -1;
                  
                  for (int l = 0; l < (int)TREEMC->truth_Upsi3muon_pt->size(); l++){
                   //std::cout << "Ron Swanson" << std::endl; 
                   TLorentzVector Upsi3_mu_truth;
                   Upsi3_mu_truth.SetPtEtaPhiM(TREEMC->truth_Upsi3muon_pt->at(l), TREEMC->truth_Upsi3muon_eta->at(l), TREEMC->truth_Upsi3muon_phi->at(l), muon_mass);
                   if (Upsi3_mu_truth.DeltaR(lepton1) < deltaRCut){
                    // std::cout << "found1_in_Upsi3 " << std::endl;
                     found1_in_Upsi3 = true;
                     found1_in_Upsi3_Index = l;
                   }
                   
                   if (Upsi3_mu_truth.DeltaR(lepton2) < deltaRCut && l != found1_in_Upsi3_Index){
                    // std::cout << "found2_in_Upsi3" << std::endl;
                     found2_in_Upsi3 = true;
                   }
                 } 
                 
                  if (found3 && found4 && found1_in_Upsi1 && found2_in_Upsi1 && !found1_in_Upsi2 && !found2_in_Upsi2 && !found1_in_Upsi3 && !found2_in_Upsi3){
                    upsi_type = 1;
                    matched = true;
                    matchedCount++;
                  
                  }
                  
                  if (found3 && found4 && found1_in_Upsi2 && found2_in_Upsi2 && !found1_in_Upsi1 && !found2_in_Upsi1 && !found1_in_Upsi3 && !found2_in_Upsi3){
                    upsi_type = 2;
                    matched = true;
                    matchedCount++;
                  }
                  
                  if (found3 && found4 && found1_in_Upsi3 && found2_in_Upsi3 && !found1_in_Upsi1 && !found2_in_Upsi1 && !found1_in_Upsi2 && !found2_in_Upsi2){
                    upsi_type = 3;
                    matched = true;
                    matchedCount++;
                  }
                  
                  if (matched){
                    
                    temp_Z_mass.push_back((lepton3 + lepton4).M());
                    temp_upsi_mass.push_back((lepton1 + lepton2).M());
           
                    //pT, eta, RAPIDITY, phi
                    temp_Z_pT.push_back((lepton3+lepton4).Pt());
                    temp_Z_eta.push_back((lepton3+lepton4).Eta());
                    temp_Z_RAPIDITY.push_back((lepton3+lepton4).Rapidity());
                    temp_Z_phi.push_back((lepton3+lepton4).Phi());
           
                    temp_upsi_pT.push_back((lepton1+lepton2).Pt());
                    temp_upsi_eta.push_back((lepton1+lepton2).Eta());
                    temp_upsi_RAPIDITY.push_back((lepton1+lepton2).Rapidity());
                    temp_upsi_phi.push_back((lepton1+lepton2).Phi());
           
                    temp_lead_pT_mu_from_Z_pT.push_back(lepton3.Pt());
                    temp_lead_pT_mu_from_Z_eta.push_back(lepton3.Eta());
                    temp_lead_pT_mu_from_Z_RAPIDITY.push_back(lepton3.Rapidity());
                    temp_lead_pT_mu_from_Z_phi.push_back(lepton3.Phi());
           
                    temp_sublead_pT_mu_from_Z_pT.push_back(lepton4.Pt());
                    temp_sublead_pT_mu_from_Z_eta.push_back(lepton4.Eta());
                    temp_sublead_pT_mu_from_Z_RAPIDITY.push_back(lepton4.Rapidity());
                    temp_sublead_pT_mu_from_Z_phi.push_back(lepton4.Phi());
           
                   temp_lead_pT_mu_from_upsi_pT.push_back(lepton1.Pt());
                   temp_lead_pT_mu_from_upsi_eta.push_back(lepton1.Eta());
                   temp_lead_pT_mu_from_upsi_RAPIDITY.push_back(lepton1.Rapidity());
                   temp_lead_pT_mu_from_upsi_phi.push_back(lepton1.Phi());
           
                   temp_sublead_pT_mu_from_upsi_pT.push_back(lepton2.Pt());
                   temp_sublead_pT_mu_from_upsi_eta.push_back(lepton2.Eta());
                   temp_sublead_pT_mu_from_upsi_RAPIDITY.push_back(lepton2.Rapidity());
                   temp_sublead_pT_mu_from_upsi_phi.push_back(lepton2.Phi());
                  
                   temp_upsi_type.push_back(upsi_type);
                  
                  
                  }
                }
                
                
                
             }
             
           }  
             
         } //Keep me, I am not related to the new stuff you are adding in the MC Truth Matching section //I should be on my own and not paired up with a bracket from the doMCTruthMatching section
    
      }
      
       
      if (TREE->pair_13_24_56->at(i) == 1){
         
         bool Z_first_upsi_phase1_second_pair_13_24_56 = false;
         bool upsi_phase1_first_Z_second_pair_13_24_56 = false; 
         //std::cout << "TREE->pair_13_24_56->at(i) == 1" << std::endl; 
         pair_13_24_56_count += 1;
         
         if ( (lepton1 + lepton3).M() > Z_mass_low && (lepton1 + lepton3).M() < Z_mass_high && (lepton2+lepton4).M()  > upsi_mass_low_phase1 && (lepton2+lepton4).M() < upsi_mass_high_phase1){
           if  ( (TREE->lepton1_charge->at(i) + TREE->lepton3_charge->at(i) == 0) && (TREE->lepton2_charge->at(i) + TREE->lepton4_charge->at(i) == 0) ) {
           
            Z_first_upsi_phase1_second_pair_13_24_56 = true;
           Z_first_upsi_phase1_second_pair_13_24_56_count +=1;
           std::cout << "Z_first_upsi_phase1_second_pair_13_24_56 = true!" <<std::endl; 
           
           }
        }
         
         if ( (lepton1 + lepton3).M() > upsi_mass_low_phase1 && (lepton1 + lepton3).M() < upsi_mass_high_phase1 && (lepton2+lepton4).M()  > Z_mass_low && (lepton2+lepton4).M() < Z_mass_high ){
            if  ( (TREE->lepton1_charge->at(i) + TREE->lepton3_charge->at(i) == 0) && (TREE->lepton2_charge->at(i) + TREE->lepton4_charge->at(i) == 0) ) {
              upsi_phase1_first_Z_second_pair_13_24_56 = true;
              upsi_phase1_first_Z_second_pair_13_24_56_count +=1;
              std::cout << "upsi_phase1_first_Z_second_pair_13_24_56 is true!" << std::endl; 
            }
      
         }
    
        
         if (Z_first_upsi_phase1_second_pair_13_24_56){
            //std::cout << "PLACEHOLDER!" << std::endl;      
            if (lepton1.Pt() < lead_mu_from_Z_pT_Cut  || lepton3.Pt() < sublead_mu_from_Z_pT_Cut){
                std::cout << "FAILED Z mu Pt Cuts" << std::endl;
                continue;
             }    
            
            if (fabs(lepton1.Eta()) > mu_from_Z_eta_Cut || fabs(lepton3.Eta()) > mu_from_Z_eta_Cut){
                std::cout << "FAILED Z mu eta Cuts!" << std::endl; 
                continue;
            }
            
            if (TREE->lepton1_isTightMuon->at(i) + TREE->lepton3_isTightMuon->at(i) != 2){  //both of them need to be tight, tight has a value of 1, 1 +1 =2 
               std::cout << "AT LEAST ONE OF THE MUS FROM A Z WAS NOT TIGHT, FAILED THE Z->MU MU BOTH MU MOST BE TIGHT CUT" << std::endl;
               continue;
           
           } 
           
            if (fabs(TREE->lepton1_impactParameterSignificance->at(i)) > mu_from_Z_3DIPSig_Cut || fabs(TREE->lepton3_impactParameterSignificance->at(i)) > mu_from_Z_3DIPSig_Cut){
               std::cout << "FAILED mu froom Z IP sig cut!" << std::endl;
               continue; 
            }
            
            if (pfIso_lep1 > pfIso_Cut_Mu_from_Z || pfIso_lep3 > pfIso_Cut_Mu_from_Z){
              std::cout << "FAILED particle flow iso cut" << std::endl;
              pfIso_Fail_Count += 1;
              continue;
           }
           
           if (TREE->dimuon1vtx->at(i) < mu_mu_from_Z_Prob_Cut){
               std::cout << "FAILED mu_mu_from_Z_Prob_Cut" << std::endl;
               continue; 
           }
           
            
            //end Z cuts
            
            //start upsi cuts 
            if (lepton2.Pt() < mu_from_upsi_pT_Cut || lepton4.Pt() < mu_from_upsi_pT_Cut){
              std::cout << "FAILED upsi mu pT cuts!" << std::endl; 
              continue;  
            }
            
            if (fabs(lepton2.Eta()) > mu_from_upsi_eta_Cut || fabs(lepton4.Eta()) > mu_from_upsi_eta_Cut){
              std::cout << "FAILED upsi mu eta cuts!" << std::endl; 
              continue; 
            }
            
            if ( (lepton2 + lepton4).M() < upsi_mass_low_phase2 || (lepton2 + lepton4).M() > upsi_mass_high_phase2){
              std::cout << "FAILED the tighter phase2 upsi mass cuts!" << std::endl;
              continue;  
            }
            
            if (TREE->lepton2_isSoftMuon->at(i) + TREE->lepton4_isSoftMuon->at(i) !=2){
                std::cout << "FAILED mu from upsi must be soft cut!" << std::endl;
                continue; 
             
            }
            
            //if //( fabs(lepton2.Rapidity()) > mu_from_upsi_RAPIDITY_Cut || fabs(lepton4.Rapidity()) > mu_from_upsi_RAPIDITY_Cut ){
            if  ( fabs((lepton2 + lepton4).Rapidity()) > upsi_RAPIDITY_Cut){
                std::cout << "FAILED upsi RAPIDITY cut!" << std::endl;
                continue; 
            
            }
            
            if (TREE->dimuon2vtx->at(i) < mu_mu_from_upsi_Prob_Cut){
               std::cout << "FAILED mu_mu_from_upsi_Prob_Cut" << std::endl;
               continue; 
           }
           
           if (pfIso_lep2 > pfIso_Cut_Mu_from_Upsi || pfIso_lep4 > pfIso_Cut_Mu_from_Upsi){
              continue;
           }
            
            if (!doMCTruthMatching){
                temp_Z_mass.push_back((lepton1 + lepton3).M());
                temp_upsi_mass.push_back((lepton2 + lepton4).M());
              
              //Pt, eta, Rapidity, Phi
                temp_Z_pT.push_back((lepton1+lepton3).Pt());
                temp_Z_eta.push_back((lepton1+lepton3).Eta());
                temp_Z_RAPIDITY.push_back((lepton1+lepton3).Rapidity());
                temp_Z_phi.push_back((lepton1+lepton3).Phi());
              
                temp_upsi_pT.push_back((lepton2+lepton4).Pt());
                temp_upsi_eta.push_back((lepton2+lepton4).Eta());
                temp_upsi_RAPIDITY.push_back((lepton2+lepton4).Rapidity());
                temp_upsi_phi.push_back((lepton2+lepton4).Phi());
              
                temp_lead_pT_mu_from_Z_pT.push_back(lepton1.Pt());
                temp_lead_pT_mu_from_Z_eta.push_back(lepton1.Eta());
                temp_lead_pT_mu_from_Z_RAPIDITY.push_back(lepton1.Rapidity());
                temp_lead_pT_mu_from_Z_phi.push_back(lepton1.Phi());
              
                temp_sublead_pT_mu_from_Z_pT.push_back(lepton3.Pt());
                temp_sublead_pT_mu_from_Z_eta.push_back(lepton3.Eta());
                temp_sublead_pT_mu_from_Z_RAPIDITY.push_back(lepton3.Rapidity());
                temp_sublead_pT_mu_from_Z_phi.push_back(lepton3.Phi());
              
                temp_lead_pT_mu_from_upsi_pT.push_back(lepton2.Pt());
                temp_lead_pT_mu_from_upsi_eta.push_back(lepton2.Eta());
                temp_lead_pT_mu_from_upsi_RAPIDITY.push_back(lepton2.Rapidity());
                temp_lead_pT_mu_from_upsi_phi.push_back(lepton2.Phi());
              
                temp_sublead_pT_mu_from_upsi_pT.push_back(lepton4.Pt());
                temp_sublead_pT_mu_from_upsi_eta.push_back(lepton4.Eta());
                temp_sublead_pT_mu_from_upsi_RAPIDITY.push_back(lepton4.Rapidity());
                temp_sublead_pT_mu_from_upsi_phi.push_back(lepton4.Phi());
              }
            
            //flag Poodles 
            if (doMCTruthMatching){
              std::cout << "Poodles 1" << std::endl; 
              
              int entriesMC = (TREEMC->fChain)->GetEntries();
              
              for(int iEntry=0; iEntry<entriesMC; iEntry++) {
                //std::cout << "Koala bear"<< std::endl; 
                (TREEMC->fChain)->GetEntry(iEntry);
                if (TREEMC->mc_event_number->at(0) == evNumThisQuad && TREEMC->mc_run_number->at(0) == runNumThisQuad && TREEMC->mc_lumi_section->at(0) == LSThisQuad ){
                  std::cout << "Poodles again!" << std::endl; 
                  
                  bool matched = false; 
                  
                  bool found1 = false; //short for lep1
                  bool found3 = false; //short for lep3
                  int found1Index = -1;
                  
                  for (int i=0; i<(int)TREEMC->truth_Zmuon_pt->size(); i++){
                    TLorentzVector Z_mu_truth;
                    Z_mu_truth.SetPtEtaPhiM(TREEMC->truth_Zmuon_pt->at(i), TREEMC->truth_Zmuon_eta->at(i), TREEMC->truth_Zmuon_phi->at(i), muon_mass);
                    if (Z_mu_truth.DeltaR(lepton1) < deltaRCut){
                      std::cout << "EUREKA" << std::endl;
                      found1Index = i;
                      found1 = true;
                    }
                    
                    if (Z_mu_truth.DeltaR(lepton3) < deltaRCut && i != found1Index){
                      std::cout << "EUREKA 2" << std::endl;
                      found3 = true;
                    }
                  }
                  
                  if (found1 && found3){
                    std::cout << "Matched Z" << std::endl;
                    matchedZCount++;
                  }
                  
                  bool found2_in_Upsi1 = false; //short for lep2
                  bool found4_in_Upsi1 = false; //short for lep4
                  int found2_in_Upsi1_Index = -1;
                  
                  for (int j = 0; j < (int)TREEMC->truth_Upsimuon_pt->size(); j++){
                    //std::cout << "Poodle" << std::endl;
                    TLorentzVector Upsi1_mu_truth;
                    Upsi1_mu_truth.SetPtEtaPhiM(TREEMC->truth_Upsimuon_pt->at(j), TREEMC->truth_Upsimuon_eta->at(j), TREEMC->truth_Upsimuon_phi->at(j), muon_mass);
                    if (Upsi1_mu_truth.DeltaR(lepton2) < deltaRCut){
                      std::cout << "EUREKA 3" << std::endl;
                      found2_in_Upsi1 = true;
                      found2_in_Upsi1_Index = j;
                    }
                    
                    if (Upsi1_mu_truth.DeltaR(lepton4) < deltaRCut && j != found2_in_Upsi1_Index){
                      std::cout << "EUREKA 4" << std::endl;
                      found4_in_Upsi1 = true;
                    }
                    
                  }
                  
                  bool found2_in_Upsi2 = false; //short for lep2
                  bool found4_in_Upsi2 = false; //short for lep4
                  int found2_in_Upsi2_Index = -1;
                  
                  for (int k = 0; k < (int)TREEMC->truth_Upsi2muon_pt->size(); k++){
                    std::cout << "Great Ape" << std::endl;
                    TLorentzVector Upsi2_mu_truth;
                    Upsi2_mu_truth.SetPtEtaPhiM(TREEMC->truth_Upsi2muon_pt->at(k), TREEMC->truth_Upsi2muon_eta->at(k), TREEMC->truth_Upsi2muon_phi->at(k), muon_mass);
                    if (Upsi2_mu_truth.DeltaR(lepton2) < deltaRCut){
                     // std::cout << "found2_in_Upsi2" << std::endl;
                      found2_in_Upsi2 = true;
                      found2_in_Upsi2_Index = k;
                    }
                    
                    if (Upsi2_mu_truth.DeltaR(lepton4) < deltaRCut && k != found2_in_Upsi2_Index){
                     // std::cout << "found4_in_Upsi2" << std::endl;
                      found4_in_Upsi2 = true;
                    }
                  }
                  
                  bool found2_in_Upsi3 = false;
                  bool found4_in_Upsi3 = false;
                  int found2_in_Upsi3_Index = -1;
                  
                  for (int l = 0; l < (int)TREEMC->truth_Upsi3muon_pt->size(); l++){
                   //std::cout << "Ron Swanson" << std::endl; 
                   TLorentzVector Upsi3_mu_truth;
                   Upsi3_mu_truth.SetPtEtaPhiM(TREEMC->truth_Upsi3muon_pt->at(l), TREEMC->truth_Upsi3muon_eta->at(l), TREEMC->truth_Upsi3muon_phi->at(l), muon_mass);
                   if (Upsi3_mu_truth.DeltaR(lepton2) < deltaRCut){
                    // std::cout << "found2_in_Upsi3 " << std::endl;
                     found2_in_Upsi3 = true;
                     found2_in_Upsi3_Index = l;
                   }
                   
                   if (Upsi3_mu_truth.DeltaR(lepton4) < deltaRCut && l != found2_in_Upsi3_Index){
                    // std::cout << "found2_in_Upsi3" << std::endl;
                     found4_in_Upsi3 = true;
                   }
                 } 
                  
                  if (found1 && found3 && found2_in_Upsi1 && found4_in_Upsi1 && !found2_in_Upsi2 && !found4_in_Upsi2 && !found2_in_Upsi3 && !found4_in_Upsi3){
                    upsi_type = 1;
                    matched = true;
                    matchedCount++;
                  }
                  
                  //Note to self 3 Dec. 2021: start here tomorrow!
                  if (found1 && found3 && found2_in_Upsi2 && found4_in_Upsi2 && !found2_in_Upsi1 && !found4_in_Upsi1 && !found2_in_Upsi3 && !found4_in_Upsi3){
                    upsi_type = 2;
                    matched = true;
                    matchedCount++;
                  }
                  
                  if (found1 && found3 && found2_in_Upsi3 && found4_in_Upsi3 && !found2_in_Upsi1 && !found4_in_Upsi1 && !found2_in_Upsi2 && !found4_in_Upsi2){
                    upsi_type = 3;
                    matched = true;
                    matchedCount++;
                  }
                  
                  if (matched){
                    std::cout << "DEBUG 1" << std::endl; 
                    temp_Z_mass.push_back((lepton1 + lepton3).M());
                    temp_upsi_mass.push_back((lepton2 + lepton4).M());
                    
                    std::cout << "DEBUG 2" << std::endl; 
                   //Pt, eta, Rapidity, Phi
                    temp_Z_pT.push_back((lepton1+lepton3).Pt());
                    temp_Z_eta.push_back((lepton1+lepton3).Eta());
                    temp_Z_RAPIDITY.push_back((lepton1+lepton3).Rapidity());
                    temp_Z_phi.push_back((lepton1+lepton3).Phi());
              
                    temp_upsi_pT.push_back((lepton2+lepton4).Pt());
                    temp_upsi_eta.push_back((lepton2+lepton4).Eta());
                    temp_upsi_RAPIDITY.push_back((lepton2+lepton4).Rapidity());
                    temp_upsi_phi.push_back((lepton2+lepton4).Phi());
              
                    temp_lead_pT_mu_from_Z_pT.push_back(lepton1.Pt());
                    temp_lead_pT_mu_from_Z_eta.push_back(lepton1.Eta());
                    temp_lead_pT_mu_from_Z_RAPIDITY.push_back(lepton1.Rapidity());
                    temp_lead_pT_mu_from_Z_phi.push_back(lepton1.Phi());
              
                    temp_sublead_pT_mu_from_Z_pT.push_back(lepton3.Pt());
                    temp_sublead_pT_mu_from_Z_eta.push_back(lepton3.Eta());
                    temp_sublead_pT_mu_from_Z_RAPIDITY.push_back(lepton3.Rapidity());
                    temp_sublead_pT_mu_from_Z_phi.push_back(lepton3.Phi());
              
                    temp_lead_pT_mu_from_upsi_pT.push_back(lepton2.Pt());
                    temp_lead_pT_mu_from_upsi_eta.push_back(lepton2.Eta());
                    temp_lead_pT_mu_from_upsi_RAPIDITY.push_back(lepton2.Rapidity());
                    temp_lead_pT_mu_from_upsi_phi.push_back(lepton2.Phi());
              
                    temp_sublead_pT_mu_from_upsi_pT.push_back(lepton4.Pt());
                    temp_sublead_pT_mu_from_upsi_eta.push_back(lepton4.Eta());
                    temp_sublead_pT_mu_from_upsi_RAPIDITY.push_back(lepton4.Rapidity());
                    temp_sublead_pT_mu_from_upsi_phi.push_back(lepton4.Phi());
                    
                    temp_upsi_type.push_back(upsi_type);
                  }
                }
              }
            }
         } //Keep me, I am separate from the truth matching stuff you are adding
    
         if (upsi_phase1_first_Z_second_pair_13_24_56) {
           
        
           //Z cuts
            if (lepton2.Pt() < lead_mu_from_Z_pT_Cut  || lepton4.Pt() < sublead_mu_from_Z_pT_Cut){
                std::cout << "FAILED Z mu Pt Cuts" << std::endl;
                continue;
             }  
            
            if (fabs(lepton2.Eta()) > mu_from_Z_eta_Cut || fabs(lepton4.Eta()) > mu_from_Z_eta_Cut){
                std::cout << "FAILED Z mu eta Cuts!" << std::endl; 
                continue;
            }
            
            if (TREE->lepton2_isTightMuon->at(i) + TREE->lepton4_isTightMuon->at(i) != 2){  //both of them need to be tight, tight has a value of 1, 1 +1 =2 
               std::cout << "AT LEAST ONE OF THE MUS FROM A Z WAS NOT TIGHT, FAILED THE Z->MU MU BOTH MU MOST BE TIGHT CUT" << std::endl;
               continue;
           
           } 
           
            if (fabs(TREE->lepton2_impactParameterSignificance->at(i)) > mu_from_Z_3DIPSig_Cut || fabs(TREE->lepton4_impactParameterSignificance->at(i)) > mu_from_Z_3DIPSig_Cut){
               std::cout << "FAILED mu froom Z IP sig cut!" << std::endl;
               continue; 
            }
            
            if (pfIso_lep2 > pfIso_Cut_Mu_from_Z || pfIso_lep4 > pfIso_Cut_Mu_from_Z){
              std::cout << "FAILED particle flow iso cut" << std::endl;
              pfIso_Fail_Count += 1;
              continue;
           }
           
           if (TREE->dimuon2vtx->at(i) < mu_mu_from_Z_Prob_Cut){
               std::cout << "FAILED mu_mu_from_Z_Prob_Cut" << std::endl;
               continue; 
            }

           
           //Upsi cuts
            
            if (lepton1.Pt() < mu_from_upsi_pT_Cut || lepton3.Pt() < mu_from_upsi_pT_Cut){
               std::cout << "FAILED upsi mu pT cuts!" << std::endl; 
               continue;
            }
            
            if (fabs(lepton1.Eta()) > mu_from_upsi_eta_Cut || fabs(lepton3.Eta()) > mu_from_upsi_eta_Cut){
               std::cout << "FAILED upsi mu eta cuts!" << std::endl;
               continue; 
            }
            
            if ( (lepton1 + lepton3).M() < upsi_mass_low_phase2 || (lepton1 + lepton3).M() > upsi_mass_high_phase2){
               std::cout << "FAILED the tighter phase2 upsi mass cuts!" << std::endl;
               continue; 
            }
            
            if (TREE->lepton1_isSoftMuon->at(i) + TREE->lepton3_isSoftMuon->at(i) != 2){
               std::cout << "FAILED mu from upsi must be soft cut" << std::endl;
               continue; 
            }
            
            if //( fabs(lepton1.Rapidity()) > mu_from_upsi_RAPIDITY_Cut || fabs(lepton3.Rapidity()) > mu_from_upsi_RAPIDITY_Cut ){
               
               (fabs ((lepton1+lepton3).Rapidity()) > upsi_RAPIDITY_Cut ){
               std::cout << "FAILED upsi RAPIDITY cut" << std::endl;
               continue; 
            }
            
            if (TREE->dimuon1vtx->at(i) < mu_mu_from_upsi_Prob_Cut){
                std::cout << "FAILED mu_mu_from_upsi_Prob_Cut" << std::endl;
                continue; 
            }
            
            if (pfIso_lep1 > pfIso_Cut_Mu_from_Upsi || pfIso_lep3 > pfIso_Cut_Mu_from_Upsi){
              continue; 
            }
            
   
            if (!doMCTruthMatching){
                temp_Z_mass.push_back((lepton2+lepton4).M());
                temp_upsi_mass.push_back((lepton1+lepton3).M());
              
              //Pt,Eta, Rapidity, Phi of Z, upsi
                temp_Z_pT.push_back((lepton2+lepton4).Pt());
                temp_Z_eta.push_back((lepton2+lepton4).Eta());
                temp_Z_RAPIDITY.push_back((lepton2+lepton4).Rapidity());
                temp_Z_phi.push_back((lepton2 + lepton4).Phi());
              
                temp_upsi_pT.push_back((lepton1+lepton3).Pt());
                temp_upsi_eta.push_back((lepton1+lepton3).Eta());
                temp_upsi_RAPIDITY.push_back((lepton1+lepton3).Rapidity());
                temp_upsi_phi.push_back((lepton1+lepton3).Phi());
              
              //pT, eta, Rapidity, Phi of daughter muons 
                temp_lead_pT_mu_from_Z_pT.push_back(lepton2.Pt());
                temp_lead_pT_mu_from_Z_eta.push_back(lepton2.Eta());
                temp_lead_pT_mu_from_Z_RAPIDITY.push_back(lepton2.Rapidity());
                temp_lead_pT_mu_from_Z_phi.push_back(lepton2.Phi());
              
                temp_sublead_pT_mu_from_Z_pT.push_back(lepton4.Pt());
                temp_sublead_pT_mu_from_Z_eta.push_back(lepton4.Eta());
                temp_sublead_pT_mu_from_Z_RAPIDITY.push_back(lepton4.Rapidity());
                temp_sublead_pT_mu_from_Z_phi.push_back(lepton4.Phi());
              
                temp_lead_pT_mu_from_upsi_pT.push_back(lepton1.Pt());
                temp_lead_pT_mu_from_upsi_eta.push_back(lepton1.Eta());
                temp_lead_pT_mu_from_upsi_RAPIDITY.push_back(lepton1.Rapidity());
                temp_lead_pT_mu_from_upsi_phi.push_back(lepton1.Phi());
              
                temp_sublead_pT_mu_from_upsi_pT.push_back(lepton3.Pt());
                temp_sublead_pT_mu_from_upsi_eta.push_back(lepton3.Eta());
                temp_sublead_pT_mu_from_upsi_RAPIDITY.push_back(lepton3.Rapidity()); 
                temp_sublead_pT_mu_from_upsi_phi.push_back(lepton3.Phi());
              }
            
            if (doMCTruthMatching){
              std::cout << "Poodles 3" << std::endl; 
              
              int entriesMC = (TREEMC->fChain)->GetEntries();
              
              for(int iEntry=0; iEntry<entriesMC; iEntry++) {
                //std::cout << "Koala bear"<< std::endl; 
                (TREEMC->fChain)->GetEntry(iEntry);
                
                if (TREEMC->mc_event_number->at(0) == evNumThisQuad && TREEMC->mc_run_number->at(0) == runNumThisQuad && TREEMC->mc_lumi_section->at(0) == LSThisQuad ){
                  std::cout << "Poodles again!" << std::endl;
                  
                  bool matched = false;
                  
                  bool found2 = false;
                  bool found4 = false;
                  int found2Index = -1;
                  
                  for (int i=0; i<(int)TREEMC->truth_Zmuon_pt->size(); i++){
                    TLorentzVector Z_mu_truth;
                    Z_mu_truth.SetPtEtaPhiM(TREEMC->truth_Zmuon_pt->at(i), TREEMC->truth_Zmuon_eta->at(i), TREEMC->truth_Zmuon_phi->at(i), muon_mass);
                    if (Z_mu_truth.DeltaR(lepton2) < deltaRCut){
                      std::cout << "EUREKA" << std::endl;
                      found2Index = i;
                      found2 = true;
                    }
                    
                    if (Z_mu_truth.DeltaR(lepton4) < deltaRCut && i != found2Index){
                      std::cout << "EUREKA 2" << std::endl;
                      found4 = true;
                    }
                  }
                  
                  if (found2 && found4){
                    std::cout << "Matched Z" << std::endl;
                    matchedZCount++;
                  }
                  
                  bool found1_in_Upsi1 = false;
                  bool found3_in_Upsi1 = false; 
                  int found1_in_Upsi1_Index = -1;
                  
                  for (int j = 0; j < (int)TREEMC->truth_Upsimuon_pt->size(); j++){
                    //std::cout << "Poodle" << std::endl;
                    TLorentzVector Upsi1_mu_truth;
                    Upsi1_mu_truth.SetPtEtaPhiM(TREEMC->truth_Upsimuon_pt->at(j), TREEMC->truth_Upsimuon_eta->at(j), TREEMC->truth_Upsimuon_phi->at(j), muon_mass);
                    if (Upsi1_mu_truth.DeltaR(lepton1) < deltaRCut){
                     // std::cout << "EUREKA 3" << std::endl;
                      found1_in_Upsi1 = true;
                      found1_in_Upsi1_Index = j;
                    }
                    
                    if (Upsi1_mu_truth.DeltaR(lepton3) < deltaRCut && j != found1_in_Upsi1_Index){
                     // std::cout << "EUREKA 4" << std::endl;
                      found3_in_Upsi1 = true;
                    }
                    
                  }
                  
                  bool found1_in_Upsi2 = false;
                  bool found3_in_Upsi2 = false;
                  int found1_in_Upsi2_Index = -1; 
                  
                  for (int k = 0; k < (int)TREEMC->truth_Upsi2muon_pt->size(); k++){
                   // std::cout << "Great Ape" << std::endl;
                    TLorentzVector Upsi2_mu_truth;
                    Upsi2_mu_truth.SetPtEtaPhiM(TREEMC->truth_Upsi2muon_pt->at(k), TREEMC->truth_Upsi2muon_eta->at(k), TREEMC->truth_Upsi2muon_phi->at(k), muon_mass);
                    if (Upsi2_mu_truth.DeltaR(lepton1) < deltaRCut){
                     // std::cout << "found2_in_Upsi2" << std::endl;
                      found1_in_Upsi2 = true;
                      found1_in_Upsi2_Index = k;
                    }
                    
                    if (Upsi2_mu_truth.DeltaR(lepton3) < deltaRCut && k != found1_in_Upsi2_Index){
                     // std::cout << "found4_in_Upsi2" << std::endl;
                      found3_in_Upsi2 = true;
                    }
                  }
                  
                  bool found1_in_Upsi3 = false;
                  bool found3_in_Upsi3 = false;
                  int found1_in_Upsi3_Index = -1;
                  
                  for (int l = 0; l < (int)TREEMC->truth_Upsi3muon_pt->size(); l++){
                   //std::cout << "Ron Swanson" << std::endl; 
                   TLorentzVector Upsi3_mu_truth;
                   Upsi3_mu_truth.SetPtEtaPhiM(TREEMC->truth_Upsi3muon_pt->at(l), TREEMC->truth_Upsi3muon_eta->at(l), TREEMC->truth_Upsi3muon_phi->at(l), muon_mass);
                   if (Upsi3_mu_truth.DeltaR(lepton1) < deltaRCut){
                    // std::cout << "found2_in_Upsi3 " << std::endl;
                     found1_in_Upsi3 = true;
                     found1_in_Upsi3_Index = l;
                   }
                   
                   if (Upsi3_mu_truth.DeltaR(lepton3) < deltaRCut && l != found1_in_Upsi3_Index){
                    // std::cout << "found2_in_Upsi3" << std::endl;
                     found3_in_Upsi3 = true;
                   }
                 } 
                 
                  if (found2 && found4 && found1_in_Upsi1 && found3_in_Upsi1 && !found1_in_Upsi2 && !found3_in_Upsi2 && !found1_in_Upsi3 && !found3_in_Upsi3){
                    upsi_type = 1;
                    matched = true;
                    matchedCount++;
                  }
                  
                  if (found2 && found4 && found1_in_Upsi2 && found3_in_Upsi2 && !found1_in_Upsi1 && !found3_in_Upsi1 && !found1_in_Upsi3 && !found3_in_Upsi3){
                    upsi_type = 2;
                    matched = true;
                    matchedCount++;
                  }
                  
                  //start here tomorrow 7 December 2021
                  if (found2 && found4 && found1_in_Upsi3 && found3_in_Upsi3 && !found1_in_Upsi1 && !found3_in_Upsi1 && !found1_in_Upsi2 && !found3_in_Upsi2){
                    upsi_type = 3;
                    matched = true;
                    matchedCount++;
                  }
                  
                  if (matched){
                   temp_Z_mass.push_back((lepton2+lepton4).M());
                   temp_upsi_mass.push_back((lepton1+lepton3).M());
              
                    //Pt,Eta, Rapidity, Phi of Z, upsi
                   temp_Z_pT.push_back((lepton2+lepton4).Pt());
                   temp_Z_eta.push_back((lepton2+lepton4).Eta());
                   temp_Z_RAPIDITY.push_back((lepton2+lepton4).Rapidity());
                   temp_Z_phi.push_back((lepton2 + lepton4).Phi());
              
                   temp_upsi_pT.push_back((lepton1+lepton3).Pt());
                   temp_upsi_eta.push_back((lepton1+lepton3).Eta());
                   temp_upsi_RAPIDITY.push_back((lepton1+lepton3).Rapidity());
                   temp_upsi_phi.push_back((lepton1+lepton3).Phi());
              
                    //pT, eta, Rapidity, Phi of daughter muons 
                   temp_lead_pT_mu_from_Z_pT.push_back(lepton2.Pt());
                   temp_lead_pT_mu_from_Z_eta.push_back(lepton2.Eta());
                   temp_lead_pT_mu_from_Z_RAPIDITY.push_back(lepton2.Rapidity());
                   temp_lead_pT_mu_from_Z_phi.push_back(lepton2.Phi());
              
                   temp_sublead_pT_mu_from_Z_pT.push_back(lepton4.Pt());
                   temp_sublead_pT_mu_from_Z_eta.push_back(lepton4.Eta());
                   temp_sublead_pT_mu_from_Z_RAPIDITY.push_back(lepton4.Rapidity());
                   temp_sublead_pT_mu_from_Z_phi.push_back(lepton4.Phi());
              
                   temp_lead_pT_mu_from_upsi_pT.push_back(lepton1.Pt());
                   temp_lead_pT_mu_from_upsi_eta.push_back(lepton1.Eta());
                   temp_lead_pT_mu_from_upsi_RAPIDITY.push_back(lepton1.Rapidity());
                   temp_lead_pT_mu_from_upsi_phi.push_back(lepton1.Phi());
              
                   temp_sublead_pT_mu_from_upsi_pT.push_back(lepton3.Pt());
                   temp_sublead_pT_mu_from_upsi_eta.push_back(lepton3.Eta());
                   temp_sublead_pT_mu_from_upsi_RAPIDITY.push_back(lepton3.Rapidity()); 
                   temp_sublead_pT_mu_from_upsi_phi.push_back(lepton3.Phi());
                   
                   temp_upsi_type.push_back(upsi_type);
                  }
                }
              }
            }
         } //Keep me separate from the truth matching stuff you are adding
    
    }
      
      if (TREE->pair_14_23_56->at(i) == 1){
       //  std::cout << "TREE->pair_14_23_56->at(i) == 1" << std::endl; 
        
         bool Z_first_upsi_phase1_second_pair_14_23_56 = false;
         bool upsi_phase1_first_Z_second_pair_14_23_56 = false; 
         pair_14_23_56_count += 1; 
         
         if ( (lepton1 + lepton4).M() > Z_mass_low && (lepton1 + lepton4).M() < Z_mass_high && (lepton2+lepton3).M()  > upsi_mass_low_phase1 && (lepton2+lepton3).M() < upsi_mass_high_phase1){
            if ( (TREE->lepton1_charge->at(i) + TREE->lepton4_charge->at(i) == 0) && (TREE->lepton2_charge->at(i) + TREE->lepton3_charge->at(i) == 0) ) {
              Z_first_upsi_phase1_second_pair_14_23_56 = true;
              Z_first_upsi_phase1_second_pair_14_23_56_count +=1;
              std::cout << "Z_first_upsi_phase1_second_pair_14_23_56 = true!" <<std::endl; 
             }
          }
         
         if ( (lepton1 + lepton4).M() > upsi_mass_low_phase1 && (lepton1 + lepton4).M() < upsi_mass_high_phase1 && (lepton2+lepton3).M()  > Z_mass_low && (lepton2+lepton3).M() < Z_mass_high ){
            if ( (TREE->lepton1_charge->at(i) + TREE->lepton4_charge->at(i) == 0) && (TREE->lepton2_charge->at(i) + TREE->lepton3_charge->at(i) == 0) ) {
              upsi_phase1_first_Z_second_pair_14_23_56 = true;
              upsi_phase1_first_Z_second_pair_14_23_56_count +=1;
              std::cout << "upsi_phase1_first_Z_second_pair_14_23_56 is true!" << std::endl; 
            }
         }
         
         if (Z_first_upsi_phase1_second_pair_14_23_56){
          
            
            if (lepton1.Pt() < lead_mu_from_Z_pT_Cut  || lepton4.Pt() < sublead_mu_from_Z_pT_Cut){
                std::cout << "FAILED Z mu Pt Cuts" << std::endl;
                continue;
             }  
            
            if (fabs(lepton1.Eta()) > mu_from_Z_eta_Cut || fabs(lepton4.Eta()) > mu_from_Z_eta_Cut){
                std::cout << "FAILED Z mu eta Cuts!" << std::endl; 
                continue;
            }
            
            if (TREE->lepton1_isTightMuon->at(i) + TREE->lepton4_isTightMuon->at(i) != 2){  //both of them need to be tight, tight has a value of 1, 1 +1 =2 
               std::cout << "AT LEAST ONE OF THE MUS FROM A Z WAS NOT TIGHT, FAILED THE Z->MU MU BOTH MU MOST BE TIGHT CUT" << std::endl;
               continue;
           
           } 
           
           if (fabs(TREE->lepton1_impactParameterSignificance->at(i)) > mu_from_Z_3DIPSig_Cut || fabs(TREE->lepton4_impactParameterSignificance->at(i)) > mu_from_Z_3DIPSig_Cut){
               std::cout << "FAILED mu froom Z IP sig cut!" << std::endl;
               continue; 
            }
            
            if (pfIso_lep1 > pfIso_Cut_Mu_from_Z || pfIso_lep4 > pfIso_Cut_Mu_from_Z){
              std::cout << "FAILED particle flow iso cut" << std::endl;
              pfIso_Fail_Count += 1;
              continue;
           }
         
           if (TREE->dimuon1vtx->at(i) < mu_mu_from_Z_Prob_Cut){
               std::cout << "FAILED mu_mu_from_Z_Prob_Cut" << std::endl;
               continue; 
           }
           
           
   
            
            //end Z cuts 
            
            //start upsi cuts 
            if (lepton2.Pt() < mu_from_upsi_pT_Cut || lepton3.Pt() < mu_from_upsi_pT_Cut){
               std::cout << "FAILED upsi mu pT cuts" <<std::endl;
               continue;
            }
            
            if (fabs(lepton2.Eta()) > mu_from_upsi_eta_Cut || fabs(lepton3.Eta()) > mu_from_upsi_eta_Cut){
               std::cout << "FAILED upsi mu eta cuts" << std::endl;
               continue;  
            }
            
            if ( (lepton2 + lepton3).M() < upsi_mass_low_phase2 || (lepton2 + lepton3).M() > upsi_mass_high_phase2 ){
               std::cout << "FAILED the tighter phase2 upsi mass cuts!" << std::endl;
               continue; 
            }
            
            if (TREE->lepton2_isSoftMuon->at(i) + TREE->lepton3_isSoftMuon->at(i) != 2){
               std::cout << "FAILED mu from upsi must be soft cut!" << std::endl;
               continue; 
            }
            
            if //( fabs(lepton2.Rapidity()) > mu_from_upsi_RAPIDITY_Cut || fabs(lepton3.Rapidity()) > mu_from_upsi_RAPIDITY_Cut ){
               ( fabs ((lepton2 + lepton3).Rapidity()) > upsi_RAPIDITY_Cut){
               std::cout << "FAILED  upsi RAPIDITY cut" << std::endl;
               continue; 
            }
            
            if (TREE->dimuon2vtx->at(i) < mu_mu_from_upsi_Prob_Cut){
               std::cout << "FAILED mu_mu_from_upsi_Prob_Cut" << std::endl;
               continue; 
           }
           
           if (pfIso_lep2 > pfIso_Cut_Mu_from_Upsi || pfIso_lep3 > pfIso_Cut_Mu_from_Upsi){
             continue; 
           }
            
 
            if (!doMCTruthMatching){
                 std::cout << "TEST AVALANCHE" << std::endl;
                 temp_Z_mass.push_back((lepton1 + lepton4).M());
                 temp_upsi_mass.push_back((lepton2+lepton3).M());
               
               //Pt, Eta, Phi, Rapidity of Z, upsi
                 temp_Z_pT.push_back((lepton1 + lepton4).Pt());
                 temp_Z_eta.push_back((lepton1+lepton4).Eta());
                 temp_Z_RAPIDITY.push_back((lepton1+lepton4).Rapidity());
                 temp_Z_phi.push_back((lepton1+lepton4).Phi());
//                
                 temp_upsi_pT.push_back((lepton2+lepton3).Pt());
                 temp_upsi_eta.push_back((lepton2+lepton3).Eta());
                 temp_upsi_RAPIDITY.push_back((lepton2+lepton3).Rapidity());
                 temp_upsi_phi.push_back((lepton2+lepton3).Phi());
               
               //Pt, Eta, phi, Rapidity of daughter muons
                 temp_lead_pT_mu_from_Z_pT.push_back(lepton1.Pt());
                 temp_lead_pT_mu_from_Z_eta.push_back(lepton1.Eta());
                 temp_lead_pT_mu_from_Z_RAPIDITY.push_back(lepton1.Rapidity());
                 temp_lead_pT_mu_from_Z_phi.push_back(lepton1.Phi());
               
                 temp_sublead_pT_mu_from_Z_pT.push_back(lepton4.Pt());
                 temp_sublead_pT_mu_from_Z_eta.push_back(lepton4.Eta());
                 temp_sublead_pT_mu_from_Z_RAPIDITY.push_back(lepton4.Rapidity());
                 temp_sublead_pT_mu_from_Z_phi.push_back(lepton4.Phi());
               
                 temp_lead_pT_mu_from_upsi_pT.push_back(lepton2.Pt());
                 temp_lead_pT_mu_from_upsi_eta.push_back(lepton2.Eta());
                 temp_lead_pT_mu_from_upsi_RAPIDITY.push_back(lepton2.Rapidity());
                 temp_lead_pT_mu_from_upsi_phi.push_back(lepton2.Phi());
               
                 temp_sublead_pT_mu_from_upsi_pT.push_back(lepton3.Pt());
                 temp_sublead_pT_mu_from_upsi_eta.push_back(lepton3.Eta());
                 temp_sublead_pT_mu_from_upsi_RAPIDITY.push_back(lepton3.Rapidity()); 
                 temp_sublead_pT_mu_from_upsi_phi.push_back(lepton3.Phi());
            }
              
            if (doMCTruthMatching){
              std::cout << "ELEPHANT" << std::endl; 
              
              int entriesMC = (TREEMC->fChain)->GetEntries();
              for(int iEntry=0; iEntry<entriesMC; iEntry++) {
                //std::cout << "Koala bear"<< std::endl; 
                (TREEMC->fChain)->GetEntry(iEntry);
                
                if (TREEMC->mc_event_number->at(0) == evNumThisQuad && TREEMC->mc_run_number->at(0) == runNumThisQuad && TREEMC->mc_lumi_section->at(0) == LSThisQuad ){
                 // std::cout << "PLACEHOLDER" << std::endl; 
                  
                  bool matched = false;
                  
                  bool found1 = false;
                  bool found4 = false;
                  int found1Index = -1;
                  
                  for (int i=0; i<(int)TREEMC->truth_Zmuon_pt->size(); i++){
                    TLorentzVector Z_mu_truth;
                    Z_mu_truth.SetPtEtaPhiM(TREEMC->truth_Zmuon_pt->at(i), TREEMC->truth_Zmuon_eta->at(i), TREEMC->truth_Zmuon_phi->at(i), muon_mass);
                    if (Z_mu_truth.DeltaR(lepton1) < deltaRCut){
                   // std::cout << "EUREKA" << std::endl;
                      found1Index = i;
                      found1 = true;
                    }
                    
                    if (Z_mu_truth.DeltaR(lepton4) < deltaRCut && i != found1Index){
                    //  std::cout << "EUREKA 2" << std::endl;
                      found4 = true;
                    }
                  }
                  
                  if (found1 && found4){
                   // std::cout << "Matched Z" << std::endl;
                    matchedZCount++;
                  }
                  
                  bool found2_in_Upsi1 = false;
                  bool found3_in_Upsi1 = false;
                  int found2_in_Upsi1_Index = -1;
                  
                  for (int j = 0; j < (int)TREEMC->truth_Upsimuon_pt->size(); j++){
                    //std::cout << "Poodle" << std::endl;
                    TLorentzVector Upsi1_mu_truth;
                    Upsi1_mu_truth.SetPtEtaPhiM(TREEMC->truth_Upsimuon_pt->at(j), TREEMC->truth_Upsimuon_eta->at(j), TREEMC->truth_Upsimuon_phi->at(j), muon_mass);
                    if (Upsi1_mu_truth.DeltaR(lepton2) < deltaRCut){
                     // std::cout << "EUREKA 3" << std::endl;
                      found2_in_Upsi1 = true;
                      found2_in_Upsi1_Index = j;
                    }
                    
                    if (Upsi1_mu_truth.DeltaR(lepton3) < deltaRCut && j != found2_in_Upsi1_Index){
                     // std::cout << "EUREKA 4" << std::endl;
                      found3_in_Upsi1 = true;
                    }
                    
                  }
                  
                  bool found2_in_Upsi2 = false;
                  bool found3_in_Upsi2 = false;
                  int found2_in_Upsi2_Index = -1;
                  
                  for (int k = 0; k < (int)TREEMC->truth_Upsi2muon_pt->size(); k++){
                   // std::cout << "Great Ape" << std::endl;
                    TLorentzVector Upsi2_mu_truth;
                    Upsi2_mu_truth.SetPtEtaPhiM(TREEMC->truth_Upsi2muon_pt->at(k), TREEMC->truth_Upsi2muon_eta->at(k), TREEMC->truth_Upsi2muon_phi->at(k), muon_mass);
                    if (Upsi2_mu_truth.DeltaR(lepton2) < deltaRCut){
                     // std::cout << "found2_in_Upsi2" << std::endl;
                      found2_in_Upsi2 = true;
                      found2_in_Upsi2_Index = k;
                    }
                    
                    if (Upsi2_mu_truth.DeltaR(lepton3) < deltaRCut && k != found2_in_Upsi2_Index){
                     // std::cout << "found4_in_Upsi2" << std::endl;
                      found3_in_Upsi2 = true;
                    }
                  }
                  
                  //start here tomorrow 8 Dec. 2021!
                  
                  bool found2_in_Upsi3 = false;
                  bool found3_in_Upsi3 = false;
                  int found2_in_Upsi3_Index = -1;
                  
                  for (int l = 0; l < (int)TREEMC->truth_Upsi3muon_pt->size(); l++){
                   //std::cout << "Ron Swanson" << std::endl; 
                   TLorentzVector Upsi3_mu_truth;
                   Upsi3_mu_truth.SetPtEtaPhiM(TREEMC->truth_Upsi3muon_pt->at(l), TREEMC->truth_Upsi3muon_eta->at(l), TREEMC->truth_Upsi3muon_phi->at(l), muon_mass);
                   if (Upsi3_mu_truth.DeltaR(lepton2) < deltaRCut){
                    // std::cout << "found2_in_Upsi3 " << std::endl;
                     found2_in_Upsi3 = true;
                     found2_in_Upsi3_Index = l;
                   }
                   
                   if (Upsi3_mu_truth.DeltaR(lepton3) < deltaRCut && l != found2_in_Upsi3_Index){
                    // std::cout << "found2_in_Upsi3" << std::endl;
                     found3_in_Upsi3 = true;
                   }
                 }
                  
                  if (found1 && found4 && found2_in_Upsi1 && found3_in_Upsi1 && !found2_in_Upsi2 && !found3_in_Upsi2 && !found2_in_Upsi3 && !found3_in_Upsi3){
                    upsi_type = 1;
                    matched = true;
                    matchedCount++;
                  }
                  
                  if (found1 && found4 && found2_in_Upsi2 && found3_in_Upsi2 && !found2_in_Upsi1 && !found3_in_Upsi1 && !found2_in_Upsi3 && !found3_in_Upsi3){
                    upsi_type = 2;
                    matched = true;
                    matchedCount++;
                  }
                  
                  if (found1 && found4 && found2_in_Upsi3 && found3_in_Upsi3 && !found2_in_Upsi1 && !found3_in_Upsi1 && !found2_in_Upsi2 && !found3_in_Upsi2){
                    upsi_type = 3;
                    matched = true;
                    matchedCount++;
                  }
                  
                  if (matched) {
                    temp_Z_mass.push_back((lepton1 + lepton4).M());
                    temp_upsi_mass.push_back((lepton2+lepton3).M());
               
                    //Pt, Eta, Phi, Rapidity of Z, upsi
                    temp_Z_pT.push_back((lepton1 + lepton4).Pt());
                    temp_Z_eta.push_back((lepton1+lepton4).Eta());
                    temp_Z_RAPIDITY.push_back((lepton1+lepton4).Rapidity());
                    temp_Z_phi.push_back((lepton1+lepton4).Phi());
//                
                    temp_upsi_pT.push_back((lepton2+lepton3).Pt());
                    temp_upsi_eta.push_back((lepton2+lepton3).Eta());
                    temp_upsi_RAPIDITY.push_back((lepton2+lepton3).Rapidity());
                    temp_upsi_phi.push_back((lepton2+lepton3).Phi());
               
                   //Pt, Eta, phi, Rapidity of daughter muons
                    temp_lead_pT_mu_from_Z_pT.push_back(lepton1.Pt());
                    temp_lead_pT_mu_from_Z_eta.push_back(lepton1.Eta());
                    temp_lead_pT_mu_from_Z_RAPIDITY.push_back(lepton1.Rapidity());
                    temp_lead_pT_mu_from_Z_phi.push_back(lepton1.Phi());
               
                    temp_sublead_pT_mu_from_Z_pT.push_back(lepton4.Pt());
                    temp_sublead_pT_mu_from_Z_eta.push_back(lepton4.Eta());
                    temp_sublead_pT_mu_from_Z_RAPIDITY.push_back(lepton4.Rapidity());
                    temp_sublead_pT_mu_from_Z_phi.push_back(lepton4.Phi());
               
                    temp_lead_pT_mu_from_upsi_pT.push_back(lepton2.Pt());
                    temp_lead_pT_mu_from_upsi_eta.push_back(lepton2.Eta());
                    temp_lead_pT_mu_from_upsi_RAPIDITY.push_back(lepton2.Rapidity());
                    temp_lead_pT_mu_from_upsi_phi.push_back(lepton2.Phi());
               
                    temp_sublead_pT_mu_from_upsi_pT.push_back(lepton3.Pt());
                    temp_sublead_pT_mu_from_upsi_eta.push_back(lepton3.Eta());
                    temp_sublead_pT_mu_from_upsi_RAPIDITY.push_back(lepton3.Rapidity()); 
                    temp_sublead_pT_mu_from_upsi_phi.push_back(lepton3.Phi());
                    
                    temp_upsi_type.push_back(upsi_type);
                  
                  
                  }
                  
                  
                }
                
              }
            }
               
               
               
         }
         
         if (upsi_phase1_first_Z_second_pair_14_23_56) {
           // std::cout << "PLACEHOLDER" << std::endl; 
             if (lepton2.Pt() < lead_mu_from_Z_pT_Cut  || lepton3.Pt() < sublead_mu_from_Z_pT_Cut){
                std::cout << "FAILED Z mu Pt Cuts" << std::endl;
                continue;
             }  
             
             if (fabs(lepton2.Eta()) > mu_from_Z_eta_Cut || fabs(lepton3.Eta()) > mu_from_Z_eta_Cut){
                std::cout << "FAILED Z mu eta Cuts!" << std::endl; 
                continue;
            }
             
             if (TREE->lepton2_isTightMuon->at(i) + TREE->lepton3_isTightMuon->at(i) != 2){  //both of them need to be tight, tight has a value of 1, 1 +1 =2 
               std::cout << "AT LEAST ONE OF THE MUS FROM A Z WAS NOT TIGHT, FAILED THE Z->MU MU BOTH MU MOST BE TIGHT CUT" << std::endl;
               continue;
           
           } 
             
             if (fabs(TREE->lepton2_impactParameterSignificance->at(i)) > mu_from_Z_3DIPSig_Cut || fabs(TREE->lepton3_impactParameterSignificance->at(i)) > mu_from_Z_3DIPSig_Cut){
               std::cout << "FAILED mu froom Z IP sig cut!" << std::endl;
               continue; 
            }
            
            if (pfIso_lep2 > pfIso_Cut_Mu_from_Z || pfIso_lep3 > pfIso_Cut_Mu_from_Z){
              std::cout << "FAILED particle flow iso cut" << std::endl;
              pfIso_Fail_Count += 1;
              continue;
            }
            
            if (TREE->dimuon2vtx->at(i) < mu_mu_from_Z_Prob_Cut){
               std::cout << "FAILED mu_mu_from_Z_Prob_Cut" << std::endl;
               continue; 
            } 
          
             
             //end Z cuts
             
             //start upsi cuts 
             
             if (lepton1.Pt() < mu_from_upsi_pT_Cut || lepton4.Pt() < mu_from_upsi_pT_Cut){
                 std::cout << "FAILED mu from  upsi pT cuts!" << std::endl;
                 continue;
             }
             
             if (fabs(lepton1.Eta()) > mu_from_upsi_eta_Cut || fabs(lepton4.Eta()) > mu_from_upsi_eta_Cut){
                 std::cout << "FAILED mu from upsi eta cuts!" << std::endl;
                 continue; 
             }
             
             if ( (lepton1 + lepton4).M() < upsi_mass_low_phase2 || (lepton1+lepton4).M() > upsi_mass_high_phase2 ){
                 std::cout << "FAILED the tighter phase2 upsi mass cuts!" << std::endl; 
                 continue; 
                
             }
             
             if (TREE->lepton1_isSoftMuon->at(i) + TREE->lepton4_isSoftMuon->at(i) != 2){
                std::cout << "FAILED mu from upsi must be soft cut!" << std::endl; 
                continue ; 
             }
             
             if //( fabs(lepton1.Rapidity()) > mu_from_upsi_RAPIDITY_Cut || fabs(lepton4.Rapidity()) > mu_from_upsi_RAPIDITY_Cut ){
                ( fabs ((lepton1 + lepton4).Rapidity()) > upsi_RAPIDITY_Cut) {
                std::cout << "FAILED upsi RAPIDITY cut" << std::endl;
                continue; 
             }
             
             if (TREE->dimuon1vtx->at(i) < mu_mu_from_upsi_Prob_Cut){
                std::cout << "FAILED mu_mu_from_upsi_Prob_Cut" << std::endl;
                continue; 
            }
            
            if (pfIso_lep1 > pfIso_Cut_Mu_from_Upsi || pfIso_lep4 > pfIso_Cut_Mu_from_Upsi){
              continue;
            }
             
          
             if (!doMCTruthMatching){
                  temp_Z_mass.push_back((lepton2+lepton3).M());
                  temp_upsi_mass.push_back((lepton1+lepton4).M());
                
                //pT, eta, phi, rapidity of Z, upsi
                  temp_Z_pT.push_back((lepton2+lepton3).Pt());
                  temp_Z_eta.push_back((lepton2+lepton3).Eta());
                  temp_Z_RAPIDITY.push_back((lepton2+lepton3).Rapidity());
                  temp_Z_phi.push_back((lepton2+lepton3).Phi());
                
                  temp_upsi_pT.push_back((lepton1+lepton4).Pt());
                  temp_upsi_eta.push_back((lepton1+lepton4).Eta());
                  temp_upsi_RAPIDITY.push_back((lepton1+lepton4).Rapidity());
                  temp_upsi_phi.push_back((lepton1+lepton4).Phi());
                
                //Pt, eta, phi, rapidity of daughter muons
                  temp_lead_pT_mu_from_Z_pT.push_back(lepton2.Pt());
                  temp_lead_pT_mu_from_Z_eta.push_back(lepton2.Eta());
                  temp_lead_pT_mu_from_Z_RAPIDITY.push_back(lepton2.Rapidity());
                  temp_lead_pT_mu_from_Z_phi.push_back(lepton2.Phi());
                
                  temp_sublead_pT_mu_from_Z_pT.push_back(lepton3.Pt());
                  temp_sublead_pT_mu_from_Z_eta.push_back(lepton3.Eta());
                  temp_sublead_pT_mu_from_Z_RAPIDITY.push_back(lepton3.Rapidity());
                  temp_sublead_pT_mu_from_Z_phi.push_back(lepton3.Phi());
                 
                  temp_lead_pT_mu_from_upsi_pT.push_back(lepton1.Pt());
                  temp_lead_pT_mu_from_upsi_eta.push_back(lepton1.Eta());
                  temp_lead_pT_mu_from_upsi_RAPIDITY.push_back(lepton1.Rapidity());
                  temp_lead_pT_mu_from_upsi_phi.push_back(lepton1.Phi());
                
                  temp_sublead_pT_mu_from_upsi_pT.push_back(lepton4.Pt());
                  temp_sublead_pT_mu_from_upsi_eta.push_back(lepton4.Eta());
                  temp_sublead_pT_mu_from_upsi_RAPIDITY.push_back(lepton4.Rapidity()); 
                  temp_sublead_pT_mu_from_upsi_phi.push_back(lepton4.Phi());
                }
             
             if (doMCTruthMatching){
               std::cout << "Doing last set of MC Truth Matching!" << std::endl; 
               
               int entriesMC = (TREEMC->fChain)->GetEntries();
               for(int iEntry=0; iEntry<entriesMC; iEntry++) {
                //std::cout << "Koala bear"<< std::endl; 
                 (TREEMC->fChain)->GetEntry(iEntry);
                 if (TREEMC->mc_event_number->at(0) == evNumThisQuad && TREEMC->mc_run_number->at(0) == runNumThisQuad && TREEMC->mc_lumi_section->at(0) == LSThisQuad ){
                  std::cout << "SUPER MARIO" << std::endl;
                  
                  bool matched = false; 
                  
                  bool found2 = false;
                  bool found3 = false;
                  int found2Index = -1;
                  
                  for (int i=0; i<(int)TREEMC->truth_Zmuon_pt->size(); i++){
                    TLorentzVector Z_mu_truth;
                    Z_mu_truth.SetPtEtaPhiM(TREEMC->truth_Zmuon_pt->at(i), TREEMC->truth_Zmuon_eta->at(i), TREEMC->truth_Zmuon_phi->at(i), muon_mass);
                    if (Z_mu_truth.DeltaR(lepton2) < deltaRCut){
                   // std::cout << "EUREKA" << std::endl;
                      found2Index = i;
                      found2 = true;
                    }
                    
                    if (Z_mu_truth.DeltaR(lepton3) < deltaRCut && i != found2Index){
                    //  std::cout << "EUREKA 2" << std::endl;
                      found3 = true;
                    }
                  }
                  
                  if (found2 && found3){
                    matchedZCount++;
                  }
                  
                  bool found1_in_Upsi1 = false;
                  bool found4_in_Upsi1 = false;
                  int found1_in_Upsi1_Index = -1;
                  
                  for (int j = 0; j < (int)TREEMC->truth_Upsimuon_pt->size(); j++){
                    //std::cout << "Poodle" << std::endl;
                    TLorentzVector Upsi1_mu_truth;
                    Upsi1_mu_truth.SetPtEtaPhiM(TREEMC->truth_Upsimuon_pt->at(j), TREEMC->truth_Upsimuon_eta->at(j), TREEMC->truth_Upsimuon_phi->at(j), muon_mass);
                    if (Upsi1_mu_truth.DeltaR(lepton1) < deltaRCut){
                     // std::cout << "EUREKA 3" << std::endl;
                      found1_in_Upsi1 = true;
                      found1_in_Upsi1_Index = j;
                    }
                    
                    if (Upsi1_mu_truth.DeltaR(lepton4) < deltaRCut && j != found1_in_Upsi1_Index){
                     // std::cout << "EUREKA 4" << std::endl;
                      found4_in_Upsi1 = true;
                    }
                    
                  }
                  
                  bool found1_in_Upsi2 = false;
                  bool found4_in_Upsi2 = false;
                  int found1_in_Upsi2_Index = -1;
                  
                  for (int k = 0; k < (int)TREEMC->truth_Upsi2muon_pt->size(); k++){
                   // std::cout << "Great Ape" << std::endl;
                    TLorentzVector Upsi2_mu_truth;
                    Upsi2_mu_truth.SetPtEtaPhiM(TREEMC->truth_Upsi2muon_pt->at(k), TREEMC->truth_Upsi2muon_eta->at(k), TREEMC->truth_Upsi2muon_phi->at(k), muon_mass);
                    if (Upsi2_mu_truth.DeltaR(lepton1) < deltaRCut){
                     // std::cout << "found2_in_Upsi2" << std::endl;
                      found1_in_Upsi2 = true;
                      found1_in_Upsi2_Index = k;
                    }
                    
                    if (Upsi2_mu_truth.DeltaR(lepton4) < deltaRCut && k != found1_in_Upsi2_Index){
                     // std::cout << "found4_in_Upsi2" << std::endl;
                      found4_in_Upsi2 = true;
                    }
                  }
                  
                  bool found1_in_Upsi3 = false;
                  bool found4_in_Upsi3 = false;
                  int found1_in_Upsi3_Index = -1;
                  
                  for (int l = 0; l < (int)TREEMC->truth_Upsi3muon_pt->size(); l++){
                   //std::cout << "Ron Swanson" << std::endl; 
                   TLorentzVector Upsi3_mu_truth;
                   Upsi3_mu_truth.SetPtEtaPhiM(TREEMC->truth_Upsi3muon_pt->at(l), TREEMC->truth_Upsi3muon_eta->at(l), TREEMC->truth_Upsi3muon_phi->at(l), muon_mass);
                   if (Upsi3_mu_truth.DeltaR(lepton1) < deltaRCut){
                    // std::cout << "found2_in_Upsi3 " << std::endl;
                     found1_in_Upsi3 = true;
                     found1_in_Upsi3_Index = l;
                   }
                   
                   if (Upsi3_mu_truth.DeltaR(lepton4) < deltaRCut && l != found1_in_Upsi3_Index){
                    // std::cout << "found2_in_Upsi3" << std::endl;
                     found4_in_Upsi3 = true;
                   }
                 }
                  
                  if (found2 && found3 && found1_in_Upsi1 && found4_in_Upsi1 && !found1_in_Upsi2 && !found4_in_Upsi2 &&  !found1_in_Upsi3 && !found4_in_Upsi3){
                    upsi_type = 1;
                    matched = true;
                    matchedCount++;  
                  }
                  
                  if (found2 && found3 && found1_in_Upsi2 && found4_in_Upsi2 && !found1_in_Upsi1 && !found4_in_Upsi1 && !found1_in_Upsi3 && !found4_in_Upsi3){
                    upsi_type = 2;
                    matched = true;
                    matchedCount++;
                  }
                  
                  if (found2 && found3 && found1_in_Upsi3 && found4_in_Upsi3 && !found1_in_Upsi1 && !found4_in_Upsi1 && !found1_in_Upsi2 && !found4_in_Upsi2){
                    upsi_type = 3;
                    matched = true;
                    matchedCount++;
                  }
                  
                  if (matched){
                    temp_Z_mass.push_back((lepton2+lepton3).M());
                    temp_upsi_mass.push_back((lepton1+lepton4).M());
                
                   //pT, eta, phi, rapidity of Z, upsi
                    temp_Z_pT.push_back((lepton2+lepton3).Pt());
                    temp_Z_eta.push_back((lepton2+lepton3).Eta());
                    temp_Z_RAPIDITY.push_back((lepton2+lepton3).Rapidity());
                    temp_Z_phi.push_back((lepton2+lepton3).Phi());
                
                    temp_upsi_pT.push_back((lepton1+lepton4).Pt());
                    temp_upsi_eta.push_back((lepton1+lepton4).Eta());
                    temp_upsi_RAPIDITY.push_back((lepton1+lepton4).Rapidity());
                    temp_upsi_phi.push_back((lepton1+lepton4).Phi());
                
                    //Pt, eta, phi, rapidity of daughter muons
                    temp_lead_pT_mu_from_Z_pT.push_back(lepton2.Pt());
                    temp_lead_pT_mu_from_Z_eta.push_back(lepton2.Eta());
                    temp_lead_pT_mu_from_Z_RAPIDITY.push_back(lepton2.Rapidity());
                    temp_lead_pT_mu_from_Z_phi.push_back(lepton2.Phi());
                
                    temp_sublead_pT_mu_from_Z_pT.push_back(lepton3.Pt());
                    temp_sublead_pT_mu_from_Z_eta.push_back(lepton3.Eta());
                    temp_sublead_pT_mu_from_Z_RAPIDITY.push_back(lepton3.Rapidity());
                    temp_sublead_pT_mu_from_Z_phi.push_back(lepton3.Phi());
                 
                    temp_lead_pT_mu_from_upsi_pT.push_back(lepton1.Pt());
                    temp_lead_pT_mu_from_upsi_eta.push_back(lepton1.Eta());
                    temp_lead_pT_mu_from_upsi_RAPIDITY.push_back(lepton1.Rapidity());
                    temp_lead_pT_mu_from_upsi_phi.push_back(lepton1.Phi());
                
                    temp_sublead_pT_mu_from_upsi_pT.push_back(lepton4.Pt());
                    temp_sublead_pT_mu_from_upsi_eta.push_back(lepton4.Eta());
                    temp_sublead_pT_mu_from_upsi_RAPIDITY.push_back(lepton4.Rapidity()); 
                    temp_sublead_pT_mu_from_upsi_phi.push_back(lepton4.Phi());
                    
                    temp_upsi_type.push_back(upsi_type);
                  }
                  
                  
                  
                  
                  
                  
                  
                  }
                
                 
               }
             }
         } //Keep me separate from the doMCTruthMatching stuff you are adding!
      }
    
    
 // Checking that we have the right syntax here      
 //     std::cout << "lepton1.Pt(): " << lepton1.Pt() << std::endl;
      
    

    //Deal with this later

      // the lines below pick up the highest in pt candidate but makes
      // the assumption that lepton 1 2 are from Z and lepton 3 4 from
      // the Upsilon. Which can be wrong.
       //QUESTION: how would you suggest we improve on this?
      //QUESTION: unless I'm mistaken here, looks like we still haven't dealt with the matching ambiguous cases, I need to review the phase 1 code and remind myself of how it works though to be sure 
      if (temp_comparison_pt_upsilon < (lepton3+lepton4).Pt()) {
        mass1_quickAndDirty = (lepton3+lepton4).M();
        temp_comparison_pt_upsilon = (lepton3+lepton4).Pt();
      }
      if (temp_comparison_pt_z < (lepton1+lepton2).Pt()) {
        mass2_quickAndDirty = (lepton1+lepton2).M();
        temp_comparison_pt_z = (lepton1+lepton2).Pt();
      }

    } // loop over the size of the leptons


//This is dirty FIX ME 
//    if (mass1_quickAndDirty > 0. && mass2_quickAndDirty > 0.) //assuming we have found good candidates for both mass1 (upsi) and mass2 (Z), book them. The "good" candidates are taken to be the highest pT ones 
//      aux->Fill();
  
   gotToEndCount += 1; 
   if (temp_Z_mass.size() > 1) {
      std::cout << "FOUND AN EVENT WITH MORE THAN ONE CANDIDATE, THROW IT AWAY! FAILED" << std::endl; 
      QuickCheckCount += 1;
      FailureCount += 1; 
//      continue;
  
   }  
 
  //the final survivor, one per event at most 
   if (temp_Z_mass.size() == 1  && temp_upsi_mass.size() == 1){
     fillCount += 1; 
     
     Z_mass =  temp_Z_mass.at(0);
     upsi_mass = temp_upsi_mass.at(0);
     
     Z_pT = temp_Z_pT.at(0);
     Z_eta = temp_Z_eta.at(0);
     Z_RAPIDITY = temp_Z_RAPIDITY.at(0);
     Z_phi = temp_Z_phi.at(0); 
     upsi_pT = temp_upsi_pT.at(0);
     upsi_eta = temp_upsi_eta.at(0);
     upsi_RAPIDITY = temp_upsi_RAPIDITY.at(0);
     upsi_phi =temp_upsi_phi.at(0);
  
     lead_pT_mu_from_Z_pT = temp_lead_pT_mu_from_Z_pT.at(0);
     lead_pT_mu_from_Z_eta = temp_lead_pT_mu_from_Z_eta.at(0);
     lead_pT_mu_from_Z_RAPIDITY = temp_lead_pT_mu_from_Z_RAPIDITY.at(0);
     lead_pT_mu_from_Z_phi = temp_lead_pT_mu_from_Z_phi.at(0); 
     sublead_pT_mu_from_Z_pT = temp_sublead_pT_mu_from_Z_pT.at(0);
     sublead_pT_mu_from_Z_eta = temp_sublead_pT_mu_from_Z_eta.at(0);
     sublead_pT_mu_from_Z_RAPIDITY =temp_sublead_pT_mu_from_Z_RAPIDITY.at(0);
     sublead_pT_mu_from_Z_phi =temp_sublead_pT_mu_from_Z_phi.at(0); 
 
     lead_pT_mu_from_upsi_pT = temp_lead_pT_mu_from_upsi_pT.at(0);
     lead_pT_mu_from_upsi_eta = temp_lead_pT_mu_from_upsi_eta.at(0);
     lead_pT_mu_from_upsi_RAPIDITY = temp_lead_pT_mu_from_upsi_RAPIDITY.at(0);
     lead_pT_mu_from_upsi_phi = temp_lead_pT_mu_from_upsi_phi.at(0); 
     sublead_pT_mu_from_upsi_pT = temp_sublead_pT_mu_from_upsi_pT.at(0);
     sublead_pT_mu_from_upsi_eta = temp_sublead_pT_mu_from_upsi_eta.at(0);
     sublead_pT_mu_from_upsi_RAPIDITY = temp_sublead_pT_mu_from_upsi_RAPIDITY.at(0); 
     sublead_pT_mu_from_upsi_phi =  temp_sublead_pT_mu_from_upsi_phi.at(0);
     
     if (doMCTruthMatching){
       upsi_type = temp_upsi_type.at(0); 
      }
      
      if (!doMCTruthMatching){
        upsi_type = -1;
      }
     aux->Fill();
    }
  //    fillCount += 1;
  //    aux->Fill();
    
    
  } // loop over the entries


std::cout << "pair_12_34_56_count: " << pair_12_34_56_count << std::endl; 
std::cout << "pair_13_24_56_count: " << pair_13_24_56_count << std::endl;
std::cout << "pair_14_23_56_count: " << pair_14_23_56_count << std::endl; 
std::cout << "pair_AMBIGUOUS_muQuad_count: " << pair_AMBIGUOUS_muQuad_count << std::endl;
std::cout << "big4MuVtx_Prob_Cut_fail_count: " << big4MuVtx_Prob_Cut_fail_count << std::endl;
std::cout << "pfIso_Fail_Count:  " << pfIso_Fail_Count << std::endl; 
std::cout << "Z_first_upsi_phase1_second_pair_12_34_56_count: " << Z_first_upsi_phase1_second_pair_12_34_56_count << std::endl;
std::cout << "upsi_phase1_first_Z_second_pair_12_34_56_count: " << upsi_phase1_first_Z_second_pair_12_34_56_count << std::endl; 
std::cout << "Z_first_upsi_phase1_second_pair_13_24_56_count: " << Z_first_upsi_phase1_second_pair_13_24_56_count << std::endl;
std::cout << "upsi_phase1_first_Z_second_pair_13_24_56_count: " << upsi_phase1_first_Z_second_pair_13_24_56_count << std::endl; 
std::cout << "Z_first_upsi_phase1_second_pair_14_23_56_count: " << Z_first_upsi_phase1_second_pair_14_23_56_count << std::endl;
std::cout << "upsi_phase1_first_Z_second_pair_14_23_56_count: " << upsi_phase1_first_Z_second_pair_14_23_56_count << std::endl; 

std::cout << "GotHereCount_Z_first_upsi_phase1_second_pair_12_34_56_Z_first_upsi_phase1_second_pair_12_34_56:  " << GotHereCount_Z_first_upsi_phase1_second_pair_12_34_56 << std::endl; 
//std::cout << "FailureCount:  " << FailureCount << std::endl;
//std::cout << "QuickCheckCount:  " << QuickCheckCount << std::endl;
std::cout << "fillCount:  " << fillCount << std::endl; 
std::cout << "gotToEndCount:  " << gotToEndCount << std::endl; 
//std::cout << "poodleCount:  " << poodleCount << std::endl; 
std::cout << "eventCounter:  " << eventCounter << std::endl;
std::cout << "mcSanityCheckCount: " << mcSanityCheckCount << std::endl; //this should be GotHereCount_Z_first_upsi_phase1_second_pair_12_34_56 times number of events in the treemc //sanity check passed
std::cout << "matchedZCount: " << matchedZCount << std::endl; 
std::cout << "matchedCount: " << matchedCount << std::endl; 
///////////////////////
//////    M C    //////
///////////////////////

  // double plot_normalization_Z = 0.;
//   double plot_normalization_Upsi = 0.;
//   int entriesMC = (TREEMC->fChain)->GetEntries();
//   for(int iEntry=0; iEntry<entriesMC; iEntry++) {
//     (TREEMC->fChain)->GetEntry(iEntry);
// 
//     for (int i=0; i<(int)TREEMC->truth_Z_mass->size(); i++) {
//       h_truth_Z_mass->Fill(TREEMC->truth_Z_mass->at(i));
//       plot_normalization_Z++;
//     }
// 
//     for (int i=0; i<(int)TREEMC->truth_Jpsi_mass->size(); i++) {
//       h_truth_Upsi_mass->Fill(TREEMC->truth_Jpsi_mass->at(i));
//       plot_normalization_Upsi++;
//     }
// 
//   }
// 
//   h_truth_Z_mass->Scale(h_reco_Z_mass->GetEntries() / plot_normalization_Z);
//   h_truth_Upsi_mass->Scale(h_reco_Upsi_mass->GetEntries() / plot_normalization_Upsi);

/////////////////////////////////////////////////////////
////////////////     P L O T T I N G     ////////////////
/////////////////////////////////////////////////////////

  TCanvas *c_masses = new TCanvas("c_masses", "c_masses", 1000, 500); c_masses->Divide(2,1);
  c_masses->cd(1); h_reco_Upsi_mass_noNewCuts->Draw("e1"); //h_truth_Upsi_mass->Draw("hesame"); //ignoring MC for the moment
  c_masses->cd(2); h_reco_Z_mass_noNewCuts->Draw("e1"); //h_truth_Z_mass->Draw("hesame"); //ignoring MC for the moment 
  h_reco_Upsi_mass_noNewCuts->Write();
  h_reco_Z_mass_noNewCuts->Write();
  c_masses->SaveAs("c_masses.pdf"); //want to save the canvas here because we split it and made it look nice 
  
  
  TCanvas *c_big4MuVtxProb_before_big4MuVtx_Prob_Cut = new TCanvas("c_big4MuVtxProb_before_big4MuVtx_Prob_Cut","c_big4MuVtxProb_before_big4MuVtx_Prob_Cut"); //last 2 are width and height
  c_big4MuVtxProb_before_big4MuVtx_Prob_Cut->cd(); h_big4MuVtxProb_before_big4MuVtx_Prob_Cut->Draw();
  h_big4MuVtxProb_before_big4MuVtx_Prob_Cut->Write();
  h_big4MuVtxProb_before_big4MuVtx_Prob_Cut->SaveAs("h_big4MuVtxProb_before_big4MuVtx_Prob_Cut.pdf");
  
  TCanvas *c_ambig_quad_count = new TCanvas("c_ambig_quad_count", "c_ambig_quad_count");
  c_ambig_quad_count->cd();
  h_ambig_quad->Draw();
  h_ambig_quad->Write();
  c_ambig_quad_count->SaveAs("c_ambig_quad_count.pdf");
  
//  TCanvas *c_dimuon_vtx = new TCanvas("c_dimuon_vtx", "c_dimuon_vtx", 1000, 500); c_dimuon_vtx->Divide(2,1); //the numbers in the canvas declaration are width then height 
//  c_dimuon_vtx->cd(1); h_dimuon_from_Z_Prob_before_Cut->Draw("e1");
//  h_dimuon_from_Z_Prob_before_Cut->Write();
 // c_dimuon_vtx->SaveAs("c_dimuon_vtx.pdf");

  TCanvas *c_cutflow_allQuadCuts = new TCanvas("c_cutflow_allQuadCuts", "c_cutflow_allQuadCuts");
  c_cutflow_allQuadCuts->cd();
  h_cutflow_allQuadCuts->Draw();
  h_cutflow_allQuadCuts->Write();
  c_cutflow_allQuadCuts->SaveAs("h_cutflow_allQuadCuts.pdf");
  
  TCanvas *c_cutflow_Z_first_upsi_phase1_second_pair_12_34_56 = new TCanvas("c_cutflow_Z_first_upsi_phase1_second_pair_12_34_56", "c_cutflow_Z_first_upsi_phase1_second_pair_12_34_56");
  c_cutflow_Z_first_upsi_phase1_second_pair_12_34_56->cd();
  h_cutflow_Z_first_upsi_phase1_second_pair_12_34_56->Draw();
  h_cutflow_Z_first_upsi_phase1_second_pair_12_34_56->Write();
  c_cutflow_Z_first_upsi_phase1_second_pair_12_34_56->SaveAs("h_cutflow_Z_first_upsi_phase1_second_pair_12_34_56.pdf");
  
  TCanvas *c_pfIso_lepN = new TCanvas("c_pfIso_lepN", "c_pfIso_lepN"); c_pfIso_lepN->Divide(2,2);
  c_pfIso_lepN->cd(1); h_pfIso_lep1->Draw();
  c_pfIso_lepN->cd(2); h_pfIso_lep2->Draw();
  c_pfIso_lepN->cd(3); h_pfIso_lep3->Draw();
  c_pfIso_lepN->cd(4); h_pfIso_lep4->Draw();
  h_pfIso_lep1->Write(); h_pfIso_lep2->Write(); h_pfIso_lep3->Write(); h_pfIso_lep4->Write();
  c_pfIso_lepN->SaveAs("c_pfIso_lepN.pdf");
  
 


ntuple->Write();
ntuple->Close();

}
