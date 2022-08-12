//How to run
// root -l -b
// .L loop_Zto4Mu.C++
// run(<fileName.root>)




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

void run(string file){

  // l o a d   t h e   t r e e s 
  root_file = new TFile(file.c_str(),"READ");
  TREE   = new tree((TTree*)root_file->Get("tree"));
  TREEMC = new treeMC((TTree*)root_file->Get("treemc"));
  
  //Announce what root_file is, thanks to Maxi for showing me how to do this
  std::cout << "////////////////////////////////////////" << std::endl;
  std::cout << "Processing file:  " << file.c_str() << std::endl;
  std::cout << "////////////////////////////////////////" << std::endl;

  //Histograms
  TH1F *h_cutflow_allQuadCuts = new TH1F("h_cutflow_allQuadCuts", "h_cutflow_allQuadCuts", 15, -0.5, 14.5); h_cutflow_allQuadCuts->SetXTitle("Cuts involving overall quad");
  h_cutflow_allQuadCuts->Sumw2();
  
  TH1F *h_pfIso_lep1 = new TH1F("h_pfIso_lep1", "h_pfIso_lep1", 60, 0, 3); h_pfIso_lep1->SetXTitle("PF Isolation for lep1");
  h_pfIso_lep1->Sumw2();
  
  TH1F *h_pfIso_lep2 = new TH1F("h_pfIso_lep2", "h_pfIso_lep2", 60, 0, 3); h_pfIso_lep2->SetXTitle("PF Isolation for lep2");
  h_pfIso_lep2->Sumw2();
  
  TH1F *h_pfIso_lep3 = new TH1F("h_pfIso_lep3", "h_pfIso_lep3", 60, 0, 3); h_pfIso_lep3->SetXTitle("PF Isolation for lep3");
  h_pfIso_lep3->Sumw2();
  
  TH1F *h_pfIso_lep4 = new TH1F("h_pfIso_lep4", "h_pfIso_lep4", 60, 0, 3); h_pfIso_lep4->SetXTitle("PF Isolation for lep4");
  h_pfIso_lep4 ->Sumw2();
  
  TH1F *h_pair_sum = new TH1F("h_pair_sum" , "h_pair_sum", 5, -0.5, 4.5); h_pair_sum->SetXTitle("Sum of pair_12_34_ZOnly, pair_13_24_ZOnly, & pair_14_23_ZOnly");
  h_pair_sum->Sumw2();
  
  TH1F *h_cart_DR_Sig = new TH1F("h_cart_DR_Sig", "h_cart_DR_Sig", 50, 0, 10); h_cart_DR_Sig->SetXTitle("Cartesian DR Significance");
  h_cart_DR_Sig->Sumw2();
  
  TH1F *h_DZ_Sig = new TH1F("h_DZ_Sig", "h_DZ_Sig", 50, 0, 10); h_DZ_Sig->SetXTitle("DZ Significance");
  h_DZ_Sig->Sumw2();
  
  // constants
  double muon_mass = 105.6583 / 1000.; //get mass in GeV
  
  //Counters
  int is_pair_12_34_ZOnly_Count = 0;
  int is_pair_13_24_ZOnly_Count = 0;
  int is_pair_14_23_ZOnly_Count = 0;
  int fillCount = 0; 
   
   //non-boolean flags

 int triggerYear = 2016; //options are 2016, 2017, 2018
//  int triggerYear = 2017;
//  int triggerYear = 2018;
  
  std::cout << "Using triggers for year:  " << triggerYear << std::endl;
  std::cout << "//////////////////////" << std::endl;
   
   //boolean flags
   
   bool doRecoToTrigMuMatching = false;
   
   if (!doRecoToTrigMuMatching){
     std::cout << "NOT including the reco muon to trigger muon matching cut" << std::endl;
     std::cout << "///////////////////////////////////////////////////////" << std::endl;
   }
   
   if (doRecoToTrigMuMatching){
     std::cout << "INCLUDING the reco muon to trigger muon matching cut" << std::endl;
     std::cout << "////////////////////////////////////////////////////" << std::endl;
   }
   
  //Cuts
  
  double pfIso_Cut = 0.35;
  
  double lepSeparationCut = 0.02;
  //double lepSeparationCut = 0.1; //this was for testing purposes, do not actually use, 
  //the value of lepSeparationCut to use from the Z->4 lepton AN is 0.02
  
  double etaCut = 2.4;
  double lepton1_pT_Cut = 20; //ordering is done in pT, i.e. lepton1 is the lepton with the leading pT, lepton2 with the subleading, etc
  double lepton2_pT_Cut = 10;
  double lepton3_pT_Cut = 5;
  double lepton4_pT_Cut = 5;
  
  double lep_3DIPSig_Cut = 4;
  double lep_dxy_Cut = 0.5;
  double lep_dz_Cut = 1;
  
  double deltaR_dimuon1vtx_dimuon2vtx_Cut = 3; //value suggested by Greg in email titled "middle of Feb. 2022 work"
  
  double offset = 0.4; //value suggested by Greg in email titled "middle of Feb. 2022 Work"
  
  double dimuonvtx_Prob_Cut = 0.05;
  
  //New skimmed root file
  double Z_mass = -99;
  double Z_pT = -99;
  double Z_eta = -99;
  double Z_phi = -99;
  
  double Z1_mass = -99;
  double Z1_pT = -99;
  double Z1_eta = -99;
  double Z1_phi = -99;
  
  double Z2_mass = -99;
  double Z2_pT = -99;
  double Z2_eta = -99;
  double Z2_phi = -99;
  
  double lepton1_pT = -99;
  double lepton2_pT = -99;
  double lepton3_pT = -99;
  double lepton4_pT = -99;
  
  double lepton1_eta = -99;
  double lepton2_eta = -99;
  double lepton3_eta = -99;
  double lepton4_eta = -99;
  
  double lepton1_phi = -99;
  double lepton2_phi = -99;
  double lepton3_phi = -99;
  double lepton4_phi = -99;
  
  TFile *ntuple = new TFile("clean_2August2022_loop_Zto4Mu_inputFileIs_12July2022_Run2016_Total_noTrigToRecoMuMatching_JanCutsRemoved.root", "RECREATE");
  TTree *aux;
  aux = new TTree("tree", "tree");
  
  //Z branches
  aux->Branch("Z_mass", &Z_mass);
  aux->Branch("Z_pT", &Z_pT);
  aux->Branch("Z_eta", &Z_eta);
  aux->Branch("Z_phi", &Z_phi);
  
  //branches for Z1, the heavier of the two pairs
  aux->Branch("Z1_mass", &Z1_mass);
  aux->Branch("Z1_pT", &Z1_pT);
  aux->Branch("Z1_eta", &Z1_eta);
  aux->Branch("Z1_phi", &Z1_phi);
  
  //branches for Z2, the lighter of the two pairs
  aux->Branch("Z2_mass", &Z2_mass);
  aux->Branch("Z2_pT", &Z2_pT);
  aux->Branch("Z2_phi", &Z2_phi);
  aux->Branch("Z2_eta", &Z2_eta);
  
  //branches for leptons, all of which are assumed to be muons with mass muon_mass
  
  //branches for lepton1 (which is assumed to be a muon)
  aux->Branch("lepton1_pT", &lepton1_pT);
  aux->Branch("lepton1_eta", &lepton1_eta);
  aux->Branch("lepton1_phi", &lepton1_phi);
  
  //branches for lepton2 (which is assumed to be a muon)
  aux->Branch("lepton2_pT", &lepton2_pT);
  aux->Branch("lepton2_eta", &lepton2_eta);
  aux->Branch("lepton2_phi", &lepton2_phi);
  
  //branches for lepton3 (which is assumed to be a muon)
  aux->Branch("lepton3_pT", &lepton3_pT);
  aux->Branch("lepton3_eta", &lepton3_eta);
  aux->Branch("lepton3_phi", &lepton3_phi);
  
  //branches for lepton4 (which is assumed to be a muon)
  aux->Branch("lepton4_pT", &lepton4_pT);
  aux->Branch("lepton4_eta", &lepton4_eta);
  aux->Branch("lepton4_phi", &lepton4_phi);
  
  
  int eventCounter = 0;
  
  int entries = (TREE->fChain)->GetEntries();
  std::cout << "number of entries:  " << entries << std::endl; 
  for(int iEntry=0; iEntry<entries; iEntry++) {
    (TREE->fChain)->GetEntry(iEntry);
     eventCounter += 1;
    if (eventCounter % 1000 == 0){
      std::cout << "Processed  " << eventCounter << "  Events" << std::endl; 
    }    
    
    std::vector<double> temp_Z_mass;
    std::vector<double> temp_Z_pT;
    std::vector<double> temp_Z_eta;
    std::vector<double> temp_Z_phi;
    
    std::vector<double> temp_Z1_mass;
    std::vector<double> temp_Z1_pT;
    std::vector<double> temp_Z1_eta;
    std::vector<double> temp_Z1_phi;
    
    std::vector<double> temp_Z2_mass;
    std::vector<double> temp_Z2_pT;
    std::vector<double> temp_Z2_eta;
    std::vector<double> temp_Z2_phi;
    
    std::vector<double> temp_lep1_pT;
    std::vector<double> temp_lep1_eta;
    std::vector<double> temp_lep1_phi;
    
    std::vector<double> temp_lep2_pT;
    std::vector<double> temp_lep2_eta;
    std::vector<double> temp_lep2_phi;
    
    std::vector<double> temp_lep3_pT;
    std::vector<double> temp_lep3_eta;
    std::vector<double> temp_lep3_phi;
    
    std::vector<double> temp_lep4_pT;
    std::vector<double> temp_lep4_eta;
    std::vector<double> temp_lep4_phi;
    
    
    temp_Z_mass.clear();
    temp_Z_pT.clear();
    temp_Z_eta.clear();
    temp_Z_phi.clear();
    
    temp_Z1_mass.clear();
    temp_Z1_pT.clear();
    temp_Z1_eta.clear();
    temp_Z1_phi.clear();
    
    temp_Z2_mass.clear();
    temp_Z2_pT.clear();
    temp_Z2_eta.clear();
    temp_Z2_phi.clear();
    
    temp_lep1_pT.clear();
    temp_lep1_eta.clear();
    temp_lep1_phi.clear();
    
    temp_lep2_pT.clear();
    temp_lep2_eta.clear();
    temp_lep2_phi.clear();
    
    temp_lep3_pT.clear();
    temp_lep3_eta.clear();
    temp_lep3_phi.clear();
    
    temp_lep4_pT.clear();
    temp_lep4_eta.clear();
    temp_lep4_phi.clear();
   
    //handle the trigger
   
    bool event_fails_trigger = true; //defaults to true
    
    //bools for 2016 triggers
    bool singleMu2016Trig1Fired = false;
    bool singleMu2016Trig2Fired = false;
    bool doubleMu2016Trig1Fired = false;
    bool doubleMu2016Trig2Fired = false;
    bool tripleMu2016Trig1Fired = false;
   
   //bools for 2017 triggers
    bool singleMu2017Trig1Fired = false;
    bool doubleMu2017Trig1Fired = false; 
    bool tripleMu2017Trig1Fired = false;
    bool tripleMu2017Trig2Fired = false; 
   
    //bools for 2018 triggers
    bool singleMu2018Trig1Fired = false;
    bool doubleMu2018Trig1Fired = false; 
    bool doubleMu2018Trig2Fired = false; 
    bool tripleMu2018Trig1Fired = false;
    bool tripleMu2018Trig2Fired = false; 
     
    //event level loop on the contents of the triggerlist
    for (int iTrig =0; iTrig < (int)TREE->triggerlist->size(); iTrig++){
     //std::cout << TREE->triggerlist->at(iTrig) << std::endl;
      std::string str (TREE->triggerlist->at(iTrig));
     //std::cout << "str:  " << str << std::endl; 
    
      if (triggerYear == 2016){
        std::string str2 ("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v"); //call this DoubleMu2016Trig1
        std::string str3 ("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v"); //call this DoubleMu2016Trig2
        std::string str4 ("HLT_IsoMu24_v"); //call this SingleMu2016Trig1
        std::string str5 ("HLT_IsoTkMu24_v"); //call this SingleMu2016Trig2
        std::string str5a ("HLT_TripleMu_12_10_5_v"); //call this TripleMu2016Trig1
         
        std::size_t foundDoubleMu2016Trig1 = str.find(str2);
        std::size_t foundDoubleMu2016Trig2 = str.find(str3);
        std::size_t foundSingleMu2016Trig1 = str.find(str4);
        std::size_t foundSingleMu2016Trig2 = str.find(str5);
        std::size_t foundTripleMu2016Trig1 = str.find(str5a);
        
        if (foundSingleMu2016Trig1 != std::string::npos){
          singleMu2016Trig1Fired = true;
         //std::cout << "singleMu2016Trig1Fired:  "  << singleMu2016Trig1Fired << std::endl; 
        }
       
        if (foundSingleMu2016Trig2 != std::string::npos){
          singleMu2016Trig2Fired = true;
        //std::cout << "singleMu2016Trig2Fired:  " << singleMu2016Trig2Fired << std::endl; 
        }
        
        if (foundDoubleMu2016Trig1 != std::string::npos){
          doubleMu2016Trig1Fired = true;
          // std::cout << "doubleMu2016Trig1Fired:  " << doubleMu2016Trig1Fired << std::endl; 
        }
        
        if (foundDoubleMu2016Trig2 != std::string::npos){
          doubleMu2016Trig2Fired = true;
        //std::cout << "doubleMu2016Trig2Fired:  " << doubleMu2016Trig2Fired << std::endl; 
        }
        
        if (foundTripleMu2016Trig1 != std::string::npos){
          tripleMu2016Trig1Fired = true;
        //   std::cout << "tripleMu2016Trig1Fired:  " << tripleMu2016Trig1Fired << std::endl; 
        }
        
        if (singleMu2016Trig1Fired || singleMu2016Trig2Fired || doubleMu2016Trig1Fired || doubleMu2016Trig2Fired || tripleMu2016Trig1Fired){
          event_fails_trigger = false;
         // std::cout << "Event passed 2016 triggers" << std::endl;
          break; 
        }
      
      }
      
      if (triggerYear == 2017){
        std::string str6 ("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v"); //call this DoubleMu2017Trig1
        std::string str7 ("HLT_IsoMu27_v"); //call this SingleMu2017Trig1
        std::string str7a ("HLT_TripleMu_12_10_5_v"); //call this TripleMu2017Trig1
        std::string str7b ("HLT_TripleMu_10_5_5_DZ_v"); //call this TripleMu2017Trig2
        
        std::size_t foundDoubleMu2017Trig1 = str.find(str6);
        std::size_t foundSingleMu2017Trig1 = str.find(str7);
        std::size_t foundTripleMu2017Trig1 = str.find(str7a);
        std::size_t foundTripleMu2017Trig2 = str.find(str7b);
      
        if (foundSingleMu2017Trig1 != std::string::npos){
          singleMu2017Trig1Fired = true;
         //  std::cout << "singleMu2017Trig1Fired:  " << singleMu2017Trig1Fired << std::endl; 
        } 
        
        if (foundDoubleMu2017Trig1 != std::string::npos){
          doubleMu2017Trig1Fired = true;
           //std::cout << "doubleMu2017Trig1Fired:  " << doubleMu2017Trig1Fired << std::endl; 
        }
        
        if (foundTripleMu2017Trig1 != std::string::npos){
          tripleMu2017Trig1Fired = true;
         //  std::cout << "tripleMu2017Trig1Fired:  " << tripleMu2017Trig1Fired << std::endl;
        }
        
        if (foundTripleMu2017Trig2 != std::string::npos){
          tripleMu2017Trig2Fired = true;
         //  std::cout << "tripleMu2017Trig2Fired:  " << tripleMu2017Trig2Fired << std::endl;
        }
        
        if (singleMu2017Trig1Fired || doubleMu2017Trig1Fired || tripleMu2017Trig1Fired || tripleMu2017Trig2Fired){
          event_fails_trigger = false;
         // std::cout << "Event passed 2017 triggers" << std::endl; 
          break;
        } 
      }
      
      if (triggerYear == 2018){
        std::string str8 ("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v"); //call this DoubleMu2018Trig1
        std::string str8a ("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v"); // call this DoubleMu2018Trig2
        std::string str9 ("HLT_IsoMu24_v"); // call this SingleMu2018Trig1
        std::string str9a ("HLT_TripleMu_12_10_5_v"); // call this TripleMu2018Trig1
        std::string str9b ("HLT_TripleMu_10_5_5_DZ_v"); // call this TripleMu2018Trig2
        
        std::size_t foundDoubleMu2018Trig1 = str.find(str8);
        std::size_t foundDoubleMu2018Trig2 = str.find(str8a);
        std::size_t foundSingleMu2018Trig1 = str.find(str9);
        std::size_t foundTripleMu2018Trig1 = str.find(str9a);
        std::size_t foundTripleMu2018Trig2 = str.find(str9b);
        
        if (foundSingleMu2018Trig1 != std::string::npos){
          singleMu2018Trig1Fired = true;
          // std::cout << "singleMu2018Trig1Fired:  " << singleMu2018Trig1Fired << std::endl; 
        }
        if (foundDoubleMu2018Trig1 != std::string::npos){
          doubleMu2018Trig1Fired = true;
       //    std::cout << "doubleMu2018Trig1Fired:  " << doubleMu2018Trig1Fired << std::endl; 
        }
        
        if (foundDoubleMu2018Trig2 != std::string::npos){
          doubleMu2018Trig2Fired = true;
        //   std::cout << "doubleMu2018Trig2Fired:  " << doubleMu2018Trig2Fired << std::endl;
        }
        
        if (foundTripleMu2018Trig1 != std::string::npos){
          tripleMu2018Trig1Fired = true;
          // std::cout << "tripleMu2018Trig1Fired:  " << tripleMu2018Trig1Fired << std::endl; 
        }
        
        if (foundTripleMu2018Trig2 != std::string::npos){
          tripleMu2018Trig2Fired = true;
         //  std::cout << "tripleMu2018Trig2Fired:  " << tripleMu2018Trig2Fired << std::endl; 
        }
        
        if (singleMu2018Trig1Fired || doubleMu2018Trig1Fired || doubleMu2018Trig2Fired || tripleMu2018Trig1Fired || tripleMu2018Trig2Fired){
          event_fails_trigger = false;
        //  std::cout << "Event passed 2018 triggers" << std::endl;
          break;
       
        }
        
      }
    }
    
    if (event_fails_trigger){
      continue;
    }
    
    for (int i=0; i<(int)TREE->lepton1_pt->size(); i++) {
      TLorentzVector lepton1, lepton2, lepton3, lepton4;
      lepton1.SetPtEtaPhiM ( TREE->lepton1_pt->at(i), TREE->lepton1_eta->at(i), TREE->lepton1_phi->at(i), muon_mass);
      lepton2.SetPtEtaPhiM ( TREE->lepton2_pt->at(i), TREE->lepton2_eta->at(i), TREE->lepton2_phi->at(i), muon_mass);
      lepton3.SetPtEtaPhiM ( TREE->lepton3_pt->at(i), TREE->lepton3_eta->at(i), TREE->lepton3_phi->at(i), muon_mass);
      lepton4.SetPtEtaPhiM ( TREE->lepton4_pt->at(i), TREE->lepton4_eta->at(i), TREE->lepton4_phi->at(i), muon_mass);
      
      //std::cout << TREE->flagZOnly->at(i) << std::endl;
      if (TREE->flagZOnly->at(i) == 0){
        continue;
      }
      
      if (TREE->flagZplusY->at(i) ==1){
        continue; 
      } //ensure we don't double count quads. Give priority to ZplusY candidates
      
      h_cutflow_allQuadCuts->Fill(1); //here are all the candidate flagZOnly and NOT flagZplusY quads
     
     //iso03 cuts here
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
    
       if (pfIso_lep1 > pfIso_Cut || pfIso_lep2 > pfIso_Cut || pfIso_lep3 > pfIso_Cut || pfIso_lep4 > pfIso_Cut){
         continue;
       }
      
      h_cutflow_allQuadCuts->Fill(2); //quads that pass the pfIso cut
      
       if (lepton1.DeltaR(lepton2) < lepSeparationCut || lepton1.DeltaR(lepton3) < lepSeparationCut || lepton1.DeltaR(lepton4) < lepSeparationCut || lepton2.DeltaR(lepton3) < lepSeparationCut || lepton2.DeltaR(lepton4) < lepSeparationCut || lepton3.DeltaR(lepton4) < lepSeparationCut){
         continue;
       } //separation cut for leptons of same flavor, as we have here since we have all muons
       
       h_cutflow_allQuadCuts->Fill(3); //quads that have their constituent muons appropriately separated and so pass the lepSeparationCut
       
       if (fabs(lepton1.Eta()) >etaCut || fabs(lepton2.Eta()) > etaCut || fabs(lepton3.Eta()) > etaCut || fabs(lepton4.Eta()) > etaCut){
         continue; 
       }
      
      h_cutflow_allQuadCuts->Fill(4); // quads that have all four muons within the bounds of  plus/minus 2.4 in eta  
     
      
      if (lepton1.Pt() < lepton1_pT_Cut || lepton2.Pt() < lepton2_pT_Cut || lepton3.Pt() < lepton3_pT_Cut || lepton4.Pt() < lepton4_pT_Cut){
        continue;
      }
      
      h_cutflow_allQuadCuts->Fill(5); //quads in which all the muons survive their appropriate pT cuts
      
      if (TREE->lepton1_isTightMuon->at(i) + TREE->lepton2_isTightMuon->at(i) + TREE->lepton3_isTightMuon->at(i) + TREE->lepton4_isTightMuon->at(i) != 4){
        continue;
      }
      
      h_cutflow_allQuadCuts->Fill(6); //quads in which all muons pass the Tight ID 
      
      if (fabs(TREE->lepton1_impactParameterSignificance->at(i)) > lep_3DIPSig_Cut || fabs(TREE->lepton2_impactParameterSignificance->at(i)) > lep_3DIPSig_Cut){
        continue;
      }
      
       h_cutflow_allQuadCuts->Fill(7); //quads in which the the lead two leptons pass the lep_3DIPSig_Cut
       
       if (fabs(TREE->lepton3_impactParameterSignificance->at(i)) > lep_3DIPSig_Cut || fabs(TREE->lepton4_impactParameterSignificance->at(i)) > lep_3DIPSig_Cut){
        continue;
      } 
      
      h_cutflow_allQuadCuts->Fill(8); //quads in which the two trailing leptons pass the lep_3DIPSig_Cut
      
      if (fabs(TREE->lepton1_dxy->at(i)) > lep_dxy_Cut || fabs(TREE->lepton2_dxy->at(i)) > lep_dxy_Cut ){
        continue;
      }
      h_cutflow_allQuadCuts->Fill(9); //quads in which the two lead leptons satisfy the lep_dxy_Cut requirements
      
      if (fabs(TREE->lepton3_dxy->at(i)) > lep_dxy_Cut || fabs(TREE->lepton4_dxy->at(i)) > lep_dxy_Cut ){
        continue;
      }
      
      h_cutflow_allQuadCuts->Fill(10); //quads in which the two trailing leptons satisfy the lep_dxy_Cut requirements
      
      if ( fabs(TREE->lepton1_dz->at(i)) > lep_dz_Cut || fabs(TREE->lepton2_dz->at(i)) > lep_dz_Cut ){
        continue;
      }
      
      h_cutflow_allQuadCuts->Fill(11); //quads in which the leading two leptons satisfy the lep_dz_Cut requirements 
      
      if ( fabs(TREE->lepton3_dz->at(i)) > lep_dz_Cut || fabs(TREE->lepton4_dz->at(i)) > lep_dz_Cut ){
        continue;
      }
      
       h_cutflow_allQuadCuts->Fill(12); //quads in which the trailing two leptons satisfy the lep_dz_Cut requirements 
      
      if (doRecoToTrigMuMatching){ 
        if (TREE->quadHasHowManyTrigMatches->at(i) < 2) {
          continue;
        }
      }
      
      h_cutflow_allQuadCuts->Fill(13); //quads that satisfy the requirement of having at least 2 reco to trigger muon matches
      
      int theSum;
      theSum = TREE->pair_12_34_ZOnly->at(i) + TREE->pair_13_24_ZOnly->at(i) + TREE->pair_14_23_ZOnly->at(i); 
      h_pair_sum->Fill(theSum);
      
      bool is_pair_12_34_ZOnly = false; bool is_pair_13_24_ZOnly = false; bool is_pair_14_23_ZOnly = false;
      
      if (theSum == 1){
      //  std::cout << "ELEPHANT" << std::endl; 
        if (TREE->pair_12_34_ZOnly->at(i) ==1){
          is_pair_12_34_ZOnly = true;
        }
        
        if (TREE->pair_13_24_ZOnly->at(i) == 1){
          is_pair_13_24_ZOnly = true;
        }
        
        if (TREE->pair_14_23_ZOnly->at(i) == 1){
          is_pair_14_23_ZOnly = true; 
        }
      
      }
      
      //theSum is usually 2
      if (theSum == 2){
        
        if(TREE->pair_12_34_ZOnly->at(i) + TREE->pair_13_24_ZOnly->at(i) == 2){
        
          double diff_12_34_ZOnly, diff_13_24_ZOnly;
        
          diff_12_34_ZOnly = fabs((lepton1 + lepton2).M() - (lepton3 + lepton4).M());
          //std::cout << "diff_12_34_ZOnly:  " << diff_12_34_ZOnly << std::endl; 
        
          diff_13_24_ZOnly = fabs((lepton1 + lepton3).M() -(lepton2 + lepton4).M());
         // std::cout << "diff_13_24_ZOnly:  " << diff_13_24_ZOnly << std::endl;
        
          if (diff_12_34_ZOnly > diff_13_24_ZOnly){
            is_pair_12_34_ZOnly = true;
           // std::cout << "is_pair_12_34_ZOnly = true" << std::endl; 
          }
          else{
            is_pair_13_24_ZOnly = true;
          //  std::cout << "is_pair_13_24_ZOnly = true" << std::endl; 
          }
        
        }
        
        if(TREE->pair_12_34_ZOnly->at(i) + TREE->pair_14_23_ZOnly->at(i) == 2){
          
          double diff_12_34_ZOnly, diff_14_23_ZOnly;
          
          diff_12_34_ZOnly = fabs((lepton1 + lepton2).M() - (lepton3 + lepton4).M());
         // std::cout << "diff_12_34_ZOnly:  " << diff_12_34_ZOnly << std::endl;
          
          diff_14_23_ZOnly = fabs( (lepton1 + lepton4).M() - (lepton2 + lepton3).M() );
          // std::cout << "diff_14_23_ZOnly:  " << diff_14_23_ZOnly << std::endl; 
          
          if (diff_12_34_ZOnly > diff_14_23_ZOnly){
            is_pair_12_34_ZOnly = true;
            //std::cout << "is_pair_12_34_ZOnly = true" << std::endl; 
          } 
          else{
            is_pair_14_23_ZOnly = true;
           // std::cout << "is_pair_14_23_ZOnly = true" << std::endl; 
          }        
         
        }
        
        if(TREE->pair_13_24_ZOnly->at(i) + TREE->pair_14_23_ZOnly->at(i) == 2) {
          
          double diff_13_24_ZOnly, diff_14_23_ZOnly;
          
          diff_13_24_ZOnly = fabs((lepton1 + lepton3).M() -(lepton2 + lepton4).M());
          //std::cout << "diff_13_24_ZOnly:  " << diff_13_24_ZOnly << std::endl; 
          
          diff_14_23_ZOnly = fabs((lepton1+lepton4).M() - (lepton2+lepton3).M());
          //std::cout << "diff_14_23_ZOnly:  " << diff_14_23_ZOnly << std::endl; 
          
          if (diff_13_24_ZOnly > diff_14_23_ZOnly){
            is_pair_13_24_ZOnly = true;
           // std::cout << "pair_13_24_ZOnly = true" << std::endl; 
          }
          else{
            is_pair_14_23_ZOnly = true;
           // std::cout << "pair_14_23_ZOnly = true" << std::endl; 
          }
        }
      
      }
      
      //We now have three mutually exclusive categories: is_pair_12_34_ZOnly, is_pair_13_24_ZOnly, and is_pair_14_23_ZOnly 
      if (is_pair_12_34_ZOnly){
        is_pair_12_34_ZOnly_Count++;
        
        //use the slot for the fourth match for the inner vector in the vector of vectors, aka the slot indexed by 3 (since indexing starts at 0)
        //i references what quad we are on (aka what we row in our vector of vectors we are at), j here is 3 because we want to look in the fourth column
        //because that's where info related to pair_12_34_ZOnly was stored in the phase1 code
        TVector3 dimuon1vtx_vec, dimuon2vtx_vec;
        dimuon1vtx_vec.SetXYZ(TREE->dimuon1vtx_xpos->at(i).at(3), TREE->dimuon1vtx_ypos->at(i).at(3), TREE->dimuon1vtx_zpos->at(i).at(3));
        //std::cout << TREE->dimuon1vtx_xpos->at(i).at(3) << std::endl; 
        
        dimuon2vtx_vec.SetXYZ(TREE->dimuon2vtx_xpos->at(i).at(3), TREE->dimuon2vtx_ypos->at(i).at(3), TREE->dimuon2vtx_zpos->at(i).at(3));
        
        //Protections against vectors whose vertex coordinates got filled with the value (-1000) that indicates that the vertex was found to be not valid in the phase1 code
        if (dimuon1vtx_vec.X() == -1000 || dimuon1vtx_vec.Y() == -1000 || dimuon1vtx_vec.Z() == -1000){
          continue; 
        }
        
        if (dimuon2vtx_vec.X() == -1000 || dimuon2vtx_vec.Y() == -1000 || dimuon2vtx_vec.Z() == -1000){
          continue; 
        }
     //    if (dimuon1vtx_vec.DeltaR(dimuon2vtx_vec) > deltaR_dimuon1vtx_dimuon2vtx_Cut){
//           //std::cout << "WATER BUFFALO" << std::endl;
//           continue; 
//         }
        
        double dR = dimuon1vtx_vec.DeltaR(dimuon2vtx_vec);
        double dZ = fabs(dimuon1vtx_vec.Z()-dimuon2vtx_vec.Z());
        
    //     if (dZ > dR + offset){
//           continue; 
//         }
//         
        //std::cout << TREE->dimuon1vtx->at(i).at(3) << std::endl; 
        
        if (TREE->dimuon1vtx->at(i).at(3) < dimuonvtx_Prob_Cut){
          continue;
        }
        
        if (TREE->dimuon2vtx->at(i).at(3) < dimuonvtx_Prob_Cut){
          continue;
        }
        //cartesian DR significance calculations
        
        double DX = fabs(dimuon1vtx_vec.X()-dimuon2vtx_vec.X());
        double DX2 = DX * DX;
   //     std::cout << "DX2:  " << DX2 << std::endl; 
        
        double error_dimu1_X2 = TREE->dimuon1vtx_xposError->at(i).at(3);
   //     std::cout << "error_dimu1_X:  " << error_dimu1_X << std::endl;
        
      //  double error_dimu1_X2 = error_dimu1_X * error_dimu1_X;
  //      std::cout << "error_dimu1_X2:  " << error_dimu1_X2 << std::endl; 
         
        double error_dimu2_X2 = TREE->dimuon2vtx_xposError->at(i).at(3);
       // double error_dimu2_X2 = error_dimu2_X * error_dimu2_X;
  //      std::cout << "error_dimu2_X2:  " << error_dimu2_X2 << std::endl;
         
        double X_err_sum_in_quad = error_dimu1_X2 + error_dimu2_X2;
 //       std::cout << "X_err_sum_in_quad:  " << X_err_sum_in_quad << std::endl; 
         
        double first_term = DX2/X_err_sum_in_quad;
 //       std::cout << "first_term:  " << first_term << std::endl; 
         
        double DY = fabs(dimuon1vtx_vec.Y() - dimuon2vtx_vec.Y());
        double DY2 = DY * DY; 
//        std::cout << "DY2:  " << DY2 << std::endl; 
         
        double error_dimu1_Y2 = TREE->dimuon1vtx_yposError->at(i).at(3);
       // double error_dimu1_Y2 = error_dimu1_Y * error_dimu1_Y; 
//        std::cout << "error_dimu1_Y2:  " << error_dimu1_Y2 << std::endl;
         
        double error_dimu2_Y2 = TREE->dimuon2vtx_yposError->at(i).at(3);
      //  double error_dimu2_Y2 = error_dimu2_Y * error_dimu2_Y;
//        std::cout << "error_dimu2_Y2:  " << error_dimu2_Y2 << std::endl; 
         
        double Y_err_sum_in_quad = error_dimu1_Y2 + error_dimu2_Y2;
//        std::cout << "Y_err_sum_in_quad:  " << Y_err_sum_in_quad << std::endl; 
         
        double second_term = DY2/Y_err_sum_in_quad;
 //       std::cout << "second_term:  " << second_term << std::endl; 
        
        double cart_DR_Sig2 = first_term + second_term;
//        std::cout << "cart_DR_Sig2:  " << cart_DR_Sig2 << std::endl; 
        
        double cart_DR_Sig = TMath::Sqrt(cart_DR_Sig2);
        std::cout << "cart_DR_Sig:  " << cart_DR_Sig << std::endl; 
        
        h_cart_DR_Sig->Fill(cart_DR_Sig);
        
        double DZ2 = dZ * dZ;
        double error_dimu1_Z2 = TREE->dimuon1vtx_zposError->at(i).at(3);
        //double error_dimu1_Z2 = error_dimu1_Z * error_dimu1_Z;
        double error_dimu2_Z2 = TREE->dimuon2vtx_zposError->at(i).at(3);
       // double error_dimu2_Z2 = error_dimu2_Z * error_dimu2_Z;
        double Z_err_sum_in_quad = error_dimu1_Z2 + error_dimu2_Z2;
        double DZ_Sig2 = DZ2/Z_err_sum_in_quad;
        double DZ_Sig = TMath::Sqrt(DZ_Sig2);
        std::cout << "DZ_Sig:  " << DZ_Sig << std::endl;
        
        h_DZ_Sig->Fill(DZ_Sig);
        
                
        //Sort out Z1, the heavier of the two pairs, vs. Z2, the lighter of the two pairs
        
        if ( (lepton1 + lepton2).M()  >  (lepton3 + lepton4).M() ){
          temp_Z1_mass.push_back( (lepton1 + lepton2).M() );
          temp_Z2_mass.push_back( (lepton3 + lepton4).M() );
          
          temp_Z1_pT.push_back( (lepton1 + lepton2).Pt() );
          temp_Z2_pT.push_back( (lepton3 + lepton4).Pt() );
          
          temp_Z1_eta.push_back( (lepton1 + lepton2).Eta() );
          temp_Z2_eta.push_back( (lepton3 + lepton4).Eta() );
          
          temp_Z1_phi.push_back( (lepton1 + lepton2).Phi() );
          temp_Z2_phi.push_back( (lepton3 + lepton4).Phi() ); 
        }
        else{
          temp_Z1_mass.push_back( (lepton3 + lepton4).M() );
          temp_Z2_mass.push_back ((lepton1 + lepton2).M() );
          
          temp_Z1_pT.push_back( (lepton3 + lepton4).Pt() );
          temp_Z2_pT.push_back( (lepton1 + lepton2).Pt() );
          
          temp_Z1_eta.push_back( (lepton3 + lepton4).Eta() );
          temp_Z2_eta.push_back( (lepton1 + lepton2).Eta() );
          
          temp_Z1_phi.push_back( (lepton3 + lepton4).Phi() );
          temp_Z2_phi.push_back( (lepton1 + lepton2).Phi() );
          
        }
      }
      
      if (is_pair_13_24_ZOnly){
        is_pair_13_24_ZOnly_Count++;
        
       
        //use the slot for the fifth match for in the inner vector in the vector of vectors, aka the slot indexed by 4 (since indexing starts at 0)
        // i references what quad we are on (aka what row in our vector of vectors we are at), j here is 4 because we want to look at the fifth column 
        //because that's where info related to the pair_13_24_ZOnly was stored in the phase1 code
        TVector3 dimuon1vtx_vec, dimuon2vtx_vec;
        dimuon1vtx_vec.SetXYZ(TREE->dimuon1vtx_xpos->at(i).at(4), TREE->dimuon1vtx_ypos->at(i).at(4), TREE->dimuon1vtx_zpos->at(i).at(4));
        dimuon2vtx_vec.SetXYZ(TREE->dimuon2vtx_xpos->at(i).at(4), TREE->dimuon2vtx_ypos->at(i).at(4), TREE->dimuon2vtx_zpos->at(i).at(4));
        
         //Protections against vectors whose vertex coordinates got filled with the value (-1000) that indicates that the vertex was found to be not valid in the phase1 code
        if (dimuon1vtx_vec.X() == -1000 || dimuon1vtx_vec.Y() == -1000 || dimuon1vtx_vec.Z() == -1000){
          continue; 
        }
        
        if (dimuon2vtx_vec.X() == -1000 || dimuon2vtx_vec.Y() == -1000 || dimuon2vtx_vec.Z() == -1000){
          continue; 
        }
        
     //    if (dimuon1vtx_vec.DeltaR(dimuon2vtx_vec) > deltaR_dimuon1vtx_dimuon2vtx_Cut){
//           //std::cout << "WATER BUFFALO 2" << std::endl;
//           continue; 
//         }
        
        double dR = dimuon1vtx_vec.DeltaR(dimuon2vtx_vec);
        double dZ = fabs(dimuon1vtx_vec.Z()-dimuon2vtx_vec.Z());
        
    //     if (dZ > dR + offset){
//           continue; 
//         }
        
         if (TREE->dimuon1vtx->at(i).at(4) < dimuonvtx_Prob_Cut){
          continue;
        }
        
        if (TREE->dimuon2vtx->at(i).at(4) < dimuonvtx_Prob_Cut){
          continue;
        }
        
                //cartesian DR significance calculations
        
        double DX = fabs(dimuon1vtx_vec.X()-dimuon2vtx_vec.X());
        double DX2 = DX * DX;
   //     std::cout << "DX2:  " << DX2 << std::endl; 
        
        double error_dimu1_X2 = TREE->dimuon1vtx_xposError->at(i).at(4);
   //     std::cout << "error_dimu1_X:  " << error_dimu1_X << std::endl;
        
      //  double error_dimu1_X2 = error_dimu1_X * error_dimu1_X;
  //      std::cout << "error_dimu1_X2:  " << error_dimu1_X2 << std::endl; 
         
        double error_dimu2_X2 = TREE->dimuon2vtx_xposError->at(i).at(4);
       // double error_dimu2_X2 = error_dimu2_X * error_dimu2_X;
  //      std::cout << "error_dimu2_X2:  " << error_dimu2_X2 << std::endl;
         
        double X_err_sum_in_quad = error_dimu1_X2 + error_dimu2_X2;
 //       std::cout << "X_err_sum_in_quad:  " << X_err_sum_in_quad << std::endl; 
         
        double first_term = DX2/X_err_sum_in_quad;
 //       std::cout << "first_term:  " << first_term << std::endl; 
         
        double DY = fabs(dimuon1vtx_vec.Y() - dimuon2vtx_vec.Y());
        double DY2 = DY * DY; 
//        std::cout << "DY2:  " << DY2 << std::endl; 
         
        double error_dimu1_Y2 = TREE->dimuon1vtx_yposError->at(i).at(4);
       // double error_dimu1_Y2 = error_dimu1_Y * error_dimu1_Y; 
//        std::cout << "error_dimu1_Y2:  " << error_dimu1_Y2 << std::endl;
         
        double error_dimu2_Y2 = TREE->dimuon2vtx_yposError->at(i).at(4);
      //  double error_dimu2_Y2 = error_dimu2_Y * error_dimu2_Y;
//        std::cout << "error_dimu2_Y2:  " << error_dimu2_Y2 << std::endl; 
         
        double Y_err_sum_in_quad = error_dimu1_Y2 + error_dimu2_Y2;
//        std::cout << "Y_err_sum_in_quad:  " << Y_err_sum_in_quad << std::endl; 
         
        double second_term = DY2/Y_err_sum_in_quad;
 //       std::cout << "second_term:  " << second_term << std::endl; 
        
        double cart_DR_Sig2 = first_term + second_term;
//        std::cout << "cart_DR_Sig2:  " << cart_DR_Sig2 << std::endl; 
        
        double cart_DR_Sig = TMath::Sqrt(cart_DR_Sig2);
        std::cout << "cart_DR_Sig Elephant:  " << cart_DR_Sig << std::endl; 
        
        h_cart_DR_Sig->Fill(cart_DR_Sig);
        
        double DZ2 = dZ * dZ;
        double error_dimu1_Z2 = TREE->dimuon1vtx_zposError->at(i).at(4);
      //  double error_dimu1_Z2 = error_dimu1_Z * error_dimu1_Z;
        double error_dimu2_Z2 = TREE->dimuon2vtx_zposError->at(i).at(4);
       // double error_dimu2_Z2 = error_dimu2_Z * error_dimu2_Z;
        double Z_err_sum_in_quad = error_dimu1_Z2 + error_dimu2_Z2;
        double DZ_Sig2 = DZ2/Z_err_sum_in_quad;
        double DZ_Sig = TMath::Sqrt(DZ_Sig2);
        std::cout << "DZ_Sig:  " << DZ_Sig << std::endl;
        
        h_DZ_Sig->Fill(DZ_Sig);
        
        //sort out Z1, the heavier of the two pairs, vs. Z2, the lighter of the two pairs
        
        if ( (lepton1 + lepton3).M() > (lepton2 + lepton4).M() ){
          temp_Z1_mass.push_back( (lepton1 + lepton3).M() );
          temp_Z2_mass.push_back( (lepton2 + lepton4).M() );
          
          temp_Z1_pT.push_back( (lepton1 + lepton3).Pt() );
          temp_Z2_pT.push_back( (lepton2 + lepton4).Pt() );
          
          temp_Z1_eta.push_back( (lepton1 + lepton3).Eta() );
          temp_Z2_eta.push_back( (lepton2 + lepton4).Eta() );
          
          temp_Z1_phi.push_back( (lepton1 + lepton3).Phi() );
          temp_Z2_phi.push_back( (lepton2 + lepton4).Phi() );
        }
        else{
          temp_Z1_mass.push_back( (lepton2 + lepton4).M() );
          temp_Z2_mass.push_back( (lepton1 + lepton3).M() );
          
          temp_Z1_pT.push_back( (lepton2 + lepton4).Pt() );
          temp_Z2_pT.push_back( (lepton1 + lepton3).Pt() );
          
          temp_Z1_eta.push_back( (lepton2 + lepton4).Eta() );
          temp_Z2_eta.push_back( (lepton1 + lepton3).Eta() );
          
          temp_Z1_phi.push_back( (lepton2 + lepton4).Phi() );
          temp_Z2_phi.push_back( (lepton1 + lepton3).Phi() );
          
          
        }
      }
      
      if (is_pair_14_23_ZOnly){
        is_pair_14_23_ZOnly_Count++;
        //use the slot for the sixth match for the inner vector in the vector of vectors, aka the slot indexed by 5 (since indexing starts at 0)
        //i references what quad we are on (aka what row in our vector of vectors we are at at), j here is 5 because we want to look at the sixth column
        //because that's where info related to the match_14_23_ZOnly was stored in the phase1 code
        TVector3 dimuon1vtx_vec, dimuon2vtx_vec;
        dimuon1vtx_vec.SetXYZ(TREE->dimuon1vtx_xpos->at(i).at(5), TREE->dimuon1vtx_ypos->at(i).at(5), TREE->dimuon1vtx_zpos->at(i).at(5));
        dimuon2vtx_vec.SetXYZ(TREE->dimuon2vtx_xpos->at(i).at(5), TREE->dimuon2vtx_ypos->at(i).at(5), TREE->dimuon2vtx_zpos->at(i).at(5));
        
         //Protections against vectors whose vertex coordinates got filled with the value (-1000) that indicates that the vertex was found to be not valid in the phase1 code
        if (dimuon1vtx_vec.X() == -1000 || dimuon1vtx_vec.Y() == -1000 || dimuon1vtx_vec.Z() == -1000){
          continue; 
        }
        
        if (dimuon2vtx_vec.X() == -1000 || dimuon2vtx_vec.Y() == -1000 || dimuon2vtx_vec.Z() == -1000){
          continue; 
        }
        
      //   if (dimuon1vtx_vec.DeltaR(dimuon2vtx_vec) > deltaR_dimuon1vtx_dimuon2vtx_Cut){
//           //std::cout << "WATER BUFFALO 3" << std::endl;
//           continue; 
//         }
        
        double dR = dimuon1vtx_vec.DeltaR(dimuon2vtx_vec);
        double dZ = fabs(dimuon1vtx_vec.Z()-dimuon2vtx_vec.Z());
        
  //       if (dZ > dR + offset){
//           continue; 
//         }
        
         if (TREE->dimuon1vtx->at(i).at(5) < dimuonvtx_Prob_Cut){
          continue;
        }
        
        if (TREE->dimuon2vtx->at(i).at(5) < dimuonvtx_Prob_Cut){
          continue;
        }
        
                //cartesian DR significance calculations
        
        double DX = fabs(dimuon1vtx_vec.X()-dimuon2vtx_vec.X());
        double DX2 = DX * DX;
   //     std::cout << "DX2:  " << DX2 << std::endl; 
        
        double error_dimu1_X2 = TREE->dimuon1vtx_xposError->at(i).at(5);
   //     std::cout << "error_dimu1_X:  " << error_dimu1_X << std::endl;
        
      //  double error_dimu1_X2 = error_dimu1_X * error_dimu1_X;
  //      std::cout << "error_dimu1_X2:  " << error_dimu1_X2 << std::endl; 
         
        double error_dimu2_X2 = TREE->dimuon2vtx_xposError->at(i).at(5);
     //   double error_dimu2_X2 = error_dimu2_X * error_dimu2_X;
  //      std::cout << "error_dimu2_X2:  " << error_dimu2_X2 << std::endl;
         
        double X_err_sum_in_quad = error_dimu1_X2 + error_dimu2_X2;
 //       std::cout << "X_err_sum_in_quad:  " << X_err_sum_in_quad << std::endl; 
         
        double first_term = DX2/X_err_sum_in_quad;
 //       std::cout << "first_term:  " << first_term << std::endl; 
         
        double DY = fabs(dimuon1vtx_vec.Y() - dimuon2vtx_vec.Y());
        double DY2 = DY * DY; 
//        std::cout << "DY2:  " << DY2 << std::endl; 
         
        double error_dimu1_Y2 = TREE->dimuon1vtx_yposError->at(i).at(5);
       // double error_dimu1_Y2 = error_dimu1_Y * error_dimu1_Y; 
//        std::cout << "error_dimu1_Y2:  " << error_dimu1_Y2 << std::endl;
         
        double error_dimu2_Y2 = TREE->dimuon2vtx_yposError->at(i).at(5);
     //   double error_dimu2_Y2 = error_dimu2_Y * error_dimu2_Y;
//        std::cout << "error_dimu2_Y2:  " << error_dimu2_Y2 << std::endl; 
         
        double Y_err_sum_in_quad = error_dimu1_Y2 + error_dimu2_Y2;
//        std::cout << "Y_err_sum_in_quad:  " << Y_err_sum_in_quad << std::endl; 
         
        double second_term = DY2/Y_err_sum_in_quad;
 //       std::cout << "second_term:  " << second_term << std::endl; 
        
        double cart_DR_Sig2 = first_term + second_term;
//        std::cout << "cart_DR_Sig2:  " << cart_DR_Sig2 << std::endl; 
        
        double cart_DR_Sig = TMath::Sqrt(cart_DR_Sig2);
        std::cout << "cart_DR_Sig Wombat:  " << cart_DR_Sig << std::endl; 
        
        h_cart_DR_Sig->Fill(cart_DR_Sig);
        
        double DZ2 = dZ * dZ;
        double error_dimu1_Z2 = TREE->dimuon1vtx_zposError->at(i).at(5);
   //     double error_dimu1_Z2 = error_dimu1_Z * error_dimu1_Z;
        double error_dimu2_Z2 = TREE->dimuon2vtx_zposError->at(i).at(5);
      //  double error_dimu2_Z2 = error_dimu2_Z * error_dimu2_Z;
        double Z_err_sum_in_quad = error_dimu1_Z2 + error_dimu2_Z2;
        double DZ_Sig2 = DZ2/Z_err_sum_in_quad;
        double DZ_Sig = TMath::Sqrt(DZ_Sig2);
        std::cout << "DZ_Sig:  " << DZ_Sig << std::endl;
        
        h_DZ_Sig->Fill(DZ_Sig);
        
        
        //sort out Z1, the heavier of the two pairs, vs. Z2, the lighter of the two pairs
        
        if ( (lepton1 + lepton4).M() > (lepton2 + lepton3).M() ){
          temp_Z1_mass.push_back( (lepton1 + lepton4).M() );
          temp_Z2_mass.push_back( (lepton2 + lepton3).M() );
          
          temp_Z1_pT.push_back( (lepton1 + lepton4).Pt() );
          temp_Z2_pT.push_back( (lepton2 + lepton3).Pt() );
          
          temp_Z1_eta.push_back( (lepton1 + lepton4).Eta() );
          temp_Z2_eta.push_back( (lepton2 + lepton3).Eta() );
          
          temp_Z1_phi.push_back( (lepton1 + lepton4).Phi() );
          temp_Z2_phi.push_back( (lepton2 + lepton3).Phi() );
        }
        else{
          temp_Z1_mass.push_back( (lepton2 + lepton3).M() );
          temp_Z2_mass.push_back( (lepton1 + lepton4).M() );
          
          temp_Z1_pT.push_back( (lepton2 + lepton3).Pt() );
          temp_Z2_pT.push_back( (lepton1 + lepton4).Pt() );
          
          temp_Z1_eta.push_back( (lepton2 + lepton3).Eta() );
          temp_Z2_eta.push_back( (lepton1 + lepton4).Eta() );
          
          temp_Z1_phi.push_back( (lepton2 + lepton3).Phi() );
          temp_Z2_phi.push_back( (lepton1 + lepton4).Phi() );
        }
      }
      //If we get here, we have a survivor
      temp_Z_mass.push_back( (lepton1 + lepton2 + lepton3 + lepton4).M() );
      temp_Z_eta.push_back( (lepton1 + lepton2 + lepton3 + lepton4).Eta() );
      temp_Z_phi.push_back ( (lepton1 + lepton2 + lepton3 + lepton4).Phi() );
      temp_Z_pT.push_back( (lepton1 + lepton2 + lepton3 + lepton4).Pt() );
      
      temp_lep1_pT.push_back(lepton1.Pt());
      temp_lep1_eta.push_back(lepton1.Eta());
      temp_lep1_phi.push_back(lepton1.Phi());
      
      temp_lep2_pT.push_back(lepton2.Pt());
      temp_lep2_eta.push_back(lepton2.Eta());
      temp_lep2_phi.push_back(lepton2.Phi());
      
      temp_lep3_pT.push_back(lepton3.Pt());
      temp_lep3_eta.push_back(lepton3.Eta());
      temp_lep3_phi.push_back(lepton3.Phi());
      
      temp_lep4_pT.push_back(lepton4.Pt());
      temp_lep4_eta.push_back(lepton4.Eta());
      temp_lep4_phi.push_back(lepton4.Phi());
      
      

    } //close loop over leptons
    if (temp_Z_mass.size() == 1){
      fillCount++;
      
      Z_mass = temp_Z_mass.at(0);
      Z_eta  = temp_Z_eta.at(0);
      Z_phi  = temp_Z_phi.at(0);
      Z_pT   = temp_Z_pT.at(0);
      
      
      
      //temp_Z_mass.size() == 1 implies that the temp_Z1, temp_Z2 quantities also are of size 1. This is true because 
      //there are no continue statements between where the temp_Z1, temp_Z2 quantities are filled and the place
      //where the temp_Z_mass quantity is filled, so whatever is filled at the temp_Z1, temp_Z2 level will
      //fall through to be filled at the temp_Z_mass level.
      //If we had an event that contained an is_pair_XY_PQ and an is_pair_XprimeYprime_PprimeQprime, the temp_Z_mass
      //size could be greater than the individual temp_Z1, temp_Z2 sizes,
      //in this case the temp_Z_mass size would be the sum of the temp_Z1 (temp_Z2) size from is_pair_XY_PQ and
      //from is_pair_XprimeYprime_PprimeQprime
        
      Z1_mass = temp_Z1_mass.at(0);
      Z2_mass = temp_Z2_mass.at(0);
      
      Z1_pT   = temp_Z1_pT.at(0);
      Z2_pT   = temp_Z2_pT.at(0);
      
      Z1_eta  = temp_Z1_eta.at(0);
      Z2_eta  = temp_Z2_eta.at(0);
      
      Z1_phi  = temp_Z1_phi.at(0);
      Z2_phi  = temp_Z2_phi.at(0);
      
      lepton1_pT = temp_lep1_pT.at(0);
      lepton1_eta = temp_lep1_eta.at(0);
      lepton1_phi = temp_lep1_phi.at(0);
      
      lepton2_pT = temp_lep2_pT.at(0);
      lepton2_eta = temp_lep2_eta.at(0);
      lepton2_phi = temp_lep2_phi.at(0);
      
      lepton3_pT = temp_lep3_pT.at(0);
      lepton3_eta = temp_lep3_eta.at(0);
      lepton3_phi = temp_lep3_phi.at(0);
      
      lepton4_pT = temp_lep4_pT.at(0);
      lepton4_eta = temp_lep4_eta.at(0);
      lepton4_phi = temp_lep4_phi.at(0);
      
      
      
      
      
      aux->Fill();
    }
  }//close loop over entries
  
  std::cout << "is_pair_12_34_ZOnly_Count:  " << is_pair_12_34_ZOnly_Count << std::endl;
  std::cout << "is_pair_13_24_ZOnly_Count:  " << is_pair_13_24_ZOnly_Count << std::endl;
  std::cout << "is_pair_14_23_ZOnly_Count:  " << is_pair_14_23_ZOnly_Count << std::endl; 
  std::cout << "fillCount: " << fillCount << std::endl; 
  /////////////////////////////////////////////////////////
////////////////     P L O T T I N G     ////////////////
/////////////////////////////////////////////////////////
  
 TCanvas *c_cutflow_allQuadCuts = new TCanvas("c_cutflow_allQuadCuts", "c_cutflow_allQuadCuts");
 c_cutflow_allQuadCuts->cd();
 h_cutflow_allQuadCuts->Draw();
 h_cutflow_allQuadCuts->Write();
 c_cutflow_allQuadCuts->SaveAs("c_cutflow_allQuadCuts_ZOnly.pdf"); 
  
 TCanvas *c_pfIso_lepN = new TCanvas("c_pfIso_lepN", "c_pfIso_lepN"); c_pfIso_lepN->Divide(2,2);
 c_pfIso_lepN->cd(1); h_pfIso_lep1->Draw();
 c_pfIso_lepN->cd(2); h_pfIso_lep2->Draw();
 c_pfIso_lepN->cd(3); h_pfIso_lep3->Draw();
 c_pfIso_lepN->cd(4); h_pfIso_lep4->Draw();
 h_pfIso_lep1->Write(); h_pfIso_lep2->Write(); h_pfIso_lep3->Write(); h_pfIso_lep4->Write();
 c_pfIso_lepN->SaveAs("c_pfIso_lepN_ZOnly.pdf");  

TCanvas *c_pair_sum = new TCanvas("c_pair_sum", "c_pair_sum");
c_pair_sum->cd();
h_pair_sum->Draw();
h_pair_sum->Write();
c_pair_sum->SaveAs("c_pair_sum_ZOnly.pdf");

TCanvas *c_cart_DR_Sig = new TCanvas("c_cart_DR_Sig", "c_cart_DR_Sig");
c_cart_DR_Sig->cd();
h_cart_DR_Sig->Draw();
h_cart_DR_Sig->Write();
c_cart_DR_Sig->SaveAs("c_cart_DR_Sig_ZOnly.pdf");

TCanvas *c_DZ_Sig = new TCanvas("c_DZ_Sig", "c_DZ_Sig");
c_DZ_Sig->cd();
h_DZ_Sig->Draw();
h_DZ_Sig->Write();
c_DZ_Sig->SaveAs("c_DZ_Sig_ZOnly.pdf"); 
  
  
  
  
  
  ntuple->Write();
  ntuple->Close();


}//close run function