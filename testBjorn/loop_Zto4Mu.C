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
  
  // v a r i a b l e s
  double muon_mass = 105.6583 / 1000.; //get mass in GeV
  
  //Counters
  
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
  
  //New skimmed root file
  double Z_mass = -99;
  double Z_pT = -99;
  double Z_eta = -99;
  double Z_phi = -99;
  
  TFile *ntuple = new TFile("18May2022_testOut_loop_Zto4Mu.root", "RECREATE");
  TTree *aux;
  aux = new TTree("tree", "tree");
  
  aux->Branch("Z_mass", &Z_mass);
  aux->Branch("Z_pT", &Z_pT);
  aux->Branch("Z_eta", &Z_eta);
  aux->Branch("Z_phi", &Z_phi);
  
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
    
    temp_Z_mass.clear();
    temp_Z_pT.clear();
    temp_Z_eta.clear();
    temp_Z_phi.clear();
    
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
      
      int theSum;
      theSum = TREE->pair_12_34_ZOnly->at(i) + TREE->pair_13_24_ZOnly->at(i) + TREE->pair_14_23_ZOnly->at(i); 
      h_pair_sum->Fill(theSum);
      
      bool is_pair_12_34_ZOnly = false; bool is_pair_13_24_ZOnly = false; bool is_pair_14_23_ZOnly = false;
      
      if (theSum == 1){
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
      
      }

    } //close loop over leptons
    aux->Fill();
  
  }//close loop over entries
  
  
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
c_pair_sum->SaveAs("c_pair_sum.pdf");
  
  
  
  
  
  
  ntuple->Write();
  ntuple->Close();


}//close run function