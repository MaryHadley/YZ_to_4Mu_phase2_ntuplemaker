//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Mar 21 08:21:13 2022 by ROOT version 6.14/09
// from TTree treemc/treemc
// found on file: ZYto4Mu_Zto4Mu_pTCut3_MC_Bjorn_17March2022_inputFileIs_MC_DPS_2016_YZ_00623223-2B20-AB42-A456-670F9B3875D5_OR_DoubleMuTrigs.root
//////////////////////////////////////////////////////////

#ifndef treeMC_h
#define treeMC_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "vector"

class treeMC {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   vector<double>  *truth_Zmuon_pt;
   vector<double>  *truth_Zmuon_eta;
   vector<double>  *truth_Zmuon_phi;
   vector<double>  *truth_Z_pt;
   vector<double>  *truth_Z_eta;
   vector<double>  *truth_Z_phi;
   vector<double>  *truth_Z_mass;
   vector<double>  *truth_Z_pdgid;
   vector<double>  *truth_Upsimuon_pt;
   vector<double>  *truth_Upsimuon_eta;
   vector<double>  *truth_Upsimuon_phi;
   vector<double>  *truth_Upsi_pt;
   vector<double>  *truth_Upsi_eta;
   vector<double>  *truth_Upsi_phi;
   vector<double>  *truth_Upsi_mass;
   vector<double>  *truth_Upsi_pdgid;
   vector<double>  *truth_Upsi2muon_pt;
   vector<double>  *truth_Upsi2muon_eta;
   vector<double>  *truth_Upsi2muon_phi;
   vector<double>  *truth_Upsi2_pt;
   vector<double>  *truth_Upsi2_eta;
   vector<double>  *truth_Upsi2_phi;
   vector<double>  *truth_Upsi2_mass;
   vector<double>  *truth_Upsi2_pdgid;
   vector<double>  *truth_Upsi3muon_pt;
   vector<double>  *truth_Upsi3muon_eta;
   vector<double>  *truth_Upsi3muon_phi;
   vector<double>  *truth_Upsi3_pt;
   vector<double>  *truth_Upsi3_eta;
   vector<double>  *truth_Upsi3_phi;
   vector<double>  *truth_Upsi3_mass;
   vector<double>  *truth_Upsi3_pdgid;
   vector<double>  *truth_Chib0_1P_UPSI_muon_pt;
   vector<double>  *truth_Chib0_1P_UPSI_muon_eta;
   vector<double>  *truth_Chib0_1P_UPSI_muon_phi;
   vector<double>  *truth_Chib0_1P_pt;
   vector<double>  *truth_Chib0_1P_eta;
   vector<double>  *truth_Chib0_1P_phi;
   vector<double>  *truth_Chib0_1P_pdgid;
   vector<double>  *truth_Chib0_1P_mass;
   vector<double>  *truth_Chib0_1P_UPSI_pt;
   vector<double>  *truth_Chib0_1P_UPSI_phi;
   vector<double>  *truth_Chib0_1P_UPSI_eta;
   vector<double>  *truth_Chib0_1P_UPSI_pdgid;
   vector<double>  *truth_Chib0_1P_UPSI_mass;
   vector<double>  *truth_Chib0_1P_photon_pt;
   vector<double>  *truth_Chib0_1P_photon_phi;
   vector<double>  *truth_Chib0_1P_photon_eta;
   vector<double>  *truth_Chib0_1P_photon_pdgid;
   vector<double>  *truth_Chib1_1P_UPSI_muon_pt;
   vector<double>  *truth_Chib1_1P_UPSI_muon_eta;
   vector<double>  *truth_Chib1_1P_UPSI_muon_phi;
   vector<double>  *truth_Chib1_1P_pt;
   vector<double>  *truth_Chib1_1P_eta;
   vector<double>  *truth_Chib1_1P_phi;
   vector<double>  *truth_Chib1_1P_pdgid;
   vector<double>  *truth_Chib1_1P_mass;
   vector<double>  *truth_Chib1_1P_UPSI_pt;
   vector<double>  *truth_Chib1_1P_UPSI_eta;
   vector<double>  *truth_Chib1_1P_UPSI_phi;
   vector<double>  *truth_Chib1_1P_UPSI_pdgid;
   vector<double>  *truth_Chib1_1P_UPSI_mass;
   vector<double>  *truth_Chib1_1P_photon_pt;
   vector<double>  *truth_Chib1_1P_photon_eta;
   vector<double>  *truth_Chib1_1P_photon_phi;
   vector<double>  *truth_Chib1_1P_photon_pdgid;
   vector<double>  *truth_Chib2_1P_UPSI_muon_pt;
   vector<double>  *truth_Chib2_1P_UPSI_muon_eta;
   vector<double>  *truth_Chib2_1P_UPSI_muon_phi;
   vector<double>  *truth_Chib2_1P_pt;
   vector<double>  *truth_Chib2_1P_eta;
   vector<double>  *truth_Chib2_1P_phi;
   vector<double>  *truth_Chib2_1P_pdgid;
   vector<double>  *truth_Chib2_1P_mass;
   vector<double>  *truth_Chib2_1P_UPSI_pt;
   vector<double>  *truth_Chib2_1P_UPSI_eta;
   vector<double>  *truth_Chib2_1P_UPSI_phi;
   vector<double>  *truth_Chib2_1P_UPSI_pdgid;
   vector<double>  *truth_Chib2_1P_UPSI_mass;
   vector<double>  *truth_Chib2_1P_photon_pt;
   vector<double>  *truth_Chib2_1P_photon_eta;
   vector<double>  *truth_Chib2_1P_photon_phi;
   vector<double>  *truth_Chib2_1P_photon_pdgid;
   vector<double>  *loop_enter_check;
   vector<double>  *mc_event_number;
   vector<double>  *mc_run_number;
   vector<unsigned int> *mc_lumi_section;
   vector<int>     *eventHasZUpsiNTo4Mu_Count;

   // List of branches
   TBranch        *b_truth_Zmuon_pt;   //!
   TBranch        *b_truth_Zmuon_eta;   //!
   TBranch        *b_truth_Zmuon_phi;   //!
   TBranch        *b_truth_Z_pt;   //!
   TBranch        *b_truth_Z_eta;   //!
   TBranch        *b_truth_Z_phi;   //!
   TBranch        *b_truth_Z_mass;   //!
   TBranch        *b_truth_Z_pdgid;   //!
   TBranch        *b_truth_Upsimuon_pt;   //!
   TBranch        *b_truth_Upsimuon_eta;   //!
   TBranch        *b_truth_Upsimuon_phi;   //!
   TBranch        *b_truth_Upsi_pt;   //!
   TBranch        *b_truth_Upsi_eta;   //!
   TBranch        *b_truth_Upsi_phi;   //!
   TBranch        *b_truth_Upsi_mass;   //!
   TBranch        *b_truth_Upsi_pdgid;   //!
   TBranch        *b_truth_Upsi2muon_pt;   //!
   TBranch        *b_truth_Upsi2muon_eta;   //!
   TBranch        *b_truth_Upsi2muon_phi;   //!
   TBranch        *b_truth_Upsi2_pt;   //!
   TBranch        *b_truth_Upsi2_eta;   //!
   TBranch        *b_truth_Upsi2_phi;   //!
   TBranch        *b_truth_Upsi2_mass;   //!
   TBranch        *b_truth_Upsi2_pdgid;   //!
   TBranch        *b_truth_Upsi3muon_pt;   //!
   TBranch        *b_truth_Upsi3muon_eta;   //!
   TBranch        *b_truth_Upsi3muon_phi;   //!
   TBranch        *b_truth_Upsi3_pt;   //!
   TBranch        *b_truth_Upsi3_eta;   //!
   TBranch        *b_truth_Upsi3_phi;   //!
   TBranch        *b_truth_Upsi3_mass;   //!
   TBranch        *b_truth_Upsi3_pdgid;   //!
   TBranch        *b_truth_Chib0_1P_UPSI_muon_pt;   //!
   TBranch        *b_truth_Chib0_1P_UPSI_muon_eta;   //!
   TBranch        *b_truth_Chib0_1P_UPSI_muon_phi;   //!
   TBranch        *b_truth_Chib0_1P_pt;   //!
   TBranch        *b_truth_Chib0_1P_eta;   //!
   TBranch        *b_truth_Chib0_1P_phi;   //!
   TBranch        *b_truth_Chib0_1P_pdgid;   //!
   TBranch        *b_truth_Chib0_1P_mass;   //!
   TBranch        *b_truth_Chib0_1P_UPSI_pt;   //!
   TBranch        *b_truth_Chib0_1P_UPSI_phi;   //!
   TBranch        *b_truth_Chib0_1P_UPSI_eta;   //!
   TBranch        *b_truth_Chib0_1P_UPSI_pdgid;   //!
   TBranch        *b_truth_Chib0_1P_UPSI_mass;   //!
   TBranch        *b_truth_Chib0_1P_photon_pt;   //!
   TBranch        *b_truth_Chib0_1P_photon_phi;   //!
   TBranch        *b_truth_Chib0_1P_photon_eta;   //!
   TBranch        *b_truth_Chib0_1P_photon_pdgid;   //!
   TBranch        *b_truth_Chib1_1P_UPSI_muon_pt;   //!
   TBranch        *b_truth_Chib1_1P_UPSI_muon_eta;   //!
   TBranch        *b_truth_Chib1_1P_UPSI_muon_phi;   //!
   TBranch        *b_truth_Chib1_1P_pt;   //!
   TBranch        *b_truth_Chib1_1P_eta;   //!
   TBranch        *b_truth_Chib1_1P_phi;   //!
   TBranch        *b_truth_Chib1_1P_pdgid;   //!
   TBranch        *b_truth_Chib1_1P_mass;   //!
   TBranch        *b_truth_Chib1_1P_UPSI_pt;   //!
   TBranch        *b_truth_Chib1_1P_UPSI_eta;   //!
   TBranch        *b_truth_Chib1_1P_UPSI_phi;   //!
   TBranch        *b_truth_Chib1_1P_UPSI_pdgid;   //!
   TBranch        *b_truth_Chib1_1P_UPSI_mass;   //!
   TBranch        *b_truth_Chib1_1P_photon_pt;   //!
   TBranch        *b_truth_Chib1_1P_photon_eta;   //!
   TBranch        *b_truth_Chib1_1P_photon_phi;   //!
   TBranch        *b_truth_Chib1_1P_photon_pdgid;   //!
   TBranch        *b_truth_Chib2_1P_UPSI_muon_pt;   //!
   TBranch        *b_truth_Chib2_1P_UPSI_muon_eta;   //!
   TBranch        *b_truth_Chib2_1P_UPSI_muon_phi;   //!
   TBranch        *b_truth_Chib2_1P_pt;   //!
   TBranch        *b_truth_Chib2_1P_eta;   //!
   TBranch        *b_truth_Chib2_1P_phi;   //!
   TBranch        *b_truth_Chib2_1P_pdgid;   //!
   TBranch        *b_truth_Chib2_1P_mass;   //!
   TBranch        *b_truth_Chib2_1P_UPSI_pt;   //!
   TBranch        *b_truth_Chib2_1P_UPSI_eta;   //!
   TBranch        *b_truth_Chib2_1P_UPSI_phi;   //!
   TBranch        *b_truth_Chib2_1P_UPSI_pdgid;   //!
   TBranch        *b_truth_Chib2_1P_UPSI_mass;   //!
   TBranch        *b_truth_Chib2_1P_photon_pt;   //!
   TBranch        *b_truth_Chib2_1P_photon_eta;   //!
   TBranch        *b_truth_Chib2_1P_photon_phi;   //!
   TBranch        *b_truth_Chib2_1P_photon_pdgid;   //!
   TBranch        *b_loop_enter_check;   //!
   TBranch        *b_mc_event_number;   //!
   TBranch        *b_mc_run_number;   //!
   TBranch        *b_mc_lumi_section;   //!
   TBranch        *b_eventHasZUpsiNTo4Mu_Count;   //!

   treeMC(TTree *tree=0);
   virtual ~treeMC();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef treeMC_cxx
treeMC::treeMC(TTree *tree) : fChain(0) 
// {
// // if parameter tree is not specified (or zero), connect the file
// // used to generate this class and read the Tree.
//    if (tree == 0) {
//       TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("ZYto4Mu_Zto4Mu_pTCut3_MC_Bjorn_17March2022_inputFileIs_MC_DPS_2016_YZ_00623223-2B20-AB42-A456-670F9B3875D5_OR_DoubleMuTrigs.root");
//       if (!f || !f->IsOpen()) {
//          f = new TFile("ZYto4Mu_Zto4Mu_pTCut3_MC_Bjorn_17March2022_inputFileIs_MC_DPS_2016_YZ_00623223-2B20-AB42-A456-670F9B3875D5_OR_DoubleMuTrigs.root");
//       }
//       TDirectory * dir = (TDirectory*)f->Get("ZYto4Mu_Zto4Mu_pTCut3_MC_Bjorn_17March2022_inputFileIs_MC_DPS_2016_YZ_00623223-2B20-AB42-A456-670F9B3875D5_OR_DoubleMuTrigs.root:/ZmuonAnalyzer");
//       dir->GetObject("treemc",tree);
// 
//    }
//    Init(tree);
// }

{

   if (tree == 0) {
      TFile *f = root_file;
      if (!f || !f->IsOpen()) {
         f = root_file;
      }
      TDirectory * dir = (TDirectory*)f->Get("ZmuonAnalyzer");
      dir->GetObject("treemc",tree);//unclear to me if this should be tree or treemc
   }
   Init(tree);
}

treeMC::~treeMC()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t treeMC::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t treeMC::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void treeMC::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   truth_Zmuon_pt = 0;
   truth_Zmuon_eta = 0;
   truth_Zmuon_phi = 0;
   truth_Z_pt = 0;
   truth_Z_eta = 0;
   truth_Z_phi = 0;
   truth_Z_mass = 0;
   truth_Z_pdgid = 0;
   truth_Upsimuon_pt = 0;
   truth_Upsimuon_eta = 0;
   truth_Upsimuon_phi = 0;
   truth_Upsi_pt = 0;
   truth_Upsi_eta = 0;
   truth_Upsi_phi = 0;
   truth_Upsi_mass = 0;
   truth_Upsi_pdgid = 0;
   truth_Upsi2muon_pt = 0;
   truth_Upsi2muon_eta = 0;
   truth_Upsi2muon_phi = 0;
   truth_Upsi2_pt = 0;
   truth_Upsi2_eta = 0;
   truth_Upsi2_phi = 0;
   truth_Upsi2_mass = 0;
   truth_Upsi2_pdgid = 0;
   truth_Upsi3muon_pt = 0;
   truth_Upsi3muon_eta = 0;
   truth_Upsi3muon_phi = 0;
   truth_Upsi3_pt = 0;
   truth_Upsi3_eta = 0;
   truth_Upsi3_phi = 0;
   truth_Upsi3_mass = 0;
   truth_Upsi3_pdgid = 0;
   truth_Chib0_1P_UPSI_muon_pt = 0;
   truth_Chib0_1P_UPSI_muon_eta = 0;
   truth_Chib0_1P_UPSI_muon_phi = 0;
   truth_Chib0_1P_pt = 0;
   truth_Chib0_1P_eta = 0;
   truth_Chib0_1P_phi = 0;
   truth_Chib0_1P_pdgid = 0;
   truth_Chib0_1P_mass = 0;
   truth_Chib0_1P_UPSI_pt = 0;
   truth_Chib0_1P_UPSI_phi = 0;
   truth_Chib0_1P_UPSI_eta = 0;
   truth_Chib0_1P_UPSI_pdgid = 0;
   truth_Chib0_1P_UPSI_mass = 0;
   truth_Chib0_1P_photon_pt = 0;
   truth_Chib0_1P_photon_phi = 0;
   truth_Chib0_1P_photon_eta = 0;
   truth_Chib0_1P_photon_pdgid = 0;
   truth_Chib1_1P_UPSI_muon_pt = 0;
   truth_Chib1_1P_UPSI_muon_eta = 0;
   truth_Chib1_1P_UPSI_muon_phi = 0;
   truth_Chib1_1P_pt = 0;
   truth_Chib1_1P_eta = 0;
   truth_Chib1_1P_phi = 0;
   truth_Chib1_1P_pdgid = 0;
   truth_Chib1_1P_mass = 0;
   truth_Chib1_1P_UPSI_pt = 0;
   truth_Chib1_1P_UPSI_eta = 0;
   truth_Chib1_1P_UPSI_phi = 0;
   truth_Chib1_1P_UPSI_pdgid = 0;
   truth_Chib1_1P_UPSI_mass = 0;
   truth_Chib1_1P_photon_pt = 0;
   truth_Chib1_1P_photon_eta = 0;
   truth_Chib1_1P_photon_phi = 0;
   truth_Chib1_1P_photon_pdgid = 0;
   truth_Chib2_1P_UPSI_muon_pt = 0;
   truth_Chib2_1P_UPSI_muon_eta = 0;
   truth_Chib2_1P_UPSI_muon_phi = 0;
   truth_Chib2_1P_pt = 0;
   truth_Chib2_1P_eta = 0;
   truth_Chib2_1P_phi = 0;
   truth_Chib2_1P_pdgid = 0;
   truth_Chib2_1P_mass = 0;
   truth_Chib2_1P_UPSI_pt = 0;
   truth_Chib2_1P_UPSI_eta = 0;
   truth_Chib2_1P_UPSI_phi = 0;
   truth_Chib2_1P_UPSI_pdgid = 0;
   truth_Chib2_1P_UPSI_mass = 0;
   truth_Chib2_1P_photon_pt = 0;
   truth_Chib2_1P_photon_eta = 0;
   truth_Chib2_1P_photon_phi = 0;
   truth_Chib2_1P_photon_pdgid = 0;
   loop_enter_check = 0;
   mc_event_number = 0;
   mc_run_number = 0;
   mc_lumi_section = 0;
   eventHasZUpsiNTo4Mu_Count = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("truth_Zmuon_pt", &truth_Zmuon_pt, &b_truth_Zmuon_pt);
   fChain->SetBranchAddress("truth_Zmuon_eta", &truth_Zmuon_eta, &b_truth_Zmuon_eta);
   fChain->SetBranchAddress("truth_Zmuon_phi", &truth_Zmuon_phi, &b_truth_Zmuon_phi);
   fChain->SetBranchAddress("truth_Z_pt", &truth_Z_pt, &b_truth_Z_pt);
   fChain->SetBranchAddress("truth_Z_eta", &truth_Z_eta, &b_truth_Z_eta);
   fChain->SetBranchAddress("truth_Z_phi", &truth_Z_phi, &b_truth_Z_phi);
   fChain->SetBranchAddress("truth_Z_mass", &truth_Z_mass, &b_truth_Z_mass);
   fChain->SetBranchAddress("truth_Z_pdgid", &truth_Z_pdgid, &b_truth_Z_pdgid);
   fChain->SetBranchAddress("truth_Upsimuon_pt", &truth_Upsimuon_pt, &b_truth_Upsimuon_pt);
   fChain->SetBranchAddress("truth_Upsimuon_eta", &truth_Upsimuon_eta, &b_truth_Upsimuon_eta);
   fChain->SetBranchAddress("truth_Upsimuon_phi", &truth_Upsimuon_phi, &b_truth_Upsimuon_phi);
   fChain->SetBranchAddress("truth_Upsi_pt", &truth_Upsi_pt, &b_truth_Upsi_pt);
   fChain->SetBranchAddress("truth_Upsi_eta", &truth_Upsi_eta, &b_truth_Upsi_eta);
   fChain->SetBranchAddress("truth_Upsi_phi", &truth_Upsi_phi, &b_truth_Upsi_phi);
   fChain->SetBranchAddress("truth_Upsi_mass", &truth_Upsi_mass, &b_truth_Upsi_mass);
   fChain->SetBranchAddress("truth_Upsi_pdgid", &truth_Upsi_pdgid, &b_truth_Upsi_pdgid);
   fChain->SetBranchAddress("truth_Upsi2muon_pt", &truth_Upsi2muon_pt, &b_truth_Upsi2muon_pt);
   fChain->SetBranchAddress("truth_Upsi2muon_eta", &truth_Upsi2muon_eta, &b_truth_Upsi2muon_eta);
   fChain->SetBranchAddress("truth_Upsi2muon_phi", &truth_Upsi2muon_phi, &b_truth_Upsi2muon_phi);
   fChain->SetBranchAddress("truth_Upsi2_pt", &truth_Upsi2_pt, &b_truth_Upsi2_pt);
   fChain->SetBranchAddress("truth_Upsi2_eta", &truth_Upsi2_eta, &b_truth_Upsi2_eta);
   fChain->SetBranchAddress("truth_Upsi2_phi", &truth_Upsi2_phi, &b_truth_Upsi2_phi);
   fChain->SetBranchAddress("truth_Upsi2_mass", &truth_Upsi2_mass, &b_truth_Upsi2_mass);
   fChain->SetBranchAddress("truth_Upsi2_pdgid", &truth_Upsi2_pdgid, &b_truth_Upsi2_pdgid);
   fChain->SetBranchAddress("truth_Upsi3muon_pt", &truth_Upsi3muon_pt, &b_truth_Upsi3muon_pt);
   fChain->SetBranchAddress("truth_Upsi3muon_eta", &truth_Upsi3muon_eta, &b_truth_Upsi3muon_eta);
   fChain->SetBranchAddress("truth_Upsi3muon_phi", &truth_Upsi3muon_phi, &b_truth_Upsi3muon_phi);
   fChain->SetBranchAddress("truth_Upsi3_pt", &truth_Upsi3_pt, &b_truth_Upsi3_pt);
   fChain->SetBranchAddress("truth_Upsi3_eta", &truth_Upsi3_eta, &b_truth_Upsi3_eta);
   fChain->SetBranchAddress("truth_Upsi3_phi", &truth_Upsi3_phi, &b_truth_Upsi3_phi);
   fChain->SetBranchAddress("truth_Upsi3_mass", &truth_Upsi3_mass, &b_truth_Upsi3_mass);
   fChain->SetBranchAddress("truth_Upsi3_pdgid", &truth_Upsi3_pdgid, &b_truth_Upsi3_pdgid);
   fChain->SetBranchAddress("truth_Chib0_1P_UPSI_muon_pt", &truth_Chib0_1P_UPSI_muon_pt, &b_truth_Chib0_1P_UPSI_muon_pt);
   fChain->SetBranchAddress("truth_Chib0_1P_UPSI_muon_eta", &truth_Chib0_1P_UPSI_muon_eta, &b_truth_Chib0_1P_UPSI_muon_eta);
   fChain->SetBranchAddress("truth_Chib0_1P_UPSI_muon_phi", &truth_Chib0_1P_UPSI_muon_phi, &b_truth_Chib0_1P_UPSI_muon_phi);
   fChain->SetBranchAddress("truth_Chib0_1P_pt", &truth_Chib0_1P_pt, &b_truth_Chib0_1P_pt);
   fChain->SetBranchAddress("truth_Chib0_1P_eta", &truth_Chib0_1P_eta, &b_truth_Chib0_1P_eta);
   fChain->SetBranchAddress("truth_Chib0_1P_phi", &truth_Chib0_1P_phi, &b_truth_Chib0_1P_phi);
   fChain->SetBranchAddress("truth_Chib0_1P_pdgid", &truth_Chib0_1P_pdgid, &b_truth_Chib0_1P_pdgid);
   fChain->SetBranchAddress("truth_Chib0_1P_mass", &truth_Chib0_1P_mass, &b_truth_Chib0_1P_mass);
   fChain->SetBranchAddress("truth_Chib0_1P_UPSI_pt", &truth_Chib0_1P_UPSI_pt, &b_truth_Chib0_1P_UPSI_pt);
   fChain->SetBranchAddress("truth_Chib0_1P_UPSI_phi", &truth_Chib0_1P_UPSI_phi, &b_truth_Chib0_1P_UPSI_phi);
   fChain->SetBranchAddress("truth_Chib0_1P_UPSI_eta", &truth_Chib0_1P_UPSI_eta, &b_truth_Chib0_1P_UPSI_eta);
   fChain->SetBranchAddress("truth_Chib0_1P_UPSI_pdgid", &truth_Chib0_1P_UPSI_pdgid, &b_truth_Chib0_1P_UPSI_pdgid);
   fChain->SetBranchAddress("truth_Chib0_1P_UPSI_mass", &truth_Chib0_1P_UPSI_mass, &b_truth_Chib0_1P_UPSI_mass);
   fChain->SetBranchAddress("truth_Chib0_1P_photon_pt", &truth_Chib0_1P_photon_pt, &b_truth_Chib0_1P_photon_pt);
   fChain->SetBranchAddress("truth_Chib0_1P_photon_phi", &truth_Chib0_1P_photon_phi, &b_truth_Chib0_1P_photon_phi);
   fChain->SetBranchAddress("truth_Chib0_1P_photon_eta", &truth_Chib0_1P_photon_eta, &b_truth_Chib0_1P_photon_eta);
   fChain->SetBranchAddress("truth_Chib0_1P_photon_pdgid", &truth_Chib0_1P_photon_pdgid, &b_truth_Chib0_1P_photon_pdgid);
   fChain->SetBranchAddress("truth_Chib1_1P_UPSI_muon_pt", &truth_Chib1_1P_UPSI_muon_pt, &b_truth_Chib1_1P_UPSI_muon_pt);
   fChain->SetBranchAddress("truth_Chib1_1P_UPSI_muon_eta", &truth_Chib1_1P_UPSI_muon_eta, &b_truth_Chib1_1P_UPSI_muon_eta);
   fChain->SetBranchAddress("truth_Chib1_1P_UPSI_muon_phi", &truth_Chib1_1P_UPSI_muon_phi, &b_truth_Chib1_1P_UPSI_muon_phi);
   fChain->SetBranchAddress("truth_Chib1_1P_pt", &truth_Chib1_1P_pt, &b_truth_Chib1_1P_pt);
   fChain->SetBranchAddress("truth_Chib1_1P_eta", &truth_Chib1_1P_eta, &b_truth_Chib1_1P_eta);
   fChain->SetBranchAddress("truth_Chib1_1P_phi", &truth_Chib1_1P_phi, &b_truth_Chib1_1P_phi);
   fChain->SetBranchAddress("truth_Chib1_1P_pdgid", &truth_Chib1_1P_pdgid, &b_truth_Chib1_1P_pdgid);
   fChain->SetBranchAddress("truth_Chib1_1P_mass", &truth_Chib1_1P_mass, &b_truth_Chib1_1P_mass);
   fChain->SetBranchAddress("truth_Chib1_1P_UPSI_pt", &truth_Chib1_1P_UPSI_pt, &b_truth_Chib1_1P_UPSI_pt);
   fChain->SetBranchAddress("truth_Chib1_1P_UPSI_eta", &truth_Chib1_1P_UPSI_eta, &b_truth_Chib1_1P_UPSI_eta);
   fChain->SetBranchAddress("truth_Chib1_1P_UPSI_phi", &truth_Chib1_1P_UPSI_phi, &b_truth_Chib1_1P_UPSI_phi);
   fChain->SetBranchAddress("truth_Chib1_1P_UPSI_pdgid", &truth_Chib1_1P_UPSI_pdgid, &b_truth_Chib1_1P_UPSI_pdgid);
   fChain->SetBranchAddress("truth_Chib1_1P_UPSI_mass", &truth_Chib1_1P_UPSI_mass, &b_truth_Chib1_1P_UPSI_mass);
   fChain->SetBranchAddress("truth_Chib1_1P_photon_pt", &truth_Chib1_1P_photon_pt, &b_truth_Chib1_1P_photon_pt);
   fChain->SetBranchAddress("truth_Chib1_1P_photon_eta", &truth_Chib1_1P_photon_eta, &b_truth_Chib1_1P_photon_eta);
   fChain->SetBranchAddress("truth_Chib1_1P_photon_phi", &truth_Chib1_1P_photon_phi, &b_truth_Chib1_1P_photon_phi);
   fChain->SetBranchAddress("truth_Chib1_1P_photon_pdgid", &truth_Chib1_1P_photon_pdgid, &b_truth_Chib1_1P_photon_pdgid);
   fChain->SetBranchAddress("truth_Chib2_1P_UPSI_muon_pt", &truth_Chib2_1P_UPSI_muon_pt, &b_truth_Chib2_1P_UPSI_muon_pt);
   fChain->SetBranchAddress("truth_Chib2_1P_UPSI_muon_eta", &truth_Chib2_1P_UPSI_muon_eta, &b_truth_Chib2_1P_UPSI_muon_eta);
   fChain->SetBranchAddress("truth_Chib2_1P_UPSI_muon_phi", &truth_Chib2_1P_UPSI_muon_phi, &b_truth_Chib2_1P_UPSI_muon_phi);
   fChain->SetBranchAddress("truth_Chib2_1P_pt", &truth_Chib2_1P_pt, &b_truth_Chib2_1P_pt);
   fChain->SetBranchAddress("truth_Chib2_1P_eta", &truth_Chib2_1P_eta, &b_truth_Chib2_1P_eta);
   fChain->SetBranchAddress("truth_Chib2_1P_phi", &truth_Chib2_1P_phi, &b_truth_Chib2_1P_phi);
   fChain->SetBranchAddress("truth_Chib2_1P_pdgid", &truth_Chib2_1P_pdgid, &b_truth_Chib2_1P_pdgid);
   fChain->SetBranchAddress("truth_Chib2_1P_mass", &truth_Chib2_1P_mass, &b_truth_Chib2_1P_mass);
   fChain->SetBranchAddress("truth_Chib2_1P_UPSI_pt", &truth_Chib2_1P_UPSI_pt, &b_truth_Chib2_1P_UPSI_pt);
   fChain->SetBranchAddress("truth_Chib2_1P_UPSI_eta", &truth_Chib2_1P_UPSI_eta, &b_truth_Chib2_1P_UPSI_eta);
   fChain->SetBranchAddress("truth_Chib2_1P_UPSI_phi", &truth_Chib2_1P_UPSI_phi, &b_truth_Chib2_1P_UPSI_phi);
   fChain->SetBranchAddress("truth_Chib2_1P_UPSI_pdgid", &truth_Chib2_1P_UPSI_pdgid, &b_truth_Chib2_1P_UPSI_pdgid);
   fChain->SetBranchAddress("truth_Chib2_1P_UPSI_mass", &truth_Chib2_1P_UPSI_mass, &b_truth_Chib2_1P_UPSI_mass);
   fChain->SetBranchAddress("truth_Chib2_1P_photon_pt", &truth_Chib2_1P_photon_pt, &b_truth_Chib2_1P_photon_pt);
   fChain->SetBranchAddress("truth_Chib2_1P_photon_eta", &truth_Chib2_1P_photon_eta, &b_truth_Chib2_1P_photon_eta);
   fChain->SetBranchAddress("truth_Chib2_1P_photon_phi", &truth_Chib2_1P_photon_phi, &b_truth_Chib2_1P_photon_phi);
   fChain->SetBranchAddress("truth_Chib2_1P_photon_pdgid", &truth_Chib2_1P_photon_pdgid, &b_truth_Chib2_1P_photon_pdgid);
   fChain->SetBranchAddress("loop_enter_check", &loop_enter_check, &b_loop_enter_check);
   fChain->SetBranchAddress("mc_event_number", &mc_event_number, &b_mc_event_number);
   fChain->SetBranchAddress("mc_run_number", &mc_run_number, &b_mc_run_number);
   fChain->SetBranchAddress("mc_lumi_section", &mc_lumi_section, &b_mc_lumi_section);
   fChain->SetBranchAddress("eventHasZUpsiNTo4Mu_Count", &eventHasZUpsiNTo4Mu_Count, &b_eventHasZUpsiNTo4Mu_Count);
   Notify();
}

Bool_t treeMC::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void treeMC::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t treeMC::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef treeMC_cxx
