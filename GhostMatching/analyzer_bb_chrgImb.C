// Studying the Spring2021 event files (Zbb): ghost matching - reco level (exclusive with exactly 2 jets [for now])
// Note: Close the event file before writing histograms to prevent seg faults
// No cuts
// status 71-79 for parton selection
// studying the charge imbalance in Zbb & Zcc samples while using GM

#include <iostream>
#include <cmath>
#include <vector>

#include <TObject.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include "ROOT/RVec.hxx"
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TSystem.h>
#include <TInterpreter.h>

using namespace std;

int main()
{
  gInterpreter->GenerateDictionary("vector<vector<int> >","vector");

  TFile *file = TFile::Open("p8_ee_Zbb_ecm91_gm7x23_auto.root");
  TTreeReader tree("events", file);
  int nEvents = tree.GetEntries();
  cout<<"Number of Events: "<<nEvents<<endl;
  
  TString histfname;
  histfname = "histZbb_gm_chrgImb.root";
  TFile *histFile = new TFile(histfname,"RECREATE");
    
  // hists for the jet loop
  TH1F* h_p_q        = new TH1F("h_p_q",       "|p| for mismatched q jets [GeV]",      100,0,50);
  TH1F* h_p_qbar     = new TH1F("h_p_qbar",    "|p| for mismatched #bar{q} jets [GeV]",100,0,50);
  TH1F* h_theta_q    = new TH1F("h_theta_q",   "#theta for mismatched q jets",         100,0,3.15);
  TH1F* h_theta_qbar = new TH1F("h_theta_qbar","#theta for mismatched #bar{q} jets",   100,0,3.15);

  TH2F* h_flvTable   = new TH2F("h_flvTable","Flavour Table",11,-5,6,11,-5,6);

  // 23 -> +5 && 7x -> +1
  TH1F* h_p_p5p1     = new TH1F("h_p_p5p1",    "|p| (23->5 && 7x->1) [GeV]", 50,0,50);
  TH1F* h_theta_p5p1 = new TH1F("h_theta_p5p1","#theta (23->5 && 7x->1)",    50,0,3.15);
  // 23 -> +5 && 7x -> -1
  TH1F* h_p_p5m1     = new TH1F("h_p_p5m1",    "|p| (23->5 && 7x->-1) [GeV]",50,0,50);
  TH1F* h_theta_p5m1 = new TH1F("h_theta_p5m1","#theta (23->5 && 7x->-1)",   50,0,3.15);
  // 23 -> +5 && 7x -> +2
  TH1F* h_p_p5p2     = new TH1F("h_p_p5p2",    "|p| (23->5 && 7x->2) [GeV]", 50,0,50);
  TH1F* h_theta_p5p2 = new TH1F("h_theta_p5p2","#theta (23->5 && 7x->2)",    50,0,3.15);
  // 23 -> +5 && 7x -> -2
  TH1F* h_p_p5m2     = new TH1F("h_p_p5m2",    "|p| (23->5 && 7x->-2) [GeV]",50,0,50);
  TH1F* h_theta_p5m2 = new TH1F("h_theta_p5m2","#theta (23->5 && 7x->-2)",   50,0,3.15);
  // 23 -> +5 && 7x -> +3
  TH1F* h_p_p5p3     = new TH1F("h_p_p5p3",    "|p| (23->5 && 7x->3) [GeV]", 50,0,50);
  TH1F* h_theta_p5p3 = new TH1F("h_theta_p5p3","#theta (23->5 && 7x->3)",    50,0,3.15);
  // 23 -> +5 && 7x -> -3
  TH1F* h_p_p5m3     = new TH1F("h_p_p5m3",    "|p| (23->5 && 7x->-3) [GeV]",50,0,50);
  TH1F* h_theta_p5m3 = new TH1F("h_theta_p5m3","#theta (23->5 && 7x->-3)",   50,0,3.15);
  // 23 -> +5 && 7x -> +4
  TH1F* h_p_p5p4     = new TH1F("h_p_p5p4",    "|p| (23->5 && 7x->4) [GeV]", 50,0,50);
  TH1F* h_theta_p5p4 = new TH1F("h_theta_p5p4","#theta (23->5 && 7x->4)",    50,0,3.15);
  // 23 -> +5 && 7x -> -4
  TH1F* h_p_p5m4     = new TH1F("h_p_p5m4",    "|p| (23->5 && 7x->-4) [GeV]",50,0,50);
  TH1F* h_theta_p5m4 = new TH1F("h_theta_p5m4","#theta (23->5 && 7x->-4)",   50,0,3.15);
  // 23 -> +5 && 7x -> +5
  TH1F* h_p_p5p5     = new TH1F("h_p_p5p5",    "|p| (23->5 && 7x->5) [GeV]", 50,0,50);
  TH1F* h_theta_p5p5 = new TH1F("h_theta_p5p5","#theta (23->5 && 7x->5)",    50,0,3.15);
  // 23 -> +5 && 7x -> -5
  TH1F* h_p_p5m5     = new TH1F("h_p_p5m5",    "|p| (23->5 && 7x->-5) [GeV]",50,0,50);
  TH1F* h_theta_p5m5 = new TH1F("h_theta_p5m5","#theta (23->5 && 7x->-5)",   50,0,3.15);
  //
  // 23 -> -5 && 7x -> +1
  TH1F* h_p_m5p1     = new TH1F("h_p_m5p1",    "|p| (23->-5 && 7x->1) [GeV]", 50,0,50);
  TH1F* h_theta_m5p1 = new TH1F("h_theta_m5p1","#theta (23->-5 && 7x->1)",    50,0,3.15);
  // 23 -> -5 && 7x -> -1
  TH1F* h_p_m5m1     = new TH1F("h_p_m5m1",    "|p| (23->-5 && 7x->-1) [GeV]",50,0,50);
  TH1F* h_theta_m5m1 = new TH1F("h_theta_m5m1","#theta (23->-5 && 7x->-1)",   50,0,3.15);
  // 23 -> -5 && 7x -> +2
  TH1F* h_p_m5p2     = new TH1F("h_p_m5p2",    "|p| (23->-5 && 7x->2) [GeV]", 50,0,50);
  TH1F* h_theta_m5p2 = new TH1F("h_theta_m5p2","#theta (23->-5 && 7x->2)",    50,0,3.15);
  // 23 -> -5 && 7x -> -2
  TH1F* h_p_m5m2     = new TH1F("h_p_m5m2",    "|p| (23->-5 && 7x->-2) [GeV]",50,0,50);
  TH1F* h_theta_m5m2 = new TH1F("h_theta_m5m2","#theta (23->-5 && 7x->-2)",   50,0,3.15);
  // 23 -> -5 && 7x -> +3
  TH1F* h_p_m5p3     = new TH1F("h_p_m5p3",    "|p| (23->-5 && 7x->3) [GeV]", 50,0,50);
  TH1F* h_theta_m5p3 = new TH1F("h_theta_m5p3","#theta (23->-5 && 7x->3)",    50,0,3.15);
  // 23 -> -5 && 7x -> -3
  TH1F* h_p_m5m3     = new TH1F("h_p_m5m3",    "|p| (23->-5 && 7x->-3) [GeV]",50,0,50);
  TH1F* h_theta_m5m3 = new TH1F("h_theta_m5m3","#theta (23->-5 && 7x->-3)",   50,0,3.15);
  // 23 -> -5 && 7x -> +4
  TH1F* h_p_m5p4     = new TH1F("h_p_m5p4",    "|p| (23->-5 && 7x->4) [GeV]", 50,0,50);
  TH1F* h_theta_m5p4 = new TH1F("h_theta_m5p4","#theta (23->-5 && 7x->4)",    50,0,3.15);
  // 23 -> -5 && 7x -> -4
  TH1F* h_p_m5m4     = new TH1F("h_p_m5m4",    "|p| (23->-5 && 7x->-4) [GeV]",50,0,50);
  TH1F* h_theta_m5m4 = new TH1F("h_theta_m5m4","#theta (23->-5 && 7x->-4)",   50,0,3.15);
  // 23 -> -5 && 7x -> +5
  TH1F* h_p_m5p5     = new TH1F("h_p_m5p5",    "|p| (23->-5 && 7x->5) [GeV]", 50,0,50);
  TH1F* h_theta_m5p5 = new TH1F("h_theta_m5p5","#theta (23->-5 && 7x->5)",    50,0,3.15);
  // 23 -> -5 && 7x -> -5
  TH1F* h_p_m5m5     = new TH1F("h_p_m5m5",    "|p| (23->-5 && 7x->-5) [GeV]",50,0,50);
  TH1F* h_theta_m5m5 = new TH1F("h_theta_m5m5","#theta (23->-5 && 7x->-5)",   50,0,3.15);
  //
  // 23 -> 0 && 7x -> +1
  TH1F* h_p_0p1     = new TH1F("h_p_0p1",    "|p| (23->0 && 7x->1) [GeV]", 50,0,50);
  TH1F* h_theta_0p1 = new TH1F("h_theta_0p1","#theta (23->0 && 7x->1)",    50,0,3.15);
  // 23 -> 0 && 7x -> -1
  TH1F* h_p_0m1     = new TH1F("h_p_0m1",    "|p| (23->0 && 7x->-1) [GeV]",50,0,50);
  TH1F* h_theta_0m1 = new TH1F("h_theta_0m1","#theta (23->0 && 7x->-1)",   50,0,3.15);
  // 23 -> 0 && 7x -> +2
  TH1F* h_p_0p2     = new TH1F("h_p_0p2",    "|p| (23->0 && 7x->2) [GeV]", 50,0,50);
  TH1F* h_theta_0p2 = new TH1F("h_theta_0p2","#theta (23->0 && 7x->2)",    50,0,3.15);
  // 23 -> 0 && 7x -> -2
  TH1F* h_p_0m2     = new TH1F("h_p_0m2",    "|p| (23->0 && 7x->-2) [GeV]",50,0,50);
  TH1F* h_theta_0m2 = new TH1F("h_theta_0m2","#theta (23->0 && 7x->-2)",   50,0,3.15);
  // 23 -> 0 && 7x -> +3
  TH1F* h_p_0p3     = new TH1F("h_p_0p3",    "|p| (23->0 && 7x->3) [GeV]", 50,0,50);
  TH1F* h_theta_0p3 = new TH1F("h_theta_0p3","#theta (23->0 && 7x->3)",    50,0,3.15);
  // 23 -> 0 && 7x -> -3
  TH1F* h_p_0m3     = new TH1F("h_p_0m3",    "|p| (23->0 && 7x->-3) [GeV]",50,0,50);
  TH1F* h_theta_0m3 = new TH1F("h_theta_0m3","#theta (23->0 && 7x->-3)",   50,0,3.15);
  // 23 -> 0 && 7x -> +4
  TH1F* h_p_0p4     = new TH1F("h_p_0p4",    "|p| (23->0 && 7x->4) [GeV]", 50,0,50);
  TH1F* h_theta_0p4 = new TH1F("h_theta_0p4","#theta (23->0 && 7x->4)",    50,0,3.15);
  // 23 -> 0 && 7x -> -4
  TH1F* h_p_0m4     = new TH1F("h_p_0m4",    "|p| (23->0 && 7x->-4) [GeV]",50,0,50);
  TH1F* h_theta_0m4 = new TH1F("h_theta_0m4","#theta (23->0 && 7x->-4)",   50,0,3.15);
  // 23 -> 0 && 7x -> +5
  TH1F* h_p_0p5     = new TH1F("h_p_0p5",    "|p| (23->0 && 7x->5) [GeV]", 50,0,50);
  TH1F* h_theta_0p5 = new TH1F("h_theta_0p5","#theta (23->0 && 7x->5)",    50,0,3.15);
  // 23 -> 0 && 7x -> -5
  TH1F* h_p_0m5     = new TH1F("h_p_0m5",    "|p| (23->0 && 7x->-5) [GeV]",50,0,50);
  TH1F* h_theta_0m5 = new TH1F("h_theta_0m5","#theta (23->0 && 7x->-5)",   50,0,3.15);
  //
  // 23 -> +5 && 7x -> 0
  TH1F* h_p_p50     = new TH1F("h_p_p50",    "|p| (23->5 && 7x->0) [GeV]", 50,0,50);
  TH1F* h_theta_p50 = new TH1F("h_theta_p50","#theta (23->5 && 7x->0)",    50,0,3.15);
  // 23 -> -5 && 7x -> 0
  TH1F* h_p_m50     = new TH1F("h_p_m50",    "|p| (23->-5 && 7x->0) [GeV]", 50,0,50);
  TH1F* h_theta_m50 = new TH1F("h_theta_m50","#theta (23->-5 && 7x->0)",    50,0,3.15);
  // 23 -> 0 && 7x -> 0
  TH1F* h_p_00      = new TH1F("h_p_00",     "|p| (23->0 && 7x->0) [GeV]", 50,0,50);
  TH1F* h_theta_00  = new TH1F("h_theta_00", "#theta (23->0 && 7x->0)",    50,0,3.15);

  //
  
  // Jets (Status 71-79)     
  TTreeReaderValue<vector<vector<int>>> jetConst7x(tree, "jetconstituents_ee_kt7x");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> jetPx7x(tree, "jets_ee_kt_px7x");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> jetPy7x(tree, "jets_ee_kt_py7x");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> jetPz7x(tree, "jets_ee_kt_pz7x");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> jetE7x(tree,  "jets_ee_kt_e7x");
  TTreeReaderValue<vector<int,ROOT::Detail::VecOps::RAdoptAllocator<int>>>     jetFlavour7x(tree,"jets_ee_kt_flavour7x");

  // Jets (Status 23)     
  TTreeReaderValue<vector<vector<int>>> jetConst23(tree, "jetconstituents_ee_kt23");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> jetPx23(tree, "jets_ee_kt_px23");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> jetPy23(tree, "jets_ee_kt_py23");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> jetPz23(tree, "jets_ee_kt_pz23");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> jetE23(tree,  "jets_ee_kt_e23");
  TTreeReaderValue<vector<int,ROOT::Detail::VecOps::RAdoptAllocator<int>>>     jetFlavour23(tree,"jets_ee_kt_flavour23");

  // event counter
  unsigned int evt = 0;

  // event loop
  while(tree.Next())
    {
      // jet loop
      if(jetE7x->size() != jetE23->size()) cout<<"Something's wrong with event no. "<<evt+1<<endl;
      int nJet = jetE23->size();
      float Px7x=0, Py7x=0, Pz7x=0, E7x=0;
      float Px23=0, Py23=0, Pz23=0, E23=0;
      TLorentzVector p_Jet7x[nJet], p_Jet23[nJet];
      for(unsigned int j=0; j<nJet; j++)
	{
	  Px7x = jetPx7x->at(j);
	  Py7x = jetPy7x->at(j);
	  Pz7x = jetPz7x->at(j);
	  E7x  = jetE7x->at(j);

	  Px23 = jetPx23->at(j);
	  Py23 = jetPy23->at(j);
	  Pz23 = jetPz23->at(j);
	  E23  = jetE23->at(j);
	  
	  p_Jet7x[j].SetPxPyPzE(Px7x,Py7x,Pz7x,E7x);
	  p_Jet23[j].SetPxPyPzE(Px23,Py23,Pz23,E23);

	  h_flvTable->Fill(jetFlavour23->at(j),jetFlavour7x->at(j));

	  // 23 -> +5 && 7x -> +1
	  if(jetFlavour23->at(j) == 5 && jetFlavour7x->at(j) == 1)
	    {
	      h_p_p5p1->Fill(p_Jet7x[j].P());
	      h_theta_p5p1->Fill(p_Jet7x[j].Theta());
	    }

	  // 23 -> +5 && 7x -> -1
	  if(jetFlavour23->at(j) == 5 && jetFlavour7x->at(j) == -1)
	    {
	      h_p_p5m1->Fill(p_Jet7x[j].P());
	      h_theta_p5m1->Fill(p_Jet7x[j].Theta());
	    }

	  // 23 -> +5 && 7x -> +2
	  if(jetFlavour23->at(j) == 5 && jetFlavour7x->at(j) == 2)
	    {
	      h_p_p5p2->Fill(p_Jet7x[j].P());
	      h_theta_p5p2->Fill(p_Jet7x[j].Theta());
	    }

	  // 23 -> +5 && 7x -> -2
	  if(jetFlavour23->at(j) == 5 && jetFlavour7x->at(j) == -2)
	    {
	      h_p_p5m2->Fill(p_Jet7x[j].P());
	      h_theta_p5m2->Fill(p_Jet7x[j].Theta());
	    }

	  // 23 -> +5 && 7x -> +3
	  if(jetFlavour23->at(j) == 5 && jetFlavour7x->at(j) == 3)
	    {
	      h_p_p5p3->Fill(p_Jet7x[j].P());
	      h_theta_p5p3->Fill(p_Jet7x[j].Theta());
	    }

	  // 23 -> +5 && 7x -> -3
	  if(jetFlavour23->at(j) == 5 && jetFlavour7x->at(j) == -3)
	    {
	      h_p_p5m3->Fill(p_Jet7x[j].P());
	      h_theta_p5m3->Fill(p_Jet7x[j].Theta());
	    }

	  // 23 -> +5 && 7x -> +4
	  if(jetFlavour23->at(j) == 5 && jetFlavour7x->at(j) == 4)
	    {
	      h_p_p5p4->Fill(p_Jet7x[j].P());
	      h_theta_p5p4->Fill(p_Jet7x[j].Theta());
	    }

	  // 23 -> +5 && 7x -> -4
	  if(jetFlavour23->at(j) == 5 && jetFlavour7x->at(j) == -4)
	    {
	      h_p_p5m4->Fill(p_Jet7x[j].P());
	      h_theta_p5m4->Fill(p_Jet7x[j].Theta());
	    }

	  // 23 -> +5 && 7x -> +5
	  if(jetFlavour23->at(j) == 5 && jetFlavour7x->at(j) == 5)
	    {
	      h_p_p5p5->Fill(p_Jet7x[j].P());
	      h_theta_p5p5->Fill(p_Jet7x[j].Theta());
	    }

	  // 23 -> +5 && 7x -> -5
	  if(jetFlavour23->at(j) == 5 && jetFlavour7x->at(j) == -5)
	    {
	      h_p_p5m5->Fill(p_Jet7x[j].P());
	      h_theta_p5m5->Fill(p_Jet7x[j].Theta());
	    }

	  //
	  
	  // 23 -> -5 && 7x -> +1
	  if(jetFlavour23->at(j) == -5 && jetFlavour7x->at(j) == 1)
	    {
	      h_p_m5p1->Fill(p_Jet7x[j].P());
	      h_theta_m5p1->Fill(p_Jet7x[j].Theta());
	    }

	  // 23 -> -5 && 7x -> -1
	  if(jetFlavour23->at(j) == -5 && jetFlavour7x->at(j) == -1)
	    {
	      h_p_m5m1->Fill(p_Jet7x[j].P());
	      h_theta_m5m1->Fill(p_Jet7x[j].Theta());
	    }

	  // 23 -> -5 && 7x -> +2
	  if(jetFlavour23->at(j) == -5 && jetFlavour7x->at(j) == 2)
	    {
	      h_p_m5p2->Fill(p_Jet7x[j].P());
	      h_theta_m5p2->Fill(p_Jet7x[j].Theta());
	    }

	  // 23 -> -5 && 7x -> -2
	  if(jetFlavour23->at(j) == -5 && jetFlavour7x->at(j) == -2)
	    {
	      h_p_m5m2->Fill(p_Jet7x[j].P());
	      h_theta_m5m2->Fill(p_Jet7x[j].Theta());
	    }

	  // 23 -> -5 && 7x -> +3
	  if(jetFlavour23->at(j) == -5 && jetFlavour7x->at(j) == 3)
	    {
	      h_p_m5p3->Fill(p_Jet7x[j].P());
	      h_theta_m5p3->Fill(p_Jet7x[j].Theta());
	    }

	  // 23 -> -5 && 7x -> -3
	  if(jetFlavour23->at(j) == -5 && jetFlavour7x->at(j) == -3)
	    {
	      h_p_m5m3->Fill(p_Jet7x[j].P());
	      h_theta_m5m3->Fill(p_Jet7x[j].Theta());
	    }

	  // 23 -> -5 && 7x -> +4
	  if(jetFlavour23->at(j) == -5 && jetFlavour7x->at(j) == 4)
	    {
	      h_p_m5p4->Fill(p_Jet7x[j].P());
	      h_theta_m5p4->Fill(p_Jet7x[j].Theta());
	    }

	  // 23 -> -5 && 7x -> -4
	  if(jetFlavour23->at(j) == -5 && jetFlavour7x->at(j) == -4)
	    {
	      h_p_m5m4->Fill(p_Jet7x[j].P());
	      h_theta_m5m4->Fill(p_Jet7x[j].Theta());
	    }

	  // 23 -> -5 && 7x -> +5
	  if(jetFlavour23->at(j) == -5 && jetFlavour7x->at(j) == 5)
	    {
	      h_p_m5p5->Fill(p_Jet7x[j].P());
	      h_theta_m5p5->Fill(p_Jet7x[j].Theta());
	    }

	  // 23 -> -5 && 7x -> -5
	  if(jetFlavour23->at(j) == -5 && jetFlavour7x->at(j) == -5)
	    {
	      h_p_m5m5->Fill(p_Jet7x[j].P());
	      h_theta_m5m5->Fill(p_Jet7x[j].Theta());
	    }

	  //
	  
	  // 23 -> 0 && 7x -> +1
	  if(jetFlavour23->at(j) == 0 && jetFlavour7x->at(j) == 1)
	    {
	      h_p_0p1->Fill(p_Jet7x[j].P());
	      h_theta_0p1->Fill(p_Jet7x[j].Theta());
	    }

	  // 23 -> 0 && 7x -> -1
	  if(jetFlavour23->at(j) == 0 && jetFlavour7x->at(j) == -1)
	    {
	      h_p_0m1->Fill(p_Jet7x[j].P());
	      h_theta_0m1->Fill(p_Jet7x[j].Theta());
	    }

	  // 23 -> 0 && 7x -> +2
	  if(jetFlavour23->at(j) == 0 && jetFlavour7x->at(j) == 2)
	    {
	      h_p_0p2->Fill(p_Jet7x[j].P());
	      h_theta_0p2->Fill(p_Jet7x[j].Theta());
	    }

	  // 23 -> 0 && 7x -> -2
	  if(jetFlavour23->at(j) == 0 && jetFlavour7x->at(j) == -2)
	    {
	      h_p_0m2->Fill(p_Jet7x[j].P());
	      h_theta_0m2->Fill(p_Jet7x[j].Theta());
	    }

	  // 23 -> 0 && 7x -> +3
	  if(jetFlavour23->at(j) == 0 && jetFlavour7x->at(j) == 3)
	    {
	      h_p_0p3->Fill(p_Jet7x[j].P());
	      h_theta_0p3->Fill(p_Jet7x[j].Theta());
	    }

	  // 23 -> 0 && 7x -> -3
	  if(jetFlavour23->at(j) == 0 && jetFlavour7x->at(j) == -3)
	    {
	      h_p_0m3->Fill(p_Jet7x[j].P());
	      h_theta_0m3->Fill(p_Jet7x[j].Theta());
	    }

	  // 23 -> 0 && 7x -> +4
	  if(jetFlavour23->at(j) == 0 && jetFlavour7x->at(j) == 4)
	    {
	      h_p_0p4->Fill(p_Jet7x[j].P());
	      h_theta_0p4->Fill(p_Jet7x[j].Theta());
	    }

	  // 23 -> 0 && 7x -> -4
	  if(jetFlavour23->at(j) == 0 && jetFlavour7x->at(j) == -4)
	    {
	      h_p_0m4->Fill(p_Jet7x[j].P());
	      h_theta_0m4->Fill(p_Jet7x[j].Theta());
	    }

	  // 23 -> 0 && 7x -> +5
	  if(jetFlavour23->at(j) == 0 && jetFlavour7x->at(j) == 5)
	    {
	      h_p_0p5->Fill(p_Jet7x[j].P());
	      h_theta_0p5->Fill(p_Jet7x[j].Theta());
	    }

	  // 23 -> 0 && 7x -> -5
	  if(jetFlavour23->at(j) == 0 && jetFlavour7x->at(j) == -5)
	    {
	      h_p_0m5->Fill(p_Jet7x[j].P());
	      h_theta_0m5->Fill(p_Jet7x[j].Theta());
	    }

	  //

	  // 23 -> +5 && 7x -> 0
	  if(jetFlavour23->at(j) == 5 && jetFlavour7x->at(j) == 0)
	    {
	      h_p_p50->Fill(p_Jet7x[j].P());
	      h_theta_p50->Fill(p_Jet7x[j].Theta());
	    }

	  // 23 -> -5 && 7x -> 0
	  if(jetFlavour23->at(j) == -5 && jetFlavour7x->at(j) == 0)
	    {
	      h_p_m50->Fill(p_Jet7x[j].P());
	      h_theta_m50->Fill(p_Jet7x[j].Theta());
	    }
	  
	  // 23 -> 0 && 7x -> 0
	  if(jetFlavour23->at(j) == 0 && jetFlavour7x->at(j) == 0)
	    {
	      h_p_00->Fill(p_Jet7x[j].P());
	      h_theta_00->Fill(p_Jet7x[j].Theta());
	    }	}
	  
      evt++;

      if(evt%100000==0)
	{
	  cout<<evt<<" events processed"<<endl;
	  cout<<"====================="<<endl<<endl;
	}
    }

  file->Close();
  cout<<"Event file closed"<<endl;

  h_p_q->Write();
  h_p_qbar->Write();
  h_theta_q->Write();
  h_theta_qbar->Write();
  h_flvTable->Write();
  ///
  h_p_p5p1->Write();
  h_theta_p5p1->Write();
  h_p_p5m1->Write();
  h_theta_p5m1->Write();
  h_p_p5p2->Write();
  h_theta_p5p2->Write();
  h_p_p5m2->Write();
  h_theta_p5m2->Write();
  h_p_p5p3->Write();
  h_theta_p5p3->Write();
  h_p_p5m3->Write();
  h_theta_p5m3->Write();
  h_p_p5p4->Write();
  h_theta_p5p4->Write();
  h_p_p5m4->Write();
  h_theta_p5m4->Write();
  h_p_p5p5->Write();
  h_theta_p5p5->Write();
  h_p_p5m5->Write();
  h_theta_p5m5->Write();
  //
  h_p_m5p1->Write();
  h_theta_m5p1->Write();
  h_p_m5m1->Write();
  h_theta_m5m1->Write();
  h_p_m5p2->Write();
  h_theta_m5p2->Write();
  h_p_m5m2->Write();
  h_theta_m5m2->Write();
  h_p_m5p3->Write();
  h_theta_m5p3->Write();
  h_p_m5m3->Write();
  h_theta_m5m3->Write();
  h_p_m5p4->Write();
  h_theta_m5p4->Write();
  h_p_m5m4->Write();
  h_theta_m5m4->Write();
  h_p_m5p5->Write();
  h_theta_m5p5->Write();
  h_p_m5m5->Write();
  h_theta_m5m5->Write();
  //
  h_p_0p1->Write();
  h_theta_0p1->Write();
  h_p_0m1->Write();
  h_theta_0m1->Write();
  h_p_0p2->Write();
  h_theta_0p2->Write();
  h_p_0m2->Write();
  h_theta_0m2->Write();
  h_p_0p3->Write();
  h_theta_0p3->Write();
  h_p_0m3->Write();
  h_theta_0m3->Write();
  h_p_0p4->Write();
  h_theta_0p4->Write();
  h_p_0m4->Write();
  h_theta_0m4->Write();
  h_p_0p5->Write();
  h_theta_0p5->Write();
  h_p_0m5->Write();
  h_theta_0m5->Write();
  //
  h_p_p50->Write();
  h_theta_p50->Write();
  h_p_m50->Write();
  h_theta_m50->Write();
  h_p_00->Write();
  h_theta_00->Write();
  
  histFile->Close();
  cout<<"Histograms written to file and file closed"<<endl;
  
  //delete jetConst;
  //jetConst = NULL;

  cout<<endl;
  return -1;
}
