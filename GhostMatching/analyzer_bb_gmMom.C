// Studying the Spring2021 event files (Zbb): ghost matching - MC level (exclusive with exactly 2 jets [for now])
// Note: Close the event file before writing histograms to prevent seg faults
// No cuts
// studying status 71-79 to understand the charge imbalance in Zbb & Zcc samples while using GM

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
  histfname = "histZbb_partonMom.root";
  TFile *histFile = new TFile(histfname,"RECREATE");
    
  // hists for the particle loop
  TH2F* h_mom7x23      = new TH2F("h_mom7x23",     "Flavour Assigning Parton's Momentum",100,0,50,100,0,50);
  TH2F* h_mom7x23_non0 = new TH2F("h_mom7x23_non0","Flavour Assigning Parton's Momentum",100,0,50,100,0,50);
  //
  TH1F* h_mom7xb = new TH1F("h_mom7xb","Flavour Assigning Parton's Momentum - 7x (b)",100,0,50);
  TH1F* h_mom7xc = new TH1F("h_mom7xc","Flavour Assigning Parton's Momentum - 7x (c)",100,0,50);
  TH1F* h_mom7xs = new TH1F("h_mom7xs","Flavour Assigning Parton's Momentum - 7x (s)",100,0,50);
  TH1F* h_mom7xu = new TH1F("h_mom7xu","Flavour Assigning Parton's Momentum - 7x (u)",100,0,50);
  TH1F* h_mom7xd = new TH1F("h_mom7xd","Flavour Assigning Parton's Momentum - 7x (d)",100,0,50);
  TH1F* h_mom7x0 = new TH1F("h_mom7x0","Flavour Assigning Parton's Momentum - 7x (0)",100,0,50);
  //
  TH1F* h_mom7xbbar = new TH1F("h_mom7xbbar","Flavour Assigning Parton's Momentum - 7x (#bar{b})",100,0,50);
  TH1F* h_mom7xcbar = new TH1F("h_mom7xcbar","Flavour Assigning Parton's Momentum - 7x (#bar{c})",100,0,50);
  TH1F* h_mom7xsbar = new TH1F("h_mom7xsbar","Flavour Assigning Parton's Momentum - 7x (#bar{s})",100,0,50);
  TH1F* h_mom7xubar = new TH1F("h_mom7xubar","Flavour Assigning Parton's Momentum - 7x (#bar{u})",100,0,50);
  TH1F* h_mom7xdbar = new TH1F("h_mom7xdbar","Flavour Assigning Parton's Momentum - 7x (#bar{d})",100,0,50);
  
  // MC Particles
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> MCpx(tree, "MC_px");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> MCpy(tree, "MC_py");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> MCpz(tree, "MC_pz");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> MCe(tree,  "MC_e");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> MCpdg(tree,"MC_pdg");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> MCstatus(tree,"MC_status");  
  // Reco Particles
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> RPe(tree,  "RP_e");

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

  // jet counter
  unsigned int bad_jet = 0;
  
  // event loop
  while(tree.Next())
    {
      ROOT::VecOps::RVec<float> pdg_gm7x(RPe->size(),0);
      ROOT::VecOps::RVec<float> pdg_gm23(RPe->size(),0);
      TLorentzVector Zero(0,0,0,0);
      ROOT::VecOps::RVec<TLorentzVector> p4_gm7x(RPe->size(),Zero);
      ROOT::VecOps::RVec<TLorentzVector> p4_gm23(RPe->size(),Zero);
      //
      // particle loop
      //float px=0, py=0, pz=0, e=0;
      TLorentzVector p4;
      for(unsigned int ip=0; ip<MCpdg->size(); ip++)
	{
	  //px = MCpx->at(ip);
	  //py = MCpy->at(ip);
	  //pz = MCpz->at(ip);
	  //e  = MCe->at(ip);
	  p4.SetPxPyPzE(MCpx->at(ip), MCpy->at(ip), MCpz->at(ip), MCe->at(ip));
	  
	  // 71-79
	  if(MCstatus->at(ip) > 70 && MCstatus->at(ip) < 80)
	    {
	      if(MCpdg->at(ip) <= 5) //only partons
		{
		  pdg_gm7x.push_back(MCpdg->at(ip));
		  p4_gm7x.push_back(p4);
		}
	    }

	  // 23
	  if(MCstatus->at(ip) == 23)
	    {
	      if(MCpdg->at(ip) <= 5) //only partons
		{
		  pdg_gm23.push_back(MCpdg->at(ip));
		  p4_gm23.push_back(p4);
		}
	    }
	}
      
      // jet constituents
      vector<int> jet1Const7x, jet2Const7x;
      if(jetConst7x->size()>=1)      jet1Const7x = jetConst7x->at(0);
      else cout<<"**No jet constituents found in event#"<<evt+1<<"**"<<endl;
      if(jetConst7x->size()>=2)      jet2Const7x = jetConst7x->at(1);
      else cout<<"**Second jet constituents not found in event#"<<evt+1<<"**"<<endl;
      //
      vector<int> jet1Const23, jet2Const23;
      if(jetConst23->size()>=1)      jet1Const23 = jetConst23->at(0);
      else cout<<"**No jet constituents found in event#"<<evt+1<<"**"<<endl;
      if(jetConst23->size()>=2)      jet2Const23 = jetConst23->at(1);
      else cout<<"**Second jet constituents not found in event#"<<evt+1<<"**"<<endl;
      
      TLorentzVector jet1P7x, jet2P7x, jet1P23, jet2P23; //deciding parton 4-momenta
      int j1result7x=0, j2result7x=0, j1result23=0, j2result23=0;
      // Status 71-79
      // JET 1
      for(int ele : jet1Const7x) 
	{
	  if(pdg_gm7x.at(ele) == 0) continue;
	  if(abs(pdg_gm7x.at(ele)) > abs(j1result7x))
	    {
	      j1result7x = pdg_gm7x.at(ele);
	      jet1P7x    = p4_gm7x.at(ele);
	    }
	}
      // JET 2
      for(int ele : jet2Const7x) 
	{
	  if(pdg_gm7x.at(ele) == 0) continue;
	  if(abs(pdg_gm7x.at(ele)) > abs(j2result7x))
	    {
	      j2result7x = pdg_gm7x.at(ele);
	      jet2P7x    = p4_gm7x.at(ele);
	    }
	}

      // Status 23
      // JET 1
      for(int ele : jet1Const23) 
	{
	  if(pdg_gm23.at(ele) == 0) continue;
	  if(abs(pdg_gm23.at(ele)) > abs(j1result23))
	    {
	      j1result23 = pdg_gm23.at(ele);
	      jet1P23    = p4_gm23.at(ele);
	    }
	}
      // JET 2
      for(int ele : jet2Const23) 
	{
	  if(pdg_gm23.at(ele) == 0) continue;
	  if(abs(pdg_gm23.at(ele)) > abs(j2result23))
	    {
	      j2result23 = pdg_gm23.at(ele);
	      jet2P23    = p4_gm23.at(ele);
	    }
	}
    
      h_mom7x23->Fill(jet1P7x.P(), jet1P23.P());
      h_mom7x23->Fill(jet2P7x.P(), jet2P23.P());
      
      if(j1result23!=0 && j1result7x!=0) h_mom7x23_non0->Fill(jet1P7x.P(), jet1P23.P());
      if(j2result23!=0 && j2result7x!=0) h_mom7x23_non0->Fill(jet2P7x.P(), jet2P23.P());
      
      if(jet1P7x.P() > jet1P23.P()) bad_jet++;
      if(jet2P7x.P() > jet2P23.P()) bad_jet++;

      if(j1result7x == 5) h_mom7xb->Fill(jet1P7x.P());
      if(j2result7x == 5) h_mom7xb->Fill(jet2P7x.P());
      if(j1result7x == 4) h_mom7xc->Fill(jet1P7x.P());
      if(j2result7x == 4) h_mom7xc->Fill(jet2P7x.P());
      if(j1result7x == 3) h_mom7xs->Fill(jet1P7x.P());
      if(j2result7x == 3) h_mom7xs->Fill(jet2P7x.P());
      if(j1result7x == 2) h_mom7xu->Fill(jet1P7x.P());
      if(j2result7x == 2) h_mom7xu->Fill(jet2P7x.P());
      if(j1result7x == 1) h_mom7xd->Fill(jet1P7x.P());
      if(j2result7x == 1) h_mom7xd->Fill(jet2P7x.P());
      if(j1result7x == 0) h_mom7x0->Fill(jet1P7x.P());
      if(j2result7x == 0) h_mom7x0->Fill(jet2P7x.P());
      if(j1result7x == -5) h_mom7xbbar->Fill(jet1P7x.P());
      if(j2result7x == -5) h_mom7xbbar->Fill(jet2P7x.P());
      if(j1result7x == -4) h_mom7xcbar->Fill(jet1P7x.P());
      if(j2result7x == -4) h_mom7xcbar->Fill(jet2P7x.P());
      if(j1result7x == -3) h_mom7xsbar->Fill(jet1P7x.P());
      if(j2result7x == -3) h_mom7xsbar->Fill(jet2P7x.P());
      if(j1result7x == -2) h_mom7xubar->Fill(jet1P7x.P());
      if(j2result7x == -2) h_mom7xubar->Fill(jet2P7x.P());
      if(j1result7x == -1) h_mom7xdbar->Fill(jet1P7x.P());
      if(j2result7x == -1) h_mom7xdbar->Fill(jet2P7x.P());

      evt++;

      if(evt%100000==0)
	{
	  cout<<evt<<" events processed"<<endl;
	  cout<<"====================="<<endl<<endl;
	}
    }

  cout<<"There are "<<bad_jet<<" inconsistent jets (p_7x > p_23 for deciding parton)"<<endl;
  
  file->Close();
  cout<<"Event file closed"<<endl;

  h_mom7x23->Write();
  h_mom7x23_non0->Write();
  h_mom7xb->Write();
  h_mom7xc->Write();
  h_mom7xs->Write();
  h_mom7xu->Write();
  h_mom7xd->Write();
  h_mom7x0->Write();
  h_mom7xbbar->Write();
  h_mom7xcbar->Write();
  h_mom7xsbar->Write();
  h_mom7xubar->Write();
  h_mom7xdbar->Write();
  //
  histFile->Close();
  cout<<"Histograms written to file and file closed"<<endl;
  
  //delete jetConst;
  //jetConst = NULL;

  cout<<endl;
  return -1;
}
