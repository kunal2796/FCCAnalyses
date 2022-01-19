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

  TFile *file1 = TFile::Open("p8_ee_Zbb_ecm91_gm7x_auto.root");
  TTreeReader tree1("events", file1);
  int nEvents = tree1.GetEntries();
  cout<<"Number of Events: "<<nEvents<<endl;
  //
  TFile *file2 = TFile::Open("p8_ee_Zbb_ecm91_gm_auto.root");
  TTreeReader tree2("events", file2);
  
  TString histfname;
  histfname = "histZbb_partonMom.root";
  TFile *histFile = new TFile(histfname,"RECREATE");
    
  // hists for the particle loop
  //TH1F* h_flv7x = new TH1F("h_flv7x","Multiplicity in Status 71-79 (No Gluons)",11,-5,6);
  
  // MC Particles
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> MCpx(tree1, "MC_px");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> MCpy(tree1, "MC_py");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> MCpz(tree1, "MC_pz");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> MCe(tree1,  "MC_e");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> MCpdg(tree1,"MC_pdg");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> MCstatus(tree1,"MC_status");  
  // Reco Particles
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> RPe(tree1,  "RP_e");

  // Jets (Status 71-79)     
  TTreeReaderValue<vector<vector<int>>> jetConst7x(tree1, "jetconstituents_ee_kt");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> jetPx7x(tree1, "jets_ee_kt_px");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> jetPy7x(tree1, "jets_ee_kt_py");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> jetPz7x(tree1, "jets_ee_kt_pz");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> jetE7x(tree1,  "jets_ee_kt_e");
  TTreeReaderValue<vector<int,ROOT::Detail::VecOps::RAdoptAllocator<int>>>     jetFlavour7x(tree1,"jets_ee_kt_flavour");

  // Jets (Status 23)     
  TTreeReaderValue<vector<vector<int>>> jetConst23(tree2, "jetconstituents_ee_kt");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> jetPx23(tree2, "jets_ee_kt_px");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> jetPy23(tree2, "jets_ee_kt_py");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> jetPz23(tree2, "jets_ee_kt_pz");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> jetE23(tree2,  "jets_ee_kt_e");
  TTreeReaderValue<vector<int,ROOT::Detail::VecOps::RAdoptAllocator<int>>>     jetFlavour23(tree2,"jets_ee_kt_flavour");

  // event counter
  unsigned int evt = 0;
  
  // event loop
  while(tree1.Next() && tree2.Next())
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
	      if(MCpdg->at(ip) > 5) continue; //only partons
	      pdg_gm7x.push_back(MCpdg->at(ip));
	      p4_gm7x.push_back(p4);
	    }

	  // 23
	  if(MCstatus->at(ip) == 23)
	    {
	      //if(MCpdg->at(ip) > 5) continue; //only partons
	      pdg_gm23.push_back(MCpdg->at(ip));
	      p4_gm23.push_back(p4);
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
      
      if(evt == 5224)
	{
	  cout<<"Reco size = "<<RPe->size()<<endl;
	  cout<<"pdg_gm7x size = "<<pdg_gm7x.size()<<endl;
	  cout<<"pdg_gm23 size = "<<pdg_gm23.size()<<endl;
	  cout<<"jet1Const7x size = "<<jet1Const7x.size()<<endl;
	  cout<<"jet2Const7x size = "<<jet2Const7x.size()<<endl;
	  cout<<"jet1Const23 size = "<<jet1Const23.size()<<endl;
	  cout<<"jet2Const23 size = "<<jet2Const23.size()<<endl;
	}
      
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

      cout<<"check 5 in "<<evt+1<<endl;
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
      cout<<"check 6"<<endl;
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

      evt++;

      if(evt%100000==0)
	{
	  cout<<evt<<" events processed"<<endl;
	  cout<<"====================="<<endl<<endl;
	}
    }

  file1->Close();
  file2->Close();
  cout<<"Event files closed"<<endl;

  //  h_flv7x->Write();
  histFile->Close();
  cout<<"Histograms written to file and file closed"<<endl;
  
  //delete jetConst;
  //jetConst = NULL;

  cout<<endl;
  return -1;
}
