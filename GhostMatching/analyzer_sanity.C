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
  histfname = "histZbb_sanity.root";
  TFile *histFile = new TFile(histfname,"RECREATE");
    
  // hists for the particle loop
  TH1F* h_momDiff = new TH1F("h_momDiff","|p|_{7x} - |p|_{23}",100,-15,15);
  
  // MC Particles
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> MCpx(tree, "MC_px");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> MCpy(tree, "MC_py");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> MCpz(tree, "MC_pz");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> MCe(tree,  "MC_e");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> MCpdg(tree,"MC_pdg");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> MCstatus(tree,"MC_status");  
  // Reco Particles
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> RPe7x(tree,  "RP_e");
  /*
  // MC Particles
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> MCpx(tree1, "MC_px");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> MCpy(tree1, "MC_py");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> MCpz(tree1, "MC_pz");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> MCe(tree1,  "MC_e");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> MCpdg(tree1,"MC_pdg");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> MCstatus(tree1,"MC_status");  
  */
  // Reco Particles
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> RPe23(tree,  "RP_e");

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
  unsigned int bad_evt = 0;
  unsigned int normal = 0, inverse = 0;
  
  // event loop
  while(tree.Next())
    {
      if(RPe7x->size() != RPe23->size()) cout<<"Different samples in event#"<<evt<<endl;

      
      if(jetE7x->size() != jetE23->size()) cout<<"Something wrong with event#"<<evt<<endl;
      TLorentzVector pJ7x[jetE7x->size()], pJ23[jetE7x->size()];
      for(unsigned int j=0; j<jetE7x->size(); j++)
	{
	  pJ7x[j].SetPxPyPzE(jetPx7x->at(j), jetPy7x->at(j), jetPz7x->at(j), jetE7x->at(j));
	  pJ23[j].SetPxPyPzE(jetPx23->at(j), jetPy23->at(j), jetPz23->at(j), jetE23->at(j));
	  //cout<<"p_jet"<<j+1<<"-> 7x: "<<pJ7x[j].P()<<", 23: "<<pJ23[j].P()<<endl;
	}
      //cout<<"-------"<<endl;

      if(pJ7x[0].P()==pJ23[0].P() && pJ7x[1].P()==pJ23[1].P())
	{
	  normal++;
	  h_momDiff->Fill(pJ7x[0].P() - pJ23[0].P());
	  h_momDiff->Fill(pJ7x[1].P() - pJ23[1].P());
	}

      if(pJ7x[0].P()==pJ23[1].P() && pJ7x[1].P()==pJ23[0].P())
	{
	  inverse++;
	  h_momDiff->Fill(pJ7x[0].P() - pJ23[1].P());
	  h_momDiff->Fill(pJ7x[1].P() - pJ23[0].P());
	}

      if(pJ7x[0].P()!=pJ23[0].P() || pJ7x[1].P()!=pJ23[1].P())
	{
	  cout<<"For event#"<<evt<<endl;
	  cout<<"Jet 1 -> 7x: "<<pJ7x[0].P()<<", 23: "<<pJ23[0].P()<<endl;
	  cout<<"Jet 2 -> 7x: "<<pJ7x[1].P()<<", 23: "<<pJ23[1].P()<<endl;
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
      
      int constSize23 = jet1Const23.size()+jet2Const23.size();
      int constSize7x = jet1Const7x.size()+jet2Const7x.size();
      if(constSize23 > constSize7x)
	{
	  //cout<<"Event#"<<evt<<"sucks!"<<endl;
	  bad_evt++;
	}
      
      evt++;

      if(evt%100000==0)
	{
	  cout<<evt<<" events processed"<<endl;
	  cout<<"====================="<<endl<<endl;
	}
    }

  cout<<"There are "<<bad_evt<<" bad events in the sample"<<endl;
  cout<<normal<<" events with normal ordering and "<<inverse<<" inverted"<<endl;
  
  file->Close();
  cout<<"Event files closed"<<endl;

  h_momDiff->Write();
  histFile->Close();
  cout<<"Histograms written to file and file closed"<<endl;
  
  //delete jetConst;
  //jetConst = NULL;

  cout<<endl;
  return -1;
}
