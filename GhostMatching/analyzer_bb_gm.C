// Studying the Spring2021 event files (Zuds files in particular): ghost matching - reco level (exclusive with exactly 2 jets [for now])
// Note: Close the event file before writing histograms to prevent seg faults
// No cuts

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
#include <TMath.h>
#include <TLorentzVector.h>
#include <TSystem.h>
#include <TInterpreter.h>

using namespace std;

int main()
{
  gInterpreter->GenerateDictionary("vector<vector<int> >","vector");

  TFile *file = TFile::Open("p8_ee_Zbb_ecm91_gm.root");
  TTreeReader tree("events", file);
  int nEvents = tree.GetEntries();
  cout<<"Number of Events: "<<nEvents<<endl;

  TString histfname;
  histfname = "histZbb_gm.root";
  TFile *histFile = new TFile(histfname,"RECREATE");
    
  // hists for the jet loop
  TH1F* h_jetFlavour = new TH1F("h_jetFlavour","Jet Flavour",11,-5,6);
  
  // reco particles                                                       
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> RPpx(tree, "RP_px");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> RPpy(tree, "RP_py");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> RPpz(tree, "RP_pz");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> RPe(tree, "RP_e");
  
  // jets and jet constituents - eekt     
  TTreeReaderValue<vector<vector<int>>> jetConst(tree, "jetconstituents_ee_kt");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> jetPx(tree, "jets_ee_kt_px");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> jetPy(tree, "jets_ee_kt_py");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> jetPz(tree, "jets_ee_kt_pz");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> jetE(tree, "jets_ee_kt_e");
  TTreeReaderValue<vector<int,ROOT::Detail::VecOps::RAdoptAllocator<int>>> jetFlavour(tree, "jets_ee_kt_flavour");

  // event counter
  int evt = 0;

  // flavour mismatch counter
  int flv_mm = 0;
  
  // event loop
  while(tree.Next())
    {
      // jet loop
      int nJet = jetE->size();
      
      for(unsigned int j=0; j<nJet; j++)
	{
	  // jet flavour
	  h_jetFlavour->Fill(jetFlavour->at(j));
	}

      if(abs(jetFlavour->at(0)) != abs(jetFlavour->at(1))) flv_mm++;
      //if(jetFlavour->at(0) == 0 || jetFlavour->at(1) == 0) flv_mm++;
      //if(jetFlavour->at(0) == 0) flv_mm += jetFlavour->at(1)/jetFlavour->at(1);
      //if(jetFlavour->at(1) == 0) flv_mm += jetFlavour->at(0)/jetFlavour->at(0);
      
      /*======================*/
      
      // jet constituents
      vector<int> jet1Const, jet2Const;
      if(jetConst->size()>=1)      jet1Const = jetConst->at(0);
      else cout<<"**No jet constituents found in event#"<<evt+1<<"**"<<endl;
      if(jetConst->size()>=2)      jet2Const = jetConst->at(1);
      else cout<<"**Second jet constituents not found in event#"<<evt+1<<"**"<<endl;
    
      // JET 1
      //cout<<"total reco particles  : "<<RPe->size()<<endl;
      //cout<<"total jet constituents: "<<jet1Const.size()+jet2Const.size()<<endl;
      //cout<<"Jet1 indices: ";
      for(int ele : jet1Const) 
	{
	  //cout<<ele<<", ";
	}
      //cout<<endl;

      // JET 2
      //cout<<"Jet2 indices: ";
      for(int ele : jet2Const) 
	{
	  //cout<<ele<<", ";
	}
      //cout<<endl;
      
      jet1Const.clear();
      jet2Const.clear();

      evt++;

      if(evt%50000==0)
	{
	  cout<<evt<<" events processed"<<endl;
	  //cout<<"Event #"<<evt<<":"<<endl;
	  //cout<<"This event has "<<RPe->size()<<" particles and "<<nJet<<" jets."<<endl;
	  cout<<"====================="<<endl<<endl;
	}
      //cout<<"====================="<<endl<<endl;
    }

  cout<<flv_mm<<" events have jets with differently assigned flavours"<<endl;
  
  file->Close();
  cout<<"Event file closed"<<endl;

  h_jetFlavour->Write();
  histFile->Close();
  cout<<"Histograms written to file and file closed"<<endl;
  
  //delete jetConst;
  //jetConst = NULL;

  cout<<endl;
  return -1;
}
