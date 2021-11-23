// Studying the Spring2021 event files (Zuds files in particular): jet flavour addignment - reco level (exclusive with exactly 2 jets)
// Note: Close the event file before writing histograms to prevent seg faults
// No cuts
// COULDN'T RECOMPILE FCCANALYSIS WITH THE NEW DEFINITION BECAUSE OF SOME ERRORS IN VERTEXFINDERACTS

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

  TFile *file = TFile::Open("p8_ee_Zuds_ecm91_reco_st23.root");
  TTreeReader tree("events", file);
  int nEvents = tree.GetEntries();
  cout<<"Number of Events: "<<nEvents<<endl;

  TString histfname;
  histfname = "histZuds_excl_st23.root";
  TFile *histFile = new TFile(histfname,"RECREATE");
  
  // hists for the jet loop
  TH1F* h_jetFlavour = new TH1F("h_jetFlavour","Jet Flavour [status 23]",6,0,6);
  TH1F* h_jetFlavour_qqbar = new TH1F("h_jetFlavour_qqbar","Jet Flavour (q-#bar{q} separation) [status 23]",7,-3,4);

  // reco particles                                                       
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> RPpx(tree, "RP_px");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> RPpy(tree, "RP_py");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> RPpz(tree, "RP_pz");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> RPe(tree, "RP_e");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> RPp(tree, "RP_p");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> RPtheta(tree, "RP_theta");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> RPmass(tree, "RP_mass");
  
  // jets and jet constituents - eekt     
  TTreeReaderValue<vector<vector<int>>> jetConst(tree, "jetconstituents_ee_kt");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> jetPx(tree, "jets_ee_kt_px");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> jetPy(tree, "jets_ee_kt_py");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> jetPz(tree, "jets_ee_kt_pz");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> jetE(tree, "jets_ee_kt_e");
  TTreeReaderValue<vector<int,ROOT::Detail::VecOps::RAdoptAllocator<int>>> jetFlavour(tree, "jets_ee_kt_flavour"); //not in the ntuple yet

  // jets and jet constituents - kt     
  //TTreeReaderValue<vector<vector<int>>> jetConst(tree, "jetconstituents_kt");
  //TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> jetPx(tree, "jets_kt_px");
  //TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> jetPy(tree, "jets_kt_py");
  //TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> jetPz(tree, "jets_kt_pz");
  //TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> jetE(tree, "jets_kt_e");
  //TTreeReaderValue<vector<int,ROOT::Detail::VecOps::RAdoptAllocator<int>>> jetFlavour(tree, "jets_kt_flavour");

  // event counter
  int evt = 0;

  // event loop
  while(tree.Next())
    {
      // jet loop
      for(unsigned int j=0; j<jetE->size(); j++)
	{
	  // jet flavour
	  h_jetFlavour->Fill(abs(jetFlavour->at(j)));

	  // jet flavour - qqbar separation
	  h_jetFlavour_qqbar->Fill(jetFlavour->at(j));
	}

      /*======================*/
      
      // jet constituents
      /*
      vector<int> jet1Const, jet2Const;
      if(jetConst->size()>=1)      jet1Const = jetConst->at(0);
      else cout<<"**No jet constituents found in event#"<<evt+1<<"**"<<endl;
      if(jetConst->size()>=2)      jet2Const = jetConst->at(1);
      else cout<<"**Second jet constituents not found in event#"<<evt+1<<"**"<<endl;
      */
      evt++;

      if(evt%50000==0)
	{
	  cout<<evt<<" events processed"<<endl;
	  //cout<<"====================="<<endl<<endl;
	}
    }

  file->Close();
  cout<<"Event file closed"<<endl;

  h_jetFlavour->Write();
  h_jetFlavour_qqbar->Write();
  histFile->Close();
  cout<<"Histograms written to file and file closed"<<endl;

  //delete jetConst;
  //jetConst = NULL;

  cout<<endl;
  return -1;
}
