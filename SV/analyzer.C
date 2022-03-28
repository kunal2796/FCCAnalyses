// Studying the Spring2021 event files (Zuds files in particular): V0 reconstruction
// Note: Close the event file before writing histograms to prevent seg faults

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
  //gInterpreter->GenerateDictionary("vector<vector<int> >","vector");

  TFile *file = TFile::Open("p8_ee_Zuds_ecm91_V0_100K.root");
  TTreeReader tree("events", file);
  int nEvents = tree.GetEntries();
  cout<<"Number of Events: "<<nEvents<<endl;
  /*
  TString histfname;
  histfname = "histZuds_V0.root";
  TFile *histFile = new TFile(histfname,"RECREATE");
  
  // hists for Ks
  TH1F* h_nKs_MC = new TH1F("h_nKs_MC","K_{S} Multiplicity from MC Truth",10,0,10);
  TH1F* h_nKs_V0 = new TH1F("h_nKs_V0","K_{S} Multiplicity from V^{0} Reconstruction",10,0,10);
  TH1F* h_nKs_diff = new TH1F("h_nKs_diff","K_{S} Multiplicity Difference (MC - V^{0})",10,0,10);

  // hists for Lambda0
  TH1F* h_nLambda_MC = new TH1F("h_nLambda_MC","#Lambda^{0} Multiplicity from MC Truth",10,0,10);
  TH1F* h_nLambda_V0 = new TH1F("h_nLambda_V0","#Lambda^{0} Multiplicity from V^{0} Reconstruction",10,0,10);
  TH1F* h_nLambda_diff = new TH1F("h_nLambda_diff","#Lambda^{0} Multiplicity Difference (MC - V^{0})",10,0,10);
  */
  // MC particles                                                       
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> MCpdg(tree, "MC_pdg");
  
  // V0
  //TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> MCpdg(tree, "MC_pdg");

  // event counter
  int evt = 0;

  int nKs = 0, nLambda0 = 0;

  // event loop
  while(tree.Next())
    {
      int nKs_evt = 0, nLambda0_evt = 0;
      
      // MC particle loop
      for(unsigned int ctr=0; ctr<MCpdg->size(); ctr++)
	{
	  if(MCpdg->at(ctr)==310)
	    {
	      nKs++;
	      nKs_evt++;
	    }

	  if(MCpdg->at(ctr)==3122)
	    {
	      nLambda0++;
	      nLambda0_evt++;
	    }
	}
	  
      evt++;

      if(evt%10000==0)
	{
	  //cout<<evt<<" events processed"<<endl;
	  cout<<"Event #"<<evt<<":"<<endl;
	  cout<<"This event has "<<nKs_evt<<" K-shorts & "<<nLambda0_evt<<" Lambda0s"<<endl;
	  cout<<"====================="<<endl<<endl;
	}
    }

  cout<<"There are "<<nKs<<" K-shorts & "<<nLambda0<<" Lambda0s in "<<nEvents<<" events"<<endl;

  file->Close();
  cout<<"Event file closed"<<endl;

  //cout<<"Histograms written to file and file closed"<<endl;

  return -1;
  cout<<"return is not the problem"<<endl;
}
