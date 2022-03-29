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
#include <TVector3.h>
#include <TSystem.h>
#include <TInterpreter.h>

using namespace std;

int main()
{
  //gInterpreter->GenerateDictionary("vector<vector<int> >","vector");

  TFile *file = TFile::Open("p8_ee_Zuds_ecm91_V0_10K.root");
  TTreeReader tree("events", file);
  int nEvents = tree.GetEntries();
  cout<<"Number of Events: "<<nEvents<<endl;
  
  TString histfname;
  histfname = "histZuds_V0.root";
  TFile *histFile = new TFile(histfname,"RECREATE");
  
  // hists for Ks
  TH1F* h_nKs_MC = new TH1F("h_nKs_MC","K_{S} Multiplicity from MC Truth",10,0,10);
  TH1F* h_nKs_V0 = new TH1F("h_nKs_V0","K_{S} Multiplicity from V^{0} Reconstruction",10,0,10);
  TH1F* h_nKs_diff = new TH1F("h_nKs_diff","K_{S} Multiplicity Difference (MC - V^{0})",13,-5,8);

  // hists for Lambda0
  TH1F* h_nLambda_MC = new TH1F("h_nLambda_MC","#Lambda^{0} Multiplicity from MC Truth",10,0,10);
  TH1F* h_nLambda_V0 = new TH1F("h_nLambda_V0","#Lambda^{0} Multiplicity from V^{0} Reconstruction",10,0,10);
  TH1F* h_nLambda_diff = new TH1F("h_nLambda_diff","#Lambda^{0} Multiplicity Difference (MC - V^{0})",10,-5,5);
  
  // MC particles                                                       
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> MCpdg(tree, "MC_pdg");
  
  // V0
  TTreeReaderValue<vector<int,ROOT::Detail::VecOps::RAdoptAllocator<int>>> V0pdg(tree, "V0_pdg");
  TTreeReaderValue<vector<double,ROOT::Detail::VecOps::RAdoptAllocator<double>>> V0invM(tree, "V0_invM");

  // event counter
  int evt = 0;

  int nKs_MC = 0, nLambda0_MC = 0;
  int nKs_V0 = 0, nLambda0_V0 = 0;

  // event loop
  while(tree.Next())
    {
      int nKs_MC_evt = 0, nLambda0_MC_evt = 0;
      int nKs_V0_evt = 0, nLambda0_V0_evt = 0;
      
      // MC particle loop
      for(unsigned int ctr=0; ctr<MCpdg->size(); ctr++)
	{
	  // Ks - MC
	  if(MCpdg->at(ctr)==310)
	    {
	      nKs_MC++;
	      nKs_MC_evt++;
	    }
	  // Lambda0 - MC
	  if(MCpdg->at(ctr)==3122)
	    {
	      nLambda0_MC++;
	      nLambda0_MC_evt++;
	    }
	}

      for(unsigned int ctr=0; ctr<V0pdg->size(); ctr++)
	{
	  // Ks - V0
	  if(V0pdg->at(ctr)==310)
	    {
	      nKs_V0++;
	      nKs_V0_evt++;
	    }
	  // Lambda0 - V0
	  if(V0pdg->at(ctr)==3122)
	    {
	      nLambda0_V0++;
	      nLambda0_V0_evt++;
	    }
	}

      h_nKs_diff->Fill(nKs_MC_evt - nKs_V0_evt);
      h_nLambda_diff->Fill(nLambda0_MC_evt - nLambda0_V0_evt);
	  
      evt++;

      if(evt%10000==0)
	{
	  //cout<<evt<<" events processed"<<endl;
	  cout<<"Event #"<<evt<<":"<<endl;
	  cout<<"This event has "<<nKs_MC_evt<<" K-shorts & "<<nLambda0_MC_evt<<" Lambda0s"<<endl;
	  cout<<"out of which "<<nKs_V0_evt<<" K-shorts & "<<nLambda0_V0_evt<<" Lambda0s were reconstructed"<<endl;
	  cout<<"====================="<<endl<<endl;
	}
    }

  cout<<"There are "<<nKs_MC<<" K-shorts & "<<nLambda0_MC<<" Lambda0s in "<<nEvents<<" events"<<endl;
  cout<<"out of which "<<nKs_V0<<" K-shorts & "<<nLambda0_V0<<" Lambda0s were reconstructed"<<endl;

  file->Close();
  cout<<"Event file closed"<<endl;

  h_nKs_diff->Write();
  h_nLambda_diff->Write();
  histFile->Close();
  cout<<"Histograms written to file and file closed"<<endl;

  return -1;
}
