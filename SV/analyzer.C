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

  //TFile *file = TFile::Open("p8_ee_Zbb_ecm91_V0_10K.root");
  //TFile *file = TFile::Open("p8_ee_Zcc_ecm91_V0_10K.root");
  TFile *file = TFile::Open("p8_ee_Zuds_ecm91_V0_10K.root");
  TTreeReader tree("events", file);
  int nEvents = tree.GetEntries();
  cout<<"Number of Events: "<<nEvents<<endl;
  
  TString histfname;
  //histfname = "histZbb_V0.root";
  //histfname = "histZcc_V0.root";
  histfname = "histZuds_V0.root";
  TFile *histFile = new TFile(histfname,"RECREATE");

  // hists for invariant mass
  TH1F* h_invM_V0      = new TH1F("h_invM_V0","V^{0} Invariant Mass",500,0.4,1.2);
  TH1F* h_invM_Ks      = new TH1F("h_invM_Ks","K_{S} Invariant Mass (MC)",100,0,1);
  TH1F* h_invM_Lambda0 = new TH1F("h_invM_Lambda0","#Lambda^{0} Invariant Mass (MC)",100,0.5,1.5);

  TH1F* h_mass_Ks      = new TH1F("h_mass_Ks","K_{S} MC Mass",100,0,1);
  TH1F* h_mass_Lambda0 = new TH1F("h_mass_Lambda0","#Lambda^{0} MC Mass",100,0.5,1.5);
  
  // hists for Ks
  //TH1F* h_nKs_MC    = new TH1F("h_nKs_MC",   "K_{S} Multiplicity from MC Truth",              10,0,10);
  //TH1F* h_nKs_V0    = new TH1F("h_nKs_V0",   "K_{S} Multiplicity from V^{0} Reconstruction",  10,0,10);
  TH1F* h_nKs_diff  = new TH1F("h_nKs_diff", "K_{S} Multiplicity Difference (MC - V^{0})",    13,-5,8);
  TH1F* h_invMKs_V0 = new TH1F("h_invMKs_V0","K_{S} Invariant Mass from V^{0} Reconstruction",100,0.493,0.503);
  TH1F* h_dKs_V0    = new TH1F("h_dKs_V0",   "K_{S} Distance",                                100,0,0.001);
  
  // hists for Lambda0
  //TH1F* h_nLambda_MC    = new TH1F("h_nLambda_MC",   "#Lambda^{0} Multiplicity from MC Truth",              10,0,10);
  //TH1F* h_nLambda_V0    = new TH1F("h_nLambda_V0",   "#Lambda^{0} Multiplicity from V^{0} Reconstruction",  10,0,10);
  TH1F* h_nLambda_diff  = new TH1F("h_nLambda_diff", "#Lambda^{0} Multiplicity Difference (MC - V^{0})",    10,-5,5);
  TH1F* h_invMLambda_V0 = new TH1F("h_invMLambda_V0","#Lambda^{0} Invariant Mass from V^{0} Reconstruction",100,1.111,1.121);
  TH1F* h_dLambda_V0    = new TH1F("h_dLambda_V0",   "#Lambda^{0} Distance",                                100,0,0.001);
  
  // MC particles                                                       
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> MCpdg(tree, "MC_pdg");
  TTreeReaderValue<vector<TLorentzVector,ROOT::Detail::VecOps::RAdoptAllocator<TLorentzVector>>> MCp4(tree, "MC_p4");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> MCmass(tree, "MC_mass");

  // Ks->pi+pi-
  TTreeReaderValue<vector<int,ROOT::Detail::VecOps::RAdoptAllocator<int>>> Ks2pipi(tree, "K0spipi_indices");

  // Lambda0->p+pi-
  TTreeReaderValue<vector<int,ROOT::Detail::VecOps::RAdoptAllocator<int>>> Lambda02ppi(tree, "Lambda0ppi_indices");

  // V0
  TTreeReaderValue<vector<int,ROOT::Detail::VecOps::RAdoptAllocator<int>>> V0pdg(tree, "V0_pdg");
  TTreeReaderValue<vector<double,ROOT::Detail::VecOps::RAdoptAllocator<double>>> V0invM(tree, "V0_invM");
  TTreeReaderValue<vector<TVector3,ROOT::Detail::VecOps::RAdoptAllocator<TVector3>>> V0pos(tree, "V0_pos");

  // event counter
  int evt = 0;

  int nKs_MC = 0, nLambda0_MC = 0;
  int nKs_V0 = 0, nLambda0_V0 = 0;
  int nKs2pipi = 0, nLambda02ppi = 0;

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
	      h_invM_Ks->Fill(MCp4->at(ctr).M());
	      //
	      h_mass_Ks->Fill(MCmass->at(ctr));
	    }
	  // Lambda0 - MC
	  else if(MCpdg->at(ctr)==3122)
	    {
	      nLambda0_MC++;
	      nLambda0_MC_evt++;
	      h_invM_Lambda0->Fill(MCp4->at(ctr).M());
	      //
	      h_mass_Lambda0->Fill(MCmass->at(ctr));
	    }
	}

      // Ks -> pi+pi- loop
      nKs2pipi += Ks2pipi->size()/3;
      //for(unsigned int ctr=0; ctr<Ks2pipi->size(); ctr++)

      // Lambda0 -> p+pi- loop
      nLambda02ppi += Lambda02ppi->size()/3;
      //for(unsigned int ctr=0; ctr<Ks2pipi->size(); ctr++)

      // V0 loop
      for(unsigned int ctr=0; ctr<V0pdg->size(); ctr++)
	{
	  TVector3 pos(V0pos->at(ctr));

	  h_invM_V0->Fill(V0invM->at(ctr));
	  
	  // Ks - V0
	  if(V0pdg->at(ctr)==310)
	    {
	      nKs_V0++;
	      nKs_V0_evt++;
	      //
	      h_invMKs_V0->Fill(V0invM->at(ctr));
	      //
	      h_dKs_V0->Fill(pos.Mag());
	    }
	  // Lambda0 - V0
	  else if(V0pdg->at(ctr)==3122)
	    {
	      nLambda0_V0++;
	      nLambda0_V0_evt++;
	      //
	      h_invMLambda_V0->Fill(V0invM->at(ctr));
	      //
	      h_dLambda_V0->Fill(pos.Mag());
	    }
	}

      h_nKs_diff->Fill(nKs_MC_evt - nKs_V0_evt);
      h_nLambda_diff->Fill(nLambda0_MC_evt - nLambda0_V0_evt);
	  
      evt++;

      if(evt%10000==0)
	{
	  //cout<<evt<<" events processed"<<endl;
	  cout<<"Event #"<<evt<<":"<<endl;
	  cout<<"This event has "<<nKs_MC_evt<<" K-shorts("<<Ks2pipi->size()/3<<" of them decay to pi+pi-) & "<<nLambda0_MC_evt<<" Lambda0s("<<Lambda02ppi->size()/3<<" of them decay to p+pi-)"<<endl;
	  cout<<"out of which "<<nKs_V0_evt<<" K-shorts & "<<nLambda0_V0_evt<<" Lambda0s were reconstructed"<<endl;
	  cout<<"====================="<<endl<<endl;
	}
    }

  cout<<"There are "<<nKs_MC<<" K-shorts("<<nKs2pipi<<" of them decay to pi+pi-) & "<<nLambda0_MC<<" Lambda0s("<<nLambda02ppi<<" of them decay to pi+pi-) in "<<nEvents<<" events"<<endl;
  cout<<"out of which "<<nKs_V0<<" K-shorts & "<<nLambda0_V0<<" Lambda0s were reconstructed"<<endl;

  file->Close();
  cout<<"Event file closed"<<endl;

  h_invM_V0->Write();
  h_invM_Ks->Write();
  h_invM_Lambda0->Write();
  h_mass_Ks->Write();
  h_mass_Lambda0->Write();
  h_nKs_diff->Write();
  h_nLambda_diff->Write();
  h_invMKs_V0->Write();
  h_invMLambda_V0->Write();
  h_dKs_V0->Write();
  h_dLambda_V0->Write();
  histFile->Close();
  cout<<"Histograms written to file and file closed"<<endl;

  return -1;
}
