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

  TFile *file1 = TFile::Open("p8_ee_Zbb_ecm91_gm7x_auto.root");
  TTreeReader tree1("events", file1);
  int nEvents = tree1.GetEntries();
  cout<<"Number of Events: "<<nEvents<<endl;
  //
  TFile *file2 = TFile::Open("p8_ee_Zbb_ecm91_gm_auto.root");
  TTreeReader tree2("events", file2);
  
  TString histfname;
  histfname = "histZbb_gm_chrgImb.root";
  TFile *histFile = new TFile(histfname,"RECREATE");
    
  // hists for the jet loop
  TH1F* h_p_q        = new TH1F("h_p_q",       "|p| for mismatched q jets [GeV]",      100,0,50);
  TH1F* h_p_qbar     = new TH1F("h_p_qbar",    "|p| for mismatched #bar{q} jets [GeV]",100,0,50);
  TH1F* h_theta_q    = new TH1F("h_theta_q",   "#theta for mismatched q jets",         100,0,3.15);
  TH1F* h_theta_qbar = new TH1F("h_theta_qbar","#theta for mismatched #bar{q} jets",   100,0,3.15);

  TH2F* h_flvTable   = new TH2F("h_flvTable","Flavour Table",11,-5,6,11,-5,6);
  
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

  // Table elements counter
  /*
  unsigned int dd=0, ddbar=0, du=0, dubar=0, ds=0, dsbar=0, dc=0, dcbar=0, db=0, dbbar=0;
  unsigned int dbard=0, dbardbar=0, dbaru=0, dbarubar=0, dbars=0, dbarsbar=0, dbarc=0, dbarcbar=0, dbarb=0, dbarbbar=0;
  unsigned int ud=0, udbar=0, uu=0, uubar=0, us=0, usbar=0, uc=0, ucbar=0, ub=0, ubbar=0;
  unsigned int dd=0, ddbar=0, du=0, dubar=0, ds=0, dsbar=0, dc=0, dcbar=0, db=0, dbbar=0;
  unsigned int dd=0, ddbar=0, du=0, dubar=0, ds=0, dsbar=0, dc=0, dcbar=0, db=0, dbbar=0;
  unsigned int dd=0, ddbar=0, du=0, dubar=0, ds=0, dsbar=0, dc=0, dcbar=0, db=0, dbbar=0;
  unsigned int dd=0, ddbar=0, du=0, dubar=0, ds=0, dsbar=0, dc=0, dcbar=0, db=0, dbbar=0;
  */
  // event loop
  while(tree1.Next() && tree2.Next())
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

  h_p_q->Write();
  h_p_qbar->Write();
  h_theta_q->Write();
  h_theta_qbar->Write();
  h_flvTable->Write();
  histFile->Close();
  cout<<"Histograms written to file and file closed"<<endl;
  
  //delete jetConst;
  //jetConst = NULL;

  cout<<endl;
  return -1;
}
