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

  // cc, uu, dd, ss event counter
  unsigned int ccEvt=0, uuEvt=0, ddEvt=0, ssEvt=0, Evt00=0;
  
  // event loop
  while(tree1.Next() && tree2.Next())
  //for(unsigned int evt=0; evt<nEvents; evt++)
    {
      //tree1.Next();
      //tree2.Next();
      /*
      if(abs(jetFlavour7x->at(0)) == 4 && abs(jetFlavour7x->at(1)) == 4) ccEvt++;
      if(abs(jetFlavour7x->at(0)) == 3 && abs(jetFlavour7x->at(1)) == 3) ssEvt++;
      if(abs(jetFlavour7x->at(0)) == 2 && abs(jetFlavour7x->at(1)) == 2) uuEvt++;
      if(abs(jetFlavour7x->at(0)) == 1 && abs(jetFlavour7x->at(1)) == 1) ddEvt++;
      */
      if(     jetFlavour7x->at(0) == 4 && jetFlavour7x->at(1) == -4) ccEvt++;
      else if(jetFlavour7x->at(0) == -4 && jetFlavour7x->at(1) == 4) ccEvt++;
      if(     jetFlavour7x->at(0) == 3 && jetFlavour7x->at(1) == -3) ssEvt++;
      else if(jetFlavour7x->at(0) == -3 && jetFlavour7x->at(1) == 3) ssEvt++;
      if(     jetFlavour7x->at(0) == 2 && jetFlavour7x->at(1) == -2) uuEvt++;
      else if(jetFlavour7x->at(0) == -2 && jetFlavour7x->at(1) == 2) uuEvt++;
      if(     jetFlavour7x->at(0) == 1 && jetFlavour7x->at(1) == -1) ddEvt++;
      else if(jetFlavour7x->at(0) == -1 && jetFlavour7x->at(1) == 1) ddEvt++;
      if(     jetFlavour7x->at(0) == 0 && jetFlavour7x->at(1) ==  0) Evt00++;
      
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

	  // case 1: unmatched jet-flavours (a: all, b: non b's)
	  //if(/*jetFlavour23->at(j) != 0 &&*/ jetFlavour7x->at(j) != jetFlavour23->at(j))
	  // case 2: all non b's
	  if(abs(jetFlavour7x->at(j)) != 5)
	  // case 3: unmatched jet-flavour for [a: differently, b: similarly] flavoured jets
	  //if(abs(jetFlavour7x->at(0)) != abs(jetFlavour7x->at(1)) && jetFlavour7x->at(j) != jetFlavour23->at(j))
	    {
	      //if(abs(jetFlavour7x->at(j)) == 5) continue;
	  
	      // momentum
	      if(jetFlavour7x->at(j) > 0) h_p_q->Fill(p_Jet7x[j].P());    // q jets
	      if(jetFlavour7x->at(j) < 0) h_p_qbar->Fill(p_Jet7x[j].P()); // qbar jets
	      
	      // polar angle
	      if(jetFlavour7x->at(j) > 0) h_theta_q->Fill(p_Jet7x[j].Theta());    // q jets
	      if(jetFlavour7x->at(j) < 0) h_theta_qbar->Fill(p_Jet7x[j].Theta()); // qbar jets	      
	    }
	}
      
      evt++;

      if(evt%100000==0)
	{
	  cout<<evt<<" events processed"<<endl;
	  cout<<"====================="<<endl<<endl;
	}
    }

  cout<<"cc events = "<<ccEvt<<endl;
  cout<<"ss events = "<<ssEvt<<endl;
  cout<<"uu events = "<<uuEvt<<endl;
  cout<<"dd events = "<<ddEvt<<endl;
  cout<<"00 events = "<<Evt00<<endl;
  
  file1->Close();
  file2->Close();
  cout<<"Event files closed"<<endl;

  h_p_q->Write();
  h_p_qbar->Write();
  h_theta_q->Write();
  h_theta_qbar->Write();
  histFile->Close();
  cout<<"Histograms written to file and file closed"<<endl;
  
  //delete jetConst;
  //jetConst = NULL;

  cout<<endl;
  return -1;
}
