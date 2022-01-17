// Studying the Spring2021 event files (Zbb): ghost matching - reco level (exclusive with exactly 2 jets [for now])
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

  TFile *file = TFile::Open("p8_ee_Zbb_ecm91.root");
  TTreeReader tree("events", file);
  int nEvents = tree.GetEntries();
  cout<<"Number of Events: "<<nEvents<<endl;
  
  TString histfname;
  histfname = "histZbb_status.root";
  TFile *histFile = new TFile(histfname,"RECREATE");
    
  // hists for the particle loop
  TH1F* h_flv7x        = new TH1F("h_flv7x",       "Multiplicity in Status 71-79 (No Gluons)",11,-5,6);
  TH1F* h_flv7x_GE     = new TH1F("h_flv7x_GE",    "q#bar{q} in Status 71-79",6,0,6);
  TH1F* h_flv7x_GEsign = new TH1F("h_flv7x_GEsign","Jet Flavour in Status 71-79",11,-5,6);
  
  // All MC particles                                                         
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> MCpx(tree, "MC_px");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> MCpy(tree, "MC_py");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> MCpz(tree, "MC_pz");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> MCe(tree, "MC_e");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> MCpdg(tree, "MC_pdg");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> MCstatus(tree, "MC_status");

  // event counter
  unsigned int evt = 0;
  
  // event loop
  while(tree.Next())
    {
      int result=0;
      vector<float> pdg_gm;
      // particle loop
      for(unsigned int ip=0; ip<MCpdg->size(); ip++)
	{
	  if(MCstatus->at(ip) < 70 || MCstatus->at(ip) > 80) continue; // 71-79
	  if(abs(MCpdg->at(ip)) == 21) continue; // remove gluons
	  h_flv7x->Fill(MCpdg->at(ip));

	  pdg_gm.push_back(MCpdg->at(ip));

	  if(abs(MCpdg->at(ip)) == 5) continue;  // remove b bbar
	  if(abs(MCpdg->at(ip)) > abs(result)) result = MCpdg->at(ip);
	}

      h_flv7x_GEsign->Fill(result);
      
      if(find(pdg_gm.begin(), pdg_gm.end(), 4) != pdg_gm.end() && find(pdg_gm.begin(), pdg_gm.end(), -4) != pdg_gm.end()) h_flv7x_GE->Fill(4);
      else if(find(pdg_gm.begin(), pdg_gm.end(), 3) != pdg_gm.end() && find(pdg_gm.begin(), pdg_gm.end(), -3) != pdg_gm.end()) h_flv7x_GE->Fill(3);
      else if(find(pdg_gm.begin(), pdg_gm.end(), 2) != pdg_gm.end() && find(pdg_gm.begin(), pdg_gm.end(), -2) != pdg_gm.end()) h_flv7x_GE->Fill(2);
      else if(find(pdg_gm.begin(), pdg_gm.end(), 1) != pdg_gm.end() && find(pdg_gm.begin(), pdg_gm.end(), -1) != pdg_gm.end()) h_flv7x_GE->Fill(1);
      
      evt++;

      if(evt%100000==0)
	{
	  cout<<evt<<" events processed"<<endl;
	  cout<<"====================="<<endl<<endl;
	}
    }

  file->Close();
  cout<<"Event files closed"<<endl;

  h_flv7x->Write();
  h_flv7x_GE->Write();
  h_flv7x_GEsign->Write();
  histFile->Close();
  cout<<"Histograms written to file and file closed"<<endl;
  
  //delete jetConst;
  //jetConst = NULL;

  cout<<endl;
  return -1;
}
