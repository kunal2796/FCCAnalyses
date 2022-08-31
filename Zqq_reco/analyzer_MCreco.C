// DIDN'T WORK
// Studying the Spring2021 event files (Zuds files in particular): association between tracks and MC final particles
// Note: Close the event file before writing histograms to prevent seg faults
// No cuts?
// Remember to edit once FCCAna is fixed and ntuples with IP sig can be made
// why is the no of MC final particles not the same as reco particles??

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

  TFile *file = TFile::Open("p8_ee_Zuds_ecm91_MCreco.root");
  TTreeReader tree("events", file);
  int nEvents = tree.GetEntries();
  cout<<"Number of Events: "<<nEvents<<endl;

  TString histfname;
  histfname = "histZuds_MCreco.root";
  TFile *histFile = new TFile(histfname,"RECREATE");

  TH1F* h_chrg = new TH1F("h_chrg","MC Final Particles Charge",20,-2,2);
  TH1F* h_PV = new TH1F("h_PV","production vtx - primary tracks",100,0,2);
  TH1F* h_nPV = new TH1F("h_nPV","production vtx - non primary tracks",100,0,2);

  // MC particles                                                       
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> MCpx(tree, "MC_px_f");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> MCpy(tree, "MC_py_f");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> MCpz(tree, "MC_pz_f");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> MCe(tree, "MC_e_f");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> MCx(tree, "MC_x_f");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> MCy(tree, "MC_y_f");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> MCz(tree, "MC_z_f");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> MCcharge(tree, "MC_charge_f");
  
  // reco particles                                                       
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> RPpx(tree, "RP_px");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> RPpy(tree, "RP_py");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> RPpz(tree, "RP_pz");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> RPe(tree, "RP_e");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> RPcharge(tree, "RP_charge");

  // track parameters
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> RPD0(tree, "RP_D0");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> RPZ0(tree, "RP_Z0");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> RPphi(tree, "RP_phi");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> RPomega(tree, "RP_omega");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> RPtanLambda(tree, "RP_tanLambda");

  // IP significance
  //TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> RPD0sig(tree, "RP_D0_sig");
  //TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> RPZ0sig(tree, "RP_Z0_sig");

  // IP covariance (until FCCAna fixed)
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> RPD0cov(tree, "RP_D0_cov");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> RPZ0cov(tree, "RP_Z0_cov");
  
  // event counter
  int evt = 0;

  int Cmc_Crp = 0;
  
  // event loop
  while(tree.Next())
    {
      // MC particle loop
      int MCchrgd = 0;
      for(unsigned int ctr=0; ctr<MCe->size(); ctr++)
	{
	  h_chrg->Fill(MCcharge->at(ctr));
	  if(abs(MCcharge->at(ctr)) != 0) MCchrgd++;
	}

      // RP particle loop
      int RPchrgd = 0;
      for(unsigned int ctr=0; ctr<RPe->size(); ctr++)
	{
	  if(abs(RPcharge->at(ctr)) != 0) RPchrgd++;
	}
      
      if(MCchrgd == RPchrgd)
	{
	  Cmc_Crp++;

	  for(unsigned int it=0; it<RPD0->size(); it++)
	    {
	      if(RPcharge->at(it) != 0)
		{
		  float vtx = sqrt(MCx->at(it)*MCx->at(it) + MCy->at(it)*MCy->at(it) + MCz->at(it)*MCz->at(it));
		  if(abs(RPD0->at(it)/sqrt(RPD0cov->at(it))) < 3) h_PV->Fill(vtx);
		  else h_nPV->Fill(vtx);
		}
	    }
	}
      
      evt++;

      if(evt%50000==0)
	{
	  //cout<<evt<<" events processed"<<endl;
	  cout<<"Event #"<<evt<<":"<<endl;
	  cout<<"There are "<<MCchrgd<<" charged MC particles and "<<RPchrgd<<" charged reco particles"<<endl;
	  cout<<"There are "<<RPe->size()<<" reco particles and "<<RPD0->size()<<" tracks"<<endl;
	  cout<<"====================="<<endl<<endl;
	}
    }

  cout<<Cmc_Crp<<" events out of "<<nEvents<<" events have same no of charged final MC particles as charged reco particles"<<endl;

  file->Close();
  cout<<"Event file closed"<<endl;

  h_chrg->Write();
  h_PV->Write();
  h_nPV->Write();
  histFile->Close();
  cout<<"Histograms written to file and file closed"<<endl;
  
  //delete jetConst;
  //jetConst = NULL;

  cout<<endl;
  return -1;
}
