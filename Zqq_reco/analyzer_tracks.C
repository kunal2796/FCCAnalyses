// Studying the Spring2021 event files (Zuds files in particular): tracks (exclusive with exactly 2 jets)
// Note: Close the event file before writing histograms to prevent seg faults
// No cuts?

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

  TFile *file = TFile::Open("p8_ee_Zuds_ecm91_reco.root");
  TTreeReader tree("events", file);
  int nEvents = tree.GetEntries();
  cout<<"Number of Events: "<<nEvents<<endl;

  TString histfname;
  histfname = "histZuds_tracks.root";
  TFile *histFile = new TFile(histfname,"RECREATE");
  
  // hists for the track loop
  TH1F* h_d0 = new TH1F("h_d0","Transverse Impact Parameter",100,0,100);
  TH1F* h_z0 = new TH1F("h_z0","Longitudinal Impact Parameter",100,0,100);
  TH1F* h_d0sig = new TH1F("h_d0sig","Transverse Impact Parameter Significance",100,-5,5);
  TH1F* h_z0sig = new TH1F("h_z0sig","Longitudinal Impact Parameter Significance",100,-5,5);

  // hist by flavour
  TH1F* h_d0_b = new TH1F("h_d0_b","Z #rightarrow b#bar{b}",100,0,100);
  TH1F* h_z0_b = new TH1F("h_z0_b","Z #rightarrow b#bar{b}",100,0,100);
  TH1F* h_d0sig_b = new TH1F("h_d0sig_b","Z #rightarrow b#bar{b}",100,-5,5);
  TH1F* h_z0sig_b = new TH1F("h_z0sig_b","Z #rightarrow b#bar{b}",100,-5,5);
  TH1F* h_d0_c = new TH1F("h_d0_c","Z #rightarrow c#bar{c}",100,0,100);
  TH1F* h_z0_c = new TH1F("h_z0_c","Z #rightarrow c#bar{c}",100,0,100);
  TH1F* h_d0sig_c = new TH1F("h_d0sig_c","Z #rightarrow c#bar{c}",100,-5,5);
  TH1F* h_z0sig_c = new TH1F("h_z0sig_c","Z #rightarrow c#bar{c}",100,-5,5);
  TH1F* h_d0_s = new TH1F("h_d0_s","Z #rightarrow s#bar{s}",100,0,100);
  TH1F* h_z0_s = new TH1F("h_z0_s","Z #rightarrow s#bar{s}",100,0,100);
  TH1F* h_d0sig_s = new TH1F("h_d0sig_s","Z #rightarrow s#bar{s}",100,-5,5);
  TH1F* h_z0sig_s = new TH1F("h_z0sig_s","Z #rightarrow s#bar{s}",100,-5,5);
  TH1F* h_d0_u = new TH1F("h_d0_u","Z #rightarrow u#bar{u}",100,0,100);
  TH1F* h_z0_u = new TH1F("h_z0_u","Z #rightarrow u#bar{u}",100,0,100);
  TH1F* h_d0sig_u = new TH1F("h_d0sig_u","Z #rightarrow u#bar{u}",100,-5,5);
  TH1F* h_z0sig_u = new TH1F("h_z0sig_u","Z #rightarrow u#bar{u}",100,-5,5);
  TH1F* h_d0_d = new TH1F("h_d0_d","Z #rightarrow d#bar{d}",100,0,100);
  TH1F* h_z0_d = new TH1F("h_z0_d","Z #rightarrow d#bar{d}",100,0,100);
  TH1F* h_d0sig_d = new TH1F("h_d0sig_d","Z #rightarrow d#bar{d}",100,-5,5);
  TH1F* h_z0sig_d = new TH1F("h_z0sig_d","Z #rightarrow d#bar{d}",100,-5,5);
  
  // reco particles                                                       
  //TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> RPpx(tree, "RP_px");
  //TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> RPpy(tree, "RP_py");
  //TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> RPpz(tree, "RP_pz");
  //TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> RPe(tree, "RP_e");
  //TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> RPp(tree, "RP_p");
  //TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> RPtheta(tree, "RP_theta");
  //TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> RPmass(tree, "RP_mass");

  // track parameters
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> RPD0(tree, "RP_D0");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> RPZ0(tree, "RP_Z0");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> RPphi(tree, "RP_phi");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> RPomega(tree, "RP_omega");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> RPtanLambda(tree, "RP_tanLambda");

  // significance
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> RPD0sig(tree, "RP_D0_sig");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> RPZ0sig(tree, "RP_Z0_sig");

  // jets and jet constituents - eekt     
  //TTreeReaderValue<vector<vector<int>>> jetConst(tree, "jetconstituents_ee_kt");
  //TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> jetPx(tree, "jets_ee_kt_px");
  //TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> jetPy(tree, "jets_ee_kt_py");
  //TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> jetPz(tree, "jets_ee_kt_pz");
  //TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> jetE(tree, "jets_ee_kt_e");
  TTreeReaderValue<vector<int,ROOT::Detail::VecOps::RAdoptAllocator<int>>> jetFlavour(tree, "jets_ee_kt_flavour");

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
      // track loop (only charged ofcourse)
      float d0=0., z0=0.;
      float d0sig=0., z0sig=0.;
      
      for(unsigned int ctr=0; ctr<RPD0->size(); ctr++)
	{	  
	  d0 = RPD0->at(ctr);
	  z0 = RPZ0->at(ctr);

	  d0sig = RPD0->at(ctr);
	  z0sig = RPZ0->at(ctr);
	  
	  // transverse IP
	  h_d0->Fill(abs(d0));
	  if(jetFlavour->size() != 0)
	    {
	      if(jetFlavour->at(0)==5 && jetFlavour->at(1)==5) h_d0_b->Fill(abs(d0)); // b
	      if(jetFlavour->at(0)==4 && jetFlavour->at(1)==4) h_d0_c->Fill(abs(d0)); // c
	      if(jetFlavour->at(0)==3 && jetFlavour->at(1)==3) h_d0_s->Fill(abs(d0)); // s
	      if(jetFlavour->at(0)==2 && jetFlavour->at(1)==2) h_d0_u->Fill(abs(d0)); // u
	      if(jetFlavour->at(0)==1 && jetFlavour->at(1)==1) h_d0_d->Fill(abs(d0)); // d
	    }
	  
	  // longitudinal IP
	  h_z0->Fill(abs(z0));
	  if(jetFlavour->size() != 0)
	    {
	      if(jetFlavour->at(0)==5 && jetFlavour->at(1)==5) h_z0_b->Fill(abs(z0)); // b
	      if(jetFlavour->at(0)==4 && jetFlavour->at(1)==4) h_z0_c->Fill(abs(z0)); // c
	      if(jetFlavour->at(0)==3 && jetFlavour->at(1)==3) h_z0_s->Fill(abs(z0)); // s
	      if(jetFlavour->at(0)==2 && jetFlavour->at(1)==2) h_z0_u->Fill(abs(z0)); // u
	      if(jetFlavour->at(0)==1 && jetFlavour->at(1)==1) h_z0_d->Fill(abs(z0)); // d
	    }
	  
	  // transverse IP significance
	  h_d0sig->Fill(d0sig);
	  if(jetFlavour->size() != 0)
	    {
	      if(jetFlavour->at(0)==5 && jetFlavour->at(1)==5) h_d0sig_b->Fill(d0sig); // b
	      if(jetFlavour->at(0)==4 && jetFlavour->at(1)==4) h_d0sig_c->Fill(d0sig); // c
	      if(jetFlavour->at(0)==3 && jetFlavour->at(1)==3) h_d0sig_s->Fill(d0sig); // s
	      if(jetFlavour->at(0)==2 && jetFlavour->at(1)==2) h_d0sig_u->Fill(d0sig); // u
	      if(jetFlavour->at(0)==1 && jetFlavour->at(1)==1) h_d0sig_d->Fill(d0sig); // d
	    }
	  
	  // longitudinal IP significance
	  h_z0sig->Fill(z0sig);
	  if(jetFlavour->size() != 0)
	    {
	      if(jetFlavour->at(0)==5 && jetFlavour->at(1)==5) h_z0sig_b->Fill(z0sig); // b
	      if(jetFlavour->at(0)==4 && jetFlavour->at(1)==4) h_z0sig_c->Fill(z0sig); // c
	      if(jetFlavour->at(0)==3 && jetFlavour->at(1)==3) h_z0sig_s->Fill(z0sig); // s
	      if(jetFlavour->at(0)==2 && jetFlavour->at(1)==2) h_z0sig_u->Fill(z0sig); // u
	      if(jetFlavour->at(0)==1 && jetFlavour->at(1)==1) h_z0sig_d->Fill(z0sig); // d
	    }
	}

      evt++;

      if(evt%50000==0)
	{
	  //cout<<evt<<" events processed"<<endl;
	  cout<<"Event #"<<evt<<":"<<endl;
	  cout<<"This event has "<<RPD0->size()<<" tracks"<<endl;
	  cout<<"====================="<<endl<<endl;
	}
    }

  file->Close();
  cout<<"Event file closed"<<endl;

  h_d0->Write();
  h_z0->Write();
  h_d0sig->Write();
  h_z0sig->Write();
  h_d0_b->Write();
  h_z0_b->Write();
  h_d0sig_b->Write();
  h_z0sig_b->Write();
  h_d0_c->Write();
  h_z0_c->Write();
  h_d0sig_c->Write();
  h_z0sig_c->Write();
  h_d0_s->Write();
  h_z0_s->Write();
  h_d0sig_s->Write();
  h_z0sig_s->Write();
  h_d0_u->Write();
  h_z0_u->Write();
  h_d0sig_u->Write();
  h_z0sig_u->Write();
  h_d0_d->Write();
  h_z0_d->Write();
  h_d0sig_d->Write();
  h_z0sig_d->Write();
  histFile->Close();
  cout<<"Histograms written to file and file closed"<<endl;

  //delete jetConst;
  //jetConst = NULL;

  cout<<endl;
  return -1;
}
