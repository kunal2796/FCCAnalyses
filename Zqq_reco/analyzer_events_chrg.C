// Studying the Spring2021 event files (Zuds files in particular): events - reco level (exclusive with exactly 2 jets)
// charged and neutrals separate
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
#include <TH2I.h>
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
  histfname = "histZuds_chrg.root";
  TFile *histFile = new TFile(histfname,"RECREATE");
  
  // hists for charged particles
  TH1F* h_n_chrg = new TH1F("h_n_chrg","Multiplicity",40,0,40);
  TH1F* h_pT_chrg = new TH1F("h_pT_chrg","p_{T} [GeV]",100,0,20);
  TH1F* h_p_chrg = new TH1F("h_p_chrg","|p| [GeV]",100,0,30);
  TH1F* h_e_chrg = new TH1F("h_e_chrg","E [GeV]",100,0,30);
  TH1F* h_theta_chrg = new TH1F("h_theta_chrg","Polar Angle (#theta)",100,0,3.15);
  TH1F* h_phi_chrg = new TH1F("h_phi_chrg","Azimuthal Angle (#phi)",100,-3.15,3.15);
  TH1F* h_invM_chrg = new TH1F("h_invM_chrg","Invariant Mass (event) [GeV]",100,0,100);
  TH1F* h_m_RP_chrg = new TH1F("h_m_RP_chrg","Mass (reco particles) [GeV]",100,0,0.14);
  
  // hists for neutral particles
  TH1F* h_n_neut = new TH1F("h_n_neut","Multiplicity",30,0,30);
  TH1F* h_pT_neut = new TH1F("h_pT_neut","p_{T} [GeV]",100,0,20);
  TH1F* h_p_neut = new TH1F("h_p_neut","|p| [GeV]",100,0,30);
  TH1F* h_e_neut = new TH1F("h_e_neut","E [GeV]",100,0,30);
  TH1F* h_theta_neut = new TH1F("h_theta_neut","Polar Angle (#theta)",100,0,3.15);
  TH1F* h_phi_neut = new TH1F("h_phi_neut","Azimuthal Angle (#phi)",100,-3.15,3.15);

  // 2D hists
  TH2I* h_CvsN = new TH2I("h_CvsN","Charged vs Neutral - Multiplicity",100,0,50, 100,0,50);

  // reco particles                                                       
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> RPpx(tree, "RP_px");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> RPpy(tree, "RP_py");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> RPpz(tree, "RP_pz");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> RPe(tree, "RP_e");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> RPp(tree, "RP_p");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> RPtheta(tree, "RP_theta");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> RPmass(tree, "RP_mass");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> RPcharge(tree, "RP_charge");
 
  // jets and jet constituents - eekt
  //TTreeReaderValue<vector<vector<int>>> jetConst(tree, "jetconstituents_ee_kt");
  //TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> jetPx(tree, "jets_ee_kt_px");
  //TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> jetPy(tree, "jets_ee_kt_py");
  //TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> jetPz(tree, "jets_ee_kt_pz");
  //TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> jetE(tree, "jets_ee_kt_e");
  //TTreeReaderValue<vector<int,ROOT::Detail::VecOps::RAdoptAllocator<int>>> jetFlavour(tree, "jets_ee_kt_flavour"); //not in the ntuple yet

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
      // particle loop
      float px=0., py=0., pz=0., p=0., e=0.;
      TLorentzVector p4, p4_chrg;
      int n_chrg=0, n_neut=0;
      
      for(unsigned int ctr=0; ctr<RPe->size(); ctr++)
	{	  
	  px = RPpx->at(ctr);
	  py = RPpy->at(ctr);
	  pz = RPpz->at(ctr);
	  p = RPp->at(ctr);
	  e = RPe->at(ctr);

	  p4.SetPxPyPzE(px, py, pz, e);

	  // charged
	  if(RPcharge->at(ctr) != 0)
	    {
	      n_chrg++;
	      p4_chrg += p4;
	      
	      // polar angle
	      h_theta_chrg->Fill(RPtheta->at(ctr));
	      
	      // azimuthal angle
	      h_phi_chrg->Fill(p4.Phi());
	      
	      // abs momenta
	      h_p_chrg->Fill(p);
	      
	      // transverse momenta
	      h_pT_chrg->Fill(p4.Pt());
	      
	      // energy
	      h_e_chrg->Fill(e);
	      
	      // mass (from FCCAnalyses function)
	      h_m_RP_chrg->Fill(RPmass->at(ctr));
	    }

	  // neutral
	  else
	    {	      
	      n_neut++;
	      
	      // polar angle
	      h_theta_neut->Fill(RPtheta->at(ctr));
	      
	      // azimuthal angle
	      h_phi_neut->Fill(p4.Phi());
	      
	      // abs momenta
	      h_p_neut->Fill(p);
	      
	      // transverse momenta
	      h_pT_neut->Fill(p4.Pt());
	      
	      // energy
	      h_e_neut->Fill(e);
	    }
	}

      // event multiplicity
      h_n_chrg->Fill(n_chrg);       // charged
      h_n_neut->Fill(n_neut);       // neutral
      h_CvsN->Fill(n_chrg, n_neut); // charged vs neutral
      
      // invariant mass - (charged particle sum)
      h_invM_chrg->Fill(p4_chrg.M());

      evt++;

      if(evt%50000==0)
	{
	  //cout<<evt<<" events processed"<<endl;
	  cout<<"Event #"<<evt<<":"<<endl;
	  cout<<"This event has "<<n_chrg<<" charged particles and "<<n_neut<<" neutral particles"<<endl;
	  cout<<"====================="<<endl<<endl;
	}
    }

  file->Close();
  cout<<"Event file closed"<<endl;

  h_n_chrg->Write();
  h_pT_chrg->Write();
  h_p_chrg->Write();
  h_e_chrg->Write();
  h_theta_chrg->Write();
  h_phi_chrg->Write();
  h_invM_chrg->Write();
  h_m_RP_chrg->Write();
  h_n_neut->Write();
  h_pT_neut->Write();
  h_p_neut->Write();
  h_e_neut->Write();
  h_theta_neut->Write();
  h_phi_neut->Write();
  h_CvsN->Write();
  histFile->Close();
  cout<<"Histograms written to file and file closed"<<endl;

  cout<<endl;
  return -1;
}
