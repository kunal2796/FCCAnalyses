// Studying the Spring2021 event files (Zbb): ghost matching - reco level (exclusive with exactly 2 jets [for now])
// Note: Close the event file before writing histograms to prevent seg faults
// No cuts
// status 71-79 for parton selection

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

  TFile *file = TFile::Open("p8_ee_Zbb_ecm91_gm7x_auto.root");
  TTreeReader tree("events", file);
  int nEvents = tree.GetEntries();
  cout<<"Number of Events: "<<nEvents<<endl;

  TString histfname;
  histfname = "histZbb_gm7x_auto.root";
  TFile *histFile = new TFile(histfname,"RECREATE");
    
  // hists for the jet loop
  TH1F* h_jetFlavour = new TH1F("h_jetFlavour","Jet Flavour",11,-5,6);
  TH1F* h_jetFlavour_diff = new TH1F("h_jetFlavour_diff","Jet Flavour (MC - GM)",21,-10,11);
  TH1F* h_jetTheta = new TH1F("h_jetTheta","Jet Axis Polar Angle",100,0,3.15);
  TH1F* h_jetPhi = new TH1F("h_jetPhi","Jet Axis Azimuthal Angle",100,-3.15,3.15);
  TH1F* h_pjet = new TH1F("h_pjet","|p| - jets [GeV]",100,0,50);
  TH1F* h_pTjet = new TH1F("h_pTjet","p_T - jets [GeV]",100,0,50);

  // MC particles                                                       
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> MCpx(tree, "MC_px");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> MCpy(tree, "MC_py");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> MCpz(tree, "MC_pz");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> MCe(tree,  "MC_e");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> MCpdg(tree,"MC_pdg");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> MCstatus(tree,"MC_status");
  
  // reco particles                                                       
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> RPpx(tree, "RP_px");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> RPpy(tree, "RP_py");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> RPpz(tree, "RP_pz");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> RPe(tree,  "RP_e");
  
  // jets and jet constituents - eekt     
  TTreeReaderValue<vector<vector<int>>> jetConst(tree, "jetconstituents_ee_kt");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> jetPx(tree, "jets_ee_kt_px");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> jetPy(tree, "jets_ee_kt_py");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> jetPz(tree, "jets_ee_kt_pz");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> jetE(tree,  "jets_ee_kt_e");
  TTreeReaderValue<vector<int,ROOT::Detail::VecOps::RAdoptAllocator<int>>> jetFlavour(tree,"jets_ee_kt_flavour");

  // event counter
  int evt = 0;

  // flavour mismatch counter
  int flv_mm = 0;
  
  // event loop
  while(tree.Next())
    {
      // jet loop
      int nJet = jetE->size();
      float Px=0, Py=0, Pz=0, E=0;
      TLorentzVector p_Jet[nJet];
      for(unsigned int j=0; j<nJet; j++)
	{
	  Px = jetPx->at(j);
	  Py = jetPy->at(j);
	  Pz = jetPz->at(j);
	  E  = jetE->at(j);
	  
	  p_Jet[j].SetPxPyPzE(Px,Py,Pz,E);

	  // jet flavour
	  h_jetFlavour->Fill(jetFlavour->at(j));

	  // jet theta
	  h_jetTheta->Fill(p_Jet[j].Theta());

	  // jet phi
	  h_jetPhi->Fill(p_Jet[j].Phi());	

	  // jet p
	  h_pjet->Fill(p_Jet[j].P());	

	  // jet pT
	  h_pTjet->Fill(p_Jet[j].Pt());	
	}

      if(abs(jetFlavour->at(0)) != abs(jetFlavour->at(1))) flv_mm++;
      //if(jetFlavour->at(0) == 0 || jetFlavour->at(1) == 0) flv_mm++;
      //if(jetFlavour->at(0) == 0) flv_mm += jetFlavour->at(1)/jetFlavour->at(1);
      //if(jetFlavour->at(1) == 0) flv_mm += jetFlavour->at(0)/jetFlavour->at(0);

      // map MC partons with the closest jet and check the diff b/n MC and assigned flavour
      int MC_flv_diff[2];
      float pxMC=0, pyMC=0, pzMC=0, eMC=0;
      TLorentzVector parton;
      for(unsigned int i=0; i<MCe->size(); i++)
	{
	  //if(MCstatus->at(i) < 70 && MCstatus->at(i) > 80) continue;
	  if(MCstatus->at(i) != 23) continue;
	  if(MCpdg->at(i) > 5) continue;

	  pxMC = MCpx->at(i);
	  pyMC = MCpy->at(i);
	  pzMC = MCpz->at(i);
	  eMC  = MCe->at(i);
	  parton.SetPxPyPzE(pxMC, pyMC, pzMC, eMC);

	  // find closer jet
	  float ang1 = parton.Angle(p_Jet[0].Vect());
	  float ang2 = parton.Angle(p_Jet[1].Vect());

	  if(ang1 < ang2) MC_flv_diff[0] = MCpdg->at(i) - jetFlavour->at(0);
	  if(ang1 > ang2) MC_flv_diff[1] = MCpdg->at(i) - jetFlavour->at(1);
	}

      h_jetFlavour_diff->Fill(MC_flv_diff[0]);
      h_jetFlavour_diff->Fill(MC_flv_diff[1]);
      
      /*======================*/
      
      // jet constituents
      vector<int> jet1Const, jet2Const;
      if(jetConst->size()>=1)      jet1Const = jetConst->at(0);
      else cout<<"**No jet constituents found in event#"<<evt+1<<"**"<<endl;
      if(jetConst->size()>=2)      jet2Const = jetConst->at(1);
      else cout<<"**Second jet constituents not found in event#"<<evt+1<<"**"<<endl;
    
      // JET 1
      //cout<<"total reco particles  : "<<RPe->size()<<endl;
      //cout<<"total jet constituents: "<<jet1Const.size()+jet2Const.size()<<endl;
      //cout<<"Jet1 indices: ";
      for(int ele : jet1Const) 
	{
	  //cout<<ele<<", ";
	}
      //cout<<endl;

      // JET 2
      //cout<<"Jet2 indices: ";
      for(int ele : jet2Const) 
	{
	  //cout<<ele<<", ";
	}
      //cout<<endl;
      
      jet1Const.clear();
      jet2Const.clear();

      evt++;

      if(evt%100000==0)
	{
	  cout<<evt<<" events processed"<<endl;
	  //cout<<"Event #"<<evt<<":"<<endl;
	  //cout<<"This event has "<<RPe->size()<<" particles and "<<nJet<<" jets."<<endl;
	  cout<<"====================="<<endl<<endl;
	}
      //cout<<"====================="<<endl<<endl;
    }

  cout<<flv_mm<<" events have jets with differently assigned flavours"<<endl;
  
  file->Close();
  cout<<"Event file closed"<<endl;

  h_jetFlavour->Write();
  h_jetFlavour_diff->Write();
  h_jetTheta->Write();
  h_jetPhi->Write();
  h_pjet->Write();
  h_pTjet->Write();
  histFile->Close();
  cout<<"Histograms written to file and file closed"<<endl;
  
  //delete jetConst;
  //jetConst = NULL;

  cout<<endl;
  return -1;
}
