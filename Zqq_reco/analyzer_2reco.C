// Studying the Spring2021 event files (Zuds files in particular): event with no jet constituents (exclusive with exactly 2 jets)
// Note: Close the event file before writing histograms to prevent seg faults
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
  /*
  TString histfname;
  histfname = "histZuds_excl.root";
  TFile *histFile = new TFile(histfname,"RECREATE");
  */
  // reco particles                                                       
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> RPpx(tree, "RP_px");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> RPpy(tree, "RP_py");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> RPpz(tree, "RP_pz");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> RPe(tree, "RP_e");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> RPp(tree, "RP_p");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> RPtheta(tree, "RP_theta");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> RPmass(tree, "RP_mass");
  
  // jets and jet constituents - eekt     
  TTreeReaderValue<vector<vector<int>>> jetConst(tree, "jetconstituents_ee_kt");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> jetPx(tree, "jets_ee_kt_px");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> jetPy(tree, "jets_ee_kt_py");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> jetPz(tree, "jets_ee_kt_pz");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> jetE(tree, "jets_ee_kt_e");
  TTreeReaderValue<vector<int,ROOT::Detail::VecOps::RAdoptAllocator<int>>> jetFlavour(tree, "jets_ee_kt_flavour"); //not in the ntuple yet

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
      evt++;
      
      if(RPe->size() != 2) continue;

      cout<<"Event #"<<evt<<" has only 2 particles"<<endl;
      cout<<"---------------------------------"<<endl;
      
      // particle loop
      float px=0., py=0., pz=0., p=0., e=0.;
      TLorentzVector p4[2];
      
      for(unsigned int ctr=0; ctr<RPe->size(); ctr++)
	{	  
	  px = RPpx->at(ctr);
	  py = RPpy->at(ctr);
	  pz = RPpz->at(ctr);
	  e = RPe->at(ctr);

	  p4[ctr].SetPxPyPzE(px, py, pz, e);
	  
	  cout<<"Particle #"<<ctr+1<<endl;
	  cout<<"px = "<<px<<", py = "<<py<<", pz = "<<pz<<", E = "<<e<<endl;
	  cout<<"Mass = "<<RPmass->at(ctr)<<endl;
	  cout<<"---------------------------------"<<endl;
	}

      float angle = p4[0].Angle(p4[1].Vect()) * (180 / 3.14159);
      cout<<"Angle between the two particles = "<<angle<<" degrees"<<endl;

      /*======================*/
      
      // jet loop
      float jPx=0., jPy=0., jPz=0., jE=0.;
      int nJet = jetE->size();
      
      for(unsigned int j=0; j<nJet; j++)
	{
	  jE = jetE->at(j);
	  jPx = jetPx->at(j);
	  jPy = jetPy->at(j);
	  jPz = jetPz->at(j);

	  cout<<"Does the jet loop start?"<<endl;							   
	}

      cout<<"====================="<<endl<<endl;
    }

  file->Close();
  cout<<"Event file closed"<<endl;

  //delete jetConst;
  //jetConst = NULL;
  
  cout<<endl;
  return -1;
}
