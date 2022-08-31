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

  //TFile *file = TFile::Open("p8_ee_Zbb_ecm91_SV_100K.root");
  TFile *file = TFile::Open("p8_ee_Zcc_ecm91_SV_100K.root");
  //TFile *file = TFile::Open("p8_ee_Zuds_ecm91_SV_10K.root");
  //TFile *file = TFile::Open("p8_ee_Zuds_ecm91_SV_100K.root");
  TTreeReader tree("events", file);
  int nEvents = tree.GetEntries();
  cout<<"Number of Events: "<<nEvents<<endl;
  
  TString histfname;
  //histfname = "histZbb_SV.root";
  histfname = "histZcc_SV.root";
  //histfname = "histZuds_SV.root";
  TFile *histFile = new TFile(histfname,"RECREATE");

  //TH1F* Zbb = new TH1F("Zbb","Z #rightarrow b#bar{b}",9,0.,9.);
  TH1F* Zcc = new TH1F("Zcc","Z #rightarrow c#bar{c}",9,0.,9.);
  //TH1F* Zuds = new TH1F("Zuds","Z #rightarrow uds",9,0.,9.);
  
  // hists for jet SVs
  TH1F* h_chi2_SVjet = new TH1F("h_chi2_SVjet","SV #chi^{2} (Per Jet)",100,0.,9.);

  // hists for event SVs
  TH1F* h_chi2_SVevt = new TH1F("h_chi2_SVevt","SV #chi^{2} (Per Event)",100,0.,9.);

  // PV
  
  // SV (per jet)
  //TTreeReaderValue<vector<TVector3,ROOT::Detail::VecOps::RAdoptAllocator<TVector3>>> SVpos_jet(tree, "SV_pos_jet");
  TTreeReaderValue<int> SVn_jet(tree, "SV_jet_n");
  TTreeReaderValue<vector<double,ROOT::Detail::VecOps::RAdoptAllocator<double>>> SVchi2_jet(tree, "SV_jet_chi2");

  // SV (per event)
  //TTreeReaderValue<vector<TVector3,ROOT::Detail::VecOps::RAdoptAllocator<TVector3>>> SVpos_evt(tree, "SV_pos_evt");
  TTreeReaderValue<vector<double,ROOT::Detail::VecOps::RAdoptAllocator<double>>> SVchi2_evt(tree, "SV_evt2_chi2");

  // event counter
  int evt = 0;

  // event loop
  while(tree.Next())
    {
      //Zbb->Fill(*SVn_jet);
      Zcc->Fill(*SVn_jet);
      //Zuds->Fill(*SVn_jet);
      
      // SV (per jet) loop
      for(unsigned int ctr=0; ctr<SVchi2_jet->size(); ctr++)
	{
	  h_chi2_SVjet->Fill(SVchi2_jet->at(ctr));
	}

      // SV (per evt) loop
      for(unsigned int ctr=0; ctr<SVchi2_evt->size(); ctr++)
	{
	  h_chi2_SVevt->Fill(SVchi2_evt->at(ctr));
	}

      evt++;

      if(evt%10000==0)
	{
	  //cout<<evt<<" events processed"<<endl;
	  cout<<"Event #"<<evt<<":"<<endl;
	  cout<<"====================="<<endl<<endl;
	}
    }

  file->Close();
  cout<<"Event file closed"<<endl;

  //Zbb->Write();
  Zcc->Write();
  //Zuds->Write();
  h_chi2_SVjet->Write();
  h_chi2_SVevt->Write();
  histFile->Close();
  cout<<"Histograms written to file and file closed"<<endl;

  return -1;
}
