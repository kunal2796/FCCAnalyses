// Studying the Spring2021 event files (Zuds files in particular): events, jets - reco level (ee-genkt)
// No kinematic cuts
// studying the flavour assignment on two leading jets in ee-genkt - currently only counting

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
  histfname = "histZuds_incl.root";
  TFile *histFile = new TFile(histfname,"RECREATE");
  
  // hists for the jet loop
  TH1F* h_njet = new TH1F("h_njet","Jet Multiplicity",20,0,20);
  TH1F* h_pjet = new TH1F("h_pjet","|p| - jets [GeV]",100,0,50);
  TH1F* h_pTjet = new TH1F("h_pTjet","p_{T} - jets [GeV]",100,0,50);
  TH1F* h_Ejet = new TH1F("h_Ejet","Energy - jet [GeV]",100,0,50);
  TH1F* h_jetTheta = new TH1F("h_jetTheta","Jet Polar Angle (#theta)",100,0,3.15);
  TH1F* h_jetPhi = new TH1F("h_jetPhi","Jet Azimuthal Angle (#phi)",100,-3.15,3.15);
  TH1F* h_invMjets = new TH1F("h_invMjets","Invariant Mass - sum of jets [GeV]",100,75,100);
  TH1F* h_invMjets_b = new TH1F("h_invMjets_b","b-jets",100,0,100);
  TH1F* h_invMjets_c = new TH1F("h_invMjets_c","c-jets",100,0,100);
  TH1F* h_invMjets_s = new TH1F("h_invMjets_s","s-jets",100,0,100);
  TH1F* h_invMjets_u = new TH1F("h_invMjets_u","u-jets",100,0,100);
  TH1F* h_invMjets_d = new TH1F("h_invMjets_d","d-jets",100,0,100);
  TH1F* h_invMjets_bHigh = new TH1F("h_invMjets_bHigh","b-jets",100,0,100);
  TH1F* h_invMjets_cHigh = new TH1F("h_invMjets_cHigh","c-jets",100,0,100);
  TH1F* h_invMjets_sHigh = new TH1F("h_invMjets_sHigh","s-jets",100,0,100);
  TH1F* h_invMjets_uHigh = new TH1F("h_invMjets_uHigh","u-jets",100,0,100);
  TH1F* h_invMjets_dHigh = new TH1F("h_invMjets_dHigh","d-jets",100,0,100);
  TH1F* h_dijet_flavour = new TH1F("h_dijet_flavour","Leading and Sub-Leading Jet Flavours",6,0,6);
  TH1F* h_invMdijets = new TH1F("h_invMdijets","Invariant Mass - sum of 2 leading jets",100,0,100);
  TH1F* h_invMdijets_b = new TH1F("h_invMdijets_b","b-jets",100,0,100);
  TH1F* h_invMdijets_c = new TH1F("h_invMdijets_c","c-jets",100,0,100);
  TH1F* h_invMdijets_s = new TH1F("h_invMdijets_s","s-jets",100,0,100);
  TH1F* h_invMdijets_u = new TH1F("h_invMdijets_u","u-jets",100,0,100);
  TH1F* h_invMdijets_d = new TH1F("h_invMdijets_d","d-jets",100,0,100);
  
  // hists for the jet constituents
  TH1F* h_angJP = new TH1F("h_angJP","Angle b/n Jet Constituents and Jet Axis",100,0,1.5);
  TH1F* h_thetaJP = new TH1F("h_thetaJP","#Delta#theta Jet Constituents & Jet Axis",100,-1.5,1.5);
  TH1F* h_phiJP = new TH1F("h_phiJP","#Delta#phi Jet Constituents & Jet Axis",100,-3.15,3.15);

  // reco particles                                                       
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> RPpx(tree, "RP_px");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> RPpy(tree, "RP_py");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> RPpz(tree, "RP_pz");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> RPe(tree, "RP_e");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> RPp(tree, "RP_p");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> RPtheta(tree, "RP_theta");
  
  // jets and jet constituents - ee-genkt     
  TTreeReaderValue<vector<vector<int>>> jetConst(tree, "jetconstituents_ee_genkt");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> jetPx(tree, "jets_ee_genkt_px");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> jetPy(tree, "jets_ee_genkt_py");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> jetPz(tree, "jets_ee_genkt_pz");
  TTreeReaderValue<vector<float,ROOT::Detail::VecOps::RAdoptAllocator<float>>> jetE(tree, "jets_ee_genkt_e");
  TTreeReaderValue<vector<int,ROOT::Detail::VecOps::RAdoptAllocator<int>>> jetFlavour(tree, "jets_ee_genkt_flavour");

  // event counter
  int evt = 0;
  
  // event loop
  while(tree.Next())
    {
      // jet loop
      float jPx=0., jPy=0., jPz=0., jE=0.;
      int nJet = jetE->size();
      TLorentzVector p_Jet[nJet], p_Jets, p_Jets_b, p_Jets_c, p_Jets_s, p_Jets_u, p_Jets_d;

      TLorentzVector p_diJet[2];
      int flavour_diJet[2];
      
      for(unsigned int j=0; j<nJet; j++)
	{
	  jE = jetE->at(j);
	  jPx = jetPx->at(j);
	  jPy = jetPy->at(j);
	  jPz = jetPz->at(j);

	  p_Jet[j].SetPxPyPzE(jPx, jPy, jPz, jE);

	  // jet theta
	  h_jetTheta->Fill(p_Jet[j].Theta());
	  
	  // jet phi
	  h_jetPhi->Fill(p_Jet[j].Phi());
	  
	  // jet |p|
	  h_pjet->Fill(p_Jet[j].P());

	  // jet pT  
	  h_pTjet->Fill(p_Jet[j].Pt());

	  // jet E
	  h_Ejet->Fill(jE);

	  // sum of jets' momenta (flavoured)
	  p_Jets += p_Jet[j];                            // all jets
	  if(jetFlavour->at(j)==5) p_Jets_b += p_Jet[j]; // b-jets
	  if(jetFlavour->at(j)==4) p_Jets_c += p_Jet[j]; // c-jets
	  if(jetFlavour->at(j)==3) p_Jets_s += p_Jet[j]; // s-jets
	  if(jetFlavour->at(j)==2) p_Jets_u += p_Jet[j]; // u-jets
	  if(jetFlavour->at(j)==1) p_Jets_d += p_Jet[j]; // d-jets

	  // leading and sub-leading jet
	  if(nJet > 1)
	    {
	      if(p_Jet[j].P() > p_diJet[0].P())
		{
		  p_diJet[1] = p_diJet[0];
		  flavour_diJet[1] = flavour_diJet[0];
		  
		  p_diJet[0] = p_Jet[j];
		  flavour_diJet[0] = jetFlavour->at(j);
		}
	      else if(p_Jet[j].P() > p_diJet[1].P())
		{
		  p_diJet[1] = p_Jet[j];
		  flavour_diJet[1] = jetFlavour->at(j);
		}
	    }
	}

      // leading & sub-leading jet flavour
      //if(nJet == 1) h_dijet_flavour->Fill(jetFlavour->at(0));
      if(nJet > 1)
	{
	  h_dijet_flavour->Fill(flavour_diJet[0]);
	  h_dijet_flavour->Fill(flavour_diJet[1]);
	}
            
      // jet multiplicity
      h_njet->Fill(nJet);
      
      // invariant mass - (jet sum)
      h_invMjets->Fill(p_Jets.M());                         // entire dataset
      if(p_Jets_b.M()!=0) h_invMjets_b->Fill(p_Jets_b.M()); // b-jets
      if(p_Jets_c.M()!=0) h_invMjets_c->Fill(p_Jets_c.M()); // c-jets
      if(p_Jets_s.M()!=0) h_invMjets_s->Fill(p_Jets_s.M()); // s-jets
      if(p_Jets_u.M()!=0) h_invMjets_u->Fill(p_Jets_u.M()); // u-jets
      if(p_Jets_d.M()!=0) h_invMjets_d->Fill(p_Jets_d.M()); // d-jets

      // invariant mass - (highest jet sum)
      if(p_Jets_b.M()>p_Jets_c.M() && p_Jets_b.M()>p_Jets_s.M() && p_Jets_b.M()>p_Jets_u.M() && p_Jets_b.M()>p_Jets_d.M()) h_invMjets_bHigh->Fill(p_Jets_b.M()); // b-jets
      if(p_Jets_c.M()>p_Jets_b.M() && p_Jets_c.M()>p_Jets_s.M() && p_Jets_c.M()>p_Jets_u.M() && p_Jets_c.M()>p_Jets_d.M()) h_invMjets_cHigh->Fill(p_Jets_c.M()); // c-jets
      if(p_Jets_s.M()>p_Jets_b.M() && p_Jets_s.M()>p_Jets_c.M() && p_Jets_s.M()>p_Jets_u.M() && p_Jets_s.M()>p_Jets_d.M()) h_invMjets_sHigh->Fill(p_Jets_s.M()); // s-jets
      if(p_Jets_u.M()>p_Jets_b.M() && p_Jets_u.M()>p_Jets_c.M() && p_Jets_u.M()>p_Jets_s.M() && p_Jets_u.M()>p_Jets_d.M()) h_invMjets_uHigh->Fill(p_Jets_u.M()); // u-jets
      if(p_Jets_d.M()>p_Jets_b.M() && p_Jets_d.M()>p_Jets_c.M() && p_Jets_d.M()>p_Jets_s.M() && p_Jets_d.M()>p_Jets_u.M()) h_invMjets_dHigh->Fill(p_Jets_d.M()); // d-jets

      // invariant mass - (leading & sub-leading jets)
      if(nJet > 1)
	{
	  TLorentzVector p_diJets = p_diJet[0] + p_diJet[1];
	  h_invMdijets->Fill(p_diJets.M());
	  if(flavour_diJet[0] == 5 && flavour_diJet[1] == 5) h_invMdijets_b->Fill(p_diJets.M()); // b-jets
	  if(flavour_diJet[0] == 4 && flavour_diJet[1] == 4) h_invMdijets_c->Fill(p_diJets.M()); // c-jets
	  if(flavour_diJet[0] == 3 && flavour_diJet[1] == 3) h_invMdijets_s->Fill(p_diJets.M()); // s-jets
	  if(flavour_diJet[0] == 2 && flavour_diJet[1] == 2) h_invMdijets_u->Fill(p_diJets.M()); // u-jets
	  if(flavour_diJet[0] == 1 && flavour_diJet[1] == 1) h_invMdijets_d->Fill(p_diJets.M()); // d-jets
	}

      /*======================*/
      
      // jet constituents
      vector<int> jet1Const, jet2Const;
      if(jetConst->size()>=1)      jet1Const = jetConst->at(0);
      else cout<<"**No jet constituents found in event#"<<evt+1<<"**"<<endl;
      if(jetConst->size()>=2)      jet2Const = jetConst->at(1);
      else cout<<"**Second jet constituents not found in event#"<<evt+1<<"**"<<endl;
      
      // JET 1
      float px_j1=0, py_j1=0, pz_j1=0, e_j1=0;
      TLorentzVector p4_j1;

      for(int ele : jet1Const) 
	{
	  px_j1 = RPpx->at(ele);
	  py_j1 = RPpy->at(ele);
	  pz_j1 = RPpz->at(ele);
	  e_j1 = RPe->at(ele);
	  
	  p4_j1.SetPxPyPzE(px_j1, py_j1, pz_j1, e_j1);

	  h_angJP->Fill(p4_j1.Angle(p_Jet[0].Vect()));     // angle b/n jet const and jet
	  h_thetaJP->Fill(p4_j1.Theta()-p_Jet[0].Theta()); // delta theta
	  h_phiJP->Fill(p4_j1.DeltaPhi(p_Jet[0]));         // delta phi
	}

      // JET 2
      float px_j2=0, py_j2=0, pz_j2=0, e_j2=0;
      TLorentzVector p4_j2;

      for(int ele : jet2Const) 
	{
	  px_j2 = RPpx->at(ele);
	  py_j2 = RPpy->at(ele);
	  pz_j2 = RPpz->at(ele);
	  e_j2 = RPe->at(ele);
	  
	  p4_j2.SetPxPyPzE(px_j2, py_j2, pz_j2, e_j2);

	  h_angJP->Fill(p4_j2.Angle(p_Jet[1].Vect()));     // angle b/n jet const and jet
	  h_thetaJP->Fill(p4_j2.Theta()-p_Jet[1].Theta()); // delta theta
	  h_phiJP->Fill(p4_j2.DeltaPhi(p_Jet[1]));         // delta phi
	}
           
      jet1Const.clear();
      jet2Const.clear();

      evt++;

      if(evt%50000==0)
	{
	  //cout<<evt<<" events processed"<<endl;
	  cout<<"Event #"<<evt<<":"<<endl;
	  cout<<"This event has "<<RPe->size()<<" particles and "<<nJet<<" jets."<<endl;
	  cout<<"====================="<<endl<<endl;
	}
    }

  file->Close();
  cout<<"Event file closed"<<endl;

  h_njet->Write();
  h_pjet->Write();
  h_pTjet->Write();
  h_Ejet->Write();
  h_jetTheta->Write();
  h_jetPhi->Write();
  h_invMjets->Write();
  h_invMjets_b->Write();
  h_invMjets_c->Write();
  h_invMjets_s->Write();
  h_invMjets_u->Write();
  h_invMjets_d->Write();
  h_invMjets_bHigh->Write();
  h_invMjets_cHigh->Write();
  h_invMjets_sHigh->Write();
  h_invMjets_uHigh->Write();
  h_invMjets_dHigh->Write();
  h_dijet_flavour->Write();
  h_invMdijets->Write();
  h_invMdijets_b->Write();
  h_invMdijets_c->Write();
  h_invMdijets_s->Write();
  h_invMdijets_u->Write();
  h_invMdijets_d->Write();
  h_angJP->Write();
  h_thetaJP->Write();
  h_phiJP->Write();
  histFile->Close();
  cout<<"Histograms written to file and file closed"<<endl;

  //delete jetConst;
  //jetConst = NULL;
  
  cout<<endl;
  return -1;
}
