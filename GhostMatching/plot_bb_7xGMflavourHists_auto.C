#include <iostream>
#include <TFile.h>
#include <TH1F.h>
#include <THStack.h>
#include <TStyle.h>
#include <TCanvas.h>

using namespace std;

int main()
{
  TFile *hists7x=TFile::Open("histZbb_gm7x_auto.root");
  TFile *hists23=TFile::Open("histZbb_gm_auto.root");
  //TFile *hists1=TFile::Open("histZuds_excl_st23.root");
  TFile *hists7xflip=TFile::Open("histZbb_gm7xflip_auto.root");

  //
  
  TH1F* jetFlavour7x = (TH1F*)hists7x->Get("h_jetFlavour");
  TH1F* jetFlavour23 = (TH1F*)hists23->Get("h_jetFlavour");
  //TH1F* jetFlavour23 = (TH1F*)hists1->Get("h_jetFlavour_qqbar");

  //

  TCanvas *c1 = new TCanvas("c1","Jet Flavour : Ghost-Matching",700,700);
  c1->SetLogy();
  c1->SetGrid();
  jetFlavour7x->GetXaxis()->SetTitle("PDG ID");
  jetFlavour7x->Draw();
  //
  c1->Print("jetFlavour_gm7x_bb_auto.pdf");

  cout<<"point 1"<<endl;
  TCanvas *c3 = new TCanvas("c3","Jet Flavour : Ghost-Matching",700,700);
  cout<<"point 2"<<endl;
  THStack *hs3 = new THStack("hs3","Flavour Assignment: GM");
  cout<<"point 3"<<endl;
  gStyle->SetPalette(kOcean);
  cout<<"point 4"<<endl;
  c3->SetLogy();
  c3->SetGrid();
  cout<<"point 5"<<endl;
  hs3->Add(jetFlavour7x);
  hs3->Add(jetFlavour23);
  cout<<"point 6"<<endl;
  //hs3->GetXaxis()->SetTitle("PDG ID");
  cout<<"point 7"<<endl;
  hs3->Draw("nostack");
  cout<<"point 8"<<endl;
  gPad->BuildLegend(0.45,0.95,0.65,0.85,"");
  cout<<"point 9"<<endl;
  //
  c3->Print("jetFlavour_gm_bb_auto.pdf");
  /*
  TCanvas *c2 = new TCanvas("c2","Jet Flavour Assignment",700,700);
  c2->Divide(2,1);
  //
  c2->cd(1);
  c2->cd(1)->SetLogy();
  c2->cd(1)->SetGrid();
  jetFlavour7x->SetTitle("Ghost Matching");
  jetFlavour7x->GetXaxis()->SetTitle("PDG ID");
  jetFlavour7x->SetMinimum(300);
  jetFlavour7x->SetMaximum(700000);
  jetFlavour7x->Draw();
  //
  c2->cd(2);
  c2->cd(2)->SetLogy();
  c2->cd(2)->SetGrid();
  jetFlavour23->SetTitle("Cone Algo with Status 23");
  jetFlavour23->GetXaxis()->SetTitle("PDG ID");
  jetFlavour23->SetMinimum(300);
  jetFlavour23->SetMaximum(700000);
  jetFlavour23->Draw();
  //
  c2->Print("jetFlavour_comparison_7xauto.pdf");  
  */
  hists7x->Close();
  hists23->Close();
  //hists1->Close();
  
  return -1;
}
