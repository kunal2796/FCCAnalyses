#include <TFile.h>
#include <TH1F.h>
#include <THStack.h>
#include <TStyle.h>
#include <TCanvas.h>

int main()
{
  TFile *hists23=TFile::Open("histZuds_gm.root");
  TFile *hists7x=TFile::Open("histZuds_gm7x.root");
  
  TH1F* jetFlavour23     = (TH1F*)hists23->Get("h_jetFlavour");
  TH1F* jetFlavour_org23 = (TH1F*)hists23->Get("h_jetFlavour_org");
  TH1F* jetTheta23       = (TH1F*)hists23->Get("h_jetTheta");
  TH1F* jetPhi23         = (TH1F*)hists23->Get("h_jetPhi");
  //
  TH1F* jetFlavour7x     = (TH1F*)hists7x->Get("h_jetFlavour");
  TH1F* jetTheta7x       = (TH1F*)hists7x->Get("h_jetTheta");
  TH1F* jetPhi7x         = (TH1F*)hists7x->Get("h_jetPhi");

  //

  TCanvas *c1 = new TCanvas("c1","Jet Flavour Assignment",710,710);
  c1->Divide(2,1);
  //
  c1->cd(1);
  c1->cd(1)->SetLogy();
  c1->cd(1)->SetGrid();
  jetFlavour23->SetTitle("GM with Status 23");
  jetFlavour23->GetXaxis()->SetTitle("PDG ID");
  jetFlavour23->SetMinimum(300);
  jetFlavour23->SetMaximum(700000);
  jetFlavour23->Draw();
  //
  c1->cd(2);
  c1->cd(2)->SetLogy();
  c1->cd(2)->SetGrid();
  jetFlavour7x->SetTitle("GM with Status 71-79");
  jetFlavour7x->GetXaxis()->SetTitle("PDG ID");
  jetFlavour7x->SetMinimum(300);
  jetFlavour7x->SetMaximum(700000);
  jetFlavour7x->Draw();
  //
  c1->Print("GM_comparison.pdf");  

  TCanvas *c2 = new TCanvas("c2","Jet-Flavour Assignment",710,710);
  c2->Divide(1,3);
  //
  c2->cd(1);
  c2->cd(1)->SetLogy();
  c2->cd(1)->SetGrid();
  jetFlavour23->SetTitle("GM with Status 23");
  jetFlavour23->GetXaxis()->SetTitle("PDG ID");
  jetFlavour23->SetMinimum(900);
  jetFlavour23->SetMaximum(700000);
  jetFlavour23->Draw();
  //
  c2->cd(2);
  c2->cd(2)->SetLogy();
  c2->cd(2)->SetGrid();
  jetFlavour7x->SetTitle("GM with Status 71-79");
  jetFlavour7x->GetXaxis()->SetTitle("PDG ID");
  jetFlavour7x->SetMinimum(900);
  jetFlavour7x->SetMaximum(700000);
  jetFlavour7x->Draw();
  //
  c2->cd(3);
  c2->cd(3)->SetLogy();
  c2->cd(3)->SetGrid();
  jetFlavour_org23->SetTitle("Original Jet Flavour");
  jetFlavour_org23->GetXaxis()->SetTitle("PDG ID");
  jetFlavour_org23->SetMinimum(900);
  jetFlavour_org23->SetMaximum(700000);
  jetFlavour_org23->Draw();
  //
  c2->Print("GM_comparison2.pdf");  
  
  hists23->Close();
  hists7x->Close();

  return -1;
}
