#include <TFile.h>
#include <TH1F.h>
#include <THStack.h>
#include <TStyle.h>
#include <TCanvas.h>

int main()
{
  TFile *hists1=TFile::Open("histZuds_excl.root");
  TFile *hists2=TFile::Open("histZuds_excl_st23.root");
  
  TH1F* jetFlavour71 = (TH1F*)hists1->Get("h_jetFlavour");
  TH1F* jetFlavour23 = (TH1F*)hists2->Get("h_jetFlavour");
  TH1F* jetFlavour23_qqbar = (TH1F*)hists2->Get("h_jetFlavour_qqbar");
  
  TCanvas *c1 = new TCanvas("c1","Jet Flavour",700,700);
  c1->Divide(2,1);
  //
  c1->cd(1);
  c1->cd(1)->SetLogy();
  c1->cd(1)->SetGrid();
  jetFlavour71->GetXaxis()->SetTitle("|PDG ID|");
  jetFlavour71->SetMinimum(300);
  jetFlavour71->SetMaximum(700000);
  jetFlavour71->Draw();
  //
  c1->cd(2);
  c1->cd(2)->SetLogy();
  c1->cd(2)->SetGrid();
  jetFlavour23->GetXaxis()->SetTitle("|PDG ID|");
  jetFlavour23->SetMinimum(300);
  jetFlavour23->SetMaximum(700000);
  jetFlavour23->Draw();
  //
  c1->Print("jetFlavour.pdf");  

  TCanvas *c2 = new TCanvas("c2","Jet Flavour - q#bar{q}",710,710);
  c2->SetLogy();
  c2->SetGrid();
  jetFlavour23_qqbar->GetXaxis()->SetTitle("PDG ID");
  jetFlavour23_qqbar->Draw();
  //
  c2->Print("jetFlavour_qqbar.pdf");  
  
  hists1->Close();
  hists2->Close();

  return -1;
}
