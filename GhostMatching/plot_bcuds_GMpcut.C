#include <TFile.h>
#include <TH1F.h>
#include <THStack.h>
#include <TStyle.h>
#include <TCanvas.h>

int main()
{
  TFile *hists_uds_pcut=TFile::Open("histZuds_gm7x_pcut.root");
  TFile *hists_uds     =TFile::Open("histZuds_gm7x_auto.root");
  TFile *hists_bb_pcut=TFile::Open("histZbb_gm7x_pcut.root");
  TFile *hists_bb     =TFile::Open("histZbb_gm7x_auto.root");
  TFile *hists_cc_pcut=TFile::Open("histZcc_gm7x_pcut.root");
  TFile *hists_cc     =TFile::Open("histZcc_gm7x_auto.root");
  
  TH1F* jetFlavour_uds     = (TH1F*)hists_uds->Get("h_jetFlavour");
  TH1F* jetFlavour_uds_pcut= (TH1F*)hists_uds_pcut->Get("h_jetFlavour");
  TH1F* jetFlavour_bb      = (TH1F*)hists_bb->Get("h_jetFlavour");
  TH1F* jetFlavour_bb_pcut = (TH1F*)hists_bb_pcut->Get("h_jetFlavour");
  TH1F* jetFlavour_cc      = (TH1F*)hists_cc->Get("h_jetFlavour");
  TH1F* jetFlavour_cc_pcut = (TH1F*)hists_cc_pcut->Get("h_jetFlavour");

  //

  TCanvas *c1 = new TCanvas("c1","Jet Flavour Assignment: Zuds",710,710);
  c1->Divide(2,1);
  //
  c1->cd(1);
  c1->cd(1)->SetLogy();
  c1->cd(1)->SetGrid();
  jetFlavour_uds->SetTitle("Z #rightarrow uds: w/o Momentum Cut");
  jetFlavour_uds->GetXaxis()->SetTitle("PDG ID");
  jetFlavour_uds->SetMinimum(300);
  jetFlavour_uds->SetMaximum(700000);
  jetFlavour_uds->Draw();
  //
  c1->cd(2);
  c1->cd(2)->SetLogy();
  c1->cd(2)->SetGrid();
  jetFlavour_uds_pcut->SetTitle("Z #rightarrow uds: w/ Momentum Cut");
  jetFlavour_uds_pcut->GetXaxis()->SetTitle("PDG ID");
  jetFlavour_uds_pcut->SetMinimum(300);
  jetFlavour_uds_pcut->SetMaximum(700000);
  jetFlavour_uds_pcut->Draw();
  //
  c1->Print("GM_uds_p_cut.pdf");  

  TCanvas *c2 = new TCanvas("c2","Jet Flavour Assignment: Zbb",720,720);
  c2->Divide(2,1);
  //
  c2->cd(1);
  c2->cd(1)->SetLogy();
  c2->cd(1)->SetGrid();
  jetFlavour_bb->SetTitle("Z #rightarrow b#bar{b}: w/o Momentum Cut");
  jetFlavour_bb->GetXaxis()->SetTitle("PDG ID");
  jetFlavour_bb->SetMinimum(300);
  jetFlavour_bb->SetMaximum(3000000);
  jetFlavour_bb->Draw();
  //
  c2->cd(2);
  c2->cd(2)->SetLogy();
  c2->cd(2)->SetGrid();
  jetFlavour_bb_pcut->SetTitle("Z #rightarrow b#bar{b}: w/ Momentum Cut");
  jetFlavour_bb_pcut->GetXaxis()->SetTitle("PDG ID");
  jetFlavour_bb_pcut->SetMinimum(300);
  jetFlavour_bb_pcut->SetMaximum(3000000);
  jetFlavour_bb_pcut->Draw();
  //
  c2->Print("GM_bb_p_cut.pdf");  

  TCanvas *c3 = new TCanvas("c3","Jet Flavour Assignment: Zcc",730,730);
  c3->Divide(2,1);
  //
  c3->cd(1);
  c3->cd(1)->SetLogy();
  c3->cd(1)->SetGrid();
  jetFlavour_cc->SetTitle("Z #rightarrow c#bar{c}: w/o Momentum Cut");
  jetFlavour_cc->GetXaxis()->SetTitle("PDG ID");
  jetFlavour_cc->SetMinimum(300);
  jetFlavour_cc->SetMaximum(3000000);
  jetFlavour_cc->Draw();
  //
  c3->cd(2);
  c3->cd(2)->SetLogy();
  c3->cd(2)->SetGrid();
  jetFlavour_cc_pcut->SetTitle("Z #rightarrow c#bar{c}: w/ Momentum Cut");
  jetFlavour_cc_pcut->GetXaxis()->SetTitle("PDG ID");
  jetFlavour_cc_pcut->SetMinimum(300);
  jetFlavour_cc_pcut->SetMaximum(3000000);
  jetFlavour_cc_pcut->Draw();
  //
  c3->Print("GM_cc_p_cut.pdf");  

  hists_uds->Close();
  hists_uds_pcut->Close();
  hists_bb->Close();
  hists_bb_pcut->Close();
  hists_cc->Close();
  hists_cc_pcut->Close();

  return -1;
}
