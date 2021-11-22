#include <TFile.h>
#include <TH1F.h>
#include <THStack.h>
#include <TStyle.h>
#include <TCanvas.h>

int main()
{
  TFile *hists=TFile::Open("histZuds_chrg.root");
  
  TH1F* pT_chrg = (TH1F*)hists->Get("h_pT_chrg");
  TH1F* p_chrg = (TH1F*)hists->Get("h_p_chrg");
  TH1F* e_chrg = (TH1F*)hists->Get("h_e_chrg");
  TH1F* theta_chrg = (TH1F*)hists->Get("h_theta_chrg");
  TH1F* phi_chrg = (TH1F*)hists->Get("h_phi_chrg");
  TH1F* invM_chrg = (TH1F*)hists->Get("h_invM_chrg");
  TH1F* m_RP_chrg = (TH1F*)hists->Get("h_m_RP_chrg");
  TH1F* n_chrg = (TH1F*)hists->Get("h_n_chrg");

  TH1F* pT_neut = (TH1F*)hists->Get("h_pT_neut");
  TH1F* p_neut = (TH1F*)hists->Get("h_p_neut");
  TH1F* e_neut = (TH1F*)hists->Get("h_e_neut");
  TH1F* theta_neut = (TH1F*)hists->Get("h_theta_neut");
  TH1F* phi_neut = (TH1F*)hists->Get("h_phi_neut");
  TH1F* n_neut = (TH1F*)hists->Get("h_n_neut");

  
  TCanvas *c1 = new TCanvas("c1","Charged Particles",700,700);
  c1->Divide(2,2);
  //
  c1->cd(1);
  n_chrg->GetXaxis()->SetTitle("#");
  n_chrg->Draw();
  //
  c1->cd(2);
  p_chrg->GetXaxis()->SetTitle("|p| [GeV]");
  c1->cd(2)->SetLogy();
  p_chrg->Draw();
  //
  c1->cd(3);
  e_chrg->GetXaxis()->SetTitle("Energy [GeV]");
  c1->cd(3)->SetLogy();
  e_chrg->Draw();
  //
  c1->cd(4);
  pT_chrg->GetXaxis()->SetTitle("p_{T} [GeV]");
  c1->cd(4)->SetLogy();
  pT_chrg->Draw();
  //
  c1->Print("chargedParticles.pdf");
  
  TCanvas *c2 = new TCanvas("c2","Charged Particles - Angular Distribution",710,710);
  c2->Divide(1,2);
  //
  c2->cd(1);
  theta_chrg->GetXaxis()->SetTitle("#theta");
  theta_chrg->Draw();
  //
  c2->cd(2);
  phi_chrg->GetXaxis()->SetTitle("#phi");
  phi_chrg->Draw();
  //
  c2->Print("charged_angDistribution.pdf");
  
  TCanvas *c3 = new TCanvas("c3","Invariant Mass (Charged)",720,720);
  c3->Divide(1,2);
  //
  c3->cd(1);
  invM_chrg->GetXaxis()->SetTitle("M_{inv} [GeV]");
  invM_chrg->Draw();
  //
  c3->cd(2);
  m_RP_chrg->GetXaxis()->SetTitle("m [GeV]");
  c3->cd(2)->SetLogy();
  m_RP_chrg->Draw();
  //
  c3->Print("invM_charged.pdf");

  TCanvas *c4 = new TCanvas("c4","Neutral Particles",730,730);
  c4->Divide(2,2);
  //
  c4->cd(1);
  n_neut->GetXaxis()->SetTitle("#");
  n_neut->Draw();
  //
  c4->cd(2);
  p_neut->GetXaxis()->SetTitle("|p| [GeV]");
  c4->cd(2)->SetLogy();
  p_neut->Draw();
  //
  c4->cd(3);
  e_neut->GetXaxis()->SetTitle("Energy [GeV]");
  c4->cd(3)->SetLogy();
  e_neut->Draw();
  //
  c4->cd(4);
  pT_neut->GetXaxis()->SetTitle("p_{T} [GeV]");
  c4->cd(4)->SetLogy();
  pT_neut->Draw();
  //
  c4->Print("neutralParticles.pdf");
  
  TCanvas *c5 = new TCanvas("c5","Neutral Particles - Angular Distribution",740,740);
  c5->Divide(1,2);
  //
  c5->cd(1);
  theta_neut->GetXaxis()->SetTitle("#theta");
  theta_neut->Draw();
  //
  c5->cd(2);
  phi_neut->GetXaxis()->SetTitle("#phi");
  phi_neut->Draw();
  //
  c5->Print("neutral_angDistribution.pdf");

  hists->Close();

  return -1;
}
