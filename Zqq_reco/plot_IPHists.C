#include <TFile.h>
#include <TH1F.h>
#include <THStack.h>
#include <TStyle.h>
#include <TCanvas.h>

int main()
{
  TFile *hists=TFile::Open("histZuds_IP.root");
  
  TH1F* d0 = (TH1F*)hists->Get("h_d0");
  TH1F* z0 = (TH1F*)hists->Get("h_z0");
  TH1F* d0sig = (TH1F*)hists->Get("h_d0sig");
  TH1F* z0sig = (TH1F*)hists->Get("h_z0sig");
  TH1F* d0_chrg = (TH1F*)hists->Get("h_d0_chrg");
  TH1F* z0_chrg = (TH1F*)hists->Get("h_z0_chrg");
  TH1F* d0sig_chrg = (TH1F*)hists->Get("h_d0sig_chrg");
  TH1F* z0sig_chrg = (TH1F*)hists->Get("h_z0sig_chrg");
  TH1F* d0_neut = (TH1F*)hists->Get("h_d0_neut");
  TH1F* z0_neut = (TH1F*)hists->Get("h_z0_neut");
  TH1F* d0sig_neut = (TH1F*)hists->Get("h_d0sig_neut");
  TH1F* z0sig_neut = (TH1F*)hists->Get("h_z0sig_neut");
  TH1F* d0_b = (TH1F*)hists->Get("h_d0_b");
  TH1F* z0_b = (TH1F*)hists->Get("h_z0_b");
  TH1F* d0sig_b = (TH1F*)hists->Get("h_d0sig_b");
  TH1F* z0sig_b = (TH1F*)hists->Get("h_z0sig_b");
  TH1F* d0_c = (TH1F*)hists->Get("h_d0_c");
  TH1F* z0_c = (TH1F*)hists->Get("h_z0_c");
  TH1F* d0sig_c = (TH1F*)hists->Get("h_d0sig_c");
  TH1F* z0sig_c = (TH1F*)hists->Get("h_z0sig_c");
  TH1F* d0_s = (TH1F*)hists->Get("h_d0_s");
  TH1F* z0_s = (TH1F*)hists->Get("h_z0_s");
  TH1F* d0sig_s = (TH1F*)hists->Get("h_d0sig_s");
  TH1F* z0sig_s = (TH1F*)hists->Get("h_z0sig_s");
  TH1F* d0_u = (TH1F*)hists->Get("h_d0_u");
  TH1F* z0_u = (TH1F*)hists->Get("h_z0_u");
  TH1F* d0sig_u = (TH1F*)hists->Get("h_d0sig_u");
  TH1F* z0sig_u = (TH1F*)hists->Get("h_z0sig_u");
  TH1F* d0_d = (TH1F*)hists->Get("h_d0_d");
  TH1F* z0_d = (TH1F*)hists->Get("h_z0_d");
  TH1F* d0sig_d = (TH1F*)hists->Get("h_d0sig_d");
  TH1F* z0sig_d = (TH1F*)hists->Get("h_z0sig_d");
  
  TCanvas *c1 = new TCanvas("c1","Reco Tracks",700,700);
  c1->Divide(2,2);
  //
  c1->cd(1);
  d0->GetXaxis()->SetTitle("|D_{0}| [mm]");
  c1->cd(1)->SetLogy();
  d0->Draw();
  //
  c1->cd(2);
  z0->GetXaxis()->SetTitle("|z_{0}| [mm]");
  c1->cd(2)->SetLogy();
  z0->Draw();
  //
  c1->cd(3);
  d0sig->GetXaxis()->SetTitle("D_{0}/#sigma_{D_{0}}");
  c1->cd(3)->SetLogy();
  d0sig->Draw();
  //
  c1->cd(4);
  z0sig->GetXaxis()->SetTitle("z_{0}/#sigma_{z_{0}}");
  c1->cd(4)->SetLogy();
  z0sig->Draw();
  //
  c1->Print("imp_para.pdf");

  TCanvas *c2 = new TCanvas("c2","Reco Tracks (stacked)",760,760);
  c2->Divide(2,2);
  //
  c2->cd(1);
  c2->cd(1)->SetLogy();
  THStack *hs1 = new THStack("hs1","Transverse Impact Parameter");
  gStyle->SetPalette(kOcean);
  hs1->Add(d0_s);
  hs1->Add(d0_u);
  hs1->Add(d0_d);
  hs1->Add(d0_c);
  hs1->Add(d0_b);
  hs1->Draw("pfc");
  hs1->GetXaxis()->SetTitle("|D_{0}| [mm]");
  gPad->BuildLegend(0.75,0.75,0.95,0.95,"");
  //
  c2->cd(2);
  c2->cd(2)->SetLogy();
  THStack *hs2 = new THStack("hs2","Longitudinal Impact Parameter");
  //gStyle->SetPalette(kOcean);
  hs2->Add(z0_s);
  hs2->Add(z0_u);
  hs2->Add(z0_d);
  hs2->Add(z0_c);
  hs2->Add(z0_b);
  hs2->Draw("pfc");
  hs2->GetXaxis()->SetTitle("|z_{0}| [mm]");
  gPad->BuildLegend(0.75,0.75,0.95,0.95,"");
  //
  c2->cd(3);
  c2->cd(3)->SetLogy();
  THStack *hs3 = new THStack("hs3","Transverse Impact Parameter Significance");
  //gStyle->SetPalette(kOcean);
  hs3->Add(d0sig_s);
  hs3->Add(d0sig_u);
  hs3->Add(d0sig_d);
  hs3->Add(d0sig_c);
  hs3->Add(d0sig_b);
  hs3->Draw("pfc");
  hs3->GetXaxis()->SetTitle("D_{0}/#sigma_{D_{0}}");
  gPad->BuildLegend(0.75,0.75,0.95,0.95,"");
  //
  c2->cd(4);
  c2->cd(4)->SetLogy();
  THStack *hs4 = new THStack("hs4","Longitudinal Impact Parameter Significance");
  //gStyle->SetPalette(kOcean);
  hs4->Add(z0sig_s);
  hs4->Add(z0sig_u);
  hs4->Add(z0sig_d);
  hs4->Add(z0sig_c);
  hs4->Add(z0sig_b);
  hs4->Draw("pfc");
  hs4->GetXaxis()->SetTitle("z_{0}/#sigma_{z_{0}}");
  gPad->BuildLegend(0.75,0.75,0.95,0.95,"");
  //
  c2->Print("imp_para_stacked_log.pdf");

  TCanvas *c3 = new TCanvas("c3","Reco Tracks (Charged)",720,720);
  c3->Divide(2,2);
  //
  c3->cd(1);
  d0_chrg->GetXaxis()->SetTitle("|D_{0}| [mm]");
  c3->cd(1)->SetLogy();
  d0_chrg->Draw();
  //
  c3->cd(2);
  z0_chrg->GetXaxis()->SetTitle("|z_{0}| [mm]");
  c3->cd(2)->SetLogy();
  z0_chrg->Draw();
  //
  c3->cd(3);
  d0sig_chrg->GetXaxis()->SetTitle("D_{0}/#sigma_{D_{0}}");
  c3->cd(3)->SetLogy();
  d0sig_chrg->Draw();
  //
  c3->cd(4);
  z0sig_chrg->GetXaxis()->SetTitle("z_{0}/#sigma_{z_{0}}");
  c3->cd(4)->SetLogy();
  z0sig_chrg->Draw();
  //
  c3->Print("imp_para_chrgd.pdf");

  TCanvas *c4 = new TCanvas("c4","Reco Tracks (Neutral)",730,730);
  c4->Divide(2,2);
  //
  c4->cd(1);
  d0_neut->GetXaxis()->SetTitle("|D_{0}| [mm]");
  c4->cd(1)->SetLogy();
  d0_neut->Draw();
  //
  c4->cd(2);
  z0_neut->GetXaxis()->SetTitle("|z_{0}| [mm]");
  c4->cd(2)->SetLogy();
  z0_neut->Draw();
  //
  c4->cd(3);
  d0sig_neut->GetXaxis()->SetTitle("D_{0}/#sigma_{D_{0}}");
  c4->cd(3)->SetLogy();
  d0sig_neut->Draw();
  //
  c4->cd(4);
  z0sig_neut->GetXaxis()->SetTitle("z_{0}/#sigma_{z_{0}}");
  c4->cd(4)->SetLogy();
  z0sig_neut->Draw();
  //
  c4->Print("imp_para_neut.pdf");
  
  hists->Close();

  return -1;
}
