#include <TFile.h>
#include <TH1F.h>
#include <THStack.h>
#include <TStyle.h>
#include <TCanvas.h>

int main()
{
  TFile *hists=TFile::Open("histZuds_incl.root");

  TH1F* njet = (TH1F*)hists->Get("h_njet");
  TH1F* pjet = (TH1F*)hists->Get("h_pjet");
  TH1F* pTjet = (TH1F*)hists->Get("h_pTjet");
  TH1F* Ejet = (TH1F*)hists->Get("h_Ejet");
  TH1F* jetTheta = (TH1F*)hists->Get("h_jetTheta");
  TH1F* jetPhi = (TH1F*)hists->Get("h_jetPhi");
  TH1F* invMjets = (TH1F*)hists->Get("h_invMjets");
  TH1F* invMjets_b = (TH1F*)hists->Get("h_invMjets_b");
  TH1F* invMjets_c = (TH1F*)hists->Get("h_invMjets_c");
  TH1F* invMjets_s = (TH1F*)hists->Get("h_invMjets_s");
  TH1F* invMjets_u = (TH1F*)hists->Get("h_invMjets_u");
  TH1F* invMjets_d = (TH1F*)hists->Get("h_invMjets_d");
  TH1F* invMjets_bHigh = (TH1F*)hists->Get("h_invMjets_bHigh");
  TH1F* invMjets_cHigh = (TH1F*)hists->Get("h_invMjets_cHigh");
  TH1F* invMjets_sHigh = (TH1F*)hists->Get("h_invMjets_sHigh");
  TH1F* invMjets_uHigh = (TH1F*)hists->Get("h_invMjets_uHigh");
  TH1F* invMjets_dHigh = (TH1F*)hists->Get("h_invMjets_dHigh");
  TH1F* dijet_flavour = (TH1F*)hists->Get("h_dijet_flavour");
  TH1F* invMdijets = (TH1F*)hists->Get("h_invMdijets");
  TH1F* invMdijets_b = (TH1F*)hists->Get("h_invMdijets_b");
  TH1F* invMdijets_c = (TH1F*)hists->Get("h_invMdijets_c");
  TH1F* invMdijets_s = (TH1F*)hists->Get("h_invMdijets_s");
  TH1F* invMdijets_u = (TH1F*)hists->Get("h_invMdijets_u");
  TH1F* invMdijets_d = (TH1F*)hists->Get("h_invMdijets_d");
  TH1F* angJP = (TH1F*)hists->Get("h_angJP");
  TH1F* thetaJP = (TH1F*)hists->Get("h_thetaJP");
  TH1F* phiJP = (TH1F*)hists->Get("h_phiJP");

  TCanvas *c3 = new TCanvas("c3","Di-Jet Flavour (ee-genkt)",720,720);
  dijet_flavour->GetXaxis()->SetTitle("PDG_{q}");
  dijet_flavour->Draw();
  c3->SetLogy();
  gPad->SetTicky();
  //
  c3->Print("dijetFlavour_ee_genkt.pdf");
  
  TCanvas *c4 = new TCanvas("c4","Jets (ee-genkt)",730,730);
  c4->Divide(2,2);
  //
  c4->cd(1);
  njet->GetXaxis()->SetTitle("#");
  njet->Draw();
  //
  c4->cd(2);
  pjet->GetXaxis()->SetTitle("|p| [GeV]");
  //c4->cd(2)->SetLogy();
  pjet->Draw();
  //
  c4->cd(3);
  Ejet->GetXaxis()->SetTitle("Energy [GeV]");
  //c4->cd(3)->SetLogy();
  Ejet->Draw();
  //
  c4->cd(4);
  pTjet->GetXaxis()->SetTitle("p_{T} [GeV]");
  //c4->cd(4)->SetLogy();
  pTjet->Draw();
  //
  c4->Print("jets_ee_genkt.pdf");

  TCanvas *c5 = new TCanvas("c5","Jets - Angular Distribution (ee-genkt)",740,740);
  c5->Divide(1,2);
  //
  c5->cd(1);
  jetTheta->GetXaxis()->SetTitle("#theta");
  //c5->cd(1)->SetLogy();
  jetTheta->Draw();
  //
  c5->cd(2);
  jetPhi->GetXaxis()->SetTitle("#phi");
  //c5->cd(2)->SetLogy();
  jetPhi->Draw();
  //
  c5->Print("jetAngularDist_ee_genkt.pdf");
  
  TCanvas *c6 = new TCanvas("c6","Invariant Mass (jet sum)",750,750);
  invMjets->GetXaxis()->SetTitle("M_{inv} [GeV]");
  invMjets->Draw();
  //
  c6->Print("invM_jetsum_ee_genkt.pdf");

  TCanvas *c7 = new TCanvas("c7","Invariant Mass - uds (stacked)",760,760);
  THStack *hs = new THStack("hs","Invariant Mass - uds");
  gStyle->SetPalette(kOcean);
  //invMjets_s->SetFillColor(kRed);
  hs->Add(invMjets_s);
  //invMjets_u->SetFillColor(kBlue);
  hs->Add(invMjets_u);
  //invMjets_d->SetFillColor(kGreen);
  hs->Add(invMjets_d);
  hs->Add(invMjets_c);
  hs->Add(invMjets_b);
  hs->Draw("pfc");
  hs->GetXaxis()->SetTitle("M_{inv} [GeV]");
  gPad->BuildLegend(0.75,0.75,0.95,0.95,"");
  //
  c7->Print("invM_stack_ee_genkt.pdf");
  
  TCanvas *c8 = new TCanvas("c8","Jet Constituents - Angular Distribution (ee_genkt)",770,770);
  c8->Divide(1,3);
  //
  c8->cd(1);
  angJP->GetXaxis()->SetTitle("[rad]");
  c8->cd(1)->SetLogy();
  angJP->Draw();
  //
  c8->cd(2);
  thetaJP->GetXaxis()->SetTitle("#Delta#theta [rad]");
  c8->cd(2)->SetLogy();
  thetaJP->Draw();
  //
  c8->cd(3);
  phiJP->GetXaxis()->SetTitle("#Delta#phi [rad]");
  c8->cd(3)->SetLogy();
  phiJP->Draw();
  //
  c8->Print("jetConstAngularDist_ee_genkt.pdf");
  
  TCanvas *c9 = new TCanvas("c9","Invariant Mass - uds (stacked)",760,760);
  THStack *hs1 = new THStack("hs1","Invariant Mass - uds [highest]");
  gStyle->SetPalette(kOcean);
  //invMjets_sHigh->SetFillColor(kRed);
  hs1->Add(invMjets_sHigh);
  //invMjets_uHigh->SetFillColor(kBlue);
  hs1->Add(invMjets_uHigh);
  //invMjets_dHigh->SetFillColor(kGreen);
  hs1->Add(invMjets_dHigh);
  hs1->Add(invMjets_cHigh);
  hs1->Add(invMjets_bHigh);
  hs1->Draw("pfc");
  hs1->GetXaxis()->SetTitle("M_{inv} [GeV]");
  gPad->BuildLegend(0.75,0.75,0.95,0.95,"");
  //
  c9->Print("invM_stack_ee_genkt_highest.pdf");

  TCanvas *c2 = new TCanvas("c2","Invariant Mass (dijet sum)",710,710);
  invMdijets->GetXaxis()->SetTitle("M_{inv} [GeV]");
  invMdijets->Draw();
  //
  c2->Print("invM_dijetsum_ee_genkt.pdf");

  TCanvas *c1 = new TCanvas("c1","Invariant Mass - uds (stacked)",700,700);
  THStack *hs2 = new THStack("hs2","Invariant Mass - uds [dijet]");
  gStyle->SetPalette(kOcean);
  //invMjets_sHigh->SetFillColor(kRed);
  hs2->Add(invMdijets_s);
  //invMjets_uHigh->SetFillColor(kBlue);
  hs2->Add(invMdijets_u);
  //invMjets_dHigh->SetFillColor(kGreen);
  hs2->Add(invMdijets_d);
  hs2->Add(invMdijets_c);
  hs2->Add(invMdijets_b);
  hs2->Draw("pfc");
  hs2->GetXaxis()->SetTitle("M_{inv} [GeV]");
  gPad->BuildLegend(0.75,0.75,0.95,0.95,"");
  //
  c1->Print("invM_stack_ee_genkt_dijet.pdf");
  
  hists->Close();

  return -1;
}
