#include <TFile.h>
#include <TH1F.h>
#include <THStack.h>
#include <TStyle.h>
#include <TCanvas.h>

int main()
{
  TFile *hists=TFile::Open("histZbb_gm_chrgImb.root");

  //
  
  TH1F* p_q    = (TH1F*)hists->Get("h_p_q");
  TH1F* p_qbar = (TH1F*)hists->Get("h_p_qbar");
  TH1F* theta_q    = (TH1F*)hists->Get("h_theta_q");
  TH1F* theta_qbar = (TH1F*)hists->Get("h_theta_qbar");
  //
  TH1F* p_p5p1     = (TH1F*)hists->Get("h_p_p5p1");
  TH1F* theta_p5p1 = (TH1F*)hists->Get("h_theta_p5p1");
  TH1F* p_p5m1     = (TH1F*)hists->Get("h_p_p5m1");
  TH1F* theta_p5m1 = (TH1F*)hists->Get("h_theta_p5m1");
  TH1F* p_p5p2     = (TH1F*)hists->Get("h_p_p5p2");
  TH1F* theta_p5p2 = (TH1F*)hists->Get("h_theta_p5p2");
  TH1F* p_p5m2     = (TH1F*)hists->Get("h_p_p5m2");
  TH1F* theta_p5m2 = (TH1F*)hists->Get("h_theta_p5m2");
  TH1F* p_p5p3     = (TH1F*)hists->Get("h_p_p5p3");
  TH1F* theta_p5p3 = (TH1F*)hists->Get("h_theta_p5p3");
  TH1F* p_p5m3     = (TH1F*)hists->Get("h_p_p5m3");
  TH1F* theta_p5m3 = (TH1F*)hists->Get("h_theta_p5m3");
  TH1F* p_p5p4     = (TH1F*)hists->Get("h_p_p5p4");
  TH1F* theta_p5p4 = (TH1F*)hists->Get("h_theta_p5p4");
  TH1F* p_p5m4     = (TH1F*)hists->Get("h_p_p5m4");
  TH1F* theta_p5m4 = (TH1F*)hists->Get("h_theta_p5m4");
  TH1F* p_p5p5     = (TH1F*)hists->Get("h_p_p5p5");
  TH1F* theta_p5p5 = (TH1F*)hists->Get("h_theta_p5p5");
  TH1F* p_p5m5     = (TH1F*)hists->Get("h_p_p5m5");
  TH1F* theta_p5m5 = (TH1F*)hists->Get("h_theta_p5m5");
  //
  TH1F* p_m5p1     = (TH1F*)hists->Get("h_p_m5p1");
  TH1F* theta_m5p1 = (TH1F*)hists->Get("h_theta_m5p1");
  TH1F* p_m5m1     = (TH1F*)hists->Get("h_p_m5m1");
  TH1F* theta_m5m1 = (TH1F*)hists->Get("h_theta_m5m1");
  TH1F* p_m5p2     = (TH1F*)hists->Get("h_p_m5p2");
  TH1F* theta_m5p2 = (TH1F*)hists->Get("h_theta_m5p2");
  TH1F* p_m5m2     = (TH1F*)hists->Get("h_p_m5m2");
  TH1F* theta_m5m2 = (TH1F*)hists->Get("h_theta_m5m2");
  TH1F* p_m5p3     = (TH1F*)hists->Get("h_p_m5p3");
  TH1F* theta_m5p3 = (TH1F*)hists->Get("h_theta_m5p3");
  TH1F* p_m5m3     = (TH1F*)hists->Get("h_p_m5m3");
  TH1F* theta_m5m3 = (TH1F*)hists->Get("h_theta_m5m3");
  TH1F* p_m5p4     = (TH1F*)hists->Get("h_p_m5p4");
  TH1F* theta_m5p4 = (TH1F*)hists->Get("h_theta_m5p4");
  TH1F* p_m5m4     = (TH1F*)hists->Get("h_p_m5m4");
  TH1F* theta_m5m4 = (TH1F*)hists->Get("h_theta_m5m4");
  TH1F* p_m5p5     = (TH1F*)hists->Get("h_p_m5p5");
  TH1F* theta_m5p5 = (TH1F*)hists->Get("h_theta_m5p5");
  TH1F* p_m5m5     = (TH1F*)hists->Get("h_p_m5m5");
  TH1F* theta_m5m5 = (TH1F*)hists->Get("h_theta_m5m5");
  //
  TH1F* p_0p1     = (TH1F*)hists->Get("h_p_0p1");
  TH1F* theta_0p1 = (TH1F*)hists->Get("h_theta_0p1");
  TH1F* p_0m1     = (TH1F*)hists->Get("h_p_0m1");
  TH1F* theta_0m1 = (TH1F*)hists->Get("h_theta_0m1");
  TH1F* p_0p2     = (TH1F*)hists->Get("h_p_0p2");
  TH1F* theta_0p2 = (TH1F*)hists->Get("h_theta_0p2");
  TH1F* p_0m2     = (TH1F*)hists->Get("h_p_0m2");
  TH1F* theta_0m2 = (TH1F*)hists->Get("h_theta_0m2");
  TH1F* p_0p3     = (TH1F*)hists->Get("h_p_0p3");
  TH1F* theta_0p3 = (TH1F*)hists->Get("h_theta_0p3");
  TH1F* p_0m3     = (TH1F*)hists->Get("h_p_0m3");
  TH1F* theta_0m3 = (TH1F*)hists->Get("h_theta_0m3");
  TH1F* p_0p4     = (TH1F*)hists->Get("h_p_0p4");
  TH1F* theta_0p4 = (TH1F*)hists->Get("h_theta_0p4");
  TH1F* p_0m4     = (TH1F*)hists->Get("h_p_0m4");
  TH1F* theta_0m4 = (TH1F*)hists->Get("h_theta_0m4");
  TH1F* p_0p5     = (TH1F*)hists->Get("h_p_0p5");
  TH1F* theta_0p5 = (TH1F*)hists->Get("h_theta_0p5");
  TH1F* p_0m5     = (TH1F*)hists->Get("h_p_0m5");
  TH1F* theta_0m5 = (TH1F*)hists->Get("h_theta_0m5");
  //
  TH1F* p_p50     = (TH1F*)hists->Get("h_p_p50");
  TH1F* theta_p50 = (TH1F*)hists->Get("h_theta_p50");
  TH1F* p_m50     = (TH1F*)hists->Get("h_p_m50");
  TH1F* theta_m50 = (TH1F*)hists->Get("h_theta_m50");
  TH1F* p_00      = (TH1F*)hists->Get("h_p_00");
  TH1F* theta_00  = (TH1F*)hists->Get("h_theta_00");

  //

  TCanvas *c1 = new TCanvas("c1","Charge Imbalance",700,700);
  c1->Divide(2,2);
  //
  c1->cd(1);
  p_q->SetTitle("q Jets");
  p_q->GetXaxis()->SetTitle("|p| [GeV]");
  p_q->Draw();
  //
  c1->cd(2);
  p_qbar->SetTitle("#bar{q} Jets");
  p_qbar->GetXaxis()->SetTitle("|p| [GeV]");
  p_qbar->Draw();
  //
  c1->cd(3);
  theta_q->SetTitle("q Jets");
  theta_q->GetXaxis()->SetTitle("#theta");
  theta_q->Draw();
  //
  c1->cd(4);
  theta_qbar->SetTitle("#bar{q} Jets");
  theta_qbar->GetXaxis()->SetTitle("#theta");
  theta_qbar->Draw();
  //
  c1->Print("chrgImb_case2_7xNon0.pdf");

  TCanvas *c2 = new TCanvas("c2","23 -> +5 && 7x -> +-1",700,700);
  c2->Divide(2,2);
  //
  c2->cd(1);
  p_p5p1->GetXaxis()->SetTitle("|p| [GeV]");
  p_p5p1->Draw();
  //
  c2->cd(2);
  p_p5m1->GetXaxis()->SetTitle("|p| [GeV]");
  p_p5m1->Draw();
  //
  c2->cd(3);
  theta_p5p1->GetXaxis()->SetTitle("#theta");
  theta_p5p1->Draw();
  //
  c2->cd(4);
  theta_p5m1->GetXaxis()->SetTitle("#theta");
  theta_p5m1->Draw();
  //
  c2->Print("p5pm1.pdf");

  TCanvas *c3 = new TCanvas("c3","23 -> +5 && 7x -> +-2",700,700);
  c3->Divide(2,2);
  //
  c3->cd(1);
  p_p5p2->GetXaxis()->SetTitle("|p| [GeV]");
  p_p5p2->Draw();
  //
  c3->cd(2);
  p_p5m2->GetXaxis()->SetTitle("|p| [GeV]");
  p_p5m2->Draw();
  //
  c3->cd(3);
  theta_p5p2->GetXaxis()->SetTitle("#theta");
  theta_p5p2->Draw();
  //
  c3->cd(4);
  theta_p5m2->GetXaxis()->SetTitle("#theta");
  theta_p5m2->Draw();
  //
  c3->Print("p5pm2.pdf");

  TCanvas *c4 = new TCanvas("c4","23 -> +5 && 7x -> +-3",700,700);
  c4->Divide(2,2);
  //
  c4->cd(1);
  p_p5p3->GetXaxis()->SetTitle("|p| [GeV]");
  p_p5p3->Draw();
  //
  c4->cd(2);
  p_p5m3->GetXaxis()->SetTitle("|p| [GeV]");
  p_p5m3->Draw();
  //
  c4->cd(3);
  theta_p5p3->GetXaxis()->SetTitle("#theta");
  theta_p5p3->Draw();
  //
  c4->cd(4);
  theta_p5m3->GetXaxis()->SetTitle("#theta");
  theta_p5m3->Draw();
  //
  c4->Print("p5pm3.pdf");

  TCanvas *c5 = new TCanvas("c5","23 -> +5 && 7x -> +-4",700,700);
  c5->Divide(2,2);
  //
  c5->cd(1);
  p_p5p4->GetXaxis()->SetTitle("|p| [GeV]");
  p_p5p4->Draw();
  //
  c5->cd(2);
  p_p5m4->GetXaxis()->SetTitle("|p| [GeV]");
  p_p5m4->Draw();
  //
  c5->cd(3);
  theta_p5p4->GetXaxis()->SetTitle("#theta");
  theta_p5p4->Draw();
  //
  c5->cd(4);
  theta_p5m4->GetXaxis()->SetTitle("#theta");
  theta_p5m4->Draw();
  //
  c5->Print("p5pm4.pdf");

  TCanvas *c6 = new TCanvas("c6","23 -> -5 && 7x -> +-1",700,700);
  c6->Divide(2,2);
  //
  c6->cd(1);
  p_m5p1->GetXaxis()->SetTitle("|p| [GeV]");
  p_m5p1->Draw();
  //
  c6->cd(2);
  p_m5m1->GetXaxis()->SetTitle("|p| [GeV]");
  p_m5m1->Draw();
  //
  c6->cd(3);
  theta_m5p1->GetXaxis()->SetTitle("#theta");
  theta_m5p1->Draw();
  //
  c6->cd(4);
  theta_m5m1->GetXaxis()->SetTitle("#theta");
  theta_m5m1->Draw();
  //
  c6->Print("m5pm1.pdf");

  TCanvas *c7 = new TCanvas("c7","23 -> -5 && 7x -> +-2",700,700);
  c7->Divide(2,2);
  //
  c7->cd(1);
  p_m5p2->GetXaxis()->SetTitle("|p| [GeV]");
  p_m5p2->Draw();
  //
  c7->cd(2);
  p_m5m2->GetXaxis()->SetTitle("|p| [GeV]");
  p_m5m2->Draw();
  //
  c7->cd(3);
  theta_m5p2->GetXaxis()->SetTitle("#theta");
  theta_m5p2->Draw();
  //
  c7->cd(4);
  theta_m5m2->GetXaxis()->SetTitle("#theta");
  theta_m5m2->Draw();
  //
  c7->Print("m5pm2.pdf");

  TCanvas *c8 = new TCanvas("c8","23 -> -5 && 7x -> +-3",700,700);
  c8->Divide(2,2);
  //
  c8->cd(1);
  p_m5p3->GetXaxis()->SetTitle("|p| [GeV]");
  p_m5p3->Draw();
  //
  c8->cd(2);
  p_m5m3->GetXaxis()->SetTitle("|p| [GeV]");
  p_m5m3->Draw();
  //
  c8->cd(3);
  theta_m5p3->GetXaxis()->SetTitle("#theta");
  theta_m5p3->Draw();
  //
  c8->cd(4);
  theta_m5m3->GetXaxis()->SetTitle("#theta");
  theta_m5m3->Draw();
  //
  c8->Print("m5pm3.pdf");

  TCanvas *c9 = new TCanvas("c9","23 -> -5 && 7x -> +-4",700,700);
  c9->Divide(2,2);
  //
  c9->cd(1);
  p_m5p4->GetXaxis()->SetTitle("|p| [GeV]");
  p_m5p4->Draw();
  //
  c9->cd(2);
  p_m5m4->GetXaxis()->SetTitle("|p| [GeV]");
  p_m5m4->Draw();
  //
  c9->cd(3);
  theta_m5p4->GetXaxis()->SetTitle("#theta");
  theta_m5p4->Draw();
  //
  c9->cd(4);
  theta_m5m4->GetXaxis()->SetTitle("#theta");
  theta_m5m4->Draw();
  //
  c9->Print("m5pm4.pdf");

  TCanvas *d1 = new TCanvas("d1","23 -> +5 && 7x -> +-5",700,700);
  d1->Divide(2,2);
  //
  d1->cd(1);
  p_p5p5->GetXaxis()->SetTitle("|p| [GeV]");
  p_p5p5->Draw();
  //
  d1->cd(2);
  p_p5m5->GetXaxis()->SetTitle("|p| [GeV]");
  p_p5m5->Draw();
  //
  d1->cd(3);
  theta_p5p5->GetXaxis()->SetTitle("#theta");
  theta_p5p5->Draw();
  //
  d1->cd(4);
  theta_p5m5->GetXaxis()->SetTitle("#theta");
  theta_p5m5->Draw();
  //
  d1->Print("p5pm5.pdf");

  TCanvas *d2 = new TCanvas("d2","23 -> -5 && 7x -> +-5",700,700);
  d2->Divide(2,2);
  //
  d2->cd(1);
  p_m5p5->GetXaxis()->SetTitle("|p| [GeV]");
  p_m5p5->Draw();
  //
  d2->cd(2);
  p_m5m5->GetXaxis()->SetTitle("|p| [GeV]");
  p_m5m5->Draw();
  //
  d2->cd(3);
  theta_m5p5->GetXaxis()->SetTitle("#theta");
  theta_m5p5->Draw();
  //
  d2->cd(4);
  theta_m5m5->GetXaxis()->SetTitle("#theta");
  theta_m5m5->Draw();
  //
  d2->Print("m5pm5.pdf");

  TCanvas *e1 = new TCanvas("e1","23 -> 0 && 7x -> +-1",700,700);
  e1->Divide(2,2);
  //
  e1->cd(1);
  p_0p1->GetXaxis()->SetTitle("|p| [GeV]");
  p_0p1->Draw();
  //
  e1->cd(2);
  p_0m1->GetXaxis()->SetTitle("|p| [GeV]");
  p_0m1->Draw();
  //
  e1->cd(3);
  theta_0p1->GetXaxis()->SetTitle("#theta");
  theta_0p1->Draw();
  //
  e1->cd(4);
  theta_0m1->GetXaxis()->SetTitle("#theta");
  theta_0m1->Draw();
  //
  e1->Print("0pm1.pdf");

  TCanvas *e2 = new TCanvas("e2","23 -> 0 && 7x -> +-2",700,700);
  e2->Divide(2,2);
  //
  e2->cd(1);
  p_0p2->GetXaxis()->SetTitle("|p| [GeV]");
  p_0p2->Draw();
  //
  e2->cd(2);
  p_0m2->GetXaxis()->SetTitle("|p| [GeV]");
  p_0m2->Draw();
  //
  e2->cd(3);
  theta_0p2->GetXaxis()->SetTitle("#theta");
  theta_0p2->Draw();
  //
  e2->cd(4);
  theta_0m2->GetXaxis()->SetTitle("#theta");
  theta_0m2->Draw();
  //
  e2->Print("0pm2.pdf");

  TCanvas *e3 = new TCanvas("e3","23 -> 0 && 7x -> +-3",700,700);
  e3->Divide(2,2);
  //
  e3->cd(1);
  p_0p3->GetXaxis()->SetTitle("|p| [GeV]");
  p_0p3->Draw();
  //
  e3->cd(2);
  p_0m3->GetXaxis()->SetTitle("|p| [GeV]");
  p_0m3->Draw();
  //
  e3->cd(3);
  theta_0p3->GetXaxis()->SetTitle("#theta");
  theta_0p3->Draw();
  //
  e3->cd(4);
  theta_0m3->GetXaxis()->SetTitle("#theta");
  theta_0m3->Draw();
  //
  e3->Print("0pm3.pdf");

  TCanvas *e4 = new TCanvas("e4","23 -> 0 && 7x -> +-4",700,700);
  e4->Divide(2,2);
  //
  e4->cd(1);
  p_0p4->GetXaxis()->SetTitle("|p| [GeV]");
  p_0p4->Draw();
  //
  e4->cd(2);
  p_0m4->GetXaxis()->SetTitle("|p| [GeV]");
  p_0m4->Draw();
  //
  e4->cd(3);
  theta_0p4->GetXaxis()->SetTitle("#theta");
  theta_0p4->Draw();
  //
  e4->cd(4);
  theta_0m4->GetXaxis()->SetTitle("#theta");
  theta_0m4->Draw();
  //
  e4->Print("0pm4.pdf");

  TCanvas *e5 = new TCanvas("e5","23 -> 0 && 7x -> +-5",700,700);
  e5->Divide(2,2);
  //
  e5->cd(1);
  p_0p5->GetXaxis()->SetTitle("|p| [GeV]");
  p_0p5->Draw();
  //
  e5->cd(2);
  p_0m5->GetXaxis()->SetTitle("|p| [GeV]");
  p_0m5->Draw();
  //
  e5->cd(3);
  theta_0p5->GetXaxis()->SetTitle("#theta");
  theta_0p5->Draw();
  //
  e5->cd(4);
  theta_0m5->GetXaxis()->SetTitle("#theta");
  theta_0m5->Draw();
  //
  e5->Print("0pm5.pdf");

  TCanvas *e6 = new TCanvas("e6","23 -> +-5 && 7x -> 0",700,700);
  e6->Divide(2,2);
  //
  e6->cd(1);
  p_p50->GetXaxis()->SetTitle("|p| [GeV]");
  p_p50->Draw();
  //
  e6->cd(2);
  p_m50->GetXaxis()->SetTitle("|p| [GeV]");
  p_m50->Draw();
  //
  e6->cd(3);
  theta_p50->GetXaxis()->SetTitle("#theta");
  theta_p50->Draw();
  //
  e6->cd(4);
  theta_m50->GetXaxis()->SetTitle("#theta");
  theta_m50->Draw();
  //
  e6->Print("pm50.pdf");

  TCanvas *e7 = new TCanvas("e7","23 -> 0 && 7x -> 0",700,700);
  e7->Divide(1,2);
  //
  e7->cd(1);
  p_00->GetXaxis()->SetTitle("|p| [GeV]");
  p_00->Draw();
  //
  e7->cd(2);
  theta_00->GetXaxis()->SetTitle("#theta");
  theta_00->Draw();
  //
  e7->Print("00.pdf");
  
  //
  
  hists->Close();
  
  return -1;
}
