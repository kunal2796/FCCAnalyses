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
  /*
  TCanvas *c2 = new TCanvas("c2","Jet Flavour Assignment",700,700);
  c2->Divide(2,1);
  //
  c2->cd(1);
  c2->cd(1)->SetLogy();
  c2->cd(1)->SetGrid();
  jetFlavour->SetTitle("Ghost Matching");
  jetFlavour->GetXaxis()->SetTitle("PDG ID");
  jetFlavour->SetMinimum(300);
  jetFlavour->SetMinimum(700000);
  jetFlavour->Draw();
  //
  c2->cd(2);
  c2->cd(2)->SetLogy();
  c2->cd(2)->SetGrid();
  jetFlavour23->SetTitle("Cone Algo with Status 23");
  jetFlavour23->GetXaxis()->SetTitle("PDG ID");
  jetFlavour23->SetMinimum(300);
  jetFlavour23->SetMinimum(700000);
  jetFlavour23->Draw();
  //
  c2->Print("jetFlavour_comparison.pdf");  
  */
  hists->Close();
  //hists1->Close();
  
  return -1;
}
