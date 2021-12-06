#include <TFile.h>
#include <TH1F.h>
#include <THStack.h>
#include <TStyle.h>
#include <TCanvas.h>

int main()
{
  TFile *hists=TFile::Open("histZqq_gm.root");
  
  TH1F* jetFlavour = (TH1F*)hists->Get("h_jetFlavour");
  

  TCanvas *c1 = new TCanvas("c1","Jet Flavour : Ghost-Matching",700,700);
  c1->SetLogy();
  c1->SetGrid();
  jetFlavour->SetMinimum(60000);
  jetFlavour->SetMaximum(350000);
  jetFlavour->GetXaxis()->SetTitle("PDG ID");
  jetFlavour->Draw();
  //
  //c1->Print("jetFlavour_gm_qq240.pdf");
  c1->Print("jetFlavour_gm_qq365.pdf");
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
  jetFlavour->SetMaximum(700000);
  jetFlavour->Draw();
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
  c2->Print("jetFlavour_comparison.pdf");  
  */
  hists->Close();
  
  return -1;
}
