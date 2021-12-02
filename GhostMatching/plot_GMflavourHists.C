#include <TFile.h>
#include <TH1F.h>
#include <THStack.h>
#include <TStyle.h>
#include <TCanvas.h>

int main()
{
  TFile *hists=TFile::Open("histZuds_gm.root");
  TFile *hists1=TFile::Open("histZuds_excl_st23.root");
  
  TH1F* jetFlavour = (TH1F*)hists->Get("h_jetFlavour");
  TH1F* jetTheta = (TH1F*)hists->Get("h_jetTheta");
  TH1F* jetPhi = (TH1F*)hists->Get("h_jetPhi");
  //
  TH1F* jetFlavour23 = (TH1F*)hists1->Get("h_jetFlavour_qqbar");
  

  TCanvas *c1 = new TCanvas("c1","Jet Flavour : Ghost-Matching",700,700);
  c1->SetLogy();
  c1->SetGrid();
  jetFlavour->GetXaxis()->SetTitle("PDG ID");
  jetFlavour->Draw();
  //
  c1->Print("jetFlavour_gm.pdf");

  TCanvas *c2 = new TCanvas("c2","Jet Flavour Assignment",710,710);
  c2->Divide(2,1);
  //
  c2->cd(1);
  c2->cd(1)->SetLogy();
  c2->cd(1)->SetGrid();
  jetFlavour->SetTitle("Ghost Matching with Status 23");
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
  
  TCanvas *c3 = new TCanvas("c3","Jet Angular Distribution (GM)",720,720);
  c3->Divide(1,2);
  //
  c3->cd(1);
  jetTheta->SetTitle("Jet Polar Angle (#theta)");
  jetTheta->GetXaxis()->SetTitle("#theta");
  jetTheta->Draw();
  //
  c3->cd(2);
  jetPhi->SetTitle("Jet Azimuthal Angle (#phi)");
  jetPhi->GetXaxis()->SetTitle("#phi");
  jetPhi->Draw();
  //
  c3->Print("jetAngularDist_gm.pdf");  
  
  hists->Close();
  hists1->Close();

  return -1;
}
