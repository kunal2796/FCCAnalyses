#include <TFile.h>
#include <TH1F.h>
#include <THStack.h>
#include <TStyle.h>
#include <TCanvas.h>

int main()
{
  TFile *hists_gm =TFile::Open("histZuds_gm.root");
  TFile *hists_reg=TFile::Open("histZuds_excl.root");
  
  TH1F* jetTheta_gm  = (TH1F*)hists_gm->Get("h_jetTheta");
  TH1F* jetPhi_gm    = (TH1F*)hists_gm->Get("h_jetPhi");
  TH1F* pjet_gm      = (TH1F*)hists_gm->Get("h_pjet");
  TH1F* pTjet_gm     = (TH1F*)hists_gm->Get("h_pTjet");
  //
  TH1F* jetTheta_reg = (TH1F*)hists_reg->Get("h_jetTheta");
  TH1F* jetPhi_reg   = (TH1F*)hists_reg->Get("h_jetPhi");
  TH1F* pjet_reg     = (TH1F*)hists_reg->Get("h_pjet");
  TH1F* pTjet_reg    = (TH1F*)hists_reg->Get("h_pTjet");
  //
  TH1F* theta_comp = new TH1F("theta_comp", "Jet #theta Comparison [w/ - w/o GM]",100,0,3.15);
  TH1F* phi_comp   = new TH1F("phi_comp", "Jet #phi Comparison [w/ - w/o GM]",100,-3.15,3.15);
  TH1F* p_comp     = new TH1F("p_comp", "Jet Momentum Comparison [w/ - w/o GM]",100,0,50);
  TH1F* pT_comp    = new TH1F("pT_comp", "Jet Transverse Momentum Comparison [w/ - w/o GM]",100,0,50);
  
  //

  // calculation of comparisons - bin by bin
  float theta_gm[100], theta_reg[100], phi_gm[100], phi_reg[100];
  float p_gm[100], p_reg[100], pT_gm[100], pT_reg[100];
  for(unsigned int i = 0; i < 100; i++)
    {
      // angular distribution
      theta_gm[i]  = jetTheta_gm->GetBinContent(i);
      theta_reg[i] = jetTheta_reg->GetBinContent(i);
      phi_gm[i]    = jetPhi_gm->GetBinContent(i);
      phi_reg[i]   = jetPhi_reg->GetBinContent(i);
      //
      theta_comp->SetBinContent(i, theta_gm[i] - theta_reg[i]);
      phi_comp->SetBinContent(i, phi_gm[i] - phi_reg[i]);

      // momenta
      p_gm[i]   = pjet_gm->GetBinContent(i);
      p_reg[i]  = pjet_reg->GetBinContent(i);
      pT_gm[i]  = pTjet_gm->GetBinContent(i);
      pT_reg[i] = pTjet_reg->GetBinContent(i);
      //
      p_comp->SetBinContent(i, p_gm[i] - p_reg[i]);
      pT_comp->SetBinContent(i, pT_gm[i] - pT_reg[i]);
    }
  
  //
  
  TCanvas *c1 = new TCanvas("c1","Jet Angular Distribution (Comparison)",700,700);
  c1->Divide(1,2);
  //
  c1->cd(1);
  theta_comp->GetXaxis()->SetTitle("#theta");
  theta_comp->Draw();
  //
  c1->cd(2);
  phi_comp->GetXaxis()->SetTitle("#phi");
  phi_comp->Draw();
  //
  c1->Print("comp_jetAngularDist.pdf");  

  TCanvas *c2 = new TCanvas("c2","Jet Momentum (Comparison)",710,710);
  c2->Divide(1,2);
  //
  c2->cd(1);
  p_comp->GetXaxis()->SetTitle("p [GeV]");
  p_comp->Draw();
  //
  c2->cd(2);
  pT_comp->GetXaxis()->SetTitle("p_T [GeV]");
  pT_comp->Draw();
  //
  c2->Print("comp_jetMomentum_comp.pdf");  

  hists_gm->Close();
  hists_reg->Close();

  return -1;
}
