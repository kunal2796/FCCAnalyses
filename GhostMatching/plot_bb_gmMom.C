#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <THStack.h>
#include <TStyle.h>
#include <TCanvas.h>

int main()
{
  TFile *hists=TFile::Open("histZbb_partonMom.root");

  //
  
  TH2F* mom7x23      = (TH2F*)hists->Get("h_mom7x23");
  TH2F* mom7x23_non0 = (TH2F*)hists->Get("h_mom7x23_non0");
  //
  TH1F* mom7xb = (TH1F*)hists->Get("h_mom7xb");
  TH1F* mom7xc = (TH1F*)hists->Get("h_mom7xc");
  TH1F* mom7xs = (TH1F*)hists->Get("h_mom7xs");
  TH1F* mom7xu = (TH1F*)hists->Get("h_mom7xu");
  TH1F* mom7xd = (TH1F*)hists->Get("h_mom7xd");
  TH1F* mom7x0 = (TH1F*)hists->Get("h_mom7x0");
  TH1F* mom7xbbar = (TH1F*)hists->Get("h_mom7xbbar");
  TH1F* mom7xcbar = (TH1F*)hists->Get("h_mom7xcbar");
  TH1F* mom7xsbar = (TH1F*)hists->Get("h_mom7xsbar");
  TH1F* mom7xubar = (TH1F*)hists->Get("h_mom7xubar");
  TH1F* mom7xdbar = (TH1F*)hists->Get("h_mom7xdbar");

  //

  TCanvas *c1 = new TCanvas("c1","Zbb : b",700,700);
  c1->Divide(1,2);
  //
  c1->cd(1);
  mom7xb->GetXaxis()->SetTitle("|p| [GeV]");
  mom7xb->Draw();
  //
  c1->cd(2);
  mom7xbbar->GetXaxis()->SetTitle("|p| [GeV]");
  mom7xbbar->Draw();
  //
  c1->Print("mom7xb.pdf");

  TCanvas *c2 = new TCanvas("c2","Zbb : c",700,700);
  c2->Divide(1,2);
  //
  c2->cd(1);
  mom7xc->GetXaxis()->SetTitle("|p| [GeV]");
  mom7xc->Draw();
  //
  c2->cd(2);
  mom7xcbar->GetXaxis()->SetTitle("|p| [GeV]");
  mom7xcbar->Draw();
  //
  c2->Print("mom7xc.pdf");

  TCanvas *c3 = new TCanvas("c3","Zbb : s",700,700);
  c3->Divide(1,2);
  //
  c3->cd(1);
  mom7xs->GetXaxis()->SetTitle("|p| [GeV]");
  mom7xs->Draw();
  //
  c3->cd(2);
  mom7xsbar->GetXaxis()->SetTitle("|p| [GeV]");
  mom7xsbar->Draw();
  //
  c3->Print("mom7xs.pdf");

  TCanvas *c4 = new TCanvas("c4","Zbb : u",700,700);
  c4->Divide(1,2);
  //
  c4->cd(1);
  mom7xu->GetXaxis()->SetTitle("|p| [GeV]");
  mom7xu->Draw();
  //
  c4->cd(2);
  mom7xubar->GetXaxis()->SetTitle("|p| [GeV]");
  mom7xubar->Draw();
  //
  c4->Print("mom7xu.pdf");

  TCanvas *c5 = new TCanvas("c5","Zbb : d",700,700);
  c5->Divide(1,2);
  //
  c5->cd(1);
  mom7xd->GetXaxis()->SetTitle("|p| [GeV]");
  mom7xd->Draw();
  //
  c5->cd(2);
  mom7xdbar->GetXaxis()->SetTitle("|p| [GeV]");
  mom7xdbar->Draw();
  //
  c5->Print("mom7xd.pdf");

  TCanvas *c6 = new TCanvas("c6","Zbb : 0",700,700);
  mom7x0->GetXaxis()->SetTitle("|p| [GeV]");
  mom7x0->Draw();
  c6->Print("mom7x0.pdf");
  
  //
  
  hists->Close();
  
  return -1;
}
