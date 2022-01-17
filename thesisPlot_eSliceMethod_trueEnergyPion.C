#include "TGraph.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TGaxis.h"
#include "TLegend.h"
#include "TFile.h"
#include "../eventSelection.h"
#include "../lambda.h"
#include "../betheBloch.h"
#include "eSlice.h"

void thesisPlot_eSliceMethod_trueEnergyPion(void)
{
  // note: this macro assumes that both graphs have the same x-axis range
  //
  
  TH1D* dEdx = new TH1D("dEdx", "", 380, 0, 380 );
  hist_bethe_mean_distance( KE_in_pion, mass_pion, dEdx);
  
  TH1* cumulative = dEdx->GetCumulative();
  cumulative->GetXaxis()->SetTitle("True track length [cm]"); cumulative->GetYaxis()->SetTitle("True deposited energy [MeV]");   
  
  TCanvas *c = new TCanvas("c", "", 600,400);
  //TPad *p1 = new TPad("p1", "", 0, 0, 1, 1);
  //p1->SetGrid();
  //TPad *p2 = new TPad("p2", "", 0, 0, 1, 1);
  //p2->SetFillStyle(4000); // will be transparent
  // p2->SetGrid();
  
  //p1->Draw();
  //p1->cd();
  cumulative->SetLineColor(kBlue);
  cumulative->Draw();
  //dEdx->Draw();
  c->Update();
  
  //scale to pad coordinates
  Float_t rightmax = 1.1*dEdx->GetMaximum();
  Float_t scale = gPad->GetUymax()/rightmax;
  dEdx->SetLineColor(kRed);
  dEdx->SetMarkerColor(kRed);
  dEdx->SetMarkerSize(0.4);
  dEdx->Scale(scale);
  //dEdx->GetYaxis()->SetRangeUser(1.5,4);
  dEdx->Draw("C SAME");


  TGaxis *axis = new TGaxis(gPad->GetUxmax(), gPad->GetUymin(), gPad->GetUxmax(), gPad->GetUymax(), 0, rightmax, 510, "+L");
  axis->SetLineColor(kRed);
  axis->SetLabelColor(kRed);
  axis->SetTitleColor(kRed);
  axis->SetLabelSize(0.03);
  axis->SetTitle("dEdx [MeV / cm]");
  axis->Draw();
  
  TLegend *leg = new TLegend(0.76, 0.45, 0.86, 0.55);
  leg->SetFillColor(0);
  leg->SetTextSize(0.036);
  leg->SetHeader("1GeV Momentum Pion", "C");
  leg->AddEntry(dEdx, "dEdx", "L");
  leg->AddEntry(cumulative, "Deposited energy", "L");
  leg->Draw();
  
  c->cd();
}

