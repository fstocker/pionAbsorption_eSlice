#include "TCanvas.h"
#include "TROOT.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TH1.h"
#include "TF1.h"
#include "THStack.h"
#include "TLegend.h"
#include "TArrow.h"
#include "TStyle.h"
#include "TColor.h"
#include "TLatex.h"
#include "TMath.h"
#include "../lambda.h"
#include "../betheBloch.h"
#include "eSlice.h"
#include <ROOT/RDataFrame.hxx>


#include <iostream>
#include <math.h>
#include <string.h>
#include <stdio.h>

#include <vector>

//using RDataFrame to cut and analyse PionTtrr

using namespace std;
using namespace ROOT::VecOps;

//Applying the eSlice Method

//***********************
//Main Function

int thesisPlot_eSliceMethod_trueClosureTest(const string mcFile){

   gInterpreter->GenerateDictionary("vector<vector<int>>", "vector");
   ROOT::RDataFrame frame(pionTree, mcFile);
   gStyle->SetNdivisions(1020);
   gStyle->SetOptStat(0);

   //output file
   //string outputNameMC = "output_eSliceMethod_trueProcess_trueE.root";
   //string outputName;
   //outputName = outputNameMC;


//   TFile *output = new TFile ( outputName.c_str() , "RECREATE");

   //access Jakes GeantFile in folder
   TFile f1("exclusive_xsec.root");
   TGraph *totInel_KE = (TGraph*)f1.Get("total_inel_KE");
   TGraph *abs_KE = (TGraph*)f1.Get("abs_KE");
   TGraph *cex_KE = (TGraph*)f1.Get("cex_KE");
   f1.Close();

   string mg_title;
   mg_title = ";True kinetic energy [MeV]; #sigma [mbarn]";

   TMultiGraph *mg = new TMultiGraph();
   mg->Add(totInel_KE);
   mg->Add(abs_KE);
   //mg->Add(cex_KE);
   mg->SetTitle(mg_title.c_str());

   //switch to output-file
   //output->cd();
   double eStart_true = 1500;
   double eEnd_true = 0;
   int nBin_true = 30;

   TH1D* hTrue_initE= new TH1D ("trueMC_initE", "MC Truth initE; True kinetic energy [MeV]; Events / 50 MeV", nBin_true, eEnd_true, eStart_true);
   TH1D* hTrue_interE= new TH1D ("trueMC_interE", "MC Truth interE; True kinetic energy [MeV]; Events / 50 MeV", nBin_true, eEnd_true, eStart_true);
   TH1D* hTrue_totInel= new TH1D ("trueMC_totInel", "MC Truth totInel; True kinetic energy [MeV]; Events / 50 MeV", nBin_true, eEnd_true, eStart_true);
   TH1D* hTrue_abs= new TH1D ("trueMC_abs", "MC Truth abs; True kinetic energy [MeV]; Events / 50 MeV", nBin_true, eEnd_true, eStart_true);
   TH1D* hTrue_incident= new TH1D ("trueMC_incident", "MC Truth incident; True kinetic energy [MeV]; Events / 50 MeV", nBin_true, eEnd_true, eStart_true);

   frame
      .Define("true_equalBin", equalBin, {"true_initKE", "true_interKE"})
      .Filter("true_beam_endZ > 0 && !true_equalBin")
      .Filter("true_beam_PDG == 211")
      .Foreach( [hTrue_initE, hTrue_interE ]( double true_initE, double true_interE)
            {

            hTrue_initE->Fill( true_initE );
            hTrue_interE->Fill( true_interE );

            },
      {"true_initKE", "true_interKE"}
   );

   frame
      .Define("true_equalBin", equalBin, {"true_initKE", "true_interKE"})
      .Filter("true_beam_endZ > 0 && !true_equalBin")
      .Filter("true_primPionInel && !isDecay")
      .Foreach( [hTrue_totInel ](double true_interE)
            {

             hTrue_totInel->Fill( true_interE );

             },
      { "true_interKE"}
   );
   frame
      .Define("true_equalBin", equalBin, {"true_initKE", "true_interKE"})
      .Filter("true_beam_endZ > 0 && !true_equalBin")
      .Filter("true_absSignal && true_pion_daughter == 0")
      .Foreach( [hTrue_abs ](double true_interE)
            {

             hTrue_abs->Fill( true_interE );

             },
      { "true_interKE"}
   );

   build_incidentHist( hTrue_initE, hTrue_interE, hTrue_incident );

   //=====================================================
   //            Prepare BetheBloch Mean for each Bin 
   //            QUESTION: take betheBloch of Pion or Muon?? Comparison to data fits better muon Bethe... 
   //            at hihger momentum ~400-800 they anway are almost the same
   //=====================================================
   TH1D* h_betheMean_pion = new TH1D("h_betheMean_pion", "Mean Energy Loss", nBin_true, eEnd_true, eStart_true);

   //fill histo with Mean dEdX of bin center
   fill_betheHisto( h_betheMean_pion, mass_pion );
   //h_betheMean_muon->Write();


   //=====================================================
   //             Computing the XS
   //=====================================================
   //
   // xs(Ebin) = (A / (Na*density*bin_size)) * dEdX(Ebin) * hInteracting / hIncident
   //
   // More Accurate use log( Ninc / (Ninc - Nint ))
   //
   //scale_factor is done in eSlice.h

   //------------------------------------------------------
   //    Absorption, Selected Interactions Reconstrucetd Energy
   //------------------------------------------------------


   TH1D* h_xs_trueE_trueAbs = new TH1D("h_xs_trueE_trueAbs", "Absorption MC", nBin_true, eEnd_true, eStart_true);
   TH1D* h_xs_trueE_trueTotInel = new TH1D("h_xs_trueE_trueTotInel", "Total Inelastic MC", nBin_true, eEnd_true, eStart_true);

   //Function to do XS with log formula

   do_XS_log( h_xs_trueE_trueAbs, hTrue_abs, hTrue_incident, h_betheMean_pion);
   
   do_XS_log_binomial_error( h_xs_trueE_trueAbs, hTrue_abs, hTrue_incident, h_betheMean_pion);

   do_XS_log( h_xs_trueE_trueTotInel, hTrue_totInel, hTrue_incident, h_betheMean_pion);
   
   do_XS_log_binomial_error( h_xs_trueE_trueTotInel, hTrue_totInel, hTrue_incident, h_betheMean_pion);   



   //=====================================================*
   //            Plotting and Style
   //=====================================================
   //
   TCanvas *c_trueE_abs = new TCanvas("c_trueE_abs", "c_trueE_abs");
   gPad->SetGrid(1,1);
   h_xs_trueE_trueAbs->SetTitle( "Absorption");
   h_xs_trueE_trueTotInel->SetTitle( "Total inelastic");
   //h_xs_trueE_trueAbs->GetXaxis()->SetRangeUser(100,1000);
   h_xs_trueE_trueAbs->GetXaxis()->SetNdivisions(1020);
   h_xs_trueE_trueAbs->GetYaxis()->SetNdivisions(1020);

   totInel_KE->SetTitle( "GEANT4 pion total inelastic;True kinetic energy [MeV]; #sigma [mb]");
   //totInel_KE->GetXaxis()->SetRangeUser(eEnd, eStart);
   totInel_KE->SetLineColor(kRed);
   totInel_KE->SetLineWidth(3);
   //totInel_KE->Draw("AC");
   abs_KE->SetTitle( "GEANT4 pion absorption;True kinetic energy [MeV]; #sigma [mb]");
   //abs_KE->GetXaxis()->SetRangeUser(eEnd, eStart);
   abs_KE->SetLineColor(kBlue);
   abs_KE->SetLineWidth(3);
   mg->GetXaxis()->SetRangeUser(0,1000);
   mg->GetYaxis()->SetRangeUser(0,1200);
   mg->Draw("AC SAME");
   h_xs_trueE_trueTotInel->SetMarkerSize(0.7);
   h_xs_trueE_trueTotInel->Draw("PE0 SAME");
   h_xs_trueE_trueAbs->SetMarkerSize(0.7);
   h_xs_trueE_trueAbs->Draw("PE0 SAME");
   c_trueE_abs->BuildLegend();

   //c_trueE_abs->Write();
   //
   TCanvas *c_initE = new TCanvas("c_initE", "c_initE");
   gPad->SetGrid(1,1);
   hTrue_initE->GetXaxis()->SetRangeUser(0,1000);
   hTrue_initE->SetLineColor(kGreen + 2);
   hTrue_initE->SetLineWidth(2);
   hTrue_initE->SetFillColorAlpha(kGreen + 2, 0.2);
   hTrue_initE->SetTitle("Pion initial distribution");
   hTrue_initE->Draw("HIST");
   
   TCanvas *c_interE = new TCanvas("c_interE", "c_interE");
   gPad->SetGrid(1,1);
   hTrue_interE->GetXaxis()->SetRangeUser(0,1000);
   hTrue_interE->SetLineColor(kGreen + 2);
   hTrue_interE->SetLineWidth(2);
   hTrue_interE->SetFillColorAlpha(kGreen + 2, 0.2);
   hTrue_interE->SetTitle("Pion interacting distribution");
   hTrue_interE->Draw("HIST");
   
   TCanvas *c_abs = new TCanvas("c_abs", "c_abs");
   gPad->SetGrid(1,1);
   hTrue_abs->GetXaxis()->SetRangeUser(0,1000);
   hTrue_abs->SetLineColor(kBlue);
   hTrue_abs->SetLineWidth(2);
   hTrue_abs->SetFillColorAlpha(kGreen + 2, 0.2);
   hTrue_abs->SetTitle("Pion absorption interacting distribution");
   hTrue_abs->Draw("HIST");
    
   TCanvas *c_incident = new TCanvas("c_incident", "c_incident");
   gPad->SetGrid(1,1);
   hTrue_incident->GetXaxis()->SetRangeUser(0,1000);
   hTrue_incident->SetLineColor(kGreen + 2);
   hTrue_incident->SetLineWidth(2);
   hTrue_incident->SetFillColorAlpha(kGreen + 2, 0.2);
   hTrue_incident->SetTitle("Pion incident distribution");
   hTrue_incident->Draw("HIST");

  // TCanvas *c_trueE_totInel = new TCanvas("c_trueE_totInel", "c_trueE_totInel");
  // gPad->SetGrid(1,1);
  // h_xs_trueE_truePion_totInel->SetTitle( "True Total Inelastic; true Kinetic Energy (MeV); #sigma (mb)");
  // h_xs_trueE_truePion_totInel->GetXaxis()->SetRangeUser(100,1100);
  // h_xs_trueE_truePion_totInel->GetXaxis()->SetNdivisions(1020);
  // h_xs_trueE_truePion_totInel->GetYaxis()->SetNdivisions(1020);

  // totInel_KE->SetTitle("Total Inelastic MC; true kinetic Energy (MeV); #sigma (mbarn)");
  // totInel_KE->GetXaxis()->SetRangeUser(eEnd, eStart);
  // totInel_KE->SetLineColor(kRed);
  // totInel_KE->SetLineWidth(3);
  // totInel_KE->Draw("AC");
  // h_xs_trueE_truePion_totInel->SetMarkerSize(0.7);
  // h_xs_trueE_truePion_totInel->Draw("PE0 SAME");

  // c_trueE_totInel->Write();
  // //output->Write();
  // //f1.Close();
  // TCanvas *c_all = new TCanvas("c_all", "c_all");
  // gPad->SetGrid(1,1);

  // //cex_KE->SetLineColor(kGreen);
  // abs_KE->SetLineColor(kBlue);
  // mg->GetXaxis()->SetRangeUser(0,1000);
  // mg->Draw("AC");
  // h_xs_trueE_truePion_totInel->Draw("PE0 SAME");
  // h_xs_trueE_trueAbs->Draw("PE0 SAME");

  // c_all->Write();
   return 0;
}


