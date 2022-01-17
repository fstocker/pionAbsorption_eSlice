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
#include "TMatrixDBase.h"
#include "TArray.h"
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


int thesisPlot_initKE_reweight(const string mcFile, const string dataFile){

   ROOT::RDataFrame inputFrame(pionTree,mcFile);
   ROOT::RDataFrame data_inputFrame(pionTree, dataFile);


   gInterpreter->GenerateDictionary("vector<vector<int>>", "vector");
   //gStyle->SetNdivisions(1020);
   gStyle->SetPaintTextFormat("3.2f");
   //gStyle->SetOptFit(0);
   //gStyle->SetOptStat(0);

   string output_name = "thesisPlot_initKE_reweight" + std::to_string((int) bin_size_int) + "MeV.root";

   TFile *output = new TFile ( output_name.c_str() , "RECREATE");

   //Selected Process and RecoE Int and Inc Histos
   //THIS is what I get from DATA
   
   auto frame = inputFrame
      //.Filter("true_beam_endZ >0 && reco_initKE <1200")
      .Filter("true_beam_endZ >0")
      //.Range(53000)
      .Define("true_equalBin", equalBin, {"true_initKE", "true_interKE"})
      .Define("reco_equalBin", equalBin, {"reco_initKE", "reco_interKE"});

   auto dataFrame1 = data_inputFrame1
      .Define("reco_equalBin", equalBin, {"reco_initKE", "reco_interKE"});


   TH1D* hReco_mc_initKE= new TH1D ("hReco_mc_initKE", "MC; Kinetic Energy [MeV]; Events/ 50 MeV", nBin_int, eEnd , eStart);
   TH1D* hReco_mc_initKE_rwData= new TH1D ("hReco_mc_initKE_rwData", "MC; Kinetic Energy [MeV]; Events/ 50 MeV", nBin_int, eEnd , eStart);
   TH1D* hReco_58XX_initKE= new TH1D ("hReco_58XX_initKE", "Data 58XX Beam Instrumentation; Kinetic Energy [MeV]; Events/ 50 MeV", nBin_int, eEnd , eStart);

   frame
      //.Filter("selected_incidentPion")
      .Foreach( [ hReco_mc_initKE, hReco_mc_initKE_rwData](double initKE, double initKE_rwData){

               hReco_mc_initKE->Fill(initKE);
               hReco_mc_initKE_rwData->Fill(initKE_rwData);
            },{"reco_initKE", "reco_initKE_rwData"}
            );

   dataFrame1
      //.Filter("selected_incidentPion")
      .Foreach( [ hReco_58XX_initKE ](double reco_initP){

               hReco_58XX_initKE->Fill(reco_initP);

            }
            ,{"reco_initKE"}
            );

   //Normalise MC to data
   hReco_mc_initKE->Scale( hReco_58XX_initKE->Integral() / hReco_mc_initKE->Integral() );

   hReco_mc_initKE->SetLineColor(kRed);
   hReco_mc_initKE->SetLineWidth(2);
   hReco_mc_initKE->SetMarkerColorAlpha(kRed,0);

   hReco_mc_initKE_rwData->Scale( hReco_58XX_initKE->Integral() / hReco_mc_initKE_rwData->Integral() );

   hReco_mc_initKE_rwData->SetLineColor(kBlue);
   hReco_mc_initKE_rwData->SetLineWidth(2);
   hReco_mc_initKE_rwData->SetMarkerColorAlpha(kBlue,0);

   hReco_58XX_initKE->SetLineColor(kBlack);
   //hReco_58XX_initKE->SetMarkerStyle(50);



   //Fit
   TF1* f1 = new TF1("f1", "gaus", 450, 1100);
   TF1* f2 = new TF1("f2", "gaus", 450, 1100);
   hReco_58XX_initKE->Fit("f1", "R Q");
   hReco_mc_initKE->Fit("f2", "R Q");

   auto legend = new TLegend(0.2,0.75,0.4,0.85);
   legend->AddEntry(hReco_58XX_initKE, "PDSP Data, Run 58XX");
   legend->AddEntry(hReco_mc_initKE, "Beam instrumentation");
   legend->SetTextSize(0.03);

   gStyle->SetOptStat(11);
   gStyle->SetOptFit(0011);

   TCanvas* c_initP_mc_data= new TCanvas("canvas_initP_mc_data","canvas_initP_mc_data", 1200,800);
   gPad->SetGrid(1,1);
   //stack->Draw("NOSTACK PE");
   hReco_58XX_initKE->SetTitle("");
   //hReco_58XX_initKE->GetYaxis()->SetRangeUser(0,9000);
   hReco_58XX_initKE->GetYaxis()->SetRangeUser(0,22000);
   hReco_58XX_initKE->Draw("PE1");
   hReco_mc_initKE->Draw("HIST SAMES");
   legend->Draw();
   //c_initP_mc_data->BuildLegend();
   c_initP_mc_data->Update();
   c_initP_mc_data->Write();

   auto legend_2 = new TLegend(0.2,0.75,0.4,0.85);
   legend_2->AddEntry(hReco_58XX_initKE, "PDSP Data, Run 58XX");
   legend_2->AddEntry(hReco_mc_initKE_rwData, "Reweight beam instrumentation");
   legend_2->AddEntry(hReco_mc_initKE, "Beam instrumentation");
   legend_2->SetTextSize(0.03);

   //gStyle->SetOptStat(11);
   //gStyle->SetOptFit(0011);

   TCanvas* c_rwData= new TCanvas("c_rwData","c_rwData", 1200,800);
   gPad->SetGrid(1,1);
   //stack->Draw("NOSTACK PE");
   hReco_mc_initKE->SetTitle("");
   //hReco_58XX_initKE->GetYaxis()->SetRangeUser(0,9000);
   hReco_mc_initKE->Draw("HIST");
   hReco_mc_initKE_rwData->Draw("HIST SAME");
   hReco_58XX_initKE->Draw("PE1 SAME");
   legend_2->Draw();
   //c_rwData->Buildlegend_2();
   c_rwData->Update();
   c_rwData->Write();

   return 0;
}

