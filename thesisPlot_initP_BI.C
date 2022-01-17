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


int thesisPlot_initP_BI(const string mcFile, const string dataFile){

   ROOT::RDataFrame inputFrame(pionTree,mcFile);
   ROOT::RDataFrame data_inputFrame(pionTree, dataFile);

   gInterpreter->GenerateDictionary("vector<vector<int>>", "vector");
   //gStyle->SetNdivisions(1020);
   gStyle->SetPaintTextFormat("3.2f");
   gStyle->SetOptFit(1);

   string output_name = "thesisPlot_initP_BI" + std::to_string((int) bin_size_int) + "MeV.root";

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


   THStack *stack = new THStack("stack_initP", "Beam Instrumentation Momentum; Momentum [GeV/c]; Events/20MeV");
   TH1D* hReco_MC_initP= new TH1D ("hReco_MC_initP", "MC; Momentum [GeV/c]; Events/20MeV", 50, 0.5 , 1.5);
   TH1D* hReco_58XX_initP= new TH1D ("hReco_58XX_initP", "Data 58XX Beam Instrumentation; Momentum [GeV/c]; Events/20MeV", 50, 0.5, 1.5);

   frame
      //.Filter("selected_incidentPion")
      .Foreach( [ hReco_MC_initP](double reco_initP){

               hReco_MC_initP->Fill(reco_initP);
            },{"beam_inst_P"}
            );

   dataFrame1
      //.Filter("selected_incidentPion")
      .Foreach( [ hReco_58XX_initP ](double reco_initP){

               hReco_58XX_initP->Fill(reco_initP);

            }
            ,{"beam_inst_P"}
            );

   //Normalise MC to data
   hReco_MC_initP->Scale( hReco_58XX_initP->Integral() / hReco_MC_initP->Integral() );

   hReco_MC_initP->SetLineColor(kRed);
   hReco_MC_initP->SetLineWidth(2);
   hReco_MC_initP->SetMarkerStyle(2);
   hReco_MC_initP->SetMarkerColor(kRed);

   hReco_58XX_initP->SetLineColor(kBlack);
   hReco_58XX_initP->SetMarkerStyle(50);


   stack->Add(hReco_58XX_initP);
   stack->Add(hReco_MC_initP);

   //Fit
   TF1* f1 = new TF1("f1", "gaus", 0.5, 1.5);
   hReco_58XX_initP->Fit("f1", "R 0Q");
   hReco_MC_initP->Fit("f1", "R 0Q");

   auto legend = new TLegend(0.2,0.75,0.4,0.85);
   legend->AddEntry(hReco_58XX_initP, "PDSP Data, Run 58XX");
   legend->AddEntry(hReco_MC_initP);
   legend->SetTextSize(0.03);

   gStyle->SetOptStat(11);
   gStyle->SetOptFit(0011);

   TCanvas* c_initP_mc_data= new TCanvas("canvas_initP_mc_data","canvas_initP_mc_data");
   gPad->SetGrid(1,1);
   //stack->Draw("NOSTACK PE");
   hReco_58XX_initP->SetTitle("");
   hReco_58XX_initP->GetYaxis()->SetRangeUser(0,9000);
   hReco_58XX_initP->Draw("PE");
   hReco_MC_initP->Draw("HIST SAMES");
   legend->Draw();
   //c_initP_mc_data->BuildLegend();
   c_initP_mc_data->Update();
   c_initP_mc_data->Write();

 
   return 0;
}

