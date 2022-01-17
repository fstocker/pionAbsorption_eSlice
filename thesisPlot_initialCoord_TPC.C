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


int thesisPlot_initialCoord_TPC(const string mcFile, const string dataFile){

   ROOT::RDataFrame inputFrame(pionTree,mcFile);
   ROOT::RDataFrame data_inputFrame(pionTree, dataFile);

   gInterpreter->GenerateDictionary("vector<vector<int>>", "vector");
   gStyle->SetNdivisions(1020);
   gStyle->SetPaintTextFormat("3.2f");
   gStyle->SetOptFit(1);

   string output_name = "thesisPlot_initP_BI" + std::to_string((int) bin_size_int) + "MeV.root";

   TFile *output = new TFile ( output_name.c_str() , "RECREATE");

   //Selected Process and RecoE Int and Inc Histos
   //THIS is what I get from DATA
   
   auto frame = inputFrame
      //.Filter("true_beam_endZ >0 && reco_initKE <1200")
      .Filter("true_beam_endZ >0")
      .Filter("primary_isBeamType");

   auto dataFrame1 = data_inputFrame1
      .Filter("primary_isBeamType");


   TH1D* hReco_MC_startX= new TH1D ("hReco_MC_startX", "MC; Reconstructed start X[cm]; % / 2 cm", 50, -80. , 20.);
   TH1D* hReco_58XX_startX= new TH1D ("hReco_58XX_startX", "Data; Reconstructed start X[cm]; % / 2 cm", 50, -80. , 20.);
   
   TH1D* hReco_MC_startY= new TH1D ("hReco_MC_startY", "MC; Reconstructed start Y[cm]; % / 2 cm", 70 , 360. , 500.);
   TH1D* hReco_58XX_startY= new TH1D ("hReco_58XX_startY", "Data; Reconstructed startY[cm]; % / 2 cm", 70, 360. , 500.);

   TH1D* hReco_MC_startZ= new TH1D ("hReco_MC_startZ", "MC; Reconstructed start Z[cm]; % / 0.2 cm", 75, -5. , 10.);
   TH1D* hReco_58XX_startZ= new TH1D ("hReco_58XX_startZ", "Data; Reconstructed start Z[cm]; % / 0.2 cm", 75, -5. , 10.);

   frame
      .Foreach( [ hReco_MC_startX, hReco_MC_startY, hReco_MC_startZ](double startX, double startY, double startZ){

               hReco_MC_startX->Fill(startX);
               hReco_MC_startY->Fill(startY);
               hReco_MC_startZ->Fill(startZ);
            
               },{"reco_beam_calo_startX", "reco_beam_calo_startY", "reco_beam_calo_startZ"}
            );

   dataFrame1
      .Foreach( [ hReco_58XX_startX, hReco_58XX_startY, hReco_58XX_startZ](double startX, double startY, double startZ){

               hReco_58XX_startX->Fill(startX);
               hReco_58XX_startY->Fill(startY);
               hReco_58XX_startZ->Fill(startZ);
            
               },{"reco_beam_calo_startX", "reco_beam_calo_startY", "reco_beam_calo_startZ"}
            );

    //Normalise MC to data
   hReco_MC_startX->Scale( 1 / hReco_MC_startX->GetEntries() );
   hReco_MC_startY->Scale( 1 / hReco_MC_startY->GetEntries() );
   hReco_MC_startZ->Scale( 1 / hReco_MC_startZ->GetEntries() );

   hReco_58XX_startX->Scale( 1 / hReco_58XX_startX->GetEntries() );
   hReco_58XX_startY->Scale( 1 / hReco_58XX_startY->GetEntries() );
   hReco_58XX_startZ->Scale( 1 / hReco_58XX_startZ->GetEntries() );

   //Cosmetics
   hReco_MC_startX->SetLineColor(kRed);
   hReco_MC_startX->SetLineWidth(2);
   hReco_MC_startX->SetMarkerStyle(8);
   hReco_MC_startX->SetMarkerColor(kRed);
   hReco_MC_startX->SetMarkerSize(0.6);

   hReco_58XX_startX->SetLineColor(kBlack);
   hReco_58XX_startX->SetMarkerStyle(8);
   hReco_58XX_startX->SetMarkerSize(0.6);

   hReco_MC_startY->SetLineColor(kRed);
   hReco_MC_startY->SetLineWidth(2);
   hReco_MC_startY->SetMarkerStyle(8);
   hReco_MC_startY->SetMarkerColor(kRed);
   hReco_MC_startY->SetMarkerSize(0.6);

   hReco_58XX_startY->SetLineColor(kBlack);
   hReco_58XX_startY->SetMarkerStyle(8);
   hReco_58XX_startY->SetMarkerSize(0.6);

   hReco_MC_startZ->SetLineColor(kRed);
   hReco_MC_startZ->SetLineWidth(2);
   hReco_MC_startZ->SetMarkerStyle(8);
   hReco_MC_startZ->SetMarkerColor(kRed);
   hReco_MC_startZ->SetMarkerSize(0.6);

   hReco_58XX_startZ->SetLineColor(kBlack);
   hReco_58XX_startZ->SetMarkerStyle(8);
   hReco_58XX_startZ->SetMarkerSize(0.6);


   //Fit
   TF1* f_mc_x = new TF1("f_mc_x", "gaus", -60, 0);
   f_mc_x->SetLineColor(kRed);
   f_mc_x->SetLineWidth(1);
   hReco_MC_startX->Fit("f_mc_x", "R Q");

   TF1* f_data_x = new TF1("f_data_x", "gaus", -60, 0);
   f_data_x->SetLineWidth(1);
   hReco_58XX_startX->Fit("f_data_x", "R Q");

   TF1* f_mc_y = new TF1("f_mc_y", "gaus", 400, 450);
   f_mc_y->SetLineColor(kRed);
   f_mc_y->SetLineWidth(1);
   hReco_MC_startY->Fit("f_mc_y", "R Q");

   TF1* f_data_y = new TF1("f_data_y", "gaus", 400, 450);
   f_data_y->SetLineWidth(1);
   hReco_58XX_startY->Fit("f_data_y", "R Q");

   TF1* f_mc_z = new TF1("f_mc_z", "gaus", -2, 3);
   f_mc_z->SetLineColor(kRed);
   f_mc_z->SetLineWidth(1);
   hReco_MC_startZ->Fit("f_mc_z", "R Q");

   TF1* f_data_z = new TF1("f_data_z", "gaus", -2, 7);
   f_data_z->SetLineWidth(1);
   hReco_58XX_startZ->Fit("f_data_z", "R Q");

   auto legend_x = new TLegend(0.2,0.75,0.45,0.85);
   legend_x->AddEntry(hReco_58XX_startX, "PDSP Data, Run 58XX");
   legend_x->AddEntry(hReco_MC_startX);
   legend_x->SetTextSize(0.03);

   auto legend_y = new TLegend(0.2,0.75,0.45,0.85);
   legend_y->AddEntry(hReco_58XX_startY, "PDSP Data, Run 58XX");
   legend_y->AddEntry(hReco_MC_startY);
   legend_y->SetTextSize(0.03);
   
   auto legend_z = new TLegend(0.6,0.25,0.83,0.35);
   legend_z->AddEntry(hReco_58XX_startZ, "PDSP Data, Run 58XX");
   legend_z->AddEntry(hReco_MC_startZ);
   legend_z->SetTextSize(0.03);

   gStyle->SetOptStat(11);
   gStyle->SetOptFit(0011);

   TCanvas* c_startX_mc_data= new TCanvas("c_startX_mc_data","c_startX_mc_data");
   gPad->SetGrid(1,1);
   hReco_58XX_startX->SetTitle("");
   hReco_58XX_startX->GetYaxis()->SetRangeUser(0,0.2);
   hReco_58XX_startX->Draw("PE1");
   hReco_MC_startX->Draw("PE1 SAMES");
   legend_x->Draw();
   c_startX_mc_data->Update();
   c_startX_mc_data->Write();

   TCanvas* c_startY_mc_data= new TCanvas("c_startY_mc_data","c_startY_mc_data");
   gPad->SetGrid(1,1);
   hReco_58XX_startY->SetTitle("");
   hReco_58XX_startY->GetYaxis()->SetRangeUser(0,0.2);
   hReco_58XX_startY->Draw("PE1");
   hReco_MC_startY->Draw("PE1 SAMES");
   legend_y->Draw();
   c_startY_mc_data->Update();
   c_startY_mc_data->Write();

   TCanvas* c_startZ_mc_data= new TCanvas("c_startZ_mc_data","c_startZ_mc_data");
   gPad->SetGrid(1,1);
   hReco_58XX_startZ->SetTitle("");
   hReco_58XX_startZ->GetYaxis()->SetRangeUser(0, 0.4);
   hReco_58XX_startZ->Draw("PE1");
   hReco_MC_startZ->Draw("PE1 SAMES");
   legend_z->Draw();
   c_startZ_mc_data->Update();
   c_startZ_mc_data->Write();

 
   return 0;
}

