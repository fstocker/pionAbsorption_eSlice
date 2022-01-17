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
#include "../eventSelection.h"
#include "eSlice.h"
#include "plotThesis_config.h"
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

int thesisPlot_interE_vs_endZ(const string mcFile, const string dataFile){

   ROOT::RDataFrame inputFrame(pionTree,mcFile);
   ROOT::RDataFrame data_inputFrame(pionTree, dataFile);

   gInterpreter->GenerateDictionary("vector<vector<int>>", "vector");
   gStyle->SetNdivisions(1020);
   gStyle->SetPaintTextFormat("3.2f");
   gStyle->SetOptFit(1);
   //gStyle->SetPalette(6, beamParticle_color);

   string output_name = "thesisPlot_APA3.root";

   TFile *output = new TFile ( output_name.c_str() , "RECREATE");

   //Selected Process and RecoE Int and Inc Histos
   //THIS is what I get from DATA
   
   auto frame = inputFrame
      //Filters
      .Filter("true_beam_endZ >0")
      .Define("reco_equalBin", equalBin, {"reco_initKE", "reco_interKE"})
      .Filter("primary_isBeamType && passBeamQuality_TPCjustPosition && !reco_equalBin");

   auto dataFrame = data_inputFrame
      .Define("reco_equalBin", equalBin, {"reco_initKE", "reco_interKE"})
      .Filter("primary_isBeamType && passBeamQuality_TPCjustPosition && !reco_equalBin");
      //Plot Variables
   
   TH2D* hReco_MC_interE_endZ= new TH2D ("hReco_MC_interE_endZ", ";  Reconstructed end Z [cm]; Reconstructed interacting kinetic Energy [MeV]", 16, 0. , 400., 24, 0., 1200.);
   TH2D* hReco_data_interE_endZ= new TH2D ("hReco_data_interE_endZ", ";  Reconstructed end Z [cm]; Reconstructed interacting kinetic Energy [MeV]", 16, 0. , 400., 24, 0., 1200.);

   frame
      .Foreach( [ hReco_MC_interE_endZ ](double endZ, double interE){

                  hReco_MC_interE_endZ->Fill( endZ, interE );
 
               },{"reco_beam_endZ", "reco_interKE_rwData"}
            );

   dataFrame
      .Foreach( [ hReco_data_interE_endZ ](double endZ, double interE ){

                  hReco_data_interE_endZ->Fill( endZ, interE );

               },{"reco_beam_endZ", "reco_interKE"}
            );

   //Normalise histos
   //
   double mc_entries = hReco_MC_interE_endZ->GetEntries();
   hReco_MC_interE_endZ->Scale( 1 / mc_entries );
   double data_entries = hReco_data_interE_endZ->GetEntries();
   hReco_data_interE_endZ->Scale( 1 / data_entries);

   gStyle->SetOptStat(0);
   gStyle->SetOptFit(0);

   TCanvas* c_mc_endZ_interE= new TCanvas("c_mc_endZ_interE","c_mc_endZ_interE", 1200, 900);
   gPad->SetGrid(1,1);
   hReco_MC_interE_endZ->Draw("colz");
   c_mc_endZ_interE->Update();
   c_mc_endZ_interE->Write();

   TCanvas* c_data_endZ_interE= new TCanvas("c_data_endZ_interE","c_data_endZ_interE", 1200, 900);
   gPad->SetGrid(1,1);
   hReco_data_interE_endZ->Draw("colz");
   c_data_endZ_interE->Update();
   c_data_endZ_interE->Write();
 
   return 0;
}

