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


int thesisPlot_vertexRes(const string mcFile){

   ROOT::RDataFrame inputFrame(pionTree, mcFile);

   gInterpreter->GenerateDictionary("vector<vector<int>>", "vector");
   gStyle->SetNdivisions(1020);
   gStyle->SetPaintTextFormat("3.2f");
   gStyle->SetOptFit(1);

   string output_name = "thesisPlot_vertexRes.root";

   TFile *output = new TFile ( output_name.c_str() , "RECREATE");

   //Selected Process and RecoE Int and Inc Histos
   //THIS is what I get from DATA
   
   auto frame = inputFrame
      //Filters
      .Filter("true_beam_endZ >0")
      .Filter("true_beam_PDG == 211 && true_primPionInel == 1 && !isDecay")
      .Filter("primary_isBeamType && reco_beam_endZ < 220")
      .Filter("reco_beam_endX != -999 && reco_beam_endY != -999 && reco_beam_endZ != -999")
      .Define("reco_true_intDist", "sqrt( pow(reco_beam_endX - true_beam_endX_SCE ,2) + pow(reco_beam_endY - true_beam_endY_SCE ,2) + pow(reco_beam_endZ - true_beam_endZ_SCE ,2) )")
      .Define("reco_true_zRes", "reco_beam_endZ - true_beam_endZ_SCE");

   
   TH1D* hReco_vertexRes = new TH1D ("hReco_vertexRes", "; Vertex #DeltaR (Reco - True) [cm]; Events / 2 cm", 25, -25. , 25.);
   TH1D* hReco_zRes = new TH1D ("hDiff_reco_true_endZ", "; Reco - True #DeltaZ [cm]; Events / 1 cm", 100 , -50. , 50.);

 

   frame
      .Foreach( [ hReco_vertexRes, hReco_zRes ]
                
            (double dist, double zRes){
               
            if(zRes < 0) hReco_vertexRes->Fill(-dist);
            else hReco_vertexRes->Fill(dist);
               
            hReco_zRes->Fill(zRes);

               },{"reco_true_intDist", "reco_true_zRes"}
            );

   
   //auto legend_z = new TLegend(0.5,0.65,0.76,0.85);
   //legend_z->AddEntry(hReco_58XX_endZ, "PDSP Data, Run 58XX");
   //legend_z->AddEntry(hReco_MC_endZ_pion, "Pion"); legend_z->AddEntry(hReco_MC_endZ_proton, "Proton"); legend_z->AddEntry(hReco_MC_endZ_beamMuon, "Beam Muon"); 
   //legend_z->AddEntry(hReco_MC_endZ_cosmic, "Cosmic"); legend_z->AddEntry(hReco_MC_endZ_decayPion, "Decay Pion"); legend_z->AddEntry(hReco_MC_endZ_other, "Other"); 
   //legend_z->SetTextSize(0.03);

   //gStyle->SetOptStat(0);
   //gStyle->SetOptFit(0);

   TCanvas* c_resVertex= new TCanvas("c_resVertex","c_resVertex", 1200, 800);
   gPad->SetGrid(1,1);
   hReco_vertexRes->Draw("HIST");
   //legend_z->Draw();
   c_resVertex->Update();
   c_resVertex->Write();
 
   TCanvas* c_resZ= new TCanvas("c_resZ","c_resZ", 1200, 800);
   gPad->SetGrid(1,1);
   hReco_zRes->Draw("HIST");
   //legend_z->Draw();
   c_resZ->Update();
   c_resZ->Write();

   return 0;
}

