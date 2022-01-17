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


int thesisPlot_vertexMichel(const string mcFile, const string dataFile){

   ROOT::RDataFrame inputFrame(pionTree,mcFile);
   ROOT::RDataFrame data_inputFrame(pionTree, dataFile);
   gInterpreter->GenerateDictionary("vector<vector<int>>", "vector");
   gStyle->SetNdivisions(510);
   gStyle->SetPaintTextFormat("3.2f");
   gStyle->SetOptFit(1);
   gStyle->SetPalette(6, beamParticle_color);

   string output_name = "thesisPlot_vertexMichel.root";

   TFile *output = new TFile ( output_name.c_str() , "RECREATE");

   //Selected Process and RecoE Int and Inc Histos
   //THIS is what I get from DATA
   
   auto frame = inputFrame
      //Filters
      .Filter("true_beam_endZ >0")
      .Filter("primary_isBeamType && passBeamQuality_TPCjustPosition && primary_ends_inAPA3")
      //Classifications
      .Define("class_true_pion", "if(true_primPionInel == 1) return true; else return false;")
      .Define("class_true_proton", "true_beam_PDG == 2212")
      .Define("class_true_cosmic", "isCosmic")
      .Define("class_true_beamMuon", "primaryMuon")
      .Define("class_true_decayPion", "isDecay && !primaryMuon && !isCosmic")
      .Define("class_true_other", "!class_true_pion && !class_true_proton && !class_true_cosmic && !class_true_beamMuon && !class_true_decayPion")
      ;

   auto dataFrame = data_inputFrame
      .Filter("primary_isBeamType && passBeamQuality_TPCjustPosition && primary_ends_inAPA3");
      //Plot Variables
   
   TH1D* hReco_MC_vertexMichel_pion= new TH1D ("hReco_MC_vertexMichel_pion", "pion; CNN MichelScore; Events / 0.02 ", 50, 0. , 1.);
   TH1D* hReco_MC_vertexMichel_proton= new TH1D ("hReco_MC_vertexMichel_proton", "proton; CNN MichelScore; Events / 0.02 ", 50, 0. , 1.);
   TH1D* hReco_MC_vertexMichel_beamMuon= new TH1D ("hReco_MC_vertexMichel_beamMuon", "beamMuon; CNN MichelScore; Events / 0.02 ", 50, 0. , 1.);
   TH1D* hReco_MC_vertexMichel_cosmic= new TH1D ("hReco_MC_vertexMichel_cosmic", "cosmic; CNN MichelScore; Events / 0.02 ", 50, 0. , 1.);
   TH1D* hReco_MC_vertexMichel_decayPion= new TH1D ("hReco_MC_vertexMichel_decayPion", "decayPion; CNN MichelScore; Events / 0.02 ", 50, 0. , 1.);
   TH1D* hReco_MC_vertexMichel_other= new TH1D ("hReco_MC_vertexMichel_other", "other; CNN MichelScore; Events / 0.02 ", 50, 0. , 1.); 

   TH1D* hReco_58XX_vertexMichel= new TH1D ("hReco_58XX_vertexMichel", "Data; CNN MichelScore; Events / 0.02 ", 50, 0. , 1.);
 

   frame
      .Foreach( [ hReco_MC_vertexMichel_pion, hReco_MC_vertexMichel_proton, hReco_MC_vertexMichel_beamMuon,
                  hReco_MC_vertexMichel_cosmic, hReco_MC_vertexMichel_decayPion, hReco_MC_vertexMichel_other]
                
            (double vertex_cnn, int vertex_hits,
             bool pion, bool proton, bool beamMuon, bool cosmic, bool decayPion, bool other){

                  double score = vertex_cnn / vertex_hits;

                  if(pion) {
                  hReco_MC_vertexMichel_pion->Fill(score);}
                  else if(proton) {
                  hReco_MC_vertexMichel_proton->Fill(score);} 
                  else if(beamMuon) {
                  hReco_MC_vertexMichel_beamMuon->Fill(score);} 
                  else if(cosmic) {
                  hReco_MC_vertexMichel_cosmic->Fill(score);} 
                  else if(decayPion) {
                  hReco_MC_vertexMichel_decayPion->Fill(score); } 
                  else if(other) {
                  hReco_MC_vertexMichel_other->Fill(score); }

               },{"reco_beam_vertex_michel_score", "reco_beam_vertex_nHits", 
                  "class_true_pion", "class_true_proton", "class_true_beamMuon", "class_true_cosmic", "class_true_decayPion", "class_true_other"}
            );

   dataFrame
      .Foreach( [ hReco_58XX_vertexMichel ](double vertex_cnn, int vertex_nHits ){

               double score = vertex_cnn / vertex_nHits;
               hReco_58XX_vertexMichel->Fill( score );
               },{"reco_beam_vertex_michel_score", "reco_beam_vertex_nHits"}
            );

   //MC Total Integral
   double totalMC =   hReco_MC_vertexMichel_pion->Integral() + hReco_MC_vertexMichel_proton->Integral()
                      + hReco_MC_vertexMichel_beamMuon->Integral() + hReco_MC_vertexMichel_cosmic->Integral()
                      + hReco_MC_vertexMichel_decayPion->Integral() + hReco_MC_vertexMichel_other->Integral();

   double scale = hReco_58XX_vertexMichel->Integral() / totalMC;

   //Normalise MC to data 
   hReco_MC_vertexMichel_pion->Scale( scale ); hReco_MC_vertexMichel_proton->Scale( scale ); hReco_MC_vertexMichel_beamMuon->Scale( scale );
   hReco_MC_vertexMichel_cosmic->Scale( scale ); hReco_MC_vertexMichel_decayPion->Scale( scale ); hReco_MC_vertexMichel_other->Scale( scale );

   //THStack
   THStack *stack_MC_vertexMichel = new THStack("stack_MC_vertexMichel", "; CNN MichelScore; Events / 0.02");
   stack_MC_vertexMichel->Add( hReco_MC_vertexMichel_pion ); stack_MC_vertexMichel->Add( hReco_MC_vertexMichel_proton ); stack_MC_vertexMichel->Add( hReco_MC_vertexMichel_beamMuon );
   stack_MC_vertexMichel->Add( hReco_MC_vertexMichel_cosmic ); stack_MC_vertexMichel->Add( hReco_MC_vertexMichel_decayPion ); stack_MC_vertexMichel->Add( hReco_MC_vertexMichel_other );

   hReco_58XX_vertexMichel->SetLineColor(kBlack);
   hReco_58XX_vertexMichel->SetMarkerStyle(8);
   hReco_58XX_vertexMichel->SetMarkerSize(0.8);
   
   auto legend = new TLegend(0.5,0.63,0.8,0.87);
   legend->AddEntry(hReco_58XX_vertexMichel, "PDSP Data, Run 58XX");
   legend->AddEntry(hReco_MC_vertexMichel_pion, "Pion"); legend->AddEntry(hReco_MC_vertexMichel_proton, "Proton"); legend->AddEntry(hReco_MC_vertexMichel_beamMuon, "Beam Muon"); 
   legend->AddEntry(hReco_MC_vertexMichel_cosmic, "Cosmic"); legend->AddEntry(hReco_MC_vertexMichel_decayPion, "Decay Pion"); legend->AddEntry(hReco_MC_vertexMichel_other, "Other"); 
   legend->SetTextSize(0.03);

   gStyle->SetOptStat(0);
   gStyle->SetOptFit(0);

   TCanvas* c_vertexMichel_mc_data= new TCanvas("c_vertexMichel_mc_data","c_vertexMichel_mc_data", 1000, 900);
   gPad->SetGrid(1,1);
   gPad->SetLogy();
   stack_MC_vertexMichel->Draw("HIST PFC PMC PLC");
   hReco_58XX_vertexMichel->Draw("PE SAME");
   legend->Draw();
   c_vertexMichel_mc_data->Update();
   c_vertexMichel_mc_data->Write();
 
   return 0;
}

