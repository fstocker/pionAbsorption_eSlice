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


int thesisPlot_APA3(const string mcFile, const string dataFile){

   ROOT::RDataFrame inputFrame(pionTree,mcFile);
   ROOT::RDataFrame data_inputFrame(pionTree, dataFile);

   gInterpreter->GenerateDictionary("vector<vector<int>>", "vector");
   gStyle->SetNdivisions(1020);
   gStyle->SetPaintTextFormat("3.2f");
   gStyle->SetOptFit(1);
   gStyle->SetPalette(6, beamParticle_color);

   string output_name = "thesisPlot_APA3.root";

   TFile *output = new TFile ( output_name.c_str() , "RECREATE");

   //Selected Process and RecoE Int and Inc Histos
   //THIS is what I get from DATA
   
   auto frame = inputFrame
      //Filters
      .Filter("true_beam_endZ >0")
      .Filter("primary_isBeamType && passBeamQuality_TPCjustPosition")
      //Classifications
      .Define("class_true_pion", "if(true_primPionInel == 1) return true; else return false;")
      .Define("class_true_proton", "true_beam_PDG == 2212")
      .Define("class_true_cosmic", "isCosmic")
      .Define("class_true_beamMuon", "primaryMuon")
      .Define("class_true_decayPion", "isDecay && !primaryMuon && !isCosmic")
      .Define("class_true_other", "!class_true_pion && !class_true_proton && !class_true_cosmic && !class_true_beamMuon && !class_true_decayPion")
      ;

   auto dataFrame = data_inputFrame
      .Filter("primary_isBeamType && passBeamQuality_TPCjustPosition");
      //Plot Variables
   
   TH1D* hReco_MC_endZ_pion= new TH1D ("hReco_MC_endZ_pion", "pion;  Reconstructed end Z [cm]; Events / 5 cm", 120, 0. , 600.);
   TH1D* hReco_MC_endZ_proton= new TH1D ("hReco_MC_endZ_proton", "proton;  Reconstructed end Z [cm]; Events / 5 cm", 120, 0. , 600.);
   TH1D* hReco_MC_endZ_beamMuon= new TH1D ("hReco_MC_endZ_beamMuon", "beamMuon;  Reconstructed end Z [cm]; Events / 5 cm", 120, 0. , 600.);
   TH1D* hReco_MC_endZ_cosmic= new TH1D ("hReco_MC_endZ_cosmic", "cosmic;  Reconstructed end Z [cm]; Events / 5 cm", 120, 0. , 600.);
   TH1D* hReco_MC_endZ_decayPion= new TH1D ("hReco_MC_endZ_decayPion", "decayPion;  Reconstructed end Z [cm]; Events / 5 cm", 120, 0. , 600.);
   TH1D* hReco_MC_endZ_other= new TH1D ("hReco_MC_endZ_other", "other;  Reconstructed end Z [cm]; Events / 5 cm", 120, 0. , 600.); 

   TH1D* hReco_58XX_endZ= new TH1D ("hReco_58XX_endZ", "Data;  Reconstructed end Z [cm]; Events / 5 cm", 120, 0. , 600.);
 

   frame
      .Foreach( [ hReco_MC_endZ_pion, hReco_MC_endZ_proton, hReco_MC_endZ_beamMuon, hReco_MC_endZ_cosmic, hReco_MC_endZ_decayPion, hReco_MC_endZ_other]
                
            (double endZ,
             bool pion, bool proton, bool beamMuon, bool cosmic, bool decayPion, bool other){
                  
                  if(pion) {
                  hReco_MC_endZ_pion->Fill(endZ);}
                  else if(proton) {
                  hReco_MC_endZ_proton->Fill(endZ);} 
                  else if(beamMuon) {
                  hReco_MC_endZ_beamMuon->Fill(endZ);} 
                  else if(cosmic) {
                  hReco_MC_endZ_cosmic->Fill(endZ);} 
                  else if(decayPion) {
                  hReco_MC_endZ_decayPion->Fill(endZ); } 
                  else if(other) {
                  hReco_MC_endZ_other->Fill(endZ); }

               },{"reco_beam_endZ", 
                  "class_true_pion", "class_true_proton", "class_true_beamMuon", "class_true_cosmic", "class_true_decayPion", "class_true_other"}
            );

   dataFrame
      .Foreach( [ hReco_58XX_endZ ](double endZ ){

               hReco_58XX_endZ->Fill( endZ );
               },{"reco_beam_endZ"}
            );

   //MC Total Integral
   double totalMC_Z =   hReco_MC_endZ_pion->Integral() + hReco_MC_endZ_proton->Integral()
                      + hReco_MC_endZ_beamMuon->Integral() + hReco_MC_endZ_cosmic->Integral()
                      + hReco_MC_endZ_decayPion->Integral() + hReco_MC_endZ_other->Integral();

   double scale_Z = hReco_58XX_endZ->Integral() / totalMC_Z;

   //Normalise MC to data 
   hReco_MC_endZ_pion->Scale( scale_Z ); hReco_MC_endZ_proton->Scale( scale_Z ); hReco_MC_endZ_beamMuon->Scale( scale_Z );
   hReco_MC_endZ_cosmic->Scale( scale_Z ); hReco_MC_endZ_decayPion->Scale( scale_Z ); hReco_MC_endZ_other->Scale( scale_Z );

   //THStack
   THStack *stack_MC_endZ = new THStack("stack_MC_endZ", "; Reconstructed end Z [cm]; Events / 5 cm");
   stack_MC_endZ->Add( hReco_MC_endZ_pion ); stack_MC_endZ->Add( hReco_MC_endZ_proton ); stack_MC_endZ->Add( hReco_MC_endZ_beamMuon );
   stack_MC_endZ->Add( hReco_MC_endZ_cosmic ); stack_MC_endZ->Add( hReco_MC_endZ_decayPion ); stack_MC_endZ->Add( hReco_MC_endZ_other );

   hReco_58XX_endZ->SetLineColor(kBlack);
   hReco_58XX_endZ->SetMarkerStyle(8);
   hReco_58XX_endZ->SetMarkerSize(0.8);
   
   auto legend_z = new TLegend(0.5,0.65,0.76,0.85);
   legend_z->AddEntry(hReco_58XX_endZ, "PDSP Data, Run 58XX");
   legend_z->AddEntry(hReco_MC_endZ_pion, "Pion"); legend_z->AddEntry(hReco_MC_endZ_proton, "Proton"); legend_z->AddEntry(hReco_MC_endZ_beamMuon, "Beam Muon"); 
   legend_z->AddEntry(hReco_MC_endZ_cosmic, "Cosmic"); legend_z->AddEntry(hReco_MC_endZ_decayPion, "Decay Pion"); legend_z->AddEntry(hReco_MC_endZ_other, "Other"); 
   legend_z->SetTextSize(0.03);

   gStyle->SetOptStat(0);
   gStyle->SetOptFit(0);

   TCanvas* c_endZ_mc_data= new TCanvas("c_endZ_mc_data","c_endZ_mc_data", 1000, 900);
   gPad->SetGrid(1,1);
   hReco_58XX_endZ->SetTitle("");
   stack_MC_endZ->Draw("HIST PFC PMC PLC");
   hReco_58XX_endZ->Draw("PE SAME");
   legend_z->Draw();
   c_endZ_mc_data->Update();
   c_endZ_mc_data->Write();
 
   return 0;
}

