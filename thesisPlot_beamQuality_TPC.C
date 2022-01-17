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


int thesisPlot_beamQuality_TPC(const string mcFile, const string dataFile){

   ROOT::RDataFrame inputFrame(pionTree,mcFile);
   ROOT::RDataFrame data_inputFrame(pionTree, dataFile);

   gInterpreter->GenerateDictionary("vector<vector<int>>", "vector");
   gStyle->SetNdivisions(510);
   gStyle->SetPaintTextFormat("3.2f");
   gStyle->SetOptFit(1);
   gStyle->SetPalette(6, beamParticle_color);

   string output_name = "thesisPlot_beamQuality_TPC.root";

   TFile *output = new TFile ( output_name.c_str() , "RECREATE");

   //Selected Process and RecoE Int and Inc Histos
   //THIS is what I get from DATA
   
   auto frame = inputFrame
      //Filters
      .Filter("true_beam_endZ >0")
      .Filter("primary_isBeamType")
      //Classifications
      .Define("class_true_pion", "if(true_primPionInel == 1) return true; else return false;")
      .Define("class_true_proton", "true_beam_PDG == 2212")
      .Define("class_true_cosmic", "isCosmic")
      .Define("class_true_beamMuon", "primaryMuon")
      .Define("class_true_decayPion", "isDecay && !primaryMuon && !isCosmic")
      .Define("class_true_other", "!class_true_pion && !class_true_proton && !class_true_cosmic && !class_true_beamMuon && !class_true_decayPion")
      //Plot Variables
      .Define("deltaX_divSigma", [](double startX){
            return (startX - mc_meanX)/mc_sigmaX;
            }, {"reco_beam_calo_startX"})
      .Define("deltaY_divSigma", [](double startY){
            return (startY - mc_meanY)/mc_sigmaY;
            }, {"reco_beam_calo_startY"})
      .Define("deltaZ_divSigma", [](double startZ){
            return (startZ - mc_meanZ)/mc_sigmaZ;
            }, {"reco_beam_calo_startZ"});

   auto dataFrame = data_inputFrame
      .Filter("primary_isBeamType")
      //Plot Variables
       .Define("deltaX_divSigma", [](double startX){
            return (startX - data_meanX)/data_sigmaX;
            }, {"reco_beam_calo_startX"})
      .Define("deltaY_divSigma", [](double startY){
            return (startY - data_meanY)/data_sigmaY;
            }, {"reco_beam_calo_startY"})
      .Define("deltaZ_divSigma", [](double startZ){
            return (startZ - data_meanZ)/data_sigmaZ;
            }, {"reco_beam_calo_startZ"});

   TH1D* hReco_MC_deltaX_pion= new TH1D ("hReco_MC_deltaX_pion", "pion; #DeltaX/#sigma_{x}; Events / 0.4", 50, -10. , 10.);
   TH1D* hReco_MC_deltaX_proton= new TH1D ("hReco_MC_deltaX_proton", "proton; #DeltaX/#sigma_{x}; Events / 0.4", 50, -10. , 10.);
   TH1D* hReco_MC_deltaX_beamMuon= new TH1D ("hReco_MC_deltaX_beamMuon", "beamMuon; #DeltaX/#sigma_{x}; Events / 0.4", 50, -10. , 10.);
   TH1D* hReco_MC_deltaX_cosmic= new TH1D ("hReco_MC_deltaX_cosmic", "cosmic; #DeltaX/#sigma_{x}; Events / 0.4", 50, -10. , 10.);
   TH1D* hReco_MC_deltaX_decayPion= new TH1D ("hReco_MC_deltaX_decayPion", "decayPion; #DeltaX/#sigma_{x}; Events / 0.4", 50, -10. , 10.);
   TH1D* hReco_MC_deltaX_other= new TH1D ("hReco_MC_deltaX_other", "other; #DeltaX/#sigma_{x}; Events / 0.4", 50, -10. , 10.); 

   TH1D* hReco_58XX_deltaX= new TH1D ("hReco_58XX_deltaX", "Data; #DeltaX/#sigma_{x}; Events / 0.4", 50, -10. , 10.);
 

   TH1D* hReco_MC_deltaY_pion= new TH1D ("hReco_MC_deltaY_pion", "pion; #DeltaY/#sigma_{y}; Events / 0.4", 50, -10. , 10.);
   TH1D* hReco_MC_deltaY_proton= new TH1D ("hReco_MC_deltaY_proton", "proton; #DeltaY/#sigma_{y}; Events / 0.4", 50, -10. , 10.);
   TH1D* hReco_MC_deltaY_beamMuon= new TH1D ("hReco_MC_deltaY_beamMuon", "beamMuon; #DeltaY/#sigma_{y}; Events / 0.4", 50, -10. , 10.);
   TH1D* hReco_MC_deltaY_cosmic= new TH1D ("hReco_MC_deltaY_cosmic", "cosmic; #DeltaY/#sigma_{y}; Events / 0.4", 50, -10. , 10.);
   TH1D* hReco_MC_deltaY_decayPion= new TH1D ("hReco_MC_deltaY_decayPion", "decayPion; #DeltaY/#sigma_{y}; Events / 0.4", 50, -10. , 10.);
   TH1D* hReco_MC_deltaY_other= new TH1D ("hReco_MC_deltaY_other", "other; #DeltaY/#sigma_{y}; Events / 0.4", 50, -10. , 10.); 

   TH1D* hReco_58XX_deltaY= new TH1D ("hReco_58XX_deltaY", "Data; #DeltaX/#sigma_{y}; Events / 0.4", 50, -10. , 10.);

   TH1D* hReco_MC_deltaZ_pion= new TH1D ("hReco_MC_deltaZ_pion", "pion; #DeltaZ/#sigma_{z}; Events / 0.4", 50, -10. , 10.);
   TH1D* hReco_MC_deltaZ_proton= new TH1D ("hReco_MC_deltaZ_proton", "proton; #DeltaZ/#sigma_{z}; Events / 0.4", 50, -10. , 10.);
   TH1D* hReco_MC_deltaZ_beamMuon= new TH1D ("hReco_MC_deltaZ_beamMuon", "beamMuon; #DeltaZ/#sigma_{z}; Events / 0.4", 50, -10. , 10.);
   TH1D* hReco_MC_deltaZ_cosmic= new TH1D ("hReco_MC_deltaZ_cosmic", "cosmic; #DeltaZ/#sigma_{z}; Events / 0.4", 50, -10. , 10.);
   TH1D* hReco_MC_deltaZ_decayPion= new TH1D ("hReco_MC_deltaZ_decayPion", "decayPion; #DeltaZ/#sigma_{z}; Events / 0.4", 50, -10. , 10.);
   TH1D* hReco_MC_deltaZ_other= new TH1D ("hReco_MC_deltaZ_other", "other; #DeltaZ/#sigma_{z}; Events / 0.4", 50, -10. , 10.); 

   TH1D* hReco_58XX_deltaZ= new TH1D ("hReco_58ZZ_deltaZ", "Data; #DeltaZ/#sigma_{z}; Events / 0.4", 50, -10. , 10.);

   frame
      .Foreach( [ hReco_MC_deltaX_pion, hReco_MC_deltaX_proton, hReco_MC_deltaX_beamMuon, hReco_MC_deltaX_cosmic, hReco_MC_deltaX_decayPion, hReco_MC_deltaX_other,
                  hReco_MC_deltaY_pion, hReco_MC_deltaY_proton, hReco_MC_deltaY_beamMuon, hReco_MC_deltaY_cosmic, hReco_MC_deltaY_decayPion, hReco_MC_deltaY_other,
                  hReco_MC_deltaZ_pion, hReco_MC_deltaZ_proton, hReco_MC_deltaZ_beamMuon, hReco_MC_deltaZ_cosmic, hReco_MC_deltaZ_decayPion, hReco_MC_deltaZ_other ]
                
            (double deltaX_sigma, double deltaY_sigma, double deltaZ_sigma,
             bool pion, bool proton, bool beamMuon, bool cosmic, bool decayPion, bool other){
                  
                  if(pion) {
                  hReco_MC_deltaX_pion->Fill(deltaX_sigma); 
                  hReco_MC_deltaY_pion->Fill(deltaY_sigma);
                  hReco_MC_deltaZ_pion->Fill(deltaZ_sigma);}
                  else if(proton) {
                  hReco_MC_deltaX_proton->Fill(deltaX_sigma); 
                  hReco_MC_deltaY_proton->Fill(deltaY_sigma);
                  hReco_MC_deltaZ_proton->Fill(deltaZ_sigma);}
                  else if(beamMuon) {
                  hReco_MC_deltaX_beamMuon->Fill(deltaX_sigma); 
                  hReco_MC_deltaY_beamMuon->Fill(deltaY_sigma);
                  hReco_MC_deltaZ_beamMuon->Fill(deltaZ_sigma);}
                  else if(cosmic) {
                  hReco_MC_deltaX_cosmic->Fill(deltaX_sigma); 
                  hReco_MC_deltaY_cosmic->Fill(deltaY_sigma);
                  hReco_MC_deltaZ_cosmic->Fill(deltaZ_sigma);}
                  else if(decayPion) {
                  hReco_MC_deltaX_decayPion->Fill(deltaX_sigma); 
                  hReco_MC_deltaY_decayPion->Fill(deltaY_sigma);
                  hReco_MC_deltaZ_decayPion->Fill(deltaZ_sigma);}
                  else if(other) {
                  hReco_MC_deltaX_other->Fill(deltaX_sigma); 
                  hReco_MC_deltaY_other->Fill(deltaY_sigma);
                  hReco_MC_deltaZ_other->Fill(deltaZ_sigma);}

               },{"deltaX_divSigma", "deltaY_divSigma", "deltaZ_divSigma", 
                  "class_true_pion", "class_true_proton", "class_true_beamMuon", "class_true_cosmic", "class_true_decayPion", "class_true_other"}
            );

   dataFrame
      .Foreach( [ hReco_58XX_deltaX, hReco_58XX_deltaY, hReco_58XX_deltaZ](double deltaX_sigma, double deltaY_sigma, double deltaZ_sigma){

               hReco_58XX_deltaX->Fill(deltaX_sigma);
               hReco_58XX_deltaY->Fill(deltaY_sigma);
               hReco_58XX_deltaZ->Fill(deltaZ_sigma);
            
               },{"deltaX_divSigma", "deltaY_divSigma", "deltaZ_divSigma"}
            );

   //MC Total Integral
   double totalMC_X =   hReco_MC_deltaX_pion->Integral() + hReco_MC_deltaX_proton->Integral()
                      + hReco_MC_deltaX_beamMuon->Integral() + hReco_MC_deltaX_cosmic->Integral()
                      + hReco_MC_deltaX_decayPion->Integral() + hReco_MC_deltaX_other->Integral();
   
   double totalMC_Y =   hReco_MC_deltaY_pion->Integral() + hReco_MC_deltaY_proton->Integral()
                      + hReco_MC_deltaY_beamMuon->Integral() + hReco_MC_deltaY_cosmic->Integral()
                      + hReco_MC_deltaY_decayPion->Integral() + hReco_MC_deltaY_other->Integral();
   
   double totalMC_Z =   hReco_MC_deltaZ_pion->Integral() + hReco_MC_deltaZ_proton->Integral()
                      + hReco_MC_deltaZ_beamMuon->Integral() + hReco_MC_deltaZ_cosmic->Integral()
                      + hReco_MC_deltaZ_decayPion->Integral() + hReco_MC_deltaZ_other->Integral();

   double scale_X = hReco_58XX_deltaX->Integral() / totalMC_X;
   double scale_Y= hReco_58XX_deltaY->Integral() / totalMC_Y;
   double scale_Z = hReco_58XX_deltaZ->Integral() / totalMC_Z;

   //Normalise MC to data
   hReco_MC_deltaX_pion->Scale( scale_X); hReco_MC_deltaX_proton->Scale( scale_X ); hReco_MC_deltaX_beamMuon->Scale( scale_X);
   hReco_MC_deltaX_cosmic->Scale( scale_X); hReco_MC_deltaX_decayPion->Scale( scale_X ); hReco_MC_deltaX_other->Scale( scale_X);
   
   hReco_MC_deltaY_pion->Scale( scale_Y); hReco_MC_deltaY_proton->Scale( scale_Y ); hReco_MC_deltaY_beamMuon->Scale( scale_Y);
   hReco_MC_deltaY_cosmic->Scale( scale_Y); hReco_MC_deltaY_decayPion->Scale( scale_Y ); hReco_MC_deltaY_other->Scale( scale_Y);
   
   hReco_MC_deltaZ_pion->Scale( scale_Z); hReco_MC_deltaZ_proton->Scale( scale_Z ); hReco_MC_deltaZ_beamMuon->Scale( scale_Z);
   hReco_MC_deltaZ_cosmic->Scale( scale_Z); hReco_MC_deltaZ_decayPion->Scale( scale_Z ); hReco_MC_deltaZ_other->Scale( scale_Z);

   //THStack
   THStack *stack_MC_deltaX = new THStack("stack_MC_deltaX", "; #Deltax / #sigma_{x}; Events / 0.4");
   stack_MC_deltaX->Add( hReco_MC_deltaX_pion ); stack_MC_deltaX->Add( hReco_MC_deltaX_proton ); stack_MC_deltaX->Add( hReco_MC_deltaX_beamMuon );
   stack_MC_deltaX->Add( hReco_MC_deltaX_cosmic ); stack_MC_deltaX->Add( hReco_MC_deltaX_decayPion ); stack_MC_deltaX->Add( hReco_MC_deltaX_other );
   
   THStack *stack_MC_deltaY = new THStack("stack_MC_deltaY", "; #Deltay / #sigma_{y}; Events / 0.4");
   stack_MC_deltaY->Add( hReco_MC_deltaY_pion ); stack_MC_deltaY->Add( hReco_MC_deltaY_proton ); stack_MC_deltaY->Add( hReco_MC_deltaY_beamMuon );
   stack_MC_deltaY->Add( hReco_MC_deltaY_cosmic ); stack_MC_deltaY->Add( hReco_MC_deltaY_decayPion ); stack_MC_deltaY->Add( hReco_MC_deltaY_other );
   
   THStack *stack_MC_deltaZ = new THStack("stack_MC_deltaZ", "; #Deltaz / #sigma_{z}; Events / 0.4");
   stack_MC_deltaZ->Add( hReco_MC_deltaZ_pion ); stack_MC_deltaZ->Add( hReco_MC_deltaZ_proton ); stack_MC_deltaZ->Add( hReco_MC_deltaZ_beamMuon );
   stack_MC_deltaZ->Add( hReco_MC_deltaZ_cosmic ); stack_MC_deltaZ->Add( hReco_MC_deltaZ_decayPion ); stack_MC_deltaZ->Add( hReco_MC_deltaZ_other );


   hReco_58XX_deltaX->SetLineColor(kBlack);
   hReco_58XX_deltaX->SetMarkerStyle(8);
   hReco_58XX_deltaX->SetMarkerSize(0.6);

   hReco_58XX_deltaY->SetLineColor(kBlack);
   hReco_58XX_deltaY->SetMarkerStyle(8);
   hReco_58XX_deltaY->SetMarkerSize(0.6);

   hReco_58XX_deltaZ->SetLineColor(kBlack);
   hReco_58XX_deltaZ->SetMarkerStyle(8);
   hReco_58XX_deltaZ->SetMarkerSize(0.6);


   auto legend_x = new TLegend(0.18,0.65,0.45,0.85);
   legend_x->AddEntry(hReco_58XX_deltaX, "PDSP Data, Run 58XX");
   legend_x->AddEntry(hReco_MC_deltaX_pion, "Pion"); legend_x->AddEntry(hReco_MC_deltaX_proton, "Proton"); legend_x->AddEntry(hReco_MC_deltaX_beamMuon, "Beam Muon"); 
   legend_x->AddEntry(hReco_MC_deltaX_cosmic, "Cosmic"); legend_x->AddEntry(hReco_MC_deltaX_decayPion, "Decay Pion"); legend_x->AddEntry(hReco_MC_deltaX_other, "Other"); 
   legend_x->SetTextSize(0.03);

   auto legend_y = new TLegend(0.18,0.65,0.45,0.85);
   legend_y->AddEntry(hReco_58XX_deltaY, "PDSP Data, Run 58XX");
   legend_y->AddEntry(hReco_MC_deltaY_pion, "Pion"); legend_y->AddEntry(hReco_MC_deltaY_proton, "Proton"); legend_y->AddEntry(hReco_MC_deltaY_beamMuon, "Beam Muon"); 
   legend_y->AddEntry(hReco_MC_deltaY_cosmic, "Cosmic"); legend_y->AddEntry(hReco_MC_deltaY_decayPion, "Decay Pion"); legend_y->AddEntry(hReco_MC_deltaY_other, "Other"); 
   legend_y->SetTextSize(0.03);
   
   auto legend_z = new TLegend(0.18,0.65,0.45,0.85);
   legend_z->AddEntry(hReco_58XX_deltaZ, "PDSP Data, Run 58XX");
   legend_z->AddEntry(hReco_MC_deltaZ_pion, "Pion"); legend_z->AddEntry(hReco_MC_deltaZ_proton, "Proton"); legend_z->AddEntry(hReco_MC_deltaZ_beamMuon, "Beam Muon"); 
   legend_z->AddEntry(hReco_MC_deltaZ_cosmic, "Cosmic"); legend_z->AddEntry(hReco_MC_deltaZ_decayPion, "Decay Pion"); legend_z->AddEntry(hReco_MC_deltaZ_other, "Other"); 
   legend_z->SetTextSize(0.03);

   gStyle->SetOptStat(0);
   gStyle->SetOptFit(0);

   TCanvas* c_deltaX_mc_data= new TCanvas("c_deltaX_mc_data","c_deltaX_mc_data",1000,800);
   gPad->SetGrid(1,1);
   hReco_58XX_deltaX->SetTitle("");
   stack_MC_deltaX->Draw("HIST PLC PFC PMC");
   hReco_58XX_deltaX->Draw("PE SAME");
   legend_x->Draw();
   c_deltaX_mc_data->Update();
   c_deltaX_mc_data->Write();

   TCanvas* c_deltaY_mc_data= new TCanvas("c_deltaY_mc_data","c_deltaY_mc_data", 1000, 800);
   gPad->SetGrid(1,1);
   hReco_58XX_deltaY->SetTitle("");
   stack_MC_deltaY->Draw("HIST PFC PLC PMC");
   hReco_58XX_deltaY->Draw("PE SAME");
   legend_y->Draw();
   c_deltaY_mc_data->Update();
   c_deltaY_mc_data->Write();

   TCanvas* c_deltaZ_mc_data= new TCanvas("c_deltaZ_mc_data","c_deltaZ_mc_data", 1000, 800);
   gPad->SetGrid(1,1);
   hReco_58XX_deltaZ->SetTitle("");
   stack_MC_deltaZ->Draw("HIST PFC PMC PLC");
   hReco_58XX_deltaZ->Draw("PE SAME");
   legend_z->Draw();
   c_deltaZ_mc_data->Update();
   c_deltaZ_mc_data->Write();
 
   return 0;
}

