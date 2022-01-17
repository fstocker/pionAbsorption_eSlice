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

int thesisPlot_secondaries(const string mcFile, const string dataFile){

   ROOT::RDataFrame inputFrame(pionTree,mcFile);
   ROOT::RDataFrame data_inputFrame(pionTree, dataFile);

   gInterpreter->GenerateDictionary("vector<vector<int>>", "vector");
   //gStyle->SetNdivisions(1020);
   gStyle->SetPaintTextFormat("3.2f");
   gStyle->SetOptFit(1);
   gStyle->SetPalette(7, secondary_color);

   //daughter types
   //1 - pion, 2- muon, 3-proton, 4-photon, 5-nucleus, 6-electron, 7-other

   string output_name = "thesisPlot_secondaries.root";

   TFile *output = new TFile ( output_name.c_str() , "RECREATE");

   //Selected Process and RecoE Int and Inc Histos
   //THIS is what I get from DATA
   
   auto frame = inputFrame
      //Filters
      .Filter("true_beam_endZ >0")
      .Filter("selected_incidentPion");
      //Classifications

   auto dataFrame = data_inputFrame
      .Filter("selected_incidentPion");
      //Plot Variables

   int nBin_dEdX = 40; double low_dEdX = 0; double high_dEdX = 8;
   int nBin_chi2 = 100; double low_chi2 = 0; double high_chi2 = 400;
   int nBin_cnn = 50; double low_cnn = 0; double high_cnn = 1;
   int nBin_nHits = 100; double low_nHits = 0; double high_nHits = 1000;

   TH1D* hReco_secondary_dEdX_pion = new TH1D("hReco_secondary_dEdX_pion", "pion; dEdX [MeV/cm];Events / 0.2 MeV", nBin_dEdX, low_dEdX, high_dEdX);
   TH1D* hReco_secondary_dEdX_muon = new TH1D("hReco_secondary_dEdX_muon", "muon; dEdX [MeV/cm];Events / 0.2 MeV", nBin_dEdX, low_dEdX, high_dEdX);
   TH1D* hReco_secondary_dEdX_proton = new TH1D("hReco_secondary_dEdX_proton", "proton; dEdX [MeV/cm];Events / 0.2 MeV", nBin_dEdX, low_dEdX, high_dEdX);
   TH1D* hReco_secondary_dEdX_photon = new TH1D("hReco_secondary_dEdX_photon", "photon; dEdX [MeV/cm];Events / 0.2 MeV", nBin_dEdX, low_dEdX, high_dEdX);
   TH1D* hReco_secondary_dEdX_nucleus = new TH1D("hReco_secondary_dEdX_nucleus", "nucleus; dEdX [MeV/cm];Events / 0.2 MeV", nBin_dEdX, low_dEdX, high_dEdX);
   TH1D* hReco_secondary_dEdX_electron = new TH1D("hReco_secondary_dEdX_electron", "electron; dEdX [MeV/cm];Events / 0.2 MeV", nBin_dEdX, low_dEdX, high_dEdX);
   TH1D* hReco_secondary_dEdX_other = new TH1D("hReco_secondary_dEdX_other", "other; dEdX [MeV/cm];Events / 0.2 MeV", nBin_dEdX, low_dEdX, high_dEdX);

   TH1D* hReco_58XX_dEdX = new TH1D("hReco_58XX_dEdX", "PDSP Data, Run 58XX; dEdX [MeV/cm];Events / 0.2 MeV", nBin_dEdX, low_dEdX, high_dEdX);
   TH1D* hReco_secondary_dEdX_help = new TH1D("hReco_secondary_dEdX_help", "; dEdX [MeV/cm]; Proton / (Proton + Pion + Muon)", nBin_dEdX, low_dEdX, high_dEdX);
   TH1D* hReco_secondary_dEdX_ratio = new TH1D("hReco_secondary_dEdX_ratio", "; dEdX [MeV/cm]; Proton / (Proton + Pion + Muon)", nBin_dEdX, low_dEdX, high_dEdX);
 
   TH1D* hReco_secondary_chi2_pion = new TH1D("hReco_secondary_chi2_pion", "pion; chi2/ndof; Events / 4 #chi^2/ndof", nBin_chi2, low_chi2, high_chi2);
   TH1D* hReco_secondary_chi2_muon = new TH1D("hReco_secondary_chi2_muon", "muon; chi2/ndof; Events / 4 #chi^2/ndof", nBin_chi2, low_chi2, high_chi2);
   TH1D* hReco_secondary_chi2_proton = new TH1D("hReco_secondary_chi2_proton", "proton; chi2/ndof; Events / 4 #chi^2/ndof", nBin_chi2, low_chi2, high_chi2);
   TH1D* hReco_secondary_chi2_photon = new TH1D("hReco_secondary_chi2_photon", "photon; chi2/ndof; Events / 4 #chi^2/ndof", nBin_chi2, low_chi2, high_chi2);
   TH1D* hReco_secondary_chi2_nucleus = new TH1D("hReco_secondary_chi2_nucleus", "nucleus; chi2/ndof; Events / 4 #chi^2/ndof", nBin_chi2, low_chi2, high_chi2);
   TH1D* hReco_secondary_chi2_electron = new TH1D("hReco_secondary_chi2_electron", "electron; chi2/ndof; Events / 4 #chi^2/ndof", nBin_chi2, low_chi2, high_chi2);
   TH1D* hReco_secondary_chi2_other = new TH1D("hReco_secondary_chi2_other", "other; chi2/ndof; Events / 4 #chi^2/ndof", nBin_chi2, low_chi2, high_chi2);

   TH1D* hReco_58XX_chi2 = new TH1D("hReco_58XX_chi2", "PDSP Data, Run 58XX; chi2/ndof; Events / 4 #chi^2/ndof", nBin_chi2, low_chi2, high_chi2);
 
   TH1D* hReco_secondary_cnn_pion = new TH1D("hReco_secondary_cnn_pion", "pion; CNN TrackScore; Events / 0.02", nBin_cnn, low_cnn, high_cnn);
   TH1D* hReco_secondary_cnn_muon = new TH1D("hReco_secondary_cnn_muon", "muon; CNN TrackScore; Events / 0.02", nBin_cnn, low_cnn, high_cnn);
   TH1D* hReco_secondary_cnn_proton = new TH1D("hReco_secondary_cnn_proton", "proton; CNN TrackScore; Events / 0.02", nBin_cnn, low_cnn, high_cnn);
   TH1D* hReco_secondary_cnn_photon = new TH1D("hReco_secondary_cnn_photon", "photon; CNN TrackScore; Events / 0.02", nBin_cnn, low_cnn, high_cnn);
   TH1D* hReco_secondary_cnn_nucleus = new TH1D("hReco_secondary_cnn_nucleus", "nucleus; CNN TrackScore; Events / 0.02", nBin_cnn, low_cnn, high_cnn);
   TH1D* hReco_secondary_cnn_electron = new TH1D("hReco_secondary_cnn_electron", "electron; CNN TrackScore; Events / 0.02", nBin_cnn, low_cnn, high_cnn);
   TH1D* hReco_secondary_cnn_other = new TH1D("hReco_secondary_cnn_other", "other; CNN TrackScore; Events / 0.02", nBin_cnn, low_cnn, high_cnn);

   TH1D* hReco_58XX_cnn = new TH1D("hReco_58XX_cnn", "PDSP Data, Run 58XX; CNN TrackScore; Events / 0.02", nBin_cnn, low_cnn, high_cnn);
 
   TH1D* hReco_secondary_nHits_pion = new TH1D("hReco_secondary_nHits_pion", "pion; nHits; Events / 10 hits", nBin_nHits, low_nHits, high_nHits);
   TH1D* hReco_secondary_nHits_muon = new TH1D("hReco_secondary_nHits_muon", "muon; nHits; Events / 10 hits", nBin_nHits, low_nHits, high_nHits);
   TH1D* hReco_secondary_nHits_proton = new TH1D("hReco_secondary_nHits_proton", "proton; nHits; Events / 10 hits", nBin_nHits, low_nHits, high_nHits);
   TH1D* hReco_secondary_nHits_photon = new TH1D("hReco_secondary_nHits_photon", "photon; nHits; Events / 10 hits", nBin_nHits, low_nHits, high_nHits);
   TH1D* hReco_secondary_nHits_nucleus = new TH1D("hReco_secondary_nHits_nucleus", "nucleus; nHits; Events / 10 hits", nBin_nHits, low_nHits, high_nHits);
   TH1D* hReco_secondary_nHits_electron = new TH1D("hReco_secondary_nHits_electron", "electron; nHits; Events / 10 hits", nBin_nHits, low_nHits, high_nHits);
   TH1D* hReco_secondary_nHits_other = new TH1D("hReco_secondary_nHits_other", "other; nHits; Events / 10 hits", nBin_nHits, low_nHits, high_nHits);

   TH1D* hReco_58XX_nHits = new TH1D("hReco_58XX_nHits", "PDSP Data, Run 58XX; nHits; Events / 10 hits", nBin_nHits, low_nHits, high_nHits);
 
 
   frame
      .Foreach( [ hReco_secondary_dEdX_pion, hReco_secondary_dEdX_muon, hReco_secondary_dEdX_proton, hReco_secondary_dEdX_photon, 
                  hReco_secondary_dEdX_nucleus, hReco_secondary_dEdX_electron, hReco_secondary_dEdX_other, 
                  hReco_secondary_chi2_pion, hReco_secondary_chi2_muon, hReco_secondary_chi2_proton, hReco_secondary_chi2_photon, 
                  hReco_secondary_chi2_nucleus, hReco_secondary_chi2_electron, hReco_secondary_chi2_other, 
                  hReco_secondary_cnn_pion, hReco_secondary_cnn_muon, hReco_secondary_cnn_proton, hReco_secondary_cnn_photon, 
                  hReco_secondary_cnn_nucleus, hReco_secondary_cnn_electron, hReco_secondary_cnn_other, 
                  hReco_secondary_nHits_pion, hReco_secondary_nHits_muon, hReco_secondary_nHits_proton, hReco_secondary_nHits_photon, 
                  hReco_secondary_nHits_nucleus, hReco_secondary_nHits_electron, hReco_secondary_nHits_other ] 
                
            (std::vector<int> &daughter_PDG,
             std::vector<double> &dEdX, std::vector<double> &chi2, std::vector<int> &ndof, 
             std::vector<double> &cnn, std::vector<int> nHits ){
            
               for(size_t i = 0; i < daughter_PDG.size(); i++){

                  int daugh_type = daughter_PDG[i];
                  
                  if(daugh_type == 1){
                     //For All Daughters Fill TrackScore Plot
                     hReco_secondary_cnn_pion->Fill( cnn[i] );
                     
                     if( cnn[i] > cut_trackScore ) { //Track-Like
                        hReco_secondary_dEdX_pion->Fill( dEdX[i] );
                        hReco_secondary_chi2_pion->Fill( chi2[i]/ndof[i] );
                     }
                     else if ( cnn[i] > 0 && cnn[i] <= cut_trackScore ) hReco_secondary_nHits_pion->Fill( nHits[i] ); //Shower Like
                  }
                  else if(daugh_type == 2){
                     //For All Daughters Fill TrackScore Plot
                     hReco_secondary_cnn_muon->Fill( cnn[i] );
                     
                     if( cnn[i] > cut_trackScore ) { //Track-Like
                        hReco_secondary_dEdX_muon->Fill( dEdX[i] );
                        hReco_secondary_chi2_muon->Fill( chi2[i]/ndof[i] );
                     }
                     else if ( cnn[i] > 0 && cnn[i] <= cut_trackScore ) hReco_secondary_nHits_muon->Fill( nHits[i] ); //Shower Like
                  }
                  else if(daugh_type == 3){
                     hReco_secondary_cnn_proton->Fill( cnn[i] );
                     
                     if( cnn[i] > cut_trackScore ) { //Track-Like
                        hReco_secondary_dEdX_proton->Fill( dEdX[i] );
                        hReco_secondary_chi2_proton->Fill( chi2[i]/ndof[i] );
                     }
                     else if ( cnn[i] > 0 && cnn[i] <= cut_trackScore ) hReco_secondary_nHits_proton->Fill( nHits[i] ); //Shower Like
                  }
                  else if(daugh_type == 4){
                     hReco_secondary_cnn_photon->Fill( cnn[i] );
                     
                     if( cnn[i] > cut_trackScore ) { //Track-Like
                        hReco_secondary_dEdX_photon->Fill( dEdX[i] );
                        hReco_secondary_chi2_photon->Fill( chi2[i]/ndof[i] );
                     }
                     else if ( cnn[i] > 0 && cnn[i] <= cut_trackScore ) hReco_secondary_nHits_photon->Fill( nHits[i] ); //Shower Like
                  }
                  else if(daugh_type == 5){
                     hReco_secondary_cnn_nucleus->Fill( cnn[i] );
                     
                     if( cnn[i] > cut_trackScore ) { //Track-Like
                        hReco_secondary_dEdX_nucleus->Fill( dEdX[i] );
                        hReco_secondary_chi2_nucleus->Fill( chi2[i]/ndof[i] );
                     }
                     else if ( cnn[i] > 0 && cnn[i] <= cut_trackScore ) hReco_secondary_nHits_nucleus->Fill( nHits[i] ); //Shower Like
                  }
                  else if(daugh_type == 6){
                     hReco_secondary_cnn_electron->Fill( cnn[i] );
                     
                     if( cnn[i] > cut_trackScore ) { //Track-Like
                        hReco_secondary_dEdX_electron->Fill( dEdX[i] );
                        hReco_secondary_chi2_electron->Fill( chi2[i]/ndof[i] );
                     }
                     else if ( cnn[i] > 0 && cnn[i] <= cut_trackScore ) hReco_secondary_nHits_electron->Fill( nHits[i] ); //Shower Like
                  }
                  else if(daugh_type == 7){
                     hReco_secondary_cnn_other->Fill( cnn[i] );
                     
                     if( cnn[i] > cut_trackScore ) { //Track-Like
                        hReco_secondary_dEdX_other->Fill( dEdX[i] );
                        hReco_secondary_chi2_other->Fill( chi2[i]/ndof[i] );
                     }
                     else if ( cnn[i] > 0 && cnn[i] <= cut_trackScore ) hReco_secondary_nHits_other->Fill( nHits[i] ); //Shower Like
                  }

               }



               },{"daughter_PDGs_types", 
                  "reco_daughter_allTrack_truncLibo_dEdX", "reco_daughter_allTrack_Chi2_proton", 
                  "reco_daughter_allTrack_Chi2_ndof", "reco_daughter_PFP_trackScore_collection", "reco_daughter_PFP_nHits"}
            );

   dataFrame
      .Foreach( [ hReco_58XX_dEdX, hReco_58XX_chi2, hReco_58XX_cnn, hReco_58XX_nHits ]
                (std::vector<double> &dEdX, std::vector<double> &chi2, std::vector<int> &ndof, 
                 std::vector<double> &cnn, std::vector<int> nHits){

               for(size_t i = 0; i < nHits.size(); i++){
                     
                     hReco_58XX_cnn->Fill( cnn[i] );
                     if( cnn[i] > cut_trackScore ){
                        hReco_58XX_dEdX->Fill( dEdX[i] );
                        hReco_58XX_chi2->Fill( chi2[i]/ndof[i] );
                     }
                     else if( cnn[i] > 0 && cnn[i] <= cut_trackScore ) hReco_58XX_nHits->Fill( nHits[i] );
                  }
               
               },{ "reco_daughter_allTrack_truncLibo_dEdX", "reco_daughter_allTrack_Chi2_proton", 
                   "reco_daughter_allTrack_Chi2_ndof", "reco_daughter_PFP_trackScore_collection", "reco_daughter_PFP_nHits"}            );

   //Ratio dEdX
   hReco_secondary_dEdX_help->Add(hReco_secondary_dEdX_proton);
   hReco_secondary_dEdX_help->Add(hReco_secondary_dEdX_pion);
   hReco_secondary_dEdX_help->Add(hReco_secondary_dEdX_muon);
   hReco_secondary_dEdX_ratio->Divide(hReco_secondary_dEdX_proton , hReco_secondary_dEdX_help );

   //MC Total Integral
   double totalMC_dEdX =   hReco_secondary_dEdX_pion->Integral() + hReco_secondary_dEdX_muon->Integral()
                      + hReco_secondary_dEdX_proton->Integral() + hReco_secondary_dEdX_photon->Integral()
                      + hReco_secondary_dEdX_nucleus->Integral() + hReco_secondary_dEdX_electron->Integral() + hReco_secondary_dEdX_other->Integral();

   double scale_dEdX = hReco_58XX_dEdX->Integral() / totalMC_dEdX;

   //MC Total Integral
   double totalMC_chi2 =   hReco_secondary_chi2_pion->Integral() + hReco_secondary_chi2_muon->Integral()
                      + hReco_secondary_chi2_proton->Integral() + hReco_secondary_chi2_photon->Integral()
                      + hReco_secondary_chi2_nucleus->Integral() + hReco_secondary_chi2_electron->Integral() + hReco_secondary_chi2_other->Integral();

   double scale_chi2 = hReco_58XX_chi2->Integral() / totalMC_chi2;
   
   //MC Total Integral
   double totalMC_cnn =   hReco_secondary_cnn_pion->Integral() + hReco_secondary_cnn_muon->Integral()
                      + hReco_secondary_cnn_proton->Integral() + hReco_secondary_cnn_photon->Integral()
                      + hReco_secondary_cnn_nucleus->Integral() + hReco_secondary_cnn_electron->Integral() + hReco_secondary_cnn_other->Integral();

   double scale_cnn = hReco_58XX_cnn->Integral() / totalMC_cnn;
   
   //MC Total Integral
   double totalMC_nHits =   hReco_secondary_nHits_pion->Integral() + hReco_secondary_nHits_muon->Integral()
                      + hReco_secondary_nHits_proton->Integral() + hReco_secondary_nHits_photon->Integral()
                      + hReco_secondary_nHits_nucleus->Integral() + hReco_secondary_nHits_electron->Integral() + hReco_secondary_nHits_other->Integral();

   double scale_nHits = hReco_58XX_nHits->Integral() / totalMC_nHits;

   //Normalise MC to data 
   hReco_secondary_dEdX_pion->Scale( scale_dEdX ); hReco_secondary_dEdX_muon->Scale( scale_dEdX ); hReco_secondary_dEdX_proton->Scale( scale_dEdX );
   hReco_secondary_dEdX_photon->Scale( scale_dEdX ); hReco_secondary_dEdX_nucleus->Scale( scale_dEdX ); hReco_secondary_dEdX_electron->Scale( scale_dEdX ); 
   hReco_secondary_dEdX_other->Scale( scale_dEdX );
   
   //Normalise MC to data 
   hReco_secondary_chi2_pion->Scale( scale_chi2 ); hReco_secondary_chi2_muon->Scale( scale_chi2 ); hReco_secondary_chi2_proton->Scale( scale_chi2 );
   hReco_secondary_chi2_photon->Scale( scale_chi2 ); hReco_secondary_chi2_nucleus->Scale( scale_chi2 ); hReco_secondary_chi2_electron->Scale( scale_chi2 ); 
   hReco_secondary_chi2_other->Scale( scale_chi2 );
   
   //Normalise MC to data 
   hReco_secondary_cnn_pion->Scale( scale_cnn ); hReco_secondary_cnn_muon->Scale( scale_cnn ); hReco_secondary_cnn_proton->Scale( scale_cnn );
   hReco_secondary_cnn_photon->Scale( scale_cnn ); hReco_secondary_cnn_nucleus->Scale( scale_cnn ); hReco_secondary_cnn_electron->Scale( scale_cnn ); 
   hReco_secondary_cnn_other->Scale( scale_cnn );
   
   //Normalise MC to data 
   hReco_secondary_nHits_pion->Scale( scale_nHits ); hReco_secondary_nHits_muon->Scale( scale_nHits ); hReco_secondary_nHits_proton->Scale( scale_nHits );
   hReco_secondary_nHits_photon->Scale( scale_nHits ); hReco_secondary_nHits_nucleus->Scale( scale_nHits ); hReco_secondary_nHits_electron->Scale( scale_nHits ); 
   hReco_secondary_nHits_other->Scale( scale_nHits );

   //THStack
   THStack *stack_secondary_dEdX = new THStack("stack_secondary_dEdX", ";dEdX [MeV/cm];Events / 0.2 MeV");
   stack_secondary_dEdX->Add( hReco_secondary_dEdX_pion ); stack_secondary_dEdX->Add( hReco_secondary_dEdX_muon ); 
   stack_secondary_dEdX->Add( hReco_secondary_dEdX_proton ); stack_secondary_dEdX->Add( hReco_secondary_dEdX_photon ); 
   stack_secondary_dEdX->Add( hReco_secondary_dEdX_nucleus ); stack_secondary_dEdX->Add( hReco_secondary_dEdX_electron );
   stack_secondary_dEdX->Add( hReco_secondary_dEdX_other );

   hReco_58XX_dEdX->SetLineColor(kBlack);
   hReco_58XX_dEdX->SetMarkerStyle(8);
   hReco_58XX_dEdX->SetMarkerSize(0.6);
   
   auto legend_dEdX = new TLegend(0.5,0.5,0.8,0.85);
   legend_dEdX->AddEntry(hReco_58XX_dEdX, "PDSP Data, Run 58XX");
   legend_dEdX->AddEntry(hReco_secondary_dEdX_pion, "Pion daughter"); legend_dEdX->AddEntry(hReco_secondary_dEdX_muon, "Muon daughter"); 
   legend_dEdX->AddEntry(hReco_secondary_dEdX_proton, "Proton daughter"); legend_dEdX->AddEntry(hReco_secondary_dEdX_photon, "Photon daughter"); 
   legend_dEdX->AddEntry(hReco_secondary_dEdX_nucleus, "Nucleus daughter"); legend_dEdX->AddEntry(hReco_secondary_dEdX_electron, "Electron daughter"); 
   legend_dEdX->AddEntry(hReco_secondary_dEdX_other, "Other"); 
   legend_dEdX->SetTextSize(0.04);

   //THStack
   THStack *stack_secondary_chi2 = new THStack("stack_secondary_chi2", "; #chi^{2}/ndof; Events / 4 #chi^{2}/ndof");
   stack_secondary_chi2->Add( hReco_secondary_chi2_pion ); stack_secondary_chi2->Add( hReco_secondary_chi2_muon ); 
   stack_secondary_chi2->Add( hReco_secondary_chi2_proton ); stack_secondary_chi2->Add( hReco_secondary_chi2_photon ); 
   stack_secondary_chi2->Add( hReco_secondary_chi2_nucleus ); stack_secondary_chi2->Add( hReco_secondary_chi2_electron );
   stack_secondary_chi2->Add( hReco_secondary_chi2_other );

   hReco_58XX_chi2->SetLineColor(kBlack);
   hReco_58XX_chi2->SetMarkerStyle(8);
   hReco_58XX_chi2->SetMarkerSize(0.8);
   
   auto legend_chi2 = new TLegend(0.5,0.6,0.8,0.85);
   legend_chi2->AddEntry(hReco_58XX_chi2, "PDSP Data, Run 58XX");
   legend_chi2->AddEntry(hReco_secondary_chi2_pion, "Pion daughter"); legend_chi2->AddEntry(hReco_secondary_chi2_muon, "Muon daughter"); 
   legend_chi2->AddEntry(hReco_secondary_chi2_proton, "Proton daughter"); legend_chi2->AddEntry(hReco_secondary_chi2_photon, "Photon daughter"); 
   legend_chi2->AddEntry(hReco_secondary_chi2_nucleus, "Nucleus daughter"); legend_chi2->AddEntry(hReco_secondary_chi2_electron, "Electron daughter"); 
   legend_chi2->AddEntry(hReco_secondary_chi2_other, "Other"); 
   legend_chi2->SetTextSize(0.03);
   
   //THStack
   THStack *stack_secondary_cnn = new THStack("stack_secondary_cnn", "; CNN TrackScore; Events / 0.02");
   stack_secondary_cnn->Add( hReco_secondary_cnn_pion ); stack_secondary_cnn->Add( hReco_secondary_cnn_muon ); 
   stack_secondary_cnn->Add( hReco_secondary_cnn_proton ); stack_secondary_cnn->Add( hReco_secondary_cnn_photon ); 
   stack_secondary_cnn->Add( hReco_secondary_cnn_nucleus ); stack_secondary_cnn->Add( hReco_secondary_cnn_electron );
   stack_secondary_cnn->Add( hReco_secondary_cnn_other );

   hReco_58XX_cnn->SetLineColor(kBlack);
   hReco_58XX_cnn->SetMarkerStyle(8);
   hReco_58XX_cnn->SetMarkerSize(0.8);
   
   auto legend_cnn = new TLegend(0.5,0.6,0.8,0.85);
   legend_cnn->AddEntry(hReco_58XX_cnn, "PDSP Data, Run 58XX");
   legend_cnn->AddEntry(hReco_secondary_cnn_pion, "Pion daughter"); legend_cnn->AddEntry(hReco_secondary_cnn_muon, "Muon daughter"); 
   legend_cnn->AddEntry(hReco_secondary_cnn_proton, "Proton daughter"); legend_cnn->AddEntry(hReco_secondary_cnn_photon, "Photon daughter"); 
   legend_cnn->AddEntry(hReco_secondary_cnn_nucleus, "Nucleus daughter"); legend_cnn->AddEntry(hReco_secondary_cnn_electron, "Electron daughter"); 
   legend_cnn->AddEntry(hReco_secondary_cnn_other, "Other"); 
   legend_cnn->SetTextSize(0.03);

   //THStack
   THStack *stack_secondary_nHits = new THStack("stack_secondary_nHits", "; Number of hits; Events / 10 hits");
   stack_secondary_nHits->Add( hReco_secondary_nHits_pion ); stack_secondary_nHits->Add( hReco_secondary_nHits_muon ); 
   stack_secondary_nHits->Add( hReco_secondary_nHits_proton ); stack_secondary_nHits->Add( hReco_secondary_nHits_photon ); 
   stack_secondary_nHits->Add( hReco_secondary_nHits_nucleus ); stack_secondary_nHits->Add( hReco_secondary_nHits_electron );
   stack_secondary_nHits->Add( hReco_secondary_nHits_other );

   hReco_58XX_nHits->SetLineColor(kBlack);
   hReco_58XX_nHits->SetMarkerStyle(8);
   hReco_58XX_nHits->SetMarkerSize(0.8);
   
   auto legend_nHits = new TLegend(0.5,0.6,0.8,0.85);
   legend_nHits->AddEntry(hReco_58XX_nHits, "PDSP Data, Run 58XX");
   legend_nHits->AddEntry(hReco_secondary_nHits_pion, "Pion daughter"); legend_nHits->AddEntry(hReco_secondary_nHits_muon, "Muon daughter"); 
   legend_nHits->AddEntry(hReco_secondary_nHits_proton, "Proton daughter"); legend_nHits->AddEntry(hReco_secondary_nHits_photon, "Photon daughter"); 
   legend_nHits->AddEntry(hReco_secondary_nHits_nucleus, "Nucleus daughter"); legend_nHits->AddEntry(hReco_secondary_nHits_electron, "Electron daughter"); 
   legend_nHits->AddEntry(hReco_secondary_nHits_other, "Other"); 
   legend_nHits->SetTextSize(0.03);

   gStyle->SetOptStat(0);
   gStyle->SetOptFit(0);

   TCanvas* c_dEdX_mc_data= new TCanvas("c_dEdX_mc_data","c_dEdX_mc_data", 1100, 900);
   c_dEdX_mc_data->Divide(1,2);
   c_dEdX_mc_data->cd(1);
   gPad->SetGrid(1,1);
   gPad->SetBottomMargin(0.01);
   stack_secondary_dEdX->Draw("HIST PFC PMC PLC");
   stack_secondary_dEdX->GetXaxis()->SetNdivisions(510);
   stack_secondary_dEdX->GetYaxis()->SetLabelSize(0.04);
   stack_secondary_dEdX->GetYaxis()->SetTitleSize(0.04);
   stack_secondary_dEdX->GetHistogram()->GetXaxis()->SetTickLength(0);
   stack_secondary_dEdX->GetHistogram()->GetXaxis()->SetLabelOffset(999);
   hReco_58XX_dEdX->Draw("PE SAME");
   legend_dEdX->Draw();
   c_dEdX_mc_data->cd(2);
   gPad->SetGrid(1,1);
   gPad->SetTopMargin(0.01);
   hReco_secondary_dEdX_ratio->Draw("HIST");
   hReco_secondary_dEdX_ratio->SetLineWidth(3);
   hReco_secondary_dEdX_ratio->GetXaxis()->SetNdivisions(510);
   hReco_secondary_dEdX_ratio->GetXaxis()->SetLabelSize(0.04);
   hReco_secondary_dEdX_ratio->GetYaxis()->SetLabelSize(0.04);
   hReco_secondary_dEdX_ratio->GetXaxis()->SetTitleSize(0.06);
   hReco_secondary_dEdX_ratio->GetYaxis()->SetTitleSize(0.04);
   c_dEdX_mc_data->Update();
   c_dEdX_mc_data->Write();

   TCanvas* c_chi2_mc_data= new TCanvas("c_chi2_mc_data","c_chi2_mc_data", 1200, 800);
   gPad->SetGrid(1,1);
   gPad->SetLogy();
   stack_secondary_chi2->Draw("HIST PFC PMC PLC");
   stack_secondary_chi2->GetXaxis()->SetNdivisions(510);
   hReco_58XX_chi2->Draw("PE SAME");
   legend_chi2->Draw();
   c_chi2_mc_data->Update();
   c_chi2_mc_data->Write();

   TCanvas* c_cnn_mc_data= new TCanvas("c_cnn_mc_data","c_cnn_mc_data", 1100, 800);
   gPad->SetGrid(1,1);
   stack_secondary_cnn->Draw("HIST PFC PMC PLC");
   stack_secondary_cnn->GetXaxis()->SetNdivisions(510);
   hReco_58XX_cnn->Draw("PE SAME");
   legend_cnn->Draw();
   c_cnn_mc_data->Update();
   c_cnn_mc_data->Write();
 
   TCanvas* c_nHits_mc_data= new TCanvas("c_nHits_mc_data","c_nHits_mc_data", 1100, 800);
   gPad->SetGrid(1,1);
   stack_secondary_nHits->Draw("HIST PFC PMC PLC");
   stack_secondary_nHits->GetXaxis()->SetNdivisions(1020);
   hReco_58XX_nHits->Draw("PE SAME");
   legend_nHits->Draw();
   c_nHits_mc_data->Update();
   c_nHits_mc_data->Write();
 
   return 0;
}

