#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
#include <iostream>
using std::cout;
using std::endl;
using namespace std;
using namespace ROOT::VecOps;

#include "TRandom.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TRatioPlot.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TH1.h"
#include "TF1.h"
#include "THStack.h"
#include "TLegend.h"
#include "TArrow.h"
#include "TStyle.h"
#include "TRandom3.h"
#include "TColor.h"
#include "TLatex.h"
#include "TMath.h"
#include "TMatrixDBase.h"
#include "TArray.h"
#include <ROOT/RDataFrame.hxx>


#include <iostream>
#include <math.h>
#include <string.h>
#include <stdio.h>
#endif

#include "../lambda.h"
#include "../betheBloch.h"
#include "eSlice.h"
#include "plotThesis_config.h"

//==============================================================================
// Global definitions
//==============================================================================

const Double_t cutdummy= -99999.0;
const Double_t dataMC_events = 70000; //root file DATA has 70121 ev, root file MC has 74149 ev


void doRatioPlot_pion( TH1D* data, TH1D* mc, TH1D* pion, TH1D* decayPion, TH1D* muon, TH1D* proton, TH1D* cosmic, TH1D* other, string c_name, string title, string png ){

   gStyle->SetPalette(5, beamParticle_noProt_color);

   double max = mc->GetBinContent(mc->GetMaximumBin());

   THStack *stack = new THStack("stack", "");
   stack->Add( pion );
   //stack->Add( proton );
   stack->Add( muon );
   stack->Add( cosmic );
   stack->Add( decayPion );
   stack->Add( other);

   TCanvas *c_ratio = new TCanvas( c_name.c_str(), c_name.c_str(), 1200,900 );
   c_ratio->SetGrid(1,1);
   gStyle->SetOptStat(0);
   //c_ratio->SetTicks(0,1);
   auto rp = new TRatioPlot( data, mc , "divsym");
   rp->SetH1DrawOpt("PE");
   rp->SetH2DrawOpt("HIST");
   rp->Draw("SAME");
   rp->GetUpperPad()->cd();
   stack->Draw("HIST PFC PMC PLC SAME");
   stack->SetTitle( title.c_str() );
   
   data->SetMarkerColor(kBlack);
   data->SetLineColor(kBlack);
   data->SetMarkerStyle(8);
   data->SetMarkerSize(0.8);
   data->SetLineWidth(2);
   data->GetYaxis()->SetRangeUser(0, max + 0.1*max);
   data->Draw("PE1 SAME");
   data->SetTitle( title.c_str() );
   data->SetTitleSize( 5);
   rp->GetUpperPad()->SetTitle( title.c_str() );
   //rp->GetUpperPad()->BuildLegend();
   rp->GetUpperPad()->SetGridx();
   rp->GetLowerRefXaxis()->SetNdivisions(1020);
   rp->GetLowerRefXaxis()->SetLabelSize(0.025);
   rp->GetLowYaxis()->SetNdivisions(505);
   rp->GetLowerRefGraph()->SetMinimum(0);
   rp->GetLowerRefGraph()->SetMaximum(2.5);
   rp->GetLowerRefYaxis()->SetTitle("Data / MC");
   rp->GetLowerRefGraph()->SetLineColor(kBlack);
   //rp->Update();
   //

   auto leg = new TLegend(0.5, 0.4, 0.9, 0.8);
   leg->AddEntry( data, "PDSP Data, Run 58XX"); leg->AddEntry( pion, "Pion"); //leg->AddEntry( proton, "Proton"); 
   leg->AddEntry( muon, "Beam muon");  leg->AddEntry( cosmic, "Cosmic"); leg->AddEntry( decayPion, "Pion decay" );
   leg->AddEntry( other, "Other");
   leg->SetTextSize(0.04);

   rp->GetUpperPad()->cd();
   leg->Draw();


   c_ratio->Update();
   //c_ratio->SaveAs( png.c_str() );

};

void doRatioPlot_pion_abs( TH1D* data, TH1D* mc, TH1D* pion_abs, TH1D* pion_cex, TH1D* pion_bg, TH1D* decayPion, 
                           TH1D* muon, TH1D* proton, TH1D* cosmic, TH1D* other, string c_name, string title, string png ){
   
   double max = mc->GetBinContent(mc->GetMaximumBin());
   gStyle->SetPalette(7, beamParticle_noProt_color_abs);
   THStack *stack = new THStack("stack", "");
   stack->Add( pion_abs );
   stack->Add( pion_cex );
   stack->Add( pion_bg );
   //stack->Add( proton );
   stack->Add( muon );
   stack->Add( cosmic );
   stack->Add( decayPion );
   stack->Add( other);

   TCanvas *c_ratio = new TCanvas( c_name.c_str(), c_name.c_str(), 1200,900 );
   c_ratio->SetGrid(1,1);
   gStyle->SetOptStat(0);
   //c_ratio->SetTicks(0,1);
   auto rp = new TRatioPlot( data, mc , "divsym");
   rp->SetH1DrawOpt("PE");
   rp->SetH2DrawOpt("HIST");
   rp->Draw("SAME");
   rp->GetUpperPad()->cd();
   stack->Draw("HIST PFC PLC PMC SAME");
   stack->SetTitle( title.c_str() );
   
   data->SetMarkerColor(kBlack);
   data->SetLineColor(kBlack);
   data->SetMarkerStyle(8);
   data->SetMarkerSize(0.8);
   data->SetLineWidth(2);
   data->Draw("PE1 SAME");
   data->GetYaxis()->SetRangeUser(0, max + 0.1*max);
   data->SetTitle( title.c_str() );
   data->SetTitleSize( 35);
   rp->GetUpperPad()->SetTitle( title.c_str() );
   //rp->GetUpperPad()->BuildLegend();
   rp->GetUpperPad()->SetGrid(1,1);
   rp->GetLowerRefXaxis()->SetNdivisions(1020);
   rp->GetLowerRefXaxis()->SetLabelSize(0.025);
   rp->GetLowYaxis()->SetNdivisions(505);
   rp->GetLowerRefGraph()->SetMinimum(0);
   rp->GetLowerRefGraph()->SetMaximum(2.5);
   rp->GetLowerRefYaxis()->SetTitle("Data / MC");
   rp->GetLowerRefGraph()->SetLineColor(kBlack);
   //rp->Update();
   //

   auto leg = new TLegend(0.5, 0.4, 0.9, 0.8);
   leg->AddEntry( data, "PDSP Data, Run58XX"); leg->AddEntry( pion_abs, "Pion absorption signal"); leg->AddEntry( pion_cex, "Pion charge exchange background");
   leg->AddEntry( pion_bg, "Pion inelastic background"); //leg->AddEntry( proton, "Proton"); 
   leg->AddEntry( muon, "Beam muon");  leg->AddEntry( cosmic, "Cosmic"); leg->AddEntry( decayPion, "Pion decay" );
   leg->AddEntry( other, "Other");

   leg->SetTextSize(0.04);

   rp->GetUpperPad()->cd();
   leg->Draw();


   c_ratio->Update();
   //c_ratio->SaveAs( png.c_str() );

};


void thesisPlot_cutByCut_endZ(const string mcFile, const string dataFile)
{
   ROOT::RDataFrame inputFrame(pionTree,mcFile);
   ROOT::RDataFrame data_inputFrame(pionTree, dataFile);

   bool doPlot = true;
   bool doCount = true;

   //Int_t pal_cutByCut[6] = { 2, 8, 9, 17, 7, 42 };
   //gStyle->SetPalette(6, beamParticle_color);

   gStyle->SetNdivisions(1020);
   //ROOT::RDataFrame data_pre(pionTree, dataPath);

   string output_name;
   //if(doMC) output_name = "unfold_wBG_mc_" + std::to_string((int) bin_size_int) + "MeV.root";
   //else 
   output_name = "thesisPlot_cutByCut.root";

   TFile *output = new TFile( output_name.c_str() , "RECREATE"); //maybe save with binning?

   //=======================  Frame Definitions ==============================================================      

   //no need to filter for only reconstructed events as with the Miss function one can take into account the not-reconstructed events
   auto frame = inputFrame
      //.Range(dataMC_events)
      .Define("class_true_pion", "if(true_primPionInel == 1) return true; else return false;")
      .Define("class_true_proton", "true_beam_PDG == 2212")
      .Define("class_true_cosmic", "isCosmic && !primaryMuon && !class_true_pion")
      .Define("class_true_beamMuon", "primaryMuon && !isCosmic && !class_true_pion")
      .Define("class_true_decayPion", "isDecay && !primaryMuon && !isCosmic && !class_true_pion")
      .Define("class_true_other", "!class_true_pion && !class_true_proton && !class_true_cosmic && !class_true_beamMuon && !class_true_decayPion")
      //.Filter("true_beam_endZ > 0")
      .Range(dataMC_events);

   auto dataFrame = data_inputFrame.Range(dataMC_events);

   //auto count_data = dataFrame.Count();
   //std::cout << "count data = " << *count_data << std::endl;

   //cout << "==================================== HISTOS =====================================" << endl;
   int nBin = 70, bin_low = -100, bin_high = 600;

   TH1D* hMC_total_pion = new TH1D("hMC_total_pion", "Pion; Reconstructed end Z [cm]; Events / 10 cm ", nBin, bin_low, bin_high);
   TH1D* hMC_total_decayPion = new TH1D("hMC_total_decayPion", "Pion decay; Reconstructed end Z [cm]; Events / 10 cm ", nBin, bin_low, bin_high);
   TH1D* hMC_total_beamMuon = new TH1D("hMC_total_beamMuon", "Beam muon; Reconstructed end Z [cm]; Events / 10 cm ", nBin, bin_low, bin_high);
   TH1D* hMC_total_proton = new TH1D("hMC_total_proton", "Proton; Reconstructed end Z [cm]; Events / 10 cm ", nBin, bin_low, bin_high);
   TH1D* hMC_total_cosmic = new TH1D("hMC_total_cosmic", "Cosmic; Reconstructed end Z [cm]; Events / 10 cm ", nBin, bin_low, bin_high);
   TH1D* hMC_total_other = new TH1D("hMC_total_other", "Other; Reconstructed end Z [cm]; Events / 10 cm ", nBin, bin_low, bin_high);

   TH1D* hMC_total = new TH1D("hMC_total", "MC total; Reconstructed end Z [cm]; Events / 10 cm ", nBin, bin_low, bin_high);
   TH1D* hData_total = new TH1D("hData_total", "Data total; Reconstructed end Z [cm]; Events / 10 cm ", nBin, bin_low, bin_high);

   TH1D* hMC_pandoraBeamType_pion = new TH1D("hMC_pandoraBeamType_pion", "Pion; Reconstructed end Z [cm]; Events / 10 cm ", nBin, bin_low, bin_high);
   TH1D* hMC_pandoraBeamType_decayPion = new TH1D("hMC_pandoraBeamType_decayPion", "Pion decay; Reconstructed end Z [cm]; Events / 10 cm ", nBin, bin_low, bin_high);
   TH1D* hMC_pandoraBeamType_beamMuon = new TH1D("hMC_pandoraBeamType_beamMuon", "Beam muon; Reconstructed end Z [cm]; Events / 10 cm ", nBin, bin_low, bin_high);
   TH1D* hMC_pandoraBeamType_proton = new TH1D("hMC_pandoraBeamType_proton", "Proton; Reconstructed end Z [cm]; Events / 10 cm ", nBin, bin_low, bin_high);
   TH1D* hMC_pandoraBeamType_cosmic = new TH1D("hMC_pandoraBeamType_cosmic", "Cosmic; Reconstructed end Z [cm]; Events / 10 cm ", nBin, bin_low, bin_high);
   TH1D* hMC_pandoraBeamType_other = new TH1D("hMC_pandoraBeamType_other", "Other; Reconstructed end Z [cm]; Events / 10 cm ", nBin, bin_low, bin_high);

   TH1D* hMC_pandoraBeamType = new TH1D("hMC_pandoraBeamType", "MC pandoraBeamType; Reconstructed end Z [cm]; Events / 10 cm ", nBin, bin_low, bin_high);
   TH1D* hData_pandoraBeamType = new TH1D("hData_pandoraBeamType", "Data pandoraBeamType; Reconstructed end Z [cm]; Events / 10 cm ", nBin, bin_low, bin_high);
   
   TH1D* hMC_beamQualityPos_pion = new TH1D("hMC_beamQualityPos_pion", "Pion; Reconstructed end Z [cm]; Events / 10 cm ", nBin, bin_low, bin_high);
   TH1D* hMC_beamQualityPos_decayPion = new TH1D("hMC_beamQualityPos_decayPion", "Pion decay; Reconstructed end Z [cm]; Events / 10 cm ", nBin, bin_low, bin_high);
   TH1D* hMC_beamQualityPos_beamMuon = new TH1D("hMC_beamQualityPos_beamMuon", "Beam muon; Reconstructed end Z [cm]; Events / 10 cm ", nBin, bin_low, bin_high);
   TH1D* hMC_beamQualityPos_proton = new TH1D("hMC_beamQualityPos_proton", "Proton; Reconstructed end Z [cm]; Events / 10 cm ", nBin, bin_low, bin_high);
   TH1D* hMC_beamQualityPos_cosmic = new TH1D("hMC_beamQualityPos_cosmic", "Cosmic; Reconstructed end Z [cm]; Events / 10 cm ", nBin, bin_low, bin_high);
   TH1D* hMC_beamQualityPos_other = new TH1D("hMC_beamQualityPos_other", "Other; Reconstructed end Z [cm]; Events / 10 cm ", nBin, bin_low, bin_high);

   TH1D* hMC_beamQualityPos = new TH1D("hMC_beamQualityPos", "MC beamQuality; Reconstructed end Z [cm]; Events / 10 cm ", nBin, bin_low, bin_high);
   TH1D* hData_beamQualityPos = new TH1D("hData_beamQualityPos", "Data beamQuality; Reconstructed end Z [cm]; Events / 10 cm ", nBin, bin_low, bin_high);

   TH1D* hMC_michelScore_pion = new TH1D("hMC_michelScore_pion", "Pion; Reconstructed end Z [cm]; Events / 10 cm ", nBin, bin_low, bin_high);
   TH1D* hMC_michelScore_decayPion = new TH1D("hMC_michelScore_decayPion", "Pion decay; Reconstructed end Z [cm]; Events / 10 cm ", nBin, bin_low, bin_high);
   TH1D* hMC_michelScore_beamMuon = new TH1D("hMC_michelScore_beamMuon", "Beam muon; Reconstructed end Z [cm]; Events / 10 cm ", nBin, bin_low, bin_high);
   TH1D* hMC_michelScore_proton = new TH1D("hMC_michelScore_proton", "Proton; Reconstructed end Z [cm]; Events / 10 cm ", nBin, bin_low, bin_high);
   TH1D* hMC_michelScore_cosmic = new TH1D("hMC_michelScore_cosmic", "Cosmic; Reconstructed end Z [cm]; Events / 10 cm ", nBin, bin_low, bin_high);
   TH1D* hMC_michelScore_other = new TH1D("hMC_michelScore_other", "Other; Reconstructed end Z [cm]; Events / 10 cm ", nBin, bin_low, bin_high);

   TH1D* hMC_michelScore = new TH1D("hMC_michelScore", "MC michelScore; Reconstructed end Z [cm]; Events / 10 cm ", nBin, bin_low, bin_high);
   TH1D* hData_michelScore = new TH1D("hData_michelScore", "Data michelScore; Reconstructed end Z [cm]; Events / 10 cm ", nBin, bin_low, bin_high);
   
   TH1D* hMC_APA3_pion = new TH1D("hMC_APA3_pion", "Pion; Reconstructed end Z [cm]; Events / 10 cm ", nBin, bin_low, bin_high);
   TH1D* hMC_APA3_decayPion = new TH1D("hMC_APA3_decayPion", "Pion decay; Reconstructed end Z [cm]; Events / 10 cm ", nBin, bin_low, bin_high);
   TH1D* hMC_APA3_beamMuon = new TH1D("hMC_APA3_beamMuon", "Beam muon; Reconstructed end Z [cm]; Events / 10 cm ", nBin, bin_low, bin_high);
   TH1D* hMC_APA3_proton = new TH1D("hMC_APA3_proton", "Proton; Reconstructed end Z [cm]; Events / 10 cm ", nBin, bin_low, bin_high);
   TH1D* hMC_APA3_cosmic = new TH1D("hMC_APA3_cosmic", "Cosmic; Reconstructed end Z [cm]; Events / 10 cm ", nBin, bin_low, bin_high);
   TH1D* hMC_APA3_other = new TH1D("hMC_APA3_other", "Other; Reconstructed end Z [cm]; Events / 10 cm ", nBin, bin_low, bin_high);

   TH1D* hMC_APA3 = new TH1D("hMC_APA3", "MC APA3; Reconstructed end Z [cm]; Events / 10 cm ", nBin, bin_low, bin_high);
   TH1D* hData_APA3 = new TH1D("hData_APA3", "Data APA3; Reconstructed end Z [cm]; Events / 10 cm ", nBin, bin_low, bin_high);

   TH1D* hMC_secondaryPi_pion_BG = new TH1D("hMC_secondaryPi_pion_BG", "Pion background; Reconstructed end Z [cm]; Events / 10 cm ", nBin, bin_low, bin_high);
   TH1D* hMC_secondaryPi_pion_abs = new TH1D("hMC_secondaryPi_pion_abs", "Pion absorption; Reconstructed end Z [cm]; Events / 10 cm ", nBin, bin_low, bin_high);
   TH1D* hMC_secondaryPi_pion_cex = new TH1D("hMC_secondaryPi_pion_cex", "Pion charge exchange; Reconstructed end Z [cm]; Events / 10 cm ", nBin, bin_low, bin_high);
   TH1D* hMC_secondaryPi_decayPion = new TH1D("hMC_secondaryPi_decayPion", "Pion decay; Reconstructed end Z [cm]; Events / 10 cm ", nBin, bin_low, bin_high);
   TH1D* hMC_secondaryPi_beamMuon = new TH1D("hMC_secondaryPi_beamMuon", "Beam muon; Reconstructed end Z [cm]; Events / 10 cm ", nBin, bin_low, bin_high);
   TH1D* hMC_secondaryPi_proton = new TH1D("hMC_secondaryPi_proton", "Proton; Reconstructed end Z [cm]; Events / 10 cm ", nBin, bin_low, bin_high);
   TH1D* hMC_secondaryPi_cosmic = new TH1D("hMC_secondaryPi_cosmic", "Cosmic; Reconstructed end Z [cm]; Events / 10 cm ", nBin, bin_low, bin_high);
   TH1D* hMC_secondaryPi_other = new TH1D("hMC_secondaryPi_other", "Other; Reconstructed end Z [cm]; Events / 10 cm ", nBin, bin_low, bin_high);

   TH1D* hMC_secondaryPi = new TH1D("hMC_secondaryPi", "MC secondaryPi; Reconstructed end Z [cm]; Events / 10 cm ", nBin, bin_low, bin_high);
   TH1D* hData_secondaryPi = new TH1D("hData_secondaryPi", "data secondaryPi; Reconstructed end Z [cm]; Events / 10 cm ", nBin, bin_low, bin_high);

   TH1D* hMC_shower_pion_BG = new TH1D("hMC_shower_pion_BG", "Pion background; Reconstructed end Z [cm]; Events / 10 cm ", nBin, bin_low, bin_high);
   TH1D* hMC_shower_pion_abs = new TH1D("hMC_shower_pion_abs", "Pion absorption; Reconstructed end Z [cm]; Events / 10 cm ", nBin, bin_low, bin_high);
   TH1D* hMC_shower_pion_cex = new TH1D("hMC_shower_pion_cex", "Pion charge exchange; Reconstructed end Z [cm]; Events / 10 cm ", nBin, bin_low, bin_high);
   TH1D* hMC_shower_decayPion = new TH1D("hMC_shower_decayPion", "Pion decay; Reconstructed end Z [cm]; Events / 10 cm ", nBin, bin_low, bin_high);
   TH1D* hMC_shower_beamMuon = new TH1D("hMC_shower_beamMuon", "Beam muon; Reconstructed end Z [cm]; Events / 10 cm ", nBin, bin_low, bin_high);
   TH1D* hMC_shower_proton = new TH1D("hMC_shower_proton", "Proton; Reconstructed end Z [cm]; Events / 10 cm ", nBin, bin_low, bin_high);
   TH1D* hMC_shower_cosmic = new TH1D("hMC_shower_cosmic", "Cosmic; Reconstructed end Z [cm]; Events / 10 cm ", nBin, bin_low, bin_high);
   TH1D* hMC_shower_other = new TH1D("hMC_shower_other", "Other; Reconstructed end Z [cm]; Events / 10 cm ", nBin, bin_low, bin_high);

   TH1D* hMC_shower = new TH1D("hMC_shower", "MC shower; Reconstructed end Z [cm]; Events / 10 cm ", nBin, bin_low, bin_high);
   TH1D* hData_shower = new TH1D("hData_shower", "Data shower; Reconstructed end Z [cm]; Events / 10 cm ", nBin, bin_low, bin_high);


   //===============================================================================================================
   //          TOTAL
   //===============================================================================================================
   
   auto mc_total = frame;

   mc_total
      .Foreach( [hMC_total, hMC_total_pion, hMC_total_decayPion, hMC_total_beamMuon, hMC_total_proton, hMC_total_cosmic, hMC_total_other]
                (double endZ, 
                 bool pion, bool proton, bool beamMuon, bool cosmic, bool decayPion, bool other)
                {
                  hMC_total->Fill( endZ );
                  if( pion ) hMC_total_pion->Fill( endZ );
                  else if( proton ) hMC_total_proton->Fill( endZ );
                  else if( beamMuon ) hMC_total_beamMuon->Fill( endZ );
                  else if( cosmic ) hMC_total_cosmic->Fill( endZ );
                  else if( decayPion ) hMC_total_decayPion->Fill( endZ );
                  else if( other ) hMC_total_other->Fill( endZ );

                },
                {"reco_beam_endZ", 
                 "class_true_pion", "class_true_proton", "class_true_beamMuon", "class_true_cosmic", "class_true_decayPion", "class_true_other"}
            );

   auto data_total = dataFrame;

   data_total
      .Foreach( [ hData_total ] (double endZ)
            {
               hData_total->Fill( endZ );
            },
            {"reco_beam_endZ"}
            );

   //===============================================================================================================
   //          pandoraBeamType
   //===============================================================================================================
   
   auto mc_pandoraBeamType = mc_total.Filter("primary_isBeamType");

   mc_pandoraBeamType
      .Foreach( [hMC_pandoraBeamType, hMC_pandoraBeamType_pion, hMC_pandoraBeamType_decayPion, hMC_pandoraBeamType_beamMuon, 
                 hMC_pandoraBeamType_proton, hMC_pandoraBeamType_cosmic, hMC_pandoraBeamType_other]
                (double endZ, 
                 bool pion, bool proton, bool beamMuon, bool cosmic, bool decayPion, bool other)
                {
                  hMC_pandoraBeamType->Fill( endZ );
                  if( pion ) hMC_pandoraBeamType_pion->Fill( endZ );
                  else if( proton ) hMC_pandoraBeamType_proton->Fill( endZ );
                  else if( beamMuon ) hMC_pandoraBeamType_beamMuon->Fill( endZ );
                  else if( cosmic ) hMC_pandoraBeamType_cosmic->Fill( endZ );
                  else if( decayPion ) hMC_pandoraBeamType_decayPion->Fill( endZ );
                  else if( other ) hMC_pandoraBeamType_other->Fill( endZ );

                },
                {"reco_beam_endZ", 
                 "class_true_pion", "class_true_proton", "class_true_beamMuon", "class_true_cosmic", "class_true_decayPion", "class_true_other"}
            );

   auto data_pandoraBeamType = data_total.Filter("primary_isBeamType");

   data_pandoraBeamType
      .Foreach( [ hData_pandoraBeamType ] (double endZ)
            {
               hData_pandoraBeamType->Fill( endZ );
            },
            {"reco_beam_endZ"}
            );

   //===============================================================================================================
   //          beamQualityPos
   //===============================================================================================================
   
   auto mc_beamQualityPos = mc_pandoraBeamType.Filter("passBeamQuality_TPCjustPosition");

   mc_beamQualityPos
      .Foreach( [hMC_beamQualityPos, hMC_beamQualityPos_pion, hMC_beamQualityPos_decayPion, hMC_beamQualityPos_beamMuon, 
                 hMC_beamQualityPos_proton, hMC_beamQualityPos_cosmic, hMC_beamQualityPos_other]
                (double endZ, 
                 bool pion, bool proton, bool beamMuon, bool cosmic, bool decayPion, bool other)
                {
                  hMC_beamQualityPos->Fill( endZ );
                  if( pion ) hMC_beamQualityPos_pion->Fill( endZ );
                  else if( proton ) hMC_beamQualityPos_proton->Fill( endZ );
                  else if( beamMuon ) hMC_beamQualityPos_beamMuon->Fill( endZ );
                  else if( cosmic ) hMC_beamQualityPos_cosmic->Fill( endZ );
                  else if( decayPion ) hMC_beamQualityPos_decayPion->Fill( endZ );
                  else if( other ) hMC_beamQualityPos_other->Fill( endZ );

                },
                {"reco_beam_endZ", 
                 "class_true_pion", "class_true_proton", "class_true_beamMuon", "class_true_cosmic", "class_true_decayPion", "class_true_other"}
            );

   auto data_beamQualityPos = data_pandoraBeamType.Filter("passBeamQuality_TPCjustPosition");

   data_beamQualityPos
      .Foreach( [ hData_beamQualityPos ] (double endZ)
            {
               hData_beamQualityPos->Fill( endZ );
            },
            {"reco_beam_endZ"}
            );

   //===============================================================================================================
   //          APA3
   //===============================================================================================================
   
   auto mc_APA3 = mc_beamQualityPos.Filter("primary_ends_inAPA3");

   mc_APA3
      .Foreach( [hMC_APA3, hMC_APA3_pion, hMC_APA3_decayPion, hMC_APA3_beamMuon, 
                 hMC_APA3_proton, hMC_APA3_cosmic, hMC_APA3_other]
                (double endZ, 
                 bool pion, bool proton, bool beamMuon, bool cosmic, bool decayPion, bool other)
                {
                  hMC_APA3->Fill( endZ );
                  if( pion ) hMC_APA3_pion->Fill( endZ );
                  else if( proton ) hMC_APA3_proton->Fill( endZ );
                  else if( beamMuon ) hMC_APA3_beamMuon->Fill( endZ );
                  else if( cosmic ) hMC_APA3_cosmic->Fill( endZ );
                  else if( decayPion ) hMC_APA3_decayPion->Fill( endZ );
                  else if( other ) hMC_APA3_other->Fill( endZ );

                },
                {"reco_beam_endZ", 
                 "class_true_pion", "class_true_proton", "class_true_beamMuon", "class_true_cosmic", "class_true_decayPion", "class_true_other"}
            );
   auto data_APA3 = data_beamQualityPos.Filter("primary_ends_inAPA3");

   data_APA3
      .Foreach( [ hData_APA3 ] (double endZ)
            {
               hData_APA3->Fill( endZ );
            },
            {"reco_beam_endZ"}
            );

   //===============================================================================================================
   //          michelScore
   //===============================================================================================================
   
   auto mc_michelScore = mc_APA3.Filter("!isPrimaryMuonCandidate");

   mc_michelScore
      .Foreach( [hMC_michelScore, hMC_michelScore_pion, hMC_michelScore_decayPion, hMC_michelScore_beamMuon, 
                 hMC_michelScore_proton, hMC_michelScore_cosmic, hMC_michelScore_other]
                (double endZ, 
                 bool pion, bool proton, bool beamMuon, bool cosmic, bool decayPion, bool other)
                {
                  hMC_michelScore->Fill( endZ );
                  if( pion ) hMC_michelScore_pion->Fill( endZ );
                  else if( proton ) hMC_michelScore_proton->Fill( endZ );
                  else if( beamMuon ) hMC_michelScore_beamMuon->Fill( endZ );
                  else if( cosmic ) hMC_michelScore_cosmic->Fill( endZ );
                  else if( decayPion ) hMC_michelScore_decayPion->Fill( endZ );
                  else if( other ) hMC_michelScore_other->Fill( endZ );

                },
                {"reco_beam_endZ", 
                 "class_true_pion", "class_true_proton", "class_true_beamMuon", "class_true_cosmic", "class_true_decayPion", "class_true_other"}
            );
   
   auto data_michelScore = data_APA3.Filter("!isPrimaryMuonCandidate");

   data_michelScore
      .Foreach( [ hData_michelScore ] (double endZ)
            {
               hData_michelScore->Fill( endZ );
            },
            {"reco_beam_endZ"}
            );


   //===============================================================================================================
   //          secondaryPi
   //===============================================================================================================
   
   auto mc_secondaryPi = mc_michelScore.Filter("has_noPion_daughter");

   mc_secondaryPi
      .Foreach( [hMC_secondaryPi, hMC_secondaryPi_pion_abs, hMC_secondaryPi_pion_cex, hMC_secondaryPi_pion_BG, 
                 hMC_secondaryPi_decayPion, hMC_secondaryPi_beamMuon, hMC_secondaryPi_proton, 
                 hMC_secondaryPi_cosmic, hMC_secondaryPi_other]
                (double endZ, int true_abs, int true_cex,
                 bool pion, bool proton, bool beamMuon, bool cosmic, bool decayPion, bool other)
                {
                  hMC_secondaryPi->Fill( endZ );
                  if( pion && true_abs ) hMC_secondaryPi_pion_abs->Fill( endZ );
                  if( pion && true_cex ) hMC_secondaryPi_pion_cex->Fill( endZ );
                  else if( pion && !true_abs && !true_cex ) hMC_secondaryPi_pion_BG->Fill( endZ );
                  else if( proton ) hMC_secondaryPi_proton->Fill( endZ );
                  else if( beamMuon ) hMC_secondaryPi_beamMuon->Fill( endZ );
                  else if( cosmic ) hMC_secondaryPi_cosmic->Fill( endZ );
                  else if( decayPion ) hMC_secondaryPi_decayPion->Fill( endZ );
                  else if( other ) hMC_secondaryPi_other->Fill( endZ );

                },
                {"reco_beam_endZ", "true_absSignal", "true_chexSignal",
                 "class_true_pion", "class_true_proton", "class_true_beamMuon", "class_true_cosmic", "class_true_decayPion", "class_true_other"}
            );
   
   auto data_secondaryPi = data_michelScore.Filter("has_noPion_daughter");

   data_secondaryPi
      .Foreach( [ hData_secondaryPi ] (double endZ)
            {
               hData_secondaryPi->Fill( endZ );
            },
            {"reco_beam_endZ"}
            );


   //===============================================================================================================
   //          shower
   //===============================================================================================================
   
   auto mc_shower = mc_secondaryPi.Filter("!has_shower_nHits_distance");

   mc_shower
      .Foreach( [hMC_shower, hMC_shower_pion_abs, hMC_shower_pion_cex, hMC_shower_pion_BG, 
                 hMC_shower_decayPion, hMC_shower_beamMuon, hMC_shower_proton, 
                 hMC_shower_cosmic, hMC_shower_other]
                (double endZ, int true_abs, int true_cex,
                 bool pion, bool proton, bool beamMuon, bool cosmic, bool decayPion, bool other)
                {
                  hMC_shower->Fill( endZ );
                  if( pion && true_abs ) hMC_shower_pion_abs->Fill( endZ );
                  if( pion && true_cex ) hMC_shower_pion_cex->Fill( endZ );
                  else if( pion && !true_abs && !true_cex ) hMC_shower_pion_BG->Fill( endZ );
                  else if( proton ) hMC_shower_proton->Fill( endZ );
                  else if( beamMuon ) hMC_shower_beamMuon->Fill( endZ );
                  else if( cosmic ) hMC_shower_cosmic->Fill( endZ );
                  else if( decayPion ) hMC_shower_decayPion->Fill( endZ );
                  else if( other ) hMC_shower_other->Fill( endZ );

                },
                {"reco_beam_endZ", "true_absSignal", "true_chexSignal",
                 "class_true_pion", "class_true_proton", "class_true_beamMuon", "class_true_cosmic", "class_true_decayPion", "class_true_other"}
            );

   auto data_shower = data_secondaryPi.Filter("!has_shower_nHits_distance");

   data_shower
      .Foreach( [ hData_shower ] (double endZ)
            {
               hData_shower->Fill( endZ );
            },
            {"reco_beam_endZ"}
            );

   //===============================================================================================================
   //          Scale MC to Data ONE FACTOR for all otherwise not COMPARABLE!!! --> Factor from Total!
   //===============================================================================================================
   double scale_total = hData_total->Integral() / hMC_total->Integral();
   hMC_total->Scale( scale_total );
   hMC_total->Sumw2(0);
   hMC_total_pion->Scale( scale_total ); hMC_total_decayPion->Scale( scale_total ); hMC_total_beamMuon->Scale( scale_total );
   hMC_total_proton->Scale( scale_total ); hMC_total_cosmic->Scale( scale_total ); hMC_total_other->Scale( scale_total );
   
   hMC_pandoraBeamType->Scale( scale_total );
   hMC_pandoraBeamType->Sumw2(0);
   hMC_pandoraBeamType_pion->Scale( scale_total ); hMC_pandoraBeamType_decayPion->Scale( scale_total ); hMC_pandoraBeamType_beamMuon->Scale( scale_total );
   hMC_pandoraBeamType_proton->Scale( scale_total ); hMC_pandoraBeamType_cosmic->Scale( scale_total ); hMC_pandoraBeamType_other->Scale( scale_total );

   hMC_beamQualityPos->Scale(scale_total );
   hMC_beamQualityPos->Sumw2(0);
   hMC_beamQualityPos_pion->Scale( scale_total ); hMC_beamQualityPos_decayPion->Scale( scale_total ); hMC_beamQualityPos_beamMuon->Scale( scale_total );
   hMC_beamQualityPos_proton->Scale( scale_total ); hMC_beamQualityPos_cosmic->Scale( scale_total ); hMC_beamQualityPos_other->Scale( scale_total );
   
   hMC_michelScore->Scale( scale_total );
   hMC_michelScore->Sumw2(0);
   hMC_michelScore_pion->Scale( scale_total ); hMC_michelScore_decayPion->Scale( scale_total ); hMC_michelScore_beamMuon->Scale( scale_total );
   hMC_michelScore_proton->Scale( scale_total ); hMC_michelScore_cosmic->Scale( scale_total ); hMC_michelScore_other->Scale( scale_total );
   
   hMC_APA3->Scale( scale_total );
   hMC_APA3->Sumw2(0);
   hMC_APA3_pion->Scale( scale_total ); hMC_APA3_decayPion->Scale( scale_total ); hMC_APA3_beamMuon->Scale( scale_total );
   hMC_APA3_proton->Scale( scale_total ); hMC_APA3_cosmic->Scale( scale_total ); hMC_APA3_other->Scale( scale_total );
    
   hMC_secondaryPi->Scale( scale_total ); 
   hMC_secondaryPi->Sumw2(0);
   hMC_secondaryPi_pion_abs->Scale( scale_total ); hMC_secondaryPi_pion_cex->Scale( scale_total ); hMC_secondaryPi_pion_BG->Scale( scale_total ); 
   hMC_secondaryPi_decayPion->Scale( scale_total ); hMC_secondaryPi_beamMuon->Scale( scale_total ); hMC_secondaryPi_proton->Scale( scale_total ); 
   hMC_secondaryPi_cosmic->Scale( scale_total ); hMC_secondaryPi_other->Scale( scale_total );
   
   hMC_shower->Scale( scale_total );
   hMC_shower->Sumw2(0);
   hMC_shower_pion_abs->Scale( scale_total ); hMC_shower_pion_cex->Scale( scale_total ); hMC_shower_pion_BG->Scale( scale_total ); 
   hMC_shower_decayPion->Scale( scale_total ); hMC_shower_beamMuon->Scale( scale_total ); hMC_shower_proton->Scale( scale_total ); 
   hMC_shower_cosmic->Scale( scale_total ); hMC_shower_other->Scale( scale_total );
   

   //===============================================================================================================
   //         Plotting 
   //===============================================================================================================
   //
   string s_total = "c_ratio_total"; 
   string s_pandoraBeamType = "c_ratio_pandoraBeamType";
   string s_beamQualityPos = "c_ratio_beamQualityPos";
   string s_APA3 = "c_ratio_APA3";
   string s_michelScore = "c_ratio_michelScore";
   string s_secondaryPi = "c_ratio_secondaryPi";
   string s_shower = "c_ratio_shower";

   string title_total = "Total"; 
   string title_pandoraBeamType = "Pandora Beam Type";
   string title_beamQualityPos = "Beam Position";
   string title_APA3 = "APA3";
   string title_michelScore = "Michel Score";
   string title_secondaryPi = "Reject Secondary Pion Candidates";
   string title_shower = "Reject Secondary Shower Candidates";

   string pdf_total = "thesisPlot/cutByCut_endZ_total.pdf"; 
   string pdf_pandoraBeamType = "thesisPlot/cutByCut_endZ_pandoraBeamType.pdf";
   string pdf_beamQualityPos = "thesisPlot/cutByCut_endZ_beamQualityPos.pdf";
   string pdf_APA3 = "thesisPlot/cutByCut_endZ_APA3.pdf";
   string pdf_michelScore = "thesisPlot/cutByCut_endZ_michelScore.pdf";
   string pdf_secondaryPi = "thesisPlot/cutByCut_endZ_secondaryPi.pdf";
   string pdf_shower = "thesisPlot/cutByCut_endZ_shower.pdf";


   if(doPlot){
   doRatioPlot_pion( hData_total, hMC_total, hMC_total_pion, hMC_total_decayPion, hMC_total_beamMuon, 
         hMC_total_proton, hMC_total_cosmic, hMC_total_other, s_total, title_total, pdf_total);
   doRatioPlot_pion( hData_pandoraBeamType, hMC_pandoraBeamType, hMC_pandoraBeamType_pion, hMC_pandoraBeamType_decayPion, hMC_pandoraBeamType_beamMuon, 
         hMC_pandoraBeamType_proton, hMC_pandoraBeamType_cosmic, hMC_pandoraBeamType_other, s_pandoraBeamType, title_pandoraBeamType, pdf_pandoraBeamType);
   doRatioPlot_pion( hData_beamQualityPos, hMC_beamQualityPos, hMC_beamQualityPos_pion, hMC_beamQualityPos_decayPion, hMC_beamQualityPos_beamMuon, 
         hMC_beamQualityPos_proton, hMC_beamQualityPos_cosmic, hMC_beamQualityPos_other, s_beamQualityPos, title_beamQualityPos, pdf_beamQualityPos);
   doRatioPlot_pion( hData_APA3, hMC_APA3, hMC_APA3_pion, hMC_APA3_decayPion, hMC_APA3_beamMuon, 
         hMC_APA3_proton, hMC_APA3_cosmic, hMC_APA3_other, s_APA3, title_APA3, pdf_APA3);
   doRatioPlot_pion( hData_michelScore, hMC_michelScore, hMC_michelScore_pion, hMC_michelScore_decayPion, hMC_michelScore_beamMuon, 
         hMC_michelScore_proton, hMC_michelScore_cosmic, hMC_michelScore_other, s_michelScore, title_michelScore, pdf_michelScore);
   
//   doRatioPlot_pion_abs( hData_secondaryPi, hMC_secondaryPi, hMC_secondaryPi_pion_abs, hMC_secondaryPi_pion_cex,
//         hMC_secondaryPi_pion_BG, hMC_secondaryPi_decayPion, hMC_secondaryPi_beamMuon, hMC_secondaryPi_proton, hMC_secondaryPi_cosmic, 
//         hMC_secondaryPi_other, s_secondaryPi, title_secondaryPi, pdf_secondaryPi);
//   
//   doRatioPlot_pion_abs( hData_shower, hMC_shower, hMC_shower_pion_abs, hMC_shower_pion_cex,
//         hMC_shower_pion_BG, hMC_shower_decayPion, hMC_shower_beamMuon, hMC_shower_proton, hMC_shower_cosmic, 
//         hMC_shower_other, s_shower, title_shower, pdf_shower);

   }

   if(doCount){

      std::cout << "----------------- Counting ----------------------------" << std::endl;
      std::cout <<"" <<std::endl;
      std::cout <<"TOTAL" << std::endl;
      std::cout << "Data          = " << *data_total.Count() << std::endl;
      std::cout << "MC            = " << *mc_total.Count() << std::endl;
      std::cout << "MC Total Pion = " << *mc_total.Filter("class_true_pion").Count() << std::endl;
      std::cout << "MC True Abs   = " << *mc_total.Filter("class_true_pion && true_absSignal").Count() << std::endl;
      std::cout << "MC True Cex   = " << *mc_total.Filter("class_true_pion && true_chexSignal").Count() << std::endl;
      std::cout << "MC True Inel  = " << *mc_total.Filter("class_true_pion && !true_absSignal && !true_chexSignal").Count() << std::endl;
      std::cout << "MC Total BG   = " << *mc_total.Filter("!class_true_pion").Count() << std::endl;
      std::cout << "MC Proton     = " << *mc_total.Filter("class_true_proton").Count() << std::endl;
      std::cout << "MC BeamMuon   = " << *mc_total.Filter("class_true_beamMuon").Count() << std::endl;
      std::cout << "MC Cosmic     = " << *mc_total.Filter("class_true_cosmic").Count() << std::endl;
      std::cout << "MC DecayPion  = " << *mc_total.Filter("class_true_decayPion").Count() << std::endl;
      std::cout << "MC Other      = " << *mc_total.Filter("class_true_other").Count() << std::endl;
      std::cout << "----------------------------------------------------------------------------" << std::endl;
      std::cout <<"Pandora" << std::endl;
      std::cout << "Data          = " << *data_pandoraBeamType.Count() << std::endl;
      std::cout << "MC            = " << *mc_pandoraBeamType.Count() << std::endl;
      std::cout << "MC Total Pion = " << *mc_pandoraBeamType.Filter("class_true_pion").Count() << std::endl;
      std::cout << "MC True Abs   = " << *mc_pandoraBeamType.Filter("class_true_pion && true_absSignal").Count() << std::endl;
      std::cout << "MC True Cex   = " << *mc_pandoraBeamType.Filter("class_true_pion && true_chexSignal").Count() << std::endl;
      std::cout << "MC True Inel  = " << *mc_pandoraBeamType.Filter("class_true_pion && !true_absSignal && !true_chexSignal").Count() << std::endl;
      std::cout << "MC Total BG   = " << *mc_pandoraBeamType.Filter("!class_true_pion").Count() << std::endl;
      std::cout << "MC Proton     = " << *mc_pandoraBeamType.Filter("class_true_proton").Count() << std::endl;
      std::cout << "MC BeamMuon   = " << *mc_pandoraBeamType.Filter("class_true_beamMuon").Count() << std::endl;
      std::cout << "MC Cosmic     = " << *mc_pandoraBeamType.Filter("class_true_cosmic").Count() << std::endl;
      std::cout << "MC DecayPion  = " << *mc_pandoraBeamType.Filter("class_true_decayPion").Count() << std::endl;
      std::cout << "MC Other      = " << *mc_pandoraBeamType.Filter("class_true_other").Count() << std::endl;
      std::cout << "----------------------------------------------------------------------------" << std::endl;
      std::cout <<"BEAM POSITION" << std::endl;
      std::cout << "Data          = " << *data_beamQualityPos.Count() << std::endl;
      std::cout << "MC            = " << *mc_beamQualityPos.Count() << std::endl;
      std::cout << "MC Total Pion = " << *mc_beamQualityPos.Filter("class_true_pion").Count() << std::endl;
      std::cout << "MC True Abs   = " << *mc_beamQualityPos.Filter("class_true_pion && true_absSignal").Count() << std::endl;
      std::cout << "MC True Cex   = " << *mc_beamQualityPos.Filter("class_true_pion && true_chexSignal").Count() << std::endl;
      std::cout << "MC True Inel  = " << *mc_beamQualityPos.Filter("class_true_pion && !true_absSignal && !true_chexSignal").Count() << std::endl;
      std::cout << "MC Total BG   = " << *mc_beamQualityPos.Filter("!class_true_pion").Count() << std::endl;
      std::cout << "MC Proton     = " << *mc_beamQualityPos.Filter("class_true_proton").Count() << std::endl;
      std::cout << "MC BeamMuon   = " << *mc_beamQualityPos.Filter("class_true_beamMuon").Count() << std::endl;
      std::cout << "MC Cosmic     = " << *mc_beamQualityPos.Filter("class_true_cosmic").Count() << std::endl;
      std::cout << "MC DecayPion  = " << *mc_beamQualityPos.Filter("class_true_decayPion").Count() << std::endl;
      std::cout << "MC Other      = " << *mc_beamQualityPos.Filter("class_true_other").Count() << std::endl;
      std::cout << "----------------------------------------------------------------------------" << std::endl;
      std::cout <<"APA 3" << std::endl;
      std::cout << "Data          = " << *data_APA3.Count() << std::endl;
      std::cout << "MC            = " << *mc_APA3.Count() << std::endl;
      std::cout << "MC Total Pion = " << *mc_APA3.Filter("class_true_pion").Count() << std::endl;
      std::cout << "MC True Abs   = " << *mc_APA3.Filter("class_true_pion && true_absSignal").Count() << std::endl;
      std::cout << "MC True Cex   = " << *mc_APA3.Filter("class_true_pion && true_chexSignal").Count() << std::endl;
      std::cout << "MC True Inel  = " << *mc_APA3.Filter("class_true_pion && !true_absSignal && !true_chexSignal").Count() << std::endl;
      std::cout << "MC Total BG   = " << *mc_APA3.Filter("!class_true_pion").Count() << std::endl;
      std::cout << "MC Proton     = " << *mc_APA3.Filter("class_true_proton").Count() << std::endl;
      std::cout << "MC BeamMuon   = " << *mc_APA3.Filter("class_true_beamMuon").Count() << std::endl;
      std::cout << "MC Cosmic     = " << *mc_APA3.Filter("class_true_cosmic").Count() << std::endl;
      std::cout << "MC DecayPion  = " << *mc_APA3.Filter("class_true_decayPion").Count() << std::endl;
      std::cout << "MC Other      = " << *mc_APA3.Filter("class_true_other").Count() << std::endl;
      std::cout << "----------------------------------------------------------------------------" << std::endl;
      std::cout <<"MICHEL SCORE" << std::endl;
      std::cout << "Data          = " << *data_michelScore.Count() << std::endl;
      std::cout << "MC            = " << *mc_michelScore.Count() << std::endl;
      std::cout << "MC Total Pion = " << *mc_michelScore.Filter("class_true_pion").Count() << std::endl;
      std::cout << "MC True Abs   = " << *mc_michelScore.Filter("class_true_pion && true_absSignal").Count() << std::endl;
      std::cout << "MC True Cex   = " << *mc_michelScore.Filter("class_true_pion && true_chexSignal").Count() << std::endl;
      std::cout << "MC True Inel  = " << *mc_michelScore.Filter("class_true_pion && !true_absSignal && !true_chexSignal").Count() << std::endl;
      std::cout << "MC Total BG   = " << *mc_michelScore.Filter("!class_true_pion").Count() << std::endl;
      std::cout << "MC Proton     = " << *mc_michelScore.Filter("class_true_proton").Count() << std::endl;
      std::cout << "MC BeamMuon   = " << *mc_michelScore.Filter("class_true_beamMuon").Count() << std::endl;
      std::cout << "MC Cosmic     = " << *mc_michelScore.Filter("class_true_cosmic").Count() << std::endl;
      std::cout << "MC DecayPion  = " << *mc_michelScore.Filter("class_true_decayPion").Count() << std::endl;
      std::cout << "MC Other      = " << *mc_michelScore.Filter("class_true_other").Count() << std::endl;
      std::cout << "----------------------------------------------------------------------------" << std::endl;
      std::cout <<"Reject Secondary Pion" << std::endl;
      std::cout << "Data          = " << *data_secondaryPi.Count() << std::endl;
      std::cout << "MC            = " << *mc_secondaryPi.Count() << std::endl;
      std::cout << "MC Total Pion = " << *mc_secondaryPi.Filter("class_true_pion").Count() << std::endl;
      std::cout << "MC True Abs   = " << *mc_secondaryPi.Filter("class_true_pion && true_absSignal").Count() << std::endl;
      std::cout << "MC True Cex   = " << *mc_secondaryPi.Filter("class_true_pion && true_chexSignal").Count() << std::endl;
      std::cout << "MC True Inel  = " << *mc_secondaryPi.Filter("class_true_pion && !true_absSignal && !true_chexSignal").Count() << std::endl;
      std::cout << "MC Total BG   = " << *mc_secondaryPi.Filter("!class_true_pion").Count() << std::endl;
      std::cout << "MC Proton     = " << *mc_secondaryPi.Filter("class_true_proton").Count() << std::endl;
      std::cout << "MC BeamMuon   = " << *mc_secondaryPi.Filter("class_true_beamMuon").Count() << std::endl;
      std::cout << "MC Cosmic     = " << *mc_secondaryPi.Filter("class_true_cosmic").Count() << std::endl;
      std::cout << "MC DecayPion  = " << *mc_secondaryPi.Filter("class_true_decayPion").Count() << std::endl;
      std::cout << "MC Other      = " << *mc_secondaryPi.Filter("class_true_other").Count() << std::endl;
      std::cout << "----------------------------------------------------------------------------" << std::endl;
      std::cout <<"Reject Secondary Shower" << std::endl;
      std::cout << "Data          = " << *data_shower.Count() << std::endl;
      std::cout << "MC            = " << *mc_shower.Count() << std::endl;
      std::cout << "MC Total Pion = " << *mc_shower.Filter("class_true_pion").Count() << std::endl;
      std::cout << "MC True Abs   = " << *mc_shower.Filter("class_true_pion && true_absSignal").Count() << std::endl;
      std::cout << "MC True Cex   = " << *mc_shower.Filter("class_true_pion && true_chexSignal").Count() << std::endl;
      std::cout << "MC True Inel  = " << *mc_shower.Filter("class_true_pion && !true_absSignal && !true_chexSignal").Count() << std::endl;
      std::cout << "MC Total BG   = " << *mc_shower.Filter("!class_true_pion").Count() << std::endl;
      std::cout << "MC Proton     = " << *mc_shower.Filter("class_true_proton").Count() << std::endl;
      std::cout << "MC BeamMuon   = " << *mc_shower.Filter("class_true_beamMuon").Count() << std::endl;
      std::cout << "MC Cosmic     = " << *mc_shower.Filter("class_true_cosmic").Count() << std::endl;
      std::cout << "MC DecayPion  = " << *mc_shower.Filter("class_true_decayPion").Count() << std::endl;
      std::cout << "MC Other      = " << *mc_shower.Filter("class_true_other").Count() << std::endl;
      std::cout << "----------------------------------------------------------------------------"<< std::endl;




   }
}

#ifndef __CINT__
int main () { thesisPlot_cutByCut_endZ(); return 0; }  // Main program when run stand-alone
#endif


