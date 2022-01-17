#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
#include <iostream>
using std::cout;
using std::endl;
using namespace std;
using namespace ROOT::VecOps;

#include "TRandom.h"
#include "TH1D.h"
#include "TCanvas.h"

#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#include "RooUnfold.h"
//#include "RooUnfoldSvd.h"
//#include "RooUnfoldTUnfold.h"
//#include "RooUnfoldIds.h"
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

//==============================================================================
// Global definitions
//==============================================================================

const Double_t cutdummy= -99999.0;
const Double_t dataMC_events = 68500; //root file DATA has 70121 ev, root file MC has 74149 ev


//==============================================================================
// Example Unfolding
//==============================================================================

void thesisPlot_meas_intIncHisto(const string mcFile, const string dataFile)
{
   ROOT::RDataFrame inputFrame(pionTree,mcFile);
   ROOT::RDataFrame data_inputFrame(pionTree, dataFile);

   gStyle->SetNdivisions(1020);
   //ROOT::RDataFrame data_pre(pionTree, dataPath);

   string output_name;
   //if(doMC) output_name = "unfold_wBG_mc_" + std::to_string((int) bin_size_int) + "MeV.root";
   //else 
   output_name = "thesisPlot_intInc_meas" + std::to_string((int) bin_size_int) + "MeV.root";

   TFile *output = new TFile( output_name.c_str() , "RECREATE"); //maybe save with binning?

   //=======================  Frame Definitions ==============================================================      

   //no need to filter for only reconstructed events as with the Miss function one can take into account the not-reconstructed events
   auto frame = inputFrame
      .Range(dataMC_events)
      .Define("reco_equalBin", equalBin, {"reco_initKE", "reco_interKE"});
   //Definitions for Response Matrix generation
   auto dataFrame = data_inputFrame
      .Range(dataMC_events)
      .Define("reco_equalBin", equalBin, {"reco_initKE", "reco_interKE"});

   //=====================================================================================   
   int lowE = 0, highE = 1200, nBin_new = highE / 50; 

   TH1D* hMeas_mc_initE= new TH1D ("meas_mc_initE", "Measured MC initE; Kinetic energy [MeV]; Events / 50 MeV", nBin_new, lowE, highE);
   TH1D* hMeas_data_initE= new TH1D ("meas_data_initE", "Measured data initE; Kinetic energy [MeV]; Events / 50 MeV", nBin_new, lowE, highE);

   TH1D* hMeas_mc_interE= new TH1D ("meas_mc_interE", "Measured MC interE; Kinetic energy [MeV]; Events / 50 MeV", nBin_new, lowE, highE);
   TH1D* hMeas_data_interE= new TH1D ("meas_data_interE", "Measured data interE; Kinetic energy [MeV]; Events / 50 MeV", nBin_new, lowE, highE);

   TH1D* hMeas_mc_totInel= new TH1D ("meas_mc_totInel", "Measured MC totInel; Kinetic energy [MeV]; Events / 50 MeV", nBin_new, lowE, highE);
   TH1D* hMeas_data_totInel= new TH1D ("meas_data_totInel", "Measured data totInel; Kinetic energy [MeV]; Events / 50 MeV", nBin_new, lowE, highE);

   TH1D* hMeas_mc_abs= new TH1D ("meas_mc_abs", "Measured MC abs; Kinetic energy [MeV]; Events / 50 MeV", nBin_new, lowE, highE);
   TH1D* hMeas_data_abs= new TH1D ("meas_data_abs", "Measured data abs; Kinetic energy [MeV]; Events / 50 MeV", nBin_new, lowE, highE);


   cout << "==================================== Building the Measured Distribution =====================================" << endl;

   frame

      .Filter("selected_incidentPion")
      .Foreach( [ hMeas_mc_initE, hMeas_mc_interE, hMeas_mc_totInel ]( double reco_initE, double reco_interE, bool recoEqualBin){

            //Measured
            if(reco_initE != -999. && reco_interE != -999. && !recoEqualBin){

            hMeas_mc_initE->Fill( reco_initE );
            hMeas_mc_interE->Fill( reco_interE );
            hMeas_mc_totInel->Fill( reco_interE );
            }


            }
            ,{"reco_initKE_rwData", "reco_interKE_rwData", "reco_equalBin"}
            //,{"reco_initKE", "reco_interKE", "reco_equalBin"}
            );

   frame

      .Filter("selected_abs")
      .Foreach( [ hMeas_mc_abs ]( double reco_initE, double reco_interE, bool recoEqualBin){

            //Measured
            if(reco_initE != -999. && reco_interE != -999. && !recoEqualBin){
            hMeas_mc_abs->Fill( reco_interE );
            }

            }
            ,{"reco_initKE_rwData", "reco_interKE_rwData", "reco_equalBin"}
            //,{"reco_initKE", "reco_interKE", "reco_equalBin"}
            );

   //====================  DATA ==========================
   dataFrame
      .Filter("selected_incidentPion")
      .Foreach( [ hMeas_data_initE, hMeas_data_interE, hMeas_data_totInel ]( double reco_initE, double reco_interE, bool recoEqualBin){

            //Measured
            if(reco_initE != -999. && reco_interE != -999. && !recoEqualBin){

            hMeas_data_initE->Fill( reco_initE );
            hMeas_data_interE->Fill( reco_interE );
            hMeas_data_totInel->Fill( reco_interE );
            }

            }
            ,{"reco_initKE", "reco_interKE", "reco_equalBin"}
            );

   dataFrame
      .Filter("selected_abs")
      .Foreach( [ hMeas_data_abs ]( double reco_initE, double reco_interE, bool recoEqualBin){

            //Measured
            if(reco_initE != -999. && reco_interE != -999. && !recoEqualBin){
            hMeas_data_abs->Fill( reco_interE );
            }

            }
            ,{ "reco_initKE", "reco_interKE", "reco_equalBin"}
            );

   hMeas_mc_initE->Write();
   hMeas_mc_interE->Write();
   hMeas_mc_totInel->Write();
   hMeas_mc_abs->Write();

   hMeas_data_initE->Write();
   hMeas_data_interE->Write();
   hMeas_data_totInel->Write();
   hMeas_data_abs->Write();

   double scale_initE = hMeas_data_initE->Integral() / hMeas_mc_initE->Integral();
   double scale_interE = hMeas_data_interE->Integral() / hMeas_mc_interE->Integral();
   double scale_totInel = hMeas_data_totInel->Integral() / hMeas_mc_totInel->Integral();
   double scale_abs = hMeas_data_abs->Integral() / hMeas_mc_abs->Integral();

   hMeas_mc_initE->Scale( scale_initE );
   hMeas_mc_interE->Scale( scale_interE );
   hMeas_mc_totInel->Scale( scale_totInel );
   hMeas_mc_abs->Scale( scale_abs );
   

   //==============================================================================================      

   auto legend_initE = new TLegend(0.5,0.65,0.76,0.85);
   legend_initE->AddEntry(hMeas_data_initE, "PDSP Data, Run 58XX");
   legend_initE->AddEntry(hMeas_mc_initE, "MC simulation");   
   legend_initE->SetTextSize(0.03);
   
   auto legend_interE = new TLegend(0.5,0.65,0.76,0.85);
   legend_interE->AddEntry(hMeas_data_interE, "PDSP Data, Run 58XX");
   legend_interE->AddEntry(hMeas_mc_interE, "MC simulation");   
   legend_interE->SetTextSize(0.03);
   
   auto legend_totInel = new TLegend(0.5,0.65,0.76,0.85);
   legend_totInel->AddEntry(hMeas_data_totInel, "PDSP Data, Run 58XX");
   legend_totInel->AddEntry(hMeas_mc_totInel, "MC simulation");   
   legend_totInel->SetTextSize(0.03);
   
   auto legend_abs = new TLegend(0.5,0.65,0.76,0.85);
   legend_abs->AddEntry(hMeas_data_abs, "PDSP Data, Run 58XX");
   legend_abs->AddEntry(hMeas_mc_abs, "MC simulation");   
   legend_abs->SetTextSize(0.03);
   

   //==============================================================================================      
   TCanvas* c_initE= new TCanvas("canvas_initE","canvas_initE", 1100,800);
   gPad->SetGrid(1,1);
   hMeas_mc_initE->SetFillColorAlpha(kBlue,0.1);
   hMeas_mc_initE->SetMarkerColorAlpha(kBlue,0);
   hMeas_mc_initE->SetLineColor(kBlue);
   hMeas_mc_initE->Draw("HIST");
   //hMeas_data_initE->SetFillColorAlpha(kBlack,0.3);
   hMeas_data_initE->Draw("SAME PE");
   legend_initE->Draw();
   c_initE->Write();
   //c_initE->SaveAs("bla.pdf");

   TCanvas* c_interE= new TCanvas("canvas_interE","canvas_interE", 1100, 800);
   gPad->SetGrid(1,1);
   hMeas_mc_interE->SetFillColorAlpha(kBlue,0.1);
   hMeas_mc_interE->SetMarkerColorAlpha(kBlue,0);
   hMeas_mc_interE->SetLineColor(kBlue);
   hMeas_mc_interE->Draw("HIST");
   //hMeas_data_interE->SetFillColorAlpha(kBlack,0.3);
   hMeas_data_interE->Draw("SAME PE1");
   legend_interE->Draw();
   c_interE->Write();

   TCanvas* c_totInel= new TCanvas("canvas_totInel","canvas_totInel", 1100, 800);
   gPad->SetGrid(1,1);
   hMeas_mc_totInel->SetFillColorAlpha(kBlue,0.1);
   hMeas_mc_totInel->SetMarkerColorAlpha(kBlue,0);
   hMeas_mc_totInel->SetLineColor(kBlue);
   hMeas_mc_totInel->Draw("HIST");
   hMeas_data_totInel->SetFillColorAlpha(kBlack,0.3);
   hMeas_data_totInel->Draw("SAME PE1");
   legend_totInel->Draw();
   c_totInel->Write();

   TCanvas* c_abs= new TCanvas("canvas_abs","canvas_abs", 1100, 800);
   gPad->SetGrid(1,1);
   hMeas_mc_abs->SetFillColorAlpha(kBlue,0.1);
   hMeas_mc_abs->SetMarkerColorAlpha(kBlue,0);
   hMeas_mc_abs->SetLineColor(kBlue);
   hMeas_mc_abs->Draw("HIST");
   //hMeas_data_abs->SetFillColorAlpha(kBlack,0.3);
   hMeas_data_abs->Draw("SAME PE1");
   legend_abs->Draw();
   c_abs->Write();

   //Build the Incident histos in order to compare

   TH1D* hMeas_mc_incident= new TH1D ("meas_mc_incident", "Measured MC incident;Kinetic energy [MeV]; Events / 50 MeV",nBin_new, lowE, highE);
   TH1D* hMeas_data_incident= new TH1D ("meas_data_incident", "Measured Data incident;Kinetic energy [MeV]; Events / 50 MeV",nBin_new, lowE, highE);

   build_incidentHist(hMeas_mc_initE, hMeas_mc_interE, hMeas_mc_incident);
   build_incidentHist(hMeas_data_initE, hMeas_data_interE, hMeas_data_incident);

   auto legend_incident = new TLegend(0.5,0.65,0.76,0.85);
   legend_incident->AddEntry(hMeas_data_incident, "PDSP Data, Run 58XX");
   legend_incident->AddEntry(hMeas_mc_incident, "MC simulation");   
   legend_incident->SetTextSize(0.03);

   TCanvas* c_incident= new TCanvas("canvas_incident","canvas_incident", 1100, 800);
   gPad->SetGrid(1,1);
   hMeas_mc_incident->SetFillColorAlpha(kBlue,0.1);
   hMeas_mc_incident->SetMarkerColorAlpha(kBlue,0);
   hMeas_mc_incident->SetLineColor(kBlue);
   hMeas_mc_incident->Draw("HIST");
   //hMeas_data_incident->SetFillColorAlpha(kBlack,0.3);
   hMeas_data_incident->Draw("SAME PE1");
   legend_incident->Draw();
   c_incident->Write();


}

#ifndef __CINT__
int main () { thesisPlot_meas_intIncHisto(mcFile, dataFile); return 0; }  // Main program when run stand-alone
#endif



/*
// Checking wether categorisation is not giving one event two contributions
// events not in event selection can not have a value and are just not contributing to unfolding
auto h_test = frame
.Define("sum", "truePion_response + truePion_miss + truePion_fake")
.Histo1D("sum");
h_test->Write();

frame
.Define("sum", "truePion_response + truePion_miss + truePion_fake")
.Foreach([](int sum, int pdg, bool selPi, bool recoBin, bool trueBin, double reco_initKE, double reco_interKE, bool res, bool miss, bool fake){

if(sum == 0){
std::cout << "PDG     = " << pdg << std::endl;
std::cout << "SelPi   = " << selPi << std::endl;
std::cout << "recoBin = " << recoBin << std::endl;
std::cout << "trueBin = " << trueBin << std::endl;
std::cout << "rInitKE = " << reco_initKE << std::endl;
std::cout << "rInterE = " << reco_interKE << std::endl;
std::cout << "response = " << res << std::endl;
std::cout << "miss = " << miss << std::endl;
std::cout << "fake = " << fake << std::endl;


}

},{"sum", "true_beam_PDG", "selected_incidentPion", "reco_equalBin", "true_equalBin", "reco_initKE", "reco_interKE", "truePion_response", "truePion_miss", "truePion_fake"});

exit;*/

/*
//==============================================================================================      
//                         Do Ratios between MC True and Reco and apply to Data Reco
//==============================================================================================      

TH1D* hScale_data_wMCratio_totInel= new TH1D ("hScale_data_wMCratio_totInel", "Data totInel with MC true / reco; Energy [MeV]; ev/bin", nBin_int, eEnd, eStart);
TH1D* hScale_data_wMCratio_abs= new TH1D ("hScale_data_wMCratio_abs", "Data abs with MC true / reco; Energy [MeV]; ev/bin", nBin_int, eEnd, eStart);

for(int i=1; i<= nBin_int; i++){

if(hMeas_mc_totInel->GetBinContent(i) != 0){
double ratio_totInel = hTrue_totInel->GetBinContent(i) / hMeas_mc_totInel->GetBinContent(i);
hScale_data_wMCratio_totInel->SetBinContent( i, hMeas_data_totInel->GetBinContent(i) * ratio_totInel );
}
else hScale_data_wMCratio_totInel->SetBinContent(i, hMeas_data_totInel->GetBinContent(i) );

if(hMeas_mc_abs->GetBinContent(i) != 0){
double ratio_abs = hTrue_abs->GetBinContent(i) / hMeas_mc_abs->GetBinContent(i);
hScale_data_wMCratio_abs->SetBinContent(i, hMeas_data_abs->GetBinContent(i) * ratio_abs );
}
else hScale_data_wMCratio_abs->SetBinContent(i, hMeas_data_abs->GetBinContent(i) );


}

//==============================================================================================      
//                         Do Ratios between MC True and Reco and apply to Data Reco
//==============================================================================================      

TH1D* hScale_data_wMCratio_incident= new TH1D ("hScale_data_wMCratio_incident", "Data incident with MC true / reco; Energy [MeV]; ev/bin", nBin_int, eEnd, eStart);

for(int i=1; i<= nBin_int; i++){

if(hMeas_mc_incident->GetBinContent(i) != 0){
double ratio_incident = hTrue_incident->GetBinContent(i) / hMeas_mc_incident->GetBinContent(i);
hScale_data_wMCratio_incident->SetBinContent( i, hMeas_data_incident->GetBinContent(i) * ratio_incident );
}
else hScale_data_wMCratio_incident->SetBinContent( i, hMeas_data_incident->GetBinContent(i) );
}

//Scale MC to Data
hMeas_mc_incident->Scale( hMeas_data_incident->Integral() / hMeas_mc_incident->Integral() );
hTrue_incident->Scale( hMeas_data_incident->Integral() / hTrue_incident->Integral() );

//==============================================================================================      
//RATIO XS
TH1D* h_wRatioXS_totInel = new TH1D("h_wRatioXS_totInel" ,"DATA wRatio totInel XS; Energy [MeV]; #sigma [mb]", nBin_int, eEnd, eStart);

do_XS_log(  h_wRatioXS_totInel, hScale_data_wMCratio_totInel, hScale_data_wMCratio_incident, h_betheMean_muon );   
do_XS_log_binomial_error( h_wRatioXS_totInel, hScale_data_wMCratio_totInel, hScale_data_wMCratio_incident, h_betheMean_muon );  
//==============================================================================================      
//RATIO XS
TH1D* h_wRatioXS_abs = new TH1D("h_wRatioXS_abs" ,"DATA wRatio abs XS; Energy [MeV]; #sigma [mb]", nBin_int, eEnd, eStart);

do_XS_log(  h_wRatioXS_abs, hScale_data_wMCratio_abs, hScale_data_wMCratio_incident, h_betheMean_muon );   
do_XS_log_binomial_error( h_wRatioXS_abs, hScale_data_wMCratio_abs, hScale_data_wMCratio_incident, h_betheMean_muon );  
*/

/*   //======= ADD Muon Content to initE and interE of pions and totInel --> FAKE!! (Extra BG) ===================================================
//Add Muons that were not selected by eventSelection in order to not double-count events
//ist that correct? probably should have many similar muons to those 8% already in the BG...(?)
frame_train
.Filter("primaryMuon && !selected_incidentPion")
.Filter("reco_initKE != -999. && reco_interKE != -999.")
.Range(muExtra_incidentPion)
.Foreach( [ &response_initE, &response_interE, &response_totInel] 
( double reco_initE, double reco_interE){

response_initE.Fake( reco_initE );
response_interE.Fake( reco_interE );
response_totInel.Fake( reco_interE );

},
{"reco_initKE", "reco_interKE"});

frame_train
.Filter("primaryMuon && !selected_abs")
.Filter("reco_initKE != -999. && reco_interKE != -999.")
.Range(muExtra_abs)
.Foreach( [ &response_abs] 
(double reco_interE){
response_abs.Fake( reco_interE );

},
{"reco_interKE"});
*/

