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
const string mcPath = "dataMC_files/eSliceMethod_mc_eventFile.root";
const string dataPath = "dataMC_files/eSliceMethod_data_eventFile.root";
const Double_t dataMC_events = 68500; //root file DATA has 70121 ev, root file MC has 74149 ev


//==============================================================================
//SYSTEMATICS Unfolding for variation of recoKE
//
// Add / remove 10% of Absorption Purity
//
// 
//==============================================================================

void unfoldSystematic_bayes_plusAbsEfficiency(const string mcFilepath, const string dataPath, bool doMC = false, bool doXS = true)
{
   ROOT::RDataFrame inputFrame(pionTree, mcFilepath);
   ROOT::RDataFrame data_inputFrame(pionTree, dataPath);

   gStyle->SetNdivisions(1020);
   //ROOT::RDataFrame data_pre(pionTree, dataPath);

   string output_name;
   //if(doMC) output_name = "unfold_wBG_mc_" + std::to_string((int) bin_size_int) + "MeV.root";
   //else 
   output_name = "unfoldSystematic_" + std::to_string((int) bin_size_int) + "MeV.root";

   TFile *output = new TFile( output_name.c_str() , "UPDATE"); //maybe save with binning?

   //=======================  Frame Definitions ==============================================================      

   //no need to filter for only reconstructed events as with the Miss function one can take into account the not-reconstructed events
   auto frame = inputFrame
      .Filter("pass_trueBeamLocation")// bool that marks in true the particles that passed the true beam location, same range in x and y as for the TPCposition cut for reco but on true var
      //.Filter("primary_isBeamType")
      .Range(dataMC_events)
      .Define("true_equalBin", equalBin, {"true_initKE", "true_interKE"})
      .Define("reco_equalBin", equalBin, {"reco_initKE", "reco_interKE"})
      //Definitions for Response Matrix generation

      .Define("truePion_response", 
            " true_beam_PDG == 211 && selected_incidentPion == 1"
            " && !true_equalBin && !reco_equalBin"
            " && reco_initKE != -999. && reco_interKE != -999.")

      .Define("truePion_miss", //either Missed Signal, or mis Reconstructed Signal
            " true_beam_PDG == 211 && !true_equalBin "
            " && (selected_incidentPion == 0 || reco_equalBin || reco_initKE == -999. || reco_interKE == -999.) ")

      .Define("truePion_fake", 
            " selected_incidentPion == 1 && (true_beam_PDG != 211 || true_equalBin)")

         //Pion Total Inelastic Distribution
      .Define("truePion_abs_response", 
            " true_absSignal && selected_abs && true_pion_daughter == 0"
            " && !true_equalBin && !reco_equalBin"
            " && reco_initKE != -999. && reco_interKE != -999.")

      .Define("truePion_abs_miss", 
            " true_absSignal && !true_equalBin && true_pion_daughter == 0"
            " && ( reco_initKE == -999. || reco_interKE == -999. || reco_equalBin || !selected_abs )")

      .Define("truePion_abs_fake", 
            " selected_abs "
            " && ( !( true_absSignal && true_pion_daughter == 0) || true_equalBin )");


   auto dataFrame = data_inputFrame
      .Range(dataMC_events)
      .Define("reco_equalBin", equalBin, {"reco_initKE", "reco_interKE"});
   //.Filter("reco_initKE < 1200");

   //auto count_data = dataFrame.Count();
   //std::cout << "count data = " << *count_data << std::endl;

   //=====================================================================================      

   //Frame to Train response
   auto frame_train = frame
      .Filter("true_beam_endZ > 0");



   cout << "==================================== TRAIN ====================================" << endl;
   //events that have reco_initBin == reco_interBin go into the Miss function as they are treated like a reco-ineff
   //events with true_initBin == true_interBin havre been filtered out already
   RooUnfoldResponse response_interE = RooUnfoldResponse(nBin_int, eEnd, eStart, 
         "response_syst_plusAbsEfficiency_interE", "");
   RooUnfoldResponse response_initE = RooUnfoldResponse(nBin_int, eEnd, eStart, 
         "response_syst_plusAbsEfficiency_initE", "");
   RooUnfoldResponse response_abs = RooUnfoldResponse(nBin_int, eEnd, eStart, 
         "response_syst_plusAbsEfficiency_abs", "");

   //If there is need to use Overflow do response_bla.UseOverflow(true)

   frame_train
      .Foreach( [ &response_initE, &response_interE](double true_initE, double true_interE, 
               double reco_initE, double reco_interE, 
               bool ev_response, bool ev_miss, bool ev_fake){

            //Signal
            if(ev_response){
            response_initE.Fill( reco_initE, true_initE);
            response_interE.Fill( reco_interE, true_interE);
            }
            //MisReconstruction
            else if( ev_miss ) {
            response_initE.Miss (true_initE);
            response_interE.Miss (true_interE);
            }
            //BackGround will be subtracted before unfolding
            else if( ev_fake ){

            response_initE.Fake(reco_initE);
            response_interE.Fake(reco_interE);
            }

      }
   ,{"true_initKE", "true_interKE", "reco_initKE_rwData", "reco_interKE_rwData", 
      "truePion_response", "truePion_miss", "truePion_fake"});

   //===================================================================================== 
   // Increase  Absorption Efficiency --> Shift from Miss to Response       
   //=====================================================================================  
   //
   double total_abs = *frame_train.Filter("truePion_abs_response").Count();
   double total_miss = *frame_train.Filter("truePion_abs_miss").Count();
   double all_abs = total_abs + total_miss;
   double eff = total_abs / all_abs;
   double new_abs =  (eff + 0.1)*all_abs;
   int delta_abs = (int) new_abs - total_abs;
   int cnt_abs = 0;

   std::cout << "Eff = " << eff << std::endl;
   std::cout << "True Abs = " << total_abs << std::endl;
   std::cout << "Miss Abs = " << total_miss << std::endl;
   std::cout << "All Abs = " << all_abs << std::endl;
   std::cout << "New Abs Sel = " << new_abs << std::endl;
   std::cout << "Delta Abs = " << delta_abs << std::endl;

   frame_train
      .Foreach( [ &response_abs, &cnt_abs, &delta_abs](double true_initE, double true_interE, 
               double reco_initE, double reco_interE, 
               bool ev_response, bool ev_miss, bool ev_fake){

            //Signal
            if(ev_response ){
            response_abs.Fill(reco_interE, true_interE);
            }
             //increase response matrix
            else if( ev_miss && cnt_abs <= delta_abs ) {
            response_abs.Fill (reco_interE, true_interE);
            ++cnt_abs;
            }
             //MisReconstruction
            else if( ev_miss && cnt_abs > delta_abs) {
            response_abs.Miss (true_interE);
            }
            //Reduce BG
            else if( ev_fake){
            response_abs.Fake(reco_interE);
            }

            }
            ,{"true_initKE", "true_interKE", "reco_initKE_rwData", "reco_interKE_rwData", 
            "truePion_abs_response", "truePion_abs_miss", "truePion_abs_fake"});





   cout << "==================================== TEST =====================================" << endl;


   TH1D* hTrue_initE= new TH1D ("trueMC_initE", "MC True; Energy [MeV]; Events / 50 MeV", nBin_int, eEnd, eStart);
   TH1D* hMeas_data_initE= new TH1D ("syst_plusAbsEfficiency_initKE", "; Energy [MeV]; Events / 50 MeV", nBin_int, eEnd, eStart);

   TH1D* hTrue_interE= new TH1D ("trueMC_interE", "MC True; Energy [MeV]; Events / 50 MeV",    nBin_int, eEnd, eStart);
   TH1D* hMeas_data_interE= new TH1D ("syst_plusAbsEfficiency_interKE", "Measured data Run 58XX interE; Energy [MeV]; Events / 50 MeV", nBin_int, eEnd, eStart);

   TH1D* hTrue_abs= new TH1D ("trueMC_abs", "MC True; Energy [MeV]; Events / 50 MeV",    nBin_int, eEnd, eStart);
   TH1D* hMeas_data_abs= new TH1D ("syst_plusAbsEfficiency_abs", "Measured data Run 58XX absorption; Energy [MeV]; Events / 50 MeV", nBin_int, eEnd, eStart);

   //Build hMeas from same prefiltered sample that was supplied to response and fake (the Miss are misReco so won't catch them in reco)
   //Build hTrue from what was part of response and Miss (the fake are not part of hTrue

   //====================  DATA ==========================
   dataFrame
      .Filter("selected_incidentPion")
      .Foreach( [ hMeas_data_initE, hMeas_data_interE ]( double reco_initE, double reco_interE, bool recoEqualBin){

            //Measured
            if(reco_initE != -999. && reco_interE != -999. && !recoEqualBin){

            hMeas_data_initE->Fill( reco_initE );
            hMeas_data_interE->Fill( reco_interE );
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

   //hMeas_data_initE->Write();
   //hMeas_data_interE->Write();
   //hMeas_data_abs->Write();


   cout << "==================================== Building the True Distribution =====================================" << endl;


   frame_train
      .Foreach( [ hTrue_initE, hTrue_interE ]( double true_initE, double true_interE, 
               bool ev_response, bool ev_miss){

            if( ev_response || ev_miss ){

            hTrue_initE->Fill( true_initE );
            hTrue_interE->Fill( true_interE );
            }

            }
            ,{"true_initKE", "true_interKE", "truePion_response", "truePion_miss"}
            );

   frame_train
      .Foreach( [ hTrue_abs ]( double true_initE, double true_interE, 
               bool ev_response, bool ev_miss){

            if( ev_response || ev_miss ){
            hTrue_abs->Fill( true_interE );
            }

            }
            ,{ "true_initKE", "true_interKE", "truePion_abs_response", "truePion_abs_miss"}
            );

   //hTrue_initE->Write();
   //hTrue_interE->Write();
   //hTrue_abs->Write();


   RooUnfoldBayes unfold_initE;
   RooUnfoldBayes unfold_interE;
   RooUnfoldBayes unfold_abs;

      cout << "==================================== UNFOLD DATA ===================================" << endl;
      unfold_initE = RooUnfoldBayes(&response_initE, hMeas_data_initE, 4);    // OR
      unfold_interE = RooUnfoldBayes(&response_interE, hMeas_data_interE, 6);    // OR
      unfold_abs = RooUnfoldBayes(&response_abs, hMeas_data_abs, 6);    // OR



   cout << "=========================== Covariance Matrices ===================================" << endl;

   TMatrixD cov_initE = unfold_initE.Ereco();
   //unfold_initE.SetMeasuredCov( cov_initE );

   TMatrixD cov_interE = unfold_interE.Ereco();
   //unfold_interE.SetMeasuredCov( cov_interE );


   TMatrixD cov_abs = unfold_abs.Ereco();
   //unfold_abs.SetMeasuredCov( cov_abs ); 

   cout << "------------------------------------ Init E -----------------------------------" << endl;
   TH1D* hUnfold_initE = (TH1D*) unfold_initE.Hreco();
   cout << "------------------------------------ Inter E -----------------------------------" << endl;
   TH1D* hUnfold_interE = (TH1D*) unfold_interE.Hreco();
   cout << "------------------------------------ abs -----------------------------------" << endl;
   TH1D* hUnfold_abs = (TH1D*) unfold_abs.Hreco();

   

   hUnfold_initE->SetNameTitle("hSystUnfold_plusAbsEfficiency_initE", "; Energy [MeV]; Events / 50 MeV");
   hUnfold_interE->SetNameTitle("hSystUnfold_plusAbsEfficiency_interE", "; Energy [MeV]; Events / 50 MeV");
   hUnfold_abs->SetNameTitle("hSystUnfold_plusAbsEfficiency_abs", "; Energy [MeV]; Events / 50 MeV");
   
   hUnfold_initE->Write();
   hUnfold_interE->Write();
   hUnfold_abs->Write();
   //print infos with PrintTable command
   cout << "==================================== UNFOLD INIT E===================================" << endl;
   unfold_initE.PrintTable (cout, hTrue_initE);
   cout << "==================================== UNFOLD INTER E===================================" << endl;
   unfold_interE.PrintTable (cout, hTrue_interE);
   cout << "==================================== UNFOLD abs E===================================" << endl;
   unfold_abs.PrintTable (cout, hTrue_abs);



   cout << "=========================== CHECK INTEGRALS INITE / INTER E===================================" << endl;
   cout << "Integral of Measured Data initE = " << hMeas_data_initE->Integral() << endl;
   cout << "Integral of Measured Data interE = " << hMeas_data_interE->Integral() << endl;
   cout << "Integral of True initE = " << hTrue_initE->Integral() << endl;
   cout << "Integral of True interE = " << hTrue_interE->Integral() << endl;
   cout << "Integral of Unfold initE = " << hUnfold_initE->Integral() << endl;
   cout << "Integral of Unfold interE = " << hUnfold_interE->Integral() << endl;

   cout << "=========================== CHECK INTEGRALS INITE / INTER E===================================" << endl;
   cout << "Integral of True abs = " << hTrue_abs->Integral() << endl;
   cout << "Integral of Unfold abs = " << hUnfold_abs->Integral() << endl;

   cout << "=========================== Correlation Matrices ===================================" << endl;
   cout << "------------------------------------ init E -----------------------------------" << endl;

   Int_t nb_initE = cov_initE.GetNrows();
   TH2D* corr_initE = new TH2D("corr_syst_plusAbsEfficiency_initE", "Correlations Init E; Pion initE [MeV];Pion initE [MeV];", nb_initE, eEnd, eStart, nb_initE, eEnd, eStart);
   correlationMatrix( cov_initE, corr_initE, nb_initE);
   corr_initE->Write();

   cout << "------------------------------------ inter E -----------------------------------" << endl;

   Int_t nb_interE = cov_interE.GetNrows();
   TH2D* corr_interE = new TH2D("corr_syst_plusAbsEfficiency_interE", "Correlations Inter E; Pion interE [MeV];Pion interE [MeV];", nb_interE, eEnd, eStart, nb_interE, eEnd, eStart);
   correlationMatrix( cov_interE, corr_interE, nb_interE);
   corr_interE->Write();

   cout << "------------------------------------ Absorption -----------------------------------" << endl;

   Int_t nb_abs = cov_abs.GetNrows();
   TH2D* corr_abs = new TH2D("corr_syst_plusAbsEfficiency_abs", "Correlations Abs; Absorption interE [MeV];Absorption interE [MeV] E", nb_abs, eEnd, eStart, nb_abs, eEnd, eStart);
   correlationMatrix( cov_abs, corr_abs, nb_abs);
   corr_abs->Write();

   cout << "=========================== Pulls  ===================================" << endl;
   //TH1D* hPull_initE= new TH1D ("hPull_initE", "Pulls initE; Energy [MeV]; Events / 50 MeV", nBin_int, eEnd, eStart);
   //TH1D* hPull_interE= new TH1D ("hPull_interE", "Pulls interE; Energy [MeV]; Events / 50 MeV",    nBin_int, eEnd, eStart);
   //TH1D* hPull_abs= new TH1D ("hPull_abs", "Pulls abs; Energy [MeV]; Events / 50 MeV",    nBin_int, eEnd, eStart);

   //pullHisto( hUnfold_initE, hTrue_initE, hPull_initE);
   //pullHisto( hUnfold_interE, hTrue_interE, hPull_interE);
   //pullHisto( hUnfold_totInel, hTrue_totInel, hPull_totInel);
   //pullHisto( hUnfold_abs, hTrue_abs, hPull_abs);

   //hPull_initE->Write(); hPull_interE->Write(); hPull_totInel->Write(); hPull_abs->Write();


   //==============================================================================================      
   auto legend_initE = new TLegend(0.5,0.65,0.76,0.85);
   legend_initE->AddEntry(hUnfold_initE, "Unfold Syst plusAbsEfficiency");
   legend_initE->AddEntry(hMeas_data_initE, "PDSP Data, Run 58XX initE");
   legend_initE->AddEntry(hTrue_initE, "MC True");
   legend_initE->SetTextSize(0.03);

   auto legend_interE = new TLegend(0.5,0.65,0.76,0.85);
   legend_interE->AddEntry(hUnfold_interE, "Unfold Syst plusAbsEfficiency");
   legend_interE->AddEntry(hMeas_data_interE, "PDSP Data, Run 58XX interE");
   legend_interE->AddEntry(hTrue_interE, "MC True");
   legend_interE->SetTextSize(0.03);

   auto legend_abs = new TLegend(0.5,0.65,0.76,0.85);
   legend_abs->AddEntry(hUnfold_abs, "Unfold Syst plusAbsEfficiency");
   legend_abs->AddEntry(hMeas_data_abs, "PDSP Data, Run 58XX absorption");
   legend_abs->AddEntry(hTrue_abs, "MC True");
   legend_abs->SetTextSize(0.03);

   //==============================================================================================      
   gStyle->SetOptFit(0);
   gStyle->SetOptStat(0);

   TCanvas* c_initE= new TCanvas("cSyst_plusAbsEfficiency_initE","cSyst_plusAbsEfficiency_initE", 1100, 800);
   gPad->SetGrid(1,1);
   hTrue_initE->SetLineColor(8);
   hTrue_initE->SetMarkerColorAlpha(kBlack,0);
   hTrue_initE->Draw("HIST ");
   hMeas_data_initE->SetMarkerColorAlpha(kBlack,0);
   hMeas_data_initE->SetFillColorAlpha(kBlack,0.2);
   hMeas_data_initE->Draw("SAME HIST");
   hUnfold_initE->Draw("PE1 SAME");
   legend_initE->Draw();
   //c_initE->Write();
   //c_initE->SaveAs("bla.pdf");

   TCanvas* c_interE= new TCanvas("cSyst_plusAbsEfficiency_interE","cSyst_plusAbsEfficiency_interE", 1100, 800);
   gPad->SetGrid(1,1);
   hTrue_interE->SetLineColor(8);
   hTrue_interE->SetMarkerColorAlpha(kBlack,0);
   hTrue_interE->Draw("HIST ");
   hMeas_data_interE->SetMarkerColorAlpha(kBlack,0);
   hMeas_data_interE->SetFillColorAlpha(kBlack,0.2);
   hMeas_data_interE->Draw("SAME HIST");
   hUnfold_interE->Draw("PE1 SAME");
   legend_interE->Draw();
   //c_interE->Write();
   //c_interE->SaveAs("bla.pdf");

   TCanvas* c_abs= new TCanvas("cSyst_plusAbsEfficiency_abs","cSyst_plusAbsEfficiency_abs", 1100, 800);
   gPad->SetGrid(1,1);
   hTrue_abs->SetLineColor(8);
   hTrue_abs->SetMarkerColorAlpha(kBlack,0);
   hTrue_abs->Draw("HIST ");
   hMeas_data_abs->SetMarkerColorAlpha(kBlack,0);
   hMeas_data_abs->SetFillColorAlpha(kBlack,0.2);
   hMeas_data_abs->Draw("SAME HIST");
   hUnfold_abs->Draw("PE1 SAME");
   legend_abs->Draw();
   //c_abs->Write();
   //c_abs->SaveAs("bla.pdf");

   //Build the Incident histos in order to compare

   TH1D* hTrue_incident= new TH1D ("trueMC_incident", "MC Truth incident; Energy [MeV]; Events / 50 MeV",    nBin_int, eEnd, eStart);
   TH1D* hMeas_data_incident= new TH1D ("meas_data_incident", "Measured Data incident; Energy [MeV]; Events / 50 MeV", nBin_int, eEnd, eStart);
   TH1D* hUnfold_incident= new TH1D ("hSystUnfold_plusAbsEfficiency_incident", "; Energy [MeV]; Events / 50 MeV", nBin_int, eEnd, eStart);

   build_incidentHist(hTrue_initE, hTrue_interE, hTrue_incident);
   build_incidentHist(hUnfold_initE, hUnfold_interE, hUnfold_incident);
   build_incidentHist(hMeas_data_initE, hMeas_data_interE, hMeas_data_incident);

   auto legend_incident = new TLegend(0.5,0.65,0.76,0.85);
   legend_incident->AddEntry(hUnfold_incident, "Unfold Syst plusAbsEfficiency");
   legend_incident->AddEntry(hMeas_data_incident, "PDSP Data, Run 58XX incident");
   legend_incident->AddEntry(hTrue_incident, "MC True");
   legend_incident->SetTextSize(0.03);


   TCanvas* c_incident= new TCanvas("cSyst_plusAbsEfficiency_incident","cSyst_plusAbsEfficiency_incident", 1100, 800);
   gPad->SetGrid(1,1);
   hTrue_incident->SetLineColor(8);
   hTrue_incident->SetMarkerColorAlpha(kBlack,0);
   hTrue_incident->Draw("HIST ");
   hMeas_data_incident->SetMarkerColorAlpha(kBlack,0);
   hMeas_data_incident->SetFillColorAlpha(kBlack,0.2);
   hMeas_data_incident->Draw("SAME HIST");
   hUnfold_incident->Draw("PE1 SAME");
   legend_incident->Draw();
   //c_incident->Write();
   //c_incident->SaveAs("bla.pdf");

   auto* R_initE = response_initE.HresponseNoOverflow();
   auto* c_response_initE = new TCanvas();
   R_initE->SetStats(0);
   R_initE->Draw("colz");
   R_initE->GetXaxis()->SetNdivisions(1020);
   R_initE->GetYaxis()->SetNdivisions(1020);
   R_initE->SetTitle("; Reco Pion InitE [MeV]; True Pion InitE [MeV]");
   gPad->SetGrid(1,1);
   //c_response_initE->Draw();

   auto* R_interE = response_interE.HresponseNoOverflow();
   auto* c_response_interE = new TCanvas();
   R_interE->SetStats(0);
   R_interE->Draw("colz");
   R_interE->GetXaxis()->SetNdivisions(1020);
   R_interE->GetYaxis()->SetNdivisions(1020);
   R_interE->SetTitle("; Reco Pion InterE [MeV]; True Pion InterE [MeV]");
   gPad->SetGrid(1,1);
   //c_response_interE->Draw();

   auto* R_abs = response_abs.HresponseNoOverflow();
   auto* c_response_abs = new TCanvas();
   R_abs->SetStats(0);
   R_abs->Draw("colz");
   R_abs->GetXaxis()->SetNdivisions(1020);
   R_abs->GetYaxis()->SetNdivisions(1020);
   R_abs->SetTitle("; Reco Abs InterE [MeV]; True Abs InterE [MeV]");
   gPad->SetGrid(1,1);
   //c_response_abs->Draw();
}

#ifndef __CINT__
int main () { unfoldSystematic_bayes_plusAbsEfficiency(mcPath, dataPath); return 0; }  // Main program when run stand-alone
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
.Filter("reco_initKE_rwData != -999. && reco_interKE_rwData != -999.")
.Range(muExtra_incidentPion)
.Foreach( [ &response_initE, &response_interE, &response_totInel] 
( double reco_initE, double reco_interE){

response_initE.Fake( reco_initE );
response_interE.Fake( reco_interE );
response_totInel.Fake( reco_interE );

},
{"reco_initKE_rwData", "reco_interKE_rwData"});

frame_train
.Filter("primaryMuon && !selected_abs")
.Filter("reco_initKE_rwData != -999. && reco_interKE_rwData != -999.")
.Range(muExtra_abs)
.Foreach( [ &response_abs] 
(double reco_interE){
response_abs.Fake( reco_interE );

},
{"reco_interKE_rwData"});
*/

