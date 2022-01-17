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

#include "lambda.h"
#include "betheBloch.h"
#include "eSlice.h"

//==============================================================================
// Global definitions
//==============================================================================

const Double_t cutdummy= -99999.0;
const string mcPath = "dataMC_files/eSliceMethod_mc_eventFile.root";
const string dataPath = "dataMC_files/eSliceMethod_data_eventFile.root";
//const string dataPath = "eSliceMethod_Prod4a_58XX_1GeV_all_09_17_21.root";
const Double_t dataMC_events = 68500; //root file DATA has 70121 ev, root file MC has 74149 ev


//==============================================================================
// Example Unfolding
//==============================================================================

void eSliceMethod_unfoldBayes_XS_dataMC(const string mcFilepath, const string dataPath, bool doMC = false, bool doXS = true)
{
   ROOT::RDataFrame inputFrame(pionTree, mcFilepath);
   ROOT::RDataFrame data_inputFrame(pionTree, dataPath);

   gStyle->SetNdivisions(1020);

   string output_name;
   output_name = "unfold_xs_data_" + std::to_string((int) bin_size_int) + "MeV.root";

   TFile *output = new TFile( output_name.c_str() , "RECREATE"); //maybe save with binning?

   //=======================  Frame Definitions ==============================================================      

   //no need to filter for only reconstructed events as with the Miss function one can take into account the not-reconstructed events
   auto frame = inputFrame
      .Filter("pass_trueBeamLocation")// bool that marks in true the particles that passed the true beam location, same range in x and y as for the TPCposition cut for reco but on true var
      //.Filter("primary_isBeamType")
      .Range(dataMC_events)
      .Define("true_equalBin", equalBin, {"true_initKE", "true_interKE"})
      .Define("reco_equalBin", equalBin, {"reco_initKE_rwData", "reco_interKE_rwData"})
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
      .Define("truePion_totInel_response", 
            " true_primPionInel && selected_incidentPion"
            " && !true_equalBin && !reco_equalBin"
            " && reco_initKE != -999. && reco_interKE != -999.")

      .Define("truePion_totInel_miss", 
            " true_primPionInel && !true_equalBin"
            " && ( reco_initKE == -999. || reco_interKE == -999. || reco_equalBin || !selected_incidentPion )")

      .Define("truePion_totInel_fake", 
            " selected_incidentPion && ( !true_primPionInel || true_equalBin ) ")

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

   //auto count_data = dataFrame.Count();
   //std::cout << "count data = " << *count_data << std::endl;

   //=====================================================================================      

   //Frame to Train response
   auto frame_train = frame
      .Filter("true_beam_endZ > 0");

   //Frame for Measurement 
   auto frame_meas = frame
      .Filter("true_beam_endZ > 0");


   cout << "==================================== TRAIN ====================================" << endl;
   //events that have reco_initBin == reco_interBin go into the Miss function as they are treated like a reco-ineff
   //events with true_initBin == true_interBin havre been filtered out already
   RooUnfoldResponse response_interE = RooUnfoldResponse(nBin_int, eEnd, eStart, 
         "response_interE", "");
   RooUnfoldResponse response_initE = RooUnfoldResponse(nBin_int, eEnd, eStart, 
         "response_initE", "");
   RooUnfoldResponse response_totInel = RooUnfoldResponse(nBin_int, eEnd, eStart, 
         "response_totInel", "");
   RooUnfoldResponse response_abs = RooUnfoldResponse(nBin_int, eEnd, eStart, 
         "response_abs", "");

   //If there is need to use Overflow do response_bla.UseOverflow(true)

   frame_train
      //.Range(half_mc)
      //.Range(53000) 
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

   //Total Inelastic
   frame_train
      //.Range(half_mc)
      //.Range(53000) 
      .Foreach( [ &response_totInel](double true_initE, double true_interE, 
               double reco_initE, double reco_interE, 
               bool ev_response, bool ev_miss, bool ev_fake){

            //Signal
            if(ev_response){
            response_totInel.Fill(reco_interE, true_interE);
            }
            //MisReconstruction
            else if( ev_miss ) {
            response_totInel.Miss (true_interE);
            }
            //BackGround will be subtracted before unfolding
            else if( ev_fake ){
            response_totInel.Fake(reco_interE);
            }

            }
            ,{"true_initKE", "true_interKE", "reco_initKE_rwData", "reco_interKE_rwData", 
            "truePion_totInel_response", "truePion_totInel_miss", "truePion_totInel_fake"});

   //Total Inelastic
   frame_train
      //.Range(half_mc)
      //.Range(53000)
      .Foreach( [ &response_abs](double true_initE, double true_interE, 
               double reco_initE, double reco_interE, 
               bool ev_response, bool ev_miss, bool ev_fake){

            //Signal
            if(ev_response){
            response_abs.Fill(reco_interE, true_interE);
            }
            //MisReconstruction
            else if( ev_miss ) {
            response_abs.Miss (true_interE);
            }
            //BackGround will be subtracted before unfolding
            else if( ev_fake ){
            response_abs.Fake(reco_interE);
            }

            }
            ,{"true_initKE", "true_interKE", "reco_initKE_rwData", "reco_interKE_rwData", 
            "truePion_abs_response", "truePion_abs_miss", "truePion_abs_fake"});

   cout << "==================================== TEST =====================================" << endl;


   TH1D* hTrue_initE= new TH1D ("trueMC_initE", "MC True; Energy [MeV]; Events / 50 MeV", nBin_int, eEnd, eStart);
   TH1D* hMeas_mc_initE= new TH1D ("meas_mc_initE", "Measured MC initE; Energy [MeV]; Events / 50 MeV", nBin_int, eEnd, eStart);
   TH1D* hMeas_data_initE= new TH1D ("meas_data_initE", "Measured data Run 58XX initE; Energy [MeV]; Events / 50 MeV", nBin_int, eEnd, eStart);

   TH1D* hTrue_interE= new TH1D ("trueMC_interE", "MC True; Energy [MeV]; Events / 50 MeV",    nBin_int, eEnd, eStart);
   TH1D* hMeas_mc_interE= new TH1D ("meas_mc_interE", "Measured MC interE; Energy [MeV]; Events / 50 MeV", nBin_int, eEnd, eStart);
   TH1D* hMeas_data_interE= new TH1D ("meas_data_interE", "Measured data Run 58XX interE; Energy [MeV]; Events / 50 MeV", nBin_int, eEnd, eStart);

   TH1D* hTrue_totInel= new TH1D ("trueMC_totInel", "MC True; Energy [MeV]; Events / 50 MeV",    nBin_int, eEnd, eStart);
   TH1D* hMeas_mc_totInel= new TH1D ("meas_mc_totInel", "Measured MC totInel; Energy [MeV]; Events / 50 MeV", nBin_int, eEnd, eStart);
   TH1D* hMeas_data_totInel= new TH1D ("meas_data_totInel", "Measured data Run 58XX totInel; Energy [MeV]; Events / 50 MeV", nBin_int, eEnd, eStart);

   TH1D* hTrue_abs= new TH1D ("trueMC_abs", "MC True; Energy [MeV]; Events / 50 MeV",    nBin_int, eEnd, eStart);
   TH1D* hMeas_mc_abs= new TH1D ("meas_mc_abs", "Measured MC abs; Energy [MeV]; Events / 50 MeV", nBin_int, eEnd, eStart);
   TH1D* hMeas_data_abs= new TH1D ("meas_data_abs", "Measured data Run 58XX absorption; Energy [MeV]; Events / 50 MeV", nBin_int, eEnd, eStart);

   //Build hMeas from same prefiltered sample that was supplied to response and fake (the Miss are misReco so won't catch them in reco)
   //Build hTrue from what was part of response and Miss (the fake are not part of hTrue

   cout << "==================================== Building the Measured Distribution =====================================" << endl;

   //==================== MC   ==========================
   frame_meas

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
            

   frame_meas
      .Filter("selected_abs")
      .Foreach( [ hMeas_mc_abs ]( double reco_initE, double reco_interE, bool recoEqualBin){

            //Measured
            if(reco_initE != -999. && reco_interE != -999. && !recoEqualBin){
            hMeas_mc_abs->Fill( reco_interE );
            }

            /*for(int i = 1; i <= nBin_int; i++){
              hMeas_mc_abs->SetBinError( i , hMeas_mc_abs->GetBinContent(i)*0.05 );

              }*/
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


   cout << "==================================== Building the True Distribution =====================================" << endl;


   frame_meas

      .Foreach( [ hTrue_initE, hTrue_interE ]( double true_initE, double true_interE, 
               bool ev_response, bool ev_miss){

            if( ev_response || ev_miss ){

            hTrue_initE->Fill( true_initE );
            hTrue_interE->Fill( true_interE );
            }

            }
            ,{"true_initKE", "true_interKE", "truePion_response", "truePion_miss"}
            );
   frame_meas

      .Foreach( [ hTrue_totInel ]( double true_initE, double true_interE, 
               bool ev_response, bool ev_miss){

            if( ev_response || ev_miss ){
            hTrue_totInel->Fill( true_interE );
            }

            }
            ,{ "true_initKE", "true_interKE", "truePion_totInel_response", "truePion_totInel_miss"}
            );

   frame_meas

      .Foreach( [ hTrue_abs ]( double true_initE, double true_interE, 
               bool ev_response, bool ev_miss){

            if( ev_response || ev_miss ){
            hTrue_abs->Fill( true_interE );
            }

            }
            ,{ "true_initKE", "true_interKE", "truePion_abs_response", "truePion_abs_miss"}
            );

   hTrue_initE->Write();
   hTrue_interE->Write();
   hTrue_totInel->Write();
   hTrue_abs->Write();


   RooUnfoldBayes unfold_initE;
   RooUnfoldBayes unfold_interE;
   RooUnfoldBayes unfold_totInel;
   RooUnfoldBayes unfold_abs;

   if(doMC){
      cout << "==================================== UNFOLD MC===================================" << endl;
      unfold_initE = RooUnfoldBayes(&response_initE, hMeas_mc_initE, 4);    // OR
      unfold_interE = RooUnfoldBayes(&response_interE, hMeas_mc_interE, 6);    // OR
      unfold_totInel = RooUnfoldBayes(&response_totInel, hMeas_mc_totInel, 6);    // OR
      unfold_abs = RooUnfoldBayes(&response_abs, hMeas_mc_abs, 6);    // OR
   }

   else{
      cout << "==================================== UNFOLD DATA ===================================" << endl;
      unfold_initE = RooUnfoldBayes(&response_initE, hMeas_data_initE, 4);    // OR
      unfold_interE = RooUnfoldBayes(&response_interE, hMeas_data_interE, 6);    // OR
      unfold_totInel = RooUnfoldBayes(&response_totInel, hMeas_data_totInel, 6);    // OR
      unfold_abs = RooUnfoldBayes(&response_abs, hMeas_data_abs, 6);    // OR
   }


   cout << "=========================== Covariance Matrices ===================================" << endl;

   TMatrixD cov_initE = unfold_initE.Ereco();
   //unfold_initE.SetMeasuredCov( cov_initE );

   TMatrixD cov_interE = unfold_interE.Ereco();
   //unfold_interE.SetMeasuredCov( cov_interE );

   TMatrixD cov_totInel = unfold_totInel.Ereco();
   //unfold_totInel.SetMeasuredCov( cov_totInel );

   TMatrixD cov_abs = unfold_abs.Ereco();
   //unfold_abs.SetMeasuredCov( cov_abs ); 

   cout << "------------------------------------ Init E -----------------------------------" << endl;
   //TH1D* hUnfold_initE = (TH1D*) unfold_initE.Hreco((RooUnfold::ErrorTreatment)RooUnfold::kCovariance);
   TH1D* hUnfold_initE = (TH1D*) unfold_initE.Hreco();
   cout << "------------------------------------ Inter E -----------------------------------" << endl;
   //TH1D* hUnfold_interE = (TH1D*) unfold_interE.Hreco((RooUnfold::ErrorTreatment)RooUnfold::kCovariance);
   TH1D* hUnfold_interE = (TH1D*) unfold_interE.Hreco();
   cout << "------------------------------------ TotInel -----------------------------------" << endl;
   //TH1D* hUnfold_totInel = (TH1D*) unfold_totInel.Hreco((RooUnfold::ErrorTreatment)RooUnfold::kCovariance);
   TH1D* hUnfold_totInel = (TH1D*) unfold_totInel.Hreco();
   cout << "------------------------------------ abs -----------------------------------" << endl;
   //TH1D* hUnfold_abs = (TH1D*) unfold_abs.Hreco((RooUnfold::ErrorTreatment)RooUnfold::kCovariance);
   TH1D* hUnfold_abs = (TH1D*) unfold_abs.Hreco();

   

   if(doMC){
      hUnfold_initE->SetNameTitle("hUnfold_initE", "Unfold MC initE; Energy [MeV]; Events / 50 MeV");
      hUnfold_interE->SetNameTitle("hUnfold_interE", "Unfold MC interE; Energy [MeV]; Events / 50 MeV");
      hUnfold_totInel->SetNameTitle("hUnfold_totInel", "Unfold MC totInel; Energy [MeV]; Events / 50 MeV");
      hUnfold_abs->SetNameTitle("hUnfold_abs", "Unfold MC abs; Energy [MeV]; Events / 50 MeV");
   }
   else{
      hUnfold_initE->SetNameTitle("hUnfold_initE", "Unfold data initE; Energy [MeV]; Events / 50 MeV");
      hUnfold_interE->SetNameTitle("hUnfold_interE", "Unfold data interE; Energy [MeV]; Events / 50 MeV");
      hUnfold_totInel->SetNameTitle("hUnfold_totInel", "Unfold data totInel; Energy [MeV]; Events / 50 MeV");
      hUnfold_abs->SetNameTitle("hUnfold_abs", "Unfold data abs; Energy [MeV]; Events / 50 MeV");
   }
   hUnfold_initE->Write();
   hUnfold_interE->Write();
   hUnfold_totInel->Write();
   hUnfold_abs->Write();
   //print infos with PrintTable command
   cout << "==================================== UNFOLD INIT E===================================" << endl;
   unfold_initE.PrintTable (cout, hTrue_initE);
   cout << "==================================== UNFOLD INTER E===================================" << endl;
   unfold_interE.PrintTable (cout, hTrue_interE);
   cout << "==================================== UNFOLD TOTINEL E===================================" << endl;
   unfold_totInel.PrintTable (cout, hTrue_totInel);
   cout << "==================================== UNFOLD abs E===================================" << endl;
   unfold_abs.PrintTable (cout, hTrue_abs);

   //cout << "=========================== SCALE MC MEas TO DATA Meas by AREA AFTER having Added MUON CONTENT===================================" << endl;
   //



   double scale_initE = hMeas_data_initE->Integral() / hMeas_mc_initE->Integral();
   //double scale_interE = hMeas_data_interE->Integral() / hMeas_mc_interE->Integral();
   //double scale_totInel = hMeas_data_totInel->Integral() / hMeas_mc_totInel->Integral();
   //double scale_abs = hMeas_data_abs->Integral() / hMeas_mc_abs->Integral();

   if(!doMC){
      hMeas_mc_initE->Scale( scale_initE );
      hMeas_mc_interE->Scale( scale_initE );
      hMeas_mc_totInel->Scale( scale_initE );
      hMeas_mc_abs->Scale( scale_initE );

      hTrue_initE->Scale( scale_initE );
      hTrue_interE->Scale( scale_initE );
      hTrue_totInel->Scale( scale_initE );
      hTrue_abs->Scale( scale_initE );
   }


   cout << "=========================== CHECK INTEGRALS INITE / INTER E===================================" << endl;
   cout << "Integral of Measured MC initE = " << hMeas_mc_initE->Integral() << endl;
   cout << "Integral of Measured MC interE = " << hMeas_mc_interE->Integral() << endl;
   cout << "Integral of Measured Data initE = " << hMeas_data_initE->Integral() << endl;
   cout << "Integral of Measured Data interE = " << hMeas_data_interE->Integral() << endl;
   cout << "Integral of True initE = " << hTrue_initE->Integral() << endl;
   cout << "Integral of True interE = " << hTrue_interE->Integral() << endl;
   cout << "Integral of Unfold initE = " << hUnfold_initE->Integral() << endl;
   cout << "Integral of Unfold interE = " << hUnfold_interE->Integral() << endl;

   cout << "=========================== CHECK INTEGRALS INITE / INTER E===================================" << endl;
   cout << "Integral of Measured totInel = " << hMeas_mc_totInel->Integral() << endl;
   cout << "Integral of Measured abs = " << hMeas_mc_abs->Integral() << endl;
   cout << "Integral of True totInel = " << hTrue_totInel->Integral() << endl;
   cout << "Integral of True abs = " << hTrue_abs->Integral() << endl;
   cout << "Integral of Unfold totInel = " << hUnfold_totInel->Integral() << endl;
   cout << "Integral of Unfold abs = " << hUnfold_abs->Integral() << endl;

   cout << "=========================== Correlation Matrices ===================================" << endl;
   cout << "------------------------------------ init E -----------------------------------" << endl;

   Int_t nb_initE = cov_initE.GetNrows();
   TH2D* corr_initE = new TH2D("corr_initE", "Correlations Init E; Pion initE [MeV];Pion initE [MeV];", nb_initE, eEnd, eStart, nb_initE, eEnd, eStart);
   correlationMatrix( cov_initE, corr_initE, nb_initE);
   corr_initE->Write();

   cout << "------------------------------------ inter E -----------------------------------" << endl;

   Int_t nb_interE = cov_interE.GetNrows();
   TH2D* corr_interE = new TH2D("corr_interE", "Correlations Inter E; Pion interE [MeV];Pion interE [MeV];", nb_interE, eEnd, eStart, nb_interE, eEnd, eStart);
   correlationMatrix( cov_interE, corr_interE, nb_interE);
   corr_interE->Write();

   cout << "------------------------------------ TotInel -----------------------------------" << endl;

   Int_t nb_totInel = cov_totInel.GetNrows();
   TH2D* corr_totInel = new TH2D("corr_totInel", "Correlations TotInel E", nb_totInel, eEnd, eStart, nb_totInel, eEnd, eStart);
   correlationMatrix( cov_totInel, corr_totInel, nb_totInel);
   corr_totInel->Write();

   cout << "------------------------------------ Absorption -----------------------------------" << endl;

   Int_t nb_abs = cov_abs.GetNrows();
   TH2D* corr_abs = new TH2D("corr_abs", "Correlations Abs; Absorption interE [MeV];Absorption interE [MeV] E", nb_abs, eEnd, eStart, nb_abs, eEnd, eStart);
   correlationMatrix( cov_abs, corr_abs, nb_abs);
   corr_abs->Write();

   cout << "=========================== Pulls  ===================================" << endl;
   TH1D* hPull_initE= new TH1D ("hPull_initE", "Pulls initE; Energy [MeV]; Events / 50 MeV", nBin_int, eEnd, eStart);
   TH1D* hPull_interE= new TH1D ("hPull_interE", "Pulls interE; Energy [MeV]; Events / 50 MeV",    nBin_int, eEnd, eStart);
   TH1D* hPull_totInel= new TH1D ("hPull_totInel", "Pulls totInel; Energy [MeV]; Events / 50 MeV",    nBin_int, eEnd, eStart);
   TH1D* hPull_abs= new TH1D ("hPull_abs", "Pulls abs; Energy [MeV]; Events / 50 MeV",    nBin_int, eEnd, eStart);

   pullHisto( hUnfold_initE, hTrue_initE, hPull_initE);
   pullHisto( hUnfold_interE, hTrue_interE, hPull_interE);
   pullHisto( hUnfold_totInel, hTrue_totInel, hPull_totInel);
   pullHisto( hUnfold_abs, hTrue_abs, hPull_abs);

   hPull_initE->Write(); hPull_interE->Write(); hPull_totInel->Write(); hPull_abs->Write();


   //==============================================================================================      
   auto legend_initE = new TLegend(0.5,0.65,0.76,0.85);
   legend_initE->AddEntry(hUnfold_initE, "Unfolded Data");
   legend_initE->AddEntry(hMeas_data_initE, "PDSP Data, Run 58XX initE");
   legend_initE->AddEntry(hTrue_initE, "MC True");
   legend_initE->SetTextSize(0.03);

   auto legend_interE = new TLegend(0.5,0.65,0.76,0.85);
   legend_interE->AddEntry(hUnfold_interE, "Unfolded Data");
   legend_interE->AddEntry(hMeas_data_interE, "PDSP Data, Run 58XX interE");
   legend_interE->AddEntry(hTrue_interE, "MC True");
   legend_interE->SetTextSize(0.03);

   auto legend_totInel = new TLegend(0.5,0.65,0.76,0.85);
   legend_totInel->AddEntry(hUnfold_totInel, "Unfolded Data");
   legend_totInel->AddEntry(hMeas_data_totInel, "PDSP Data, Run 58XX totInel");
   legend_totInel->AddEntry(hTrue_totInel, "MC True");
   legend_totInel->SetTextSize(0.03);

   auto legend_abs = new TLegend(0.5,0.65,0.76,0.85);
   legend_abs->AddEntry(hUnfold_abs, "Unfolded Data");
   legend_abs->AddEntry(hMeas_data_abs, "PDSP Data, Run 58XX absorption");
   legend_abs->AddEntry(hTrue_abs, "MC True");
   legend_abs->SetTextSize(0.03);

   //==============================================================================================      
   gStyle->SetOptFit(0);
   gStyle->SetOptStat(0);

   TCanvas* c_initE= new TCanvas("canvas_initE","canvas_initE", 1100, 800);
   gPad->SetGrid(1,1);
   hTrue_initE->SetLineColor(8);
   hTrue_initE->SetMarkerColorAlpha(kBlack,0);
   hTrue_initE->Draw("HIST ");
   hMeas_data_initE->SetMarkerColorAlpha(kBlack,0);
   hMeas_data_initE->SetFillColorAlpha(kBlack,0.2);
   hMeas_data_initE->Draw("SAME HIST");
   hUnfold_initE->Draw("PE1 SAME");
   //hMeas_mc_initE->SetFillColorAlpha(kBlue,0.1);
   //hMeas_mc_initE->SetLineColor(kBlue);
   //hMeas_mc_initE->Draw("SAME HIST");
   legend_initE->Draw();
   c_initE->Write();
   //c_initE->SaveAs("bla.pdf");

   TCanvas* c_interE= new TCanvas("canvas_interE","canvas_interE", 1100, 800);
   gPad->SetGrid(1,1);
   hTrue_interE->SetLineColor(8);
   hTrue_interE->SetMarkerColorAlpha(kBlack,0);
   hTrue_interE->Draw("HIST ");
   hMeas_data_interE->SetMarkerColorAlpha(kBlack,0);
   hMeas_data_interE->SetFillColorAlpha(kBlack,0.2);
   hMeas_data_interE->Draw("SAME HIST");
   hUnfold_interE->Draw("PE1 SAME");
   //hMeas_mc_interE->SetFillColorAlpha(kBlue,0.1);
   //hMeas_mc_interE->SetLineColor(kBlue);
   //hMeas_mc_interE->Draw("SAME HIST");
   legend_interE->Draw();
   c_interE->Write();
   //c_interE->SaveAs("bla.pdf");

   TCanvas* c_totInel= new TCanvas("canvas_totInel","canvas_totInel", 1100, 800);
   gPad->SetGrid(1,1);
   hTrue_totInel->SetLineColor(8);
   hTrue_totInel->SetMarkerColorAlpha(kBlack,0);
   hTrue_totInel->Draw("HIST ");
   hMeas_data_totInel->SetMarkerColorAlpha(kBlack,0);
   hMeas_data_totInel->SetFillColorAlpha(kBlack,0.2);
   hMeas_data_totInel->Draw("SAME HIST");
   hUnfold_totInel->Draw("PE1 SAME");
   //hMeas_mc_totInel->SetFillColorAlpha(kBlue,0.1);
   //hMeas_mc_totInel->SetLineColor(kBlue);
   //hMeas_mc_totInel->Draw("SAME HIST");
   legend_totInel->Draw();
   c_totInel->Write();
   //c_totInel->SaveAs("bla.pdf");

   TCanvas* c_abs= new TCanvas("canvas_abs","canvas_abs", 1100, 800);
   gPad->SetGrid(1,1);
   hTrue_abs->SetLineColor(8);
   hTrue_abs->SetMarkerColorAlpha(kBlack,0);
   hTrue_abs->Draw("HIST ");
   hMeas_data_abs->SetMarkerColorAlpha(kBlack,0);
   hMeas_data_abs->SetFillColorAlpha(kBlack,0.2);
   hMeas_data_abs->Draw("SAME HIST");
   hUnfold_abs->Draw("PE1 SAME");
   //hMeas_mc_abs->SetFillColorAlpha(kBlue,0.1);
   //hMeas_mc_abs->SetLineColor(kBlue);
   //hMeas_mc_abs->Draw("SAME HIST");
   legend_abs->Draw();
   c_abs->Write();
   //c_abs->SaveAs("bla.pdf");

   //Build the Incident histos in order to compare

   TH1D* hTrue_incident= new TH1D ("trueMC_incident", "MC Truth incident; Energy [MeV]; Events / 50 MeV",    nBin_int, eEnd, eStart);
   TH1D* hMeas_mc_incident= new TH1D ("meas_mc_incident", "Measured MC incident; Energy [MeV]; Events / 50 MeV", nBin_int, eEnd, eStart);
   TH1D* hMeas_data_incident= new TH1D ("meas_data_incident", "Measured Data incident; Energy [MeV]; Events / 50 MeV", nBin_int, eEnd, eStart);
   TH1D* hUnfold_incident= new TH1D ("hUnfold_incident", "Unfolded incident; Energy [MeV]; Events / 50 MeV", nBin_int, eEnd, eStart);

   build_incidentHist(hTrue_initE, hTrue_interE, hTrue_incident);
   build_incidentHist(hMeas_mc_initE, hMeas_mc_interE, hMeas_mc_incident);
   build_incidentHist(hUnfold_initE, hUnfold_interE, hUnfold_incident);
   build_incidentHist(hMeas_data_initE, hMeas_data_interE, hMeas_data_incident);

   auto legend_incident = new TLegend(0.5,0.65,0.76,0.85);
   legend_incident->AddEntry(hUnfold_incident, "Unfolded Data");
   legend_incident->AddEntry(hMeas_data_incident, "PDSP Data, Run 58XX incident");
   legend_incident->AddEntry(hTrue_incident, "MC True");
   legend_incident->SetTextSize(0.03);

   //==============================================================================================      
   //build_incidentHist(hTrue_initE, hTrue_interE, hTrue_incident);

   //USE TRUE-INITE to build INCIDENT E histo for MEAS and UNFOLD
   //hTrue_initE->Scale( hMeas_mc_interE->Integral() / hTrue_initE->Integral() );
   //build_incidentHist(hTrue_initE, hMeas_mc_interE, hMeas_mc_incident);


   //hTrue_initE->Scale( hUnfold_interE->Integral() / hTrue_initE->Integral() );
   //build_incidentHist(hTrue_initE, hUnfold_interE, hUnfold_incident);

   TCanvas* c_incident= new TCanvas("canvas_incident","canvas_incident", 1100, 800);
   gPad->SetGrid(1,1);
   hTrue_incident->SetLineColor(8);
   hTrue_incident->SetMarkerColorAlpha(kBlack,0);
   hTrue_incident->Draw("HIST ");
   hMeas_data_incident->SetMarkerColorAlpha(kBlack,0);
   hMeas_data_incident->SetFillColorAlpha(kBlack,0.2);
   hMeas_data_incident->Draw("SAME HIST");
   hUnfold_incident->Draw("PE1 SAME");
   //hMeas_mc_incident->SetFillColorAlpha(kBlue,0.1);
   //hMeas_mc_incident->SetLineColor(kBlue);
   //hMeas_mc_incident->Draw("SAME HIST");
   legend_incident->Draw();
   c_incident->Write();
   //c_incident->SaveAs("bla.pdf");

   auto* R_initE = response_initE.HresponseNoOverflow();
   auto* c_response_initE = new TCanvas();
   R_initE->SetStats(0);
   R_initE->Draw("colz");
   R_initE->GetXaxis()->SetNdivisions(1020);
   R_initE->GetYaxis()->SetNdivisions(1020);
   R_initE->SetTitle("; Reco Pion InitE [MeV]; True Pion InitE [MeV]");
   gPad->SetGrid(1,1);
   c_response_initE->Draw();

   auto* R_interE = response_interE.HresponseNoOverflow();
   auto* c_response_interE = new TCanvas();
   R_interE->SetStats(0);
   R_interE->Draw("colz");
   R_interE->GetXaxis()->SetNdivisions(1020);
   R_interE->GetYaxis()->SetNdivisions(1020);
   R_interE->SetTitle("; Reco Pion InterE [MeV]; True Pion InterE [MeV]");
   gPad->SetGrid(1,1);
   c_response_interE->Draw();

   auto* R_totInel = response_totInel.HresponseNoOverflow();
   auto* c_response_totInel = new TCanvas();
   R_totInel->SetStats(0);
   R_totInel->Draw("colz");
   R_totInel->GetXaxis()->SetNdivisions(1020);
   R_totInel->GetYaxis()->SetNdivisions(1020);
   R_totInel->SetTitle("Response Matrix Pion Total Inelastic Energy Distribution");
   gPad->SetGrid(1,1);
   c_response_totInel->Draw();

   auto* R_abs = response_abs.HresponseNoOverflow();
   auto* c_response_abs = new TCanvas();
   R_abs->SetStats(0);
   R_abs->Draw("colz");
   R_abs->GetXaxis()->SetNdivisions(1020);
   R_abs->GetYaxis()->SetNdivisions(1020);
   R_abs->SetTitle("; Reco Abs InterE [MeV]; True Abs InterE [MeV]");
   gPad->SetGrid(1,1);
   c_response_abs->Draw();

   if(doXS){
      //access Jakes GeantFile in folder
      TFile f1("exclusive_xsec.root");
      TGraph *totInel_KE = (TGraph*)f1.Get("total_inel_KE");
      TGraph *abs_KE = (TGraph*)f1.Get("abs_KE");
      TGraph *cex_KE = (TGraph*)f1.Get("cex_KE");
      f1.Close();

      output->cd();

      TH1D* h_betheMean_muon = new TH1D("h_betheMean_muon", "Mean Energy Loss", nBin_int, eEnd, eStart);

      //fill histo with Mean dEdX of bin center
      for(int i = 1; i <= nBin_int; i++){
         h_betheMean_muon->SetBinContent(i , betheBloch( eEnd + (i - 0.5)*bin_size_int  , mass_muon) );
      };


      TH1D* h_unfoldXS_totInel = new TH1D("h_unfoldXS_totInel" ,"DATA Unfold Total Inelastic XS; Kinetic energy [MeV]; #sigma [mbarn]", nBin_int, eEnd, eStart);

      do_XS_log(  h_unfoldXS_totInel, hUnfold_totInel, hUnfold_incident, h_betheMean_muon );   
      do_XS_log_binomial_error( h_unfoldXS_totInel, hUnfold_totInel, hUnfold_incident, h_betheMean_muon );   

      TH1D* h_TrueXS_totInel = new TH1D("h_TrueXS_totInel" ,"True Total Inelastic XS; Kinetic energy [MeV]; #sigma [mbarn]", nBin_int, eEnd, eStart);

      do_XS_log(  h_TrueXS_totInel, hTrue_totInel, hTrue_incident, h_betheMean_muon );   
      do_XS_log_binomial_error( h_TrueXS_totInel, hTrue_totInel, hTrue_incident, h_betheMean_muon );   

      TH1D* h_rawXS_totInel = new TH1D("h_rawXS_totInel" ,"DATA Raw totInel XS; Kinetic energy [MeV]; #sigma [mbarn]", nBin_int, eEnd, eStart);

      do_XS_log(  h_rawXS_totInel, hMeas_data_totInel, hMeas_data_incident, h_betheMean_muon );   
      do_XS_log_binomial_error( h_rawXS_totInel, hMeas_data_totInel, hMeas_data_incident, h_betheMean_muon );   


      TCanvas *c_unfold_totInel = new TCanvas("c_unfold_totInel", "c_unfold_totInel", 1100, 800);
      gPad->SetGrid(1,1);
      h_unfoldXS_totInel->GetXaxis()->SetRangeUser(400,900);
      h_unfoldXS_totInel->GetXaxis()->SetNdivisions(1020);
      h_unfoldXS_totInel->GetYaxis()->SetNdivisions(1020);

      h_TrueXS_totInel->GetXaxis()->SetRangeUser(300,1000);

      totInel_KE->SetTitle( "Geant Total Inelasitc XS; Kinetic energy [MeV]; #sigma [mbarn]");
      totInel_KE->GetXaxis()->SetRangeUser(eEnd, eStart);
      totInel_KE->SetLineColor(kRed);
      totInel_KE->SetLineWidth(3);
      totInel_KE->Draw("AC");
      h_unfoldXS_totInel->SetMarkerSize(0.7);
      h_unfoldXS_totInel->Draw("PE0 SAME");
      //h_TrueXS_totInel->SetMarkerColor(8);
      //h_TrueXS_totInel->SetMarkerStyle(22);
      //h_TrueXS_totInel->Draw("PE0 SAME");
      //h_rawXS_totInel->SetMarkerColorAlpha(kOrange, 0.3);
      //h_rawXS_totInel->Draw("P SAME");
      c_unfold_totInel->BuildLegend();

      c_unfold_totInel->Write();

      TH1D* h_unfoldXS_abs = new TH1D("h_unfoldXS_abs" ,"DATA Unfold Absorption XS; Kinetic energy [MeV]; #sigma [mbarn]", nBin_int, eEnd, eStart);

      do_XS_log(  h_unfoldXS_abs, hUnfold_abs, hUnfold_incident, h_betheMean_muon );   
      do_XS_log_binomial_error( h_unfoldXS_abs, hUnfold_abs, hUnfold_incident, h_betheMean_muon );   
      
      h_unfoldXS_abs->Write();

      TH1D* h_TrueXS_abs = new TH1D("h_TrueXS_abs" ,"True Absorption XS; Kinetic energy [MeV]; #sigma [mbarn]", nBin_int, eEnd, eStart);

      do_XS_log(  h_TrueXS_abs, hTrue_abs, hTrue_incident, h_betheMean_muon );   
      do_XS_log_binomial_error( h_TrueXS_abs, hTrue_abs, hTrue_incident, h_betheMean_muon );   

      TH1D* h_rawXS_abs = new TH1D("h_rawXS_abs" ,"DATA Raw abs XS; Kinetic energy [MeV]; #sigma [mbarn]", nBin_int, eEnd, eStart);

      do_XS_log(  h_rawXS_abs, hMeas_data_abs, hMeas_data_incident, h_betheMean_muon );   
      do_XS_log_binomial_error( h_rawXS_abs, hMeas_data_abs, hMeas_data_incident, h_betheMean_muon );   

      TCanvas *c_unfold_abs = new TCanvas("c_unfold_abs", "c_unfold_abs", 1100, 800);
      gPad->SetGrid(1,1);
      h_unfoldXS_abs->GetXaxis()->SetRangeUser(450,1100);
      h_unfoldXS_abs->GetXaxis()->SetNdivisions(1020);
      h_unfoldXS_abs->GetYaxis()->SetNdivisions(1020);

      h_TrueXS_abs->GetXaxis()->SetRangeUser(450,1100);

      abs_KE->SetTitle( "Geant Absorption XS; Kinetic energy [MeV]; #sigma [mbarn]");
      abs_KE->GetXaxis()->SetRangeUser(0, eStart);
      abs_KE->SetLineColor(kBlue);
      abs_KE->SetLineWidth(3);
      abs_KE->Draw("AC");
      h_unfoldXS_abs->SetMarkerSize(0.7);
      h_unfoldXS_abs->Draw("PE0 SAME");
      //h_TrueXS_abs->SetMarkerColor(8);
      //h_TrueXS_abs->SetMarkerStyle(22);
      //h_TrueXS_abs->Draw("PE0 SAME");
      //h_rawXS_abs->GetXaxis()->SetRangeUser(400,1100);
      //h_rawXS_abs->SetMarkerColorAlpha(kOrange, 0.3);
      //h_rawXS_abs->Draw("P SAME");
      c_unfold_abs->BuildLegend();

      c_unfold_abs->Write();

   }



}

#ifndef __CINT__
int main () { eSliceMethod_unfoldBayes_XS_dataMC(mcPath, dataPath); return 0; }  // Main program when run stand-alone
#endif




