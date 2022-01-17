#include "TCanvas.h"
#include "TROOT.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TGraphMultiErrors.h"
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
#include "eSlice.h"
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


int thesisPlot_systematic_xs(){

      TFile f2("unfold_xs_data_50MeV.root", "READ");
      TH1D *hNominalUnfold_initE = (TH1D*)f2.Get("hUnfold_initE");
      TH1D *hNominalUnfold_interE = (TH1D*)f2.Get("hUnfold_interE");
      TH1D *hNominalUnfold_abs = (TH1D*)f2.Get("hUnfold_abs");
      TH1D *hNominalUnfold_xs = (TH1D*)f2.Get("h_unfoldXS_abs");

      hNominalUnfold_initE->SetDirectory(0); hNominalUnfold_interE->SetDirectory(0); hNominalUnfold_abs->SetDirectory(0);
      f2.Close();

      TFile f1("unfoldSystematic_50MeV.root", "READ");
      TH1D *hSyst_plusEnergy_initE = (TH1D*)f1.Get("hSystUnfold_plusEnergy_initE");
      TH1D *hSyst_plusEnergy_interE = (TH1D*)f1.Get("hSystUnfold_plusEnergy_interE");
      TH1D *hSyst_plusEnergy_abs = (TH1D*)f1.Get("hSystUnfold_plusEnergy_abs");
      TH1D *hSyst_minusEnergy_initE = (TH1D*)f1.Get("hSystUnfold_minusEnergy_initE");
      TH1D *hSyst_minusEnergy_interE = (TH1D*)f1.Get("hSystUnfold_minusEnergy_interE");
      TH1D *hSyst_minusEnergy_abs = (TH1D*)f1.Get("hSystUnfold_minusEnergy_abs");
      
      TH1D *hSyst_plusAbsEvent_initE = (TH1D*)f1.Get("hSystUnfold_plusAbsEvent20_initE");
      TH1D *hSyst_plusAbsEvent_interE = (TH1D*)f1.Get("hSystUnfold_plusAbsEvent20_interE");
      TH1D *hSyst_plusAbsEvent_abs = (TH1D*)f1.Get("hSystUnfold_plusAbsEvent20_abs");
      TH1D *hSyst_minusAbsEvent_initE = (TH1D*)f1.Get("hSystUnfold_minusAbsEvent20_initE");
      TH1D *hSyst_minusAbsEvent_interE = (TH1D*)f1.Get("hSystUnfold_minusAbsEvent20_interE");
      TH1D *hSyst_minusAbsEvent_abs = (TH1D*)f1.Get("hSystUnfold_minusAbsEvent20_abs");
      
      TH1D *hSyst_plusAbsPurity_initE = (TH1D*)f1.Get("hSystUnfold_plusAbsPurity_initE");
      TH1D *hSyst_plusAbsPurity_interE = (TH1D*)f1.Get("hSystUnfold_plusAbsPurity_interE");
      TH1D *hSyst_plusAbsPurity_abs = (TH1D*)f1.Get("hSystUnfold_plusAbsPurity_abs");
      TH1D *hSyst_minusAbsPurity_initE = (TH1D*)f1.Get("hSystUnfold_minusAbsPurity_initE");
      TH1D *hSyst_minusAbsPurity_interE = (TH1D*)f1.Get("hSystUnfold_minusAbsPurity_interE");
      TH1D *hSyst_minusAbsPurity_abs = (TH1D*)f1.Get("hSystUnfold_minusAbsPurity_abs");

      TH1D *hSyst_plusAbsEfficiency_initE = (TH1D*)f1.Get("hSystUnfold_plusAbsEfficiency_initE");
      TH1D *hSyst_plusAbsEfficiency_interE = (TH1D*)f1.Get("hSystUnfold_plusAbsEfficiency_interE");
      TH1D *hSyst_plusAbsEfficiency_abs = (TH1D*)f1.Get("hSystUnfold_plusAbsEfficiency_abs");
      TH1D *hSyst_minusAbsEfficiency_initE = (TH1D*)f1.Get("hSystUnfold_minusAbsEfficiency_initE");
      TH1D *hSyst_minusAbsEfficiency_interE = (TH1D*)f1.Get("hSystUnfold_minusAbsEfficiency_interE");
      TH1D *hSyst_minusAbsEfficiency_abs = (TH1D*)f1.Get("hSystUnfold_minusAbsEfficiency_abs");
      
      hSyst_plusEnergy_initE->SetDirectory(0);        hSyst_plusAbsEvent_initE->SetDirectory(0); 
      hSyst_plusEnergy_interE->SetDirectory(0);       hSyst_plusAbsEvent_interE->SetDirectory(0);
      hSyst_plusEnergy_abs->SetDirectory(0);          hSyst_plusAbsEvent_abs->SetDirectory(0);
      hSyst_minusEnergy_initE->SetDirectory(0);       hSyst_minusAbsEvent_initE->SetDirectory(0);
      hSyst_minusEnergy_interE->SetDirectory(0);      hSyst_minusAbsEvent_interE->SetDirectory(0);
      hSyst_minusEnergy_abs->SetDirectory(0);         hSyst_minusAbsEvent_abs->SetDirectory(0);      
      
      hSyst_plusAbsPurity_initE->SetDirectory(0);     hSyst_plusAbsEfficiency_initE->SetDirectory(0);
      hSyst_plusAbsPurity_interE->SetDirectory(0);    hSyst_plusAbsEfficiency_interE->SetDirectory(0);
      hSyst_plusAbsPurity_abs->SetDirectory(0);       hSyst_plusAbsEfficiency_abs->SetDirectory(0);
      hSyst_minusAbsPurity_initE->SetDirectory(0);    hSyst_minusAbsEfficiency_initE->SetDirectory(0);
      hSyst_minusAbsPurity_interE->SetDirectory(0);   hSyst_minusAbsEfficiency_interE->SetDirectory(0);
      hSyst_minusAbsPurity_abs->SetDirectory(0);      hSyst_minusAbsEfficiency_abs->SetDirectory(0);

      f1.Close();
      //
      //Build the Incidents
      
      TH1D* hNominalUnfold_incident= new TH1D ("hNominalUnfold_incident", "; Energy [MeV]; Events / 50 MeV", nBin_int, eEnd, eStart);
      TH1D* hSyst_plusEnergy_incident= new TH1D ("hSyst_plusEnergy_incident", "; Energy [MeV]; Events / 50 MeV", nBin_int, eEnd, eStart);
      TH1D* hSyst_minusEnergy_incident= new TH1D ("hSyst_minusEnergy_incident", "; Energy [MeV]; Events / 50 MeV", nBin_int, eEnd, eStart);
      TH1D* hSyst_plusAbsEvent_incident= new TH1D ("hSyst_plusAbsEvent_incident", "; Energy [MeV]; Events / 50 MeV", nBin_int, eEnd, eStart);
      TH1D* hSyst_minusAbsEvent_incident= new TH1D ("hSyst_minusAbsEvent_incident", "; Energy [MeV]; Events / 50 MeV", nBin_int, eEnd, eStart);
      TH1D* hSyst_plusAbsPurity_incident= new TH1D ("hSyst_plusAbsPurity_incident", "; Energy [MeV]; Events / 50 MeV", nBin_int, eEnd, eStart);
      TH1D* hSyst_minusAbsPurity_incident= new TH1D ("hSyst_minusAbsPurity_incident", "; Energy [MeV]; Events / 50 MeV", nBin_int, eEnd, eStart);
      TH1D* hSyst_plusAbsEfficiency_incident= new TH1D ("hSyst_plusAbsEfficiency_incident", "; Energy [MeV]; Events / 50 MeV", nBin_int, eEnd, eStart);
      TH1D* hSyst_minusAbsEfficiency_incident= new TH1D ("hSyst_minusAbsEfficiency_incident", "; Energy [MeV]; Events / 50 MeV", nBin_int, eEnd, eStart);

      build_incidentHist( hNominalUnfold_initE , hNominalUnfold_interE , hNominalUnfold_incident );
      build_incidentHist( hSyst_plusEnergy_initE , hSyst_plusEnergy_interE , hSyst_plusEnergy_incident );
      build_incidentHist( hSyst_minusEnergy_initE , hSyst_minusEnergy_interE , hSyst_minusEnergy_incident );
      build_incidentHist( hSyst_plusAbsEvent_initE , hSyst_plusAbsEvent_interE , hSyst_plusAbsEvent_incident );
      build_incidentHist( hSyst_minusAbsEvent_initE , hSyst_minusAbsEvent_interE , hSyst_minusAbsEvent_incident );
      build_incidentHist( hSyst_plusAbsPurity_initE , hSyst_plusAbsPurity_interE , hSyst_plusAbsPurity_incident );
      build_incidentHist( hSyst_minusAbsPurity_initE , hSyst_minusAbsPurity_interE , hSyst_minusAbsPurity_incident );
      build_incidentHist( hSyst_plusAbsEfficiency_initE , hSyst_plusAbsEfficiency_interE , hSyst_plusAbsEfficiency_incident );
      build_incidentHist( hSyst_minusAbsEfficiency_initE , hSyst_minusAbsEfficiency_interE , hSyst_minusAbsEfficiency_incident ); 
      
      //Bethe Bloch
      TH1D* h_betheMean_muon = new TH1D("h_betheMean_muon", "Mean Energy Loss", nBin_int, eEnd, eStart);
      for(int i = 1; i <= nBin_int; i++){
         h_betheMean_muon->SetBinContent(i , betheBloch( eEnd + (i - 0.5)*bin_size_int  , mass_muon) );
      };
      
      //================================================================================

      TH1D* hXS_nominalUnfold = new TH1D("hXS_nominalUnfold" ,"; Kinetic energy [MeV]; #sigma [mbarn]", nBin_int, eEnd, eStart);
      TH1D* hXS_syst_plusEnergy = new TH1D("hXS_syst_plusEnergy" ,"; Kinetic energy [MeV]; #sigma [mbarn]", nBin_int, eEnd, eStart);
      TH1D* hXS_syst_minusEnergy = new TH1D("hXS_syst_minusEnergy" ,"; Kinetic energy [MeV]; #sigma [mbarn]", nBin_int, eEnd, eStart);
      TH1D* hXS_syst_plusAbsEvent = new TH1D("hXS_syst_plusAbsEvent" ,"; Kinetic energy [MeV]; #sigma [mbarn]", nBin_int, eEnd, eStart);
      TH1D* hXS_syst_minusAbsEvent = new TH1D("hXS_syst_minusAbsEvent" ,"; Kinetic energy [MeV]; #sigma [mbarn]", nBin_int, eEnd, eStart);
      TH1D* hXS_syst_plusAbsPurity = new TH1D("hXS_syst_plusAbsPurity" ,"; Kinetic energy [MeV]; #sigma [mbarn]", nBin_int, eEnd, eStart);
      TH1D* hXS_syst_minusAbsPurity = new TH1D("hXS_syst_minusAbsPurity" ,"; Kinetic energy [MeV]; #sigma [mbarn]", nBin_int, eEnd, eStart);
      TH1D* hXS_syst_plusAbsEfficiency = new TH1D("hXS_syst_plusAbsEfficiency" ,"; Kinetic energy [MeV]; #sigma [mbarn]", nBin_int, eEnd, eStart);
      TH1D* hXS_syst_minusAbsEfficiency = new TH1D("hXS_syst_minusAbsEfficiency" ,"; Kinetic energy [MeV]; #sigma [mbarn]", nBin_int, eEnd, eStart);
      
      do_XS_log(  hXS_nominalUnfold, hNominalUnfold_abs, hNominalUnfold_incident, h_betheMean_muon );
      do_XS_log_binomial_error(  hXS_nominalUnfold, hNominalUnfold_abs, hNominalUnfold_incident, h_betheMean_muon );
      
      do_XS_log(  hXS_syst_plusEnergy, hSyst_plusEnergy_abs, hSyst_plusEnergy_incident, h_betheMean_muon );
      do_XS_log(  hXS_syst_minusEnergy, hSyst_minusEnergy_abs, hSyst_minusEnergy_incident, h_betheMean_muon );
      do_XS_log(  hXS_syst_plusAbsEvent, hSyst_plusAbsEvent_abs, hSyst_plusAbsEvent_incident, h_betheMean_muon );
      do_XS_log(  hXS_syst_minusAbsEvent, hSyst_minusAbsEvent_abs, hSyst_minusAbsEvent_incident, h_betheMean_muon );
      do_XS_log(  hXS_syst_plusAbsPurity, hSyst_plusAbsPurity_abs, hSyst_plusAbsPurity_incident, h_betheMean_muon );
      do_XS_log(  hXS_syst_minusAbsPurity, hSyst_minusAbsPurity_abs, hSyst_minusAbsPurity_incident, h_betheMean_muon );
      do_XS_log(  hXS_syst_plusAbsEfficiency, hSyst_plusAbsEfficiency_abs, hSyst_plusAbsEfficiency_incident, h_betheMean_muon );
      do_XS_log(  hXS_syst_minusAbsEfficiency, hSyst_minusAbsEfficiency_abs, hSyst_minusAbsEfficiency_incident, h_betheMean_muon );
      
      //================================================================================
      // Errors!
      // high --> higher XS, err[0] at 450 MeV, err[nBin_int -1] at 1000MeV
      double err_energy_high[nBin_int], err_energy_low[nBin_int];
      double err_eff_high[nBin_int], err_eff_low[nBin_int];
      double err_pur_high[nBin_int], err_pur_low[nBin_int];
      //symmetric error, xs variation, xs calc error
      double err_xs_var_high[nBin_int], err_xs_var_low[nBin_int];
      double xs_calc[nBin_int];
      double err_xs_calc[nBin_int];
      double err_sqrtSumSquare_high[nBin_int], err_sqrtSumSquare_low[nBin_int];
      double err_justSyst_high[nBin_int], err_justSyst_low[nBin_int];

      for(int i=0; i < nBin_int; i++){

         double nominal = hXS_nominalUnfold->GetBinContent(i+1);
         xs_calc[i] = nominal;

         //std::cout << "Bin = " << i+ 1 << "Value = " << nominal << std::endl;

         err_energy_high[i] = hXS_syst_plusEnergy->GetBinContent(i+1) - nominal;
         err_energy_low[i] = nominal - hXS_syst_minusEnergy->GetBinContent(i+1);
         err_xs_var_high[i] = hXS_syst_plusAbsEvent->GetBinContent(i+1) - nominal;
         err_xs_var_low[i] = nominal - hXS_syst_minusAbsEvent->GetBinContent(i+1);
         err_pur_high[i] = hXS_syst_plusAbsPurity->GetBinContent(i+1) - nominal;
         err_pur_low[i] = nominal - hXS_syst_minusAbsPurity->GetBinContent(i+1);
         err_eff_high[i] = hXS_syst_plusAbsEfficiency->GetBinContent(i+1) - nominal;
         err_eff_low[i] = nominal - hXS_syst_minusAbsEfficiency->GetBinContent(i+1);
         err_xs_calc[i] = hXS_nominalUnfold->GetBinError(i+1);

      }

      for(int i=0; i <nBin_int-1; i++){ //nBin_int -1 because last bin 950-1000 does have xs=0 and errors are 0 / nan

         double highSquare = 0, lowSquare = 0;

         if(err_energy_high[i] > 0) highSquare = highSquare + pow( err_energy_high[i], 2);
         else lowSquare = lowSquare + pow( err_energy_high[i], 2); //if error is negative contribution
         if(err_energy_low[i] > 0) lowSquare = lowSquare + pow( err_energy_low[i], 2);
         else highSquare = highSquare + pow( err_energy_low[i], 2); //if error is negative contribution
         
         if(err_xs_var_high[i] > 0) highSquare = highSquare + pow( err_xs_var_high[i], 2);
         else lowSquare = lowSquare + pow( err_xs_var_high[i], 2); //if error is negative contribution
         if(err_xs_var_low[i] > 0) lowSquare = lowSquare + pow( err_xs_var_low[i], 2);
         else highSquare = highSquare + pow( err_xs_var_low[i], 2); //if error is negative contribution

         if(err_pur_high[i] > 0) highSquare = highSquare + pow( err_pur_high[i], 2);
         else lowSquare = lowSquare + pow( err_pur_high[i], 2); //if error is negative contribution
         if(err_pur_low[i] > 0) lowSquare = lowSquare + pow( err_pur_low[i], 2);
         else highSquare = highSquare + pow( err_pur_low[i], 2); //if error is negative contribution
         
         if(err_eff_high[i] > 0) highSquare = highSquare + pow( err_eff_high[i], 2);
         else lowSquare = lowSquare + pow( err_eff_high[i], 2); //if error is negative contribution
         if(err_eff_low[i] > 0) lowSquare = lowSquare + pow( err_eff_low[i], 2);
         else highSquare = highSquare + pow( err_eff_low[i], 2); //if error is negative contribution

         //highSquare = highSquare + pow(err_xs_calc[i],2);
         //lowSquare = lowSquare + pow(err_xs_calc[i],2);

         err_sqrtSumSquare_high[i] = sqrt( highSquare );
         err_sqrtSumSquare_low[i] = sqrt( lowSquare );
         err_justSyst_high[i] = sqrt( highSquare );
         err_justSyst_low[i] = sqrt( lowSquare );
         
         std::cout << "Bin 450 + 50*" << i << "  ============" <<std::endl;
         std::cout << "XS central value = " <<  xs_calc[i] << std::endl;
         std::cout << "Plus Energy Error = " << err_energy_high[i] << std::endl;
         std::cout << "Minus Energy Error = " << err_energy_low[i] << std::endl;
         std::cout << "Plus XS Error = " << err_xs_var_high[i] << std::endl;
         std::cout << "Minus XS Error = " << err_xs_var_low[i] << std::endl;
         std::cout << "Plus Purity Error = " << err_pur_high[i] << std::endl;
         std::cout << "Minus Purity Error = " << err_pur_low[i] << std::endl;
         std::cout << "Plus Efficiency Error = " << err_eff_high[i] << std::endl;
         std::cout << "Minus Efficiency Error = " << err_eff_low[i] << std::endl;
         std::cout << "XS calculation Error = " << err_xs_calc[i] << std::endl;
         std::cout << "============================================" << std::endl;
         std::cout << "Total Error High = " << err_sqrtSumSquare_high[i] << std::endl;
         std::cout << "Total Error Low = " << err_sqrtSumSquare_low[i] << std::endl;
         std::cout << "============================================" << std::endl;


      }

      //adding stat error to syst error
      for(int i=0; i < nBin_int -1; i++){

         err_sqrtSumSquare_high[i] = sqrt( pow(err_sqrtSumSquare_high[i], 2) + pow(err_xs_calc[i],2));
         err_sqrtSumSquare_low[i] = sqrt( pow(err_sqrtSumSquare_low[i], 2) + pow(err_xs_calc[i],2));
         
         std::cout << "Bin 450 + 50*" << i << "  ============" <<std::endl;
         std::cout << "============================================" << std::endl;
         std::cout << "Total Error High = " << err_sqrtSumSquare_high[i] << std::endl;
         std::cout << "Total Error Low = " << err_sqrtSumSquare_low[i] << std::endl;
         std::cout << "============================================" << std::endl;

      }

      double kineticEnergy[nBin_int], pseudoErr[nBin_int];
      for(int i = 0; i < nBin_int; i++){ 
         kineticEnergy[i] = 475 + 50*i;
         pseudoErr[i] = 25;
      }

      auto xsAsymmErr = new TGraphAsymmErrors(nBin_int-1, kineticEnergy, xs_calc, pseudoErr, pseudoErr, err_sqrtSumSquare_low, err_sqrtSumSquare_high);
      xsAsymmErr->SetFillColor(kRed);
      xsAsymmErr->SetMarkerColor(kRed);
      xsAsymmErr->SetLineColor(kRed);
      xsAsymmErr->SetFillStyle(3002);

      //auto xsMultiError = new TGraphMultiErrors(nBin_int-1, kineticEnergy, xs_calc, pseudoErr, pseudoErr, err_xs_calc, err_xs_calc ); //errors on y are stat
      //xsMultiError->AddYError(err_sqrtSumSquare_low, err_sqrtSumSquare_high);

      auto xsErr_stat = new TGraphAsymmErrors(nBin_int-1, kineticEnergy, xs_calc, 0, 0, err_xs_calc, err_xs_calc);
      xsErr_stat->SetMarkerColorAlpha(kRed,0);
      xsErr_stat->SetLineColor(kBlack);
      xsErr_stat->SetLineWidth(2);
      //xsErr_stat->SetFillStyle(3002);

      auto xsErr_syst = new TGraphAsymmErrors(nBin_int-1, kineticEnergy, xs_calc, 0, 0, err_justSyst_low, err_justSyst_high);
      //xsErr_syst->SetFillColor(kRed);
      xsErr_syst->SetMarkerColorAlpha(kRed,0);
      xsErr_syst->SetLineColor(kBlue);
      //xsErr_syst->SetFillStyle(3002);

      //================================================================================
      //LADS results
      double lads_xs[5] = {180,320,350,280,220};
      double lads_xs_err[5] = {40,60,40,30,10}; //from jake cm 2021 presentation. LADS data
      double lads_energy[5] = {70,118,162,230,330};
      double lads_energy_err[5] = {0,0,0,0,0};

      TGraphErrors *lads_graph = new TGraphErrors(5,lads_energy,lads_xs, lads_energy_err, lads_xs_err);
      //================================================================================

      //access Jakes GeantFile in folder
      TFile f3("exclusive_xsec.root");
      TGraph *abs_KE = (TGraph*)f3.Get("abs_KE");
      f3.Close();
   
      string output_name;
      output_name = "systematic_XS_" + std::to_string((int) bin_size_int) + "MeV.root";
      TFile *output = new TFile( output_name.c_str() , "RECREATE"); //maybe save with binning?
      output->cd();

      hXS_nominalUnfold->Write();
      hXS_syst_plusEnergy->Write();
      hXS_syst_minusEnergy->Write();
      hXS_syst_plusAbsEvent->Write();
      hXS_syst_minusAbsEvent->Write();
      hXS_syst_plusAbsPurity->Write();
      hXS_syst_minusAbsPurity->Write();
      hXS_syst_plusAbsEfficiency->Write(); 
      hXS_syst_minusAbsEfficiency->Write();

      hXS_nominalUnfold->SetMarkerColor(kBlack);                  hXS_nominalUnfold->SetMarkerStyle(20);
      hXS_syst_plusEnergy->SetMarkerColor(kRed + 3);              hXS_syst_plusEnergy->SetMarkerStyle(54);
      hXS_syst_minusEnergy->SetMarkerColor(kRed - 3);             hXS_syst_minusEnergy->SetMarkerStyle(54);
      hXS_syst_plusAbsEvent->SetMarkerColor(kBlue + 3);           hXS_syst_plusAbsEvent->SetMarkerStyle(57);
      hXS_syst_minusAbsEvent->SetMarkerColor(kBlue - 9);          hXS_syst_minusAbsEvent->SetMarkerStyle(57);
      hXS_syst_plusAbsPurity->SetMarkerColor(kGreen + 3);         hXS_syst_plusAbsPurity->SetMarkerStyle(63);
      hXS_syst_minusAbsPurity->SetMarkerColor(kGreen -6);         hXS_syst_minusAbsPurity->SetMarkerStyle(63);
      hXS_syst_plusAbsEfficiency->SetMarkerColor(kYellow + 3);    hXS_syst_plusAbsEfficiency->SetMarkerStyle(66); 
      hXS_syst_minusAbsEfficiency->SetMarkerColor(kYellow - 6);   hXS_syst_minusAbsEfficiency->SetMarkerStyle(66);
      
      hXS_nominalUnfold->SetLineColorAlpha(kBlack,0);                      hXS_nominalUnfold->SetMarkerSize(1.5);
      hXS_syst_plusEnergy->SetLineColorAlpha(kBlack,0);                    hXS_syst_plusEnergy->SetMarkerSize(1.5);
      hXS_syst_minusEnergy->SetLineColorAlpha(kBlack,0);                   hXS_syst_minusEnergy->SetMarkerSize(1.5);
      hXS_syst_plusAbsEvent->SetLineColorAlpha(kBlack,0);                  hXS_syst_plusAbsEvent->SetMarkerSize(1.5);
      hXS_syst_minusAbsEvent->SetLineColorAlpha(kBlack,0);                 hXS_syst_minusAbsEvent->SetMarkerSize(1.5);
      hXS_syst_plusAbsPurity->SetLineColorAlpha(kBlack,0);                 hXS_syst_plusAbsPurity->SetMarkerSize(1.5);
      hXS_syst_minusAbsPurity->SetLineColorAlpha(kBlack,0);                hXS_syst_minusAbsPurity->SetMarkerSize(1.5);
      hXS_syst_plusAbsEfficiency->SetLineColorAlpha(kBlack,0);             hXS_syst_plusAbsEfficiency->SetMarkerSize(1.5); 
      hXS_syst_minusAbsEfficiency->SetLineColorAlpha(kBlack,0);            hXS_syst_minusAbsEfficiency->SetMarkerSize(1.5);

      auto legend = new TLegend();
      legend->AddEntry(abs_KE, "Geant Absorption Cross-Section"); legend->AddEntry(hXS_nominalUnfold, "Data unfold nominal");
      legend->AddEntry(hXS_syst_plusEnergy, "Energy +#Delta E");
      legend->AddEntry(hXS_syst_minusEnergy, "Energy -#Delta E");
      legend->AddEntry(hXS_syst_plusAbsEvent, "Simulated cross-section +20%");
      legend->AddEntry(hXS_syst_minusAbsEvent, "Simulated cross-section -20%");
      legend->AddEntry(hXS_syst_plusAbsPurity, "Event selection purity +10%");
      legend->AddEntry(hXS_syst_minusAbsPurity, "Event selection purity -10%");
      legend->AddEntry(hXS_syst_plusAbsEfficiency, "Event selection efficiency +10%");
      legend->AddEntry(hXS_syst_minusAbsEfficiency, "Event selection efficiency -10%");
      legend->SetTextSize(0.03);
 
      auto legend_energy = new TLegend();
      legend_energy->AddEntry(abs_KE, "Simulated Cross-Section"); legend_energy->AddEntry(hXS_nominalUnfold, "Data unfold nominal");
      legend_energy->AddEntry(hXS_syst_plusEnergy, "Energy +#Delta E");
      legend_energy->AddEntry(hXS_syst_minusEnergy, "Energy -#Delta E");
      legend_energy->SetTextSize(0.03);
      
      auto legend_absEvent = new TLegend();
      legend_absEvent->AddEntry(abs_KE, "Simulated Cross-Section"); legend_absEvent->AddEntry(hXS_nominalUnfold, "Data unfold nominal");
      legend_absEvent->AddEntry(hXS_syst_plusAbsEvent, "Geant cross-section +20%");
      legend_absEvent->AddEntry(hXS_syst_minusAbsEvent, "Geant cross-section -20%");
      legend_absEvent->SetTextSize(0.03);
      
      auto legend_absPurity = new TLegend();
      legend_absPurity->AddEntry(abs_KE, "Simulated Cross-Section"); legend_absPurity->AddEntry(hXS_nominalUnfold, "Data unfold nominal");
      legend_absPurity->AddEntry(hXS_syst_plusAbsPurity, "Event selection purity +10%");
      legend_absPurity->AddEntry(hXS_syst_minusAbsPurity, "Event selection purity -10%");
      legend_absPurity->SetTextSize(0.03);

      auto legend_absEfficiency = new TLegend();
      legend_absEfficiency->AddEntry(abs_KE, "Simulated Cross-Section"); legend_absEfficiency->AddEntry(hXS_nominalUnfold, "Data unfold nominal");
      legend_absEfficiency->AddEntry(hXS_syst_plusAbsEfficiency, "Event selection efficiency +10%");
      legend_absEfficiency->AddEntry(hXS_syst_minusAbsEfficiency, "Event selection efficiency -10%");
      legend_absEfficiency->SetTextSize(0.03);
 
      auto legend_absGraph = new TLegend();
      legend_absGraph->AddEntry(xsAsymmErr, "Unfolded Pion Absorption Cross-Section, ProtoDUNE-SP"); 
      legend_absGraph->AddEntry(lads_graph, "LADS Pion Absortpion Data");
      legend_absGraph->AddEntry(abs_KE, "Simulated Absorption Cross-Section"); 
      legend_absGraph->SetTextSize(0.03);

      auto legend_absStatSys = new TLegend();
      legend_absStatSys->AddEntry(xsAsymmErr, "Unfolded Pion Absorption Cross-Section, ProtoDUNE-SP"); 
      legend_absStatSys->AddEntry(xsErr_stat, "#sigma_{stat}"); 
      legend_absStatSys->AddEntry(xsErr_syst, "#sigma_{syst}"); 
      legend_absStatSys->AddEntry(lads_graph, "LADS Pion Absortpion Data");
      legend_absStatSys->AddEntry(abs_KE, "Simulated Absorption Cross-Section"); 
      legend_absStatSys->SetTextSize(0.03);


      gStyle->SetOptStat(0);
      //gStyle->SetErrorX(0);

      TCanvas *c_syst_absAll = new TCanvas("c_syst_absAll", "c_syst_absAll", 1100, 800);
      gPad->SetGrid(1,1);
      abs_KE->SetTitle( "Cross-section variations; Kinetic energy [MeV]; #sigma [mbarn]");
      abs_KE->GetXaxis()->SetRangeUser(0,1000);
      abs_KE->GetXaxis()->SetNdivisions(1020);
      abs_KE->GetYaxis()->SetRangeUser(0,700);
      abs_KE->SetMarkerColorAlpha(kBlue, 0);
      abs_KE->SetLineColor(kBlue);
      abs_KE->SetLineWidth(3);
      abs_KE->Draw("AC");
      hXS_syst_plusEnergy->Draw("HIST P SAME");
      hXS_syst_minusEnergy->Draw("HIST P SAME");
      hXS_syst_plusAbsEvent->Draw("HIST P SAME");
      hXS_syst_minusAbsEvent->Draw("HIST P SAME");
      hXS_syst_plusAbsPurity->Draw("HIST P SAME");
      hXS_syst_minusAbsPurity->Draw("HIST P SAME");
      hXS_syst_plusAbsEfficiency->Draw("HIST P SAME");
      hXS_syst_minusAbsEfficiency->Draw("HIST P SAME");
      hXS_nominalUnfold->Draw("HIST P SAME");
      legend->Draw();
      c_syst_absAll->Update();
      c_syst_absAll->Write();

      TCanvas *c_syst_absEnergy = new TCanvas("c_syst_absEnergy", "c_syst_absEnergy", 1100, 800);
      gPad->SetGrid(1,1);
      //abs_KE->GetXaxis()->SetRangeUser(400,1000);
      //abs_KE->GetYaxis()->SetRangeUser(0,350);
      abs_KE->Draw("AC");
      hXS_nominalUnfold->Draw("HIST P SAME");
      hXS_syst_plusEnergy->Draw("HIST P SAME");
      hXS_syst_minusEnergy->Draw("HIST P SAME");
      legend_energy->Draw();
      c_syst_absEnergy->Update();
      c_syst_absEnergy->Write();

      TCanvas *c_syst_absEvent = new TCanvas("c_syst_absEvent", "c_syst_absEvent", 1100, 800);
      gPad->SetGrid(1,1);
      abs_KE->Draw("AC");
      hXS_syst_plusAbsEvent->Draw("HIST P SAME");
      hXS_syst_minusAbsEvent->Draw("HIST P SAME");
      hXS_nominalUnfold->Draw("HIST P SAME");
      legend_absEvent->Draw();
      c_syst_absEvent->Update();
      c_syst_absEvent->Write();

      TCanvas *c_syst_absPurity = new TCanvas("c_syst_AbsPurity", "c_syst_AbsPurity", 1100, 800);
      gPad->SetGrid(1,1);
      abs_KE->Draw("AC");
      hXS_syst_plusAbsPurity->Draw("HIST P SAME");
      hXS_syst_minusAbsPurity->Draw("HIST P SAME");
      hXS_nominalUnfold->Draw("HIST P SAME");
      legend_absPurity->Draw();
      c_syst_absPurity->Update();
      c_syst_absPurity->Write();

      TCanvas *c_syst_absEfficiency = new TCanvas("c_syst_absEfficiency", "c_syst_absEfficiency", 1100, 800);
      gPad->SetGrid(1,1);
      abs_KE->Draw("AC");
      hXS_syst_plusAbsEfficiency->Draw("HIST P SAME");
      hXS_syst_minusAbsEfficiency->Draw("HIST P SAME");
      hXS_nominalUnfold->Draw("HIST P SAME");
      legend_absEfficiency->Draw();
      c_syst_absEfficiency->Update();
      c_syst_absEfficiency->Write();

      TMultiGraph *mg = new TMultiGraph();
      mg->Add(abs_KE, "AC");
      mg->Add(lads_graph, "A*");
      mg->Add(xsAsymmErr, "A2");
      mg->Add(xsAsymmErr, "AP");
      mg->Add(xsErr_stat, "[]");

      TCanvas *c_syst_final = new TCanvas("c_syst_final", "c_syst_final", 1100, 800);
      gPad->SetGrid(1,1);
      mg->GetXaxis()->SetRangeUser(0,1000);
      mg->GetYaxis()->SetRangeUser(0,700);
      mg->GetXaxis()->SetNdivisions(1020);
      mg->GetXaxis()->SetTitle("Kinetic energy [MeV]");
      mg->GetYaxis()->SetTitle("#sigma [mbarn]");
      mg->Draw("AP");
      //hXS_nominalUnfold->Draw("HIST P SAME");
      legend_absGraph->Draw();


      TMultiGraph *mg_stat_sys = new TMultiGraph();
      mg_stat_sys->Add(abs_KE, "AC");
      mg_stat_sys->Add(lads_graph, "A*");
      mg_stat_sys->Add(xsAsymmErr, "AP2");
      //mg_stat_sys->Add(xsAsymmErr, "AP");
      mg_stat_sys->Add(xsErr_stat, "||");
      mg_stat_sys->Add(xsErr_syst, "[]");

      TCanvas *c_syst_final_statSys = new TCanvas("c_syst_final_statSys", "c_syst_final_statSys", 1100, 800);
      gPad->SetGrid(1,1);
      mg_stat_sys->GetXaxis()->SetRangeUser(0,1000);
      mg_stat_sys->GetYaxis()->SetRangeUser(0,700);
      mg_stat_sys->GetXaxis()->SetNdivisions(1020);
      mg_stat_sys->GetXaxis()->SetTitle("Kinetic energy [MeV]");
      mg_stat_sys->GetYaxis()->SetTitle("#sigma [mbarn]");
      mg_stat_sys->Draw("AP");
      //hXS_nominalUnfold->Draw("HIST P SAME");
      legend_absStatSys->Draw();


   return 0;
}

