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

const string path = "eSliceMethod_Prod4a_mc_1GeV_all_08_22_21.root";
const string dataPath1 = "eSliceMethod_Prod4a_58XX_1GeV_all_08_22_21.root";

int thesisPlot_geantXS(){

      TFile f1("exclusive_xsec.root");
      TGraph *totInel_mom = (TGraph*)f1.Get("total_inel_KE");
      TGraph *abs_mom = (TGraph*)f1.Get("abs_KE");
      TGraph *cex_mom = (TGraph*)f1.Get("cex_KE");
      TGraph *prod_mom = (TGraph*)f1.Get("prod_KE");
      TGraph *inel_mom = (TGraph*)f1.Get("inel_KE");
      f1.Close();

      double lads_xs[5] = {180,320,350,280,220};
      double lads_xs_err[5] = {40,60,40,30,10}; //from jake cm 2021 presentation. LADS data
      double lads_energy[5] = {70,118,162,230,330};
      double lads_energy_err[5] = {0,0,0,0,0};

      TGraphErrors *lads_graph = new TGraphErrors(5,lads_energy,lads_xs, lads_energy_err, lads_xs_err);

      auto legend = new TLegend(0.2,0.75,0.4,0.85);
      legend->AddEntry(lads_graph, "LADS data");
      legend->AddEntry(totInel_mom, "Total inealstic");
      legend->AddEntry(inel_mom, "Inelastic");
      legend->AddEntry(abs_mom, "Absorption");
      legend->AddEntry(cex_mom, "Charge Exchange");
      legend->AddEntry(prod_mom, "Pion Production");

      totInel_mom->SetLineWidth(3);
      totInel_mom->SetLineColor(kRed);

      abs_mom->SetLineWidth(3);
      abs_mom->SetLineColor(kBlue);
      cex_mom->SetLineWidth(3);
      cex_mom->SetLineColor(kGreen);
      prod_mom->SetLineWidth(3);
      prod_mom->SetLineColor(kMagenta);
      inel_mom->SetLineWidth(3);
      inel_mom->SetLineColor(kOrange);

      TMultiGraph *mg = new TMultiGraph();
      mg->Add(lads_graph, "AP");
      mg->Add(totInel_mom, "AC");
      mg->Add(inel_mom, "AC");
      mg->Add(abs_mom, "AC");
      mg->Add(cex_mom, "AC");
      mg->Add(prod_mom, "AC");
      
      mg->GetXaxis()->SetRangeUser(0,2000);
      mg->GetXaxis()->SetTitle("Pion kinetic Energy [MeV]");
      mg->GetYaxis()->SetTitle("#sigma [mbarn]");
      
      TCanvas* c_geant= new TCanvas("geant","geant");
      gPad->SetGrid(1,1);
      mg->Draw("A");
      legend->Draw();
   c_geant->Update();
   //c_geant->Write();

 
   return 0;
}

