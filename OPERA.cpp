#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif

#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooGaussian.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH2.h"
#include "TStyle.h"

#include <iostream>

using namespace RooFit;
using namespace RooStats;


// Choose between single-channel and 4-independent channels mode, default single channel
int OPERA(int channel_mode = 1) {

   if (channel_mode == 1) {
    RooConstVar n_signal("n_signal", "nominal_signal", 2.64);
    RooConstVar sigma_b("sigma_b", "background_error", 0.05);
    RooRealVar nobs("nobs", "observed_events", 0., 10.);
    RooRealVar b("b", "background_events", 0.25, 0., 0.5);
    RooRealVar sig_strength("sig_strength", "#mu", 0.5, 0., 1.);

    RooFormulaVar nexp("nexp", "@0 * @1 + @2",
                       RooArgSet(sig_strength, n_signal, b));

    RooPoisson model("model", "expected_events_distr", nobs, nexp);

    RooRealVar b0("b0", "best_bg_estimate", 0.25, 0., 1.);

    RooGaussian bg_constraint("bg_gconstraint", "bg_gaussian_constraint", b0, b, sigma_b);

    RooProdPdf model_constraint("model_constraint", "model with constraint", RooArgSet(model, bg_constraint));

    RooWorkspace w("w");

    w.import(model_constraint);

    ModelConfig mc("ModelConfig", &w);
    mc.SetPdf(*w.pdf("model_constraint"));
    mc.SetParametersOfInterest(*w.var("sig_strength"));
    mc.SetObservables(*w.var("nobs"));
    mc.SetSnapshot(*w.var("sig_strength")); 
    mc.SetGlobalObservables("b0");
    w.var("b0")->setConstant(true);   
    w.import(mc);

    RooDataSet data("data", "", nobs);
    nobs.setVal(5);
    data.add(nobs);
    w.import(data);
    w.writeToFile("OPERA_1.root", true);
  }

  else if (channel_mode == 4) {

   RooRealVar sig_strength("sig_strength", "#mu", 0.5, 0., 1.); // one for all

   // first channel
    RooConstVar n_signal_1("n_signal_1", "nominal_signal", 0.52);
    RooConstVar sigma_b_1("sigma_b_1", "background_error", 0.01);
    RooRealVar nobs_1("nobs_1", "observed_events", 0., 10.);
    RooRealVar b_1("b_1", "background_events", 0.04, 0., 0.08);

    RooFormulaVar nexp_1("nexp_1", "@0 * @1 + @2", RooArgSet(sig_strength, n_signal_1, b_1));

    RooPoisson model_1("model_1", "expected_events_distr", nobs_1, nexp_1);

    RooRealVar b0_1("b0_1", "best_bg_estimate", 0.04, 0., 0.08);

    RooGaussian bg_constraint_1("bg_gconstraint_1", "bg_gaussian_constraint", b0_1, b_1, sigma_b_1);


   // second channel
    RooConstVar n_signal_2("n_signal_2", "nominal_signal", 0.73);
    RooConstVar sigma_b_2("sigma_b_2", "background_error", 0.03);
    RooRealVar nobs_2("nobs_2", "observed_events", 0., 10.);
    RooRealVar b_2("b_2", "background_events", 0.17, 0., 0.34);

    RooFormulaVar nexp_2("nexp_2", "@0 * @1 + @2", RooArgSet(sig_strength, n_signal_2, b_2));

    RooPoisson model_2("model_2", "expected_events_distr", nobs_2, nexp_2);

    RooRealVar b0_2("b0_2", "best_bg_estimate", 0.17, 0., 0.34);

    RooGaussian bg_constraint_2("bg_gconstraint_2", "bg_gaussian_constraint", b0_2, b_2, sigma_b_2);


   // third channel
    RooConstVar n_signal_3("n_signal_3", "nominal_signal", 0.61);
    RooConstVar sigma_b_3("sigma_b_3", "background_error", 0.001);
    RooRealVar nobs_3("nobs_3", "observed_events", 0., 10.);
    RooRealVar b_3("b_3", "background_events", 0.004, 0., 0.008);

    RooFormulaVar nexp_3("nexp_3", "@0 * @1 + @2", RooArgSet(sig_strength, n_signal_3, b_3));

    RooPoisson model_3("model_3", "expected_events_distr", nobs_3, nexp_3);

    RooRealVar b0_3("b0_3", "best_bg_estimate", 0.004, 0., 0.008);

    RooGaussian bg_constraint_3("bg_gconstraint_3", "bg_gaussian_constraint", b0_3, b_3, sigma_b_3);


   // fourth channel
    RooConstVar n_signal_4("n_signal_4", "nominal_signal", 0.78);
    RooConstVar sigma_b_4("sigma_b_4", "background_error", 0.01);
    RooRealVar nobs_4("nobs_4", "observed_events", 0., 10.);
    RooRealVar b_4("b_4", "background_events", 0.03, 0., 0.06);

    RooFormulaVar nexp_4("nexp_4", "@0 * @1 + @2", RooArgSet(sig_strength, n_signal_4, b_4));

    RooPoisson model_4("model_4", "expected_events_distr", nobs_4, nexp_4);

    RooRealVar b0_4("b0_4", "best_bg_estimate", 0.03, 0., 0.06);

    RooGaussian bg_constraint_4("bg_gconstraint_4", "bg_gaussian_constraint", b0_4, b_4, sigma_b_4);


    RooProdPdf model_constraint("model_constraint", "model with constraint", RooArgSet(model_1, model_2, model_3, model_4,
                                                             bg_constraint_1, bg_constraint_2, bg_constraint_3, bg_constraint_4));


    RooWorkspace w("w");

    w.import(model_constraint);

    ModelConfig mc("ModelConfig", &w);
    mc.SetPdf(*w.pdf("model_constraint"));
    mc.SetParametersOfInterest(*w.var("sig_strength"));
    mc.SetObservables(RooArgSet(*w.var("nobs_1"), *w.var("nobs_2"), *w.var("nobs_3"), *w.var("nobs_4")));
    mc.SetSnapshot(*w.var("sig_strength"));   // needed for hypothesis test
    mc.SetGlobalObservables(RooArgSet(b0_1, b0_2, b0_3, b0_4)); // b0 should be treated as auxiliary obs in frequenist stat and varied in pseudo exp
    w.var("b0_1")->setConstant(true);    // it has a range as glob obs, but it is const //potrebbe non servire perché è già const
    w.var("b0_2")->setConstant(true);   
    w.var("b0_3")->setConstant(true);    
    w.var("b0_4")->setConstant(true);    
    w.import(mc); 

    RooDataSet data("data", "", RooArgSet(nobs_1, nobs_2, nobs_3, nobs_4));
    nobs_1.setVal(3);
    nobs_2.setVal(1);
    nobs_3.setVal(1);
    nobs_4.setVal(0);
    data.add(RooArgSet(nobs_1, nobs_2, nobs_3, nobs_4));
    w.import(data);
    w.writeToFile("OPERA_4.root", true);
  }

   else {
    std::cout << "Entry is not valid, select 1 or 4\n";
  }

  return 0;
}