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


// Choose between single-channel and 4-independent channels mode
int OPERA(Int_t channel_mode) {

  //Int_t n_channel = channel_mode;

//   if ((channel_mode != 1) || (channel_mode != 4)) {
//     std::cout << "Entry is not valid, select 1 or 4\n";
//     return 0;
//   }

   if (channel_mode == 1) {
    RooConstVar n_signal("n_signal", "nominal_signal", 2.64);
    RooConstVar sigma_b("sigma_b", "background_error", 0.05);
    RooRealVar nobs("nobs", "observed_events", 0., 10.);
    RooRealVar b("b", "background_events", 0.25, 0., 0.5);
    RooRealVar sig_strength("sig_strength", "#mu", 0.5, 0., 1.);

    RooFormulaVar nexp("nexp", "@0 * @1 + @2",
                       RooArgSet(sig_strength, n_signal, b));

    RooPoisson model("model", "expected_events_distr", nobs, nexp);

    RooRealVar b0("b0", "best_bg_estimation", 0.25, 0., 1.);

    RooGaussian bg_constraint("bg_gconstraint", "bg_gaussian_constraint", b0, b, sigma_b);

    RooProdPdf model_constraint("model_constraint", "model with constraint", RooArgSet(model, bg_constraint));

    RooWorkspace w("w");

    w.import(model_constraint);

    ModelConfig mc("ModelConfig", &w);
    mc.SetPdf(*w.pdf("model"));
    mc.SetParametersOfInterest(*w.var("sig_strength"));
    mc.SetObservables(*w.var("nobs"));
    mc.SetSnapshot(*w.var("sig_strength"));   // needed for hypothesis test
    mc.SetGlobalObservables("b0"); // b0 should be treated as auxiliary obs in
                                   // frequenist stat and varied in pseudo exp
    w.var("b0")->setConstant(true);    // it has a range as glob obs, but it is const //potrebbe non servire perché è già const
    w.import(mc); // va importato prima di data, altrimenti non lo trova

    RooDataSet data("data", "", nobs);
    nobs.setVal(5);
    data.add(nobs);
    w.import(data);
    w.writeToFile("OPERA.root", true);
  }

  else if (channel_mode == 4) {
    return 0;
  }

  return 0;
}