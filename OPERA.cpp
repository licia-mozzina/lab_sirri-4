#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif

#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "RooFitResult.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TStyle.h"
#include "TH2.h"

using namespace RooFit;
using namespace RooStats;

void OPERA( /* int nobs = 3,           // number of observed events
                     double b = 1,           // number of background events
                     double sigmab = 0.2 */)   // relative uncertainty in b
{

   RooConstVar n_signal ("n_signal", "nominal_signal", 2.64);
   RooConstVar bg_err ("bg_err", "background_error", 0.05);
   RooRealVar nobs ("nobs", "observed_events", 0., 10.);
   RooRealVar bg ("bg", "background_events", 0.25, 0., 0.5);
   RooRealVar sig_strength ("sig_strength", "#mu", 0.5, 0., 1.);

   RooFormulaVar nexp ("nexp", "@0 * @1 + @2", RooArgSet(sig_strength, n_signal, bg));

   RooPoisson model("model", "expected_events_distr", nobs, nexp);

   RooRealVar bg_est ("bg_est", "best_bg_estimation", 0.25, 0., 1.);

   RooGaussian bg_constraint ("bg_gconstraint", "bg_gaussian_constraint", bg_est, bg, bg_err);

   RooProdPdf model_constraint("model_constraint","model with constraint",RooArgSet(model,bg_constraint));

   RooDataSet data("data","", nobs);
   nobs.setVal(5);
   data.add(nobs);

   RooWorkspace w("w");

   w.import(model_constraint);
   w.import(data);
   w.writeToFile("OPERA.root", true);

   // TFile f ("OPERA.root", "recreate");
   // w.Write();
   // f.Close();
   
// // make Poisson model * Gaussian constraint
//    w.factory("sum:nexp(s[3,0,15],b[1,0,10])");
// // Poisson of (n | s+b)
//    w.factory("Poisson:pdf(nobs[0,50],nexp)");
//    w.factory("Gaussian:constraint(b0[0,10],b,sigmab[1])");
//    w.factory("PROD:model(pdf,constraint)");


//    w.var("b0")->setVal(b);
//    w.var("b0")->setConstant(true); // needed for being treated as global observables
//    w.var("sigmab")->setVal(sigmab*b);  
   

//    ModelConfig mc("ModelConfig",&w);
//    mc.SetPdf(*w.pdf("model"));
//    mc.SetParametersOfInterest(*w.var("s"));
//    mc.SetObservables(*w.var("nobs"));
//    mc.SetNuisanceParameters(*w.var("b"));

//    // these are needed for the hypothesis tests
//    mc.SetSnapshot(*w.var("s"));
//    mc.SetGlobalObservables(*w.var("b0"));

//    mc.Print();
//    // import model in the workspace 
//    w.import(mc);

//    // make data set with the namber of observed events
//    RooDataSet data("data","", *w.var("nobs"));
//    w.var("nobs")->setVal(3);
//    data.add(*w.var("nobs") );
//    // import data set in workspace and save it in a file
//    w.import(data);

//    w.Print();

//    TString fileName = "CountingModel.root"; 

//    // write workspace in the file (recreate file if already existing)
//    w.writeToFile(fileName, true);

}