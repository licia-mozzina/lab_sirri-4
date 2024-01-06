#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif

#include "RooDataHist.h"
#include "RooHistPdf.h"
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

int higgs_CMS() {
    //Declaring the observable
    RooRealVar inv_mass ("inv_mass", "4_leptons_invariant_mass_distr", 70., 180.);
    inv_mass.setBins(37);

    // Importing the data
    RooDataHist inv_mass_dh ("inv_mass_dh", "histo_inv_mass", inv_mass);
    RooDataHist ttbar_dh("ttbar_dh", "histo_ttbar", inv_mass);  
    RooDataHist ZZ_dh("ZZ_dh", "histo_ZZ", inv_mass);  
    RooDataHist DY_dh("DY_dh", "histo_DY", inv_mass);  

    ifstream f1("_cms_higgs_data.txt");
    Double_t v, w; // basterebbe una dichiarazione, poi cambiano sempre valore, anche per nome file
    while(!f1.eof()) {
        f1 >> v >> w;
        inv_mass.setVal(v);
        inv_mass_dh.set(inv_mass, w);
    }

    ifstream f2("_cms_higgs_TTbarto4l.txt");
    while(!f2.eof()) {
        f2 >> v >> w;
        inv_mass.setVal(v);
        ttbar_dh.set(inv_mass, w);
    }
  
    ifstream f3("_cms_higgs_ZZto4l.txt");
    while(!f3.eof()) {
        f3 >> v >> w;
        inv_mass.setVal(v);
        ZZ_dh.set(inv_mass, w);
    }

    ifstream f4("_cms_higgs_DYto4l.txt");
    while(!f4.eof()) {
        f4 >> v >> w;
        inv_mass.setVal(v);
        DY_dh.set(inv_mass, w);
    }
    
    //Building background distributions from data
    RooHistPdf hpdf_ttbar("hpdf_ttbar", "histpdf_ttbar", inv_mass, ttbar_dh);
    RooHistPdf hpdf_ZZ("hpdf_ZZ", "histpdf_ZZ", inv_mass, ZZ_dh);
    RooHistPdf hpdf_DY("hpdf_DY", "histpdf_DY", inv_mass, DY_dh);

    // Building the composite not extended model for background
    // Defining the fraction to be proportional to the relative weight of the histogram over the obs range
    // via histo integral
    RooRealVar f_ttbar ("f_ttbar", "fraction_bkg_ttbar", ttbar_dh.sum(kTRUE) / (180. - 70.), 0., 1.); 
    f_ttbar.setConstant();
    RooRealVar f_DY ("f_DY", "fraction_bkg_DY", DY_dh.sum(kTRUE) / (180. - 70.), 0., 1.);
    f_DY.setConstant();


    RooAddPdf bkg_model ("bkg_model", "background_model", RooArgList(hpdf_ttbar, hpdf_DY, hpdf_ZZ), RooArgList(f_ttbar, f_DY));
    RooPlot *higgs1 = inv_mass.frame();
    ttbar_dh.plotOn(higgs1, Name("Data"));
    bkg_model.plotOn(higgs1, Name("Model Fit")); 
    bkg_model.plotOn(higgs1, Components(hpdf_ttbar), LineColor(kViolet-1), Name("tt"));
    bkg_model.plotOn(higgs1, Components(hpdf_DY), LineColor(kViolet-1), Name("dy"));
    bkg_model.plotOn(higgs1, Components(hpdf_ZZ), LineColor(kMagenta), Name("zz"));

    TCanvas *c = new TCanvas("c", "HIGGS", 1600, 800);
    higgs1->Draw();

    // Building signal + background model
    RooRealVar m_higgs ("m_higgs", "Higgs_mass", 125., 110., 140);
    RooConstVar sigma_higgs ("sigma_higgs", "Higgs_mass_sigma", 3.); //fixed to 3 GeV
    RooRealVar f_s ("f_s", "signal_fraction", 0.5, 0., 1.);
    RooGaussian gaus_higgs ("gaus_higgs", "higgs_signal_distr", inv_mass, m_higgs, sigma_higgs);
    RooAddPdf model ("model", "Higgs_composite_model", RooArgList(gaus_higgs, bkg_model), f_s);

    // Fitting, saving fit values and plotting
    RooFitResult * result = model.fitTo(inv_mass_dh, Save());

    RooPlot *higgs = inv_mass.frame();
    inv_mass_dh.plotOn(higgs, Name("Data"));
    model.plotOn(higgs, Components(gaus_higgs), LineColor(kRed), Name("gaus"));
    model.plotOn(higgs, Components(bkg_model), LineColor(kBlack), Name("bkg"));
    model.plotOn(higgs, Name("Model Fit")); 

    //Best fit parameters
    model.paramOn(higgs, Label("fit result"), Format("NE",AutoPrecision(1)));

    TCanvas *c1 = new TCanvas("c1", "HIGGS", 1600, 800);
    higgs->Draw();

    TLegend *leg1 = new TLegend(0.70, 0.70, 0.85, 0.87);
    leg1->SetFillColor(kWhite);
    leg1->SetLineColor(kBlack);
    leg1->AddEntry("Data", "Data", "P");
    leg1->AddEntry("Model Fit", "Fit function", "LP");
   
    leg1->Draw("SAME");
    c1->Print("higgs.png");    

   // result->Print("v"); 

    //Writing fit results to txt file 
    ofstream myfile;
    myfile.open("minos_fit.txt");
    result->printMultiline(myfile, 1111, kTRUE);
    myfile.close();



    // Creating the ModelConfig for the 95% confidence interval
    RooWorkspace ws("ws");

    ModelConfig mc("ModelConfig", &ws);

    ws.import(model);
    //mc.SetParametersOfInterest(RooArgList(*ws.var("model"), *ws.var("inv_mass"), *ws.var("m_higgs")));
   // mc.SetNuisanceParameters(*ws.var("f_s"));
    ws.import(mc);

    // Direi che non servono con macro root apposta
    // RooAbsReal* ll = model.createNLL(inv_mass_dh);
    // RooAbsReal* pll = ll->createProfile(m_higgs);

    ws.writeToFile("higgs_CMS.root", true);

    return 0;
}