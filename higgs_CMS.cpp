#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif

#include "RooDataHist.h"
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
    //Declaring the observables
    RooRealVar inv_mass ("inv_mass", "4_leptons_invariant_mass_distr", 70., 180.);
    inv_mass.setBins(37);

    RooRealVar ttbar ("ttbar", "ttbar_decay", 70., 180.);
    ttbar.setBins(37);

    RooRealVar ZZ ("ZZ", "ZZ_decay", 70., 180.);
    ZZ.setBins(37);

    RooRealVar DY ("DY", "DY_decay", 70., 180.);
    DY.setBins(37);

    // Importing the data
    RooDataHist inv_mass_dh("inv_mass_dh", "histo_inv_mass", inv_mass);  
    RooDataHist ttbar_dh("ttbar_dh", "histo_ttbar", ttbar);  
    RooDataHist ZZ_dh("ZZ_dh", "histo_ZZ", ZZ);  
    RooDataHist DY_dh("DY_dh", "histo_DY", DY);  

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
        ttbar.setVal(v);
        ttbar_dh.set(ttbar, w);
    }
  
    ifstream f3("_cms_higgs_ZZto4l.txt");
    while(!f3.eof()) {
        f3 >> v >> w;
        ZZ.setVal(v);
        ZZ_dh.set(ZZ, w);
    }

    ifstream f4("_cms_higgs_DYto4l.txt");
    while(!f4.eof()) {
        f4 >> v >> w;
        DY.setVal(v);
        DY_dh.set(DY, w);
    }
    
    //Building background distributions from data
    RooHistPdf hpdf_ttbar("hpdf_ttbar", "histpdf_ttbar", ttbar, ttbar_dh, 0);
    RooHistPdf hpdf_ZZ("hpdf_ZZ", "histpdf_ZZ", ZZ, ZZ_dh, 0);
    RooHistPdf hpdf_DY("hpdf_DY", "histpdf_DY", DY, DY_dh, 0);

    // Building the composite not extended model for background
    // Defining the fraction to be proportional to the relative weight of the histogram over the obs range
    // via histo integral
    RooRealVar f_ttbar ("f_ttbar", "fraction_bkg_ttbar", ttbar_dh.sum(kTRUE) / (180. - 70.), 0., 1.); //vedo se mettere intervallo
    f_ttbar.setConstant();
    RooRealVar f_DY ("f_DY", "fraction_bkg_DY", DY_dh.sum(kTRUE) / (180. - 70.), 0., 1.);
    f_DY.setConstant();

    RooAddPdf bkg_model ("bkg_model", "background_model", RooArgList(hpdf_ttbar, hpdf_DY, hpdf_ZZ), RooArgList(f_ttbar, f_DY));

    // Building signal + background model
    RooRealVar m_higgs ("m_higgs", "Higgs_mass", 125., 110., 140);
    RooConstVar sigma_higgs ("sigma_higgs", "Higgs_mass_sigma", 3.); //fixed to 3 GeV
    RooRealVar f_s ("f_s", "signal_fraction", 0.5, 0., 1.);
    RooGaussian gaus_higgs ("gaus_higgs", "higgs_signal_distr", inv_mass, m_higgs, sigma_higgs);
    RooAddPdf model ("model", "Higgs_composite_model", RooArgList(gaus_higgs, bkg_model), f_s);

    // Fitting and plotting
    //model.fitTo(inv_mass_dh);

    RooPlot *higgs = inv_mass.frame();
    inv_mass_dh.plotOn(higgs, Name("Data"));
    model.plotOn(higgs, Name("Model Fit")); 
    model.plotOn(higgs, Components(gaus_higgs), LineColor(kViolet-1), Name("gaus"));
    model.plotOn(higgs, Components(bkg_model), LineColor(kMagenta), Name("bkg"));

    TCanvas *c1 = new TCanvas("c1", "HIGGS", 1600, 800);
    higgs->Draw();
    c1->Print("higgs.png");       



    return 0;
}