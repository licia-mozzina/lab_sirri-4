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

int prova() {
    //Declaring the observable
    RooRealVar inv_mass ("inv_mass", "4_leptons_invariant_mass_distr", 70., 180.);
    inv_mass.setBins(37);
    
    // Importing the data
    RooDataHist inv_mass_dh ("inv_mass_dh", "histo_inv_mass", inv_mass);
    RooDataHist ttbar_dh("ttbar_dh", "histo_ttbar", inv_mass);  //STESSA VARIABILE PER TUTTI
    RooDataHist ZZ_dh("ZZ_dh", "histo_ZZ", inv_mass);  
    RooDataHist DY_dh("DY_dh", "histo_DY", inv_mass); 

    //ttbar_dh.SumW2Error(kTRUE);

    ifstream f1("_cms_higgs_data.txt");
    Double_t v, w;
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

    RooPlot *higgs = inv_mass.frame();
    bkg_model.plotOn(higgs);

    higgs->Draw();
    
    return 0;
}