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

#include "RooStats/ModelConfig.h"
#include "RooStats/FeldmanCousins.h"

#include <iostream>

using namespace RooFit;
using namespace RooStats;

// Choose between profile likelihood calculator, 
int higgs_CMS2(Int_t mode = 1) {
    //Declaring the observable
    Double_t low_lim = 70.;
    Double_t upp_lim = 180.;
    RooRealVar inv_mass ("inv_mass", "4_leptons_invariant_mass_distr", low_lim, upp_lim);
    Int_t n_bins = 37;
    inv_mass.setBins(n_bins);

    // Importing the data
    RooDataHist inv_mass_dh ("inv_mass_dh", "histo_inv_mass", inv_mass);
    RooDataHist ttbar_dh("ttbar_dh", "histo_ttbar", inv_mass);  
    RooDataHist ZZ_dh("ZZ_dh", "histo_ZZ", inv_mass);  
    RooDataHist DY_dh("DY_dh", "histo_DY", inv_mass);  

    ifstream f1("_cms_higgs_data.txt");
    Double_t v, w; 
    while(!f1.eof()) {
        f1 >> v >> w;
        if (v > low_lim && v < upp_lim) { // double-checking for underflow/overflow values 
                                            // anche perchè altrimenti per n_bin < 37 eof salta ad ultimo valore  
                                            // e lo assegna due volte a v prima di uscire dal loop
            inv_mass.setVal(v);
            inv_mass_dh.set(inv_mass, w);
        } else {
            continue;
        }
    }

    ifstream f2("_cms_higgs_TTbarto4l.txt");
    while(!f2.eof()) {
        f2 >> v >> w;
        if (v > low_lim && v < upp_lim) {
            inv_mass.setVal(v);
            ttbar_dh.set(inv_mass, w);
        } else {
            continue;
        }
    }
  
    ifstream f3("_cms_higgs_ZZto4l.txt");
    while(!f3.eof()) {
        f3 >> v >> w;
        if (v > low_lim && v < upp_lim) {
            inv_mass.setVal(v);
            ZZ_dh.set(inv_mass, w);
        } else {
            continue;
        }
    }

    ifstream f4("_cms_higgs_DYto4l.txt");
    while(!f4.eof()) {
        f4 >> v >> w;
        if (v > low_lim && v < upp_lim) {
            inv_mass.setVal(v);
            DY_dh.set(inv_mass, w);
        } else {
            continue;
        }
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

    // RooPlot *higgs1 = inv_mass.frame();
    // ttbar_dh.plotOn(higgs1, Name("Data"));
    // bkg_model.plotOn(higgs1, Name("Model Fit")); 
    // bkg_model.plotOn(higgs1, Components(hpdf_ttbar), LineColor(kViolet-1), Name("tt"));
    // bkg_model.plotOn(higgs1, Components(hpdf_DY), LineColor(kViolet-1), Name("dy"));
    // bkg_model.plotOn(higgs1, Components(hpdf_ZZ), LineColor(kMagenta), Name("zz"));

    // TCanvas *c = new TCanvas("c", "HIGGS", 1600, 800);
    // higgs1->Draw();

    // Building signal + background model
    RooRealVar m_higgs ("m_higgs", "Higgs_mass", 125., 110., 140.);
    RooConstVar sigma_higgs ("sigma_higgs", "Higgs_mass_sigma", 3.); //fixed to 3 GeV
    RooRealVar f_s ("f_s", "signal_fraction", 0.5, 0., 1.);
    RooFormulaVar f_b ("f_b", "1 -@0", f_s);
    RooGaussian gaus_higgs ("gaus_higgs", "higgs_signal_distr", inv_mass, m_higgs, sigma_higgs);
    RooAddPdf model ("model", "Higgs_composite_model", RooArgList(gaus_higgs, bkg_model), RooArgSet(f_s, f_b)); //SERVE IL COEFFICIENTE PER ENTRAMBI

    // Fitting, saving fit values and plotting
    RooFitResult * result = model.fitTo(inv_mass_dh, Save());

    RooPlot *higgs = inv_mass.frame();
    higgs->SetTitle("4 leptons invariant mass");
    higgs->GetXaxis()->SetTitle("m_{4l} [GeV/c^{2}]");
    higgs->GetYaxis()->SetTitle("Events / 3 GeV");
    inv_mass_dh.plotOn(higgs, Name("Data"));
    model.plotOn(higgs, Components(gaus_higgs), LineColor(kRed), Name("Signal"));
    model.plotOn(higgs, Components(bkg_model), LineColor(kBlack), Name("Bkg"));
    model.plotOn(higgs, Name("Model Fit")); 

    //Best fit parameters
    model.paramOn(higgs, Label("fit result"), Format("NE",AutoPrecision(1)), Layout(0.69, 0.75, 0.65));

    TCanvas *c1 = new TCanvas("c1", "HIGGS", 1600, 800);
    higgs->Draw();

    TLegend *leg1 = new TLegend(0.70, 0.70, 0.85, 0.87);
    leg1->SetFillColor(kWhite);
    leg1->SetLineColor(kBlack);
    leg1->AddEntry("Data", "Data", "P");
    leg1->AddEntry("Signal", "Signal", "L");
    leg1->AddEntry("Bkg", "Background", "L");
    leg1->AddEntry("Model Fit", "Signal + Background Fit", "LP");
   
    leg1->Draw("SAME");
    c1->Print("higgs2.png");  

   // result->Print("v"); 

    //Writing fit results to txt file 
    ofstream myfile;
    myfile.open("higgs_CMS_fit2.txt");
    result->printMultiline(myfile, 1111, kTRUE);
    myfile.close();



    // Creating the ModelConfig 
    RooWorkspace ws ("ws", "Workspace");

    ModelConfig mc ("mc", "ModelConfig", &ws);

    ws.import(model);
    ws.import(inv_mass_dh);
//    mc.SetParametersOfInterest(RooArgSet(*ws.pdf("model"), *ws.var("inv_mass"), *ws.var("m_higgs")));
//    mc.SetNuisanceParameters(*ws.var("f_s"));

    mc.SetPdf(*ws.pdf("model"));
    ws.defineSet("POI","m_higgs");
    ws.defineSet("nuisP","f_s"); 
    mc.SetParametersOfInterest(*ws.var("m_higgs"));
    mc.SetNuisanceParameters(*ws.var("f_s"));
    mc.SetObservables("inv_mass");

    ws.import(mc);

    ws.writeToFile("higgs_CMS2.root", true);


    //Profile Likelihood 95% confidence interval  
    ProfileLikelihoodCalculator plc(inv_mass_dh,mc);
    plc.SetConfidenceLevel(0.95); // 90% interval
    LikelihoodInterval* interval = plc.GetInterval();

    RooRealVar* firstPOI = (RooRealVar*) mc.GetParametersOfInterest()->first();
    double lowerLimit = interval->LowerLimit(*firstPOI);
    double upperLimit = interval->UpperLimit(*firstPOI);

    cout << "\n95 % interval on " <<firstPOI->GetName()<<" is : ["<<     lowerLimit << ", "<<     upperLimit <<"] "<<endl; 

    //FeldmanCousins
    FeldmanCousins fc (inv_mass_dh, mc);
    //fc.SetConfidenceLevel(0.90);        // quello sopra lo definisce per tutti 
    fc.SetPdf(model);
    fc.SetData(inv_mass_dh); 
    fc.SetParameters(m_higgs);
    fc.UseAdaptiveSampling(true);
    fc.FluctuateNumDataEntries(false);
    fc.SetNBins(100); // number of points to test per parameter
    fc.SetTestSize(.1);
    ConfInterval * fcint = fc.GetInterval();

    



    return 0;
}