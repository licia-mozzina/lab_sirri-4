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

int higgs_CMS(Int_t mode = 1) {
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


    // Building signal + background model
    RooRealVar m_higgs ("m_higgs", "Higgs_mass", 125., 110., 140.);
    RooConstVar sigma_higgs ("sigma_higgs", "Higgs_mass_sigma", 3.); //fixed to 3 GeV
    RooRealVar f_s ("f_s", "signal_fraction", 0.5, 0., 1.);
    RooFormulaVar f_b ("f_b", "1 -@0", f_s);
    RooGaussian gaus_higgs ("gaus_higgs", "higgs_signal_distr", inv_mass, m_higgs, sigma_higgs);
    RooAddPdf model ("model", "Higgs_composite_model", RooArgList(gaus_higgs, bkg_model), RooArgSet(f_s, f_b)); 

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

    c1->Print("higgs_CMS.png"); 

    //Writing fit results to txt file 
    ofstream myfile;
    myfile.open("higgs_CMS_fit.txt");
    result->printMultiline(myfile, 1111, kTRUE);
    myfile.close();



    // Creating the the workspace and the ModelConfig 
    RooWorkspace ws ("ws", "Workspace");

    ModelConfig mc ("mc", "ModelConfig", &ws);

    ws.import(model);
    ws.import(inv_mass_dh);

    // Defining pdf, parameters of interest and nuisance parameters
    mc.SetPdf(*ws.pdf("model"));
    ws.defineSet("POI","m_higgs");
    ws.defineSet("nuisP","f_s"); 
    // ws.defineSet("POI","f_s");           // for computing p-values and significance of f_s as signal strength, null hypothesis f_s = 0
    // ws.defineSet("nuisP","m_higgs"); 
    mc.SetParametersOfInterest(*ws.set("POI"));
    mc.SetNuisanceParameters(*ws.set("nuisP"));
    mc.SetObservables("inv_mass");

    ws.import(mc);

    ws.writeToFile("higgs_CMS.root", true);  


    //***************************************************
    // Computing significance (p-value) as function of the signal mass
    //***************************************************

    // ModelConfig corresponds to signal + background model
    mc.SetName("S+B Model");      
    RooRealVar* poi = (RooRealVar*) mc.GetParametersOfInterest()->first(); // m_higgs
    
    // Setting POI snapshot in S+B model for expected significance
    poi->setVal(200);  
    mc.SetSnapshot(*poi);

    // Creating the B model from the S+B model
    ModelConfig * bModel = (ModelConfig*) mc.Clone();
    bModel->SetName("B Model");      
    poi->setVal(0);
    bModel->SetSnapshot(*poi);

    // Storing for mass values
    vector<Double_t> masses;
    vector<Double_t> p0values;
    vector<Double_t> p0valuesExpected;

    // Setting lower/upper bound for mass values
    Double_t m_higgs_min = 112;
    Double_t m_higgs_max = 158;

    // Looping on mass values 
    Int_t npoints = 30;

    for( Double_t mass=m_higgs_min; mass<=m_higgs_max; mass += (m_higgs_max-m_higgs_min)/Double_t(npoints) )   {
        cout << endl << endl << "Running for mass: " << mass << endl << endl;      
        ws.var("m_higgs")->setVal(mass);

        // Asymptotic calculator for the likelihood
        AsymptoticCalculator *  ac = new AsymptoticCalculator(inv_mass_dh, mc, *bModel);
        ac->SetOneSidedDiscovery(true);  // for one-side discovery test                                      
        AsymptoticCalculator::SetPrintLevel(-1);


        HypoTestResult* asymCalcResult = ac->GetHypoTest();
 
        asymCalcResult->Print();
     
        masses.push_back( mass );

        // P-values for null hypothesis
        p0values.push_back( asymCalcResult->NullPValue() );

        // Retrieving expected p-values
        Double_t expectedP0 = AsymptoticCalculator::GetExpectedPValues(  asymCalcResult->NullPValue(),  asymCalcResult->AlternatePValue(), 0, false);
        p0valuesExpected.push_back( expectedP0 );
        std::cout << "expected p0 = " << expectedP0 << std::endl;   
     }   
     
    //Plotting p-values
    TCanvas *c2 = new TCanvas("c2", "HIGGS", 1600, 800);

    TGraph * graph1  = new TGraph(masses.size(),&masses[0],&p0values[0]);   
    TGraph * graph2  = new TGraph(masses.size(),&masses[0],&p0valuesExpected[0]);   
    graph1->SetMarkerStyle(20);
    graph1->Draw("APC");
    graph2->SetLineStyle(2);
    graph2->Draw("SAME C");
    graph1->GetXaxis()->SetTitle("Mass [GeV/c^{2}]");
    graph1->GetYaxis()->SetTitle("p0 value");
    graph1->GetYaxis()->SetLimits(-1, 0.2);
    graph1->SetTitle("Significance vs Mass");
    graph1->SetMinimum(graph2->GetMinimum());
    graph1->SetLineColor(kBlue);
    graph2->SetLineColor(kRed);
    graph2->SetLineWidth(3);
    gPad->SetLogy(true);

    TLegend *leg2 = new TLegend(0.70, 0.67, 0.85, 0.82);
    leg2->SetFillColor(kWhite);
    leg2->SetLineColor(kBlack);
    leg2->AddEntry(graph1, "Observed p-values", "P");
    leg2->AddEntry(graph2, "Expected p-values", "L");
   
    leg2->Draw("SAME"); 

    c2->Print("p0_values.png");

    return 0;
}