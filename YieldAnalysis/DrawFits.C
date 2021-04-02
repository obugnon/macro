/*
 *  DrawFits.C
 *
 *  Created by Ophelie Bugnon on 01/11/20.
 *
 */
#include <vector>
#include <iostream>
#include <algorithm>
#include <stdio.h>

#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TString.h"
#include "TMath.h"
#include "TLatex.h"

#include "FitYield.C"
// #include </Users/obugnon/Documents/ALICE/AnalyseJPsi/macro/SetRangeAndNameTest.C>


//___________________________________________________________________________________________________________
//___________________________________________________________________________________________________________
void DrawMultipleFits(Int_t minCent, Int_t maxCent)
{
    //Data
    TCanvas *cYield = new TCanvas(Form("%d-%d%%",minCent, maxCent),"J/#psi Raa");
    TH1D* graphYield = GetYieldGraph(minCent, maxCent, kFALSE, kFALSE);
    graphYield->SetTitle(Form("%d-%d%%", minCent, maxCent));
    graphYield->SetMarkerColor(kBlack);
    graphYield->SetLineColor(kBlack);
    graphYield->SetFillColor(kBlack);
    graphYield->SetMarkerStyle(8);
    cYield->SetLogy();
    graphYield->Draw();
    gStyle->SetOptStat(0000);

    //TF1 file
    TFile *outputFile = new TFile("FitFunctionsYield.root", "UPDATE");
    TDirectory *outputList ;
    std::vector<std::vector<Double_t>> *vect;
    std::vector<std::vector<Double_t>> parameters;
    Double_t chi2;

    TString sRange;
    TString sTest;

    //Tests
    int nTests;
    eFunction fFit=kPowLawFree;
    std::vector<Double_t> tMinFit;

    Int_t iCent; 
    if(minCent == 30 && maxCent == 50) iCent = 2;
    else if(minCent == 10 && maxCent == 30) iCent = 1; 
    else iCent = 0 ;

    if(minCent<50) 
    {
        nTests = 4;
        tMinFit={0, 0.3, 0.65, 1};
    }    
    else {
        nTests = 3;
        tMinFit={0.3, 0.65, 1};
    }



    
    Int_t lineColor[5]={kRed+1, kAzure-1, kGreen+2, kOrange-3};
    Int_t lineStyle[5]={1, 10, 9, 7};

    TLegend* legend = new TLegend(0.4,0.61,0.89,0.89);
    legend->SetBorderSize(0);
    legend->SetHeader("J/#psi Yield","C");

    TF1* currentFunction;
    for(int i=0; i<nTests; i++)
    {
       
        sRange = SetRangeFunction(-4,-2.5, minCent, maxCent); 
        cout << sRange.Data() << endl;
        sTest = SetNameTest(kInvariantYield, fFit, kFALSE, kFALSE, tMinFit[i], 15, 0);
        cout << sTest.Data() << endl;

        outputList = (TDirectory*) outputFile->Get(Form("%s_%s", sRange.Data(), sTest.Data()));
        outputList->GetObject("parameters", vect);
        parameters = *vect;
        chi2=parameters[0][parameters[0].size()-1];

        currentFunction = (TF1*)outputList->Get("function");
        currentFunction->SetLineColor(lineColor[i]);
        currentFunction->SetLineStyle(lineStyle[i]);
        currentFunction->Draw("SAME");

        legend->AddEntry(currentFunction, Form("Power Law, fitting range %.2f-15 GeV/c, #chi^{2}/NDF=%.2f", tMinFit[i], chi2), "l");
    
    }
    legend->Draw();
}
