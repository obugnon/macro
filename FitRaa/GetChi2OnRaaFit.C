/*
 *  GetChi2OnRaaFit.C
 *
 *  Created by Ophelie Bugnon on 28/01/21.
 *
 */
#include <vector>
#include <iostream>
#include <algorithm>
#include <stdio.h>

#include <TStyle.h>
#include "TFile.h"
#include "TH1.h"
#include "TGraphErrors.h"
#include "TLine.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TString.h"
#include "TMath.h"
#include "TLatex.h"
#include "TMatrixD.h"
#include "TFitResult.h"
#include "TLegend.h"

#include "FitFunctions.C"
#include "AliPWGFunc.h"

#include </Users/obugnon/Documents/ALICE/AnalyseJPsi/macro/ResultsRaa.C>
#include </Users/obugnon/Documents/ALICE/AnalyseJPsi/macro/ResultsPPCrossSection.C>
#include </Users/obugnon/Documents/ALICE/AnalyseJPsi/macro/ResultsAcceptanceEfficiency.C>
#include </Users/obugnon/Documents/ALICE/AnalyseJPsi/macro/SetRangeAndNameTest.C>


#ifndef JPSI_MASS
#define JPSI_MASS  3.096916
#define PSI2S_MASS 3.686109
#endif


// Double_t minFit = 0.;
// Double_t maxFit = 15.;

std::vector<std::vector<Double_t>> vectPt0 = {
     //0 pour free, mass Jpsi, mass Jpsi/2, meanPt, 2*meanPt
     {0, 3.097, 1.549, 2.09, 4.18},
     {0, 3.097, 1.549, 2.20, 4.40},
     {0, 3.097, 1.549, 2.35, 4.70}
    //meanPt pour 50-70 : 2.44 et pour 70-90 : 2.34
};

//___________________________________________________________________________________________________________
//___________________________________________________________________________________________________________
TH1D* GetRAAGraph(Int_t minCent, Int_t maxCent, Bool_t isErrStatOnly, Bool_t isErrSystOnly)
{
    //tableaux de la valeur moyenne du bin, de son ecart type et de sa valeur inférieure
    std::vector<std::vector<Double_t>> xRange;
    xRange.push_back({0.15, 0.475, 0.825, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 11, 13.5});
    xRange.push_back({0.15, 0.475, 0.825, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 9, 11, 13.5});
    xRange.push_back({0.15, 0.475, 0.825, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 9, 11});

    std::vector<std::vector<Double_t>> dxRange;
    dxRange.push_back({0.15, 0.175, 0.175, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1, 1.5});
    dxRange.push_back({0.15, 0.175, 0.175, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1, 1, 1.5});
    dxRange.push_back({0.15, 0.175, 0.175, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1, 1});

    std::vector<std::vector<Double_t>> xLowEdge;
    xLowEdge.push_back({0., 0.3, 0.65, 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 12., 15.});
    xLowEdge.push_back({0., 0.3, 0.65, 1., 2., 3., 4., 5., 6., 7., 8., 10., 12., 15.});
    xLowEdge.push_back({0., 0.3, 0.65, 1., 2., 3., 4., 5., 6., 7., 8., 10., 12.});

    Int_t iCent; // indice permettant d'acceder à la bonne ligne du tableau
    if(minCent == 70 && maxCent == 90) iCent = 2;
    else if(minCent == 50 && maxCent == 70) iCent = 1; 
    else iCent = 0 ; 

    Int_t nBins = xRange[iCent].size();

    Double_t* tLowEdge;
    tLowEdge = &xLowEdge[iCent][0];

    TH1D* graph = new TH1D("J/psi Raa", "", nBins, tLowEdge);
    graph->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    graph->GetYaxis()->SetTitle("R_{AA}^{J/#psi}");
    graph->SetMarkerColor(kBlue);
    graph->SetLineColor(kBlue);
    graph->SetFillColor(kBlue);
    graph->SetFillStyle(0);
    

    //Remplissage de l'histogramme
    std::vector<Double_t> vectRaa;
    std::vector<Double_t> vectCRpp;
    std::vector<Double_t> vectAccEff;
    Double_t MCinput_correlations;

    Double_t nRaa; 
    Double_t errRaa;
    Double_t errSystRaa;

    for(int i=0; i<nBins; i++)
    {
        vectRaa.clear();
        vectRaa = RaaResultsVSpT[SetRangeValue(kRaa, -4, -2.5, xRange[iCent][i]-dxRange[iCent][i], xRange[iCent][i]+dxRange[iCent][i], minCent, maxCent)];

        vectCRpp.clear();
        vectCRpp = ppCrossSection[SetRangeValue(kCRpp, -4, -2.5, xRange[iCent][i]-dxRange[iCent][i], xRange[iCent][i]+dxRange[iCent][i])];
        
        vectAccEff.clear();
        vectAccEff = AccEffHadro[SetRangeValue(kAccEff, -4, -2.5, xRange[iCent][i]-dxRange[iCent][i], xRange[iCent][i]+dxRange[iCent][i], minCent, maxCent)];

        nRaa=vectRaa[kRaaValue];

        // if(isErrStatOnly && !isErrSystOnly) errRaa=vectRaa[kRaaStat];
        // else if(isErrSystOnly && !isErrStatOnly) errRaa=vectRaa[kRaaSyst];
        // else errRaa=TMath::Sqrt(TMath::Power((vectRaa[kRaaStat]/nRaa),2) + TMath::Power((vectRaa[kRaaSyst]/nRaa),2))*nRaa;

        if(isErrStatOnly && !isErrSystOnly) errRaa = TMath::Sqrt(TMath::Power((vectRaa[kRaaStat]/nRaa),2) - TMath::Power((vectAccEff[kAccEffStatError]/vectAccEff[kAccEffValue]), 2)  -  TMath::Power((vectCRpp[kCRppStatError]/vectCRpp[kCRppValue]), 2))*nRaa ;

        else if(isErrSystOnly && !isErrStatOnly) errRaa=TMath::Sqrt(TMath::Power((vectRaa[kRaaSyst]/nRaa),2) - TMath::Power((vectCRpp[kCRppSystError]/vectCRpp[kCRppValue]), 2)) * nRaa;

        else errRaa=TMath::Sqrt(TMath::Power((vectRaa[kRaaStat]/nRaa),2) - TMath::Power((vectAccEff[kAccEffStatError]/vectAccEff[kAccEffValue]), 2)  -  TMath::Power((vectCRpp[kCRppStatError]/vectCRpp[kCRppValue]), 2) + TMath::Power((vectRaa[kRaaSyst]/nRaa),2) - TMath::Power((vectCRpp[kCRppSystError]/vectCRpp[kCRppValue]), 2))*nRaa;

        // printf("Uncertainty considered in %.2f-%.2f is %.4f%\n", xRange[iCent][i]-dxRange[iCent][i], xRange[iCent][i]+dxRange[iCent][i], errRaa);
        
        graph->SetBinContent(i+1, nRaa);
        graph->SetBinError(i+1, errRaa);
    }

    return graph;
}
//___________________________________________________________________________________________________________
//___________________________________________________________________________________________________________
void GetChi2OnRaaFit(eFunction fitFunction, Int_t minCent, Int_t maxCent, Double_t minFit, Double_t maxFit, Double_t nPT0, Bool_t isErrStatOnly=kFALSE, Bool_t isErrSystOnly=kFALSE)
{
    TString sRange = SetRangeFunction(-4,-2.5, minCent, maxCent);  
    TString sTest = SetNameTest(kRaa, fitFunction, isErrStatOnly, isErrSystOnly, minFit, maxFit, nPT0);

    TH1D* graphRaa = GetRAAGraph(minCent, maxCent, isErrStatOnly, isErrSystOnly);
    
    //TF1 file
    TFile *inputFile = new TFile("FitFunctionsRaa.root");
    TDirectory *inputList;
    TF1* currentFunction;

    inputList = (TDirectory*) inputFile->Get(Form("%s_%s", sRange.Data(), sTest.Data()));
    printf("Test %s_%s\n", sRange.Data(), sTest.Data());
    currentFunction = (TF1*)inputList->Get("function");

    // Compute chi2
    Double_t chi2=0;
    // Double_t val_pT[7]={0.15, 0.475, 0.825, 1.5, 2.5, 3.5, 4.5};    
    // Double_t val_pT[6]={0.475, 0.825, 1.5, 2.5, 3.5, 4.5};
    // Double_t val_pT[11]={1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 11, 13.5};
    Double_t val_pT[3]={0.15, 0.475, 0.825};    


    Double_t exp_value, exp_err, fit_value;

    for(int i=0; i<3; i++)
    {
        exp_value = graphRaa->GetBinContent(graphRaa->GetXaxis()->FindBin(val_pT[i]));
        exp_err = graphRaa->GetBinError(graphRaa->GetXaxis()->FindBin(val_pT[i]));
        fit_value = currentFunction->Eval(val_pT[i]);
        // printf("exp_value=%.2f and fit_value=%.2f ---> %.5f\n", exp_value, fit_value, TMath::Power(exp_value-fit_value, 2)/TMath::Power(exp_err,2));
        printf("exp_value=%.2f and fit_value=%.2f ---> %.5f\n", exp_value, fit_value, TMath::Abs(exp_value-fit_value)/exp_err,2);
        chi2 += TMath::Power(exp_value-fit_value, 2)/TMath::Power(exp_err,2);
    }
    Int_t NDF;
    if(fitFunction==kWoodSaxonFree) NDF=3;
    else NDF=3;
    printf("Test %s : Chi2=%.2f ---> Chi2/NDF = %.2f\n\n", sTest.Data() ,chi2, chi2/NDF);
    // printf("Chi2=%.2f and NDF=%i given chi2/NDF=%.2f\n", chi2, NDF, chi2/NDF);
}

//___________________________________________________________________________________________________________
//___________________________________________________________________________________________________________
void MultipleTests(Int_t minCent, Int_t maxCent)
{
    //Tests
    int nTests;
    std::vector<eFunction> fFit;
    std::vector<Double_t> tabPt0;
    std::vector<Double_t> minFit;

    Int_t iCent; // indice permettant d'acceder à la bonne ligne du tableau
    if(minCent == 30 && maxCent == 50) iCent = 2;
    else if(minCent == 10 && maxCent == 30) iCent = 1; 
    else iCent = 0 ;

    if(minCent<50) 
    {
        nTests = 5;
        fFit={kWoodSaxonFree, kWoodSaxonFixed, kWoodSaxonFixed, kWoodSaxonFixed, kWoodSaxonFixed};
        tabPt0=vectPt0[iCent];
    }    
    else {
        nTests = 2;
        fFit={kPol1, kConst};
        tabPt0={0, 0};
    }
    minFit={0.65, 1};
    // minFit={0};

    for(int i=0; i<fFit.size(); i++)
    {
        for(int j=0; j<minFit.size(); j++)
        {
            GetChi2OnRaaFit(fFit[i], minCent, maxCent, minFit[j], 15, tabPt0[i], kFALSE, kFALSE);
        }
    }

}