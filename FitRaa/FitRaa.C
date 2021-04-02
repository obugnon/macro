/*
 *  FitRaa.C
 *
 *  Created by Ophelie Bugnon on 21/07/20.
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


 std::vector<std::vector<Double_t>> vectPt0 = {
     //0 pour free, mass Jpsi, mass Jpsi/2, meanPt, 2*meanPt
     {0, 3.097, 1.549, 2.10, 4.20},
     {0, 3.097, 1.549, 2.19, 4.38},
     {0, 3.097, 1.549, 2.34, 4.68}
    //meanPt pour 50-70 : 2.43 et pour 70-90 : 2.34
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

    Int_t iCent; 
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
void FitRaa(eFunction fitFunction, Int_t minCent, Int_t maxCent, Double_t minFit, Double_t maxFit, Double_t nPT0, Bool_t isSaved=kTRUE, Bool_t isErrStatOnly=kFALSE, Bool_t isErrSystOnly=kFALSE)
{
    TString sRange = SetRangeFunction(-4,-2.5, minCent, maxCent);  
    TString sTest = SetNameTest(kRaa, fitFunction, isErrStatOnly, isErrSystOnly, minFit, maxFit, nPT0);

    TCanvas *cRaa = new TCanvas(Form("%d-%d%%_%s", minCent, maxCent, sTest.Data()),"J/#psi Raa");
    TH1D* graphRaa = GetRAAGraph(minCent, maxCent, isErrStatOnly, isErrSystOnly);
    graphRaa->Draw();
    gStyle->SetOptFit(1111);
    graphRaa->SetTitle(Form("%d-%d%%_%s", minCent, maxCent, sTest.Data()));


    //Fit
    TF1* fFit; 
    TFitResultPtr fitStatus;
    Int_t nPar;

    if (fitFunction == kWoodSaxonFree)
    {
        nPar = 4;
        fFit = GetWoodSaxonLikeFunction(0.25, 3.097, 0, 15, kFALSE);
        fitStatus=graphRaa->Fit("WoodSaxonLikeFunction","EIS","", minFit, maxFit);
    }
    else if (fitFunction == kWoodSaxonFixed)
    {
        nPar = 4;
        fFit = GetWoodSaxonLikeFunction(0.2, nPT0, 0, 15, kTRUE);
        fitStatus=graphRaa->Fit("WoodSaxonLikeFunction","EIS","",  minFit, maxFit);
    }
    else if (fitFunction == kPol1)
    {
        nPar = 2;
        fFit = GetPol1Function(1, -0.01, 0, 15);
        fitStatus=graphRaa->Fit("Pol1Function","EIS","", minFit, maxFit);
    }
    else if (fitFunction == kPol3)
    {
        nPar = 4;
        fFit = GetPol3Function(1, 0, 15);
        fitStatus=graphRaa->Fit("Pol3Function","EIS","", minFit, maxFit);
    }
    else if (fitFunction == kConst)
    {
        nPar = 2;
        fFit = GetPol1Function(1, 0, 0, 15, kTRUE);
        fitStatus=graphRaa->Fit("Pol1Function","EIS","", minFit, maxFit);
    }

    //Quality of the fit 
    Double_t chi2 = fFit->GetChisquare()/fFit->GetNDF();
    Double_t covMatrixStatus = fitStatus->CovMatrixStatus();

    cout << "Fit test " << sTest.Data() << " : fits status = " << fitStatus << ", cov matrix status = " << covMatrixStatus << " and chi2/NDF = " << chi2 << endl;

    if(fitStatus!=0 || chi2 > 2.5 || covMatrixStatus!=3){
        printf("--------------------------------------\n Fit didn't converge !!!!! \n --------------------------------------\n");
        // return;
    }


    //Draw extrapolated function
    fFit->SetLineStyle(2);
    fFit->Draw("SAME");

    //Save function results 
    Double_t params[nPar];
    fFit->GetParameters(params);
    const Double_t *param_errors=fFit->GetParErrors();

    std::vector<std::vector<Double_t>> fitParameters;
    std::vector<Double_t> vParams;
    vParams.assign(params, params+nPar);
    vParams.push_back(chi2);
    std::vector<Double_t> vParam_errors;
    vParam_errors.assign(param_errors, param_errors+nPar);
    vParam_errors.push_back(0);
    fitParameters.push_back(vParams);
    fitParameters.push_back(vParam_errors);

    TMatrixDSym covarianceMatrix = fitStatus->GetCovarianceMatrix();

    if (isSaved)
    {
        //Save Histogram
        cRaa->SaveAs(".pdf");
        TFile *outputFile = new TFile("FitFunctionsRaa.root", "UPDATE");
        TDirectory *outputList = (TDirectory*) outputFile->Get(Form("%s_%s", sRange.Data(), sTest.Data()));
        if(!outputList) 
        {
            outputList = outputFile->mkdir(Form("%s_%s", sRange.Data(), sTest.Data()));
        }
        outputList->cd();
        fFit->Write("function");
        outputList->WriteObject(&fitParameters,"parameters");
        outputList->WriteObject(&covarianceMatrix,"covarianceMatrix");
        outputFile->Close();
    }
    
}
//___________________________________________________________________________________________________________
//___________________________________________________________________________________________________________
void MultipleFitRaa(Int_t minCent, Int_t maxCent, Bool_t isSaved)
{
    //Tests
    int nTests;
    std::vector<eFunction> fFit;
    std::vector<Double_t> tabPt0;
    std::vector<Double_t> minFit;

    Int_t iCent; 
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
            FitRaa(fFit[i], minCent, maxCent, minFit[j], 15, tabPt0[i], isSaved, kFALSE, kFALSE);
        }
    }

}
//___________________________________________________________________________________________________________
//___________________________________________________________________________________________________________
void DrawMultipleFits(Int_t minCent, Int_t maxCent)
{
    //Data
    TCanvas *cRaa = new TCanvas(Form("%d-%d%%",minCent, maxCent),"J/#psi Raa");
    TH1D* graphRaa = GetRAAGraph(minCent, maxCent, kFALSE, kFALSE);
    graphRaa->SetTitle(Form("%d-%d%%", minCent, maxCent));
    graphRaa->SetMarkerColor(kBlack);
    graphRaa->SetLineColor(kBlack);
    graphRaa->SetFillColor(kBlack);
    graphRaa->SetMarkerStyle(8);
    graphRaa->Draw();
    graphRaa->GetYaxis()->SetRangeUser(0.1, 1);
    gStyle->SetOptStat(0000);

    //TF1 file
    TFile *inputFile = new TFile("FitFunctionsRaa.root");
    TDirectory *inputList ;
    std::vector<std::vector<Double_t>> *vect;
    std::vector<std::vector<Double_t>> parameters;
    Double_t chi2;

    TString sRange;
    TString sTest;

    //Tests
    int nTests;
    std::vector<eFunction> fFit;
    std::vector<Double_t> tabPt0;

    Int_t iCent; 
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

    
    Int_t color[5]={kRed+1, kAzure-1, kGreen+2, kViolet-6, kOrange-3};
    Double_t tMinFit[2]={0.65, 1};
    // Double_t tMinFit[3]={0, 0.65, 1};
    Int_t lineStyle[3]={1, 10, 4};

    TLegend* legend = new TLegend(0.4,0.6,0.89,0.89);
    legend->SetBorderSize(0);
    legend->SetHeader("J/#psi R_{AA}","C");
    TString sLegend;

    TF1* currentFunction;
    for(int i=0; i<2; i++)
    {
        for(int j=0; j<nTests; j++)
        {
            sRange = SetRangeFunction(-4,-2.5, minCent, maxCent); 
            sTest = SetNameTest(kRaa, fFit[j], kFALSE, kFALSE, tMinFit[i], 15, tabPt0[j]);
            //if(j==2 || (j>3 && tMinFit[i]==1)) continue;

            inputList = (TDirectory*) inputFile->Get(Form("%s_%s", sRange.Data(), sTest.Data()));
            inputList->GetObject("parameters", vect);
            parameters = *vect;
            chi2=parameters[0][parameters[0].size()-1];
            // printf("chi2/NDF = %.3f\n", chi2);

            currentFunction = (TF1*)inputList->Get("function");
            currentFunction->SetLineColor(color[j]);
            currentFunction->SetLineStyle(lineStyle[i]);

            currentFunction->Draw("SAME");

            switch (fFit[j])
            {
                case kWoodSaxonFixed:
                sLegend = Form("WS with p_{T0} = %.3f", tabPt0[j]);
                break;
                case kWoodSaxonFree:
                sLegend = "WS free parameters";
                break;
                case kPol1:
                sLegend = "Linear fit";
                break;
                case kPol3:
                sLegend = "Third degree polynomial";
                break;
                case kConst:
                sLegend = "Constant fit";
                break;
            }
            legend->AddEntry(currentFunction, Form("%s, fitting range %.2f-15 GeV/c, #chi^{2}=%.2f", sLegend.Data(), tMinFit[i], chi2), "l");
        }
    }
    legend->Draw();

 
}