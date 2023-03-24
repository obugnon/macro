/*
 *  FitAcceptanceEfficiency.C
 *
 *  Created by Ophelie Bugnon on 21/07/20.
 *
 */
#include <vector>
#include <iostream>
#include <algorithm>
#include <stdio.h>

#include <TStyle.h>
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
#include "TFile.h"

#include "FitFunctions.C"
#include "AliPWGFunc.h"
#include </Users/obugnon/Documents/ALICE/AnalyseJPsi/macro/ResultFiles/ResultsAcceptanceEfficiency.C>
#include </Users/obugnon/Documents/ALICE/AnalyseJPsi/macro/SetRangeAndNameTest.C>


#ifndef JPSI_MASS
#define JPSI_MASS   3.096916
#define PSI2S_MASS 3.686109
#endif

Double_t minFit = 0.;
Double_t maxFit = 15.;


TH1D* GetAccEffGraph(Double_t minCent, Double_t maxCent, Bool_t isErrStatOnly, Bool_t isErrSystOnly, Double_t minY, Double_t maxY)
{
    // ranges pT
    // tableaux de la valeur moyenne du bin, de son ecart type et de sa valeur inférieure
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

    Int_t centClase;
    if (minCent==0) centClase=0;
    else if (minCent==10) centClase=1;
    else if (minCent==30) centClase=2;
    else if (minCent==50) centClase=3;
    else if (minCent==70) centClase=4;

    std::vector<Double_t> vectAccEff;
    Double_t nAccEff; 
    Double_t errAccEff;
    Double_t errSystAccEff;

    Double_t EffresponseMTR;
    Double_t MCinput_correlations;

    Double_t* tLowEdge;
    tLowEdge = &xLowEdge[iCent][0];

    TH1D* graph = new TH1D("J/psi Acceptance Efficiency", "", nBins, tLowEdge);
    graph->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    graph->GetYaxis()->SetTitle("A#times#epsilon (J/#psi)");
    graph->SetMarkerColor(kBlue);
    graph->SetLineColor(kBlue);
    graph->SetFillColor(kBlue);
    graph->SetFillStyle(0);

    for(int i=0; i<nBins; i++)
    {
        vectAccEff.clear();
        vectAccEff = AccEffHadro[SetRangeValue(kAccEff, minY, maxY, xRange[iCent][i]-dxRange[iCent][i], xRange[iCent][i]+dxRange[iCent][i], minCent, maxCent)];

        EffresponseMTR = SystEffMTR[SetNameSystAccEff(minY, maxY, xRange[iCent][i]-dxRange[iCent][i], xRange[iCent][i]+dxRange[iCent][i])][0];
        MCinput_correlations = SystMCinput_correlations[SetNameSystAccEff(minY, maxY, xRange[iCent][i]-dxRange[iCent][i], xRange[iCent][i]+dxRange[iCent][i])][0];

        nAccEff=vectAccEff[kAccEffValue];
        errSystAccEff = TMath::Sqrt(TMath::Power(MC_inputVSpt_stat/100, 2) + TMath::Power(MCinput_correlations/100, 2) + TMath::Power((EffMCH/100), 2) + TMath::Power((EffresponseMTR/100), 2) + TMath::Power((EffMTR/100), 2) + TMath::Power(0.01, 2) )*nAccEff;

        if(isErrStatOnly && !isErrSystOnly) errAccEff=vectAccEff[kAccEffStatError];

        else if(isErrSystOnly && !isErrStatOnly) errAccEff = errSystAccEff;

        else errAccEff=TMath::Sqrt(TMath::Power(vectAccEff[kAccEffStatError]/vectAccEff[kAccEffValue],2) + TMath::Power(errSystAccEff/nAccEff, 2))*vectAccEff[kAccEffValue];

        graph->SetBinContent(i+1, nAccEff);
        graph->SetBinError(i+1, errAccEff);
    }


    return graph;
}

void FitAcceptanceEfficiency(eFunction fitFunction, Int_t minCent, Int_t maxCent, Double_t minFit, Double_t maxFit, Bool_t isSaved=kFALSE, Bool_t isErrStatOnly=kTRUE, Bool_t isErrSystOnly=kFALSE, Double_t minY=-4, Double_t maxY=-2.5)
{
    TString sTest = SetNameTest(kAccEff, fitFunction, isErrStatOnly, isErrSystOnly, minFit, maxFit);
    TString sRangeFunction = SetRangeFunction(minY, maxY, minCent, maxCent);
    TCanvas *cAccEff = new TCanvas(Form("%s-%s",sRangeFunction.Data(), sTest.Data()),"J/#psi Acceptance Efficiency");
    TH1D* graphAccEff = GetAccEffGraph(minCent, maxCent, isErrStatOnly, isErrSystOnly, minY, maxY);
    graphAccEff->Draw();
    graphAccEff->GetXaxis()->SetRangeUser(minFit, maxFit+2);
    gStyle->SetOptFit(1111);
    graphAccEff->SetTitle(Form("%s-%s", sRangeFunction.Data(), sTest.Data()));

    //Fit
    TF1* fFit; 
    TFitResultPtr fitStatus;
    Int_t nPar;
    Double_t chi2;
    Double_t covMatrixStatus;
    Int_t iter=0;
    // ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(10000); 
    // ROOT::Math::MinimizerOptions::SetDefaultTolerance(0.5); 
    do{
        if (fitFunction == kRatioLevy)
        {
            if(iter==0)
            {
                nPar = 7;
                fFit = GetRatioLevyFunction(0.12, minFit, maxFit);
            }
           fitStatus=graphAccEff->Fit("RatioLevyFunction","RES");
        }
        else if (fitFunction == kPol3)
        {
            if(iter==0)
            {
                nPar = 4;
                fFit = GetPol3Function(0.15, minFit, maxFit);
            }    
            fitStatus=graphAccEff->Fit("Pol3Function","ERS");
        }
        else if (fitFunction == kRatioPowerLaw)
        {
            if(iter==0)
            {
                nPar = 7;
                fFit = GetRatioPowerLaw(0.15, minFit, maxFit);
            }    
            fitStatus=graphAccEff->Fit("RatioPowerLaw","RES");
        }
        
        //Quality of the fit 
        chi2 = fFit->GetChisquare()/fFit->GetNDF();
        covMatrixStatus = fitStatus->CovMatrixStatus();

        iter++;
        if (iter > 10)
        {
            cout << "______________________________" << endl;
            cout << " The fit " << sRangeFunction << "-" << sTest << "has not converged "    << endl;    
            cout << "______________________________" << endl;
            break;
        }
    }
    // while (fitStatus !=0 ||  chi2 > 3 || covMatrixStatus != 3);
    while (fitStatus !=0 || covMatrixStatus != 3);
    cout << "Fit test " << sTest.Data() << " : fits status = " << fitStatus << ", cov matrix status = " << covMatrixStatus << " and chi2/NDF = " << chi2 << endl;
    

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

    if (isSaved)// 
    {
        //Save Histogram
        cAccEff->SaveAs(".pdf");
        TFile *outputFile = new TFile("$LOWPT/macro/ResultFiles/FitFunctionsAcceptanceEfficiency.root", "UPDATE");
        TDirectory *outputList = (TDirectory*) outputFile->Get(Form("%s_%s", sRangeFunction.Data(), sTest.Data()));
        if(!outputList) 
        {
            outputList = outputFile->mkdir(Form("%s_%s", sRangeFunction.Data(), sTest.Data()));
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
void MultipleFitAcceptanceEfficiency(Bool_t isSaved)
{
    //Tests
    int nTests=2;
    std::vector<eFunction> fFit={kRatioLevy, kPol3};
    Int_t centRange[6]={0, 10, 30, 50, 70, 90};

    for(int i=0; i<5; i++)
    {
        for(int j=0; j<2; j++)
        {
            FitAcceptanceEfficiency(fFit[j], centRange[i], centRange[i+1], 0., 10, isSaved);
        }
    }

}
//___________________________________________________________________________________________________________
//___________________________________________________________________________________________________________
void DrawMultipleFits(Int_t minCent, Int_t maxCent)
{
    //Data
    TCanvas *cAccEff = new TCanvas("","J/#psi Acceptance Efficiency");
    cAccEff->SetRightMargin(0.03);
    cAccEff->SetTopMargin(0.03);
    TH1D* graphAccEff = GetAccEffGraph(minCent, maxCent, kTRUE, kFALSE, -4, -2.5);
    graphAccEff->SetMarkerColor(kBlack);
    graphAccEff->SetLineColor(kBlack);
    graphAccEff->SetFillColor(kBlack);
    graphAccEff->SetMarkerStyle(8);
    graphAccEff->Draw();
    graphAccEff->GetYaxis()->SetRangeUser(0.09, 0.5);
    graphAccEff->GetXaxis()->SetRangeUser(0, 10);
    gStyle->SetOptStat(0000);

    TLatex * text = new TLatex (0.5,0.45,"Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV");
    text->Draw("SAME");
    text->SetTextSizePixels(20);

    TLatex * text1 = new TLatex (0.5,0.405,"J/#psi #rightarrow #mu^{+}#mu^{-}, 2.5 < #it{y} < 4");
    text1->SetTextSizePixels(18);
    text1->Draw("SAME");

    TLatex * text2 = new TLatex (0.5,0.375,Form("%i-%i %%",minCent,maxCent));
    text2->SetTextSizePixels(18);
    text2->Draw("SAME");

    //TF1 file
    TFile *inputFile = new TFile("$LOWPT/macro/ResultFiles/FitFunctionsAcceptanceEfficiency.root");
    TDirectory *inputList ;
    std::vector<std::vector<Double_t>> *vect;
    std::vector<std::vector<Double_t>> vectParam;
    Double_t chi2;

    TString sRange;
    TString sTest;

    //Tests
    int nTests;
    std::vector<eFunction> fFit;

    Int_t iCent; 
    if(minCent == 30 && maxCent == 50) iCent = 2;
    else if(minCent == 10 && maxCent == 30) iCent = 1; 
    else iCent = 0 ;

    nTests = 2;
    fFit={kPol3, kRatioLevy};

    Int_t color[2]={kAzure-1, kRed+1};
    Int_t style[2]={7, 10};

    TLegend* legend = new TLegend(0.4,0.6,0.89,0.89);
    legend->SetBorderSize(0);
    // legend->SetHeader("J/#psi R_{AA}","C");
    TString sLegend;

    TF1* currentFunction;

    for(int j=0; j<nTests; j++)
    {
        sRange = SetRangeFunction(-4,-2.5, minCent, maxCent); 
        sTest = SetNameTest(kAccEff, fFit[j], kTRUE, kFALSE, 0, 10);
        

        inputList = (TDirectory*) inputFile->Get(Form("%s_%s", sRange.Data(), sTest.Data()));
        inputList->GetObject("parameters", vect);
        vectParam = *vect;
    
        Int_t nPar = vectParam[0].size();
        Double_t parameters[nPar-1];
        Double_t parametersErrors[nPar-1];

        for(int i=0; i<vectParam[0].size()-1; i++)
        {
            parameters[i]=vectParam[0][i];
            parametersErrors[i]=vectParam[1][i];

        }
        chi2=vectParam[0][vectParam[0].size()-1];
        // printf("chi2/NDF = %.3f\n", chi2);

        if (fFit[j] == kRatioLevy)
        {
            currentFunction = GetRatioLevyFunction(0.12, minFit, maxFit);
        }
        else if (fFit[j] == kPol3)
        {
            currentFunction = GetPol3Function(0.15, minFit, maxFit);
        }
        
        currentFunction->SetParameters(parameters);
        currentFunction->SetParErrors(parametersErrors);
        // currentFunction = (TF1*)inputList->Get("function");
        currentFunction->SetLineColor(color[j]);
        currentFunction->SetLineStyle(style[j]);
        currentFunction->SetLineWidth(3);

        currentFunction->Draw("SAME");

        switch (fFit[j])
        {
            case kPol3:
            sLegend = "3rd degree polynomial";
            break;
            case kRatioLevy:
            sLegend = "Ratio of Levy functions";
            break;
        }
        legend->AddEntry(currentFunction, Form("%s, [0-10] GeV/#it{c}, #chi^{2}/NDF = %.2f", sLegend.Data(),chi2), "l");
        
    }
    legend->Draw();

 
}