/*
 *  FitCrossSection.C
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

#include "FitFunctions.C"
#include "AliPWGFunc.h"
#include </Users/obugnon/Documents/ALICE/AnalyseJPsi/macro/ResultFiles/ResultsPPCrossSection.C>
#include </Users/obugnon/Documents/ALICE/AnalyseJPsi/macro/SetRangeAndNameTest.C>


#ifndef JPSI_MASS
#define JPSI_MASS  3.096916
#define PSI2S_MASS 3.686109
#endif

Double_t minFit = 0.;
Double_t maxFit = 15.;


TH1D* GetppGraph(Bool_t isErrStatOnly, Bool_t isErrSystOnly)
{
    //ranges pT
    const int nBins = 13;
    Double_t x[nBins]={0.15, 0.65, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 11, 13.5};
    Double_t dx[nBins]={0.15, 0.35, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1, 1.5};
    Float_t xLowEdge[nBins+1]={0., 0.3, 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 12., 15.};

    std::vector<Double_t> vectCRpp;
    Double_t nCrossSect; 
    Double_t errCrossSect;

    TH1D* graph = new TH1D("J/psi differential pp cross section", "", nBins, xLowEdge);
    graph->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    graph->GetYaxis()->SetTitle("d^{2}#sigma_{J/#psi}/d#it{p}_{T}d#it{y} (#mub/(GeV/#it{c}))");
    graph->SetMarkerColor(kBlue);
    graph->SetLineColor(kBlue);
    graph->SetFillColor(kBlue);
    graph->SetFillStyle(0);

    for(int i=0; i<nBins; i++)
    {
        vectCRpp.clear();
        vectCRpp = ppCrossSection[SetRangeValue(kCRpp, -4, -2.5, x[i]-dx[i], x[i]+dx[i])];
        nCrossSect=vectCRpp[kCRppValue];
        if(isErrStatOnly && !isErrSystOnly) errCrossSect=vectCRpp[kCRppStatError];
        else if(isErrSystOnly && !isErrStatOnly) errCrossSect=vectCRpp[kCRppSystError];
        else errCrossSect=TMath::Sqrt(TMath::Power((vectCRpp[kCRppStatError]/vectCRpp[kCRppValue]),2) + TMath::Power((vectCRpp[kCRppSystError]/vectCRpp[kCRppValue]),2))*vectCRpp[kCRppValue];
        graph->SetBinContent(i+1, nCrossSect);
        graph->SetBinError(i+1, errCrossSect);
    }


    return graph;
}

std::vector<Double_t> FitCrossSection(eFunction fitFunction, Bool_t isErrStatOnly, Bool_t isErrSystOnly, Double_t minPt, Double_t maxPt, Bool_t isSaved)
{
    TString sTest = SetNameTest(kCRpp, fitFunction, isErrStatOnly, isErrSystOnly);
    TCanvas *cPPcrosssection = new TCanvas(Form("%s",sTest.Data()),"J/#psi differential pp cross section");
    TH1D* graphPP = GetppGraph(isErrStatOnly, isErrSystOnly);
    graphPP->Draw();
    gStyle->SetOptFit(1111);
    cPPcrosssection->SetLogy();
    graphPP->SetTitle(Form("%s", sTest.Data()));

    //Fit
    TF1* fFit=nullptr; 
    TFitResultPtr fitStatus;
    Int_t nPar=0;
    if (fitFunction == kPowLawFree)
    {
        nPar = 4;
        fFit = GetPowerLawFunction(0.65, minFit, maxFit, kFALSE);
        fitStatus=graphPP->Fit("PowerLawFunction","EIRS");
    }
    else if (fitFunction == kPowLawFixed)
    {
        nPar = 4;
        fFit = GetPowerLawFunction(0.65, minFit, maxFit, kTRUE);
        fitStatus=graphPP->Fit("PowerLawFunction","EIRS");
    }
    else if (fitFunction == kLevy)
    {
        nPar = 4;
        fFit = GetLevyFunction(0.65, JPSI_MASS, minFit, maxFit);
        fitStatus=graphPP->Fit("LevyFunction","EIRS");
    }
    else if (fitFunction == kUA1)
    {
        nPar = 6;
        AliPWGFunc* AliFunc = new AliPWGFunc();
        fFit = AliFunc->GetUA1(JPSI_MASS, 2.70, 17.15, 19.89, 0.66, 5.53, "UA1Function");
        fFit->SetParLimits(1, 0.01, 5);
        fitStatus=graphPP->Fit("UA1Function","EIS",0,15);
    }
    //Quality of the fit 
    Double_t chi2 = fFit->GetChisquare()/fFit->GetNDF();
    Double_t covMatrixStatus = fitStatus->CovMatrixStatus();

    cout << "Fit test " << sTest.Data() << " : fits status = " << fitStatus << ", cov matrix status = " << covMatrixStatus << " and chi2/NDF = " << chi2 << endl;

    //Integration
    Double_t params[nPar];
    fFit->GetParameters(params);
    const Double_t *param_errors=fFit->GetParErrors();
    Double_t covmat[nPar][nPar];
    for (int k=0;k<nPar;k++)
	{
		for( int t=0;t<nPar;t++)
		{ 
		    covmat[k][t]=(fitStatus->GetCovarianceMatrix())(k,t);
        }    
	}
    Double_t N_CrossSection = fFit->Integral(minPt, maxPt)/(maxPt-minPt); 
    Double_t Err_CrossSection = fFit->IntegralError(minPt, maxPt, &params[0], &covmat[0][0], 0.1)/(maxPt-minPt);
    cout << "Cross section in " << minPt << "-" << maxPt << " is " << N_CrossSection << " pm " << Err_CrossSection << endl;

    //SetResults
    std::vector<Double_t> fitResult;
    fitResult.push_back(N_CrossSection);
    fitResult.push_back(Err_CrossSection);
    fitResult.push_back(chi2);
    fitResult.push_back(0);
    fitResult.push_back(fitStatus);
    fitResult.push_back(0);
    fitResult.push_back(covMatrixStatus);
    fitResult.push_back(0);

    //Save function results 
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

    if(isSaved)
    {
        cPPcrosssection->SaveAs(".pdf");
        TFile *outputFile = new TFile("$LOWPT/macro/ResultFiles/FitFunctionsCRpp.root", "UPDATE");
        TDirectory *outputList = (TDirectory*) outputFile->Get(Form("%s", sTest.Data()));
        if(!outputList) 
        {
            outputList = outputFile->mkdir(Form("%s", sTest.Data()));
        }
        outputList->cd();
        fFit->Write("function");
        outputList->WriteObject(&fitParameters,"parameters");
        outputList->WriteObject(&covarianceMatrix,"covarianceMatrix");
        outputFile->Close();
    }

    return fitResult;
}

//___________________________________________________________________________________________________________
//___________________________________________________________________________________________________________
void DrawMultipleFits()
{
    //Data
    TCanvas *cPPcrosssection = new TCanvas("","J/#psi differential pp cross section");
    cPPcrosssection->SetRightMargin(0.03);
    cPPcrosssection->SetTopMargin(0.03);
    TH1D* graphPP = GetppGraph(kFALSE, kFALSE);
    // graphRaa->SetTitle(Form("%d-%d%%", minCent, maxCent));
    graphPP->SetMarkerColor(kBlack);
    graphPP->SetLineColor(kBlack);
    graphPP->SetFillColor(kBlack);
    graphPP->SetMarkerStyle(8);
    graphPP->Draw();
    graphPP->GetYaxis()->SetRangeUser(0.001, 2);
    gStyle->SetOptStat(0000);

    TLatex * text = new TLatex (6,1.3,"pp #sqrt{#it{s}_{NN}} = 5.02 TeV");
    text->Draw("SAME");
    text->SetTextSizePixels(20);

    TLatex * text1 = new TLatex (6,1.18,"J/#psi #rightarrow #mu^{+}#mu^{-}, 2.5 < #it{y} < 4");
    text1->SetTextSizePixels(18);
    text1->Draw("SAME");


    //TF1 file
    TFile *inputFile = new TFile("$LOWPT/macro/ResultFiles/FitFunctionsCRpp.root");
    TDirectory *inputList ;
    std::vector<std::vector<Double_t>> *vect;
    std::vector<std::vector<Double_t>> parameters;
    Double_t chi2;

    TString sRange;
    TString sTest;

    //Tests
    int nTests;
    std::vector<eFunction> fFit;

    nTests = 3;
    fFit={kPowLawFree, kPowLawFixed, kLevy};

    Int_t color[5]={kRed+1, kAzure-1, kGreen+2};

    TLegend* legend = new TLegend(0.4,0.6,0.89,0.89);
    legend->SetBorderSize(0);
    TString sLegend;

    TF1* currentFunction;

    for(int j=0; j<nTests; j++)
    {
        TString sTest = SetNameTest(kCRpp, fFit[j], kFALSE, kFALSE);

        inputList = (TDirectory*) inputFile->Get(Form("%s", sTest.Data()));
        inputList->GetObject("parameters", vect);
        parameters = *vect;
        chi2=parameters[0][parameters[0].size()-1];
        // printf("chi2/NDF = %.3f\n", chi2);

        currentFunction = (TF1*)inputList->Get("function");
        currentFunction->SetLineColor(color[j]);

        currentFunction->Draw("SAME");

        switch (fFit[j])
        {
            case kPowLawFixed:
            sLegend = "Powerlaw D = 2";
            break;
            case kPowLawFree:
            sLegend = "Powerlaw D free";
            break;
            case kLevy:
            sLegend = "Levy function";
            break;
        }
        legend->AddEntry(currentFunction, Form("%s, [0-15] GeV/#it{c}, #chi^{2}/NDF = %.2f", sLegend.Data(), chi2), "l");
        
    }
    legend->Draw();

 
}