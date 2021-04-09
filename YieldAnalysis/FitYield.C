/*
 *  FitYield.C
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
#include "TList.h"
#include "TDirectory.h"
#include "TObjArray.h"
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

#include </Users/obugnon/Documents/ALICE/AnalyseJPsi/macro/SetRangeAndNameTest.C>
#include </Users/obugnon/Documents/ALICE/AnalyseJPsi/macro/ResultFiles/ResultsYield.C>

#ifndef JPSI_MASS
#define JPSI_MASS  3.096916
#define PSI2S_MASS 3.686109
#endif

TH1D* GetYieldGraph(Int_t minCent, Int_t maxCent, Bool_t isErrStatOnly, Bool_t isErrSystOnly)
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
    if (minCent == 70 && maxCent == 90) iCent = 2;
    else if (minCent == 50 && maxCent == 70) iCent = 1; 
    else iCent = 0 ; 

    Int_t nBins = xRange[iCent].size();

    Double_t* tLowEdge;
    tLowEdge = &xLowEdge[iCent][0];

    TH1D* graph = new TH1D("J/psi Yiled", "", nBins, tLowEdge);
    graph->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    graph->GetYaxis()->SetTitle("d^{2}Y_{J/#psi}/dy#dotdp_{T}");
    graph->SetMarkerColor(kBlue);
    graph->SetLineColor(kBlue);
    graph->SetFillColor(kBlue);
    graph->SetFillStyle(0);
    

    //Remplissage de l'histogramme
    std::vector<Double_t> vectYield;
    Double_t nYield; 
    Double_t errYield;

    for(int i=0; i<nBins; i++)
    {
        vectYield.clear();
        vectYield = YieldResults[SetRangeValue(kInvariantYield, -4, -2.5, xRange[iCent][i]-dxRange[iCent][i], xRange[iCent][i]+dxRange[iCent][i], minCent, maxCent)];
        nYield=vectYield[kYieldValue];
        if(isErrStatOnly && !isErrSystOnly) errYield=vectYield[kYieldStat];
        else if(isErrSystOnly && !isErrStatOnly) errYield=vectYield[kYieldSyst];
        else errYield=TMath::Sqrt(TMath::Power((vectYield[kYieldStat]/vectYield[kYieldValue]),2) + TMath::Power((vectYield[kYieldSyst]/vectYield[kYieldValue]),2))*vectYield[kYieldValue];
        graph->SetBinContent(i+1, nYield);
        graph->SetBinError(i+1, errYield);
    }

    return graph;
}


std::vector<Double_t> FitYield(eFunction fitFunction, Int_t minCent, Int_t maxCent, Bool_t isErrStatOnly, Bool_t isErrSystOnly, Double_t minPt, Double_t maxPt, Bool_t isSaved, Double_t minFit, Double_t maxFit)
{
    TString sRange = SetRangeFunction(-4, -2.5, minCent, maxCent); 
    TString sTest = SetNameTest(kInvariantYield, fitFunction, isErrStatOnly, isErrSystOnly, minFit, maxFit);

    TCanvas *cYield = new TCanvas(Form("%d-%d%%_%s", minCent, maxCent, sTest.Data()),"J/#psi Yield");
    TH1D* graphYield = GetYieldGraph(minCent, maxCent, isErrStatOnly, isErrSystOnly);
    graphYield->Draw();
    gStyle->SetOptFit(1111);
    graphYield->SetTitle(Form("%d-%d%%_%s", minCent, maxCent, sTest.Data()));


    //Fit
    TF1* fFit=nullptr; 
    TFitResultPtr fitStatus;
    Int_t nPar=0;

    if (fitFunction == kPowLawFree)
    {
        nPar = 4;
        fFit = GetPowerLawFunction(0.025, 0, 15, kFALSE);
        fitStatus=graphYield->Fit("PowerLawFunction","EIS","", minFit, maxFit);
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
    
    Double_t N_Yield = fFit->Integral(minPt, maxPt)/(maxPt-minPt); 
    Double_t Err_Yield = fFit->IntegralError(minPt, maxPt, &params[0], &covmat[0][0], 0.1)/(maxPt-minPt);
    cout << "Yield in " << minPt << "-" << maxPt << " is " << N_Yield << " pm " << Err_Yield << endl;

    Double_t meanPt = fFit->Mean(0.,maxFit);
    Double_t meanPtSquare = fFit->Moment(2,0.,maxFit);

    cout << "Mean pT " << meanPt << " and mean pt square " << meanPtSquare << endl;

    //Draw extrapolated function
    fFit->SetLineStyle(2);
    fFit->Draw("SAME");
    //Save Histogram
    if(isSaved) cYield->SaveAs(".pdf");
    
    //SetResults
    std::vector<Double_t> fitResult;
    fitResult.push_back(N_Yield);
    fitResult.push_back(Err_Yield);
    fitResult.push_back(meanPt);
    fitResult.push_back(0);
    fitResult.push_back(meanPtSquare);
    fitResult.push_back(0);
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

    if (isSaved)
    {
        TFile *outputFile = new TFile("$LOWPT/macro/ResultFiles/FitFunctionsYield.root", "UPDATE");
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
    
    return fitResult;
}