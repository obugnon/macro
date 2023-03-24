/*
 *  FitCentralityLimits.C
 *
 *  Created by Ophelie Bugnon on 25/11/20.
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

Double_t ratioV0[3]={0.99684, 1.00193, 1.01368};
Double_t errratioV0[3]={0.01193, 0.01316, 0.01702};
Double_t centRange[3]={10, 30, 65};
Double_t errcentRange[3]={10, 10, 25};

Double_t Pol1Function(Double_t *x, Double_t *par)
{
	//par[0] = Normalization
	//par[1] = slope
	return  par[0] + par[1]*x[0];
}
TF1* GetPol1Function(Double_t initNorm, Double_t initSlope, Double_t xmin, Double_t xmax, Bool_t isFixedPar1=kFALSE)
{
	TF1* funct = new TF1("Pol1Function", Pol1Function, xmin, xmax, 2);
	funct->SetParameters(initNorm,initSlope);
	funct->SetParNames("Norm", "slope");
	if(isFixedPar1) funct->FixParameter(1, 0);
	return funct;
}
void FitCentralityLimits(Bool_t isIntegral)
{
  TCanvas *cRawNumber = new TCanvas("c1","Jpsi coherent photoproduction cross section");
  TH1F *fOption = new TH1F("Ratio","", 90, 0, 90);
  fOption->GetYaxis()->SetRangeUser(0.95, 1.05);  
  fOption->GetYaxis()->SetTitle("N_{J/#psi}^{V0Mplus05}/N_{J/#psi}^{V0Mminus05}/2");
  fOption->GetXaxis()->SetTitle("centrality (%)");
  gStyle->SetOptFit(111);
  // gStyle->SetOptStat(0);
	fOption->Draw(); 
    
  TGraphErrors* gr_syst = new TGraphErrors(3);
      
  for (int j=0; j<3; j++)
  {
    gr_syst->SetPoint(j, centRange[j], ratioV0[j]);
    gr_syst->SetPointError(j, errcentRange[j], errratioV0[j]);
  } 

  TF1* fFit = GetPol1Function(0.99, 0, 0, 90);

  gr_syst->Draw("PSAME");
  TFitResultPtr fitStatus=gr_syst->Fit("Pol1Function","IEFS");
  Double_t chi2 = fFit->GetChisquare()/fFit->GetNDF();
  Double_t covMatrixStatus = fitStatus->CovMatrixStatus();
  cout << "Prob " << fFit->GetProb()<<endl;
  cout << "Fit status = " << fitStatus << ", cov matrix status = " << covMatrixStatus << " and chi2/NDF = " << chi2 << endl;

    Double_t Integral[5];
    Double_t meanCent[5]={5, 20, 40, 60, 80};
    Double_t edgeCent[5]={0, 30, 50, 70, 90};
    Double_t newCentRange[6]={0, 10, 30, 50, 70, 90};

    for(int i=0; i<5; i++)
    {
        if(isIntegral)
        {
            Integral[i] = fFit->Integral(newCentRange[i], newCentRange[i+1])/(newCentRange[i+1] - newCentRange[i]);
        }
        else Integral[i] = fFit->Eval(edgeCent[i]);
        printf("Centrality limits in %.f-%.f%% is %.5f for uncertainty of %.2f\n", newCentRange[i], newCentRange[i+1], Integral[i], TMath::Abs(1-Integral[i])*100);
    }
    TGraph* gr_newValues;
    if(isIntegral)  gr_newValues = new TGraph(5, meanCent, Integral);
    else gr_newValues = new TGraph(5, edgeCent, Integral);
    gr_newValues->SetMarkerColor(9);
    gr_newValues->SetMarkerSize(3);
    gr_newValues->SetMarkerStyle(3);

    gr_newValues->Draw("PSAME");
}
