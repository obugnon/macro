/*
 *  DrawSignalExtractionResults.C
 *
 *  Created by Ophelie Bugnon on 05/03/20.
 *
 */
#include <vector>
#include <iostream>
#include <algorithm>
#include <stdio.h>

#include "TFile.h"
#include "TH1.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TString.h"
#include "TMath.h"
#include "TLatex.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TLegend.h"


#include </Users/obugnon/Documents/ALICE/AnalyseJPsi/macro/ResultFiles/ResultsYield.C>
#include </Users/obugnon/Documents/ALICE/AnalyseJPsi/macro/SetRangeAndNameTest.C>


void DrawYieldResults()
{
  Double_t minY = -4;
  Double_t maxY = -2.5;

  std::vector<std::vector<Double_t>> xRange;
  xRange.push_back({0.15, 0.475, 0.825, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 11, 13.5});
  xRange.push_back({0.15, 0.475, 0.825, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 9, 11, 13.5});
  xRange.push_back({0.15, 0.475, 0.825, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 9, 11});

  std::vector<std::vector<Double_t>> dxRange;
  dxRange.push_back({0.15, 0.175, 0.175, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1, 1.5});
  dxRange.push_back({0.15, 0.175, 0.175, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1, 1, 1.5});
  dxRange.push_back({0.15, 0.175, 0.175, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1, 1});

  //initializers
  TCanvas *cRawNumber;
  TH1F *fOption;
  TGraphErrors* grCurrent_syst;
  TGraphErrors* grCurrent_stat;
  Int_t color[5]{12, 9, 32, 42, 46};
  Int_t marker[5]{2, 23, 22, 21, 20};

  cRawNumber = new TCanvas("c2","Raw number of Jpsi");
  fOption = new TH1F("Raw number of Jpsi","", 15, 0, 15);
  fOption->GetYaxis()->SetTitle("d^{2}Y_{J/#psi} /(dp_{T}dy)  (GeV/c)^{-1}");

  fOption->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  fOption->SetMaximum(0.05);
  fOption->SetMinimum(0.0000001);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  fOption->Draw();

  //legend
  TLegend* legend = new TLegend(0.1,0.7,0.48,0.9);
  TLatex * text = new TLatex (8,0.018,"ALICE, Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV");
  text->Draw("SAME");
  text->SetTextSizePixels(22);

  TLatex * text1 = new TLatex (8,0.015,"J/#psi #rightarrow #mu^{+}#mu^{-}, 2.5 < #it{y} < 4");
  text1->SetTextSizePixels(20);
  text1->Draw("SAME");

  TLatex * text2 = new TLatex (1,0.0015,"Centrality");
  text2->SetTextSizePixels(20);
  text2->Draw("SAME");
  
  
  //values
  Int_t iCent;
  Int_t nBins;
  Double_t minCent[5]={0, 10, 30, 50, 70}; 
  Double_t maxCent[5]={10, 30, 50, 70, 90}; 

  std::vector<Double_t> vectNJpsi;
  std::vector<Double_t> nJpsi; 
  std::vector<Double_t> errJpsi_stat;
  std::vector<Double_t> errJpsi_syst;

  Double_t* tJpsi;
  Double_t* terrJpsi_stat;
  Double_t* terrJpsi_syst;
  Double_t* txRange;
  Double_t* tdxRange;

  //fits
  TFile *outputFile = new TFile("$LOWPT/macro/ResultFiles/FitFunctionsYield.root", "UPDATE");
  TDirectory *outputList ;
  std::vector<std::vector<Double_t>> *vect;
  std::vector<std::vector<Double_t>> parameters;
  eFunction fFit=kPowLawFree;
  TString sRange;
  TString sTest;
  TF1* currentFunction;
  TF1* extrapolatedFunction;

  for (int i=0; i<5; i++)
  {
    if (minCent[i] == 70 && maxCent[i] == 90) iCent = 2;
    else if (minCent[i] == 50 && maxCent[i] == 70) iCent = 1; 
    else iCent = 0 ;
    nBins = xRange[iCent].size();

    nJpsi.clear();
    errJpsi_stat.clear();
    errJpsi_syst.clear();
    for (int j = 0; j < nBins; j++)
    {
      vectNJpsi = YieldResults[SetRangeValue(kInvariantYield, minY, maxY, xRange[iCent][j]-dxRange[iCent][j], xRange[iCent][j]+dxRange[iCent][j], minCent[i], maxCent[i])];

      nJpsi.push_back(vectNJpsi[kYieldValue]);
      errJpsi_stat.push_back(vectNJpsi[kYieldStat]);
      errJpsi_syst.push_back(vectNJpsi[kYieldSyst]);

    }

    tJpsi = &nJpsi[0];
    terrJpsi_stat = &errJpsi_stat[0];
    terrJpsi_syst = &errJpsi_syst[0];
    txRange = &xRange[iCent][0];
    tdxRange = &dxRange[iCent][0];

    grCurrent_syst = new TGraphErrors(nBins,txRange,tJpsi,tdxRange,terrJpsi_syst);
    grCurrent_syst->SetLineColor(color[i]);
    grCurrent_syst->SetFillColor(color[i]);
    grCurrent_syst->SetFillStyle(0);
    grCurrent_syst->Draw("5same");

    grCurrent_stat = new TGraphErrors(nBins,txRange,tJpsi,0,terrJpsi_stat);
    grCurrent_stat->SetMarkerColor(color[i]);
    grCurrent_stat->SetMarkerStyle(marker[i]);
    grCurrent_stat->SetLineColor(color[i]);
    grCurrent_stat->Draw("psame");

    legend->AddEntry(grCurrent_stat,Form("%.f-%.f %%",minCent[i], maxCent[i]),"lep");

    sRange = SetRangeFunction(-4,-2.5, minCent[i], maxCent[i]); 
    cout << sRange.Data() << endl;
    sTest = SetNameTest(kInvariantYield, fFit, kFALSE, kFALSE, 1, 15, 0);
    cout << sTest.Data() << endl;

    outputList = (TDirectory*) outputFile->Get(Form("%s_%s", sRange.Data(), sTest.Data()));
    outputList->GetObject("parameters", vect);
    parameters = *vect;

    currentFunction = (TF1*)outputList->Get("function");
    currentFunction->SetLineColor(color[i]);
    currentFunction->SetLineStyle(1);
    currentFunction->SetLineWidth(3);
    currentFunction->SetRange(1, 15);
    currentFunction->Draw("SAME");
    extrapolatedFunction = (TF1*)outputList->Get("function");
    extrapolatedFunction->SetLineColor(color[i]);
    extrapolatedFunction->SetLineStyle(7);
    extrapolatedFunction->SetLineWidth(3);
    extrapolatedFunction->SetRange(0, 1);
    extrapolatedFunction->Draw("SAME");
  }  
  legend->Draw();
}   