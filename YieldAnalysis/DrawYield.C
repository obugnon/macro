/*
 *  DrawYield.C
 *
 *  Created by Ophelie Bugnon on 19/05/20.
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
#include "TF1.h"

#include </Users/obugnon/Documents/ALICE/AnalyseJPsi/macro/ResultsYield.C>
#include </Users/obugnon/Documents/ALICE/AnalyseJPsi/macro/SetRangeAndNameTest.C>


void DrawYieldVSpt(Int_t minCent, Int_t maxCent)
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

  std::vector<std::vector<Double_t>> dxRange_syst;
  dxRange_syst.push_back({0.05, 0.058, 0.058, 0.17, 0.17, 0.17, 0.17, 0.17, 0.17, 0.17, 0.17, 0.17, 0.33, 0.5});
  dxRange_syst.push_back({0.05, 0.058, 0.058, 0.17, 0.17, 0.17, 0.17, 0.17, 0.17, 0.17, 0.33, 0.33, 0.5});
  dxRange_syst.push_back({0.05, 0.058, 0.058, 0.17, 0.17, 0.17, 0.17, 0.17, 0.17, 0.17, 0.33, 0.33});


  Int_t iCent; 
  if (minCent == 70 && maxCent == 90) iCent = 2;
  else if (minCent == 50 && maxCent == 70) iCent = 1; 
  else iCent = 0 ; 

  Int_t nBins = xRange[iCent].size();
  
  //Fill Histo
  std::vector<Double_t> vectYield;
  std::vector<Double_t> nYield; 
  std::vector<Double_t> errYield_stat;
  std::vector<Double_t> errYield_syst;

  for(int i=0; i<nBins; i++)
  {
    vectYield.clear();
    vectYield = YieldResults[SetRangeValue(kInvariantYield, -4, -2.5, xRange[iCent][i]-dxRange[iCent][i], xRange[iCent][i]+dxRange[iCent][i], minCent, maxCent)];
    nYield.push_back(vectYield[kYieldValue]);
    errYield_stat.push_back(vectYield[kYieldStat]);
    errYield_syst.push_back(vectYield[kYieldSyst]);
  }

  Double_t* txRange;
  txRange = &xRange[iCent][0];
  Double_t* tdxRange;
  tdxRange = &dxRange[iCent][0];
  Double_t* tdxRange_syst;
  tdxRange_syst = &dxRange_syst[iCent][0];

  Double_t* tYield;
  tYield = &nYield[0];
  Double_t* terrYield_stat;
  terrYield_stat = &errYield_stat[0];
  Double_t* terrYield_syst;
  terrYield_syst = &errYield_syst[0];

  TCanvas *cRawNumber = new TCanvas("c1","Yield of Jpsi");
  TH1F *fOption = new TH1F("Jpsi Yield","", 15, 0, 15);
  fOption->GetYaxis()->SetTitle("d^{2}Y_{J/#psi}/dy#dotdp_{T} ");
  fOption->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  fOption->SetMaximum(12);
  cRawNumber->SetLogy();
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  fOption->Draw();

  TGraphErrors* gr1_syst = new TGraphErrors(nBins,txRange,tYield,tdxRange_syst,terrYield_syst);
  gr1_syst->SetMarkerColor(9);
  gr1_syst->SetMarkerStyle(2);
  // gr1_syst->SetMarkerSize(2);
  gr1_syst->SetLineColor(9);
  // gr1_syst->SetLineWidth(2);
  gr1_syst->Draw("5same");

  TGraphErrors* gr1_stat = new TGraphErrors(nBins,txRange,tYield,tdxRange,terrYield_stat);
  gr1_stat->SetLineColor(9);
  gr1_stat->SetFillColor(9);
  gr1_stat->SetFillStyle(0);
  // gr1_stat->SetLineWidth(2);
  gr1_stat->Draw("psame");
  
  fOption->GetYaxis()->SetRangeUser(0.000001, 0.1);  
  //legend
  TLegend* dlegend = new TLegend(0.1,0.7,0.48,0.9);
  dlegend->AddEntry((TObject*)0, "Pb-Pb collisions #sqrt{s_{NN} }= 5.02 TeV", "");
  dlegend->AddEntry((TObject*)0,Form("%i-%i %%",minCent,maxCent),"");
  dlegend->AddEntry((TObject*)0,"-4 < y < -2.5","");
  dlegend->Draw();
  gStyle->SetLegendBorderSize(0);
}
