/*
 *  DrawRaaVSpT.C
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
#include "TLegend.h"
#include "TROOT.h"
#include "TStyle.h"

#include </Users/obugnon/Documents/ALICE/AnalyseJPsi/macro/ResultFiles/ResultsRaa.C>
#include </Users/obugnon/Documents/ALICE/AnalyseJPsi/macro/SetRangeAndNameTest.C>


void DrawRaaVSpt(Int_t minCent, Int_t maxCent)
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

  Int_t iGlobalError; 
  if (minCent == 70) iGlobalError = 4;
  else if (minCent == 50) iGlobalError = 3; 
  else if (minCent == 30) iGlobalError = 2; 
  else if (minCent == 10) iGlobalError = 1; 
  else if (minCent == 0) iGlobalError = 0; 

  Int_t nBins = xRange[iCent].size();
  
  //Remplissage de l'histogramme
  std::vector<Double_t> vectRaa;
  std::vector<Double_t> nRaa; 
  std::vector<Double_t> errRaa_stat;
  std::vector<Double_t> errRaa_syst;

  for(int i=0; i<nBins; i++)
  {
    vectRaa.clear();
    vectRaa = RaaResultsVSpT[SetRangeValue(kRaa, -4, -2.5, xRange[iCent][i]-dxRange[iCent][i], xRange[iCent][i]+dxRange[iCent][i], minCent, maxCent)];
    nRaa.push_back(vectRaa[kRaaValue]);
    errRaa_stat.push_back(vectRaa[kRaaStat]);
    errRaa_syst.push_back(vectRaa[kRaaSyst]);
  }


  Double_t* txRange;
  txRange = &xRange[iCent][0];
  Double_t* tdxRange;
  tdxRange = &dxRange[iCent][0];
  Double_t* tdxRange_syst;
  tdxRange_syst = &dxRange_syst[iCent][0];

  Double_t* tRaa;
  tRaa = &nRaa[0];
  Double_t* terrRaa_stat;
  terrRaa_stat = &errRaa_stat[0];
  Double_t* terrRaa_syst;
  terrRaa_syst = &errRaa_syst[0];

  TCanvas *cRawNumber = new TCanvas("c1","Raa of Jpsi");
  TH1F *fOption = new TH1F("Jpsi Raa","", 15, 0, 15);
  fOption->GetYaxis()->SetTitle("J/#psi R_{AA}");
  fOption->GetYaxis()->SetRangeUser(0.1, 1.4);  
  fOption->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  fOption->SetMaximum(12);
  // cRawNumber->SetLogy();
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  fOption->Draw();

  TGraphErrors* gr1_syst = new TGraphErrors(nBins,txRange,tRaa,tdxRange_syst,terrRaa_syst);
  gr1_syst->SetMarkerColor(9);
  gr1_syst->SetMarkerStyle(2);
  // gr1_syst->SetMarkerSize(2);
  gr1_syst->SetLineColor(9);
  // gr1_syst->SetLineWidth(2);
  gr1_syst->Draw("5same");

  TGraphErrors* gr1_stat = new TGraphErrors(nBins,txRange,tRaa,tdxRange,terrRaa_stat);
  gr1_stat->SetLineColor(9);
  gr1_stat->SetFillColor(9);
  gr1_stat->SetFillStyle(0);
  // gr1_stat->SetLineWidth(2);
  gr1_stat->Draw("psame");
  

  TF1 *mean = new TF1("mean","[0]",0,16);
  mean->SetParameter(0, 1);
  mean->SetLineColor(1);
  mean->Draw("SAME");
  fOption->GetYaxis()->SetRangeUser(0.1, 1.4);  


  
  TGraphErrors *globalError= new TGraphErrors(1);
  globalError->SetPoint(0,14.5,1);
  globalError->SetPointError(0,0.5,errSystUncorrVSpt_Raa[iGlobalError]/100);
  globalError->SetLineColorAlpha(9,0);
  globalError->SetFillColorAlpha(9, 0.75);
  globalError->SetFillStyle(1001);
  globalError->Draw("2same");

   
  //legend
  TLegend* dlegend = new TLegend(0.4,0.7,0.8,0.89);
  dlegend->SetBorderSize(0);
  dlegend->AddEntry((TObject*)0, "Pb-Pb collisions #sqrt{s_{NN} }= 5.02 TeV", "");
  dlegend->AddEntry((TObject*)0,"-4 < y < -2.5","");
  dlegend->AddEntry((TObject*)0,Form("%i-%i %%",minCent,maxCent),"");
  dlegend->AddEntry((TObject*)0,"","");

  dlegend->AddEntry(gr1_stat,"L_{int} = 750 #mub^{-1}" ,"lep");
  dlegend->AddEntry(gr1_syst,"uncorrelated systematic uncertainty","f");

  
  dlegend->Draw();
  gStyle->SetLegendBorderSize(0);
}
