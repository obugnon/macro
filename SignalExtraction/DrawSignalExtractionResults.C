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

#include </Users/obugnon/Documents/ALICE/AnalyseJPsi/macro/ResultsSignalExtraction.C>
#include </Users/obugnon/Documents/ALICE/AnalyseJPsi/macro/SetRangeAndNameTest.C>

void ResultsLowPt()
{
  Double_t minY = -4;
  Double_t maxY = -2.5;

  Double_t xRange[5]={5, 20, 40, 60, 80};
  Double_t dxRange[5]={5, 10, 10, 10, 10};

  Double_t ptRange[4]={0, 0.3, 1, 8};

  TCanvas *cCrossSection = new TCanvas("c1","Raw number of Jpsi");
  TH1F *fOption = new TH1F("Cross section distribution","", 9, 0, 90);
  fOption->GetYaxis()->SetTitle("Nb of J/#psi");
  fOption->GetXaxis()->SetTitle("cent (%)");
  fOption->SetMaximum(350000);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  fOption->Draw();

  Int_t color[3]={46, 36, 9};
  Int_t marker[3]={22, 21, 20};

  //Pour le remplissage de l'histogramme
  TGraphErrors* grCurrent_syst;
  TGraphErrors* grCurrent_stat;

  std::vector<Double_t> vectNJpsi;
  std::vector<Double_t> nJpsi; 
  std::vector<Double_t> errJpsi_stat;
  std::vector<Double_t> errJpsi_syst;

  Double_t* tJpsi;
  Double_t* terrJpsi_stat;
  Double_t* terrJpsi_syst;

  TLegend* legend = new TLegend(0.1,0.7,0.48,0.9);
  legend->AddEntry((TObject*)0, "Pb-Pb collisions #sqrt{s_{NN} }= 5.02 TeV", "");

  for (int i=0; i<3; i++)
  {
    nJpsi.clear();
    errJpsi_stat.clear();
    errJpsi_syst.clear();
    for (int j = 0; j < 5; j++)
    {
      vectNJpsi.clear();
      vectNJpsi = SignalResultsLowPt[SetRangeValue(kRawYield, minY, maxY, ptRange[i], ptRange[i+1], (Int_t)(xRange[j]-dxRange[j]), (Int_t)(xRange[j]+dxRange[j]))]; 
      nJpsi.push_back(vectNJpsi[kJpsiMean]);
      errJpsi_stat.push_back(vectNJpsi[kJpsiStatError]);
      errJpsi_syst.push_back(vectNJpsi[kJpsiSysError]);
    }
    
    tJpsi = &nJpsi[0];
    terrJpsi_stat = &errJpsi_stat[0];
    terrJpsi_syst = &errJpsi_syst[0];

    grCurrent_syst = new TGraphErrors(5,xRange,tJpsi,dxRange,terrJpsi_syst);
    grCurrent_syst->SetLineColor(color[i]);
    grCurrent_syst->SetFillColor(color[i]);
    grCurrent_syst->SetFillStyle(0);
    grCurrent_syst->Draw("5same");

    grCurrent_stat = new TGraphErrors(5,xRange,tJpsi,0,terrJpsi_stat);
    grCurrent_stat->SetMarkerColor(color[i]);
    grCurrent_stat->SetMarkerStyle(marker[i]);
    grCurrent_stat->SetLineColor(color[i]);
    grCurrent_stat->Draw("psame");

    legend->AddEntry(grCurrent_stat,Form("%.1f - %.1f GeV/c",ptRange[i], ptRange[i+1]),"lep");
      
  }
  legend->Draw();
  

}
//_______________________________________________________________________________________
//_______________________________________________________________________________________
Double_t ComputeDNDpT(Double_t nJpsi, Double_t dpT)
{
  Double_t result = nJpsi/(2*dpT); 
  result = result + 0.5 - (result<0);
  return result;
}

//_______________________________________________________________________________________
//_______________________________________________________________________________________
void ResultsHadro(Bool_t isDraw_dResults)
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

  if(isDraw_dResults)
  {
    cRawNumber = new TCanvas("c1","dN_{J/#psi}/dp_{T}");
    fOption = new TH1F("dN_{J/#psi}/dp_{T}","", 15, 0, 15);
  	fOption->GetYaxis()->SetTitle("dN_{J/#psi}/dp_{T}");
  }
  else
  {
    cRawNumber = new TCanvas("c2","Raw number of Jpsi");
    fOption = new TH1F("Raw number of Jpsi","", 15, 0, 15);
  	fOption->GetYaxis()->SetTitle("Nb of J/#psi");
  }  
  fOption->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  fOption->SetMaximum(350000);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  fOption->Draw();

  //legend
  TLegend* legend = new TLegend(0.1,0.7,0.48,0.9);
  legend->AddEntry((TObject*)0, "Pb-Pb collisions #sqrt{s_{NN} }= 5.02 TeV", "");

  //values
  Int_t iCent; // indice permettant d'acceder à la bonne ligne du tableau
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
      vectNJpsi = SignalResultsHadro[SetRangeValue(kRawYield, minY, maxY, xRange[iCent][j]-dxRange[iCent][j], xRange[iCent][j]+dxRange[iCent][j], minCent[i], maxCent[i])];

      if(isDraw_dResults)
      {
        nJpsi.push_back(ComputeDNDpT(vectNJpsi[kJpsiMean], dxRange[iCent][j]));
        errJpsi_stat.push_back(ComputeDNDpT(vectNJpsi[kJpsiStatError], dxRange[iCent][j]));
        errJpsi_syst.push_back(ComputeDNDpT(vectNJpsi[kJpsiSysError], dxRange[iCent][j]));
      }
      else
      { 
      nJpsi.push_back(vectNJpsi[kJpsiMean]);
      errJpsi_stat.push_back(vectNJpsi[kJpsiStatError]);
      errJpsi_syst.push_back(vectNJpsi[kJpsiSysError]);
      } 
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
  }  
  legend->Draw();
}   