/*
 *  DrawRaaVScent.C
 *
 *  Created by Ophelie Bugnon on 02/07/20.
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


int minCent[5]={70,50,30,10,0};
int maxCent[5]={90,70,50,30,10};

Double_t nPart[5]={11.35,42.66,108.97,224.95, 357.3};

Double_t minPt[3]={1, 0.3, 0};
Double_t maxPt[3]={2, 1, 0.3};

int nColor[3]={kBlack, kBlue, kRed};
int nMarker[3]={22, 21, 20};


void DrawExcessRaaVScentrality()
{
  TCanvas *cRawNumber = new TCanvas("c1","Raa of Jpsi");
    TH1F *fOption = new TH1F("Raa of Jpsi","", 40, 0, 400);
  	fOption->GetYaxis()->SetTitle("J/#psi R_{AA}");
    // fOption->GetYaxis()->SetRangeUser(0., 1.4);  
  	fOption->GetXaxis()->SetTitle("<N_{part}>");
  	fOption->SetMaximum(12);
    gStyle->SetOptFit(0);
    gStyle->SetOptStat(0);
	fOption->Draw();

    //legend
    TLegend* dlegend = new TLegend(0.1,0.7,0.48,0.9);
    gStyle->SetLegendBorderSize(0);
    // dlegend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
    TLatex * text = new TLatex (1.2,1.32,"Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV");
    text->Draw("SAME");
    text->SetTextSizePixels(20);

    TLatex * text1 = new TLatex (8,1.32,"J/#psi #rightarrow #mu^{+}#mu^{-}, 2.5 < #it{y} < 4");
    text1->SetTextSizePixels(18);
    text1->Draw("SAME");


    for(int i=0; i<3; i++)
    {
      Int_t iGlobalError; 
      if (i==0) iGlobalError = 2;
      else if (i==1) iGlobalError = 1; 
      else if (i==2) iGlobalError = 0; 

    
      TGraphErrors* gr_syst = new TGraphErrors(5);
      
      TGraphErrors* gr_stat = new TGraphErrors(5);

      TGraphErrors *gr_global = new TGraphErrors(1);
      
      for (int j=0; j<5; j++)
      {
        gr_syst->SetPoint(j, nPart[j], RaaResultsVSCent[SetRangeValue(kRaa, -4, -2.5, minPt[i], maxPt[i], minCent[j], maxCent[j])][kRaaValue]);
        cout << RaaResultsVSCent[SetRangeValue(kRaa, -4, -2.5, minPt[i], maxPt[i], minCent[j], maxCent[j])][kRaaValue] << endl;
        gr_syst->SetPointError(j,2, RaaResultsVSCent[SetRangeValue(kRaa, -4, -2.5, minPt[i], maxPt[i], minCent[j], maxCent[j])][kRaaSyst]);

        gr_stat->SetPoint(j, nPart[j], RaaResultsVSCent[SetRangeValue(kRaa, -4, -2.5, minPt[i], maxPt[i], minCent[j], maxCent[j])][kRaaValue]);
        gr_stat->SetPointError(j,0, RaaResultsVSCent[SetRangeValue(kRaa, -4, -2.5, minPt[i], maxPt[i], minCent[j], maxCent[j])][kRaaStat]); 
      }

      gr_syst->SetLineColor(nColor[i]);
      gr_syst->SetFillStyle(0);
      gr_syst->Draw("5same");

      gr_stat->SetMarkerColor(nColor[i]);
      gr_stat->SetMarkerSize(0.5);
      gr_stat->SetMarkerStyle(nMarker[i]);
      gr_stat->SetLineColor(nColor[i]);
      gr_stat->SetFillColor(nColor[i]);
      gr_stat->SetFillStyle(0);
      gr_stat->Draw("psame");

      gr_global->SetPoint(0,395-i*10,1);
      gr_global->SetPointError(0,5,errSystUncorrVScent_Raa[iGlobalError]/100);
      gr_global->SetLineColorAlpha(nColor[i],0);
      gr_global->SetFillColorAlpha(nColor[i], 0.75);
      gr_global->SetFillStyle(1001);
      gr_global->Draw("2same");
  
      dlegend->AddEntry(gr_stat,Form("%.1f < p_{T} < %.1f GeV/c",minPt[i], maxPt[i]) ,"lep");

    }
    TF1 *mean = new TF1("mean","[0]",0,450);
    mean->SetParameter(0, 1);
    mean->SetLineColor(1);
    mean->Draw("SAME");    
    dlegend->Draw();
}
