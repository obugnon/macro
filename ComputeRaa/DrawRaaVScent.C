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

#include </Users/obugnon/Documents/ALICE/AnalyseJPsi/macro/ResultsRaa.C>
#include </Users/obugnon/Documents/ALICE/AnalyseJPsi/macro/figTemplate.C>
#include </Users/obugnon/Documents/ALICE/AnalyseJPsi/macro/SetRangeAndNameTest.C>


int minCent[5]={70,50,30,10,0};
int maxCent[5]={90,70,50,30,10};

Double_t nPart[5]={11.35,42.66,108.97,224.95, 357.3};

Double_t minPt[3]={1,0.3,0};
Double_t maxPt[3]={2,1, 0.3};

int nColor[3]={kBlack, kBlue, kRed};
int nMarker[3]={22, 21, 20};


void ResultsHadro()
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
    dlegend->AddEntry((TObject*)0, "Asked for preliminary", "");
    dlegend->AddEntry((TObject*)0, "", "");
    dlegend->AddEntry((TObject*)0, "Pb-Pb collisions #sqrt{s_{NN} }= 5.02 TeV", "");
    dlegend->AddEntry((TObject*)0,"-4 < y < -2.5",""); 

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

//For Preliminary Figure
void ResultsHadro_ForPreliminary()
{
  // Load necessary libraries
  LoadLibs();
  // Set the default style
  SetStyle();

  // Prepare Figure, please stick to the default canvas size(s) unless absolutely necessary in your case
  // Rectangular
  TCanvas *cfig = new TCanvas("cfig", "Alice Figure Template", 800, 600); 
  // Square
  //TCanvas *cfig = new TCanvas("cfig", "Alice Figure Template", 800, 800); 
  cfig->SetLogy();
  // Set Titles etc..
  TH1 * h = cfig->DrawFrame(0,0.4,150,12);

  // === Commonly used x/ titles: ===
    // pt invariant yields
    const char *  texPtX="#it{p}_{T} (GeV/#it{c})";
    const char *  texPtY="1/#it{N}_{ev} 1/(2#pi#it{p}_{T}) d#it{N}/(d#it{p}_{T}d#it{y}) ((GeV/#it{c})^{-2})";
    // mt invariant yields
    const char *  texMtX="#it{m}_{T} (GeV/#it{c}^{2})";
    const char *  texMtY="1/#it{N}_{ev} 1/(2#pi#it{m}_{T}) d#it{N}/(d#it{m}_{T}d#it{y}) ((GeV/#it{c}^{2})^{-2}) "; 
    // Invariant mass with decay products K and pi
    const char *  texMassX="#it{M}_{K#pi} (GeV/#it{c}^{2})";
    const char *  texMassY="d#it{N}/(d#it{M}_{K#pi})";
    // <pt>, npart
    const char * texMeanPt    = "#LT#it{p}_{T}#GT (GeV/#it{c})";
    const char * texMeanNpart = "#LT#it{N}_{part}#GT";
    // Raa
    const char *  texRaa="#it{R}_{AA}";

    // Set titles
    h->SetXTitle(texMeanNpart);
    // Please be consistent on the y label
    h->SetYTitle(texRaa);

    // Draw the logo   
    //  0: Just "ALICE" (for final data), to be added only if ALICE does not appear otherwise (e.g. in the legend)
    //  >0: ALICE Preliminary
    DrawLogo(1, 0.59, 0.81);

    // You should always specify the colliding system
    // NOTATION: pp, p-Pb, Pb-Pb. 
    // Don't forget to use #sqrt{s_{NN}} for p-Pb and Pb-Pb
    // You can change the position of this with
    // TLatex * text = new TLatex (202,6.16,"Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV");
    // text->Draw("SAME");
    // text->SetTextSizePixels(20);

    // TLatex * text2 = new TLatex (202.6,4.79,"J/#psi #rightarrow #mu^{+}#mu^{-}, 2.5 < #it{y} < 4");
    // text2->SetTextSizePixels(20);
    // text2->Draw("SAME");
    TLatex * text = new TLatex (102,6.16,"Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV");
    text->Draw("SAME");
    text->SetTextSizePixels(20);

    TLatex * text2 = new TLatex (102.6,4.79,"J/#psi #rightarrow #mu^{+}#mu^{-}, 2.5 < #it{y} < 4");
    text2->SetTextSizePixels(20);
    text2->Draw("SAME");
  
    //Legend, if needed
    TLegend * leg = new TLegend( 0.56,  0.48,  0.86, 0.66);

    TF1 *mean = new TF1("mean","[0]",0,450);
    mean->SetParameter(0, 1);
    mean->SetLineColor(1);
    mean->Draw("SAME");  
    

    for(int i=0; i<3; i++)
    {
      Int_t iGlobalError; 
      if (i==0) iGlobalError = 2;
      else if (i==1) iGlobalError = 1; 
      else if (i==2) iGlobalError = 0; 
    
      TGraphErrors* gr_syst = new TGraphErrors(5);
      
      TGraphErrors* gr_stat = new TGraphErrors(5);

      TGraphErrors *gr_global = new TGraphErrors(1);
      
      for (int j=0; j<3; j++)
      {
        if(minPt[i]==0 && j > 2) continue;
        gr_syst->SetPoint(j, nPart[j], RaaResultsVSCent[SetRangeValue(kRaa, -4, -2.5, minPt[i], maxPt[i], minCent[j], maxCent[j])][kRaaValue]);
        cout << RaaResultsVSCent[SetRangeValue(kRaa, -4, -2.5, minPt[i], maxPt[i], minCent[j], maxCent[j])][kRaaValue] << endl;
        gr_syst->SetPointError(j,2, RaaResultsVSCent[SetRangeValue(kRaa, -4, -2.5, minPt[i], maxPt[i], minCent[j], maxCent[j])][kRaaSyst]);

        gr_stat->SetPoint(j, nPart[j], RaaResultsVSCent[SetRangeValue(kRaa, -4, -2.5, minPt[i], maxPt[i], minCent[j], maxCent[j])][kRaaValue]);
        gr_stat->SetPointError(j,0, RaaResultsVSCent[SetRangeValue(kRaa, -4, -2.5, minPt[i], maxPt[i], minCent[j], maxCent[j])][kRaaStat]); 
      }

      gr_syst->SetLineColor(colors[i]);
      gr_syst->SetFillStyle(0);
      gr_syst->Draw("5same");

      gr_stat->SetMarkerColor(colors[i]);
      gr_stat->SetMarkerSize(0.8);
      gr_stat->SetMarkerStyle(markers[i]);
      gr_stat->SetLineColor(colors[i]);
      gr_stat->SetFillColor(colors[i]);
      gr_stat->SetFillStyle(0);
      gr_stat->Draw("psame");

      gr_global->SetPoint(0,147-i*6,1);
      gr_global->SetPointError(0,3,errSystUncorrVScent_Raa[iGlobalError]/100);
      gr_global->SetLineColorAlpha(fillColors[i],0);
      gr_global->SetFillColorAlpha(fillColors[i], 0.75);
      gr_global->SetFillStyle(1001);
      gr_global->Draw("2same");
  
      if(minPt[i]==0.3) leg->AddEntry(gr_stat,Form("%.1f < #it{p}_{T} < %.f GeV/#it{c}",minPt[i], maxPt[i]) ,"lp");
      else if(maxPt[i]==0.3) leg->AddEntry(gr_stat,Form("%.f < #it{p}_{T} < %.1f GeV/#it{c}",minPt[i], maxPt[i]) ,"lp");
      else leg->AddEntry(gr_stat,Form("%.f < #it{p}_{T} < %.f GeV/#it{c}",minPt[i], maxPt[i]) ,"lp");

    }
    leg->SetFillColor(0);
    leg->SetTextSize(gStyle->GetTextSize()*0.6);
    leg->Draw();


      
}