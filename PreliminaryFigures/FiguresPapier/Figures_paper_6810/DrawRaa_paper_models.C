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

#include </Users/obugnon/Documents/ALICE/AnalyseJPsi/macro/ResultFiles/ResultsRaa.C>
#include </Users/obugnon/Documents/ALICE/AnalyseJPsi/macro/figTemplate.C>
#include </Users/obugnon/Documents/ALICE/AnalyseJPsi/macro/SetRangeAndNameTest.C>


int minCent[5]={70,50,30,10,0};
int maxCent[5]={90,70,50,30,10};

Double_t nPart[5]={11.35,42.66,108.97,224.95, 357.3};

Double_t minPt[3]={1,0.3,0};
Double_t maxPt[3]={2,1, 0.3};

int nColor[3]={kBlack, kBlue, kRed};
int nMarker[3]={22, 21, 20};
int nFillStyle[3]={3003, 3004, 3005};

Double_t impactParameterRange[26]={0, 0.6, 1.2, 1.8, 2.4, 3, 3.6, 4.2, 4.8, 5.4, 6, 6.6, 7.2, 7.8, 8.4, 9, 9.6, 10.2, 10.8, 11.4, 12, 12.6, 13.2, 13.8, 14.4, 15};

//For Preliminary Figure
void DrawRaa_paper(Bool_t isModels)
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
  cfig->SetBottomMargin(0.15);
  cfig->SetLeftMargin(0.12);
  cfig->SetRightMargin(0.08);
  cfig->SetTopMargin(0.08);
  cfig->SetLogy();
  // Set Titles etc..
  TH1 * h = cfig->DrawFrame(0,0.4,400,12);

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
   //  DrawLogo(0, 0.59, 0.81);

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
    TLatex * text = new TLatex (102,6.16,"ALICE, Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV");
    text->Draw("SAME");
    text->SetTextSizePixels(28);

    TLatex * text2 = new TLatex (102.6,4.79,"J/#psi #rightarrow #mu^{+}#mu^{-}, 2.5 < #it{y} < 4");
    text2->SetTextSizePixels(28);
    text2->Draw("SAME");

    //Legend, if needed
    TLegend * leg = new TLegend( 0.56,  0.48,  0.86, 0.66);
    if(isModels)
      {
        leg->SetNColumns(3);
        // leg->AddEntry((TObject*)0," ","");
        // leg->AddEntry((TObject*)0,"Data" ,"");
        // leg->AddEntry((TObject*)0,"Transport" ,"");

      }  
    
    TF1 *mean = new TF1("mean","[0]",0,450);
    mean->SetParameter(0, 1);
    mean->SetLineColor(1);
    mean->Draw("SAME");  
    
    TGraphErrors* gr_syst[3];
    TGraphErrors* gr_stat[3];
    TGraphErrors* gr_global[3];
    TGraph *gr_filled[3];

    for(int i=0; i<3; i++)
    {
      Int_t iGlobalError; 
      if (i==0) iGlobalError = 2;
      else if (i==1) iGlobalError = 1; 
      else if (i==2) iGlobalError = 0; 
    
      gr_syst[i] = new TGraphErrors(5);
      
      gr_stat[i] = new TGraphErrors(5);

      gr_global[i] = new TGraphErrors(1);

      for (int j=0; j<5; j++)
      {
      //   if(minPt[i]==0 && j > 2) continue;
        gr_syst[i]->SetPoint(j, nPart[j], RaaResultsVSCent[SetRangeValue(kRaa, -4, -2.5, minPt[i], maxPt[i], minCent[j], maxCent[j])][kRaaValue]);
        // cout << RaaResultsVSCent[SetRangeValue(kRaa, -4, -2.5, minPt[i], maxPt[i], minCent[j], maxCent[j])][kRaaValue] << endl;
        gr_syst[i]->SetPointError(j,2, RaaResultsVSCent[SetRangeValue(kRaa, -4, -2.5, minPt[i], maxPt[i], minCent[j], maxCent[j])][kRaaSyst]);

        gr_stat[i]->SetPoint(j, nPart[j], RaaResultsVSCent[SetRangeValue(kRaa, -4, -2.5, minPt[i], maxPt[i], minCent[j], maxCent[j])][kRaaValue]);
        gr_stat[i]->SetPointError(j,0, RaaResultsVSCent[SetRangeValue(kRaa, -4, -2.5, minPt[i], maxPt[i], minCent[j], maxCent[j])][kRaaStat]); 
      }

      gr_syst[i]->SetLineColor(colors[i]);
      gr_syst[i]->SetFillStyle(0);
      gr_syst[i]->Draw("5same");

      gr_stat[i]->SetMarkerColor(colors[i]);
      gr_stat[i]->SetMarkerSize(0.8);
      gr_stat[i]->SetMarkerStyle(markers[i]);
      gr_stat[i]->SetLineColor(colors[i]);
      gr_stat[i]->SetFillColor(colors[i]);
      gr_stat[i]->SetFillStyle(0);
      gr_stat[i]->Draw("psame");

      gr_global[i]->SetPoint(0,393-i*14,1);
      gr_global[i]->SetPointError(0,7,errSystUncorrVScent_Raa[iGlobalError]/100);
      gr_global[i]->SetLineColorAlpha(fillColors[i],0);
      gr_global[i]->SetFillColorAlpha(fillColors[i], 0.75);
      gr_global[i]->SetFillStyle(1001);
      gr_global[i]->Draw("2same");
  
      //for models
      if (isModels)
      {
        Int_t n=26;
        Double_t nPartRange[26];
        for (int j=0; j<n; j++)
        {
          nPartRange[j]= Raa_Prompt_Prediction[SetRangeValue(kRaaModels, -4, -2.5, minPt[2], maxPt[2], 0, 0, impactParameterRange[j])][0];
        }
        TGraph* gr_min = new TGraph(n);
        TGraph* gr_max = new TGraph(n);
        gr_filled[i]= new TGraph(2*n);
        for (int j=0; j<n; j++)
        {
          gr_min->SetPoint(j, nPartRange[j], Raa_Prompt_Prediction[SetRangeValue(kRaaModels, -4, -2.5, minPt[i], maxPt[i], 0, 0, impactParameterRange[j])][1]);
          gr_max->SetPoint(j, nPartRange[j], Raa_Prompt_Prediction[SetRangeValue(kRaaModels, -4, -2.5, minPt[i], maxPt[i], 0, 0, impactParameterRange[j])][2]);
          gr_filled[i]->SetPoint(j, nPartRange[j], Raa_Prompt_Prediction[SetRangeValue(kRaaModels, -4, -2.5, minPt[i], maxPt[i], 0, 0, impactParameterRange[j])][2]);
          gr_filled[i]->SetPoint(j+n, nPartRange[n-j-1], Raa_Prompt_Prediction[SetRangeValue(kRaaModels, -4, -2.5, minPt[i], maxPt[i], 0, 0, impactParameterRange[n-j-1])][1]);
        }
        gr_min->SetLineColorAlpha(colors[i], 0.5);
        gr_max->SetLineColorAlpha(colors[i], 0.5);
        gr_min->SetLineWidth(1);
        gr_max->SetLineWidth(1);
        gr_min->Draw("lsame");
        gr_max->Draw("lsame");
        gr_filled[i]->SetFillStyle(nFillStyle[i]);
        gr_filled[i]->SetLineColorAlpha(colors[i], 0.4);
        gr_filled[i]->SetFillColorAlpha(colors[i], 0.4);
        gr_filled[i]->Draw("f");
        // leg->AddEntry(gr_filled, "Model" ,"f");
      }
      
    }

    if(isModels)
    {
    leg->AddEntry((TObject*)0,"#it{p}_{T} < 0.3 GeV/#it{c}" ,"");
    leg->AddEntry(gr_stat[2],"Data" ,"ep");
    leg->AddEntry(gr_filled[2], "Model" ,"f");
    leg->AddEntry((TObject*)0,"0.3 < #it{p}_{T} < 1 GeV/#it{c}" ,"");
    leg->AddEntry(gr_stat[1],"Data" ,"ep");
    leg->AddEntry(gr_filled[1], "Model" ,"f");
    leg->AddEntry((TObject*)0,"1< #it{p}_{T} < 2 GeV/#it{c}" ,"");
    leg->AddEntry(gr_stat[0],"Data" ,"ep");
    leg->AddEntry(gr_filled[0], "Model" ,"f");
    }
    else
    {
    leg->AddEntry(gr_stat[2],"#it{p}_{T} < 0.3 GeV/#it{c}" ,"ep");
    leg->AddEntry(gr_stat[1],"0.3 < #it{p}_{T} < 1 GeV/#it{c}" ,"ep");
    leg->AddEntry(gr_stat[0],"1< #it{p}_{T} < 2 GeV/#it{c}" ,"ep");
    }
    leg->SetFillColor(0);
    leg->SetTextSize(gStyle->GetTextSize()*0.9);
    leg->Draw();

}