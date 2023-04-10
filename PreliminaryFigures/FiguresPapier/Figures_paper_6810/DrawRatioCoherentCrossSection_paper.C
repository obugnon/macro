#include <vector>
#include <iostream>
#include <algorithm>
#include <stdio.h>

#include "TFile.h"
#include "TH1.h"
#include "TF1.h"
#include "TMultiGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TString.h"
#include "TMath.h"
#include "TLatex.h"
#include "TFitResult.h"

#include </Users/obugnon/Documents/ALICE/AnalyseJPsi/macro/SetRangeAndNameTest.C>
#include </Users/obugnon/Documents/ALICE/AnalyseJPsi/macro/ResultFiles/ResultsCoherentCrossSection.C>
#include </Users/obugnon/Documents/ALICE/AnalyseJPsi/macro/figTemplate.C>


int minCent[5]={0, 10, 30,50,70};
int maxCent[5]={10, 30, 50,70,90};
Double_t cent[5]={5, 20, 40, 60, 80};
Double_t dcent[5]={5, 10, 10, 10, 10};


//Npart 5.02 TeV
Double_t nPart5[5]={11.35,42.66,108.97,224.95, 357.3};
Double_t nPart2[5]={11.41, 42.66, 108.4, 222.9, 354.7};
int minCentInverse[5]={70,50,30, 10, 0};
int maxCentInverse[5]={90,70,50, 30, 10};
Double_t dnPart5[5]={4, 4, 4, 4, 5};
Double_t dnPart2[5]={5, 5, 5, 5, 5};

// Preferred colors and markers
// const Int_t fillColors[] = {kGray+1, kAzure-4,  kRed-6, kGreen-8, kViolet-8, kOrange-9}; // for syst bands
// const Int_t colors[] = {kBlack, kAzure-1, kRed+1, kGreen+2, kViolet-6, kOrange-3};
// const Int_t markers[] = {kFullCircle, kFullSquare,kOpenCircle,kOpenSquare,kFullDiamond,kOpenDiamond};

//Values from paper on excess in peripheral collisions at 2.76TeV
Double_t fI=0.14;
Double_t fD=0.10;
Double_t errfD=0.06; //absolute uncertainty 
Double_t errfI_plus=0.16; //absolute uncertainty 
Double_t errfI_minus=0.05; //absolute uncertainty

void DrawCoherentCrossSectionForPreliminary(Bool_t isGlobalSystematicUncluded, Bool_t isModels, Int_t nScenario, Bool_t isUPC)
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
//   cfig->SetLogy();
  // Set Titles etc..


  Double_t minXaxis = 29;
  Double_t maxXaxis;
  if (isUPC) maxXaxis = 120;
  else  maxXaxis = 91;
  Double_t minYaxis = 0;
  Double_t maxYaxis = 8;
  TH1 * h = cfig->DrawFrame(minXaxis,minYaxis,maxXaxis,maxYaxis);

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
  const char * texCent = "Centrality (%)";
  // Raa
  const char *  texRaa="#it{R}_{AA}";
  const char * texCrossSection="d#sigma_{coh}^{J/#psi}/dy (#sqrt{s_{NN}} = 5.02 TeV) / d#sigma_{coh}^{J/#psi}/dy (#sqrt{s_{NN}} = 2.76 TeV) ";

  // Set titles
  h->SetXTitle(texCent);
  // Please be consistent on the y label
  h->SetYTitle(texCrossSection);

  // Draw the logo   
  //  0: Just "ALICE" (for final data), to be added only if ALICE does not appear otherwise (e.g. in the legend)
  //  >0: ALICE Preliminary
//   DrawLogo(0, 0.59, 0.81);

  // You should always specify the colliding system
  // NOTATION: pp, p-Pb, Pb-Pb. 
  // Don't forget to use #sqrt{s_{NN}} for p-Pb and Pb-Pb
  // You can change the position of this with

  TLatex * text = new TLatex (minXaxis+7,maxYaxis-1,"ALICE, Pb-Pb");
  text->SetTextFont(42);
  // text->SetTextSize(0.03652174);
  text->SetLineWidth(2);
  text->Draw("SAME");
  text->SetTextSizePixels(24);

  TLatex * text2 = new TLatex (minXaxis+7,maxYaxis-2,"J/#psi #rightarrow #mu^{+}#mu^{-}, 2.5 < #it{y} < 4");
  text2->SetTextFont(42);
  // text2->SetTextSize(0.03478261);
  text2->SetLineWidth(2);
  text2->Draw();
  text->SetTextSizePixels(22);
  TLatex * text3 = new TLatex (minXaxis+7,maxYaxis-2.5,"#it{p}_{T} < 0.3 GeV/#it{c}");
  text3->SetTextFont(42);;
  // text3->SetTextSize(0.03478261);
  text3->SetLineWidth(2);
  text->SetTextSizePixels(22);
  text3->Draw("SAME");

  TLatex * text4 = new TLatex (minXaxis+18,maxYaxis-1,"Models");
  text4->SetTextFont(42);;
  text4->SetTextSize(0.03478261);
  text4->SetLineWidth(2);
  // text4->Draw("SAME");
 
 //Legend, if needed
  TLegend * leg = new TLegend( 0.56,  0.48,  0.86, 0.66);
  TLegend* dlegend2 = new TLegend(0.45,0.586,0.8,0.786);


//Models 5 TeV
    Double_t nValues;
    if(isUPC) nValues=4;
    else  nValues=3;
    TGraphErrors* gr_IIM = new TGraphErrors(3);
    TGraphErrors* gr_GBW = new TGraphErrors(3);
    TGraphErrors* gr_GG = new TGraphErrors(nValues-2);
    TGraphErrors* gr_GS = new TGraphErrors(nValues-2);
    TGraphErrors* gr_Glaub = new TGraphErrors(3);
    TGraphErrors* gr_Geom = new TGraphErrors(3);

  //-------------------------------------------------------------------
  if(isModels)
  {
    Double_t nCoh_IIM, nCoh_GBW, nCoh_GG, nCoh_GS, nCoh_Glaub, nCoh_Geom;


    for (int j=2; j<2+nValues; j++)
    {
      nCoh_IIM = Models5TeV[Form("Ducati_IIM%d_%d_%d", nScenario,  minCent[j], maxCent[j])]/Models2TeV[Form("Ducati_IIM%d_%d_%d", nScenario,  minCent[j], maxCent[j])]; 
      nCoh_GBW = Models5TeV[Form("Ducati_GBW%d_%d_%d", nScenario,  minCent[j], maxCent[j])]/Models2TeV[Form("Ducati_GBW%d_%d_%d", nScenario,  minCent[j], maxCent[j])];

      cout << "Rapport pour GBW = " << nCoh_GBW << endl;
      cout << "Rapport pour IIM = " << nCoh_IIM << endl;

      gr_IIM->SetPoint(j, cent[j], nCoh_IIM);
      gr_IIM->SetPointError(j,2.5, 0);

      gr_GBW->SetPoint(j, cent[j], nCoh_GBW);
      gr_GBW->SetPointError(j,2.5, 0);

      nCoh_Glaub = Models5TeV[Form("Klusek_Glauber_%d_%d",  minCent[j], maxCent[j])]/Models2TeV[Form("Klusek_Glauber_%d_%d",  minCent[j], maxCent[j])];
      nCoh_Geom = Models5TeV[Form("Klusek_Geometrical_%d_%d",  minCent[j], maxCent[j])]/Models2TeV[Form("Klusek_Geometrical_%d_%d",  minCent[j], maxCent[j])];

      gr_Glaub->SetPoint(j, cent[j], nCoh_Glaub);
      gr_Glaub->SetPointError(j,2.5, 0);

      gr_Geom->SetPoint(j, cent[j], nCoh_Geom);
      gr_Geom->SetPointError(j,2.5, 0);

      // if(j==4 || j==5)
      if(j==4)
      {
        nCoh_GS = Models5TeV[Form("Guillermo_GS_%d_%d",  minCent[j], maxCent[j])]/Models2TeV[Form("Guillermo_GS_%d_%d",  minCent[j], maxCent[j])];
        nCoh_GG = Models5TeV[Form("Guillermo_GG_%d_%d",  minCent[j], maxCent[j])]/Models2TeV[Form("Guillermo_GG_%d_%d",  minCent[j], maxCent[j])];
      
        gr_GG->SetPoint(j, cent[j], nCoh_GG);
        gr_GG->SetPointError(j,2.5, 0);

        gr_GS->SetPoint(j, cent[j], nCoh_GS);
        gr_GS->SetPointError(j,2.5, 0);
      }  
      if(j==5)
      {
        nCoh_GS = Models5TeV["Guillermo_GS_UPC"]/Models2TeV["Guillermo_GS_UPC"];
        nCoh_GG = Models5TeV["Guillermo_GG_UPC"]/Models2TeV["Guillermo_GG_UPC"];
      
        gr_GG->SetPoint(j, 110, nCoh_GG);
        gr_GG->SetPointError(j,2.5, 0);

        gr_GS->SetPoint(j, 110, nCoh_GS);
        gr_GS->SetPointError(j,2.5, 0);
      }  

    }
    //settings IIM
    gr_IIM->SetLineWidth(3);
    gr_IIM->SetLineColor(colors[2]);
    gr_IIM->SetLineStyle(1);
    gr_IIM->Draw("psame");
    //settings GBW
    gr_GBW->SetLineWidth(3);
    gr_GBW->SetLineColor(colors[3]);
    gr_GBW->SetLineStyle(4);
    gr_GBW->Draw("psame");
    //settings GG
    gr_GG->SetLineWidth(3);
    gr_GG->SetLineColor(colors[1]);
    gr_GG->SetLineStyle(2);
    gr_GG->Draw("psame");
    //settings GS
    gr_GS->SetLineWidth(3);
    gr_GS->SetLineColor(colors[3]);
    gr_GS->SetLineStyle(8);
    // gr_GS->Draw("psame");
    //setting Glauber
    gr_Glaub->SetLineWidth(3);
    gr_Glaub->SetLineColor(colors[5]);
    gr_Glaub->SetLineStyle(1);
    gr_Glaub->Draw("psame");
    //setting Geometric
    gr_Geom->SetLineWidth(3);
    gr_Geom->SetLineColor(colors[4]);
    gr_Geom->SetLineStyle(2);
    // gr_Geom->Draw("psame");
  } 

//My data
//-------------------------------------------------------------------
  TGraphAsymmErrors* gr_syst = new TGraphAsymmErrors(nValues);
  TGraphAsymmErrors* gr_globsyst = new TGraphAsymmErrors(nValues);
  TGraphErrors* gr_stat = new TGraphErrors(nValues);

  Double_t nCoh_5, errCoh_stat_5, errCoh_syst_5, errCoh_glob_5, errCoh_syst_plus_5, errCoh_syst_minus_5;
  Double_t nCoh_276, errCoh_stat_276, errCoh_syst_plus_276, errCoh_syst_minus_276, errCoh_glob_276;
  Double_t nCoh, errCoh_stat, errCoh_syst_plus, errCoh_syst_minus, errCoh_glob_plus, errCoh_glob_minus;

  for (int j=2; j<5; j++)
  {

    nCoh_5 = ResultsCRThisAnalysis[SetRangeValue(kRawYield, -4, -2.5, 0, 0.3, minCent[j], maxCent[j])][kCrossSection]; 
    errCoh_stat_5 = ResultsCRThisAnalysis[SetRangeValue(kRawYield, -4, -2.5, 0, 0.3, minCent[j], maxCent[j])][kCrossSectionStat];
    errCoh_syst_5 = ResultsCRThisAnalysis[SetRangeValue(kRawYield, -4, -2.5, 0, 0.3, minCent[j], maxCent[j])][kCrossSectionSyst];
    errCoh_glob_5 = ResultsCRThisAnalysis[SetRangeValue(kRawYield, -4, -2.5, 0, 0.3, minCent[j], maxCent[j])][kCrossSectionSystCor];

    nCoh_276 = ResultsCRRun1Data[SetRangeValue(kRawYield, -4, -2.5, 0, 0.3, minCent[j], maxCent[j])][kCrossSection]; 
    errCoh_stat_276 = ResultsCRRun1Data[SetRangeValue(kRawYield, -4, -2.5, 0, 0.3, minCent[j], maxCent[j])][kCrossSectionStat];
    errCoh_syst_plus_276 = ResultsCRRun1Data[SetRangeValue(kRawYield, -4, -2.5, 0, 0.3, minCent[j], maxCent[j])][2];
    errCoh_syst_minus_276 = ResultsCRRun1Data[SetRangeValue(kRawYield, -4, -2.5, 0, 0.3, minCent[j], maxCent[j])][3];
    errCoh_glob_276 = ResultsCRRun1Data[SetRangeValue(kRawYield, -4, -2.5, 0, 0.3, minCent[j], maxCent[j])][4];

    errCoh_syst_plus_276 = TMath::Sqrt(TMath::Power(errCoh_syst_plus_276/nCoh_276, 2) - TMath::Power(errfI_minus/(1+fI+fD),2) - TMath::Power(errfD/(1+fI+fD), 2))*nCoh_276;
    errCoh_syst_minus_276 = TMath::Sqrt(TMath::Power(errCoh_syst_minus_276/nCoh_276, 2) - TMath::Power(errfI_plus/(1+fI+fD),2) - TMath::Power(errfD/(1+fI+fD), 2))*nCoh_276;


    nCoh = nCoh_5/nCoh_276;
    errCoh_stat = TMath::Sqrt(TMath::Power(errCoh_stat_5/nCoh_5, 2) + TMath::Power(errCoh_stat_276/nCoh_276, 2))*nCoh;
    
    //Avec séparation de l'erreur systematique
    errCoh_syst_plus = TMath::Sqrt(TMath::Power(errCoh_syst_5/nCoh_5, 2) + TMath::Power(errCoh_syst_plus_276/nCoh_276, 2))*nCoh;
    errCoh_syst_minus = TMath::Sqrt(TMath::Power(errCoh_syst_5/nCoh_5, 2) + TMath::Power(errCoh_syst_minus_276/nCoh_276, 2))*nCoh;
    errCoh_glob_plus = TMath::Sqrt(TMath::Power(errCoh_glob_5/nCoh_5, 2) + TMath::Power(errCoh_glob_276/nCoh_276, 2) + TMath::Power(errfI_minus/(1+fI+fD),2) + TMath::Power(errfD/(1+fI+fD), 2) - TMath::Power(1/100, 2) - TMath::Power(0.5/100, 2) - 2*TMath::Power(1/100, 2))*100;
    errCoh_glob_minus = TMath::Sqrt(TMath::Power(errCoh_glob_5/nCoh_5, 2) + TMath::Power(errCoh_glob_276/nCoh_276, 2) + TMath::Power(errfI_plus/(1+fI+fD),2) + TMath::Power(errfD/(1+fI+fD), 2) - TMath::Power(1/100, 2) - TMath::Power(0.5/100, 2) - 2*TMath::Power(1/100, 2))*100;

    //Sans séparation de l'erreur systematique
    // errCoh_syst_plus = TMath::Sqrt(TMath::Power(errCoh_syst_5/nCoh_5, 2) + TMath::Power(errCoh_syst_plus_276/nCoh_276, 2) + TMath::Power(errCoh_glob_5/nCoh_5, 2) + TMath::Power(errCoh_glob_276/nCoh_276, 2) + TMath::Power(errfI_minus/(1+fI+fD),2) + TMath::Power(errfD/(1+fI+fD), 2) - TMath::Power(1/100, 2) - TMath::Power(0.5/100, 2) - 2*TMath::Power(1/100, 2))*nCoh;

    // errCoh_syst_minus = TMath::Sqrt(TMath::Power(errCoh_syst_5/nCoh_5, 2) + TMath::Power(errCoh_syst_minus_276/nCoh_276, 2) + TMath::Power(errCoh_glob_5/nCoh_5, 2) + TMath::Power(errCoh_glob_276/nCoh_276, 2) + TMath::Power(errfI_plus/(1+fI+fD),2) + TMath::Power(errfD/(1+fI+fD), 2) - TMath::Power(1/100, 2) - TMath::Power(0.5/100, 2) - 2*TMath::Power(1/100, 2))*nCoh;

    gr_syst->SetPoint(j-2, cent[j], nCoh);
    gr_syst->SetPointError(j-2,1.3, 1.3, TMath::Max(errCoh_syst_plus, errCoh_syst_minus), TMath::Max(errCoh_syst_plus, errCoh_syst_minus));

    gr_globsyst->SetPoint(j-2, cent[j], nCoh);
    gr_globsyst->SetPointError(j-2, 1, 1, errCoh_glob_minus/100*nCoh, errCoh_glob_plus/100*nCoh);

    gr_stat->SetPoint(j-2, cent[j], nCoh);
    gr_stat->SetPointError(j-2,dcent[j], errCoh_stat); 
    
    printf("In the centrality %i-%i%% the ratio is %.2f pm %.2f pm + %.3f - %.3f pm + %.3f - %.3f\n", minCent[j], maxCent[j], nCoh, errCoh_stat, errCoh_syst_plus, errCoh_syst_minus, errCoh_glob_plus/100*nCoh, errCoh_glob_minus/100*nCoh);
  }

  if(isUPC)
  {
    nCoh_5 = ResultsUPC["rapidity-4.0--2.5_5.02"][kCrossSection]; 
    errCoh_stat_5 = ResultsUPC["rapidity-4.0--2.5_5.02"][kCrossSectionStat];
    errCoh_syst_plus_5 = ResultsUPC["rapidity-4.0--2.5_5.02"][2];
    errCoh_syst_minus_5 = ResultsUPC["rapidity-4.0--2.5_5.02"][3];

    nCoh_276 = ResultsUPC["rapidity-4.0--2.5_2.76"][kCrossSection]; 
    errCoh_stat_276 = ResultsUPC["rapidity-4.0--2.5_2.76"][kCrossSectionStat];
    errCoh_syst_plus_276 = ResultsUPC["rapidity-4.0--2.5_2.76"][2];
    errCoh_syst_minus_276 = ResultsUPC["rapidity-4.0--2.5_2.76"][3];

    nCoh = nCoh_5/nCoh_276;
    errCoh_stat = TMath::Sqrt(TMath::Power(errCoh_stat_5/nCoh_5, 2) + TMath::Power(errCoh_stat_276/nCoh_276, 2))*nCoh;
    errCoh_syst_plus = TMath::Sqrt(TMath::Power(errCoh_syst_5/nCoh_5, 2) + TMath::Power(errCoh_syst_plus_276/nCoh_276, 2))*nCoh;
    errCoh_syst_minus = TMath::Sqrt(TMath::Power(errCoh_syst_5/nCoh_5, 2) + TMath::Power(errCoh_syst_minus_276/nCoh_276, 2))*nCoh;

    gr_syst->SetPoint(4, 110, nCoh);
    gr_syst->SetPointError(4, 1, 1, errCoh_syst_minus, errCoh_syst_plus);

    gr_stat->SetPoint(4, 110, nCoh);
    gr_stat->SetPointError(4,10, errCoh_stat); 
  }

  if(isGlobalSystematicUncluded)
  {
    gr_globsyst->SetFillColorAlpha(fillColors[0], 0.5);
    gr_globsyst->SetLineColor(fillColors[0]);
    gr_globsyst->SetFillStyle(3001);
    gr_globsyst->Draw("5same");
  }
  gr_syst->SetLineColor(colors[0]);
  gr_syst->SetFillStyle(0);
  gr_syst->Draw("5same");

  gr_stat->SetMarkerColor(colors[0]);
  gr_stat->SetMarkerSize(1);
  gr_stat->SetMarkerStyle(markers[1]);
  gr_stat->SetLineColor(colors[0]);
  gr_stat->SetFillColor(colors[0]);
  gr_stat->SetFillStyle(0);
  gr_stat->Draw("psame");

  if(!isGlobalSystematicUncluded) 
  {
    TLatex * text4 = new TLatex (minXaxis+5,maxYaxis-3,Form("Cent. corr. syst. uncert. = ^{%.1f%%}_{%.1f%%}", errCoh_glob_plus, errCoh_glob_minus));
    text4->Draw("SAME");
    text4->SetTextSizePixels(20);
  }  
  else
  {
    leg->AddEntry(gr_globsyst,"Cent. corr. syst. uncert.","f");
  }

  if(isModels)
  {
    dlegend2->SetNColumns(1);
    dlegend2->AddEntry(gr_stat,"Data","pfe");
    // dlegend2->AddEntry((TObject*)0, "M.B. Gay-Ducati et al., Phys. Rev. D97 (2018) 116013", "");
    // dlegend2->AddEntry((TObject*)0, "M.B. Gay-Ducati et al.", "");
    // dlegend2->AddEntry((TObject*)0, "Phys. Rev. D97 (2018) 116013", "");
    dlegend2->AddEntry(gr_IIM,"IIM","l");
    dlegend2->AddEntry(gr_GBW,"GBW","l");
    // dlegend2->AddEntry((TObject*)0,"","");
    // dlegend2->AddEntry((TObject*)0, "J. Cepila et al., Phys. Rev. C97 (2018) 024901", "");
    // dlegend2->AddEntry((TObject*)0, "J. Cepila et al.", "");
    // dlegend2->AddEntry((TObject*)0, "Phys. Rev. C97 (2018) 024901", "");
    dlegend2->AddEntry(gr_GG,"GG-hs ","l");
    // dlegend2->AddEntry(gr_GS,"GS-hs","l");

    // dlegend2->AddEntry((TObject*)0,"","");
    // dlegend2->AddEntry((TObject*)0, "M. Klusek-Gawenda et al.", "");
    // dlegend2->AddEntry((TObject*)0, "Phys. Rev. C93 (2016) 044912", "");
    // dlegend2->AddEntry(gr_Geom,"Geometrical ","l");
    dlegend2->AddEntry(gr_Glaub,"VDM","l");
  }
  leg->SetFillColor(0);
  leg->SetTextSize(gStyle->GetTextSize()*0.55);
  if(isGlobalSystematicUncluded)leg->Draw();
  dlegend2->SetFillColor(0);
  dlegend2->SetTextSize(gStyle->GetTextSize()*0.55);
  dlegend2->Draw();


  if(isUPC) 
  {
    TLegend * dlegAxis = new TLegend( 0.56,  0.48,  0.86, 0.66);
    dlegAxis->AddEntry((TObject*)0, "   #plus#infty", "");
    dlegAxis->SetFillColor(0);
    dlegAxis->SetTextSize(gStyle->GetTextSize()*4);
    dlegAxis->Draw();
    }
}  