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
int minCentInverse[5]={70,50,30, 10, 0};
int maxCentInverse[5]={90,70,50, 30, 10};

//Npart 2.76 TeV
Double_t nPart2[5]={11.41, 42.66, 108.4, 222.9, 354.7};
// Double_t nPart2[3]={11.41, 42.66, 108.4};
Double_t dnPart2[5]={5.5, 5.5, 5.5, 5.5, 5.5};


// Preferred colors and markers
// const Int_t fillColors[] = {kGray+1, kAzure-4,  kRed-6, kGreen-8, kViolet-8, kOrange-9}; // for syst bands
// const Int_t colors[] = {kBlack, kAzure-1, kRed+1, kGreen+2, kViolet-6, kOrange-3};
// const Int_t markers[] = {kFullCircle, kFullSquare,kOpenCircle,kOpenSquare,kFullDiamond,kOpenDiamond};

void DrawCoherentCrossSectionForPreliminary(Bool_t isGlobalSystematicUncluded, Bool_t isModels)
{
// Load necessary libraries
  LoadLibs();
  // Set the default style
  SetStyle();

  // Prepare Figure, please stick to the default canvas size(s) unless absolutely necessary in your case
  // Rectangular
  TCanvas *cfig = new TCanvas("cfig", "Alice Figure Template", 800, 600); 
  cfig->SetTopMargin(0.05);
  cfig->SetRightMargin(0.04);
  // Square
  //TCanvas *cfig = new TCanvas("cfig", "Alice Figure Template", 800, 800); 
//   cfig->SetLogy();
  // Set Titles etc..
  TH1 * h = cfig->DrawFrame(0,0,400,550);

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
  const char * texCrossSection="d#sigma_{coh}^{J/#psi}/d#it{y} (#mub)";

  // Set titles
  h->SetXTitle(texMeanNpart);
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

  TLatex * text = new TLatex (25,475,"ALICE, Pb-Pb #sqrt{#it{s}_{NN}} = 2.76 TeV");
  text->Draw("SAME");
  text->SetTextSizePixels(24);

  TLatex * text2 = new TLatex (25,396,"J/#psi #rightarrow #mu^{+}#mu^{-}, 2.5 < #it{y} < 4");
  text2->SetTextSizePixels(22);
  text2->Draw("SAME");
  TLatex * text3 = new TLatex (25,342.6,"#it{p}_{T} < 0.3 GeV/#it{c}");
  text3->SetTextSizePixels(22);
  text3->Draw("SAME");

    TLatex * text4 = new TLatex (225,625,"Models");
  // text4->Draw("SAME");
  text4->SetTextSizePixels(22);
 
 //Legend, if needed
  TLegend * leg = new TLegend( 0.56,  0.48,  0.86, 0.66);
  TLegend* dlegend2 = new TLegend(0.58,0.64,0.85,0.90);


//Models 5 TeV
    TGraphErrors* gr_IIM2 = new TGraphErrors(5);
    TGraphErrors* gr_GBW2 = new TGraphErrors(5);
    TGraphErrors* gr_IIM3 = new TGraphErrors(5);
    TGraphErrors* gr_GBW3 = new TGraphErrors(5);
    TGraphErrors* gr_GG = new TGraphErrors(1);
    TGraphErrors* gr_GS = new TGraphErrors(1);
    TGraphErrors* gr_Glaub = new TGraphErrors(5);
    TGraphErrors* gr_Geom = new TGraphErrors(5);

  //-------------------------------------------------------------------
  if(isModels)
  {
    Double_t nCoh_IIM2, nCoh_GBW2, nCoh_IIM3, nCoh_GBW3, nCoh_GG, nCoh_GS, nCoh_Glaub, nCoh_Geom;


    for (int j=0; j<5; j++)
    {
      nCoh_IIM2 = Models2TeV[Form("Ducati_IIM2_%d_%d", minCentInverse[j], maxCentInverse[j])]; 
      nCoh_GBW2 = Models2TeV[Form("Ducati_GBW2_%d_%d",  minCentInverse[j], maxCentInverse[j])];
      nCoh_IIM3 = Models2TeV[Form("Ducati_IIM3_%d_%d", minCentInverse[j], maxCentInverse[j])]; 
      nCoh_GBW3 = Models2TeV[Form("Ducati_GBW3_%d_%d",  minCentInverse[j], maxCentInverse[j])];

      gr_IIM2->SetPoint(j, nPart2[j], nCoh_IIM2);
      gr_IIM2->SetPointError(j,dnPart2[j]*2, 0);

      gr_GBW2->SetPoint(j, nPart2[j], nCoh_GBW2);
      gr_GBW2->SetPointError(j,dnPart2[j]*2, 0);

      gr_IIM3->SetPoint(j, nPart2[j], nCoh_IIM3);
      gr_IIM3->SetPointError(j,dnPart2[j]*2, 0);

      gr_GBW3->SetPoint(j, nPart2[j], nCoh_GBW3);
      gr_GBW3->SetPointError(j,dnPart2[j]*2, 0);

      nCoh_Glaub = Models2TeV[Form("Klusek_Glauber_%d_%d", minCentInverse[j], maxCentInverse[j])]; 
      nCoh_Geom = Models2TeV[Form("Klusek_Geometrical_%d_%d", minCentInverse[j], maxCentInverse[j])];

      gr_Glaub->SetPoint(j, nPart2[j], nCoh_Glaub);
      gr_Glaub->SetPointError(j,dnPart2[j]*2, 0);

      gr_Geom->SetPoint(j, nPart2[j], nCoh_Geom);
      gr_Geom->SetPointError(j,dnPart2[j]*2, 0);

      if(j==0)
      {
        nCoh_GS = Models2TeV[Form("Guillermo_GS_%d_%d",  minCentInverse[j], maxCentInverse[j])];
        nCoh_GG = Models2TeV[Form("Guillermo_GG_%d_%d",  minCentInverse[j], maxCentInverse[j])];
      
        gr_GG->SetPoint(j, nPart2[j], nCoh_GG);
        gr_GG->SetPointError(j,dnPart2[j]*2, 0);

        gr_GS->SetPoint(j, nPart2[j], nCoh_GS);
        gr_GS->SetPointError(j,dnPart2[j]*2, 0);
      }  

    }
    //settings IIM2
    gr_IIM2->SetLineWidth(3);
    gr_IIM2->SetLineColor(colors[2]);
    gr_IIM2->SetLineStyle(7);
    gr_IIM2->Draw("psame");
    //settings IIM3
    gr_IIM3->SetLineWidth(3);
    gr_IIM3->SetLineColor(colors[2]);
    gr_IIM3->SetLineStyle(1);
    gr_IIM3->Draw("psame");
    //settings GBW2
    gr_GBW2->SetLineWidth(3);
    gr_GBW2->SetLineColor(colors[3]);
    gr_GBW2->SetLineStyle(7);
    gr_GBW2->Draw("psame");
    //settings GBW2
    gr_GBW3->SetLineWidth(3);
    gr_GBW3->SetLineColor(colors[3]);
    gr_GBW3->SetLineStyle(1);
    gr_GBW3->Draw("psame");
    //setting Glauber
    gr_Glaub->SetLineWidth(3);
    gr_Glaub->SetLineColor(colors[5]);
    gr_Glaub->SetLineStyle(1);
    gr_Glaub->Draw("psame");
    //setting Geometric
    gr_Geom->SetLineWidth(3);
    gr_Geom->SetLineColor(colors[4]);
    gr_Geom->SetLineStyle(1);
    // gr_Geom->Draw("psame");
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
    
  } 

//My data
//-------------------------------------------------------------------
  TGraphAsymmErrors* gr_syst = new TGraphAsymmErrors(3);
  TGraphErrors* gr_globsyst = new TGraphErrors(3);
  TGraphErrors* gr_stat = new TGraphErrors(3);
  TGraphErrors* gr_upper_limit = new TGraphErrors(2);
  TObjArray arrows;

  Double_t nCoh, errCoh_stat, errCoh_syst_inf, errCoh_syst_sup, errCoh_glob, nUpperLimit;

  Int_t tminCentConsidered[5];
  Int_t tmaxCentConsidered[5];
  Double_t xRange[5];
  Double_t dxRange[5];

  std::copy(std::begin(minCentInverse), std::end(minCentInverse), std::begin(tminCentConsidered));
  std::copy(std::begin(maxCentInverse), std::end(maxCentInverse), std::begin(tmaxCentConsidered));
  std::copy(std::begin(nPart2), std::end(nPart2), std::begin(xRange));
  std::copy(std::begin(dnPart2), std::end(dnPart2), std::begin(dxRange));

  for (int j=0; j<5; j++)
  {

    nCoh = ResultsCRRun1Data[SetRangeValue(kRawYield, -4, -2.5, 0, 0.3, tminCentConsidered[j], tmaxCentConsidered[j])][kCrossSection]; 
    printf("Value %i is %.2f\n", j, nCoh);
    errCoh_stat = ResultsCRRun1Data[SetRangeValue(kRawYield, -4, -2.5, 0, 0.3, tminCentConsidered[j], tmaxCentConsidered[j])][kCrossSectionStat];
    errCoh_syst_inf = ResultsCRRun1Data[SetRangeValue(kRawYield, -4, -2.5, 0, 0.3, tminCentConsidered[j], tmaxCentConsidered[j])][3];
    errCoh_syst_sup = ResultsCRRun1Data[SetRangeValue(kRawYield, -4, -2.5, 0, 0.3, tminCentConsidered[j], tmaxCentConsidered[j])][2];
    errCoh_glob = ResultsCRRun1Data[SetRangeValue(kRawYield, -4, -2.5, 0, 0.3, tminCentConsidered[j], tmaxCentConsidered[j])][4];
    printf("Get %i value ok\n", j);
    
    if(j<3)
    {
      gr_syst->SetPoint(j, xRange[j], nCoh);
      gr_syst->SetPointError(j,dxRange[j],dxRange[j], errCoh_syst_inf, errCoh_syst_sup);

      gr_globsyst->SetPoint(j, xRange[j], nCoh);
      gr_globsyst->SetPointError(j,dxRange[j], errCoh_glob);

      gr_stat->SetPoint(j, xRange[j], nCoh);
      gr_stat->SetPointError(j,0, errCoh_stat); 
      printf("Fill histo %i ok\n", j);
    }
    else
    {
      gr_upper_limit->SetPoint(j-3, xRange[j], nCoh);
      gr_upper_limit->SetPointError(j-3, dxRange[j], 0);
      arrows.Add(new TArrow(xRange[j], nCoh, xRange[j], 0, 0.015, "|>"));
      dynamic_cast<TArrow*>(arrows.Last())->Draw("same|>");
    }
    

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

  gr_upper_limit->SetMarkerColor(colors[0]);
  gr_upper_limit->SetMarkerSize(1);
  gr_upper_limit->SetMarkerStyle(1);
  gr_upper_limit->SetLineColor(colors[0]);
  gr_upper_limit->SetFillColor(colors[0]);
  gr_upper_limit->SetFillStyle(0);
  gr_upper_limit->Draw("pezsame");
  
  if(!isGlobalSystematicUncluded) 
  {
    TLatex * text4 = new TLatex (25,291.5,"Cent. corr. syst. uncert. = ^{+13.7%}_{-14.3%}");
    text4->Draw("SAME");
    text4->SetTextSizePixels(20);
  }  

  if(isModels)
  {
    dlegend2->SetNColumns(2);
    dlegend2->AddEntry(gr_stat,"Data","pfe");
    dlegend2->AddEntry((TObject*)0, "", "");
    // dlegend2->AddEntry((TObject*)0, "M.B. Gay-Ducati et al., Phys. Rev. D97 (2018) 116013", "");
    // dlegend2->AddEntry((TObject*)0, "M.B. Gay-Ducati et al.", "");
    // dlegend2->AddEntry((TObject*)0, "", "");
    // dlegend2->AddEntry((TObject*)0, "Phys. Rev. D97 (2018) 116013", "");
    // dlegend2->AddEntry((TObject*)0, "", "");
    // dlegend2->SetNColumns(2);

    dlegend2->AddEntry(gr_IIM2,"IIM S2","l");
    dlegend2->AddEntry(gr_IIM3,"IIM S3","l");
    dlegend2->AddEntry(gr_GBW2,"GBW S2","l");
    dlegend2->AddEntry(gr_GBW3,"GBW S3","l");
    // dlegend2->SetNColumns(1);
    // dlegend2->AddEntry((TObject*)0,"","");
    // dlegend2->AddEntry((TObject*)0,"","");
    // dlegend2->AddEntry((TObject*)0, "J. Cepila et al., Phys. Rev. C97 (2018) 024901", "");
    // dlegend2->AddEntry((TObject*)0, "J. Cepila et al.", "");
    // dlegend2->AddEntry((TObject*)0,"","");
    // dlegend2->AddEntry((TObject*)0, "Phys. Rev. C97 (2018) 024901", "");
    // dlegend2->AddEntry((TObject*)0,"","");
    dlegend2->AddEntry(gr_GG,"GG-hs ","l");
    // dlegend2->AddEntry(gr_GS,"GS-hs","l");
    dlegend2->AddEntry((TObject*)0,"","");

    // dlegend2->AddEntry((TObject*)0,"","");
    // dlegend2->AddEntry((TObject*)0,"","");
    // dlegend2->AddEntry((TObject*)0, "M. Klusek-Gawenda et al.", "");
    // dlegend2->AddEntry((TObject*)0,"","");
    // dlegend2->AddEntry((TObject*)0, "Phys. Rev. C93 (2016) 044912", "");
    // dlegend2->AddEntry((TObject*)0,"","");
    // dlegend2->AddEntry(gr_Geom,"Geometrical ","l");
    dlegend2->AddEntry(gr_Glaub,"VDM","l");
  }
  leg->SetFillColor(0);
  leg->SetTextSize(gStyle->GetTextSize()*0.55);
  // leg->Draw();
  dlegend2->SetFillColor(0);
  dlegend2->SetTextSize(gStyle->GetTextSize()*0.7);
  dlegend2->Draw();
}