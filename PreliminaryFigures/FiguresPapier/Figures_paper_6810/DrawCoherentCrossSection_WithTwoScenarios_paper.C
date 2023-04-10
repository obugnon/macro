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
int minCentInverse[5]={70,50,30, 10, 0};
int maxCentInverse[5]={90,70,50, 30, 10};
// Double_t nPart5[3]={11.35,42.66,108.97};
// Double_t dnPart5[5]={1.5, 1.5, 1.5,1.5, 1.5};
Double_t dnPart5[5]={5.5, 5.5, 5.5, 5.5, 5.5};


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
  TH1 * h = cfig->DrawFrame(0,0,400,700);

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

  TLatex * text = new TLatex (25,625,"ALICE, Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV");
  text->Draw("SAME");
  text->SetTextSizePixels(24);

  TLatex * text2 = new TLatex (25,546,"J/#psi #rightarrow #mu^{+}#mu^{-}, 2.5 < #it{y} < 4");
  text2->SetTextSizePixels(22);
  text2->Draw("SAME");
  TLatex * text3 = new TLatex (25,492.6,"#it{p}_{T} < 0.3 GeV/#it{c}");
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
      nCoh_IIM2 = Models5TeV[Form("Ducati_IIM2_%d_%d", minCentInverse[j], maxCentInverse[j])]; 
      nCoh_GBW2 = Models5TeV[Form("Ducati_GBW2_%d_%d", minCentInverse[j], maxCentInverse[j])];

      nCoh_IIM3 = Models5TeV[Form("Ducati_IIM3_%d_%d", minCentInverse[j], maxCentInverse[j])]; 
      nCoh_GBW3 = Models5TeV[Form("Ducati_GBW3_%d_%d", minCentInverse[j], maxCentInverse[j])];

      gr_IIM2->SetPoint(j, nPart5[j], nCoh_IIM2);
      gr_IIM2->SetPointError(j,dnPart5[j]*2, 0);

      gr_GBW2->SetPoint(j, nPart5[j], nCoh_GBW2);
      gr_GBW2->SetPointError(j,dnPart5[j]*2, 0);

      gr_IIM3->SetPoint(j, nPart5[j], nCoh_IIM3);
      gr_IIM3->SetPointError(j,dnPart5[j]*2, 0);

      gr_GBW3->SetPoint(j, nPart5[j], nCoh_GBW3);
      gr_GBW3->SetPointError(j,dnPart5[j]*2, 0);
      
      nCoh_Glaub = Models5TeV[Form("Klusek_Glauber_%d_%d", minCentInverse[j], maxCentInverse[j])]; 
      nCoh_Geom = Models5TeV[Form("Klusek_Geometrical_%d_%d", minCentInverse[j], maxCentInverse[j])];

      gr_Glaub->SetPoint(j, nPart5[j], nCoh_Glaub);
      gr_Glaub->SetPointError(j,dnPart5[j]*2, 0);

      gr_Geom->SetPoint(j, nPart5[j], nCoh_Geom);
      gr_Geom->SetPointError(j,dnPart5[j]*2, 0);

      if(j==0)
      {
        nCoh_GS = Models5TeV[Form("Guillermo_GS_%d_%d",  minCentInverse[j], maxCentInverse[j])];
        nCoh_GG = Models5TeV[Form("Guillermo_GG_%d_%d",  minCentInverse[j], maxCentInverse[j])];
      
        gr_GG->SetPoint(j, nPart5[j], nCoh_GG);
        gr_GG->SetPointError(j,dnPart5[j]*2, 0);

        gr_GS->SetPoint(j, nPart5[j], nCoh_GS);
        gr_GS->SetPointError(j,dnPart5[j]*2, 0);
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
    // gr_Geom->SetLineStyle(1);
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
  TGraphErrors* gr_syst = new TGraphErrors(4);
  TGraphErrors* gr_globsyst = new TGraphErrors(4);
  TGraphErrors* gr_stat = new TGraphErrors(4);
  TGraphErrors* gr_upper_limit = new TGraphErrors(1);
  TObjArray arrows;

  Double_t nCoh, errCoh_stat, errCoh_syst, errCoh_glob, nUpperLimit;

  Int_t tminCentConsidered[5];
  Int_t tmaxCentConsidered[5];
  Double_t xRange[5];
  Double_t dxRange[5];

  std::copy(std::begin(minCentInverse), std::end(minCentInverse), std::begin(tminCentConsidered));
  std::copy(std::begin(maxCentInverse), std::end(maxCentInverse), std::begin(tmaxCentConsidered));
  std::copy(std::begin(nPart5), std::end(nPart5), std::begin(xRange));
  std::copy(std::begin(dnPart5), std::end(dnPart5), std::begin(dxRange));

  for (int j=0; j<5; j++)
  {

    nCoh = ResultsCRThisAnalysis[SetRangeValue(kRawYield, -4, -2.5, 0, 0.3, tminCentConsidered[j], tmaxCentConsidered[j])][kCrossSection]; 
    errCoh_stat = ResultsCRThisAnalysis[SetRangeValue(kRawYield, -4, -2.5, 0, 0.3, tminCentConsidered[j], tmaxCentConsidered[j])][kCrossSectionStat];
    errCoh_syst = ResultsCRThisAnalysis[SetRangeValue(kRawYield, -4, -2.5, 0, 0.3, tminCentConsidered[j], tmaxCentConsidered[j])][kCrossSectionSyst];
    errCoh_glob = ResultsCRThisAnalysis[SetRangeValue(kRawYield, -4, -2.5, 0, 0.3, tminCentConsidered[j], tmaxCentConsidered[j])][kCrossSectionSystCor];
    
    if(j<4)
    {
      gr_syst->SetPoint(j, xRange[j], nCoh);
      gr_syst->SetPointError(j,dxRange[j], errCoh_syst);

      gr_globsyst->SetPoint(j, xRange[j], nCoh);
      gr_globsyst->SetPointError(j,dxRange[j], errCoh_glob);

      gr_stat->SetPoint(j, xRange[j], nCoh);
      gr_stat->SetPointError(j,0, errCoh_stat); 
    }
    else
    {
      nUpperLimit=TMath::Sqrt((errCoh_stat*errCoh_stat)+(errCoh_syst*errCoh_syst))*2+nCoh;
      printf("Upper Limit is %.2f\n", nUpperLimit);
      gr_upper_limit->SetPoint(0, xRange[j], nUpperLimit);
      gr_upper_limit->SetPointError(0, dxRange[j], 0);
      arrows.Add(new TArrow(xRange[j], nUpperLimit, xRange[j], 0, 0.015, "|>"));

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
  
  dynamic_cast<TArrow*>(arrows.Last())->Draw();


  if(!isGlobalSystematicUncluded) 
  {
    TLatex * text4 = new TLatex (25,441.5,Form("Cent. corr. syst. uncert. = %.1f%%", errGlob_ThisAnalysis));
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