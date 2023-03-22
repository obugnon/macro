/*
This macro compare the pt and rapidity shapes obtained for various iteration and centrality bins. It is not an intermediate step but rather used after the last iteration to produce a summary plot
*/

//--------------------------------------------------------------------------------------------//
//Compare for each centrality the different input shapes of the different iteration
void CompareIterationsShape(int centMin = 0, int centMax = 10, int lastIteration = 3)
{

  if (lastIteration == -1)
  {
    //Get the last iteration for the first multiplicity bin (for now assume it is the same for all the bins):
    TString strLastIter = GetLastIteration(centMin, centMax, "PtShapeIterations");
    lastIteration = ((TString)strLastIter(strLastIter.Length() - 1, strLastIter.Length() - 1)).Atoi();
  }

  int numberOfRapidityBins;
  const Double_t *arrayRapidityBins = GetBinsArray(centMin, centMax, "Rap", numberOfRapidityBins);

  int numberOfPtBins;
  const Double_t *arrayPtBins = GetBinsArray(centMin, centMax, "Pt", numberOfPtBins);

  TCanvas *canRapShapes = new TCanvas(Form("canRapShapes_Cent%dto%d", centMin, centMax), "", 900, 800);
  SetCanvasStyle(canRapShapes);
  TH1F *histoEmptyRap = new TH1F(Form("histoEmptyRap_Cent%dto%d", centMin, centMax), "", 10, arrayRapidityBins[0], arrayRapidityBins[numberOfRapidityBins]);
  SetHistoStyle(histoEmptyRap, kWhite, kFullCircle, 0.8, "#it{y}", "dN/d#it{y}", 1, kFALSE);
  histoEmptyRap->GetXaxis()->SetLabelSize(0.045);
  histoEmptyRap->GetXaxis()->SetTitleSize(0.045);
  histoEmptyRap->GetYaxis()->SetRangeUser(0.35, 1);
  histoEmptyRap->Draw();

  TLegend *legRap = new TLegend(0.639, 0.185, 0.889, 0.536);
  legRap->SetFillStyle(0);
  legRap->SetLineColorAlpha(0, 0);
  legRap->SetTextColor(kBlack);
  legRap->SetMargin(0.1);
  legRap->AddEntry("NULL", Form("centrality: %d-%d %%", centMin, centMax), "");

  TCanvas *canPtShapes = new TCanvas(Form("canPtShapes_Cent%dto%d", centMin, centMax), "", 900, 800);
  SetCanvasStyle(canPtShapes);
  canPtShapes->SetLogy();
  TH1F *histoEmptyPt = new TH1F(Form("histoEmptyPt_Cent%dto%d", centMin, centMax), "", 10, arrayPtBins[0], arrayPtBins[numberOfPtBins]);
  SetHistoStyle(histoEmptyPt, kWhite, kFullCircle, 0.8, ptAxisTitle, "dN/d#it{p}_{T}", 1, kFALSE);
  histoEmptyPt->GetXaxis()->SetLabelSize(0.045);
  histoEmptyPt->GetXaxis()->SetTitleSize(0.045);
  histoEmptyPt->GetYaxis()->SetRangeUser(0.00001, 0.5);
  histoEmptyPt->Draw();

  TLegend *legPt = new TLegend(0.25, 0.185, 0.5, 0.536);
  legPt->SetFillStyle(0);
  legPt->SetLineColorAlpha(0, 0);
  legPt->SetTextColor(kBlack);
  legPt->SetMargin(0.1);
  legPt->AddEntry("NULL", Form("centrality: %d-%d %%", centMin, centMax), "");

  //For each iteration get the rapidity and pt shapes
  for (int iIteration = 0; iIteration <= lastIteration; iIteration++)
  {
    canRapShapes->cd();
    TF1 *fRap = GetRapDistribution(centMin, centMax, iIteration);
    fRap->SetLineColor((arrayOfColors[iIteration]));
    fRap->SetParameter(0, fRap->GetParameter(0) / fRap->Integral(arrayRapidityBins[0], arrayRapidityBins[numberOfRapidityBins]));
    fRap->Draw("same");
    legRap->AddEntry(fRap, Form("iteration %d", iIteration), "l");

    canPtShapes->cd();
    TF1 *fPt = GetPtDistribution(centMin, centMax, iIteration, arrayPtBins[0], arrayPtBins[numberOfPtBins]);
    fPt->SetLineColor((arrayOfColors[iIteration]));
    fPt->SetParameter(0, fPt->GetParameter(0) / fPt->Integral(arrayPtBins[0], arrayPtBins[numberOfPtBins]));
    fPt->Draw("same");
    legPt->AddEntry(fPt, Form("iteration %d", iIteration), "l");
  }
  canRapShapes->cd();
  legRap->Draw();
  canRapShapes->SaveAs(Form("Cent-%dto%d/FinalPlots/RapShapesVsIter.png", centMin, centMax));

  canPtShapes->cd();
  legPt->Draw();
  canPtShapes->SaveAs(Form("Cent-%dto%d/FinalPlots/PtShapesVsIter.png", centMin, centMax));
  //gDirectory->GetList()->Delete();
}
//--------------------------------------------------------------------------------------------//
