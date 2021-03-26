/*
This macro compare the pt and rapidity Acc$effi obtained for various iteration and centrality bins. It is not an intermediate step but rather used after the last iteration to produce a summary plot
*/

//Compare for each centrality the different input shapes of the different iteration

std::vector<std::vector<TString>> CompareEffi(int centMin = 0, int centMax = 10, int lastIteration = 3, TString rapOrPt = "Rap")
{

  TCanvas *canAccEffi = new TCanvas(Form("canIterationsAccEffiVs%s_Cent-%dto%d", rapOrPt.Data(), centMin, centMax), "", 900, 800);
  SetCanvasStyle(canAccEffi);

  TLegend *leg;
  if (rapOrPt == "Pt")
  {
    leg = new TLegend(0.13, 0.54, 0.38, 0.88);
  }
  else
  {
    leg = new TLegend(0.37, 0.18, 0.62, 0.49);
  }
  leg->SetMargin(0.1);
  leg->SetFillStyle(0);
  leg->SetLineColorAlpha(0, 0);
  leg->SetTextColor(kBlack);
  leg->AddEntry("NULL", Form("centrality %d-%d %%", centMin, centMax), "");

  TPad *padMain = padMain = new TPad("padMain", "padMain", 0, 0.3, 1, 1, 0);
  padMain->SetBottomMargin(0.);
  padMain->Draw();

  TPad *padRatio = new TPad("padRatio", "padRatio", 0, 0, 1, 0.3, 0);
  padRatio->SetBottomMargin(0.4);
  padRatio->SetTopMargin(0.);
  padRatio->Draw();

  TFile *inputAccEffFile_0 = TFile::Open(Form("Cent-%dto%d/AccEffi/iter-0/AccEffiValues.root", centMin, centMax));
  TH1F *histoAccEffi_0 = ((TH1F *)inputAccEffFile_0->Get(Form("histoAccEffiVs%s", rapOrPt.Data())));
  SetHistoStyle(histoAccEffi_0, arrayOfColors[0], kFullCircle, 0.8, rapOrPt == "Rap" ? yAxisTitle : ptAxisTitle, accEffiAxisTitle, 1, kFALSE);
  padMain->cd();
  histoAccEffi_0->Draw();
  leg->AddEntry(histoAccEffi_0, Form("Iteration: Step-0"), "p");
  leg->Draw();

  // assign a string for each bin (it will be used later for html/latex table drawing)
  std::vector<std::vector<TString>> accEffiVector;
  std::vector<TString> headerVector;
  headerVector.push_back(rapOrPt == "Rap" ? "y" : "pt (GeV/c)");
  headerVector.push_back("AccEffi-0");
  accEffiVector.push_back(headerVector);
  int numberOfBins = histoAccEffi_0->GetNbinsX();
  for (int iBin = 1; iBin <= numberOfBins; iBin++)
  {
    std::vector<TString> binVector;

    binVector.push_back(Form("%2.2f,%2.2f", histoAccEffi_0->GetXaxis()->GetBinLowEdge(iBin), histoAccEffi_0->GetXaxis()->GetBinUpEdge(iBin)));

    binVector.push_back(Form("%2.2f%s%4.4f", histoAccEffi_0->GetBinContent(iBin), plusOrMinusSymbol.Data(), histoAccEffi_0->GetBinError(iBin)));

    accEffiVector.push_back(binVector);
  }

  for (int iteration = 1; iteration <= lastIteration; iteration++)
  {
    TFile *inputAccEffFile = TFile::Open(Form("Cent-%dto%d/AccEffi/iter-%d/AccEffiValues.root", centMin, centMax, iteration));
    TH1F *histoAccEffi = ((TH1F *)inputAccEffFile->Get(Form("histoAccEffiVs%s", rapOrPt.Data())));
    SetHistoStyle(histoAccEffi, arrayOfColors[iteration], kFullCircle, 0.8, rapOrPt == "Rap" ? yAxisTitle : ptAxisTitle, accEffiAxisTitle, 1, kFALSE);

    padMain->cd();
    histoAccEffi->Draw("same");
    leg->AddEntry(histoAccEffi, Form("Iteration: Step-%d", iteration), "p");

    // ratio:
    padRatio->cd();
    TH1F *histoRatio = new TH1F(*histoAccEffi);
    histoRatio->Divide(histoAccEffi_0);
    SetHistoStyle(histoRatio, arrayOfColors[iteration], kFullCircle, 0.8, rapOrPt == "Rap" ? yAxisTitle : ptAxisTitle, "ratio", 1, kTRUE);
    std::vector<double> userRange = GetAVisuallyGoodYRangeUser(histoRatio->GetMinimum() - histoRatio->GetBinError(histoRatio->GetNbinsX()), histoRatio->GetMaximum() + histoRatio->GetBinError(histoRatio->GetNbinsX()), 1);
    histoRatio->GetYaxis()->SetRangeUser(userRange[0], userRange[1]);
    histoRatio->GetYaxis()->SetNdivisions(505);

    histoRatio->Draw(iteration == 1 ? "" : "same");
    if (iteration == 1)
    {
      TLine *lineAtOne = new TLine(histoRatio->GetXaxis()->GetXmin(), 1, histoRatio->GetXaxis()->GetXmax(), 1);
      lineAtOne->SetLineStyle(kDashed);
      lineAtOne->Draw("same");
    }

    // push to header:
    accEffiVector[0].push_back(Form("AccEffi-%d", iteration));
    accEffiVector[0].push_back(Form("|(%d-0)/0| (%%)", iteration));
    for (int iBin = 1; iBin <= numberOfBins; iBin++)
    {
      accEffiVector[iBin].push_back(Form("%2.2f%s%4.4f", histoAccEffi->GetBinContent(iBin), plusOrMinusSymbol.Data(), histoAccEffi->GetBinError(iBin)));
      accEffiVector[iBin].push_back(Form("%2.2f", 100 * TMath::Abs(histoRatio->GetBinContent(iBin) - 1)));
    }
  }

  canAccEffi->SaveAs(Form("Cent-%dto%d/FinalPlots/%sAccEffiVsIter.png", centMin, centMax, rapOrPt.Data()));
  //gDirectory->GetList()->Delete();
  return accEffiVector;
}

void CompareIterationsAccEffi(int centMin, int centMax, int lastIteration, std::vector<std::vector<TString>> &accEffiVsRap, std::vector<std::vector<TString>> &accEffiVsPt)
{
  if (lastIteration == -1)
  {
    //Get the last iteration for the first multiplicity bin (for now assume it is the same for all the bins):
    TString strLastIter = GetLastIteration(centMin, centMax, "PtShapeIterations");
    lastIteration = ((TString)strLastIter(strLastIter.Length() - 1, strLastIter.Length() - 1)).Atoi();
  }

  accEffiVsRap = CompareEffi(centMin, centMax, lastIteration, "Rap");
  accEffiVsPt = CompareEffi(centMin, centMax, lastIteration, "Pt");
}
