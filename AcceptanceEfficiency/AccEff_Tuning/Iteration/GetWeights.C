/*
This macro return the weight as a function of y for a given iteration step. The weight is given by the ratio of the normalised y fit from the current step and the one before (for the step-0 take the generated input shape).
The macro is basically the one sent by Laure (weight.C)
*/
#define UNUSED(x) (void)(x)

//--------------------------------------------------------------------------------------------//
//Define the two rapidity distributions and the shape descrirapions
TF1 *rapDistributionOld = NULL;
TF1 *rapDistributionNew = NULL;
TF1 *rapDistributionNewScaled = NULL;
TF1 *rapWeight = NULL;
Double_t rapIntegralOld = 0;
Double_t rapIntegralNew = 0;

Double_t RescaleRap(Double_t *x, Double_t *par)
{
  UNUSED(par);
  const Double_t xx = x[0];
  return (rapDistributionNew->Eval(x[0]) * rapIntegralOld / rapIntegralNew);
}

Double_t rapdist(Double_t *x, Double_t *par)
{
  return (par[0] * TMath::Exp(-0.5 * TMath::Power(x[0] / par[1], 2)));
}

Double_t DivideRap(Double_t *x, Double_t *par)
{
  UNUSED(par);
  return (rapDistributionNewScaled->Eval(x[0])) / (rapDistributionOld->Eval(x[0]));
}
//--------------------------------------------------------------------------------------------//

//--------------------------------------------------------------------------------------------//
//Simle funtion to get a TF1 from a given step
TF1 *GetRapDistribution(Int_t centMin = 0, Int_t centMax = 10, Int_t iteration = 1, Double_t rapMin = -4, Double_t rapMax = -2.5)
{
  TString inputData = gSystem->GetFromPipe(Form("cat Cent-%dto%d/RapShapeIterations/iter-%d/values.txt", centMin, centMax, iteration));
  TObjArray *objInputData = inputData.Tokenize("\n");

  TF1 *rapDistribution = new TF1(Form("RapDistribution_%d", iteration), rapdist, rapMin, rapMax, 2);
  rapDistribution->FixParameter(0, ((TObjString *)objInputData->UncheckedAt(0))->String().Atof());
  rapDistribution->FixParameter(1, ((TObjString *)objInputData->UncheckedAt(2))->String().Atof());
  return rapDistribution;
}
//--------------------------------------------------------------------------------------------//

//--------------------------------------------------------------------------------------------//
//The main function. It return a histogram with the weight values (TODO: return the weight as TF1 without converting it to TH1. This has shown some problems when used for weighting)
TH1 *GetRapWeight(Int_t centMin = 0, Int_t centMax = 10, Int_t iteration = 1, Bool_t drawIt = kFALSE)
{
  PrintMessage("Getting y weight", "-");

  if (iteration < 1)
  {
    cout << "This macro needs at least one previous iteration." << endl;
    return NULL;
  }

  int numberOfRapidityBins;
  const Double_t *arrayRapidityBins = GetBinsArray(centMin, centMax, "Rap", numberOfRapidityBins);

  rapDistributionOld = GetRapDistribution(centMin, centMax, iteration - 1, arrayRapidityBins[0], arrayRapidityBins[numberOfRapidityBins]);
  rapIntegralOld = rapDistributionOld->Integral(arrayRapidityBins[0], arrayRapidityBins[numberOfRapidityBins]);

  rapDistributionNew = GetRapDistribution(centMin, centMax, iteration, arrayRapidityBins[0], arrayRapidityBins[numberOfRapidityBins]);
  rapIntegralNew = rapDistributionNew->Integral(arrayRapidityBins[0], arrayRapidityBins[numberOfRapidityBins]);

  rapDistributionNewScaled = new TF1("rapgenrescaled", RescaleRap, arrayRapidityBins[0], arrayRapidityBins[numberOfRapidityBins], 0);
  rapWeight = new TF1("rapDivide", DivideRap, arrayRapidityBins[0], arrayRapidityBins[numberOfRapidityBins], 0);

  TH1 *histo = rapWeight->CreateHistogram();

  if (drawIt)
  {
    TCanvas *canWeightRap = new TCanvas(Form("canWeightRap_%d_Cent-%dto%d", iteration, centMin, centMax), "", 900, 800);
    SetCanvasStyle(canWeightRap);

    TPad *padMain, *padRatio;
    padMain = new TPad("padMain", "padMain", 0, 0.3, 1, 1, 0);
    padMain->SetBottomMargin(0.);
    padMain->Draw();

    padRatio = new TPad("padRatio", "padRatio", 0, 0, 1, 0.3, 0);
    padRatio->SetBottomMargin(0.4);
    padRatio->SetTopMargin(0.);
    padRatio->Draw();

    padMain->cd();
    // gPad->SetLogy();
    rapDistributionNewScaled->SetLineColor(myColorsMap["red"]);
    rapDistributionNewScaled->Draw();
    rapDistributionOld->SetLineColor(myColorsMap["blue"]);
    rapDistributionOld->Draw("SAME");

    TLegend *legMain = new TLegend(0.37, 0.59, 0.59, 0.83);
    legMain->SetMargin(0.1);
    legMain->SetFillStyle(0);
    legMain->SetLineColorAlpha(0, 0);
    legMain->SetTextColor(kBlack);
    legMain->AddEntry("NULL", Form("Centrality %d-%d %%", centMin, centMax), "");
    legMain->AddEntry(rapDistributionOld, Form("Iteration: Step-%d", iteration - 1), "l");
    legMain->AddEntry(rapDistributionNewScaled, Form("Iteration: Step-%d", iteration), "l");
    legMain->Draw();

    padRatio->cd();
    TH1F *histoEmpty = new TH1F("empty_rap", "", 10, arrayRapidityBins[0], arrayRapidityBins[numberOfRapidityBins]);
    SetHistoStyle(histoEmpty, kBlack, kFullCircle, 0.8, yAxisTitle, "", 1, kTRUE);
    std::vector<double> rangeUser = GetAVisuallyGoodYRangeUser(histo->GetMinimum(), histo->GetMaximum(), 0.5);
    histoEmpty->GetYaxis()->SetRangeUser(rangeUser[0], rangeUser[1]);
    histoEmpty->GetYaxis()->SetNdivisions(505);
    histoEmpty->Draw("");
    rapWeight->SetLineColor(myColorsMap["red"]);
    rapWeight->Draw("same");
    //Draw line at one:
    TLine *lineAtOne = new TLine(histoEmpty->GetXaxis()->GetXmin(), 1, histoEmpty->GetXaxis()->GetXmax(), 1);
    lineAtOne->SetLineStyle(kDashed);
    lineAtOne->Draw();
    canWeightRap->SaveAs(Form("Cent-%dto%d/RapShapeIterations/iter-%d/RapWeight.png", centMin, centMax, iteration));
  }

  return histo;
}
//--------------------------------------------------------------------------------------------//
/*
This macro return the weight as a function of pT for a given iteration step. The weight is given by the ratio of the normalised Pt fit from the current step and the one before (for the step-0 take the generated input shape).
The macro is basically the one sent by Laure (weight.C)
*/

//--------------------------------------------------------------------------------------------//
//Define the two pt distributions and the shape descriptions
TF1 *ptDistributionOld = NULL;
TF1 *ptDistributionNew = NULL;
TF1 *ptDistributionNewScaled = NULL;
Double_t ptIntegralOld = 0;
Double_t ptIntegralNew = 0;

Double_t RescalePt(Double_t *x, Double_t *par)
{
  UNUSED(par);
  const Double_t xx = x[0];
  return (ptDistributionNew->Eval(x[0]) * ptIntegralOld / ptIntegralNew);
}

Double_t DividePt(Double_t *x, Double_t *par)
{
  UNUSED(par);
  return (ptDistributionNewScaled->Eval(x[0])) / (ptDistributionOld->Eval(x[0]));
}
//--------------------------------------------------------------------------------------------//

//--------------------------------------------------------------------------------------------//
//Simle funtion to get a TF1 from a given step
TF1 *GetPtDistribution(Int_t centMin = 0, Int_t centMax = 10, Int_t iteration = 1, Double_t ptMin = 0.3, Double_t ptMax = 15)
{
  TString inputData = gSystem->GetFromPipe(Form("cat Cent-%dto%d/PtShapeIterations/iter-%d/values.txt", centMin, centMax, iteration));
  TObjArray *objInputData = inputData.Tokenize("\n");

  TF1 *ptDistribution = new TF1(Form("PtDistribution_%d", iteration), pT_shape, ptMin, ptMax, 4);
  for (Int_t iParam = 0; iParam < 4; iParam++)
  {
    ptDistribution->FixParameter(iParam, ((TObjString *)objInputData->UncheckedAt(iParam))->String().Atof());
  }

  return ptDistribution;
}
//--------------------------------------------------------------------------------------------//

//--------------------------------------------------------------------------------------------//
//The main function. It return a histogram with the weight values (TODO: return the weight as TF1 without converting it to TH1. This has shown some problems when used for weighting)
TH1 *GetPtWeight(Int_t centMin = 0, Int_t centMax = 10, Int_t iteration = 1, Bool_t drawIt = kFALSE)
{
  PrintMessage("Getting Pt weight", "-");

  int numberOfPtBins;
  const Double_t *arrayPtBins = GetBinsArray(centMin, centMax, "Pt", numberOfPtBins);

  ptDistributionOld = GetPtDistribution(centMin, centMax, iteration - 1, arrayPtBins[0], arrayPtBins[numberOfPtBins]);

  ptIntegralOld = ptDistributionOld->Integral(arrayPtBins[0], arrayPtBins[numberOfPtBins]);

  ptDistributionNew = GetPtDistribution(centMin, centMax, iteration, arrayPtBins[0], arrayPtBins[numberOfPtBins]);
  ptIntegralNew = ptDistributionNew->Integral(arrayPtBins[0], arrayPtBins[numberOfPtBins]);

  ptDistributionNewScaled = new TF1("pTgenrescaled", RescalePt, arrayPtBins[0], arrayPtBins[numberOfPtBins], 0);
  TF1 *ptWeight = new TF1(Form("ptWeight_%d", iteration), DividePt, arrayPtBins[0], arrayPtBins[numberOfPtBins], 0);
  TH1 *histo = ptWeight->CreateHistogram();

  if (drawIt)
  {
    TCanvas *canWeightPt = new TCanvas(Form("canWeightPt_%d_Cent-%dto%d", iteration, centMin, centMax), "", 900, 800);
    SetCanvasStyle(canWeightPt);
    canWeightPt->SetLogy();

    TPad *padMain, *padRatio;
    padMain = new TPad("padMain", "padMain", 0, 0.3, 1, 1, 0);
    padMain->SetBottomMargin(0.);
    padMain->Draw();

    padRatio = new TPad("padRatio", "padRatio", 0, 0, 1, 0.3, 0);
    padRatio->SetBottomMargin(0.4);
    padRatio->SetTopMargin(0.);
    padRatio->Draw();

    padMain->cd();
    gPad->SetLogy();
    ptDistributionNewScaled->SetLineColor(myColorsMap["red"]);
    ptDistributionNewScaled->Draw();
    ptDistributionOld->SetLineColor(myColorsMap["blue"]);
    ptDistributionOld->Draw("SAME");

    TLegend *legMain = new TLegend(0.57, 0.58, 0.82, 0.83);
    legMain->SetMargin(0.1);
    legMain->SetFillStyle(0);
    legMain->SetLineColorAlpha(0, 0);
    legMain->SetTextColor(kBlack);
    legMain->AddEntry("NULL", Form("Centrality %d-%d %%", centMin, centMax), "");
    legMain->AddEntry(ptDistributionOld, Form("Iteration: Step-%d", iteration - 1), "l");
    legMain->AddEntry(ptDistributionNewScaled, Form("Iteration: Step-%d", iteration), "l");
    legMain->Draw();

    padRatio->cd();
    TH1F *histoEmpty = new TH1F("empty", "", 10, arrayPtBins[0], arrayPtBins[numberOfPtBins]);
    SetHistoStyle(histoEmpty, kBlack, kFullCircle, 0.8, ptAxisTitle, "", 1, kTRUE);
    std::vector<double> rangeUser = GetAVisuallyGoodYRangeUser(histo->GetMinimum(), histo->GetMaximum(), 0.5);
    histoEmpty->GetYaxis()->SetRangeUser(rangeUser[0], rangeUser[1]);
    histoEmpty->GetYaxis()->SetNdivisions(505);
    histoEmpty->Draw("");
    ptWeight->SetLineColor(myColorsMap["red"]);
    ptWeight->Draw("same");
    //Draw line at one:
    TLine *lineAtOne = new TLine(histoEmpty->GetXaxis()->GetXmin(), 1, histoEmpty->GetXaxis()->GetXmax(), 1);
    lineAtOne->SetLineStyle(kDashed);
    lineAtOne->Draw();
    canWeightPt->SaveAs(Form("Cent-%dto%d/PtShapeIterations/iter-%d/PtWeight.png", centMin, centMax, iteration));
  }

  return histo;
}
//--------------------------------------------------------------------------------------------//