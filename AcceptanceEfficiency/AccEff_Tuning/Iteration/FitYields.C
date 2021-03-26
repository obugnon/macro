/*
For a given iteration step, this macro calculate the yield of J/psi vs pt using the Acc*Eff values obtained from the previous step. It then fits the distribution and store the fit parameters to be used in the next step. The fit function is a power law plus an exponential for high pt. For the zeroth iteration the fit is performed on the generated jpsi in a given centrality which must result in parameters similar to the one used for the MC generation (with no centrality dependence: THIS HAS TO BE CHECKED)
*/

//--------------------------------------------------------------------------------------------//
// TFile *inputRawYieldFile = TFile::Open("./RawYields/NJpsi_Raw.root");

//--------------------------------------------------------------------------------------------//

int FitPtYield(TH1 *histoYiledVsPt, bool draw = false)
{
  TF1 *fPtDistribution = new TF1("fPtDistribution", pT_shape, histoYiledVsPt->GetXaxis()->GetXmin(), histoYiledVsPt->GetXaxis()->GetXmax(), 4);
  fPtDistribution->SetLineColor(myColorsMap["red"]);
  fPtDistribution->SetParameter(0, histoYiledVsPt->GetMaximum());
  fPtDistribution->SetParameter(1, 4);
  fPtDistribution->SetParameter(2, 2);
  fPtDistribution->SetParameter(3, 3);

  TString fitOpt = "SRIQ";
  if (!draw)
  {
    fitOpt.Append("0");
  }

  double initParam = 3;
  double paramStep = -0.01;
  int customStepCounter = 0;
  TFitResultPtr fitResultPtr;
  double chi2OverNDF;
  do
  {
    fPtDistribution->SetParameter(1, initParam + paramStep * customStepCounter);
    fitResultPtr = histoYiledVsPt->Fit(fPtDistribution, fitOpt, "", ptFitRangeMin, ptFitRangeMax);
    chi2OverNDF = fPtDistribution->GetChisquare() / fPtDistribution->GetNDF();
    customStepCounter++;
  } while (fitResultPtr->CovMatrixStatus() != 3 || !fitResultPtr->IsValid() || customStepCounter > 100);
  return customStepCounter;
}

//--------------------------------------------------------------------------------------------//
void FitPtYield(int centMin = 0, int centMax = 10, Int_t iteration = 1, bool fitOnly = false)
{
  PrintMessage("Fitting pt yield", "-");
  TCanvas *canYieldPt = nullptr;
  if (!fitOnly)
  {
    canYieldPt = new TCanvas(Form("canYieldPt_%d_Cent-%dto%d", iteration, centMin, centMax), "", 900, 800);
    SetCanvasStyle(canYieldPt);
    canYieldPt->SetLogy();
  }

  //Get the corresponding Acc*Effi file
  TFile *inputAccEffFile = TFile::Open(Form("Cent-%dto%d/AccEffi/iter-%d/AccEffiValues.root", centMin, centMax, iteration - 1));
  TH1F *histoAccEffiVsPt = ((TH1F *)inputAccEffFile->Get("histoAccEffiVsPt"));

  TH1F *histoJpsiNumbersVsPt = ((TH1F *)inputFileRawYield->Get(Form("histoPtYield_cent%dto%d", centMin, centMax)));
  TH1F *histoYiledVsPt = new TH1F(*histoJpsiNumbersVsPt);
  histoYiledVsPt->Divide(histoAccEffiVsPt);
  SetHistoStyle(histoYiledVsPt, kBlack, kFullCircle, 0.9, ptAxisTitle, "d#it{N}/d#it{p}_{T}", 1, kFALSE);
  int customStepCounter = FitPtYield(histoYiledVsPt, true);
  TF1 *fPtDistribution = static_cast<TF1 *>(histoYiledVsPt->GetListOfFunctions()->FindObject("fPtDistribution"));
  histoYiledVsPt->GetYaxis()->SetTitleOffset(1.4);
  if (fitOnly)
  {
    return;
  }
  TLegend *legMain = new TLegend(0.54, 0.65, 0.86, 0.87);
  legMain->SetMargin(0.1);
  legMain->SetFillStyle(0);
  legMain->SetLineColorAlpha(0, 0);
  legMain->SetTextColor(customStepCounter < 100 ? kBlack : kRed);
  legMain->AddEntry("NULL", Form("Centrality: %d-%d %%", centMin, centMax), "");
  legMain->AddEntry("NULL", Form("Iteration: Step-%d", iteration), "");
  legMain->AddEntry(fPtDistribution, Form("C#times#frac{#it{p}_{T}}{  (1+ (#it{p}_{T}/%2.2f)^{%2.2f})^{%2.2f}  }", fPtDistribution->GetParameter(1), fPtDistribution->GetParameter(2), fPtDistribution->GetParameter(3)), "l");
  legMain->Draw();
  //--------------------------------------------------------------------------------------------//

  //--------------------------------------------------------------------------------------------//
  //Store the values of the fit. Use 1 (arbitrary) for the normalisation parameter (1st) as a normalisation will be done in later steps.
  gSystem->Exec(Form("mkdir -p Cent-%dto%d/PtShapeIterations/iter-%d", centMin, centMax, iteration));
  canYieldPt->SaveAs(Form("Cent-%dto%d/PtShapeIterations/iter-%d/PtYield.png", centMin, centMax, iteration));
  ofstream outputFile(Form("Cent-%dto%d/PtShapeIterations/iter-%d/values.txt", centMin, centMax, iteration), std::ofstream::out);
  outputFile << "1" << endl;
  outputFile << fPtDistribution->GetParameter(1) << endl;
  outputFile << fPtDistribution->GetParameter(2) << endl;
  outputFile << fPtDistribution->GetParameter(3) << endl;
  outputFile.close();
  //gDirectory->GetList()->Delete();
  //--------------------------------------------------------------------------------------------//
}
//--------------------------------------------------------------------------------------------//

/*
For a given iteration step, this macro calculate the yield of j/psi vs rapidity using the Acc*Eff values obtained from the previous step. It then fits the distribution and store the fit parameters to be used in the next step. The fit function is a 0-centered gaussian.
*/

int FitRapYield(TH1 *histoYiledVsRapidity, bool draw = false)
{
  TF1 *fRapDistribution = new TF1("fRapDistribution", "[0]*exp(-0.5*(x/[1])**2)", histoYiledVsRapidity->GetXaxis()->GetXmin(), histoYiledVsRapidity->GetXaxis()->GetXmax());
  fRapDistribution->SetParameter(0, histoYiledVsRapidity->GetMaximum());

  TString fitOpt = "SRIQ";
  if (!draw)
  {
    fitOpt.Append("0");
  }

  fRapDistribution->SetLineColor(myColorsMap["red"]);
  double initWidth = 3;
  double widthStep = -0.01;
  int customStepCounter = 0;
  TFitResultPtr fitResultPtr;
  double chi2OverNDF;
  do
  {
    fRapDistribution->SetParameter(1, initWidth + widthStep * customStepCounter);
    fitResultPtr = histoYiledVsRapidity->Fit(fRapDistribution, fitOpt.Data(), "", rapFitRangeMin, rapFitRangeMax);
    chi2OverNDF = fRapDistribution->GetChisquare() / fRapDistribution->GetNDF();
  } while (fitResultPtr->CovMatrixStatus() != 3 || !fitResultPtr->IsValid() || customStepCounter > 100);
  return customStepCounter;
}

//--------------------------------------------------------------------------------------------//
void FitRapYield(int centMin = 0, int centMax = 10, Int_t iteration = 1, bool fitOnly = false)
{
  PrintMessage("Fitting y yield", "-");

  TCanvas *canYieldRap = nullptr;
  if (!fitOnly)
  {
    canYieldRap = new TCanvas(Form("canYieldRap_%d_Cent-%dto%d", iteration, centMin, centMax), "", 900, 800);
    SetCanvasStyle(canYieldRap);
  }

  //Get the corresponding Acc*Effi file
  TFile *inputAccEffFile = TFile::Open(Form("Cent-%dto%d/AccEffi/iter-%d/AccEffiValues.root", centMin, centMax, iteration - 1));
  TH1F *histoAccEffiVsRap = ((TH1F *)inputAccEffFile->Get("histoAccEffiVsRap"));

  TH1F *histoJpsiNumbersVsRap = ((TH1F *)inputFileRawYield->Get(Form("histoRapYield_cent%dto%d", centMin, centMax)));
  TH1F *histoYiledVsRapidity = new TH1F(*histoJpsiNumbersVsRap);
  histoYiledVsRapidity->Divide(histoAccEffiVsRap);
  SetHistoStyle(histoYiledVsRapidity, kBlack, kFullCircle, 0.9, yAxisTitle, "d#it{N}/d#it{y}", 1, kFALSE);
  int customStepCounter = FitRapYield(histoYiledVsRapidity, true);
  TF1 *fRapDistribution = static_cast<TF1 *>(histoYiledVsRapidity->GetListOfFunctions()->FindObject("fRapDistribution"));
  if (fitOnly)
  {
    return;
  }
  TLegend *legMain = new TLegend(0.20, 0.61, 0.50, 0.87);
  legMain->SetMargin(0.1);
  legMain->SetFillStyle(0);
  legMain->SetLineColorAlpha(0, 0);
  legMain->SetTextColor(customStepCounter < 100 ? kBlack : kRed);
  legMain->AddEntry("NULL", Form("Centrality: %d-%d %%", centMin, centMax), "");
  legMain->AddEntry("NULL", Form("Iteration: Step-%d", iteration), "");
  legMain->AddEntry(histoYiledVsRapidity, "Data", "p");
  legMain->AddEntry(fRapDistribution, Form("#it{Gaussian}(0,%2.2f)", fRapDistribution->GetParameter(1)), "l");
  legMain->Draw();
  //--------------------------------------------------------------------------------------------//

  //--------------------------------------------------------------------------------------------//
  //Store the values of the fit. Use 1 (arbitrary) for the normalisation parameter (1st) as a normalisation will be done in later steps.
  gSystem->Exec(Form("mkdir -p Cent-%dto%d/RapShapeIterations/iter-%d", centMin, centMax, iteration));
  canYieldRap->SaveAs(Form("Cent-%dto%d/RapShapeIterations/iter-%d/RapYield.png", centMin, centMax, iteration));
  ofstream outputFile(Form("Cent-%dto%d/RapShapeIterations/iter-%d/values.txt", centMin, centMax, iteration), std::ofstream::out);
  outputFile << "1" << endl;
  outputFile << "0" << endl;
  outputFile << fRapDistribution->GetParameter(1) << endl;
  outputFile.close();
  //gDirectory->GetList()->Delete();
}
//--------------------------------------------------------------------------------------------//