/*
This macro compare the pt and rapidity dependence of the Acc*Eff that correposnds to two different steps A and B.  
*/
//--------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------//
Double_t DrawOrCompareEffi(int centMin = 0, int centMax = 10, Int_t iteration_A = 2, Int_t iteration_B = 3, TString rapOrPt = "Rap")
{

  TFile *inputAccEffFile_A = TFile::Open(Form("Cent-%dto%d/AccEffi/iter-%d/AccEffiValues.root", centMin, centMax, iteration_A));
  TH1F *histoAccEffi_A = ((TH1F *)inputAccEffFile_A->Get(Form("histoAccEffiVs%s", rapOrPt.Data())));

  TH1F *histoAccEffi_B = nullptr;
  if (iteration_A != iteration_B)
  {
    TFile *inputAccEffFile_B = TFile::Open(Form("Cent-%dto%d/AccEffi/iter-%d/AccEffiValues.root", centMin, centMax, iteration_B));
    histoAccEffi_B = ((TH1F *)inputAccEffFile_B->Get(Form("histoAccEffiVs%s", rapOrPt.Data())));
  }

  //--------------------------------------------------------------------------------------------//
  TCanvas *canAccEffi = new TCanvas(Form("canAccEffiVs%s_%dto%d_Cent-%dto%d", rapOrPt.Data(), iteration_A, iteration_B, centMin, centMax), "", 900, 800);
  SetCanvasStyle(canAccEffi);

  TPad *padMain, *padRatio = nullptr;
  if (iteration_A != iteration_B)
  {
    padMain = new TPad("padMain", "padMain", 0, 0.3, 1, 1, 0);
    padMain->SetBottomMargin(0.);
    padMain->Draw();

    padRatio = new TPad("padRatio", "padRatio", 0, 0, 1, 0.3, 0);
    padRatio->SetBottomMargin(0.4);
    padRatio->SetTopMargin(0.);
    padRatio->Draw();

    padMain->cd();
  }

  SetHistoStyle(histoAccEffi_A, kBlue, kFullCircle, 0.8, rapOrPt == "Rap" ? yAxisTitle : ptAxisTitle, accEffiAxisTitle, 1, kFALSE);
  histoAccEffi_A->Draw();
  if (iteration_A != iteration_B)
  {
    SetHistoStyle(histoAccEffi_B, kRed, kFullCircle, 0.8, "", "", 1, kFALSE);
    histoAccEffi_B->Draw("SAME");
  }

  TLegend *leg;
  if (rapOrPt == "Pt")
  {
    leg = new TLegend(0.21, 0.70, 0.46, 0.85);
  }
  else
  {
    leg = new TLegend(0.42, 0.27, 0.67, 0.42);
  }
  leg->SetMargin(0.1);
  leg->SetFillStyle(0);
  leg->SetLineColorAlpha(0, 0);
  leg->SetTextColor(kBlack);
  leg->AddEntry("NULL", Form("centrality %d-%d %%", centMin, centMax), "");
  leg->AddEntry(histoAccEffi_A, Form("Iteration: Step-%d", iteration_A), "ep");
  if (iteration_A != iteration_B)
  {
    leg->AddEntry(histoAccEffi_B, Form("Iteration: Step-%d", iteration_B), "ep");
  }
  leg->Draw();

  Double_t maxDifference = 0;
  if (iteration_A != iteration_B)
  {
    padRatio->cd();
    TH1F *histoRatio = new TH1F(*histoAccEffi_B);
    histoRatio->Divide(histoAccEffi_A);
    SetHistoStyle(histoRatio, kBlack, kFullCircle, 0.8, rapOrPt == "Rap" ? yAxisTitle : ptAxisTitle, "ratio", 1, kTRUE);
    std::vector<double> userRange = GetAVisuallyGoodYRangeUser(histoRatio->GetMinimum() - histoRatio->GetBinError(histoRatio->GetNbinsX()), histoRatio->GetMaximum() + histoRatio->GetBinError(histoRatio->GetNbinsX()), 1);
    histoRatio->GetYaxis()->SetRangeUser(userRange[0], userRange[1]);
    histoRatio->GetYaxis()->SetNdivisions(505);
    histoRatio->Draw();
    //Draw line at one:
    TLine *lineAtOne = new TLine(histoRatio->GetXaxis()->GetXmin(), 1, histoRatio->GetXaxis()->GetXmax(), 1);
    lineAtOne->SetLineStyle(kDashed);
    lineAtOne->Draw();
    Double_t histoMax = histoRatio->GetBinContent(histoRatio->GetMaximumBin());
    Double_t histoMin = histoRatio->GetBinContent(histoRatio->GetMinimumBin());
    maxDifference = TMath::Max(TMath::Abs(histoMax - 1), TMath::Abs(1 - histoMin));
  }

  if (iteration_A != iteration_B)
  {
    canAccEffi->SaveAs(Form("Cent-%dto%d/AccEffi/iter-%d/AccEffVs%s_ComparisonToIter-%d.png", centMin, centMax, iteration_A, rapOrPt.Data(), iteration_B));
  }
  else
  {
    canAccEffi->SaveAs(Form("Cent-%dto%d/AccEffi/iter-%d/AccEffVs%s.png", centMin, centMax, iteration_A, rapOrPt.Data()));
  }
  //gDirectory->GetList()->Delete();
  return maxDifference;
  //--------------------------------------------------------------------------------------------//
}

double DrawAndCompareAccEffi(int centMin, int centMax, Int_t iteration)
{
  // Draw it alone for all cases
  DrawOrCompareEffi(centMin, centMax, iteration, iteration, "Rap");
  DrawOrCompareEffi(centMin, centMax, iteration, iteration, "Pt");
  PrintMessage(Form("Calculating Acc*Effi"), "-");
  if (iteration != 0)
  {
    PrintMessage(Form("Comparing Acc*Effi with step %d", iteration - 1), "-");
    // Draw it vs iteration -1
    Double_t maxDiffVsRap = DrawOrCompareEffi(centMin, centMax, iteration, iteration - 1, "Rap");
    Double_t maxDiffVsPt = DrawOrCompareEffi(centMin, centMax, iteration, iteration - 1, "Pt");
    return TMath::Max(maxDiffVsRap, maxDiffVsPt);
  }
  return 0;
}