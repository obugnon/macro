/*
For the zeroth iteration the fit is performed on the generated jpsi in a given centrality which must result in parameters similar to the one used for the MC generation (with no centrality dependence: THIS HAS TO BE CHECKED)
TODO: Incorporate this into a more general FitYields.C that can fit either corrected yield or generated distributions
*/

//--------------------------------------------------------------------------------------------//
void FitPtGen(int centMin = 0, int centMax = 10)
{
    PrintMessage("Fitting generated Pt", "-");
    TCanvas *canPtGen = new TCanvas(Form("canPtGen_cent%dto%d", centMin, centMax), "", 900, 800);
    SetCanvasStyle(canPtGen);
    canPtGen->SetLogy();

    //Get the corresponding Acc*Effi file
    TFile *inputMCFile = TFile::Open(Form("Cent-%dto%d/AccEffi/iter-0/AccEffiValues.root", centMin, centMax));
    TH1F *histoGenPt = ((TH1F *)inputMCFile->Get("histoGenPt"));

    // Divide each bin by its width:
    Int_t numberOfBins = histoGenPt->GetNbinsX();
    for (int iBin = 1; iBin <= numberOfBins; iBin++)
    {
        histoGenPt->SetBinContent(iBin, histoGenPt->GetBinContent(iBin) / histoGenPt->GetBinWidth(iBin));
        histoGenPt->SetBinError(iBin, histoGenPt->GetBinError(iBin) / histoGenPt->GetBinWidth(iBin));
    }
    SetHistoStyle(histoGenPt, kBlack, kFullCircle, 0.9, ptAxisTitle, "d#it{N}/d#it{p}_{T}", 1, kFALSE);

    TF1 *fPtDistribution = new TF1("fPtDistribution", pT_shape, histoGenPt->GetXaxis()->GetXmin(), histoGenPt->GetXaxis()->GetXmax(), 4);
    fPtDistribution->SetLineColor(myColorsMap["red"]);
    fPtDistribution->SetParameter(0, histoGenPt->GetMaximum());
    fPtDistribution->SetParameter(1, 4);
    fPtDistribution->SetParameter(2, 2);
    fPtDistribution->SetParameter(3, 3);

    double initParam = 3;
    double paramStep = -0.01;
    int customStepCounter = 0;
    TFitResultPtr fitResultPtr;
    double chi2OverNDF;
    do
    {
        fPtDistribution->SetParameter(1, initParam + paramStep * customStepCounter);
        fitResultPtr = histoGenPt->Fit(fPtDistribution, "SRIQ", "", ptFitRangeMin, ptFitRangeMax);
        chi2OverNDF = fPtDistribution->GetChisquare() / fPtDistribution->GetNDF();
        customStepCounter++;
    } while (fitResultPtr->CovMatrixStatus() != 3 || !fitResultPtr->IsValid() || customStepCounter > 100);

    histoGenPt->GetYaxis()->SetTitleOffset(1.4);
    TLegend *legMain = new TLegend(0.20, 0.27, 0.50, 0.48);
    legMain->SetMargin(0.1);
    legMain->SetFillStyle(0);
    legMain->SetLineColorAlpha(0, 0);
    legMain->SetTextColor(customStepCounter < 100 ? kBlack : kRed);
    legMain->AddEntry("NULL", Form("Centrality: %d-%d %%", centMin, centMax), "");
    legMain->AddEntry(histoGenPt, "MC (generated)", "p");
    legMain->AddEntry(fPtDistribution, Form("C#times#frac{#it{p}_{T}}{  (1+ (#it{p}_{T}/%2.2f)^{%2.2f})^{%2.2f}  }", fPtDistribution->GetParameter(1), fPtDistribution->GetParameter(2), fPtDistribution->GetParameter(3)), "l");
    legMain->Draw();
    //--------------------------------------------------------------------------------------------//

    //--------------------------------------------------------------------------------------------//
    //Store the values of the fit. Use 1 (arbitrary) for the normalisation parameter (1st) as a normalisation will be done in later steps.
    gSystem->Exec(Form("mkdir -p Cent-%dto%d/PtShapeIterations/iter-0", centMin, centMax));
    canPtGen->SaveAs(Form("Cent-%dto%d/PtShapeIterations/iter-0/PtGen.png", centMin, centMax));
    ofstream outputFile(Form("Cent-%dto%d/PtShapeIterations/iter-0/values.txt", centMin, centMax), std::ofstream::out);
    outputFile << "1" << endl;
    outputFile << fPtDistribution->GetParameter(1) << endl;
    outputFile << fPtDistribution->GetParameter(2) << endl;
    outputFile << fPtDistribution->GetParameter(3) << endl;
    outputFile.close();
    //--------------------------------------------------------------------------------------------//
    //gDirectory->GetList()->Delete();
}
//--------------------------------------------------------------------------------------------//

//--------------------------------------------------------------------------------------------//
void FitRapGen(int centMin = 0, int centMax = 10)
{
    PrintMessage("Fitting generated y", "-");

    TCanvas *canRapGen = new TCanvas(Form("canRapGen_cent%dto%d", centMin, centMax), "", 900, 800);
    SetCanvasStyle(canRapGen);

    //Get the corresponding Acc*Effi file
    TFile *inputMCFile = TFile::Open(Form("Cent-%dto%d/AccEffi/iter-0/AccEffiValues.root", centMin, centMax));
    TH1F *histoGenRap = ((TH1F *)inputMCFile->Get("histoGenRapidity"));

    // Divide each bin by its width:
    Int_t numberOfBins = histoGenRap->GetNbinsX();
    for (int iBin = 1; iBin <= numberOfBins; iBin++)
    {
        histoGenRap->SetBinContent(iBin, histoGenRap->GetBinContent(iBin) / histoGenRap->GetBinWidth(iBin));
        histoGenRap->SetBinError(iBin, histoGenRap->GetBinError(iBin) / histoGenRap->GetBinWidth(iBin));
    }
    SetHistoStyle(histoGenRap, kBlack, kFullCircle, 0.9, yAxisTitle, "d#it{N}/d#it{y}", 1, kFALSE);

    TF1 *fRapDistribution = new TF1("fRapDistribution", "[0]*exp(-0.5*(x/[1])**2)", histoGenRap->GetXaxis()->GetXmin(), histoGenRap->GetXaxis()->GetXmax());
    fRapDistribution->SetLineColor(myColorsMap["red"]);
    fRapDistribution->SetParameter(0, histoGenRap->GetMaximum());
    fRapDistribution->SetParameter(1, 3);

    double initWidth = 3;
    double widthStep = -0.01;
    int customStepCounter = 0;
    TFitResultPtr fitResultPtr;
    double chi2OverNDF;
    do
    {
        fRapDistribution->SetParameter(1, initWidth + widthStep * customStepCounter);
        fitResultPtr = histoGenRap->Fit(fRapDistribution, "RISQ", "", rapFitRangeMin, rapFitRangeMax);
        chi2OverNDF = fRapDistribution->GetChisquare() / fRapDistribution->GetNDF();
    } while ((fitResultPtr->CovMatrixStatus() != 3) || (!fitResultPtr->IsValid()) || (customStepCounter > 100));

    TLegend *legMain = new TLegend(0.235, 0.745, 0.455, 0.88);
    legMain->SetMargin(0.1);
    legMain->SetFillStyle(0);
    legMain->SetLineColorAlpha(0, 0);
    legMain->SetTextColor(customStepCounter < 100 ? kBlack : kRed);
    legMain->AddEntry("NULL", Form("Centrality: %d-%d %%", centMin, centMax), "");
    legMain->AddEntry(histoGenRap, "MC (generated)", "p");
    legMain->AddEntry(fRapDistribution, Form("#it{Gaussian}(0,%2.2f)", fRapDistribution->GetParameter(1)), "l");
    legMain->Draw();
    //--------------------------------------------------------------------------------------------//

    //--------------------------------------------------------------------------------------------//
    //Store the values of the fit. Use 1 (arbitrary) for the normalisation parameter (1st) as a normalisation will be done in later steps.
    gSystem->Exec(Form("mkdir -p Cent-%dto%d/RapShapeIterations/iter-0", centMin, centMax));
    canRapGen->SaveAs(Form("Cent-%dto%d/RapShapeIterations/iter-0/RapGen.png", centMin, centMax));
    ofstream outputFile(Form("Cent-%dto%d/RapShapeIterations/iter-0/values.txt", centMin, centMax), std::ofstream::out);
    outputFile << "1" << endl;
    outputFile << "0" << endl;
    outputFile << fRapDistribution->GetParameter(1) << endl;
    outputFile.close();
    //gDirectory->GetList()->Delete();
}
//--------------------------------------------------------------------------------------------//