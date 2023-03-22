
std::vector<vector<int>> TxtToVector(TString strEventsPerRun = "0 0", TString outerSplitter = "\n", TString innerSplitter = ",")
{
    std::vector<vector<int>> vectorBins;
    std::vector<int> vectorSingleBin;

    TObjArray *objBins = strEventsPerRun.Tokenize(outerSplitter);
    TIter nextBin(objBins);
    TObjString *strBin;
    while ((strBin = (TObjString *)nextBin()))
    {
        TString currentBin = strBin->GetName();
        vectorSingleBin.clear();

        TObjArray *objValues = currentBin.Tokenize(innerSplitter);
        vectorSingleBin.push_back(((TObjString *)objValues->UncheckedAt(0))->String().Atoi());
        vectorSingleBin.push_back(((TObjString *)objValues->UncheckedAt(1))->String().Atoi());
        vectorBins.push_back(vectorSingleBin);
    }
    return vectorBins;
}
TH1I *FillNCMULEventsPerRun(TString histoName = "histoCMULEventPerRun", TString runLitsPath = "./CMULtot.txt")
{
    std::vector<vector<int>> vectorOfEventsPerRun = TxtToVector(gSystem->GetFromPipe(Form("cat %s", runLitsPath.Data())), "\n", ",");
    int numberOfRuns = (int)vectorOfEventsPerRun.size();
    cout << Form("Filling CMUL per run histogram for %d runs", numberOfRuns) << endl;

    // Hardcoded as you might want to have fix-sized histograms regradless of the input
    const int firstRun = 244918;
    const int lastRun = 297595;
    int maxNumberOfRuns = lastRun - firstRun;

    TH1I *histo = new TH1I(histoName, "", maxNumberOfRuns, firstRun, lastRun);
    for (int iRun = 0; iRun < numberOfRuns; iRun++)
    {
        histo->SetBinContent(vectorOfEventsPerRun[iRun][0] - firstRun + 1, vectorOfEventsPerRun[iRun][1]);
    }
    return histo;
}
TH1F *FillNJpsiVsCentrality(TString histoName = "histoNJpsiVsCentrality")
{
    Double_t arrayCentralityBins[] = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90};
    int numberOfCentralityBins = sizeof(arrayCentralityBins) / sizeof(arrayCentralityBins[0]) - 1;

    Double_t arrayNJpsi[][2] = {
        {361984., 3262.},
        {227930., 2325.},
        {146050., 1664.},
        {82286., 1001.},
        {48791., 642.},
        {25160., 364.},
        {12897., 212.},
        {5712., 120.},
        {2250., 64.}};
    TH1F *histo = new TH1F(histoName, "", numberOfCentralityBins, &arrayCentralityBins[0]);
    for (int iBin = 1; iBin <= numberOfCentralityBins; iBin++)
    {
        histo->SetBinContent(iBin, arrayNJpsi[iBin - 1][0]);
        histo->SetBinError(iBin, arrayNJpsi[iBin - 1][1]);
    }
    return histo;
}
void FillHistograms(TString outputFileName = "CMUL_Centrality_Weights.root")
{
    TFile *outputFile = new TFile(outputFileName, "recreate");
    outputFile->cd();

    FillNCMULEventsPerRun();
    FillNJpsiVsCentrality();
    outputFile->Write();
    delete outputFile;
}