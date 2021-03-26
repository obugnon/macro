/*
This macro takes as inputs a tree that contains dimuon information (generated and reconstructed) and fill them in histograms. The Acc*Effi is calculated as well.
This macros has been modified to weight the pt and rapidity shape via an iterative procedure. It also need the centrality bin as an argument.
*/

enum enumDimuon
{
  kDimuonInvMass,
  kDimuonPt,
  kDimuonRapidity
};

//-------------------------------------------------------------------------------------------------------------------------------------//
void FillHistoAndGetAccEffi(Int_t centMin = 0, Int_t centMax = 10, Int_t iteration = 0, Bool_t isSpecifBins=kFALSE)//modified
{
  Double_t specificArray[4]={0, 0.3, 1, 8}; //modified

  //-------------------------------------------------------------------------------------------------------------------------------------//

  int numberOfRapidityBins;
  const Double_t *arrayRapidityBins = GetBinsArray(centMin, centMax, "Rap", numberOfRapidityBins);

  int numberOfPtBins;
  const Double_t *arrayPtBins;
  if(isSpecifBins){ 
    numberOfPtBins=3;
    arrayPtBins = specificArray;
  } 
  else arrayPtBins = GetBinsArray(centMin, centMax, "Pt", numberOfPtBins);

  //-------------------------------------------------------------------------------------------------------------------------------------//

  //-------------------------------------------------------------------------------------------------------------------------------------//
  //Input files
  TString arrayOfPeriods[] = {"LHC16e2", "LHC16e2plus", "LHC19a2_18q", "LHC19a2_18r"};
  int numberOfPeriods = sizeof(arrayOfPeriods) / sizeof(arrayOfPeriods[0]);
  TChain *mergedTree = new TChain("eventsTree");
  //Loop over the priods and for each one, access the different .root files and attach them to the chain:
  for (int iPeriod = 0; iPeriod < numberOfPeriods; iPeriod++)
  {
    TString periodDirPath;
    periodDirPath.Form("%s/%s/", inputMCFilesDir.Data(), arrayOfPeriods[iPeriod].Data());
    //Loop over files in subdirectory
    TSystemDirectory periodDir(periodDirPath, periodDirPath);
    TList *listOfFilesPerPeriod = periodDir.GetListOfFiles();
    if (listOfFilesPerPeriod)
    {
      TSystemFile *fileInPeriod;
      TString fileInPeriodName;
      TIter next(listOfFilesPerPeriod);
      while ((fileInPeriod = (TSystemFile *)next()))
      {
        fileInPeriodName = fileInPeriod->GetName();
        if (!fileInPeriod->IsDirectory() && fileInPeriodName.EndsWith("New.root"))
        {
          if (verbose != 0)
          {
            cout << Form("Attaching %s of %s", fileInPeriodName.Data(), arrayOfPeriods[iPeriod].Data()) << endl;
          }
          TString fileInPeriodPath;
          fileInPeriodPath.Form("%s/%s", periodDirPath.Data(), fileInPeriodName.Data());
          mergedTree->Add(fileInPeriodPath);
        }
      }
    }
    //end of loop over files in directory
  }
  //------------------------------------------------------------------------------------------------------//

  Int_t nEntries = mergedTree->GetEntries();

  std::vector<double> *tempoVectorGenDimuon = 0;
  std::vector<double> *tempoVectorRecDimuon = 0;
  Int_t runNumber;
  Float_t centrality;
  mergedTree->SetBranchAddress("GenDimuon", &tempoVectorGenDimuon);
  mergedTree->SetBranchAddress("RecDimuon", &tempoVectorRecDimuon);
  mergedTree->SetBranchAddress("runNumber", &runNumber);
  mergedTree->SetBranchAddress("centrality", &centrality);
  //-------------------------------------------------------------------------------------------------------------------------------------//

  //-------------------------------------------------------------------------------------------------------------------------------------//

  //Histograms for generated
  TH1F *histoGenPt = new TH1F("histoGenPt", "", numberOfPtBins, &arrayPtBins[0]);
  histoGenPt->Sumw2();

  TH1F *histoGenRapidity = new TH1F("histoGenRapidity", "", numberOfRapidityBins, &arrayRapidityBins[0]);
  histoGenRapidity->Sumw2();

  //Histograms for reconstructed

  TH1F *histoRecPt = new TH1F("histoRecPt", "", numberOfPtBins, &arrayPtBins[0]);
  histoRecPt->Sumw2();

  TH1F *histoRecRapidity = new TH1F("histoRecRapidity", "", numberOfRapidityBins, &arrayRapidityBins[0]);
  histoRecRapidity->Sumw2();

  //-------------------------------------------------------------------------------------------------------------------------------------//

  //-------------------------------------------------------------------------------------------------------------------------------------//
  //Pt and y Weight histograms (for iteration step > 0)
  TH1 *histoPtWeight = NULL;
  TH1 *histoRapWeight = NULL;

  if (iteration > 0)
  {
    histoPtWeight = GetPtWeight(centMin, centMax, iteration, 1);
    histoRapWeight = GetRapWeight(centMin, centMax, iteration, 1);
    for (int iIteration = iteration - 1; iIteration > 0; iIteration--)
    {
      histoPtWeight->Multiply(GetPtWeight(centMin, centMax, iIteration));
      histoRapWeight->Multiply(GetRapWeight(centMin, centMax, iIteration));
    }
  }

  if (!histoNeventPerRunGen)
  {
    PrintMessage("Finding the number of simulated events per run and centrality", "-");
    histoNeventPerRunGen = static_cast<TH1I *>(histoCMULEventPerRun->Clone("histoNEventsPerRunGen"));
    histoNeventPerRunGen->Reset();
    histoNJpsiVsCentralityGen = static_cast<TH1F *>(histoNJpsiVsCentrality->Clone("histoNJpsiVsCentralityGen"));
    histoNJpsiVsCentralityGen->Reset();
    for (Int_t iEvent = 0; iEvent < nEntries; iEvent++)
    {
      if (iEvent % 100 == 0)
      {
        printProgress((double)iEvent / nEntries);
      }
      mergedTree->GetEntry(iEvent);
      histoNeventPerRunGen->Fill(runNumber);
      histoNJpsiVsCentralityGen->Fill(centrality);
    }
    cout << endl;

    histoNeventPerRunGen->Scale(totalCMULEvents / histoNeventPerRunGen->Integral());
    histoCMULEventPerRun->Divide(histoNeventPerRunGen);
    histoNJpsiVsCentralityGen->Scale(totalNumberOfJpsi / histoNJpsiVsCentralityGen->Integral());
    histoNJpsiVsCentrality->Divide(histoNJpsiVsCentralityGen);
  }

  //-------------------------------------------------------------------------------------------------------------------------------------//
  PrintMessage("Filling generated and reconstructed histograms", "-");
  //-------------------------------------------------------------------------------------------------------------------------------------//
  for (Int_t iEvent = 0; iEvent < nEntries; iEvent++)
  {

    if (iEvent % 100 == 0)
    {
      printProgress((double)iEvent / nEntries);
    }
    mergedTree->GetEntry(iEvent);

    if (centrality < centMin || centrality >= centMax)
    {
      continue;
    }

    // Double_t weightDimuon = 1;
    //Weight CMUL per run and CMUL vs centrality
    // Double_t cmulFractionPerRun = histoCMULEventPerRun->GetBinContent(histoCMULEventPerRun->GetXaxis()->FindBin(runNumber)) / totalCMULEvents;
    // Double_t njpsiVsCentralityFraction = histoNJpsiVsCentrality->GetBinContent((int)(centrality / 10.) + 1) / totalNumberOfJpsi;
    // weightDimuon = weightDimuon / (cmulFractionPerRun * njpsiVsCentralityFraction);
    double weightDimuon = histoCMULEventPerRun->GetBinContent(histoCMULEventPerRun->GetXaxis()->FindBin(runNumber)) * histoNJpsiVsCentrality->GetBinContent((int)(centrality / 10.) + 1);

    if (tempoVectorGenDimuon->at(kDimuonRapidity) < arrayRapidityBins[numberOfRapidityBins] && tempoVectorGenDimuon->at(kDimuonRapidity) > arrayRapidityBins[0] && tempoVectorGenDimuon->at(kDimuonPt) > arrayPtBins[0] && tempoVectorGenDimuon->at(kDimuonPt) < arrayPtBins[numberOfPtBins])
    {
      if (iteration > 0)
      {
        weightDimuon *= histoPtWeight->GetBinContent(histoPtWeight->GetXaxis()->FindBin(tempoVectorGenDimuon->at(kDimuonPt))) * histoRapWeight->GetBinContent(histoRapWeight->GetXaxis()->FindBin(tempoVectorGenDimuon->at(kDimuonRapidity)));
      }

      histoGenPt->Fill(tempoVectorGenDimuon->at(kDimuonPt), weightDimuon);
      histoGenRapidity->Fill(tempoVectorGenDimuon->at(kDimuonRapidity), weightDimuon);
    }
    if (tempoVectorRecDimuon->size() > 0 && weightDimuon > 0)
    {
      if (tempoVectorRecDimuon->at(kDimuonRapidity) < arrayRapidityBins[numberOfRapidityBins] && tempoVectorRecDimuon->at(kDimuonRapidity) > arrayRapidityBins[0])
      {

        histoRecPt->Fill(tempoVectorRecDimuon->at(kDimuonPt), weightDimuon);
        histoRecRapidity->Fill(tempoVectorRecDimuon->at(kDimuonRapidity), weightDimuon);
      }
    }
  }
  //-------------------------------------------------------------------------------------------------------------------------------------//

  //-------------------------------------------------------------------------------------------------------------------------------------//
  //Acc*Effi calculation

  TH1F *histoAccEffiVsPt = new TH1F(*histoRecPt);
  histoAccEffiVsPt->SetName("histoAccEffiVsPt");
  histoAccEffiVsPt->Divide(histoGenPt);

  TH1F *histoAccEffiVsRap = new TH1F(*histoRecRapidity);
  histoAccEffiVsRap->SetName("histoAccEffiVsRap");
  histoAccEffiVsRap->Divide(histoGenRapidity);
  cout << endl
       << endl;

  //Output File
  gSystem->Exec(Form("mkdir -p Cent-%dto%d/AccEffi/iter-%d", centMin, centMax, iteration));
  TFile *outputFile;
  if(isSpecifBins) outputFile = TFile::Open(Form("Cent-%dto%d/AccEffi/iter-%d/AccEffiValuesOtherRange.root", centMin, centMax, iteration), "recreate");
  else outputFile = TFile::Open(Form("Cent-%dto%d/AccEffi/iter-%d/AccEffiValues.root", centMin, centMax, iteration), "recreate");
  histoGenPt->Write();
  histoGenRapidity->Write();
  histoRecPt->Write();
  histoRecRapidity->Write();
  histoAccEffiVsPt->Write();
  histoAccEffiVsRap->Write();
  delete outputFile;
  //-------------------------------------------------------------------------------------------------------------------------------------//
}
