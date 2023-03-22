#include "RootIncluders.C"
#include "Common.C"

enum enumDimuon
{
  kDimuonInvMass,
  kDimuonPt,
  kDimuonRapidity
};

//-------------------------------------------------------------------------------------------------------------------------------------//
void FillHistoAndGetAccEffi()
{
  //-------------------------------------------------------------------------------------------------------------------------------------//
  //Input files
  TString arrayOfPeriods[] = {"LHC20g13", "LHC20j4_18q", "LHC20j4_18r"};
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
  //Output File
  gSystem->Exec("mkdir -p AccEffi_withPtCuts_rec");
  TFile *outputFile = TFile::Open("AccEffi_withPtCuts_rec/AccEffiValues.root", "recreate");

  // gSystem->Exec("mkdir -p AccEffi_withPtCuts_all");
  // TFile *outputFile = TFile::Open("AccEffi_withPtCuts_all/AccEffiValues.root", "recreate");

  // gSystem->Exec("mkdir -p AccEffi_withoutPtCuts_all");
  // TFile *outputFile = TFile::Open("AccEffi_withoutPtCuts_all/AccEffiValues.root", "recreate");

  const Int_t numberOfCentBins=5;
  const Double_t arrayCentBins[numberOfCentBins+1]={0, 10, 30, 50, 70, 90};
  //Histograms for generated
  TH1F *histoGenCent = new TH1F("histoGenCent", "", numberOfCentBins, &arrayCentBins[0]);
  histoGenCent->Sumw2();

  //Histograms for reconstructed
  TH1F *histoRecCent = new TH1F("histoRecCent", "",  numberOfCentBins, &arrayCentBins[0]);
  histoRecCent->Sumw2();

  //Histograms for reconstructed
  TH1F *histoCent = new TH1F("histoCent", "",  90, 0, 90);
  histoCent->Sumw2();

  //-------------------------------------------------------------------------------------------------------------------------------------//

  //-------------------------------------------------------------------------------------------------------------------------------------//
  //Run Weight histograms

  //Input files for CMUL weighting:
  TH1I *histoNeventPerRunGen = nullptr;
  TFile *fileData = TFile::Open("./DataWeights/CMUL_Centrality_Weights.root");
  TH1I *histoCMULEventPerRun = ((TH1I *)fileData->Get("histoCMULEventPerRun"));
  Double_t totalCMULEvents = histoCMULEventPerRun->Integral();
  Double_t totalGenEvents;
  if (!histoNeventPerRunGen)
  {
    PrintMessage("Finding the number of events per run in the simulations", "-");
    histoNeventPerRunGen = static_cast<TH1I *>(histoCMULEventPerRun->Clone("histoNEventsPerRunGen"));
    histoNeventPerRunGen->Reset();
    
    for (Int_t iEvent = 0; iEvent < nEntries; iEvent++)
    {
      if (iEvent % 100 == 0)
      {
        printProgress((double)iEvent / nEntries);
      }
      mergedTree->GetEntry(iEvent);
      histoNeventPerRunGen->Fill(runNumber);
    }
    totalGenEvents = histoNeventPerRunGen->Integral();
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

    histoCent->Fill(centrality);
    
    Double_t nGenEvent = histoNeventPerRunGen->GetBinContent(histoNeventPerRunGen->GetXaxis()->FindBin(runNumber));
    Double_t nCMULEvent = histoCMULEventPerRun->GetBinContent(histoCMULEventPerRun->GetXaxis()->FindBin(runNumber));

    Double_t weightDimuon = 1;
    weightDimuon = totalGenEvents/nGenEvent*nCMULEvent/totalCMULEvents;

    // if (tempoVectorGenDimuon->at(kDimuonPt) > 0. && tempoVectorGenDimuon->at(kDimuonPt) < 0.3 && tempoVectorGenDimuon->at(kDimuonRapidity) > -4. && tempoVectorGenDimuon->at(kDimuonRapidity) < -2.5)
    if (tempoVectorGenDimuon->at(kDimuonPt) > 0. && tempoVectorGenDimuon->at(kDimuonRapidity) > -4. && tempoVectorGenDimuon->at(kDimuonRapidity) < -2.5)
    {
      histoGenCent->Fill(centrality, weightDimuon);
    }
    if (tempoVectorRecDimuon->size() > 0 && weightDimuon > 0)
    {
      if (tempoVectorRecDimuon->at(kDimuonPt) > 0. && tempoVectorRecDimuon->at(kDimuonPt) < 0.3 && tempoVectorRecDimuon->at(kDimuonRapidity) > -4. && tempoVectorRecDimuon->at(kDimuonRapidity) < -2.5)
      // if (tempoVectorRecDimuon->at(kDimuonPt) > 0. && tempoVectorRecDimuon->at(kDimuonRapidity) > -4. && tempoVectorRecDimuon->at(kDimuonRapidity) < -2.5)
      {
        histoRecCent->Fill(centrality, weightDimuon);
      }
    }
  }
  //-------------------------------------------------------------------------------------------------------------------------------------//

  //-------------------------------------------------------------------------------------------------------------------------------------//
  //Acc*Effi calculation
  outputFile->cd();
  TH1F *histoAccEffiVsCent = new TH1F(*histoRecCent);
  histoAccEffiVsCent->SetName("histoAccEffiCent");
  histoAccEffiVsCent->Divide(histoGenCent);

  cout << endl
       << endl;
  //-------------------------------------------------------------------------------------------------------------------------------------//

  outputFile->Write();
  outputFile->Close();
}
