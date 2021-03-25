/*
 *  ComputeLuminosity.C
 *
 *  Created by Ophelie Bugnon on 24/02/21.
 *
 */
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <iostream>

#include "TFile.h"
#include "TList.h"
#include "THnSparse.h"
#include "TH1.h"
#include "TF1.h"
#include "TString.h"
#include "TMath.h"

Double_t sigmaV0[5] = {4.708, 4.105,  4.073, 3.970, 3.933}; // 2015 puis 295585 - 295589, 295612 - 295615, 295665 - 296198, 296240 - 297624

TString fileLocation = "/Users/obugnon/Documents/ALICE/AnalyseJPsi/AnalysisResults_5TeV";

//__________________________________________________________________________________________________
//__________________________________________________________________________________________________
std::vector<int> runListToVector(TString runList)
{
  std::vector<int> vectorRuns;
  TString strRunList = gSystem->GetFromPipe(Form("cat %s",runList.Data()));
  TObjArray *objRunList = strRunList.Tokenize("\n");
  TIter nextRun(objRunList);
  TObjString *strRun;
  while ((strRun=(TObjString*)nextRun())) 
  {
    TString currentRun = strRun->GetName();
    int runNumber = currentRun.Atoi();
    vectorRuns.push_back(runNumber);
  }
  return vectorRuns;
}
//__________________________________________________________________________________________________
//__________________________________________________________________________________________________
void ComputeLuminosity()
{
    std::vector<int> vecRun = runListToVector(Form("%s/runList_tot.txt", fileLocation.Data()));

    TFile* analysis;
    TList* lCMULhistos; 
    TH1* hCMUL; 
    TList* lCINThistos; 
    TH1* hCINT; 
    TH1* hV0MinCINT; 
    TH1* hMULinCINT; 

    Double_t nCMUL, nV0MinCINT, nMULinCINT, nCMUL_tot;
    Double_t Lumi, Lumi_tot, FNorm, FNorm_tot;
    int j;

    nCMUL_tot=0;
    FNorm_tot=0;
    Lumi_tot=0; 

    for(int i=0; i<vecRun.size(); i++)
    {
      if(vecRun[i] <= 246994) analysis = TFile::Open(Form("%s/AnalysisResults_LHC15o_forFnorm.root", fileLocation.Data()));
      else if(vecRun[i] >= 295584 && vecRun[i] <= 296623) analysis = TFile::Open(Form("%s/AnalysisResults_LHC18q_forFnorm.root", fileLocation.Data()));
      else if(vecRun[i] >= 296690 && vecRun[i] <= 297595) analysis = TFile::Open(Form("%s/AnalysisResults_LHC18r_forFnorm.root", fileLocation.Data()));

      lCMULhistos = (TList*)analysis->Get("EventHistos_CMUL7");
      hCMUL = (TH1*) lCMULhistos->FindObject("fHistoPSEventsPerRun");

      lCINThistos = (TList*)analysis->Get("EventHistos_CINT7");
      hCINT = (TH1*) lCINThistos->FindObject("fHistoPSEventsPerRun");
      hV0MinCINT = (TH1*) lCINThistos->FindObject("fHisto0V0MEventsInCINT7");
      // hMULinCINT = (TH1*) lCINThistos->FindObject("fHisto0MULand0V0MEventsInCINT7");
      // hMULinCINT = (TH1*) lCINThistos->FindObject("fHistoCMULEventsInCINT7");
      hMULinCINT = (TH1*) lCINThistos->FindObject("fHisto0MULEventsInCINT7");

      nCMUL = hCMUL->GetBinContent(hCMUL->FindBin(vecRun[i]));
      nCMUL_tot += nCMUL;
      nV0MinCINT = hV0MinCINT->GetBinContent(hV0MinCINT->FindBin(vecRun[i]));
      nMULinCINT = hMULinCINT->GetBinContent(hMULinCINT->FindBin(vecRun[i]));

      if(nCMUL==0 || nV0MinCINT==0 || nMULinCINT==0) printf("------------Wrong value found for run %i------------\n", vecRun[i]);
      
      if(vecRun[i] <= 246994) j=0;
      else if (vecRun[i] >= 295584 && vecRun[i] <= 295589) j=1;
      else if (vecRun[i] >= 295612 && vecRun[i] <= 295615) j=2;
      else if (vecRun[i] >= 295665 && vecRun[i] <= 296198) j=3;
      else if (vecRun[i] >= 296240 && vecRun[i] <= 297624) j=4;
      
      FNorm = nV0MinCINT/nMULinCINT;
      FNorm_tot += nCMUL*FNorm;

      Lumi = nCMUL*FNorm/(sigmaV0[j]*TMath::Power(10,6));
      Lumi_tot += Lumi; 
    }

    FNorm_tot=FNorm_tot/nCMUL_tot;
    printf("Fnorm for the period is = %.2f \n", FNorm_tot);
    printf("Lumi = %.2f \n", Lumi_tot);
    return;
}
