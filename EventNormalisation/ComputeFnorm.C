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

// Double_t sigmaV0 = 4.706;//2015
// Double_t sigmaV0 = 3.933;//2018
Double_t sigmaV0[5] = {4.708, 4.105,  4.073, 3.970, 3.933}; // 2015 puis 295585 - 295589, 295612 - 295615, 295665 - 296198, 296240 - 297624
Double_t dsigmaV0 = 0.004;

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
void ComputFnorm(TString period)
{
    std::vector<int> vecRun = runListToVector(Form("%s/runList_%s.txt", fileLocation.Data(), period.Data()));

    TCanvas* cCINTRatio = new TCanvas("cCINTratio","");
    TH1D* hCINTRatio = new TH1D("hCINTRatio","",vecRun.size(),1,vecRun.size());
    hCINTRatio->Sumw2();
    hCINTRatio->GetXaxis()->SetTitle("run");
    hCINTRatio->SetLabelSize(0.03,"X");
    hCINTRatio->SetTitleOffset(1.60, "X");
    hCINTRatio->SetLabelSize(0.03,"Y");
    hCINTRatio->GetYaxis()->SetTitle("N_{CINT}/N_{CINT&V0M}");
    hCINTRatio->SetLineColor(9);


    TFile* analysis = TFile::Open(Form("%s/AnalysisResults_%s_forFnorm.root", fileLocation.Data(), period.Data()));

    TList* lCMULhistos = (TList*)analysis->Get("EventHistos_CMUL7");
    TH1* hCMUL = (TH1*) lCMULhistos->FindObject("fHistoPSEventsPerRun");

    TList* lCINThistos = (TList*)analysis->Get("EventHistos_CINT7");
    TH1* hCINT = (TH1*) lCINThistos->FindObject("fHistoPSEventsPerRun");
    TH1* hV0MinCINT = (TH1*) lCINThistos->FindObject("fHisto0V0MEventsInCINT7");
    // TH1* hMULinCINT = (TH1*) lCINThistos->FindObject("fHisto0MULand0V0MEventsInCINT7");
    // TH1* hMULinCINT = (TH1*) lCINThistos->FindObject("fHistoCMULEventsInCINT7");
    TH1* hMULinCINT = (TH1*) lCINThistos->FindObject("fHisto0MULEventsInCINT7");

    Double_t nCMUL, nV0MinCINT, nMULinCINT, nCMUL_tot, nCINT, errCMUL, errV0MinCINT, errMULinCINT, errCMUL_tot, errCINT;
    Double_t Lumi_tot, FNorm, errFNorm, FNorm_tot, errFNorm_tot;

    Double_t FNorm_test, errFNorm_test;
    Double_t FNorm_transition, errFNorm_transition, FNorm_tot_transition;
    nCMUL_tot=0;
    FNorm_tot=0;
    errFNorm_tot=0;
    Lumi_tot=0; 

    for(int i=0; i<vecRun.size(); i++)
    {
      nCMUL = hCMUL->GetBinContent(hCMUL->FindBin(vecRun[i]));
      errCMUL = hCMUL->GetBinError(hCMUL->FindBin(vecRun[i]));
      nCMUL_tot += nCMUL;
      // printf("Run %i : %i CMUL events ---> Total of %i \n", vecRun[i], nCMUL, nCMUL_tot);
      nV0MinCINT = hV0MinCINT->GetBinContent(hV0MinCINT->FindBin(vecRun[i]));
      errV0MinCINT = hV0MinCINT->GetBinError(hV0MinCINT->FindBin(vecRun[i]));
      nMULinCINT = hMULinCINT->GetBinContent(hMULinCINT->FindBin(vecRun[i]));
      errMULinCINT = hMULinCINT->GetBinError(hMULinCINT->FindBin(vecRun[i]));
      nCINT = hCINT->GetBinContent(hCINT->FindBin(vecRun[i]));
      errCINT = hCINT->GetBinError(hCINT->FindBin(vecRun[i]));
      
      //Vrai FNorm
      FNorm = nV0MinCINT/nMULinCINT;
      errFNorm = TMath::Sqrt((errV0MinCINT/nV0MinCINT)*(errV0MinCINT/nV0MinCINT)+(errMULinCINT/nMULinCINT)*(errMULinCINT/nMULinCINT))*FNorm;
      FNorm_tot += nCMUL*FNorm;
      errFNorm_tot += nCMUL*errFNorm;

      FNorm_test = nV0MinCINT/nCINT;
      errFNorm_test = TMath::Sqrt((errV0MinCINT/nV0MinCINT)*(errV0MinCINT/nV0MinCINT)+(errCINT/nCINT)*(errCINT/nCINT))*FNorm_test;
      
      FNorm_transition = nCINT/nV0MinCINT;
      FNorm_tot_transition += nCMUL*FNorm_transition;

      hCINTRatio->SetBinContent(i+1,FNorm_test);
      hCINTRatio->SetBinError(i+1,errFNorm_test);
   		hCINTRatio->GetXaxis()->SetBinLabel(i+1, Form("%i",vecRun[i]));

    }

    hCINTRatio->Draw();

    FNorm_tot=FNorm_tot/nCMUL_tot;
    FNorm_tot_transition=FNorm_tot_transition/nCMUL_tot;
    printf("Fnorm for the period is = %.2f \n", FNorm_tot);
    printf("Test Fnorm for the period is = %.2f \n", FNorm_tot_transition);

    Lumi_tot = nCMUL_tot*FNorm_tot/(sigmaV0*TMath::Power(10,6));
    printf("Lumi = %.2f \n", Lumi_tot);
}
//__________________________________________________________________________________________________
//__________________________________________________________________________________________________
void VerifieNumberOfCMUL(TString period)
{
    TFile* analysis = TFile::Open(Form("%s/AnalysisResults_%s_forFnorm.root", fileLocation.Data(), period.Data()));
    TFile* analysis2 = TFile::Open(Form("%s/Individual_data_prod/AnalysisResults_15o_AOD229.root", fileLocation.Data()));

    TList* lCMULhistos = (TList*)analysis->Get("EventHistos_CMUL7");
    TH1* hCMUL = (TH1*) lCMULhistos->FindObject("fHistoPSEventsPerRun");
    TList* lCMULhistos2 = (TList*)analysis2->Get("EventHistos_CMUL7");
    TH1* hCMUL2 = (TH1*) lCMULhistos2->FindObject("fHistoPSEventsPerRun");

    std::vector<int> vecRun = runListToVector(Form("%s/runList_%s.txt", fileLocation.Data(), period.Data()));

    Int_t nCMUL, nCMUL_2, nCMUL_tot;
    Double_t diff;
    

    for(int i=0; i<vecRun.size(); i++)
    {
        nCMUL = hCMUL->GetBinContent(hCMUL->FindBin(vecRun[i]));
        nCMUL_2 = hCMUL2->GetBinContent(hCMUL2->FindBin(vecRun[i]));
        diff = TMath::Abs(nCMUL-nCMUL_2);
        if(nCMUL < nCMUL_2) printf("In run %i there is %i CMUL events in first file and %i in the second --> diff is %.2f %%\n", vecRun[i], nCMUL_2, nCMUL, diff/nCMUL_2*100);
        // if(nCMUL < nCMUL_2) printf("%i ", vecRun[i]);
        nCMUL_tot+=nCMUL;
    }
    printf("%i ", nCMUL_tot);
}
//__________________________________________________________________________________________________
//__________________________________________________________________________________________________
void ComputStandartFnorm(TString period)
{
    TFile* analysis = TFile::Open(Form("%s/AnalysisResults_%s_forFnorm.root", fileLocation.Data(), period.Data()));

    TList* lCMULhistos = (TList*)analysis->Get("EventHistos_CMUL7");
    if(!lCMULhistos) printf("Cannot find CMUL histogramms\n"); 
    TH1* hCMUL = (TH1*) lCMULhistos->FindObject("fHistoPSEventsPerRun");
    if(!hCMUL) printf("Cannot find CMUL histo\n"); 

    TList* lCINThistos = (TList*)analysis->Get("EventHistos_CINT7");
    if(!lCINThistos) printf("Cannot find CINT histogramms\n"); 
    TH1* hMSLinCINT = (TH1*) lCINThistos->FindObject("fHisto0MSLEventsInCINT7");
    TH1* hCINT = (TH1*) lCINThistos->FindObject("fHistoPSEventsPerRun");
    if(!hMSLinCINT || !hCINT) printf("Cannot find CINT histo\n"); 

    TList* lCMSLhistos = (TList*)analysis->Get("EventHistos_CMSL7");
    if(!lCMSLhistos) printf("Cannot find CMSL histogramms\n"); 
    TH1* hMULinCMSL = (TH1*) lCMSLhistos->FindObject(" fHisto0MULEventsInCMSL");
    TH1* hCMSL = (TH1*) lCMSLhistos->FindObject("fHistoPSEventsPerRun");
    if(!hMULinCMSL) printf("Cannot find MULinCMSL histo\n"); 
    if(!hCMSL) printf("Cannot find CMSL histo\n"); 

    std::vector<int> vecRun = runListToVector(Form("%s/runList_%s.txt", fileLocation.Data(), period.Data()));

    Double_t nCMUL, nMSLinCINT, nCINT, nMULinCMSL, nCMSL, nCMUL_tot;
    Double_t FNorm, FNorm_tot;
    nCMUL_tot=0;
    FNorm_tot=0;
    
    for(int i=0; i<vecRun.size(); i++)
    {
        nCMUL = hCMUL->GetBinContent(hCMUL->FindBin(vecRun[i]));
        nCMUL_tot += nCMUL;
        // printf("Run %i : %i CMUL events ---> Total of %i \n", vecRun[i], nCMUL, nCMUL_tot);
        nMSLinCINT = hMSLinCINT->GetBinContent(hMSLinCINT->FindBin(vecRun[i]));
        nCINT = hCINT->GetBinContent(hCINT->FindBin(vecRun[i]));

        nMULinCMSL = hMULinCMSL->GetBinContent(hMULinCMSL->FindBin(vecRun[i]));
        nCMSL = hCMSL->GetBinContent(hCMSL->FindBin(vecRun[i]));
        
        FNorm =  ((nCINT/nMSLinCINT)*(nCMSL/nMULinCMSL));
        // printf("Fnorm = %.2.f \n", FNorm);
        FNorm_tot += nCMUL*FNorm;
    }
    FNorm_tot=FNorm_tot/nCMUL_tot;
    printf("Fnorm for the period is = %.2f \n", FNorm_tot);
    Double_t Lumi;
    Lumi = nCMUL_tot*11.88/(7.67*TMath::Power(10,6));
    printf("Lumi = %.2f \n", Lumi);

}
// // __________________________________________________________________________________________________
// // __________________________________________________________________________________________________
// void PileUpFactor(Int_t runNumber, const char* period, Double_t& value, Double_t& error)
// {
//   TFile* analysis = TFile::Open(Form("%s/AnalysisResults_%s_forFnorm.root", fileLocation.Data(), period.Data()));
//   TList* lCINThistos = (TList*)analysis->Get("EventHistos_CINT7");

//   AliAnalysisTriggerScalers* triggerScalers = new AliAnalysisTriggerScalers(Form("%s/runList_%s.txt", fileLocation.Data(), period.Data()),"alien://folder=/alice/data/2015/OCDB");
  
//   //Get RunList information
//   std::vector<int> vecRun = runListToVector(Form("%s/runList_%s.txt", fileLocation.Data(), period.Data()));
  
//   Double_t nPurity, errorPurity, nPileup, errorPileup;
//   GetPurityFactor(eventHistos_CINT7, runNumber, nPurity, errorPurity);
//   triggerScalers->GetPileUpFactor(runNumber, "CINT7-B-NOPF-MUFAST", nPurity, nPileup, errorPileup);
  
// }
// //__________________________________________________________________________________________________
// //__________________________________________________________________________________________________
// void GetPurityFactor(TString period, Int_t runNumber, Double_t& value, Double_t& error)
// {
//   TFile* analysis = TFile::Open(Form("%s/AnalysisResults_%s_forFnorm.root", fileLocation.Data(), period.Data()));
//   value = error = 0.0;
  
//   TH1* hEventsPS = (TH1*)eventHistos->FindObject("fHistoPSEventsPerRun");
//   TH1* hEventsALL = (TH1*)eventHistos->FindObject("fHistoEventsBeforePSPerRun");

//   Double_t nPS = hEventsPS->GetBinContent(hEventsPS->FindBin(runNumber));
//   if (!nPS) {
//     printf(Form("Cannot work without data for run %d", runNumber));
//     return;
//   }
//   Double_t errorPS = hEventsPS->GetBinError(hEventsPS->FindBin(runNumber));
  
//   Double_t nALL = hEventsALL->GetBinContent(hEventsALL->FindBin(runNumber));
//   if (!nALL) {
//     printf(Form("Cannot work without data for run %d", runNumber));
//     return;
//   }
//   Double_t errorALL = hEventsALL->GetBinError(hEventsALL->FindBin(runNumber));
   
//   value = (nPS/nALL);
//   error = TMath::Sqrt(errorPS*errorPS/(nPS*nPS) + errorALL*errorALL/(nALL*nALL))*(nPS/nALL);
   
//   return;
// }