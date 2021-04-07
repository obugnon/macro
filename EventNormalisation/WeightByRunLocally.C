/*
 *  WeightByRunLocally.C
 *
 *  Created by Ophelie Bugnon on 11/03/20.
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
#include "TSystem.h"


void FiletoVector(TString file, std::vector<Int_t> &tab1, std::vector<Int_t> &tab2)
{
    FILE* fichier = NULL;
    Double_t entry;

    fichier = fopen(Form("%s", file.Data()), "r");
    if (fichier == NULL) 
    {
        printf("No file found for run list");
        return;
    }    
    char nextChar = 'a';
    TString nextEntry = "";

    do
    {
        nextChar = fgetc(fichier);
        if(nextChar == ' ')
        {
            //nextEntry to double dans vector 1
            entry = nextEntry.Atoi();
            tab1.push_back(entry);
            nextEntry = "";
        }
        else if(nextChar == '\n')
        {
            //nextEntry to double dans vector 2
            entry = nextEntry.Atoi();
            tab2.push_back(entry);
            nextEntry = "";
        }
        else
        {
            nextEntry += nextChar;
        } 
    } while (nextChar != EOF);
    
    fclose(fichier);
    return;
} 


void WeightByRunLocally(TString fileLocation="~/Documents/ALICE/Data/2016/LHC16e2tot_weightedByRun", TString fileName = "AnalysisResults.root", TString runListandEvents = "CMULevents_LHC15o.txt")
{
    std::vector<Int_t> runList;
    std::vector<Int_t> nCMUL;
    Int_t totalCMUL = 400381507; //LHC15o=126395768, LHC18q=111006434, LHC18r=162979305, Total=400381507
    Int_t totalGeneratedEvents = 13464798; //LHC16e2=5417445, LHC16e2plus=8047353, LHC16e2tot=13464798, LHC19a2=101912143, LHC18c2incoh=11610906, LHC18c2coh=11620915, LHC20a6coh=26539830, LHC20a6incoh=26259793
    Int_t nGeneratedEvents; 
    Double_t weight;

  
    FiletoVector(Form("%s", runListandEvents.Data()), runList, nCMUL);
    if (runList.size() != nCMUL.size()) 
    {
        printf("run list size is different than nCMUL size");
        return;
    }    
  
    // loop over runs
    Int_t nMissingFiles = 0;
    TList* generatedHistos;
    TList* reconstructedHistos;
    TList* eventHistos;
    THnSparse* hDimuon;
    THnSparse* hSingleMuon;
    TH1I* hEventsPerRun;

    TList* generatedHistosLight;
    TList* reconstructedHistosLight;
    THnSparse* hDimuonGenerated;
    THnSparse* hSingleMuonGenerated;
    THnSparse* hDimuonReconstructed;
    THnSparse* hSingleReconstructed;

    int n=0;

    for(int i=0; i<runList.size(); i++)
    {
        TFile* analysis = new TFile(Form("%s/runs/%i/%s", fileLocation.Data(), runList[i], fileName.Data()), "UPDATE");
        if (!analysis) 
        {
            printf(Form("Missing file run %i",runList[i])); 
            n++;
            continue;
        }    
    
        //compute weight
        eventHistos = (TList*)analysis->Get("EventHistos_CAny");
        // hEventsPerRun = (TH1I*)eventHistos->FindObject("fHistoTotalEventsPerRun");
        hEventsPerRun = (TH1I*)eventHistos->FindObject("fHistoPSEventsPerRun");
        nGeneratedEvents=hEventsPerRun->GetEntries();
        weight = (Double_t)totalGeneratedEvents/nGeneratedEvents*nCMUL[i]/totalCMUL;
        cout << "ncmul = " << nCMUL[i] << " and ncany = " << nGeneratedEvents << " for weight = " << weight << endl;

        // slim the THnSparse
        const Int_t axis_dimuon[4] = {0,1,2,3};
        const Int_t axis_singlemuon[5] = {0,1,2,3,4};
        // const Int_t axis_dimuon[3] = {0,1,2};
        // const Int_t axis_singlemuon[4] = {0,1,2,3};
        TList* generatedHistosLight = new TList();
        generatedHistosLight->SetName("GeneratedHistos_CAny");
        TList* reconstructedHistosLight = new TList();
        reconstructedHistosLight->SetName("ReconstructedHistos_CAny");

        //rescale generated histo
        // generatedHistos = (TList*)analysis->Get("GeneratedHistos_CAny");//STARLIGHT
        generatedHistos = (TList*)analysis->Get("SingleMuonHistos_CAny");//MCemb
        hDimuon = (THnSparse*)generatedHistos->FindObject("fHistoJPsiGenerated");
        hSingleMuon = (THnSparse*)generatedHistos->FindObject("fHistoSingleMuonGenerated");
        hDimuon->Scale(weight);
        hSingleMuon->Scale(weight);

        //slim generated histo and add to list 
        hDimuonGenerated = hDimuon->Projection(4, axis_dimuon,"e");
        hDimuonGenerated->SetName("fHistoJPsiGenerated");
        hSingleMuonGenerated = hSingleMuon->Projection(5, axis_singlemuon,"e");
        hSingleMuonGenerated->SetName("fHistoSingleMuonGenerated");
        generatedHistosLight->Add(hDimuonGenerated);
        generatedHistosLight->Add(hSingleMuonGenerated);
        

        //rescale reconstructed histo
        // reconstructedHistos = (TList*)analysis->Get("ReconstructedHistos_CAny");//STARLIGHT
        reconstructedHistos = (TList*)analysis->Get("DiMuonHistos_CAny");//MCemb
        hDimuon = (THnSparse*)reconstructedHistos->FindObject("fHistoDiMuonOS");
        hSingleMuon = (THnSparse*)reconstructedHistos->FindObject("fHistoSingleMuon");
        hDimuon->Scale(weight);
        hSingleMuon->Scale(weight);

        //slim generated histo and add to list 
        hDimuonReconstructed = hDimuon->Projection(4, axis_dimuon,"e");
        hDimuonReconstructed->SetName("fHistoDiMuonOS");
        hSingleReconstructed = hSingleMuon->Projection(5, axis_singlemuon,"e");
        hSingleReconstructed->SetName("fHistoSingleMuon");
        reconstructedHistosLight->Add(hDimuonReconstructed);
        reconstructedHistosLight->Add(hSingleReconstructed);
        
        analysis->ls();

        //analysis->Write("",TObject::kOverwrite);
        analysis->Close();
        TFile* savedFile = new TFile(Form("%s/runs/%i/AnalysisResultsWeighted.root", fileLocation.Data(), runList[i]), "RECREATE");
        generatedHistosLight->Write("GeneratedHistos_CAny",TObject::kSingleKey);
        reconstructedHistosLight->Write("ReconstructedHistos_CAny",TObject::kSingleKey);
        savedFile->Close();
        savedFile->~TFile();
        generatedHistosLight->~TList();
        reconstructedHistosLight->~TList();
    }

    printf(Form("Number of missing files %i", n));
    // generatedHistos->~TList();
    // reconstructedHistos->~TList();
    // eventHistos->~TList();
    // hDimuon->~THnSparse();
    // hSingleMuon->~THnSparse();
    // hEventsPerRun->~TH1I();
    
}
    