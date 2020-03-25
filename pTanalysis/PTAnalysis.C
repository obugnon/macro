
#include <vector>
#include <iostream>
#include <algorithm>
#include <stdio.h>

#include "TFile.h"
#include "TList.h"
#include "TLegend.h"
#include "THnSparse.h"
#include "TH1.h"
#include "TLine.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TString.h"
#include "TMath.h"
#include "TLatex.h"


//TH1* GetDimuonVSpT(const char* file, Double_t minMas, Double_t maxMass, Double_t minPt, Double_t maxPt, Double_t minY, Double_t maxY,  Int_t minCent, Int_t maxCent)
//void GetDimuonVSpT(const char* file, Double_t minMass = 2.85, Double_t maxMass = 3.35, Double_t minPt = 0., Double_t maxPt = 2.7, Double_t minY = -4, Double_t maxY = -2.5,  Int_t minCent = 70, Int_t maxCent=90)
void GetDimuonVSpT(const char* file, Double_t minMass, Double_t maxMass, Double_t minPt, Double_t maxPt, Double_t minY, Double_t maxY,  Int_t minCent, Int_t maxCent)
{
    TCanvas* cNmumuVSpT = new TCanvas("cNmumuVSpT","");
	    gStyle->SetOptStat(0);
        gStyle->SetHistLineColor(4);
        //cNmumuVSpT->SetLogy();

    TFile* analysis = TFile::Open(Form("%s", file));
        TList* diMuonHistos = (TList*)analysis->Get("DiMuonHistos_CMUL7");
        THnSparse* hOppositeSign = (THnSparse*)diMuonHistos->FindObject("fHistoDiMuonOS");

        hOppositeSign->GetAxis(0)->SetRangeUser(minMass, maxMass);
        hOppositeSign->GetAxis(1)->SetRangeUser(minPt,maxPt);
        hOppositeSign->GetAxis(2)->SetRangeUser(minY,maxY);
        hOppositeSign->GetAxis(3)->SetRangeUser(minCent,maxCent);

    TH1* invMassDist = hOppositeSign->Projection(1);
        invMassDist->GetXaxis()->SetTitle("Dimuon p_{T} GeV/c");
        invMassDist->GetYaxis()->SetTitle("Counts per 100 MeV/c");
        invMassDist->SetTitle(" ");
        invMassDist->Rebin(2);

    invMassDist->Draw();

    //Set Legend
        TLatex* text = new TLatex();
        text->SetNDC();
        text->SetTextAlign(13);
	    text->SetTextFont(43);
	    text->SetTextSize(20);
        text->DrawLatex(.60, .85, "Pb-Pb #sqrt{s_{NN}} = 5.02 TeV");
        text->DrawLatex(.60, .80, Form("%i-%i %%",minCent,maxCent));
        text->DrawLatex(.60, .75, Form("%.2f < y < %.2f",minY,maxY));
        text->DrawLatex(.60, .70, Form("%.2f < m_{#mu#mu} < %.2f  GeV/c^{2}",minMass,maxMass));
          
    //return invMassDist;
}

void GetSTARLIGHTDimuon(Double_t minMass=1, Double_t maxMass=8, Double_t minPt=0., Double_t maxPt=8, Double_t minY=-4, Double_t maxY=-2.5)
{
    TCanvas* cNmumuVSpT = new TCanvas("cNmumuVSpT","");
	    gStyle->SetOptStat(0);
        //gStyle->SetHistLineColor(4);
        cNmumuVSpT->SetLogy();

        TFile* analysisCoh = TFile::Open("~/Documents/ALICE/AnalyseJPsi/AnalysisResults_5TeV/AnalysisResults_cohJpsiSTARLIGHT.root");
        TList* cohHistos = (TList*)analysisCoh->Get("ReconstructedHistos_CAny");
        THnSparse* hOScoh = (THnSparse*)cohHistos->FindObject("fHistoDiMuonOS");

        TFile* analysisIncoh = TFile::Open("~/Documents/ALICE/AnalyseJPsi/AnalysisResults_5TeV/AnalysisResults_incohJpsiSTARLIGHT.root");
        TList* incohHistos = (TList*)analysisIncoh->Get("ReconstructedHistos_CAny");
        THnSparse* hOSincoh = (THnSparse*)incohHistos->FindObject("fHistoDiMuonOS");

        hOScoh->GetAxis(0)->SetRangeUser(minMass, maxMass);
        hOScoh->GetAxis(1)->SetRangeUser(minPt,maxPt);
        hOScoh->GetAxis(2)->SetRangeUser(minY,maxY);
        

        hOSincoh->GetAxis(0)->SetRangeUser(minMass, maxMass);
        hOSincoh->GetAxis(1)->SetRangeUser(minPt,maxPt);
        hOSincoh->GetAxis(2)->SetRangeUser(minY,maxY);
       
        TH1* invMassDistCoh = hOScoh->Projection(1);
        invMassDistCoh->GetXaxis()->SetTitle("Dimuon p_{T} GeV/c");
        invMassDistCoh->GetYaxis()->SetTitle("Counts per 100 MeV/c");
        invMassDistCoh->SetTitle("OS dimuon distribution from STARLIGHT simulation");
        invMassDistCoh->SetLineColor(46);
        invMassDistCoh->Rebin(2);
        invMassDistCoh->Draw();

        TH1* invMassDistIncoh = hOSincoh->Projection(1);
        invMassDistIncoh->Rebin(2);
        invMassDistIncoh->SetLineColor(9);
        invMassDistIncoh->Draw("SAME");

    //Set Legend
       TLegend* legend = new TLegend(0.1,0.7,0.48,0.9);
        // legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
        legend->AddEntry(invMassDistCoh,"J/#psi coherent photoproduction","l");
        legend->AddEntry(invMassDistIncoh,"J/#psi incoherent photoproduction","l");
        
        legend->Draw();

        // TLatex* text = new TLatex();
        // text->SetNDC();
        // text->SetTextAlign(13);
	    // text->SetTextFont(43);
	    // text->SetTextSize(20);
        // text->DrawLatex(.60, .85, "Pb-Pb #sqrt{s_{NN}} = 5.02 TeV");
        // text->DrawLatex(.60, .75, Form("%.2f < y < %.2f",minY,maxY));
        // text->DrawLatex(.60, .70, Form("%.2f < m_{#mu#mu} < %.2f  GeV/c^{2}",minMass,maxMass));
          
    //return invMassDist;
}