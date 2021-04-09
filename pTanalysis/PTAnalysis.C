
/*
 *  PTAnalysis.C
 *
 *  Created by Ophelie Bugnon on 13/05/19.
 *
 */
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
#include "TStyle.h"

const char* fileLocation ="$LOWPT/AnalysisResults_5TeV";
const char* fileNameData = "AnalysisResults_DataMerged.root";
const char* fileNameData2015 = "Individual_data_prod/AnalysisResults_15o_AOD229.root";
const char* fileNameSim = "AnalysisResults_cohJpsi_weighted.root";


void GetDimuonVSpT(const char* file=fileNameData, Double_t minPt=0, Double_t maxPt=2, Double_t minY=-4, Double_t maxY=-2.5,  Int_t minCent=70, Int_t maxCent=90)
{
    TCanvas* cNmumuVSpT = new TCanvas("cNmumuVSpT","");
	    gStyle->SetOptStat(0);
        gStyle->SetHistLineColor(4);
        //cNmumuVSpT->SetLogy();

    TFile* analysis = TFile::Open(Form("%s/%s", fileLocation,fileNameData));
        TList* diMuonHistos = (TList*)analysis->Get("DiMuonHistos_CMUL7");
        THnSparse* hOppositeSign = (THnSparse*)diMuonHistos->FindObject("fHistoDiMuonOS");

// invariant mass corresponding to the J/psi
        hOppositeSign->GetAxis(0)->SetRangeUser(2.8, 3.4);
        hOppositeSign->GetAxis(1)->SetRangeUser(minPt,maxPt);
        hOppositeSign->GetAxis(2)->SetRangeUser(minY,maxY);
        hOppositeSign->GetAxis(3)->SetRangeUser(minCent,maxCent);

    TH1* invMassDist = hOppositeSign->Projection(1);
        invMassDist->GetXaxis()->SetTitle("p_{T} (GeV/c)");
        invMassDist->GetYaxis()->SetTitle("d^{2}N/ dp_{T}dm (GeV^{-2}.c^{3})");
        invMassDist->SetTitle(" ");
        invMassDist->SetLineColor(46);
        invMassDist->SetMarkerStyle(20);
        invMassDist->SetMarkerColor(46);
        invMassDist->SetMarkerSize(0.7);        
        // invMassDist->Rebin(2);

        //Normalisation
        Double_t widthPt = invMassDist->GetXaxis()->GetBinWidth(1);
        invMassDist->Scale(widthPt*(3.4-2.8));

    invMassDist->Draw();

// invariant mass before to the J/psi
    THnSparse* hOppositeSignLower = (THnSparse*)diMuonHistos->FindObject("fHistoDiMuonOS");

        hOppositeSignLower->GetAxis(0)->SetRangeUser(1.5, 2.8);
        hOppositeSignLower->GetAxis(1)->SetRangeUser(minPt,maxPt);
        hOppositeSignLower->GetAxis(2)->SetRangeUser(minY,maxY);
        hOppositeSignLower->GetAxis(3)->SetRangeUser(minCent,maxCent);

    TH1* invMassDist2 = hOppositeSignLower->Projection(1);
        invMassDist2->SetTitle(" ");
        invMassDist2->SetLineColor(9);
        invMassDist2->SetMarkerStyle(21);
        invMassDist2->SetMarkerColor(9);        
        invMassDist2->SetMarkerSize(0.7);        
        // invMassDist2->Rebin(2);

        invMassDist2->Scale(widthPt*(2.8-1.5));

        Double_t nEntriesCoh= invMassDist->Integral(invMassDist->FindBin(0.65),invMassDist->FindBin(2));
        Double_t nEntriesInclusive= invMassDist2->Integral(invMassDist2->FindBin(0.65),invMassDist2->FindBin(2));
        Double_t scale = nEntriesCoh/nEntriesInclusive;
        printf("scale = %.2f\n", scale);
        invMassDist2->Scale(scale);
        invMassDist2->Draw("SAME");


    //legend
    TLegend* dlegend = new TLegend(0.5,0.66,0.88,0.88);
    dlegend->SetBorderSize(0);
    // dlegend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
    dlegend->AddEntry((TObject*)0, "Pb-Pb collisions #sqrt{s_{NN} }= 5.02 TeV", "");
    dlegend->AddEntry((TObject*)0,Form("%i-%i %%",minCent,maxCent),"");
    dlegend->AddEntry((TObject*)0,Form("%.2f < y < %.2f",minY,maxY),"");
    dlegend->AddEntry((TObject*)0,"","");

    dlegend->AddEntry(invMassDist,Form("#mu^{+}#mu^{-}, %.2f < m_{#mu#mu} < %.2f  GeV/c^{2}",2.8,3.4),"lep");
    dlegend->AddEntry(invMassDist2,Form("#mu^{+}#mu^{-} (#times %.2f), %.2f < m_{#mu#mu} < %.2f  GeV/c^{2}",scale, 1.5,2.8),"lep");

    dlegend->Draw();

}

void GetSTARLIGHTDimuon(Double_t minMass=2.8, Double_t maxMass=3.4, Double_t minPt=0., Double_t maxPt=2, Double_t minY=-4, Double_t maxY=-2.5)
{
    TCanvas* cNmumuVSpT = new TCanvas("cNmumuVSpT","");
	    gStyle->SetOptStat(0);
        //gStyle->SetHistLineColor(4);
        cNmumuVSpT->SetLogy();

        TFile* analysisCoh = TFile::Open(Form("%s/AnalysisResults_cohJpsiSTARLIGHT.root", fileLocation));
        TList* cohHistos = (TList*)analysisCoh->Get("ReconstructedHistos_CAny");
        THnSparse* hOScoh = (THnSparse*)cohHistos->FindObject("fHistoDiMuonOS");

        TFile* analysisIncoh = TFile::Open(Form("%s/AnalysisResults_incohJpsiSTARLIGHT.root", fileLocation));
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

void compareDataVsSim(const char* filedata=fileNameData, const char* filesim=fileNameSim, Double_t minMass=2.8, Double_t maxMass=3.4, Double_t minPt=0, Double_t maxPt=2, Double_t minY=-4, Double_t maxY=-2.5,  Int_t minCent=70, Int_t maxCent=90)
{
    TCanvas* cNmumuVSpT = new TCanvas("cNmumuVSpT","");
	    gStyle->SetOptStat(0);
        gStyle->SetHistLineColor(4);
        //cNmumuVSpT->SetLogy();


// Ful data sample
    TFile* analysis = TFile::Open(Form("%s/%s", fileLocation,fileNameData));
        TList* diMuonHistos = (TList*)analysis->Get("DiMuonHistos_CMUL7");
        THnSparse* hOppositeSign = (THnSparse*)diMuonHistos->FindObject("fHistoDiMuonOS");

        hOppositeSign->GetAxis(0)->SetRangeUser(minMass, maxMass);
        hOppositeSign->GetAxis(1)->SetRangeUser(minPt,maxPt);
        hOppositeSign->GetAxis(2)->SetRangeUser(minY,maxY);
        hOppositeSign->GetAxis(3)->SetRangeUser(minCent,maxCent);

    TH1* invMassDist = hOppositeSign->Projection(1);
        invMassDist->GetXaxis()->SetTitle("p_{T} (GeV/c)");
        invMassDist->GetYaxis()->SetTitle("d^2N/ dp_{T}dm (GeV^{-2}.c^3)");
        invMassDist->SetTitle(" ");
        invMassDist->Rebin(2);
        invMassDist->Draw();

        TFile* analysisCoh = TFile::Open(Form("%s/%s", fileLocation,fileNameSim));
        TList* cohHistos = (TList*)analysisCoh->Get("ReconstructedHistos_CAny");
        THnSparse* hOScoh = (THnSparse*)cohHistos->FindObject("fHistoDiMuonOS");

        hOScoh->GetAxis(0)->SetRangeUser(minMass, maxMass);
        hOScoh->GetAxis(1)->SetRangeUser(minPt,maxPt);
        hOScoh->GetAxis(2)->SetRangeUser(minY,maxY);
       
       
        TH1* invMassDistCoh = hOScoh->Projection(1);
        // invMassDistCoh->GetXaxis()->SetTitle("Dimuon p_{T} GeV/c");
        // invMassDistCoh->GetYaxis()->SetTitle("Counts per 100 MeV/c");
        // invMassDistCoh->SetTitle("OS dimuon distribution from STARLIGHT simulation");
    invMassDistCoh->Rebin(2);
        Double_t scale;
        Double_t nFirstBinData= invMassDist->GetBinContent(1);
        Double_t nFirstBinSim= invMassDistCoh->GetBinContent(1);
        scale = 3./4.*nFirstBinData/nFirstBinSim;
        cout << scale << " and " << nFirstBinSim << " = " << scale*nFirstBinSim << endl;
        invMassDistCoh->SetLineColor(kRed);
        invMassDistCoh->Scale(scale);
        invMassDistCoh->Draw("HIST SAME");



    //legend
    TLegend* dlegend = new TLegend(0.1,0.7,0.48,0.9);
    // dlegend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
    dlegend->AddEntry((TObject*)0, "Pb-Pb collisions #sqrt{s_{NN} }= 5.02 TeV", "");
    dlegend->AddEntry((TObject*)0,Form("%i-%i %%",minCent,maxCent),"");
    dlegend->AddEntry((TObject*)0,Form("%.2f < y < %.2f",minY,maxY),"");
    dlegend->AddEntry((TObject*)0,Form("%.2f < m_{#mu#mu} < %.2f  GeV/c^{2}",minMass,maxMass),"");
    dlegend->AddEntry((TObject*)0,"","");
    dlegend->AddEntry(invMassDist,"OS dimuon (data)","lep");
    dlegend->AddEntry(invMassDistCoh,"J/#psi coherent photoproduction (simulation)","l");
    dlegend->Draw();

}

void compare2DataSet(const char* filedata1=fileNameData, const char* filedata2=fileNameData2015, const char* filesim=fileNameSim, Double_t minMass=2.8, Double_t maxMass=3.4, Double_t minPt=0, Double_t maxPt=2, Double_t minY=-4, Double_t maxY=-2.5,  Int_t minCent=70, Int_t maxCent=90)
{
    TCanvas* cNmumuVSpT = new TCanvas("cNmumuVSpT","");
	    gStyle->SetOptStat(0);
        gStyle->SetHistLineColor(4);
        //cNmumuVSpT->SetLogy();

// Full data set
    TFile* analysis = TFile::Open(Form("%s/%s", fileLocation,fileNameData));
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
        invMassDist->SetMarkerColor(9);
        invMassDist->SetLineColor(9);
        invMassDist->Draw();

//2015 data only
        TFile* analysis2 = TFile::Open(Form("%s/%s", fileLocation,fileNameData2015));
        TList* diMuonHistos2 = (TList*)analysis2->Get("DiMuonHistos_CMUL7");
        THnSparse* hOppositeSign2 = (THnSparse*)diMuonHistos2->FindObject("fHistoDiMuonOS");

        hOppositeSign2->GetAxis(0)->SetRangeUser(minMass, maxMass);
        hOppositeSign2->GetAxis(1)->SetRangeUser(minPt,maxPt);
        hOppositeSign2->GetAxis(2)->SetRangeUser(minY,maxY);
        hOppositeSign2->GetAxis(3)->SetRangeUser(minCent,maxCent);

    TH1* invMassDist2 = hOppositeSign2->Projection(1);
        invMassDist2->GetXaxis()->SetTitle("Dimuon p_{T} GeV/c");
        invMassDist2->GetYaxis()->SetTitle("Counts per 100 MeV/c");
        invMassDist2->SetTitle(" ");
        invMassDist2->Rebin(2);
        invMassDist2->SetMarkerColor(46);
        invMassDist2->SetLineColor(46);
        invMassDist2->Draw("SAME");
 
 //Sim
        TFile* analysisCoh = TFile::Open(Form("%s/%s", fileLocation,fileNameSim));
        TList* cohHistos = (TList*)analysisCoh->Get("ReconstructedHistos_CAny");
        THnSparse* hOScoh = (THnSparse*)cohHistos->FindObject("fHistoDiMuonOS");

        hOScoh->GetAxis(0)->SetRangeUser(minMass, maxMass);
        hOScoh->GetAxis(1)->SetRangeUser(minPt,maxPt);
        hOScoh->GetAxis(2)->SetRangeUser(minY,maxY);
       
       
        TH1* invMassDistCoh = hOScoh->Projection(1);
        // invMassDistCoh->GetXaxis()->SetTitle("Dimuon p_{T} GeV/c");
        // invMassDistCoh->GetYaxis()->SetTitle("Counts per 100 MeV/c");
        // invMassDistCoh->SetTitle("OS dimuon distribution from STARLIGHT simulation");
    invMassDistCoh->Rebin(2);
        Double_t scale;
        Double_t nFirstBinData= invMassDist->GetBinContent(1);
        Double_t nFirstBinSim= invMassDistCoh->GetBinContent(1);
        scale = 0.9*nFirstBinData/nFirstBinSim;
        cout << scale << " and " << nFirstBinSim << " = " << scale*nFirstBinSim << endl;
        invMassDistCoh->SetLineColor(9);
        invMassDistCoh->Scale(scale);
        invMassDistCoh->Draw("HIST SAME");



    //legend
    TLegend* dlegend = new TLegend(0.1,0.7,0.48,0.9);
    // dlegend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
    dlegend->AddEntry((TObject*)0, "Pb-Pb collisions #sqrt{s_{NN} }= 5.02 TeV", "");
    dlegend->AddEntry((TObject*)0,Form("%i-%i %%",minCent,maxCent),"");
    dlegend->AddEntry((TObject*)0,Form("%.2f < y < %.2f",minY,maxY),"");
    dlegend->AddEntry((TObject*)0,Form("%.2f < m_{#mu#mu} < %.2f  GeV/c^{2}",minMass,maxMass),"");
    dlegend->AddEntry((TObject*)0,"","");
    dlegend->AddEntry(invMassDist,"OS dimuon (Full data set)","lep");
    dlegend->AddEntry(invMassDist2,"OS dimuon (2015 data set)","lep");
    dlegend->AddEntry(invMassDistCoh,"J/#psi coherent photoproduction (simulation for full data set)","l");
    dlegend->Draw();

}