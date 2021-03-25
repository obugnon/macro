/*
 *  TailsExtraction.C
 *
 *  Created by Ophelie Bugnon on 12/12/19.
 *
 */
#include <vector>
#include <iostream>
#include <algorithm>
#include <stdio.h>
#include <stdarg.h>

#include "TFile.h"
#include "TList.h"
#include "THnSparse.h"
#include "TH1.h"
#include "TLine.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TString.h"
#include "TMath.h"
#include "TLatex.h"
#include "TVector.h"
#include "TFitResult.h"
#include "TLegend.h"

#include "InitializationTailsExtraction.C"
#include "SetRangeAndNameTest.C"


#define JPSI_MASS   3.096916
#define PSI2S_MASS 3.686109

TH1* GetInvMass(const char* file, Etails pTails, Double_t minMas, Double_t maxMass, Double_t minPt, Double_t maxPt, Double_t minY, Double_t maxY);

// const char* fileName = "AnalysisResults_embedding_weighted.root";
//const char* fileName = "AnalysisResults_pp.root";
const char* fileName = "AnalysisResults_cohJpsi_weighted.root";
// const char* fileName = "AnalysisResults_incohJpsi_weighted.root";
const char* fileLocation = "~/Documents/ALICE/AnalyseJPsi/AnalysisResults_5TeV";


Double_t min_y = -4.;
Double_t max_y = -2.5;
Double_t min_pt = 0.;
Double_t max_pt = 0.3;
Double_t minMass = 1.;
Double_t maxMass = 5.;
Double_t min_fit = 1;
Double_t max_fit = 5;
Efunction f_BackGround = kVWGQuadratic;
Efunction f_Signal = kCBExtended;
Etails p_tails = kSTARLIGHTcoh;

//___________________________________________________________________________________________________________
//___________________________________________________________________________________________________________
void TailsExtraction_signal(const char* file=fileName, Double_t minY=min_y, Double_t maxY=max_y, Double_t minPt=min_pt, Double_t maxPt=max_pt, Efunction fSignal=f_Signal, Etails pTails=p_tails, Double_t minFit=min_fit, Double_t maxFit=max_fit, Bool_t isSaved = kFALSE, Bool_t isComparedToIntegratedJPsi=kFALSE)
{
    std::vector<Double_t> results;
    //TVector results;
    Int_t nb_sig = GetNPar(fSignal);
    TString nameFunction;
    if (fSignal==kCBExtended) nameFunction = "Crystal Ball Extended";
    else if (fSignal==kNA60) nameFunction = "NA60";

    //Get Histogram
        TString nameTails = SetTailsName(pTails, fSignal, f_BackGround, minY, maxY, minPt, maxPt, minFit, maxFit);
        TH1* invMassDist = GetInvMass(Form("%s/%s",fileLocation, file), pTails, minMass, maxMass, minPt, maxPt, minY, maxY);
        invMassDist->SetTitle(Form("%s in %.1f-%.1f GeV/c", nameFunction.Data(), minPt, maxPt));

        TCanvas* cInvMassDist = new TCanvas(Form("c_%s", nameTails.Data()), Form("Invariant Mass distribution"));
        gStyle->SetOptFit(1111);
        // gStyle->SetOptFit(0);
        // gStyle->SetOptStat(0);
        gStyle->SetHistLineColor(1);

        invMassDist->Draw();

	//Fit Invariant mass distribution
        TF1 *signalFunction = SignalFunction(fSignal,pTails, kJPsi, minFit, maxFit, kFALSE,kFALSE);
        Double_t n_jpsi = invMassDist->GetBinContent(invMassDist->GetXaxis()->FindBin(JPSI_MASS));

        signalFunction->SetParameter(0,n_jpsi);

        int secur = 0;
        TFitResultPtr fitStatus;
        Double_t chi2 = 3;
        int covMatrixStatus = 2;

        do
        {
            fitStatus = invMassDist->Fit("fitSignal","LS","",minFit,maxFit);
            chi2 = signalFunction->GetChisquare()/signalFunction->GetNDF();
            covMatrixStatus = fitStatus->CovMatrixStatus();

            secur++;
            if (secur > 15){
                cout << "______________________________" << endl;
                cout << " The fit has not converged "    << endl;    
                cout << "______________________________" << endl;
                break;
            }
            cout << "FitStatus = " << fitStatus << "     chi2/NDF = " << chi2 << "     cov matrix status = " << covMatrixStatus << endl;
        }
        // while (fitStatus !=0 ||  chi2 > 2.5 || covMatrixStatus != 3);
        while (fitStatus !=0 || covMatrixStatus != 3);

    
	//Set Results
        Double_t para[nb_sig];
        signalFunction->GetParameters(para);
        for(int i=3; i<nb_sig; i++) {
            cout << para[i] << endl;
            results.push_back(para[i]);
        }

    //Save tails
    if (isSaved == kTRUE)
    {
        TFile *fout = TFile::Open("Tails_PbPb_5TeV_weighted.root", "UPDATE");
        fout->WriteObject(&results, Form("%s", nameTails.Data()));
        fout->Close();
    }

    if(isComparedToIntegratedJPsi)
    {
        //Compute Number of JPsi
	    double dx = invMassDist->GetBinWidth(1);
	    Double_t N_JPsi = signalFunction->Integral(1,5)/dx;
        Double_t Err_JPsi_intHisto;
        Double_t N_JPsi_intHisto = invMassDist->IntegralAndError(invMassDist->GetXaxis()->FindBin(1.), invMassDist->GetXaxis()->FindBin(5), Err_JPsi_intHisto);

	    Double_t covmat[nb_sig][nb_sig];
	    for (int k=0;k<nb_sig;k++)
	    {
		    for( int t=0;t<nb_sig;t++)
		    { 
			    covmat[k][t]=(fitStatus->GetCovarianceMatrix())(k,t);
            }    
	    }

	    Double_t Err_JPsi = signalFunction->IntegralError(1,5,&para[0],&covmat[0][0], 0.05)/(dx);
	    cout << "Number of J/psi " << N_JPsi << " with associated error " << Err_JPsi << " and number of J/psi in histo " << N_JPsi_intHisto << " with associated error " << Err_JPsi_intHisto <<endl;
        TLatex* text = new TLatex();
            text->SetNDC();
            text->SetTextAlign(12);
            text->SetTextFont(43);
            text->SetTextSize(15);

            text->DrawLatex(0.62, 0.85, Form("Integral of the functionnal form NJpsi = %.f #pm %.f", N_JPsi, Err_JPsi));
            text->DrawLatex(0.62, 0.80, Form("Integral of the histogram NJpsi = %.f #pm %.f", N_JPsi_intHisto, Err_JPsi_intHisto));

    }
    return;    
}

// Set manually a set of tails
//___________________________________________________________________________________________________________
//___________________________________________________________________________________________________________
void AddTailsSet(Efunction fSignal, Etails pTails, Efunction fBackground, Double_t minY, Double_t maxY, Double_t minPt, Double_t maxPt, Double_t minFit, Double_t maxFit)
{
    std::vector<Double_t> results;
    results.push_back(0.891);
    results.push_back(8.345);
    results.push_back(1.833);
    results.push_back(20.105);

    TString nameTails = SetTailsName(pTails, fSignal, fBackground, minY, maxY, minPt, maxPt, minFit, maxFit);
    TFile *fout = TFile::Open("Tails_PbPb_5TeV_weighted.root", "UPDATE");
    fout->WriteObject(&results, Form("%s", nameTails.Data()));
    fout->Close();
}

//___________________________________________________________________________________________________________
//___________________________________________________________________________________________________________
TH1* GetInvMass(const char* file, Etails pTails, Double_t minMass, Double_t maxMass, Double_t minPt, Double_t maxPt, Double_t minY, Double_t maxY)
{
    TFile* analysis = TFile::Open(Form("%s", file));
    if(!analysis) {
            Error("GetInvMass",Form("Cannot open Analysis File %s", file));
            return results;
        }
    TH1* invMassDist;
    if(pTails == kPP)
    {
        THnSparse* hOppositeSign = static_cast<THnSparse*>(analysis->FindObjectAny("CMUL_hDimuon"));

        hOppositeSign->GetAxis(2)->SetRangeUser(minMass, maxMass);
        hOppositeSign->GetAxis(0)->SetRangeUser(minPt,maxPt);
        hOppositeSign->GetAxis(1)->SetRangeUser(minY,maxY);
        invMassDist = hOppositeSign->Projection(2,"e");
   }
    else if(pTails == kSTARLIGHTcoh || pTails == kSTARLIGHTincoh || pTails == kEMB)
    {
        TList* diMuonHistos = (TList*)analysis->Get("ReconstructedHistos_CAny");
        THnSparse* hOppositeSign = (THnSparse*)diMuonHistos->FindObject("fHistoDiMuonOS");

        hOppositeSign->GetAxis(0)->SetRangeUser(minMass, maxMass);
        hOppositeSign->GetAxis(1)->SetRangeUser(minPt,maxPt);
        hOppositeSign->GetAxis(2)->SetRangeUser(minY,maxY);
        invMassDist = hOppositeSign->Projection(0,"e");
    } 
    invMassDist->GetXaxis()->SetTitle("m_{#mu#mu} GeV/c^{2}");
    invMassDist->GetYaxis()->SetTitle("Counts per 20 MeV/c^{2}");
    invMassDist->SetTitle(" ");
    return invMassDist;
}
