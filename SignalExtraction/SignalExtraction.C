#include <vector>
#include <iostream>
#include <algorithm>
#include <stdio.h>

#include <TStyle.h>
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
#include "TMatrixD.h"
#include "TFitResult.h"

#include "InitializeFunctions.C"
#include "SetRangeAndNameTest.C"
// #include "SetAndGetFitResults.C"

TH1* GetInvMassHisto(TFile* file, Double_t minMas, Double_t maxMass, Double_t minPt, Double_t maxPt, Double_t minY, Double_t maxY,  Int_t minCent, Int_t maxCent);
Bool_t isTailsParameter(TString nameTails);

#define JPSI_MASS   3.096916
#define PSI2S_MASS 3.686109

//___________________________________________________________________________________________________________
//___________________________________________________________________________________________________________
std::vector<Double_t> SignalExtraction(const char* file, Double_t minY, Double_t maxY, Double_t minPt, Double_t maxPt, Int_t minCent, Int_t maxCent, Double_t minMass, Double_t maxMass, Efunction fBackGround, Efunction fSignal, Etails pTails, Double_t minFit, Double_t maxFit, Bool_t isDisplayed, Bool_t isSaved)
{
    std::vector<Double_t> results;
    Int_t nb_bg = GetNPar(fBackGround);
    Int_t nb_sig = GetNPar(fSignal);

    //Get Histogram
        TFile* analysis = TFile::Open(Form("%s", file));
        TString nameTest = SetNameTest(fBackGround, fSignal, pTails, minFit, maxFit);
        TString rangeName = SetRangeName(minY, maxY, minPt, maxPt, minCent, maxCent);
        TString tailsName = SetTailsName(pTails, fSignal, fBackGround, minY, maxY, minPt, maxPt, minFit, maxFit);
        if(! isTailsParameter(tailsName))return results;

        TH1* invMassDist = GetInvMassHisto(analysis, minMass, maxMass, minPt, maxPt, minY, maxY, minCent, maxCent);
        invMassDist->SetTitle(Form("test %s in range %s", nameTest.Data(), rangeName.Data()));

	//Fit Background
        TF1* BGFunction =  BackGroundFunction(fBackGround, 0., 10.);
        
        reject = kTRUE;
        Double_t n_BG = invMassDist->GetBinContent(invMassDist->GetXaxis()->FindBin(minFit));
        BGFunction->SetParameter(0,n_BG);
        
        int secur1 = 0;
        TFitResultPtr fitStatus1;
        Double_t chi21 = 3;
        int covMatrixStatus1 = 2;

        do
        {
            fitStatus1 = invMassDist->Fit("fitBG","LSN","", minFit, maxFit);
            chi21 = BGFunction->GetChisquare()/BGFunction->GetNDF();
            covMatrixStatus1 = fitStatus1->CovMatrixStatus();
              
            secur1++;
            if (secur1 > 3)  break;
            cout << "FitStatus = " << fitStatus1 << "     chi2/NDF = " << chi21 << "     cov matrix status = " << covMatrixStatus1 << endl;
        }
        // while (fitStatus1 !=0 ||  chi21 > 2.5  || covMatrixStatus1 != 3);
        while (fitStatus1 !=0 || covMatrixStatus1 != 3);

        reject = kFALSE;


	//Fit Invariant mass distribution
        TF1* fitFunction =  DistributionFunction(fBackGround, fSignal, tailsName, 0., 10.);
        Double_t n_jpsi = invMassDist->GetBinContent(invMassDist->GetXaxis()->FindBin(JPSI_MASS));
        n_jpsi -= BGFunction->Eval(JPSI_MASS);

        for(int i = 0; i < nb_bg; i++)
        {
            fitFunction->SetParameter(i,BGFunction->GetParameter(i));
        }  
        fitFunction->SetParameter(nb_bg,n_jpsi);

        int secur2 = 0;
        TFitResultPtr fitStatus;
        Double_t chi22 = 3;
        int covMatrixStatus2 = 2;

        do
        {
            fitStatus = invMassDist->Fit("fitDistrib","LSN","",minFit,maxFit);
            chi22 = fitFunction->GetChisquare()/fitFunction->GetNDF();
            covMatrixStatus2 = fitStatus->CovMatrixStatus();
            cout << "FitStatus = " << fitStatus << "     chi2/NDF = " << chi22 << "     cov matrix status = " << covMatrixStatus2 << endl;

            secur2++;
            if (secur2 > 2){
            cout << "______________________________" << endl;
            cout << " The fit has not converged "    << endl;    
            cout << "______________________________" << endl;
            break;
            }
        }
        while (fitStatus !=0 ||  chi22 > 2.5 || covMatrixStatus2 != 3);
        
	//Draw Function
	    TF1 *JPsiFunction = SignalFunction(fSignal, tailsName, kJPsi, 0., 10.);
        TF1 *Psi2SFunction = SignalFunction(fSignal, tailsName, kPsi2S, 0., 10.);
        Double_t para[nb_bg+nb_sig+1];
       
        fitFunction->GetParameters(para);
	    BGFunction->SetParameters(para);
	    JPsiFunction->SetParameters(&(para[nb_bg]));
        Psi2SFunction->SetParameters(&(para[nb_bg]));
        Psi2SFunction->SetParameter(0, para[nb_bg+nb_sig]);
        Psi2SFunction->SetParameter(1, para[nb_bg+1]+3.686109 - 3.096916);
        Psi2SFunction->SetParameter(2, para[nb_bg+2]*1.05);

        Double_t paraJPsi[nb_sig];
        JPsiFunction->GetParameters(paraJPsi);
        Double_t paraPsi2S[nb_sig];
        Psi2SFunction->GetParameters(paraPsi2S);

    //Compute Number of JPsi
	    double dx = invMassDist->GetBinWidth(1);
	    Double_t N_JPsi = JPsiFunction->Integral(0,10)/dx;
        Double_t mean_JPsi = fitFunction->GetParameter(nb_bg+1);
        Double_t Errmean_JPsi = fitFunction->GetParError(nb_bg+1);
        Double_t sigma_JPsi = fitFunction->GetParameter(nb_bg+2);
        Double_t Errsigma_JPsi = fitFunction->GetParError(nb_bg+2);

	    Double_t covmat[nb_sig][nb_sig];
	    for (int k=0;k<nb_sig;k++)
	    {
		    for( int t=0;t<nb_sig;t++)
		    { 
			    covmat[k][t]=(fitStatus->GetCovarianceMatrix())(k+nb_bg,t+nb_bg);
            }    
	    }

	    Double_t Err_JPsi = JPsiFunction->IntegralError(0,10,&paraJPsi[0],&covmat[0][0], 0.05)/(dx);
	    cout << "Number of J/psi " << N_JPsi << " with associated error " << Err_JPsi << endl;
        
    //Compute Number of Psi2S
	    Double_t N_Psi2S = Psi2SFunction->Integral(0,9)/dx;

	    Double_t covmatPsi2S[nb_sig][nb_sig];
        int m=0;
        int n=0;
	    for (int k=0;k<nb_sig;k++)
	    {
		    for( int t=0;t<nb_sig;t++)
		    { 
            if(k==0) m=nb_sig;
            else m=k;
            if(t==0) n=nb_sig;
            else n=t;
               
			    covmatPsi2S[k][t]=(fitStatus->GetCovarianceMatrix())(m+nb_bg,n+nb_bg);
		    }
	    }
        
	    Double_t Err_Psi2S = Psi2SFunction->IntegralError(0,9.,&paraPsi2S[0],&covmatPsi2S[0][0], 0.05)/(dx);
	    cout << "Number of Psi2S " << N_Psi2S << " with associated error " << Err_Psi2S << endl;

    //Set Results
        results.push_back(N_JPsi);
        results.push_back(Err_JPsi);
        results.push_back(N_Psi2S);
        results.push_back(Err_Psi2S);
        results.push_back(mean_JPsi);
        results.push_back(Errmean_JPsi);
        results.push_back(sigma_JPsi);
        results.push_back(Errsigma_JPsi);
        results.push_back(chi22);
        results.push_back(0.);
        results.push_back(fitStatus);
        results.push_back(covMatrixStatus2);

        
    //Display 
        if(isDisplayed==kTRUE)
        {
            TCanvas* cInvMassDist = new TCanvas(Form("%s_%s",rangeName.Data(), nameTest.Data()), Form("Range %s for test %s",rangeName.Data(), nameTest.Data()));
            //gStyle->SetOptFit(1111);
            gStyle->SetOptFit(0);
            gStyle->SetOptStat(0);
            gStyle->SetHistLineColor(1);

            invMassDist->Draw();
            fitFunction->SetLineColor(kBlue);
            fitFunction->Draw("SAME");
            BGFunction->SetLineColor(kBlue);
            BGFunction->SetLineStyle(7);
            BGFunction->Draw("SAME");
            JPsiFunction->Draw("SAME");
            Psi2SFunction->SetLineColor(8);
            Psi2SFunction->Draw("SAME");

        //Set Legend
            TLatex* text = new TLatex();
            text->SetNDC();
            text->SetTextAlign(12);
            text->SetTextFont(43);
            text->SetTextSize(15);

            // text->DrawLatex(0.62, 0.80, Form("%i-%i %%",minCent, maxCent));
            // text->DrawLatex(0.62, 0.76, Form("%.2f < y < %.2f",minY, maxY));
            // text->DrawLatex(0.62, 0.72, Form("p_{T} < %.2f GeV/c", maxPt));
            //text->DrawLatex(0.62, 0.8, Form("N_{J/#psi} = %i #pm %i",TMath::Nint(N_JPsi),TMath::Nint(Err_JPsi)));

            
            text->DrawLatex(0.65, 0.88, Form("N_{J/#psi} = %i #pm %i",TMath::Nint(N_JPsi),TMath::Nint(Err_JPsi)));
            for (int np = 0; np < nb_sig + 1; np++)
            {
               text->DrawLatex(0.65, 0.84-np*0.04, Form("%s = %.3f #pm %.3f", fitFunction->GetParName(np+nb_bg), fitFunction->GetParameter(np+nb_bg), fitFunction->GetParError(np+nb_bg)));
            }
            text->DrawLatex(0.65, 0.84-(nb_bg + nb_sig-2)*0.04, Form("#chi^{2}/NDF = %.3f", chi22) );
            text->DrawLatex(0.65, 0.84-(nb_bg + nb_sig-1)*0.04, Form("Fit status = %i", (int)fitStatus) );
            text->DrawLatex(0.65, 0.84-(nb_bg + nb_sig)*0.04, Form("Cov matrix status = %i", covMatrixStatus2) );

            if(isSaved) cInvMassDist->SaveAs(".pdf");
        }    
    analysis->Close();
    return results;

}

//___________________________________________________________________________________________________________
//___________________________________________________________________________________________________________
TH1* GetInvMassHisto(TFile* file, Double_t minMas, Double_t maxMass, Double_t minPt, Double_t maxPt, Double_t minY, Double_t maxY,  Int_t minCent, Int_t maxCent)
{
    // TFile* analysis = TFile::Open(Form("%s", file));
    TList* diMuonHistos = (TList*)file->Get("DiMuonHistos_CMUL7");
    THnSparse* hOppositeSign = (THnSparse*)diMuonHistos->FindObject("fHistoDiMuonOS");

    hOppositeSign->GetAxis(0)->SetRangeUser(minMas, maxMass);
    hOppositeSign->GetAxis(1)->SetRangeUser(minPt,maxPt);
    hOppositeSign->GetAxis(2)->SetRangeUser(minY,maxY);
    hOppositeSign->GetAxis(3)->SetRangeUser(minCent,maxCent);

    TH1* invMassDist = hOppositeSign->Projection(0,"e");
    invMassDist->GetXaxis()->SetTitle("m_{#mu#mu} GeV/c^{2}");
    invMassDist->GetYaxis()->SetTitle("Counts per 20 MeV/c^{2}");
    invMassDist->SetTitle(" ");

    return invMassDist;
}
//___________________________________________________________________________________________________________
//___________________________________________________________________________________________________________
Bool_t isTailsParameter(TString nameTails)
{
    TFile* analysis = TFile::Open("~/Documents/ALICE/AnalyseJPsi/macro/SignalExtraction/Tails_PbPb_5TeV.root");
    std::vector<Double_t> *vect;
    analysis->GetObject(Form("%s", nameTails.Data()), vect);
    if(vect) return kTRUE;
    else return kFALSE;

}