#include <vector>
#include <iostream>
#include <algorithm>
#include <stdio.h>
#include "TVector.h"
#include "TF1.h"
#include "TString.h"
#include "TMath.h"
#include "TFile.h"

#include "FitFunctions.C"
//Initial parameters to fit background
	Double_t par_VWG[4]={1., 1.235, 1.600, 2.928};
	TString name_VWG[4]={"Norm_{VWG}", "#mu_{VWG}", "#alpha_{VWG}", "#beta_{VWG}"};

	// Double_t par_VWGQuadratic[5]={1., 2.15, 0.5, 0.14, 0.075}; //0-90
    // Double_t par_VWGQuadratic[5]={1., 0.59, 1.09, 0.023, 0.003}; //0-10 high pt
	//Double_t par_VWGQuadratic[5]={1., 1.235, 1.600, 2.928, -2.894}; //50-70// 30-50
	Double_t par_VWGQuadratic[5]={1., 2.3, 0.1, 1.82, -0.79};//70-90// 50-70 70-90
	//Double_t par_VWGQuadratic[5]={1., 2.3, 0.1, 1.82, -0.79};//70-90
	// Double_t par_VWGQuadratic[5]={1., 2.3, 0.5, 0.3, -0.79}; //0-10, 10-30, 30-50// 0-10, 10-30
	//Double_t par_VWGQuadratic[5]={1., 2., 0.4, 3, -2.}; // //80-90 
	//Double_t par_VWGQuadratic[5]={1., 1., 0.3, 0.5, -0.1}; //80-90 //
	TString name_VWGQuadratic[5]={"Norm_{VWG}", "#mu_{VWG}", "#alpha_{VWG}", "#beta_{VWG}", "#gamma_{VWG}"};

	
	// Double_t par_Pol2Pol3[6]={1., 1., 1., 1., 1., 1.};
    // Double_t par_Pol2Pol3[6]={1., -0.37, 0.035, 46.43, -32.94, 6.57}; //0-90
    // Double_t par_Pol2Pol3[6]={13., -0.42, 0.85, 20.50, -15.86, 3.3}; //70-90 large range
    // Double_t par_Pol2Pol3[6]={8, 2, 0, 22., -15., 3.}; //50-70 && 70-90 small rangZ
    // Double_t par_Pol2Pol3[6]={13., 5.01901, 0.106750, 22.2812, -15.6848, 3.13765}; //70-90 small range CB-coh CB-EMB
    // Double_t par_Pol2Pol3[6]={8.96996, -0.09, 0.666, 15.05, -11.52, 2.35}; //70-90 small range NA60
	// Double_t par_Pol2Pol3[6]={1., -0.35, 0.035, 6.5, -5.1, 1.04};
    Double_t par_Pol2Pol3[6]={1., 1.25, 0.45, 3.82, -3.11, 0.64};//70-90 starlight low pt
    // Double_t par_Pol2Pol3[6]={1., -0.28, 0.013, 21.4, -12.3, 3.};//10-30
	TString name_Pol2Pol3[6]={"Norm", "a1", "a2", "b1", "b2", "b3"};

//Initial parameters to fit signal
	Double_t par_jpsi[3]={1., 3.1, 0.065};
    Double_t par_psi[3]={1.7, par_jpsi[1]+3.686109-3.096916, par_jpsi[2]*1.05};

	TString name_CBjpsi[7]={"Norm_{J/#psi}", "#mu_{J/#psi}", "#sigma_{J/#psi}", "#alpha^{L}_{J/#psi}","n^{L}_{J/#psi}", "#alpha^{R}_{J/#psi}", "n^{R}_{J/#psi}"};
    TString name_CBpsi[3]={"Norm_{#psi'}", "#mu_{#psi'}", "#sigma_{#psi'}"};
	TString name_CBpsi_full[7]={"Norm_{#psi'}", "#mu_{#psi'}", "#sigma_{#psi'}", "#alpha^{L}_{#psi'}","n^{L}_{#psi'}", "#alpha^{R}_{#psi'}", "n^{R}_{#psi'}"};
	
	TString name_NA60jpsi[11]={"Norm_{J/#psi}", "#mu_{J/#psi}", "#sigma_{J/#psi}", "#alpha^{L}_{J/#psi}", "p1^{L}_{J/#psi}", "p2^{L}_{J/#psi}", "p3^{L}_{J/#psi}", "#alpha^{R}_{J/#psi}", "p1^{R}_{J/#psi}", "p2^{R}_{J/#psi}", "p3^{R}_{J/#psi}"};
    TString name_NA60psi[3]={"Norm_{#psi'}", "#mu_{#psi'}", "#sigma_{#psi'}"};
	TString name_NA60psi_full[11]={"Norm_{#psi'}", "#mu_{#psi'}", "#sigma_{#psi'}", "#alpha^{L}_{#psi'}", "p1^{L}_{#psi'}", "p2^{L}_{#psi'}", "p3^{L}_{#psi'}", "#alpha^{R}_{#psi'}", "p1^{R}_{#psi'}", "p2^{R}_{#psi'}", "p3^{R}_{#psi'}"};
	


//___________________________________________________________________________________________________________
//___________________________________________________________________________________________________________
//Get tails parameter from a file that contains vectors
std::vector<Double_t> GetTailsParameter(TString nameTails)
{
    TFile* analysis = TFile::Open("~/Documents/ALICE/AnalyseJPsi/macro/SignalExtraction/Tails_PbPb_5TeV.root");
    std::vector<Double_t> *vect;
    analysis->GetObject(Form("%s", nameTails.Data()), vect);
    std::vector<Double_t> param_tails = *vect;
    
    return param_tails;

}
//___________________________________________________________________________________________________________
//___________________________________________________________________________________________________________
TF1* BackGroundFunction(Efunction fName, Double_t xmin, Double_t xmax)
{
	Double_t* par;
	TString* name_par;
	TF1* BGFunction;
    Int_t nbPar = GetNPar(fName);
	
	switch (fName)
	{
	    case kVWG:
		    par = par_VWG;
		    name_par = name_VWG;
		    BGFunction =  new TF1("fitBG", VWG, xmin, xmax, GetNPar(kVWG));
		break;
	
	    case kVWGQuadratic:
		    par = par_VWGQuadratic;
		    name_par = name_VWGQuadratic;
		    BGFunction =  new TF1("fitBG", VWGQuadratic, xmin, xmax, GetNPar(kVWGQuadratic));
	    break;

        case kPol2Pol3:
            par = par_Pol2Pol3;
		    name_par = name_Pol2Pol3;
		    BGFunction =  new TF1("fitBG", Pol2Pol3, xmin, xmax, GetNPar(kPol2Pol3));
		break;    
	}

	for (int i=0; i<nbPar; i++)
	{
		BGFunction->SetParameter(i,par[i]);
		BGFunction->SetParName(i,name_par[i]);
	}
	return BGFunction;
}

//___________________________________________________________________________________________________________
//___________________________________________________________________________________________________________
TF1* SignalFunction(Efunction fSig, TString nTails, Epart fPart, Double_t xmin, Double_t xmax)
{
    Double_t* par_signal;
    TString* name_par_signal;
    std::vector<Double_t> par_tails = GetTailsParameter(nTails);
    

    Int_t nb_par_signal = GetNPar(fSig);

	TF1* fSignal;
    
    switch (fSig)
    {
        case kCBExtended:
            if (fPart == kJPsi) 
            {
                par_signal = par_jpsi;
                name_par_signal = name_CBjpsi;
            }
            else if (fPart == kPsi2S)
            {
                par_signal = par_psi;
                name_par_signal = name_CBpsi_full;
            }

            fSignal = new TF1("fitSignal", CrystalBallExtended, xmin, xmax, GetNPar(kCBExtended));
        break;
    
        case kNA60:
            if (fPart == kJPsi) 
            {
                par_signal = par_jpsi;
                name_par_signal = name_NA60jpsi;
            }
            else if (fPart == kPsi2S)
            {
                par_signal = par_psi;
                name_par_signal = name_NA60psi_full;
            }
    
            fSignal = new TF1("fitSignal", NA60, xmin, xmax, GetNPar(kNA60));

        break;
    }

	for (int i=0; i<nb_par_signal; i++)
	{
        //Norm, mean and sigma JPsi
        if (i < 3) fSignal->SetParameter(i, par_signal[i]);
        //Tails Jpsi
        else fSignal->FixParameter(i, par_tails[i-3]);
        
        fSignal->SetParName(i, name_par_signal[i]);
    }
    
	return fSignal;
}

//___________________________________________________________________________________________________________
//___________________________________________________________________________________________________________
TF1* DistributionFunction(Efunction fBG, Efunction fSig, TString nTails, Double_t xmin, Double_t xmax)
{
	Double_t* par_bg;
    Double_t* par_signal1;
    Double_t* par_signal2;
    std::vector<Double_t> par_tails = GetTailsParameter(nTails);

    TString* name_par_bg;
    TString* name_par_signal1;
    TString* name_par_signal2;

    Int_t nb_par_bg = GetNPar(fBG);
    Int_t nb_par_signal = GetNPar(fSig);
    Int_t nb_par = nb_par_bg + nb_par_signal + 1;

    Int_t fchoice = fBG + fSig;

	TF1* fDistrib;
    
	
	switch (fBG)
	{
	    case kVWGQuadratic:
		    par_bg = par_VWGQuadratic;
		    name_par_bg = name_VWGQuadratic;
		break;
	
	    case kPol2Pol3:
		    par_bg = par_Pol2Pol3;
		    name_par_bg = name_Pol2Pol3;
		break;
	}

    switch (fSig)
    {
        case kCBExtended:
            par_signal1 = par_jpsi;
            par_signal2 = par_psi;
            name_par_signal1 = name_CBjpsi;
            name_par_signal2 = name_CBpsi;
        break;
    
        case kNA60:
            par_signal1 = par_jpsi;
            par_signal2 = par_psi;
            name_par_signal1 = name_NA60jpsi;
            name_par_signal2 = name_NA60psi;
        break;
    }

    switch (fchoice)
    {
        case 220:
            fDistrib =  new TF1("fitDistrib", VWGquad_DoubleCBext, xmin, xmax, GetNPar(kVWGQuadratic)+GetNPar(kCBExtended)+1);
        break;

        case 230:
            fDistrib =  new TF1("fitDistrib", VWGquad_DoubleNA60, xmin, xmax, GetNPar(kVWGQuadratic)+GetNPar(kNA60)+1);
        break;

        case 320:
            fDistrib =  new TF1("fitDistrib", Pol2Pol3_DoubleCBext, xmin, xmax, GetNPar(kPol2Pol3)+GetNPar(kCBExtended)+1);
        break;

        case 330:
            fDistrib =  new TF1("fitDistrib", Pol2Pol3_DoubleNA60, xmin, xmax, GetNPar(kPol2Pol3)+GetNPar(kNA60)+1);
        break;
    }

    //Background
	for (int i=0; i<nb_par_bg; i++)
	{
        fDistrib->SetParameter(i,par_bg[i]);
		fDistrib->SetParName(i,name_par_bg[i]);
	}
    //Signal
    for (int j = 0; j < nb_par_signal; j++)
    {
        //Norm, mean and sigma JPsi
        if (j < 3) 
        {
            fDistrib->SetParameter(j+nb_par_bg, par_signal1[j]);
            if (j == 1 )   
            {
                fDistrib->SetParLimits(j+nb_par_bg,2.9,3.3);
            }
            if (j == 2) 
            {
                fDistrib->SetParLimits(j+nb_par_bg,0.04,0.1);  
            }
        }
        //Tails JPsi
        else fDistrib->FixParameter(j+nb_par_bg, par_tails[j-3]);
        fDistrib->SetParName(j+nb_par_bg, name_par_signal1[j]);
    }
    //psi2s
    fDistrib->SetParameter(nb_par_bg+nb_par_signal,par_signal2[0]);
    fDistrib->SetParName(nb_par_bg+nb_par_signal,name_par_signal2[0]);
	return fDistrib;
}