/*
 *  InitializeFunctions.C
 *
 *  Created by Ophelie Bugnon on 06/12/19.
 *
 */

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

    Double_t par_VWGQuadratic[5]={1., 2.3, 0.5, 0.3, -0.79}; //low pT : 0-10, 10-30, 30-50, 0-90
    // Double_t par_VWGQuadratic[5]={1., 2.3, 0., 1.82, -0.79};//low pT : 50-70 70-90
    // Double_t par_VWGQuadratic[5]={1., 1.9, 0.47, 0.47, 0.096};//high pT : 50-70 70-90

    // Double_t par_VWGQuadratic[5]={1., 2.2, 0.22, 1.28, -3.8};//50-70 70-90
	// Double_t par_VWGQuadratic[5]={1., 2.15, 0.5, 0.14, 0.075}; //0-90
	// Double_t par_VWGQuadratic[5]={1., 2.3, 0.1, 1.82, -0.79};//70-90
	// Double_t par_VWGQuadratic[5]={1., 2., 0.4, 3, -2.}; // //80-90 
	// Double_t par_VWGQuadratic[5]={1., 1., 0.3, 0.5, -0.1}; //Very high pT : 30-50, 70-90
	// Double_t par_VWGQuadratic[5]={1., 2.3, 0.032, 1.9, -0.1}; //Very high pT : mid cent
    // Double_t par_VWGQuadratic[5]={1., 2.9, 0.97, 0.45, 0.058}; //Very high pT : 50-70


	TString name_VWGQuadratic[5]={"Norm_{VWG}", "#mu_{VWG}", "#alpha_{VWG}", "#beta_{VWG}", "#gamma_{VWG}"};

	
    Double_t par_Pol2Pol3[6]={1., -0.35, 0.035, 6.5, -5.1, 1.04}; //low pT : 0-10, 10-30, 30-50, 0-90
    // Double_t par_Pol2Pol3[6]={1., 1.25, 0.45, 3.82, -3.11, 0.64};//low pT : 50-70, 70-90
    // Double_t par_Pol2Pol3[6]={1., 1050, -250, 116, -475, 232};//high pT : 50-70, 70-90
    // Double_t par_Pol2Pol3[6]={8.96996, -0.09, 0.666, 15.05, -11.52, 2.35}; //low pT : 70-90 
    // Double_t par_Pol2Pol3[6]={8, 2, 0, 22., -15., 3.}; //high pT:  70-90 small range


    // Double_t par_Pol2Pol3[6]={13., 5.01901, 0.106750, 22.2812, -15.6848, 3.13765}; //low pT : 70-90 Starlight
    // Double_t par_Pol2Pol3[6]={1., -0.37, 0.035, 46.43, -32.94, 6.57}; //0-90
    // Double_t par_Pol2Pol3[6]={1., 100, -19, 430, 100, -27}; //0-90
    // Double_t par_Pol2Pol3[6]={13., -0.42, 0.85, 20.50, -15.86, 3.3}; // Very high pT : 10-30, 70-90 
    // Double_t par_Pol2Pol3[6]={13., -0.68, 0.813, 1.532, -1.27, 0.25}; // Very high pT : 50-70 

    TString name_Pol2Pol3[6]={"Norm", "a1", "a2", "b1", "b2", "b3"};

    Double_t par_DoubleExp[4]={1265, -1.4, -0.04, -0.04}; 
    TString name_DoubleExp[4]={"Norm1_{DoubleExp}", "#alpha1__{DoubleExp}", "Norm2_{DoubleExp}", "#alpha2__{DoubleExp}"};

    Double_t par_Exp[2]={1, 1}; 
    TString name_Exp[2]={"Norm1_{DoubleExp}", "#alpha1__{DoubleExp}"};


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
    TFile* analysis = TFile::Open("$LOWPT/macro/ResultFiles/Tails_PbPb_5TeV.root");
    std::vector<Double_t> param_tails;
    if(!analysis)
    {
        Error("InitializeFunctions","Cannot open Analysis File Tails_PbPb_5TeV.root");
        return param_tails;
    }
    std::vector<Double_t> *vect;
    analysis->GetObject(Form("%s", nameTails.Data()), vect);
    param_tails = *vect;
    
    return param_tails;

}
//___________________________________________________________________________________________________________
//___________________________________________________________________________________________________________
TF1* BackGroundFunction(Efunction fName, Double_t xmin, Double_t xmax)
{
	Double_t* par=nullptr;
	TString* name_par=nullptr;
	TF1* BGFunction=nullptr;
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

        case kDoubleExp:
            par = par_DoubleExp;
            name_par = name_DoubleExp;
            BGFunction = new TF1("fitBG", DoubleExp, xmin, xmax, GetNPar(kDoubleExp));
        break;

        case kExp:
            par = par_Exp;
            name_par = name_Exp;
            BGFunction = new TF1("fitBG", Exp, xmin, xmax, GetNPar(kExp));
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
    Double_t* par_signal=nullptr;
    TString* name_par_signal=nullptr;
    std::vector<Double_t> par_tails = GetTailsParameter(nTails);
    

    Int_t nb_par_signal = GetNPar(fSig);

	TF1* fSignal=nullptr;
    
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
	Double_t* par_bg=nullptr;
    Double_t* par_signal1=nullptr;
    Double_t* par_signal2=nullptr;
    std::vector<Double_t> par_tails = GetTailsParameter(nTails);

    TString* name_par_bg=nullptr;
    TString* name_par_signal1=nullptr;
    TString* name_par_signal2=nullptr;

    Int_t nb_par_bg = GetNPar(fBG);
    Int_t nb_par_signal = GetNPar(fSig);
    Int_t nb_par = nb_par_bg + nb_par_signal + 1;

    Int_t fchoice = fBG + fSig;

	TF1* fDistrib=nullptr;
    
	
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

        case kDoubleExp:
            par_bg = par_DoubleExp;
            name_par_bg = name_DoubleExp;
        break;

        case kExp:
            par_bg = par_Exp;
            name_par_bg = name_Exp;
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

        case 420:
            fDistrib =  new TF1("fitDistrib", DoubleExp_DoubleCBext, xmin, xmax, GetNPar(kDoubleExp)+GetNPar(kCBExtended)+1);
        break;

        case 430:
            fDistrib =  new TF1("fitDistrib", DoubleExp_DoubleNA60, xmin, xmax, GetNPar(kDoubleExp)+GetNPar(kNA60)+1);
        break;

        case 520:
            fDistrib =  new TF1("fitDistrib", Exp_DoubleCBext, xmin, xmax, GetNPar(kExp)+GetNPar(kCBExtended)+1);
        break;

        case 530:
            fDistrib =  new TF1("fitDistrib", Exp_DoubleNA60, xmin, xmax, GetNPar(kExp)+GetNPar(kNA60)+1);
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
    fDistrib->SetParLimits(nb_par_bg+nb_par_signal,0, 200000000);
    fDistrib->SetParName(nb_par_bg+nb_par_signal,name_par_signal2[0]);
	return fDistrib;
}