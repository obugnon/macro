#include <vector>
#include <iostream>
#include <algorithm>
#include <stdio.h>
#include "TVector.h"
#include "TF1.h"
#include "TString.h"
#include "TMath.h"

#include "FitFunctions.C"


//___________________________________________________________________________________________________________
// Tails parameter for LHC15o name+source_function
	//Double_t emb_NA60[8]={-0.45,	0.01,	 0.5, 	0.26, 		2.7,	 0.0,	 0.58, 0.3}; //initial
	Double_t emb_NA60[8]={-0.446, 0.005, 0.529, 0.256, 2.653, 0.005, 0.583, 0.254};
	Double_t emb_CB[4]={0.87, 4.5, 2.5, 2.7}; //alphaL,nL,alphaR,nR
	Double_t ppData_CB[4]={0.94, 36.8, 9.40, 5.3}; //initial one
	//Double_t ppData_CB_qVWG_fitrangelarge[4] = {0.95, 21.6, 12., 1.6};
	//Double_t ppData_CB_qVWG_fitrangesmall[4] = {0.94,36.8, 9.4, 5.3};
	//Double_t ppData_CB_Pol2Pol3_fitrangelarge[4] = {0.94, 27.3, 8.85, 5.1};
	//Double_t ppData_CB_Pol2Pol3_fitrangesmall[4] = {0.93, 50.2, 5.95, 17.7};
	Double_t coherent_CB[4] = {0.803, 120, 2.32, 7.2};
	Double_t incoherent_CB[4] = {0.86, 11.9, 2.36, 4.1};
	Double_t geant4_CB[4]={0.95, 3.8, 2.3, 5.};
	Double_t geant4_NA60[8]={-0.072, 0.008, 0.763,  0.347, 2.201, 0.034, 1.157, 0.453};
	Double_t coherent_NA60[8] = {0.018, 0.763, 0.333, 0.104, 1.686, 0.454, -0.382, 2.065};
	Double_t incoherent_NA60[8] = {0.009, 0.628, 0.285, 0.012, 0.813, 0.336, -0.553, 2.449};
	

//Initial parameters to fit background
	Double_t par_VWG[4]={1., 1.235, 1.600, 2.928};
	TString name_VWG[4]={"Norm_{VWG}", "#mu_{VWG}", "#alpha_{VWG}", "#beta_{VWG}"};

	Double_t par_VWGQuadratic[5]={1., 2.15, 0.5, 0.14, 0.075};
	//Double_t par_VWGQuadratic[5]={1., 1.235, 1.600, 2.928, -2.894}; //50-70// 30-50
	//Double_t par_VWGQuadratic[5]={1., 2.3, 0.1, 1.82, -0.79};//70-90// 50-70 70-90
	//Double_t par_VWGQuadratic[5]={1., 2.3, 0.1, 1.82, -0.79};//70-90
	//Double_t par_VWGQuadratic[5]={1., 2.3, 0.5, 0.3, -0.79}; //0-10, 10-30, 30-50// 0-10, 10-30
	//Double_t par_VWGQuadratic[5]={1., 2., 0.4, 3, -2.}; // //80-90 
	//Double_t par_VWGQuadratic[5]={1., 1., 0.3, 0.5, -0.1}; //80-90 //
	TString name_VWGQuadratic[5]={"Norm_{VWG}", "#mu_{VWG}", "#alpha_{VWG}", "#beta_{VWG}", "#gamma_{VWG}"};

	
	Double_t par_Pol2Pol3[6]={1., 1., 1., 1., 1., 1.};
	//Double_t par_Pol2Pol3[6]={1., -0.35, 0.035, 6.5, -5.1, 1.04};
	TString name_Pol2Pol3[6]={"Norm", "a1", "a2", "b1", "b2", "b3"};

//Initial parameters to fit signal
	Double_t par_CBjpsi[3]={1., 3.1, 0.065};
	TString name_CBjpsi[7]={"Norm_{J/#psi}", "#mu_{J/#psi}", "#sigma_{J/#psi}", "#alpha^{L}_{J/#psi}","n^{L}_{J/#psi}", "#alpha^{R}_{J/#psi}", "n^{R}_{J/#psi}"};

	Double_t par_CBpsi[3]={1.7, par_CBjpsi[1]+3.686109-3.096916, par_CBjpsi[2]*1.05};
	TString name_CBpsi[7]={"Norm_{#psi'}", "#mu_{#psi'}", "#sigma_{#psi'}", "#alpha^{L}_{#psi'}","n^{L}_{#psi'}", "#alpha^{R}_{#psi'}", "n^{R}_{#psi'}"};
	//TString name_CBpsi[3]={"Norm_{#psi'}", "#mu_{#psi'}", "#sigma_{#psi'}"};

	Double_t par_NA60jpsi[3]={1.,3.1, 0.065};
	TString name_NA60jpsi[11]={"Norm_{J/#psi}", "#mu_{J/#psi}", "#sigma_{J/#psi}", "#alpha^{L}_{J/#psi}", "p1^{L}_{J/#psi}", "p2^{L}_{J/#psi}", "p3^{L}_{J/#psi}", "#alpha^{R}_{J/#psi}", "p1^{R}_{J/#psi}", "p2^{R}_{J/#psi}", "p3^{R}_{J/#psi}"};
		
	Double_t par_NA60psi[3]={1.7, par_NA60jpsi[1]+3.686109-3.096916, par_NA60jpsi[2]*1.05};
	TString name_NA60psi[11]={"Norm_{#psi'}", "#mu_{#psi'}", "#sigma_{#psi'}", "#alpha^{L}_{#psi'}", "p1^{L}_{#psi'}", "p2^{L}_{#psi'}", "p3^{L}_{#psi'}", "#alpha^{R}_{#psi'}", "p1^{R}_{#psi'}", "p2^{R}_{#psi'}", "p3^{R}_{#psi'}"};
	//TString name_NA60psi[3]={"Norm_{#psi'}", "#mu_{#psi'}", "#sigma_{#psi'}"};

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
TF1* SignalFunction(Efunction fSig, Etails fTails, Epart fPart, Double_t xmin, Double_t xmax, Bool_t isLeftFixed = kFALSE, Bool_t isRightFixed = kFALSE)
{
    Double_t* par_signal;
    Double_t* par_tails;
    TString* name_par_signal;

    Int_t nb_par_signal = GetNPar(fSig);

	TF1* fSignal;
    
    switch (fSig)
    {
        case kCBExtended:
            if (fPart == kJPsi) 
            {
                par_signal = par_CBjpsi;
                name_par_signal = name_CBjpsi;
            }
            else if (fPart == kPsi2S)
            {
                par_signal = par_CBpsi;
                name_par_signal = name_CBpsi;
            }
            if (fTails == kEMB) par_tails = emb_CB;
            else if (fTails == kPP) par_tails = ppData_CB;
            else if (fTails == kSTARLIGHTcoh) par_tails = coherent_CB;
            else if (fTails == kSTARLIGHTincoh) par_tails = incoherent_CB;

            fSignal = new TF1("fitSignal", CrystalBallExtended, xmin, xmax, GetNPar(kCBExtended));
        break;
    
        case kNA60:
            if (fPart == kJPsi) 
            {
                par_signal = par_NA60jpsi;
                name_par_signal = name_NA60jpsi;
            }
            else if (fPart == kPsi2S)
            {
                par_signal = par_NA60psi;
                name_par_signal = name_NA60psi;
            }
        
            if (fTails == kEMB) par_tails = emb_NA60;
            else if (fTails == kSTARLIGHTcoh) par_tails = coherent_NA60;
            else if (fTails == kSTARLIGHTincoh) par_tails = incoherent_NA60;

            fSignal = new TF1("fitSignal", NA60, xmin, xmax, GetNPar(kNA60));

        break;
    }

	for (int i=0; i<nb_par_signal; i++)
	{
        //Norm, mean and sigma JPsi
        if (i < 3) fSignal->SetParameter(i, par_signal[i]);
        //Tails Jpsi left
        else if ( i < 3+(nb_par_signal-3)/2)
        {
            if(isLeftFixed == kFALSE) fSignal->SetParameter(i, par_tails[i-3]);
            else fSignal->FixParameter(i, par_tails[i-3]);
        }
        //Tails Jpsi right
        else
        {
            if(isRightFixed == kFALSE) fSignal->SetParameter(i, par_tails[i-3]);
            else fSignal->FixParameter(i, par_tails[i-3]);
        }
        fSignal->SetParName(i, name_par_signal[i]);
    }
    
	return fSignal;
}

//___________________________________________________________________________________________________________
//___________________________________________________________________________________________________________
TF1* DistributionFunction(Efunction fBG, Efunction fSig, Etails fTails, Double_t xmin, Double_t xmax, Bool_t isLeftFixed = kFALSE, Bool_t isRightFixed = kFALSE)
{
	Double_t* par_bg;
    Double_t* par_signal1;
    Double_t* par_signal2;
    Double_t* par_tails;
	TString* name_par_bg;
    TString* name_par_signal1;
    TString* name_par_signal2;

    Int_t nb_par_bg = GetNPar(fBG);
    Int_t nb_par_signal = GetNPar(fSig);
    //Int_t nb_par = nb_par_bg + 2*nb_par_signal;
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
            par_signal1 = par_CBjpsi;
            par_signal2 = par_CBpsi;
            name_par_signal1 = name_CBjpsi;
            name_par_signal2 = name_CBpsi;
        
            if (fTails == kEMB) par_tails = emb_CB;
            else if (fTails == kPP) par_tails = ppData_CB;
            else if (fTails == kSTARLIGHTcoh) par_tails = coherent_CB;
            else if (fTails == kSTARLIGHTincoh) par_tails = incoherent_CB;
        break;
    
        case kNA60:
            par_signal1 = par_NA60jpsi;
            par_signal2 = par_NA60psi;
            name_par_signal1 = name_NA60jpsi;
            name_par_signal2 = name_NA60psi;
        
            if (fTails == kEMB) par_tails = emb_NA60;
            else if (fTails == kSTARLIGHTcoh) par_tails = coherent_NA60;
            else if (fTails == kSTARLIGHTincoh) par_tails = incoherent_NA60;
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
        //Tails JPsi left
        else if ( j < 3+(nb_par_signal-3)/2)
        {
            if(isLeftFixed == kFALSE) fDistrib->SetParameter(j+nb_par_bg, par_tails[j-3]);
            else fDistrib->FixParameter(j+nb_par_bg, par_tails[j-3]);
        }
        //Tails JPsi right
        else
        {
            if(isRightFixed == kFALSE) fDistrib->SetParameter(j+nb_par_bg, par_tails[j-3]);
            else fDistrib->FixParameter(j+nb_par_bg, par_tails[j-3]);
        }
        fDistrib->SetParName(j+nb_par_bg, name_par_signal1[j]);
    }
    //psi2s
    fDistrib->SetParameter(nb_par_bg+nb_par_signal,par_signal2[0]);
    fDistrib->SetParLimits(nb_par_bg+nb_par_signal,0., 10000);
    fDistrib->SetParName(nb_par_bg+nb_par_signal,name_par_signal2[0]);
    
	return fDistrib;
}