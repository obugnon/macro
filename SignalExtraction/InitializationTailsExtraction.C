/*
 *  InitializationTailsExtraction.C
 *
 *  Created by Ophelie Bugnon on 28/01/20.
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
#include "TStyle.h"

#include "FitFunctions.C"


//___________________________________________________________________________________________________________
// Tails parameter for LHC15o name+source_function
	Double_t emb_NA60[8]={-1.12,	0.004,	 0.5, 	0.24, 1.93,	 0.004,	 0.62, 0.3}; //initial
	// Double_t emb_NA60[8]={-0.446, 0.005, 0.529, 0.256, 2.653, 0.005, 0.583, 0.254};
	Double_t emb_CB[4]={0.87, 4.5, 2.5, 2.7}; //low pt : alphaL,nL,alphaR,nR
    // Double_t emb_CB[4]={1.2, 3.27, 2., 2.9}; //high pt : alphaL,nL,alphaR,nR
    // Double_t emb_CB[4]={1.2, 3.8, 1.7, 3.42}; //differential rapidity : alphaL,nL,alphaR,nR

	Double_t ppData_CB[4]={0.94, 36.8, 9.40, 5.3}; //initial one
	//Double_t ppData_CB_qVWG_fitrangelarge[4] = {0.95, 21.6, 12., 1.6};
	//Double_t ppData_CB_qVWG_fitrangesmall[4] = {0.94,36.8, 9.4, 5.3};
	//Double_t ppData_CB_Pol2Pol3_fitrangelarge[4] = {0.94, 27.3, 8.85, 5.1};
	//Double_t ppData_CB_Pol2Pol3_fitrangesmall[4] = {0.93, 50.2, 5.95, 17.7};
	Double_t coherent_CB[4] = {0.803, 120, 2.32, 7.2};
	Double_t incoherent_CB[4] = {0.86, 11.9, 2.36, 4.1};
	Double_t geant4_CB[4]={0.95, 3.8, 2.3, 5.};
	Double_t geant4_NA60[8]={-0.072, 0.008, 0.763,  0.347, 2.201, 0.034, 1.157, 0.453};
	// Double_t coherent_NA60[8] = {0.018, 0.763, 0.333, 0.104, 1.686, 0.454, -0.382, 2.065};
    Double_t coherent_NA60[8] = {-0.342717, 0.186506, 1.25663, 0.15361, 2.1523, 0.0699892, 1.50153, 0.489294};
	// Double_t incoherent_NA60[8] = {0.009, 0.628, 0.285, 0.012, 0.813, 0.336, -0.553, 2.449};
	Double_t incoherent_NA60[8] = {-0.6001, 0.20903, 1.0602, 0.0311, 2.290, 0.1790 , 1.359, 0.106};


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
TF1* SignalFunction(Efunction fSig, Etails fTails, Epart fPart, Double_t xmin, Double_t xmax, Bool_t isLeftFixed = kFALSE, Bool_t isRightFixed = kFALSE)
{
    Double_t* par_signal=nullptr;
    Double_t* par_tails=nullptr;
    TString* name_par_signal=nullptr;

    Int_t nb_par_signal = GetNPar(fSig);

	TF1* fSignal=nullptr;
    
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
