/*
 *  FitFunctions.C
 *
 *  Created by Ophelie Bugnon on 20/07/20.
 *
 */
#include <vector>
#include <stdio.h>
#include <iostream>
#include <algorithm>
#include "TF1.h"

Double_t PowerLawFunction(Double_t *x, Double_t *par)
{
    //par[0] = Normalization
	//par[1] = pT0
	//par[2] = n1
	//par[3] = n2
    return par[0]*x[0]/TMath::Power((1+TMath::Power(x[0]/par[1],par[3])),par[2]);
}

Double_t LevyFunction(Double_t *x, Double_t *par)
{
	//par[0] = Normalization
	//par[1] = n
	//par[2] = T
	//par[3] = mass
	return  x[0]*par[0]*(par[1]-1)*(par[1]-2)/( par[1]*par[2]*(par[1]*par[2]+par[3]*(par[1]-2) ) ) * TMath::Power(( 1 + (TMath::Sqrt(par[3]*par[3]+x[0]*x[0]) - par[3])/(par[1]*par[2])),-par[1]);
}

TF1* GetPowerLawFunction(Double_t initNorm, Double_t xmin, Double_t xmax, Bool_t isFixedPar3=kFALSE)
{
	TF1* funct = new TF1("PowerLawFunction", PowerLawFunction, xmin, xmax, 4);
	funct->SetParameters(initNorm, 2., 2., 2.);
	funct->SetParNames("Norm", "pt0", "n1", "n2");
	if(isFixedPar3) funct->FixParameter(3, 2.);
	return funct;
}

TF1* GetLevyFunction(Double_t initNorm, Double_t mass, Double_t xmin, Double_t xmax)
{
	TF1* funct = new TF1("LevyFunction", LevyFunction, xmin, xmax, 4);
	funct->SetParameters(initNorm,5., 0.6, mass);
	funct->SetParNames("Norm", "n", "T", "mass");
	funct->FixParameter(3, mass);
	return funct;
}

