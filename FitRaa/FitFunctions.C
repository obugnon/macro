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


Double_t WoodSaxonLikeFunction(Double_t *x, Double_t *par)
{
    //par[0] = Normalization
	//par[1] = DeltaRaa
	//par[2] = pT0
	//par[3] = sigmapT
    return par[0]*(1+ par[1]/(1+TMath::Exp((x[0]-par[2])/par[3])));
}

Double_t Pol1Function(Double_t *x, Double_t *par)
{
	//par[0] = Normalization
	//par[1] = slope
	return  par[0] + par[1]*x[0];
}

Double_t Pol3Function(Double_t *x, Double_t *par)
{
	//par[0] = Normalization
	//par[1] = slope
	return  par[0] + par[1]*x[0] + par[2]*x[0]*x[0] + par[3]*x[0]*x[0]*x[0];
}

TF1* GetWoodSaxonLikeFunction(Double_t initNorm, Double_t initPT0, Double_t xmin, Double_t xmax, Bool_t isFixedPar2=kFALSE)
{
	TF1* funct = new TF1("WoodSaxonLikeFunction", WoodSaxonLikeFunction, xmin, xmax, 4);
	funct->SetParameters(initNorm, 2., initPT0, 2.);
	funct->SetParNames("Norm", "Delta", "pT0", "sigma");
	funct->SetParLimits(2, 0., 20.);
	if(isFixedPar2) funct->FixParameter(2, initPT0);
	return funct;
}

TF1* GetPol1Function(Double_t initNorm, Double_t initSlope, Double_t xmin, Double_t xmax, Bool_t isFixedPar1=kFALSE)
{
	TF1* funct = new TF1("Pol1Function", Pol1Function, xmin, xmax, 2);
	funct->SetParameters(initNorm,initSlope);
	funct->SetParNames("Norm", "slope");
	if(isFixedPar1) funct->FixParameter(1, 0);
	return funct;
}


TF1* GetPol3Function(Double_t initNorm, Double_t xmin, Double_t xmax)
{
	TF1* funct = new TF1("Pol3Function", Pol3Function, xmin, xmax, 4);
	funct->SetParameters(initNorm,-0.4,4,-4);
	funct->SetParNames("Norm", "a1", "a2", "a3");
	return funct;
}
