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

Double_t RatioLevyFunction(Double_t *x, Double_t *par)
{
	
	//par[0] = dNdy
	//par[1] = n1
	//par[2] = T1
	//par[3] = mass1
    // Int_t param1 = par[1];
	Double_t bigCoef1=((par[1]-1)*(par[1]-2))/(par[1]*par[2]*(par[1]*par[2]+par[3]*(par[1]-2)));
	Double_t inPower1 = 1.0 + (TMath::Sqrt(x[0]*x[0]+par[3]*par[3])-par[3]) /(par[1]*par[2]);

	//par[4] = n2
	//par[5] = T2
	//par[6] = mass2

	Double_t bigCoef2 = ((par[4]-1)*(par[4]-2))/(par[4]*par[5]*(par[4]*par[5]+par[6]*(par[4]-2)));
	Double_t inPower2 = 1.0 + (TMath::Sqrt(x[0]*x[0]+par[6]*par[6])-par[6]) /(par[4]*par[5]);

	if((x[0] *bigCoef2 * TMath::Power(inPower2,(-1.0)*par[4]))==0) return 0;
	else return par[0] * (x[0] * bigCoef1 * TMath::Power(inPower1,(-1.0)*par[1]))/(x[0] *bigCoef2 * TMath::Power(inPower2,(-1.0)*par[4]));
}


Double_t Pol3Function(Double_t *x, Double_t *par)
{
	//par[0] = Normalization
	//par[1] = slope
	return  par[0] + par[1]*x[0] + par[2]*x[0]*x[0] + par[3]*x[0]*x[0]*x[0];
}

Double_t RatioPowerLaw(Double_t *x, Double_t *par)
{
	//par[0] = Normalization
	//par[1] = pT01
	//par[2] = n11
	//par[3] = n12

	//par[4] = pT01
	//par[5] = n21
	//par[6] = n22
	Double_t PowerLaw1 = x[0]/TMath::Power((1+TMath::Power(x[0]/par[1],par[3])),par[2]);
	Double_t PowerLaw2 = x[0]/TMath::Power((1+TMath::Power(x[0]/par[4],par[6])),par[5]);
    return par[0]*PowerLaw1/PowerLaw2;
}

TF1* GetRatioLevyFunction(Double_t initNorm, Double_t xmin, Double_t xmax)
{
	TF1* funct = new TF1("RatioLevyFunction", RatioLevyFunction, xmin, xmax, 7);
	// funct->SetParameters(initNorm, 9.35, 0.77, 1.03, 5.32, 0.29, 4.9); //work for 70-90  
	// funct->SetParLimits(0, 0., 1);
	// funct->SetParLimits(1, 5., 20);
	// funct->SetParLimits(2, 0.001, 1);
	// funct->SetParLimits(3, 0.5, 5);
	// funct->SetParLimits(4, 4., 10);
	// funct->SetParLimits(5, 0.001, 1);
	// funct->SetParLimits(6, 0.5, 20);
	// funct->SetParameters(initNorm, 9.54, 0.77, 1.03, 5.35, 0.29, 4.9);// work for 0-10
	// funct->SetParLimits(0, 0., 1);
	// funct->SetParLimits(1, 1., 200);
	// funct->SetParLimits(2, 0.001, 1);
	// funct->SetParLimits(3, 0.5, 5);
	// funct->SetParLimits(4, 1., 10);
	// funct->SetParLimits(5, 0.001, 1);
	// funct->SetParLimits(6, 0.5, 20);
	// funct->SetParameters(initNorm, 48.6, 0.97, 0.88, 8.63, 0.39, 4.6);// work for 50-70, 70-90, 10-30 and 0-10
	// funct->SetParLimits(0, 0., 1);
	// funct->SetParLimits(1, 1., 200);
	// funct->SetParLimits(2, 0.001, 1);
	// funct->SetParLimits(3, 0.5, 5);
	// funct->SetParLimits(4, 1., 10);
	// funct->SetParLimits(5, 0.001, 1);
	// funct->SetParLimits(6, 0.5, 20);
	funct->SetParameters(initNorm, 14, 0.92, 0.90, 5.66, 0.31, 5.27);// work for 50-70, 70-90, 30-50 and 0-10
	funct->SetParLimits(0, 0., 1);
	funct->SetParLimits(1, 1., 200);
	funct->SetParLimits(2, 0.001, 1);
	funct->SetParLimits(3, 0.5, 5);
	funct->SetParLimits(4, 1., 10);
	funct->SetParLimits(5, 0.001, 1);
	funct->SetParLimits(6, 0.5, 20);
	// funct->SetParameters(initNorm, 11.85, 0.82, 0.99, 5.54, 0.29, 5.29);// 70-90 from only stat
	// funct->SetParameters(initNorm, 9.63, 0.71, 1.09, 5.52, 0.27, 5.03);// 50-70 from only stat
	// funct->SetParameters(initNorm, 11.66, 0.84, 0.96, 4.56, 0.19, 7.95);// 30-50 from only stat
	// funct->SetParameters(initNorm, 11.98, 0.76, 1.08, 5.79, 0.27, 5.29);// 10-30 from only stat
	// funct->SetParameters(initNorm, 9.35, 0.77, 1.03, 5.32, 0.29, 4.92);// 0-10 from only stat
	// funct->SetParameters(initNorm, 18.03, 0.97, 0.92, 5.76, 0.29, 6.04);// 0-10 from only stat
	// funct->SetParLimits(0, 0., 1);
	// funct->SetParLimits(1, 5, 20);
	// funct->SetParLimits(2, 0.1, 1.5);
	// funct->SetParLimits(3, 0.1, 5);
	// funct->SetParLimits(4, 1, 10);
	// funct->SetParLimits(5, 0.1, 1);
	// funct->SetParLimits(6, 1., 20);
	funct->SetParNames("dN/dy", "n1", "T1", "m1", "n2", "T2", "m2");
	return funct;
}

TF1* GetRatioPowerLaw(Double_t initNorm, Double_t xmin, Double_t xmax)
{
	TF1* funct = new TF1("RatioPowerLaw", RatioPowerLaw, xmin, xmax, 7);
	funct->SetParameters(initNorm, 6, 5, 1.5, 3.4, 2, 2.7);//0-10, 10-30
	funct->SetParLimits(1, 1., 10);
	funct->SetParLimits(2, 1, 7.5);
	funct->SetParLimits(3, 0.5, 5);
	// funct->SetParameters(initNorm, 5, 7, 1.8, 3.4, 5, 2.7);// 
	funct->SetParNames("dN/dy", "pT01", "n11", "n12", "pT02", "n21", "n22");
	return funct;
}

TF1* GetPol3Function(Double_t initNorm, Double_t xmin, Double_t xmax)
{
	TF1* funct = new TF1("Pol3Function", Pol3Function, xmin, xmax, 4);
	funct->SetParameters(initNorm,-0.29,0.09,-0.004);
	// funct->SetParLimits(0, 0., 1000000);
	funct->SetParLimits(1, -100, -0.000000001);
	funct->SetParLimits(2, 0., 100);
	funct->SetParLimits(3, -100, -0.000000001);
	funct->SetParNames("Norm", "a1", "a2", "a3");
	return funct;
}