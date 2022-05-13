#define _USE_MATH_DEFINES
#include "Header.h"
#include "AnalyticPhi23.h"
#include "CalDerivation.h"
#include <complex>
#include <vector>
#include <cmath>
using namespace std;
using namespace ResFit;

const complex<double> M_I = complex<double>{ 0,1 };

double FunctionPhaseSpace(double sqrts) {
	return AnalyticPhi23(sqrts);
}

double FirstDerPhaseSpace(double sqrts) {
	return AnalyticPhi23_derv1(sqrts) / AnalyticPhi23(sqrts);
}

double SecondDerPhaseSpace(double sqrts) {
	return AnalyticPhi23_derv2(sqrts) / AnalyticPhi23(sqrts);
}

complex<double> FunctionBwigner(double sqrts, resonance para) {
	return 1. / (sqrts * sqrts - para._pole);
}

complex<double> FirstDerPhase(double sqrts, resonance para) {
	return M_I;
}

complex<double> FirstDerBrandRatio(double sqrts, resonance para) {
	return 1. / (2. * para._branchratio);
}

complex<double> FirstDerWidth(double sqrts, resonance para) {
	return 1. / (2. * para._width) - M_I * para._mass * FunctionBwigner(sqrts, para);
}

complex<double> FirstDerMass(double sqrts, resonance para) {
	return (para._mass + sqrts * sqrts / para._mass) * FunctionBwigner(sqrts, para) - FirstDerPhaseSpace(para._mass) / 2.;
}

complex<double> FirstDerMass(double sqrts, resonance para, double PhaseSpace, double PhaseSpace_1st) {
        return (para._mass + sqrts * sqrts / para._mass) * FunctionBwigner(sqrts, para) - (PhaseSpace_1st/PhaseSpace) / 2.;
}

complex<double> SecondDerWidth(double sqrts, resonance para) {
	return -1. / (2. * para._width * para._width) - pow((para._mass * FunctionBwigner(sqrts, para)), 2);
}

complex<double> SecondDerMass(double sqrts, resonance para) {
	return (pow(para._mass, 4) + 4 * pow(para._mass * sqrts, 2) - pow(sqrts, 4) - 2. * M_I * pow(sqrts, 2) * para._mass * para._width) * pow(FunctionBwigner(sqrts, para) / para._mass, 2) \
		+ pow(FirstDerPhaseSpace(para._mass), 2) / 2. - SecondDerPhaseSpace(para._mass) / 2.;
}

complex<double> SecondDerMass(double sqrts, resonance para, double PhaseSpace, double PhaseSpace_1st, double PhaseSpace_2nd) {
        return (pow(para._mass, 4) + 4 * pow(para._mass * sqrts, 2) - pow(sqrts, 4) - 2. * M_I * pow(sqrts, 2) * para._mass * para._width) * pow(FunctionBwigner(sqrts, para) / para._mass, 2) \
                + pow((PhaseSpace_1st/PhaseSpace), 2) / 2. - (PhaseSpace_2nd/PhaseSpace) / 2.;
}

complex<double> SecondDerMassWidth(double sqrts, resonance para) {
	return -M_I * (para._mass * para._mass + sqrts * sqrts) * pow(FunctionBwigner(sqrts, para), 2);
}

complex<double> SecondDerFr(double sqrts, resonance para){
        return -1.0/(2*para._branchratio*para._branchratio);
}

complex<double> FunctionAmplitude(double sqrts, resonance para) {
	return para._mass / sqrts * sqrt(12. * M_PI * para._branchratio * para._width) * FunctionBwigner(sqrts, para) * \
		sqrt(FunctionPhaseSpace(sqrts) / FunctionPhaseSpace(para._mass)) * pow(M_E, (M_I * para._phase));
}

complex<double> FunctionFirstDerivative(double sqrts, resonance para, int idx) {
	switch(idx){
	case 0: return FirstDerMass(sqrts, para);
	case 1: return FirstDerWidth(sqrts, para);
	case 2: return FirstDerPhase(sqrts, para);
	case 3: return FirstDerBrandRatio(sqrts, para);
	}
 	return 0.;
}

complex<double> FunctionFirstDerivative(double sqrts, resonance para, int idx,double phasespace, double phasespace_1stDerv) {
        switch(idx){
        case 0: return FirstDerMass(sqrts, para, phasespace, phasespace_1stDerv);
        case 1: return FirstDerWidth(sqrts, para);
        case 2: return FirstDerPhase(sqrts, para);
        case 3: return FirstDerBrandRatio(sqrts, para);
        }
        return 0.;
}

complex<double> FunctionSecondDerivative(double sqrts, resonance para, int idx, int jdx) {
	if (idx == 0 && jdx == 0)
		return SecondDerMass(sqrts, para);
	else if ((idx == 0 && jdx == 1) || (idx == 1 && jdx == 0))
		return SecondDerMassWidth(sqrts, para);
	else if (idx == 1 && jdx == 1)
		return SecondDerWidth(sqrts, para);
        else if(idx==3 && jdx ==3)
                return SecondDerFr(sqrts,para);
	else
		return 0.;
}

complex<double> FunctionSecondDerivative(double sqrts, resonance para, int idx, int jdx, double phasespace, double phasespace_1stDerv, double phasespace_2ndDerv) {
        if (idx == 0 && jdx == 0)
                return SecondDerMass(sqrts, para, phasespace,phasespace_1stDerv,phasespace_2ndDerv);
        else if ((idx == 0 && jdx == 1) || (idx == 1 && jdx == 0))
                return SecondDerMassWidth(sqrts, para);
        else if (idx == 1 && jdx == 1)
                return SecondDerWidth(sqrts, para);
        else if(idx==3 && jdx ==3)
                return SecondDerFr(sqrts,para);
        else
                return 0.;
}
