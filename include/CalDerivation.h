#pragma once
#ifndef HEADER_CALDERIVATION
#define HEADER_CALDERIVATION

#include "Header.h"
#include <complex>
#include <vector>

using namespace std;
using namespace ResFit;

inline void SetDataPhaseSpace(Bdata & data){
    data.pphasespace = AnalyticPhi23(data.penergy);
}

double FunctionPhaseSpace(double sqrts);
double FirstDerPhaseSpace(double sqrts);
double SecondDerPhaseSpace(double sqrts);
complex<double> FunctionBwigner(double sqrts, resonance para);
complex<double> FirstDerPhase(double sqrts, resonance para);
complex<double> FirstDerBrandRatio(double sqrts, resonance para);
complex<double> FirstDerWidth(double sqrts, resonance para);
complex<double> FirstDerMass(double sqrts, resonance para);
complex<double> FirstDerMass(double sqrts, resonance para, double phasespace, double phasespace_1stDerv);
complex<double> SecondDerWidth(double sqrts, resonance para);
complex<double> SecondDerMass(double sqrts, resonance para);
complex<double> SecondDerFr(double sqrts, resonance para);
complex<double> SecondDerMassWidth(double sqrts, resonance para);
complex<double> FunctionAmplitude(double sqrts, resonance para);
complex<double> FunctionFirstDerivative(double sqrts, resonance para, int idx);
complex<double> FunctionFirstDerivative(double sqrts, resonance para, int idx, double phasespace, double phasespace_1stDerv);
complex<double> FunctionSecondDerivative(double sqrts, resonance para, int idx, int jdx);
complex<double> FunctionSecondDerivative(double sqrts, resonance para, int idx, int jdx, double phasespace, double phasespace_1stDerv,double phasespace_2ndDerv);
#endif
