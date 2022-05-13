#pragma once

#ifndef HEADER_ANALYTICPHI23
#define HEADER_ANALYTICPHI23

double AnalyticPhi23(double x, double m1 = 0.1395706, double m2 = 0.1395706, double m3 = 3.0969);
double AnalyticPhi23_derv1(double x, double m1 = 0.1395706, double m2 = 0.1395706, double m3 = 3.0969);
double AnalyticPhi23_derv2(double x, double m1 = 0.1395706, double m2 = 0.1395706, double m3 = 3.0969);
double Coeff(double x, double m1, double m2, double m3);
double Coeff1(double x, double m1, double m2, double m3);
double Coeff2(double x, double m1, double m2, double m3);
double Coeff3(double x, double m1, double m2, double m3);
double Coeff4(double x, double m1, double m2, double m3);
double Coeff_derv1(double x, double m1, double m2, double m3);
double Coeff1_derv1(double x, double m1, double m2, double m3);
double Coeff2_derv1(double x, double m1, double m2, double m3);
double Coeff3_derv1(double x, double m1, double m2, double m3);
double Coeff4_derv1(double x, double m1, double m2, double m3);
double Coeff_derv2(double x, double m1, double m2, double m3);
double Coeff1_derv2(double x, double m1, double m2, double m3);
double Coeff2_derv2(double x, double m1, double m2, double m3);
double Coeff3_derv2(double x, double m1, double m2, double m3);
double Coeff4_derv2(double x, double m1, double m2, double m3);
double comp_ellint_1_derv1(double x);
double comp_ellint_1_derv2(double x);
double comp_ellint_2_derv1(double x);
double comp_ellint_2_derv2(double x);
double comp_ellint_3_derv1_n(double x);
double comp_ellint_3_derv1_k(double x);
double comp_ellint_3_derv2_kk(double y, double x);
double comp_ellint_3_derv2_kn(double y, double x);
double comp_ellint_3_derv2_nk(double y, double x);
double comp_ellint_3_derv2_nn(double y, double x);
double QPlus_derv1(double x, double m1, double m2, double m3);
double QMinus_derv1(double x, double m1, double m2, double m3);
double QPlus_derv2(double x, double m1, double m2, double m3);
double QMinus_derv2(double x, double m1, double m2, double m3);
double K_derv1(double x, double m1, double m2, double m3);
double K_derv2(double x, double m1, double m2, double m3);
double A1_derv1(double x, double m1, double m2, double m3);
double A1_derv2(double x, double m1, double m2, double m3);
double A2_derv1(double x, double m1, double m2, double m3);
double A2_derv2(double x, double m1, double m2, double m3);
double comp_ellint_1_derv1_over_X(double x, double m1, double m2, double m3);
double comp_ellint_1_derv2_over_X(double x, double m1, double m2, double m3);
double comp_ellint_2_derv1_over_X(double x, double m1, double m2, double m3);
double comp_ellint_2_derv2_over_X(double x, double m1, double m2, double m3);
double comp_ellint_3_derv1_A1_over_X(double x, double m1, double m2, double m3);
double comp_ellint_3_derv2_A1_over_X(double x, double m1, double m2, double m3);
double comp_ellint_3_derv1_A2_over_X(double x, double m1, double m2, double m3);
double comp_ellint_3_derv2_A2_over_X(double x, double m1, double m2, double m3);

#endif