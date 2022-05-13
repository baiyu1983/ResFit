#ifndef TRANSFORMER_H
#define TRANSFORMER_H
#include "TMatrixD.h"
#include "TComplex.h"
#include "Header.h"
#include <string>
#include <vector>

std::vector<TComplex> GetDVector(std::vector<TComplex> Coeffs, std::vector<TComplex> Poles);
std::vector<TComplex> GetZeros(std::vector<TComplex> Coeffs, std::vector<TComplex> Poles);
std::vector<TComplex> GetZerosFromD(std::vector<TComplex> DVector,std::vector<TComplex> Poles);//Get minus Zeros, which are alphas
std::vector<TComplex> GetSolution(std::vector<TComplex> Zeros, std::vector<TComplex> Coeffs, std::vector<TComplex> Poles,int idx);
TMatrixD GetErrorMatrix(TMatrixD inputMatrix, std::vector<TComplex> Zeros, std::vector<TComplex> Coeffs1, std::vector<TComplex> Coeffs2, std::vector<TComplex> Poles, int idx);
//TMatrixD GetA_Over_R_Matrix_AP(std::vector<TComplex> Coeffs, std::vector<TComplex> Poles,  std::vector<PhasespacePoint>); 
TMatrixD GetA_Over_R_Matrix_AP(std::vector<TComplex> Coeffs, std::vector<TComplex> Poles); //Partial(A,Phi)/Partial(R,Phi), Coeffs in A-Phi represention
TMatrixD GetA_Over_R_Matrix_RP(std::vector<TComplex> Coeffs, std::vector<TComplex> Poles); //Partial(A,Phi)/Partial(R,Phi), Coeffs in R-Phi represention
//TMatrixD GetRPhi_Over_XY_Matrix(std::vector<TComplex> Coeffs_r,std::vector<TComplex> Poles,double *);//
TMatrixD GetRPhi_Over_XY_Matrix(std::vector<TComplex> Coeffs_r,std::vector<TComplex> Poles);//Modified 6-1
TMatrixD GetKernelMatrix(std::vector<TComplex> Coeffs_r, std::vector<TComplex> Zeros, std::vector<TComplex> Poles);  //Transfer from Zx Zy to A0 and Alpha
TMatrixD GetHigherOrderCorrections(std::vector<TComplex> Coeffs, std::vector<TComplex> Poles, int target_solution_idx, double * Covariance_3Point, double * Covariance_4Point,std::vector<double> epsilons_2d, std::vector<double> epsilons_3d);
TMatrixD GetJacobianCombined(std::vector<TComplex> Zeros, std::vector<TComplex> Coeffs1, std::vector<TComplex> Coeffs2, std::vector<TComplex> Poles, int idx);
void FillSecondDerivation(double * result, std::vector<TComplex> Coeffs, std::vector<TComplex> Poles,  int target_solution_idx, std::vector<double> epsilons_2d);
void FillThirdDerivation(double * result , std::vector<TComplex> Coeffs, std::vector<TComplex> Poles,  int target_solution_idx, std::vector<double> epsilons_2d, std::vector<double> epsilons_3d);
#endif

