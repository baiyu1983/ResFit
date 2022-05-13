#ifndef SETMATRIX_H
#define SETMATRIX_H

#include "Header.h"
#include <vector>
using namespace std;
void SetFirstDMatrix(TMatrixD *output, std::vector<resonance> v_res, std::vector<Ddata> v_data, double M_MIN, double M_MAX);
void SetSecondDMatrix(TMatrixD *smatrix, vector<resonance> v_res, vector<Ddata> v_data, double M_MIN, double M_MAX);
double GetChi2(vector<resonance> v_res, vector<Ddata> v_data, double M_MIN, double M_MAX);
#endif
