#include "Header.h"
#include "TTree.h"
#include "TFile.h"
#include <string>
#include <iostream>
#include <time.h>
#include "TMatrixD.h"
#include "SetMatrix.h"
#define EPSILON 1.0E-4
#define DELTA 1.0

vector<resonance> Initialization(double Chi2_Min, double Delta, TMatrixD * MDelta, std::vector<resonance> v_res, std::vector<Ddata> v_data, int iipar);
vector<resonance> UpdateOnce(double Chi2_Min, double Delta, std::vector<resonance> v_res, std::vector<Ddata> v_data, int iipar);
