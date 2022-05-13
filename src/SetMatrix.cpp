#include <complex>
#include <iostream>
#include <vector>
#include <cmath>
#include "CalDerivation.h"
#include "TMatrixD.h"

#define UNIT 0.38936951

using namespace std;
double GetChi2(vector<resonance> v_res, vector<Ddata> v_data, double M_MIN, double M_MAX){
    const int NData = v_data.size(); const int NRes = v_res.size();
    double Chi2 = 0;
    complex<double> amplitude[NData*NRes];
    complex<double> sumAmplitude[NData];
    for (int s = 0; s < NData; s++) {
         double sqrts = M_MIN + (s + 0.5) * (M_MAX - M_MIN) / ((double)NData);
         sumAmplitude[s] = 0;
         for (int i = 0; i < NRes; i++) {
                 amplitude[s * NRes + i] = FunctionAmplitude(sqrts, v_res.at(i));
                 sumAmplitude[s] += amplitude[s * NRes + i];
/*
                 for(int j=0;j<4;j++){
                        firstDer[s * NRes * 4 + i * 4 + j] = FunctionFirstDerivative(sqrts, v_res.at(i), j);
                 }
*/
         }
//         firstDer[s * NRes * 4 + 2] = 0;
         double tmp = abs(sumAmplitude[s]);
         double data_exp = tmp*tmp*UNIT;
         double data_deviation = data_exp - v_data.at(s).pdata;
         Chi2 += data_deviation*data_deviation/(v_data.at(s).pdataerror*v_data.at(s).pdataerror);
    }   
    return Chi2;
}

void SetFirstDMatrix(TMatrixD *output, vector<resonance> v_res, vector<Ddata> v_data, double M_MIN, double M_MAX){
    const int NData = v_data.size(); const int NRes = v_res.size();
    for (int i = 0; i < 4 * NRes; i++) {
                (*output)[i][0] = 0.0;
    }
/*
    complex<double> * amplitude = new complex<double>[NData*NRes];
    complex<double> * sumAmplitude = new complex<double>[NData];
    complex<double>* firstDer = new complex<double>[NData * 4 * NRes];;
*/
    complex<double> amplitude[NData*NRes];
    complex<double> sumAmplitude[NData];
    complex<double> firstDer[NData * 4 * NRes];

    for (int s = 0; s < NData; s++) {
         double sqrts = M_MIN + (s + 0.5) * (M_MAX - M_MIN) / ((double)NData);
         sumAmplitude[s] = 0;
         for (int i = 0; i < NRes; i++) {
                 amplitude[s * NRes + i] = FunctionAmplitude(sqrts, v_res.at(i));
                 sumAmplitude[s] += amplitude[s * NRes + i];
                 for(int j=0;j<4;j++){
                        firstDer[s * NRes * 4 + i * 4 + j] = FunctionFirstDerivative(sqrts, v_res.at(i), j);
                 }
         }
         firstDer[s * NRes * 4 + 2] = 0;

    }
/*
    for (int s = 0; s < NData; s++) {
          sqrts = M_MIN + (s + 0.5) * (M_MAX - M_MIN) / numofData;
          for (int i = 0; i < numofPara; i++) {
                  for (int j = 0; j < 4; j++) {
                          firstDer[s * numofPara * 4 + i * 4 + j] = FunctionFirstDerivative(sqrts, mpara[i], j);
                  }
          }
          firstDer[s * numofPara * 4 + 2] = 0;
    }
*/
    for (int l = 0; l < NData; l++) {
          for (int i = 0; i < 4 * NRes; i++){
                  (*output)[i][0] += ((UNIT * 4) / pow(v_data.at(l).pdataerror, 2)\
                          * (UNIT * pow(abs(sumAmplitude[l]), 2) - v_data.at(l).pdata)\
                          * real(firstDer[l * 4 * NRes + i] * amplitude[l * NRes + i / 4] * conj(sumAmplitude[l])));
          }
    }
//    delete  [] amplitude; delete [] sumAmplitude; delete [] firstDer;
    return;
}


void SetSecondDMatrix(TMatrixD *smatrix, vector<resonance> v_res, vector<Ddata> v_data, double M_MIN, double M_MAX){
     const int NData = v_data.size(); const int NRes = v_res.size();
     for (int i = 0; i < 4 * NRes; i++) {
           for (int j = 0; j < 4 * NRes; j++) {
                   (*smatrix)[i][j] = 0.0;
           }
     }

   
     complex<double> amplitude[NData*NRes];
     complex<double> sumAmplitude[NData];
     complex<double> firstDer[NData * 4 * NRes];
     complex<double> secondDer[NData * (4 * NRes) * (4 * NRes)];

     for (int s = 0; s < NData; s++) {
         double sqrts = M_MIN + (s + 0.5) * (M_MAX - M_MIN) / ((double)NData);
//         std::cout<<"sqrt["<<s<<"] = "<<sqrts<<std::endl;
         sumAmplitude[s] = 0;
         for (int i = 0; i < NRes; i++) {
                 amplitude[s * NRes + i] = FunctionAmplitude(sqrts, v_res.at(i));
                 sumAmplitude[s] += amplitude[s * NRes + i];
                 for(int j=0;j<4;j++){
                        firstDer[s * NRes * 4 + i * 4 + j] = FunctionFirstDerivative(sqrts, v_res.at(i), j);
                 }
         }
         firstDer[s * NRes * 4 + 2] = 0;

     }

     for (int s = 0; s < NData; s++) {
                double sqrts = M_MIN + (s + 0.5) * (M_MAX - M_MIN) / ((double)NData);
                for (int i = 0; i < NRes; i++) {
                        for (int j = 0; j < NRes; j++) {
                                for (int k = 0; k < 4; k++) {
                                        for (int m = 0; m < 4; m++) {
                                                if (i != j)
                                                        secondDer[s * NRes * 4 * NRes * 4 + i * 4 * NRes * 4 + k * NRes * 4 + j * 4 + m] = 0;
                                                else {
                                                        secondDer[s * NRes * 4 * NRes * 4 + i * 4 * NRes * 4 + k * NRes * 4 + j * 4 + m] \
                                                                = firstDer[s * NRes * 4 + i * 4 + k] * firstDer[s * NRes * 4 + j * 4 + m] + FunctionSecondDerivative(sqrts, v_res.at(i), k, m);
                                                }
                                        }
                                }
                        }
                }
     }

     for (int i = 0; i < 4 * NRes; i++) {
              for (int j = 0; j < 4 * NRes; j++) {
                      for (int l = 0; l < NData; l++) {
                              (*smatrix)[i][j] += ((UNIT * 4) / pow(v_data.at(l).pdataerror, 2)) * \
                                      ((real(firstDer[l * NRes * 4 + i] * conj(firstDer[l * NRes * 4 + j]) * amplitude[l * NRes + i / 4] * conj(amplitude[l * NRes + j / 4])) \
                                              + real(secondDer[l * 4 * NRes * 4 * NRes + i * 4 * NRes + j] * amplitude[l * NRes + i / 4] * conj(sumAmplitude[l])))\
                                              * (UNIT * pow(abs(sumAmplitude[l]), 2) - v_data.at(l).pdata)\
                                              + real(firstDer[l * NRes * 4 + i] * amplitude[l * NRes + i / 4] * conj(sumAmplitude[l]))\
                                              * real(firstDer[l * NRes * 4 + j] * amplitude[l * NRes + j / 4] * conj(sumAmplitude[l])) * 2 * UNIT);
                      }
              }
      }

      return;
}

