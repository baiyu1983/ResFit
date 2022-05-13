#include "transformer.h"
#include "Header.h"
#include "TComplex.h"
#include <vector>
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <stdlib.h>
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TMatrixD.h"
#include "TH2D.h"
#include "AnalyticPhi23.h"
#include <Python.h>
#include "SetMatrix.h"
using namespace std;

vector<TComplex> VTransAmpToR(vector<TComplex> vCoeff, vector<TComplex> vPole){ //Transform the paramter in amplitude/phi to coeff phi
      vector<TComplex> v_result; v_result.clear();
      if(vPole.size()==vCoeff.size() && vPole.size()>1){
            for(size_t i=0; i< vCoeff.size();i++){
                  double tmp_K = sqrt((-vPole.at(i).Im())*sqrt(vPole.at(i).Re())/AnalyticPhi23(sqrt(vPole.at(i).Re())));
                  TComplex r_tmp = TComplex(tmp_K*sqrt(vCoeff.at(i).Rho()),vCoeff.at(i).Theta(),kTRUE);
                  v_result.push_back(r_tmp);
            }
      }else{
         cout<<"Wring input vectors"<<endl;
      }
      return v_result;
}

vector<TComplex> GetDVector(vector<TComplex> Coeffs, vector<TComplex> Poles){
      vector<TComplex> v_result; v_result.clear();
      if(Coeffs.size()==Poles.size()&&Coeffs.size()>1){
             TComplex A[Coeffs.size()]; A[0] = TComplex(0,0);
             for(int i=0;i<Coeffs.size();i++){
                  A[0]+=Coeffs.at(i);
                  if(i>=1){A[i]=Coeffs.at(i)*(Poles.at(i)-Poles.at(0));}
             }
             for(int k=0;k<Coeffs.size()-1;k++){
                  TComplex R=A[k+1]/A[0];
                  for(int i=0;i<Coeffs.size()-1;i++){
                        if(i!=k) R=R*(Poles.at(k+1)-Poles.at(i+1));
                  }
                  v_result.push_back(R);
             }
      }
      return v_result;
}

vector<TComplex> GetZerosFromD(vector<TComplex> DVector,vector<TComplex> Poles){
      vector<TComplex> v_result; v_result.clear();
      if(DVector.size()-Poles.size()==-1&&Poles.size()>=2){
             const int n = Poles.size()-1;
             TComplex PMatrix[n][n];
             TComplex PVector[n];
             for(int i=0;i<n;i++){
                for(int j=0;j<n;j++){
                   if(j==0)PMatrix[i][j] = TComplex(1,0);
                   else PMatrix[i][j] = PMatrix[i][j-1]*Poles.at(i+1);
                }
                PVector[i] = PMatrix[i][Poles.size()-2]*Poles.at(i+1);
             }
             TComplex RightSide[n];//Rightside of equation A3-p^n in the left
             for(int i=0;i<n;i++){
                   RightSide[i]= DVector.at(i)-PVector[i];
             }

             double PMatrix_Re[n][n], PMatrix_Im[n][n];  double RightSide_Re[n], RightSide_Im[n];
             for(int i=0; i<n;i++){
                 for(int j=0;j<n;j++){
                      PMatrix_Re[i][j] = PMatrix[i][j].Re();
                      PMatrix_Im[i][j] = PMatrix[i][j].Im();
                 }
                 RightSide_Re[i] = RightSide[i].Re();
                 RightSide_Im[i] = RightSide[i].Im();
             }
             PyObject *pModule,*pFunc;
             PyObject *pArgs, *pValue;
             PyObject *pMatrix_Re = PyTuple_New(n*n);
             PyObject *pMatrix_Im = PyTuple_New(n*n);
             PyObject *pVector_Re = PyTuple_New(n);
             PyObject *pVector_Im = PyTuple_New(n);
             PyRun_SimpleString("print(os.listdir())");
             pModule = PyImport_Import(PyUnicode_FromString("solution"));
//             pModule = PyImport_Import(PyString_FromString("solution"));
             pFunc = PyObject_GetAttrString(pModule,"getzeros");
             for(int i=0;i<n;i++){
                  PyTuple_SetItem(pVector_Re,i,PyFloat_FromDouble(RightSide_Re[i]));
                  PyTuple_SetItem(pVector_Im,i,PyFloat_FromDouble(RightSide_Im[i]));
                  for(int j=0;j<n;j++){
                       PyTuple_SetItem(pMatrix_Re,i*n+j,PyFloat_FromDouble(PMatrix_Re[i][j]));
                       PyTuple_SetItem(pMatrix_Im,i*n+j,PyFloat_FromDouble(PMatrix_Im[i][j]));
                  }
             }
             pArgs = PyTuple_New(4);
             PyTuple_SetItem(pArgs,0, pMatrix_Re);
             PyTuple_SetItem(pArgs,1, pMatrix_Im);
             PyTuple_SetItem(pArgs,2, pVector_Re);
             PyTuple_SetItem(pArgs,3, pVector_Im);

             pValue = PyObject_CallObject(pFunc, pArgs);
             for(size_t i=0;i<n;i++){
                  v_result.push_back(-TComplex(PyFloat_AsDouble(PyTuple_GetItem(pValue,2*i)),PyFloat_AsDouble(PyTuple_GetItem(pValue,2*i+1))));//return minus zeros, which are alphas
             }
      }
      return v_result;
}

vector<TComplex> GetZeros(vector<TComplex> input_parameter, vector<TComplex> input_poles){  // in r phi representation
     vector<TComplex> v_result; v_result.clear();
     int NRes = input_poles.size();
     vector<TComplex> DVector = GetDVector(input_parameter,input_poles);
     if(NRes>1){
         v_result = GetZerosFromD(DVector,input_poles);
     }
     return v_result;
}

std::vector<TComplex> GetSolution(std::vector<TComplex> Zeros, std::vector<TComplex> Coeffs, std::vector<TComplex> Poles,int idx){
      int NSol = 1;//Number of Solutions
      vector<int> v_mask;
      for(size_t i =0;i<Zeros.size();i++){
          v_mask.push_back((idx&NSol)>>i);
          NSol *=2;
      }



      TComplex A0(0,0);
      for(size_t i=0; i<Coeffs.size();i++){
          A0 += Coeffs.at(i);
      }

      vector<TComplex> output; output.clear();
      output.push_back(A0);
      for(size_t i=1; i<Coeffs.size();i++){
           TComplex z = Coeffs.at(i);
           for(size_t j=0;j<v_mask.size();j++){
                 if(v_mask[j] ==0)continue;
                 else{
                      z = z*(Poles.at(i)+TComplex::Conjugate(Zeros.at(j)))/(Zeros.at(j)+Poles.at(i));
                 }
           }
           output.push_back(z);
      }

      for(size_t i=1;i<output.size();i++){
           output.at(0) -= output.at(i);
      }
      return output;
}

TMatrixD GetA_Over_R_Matrix_RP(vector<TComplex> Coeffs, vector<TComplex> Poles){
      const int n = 4*Poles.size();
      TMatrixD outputM(n,n);
      for(int i=0;i<n;i++){
           double M = sqrt(Poles.at(int(i/4)).Re());
           double W = -Poles.at(int(i/4)).Im()/M;
/*
           double PhaseSpace = GetPhaseSpace(M,PS);
           double PhaseSpaceDerivation = GetPhaseSpaceDerivation(M,PS);
*/
           double PhaseSpace = AnalyticPhi23(M);
           double PhaseSpaceDerivation = AnalyticPhi23_derv1(M);
           double A = Coeffs.at(int(i/4)).Rho2()*PhaseSpace/(M*M*W);
           for(int j=0;j<n;j++){
                int iv = i%4; int jv = j%4;
                if(i/4!=j/4)outputM[i][j] = 0;
                else if(iv!=3){
                      if(i!=j)outputM[i][j] =0;
                      else outputM[i][j] = 1.0;
                }else{
                      if(jv == 0){
                           outputM[i][j] = (-2.0/M+PhaseSpaceDerivation/PhaseSpace)*A;
                      }else if(jv == 1){
                           outputM[i][j] = -A/W;
                      }else if(jv == 2){
                           outputM[i][j] = 0;
                      }else{
                           outputM[i][j] = 2*A / Coeffs.at(int(i/4)).Rho();
                      }
                }
           }
      }
      outputM.Invert();
      return outputM;
}

TMatrixD GetRPhi_Over_XY_Matrix(vector<TComplex> Coeffs, vector<TComplex> Poles){//Modified 6-1
      const int n =4* Coeffs.size();
      TMatrixD outputM(n,n);
      
      for(int i=0;i<n;i++){
         int ir = i/4; int iv = i%4;//iv= {Mass*Mass, -Mass*Width, X,Y}
         double M = sqrt(Poles.at(ir).Re());
         double W = -Poles.at(ir).Im()/M;
/*
         cout<<"In GetRPhi_OVer_XY_Matrix , Resonance["<<ir<<"] , Re(Coeff) = " <<Coeffs.at(ir).Re()<<", Im(Coeff) = "<<Coeffs.at(ir).Im()<<endl;
*/    
         for(int j=0;j<n;j++){
             int jr = j/4; int jv = j%4; //jv = {Mass, Width, R,Phi}
           
             outputM[i][j] = 0; double norm(n/4);
             if(jr>0 && (jv==2)){
                  double norm1(-(n/4));double norm2 = 1+1.0/norm1; norm = (ir == jr)?1.0/norm2:norm1;
             }
             if(iv==0){ //----M^2 
                  if(jv ==0) outputM[i][j] = (ir==jr)?2*M:0;
             }else if(iv == 1){ //----- -MW
                  if(jv == 0) outputM[i][j] = (ir==jr)?-W:0;
                  if(jv == 1)outputM[i][j] = (ir==jr)?-M:0;
             }else if(iv ==2){//-------- X
                  if(jv ==2 ) outputM[i][j] = (ir==jr)?-Coeffs.at(ir).Im():0;//Modified 6-1
                  if(jv == 3) outputM[i][j] = (ir==jr)?Coeffs.at(ir).Re()/Coeffs.at(ir).Rho():0;
             }else if(iv==3){//--------Y
                  if(jv == 2) {
                         outputM[i][j] = (ir==jr)?Coeffs.at(ir).Re():0;//Modified 6-1
                  }
                  if(jv == 3) outputM[i][j] = (ir==jr)?Coeffs.at(ir).Im()/Coeffs.at(ir).Rho():0;
             }
         }
      }

      return outputM;
}

TMatrixD GetKernelMatrix(vector<TComplex> Coeffs_r, vector<TComplex> Zeros, vector<TComplex> Poles){
      int n = Zeros.size()+1;
      double T_alpha_Re[4*n*n], T_alpha_Im[4*n*n] ,T_z_Re[4*n*n], T_z_Im[4*n*n];
      TComplex A0(0,0);
      for(int i=0;i<n;i++){
          A0=A0+Coeffs_r.at(i);
      }

      TComplex InverseA0 = TComplex(1.0,0.0)/A0;

      for(int i =0; i<2*n; i++){
            int ir = i/2; int iv = i%2; //ir is the index of resonances, iv is [pole, coeff]
            TComplex pksum(0,0);//Sum of 1/(pi-pj), j!=i
            for(int m =0; m<n;m++){
                  if(ir == m)continue;
                  pksum += TComplex(1,0)/(Poles.at(ir)-Poles.at(m));
            }
            TComplex alphapksum(0,0);// Sume of 1/(alpha_i+p_k)
            for(int m =0;m<Zeros.size();m++){
                 alphapksum += TComplex(1,0)/(Zeros.at(m)+Poles.at(ir));
            }
            for(int j =0;j<2*n;j++){
                 int jr = j/2;int jv = j%2;
                 if(iv == 0){
                     if(i!=j) {T_z_Re[i*2*n+j] =0; T_z_Im[i*2*n+j] =0; T_alpha_Re[i*2*n+j] =0 ; T_alpha_Im[i*2*n+j] =0;}
                     else {T_z_Re[i*2*n+j] =1; T_z_Im[i*2*n+j] =0; T_alpha_Re[i*2*n+j] =1 ; T_alpha_Im[i*2*n+j] =0;}
                 }else{
                     if(ir == 0){
                          T_z_Re[i*2*n+j] = (jv==1)?1.0:0.0; T_z_Im[i*2*n+j] = 0.0;
                          T_alpha_Re[i*2*n+j] = (j==1)?1.0:0.0;  T_alpha_Im[i*2*n+j] = 0;
                     }
                     else{
                          TComplex tmp_talpha(0,0);
                          if(jr>0){
                              if(jv == 1){
                                 tmp_talpha = TComplex(TComplex(1.0,0)/(Zeros.at(jr-1)+Poles.at(ir)));
                              }
                              else if(jr == ir){
                                 tmp_talpha = alphapksum;
                              }
                          }else{
                              tmp_talpha = (jv==1)?InverseA0:TComplex(0,0);
                          }
                          T_alpha_Re[2*n*i+j] = tmp_talpha.Re();T_alpha_Im[2*n*i+j] = tmp_talpha.Im();
                          TComplex tmp_tz(0,0);
                          if(jv == 1){
                              tmp_tz = (i==j)?TComplex(1.0)/Coeffs_r.at(ir):TComplex(0,0);
                          }else{
                              tmp_tz = (jr == ir)?pksum:TComplex(-1.0,0)/(Poles.at(ir)-Poles.at(jr));
                          }
                          T_z_Re[2*n*i+j] = tmp_tz.Re();  T_z_Im[2*n*i+j] = tmp_tz.Im();
                     }
                 }
            }
      }

      PyObject *pModule,*pFunc;
      PyObject *pArgs, *pValue;
      PyObject *pTAlpha_Re = PyTuple_New(4*n*n);
      PyObject *pTAlpha_Im = PyTuple_New(4*n*n);
      PyObject *pTz_Re = PyTuple_New(4*n*n);
      PyObject *pTz_Im = PyTuple_New(4*n*n);
      PyRun_SimpleString("print(os.listdir())");
      pModule = PyImport_Import(PyUnicode_FromString("solution"));
//      pModule = PyImport_Import(PyString_FromString("solution"));
      pFunc = PyObject_GetAttrString(pModule,"kernelmatrix");

      for(int i=0;i<2*n;i++){
           for(int j=0;j<2*n;j++){
                PyTuple_SetItem(pTAlpha_Re,i*2*n+j,PyFloat_FromDouble(T_alpha_Re[i*2*n+j]));
                PyTuple_SetItem(pTAlpha_Im,i*2*n+j,PyFloat_FromDouble(T_alpha_Im[i*2*n+j]));
                PyTuple_SetItem(pTz_Re,i*2*n+j,PyFloat_FromDouble(T_z_Re[i*2*n+j]));
                PyTuple_SetItem(pTz_Im,i*2*n+j,PyFloat_FromDouble(T_z_Im[i*2*n+j]));
           }
      }

      pArgs = PyTuple_New(4);
      PyTuple_SetItem(pArgs,0, pTAlpha_Re);
      PyTuple_SetItem(pArgs,1, pTAlpha_Im);
      PyTuple_SetItem(pArgs,2, pTz_Re);
      PyTuple_SetItem(pArgs,3, pTz_Im);

      pValue = PyObject_CallObject(pFunc, pArgs);
      TMatrixD outputM(4*n,4*n);

      int idx = 0;
      for(size_t i=0;i<2*n;i++){
             for(size_t j=0;j<2*n;j++){
                  outputM[2*i][2*j] = PyFloat_AsDouble(PyTuple_GetItem(pValue,2*idx));
                  outputM[2*i][2*j+1] = -PyFloat_AsDouble(PyTuple_GetItem(pValue,2*idx+1));
                  outputM[2*i+1][2*j] = PyFloat_AsDouble(PyTuple_GetItem(pValue,2*idx+1));
                  outputM[2*i+1][2*j+1] = PyFloat_AsDouble(PyTuple_GetItem(pValue,2*idx));
                  idx+=1;
             }
      }

      outputM.Invert(); //Modified 6-1
      return outputM;
}


TMatrixD GetJacobianCombined(vector<TComplex> Zeros, vector<TComplex> Coeffs1, vector<TComplex> Coeffs2, vector<TComplex> Poles, int idx){
     int NSol = 1;//Number of Solutions
     vector<int> v_mask;
     for(size_t i =0;i<Zeros.size();i++){
          v_mask.push_back((idx&NSol)>>i);
          NSol *=2;
     }

     idx = idx%NSol;
     int NRes = Poles.size();

     TMatrixD ARTransInput = GetA_Over_R_Matrix_RP(Coeffs1, Poles);

     TMatrixD RPhiToXYInput = GetRPhi_Over_XY_Matrix(Coeffs1,Poles); //Modified 6-1
     TMatrixD KernelInput = GetKernelMatrix(Coeffs1,Zeros,Poles);
     TMatrixD ARTransOutput = GetA_Over_R_Matrix_RP(Coeffs2, Poles);

     TMatrixD RPhiToXYOutput = GetRPhi_Over_XY_Matrix(Coeffs2,Poles); //Modified 6-1
     vector<TComplex> Zeros2; Zeros2.clear();
     for(int i=0;i<v_mask.size();i++){
           if(v_mask.at(i)!=0){Zeros2.push_back(TComplex::Conjugate(Zeros.at(i)));}
           else{Zeros2.push_back(Zeros.at(i));};
     }

     TMatrixD KernelOutput = GetKernelMatrix(Coeffs2,Zeros2,Poles);

     TMatrixD MMInput = KernelInput*RPhiToXYInput*ARTransInput;
     TMatrixD MMOutput = KernelOutput*RPhiToXYOutput*ARTransOutput;
     TMatrixD MRotation(4*Zeros.size()+4,4*Zeros.size()+4);

     MMOutput.Invert();
     TMatrixD TAlpha(4*Zeros.size()+4,4*Zeros.size()+4);

     for(int i =0;i<4*Zeros.size()+4;i++){
          for(int j =0;j<4*Zeros.size()+4;j++){
                if(i!=j)TAlpha[i][j] = 0;
                else if(i%4==3 && i>=4){ //Not sure if this always works properly
                     if(v_mask.at(i/4-1)!=0)
                         TAlpha[i][j] = -1;
                     else
                         TAlpha[i][j] = 1;
                }else{
                     TAlpha[i][j] = 1;
                }
          }
     }

     for(int i =0;i<4*Zeros.size()+4;i++){
          for(int j =0;j<4*Zeros.size()+4;j++){
              if(i==2)MRotation[i][j]=0;
              else if(i%4==2){
                      if(j==2)MRotation[i][j]=-1;
                      else if(i==j)MRotation[i][j]=1;
                      else MRotation[i][j]=0;
              }else{
                      MRotation[i][j] = (i==j)?1:0;
              }
          }
     }

     TMatrixD MTM = MRotation*MMOutput*TAlpha*MMInput;
     return MTM;
}

TMatrixD GetErrorMatrix(TMatrixD InputMatrix, vector<TComplex> Zeros, vector<TComplex> Coeffs1, vector<TComplex> Coeffs2, vector<TComplex> Poles, int idx){
     int NRes = Poles.size();
     TMatrixD inputMatrix(4*NRes,4*NRes);
     for(int i=0; i< 4*NRes;i++){
         for(int j =0; j< 4*NRes;j++){
              if(i==2||j==2)inputMatrix[i][j] = 0;
              else{
                  int idx = (i>1)?(i-1):i;
                  int jdx = (j>1)?(j-1):j;
                  inputMatrix[i][j] = InputMatrix[idx][jdx];
              }
         }
     }



     TMatrixD MTM = GetJacobianCombined(Zeros, Coeffs1, Coeffs2, Poles, idx);
     TMatrixD MTMT = MTM;  MTMT.Transpose(MTM);

     TMatrixD outputMatrix = MTM*inputMatrix*MTMT; //Modified 6-1

     TMatrixD OutputMatrix(4*NRes-1,4*NRes-1);

     for(int i=0;i < 4*NRes-1;i++){
         for(int j=0; j<4*NRes-1;j++){
                     int idx_j = (j>1)?j+1:j;
                     int idx_i = (i>1)?i+1:i;
                     OutputMatrix[i][j] = outputMatrix[idx_i][idx_j];
         }
     }
     return OutputMatrix;
}

int main(int argc,char *argv[]){
     if(argc!=6){cout<<"Error.Usage: ./transformer number_of_resonances input_file input_tree_name entry_idx solution_idx"<<endl;}
     const int _nresonance = atoi(argv[1]);
     const int _solution_idx = atoi(argv[5]);
     int entry_idx= atoi(argv[4]);
     TFile finput(argv[2],"read");
     TTree *tinput = (TTree*)finput.Get(argv[3]);
     
     const int _npar= 4*_nresonance;
     int nsolution = 1;
     for(int i=1;i<_nresonance;i++){
          nsolution = 2*nsolution;
     }
     const int _nsolution = nsolution;
 
     double Chi2[_nsolution];

     double Value[_nsolution][_nresonance][4];
     double Error[_nsolution][_nresonance][4];
     double Covariance[_nsolution][_nresonance][_nresonance][4][4];
     double Data[100];
     double DataError[100];
     
/*
     double Value[_nsolution*_nresonance*_npar];   
     double Covariance[_nsolution][_nresonance][_nresonance][_npar][_npar];
*/
     tinput->SetBranchAddress("Chi2",Chi2);
     tinput->SetBranchAddress("Value",Value);
     tinput->SetBranchAddress("Error",Error);
     tinput->SetBranchAddress("Covariance",Covariance);
     tinput->SetBranchAddress("Data",Data);
     tinput->SetBranchAddress("DataError",DataError);

     Int_t nentries = tinput->GetEntries();

     Py_Initialize();
     PyRun_SimpleString("import sys; import os; sys.path.append(os.getcwd())");

     vector<TComplex> poles; poles.clear();
     vector<TComplex> coeffs; coeffs.clear();

     tinput->GetEntry(entry_idx);
     vector<resonance> v_res_0;
     for(int i=0; i< _nresonance; i++){

             double tmp_M = Value[0][i][0]; double tmp_W = Value[0][i][1]; 
             double tmp_Phi = Value[0][i][2]; double tmp_Fr = Value[0][i][3];

/*
             double tmp_M = *(Value+i*4); double tmp_W = *(Value+i*4+1); 
             double tmp_Phi = *(Value+i*4+2); double tmp_Fr = *(Value+i*4+3);
*/
             TComplex tmp_poles = TComplex( tmp_M*tmp_M,-tmp_M*tmp_W);
             TComplex tmp_coeffs(0,0);
             tmp_coeffs = TComplex(tmp_Fr,tmp_Phi,kTRUE);
             std::cout<<"M["<<i<<"] = "<<tmp_M<<"\t"<<"W["<<i<<"] = "<<tmp_W<<"\t"<<"Phi["<<i<<"] = "<<tmp_Phi<<"\t"<<"FR["<<i<<"] = "<<tmp_Fr<<std::endl;
             poles.push_back(tmp_poles);  coeffs.push_back(tmp_coeffs);
             resonance tmp_res; tmp_res.mass=tmp_M; tmp_res.width = tmp_W; tmp_res.branchratio = tmp_Fr; tmp_res.phase = tmp_Phi; v_res_0.push_back(tmp_res);
     }
     
     vector<Ddata> v_data;
     for(int i=0;i<100;i++){
          Ddata tmp;
          tmp.pdata = Data[i];
          tmp.pdataerror = DataError[i];
          v_data.push_back(tmp);
     }

     TMatrixD InitialCovarianceMatrix(_npar-1,_npar-1);
     int idx = 0; int jdx = 0;
     for(int i=0;i<_npar;i++){
             if(i == 2) continue; //Phi[0]-Row
             jdx = 0;
             for(int j=0;j<_npar;j++){
                 if(j==2) continue; //Phi[0]-Columne
                 int in1 = i/4; int in2 = j/4; int ii1 = i%4; int ii2=j%4;

                 InitialCovarianceMatrix[idx][jdx] = Covariance[0][in1][in2][ii1][ii2];
/*
               InitialCovarianceMatrix[idx][jdx] = *(Covariance+in1*_nresonance*4*4+in2*4*4+ii1*4+ii2);
*/
                 jdx += 1;
             }
             idx += 1;
        } //End of Reading ErroMatrix

     cout<<"Input Error Matrix MIGRAD, Det = "<<InitialCovarianceMatrix.Determinant()<<endl;
     for(int i =0 ; i<4*_nresonance-1;i++){
           for(int j =0;j <4*_nresonance-1;j++){
                 cout<<InitialCovarianceMatrix[i][j]<<"    ";
           }
           cout<<"//"<<endl;
     }
     cout<<"Input Error MIGRAD"<<endl;
     for(int i=0;i<_npar;i++){
         if(i==2)continue;
         cout<<"Error Migrad["<<i<<"] : "<<"\t"<<Error[0][i/4][i%4]<<std::endl;
     }

     cout<<"Output Error MIGRAD"<<endl;
     for(int i=0;i<_npar;i++){
         if(i==2)continue;
         cout<<"Error Migrad["<<i<<"] : "<<"\t"<<Error[1][i/4][i%4]<<std::endl;
     }


     vector<TComplex> coeffs_r = VTransAmpToR(coeffs,poles);

     vector<TComplex> DVector = GetDVector(coeffs_r,poles);
 
     vector<TComplex> myZeros = GetZerosFromD(DVector,poles);
 
     vector<TComplex> coeffs_r2 = GetSolution(myZeros,coeffs_r,poles,_solution_idx);
 
     TMatrixD error_matrix_2 = GetErrorMatrix(InitialCovarianceMatrix, myZeros, coeffs_r, coeffs_r2, poles,_solution_idx);
 
     cout<<"Output Matrix From Prapation MIGRAD: "<<endl;
     for(int i=0; i<_npar-1;i++){
          for(int j=0; j<_npar-1; j++){
               cout<<error_matrix_2[i][j]<<"\t";
          }
          cout<<"//"<<endl;
     }
 
     cout<<"Output Matrix From Fitting MIGRAD"<<endl;
     for(int i=0;i<_npar;i++){
          for(int j=0;j<_npar;j++){
               cout<<Covariance[_solution_idx][i/4][j/4][i%4][j%4]<<"\t";
          }
          cout<<"//"<<endl;
     }

     cout<<"Printing Errors Prapation MIGRAD: "<<endl;
    for(int i=0;i<_npar-1;i++){
              cout<<"Error["<<i<<"]: = "<<"\t"<<sqrt(error_matrix_2[i][i])<<endl;
    }
 
    TMatrixD * secondderv_0 = new TMatrixD(_npar, _npar);
    SetSecondDMatrix(secondderv_0, v_res_0, v_data, 3.4, 5.1);
    TMatrixD *secondderv_0_reduced = new TMatrixD(_npar-1,_npar-1);    
    for(int i=0;i<_npar;i++){
         if(i==2)continue;
         int ii = (i<2)?i:i-1;
         for(int j=0; j<_npar;j++){
             if(j==2)continue;
             int jj = (j<2)?j:j-1;
             (*secondderv_0_reduced)[ii][jj] = (*secondderv_0)[i][j]/2.0;
         }
    }

    secondderv_0_reduced->Invert();

    cout<<"Input Matrix HESS: "<<endl;
    for(int i=0;i<_npar-1;i++){
       for(int j=0; j<_npar-1;j++){
            cout<<1.0E6*(*secondderv_0_reduced)[i][j]<<"\t";
       }
       cout<<"//"<<endl;
    }
    cout<<"Input Error HESS: "<<endl;
    for(int i=0;i<_npar-1;i++){
       cout<<"HESS Error["<<i<<"] : "<<"\t"<<1.0E3*sqrt((*secondderv_0_reduced)[i][i])<<endl;
    }

    vector<resonance> v_res_1;
    for(int i=0; i< _nresonance; i++){
             double tmp_M = Value[_solution_idx][i][0]; double tmp_W = Value[_solution_idx][i][1];
             double tmp_Phi = Value[_solution_idx][i][2]; double tmp_Fr = Value[_solution_idx][i][3];
//             poles.push_back(tmp_poles);  coeffs.push_back(tmp_coeffs);
             resonance tmp_res; tmp_res.mass=tmp_M; tmp_res.width = tmp_W; tmp_res.branchratio = tmp_Fr; tmp_res.phase = tmp_Phi; v_res_1.push_back(tmp_res);
     }
     TMatrixD * secondderv_1 = new TMatrixD(_npar, _npar);
     SetSecondDMatrix(secondderv_1, v_res_1, v_data, 3.4, 5.1);
     TMatrixD *secondderv_1_reduced = new TMatrixD(_npar-1,_npar-1);
     TMatrixD InitialCovarianceMatrix_1(_npar-1,_npar-1);


     for(int i=0;i<_npar;i++){
         if(i==2)continue;
         int ii = (i<2)?i:i-1;
         for(int j=0; j<_npar;j++){
             if(j==2)continue;
             int jj = (j<2)?j:j-1;
             (*secondderv_1_reduced)[ii][jj] = (*secondderv_1)[i][j]/2.0;
             InitialCovarianceMatrix_1[ii][jj] = (*secondderv_0)[i][j]/2.0E6;
         }
      }
   
      secondderv_1_reduced->Invert();
      InitialCovarianceMatrix_1.Invert();
 
      cout<<"Output Matrix HESS: "<<endl;
      for(int i=0;i<_npar-1;i++){
         for(int j=0; j<_npar-1;j++){
              cout<<1.0E6*(*secondderv_1_reduced)[i][j]<<"\t";
         }
         cout<<"//"<<endl;
      }
      cout<<"Output Error HESS: "<<endl;
      for(int i=0;i<_npar-1;i++){
         cout<<"HESS Error["<<i<<"] : "<<"\t"<<1.0E3*sqrt((*secondderv_1_reduced)[i][i])<<endl;
      }

      TMatrixD error_matrix_2_hess = GetErrorMatrix(InitialCovarianceMatrix_1, myZeros, coeffs_r, coeffs_r2, poles,_solution_idx);
      cout<<"Output Matrix From Prapation HESS: "<<endl;
      for(int i=0; i<_npar-1;i++){
           for(int j=0; j<_npar-1; j++){
                cout<<error_matrix_2_hess[i][j]<<"\t";
           }
           cout<<"//"<<endl;
      }
      cout<<"Outpu Error Prapation HESS:"<<endl;
      for(int i=0; i<_npar-1;i++){
           cout<<"HESS Error Prapation["<<i<<"]"<<"\t"<<sqrt(error_matrix_2_hess[i][i])<<endl;
      }

      Py_Finalize();
      delete secondderv_0;
      delete secondderv_0_reduced;
      delete secondderv_1;
      delete secondderv_1_reduced;
}
