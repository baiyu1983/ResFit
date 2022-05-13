#include "TTree.h"
#include "TFile.h"
#include "TComplex.h"
#include "TF1.h"
#include <string>
#include <iostream>
#include <algorithm>
#include <functional>
#include <time.h>
//#include "TMatrixD.h"
#include "Header.h"
//#include "SetMatrix.h"
#include "AnalyticPhi23.h"
#include "CalDerivation.h"
#define EPSILON 1.0E-4
#define DELTA 1.0
#define UNIT 0.38936951

using namespace std;
/*
void chi2::AddRes(resonance r){
    _states.push_back(r);_nres+=1;
    _poles.push_back(TComplex(r.mass*r.mass,-r.mass*r.width));
    TComplex tmp_coeffs(1,0);
    tmp_coeffs *= r.mass;
    tmp_coeffs *= sqrt(r.width*r.branchratio);
    tmp_coeffs *= TComplex::Exp(TComplex::I()*r.phase);
    tmp_coeffs *= sqrt(1.0/AnalyticPhi23(r.mass));
    _coeffs.push_back(tmp_coeffs);
}
*/
/*
template<typename T1,typename T2>
bool Amplitude<T1,T2>::operator() (T1 states, T2 data){
    TComplex this_amp(1,0);
    this_amp *= states.mass;
    this_amp *= sqrt(states.width*states.branchratio);
    this_amp *= TComplex::Exp(TComplex::I()*states.phase);
    this_amp *= sqrt(1.0/AnalyticPhi23(states.mass));
    TComplex BW = TComplex::One()*(data.penergy*data.penergy-states.mass*states.mass)+TComplex::I()*(states.mass*states.width);
    _amplitude += this_amp/BW;
    return true;
}
*/
/*
void chi2::operator() (Ddata data){
//    Amplitude<resonance,Ddata> ampsum;
//    ampsum = for_each((this->_states).begin(),(this->_states).end(),ampsum);
    //Amplitude<resonance,Ddata,void> ampsum = 
    //TComplex amplitude = 
//    for_each((this->_states).begin(),(this->_states).end(),bind2nd(ref(ampsum),data));
//    double xs = GetXS(ampsum._getAmplitude(),data.penergy);
    TComplex s=TComplex(data.penergy*data.penergy,0);
    TComplex ampsum(0,0);
    for(int i=0;i<this->_nres;i++){
        ampsum += (this->_coeffs).at(i)/(s-(this->_poles).at(i));
    }
    double xs = GetXS(ampsum,data.penergy);
    _chi2 += TMath::Power((xs-data.pdata)/data.pdataerror,2);
}
*/
/*
TComplex amp(vector<resonance> states, double sqrts){
       TComplex result;
       for(size_t i =0; i<states.size();i++ ){
            TComplex this_amp = TComplex::One();
            double M = states.at(i).mass;
            double WT = states.at(i).width;
            double Phase = states.at(i).phase;
            double BR = states.at(i).branchratio;
            this_amp *= M;
            this_amp *= sqrt(WT*BR);
            this_amp *= TComplex::Exp(TComplex::I()*Phase);
*/
/*
            if(M==4.0)
                this_amp *= sqrt(1.0/PS40);
            else if(M==4.2)
                this_amp *= sqrt(1.0/PS42);
*/
//            std::cout<<"M["<<i<<"] = "<<M<<std::endl;       
/*
            this_amp *= sqrt(1.0/AnalyticPhi23(M));
            TComplex BW = TComplex::One()*(sqrts*sqrts-M*M)+TComplex::I()*(M*WT);
            result += this_amp/BW;
        }
        return result;
}
*/
/*
double amp2TF(double* x, double * mypar){
    const int nres(*(mypar+0));
    if(nres ==0) return x[0]*mypar[4*nres+2]-mypar[4*nres+1];
    vector<resonance> v_state; v_state.clear();
    double PA[nres];
//    double Phase[nres];
    for(size_t i=0;i<nres;i++){
       PA[i] = mypar[4*i+3];
    }
*/
/*
    if(nres == 1){
        Phase[0] = PA[0];
    }else{
        Phase[0] = PA[0];
        for(size_t i=1;i<nres;i++){
             Phase[0] -= PA[i];
        }
        Phase[0] = Phase[0]/double(nres);
        for(size_t i=1;i<nres;i++){
             Phase[i] = PA[i]+Phase[0];
        }
    }
*/
/*
    for(size_t i=0; i<nres;i++){
        resonance tmp_state;
        tmp_state.mass = mypar[4*i+1];
        tmp_state.width = mypar[4*i+2];

//        tmp_state.phase = Phase[i]; // the input phase angle parameter needs to be altered

        tmp_state.phase = PA[i]; //Phase Anlge of resoances
        tmp_state.branchratio = mypar[4*i+4];
        v_state.push_back(tmp_state);
    }
    return UNIT*(amp(v_state,x[0])+x[0]*mypar[4*nres+2]+mypar[4*nres+1]).Rho2()*12.0*TMath::Pi()*AnalyticPhi23(x[0])/(x[0]*x[0]); //Interferance background
}
*/
/*
double GetChi2(chi2 &chi2_tmp, vector<Ddata> v_data){
    for_each(v_data.begin(),v_data.end(),chi2_tmp);
    return chi2_tmp._chi2;
}
*/
void Chi2Test(string inputfile, string treename, int NRES, int NData, Long64_t jentry){
     int NPar = NRES*4;
     int NSol = 1;
     for(int i=0;i<NRES;i++){
           NSol = NSol*2;
     }

     NSol /= 2;

     double Chi2[NSol];
     double value[NSol][NRES][4]; //[NSolution][NResonance][NVariable];
     double error[NSol][NRES][4];
     int ErrorStatus[NSol];
     double Data[NData];
     double DataError[NData];

     TFile finput(inputfile.c_str(),"read");
     TTree *tinput = (TTree*)finput.Get(treename.c_str());

     tinput->SetBranchAddress("Chi2", Chi2);
     tinput->SetBranchAddress("Value",value);
     tinput->SetBranchAddress("Error",error);
     tinput->SetBranchAddress("ErrorStatus",ErrorStatus);
     tinput->SetBranchAddress("Data",Data);
     tinput->SetBranchAddress("DataError",DataError);

     Int_t nentries = tinput->GetEntries();
     tinput->GetEntry(jentry);
//     vector<resonance> v_res_central; v_res_central.clear(); v_res_central.resize(0);
//     TF1 * fit_resonance = new TF1("fit_resonance",amp2TF,3.4,5.1,4*NRES+3);
//     fit_resonance->FixParameter(0,3);

//     ResFit::Chi2_Gradient chi2test;
     ResFit::Chi2_Hess chi2test;
/*   
     for(int i=0; i<NRES ;i++){

          resonance tmp;
          tmp.mass  = value[0][i][0];
          tmp.width = value[0][i][1];
          tmp.phase = value[0][i][2];
          tmp.branchratio = value[0][i][3];
          chi2test.AddRes(tmp);
     }
*/
/*
     fit_resonance->FixParameter(4*NRES+1,0);
     fit_resonance->FixParameter(4*NRES+2,0);
*/
     vector<ResFit::Bdata> v_data; 
     double binwidth = (5.1-3.4)/double(NData);
     for(int i=0;i<NData;i++){
          ResFit::Bdata tmp;
          tmp.pdata = Data[i];
          tmp.pdataerror = DataError[i];
          tmp.penergy = 3.4+double(i)*binwidth+0.5*binwidth;
          v_data.push_back(tmp);
     }
/*
     struct timespec rup;
     long time_minos = 0;
     int gt = clock_gettime(CLOCK_MONOTONIC, &rup);
     if(gt == -1){
           printf("Cannot get system reource usage information! \n");
           pthread_exit(NULL);   //terminate calling thread!
     }
     long mark_start_s = rup.tv_sec;
     long mark_start_ns = rup.tv_nsec;
*/

     for(int i=0; i<NRES ;i++){
          ResFit::resonance tmp;
          string resname("res_");
          resname.append(to_string(i));
          tmp._mass  = value[0][i][0];
          tmp._width = value[0][i][1];
          tmp._phase = value[0][i][2];
          tmp._branchratio = value[0][i][3];
          tmp._parameterized1 = true;
          tmp._parameterized2 = false;
          chi2test.AddRes(resname,tmp);
     }
//     chi2test.FixPhaseAngle("res_1");
     chi2test.FixPhaseAngle();
     std::map<std::string,ResFit::resonance> v_res = chi2test.GetResonance();
     int counter =0;

     for(std::map<std::string,ResFit::resonance>::const_iterator itr_res=v_res.begin(); itr_res!=v_res.end(); ++itr_res){
           std::cout<<"Resonance["<<counter<<"] , name = "<<itr_res->first<<", Mass = "<<(itr_res->second)._mass<<", Width = "<<(itr_res->second)._width<<", phaseangle = "<<(itr_res->second)._phase<<", fr = "<<(itr_res->second)._branchratio<<", pole = "<<(itr_res->second)._pole<<", coeff = "<<(itr_res->second)._coeff<<std::endl;
           counter+=1;
     }

     chi2test.ReSet();
//     for_each(v_data.begin(),v_data.end(),ResFit::SetDataPhaseSpace);
     struct timespec rup;
     long time_minos = 0;
     int gt = clock_gettime(CLOCK_MONOTONIC, &rup);
     if(gt == -1){
           printf("Cannot get system reource usage information! \n");
           pthread_exit(NULL);   //terminate calling thread!
     }
     long mark_start_s = rup.tv_sec;
     long mark_start_ns = rup.tv_nsec;
//     double myChi2 = GetChi2(chi2test, v_data);
     for_each(v_data.begin(),v_data.end(),SetDataPhaseSpace);
/*
     for(int i=0;i < v_data.size();i++){
           std::cout<<"data["<<i<<"] , energy = "<<v_data.at(i).penergy<<", phasespace = "<<v_data.at(i).pphasespace<<std::endl;
     }
*/
     for_each(v_data.begin(),v_data.end(),ref(chi2test));

     std::cout<<"Chi2 = "<<chi2test.GetChi2()<<std::endl;
     TMatrixD * chi2_grad = chi2test.GetChi2Gradient();
     for(int i=0;i<4*NRES;i++){
           std::cout<<"Gradient["<<i<<"] = "<<(*chi2_grad)[i][0]<<std::endl;
     }
     TMatrixD *hess_reduced = chi2test.GetHess_Reduced();
     hess_reduced->Invert();
     gt = clock_gettime(CLOCK_MONOTONIC, &rup);
     assert(gt != -1);
     long mark_end_s = rup.tv_sec;
     long mark_end_ns = rup.tv_nsec;
 /*  
     TMatrixD *hess = chi2test.GetHess();
     for(int i=0; i< 4*NRES; i++){
          for(int j=0;j<4*NRES;j++){
              std::cout<<(*(hess))[i][j]<<"\t";
          }std::cout<<endl;
     }

     TMatrixD *hess_reduced = chi2test.GetHess_Reduced();
     for(int i=0; i< 4*NRES-1; i++){
          for(int j=0;j<4*NRES-1;j++){
              std::cout<<(*(hess_reduced))[i][j]<<"\t";
          }std::cout<<endl;
     }
*/
     std::cout<<"Printing Covariance Matrix : "<<std::endl;
//     hess_reduced->Invert();
     for(int i=0; i< 4*NRES-1; i++){
          for(int j=0;j<4*NRES-1;j++){
              std::cout<<2.0E6*(*(hess_reduced))[i][j]<<"\t";
          }std::cout<<endl;
     }

     time_minos += ((mark_end_s-mark_start_s)*1000000000+(mark_end_ns-mark_start_ns));
     std::cout<<"Time : "<<time_minos*1.0E-9<<std::endl;

}

int main(int argc,char *argv[]){
     if(argc!=5){std::cout<<"Error.Usage: ./chi2_test number_of_resonances input_file input_tree_name entry_idx" <<std::endl; return -1;}
     const int _nresonance = atoi(argv[1]);
//     const int _solution_idx = atoi(argv[5]);
     int entry_idx= atoi(argv[4]);
/*
     TFile finput(argv[2],"read");
     TTree *tinput = (TTree*)finput.Get(argv[3]);
*/ 
     Chi2Test(string(argv[2]), string(argv[3]), _nresonance, 100, entry_idx);
     return 0;
}
