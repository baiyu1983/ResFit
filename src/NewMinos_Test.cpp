#include "TTree.h"
#include "TFile.h"
#include <string>
#include <iostream>
#include <time.h>
#include "TMatrixD.h"
#include "Header.h"
//#include "SetMatrix.h"
#include "CalDerivation.h"
#define EPSILON 1.0E-4
#define DELTA 1.0

using namespace std;
using namespace ResFit;

void NewMinos_Test(string inputfile, string treename,  int NRES, int NData, Long64_t jentry, int iipar){
     int NPar = NRES*4;
     int NSol = 1;
     for(int i=0;i<NRES;i++){
           NSol = NSol*2;
     }

     NSol /= 2;
     int ipar = (iipar>2)?iipar-1:iipar;

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

     Chi2_Fitter test_fitter;

     for(int i=0; i<NRES ;i++){
          string resname("res_");
          resname.append(to_string(i));
          resonance tmp;
          tmp._mass  = value[0][i][0];
          tmp._width = value[0][i][1];
          tmp._phase = value[0][i][2];
          tmp._branchratio = value[0][i][3];
          test_fitter.AddRes(resname,tmp);
     }

//     vector<Bdata> v_data;
     double binwidth = (5.1-3.4)/double(NData);
     for(int i=0;i<NData;i++){
          Bdata tmp;
          tmp.pdata = Data[i];
          tmp.pdataerror = DataError[i];
          tmp.penergy = 3.4+double(i)*binwidth+0.5*binwidth;
          test_fitter.AddData(tmp);
     }
     struct timespec rup;
     long time_minos = 0;
     int gt = clock_gettime(CLOCK_MONOTONIC, &rup);
     if(gt == -1){
           printf("Cannot get system reource usage information! \n");
           pthread_exit(NULL);   //terminate calling thread!
     }
     long mark_start_s = rup.tv_sec;
     long mark_start_ns = rup.tv_nsec;

     test_fitter.FixPhaseAngle();
     test_fitter.ReSet(); 
     test_fitter.UpdateHess();
     double chi2_min = test_fitter.GetChi2Min();
     std::cout<<"chi2_min = "<<chi2_min<<std::endl;
     double Delta=3.0;
     test_fitter.SetLimit(std::pair<std::string,ResFit::respara>("res_0",Branchratio),false,Delta);
/*
     test_fitter.InitialMinos(string("res_0"),Width,false,Delta);
     bool finished=false;
     while(!finished){
        finished = test_fitter.UpdateLimitOnce(std::pair<std::string,ResFit::respara>("res_0_up",Width),chi2_min,false,Delta);
     }
*/
     gt = clock_gettime(CLOCK_MONOTONIC, &rup);
     assert(gt != -1);
     long mark_end_s = rup.tv_sec;
     long mark_end_ns = rup.tv_nsec;
     time_minos += ((mark_end_s-mark_start_s)*1000000000+(mark_end_ns-mark_start_ns));
     std::cout<<"Time : "<<time_minos*1.0E-9<<std::endl;
}

int main(int argc, char * argv[]){
     if(argc!=6){std::cout<<"Error.Usage: ./chi2_test number_of_resonances input_file input_tree_name entry_idx para_idx" <<std::endl; return -1;}
     const int _nresonance = atoi(argv[1]);
     int entry_idx= atoi(argv[4]);
     int para_idx = atoi(argv[5]);
     NewMinos_Test(string(argv[2]), string(argv[3]), _nresonance, 100,  entry_idx, para_idx);
     return 0;
}
