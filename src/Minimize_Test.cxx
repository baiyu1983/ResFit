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

void NewMinos_Test(string inputfile, string treename,  int NRES, int NData, Long64_t jentry, Long64_t NEvents){
     int NPar = NRES*4;
     int NSol = 1;
     for(int i=0;i<NRES;i++){
           NSol = NSol*2;
     }

     NSol /= 2;
     int ipar = (iipar>2)?iipar-1:iipar;
     double Data[NData];
     double DataError[NData];

     TFile finput(inputfile.c_str(),"read");
     TTree *tinput = (TTree*)finput.Get(treename.c_str());
   
     tinput->SetBranchStatus("Chi2", 0);
     tinput->SetBranchStatus("Value",0);
     tinput->SetBranchStatus("Error",0);
     tinput->SetBranchStatus("ErrorStatus",0);
     
     Int_t nentries = tinput->GetEntries();
     
     struct timespec rup;
     long time_minos = 0;
     int gt = clock_gettime(CLOCK_MONOTONIC, &rup);
     if(gt == -1){
           printf("Cannot get system reource usage information! \n");
           pthread_exit(NULL);   //terminate calling thread!
     }
     long mark_start_s = rup.tv_sec;
     long mark_start_ns = rup.tv_nsec;

     for(Long64_t i=jentry; i<jentry+NEvents;i++){
         tinput->GetEntry(i);
         Chi2_Fitter test_fitter;
         //Load Data
         double binwidth = (5.1-3.4)/double(NData);
         for(int i=0;i<NData;i++){
             Bdata tmp;
             tmp.pdata = Data[i];
             tmp.pdataerror = DataError[i];
             tmp.penergy = 3.4+double(i)*binwidth+0.5*binwidth;
             test_fitter.AddData(tmp);
         }
         
         test_fitter.AddRes("X_4000",4.0,0.1,0.0,1.0);
         test_fitter.AddRes("X_4200",4.2,0.15,1.5,1.0);
         test_fitter.AddRes("X_4600",4.6,0.25,2.5,1.0);
//         test_fitter.Minimize();
     }

     gt = clock_gettime(CLOCK_MONOTONIC, &rup);
     assert(gt != -1);
     long mark_end_s = rup.tv_sec;
     long mark_end_ns = rup.tv_nsec;
     time_minos += ((mark_end_s-mark_start_s)*1000000000+(mark_end_ns-mark_start_ns));
     std::cout<<"Time : "<<time_minos*1.0E-9<<std::endl;

}

int main(int argc, char * argv[]){
     if(argc!=6){std::cout<<"Error.Usage: ./chi2_test number_of_resonances input_file input_tree_name first_entry n_entry" <<std::endl; return -1;}
     const int _nresonance = atoi(argv[1]);
     int entry_idx= atoi(argv[4]);
     int para_idx = atoi(argv[5]);
     NewMinos_Test(string(argv[2]), string(argv[3]), _nresonance, 100,  entry_idx, para_idx);
     return 0;
}
