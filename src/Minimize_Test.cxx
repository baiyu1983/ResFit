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

void Minimize_Test(string inputfile, string treename,  int NRES, int NData, Long64_t jentry, Long64_t NEvents){
     int NPar = NRES*4;
     int NSol = 1;
     for(int i=0;i<NRES;i++){
           NSol = NSol*2;
     }

     NSol /= 2;
//     int ipar = (iipar>2)?iipar-1:iipar;
     double Data[NData];
     double DataError[NData];

     TFile finput(inputfile.c_str(),"read");
     TTree *tinput = (TTree*)finput.Get(treename.c_str());
   
     tinput->SetBranchStatus("Chi2", 0);
     tinput->SetBranchStatus("Value",0);
     tinput->SetBranchStatus("Error",0);
     tinput->SetBranchStatus("ErrorStatus",0);
     tinput->SetBranchAddress("Data",Data);
     tinput->SetBranchAddress("DataError",DataError);
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

         test_fitter.FixPhaseAngle();
/*
         test_fitter.ReSet();
         test_fitter.UpdateHess();
         std::cout<<"Chi2 = "<<test_fitter.GetChi2Min()<<std::endl;
*/
         test_fitter.Minimize(); 
//         std::cout<<"Chi2 = "<<test_fitter.GetChi2Min()<<std::endl;
         TMatrixD * covariance = test_fitter.GetCovariance_Reduced();
         TMatrixD * gradient = test_fitter.GetGradient_Reduced();
         std::cout<<"M(X4000) :  "<<test_fitter.GetParameter("X_4000",ResFit::Mass)<<std::endl;
         std::cout<<"W(X4000) :  "<<test_fitter.GetParameter("X_4000",ResFit::Width)<<std::endl;
         std::cout<<"P(X4000) :  "<<test_fitter.GetParameter("X_4000",ResFit::Phase)<<std::endl;
         std::cout<<"B(X4000) :  "<<test_fitter.GetParameter("X_4000",ResFit::Branchratio)<<std::endl;
         std::cout<<"M(X4200) :  "<<test_fitter.GetParameter("X_4200",ResFit::Mass)<<std::endl;
         std::cout<<"W(X4200) :  "<<test_fitter.GetParameter("X_4200",ResFit::Width)<<std::endl;
         std::cout<<"P(X4200) :  "<<test_fitter.GetParameter("X_4200",ResFit::Phase)<<std::endl;
         std::cout<<"B(X4200) :  "<<test_fitter.GetParameter("X_4200",ResFit::Branchratio)<<std::endl;
         std::cout<<"M(X4600) :  "<<test_fitter.GetParameter("X_4600",ResFit::Mass)<<std::endl;
         std::cout<<"W(X4600) :  "<<test_fitter.GetParameter("X_4600",ResFit::Width)<<std::endl;
         std::cout<<"P(X4600) :  "<<test_fitter.GetParameter("X_4600",ResFit::Phase)<<std::endl;
         std::cout<<"B(X4600) :  "<<test_fitter.GetParameter("X_4600",ResFit::Branchratio)<<std::endl;
 
         for(int i = 0; i<4*NRES-1;i++){
               std::cout<<"Gradient["<<i<<"]  = " <<(*gradient)[i][0]<<std::endl;
         }
     
         for(int i =0;i<4*NRES-1;i++){
              for(int j=0;j<4*NRES-1;j++){
                  std::cout<<(*covariance)[i][j]<<"\t";
              }std::cout<<std::endl;
         }
     }

     gt = clock_gettime(CLOCK_MONOTONIC, &rup);
     assert(gt != -1);
     long mark_end_s = rup.tv_sec;
     long mark_end_ns = rup.tv_nsec;
     time_minos += ((mark_end_s-mark_start_s)*1000000000+(mark_end_ns-mark_start_ns));
/*
     std::cout<<"Chi2 = "<<test_fitter.GetChi2Min()<<std::endl;
     TMatrxiD * covariance = test_fitter.GetCovariance_Reduced();
     TMatrxiD * graident = test_fitter.GetGetGradient_Reduced();
     for(int i = 0; i<4*NRES-1;i++){
           std::cout<<"Gradient["<<i<<"] " = (*graident)[i][0];
     }

     for(int i =0;i<4*NRE1-1;i++){
          for(int j=0;j<4*NRES-1;j++){
              std::cout<<(*covariance)[i][j]<<"\t";
          }std::cout<<std::endl;
     }
*/
     std::cout<<"Time : "<<time_minos*1.0E-9<<std::endl;

}

int main(int argc, char * argv[]){
     if(argc!=6){std::cout<<"Error.Usage: ./chi2_test number_of_resonances input_file input_tree_name first_entry n_entry" <<std::endl; return -1;}
     const int _nresonance = atoi(argv[1]);
     Long64_t entry_idx= atoi(argv[4]);
     Long64_t NEvents = atoi(argv[5]);
//     Minimize_Test(string(argv[2]), string(argv[3]), _nresonance, 100,  entry_idx, para_idx);
     Minimize_Test( std::string(argv[2]),  string(argv[3]), _nresonance,100 ,entry_idx, NEvents);
     return 0;
}
