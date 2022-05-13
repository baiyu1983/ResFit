#include "TTree.h"
#include "TFile.h"
#include "TLine.h"
#include "TAxis.h"
#include <string>
#include <iostream>
#include <time.h>
#include "TMatrixD.h"
#include "Header.h"
#include "Contour.h"
#include "TGraph.h"
#include "TCanvas.h"
//#include "SetMatrix.h"
#include "CalDerivation.h"
#define EPSILON 1.0E-4
#define DELTA 1.0

using namespace std;
using namespace ResFit;

void Contour_Test(string inputfile, string treename,  int NRES, int NData, Long64_t jentry, string name_x, string name_y, respara val_x, respara val_y, int np){
     int NPar = NRES*4;
     int NSol = 1;
     for(int i=0;i<NRES;i++){
           NSol = NSol*2;
     }

     NSol /= 2;
//     int ipar = (iipar>2)?iipar-1:iipar;

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

     double binwidth = (5.1-3.4)/double(NData);
     for(int i=0;i<NData;i++){
          Bdata tmp;
          tmp.pdata = Data[i];
          tmp.pdataerror = DataError[i];
          tmp.penergy = 3.4+double(i)*binwidth+0.5*binwidth;
          test_fitter.AddData(tmp);
     }

     double x_central[1] ={ test_fitter.GetParameter(name_x,val_x)};
     double y_central[1] ={ test_fitter.GetParameter(name_y,val_y)};
     std::cout<<"x_central = "<<x_central[0]<<", y_central = "<<y_central[0]<<std::endl;

     Contour cont_test_05;
     Contour cont_ell_05;
     cont_ell_05.SetName("cont_ell_05");
     cont_test_05.SetName("cont_test_05");

     Contour cont_test_1;
     Contour cont_ell_1;
     cont_ell_1.SetName("cont_ell_1");
     cont_test_1.SetName("cont_test_1");
     
     Contour cont_test_15;
     Contour cont_ell_15;
     cont_ell_15.SetName("cont_ell_15");
     cont_test_15.SetName("cont_test_15");

     Contour cont_test_20;
     Contour cont_ell_20;
     cont_ell_20.SetName("cont_ell_20");
     cont_test_20.SetName("cont_test_20");

     test_fitter.FixPhaseAngle();
     test_fitter.ReSet();
     test_fitter.UpdateHess(); 
//     Chi2_Hess par_central = test_fitter.GetCentral();


     test_fitter.FillContour(cont_ell_05,name_x,name_y,Phase,Branchratio,np-1,false,0.25);
     test_fitter.FillContour(cont_test_05,name_x,name_y,Phase,Branchratio,np-1,true,0.25);

     test_fitter.ReSet();
     test_fitter.UpdateHess(); 

     test_fitter.FillContour(cont_ell_1,name_x,name_y,Phase,Branchratio,np-1,false,1);
     test_fitter.FillContour(cont_test_1,name_x,name_y,Phase,Branchratio,np-1,true,1);

     test_fitter.ReSet();
     test_fitter.UpdateHess();
  
     test_fitter.FillContour(cont_ell_15,name_x,name_y,Phase,Branchratio,np-1,false,2.25);
     test_fitter.FillContour(cont_test_15,name_x,name_y,Phase,Branchratio,np-1,true,2.25);

     test_fitter.ReSet();
     test_fitter.UpdateHess();

     test_fitter.FillContour(cont_ell_20,name_x,name_y,Phase,Branchratio,np-1,false,4);
     test_fitter.FillContour(cont_test_20,name_x,name_y,Phase,Branchratio,np-1,true,4);

     std::map<int,std::pair<double,double> > cont_data_05 = cont_test_05.GetData();
     std::map<int,std::pair<double,double> > cont_data_ell_05 = cont_ell_05.GetData();
     std::map<int,std::pair<double,double> > cont_data_1 = cont_test_1.GetData();
     std::map<int,std::pair<double,double> > cont_data_ell_1 = cont_ell_1.GetData();
     std::map<int,std::pair<double,double> > cont_data_15 = cont_test_15.GetData();
     std::map<int,std::pair<double,double> > cont_data_ell_15 = cont_ell_15.GetData();
     std::map<int,std::pair<double,double> > cont_data_20 = cont_test_20.GetData();
     std::map<int,std::pair<double,double> > cont_data_ell_20 = cont_ell_20.GetData();

     double x_05[np+1],y_05[np+1];
     double x_ell_05[np+1], y_ell_05[np+1];
     double x_10[np+1],y_10[np+1];
     double x_ell_10[np+1], y_ell_10[np+1];
     double x_15[np+1],y_15[np+1];
     double x_ell_15[np+1], y_ell_15[np+1];
     double x_20[np+1],y_20[np+1];
     double x_ell_20[np+1], y_ell_20[np+1];

     int counter =0;
     for(std::map<int,std::pair<double,double> >::const_iterator itr = cont_data_05.begin();itr!=cont_data_05.end();++itr){
//          std::cout<<itr->first<<"\t"<<(itr->second).first-x_central<<"\t"<<(itr->second).second-y_central<<std::endl;
          x_05[counter] = (itr->second).first; y_05[counter] = (itr->second).second; counter+=1;
     }    
     counter = 0;

     for(std::map<int,std::pair<double,double> >::const_iterator itr = cont_data_ell_05.begin();itr!=cont_data_ell_05.end();++itr){
//          std::cout<<itr->first<<"\t"<<(itr->second).first-x_central<<"\t"<<(itr->second).second-y_central<<std::endl;
          x_ell_05[counter] = (itr->second).first; y_ell_05[counter] = (itr->second).second; counter+=1;
     }  
     counter =0;
     for(std::map<int,std::pair<double,double> >::const_iterator itr = cont_data_1.begin();itr!=cont_data_1.end();++itr){
//          std::cout<<itr->first<<"\t"<<(itr->second).first-x_central<<"\t"<<(itr->second).second-y_central<<std::endl;
          x_10[counter] = (itr->second).first; y_10[counter] = (itr->second).second; counter+=1;
     }
     counter = 0;
     
     for(std::map<int,std::pair<double,double> >::const_iterator itr = cont_data_ell_1.begin();itr!=cont_data_ell_1.end();++itr){
//          std::cout<<itr->first<<"\t"<<(itr->second).first-x_central<<"\t"<<(itr->second).second-y_central<<std::endl;
          x_ell_10[counter] = (itr->second).first; y_ell_10[counter] = (itr->second).second; counter+=1;
     }
     counter = 0;
     for(std::map<int,std::pair<double,double> >::const_iterator itr = cont_data_15.begin();itr!=cont_data_15.end();++itr){
//          std::cout<<itr->first<<"\t"<<(itr->second).first-x_central<<"\t"<<(itr->second).second-y_central<<std::endl;
          x_15[counter] = (itr->second).first; y_15[counter] = (itr->second).second; counter+=1;
     }
     counter = 0;
     
     for(std::map<int,std::pair<double,double> >::const_iterator itr = cont_data_ell_15.begin();itr!=cont_data_ell_15.end();++itr){
//          std::cout<<itr->first<<"\t"<<(itr->second).first-x_central<<"\t"<<(itr->second).second-y_central<<std::endl;
          x_ell_15[counter] = (itr->second).first; y_ell_15[counter] = (itr->second).second; counter+=1;
     }
     counter =0 ;
     for(std::map<int,std::pair<double,double> >::const_iterator itr = cont_data_20.begin();itr!=cont_data_20.end();++itr){
//          std::cout<<itr->first<<"\t"<<(itr->second).first-x_central<<"\t"<<(itr->second).second-y_central<<std::endl;
          x_20[counter] = (itr->second).first; y_20[counter] = (itr->second).second; counter+=1;
     }
     counter = 0;
     
     for(std::map<int,std::pair<double,double> >::const_iterator itr = cont_data_ell_20.begin();itr!=cont_data_ell_20.end();++itr){
//          std::cout<<itr->first<<"\t"<<(itr->second).first-x_central<<"\t"<<(itr->second).second-y_central<<std::endl;
          x_ell_20[counter] = (itr->second).first; y_ell_20[counter] = (itr->second).second; counter+=1;
     }

     TGraph * g_cont_05 = new TGraph(np+1,x_05,y_05);
     TGraph * g_cont_ell_05 = new TGraph(np+1,x_ell_05,y_ell_05);

     TGraph * g_cont_10 = new TGraph(np+1,x_10,y_10);
     TGraph * g_cont_ell_10 = new TGraph(np+1,x_ell_10,y_ell_10);
     
     TGraph * g_cont_15 = new TGraph(np+1,x_15,y_15);
     TGraph * g_cont_ell_15 = new TGraph(np+1,x_ell_15,y_ell_15);

     TGraph * g_cont_20 = new TGraph(np+1,x_20,y_20);
     TGraph * g_cont_ell_20 = new TGraph(np+1,x_ell_20,y_ell_20);

     TCanvas *c1 = new TCanvas("c1","",800,600);
     c1->cd();
     g_cont_20->SetTitle("Contour plot of two variables");
     g_cont_20->GetXaxis()->SetTitle("Resonance 2 Phase Angle [Rad]");
     g_cont_20->GetYaxis()->SetTitle("Resonance 2 #Gamma_ee #times Br[eV]");
     g_cont_20->Draw("AL");

     g_cont_ell_20->SetLineColor(kRed);
     g_cont_ell_20->SetLineStyle(2);
     g_cont_ell_20->Draw("L same");

     g_cont_10->Draw("L same");
     
     g_cont_ell_10->SetLineColor(kRed);
     g_cont_ell_10->SetLineStyle(2);
     g_cont_ell_10->Draw("L same");

     g_cont_15->Draw("L same");

     g_cont_ell_15->SetLineColor(kRed);
     g_cont_ell_15->SetLineStyle(2);
     g_cont_ell_15->Draw("L same");

     g_cont_05->Draw("L same");

     g_cont_ell_05->SetLineColor(kRed);
     g_cont_ell_05->SetLineStyle(2);
     g_cont_ell_05->Draw("L same");

     TGraph * g_center = new TGraph(1,x_central,y_central);
     g_center->SetMarkerStyle(43);
     g_center->SetMarkerSize(2.0);
     g_center->Draw("P same");

     TLine * line_20_plus = new TLine(x_20[0],0.62,x_20[0],1.22); 
     TLine * line_20_mins = new TLine(x_20[np/2],0.62,x_20[np/2],1.22);
     TLine * line_ell_20_plus = new TLine(x_ell_20[0],0.62,x_ell_20[0],1.22);
     TLine * line_ell_20_mins = new TLine(x_ell_20[np/2],0.62,x_ell_20[np/2],1.22);
     
     line_ell_20_plus->SetLineColor(2); line_ell_20_plus->SetLineColor(kRed);
     line_ell_20_mins->SetLineColor(2); line_ell_20_mins->SetLineColor(kRed);

     TLine * line_20_down = new TLine(1.092,y_20[np/4],2.1,y_20[np/4]);
     TLine * line_20_up = new TLine(1.092,y_20[3*np/4],2.1,y_20[3*np/4]);

     TLine * line_ell_20_down = new TLine(1.092,y_ell_20[np/4],2.1,y_ell_20[np/4]);
     TLine * line_ell_20_up = new TLine(1.092,y_ell_20[3*np/4],2.1,y_ell_20[3*np/4]);
     line_ell_20_down->SetLineColor(2); line_ell_20_down->SetLineColor(kRed);
     line_ell_20_up->SetLineColor(2); line_ell_20_up->SetLineColor(kRed);
//     line_20_plus->DrawLine(x_20[0],6.5,x_20[0],12.5);
     line_20_plus->Draw("L same");
     line_20_mins->Draw("L same");
     line_ell_20_plus->Draw("L same");
     line_ell_20_mins->Draw("L same");
     line_20_down->Draw("L same");
     line_20_up->Draw("L same");
     line_ell_20_down->Draw("L same");
     line_ell_20_up->Draw("L same");
     c1->Print("cont_test.pdf");
}


int main(int argc, char * argv[]){
     if(argc!=10){std::cout<<"Error.Usage: ./cont_test number_of_resonances input_file input_tree_name entry_idx name_x name_y val_x val_y np" <<std::endl; return -1;}
     const int _nresonance = atoi(argv[1]);
     int entry_idx= atoi(argv[4]);
     respara val_x = ResFit::respara(atoi(argv[7]));
     respara val_y = ResFit::respara(atoi(argv[8]));
     int np = atoi(argv[9]);
     std::string name_x = std::string(argv[5]);
     std::string name_y = std::string(argv[6]);
     Contour_Test(string(argv[2]), string(argv[3]), _nresonance, 100,  entry_idx, name_x, name_y, val_x, val_y, np);
     return 0;
}
