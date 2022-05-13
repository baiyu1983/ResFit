#include "Header.h"
#include "AnalyticPhi23.h"
#include <algorithm>
#include <iostream>
#include "TMatrixD.h"
#include "CalDerivation.h"

namespace ResFit{ 
/*
   resonance Model::GetResonance(string resname){
        std::map<std::string,resonance>::iterator itr = FindRes(resname);
        if( itr != _vres.end()){
           return (itr->second);
        }else{
           return ;
        }
   }
*/
   void ResFit::Model::AddRes(std::string resname, double mass, double width , double phase, double branchratio){
           resonance tmp; tmp._mass = mass; tmp._width = width; tmp._phase = phase; tmp._branchratio = branchratio;
           this->AddRes(resname, tmp);
   }
   
   void ResFit::Model::AddRes(std::string resname, ResFit::resonance res){
       if((this->_vres).count(resname)){
            std::cout<<"WARNING: Resonance with name "<<resname<<" already exists!"<<std::endl;
            return;
       }else if(resname==""){
            std::cout<<"ERROR: Resonance must have a name!"<<std::endl;
            return;
       }
       else{
/*
            std::cout<<"Adding resonance "<<resname<<std::endl;
            std::cout<<"Mass = "<<res._mass<<std::endl;
            std::cout<<"Width = "<<res._width<<std::endl;
            std::cout<<"Phase = "<<res._phase<<std::endl;
            std::cout<<"BranchRatio = "<<res._branchratio<<std::endl;
*/
            res._pole=std::complex<double>(res._mass*res._mass,-res._width*res._mass);
            std::complex<double> tmp_coeffs(1,0);
            tmp_coeffs *= res._mass;
            tmp_coeffs *= sqrt(res._width*res._branchratio);
            tmp_coeffs *= std::exp(std::complex<double>(0,res._phase));
            tmp_coeffs *= sqrt(1.0/AnalyticPhi23(res._mass));
            res._coeff=tmp_coeffs;
            res._parameterized2 = true;
            (this->_vres).insert(std::pair<std::string,ResFit::resonance>(resname,res));            
            _nres = (this->_vres).size();
       }
   }  

   int ResFit::Model::GetVarIndex (std::string resname, ResFit::respara var) const{
        int counter = 0;
        bool find_resonance =false;
        int result = -1;
        for(std::map<std::string,ResFit::resonance>::const_iterator itr= _vres.begin(); itr !=_vres.end();++itr){
             if(itr->first != resname){
                  counter += 1; continue;
             }else{
                  find_resonance = true;
                  result = 4*counter+var;
                  break;
             }
        }
        if(!find_resonance || (var>=4)){
             std::cout<<"Error::Variable["<<var<<"] not found for resonance "<<resname<<std::endl;
             return -1;
        }
        else return result;
   }

   int ResFit::Model::GetVarIndex_Reduced (std::string resname, ResFit::respara var) const{
        int counter = 0;
        bool find_resonance = false;
        int result = -1;
        if(resname == _fixed_phase_resonance && var == Phase){std::cout<<"Error: Fixed phase angle doesn't have reduced index! resname = "<<resname<<", variable = "<<var<<std::endl;return -1;}
        for(std::map<std::string,ResFit::resonance>::const_iterator itr= _vres.begin(); itr !=_vres.end();++itr){
             if(itr->first != resname){
                  counter += ((itr->first == _fixed_phase_resonance)?3:4); continue;
             }else{
                  find_resonance = true;
                  if(itr->first != _fixed_phase_resonance){
                     result = counter+var;
                     break;
                  }else{
                       result = (var==Branchratio)?counter+var-1:counter+var;
                       break;
                  }
             }
        }
        if(!find_resonance || var>=4){
           std::cout<<"Error::Variable["<<var<<"] not found for resonance "<<resname<<std::endl;
           return -1;
        }
        else return result;
   }


   void ResFit::Model::UpdateRes(std::string resname, ResFit::resonance res){
       if((this->_vres).count(resname) == 0){
            std::cout<<"WARNING: Resonance with name "<<resname<<" not exists!"<<std::endl;
            return;
       }else if(resname==""){
            std::cout<<"ERROR: Resonance must have a name!"<<std::endl;
            return;
       }else{
            (this->_vres)[resname]._mass = res._mass;
            (this->_vres)[resname]._width = res._width;
            (this->_vres)[resname]._phase = res._phase;
            (this->_vres)[resname]._branchratio = res._branchratio;
            (this->_vres)[resname]._pole = std::complex<double>(res._mass*res._mass,-res._width*res._mass);
            std::complex<double> tmp_coeffs(1,0);
            tmp_coeffs *= res._mass;
            tmp_coeffs *= sqrt(res._width*res._branchratio);
            tmp_coeffs *= std::exp(std::complex<double>(0,res._phase));
            tmp_coeffs *= sqrt(1.0/AnalyticPhi23(res._mass));
            (this->_vres)[resname]._coeff = tmp_coeffs;
            (this->_vres)[resname]._parameterized1 = true;
            (this->_vres)[resname]._parameterized2 = true; 

       }
   }

   void ResFit::Model::DelRes(std::string resname){
       (this->_vres).erase(resname);
       _nres = (this->_vres).size();
   }

   void ResFit::Model::CopyRes(ResFit::Model v_source){
       _vres.clear(); _nres=0;
       std::map<std::string, ResFit::resonance>  v_resonance = v_source.GetResonance();
       for(std::map<std::string, ResFit::resonance>::iterator itr= v_resonance.begin(); itr!=v_resonance.end();++itr){
              this->AddRes(itr->first,itr->second);
       }
   }

   void ResFit::Model::FixPhaseAngle(std::string resname){
       std::map<std::string,ResFit::resonance>::iterator itr =  (this->_vres).find(resname);
       if(itr == (this->_vres).end()){
            std::cout<<"WARNING : Resonance "<<resname<<" not founded, the phase angle of the first resonance will be set to 0!"<<std::endl;
            this->FixPhaseAngle();
            return;
       }

       _fixed_phase_resonance = (itr->first);
       double phase_to_fix = (itr->second)._phase;
       for(itr = (this->_vres).begin(); itr!=(this->_vres).end();++itr){
             ((itr->second)._phase) -= phase_to_fix;
             ((itr->second)._coeff) *= std::exp(std::complex<double>(0,-phase_to_fix));
       }   
       _fixed_phase_paridx = this->GetVarIndex(resname,Phase);
   }

   void ResFit::Model::FixPhaseAngle(){
       std::map<std::string,ResFit::resonance>::iterator itr = (this->_vres).begin();
       _fixed_phase_resonance = itr->first;
       double phase_to_fix = (itr->second)._phase;
       for(itr = (this->_vres).begin(); itr!=(this->_vres).end();++itr){
             ((itr->second)._phase) -= phase_to_fix;
             ((itr->second)._coeff) *= std::exp(std::complex<double>(0,-phase_to_fix));
       }
       _fixed_phase_paridx = this->GetVarIndex(_fixed_phase_resonance,Phase);
   }

   void ResFit::Chi2::operator() (Bdata data){
       double s = data.penergy*data.penergy;
       std::complex<double> ampsum(0,0);
       for(std::map<std::string,ResFit::resonance>::iterator itr = (this->_vres).begin();itr!=(this->_vres).end();++itr){
           ampsum += (itr->second)._coeff/(s-(itr->second)._pole);
       }
/*
       double phasespace = AnalyticPhi23(data.penergy);
*/
       double xs = GetXS(ampsum,data.penergy,data.pphasespace);
       _chi2 += TMath::Power((xs-data.pdata)/data.pdataerror,2);
   }

   void ResFit::Chi2_Gradient::ReSet(){
       _chi2 = 0; _chi2_gradient = new TMatrixD(4*(this->_nres),1);
       _chi2_gradient_reduced = new TMatrixD(4*(this->_nres)-1,1);
       for(int i=0;i<4*(this->_nres)-1;i++){
          (*_chi2_gradient)[i][0] = 0;
          (*_chi2_gradient_reduced)[i][0] = 0;
       }
       (*_chi2_gradient)[4*(this->_nres)-1][0] = 0;
       
       _phasespace.clear();
       _phasespace_1stder.clear(); //_phasespace_1stder.resize(0);
       for(std::map<std::string,ResFit::resonance>::const_iterator itr= (this->_vres).begin(); itr!=(this->_vres).end();++itr){
            _phasespace.insert(std::pair<std::string,double>(itr->first,AnalyticPhi23((itr->second)._mass)));
            _phasespace_1stder.insert(std::pair<std::string,double>(itr->first,AnalyticPhi23_derv1((itr->second)._mass)));
       }
   }


   void ResFit::Chi2_Gradient::operator()(Bdata data){
       double s = data.penergy*data.penergy;
/*
       double phasespace = AnalyticPhi23(data.penergy);
*/
       std::complex<double> ampsum(0,0);
/*     //This chech should be done before implementing this functor
       if((this->_vres).size()!=(this->_nres)){
           std::cout<<"ERROR : Resonance size : "<<(this->_vres).size()<<" not equal to nres : "<<this->_nres<<"!"<<std::endl;
           return;
       }
*/
       const int NRES = (this->_nres);
       std::complex<double> amp[NRES] ;

       int counter =0;
       for(std::map<std::string,ResFit::resonance>::iterator itr = (this->_vres).begin();itr!=(this->_vres).end();++itr){
           amp[counter] = (itr->second)._coeff/(s-(itr->second)._pole);
           ampsum +=  amp[counter];
           counter +=1;
       }
       double xs = GetXS(ampsum,data.penergy,data.pphasespace);
       _chi2 += TMath::Power((xs-data.pdata)/data.pdataerror,2);
       double delta = (xs - data.pdata)/(data.pdataerror*data.pdataerror);      
       counter = 0;
       double KCoeff =GetXSConst(data.penergy, data.pphasespace) ;
       bool passfixed = false;
       for(std::map<std::string,ResFit::resonance>::iterator itr = (this->_vres).begin();itr!=(this->_vres).end();++itr){
           std::complex<double> tmp_1stDerv[4]={0,0,0,0};
           for(int i=0;i<4;i++){
                if(this->_fixed_phase_resonance == itr->first && i==2){passfixed = true ; continue;}
                tmp_1stDerv[i] = FunctionFirstDerivative(data.penergy, itr->second, i, (this->_phasespace)[itr->first] ,(this->_phasespace_1stder)[itr->first]);
                (*_chi2_gradient)[4*counter+i][0] += 4*KCoeff*std::real(tmp_1stDerv[i]*amp[counter]*std::conj(ampsum))*delta;
                if(!passfixed){
                     (*_chi2_gradient_reduced)[4*counter+i][0] = (*_chi2_gradient)[4*counter+i][0];
                }else if(this->_fixed_phase_resonance == itr->first){
                     if(i<2)(*_chi2_gradient_reduced)[4*counter+i][0] = (*_chi2_gradient)[4*counter+i][0];
                     else (*_chi2_gradient_reduced)[4*counter+i-1][0] = (*_chi2_gradient)[4*counter+i][0];
                }else{
                     (*_chi2_gradient_reduced)[4*counter+i-1][0] = (*_chi2_gradient)[4*counter+i][0];
                }
           }
           counter+=1;
       }
   }

   void ResFit::Chi2_Hess::ReSet(){
       _chi2 = 0;
       _chi2_gradient = new TMatrixD(4*(this->_nres),1);
       _chi2_gradient_reduced = new TMatrixD(4*(this->_nres)-1,1);
       for(int i=0;i<4*(this->_nres)-1;i++){
          (*_chi2_gradient)[i][0] = 0;
          (*_chi2_gradient_reduced)[i][0] = 0;
       }
       (*_chi2_gradient)[4*(this->_nres)-1][0] = 0;


       _phasespace.clear();
       _phasespace_1stder.clear(); //_phasespace_1stder.resize(0);
       _phasespace_2ndder.clear();
       for(std::map<std::string,ResFit::resonance>::const_iterator itr= (this->_vres).begin(); itr!=(this->_vres).end();++itr){
            _phasespace.insert(std::pair<std::string,double>(itr->first,AnalyticPhi23((itr->second)._mass)));
            _phasespace_1stder.insert(std::pair<std::string,double>(itr->first,AnalyticPhi23_derv1((itr->second)._mass)));
            _phasespace_2ndder.insert(std::pair<std::string,double>(itr->first,AnalyticPhi23_derv2((itr->second)._mass)));
       }

       if(_chi2_hess){
            delete _chi2_hess;
       }
       _chi2_hess = new TMatrixD(4*(this->_vres).size(),4*(this->_vres).size());
   
       for(int i=0;i<4*(this->_vres).size();i++){
             for(int j=0; j<4*(this->_vres).size();j++){
                (*_chi2_hess)[i][j] =0;
             }
       }

       if(_chi2_hess_reduced){
            delete _chi2_hess_reduced;
       }
       _chi2_hess_reduced = new TMatrixD(4*(this->_vres).size()-1,4*(this->_vres).size()-1);
      
       for(int i=0;i<4*(this->_vres).size()-1;i++){
             for(int j=0; j<4*(this->_vres).size()-1;j++){
                (*_chi2_hess_reduced)[i][j] =0;
             }
       }
   }


   void ResFit::Chi2_Hess::operator()(Bdata data){
       double s = data.penergy*data.penergy;
       std::complex<double> ampsum(0,0);
       
       const int NRES = (this->_nres);
       std::complex<double> amp[NRES] ;

       int counter =0;
       for(std::map<std::string,ResFit::resonance>::iterator itr = (this->_vres).begin();itr!=(this->_vres).end();++itr){
           amp[counter] = (itr->second)._coeff/(s-(itr->second)._pole);
           ampsum +=  amp[counter];
           counter +=1;
       }
       double xs = GetXS(ampsum,data.penergy,data.pphasespace);
       _chi2 += TMath::Power((xs-data.pdata)/data.pdataerror,2);

       double delta = (xs - data.pdata)/(data.pdataerror*data.pdataerror);

       counter = 0;
       double KCoeff =GetXSConst(data.penergy, data.pphasespace) ;
       std::complex<double> tmp_1stDerv[4*NRES];
       double xs_1stDerv[4*NRES];
       bool ipass_fixed = false;
       for(std::map<std::string,ResFit::resonance>::iterator itr = (this->_vres).begin();itr!=(this->_vres).end();++itr){
           for(int i=0;i<4;i++){
                if(this->_fixed_phase_resonance == itr->first && i==2){ipass_fixed = true; continue;}
                tmp_1stDerv[counter*4+i] = FunctionFirstDerivative(data.penergy, itr->second, i, (this->_phasespace)[itr->first] ,(this->_phasespace_1stder)[itr->first]);
                xs_1stDerv[counter*4+i] = 2*KCoeff*std::real(tmp_1stDerv[counter*4+i]*amp[counter]*std::conj(ampsum));
                (*_chi2_gradient)[4*counter+i][0] += xs_1stDerv[counter*4+i]*2*delta;
                if(!ipass_fixed){
//                     std::cout<<"DEBUG : (index1,index2) = "<<"("<<4*counter+i<<","<<4*counter+i<<")"<<std::endl;
                     (*_chi2_gradient_reduced)[4*counter+i][0] = (*_chi2_gradient)[4*counter+i][0];
                }else if(this->_fixed_phase_resonance == itr->first){
                     if(i<2){(*_chi2_gradient_reduced)[4*counter+i][0] = (*_chi2_gradient)[4*counter+i][0];
 //                            std::cout<<"DEBUG : (index1,index2) = "<<"("<<4*counter+i<<","<<4*counter+i<<")"<<std::endl;
                     }
                     else{(*_chi2_gradient_reduced)[4*counter+i-1][0] = (*_chi2_gradient)[4*counter+i][0];
//                             std::cout<<"DEBUG : (index1,index2) = "<<"("<<4*counter+i-1<<","<<4*counter+i<<")"<<std::endl;
                     }
                }else{
                     (*_chi2_gradient_reduced)[4*counter+i-1][0] = (*_chi2_gradient)[4*counter+i][0];
//                     std::cout<<"DEBUG : (index1,index2) = "<<"("<<4*counter+i-1<<","<<4*counter+i<<")"<<std::endl;
                }
           }
           counter+=1;
       }
//       std::cout<<"DEBUG: gradient calculated"<<std::endl;
       counter =0;
       int icounter =0; 
       ipass_fixed = false;
       for(std::map<std::string,ResFit::resonance>::iterator itr = (this->_vres).begin();itr!=(this->_vres).end();++itr){
           for(int i=0;i<4;i++){
                 if((!ipass_fixed) && i==2){ //Check if the phase angle need to be fit in this resonance i
                     if(itr->first == (this->_fixed_phase_resonance))ipass_fixed = true;
                 }
                 int jcounter = 0; int idx = 4*icounter+i; int idx_reduced = ipass_fixed?(idx-1):idx;
                 bool jpass_fixed = false;
                 for(std::map<std::string,ResFit::resonance>::iterator jtr = (this->_vres).begin();jtr!=(this->_vres).end();++jtr){
                      for(int j=0;j<4;j++){
                           if((!jpass_fixed) && j==2){  //Check if the phase angle need to be fit in this resonance j
                                 if(jtr->first == (this->_fixed_phase_resonance))jpass_fixed = true;
                           }    
                           int jdx = 4*jcounter+j;  int jdx_reduced = jpass_fixed?(jdx-1):jdx;
                           if(jdx<idx)continue;
                           if((itr->first == (this->_fixed_phase_resonance) && i==2) || (jtr->first == (this->_fixed_phase_resonance) && j==2)  )continue;
                           double tmp_elem = 0;
                           tmp_elem += 2*std::real(std::conj(tmp_1stDerv[idx]*amp[icounter])*(tmp_1stDerv[jdx]*amp[jcounter])); //Product of first derv
                           if(icounter == jcounter){            //Adding 2nd Derv
                                 tmp_elem += 2*std::real(std::conj(ampsum)*amp[icounter]*(tmp_1stDerv[idx]*tmp_1stDerv[jdx]+FunctionSecondDerivative(data.penergy,itr->second,i,j,(this->_phasespace)[itr->first] ,(this->_phasespace_1stder)[itr->first],(this->_phasespace_2ndder)[itr->first] )));
                           }
                           tmp_elem *= (2*delta);
                           tmp_elem *= KCoeff;
                           tmp_elem += 2*xs_1stDerv[idx]*xs_1stDerv[jdx]/(data.pdataerror*data.pdataerror);  
                           (*_chi2_hess)[idx][jdx] += tmp_elem;
                           if(idx!=jdx)
                               (*_chi2_hess)[jdx][idx] += tmp_elem;
//                           if((itr->first == (this->_fixed_phase_resonance) && i==2) || (jtr->first == (this->_fixed_phase_resonance) && j==2)  )continue;
                           (*_chi2_hess_reduced)[idx_reduced][jdx_reduced] += tmp_elem;
                           if(idx_reduced != jdx_reduced)
                               (*_chi2_hess_reduced)[jdx_reduced][idx_reduced] += tmp_elem;
                      }//End of j , value of 2nd resonance
                      jcounter += 1;
                 }//End of jtr(2nd resonace)
           } //End of i, value of 1st resonance
           icounter+=1;
       }//Endf of itr(1st resonance)

   }//Endf of Definition ResFit::Chi2_Hess


   void ResFit::Chi2_Hess::SetCovariance(){
       int NRES=this->_nres; 
       _chi2_covariance = new TMatrixD(4*NRES,4*NRES);
       _chi2_covariance_reduced = new TMatrixD(4*NRES-1,4*NRES-1);
       for(int i=0;i<4*NRES-1;i++){
           for(int j=0;j<4*NRES-1;j++){
               (*_chi2_covariance_reduced)[i][j] = (*_chi2_hess_reduced)[i][j]*0.5;
           }
       }
       _chi2_covariance_reduced->Invert();
       for(int i=0;i<4*NRES;i++){ 
            for(int j=0;j<4*NRES;j++){
               if(i==_fixed_phase_paridx || j==_fixed_phase_paridx){
                    (*_chi2_covariance)[i][j] = 0;continue;
               }else{
                    int ii = (i<_fixed_phase_paridx)?i:i-1;
                    int jj = (j<_fixed_phase_paridx)?j:j-1;
/*
                    std::cout<<"NRES = "<<NRES<<", fixed phase idx = "<<_fixed_phase_paridx<<std::endl;
                    std::cout<<"DEBUG: (i,j) = ("<<i<<", "<<j<<"), (ii,jj) = ("<<ii<<", "<<jj<<")"<<std::endl;
*/
                    (*_chi2_covariance)[i][j]  = (*_chi2_covariance_reduced)[ii][jj] ;
               }
            }
       }      
   }
} //End of ResFit
