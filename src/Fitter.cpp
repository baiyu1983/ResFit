#include "Header.h"
#include "Contour.h"
#include <iostream>
#include "AnalyticPhi23.h"

namespace ResFit{
/*
    Chi2_Fitter::~Chi2_Fitter(){
         _para_central.~Chi2_Hess();
         _contour_x.clear(); _contour_y.clear();
         _para_limits.clear();
    }
*/
void Chi2_Fitter::AddData(Bdata data){
     double phasespace = AnalyticPhi23(data.penergy);
     Bdata tmp_data;
     tmp_data.penergy = data.penergy;
     tmp_data.pdata = data.pdata;
     tmp_data.pdataerror = data.pdataerror;
     tmp_data.pphasespace = phasespace;
     _data.push_back(tmp_data);
}

void Chi2_Fitter::AddRes(std::string resname, double m, double w, double p,  double b){
     resonance tmp; tmp._mass = m; tmp._width = w; tmp._phase = p; tmp._branchratio = b;
     this->AddRes(resname, tmp);
}

void Chi2_Fitter::FixPhaseAngle(std::string resname){
     _para_central.FixPhaseAngle(resname);
     if(_para_limits.empty())return;
     for(std::map<std::pair<std::string,ResFit::respara >, ResFit::Chi2_Hess >::iterator itr = _para_limits.begin(); itr!=_para_limits.end();++itr){
           (itr->second).FixPhaseAngle(resname);
     }
}

void Chi2_Fitter::FixPhaseAngle(){
     _para_central.FixPhaseAngle();
     if(_para_limits.empty())return;
     for(std::map<std::pair<std::string,ResFit::respara >, ResFit::Chi2_Hess >::iterator itr = _para_limits.begin(); itr!=_para_limits.end();++itr){
           (itr->second).FixPhaseAngle();
     }
}

void Chi2_Fitter::ReSet(){
    _para_central.ReSet();
    _para_limits.clear(); 
}

void Chi2_Fitter::Minimize() {
    int zerocounter;  //counts the number of accept zero of Chi2_Gradient
    int const_mu = 2;
    int NRES = _para_central.GetNRes();
    double lambda = _lambda;
    this->ReSet();
    this->UpdateHess();
    TMatrixD * chi2_gradient_reduced = _para_central.GetChi2Gradient_Reduced();
/*
    std::map<std::string, ResFit::resonance> res_init = _para_central.GetResonance();
    for(std::map<std::string, ResFit::resonance>::iterator itr = res_init.begin(); itr!=res_init.end();++itr){
           std::cout<<itr->first<<"\t"<<(itr->second)[Mass]<<"\t"<<(itr->second)[Width]<<"\t"<<(itr->second)[Phase]<<"\t"<<(itr->second)[Branchratio]<<std::endl;
    }
*/    
    for (int circounter = 0; circounter <= 10; circounter++) {  //Converge to the desired precison within 100 times.
        chi2_gradient_reduced = _para_central.GetChi2Gradient_Reduced();
        for (zerocounter = 0; zerocounter < 4 * NRES - 1; zerocounter++)
            if ((isnan((*chi2_gradient_reduced)[zerocounter][0]) || isinf((*chi2_gradient_reduced)[zerocounter][0])) == 1) return;        //Failure or error to converge.
        for (zerocounter = 0; zerocounter < 4 * NRES - 1; zerocounter++)
            if (abs((*chi2_gradient_reduced)[zerocounter][0]) > 10E-5) break;
        if (zerocounter == 4 * NRES - 1) return;    //Sufficient accuracy!
        ResFit::Chi2 tmp_chi2_lamda, tmp_chi2_lamda_mu;
        ParameterUpdateCal(tmp_chi2_lamda, lambda);
        ParameterUpdateCal(tmp_chi2_lamda_mu, lambda / const_mu);
        std::map<std::string,ResFit::resonance> v_res_tmp_chi2_lamda_mu = tmp_chi2_lamda_mu.GetResonance();
        std::map<std::string,ResFit::resonance> v_res_tmp_chi2_lamda = tmp_chi2_lamda.GetResonance();
        if (tmp_chi2_lamda.GetChi2() <= _para_central.GetChi2()) {
            if (tmp_chi2_lamda_mu.GetChi2() <= _para_central.GetChi2()) {
                UpdateParameter(v_res_tmp_chi2_lamda_mu);
                lambda = lambda / const_mu;
            }
            else UpdateParameter(v_res_tmp_chi2_lamda);
        }
        else lambda = lambda * const_mu;
    }
//    std::for_each(_data.begin(),_data.end(), std::ref(_para_central));
//    this->UpdateHess();
}

void Chi2_Fitter::UpdateParameter(std::map<std::string, ResFit::resonance>  changed_res) {
    int NRES = changed_res.size();
/*
    for (int i = 0; i < NRES; i++) {
        ResFit::resonance tmp;
        std::string resname("res_");
        resname.append(std::to_string(i));
        _para_central.UpdateRes(resname, changed_res[resname]);
    }
*/
    for(std::map<std::string, ResFit::resonance>::iterator itr = changed_res.begin(); itr!=changed_res.end();++itr){
        _para_central.UpdateRes(itr->first, changed_res[itr->first]);
    }

//    _para_central.FixPhaseAngle();
/*
    _para_central.ReSet();
    std::for_each(_data.begin(),_data.end(),std::ref(_para_central));
    _para_central.SetCovariance();
*/
    this->UpdateHess();
}

void Chi2_Fitter::ParameterUpdateCal(ResFit::Chi2 & _tmp_chi2, double lamda_mu) {
    int NRES = _para_central.GetNRes();
    TMatrixD calculateMatrix(4 * NRES - 1, 4 * NRES - 1);
    TMatrixD * chi2_hess_reduced =_para_central.GetHess_Reduced() ;
    TMatrixD * chi2_gradient_reduced = _para_central.GetChi2Gradient_Reduced();
    for (int i = 0; i < 4 * NRES - 1; i++) {
        for (int j = 0; j < 4 * NRES - 1; j++) {
            calculateMatrix[i][j] = (*chi2_hess_reduced)[i][j] + lamda_mu * (i == j);      //i==j stand for IdentityMatrix
        }
    }
    calculateMatrix.Invert();

    TMatrixD parameter_matrix(4 * NRES, 1);
    TMatrixD tmp_parameter_matrix = calculateMatrix * (*chi2_gradient_reduced);
    for (int i = 0; i < 4 * NRES; i++) {
        if (i == 2) parameter_matrix[2][0] = 0.;
        else parameter_matrix[i][0] = tmp_parameter_matrix[i - (i >= 2)][0];
    }

    std::map<std::string, ResFit::resonance> _variable_res = _para_central.GetResonance();
/*
    for (int i = 0; i < NRES; i++) {
        ResFit::resonance tmp;
        std::string resname("res_");
        resname.append(std::to_string(i));
        tmp._mass = _variable_res[resname]._mass - parameter_matrix[4 * i][0];
        tmp._width = _variable_res[resname]._width - parameter_matrix[4 * i + 1][0];
        tmp._phase = _variable_res[resname]._phase - parameter_matrix[4 * i + 2][0];
        tmp._branchratio = _variable_res[resname]._branchratio - parameter_matrix[4 * i + 3][0];
        _tmp_chi2.AddRes(resname, tmp);
    }
*/
    int res_counter = 0;
    for(std::map<std::string, ResFit::resonance>::iterator itr = _variable_res.begin(); itr!=_variable_res.end();++itr){
        ResFit::resonance tmp;
        tmp._mass = _variable_res[itr->first]._mass - parameter_matrix[4 * res_counter][0];
        tmp._width = _variable_res[itr->first]._width - parameter_matrix[4 * res_counter + 1][0];
        tmp._phase = _variable_res[itr->first]._phase - parameter_matrix[4 * res_counter + 2][0];
        tmp._branchratio = _variable_res[itr->first]._branchratio - parameter_matrix[4 * res_counter + 3][0];
        _tmp_chi2.AddRes(itr->first,tmp);
        res_counter+=1;
    }
    _tmp_chi2.FixPhaseAngle(_para_central.GetFixedPhaseRes());
    _tmp_chi2.ReSet();
    std::for_each(_data.begin(), _data.end(), std::ref(_tmp_chi2));
}

void Chi2_Fitter::UpdateHess(){
    _para_central.ReSet();
    for_each(_data.begin(),_data.end(),std::ref(_para_central));
//    std::cout<<"Hessina Matrix Calculated "<<std::endl;
    _para_central.SetCovariance();
}
/*
    void InitalMinos(std::string resname_y, ResFit::respara variable_y, std::string resname_x, ResFit::respara variable_x, double x_val, double Delta, bool updatehess){   //x is fixed, initial for variable y
        if(updatehess)
           this->UpdateHess();
        TMatrixD * hess_reduced = v_res_central.GetHess_Reduced();
        hess_reduced->Invert();
        int NRES =_para_central.GetNRes() ;
        double delta0[1][2];
        double kelta0[1][2];

        int var_idx_reduced =  v_res_central.GetVarIndex_Reduced(resname, variable);
        
    }
*/
void Chi2_Fitter::InitialMinos(std::string resname, ResFit::respara variable, bool updatehessstd, double Delta){
    if(updatehessstd)
        this->UpdateHess();
    std::map<std::string,ResFit::resonance> v_res_central = _para_central.GetResonance();
    TMatrixD * hess_reduced = _para_central.GetCovariance_Reduced();
//        hess_reduced->Invert();//Need to be re-inverted later
    int fixed_phase_index = _para_central.GetVarIndex(_para_central.GetFixedPhaseRes(),ResFit::Phase);
    int NRES =_para_central.GetNRes() ;
    double delta0[1][2];
    double kelta0[1][2]; 

    int var_idx_reduced =  _para_central.GetVarIndex_Reduced(resname, variable);
    int var_idx =  _para_central.GetVarIndex(resname, variable);
    if(var_idx_reduced == -1){
         std::cout<<"Initialize Minos Error : variable index abnormal"<<std::endl;
//             hess_reduced->Invert();
         return;
    }

    delta0[0][0] = sqrt(Delta*((*hess_reduced)[var_idx_reduced][var_idx_reduced]));
    delta0[0][1] = -sqrt(Delta*((*hess_reduced)[var_idx_reduced][var_idx_reduced]));
    kelta0[0][0] = 2.0*Delta/delta0[0][0];
    kelta0[0][1] = 2.0*Delta/delta0[0][1];
/*
        if(resname == _para_central.GetFixedPhaseRes()){
            if(variable == Phase){delta0[0][0] = 0; delta0[0][1] = 0; kelta0[0][0] = 0; kelta0[0][1] = 0;}
        }else{
            ResFit::resonance v_res_central = _para_central.GetResonance();
            int icounter = 0; bool passed_fixed = false;
            for(std::map<std::string, ResFit::resonance>::const_iterator itr = v_res_central.begin(); itr!=v_res_central.end();++itr){
                 if(itr == _para_central.GetFixedPhaseRes() && variable == ResFit::Phase){
                      delta0[0][0] = 0; delta0[0][1] = 0; kelta0[0][0] = 0; kelta0[0][1] = 0;
                      break;
                 }
                 else{
                     if(resname == itr->first ){                 
                        delta0[0][0] = sqrt(2*Delta*((*hess_reduced)[ii][ii]));
                        delta0[0][1] = -sqrt(2*Delta*((*hess_reduced)[ii][ii]));
                        kelta0[0][0] = 2.0*Delta/delta0[i][0];
                        kelta0[0][1] = 2.0*Delta/delta0[i][1];
                     }
                 }
            }    
        }
*/
/*
    std::cout<<"Delta = "<<delta0[0][0]<<std::endl;
    std::cout<<"Kelta = "<<kelta0[0][0]<<std::endl;
    std::cout<<"Fixed phase index = "<<fixed_phase_index<<std::endl;
*/
    TMatrixD * KM0_up_Redu = new TMatrixD(4*NRES-1,1);
    TMatrixD * KM0_down_Redu = new TMatrixD(4*NRES-1,1);
    
    for(int i=0;i<4*NRES-1;i++){
        int ii = (i<fixed_phase_index)?i:i+1;
        (*KM0_up_Redu)[i][0] = (ii==var_idx)?kelta0[0][0]:0;
        (*KM0_down_Redu)[i][0] = (ii==var_idx)?kelta0[0][1]:0;
    }

    TMatrixD * DM0_up_Redu = new TMatrixD(4*NRES-1,1);
    TMatrixD * DM0_down_Redu = new TMatrixD(4*NRES-1,1);
    (*DM0_up_Redu) = 0.5*(*hess_reduced)*(*KM0_up_Redu);
    (*DM0_down_Redu) = 0.5*(*hess_reduced)*(*KM0_down_Redu);    

    Chi2_Hess chi2_hess_up;
    Chi2_Hess chi2_hess_down;
/*
    for(int i=0;i<4*NRES-1;i++){
           std::cout<<"DEBUG::KMatrix_up["<<i<<"][0] = "<<(*KM0_up_Redu)[i][0]<<std::endl;
    } 

    for(int i=0;i<4*NRES-1;i++){
           std::cout<<"DEBUG::DMatrix_up["<<i<<"][0] = "<<(*DM0_up_Redu)[i][0]<<std::endl;
    }

    for(int i=0;i<4*NRES-1;i++){
           std::cout<<"DEBUG::KMatrix_down["<<i<<"][0] = "<<(*KM0_down_Redu)[i][0]<<std::endl;
    }

    for(int i=0;i<4*NRES-1;i++){
           std::cout<<"DEBUG::DMatrix_down["<<i<<"][0] = "<<(*DM0_down_Redu)[i][0]<<std::endl;
    }
*/
    bool passfix = false;
    int icounter = 0;
    for(std::map<std::string, ResFit::resonance>::iterator itr = v_res_central.begin(); itr!=v_res_central.end();++itr){
          resonance tmp_up;
          resonance tmp_down;
          if(itr->first == _para_central.GetFixedPhaseRes()){ 
                tmp_up._phase =0;
                tmp_up._mass = (itr->second)[Mass]+(*DM0_up_Redu)[4*icounter][0];
                tmp_up._width = (itr->second)[Width]+(*DM0_up_Redu)[4*icounter+1][0];
                tmp_up._branchratio = (itr->second)[Branchratio]+(*DM0_up_Redu)[4*icounter+2][0];
                tmp_down._mass = (itr->second)[Mass]+(*DM0_down_Redu)[4*icounter][0];
                tmp_down._width = (itr->second)[Width]+(*DM0_down_Redu)[4*icounter+1][0];
                tmp_down._branchratio = (itr->second)[Branchratio]+(*DM0_down_Redu)[4*icounter+2][0];
                tmp_down._phase = 0;
                passfix = true;
          }else if(!passfix){
                tmp_up._mass  = (itr->second)[Mass] + (*DM0_up_Redu)[4*icounter][0];
                tmp_up._width = (itr->second)[Width] + (*DM0_up_Redu)[4*icounter+1][0];
                tmp_up._phase = (itr->second)[Phase] + (*DM0_up_Redu)[4*icounter+2][0];
                tmp_down._mass  = (itr->second)[Mass] + (*DM0_down_Redu)[4*icounter][0];
                tmp_down._width = (itr->second)[Width] + (*DM0_down_Redu)[4*icounter+1][0];
                tmp_down._phase = (itr->second)[Phase] + (*DM0_down_Redu)[4*icounter+2][0];
                tmp_up._branchratio = (itr->second)[Branchratio]+(*DM0_up_Redu)[4*icounter+3][0];
                tmp_down._branchratio = (itr->second)[Branchratio]+(*DM0_down_Redu)[4*icounter+3][0];
          }else{
                tmp_up._mass  = (itr->second)[Mass] + (*DM0_up_Redu)[4*icounter-1][0];
                tmp_up._width = (itr->second)[Width] + (*DM0_up_Redu)[4*icounter][0];
                tmp_up._phase = (itr->second)[Phase] + (*DM0_up_Redu)[4*icounter+1][0];
                tmp_down._mass  = (itr->second)[Mass] + (*DM0_down_Redu)[4*icounter-1][0];
                tmp_down._width = (itr->second)[Width] + (*DM0_down_Redu)[4*icounter][0];
                tmp_down._phase = (itr->second)[Phase] + (*DM0_down_Redu)[4*icounter+1][0];
                tmp_up._branchratio = (itr->second)[Branchratio]+(*DM0_up_Redu)[4*icounter+2][0];
                tmp_down._branchratio = (itr->second)[Branchratio]+(*DM0_down_Redu)[4*icounter+2][0];
          }
          chi2_hess_up.AddRes(itr->first, tmp_up);
          chi2_hess_down.AddRes(itr->first, tmp_down);
/*
              std::cout<<"DEBUG::In function MinosInitialization Resonance["<<itr->first<<" Up : Mass = "<<tmp_up[Mass]<<std::endl;
              std::cout<<"DEBUG::In function MinosInitialization Resonance["<<itr->first<<" Up : Width = "<<tmp_up[Width]<<std::endl;
              std::cout<<"DEBUG::In function MinosInitialization Resonance["<<itr->first<<" Up : Phase = "<<tmp_up[Phase]<<std::endl;
              std::cout<<"DEBUG::In function MinosInitialization Resonance["<<itr->first<<" Up : Branchratio = "<<tmp_up[Branchratio]<<std::endl;
*/
              icounter += 1;
   }
   chi2_hess_up.FixPhaseAngle(_para_central.GetFixedPhaseRes());
   chi2_hess_down.FixPhaseAngle(_para_central.GetFixedPhaseRes());
   std::string name_up(resname);
   std::string name_down(resname);
   name_up.append("_up"); name_down.append("_down");

   _para_limits.insert(std::pair<std::pair<std::string,ResFit::respara>, ResFit::Chi2_Hess>(std::pair<std::string,ResFit::respara>(name_up,variable), chi2_hess_up));
   _para_limits.insert(std::pair<std::pair<std::string,ResFit::respara>, ResFit::Chi2_Hess>(std::pair<std::string,ResFit::respara>(name_down,variable), chi2_hess_down));
 /*  
   for(int i=0;i<4*NRES-1;i++){
        std::cout<<"Initialize Minos : Error["<<i<<"] = "<<sqrt((*hess_reduced)[i][i])<<std::endl;
   }
*/
//        hess_reduced->Invert();
}

bool Chi2_Fitter::UpdateLimitOnce(std::pair<std::string,ResFit::respara> par_idx, double Chi2_Min,bool updatehessstd, double Delta){
    std::map<std::pair<std::string,ResFit::respara>,ResFit::Chi2_Hess>::iterator itr_pl = _para_limits.find(par_idx);
    if(itr_pl== _para_limits.end()){
        std::cout<<"WARNING : limit "<<par_idx.first<<" , value "<<par_idx.second<<" not found, no limits will be updated!"<<std::endl;
        return false;
    } 

    TMatrixD * hess_central = _para_central.GetCovariance_Reduced();
//        hess_central->Invert(); 

    ResFit::Chi2_Hess par_limit = itr_pl->second;
    std::map<std::string,ResFit::resonance> v_resonance_tmp = par_limit.GetResonance();
    (par_limit).ReSet();
/*
    for(std::map<std::string,ResFit::resonance>::iterator itr=v_resonance_tmp.begin(); itr!=v_resonance_tmp.end();++itr){
         std::cout<<"Before Chi2 Calculation , Resonance["<<itr->first<<"] Mass : "<<(itr->second)[Mass]<<std::endl;
         std::cout<<"Before Chi2 Calculation , Resonance["<<itr->first<<"] Width : "<<(itr->second)[Width]<<std::endl;
         std::cout<<"Before Chi2 Calculation , Resonance["<<itr->first<<"] Phase : "<<(itr->second)[Phase]<<std::endl;
         std::cout<<"Before Chi2 Calculation , Resonance["<<itr->first<<"] Branchratio : "<<(itr->second)[Branchratio]<<std::endl;
    }
    std::cout<<"Chi2 before chi2 calculation : "<<par_limit.GetChi2()<<std::endl;
*/
    for_each(_data.begin(),_data.end(),std::ref(par_limit));
    int NRES =_para_central.GetNRes() ;

    double eta = Chi2_Min+Delta-(par_limit).GetChi2();
    TMatrixD * MK0 = new TMatrixD(4*NRES,1);
    TMatrixD * firstderv = (par_limit).GetChi2Gradient();
    
    std::map<std::string,ResFit::resonance> v_resonance = par_limit.GetResonance();
    int fixed_phase_index = par_limit.GetVarIndex(par_limit.GetFixedPhaseRes(),ResFit::Phase);
    if(par_idx.first.length()!=5+par_idx.first.rfind("_down") && par_idx.first.length()!=3+par_idx.first.rfind("_up")){
         std::cout<<"ERROR : Illegal limit name "<<par_idx.first<<", must be ended with \'_up\' or \'_down\'!"<<std::endl;
         return false;
    }
    std::string var_name("");
    if(par_idx.first.length()==5+par_idx.first.rfind("_down") ){
         var_name = (par_idx.first).substr(par_idx.first.find_first_not_of("_down"),par_idx.first.length()-5);
    }else{
         var_name = (par_idx.first).substr(par_idx.first.find_first_not_of("_up"),par_idx.first.length()-3);
    }
    int ipar =  _para_central.GetVarIndex_Reduced(var_name, par_idx.second);

//    std::cout<<"DEBUG: IN function UpdateLimitOnce, reduced parameter idx = "<<ipar<<std::endl;
    /*    
        for(int i=0;i<4*NRES;i++){
          (*firstderv)[i][0] = v_firstderv.at(i);
        }
*/
    for(int i=0;i<4*NRES-1;i++){
       (*MK0)[i+1][0] = -((i<fixed_phase_index)?(*firstderv)[i][0]:(*firstderv)[i+1][0]);//To Be Checked!          
    }
    (*MK0)[0][0] = eta;
//    std::cout<<"DEBUG: IN function UpdateLimitOnce, eta = "<<eta<<std::endl;
    TMatrixD * secondderv = par_limit.GetHess();
    TMatrixD * MH0 = new TMatrixD(4*NRES,4*NRES);

    for(int i=0;i<4*NRES;i++){
       if(i==0){
           for(int j=0;j<4*NRES;j++){
               if(j==fixed_phase_index)continue;
               (*MH0)[i][((j<fixed_phase_index)?j:j-1)] = (*firstderv)[j][0];
           }
           (*MH0)[i][4*NRES-1] = 0;
           continue;
       }
       for(int j=0;j<4*NRES-1;j++){
           int jj = (j<fixed_phase_index)?j:j+1;
           int ii = (i-1<fixed_phase_index?i-1:i);
           (*MH0)[i][j] = (*secondderv)[ii][jj];
       }
       if(i-1 != ipar)(*MH0)[i][4*NRES-1] = 0;
       else (*MH0)[i][4*NRES-1] = -1;
    }

    MH0->Invert();
    TMatrixD * Epsilon = new TMatrixD(4*NRES,1);
    (*Epsilon) = (*MH0)*(*MK0);
    TMatrixD *Epsilon_Aux = new TMatrixD(4*NRES,1);

    for(int i=0;i<4*NRES;i++){
          if(i<fixed_phase_index)
             (*Epsilon_Aux)[i][0] = (*Epsilon)[i][0];
          else if(i==fixed_phase_index)
             (*Epsilon_Aux)[i][0]  = 0;
          else (*Epsilon_Aux)[i][0] = (*Epsilon)[i-1][0];
    }
/*
        for(std::map<std::string,ResFit::resonance>::iterator itr=v_resonance.begin(); itr!=v_resonance.end();++itr){
             std::cout<<"Before Update , Resonance["<<itr->first<<"] Mass : "<<(itr->second)[Mass]<<std::endl;
             std::cout<<"Before Update , Resonance["<<itr->first<<"] Width : "<<(itr->second)[Width]<<std::endl;
             std::cout<<"Before Update , Resonance["<<itr->first<<"] Phase : "<<(itr->second)[Phase]<<std::endl;
             std::cout<<"Before Update , Resonance["<<itr->first<<"] Branchratio : "<<(itr->second)[Branchratio]<<std::endl;
        }
*/
/*
    for(int i=0;i<4*NRES;i++){
         std::cout<<"Epsilon["<<i<<"] = "<<(*Epsilon_Aux)[i][0]<<std::endl;
    }
*/
    int icounter =0;
    for(std::map<std::string,ResFit::resonance>::iterator itr=v_resonance.begin(); itr!=v_resonance.end();++itr){
         resonance tmp_res;
         tmp_res._mass = (itr->second)[Mass]+(*Epsilon_Aux)[icounter*4][0];
         tmp_res._width = (itr->second)[Width]+(*Epsilon_Aux)[icounter*4+1][0];
         tmp_res._phase = (itr->second)[Phase]+(*Epsilon_Aux)[icounter*4+2][0];
         tmp_res._branchratio = (itr->second)[Branchratio]+(*Epsilon_Aux)[icounter*4+3][0];
         //par_limit.UpdateRes(itr->first, tmp_res); icounter+=1; //Why doesn't it work?
         (itr_pl->second).UpdateRes(itr->first, tmp_res); icounter +=1;
    }
 /*  
        for(std::map<std::string,ResFit::resonance>::iterator itr=v_resonance.begin(); itr!=v_resonance.end();++itr){
             std::cout<<"After Update , Resonance["<<itr->first<<"] Mass : "<<(itr->second)[Mass]<<std::endl;
             std::cout<<"After Update , Resonance["<<itr->first<<"] Width : "<<(itr->second)[Width]<<std::endl;
             std::cout<<"After Update , Resonance["<<itr->first<<"] Phase : "<<(itr->second)[Phase]<<std::endl;
             std::cout<<"After Update , Resonance["<<itr->first<<"] Branchratio : "<<(itr->second)[Branchratio]<<std::endl;
        }
*/
    bool result = true;
    for(int i =0; i<4*NRES-1;i++){
//         std::cout<<"Error["<<i<<"] = "<<sqrt((*hess_central)[i][i])<<std::endl;
         if(fabs((*Epsilon)[i][0])>sqrt((*hess_central)[i][i])*_converge_criteria){
             result = false; break;
         }
    }
//        hess_central->Invert();
//        delete firstderv;
    delete MH0;
    delete MK0;
    delete Epsilon;
    delete Epsilon_Aux;
    return result;

}

void Chi2_Fitter::SetLimit(std::pair<std::string,ResFit::respara>  par_limit_idx, bool updatehessstd, double Delta){
    if(updatehessstd)
       this->UpdateHess();
    this->InitialMinos(par_limit_idx.first,par_limit_idx.second,false,Delta);
    bool finished=false;
    double chi2_min = _para_central.GetChi2();
    std::string name_up(par_limit_idx.first); name_up.append("_up");
    std::string name_down(par_limit_idx.first);  name_down.append("_down");

    while(!finished){
        finished = this->UpdateLimitOnce(std::pair<std::string,ResFit::respara>(name_up,par_limit_idx.second),chi2_min,false,Delta);
    }


    finished = false;
    while(!finished){
        finished = this->UpdateLimitOnce(std::pair<std::string,ResFit::respara>(name_down,par_limit_idx.second),chi2_min,false,Delta);
    }

}

void Chi2_Fitter::SetLimit(Chi2_Hess & chi2_hess, std::string resname_x, std::string resname_y, respara para_x, respara para_y, std::string res_with_fixed_angle, double chi2_min, double Delta){   
//    std::cout<<"Mark N1"<<std::endl;
    int idx_fixed_angle = chi2_hess.GetFixedIndex();
//    std::cout<<"Mark N2"<<std::endl;
    int idx = chi2_hess.GetVarIndex(resname_x,para_x);
//    std::cout<<"Mark N3"<<std::endl;
    int idx_reduced = chi2_hess.GetVarIndex_Reduced(resname_x,para_x);
//    std::cout<<"Mark N4"<<std::endl;
    int idy = chi2_hess.GetVarIndex(resname_y,para_y);
//    std::cout<<"Mark N5"<<std::endl;
    int idy_reduced = chi2_hess.GetVarIndex_Reduced(resname_y,para_y);
//    std::cout<<"Mark N6"<<std::endl;
    int NRES = chi2_hess.GetNRes() ;
//    std::cout<<"Mark N7, NRES = "<<NRES<<", idx = "<<idx_reduced<<", idy = "<<idy_reduced<<std::endl;
    TMatrixD * H_expanded = new TMatrixD(4*NRES+1,4*NRES+1);
    TMatrixD* v_vector = new TMatrixD(4*NRES+1,1);
    TMatrixD* u_vector = new TMatrixD(4*NRES+1,1);//Hu=v, resonance+=u
//    std::cout<<"Mark N8"<<std::endl;
    bool finished = false;
    int loop_counter = 0;
    while(!finished){
        chi2_hess.ReSet();
//        std::cout<<"Mark A1"<<std::endl;
        for_each(_data.begin(), _data.end(), std::ref(chi2_hess));
//        std::cout<<"Mark A2"<<std::endl;
        double eta = chi2_min+Delta-chi2_hess.GetChi2();
//         std::cout<<"Mark A3"<<std::endl;
        TMatrixD * gradient_reduced = chi2_hess.GetChi2Gradient_Reduced();
//         std::cout<<"Mark A4"<<std::endl;
//Set v vector
//        std::cout<<"Mark A"<<std::endl;
        (*v_vector)[0][0] = eta;
        for(int i =1;i<4*NRES;i++){
            (*v_vector)[i][0] = -(*gradient_reduced)[i-1][0];
        }
        (*v_vector)[4*NRES][0] = 0 ;
//        std::cout<<"Mark B"<<std::endl;
//Set H_expanded
        TMatrixD * hess_reduced = chi2_hess.GetHess_Reduced();
        for(int j=0;j<4*NRES+1;j++){
            if(j<4*NRES-1)(*H_expanded)[0][j] = (*gradient_reduced)[j][0];
            else (*H_expanded)[0][j] = 0;
        }
//        std::cout<<"Mark C"<<std::endl;   
        for(int i=1;i<4*NRES;i++){
           for(int j=0;j<4*NRES-1;j++){
               (*H_expanded)[i][j] = (*hess_reduced)[i-1][j];
           }
           (*H_expanded)[i][4*NRES-1] = (i-1==idx_reduced)?-1:0;
           (*H_expanded)[i][4*NRES] = (i-1==idy_reduced)?-1:0;
        }
//        std::cout<<"Mark D"<<std::endl;
        for(int j=0;j<4*NRES+1;j++){
           (*H_expanded)[4*NRES][j] = (j==idx_reduced)?1:0;
        }
//        std::cout<<"Mark E"<<std::endl;
//Now get u vector
/*
        std::cout<<"Loop "<<loop_counter<<std::endl;
        std::cout<<"eta = "<<eta<<std::endl;
        std::cout<<"Printing Gradient"<<std::endl;
        for(int i =0;i<4*NRES-1;i++){
              std::cout<<(*gradient_reduced)[i][0]<<"\t";
        }std::cout<<std::endl;
        std::cout<<"Printing expanded H matrix"<<std::endl;
        for(int i=0;i<4*NRES+1;i++){
            for(int j=0; j<4*NRES+1; j++){
                std::cout<<(*H_expanded)[i][j]<<"\t";
            }std::cout<<std::endl;
        }
*/
        H_expanded->Invert();
//        std::cout<<"Mark F"<<std::endl;
        (*u_vector) = (*H_expanded)*(*v_vector);
//        std::cout<<"Mark G"<<std::endl;
//Upate Resonance
        int counter = 0;
        std::map<std::string,ResFit::resonance> v_resonance = chi2_hess.GetResonance();
        for(std::map<std::string,ResFit::resonance>::iterator itr = v_resonance.begin(); itr!= v_resonance.end(); ++itr){
             ResFit::resonance tmp_res;
             if(itr->first == _para_central.GetFixedPhaseRes()){
                 tmp_res._mass = (itr->second)[Mass]+(*u_vector)[counter+Mass][0];
                 tmp_res._width = (itr->second)[Width]+(*u_vector)[counter+Width][0];
                 tmp_res._phase = 0;
                 tmp_res._branchratio = (itr->second)[Branchratio]+(*u_vector)[counter+Branchratio-1][0];
             }else{
                 tmp_res._mass = (itr->second)[Mass]+(*u_vector)[counter+Mass][0];
                 tmp_res._width = (itr->second)[Width]+(*u_vector)[counter+Width][0];
                 tmp_res._phase =  (itr->second)[Phase]+(*u_vector)[counter+Phase][0];
                 tmp_res._branchratio = (itr->second)[Branchratio]+(*u_vector)[counter+Branchratio][0];
             }
             chi2_hess.UpdateRes(itr->first,tmp_res);            
             if(itr->first == _para_central.GetFixedPhaseRes())counter+=3;
             else counter+=4;
        }
        loop_counter += 1;
//Determine the to continue the loop or not
/*
        std::cout<<"Printing V matrix"<<std::endl;
        for(int i =0;i<4*NRES+1;i++){
              std::cout<<(*v_vector)[i][0]<<"\t";
        }std::cout<<std::endl;
        std::cout<<"Printing U matrix"<<std::endl;
        for(int i =0;i<4*NRES+1;i++){
              std::cout<<(*u_vector)[i][0]<<"\t";
        }std::cout<<std::endl;
*/
        bool finished_in_this_step = true;
        TMatrixD * covariance_central = _para_central.GetCovariance_Reduced();
        for(int i=0;i<4*NRES-1;i++){
             if(fabs((*u_vector)[i][0])>_converge_criteria*sqrt((*covariance_central)[i][i])){
                  finished_in_this_step = false; continue;
             }
        }
        finished = finished_in_this_step;
    }


    delete H_expanded;
    delete v_vector;
    delete u_vector;
}

void Chi2_Fitter::PrintLimits(){
    if(_para_limits.empty())return;
    for(std::map<std::pair<std::string, ResFit::respara>,Chi2_Hess>::const_iterator itr = _para_limits.begin();itr!=_para_limits.end();++itr){
//             int idx = (itr->second).GetVarIndex();
         std::cout<<(itr->first).first<<"\t"<<(itr->first).second<<"\t"<<std::endl;
    }
}

void Chi2_Fitter::FillContour(Contour & c, std::string name_x, std::string name_y, respara idx, respara idy, int NP, bool update,double delta){
     int fixed_phase_index = _para_central.GetVarIndex(_para_central.GetFixedPhaseRes(),ResFit::Phase);
     int NRES =_para_central.GetNRes() ;

     double chi2_min = this->GetChi2Min();
     int var_x_idx_reduced =  _para_central.GetVarIndex_Reduced(name_x, idx);
     int var_x_idx =  _para_central.GetVarIndex(name_x, idx);
     int var_y_idx_reduced =  _para_central.GetVarIndex_Reduced(name_y, idy);
     int var_y_idx =  _para_central.GetVarIndex(name_y, idy);
//Clear Data in contour
     c.ClearData();

//Check List NP Number,at least 8+1=9 points
     if(NP%4!=0)NP+= (4-NP%4);
     if(NP<=4)NP+=4;

     this->ReSet();
     this->UpdateHess();

     double x_central = _para_central.GetParameter(name_x, idx);
     double y_central = _para_central.GetParameter(name_y, idy); 
//     std::cout<<"var_x_idx_reduced = "<<var_x_idx_reduced<<",\t var_y_idx_reduced = "<<var_y_idx_reduced<<std::endl;
     TMatrixD * hess_reduced = _para_central.GetCovariance_Reduced();
     TMatrixD *delta0 =new TMatrixD(2,2);
     (*delta0)[0][0] = ((*hess_reduced)[var_x_idx_reduced][var_x_idx_reduced]);
     (*delta0)[1][1] = ((*hess_reduced)[var_y_idx_reduced][var_y_idx_reduced]);
     (*delta0)[0][1] = ((*hess_reduced)[var_y_idx_reduced][var_x_idx_reduced]);
     (*delta0)[1][0] = ((*hess_reduced)[var_x_idx_reduced][var_y_idx_reduced]);
     delta0->Invert();
     int nq = NP/4.0;
     double x_min_el = -sqrt(delta*((*hess_reduced)[var_x_idx_reduced][var_x_idx_reduced]));  //With ellpitcal hypothesis
     double x_max_el = sqrt(delta*((*hess_reduced)[var_x_idx_reduced][var_x_idx_reduced]));
     double y_x_min_el =0;   //y value at x_min with ellplitcal hypothesis 
     double y_x_max_el = 0;   //y value at x_max with ellplitical hypothesis
//Filling elliptical curve
     if(!update){    
          double x_min = x_min_el;
          double x_max = x_max_el;
          double step = (x_max-x_min)/double(2*nq);
          double y[2];
          Solve_y_by_x(y,x_min,delta,delta0);
          c.Insert(0,x_min+x_central,y[0]+y_central);
          y_x_min_el = 0;
          for(int i=1;i<=2*nq;i++){
             double x = x_min+double(i)*step;
             Solve_y_by_x(y,x,delta,delta0);
             c.Insert(-i,x+x_central,y[0]+y_central);
             c.Insert(i,x+x_central,y[1]+y_central);
             if(i==2*nq){
                  y_x_max_el = y[0];
             }
          }
          c.SetNFP(4*nq+1);
          delete delta0;
          return;
     }

     double y_tmp[2];
     Solve_y_by_x(y_tmp,x_min_el,delta,delta0);
     y_x_min_el = y_tmp[0];
//     std::cout<<"at x_min, y_tmp[0] = "<<y_tmp[0]<<", y_tmp[1] = "<<y_tmp[1]<<std::endl;
     Solve_y_by_x(y_tmp,x_max_el,delta,delta0);
     y_x_max_el = y_tmp[0];
//     std::cout<<"at x_max, y_tmp[0] = "<<y_tmp[0]<<", y_tmp[1] = "<<y_tmp[1]<<std::endl;

     TMatrixD * K_x_max_el21 = new TMatrixD(2,1);
     TMatrixD * K_x_min_el21 = new TMatrixD(2,1);
     TMatrixD * K_x_max_el = new TMatrixD(4*NRES-1,1);//Gradient 
     TMatrixD * K_x_min_el = new TMatrixD(4*NRES-1,1);
     TMatrixD * u_x_max_el = new TMatrixD(4*NRES-1,1); //INitial value with eplli hypothesis, values are centered at zero
     TMatrixD * u_x_min_el = new TMatrixD(4*NRES-1,1);
     
     TMatrixD * v_x_max_el21 = new TMatrixD(2,1);
     TMatrixD * v_x_min_el21 = new TMatrixD(2,1);

     TMatrixD * u_x_max = new TMatrixD(4*NRES-1,1); //Inital value by Minos, centered at central value
     TMatrixD * u_x_min = new TMatrixD(4*NRES-1,1);

     (*v_x_max_el21)[0][0] = x_max_el;  (*v_x_max_el21)[1][0] = y_x_max_el;
     (*v_x_min_el21)[0][0] = x_min_el;  (*v_x_min_el21)[1][0] = y_x_min_el;

     (*K_x_max_el21) = 2.0*(*delta0)*(*v_x_max_el21);
     (*K_x_min_el21) = 2.0*(*delta0)*(*v_x_min_el21);

     for(int i=0;i<4*NRES-1;i++){
          if(i==var_x_idx_reduced){
             (*K_x_max_el)[i][0] = (*K_x_max_el21)[0][0];
             (*K_x_min_el)[i][0] = (*K_x_min_el21)[0][0];
          }else if(i == var_y_idx_reduced){
              (*K_x_max_el)[i][0] = (*K_x_max_el21)[1][0];
              (*K_x_min_el)[i][0] = (*K_x_min_el21)[1][0];
          }else{
              (*K_x_max_el)[i][0] = 0;
              (*K_x_min_el)[i][0] = 0;
          }
     }

     (*u_x_max_el)=0.5*(*hess_reduced)*(*K_x_max_el);
     (*u_x_min_el)=0.5*(*hess_reduced)*(*K_x_min_el);

//Iteration for points filling, now get min and max value of x         
     this->SetLimit(std::pair<std::string,ResFit::respara>(name_x,idx),false,delta);
     std::string name_x_up(name_x); std::string name_x_down(name_x); name_x_up.append("_up"); name_x_down.append("_down");
     std::pair<std::string, respara> par_x_up = std::pair<std::string,respara>(name_x_up, idx);
     std::pair<std::string, respara> par_x_down = std::pair<std::string,respara>(name_x_down, idx);
/*
     std::cout<<"par_x_up: key = "<<par_x_up.first<<"\t, par_x_up: vlaue = "<<idx<<std::endl;
     std::cout<<"par_y_down: key = "<<par_x_down.first<<"\t, par_y_down: vlaue = "<<idy<<std::endl;
*/
     double x_max = _para_limits[par_x_up].GetParameter(name_x,idx);
     double x_min = _para_limits[par_x_down].GetParameter(name_x,idx);
     double y_x_max = _para_limits[par_x_up].GetParameter(name_y,idy);
     double y_x_min = _para_limits[par_x_down].GetParameter(name_y,idy);

     std::map<std::string,resonance> res_limit_up = _para_limits[par_x_up].GetResonance();
     std::map<std::string,resonance> res_limit_down = _para_limits[par_x_down].GetResonance();

     int counter = 0;
     for(std::map<std::string,resonance>::iterator itr = res_limit_up.begin(); itr!=res_limit_up.end();++itr){
          if(itr->first==_para_central.GetFixedPhaseRes()){
               (*u_x_max)[counter+Mass][0] = (itr->second)[Mass];
               (*u_x_max)[counter+Width][0] = (itr->second)[Width];
               (*u_x_max)[counter+Branchratio-1][0] = (itr->second)[Branchratio];
               counter += 3;               
          }else{
               (*u_x_max)[counter+Mass][0] = (itr->second)[Mass];
               (*u_x_max)[counter+Width][0] = (itr->second)[Width];
               (*u_x_max)[counter+Phase][0] = (itr->second)[Phase];
               (*u_x_max)[counter+Branchratio][0] = (itr->second)[Branchratio];
               counter += 4;
          }
     }

     counter = 0;
     for(std::map<std::string,resonance>::iterator itr = res_limit_down.begin(); itr!=res_limit_down.end();++itr){
          if(itr->first==_para_central.GetFixedPhaseRes()){
               (*u_x_min)[counter+Mass][0] = (itr->second)[Mass];
               (*u_x_min)[counter+Width][0] = (itr->second)[Width];
               (*u_x_min)[counter+Branchratio-1][0] = (itr->second)[Branchratio];
               counter += 3;
          }else{
               (*u_x_min)[counter+Mass][0] = (itr->second)[Mass];
               (*u_x_min)[counter+Width][0] = (itr->second)[Width];
               (*u_x_min)[counter+Phase][0] = (itr->second)[Phase];
               (*u_x_min)[counter+Branchratio][0] = (itr->second)[Branchratio];
               counter += 4;
          }
     }

     TMatrixD *KK = new TMatrixD(4*NRES-1,4*NRES-1); //Slope, minos_v =kk*ellpitcal_v+bb
     TMatrixD *BB = new TMatrixD(4*NRES-1,1);//Intercept

     for(int i=0;i < 4*NRES-1;i++){
         (*BB)[i][0] = 0.5*((*u_x_max)[i][0]+(*u_x_min)[i][0]);
         for(int j=0;j<4*NRES-1;j++){
              if(i!=j)(*KK)[i][j] = 0;
              else {
                     (*KK)[i][j] = ((*u_x_max)[i][0]-(*u_x_min)[i][0])/((*u_x_max_el)[i][0]-(*u_x_min_el)[i][0]);
/*
                     std::cout<<"u_x_max["<<i<<"] = "<<(*u_x_max)[i][0]<<",\t u_x_min["<<i<<"] = "<<(*u_x_min)[i][0]<<std::endl;
                     std::cout<<"u_x_max_el["<<i<<"] = "<<(*u_x_max_el)[i][0]<<",\t u_x_min_el["<<i<<"] = "<<(*u_x_min_el)[i][0]<<std::endl;
*/
                   }
         }
     }
 
     double step = (x_max_el-x_min_el)/double(2*nq);
     c.Insert(0,x_min,y_x_min);
     TMatrixD *tmp_el_up = new TMatrixD(4*NRES-1,1);
     TMatrixD *tmp_el_down = new TMatrixD(4*NRES-1,1);

     TMatrixD *tmp_up = new TMatrixD(4*NRES-1,1);
     TMatrixD *tmp_down = new TMatrixD(4*NRES-1,1);
     
     std::map<std::string,resonance> res_cent = _para_central.GetResonance();

     for(int i=1;i<2*nq;i++){
//Initialization
         Chi2_Hess tmp_chi2_up, tmp_chi2_down;
         double y[2];//y[0]down, y[1]up
         double x = x_min_el+step*double(i);
         Solve_y_by_x(y,x,delta,delta0); 
         Solve_u_by_x(tmp_el_up, hess_reduced, x,  y[1],  delta, delta0 , NRES, var_x_idx_reduced, var_y_idx_reduced);
         Solve_u_by_x(tmp_el_down, hess_reduced, x,  y[0],  delta, delta0 , NRES, var_x_idx_reduced, var_y_idx_reduced);
         (*tmp_up) = (*KK)*(*tmp_el_up)+(*BB);
         (*tmp_down) = (*KK)*(*tmp_el_down)+(*BB);
         counter = 0;

         for(std::map<std::string, resonance>::iterator itr = res_cent.begin();itr!=res_cent.end();++itr){                                
              if((itr->first)==_para_central.GetFixedPhaseRes()){
                 tmp_chi2_up.AddRes(itr->first,(*tmp_up)[counter+Mass][0],(*tmp_up)[counter+Width][0],0,(*tmp_up)[counter+Branchratio-1][0]);
                 counter+=3;
              }else{
                 tmp_chi2_up.AddRes(itr->first,(*tmp_up)[counter+Mass][0],(*tmp_up)[counter+Width][0],(*tmp_up)[counter+Phase][0],(*tmp_up)[counter+Branchratio][0]);
                 counter+=4;
              }
         }

         counter = 0;
         for(std::map<std::string, resonance>::iterator itr = res_cent.begin();itr!=res_cent.end();++itr){
              if((itr->first)==_para_central.GetFixedPhaseRes()){
                 tmp_chi2_down.AddRes(itr->first,(*tmp_down)[counter+Mass][0],(*tmp_down)[counter+Width][0],0,(*tmp_down)[counter+Branchratio-1][0]);
                 counter+=3;
              }else{
                 tmp_chi2_down.AddRes(itr->first,(*tmp_down)[counter+Mass][0],(*tmp_down)[counter+Width][0],(*tmp_down)[counter+Phase][0],(*tmp_down)[counter+Branchratio][0]);
                 counter+=4;
              }
         }

//Update 
/*
         std::cout<<"======================================Point indx : "<<i<<"==============================================="<<std::endl;
         std::cout<<"x = "<<x<<", y[0] = "<<y[0]<<", y[1] = "<<y[1]<<std::endl;
         std::cout<<"Printing KK"<<std::endl;
         for(int i= 0; i<4*NRES-1;i++){
              std::cout<<(*KK)[i][0]<<"\t";
         }std::cout<<std::endl;
         std::cout<<"Printing BB"<<std::endl;
         for(int i= 0; i<4*NRES-1;i++){
              std::cout<<(*BB)[i][0]<<"\t";
         }std::cout<<std::endl;
         std::cout<<"Printing tmp up el resonance"<<std::endl;
         for(int i= 0; i<4*NRES-1;i++){
              std::cout<<(*tmp_el_up)[i][0]<<"\t";
         }std::cout<<std::endl;
         std::cout<<"Printing tmp down el resonance"<<std::endl;
         for(int i= 0; i<4*NRES-1;i++){
              std::cout<<(*tmp_el_down)[i][0]<<"\t";
         }std::cout<<std::endl;
         std::cout<<"Printing tmp up resonance"<<std::endl;
         for(int i= 0; i<4*NRES-1;i++){
              std::cout<<(*tmp_up)[i][0]<<"\t";
         }std::cout<<std::endl;
         std::cout<<"Printing tmp down resonance"<<std::endl;
         for(int i= 0; i<4*NRES-1;i++){
              std::cout<<(*tmp_down)[i][0]<<"\t";
         }std::cout<<std::endl;
*/
         tmp_chi2_up.FixPhaseAngle(_para_central.GetFixedPhaseRes());  tmp_chi2_down.FixPhaseAngle(_para_central.GetFixedPhaseRes());
         SetLimit(tmp_chi2_up , name_x, name_y, idx, idy, _para_central.GetFixedPhaseRes(),  chi2_min,  delta);
         SetLimit(tmp_chi2_down , name_x, name_y, idx, idy, _para_central.GetFixedPhaseRes(),  chi2_min,  delta);
         c.Insert(i,tmp_chi2_up.GetParameter(name_x,idx),tmp_chi2_up.GetParameter(name_y,idy));
         c.Insert(-i,tmp_chi2_down.GetParameter(name_x,idx),tmp_chi2_down.GetParameter(name_y,idy));
     }
     c.Insert(2*nq,x_max,y_x_max);
     c.Insert(-2*nq,x_max,y_x_max);
     delete K_x_max_el21;
     delete K_x_min_el21;
     delete K_x_max_el;
     delete K_x_min_el;
     delete u_x_max_el;
     delete u_x_min_el;
     delete v_x_max_el21;
     delete v_x_min_el21;
     delete u_x_max;
     delete u_x_min;
     delete KK; delete BB;
     delete tmp_el_up;
     delete tmp_el_down;
     delete tmp_up;
     delete tmp_down;
}

void Solve_y_by_x(double y[2], double x, double delta, TMatrixD* inverse_cov){
     double A =  (*inverse_cov)[0][0] ; double B= 0.5*((*inverse_cov)[0][1] + (*inverse_cov)[1][0]); double C= (*inverse_cov)[1][1] ;
     double D = B*B*x*x-C*(A*x*x-delta);
//     std::cout<<"In function Solve_y_by_x : x = "<<x<<", delta = "<<delta<<", Delta = "<<B*B*x*x-C*(A*x*x-delta)<<std::endl;
     if(!((A*C-B*B)>0)){
          std::cout<<"ERROR::2*2 Matrix not positive determined"<<std::endl;
          y[0] = 0;
          y[1] = 0; return;
     }
     if(D<0){
          y[0] = -B*x/(C);
          y[1] = -B*x/(C);
          return; 
     }else{
         y[0] = (-B*x-sqrt(D))/(C);
         y[1] = (-B*x+sqrt(D))/(C);
         return;
     }
}

void Solve_u_by_x(TMatrixD* u, TMatrixD * central_covariance ,double x, double y, double delta, TMatrixD * inverse_cov, int NRES, int idx_reduced, int idy_reduced){

     TMatrixD * K_x_el21 = new TMatrixD(2,1);
     TMatrixD * K_x_el = new TMatrixD(4*NRES-1,1);
     TMatrixD * v_x_el21 = new TMatrixD(2,1);

     (*v_x_el21)[0][0] = x; (*v_x_el21)[1][0] = y;
     (*K_x_el21)=2.0*(*inverse_cov)*(*v_x_el21);

     for(int i=0;i<4*NRES-1;i++){
          if(i==idx_reduced){
             (*K_x_el)[i][0] = (*K_x_el21)[0][0];
          }else if(i == idy_reduced){
             (*K_x_el)[i][0] = (*K_x_el21)[1][0];
          }else{
             (*K_x_el)[i][0] = 0;
          }
     }

     (*u) = 0.5*(*central_covariance)*(*K_x_el);
     delete K_x_el21;
     delete K_x_el;
     delete v_x_el21;
     return;
}

}
