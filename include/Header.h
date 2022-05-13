#ifndef HEAD_H
#define HEAD_H
#include <cmath>
#include <complex>
#include <vector>
#include <map>
#include <algorithm>
//#include <functional>
//#include "TGraph.h"
//#include "TComplex.h"
#include "TMath.h"
#include "TMatrixD.h"
#include "Contour.h"
#include "AnalyticPhi23.h"
#define UNIT 0.38936951

namespace ResFit{
 
enum respara {Mass, Width, Phase, Branchratio};

struct resonance{    //Resonance class
    double _mass;   //if mass, width, branchratio and phase are set, the resonance is parametreized_1 = true
    double _width;
    double _branchratio;
    double _phase;
    std::complex<double> _pole;  //if mass, width, branchratio and phase are set, the resonance is parametreized_2 = true, either parametreized_1 or parametreized_2, or both are true
    std::complex<double> _coeff;
    bool _parameterized1;
    bool _parameterized2;
    double operator[] (respara ipar){
        double return_value = 0;
        switch(ipar){  
           case Mass : {return_value = _mass;break;}
           case Width : {return_value = _width;break;}
           case Branchratio : {return_value = _branchratio; break;}
           case Phase : {return_value = _phase; break;}
           default: break;
        }
        return return_value;
    }
};

struct Bdata {       //Binned Data
    double penergy;
    double pdata;
    double pdataerror;
    double pphasespace;
};
/*
inline void SetDataPhaseSpace(Bdata & data){
    data.pphasespace = AnalyticPhi23(data.penergy);
}
*/
inline double maxabs(double a, double b){
     return (fabs(a)>fabs(b)?a:b);
}

inline double GetXS(std::complex<double> amp, double energy, double phasespace){
     return  UNIT*std::norm(amp)*12.0*TMath::Pi()*phasespace/(energy*energy); //Interferance background
}

inline double GetXSConst(double energy, double phasespace){
    return UNIT*12.0*TMath::Pi()*phasespace/(energy*energy);
}

class Model{  //Base Class of Chi2, HESS, Minimize, Minos, Contour
     protected:
        int _nres;
        std::string _fixed_phase_resonance;
        int _fixed_phase_paridx; //index of parameter of fixed phase angle
        std::map<std::string,ResFit::resonance> _vres;
     public:
        Model():_nres(0), _fixed_phase_resonance(""),_fixed_phase_paridx(-1), _vres(std::map<std::string,ResFit::resonance>()){}
        void AddRes(std::string resname, ResFit::resonance r);
        void AddRes(std::string resname, double mass, double width, double phase, double branchratio);
        void DelRes(std::string resname);
        void UpdateRes(std::string resname, ResFit::resonance r);
        void FixPhaseAngle(std::string resname); //Fix the PhaseAngle of resonance with name resname to 0
        void FixPhaseAngle(); //Fix the first PhaseAngle of resonance to 0, either FixPhaseAngle and FixPhaseAngle should be called in case of modificantion of _vres
        void CopyRes(ResFit::Model v_source); //Copy the basic information from one resonance list to this Model
/*
        resonance GetResonance(std::string resname)();
        double GetMass(std::string resname);
        double GetWidth(std::string resname);
        double GetFR(std::string resname);
        double GetPhi(std::string resname);
        TComplex GetPole(std::string resname);
        TComplex GetCoeff(std::string resname);
        bool GetParameterizedStatus1(std::string resname);
        bool GetParameterizedStatus2(std::string resname);
        void SetMass(std::string resname, double mass);
        void SetWidth(std::string resname, double width);
        void SetFR(std::string resname, double FR);
        void SetPhi(std::string resname, double Phi);
        void SetPole(std::string resname, TComplex pole);
        void AddRes(std::string, resonance);
*/
        std::map<std::string,ResFit::resonance>::iterator FindRes(std::string resname){return _vres.find(resname);}
        std::map<std::string,ResFit::resonance> GetResonance() const{return _vres;}
        int GetNRes() const{return _nres;}
        std::string GetFixedPhaseRes() const{return _fixed_phase_resonance;}
        int GetFixedIndex() const{return _fixed_phase_paridx;}
        int GetVarIndex(std::string resname, ResFit::respara) const; //Get The row or column index of a resonance variable
        int GetVarIndex_Reduced(std::string resname, ResFit::respara) const;   //Get The row or column index of a resonance variable, excluding fixed phase angle
        double GetParameter(std::string resname, ResFit::respara val){return (_vres[resname])[val];}
};

class Chi2 : public ResFit::Model{
     protected :      
        double _chi2;
     public : 
        Chi2():_chi2(0){}
        double GetChi2() const{return _chi2;}
        void ReSet(){_chi2=0;}     //Necessary when data is updated or resonance parameters changed.
        void operator()(ResFit::Bdata); //Recommand after ReSet
};

class Chi2_Gradient : public ResFit::Chi2{
     protected : 
/*
        std::vector<double> _chi2_gradient;
*/
        TMatrixD * _chi2_gradient;
        TMatrixD * _chi2_gradient_reduced;
        std::map<std::string,double> _phasespace;
        std::map<std::string,double> _phasespace_1stder;
     public :  
        Chi2_Gradient():_chi2_gradient(NULL),_chi2_gradient_reduced(NULL), _phasespace(std::map<std::string,double>()), _phasespace_1stder(std::map<std::string,double>()){}
/*
        std::vector<double> GetChi2Gradient(){return _chi2_gradient;}
*/
        TMatrixD* GetChi2Gradient() const{return _chi2_gradient;}
        TMatrixD* GetChi2Gradient_Reduced() const{return _chi2_gradient_reduced;}
        void ReSet(); //ReSet Chi2 and Gradient and firsDerv of phasespace according to Resonance, necessary when resonance parameters are updated, or data modified
        void operator()(ResFit::Bdata); //Calcuate chi2 and chi2 gradient, recommand after ReSet
};


class Chi2_Hess : public ResFit::Chi2_Gradient{
     protected:
        std::map<std::string,double> _phasespace_2ndder;
        TMatrixD *_chi2_hess;
        TMatrixD *_chi2_hess_reduced;
        TMatrixD *_chi2_covariance;
        TMatrixD *_chi2_covariance_reduced;
     public:
        Chi2_Hess():_chi2_hess(NULL),_chi2_hess_reduced(NULL),_chi2_covariance_reduced(NULL),_chi2_covariance(NULL),_phasespace_2ndder(std::map<std::string,double>()){}
        ~Chi2_Hess(){if(_chi2_hess){delete _chi2_hess;}}
        void  SetCovariance();
        TMatrixD* GetHess() const{return _chi2_hess;}
        TMatrixD* GetHess_Reduced() const{return _chi2_hess_reduced;}
        TMatrixD* GetCovariance() const{return _chi2_covariance;}
        TMatrixD* GetCovariance_Reduced() const{return _chi2_covariance_reduced;}
        void ReSet();
        void operator()(ResFit::Bdata);
};

class Chi2_Fitter{
    private:
       Chi2_Hess _para_central;  //Including central value of parameters, and hessian matrix. It will be updated in minimization
       std::map<std::pair<std::string, ResFit::respara>,Chi2_Hess>  _para_limits;
       int _nfp; // number of points to plot contour
       double _error_definition;   //Error Definition, default value is 1
       double _converge_criteria; //Ratio of update comparing to error
 //      std::map<std::pair<std::pair<std::string, ResFit::respara>,std::pair<std::string, ResFit::respara> >,Contour> _contour;
       std::vector<Bdata> _data;
    public:
       Chi2_Fitter():_nfp(200),_error_definition(1.0),_converge_criteria(1.0E-5),_data(std::vector<Bdata>(0)),_para_central(ResFit::Chi2_Hess()),_para_limits(std::map<std::pair<std::string,ResFit::respara>,ResFit::Chi2_Hess>()){}
       ~Chi2_Fitter(){}
       void AddData(Bdata);
//       void DelData(Delet Criteria);
       void AddRes(std::string resname, resonance res){_para_central.AddRes(resname,res);}
       void AddRes(std::string resname, double mass, double width, double phase, double branchratio);
       void DelRes(std::string resname){_para_central.DelRes(resname);}
       void SetErrorDef(double x){_error_definition = x;}
       void SetPrecision(double x){_converge_criteria = x;}
       void SetNFP(int n){_nfp = n;}
       void ReSet();
/*
       void Minimize();
       void UpdateParameter(std::map<std::string, ResFit::resonance> changed_res);
       void ParameterUpdateCal(ResFit::Chi2 _tmp_chi2, double lamda_mu);
*/
       void UpdateHess();  //Update Hessian matrix of central value
       void InitialMinos(std::string, ResFit::respara, bool updatehess=true, double Delta=1.0);   //Initialize for one parameter
//       void InitialMinos(std::string, std::String, ResFit::respara, ResFit::respara, bool updatehess=true); //
//       void InitialMinos();   //Initialize for all parameters
//       void UpdateMinos();
//       void UpdateContour();
       bool UpdateLimitOnce(std::pair<std::string,ResFit::respara>  par_limit_idx, double Chi2_Min, bool updatehessstd=true, double Delta=1.0); //return false if iteration is not finished, true iteration is ceased.
       void SetLimit(std::pair<std::string,ResFit::respara>  par_limit_idx, bool updatehessstd=true, double Delta=1.0);
       void SetLimit(Chi2_Hess & , std::string resname_x, std::string resname_y, respara para_x, respara para_y, std::string res_with_fixed_angle, double chi2_min, double Delta);//Set limit with fixed x and float y, Chi2_Hess must be initialized
       void FixPhaseAngle();
       void FixPhaseAngle(std::string);
       void PrintLimits();
//       void SetHess();
//       void SetMinos();
       void FillContour(Contour & c, std::string res_name_x, std::string res_name_2, respara value_x, respara value_y, int NP, bool update = false, double delta=1.0);//Fill data in Contour C
       double GetErrorDef() const{return _error_definition;}
       double GetPrecision() const{return _converge_criteria;}
       double GetNFP() const{return _nfp;}
       double GetChi2Min() const{return _para_central.GetChi2();}
       TMatrixD* GetHess() const{return _para_central.GetHess();}
       TMatrixD* GetHess_Reduced() const{return _para_central.GetHess_Reduced();}
       TMatrixD* GetCovariance() const{return _para_central.GetCovariance();}
       TMatrixD* GetCovariance_Reduced() const{return _para_central.GetCovariance_Reduced();}
       double GetParameter(std::string resname, ResFit::respara val) {return _para_central.GetParameter(resname,val);}
//       Chi2_Hess GetCentral() const{return _para_central;}
       double FillContour();
       double SetContourPoint();
             
//       GetHess();
//       GetHess_Reduced();
//       GetCovariance();
//       GetCovariance_Reduced();
//       GetCentralValue();
};

void Solve_y_by_x(double y[2], double x, double delta, TMatrixD* inverse_cov);
void Solve_u_by_x(TMatrixD* u, TMatrixD * central_covariance ,double x, double y, double delta, TMatrixD * inverse_cov, int NRES, int idx_reduced, int idy_reduced);
} //End of namespace Res::Fit

/*
struct ResModel{
     vector<resonance> _ResVector;
     double _Chi2;
     TMatrixD _Covariance;
     double _Gradient;
     TMatrixD _EdgePoint;
     TGraph * _Contour; 
};
*/

#endif
