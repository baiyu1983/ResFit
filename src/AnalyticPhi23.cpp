#define MPI 0.1395706
#define MJPSI 3.0969
#define ME 0.000511
#include "TMath.h"
#include "Math/SpecFuncMathMore.h"
#include "AnalyticPhi23.h"
#include <iostream>

using namespace std;

double AnalyticPhi23(double x, double m1,double m2,double m3){
      double Q_Plus  = (x+m1+m2+m3)*(x+m1-m2-m3)*(x-m1+m2-m3)*(x-m1-m2+m3);
      double Q_Minus = (x-m1-m2-m3)*(x-m1+m2+m3)*(x+m1-m2+m3)*(x+m1+m2-m3);
      double K=sqrt(Q_Minus/Q_Plus);
      double A1 = ( (x-m3)*(x-m3)-(m1+m2)*(m1+m2))/((x-m3)*(x-m3)-(m1-m2)*(m1-m2));
      double A2 = (m1-m2)*(m1-m2)*A1/((m1+m2)*(m1+m2));

      double Ek = ROOT::Math::comp_ellint_2(K);
      double Kk = ROOT::Math::comp_ellint_1(K);
      double PIk1 = ROOT::Math::comp_ellint_3(A1,K);
      double PIk2 = ROOT::Math::comp_ellint_3(A2,K);

      double result = 0.5*Q_Plus*(m1*m1+m2*m2+m3*m3+x*x)*Ek;
      result += 4*m1*m2*( (x-m3)*(x-m3) -(m1-m2)*(m1-m2) )*( (x+m3)*(x+m3) -m3*x+m1*m2  )*Kk;
      result += 8*m1*m2*( (m1*m1+m2*m2)*(x*x+m3*m3) -2*m1*m1*m2*m2 -2*m3*m3*x*x ) *PIk1;
      result -= 8*m1*m2*(x*x-m3*m3)*(x*x-m3*m3)*PIk2;
      result /= sqrt(Q_Plus);

      return result/(256*TMath::Pi()*TMath::Pi()*TMath::Pi()*x*x*x*sqrt(x*x-4*ME*ME));
}

double Coeff(double x, double m1, double m2, double m3){
      double Q_Plus  = (x+m1+m2+m3)*(x+m1-m2-m3)*(x-m1+m2-m3)*(x-m1-m2+m3);
      return 1.0/(256*TMath::Pi()*TMath::Pi()*TMath::Pi()*x*x*x*sqrt(x*x-4*ME*ME)*sqrt(Q_Plus));
}

double Coeff1(double x, double m1, double m2, double m3){  //Coefficient of 1-st type elliptic integral function
      return 4*m1*m2*( (x-m3)*(x-m3) -(m1-m2)*(m1-m2) )*( (x+m3)*(x+m3) -m3*x+m1*m2  );
}

double Coeff2(double x, double m1, double m2, double m3){ //Coefficient of 2-nd type elliptic integral function
      double Q_Plus = (x+m1+m2+m3)*(x+m1-m2-m3)*(x-m1+m2-m3)*(x-m1-m2+m3);
      return 0.5*Q_Plus*(m1*m1+m2*m2+m3*m3+x*x);
}

double Coeff3(double x, double m1, double m2, double m3){ ////Coefficient of 3-nd type elliptic integral function, alpha1
      double A1 = ( (x-m3)*(x-m3)-(m1+m2)*(m1+m2))/((x-m3)*(x-m3)-(m1-m2)*(m1-m2));
      double A2 = (m1-m2)*(m1-m2)*A1/((m1+m2)*(m1+m2));
      return 8*m1*m2*( (m1*m1+m2*m2)*(x*x+m3*m3) -2*m1*m1*m2*m2 -2*m3*m3*x*x );
}

double Coeff4(double x, double m1, double m2, double m3){ ////Coefficient of 3-nd type elliptic integral function, alpha2
      double A1 = ( (x-m3)*(x-m3)-(m1+m2)*(m1+m2))/((x-m3)*(x-m3)-(m1-m2)*(m1-m2));
      double A2 = (m1-m2)*(m1-m2)*A1/((m1+m2)*(m1+m2));
      return -8*m1*m2*(x*x-m3*m3)*(x*x-m3*m3);
}

double Coeff_derv1(double x, double m1, double m2, double m3){
      double Q_Plus  = (x+m1+m2+m3)*(x+m1-m2-m3)*(x-m1+m2-m3)*(x-m1-m2+m3);
      double result = 3.0/(x*x)+1.0/(x*x-4*ME*ME)+0.5*QPlus_derv1(x,m1,m2,m3)/(x*Q_Plus);
      result /= (sqrt(Q_Plus)*x*x*sqrt(x*x-4*ME*ME));
      result /= -256*TMath::Pi()*TMath::Pi()*TMath::Pi();
      return result;
}

double Coeff1_derv1(double x, double m1, double m2, double m3){
      double result = 2*(x-m3)*((x+m3)*(x+m3)-m3*x+m1*m2);
      result += ((x-m3)*(x-m3)-(m1-m2)*(m1-m2))*(2*x+m3);
      return result*4*m1*m2;
}

double Coeff2_derv1(double x, double m1, double m2, double m3){
      double Q_Plus = (x+m1+m2+m3)*(x+m1-m2-m3)*(x-m1+m2-m3)*(x-m1-m2+m3);
      double QP_d1 = QPlus_derv1(x,m1,m2,m3);
      return 0.5*(Q_Plus*2*x+QP_d1*(x*x+m1*m1+m2*m2+m3*m3));
}

double Coeff3_derv1(double x, double m1, double m2, double m3){
      double result = 2*x*(m1*m1+m2*m2);
      result -= 4*m3*m3*x;
      return result*8*m1*m2;
}

double Coeff4_derv1(double x, double m1, double m2, double m3){
      return -32*m1*m2*x*(x*x-m3*m3);
}


double Coeff_derv2(double x, double m1, double m2, double m3){
      double Q_Plus  = (x+m1+m2+m3)*(x+m1-m2-m3)*(x-m1+m2-m3)*(x-m1-m2+m3);
      double QPd1 = QPlus_derv1(x,m1,m2,m3);
      double result = (2/(x*x) + 1/(x*x-4*ME*ME)+QPd1/(2*x*Q_Plus));
      result *= (3/(x*x)+1.0/(x*x-4*ME*ME)+QPd1/(2*x*Q_Plus));
      result += (6/(x*x*x*x)+2/((x*x-4*ME*ME)*(x*x-4*ME*ME)) -QPlus_derv2(x,m1,m2,m3)/(2*x*x*Q_Plus) + QPd1*(1/(2*x*x*x*Q_Plus)+QPd1/(2*x*x*Q_Plus*Q_Plus)));
      result /= (sqrt(Q_Plus)*x*sqrt(x*x-4*ME*ME));
      return result/(256*TMath::Pi()*TMath::Pi()*TMath::Pi());
}

double Coeff1_derv2(double x, double m1, double m2, double m3){
      double result =0;
      result += 2*( (x+m3)*(x+m3) -m3*x+m1*m2 );
      result += 2*(x-m3)*(2*x+m3);
      result += 2*(x-m3)*(2*x+m3);
      result += 2*((x-m3)*(x-m3)-(m1-m2)*(m1-m2));
      return result*4*m1*m2;
}

double Coeff2_derv2(double x, double m1, double m2, double m3){
      double Q_Plus = (x+m1+m2+m3)*(x+m1-m2-m3)*(x-m1+m2-m3)*(x-m1-m2+m3);
      double result = 0;
      result += QPlus_derv2(x,m1,m2,m3)*(x*x+m1*m1+m2*m2+m3*m3);
      result += 2*QPlus_derv1(x,m1,m2,m3)*(2*x);
      result += Q_Plus*2;
      return result/2.0;
}

double Coeff3_derv2(double x, double m1, double m2, double m3){
     double result =0;
     result += (m1*m1+m2*m2)*2;
     result -= 2*m3*m3*2;
     return result*8*m1*m2;
}

double Coeff4_derv2(double x, double m1, double m2, double m3){
     return -8*m1*m2*(12*x*x-4*m3*m3);
}

double AnalyticPhi23_derv1(double x, double m1, double m2, double m3){
      double C =  Coeff(x, m1, m2, m3);
      double C1 = Coeff1(x,m1,m2,m3);
      double C2 = Coeff2(x,m1,m2,m3);
      double C3 = Coeff3(x,m1,m2,m3);
      double C4 = Coeff4(x,m1,m2,m3);
      double C_d1 =  Coeff_derv1(x,m1,m2,m3 );
      double C1_d1 = Coeff1_derv1(x,m1,m2,m3 ); 
      double C2_d1 = Coeff2_derv1(x,m1,m2,m3 );
      double C3_d1 = Coeff3_derv1(x,m1,m2,m3 );
      double C4_d1 = Coeff4_derv1(x,m1,m2,m3 );
      double E1_d1 = comp_ellint_1_derv1_over_X(x,m1,m2,m3);
      double E2_d1 = comp_ellint_2_derv1_over_X(x,m1,m2,m3);
      double E3_d1 = comp_ellint_3_derv1_A1_over_X(x,m1,m2,m3);
      double E4_d1 = comp_ellint_3_derv1_A2_over_X(x,m1,m2,m3);

      double Q_Plus  = (x+m1+m2+m3)*(x+m1-m2-m3)*(x-m1+m2-m3)*(x-m1-m2+m3);
      double Q_Minus = (x-m1-m2-m3)*(x-m1+m2+m3)*(x+m1-m2+m3)*(x+m1+m2-m3);
      double K=sqrt(Q_Minus/Q_Plus);
      double A1 = ( (x-m3)*(x-m3)-(m1+m2)*(m1+m2))/((x-m3)*(x-m3)-(m1-m2)*(m1-m2));
      double A2 = (m1-m2)*(m1-m2)*A1/((m1+m2)*(m1+m2));

      double E2 = ROOT::Math::comp_ellint_2(K);
      double E1 = ROOT::Math::comp_ellint_1(K);
      double E3 = ROOT::Math::comp_ellint_3(A1,K);
      double E4 = ROOT::Math::comp_ellint_3(A2,K);
      /*
      cout<<"C_d = "<<C_d1<<endl;
      cout<<"C1_d = "<<C1_d1<<endl;
      cout<<"C2_d = "<<C2_d1<<endl;
      cout<<"C3_d = "<<C3_d1<<endl;
      cout<<"C4_d = "<<C4_d1<<endl;
      cout<<"E1_d = "<<E1_d1<<endl;
      cout<<"E2_d = "<<E2_d1<<endl;
      cout<<"E3_d = "<<E3_d1<<endl;
      cout<<"E4_d = "<<E4_d1<<endl;
      cout<<"Phasespace = "<<C*(C1*E1+C2*E2+C3*E3+C4*E4)<<endl;
      */
      double result =  C_d1*(E1*C1+E2*C2+E3*C3+E4*C4);
      result += C*(E1*C1_d1+E1_d1*C1);
      result += C*(E2*C2_d1+E2_d1*C2);
      result += C*(E3*C3_d1+E3_d1*C3);
      result += C*(E4*C4_d1+E4_d1*C4);
      return result;
}

double AnalyticPhi23_derv2(double x, double m1, double m2, double m3){
      double C =  Coeff(x, m1, m2, m3);
      double C1 = Coeff1(x,m1,m2,m3);
      double C2 = Coeff2(x,m1,m2,m3);
      double C3 = Coeff3(x,m1,m2,m3);
      double C4 = Coeff4(x,m1,m2,m3);
      double C_d1 =  Coeff_derv1(x,m1,m2,m3 );
      double C1_d1 = Coeff1_derv1(x,m1,m2,m3 );
      double C2_d1 = Coeff2_derv1(x,m1,m2,m3 );
      double C3_d1 = Coeff3_derv1(x,m1,m2,m3 );
      double C4_d1 = Coeff4_derv1(x,m1,m2,m3 );
      double E1_d1 = comp_ellint_1_derv1_over_X(x,m1,m2,m3);
      double E2_d1 = comp_ellint_2_derv1_over_X(x,m1,m2,m3);
      double E3_d1 = comp_ellint_3_derv1_A1_over_X(x,m1,m2,m3);
      double E4_d1 = comp_ellint_3_derv1_A2_over_X(x,m1,m2,m3);

      double Q_Plus  = (x+m1+m2+m3)*(x+m1-m2-m3)*(x-m1+m2-m3)*(x-m1-m2+m3);
      double Q_Minus = (x-m1-m2-m3)*(x-m1+m2+m3)*(x+m1-m2+m3)*(x+m1+m2-m3);
      double K=sqrt(Q_Minus/Q_Plus);
      double A1 = ( (x-m3)*(x-m3)-(m1+m2)*(m1+m2))/((x-m3)*(x-m3)-(m1-m2)*(m1-m2));
      double A2 = (m1-m2)*(m1-m2)*A1/((m1+m2)*(m1+m2));

      double E2 = ROOT::Math::comp_ellint_2(K);
      double E1 = ROOT::Math::comp_ellint_1(K);
      double E3 = ROOT::Math::comp_ellint_3(A1,K);
      double E4 = ROOT::Math::comp_ellint_3(A2,K);

      double C_d2 =  Coeff_derv2(x,m1,m2,m3 );
      double C1_d2 = Coeff1_derv2(x,m1,m2,m3 );
      double C2_d2 = Coeff2_derv2(x,m1,m2,m3 );
      double C3_d2 = Coeff3_derv2(x,m1,m2,m3 );
      double C4_d2 = Coeff4_derv2(x,m1,m2,m3 );
      double E1_d2 = comp_ellint_1_derv2_over_X(x,m1,m2,m3);
      double E2_d2 = comp_ellint_2_derv2_over_X(x,m1,m2,m3);
      double E3_d2 = comp_ellint_3_derv2_A1_over_X(x,m1,m2,m3);
      double E4_d2 = comp_ellint_3_derv2_A2_over_X(x,m1,m2,m3);

      double result = 0;
      result += C_d2*(E1*C1+E2*C2+E3*C3+E4*C4);
      result += 2*C_d1*(E1_d1*C1+E1*C1_d1);
      result += 2*C_d1*(E2_d1*C2+E2*C2_d1);
      result += 2*C_d1*(E3_d1*C3+E3*C3_d1);
      result += 2*C_d1*(E4_d1*C4+E4*C4_d1);
      result += C*(E1_d2*C1+2*E1_d1*C1_d1+E1*C1_d2);
      result += C*(E2_d2*C2+2*E2_d1*C2_d1+E2*C2_d2);
      result += C*(E3_d2*C3+2*E3_d1*C3_d1+E3*C3_d2);
      result += C*(E4_d2*C4+2*E4_d1*C4_d1+E4*C4_d2);
      return result;
}

double QPlus_derv1(double x, double m1, double m2, double m3){
      double QPlus = (x+m1+m2+m3)*(x+m1-m2-m3)*(x-m1+m2-m3)*(x-m1-m2+m3);
      double r = 1.0/(x+m1+m2+m3)+1.0/(x+m1-m2-m3)+1.0/(x-m1+m2-m3)+1.0/(x-m1-m2+m3);
      return r*QPlus;
}

double QPlus_derv2(double x, double m1, double m2, double m3){
      double q1 = x+m1+m2+m3;  double q2 = x+m1-m2-m3;  double q3 = x-m1+m2-m3; double q4 = x-m1-m2+m3;
      return 2*(q1*q2+q1*q3+q1*q4+q2*q3+q2*q4+q3*q4);
}

double QMinus_derv1(double x, double m1, double m2, double m3){
      double QMinus = (x-m1-m2-m3)*(x-m1+m2+m3)*(x+m1-m2+m3)*(x+m1+m2-m3);
      double r = 1.0/(x-m1-m2-m3)+1.0/(x-m1+m2+m3)+1.0/(x+m1-m2+m3)+1.0/(x+m1+m2-m3);
      return r*QMinus;
}

double QMinus_derv2(double x, double m1, double m2, double m3){
      double q1 = x-m1-m2-m3;  double q2 = x-m1+m2+m3; double q3 = x+m1-m2+m3; double q4 = x+m1+m2-m3;
      return 2*(q1*q2+q1*q3+q1*q4+q2*q3+q2*q4+q3*q4);
}

double K_derv1(double x, double m1, double m2, double m3){
      double Q_Plus  = (x+m1+m2+m3)*(x+m1-m2-m3)*(x-m1+m2-m3)*(x-m1-m2+m3);
      double Q_Minus = (x-m1-m2-m3)*(x-m1+m2+m3)*(x+m1-m2+m3)*(x+m1+m2-m3);
      return 0.5*sqrt(Q_Plus/Q_Minus)*(QMinus_derv1(x,m1,m2,m3)*Q_Plus-QPlus_derv1(x,m1,m2,m3)*Q_Minus)/(Q_Plus*Q_Plus);
}

double K_derv2(double x, double m1, double m2, double m3){
      double Q_Plus  = (x+m1+m2+m3)*(x+m1-m2-m3)*(x-m1+m2-m3)*(x-m1-m2+m3);
      double Q_Minus = (x-m1-m2-m3)*(x-m1+m2+m3)*(x+m1-m2+m3)*(x+m1+m2-m3);
      double q_plus_derv1 = QPlus_derv1(x,m1,m2,m3);
      double q_minus_derv1 = QMinus_derv1(x,m1,m2,m3);
      double q_plus_derv2 = QPlus_derv2(x,m1,m2,m3);
      double q_minus_derv2 = QMinus_derv2(x,m1,m2,m3);
      double K=sqrt(Q_Minus/Q_Plus);      
      double result = 0.5*K*(q_plus_derv1*Q_Minus-q_minus_derv1*Q_Plus)*(q_minus_derv1*Q_Plus-q_plus_derv1*Q_Minus)/(Q_Plus*Q_Plus*Q_Minus*Q_Minus);
      result += ( (q_minus_derv2*Q_Plus-q_plus_derv2*Q_Minus)*(Q_Plus*Q_Plus) -2*Q_Plus*q_plus_derv1*(q_minus_derv1*Q_Plus-q_plus_derv1*Q_Minus))/(Q_Plus*Q_Plus*Q_Plus*Q_Plus*K);
      return result/2.0;
}

double A1_derv1(double x, double m1, double m2, double m3){ 
      double d = (x-m3)*(x-m3)-(m1-m2)*(m1-m2);
      return 2*4*m1*m2*(x-m3)/(d*d);
}

double A1_derv2(double x, double m1, double m2, double m3){
       double d = (x-m3)*(x-m3)-(m1-m2)*(m1-m2);
       double n1 = 2*d*d;
       double n2 = 2*(x-m3)*2*d*2*(x-m3);
       return 4*m1*m2*(n1-n2)/(d*d*d*d);
}

double A2_derv1(double x, double m1, double m2, double m3){
      double r = (m1-m2)*(m1-m2)/((m1+m2)*(m1+m2));
      if(r==0)return 0;
      else return r*A1_derv1(x,m1,m2,m3);
}

double A2_derv2(double x, double m1, double m2, double m3){
      double r = (m1-m2)*(m1-m2)/((m1+m2)*(m1+m2));
      if(r==0)return 0;
      else return r*A1_derv2(x,m1,m2,m3);
}
	
double comp_ellint_1_derv1(double x){
      return ROOT::Math::comp_ellint_2(x)/(x*(1-x*x))-ROOT::Math::comp_ellint_1(x)/x;
}

double comp_ellint_1_derv2(double x){
      return (x*ROOT::Math::comp_ellint_1(x)-(1-3*x*x)*comp_ellint_1_derv1(x))/(x*(1-x*x));
}	

double comp_ellint_2_derv1(double x){
      return (ROOT::Math::comp_ellint_2(x) - ROOT::Math::comp_ellint_1(x))/x;
}	

double comp_ellint_2_derv2(double x){
      return ( x*ROOT::Math::comp_ellint_2(x)/(x*x-1) -comp_ellint_2_derv1(x)  )/x;
}	

double comp_ellint_3_derv1_n(double y, double x){
      return (  ROOT::Math::comp_ellint_2(x) +(x*x-y)*ROOT::Math::comp_ellint_1(x)/y + (y*y-x*x)*ROOT::Math::comp_ellint_3(y,x)/y  )/( 2*(x*x-y) * (y-1) )	;
}

double comp_ellint_3_derv1_k(double y, double x){
      return (ROOT::Math::comp_ellint_2(x)/(x*x-1) + ROOT::Math::comp_ellint_3(y,x))*x/(y-x*x);	 
}

double comp_ellint_3_derv2_kk(double y, double x){
      double z1 = (y+x*x)*(ROOT::Math::comp_ellint_2(x)/(x*x-1) + ROOT::Math::comp_ellint_3(y,x))/((y-x*x)*(y-x*x));
      double z2 = comp_ellint_2_derv1(x)/(x*x-1) -ROOT::Math::comp_ellint_2(x)*2*x/((x*x-1)*(x*x-1));
      double z3 = comp_ellint_3_derv1_k(y, x);	      
      return z1+x*(z2+z3)/(y-x*x);
}

double comp_ellint_3_derv2_kn(double y, double x){ 
      return -comp_ellint_3_derv1_k(y,x)/(y-x*x)+x*comp_ellint_3_derv1_n(y,x)/(y-x*x);	
}

double comp_ellint_3_derv2_nk(double y, double x){
      return comp_ellint_3_derv2_kn(y,x);	
}

double comp_ellint_3_derv2_nn(double y, double x){
      double c1 = ROOT::Math::comp_ellint_3(y,x);
      double d1 = comp_ellint_3_derv1_n(y,x);
      double r1 = 1.0/(x*x-y)-1.0/(y-1);
      double z1 = d1*r1;
      double z2 = -x*x*ROOT::Math::comp_ellint_1(x)/(y*y);
      double r3 = (x*x+y*y)/(y*y);
      return z1+(z2+c1*r3+d1*(y*y-x*x)/y)/(2*(x*x-y)*(y-1)); 	      
}

double comp_ellint_1_derv1_over_X(double x, double m1, double m2, double m3){
      double Q_Plus  = (x+m1+m2+m3)*(x+m1-m2-m3)*(x-m1+m2-m3)*(x-m1-m2+m3);
      double Q_Minus = (x-m1-m2-m3)*(x-m1+m2+m3)*(x+m1-m2+m3)*(x+m1+m2-m3);
      double K=sqrt(Q_Minus/Q_Plus);
      return comp_ellint_1_derv1(K)*K_derv1(x,m1,m2,m3);
}

double comp_ellint_1_derv2_over_X(double x, double m1, double m2, double m3){
      double Q_Plus  = (x+m1+m2+m3)*(x+m1-m2-m3)*(x-m1+m2-m3)*(x-m1-m2+m3);
      double Q_Minus = (x-m1-m2-m3)*(x-m1+m2+m3)*(x+m1-m2+m3)*(x+m1+m2-m3);
      double K=sqrt(Q_Minus/Q_Plus);
      double kderv1 = K_derv1(x,m1,m2,m3);
      return comp_ellint_1_derv2(K)*kderv1*kderv1+comp_ellint_1_derv1(K)*K_derv2(x,m1,m2,m3);
}

double comp_ellint_2_derv1_over_X(double x, double m1, double m2, double m3){
      double Q_Plus  = (x+m1+m2+m3)*(x+m1-m2-m3)*(x-m1+m2-m3)*(x-m1-m2+m3);
      double Q_Minus = (x-m1-m2-m3)*(x-m1+m2+m3)*(x+m1-m2+m3)*(x+m1+m2-m3);
      double K=sqrt(Q_Minus/Q_Plus);
      return comp_ellint_2_derv1(K)*K_derv1(x,m1,m2,m3);
}

double comp_ellint_2_derv2_over_X(double x, double m1, double m2, double m3){
      double Q_Plus  = (x+m1+m2+m3)*(x+m1-m2-m3)*(x-m1+m2-m3)*(x-m1-m2+m3);
      double Q_Minus = (x-m1-m2-m3)*(x-m1+m2+m3)*(x+m1-m2+m3)*(x+m1+m2-m3);
      double K=sqrt(Q_Minus/Q_Plus);
      double kderv1 = K_derv1(x,m1,m2,m3);
      return comp_ellint_2_derv2(K)*kderv1*kderv1+comp_ellint_2_derv1(K)*K_derv2(x,m1,m2,m3);
}

double comp_ellint_3_derv1_A1_over_X(double x, double m1, double m2, double m3){
      double Q_Plus  = (x+m1+m2+m3)*(x+m1-m2-m3)*(x-m1+m2-m3)*(x-m1-m2+m3);
      double Q_Minus = (x-m1-m2-m3)*(x-m1+m2+m3)*(x+m1-m2+m3)*(x+m1+m2-m3);
      double K=sqrt(Q_Minus/Q_Plus);
      double A1 = ( (x-m3)*(x-m3)-(m1+m2)*(m1+m2))/((x-m3)*(x-m3)-(m1-m2)*(m1-m2));
      return comp_ellint_3_derv1_n(A1,K)*A1_derv1(x,m1,m2,m3)+comp_ellint_3_derv1_k(A1,K)*K_derv1(x,m1,m2,m3);
}

double comp_ellint_3_derv1_A2_over_X(double x, double m1, double m2, double m3){ 
      double Q_Plus  = (x+m1+m2+m3)*(x+m1-m2-m3)*(x-m1+m2-m3)*(x-m1-m2+m3);
      double Q_Minus = (x-m1-m2-m3)*(x-m1+m2+m3)*(x+m1-m2+m3)*(x+m1+m2-m3);
      double K=sqrt(Q_Minus/Q_Plus);
      double A1 = ( (x-m3)*(x-m3)-(m1+m2)*(m1+m2))/((x-m3)*(x-m3)-(m1-m2)*(m1-m2));
      double A2 = (m1-m2)*(m1-m2)*A1/((m1+m2)*(m1+m2));
      /*
      cout<<"In function comp_ellint_3_derv1_A2_over_X, A2 = "<<A2<<endl;
      */ 
      if(A2==0){
	      /*
	       cout<<"In function comp_ellint_3_derv1_A2_over_X, A2_derv1 = "<<A2_derv1(x,m1,m2,m3)<<endl;
	       cout<<" In function comp_ellint_3_derv1_A2_over_X, comp_ellint_3_derv1_n(A2,K) = "<<comp_ellint_3_derv1_n(0,K)<<endl;
	      */ 
	       return comp_ellint_3_derv1_k(0,K)*K_derv1(x,m1,m2,m3);
      }
      else return comp_ellint_3_derv1_n(A2,K)*A2_derv1(x,m1,m2,m3)+comp_ellint_3_derv1_k(A2,K)*K_derv1(x,m1,m2,m3);
    
}

double comp_ellint_3_derv2_A1_over_X(double x, double m1, double m2, double m3){
      double Q_Plus  = (x+m1+m2+m3)*(x+m1-m2-m3)*(x-m1+m2-m3)*(x-m1-m2+m3);
      double Q_Minus = (x-m1-m2-m3)*(x-m1+m2+m3)*(x+m1-m2+m3)*(x+m1+m2-m3);
      double K=sqrt(Q_Minus/Q_Plus);
      double A1 = ( (x-m3)*(x-m3)-(m1+m2)*(m1+m2))/((x-m3)*(x-m3)-(m1-m2)*(m1-m2));
      double DKDx = K_derv1(x,m1,m2,m3);
      double DA1Dx = A1_derv1(x, m1, m2, m3);
      return comp_ellint_3_derv1_n(A1,K)*A1_derv2(x,m1,m2,m3) + comp_ellint_3_derv1_k(A1,K)*K_derv2(x,m1,m2,m3) + comp_ellint_3_derv2_nn(A1,K)*DA1Dx*DA1Dx + \
	      comp_ellint_3_derv2_kk(A1,K)*DKDx*DKDx+2*comp_ellint_3_derv2_kn(A1,K)*DKDx*DA1Dx; 
}

double comp_ellint_3_derv2_A2_over_X(double x, double m1, double m2, double m3){
      double Q_Plus  = (x+m1+m2+m3)*(x+m1-m2-m3)*(x-m1+m2-m3)*(x-m1-m2+m3);
      double Q_Minus = (x-m1-m2-m3)*(x-m1+m2+m3)*(x+m1-m2+m3)*(x+m1+m2-m3);
      double K=sqrt(Q_Minus/Q_Plus);
      double A1 = ( (x-m3)*(x-m3)-(m1+m2)*(m1+m2))/((x-m3)*(x-m3)-(m1-m2)*(m1-m2));
      double A2 = (m1-m2)*(m1-m2)*A1/((m1+m2)*(m1+m2));
      double DKDx = K_derv1(x,m1,m2,m3);
      double DA2Dx = A2_derv1(x, m1, m2, m3); 
      if(A2 ==0)return comp_ellint_3_derv1_k(0,K)*K_derv2(x,m1,m2,m3) + comp_ellint_3_derv2_kk(0,K)*DKDx*DKDx;
      else return comp_ellint_3_derv1_n(A2,K)*A2_derv2(x,m1,m2,m3) + comp_ellint_3_derv1_k(A2,K)*K_derv2(x,m1,m2,m3) + comp_ellint_3_derv2_nn(A2,K)*DA2Dx*DA2Dx + \
              comp_ellint_3_derv2_kk(A2,K)*DKDx*DKDx+2*comp_ellint_3_derv2_kn(A2,K)*DKDx*DA2Dx;
}


