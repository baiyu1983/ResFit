#ifndef CONTOUR_H
#define CONTOUR_H
#include "Header.h"

namespace ResFit{

class Contour{
      private:
         int _nfp;
         std::map<int,std::pair<double,double> > _data;
         std::string _name;
         std::string _x_name;
         std::string _y_name;
         double _z;
      public:
         Contour():_nfp(0),_data(std::map<int,std::pair<double,double>  > ()),_name(""),_x_name(""),_y_name(""){}
         Contour(const ResFit::Contour & a_contour); //Copy Reference Constructor, name of contour not set
         void operator = (ResFit::Contour );
         void SetNFP(int nfp){_nfp = nfp;}
         double GetZ()const{return _z;}
         int GetNFP() const{return _nfp;}
         void SetName(std::string name){_name=name;}
         std::string GetName() const{return _name;}
         void SetXName(std::string name){_x_name = name;}
         void SetYName(std::string name){_y_name = name;}
         void SetZ(double z){_z = z;}
         std::string GetXName() const{return _x_name;}
         std::string GetYName() const{return _y_name;}
         void Insert(int n, double x, double y);
//         std::pair<double,double> GetData(int n) const{return _data[n];}
         std::map<int,std::pair<double,double> > GetData() const{return _data;}
         void ClearData(){_data.clear();}
         void SetData(int n,double x, double y);
         void Sort();
         void Draw();
};

}
#endif
