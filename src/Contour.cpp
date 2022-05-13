#include "Header.h"
#include "Contour.h"
namespace ResFit{

ResFit::Contour::Contour(const ResFit::Contour & c){ //Copy Constructor, name not set;
    _nfp = c.GetNFP();
    _x_name = c.GetXName();
    _y_name = c.GetYName();
    std::map<int,std::pair<double,double> > data=c.GetData();
    for(std::map<int,std::pair<double,double> >::const_iterator itr= data.begin(); itr!=data.end();++itr){
         _data.insert(std::pair<int,std::pair<double,double> >((itr->first),std::pair<double,double>((itr->second).first,(itr->second).second)));
    }
}

void ResFit::Contour::operator=(ResFit::Contour c){
    _nfp = c.GetNFP();
    _x_name = c.GetXName();
    _y_name = c.GetYName();
    std::map<int,std::pair<double,double> > data=c.GetData();
    for(std::map<int,std::pair<double,double> >::const_iterator itr= data.begin(); itr!=data.end();++itr){
         _data.insert(std::pair<int,std::pair<double,double> >((itr->first),std::pair<double,double>((itr->second).first,(itr->second).second)));
    }
}

void ResFit::Contour::Insert(int n, double x, double y){
    _data.insert(std::pair<int, std::pair<double,double> >(n,std::pair<double,double>(x,y)) );
}

void ResFit::Contour::SetData(int n,double x,double y){
    std::map<int,std::pair<double,double> >::const_iterator itr = _data.find(n);
    if(itr == _data.end()){
         Insert(n,x,y);
    }else{
         _data[n] = std::pair<double,double>(x,y);
    }
}

}
