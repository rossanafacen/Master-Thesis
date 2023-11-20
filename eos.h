#include <string>
#include <iomanip>
#include <cmath>
#include <vector>
//#include <boost/math/interpolators/barycentric_rational.hpp>
#include <iostream>
#include "spline.h"

using std::string;
using std::scientific;
using std::setw;
using std::setprecision;
using std::endl;
using std::ofstream;
using std::ostringstream;


struct discrete{
    
    const std::size_t NumberOfTemperatureGridPoints=100;
    const double Tfin=1.0;
    
    std::vector< double> discreteTemperature, discreteEntropy, discreteEnergy;
    
    
     double Power( double, int);
    
    
    
    discrete();
    
    double energy( double);
    double entropy( double);
    
    
    
};


class EquationOfState {
    
public: 
    
    
    EquationOfState();
  //  ~EquationOfState();
    
     double energy( double);
     double entropy( double);
    
    
    
     double StoE( double);
     double EtoS( double);

private:
    
    discrete temporary;
    
    tk::spline S2E;
    tk::spline E2S;
    //boost::math::barycentric_rational< double>  S2E;
    //boost::math::barycentric_rational< double>  E2S;
    
 
};
