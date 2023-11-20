#include <string>
#include <iomanip>
#include <cmath>
#include <vector>
#include "eos.h"

using std::string;
using std::scientific;
using std::setw;
using std::setprecision;
using std::endl;
using std::ofstream;
using std::ostringstream;


EquationOfState:: EquationOfState():
    temporary(),
    //S2E(std::move(temporary.discreteEntropy), std::move(temporary.discreteEnergy)),
    //E2S(std::move(temporary.discreteEnergy), std::move(temporary.discreteEntropy))
S2E(temporary.discreteEntropy, temporary.discreteEnergy),
E2S(temporary.discreteEnergy, temporary.discreteEntropy)
{
    
}

  
//EquationOfState::~EquationOfState(){
    
//}
    

double EquationOfState::EtoS( double E)
{
    
    return E2S(E);
    
}

double EquationOfState::StoE( double S)
{
    return S2E(S);
    
}

double EquationOfState::energy (double T)
{
    return temporary.energy(T);
}

double EquationOfState::entropy (double T)
{
    return temporary.entropy(T);
}

discrete::discrete():discreteTemperature(NumberOfTemperatureGridPoints,0), discreteEntropy(NumberOfTemperatureGridPoints,0), discreteEnergy(NumberOfTemperatureGridPoints,0){
    for (std::size_t i=0; i< NumberOfTemperatureGridPoints ; i++) {
        
        discreteTemperature[i]= Tfin*(1.0*i)/((1.0*NumberOfTemperatureGridPoints) -1.0);
    
        discreteEnergy[i]=energy(discreteTemperature[i]);
       
    
        discreteEntropy[i]=entropy(discreteTemperature[i]);
        
       
    }
    
};


double discrete::energy( double T)
{
    //double E=exp(1);
    double result=0;
     //double Pi=acos(-1); // CHANGE
     if ( T==0 ){
         
          result=0;
    } else{
        result= (exp((-0.00009123176015137558 -0.1095231425037804*T)/pow(T,2))*pow(T,2)*(3.831255349484765e-8 + 0.000021400588024804855*T - 0.00029690569178115957*pow(T,2) -0.009695660474441342*pow(T,3) + 0.3647383055723698*pow(T,4) - 5.404720030734935*pow(T,5) + 46.26425935097975*pow(T,6) -251.70725781035893*pow(T,7) + 889.2684561328422*pow(T,8) - 1932.746777393406*pow(T,9) + 
2034.1450534655232*pow(T,10)))/pow(0.0010461910330715387 - 
0.016997351850612453*T + 0.12595528994069613*pow(T,2) - 
0.510477729039857*pow(T,3) + pow(T,4),2);
        
        
//        (exp((-0.0000912318 - 
//  0.109523*T)/pow(T,2))*pow(T,2)*(3.83126e-8 + 0.0000214006*T - 
 //  0.000296906*pow(T,2) - 0.00969566*pow(T,3) + 0.364738*pow(T,4) - 5.40472*pow(T,5) + 
  // 46.2643*pow(T,6) - 251.707*pow(T,7) + 889.268*pow(T,8) - 1932.75*pow(T,9) + 
//   2034.15*pow(T,10)))/pow((0.00104619 - 0.0169974*T + 0.125955*pow(T,2) - 
 // 0.510478*pow(T,2) + pow(T,4)),2);
        
    }
    return(result);
    

}


 double discrete::entropy( double T)
{
     double result=0;
   //  double E=exp(1);
   //  double Pi=acos(-1); // CHANGE
    
    if ( T==0 ){
         
          result=0;
    } else{
        result=(exp((-0.00009123176015137558 - 
0.1095231425037804*T)/pow(T,2))*T*(3.831255349484765e-8 + 
0.000021400588024804855*T - 0.00008693192041427013*pow(T,2) - 
0.018444802527704172*pow(T,3) + 0.5368632648594897*pow(T,4) - 
7.45225658178704*pow(T,5) + 62.30776783074574*pow(T,6) - 
336.8744183413885*pow(T,7) + 1191.0264050273945*pow(T,8) - 
2590.122576707207*pow(T,9) + 
2712.1934046206975*pow(T,10)))/pow(0.0010461910330715387 - 
0.016997351850612453*T + 0.12595528994069613*pow(T,2) - 
0.510477729039857*pow(T,3) + pow(T,4),2);
        
        
//        (Power(E,(-0.00009123176015137558 - 0.1095231425037804*T)/Power(T,2))*T*
  //   (0.00000003831255349484765 + 0.000021400588024804855*T - 0.00008693192041427013*Power(T,2) -
  //     0.018444802527704172*Power(T,3) + 0.5368632648594897*Power(T,4) - 7.45225658178704*Power(T,5) +
//       62.30776783074574*Power(T,6) - 336.8744183413885*Power(T,7) + 1191.0264050273945*Power(T,8) -
//       2590.122576707207*Power(T,9) + 2712.1934046206975*Power(T,10)))/
 //  Power(0.0010461910330715387 - 0.016997351850612453*T + 0.12595528994069613*Power(T,2) - 0.510477729039857*Power(T,3) +
   //  Power(T,4),2);
        
    }
    return(result);
    
}

 double discrete::Power( double a, int b)
{
    
    return(pow(a,b));
    
}
