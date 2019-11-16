#include <stdio.h>
#include <math.h> 
using namespace std;

double derivativeX1(double x1, double x2){

  double response;

  double dividend;

  dividend = pow(x1,6)/6 - pow(x1,4) + 2*pow(x1,2) + x1*x2 + pow(x2,2);

  response = (pow(x1,5) - 4*pow(x1,3) + 4*x1 + x2)/(2*sqrt(dividend));

  return response;
}

double derivativeX2(double x1, double x2){

  double response;

  double dividend;

  dividend = pow(x1,6)/6 - pow(x1,4) + 2*pow(x1,2) + x1*x2 + pow(x2,2);

  response = (x1 + 2*x2 )/(2*sqrt(dividend));

  return response;
}
