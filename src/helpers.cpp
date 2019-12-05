#include <stdio.h>
#include <iostream>
#include <math.h> 
#include <vector>
#include "helpers.h"

using namespace std;


double function(double x1, double x2)
{
  // cout << endl << "func value: " << sqrt(pow(x1,6)/6 - pow(x1,4) + 2*pow(x1,2) + x1*x2 + pow(x2,2)) << endl;
  return sqrt(pow(x1,6)/6 - pow(x1,4) + 2*pow(x1,2) + x1*x2 + pow(x2,2));
}

double phi(double x1, double x2, double d1, double d2, double t)
{
  return function(x1 + d1*t, x2 + d2*t);
}

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

  response = (x1 + 2*x2)/(2*sqrt(dividend));

  //cout << "dx2 " << response << endl;

  return response;
}

double hessianX1X1(double x1, double x2){
  double numerator1, numerator2;

  double dividend1, dividend2;

  numerator1 = 5*pow(x1,4) - 12*pow(x1,2) + 4;

  numerator2 = pow((pow(x1,5) - 4*pow(x1,3) + 4*x1 + x2),2);

  dividend1 = 2*(function(x1,x2));

  dividend2 = 4*pow(function(x1,x2) ,(3));

  // cout << "hessianx1x1: " << numerator1/dividend1 - numerator2/dividend2 << endl;

  return numerator1/dividend1 - numerator2/dividend2;
}

double hessianX1X2(double x1, double x2){
  double numerator1, numerator2;

  double dividend1, dividend2;

  numerator1 = 1;

  numerator2 = (pow(x1,5) - 4*pow(x1,3) + 4*x1 + x2)*(x1 + 2*x2);

  dividend1 = 2*(function(x1,x2));

  dividend2 = 4*pow(function(x1,x2) ,(3));

  // cout << "hessianx1x2: " << numerator1/dividend1 - numerator2/dividend2 << endl;

  return numerator1/dividend1 - numerator2/dividend2;
}

double hessianX2X1(double x1, double x2){
  double numerator1, numerator2;

  double dividend1, dividend2;
  
  numerator1 = 1;

  numerator2 = (pow(x1,5) - 4*pow(x1,3) + 4*x1 + x2)*(x1 + 2*x2);

  dividend1 = 2*(function(x1,x2));

  dividend2 = 4*pow(function(x1,x2) ,(3));

  // cout << "hessianx2x1: " << numerator1/dividend1 - numerator2/dividend2 << endl;

  return numerator1/dividend1 - numerator2/dividend2;
}

double hessianX2X2(double x1, double x2){
  double numerator1, numerator2;

  double dividend1, dividend2;

  numerator1 = 1;

  numerator2 = pow((x1 + 2*x2),2);

  dividend1 = (function(x1,x2));

  dividend2 = 4*pow(function(x1,x2) ,(3));

  // cout << "hessianx2x2: " << numerator1/dividend1 - numerator2/dividend2 << endl;

  return numerator1/dividend1 - numerator2/dividend2;
}

vector<double> dnewton(double x1, double x2,double gradX1, double gradX2){
  double x1x2 = hessianX1X2(x1,x2);
  double x1x1 = hessianX1X1(x1,x2);
  double x2x2 = hessianX2X2(x1,x2);
  double x2x1 = hessianX2X1(x1,x2);

  double hessian[2][2] = {{x1x1,x1x2},{x2x1,x2x2}};

  double detHessian = hessian[0][0]*hessian[1][1] - (hessian[1][0]*hessian[0][1]);

  double inv_hessian[2][2] = {{hessian[1][1]/detHessian,-hessian[0][1]/detHessian},{-hessian[1][0]/detHessian, hessian[0][0]/detHessian}};

  double gradient[2] = {derivativeX1(x1,x2) ,derivativeX2(x1,x2)};

  vector<double> d = {-(inv_hessian[0][0]*gradX1 + inv_hessian[0][1]*gradX2),-(inv_hessian[1][0]*gradX1 + inv_hessian[1][1]*gradX2)};
  
  return d;
}

vector<double> p(double x1,double x2,double x1_1, double x2_1){
  vector<double> p = {x1_1 - x1, x2_1 - x2};
  return p;
}

vector<double> q(double x1,double x2,double x1_1, double x2_1){
  vector<double> q = { derivativeX1(x1_1,x2_1) - derivativeX1(x1,x2),derivativeX2(x1_1,x2_1) - derivativeX2(x1,x2) };
  return q;
}

vector<vector<double>> BFGS(vector<vector<double>> H, double x1,double x2,double x1_1, double x2_1){
  // Hk+1 = firstpart + secondpart + thirdpart
  // firstpart = Hk

  double tolerance = 0.0000001;

  vector<vector<double>> secondpart(2,vector<double>(2, 0));
  vector<vector<double>> thirdpart(2,vector<double>(2, 0));
  vector<vector<double>> result(2,vector<double>(2, 0));

  vector<double> pk = p(x1,x2,x1_1,x2_1);
  vector<double> qk = q(x1,x2,x1_1,x2_1);

  double dividend = pk[0]*qk[0] + pk[1]*qk[1];

  if (dividend < tolerance) {
    return H;
  }

  // double parantheses = 1 + (H[0][0]*qk[0]*qk[0] + H[1][0]*qk[0]*qk[0] + H[0][1]*qk[1]*qk[1] + H[1][1]*qk[1]*qk[1])/dividend;
  double parantheses = 1 + ((qk[0] * H[0][0] + qk[1] * H[1][0]) * qk[0] + (qk[0] * H[0][1] + qk[1] * H[1][1])*qk[1])/dividend;

  secondpart[0][0] = (parantheses/dividend)*(pk[0]*pk[0]);
  secondpart[0][1] = (parantheses/dividend)*(pk[0]*pk[1]);
  secondpart[1][0] = (parantheses/dividend)*(pk[1]*pk[0]);
  secondpart[1][1] = (parantheses/dividend)*(pk[1]*pk[1]);

  // cout << "dividend: " << dividend << endl;
  thirdpart[0][0] = (pk[0]*qk[0]*H[0][0] + pk[0]*qk[1]*H[1][0] + (H[0][0]*qk[0] + H[0][1]*qk[1])*pk[0])/dividend;
  thirdpart[0][1] = (pk[0]*qk[0]*H[0][1] + pk[0]*qk[1]*H[1][1] + (H[0][0]*qk[0] + H[0][1]*qk[1])*pk[1])/dividend;
  thirdpart[1][0] = (pk[1]*qk[0]*H[0][0] + pk[1]*qk[1]*H[1][0] + (H[1][0]*qk[0] + H[1][1]*qk[1])*pk[0])/dividend;
  thirdpart[1][1] = (pk[1]*qk[0]*H[0][1] + pk[1]*qk[1]*H[1][1] + (H[1][0]*qk[0] + H[1][1]*qk[1])*pk[1])/dividend;

  
  result[0][0] = H[0][0] + secondpart[0][0] - thirdpart[0][0];
  result[0][1] = H[0][1] + secondpart[0][1] - thirdpart[0][1];
  result[1][0] = H[1][0] + secondpart[1][0] - thirdpart[1][0];
  result[1][1] = H[1][1] + secondpart[1][1] - thirdpart[1][1];

  return result;
}

vector<vector<double>> inverse(vector<vector<double>> H){

  vector<vector<double>> result(2,vector<double>(2, 0));
  
  double det = H[0][0]*H[1][1] - H[1][0]*H[0][1];

  result[0][0] = H[1][1]/det;
  result[0][1] = 0 - H[0][1]/det;
  result[1][0] = 0 - H[1][0]/det;
  result[1][1] = H[0][0]/det;
  
  return result;

}


