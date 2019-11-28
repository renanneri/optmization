#include <stdio.h>
#include <iostream>
#include <math.h> 
#include <vector>
#include "otm.h"

using namespace std;

Optimization::Optimization(){

}

Optimization::~Optimization(){

}

double Optimization::function(double x1, double x2)
{
  return sqrt(pow(x1,6)/6 - pow(x1,4) + 2*pow(x1,2) + x1*x2 + pow(x2,2));
}

double Optimization::phi(double x1, double x2, double d1, double d2, double t)
{
  return function(x1 + d1*t, x2 + d2*t);
}

double Optimization::derivativeX1(double x1, double x2){

  double response;

  double dividend;

  dividend = pow(x1,6)/6 - pow(x1,4) + 2*pow(x1,2) + x1*x2 + pow(x2,2);
 
  response = (pow(x1,5) - 4*pow(x1,3) + 4*x1 + x2)/(2*sqrt(dividend));

  return response;
}

double Optimization::derivativeX2(double x1, double x2){

  double response;

  double dividend;

  dividend = pow(x1,6)/6 - pow(x1,4) + 2*pow(x1,2) + x1*x2 + pow(x2,2);

  response = (x1 + 2*x2)/(2*sqrt(dividend));

  //cout << "dx2 " << response << endl;

  return response;
}

double Optimization::hessianX1X1(double x1, double x2){
  double response;

  double dividend;

  dividend = pow(((x1,6) - 6*pow(x1,4) + 12*pow(x1,2) + 6*x1*x2 + 6*pow(x2,2)),(3/2)); 

  response = (sqrt(3/2)*((5*pow(x1,4) - 12*pow(x1,2) + 4)*(pow(x1,6) - 6*pow(x1,4) + 12*pow(x1,2) + 6*x1*x2 + 6*pow(x2,2)) - 3*pow((pow(x1,5) - 4*pow(x1,3) + 4*x1 + x2),2)))/dividend;

  return response;
}
double Optimization::hessianX1X2(double x1, double x2){
  double response;

  double dividend;

  dividend = pow(((x1,6) - 6*pow(x1,4) + 12*pow(x1,2) + 6*x1*x2 + 6*pow(x2,2)),(3/2));

  response = (sqrt(3/2)*(-2*pow(x1,6) - 6*pow(x1,5)*x2 + 6*pow(x1,4) + 24*pow(x1,3)*x2 - 21*x1*x2 ))/dividend;

  return response;
}
double Optimization::hessianX2X1(double x1, double x2){
  double response;

  double dividend;

  dividend = pow(((x1,6) - 6*pow(x1,4) + 12*pow(x1,2) + 6*x1*x2 + 6*pow(x2,2)),(3/2));

  response = (sqrt(3/2)*(-2*pow(x1,6) - 6*pow(x1,5)*x2 + 6*pow(x1,4) + 24*pow(x1,3)*x2 - 21*x1*x2 ))/dividend;

  return response;
}
double Optimization::hessianX2X2(double x1, double x2){
  double response;

  double dividend;

  dividend = pow(((x1,6) - 6*pow(x1,4) + 12*pow(x1,2) + 6*x1*x2 + 6*pow(x2,2)),(3/2));

  response = (sqrt(3/2)*(pow(x1,2)*(2*pow(x1,4) - 12*pow(x1,2) + 21)))/dividend;

  return response;
}

double Optimization::goldenSection(double x1, double x2, double eps, double p, double d1, double d2)
{
  double theta1 = (3 - sqrt(5))/2.0;
  double theta2 = 1 - theta1;

  // Obtenção do intervalo [a, b]
  double a = 0.0;
  double s = p;
  double b = 2*p;

  // cout << "func " << derivativeX1(x1,x2)*d1 + derivativeX2(x1, x2)*d2 << endl;


  while (function(x1+b*d1, x2+b*d2) < function(x1+s*d1, x2+s*d2)){
    a = s;
    s = b;
    b = 2*b;
    // cout << "f: " << function(x1, x2) << endl;
    // cout << "phi1: " << phi(x1, x2, d1, d2, b) << endl;
    // cout << "phi2: " << phi(x1, x2, d1, d2, s) << endl;
    // cout << "a:" << a << " s:" << s << " b:" << b << endl;
  }

  // Obtenção de t*
  double u = a + theta1*(b - a);
  double v = a + theta2*(b - a);
  int cont = 0;
  while (abs(b - a) > eps){
    // cout << "b - a: " << b-a << endl;
    // if (cont > 1000) break;
    if (function(x1+u*d1, x2+u*d2) < function(x1+v*d1, x2+v*d2)){
      b = v;
      v = u;
      u = a + theta1*(b - a);
    } else {
      a = u;
      u = v;
      v = a + theta2*(b - a);
    }
    cont++;
  }

  // cout << endl << "t: " << (u+v)/2 << endl;
  return (u + v)/2;
}

double Optimization::armijo(double x1, double x2, double d1, double d2, double gama, double eta)
{
  double t = 1;
  double functionStep = function(x1 + t*d1, x2 + t*d2);
  double functionImage = function(x1, x2);
  double regularizer = eta*t*(derivativeX1(x1, x2)*d1 + derivativeX2(x1, x2)*d2);

  // cout << endl << "func step: " << functionStep << endl << "func image: " << functionImage << endl << "reg: " << regularizer << endl; 

  while (functionStep > functionImage + regularizer){
    t = gama*t;
    functionStep = function(x1 + t*d1, x2 + t*d2);
    regularizer = eta*t*(derivativeX1(x1, x2)*d1 + derivativeX2(x1, x2)*d2);

    // cout << endl << "regularizer: " << regularizer << endl;

  }

  return t;
}

double Optimization::gradient(double x1, double x2){
  int k = 0;
  double d1, d2, t;
  double gradX1 = derivativeX1(x1, x2);
  double gradX2 = derivativeX2(x1, x2);

  double bestX1, bestX2;

  bestX1 = x1;
  bestX2 = x2;
  int it = 0;
  double tolerance = 0.001;

  while (gradX1 != 0 && gradX2 != 0){
    if (k == 1000) break;
    d1 = -gradX1;
    d2 = -gradX2;
    t = 1;
    x1 = x1 + t*d1;
    x2 = x2 + t*d2;
    k++;

    if(function(bestX1,bestX2) > function(x1,x2) ){
      bestX1 = x1;
      bestX2 = x2;
      it = k;
    }

    gradX1 = derivativeX1(x1, x2);
    gradX2 = derivativeX2(x1, x2);

    // cout << endl << "Grad x1: " << gradX1 << endl << "Grad x2: " << gradX2 << endl;
  }
}

vector<double> Optimization::dnewton(double x1, double x2,double gradX1, double gradX2){
  double x1x2 = hessianX1X2(x1,x2);
  double x1x1 = hessianX1X1(x1,x2);
  double x2x2 = hessianX2X2(x1,x2);
  double x2x1 = hessianX2X1(x1,x2);

  double hessian[2][2] = {{x1x1,x1x2},{x2x1,x2x2}};

  double detHessian = hessian[0][0]*hessian[1][1] - (hessian[1][0]*hessian[0][1]);

  double inv_hessian[2][2] = {{hessian[1][1]/detHessian,-hessian[0][1]/detHessian},{-hessian[1][0]/detHessian, hessian[0][0]/detHessian}};

  double gradient[2] = {derivativeX1(x1,x2) ,derivativeX2(x1,x2)};

  vector<double> d = {-(inv_hessian[0][0]*gradX1 + inv_hessian[1][0]*gradX2),-(inv_hessian[1][0]*gradX1 + inv_hessian[1][1]*gradX2)};
  
  return d;
}

double Optimization::newton(double x1, double x2){
  int k = 0;
  double d1, d2, t;
  double gradX1 = derivativeX1(x1, x2);
  double gradX2 = derivativeX2(x1, x2);

  double bestX1, bestX2;

  bestX1 = x1;
  bestX2 = x2;
  int it = 0;
  vector<double> d;

  while (gradX1 != 0 && gradX2 != 0){
    if (k == 100) break;
    t = 1;
    d = dnewton(x1,x2,gradX1,gradX2);
    cout << "x1: " << x1 << "   x2: " << x2 << endl;
    cout << "D1: " << d[0] << "   D2: " << d[1] << endl;
    cout << "gradx1: " << gradX1 << endl;
    cout << "gradx2: " << gradX2 << endl << endl << endl;
    x1 = x1 + t*d[0];
    x2 = x2 + t*d[1];
    k++;
    if(function(bestX1,bestX2) > function(x1,x2) ){
      bestX1 = x1;
      bestX2 = x2;
      it = k;
    }
    gradX1 = derivativeX1(x1, x2);
    gradX2 = derivativeX2(x1, x2);
  }
    cout << "Best x: " << bestX1 << endl << "Best y: " << bestX2 << endl << "Interacao: " << it << endl;
};

double Optimization::quasiNewton(double x1, double x2){
  int k = 0;
  double d1, d2, t, old_x1, old_x2;
  vector<vector<double>> H = { { 1, 0}, { 0, 1}};
  double gradX1 = derivativeX1(x1, x2);
  double gradX2 = derivativeX2(x1, x2);

  double bestX1, bestX2;

  bestX1 = x1;
  bestX2 = x2;
  int it = 0;

  while (gradX1 != 0 && gradX2 != 0){
    if (k == 10) break;
    t = 1;

    d1 = - (H[0][0]*gradX1 + H[0][1]*gradX2);
    d2 = - (H[1][0]*gradX1 + H[1][1]*gradX2);
    
    cout << "x1: " << x1 << "   x2: " << x2 << endl;
    cout << "D1: " << d1 << "   D2: " << d2 << endl;
    cout << "gradx1: " << gradX1 << "   gradx2: " << gradX2 << endl << endl << endl;

    old_x1 = x1;
    old_x2 = x2;

    x1 = x1 + t*d1;
    x2 = x2 + t*d2;
    k++;
    if(function(bestX1,bestX2) > function(x1,x2) ){

      bestX1 = x1;
      bestX2 = x2;
      
      it = k;

    }
    gradX1 = derivativeX1(x1, x2);
    gradX2 = derivativeX2(x1, x2);
    H = BFGS(H,old_x1,old_x2,x1,x2);
  }
  cout << "Best x: " << bestX1 << endl << "Best y: " << bestX2 << endl << "Interacao: " << it << endl;
}

vector<double> Optimization::p(double x1,double x2,double x1_1, double x2_1){
  vector<double> p = {x1_1 - x1, x2_1 - x2};
  return p;
}

vector<double> Optimization::q(double x1,double x2,double x1_1, double x2_1){
  vector<double> q = { derivativeX1(x1_1,x2_1) - derivativeX1(x1,x2),derivativeX2(x1_1,x2_1) - derivativeX2(x1,x2) };
  return q;
}

vector<vector<double>> Optimization::BFGS(vector<vector<double>> H, double x1,double x2,double x1_1, double x2_1){
  // Hk+1 = firstpart + secondpart + thirdpart
  // firstpart = Hk

  vector<vector<double>> secondpart(2,vector<double>(2));
  vector<vector<double>> thirdpart(2,vector<double>(2));
  vector<vector<double>> result(2,vector<double>(2));
  
  vector<double> pk = p(x1,x2,x1_1,x2_1);
  vector<double> qk = q(x1,x2,x1_1,x2_1);

  double dividend = pk[0]*qk[0] + pk[1]*qk[1];


  double parantheses = 1 + (H[0][0]*qk[0]*qk[0] + H[1][0]*qk[0]*qk[1] + H[0][1]*qk[0]*qk[1] + H[1][1]*qk[1]*qk[1])/dividend;


  secondpart[0][0] = (parantheses/dividend)*(pk[0]*pk[0]);
  secondpart[0][1] = (parantheses/dividend)*(pk[1]*pk[0]);
  secondpart[1][0] = (parantheses/dividend)*(pk[1]*pk[0]);
  secondpart[1][1] = (parantheses/dividend)*(pk[1]*pk[1]);


  thirdpart[0][0] = (pk[0]*qk[0]*H[0][0] + (H[0][0]*qk[0] + H[0][1]*qk[1])*pk[0])/dividend;
  thirdpart[0][1] = (pk[0]*qk[1]*H[0][1] + (H[0][0]*qk[0] + H[0][1]*qk[1])*pk[1])/dividend;
  thirdpart[1][0] = (pk[1]*qk[0]*H[1][0] + (H[1][0]*qk[0] + H[1][1]*qk[1])*pk[0])/dividend;
  thirdpart[1][1] = (pk[1]*qk[1]*H[1][1] + (H[1][0]*qk[0] + H[1][1]*qk[1])*pk[1])/dividend;

  
  result[0][0] = H[0][0] + secondpart[0][0] - thirdpart[0][0];
  result[0][1] = H[0][1] + secondpart[0][1] - thirdpart[0][1];
  result[1][0] = H[1][0] + secondpart[1][0] - thirdpart[1][0];
  result[1][1] = H[1][1] + secondpart[1][1] - thirdpart[1][1];

  
  return result;
}


