#include <stdio.h>
#include <iostream>
#include <math.h> 
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

  // cout << "dx1 " << response << endl;

  return response;
}

double Optimization::derivativeX2(double x1, double x2){

  double response;

  double dividend;

  dividend = pow(x1,6)/6 - pow(x1,4) + 2*pow(x1,2) + x1*x2 + pow(x2,2);

  response = (x1 + 2*x2)/(2*sqrt(dividend));

  // cout << "dx2 " << response << endl;

  return response;
}

double Optimization::goldenSection(double x1, double x2, double eps, double p, double d1, double d2)
{
  double theta1 = (3 - sqrt(5))/2;
  double theta2 = 1 - theta1;

  // Obtenção do intervalo [a, b]
  double a = 0.0;
  double s = p;
  double b = 2*p;

  // cout << "func " << function(0.32, 1.28) << endl;


  while (phi(x1, x2, d1, d2, b) < phi(x1, x2, d1, d2, s)){
    a = s;
    s = b;
    b = 2*b;
    cout << "f: " << function(x1, x2) << endl;
    cout << "phi1: " << phi(x1, x2, d1, d2, b) << endl;
    cout << "phi2: " << phi(x1, x2, d1, d2, s) << endl;
    // cout << "a:" << a << " s:" << s << " b:" << b << endl;
  }

  // Obtenção de t*
  double u = a + theta1*(b - a);
  double v = a + theta2*(b - a);
  int cont = 0;
  while ((b - a) > eps){
    if (cont > 1000) break;
    if (phi(x1, x2, d1, d2, u) < phi(x1, x2, d1, d2, v)){
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

  return (u + v)/2;
}

double Optimization::gradient(double x1, double x2){
  int k = 0;

  double d1, d2, t;

  double gradX1 = derivativeX1(x1, x2);
  double gradX2 = derivativeX2(x1, x2);


  while (gradX1 != 0 && gradX2 != 0){
    // if (k == 10000) break;
    // cout << "gradx1: " << gradX1 << endl;
    // cout << "gradx2: " << gradX2 << endl;
    d1 = -gradX1;
    d2 = -gradX2;

    t = 1;// goldenSection(x1, x2, 0.000001, 0.000001, d1, d2);
    x1 = x1 + t*d1;
    x2 = x2 + t*d2;
    k++;
    gradX1 = derivativeX1(x1, x2);
    gradX2 = derivativeX2(x1, x2);
    // cout << k << endl;
  }
  cout << "x1 " << x1 << endl;
  cout << "x2 " << x2 << endl;
}


