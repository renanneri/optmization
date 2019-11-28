#include <stdio.h>
#include <iostream>
#include <math.h> 
#include "otm.h"

using namespace std;

int main()
{
    Optimization opt = Optimization();
    // opt.derivativeX1(1,1);
    // opt.derivativeX2(1,1);
    cout << "Newton: " << endl;
    opt.newton(1, 1);
    cout << "Quase Newton: " << endl;
    return 0;
}