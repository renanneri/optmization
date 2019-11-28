#include <stdio.h>
#include <iostream>
#include <math.h> 
#include "otm.h"

int main()
{
    Optimization opt = Optimization();
    opt.gradient(0.1,0.1);
    return 0;
}