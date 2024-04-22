#include <stdio.h>
#include <iostream>
#include "fileread.h"

#define UssFov 120
int main()
{
    TPointf64 e1focus0 = {-3,0};
    TPointf64 e1focus1 = {3,0};
    float64 e1semiaxial = 5;
    Tflag e1type = ELLIPSE;

    TPointf64 e2focus0 = {2,4};
    TPointf64 e2focus1 = {3,8};
    float64 e2semiaxial = 4;
    Tflag e2type = CIRCLE;

    stIntersection solution =  getellipseintersection(e1focus0,e1focus1,e1semiaxial,UssFov,e1type,
                                                   e2focus0,e2focus1,e2semiaxial,UssFov,e2type);

     return 0;

}
