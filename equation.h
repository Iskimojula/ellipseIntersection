#ifndef EQUATION_H
#define EQUATION_H
#include "matrix.h"


#define QUARTIC    4
#define CUBIC      3
#define QUADRATIC  2
#define SIMPLE     1
typedef  struct
{
    float64 sol[QUARTIC];
    int maxsol;
    int realnum;

}stEquSolution;

stEquSolution QuarticEquationSolving(const float64 a,
                                     const float64 b,
                                     const float64 c,
                                     const float64 d,
                                     const float64 e);//天珩公式

stEquSolution CubicEquationSolving(const float64 a,//盛金公式
                                   const float64 b,
                                   const float64 c,
                                   const float64 d);

//存在两个相同根时,输出的个数为1
stEquSolution QuadraticEquationSolving(const float64 a,
                                       const float64 b,
                                       const float64 c);

stEquSolution SimpleEquationSolving(const float64 a, const float64 b);

#endif // EQUATION_H
