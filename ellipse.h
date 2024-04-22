#ifndef ELLIPSE_H
#define ELLIPSE_H
#include "matrix.h"
#include "equation.h"
#define SPATIALDIMENSIONS 3
#define POLYOMIAL         6

#define CIRCLE  0
#define ELLIPSE 1

typedef struct
{
    TPointf64 point[4];
    int realnum;

}stIntersection;

typedef  struct
{
    float64 X[4];
    float64 Y[4];
    int realnum;
}PyIntersection;

typedef struct
{
   float64 A;
   float64 B;
   float64 C;
   float64 D;
   float64 E;
   float64 F;
}stPC;

typedef struct
{
    TPointf64 center;
    float64  a;
    float64  b;
    float64  theta;
    vec2Df64 Vecfov;
    stPC  homogeneous;
    matrixf64  M[SPATIALDIMENSIONS][SPATIALDIMENSIONS];
}stEllipse;

/*public*/
#ifdef __cplusplus
extern "C" {
#endif

stIntersection getellipseintersection(const TPointf64 e1focus0,const TPointf64 e1focus1 ,const float64 e1semiaxial,const float64 ussFov1,const Tflag e1type,
                                           const TPointf64 e2focus0 , const TPointf64 e2focus1 ,const float64 e2semiaxial,const float64 ussFov2,const Tflag e2type);
#ifdef __cplusplus
}
#endif

/*private*/
stEllipse EllipseGenerating(const float64 centerx , const float64 centery ,const float64 a_,const float64 b_,const float64 theta_);
Tflag    CheckIsPoint(const stEllipse *ellipse, const TPointf64 apoint);
stIntersection ellipse_intersection(const stEllipse *ellipse1,const stEllipse *ellipse2);
stEllipse EllipseGenerating(const float64 centerx , const float64 centery ,const float64 a_,const float64 b_,const float64 theta_);
float64  LinearCombinationSolve(const stEllipse *el1,const stEllipse *el2);
void     LinearCombinationMatrix(const stEllipse *el1,const stEllipse *el2,const float64 coeff,matrixf64 * s);
Tflag    CheckIsPoint(const stEllipse *ellipse, const TPointf64 apoint);

Tflag ExpressionConversion(float64 *property, const float64 x0,const float64 y0 , const float64 x1, const float64 y1, const float64 a, const Tflag type);
int  Line_Ellipse_intersection(const stEllipse *ellipse, const float64 a,const float64 b, const float64 c, TPointf64* point);
stIntersection ellipse_ellipse_intersection(const stEllipse *ell1,const stEllipse *ell2, const matrixf64 Mat[3][3]);
Tflag    CheckIsFov(const TPointf64 secpoint, const stEllipse *ellipse, const float64 fov);
#endif // ELLIPSE_H
