#ifndef ALGOHEADER_H
#define ALGOHEADER_H

#include <assert.h>
#include <stdio.h>

//#undef
#define NDEBUG
typedef unsigned char  Tflag;
typedef float float32;
typedef  float32 float32_t;
typedef  float32_t  matrixf32;
typedef  double float64;

typedef struct
{
    float64 x;
    float64 y;
}TPointf64;

typedef TPointf64 vec2Df64;
typedef float64 matrixf64;

#ifndef PI
#define PI 3.14159265358979323846
#endif

#ifndef eps
#define eps  1e-5
#endif

#ifndef eps1
#define eps1 1e-10
#endif

#ifndef eps2
#define eps2 1e-12
#endif

#ifndef sgn
#define sgn(x)   ((x>0)?1:-1)
#endif

#ifndef IsFloatEqual0
#define IsFloatEqual0(x) (((x> -eps) && (x< eps))?1:0)
#endif

#ifndef MaxC
#define MaxC(a,b) (((a) > (b)) ? (a) : (b))
#endif

#ifndef MinC
#define MinC(a,b) (((a) < (b)) ? (a) : (b))
#endif

#ifndef Deg2Rad
#define Deg2Rad(x) (x*PI / 180.0f)
#endif
#endif // ALGOHEADER_H
