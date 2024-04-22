#include "ellipse.h"
stEllipse EllipseGenerating(const float64 centerx ,
                       const float64 centery ,
                       const float64 a_,
                       const float64 b_,
                       const float64 theta_)
{

    assert(!IsFloatEqual0(a_) && !IsFloatEqual0(b_));

    float64 A = (sin(theta_)/b_)*(sin(theta_)/b_) + (cos(theta_)/a_)*(cos(theta_)/a_);
    float64 B = 2*((1/a_)*(1/a_) - (1/b_)*(1/b_))*sin(theta_)*cos(theta_);
    float64 C = (cos(theta_)/b_)*(cos(theta_)/b_) + (sin(theta_)/a_)*(sin(theta_)/a_);
    float64 D = -(2*A*centerx + B*centery);
    float64 E = -(2*C*centery + B*centerx);
    float64 F = -(D*centerx + E*centery)/2 - 1;

    stEllipse ellipse ={
        .center = {
            .x = centerx,
            .y = centery,
        },

        .a = a_,
        .b = b_,
        .theta = theta_,
        .homogeneous = {
            .A= A,
            .B= B,
            .C= C,
            .D= D,
            .E= E,
            .F= F,},
        .M = {{A, B/2, D/2},
              {B/2, C, E/2},
              {D/2, E/2, F}},
    };

    ellipse.Vecfov.x = cos(ellipse.theta);
    ellipse.Vecfov.y = sin(ellipse.theta);

    Vec_rotate(&ellipse.Vecfov,PI/2);

    return ellipse;

}
float64  LinearCombinationSolve(const stEllipse *el1,const stEllipse *el2)
{
    /*计算 x^3的系数*/
    float64 a = Matrix_det((matrixf64 *)el1->M,3,3);

    /*计算 x^2的系数*/
    matrixf64 M1minor[3][3];
    memset(&M1minor[0][0],0,sizeof (matrixf64)*3*3);
    Matrix_cofactor((matrixf64 *)el1->M,3,3,(matrixf64 *)M1minor);

    matrixf64 M1minor_mult_M2[3][3];
    memset(&M1minor_mult_M2[0][0],0,sizeof (matrixf64)*3*3);
    Matrix_left_mul((matrixf64 *)M1minor,(matrixf64 *)el2->M,(matrixf64 *)M1minor_mult_M2,3,3,3,3);
    float64 b = Matrix_trace((matrixf64 *)M1minor_mult_M2,3,3);

    /*计算x的系数*/
    matrixf64 M2minor[3][3];
    memset(&M2minor[0][0],0,sizeof (matrixf64)*3*3);
    Matrix_cofactor((matrixf64 *)el2->M,3,3,(matrixf64 *)M2minor);

    matrixf64 M1_mult_M2minor[3][3];
    memset(&M1_mult_M2minor[0][0],0,sizeof (matrixf64)*3*3);
    Matrix_left_mul((matrixf64 *)el1->M,(matrixf64 *)M2minor,(matrixf64 *)M1_mult_M2minor,3,3,3,3);
    float64 c = Matrix_trace((matrixf64 *)M1_mult_M2minor,3,3);

    /*计算常数项*/
    float64 d = Matrix_det((matrixf64 *)el2->M,3,3);

    /*解一元三次方程*/
    stEquSolution solve = CubicEquationSolving(a,b,c,d);

      return solve.sol[0];


}

void  LinearCombinationMatrix(const stEllipse *el1,const stEllipse *el2,const float64 coeff,matrixf64 * s)
{
    Matrix_assign((matrixf64 *)s,(matrixf64 *)el1->M,3,3);

    Matrix_dot_mul((matrixf64 *)s,coeff,3,3);

    Matrix_adder((matrixf64 *)s,(matrixf64 *)el2->M,3,3);

}

Tflag  CheckIsPoint(const stEllipse *ellipse, const TPointf64 apoint)
{
          return IsFloatEqual0( ellipse->homogeneous.A * apoint.x * apoint.x
            + ellipse->homogeneous.B * apoint.x * apoint.y
            + ellipse->homogeneous.C * apoint.y * apoint.y
            + ellipse->homogeneous.D * apoint.x
            + ellipse->homogeneous.E * apoint.y
            + ellipse->homogeneous.F);
}

Tflag    CheckIsFov(const TPointf64 secpoint, const stEllipse *ellipse, const float64 fov)
{
    vec2Df64 secvec ;
    secvec.x = (secpoint.x - ellipse->center.x);
    secvec.y = (secpoint.y - ellipse->center.y);

    return (fabs(Vec_angle(secvec,ellipse->Vecfov))< fov);
}

int  Line_Ellipse_intersection(const stEllipse *ellipse, const float64 a,const float64 b, const float64 c, TPointf64* point)
{
    if(a == 0 && b == 0)
    {
        return 0;
    }

    if(fabs(a) > fabs(b))
    {

        float64 d = ellipse->M[0][0]*b*b - 2*ellipse->M[0][1]*a*b + ellipse->M[1][1]*a*a;
        float64 e = 2*(ellipse->M[1][2]*a*a + ellipse->M[0][0]*b*c - ellipse->M[0][1]*a*c - ellipse->M[0][2]*a*b);
        float64 f = ellipse->M[0][0]*c*c - 2*ellipse->M[0][2]*a*c + ellipse->M[2][2]*a*a;
        stEquSolution sol = QuadraticEquationSolving(d,e,f);

        for(int i =0;i<sol.realnum;i++)
        {
            point[i].y = sol.sol[i];
            point[i].x = (-c - b* point[i].y)/a;
        }

        return sol.realnum;

    }

    float64 d = ellipse->M[0][0]*b*b - 2*ellipse->M[0][1]*a*b + ellipse->M[1][1]*a*a;
    float64 e = 2*(ellipse->M[0][2]*b*b + ellipse->M[1][1]*a*c - ellipse->M[0][1]*b*c - ellipse->M[1][2]*a*b);
    float64 f = ellipse->M[1][1]*c*c - 2*ellipse->M[1][2]*b*c + ellipse->M[2][2]*b*b;

    stEquSolution sol = QuadraticEquationSolving(d,e,f);
    for(int i =0;i<sol.realnum;i++)
    {
        point[i].x = sol.sol[i];
        point[i].y = (-c - a* point[i].x)/b;
    }

    return sol.realnum;

}

stIntersection ellipse_ellipse_intersection(const stEllipse *ell1,const stEllipse *ell2, const matrixf64 Mat[3][3])
{
    float64 a = Mat[0][0]; float64 b =Mat[0][1] ; float64 c =Mat[1][1] ; float64 d = Mat[0][2] ; float64 e = Mat[1][2]; float64 f = Mat[2][2];
    float64 a33 = Mat[0][0] * Mat[1][1] - Mat[0][1]* Mat[1][0];

    stIntersection points;
    memset(&points,0,sizeof (stIntersection));
    //need check
    if(a33 > eps2)//    2 1 90 -1 -1//    2 1 180 3 -3
    {
        points.point[0].x = (b*e - c*d) / a33;
        points.point[0].y = (b*d - a*e) / a33;

        points.realnum = (CheckIsPoint(ell1,points.point[0]) && CheckIsPoint(ell2,points.point[0]))?
                          1:0;
        return points;
    }

    if(a33 < -eps2)
    {
        a33 = sqrt(-a33);
        c = (b*d - a*e) / a33;
        TPointf64 tmp[4];
        memset(tmp,0,sizeof (TPointf64));

        int num1 = Line_Ellipse_intersection(ell2,  a, b-a33, d-c, &tmp[0]);
        int num2 = Line_Ellipse_intersection(ell2,  a, b+a33, d+c, &tmp[num1]);

        for(int i =0;i< (num1 + num2);i++)
        {
            if(CheckIsPoint(ell1,tmp[i]))
            {
                points.point[points.realnum] = tmp[i];
                points.realnum ++ ;

            }
        }

        return points;
    }

    if(fabs(b) < eps2)
    {
        if(fabs(a) > MaxC(fabs(c),eps2))
        {
            stEquSolution sol =  QuadraticEquationSolving(a,2*d,f);
            int num1 = 0,num2=0;
            TPointf64 tmp[4];
            memset(tmp,0,sizeof (TPointf64));
            if(sol.realnum > 0)
            {
               num1 = Line_Ellipse_intersection(ell2,  -1,0,sol.sol[0],&tmp[0]);
            }

            if(sol.realnum > 1)
            {
                num2 = Line_Ellipse_intersection(ell2,  -1,0,sol.sol[1],&tmp[num1]);
            }

            for(int i =0;i< (num1 + num2);i++)
            {
                if(CheckIsPoint(ell1,tmp[i]))
                {
                    points.point[points.realnum] = tmp[i];
                    points.realnum ++ ;

                }
            }
            return points;
        }

        if(fabs(c) > MaxC(fabs(a), eps2))//    [4.0, 1.0, 270.0, 0.0, 1.0] //    [3.0, 1.0, 360.0, 0.0, -1.0]
        {
            stEquSolution sol =  QuadraticEquationSolving(c, 2*e, f);
            int num1 = 0,num2=0;

            TPointf64 tmp[4];
            memset(tmp,0,sizeof (TPointf64));

            if(sol.realnum > 0)
            {
               num1 = Line_Ellipse_intersection(ell2,  0,-1,sol.sol[0],&tmp[0]);
            }

            if(sol.realnum > 1)
            {
                num2 = Line_Ellipse_intersection(ell2,  0,-1,sol.sol[1],&tmp[num1]);
            }

            for(int i =0;i< (num1 + num2);i++)
            {
                if(CheckIsPoint(ell1,tmp[i]))
                {
                    points.point[points.realnum] = tmp[i];
                    points.realnum ++ ;

                }
            }
            return points;;
        }

        TPointf64 tmp[4];
        memset(tmp,0,sizeof (TPointf64));
        int num = Line_Ellipse_intersection(ell2,2*d,2*e,f,&tmp[0]);

        for(int i =0;i< num;i++)
        {
            if(CheckIsPoint(ell1,tmp[i]))
            {
                points.point[points.realnum] = tmp[i];
                points.realnum ++ ;

            }
        }
        return points;
    }

    if(a<0)
    {
        a = -a; b = -b; c = -c; d = -d; e = -e; f = -f;
    }

    float64 s = d*d - a*f;

    if (b*d*e < -eps2 || fabs(a*e*e - c*d*d) > eps2 || s < -eps2)
    {
        return points;
    }

    if(s< eps)
    {
        TPointf64 tmp[4];
        memset(tmp,0,sizeof (TPointf64));
        int num = Line_Ellipse_intersection(ell2,a,b,d,&tmp[0]);
        for(int i = 0;i< num;i++)
        {
            if(CheckIsPoint(ell1,tmp[i]))
            {
                points.point[points.realnum] = tmp[i];
                points.realnum ++ ;

            }
        }
        return points;
    }

    s = sqrt(s);

    int num1 = 0,num2=0;
    TPointf64 tmp[4];
    memset(tmp,0,sizeof (TPointf64));
    num1 = Line_Ellipse_intersection(ell2,  a,b,d-s,&tmp[0]);
    num2 = Line_Ellipse_intersection(ell2,  a,b,d+s,&tmp[num1]);
    for(int i =0;i< (num1 + num2);i++)
    {
        if(CheckIsPoint(ell1,tmp[i]))
        {
            points.point[points.realnum] = tmp[i];
            points.realnum ++ ;

        }
    }
    return points;
}

Tflag ExpressionConversion(float64 *property, const float64 x0,const float64 y0 , const float64 x1, const float64 y1, const float64 semiaxial, const Tflag type)
{
    if(type)
    {
        float64 centerx = (float64) (x1 + x0)/2.0f;
        float64 centery = (float64) (y1 + y0)/2.0f;

        float64 a = semiaxial;
        float64 c = sqrt(pow((x1-x0),2) + pow((y1-y0),2))/2;
        float64 b = sqrt(a*a - c*c);

        if(isnan(b))
        {
            return 0;
        }

        float64 theta = (fabs(x1-x0) < eps1)?(PI/2):atan((y1-y0)/(x1-x0));

        property[0] = centerx;
        property[1] = centery;
        property[2] = a;
        property[3] = b;
        property[4] = theta;

        return 1;

    }
    else //circle
    {
        float64 centerx = (float64) x0;
        float64 centery = (float64) y0;

        float64 a = semiaxial;
        float64 b = a;
        float64 theta = (fabs(x1-x0) < eps1)?(PI/2):atan((y1-y0)/(x1-x0));

        property[0] = centerx;
        property[1] = centery;
        property[2] = a;
        property[3] = b;
        property[4] = theta;

        return 1;
    }
}
stIntersection ellipse_intersection(const stEllipse *ellipse1,const stEllipse *ellipse2)
{
    float64 coefficent = LinearCombinationSolve(ellipse1,ellipse2);

    matrixf64 S[3][3];
    memset((matrixf64*)S,0,sizeof (matrixf64)*3*3);
    LinearCombinationMatrix(ellipse1,ellipse2,coefficent,(matrixf64*)S);

    return ellipse_ellipse_intersection(ellipse1,ellipse2,S);
}
stIntersection getellipseintersection(const TPointf64 e1focus0,const TPointf64 e1focus1 ,const float64 e1semiaxial,const float64 ussFov1,const Tflag e1type,
                                           const TPointf64 e2focus0 , const TPointf64 e2focus1 ,const float64 e2semiaxial,const float64 ussFov2,const Tflag e2type)
{
    float64 property1[5];
    float64 property2[5];
    memset(property1,0,sizeof (float64)*5);
    memset(property2,0,sizeof (float64)*5);
    ExpressionConversion(property1,e1focus0.x,e1focus0.y,e1focus1.x,e1focus1.y,e1semiaxial,e1type);
    ExpressionConversion(property2,e2focus0.x,e2focus0.y,e2focus1.x,e2focus1.y,e2semiaxial,e2type);

    stEllipse Ellipse1 = EllipseGenerating(  property1[0] ,   property1[1] , property1[2],  property1[3],  property1[4]);
    stEllipse Ellipse2 = EllipseGenerating(  property2[0] ,   property2[1] , property2[2],  property2[3],  property2[4]);

    stIntersection solution = ellipse_intersection(&Ellipse1,&Ellipse2);
    stIntersection solutionfov;
    memset(&solutionfov,0,sizeof (stIntersection));

    for(int i = 0; i< solution.realnum;i++)
    {
        if(CheckIsFov(solution.point[i],&Ellipse1,Deg2Rad(ussFov1/2)) &&
                CheckIsFov(solution.point[i],&Ellipse2,Deg2Rad(ussFov2/2)))
        {
            solutionfov.point[solutionfov.realnum] = solution.point[i];
            solutionfov.realnum ++;
        }
    }

    return  solutionfov;


}
