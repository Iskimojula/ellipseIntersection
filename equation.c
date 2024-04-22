#include "equation.h"

stEquSolution QuarticEquationSolving(const float64 a,
                                     const float64 b,
                                     const float64 c,
                                     const float64 d,
                                     const float64 e)
{
     stEquSolution sol;
     memset(&sol,0,sizeof (stEquSolution));
     sol.maxsol = QUARTIC;
     float64 D = 3*pow(b,2) - 8*a*c; //3b^2 - 8ac
     float64 E = -pow(b,3) + 4*a*b*c-8*pow(a,2)*d;// -b^3 + 4*a*b*Sc - 8*a^2*d
     float64 F= 3*pow(b,4)+16*pow(a,2)*pow(c,2)-16*a*pow(b,2)*c + 16*pow(a,2)*b*d - 64*pow(a,3)*e;

     float64 A    = pow(D,2)-3*F;
     float64 B    = D*F-9*pow(E,2);
     float64 C    = pow(F,2) - 3*D*pow(E,2);
     float64 DETA = pow(B,2) - 4*A*C;

     if(IsFloatEqual0(D) && IsFloatEqual0(E) && IsFloatEqual0(F))
     {
        for(int i = 0;i< QUARTIC;i++)
        {
            sol.sol[i] = -b/(4*a);
        }
        sol.realnum = QUARTIC;
        return sol;
     }

     if(!IsFloatEqual0(D*E*F) && IsFloatEqual0(A) && IsFloatEqual0(B) && IsFloatEqual0(C))
     {
         sol.sol[0] = (-b*D + 9*E)/(4*a*D);
         sol.sol[1] = (-b*D - 3*E)/(4*a*D);
         sol.sol[2] = (-b*D - 3*E)/(4*a*D);
         sol.sol[3] = (-b*D - 3*E)/(4*a*D);
         sol.realnum = QUARTIC;
         return sol;
     }

     if(IsFloatEqual0(E) && IsFloatEqual0(F) && D>0)
     {
         sol.sol[0] = (-b +sqrt(D))/(4*a);
         sol.sol[1] = (-b +sqrt(D))/(4*a);

         sol.sol[2] = (-b -sqrt(D))/(4*a);
         sol.sol[3] = (-b -sqrt(D))/(4*a);
         sol.realnum = QUARTIC;
         return sol;
     }

     if(!IsFloatEqual0(A*B*C) && IsFloatEqual0(DETA) && A*B>0)
     {
        float64 tmpvar1 = 2*A*E/B;
        float64 tmpvar2 = sqrt(2*B/A);

        sol.sol[0] = (-b +tmpvar1 + tmpvar2)/(4*a);
        sol.sol[1] = (-b +tmpvar1 - tmpvar2)/(4*a);

        sol.sol[2] = (-b -tmpvar1)/(4*a);
        sol.sol[3] = (-b -tmpvar1)/(4*a);
        sol.realnum = QUARTIC;
        return sol;
     }

     if(!IsFloatEqual0(A*B*C) && IsFloatEqual0(DETA) && A*B<0)
     {
         float64 tmpvar1 = 2*A*E/B;
         sol.sol[0] = (-b -tmpvar1)/(4*a);
         sol.sol[1] = (-b -tmpvar1)/(4*a);
         sol.realnum = QUARTIC - 2;
         return sol;
     }

     if(DETA > 0)
     {
         float64 z1 = A*D+3*(-B + sqrt(DETA))/2;
         float64 z2 = A*D+3*(-B - sqrt(DETA))/2;
         float64 z = pow(D,2)-D*(cbrt(z1)+cbrt(z2)) + pow((cbrt(z1)+cbrt(z2)),2)-3*A;

         float64 tmp1 = sgn(E)*sqrt((D+cbrt(z1)+cbrt(z2))/3);
         float64 tmp2 = sqrt((2*D-(cbrt(z1)+cbrt(z2))+2*sqrt(z))/3);

         sol.sol[0] = (-b + tmp1 + tmp2)/(4*a);
         sol.sol[1] = (-b + tmp1 - tmp2)/(4*a);
         sol.realnum = QUARTIC - 2;
         return sol;
     }

     if(DETA < 0 && D>0 && F>0)
     {
         float64 theta120 = (2*PI)/3;
         float64 theta = acos((3*B-2*A*D)/(2*A*sqrt(A)));
         float64 y1 = (D-2*sqrt(A)*cos(theta/3))/3;
         float64 y2 = (D-2*sqrt(A)*cos(theta/3 + theta120))/3;
         float64 y3 = (D-2*sqrt(A)*cos(theta/3 - theta120))/3;

         if(IsFloatEqual0(E))
         {
             float64 tmp1 = sqrt(D+2*sqrt(F));
             float64 tmp2 = sqrt(D-2*sqrt(F));

             assert(!isnan(tmp2));
             sol.sol[0] = (-b +tmp1 )/(4*a);
             sol.sol[1] = (-b -tmp1 )/(4*a);

             sol.sol[2] = (-b +tmp2)/(4*a);
             sol.sol[3] = (-b -tmp2)/(4*a);
             sol.realnum = QUARTIC;
             return sol;
         }
         else
         {
             sol.sol[0] = (-b +sgn(E)*sqrt(y1) + (sqrt(y2)+sqrt(y3)))/(4*a);
             sol.sol[1] = (-b +sgn(E)*sqrt(y1) - (sqrt(y2)+sqrt(y3)))/(4*a);

             sol.sol[2] = (-b -sgn(E)*sqrt(y1) + (sqrt(y2)-sqrt(y3)))/(4*a);
             sol.sol[3] = (-b -sgn(E)*sqrt(y1) + (sqrt(y2)-sqrt(y3)))/(4*a);
             sol.realnum = QUARTIC;
             return sol;
         }

     }

              return sol;
}
stEquSolution CubicEquationSolving(const float64 a,
                                   const float64 b,
                                   const float64 c,
                                   const float64 d)
{
      stEquSolution sol;
      memset(&sol,0,sizeof (stEquSolution));
      sol.maxsol = CUBIC;

      float64 A = pow(b,2)-3*a*c;
      float64 B = b*c - 9*a*d;
      float64 C = pow(c,2)-3*b*d;
      float64 DETA = pow(B,2)-4*A*C;

      if( A == 0 && B == 0)
      {
          sol.sol[0] = -b/(3*a);
          sol.sol[1] = -b/(3*a);
          sol.sol[2] = -b/(3*a);

          sol.realnum = CUBIC -2;
          return sol;
      }

      if(DETA > 0 )
      {
          float64 Y1 = A*b + (3*a/2)*(-B + sqrt(DETA));
          float64 Y2 = A*b + (3*a/2)*(-B - sqrt(DETA));

          sol.sol[0] = - (b+cbrt(Y1)+cbrt(Y2))/(3*a);
          sol.realnum = CUBIC -2;
          return sol;
      }

      if(DETA > -eps2)
      {
          sol.sol[0] = -B/(2*A);
          sol.sol[1] = -B/(2*A);
          sol.sol[2] =  (a*B-A*b)/(a*A);
          sol.realnum = CUBIC;
          return sol;
      }

      if(DETA < 0)
      {
          float64 theta = acos((2*A*b - 3*a*B)/(2*A*sqrt(A)));
           for(int i =0;i<CUBIC;i++)
           {
                sol.sol[i] = -(b+2*sqrt(A)*cos((theta + 2*i*PI)/3))/(3*a);
           }
           sol.realnum = CUBIC;
           return sol;
      }
      return sol;

}
stEquSolution QuadraticEquationSolving(const float64 a,
                                       const float64 b,
                                       const float64 c)
{
    stEquSolution sol;
    memset(&sol,0,sizeof (stEquSolution));
    sol.maxsol = QUADRATIC;

    double A = pow(b,2) - 4*a*c;

    if(A<0)
    {
        if(A< -eps1)
        {
            return  sol;
        }
        A = 0;
    }


    if(sqrt(A) == 0 || sqrt(A) < a*eps1)
    {
        sol.sol[0] = -b/(2*a);
        sol.sol[1] = -b/(2*a);

        sol.realnum = QUADRATIC-1;
        return sol;
    }

    if(A > 0)
    {
        sol.sol[1] = (-b +sqrt(A))/(2*a);
        sol.sol[0] = (-b -sqrt(A))/(2*a);
        sol.realnum = QUADRATIC;
        return sol;
    }

    return sol;
}
stEquSolution SimpleEquationSolving(const float64 a, const float64 b)
{
    stEquSolution sol;
    memset(&sol,0,sizeof (stEquSolution));
    sol.maxsol = SIMPLE;

    if(!IsFloatEqual0(a))
    {
        sol.sol[0] = -b/a;
        sol.realnum = SIMPLE;
        return sol;
    }

    return sol;
}
