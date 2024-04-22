#include "matrix.h"

static matrixf64 MatDataPool[MAXMATRIXSIZE*MAXMATRIXSIZE] = {0};
void Matrix_left_mul(const matrixf64 *left,const matrixf64 *right,matrixf64 *rest,const int row1,const int col1,const int row2,const int col2)
{
    assert(col1==row2);
    for (int i = 0; i < row1; i++)//从第i行开始
    {
        for (int j = 0; j < col2; j++)//从第j列开始
        {
            for (int k = 0; k < col1; k++)//i行元素和j列元素相乘，结果累加
            {
               *(rest+i*col2+j) += (*(left+i*col1+k)) * (*(right+k*col2+j));
            }
        }
    }
}

void Matrix_transposition(const  matrixf64 *A,matrixf64 *A_,const int row,const int col)
{
    for (int i = 0; i < row; i++)//从第i行开始
    {
        for (int j = 0; j < col; j++)//从第j列开始
        {
               *(A_+j*row+i)  = *(A+i*col+j);
        }
    }
}

void Matrix_adder( matrixf64 *augend,const matrixf64 *addend,const int row,const int col)
{
    for (int i = 0; i < row; i++)//从第i行开始
    {
        for (int j = 0; j < col; j++)//从第j列开始
        {
               *(augend+i*col+j) = *(augend+i*col+j) + *(addend+i*col+j);

        }
    }
}


void Matrix_dot_mul( matrixf64 *multiplicand,const float64 multiplier,const int row,const int col)
{
    for (int i = 0; i < row; i++)//从第i行开始
    {
        for (int j = 0; j < col; j++)//从第j列开始
        {
               *(multiplicand+i*col+j) = (*(multiplicand+i*col+j)) * multiplier;

        }
    }
}
void Matrix_assign( matrixf64 *des,const matrixf64 *str,const int row,const int col)
{
    for (int i = 0; i < row; i++)//从第i行开始
    {
        for (int j = 0; j < col; j++)//从第j列开始
        {
               *(des+i*col+j) = *(str+i*col+j);

        }
    }
}
void Matrix_minor(const matrixf64 *A_,const int i_,const int j_,const int row,const int col,matrixf64 *res )
{

    matrixf64 * ptr = (matrixf64 * ) MatDataPool;
    memset(ptr,0,sizeof (matrixf64)*MAXMATRIXSIZE*MAXMATRIXSIZE);

    for(int i = 0;i<row;i++ )
    {
        if(i == i_)
        {
            continue;
        }

        for(int j = 0; j<col;j++)
        {
           if(j==j_)
           {
               continue;
           }
           *ptr = *(A_ +i*col+ j);
            ptr ++;
        }
    }

    int i = 0,j =0;
    for(int k = 0;k< (col-1)*(row -1);k++)
    {
        if(!k%(col -1) && k!=0)
        {
            i++;
            j=0;
        }

        *(res + i*col + j) = MatDataPool[k];
        j++;
    }

}

float64 Matrix_det(const matrixf64 *A,const int row, const int col)
{
      if(row == 2 && col == 2)
      {
          return (*(A + 0*col + 0) )*(*(A + 1*col + 1)) - (*(A + 0*col + 1) )*(*(A + 1*col + 0));
      }

      else
      {
          float64 det = 0;
          for(int j = 0; j< col ;j++)
          {
              matrixf64 newA[MAXMATRIXSIZE][MAXMATRIXSIZE] = {0};
              Matrix_minor((matrixf64*)A,0,j, row, col,(matrixf64 *)newA);
              det += pow(-1,(0+j))*(*(A + 0*col + j))*Matrix_det((matrixf64 *)newA, row -1, col -1 );
          }
        return det;
      }
}

float64 Matrix_trace(const matrixf64 *A,const int row,const int col)
{
    assert(row == col);

    float64 trace = 0;

    for(int i =0;i< col;i++)
    {
        trace += *(A + i*col + i);
    }
    return trace;
}

void Matrix_cofactor(const matrixf64 *A,const int row,const int col,matrixf64* res)
{
    for(int i =0;i<row;i++)
    {
        for(int j = 0;j<col;j++)
        {

            matrixf64 newA[MAXMATRIXSIZE][MAXMATRIXSIZE] = {0};
            Matrix_minor((matrixf64*)A,i,j, row, col,(matrixf64 *)newA);
            *(res + i*col + j) = pow(-1,(i+j))*Matrix_det( (matrixf64 *)newA, row-1, col-1);
        }
    }
}

void Matrix_inverse(const matrixf64 *A,const int row,const int col,matrixf64 *res)
{
    matrixf64 tmpArray[MAXMATRIXSIZE][MAXMATRIXSIZE] = {0};
    matrixf64 *tmp = (matrixf64 *) tmpArray;
    float64 det = Matrix_det((matrixf64 *)A,3,3);
    assert(!IsFloatEqual0(det));
    Matrix_cofactor(A,row,col,tmp);
    Matrix_transposition((matrixf64 *)tmp,(matrixf64 *)res,3,3);
    Matrix_dot_mul((matrixf64 *)res,(float64)1/det,3,3);
}

void Vec_rotate(vec2Df64 *vec,const float64 theta_)
{
    float64 x0 = vec->x;
    float64 y0 = vec->y;

    vec->x = x0 * cos(theta_) - y0 * sin(theta_);
    vec->y = x0 * sin(theta_) - y0 * cos(theta_);

}

float64 Vec_angle(const vec2Df64 vec1, const vec2Df64 vec2)
{
    float64 angle =0.0f;

    float64 dot = vec1.x * vec2.x + vec1.y * vec2.y;
    float64 mod1 = sqrt(vec1.x * vec1.x +  vec1.y*vec1.y);
    float64 mod2 = sqrt(vec2.x * vec2.x +  vec2.y*vec2.y);
    angle =  acos(dot/(mod1*mod2));
    assert(!isnan(angle));
    return angle;
}
