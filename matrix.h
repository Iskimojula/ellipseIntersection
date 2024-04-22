#ifndef MATRIX_H
#define MATRIX_H
#include <math.h>


#include "algoheader.h"
#define MAXMATRIXSIZE     5
//矩阵左乘
void Matrix_left_mul(const matrixf64 *left,
                     const matrixf64 *right,
                     matrixf64 *rest,
                     const int row1,const int col1,
                     const int row2,const int col2);
//矩阵数乘
void Matrix_dot_mul( matrixf64 *multiplicand,
                     const float64 multiplier,
                     const int row1,const int col1);
//转置
void Matrix_transposition(const  matrixf64 *A,
                           matrixf64 *A_,
                           const int row,
                           const int col);
//加法
void Matrix_adder( matrixf64 *augend,
                     const matrixf64 *addend,
                     const int row,const int col);
//复制
void Matrix_assign( matrixf64 *des,
                     const matrixf64 *str,
                     const int row,const int col);
/*任意维度矩阵的行列式，通过修改MAXMATRIXSIZE参数修改最大维度*/
float64 Matrix_det(const matrixf64 *A,
                   const int row,
                   const int col);
/*矩阵A i_,j_ 元素的余子矩阵，row,col为矩阵A的维度，索引从0开始*/
void Matrix_minor(const matrixf64 *A_,
                     const int i_,
                     const int j_,
                     const int row,
                     const int col,
                     matrixf64 *res);
/*求矩阵的迹*/
float64 Matrix_trace(const matrixf64 *A,
                     const int row,
                     const int col);
/*求矩阵的协因数阵*/
void Matrix_cofactor(const matrixf64 *A,
                        const int row,
                        const int col,
                        matrixf64* res);
/*求矩阵的逆矩阵*/
void Matrix_inverse(const matrixf64 *A,
                    const int row,
                    const int col,
                    matrixf64 *res);


/*vector*/
//二维向量逆时针旋转一个角度theta
void Vec_rotate(vec2Df64 *vec,const float64 theta_);
float64 Vec_angle(const vec2Df64 vec1, const vec2Df64 vec2);

#endif // MATRIX_H

