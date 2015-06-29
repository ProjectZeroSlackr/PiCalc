#ifndef FIXPOINT_H_
#define FIXPOINT_H_

#include "sys_cfg.h"

#define FIXPOINT_LEN 80
typedef struct {UINT32 Data[FIXPOINT_LEN];int Sign;} FIXPOINT;

/*
** The Fixed point length is based upon the number of digits you
** are putting into each element. Let's put 256 decimals into it.
** 
** Knuth gives the formula:
**   2^(3K+2L+2.5+lg(K)-M)<=0.5,
**   with: K=Log2(fft_length)
**         L=Log2(number_base)
**         M=data_type_bits.
** 
** In our example, K=30, L=851.  So, based on that, 'M' needs to be at
** least 1800 bits, preferably higher.  Let's say 1820.  Since this code
** works with base 10, converting bases that means 548 decimals.  Since
** the code is putting 8 decimals into each UINT32 data type, that means
** 69 UINT32's.  Things work better when aligned, so let's say 70 words.
** The data format needs at least 1 extra word for the integer part,
** and I really like the idea of a few extra for insurance for the trig,
** so, let's say 80.
**
** As you can see, although there is a formula, at least part of it
** is personal feeling etc.  Above, I say 80, but 72 would probably
** have worked just as well.
*/


void F_Abs(FIXPOINT *Num);
void F_Add(FIXPOINT *Sum,FIXPOINT *Num1, FIXPOINT *Num2);
void F_AddOne(FIXPOINT *Num);
void F_Clear(FIXPOINT *Num);
void F_DivBy2(FIXPOINT *Quot,FIXPOINT *Dividend);
void F_DivByFloat(FIXPOINT *Quot, FIXPOINT *Num, double Val);
void F_Dump(char *Str, FIXPOINT *Num);
void F_Mul(FIXPOINT *Prod,FIXPOINT *Num1, FIXPOINT *Num2);
void F_Negate(FIXPOINT *Num);
void F_Set(FIXPOINT *Num,UINT32 V1,UINT32 V2);
void F_SetSign(FIXPOINT *Num,int Sign);
void F_ShiftR(FIXPOINT *Num);
void F_Sub(FIXPOINT *Dif,FIXPOINT *Num1, FIXPOINT *Num2);
void InitFixPoint(void);
void DeInitFixPoint(void);


/* Not general purpose.  For trig generation only */
void F_Sqrt(FIXPOINT *Root, FIXPOINT *Num);
void F_Divide(FIXPOINT *D, FIXPOINT *Num, FIXPOINT *Denom);

void
SimpleFFTMul(UINT32 *Prod, UINT32 *Num1, UINT32 *Num2, int NumLen,
             int FFTLen, double *FFTNum1, double *FFTNum2);


#endif

