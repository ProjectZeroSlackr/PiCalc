#ifndef FIXPOINT_H_
#define FIXPOINT_H_

#include "sys_cfg.h"

#define FIXPOINT_LEN 14
typedef struct {UINT32 Data[FIXPOINT_LEN];int Sign;} FIXPOINT;
/* The length needs to be:


**   Digits*2
** + Log10(transform length) (for summation)
** + 1 for integer part of the data format
** + 2*Log10(transform length) (for trig)
** + Log2(pow-2 of transform length)
** + 3.5 for good measure
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
void F_SetStr(FIXPOINT *Num,char *Str);


#endif

