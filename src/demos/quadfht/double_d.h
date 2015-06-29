#ifndef DOUBLE_DOUBLE_H_
#define DOUBLE_DOUBLE_H_

//#include "sys_cfg.h"

typedef struct {double hi,lo;} Quad;

void DD_Cmpr(Quad *a, Quad *b, int *ic);
#if 0
void DD_Set1(double a,Quad *b);
void DD_Set2(double a,double b,Quad *c);
void DD_Set(Quad *a, Quad *b);
void DD_GetD(Quad *a,double *b);
void DD_Add(Quad *dda,Quad *ddb, Quad *ddc);
void DD_Sub(Quad *dda, Quad *ddb, Quad *ddc);
void DD_Mul(Quad *dda, Quad *ddb, Quad *ddc);
#endif
void DD_MulD(Quad *dda, double db, Quad *ddc);
void DD_MulDD(double da, double db, Quad *ddc);
void DD_Div(Quad *dda, Quad *ddb, Quad *ddc);
void DD_DivD(Quad *dda, double db, Quad *ddc);
void DD_IntFrac(Quad *a, Quad *b, Quad *c);
void DD_Round(Quad *a, Quad *b);
void DD_Sqrt(Quad *a, Quad *b);
void DD_NPwr(Quad *a, int n, Quad*b);
void DD_OutC(Quad *a,char *b);
void DD_Init(void);
void DD_DeInit(void);
void DD_Dump(Quad *Num);
void DD_SetSign(Quad *Num,int Sign);

#if 1
static inline void
DD_Set1(double a,Quad *b)
/*
** This routine converts the DP number A to DD form in B.  All bits of
** A are recovered in B.  However, note for example that if A = 0.1D0 and N
** is 0, then B will NOT be the DD equivalent of 1/10.
*/
{
b->hi=a;
b->lo=0.0;
}

static inline void
DD_Set2(double a,double b,Quad *c)
/*
** This routine converts the DP number A to DD form in B.  All bits of
** A are recovered in B.  However, note for example that if A = 0.1D0 and N
** is 0, then B will NOT be the DD equivalent of 1/10.
*/
{
c->hi=a;
c->lo=b;
}

static inline void
DD_Set(Quad *a, Quad *b)
{
b->hi=a->hi;
b->lo=a->lo;
}

static inline void
DD_GetD(Quad *a,double *b)
{
*b = a->hi;
}


static inline void
DD_Add(Quad *dda,Quad *ddb, Quad *ddc)
{double e,t1,t2;
t1=dda->hi+ddb->hi;
e=t1-dda->hi;
t2=((ddb->hi-e)+(dda->hi-(t1-e)))+dda->lo+ddb->lo;
ddc->hi=t1+t2;
ddc->lo=t2-(ddc->hi-t1);
}

static inline void
DD_Sub(Quad *dda, Quad *ddb, Quad *ddc)
{double e,t1,t2;
t1=dda->hi - ddb->hi;
e=t1-dda->hi;
t2=((-ddb->hi - e) + (dda->hi - (t1 - e))) + dda->lo - ddb->lo;
ddc->hi=t1+t2;
ddc->lo=t2-(ddc->hi - t1);
}

static inline void
DD_Mul(Quad *dda, Quad *ddb, Quad *ddc)
{double a1,a2,b1,b2,cona,conb,c11,c21,c2,e,split,t1,t2;
split=134217729.0; /* 1+2^27 */
cona = dda->hi * split;
conb = ddb->hi * split;
a1 = cona - (cona - dda->hi);
b1 = conb - (conb - ddb->hi);
a2 = dda->hi - a1;
b2 = ddb->hi - b1;
c11= dda->hi * ddb->hi;
c21= (((a1 * b1 - c11) + a1 * b2) + a2 * b1) + a2 * b2;
c2 = dda->hi * ddb->lo + dda->lo * ddb->hi;
t1 = c11+c2;
e  = t1 - c11;
t2 = ((c2 -e) + (c11 - (t1 -e))) + c21 + dda->lo * ddb->lo;
ddc->hi = t1 + t2;
ddc->lo = t2 - (ddc->hi - t1);
}
#endif

/*
On systems with a 'fused multiply add', the following Fortran code can
be coded instead of the generic method above.

--------------------

!   Multiply dda(1) * ddb(1), assuming a fused multiply-add.

c11 = dda(1) * ddb(1)
c21 = dda(1) * ddb(1) - c11

!   Compute dda(1) * ddb(2) + dda(2) * ddb(1) (only high-order word is needed).

c2 = dda(1) * ddb(2) + dda(2) * ddb(1)

!   Compute (c11, c21) + c2 using Knuth's trick, also adding low-order product.

t1 = c11 + c2
e = t1 - c11
t2 = ((c2 - e) + (c11 - (t1 - e))) + c21 + dda(2) * ddb(2)

!   The result is t1 + t2, after normalization.

ddc(1) = t1 + t2
ddc(2) = t2 - (ddc(1) - t1)
*/


#define F_Add(Sum,Num1,Num2)         DD_Add(Num1,Num2,Sum);
#define F_Clear(Num)                 DD_Set1(0.0,Num)
#define F_DivBy2(Quot,Dividend)      DD_DivD(Dividend,2.0,Quot)
#define F_DivByFloat(Quot, Num, Val) DD_DivD(Num,Val,Quot)
#define F_Mul(Prod,Num1,Num2)        DD_Mul(Num1,Num2,Prod)
#define F_Set(Num,V1,V2)             DD_Set2(V1,V2,Num)
//#define F_SetSign(Num,Sign)          (Num)->hi=fabs((Num)->hi)*Sign;(Num)->lo=fabs((Num)->lo)*Sign
#define F_SetSign(Num,Sign)          (Num)->hi=fabs((Num)->hi)*Sign
#define F_Sub(Dif,Num1,Num2)         DD_Sub(Num1,Num2,Dif)
#define InitFixPoint()               DD_Init()
#define DeInitFixPoint()             DD_DeInit()

/*
#define F_Negate(Num)                (Num)->hi=-((Num)->hi);(Num)->lo=-((Num)->lo)
#define F_Abs(Num)                   {(Num)->hi=fabs((Num)->hi);(Num)->lo=fabs((Num)->lo);}
#define F_AddOne(FIXPOINT *Num);
#define F_ShiftR(FIXPOINT *Num);
#define F_SetStr(FIXPOINT *Num,char *Str);
#define F_Dump(char *Str, FIXPOINT *Num);
*/

#endif


