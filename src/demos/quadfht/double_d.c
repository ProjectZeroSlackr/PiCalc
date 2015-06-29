#include "pi.h"
#include <float.h>
#include <math.h>
#include "double_d.h"

static int OldFPU=0;

//#define VDouble volatile double
#define VDouble double

void DD_Init(void)
{
OldFPU=_control87(PC_53|RC_NEAR,MCW_PC|MCW_RC);
}

void DD_DeInit(void)
{
if (OldFPU) _control87(OldFPU,~0);
}

#define inline

void
DD_SetSign(Quad *Num,int Sign)
{
Num->hi=fabs(Num->hi)*Sign;
}

inline void
DD_Cmpr(Quad *a, Quad *b, int *ic)
{
if (a->hi < b->hi) *ic=-1;
else if (a->hi == b->hi)
  {
   if (a->lo < b->lo) *ic=-1;
   else if (a->lo == b->lo) *ic=0;
   else *ic=1;
  }
else *ic=1;
}

#if 0
inline void
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

inline void
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

inline void
DD_Set(Quad *a, Quad *b)
{
b->hi=a->hi;
b->lo=a->lo;
}

inline void
DD_GetD(Quad *a,double *b)
{
*b = a->hi;
}


inline void
DD_Add(Quad *dda,Quad *ddb, Quad *ddc)
{VDouble e,t1,t2;
t1=dda->hi+ddb->hi;
e=t1-dda->hi;
t2=((ddb->hi-e)+(dda->hi-(t1-e)))+dda->lo+ddb->lo;
ddc->hi=t1+t2;
ddc->lo=t2-(ddc->hi-t1);
}

inline void
DD_Sub(Quad *dda, Quad *ddb, Quad *ddc)
{VDouble e,t1,t2;
t1=dda->hi - ddb->hi;
e=t1-dda->hi;
t2=((-ddb->hi - e) + (dda->hi - (t1 - e))) + dda->lo - ddb->lo;
ddc->hi=t1+t2;
ddc->lo=t2-(ddc->hi - t1);
}

inline void
DD_Mul(Quad *dda, Quad *ddb, Quad *ddc)
{VDouble a1,a2,b1,b2,cona,conb,c11,c21,c2,e,split,t1,t2;
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

inline void
DD_MulD(Quad *dda, double db, Quad *ddc)
{VDouble a1,a2,b1,b2,cona,conb,c11,c21,c2,e,split,t1,t2;
split=134217729.0; /* 1+2^27 */

cona = dda->hi * split;
conb = db * split;
a1   = cona - (cona - dda->hi);
b1   = conb - (conb - db);
a2   = dda->hi - a1;
b2   = db - b1;
c11  = dda->hi * db;
c21  = (((a1 * b1 - c11) + a1 * b2) + a2 * b1) + a2 * b2;
c2   = dda->lo*db;
t1   = c11+c2;
e    = t1-c11;
t2   = ((c2 - e) + (c11 - (t1 - e))) + c21;
ddc->hi = t1 + t2;
ddc->lo = t2 - (ddc->hi - t1);
}

inline void
DD_MulDD(double da, double db, Quad *ddc)
{VDouble a1,a2,b1,b2,cona,conb,split,s1,s2;
split=134217729.0; /* 1+2^27 */
cona = da * split;
conb = db * split;
a1 = cona - (cona - da);
b1 = conb - (conb - db);
a2 = da - a1;
b2 = db - b1;

s1 = da * db;
s2 = (((a1 * b1 - s1) + a1 * b2) + a2 * b1) + a2 * b2;
ddc->hi=s1;
ddc->lo=s2;
}

inline void
DD_Div(Quad *dda, Quad *ddb, Quad *ddc)
{VDouble a1,a2,b1,b2,cona,conb,c11,c2,c21,e,split,s1,s2;
 VDouble t1,t2,t11,t12,t21,t22;
split=134217729.0; /* 1+2^27 */

s1 = dda->hi / ddb->hi;

cona = s1 * split;
conb = ddb->hi * split;
a1 = cona - (cona - s1);
b1 = conb - (conb - ddb->hi);
a2 = s1 - a1;
b2 = ddb->hi - b1;

c11 = s1 * ddb->hi;
c21 = (((a1 * b1 - c11) + a1 * b2) + a2 * b1) + a2 * b2;
c2  = s1 * ddb->lo;
t1  = c11 + c2;
e   = t1 - c11;
t2  = ((c2 - e) + (c11 - (t1 - e))) + c21;

t12 = t1 + t2;
t22 = t2 - (t12 - t1);

t11 = dda->hi - t12;
e = t11 - dda->hi;
t21 = ((-t12 - e) + (dda->hi - (t11 - e))) + dda->lo - t22;

s2 = (t11 + t21) / ddb->hi;
ddc->hi = s1+s2;
ddc->lo = s2 - (ddc->hi - s1);
}

inline void
DD_DivD(Quad *dda, double db, Quad *ddc)
{VDouble a1,a2,b1,b2,cona,conb,e,split;
 VDouble t1,t2,t11,t12,t21,t22;
split=134217729.0; /* 1+2^27 */

t1=dda->hi / db;
cona = t1 * split;
conb = db * split;
a1 = cona - (cona - t1);
b1 = conb - (conb - db);
a2 = t1 - a1;
b2 = db - b1;

t12 = t1 * db;
t22 = (((a1 * b1 - t12) + a1 * b2) + a2 * b1) + a2 * b2;

t11 = dda->hi - t12;
e = t11 - dda->hi;
t21 = ((-t12 - e) + (dda->hi - (t11 - e))) + dda->lo - t22;

t2 = (t11 + t21) / db;

ddc->hi = t1 + t2;
ddc->lo = t2 - (ddc->hi - t1);
}

inline void
DD_IntFrac(Quad *a, Quad *b, Quad *c)
/*
** Sets B to the integer part of the DD number A and sets C equal to the
** fractional part of A.  Note that if A = -3.3, then B = -3 and C = -0.3.
*/
{Quad con,f,s0,s1;
 VDouble t105,t52;
 int ic;
 t52 =67108864.0*67108864.0; /* 2^52 */
 t105=67108864.0*67108864.0*67108864.0*67108864.0*2.0; /* 2^105 */


if (a->hi == 0.0) {b->hi=b->lo=c->hi=c->lo=0.0;return;}

if (a->hi >= t105) FatalError("DD_IntFrac: Argument is too large.\n");
f.hi=1.0;
f.lo=0.0;
con.hi=t105;
con.lo=t52;
if (a->hi > 0.0)
  {
   DD_Add(a,&con,&s0);
   DD_Sub(&s0,&con,b);
   DD_Cmpr(a,b,&ic);
   if (ic >= 0) DD_Sub(a,b,c);
   else
     {
      DD_Sub(b,&f,&s1);
      b->hi=s1.hi;
      b->lo=s1.lo;
      DD_Sub(a,b,c);
     }
  }
else
  {
   DD_Sub(a,&con,&s0);
   DD_Add(&s0,&con,b);
   DD_Cmpr(a,b,&ic);
   if (ic <= 0) DD_Sub(a,b,c);
   else
     {
      DD_Add(b,&f,&s1);
      b->hi=s1.hi;
      b->lo=s1.lo;
      DD_Sub(a,b,c);
     }
  }
}

inline void
DD_Round(Quad *a, Quad *b)
{Quad con,s0;
 VDouble t105,t52;
 t52 =67108864.0*67108864.0; /* 2^52 */
 t105=67108864.0*67108864.0*67108864.0*67108864.0*2.0; /* 2^105 */
con.hi=t105;
con.lo=t52;

if (a->hi==0.0) {b->hi=b->lo=0.0;return;}
if (a->hi >= t105) FatalError("DD_Floor: Argument is too large.\n");

if (a->hi > 0.0)
  {
   DD_Add(a,&con,&s0);
   DD_Sub(&s0,&con,b);
  }
else
  {
   DD_Sub(a,&con,&s0);
   DD_Add(&s0,&con,b);
  }
}

inline void
DD_Sqrt(Quad *a, Quad *b)
{VDouble t1,t2,t3;
 Quad s0,s1;

t1 =1.0/sqrt(a->hi);
t2 = a->hi*t1;
DD_MulDD(t2,t2,&s0);
DD_Sub(a,&s0,&s1);
t3 = 0.5 * s1.hi * t1;
s0.hi=t2;
s0.lo=0.0;
s1.hi=t3;
s1.lo=0.0;
DD_Add(&s0,&s1,b);
}

void
DD_NPwr(Quad *a, int n, Quad*b)
{
/*
 This computes the N-th power of the DD number A and returns the DD result
 in B.  When N is zero, 1 is returned.  When N is negative, the reciprocal
 of A ^ |N| is returned.

 This routine employs the binary method for exponentiation.
*/
int j,kk,kn,mn,nn;
VDouble cl2,t1;
Quad s0,s1,s2;

cl2=1.4426950408889633;
if (a->hi == 0.0)
  {
   if (n >= 0) {b->hi=b->lo=0.0;return;}
   else FatalError("DD_NPwr: Argument is zero and N is negative or zer.\n");
  }

nn = abs (n);
if (nn==0) {b->hi=1.0;b->lo=0.0;return;}
else if (nn==1) {s2.hi=a->hi;s2.lo=a->lo;}
else if (nn==2) {DD_Mul(a,a,&s2);}
else
  {
   t1 = nn;
   //mn = cl2 * log (t1) + 1.d0 + 1.d-14;
   mn = cl2 * log(t1) + 1.0 + 1.0e-14;
   s0.hi = a->hi;
   s0.lo = a->lo;
   s2.hi = 1.0;
   s2.lo = 0.0;
   kn = nn;

   for (j=1;j<=mn;j++)
     {
      kk=kn/2;
      if (kn != 2*kk) {DD_Mul(&s2,&s0,&s1);s2.hi=s1.hi;s2.lo=s1.lo;}
      kn=kk;
      if (j < mn) {DD_Mul(&s0,&s0,&s1);s0.hi=s1.hi;s0.lo=s1.lo;}
     }
  }

if (n < 0)
  {
   s1.hi=1.0;s1.lo=0.0;
   DD_Div(&s1,&s2,&s0);
   s2.hi=s0.hi;s2.lo=s0.lo;
  }

b->hi=s2.hi;
b->lo=s2.lo;
}

void
DD_OutC(Quad *a,char *b)
{int i,ii,ln,nx;
 int ib[41];
 VDouble t1;
 Quad f,s0,s1;
 char ca[20],*digits;

b--; /* allow for fortran indexing */

ln=40;digits="*0123456789abcdefghijklmnopqrstuvwxyz";

f.hi=10.0;f.lo=0.0;
for (i=1;i<=ln;i++) ib[i]=0;

if (a->hi != 0.0)
  {
   t1=log10(fabs(a->hi));
   if (t1 >= 0.0) nx=t1;
   else           nx=t1-1.0;
   DD_NPwr(&f,nx,&s0);
   DD_Div(a,&s0,&s1);
   if (s1.hi < 0.0) {s1.hi = -s1.hi; s1.lo = -s1.lo;}

   while (1)
     {
      if (s1.hi < 1.0)
        {
         nx--;
         DD_MulD(&s1,10.0,&s0);
         s1.hi=s0.hi;s1.lo=s0.lo;
        }
      else if (s1.hi >= 10.0)
        {
         nx++;
         DD_DivD(&s1,10.0,&s0);
         s1.hi=s0.hi;s1.lo=s0.lo;
        }
      else break;
     }
  }
else
  {
   nx=0;
   s1.hi=s1.lo=0.0;
  }
for (i=1;i<= ln-8;i++)
  {
   ii=s1.hi;
   ib[i]=ii;
   f.hi=ii;
   DD_Sub(&s1,&f,&s0);
   DD_MulD(&s0,10.0,&s1);
  }

for (i=ln-8;i>=2;i--)
  if (ib[i] < 0)
    {
     ib[i]+=10;
     ib[i-1]--;
    }

if (ib[1] < 0) FatalError("DD_OutC: negative leading digit.");

if (ib[ln-8] >= 5)
  {
   ib[ln-9]++;

   for (i=ln-9;i>=2;i--)
     if (ib[i]==10)
       {
        ib[i]=0;
        ib[i-1]++;
       }

   if (ib[1] == 10)
     {
      ib[1]=1;
      nx++;
     }
  }

b[1]=b[2]=b[3]=' ';

if (a->hi < 0.0) b[3]='-';

ii = ib[1];
b[4] = digits[ii+1];
b[5] = '.';
b[ln] = ' ';
b[ln+1]=0;

for (i=2;i<=ln-9;i++)
  {
   ii=ib[i];
   b[i+4]=digits[ii+1];
  }

sprintf(ca,"%+04d",nx);
b[ln-4]='e';
for (i=1;i<=4;i++) b[ln-4+i]=ca[i-1];
}

void DD_Dump(Quad *Num)
{char Str[80];
DD_OutC(Num,Str);
printf("%s\n",Str);
}

