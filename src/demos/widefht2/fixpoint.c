/*
** The fixpoint format is:
** Big endian.  Base 1e8.  Signed.
** First element integer, remainder fractional parts.  (The integer
**  part is really just for internal use.  For the sqrt & Div.)
**
** However, at times I manipulate directly into the data strucute.
** It's easier to do that than to write wrapper functions.
**
** The 'main' part of the format is the fractional part.  Since we
** don't really have an integer part to use, we have to scale
** the data we put into the FFT elements.  (ie: Knuth.)  When we
** put the data into the FIXPOINT, we scale it by number of 'digits'
** we put in, plus the power of two for the transform.  When we
** later extract, we scale by the square of those, of course.
*/
#include "pi.h"
#include "fixpoint.h"
#include "block.h"


//static UINT32 WorkD[FIXPOINT_LEN*2];
static UINT32 WorkD[FIXPOINT_LEN*4];

void
InitFixPoint(void)
{
}

void
DeInitFixPoint(void)
{
}

void
F_Abs(FIXPOINT *Num)
{
Num->Sign=1;
}

void
F_Add(FIXPOINT *Sum,FIXPOINT *Num1, FIXPOINT *Num2)
{int Overflow;
if (Num1->Sign==Num2->Sign)
  {
   Overflow=BlockAdd(Sum->Data,Num1->Data,Num2->Data,0,FIXPOINT_LEN);
   Sum->Sign=Num1->Sign;
  }
else
  {
   if (Num1->Sign < 0)
     {
      Overflow=BlockSub(Sum->Data,Num2->Data,Num1->Data,0,FIXPOINT_LEN);
      Sum->Sign=1;
      if (Overflow)
        {
         Sum->Sign=-1;
         BlockNegate(Sum->Data,0,FIXPOINT_LEN,0);
        }
     }
   else
     {
      Overflow=BlockSub(Sum->Data,Num1->Data,Num2->Data,0,FIXPOINT_LEN);
      Sum->Sign=1;
      if (Overflow)
        {
         Sum->Sign=-1;
         BlockNegate(Sum->Data,0,FIXPOINT_LEN,0);
        }
     }
  }
}

void
F_AddOne(FIXPOINT *Num)
{
Num->Data[0]++;
}

void
F_Clear(FIXPOINT *Num)
{size_t x;
for (x=0;x<FIXPOINT_LEN;x++) Num->Data[x]=0;
Num->Sign=1;
}

void
F_DivBy2(FIXPOINT *Quot,FIXPOINT *Dividend)
{
BlockDivBy(Quot->Data,Dividend->Data,2,0,FIXPOINT_LEN);
}

/*
** Danger... NOT general purpose.  Only powers of two.
*/
void
F_DivByFloat(FIXPOINT *Quot, FIXPOINT *Num, double Val)
{
while (Val > 1)
  {
   F_DivBy2(Quot,Num);
   Val/=2.0;
  }
}

void
F_Dump(char *Str, FIXPOINT *Num)
{size_t x;

  printf("%s : ", Str);
  if (Num->Sign > 0) printf("+"); else printf("-");
  printf("%u.", Num->Data[0]);
  for (x = 1; x < FIXPOINT_LEN/2; x++)
    printf("%08u", Num->Data[x]);
  printf(":");
  for (x = FIXPOINT_LEN/2; x < FIXPOINT_LEN; x++)
    printf("%08u", Num->Data[x]);
  printf("\n");
}

static void
F_FFTMul(UINT32 *Prod, UINT32 *Num1, UINT32 *Num2);

void
F_Mul(FIXPOINT *Prod,FIXPOINT *Num1, FIXPOINT *Num2)
{
BlockClear(WorkD,&WorkD[FIXPOINT_LEN*2-1]);
//BlockSlowMul(WorkD,Num1->Data,Num2->Data,FIXPOINT_LEN);
F_FFTMul(WorkD,Num1->Data,Num2->Data);
BlockCopy(Prod->Data,WorkD+1,FIXPOINT_LEN);
if (Num1->Sign==Num2->Sign) Prod->Sign=1; else Prod->Sign=-1;
}

void
F_Negate(FIXPOINT *Num)
{
BlockNegate(Num->Data,0,FIXPOINT_LEN,0);
if (Num->Sign > 0) Num->Sign=-1; else Num->Sign=1;
}

void
F_Set(FIXPOINT *Num,UINT32 V1,UINT32 V2)
{
F_Clear(Num);
Num->Data[0]=V1;
Num->Data[1]=V2;
}

void
F_SetSign(FIXPOINT *Num,int Sign)
{
Num->Sign=Sign;
}

void
F_ShiftR(FIXPOINT *Num)
{int x;
for (x=FIXPOINT_LEN-1;x>0;x--) Num->Data[x]=Num->Data[x-1];
Num->Data[0]=0;
}

void
F_Sub(FIXPOINT *Dif,FIXPOINT *Num1, FIXPOINT *Num2)
{int Underflow;

if (Num1->Sign!=Num2->Sign)
  {
   BlockAdd(Dif->Data,Num1->Data,Num2->Data,0,FIXPOINT_LEN);
   Dif->Sign=Num1->Sign;
  }
else if ((Num1->Sign > 0) && (Num2->Sign > 0))
  {
   Underflow=BlockSub(Dif->Data,Num1->Data,Num2->Data,0,FIXPOINT_LEN);
   Dif->Sign=1;
   if (Underflow)
     {
      Dif->Sign= -1;
      BlockNegate(Dif->Data,0,FIXPOINT_LEN,0);
     }
  }
else
  {
   Underflow=BlockSub(Dif->Data,Num2->Data,Num1->Data,0,FIXPOINT_LEN);
   Dif->Sign=1;
   if (Underflow)
     {
      Dif->Sign= -1;
      BlockNegate(Dif->Data,0,FIXPOINT_LEN,0);
     }
  }
}

void
F_SetStr(FIXPOINT *Num,char *Str)
/* Simple minded.  Few safety checks. */
{int x;char Dig[10];
if ((strlen(Str) % 8) != 0)
  FatalError("F_SetStr called with odd length.\n%s\n",Str);

F_Clear(Num);
x=0;
while (*Str)
  {
   memset(Dig,0,10);strncpy(Dig,Str,8);
   Num->Data[x]=atoi(Dig);
   Str+=8;
   x++;
   if (x>=FIXPOINT_LEN) return; /* ignore extra */
  }
}

void
SimpleFFTMul(UINT32 *Prod, UINT32 *Num1, UINT32 *Num2, int NumLen,
             int FFTLen, double *FFTNum1, double *FFTNum2);


static void
F_FFTMul(UINT32 *Prod, UINT32 *Num1, UINT32 *Num2)
{
  static double *FFTNum1=NULL;
  static double *FFTNum2=NULL;
  static int FFTLen=0;
  static int NumLen=0;
  static int NumLen2=0;

  StartTimer(CRTTime);
  if (FFTNum1==NULL)
    {
     FFTLen=FIXPOINT_LEN;
     NumLen=FIXPOINT_LEN;
     NumLen2=FIXPOINT_LEN*2;
     FFTLen=1;
     while (FFTLen < NumLen) FFTLen*=2;
     FFTLen*=2;                       /* for 8 digits in num, 4 in fft */
     FFTLen*=2;                       /* for zero padding */
     FFTNum1=(double*)malloc(FFTLen*sizeof(double));
     FFTNum2=(double*)malloc(FFTLen*sizeof(double));
     if ((FFTNum1==NULL) || (FFTNum2==NULL))
       FatalError("Unable to allocate memory for fixedpoint multiplication.\n");
    }

  SimpleFFTMul(Prod,Num1,Num2,NumLen,FFTLen,FFTNum1,FFTNum2);
  StopTimer(CRTTime);
}





