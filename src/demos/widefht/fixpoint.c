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

static void F_FFTMul(UINT32 *Prod, UINT32 *Num1, UINT32 *Num2);

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



static void
F_RevSubInt(INT32 Val,FIXPOINT *Num)
/* Hardwired just for sqrt.  No signs. */
{INT32 Sign;
  Sign = BlockNegate(Num->Data,0,FIXPOINT_LEN,Val);
  if (Sign) BlockNegate(Num->Data,0,FIXPOINT_LEN,0);
}

/*
** d = a/b by computing the reciprocal of b and then multiplying
** that by a.
**
** r = r*(2-br)
** d = a * r
*/
void
F_Divide(FIXPOINT *D, FIXPOINT *Num, FIXPOINT *Denom)
{ size_t SubLen;
  size_t Redo=FIXPOINT_LEN;
  FIXPOINT Work;

  {double Prep;
   Prep=Denom->Data[0]+((double)Denom->Data[1])/1.0e8;
   Prep=1.0/Prep;
   SubLen=2;
   F_Clear(D);
   D->Data[0]=Prep;
   Prep-=floor(Prep);Prep*=1.0e8;
   D->Data[1]=Prep;
  }

  while (SubLen < FIXPOINT_LEN)
    {
      SubLen *= 2;if (SubLen > FIXPOINT_LEN) SubLen = FIXPOINT_LEN;

      F_Mul(&Work, D, Denom);
      F_RevSubInt(2,&Work);
      F_Mul(D, D, &Work);
      if (SubLen == Redo) {SubLen/=2;Redo=0;}
    }

  F_Mul(D, Num, D);
}

/*
** Computes the square root of 'n' by computing the reciprocal and then
** multiplying that by the original number.
** r = r*(3-n*r^2)/2
*/
void
F_Sqrt(FIXPOINT *Root, FIXPOINT *Num)
{
  size_t SubLen=1;
  size_t Redo=FIXPOINT_LEN;
  FIXPOINT Work;

  {double Prep;
   Prep=Num->Data[0]+((double)Num->Data[1])/1.0e8;
   Prep=sqrt(1.0/Prep);
   SubLen=2;
   F_Clear(Root);
   Root->Data[0]=Prep;
   Prep-=floor(Prep);Prep*=1.0e8;
   Root->Data[1]=Prep;
  }

  while (SubLen < FIXPOINT_LEN)
    {
      SubLen *= 2;if (SubLen > FIXPOINT_LEN) SubLen = FIXPOINT_LEN;
      F_Mul(&Work, Root, Root);
      F_Mul(&Work, Num, &Work);
      F_RevSubInt(3,&Work);
      F_Mul(Root,Root,&Work);
      F_DivBy2(Root, Root);
      if (SubLen == Redo) {SubLen/=2;Redo=0;}
    }

  F_Mul(Root,Num,Root);
}

/*
************************************************************
** special FFT based multiplication for FIXPOINT numbers. **
************************************************************
*/

/*
** Reorder complex data array by bit reversal rule.
*/
static void
F_FFTReOrder(double *point, int Len)
{
  int Index, xednI, x;
  double temp;

  xednI = 0;
  for (Index = 0; Index < Len; Index += 2)
    {
      if (Index < xednI)
        {
          temp = point[Index];
          point[Index] = point[xednI];
          point[xednI] = temp;
          temp = point[Index + 1];
          point[Index + 1] = point[xednI + 1];
          point[xednI + 1] = temp;
        }
      x = Len / 2;
      while ((x > 0) && (x <= xednI))
        {
          xednI -= x;
          x /= 2;
        }
      xednI += x;
    }
}


/*
**  This is fairly much a standard Fast Fourier Transform, except
**  I don't do the 'normalization' at the end of the inverse FFT.
**  Due to memory access times, it's faster to do it later, when
**  we are accessing the memory anyway.  It has been tweaked a bit
**  for my compiler, but nothing that should cause your compiler
**  any problems.  It wouldn't be hard to convert it back to array
**  indexing instead of pointers.
*/
static void
F_SpecialFFT(double *point, int Len, int Dir)
{
  int Step, HalfStep;
  double tcos, tsin;
  double Cosine, Sine;
  double *LastPoint, *ep, *bp;
  double *ap, *fp;
  double Theta;

/*
** We are 'Len' of 'complex' data points.  Since each complex data
** point is actually two doubles, that means the 'double' length is
** Len*2
*/
  Len *= 2;

/*
** Reorder complex data vector by bit reversal rule
*/
  F_FFTReOrder(point, Len);

  LastPoint = &point[Len];
  Theta = atan(1.0) * 4.0;
  if (Dir < 0)
    Theta = -Theta;

/*
** Do the transform
*/
  Step = 2;
  while (Step < Len)
    {
      HalfStep = Step;
      Step *= 2;
      ep = &point[HalfStep];

/*
** The first pass of this loop down below will always have
** Cosine=1 and Sine=0, so we can optimize this a bit.
*/
      for (ap = &point[0]; ap < LastPoint; ap += Step)
        {
          double tr, ti, t;
          fp = ap + HalfStep;
          tr = *fp;
          *fp = (t = *ap) - tr;
          *ap = t + tr;

          ti = *(fp + 1);
          *(fp + 1) = (t = *(ap + 1)) - ti;
          *(ap + 1) = t + ti;
        }

      Cosine = 1.0;
      Sine = 0.0;
      tsin = sin(Theta);
      tcos = sin(Theta *= 0.5);
      tcos = -2.0 * tcos * tcos;
      for (bp = &point[2]; bp < ep; bp += 2)
        {
          {
            double temp;
            temp = Cosine;
            Cosine = Cosine * tcos - Sine * tsin + Cosine;
            Sine = Sine * tcos + temp * tsin + Sine;
          }
          for (ap = bp; ap < LastPoint; ap += Step)
            {
              double tr, ti, t;
              fp = ap + HalfStep;
              tr = (ti = *fp) * Cosine - *(fp + 1) * Sine;
              *fp = (t = *ap) - tr;
              *ap = t + tr;

              ti = ti * Sine + *(fp + 1) * Cosine;
              *(fp + 1) = (t = *(ap + 1)) - ti;
              *(ap + 1) = t + ti;
            }
        }
    }
/*
** We would normally 'normalize' the result, but for us, it's
** faster to do it later, when we are accessing the memory anyway.
*/
}

/*
** This is the routine that lets us use a 'Complex' FFT with our
** Real only data.
*/
static void
F_RFFT(double *data, int Len, int Dir)
{
  double Cosine, Sine, tcos, tsin, temp, theta;
  int i, j;
  double NegHalf, PosHalf;

  Len /= 2;
  theta = -(4.0 * atan(1.0)) / Len;
  NegHalf = -0.5;
  PosHalf = 0.5;
  if (Dir > 0)
    {
      NegHalf = 0.5;
      PosHalf = -0.5;
      theta = -theta;
      F_SpecialFFT(data, Len, 1);
    }

  tcos = sin(theta / 2);
  tcos = -2.0 * tcos * tcos;
  tsin = sin(theta);
  Cosine = 1 + tcos;
  Sine = tsin;
  for (i = 2, j = (2 * Len) - i; i < Len; i += 2, j -= 2)
    {
      double t1r, t1i, t2r, t2i, a, b;
      t1r = 0.5 * ((a = data[i]) + (b = data[j]));
      t2i = PosHalf * (a - b);
      t2r = NegHalf * ((a = data[i + 1]) + (b = data[j + 1]));
      t1i = 0.5 * (a - b);
      a = Cosine * t2r;
      b = Sine * t2i;
      data[i] = t1r + a - b;
      data[j] = t1r - a + b;
      a = Cosine * t2i;
      b = Sine * t2r;
      data[i + 1] = a + b + t1i;
      data[j + 1] = a + b - t1i;
      temp = Cosine;
      Cosine = Cosine * tcos - Sine * tsin + Cosine;
      Sine = Sine * tcos + temp * tsin + Sine;
    }

  temp = data[0];
  data[0] = temp + data[1];
  data[1] = temp - data[1];

  if (Dir < 0)
    {
      data[0] /= 2.0;
      data[1] /= 2.0;
      F_SpecialFFT(data, Len, -1);
    }
}

/*
** Do a Fast Fourier Transform based multiplication.  To save some
** time, we are using the fact that our numbers are purely real and
** doing a half sized real FFT for each number, followed by the
** convolution and then the half sized real FFT inverse.
*/
void
SimpleFFTMul(UINT32 *Prod, UINT32 *Num1, UINT32 *Num2, int NumLen,
             int FFTLen, double *FFTNum1, double *FFTNum2)
{
  int x,y;
  double Carry, RawPyramid, Pyramid, PyramidError;
  double inv;
  double MaxFFTError = 0.0;
  int NumLen2=NumLen*2;

  for (y=0,x=NumLen-1;x>=0;x--,y+=2)
    {
     FFTNum1[y  ]=Num1[x] % 10000;
     FFTNum1[y+1]=Num1[x] / 10000;
    }
  for (;y<FFTLen;y++) FFTNum1[y]=0.0;
  F_RFFT(FFTNum1, FFTLen, 1);

/*
** If we are squaring a number, we can save the cost
** of a FFT.
*/
  if (Num1 != Num2)
    {
      for (y=0,x=NumLen-1;x>=0;x--,y+=2)
        {
         FFTNum2[y  ]=Num2[x] % 10000;
         FFTNum2[y+1]=Num2[x] / 10000;
        }
      for (;y<FFTLen;y++) FFTNum2[y]=0.0;
      F_RFFT(FFTNum2, FFTLen, 1);

      /*
      ** Now do the convolution
      */
      FFTNum1[0] = FFTNum1[0] * FFTNum2[0];
      FFTNum1[1] = FFTNum1[1] * FFTNum2[1];
      for (x = 2; x < FFTLen; x += 2)
        {
          double a, b, c, d;
          FFTNum1[x] = (a = FFTNum1[x])     * (b = FFTNum2[x]) -
                       (c = FFTNum1[x + 1]) * (d = FFTNum2[x + 1]);
          FFTNum1[x + 1] = a * d + c * b;
        }
    }
  else
    {
      /*
      ** Now do the convolution
      */
      FFTNum1[0] = FFTNum1[0] * FFTNum1[0];
      FFTNum1[1] = FFTNum1[1] * FFTNum1[1];
      for (x = 2; x < FFTLen; x += 2)
        {
          double a, c;
          a = FFTNum1[x];
          c = FFTNum1[x + 1];
          FFTNum1[x] = a * a - c * c;
          FFTNum1[x + 1] = a * c + c * a;
        }
    }

/*
** Now do an Inverse FFT
*/
  F_RFFT(FFTNum1, FFTLen, -1);

/*
** Now round the results, and release our carries, and store
** the results in the 'prod' array.  Also do the normalization
** we didn't do in the FFT.  And, as usual, it's slightly faster
** to multiply by the reciprocal, instead of doing the division.
*/
  Carry = 0.0;
  inv = 2.0 / FFTLen;
  y=NumLen2-1;
  x=0;
  while (x<FFTLen)
    {UINT32 D1,D2;
      RawPyramid = FFTNum1[x++] * inv + Carry;
      Pyramid = floor(RawPyramid + 0.5);
      PyramidError = fabs(RawPyramid - Pyramid);
      if (PyramidError > MaxFFTError) MaxFFTError = PyramidError;
      Carry = floor(Pyramid / 10000);
      D1 = Pyramid - Carry * 10000;

      RawPyramid = FFTNum1[x++] * inv + Carry;
      Pyramid = floor(RawPyramid + 0.5);
      PyramidError = fabs(RawPyramid - Pyramid);
      if (PyramidError > MaxFFTError) MaxFFTError = PyramidError;
      Carry = floor(Pyramid / 10000);
      D2 = Pyramid - Carry * 10000;

      Prod[y--]=D2*10000+D1;
    }

  /*
  ** Do a bit of 'sanity' error checking.  This value
  ** is based on some testing with a 'test jig' and
  ** personal opinion.  Should be good enough to catch
  ** systems with poor FPUs.  However, if the value ever
  ** actually reaches 0.5, then you know for a FACT that
  ** the FFTMul() has failed.  David Bailey suggested
  ** a value of 0.001, but that's far far to restrictive.
  */
  if (MaxFFTError >= 0.125)
    {
      printf("**WARNING** Max FFT Error: %f\n", MaxFFTError);
      puts("The FFT is approaching its limit.  As long as");
      puts("the FFTLimit is below 8,388,608 the results are");
      puts("probably correct, but you may want to verify");
      puts("the results with an independant value.  If the");
      puts("value is at least 0.5, then the FFTMul() has");
      puts("definetly failed and the results are incorrect.");
      exit(1);
    }
}

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




