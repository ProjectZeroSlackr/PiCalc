/*
** This pi calculation program implements the Salamin / Brent / Gauss
** Arithmetic-Geometric Mean pi formula.
**
** Version: v1.2.5, March 15, 1999.
**
** This version is a 'cleaned up' version of my much older v1.2.0.1,
** that I placed into the public domain on December 18, 1996.
**
** Although I made a number of textual changes (spacing, etc.), the
** only functional changes were:
**
** 1) I removed the program's ability to handle non-powers of two.
**
** 2) I removed FractalMul's ability to handle lengths that
**    were not powers of two.  (Although I did leave in the
**    code, in case you were curious about how to do it.)
**
** 3) I changed the divide and square root routines from checking
**    the last level and possibly redoing it, to simply always
**    redoing the 'nth' level.  This is better, especially since
**    I no longer need to handle non-power of two lengths.
**
** I barely managed to resist fixing 'all' the other inefficiencies
** in here.  It could certainly use it, but that would destroy its
** usefulness as part of my pi tutorial.
**
** This program is placed into the public domain by its author,
** Carey Bloodworth, on March 15, 1999.
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <limits.h>

/*
** Set the FFTLimit here. See comments at end of program for details.
*/
#ifndef FFTLimit
#error You need to set the FFTLimit for your system.  See program for details.
#endif

/*
** Warn the user if he tries to run it on a 16 bit compiler.
*/
#if SHRT_MAX == INT_MAX
#error Forget it!  A 16 bit compiler would not have enough memory.
#endif

/* Which iteration to redo for accuracy */
#define REDO_LEN (Len/16)

typedef short int Short;
typedef long int Long;

Short *AGM_A, *AGM_B, *AGM_Sum, *RRSWork, *MulProd, *AGMWork, *OldRoot;
double *FFTNum1, *FFTNum2;
int MaxIter;

/*
** RRSWork is explicitly used in Reciprocal() and RSqrt()
**   It is implicitly used in Divide() and Sqrt()
**   RRSWork stands for _R_eciprocal() & _RS_qrt() _Work_ var.
** OldRoot is used in the Sqrt() to hold the previous square
**   root.  This saves us some time in RSqrt().
**   It's also used in Divide(), since we are done with it by then.
** MulProd is where the multiplication puts its double
**   length answer.  It's used in lots of places, of course.
**   It's also big enough for all of FractalMul()'s needs.
** AGMWork is used in the AGM iteration as a temp variable.
** FFTNum1 and FFTNum2 are the two 'massive' arrays used in the
**   Fast Fourier Transform used to multiply the numbers.
** MaxIter is the maximum 'reasonable' number of iterations that
**   the reciprocal, square root, and agm iteration could take.
**   Anything beyond that is considered a fatal error.
*/

/*
** Just a predeclaration.
*/
void Mul(Short * prod, Short * a, Short * b, int Len);

/*
*******************************
**      Low level stuff      **
*******************************
*/

/*
** The lower of two int numbers
*/
int
Min(int Num1, int Num2)
{
  if (Num1 < Num2)
    return Num1;
  return Num2;
}

/*
** a=b
*/
void
Copy(Short *Dest, Short *Src, int Len)
{
  while (Len--)
    Dest[Len] = Src[Len];
}

/*
** Find where two numbers differ.  Useful for checking convergence.
*/
int
DiffWhere(Short *Num1, Short *Num2, int Len)
{
  int x;

  for (x = 0; x < Len; x++)
    if (Num1[x] != Num2[x])
      return x;
  return Len;
}

/*
** Is 'Num' exactly zero?
*/
int
IsZero(Short *Num, int Len)
{
  while (Len--)
    if (*Num++ != 0)
      return 0;
  return 1;
}

/*
** This rounds and copies a 2x mul result into a normal result
** Our number format will never have more than one unit of integer,
** and after a mul, we have two, so we need to fix that.
*/
void
Round2x(Short *Dest, Short *Src, int Len)
{
  int x;
  Short carry;

  carry = 0;
  if (Src[Len + 1] >= 5000)
    carry = 1;
  for (x = Len; x > 0; x--)
    {
      carry += Src[x];
      Dest[x - 1] = (Short)(carry % 10000);
      carry = (Short)(carry / 10000);
    }
}

/*
** n = n / 2
*/
void
DivBy2(Short *Num, int Len)
{
  int x;
  Short temp;

  temp = 0;
  for (x = 0; x < Len; x++)
    {
      temp = (Short)(Num[x] + temp * 10000);
      Num[x] = (Short)(temp / 2);
      temp = (Short)(temp % 2);
    }
}

int
IsPow2(int Num)
{
  return ((Num & -Num) == Num);
}


/*
*******************************
**      Output routines      **
*******************************
*/

/*
** Print out a simple straight representation of our big numbers
** Used only in debugging sessions.
*/
void PrintShort2(char *Str, Short *Num, int Len)
{
  int x;

  printf("%s ", Str);
  printf("%u.", Num[0]);
  for (x = 1; x < Len; x++)
    printf("%04u", Num[x]);
  printf("\n");
}


/*
** Print out a formated representation of our big numbers
**
** The formatting is based on how the Gutenberg PIMIL10.TXT
** is formatted, allowing easy comparision of results.
*/
void
PrintShort(char *Str, Short *Num, int Len)
{
  int x;
  int printed = 0;

  printf("%s", Str);

  printf("%u.\n\n", Num[0]);

  for (x = 1; x < Len; x++)
    {
      printf("%02u", Num[x] / 100);
      printed += 2;
      if ((printed % 1000) == 0)    {printf("\n\n\n\n");printed = 0;}
      else if ((printed % 50) == 0)  printf("\n");
      else if ((printed % 10) == 0)  printf(" ");
      printf("%02u", Num[x] % 100);
      printed += 2;
      if ((printed % 1000) == 0)    {printf("\n\n\n\n");printed = 0;}
      else if ((printed % 50) == 0)  printf("\n");
      else if ((printed % 10) == 0)  printf(" ");
    }
  printf("\n");

}


/*
*******************************
**         Basic math        **
*******************************
*/


/*
** Sum = Num1 + Num2
** Just like a regular add, except the carry can ripple on up higher,
** beyond where we end our addition.
*/
Short
RippleAdd(Short *Limit, Short *Sum, Short *Num1, Short *Num2, int Len)
{
  Short temp, Carry;

  Carry = 0;
  Num1 += Len;
  Num2 += Len;
  Sum += Len;
  while (Len--)
    {
      temp     = (Short)(*(--Num1) + *(--Num2) + Carry);
      *(--Sum) = (Short)(temp % 10000);
      Carry    = (Short)(temp / 10000);
    }

  while ((Sum > Limit) && (Carry != 0))
    {
      temp  = (Short)(*(--Sum) + Carry);
      *Sum  = (Short)(temp % 10000);
      Carry = (Short)(temp / 10000);
    }
  return Carry;
}

/*
** Sum = Num1 + Num2
*/
Short
Add(Short *Sum, Short *Num1, Short *Num2, int Len)
{
  return RippleAdd(Sum, Sum, Num1, Num2, Len);
}

/*
** Dif = Num1-Num2
** Just like a regular sub, except the borrow can ripple on up higher,
** beyond where we end our subtraction.
*/
Short
RippleSub(Short *Limit, Short *Dif, Short *Num1, Short *Num2, int Len)
{
  Short temp, Borrow;

  Borrow = 0;
  Num1 += Len;
  Num2 += Len;
  Dif += Len;
  while (Len--)
    {
      temp = (Short)(*(--Num1) - *(--Num2) - Borrow);
      Borrow = 0;
      if (temp < 0)
        {
          Borrow = 1;
          temp = (Short)(temp + 10000);
        }
      *(--Dif) = temp;
    }

  while ((Dif > Limit) && (Borrow != 0))
    {
      temp = (Short)(*(--Dif) - Borrow);
      Borrow = 0;
      if (temp < 0)
        {
          Borrow = 1;
          temp = (Short)(temp +10000);
        }
      *Dif = temp;
    }
  return Borrow;
}

/*
** Dif = Num1-Num2
*/
Short
Sub(Short *Dif, Short *Num1, Short *Num2, int Len)
{
  return RippleSub(Dif, Dif, Num1, Num2, Len);
}

/*
** Num = 0 - Num
*/
void
Negate(Short *Num, int Len)
{
  int x;
  Short d, Borrow;

  Borrow = 0;
  for (x = Len - 1; x >= 0; x--)
    {
      d = (Short)(0 - Num[x] - Borrow);
      Borrow = 0;
      if (d < 0)
        {
          Borrow = 1;
          d = (Short)(d + 10000);
        }
      Num[x] = d;
    }
}

/*
*******************************
**       Multiplication      **
*******************************
*/

/*
** Prod = Num1*Num2.
** Simple, slow multiplication method.  The space for 'prod'
** should of course be double length.
*/
void
SlowMul(Short *Prod, Short *Num1, Short *Num2, int Len)
{
  Long p;
  int Ndx1, Ndx2, NdxP;

  for (NdxP = 0; NdxP < Len * 2; NdxP++)
    Prod[NdxP] = 0;
  for (Ndx2 = Len - 1; Ndx2 >= 0; Ndx2--)
    {
      if (Num2[Ndx2] == 0)
        continue;
      NdxP = Ndx2 + Len;
      p = 0;
      for (Ndx1 = Len - 1; Ndx1 >= 0; Ndx1--)
        {
          p = (Long) Num1[Ndx1] * (Long) Num2[Ndx2] + p + (Long) Prod[NdxP];
          Prod[NdxP] = (Short) (p % 10000);
          p = p / 10000;
          NdxP--;
        }
      while ((p) && (NdxP >= 0))
        {
          p += (Long) Prod[NdxP];
          Prod[NdxP] = (Short) (p % 10000);
          p = p / 10000;
          NdxP--;
        }
    }
}

#if 0
/*
** This is based on the simple Fractal / Divide and Conquer
** O(n^1.585) method.
**
** It's fairly simple.  a*b is: a1b1(B^2+B)+(a1-a2)(b2-b1)B+a2b2(B+1)
**
** The simple formula above really only works when the length of the
** numbers is a power of two.  The implementation below can deal
** with numbers of any length, as long as both 'a' and 'b' are the
** same length.  Doing different lengths would have complicated it
** more, and wasn't needed anyway.
**
** For lengths of powers of two, you need 4*SIZE storage.  For other
** sizes, you need 5*SIZE.  This is guaranteed enough by the nature
** of the formula, and the way I do the storage.
**
** Also, the program attempts to divide the numbers so that at
** least one of them will be a power of two.  This lets the FFT
** handle at least this part efficiently.  Overall, it's better to
** calculate digits of pi that are powers of two.
*/
void
FractalMul(Short *prod, Short *a, Short *b, Indexer Len)
{
  int SignA, SignB;
  Short *offset;
  Short *WorkArea;
  Short *a1, *a2, *b1, *b2;
  Indexer x, HLen1, HLen2, OddLen;
  /*
  ** HLen1 is the length of the first sub-part of the numbers
  ** HLen2 is the length of the second part.  HLen1 is always
  ** the longest if the length isn't even.  HLen stands for
  ** 'Half Length'.
  ** OddLen keeps track of whether it was an even length.
  **
  ** The a1,a2,b1,b2 are used just so you can relate to the
  ** formula given above.
  */

  /*
  ** Figure out how to divide the number.  Preferably, so that
  ** at least one of the numbers will be a power of two.
  ** We also don't want the second number to be either 0 or 1.
  */
  HLen1 = Len;
  /*
  ** Make it a power of two by stripping off all
  ** but the high bit.
  */
  while (!IsPow2(HLen1))
    HLen1 ^= (HLen1 & -HLen1);
  if (HLen1 == Len)           HLen1 /= 2;
  else if (HLen1 == Len - 1)  HLen1 = Len - HLen1 / 2;
  HLen2 = Len - HLen1;
  OddLen = HLen1 - HLen2;
  WorkArea = prod + 2 * Len;
  a1 = a;
  a2 = a1 + HLen1;
  b1 = b;
  b2 = b1 + HLen1;

  /*
  ** Now do the (b2-b1) part of the formula.
  */
  for (x = 0; x < OddLen; x++) prod[x] = 0;
  Copy(prod + OddLen, b2, HLen2);
  SignB = Sub(prod, prod, b1, HLen1);
  if (SignB) Negate(prod, HLen1);

  /*
  ** Now do the (a1-a2) part of the formula.
  ** I don't like damaging the number, even though we'll
  ** fix it later.  It does however reduce the working space
  ** you need.
  */
  SignA = RippleSub(a1, a1 + OddLen, a1 + OddLen, a2, HLen2);
  if (SignA) Negate(a1, HLen1);

  /*
  ** Now multiply those two differences (a1-a2)*(b2-b1)
  ** and then, depending on the signs of the two subtractions
  ** either add or subtract it into our answer, which so
  ** far is still zero.  By doing the middle part of the
  ** equation first, I can save some memory.
  */
  Mul(WorkArea, prod, a1, HLen1);
  for (x = 0; x < 2 * Len; x++) prod[x] = 0;
  offset = prod + HLen1 - OddLen;
  if (SignA == SignB) RippleAdd(prod, offset, offset, WorkArea, HLen1 * 2);
  else                RippleSub(prod, offset, offset, WorkArea, HLen1 * 2);

  /*
  ** Now we have to restore the 'a1' we damaged.
  */
  if (SignA) Negate(a1, HLen1);
  RippleAdd(a1, a1 + OddLen, a1 + OddLen, a2, HLen2);

  /*
  ** Now we do the first part of the formula, a1*b1
  ** We then add that product into our answer.
  ** We have to add it in two places, remember?
  */
  Mul(WorkArea, a1, b1, HLen1);
  offset = prod + HLen2;
  RippleAdd(prod, offset, offset, WorkArea, HLen1 * 2);
  Add(prod, prod, WorkArea, HLen1 * 2);

  /*
  ** Now we do the last part of the formula, a2*b2
  ** We then add that product into our answer.
  ** We have to add it in two places, remember?
  */
  Mul(WorkArea, a2, b2, HLen2);
  offset = prod + HLen1 + OddLen;
  RippleAdd(prod, offset, offset, WorkArea, HLen2 * 2);
  offset = prod + HLen1 * 2;
  RippleAdd(prod, offset, offset, WorkArea, HLen2 * 2);
}
#endif


/*
** This is based on the simple Fractal / Divide and Conquer / Karatsuba
** O(n^1.585) method.
**
** It's fairly simple.  a*b is: a1b1(B^2+B)+(a1-a2)(b2-b1)B+a2b2(B+1)
**
** You need 4*SIZE storage for the product and the working space.
*/
void
FractalMul(Short *Prod, Short *Num1, Short *Num2, int Len)
{int x, HLen, Sign1, Sign2;
 Short *offset;
 Short *Work=Prod+2*Len;

  HLen = Len/2;

  /* Do x=Right(Num1)-Left(Num1) y=Left(Num2)-Right(Num2) */
  Sign1 = Sub(Prod, Num1+HLen, Num1, HLen);     if (Sign1) Negate(Prod, HLen);
  Sign2 = Sub(Prod+HLen, Num2, Num2+HLen, HLen);if (Sign2) Negate(Prod+HLen, HLen);

  Mul(Work, Prod, Prod+HLen, HLen);

  for (x=0;x<Len*2;x++) Prod[x]=0;
  offset = Prod + HLen;
  if (Sign1 == Sign2)  RippleAdd(Prod,offset,offset,Work,Len);
  else                 RippleSub(Prod,offset,offset,Work,Len);

  Mul(Work, Num1, Num2, HLen);
  offset = Prod + HLen;RippleAdd(Prod,offset,offset,Work,Len);
  Add(Prod, Prod, Work, Len);

  Mul(Work, Num1 + HLen, Num2 + HLen, HLen);
  offset = Prod + HLen;RippleAdd(Prod,offset,offset,Work,Len);
  offset = Prod + Len; RippleAdd(Prod,offset,offset,Work,Len);
}



/*
*******************************
**     FFT Multiplication    **
*******************************
*/


/*
** Reorder complex data array by bit reversal rule.
*/
void
FFTReOrder(double *point, int Len)
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
void
SpecialFFT(double *point, int Len, int Dir)
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
  FFTReOrder(point, Len);

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
void
RFFT(double *data, int Len, int Dir)
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
      SpecialFFT(data, Len, 1);
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
      SpecialFFT(data, Len, -1);
    }
}

/*
** Do a Fast Fourier Transform based multiplication.  To save some
** time, we are using the fact that our numbers are purely real and
** doing a half sized real FFT for each number, followed by the
** convolution and then the half sized real FFT inverse.
*/
void
FFTMul(Short *Prod, Short *Num1, Short *Num2, int Len)
{
  int x, Len2 = Len * 2;
  double Carry, RawPyramid, Pyramid, PyramidError;
  double inv;
  double MaxFFTError = 0.0;

/*
** This is a radix-2 FFT, so the length has to be a power of two.
*/
  if (!IsPow2(Len))
    {
      printf("The FFT size is not a power of two\n");
      exit(EXIT_FAILURE);
    }

/*
** Make sure we aren't trying to do too big of a FFT.
** Note that this is 'element' size, not digits, like FFTLimit.
*/
  if ((long) Len > 2097152L)
    {
      printf("Too large of a FFT.  I fail over 8,388,608.\n");
      exit(EXIT_FAILURE);
    }

  for (x = 0; x < Len; x++)    FFTNum1[x] = Num1[Len - x - 1];
  for (x = Len; x < Len2; x++) FFTNum1[x] = 0.0;
  RFFT(FFTNum1, Len2, 1);

/*
** If we are squaring a number, we can save the cost
** of a FFT.
*/
  if (Num1 != Num2)
    {
      for (x = 0; x < Len; x++)    FFTNum2[x] = Num2[Len - x - 1];
      for (x = Len; x < Len2; x++) FFTNum2[x] = 0.0;
      RFFT(FFTNum2, Len2, 1);

      /*
      ** Now do the convolution
      */
      FFTNum1[0] = FFTNum1[0] * FFTNum2[0];
      FFTNum1[1] = FFTNum1[1] * FFTNum2[1];
      for (x = 2; x < Len2; x += 2)
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
      for (x = 2; x < Len2; x += 2)
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
  RFFT(FFTNum1, Len2, -1);

/*
** Now round the results, and release our carries, and store
** the results in the 'prod' array.  Also do the normalization
** we didn't do in the FFT.  And, as usual, it's slightly faster
** to multiply by the reciprocal, instead of doing the division.
*/
  Carry = 0.0;
  inv = 1.0 / Len;
  for (x = Len2; x > 0; x--)
    {
      RawPyramid = FFTNum1[Len2 - x] * inv + Carry;
      Pyramid = floor(RawPyramid + 0.5);
      PyramidError = fabs(RawPyramid - Pyramid);
      if (PyramidError > MaxFFTError)
        {
          MaxFFTError = PyramidError;
/*        printf("New Max FFT Error=%f\n",MaxFFTError);*/
        }
      Carry = floor(Pyramid / 10000);
      Prod[x - 1] = Pyramid - Carry * 10000;
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
    }
}


/*
** FINALLY we get finished with the three specialized multiply
** routines and get to the generic one we actually call!
*/

/*
** This routine is the master multiplication routine.  It decides
** which multiplication routine is most appropriate for the size
** of the numbers.  The regular slow O(N^2) method, the fractal
** O(N^1.585) method, or the FFT O(N Log2(N)) method.
*/
void
Mul(Short *Prod, Short *Num1, Short *Num2, int Len)
{
  if (Len <= 16)
      SlowMul(Prod, Num1, Num2, Len);
  /*
  ** Do the FFT based multiplication, if the size is within
  ** what we can handle, and the length is a power of two.
  ** Note: We have to convert our FFTLimit of digits to the
  ** number of array elements used.  Simple division by 4.
  */
  else if ((Len <= FFTLimit / 4) && IsPow2(Len))
    FFTMul(Prod, Num1, Num2, Len);
  else
    FractalMul(Prod, Num1, Num2, Len);
}

/*
*******************************
**  Division and square root **
*******************************
*/


/*
** Computes the square root of 'n' by computing the reciprocal and then
** multiplying that by the original number.
** r = r*(3-n*r^2)/2
*/
void
Sqrt(Short *Root, Short *Num, int Len)
{
  int x, SubLen=4;
  static int PreSubLen = 2;
  int Redo=REDO_LEN;

  if      ((Num[0]== 0) && (Num[1]==5000) && (Num[2]==   0) && (Num[3]==   0))
    {OldRoot[0]=1;OldRoot[1]=4142;OldRoot[2]=1356;OldRoot[3]=2373;}
  else if ((Num[0]== 0) && (Num[1]==7071) && (Num[2]== 678) && (Num[3]==1186))
    {OldRoot[0]=1;OldRoot[1]=1892;OldRoot[2]= 711;OldRoot[3]=5002;}
  else if ((Num[0]== 0) && (Num[1]==7177) && (Num[2]==4998) && (Num[3]==6377))
    {OldRoot[0]=1;OldRoot[1]=1803;OldRoot[2]=5706;OldRoot[3]=4195;}
  else if ((Num[0]== 0) && (Num[1]==7177) && (Num[2]==7001) && (Num[3]== 976))
    {OldRoot[0]=1;OldRoot[1]=1803;OldRoot[2]=4059;OldRoot[3]=9073;}
  else if ((Num[0]== 0) && (Num[1]==7177) && (Num[2]==7001) && (Num[3]==1046))
    {
     OldRoot[0]=1;OldRoot[1]=1803;OldRoot[2]=4059;OldRoot[3]=9016;
     SubLen=PreSubLen;PreSubLen*=2;
    }
  else
    {
     PrintShort("Unknown Square root: ",Num,16);
     exit(EXIT_FAILURE);
    }

  while (SubLen < Len)
    {
      fputc('S', stderr);
      SubLen *= 2;if (SubLen > Len) SubLen = Len;

      Mul(MulProd, OldRoot, OldRoot, SubLen);
      Round2x(RRSWork, MulProd, SubLen);
      Mul(MulProd, Num, RRSWork, SubLen);
      Round2x(RRSWork, MulProd, SubLen);

      MulProd[0] = 3;
      for (x = 1; x < SubLen; x++) MulProd[x] = 0;
      Sub(RRSWork, MulProd, RRSWork, SubLen);

      Mul(MulProd, OldRoot, RRSWork, SubLen);
      Round2x(OldRoot, MulProd, SubLen);
      DivBy2(OldRoot, SubLen);
      if (SubLen == Redo) {SubLen/=2;Redo=0;}
    }
  Mul(MulProd, Num, OldRoot, Len);
  Round2x(Root, MulProd, Len);
}

/*
** d = a/b by computing the reciprocal of b and then multiplying
** that by a.
**
** r = r*(2-br)
** d = a * r
*/
void
Divide(Short *d, Short *Num, Short *Denom, int Len)
{ Short *R=OldRoot;
  int x, SubLen;
  int Redo=REDO_LEN;

  /*
  ** Estimate our reciprocal.  I can cheat because
  ** I already know what it's supposed to be.
  */
  R[0] = 1; R[1] = 942; R[2] = 1980; R[3] = 7613; SubLen = 4;

  while (SubLen < Len)
    {
      fputc('R', stderr);
      SubLen *= 2;if (SubLen > Len) SubLen = Len;

      Mul(MulProd, R, Denom, SubLen);
      Round2x(RRSWork, MulProd, SubLen);

      MulProd[0] = 2;
      for (x = 1; x < SubLen; x++) MulProd[x] = 0;
      Sub(RRSWork, MulProd, RRSWork, SubLen);

      Mul(MulProd, R, RRSWork, SubLen);
      Round2x(R, MulProd, SubLen);
      if (SubLen == Redo) {SubLen/=2;Redo=0;}
    }

  Mul(MulProd, Num, R, Len);
  Round2x(d, MulProd, Len);
}


/*
*******************************
**      The AGM itself       **
*******************************
*/

void
AGM(int Len)
{
  int x;
  double Pow2, tempd, carryd;
  int AGMIter, loops;
  time_t LoopTime;

  Pow2 = 4.0;

  fprintf(stderr, "Initialization\n");
  LoopTime = time(NULL);
  for (x = 0; x < Len; x++)
    AGM_A[x] = AGM_B[x] = AGM_Sum[x] = AGMWork[x] = OldRoot[x] = 0;
  AGM_A[0] = 1;
  AGM_Sum[0] = 1;
  AGMWork[1] = 5000;
  Sqrt(AGM_B, AGMWork, Len);
  fprintf(stderr, " Time= %0.0f\n", difftime(time(NULL), LoopTime));

  /*
  printf("AGMIter %d\n", AGMIter);
  PrintShort2("a AGM: ", a, Len);
  PrintShort2("b AGM: ", b, Len);
  PrintShort2("C sum: ", c, Len);
  */

  loops = AGMIter = MaxIter;
  do
    {
      fprintf(stderr, "Pass %d\n", AGMIter - 5);
      LoopTime = time(NULL);
      /* w = (a-b)/2      */
      Sub(AGMWork, AGM_A, AGM_B, Len);
      DivBy2(AGMWork, Len);

      /* m = w*w          */
      Mul(MulProd, AGMWork, AGMWork, Len);
      Round2x(AGMWork, MulProd, Len);

      /* m = m* w^(J+1)   */
      carryd = 0.0;
      for (x = Len - 1; x >= 0; x--)
        {
          tempd = Pow2 * AGMWork[x] + carryd;
          carryd = floor(tempd / 10000.0);
          AGMWork[x] = tempd - carryd * 10000.0;
        }
      Pow2 *= 2.0;
      /* c = c - m        */
      Sub(AGM_Sum, AGM_Sum, AGMWork, Len);

      /* See if it's done */
      if (IsZero(AGMWork, Len))
        AGMIter = 0;

      /* Save some time?  */
      if (AGMIter != 0)
        {
          /* m = a*b          */
          Mul(MulProd, AGM_A, AGM_B, Len);
        }
      /* a = (a+b)/2      */
      Add(AGM_A, AGM_A, AGM_B, Len);
      DivBy2(AGM_A, Len);

      /* Save some time?  */
      if (AGMIter != 0)
        {
          Round2x(AGMWork, MulProd, Len);
          /* b = sqrt(a*b)    */
          Sqrt(AGM_B, AGMWork, Len);
        }

      /*
      printf("AGMIter %d\n", AGMIter);
      PrintShort2("a AGM: ", AGM_A, Len);
      PrintShort2("b AGM: ", AGM_B, Len);
      PrintShort2("C sum: ", AGM_Sum, Len);
      */

      loops--;
      fprintf(stderr, " Time= %0.0f\n", difftime(time(NULL), LoopTime));
    }
  while (--AGMIter > 0);

  if (loops == 0)
    {
      printf("AGM iter failure.  Loops=0\n");
      PrintShort("a AGM: ", AGM_A, Len);
      PrintShort("b AGM: ", AGM_B, Len);
      PrintShort("c Sum: ", AGM_Sum, Len);
      exit(EXIT_FAILURE);
    }

  Mul(MulProd, AGM_A, AGM_A, Len);
  Round2x(AGM_A, MulProd, Len);
  carryd = 0.0;
  for (x = Len - 1; x >= 0; x--)
    {
      tempd = 4.0 * AGM_A[x] + carryd;
      carryd = floor(tempd / 10000.0);
      AGM_A[x] = tempd - carryd * 10000.0;
    }

  Divide(AGM_B, AGM_A, AGM_Sum, Len);
}

void
usage(void)
{
  puts("This program requires the number of digits to calculate.");
  puts("This should be at least 16 digits.  The length must be a");
  puts("power of two.");
  exit(EXIT_FAILURE);
}

int
main(int argc, char *argv[])
{
  time_t StartTime, EndTime;
  int Len, Digits;

  if (argc < 2) usage();

  Digits = atoi(argv[1]);
  if (!IsPow2(Digits)) usage();
  Len = (Digits + 3) / 4;

  if (Len < 4) usage();
  /*
  fprintf(stderr,"data Len is %d, # digits is %d\n",Len,Digits);
  */

  /*
  ** No iteration (AGM, RSqrt, or Reciprocal) should ever
  ** take more than this.  If so... then it's an error!
  */
  MaxIter = log(Digits) / log(2.0) + 5;

  AGM_A = (Short *) calloc(Len, sizeof(Short));
  AGM_B = (Short *) calloc(Len, sizeof(Short));
  AGM_Sum = (Short *) calloc(Len, sizeof(Short));
  RRSWork = (Short *) calloc(Len, sizeof(Short));
  AGMWork = (Short *) calloc(Len, sizeof(Short));
  OldRoot = (Short *) calloc(Len, sizeof(Short));
  if (IsPow2(Len) && (Digits <= FFTLimit))
    MulProd = (Short *) calloc(Len * 2, sizeof(Short));
  else
    MulProd = (Short *) calloc(Len * 4, sizeof(Short));

  FFTNum1 = (double *) calloc(Min(FFTLimit / 2, Len * 2), sizeof(double));
  FFTNum2 = (double *) calloc(Min(FFTLimit / 2, Len * 2), sizeof(double));

  if ((AGM_A == NULL) || (AGM_B == NULL) || (AGM_Sum == NULL) ||
      (RRSWork == NULL) || (AGMWork == NULL) ||
      (OldRoot == NULL) || (MulProd == NULL) ||
      (FFTNum1 == NULL) || (FFTNum2 == NULL))
    {
      printf("Unable to allocate that much memory\n");
      exit(EXIT_FAILURE);
    }

  StartTime = time(NULL);

  /*
  ** Now actually do the Arithmetic-Geometric-Mean
  */
  AGM(Len);

  EndTime = time(NULL);

  PrintShort("", AGM_B, Len);

  printf("\nTotal Execution time: %0.0f\n",
         difftime(EndTime, StartTime));


  free(AGM_A);
  free(AGM_B);
  free(AGM_Sum);
  free(RRSWork);
  free(AGMWork);
  free(OldRoot);
  free(MulProd);
  free(FFTNum1);
  free(FFTNum2);

  return EXIT_SUCCESS;
}


/*
** The FFTLimit is the largest FFT based multiplication you can do
** in physical memory.  You do *NOT* want disk based virtual memory
** to kick in in the middle of a FFT, although virtual memory usage
** elsewhere is some what tolerable.
**
** This is a matter of obtaining the best performance possible, not
** of program ability.  Even if you set this too low or too high
** you can still compute a large number of digits of pi.
**
** You need to set the FFTLimit according to the following guide:
**
** With virtual memory (either OS or from the compiler), basically,
** you should divide the physical memory given to your program by
** 4, then round down to the next power of two, and that's the number
** of digits that can be multiplied together in a single 'chunk'.
** Numbers larger than that will first be broken into smaller chunks
** that the FFT can handle.
**
** For a 4meg (4194304) computer, you'd have 4194304/4=1,048,576
** But, even though that's a power of two, you have to gown down to
** the next one, so you'll still have enough room for your OS and
** the program.  So it'd be 524,288
**
** Assuming virtual memory:
** A     4 meg computer can manage a threshold of   524,288.
** A  5- 8 meg computer can manage a threshold of 1,048,576.
** A  9-16 meg computer can manage a threshold of 2,097,152.
** A 17-32 meg computer can manage a threshold of 4,194,304.
** A 33-?? meg computer can manage a threshold of 8,388,608.
**
** If you don't have virtual memory, then you just have to determine
** yourself just how much you can spare to the FFT.  16 megs of
** memory is enough to do 1,048,576 digits of pi entirely within
** physical memory and with the FFTLimit set at 1,048,576.  You could
** just set it extremely high, and if you try to compute too many,
** it'll tell you when you try.  However, since the FFT uses most of
** the memory, reducing this to a level smaller than the number of
** digits you are computing will actually let you compute more digits.
** Reducing the FFTLimit to '0' would turn off the FFT (and it's large
** memory usage) and let you compute even more digits, but increase
** the run time considerably.  Without the FFT, the program isn't even
** really practical.  Other programs, such as a classic arctangent,
** like in the C Snippets, would be better.  I can't help you with
** this.  You'll have to decide for yourself.
**
** If you've got enough memory, it is definetly best to set it as
** high as possible.  This FFT based multiplication implementation has
** a *MAX* size of 8,388,608.  Beyond that, it fails due to the limited
** precision of an 8 byte 'double' (ie: it overflows.)
**
** Examples:
#define FFTLimit   524288
#define FFTLimit  1048576
#define FFTLimit  2097152
#define FFTLimit  4194304
#define FFTLimit  8388608
**
** Total Memory usage:
** There are 3 actual formula variables of length   SIZE
** There are 3 temporary / work variables of length SIZE
** There is 1 multiplication variable of length     SIZE*5
**   (Could be SIZE*4 if SIZE is a power of two, or just SIZE*2
**    if it's power of two, and it fits entirely within the FFT.)
**
** SIZE is (DIGITS_WANTED+3)/4.
** Each array element is a short 16 bit integer.
**
** So, to compute 1024k digits totally within memory and the numbers
** fitting entirely within the FFT, SIZE would be 256k, you'd have 6
** 'SIZE' variables, plus a 2*SIZE multiply product array.  That's
** a total of 8*256k*2 (2 bytes per short)=4096k.  The FFTLen would
** of course be 1024k, so that'd be 1024k*8=8192k, plus the 4096k
** we need for the algorithm itself, for a total of 12meg.
**
** By reducing the size of the FFTLen, you could reduce the memory
** used so it would fit in an 8meg machine.  For a 4meg computer,
** you are going to have to use a virtual memory compiler if you
** actually want one million digits of pi.  On my computer, I can
** get 1 million digits in slightly under 3 hours.  But, using only
** 4 meg and virtual memory, and a reduced FFTLimit, it takes about
** 17 hours!  But, at least it is indeed possible.  Ohh, one last note,
** if you do use virtual memory (ie: disk), I'd really recommend not
** using any disk compression (stacker, doublespace, etc.) where the
** swap file is stored.  It wouldn't hurt anything, but it sure would
** be a waste of time.
*/

