#include "pi.h"
#include "block.h"
/*
** This file contains the 'core' math routines.  They only work
** with in-memory blocks of data.
**
** It doesn't know about BigInt, disk, etc. etc.
*/

void
BlockClear(INT32 *Beg,INT32 *End)
{
if (Beg > End) while (End!=Beg) *End++=0;
else           while (Beg!=End) *Beg++=0;
}

void
BlockCopy(INT32 *Dest,INT32 *Src,size_t Len)
{
  memmove(Dest,Src,Len*sizeof(INT32));
}

void
BlockPack2(INT32 *Data,size_t Len)
/*
** Pack a sequence of bytes holding 2 digits into a
** 'Len' sequence of INT32 words holding 8 digits.
*/
{unsigned char *Byte=(unsigned char*)Data;
  while (Len)
   {INT32 D;
    D=Byte[0];
    D=(D*100)+Byte[1];
    D=(D*100)+Byte[2];
    D=(D*100)+Byte[3];
    *Data=D;
    Byte+=4;
    Data+=1;
    Len--;
   }
}

void
BlockPack4(INT32 *Data,size_t Len)
/*
** Pack a sequence of short ints holding 4 digits into a
** 'Len' sequence of INT32 words holding 8 digits.
*/
{unsigned short int *Word=(unsigned short int*)Data;
  while (Len)
   {INT32 D;
    D=Word[0];
    D=(D*10000)+Word[1];
    *Data=D;
    Word+=2;
    Data+=1;
    Len--;
   }
}

void
BlockUnpack2(INT32 *Data,size_t Len)
/*
** UnPack a sequence of 'Len' INT32 words holding 8 digits.
** into a sequence of bytes holding 2 digits.
*/
{unsigned char *Byte=(unsigned char*)Data;
  while (Len)
   {INT32 D=*Data;
    INT32 A,B;
    A=(D / 10000);
    B=(D % 10000);
    Byte[0]=(unsigned char)(A / 100);
    Byte[1]=(unsigned char)(A % 100);
    Byte[2]=(unsigned char)(B / 100);
    Byte[3]=(unsigned char)(B % 100);
/*
    Byte[0]=(D / 1000000) % 100;
    Byte[1]=(D / 10000) % 100;
    Byte[2]=(D / 100) % 100;
    Byte[3]=(D) % 100;
*/
    Byte+=4;
    Data+=1;
    Len--;
   }
}

void
BlockUnpack4(INT32 *Data,size_t Len)
/*
** UnPack a sequence of 'Len' INT32 words holding 8 digits.
** into a sequence of short ints holding 4 digits.
*/
{unsigned short int *Word=(unsigned short int*)Data;
  while (Len)
   {INT32 D=*Data;
    Word[0]=(unsigned short int)(D / 10000);
    Word[1]=(unsigned short int)(D % 10000);
    Word+=2;
    Data+=1;
    Len--;
   }
}

size_t
BlockCountZeros(INT32 *Num,size_t Len)
{size_t x,count;
  count=0;
  for (x=0;x < Len;x++) if (Num[x] == 0) count++;
  return count;
}

INT32
BlockDivBy(INT32 *Result, INT32 *Num,INT32 Val,INT32 Remain,size_t Len)
{size_t x;
 INT32 Temp=(INT32)Remain;

if (Val > 10) FatalError("BlockDivBy called with %d\n",Val);

  if (Val==2)
    for (x = 0; x < Len; x++)
      {
        Temp = ((INT32)Num[x]) + Temp * 100000000;
        Result[x] =(INT32) (Temp / 2);
        Temp      = Temp % 2;
      }
  else if (Val==4)
    for (x = 0; x < Len; x++)
      {
        Temp = ((INT32)Num[x]) + Temp * 100000000;
        Result[x] =(INT32) (Temp / 4);
        Temp      = Temp % 4;
      }
  else
    for (x = 0; x < Len; x++)
      {
        Temp = ((INT32)Num[x]) + Temp * 100000000;
        Result[x] =(INT32) (Temp / Val);
        Temp      = Temp % Val;
      }
return (INT32)Temp;
}

double
BlockMulByFloat(INT32 *Num, double Val, double Carry,size_t Len)
{double Temp;

 while (Len)
   {
    Len--;
    Temp     = Num[Len] * Val + Carry;
    Carry    = floor(Temp /   100000000.0);
    Num[Len] = Temp - Carry * 100000000.0;
   }
return Carry;
}

INT32
BlockMulBy(INT32 *Result,INT32 *Num, INT32 Val, INT32 Carry, size_t Len)
{INT32 Prod=Carry;

  if (Val > 10) /* too big to fit into 32 bit integer */
    {
     if (Result != Num) BlockCopy(Result,Num,Len);
     return (INT32)BlockMulByFloat(Result,Val,Carry,Len);
    }

  while (Len)
    {
     Len--;
     Prod = (INT32)Num[Len] * Val + Prod;
     Result[Len] = (INT32) (Prod % 100000000);
     Prod        = (Prod / 100000000);
    }

return (INT32)Prod;
}

INT32
BlockAdd(INT32 *Sum, INT32 *Num1, INT32 *Num2, INT32 Carry, size_t Len)
{
  while (Len)
    {
     Len--;
     Carry    = (INT32)(Num1[Len] + Num2[Len] + Carry);
     Sum[Len] = (INT32)(Carry % 100000000);
     Carry    = (INT32)(Carry / 100000000);
    }
  return Carry;
}

INT32
BlockSub(INT32 *Dif, INT32 *Num1, INT32 *Num2, INT32 Borrow,size_t Len)
{INT32 temp;

  while (Len)
    {
     Len--;
     temp = (INT32)(Num1[Len] - Num2[Len] - Borrow);
     Borrow = 0;
     if (temp < 0) { Borrow = 1; temp = (INT32)(temp+100000000); }
     Dif[Len] = temp;
    }
  return Borrow;
}

INT32
BlockNegate(INT32 *Num, INT32 Borrow,size_t Len,INT32 Val0)
{INT32 Temp;

  while (Len > 1)
    {
     Len--;
     Temp = (INT32)(0 - Num[Len] - Borrow);
     Borrow = 0;
     if (Temp < 0) {Borrow = 1; Temp =(INT32)(Temp+100000000);}
     Num[Len] = Temp;
    }
/*
** Do the first digit seperately.  Allows us to do more than simple
** negation.  We can do Num = x.0000 - Num
*/
  Temp = (INT32)(Val0 - Num[0] - Borrow);
  Borrow = 0;
  if (Temp < 0) {Borrow = 1; Temp =(INT32)(Temp+100000000);}
  Num[0] = Temp;

  return Borrow;
}


INT32
BlockSlowMul(INT32 *Prod, INT32 *Num1, INT32 *Num2, size_t Len)
/*
** Simple, slow schoolboy multiplication.
** We don't have any way to multiple two 1e8 numbers and get
** a 1e16 result, so we have to break the numbers down into
** Len*2 1e4 numbers, and then later recombine them into 1e8.
**
** The packing & unpacking isn't the most efficient, but hey,
** this whole mul is inefficient anyway!  This routine is only
** going to be used for small multiplications, so none of it
** matters as long as it gets the right answer.
*/
{size_t Ndx1, Ndx2, PNdx;
 UINT32 Carry=0;
 unsigned short int *P,*N1,*N2;

  if ((Num1==Prod) || (Num2==Prod))
    FatalError("BlockSlowMul can't have the output as one of the inputs.\n");

  BlockUnpack4(Prod,Len*2);
  BlockUnpack4(Num1,Len);
  if (Num1 != Num2) BlockUnpack4(Num2,Len);
  Len*=2;
  P=(unsigned short int *)Prod;
  N1=(unsigned short int *)Num1;
  N2=(unsigned short int *)Num2;

  Ndx2=Len;
  while (Ndx2)
    {
      Ndx2--;
      PNdx = Ndx2+Len;
      if (N2[Ndx2] != 0)
        {
         Ndx1=Len;
         while (Ndx1)
           {
             Ndx1--;
             Carry   = ((UINT32) N1[Ndx1]) * ((UINT32) N2[Ndx2]) +
                       Carry + ((UINT32) P[PNdx]);
             P[PNdx] = (unsigned short int)(Carry % 10000);
             Carry   =          Carry / 10000;
             PNdx--;
           }
        }
      P[Ndx2]+=(unsigned short int)Carry;
      Carry=0;
    }
Carry=P[0] / 10000;
P[0] = (unsigned short int)(P[0] % 10000);

Len/=2;
BlockPack4(Prod,Len*2);
BlockPack4(Num1,Len);
if (Num1 != Num2) BlockPack4(Num2,Len);

return (INT32)Carry;
}


