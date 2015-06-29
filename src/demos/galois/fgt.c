// A lot of modint's on the stack.  Fix.
#include "pi.h"
#include "primes.h"
#include "fgt.h"
#include "modmath.h"

ModInt Prime;
CmplxModInt PrimvRoot,NthRoot,NthRoot1;
ModInt MulInv; /* easier as a single value */
int PrimeQ=0;

static CmplxModInt FGTPowRoots[50];

static void
FGTReOrder(CmplxModInt *Data, size_t Len)
{size_t Index,xednI,k;
 CmplxModInt Temp;
xednI=0;
for (Index=0;Index<Len;Index++)
  {
   if (xednI > Index)
     {
      Temp=Data[xednI];
      Data[xednI]=Data[Index];
      Data[Index]=Temp;
     }
   k=Len/2;
   while ((k <= xednI) && (k >=1)) {xednI-=k;k/=2;}
   xednI+=k;
  }
}

static void
FGT_T(CmplxModInt *Data, size_t Len)
/*
** A Radix 2, decimation in time, iterative FGT.
*/
{size_t j,Step,HalfStep,z;
 CmplxModInt u,w;

StartTimer(FFTITime);

Step=1;z=1;

while (Step < Len)
  {CmplxModInt *EndP=&Data[Len];
   HalfStep=Step;
   Step*=2;

   CmplxModSet(&u,1,0);
   w=FGTPowRoots[z++];
   {CmplxModInt *P1,*P2;
    P1=Data;P2=&Data[HalfStep];
    while (P1 < EndP)
      {
       CmplxModAddSub(P1,P2);
       P1+=Step;P2+=Step;
      }
   }
   for (j=1;j<HalfStep;j++)
     {CmplxModInt *P1,*P2;
      P1=&Data[j];P2=P1+HalfStep;
      CmplxModMul(&u,&u,&w);
      while (P1 < EndP)
        {
         CmplxModMul(P2,P2,&u);
         CmplxModAddSub(P1,P2);
         P1+=Step;P2+=Step;
        }
     }
  }
StopTimer(FFTITime);
}

static void
RFGT_T(CmplxModInt *Data,size_t Len)
/* Recursive decimation in time */
{size_t x,Len2;
 CmplxModInt P,Nth;
 CmplxModInt *Left,*Right;

Len2=Len/2;

if (Len<=(CPU_CACHE/sizeof(CmplxModInt))) {FGT_T(Data,Len);return;}
if (Len2 >= 2) {RFGT_T(Data,Len2);RFGT_T(Data+Len2,Len2);}

StartTimer(FFTRTime);
CmplxModSet(&P,1,0);
Nth=FGTPowRoots[Log2(Len)];Left=Data;Right=Data+Len2;
for (x=0;x<Len2;x++)
  {
   CmplxModMul(Right,Right,&P);
   CmplxModAddSub(Left,Right);
   Left++;Right++;
   CmplxModMul(&P,&P,&Nth);
  }
StopTimer(FFTRTime);
}

static void
DoFGT(CmplxModInt *Data,size_t Len, int Dir)
{size_t x=0,Power=Len;
 CmplxModInt Root;

StartTimer(FFTTime);

if (Dir > 0) Root=NthRoot;
else         Root=NthRoot1;
while (Power)
  {
   CmplxModIPow(&FGTPowRoots[x++],&Root,Power);
   Power/=2;
  }

FGTReOrder(Data,Len);
RFGT_T(Data,Len);
//FGT_T(Data,Len);

StopTimer(FFTTime);
}

static void
LoadNumBlockIntoFGT(size_t NX, INT32* NBuf, CmplxModInt *FGTData)
{
while (NX)
  {
//   NX-=RawCmplxModInts;
   NX-=2;
   CmplxModSet(FGTData,NBuf[NX+1],NBuf[NX]);
   FGTData++;
  }
}

void
DoFwdTransforms(CmplxModInt *FGTData,BigInt Num, size_t NumLen)
{size_t x,NX;
 size_t FGTLen=CalcFGTLen(NumLen);
 size_t NLen=NumLen;
 size_t FGTPos=0;
 INT32 *NBuf=(INT32*)FixedBuf;

  StartTimer(LoadTime);
  while (NLen)
    {
     NX=Min(FIXEDBUF_SIZE/sizeof(INT32),NLen);
     NLen-=NX;
     ReadNumIntoBuf(Num+NLen,NBuf,NX);
     LoadNumBlockIntoFGT(NX,NBuf,FGTData+FGTPos);
//     FGTPos+=(NX/RawCmplxModInts);
     FGTPos+=(NX/2);
    }

  for (x = FGTLen/2; x < FGTLen; x++)
    CmplxModSet(&FGTData[x],0,0);
  StopTimer(LoadTime);

  DoFGT(FGTData, FGTLen, 1);
}

void
DoRevTransforms(CmplxModInt *FGTData, size_t NumLen)
{size_t FGTLen=CalcFGTLen(NumLen);

 DoFGT(FGTData, FGTLen, -1);
}

void
SetFGTSize(size_t NumLen)
/*
** Setup for a specific size of FGT.
*/
{size_t FGTLen=CalcFGTLen(NumLen);
 ModInt Temp;

IModPow(&Temp,2,PrimeQ+1-Log2(FGTLen));
CmplxModPow(&NthRoot,&PrimvRoot,&Temp);
NthRoot1=NthRoot;CmplxConj(&NthRoot1);
IModPow(&MulInv,2,PrimeQ-Log2(FGTLen)-2); /* -2 is for the /4 the FGT needs */
}


void
InitFFT(size_t Size)
{ModInt Temp1,Temp2;


HexStr2Mod(&Prime,"1FFFFFFFFFFFFFFFFFFFFFF"); /* 2^89-1 */
PrimeQ=89;

IModPow(&Temp1,2,PrimeQ-2);
ModSet(&Temp2,2);
ModPow(&PrimvRoot.r,&Temp2,&Temp1);

IModPow(&Temp1,2,PrimeQ-2);
Temp2=Prime;ModSubInt(&Temp2,3);
ModPow(&PrimvRoot.i,&Temp2,&Temp1);
}

void
DeInitFFT(void)
{
}


