#include "pi.h"
#include "primes.h"
#include "ntt.h"
#include "modmath.h"

ModInt Prime,PrimvRoot,NthRoot,NthRoot1,MulInv;

static ModInt NTTPowRoots[50];
static int DoZeroPadCuts=0;

#define DO_ZERO_PAD_CUT

static ModInt
ModPow(ModInt Base,ModInt Expon)
{ModInt prod,b;

if (Expon<=0) return 1;

b=Base;
while (!(Expon&1)) {b=ModMul(b,b);Expon>>=1;}
prod=b;

while (Expon>>=1)
  {
   b=ModMul(b,b);
   if (Expon&1) prod=ModMul(prod,b);
  }
return prod;
}

static ModInt
FindInverse(ModInt Num)
{ModInt i;
i=ModPow(Num,Prime-2);
/*
** Num*3 can overflow causing the check to fail.
if (ModMul(Num*3,i) != 3)
  FatalError("Unable to find Mul inverse for %u mod %u\n",Num,Modulus);
*/
return i;
}

/* ========== End low level modular math stuff ========== */

#if 0
static void
NTTReOrder(ModInt *Data, size_t Len)
{size_t Index,xednI,k;
xednI=0;
for (Index=0;Index<Len;Index++)
  {
   if (xednI > Index)
     {ModInt Temp;
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
NTT_T(ModInt *Data, size_t Len)
/*
** A Radix 2, decimation in time, iterative NTT.
*/
{size_t j,Step,HalfStep,z;
 ModInt u,w;

StartTimer(FFTITime);

if (DoZeroPadCuts) {Step=2;z=2;}
else {Step=1;z=1;}

while (Step < Len)
  {ModInt *EndP=&Data[Len];
   HalfStep=Step;
   Step*=2;

   u=1;
   w=NTTPowRoots[z++];
   {ModInt *P1,*P2;
    P1=Data;P2=&Data[HalfStep];
    while (P1 < EndP)
      {ModInt D1=*P1,D2=*P2;
       *P2=ModSub(D1,D2);
       *P1=ModAdd(D1,D2);
       P1+=Step;P2+=Step;
      }
   }
   for (j=1;j<HalfStep;j++)
     {ModInt *P1,*P2;
      P1=&Data[j];P2=P1+HalfStep;
      u=ModMul(u,w);
      while (P1 < EndP)
        {ModInt D1,D2;
         D1=*P1;D2=*P2;
         D2=ModMul(D2,u);
         *P2=ModSub(D1,D2);
         *P1=ModAdd(D1,D2);
         P1+=Step;P2+=Step;
        }
     }
  }
StopTimer(FFTITime);
}

static void
RNTT_T(ModInt *Data,size_t Len)
/* Recursive decimation in time */
{size_t x,Len2;
 ModInt P=1,Nth;
 ModInt *Left,*Right;

Len2=Len/2;

if (Len<=(CPU_CACHE/sizeof(ModInt))) {NTT_T(Data,Len);return;}
if (Len2 >= 2) {RNTT_T(Data,Len2);RNTT_T(Data+Len2,Len2);}

StartTimer(FFTRTime);
Nth=NTTPowRoots[Log2(Len)];Left=Data;Right=Data+Len2;
for (x=0;x<Len2;x++)
  {ModInt b,c;
   b=Left[x];c=ModMul(Right[x],P);
   Left[x] = ModAdd(b,c);
   Right[x]= ModSub(b,c);
   P=ModMul(P,Nth);
  }
StopTimer(FFTRTime);
}
#endif

#if 1
static void
NTT_T(ModInt *Data, size_t Len)
/*
** A radix 2, decimation in time, iterative NTT.
*/
{size_t j,Step,HalfStep,z;
 ModInt u,w;
 ModInt *EndP=&Data[Len];

StartTimer(FFTITime);

Step=1;z=1;

while (Step < Len/2)
  {
   HalfStep=Step;
   Step*=2;

   u=1;
   w=NTTPowRoots[z++];

   {ModInt *P1,*P2;
    P1=Data;P2=&Data[HalfStep];
    while (P1 < EndP)
      {ModInt D1=*P1;
       ModInt D2=*P2;
       *P2=ModSub(D1,D2);
       *P1=ModAdd(D1,D2);
       P1+=Step;
       P2+=Step;
      }
    for (j=1;j<HalfStep;j++)
      {
       P1=&Data[j];P2=P1+HalfStep;
       u=ModMul(u,w);
       while (P1 < EndP)
         {ModInt D1,D2;
          D1=*P1;D2=*P2;
          D2=ModMul(D2,u);
          *P2=ModSub(D1,D2);
          *P1=ModAdd(D1,D2);
          P1+=Step;P2+=Step;
         }
      }
   }
  }
/* do the last pass to save some unneeded loop overhead. */
{ModInt *L,*R;
   HalfStep=Step;
   Step*=2;
   L=Data;R=Data+HalfStep;
   u=1;
   w=NTTPowRoots[z++];
    for (j=0;j<HalfStep;j++)
      {ModInt D1,D2;
       D1=*L;D2=*R;
       D2=ModMul(D2,u);
       *R=ModSub(D1,D2);
       *L=ModAdd(D1,D2);
       L++;R++;
       u=ModMul(u,w);
      }
 }
StopTimer(FFTITime);
}

static void
RNTT_T(ModInt *Data,size_t Len)
{size_t x,Len2;
 ModInt P=1,Nth;
 ModInt *Left,*Right;

Len2=Len/2;

if (Len<=(CPU_CACHE/sizeof(ModInt))) {NTT_T(Data,Len);return;}
if (Len2 >= 2) {RNTT_T(Data,Len2);RNTT_T(Data+Len2,Len2);}

StartTimer(FFTRTime);
Left=Data;Right=Data+Len2;
Nth=NTTPowRoots[Log2(Len)];

for (x=0;x<Len2;x+=4)
  {ModInt C0,C1,C2,C3;
   C0=ModMul(Right[0],P);P=ModMul(P,Nth);
   C1=ModMul(Right[1],P);P=ModMul(P,Nth);
   C2=ModMul(Right[2],P);P=ModMul(P,Nth);
   C3=ModMul(Right[3],P);P=ModMul(P,Nth);
   Right[0]=ModSub(Left[0],C0);
   Left[0] =ModAdd(Left[0],C0);
   Right[1]=ModSub(Left[1],C1);
   Left[1] =ModAdd(Left[1],C1);
   Right[2]=ModSub(Left[2],C2);
   Left[2] =ModAdd(Left[2],C2);
   Right[3]=ModSub(Left[3],C3);
   Left[3] =ModAdd(Left[3],C3);
   Left+=4;Right+=4;
  }

StopTimer(FFTRTime);
}
#endif

static void
NTT_F(ModInt *Data, size_t Len)
/* Iterative decimation in frequency */
{size_t j,Step,HalfStep;
 ModInt w,Nth;
 ModInt *EndP=&Data[Len];

StartTimer(FFTITime);
Step=Len;
while (Step>1)
  {
   HalfStep=Step/2;

   Nth=1;
   w=NTTPowRoots[Log2(Step)];
   {ModInt *L=Data,*R=Data+HalfStep;
    while (L < EndP)
      {ModInt u=*L,v=*R;
       *L=ModAdd(u,v);
       *R=ModSub(u,v);
       L+=Step;R+=Step;
      }
   }
   for (j=1;j<HalfStep;j++)
     {ModInt *L=&Data[j],*R=&Data[j+HalfStep];
      Nth=ModMul(Nth,w);
      while (L < EndP)
        {ModInt u=*L,v=*R;
         *L=ModAdd(u,v);
         *R=ModMul(ModSub(u,v),Nth);
         L+=Step;R+=Step;
        }
     }
   Step=HalfStep;
  }
StopTimer(FFTITime);
}

static void
RNTT_F(ModInt *Data,size_t Len)
/* Recursive decimation in frequency */
{size_t k,Len2;
 ModInt P=1;
 ModInt *Left,*Right;
 ModInt Nth;

StartTimer(FFTRTime);
Len2=Len/2;
Left=Data;Right=Data+Len2;
Nth=NTTPowRoots[Log2(Len)];
if (DoZeroPadCuts)
  for (k=0;k<Len2;k++)
    {
     Right[k]=ModMul(Left[k],P);
     P=ModMul(P,Nth);
    }
else
  for (k=0;k<Len2;k++)
    {ModInt x,y;
     x=Left[k];y=Right[k];
     Left[k] =ModAdd(x,y);
     Right[k]=ModMul(ModSub(x,y),P);
     P=ModMul(P,Nth);
    }
DoZeroPadCuts=0;
StopTimer(FFTRTime);

if (Len2<=(CPU_CACHE/sizeof(ModInt)))
  {NTT_F(Data,Len2);NTT_F(Data+Len2,Len2);return;}
RNTT_F(Data,Len2);
RNTT_F(Data+Len2,Len2);
}


static void
DoNTT(ModInt *Data,size_t Len, int Dir)
{size_t x=0,Power=Len;
 ModInt Root;

StartTimer(FFTTime);
DoZeroPadCuts=0;

if (Dir > 0) Root=NthRoot;
else         Root=NthRoot1;
while (Power)
  {
   NTTPowRoots[x++]=ModPow(Root,Power);
   Power/=2;
  }

#ifdef DO_ZERO_PAD_CUT
if (Dir > 0) DoZeroPadCuts=1;
#endif

/*NTTReOrder(Data,Len);RNTT_T(Data,Len);*/
if (Dir > 0) RNTT_F(Data,Len);
else         RNTT_T(Data,Len);

StopTimer(FFTTime);
}

static void
LoadNumBlockIntoNTT(size_t NX, INT32* NBuf, ModInt *NTTData)
{size_t x;
/*
** HARDWIRED
*/
for (x = 0; x < NX; x++)
  {
   NTTData[x*2]  =NBuf[NX - x - 1] % 10000;
   NTTData[x*2+1]=NBuf[NX - x - 1] / 10000;
  }
}

void
DoFwdTransforms(ModInt *NTTData,BigInt Num, size_t NumLen)
{size_t x,NX;
 size_t NTTLen=CalcNTTLen(NumLen);
 size_t NLen=NumLen;
 size_t NTTPos=0;
 INT32 *NBuf=(INT32*)FixedBuf;

  StartTimer(LoadTime);
  while (NLen)
    {
     NX=Min(FIXEDBUF_SIZE/sizeof(INT32),NLen);
     NLen-=NX;
     ReadNumIntoBuf(Num+NLen,NBuf,NX);
     LoadNumBlockIntoNTT(NX,NBuf,NTTData+NTTPos);
     NTTPos+=(NX*2); // HARDWIRED
    }

#ifdef DO_ZERO_PAD_CUT
  for (x = NTTLen/2; x < NTTLen; x++) NTTData[x] = NTTData[x-NTTLen/2];
#else
  for (x = NTTLen/2; x < NTTLen; x++) NTTData[x] = 0;
#endif
  StopTimer(LoadTime);

  DoNTT(NTTData, NTTLen, 1);
}

void
DoRevTransforms(ModInt *NTTData, size_t NumLen)
{size_t NTTLen=CalcNTTLen(NumLen);

 DoNTT(NTTData, NTTLen, -1);
}

void
SetNTTSize(size_t NumLen)
/*
** Setup for a specific size of NTT.
*/
{size_t NTTLen=CalcNTTLen(NumLen);

 NthRoot  =ModPow(PrimvRoot,Prime-1-(Prime-1)/NTTLen);
 NthRoot1 =ModPow(PrimvRoot,(Prime-1)/NTTLen);
 MulInv   =FindInverse(NTTLen);
}


void
InitFFT(size_t Size)
{
Prime    =PrimeList[0];
PrimvRoot=PrimvRootList[0];
}

void
DeInitFFT(void)
{
}


