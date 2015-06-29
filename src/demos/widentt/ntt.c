// A lot of modint's on the stack.  Fix.
#include "pi.h"
#include "primes.h"
#include "ntt.h"
#include "modmath.h"

ModInt Prime,PrimvRoot,NthRoot,NthRoot1,MulInv;

static ModInt NTTPowRoots[50];
static int DoZeroPadCuts=0;

#define DO_ZERO_PAD_CUT

#if 0
static void
NTTReOrder(ModInt *Data, size_t Len)
{size_t Index,xednI,k;
 ModInt Temp;
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
#endif

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

   ModSet1(&u);
   w=NTTPowRoots[z++];
   {ModInt *P1,*P2;
    P1=Data;P2=&Data[HalfStep];
    while (P1 < EndP)
      {
       ModAddSub(P1,P2);
       P1+=Step;P2+=Step;
      }
   }
   for (j=1;j<HalfStep;j++)
     {ModInt *P1,*P2;
      P1=&Data[j];P2=P1+HalfStep;
      ModMul(&u,&u,&w);
      while (P1 < EndP)
        {
         ModMul(P2,P2,&u);
         ModAddSub(P1,P2);
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
 ModInt P,Nth;
 ModInt *Left,*Right;

Len2=Len/2;

if (Len<=(CPU_CACHE/sizeof(ModInt))) {NTT_T(Data,Len);return;}
if (Len2 >= 2) {RNTT_T(Data,Len2);RNTT_T(Data+Len2,Len2);}

StartTimer(FFTRTime);
ModSet1(&P);
Nth=NTTPowRoots[Log2(Len)];Left=Data;Right=Data+Len2;
for (x=0;x<Len2;x++)
  {
   ModMul(Right,Right,&P);
   ModAddSub(Left,Right);
   Left++;Right++;
   ModMul(&P,&P,&Nth);
  }
StopTimer(FFTRTime);
}

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

   ModSet1(&Nth);
   w=NTTPowRoots[Log2(Step)];
   {ModInt *L=Data,*R=Data+HalfStep;
    while (L < EndP)
      {
       ModAddSub(L,R);
       L+=Step;R+=Step;
      }
   }
   for (j=1;j<HalfStep;j++)
     {ModInt *L=&Data[j],*R=&Data[j+HalfStep];
      ModMul(&Nth,&Nth,&w);
      while (L < EndP)
        {
         ModAddSub(L,R);
         ModMul(R,R,&Nth);
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
 ModInt *Left,*Right;
 ModInt P,Nth;

StartTimer(FFTRTime);
Len2=Len/2;
Left=Data;Right=Data+Len2;
Nth=NTTPowRoots[Log2(Len)];
ModSet1(&P);
if (DoZeroPadCuts)
  for (k=0;k<Len2;k++)
    {
     ModMul(Right,Left,&P);
     ModMul(&P,&P,&Nth);
     Left++;Right++;
    }
else
  for (k=0;k<Len2;k++)
    {
     ModAddSub(Left,Right);
     ModMul(Right,Right,&P);
     Left++;Right++;
     ModMul(&P,&P,&Nth);
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
   ModIPow(&NTTPowRoots[x++],&Root,Power);
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
{size_t z;
while (NX)
  {
   NX-=RawModInts;
   ModClear(NTTData);
   for (z=0;z<RawModInts;z++)
     {
      CB_MulInt(NTTData,100000000);
      CB_AddInt(NTTData,NBuf[NX+z]);
     }
   NTTData++;
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
     NTTPos+=(NX/RawModInts);
    }

  for (x = NTTLen/2; x < NTTLen; x++)
#ifdef DO_ZERO_PAD_CUT
    NTTData[x]=NTTData[x-NTTLen/2];
#else
    ModClear(&NTTData[x]);
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
 ModInt Temp;

  ModSet(&Temp,NTTLen);NthRoot=Prime;ModSub1(&NthRoot);ModSub1(&NthRoot);
 ModPow(&MulInv,&Temp,&NthRoot);

  Temp=Prime;ModSub1(&Temp);ModDivInt(&Temp,NTTLen);
  ModSub(&Temp,&Prime,&Temp);ModSub1(&Temp);
 ModPow(&NthRoot, &PrimvRoot,&Temp);

  Temp=Prime;ModSub1(&Temp);ModDivInt(&Temp,NTTLen);
 ModPow(&NthRoot1, &PrimvRoot,&Temp);

// ModPow(NthRoot, PrimvRoot,Prime-1-(Prime-1)/NTTLen);
// ModPow(Nthroot1,PrimvRoot,(Prime-1)/NTTLen);
// MulInv   =FindInverse(NTTLen);
}


void
InitFFT(size_t Size)
{
#if 0
HexStr2Mod(&Prime,"CF00000000000000000000000000000000000000000000000000000000000001"); /* 207*2^248+1 */
HexStr2Mod(&PrimvRoot,"5"); /* primitive root=5 */
#endif

#if 1
HexStr2Mod(&Prime,"4083000000000000000000000000000000000000000000000000000000000001"); /* 16515*2^240+1 */
HexStr2Mod(&PrimvRoot,"d"); /* primitive root=13 */
#endif

#if 0
HexStr2Mod(&Prime,"1FE00000000000000000000000000000000000000000000000000000000001"); /* 255*2^237+1 */
HexStr2Mod(&PrimvRoot,"13"); /* primitive root=19 ($13) */
#endif

#if 0
HexStr2Mod(&Prime,"ffffffff00000001"); // 18446742974197923841
HexStr2Mod(&PrimvRoot,"7");
#endif
}

void
DeInitFFT(void)
{
}


