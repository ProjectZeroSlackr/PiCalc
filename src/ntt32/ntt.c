#include "pi.h"
#include "primes.h"
#include "ntt.h"
#include "modmath.h"
#include "vector.h"


ModInt Prime;
double RecipPrime=0.0;
ModInt Inverses[NPrimes][NPrimes];

#if defined(VECTOR_NTT)
ModInt* TrigTable;
ModInt* TrigTables[2][NPrimes];
#endif

ConstList Consts[NPrimes]; /* This should be static, but I need to cheat for Knuth. */
ModInt NTTPowRoots[50];

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


#if defined(VECTOR_NTT)
/*
** My trig table generator is faster than Jason's.  Fast enough
** I can just use it 'on the fly'.
*/
static void
GenTrigTable_T(size_t Len)
{size_t j,Step,HalfStep,z;
 UINT32 u,w;
 size_t TNdx=0;

Step=2;z=2;
while (Step < Len)
  {
   HalfStep=Step;Step*=2;

   u=1; w=NTTPowRoots[z++];
   for (j=0;j<HalfStep;j++) {TrigTable[TNdx++]=u;u=ModMul(u,w);}
  }
}

void
GenTrigTable_F(size_t Len)
/* non-recursive decimation in frequency */
{size_t j,step,halfstep;
 ModInt w,Nth;
 size_t TNdx=0;

step=Len;
while (step)
  {
   halfstep=step/2;

   Nth=1; w=NTTPowRoots[Log2(step)];
   for (j=0;j<halfstep;j++) {TrigTable[TNdx++]=Nth;Nth=ModMul(Nth,w);}
   step=halfstep;
  }
}

static void
NTT_T( ModInt *Data, size_t Len )
/*
** Do a table based, vector style Decimation in Time NTT
**
** This routine was originally written by Jason P. and is public domain.
*/
{
 ModInt *Left, *Right, *roots=TrigTable+2;
 size_t blocks, size;
 size_t Half=Len/2;

StartTimer(FFTITime);

PrepVector(Prime,Len);
VectorNTT_First2( Data, TrigTable[1], Len );

for (size=4;size < Len;size*=2)
  {
   blocks = Half/size;
   Left  = Data;
   Right = Data + size;

   while ( blocks )
     {
      VectorModMul( Right, roots, size );
      VectorModButterfly( Left, Right, size );
      Left  += 2*size;
      Right += 2*size;
      blocks--;
     }

   roots += size;
  }

StopTimer(FFTITime);
}

static void
NTT_F( ModInt *Data, size_t Len )
/*
** Do a table based, vector style Decimation in Frequency NTT
**
** This routine was based on one originally written by
** Jason P. and is public domain.
*/
{
 ModInt *Left, *Right, *roots;
 size_t blocks, size;
 size_t Half=Len/2;

StartTimer(FFTITime);

roots=TrigTable + ((CPU_CACHE/sizeof(ModInt)) - Len);
PrepVector(Prime,Len);
for (size=Len/2;size > 2;size/=2)
  {
   blocks = Half/size;
   Left  = Data;
   Right = Data + size;

   while ( blocks )
     {
      VectorModButterfly( Left, Right, size );
      VectorModMul( Right, roots, size );
      Left  += 2*size;
      Right += 2*size;
      blocks--;
     }

   roots += size;
  }

VectorNTT_Last2( Data, roots[1], Len );

StopTimer(FFTITime);
}

#else
/* Not vector iterative NTTs */

static void
NTT_T(ModInt *Data, size_t Len)
/*
** A radix 2, decimation in time, iterative NTT.
*/
{size_t j,Step,HalfStep,z;
 ModInt u,w;
 ModInt *EndP=&Data[Len];

StartTimer(FFTITime);

Step=1;
z=1;

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
   L=Data;R=Data+HalfStep;
   u=1;
   w=NTTPowRoots[z];
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

/* End of vector & non-vector iterative NTTs */
#endif



#ifdef VECTOR_RNTT
ModInt PList[1024];
#endif

static void
RNTT_T(ModInt *Data,size_t Len)
{size_t x,Len2;
 ModInt P=1,Nth;
 ModInt *Left,*Right;
#ifdef VECTOR_RNTT
 size_t BlockLen,y;
#endif

Len2=Len/2;

if (Len<=(CPU_CACHE/sizeof(ModInt))) {NTT_T(Data,Len);return;}
if (Len2 >= 2) {RNTT_T(Data,Len2);RNTT_T(Data+Len2,Len2);}

StartTimer(FFTRTime);
Left=Data;Right=Data+Len2;
Nth=NTTPowRoots[Log2(Len)];

#ifdef VECTOR_RNTT
BlockLen=Min(1024,(CPU_CACHE/sizeof(ModInt)));
if (BlockLen*VECTOR_RNTT <= Len2)
  {/* Only efficient if we can do at least 4 */ /* Sub-optimal???? */
   for (y=0;y<BlockLen;y++) {PList[y]=P;P=ModMul(P,Nth);}

   PrepVector(Prime,BlockLen);
   for (x=0;x<Len2;x+=BlockLen)
     {
      VectorModMul( Right, PList,BlockLen );
      VectorModButterfly( Left, Right, BlockLen );
      if (x+BlockLen != Len2)
        VectorModMulC( PList, P ,BlockLen );
      Left+=BlockLen;Right+=BlockLen;
     }
  }
else
#endif
  {
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
  }

StopTimer(FFTRTime);
}

static void
RNTT_F(ModInt *Data,size_t Len)
/* Recursive decimation in frequency */
{size_t k,Len2;
 ModInt P=1;
 ModInt *Left,*Right;
 ModInt Nth;
#ifdef VECTOR_RNTT
 size_t BlockLen,x,y;
#endif

if (Len<=(CPU_CACHE/sizeof(ModInt))) {NTT_F(Data,Len);return;}

StartTimer(FFTRTime);
Len2=Len/2;
Left=Data;Right=Data+Len2;
Nth=NTTPowRoots[Log2(Len)];
#ifdef VECTOR_RNTT
BlockLen=Min(1024,(CPU_CACHE/sizeof(ModInt)));
if (BlockLen*VECTOR_RNTT <= Len2)
  {/* Only efficient if we can do at least 4 */ /* Sub-optimal???? */
   for (y=0;y<BlockLen;y++) {PList[y]=P;P=ModMul(P,Nth);}

   PrepVector(Prime,BlockLen);
   for (x=0;x<Len2;x+=BlockLen)
     {
      VectorModButterfly( Left, Right, BlockLen );
      VectorModMul( Right, PList,BlockLen );
      if (x+BlockLen != Len2)
        VectorModMulC( PList, P ,BlockLen );
      Left+=BlockLen;Right+=BlockLen;
     }
  }
else
#endif
  {
   for (k=0;k<Len2;k++)
     {ModInt x,y;
      x=Left[k];y=Right[k];
      Left[k] =ModAdd(x,y);
      Right[k]=ModMul(ModSub(x,y),P);
      P=ModMul(P,Nth);
     }
  }

StopTimer(FFTRTime);

RNTT_F(Data,Len2);
RNTT_F(Data+Len2,Len2);
}


static void
DoNTT(ModInt *Data,size_t Len, int Dir, size_t Pass)
{size_t x=0,Power=Len;
 ModInt Root;

SetModPrime(Pass);
StartTimer(FFTTime);

if (Dir > 0) Root=Consts[Pass].NthRoot;
else         Root=Consts[Pass].NthRoot1;
while (Power)
  {
   NTTPowRoots[x++]=ModPow(Root,Power);
   Power/=2;
  }

#if defined(VECTOR_NTT)
 if (Dir > 0) TrigTable=TrigTables[1][Pass];
 else         TrigTable=TrigTables[0][Pass];
#endif

if (Dir > 0) RNTT_F(Data,Len);
else         RNTT_T(Data,Len);

StopTimer(FFTTime);
}

static void
LoadNumBlockIntoNTT(size_t NX,INT32* NBuf,size_t Pass, ModInt *NTTData)
{

SetModPrime(Pass);

while (NX)
  {ModInt R,D,Z;
   NX-=RawModInts;
   R=((ModInt)NBuf[NX+0]);
   Z=((ModInt)NBuf[NX+1]);D=ModMul(R,BaseVar);R=ModAdd(Z,D);
#if RawModInts==4
   Z=((ModInt)NBuf[NX+2]);D=ModMul(R,BaseVar);R=ModAdd(Z,D);
   Z=((ModInt)NBuf[NX+3]);D=ModMul(R,BaseVar);R=ModAdd(Z,D);
#endif
   *NTTData=R;NTTData++;
  }
}

void
DoFwdTransforms(ModInt *NTTData,BigInt Num, size_t NumLen,size_t StartP,size_t EndP)
{size_t Pass,x,NX;
 size_t NTTLen=CalcNTTLen(NumLen);
 ModInt *Data;
 size_t NLen=NumLen;
 size_t NTTPos=0;
 INT32 *NBuf=(INT32*)FixedBuf;

  StartTimer(LoadTime);
  while (NLen)
    {
     NX=Min(FIXEDBUF_SIZE/sizeof(INT32),NLen);
     NLen-=NX;
     ReadNumIntoBuf(Num+NLen,NBuf,NX);
     for (Pass=StartP;Pass < EndP;Pass++)
       LoadNumBlockIntoNTT(NX,NBuf,Pass,NTTData+NTTLen*(Pass-StartP) +NTTPos);
     NTTPos+=(NX/RawModInts);
    }

  for (Pass=StartP;Pass < EndP;Pass++)
    {
     Data=NTTData+NTTLen*(Pass-StartP);
     for (x = NTTLen/2; x < NTTLen; x++) Data[x] = 0;
    }
  StopTimer(LoadTime);

  for (Pass=StartP;Pass < EndP; Pass++)
     DoNTT(NTTData+NTTLen*(Pass-StartP), NTTLen, 1, Pass);
}

void
DoRevTransforms(ModInt *NTTData, size_t NumLen, size_t StartP, size_t EndP)
{size_t Pass;
 size_t NTTLen=CalcNTTLen(NumLen);

  for (Pass=StartP; Pass < EndP; Pass++)
     DoNTT(NTTData+NTTLen*(Pass-StartP), NTTLen, -1, Pass);
}

void
SetNTTSize(size_t NumLen)
/*
** Setup for a specific size of NTT.
*/
{size_t x,Pass;
 size_t NTTLen=CalcNTTLen(NumLen);

  for (x=0;x<NPrimes;x++)
    {ModInt PrimvRoot;
     Prime    =Consts[x].Prime=PrimeList[x];
#ifdef NEED_RECIP_PRIME
     RecipPrime=Consts[x].RecipPrime= 1.0 / Prime;
#endif
     PrimvRoot=Consts[x].PrimvRoot=PrimvRootList[x];
               Consts[x].NthRoot=ModPow(PrimvRoot,Prime-1-(Prime-1)/NTTLen);
               Consts[x].NthRoot1=ModPow(PrimvRoot,(Prime-1)/NTTLen);
               Consts[x].MulInv=FindInverse(NTTLen);
    }

  for (Pass=0;Pass < NPrimes; Pass++)
     {
      SetModPrime(Pass);
      for (x=0;x < Pass;x++)
        Inverses[x][Pass]=FindInverse(PrimeList[x]);
     }
}

void
InitFFT(size_t Size)
#if defined(VECTOR_NTT)
{int d,p,x;size_t Power;ModInt Root;

InitVectorModMath(Size);

for (d=0;d<2;d++)
  for (p=0;p<NPrimes;p++)
    {
     Power=Min(CalcNTTLen(Size),CPU_CACHE/sizeof(ModInt));
     Power=CPU_CACHE/sizeof(ModInt);
     TrigTables[d][p]=(ModInt*)AlignedMalloc(Power*sizeof(ModInt));
     if (TrigTables[d][p]==NULL)
       FatalError("Unable to allocate eough memory for the trig tables. %d %d\n",d,p);
     SetNTTSize(CalcNumLen(Power));
     SetModPrime(p);
     TrigTable=TrigTables[d][p];

     if (d > 0) Root=Consts[p].NthRoot;
     else       Root=Consts[p].NthRoot1;
     x=0;
     while (Power)
       {
        NTTPowRoots[x++]=ModPow(Root,Power);
        Power/=2;
       }
     if (d > 0) GenTrigTable_F(CPU_CACHE/sizeof(ModInt));
     else       GenTrigTable_T(CPU_CACHE/sizeof(ModInt));
    }
}
#else
{
InitVectorModMath(Size);
}
#endif

void
DeInitFFT(void)
#if defined(VECTOR_NTT)
{int d,p;
for (d=0;d<2;d++)
  for (p=0;p<NPrimes;p++)
     AlignedFree(TrigTables[d][p]);
DeInitVectorModMath();
}
#else
{
DeInitVectorModMath();
}
#endif

