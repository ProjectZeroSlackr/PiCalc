#include "pi.h"
#include "primes.h"
#include "ntt.h"
#include "modmath.h"
#include "vector.h"


ModInt Prime;
double RecipPrime=0.0;
ModInt Inverses[NPrimes][NPrimes];
float chopper64 = 9223372036854775808.0;        // 2^63
double dmodulus, imodulus;
int modulus;

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

#include "sixstep.c"

static void
DoNTT(ModInt *Data,size_t Len, int Dir, size_t Pass)
{size_t x=0,Power=Len;
 ModInt Root;
 int OldFPU=_control87(0,0);

SetModPrime(Pass);
StartTimer(FFTTime);

if (Dir > 0) Root=Consts[Pass].NthRoot;
else         Root=Consts[Pass].NthRoot1;
while (Power)
  {
   NTTPowRoots[x++]=ModPow(Root,Power);
   Power/=2;
  }

/* Only the APFloat six-step needs these vars, so we can set them here. */
 modulus   = Prime;
 dmodulus  = Prime;
 imodulus  = 1.0/Prime;

_control87(RC_CHOP,MCW_RC);
if (Dir > 0) FSixStep(Data,Consts[Pass].PrimvRoot,Root, 1,Len);
else         RSixStep(Data,Consts[Pass].PrimvRoot,Root,-1,Len);
_control87(OldFPU,MCW_RC);

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
#if defined(VECTOR_NTT) || defined(TRIG_NTT)
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
     GenTrigTable(CPU_CACHE/sizeof(ModInt));
    }
}
#else
{
InitVectorModMath(Size);
}
#endif

void
DeInitFFT(void)
#if defined(VECTOR_NTT) || defined(TRIG_NTT)
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

