#include "pi.h"
#include "ntt.h"
#include "modmath.h"

/*
** In the spirit of Jason P.'s ntt586 (which these routines were derived
** from), this file is also placed into the Public Domain.
*/

/*
** This routine is C code that targets the Pentium processor.
**
** Under DJGPP 2.7.2.1, it mostly generates very good code and the
** total program speed (including C versions of crt.h and modmath.h)
** is actually only about 60% slower than the asm.
**
** Some routines could use some improvement.
**
** How your compiler will deal with this code is the biggest
** unknown.  And I can't help you here.  On the other hand, since
** this is just simple C, you can tweak it as needed.
*/

/*
** These routines _require_ that the FPU be in full 'long double'
** mode and that it has rounding enabled.  With many PC compilers,
** you use the function _control87(), but it does vary, as does
** which header that function is defined in.  If you need to
** set it, you could do so in "InitVector()", but if you _do_ need
** to set it, the compiler might not appreciate it since it is
** probably depending upon it to be a certain way.
*/

/*
** GCC v2.8.1 has absolutely terrible data alignment problems.
** What's worse, the method GCC tells you (in GCC.I9 doc file)
** to use to fix alignment problems, doesn't work!
**
** This means you are going to have to fix this yourself!  *AND*
** this will change depending upon whether you are using disk numbers
** or virtual memory, or the Pentium timers.  Just move the alignment
** variable out of the comments.  You may need to adjust the number
** of elements in the array, too.

** The BigInt is there since that data type changes size.  This allows
** a resonably consistant alignment.  The int array is there in case
** you need to align more.
**
** **HOWEVER** there is no guarantee that the assembler or linker
** will actually put the alignment variable where it will do any
** good!!  There is absolutely no way around this.  Sorry.  I really
** wish there was a better generic solution.
*/
BigInt VectorAlign1;
int VectorAlign2[1];
ModInt JPrime;
double VectStuff[9];
union {double D;long int i[2];} ChopShop[4];
/*
** 0-3 For chopping double's to int
** 4   MAGIC1=3.0*2147483648.0*2147483648.0;
** 5   MAGIC2=3.0*2147483648.0*2147483648.0;
** 6   JUSTIFY=6755399441055744.0;
** 7   Prime
** 8   Reciprocal prime.
*/
typedef long double LDouble;

/*
** Since the x87 has a terrible round to integer, we have to cheat.
** While we are at it, also take care of the possible underflow.
*/
#define Chop(Dest,z)                 \
  {register long int ALU;            \
   ALU = ChopShop[z].i[0];           \
   ALU = ALU + ((ALU>>31) & JPrime); \
   Dest = ALU;                       \
  }

void
VectorModAdd(ModInt *Num1, ModInt *Num2, size_t Len)
{
if ((Len < 4) || (Len & 3))
  FatalError("VectorModAdd called with incorrect length: %lu\n",(ULINT)Len);

while (Len)
  {
   Num1[0]=ModAdd(Num1[0],Num2[0]);
   Num1[1]=ModAdd(Num1[1],Num2[1]);
   Num1[2]=ModAdd(Num1[2],Num2[2]);
   Num1[3]=ModAdd(Num1[3],Num2[3]);
   Num1+=4;Num2+=4;Len-=4;
  }
}

void
VectorModSub(ModInt *Num1, ModInt *Num2, size_t Len)
{
if ((Len < 4) || (Len & 3))
  FatalError("VectorModSub called with incorrect length: %lu\n",(ULINT)Len);

while (Len)
  {
   Num1[0]=ModSub(Num1[0],Num2[0]);
   Num1[1]=ModSub(Num1[1],Num2[1]);
   Num1[2]=ModSub(Num1[2],Num2[2]);
   Num1[3]=ModSub(Num1[3],Num2[3]);
   Num1+=4;Num2+=4;Len-=4;
  }
}

void
VectorModMul( ModInt *a, ModInt *b, size_t size )
/*  Num1[0...len] <- Num1[0...len] * Num2[0...len] mod n */
{
   do {
       register LDouble f0a, f0b;
       register LDouble f1a, f1b;

      f0a = (double)a[0];
      f0b = (double)b[0];
      f1a = (double)a[1];
      f1b = (double)b[1];
      f0a = f0a * f0b;
      f0b = f0a;
      f1a = f1a * f1b;
      f1b = f1a;

/* f0a = quotient, f0b = remainder */
      f0a = f0a * VectStuff[8];
      f1a = f1a * VectStuff[8];

/* round to nearest by justifying mantissa */
      f0a = f0a + VectStuff[4];
      f1a = f1a + VectStuff[4];
      f0a = f0a - VectStuff[5];
      f1a = f1a - VectStuff[5];

/* compute a*b-n*q */
      f0a = f0a * VectStuff[7];
      f1a = f1a * VectStuff[7];
      f0b = f0b - f0a;
      f1b = f1b - f1a;

/* push answer to bottom of 53-bit mantissa */
      f0b = f0b + VectStuff[6];
      f1b = f1b + VectStuff[6];
/* store to ChopShop temporary space */
      ChopShop[0].D=f0b;
      ChopShop[1].D=f1b;

      Chop(a[0],0);
      Chop(a[1],1);

      a+=2;b+=2;size-=2;
   } while(size);
}

void
VectorModMulC( ModInt *a, ModInt b, size_t size )
/*  Num1[0...len] <- Num1[0...len] * Num2 mod n */
{
   do {
       register LDouble f0a, f0b;
       register LDouble f1a, f1b;

      f0a = (double)a[0];
      f0b = (double)b;
      f1a = (double)a[1];
      f1b = (double)b;
      f0a = f0a * f0b;
      f0b = f0a;
      f1a = f1a * f1b;
      f1b = f1a;

/* f0a = quotient, f0b = remainder */
      f0a = f0a * VectStuff[8];
      f1a = f1a * VectStuff[8];

/* round to nearest by justifying mantissa */
      f0a = f0a + VectStuff[4];
      f1a = f1a + VectStuff[4];
      f0a = f0a - VectStuff[5];
      f1a = f1a - VectStuff[5];

/* compute a*b-n*q */
      f0a = f0a * VectStuff[7];
      f1a = f1a * VectStuff[7];
      f0b = f0b - f0a;
      f1b = f1b - f1a;

/* push answer to bottom of 53-bit mantissa */
      f0b = f0b + VectStuff[6];
      f1b = f1b + VectStuff[6];
/* store to ChopShop temporary space */
      ChopShop[0].D=f0b;
      ChopShop[1].D=f1b;

      Chop(a[0],0);
      Chop(a[1],1);

      a+=2;size-=2;
   } while(size);
}

void
VectorModButterfly(ModInt *Num1, ModInt *Num2, size_t Len)
/*
** Do a simple vector modular 'butterfly'.
**
** N1=Num1[x];
** Num1[x]=N1 + Num2[x] mod Prime;
** Num2[x]=N1 - Num2[x] mod Prime;
*/
{
if ((Len < 4) || (Len & 3))
  FatalError("VectorModButterfly called with incorrect length: %lu\n",(ULINT)Len);

while (Len)
  {ModInt L,R;

   L=Num2[0];
   R=Num1[0];
   L=L+Num1[0];
   R=R-Num2[0];
   L=L-JPrime;
   {INT32 TR,TL;
    TR=R>>31;
    TL=L>>31;
    TR=TR & JPrime;
    TL=TL & JPrime;
    R=R+TR;
    L=L+TL;
    Num1[0]=L;
    Num2[0]=R;
   }
   Num1+=1;Num2+=1;Len-=1;
  }
}

void
VectorNTT_First2(ModInt *Data, ModInt Trig, size_t Len)
/*
** Perform the first two passes of the DiT NTT
** It's a little more complicated than the others.
** It needs to be done because the vector code works on a granularity
** of 4 (meaning all operations need to be a multiple of 4)
** and the first two passes of the NTT would be smaller than that.
**
** Also, the code has been aranged to (strongly) encourage parallel
** behavior between the integer unit and the FPU unit while doing
** the one required modular multiply.
*/
{size_t x;
 ModInt D0,D1,D2,D3;
 double DTrig=Trig;

if ((Len < 4) || (Len & 3))
  FatalError("VectorNTT_First2 called with incorrect length: %lu\n",(ULINT)Len);

for (x=0;x<Len;x+=4)
  {ModInt T;
   register LDouble f0,f1;

   D2=Data[2];D3=Data[3];
   T=D3;
   D3=ModSub(D2,D3);        f0=D3;f0=f0*DTrig;
   D2=ModAdd(D2,T);         f1=f0;f0=f0*VectStuff[8];
   D0=Data[0];              f0=f0+VectStuff[4];
   D1=Data[1];              f0=f0-VectStuff[5];
   T=D0;                    f0=f0*VectStuff[7];
   D0=ModAdd(T,D1);         f1=f1-f0;
   D1=ModSub(T,D1);         f1=f1+VectStuff[6];

   T=D0;                    ChopShop[0].D=f1;
   D0=ModAdd(T,D2);
   D2=ModSub(T,D2);
   Chop(D3,0);

   T=D1;
   D1=ModAdd(T,D3);
   D3=ModSub(T,D3);
   Data[0]=D0;
   Data[1]=D1;
   Data[2]=D2;
   Data[3]=D3;

   Data+=4;
  }
}

void
VectorNTT_Last2(ModInt *Data, ModInt Trig, size_t Len)
/*
** Perform the Last two passes of the DiF NTT
** It's a little more complicated than the others.
** It needs to be done because the vector code works on a granularity
** of 4 (meaning all operations need to be a multiple of 4)
** and the last two passes of the NTT would be smaller than that.
*/
/*
** Yes, I know this isn't the 'c586' style.  I'll fix it later.
*/
{size_t x;
 ModInt D0,D1,D2,D3;

#if 0
if ((Len < 4) || (Len & 3))
  FatalError("VectorNTT_Last2 called with incorrect length: %lu\n",(ULINT)Len);
#endif

for (x=0;x<Len;x+=4)
  {ModInt T;
   D0=Data[0];D1=Data[1];D2=Data[2];D3=Data[3];

   T=D0;
   D0=ModAdd(T,D2);
   D2=ModSub(T,D2);

   T=D1;
   D1=ModAdd(T,D3);
   D3=ModSub(T,D3);
   D3=ModMul(D3,Trig);

   T=D0;
   D0=ModAdd(T,D1);
   D1=ModSub(T,D1);
   T=D2;
   D2=ModAdd(T,D3);
   D3=ModSub(T,D3);

   Data[0]=D0;Data[1]=D1;Data[2]=D2;Data[3]=D3;
   Data+=4;
  }
}


void
PrepVector(UINT32 Prime, size_t Len)
/*
** Prepare for a vector of length 'Len' with prime 'Prime'.  This
** will be valid until this routine is called again, not just for
** the next vector opeartion.
*/
{
 JPrime=Prime;
 VectStuff[7]=Prime;
 VectStuff[8]=1.0/Prime;
}

/*
** Initialize for vector mod math.  This also includes checking the
** rounding mode, which must be in 'round to nearest' mode.  If it
** isn't, then setting it will be system & compiler dependant.
**
** Even though FLT_ROUNDS can be a value (rather than a constant) and
** can change during the program run, rather than checking it in
** every VectorMod...() call, I'll assume that it's not going to change
** during the program run.
**
** Also see the docs and version.txt
*/
void
InitVectorModMath(size_t Len)
{volatile long double LD;

LD=1.0;
LD+=LDBL_EPSILON;
if (LD==1.0)
  FatalError("It seems that you can't do long double math.  Sorry.  See vector.c\n");

/*
** This assumes that the radix is 2!!  Too much trouble to try and code it
** for any other radix.
*/
if (LDBL_MANT_DIG < 64)
  FatalError("This version needs long doubles to be at least 64 bits.  See vector.c\n");

 if (FLT_ROUNDS != 1)
   FatalError("It appears that the FPU isn't rounding to nearest.\nSee docs, version.txt, and vector.c");

CHECK_ALIGNMENT(VectStuff);

VectStuff[4]=3.0*2147483648.0*2147483648.0;
VectStuff[5]=3.0*2147483648.0*2147483648.0;
VectStuff[6]=6755399441055744.0;
}

void
DeInitVectorModMath(void)
{
}


