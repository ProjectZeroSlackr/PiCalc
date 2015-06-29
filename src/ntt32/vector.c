#include "pi.h"
#include "ntt.h"
#include "modmath.h"

void
VectorModAdd(ModInt *Num1, ModInt *Num2, size_t Len)
/*
** Do a simple vector modular addition.
**
** Num1[x]=Num1[x] + Num2[x] mod Prime;
*/
{
#if 0
if ((Len < 4) || (Len & 3))
  FatalError("VectorModAdd called with incorrect length: %lu\n",(ULINT)Len);
#endif

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
/*
** Do a simple vector modular subtraction.
**
** Num1[x]=Num1[x] - Num2[x] mod Prime;
*/
{
#if 0
if ((Len < 4) || (Len & 3))
  FatalError("VectorModSub called with incorrect length: %lu\n",(ULINT)Len);
#endif

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
VectorModMul(ModInt *Num1, ModInt *Num2, size_t Len)
/*
** Do a simple vector modular multiplication.
**
** Num1[x]=Num1[x] + Num2[x] mod Prime;
*/
{
#if 0
if ((Len < 4) || (Len & 3))
  FatalError("VectorModMul called with incorrect length: %lu\n",(ULINT)Len);
#endif

while (Len)
  {
   Num1[0]=ModMul(Num1[0],Num2[0]);
   Num1[1]=ModMul(Num1[1],Num2[1]);
   Num1[2]=ModMul(Num1[2],Num2[2]);
   Num1[3]=ModMul(Num1[3],Num2[3]);
   Num1+=4;Num2+=4;Len-=4;
  }
}

void
VectorModMulC(ModInt *Num1, ModInt Num2, size_t Len)
/*
** Do a simple array modular multiplication by a constant
**
** Num1[x]=Num1[x] + Num2 mod Prime;
*/
{
#if 0
if ((Len < 4) || (Len & 3))
  FatalError("VectorModMulC called with incorrect length: %lu\n",(ULINT)Len);
#endif

while (Len)
  {
   Num1[0]=ModMul(Num1[0],Num2);
   Num1[1]=ModMul(Num1[1],Num2);
   Num1[2]=ModMul(Num1[2],Num2);
   Num1[3]=ModMul(Num1[3],Num2);
   Num1+=4;Len-=4;
  }
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
#if 0
if ((Len < 4) || (Len & 3))
  FatalError("VectorModButterfly called with incorrect length: %lu\n",(ULINT)Len);
#endif

while (Len)
  {ModInt L,R;
   L=*Num1;R=*Num2;
   *Num1=ModAdd(L,R);
   *Num2=ModSub(L,R);
   Num1++;Num2++;Len--;
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
*/
{size_t x;
 ModInt D0,D1,D2,D3;

#if 0
if ((Len < 4) || (Len & 3))
  FatalError("VectorNTT_First2 called with incorrect length: %lu\n",(ULINT)Len);
#endif

for (x=0;x<Len;x+=4)
  {ModInt T;
   D0=Data[0];D1=Data[1];D2=Data[2];D3=Data[3];
   T=D0;
   D0=ModAdd(T,D1);
   D1=ModSub(T,D1);
   T=D2;
   D2=ModAdd(T,D3);
   D3=ModSub(T,D3);

   T=D0;
   D0=ModAdd(T,D2);
   D2=ModSub(T,D2);

   D3=ModMul(D3,Trig);
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
}

void
InitVectorModMath(size_t Len)
/*
** Do whatever you need to initialize for the later vector
** operations.  This includes activating any hardware or setting
** any variables.
*/
{
}

void
DeInitVectorModMath(void)
/*
** Do what ever clean up you need to do at the end of the program.
*/
{
}


