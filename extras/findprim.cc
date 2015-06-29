/*
** This program is a stand alone program that finds all the suitable
** primes to use in a Number Theoretic Transform.
*/
/*
** This program uses GNU C++'s 'integer.h' package, which allows
** very large integers.  This lets you (slowly) find primes
** larger than your cpu word size.
*/

/*
** Set the maximun length of the transform you are wanting to do.
** This includes the doubling for the product length.
*/
#define MAX_TRANSFORM_LENGTH 8388608UL


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <integer.h>
#include <string.h>

Integer Prime;

int isprime(Integer Num)
/*
** Is 'Num' a prime number?
*/
{Integer d, q, r;

if (Num==2) return 1;
if (Num==3) return 1;
if ((Num % 2) == 0) return 0;
d = 1;
do
  {
    d+=2;
    q = Num / d;r = Num % d;
    if (r == 0)  return 0;
  } while (d < q);
return 1;
}

Integer NextPrime(Integer prime)
{
while (prime > 0) if (isprime(++prime)) return prime;
return 0;
}

Integer FindSpecialPrime(Integer num,Integer Max)
{Integer p,prevp,k;
k=(Max-1)/num;
while (k!=0)
  {p=k*num+1L;
/* I think you can use almost any prime, but I'm not sure. */
//   if (k&1) // if it's even, we've already done it higher up.
//   if ((k==2) || ((k&1)==1))
     if (isprime(p)) return p;
   --k;
  }
return 0;
}

Integer Mul(Integer a,Integer b)
{Integer p;
p=a;
p*=b;
return p % Prime;
}

Integer Add(Integer a,Integer b)
{Integer s;
s=a;
s+=b;
while (s >= Prime) s-=Prime;
return (Integer) s;
}

Integer Sub(Integer a,Integer b)
{Integer d;
d=a;
d-=b;
while (d<0) d+=Prime;
return (Integer)d;
}

Integer ModPow(Integer base,Integer expon)
{Integer prod,b;

if (expon<=0) return 1;

b=base;
while ((expon&1)==0)
  {
   b=Mul(b,b);
   expon>>=1;
  }
prod=b;

while ((expon>>=1)!=0)
  {
   b=Mul(b,b);
   if ((expon&1)==1) prod=Mul(prod,b);
  }
return prod;
}

int PrimeFactorize(Integer pf[],Integer num)
{Integer p;int x;
p=1;x=0;
do
  {
    p=NextPrime(p);
    if ((num % p)==0) {pf[x]=p;++x;}
  } while (num/p > p);
return x;
}

Integer FindPrimvRoot(Integer len,Integer modulus)
{Integer m,f,prod,root;
 int x,np;
 Integer pf[50];

m=modulus-1;

np=PrimeFactorize(pf,m);

root=2;
do
  {
   root=NextPrime(root);
   if ((root > modulus) || (root==0))
//     {printf("failed.\n");return 0;}
     return 0;

   for (x=0;x<np;x++)
     if (ModPow(root,m/pf[x])==1) break;
  } while (x < np);
return root;
}

Integer FindInverse(Integer Len,Integer Modulus)
{Integer i;
i=ModPow(Len,Modulus-2L);
/*if (Mul(Len*3,i)!=3L) printf("Wrong. %lu\n",Mul(Len*3,i));*/
return i;
}

char *
Bits(Integer Num,int NBits)
{static char str[80];
 Integer x;
str[0]=0;
x=1;while (NBits--) x*=2;
x/=2;
while (x!=0)
  {
   if ((x & Num)!=0) strcat(str,"1");
   else              strcat(str,"0");
   x/=2;
  }
return str;
}

int main(void)
{Integer num,temp;
 Integer BiggestPrime,k;
 Integer Length;
 Integer PrimvRoot,MulInverse,OtherRoot;
 Integer MAX;
 int NumBits=32,x;


for (Length=1073741824UL; Length >=MAX_TRANSFORM_LENGTH;Length/=2)
  {
   cout << "Max Len FFT: " << Length << "\n";
   MAX=ULONG_MAX;
   MAX=LONG_MAX; /* primes need to be 'signed' */
//   MAX =65536;MAX*=65536; /* 32 bits */
//   MAX*=65536;MAX*=65536; /* 64 bits */
   MAX=1;for (x=0;x<NumBits;x++) MAX*=2;
   while (1)
     {
      Prime=BiggestPrime=FindSpecialPrime(Length,MAX);
      if (Prime==0) break;
      PrimvRoot=FindPrimvRoot(Length,BiggestPrime);
      if (PrimvRoot!=0)
        {
         MulInverse=FindInverse(Length,BiggestPrime);
         cout << "Prime=" << BiggestPrime << " ";
//         printf("Prime=%10u ",BiggestPrime);
//         printf("k=%4u ",(BiggestPrime-1)/Length);
         cout << "PRoot=" << PrimvRoot << " ";
//         printf("PRoot=%4u ",PrimvRoot);
//         printf("W=%10u ",ModPow(PrimvRoot,BiggestPrime-1-(BiggestPrime-1)/Length));
//         printf("W^-1=%10u ",ModPow(PrimvRoot,(BiggestPrime-1)/Length));
//         printf("Inv=%10u ",MulInverse);
//         printf("Bits=%5.2f ",log(BiggestPrime)/log(2.0)-0.02);
//         printf("%s",Bits(BiggestPrime,32));
         cout << Bits(BiggestPrime,NumBits) << " ";
//         printf("\n");
         cout << "\n";
        }
      MAX=Prime-1;
     }
  }

return 0;
}


