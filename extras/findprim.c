/*
** This program is a stand alone program that finds all the suitable
** primes to use in a Number Theoretic Transform.
*/
/*
** This program depends on unsigned numbers rolling over at its limit
** to a lower value.
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

typedef unsigned long int USingle;
typedef unsigned long long int UDouble;

USingle Prime;

int isprime(USingle Num)
/*
** Is 'Num' a prime number?
*/
{USingle d, q, r;

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

USingle NextPrime(USingle prime)
{
while (prime > 0) if (isprime(++prime)) return prime;
return 0;
}

USingle FindSpecialPrime(USingle num,USingle Max)
{USingle p,prevp,k;
k=(Max-1)/num;
/*k=(LONG_MAX-1)/num; */
while (k!=0)
  {p=k*num+1L;
/* I think you can use almost any prime, but I'm not sure. */
//   if (k&1) // if it's even, we've already done it higher up.
   if ((k==2) || (k&1))
     if (isprime(p)) return p;
   k--;
  }
return 0;
}

USingle GetFactor(USingle r,USingle m)
{USingle quot,rem;
while (r)
  {
   if (isprime(++r))
     {
      quot=m/r;
      rem=m%r;
      if (rem==0) return r;
      if (quot < r) return 0;
     }
  }
return 0;
}

USingle Mul(USingle a,USingle b)
{UDouble p;
p=a;
p*=b;
return p % Prime;
}

USingle Add(USingle a,USingle b)
{UDouble s;
s=a;
s+=b;
while (s >= Prime) s-=Prime;
return (USingle) s;
}

USingle Sub(USingle a,USingle b)
{UDouble d;
d=a;
d-=b;
while (d<0) d+=Prime;
return (USingle)d;
}

USingle ModPow(USingle base,USingle expon)
{USingle prod,b;

if (expon<=0) return 1;

b=base;
while (!(expon&1))
  {
   b=Mul(b,b);
   expon>>=1;
  }
prod=b;

while (expon>>=1)
  {
   b=Mul(b,b);
   if (expon&1) prod=Mul(prod,b);
  }
return prod;
}

int PrimeFactorize(USingle pf[],USingle num)
{USingle p,x;
p=1;x=0;
do
  {
    p=NextPrime(p);
    if ((num % p)==0) pf[x++]=p;
  } while (num/p > p);
return x;
}

USingle FindPrimvRoot(USingle len,USingle modulus)
{USingle m,f,prod,root;
 int x,np;
 USingle pf[50];

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

USingle FindInverse(USingle Len,USingle Modulus)
{USingle i;
i=ModPow(Len,Modulus-2L);
/*if (Mul(Len*3,i)!=3L) printf("Wrong. %lu\n",Mul(Len*3,i));*/
return i;
}

char *
Bits(USingle Num)
{static char str[80];
 unsigned int x;
str[0]=0;
x=1<<31;
while (x)
  {
   if (x & Num) strcat(str,"1");
   else         strcat(str,"0");
   x/=2;
  }
return str;
}

int Log2(unsigned int Num)
{int x=-1;
  if (Num==0) return 0;
  while (Num) {x++;Num/=2;}
  return x;
}


int main(void)
{USingle num,temp;
 USingle BiggestPrime,k;
 USingle Length;
 USingle PrimvRoot,MulInverse,OtherRoot;
 USingle MAX;

for (Length=1073741824UL; Length >=MAX_TRANSFORM_LENGTH;Length/=2)
  {
   printf("Max Len FFT: %u Pow2=%d\n",Length,Log2(Length));
   MAX=LONG_MAX; /* primes need to be 'signed' */
   MAX=ULONG_MAX;
   while (1)
     {
      Prime=BiggestPrime=FindSpecialPrime(Length,MAX);
      if (Prime==0) break;
      PrimvRoot=FindPrimvRoot(Length,BiggestPrime);
      if (PrimvRoot)
        {
         MulInverse=FindInverse(Length,BiggestPrime);
         printf("Prime=%10u ",BiggestPrime);
         printf("k=%4u ",(BiggestPrime-1)/Length);
         printf("PRoot=%4u ",PrimvRoot);
/*         printf("W=%10u ",ModPow(PrimvRoot,BiggestPrime-1-(BiggestPrime-1)/Length));*/
/*         printf("W^-1=%10u ",ModPow(PrimvRoot,(BiggestPrime-1)/Length));*/
/*         printf("Inv=%10u ",MulInverse); */
         printf("Bits=%5.2f ",log(BiggestPrime)/log(2.0)-0.02);
/*         printf("%s",Bits(BiggestPrime));*/
         printf("\n");
        }
      MAX=Prime-1;
     }
  }

return 0;
}


