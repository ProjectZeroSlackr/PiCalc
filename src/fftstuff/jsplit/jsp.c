/*
** Please note:
**
** This code was written to link into jason's C based split radix FFT
** pair.  He does the signs differently from what I do.  His DiF is a
** negative sign, and his DiT is positive.  Opposite from what I do.
** This also effects the sign in the real wrapper.
*/

typedef struct {FFT_DATA_TYPE r,i;} Cmplx;

#define POS (+1.0)
#define NEG (-1.0)

/* A couple of Konstants */
double K_PI_    =3.14159265358979323846L;
double K_SQRT05_=0.70710678118654752440L;

Cmplx *Twiddles=NULL;

#include "jsp.h"

#define RFButterfly(Ndx,Sr,Pr,Si,Pi)                        \
    {double temp_r,temp_i,t;size_t z=(Ndx);                 \
     Left[z].r  = (temp_r = Left[z].r) + (t = Right[z].r);  \
     temp_r  -= t;                                          \
     Left[z].i  = (temp_i = Left[z].i) + (t = Right[z].i);  \
     temp_i  -= t;                                          \
     Right[z].r = temp_r*(Pr)*(Sr) - temp_i*(Pi)*(Si);      \
     Right[z].i = temp_i*(Pr)*(Sr) + temp_r*(Pi)*(Si);      \
    }

void FwdRFFT_F(Cmplx *Data,int Len)
/*
** Recursive Decimation in frequency FFT
*/
{int x;
 Cmplx *Left,*Right;
 TRIG_VARS
 int Len2,Len4,Len8;

if (Len<=(CPU_CACHE/sizeof(Cmplx)))
  {size_t TwidLen=CPU_CACHE/sizeof(Cmplx);
   Cmplx *Trig=Twiddles;
   while (Len < TwidLen) {TwidLen/=2;Trig+=TwidLen;}
   single_fft(Data,Trig,Len);
   return;
  }

Len2=Len/2;Len4=Len/4;Len8=Len/8;
Left=&Data[0];Right=&Data[Len2];

INIT_TRIG(Len2,POS); // Yes, POS because we do the signs below.
RFButterfly(0,POS,1.0,NEG,0.0);
for (x=1;x<Len8;x++)
  {
    NEXT_TRIG_POW(POS);
    RFButterfly(x,     POS,Pow_r,NEG,Pow_i);
    RFButterfly(Len2-x,NEG,Pow_r,NEG,Pow_i);
    RFButterfly(Len4-x,POS,Pow_i,NEG,Pow_r);
    RFButterfly(Len4+x,NEG,Pow_i,NEG,Pow_r);
  }

if (Len8)
  {double sq=K_SQRT05_;
   RFButterfly(Len8,     POS,sq,NEG,sq);
   RFButterfly(Len2-Len8,NEG,sq,NEG,sq);
  }

if (Len4)
   RFButterfly(Len4,POS,0.0,NEG,1.0);

if (Len >= 4)
  {
   FwdRFFT_F(Data,     Len2);
   FwdRFFT_F(Data+Len2,Len2);
  }
}

/* The butterfly for the recursive FFTs */
#define RButterfly(Ndx,Sr,Pr,Si,Pi)                   \
   {int z=(Ndx);                                      \
    {double tr, ti, t;                                \
     tr = (ti = Right[z].r) * (Pr) * (Sr) -           \
          Right[z].i * (Pi) * (Si);                   \
     Right[z].r = (t = Left[z].r) - tr;               \
     Left[z].r  = t + tr;                             \
     ti = ti * (Pi) * (Si) +                          \
          Right[z].i * (Pr) *(Sr);                    \
     Right[z].i = (t = Left[z].i) - ti;               \
     Left[z].i  = t + ti;                             \
    }                                                 \
   }

static void
RevRFFT_T(Cmplx *Data,int Len)
/*
** Recursive Reverse Decimation in Time
*/
{int x,Len2,Len4,Len8;
 Cmplx *Left,*Right;
 TRIG_VARS;

Len2 = Len/2;Len4 = Len/4;Len8 = Len/8;

if (Len<=(CPU_CACHE/sizeof(Cmplx)))
  {size_t TwidLen=CPU_CACHE/sizeof(Cmplx);
   Cmplx *Trig=Twiddles;
   while (Len < TwidLen) {TwidLen/=2;Trig+=TwidLen;}
   single_ifft(Data,Trig,Len);
   return;
  }
if (Len2 >= 2) {RevRFFT_T(Data,Len2);RevRFFT_T(Data+Len2,Len2);}

INIT_TRIG(Len2,POS);
Left=&Data[0];
Right=&Data[Len2];
RButterfly(0,POS,1.0,POS,0.0);
for (x=1;x<Len8;x++)
  {
    NEXT_TRIG_POW(POS);
    RButterfly(x,     POS,Pow_r,POS,Pow_i);
    RButterfly(Len2-x,NEG,Pow_r,POS,Pow_i);
    RButterfly(Len4-x,POS,Pow_i,POS,Pow_r);
    RButterfly(Len4+x,NEG,Pow_i,POS,Pow_r);
  }

if (Len8)
   {/*long*/double sq=K_SQRT05_;
    RButterfly(Len8,     POS,sq,POS,sq);
    RButterfly(Len2-Len8,NEG,sq,POS,sq);
   }
if (Len4)
    RButterfly(Len4,POS,0,POS,1.0);
}


void
InitFFT(size_t Len)
{int TLen=(CPU_CACHE/sizeof(Cmplx));
 long i,inc=1;
 double root=2*K_PI_/TLen;
 Cmplx* Twid;

Twiddles = (Cmplx *)malloc( (TLen+2) * sizeof(Cmplx) );
if (Twiddles==NULL) FatalError("Unable to allocate memory for trig twiddles.\n");

Twid=Twiddles;

/*
** fills the array x[] with the twiddle factors needed
** for a split-radix FFT. The roots are stored in the order
** they're used; for simplicity, make x of size n even though
** this wastes a little memory.
*/
 while( TLen>2 )
   {
    for(i=0; i<TLen/4; i++)
      {
       Twid->r = cos(1*i*inc*root); Twid->i = -sin(1*i*inc*root); Twid++;
       Twid->r = cos(3*i*inc*root); Twid->i = -sin(3*i*inc*root); Twid++;
      }
    TLen = TLen/2;
    inc = inc*2;
   }
}

void
DeInitFFT(void)
{
free(Twiddles);
}


