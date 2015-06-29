#include "pi.h"
#include "fft.h"
#include "cache.h"

#include <math.h>

#include "fht_f32.h"
#include "fht_t32.h"
#define FXT_FHT_F fht_f32
#define FXT_FHT_T fht_t32
#define HARDWIRED 32


/*
** Precalculated trig in long double (64 bits of mantissa).
** Perhaps not accurate to the last bit, but it's good enough.
*/
static FFT_DATA_TYPE SineTable[]={
/* sin(MY_PI     ) */      0.000000000000000000000000000000000000000e+00,
/* sin(MY_PI/1   ) */      0.000000000000000000000000000000000000000e+00,
/* sin(MY_PI/2   ) */      1.000000000000000000000000000000000000000e+00,
/* sin(MY_PI/4   ) */      7.071067811865475027412186737052479656995e-01,
/* sin(MY_PI/8   ) */      3.826834323650897575798401906155277174548e-01,
/* sin(MY_PI/16  ) */      1.950903220161282603309360617060974618653e-01,
/* sin(MY_PI/32  ) */      9.801714032956059818174621156572356994729e-02,
/* sin(MY_PI/64  ) */      4.906767432741801234115375240918410781887e-02,
/* sin(MY_PI/128 ) */      2.454122852291228707409531661909340982675e-02,
/* sin(MY_PI/256 ) */      1.227153828571992560101701352781589093865e-02,
/* sin(MY_PI/512 ) */      6.135884649154475120776119911880641666357e-03,
/* sin(MY_PI/1k  ) */      3.067956762965976150181468540267815114930e-03,
/* sin(MY_PI/2k  ) */      1.533980186284765552615777517431183696317e-03,
/* sin(MY_PI/4k  ) */      7.669903187427044969921158257264437452250e-04,
/* sin(MY_PI/8k  ) */      3.834951875713955740925670268026692610874e-04,
/* sin(MY_PI/16k ) */      1.917475973107032999678822626776764082024e-04,
/* sin(MY_PI/32k ) */      9.587379909597734213219655252657958044438e-05,
/* sin(MY_PI/64k ) */      4.793689960306688268003999509048185245774e-05,
/* sin(MY_PI/128k) */      2.396844980841821779338901565736819065933e-05,
/* sin(MY_PI/256k) */      1.198422490506970595472088780830688392598e-05,
/* sin(MY_PI/512k) */      5.992112452642427609071640315363538320526e-06,
/* sin(MY_PI/1m  ) */      2.996056226334660633689455089267994480906e-06,
/* sin(MY_PI/2m  ) */      1.498028113169011170304617541759739651752e-06,
/* sin(MY_PI/4m  ) */      7.490140565847156918673210856951527603087e-07,
/* sin(MY_PI/8m  ) */      3.745070282923841092784580930619142691285e-07,
/* sin(MY_PI/16m ) */      1.872535141461953375578708413939921229030e-07,
/* sin(MY_PI/32m ) */      9.362675707309807915379451515036635100842e-08,
/* sin(MY_PI/64m ) */      4.681337853654909086833363351942693952878e-08,
/* sin(MY_PI/128m) */      2.340668926827455184830686918395770135248e-08,
/* sin(MY_PI/256m) */      1.170334463413727672429463788628112297374e-08,
/* sin(MY_PI/512m) */      5.851672317068638462436019898049721632560e-09,
/* sin(MY_PI/1g  ) */      2.925836158534319244011595584353813137568e-09,
/* sin(MY_PI/2g  ) */      1.462918079267159623523680833656612776394e-09
};

#define K_PI_      3.1415926535897932384626433832
#define MY_PI      3.1415926535897932384626433832
#define MY_SQRT_2  0.7071067811865475244008443621
#define MY_SQRT2   1.4142135623730950488016887242

#define USE_RECURSIVE 1
#ifdef USE_RECURSIVE

#ifdef  FFT_TRIG_SIN
/* use explicit sin() calls.  Works on all processors & compilers */
#define TRIG_VARS                   \
 LDouble Theta,Angle;               \
 size_t TLen,TNdx;                  \
 LDouble Pow_r,Pow_i;

#define INIT_TRIG(LENGTH,DIR)       \
 Pow_r=1.0;Pow_i=0.0;               \
 Theta=K_PI_/(LENGTH);              \
 TLen=LENGTH;TNdx=0;                \
 Angle=0.0;

#define NEXT_TRIG_POW(DIR)          \
 Angle=(K_PI_*(++TNdx))/TLen;       \
 Pow_r=SINE(Angle*0.5);             \
 Pow_r=1.0-2.0*Pow_r*Pow_r;         \
 Pow_i=SINE(Angle)*(DIR);
#endif

#ifdef  FFT_TRIG_SIN_RECUR
/* Slightly faster than explicit trig.  Only occasionally calls sin() */
#define TRIG_VARS                   \
 size_t TLen,TNdx;                  \
 LDouble Pow_r,Pow_i,Nth_r,Nth_i;

#define INIT_TRIG(LENGTH,DIR)       \
 TNdx=0;TLen=LENGTH;                \
 Pow_r=1.0;Pow_i=0.0;               \
 Nth_r=SINE(K_PI_/((LENGTH)*2));    \
 Nth_r=-2.0*Nth_r*Nth_r;            \
 Nth_i=SINE(K_PI_/(LENGTH))*(DIR);

#define NEXT_TRIG_POW(DIR)          \
 if (((++TNdx)&7)==0)               \
   {LDouble Angle;                  \
    Angle=(K_PI_*(TNdx))/TLen;      \
    Pow_r=SINE(Angle*0.5);          \
    Pow_r=1.0-2.0*Pow_r*Pow_r;      \
    Pow_i=SINE(Angle)*(DIR);        \
   }                                \
 else                               \
   {LDouble temp;                   \
     temp = Pow_r;                  \
     Pow_r = Pow_r * Nth_r - Pow_i  \
             * Nth_i + Pow_r;       \
     Pow_i = Pow_i * Nth_r + temp   \
             * Nth_i + Pow_i;       \
   }
#endif

#ifdef  FFT_WIDE_FPU_RECUR
/*
** Use a simple trig recurance and only do sin() at init.
** Requires 64+ bit FPU registers and a compiler that will keep the vars
** in the FPU registers.
*/
#define TRIG_VARS                   \
 LDouble Nth_r,Nth_i,Pow_r,Pow_i;

/*
#define INIT_TRIG(LENGTH,DIR)       \
 Pow_r=1.0;Pow_i=0.0;               \
 Nth_r=SINE(K_PI_/((LENGTH)*2));    \
 Nth_r=-2.0*Nth_r*Nth_r;            \
 Nth_i=SINE(K_PI_/(LENGTH))*(DIR);
*/
#define INIT_TRIG(LENGTH,DIR)       \
 {size_t x=Log2(LENGTH)+1;          \
  Pow_r=1.0;Pow_i=0.0;              \
  Nth_i=SineTable[x];               \
  Nth_r=SineTable[x+1];             \
  Nth_r=-2.0*Nth_r*Nth_r;           \
 }


#define NEXT_TRIG_POW(DIR)          \
 {LDouble temp;                     \
   temp = Pow_r;                    \
   Pow_r = Pow_r * Nth_r - Pow_i    \
           * Nth_i + Pow_r;         \
   Pow_i = Pow_i * Nth_r + temp     \
           * Nth_i + Pow_i;         \
 }
#endif


/* Macro for the DiF FHT */
#define FHT_F2Butterfly(N1,N2,C,S)           \
 {FFT_DATA_TYPE D1,D2;                       \
  int i1=N1, i2=N2;                          \
  D1=Left[i1];D2=Left[i2];                   \
  {FFT_DATA_TYPE temp;                       \
   Left[i1] =D1+(temp=Right[i1]);D1=D1-temp; \
   Left[i2] =D2+(temp=Right[i2]);D2=D2-temp; \
  }                                          \
  Right[i1]=D1*(C)+D2*(S);                   \
  Right[i2]=D1*(S)-D2*(C);                   \
 }

void RFHT_F(FFT_DATA_TYPE *Data,int Len)
/*
** Recursive Decimation in Frequency style Fast Hartley Transform
*/
{int x,Len2,Len4;
 TRIG_VARS
 FFT_DATA_TYPE *Left,*Right;

#ifdef USE_HARD_WIRED
if (Len==HARDWIRED) {FXT_FHT_F(Data);return;}
#endif
Len/=2;Left=&Data[0];Right=&Data[Len];
if (Len==2)
  {FFT_DATA_TYPE d0=Data[0]; FFT_DATA_TYPE d1=Data[1];
   FFT_DATA_TYPE d2=Data[2]; FFT_DATA_TYPE d3=Data[3];
   {FFT_DATA_TYPE d02=d0+d2; FFT_DATA_TYPE d13=d1+d3;
    Data[0]=d02+d13; Data[1]=d02-d13;
   }
   {FFT_DATA_TYPE d02=d0-d2; FFT_DATA_TYPE d13=d1-d3;
    Data[2]=d02+d13; Data[3]=d02-d13;
   }
   return;
  }

{FFT_DATA_TYPE t1,t2;
 t1=Left[0];t2=Right[0];
 Left[0]=t1+t2;Right[0]=t1-t2;
 t1=Left[Len/2];t2=Right[Len/2];
 Left[Len/2]=t1+t2;Right[Len/2]=t1-t2;
}

INIT_TRIG(Len,1.0);

Len2=Len/2;
Len4=Len/4;
for (x=1;x<Len4;x++)
  {
   NEXT_TRIG_POW(1.0);
   FHT_F2Butterfly(x,Len-x,Pow_r,Pow_i);
   FHT_F2Butterfly(Len2-x,Len2+x,Pow_i,Pow_r);
  }

/* Now do the two Len/4 points the loop missed */
if (Len4)
  {FFT_DATA_TYPE sq=MY_SQRT_2; /* variable allows optimizations */
   FHT_F2Butterfly(Len4,Len-Len4,sq,sq);
  }

if (Len>=2) RFHT_F(Left, Len);
if (Len>=2) RFHT_F(Right,Len);
}


/* Macro for the DiT FHT */
#define FHT_T2Butterfly(N1,N2,C,S) \
 {FFT_DATA_TYPE Rx,Ri;           \
  int i1=N1,i2=N2;               \
  Rx=Right[i1];Ri=Right[i2];     \
  {FFT_DATA_TYPE cas1,Lx;        \
   cas1=Rx*(C)+Ri*(S);           \
   Lx=Left[i1];                  \
   Left[i1]  = Lx+cas1;          \
   Right[i1] = Lx-cas1;          \
  }                              \
  {FFT_DATA_TYPE cas2,Li;        \
   cas2=Rx*(S)-Ri*(C);           \
   Li=Left[i2];                  \
   Left[i2]  = Li+cas2;          \
   Right[i2] = Li-cas2;          \
  }                              \
 }

/* Macro for the DiT FHT */
#define FHT_T1Butterfly(N1,N2,C,S)         \
 {int i1=N1,i2=N2;                         \
  FFT_DATA_TYPE cas1=Right[i1]*(C)+Right[i2]*(S); \
  FFT_DATA_TYPE temp=Left[i1];                    \
  Left[i1] = temp + cas1;                  \
  Right[i2]= temp - cas1;                  \
 }

void RFHT_T(FFT_DATA_TYPE *Data,int Len)
/*
** recursive Decimation in Time style Fast Hartley Transform
*/
{int x,Len2,Len4;
 TRIG_VARS
 FFT_DATA_TYPE *Left,*Right;

#ifdef USE_HARD_WIRED
if (Len==HARDWIRED) {FXT_FHT_T(Data);return;}
#endif
Len/=2;Right=&Data[Len];Left=&Data[0];
if (Len==4)
  {FFT_DATA_TYPE d45,d67,sd0123,dd0123;
   {FFT_DATA_TYPE ss0123,ds0123,ss4567,ds4567;
    {FFT_DATA_TYPE s01,s23,d01,d23;
     d01 = Data[0] - Data[1];
     s01 = Data[0] + Data[1];
     d23 = Data[2] - Data[3];
     s23 = Data[2] + Data[3];
     ds0123 = (s01 - s23);
     ss0123 = (s01 + s23);
     dd0123 = (d01 - d23);
     sd0123 = (d01 + d23);
    }
    {FFT_DATA_TYPE s45,s67;
     s45 = Data[4] + Data[5];
     s67 = Data[6] + Data[7];
     d45 = Data[4] - Data[5];
     d67 = Data[6] - Data[7];
     ds4567 = (s45 - s67);
     ss4567 = (s45 + s67);
    }
    Data[4] = ss0123 - ss4567;
    Data[0] = ss0123 + ss4567;
    Data[6] = ds0123 - ds4567;
    Data[2] = ds0123 + ds4567;
   }
   d45 *= MY_SQRT2;
   d67 *= MY_SQRT2;
   Data[5] = sd0123 - d45;
   Data[1] = sd0123 + d45;
   Data[7] = dd0123 - d67;
   Data[3] = dd0123 + d67;
   return;
  }

RFHT_T(&Left[0], Len);
RFHT_T(&Right[0],Len);

/* Do the special x=0 loop below. */
FHT_T1Butterfly(0,0,1.0,0.0);

INIT_TRIG(Len,1.0);

Len2=Len/2;
Len4=Len/4;
for (x=1;x<Len4;x++)
  {
   NEXT_TRIG_POW(1.0);
   FHT_T2Butterfly(x,Len-x,Pow_r,Pow_i);
   FHT_T2Butterfly(Len2-x,Len2+x,Pow_i,Pow_r);
  }

/* Now do the two Len/4 points the loop missed */
if (Len4)
  {FFT_DATA_TYPE sq=MY_SQRT_2; /* variable allows optimizations */
   FHT_T2Butterfly(Len4,Len-Len4,sq,sq);
  }
/* Now do the Len/2 point the loop couldn't do. */
if (Len2)
  FHT_T1Butterfly(Len2,Len2,0.0,1.0);
}

#else
/* use iterative, but give them the same name. */

#ifndef M_SQRT2
#define M_SQRT2   1.4142135623730950488016887242
#endif
#ifndef M_SQRT_2
#define M_SQRT_2  0.7071067811865475244008443621
#endif


#define NextTrigPow(Pr,Pi,Nr,Ni)    \
 {/*long*/ double temp;             \
   temp = Pr;                       \
   Pr = Pr * Nr - Pi   * Ni + Pr;   \
   Pi = Pi * Nr + temp * Ni + Pi;   \
 }

#define FHT_T1Butterfly(N1,N2,C,S)         \
 {int i1=N1,i2=N2;                         \
  double cas1=Right[i1]*(C)+Right[i2]*(S); \
  double temp=Left[i1];                    \
  Left[i1] = temp + cas1;                  \
  Right[i2]= temp - cas1;                  \
 }

#define FHT_T2Butterfly(N1,N2,C,S) \
 {double Rx,Ri;                  \
  int i1=N1,i2=N2;               \
  Rx=Right[i1];Ri=Right[i2];     \
  {double cas1,Lx;               \
   cas1=Rx*(C)+Ri*(S);           \
   Lx=Left[i1];                  \
   Left[i1]  = Lx+cas1;          \
   Right[i1] = Lx-cas1;          \
  }                              \
  {double cas2,Li;               \
   cas2=Rx*(S)-Ri*(C);           \
   Li=Left[i2];                  \
   Left[i2]  = Li+cas2;          \
   Right[i2] = Li-cas2;          \
  }                              \
 }

void
RFHT_T(double *Data, int Len)
/* Decimation in Time Hartley transform */
{
  int Step, Step2, Step4;
 /* long */ double Sin0,Sin,Cos0,Cos;
  int a,b;
  int TrigIndex=1;
  double *Right,*Left;
  double *Last=&Data[Len];

  TrigIndex=1;
  Step = 1;
  while (Step < Len)
    {
      Step *= 2;
      Step2 = Step/2;
      Step4 = Step/4;
      Left  = &Data[0];
      Right = &Data[Step2];

      Sin0 = SineTable[TrigIndex++];
      Cos0 = SineTable[TrigIndex];Cos0 = -2.0 * Cos0 * Cos0;
      Sin=Sin0;Cos=1.0+Cos0;

      /* Do the special b=0 loop below */
      for (a = 0; a < Len; a+=Step)
         {
          FHT_T1Butterfly(a,a,1.0,0.0);
          if (Step4)
            FHT_T1Butterfly(a+Step4,a+Step4,1.0,0.0);
         }

      for (b = 1; b < Step4; b++)
        {int I1,I2;
         for (I1=b,I2=Step2-b; I1 < Len; I1+=Step,I2+=Step)
            FHT_T2Butterfly(I1,I2,Cos,Sin);
         NextTrigPow(Cos,Sin,Cos0,Sin0);
        }
    }
}

/* Macro for the DiF FHT */
#define FHT_F1Butterfly(N1,N2,C,S)           \
 {double D1,D2;                              \
  int i1=N1, i2=N2;                          \
  D1=Left[i1];D2=Left[i2];                   \
  {double temp;                              \
   Left[i1] =D1+(temp=Right[i2]);D1=D1-temp; \
  }                                          \
  Right[i2]=D1*(C)+D2*(S);                   \
 }

#define FHT_F2Butterfly(N1,N2,C,S)           \
 {double D1,D2;                              \
  int i1=N1, i2=N2;                          \
  D1=Left[i1];D2=Left[i2];                   \
  {double temp;                              \
   Left[i1] =D1+(temp=Right[i1]);D1=D1-temp; \
   Left[i2] =D2+(temp=Right[i2]);D2=D2-temp; \
  }                                          \
  Right[i1]=D1*(C)+D2*(S);                   \
  Right[i2]=D1*(S)-D2*(C);                   \
 }


void
RFHT_F(double *Data, int Len)
/* Decimation in frequency Hartley transform */
{
  int Step, Step2, Step4;
 /* long */ double Sin0,Sin,Cos0,Cos;
  int a,b;
  int TrigIndex;
  double *Right,*Left;

  TrigIndex=Log2(Len);
  Step = Len;
  while (Step > 1)
    {
      Step2 = Step/2;
      Step4 = Step/4;
      Left  = &Data[0];
      Right = &Data[Step2];

      Sin0 = SineTable[TrigIndex];
      Cos0 = SineTable[TrigIndex+1];Cos0 = -2.0 * Cos0 * Cos0;
      Sin=Sin0;Cos=1.0+Cos0;
      TrigIndex--;

      /* Do the special b=0 loop below */
      for (a = 0; a < Len; a+=Step)
        {
         FHT_F1Butterfly(a,a,1.0,0.0);
         if (Step4)
           FHT_F1Butterfly(a+Step4,a+Step4,1.0,0.0);
        }

      for (b = 1; b < Step4; b++)
        {int I1,I2;
         for (I1=b,I2=Step2-b; I1 < Len; I1+=Step,I2+=Step)
            FHT_F2Butterfly(I1,I2,Cos,Sin);
         NextTrigPow(Cos,Sin,Cos0,Sin0);
        }
      Step /= 2;
    }
}


#endif

/*
** For loading our numbers into the FFTNum array.
*/
static void
PutNumIntoFFTNum(FFT_DATA_TYPE *FFTNum, BigInt Num, size_t NumLen)
{size_t x;

 StartTimer(LoadTime);
 if (NumLen <8) FatalError("Too small of a FFT.\n");

 {INT32 *N;
  size_t FFTLen=CalcFFTLen(NumLen);
  N=(INT32*)(FFTNum+FFTLen);
  ReadNumIntoBuf(Num,N-NumLen,NumLen);
#if FFT_DIGITS==4
  for (x=0; x < FFTLen/2; x+=2)
    {
     --N;
     FFTNum[x  ]=*N % 10000;
     FFTNum[x+1]=*N / 10000;
    }
#else
  for (x=0; x < FFTLen/2; x+=4)
    {INT32 D;
     --N;
     D=*N;
     FFTNum[x  ]= D % 100;D=D / 100;
     FFTNum[x+1]= D % 100;D=D / 100;
     FFTNum[x+2]= D % 100;D=D / 100;
     FFTNum[x+3]= D;
    }
#endif
  for (x = FFTLen/2; x < FFTLen; x++) FFTNum[x] = 0.0;
 }
  StopTimer(LoadTime);
}


/*
** Do a convolution.  The data is scrambled.
*/
static void
DoScrambledFHTConv(FFT_DATA_TYPE *FFTNum1, FFT_DATA_TYPE *FFTNum2,size_t Len2)
{size_t Step=2;
 size_t Step2=Step*2;
 FFT_DATA_TYPE dn=0.5 / Len2;
 size_t x,y;

 while (Step < Len2)
  {
   for (x=Step,y=Step2-1;x<Step2;x+=2,y-=2)
     {FFT_DATA_TYPE h1p,h1m,h2p,h2m;
      FFT_DATA_TYPE s1,d1;
      h1p=FFTNum1[x];
      h1m=FFTNum1[y];
      s1=h1p+h1m;
      d1=h1p-h1m;
      h2p=FFTNum2[x];
      h2m=FFTNum2[y];
      FFTNum1[x]=(h2p*s1+h2m*d1)*dn;
      FFTNum1[y]=(h2m*s1-h2p*d1)*dn;
     }
   Step*=2;
   Step2*=2;
  }
 FFTNum1[0]   = FFTNum1[0]   * 2.0 * dn * FFTNum2[0];
 FFTNum1[1]   = FFTNum1[1]   * 2.0 * dn * FFTNum2[1];
}

void
DoConvolution(FFT_DATA_TYPE *FFTNum,int Cache,size_t Len2)
{
StartTimer(ConvTime);
if (Cache==-2) /* In memory convolution of both nums */
  {FFT_DATA_TYPE *FFTNum2=FFTNum+Len2;
   DumpDebug("Cm...");
   DoScrambledFHTConv(FFTNum,FFTNum2,Len2);
  }
else if (Cache==-1)  /* In memory self convolution. */
  {
   DumpDebug("C");
   DoScrambledFHTConv(FFTNum,FFTNum,Len2);
  }
#ifdef VIRTUAL_CACHE
else /* Virtual Mem based convolution. */
  {FFT_DATA_TYPE *Num1=FFTNum;
   FFT_DATA_TYPE *Num2=FFTCash[Cache].Mem;

   if (Num2==NULL)
     FatalError("Cache %d doesn't exist.\n",Cache);

   DumpDebug("C%d...",Cache);
   DoScrambledFHTConv(Num1,Num2,Len2);
  }
#else
/*
** Eyuk!  What a kludge!
*/
else
  {FILE *f;FFT_DATA_TYPE h1p,h1m;
   size_t Step=2;
   size_t Step2=Step*2;
   FFT_DATA_TYPE dn=0.5 / Len2;
   size_t x,y;
   size_t BufSize=(FIXEDBUF_SIZE/2);
   FFT_DATA_TYPE *Bufx=(FFT_DATA_TYPE*)FixedBuf;
   FFT_DATA_TYPE *Bufy=(FFT_DATA_TYPE*)(FixedBuf+BufSize);

   DumpDebug("C%d...",Cache);
   f=fopen(FFTCash[Cache].Name,"rb");
   if (f==NULL)
     FatalError("Unable to open '%s' for convolution.\n",FFTCash[Cache].Name);
   while (Step < Len2)
    {x=Step;y=Step2-1;
     while (x < Step2)
       {size_t Sz,xx,yy;
        Sz=Min(BufSize/sizeof(FFT_DATA_TYPE),Step2-x);
        fseek(f,x*sizeof(FFT_DATA_TYPE),SEEK_SET);
        fread(Bufx,sizeof(FFT_DATA_TYPE),Sz,f);
        fseek(f,(y-Sz+1)*sizeof(FFT_DATA_TYPE),SEEK_SET);
        fread(Bufy,sizeof(FFT_DATA_TYPE),Sz,f);
        xx=0;yy=Sz-1;
        while (xx < Sz)
          {FFT_DATA_TYPE h1p,h1m,h2p,h2m;
           FFT_DATA_TYPE s1,d1;
            h1p=Bufx[xx];h1m=Bufy[yy];
            s1=h1p+h1m;
            d1=h1p-h1m;
            h2p=FFTNum[x];
            h2m=FFTNum[y];
            FFTNum[x]=(h2p*s1+h2m*d1)*dn;
            FFTNum[y]=(h2m*s1-h2p*d1)*dn;
            x+=2;xx+=2;
            y-=2;yy-=2;
          }
       }
     Step*=2;
     Step2*=2;
     }
   fseek(f,0*sizeof(FFT_DATA_TYPE),SEEK_SET);fread(&h1p,sizeof(FFT_DATA_TYPE),1,f);
   fseek(f,1*sizeof(FFT_DATA_TYPE),SEEK_SET);fread(&h1m,sizeof(FFT_DATA_TYPE),1,f);
   FFTNum[0] = h1p * 2.0 * dn * FFTNum[0];
   FFTNum[1] = h1m * 2.0 * dn * FFTNum[1];
   if (ferror(f)) FatalError("Error convoluting '%s'\n",FFTCash[Cache].Name);
   fclose(f);
  }
#endif
StopTimer(ConvTime);
}


void
FwdTransform(FFT_DATA_TYPE *ddata,BigInt Num,size_t NumLen)
{
PutNumIntoFFTNum(ddata, Num, NumLen);
StartTimer(FFTTime);
RFHT_F(ddata,CalcFFTLen(NumLen));
StopTimer(FFTTime);
}

void
RevTransform(FFT_DATA_TYPE *ddata,size_t NumLen)
{
StartTimer(FFTTime);
RFHT_T(ddata,CalcFFTLen(NumLen));
StopTimer(FFTTime);
}

void
InitFFT(size_t Len)
/*
** If you need to use your own trig tables, you can calculate
** them here.
*/
{
#if (LDBL_DIG > 18) || (LDBL_MANT_DIG > 64)
int x;unsigned int P;
SineTable[0]=-9.0;  /* Just a wild value that's sure to cause an error if used. */
x=1;P=1;
while (P<=Len*4)
  {
   SineTable[x]=sin(MY_PI/P);
/*   printf("%.39Le\n",(long double)sin(MY_PI/P));*/
   P*=2;
   x++;
  }
#endif
}

void
DeInitFFT(void)
{
}



