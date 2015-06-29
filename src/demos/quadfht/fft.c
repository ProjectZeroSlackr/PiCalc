#include "pi.h"
#include "fft.h"
#include "cache.h"

#include <math.h>

static FFT_DATA_TYPE SineTable[50];
//static FFT_DATA_TYPE CoSineTable[50];
static FFT_DATA_TYPE Sqrt2;
static FFT_DATA_TYPE Sqrt05;

#define TRIG_VARS                   \
 FFT_DATA_TYPE Nth_r,Nth_i,Pow_r,Pow_i;

#define INIT_TRIG(LENGTH)           \
 {size_t x=Log2(LENGTH)+1;          \
  F_Set(&Pow_r,1,0);                \
  F_Set(&Pow_i,0,0);                \
  Nth_i=SineTable[x];               \
  Nth_r=SineTable[x+1];             \
  F_Mul(&Nth_r,&Nth_r,&Nth_r);      \
  DD_MulD(&Nth_r,-2.0,&Nth_r);      \
 }

#define NEXT_TRIG_POW               \
 {FFT_DATA_TYPE temp,T2;            \
   temp=Pow_r;                      \
   F_Mul(&Pow_r,&Pow_r,&Nth_r);     \
   F_Mul(&T2,&Pow_i,&Nth_i);        \
   F_Add(&Pow_r,&Pow_r,&temp);      \
   F_Sub(&Pow_r,&Pow_r,&T2);        \
   F_Mul(&T2,&Pow_i,&Nth_r);        \
   F_Add(&Pow_i,&Pow_i,&T2);        \
   F_Mul(&T2,&temp,&Nth_i);         \
   F_Add(&Pow_i,&Pow_i,&T2);        \
 }

// C=1, S=0 N1==N2
#define FHT_T1Butterfly(N1,N2,C,S)                \
 {size_t i1=N1,i2=N2;                             \
  FFT_DATA_TYPE cas1,temp;                        \
  cas1=Right[i1];                                 \
  temp=Left[i1];                                  \
  F_Add(&Left[i1],&temp,&cas1);                   \
  F_Sub(&Right[i2],&temp,&cas1);                  \
 }

#define FHT_T2Butterfly(N1,N2,C,S) \
 {FFT_DATA_TYPE Rx,Ri;             \
  size_t i1=N1,i2=N2;              \
  Rx=Right[i1];Ri=Right[i2];       \
  {FFT_DATA_TYPE cas1,Lx;          \
   F_Mul(&cas1,&Rx,&C);            \
   F_Mul(&Lx,&Ri,&S);              \
   F_Add(&cas1,&cas1,&Lx);         \
   Lx=Left[i1];                    \
   F_Add(&Left[i1],&Lx,&cas1);     \
   F_Sub(&Right[i1],&Lx,&cas1);    \
  }                                \
  {FFT_DATA_TYPE cas2,Li;          \
   F_Mul(&cas2,&Rx,&S);            \
   F_Mul(&Li,&Ri,&C);              \
   F_Sub(&cas2,&cas2,&Li);         \
   Li=Left[i2];                    \
   F_Add(&Left[i2],&Li,&cas2);     \
   F_Sub(&Right[i2],&Li,&cas2);    \
  }                                \
 }

/* Macro for the DiF FHT */
// C=1, S=0 N1==N2;
#define FHT_F1Butterfly(N1,N2,C,S)           \
 {FFT_DATA_TYPE D1,D2;                       \
  size_t i1=N1, i2=N2;                       \
  D1=Left[i1];D2=Left[i2];                   \
  {FFT_DATA_TYPE temp;                       \
   temp=Right[i2];                           \
   F_Add(&Left[i1],&D1,&temp);               \
   F_Sub(&Right[i2],&D1,&temp);              \
  }                                          \
 }

#define FHT_F2Butterfly(N1,N2,C,S)           \
 {FFT_DATA_TYPE D1,D2,temp;                  \
  size_t i1=N1, i2=N2;                       \
  D1=Left[i1];D2=Left[i2];                   \
   temp=Right[i1];                           \
   F_Add(&Left[i1],&D1,&temp);               \
   F_Sub(&D1,&D1,&temp);                     \
   temp=Right[i2];                           \
   F_Add(&Left[i2],&D2,&temp);               \
   F_Sub(&D2,&D2,&temp);                     \
  F_Mul(&temp,&D1,&C);                       \
  F_Mul(&Right[i1],&D2,&S);                  \
  F_Add(&Right[i1],&Right[i1],&temp);        \
  F_Mul(&temp,&D1,&S);                       \
  F_Mul(&Right[i2],&D2,&C);                  \
  F_Sub(&Right[i2],&temp,&Right[i2]);        \
 }

void
FHT_T(FFT_DATA_TYPE *Data, size_t Len)
/* Decimation in Time Hartley transform */
{
  size_t Step, Step2, Step4;
  TRIG_VARS
  size_t a,b;
  FFT_DATA_TYPE *Right,*Left;
StartTimer(FFTRTime);

  Step = 1;
  while (Step < Len)
    {
      Step *= 2;
      Step2 = Step/2;
      Step4 = Step/4;
      Left  = &Data[0];
      Right = &Data[Step2];

      INIT_TRIG(Step2);

      /* Do the special b=0 loop below */
      for (a = 0; a < Len; a+=Step)
         {
          FHT_T1Butterfly(a,a,1.0,0.0);
          if (Step4)
            FHT_T1Butterfly(a+Step4,a+Step4,1.0,0.0);
         }

      for (b = 1; b < Step4; b++)
        {size_t I1,I2;
         NEXT_TRIG_POW;
         for (I1=b,I2=Step2-b; I1 < Len; I1+=Step,I2+=Step)
            FHT_T2Butterfly(I1,I2,Pow_r,Pow_i);
        }
    }
StopTimer(FFTRTime);
}

void
FHT_F(FFT_DATA_TYPE *Data, size_t Len)
/* Decimation in frequency Hartley transform */
{
  size_t Step, Step2, Step4;
  TRIG_VARS
  size_t a,b;
  FFT_DATA_TYPE *Right,*Left;

StartTimer(FFTITime);
  Step = Len;
  while (Step > 1)
    {
      Step2 = Step/2;
      Step4 = Step/4;
      Left  = &Data[0];
      Right = &Data[Step2];

      INIT_TRIG(Step2);

      /* Do the special b=0 loop below */
      for (a = 0; a < Len; a+=Step)
        {
         FHT_F1Butterfly(a,a,1.0,0.0);
         if (Step4)
           FHT_F1Butterfly(a+Step4,a+Step4,1.0,0.0);
        }

      for (b = 1; b < Step4; b++)
        {size_t I1,I2;
         NEXT_TRIG_POW;
         for (I1=b,I2=Step2-b; I1 < Len; I1+=Step,I2+=Step)
            FHT_F2Butterfly(I1,I2,Pow_r,Pow_i);
        }
      Step /= 2;
    }
StopTimer(FFTITime);
}


void RFHT_F(FFT_DATA_TYPE *Data,int Len)
/*
** Recursive Decimation in Frequency style Fast Hartley Transform
*/
{int x,Len2,Len4;
 TRIG_VARS
 FFT_DATA_TYPE *Left,*Right;

Len/=2;Left=&Data[0];Right=&Data[Len];
if (Len==2)
  {FFT_DATA_TYPE d0=Data[0]; FFT_DATA_TYPE d1=Data[1];
   FFT_DATA_TYPE d2=Data[2]; FFT_DATA_TYPE d3=Data[3];
StartTimer(FFTITime);
StartTimer(CRTLoad);
   {FFT_DATA_TYPE d02; FFT_DATA_TYPE d13;
    DD_Add(&d0,&d2,&d02);
    DD_Add(&d1,&d3,&d13);
    DD_Add(&d02,&d13,&Data[0]);
    DD_Sub(&d02,&d13,&Data[1]);
   }
   {FFT_DATA_TYPE d02,d13;
    DD_Sub(&d0,&d2,&d02);
    DD_Sub(&d1,&d3,&d13);
    DD_Add(&d02,&d13,&Data[2]);
    DD_Sub(&d02,&d13,&Data[3]);
   }
StopTimer(FFTITime);
StopTimer(CRTLoad);
   return;
  }

if (Len<=(CPU_CACHE/sizeof(FFT_DATA_TYPE))) {FHT_F(Data,Len*2);return;}

StartTimer(FFTITime);

{FFT_DATA_TYPE t1,t2;
 t1=Left[0];t2=Right[0];
 DD_Add(&t1,&t2,&Left[0]);
 DD_Sub(&t1,&t2,&Right[0]);
 t1=Left[Len/2];t2=Right[Len/2];
 DD_Add(&t1,&t2,&Left[Len/2]);
 DD_Sub(&t1,&t2,&Right[Len/2]);
}

INIT_TRIG(Len);

Len2=Len/2;
Len4=Len/4;
for (x=1;x<Len4;x++)
  {
   NEXT_TRIG_POW;
   FHT_F2Butterfly(x,Len-x,Pow_r,Pow_i);
   FHT_F2Butterfly(Len2-x,Len2+x,Pow_i,Pow_r);
  }

/* Now do the two Len/4 points the loop missed */
if (Len4)
//  {FFT_DATA_TYPE sq=MY_SQRT_2; /* variable allows optimizations */
//   FHT_F2Butterfly(Len4,Len-Len4,sq,sq);
  {
   FHT_F2Butterfly(Len4,Len-Len4,Sqrt05,Sqrt05);
  }
StopTimer(FFTITime);

if (Len>=2) RFHT_F(Left, Len);
if (Len>=2) RFHT_F(Right,Len);
}

void RFHT_T(FFT_DATA_TYPE *Data,int Len)
/*
** recursive Decimation in Time style Fast Hartley Transform
*/
{int x,Len2,Len4;
 TRIG_VARS
 FFT_DATA_TYPE *Left,*Right;

Len/=2;Right=&Data[Len];Left=&Data[0];
if (Len<4) FatalError("OOps, how did rfht get called with %d\n",Len*2);
if (Len==4)
  {FFT_DATA_TYPE d45,d67,sd0123,dd0123;
StartTimer(FFTRTime);
StartTimer(CRTTime);
   {FFT_DATA_TYPE ss0123,ds0123,ss4567,ds4567;
    {FFT_DATA_TYPE s01,s23,d01,d23;
     DD_Sub(&Data[0],&Data[1],&d01);
     DD_Add(&Data[0],&Data[1],&s01);
     DD_Sub(&Data[2],&Data[3],&d23);
     DD_Add(&Data[2],&Data[3],&s23);
     DD_Sub(&s01,&s23,&ds0123);
     DD_Add(&s01,&s23,&ss0123);
     DD_Sub(&d01,&d23,&dd0123);
     DD_Add(&d01,&d23,&sd0123);
    }
    {FFT_DATA_TYPE s45,s67;
     DD_Add(&Data[4],&Data[5],&s45);
     DD_Add(&Data[6],&Data[7],&s67);
     DD_Sub(&Data[4],&Data[5],&d45);
     DD_Sub(&Data[6],&Data[7],&d67);
     DD_Sub(&s45,&s67,&ds4567);
     DD_Add(&s45,&s67,&ss4567);
    }
    DD_Sub(&ss0123,&ss4567,&Data[4]);
    DD_Add(&ss0123,&ss4567,&Data[0]);
    DD_Sub(&ds0123,&ds4567,&Data[6]);
    DD_Add(&ds0123,&ds4567,&Data[2]);
   }
   DD_Mul(&d45,&Sqrt2,&d45);
   DD_Mul(&d67,&Sqrt2,&d67);
   DD_Sub(&sd0123,&d45,&Data[5]);
   DD_Add(&sd0123,&d45,&Data[1]);
   DD_Sub(&dd0123,&d67,&Data[7]);
   DD_Add(&dd0123,&d67,&Data[3]);
StopTimer(FFTRTime);
StopTimer(CRTTime);
   return;
  }
if (Len<=(CPU_CACHE/sizeof(FFT_DATA_TYPE))) {FHT_T(Data,Len*2);return;}

RFHT_T(&Left[0], Len);
RFHT_T(&Right[0],Len);

StartTimer(FFTRTime);
/* Do the special x=0 loop below. */
FHT_T1Butterfly(0,0,1.0,0.0);

INIT_TRIG(Len);

Len2=Len/2;
Len4=Len/4;
for (x=1;x<Len4;x++)
  {
   NEXT_TRIG_POW;
   FHT_T2Butterfly(x,Len-x,Pow_r,Pow_i);
   FHT_T2Butterfly(Len2-x,Len2+x,Pow_i,Pow_r);
  }

/* Now do the two Len/4 points the loop missed */
if (Len4)
//  {FFT_DATA_TYPE sq=MY_SQRT_2; /* variable allows optimizations */
//   FHT_T2Butterfly(Len4,Len-Len4,sq,sq);
  {
   FHT_T2Butterfly(Len4,Len-Len4,Sqrt05,Sqrt05);
  }
/* Now do the Len/2 point the loop couldn't do. */
if (Len2)
  FHT_T1Butterfly(Len2,Len2,0.0,1.0);
StopTimer(FFTRTime);
}

/*
** For loading our numbers into the FFTNum array.
*/
static void
LoadNumBlockIntoNTT(size_t NX, INT32* NBuf, FFT_DATA_TYPE *FFTData)
{
while (NX)
  {
   NX-=RAW_FFT_DIG;
   DD_Set1((double)NBuf[NX],FFTData);
   FFTData++;
  }
}

static void
PutNumIntoFFTNum(FFT_DATA_TYPE *FFTNum, BigInt Num, size_t NumLen)
{size_t x,NX;
 size_t FFTLen=CalcFFTLen(NumLen);
 size_t NLen=NumLen;
 size_t FFTPos=0;
 INT32 *NBuf=(INT32*)FixedBuf;
 StartTimer(LoadTime);
 if (NumLen <8) FatalError("Too small of a FFT.\n");

  while (NLen)
    {
     NX=Min(FIXEDBUF_SIZE/sizeof(INT32),NLen);
     NLen-=NX;
     ReadNumIntoBuf(Num+NLen,NBuf,NX);
     LoadNumBlockIntoNTT(NX,NBuf,FFTNum+FFTPos);
     FFTPos+=(NX/RAW_FFT_DIG);
    }
  for (x = FFTLen/2; x < FFTLen; x++) F_Clear(&FFTNum[x]);

  StopTimer(LoadTime);
}

/*
** Do a convolution.  The data is scrambled.
*/
static void
DoScrambledFHTConv(FFT_DATA_TYPE *FFTNum1, FFT_DATA_TYPE *FFTNum2,size_t Len2)
{size_t Step=2;
 size_t Step2=Step*2;
 size_t x,y;
 double R=1.0/(2.0*Len2);

 while (Step < Len2)
  {
   for (x=Step,y=Step2-1;x<Step2;x+=2,y-=2)
     {FFT_DATA_TYPE h1p,h1m,h2p,h2m;
      FFT_DATA_TYPE s1,d1;
      h1p=FFTNum1[x];
      h1m=FFTNum1[y];
      F_Add(&s1,&h1p,&h1m);
      F_Sub(&d1,&h1p,&h1m);
      h2p=FFTNum2[x];
      h2m=FFTNum2[y];
      F_Mul(&h1p,&h2p,&s1);F_Mul(&h1m,&h2m,&d1);
      F_Add(&FFTNum1[x],&h1p,&h1m);
//      F_DivByFloat(&FFTNum1[x],&FFTNum1[x],2.0*Len2);
      DD_MulD(&FFTNum1[x],R,&FFTNum1[x]);
      F_Mul(&h1m,&h2m,&s1);F_Mul(&h1p,&h2p,&d1);
      F_Sub(&FFTNum1[y],&h1m,&h1p);
//      F_DivByFloat(&FFTNum1[y],&FFTNum1[y],2.0*Len2);
      DD_MulD(&FFTNum1[y],R,&FFTNum1[y]);
     }
   Step*=2;
   Step2*=2;
  }
 F_Mul(&FFTNum1[0],&FFTNum1[0],&FFTNum2[0]);
 F_Add(&FFTNum1[0],&FFTNum1[0],&FFTNum1[0]);
// F_DivByFloat(&FFTNum1[0],&FFTNum1[0],2.0*Len2);
 DD_MulD(&FFTNum1[0],R,&FFTNum1[0]);
 F_Mul(&FFTNum1[1],&FFTNum1[1],&FFTNum2[1]);
 F_Add(&FFTNum1[1],&FFTNum1[1],&FFTNum1[1]);
// F_DivByFloat(&FFTNum1[1],&FFTNum1[1],2.0*Len2);
 DD_MulD(&FFTNum1[1],R,&FFTNum1[1]);
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
#error
#if 0
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
#endif
StopTimer(ConvTime);
}


void
FwdTransform(FFT_DATA_TYPE *ddata,BigInt Num,size_t NumLen)
{
PutNumIntoFFTNum(ddata, Num, NumLen);
StartTimer(FFTTime);
RFHT_F(ddata,CalcFFTLen(NumLen));
//FHT_F(ddata,CalcFFTLen(NumLen));
StopTimer(FFTTime);
}

void
RevTransform(FFT_DATA_TYPE *ddata,size_t NumLen)
{
StartTimer(FFTTime);
//FHT_T(ddata,CalcFFTLen(NumLen));
RFHT_T(ddata,CalcFFTLen(NumLen));
StopTimer(FFTTime);
}

void
InitFFT(size_t Len)
/*
** If you need to use your own trig tables, you can calculate
** them here.
*/
{int x;double P;
 FFT_DATA_TYPE cosine,sine,One,C1,S1;
InitFixPoint();

DD_Set2(0.0,0.0,&SineTable[0]);
DD_Set2(0.0,0.0,&SineTable[1]);
DD_Set2(1.0,0.0,&SineTable[2]);

DD_Set2(0.0,0.0,&cosine);
DD_Set2(1.0,0.0,&sine);
DD_Set2(1.0,0.0,&One);

DD_Set2(2.0,0.0,&Sqrt2);DD_Sqrt(&Sqrt2,&Sqrt2);
DD_Set2(0.5,0.0,&Sqrt05);DD_Sqrt(&Sqrt05,&Sqrt05);

x=3;P=4;
while (P<=4.0*CalcFFTLen((double)Len) )
  {
   if (x>39) break;
/*   SineTable[x]=sin(M_PI/P);*/

   DD_Add(&One,&cosine,&cosine);
   DD_DivD(&cosine,2.0,&cosine);
//printf("Orig: ");DD_Dump(&cosine);
   DD_Sqrt(&cosine,&C1);
//printf("Sqrt: ");DD_Dump(&C1);
//{Quad z;DD_Mul(&C1,&C1,&z);printf("Sqr : ");DD_Dump(&z);}
   DD_Set(&C1,&cosine);
//printf("Cosine x=%2d ",x);DD_Dump(&cosine);

   DD_Add(&C1,&C1,&C1);
   DD_Div(&sine,&C1,&S1);
   DD_Set(&S1,&sine);
   DD_Set(&sine,&SineTable[x]);
/*
{Quad misc,check;
   DD_Mul(&sine,&sine,&misc);
   DD_Add(&misc,&misc,&misc);
   DD_Set2(1.0,0.0,&check);
   DD_Sub(&check,&misc,&check);
printf("Check  x=%2d ",x-1);DD_Dump(&check);
printf("\n");
}
*/
/*
printf("Sin()          %.20e\n",sin(M_PI/P));
printf("Sine   x=%2d ",x);DD_Dump(&sine);
//printf("\n");

printf("Cos()          %.20e\n",cos(M_PI/P));
printf("Cosine x=%2d ",x);DD_Dump(&cosine);
*/

   P*=2;
   x++;
  }

}

void
DeInitFFT(void)
{
DeInitFixPoint();
}



