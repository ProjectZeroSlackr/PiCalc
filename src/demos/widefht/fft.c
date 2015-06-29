#include "pi.h"
#include "fft.h"
#include "cache.h"

#include <math.h>

static FFT_DATA_TYPE SineTable[40];

#define TRIG_VARS                   \
 FFT_DATA_TYPE Nth_r,Nth_i,Pow_r,Pow_i;

#define INIT_TRIG(LENGTH)           \
 {size_t x=Log2(LENGTH)+1;          \
  F_Set(&Pow_r,1,0);                \
  F_Set(&Pow_i,0,0);                \
  Nth_i=SineTable[x];               \
  Nth_r=SineTable[x+1];             \
  F_Mul(&Nth_r,&Nth_r,&Nth_r);      \
  F_Add(&Nth_r,&Nth_r,&Nth_r);      \
  F_SetSign(&Nth_r,-1);             \
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
}

void
FHT_F(FFT_DATA_TYPE *Data, size_t Len)
/* Decimation in frequency Hartley transform */
{
  size_t Step, Step2, Step4;
  TRIG_VARS
  size_t a,b;
  FFT_DATA_TYPE *Right,*Left;

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
}

/*
** For loading our numbers into the FFTNum array.
** Remember, this is a fixed point FFT.  We need to scale our data.
** (See Knuth, Vol2, 3rd, p367)  To allow for a max transform
** length of 2^40, we also scale it by that amount.  Since 2^40
** isn't our base, we can instead scale it by 1.1e12.  But, it's
** easier to scale by 1.0e16 because that's just two words.  Since
** the first fixed point element is the integer part,
*/
static void
LoadNumBlockIntoNTT(size_t NX, INT32* NBuf, FFT_DATA_TYPE *FFTData)
{int z;
while (NX)
  {
   NX-=RAW_FFT_DIG;
   F_Clear(FFTData);
   for (z=0;z<RAW_FFT_DIG;z++)
     {
      FFTData->Data[3+z]=NBuf[NX+z];
     }
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
      F_DivByFloat(&FFTNum1[x],&FFTNum1[x],2.0*Len2);
      F_Mul(&h1m,&h2m,&s1);F_Mul(&h1p,&h2p,&d1);
      F_Sub(&FFTNum1[y],&h1m,&h1p);
      F_DivByFloat(&FFTNum1[y],&FFTNum1[y],2.0*Len2);
     }
   Step*=2;
   Step2*=2;
  }
 F_Mul(&FFTNum1[0],&FFTNum1[0],&FFTNum2[0]);
 F_Add(&FFTNum1[0],&FFTNum1[0],&FFTNum1[0]);
 F_DivByFloat(&FFTNum1[0],&FFTNum1[0],2.0*Len2);
 F_Mul(&FFTNum1[1],&FFTNum1[1],&FFTNum2[1]);
 F_Add(&FFTNum1[1],&FFTNum1[1],&FFTNum1[1]);
 F_DivByFloat(&FFTNum1[1],&FFTNum1[1],2.0*Len2);
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
FHT_F(ddata,CalcFFTLen(NumLen));
StopTimer(FFTTime);
}

void
RevTransform(FFT_DATA_TYPE *ddata,size_t NumLen)
{
StartTimer(FFTTime);
FHT_T(ddata,CalcFFTLen(NumLen));
StopTimer(FFTTime);
}

void
InitFFT(size_t Len)
/*
** If you need to use your own trig tables, you can calculate
** them here.
*/
{int x;double P;
 FIXPOINT C1,S1,sine,cosine;

InitFixPoint();

F_Set(&SineTable[0],0,0);
F_Set(&SineTable[1],0,0);
F_Set(&SineTable[2],1,0);
F_Set(&cosine,0,0);F_Set(&sine,1,0);
x=3;P=4;
while (P<=4.0*CalcFFTLen((double)Len) )
  {
   if (x>39) break;
/*   SineTable[x]=sin(M_PI/P);*/
/*   printf("%.39Le\n",(long double)sin(M_PI/P));*/

   F_AddOne(&cosine);F_DivBy2(&cosine,&cosine);
   F_Sqrt(&C1,&cosine);
   cosine=C1;
/*   F_Dump("cosine",&cosine);*/
   F_Add(&C1,&C1,&C1);
   F_Divide(&S1,&sine,&C1);
   sine=S1;
   SineTable[x]=sine;
/* Compute checking cosine
   printf("x=%d P=%f ",x,P);F_Dump("  sine",&sine);
   F_Mul(&misc,&sine,&sine);
   F_Add(&misc,&misc,&misc);
   F_Set(&pi,1,0);
   F_Sub(&pi,&pi,&misc);
   F_Dump("Check ",&pi);
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



