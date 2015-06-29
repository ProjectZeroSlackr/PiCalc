#include "pi.h"
#include "bigint.h"
#include "agm.h"
#include "borwein.h"

char *PROG_VERSION="Carey Bloodworth's AGM/Borwein PI program version 2.4\n";

size_t DigitsToCompute=0;
size_t LenToCompute=0;
int Num1IsCached=0, Num2IsCached=0;
int SaveNum1FFT=0, SaveNum2FFT=0;

TIMER_TYPE FFTMulTime, FFTITime;
TIMER_TYPE FFTRTime, MulTime, ConvTime, SaveTime, CarryTime;
TIMER_TYPE FFTTime, LoadTime, CRTLoad, CRTTime;
TIMER_TYPE DiskIOTime,FFTDisk;
TIMER_TYPE UnusedTimers[3]; /* for alignment & future expansion. */

UINT32 BaseVar   =100000000;     /* our base.  For convenience. */
INT32 BI_One     = 10000000;
INT32 BI_Two     = 20000000;
INT32 BI_Three   = 30000000;
INT32 BI_OneHalf =  5000000;
BigInt DSWork, OldRoot;
size_t FatalPasses;
char *PI_CFG_FILE="pi.ini";

struct CfgStruct Cfg;



#ifdef USE_MSR

volatile unsigned long long MSR;

void ReadMSRReg(int r)
{
 asm volatile ("movl %1,%%ecx ; .byte 15,50 ; movl %%eax,%0 ; movl %%edx,4+%0 "
      : "=o" (MSR)
      : "g" (r)
      : "%eax", "%ebx", "%edx", "%ecx", "cc");
}

void WriteMSRReg(int r)
{
 asm volatile ("movl %1,%%ecx ; movl %0,%%eax ; movl 4+%0,%%edx ; .byte 15,48 "
      : "=o" (MSR)
      : "g" (r)
      : "%eax", "%ebx", "%edx", "%ecx", "cc");
}

void PrintMSRRegs(void)
{int Reg;
for (Reg=0x10;Reg<=0x13;Reg++)
  {
   if (Reg==3) continue;
   if (Reg==10) continue;
   if (Reg==15) continue;
   if (Reg==20) continue;
   printf("Reg 0x%02x=",Reg);fflush(stdout);
   ReadMSRReg(Reg);
   printf("%016llx\n",MSR);
  }
}

#endif

size_t
DetectCacheSize(void *Core,size_t CoreSize)
{size_t x;
 long int sum;
 size_t Size,BestSize;
 double Gain,BestGain;
 size_t Loops,PrevLoops;
 clock_t Start,ETime;
 long int volatile *Mem=(long int *)Core;

fprintf(stderr,"Calculating cache size.  This may take several minutes.\nPlease wait...\n");
CoreSize/=sizeof(long int);
if (CoreSize > 4194304) CoreSize=4194304;
BestSize=CoreSize;
BestGain=0.0;
PrevLoops=0;
for (Size=CoreSize;Size>=1024;Size/=2)
  {
   sum=0;x=0;Loops=0;
   ETime=clock();
   do {Start=clock();} while (Start==ETime);
   do
     {
      sum += Mem[x+30]; sum += Mem[x+ 1];
      sum += Mem[x+28]; sum += Mem[x+ 3];
      sum += Mem[x+26]; sum += Mem[x+ 5];
      sum += Mem[x+24]; sum += Mem[x+ 7];
      sum += Mem[x+22]; sum += Mem[x+ 9];
      sum += Mem[x+20]; sum += Mem[x+11];
      sum += Mem[x+18]; sum += Mem[x+13];
      sum += Mem[x+16]; sum += Mem[x+15];
      sum += Mem[x+14]; sum += Mem[x+17];
      sum += Mem[x+12]; sum += Mem[x+19];
      sum += Mem[x+10]; sum += Mem[x+21];
      sum += Mem[x+ 8]; sum += Mem[x+23];
      sum += Mem[x+ 6]; sum += Mem[x+25];
      sum += Mem[x+ 4]; sum += Mem[x+27];
      sum += Mem[x+ 2]; sum += Mem[x+29];
      sum += Mem[x+ 0]; sum += Mem[x+31];
      Mem[0]=sum;
      x=(x+32) % Size;
      Loops++;
      ETime=clock();
     } while (ETime<=Start+(5*CLOCKS_PER_SEC));

   if (Size==CoreSize) PrevLoops=Loops;
   else
     {
      Gain=((double)Loops)/((double)PrevLoops);
      if (Gain > BestGain)
        {
         BestGain=Gain;
         BestSize=Size;
        }
      PrevLoops=Loops;
     }
  }

/* Return half the size of what we estimated the cache size to be. */
return BestSize*sizeof(long int) / 2;
}



char *
Num2Str(size_t Num)
{static char Str[80];
if      (Num < 1024)         sprintf(Str,"%lu", (ULINT)Num);
else if (Num < 1048576)      sprintf(Str,"%luk",(ULINT)Num/1024UL);
else if (Num < 1073741824UL) sprintf(Str,"%lum",(ULINT)Num/1048576UL);
else                         sprintf(Str,"%lug",(ULINT)Num/1073741824UL);
return Str;
}

void
BackSpace(int b)
{int x;
for (x=0;x<b;x++) fprintf(stderr,"\b");
for (x=0;x<b;x++) fprintf(stderr," ");
for (x=0;x<b;x++) fprintf(stderr,"\b");
}

void
ExitPrg(int EVal)
{
if (Cfg.Macintosh)
   printf("Select Quit from the File menu to quit\n");
exit(EVal);
}

/*
** A fatal error of some sort occured.  Print out the info and
** abort the program.
*/
void
FatalError(char *fmt,...)
{va_list argptr;
 int Err=errno;

  fprintf(stderr,"\a");
  va_start(argptr, fmt);
  vfprintf(stderr, fmt, argptr);
  fflush(stderr);
  va_end(argptr);

  fprintf(stdout,"\a");
  va_start(argptr, fmt);
  vfprintf(stdout, fmt, argptr);
  fflush(stdout);
  va_end(argptr);

  fprintf(stderr,"Errno=%d\n",Err);
  fprintf(stdout,"Errno=%d\n",Err);

  ExitPrg(EXIT_FAILURE);
}


/*
** Conditionally printf some stuff to stderr
*/
void
DumpDebug(char *fmt,...)
{va_list argptr;

if (Cfg.DumpDebug)
  {
   va_start(argptr, fmt);
   vfprintf(stderr, fmt, argptr);
   fflush(stderr);
   va_end(argptr);
  }
}


/*
** Misc. minor support routines.
*/

/*
** The lower of two int numbers
*/
int
Min(size_t Num1, size_t Num2)
{
  if (Num1 < Num2) return Num1;
  return Num2;
}

/*
** The higher of two int numbers
*/
int
Max(size_t Num1, size_t Num2)
{
  if (Num1 > Num2) return Num1;
  return Num2;
}

/*
** Is a number a power of two?
*/
int
IsPow2(size_t Num)
{
  return ((Num & -Num) == Num);
}

/*
** What integer power of two is the number?
*/
int
Log2(size_t Num)
{int x=-1;
  if (Num==0) return 0;
  while (Num) {x++;Num/=2;}
  return x;
}

/*
** This gets 4 BigInt digit (16 regular digits) and converts
** it to a string for easy comparasion for checkings etc.
**
** It returns _static_ storage.
*/
char *
GetCheckStr(BigInt Num)
{static char Str[32];
sprintf(Str,"%08d%08d",GetBigIntDigit(Num,0),GetBigIntDigit(Num,1));
return Str;
}

#ifdef DO_PENTIUM_TIMINGS
char *comma(unsigned long long val)
{
  static char buf2[80];
  char buf1[80];
  int i, j, k;
  sprintf(buf2,"%llu", val);
  j=79;buf1[79]='\0';
  k=0;
  for (i=strlen(buf2)-1;i>=0;i--)
    {
     buf1[--j]=buf2[i];
     k++;
     if ((k==3) && (i!=0)) {buf1[--j]=',';k=0;}
    }
strcpy(buf2,&buf1[j]);
return buf2;
}
#endif

static void
DumpTimings2(double Total,FILE *F)
{
#ifdef DO_TIMINGS

#ifdef DO_PENTIUM_TIMINGS
#ifdef USE_MSR
fprintf(F,"Total Time=%.0f seconds.  Timings below are MSR events.\n",
       (Total);
#else

fprintf(F,"Total Time=%.0f seconds.  Timings below are CPU clock cycles.\n",
       (Total);
#endif

fprintf(F,"Mul    =%16s ",comma(MulTime));
fprintf(F,"FFTMul =%16s\n",comma(FFTMulTime));
fprintf(F,"FFT    =%16s ",comma(FFTTime));
fprintf(F,"FFTI   =%16s ",comma(FFTITime));
fprintf(F,"FFTR   =%16s\n",comma(FFTRTime));
fprintf(F,"Conv   =%16s ",comma(ConvTime));
fprintf(F,"Carry  =%16s ",comma(CarryTime));
fprintf(F,"CRT    =%16s\n",comma(CRTTime));
fprintf(F,"Load   =%16s ",comma(LoadTime));
fprintf(F,"Save   =%16s ",comma(SaveTime));
fprintf(F,"CRTLoad=%16s\n",comma(CRTLoad));
fprintf(F,"FFTDisk=%16s ",comma(FFTDisk));
fprintf(F,"DiskIO =%16s\n",comma(DiskIOTime));

#else /* Regular timings */

#ifdef __DJGPP__
fprintf(F,"Total Time=%.0f seconds.  Timings below are 18.2 ticks.\n",Total);
#else
fprintf(F,"Total Time=%.0f seconds.\n",Total);
#endif
fprintf(F,"MulTime=%.0f FFTMulTime=%.0f\n",
        (double)MulTime,(double)FFTMulTime);
fprintf(F,"FFTTime=%.0f FFTITime=%.0f FFTRTime=%.0f FFTDisk=%.0f\n",
        (double)FFTTime,(double)FFTITime,(double)FFTRTime,(double)FFTDisk);
fprintf(F,"ConvTime=%.0f CarryTime=%.0f CRTTime=%.0f DiskIO=%.0f\n",
        (double)ConvTime,(double)CarryTime,(double)CRTTime,(double)DiskIOTime);
fprintf(F,"LoadTime=%.0f SaveTime=%.0f CLoad=%.0f\n",
        (double)LoadTime,(double)SaveTime,(double)CRTLoad);
#endif

#endif
}

void
DumpTimings(double Total)
{
DumpTimings2(Total,stderr);
/*DumpTimings2(Total,stdout);*/
}

void
ClearTimings(void)
{
MulTime=0;
FFTMulTime=0;
FFTTime=0;FFTRTime=0;FFTITime=0;FFTDisk=0;
ConvTime=0;
CarryTime=0;
SaveTime=0;
LoadTime=0;
CRTLoad=0;
CRTTime=0;
DiskIOTime=0;
}

static struct CfgNameStruct SettingsVarList[]={
  {"Settings","UseCommandLine" ,&Cfg.UseCommandLine ,Cfg_Boolean  ,1},
  {"Settings","DeleteSaveFile" ,&Cfg.DeleteSaveFile ,Cfg_Boolean  ,1},
  {"Settings","Macintosh"      ,&Cfg.Macintosh      ,Cfg_Boolean  ,1},
  {"Settings","DumpDebug"      ,&Cfg.DumpDebug      ,Cfg_Boolean  ,1},
  {"Settings","AlwaysSaveData" ,&Cfg.AlwaysSaveData ,Cfg_Boolean  ,1},
  {"Settings","UseFastAGM"     ,&Cfg.UseFastAGM     ,Cfg_Boolean  ,1},
  {"Settings","UseCRC"         ,&Cfg.UseCRC         ,Cfg_Boolean  ,1},
  {"Settings","SavePiToFile"   ,&Cfg.SavePiToFile   ,Cfg_Boolean  ,1},
  {"Settings","OutputFormat"   ,&Cfg.OutputFormat   ,Cfg_Integer  ,1},
  {"Settings","AllowFractalMul",&Cfg.AllowFractalMul,Cfg_Boolean  ,1},
  {"Settings","AGMSelfCheck"   ,&Cfg.AGMSelfCheck   ,Cfg_Integer  ,1},
  {"RunTime" ,"PiFormulaToUse" ,&Cfg.PiFormulaToUse ,Cfg_Integer  ,1},
  {"RunTime" ,"PassesToCompute",&Cfg.PassesToCompute,Cfg_Integer  ,1},
  {"RunTime" ,"DigitsToCompute",&Cfg.DigitsToCompute,Cfg_ULInteger,1},
  {NULL,NULL,NULL,Cfg_String,0}};

void
ReadConfig(void)
{char Str[80];
if (ReadCfgItem(PI_CFG_FILE,"Program","Title",Str,Cfg_String,78)==0)
    fprintf(stderr,"\a\n** Unable to locate Program identification string in %s **\n\n",PI_CFG_FILE);

#ifdef NO_COMMAND_LINE
Cfg.UseCommandLine=FALSE;
#endif
Cfg.UseCommandLine=FALSE;
Cfg.DeleteSaveFile=FALSE;
Cfg.Macintosh=FALSE;
Cfg.DumpDebug=FALSE;
Cfg.AlwaysSaveData=FALSE;
Cfg.UseFastAGM=TRUE;
Cfg.UseCRC=FALSE;
Cfg.SavePiToFile=TRUE;
Cfg.AGMSelfCheck=0;
Cfg.OutputFormat=3;
Cfg.DigitsToCompute=0;
/*Cfg.LenToCompute=0;*/ /* Not read in pi.ini */
/*Cfg.HalfMuls=FALSE;*/ /* Not read in pi.ini */
Cfg.AllowFractalMul=TRUE;
Cfg.PiFormulaToUse=0;
Cfg.PassesToCompute=0;

if (!LoadConfigStruct(PI_CFG_FILE,SettingsVarList))
  FatalError("Error while reading the %s file.\n",PI_CFG_FILE);

if (Cfg.Macintosh) Cfg.UseCommandLine=FALSE;

#ifdef NO_COMMAND_LINE
Cfg.UseCommandLine=FALSE;
#endif
}

static void
ReadLine(char *Str,FILE *f)
{char *s;int x;
s=fgets(Str,LINE_INPUT_MAX-2,f);
if (s==NULL) {Str[0]='\0';return;}
x=strlen(Str);
if (x)
  if (Str[x-1]=='\n') Str[x-1]='\0';
}

static void
Usage(void)
{
  puts("This program computes large number of digits of pi.");
  puts("");
  puts("The first option should be the number of digits to compute.  This");
  puts("  must be a power of two (such as 1024 2048 4096 8192 16384 32768");
  puts("  65536 131072 262144 524288 1048576 2097152 4194304, and so on.)");
  puts("  This can be abreivated as 8k or 1m for 8192 and 1048576 etc.");
  puts("");
  puts("The second option specifies which pi formula to use.");
  puts("  1 or 'C' means to use the memory frugal classic AGM. (slowest)");
  puts("  2 or 'A' means to use the traditional 1/sqrt(2) AGM. (fastest)");
  puts("  3 or 'S' means to use the sqrt(3) AGM.               (fast)");
  puts("  4 or 'B' means to use the Borwein Quartic formula.   (slow)");
  puts("  Anything else will use the traditional Salamin 1/sqrt(2) AGM.");
  puts("");
  puts("The third option specifies how many iterations of the specified");
  puts("  formula to compute.  If you specify this, you must specify");
  puts("  the second option too.  If the program doesn't complete the");
  puts("  calculation, it will save the data to disk for later resumption.");
  puts("  The resumption will be automatic when you later run the program");
  puts("  with the same arguments.");
  puts("");
  puts("");

  ExitPrg(EXIT_FAILURE);
}

int
main(int argc, char *argv[])
{int DoTimings=0;
#ifdef USE_MSR
 TIMER_TYPE Total_MSR=0;
#endif

  ReadConfig();

  if ((Cfg.DigitsToCompute!=0) &&
      (Cfg.PiFormulaToUse!=0)  &&
      (Cfg.PassesToCompute!=0))
    {
     Cfg.LenToCompute=(Cfg.DigitsToCompute+RawIntDigits-1)/RawIntDigits;
     FatalPasses=Log2(Cfg.LenToCompute)+4;
     if (Cfg.PassesToCompute==-1) Cfg.PassesToCompute=FatalPasses;
    }
  else
    {
     if (!Cfg.UseCommandLine)
       {char InputStr[LINE_INPUT_MAX];

        printf("#============================================================#\n");
        printf("I                                                            I\n");
        printf("I  PPPPPP  IIIIIII            AAAAA    GGGGG  MMM         M  I\n");
        printf("I  PPP   P   III             AAA   A  GGG   G MMMM       MM  I\n");
        printf("I  PPP   P   III             AAA   A  GGG     MMM M     M M  I\n");
        printf("I  PPPPPP    III     ãããã    AAAAAAA  GGG     MMM  M   M  M  I\n");
        printf("I  PPP       III             AAA   A  GGG  GG MMM   M M   M  I\n");
        printf("I  PPP       III             AAA   A  GGG   G MMM    M    M  I\n");
        printf("I  PPP     IIIIIII           AAA   A   GGGGG  MMM         M  I\n");
        printf("I                                                            I\n");
        printf("I Pi-AGM v2.4 by Carey Bloodworth                            I\n");
        printf("#============================================================#\n");

        do
          {int Digits,Len;
           printf("\nHow many digits of Pi do you want to calculate?  It should be\n");
           printf("a power of two, and can be followed by a k meaning *1024, or\n");
           printf("followed by a m meaning *1048576.  Ex: 1k 1024 2k 8k 1m\n");
           printf("Digits of pi: ");fflush(stdout);
           ReadLine(InputStr,stdin);
           Digits = atoi(InputStr);
           if (Digits == 0) ExitPrg(EXIT_SUCCESS);
           Len=strlen(InputStr);
           if (Len)
             {
              if (tolower(InputStr[Len-1])=='k') Digits*=1024;
              else if (tolower(InputStr[Len-1])=='m') Digits*=1048576L;
              else if (tolower(InputStr[Len-1])=='g') Digits*=1073741824L;
             }
           if (Digits < 0) DoTimings=Digits=-Digits;
           Len = (Digits + RawIntDigits-1) / RawIntDigits;
           Cfg.DigitsToCompute=Digits;
           Cfg.LenToCompute=Len;
          } while ((Cfg.LenToCompute < 64) || (!IsPow2(Cfg.LenToCompute)));
        FatalPasses=Cfg.PassesToCompute=Log2(Cfg.LenToCompute)+4;

        printf("\nWould formula would you like to use?\n");
        printf("1) The classic, memory frugal Salamin AGM\n");
        printf("2) The fast Salamin AGM\n");
        printf("3) The 'sqrt(3.0)' AGM\n");
        printf("4) The Borwein Quartic formula\n");
        printf("(1-4)? ");fflush(stdout);
        ReadLine(InputStr,stdin);
        Cfg.PiFormulaToUse=2;
        if (InputStr[0]=='1') Cfg.PiFormulaToUse=1;
        if (InputStr[0]=='2') Cfg.PiFormulaToUse=2;
        if (InputStr[0]=='3') Cfg.PiFormulaToUse=3;
        if (InputStr[0]=='4') Cfg.PiFormulaToUse=4;
        if (tolower(InputStr[0])=='b') Cfg.PiFormulaToUse=4;
        if (tolower(InputStr[0])=='a') Cfg.PiFormulaToUse=2;
        if (tolower(InputStr[0])=='c') Cfg.PiFormulaToUse=1;
        if (tolower(InputStr[0])=='s') Cfg.PiFormulaToUse=3;

        printf("\nThere will be an estimated %lu passes.  How many of\n",(ULINT)FatalPasses-4);
        printf("those passes do you wish to do (-1 for all)? ");fflush(stdout);
        ReadLine(InputStr,stdin);
        if ((InputStr[0]) && (InputStr[0]!='\n'))
          Cfg.PassesToCompute=atoi(InputStr);
        if (Cfg.PassesToCompute==-1) Cfg.PassesToCompute=FatalPasses;

        printf("\n\n");
       } /* Not predefined, no command line */
     else /* Not predefined, use command line */
       {int Digits,Len;
        if (argc < 2) Usage();

        if ((strcmp(argv[1],"-?")==0)   ||
            (strcmp(argv[1],"?")==0)    ||
            (strcmp(argv[1],"/?")==0)   ||
            (strcmp(argv[1],"\\?")==0)  ||
            (strcmp(argv[1],"-h")==0)   ||
            (strcmp(argv[1],"-H")==0)   ||
            (strcmp(argv[1],"help")==0) ||
            (strcmp(argv[1],"HELP")==0) ||
            (strcmp(argv[1],"/h")==0)   ||
            (strcmp(argv[1],"/H")==0)   ||
            (strcmp(argv[1],"\\h")==0)  ||
            (strcmp(argv[1],"\\H")==0)) Usage();

        Digits = atoi(argv[1]);
        Len=strlen(argv[1]);
        if (Len)
          {
           if (tolower(argv[1][Len-1])=='k')      Digits*=1024;
           else if (tolower(argv[1][Len-1])=='m') Digits*=1048576L;
           else if (tolower(argv[1][Len-1])=='g') Digits*=1073741824L;
          }
        if (Digits < 0) DoTimings=Digits=-Digits;
        Len = (Digits + RawIntDigits-1) / RawIntDigits;
        Cfg.DigitsToCompute=Digits;
        Cfg.LenToCompute=Len;

        if (Cfg.LenToCompute < 16) Usage();
        if (!IsPow2(Cfg.LenToCompute)) Usage();

        Cfg.PiFormulaToUse=2;
        if (argc > 2)
          {
           if (argv[2][0]=='1') Cfg.PiFormulaToUse=1;
           if (argv[2][0]=='2') Cfg.PiFormulaToUse=2;
           if (argv[2][0]=='3') Cfg.PiFormulaToUse=3;
           if (argv[2][0]=='4') Cfg.PiFormulaToUse=4;
           if (tolower(argv[2][0]=='b')) Cfg.PiFormulaToUse=4;
           if (tolower(argv[2][0]=='a')) Cfg.PiFormulaToUse=2;
           if (tolower(argv[2][0]=='c')) Cfg.PiFormulaToUse=1;
           if (tolower(argv[2][0]=='s')) Cfg.PiFormulaToUse=3;
           /* secret cheat to force fast/slow AGM */
           if (argv[2][1]=='f') Cfg.UseFastAGM=TRUE;
           if (argv[2][1]=='s') Cfg.UseFastAGM=FALSE;
          }

       FatalPasses=Cfg.PassesToCompute=Log2(Cfg.LenToCompute)+4;
       if (argc > 3) Cfg.PassesToCompute=atoi(argv[3]);
      } /* Not predefined, use command line. */
    } /* Not predefined. */

  if (Cfg.DigitsToCompute > 1073741824UL)
    FatalError("The program has an upper limit of 1 billion (1024m) digits.\n");

  setbuf(stderr,NULL);

  Cfg.HalfMuls=FALSE;
  if ((Cfg.UseFastAGM) && (!DoTimings))
    if ((Cfg.PiFormulaToUse==2) || (Cfg.PiFormulaToUse==3))
      Cfg.HalfMuls=TRUE;

  CheckDiskCache();
  InitBigIntPkg(Cfg.LenToCompute);

  DSWork  = CreateBigInt(Cfg.LenToCompute);
  OldRoot = CreateBigInt(Cfg.LenToCompute);

  ClearTimings();

#ifdef USE_MSR
MSR=0;WriteMSRReg(0x12);
ReadMSRReg(0x11);
MSR&= ~511;
MSR|= (3 << 6);MSR|= (USE_MSR);
WriteMSRReg(0x11);
MSR=0;WriteMSRReg(0x12);
StartTimer(Total_MSR);
#endif

  /*
  ** Now actually do the pi formula
  */
  if (DoTimings)                  TestMath(Cfg.LenToCompute);
  else if (Cfg.PiFormulaToUse==4) ComputeBorwein(Cfg.LenToCompute,Cfg.PassesToCompute);
  else                            ComputeAGM(Cfg.LenToCompute,Cfg.PassesToCompute);


#ifdef USE_MSR
StopTimer(Total_MSR);
printf("Total MSR Events=%s\n",comma(Total_MSR));
ReadMSRReg(0x11);
MSR&= ~511;
WriteMSRReg(0x11);
#endif

  DeInitBigIntPkg();

  ExitPrg(EXIT_SUCCESS);
  return EXIT_SUCCESS;
}


