#ifndef PI_H_
#define PI_H_ 1

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdarg.h>
#include <stddef.h>
#include <math.h>
#include <time.h>
#include <limits.h>
#include <float.h>
#include <errno.h>

typedef unsigned long int ULINT;

/*
** The standard C macro FILENAME_MAX might only be large enough for
** just the filename, and not any path.<disgust>  I need to allow
** enough for paths, too.
*/
#define MAX_FILENAME 128

#define LINE_INPUT_MAX 256

#include "config.h"
#include "sys_cfg.h"

/*
** Accumulate MSR counts, instead of RDTSR cycle counts
** Pentium & Pentium/MMX specific!  NOT Pentium II or III compatable.
** DJGPP specific.
**
** ##** MUST **## be running in 'ring 0' privelege level.  IE: Operating
** system privilege level.  Under DOS, that means using the cwsdpr0.exe
** instead of the regular DPMI, and that means no virtual memory will
** be available.
**
** You should also remember that these counters will count any event,
** even if it occurs in an interrupt handler, DPMI server, DOS, etc.
** So, these counts will not be 100% accurate, and will fluctuate
** some between runs.
**
** 0x3  Data read misses
** 0xa  Bank conflicts.
** 0xb  Misaligned data accesses
** 0x16 Instructions executed.
** 0x1a pipeline stall due to data read delays (cycles duration)
** 0x22 FLOPs.  Inlcudes int muls, int->float, etc.
** 0x32 FPU stalls due to read delays.  (MMX) (cycles duration)
#define USE_MSR 0x22
*/
#ifndef DO_PENTIUM_TIMINGS
#undef USE_MSR
#endif

#include "port.h"
#include "malloc.h"
#include "ini.h"
#include "block.h"
#include "bigint.h"

#define CHECK_ALIGNMENT(Var) if (((UINT32)&Var) % 8) \
    fprintf(stderr,"\n**Warning** Misalignment [%p] in file: %s, Line: %d\n",&Var,__FILE__,__LINE__);

//    fprintf(stderr,"\n**Warning** Misalignment [%p] at file: %s, Line: %d\n",((UINT32)(Var)),__FILE__,__LINE__);

#ifndef TRUE
#define TRUE 1
#endif
#ifndef FALSE
#define FALSE 0
#endif
#ifndef YES
#define YES 1
#endif
#ifndef NO
#define NO 0
#endif

#ifdef USE_DISK_NUM
#define NO_NUM 0
#else
#define NO_NUM NULL
#endif

/* how many decimal digits per RawInt (ie: INT32).  Just convenience. */
#define RawIntDigits 8

struct CfgStruct
  {
    int    UseCommandLine;
    int    DeleteSaveFile;
    int    Macintosh;
    int    DumpDebug;
    int    AlwaysSaveData;
    int    UseFastAGM;
    int    UseCRC;
    int    SavePiToFile;
    int    AGMSelfCheck;
    int    OutputFormat;
    size_t DigitsToCompute;
    size_t LenToCompute; /* Just needs to be here.  You shouldn't set it in the pi.ini */
    int    HalfMuls;     /* Just needs to be here.  You shouldn't set it in the pi.ini */
    int    AllowFractalMul;
    int    PiFormulaToUse;
    int    PassesToCompute;
  };

extern struct CfgStruct Cfg;

void  ClearTimings(void);
size_t DetectCacheSize(void *Core,size_t CoreSize);
void  DumpDebug(char *fmt,...);
void  DumpTimings(double Total);
void  ExitPrg(int EVal);
void  FatalError(char *fmt,...);
char *GetCheckStr(BigInt Num);
int   IsPow2(size_t Num);
int   Log2(size_t Num);
int   Max(size_t Num1, size_t Num2);
int   Min(size_t Num1, size_t Num2);
int   TestKeyboard(void);
char *Num2Str(size_t Num);
void  BackSpace(int b);
void  TestMath(size_t Digits);

extern char *PROG_VERSION;
extern char *PI_CFG_FILE;
extern BigInt DSWork, OldRoot;
extern size_t FatalPasses;

extern int Num1IsCached, Num2IsCached;
extern int SaveNum1FFT, SaveNum2FFT;

extern TIMER_TYPE FFTMulTime, FFTITime;
extern TIMER_TYPE FFTRTime, MulTime, ConvTime, SaveTime, CarryTime;
extern TIMER_TYPE FFTTime, LoadTime, CRTLoad, CRTTime;
extern TIMER_TYPE DiskIOTime,FFTDisk;

extern UINT32 BaseVar;
extern INT32 BI_One,BI_Two,BI_Three,BI_OneHalf;

/*
** Which iteration to redo for accuracy.  Should be divisible by
** 4 because of the Borwein quartic formula.  (Len/16) works very
** well, with good accuracy and little extra runtime.  You shouldn't
** really need to bother this.  It's just here because it needs to
** be somewhere.
*/
#define REDO_LEN (Len/16)


#endif

