#include "pi.h"
#include "cache.h"

/*
** FFT caching
*/
struct FFTHdr *FFTCash=NULL;
int    MaxFFTCacheLines=10;

/*
** Flush our secondary FFT Cache.  Flushes all _but_ the 'Tag' one.
** A 'tag' less than or equal to 0 means flush them all.
*/
void
FlushFFTCache(int Tag)
{int x;
for (x=0;x<MaxFFTCacheLines;x++)
  {
   if ((Tag>0) && (FFTCash[x].Tag==Tag)) continue;
#ifndef VIRTUAL_CACHE
   if (FFTCash[x].FFTLen > 0) remove(FFTCash[x].Name);
#endif
   FFTCash[x].NumLen=0;
   FFTCash[x].FFTLen=0;
   FFTCash[x].ChkNum=NO_NUM;
   FFTCash[x].Tag=-1;
   FFTCash[x].Pass=-1;
#ifdef VIRTUAL_CACHE
   if (FFTCash[x].Mem) AlignedFree(FFTCash[x].Mem);
   FFTCash[x].Mem=NULL;
#endif
  }
}

/*
** Delete a particular cache.
*/
void
DeleteFFTCache(int Tag)
{int x;
for (x=0;x<MaxFFTCacheLines;x++)
  {
   if (Tag != FFTCash[x].Tag) continue;
#ifndef VIRTUAL_CACHE
   if (FFTCash[x].FFTLen > 0) remove(FFTCash[x].Name);
#endif
   FFTCash[x].NumLen=0;
   FFTCash[x].FFTLen=0;
   FFTCash[x].ChkNum=NO_NUM;
   FFTCash[x].Tag=-1;
   FFTCash[x].Pass=-1;
#ifdef VIRTUAL_CACHE
   if (FFTCash[x].Mem) AlignedFree(FFTCash[x].Mem);
   FFTCash[x].Mem=NULL;
#endif
   break;
  }
}

/*
** Check to see if a number has already been cached.
*/
int
CheckFFTCache(BigInt Num, size_t NumLen, int Tag,int Pass)
{int x;

/* Never match with entry 0 because that's special */
for (x=1;x<=MaxFFTCacheLines;x++)
  {
   if ((FFTCash[x].Tag==Tag) &&
       (FFTCash[x].Pass==Pass) &&
       (FFTCash[x].NumLen==NumLen) &&
       (FFTCash[x].ChkNum==Num))
     return x;
  }
return 0;
}

void
LoadFFTFromCache(int Cache,void *FFTNum)
{size_t FFTLen;
#ifndef VIRTUAL_CACHE
 FILE *f;
 size_t x;
#endif

if (Cache > MaxFFTCacheLines)
  FatalError("Tried to load too high of a cache line: %d\n",Cache);
FFTLen=FFTCash[Cache].FFTLen;

DumpDebug("Load%d (%d)...",Cache,FFTCash[Cache].Tag);
#ifdef VIRTUAL_CACHE
if (FFTCash[Cache].Mem)
  memcpy(FFTNum,FFTCash[Cache].Mem,sizeof(FFT_DATA_TYPE)*FFTLen);
#else
StartTimer(SaveTime);
StartTimer(DiskIOTime);
f=fopen(FFTCash[Cache].Name,"rb");
StopTimer(DiskIOTime);
StopTimer(SaveTime);
if (f==NULL)
   FatalError("Unable to open '%s' for reading.\n",FFTCash[Cache].Name);

StartTimer(SaveTime);
StartTimer(DiskIOTime);
x=fread(FFTNum,sizeof(FFT_DATA_TYPE),(size_t)FFTLen,f);
StopTimer(DiskIOTime);
StopTimer(SaveTime);
if ((x!=FFTLen) || (ferror(f)))
   FatalError("Error reading '%s'.\n",FFTCash[Cache].Name);
fclose(f);
#endif
}


static int
SaveFFTIntoSpecifiedCache(void *FFTNum, size_t FFTLen,
                          BigInt Num, size_t NumLen,
                          int Tag, int Pass, int Line)
/*
** Save our FFT into the specified cache line and return the cache num.
** If we can't save it, return -1.  Do NOT abort if it fails.
**
** Num  can be NO_NUM (ie: NULL)
*/
{
#ifndef VIRTUAL_CACHE
 FILE *f;
 size_t x;
#endif

  if (Line >= MaxFFTCacheLines) return -1;

#ifdef VIRTUAL_CACHE
  if (FFTCash[Line].Mem) AlignedFree(FFTCash[Line].Mem);
  FFTCash[Line].Mem=AlignedMalloc(sizeof(FFT_DATA_TYPE)*FFTLen);
  if (FFTCash[Line].Mem==NULL) return -1;
  memcpy(FFTCash[Line].Mem,FFTNum,sizeof(FFT_DATA_TYPE)*FFTLen);
#else
  StartTimer(SaveTime);
  StartTimer(DiskIOTime);
  f=fopen(FFTCash[Line].Name,"wb");
  StopTimer(DiskIOTime);
  StopTimer(SaveTime);
  if (f==NULL) return -1;

  StartTimer(SaveTime);
  StartTimer(DiskIOTime);
  x=fwrite(FFTNum,sizeof(FFT_DATA_TYPE),(size_t)FFTLen,f);
  StopTimer(SaveTime);
  StopTimer(DiskIOTime);
  if ( (x!=FFTLen) || (ferror(f)) )
    {fclose(f);remove(FFTCash[Line].Name);return -1;}
  fclose(f);
#endif

  FFTCash[Line].NumLen=NumLen;
  FFTCash[Line].FFTLen=FFTLen;
  FFTCash[Line].ChkNum=Num;
  FFTCash[Line].Tag=Tag;
  FFTCash[Line].Pass=Pass;

return Line;
}

int
SaveFFTIntoCache(void *FFTNum, size_t FFTLen,BigInt Num, size_t NumLen,
                 int Tag, int Pass)
/*
** Save our FFT into a cache and return the cache num.
** If we can't save it, return 0.  Do NOT abort if it fails.
*/
{int x;
 int Line=0;

  while (Line < MaxFFTCacheLines)
    if (FFTCash[++Line].NumLen==0) break;

  if (Line >= MaxFFTCacheLines) return 0;

  DumpDebug("Save%d (%d)...",Line,Tag);

  x=SaveFFTIntoSpecifiedCache(FFTNum,FFTLen,Num,NumLen,Tag,Pass,Line);
  if (x==-1) return 0;
  return Line;
}

void
SaveFFTIntoCache0(void *FFTNum, size_t FFTLen)
{int Saved;

 DumpDebug("Save0...");

 Saved=SaveFFTIntoSpecifiedCache(FFTNum,FFTLen,NO_NUM,0,1,0,0);
 if (Saved==0) return;  /* Success */

/*
** If we've reached here, the above save failed.
** It's time to panic!
*/
 fprintf(stderr,"\n** NOTICE **  I had to delete non-critical, performance improving\n");
 fprintf(stderr,"FFT data to be able to operate.  More disk space is recommended.\n");
 fprintf(stderr,"This is not a failure.  I just didn't have enough space to keep\n");
 fprintf(stderr,"all the optional, performance improving, caches as well as the\n");
 fprintf(stderr,"data that I do actually need.\n");
 FlushFFTCache(0); /* delete everything. */

 Saved=SaveFFTIntoSpecifiedCache(FFTNum,FFTLen,NO_NUM,0,1,0,0);
 if (Saved==0) return;  /* Success */

/*
** Even after deleting everything, we _still_ couldn't save it!
** Nothing else to do except abort the program.
*/
 printf("**** PROGRAM FAILURE ****\n");
 printf("I am unable to save a number %ld bytes long\n",(long int)FFTLen*sizeof(FFT_DATA_TYPE));
 printf("I need it to be able to multiply numbers, and since I can't save\n");
 printf("it, I can't multiply.  I tried deleting all the FFT caches, which\n");
 printf("can take a lot of disk space, but that still wasn't enough.\n\n");
 printf("You are just going to have to clean off your disk or put the\n");
 printf("caches (espc. cache #0) where I can have more space to work in.\n");
 ExitPrg(EXIT_FAILURE);
}

void
InitFFTCache(size_t Len)
{int x;
 char Str[MAX_FILENAME+8];

if (ReadCfgItem(PI_CFG_FILE,"FFT-Cache","MaxCaches",&x,Cfg_ULInteger,1)==1)
  MaxFFTCacheLines=x;

MaxFFTCacheLines++; /* Allow for 0 being the special convolution cache. */

FFTCash=(struct FFTHdr*)calloc((size_t)MaxFFTCacheLines,sizeof(struct FFTHdr));
if (FFTCash==NULL)
  FatalError("Unable to allocate memory for FFT cache.\n");

for (x=0;x<MaxFFTCacheLines;x++)
  {
   FFTCash[x].Name=(char*)malloc(MAX_FILENAME+1);
   if (FFTCash[x].Name==NULL)
     FatalError("Unable to allocate memory for FFT cache.\n");
  }

if (ReadCfgItem(PI_CFG_FILE,"Files","Convolution",FFTCash[0].Name,
               Cfg_String,MAX_FILENAME)==0)
   strcpy(FFTCash[0].Name,"/opt/Tools/PiCalc/Output/piconvl.tmp");

for (x=1;x<MaxFFTCacheLines;x++)
  {
   sprintf(Str,"Cache%d",x);
   if (ReadCfgItem(PI_CFG_FILE,"FFT-Cache",Str,FFTCash[x].Name,
                  Cfg_String,MAX_FILENAME)==0)
      sprintf(FFTCash[x].Name,"/opt/Tools/PiCalc/Output/cache%d.tmp",x);
  }

for (x=0;x<MaxFFTCacheLines;x++)
  {
#ifdef VIRTUAL_CACHE
   FFTCash[x].Mem=NULL;
#endif
   FFTCash[x].FFTLen=0;
   FFTCash[x].Tag=0;
  }

FlushFFTCache(0);
}

void
DeInitFFTCache(void)
{int x;

  FlushFFTCache(0);

  for (x=0;x<MaxFFTCacheLines;x++)
    remove(FFTCash[x].Name);

  for (x=0;x<MaxFFTCacheLines;x++)
    free(FFTCash[x].Name);

  free(FFTCash);
}

