#include "pi.h"
#include "block.h"
#include "bigint.h"
#include "bigmul.h"
#include "cache.h"

/* Just for debugging of the internal stuff
#define USE_BLOCK_MEM 1
#undef USE_DISK_NUM
#undef USE_VIRT_MEM
*/

/*
** Just so you can later scan the executable to see what type
** you compiled it for.
*/
#if  defined(USE_DISK_NUM)
char *ProgramMemoryType="Program memory type is: Disk based numbers";
#elif defined(USE_VIRT_MEM)
char *ProgramMemoryType="Program memory type is: Virtual memory";
#else
char *ProgramMemoryType="Program memory type is: Debugging Block based numbers";
#endif


/*
** A fixed size buffer of FIXEDBUF_SIZE bytes.
** It is available to every file and module that needs to buffer some
** BigInt data.
*/
void *FixedBuf=NULL;
size_t FIXEDBUF_SIZE;


#if defined(USE_BLOCK_MEM) || defined(USE_DISK_NUM)
#define USE_BLOCKS 1
#define SMALL_BLOCK_SIZE 16384
#else
#undef USE_BLOCKS
#endif

size_t CoreMemAvail;
void  *CoreMemPtr=NULL;
size_t CPU_CACHE;

static char SaveFileName[MAX_FILENAME+8];


#ifdef USE_VIRT_MEM
static BigInt VarPointers[64];
static size_t VarsAllocated=0;
#else

struct VarPkg {double Start,End;
               size_t Len;
               FILE *fp;
               char *FileName;
               UINT32 *Buffer;
               int BufStat;
              };

static struct VarPkg DiskVars[64];
static BigInt NextBIAddr=0;
static size_t VarsAllocated=0;
static size_t DiskBufSize=0;
#endif

int
CheckPhysicalMem(size_t Bytes)
{size_t pm, vm;

  vm = GetVirtualMemSize();
  pm = GetPhysicalMemSize();

  fprintf(stderr,"Phys mem left: %lu  ", (ULINT)pm);
  fprintf(stderr,"Virt mem left: %lu\n", (ULINT)vm);
  if ((pm < Bytes) && (vm >= Bytes))
    {
     fprintf(stderr,"\aWarning.  The program is wanting %lu bytes but you\n",(ULINT)Bytes);
     fprintf(stderr,"only have %lu bytes of physical memory available.\n",(ULINT)pm);
     fprintf(stderr,"You do, however, have %lu bytes of virtual memory\n",(ULINT)vm);
     fprintf(stderr,"so the allocation will succeed.  Be careful!\n");
     return 0;
    }

return 1;
}

#ifdef USE_SBRK
/* override the globally defined functions of the same name. */
#include <unistd.h>

static void *
AlignedMalloc(size_t Bytes)
{void *Buf;
Buf=sbrk(Bytes);
if ( ((INT32)Buf) == -1) Buf=NULL;
return Buf;
}

static void
AlignedFree(void *ptr)
{
/* can't free sbrk() memory. */
}
#endif

void
InitBigIntPkg(size_t Len)
{size_t CoreMemWanted;

  CoreMemWanted=CheckFFTMemReq(Len);
  CoreMemPtr=NULL;

  if (ReadCfgItem(PI_CFG_FILE,"Memory","Physical",&CoreMemAvail,Cfg_ULInteger,1)==0)
    CoreMemAvail=ULONG_MAX;
  if (CoreMemAvail < 1024*1024)
    FatalError("Come on... Surely you can give at least a meg of physical memory!\n");

  fprintf(stderr,"Mem specified: %lu  Mem wanted: %lu\n",
          (ULINT)CoreMemAvail,(ULINT)CoreMemWanted);
  CoreMemAvail=Min(CoreMemAvail,CoreMemWanted);
  do
    {
     CheckPhysicalMem(CoreMemAvail);
     CoreMemPtr=AlignedMalloc(CoreMemAvail);
     if (CoreMemPtr!=NULL) break;
     CoreMemAvail/=2;
    } while ((CoreMemPtr==NULL) && (CoreMemAvail >= 64*1024));

  if (CoreMemPtr==NULL)
    FatalError("\nUnable to allocate even %lu bytes for internal working space.\n",
               (ULINT)CoreMemAvail);

#ifdef USE_DISK_NUM
  if (ReadCfgItem(PI_CFG_FILE,"Memory","DiskBuffer",
                  &DiskBufSize,Cfg_ULInteger,1)==0)
    DiskBufSize=0;
  DiskBufSize/=8; /* allow for 8 buffers */
/*
printf("DiskBuffer size=%d\n",DiskBufSize);
*/
#endif

  if (ReadCfgItem(PI_CFG_FILE,"Memory","Cache",&CPU_CACHE,Cfg_ULInteger,1)==0)
    CPU_CACHE=128*1024;;
/*  if (CPU_CACHE==0) CPU_CACHE=CoreMemAvail;*/
  if (CPU_CACHE==0)
    {
     CPU_CACHE=DetectCacheSize(CoreMemPtr,CoreMemAvail);
     UpdateCfgItem(PI_CFG_FILE,"Memory","Cache",&CPU_CACHE,Cfg_ULInteger,1);
     fprintf(stderr,"Calculated good cache setting as %lu bytes.\n",(ULINT)CPU_CACHE);
    }

  if (ReadCfgItem(PI_CFG_FILE,"Memory","Buffer",
                  &FIXEDBUF_SIZE,Cfg_ULInteger,1)==0)
    FIXEDBUF_SIZE=128*1024;
  if (FIXEDBUF_SIZE < 32768) FIXEDBUF_SIZE=32768;
  FixedBuf=AlignedMalloc((FIXEDBUF_SIZE/sizeof(INT32))*sizeof(INT32));
  if (FixedBuf==NULL)
    FatalError("Unable to allocate %lu bytes for internal misc. buffers.\n",(ULINT)FIXEDBUF_SIZE);

  fprintf(stderr,"Allocated %lu bytes of core and %lu bytes of Buffer.\n",
         (ULINT)CoreMemAvail,(ULINT)FIXEDBUF_SIZE);
  VarsAllocated=0;

  if (ReadCfgItem(PI_CFG_FILE,"Files","SaveFile",
                  SaveFileName,Cfg_String,MAX_FILENAME)==0)
     strcpy(SaveFileName,"piceb.sav");

  InitFFT(Len);
  InitFFTMul(Len);
  InitFFTCache(Len);
}

void
DeInitBigIntPkg(void)
{int x;
  DeInitFFT();
  DeInitFFTMul();
  DeInitFFTCache();
#ifdef USE_DISK_NUM
  for (x=0;x<VarsAllocated;x++)
    {
     if (DiskVars[x].fp) fclose(DiskVars[x].fp);
     if (DiskVars[x].FileName) remove(DiskVars[x].FileName);
     if (DiskVars[x].FileName) free(DiskVars[x].FileName);
     if (DiskVars[x].Buffer) free(DiskVars[x].Buffer);
    }
#else
  for (x=0;x<VarsAllocated;x++)
    AlignedFree(VarPointers[x]);
#endif
  AlignedFree(FixedBuf);
  AlignedFree(CoreMemPtr);
}

BigInt
CreateBigInt(size_t Len)
{BigInt Num;

  if (VarsAllocated >= 50) FatalError("Too many BigInts.  I can only do 50.\n");

  if (Len==0)
    FatalError("Why are you trying to create a zero length BigInt?\n");

 if (((double)Len)*sizeof(INT32) > 2147483648.0)
   FatalError("Unable to create a BigInt that is more than 2g bytes. %lu\n",(ULINT)Len);

#ifdef USE_DISK_NUM
 {struct VarPkg *vp;
  char Str[MAX_FILENAME+8];

  vp=&DiskVars[VarsAllocated];
  vp->FileName=(char*)malloc(MAX_FILENAME+1);
  if (vp->FileName==NULL)
    FatalError("Unable to allocate memory for BigInt variable filename.\n");

  sprintf(Str,"Var%ld",VarsAllocated+1);
  if (ReadCfgItem(PI_CFG_FILE,"DiskNumbers",Str,
                  vp->FileName,Cfg_String,MAX_FILENAME)==0)
     sprintf(vp->FileName,"pivar%ld.tmp",VarsAllocated+1);

  vp->fp=fopen(vp->FileName,"wb+");
  if (vp->fp==NULL)
    FatalError("Unable to open var file %d %s\n",VarsAllocated+1,vp->FileName);

  if (ferror(vp->fp))
    FatalError("Error creating bigint.\n");

  DiskBufSize=NULL;
  if (DiskBufSize)
    vp->Buffer=(UINT32*)malloc(DiskBufSize);
/* It's okay if its NULL */
  vp->BufStat=0;
//if (vp->Buffer==NULL) printf("Var %d doesn't have a buffer.\n",VarsAllocated);

/*
** Put a 1024 element header and trailer around the address of
** the number.  This protects 'NULL' and a few addresses below
** and above what you are supposed to use.  (ie: a simple form of
** address protection.)  No actual space is allocated for this.
*/
  Num=1024+NextBIAddr;
  vp->Start=Num;
  vp->Len=Len;
  vp->End=Num+Len;
  NextBIAddr+=(1024+Len+1024);
  fseek(vp->fp,Len*sizeof(INT32)-1,SEEK_SET);fputc(' ',vp->fp);
  fflush(vp->fp);
  if (ferror(vp->fp)) FatalError("Error creating a %lu long BigInt.\n",(ULINT)Len);
  VarsAllocated++;
 }
#else
  Num=(BigInt)AlignedMalloc(Len*sizeof(INT32));
  if (Num==NULL) FatalError("Unable to create a %lu long BigInt",(ULINT)Len);
  VarPointers[VarsAllocated++]=Num;
#endif
  return Num;
}


/**************************************************
** The four routines that actually access the data.
**************************************************/

/*
** Read a block of data from a number into a buffer in main memory.
** The buffer must be large enough to hold the data.
*/
#ifdef USE_DISK_NUM

UINT32
ReadFromDiskBuffer(INT32 *Buf,BigInt Num,size_t Len)
{int x;
 int v=-1;
 size_t Ndx,ENdx;

for (x=0;x<VarsAllocated;x++)
  if ((Num >=DiskVars[x].Start) &&
      (Num < DiskVars[x].End))
    {v=x;break;}

if (v==-1) return 0;

if (DiskVars[v].Buffer==NULL) return 0;
/*
if (DiskVars[v].BufStat==0) return 0;
*/

Ndx=(size_t)(Num - DiskVars[v].Start);
if (Ndx >= (DiskBufSize/sizeof(INT32))) return 0;
ENdx=Ndx+Len;if (ENdx >= (DiskBufSize/sizeof(INT32))) ENdx=DiskBufSize/sizeof(INT32);
//printf("R: v=%d %d %d %d %d\n",v,Ndx,Len,ENdx,ENdx-Ndx);
memcpy(Buf,&DiskVars[v].Buffer[Ndx],(ENdx-Ndx)*sizeof(INT32));
return (ENdx-Ndx);
}

UINT32
SaveIntoDiskBuffer(BigInt Num,INT32 *Buf,size_t Len)
{int x;
 int v=-1;
 size_t Ndx,ENdx;

for (x=0;x<VarsAllocated;x++)
  if ((Num >=DiskVars[x].Start) &&
      (Num < DiskVars[x].End))
    {v=x;break;}

if (v==-1) return 0;

if (DiskVars[v].Buffer==NULL) return 0;
/*
if (DiskVars[v].BufStat==0) return 0;
*/

Ndx=(size_t)(Num - DiskVars[v].Start);
if (Ndx >= (DiskBufSize/sizeof(INT32))) return 0;
ENdx=Ndx+Len;if (ENdx >= (DiskBufSize/sizeof(INT32))) ENdx=DiskBufSize/sizeof(INT32);
//printf("W: v=%d %d %d %d %d\n",v,Ndx,Len,ENdx,ENdx-Ndx);
memcpy(&DiskVars[v].Buffer[Ndx],Buf,(ENdx-Ndx)*sizeof(INT32));
return (ENdx-Ndx);
}

static FILE *
PositionFile(BigInt Addr,size_t Len)
{int x;
 int v=-1;
 long int offset;

for (x=0;x<VarsAllocated;x++)
  if ((Addr >=DiskVars[x].Start) &&
      (Addr < DiskVars[x].End))
    {v=x;break;}

if (v==-1) FatalError("Disk Var %f doesn't relate to any current bigint.\n",Addr);

if (Addr+Len > DiskVars[v].End)
  FatalError("Attempting to access a BigInt beyond its allocated space.\n");

if (((((double)Addr)-DiskVars[v].Start)*sizeof(INT32)) > (double)LONG_MAX)
  FatalError("Some how, in PositionFile, the offset is way too big.\n%f %f\n",
             Addr,DiskVars[v].Start);

offset=(long int)(Addr - DiskVars[v].Start)*sizeof(INT32);
if (fseek(DiskVars[v].fp,offset,SEEK_SET) != 0)
  FatalError("Unable to fseek to position %u for reading.\n",offset);

if (ferror(DiskVars[v].fp))
  FatalError("There was an error of some sort in PositionFile. %f\n",Addr);

return DiskVars[v].fp;
}

void
ReadNumIntoBuf(BigInt Num,INT32 *Buf,size_t Len)
{size_t x;
 FILE *f;

if (Len==0) return;

x=ReadFromDiskBuffer(Buf,Num,Len);
Num+=x;
Len-=x;
Buf+=x;

if (Len==0) return;

StartTimer(DiskIOTime);
f=PositionFile(Num,Len);

x=fread(Buf,sizeof(INT32),Len,f);
if (ferror(f)) FatalError("Error reading bigint file.\n");
if (x!=Len) FatalError("Failure reading %u INT32s.  Instead: %u.\n",Len,x);

StopTimer(DiskIOTime);
}

/*
** Return a specific index'ed digit.
** Index is zero based.
*/
INT32
GetBigIntDigit(BigInt Num,size_t Ndx)
{INT32 Buf[1];
ReadNumIntoBuf(Num+Ndx,Buf,1);
return Buf[0];
}

/*
** Wite a block of data in main memory to a number.
*/
void
WriteBufIntoNum(INT32 *Buf,BigInt Num,size_t Len)
{size_t x;
 FILE *f;

if (Len==0) return;

x=SaveIntoDiskBuffer(Num,Buf,Len);
Num+=x;
Len-=x;
Buf+=x;

if (Len==0) return;

StartTimer(DiskIOTime);
f=PositionFile(Num,Len);

x=fwrite(Buf,sizeof(INT32),Len,f);
fflush(f);
if (ferror(f)) FatalError("Error writing bigint file.\n");
if (x!=Len) FatalError("Failure writing %u INT32s.  Instead: %u.\n",Len,x);

StopTimer(DiskIOTime);
}

/*
** Sets a specific index'ed digit.
** Index is zero based.
*/
void
SetBigIntDigit(BigInt Num,size_t Ndx,INT32 Val)
{INT32 Buf[1];
Buf[0]=Val;
WriteBufIntoNum(Buf,Num+Ndx,1);
}


/* virtual memory based numbers */
#else

/*
** Read a block of data from a number into a buffer in main memory.
** The buffer must be large enough to hold the data.
*/
void
ReadNumIntoBuf(BigInt Num,INT32 *Buf,size_t Len)
{
if (Len==0) return;
memcpy(Buf,Num,Len*sizeof(INT32));
}

/*
** Wite a block of data to a number.
*/
void
WriteBufIntoNum(INT32 *Buf,BigInt Num,size_t Len)
{
if (Len==0) return;
memcpy(Num,Buf,Len*sizeof(INT32));
}

/*
** Return a specific index'ed digit.
** Index is zero based.
*/
INT32
GetBigIntDigit(BigInt Num,size_t Ndx)
{
return Num[Ndx];
}

/*
** Sets a specific index'ed digit.
** Index is zero based.
*/
void
SetBigIntDigit(BigInt Num,size_t Ndx,INT32 Val)
{
Num[Ndx]=Val;
}
#endif


/***************************************
** The actual routines that do the math.
***************************************/


/*
** Does 'Num' have the 2 values as it's 2 most significant digits?
** Used in Newton to preset the starting values.
*/
int
NumIs(BigInt Num, INT32 Val1, INT32 Val2)
{
  if ((GetBigIntDigit(Num,0)==Val1) &&
      (GetBigIntDigit(Num,1)==Val2)) return 1;
  return 0;
}

/*
** Zero out all or part of a BigInt
*/
void
ClearBigInt(BigInt Num,size_t Len)
#ifdef USE_BLOCKS
{size_t BufSz=CoreMemAvail/sizeof(INT32);
 INT32 *Buf=(INT32*)CoreMemPtr;

  BlockClear(Buf,Buf+Min(Len,BufSz));
  while (Len)
    {int x;
     x=Min(BufSz,Len);
     WriteBufIntoNum(Buf,Num,x);
     Num+=x;
     Len-=x;
    }
}
#else
{
  BlockClear(Num,&Num[Len]);
}
#endif

/*
** Initialize a number to zero and explicitly set the 2 MSD.
*/
void
SetNum(BigInt Num,size_t Len, INT32 Val1, INT32 Val2)
{
  ClearBigInt(Num,Len);
  SetBigIntDigit(Num,0,Val1);
  SetBigIntDigit(Num,1,Val2);
}

/*
** Add a small integer to Num
*/
void
AddInt(BigInt Num,INT32 Val)
{
  SetBigIntDigit(Num,0,(INT32)(GetBigIntDigit(Num,0)+Val));
}

/*
** Subtract a small integer from Num.  Num = Num - int.0000
** Int is 8 digits, with high digit the integer.
*/
void
SubInt(BigInt Num, INT32 Val)
{
  SetBigIntDigit(Num,0,(INT32)(GetBigIntDigit(Num,0)-Val));
}

/*
** Yuk.  Rework.
*/
/*
** Num = int.00000 - Num
** Int is 8 digits, with high digit the integer.
*/
INT32
RevSubInt(INT32 Val,BigInt Num, size_t Len)
#ifdef USE_BLOCKS
{size_t L;
 INT32 Borrow=0;
 BigInt N;
 size_t BufSz=CoreMemAvail/sizeof(INT32);
 INT32 *Buf=(INT32*)CoreMemPtr;

  N=Num+Len;
  L=Len;
  while (L)
    {int x;
     x=Min(BufSz,L);
     N-=x;L-=x;
     ReadNumIntoBuf(N,Buf,x);
     if (L) Borrow=BlockNegate(Buf,Borrow,x,0);
     else   Borrow=BlockNegate(Buf,Borrow,x,Val);
     WriteBufIntoNum(Buf,N,x);
    }
  if (Borrow)
    {
     N=Num+Len;
     L=Len;
     Borrow=0;
     while (L)
       {size_t x;
        x=Min(BufSz,L);
        N-=x;L-=x;
        ReadNumIntoBuf(N,Buf,x);
        Borrow=BlockNegate(Buf,Borrow,x,0);
        WriteBufIntoNum(Buf,N,x);
       }
     return 1;
    }
  return 0;
}
#else
{INT32 Sign;
  Sign = BlockNegate(Num,0,Len,Val);
  if (Sign) BlockNegate(Num,0,Len,0);
  return Sign;
}
#endif

/*
** Num1=Num2
*/
void
Copy(BigInt Dest, BigInt Src, size_t Len)
#ifdef USE_BLOCKS
{size_t BufSz=CoreMemAvail/sizeof(INT32);
 INT32 *Buf=(INT32*)CoreMemPtr;
if ((Src < Dest) && (Dest < Src+Len))
  {
   Src+=Len;Dest+=Len;
   while (Len)
     {size_t x;
      x=Min(BufSz,Len);
      Src-=x;Dest-=x;
      ReadNumIntoBuf(Src,Buf,x);
      WriteBufIntoNum(Buf,Dest,x);
      Len-=x;
     }
  }
else
  while (Len)
    {size_t x;
     x=Min(BufSz,Len);
     ReadNumIntoBuf(Src,Buf,x);
     WriteBufIntoNum(Buf,Dest,x);
     Src+=x;Dest+=x;
     Len-=x;
    }
}
#else
{
  BlockCopy(Dest,Src,Len);
}
#endif

/*
** Is 'Num' exactly zero?
*/
int
IsZero(BigInt Num, size_t Len)
#ifdef USE_BLOCKS
{size_t BufSz=CoreMemAvail/sizeof(INT32);
 INT32 *Buf=(INT32*)CoreMemPtr;

  while (Len)
    {size_t x;
     x=Min(Min(BufSz,SMALL_BLOCK_SIZE),Len);
     ReadNumIntoBuf(Num,Buf,x);
     if (BlockCountZeros(Buf,x) != x) return 0;
     Num+=x;
     Len-=x;
    }
return 1;
}
#else
{
  return BlockCountZeros(Num,Len) == Len;
}
#endif

/*
** Find first non zero digit, else return the Length.
*/
size_t
FindFirstNonZero(BigInt Num, size_t Len)
#ifdef USE_BLOCKS
{size_t First;
 size_t BufSz=CoreMemAvail/sizeof(INT32);
 INT32 *Buf=(INT32*)CoreMemPtr;

  First=0;

  while (Len)
    {size_t x,y;
     x=Min(Min(BufSz,SMALL_BLOCK_SIZE),Len);
     ReadNumIntoBuf(Num,Buf,x);
     for (y=0;y<x;y++)
       if (Buf[y]) return y+First;
     Num+=x;
     Len-=x;
     First+=x;
    }
return First;
}
#else
{size_t x;
 for (x=0;x<Len;x++) if (Num[x]) return x;
 return Len;
}
#endif

/*
** Count the number of zeros in the number.  Knowing that can let
** us do a special type of multiplication, which is faster than
** the normal one.
**
** Left is the most significant half, right is the lest significant
** half.  If those variables are NULL, the results wont be returned.
** The return value itself is the total number of zeros.
*/
size_t
CountZeros(BigInt Num, size_t *Left,size_t *Right,size_t Len)
#ifdef USE_BLOCKS
{size_t Half=Len/2;
 size_t LenL=Half,LenR=Half;
 size_t L=0,R=0;
 size_t BufSz=CoreMemAvail/sizeof(INT32);
 INT32 *Buf=(INT32*)CoreMemPtr;

  while (LenL)
    {size_t x;
     x=Min(BufSz,LenL);
     ReadNumIntoBuf(Num,Buf,x);
     L+=BlockCountZeros(Buf,x);
     Num+=x;LenL-=x;
    }
  while (LenR)
    {size_t x;
     x=Min(BufSz,LenR);
     ReadNumIntoBuf(Num,Buf,x);
     R+=BlockCountZeros(Buf,x);
     Num+=x;LenR-=x;
    }
  if (Left)  *Left=L;
  if (Right) *Right=R;
  return L+R;
}
#else
{size_t Half=Len/2;
 size_t L,R;

  L = BlockCountZeros(&Num[0],Half);
  R = BlockCountZeros(&Num[Half],Half);
  if (Left)  *Left=L;
  if (Right) *Right=R;
  return L+R;
}
#endif

/*
** Result = Num / Val
*/
void
DivBy(BigInt Result, BigInt Num,INT32 Val, size_t Len)
#ifdef USE_BLOCKS
{INT32 Remain=0;
 size_t BufSz=CoreMemAvail/sizeof(INT32);
 INT32 *Buf=(INT32*)CoreMemPtr;

  while (Len)
    {size_t x;
     x=Min(BufSz,Len);
     ReadNumIntoBuf(Num,Buf,x);
     Remain=BlockDivBy(Buf,Buf,Val,Remain,x);
     WriteBufIntoNum(Buf,Result,x);
     Num+=x;
     Result+=x;
     Len-=x;
    }
}
#else
{
  BlockDivBy(Result,Num,Val,0,Len);
}
#endif

/*
** Num=Num*dig
**
** Simple multiplication of 'Num' by any value.  It doesn't do
** shift or overflow, etc.  It's just a simple, basic single
** value multiplication routine.  It trusts you to know when
** to properly use it.
**
** It's a bit slow (since it uses floating point), so it's not
** something you should use very often.  It's only really useful
** in the places where we need to multiply by a power of two.
*/
double
MulByFloat(BigInt Num, double Val, size_t Len)
#ifdef USE_BLOCKS
{double Carry=0.0;
 size_t BufSz=CoreMemAvail/sizeof(INT32);
 INT32 *Buf=(INT32*)CoreMemPtr;

  Num+=Len;
  while (Len)
    {size_t x;
     x=Min(BufSz,Len);
     Num-=x;
     ReadNumIntoBuf(Num,Buf,x);
     Carry=BlockMulByFloat(Buf,Val,Carry,x);
     WriteBufIntoNum(Buf,Num,x);
     Len-=x;
    }
return Carry;
}
#else
{
return BlockMulByFloat(Num,Val,0.0,Len);
}
#endif

/*
** Num=Num*dig
**
** Simple multiplication of 'Num' by any value.  It doesn't do
** shift or overflow, etc.  It's just a simple, basic single
** value multiplication routine.  It trusts you to know when
** to properly use it.
*/
INT32
MulBy(BigInt Result,BigInt Num, INT32 Val, size_t Len)
#ifdef USE_BLOCKS
{INT32 Carry=0;
 size_t BufSz=CoreMemAvail/sizeof(INT32);
 INT32 *Buf=(INT32*)CoreMemPtr;

  Num+=Len;
  Result+=Len;
  while (Len)
    {size_t x;
     x=Min(BufSz,Len);
     Num-=x;Result-=x;
     Len-=x;
     ReadNumIntoBuf(Num,Buf,x);
     Carry=BlockMulBy(Buf,Buf,Val,Carry,x);
     WriteBufIntoNum(Buf,Result,x);
    }
return Carry;
}
#else
{INT32 Carry;
Carry=BlockMulBy(Result,Num,Val,0,Len);
return Carry;
}
#endif

/*
** A BigInt library support routine to ripple a carry up
** a number.
*/
static INT32
RippleCarry(BigInt Limit,BigInt Sum,INT32 Carry)
#ifdef USE_BLOCKS
{size_t BufSz=CoreMemAvail/sizeof(INT32);
 INT32 *Buf=(INT32*)CoreMemPtr;

  while ((Sum > Limit) && (Carry))
    {size_t Len,x;
      Len = Min(Min(BufSz,SMALL_BLOCK_SIZE),(size_t)(Sum-Limit));
      Sum-=Len;
      ReadNumIntoBuf(Sum,Buf,Len);
      x=Len;
      while ((x) && (Carry))
        {
         x--;
         Carry = (INT32)(Buf[x] + Carry);
         Buf[x] = (INT32)(Carry % BaseVar);
         Carry  = (INT32)(Carry / BaseVar);
        }
      WriteBufIntoNum(Buf,Sum,Len);
    }
return Carry;
}
#else
{
  while ((Sum > Limit) && (Carry))
    {
      Carry = (INT32)(*(--Sum) + Carry);
      *Sum  = (INT32)(Carry % BaseVar);
      Carry = (INT32)(Carry / BaseVar);
    }
  return Carry;
}
#endif

/*
** Sum = Num1 + Num2
**
** Very much like the regular ADD routine, except the carries
** can ripple on up beyond where we end our addition.
*/
INT32
RippleAdd(BigInt Limit,BigInt Sum,BigInt Num1,BigInt Num2,size_t Len)
#ifdef USE_BLOCKS
{INT32 Carry=0;
 size_t BufSz=CoreMemAvail/sizeof(INT32);
 INT32 *Buf1=(INT32*)CoreMemPtr;
 INT32 *Buf2=NULL;

  if (Num1!=Num2) {BufSz/=2;Buf2=Buf1+BufSz;}
  Num1+=Len;
  Num2+=Len;
  Sum+=Len;
  while (Len)
    {size_t x;
     x=Min(BufSz,Len);
     Num1-=x;Num2-=x;Sum-=x;
     ReadNumIntoBuf(Num1,Buf1,x);
     if (Num1==Num2)                   Carry = BlockAdd(Buf1,Buf1,Buf1,Carry,x);
     else {ReadNumIntoBuf(Num2,Buf2,x);Carry = BlockAdd(Buf1,Buf1,Buf2,Carry,x);}
     WriteBufIntoNum(Buf1,Sum,x);
     Len-=x;
    }
  Carry = RippleCarry(Limit,Sum,Carry);
  return Carry;
}
#else
{INT32 Carry;
  Carry = BlockAdd(Sum,Num1,Num2,0,Len);
  Carry = RippleCarry(Limit,Sum,Carry);
  return Carry;
}
#endif

/*
** Sum = Num1 + Num2
**
** Returns whether a carry occured off the end.
*/
INT32
Add(BigInt Sum, BigInt Num1, BigInt Num2, size_t Len)
{
  return RippleAdd(Sum,Sum,Num1,Num2,Len);
}

/*
** Dif = Num1 - Num2
** Just like a regular sub, except the borrow can ripple on up higher,
** beyond where we end our subtraction.
*/
INT32
RippleSub(BigInt Limit,BigInt Dif,BigInt Num1,BigInt Num2,size_t Len)
#ifdef USE_BLOCKS
{INT32 Borrow=0;
 size_t BufSz=CoreMemAvail/sizeof(INT32);
 INT32 *Buf1=(INT32*)CoreMemPtr;
 INT32 *Buf2=NULL;

  if (Num1!=Num2) {BufSz/=2;Buf2=Buf1+BufSz;}
  Num1+=Len;
  Num2+=Len;
  Dif+=Len;
  while (Len)
    {size_t x;
     x=Min(BufSz,Len);
     Num1-=x;Num2-=x;Dif-=x;
     ReadNumIntoBuf(Num1,Buf1,x);
     if (Num1==Num2)                   Borrow = BlockSub(Buf1,Buf1,Buf1,Borrow,x);
     else {ReadNumIntoBuf(Num2,Buf2,x);Borrow = BlockSub(Buf1,Buf1,Buf2,Borrow,x);}
     WriteBufIntoNum(Buf1,Dif,x);
     Len-=x;
    }

  while ((Dif > Limit) && (Borrow))
    {size_t Len,x;
      Len = Min(Min(BufSz,SMALL_BLOCK_SIZE),Dif-Limit);
      Dif-=Len;
      ReadNumIntoBuf(Dif,Buf1,Len);
      x=Len;
      while (x)
        {INT32 temp;
         x--;
         temp = (INT32)(Buf1[x] - Borrow);
         Borrow = 0;
         if (temp < 0) { Borrow = 1; temp =(INT32)(temp + BaseVar); }
         Buf1[x] = temp;
        }
      WriteBufIntoNum(Buf1,Dif,Len);
    }
  return Borrow;
}
#else
{INT32 temp, Borrow;

  Borrow = BlockSub(Dif,Num1,Num2,0,Len);

  while ((Dif > Limit) && (Borrow))
    {
      temp = (INT32)(*(--Dif) - Borrow);
      Borrow = 0;
      if (temp < 0) { Borrow = 1; temp =(INT32)(temp+BaseVar); }
      *Dif = temp;
    }

  return Borrow;
}
#endif

/*
** Dif = Num1 - Num2
**
** Returns whether a borrow occurred off the end.
*/
INT32
Sub(BigInt Dif, BigInt Num1, BigInt Num2, size_t Len)
{
  return RippleSub(Dif,Dif,Num1,Num2,Len);
}

/*
** Num = 0 - Num
*/
void
Negate(BigInt Num, size_t Len)
#ifdef USE_BLOCKS
{INT32 Borrow=0;
 size_t BufSz=CoreMemAvail/sizeof(INT32);
 INT32 *Buf=(INT32*)CoreMemPtr;

  Num+=Len;
  while (Len)
    {size_t x;
     x=Min(BufSz,Len);
     Num-=x;
     ReadNumIntoBuf(Num,Buf,x);
     Borrow=BlockNegate(Buf,Borrow,x,0);
     WriteBufIntoNum(Buf,Num,x);
     Len-=x;
    }
}
#else
{
  BlockNegate(Num,0,Len,0);
}
#endif

/*
** Special AGM formula.
**
** When I have enough memory to hold two full variables at once,
** I can easily, and efficiently, do the formula:
**
** AGM_B2=(AGM_B2+AGM_A2+4*AGM_C2)/2
**
** When I don't, I can do it the slow, individual way.
*/
int
SpecialAGMFunc1(BigInt AGM_B2,BigInt AGM_A2,BigInt AGM_C2, size_t Len)
{int Sign;
#ifdef USE_BLOCKS
if (Len*2*sizeof(INT32) <= CoreMemAvail)
  {size_t BufSz=CoreMemAvail/(sizeof(INT32)*2);
   INT32 *Buf1=(INT32*)CoreMemPtr;
   INT32 *Buf2=Buf1+BufSz;

   ReadNumIntoBuf(AGM_B2,Buf1,Len);
   ReadNumIntoBuf(AGM_A2,Buf2,Len);
   BlockAdd(Buf1,Buf1,Buf2,0,Len);
   ReadNumIntoBuf(AGM_C2,Buf2,Len);
   BlockMulBy(Buf2,Buf2,4,0,Len);
   Sign=BlockSub(Buf1,Buf1,Buf2,0,Len);
   BlockDivBy(Buf1,Buf1,2,0,Len);
   WriteBufIntoNum(Buf1,AGM_B2,Len);
  }
else
#endif
  {
   Add(AGM_B2,AGM_B2,AGM_A2,Len);
   MulBy(AGM_A2,AGM_C2,4,Len);
   Sign=Sub(AGM_B2,AGM_B2,AGM_A2,Len);
   DivBy(AGM_B2,AGM_B2,2,Len);
  }
return Sign;
}

/*
** Do Dif=(Num1-Num2)/2
**
** If we have enough core memory, we can do it efficiently.
** If not, just call the regular sub and div.
*/
INT32
HalfDiff(BigInt Dif, BigInt Num1, BigInt Num2, size_t Len)
{INT32 Sign;
#ifdef USE_BLOCKS
if (Len*2*sizeof(INT32) <= CoreMemAvail)
  {size_t BufSz=CoreMemAvail/(sizeof(INT32)*2);
   INT32 *Buf1=(INT32*)CoreMemPtr;
   INT32 *Buf2=Buf1+BufSz;

   ReadNumIntoBuf(Num1,Buf1,Len);
   ReadNumIntoBuf(Num2,Buf2,Len);
   Sign=BlockSub(Buf1,Buf1,Buf2,0,Len);
   BlockDivBy(Buf1,Buf1,2,0,Len);
   WriteBufIntoNum(Buf1,Dif,Len);
  }
else
#endif
  {
   Sign=Sub(Dif,Num1,Num2,Len);
   DivBy(Dif,Dif,2,Len);
  }
return Sign;
}

/*
** Plain old fashioned full multiplication.
*/
void
FullMul(BigInt Prod, BigInt Num1, BigInt Num2, size_t Len)
{
  StartTimer(MulTime);
  FFTMul(Prod,Num1,Num2,Len,Len,10);
  StopTimer(MulTime);
  Num1IsCached=Num2IsCached=SaveNum1FFT=SaveNum2FFT=0;
}

void
SpecialSquare(BigInt Prod, BigInt Num, size_t Len, BigInt Work)
/*
** Since we are squaring, we can save a tinsy bit of time
** A half sized FFT cots 98 units.  A full size costs 209 units.
** (those numbers came from actual old fft measurments.)
**
** A regular square would cost 2*209=418.
** Doing it this way, we are doing 4 half sized, 4*98=392.
** Allowing for the other overhead, that's a _slight_ savings.
** This routine can only be safely used in a few places.
*/
{size_t Half=Len/2;

  StartTimer(MulTime);
  if (Len <= 64)
      FullMul(Prod,Num,Num,Len);
  else
    {
/* Work= Upper(Num) * Lower(Num) */
     SaveNum1FFT=1;
     FFTMul(Work,Num,Num+Half,Half,Half,0);
     Num1IsCached=1;
/* Prod=Upper(Num)^2 */
     FFTMul(Prod,Num,Num,Half,Len,0);
/* Prod=Prod + 2*Shifted(Work) */
     RippleAdd(Prod,Prod+Half,Prod+Half,Work,Half);
     RippleAdd(Prod,Prod+Half,Prod+Half,Work,Half);
     MulBy(Prod,Prod,10,Len);
     DeleteFFTCache(1);
    }
  StopTimer(MulTime);
  Num1IsCached=Num2IsCached=SaveNum1FFT=SaveNum2FFT=0;
}

void
SpecialFullMul(BigInt Prod, BigInt Num1, BigInt Num2, size_t Len, BigInt Work)
/*
** This FullMul() is specifically designed to multiply two
** full sized numbers using only half sized multiplies.  This
** does, of course, take more time.  The only reason we do this
** is because in the AGM, we need to do one full sized multiplication
** while everywhere else in the program can get by with only half
** sized muls.  That one full sized mul reduces the total capacity of
** the multiplication routine.  By breaking it into half sized muls,
** we can let the program go higher than it otherwise would be
** capable of.
**
** Also, by cheating and doing only three of the four muls, we save
** save time, but we also let a little bit of inaccuracy creep in.
** In this case, though, it's okay.
**
** A half sized FFT cots 98 units.  A full size costs 209 units.
** (those numbers came from actual old fft measurments.)
*/
{size_t Half=Len/2;

  StartTimer(MulTime);
  if (Len <= 64)
      FullMul(Prod,Num1,Num2,Len);
  else
    {
     SaveNum1FFT=1;SaveNum2FFT=2;
     FFTMul(Prod,Num1,Num2,Half,Len,0);
     Num1IsCached=1;
     FFTMul(Work,Num1,Num2+Half,Half,Len,0);
     RippleAdd(Prod,Prod+Half,Prod+Half,Work,Half);
     Num2IsCached=2;
     FFTMul(Work,Num1+Half,Num2,Half,Len,0);
     RippleAdd(Prod,Prod+Half,Prod+Half,Work,Half);
     MulBy(Prod,Prod,10,Len);
     DeleteFFTCache(1);
     DeleteFFTCache(2);
    }
  StopTimer(MulTime);
  Num1IsCached=Num2IsCached=SaveNum1FFT=SaveNum2FFT=0;
}


static void
HalfMul0(BigInt Prod, BigInt Num1, BigInt Num2, size_t Len,INT32 Scaler)
{size_t z1,z2;
 size_t Half=Len/2;

  if (Len <= 64)
     {
      StartTimer(MulTime);
      FFTMul(Prod,Num1,Num2,Len,Len,Scaler);
      StopTimer(MulTime);
     }
  else
    {INT32 Carry;
     StartTimer(MulTime);
     z2=FindFirstNonZero(Num2,Half);
     if (Num1==Num2) z1=z2;
     else z1=FindFirstNonZero(Num1,Half);

     Carry=FFTMul(Prod+z1+z2,Num1+z1,Num2+z2,Half,Len-z1-z2,Scaler);
     ClearBigInt(Prod,z1+z2);
     if ((z1+z2) && (Carry))
       SetBigIntDigit(Prod,z1+z2-1,Carry);
     StopTimer(MulTime);
    }

  Num1IsCached=Num2IsCached=SaveNum1FFT=SaveNum2FFT=0;
}

/*
** Multiply two half precision numbers together and get a full
** precision result.  Just an optimization for the division and
** root routines.  Part of "Karp's trick".
*/
void
HalfMul(BigInt Prod, BigInt Num1, BigInt Num2, size_t Len)
{
HalfMul0(Prod,Num1,Num2,Len,10);
}

/*
** Like regular HalfMul(), except we also divide by 2.  Used
** in a few places and it's faster to work it in like this,
** rather than doing additional disk load & stores.
*/
void
HalfMul2(BigInt Prod, BigInt Num1, BigInt Num2, size_t Len)
{
HalfMul0(Prod,Num1,Num2,Len,5);
}


/*
** Right half of Num1 is zero.  Not much of an optimization, but I
** might as well do it.
*/
void
N1R0Mul(BigInt Prod, BigInt Num1, BigInt Num2, BigInt Work,size_t Len)
{size_t Half=Len/2;

  if (Len <= 64)
     {
      StartTimer(MulTime);
      FFTMul(Prod,Num1,Num2,Len,Len,10);
      StopTimer(MulTime);
      Num1IsCached=Num2IsCached=SaveNum1FFT=SaveNum2FFT=0;
      return;
     }

  SaveNum1FFT=1;
  StartTimer(MulTime);
  FFTMul(Work,Num1,Num2+Half,Half,Len,0);
  Num1IsCached=1;
  FFTMul(Prod,Num1,Num2     ,Half,Len,0);
  RippleAdd(Prod,Prod+Half,Prod+Half,Work,Half);
  MulBy(Prod,Prod,10,Len);
  StopTimer(MulTime);
  DeleteFFTCache(1);
  Num1IsCached=Num2IsCached=SaveNum1FFT=SaveNum2FFT=0;
}


/*
** Provide for some CRC-32 checking.
*/
static UINT32 CRCTable[256];

static void
BuildCRCTable(void)
{int i,j;UINT32 crc;
for (i=0;i<256;i++)
  {
   crc=i;
   for (j=8;j>0;j--)
     {
      if (crc & 1) crc=(crc >> 1) ^ 0xEDB88320UL;
      else         crc >>= 1;
     }
   CRCTable[i]=crc;
  }
}

static UINT32
CalcCRC(size_t count, UINT32 CRC, INT32 *Buf)
{UINT32 t1,t2;
if (Cfg.UseCRC)
  while (count)
    {unsigned char Byte;
     Byte=(unsigned char)((*Buf >> 24) & 0xff);
     t1=(CRC >> 8) & 0x00FFFFFFUL;
     t2=CRCTable[(CRC ^ Byte) & 0xff];
     CRC=t1 ^ t2;

     Byte=(unsigned char)((*Buf >> 16) & 0xff);
     t1=(CRC >> 8) & 0x00FFFFFFUL;
     t2=CRCTable[(CRC ^ Byte) & 0xff];
     CRC=t1 ^ t2;

     Byte=(unsigned char)((*Buf >> 8) & 0xff);
     t1=(CRC >> 8) & 0x00FFFFFFUL;
     t2=CRCTable[(CRC ^ Byte) & 0xff];
     CRC=t1 ^ t2;

     Byte=(unsigned char)((*Buf) & 0xff);
     t1=(CRC >> 8) & 0x00FFFFFFUL;
     t2=CRCTable[(CRC ^ Byte) & 0xff];
     CRC=t1 ^ t2;

     count--;Buf++;
    }
return CRC;
}

/***************************
** Do raw disk BigInt I/O **
***************************/

void
DeleteSaveFile(void)
{
remove(SaveFileName);
}

static void
LoadOneVar(FILE *f,BigInt Var,size_t Len,int Which)
{UINT32 CRC=0xFFFFFFFFUL;

#ifdef USE_BLOCKS
 INT32 *Buf=(INT32*)CoreMemPtr;
 size_t BufSize=CoreMemAvail/sizeof(INT32);

if (Var)
  {size_t NLen=Len;size_t x,y;
    while (NLen)
      {x=Min(BufSize,NLen);
       y=fread(Buf,sizeof(INT32),x,f);
       if (y!=x) FatalError("fread failure in LoadOneVar.\n");
       WriteBufIntoNum(Buf,Var,x);
       CRC=CalcCRC(x,CRC,Buf);
       Var+=x;NLen-=x;
      }
   if (Cfg.UseCRC)
     fprintf(stderr,"Var %u crc32=%08x\n",Which,~CRC);
  }
#else
if (Var!=NO_NUM)
  {size_t y;
   y=fread(Var,sizeof(INT32),Len,f);
   if (y!=Len) FatalError("fread failure in LoadOneVar.\n");
   CRC=CalcCRC(Len,CRC,Var);
   if (Cfg.UseCRC)
     fprintf(stderr,"Var %u crc32=%08x\n",Which,~CRC);
  }
#endif
}

/*
** Try loading the 'save data' file.  If it doesn't exist,
** just return with a '0' so we can know that we didn't load
** anything.  If an error occurs during loading, abort the
** program.  If everything loads okay, return a '1'.
*/
int
LoadData(BigInt Var1, BigInt Var2, BigInt Var3, BigInt Var4,
         BigInt Var5, BigInt Var6,
         clock_t *StartTime, double *Pow2, size_t *Pass, size_t Len,
         size_t Algor)
{FILE *f;
 size_t x;
 double elap;
 char str[80];

BuildCRCTable();
f=fopen(SaveFileName,"rb");
if (f==NULL) return 0;

fprintf(stderr,"Loading vars");
fread(str,strlen(PROG_VERSION),1,f);
if (memcmp(str,PROG_VERSION,strlen(PROG_VERSION))!=0)
   FatalError("Wrong save file version.\n");

fread(&x,sizeof(size_t),1,f);
if (x!=Algor) FatalError("Wrong algorithm.  The save is for %u\n",x);

fread(&x,sizeof(size_t),1,f);
if (x!=Len) FatalError("The save is for %u digits, not %u.\n",
                       x*RawIntDigits,Len*RawIntDigits);

fread(&elap,sizeof(double),1,f);
fread(Pow2,sizeof(double),1,f);
fread(Pass,sizeof(size_t),1,f);

LoadOneVar(f,Var1,Len,1);
LoadOneVar(f,Var2,Len,2);
LoadOneVar(f,Var3,Len,3);
LoadOneVar(f,Var4,Len,4);
LoadOneVar(f,Var5,Len,5);
LoadOneVar(f,Var6,Len,6);

if (ferror(f)!=0)
  {
   fclose(f);
   FatalError("Error loading data from %s.\n",SaveFileName);
  }

BackSpace(12);
fclose(f);
if (Cfg.DeleteSaveFile)
  DeleteSaveFile();
*StartTime=elap*(float)CLOCKS_PER_SEC;
printf("Time: %li\n",StartTime);
return 1;
}

static void
SaveOneVar(FILE *f,BigInt Var,size_t Len,int Which)
{UINT32 CRC=0xFFFFFFFFUL;

#ifdef USE_BLOCKS
 INT32 *Buf=(INT32*)CoreMemPtr;
 size_t BufSize=CoreMemAvail/sizeof(INT32);

if (Var)
  {size_t NLen=Len;size_t x,y;
    while (NLen)
      {x=Min(BufSize,NLen);
       ReadNumIntoBuf(Var,Buf,x);
       y=fwrite(Buf,sizeof(INT32),x,f);
       if (y!=x) FatalError("Failure writing in SaveOneVar.\n");
       CRC=CalcCRC(x,CRC,Buf);
       Var+=x;NLen-=x;
      }
   if (Cfg.UseCRC)
     fprintf(stderr,"Var %u crc32=%08x\n",Which,~CRC);
  }
#else
if (Var!=NO_NUM)
  {size_t y;
   y=fwrite(Var,sizeof(INT32),Len,f);
   if (y!=Len) FatalError("Failure writing in SaveOneVar.\n");
   CRC=CalcCRC(Len,CRC,Var);
   if (Cfg.UseCRC)
     fprintf(stderr,"Var %u crc32=%08x\n",Which,~CRC);
  }
#endif
}

/*
** Save all the data that we need to later continue our
** computation.
*/
void
SaveData(BigInt Var1, BigInt Var2, BigInt Var3, BigInt Var4,
         BigInt Var5, BigInt Var6,
         double StartTime, double Pow2, size_t Pass, size_t Len,
         size_t Algor)
{FILE *f;
 double elap;

fprintf(stderr,"Saving vars");
BuildCRCTable();
DeleteSaveFile();
f=fopen(SaveFileName,"wb");
if (f==NULL)
   FatalError("Unable to open %s for saving data.\n",SaveFileName);

elap=StartTime;
fwrite(PROG_VERSION,strlen(PROG_VERSION),1,f);
fwrite(&Algor,sizeof(size_t),1,f);
fwrite(&Len,sizeof(size_t),1,f);
fwrite(&elap,sizeof(double),1,f);
fwrite(&Pow2,sizeof(double),1,f);
fwrite(&Pass,sizeof(size_t),1,f);
SaveOneVar(f,Var1,Len,1);
SaveOneVar(f,Var2,Len,2);
SaveOneVar(f,Var3,Len,3);
SaveOneVar(f,Var4,Len,4);
SaveOneVar(f,Var5,Len,5);
SaveOneVar(f,Var6,Len,6);
if (ferror(f)!=0)
  {
   fclose(f);
   FatalError("Error saving data to %s.\n",SaveFileName);
  }

BackSpace(11);
fclose(f);
}


/*********************************
** Do formatted output of BigInt's
*********************************/


/*
** Print out a simple straight representation of our big numbers
** Used only in debugging sessions.
*/
void
DumpBigInt(char *Str, BigInt Num, size_t Len)
{size_t x;

/*  printf("%10.10s : ", Str);*/
  printf("%s : ", Str);
  printf("%u.", GetBigIntDigit(Num,0)/10000000);
  printf("%07u",GetBigIntDigit(Num,0)%10000000);
  for (x = 1; x < Len/2; x++)
    printf("%08u", GetBigIntDigit(Num,x));
  printf(":");
  for (x = Len/2; x < Len; x++)
    printf("%08u", GetBigIntDigit(Num,x));
  printf("\n");
}

/*
** Print out a formated representation of our big numbers
*/
static FILE *OutPutFile;/*=stdout;*/

/* Local function for PrintFormattedPi() */
static void
PrintDig(int Dig,INT32 printed)
{
if (Cfg.OutputFormat==1)
  {
   /* My old Guttenberg output format. */
   fprintf(OutPutFile,"%d",Dig);
   if ((printed % 1000) == 0)    fprintf(OutPutFile,"\n\n\n\n");
   else if ((printed % 50) == 0) fprintf(OutPutFile,"\n");
   else if ((printed % 10) == 0) fprintf(OutPutFile," ");
  }
else if (Cfg.OutputFormat==2)
  {/* SuperPi */
   fprintf(OutPutFile,"%d",Dig);
   if ((printed % 1000) == 0)    fprintf(OutPutFile,"\n\n");
   else if ((printed % 50) == 0) fprintf(OutPutFile,"\n");
   else if ((printed % 10) == 0) fprintf(OutPutFile," ");
  }
else /*if (Cfg.OutputFormat==3) Default to this format */
  {/* Japanese archive format */
   fprintf(OutPutFile,"%d",Dig);
   if ((printed % 100) == 0) fprintf(OutPutFile,"\n");
   else if ((printed % 10) == 0) fprintf(OutPutFile," ");
  }
}

static void
DumpUnformattedPi(char *Str, double ETime, BigInt Num, size_t Len)
{INT32 *Buf=(INT32*)FixedBuf;
 char Mask[MAX_FILENAME+1];

  OutPutFile=stdout;
  if (Cfg.SavePiToFile)
    {char str[80];

    if (ReadCfgItem(PI_CFG_FILE,"Files","Pi_Outfile_Mask",
                    Mask,Cfg_String,MAX_FILENAME)==0)
       strcpy(Mask,"pi%s.txt");

     sprintf(str,Mask,Num2Str(Len*RawIntDigits));
     OutPutFile=fopen(str,"wb");
     if (OutPutFile==NULL)
       FatalError("Unable to open file \'%s\' for saving output.\n",str);
    }

  while (Len)
    {size_t x,y;
     x=Min(FIXEDBUF_SIZE/sizeof(INT32),Len);
     ReadNumIntoBuf(Num,Buf,x);
     Num+=x;Len-=x;
     BlockUnpack2(Buf,x);
/* Make it BCD
     {char  *CBuf=(char*)FixedBuf;
      for (y=0;y<x*4;y++) CBuf[y]=(CBuf[y]/10)*16+(CBuf[y]%10);
     }
*/
     y=fwrite(Buf,sizeof(char),x*4,OutPutFile);
     if (y != x*4)
       FatalError("Error during Dumping of pi.\n");
    }

  if (Cfg.SavePiToFile)
    fclose(OutPutFile);
}

void
PrintFormattedPi(char *Str, double ETime, BigInt Num, size_t Len)
{INT32 Dig;
 int printed = 0;
 INT32 *Buf=(INT32*)FixedBuf;
 int TotalDigits=Len*RawIntDigits;
 char Mask[MAX_FILENAME+1];

  OutPutFile=stdout;
  if (Cfg.OutputFormat==0)
    {
     DumpUnformattedPi(Str,ETime,Num,Len);
     return;
    }

  if (Cfg.SavePiToFile)
    {char str[80];
    if (ReadCfgItem(PI_CFG_FILE,"Files","Pi_Outfile_Mask",
                    Mask,Cfg_String,MAX_FILENAME)==0)
       strcpy(Mask,"pi%s.txt");

     sprintf(str,Mask,Num2Str(Len*RawIntDigits));
     OutPutFile=fopen(str,"w");
     if (OutPutFile==NULL)
       FatalError("Unable to open file \'%s\' for saving output.\n",str);
    }

  while (Len)
    {size_t x,y;
     x=Min(FIXEDBUF_SIZE/sizeof(INT32),Len);
     ReadNumIntoBuf(Num,Buf,x);
     Num+=x;Len-=x;
     for (y=0;y<x;y++)
       {
        Dig=Buf[y];
        PrintDig((Dig/10000000)%10,++printed);
        if (printed==1) {printed=0;fprintf(OutPutFile,".\n\n");}
        PrintDig((Dig/1000000)%10,++printed);
        PrintDig((Dig/100000)%10,++printed);
        PrintDig((Dig/10000)%10,++printed);
        PrintDig((Dig/1000)%10,++printed);
        PrintDig((Dig/100)%10,++printed);
        PrintDig((Dig/10)%10,++printed);
        PrintDig(Dig%10,++printed);
       }
    }

  fprintf(OutPutFile,"\n\n%s computed %u digits in %0.0f seconds.\n\n", Str,TotalDigits,ETime);

  if (Cfg.SavePiToFile)
    fclose(OutPutFile);
}


