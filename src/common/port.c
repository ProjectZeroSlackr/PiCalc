#include "pi.h"

/*
** Test the keyboard to see if any key has been pressed (0=no, 1=yes)
** Used to allow the user to exit a pi calculation with the state saved.
*/
#ifdef __DJGPP__
#include <pc.h>
int TestKeyboard(void)
{
return kbhit();
}

#else

int TestKeyboard(void)
{
return 0;
}
#endif


#if defined(__DJGPP__) && defined(USE_DISK_NUM)
/*
** DOS specific.  Just to warn the user if disk write caching
** is enabled.
*/

#include <dos.h>

union REGS r;

void CheckDiskCache(void)
{int Drive,WC;

WC=0;

r.x.ax=0x4a10;
r.x.bx=0;
r.x.cx=0xebab;
int86(0x2f,&r,&r);

if (r.x.ax != 0xbabe) return;

for (Drive=0;Drive<10;Drive++)
  {
   r.x.ax=0x4a10;
   r.x.bx=3;
   r.h.bp=Drive;
   r.h.dl=0;
   int86(0x2f,&r,&r);
   if ((r.x.ax == 0xbabe) && (r.h.dl != 0xff))
     {
      if ((r.h.dl & 0x40)==0) {WC++;break;}
/*
      if (r.h.dl & 0x80) printf("Drive %c is not cached.\n",'A'+Drive);
      else               printf("Drive %c is read cached.\n",'A'+Drive);
      if (r.h.dl & 0x40) printf("Drive %c is write through.\n",'A'+Drive);
      else               printf("Drive %c is write cached.\n",'A'+Drive);
*/
     }
  }

if (WC)
  {
   fprintf(stderr,"Notice:  It appears that at least one drive has write caching\n");
   fprintf(stderr,"enabled.  Although this will not prevent the program from running,\n");
   fprintf(stderr,"it is strongly discouraged because it greatly increases the number\n");
   fprintf(stderr,"of disk head movements.  For SmartDrv, you can disable all disk\n");
   fprintf(stderr,"write caching with the command: smartdrv /x\n\n");
  }
}

#else
void CheckDiskCache(void)
{
}
#endif


/*
** Return the amount of physical and virtual memory currently available.
*/

#if defined(__DJGPP__)
#include <go32.h>
#include <dpmi.h>

size_t
GetPhysicalMemSize(void)
{
 return (size_t)_go32_dpmi_remaining_physical_memory();
}

size_t
GetVirtualMemSize(void)
{
 return (size_t)_go32_dpmi_remaining_virtual_memory();
}

#else
/*
** No way to know how much memory is actually available, so assume we
** have an infinite amount.  Of course, since I don't know whether size_t
** is an unsigned integer or unsigned long integer, I can't use the
** proper macro.  But, since 32 bit computers will have long and int the
** same size, it actually works out okay.
*/

size_t
GetPhysicalMemSize(void)
{
/* return (size_t)ULONG_MAX;*/
 return (size_t)UINT_MAX;
}

size_t
GetVirtualMemSize(void)
{
/* return (size_t)ULONG_MAX;*/
 return (size_t)UINT_MAX;
}

#endif


