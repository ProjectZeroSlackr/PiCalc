/*
** program to extract 1,000,000 digit markers from Guttenberg pi file.
*/
#include <stdio.h>
#include <errno.h>

/* whether line feeds on disk are 1 or 2 chars */
#define LF 2

/* Print out 100 every ??? digits */
#define INC 1000000

char Str[800];
FILE *fin;
int count=0;


int main(void)
{int pos; /* current pos */
 int eop; /* end of file */
 int iop; /* initial pos */

fin=fopen("f:\\pi512m.txt","r");
if (fin==NULL) return 0;

fseek(fin,0,SEEK_END);
eop=ftell(fin);
fseek(fin,0,SEEK_SET);
iop=pos=2+2*LF;

while (pos < eop)
  {
   fseek(fin,pos,SEEK_SET);
   fgets(Str,200,fin);
    printf(       "%3d million %s",count/INC,Str);
   fprintf(stderr,"%3d million %s",count/INC,Str);
   fgets(Str,200,fin);
    printf(       "            %s",Str);
   fprintf(stderr,"            %s",Str);
   count+=INC;
   {int line;
    line=count/50;
    pos=iop+(line*(54+LF));
    pos+=(line/20)*(3*LF);
   }
   if (ferror(fin))
     {
      fprintf(stderr,"Error reading file. errno=%d\n",errno);
      break;
     }
  }

fclose(fin);
return 0;
}


