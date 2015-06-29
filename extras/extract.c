/*
** Simple program to extract 1m digit markers from Guttenberg pi file.
*/
#include <stdio.h>
#include <errno.h>

char Str[800];
FILE *fin;
int count=0;


int main(void)
{
fin=fopen("pi128m.txt","r");
if (fin==NULL) return 0;

fgets(Str,80,fin);
fgets(Str,80,fin);

while (!feof(fin))
  {
   fgets(Str,80,fin);
   if (Str[0]=='\n') continue;
   if ((count % 1000000)==0)  printf("%3d million %s",count/1000000,Str);
   if ((count % 1000000)==50) printf("            %s",Str);
   count+=50;
   if (ferror(fin))
     {
      fprintf(stderr,"Error reading file. errno=%d\n",errno);
      break;
     }
  }

fclose(fin);
return 0;
}

