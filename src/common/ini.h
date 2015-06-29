/*
** This .ini management code is placed into the Public Domain
** by Carey Bloodworth on May 26, 1998.
*/
#ifndef INI_H_
#define INI_H_

enum CfgVarTypes {Cfg_String=0,  /* A character string           */
                  Cfg_Integer,   /* 16/32 bit signed int         */
                  Cfg_UInteger,  /* 16/32 bit unsigned int       */
                  Cfg_SInteger,  /* Short int                    */
                  Cfg_USInteger, /* Unsigned short int           */
                  Cfg_LInteger,  /* signed Long int              */
                  Cfg_ULInteger, /* Unsigned long int            */
                  Cfg_Boolean,   /* Boolean is 'int' N/Y F/T 0/1 */
                  Cfg_Bytes      /* A 'char' but not a string    */
                 };

struct CfgNameStruct
       {char *Section;
        char *Name;
        void *DataPtr;
        enum CfgVarTypes VarType;
        unsigned int Max;
       };


unsigned int ReadCfgItem(char *FileName,
                char *SectionName,
                char *VarName,
                void *DataPtr,
                enum CfgVarTypes DataType,
                unsigned int MaxItems);

unsigned int UpdateCfgItem(char *FileName,
                  char *SectionName,
                  char *VarWanted,
                  void *DataPtr,
                  enum CfgVarTypes DataType,
                  unsigned int HowMany);

int
LoadConfigStruct(char *FileName, struct CfgNameStruct *mv);

#endif

