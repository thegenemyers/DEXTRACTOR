/************************************************************************************\
*                                                                                    *
* Copyright (c) 2014, Dr. Eugene W. Myers (EWM). All rights reserved.                *
*                                                                                    *
* Redistribution and use in source and binary forms, with or without modification,   *
* are permitted provided that the following conditions are met:                      *
*                                                                                    *
*  · Redistributions of source code must retain the above copyright notice, this     *
*    list of conditions and the following disclaimer.                                *
*                                                                                    *
*  · Redistributions in binary form must reproduce the above copyright notice, this  *
*    list of conditions and the following disclaimer in the documentation and/or     *
*    other materials provided with the distribution.                                 *
*                                                                                    *
*  · The name of EWM may not be used to endorse or promote products derived from     *
*    this software without specific prior written permission.                        *
*                                                                                    *
* THIS SOFTWARE IS PROVIDED BY EWM ”AS IS” AND ANY EXPRESS OR IMPLIED WARRANTIES,    *
* INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND       *
* FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL EWM BE LIABLE   *
* FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES *
* (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS  *
* OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY      *
* THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING     *
* NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN  *
* IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                                      *
*                                                                                    *
* For any issues regarding this software and its use, contact EWM at:                *
*                                                                                    *
*   Eugene W. Myers Jr.                                                              *
*   Bautzner Str. 122e                                                               *
*   01099 Dresden                                                                    *
*   GERMANY                                                                          *
*   Email: gene.myers@gmail.com                                                      *
*                                                                                    *
\************************************************************************************/

/*******************************************************************************************
 *
 *  Compresses a .fasta file into a 2-bit per base .dexta file
 *
 *  Author:  Gene Myers
 *  Date  :  January 12, 2014
 *
 ********************************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <unistd.h>
#include <sys/stat.h>

#include "DB.h"

static char *Usage = "[-vk] <path:fasta> ...";

#define MAX_BUFFER 100000

//  Compress read into 2-bits per base (from [0-3] per byte representation

int main(int argc, char *argv[])
{ int     VERBOSE;
  int     KEEP;

  { int i, j, k;
    int flags[128];

    ARG_INIT("dexta")

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        { ARG_FLAGS("vk") }
      else
        argv[j++] = argv[i];
    argc = j;

    VERBOSE = flags['v'];
    KEEP    = flags['k'];

    if (argc == 1)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        exit (1);
      }
  }

  // For each fasta file do:

  { char   *read;
    int     rmax;
    int     i;

    rmax  = MAX_BUFFER + 30000;
    read  = (char *) Malloc(rmax+1,"Allocating read buffer");
    if (read == NULL)
      exit (1);

    for (i = 1; i < argc; i++)

      { char *pwd, *root;
        FILE *input, *output;
        int   eof;

        // Open fasta file

        pwd   = PathTo(argv[i]);
        root  = Root(argv[i],".fasta");
        input = Fopen(Catenate(pwd,"/",root,".fasta"),"r");
        if (input == NULL)
          exit (1);
        output = Fopen(Catenate(pwd,"/",root,".dexta"),"w");
        if (output == NULL)
          exit (1);

        if (VERBOSE)
          { fprintf(stderr,"Processing '%s' ...\n",root);
            fflush(stderr);
          }

        // Read the first header and output the endian key and short name

        { char  *slash;
          uint16 half;
          int    x;

          eof = (fgets(read,MAX_BUFFER,input) == NULL);
          if (read[strlen(read)-1] != '\n')
            { fprintf(stderr,"Line 1: Fasta line is too long (> %d chars)\n",MAX_BUFFER-2);
              exit (1);
            }
          if (!eof && read[0] != '>')
            { fprintf(stderr,"Line 1: First header in fasta file is missing\n");
              exit (1);
            }

          slash = index(read,'/');
          if (slash == NULL)
            { fprintf(stderr,"%s: Header line incorrectly formatted ?\n",Prog_Name);
              exit (1);
            }

          half = 0x33cc;
          fwrite(&half,sizeof(uint16),1,output);

          x = slash-read;
          fwrite(&x,sizeof(int),1,output);
          fwrite(read,1,slash-read,output);
        }

        //  For each read do

        { int  nline, rlen, lwell;

          nline = 1;
          rlen  = 0;
          lwell = 0;
          while (!eof)
            { int    well, beg, end, qv;
              char  *slash;
              uint16 half;
              uint8  byte;

              //  Next header is always at read+(rlen+1).  Interpret its fields

              slash = index(read+(rlen+1),'/');
              if (slash == NULL)
                { fprintf(stderr,"%s: Header line incorrectly formatted ?\n",Prog_Name);
    		    exit (1);
                }
              if (sscanf(slash+1,"%d/%d_%d RQ=0.%d\n",&well,&beg,&end,&qv) != 4)
                { fprintf(stderr,"%s: Header line incorrectly formatted ?\n",Prog_Name);
                  exit (1);
                }

              //  Read fasta sequence (@read) and stop at eof or after having read next header

              rlen = 0;
              while (1)
                { int x;

                  eof = (fgets(read+rlen,MAX_BUFFER,input) == NULL);
                  nline += 1;
                  x      = strlen(read+rlen)-1;
                  if (read[rlen+x] != '\n')
                    { fprintf(stderr,"Line %d: Fasta line is too long (> %d chars)\n",
                                     nline,MAX_BUFFER-2);
                      exit (1);
                    }
                  if (eof || read[rlen] == '>')
                    break;
                  rlen += x;
                  if (rlen + MAX_BUFFER > rmax)
                    { rmax = ((int) (1.2 * rmax)) + 1000 + MAX_BUFFER;
                      read = (char *) Realloc(read,rmax+1,"Reallocaing read buffer");
                      if (read == NULL)
                        exit (1);
                    }
                }
              read[rlen] = '\0';

              //  Compress the header fields and output (except for short name, only output once)

              while (well - lwell >= 255)
                { byte = 0xff;
                  fwrite(&byte,1,1,output);
                  lwell += 255;
                }
              byte = (uint8) (well-lwell);
              fwrite(&byte,1,1,output);
              lwell = well;

              half = (uint16) beg;
              fwrite(&half,sizeof(uint16),1,output);
              half = (uint16) end;
              fwrite(&half,sizeof(uint16),1,output);
              half = (uint16) qv;
              fwrite(&half,sizeof(uint16),1,output);

              //  Compress read and output

              Number_Read(read);
              Compress_Read(rlen,read);
              fwrite(read,1,COMPRESSED_LEN(rlen),output);
            }
        }

        if (!KEEP)
          unlink(Catenate(pwd,"/",root,".fasta"));
        free(root);
        free(pwd);

        if (VERBOSE)
          { fprintf(stderr,"Done\n");
            fflush(stderr);
          }
      }

    free(read);
  }

  exit (0);
}
