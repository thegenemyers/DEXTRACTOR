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
 *  Uncompresses a .dexta file (2-bit per base compression) back to a .fasta file
 *
 *  Author:  Gene Myers
 *  Date  :  January 12, 2014
 ********************************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <unistd.h>
#include <sys/stat.h>

#include "DB.h"

static char *Usage = "[-vk] <path:dexta> ...";

#define MAX_BUFFER 100000

//  Uncompress read from 2-bits per base into [0-3] per byte representation

static void flip_long(void *w)
{ uint8 *v = (uint8 *) w;
  uint8  x;

  x    = v[0];
  v[0] = v[3];
  v[3] = x;
  x    = v[1];
  v[1] = v[2];
  v[2] = x;
}

static void flip_short(void *w)
{ uint8 *v = (uint8 *) w;
  uint8  x;

  x    = v[0];
  v[0] = v[1];
  v[1] = x;
}

int main(int argc, char *argv[])
{ int     VERBOSE;
  int     KEEP;

  { int i, j, k;
    int flags[128];

    ARG_INIT("undexta")

    VERBOSE = 1;
    KEEP    = 0;

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

  // For each .dexta file do

  { char   *read;
    int     rmax;
    int     i;

    rmax  = MAX_BUFFER + 30000;
    read  = (char *) Malloc(rmax+1,"Allocating read buffer");
    for (i = 1; i < argc; i++)
      { char *pwd, *root;
        FILE *input, *output;

        // Open dexta file

        pwd   = PathTo(argv[i]);
        root  = Root(argv[i],".dexta");
        input = Fopen(Catenate(pwd,"/",root,".dexta"),"r");
        if (input == NULL)
          exit (1);
        output = Fopen(Catenate(pwd,"/",root,".fasta"),"w");
        if (output == NULL)
          exit (1);

        if (VERBOSE)
          { fprintf(stderr,"Processing '%s' ...\n",root);
            fflush(stderr);
          }

        { char *name;
          int   well, flip;

          // Read endian key and short name common to all headers

          { uint16 half;

            if (fread(&half,sizeof(uint16),1,input) != 1)
              SYSTEM_ERROR
            flip = (half != 0x33cc);

            if (fread(&well,sizeof(int),1,input) != 1)
              SYSTEM_ERROR
            if (flip) flip_long(&well);
            name = (char *) Malloc(well+1,"Allocating header prefix");
            if (well > 0)
              { if (fread(name,well,1,input) != 1)
                  SYSTEM_ERROR
              }
            name[well] = '\0';
          }

          // For each encoded entry do

          well = 0;
          while (1)
            { int    rlen, beg, end, qv;
              int    clen;
              uint16 half;
              uint8  byte;

              //  Read and decompress header and output

              if (fread(&byte,1,1,input) < 1) break;
              while (byte == 255)
                { well += 255;
                  if (fread(&byte,1,1,input) != 1)
                    SYSTEM_ERROR
                }
              well += byte;

              if (flip)
                { if (fread(&half,sizeof(uint16),1,input) != 1)
                    SYSTEM_ERROR
                  flip_short(&half);
                  beg = half;
                  if (fread(&half,sizeof(uint16),1,input) != 1)
                    SYSTEM_ERROR
                  flip_short(&half);
                  end = half;
                  if (fread(&half,sizeof(uint16),1,input) != 1)
                    SYSTEM_ERROR
                  flip_short(&half);
                  qv = half;
                }
              else
                { if (fread(&half,sizeof(uint16),1,input) != 1)
                    SYSTEM_ERROR
                  beg = half;
                  if (fread(&half,sizeof(uint16),1,input) != 1)
                    SYSTEM_ERROR
                  end = half;
                  if (fread(&half,sizeof(uint16),1,input) != 1)
                    SYSTEM_ERROR
                  qv = half;
                }

              fprintf(output,"%s/%d/%d_%d RQ=0.%d\n",name,well,beg,end,qv);

              //  Read compressed sequence (into buffer big enough for uncompressed sequence)
              //  Uncompress and output 80 symbols to a line

              rlen = end-beg;
              if (rlen > rmax)
                { rmax = ((int) (1.2 * rmax)) + 1000 + MAX_BUFFER;
                  read = (char *) Realloc(read,rmax+1,"Allocating read buffer");
                }
              clen = COMPRESSED_LEN(rlen);
              if (clen > 0)
                { if (fread(read,clen,1,input) != 1)
                    SYSTEM_ERROR
                }
              Uncompress_Read(rlen,read);
              Lower_Read(read);

              { int j;

                for (j = 0; j < rlen; j += 80)
                  if (j+80 > rlen)
                    fprintf(output,"%.*s\n", rlen-j, read+j);
                  else
                    fprintf(output,"%.80s\n", read+j);
              }
            }

          free(name);
        }

        if (!KEEP)
          unlink(Catenate(pwd,"/",root,".dexta"));
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
