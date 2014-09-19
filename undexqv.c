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
 *  Uncompressor for .dexqv files
 *
 *  Author:  Gene Myers
 *  Date:    Jan 18, 2014
 *
 ********************************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

#include "DB.h"

static char *Usage = "[-vk] <path:dexqv> ...";

static void flip_short(void *w)
{ uint8 *v = (uint8 *) w;
  uint8  x;

  x    = v[0];
  v[0] = v[1];
  v[1] = x;
}

int main(int argc, char* argv[])
{ int VERBOSE;
  int KEEP;

  { int i, j, k;
    int flags[128];

    ARG_INIT("undexqv")

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

  //  For each .dexqv file to be decompressed

  { int   i;
    char *entry[5] = { NULL, NULL, NULL, NULL, NULL };
    int   emax     = -1;
    
    for (i = 1; i < argc; i++)
      { char     *pwd, *root;
        FILE     *input, *output;
        QVcoding *coding;

        //   Open it and the appropriately named .quiva file

        pwd   = PathTo(argv[i]);
        root  = Root(argv[i],".dexqv");
        input = Fopen(Catenate(pwd,"/",root,".dexqv"),"r");
        if (input == NULL)
          exit (1);
        output = Fopen(Catenate(pwd,"/",root,".quiva"),"w");
        if (output == NULL)
          exit (1);

        if (VERBOSE)
          { fprintf(stderr,"Processing '%s' ...\n",root);
            fflush(stderr);
          }

        // Read in compression scheme

        coding = Read_QVcoding(input);

        //  For each compressed entry do

        { int well;

          well = 0;
          while (1)
            { int    beg, end, qv, rlen;
              uint16 half;
              uint8  byte;
              int    e;

              //  Decode the compressed header and write it out

              if (fread(&byte,1,1,input) < 1) break;
              while (byte == 255)
                { well += 255;
                  if (fread(&byte,1,1,input) != 1)
                    SYSTEM_ERROR
                }
              well += byte;

              if (coding->flip)
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

              fprintf(output,"%s/%d/%d_%d RQ=0.%d\n",coding->prefix,well,beg,end,qv);

              //  Decode the QV entry and write it out

              rlen = end-beg;
              if (rlen > emax)
                { emax = ((int) (1.2*rlen)) + 1000;
                  entry[0] = (char *) Realloc(entry[0],5*emax,"Reallocating QV entry buffer");
                  if (entry[0] == NULL)
                    exit (1);
                  for (e = 1; e < 5; e++)
                    entry[e] = entry[e-1] + emax;
                }

              Uncompress_Next_QVentry(input,entry,coding,rlen);

              for (e = 0; e < 5; e++)
                fprintf(output,"%.*s\n",rlen,entry[e]);
            }
	}

        //  Clean up for the next file

	Free_QVcoding(coding);

        fclose(input);
        fclose(output);

        if (!KEEP)
          unlink(Catenate(pwd,"/",root,".dexqv"));
        free(root);
        free(pwd);

        if (VERBOSE)
          { fprintf(stderr,"Done\n");
            fflush(stderr);
          }
      }
  }

  free(QVentry());

  exit (0);
}
