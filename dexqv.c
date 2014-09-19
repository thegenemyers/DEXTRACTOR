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
 *  Compressor for .quiv files, customized Huffman codes for each stream based on the
 *    histogram of values occuring in the given file.  The two low complexity streams
 *    (deletionQV and substitutionQV) use a Huffman coding of the run length of the prevelant
 *    character.
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

static char *Usage = "[-vkl] <path:quiva> ...";

int main(int argc, char* argv[])
{ int        VERBOSE;
  int        KEEP;
  int        LOSSY;

  { int   i, j, k;
    int   flags[128];

    ARG_INIT("dexqv")

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        { ARG_FLAGS("vkl") }
      else
        argv[j++] = argv[i];
    argc = j;

    VERBOSE = flags['v'];
    KEEP    = flags['k'];
    LOSSY   = flags['l'];

    if (argc == 1)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        exit (1);
      }
  }

  // For each .quiva file to be compressed:

  { int       i;

    for (i = 1; i < argc; i++)
      { char     *pwd, *root;
        FILE     *input, *output;
        QVcoding *coding;
        
        pwd   = PathTo(argv[i]);
        root  = Root(argv[i],".quiva");
        input = Fopen(Catenate(pwd,"/",root,".quiva"),"r");
        if (input == NULL)
          exit (1);
        output = Fopen(Catenate(pwd,"/",root,".dexqv"),"w");
        if (output == NULL)
          exit (1);

        if (VERBOSE)
          { fprintf(stderr,"Processing '%s' ...\n",root);
            fflush(stderr);
          }

        //  Scan the file collecting statistics for Huffman schemes
        
        QVcoding_Scan(input);

        //  Create and output the encoding schemes

        coding = Create_QVcoding(LOSSY);

        { char *slash, *read;

          rewind (input);
          Read_Lines(input,1);
          read = QVentry();

          slash = index(read+1,'/');
          if (slash != NULL)
            { coding->prefix = (char *) malloc((slash-read)+1);
              if (coding->prefix == NULL)
                 { fprintf(stderr,"%s: Out of memory (Allocating header prefix)\n",Prog_Name);
                   exit (1);
                 }
              *slash = '\0';
              if (strcpy(coding->prefix,read) == NULL)
                SYSTEM_ERROR
              *slash = '/';
            }
        }

        Write_QVcoding(output,coding);

        //  For each entry do

        { int lwell;

          rewind (input);

          lwell = 0;
          while (Read_Lines(input,1) > 0)
            { int    well, beg, end, qv;
              char  *slash;
	      uint16 half;
              uint8  byte;

              //  Interpret the header, encode and write out the fields

              slash = index(QVentry(),'/');
              sscanf(slash+1,"%d/%d_%d RQ=0.%d\n",&well,&beg,&end,&qv);

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

              Compress_Next_QVentry(input,output,coding,LOSSY);
            }
        }

        //  Clean up for the next file

        Free_QVcoding(coding);

        fclose(input);
        fclose(output);

        if (!KEEP)
          unlink(Catenate(pwd,"/",root,".quiva"));
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
