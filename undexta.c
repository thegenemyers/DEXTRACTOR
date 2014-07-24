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
        root  = Root(argv[i],".fasta");
        input = Fopen(Catenate(pwd,"/",root,".dexta"),"r");
        if (input == NULL)
          exit (1);
        output = Fopen(Catenate(pwd,"/",root,".fasta"),"w");
        if (output == NULL)
          exit (1);

        if (VERBOSE)
          { fprintf(stderr,"Processing '%s' ...",root);
            fflush(stderr);
          }

        { char *name;
          int   well, flip;

          // Read endian key and short name common to all headers

          { uint16 half;

            fread(&half,sizeof(uint16),1,input);
            flip = (half != 0x33cc);

            fread(&well,sizeof(int),1,input);
            if (flip) flip_long(&well);
            name = (char *) Malloc(well+1,"Allocating header prefix");
            fread(name,1,well,input);
            name[well] = '\0';
          }

          // For each encoded entry do

          well = 0;
          while (1)
            { int    rlen, beg, end, qv;
              uint16 half;
              uint8  byte;

              //  Read and decompress header and output

              if (fread(&byte,1,1,input) < 1) break;
              while (byte == 255)
                { well += 255;
                  fread(&byte,1,1,input);
                }
              well += byte;

              if (flip)
                { fread(&half,sizeof(uint16),1,input);
                  flip_short(&half);
                  beg = half;
                  fread(&half,sizeof(uint16),1,input);
                  flip_short(&half);
                  end = half;
                  fread(&half,sizeof(uint16),1,input);
                  flip_short(&half);
                  qv = half;
                }
              else
                { fread(&half,sizeof(uint16),1,input);
                  beg = half;
                  fread(&half,sizeof(uint16),1,input);
                  end = half;
                  fread(&half,sizeof(uint16),1,input);
                  qv = half;
                }

              fprintf(output,"%s/%d/%d_%d RQ=0.%d\n",name,well,beg,end,qv);

              //  Read compressed sequence (into buffer big enough for uncompressed sequence)
              //  Uncompress and output 80 symbols to a line

              rlen = end-beg;
              if (rlen > rmax)
                { rmax = 1.2 * rmax + 1000 + MAX_BUFFER;
                  read = (char *) Realloc(read,rmax+1,"Allocating read buffer");
                }
              fread(read,1,COMPRESSED_LEN(rlen),input);
              Uncompress_Read(rlen,read);
              Lower_Read(read);

              for (i = 0; i < rlen; i += 80)
                if (i+80 > rlen)
                  fprintf(output,"%.*s\n", rlen-i, read+i);
                else
                  fprintf(output,"%.80s\n", read+i);
            }

          free(name);
        }

        if (!KEEP)
          unlink(Catenate(pwd,"/",root,".dexta"));
        free(root);
        free(pwd);

        if (VERBOSE)
          { fprintf(stderr," Done\n");
            fflush(stderr);
          }
      }

    free(read);
  }

  exit (0);
}
