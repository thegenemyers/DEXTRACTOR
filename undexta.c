/*******************************************************************************************
 *
 *  Uncompresses a .dexta file (2-bit per base compression) back to a .fasta file
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

#define MAX_BUFFER 100000

static char *Usage = "[-Sk] <path:dexta> ...";

#define UNDEXTA

#include "shared.c"

//  Uncompress read form 2-bits per base into [0-3] per byte representation

int main(int argc, char *argv[])
{ int     VERBOSE;
  int     KEEP;

  Program_Name = argv[0];

  { int   i, j, k;

    VERBOSE = 1;
    KEEP    = 0;

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        for (k = 1; argv[i][k] != '\0'; k++)
          if (argv[i][k] == 'S')
            VERBOSE = 0;
          else if (argv[i][k] == 'k')
            KEEP = 1;
          else
            { fprintf(stderr,"%s: -%c is an illegal option\n",Program_Name,argv[i][k]);
              exit (1);
            }
      else
        argv[j++] = argv[i];
    argc = j;

    if (argc == 1)
      { fprintf(stderr,"Usage: %s %s\n",Program_Name,Usage);
        exit (1);
      }
  }

  // For each .dexta file do

  { char   *read;
    int     rmax;
    int     i;

    rmax  = MAX_BUFFER + 30000;
    read  = (char *) Guarded_Alloc(NULL,rmax+1);
    for (i = 1; i < argc; i++)
      { FILE   *input, *output;
        char   *full;

        // Open dexta file

        { char *path;
          int   epos;

          path = argv[i];
          full = (char *) Guarded_Alloc(NULL,strlen(path)+20);
          epos = strlen(path);
          if (epos >= 6 && strcasecmp(path+(epos-6),".dexta") == 0)
            strcpy(full,path);
          else
            { epos += 6;
              sprintf(full,"%s.dexta",path);
            }

          input = Guarded_Fopen(full,"r");

          if (VERBOSE)
            { fprintf(stderr,"Processing '%s' ...",full);
              fflush(stderr);
            }

          strcpy(full+(epos-6),".fasta");
          output = Guarded_Fopen(full,"w");

          strcpy(full+(epos-6),".dexta");
        }

        { char *name;
          int   well, flip;

          // Read endian key and short name common to all headers

          { uint16 half;

            fread(&half,sizeof(uint16),1,input);
            flip = (half != 0x33cc);

            fread(&well,sizeof(int),1,input);
            if (flip) flip_long(&well);
            name = (char *) Guarded_Alloc(NULL,well+1);
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
              //  Uncompress and output 70 symbols to a line

              rlen = end-beg;
              if (rlen > rmax)
                { rmax = 1.2 * rmax + 1000 + MAX_BUFFER;
                  read = (char *) Guarded_Alloc(read,rmax+1);
                }
              fread(read,1,COMPRESSED_LEN(rlen),input);
              Uncompress_Read(rlen,read);

              for (i = 0; i < rlen; i += 70)
                if (i+70 > rlen)
                  fprintf(output,"%.*s\n", rlen-i, read+i);
                else
                  fprintf(output,"%.70s\n", read+i);
            }

          free(name);
        }

        if (!KEEP)
          unlink(full);
        free(full);

        if (VERBOSE)
          { fprintf(stderr," Done\n");
            fflush(stderr);
          }
      }

    free(read);
  }

  exit (0);
}
