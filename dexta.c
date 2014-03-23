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

#define MAX_BUFFER 100000

static char *Usage = "[-Sk] <path:fasta> ...";

#define DEXTA

#include "shared.c"

//  Compress read into 2-bits per base (from [0-3] per byte representation

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

  // For each fasta file do:

  { char   *read;
    int     rmax;
    int     i;

    rmax  = MAX_BUFFER + 30000;
    read  = (char *) Guarded_Alloc(NULL,rmax+1);
    for (i = 1; i < argc; i++)

      { FILE   *input, *output;
        char   *full;
        int     eof;

        // Open fasta file

        { char *path;
          int   epos;

          path = argv[i];
          full = (char *) Guarded_Alloc(NULL,strlen(path)+20);
          epos = strlen(path);
          if (epos >= 6 && strcasecmp(path+(epos-6),".fasta") == 0)
            strcpy(full,path);
          else
            { epos += 6;
              sprintf(full,"%s.fasta",path);
            }

          input = Guarded_Fopen(full,"r");

          if (VERBOSE)
            { fprintf(stderr,"Processing '%s' ...",full);
              fflush(stderr);
            }

          strcpy(full+(epos-6),".dexta");
          output = Guarded_Fopen(full,"w");

          strcpy(full+(epos-6),".fasta");
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
            { fprintf(stderr,"%s: Header line incorrectly formatted ?\n",Program_Name);
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
                { fprintf(stderr,"%s: Header line incorrectly formatted ?\n",Program_Name);
    		    exit (1);
                }
              if (sscanf(slash+1,"%d/%d_%d RQ=0.%d\n",&well,&beg,&end,&qv) != 4)
                { fprintf(stderr,"%s: Header line incorrectly formatted ?\n",Program_Name);
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
                    { rmax = 1.2 * rmax + 1000 + MAX_BUFFER;
                      read = (char *) Guarded_Alloc(read,rmax+1);
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

              Compress_Read(rlen,read);
              fwrite(read,1,COMPRESSED_LEN(rlen),output);
            }
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
