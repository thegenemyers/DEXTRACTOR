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

#undef  DEBUG

#define MIN_BUFFER 1000

static char *Usage = "[-Sk] <path:dexqv> ...";

#define UNDEXQV

#include "shared.c"

static int LittleEndian;  //  Little-endian machine ?
static int Flip;          //  Flip endian of all input shorts and ints


/*******************************************************************************************
 *
 *  Huffman Decoding Routines
 *
 ********************************************************************************************/

typedef struct
  { int    type;             //  0 => normal, 2 => truncated
    uint32 codebits[256];    //  If type = 2, then code 255 is the special code for
    int    codelens[256];    //    non-Huffman exceptions
    int    lookup[0x10000];
  } HScheme;

#ifdef DEBUG

static void Print_Table(HScheme *scheme)
{ uint32 specval, mask, code, *bits;
  int    speclen, clen, *lens;
  int    i, k;

  bits = scheme->codebits;
  lens = scheme->codelens;
  if (scheme->type == 2)
    { specval = bits[255];
      speclen = lens[255];
    }
  else
    specval = speclen = 0x7fffffff;

  printf("\nCode Table:\n");
  for (i = 0; i < 256; i++)
    if (lens[i] > 0)
      { clen = lens[i];
        mask = (1 << clen);
        code = bits[i];
        printf(" %3d: %2d ",i,clen);
        for (k = 0; k < clen; k++)
          { mask >>= 1;
            if (code & mask)
              printf("1");
            else
              printf("0");
          }
        if (code == specval && clen == speclen)
          printf(" ***");
        printf("\n");
      }
}

#endif
 
  //  Allocate and read a code table from in, and return a pointer to it.

static HScheme *Read_Scheme(FILE *in)
{ HScheme *scheme;
  int     *look, *lens;
  uint32  *bits, base;
  int      i, j, powr;
  uint8    x;

  scheme = (HScheme *) Guarded_Alloc(NULL,sizeof(HScheme));

  lens = scheme->codelens;
  bits = scheme->codebits;
  look = scheme->lookup;

  fread(&x,1,1,in);
  scheme->type = x;
  for (i = 0; i < 256; i++)
    { fread(&x,1,1,in);
      lens[i] = x;
      if (x > 0)
        fread(bits+i,sizeof(uint32),1,in);
      else
        bits[i] = 0;
    }

  if (Flip)
    { for (i = 0; i < 256; i++)
        flip_long(bits+i);
    }

  for (i = 0; i < 256; i++)
    { if (lens[i] > 0)
        { base = (bits[i] << (16-lens[i]));
          powr = (1 << (16-lens[i]));
          for (j = 0; j < powr; j++)
            look[base+j] = i;
        }
    }

  return (scheme);
}

  //  Read and decode from in, the next rlen symbols into read according to scheme

static void Decode(HScheme *scheme, FILE *in, char *read, int rlen)
{ int    *look, *lens;
  int     signal, ilen;
  uint64  icode;
  uint32 *ipart;
  uint16 *xpart;
  uint8  *cpart;
  int     j, n, c;

  if (LittleEndian)
    { ipart = ((uint32 *) (&icode));
      xpart = ((uint16 *) (&icode)) + 2;
      cpart = ((uint8  *) (&icode)) + 5;
    }
  else
    { ipart = ((uint32 *) (&icode)) + 1;
      xpart = ((uint16 *) (&icode)) + 1;
      cpart = ((uint8  *) (&icode)) + 2;
    }

  if (scheme->type == 2)
    signal  = 255;
  else
    signal  = 256;
  lens = scheme->codelens;
  look = scheme->lookup;

#define GET					\
  if (n > ilen)					\
    { icode <<= ilen;				\
      fread(ipart,sizeof(uint32),1,in);		\
      ilen    = n-ilen;				\
      icode <<= ilen;				\
      ilen    = 32-ilen;			\
    }						\
  else						\
    { icode <<= n;				\
      ilen   -= n;				\
    }

#define GETFLIP					\
  if (n > ilen)					\
    { icode <<= ilen;				\
      fread(ipart,sizeof(uint32),1,in);		\
      flip_long(ipart);				\
      ilen    = n-ilen;				\
      icode <<= ilen;				\
      ilen    = 32-ilen;			\
    }						\
  else						\
    { icode <<= n;				\
      ilen   -= n;				\
    }

  n    = 16;
  ilen = 0;
  if (Flip)
    for (j = 0; j < rlen; j++)
      { GETFLIP
        c = look[*xpart];
        n = lens[c];
        if (c == signal)
          { GETFLIP
            c = *cpart;
            n = 8;
          }
        read[j] = c;
      }
  else
    for (j = 0; j < rlen; j++)
      { GET
        c = look[*xpart];
        n = lens[c];
        if (c == signal)
          { GET
            c = *cpart;
            n = 8;
          }
        read[j] = c;
      }
}

  //  Read and decode from in, the next rlen symbols into read according to non-rchar scheme
  //    neme, and the rchar runlength shceme reme

static void Decode_Run(HScheme *neme, HScheme *reme, FILE *in, char *read,
                       int rlen, int rchar)
{ int    *nlook, *nlens;
  int    *rlook, *rlens;
  int     nsignal, rsignal, ilen;
  uint64  icode;
  uint32 *ipart;
  uint16 *xpart;
  uint8  *cpart;
  int     j, n, c, k;

  if (LittleEndian)
    { ipart = ((uint32 *) (&icode));
      xpart = ((uint16 *) (&icode)) + 2;
      cpart = ((uint8  *) (&icode)) + 5;
    }
  else
    { ipart = ((uint32 *) (&icode)) + 1;
      xpart = ((uint16 *) (&icode)) + 1;
      cpart = ((uint8  *) (&icode)) + 2;
    }

  if (neme->type == 2)
    nsignal = 255;
  else
    nsignal = 256;
  nlens = neme->codelens;
  nlook = neme->lookup;

  if (reme->type == 2)
    rsignal = 255;
  else
    rsignal = 256;
  rlens = reme->codelens;
  rlook = reme->lookup;

  n    = 16;
  ilen = 0;
  if (Flip)
    for (j = 0; j < rlen; j++)
      { GETFLIP
        c = rlook[*xpart];
        n = rlens[c];
        if (c == rsignal)
          { GETFLIP
            c = *xpart;
            n = 16;
          }
        for (k = 0; k < c; k++)
          read[j++] = rchar;

        if (j < rlen)
          { GETFLIP
            c = nlook[*xpart];
            n = nlens[c];
            if (c == nsignal)
              { GETFLIP
                c = *cpart;
                n = 8;
              }
            read[j] = c;
          }
      }
  else
    for (j = 0; j < rlen; j++)
      { GET
        c = rlook[*xpart];
        n = rlens[c];
        if (c == rsignal)
          { GET
            c = *xpart;
            n = 16;
          }
        for (k = 0; k < c; k++)
          read[j++] = rchar;

        if (j < rlen)
          { GET
            c = nlook[*xpart];
            n = nlens[c];
            if (c == nsignal)
              { GET
                c = *cpart;
                n = 8;
              }
            read[j] = c;
          }
      }
}


/*******************************************************************************************
 *
 *  Tag decompression routines
 *
 ********************************************************************************************/

  //  Count the # of non-rchar symbols in qvs[0..rlen-1]

static int Packed_Length(char *qvs, int rlen, int rchar)
{ int k, clen;

  clen = 0;
  for (k = 0; k < rlen; k++)
    if (qvs[k] != rchar)
      clen += 1;
  return (clen);
}

  //  Unpack tags by moving its i'th char to position k where qvs[k] is the i'th non-rchar
  //    symbol in qvs.  All other chars are set to rchar.  rlen is the length of qvs and
  //    the unpacked result, clen is the initial length of tags.

static void Unpack_Tag(char *tags, int clen, char *qvs, int rlen, int rchar)
{ int j, k;

  j = clen-1;
  for (k = rlen-1; k >= 0; k--)
    { if (qvs[k] == rchar)
        tags[k] = 'N';
      else
        tags[k] = tags[j--];
    }
}

/*******************************************************************************************
 *
 *  Main routine
 *
 ********************************************************************************************/

int main(int argc, char* argv[])
{ int       VERBOSE;
  int       KEEP;

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

  { uint32 x = 3;
    uint8 *b = (uint8 *) (&x);

    LittleEndian = (b[0] == 3);
  }

  //  For each .dexqv file to be decompressed

  { char   *read;
    int     rmax;
    int     i;

    rmax  = 0;
    read  = NULL;

    for (i = 1; i < argc; i++)
      { FILE *input, *output;
        char *full;
  
        //   Open it and the appropriately named .quiva file
  
        { char *path;
          int   epos;
  
          path = argv[i];
          full = (char *) Guarded_Alloc(NULL,strlen(path)+20);
          epos = strlen(path);
          if (epos >= 6 && strcasecmp(path+(epos-6),".dexqv") == 0)
            strcpy(full,path);
          else
            { epos += 6;
              sprintf(full,"%s.dexqv",path);
            }
  
          input = Guarded_Fopen(full,"r");
  
          if (VERBOSE)
            { fprintf(stderr,"Processing '%s' ... ",full);
              fflush(stderr);
            }
  
          strcpy(full+(epos-6),".quiva");
          output = Guarded_Fopen(full,"w");
  
          strcpy(full+(epos-6),".dexqv");
        }
  
        //  Read the short name common to all headers

        { HScheme  *delScheme, *insScheme, *mrgScheme, *subScheme;
          HScheme  *dRunScheme, *sRunScheme;
          int       delChar, subChar;
          char     *name;
          int       well;

          // Read endian key, run chars, and short name common to all headers

          { uint16 half;

            fread(&half,sizeof(uint16),1,input);
            Flip = (half != 0x33cc);

            fread(&half,sizeof(uint16),1,input);
            if (Flip) flip_short(&half);
            delChar = half;
            if (delChar >= 256)
              delChar = -1;

            fread(&half,sizeof(uint16),1,input);
            if (Flip) flip_short(&half);
            subChar = half;
            if (subChar >= 256)
              subChar = -1;

            fread(&well,sizeof(int),1,input);
            if (Flip) flip_long(&well);
            name = (char *) Guarded_Alloc(NULL,well+1);
            fread(name,1,well,input);
            name[well] = '\0';
          }
  
          //  Read the Huffman schemes used to compress the data
    
          delScheme  = Read_Scheme(input);
          if (delChar >= 0)
            dRunScheme = Read_Scheme(input);
          insScheme  = Read_Scheme(input);
          mrgScheme  = Read_Scheme(input);
          subScheme  = Read_Scheme(input);
          if (subChar >= 0)
            sRunScheme = Read_Scheme(input);
  
          //  For each compressed entry do
  
          well = 0;
          while (1)
            { int    rlen, clen, beg, end, qv;
              uint16 half;
              uint8  byte;
    
              //  Decode the compressed header and write it out
    
              if (fread(&byte,1,1,input) < 1) break;
              while (byte == 255)
                { well += 255;
                  fread(&byte,1,1,input);
                }
              well += byte;
    
              if (Flip)
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
    
              //  Make sure read buffer is big enough for the uncompressed streams
    
              rlen = end-beg;
              if (rlen > rmax)
                { rmax = 1.2 * rlen + 1000;
                  read = (char *) Guarded_Alloc(read,2*rmax+1);
                }
    
              //  Decode each stream and write to output
    
              if (delChar < 0)
                { Decode(delScheme, input, read, rlen);
                  clen = rlen;
                }
              else
                { Decode_Run(delScheme, dRunScheme, input, read, rlen, delChar);
                  clen = Packed_Length(read,rlen,delChar);
                }
              fprintf(output,"%.*s\n",rlen,read);
    
              fread(read+rmax,1,COMPRESSED_LEN(clen),input);
              Uncompress_Read(clen,read+rmax);
              if (delChar >= 0)
                Unpack_Tag(read+rmax,clen,read,rlen,delChar);
              fprintf(output,"%.*s\n",rlen,read+rmax);
    
              Decode(insScheme, input, read, rlen);
              fprintf(output,"%.*s\n",rlen,read);
    
              Decode(mrgScheme, input, read, rlen);
              fprintf(output,"%.*s\n",rlen,read);
    
              if (subChar < 0)
                Decode(subScheme, input, read, rlen);
              else
                Decode_Run(subScheme, sRunScheme, input, read, rlen, subChar);
              fprintf(output,"%.*s\n",rlen,read);
            }
  
          free(delScheme);
          free(insScheme);
          free(mrgScheme);
          free(subScheme);
          if (delChar >= 0)
            free(dRunScheme);
          if (subChar >= 0)
            free(sRunScheme);

          free(name);
        }
  
        fclose(input);
        fclose(output);
  
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
