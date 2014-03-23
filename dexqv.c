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

#undef DEBUG

#define MIN_BUFFER 1000

#define HUFF_CUTOFF  16   //  This cannot be larger than 16 !

static char *Usage = "[-Skl] <path:quiva> ...";

#define DEXQV

#include "shared.c"

/*******************************************************************************************
 *
 *  Huffman Encoding Routines
 *
 ********************************************************************************************/

typedef struct
  { int    type;             //  0 => normal, 1 => normal but has long codes, 2 => truncated
    uint32 codebits[256];    //  If type = 2, then code 255 is the special code for
    int    codelens[256];    //    non-Huffman exceptions
  } HScheme;

typedef struct _HTree
  { struct _HTree *lft, *rgt; 
    uint64         count;
  } HTree;

  //  Establish heap property from node s down (1 is root, siblings of n are 2n and 2n+1)
  //    assuming s is the only perturbation in the tree.

static void Reheap(int s, HTree **heap, int hsize)
{ int      c, l, r;
  HTree   *hs, *hr, *hl;

  c  = s;
  hs = heap[s];
  while ((l = 2*c) <= hsize)
    { r  = l+1;
      hl = heap[l];
      hr = heap[r];
      if (r > hsize || hr->count > hl->count)
        { if (hs->count > hl->count)
            { heap[c] = hl;
              c = l;
            }
          else
            break;
        }
      else
        { if (hs->count > hr->count)
            { heap[c] = hr;
              c = r;
            }
          else
            break;
        }
    }
  if (c != s)
    heap[c] = hs;
}

  //  Given Huffman tree build a table of codes from it, the low-order codelens[s] bits
  //    of codebits[s] contain the code for symbol s.

static void Build_Table(HTree *node, int code, int len, uint32 *codebits, int *codelens)
{ if (node->rgt == NULL)
    { uint64 symbol = (uint64) (node->lft);
      codebits[symbol] = code;
      codelens[symbol] = len;
    }
  else
    { code <<= 1;
      len   += 1;
      Build_Table(node->lft,code,len,codebits,codelens);
      Build_Table(node->rgt,code+1,len,codebits,codelens);
    }
}

#ifdef DEBUG

  //  For debug, show the coding table

static void Print_Table(HScheme *scheme, uint64 *hist, int infosize)
{ uint64 total_bits;
  uint32 specval, mask, code, *bits;
  int    speclen, clen, *lens;
  int    i, k;

  total_bits = 0;
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
          { printf(" ***");
            total_bits += (clen+infosize)*hist[i];
          }
        else
          total_bits += clen*hist[i];
        printf("\n");
      }
  printf("\nTotal Bytes = %lld\n",(total_bits-1)/8+1);
}

  //  For debug, show the histogram

static void Print_Histogram(uint64 *hist)
{ int    i, low, hgh;
  uint64 count;

  for (hgh = 255; hgh >= 0; hgh--)
    if (hist[hgh] != 0)
      break;
  for (low = 0; low < 256; low++)
    if (hist[low] != 0)
      break;
  count = 0;
  for (i = low; i <= hgh; i++)
    count += hist[i];

  for (i = hgh; i >= low; i--)
    printf("    %3d: %8llu %5.1f%%\n",i,hist[i],(hist[i]*100.)/count);
}

#endif

  // For the non-zero symbols in hist, compute a huffman tree over them, and then
  //   build a table of the codes.  If inscheme is not NULL, then place all symbols
  //   with code 255 or with more than HUFF_CUTOFF bits in the encoding by inscheme
  //   as a single united entity, whose code signals that the value of these symbols
  //   occur explicitly in 8 (values) or 16 (run lengths) bits following the code.
  //   All the symbols in this class will have the same entry in the code table and
  //   255 is always in this class.

static HScheme *Huffman(uint64 *hist, HScheme *inscheme)
{ HScheme *scheme;
  HTree   *heap[257];
  HTree    node[512];
  int      hsize;
  HTree   *lft, *rgt;
  int      value, range;
  int     i;

  scheme = (HScheme *) Guarded_Alloc(NULL,sizeof(HScheme));

  hsize = 0;                        //  Load heap
  value = 0;
  if (inscheme != NULL)
    { node[0].count = 0;
      node[0].lft   = (HTree *) (uint64) 255;
      node[0].rgt   = NULL;
      heap[++hsize] = node+(value++);
    }
  for (i = 0; i < 256; i++)
    if (hist[i] > 0)
      { if (inscheme != NULL && (inscheme->codelens[i] > HUFF_CUTOFF || i == 255))
          node[0].count += hist[i];
        else
          { node[value].count = hist[i];
            node[value].lft   = (HTree *) (uint64) i;
            node[value].rgt   = NULL;
            heap[++hsize] = node+(value++);
          }
      }

  for (i = hsize/2; i >= 1; i--)    //  Establish heap property
    Reheap(i,heap,hsize);

  range = value;                    //   Merge pairs with smallest count until have a tree
  for (i = 1; i < value; i++)
    { lft = heap[1];
      heap[1] = heap[hsize--];
      Reheap(1,heap,hsize);
      rgt = heap[1];
      node[range].lft = lft;
      node[range].rgt = rgt;
      node[range].count = lft->count + rgt->count;
      heap[1] = node+(range++);
      Reheap(1,heap,hsize);
    }

  for (i = 0; i < 256; i++)        //  Build the code table
    { scheme->codebits[i] = 0;
      scheme->codelens[i] = 0;
    }

  Build_Table(node+(range-1),0,0,scheme->codebits,scheme->codelens);

  if (inscheme != NULL)            //  Set scheme type and if truncated (2), map truncated codes
    { scheme->type = 2;            //    to code and length for 255
      for (i = 0; i < 255; i++)
        if (inscheme->codelens[i] > HUFF_CUTOFF || scheme->codelens[i] > HUFF_CUTOFF)
          { scheme->codelens[i] = scheme->codelens[255];
            scheme->codebits[i] = scheme->codebits[255];
          }
    }
  else
    { scheme->type = 0;
      for (i = 0; i < 256; i++)
        { if (scheme->codelens[i] > HUFF_CUTOFF)
            scheme->type = 1;
        }
    }

  return (scheme);
}

  //  Write the code table to out.

static void Write_Scheme(HScheme *scheme, FILE *out)
{ int     i;
  uint8   x;
  uint32 *bits;
  int    *lens;

  lens = scheme->codelens;
  bits = scheme->codebits;
  x = scheme->type;
  fwrite(&x,1,1,out);
  for (i = 0; i < 256; i++)
    { x = (uint8) (lens[i]);
      fwrite(&x,1,1,out);
      if (x > 0)
        fwrite(bits+i,sizeof(uint32),1,out);
    }
}

  //  Encode read[0..rlen-1] according to scheme and write to out

static void Encode(HScheme *scheme, FILE *out, uint8 *read, int rlen)
{ uint32  x, c, ocode;
  int     n, k, olen, llen;
  int    *nlens;
  uint32 *nbits;
  uint32  nspec;
  int     nslen;

  nlens = scheme->codelens;
  nbits = scheme->codebits;

  if (scheme->type == 2)
    { nspec = nbits[255];
      nslen = nlens[255];
    }
  else
    nspec = nslen = 0x7fffffff;

#define OCODE(L,C)				\
{ int    len  = olen + (L);			\
  uint32 code = (C);				\
						\
  llen = olen;					\
  if (len >= 32)				\
    { olen   = len-32;				\
      ocode |= (code >> olen);			\
      fwrite(&ocode,sizeof(uint32),1,out);	\
      if (olen > 0)				\
        ocode = (code << (32-olen));		\
      else					\
        ocode = 0;				\
    } 						\
  else						\
    { olen   = len;				\
      ocode |= (code << (32-olen));;		\
    }						\
}

  llen  = 0;
  olen  = 0;
  ocode = 0;
  for (k = 0; k < rlen; k++)
    { x = read[k];
      n = nlens[x];
      c = nbits[x];
      OCODE(n,c);
      if (c == nspec && n == nslen)
        OCODE(8,x);
    }

  if (olen > 0)                              //  Tricky: must pad so decoder does not read past
    { fwrite(&ocode,sizeof(uint32),1,out);   //    last integer int the coded output.
      if (llen > 16 && olen > llen)
        fwrite(&ocode,sizeof(uint32),1,out);
    }
  else if (llen > 16)
    fwrite(&ocode,sizeof(uint32),1,out);
}

  //  Encode read[0..rlen-1] according to non-rchar table neme, and run-length table reme for
  //    runs of rchar characters.  Write to out.

static void Encode_Run(HScheme *neme, HScheme *reme, FILE *out, uint8 *read, int rlen, int rchar)
{ uint32  x, c, ocode;
  int     n, h, k, olen, llen;
  int    *nlens, *rlens;
  uint32 *nbits, *rbits;
  uint32  nspec, rspec;
  int     nslen, rslen;

  nlens = neme->codelens;
  nbits = neme->codebits;
  rlens = reme->codelens;
  rbits = reme->codebits;

  if (neme->type == 2)
    { nspec = nbits[255];
      nslen = nlens[255];
    }
  else
    nspec = nslen = 0x7fffffff;

  rspec = rbits[255];
  rslen = rlens[255];

  llen  = 0;
  olen  = 0;
  ocode = 0;
  k     = 0;
  while (k < rlen)
    { h = k;
      while (k < rlen && read[k] == rchar)
        k += 1;
      if (k-h >= 255)
        x = 255;
      else
        x = k-h;
      n = rlens[x];
      c = rbits[x];
      OCODE(n,c);
      if (c == rspec && n == rslen)
        OCODE(16,k-h);
      if (k < rlen)
        { x = read[k];
          n = nlens[x];
          c = nbits[x];
          OCODE(n,c);
          if (c == nspec && n == nslen)
            OCODE(8,x);
          k += 1;
        }
    }

  if (olen > 0)
    { fwrite(&ocode,sizeof(uint32),1,out);
      if (llen > 16 && olen > llen)
        fwrite(&ocode,sizeof(uint32),1,out);
    }
  else if (llen > 16)
    fwrite(&ocode,sizeof(uint32),1,out);
}


/*******************************************************************************************
 *
 *  Reader, Run Histogrammmer, and DelTag Compacter routines
 *
 ********************************************************************************************/

static char  *Read;
static int    Rmax;
static int    Nline;

//  If nlines == 1 trying to read a single header, nlines = 5 trying to read 5 QV/fasta lines
//    for a sequence.  Place line j at Read+j*Rmax and the length of every line is returned
//    unless eof occurs while trying to read a header, in which case return -1.

static int Read_Seq(FILE *input, int nlines)
{ int   i, rlen;
  char *other;

  Nline += 1;
  if (fgets(Read,Rmax,input) == NULL)
    { if (nlines == 1)
        return (-1);
      else
        { fprintf(stderr,"Line %d: incomplete last entry of .quiv file\n",Nline);
          exit (1);
        }
    }
  rlen = strlen(Read);
  while (Read[rlen-1] != '\n')
    { Rmax = 1.4*Rmax + MIN_BUFFER;
      Read = (char *) Guarded_Alloc(Read,5*Rmax);
      if (fgets(Read+rlen,Rmax-rlen,input) == NULL)
        { fprintf(stderr,"Line %d: Last line does not end with a newline !\n",Nline);
          exit (1);
        }
      rlen += strlen(Read+rlen);
    }
  other = Read;
  for (i = 1; i < nlines; i++)
    { other += Rmax;
      Nline += 1;
      if (fgets(other,Rmax,input) == NULL)
        { fprintf(stderr,"Line %d: incomplete last entry of .quiv file\n",Nline);
          exit (1);
        }
      if (rlen != (int) strlen(other))
        { fprintf(stderr,"Line %d: Lines for an entry are not the same length\n",Nline);
          exit (1);
        }
    }
  return (rlen-1);
}

//  Histogram runlengths of symbol runChar in stream[0..rlen-1] into run.

static void Histogram_Seqs(uint64 *hist, uint8 *stream, int rlen)
{ int k;

  for (k = 0; k < rlen; k++)
    hist[stream[k]] += 1;
}

static void Histogram_Runs(uint64 *run, uint8 *stream, int rlen, int runChar)
{ int k, h;

  k = 0;
  while (k < rlen)
    { h = k;
      while (k < rlen && stream[k] == runChar)
        k += 1;
      if (k-h >= 256)
        run[255] += 1;
      else
        run[k-h] += 1;
      if (k < rlen)
        k += 1;
    }
}

//  Keep only the symbols in tags[0..rlen-1] for which qvs[k] != rchar and
//    return the # of symbols kept.

static int Pack_Tag(char *tags, char *qvs, int rlen, int rchar)
{ int j, k;

  j = 0;
  for (k = 0; k < rlen; k++)
    { if (qvs[k] != rchar)
        tags[j++] = tags[k];
    }
  return (j);
}

/*******************************************************************************************
 *
 *  Main routine
 *
 ********************************************************************************************/

int main(int argc, char* argv[])

{ uint64    *delHist, *insHist, *mrgHist, *subHist, *delRun, *subRun;
  int        VERBOSE;
  int        KEEP;
  int        LOSSY;

  Program_Name = argv[0];

  { int   i, j, k;

    VERBOSE = 1;
    KEEP    = 0;
    LOSSY   = 0;

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        for (k = 1; argv[i][k] != '\0'; k++)
          if (argv[i][k] == 'S')
            VERBOSE = 0;
          else if (argv[i][k] == 'k')
            KEEP = 1;
          else if (argv[i][k] == 'l')
            LOSSY = 1;
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

  // Create histograms needed to determine the various huffman schemes

  delHist = (uint64 *) Guarded_Alloc(NULL,sizeof(uint64)*256);
  delRun  = (uint64 *) Guarded_Alloc(NULL,sizeof(uint64)*256);
  insHist = (uint64 *) Guarded_Alloc(NULL,sizeof(uint64)*256);
  mrgHist = (uint64 *) Guarded_Alloc(NULL,sizeof(uint64)*256);
  subHist = (uint64 *) Guarded_Alloc(NULL,sizeof(uint64)*256);
  subRun  = (uint64 *) Guarded_Alloc(NULL,sizeof(uint64)*256);

  // For each .quiva file to be compressed:

  { int i;

    Rmax = MIN_BUFFER;
    Read = (char *) Guarded_Alloc(NULL,Rmax);

    for (i = 1; i < argc; i++)

      { HScheme *delScheme, *insScheme, *mrgScheme, *subScheme;
        HScheme *dRunScheme, *sRunScheme;
        FILE    *input, *output;
        int      delChar, subChar;
        uint64   totChar;
        char    *full;

        //  Open the .quiva file for reading and create & open .dexqv file for writing

        { char *path;
          int   epos;

          path = argv[i];
          full = (char *) Guarded_Alloc(NULL,strlen(path)+20);
          epos = strlen(path);
          if (epos >= 6 && strcasecmp(path+(epos-6),".quiva") == 0)
            strcpy(full,path);
          else
            { epos += 6;
              sprintf(full,"%s.quiva",path);
            }

          input = Guarded_Fopen(full,"r");

          if (VERBOSE)
            { fprintf(stderr,"Processing '%s' ...",full);
              fflush(stderr);
            }

          strcpy(full+(epos-6),".dexqv");
          output = Guarded_Fopen(full,"w");

          strcpy(full+(epos-6),".quiva");
        }

        //  Zero the histograms

        bzero(delHist,sizeof(uint64)*256);
        bzero(delRun,sizeof(uint64)*256);
        bzero(mrgHist,sizeof(uint64)*256);
        bzero(insHist,sizeof(uint64)*256);
        bzero(subHist,sizeof(uint64)*256);
        bzero(subRun,sizeof(uint64)*256);

        //  Make a sweep through the .quiva entries, histogramming the relevant things
        //    and figuring out the run chars for the deletion and substition streams

        totChar = 0;
        delChar = -1;
        subChar = -1;
        Nline   = 0;
        while (1)
          { int    rlen, well, beg, end, qv;
            char  *slash;

            rlen = Read_Seq(input,1);
            if (rlen < 0)
              break;
            if (rlen == 0 || Read[0] != '@')
              { fprintf(stderr,"Line %d: Header in quiv file is missing\n",Nline);
                exit (1);
              }
            slash = index(Read+1,'/');
            if (slash == NULL)
  	    { fprintf(stderr,"%s: Line %d: Header line incorrectly formatted ?\n",
                               Program_Name,Nline);
                exit (1);
              }
            if (sscanf(slash+1,"%d/%d_%d RQ=0.%d\n",&well,&beg,&end,&qv) != 4)
              { fprintf(stderr,"%s: Line %d: Header line incorrectly formatted ?\n",
                               Program_Name,Nline);
                exit (1);
              }

            rlen = Read_Seq(input,5);

            Histogram_Seqs(delHist,(uint8 *) (Read),rlen);
            Histogram_Seqs(insHist,(uint8 *) (Read+2*Rmax),rlen);
            Histogram_Seqs(mrgHist,(uint8 *) (Read+3*Rmax),rlen);
            Histogram_Seqs(subHist,(uint8 *) (Read+4*Rmax),rlen);

            if (delChar < 0)
              { int   k;
                char *del = Read+Rmax;

                for (k = 0; k < rlen; k++)
                  if (del[k] == 'N')
                    { delChar = Read[k];
                      break;
                    }
              }
            if (delChar >= 0)
              Histogram_Runs( delRun,(uint8 *) (Read),rlen,delChar);
            totChar += rlen;
            if (subChar < 0)
              { if (totChar >= 100000)
                  { int k;

                    subChar = 0;
                    for (k = 1; k < 256; k++)
                      if (subHist[k] > subHist[delChar])
                        subChar = k;
                  }
              }
            if (subChar >= 0)
              Histogram_Runs( subRun,(uint8 *) (Read+4*Rmax),rlen,subChar);
          }

        //  Check whether using a subtitution run char is a win

        if (totChar < 200000 || subHist[subChar] < .5*totChar)
          subChar = -1;

        //  If lossy encryption is enabled then scale insertions and merge QVs.

        if (LOSSY)
          { int k;

            for (k = 0; k < 256; k += 2)
              { insHist[k] += insHist[k+1];
                insHist[k+1] = 0;
              }

            for (k = 0; k < 256; k += 4)
              { mrgHist[k] += mrgHist[k+1];
                mrgHist[k] += mrgHist[k+2];
                mrgHist[k] += mrgHist[k+3];
                mrgHist[k+1] = 0;
                mrgHist[k+2] = 0;
                mrgHist[k+3] = 0;
              }
          }

        //  Build a Huffman scheme for each stream entity from the histograms

#define SCHEME_MACRO(meme,hist,label,bits)	\
  scheme = Huffman( (hist), NULL);		\
  if (scheme->type)				\
    { (meme) = Huffman( (hist), scheme);	\
      free(scheme);				\
    }						\
  else						\
    (meme) = scheme;

#ifdef DEBUG

#define MAKE_SCHEME(meme,hist,label,bits)	\
  SCHEME_MACRO(meme,hist,label,bits)		\
  printf("\n%s\n", (label) );			\
  Print_Histogram( (hist));			\
  Print_Table( (meme), (hist), (bits));	

#else

#define MAKE_SCHEME(meme,hist,label,bits)	\
  SCHEME_MACRO(meme,hist,label,bits)

#endif

        { HScheme *scheme;

          if (delChar < 0)
            { MAKE_SCHEME(delScheme,delHist, "Hisotgram of Deletion QVs", 8); }
          else
            { delHist[delChar] = 0;
              MAKE_SCHEME(delScheme,delHist, "Hisotgram of Deletion QVs less run char", 8);
              MAKE_SCHEME(dRunScheme,delRun, "Histogram of Deletion Runs QVs", 16);
#ifdef DEBUG
              printf("\nRun char is '%c'\n",delChar);
#endif
            }

#ifdef DEBUG
          { int    k;
            uint64 count;

            count = 0;
            for (k = 0; k < 256; k++)
              count += delHist[k];
            printf("\nDelTag will require %lld bytes\n",count/4);
          }
#endif

          MAKE_SCHEME(insScheme,insHist, "Hisotgram of Insertion QVs", 8);
          MAKE_SCHEME(mrgScheme,mrgHist, "Hisotgram of Merge QVs", 8);

          if (subChar < 0)
            { MAKE_SCHEME(subScheme,subHist, "Hisotgram of Subsitution QVs", 8); }
          else
            { subHist[subChar] = 0;
              MAKE_SCHEME(subScheme,subHist, "Hisotgram of Subsitution QVs less run char", 8);
              MAKE_SCHEME(sRunScheme,subRun, "Histogram of Substitution Run QVs", 16);
#ifdef DEBUG
              printf("\nRun char is '%c'\n",subChar);
#endif
            }
        }

        //  In a second pass ...

        rewind(input);

        //   Get the first header and write out the endian key, run chars, and short name

        { uint16 half;
          int    spos;
          char  *slash;

          half = 0x33cc;
          fwrite(&half,sizeof(uint16),1,output);

          if (delChar < 0)
            half = 256;
          else
            half = delChar;
          fwrite(&half,sizeof(uint16),1,output);

          if (subChar < 0)
            half = 256;
          else
            half = subChar;
          fwrite(&half,sizeof(uint16),1,output);

          Read_Seq(input,1);
          slash = index(Read,'/');
          spos = slash-Read;
          fwrite(&spos,sizeof(int),1,output);
          fwrite(Read,1,spos,output);
        }

        //   Write out the scheme tables

        Write_Scheme(delScheme,output);
        if (delChar >= 0)
          Write_Scheme(dRunScheme,output);
        Write_Scheme(insScheme,output);
        Write_Scheme(mrgScheme,output);
        Write_Scheme(subScheme,output);
        if (subChar >= 0)
          Write_Scheme(sRunScheme,output);

        //   For each entry

        { int lwell;

          lwell = 0;
          do
            { int    rlen, clen;
              int    well, beg, end, qv;
              char  *slash;
              uint16 half;
              uint8  byte;

              //  Interpret the header, encode and write out the fields

              slash = index(Read,'/');
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

              //  Get all 5 streams, compress each with its scheme, and output

              rlen = Read_Seq(input,5);
              if (delChar < 0)
                { Encode(delScheme, output, (uint8 *) Read, rlen);
                  clen = rlen;
                }
              else
                { Encode_Run(delScheme, dRunScheme, output, (uint8 *) Read, rlen, delChar);
                  clen = Pack_Tag(Read+Rmax,Read,rlen,delChar);
                }
              Compress_Read(clen,Read+Rmax);
              fwrite(Read+Rmax,1,COMPRESSED_LEN(clen),output);

              if (LOSSY)
                { uint8 *insert = (uint8 *) (Read+2*Rmax);
                  uint8 *merge  = (uint8 *) (Read+3*Rmax);
                  int    k;

                  for (k = 0; k < rlen; k++)
                    { insert[k] = ((insert[k] >> 1) << 1);
                      merge[k]  = (( merge[k] >> 2) << 2);
                    }
                }

              Encode(insScheme, output, (uint8 *) (Read+2*Rmax), rlen);
              Encode(mrgScheme, output, (uint8 *) (Read+3*Rmax), rlen);
              if (subChar < 0)
                Encode(subScheme, output, (uint8 *) (Read+4*Rmax), rlen);
              else
                Encode_Run(subScheme, sRunScheme, output, (uint8 *) (Read+4*Rmax), rlen, subChar);
            }
          while (Read_Seq(input,1) >= 0);
        }

        //  Clean up for the next iteration

        if (subChar >= 0)
          free(sRunScheme);
        free(subScheme);
        free(mrgScheme);
        free(insScheme);
        if (delChar >= 0)
          free(dRunScheme);
        free(delScheme);

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

    free(Read);
  }

  free(subRun);
  free(subHist);
  free(mrgHist);
  free(insHist);
  free(delRun);
  free(delHist);

  exit (0);
}
