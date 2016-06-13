/*******************************************************************************************
 *
 *  Dextractor: pullls requested info out of .bax.h5 files produced by Pacbio
 *
 *
 *  Author:  Martin Pippel
 *  Date  :  Dec 12, 2013
 *
 *  Author:  Gene Myers
 *  Date:    Jan 8, 2014, redesign of the modes of operation and flags, and also the
 *               logic for extraction in writeBaxReads
 *
 ********************************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include <ctype.h>
#include <sys/stat.h>

#include <hdf5.h>

#include "DB.h"

#define LOWER_OFFSET 32
#define PHRED_OFFSET 33

static char *Usage = "[-vq] [-o[<path>]] [-l<int(500)>] [-s<int(750)>] <input:bax_h5> ...";

#define DEXTRACT

// Exception codes

#define CANNOT_OPEN_BAX_FILE   1
#define BAX_BASECALL_ERR       2
#define BAX_DEL_ERR            3
#define BAX_TAG_ERR            4
#define BAX_INS_ERR            5
#define BAX_MRG_ERR            6
#define BAX_SUB_ERR            7
#define BAX_QV_ERR             8
#define BAX_NR_EVENTS_ERR      9
#define BAX_REGION_ERR        10
#define BAX_HOLESTATUS_ERR    11

typedef struct
  { char   *fullName;      // full file path
    char   *shortName;     // without path and file extension (used in header line)
    int     fastq;         // if non-zero produce a fastq file instead of a fasta file
    int     quivqv;        // if non-zero produce a quiv file

    hsize_t numBP;         // sum of all raw read lengths
    char   *baseCall;      // 7 streams that may be extracted dependent on flag settings
    char   *delQV;
    char   *delTag;
    char   *insQV;
    char   *mergeQV;
    char   *subQV;
    char   *fastQV;

    hsize_t numZMW;        // number of wells/holes
    int    *readLen;       // length of each read in events
    char   *holeType;      // Hole type, only SEQUENCING holes are extracted

    hsize_t numHQR;        // number of regions
    int    *regions;       // region information (5 ints per entry)

    int     delLimit;     //  The Del QV associated with N's in the Del Tag

  } BaxData;

//  Initialize *the* BaxData structure

static void initBaxData(BaxData *b, int fastq, int quivqv)
{ b->fullName  = NULL;
  b->shortName = NULL;
  b->fastq     = fastq;
  b->quivqv    = quivqv;
  b->baseCall  = NULL;
  b->delQV     = NULL;
  b->delTag    = NULL;
  b->insQV     = NULL;
  b->mergeQV   = NULL;
  b->subQV     = NULL;
  b->fastQV    = NULL;
  b->readLen   = NULL;
  b->holeType  = NULL;
  b->regions   = NULL;
  b->delLimit  = 0;
}

//  Record the names of the next bax file and reset the memory buffer high-water mark

static void initBaxNames(BaxData *b, char *fname, char *hname)
{ b->fullName  = fname;
  b->shortName = hname;
  b->numBP     = 0;
  b->numZMW    = 0;
  b->numHQR    = 0;
}

//  Check if memory needed is above highwater mark, and if so allocate

static void ensureBases(BaxData *b, hsize_t len)
{ static hsize_t smax = 0;

  b->numBP = len;
  if (smax < len)
    { smax = 1.2*len + 10000;
      b->baseCall = (char *) Realloc(b->baseCall, smax, "Allocating basecall vector");
      if (b->fastq)
        b->fastQV = (char *) Realloc(b->fastQV, smax, "Allocating fastq vector");
      if (b->quivqv)
        { b->delQV   = (char *) Realloc(b->delQV, 5ll*smax, "Allocating 5 QV vectors");
          b->delTag  = b->delQV   + smax;
          b->insQV   = b->delTag  + smax;
          b->mergeQV = b->insQV   + smax;
          b->subQV   = b->mergeQV + smax;
        }
    }
}

static void ensureZMW(BaxData *b, hsize_t len)
{ static hsize_t smax = 0;

  b->numZMW = len;
  if (smax < len)
    { smax = 1.2*len + 10000;
      b->holeType = (char *) Realloc(b->holeType, smax, "Allocating hole vector");
      b->readLen  = (int *) Realloc(b->readLen , smax * sizeof(int), "Allocating event vector");
    }
}

static void ensureHQR(BaxData *b, hsize_t len)
{ static hsize_t smax = 0;

  b->numHQR = len;
  if (smax < len)
    { smax = 1.2*len + 10000;
      b->regions = (int *) Realloc(b->regions, (5ll*smax+1)*sizeof(int), "Allocating region vector");
    }
}

// Fetch the relevant contents of the current bax.h5 file and return the H5 file id.

static int getBaxData(BaxData *b)
{ hid_t   field_space;
  hid_t   field_set;
  hsize_t field_len[2];
  hid_t   file_id;
  herr_t  stat;
  int     ecode;

  H5Eset_auto(H5E_DEFAULT,0,0); // silence hdf5 error stack

  file_id = H5Fopen(b->fullName, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file_id < 0)
    return (CANNOT_OPEN_BAX_FILE);

#ifdef DEBUG
  printf("PROCESSING %s, file_id: %d\n", baxFileName, file_id);
#endif

#define GET_SIZE(path,error)									\
  { ecode = error;										\
    if ((field_set = H5Dopen2(file_id, path, H5P_DEFAULT)) < 0) goto exit0;			\
    if ((field_space = H5Dget_space(field_set)) < 0) goto exit1;				\
    H5Sget_simple_extent_dims(field_space, field_len, NULL);					\
  }

#define FETCH(field,type)									\
  { stat = H5Dread(field_set, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, b->field);			\
    H5Sclose(field_space);									\
    H5Dclose(field_set);									\
    if (stat < 0) goto exit0;									\
  }

#define CHECK_FETCH(path,error,field,type,cntr)							\
  { GET_SIZE(path,error)									\
    if (b->cntr != field_len[0]) goto exit2;							\
    FETCH(field,type)										\
  }

  GET_SIZE("/PulseData/BaseCalls/Basecall",BAX_BASECALL_ERR)
  ensureBases(b,field_len[0]);
  FETCH(baseCall,H5T_NATIVE_UCHAR)
  if (b->fastq)
    CHECK_FETCH("/PulseData/BaseCalls/QualityValue",BAX_QV_ERR,fastQV,H5T_NATIVE_UCHAR,numBP)
  if (b->quivqv)
    { CHECK_FETCH("/PulseData/BaseCalls/DeletionQV",    BAX_DEL_ERR,delQV,  H5T_NATIVE_UCHAR,numBP)
      CHECK_FETCH("/PulseData/BaseCalls/DeletionTag",   BAX_TAG_ERR,delTag, H5T_NATIVE_UCHAR,numBP)
      CHECK_FETCH("/PulseData/BaseCalls/InsertionQV",   BAX_INS_ERR,insQV,  H5T_NATIVE_UCHAR,numBP)
      CHECK_FETCH("/PulseData/BaseCalls/MergeQV",       BAX_MRG_ERR,mergeQV,H5T_NATIVE_UCHAR,numBP)
      CHECK_FETCH("/PulseData/BaseCalls/SubstitutionQV",BAX_SUB_ERR,subQV,  H5T_NATIVE_UCHAR,numBP)
    }

  GET_SIZE("/PulseData/BaseCalls/ZMW/HoleStatus",BAX_HOLESTATUS_ERR)
  ensureZMW(b,field_len[0]);
  FETCH(holeType,H5T_NATIVE_UCHAR)
  CHECK_FETCH("/PulseData/BaseCalls/ZMW/NumEvent",BAX_NR_EVENTS_ERR,readLen,H5T_NATIVE_INT,numZMW)

  GET_SIZE("/PulseData/Regions",BAX_REGION_ERR)
  ensureHQR(b,field_len[0]);
  FETCH(regions,H5T_NATIVE_INT)

  //  Find the Del QV associated with N's in the Del Tag

  if (b->quivqv)
    { hsize_t  i;
  
      for (i = 0; i < b->numBP; i++)
        if (b->delTag[i] == 'N')
          { b->delLimit = b->delQV[i];
            break;
          }
    }

  return (0);

exit2:
  H5Sclose(field_space);
exit1:
  H5Dclose(field_set);
exit0:
  H5Fclose(file_id);
  return (ecode);
}

// Find the good read invervals of the baxfile b(FileID), output the reads of length >= minLen and
//   score >= minScore to output (for the fasta or fastq part) and qvquiv (if b->quivqv is set)

static char *fasta_header = ">%s/%d/%d_%d RQ=0.%d\n";
static char *fastq_header = "@%s/%d/%d_%d RQ=0.%d\n";

static void writeBaxReads(BaxData *b, int minLen, int minScore, FILE *output, FILE* qvquiv)
{ int   nreads, *rlen;
  int   roff, *hlen, *cur, h, w;
  int   tolower;
  char *header;

  char   *baseCall;
  char   *delQV;
  char   *delTag;
  char   *insQV;
  char   *mergeQV;
  char   *subQV;
  char   *fastQV;

  baseCall = b->baseCall;
  delQV    = b->delQV;
  delTag   = b->delTag;
  insQV    = b->insQV;
  mergeQV  = b->mergeQV;
  subQV    = b->subQV;
  fastQV   = b->fastQV;

#ifdef DEBUG
  printf("printSubreadFields\n");
#endif

#define HOLE   0
#define TYPE   1
#define    ADAPTER_REGION 0
#define    INSERT_REGION  1
#define    HQV_REGION     2
#define START  2
#define FINISH 3
#define SCORE  4

  //  Find the HQV regions and output as reads according to the various output options

  tolower = isupper(b->baseCall[0]);
  if (b->fastq)
    header = fastq_header;
  else
    header = fasta_header;

  rlen    = b->readLen;
  roff    = 0;
  cur     = b->regions;
  nreads  = b->numZMW + cur[HOLE];
  hlen    = rlen - cur[HOLE];
  cur[5*b->numHQR] = nreads; 

  for (h = cur[HOLE], w = 0; h < nreads; h++, w++)
    { int *bot, *top, *hqv, *r;
      int hbeg, hend, qv;
      int ibeg, iend;

      if (hlen[h] >= minLen)
        { while (cur[HOLE] < h)
            cur += 5;
          bot = hqv = cur;
          while (cur[HOLE] <= h)
            { if (cur[TYPE] == HQV_REGION)
                hqv = cur;
              cur += 5;
            }
          top = cur-5;

          qv = hqv[SCORE];
          if (qv >= minScore)
            { hbeg = hqv[START];
              hend = hqv[FINISH];
              for (r = bot; r <= top; r += 5)
                { if (r[TYPE] != INSERT_REGION)
                    continue;

                  ibeg = r[START];
                  iend = r[FINISH];

                  if (ibeg < hbeg)
                    ibeg = hbeg;
                  if (iend > hend)
                    iend = hend;
                  if (iend - ibeg < minLen || b->holeType[w] > 0)
                    continue;

                  fprintf(output,header,b->shortName,h,ibeg,iend,qv);

                  ibeg += roff;
                  iend += roff;

                  if (tolower)
                    { int a;

                      for (a = ibeg; a < iend; a++)
                        baseCall[a] += LOWER_OFFSET;
                      if (b->quivqv)
                        for (a = ibeg; a < iend; a++)
                          delTag[a] += LOWER_OFFSET;
                    }

                  if (b->fastq)
                    { int a;

                      fprintf(output,"%.*s\n", iend-ibeg, baseCall + ibeg);
                      fprintf(output,"+\n");
                      for (a = ibeg; a < iend; a++)
                        fputc(fastQV[a]+PHRED_OFFSET,output);
                      fputc('\n',output);
                    }
                  else
                    { int a;

                      for (a = ibeg; a < iend; a += 80)
                        if (a+80 > iend)
                          fprintf(output,"%.*s\n", iend-a, baseCall + a);
                        else
                          fprintf(output,"%.80s\n", baseCall + a);
                    }

                  if (b->quivqv)
                    { int a, d;

                      fprintf(qvquiv,"@%s/%d/%d_%d RQ=0.%d\n",
                                     b->shortName,h,ibeg-roff,iend-roff,qv);

                      d = b->delLimit;
                      for (a = ibeg; a < iend; a++)
                        { if (delQV[a] == d)
                            delTag[a] = 'n';
                          delQV[a]   += PHRED_OFFSET;
                          insQV[a]   += PHRED_OFFSET;
                          mergeQV[a] += PHRED_OFFSET;
                          subQV[a]   += PHRED_OFFSET;
                        }

                      iend -= ibeg;
                      fprintf (qvquiv, "%.*s\n", iend, delQV + ibeg);
                      fprintf (qvquiv, "%.*s\n", iend, delTag + ibeg);
                      fprintf (qvquiv, "%.*s\n", iend, insQV + ibeg);
                      fprintf (qvquiv, "%.*s\n", iend, mergeQV + ibeg);
                      fprintf (qvquiv, "%.*s\n", iend, subQV + ibeg);
                    }
                }
            }
        }
      roff += hlen[h];
    }
}

//  Print an error message

static void printBaxError(int ecode)
{ fprintf(stderr,"  *** Warning ***: ");
  switch (ecode)
    { case CANNOT_OPEN_BAX_FILE:
        fprintf(stderr,"Cannot open bax file:\n");
        break;
      case BAX_BASECALL_ERR:
        fprintf(stderr,"Cannot parse /PulseData/BaseCalls/Basecall from file:\n");
        break;
      case BAX_DEL_ERR:
        fprintf(stderr,"Cannot parse /PulseData/BaseCalls/DeletionQV from file:\n");
        break;
      case BAX_TAG_ERR:
        fprintf(stderr,"Cannot parse /PulseData/BaseCalls/DeletionTag from file:\n");
        break;
      case BAX_INS_ERR:
        fprintf(stderr,"Cannot parse /PulseData/BaseCalls/InsertionQV from file:\n");
        break;
      case BAX_MRG_ERR:
        fprintf(stderr,"Cannot parse /PulseData/BaseCalls/MergeQV from file:\n");
        break;
      case BAX_SUB_ERR:
        fprintf(stderr,"Cannot parse /PulseData/BaseCalls/SubstitutionQV from file:\n");
        break;
      case BAX_QV_ERR:
        fprintf(stderr,"Cannot parse /PulseData/BaseCalls/QualityValue from file:\n");
        break;
      case BAX_NR_EVENTS_ERR:
        fprintf(stderr,"Cannot parse /PulseData/BaseCalls/ZMW/NumEvent from file:\n");
        break;
      case BAX_REGION_ERR:
        fprintf(stderr,"Cannot parse /PulseData/Regions from file:\n");
        break;
      case BAX_HOLESTATUS_ERR:
        fprintf(stderr,"Cannot parse /PulseData/BaseCalls/ZMW/HoleStatus from file:\n");
        break;
      default: 
        fprintf(stderr,"Cannot parse bax file:\n");
        break;
    }
  fflush(stderr);
}

//  Free *the* bax data structure

static void freeBaxData(BaxData *b)
{ free(b->baseCall);
  free(b->delQV);
  free(b->fastQV);
  free(b->holeType);
  free(b->readLen);
  free(b->regions);
}

int main(int argc, char* argv[])
{ char *output;
  FILE *fileOut;
  FILE *fileQuiv;

  int FASTQ;
  int QUIVQV;
  int MIN_LEN;
  int MIN_SCORE;
  int VERBOSE;

  BaxData b;

  //  Check that zlib library is present

  if ( ! H5Zfilter_avail(H5Z_FILTER_DEFLATE))
    { fprintf(stderr,"%s: zlib library is not present, check build/installation\n",Prog_Name);
      exit (1);
    }

  { int   i, j, k;
    int   flags[128];
    char *eptr;

    ARG_INIT("dextract")

    MIN_LEN   = 500;
    MIN_SCORE = 750;
    QUIVQV    = 0;
    output    = NULL;

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("qv")
            break;
          case 's':
            ARG_NON_NEGATIVE(MIN_SCORE,"Subread score threshold")
            break;
          case 'l':
            ARG_NON_NEGATIVE(MIN_LEN,"Minimum length threshold")
            break;
          case 'o':
            QUIVQV = 1;
            output = argv[i]+2;
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    VERBOSE = flags['v'];
    FASTQ   = flags['q'];

    if (argc == 1)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        exit (1);
      }
  }

  fileQuiv = NULL;
  if (QUIVQV)
    { int explicit;

      explicit = (*output != '\0');
      if ( ! explicit)
        output = Root(argv[1],NULL);

      if (FASTQ)
        fileOut = Fopen(Catenate("","",output,".fastq"), "w");
      else
        fileOut = Fopen(Catenate("","",output,".fasta"), "w");
      fileQuiv = Fopen(Catenate("","",output,".quiva"), "w");
      if (fileOut == NULL || fileQuiv == NULL)
        exit (1);

      if (explicit)
        output = Root(output,NULL);
    }
  else
    fileOut = stdout;

  if (VERBOSE)
    { fprintf(stderr, "Minimum length: %d\n", MIN_LEN);
      fprintf(stderr, "Minimum score : %d\n", MIN_SCORE);
    }

  initBaxData(&b,FASTQ,QUIVQV);

  { int i;

    for (i = 1; i < argc; i++)
      { char *root, *full, *input;
        int   ecode;

        { char *pwd;
          FILE *in;

          pwd   = PathTo(argv[i]);
          root  = Root(argv[i],".bax.h5");
          full  = Strdup(Catenate(pwd,"/",root,".bax.h5"),"Allocating full name");
          input = Root(argv[i],NULL);
          
          free(pwd);

          if ((in = fopen(full,"r")) == NULL)
            { fprintf(stderr,"%s: Cannot find %s !\n",Prog_Name,input);
              exit (1);
            }
          else
            fclose(in);
	}

        if (QUIVQV)
          initBaxNames(&b,full,output);
        else
          initBaxNames(&b,full,input);

        if (VERBOSE)
          { fprintf(stderr, "Fetching file : %s ...\n", root); fflush(stderr); }

        if ((ecode = getBaxData(&b)) == 0)
          { if (VERBOSE)
              { fprintf(stderr, "Extracting subreads ...\n"); fflush(stderr); }
            writeBaxReads(&b, MIN_LEN, MIN_SCORE, fileOut, fileQuiv);
            if (VERBOSE)
              { fprintf(stderr, "Done\n"); fflush(stdout); }
          }
        else
          { if (VERBOSE)
              fprintf(stderr, "Skipping due to failure\n"); 
            else
              fprintf(stderr, "Skipping %s due to failure\n",root); 
            printBaxError(ecode);
          }

        free(root);
        free(full);
        free(input);
      }
  }

  freeBaxData(&b);
  if (fileOut != stdout)
    fclose(fileOut);
  if (fileQuiv != NULL)
    fclose(fileQuiv);
  if (QUIVQV)
    free(output);

  exit (0);
}
