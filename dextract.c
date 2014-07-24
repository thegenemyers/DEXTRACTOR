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

#define CANNOT_OPEN_BAX_FILE   -1
#define BAX_BASECALL_ERR       -2
#define BAX_DELETIONQV_ERR     -3
#define BAX_DELETIONTAG_ERR    -4
#define BAX_INSERTIONQV_ERR    -5
#define BAX_MERGEQV_ERR        -6
#define BAX_SUBSTITUTIONQV_ERR -7
#define BAX_QV_ERR             -8
#define BAX_NR_EVENTS_ERR      -9
#define BAX_REGION_ERR        -10

typedef struct
  { char   *fullName;      // full file path
    char   *shortName;     // without path and file extension (used in header line)
    hsize_t length;        // sum of all raw read lengths
    int     fastq;         // if non-zero produce a fastq file instead of a fasta file
    int     quivqv;        // if non-zero produce a quiv file
    char   *baseCall;      // 7 streams that may be extracted dependent on flag settings
    char   *delQV;
    char   *delTag;
    char   *insQV;
    char   *mergeQV;
    char   *subQV;
    char   *fastQV;
  } BaxData;

//  Initialize *the* BaxData structure

static void initBaxData(BaxData *b, int fastq, int quivqv)
{ b->fullName  = NULL;
  b->shortName = NULL;
  b->fastq     = fastq;
  b->quivqv    = quivqv;
  b->baseCall  = NULL;
  b->delQV     = NULL;
  b->fastQV    = NULL;
}

//  Record the names of the next bax file and reset the memory buffer high-water mark

static void initBaxNames(BaxData *b, char *fname, char *hname)
{ b->fullName  = fname;
  b->shortName = hname;
  b->length    = 0;
}

//  Check if memory needed is above highwater mark, and if so allocate

static void ensureCapacity(BaxData *b)
{ static hsize_t smax = 0;

  if (smax < b->length)
    { smax = 1.2*b->length + 10000;
      b->baseCall = (char *) Realloc(b->baseCall, smax, "Allocating basecall vector");
      if (b->fastq)
        b->fastQV = (char *) Realloc(b->fastQV, smax, "Allocating fastq vector");
      if (b->quivqv)
        { b->delQV   = (char *) Realloc(b->delQV, 5ll*smax, "Allocating QV vector");
          b->delTag  = b->delQV + smax;
          b->insQV   = b->delTag + smax;
          b->mergeQV = b->insQV + smax;
          b->subQV   = b->mergeQV + smax;
        }
    }
}

// Fetch the relevant contents of the current bax.h5 file and return the H5 file id.

static hid_t getBaxData(BaxData *b)
{ hid_t   field_space;
  hid_t   field_set;
  hsize_t field_len[1];
  hid_t   file_id;
  herr_t  stat;

  H5Eset_auto(H5E_DEFAULT,0,0); // silence hdf5 error stack

  file_id = H5Fopen(b->fullName, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file_id < 0)
    return (CANNOT_OPEN_BAX_FILE);

#if DEBUG
  printf("PROCESSING %s, file_id: %d\n", baxFileName, file_id);
#endif

#define FETCH(field,path,error)									\
  { field_set   = H5Dopen2(file_id, path, H5P_DEFAULT);						\
    field_space = H5Dget_space(field_set);							\
    if (field_set < 0 || field_space < 0)							\
      { H5Fclose(file_id);									\
        return (error);										\
      }												\
    H5Sget_simple_extent_dims(field_space, field_len, NULL);					\
    if (b->length == 0)										\
      { b->length = field_len[0];								\
        ensureCapacity(b);									\
      }												\
    else											\
      { if (b->length != field_len[0])								\
          return (error);									\
      }												\
    stat = H5Dread(field_set, H5T_NATIVE_UCHAR, H5S_ALL, H5S_ALL, H5P_DEFAULT, b->field);	\
    if (stat < 0)										\
      return (error);										\
    H5Sclose(field_space);									\
    H5Dclose(field_set);									\
  }

  FETCH(baseCall,"/PulseData/BaseCalls/Basecall",BAX_BASECALL_ERR)
  if (b->fastq)
    FETCH(fastQV,"/PulseData/BaseCalls/QualityValue",BAX_QV_ERR)
  if (b->quivqv)
    { FETCH(delQV,"/PulseData/BaseCalls/DeletionQV",BAX_DELETIONQV_ERR)
      FETCH(delTag,"/PulseData/BaseCalls/DeletionTag",BAX_DELETIONTAG_ERR)
      FETCH(insQV,"/PulseData/BaseCalls/InsertionQV",BAX_INSERTIONQV_ERR)
      FETCH(mergeQV,"/PulseData/BaseCalls/MergeQV",BAX_MERGEQV_ERR)
      FETCH(subQV,"/PulseData/BaseCalls/SubstitutionQV",BAX_SUBSTITUTIONQV_ERR)
    }

  return (file_id);
}

// Find the good read invervals of the baxfile b(FileID), output the reads of length >= minLen and
//   score >= minScore to output (for the fasta or fastq part) and qvquiv (if b->quivqv is set)

static char *fasta_header = ">%s/%d/%d_%d RQ=0.%d\n";
static char *fastq_header = "@%s/%d/%d_%d RQ=0.%d\n";

static int writeBaxReads(BaxData *b, hid_t FileID,
                         int minLen, int minScore,
                         FILE *output, FILE* qvquiv)
{ hid_t   read_set, read_space;
  hid_t   region_set, region_space;
  hsize_t read_num[1],  region_num[2];

#if DEBUG
  printf("printSubreadFields\n");
#endif

  if(FileID < 0)
    return (CANNOT_OPEN_BAX_FILE);

  // Get read lengths out of the bax file

  read_set   = H5Dopen2(FileID, "/PulseData/BaseCalls/ZMW/NumEvent", H5P_DEFAULT);
  read_space = H5Dget_space(read_set);
  if (read_space < 0 || read_set < 0)
    return (BAX_NR_EVENTS_ERR);
  if (H5Sget_simple_extent_dims(read_space, read_num, NULL) < 0)
    { H5Sclose(read_space);
      H5Dclose(read_set);
      return (BAX_NR_EVENTS_ERR);
    }

  { int nreads = read_num[0];
    int rlen[nreads];

    if (H5Dread(read_set, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, rlen) < 0)
      { H5Sclose(read_space);
        H5Dclose(read_set);
        return (BAX_NR_EVENTS_ERR);
      }
    H5Sclose(read_space);
    H5Dclose(read_set);

    // Get the region annotations out of the bax file

    region_set   = H5Dopen2(FileID, "/PulseData/Regions", H5P_DEFAULT);
    region_space = H5Dget_space(region_set);
    if (region_set < 0 || region_space < 0)
      return (BAX_REGION_ERR);
    if (H5Sget_simple_extent_dims(region_space, region_num, NULL)<0)
      { H5Sclose(region_space);
        H5Dclose(region_set);
        return (BAX_REGION_ERR);
      }

    { int   nregions = region_num[0];
      int   region[5*nregions+1];
      int   roff, *hlen, *cur, h;
      int   tolower;
      char *header;

      if (H5Dread(region_set, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, region) < 0)
        { H5Sclose(region_space);
          H5Dclose(region_set);
          return (BAX_REGION_ERR);
        }
      H5Sclose(region_space);
      H5Dclose(region_set);

#define HOLE   0
#define TYPE   1
#define    ADAPTER_REGION 0
#define    INSERT_REGION  1
#define    HQV_REGION     2
#define START  2
#define FINISH 3
#define SCORE  4

      //  Find the HQV regions and output as reads to the various output options

      tolower = isupper(b->baseCall[0]);
      if (b->fastq)
        header = fastq_header;
      else
        header = fasta_header;

      roff    = 0;
      cur     = region;
      nreads += cur[HOLE];
      hlen    = rlen - cur[HOLE];
      region[5*nregions] = nreads; 
      for (h = cur[HOLE]; h < nreads; h++)
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
                      if (iend - ibeg < minLen)
                        continue;

                      fprintf(output,header,b->shortName,h,ibeg,iend,qv);

                      ibeg += roff;
                      iend += roff;

                      if (tolower)
                        { int a;

                          for (a = ibeg; a < iend; a++)
                            b->baseCall[a] += LOWER_OFFSET;
                          if (b->quivqv)
                            for (a = ibeg; a < iend; a++)
                              b->delTag[a] += LOWER_OFFSET;
                        }

                      if (b->fastq)
                        { int a;

                          fprintf(output,"%.*s\n", iend-ibeg, b->baseCall + ibeg);
                          fprintf(output,"+\n");
                          for (a = ibeg; a < iend; a++)
                            fputc(b->fastQV[a]+PHRED_OFFSET,output);
                          fputc('\n',output);
                        }
                      else
                        { int a;

                          for (a = ibeg; a < iend; a += 80)
                            if (a+80 > iend)
                              fprintf(output,"%.*s\n", iend-a, b->baseCall + a);
                            else
                              fprintf(output,"%.80s\n", b->baseCall + a);
                        }

                      if (b->quivqv)
                        { int a;

                          fprintf(qvquiv,"@%s/%d/%d_%d RQ=0.%d\n",
                                         b->shortName,h,ibeg-roff,iend-roff,qv);

                          for (a = ibeg; a < iend; a++)
                            fputc(b->delQV[a]+PHRED_OFFSET,qvquiv);
                          fputc('\n',qvquiv);

                          fprintf (qvquiv, "%.*s\n", iend-ibeg, b->delTag + ibeg);

                          for (a = ibeg; a < iend; a++)
                            fputc(b->insQV[a]+PHRED_OFFSET,qvquiv);
                          fputc('\n',qvquiv);

                          for (a = ibeg; a < iend; a++)
                            fputc(b->mergeQV[a]+PHRED_OFFSET,qvquiv);
                          fputc('\n',qvquiv);

                          for (a = ibeg; a < iend; a++)
                            fputc(b->subQV[a]+PHRED_OFFSET,qvquiv);
                          fputc('\n',qvquiv);
                        }
                    }
                }
            }
          roff += hlen[h];
        }
    }
  }

  return (1);
}

//  Print an error message

static void printBaxError(int errorCode)
{ fprintf(stderr,"  *** Warning ***: ");
  switch (errorCode)
    { case CANNOT_OPEN_BAX_FILE:
        fprintf(stderr,"Cannot open bax file:\n");
        break;
      case BAX_BASECALL_ERR:
        fprintf(stderr,"Cannot parse /PulseData/BaseCalls/Basecall from file:\n");
        break;
      case BAX_DELETIONQV_ERR:
        fprintf(stderr,"Cannot parse /PulseData/BaseCalls/DeletionQV from file:\n");
        break;
      case BAX_DELETIONTAG_ERR:
        fprintf(stderr,"Cannot parse /PulseData/BaseCalls/DeletionTag from file:\n");
        break;
      case BAX_INSERTIONQV_ERR:
        fprintf(stderr,"Cannot parse /PulseData/BaseCalls/InsertionQV from file:\n");
        break;
      case BAX_MERGEQV_ERR:
        fprintf(stderr,"Cannot parse /PulseData/BaseCalls/MergeQV from file:\n");
        break;
      case BAX_SUBSTITUTIONQV_ERR:
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
      default: 
        fprintf(stderr,"Cannot parse bax file:\n");
        break;
    }
  fflush(stderr);
}

//  Free *the* bax data structure

static void freeBaxData(BaxData *b)
{ if (b->baseCall != NULL)
    free(b->baseCall);
  if (b->delQV != NULL)
    free(b->delQV);
  if (b->fastQV != NULL)
    free(b->fastQV);
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
    { if (*output == '\0')
        output = Root(argv[1],NULL);
      else
        output = Strdup(output,"Allocating output file root");
      if (FASTQ)
        fileOut = Fopen(Catenate("","",output,".fastq"), "w");
      else
        fileOut = Fopen(Catenate("","",output,".fasta"), "w");
      fileQuiv = Fopen(Catenate("","",output,".quiva"), "w");
      if (fileOut == NULL || fileQuiv == NULL)
        exit (1);
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
        hid_t fileID;
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
          { fprintf(stderr, "Fetching file : %s ...", root); fflush(stderr); }
        ecode = fileID = getBaxData(&b);

        if (ecode > 0)
          { if (VERBOSE)
              { fprintf(stderr, " Extracting subreads ..."); fflush(stderr); }
            ecode = writeBaxReads(&b, fileID, MIN_LEN, MIN_SCORE, fileOut, fileQuiv);

            H5Fclose(fileID);
          }

        if (ecode < 0)
          { if (VERBOSE)
              fprintf(stderr, " Skipping due to failure\n"); 
            else
              fprintf(stderr, " Skipping %s due to failure\n",root); 
            printBaxError(ecode);
          }
        else
          { if (VERBOSE)
              { fprintf(stderr, " Done\n"); fflush(stdout); }
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
