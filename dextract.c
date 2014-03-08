/*******************************************************************************************
 *
 *  Dextractor: pullls requested info out of .bax.h5 files produced by Pacbio
 *
 *
 *  Author:  Martin Pippel
 *  Date  :  Dec 12, 2013
 *
 *  Author:  Gene Myers
 *  Date:    Jan 8, 2014, rewrite of the entire code base
 *
 ********************************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include <sys/stat.h>

#include <hdf5.h>

#define PHRED_OFFSET  33

static char *Usage = "[-qS] [-o[<path>]] [-l<int(500)>] [-s<int(750)>] <input:bax_h5> ...";

#define DEXTRACT

#include "shared.c"

// Exception codes

#define CANNOT_OPEN_BAX_FILE 	-1
#define BAX_BASECALL_ERR 	-2
#define BAX_DELETIONQV_ERR 	-3
#define BAX_DELETIONTAG_ERR 	-4
#define BAX_INSERTIONQV_ERR 	-5
#define BAX_MERGEQV_ERR 	-6
#define BAX_SUBSTITUTIONQV_ERR 	-7
#define BAX_QV_ERR 		-8
#define BAX_NR_EVENTS_ERR 	-9
#define BAX_REGION_ERR		-10

typedef struct
  { char   *fullName;      // full path
    char   *shortName;     // without path and file extension (used for output: fasta, fastq)
    hsize_t length;        // sum of all raw reads
    int     id;            // arbitrary file ID, used within DB creation
    int     fastq;         // if set produce a fastq file instead of a fasta file
    int     quivqv;        // if set produce a quiv file
    char   *baseCall;      // fields, that can be extracted
    char   *delQV;
    char   *delTag;
    char   *insQV;
    char   *mergeQV;
    char   *subQV;
    char   *qualityValue;
  } BaxData;

//  Initialize *the* BaxData structure

static void initBaxData(BaxData *b, int fastq, int quivqv)
{ b->fullName     = NULL;
  b->shortName    = NULL;
  b->length       = 0;
  b->id           = -1;
  b->fastq        = fastq;
  b->quivqv       = quivqv;
  b->baseCall     = NULL;
  b->delQV        = NULL;
  b->qualityValue = NULL;
}

//  Record the names of the next bax file and reset the memory buffer high-water mark

static void initBaxNames(BaxData *b, char *fname, char *hname)
{ char *beg, *end;

  beg = strrchr(hname, '/' );
  if (beg != NULL)
    beg += 1;
  else
    beg = hname;
  end = strchr(beg,'.');
  if (end != NULL)
    *end = '\0';

  free(b->shortName);

  b->id       += 1;
  b->fullName  = fname;
  b->shortName = Guarded_Strdup(beg);
  b->length    = 0;

  if (end != NULL)
    *end = '.';
}

//  Check if memory needed is above highwater mark, and if so allocate

static void ensureCapacity(BaxData *b, hsize_t length)
{ static hsize_t smax = 0;

  if (smax < length)
    { smax = 1.2*length + 10000;
      b->baseCall = (char *) Guarded_Alloc(b->baseCall, smax);
      if (b->fastq)
        b->qualityValue = (char *) Guarded_Alloc(b->qualityValue, smax);
      if (b->quivqv)
        { b->delQV   = (char *) Guarded_Alloc(b->delQV, 5ll*smax);
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
 
  H5Eset_auto(H5E_DEFAULT,0,0); // silence hdf5 error stack

  file_id = H5Fopen(b->fullName, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file_id < 0)
    return (CANNOT_OPEN_BAX_FILE);

#if DEBUG
  printf("PROCESSING %s, file_id: %d\n", baxFileName, file_id);
#endif

#define FETCH(field,path,error)								\
  { field_set   = H5Dopen2(file_id, path, H5P_DEFAULT);					\
    field_space = H5Dget_space(field_set);						\
    if (field_set < 0 || field_space < 0)						\
      { H5Fclose(file_id);								\
        return (error);									\
      }											\
    if (b->length == 0)									\
      { H5Sget_simple_extent_dims(field_space, field_len, NULL);			\
        b->length = field_len[0];							\
        ensureCapacity(b, field_len[0]);						\
      }											\
    H5Dread(field_set, H5T_NATIVE_UCHAR, H5S_ALL, H5S_ALL, H5P_DEFAULT, b->field);	\
    H5Sclose(field_space);								\
    H5Dclose(field_set);								\
  }

  FETCH(baseCall,"/PulseData/BaseCalls/Basecall",BAX_BASECALL_ERR)
  if (b->fastq)
    FETCH(qualityValue,"/PulseData/BaseCalls/QualityValue",BAX_QV_ERR)
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

    { int nregions = region_num[0];
      int region[5*nregions+1];
      int roff, *hlen, *cur, h;

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
              
                      if (b->fastq)
                        { int a;
  
                          fprintf(output,"@%s/%d/%d_%d RQ=0.%d\n",b->shortName,h,ibeg,iend,qv);
  
                          ibeg += roff;
                          iend += roff;
  
                          fprintf(output,"%.*s\n", iend-ibeg, b->baseCall + ibeg);
                          fprintf(output,"+\n");
                          for (a = ibeg; a < iend; a++)
                            fputc(b->qualityValue[a]+PHRED_OFFSET,output);
                          fputc('\n',output);
                        }
                      else
                        { int a;

                          fprintf(output,">%s/%d/%d_%d RQ=0.%d\n",b->shortName,h,ibeg,iend,qv);
  
                          ibeg += roff;
                          iend += roff;
                           
                          for (a = ibeg; a < iend; a += 70)
                            if (a+70 > iend)
                              fprintf(output,"%.*s\n", iend-a, b->baseCall + a);
                            else
                              fprintf(output,"%.70s\n", b->baseCall + a);
                        }
  
                      if (b->quivqv)
                        { int a;
  
                          ibeg -= roff;
                          iend -= roff;
  
                          fprintf(qvquiv,"@%s/%d/%d_%d RQ=0.%d\n",b->shortName,h,ibeg,iend,qv);
  
                          ibeg += roff;
                          iend += roff;
  
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
{ b->length = 0;
  if (b->baseCall != NULL)
    free(b->baseCall);
  if (b->delQV != NULL)
    free(b->delQV);
  if (b->qualityValue != NULL)
    free(b->qualityValue);
  free(b->shortName);
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

  FASTQ     = 0;
  QUIVQV    = 0;
  VERBOSE   = 1;
  MIN_LEN   = 500;
  MIN_SCORE = 750;
  output    = NULL;

  Program_Name = argv[0];

  { int   i, j, k;
    char *eptr;

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            for (k = 1; argv[i][k] != '\0'; k++)
              if (argv[i][k] == 'S')
                VERBOSE = 0;
              else if (argv[i][k] == 'q')
                FASTQ = 1;
              else
                { fprintf(stderr,"%s: -%c is an illegal option\n",Program_Name,argv[i][k]);
                  exit (1);
                }
            break;
          case 's':
            MIN_SCORE = strtol(argv[i]+2,&eptr,10);
            if (*eptr != '\0' || argv[i][2] == '\0')
              { fprintf(stderr,"%s: -s argument is not an integer\n",Program_Name);
                exit (1);
              }
            if (MIN_SCORE < 0)
              { fprintf(stderr,"%s: Subread score must be non-negative (%d)\n",
                               Program_Name,MIN_SCORE);
                exit (1);
              }
            break;
          case 'l':
            MIN_LEN = strtol(argv[i]+2,&eptr,10);
            if (*eptr != '\0' || argv[i][2] == '\0')
              { fprintf(stderr,"%s: -l argument is not an integer\n",Program_Name);
                exit (1);
              }
            if (MIN_LEN < 0)
              { fprintf(stderr,"%s: Subread length must be non-negative (%d)\n",
                               Program_Name,MIN_LEN);
                exit (1);
              }
            break;
          case 'o':
            QUIVQV = 1;
            output = argv[i]+2;
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    if (argc == 1)
      { fprintf(stderr,"Usage: %s %s\n",Program_Name,Usage);
        exit (1);
      }
  }

  fileQuiv = NULL;
  if (QUIVQV)
    { char *name;
      char *eos;

      if (*output == '\0')
        { output = Guarded_Strdup(argv[1]);   //  will not be freed, OK as only once

          eos = strrchr(output,'/');
          if (eos == NULL)
            eos = output;
          eos = strchr(eos+1,'.');
          if (eos != NULL)
            *eos   = '\0';
        }
 
      name = (char *) Guarded_Alloc(NULL,strlen(output) + 20);

      if (FASTQ)
        sprintf(name,"%s.fastq",output);
      else
        sprintf(name,"%s.fasta",output);
      fileOut = Guarded_Fopen(name, "w");

      if (QUIVQV)
        { sprintf(name,"%s.quiva",output);
          fileQuiv = Guarded_Fopen(name, "w");
        }

      free(name);
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
      { char *full;
        hid_t fileID;
        int   ecode;
  
        { int   epos;
          char *input;
  
          input = argv[i];
          full = (char *) Guarded_Alloc(NULL,strlen(input) + 20);
          epos = strlen(input);
          if (epos >= 7 && strcasecmp(input+(epos-7),".bax.h5") == 0)
            strcpy(full,input);
          else
            { FILE *in;
    
              sprintf(full,"%s.bax.h5",input);
              if ((in = fopen(full,"r")) == NULL)
                { fprintf(stderr,"%s: Cannot find %s !\n",Program_Name,input);
                  exit (1);
                }
              else
                fclose(in);
            }
  
          if (QUIVQV)
            initBaxNames(&b,full,output);
          else
            initBaxNames(&b,full,input);
        }
    
        if (VERBOSE)
          { fprintf(stderr, "Fetching file : %s ...", full); fflush(stderr); }
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
              fprintf(stderr, " Skipping %s due to failure\n",full); 
            printBaxError(ecode);
          }
        else
          { if (VERBOSE)
              { fprintf(stderr, " Done\n"); fflush(stdout); }
          }
      }
  }

  freeBaxData(&b);
  if (fileOut != stdout)
    fclose(fileOut);
  if (fileQuiv != NULL)
    fclose(fileQuiv);

  exit (0);
}
