/*******************************************************************************************
 *
 *  SAM/BAM reader & pacbio extractor
 *    Uses zlib file IO routines to read either SAM or BAM encoding and extracts just the
 *    information needed for by the Dazzler (sequence, fasta header, well, beg & end pulse,
 *    per base snr, and pulse width sequence).
 *
 *  Author:  Gene Myers
 *  Date  :  Oct. 9, 2016
 *
 ********************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <zlib.h>

#include "sam.h"
#include "DB.h"

 // Big to little endian converters

#define SWAP(i,j)  (x = v[i], v[i] = v[j], v[j] = x)

static void flip_short(void *w)
{ uint8 *v = (uint8 *) w;
  uint8  x;

  SWAP(0,1);
}

static void flip_int(void *w)
{ uint8 *v = (uint8 *) w;
  uint8  x;

  SWAP(0,3);
  SWAP(1,2);
}

static void flip_double(void *w)
{ uint8 *v = (uint8 *) w;
  uint8  x;

  SWAP(0,7);
  SWAP(1,6);
  SWAP(2,5);
  SWAP(3,4);
}

 // Universal data buffer for all routines

static uint8 *data = NULL;
static int    dmax = 0;

static int make_room(int len)
{ if (len > dmax)
    { if (dmax == 0)
        dmax = 8192;
      else
        dmax *= 2;
      if (dmax >= 0x7fffffff)
        { fprintf(stderr,"%s: More than MAX_INT memory requested ? Corrupt BAM?\n",Prog_Name);
          return (1);
        }
      while (len > dmax)
        dmax *= 2;
      if (dmax >= 0x7fffffff)
        { fprintf(stderr,"%s: More than MAX_INT memory requested ? Corrupt BAM?\n",Prog_Name);
          return (1);
        }
      data = realloc(data,dmax);
      if (data == NULL)
        { fprintf(stderr,"%s: Could not (re)allocate %d bytes of memory\n",Prog_Name,dmax);
          return (1);
        }
    }
  return (0);
}

 // Universal sam Record for all routines

static samRecord theR;
static int       rmax = 0;  //  maximum length for a sequence

static int init_record(int len)
{ if (len > rmax)
    { if (rmax == 0)
        theR.seq = NULL;
      rmax = 1.2*len + 1000;
      theR.seq = realloc(theR.seq,rmax);
      if (theR.seq == NULL)
        { fprintf(stderr,"%s: Could not (re)allocate %d bytes of memory\n",Prog_Name,2*rmax);
          return (1);
        }
    }
  return (0);
}


/*******************************************************************************************
 *
 *  FILE HANDLING
 *
 ********************************************************************************************/

samFile *sam_open(char *name)
{ samFile  *sf;
  gzFile    ptr = NULL;
  int       one  = 1;

  sf = (samFile *) malloc(sizeof(samFile));
  if (sf == NULL)
    return (NULL);

  sf->name = strdup(name);
  if (sf->name == NULL)
    goto error;

  if (strcmp(name,"-") == 0)
    ptr = gzdopen(fileno(stdin),"r");
  else
    ptr = gzopen(name,"r");
  if (ptr == NULL)
    goto error;
  gzbuffer(ptr,0x10000);

  if (gzdirect(ptr))
    sf->format = sam;
  else
    sf->format = bam;

  sf->ptr    = ptr;
  sf->is_big = ( *((char *) (&one)) == 0);
  sf->nline  = 0;

  return (sf);

error:
  free(sf->name);
  free(sf);
  if (ptr != NULL)
    gzclose(ptr);
  return (NULL);
}

int sam_eof(samFile *sf)
{ int c;

  c = gzgetc(sf->ptr);
  gzungetc(c,sf->ptr);
  return (c < 0);
}

int sam_close(samFile *sf)
{ int ret;

  ret = gzclose(sf->ptr);
  free(sf->name);
  free(sf);
  return (ret != Z_OK);
}

static int sam_getline(samFile *sf, int clen)
{ gzFile file = sf->ptr;

  ++sf->nline;

  if (make_room(8192))
    return (-1);
  while (gzgets(file,(char *) (data+clen),dmax-clen) != NULL)
    { clen += strlen((char *) (data+clen));
      if (data[clen-1] == '\n')
        return (clen);
      dmax *= 2;
      data = realloc(data,dmax);
      if (data == NULL)
        { fprintf(stderr,"%s: Could not (re)allocate %d bytes of memory\n",Prog_Name,dmax);
          return (-1);
        }
    }
  if (gzeof(file))
    return (0);
  else
    { fprintf(stderr,"%s: Could not get a line from %s\n",Prog_Name,sf->name);
      return (-1);
    }
}


/*******************************************************************************************
 *
 *  HEADER PROCESSING
 *
 ********************************************************************************************/

static int bam_header_read(samFile *sf)
{ gzFile file = sf->ptr;
  int     nlen, ncnt, tlen;
  int     i;

  // read "BAM\1"

  { int  ret;
    char buf[4];

    ret = gzread(file, buf, 4);
    if (ret != 4 || strncmp(buf, "BAM\1", 4) != 0)
      { fprintf(stderr, "%s: Corrupted BAM header\n",Prog_Name);
        return (1);
      }
  }

  // read plain text

  if (gzread(file, &tlen, 4) != 4)
    goto IO_error;
  if (sf->is_big)
    flip_int(&tlen);

  if (make_room(tlen+1))
    return (1);
  if (gzread(file, data, tlen) != tlen)
    goto IO_error;
  data[tlen++] = 0;              // make sure it is NULL terminated

  //  read through number of reference sequences

  if (gzread(file, &ncnt, 4) != 4)
    goto IO_error;
  if (sf->is_big)
    flip_int(&ncnt);
  if (ncnt < 0)
    goto corrupted;

  // read through reference sequence names and lengths

  for (i = 0; i < ncnt; i++)
    { if (gzread(file, &nlen, 4) != 4)
        goto IO_error;
      if (sf->is_big)
        flip_int(&nlen);
      if (nlen <= 0)
        goto corrupted;
      if (make_room(tlen+nlen+5))
        return (1);
      if (gzread(file, data+(tlen+1), nlen+4) != nlen+4)
        goto IO_error;
    }
  return (0);

IO_error:
  fprintf(stderr, "%s: truncated bam header\n",Prog_Name);
  return (1);

corrupted:
  fprintf(stderr, "%s: invalid BAM binary header\n",Prog_Name);
  return (1);
}

static int sam_header_read(samFile *sf)
{ int dlen, c;

  dlen = 0;
  while (1)
    { c = gzgetc(sf->ptr);
      if (c == EOF)
        break;
      gzungetc(c,sf->ptr);
      if (c != '@')
        break;
      dlen = sam_getline(sf,dlen);
      if (dlen < 0)
        return (1);
    }
  if (dlen == 0)
    { if (make_room(1))
        return (1);
      data[0] = '\0';
    }
  return (0);
}

 //  Sequence conversion tables

static char *INT_2_IUPAC      = "=acmgrsvtwyhkdbn";
static char  INT_2_NUMBER[16] = { 0, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0 };

static char  IUPAC_2_NUMBER[256] =
  { 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0,
    0, 1, 2, 3, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0,

    0, 0, 1, 1, 0, 0, 0, 2,   0, 0, 0, 2, 0, 0, 0, 0,
    0, 0, 0, 1, 3, 0, 0, 0,   0, 1, 0, 0, 0, 0, 0, 0,
    0, 0, 1, 1, 0, 0, 0, 2,   0, 0, 0, 2, 0, 0, 0, 0,
    0, 0, 0, 1, 3, 0, 0, 0,   0, 1, 0, 0, 0, 0, 0, 0,
  };

static char  IUPAC_2_DNA[256] =
  { 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a',   'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a',
    'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a',   'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a',
    'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a',   'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a',
    'a', 'c', 'g', 't', 'a', 'a', 'a', 'a',   'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a',

    'a', 'a', 'c', 'c', 'a', 'a', 'a', 'g',   'a', 'a', 'a', 'g', 'a', 'a', 'a', 'a',
    'a', 'a', 'a', 'c', 't', 'a', 'a', 'a',   'a', 'c', 'a', 'a', 'a', 'a', 'a', 'a',
    'a', 'a', 'c', 'c', 'a', 'a', 'a', 'g',   'a', 'a', 'a', 'g', 'a', 'a', 'a', 'a',
    'a', 'a', 'a', 'c', 't', 'a', 'a', 'a',   'a', 'c', 'a', 'a', 'a', 'a', 'a', 'a',
  };
 

static char *SEQ_CONVERT;  //  1 of 4 tables abave: sam vs. bam, numberic vs. alpha
static int   ARR_CONVERT;

int sam_header_process(samFile *sf, int numeric)
{ char *desc, *eol, *eod, *subs, *pw;
  int   status;

  if (sf->format == bam)
    { if (bam_header_read(sf))
        return (-1);
    }
  else
    { if (sam_header_read(sf))
        return (-1);
    }

  if (numeric)
    { if (sf->format == sam)
        SEQ_CONVERT = IUPAC_2_NUMBER;
      else
        SEQ_CONVERT = INT_2_NUMBER;
      ARR_CONVERT = -1;
    }
  else
    { if (sf->format == sam)
        SEQ_CONVERT = IUPAC_2_DNA;
      else
        SEQ_CONVERT = INT_2_IUPAC;
      ARR_CONVERT = '0';
    }

  return (status);
}


/*******************************************************************************************
 *
 *  RECORD PROCESSING
 *
 ********************************************************************************************/

// Looking for dq, dt, iq, mq, sq

 //  SAM Type letter to data size, save Z, H -> 9, and B -> 10

static int is_integer[128] =
  { 0, 0,  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0,
    0, 0,  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0,
    0, 0,  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0,
    0, 0,  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0,

    0, 0,  0, 1, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0, 0, 0,
    0, 0,  0, 1, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0,
    0, 0,  0, 1, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0, 0, 0,
    0, 0,  0, 1, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0,
  };

static int bam_tag_size[128] =
  { 0, 0,  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0,
    0, 0,  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0,
    0, 0,  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0,
    0, 0,  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0,

    0, 1, 10, 1, 0, 0, 0, 0,   9, 4, 0, 0, 0, 0, 0, 0,
    0, 0,  0, 2, 0, 0, 0, 0,   0, 0, 9, 0, 0, 0, 0, 0,
    0, 0,  0, 1, 8, 0, 4, 0,   0, 4, 0, 0, 0, 0, 0, 0,
    0, 0,  0, 2, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0,
  };

  // Convert big endian data to little endian

static void flip_auxilliary(uint8 *s, uint8 *e)
{ int  i, n, x;

  while (s < e)
    { s += 2;
      x = bam_tag_size[*s++];
      switch (x)
      { case 1:
          s += 1;
          break;
        case 2:
          flip_short(s);
          s += 2;
          break;
        case 4:
          if (is_integer[s[-1]])
            flip_int(s);
          s += 4;
          break;
        case 8:
          flip_double(s);
          s += 8;
          break;
        case 9:   //  Z or H
          while (*s != 0)
            s += 1;
          s += 1;
          break;
        case 10:  //  B
          x = bam_tag_size[*s++];
          flip_int(s);
          n = *((int *) s);
          s += 4;
          switch (x)
          { case 1:
              s += n;
              break;
            case 2:
              for (i = 0; i < n; i++, s += 2)
                flip_short(s);
              break;
            case 4:
              for (i = 0; i < n; i++, s += 4)
                flip_int(s);
              break;
            case 8:
              for (i = 0; i < n; i++, s += 8)
                flip_double(s);
              break;
          }
          break;
       }
    } 
}

static int bam_record_read(samFile *sf)
{ int ldata, lname, lcigar, lseq, aux;

  { int    ret;      //  read next block
    uint32 x[9];

    if ((ret = gzread(sf->ptr, x, 36)) != 36)
      { if (ret == 0)
          return (0);   // normal end-of-file
        else
          { fprintf(stderr,"%s: Unexpected end of input file\n",Prog_Name);
            return (-1);
          }
      }

    if (sf->is_big)
      { flip_int(x);
        flip_int(x + 3);
        flip_int(x + 4);
        flip_int(x + 5);
      }

    ldata  = ((int) x[0]) - 32;
    lname  = (x[3] & 0xff);
    lcigar = (x[4] & 0xffff);
    lseq   = x[5];

    if (ldata < 0 || lseq < 0 || lname < 1)
      return (-1);

    aux = lname + (lcigar<<2) + ((lseq + 1)>>1) + lseq;
    if (aux > ldata)
      { fprintf(stderr,"%s: Non-sensical SAM record, file corrupted?\n",Prog_Name);
        return (-1);
      }

    if (make_room(ldata))
      return (-1);

    if (gzread(sf->ptr, data, ldata) != ldata)
      { fprintf(stderr,"%s: Unexpected end of input file\n",Prog_Name);
        return (-1);
      }

    if (sf->is_big)
      flip_auxilliary(data+aux, data+ldata);
  }
  
  { uint8 *t;     //  Load header and sequence from required fields
    int    i, e;
    char  *seq, *eoh;

    if (lseq <= 0)
      { fprintf(stderr,"%s: no sequence for subread !?\n",Prog_Name);
        return (-1);
      }

    if (init_record(lseq))
      return (-1);

    theR.header = (char *) data;
    theR.len    = lseq;
    seq         = theR.seq;

    eoh = index(theR.header,'/');
    if (eoh != NULL)
      *eoh = 0;

    t = data + (lname + (lcigar<<2)); 
    for (e = i = 0; i < lseq-1; i += 2, e++)
      { seq[i]   = SEQ_CONVERT[t[e] >> 4];
        seq[i+1] = SEQ_CONVERT[t[e] & 0xf];
      }
    if (i < lseq)
      seq[i] = SEQ_CONVERT[t[e] >> 4];
  }

  return (1);
}

#define CHECK(cond, msg)				\
{ if ((cond))						\
    { fprintf(stderr, "%s: %s\n", Prog_Name, msg);	\
       return (-1);					\
    }							\
}

#define NEXT_ITEM(b,e)					\
{ b = e;						\
  while (*e != '\0' && *e != '\t')			\
    e++;						\
  CHECK( *e == '\0', "Missing one or more fields")	\
  *e = 0;						\
}

static int sam_record_read(samFile *sf)
{ char  *p;
  int    qlen, ret;

  //  read next line

  ret = sam_getline(sf,0);
  if (ret <= 0)
    return (ret);

  p = (char *) data;

  { char *q, *seq;     //  Load header and sequence from required fields
    int   i;

    NEXT_ITEM(q,p)
    qlen = p-q;
    CHECK( qlen <= 1, "Empty header name")
    CHECK( qlen > 255, "Header is too long")

    theR.header = q;
    q = index(q,'/');         // Truncate pacbio well & pulse numbers
    if (q != NULL && q < p)
      *q = 0;

    for (i = 0; i < 8; i++)   // Skip next 8 required fields
      { p = index(p+1,'\t');
        CHECK( p == NULL, "Too few required fields in SAM record, file corrupted?")
      }
    p += 1;

    NEXT_ITEM(q,p)
    qlen = p-q;
    CHECK (*q == '*', "No sequence for read?");

    if (init_record(qlen))
      return (-1);

    seq = theR.seq;
    theR.len = qlen;
    for (i = 0; i < qlen; i++)
      seq[i] = SEQ_CONVERT[(int) (*q++)];

    p = index(p+1,'\t');  // Skip qual
    CHECK( p == NULL, "No auxilliary tags in SAM record, file corrupted?")
  }

  return (1);
}

static samRecord _SAM_EOF;
samRecord *SAM_EOF = &_SAM_EOF;

samRecord *sam_record_extract(samFile *sf)
{ int64 ret;

  if (sf->format == bam)
    ret = bam_record_read(sf);
  else
    ret = sam_record_read(sf);

  if (ret <= 0)
    return (SAM_EOF);

  return (&theR);
}
