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

#ifndef _SAM_BAM
#define _SAM_BAM

#include <stdint.h>
#include <zlib.h>

#include "DB.h"

typedef enum { sam, bam } samFormat;

typedef struct
  { samFormat format;  //  sam or bam
    int       is_big;  //  endian (bam only)
    int       nline;   //  current line number (sam only)
    char     *name;    //  file name
    gzFile    ptr;     //  pointer to file descriptor
} samFile;

typedef struct
  { int   len;
    char *header;
    char *seq;
  } samRecord;

  // sam_open: NULL => error, open file otherwise
  // sam_close: 1 => error, 0 otherwise OK
  // sam_eof: 1 => eof or error, 0 otherwise 
  //   error message *will not* have been sent to stderr.

samFile *sam_open(char *sf);          //   Open a SAM/BAM file for reading
int      sam_close(samFile *sf);      //   Close an open SAM/BAM file
int      sam_eof(samFile *sf);        //   Return non-zero if at eof

  // sam_header_process: 0 => no codec, 1 => codec, -1 => error
  // sam_record_extract: NULL => error, EOF => error, filled in sam record otherwise
  //   error message *will* have been sent to stderr.

extern samRecord *SAM_EOF;

int        sam_header_process(samFile *sf, int numeric);
samRecord *sam_record_extract(samFile *sf)

#endif // _SAM_BAM
