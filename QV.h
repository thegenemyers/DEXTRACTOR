/********************************************************************************************
 *
 *  Recreate all the .fasta files that have been loaded into a specified database.
 *
 *  Author:  Gene Myers
 *  Date  :  May 2014
 *
 ********************************************************************************************/

#ifndef _QV_COMPRESSOR

#define _QV_COMPRESSOR

#include "DB.h"

  //  A PacBio compression scheme

typedef struct
  { void    *delScheme;   //  Huffman scheme for deletion QVs
    void    *insScheme;   //  Huffman scheme for insertion QVs
    void    *mrgScheme;   //  Huffman scheme for merge QVs
    void    *subScheme;   //  Huffman scheme for substitution QVs
    void    *dRunScheme;  //  Huffman scheme for deletion run lengths (if delChar > 0)
    void    *sRunScheme;  //  Huffman scheme for substitution run lengths (if subChar > 0)
    int      delChar;     //  If > 0, run-encoded deletion value
    int      subChar;     //  If > 0, run-encoded substitution value
    int      flip;        //  Need to flip multi-byte integers
    char    *prefix;      //  Header line prefix
  } QVcoding;

  // Read the next nlines of input, and QVentry returns a pointer to the first line if needed.

int       Read_Lines(FILE *input, int nlines);
char     *QVentry();

  // Read the .quiva file on input and record frequency statistics.  If zero is set then
  //   restart statistics gathering, otherwise accumulate statisttics.

uint64    QVcoding_Scan(FILE *input, int zero);

  // Given QVcoding_Scan has been called at least once, create an encoding scheme based on
  //   the accumulated statistics and return a pointer to it.  The returned encoding object
  //   is *statically allocated within the routine.  If lossy is set then use a lossy scaling
  //   for the insertion and merge streams.

QVcoding *Create_QVcoding(int lossy);

  //  Read/write a coding scheme to input/output.  The encoding object returned by the reader
  //    is *statically* allocated within the routine.

QVcoding *Read_QVcoding(FILE *input);
void      Write_QVcoding(FILE *output, QVcoding *coding);

  //  Free all the auxiliary storage associated with coding (but not the object itself!)

void      Free_QVcoding(QVcoding *coding);

  //  Assuming the file pointer is positioned just beyond an entry header line, read the
  //    next set of 5 QV lines, compress them according to 'coding', and output.  If lossy
  //    is set then the scheme is a lossy one.

void      Compress_Next_QVentry(FILE *input, FILE *output, QVcoding *coding, int lossy);

  //  Assuming the input is position just beyond the compressed encoding of an entry header,
  //    read the set of compressed encodings for the ensuing 5 QV vectors, decompress them,
  //    and write their decompressed values to output.  The parameter rlen computed from the
  //    preceeding header line, critically provides the length of each of the 5 vectors.

void      Uncompress_Next_QVentry(FILE *input, FILE *output, QVcoding *coding, int rlen);

#endif // _QV_COMPRESSOR
