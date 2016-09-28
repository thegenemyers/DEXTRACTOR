# The Dextractor and Compression Command Library

## _Authors: Gene Myers_

For typeset documentation, examples of use, and design philosophy please go to
my [blog](https://dazzlerblog.wordpress.com/command-guides/dextractor-command-guide).

The Dextractor commands allow one to pull exactly and only the information needed for
assembly and reconstruction from the source .bax.h5 HDF5 files produced by the PacBio
RS II sequencer.  Generally speaking, this information is the sequence of all the reads
coded in the .bax.h5 file and a number of quality value (QV) streams needed by Quiver
to produce a highly accurate consensus sequence as the last step in the assembly
process.  The Dextractor therefore produces a .fasta file of the sequence of all the
reads, and a .quiva file containing the QV stream information in a .fastq readable
format.  For each of these two file types the library contains commands to compress the
given file type, and to decompress it, which is a reversible process delivering the
original uncompressed file.  In this way, users of a PacBio can keep the data needed
for assembly spooled up on disk in 1/14th the space occupied by the .bax.h5 files which
can be archived to a cheap backup medium such as tape, should the raw data ever need to
be consulted again (we expect never unless the spooled up data is compromised or lost
in some way).  The compressor/decompressor pairs are endian-aware so moving compressed
files between machines is possible.

```
1. dextract  [-vq] [-o[<path>]] [-l<int(500)>] [-s<int(750)>] <input:bax.h5> ...
```

The dextract'or takes the .bax.h5 files produced for a given SMRT cell as
input and:

1. if the -o option is set, then the information needed for Quiver is extracted
and put in a file named \<path\>.quiva.  If the -q option is not set, then the
sequence of each read is placed in a file named \<path\>.fasta, otherwise a
.fastq file of the sequence and the imputed "quality values" for each base in
the sequence is placed in a file named \<path\>.fastq.  We personally do not
find these values useful and so never set -q but we give you the option in
case your downstream processes use such values.

  If \<path\> is missing, then the path of the first .bax.h5 file is used for the
  output file name, less any suffixes which are replaced by .fasta and .quiva.
  E.G., the call "dextract -o EColi.1.bax.h5 EColi.2.bax.h5 Ecoli.3.bax.h5"
  will result in the files EColi.fasta and Ecoli.quiva.

2. if the -o option is not set, then if the -q option is also not set, then a
.fasta file of the sequence of each read is written to the standard output.
Otherwise a .fastq file is written to the standard output.

If the -v option is set then the program reports the processing of each .bax.h5
file, otherwise it runs silently.  The parameter -l determines the shortest read
length to be extracted (default 500) and the -s parameter determines the minimum
quality/score of reads to be extracted (default 750 = 75%).


```
2. dexta   [-vk] ( -i | <path:fasta> .. .)
   undexta [-vkU] [-w<int(80)>] ( -i | <path:dexta> ... )
```

Dexta compresses a set of .fasta files (produced by either Pacbio's software or
dextract) and replaces them with new files with a .dexta extension.  That is,
submitting G.fasta will result in a compressed image G.dexta, and G.fasta
will no longer exist.  With the -k option the .fasta source is *not* removed.  If
-v is set, then the program reports its progress on each file.  Otherwise it runs
completely silently (good for batch jobs to an HPC cluster).  The compression
factor is always slightly better than 4.0.  Undexta reverses the compression of
dexta, replacing the uncompressed image of G.dexta with G.fasta.  By default the
sequences output by undexta are in lower case and 80 chars per line.  The -U
option specifies upper case should be used, and the characters per line, or line
width, can be set to any positive value with the -w option

With the -i option set, the programs run as UNIX pipes that take .fasta (.dexta)
input from the standard input and write .dexta (.fasta) to the standard output.
In this case the -k option has no effect.


```
3. dexqv   [-vkl] <path:quiva> ...
   undexqv [-vkU] <path:dexqv> ...
```

Dexqv compresses a set of .quiva files (produced by dextract) into new files with a
.dexqv extension.  That is, submitting G.quiva will result in a compressed image
G.dexqv, and G.quiva will not longer exist.  The -k flag prevents the removal
of G.quiva.   With -v set progress is reported, otherwise the command runs
silently.  If slightly more compression is desired at the expense of being a bit
"lossy" then set the -l option.  This option is experimental in that it remains to
be seen if Quiver gives the same results with the scaled values responsible for the
loss.  Undexqv reverses the compression of dexqv, replacing the uncompressed image
of G.dexqv with G.quiva.  The flags are analgous to the v & k flags for dexqv.
The compression factor is typically 3.4 or so (4.0 or so with -l set).  By .fastq
convention each QV vector is output by undexqv as a line without intervening
new-lines, and by default the Deletion Tag vector is in lower case letters. The -U
option specifies upper case letters should instead be used for said vector.


To compile the programs you must have the HDF5 library installed on your system and
the library and include files for said must be on the appropriate search paths.  The
HDR5 library in turn depends on the presence of zlib, so make sure it is also installed
on your system.  The most recent version of the source for the HDF5 library can be
obtained [here](https://support.hdfgroup.org/downloads/index.html).
