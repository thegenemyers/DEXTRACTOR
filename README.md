# The Dextractor and Compression Command Library

## _Authors: Gene Myers_

For typeset documentation, examples of use, and design philosophy please go to
my [blog](https://dazzlerblog.wordpress.com/command-guides/dextractor-command-guide).

The Dextractor commands allow one to pull exactly and only the information needed for
assembly and reconstruction from the source HDF5 files produced by the PacBio
RS II sequencer, or from the source BAM files produced by the PacBio Sequel
sequencer.  The program **dextract** can pull information from either HDF5 .bax.h5 files
or BAM/SAM .subreads.[bs]am files, where it determines the file type based on the file
suffix.  Both file
formats contain the basic sequence information for each subread and contain either the
quality streams needed by the Quiver consensus programs, or the pulse width and signal-
to-noise ratio information needed by the newer Arrow consensus program, or both. 
Dextract can be directed to pull the sequence information into a .fasta file by
setting the -f option, to pull the Quiver information (if present) into a .quiva file
by setting the -q option, and to pull the Arrow information (if present) into a .arrow
file by setting the -a option.  As detailed in the command description below, the
program can produce any combination of the three output files (if the information is
present).

For each of the three extracted file types -- .fasta, .quiva, and .arrow --
the library contains commands to compress the given file type, and to decompress
it, which is a reversible process delivering the original uncompressed file.   The
compressed .fasta files, with the extension .dexta, consume 1/4 byte per base.
The compressed .quiva files, with the extension .dexqv, consume 1.5 bytes per base on
average, and the compressed .arrow files, with the extension .dexar, consume 1/4 byte
per base.  In this way, users can keep just the data needed for assembly spooled up on disk
in much less space than the source files and keep said source files archived to a cheap
backup medium such as tape, should the raw data ever need to be consulted again
(we expect never unless the spooled up data is compromised or lost in some way).
The compressor/decompressor pairs are endian-aware so moving compressed files between
machines is possible.

Finally, one can also directly produce a Dazzler database from the source files without
explicitly extracting the uncompressed file types.  Since a Dazzler database keeps the
extracted information in compressed form, it is actually recommended that you directly
build a DB from your source data and then archive the raw data.  Doing so saves time
and there is no need for any large intermediate .fasta, .quiva, or .arrow files.  But
we continue to provide dextract and sequel for those who which to build their own
assembly pipelines and not use our DBs as an organizing principle.

```
1. dextract [-vfaq] [-o[<path>]] [-e<expr(ln>=500 && rq>=750)>] <input:pacbio> ...
```

Dextract takes a series of .bax.h5 or .subreads.[bs]am files as input, and depending on
the option flags settings produces:

1. (-f) a .fasta file containing subread sequences, each with a "standard" Pacbio header
consisting of the movie name, well number, pulse range, and read quality value.

2. (-a) a FASTA format .arrow file containing the pulse width stream for each subread, with
a header that contains the movie name and the 4 channel SNR values.

3. (-q) a FASTQ-like .quiva file containing for each subread the same header as the
.fasta file above, save that it starts with an @-sign, followed by the 5 quality
value streams used by Quiver, one per line, where the order of the streams is:
deletion QVs, deletion Tags, insertion QVs, merge QVs, and last substitution QVs. 

If the -v option is set then the program reports the processing of each PacBio input
file, otherwise it runs silently.  If none of the -f, -a, or -q flags is set, then by
default -f is assumed.  The destination of the extracted information is controlled
by the -o parameter as follows:

1. If -o is absent, then for each input file X.bax.h5 or X.subreads.[bs]am, dextract
will produce X.fasta, X.arrow, and/or X.quiva as per the option flags.

2. If -o is present and followed by a path Y, then the concatenation of the output for
the input files is placed in Y.fasta, Y.arrow, and/or Y.quiva as per the option flags.

3. If -o is present but with no following path, then the output is sent to the standard
output (to enable a UNIX pipe if desired).  In this case only one of the flags -f, -a,
or -q can be set.

One can select which subreads are transferred to the output with the -e option that
contains a C-style boolean expression over integer constants and the 8 variable
names that for a given subread designate:

```
zm  - well number
ln  - length of the subread
rq  - quality value of the subread (normalized to [0,1000])
bc1 - # of the first barcode
bc2 - # of the second barcode
bq  - quality of the barcode detection (normalized to [0,100])
np  - number of passes producing subread
qs  - start pulse of subread
```

For each subread, the expression is evaluated and the subread is output only if the
expression is true.
The default filter is "ln >= 500 && rq >= 750", that is, only subreads longer than
500bp with a quality score of .75 or better will be output.  If a variable is undefined
for a subread (e.g. bar codes are often not present), the value of the variable will be -1.

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

With the -i option set, the program runs as a UNIX pipe that takes .fasta (.dexta)
input from the standard input and writes .dexta (.fasta) to the standard output.
In this case the -k option has no effect.

```
3. dexar   [-vk] ( -i | <path:arrow> .. .)
   undexar [-vk] [-w<int(80)>] ( -i | <path:dexar> ... )
```

Dexar compresses a set of .arrow files
and replaces them with new files with a .dexar extension.  That is,
submitting G.arrow will result in a compressed image G.dexar, and G.arrow
will no longer exist.  With the -k option the .arrow source is *not* removed.  If
-v is set, then the program reports its progress on each file.  Otherwise it runs
completely silently (good for batch jobs to an HPC cluster).  The compression
factor is always slightly better than 4.0.  Undexar reverses the compression of
dexar, replacing the uncompressed image of G.dexar with G.arrow.  By default the
sequences output by undexar are 80 chars per line.  The characters per line, or
line width, can be set to any positive value with the -w option

With the -i option set, the program runs as a UNIX pipe that takes .arrow (.dexar)
input from the standard input and writes .dexar (.arrow) to the standard output.
In this case the -k option has no effect.


```
4. dexqv   [-vkl] <path:quiva> ...
   undexqv [-vkU] <path:dexqv> ...
```

Dexqv compresses a set of .quiva files into new files with a
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

```
5. dex2DB [-vlaq] [-e<expr(ln>=500 && rq>=750)>] 
              <path:db> ( -f<file> | <input:pacbio> ... )
```

Builds an initial data base, or adds to an existing database, *directly* from either
(a) the list of .bax.h5 or .subreads.[bs]am files following the database name argument,
or (b) the list of PacBio source files in \<file\> if the -f option is used.
One can filter which reads are added to the DB with the -e option (see dextract above).

On a first call to dex2DB, i.e. one that creates the database, the settings of the
-a and -q flags, determine the type of the DB as follows.  If the -a option is set,
then Arrow information is added to the DB and the DB is an Arrow-DB (A-DB).  If
the -q option is set, then Quiver information is added to the DB and the DB is a
Quiver-DB (Q-DB).  If neither flag is set, then the DB is a simple Sequence-DB (S-DB).
One can never specify
both the -a and the -q flags together.
After the first additions to the database, the settings of the -a and -q flags must be
consistent with the type of the database call, or left unset in which case the type
of the database determines which information is pulled from the inputs.

If the files are being added to an existing database, and the
partition settings of the DB have already been set (see DBsplit below), then the
partitioning of the database is updated to include the new data.
If the -l option is set, then either the -q option must be set or the DB must already be
established as a Q-DB.  If it is set then the Quiver streams are lossy compressed
(see quiva2DB in the DAZZ_DB module).
