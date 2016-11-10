PATH_HDF5 = /sw/apps/hdf5/current
PATH_HDF5 = /usr/local/hdf5

DEST_DIR  = ~/bin

CFLAGS = -O3 -Wall -Wextra -Wno-unused-result -fno-strict-aliasing

ALL = dextract dexta undexta dexar undexar dexqv undexqv dex2DB

all: $(ALL)

dextract: dextract.c sam.c bax.c expr.c expr.h bax.h DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -I$(PATH_HDF5)/include -L$(PATH_HDF5)/lib -o dextract dextract.c sam.c bax.c expr.c DB.c QV.c -lhdf5 -lz

dexta: dexta.c DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o dexta dexta.c DB.c QV.c

undexta: undexta.c DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o undexta undexta.c DB.c QV.c

dexar: dexar.c DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o dexar dexar.c DB.c QV.c

undexar: undexar.c DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o undexar undexar.c DB.c QV.c

dexqv: dexqv.c DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o dexqv dexqv.c DB.c QV.c

undexqv: undexqv.c DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o undexqv undexqv.c DB.c QV.c

dex2DB: dex2DB.c sam.c bax.c DB.c QV.c bax.h DB.h QV.h
	gcc $(CFLAGS) -I$(PATH_HDF5)/include -L$(PATH_HDF5)/lib -o dex2DB dex2DB.c sam.c bax.c DB.c QV.c -lhdf5 -lz

clean:
	rm -f $(ALL)
	rm -fr *.dSYM
	rm -f dextract.tar.gz

install:
	cp $(ALL) $(DEST_DIR)

package:
	make clean
	tar -zcf dextract.tar.gz README.md Makefile *.h *.c
