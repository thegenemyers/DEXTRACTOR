PATH_HDF5 = /sw/apps/hdf5/current
PATH_HDF5 = /usr/local/hdf5

DEST_DIR  = ~/bin

CFLAGS = -O3 -Wall -Wextra -Wno-unused-result -fno-strict-aliasing

ALL = dextract dexta undexta dexqv undexqv

all: $(ALL)

dextract: dextract.c DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -I$(PATH_HDF5)/include -L$(PATH_HDF5)/lib -o dextract dextract.c DB.c QV.c -lhdf5

dexta: dexta.c DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o dexta dexta.c DB.c QV.c

undexta: undexta.c DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o undexta undexta.c DB.c QV.c

dexqv: dexqv.c DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o dexqv dexqv.c DB.c QV.c

undexqv: undexqv.c DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o undexqv undexqv.c DB.c QV.c

clean:
	rm -f $(ALL)
	rm -fr *.dSYM
	rm -f dextract.tar.gz

install:
	cp $(ALL) $(DEST_DIR)

package:
	make clean
	tar -zcf dextract.tar.gz README Makefile *.h *.c
