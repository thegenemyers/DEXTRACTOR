PATH_HDF5 = /sw/apps/hdf5/current
PATH_HDF5 = /usr/local/hdf5
CFLAGS = -O3 -Wall -Wextra -fno-strict-aliasing

all: dextract dexta undexta dexqv undexqv

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
	rm -f dextract dexta undexta dexqv undexqv dextract.tar.gz
	rm -fr *.dSYM
	rm -f dextract.tar.gz

install:
	cp dextract dexta undexta dexqv undexqv ~/bin

package:
	make clean
	tar -zcf dextract.tar.gz README Makefile *.c *.h
