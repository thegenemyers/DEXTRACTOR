CFLAGS = -O3 -Wall -Wextra -fno-strict-aliasing
% CFLAGS = -O3 -Wall -Wextra -fno-strict-aliasing -L/sw/apps/hdf5/current/lib -I/sw/apps/hdf5/current/include

all: dextract dexta undexta dexqv undexqv

dextract: dextract.c DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o dextract dextract.c DB.c QV.c -lhdf5

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

install:
	cp dextract dexta undexta dexqv undexqv ~/bin

package:
	tar -zcf dextract.tar.gz README Makefile *.c *.h
