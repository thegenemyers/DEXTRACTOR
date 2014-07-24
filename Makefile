CFLAGS = -O4 -Wall -Wextra

all: dextract dexta undexta dexqv undexqv

dextract: dextract.c DB.c DB.h
	gcc $(CFLAGS) -o dextract dextract.c DB.c -lhdf5

dexta: dexta.c DB.c DB.h
	gcc $(CFLAGS) -o dexta dexta.c DB.c

undexta: undexta.c DB.c DB.h
	gcc $(CFLAGS) -o undexta undexta.c DB.c

dexqv: dexqv.c DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o dexqv dexqv.c DB.c QV.c

undexqv: undexqv.c DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -fno-strict-aliasing -o undexqv undexqv.c DB.c QV.c

clean:
	rm -f dextract dexta undexta dexqv undexqv dextract.tar.gz

install:
	cp dextract dexta undexta dexqv undexqv ~/bin

package:
	tar -zcf dextract.tar.gz README Makefile *.c *.h
