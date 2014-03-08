CFLAGS = -O4 -Wall -Wextra

all: dextract dexta undexta dexqv undexqv

dextract: dextract.c shared.c
	gcc $(CFLAGS) -o dextract dextract.c -lhdf5

dexta: dexta.c shared.c
	gcc $(CFLAGS) -o dexta dexta.c

undexta: undexta.c shared.c
	gcc $(CFLAGS) -o undexta undexta.c

dexqv: dexqv.c shared.c
	gcc $(CFLAGS) -o dexqv dexqv.c

undexqv: undexqv.c shared.c
	gcc $(CFLAGS) -fno-strict-aliasing -o undexqv undexqv.c

clean:
	rm -f dextract dexta undexta dexqv undexqv dextract.tar.gz

install:
	cp dextract dexta undexta dexqv undexqv ~/bin

package:
	tar -zcf dextract.tar.gz README *.c Makefile
