CC ?= gcc
CFLAGS  += -O4 -Wall -Wextra 
LDFLAGS += -lhdf5 
AUTOTARGETS = dextract dexta undexta dexqv 

all: $(AUTOTARGETS) undexqv 

$(AUTOTARGETS): % : %.c
	$(CC) $(CFLAGS) $(LDFLAGS) $^ -o $@ 

undexqv: undexqv.c 
	$(CC) $(CFLAGS) -fno-strict-aliasing -o undexqv undexqv.c

clean:
	rm -f dextract dexta undexta dexqv undexqv dextract.tar.gz

install:
	cp dextract dexta undexta dexqv undexqv ~/bin

package:
	tar -zcf dextract.tar.gz README *.c Makefile
