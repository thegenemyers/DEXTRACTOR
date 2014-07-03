CC ?= gcc

CFLAGS ?= -O3 -march=native #Default flags in Release mode
CFLAGS +=  -Wall -Wextra 

LDFLAGS += -lhdf5 
AUTOTARGETS = dextract dexta undexta dexqv 
TARGETS = $(AUTOTARGETS) undexqv 

all: $(TARGETS)

$(AUTOTARGETS): % : %.c
	$(CC) $(CFLAGS) $(LDFLAGS) $^ -o $@ 

undexqv: undexqv.c 
	$(CC) $(CFLAGS) -fno-strict-aliasing -o $@ $<

clean:
	rm -f $(TARGETS) dextract.tar.gz

install:
	cp $(TARGETS) ~/bin

package:
	tar -zcf dextract.tar.gz README *.c Makefile
