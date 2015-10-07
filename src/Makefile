DEFS = -DPKGDATADIR=\"$(R_PACKAGE_DIR)/extdata/\"

CFILES = \
	SeqLibR.c  \
	sequence.c  \
	energy.c  \
	util.c  \
	getopt.c  \
	xmalloc.c  \
	hybrid-ss-min.c

OFILES = $(CFILES:.c=.o)

all: SeqLibR.so

.c.o:
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) $(DEFS) -c $< -o $@


# target and dependecies
SeqLibR.so : $(OFILES) Makefile
	$(CC) -shared $(CFLAGS) -o SeqLibR.so $(OFILES) $(LDFLAGS) $(ALLLIBS)  

clean :
	rm -rf *.o core SeqLibR.so *~

SeqLibR.so : $(OFILES)
SeqLibR.o : SeqLibR.c sequence.h
sequence.o : sequence.c sequence.h hybrid-ss-min.h
energy.o : energy.c energy.h xmalloc.h util.h
getopt.o : getopt.h
xmalloc.o: xmalloc.c xmalloc.h
hybrid-ss-min.o: util.h hybrid-ss-min.c

