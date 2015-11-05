CC = gcc
OBJS_UTIL = readfasta.o readtable.o util.o gv.o instance.o
OBJS = cgaln.o aln_block.o aln_colony.o aln_base.o iterativealign.o dp.o output.o $(OBJS_UTIL)
OBJS2 = maketable.o $(OBJS_UTIL)
HEADS = seq.h util.h readfasta.h cgaln.h BA.h NA.h dp.h
HEADS2 = seq.h util.h readfasta.h
CFLAGS += -Wall -O3 -W
CFLAGS2 = $(CFLAGS) -lm

all: Cgaln maketable

Cgaln: $(OBJS)
	$(CC) -o $@ $^ $(CFLAGS2) 

maketable: $(OBJS2)
	$(CC) -o $@ $^ $(CFLAGS2) 

.SUFFIXES: .o .c
.c.o:
	gcc -c $< $(CFLAGS)

clean:
	rm $(OBJS) maketable.o
	rm *~

$(OBJS): $(HEADS) Makefile
$(OBJS2): $(HEADS2) Makefile