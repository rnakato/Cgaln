/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of Cgaln-1.*.* sources.
 */
#ifndef CGALN_H_
#define CGALN_H_

#include "readfasta.h"
#include "util.h"

#define SEED_RATIO_CONS 20
#define XDROP_DEFAULT 4000      /* default XDROP*/
#define X_D_RATIO_DEFAULT 1500  /* default ration between XDROP and gap penalty */ 

extern option opt;
extern int reverse, block, iterative, ext_gappeddp, short_reverse, printregion, filtercln, consistent;

extern TYPE_SEEDT *seedtable_a, *seedtable_b, *seedtable_a_revcom;
extern TYPE_BLKT **blktable_a, **blktable_a_revcom;
extern struct poistable *poistable_a, *poistable_b;
extern outputalignment aln_for, aln_rev;

/* readtable.c */
void read_table();
void read_table_revcom();
void get_prefix(char *, char *);
void read_summary(inputdata *);

/* output.c */
void output_BAresult();
void print_HSP(struct hsp *, Fasta *, Fasta *, int, int, int);
void print_DPresult(int, int, int, Fasta *, Fasta *, int);
void print_region(struct region *, int);
void output_fastaboundary();
#endif /*CGALN_H_*/
