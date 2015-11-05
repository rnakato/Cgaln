/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of Cgaln-1.*.* sources.
 */
#ifndef NA_H_
#define NA_H_

#include <stdio.h>
#include <stdlib.h>
#include "cgaln.h"

/* aln_base.c */
#define HSPCONS 2
#define X_EXTENDHSP 50
/* iterativealign.c */
#define SEEDWEIGHT_2 7
#define HSP_SIZEMAX 1000
#define HSP_CHAIN 20
#define PENALTY_WEIGHT 2000
#define PENALTY_WEIGHT2 2
#define HSP_UNITE_THRE 10
#define IA_THRE1 50 
#define IA_THRE2 200

struct diag_table{
  int len;
  int x_e, sy;
  int S;
  int diag;
  int m_end;
};

struct hsplist{
  struct hsp *hsp;
  int num;
  int nummax;
  int last_i;
  int seedweight;
  int penaweight;
};

enum{FIRST, SECOND};

/* aln_base.c */
void NA(FILE *, Fasta *, Fasta *, outputalignment *, int, int *);
void renewHSP(struct hsplist *, struct diag_table *, Fasta *, Fasta *, int, int, int, int, int);
void hspcopy_from_diagtable(struct hsplist *, struct diag_table *, int);
void scan_diagtable(struct hsplist *, struct diag_table *, int, int);
struct hsplist *hsplist_new(int, int);
struct diag_table *diag_table_new(int, int);
void subst_region(struct region *, int, int, int, int);
/* iterativealign.c */
void chain_and_IA(struct hsplist *, struct region *, Fasta *, Fasta *, int, int, int);
int ChainHSP(struct hsplist *);
void maketable4region(Fasta *, int **, int **, int, int, int, int);
int DP4smallregion(struct region *, Fasta *, Fasta *, int, int);
/* dp.c */
void DP(struct region *, Fasta *, Fasta *, int);
void DP_xdrop(struct region *, Fasta *, Fasta *, int, int);
#endif /*NA_H_*/
