/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of Cgaln-1.*.* sources.
 */
#ifndef BA_H_
#define BA_H_
#include "cgaln.h"

/* aln_block.c */
#define ALLOW_RAND_CELL 5
#define SIZE_BL_DEFAULT 1000
#define CLNMAX_DEFAULT 10000
#define CONS_FC 1

struct colony{
  TYPE_BLKT ini_i,ini_j;
  double m_value; 
  TYPE_BLKT m_i,m_j;
  TYPE_BLKT_SIGHNED i_j_plus;
  TYPE_BLKT_SIGHNED i_j_minus; 
  struct clnarray *clnarray;
  int num;
};

struct colonyset{
  struct colony *cln;
  int num;
  int nummax;
};

struct column{
  double score;
  int id;
};


/* aln_block.c */
void BA(outputalignment *, int **, unsigned short **, int);
void Score_Blocksim(int *, unsigned short *, double *, TYPE_SEEDT *, TYPE_BLKT **, int ,int);
/* aln_colony.c */
void Find_colony(struct colonyset *, double *, struct column *, int);
void Align_colony(struct colony *, int **, unsigned short **, TYPE_SEEDT *, TYPE_BLKT **);

#endif /*BA_H_*/
