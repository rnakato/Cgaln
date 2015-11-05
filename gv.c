/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of Cgaln-1.*.* sources.
 */
#include "seq.h"

int seedfigure[][SEEDWEIGHT_MAX+1] = {
  {},                                /* 0 */
  {0},                               /* 1 */
  {0,2},                             /* 2 */
  {0,1,3},                           /* 3 */
  {0,1,2,4},                         /* 4 */
  {0,1,2,4,6},                       /* 5 */
  {0,1,4,5,6,8},                     /* 6 */
  {0,1,3,4,5,7,9},                   /* 7 */
  {0,1,3,6,7,9,10,11},               /* 8 */
  {0,1,2,4,5,8,10,12,13},            /* 9 */
  {0,1,3,4,7,9,10,12,13,14},         /* 10 */
  {0,1,2,4,7,9,12,13,15,16,17},      /* 11 */
  {0,1,2,3,5,7,10,11,14,16,17,18},   /* 12 */
  {0,1,2,3,5,8,9,12,13,15,17,18,19}  /* 13 */
};
int window[] = {0,1,3,4,5,7,9,10,12,14,15,18,19,20};  /* length of spaced seed */
int size[] = {1,4,16,64,256,1024,4096,16384,65536,262144,1048576,4194304,16777216,67108864};   /* 4^n */
int reverse, block, iterative, ext_gappeddp, short_reverse, printregion, filtercln, consistent;

option opt;
param par;
inputdata *idata_a, *idata_b;
outputalignment aln_for, aln_rev;

TYPE_SEEDT *seedtable_a, *seedtable_b, *seedtable_a_revcom;
TYPE_BLKT **blktable_a, **blktable_a_revcom;
struct poistable *poistable_a, *poistable_b;
