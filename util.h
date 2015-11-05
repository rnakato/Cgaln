/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of Cgaln-1.*.* sources.
 */
#ifndef UTIL_H_
#define UTIL_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "seq.h"

/* util.c */
int SeedScoring_masked(Fasta *, int, int);
int SeedScoring_notmasked(Fasta *, int, int);
char convert_num2base(char);
char define_base_nomask(char);
char convert_complement(char);
char *checkfilename(char *);
FILE *my_fopen_r(char *);
FILE *my_fopen_w(char *);
FILE *my_fopen_ab(char *);
void *my_malloc(size_t, char *);
void *my_calloc(size_t, size_t, char *);
void *my_realloc(void *, size_t, char *);

#endif /*UTIL_H_*/
