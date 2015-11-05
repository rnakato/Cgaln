/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of Cgaln-1.*.* sources.
 */
#ifndef READFASTA_H_
#define READFASTA_H_

#include <ctype.h>
#include "util.h"

#define SIZE_HEAD 1000
#define SIZE_BODY_DEFAULT 3000000
#define BUF 65536

/* readfasta.c*/
void read_multifasta(FILE *, Fasta *, int, int *);
int countfasta(char *);

#endif /*READFASTA_H_*/
