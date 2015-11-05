/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of Cgaln-1.*.* sources.
 */
#include "util.h"

int SeedScoring_masked(Fasta *fst, int posi, int seedweight){
  int nuc, seedvalue = 0;
  int *sp = seedfigure[seedweight];
  while(seedweight--){
    nuc = fst->body[posi + *sp++];
    if(nuc <= -1) return(-1);
    if(nuc == 4) return(-2);
    seedvalue = (seedvalue << 2) + nuc;
  }
  return seedvalue;
}

int SeedScoring_notmasked(Fasta *fst, int posi, int seedweight){
  int nuc, seedvalue = 0;
  int *sp = seedfigure[seedweight];
  while(seedweight--){
    nuc = fst->body[posi + *sp++];
    if(nuc == -1) return(-1);
    if(nuc == 4) return(-2);
    seedvalue = (seedvalue << 2) + define_base_nomask(nuc);
  }
  return seedvalue;
}

char define_base_nomask(char value){
  char x=-1;
  if(value<-1){
    switch(value){
    case -2: x=0; break;
    case -3: x=1; break;
    case -4: x=2; break;
    case -5: x=3; break;
    }
  }else x=value;
  return x;
}

char convert_complement(char value){
  char x;
  switch(value){
  case 3: x=0; break;
  case 2: x=1; break;
  case 1: x=2; break;
  case 0: x=3; break;
  case -2: x=-5; break;
  case -3: x=-4; break;
  case -4: x=-3; break;
  case -5: x=-2; break;
  default: x=value; break;
  }
  return x;
}

char convert_num2base(char value){
  char c=0;
  switch(value){
  case  0: c='A'; break;
  case  1: c='C'; break;
  case  2: c='G'; break;
  case  3: c='T'; break;
  case -2: c='a'; break;
  case -3: c='c'; break;
  case -4: c='g'; break;
  case -5: c='t'; break;
  case -1: c='N'; break;
  case  4: c='$'; break;
  }
  return c;
}

char *checkfilename(char *inputname){
  char *p = strrchr(inputname, '/');
  if(p) return p+1;
  else return inputname;
}

FILE *my_fopen_r(char *filename){
  FILE *IN;
  if((IN = fopen(filename, "r"))==NULL){
    fprintf(stderr,"[E] Cannot open <%s>.\n", filename); 
    exit(1);
  }
  return IN;
}

FILE *my_fopen_w(char *filename){
  FILE *IN;
  if((IN = fopen(filename, "w"))==NULL){
    fprintf(stderr,"[E] Cannot open <%s>.\n", filename); 
    exit(1);
  }
  return IN;
}

FILE *my_fopen_ab(char *filename){
  FILE *IN;
  if((IN = fopen(filename, "ab"))==NULL){
    fprintf(stderr,"[E] Cannot open <%s>.\n", filename); 
    exit(1);
  }
  return IN;
}

void *my_malloc(size_t s, char *name){
  void *p;
  p = malloc(s);
  if(!p){
    fprintf(stderr, "[E] failed malloc: %s\n", name);
    exit(1);
  }
  return p;
}

void *my_calloc(size_t n, size_t s, char *name){
  void *p;
  p = calloc(n,s);
  if(!p){
    fprintf(stderr,"[E]failed calloc: %s\n", name); 
    exit(1);
  }
  return p;
}

void *my_realloc(void *p, size_t s, char *name){
  p = realloc(p,s);
  if(!p){
    fprintf(stderr,"[E]failed realloc: %s\n", name); 
    exit(1);
  }
  return p;
}
