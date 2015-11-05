/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of Cgaln-1.*.* sources.
 */
#include "readfasta.h"

void read_multifasta(FILE *fp, Fasta *fst, int strand, int *state){
  int i, numh=0, numb=0;
  char c;
  char *bf_head, *bf_body;
  int bodysize = SIZE_BODY_DEFAULT;
  bf_head = (char *)my_malloc(sizeof(char)*SIZE_HEAD, "bf_head");
  bf_body = (char *)my_malloc(sizeof(char)*bodysize, "bf_body");

  while((c = fgetc(fp)) != EOF){
    switch(*state){
    case 0:
      if(c=='>') *state=1;
      break;
    case 1:              /*header*/
      if(c=='\n'){
	*state=2;
	bf_head[numh] = '\0';
      }else{
	bf_head[numh] = c; 	
	numh++;
      }
      break;
    case 2:              /*body*/
      if(c=='>'){
	bf_body[numb] = '\0'; 
	*state=1;
	goto final;
      }else if(isalpha(c)){
	switch(c){
	case 'A': if(!strand) bf_body[numb] =  0; else bf_body[numb] =  3; break;
	case 'C': if(!strand) bf_body[numb] =  1; else bf_body[numb] =  2; break;
	case 'G': if(!strand) bf_body[numb] =  2; else bf_body[numb] =  1; break;
	case 'T': if(!strand) bf_body[numb] =  3; else bf_body[numb] =  0; break;
	case 'a': if(!strand) bf_body[numb] = -2; else bf_body[numb] = -5; break;
	case 'c': if(!strand) bf_body[numb] = -3; else bf_body[numb] = -4; break;
	case 'g': if(!strand) bf_body[numb] = -4; else bf_body[numb] = -3; break;
	case 't': if(!strand) bf_body[numb] = -5; else bf_body[numb] = -2; break;
	case 'N': bf_body[numb] = -1; break;
	case '$': bf_body[numb] =  4; break;
	default:  bf_body[numb] = -1; break;
	}
	numb++;
	if(numb >= bodysize){
	  bodysize += SIZE_BODY_DEFAULT;
	  bf_body = my_realloc(bf_body, sizeof(char)*bodysize, "bf_body");
	}
      }else{
      }
      break;
    }
  }
  
  final:
  /*----- substitute for fst -----*/
  fst->length = numb;
  fst->head = (char *)my_malloc(sizeof(char)*(numh+1), "fst->head");
  strcpy(fst->head, bf_head);
  fst->body = (char *)my_malloc(sizeof(char)*(numb+1), "fst->body");

  for(i=0; i<numb; i++) fst->body[i] = bf_body[i + strand*(numb-1-2*i)];  /* if(strand==1) fst->body[i] = bf_body[numb-1-i] */
  
  free(bf_head);
  free(bf_body);
  return;
}


int countfasta(char *seqname){
  int i, len, count=0, state=2;
  char c, buf[BUF];
  FILE *IN = my_fopen_r(seqname);

  while((len=fread(buf,sizeof(char),BUF,IN))>0){
    for(i=0; i<len; i++){
      c=buf[i];
      switch(state){
      case 1:   if(c=='\n') state = 2;
	break;
      case 2:   if(c=='>'){
	  state=1;
	  count++;
	}
	break;
      }
    }
  }

  fclose(IN);
  return count;
}
