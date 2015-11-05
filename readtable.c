/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of Cgaln-1.*.* sources.
 */
#include "seq.h"
#include "cgaln.h"
struct elem{
  char str[12800];
};

static void *read_table_each(size_t size, int tablesize, inputdata *idata, char *postfix){
  FILE *IN;
  char filename[INPUTFILENAME_MAX];
  int i, strand=0;
  if(strstr(postfix, "-revcom")) strand=1;
  unsigned long seeksize=0;

  if(!strcmp(postfix, ".poistable")){
    sprintf(filename, "%s/%s%s-%d", par.tabledir, checkfilename(idata->seqname), postfix, par.blocksize);
  }else if(!strcmp(postfix, ".seedtable") || !strcmp(postfix, ".seedtable-revcom")){
    sprintf(filename, "%s/%s%s-%d", par.tabledir, checkfilename(idata->seqname), postfix, par.kmer_ba);    
  }else{
    sprintf(filename, "%s/%s%s-%d-%d", par.tabledir, checkfilename(idata->seqname), postfix, par.kmer_ba, par.blocksize);
  }

  if(!strcmp(postfix, ".blktable") || !strcmp(postfix, ".blktable-revcom")){
    if(idata->cnt){
      for(i=0; i<idata->cnt; i++) seeksize += size*idata->ksum[strand][i];
    }
  }else seeksize = size*tablesize*idata->cnt;

  //  printf("%s, cnt %d, seeksize=%ld, tablesize=%d\n", filename, idata->cnt, seeksize, (int)(size*tablesize));

  void *table_p = my_malloc(size*tablesize, postfix);
  IN = my_fopen_r(filename);
  if(seeksize) fseek(IN, seeksize, SEEK_SET);
  if(fread(table_p, size*tablesize, 1, IN)!=1){
    if(ferror(IN)!=0) printf("ferror\n");
    if(feof(IN)!=0) printf("feof\n");
    fprintf(stderr,"[E]fread error:%s\n", filename);
    exit(1);
  }
  fclose(IN);
  return table_p;
}

static TYPE_BLKT **alloc_blktable(TYPE_BLKT *blktable_sum, TYPE_SEEDT *seedtable, int size){
  int i,j;
  TYPE_BLKT **blktable = (TYPE_BLKT **)my_malloc(size * sizeof(TYPE_BLKT *), "blktable");
  for(i=0; i<size; i++){
    if(seedtable[i]){  
      blktable[i] = (TYPE_BLKT *)my_malloc(seedtable[i] * sizeof(TYPE_BLKT), "blktable[i]");
      for(j=0; j<seedtable[i]; j++) blktable[i][j] = *(blktable_sum++);
    }else{
      blktable[i]=NULL;
    }
  }
  return blktable;
}

void read_table(){
  int i;
  /*----- seedtable_a -----*/
  seedtable_a = (TYPE_SEEDT *)read_table_each(sizeof(TYPE_SEEDT), size[par.kmer_ba], idata_a, ".seedtable");
  /*----- seedtable_b -----*/
  seedtable_b = (TYPE_SEEDT *)read_table_each(sizeof(TYPE_SEEDT), size[par.kmer_ba], idata_b, ".seedtable");
  /* calculate threshold of seed occurrence num */
  par.num_a=0; par.num_b=0;
  for(i=0; i<size[par.kmer_ba]; i++){
    if(seedtable_a[i]) par.num_a++;
    if(seedtable_b[i]) par.num_b++;
  }
  par.seedthre_a = idata_a->length_each[idata_a->cnt] / par.num_a * SEED_RATIO_CONS;
  par.seedthre_b = idata_b->length_each[idata_b->cnt] / par.num_b * SEED_RATIO_CONS;

  /*----- blktable_a -----*/
  TYPE_BLKT *blktable_sum = (TYPE_BLKT *)read_table_each(sizeof(TYPE_BLKT), idata_a->ksum[0][idata_a->cnt], idata_a, ".blktable");
  blktable_a = alloc_blktable(blktable_sum, seedtable_a, size[par.kmer_ba]);
  free(blktable_sum);

  /*----- poistable_a -----*/
  poistable_a = (struct poistable *)read_table_each(sizeof(struct poistable), KMER_NUMMAX*POISSON_MAX, idata_a, ".poistable");
  /*----- poistable_b -----*/
  poistable_b = (struct poistable *)read_table_each(sizeof(struct poistable), KMER_NUMMAX*POISSON_MAX, idata_b, ".poistable");

  return;
}


void read_table_revcom(){
  /*----- seedtable_a_revcom -----*/
  seedtable_a_revcom  = (TYPE_SEEDT *)read_table_each(sizeof(TYPE_SEEDT), size[par.kmer_ba], idata_a, ".seedtable-revcom");
  /*----- blktable_a_revcom -----*/
  TYPE_BLKT *blktable_sum = (TYPE_BLKT *)read_table_each(sizeof(TYPE_BLKT), idata_a->ksum[1][idata_a->cnt], idata_a, ".blktable-revcom");
  blktable_a_revcom = alloc_blktable(blktable_sum, seedtable_a_revcom, size[par.kmer_ba]);
  free(blktable_sum);
  return;
}

static void chomp(char *str){
  char *p = strchr(str, '\n');
  if(p) p[0]='\0';
  return;
}

int ParseLine(char *str, struct elem clm[]){
  int i, j=0, num=0, len=strlen(str);
  char *strtemp = (char *)my_calloc(len, sizeof(char), "ParseLine");
  for(i=0; i<=len; i++){
    if(str[i]=='\0' || str[i]=='\n'){
      strtemp[j]='\0';
      strcpy(clm[num].str, strtemp);
      free(strtemp);
      return ++num;
    }
    if(str[i]=='\t'){
      strtemp[j]='\0';
      strcpy(clm[num].str, strtemp);
      num++; 
      j=0;
    }else{
      strtemp[j]=str[i];
      j++;
    }
  }
  free(strtemp);
  return num;
}

void read_summary(inputdata *idata){
  int i;
  struct elem clm[10];
  char filename[INPUTFILENAME_MAX], str[256];
  sprintf(filename, "%s/%s-%d-%d.txt", par.tabledir, checkfilename(idata->seqname), par.kmer_ba, par.blocksize);
  FILE *IN = my_fopen_r(filename);

  if(fgets(str, 256, IN)){ /* first line */
    chomp(str);
    if(strcmp(str, idata->seqname)){
      fprintf(stderr, "wrong table: %s <-> %s\n", str, idata->seqname);
      exit(0);
    }
  }

  while(fgets(str, 256, IN)){
    if(!strcmp(str, "\n")) continue;
    chomp(str);
    if(!strncmp(str, "fasta num: ", 11)){
      idata->fstnum      = atoi(str+11);
      idata->blocknum    = (int *)my_malloc(idata->fstnum * sizeof(int), "idata->blocknum");
      idata->length_each = (int *)my_malloc(idata->fstnum * sizeof(int), "idata->length_each");
      for(i=0; i<2; i++) idata->ksum[i] = (int *)my_malloc(idata->fstnum * sizeof(int), "idata->ksum[i]");
    }else if(!strncmp(str, "total length: ", 14)) idata->length_total = atol(str+14);
    else{
      ParseLine(str, clm);
      idata->length_each[atoi(clm[0].str)] = atoi(clm[1].str);
      idata->blocknum[atoi(clm[0].str)]    = atoi(clm[2].str);
      idata->ksum[0][atoi(clm[0].str)]     = atoi(clm[3].str);
      idata->ksum[1][atoi(clm[0].str)]     = atoi(clm[4].str);
    }
  }
  fclose(IN);
}
