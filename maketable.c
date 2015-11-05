/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of Cgaln-1.*.* sources.
 */
#include "util.h"
#include "readfasta.h"
#include <math.h>
#include <sys/stat.h>

static TYPE_SEEDT *seedtable; 
static struct poistable *poistable;
static TYPE_BLKT **blktable;
static TYPE_BLKT *blktable_sum;

static int make_seedtable(Fasta *, int);
static int write_table(Fasta *, char *, int);
static void write_table_each(void *, size_t, int, char *, char *);
static double pois(int, double);
static double pois2(int, double *);
static void argv_init(int, char **, char *);
static void Output_summary(char *, int **);

static const char Usage[] =
"./maketable [options] <sequence>\n\
       -o specify output directory (default: CgalnTable)\n\
       -K# specify k-mer size (11 as default)\n\
       -BS# specify blocksize (10000 as default)\n\n";

static void rm_table(char *head){
  struct stat st;
  char filename[INPUTFILENAME_MAX];
  sprintf(filename, "%s.poistable-%d", head, par.blocksize);
  if(!stat(filename, &st)) remove(filename);
  sprintf(filename, "%s.seedtable-%d", head, par.kmer_ba);
  if(!stat(filename, &st)) remove(filename);
  sprintf(filename, "%s.seedtable-revcom-%d", head, par.kmer_ba);
  if(!stat(filename, &st)) remove(filename);
  sprintf(filename, "%s.blktable-%d-%d", head, par.kmer_ba, par.blocksize);
  if(!stat(filename, &st)) remove(filename);
  sprintf(filename, "%s.blktable-revcom-%d-%d", head, par.kmer_ba, par.blocksize);
  if(!stat(filename, &st)) remove(filename);
}

int main(int argc, char *argv[]){
  FILE *IN;
  int i, strand, state;
  char seqname[INPUTFILENAME_MAX], prefix[INPUTFILENAME_MAX];

  argv_init(argc, argv, seqname);
  sprintf(prefix, "%s/%s", par.tabledir, checkfilename(seqname));
  idata_a->length_each = (int *)my_calloc(idata_a->fstnum, sizeof(int), "idata.length_each");
  int *ksum[2];
  ksum[0] = (int *)my_calloc(idata_a->fstnum, sizeof(int), "ksum");
  ksum[1] = (int *)my_calloc(idata_a->fstnum, sizeof(int), "ksum_rev");
  Fasta *fst = fasta_new();

  mkdir(par.tabledir, 0755);
  rm_table(prefix);

  /*------ read fastafile ------*/
  for(strand=0; strand<2; strand++){
    if(!strand) printf("<forward strand>\n"); else printf("<reverse strand>\n");
    IN = my_fopen_r(seqname);
    state=0;

    for(i=0; i<idata_a->fstnum; i++){
      read_multifasta(IN, fst, strand, &state);
      printf("sequence %d\n", i+1);
      printf(">%s\n", fst->head);
      printf("length: %d (%d blocks)\n\n", fst->length, fst->length/par.blocksize +1);
      idata_a->length_each[i] = fst->length;
      
      ksum[strand][i] = write_table(fst, prefix, strand);

      free(fst->head);
      free(fst->body);
    }
    fclose(IN);
  }
  free(fst);

  Output_summary(seqname, ksum);
  free(ksum[0]); free(ksum[1]);
  idata_delete(idata_a);

  return 0;
}

static void Output_summary(char *seqname, int **ksum){
  int i;
  FILE *OUT;
  char filename[INPUTFILENAME_MAX];
  sprintf(filename, "%s/%s-%d-%d.txt", par.tabledir, checkfilename(seqname), par.kmer_ba, par.blocksize);
  OUT = my_fopen_w(filename);
  fprintf(OUT, "%s\n", seqname);
  fprintf(OUT, "fasta num: %d\n", idata_a->fstnum);
  for(i=0; i<idata_a->fstnum; i++) idata_a->length_total += idata_a->length_each[i];
  fprintf(OUT, "total length: %ld\n", idata_a->length_total);
  for(i=0; i<idata_a->fstnum; i++){
    fprintf(OUT, "%d\t%d\t%d\t%d\t%d\n", i, idata_a->length_each[i], idata_a->length_each[i]/par.blocksize+1, ksum[0][i], ksum[1][i]);
  }
  fclose(OUT);
}

static int write_table(Fasta *fst, char *head, int strand){
  int ksum = make_seedtable(fst, strand);
  if(!strand){
    write_table_each(seedtable,    sizeof(TYPE_SEEDT), size[par.kmer_ba], head, ".seedtable");
    write_table_each(blktable_sum, sizeof(TYPE_BLKT), ksum, head, ".blktable");
    write_table_each(poistable,    sizeof(struct poistable), KMER_NUMMAX*POISSON_MAX, head, ".poistable");
  }else{
    write_table_each(seedtable,    sizeof(TYPE_SEEDT), size[par.kmer_ba], head, ".seedtable-revcom");
    write_table_each(blktable_sum, sizeof(TYPE_BLKT), ksum, head, ".blktable-revcom");
  }
  return ksum;
}


static void write_table_each(void *table, size_t size, int tablesize, char *head, char *postfix){
  char filename[INPUTFILENAME_MAX];
  if(!strcmp(postfix, ".poistable")){
    sprintf(filename, "%s%s-%d", head, postfix, par.blocksize);
  }else if(!strcmp(postfix, ".seedtable") || !strcmp(postfix, ".seedtable-revcom")){
    sprintf(filename, "%s%s-%d", head, postfix, par.kmer_ba);
  }else{
    sprintf(filename, "%s%s-%d-%d", head, postfix, par.kmer_ba, par.blocksize);
  }
  //printf("%s, tablesize=%ld\n", filename, (long)(size*tablesize));

  FILE *IN = my_fopen_ab(filename);
  if(fwrite(table, size*tablesize, 1, IN)!=1){
    fprintf(stderr,"[E] fwrite error:%s\n", filename); 
    exit(1);
  }
  fclose(IN);
  free(table);
  return;
}

static int make_seedtable(Fasta *fst, int strand){
  int i,j, seedvalue, ksum=0, num=0;
  int end = fst->length - window[par.kmer_ba];
  seedtable = (TYPE_SEEDT *)my_calloc(size[par.kmer_ba], sizeof(TYPE_SEEDT), "seedtable");
  blktable  = (TYPE_BLKT **)my_malloc(size[par.kmer_ba] * sizeof(TYPE_BLKT *), "blktable");

  for(i=0; i<=end; i++){
    seedvalue = SeedScoring_masked(fst, i, par.kmer_ba);
    if(seedvalue != -1 && seedtable[seedvalue] != -1){
      seedtable[seedvalue]++;
      if(seedtable[seedvalue] >=KMER_NUMMAX) seedtable[seedvalue] = -1;
    }
  }

  for(i=0; i<size[par.kmer_ba]; i++){
    if(seedtable[i]>0){
      blktable[i] = (TYPE_BLKT *)my_malloc(seedtable[i] * sizeof(TYPE_BLKT), "blktable[i]");
      ksum += seedtable[i];
    }else{
      blktable[i] = NULL;
    }
    seedtable[i] = 0; /* reinitialize*/
  }
  for(i=0; i<=end; i++){
    seedvalue = SeedScoring_masked(fst, i, par.kmer_ba);
    if(seedvalue != -1){
      if(blktable[seedvalue]){
	blktable[seedvalue][seedtable[seedvalue]] = i/par.blocksize;
	seedtable[seedvalue]++;
      }
    }
  }

  blktable_sum = (TYPE_BLKT *)my_malloc(ksum * sizeof(TYPE_BLKT), "blktable_sum");
  for(i=0; i<size[par.kmer_ba]; i++){
    for(j=0; j<seedtable[i]; j++){
      blktable_sum[num++] = blktable[i][j];
    }
  }
  blktable_delete(blktable, size[par.kmer_ba]);

  /*--- poisson table (forward only) ---*/
  if(!strand){
    poistable = (struct poistable *)my_calloc(KMER_NUMMAX*POISSON_MAX, sizeof(struct poistable), "poistable");
    int M_num = fst->length/par.blocksize;
    double myu, temp[POISSON_MAX];
    for(i=0; i<KMER_NUMMAX; i++){
      myu = i / (double)M_num;
      for(j=0; j<POISSON_MAX; j++){
	temp[j] = pois(j, myu);
	poistable[i*POISSON_MAX +j].first  = TEMP_TERM * log(temp[j]);
	poistable[i*POISSON_MAX +j].second = TEMP_TERM * pois2(j, temp);
      }
    }
  }

  return ksum;
}

/* Poisson */
static double pois(int k, double myu){
  int i;
  double term1=1, term2=1, term3=1/exp(myu);
  if(!k) return (term3);
  for(i=1; i<=k; i++){
    term1 = term1 * myu;
    term2 = term2 * i;
  }
  return (term1*term2*term3);
}

static double pois2(int k, double *temp){
  int i;
  double p=1;
  for(i=0; i<k; i++) p -= temp[i];
  return (log(p));
}

static void dump_init(char *seqname){
  printf("============================\n");
  printf("inputfile: %s\n", seqname);
  printf("number of fasta sequences: %d\n", idata_a->fstnum);
  printf("output directory: %s\n", par.tabledir);
  printf("seed length: %d-mer\n", par.kmer_ba);
  printf("seed size: %d\n", size[par.kmer_ba]);
  printf("block size: %.1fkbp\n", par.blocksize/(double)1000);
  printf("============================\n\n");
}

static void global_variable_initialization(){
  par.kmer_ba = SEEDWEIGHT_DEFAULT;  
  par.blocksize = BLOCKSIZE_DEFAULT;
  strcpy(par.tabledir, TABLE_PATH);
  idata_a = idata_new();
}

static void argv_init(int argc, char **argv, char *seqname){
  int i;
  int nseq=0;

  global_variable_initialization();

  for(i=1; i<argc; i++){
    if(!strncmp(argv[i], "-BS", 3))     par.blocksize = atoi(argv[i]+3);
    else if(!strncmp(argv[i], "-K", 2)) par.kmer_ba = atoi(argv[i]+2);
    else if(!strcmp(argv[i], "-o"))     strcpy(par.tabledir, argv[++i]);
    else{
      nseq++;
      if(nseq > 1) goto err;
      strcpy(seqname, argv[i]);
    }
  }
  if(nseq!=1) goto err;

  idata_a->fstnum = countfasta(seqname);
  dump_init(seqname);
  return;

 err:
  printf("irregal inputfile number:%d\n", nseq);
  printf(Usage);
  exit(0);
}
