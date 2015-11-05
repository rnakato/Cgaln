/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of Cgaln-1.*.* sources.
 */
#include "BA.h"
#include "NA.h"
#include <time.h>

#define D_PRIME 37  /* for double hash */

static void CGAT();
static void Make_SeedTable(Fasta *, int **, unsigned short **);
static void argv_init(int, char **);
static void define_xdrop_and_shiftp(int, int);
static void dump_init();

static const char Usage[] =
"./Cgaln <sequence1> <sequence2> -o <outputfile>\n\n\
       -t: table directory (default: CgalnTable)\n\
       -r: both strand (default: forward strand)\n\
       -b: output block-level alignment and exit (gnuplot format only) \n\
       -nc: no chaining (default: with chaining)\n\
       -ia: with iterative alignment (default: off)\n\
       -pr: print squares of regions for iterative alignment (default: off)\n\
       -fc: filter colony to extract consistent set (default: off)\n\
       -cons: filter inconsistent HSPs at the HSP-chaining (defaults: off)\n\
       -noext: (when -ia is on) omit extending HSPs with gapped DP (default: do not omit)\n\
       -sr: consider short-reverse complement (default: off)\n\
       -k: seed weight; -k1: 11of18; -k2: 12of19; -k3: 13of20 (default: 11of18)\n\
       -BS(10000): block size\n\
       -X(4000): X drop-off at block-level\n\
       -R: ratio between X drop-off score and gap penalty at block-level (default: 1500)\n\
       -nl: don't show the lines of fasta-end for multifasta file (default: on)\n\
       -otype<int>: output format; 0: gnuplot; 1: ensfile; 2: fasta (default: gnuplot)\n\
       -debug: debug mode\n\n";

int main(int argc, char *argv[]){
  clock_t start, end;
  start = clock();

  /*--- initialization ---*/
  argv_init(argc, argv);
  printf("genomeA: %d fasta, genomeB: %d fasta\n", idata_a->fstnum, idata_b->fstnum);

  /*-------- CGAT algorithm ---------*/
  CGAT();

  idata_delete(idata_a);
  idata_delete(idata_b);

  end = clock();
  if(opt.debug) printf("total time: %.2f sec.\n", (double)(end-start)/CLOCKS_PER_SEC);

  return 0;
}

/* CGAT
 * read and make tables
 * BA -> NA 
 */
static void CGAT(){
  FILE *IN=NULL, *IN_b=NULL, *IN_rev=NULL;
  int **table_value=NULL;
  unsigned short **table_num=NULL;
  int state_a2, state_b=0, state_rev;
  Fasta *fst1 = fasta_new(), *fst2 = fasta_new(), *fst_rev = fasta_new();
  clock_t start, start1, end0, end1, end2;

  /*----- do BA and NA for each fasta-pair-----*/
  /* genomeB file open */
  IN_b = my_fopen_r(idata_b->seqname);
  par.OUT = my_fopen_w(par.outputfile);

  for(idata_b->cnt=0; idata_b->cnt < idata_b->fstnum; idata_b->cnt++){
    start = clock();

    /* MakeTable */
    read_multifasta(IN_b, fst2, FORWARD, &state_b);
    table_value = (int **)my_malloc(idata_b->blocknum[idata_b->cnt] * sizeof(int *), "table_b_value");
    table_num = (unsigned short **)my_malloc(idata_b->blocknum[idata_b->cnt] * sizeof(unsigned short *), "table_b_num");
    Make_SeedTable(fst2, table_value, table_num);

    end0 = clock();
    if(opt.debug) printf("MakeTable time: %.2f sec.\n", (double)(end0-start)/CLOCKS_PER_SEC);

    IN     = my_fopen_r(idata_a->seqname);
    IN_rev = my_fopen_r(idata_a->seqname);
    state_a2=0; state_rev=0;

    for(idata_a->cnt=0; idata_a->cnt < idata_a->fstnum; idata_a->cnt++){      
      printf("\ngenomeA-fasta%d (%d blocks) - genomeB-fasta%d (%d blocks)\n",
	     idata_a->cnt+1, idata_a->blocknum[idata_a->cnt], idata_b->cnt+1, idata_b->blocknum[idata_b->cnt]);
      start1 = clock();

      /*--- BA: the results are stored in aln_for/rev ---*/
      BA(&aln_for, table_value, table_num, FORWARD);
      if(reverse) BA(&aln_rev, table_value, table_num, REVERSE);
      if(idata_a->cnt == idata_a->fstnum -1) table_b_delete(table_value, table_num, idata_b->blocknum[idata_b->cnt]);

      end1 = clock();
      if(opt.debug) printf("BA time: %.2f sec.\n", (double)(end1-start1)/CLOCKS_PER_SEC);
      /*--- (if -b is on) output BA result and skip NA ---*/
      if(block){
	output_BAresult();
	continue;
      }

      /*--- NA: detailed alignmend within colonies in bl ---*/
      NA(IN, fst1, fst2, &aln_for, FORWARD, &state_a2);
      if(reverse) NA(IN_rev, fst_rev, fst2, &aln_rev, REVERSE, &state_rev);

      end2 = clock();
      if(opt.debug) printf("NA time: %.2f sec.\n", (double)(end2-end1)/CLOCKS_PER_SEC);
    }

    free(fst2->head);
    free(fst2->body);
    fclose(IN);
    fclose(IN_rev);
  }

  if(opt.boundary) output_fastaboundary();

  free(fst1); free(fst2); free(fst_rev);
  fclose(IN_b);
  fclose(par.OUT);
}

static int hashing(int seedvalue, struct seed *hash, int hashsize){
  int h, h2 = D_PRIME;
  int cnt = 0; 
  while(1){
    h = (seedvalue + h2*cnt) % hashsize;
    if(hash[h].value == -1){
      hash[h].value = seedvalue;
      hash[h].num = 1;
      break;
    }else if(hash[h].value == seedvalue){
      hash[h].num++;
      break;
    }else{
      cnt++;
    }
  }
  return 0;
}

static int primes(int max){ 
  int m, i;
  for(m=max; m>=2; m--){
    for(i=3; i<max/2; i+=2){
      if(!(m%i)) break;
    }
    if(i>=max/2) break;
  }
  return m;
}

static void Make_SeedTable(Fasta *fst, int **table_value, unsigned short **table_num){
  int i,j;
  int hashsize = primes((int)(par.blocksize*1.2));
  int start, end, num, seedvalue, M_b=idata_b->blocknum[idata_b->cnt];
  struct seed *hash = (struct seed *)my_malloc(hashsize * sizeof(struct seed), "hash");
  struct seed *seedlist_temp = (struct seed *)my_malloc(par.blocksize * sizeof(struct seed), "seedlist_temp");

  for(i=0; i<M_b; i++){   
    for(j=0; j<hashsize; j++) hash[j].value = -1;
    start = i*par.blocksize;
    if(i==M_b-1) end = fst->length - window[par.kmer_ba];
    else end = (i+1)*par.blocksize;

    for(j=start; j<end; j++){
      seedvalue = SeedScoring_masked(fst, j, par.kmer_ba);
      if(seedvalue != -1) hashing(seedvalue, hash, hashsize);
    }   

    num=0;
    for(j=0; j<hashsize; j++){
      if(hash[j].value != -1) seedlist_temp[num++] = hash[j];
    }
    table_value[i] = (int *)my_malloc((num+1) * sizeof(int), "table_b_value[i]");
    table_num[i] = (unsigned short *)my_malloc((num+1) * sizeof(unsigned short), "table_b_num[i]");
    for(j=0; j<num; j++){
      table_value[i][j] = seedlist_temp[j].value;
      table_num[i][j] = seedlist_temp[j].num;
    }
    table_value[i][num] = -1;    /* guard */
  }

  free(seedlist_temp);
  free(hash);
  return;
}


static void global_variable_initialization(){  
  reverse=0; block=0; iterative=0; ext_gappeddp=1; short_reverse=0; 
  printregion=0; filtercln=0; consistent=0;
  opt.debug=0;
  opt.boundary=1;
  opt.chaining=1;
  opt.otype=0;

  par.kmer_ba = SEEDWEIGHT_DEFAULT;
  par.kmer_na = SEEDWEIGHT_DEFAULT;
  par.blocksize = BLOCKSIZE_DEFAULT;
  par.xdrop = XDROP_DEFAULT;
  par.x_d_ratio = X_D_RATIO_DEFAULT;
  strcpy(par.tabledir, TABLE_PATH);
  strcpy(par.outputfile, "");
  par.OUT=NULL;

  idata_a = idata_new();
  idata_b = idata_new();

  aln_for.bl=NULL;
  aln_rev.bl=NULL;
}

static void argv_init(int argc, char **argv){
  int i;
  int nseq=0;
  global_variable_initialization();

  for(i=1; i<argc; i++){
    if(!strcmp(argv[i], "-r"))            reverse = 1;
    else if(!strcmp(argv[i], "-b"))       block = 1;
    else if(!strcmp(argv[i], "-debug"))   opt.debug = 1;
    else if(!strcmp(argv[i], "-nc"))      opt.chaining = 0;
    else if(!strcmp(argv[i], "-ia"))      iterative = 1;
    else if(!strcmp(argv[i], "-noext"))   ext_gappeddp = 0;
    else if(!strcmp(argv[i], "-sr"))      short_reverse = 1;
    else if(!strcmp(argv[i], "-pr"))      printregion = 1;
    else if(!strcmp(argv[i], "-fc"))      filtercln = 1;
    else if(!strcmp(argv[i], "-cons"))    consistent = 1;
    else if(!strcmp(argv[i], "-k1"))      par.kmer_ba = 11;
    else if(!strcmp(argv[i], "-k2"))      par.kmer_ba = 12;
    else if(!strcmp(argv[i], "-k3"))      par.kmer_ba = 13;
    else if(!strcmp(argv[i], "-t"))       strcpy(par.tabledir, argv[++i]);
    else if(!strcmp(argv[i], "-o"))       strcpy(par.outputfile, argv[++i]);
    else if(!strcmp(argv[i], "-nl"))      opt.boundary = 0;
    else if(!strncmp(argv[i], "-K", 2))   par.kmer_ba   = atoi(argv[i]+2);
    else if(!strncmp(argv[i], "-BS", 3))  par.blocksize = atoi(argv[i]+3);
    else if(!strncmp(argv[i], "-X", 2))   par.xdrop     = atoi(argv[i]+2);
    else if(!strncmp(argv[i], "-R", 2))   par.x_d_ratio = atoi(argv[i]+2);
    else if(!strncmp(argv[i], "-otype", 6)) opt.otype   = atoi(argv[i]+6);
    else{
      if(nseq==0){
	strcpy(idata_a->seqname, argv[i]);
	nseq++;
      }else if(nseq==1){
	strcpy(idata_b->seqname, argv[i]);
	nseq++;
      }else{
	nseq++;
      }
    }
  }
  if(nseq!=2){
    fprintf(stderr, "irregal inputfile number:%d\n", nseq);
    goto err;
  }
  if(!strcmp(par.outputfile, "")){
    fprintf(stderr, "please specify outputfile name.\n");
    goto err;
  }
  define_xdrop_and_shiftp(par.kmer_ba, par.blocksize);

  if(block) opt.otype=0;

  /*--- read summary ---*/
  read_summary(idata_a);
  read_summary(idata_b);
  dump_init();
  return;

 err:
  printf(Usage);
  exit(0);
}

static char *str_boolean[2]={ "OFF", "ON"};
static char *str_output[3]={ "gnuplot", "ensfile", "fasta"};
static void dump_init(){
  printf("============================\n");
  printf("genome A: %s (%d fasta)\n", idata_a->seqname, idata_a->fstnum);
  printf("genome B: %s (%d fasta)\n", idata_b->seqname, idata_b->fstnum);
  printf("output file: %s\n", par.outputfile);
  printf("table directory: %s\n", par.tabledir);
  printf("reverse complement: %s\n", str_boolean[reverse]);
  printf("<for BA (block-level alignment) option>\n");
  printf("\tBA seed length: %d-mer\tseed size: %d\n", par.kmer_ba, size[par.kmer_ba]);
  printf("\tblock size: %.1fkbp\n",   par.blocksize/(double)1000);
  printf("\tfilter inconsistent colonies: %s\n", str_boolean[filtercln]);
  if(block) printf("output block level result.\n");
  else{
    printf("<for NA (nucleotide-level alignment) option>\n");
    printf("\tNA seed length: %d-mer\tseed size: %d\n", par.kmer_na, size[par.kmer_na]);
    printf("\tHSP chaining: %s\n", str_boolean[opt.chaining]);
    if(opt.chaining){
      printf("\tfilter inconsistent HSPs: %s\n", str_boolean[consistent]);
      printf("\titerative alignment: %s\n", str_boolean[iterative]);
      if(iterative) printf("\tshort-reverse: %s\n", str_boolean[short_reverse]);
    }
  }
  printf("output format: %s\n", str_output[opt.otype]);
  printf("============================\n\n");
}


static void define_xdrop_and_shiftp(int kmer, int blocksize){
  switch(kmer){
  case 11: break;
  case 12: par.xdrop = par.xdrop * 11/25; break;
  case 13: par.xdrop = par.xdrop * 4/25;  break;
  }

  double bsratio = blocksize / (double)10000;
  double a = 3/14.0;
  double b = 11/14.0;
  par.xdrop = par.xdrop * (a * bsratio* bsratio + b * bsratio);
  par.shiftp = par.xdrop * 100 / par.x_d_ratio;

  return;
}
