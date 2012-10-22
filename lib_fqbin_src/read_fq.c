
#include "lib_fqbin.h"
#include <stdio.h>
#include <ctype.h>

#include <unistd.h>


// process a fastq file adding it to fbin file
int find_in_fastq(char *fname, char *seq_name)
{

    // allocate strings
    char *name;
    if ((name = malloc(MAXSEQNAME)) == NULL) {
        puts("Memory allocation error!");
        return EXIT_FAILURE;
    }
    
    char *fasta;
    if ((fasta = malloc(MAXSEQLENGTH)) == NULL) {
        puts("Memory allocation error!");
        return EXIT_FAILURE;
    }
    
    char *qual;
    if ((qual = malloc(MAXSEQLENGTH)) == NULL) {
        puts("Memory allocation error!");
        return EXIT_FAILURE;
    }
    
    char *comments;
    if ((comments = malloc(MAXSEQLENGTH)) == NULL) {
        puts("Memory allocation error!");
        return EXIT_FAILURE;
    }
    
    static time_t curr_time=0;
    static time_t prev_time=0;
    
    prev_time=time(NULL);
    
    FILE *fastq_file=NULL;
    
    
    int valid=0;
    int res=0;
    int r=0;
    // Open fasta and qual files
    if (strcmp(fname,"-")==0){
        fastq_file=stdin;
    }else{

        open_file(fname,&fastq_file);
    }
    
    if (fastq_file==NULL){
        printf("TRESb\n");
    }
        
    // for each sequence on fastq file
    while (valid=get_next_seq_fastq(fastq_file,&name,&fasta,&qual,&comments)){
        if(valid==1)
        {
            r++;
            if (strcmp(name,seq_name)==0)
            {
                printf("@%s %s\n%s\n", name,comments, fasta);
                printf("+%s\n%s\n",name,qual);
                break;
            }

        }else{
            fprintf(stderr,"Invalid sequence found %s. Aborting import.\n",name);
            res=-1;
            break;
        }
        
        // if ((r%10000)==0) {
        // }
        
    }

    curr_time=time(NULL);    
    printf("\nEnd fastq processing. %d seqs in %.0f s. Rate: %8.2f seqs/s\n",r,difftime(curr_time,prev_time),r/difftime(curr_time,prev_time));

    // free mem
    free(name);
    free(fasta);
    free(qual);
    free(comments);
    
    // close files
    fclose(fastq_file);
    
    return res;
}

void usage(){
    printf("Usage: read_fq fastq_file seq_name\n\n");
    
    exit(-1);
    
}



/*******************************************************/
/* main                                                */
/*******************************************************/
int main(int argc, char *argv[])
{
    
    int ch;

    while ((ch = getopt(argc, argv, "h")) != -1) {
            switch (ch) {

            case 'h':
                    usage();
                    break;
            case '?':
            default:
                    usage();
            }
    }

    argc -= optind;
    argv += optind;
    
    
  // check params
  if (argc<2)
  {
    usage();
    return -1;
  }

  int res=find_in_fastq(argv[0],argv[1]);

  return res;
}

