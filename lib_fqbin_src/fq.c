
#include "libfbin.h"
#include <stdio.h>
#include <ctype.h>

#include <unistd.h>


// process a fastq file adding it to fbin file
int iterate_fastq(char *fname, int only_extras, int output_fasta, int output_extras)
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

                if (!only_extras){
                     if (output_fasta){
                         printf(">%s %s\n%s\n", name, comments, fasta);
                     }else{
                         printf("@%s %s\n%s\n", name,comments, fasta);
                       printf("+%s\n%s\n",name,qual);
                     }
                }
                
                // if ((extras!=NULL) && (output_extras)) printf ("EXTRAS:%s\n",extras);
            
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
    printf("Usage: fq [-f][-e][-E] fbin_file seq_name\n\n");
    printf("    -f    Output sequence in fasta format\n");
    printf("    -e    Output extras for sequence\n");
    printf("    -E    Output only extras for sequence\n");
    
    exit(-1);
    
}



/*******************************************************/
/* main                                                */
/*******************************************************/
int main(int argc, char *argv[])
{
    
    int ch;

    int output_fasta = 0;
    int output_extras = 0;
    int only_extras = 0;

    while ((ch = getopt(argc, argv, "feEh")) != -1) {
            switch (ch) {
            case 'f':
                    output_fasta = 1;
                    break;
            case 'e':
                    output_extras = 1;
                    break;
            case 'E':
                    output_extras = 1;
                    only_extras = 1;
                    break;

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
  if (argc<1)
  {
    usage();
    return -1;
  }

  int res=iterate_fastq(argv[0],only_extras, output_fasta, output_extras);

  return res;
}

