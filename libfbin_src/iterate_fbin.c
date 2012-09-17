#include "libfbin.h"
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>


void usage(){
    
    printf("Usage: iterate_fbin [-F|-q] fbin_file\n\n");
    printf("By default outputs in fastq format\n\n");
    printf("    -F    Output only sequence in fasta format\n");
    printf("    -q    Output only qualities in phred format\n\n");
    printf("    -e    Output only extras\n\n");    
    
    exit(-1);
    
}

int print_file(struct file_data *filed, int only_fasta, int only_qual, int only_extras){
    char *sname=NULL;
    char *fasta=NULL;
    char *qual=NULL;
    char *extras=NULL;
    
    int len=0;
    int i=0;
    // int size=5000;
    int res=0;
    
    
    while ((res=read_data_sequential(filed, &sname, &fasta, &qual, &extras))==0)
    {

            if (only_fasta){
                
                printf(">%s %s\n", sname, extras);
                len=strlen(fasta);
                i=0;
                for(i = 0; i < len; i+=70)
                {
                    printf("%.70s\n", fasta+i);
                }
            }else if (only_qual){
                printf(">%s %s\n", sname, extras);
                if (qual!=NULL){
                    len=strlen(qual);
                    i=0;
                    for(i = 0; i < len; i++)
                    {
                        printf("%02d ", qual[i]-33);
                        if (((i+1)%30 == 0) || (i==len-1)) printf("\n");
                    }
                }
                
            }else if (only_extras){
                printf(">%s %s\n", sname, extras);
                if (extras!=NULL) printf ("%s\n",extras);
            }else{
                printf("@%s %s\n%s\n", sname, extras, fasta);
                         printf("+\n%s\n",qual);
                            // printf("+%s\n%s\n",sname,qual);
            }
            

    	if ( fasta!=NULL ) {free(fasta);fasta=NULL;}
    	if ( qual!=NULL ) {free(qual);qual=NULL;}
    	if ( extras!=NULL ) {free(extras);extras=NULL;}
    }
    
    return res;
}



/*******************************************************/
/* main                                                */
/*******************************************************/
int main(int argc, char *argv[])
{

    //gzFile gzf_bin;
    // struct file_data filed;

    struct file_data *filed=NULL;

    int ch;

    int output_fasta = 0;
    int output_qual = 0;
    int output_extras = 0;

    while ((ch = getopt(argc, argv, "Fqeh")) != -1) {
            switch (ch) {
            case 'F':
                    output_fasta = 1;
                    break;
            case 'q':
                    output_qual=1;
                    break;
            case 'e':
                    output_extras=1;
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
    // printf("argc: %d", argc);
    // printf("argv: %s", argv[0]);

    if (argc!=1)
    {
        usage();
    }

    // initialize reads
    if (initialize_sequential_reads(&filed, argv[0])==-1){
        printf("File %s does not exists",argv[0]);
        exit(-1);
    }

    int res=print_file(filed,output_fasta,output_qual,output_extras);

    close_sequential_reads(filed);

    return res;
}

