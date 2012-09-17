#include "libfbin.h"
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>



#include <unistd.h>


void usage(){
    printf("Usage: fbin_generate_index fbin_file\n\n");
    // printf("    -f    Output sequence in fasta format\n\n");
    
    exit(-1);
    
}

/*******************************************************/
/* main                                                */
/*******************************************************/
int main(int argc, char *argv[])
{

    char *fasta=NULL;
    char *qual=NULL;
    char *extras=NULL;
    int size=5000;
    int res=0;

    int ch;

    int output_fasta = 0;
    int output_qual = 0;

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
    // printf("argc: %d", argc);
    // printf("argv: %s", argv[0]);

    if (argc!=1)
    {
        usage();
    }


    if (regenerate_index(argv[0])==-1){
        printf("File %s does not exists",argv[0]);
        exit(-1);
    }

    exit(0);

}

