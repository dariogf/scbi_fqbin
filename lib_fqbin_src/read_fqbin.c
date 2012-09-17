#include "libfbin.h"
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>

#include <unistd.h>


void usage(){
    printf("Usage: read_fbin [-f][-e][-E] fbin_file seq_name\n\n");
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

	char *fasta=NULL;
	char *qual=NULL;
	char *extras=NULL;
	int size=5000;
	int res=0;

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
    // printf("argc: %d", argc);
    // printf("argv: %s", argv[0]);


	if (argc!=2)
	{
        usage();
	  return -1;
	}

res=read_seq(argv[0],argv[1], &fasta, &qual, &extras);



  if (res==0){
  
      if (!only_extras){
           if (output_fasta){
               printf(">%s\n%s\n", argv[1], fasta);
           }else{
               printf("@%s\n%s\n", argv[1], fasta);
             printf("+%s\n%s\n",argv[1],qual);
           }
   }
      if ((extras!=NULL) && (output_extras)) printf ("EXTRAS:%s\n",extras);
   
  }

	//res=read_seq("prueba2gz.fbin","F143CJN01EN6AH", &fasta, &qual, &extras);
	if ( fasta!=NULL ) {free(fasta);fasta=NULL;}
	if ( qual!=NULL ) {free(qual);qual=NULL;}
	if ( extras!=NULL ) {free(extras);extras=NULL;}
	return res;
}

