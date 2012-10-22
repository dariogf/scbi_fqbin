#include "lib_fqbin.h"
#include <stdio.h>
#include <ctype.h>

#include <unistd.h>


void usage(){
  
    // printf("Usage: mk_fqbin [-i] [-f flatten_qual] [-d discretize_qual] [-e extras_file] -o output_file [fastq_file_input_file]\n\n");
    
    printf("mk_fqbin converts a fastQ input file or STDIN stream (use no filename or '-') to compressed fqbin format.\n\n");
        
    printf("Usage: mk_fqbin [OPTIONS] -o output_file [fastq_file_input_file]\n\n");
    
    printf("Options:\n");
    printf("    -i create random access index\n");
    printf("    -d discretize_qual: quality values are discretized in groups of size discretize_qual. This way less quality values are used and a better compression is obtained\n");
    printf("    -f flatten_qual: quality values over flatten_qual (use phred scale) will be set to flatten_qual value in order to achieve a better compression\n");
    printf("    -e extras_file: a file with extra metadata for each sequence if standard FASTA format\n");

    printf("Mandatory parameters:\n");
    printf("    -o output_file: output fqbin file\n");
    printf("    -F input is in fasta format (will look for filename.qual for qualities) \n");
    
    printf("\nSCBI - Supercomputación y Bioinformática. University of Malaga. http://www.scbi.uma.es. Copyright 2011\n\n");
    
    exit(-1);
    
}

/*******************************************************/
/* main                                                */
/*******************************************************/
int main(int argc, char *argv[])
{
  
  int res = 0;
  
  int ch;

  int output_fasta = 0;
  int output_qual = 0;
  
  int flatten_qual=0;
  int discretize_qual=0;
  
  char *extras_file=NULL;
  char *output_file=NULL;
  int create_index=0;
  int input_in_fasta=0;

  while ((ch = getopt(argc, argv, "o:e:d:f:iFh")) != -1) {
      // printf("opt %d\n",ch);
          switch (ch) {
          case 'e':
                  // strcopy(extras_file,optarg);
                  extras_file=optarg;
                  break;
          case 'f':
                  flatten_qual = atoi(optarg)+33;
                  break;
          case 'd':
                  discretize_qual = atoi(optarg);
                  
                  if(discretize_qual<2)
                  {
                    discretize_qual=0;
                  }
                  break;
          case 'o':
                  output_file=optarg;
                  break;
          case 'i':
                  create_index=1;
                  break;
          case 'F':
                  input_in_fasta=1;
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
  
  if(output_file==NULL)
  {
      printf("Output file is a mandatory option. Provide one with -o filename\n"); 
      usage();
      exit -1;
  }
  
  printf("Extra metadata file: %s\nFlattenning qual over:%d (%c char)\nDistretizing qual in groups of:%d\n",extras_file,flatten_qual,flatten_qual,discretize_qual);
  
  if (create_index){printf("Creating random access index\n");}
  

   if(input_in_fasta)
   {

       // check remaining params
       if (argc==1){
         res=process_fasta(argv[0],extras_file,output_file,discretize_qual,flatten_qual,create_index);
       }
       else if(argc==0)
       {
           res=process_fasta("-",extras_file,output_file,discretize_qual,flatten_qual,create_index);
       }
       else{
         usage();
       }
       
   }else{
       // check remaining params
       if (argc==1){
         res=process_fastq(argv[0],extras_file,output_file,discretize_qual,flatten_qual,create_index);
       }
       else if(argc==0)
       {
           res=process_fastq("-",extras_file,output_file,discretize_qual,flatten_qual,create_index);
       }
       else{
         usage();
       }
       
   }
 
  return res;
}

