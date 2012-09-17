#include "libreria_gz.h"
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>


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

if (argc!=3)
{
  printf("Usage %s fbin_file seq_name\n\n",argv[0]);
  return -1;
}
 


res=read_seq(argv[1],argv[2], &fasta, &qual, &extras);

if (res==0){
	
	printf(">%s\n%s\n", argv[2], fasta);
	printf("%s\n",qual);
	if (extras!=NULL) printf ("extras:%s\n",extras);
}

//res=read_seq("prueba2gz.fbin","F143CJN01EN6AH", &fasta, &qual, &extras);
if ( fasta!=NULL ) {free(fasta);fasta=NULL;}
if ( qual!=NULL ) {free(qual);qual=NULL;}
if ( extras!=NULL ) {free(extras);extras=NULL;}

return res;
}

