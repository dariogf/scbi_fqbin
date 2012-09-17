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

//gzFile gzf_bin;
// struct file_data filed;

struct file_data *filed=NULL;

if (argc!=2)
{
  printf("Usage %s fbin_file\n\n",argv[0]);
  return -1;
}

initialize_sequential_reads(&filed, argv[1]);

char *sname=NULL;

while ((res=read_data_sequential(filed, &sname, &fasta, &qual, &extras))==0)
{
	// printf("res:%d\n",res);
	if (res==0){
		
		printf(">%s\n%s\n", sname, fasta);
		printf("%s\n",qual);
		if (extras!=NULL) printf ("extras:%s\n",extras);
	}
	
	if ( fasta!=NULL ) {free(fasta);fasta=NULL;}
	if ( qual!=NULL ) {free(qual);qual=NULL;}
	if ( extras!=NULL ) {free(extras);extras=NULL;}
}

close_sequential_reads(filed);

return res;
}

