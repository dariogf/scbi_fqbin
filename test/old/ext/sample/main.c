#include "libreria_gz.h"
#include <stdio.h>
#include <ctype.h>


/*******************************************************/
/* main                                                */
/*******************************************************/
int main(int argc, char *argv[])
{

char *fasta=NULL;
char *qual=NULL;
char *extras=NULL;
int size=5000;
            
if (argc!=4)
{
  printf("Usage %s fasta_file qual_file output_file\n\n",argv[0]);
  exit(-1);
}
 


// prueba de lectura
printf ("Comienzo\n");

//int res=read_seq("data/F143CJN01.fbin","F143CJN01D96I9", &fasta, &qual, &extras);
//int res=read_seq("fasta_filt_w_dict.fbin","F143CJN01D96I9", &fasta, &qual, &extras);

// prueba de escritura
/*int res=write_seq("borrar","S143CJN01D96I9", "primero","QUANTAQUANTAQUANTAQUANTAQUANTA","");
    res=write_seq("borrar","3CJN01D96I9", "el de enmedio","QUANTAQUANTAQUANTAQUANTAQUANTA","");
    res=write_seq("borrar","3CJsdasdaN01D96I9", "ultimo","QUANTAQUANTAQUANTAQUANTAQUANTA","");
*/
//init_dicts(dict_fasta,dict_qual,32767);
//int res=process_biofile("data/F143CJN01.fasta","/tmp/prueba2gz.fbin");
int res=process_biofile(argv[1],argv[3]);

//res=read_seq("fasta_filt_w_dict.fbin","F143CJN01D96I9", &fasta, &qual, &extras);
/*
int res=read_seq("prueba2.fbin","F143CJN01DI5MZ", &fasta, &qual, &extras);
printf ("-------------------------------------------------------------\n");
printf ("RES of read_seq1 call is :%d\n",res);
if ( res==0 ) printf ("fasta:%s\n size:%d\n",fasta,sizeof(fasta));
if ( res==0 ) printf ("qual:%s\n",qual);
if (( res==0 )&& (extras!=NULL)) printf ("extras:%s\n",extras);
if ( fasta!=NULL ) {free(fasta);fasta=NULL;}
if ( qual!=NULL ) {free(qual);qual=NULL;}
if ( extras!=NULL ) {free(extras);extras=NULL;}
*/
//int res=read_seq("/tmp/prueba2gz.fbin","F143CJN01BO14N", &fasta, &qual, &extras);
//int res=read_seq("/tmp/prueba2gz.fbin","F143CJN01D2X26", &fasta, &qual, &extras);
//res=read_seq("prueba2gz.fbin","F143CJN01EBIJN", &fasta, &qual, &extras);
res=read_seq(argv[3],"F143CJN01DZW7L", &fasta, &qual, &extras);
//res=read_seq("prueba2gz.fbin","F143CJN01EN6AH", &fasta, &qual, &extras);
printf ("-------------------------------------------------------------\n");
printf ("RES of read_seq2 call is :%d\n",res);
if ( res==0 ) printf ("fasta:%s\n fasta size:%d\n",fasta,strlen(fasta));
if ( res==0 ) printf ("qual:%s\n",qual);
if (( res==0 )&& (extras!=NULL)) printf ("extras:%s\n",extras);
if ( fasta!=NULL ) {free(fasta);fasta=NULL;}
if ( qual!=NULL ) {free(qual);qual=NULL;}
if ( extras!=NULL ) {free(extras);extras=NULL;}

printf ("***************************\n");
printf ("Sequential reads\n");

initialize_sequential_reads(argv[3]);
char *sname=NULL;
while (read_data_sequential(&sname, &fasta, &qual, &extras)==0)
{
	printf ("***************************\n");
	printf ("RES of read_seq2 call is :%d, sname:%s\n",res,sname);
	if ( res==0 ) printf ("fasta:%s fasta size:%d\n",fasta,strlen(fasta));
	if ( res==0 ) printf ("qual:%s",qual);
	if (( res==0 )&& (extras!=NULL)) printf ("extras:%s\n",extras);
	if ( fasta!=NULL ) {free(fasta);fasta=NULL;}
	if ( qual!=NULL ) {free(qual);qual=NULL;}
	if ( extras!=NULL ) {free(extras);extras=NULL;}
}
close_sequential_reads();

return res;
}

