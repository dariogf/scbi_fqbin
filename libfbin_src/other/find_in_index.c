
#include <stdio.h>
#include <string.h>
#include <time.h>


#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <errno.h>

#include <zlib.h>
#include <zlib.h>
#include <stdlib.h>

// Maximum file name (including .idx)
#define MAXFNAME 512

// Maximum lenght of the name of a sequence
#define MAXSEQNAME 1024
#define MAXSEQLENGTH 150000000
#define DEBUG 1
#define FALSE 0 
#define TRUE 1 

#define INVALID_FASTQ_FORMAT -5
#define INVALID_FASTA_FORMAT -6

#define SEQ_METADATA 10000


// int mystrcmp(const char *a,const char *b)
// {
//   return strlen(a)-strlen(b)?strlen(a)-strlen(b):strcmp(a,b);
// }


int first_line(gzFile file, line){
    
    
    
}



long long find_seq_in_hash(char *filename,char *sname)
{

  char file_name[MAXFNAME];
  // char indexname[MAXFNAME];
  int error;
  char sname1[MAXSEQNAME];// sequence name
  char sname2[MAXSEQNAME];// sequence name
  long long gz_chunk=0;
  char tmp[SEQ_METADATA];
  long long res=-1;

  // to save min, max sequences and current chunk
  char min_name[MAXSEQNAME];
  char max_name[MAXSEQNAME];
  long long current_chunk=0;
  

  strcpy(min_name,"");
  strcpy(max_name,"");

  // calc index and hash name
  // snprintf(indexname,MAXFNAME,"%s.index",filename);
  snprintf(file_name,MAXFNAME,"%s.index",filename);
  
  // open index and hash file
  gzFile gzhash_file=gzopen(file_name,"r");
  
  if (gzhash_file==NULL) {
  	fprintf(stderr,"error opening gzhash_file :%s\n",gzerror(gzhash_file,&error));
  	return -2;
  }
  
  // repeat until EOF
  while ( gzgets(gzhash_file,tmp,sizeof(tmp))!=Z_NULL ) {
    
      // printf("%s\n",tmp);
    // parse string
    int reads=sscanf(tmp,"%s %s %lld",sname1,sname2,&gz_chunk);

    if(reads==3) // valid index line
    {
        // 
        if((mystrcmp(sname,sname1)>=0) && (mystrcmp(sname,sname2)<=0))
        {
            printf("%s in [%s,%s]\n",sname,sname1,sname2);
            res = gz_chunk;
            break;
        }else{
            printf("%s NOT IN [%s,%s]\n",sname,sname1,sname2);
        }
        
    }
        
  }
  
  // close files
  gzclose(gzhash_file);
  
  return res;
}


/*******************************************************/
/* main                                                */
/*******************************************************/
int main(int argc, char *argv[])
{
  // check params
  if (argc!=3)
  {
    printf("Usage %s fbin_index_file seq_name\n\n",argv[0]);
    return -1;
  }
  
  int c1=mystrcmp("SRR314795.1","SRR314795.1000000");
  int c2=mystrcmp("SRR314795.1000000","SRR314795.9");
  printf("RES: %d,%d\n",c1,c2);
  
  long long chunk=find_seq_in_hash(argv[1],argv[2]);
  
  printf("Chunk: %lld\n",chunk);
  
}

