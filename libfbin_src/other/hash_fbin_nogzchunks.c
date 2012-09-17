
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

// creates a hash from an index file with the desired chunk size. Chunk size can be 
// adjusted to fit a good compromise between access speed and used space.
int hash_index_file(char *filename, int chunk_size, int skip_sort)
{

  char hash_file_name[MAXFNAME];
  char indexname[MAXFNAME];
  
  char sname[MAXSEQNAME];// sequence name
  long long beginH, gz_chunk=0;
  char tmp[SEQ_METADATA];
  int res=0;
  int error;

  // to save min, max sequences and current chunk
  char min_name[MAXSEQNAME];
  char max_name[MAXSEQNAME];
  long long current_chunk=0;
  long long count=0;
  

  strcpy(min_name,"");
  strcpy(max_name,"");

	// calc index and hash name
  snprintf(indexname,MAXFNAME,"%s.index",filename);
  snprintf(hash_file_name,MAXFNAME,"%s.index.hash",filename);
  
  // sort index file by external  command
  if(skip_sort==0)
  {
      char cmd[10000];
      snprintf(cmd,10000,"sort_index.sh %s",indexname);
      system(cmd);
  }

  // use sorted index
  // snprintf(indexname,MAXFNAME,"%s.index.sort",filename);
  
  // open index and hash file
  gzFile gzhash_file=gzopen(hash_file_name,"wb");
  gzFile gzfile_index=gzopen(indexname,"r");
  
  if (gzfile_index==NULL) {
	fprintf(stderr,"error opening gzfile_index :%s\n",gzerror(gzfile_index,&error));
	return -2;
  }
  
  if (gzhash_file==NULL) {
  	fprintf(stderr,"error opening gzhash_file :%s\n",gzerror(gzhash_file,&error));
  	return -2;
  }
  
  // repeat until EOF
  while ( gzgets(gzfile_index,tmp,sizeof(tmp))!=Z_NULL ) {
    
    // parse string
    sscanf(tmp,"%s %lld %lld",sname,&gz_chunk,&beginH);

    if(strcmp(sname,"UMACOMPRESSEDFORMAT")!=0) // valid index line
    {
    
        // clear chunk_data if any
        // if (gz_chunk!=current_chunk){
        if ((count%chunk_size)==0){
            if (strcmp(min_name,"")!=0){
                // there are data to write
                res=gzprintf(gzhash_file,"%s %s %lld\n",min_name,max_name,current_chunk);
            }
            strcpy(min_name,"");
            strcpy(max_name,"");
            // current_chunk=gz_chunk;
            current_chunk = gztell(gzfile_index);
        }

        // save min_name
    	if((strcmp(min_name,"")==0) || (strcmp(sname,min_name)<0))
    	{
            // replace min_name
            strcpy(min_name,sname);
    	}
	
    	//save max_name
    	if((strcmp(max_name,"")==0) || (strcmp(sname,max_name)>0))
    	{
            strcpy(max_name,sname);
      }
        
      count++;
    }
        
  }
  
  if (strcmp(min_name,"")!=0){
      // there are data to write
      res=gzprintf(gzhash_file,"%s %s %lld\n",min_name,max_name,current_chunk);
  }
  
  // close files
  gzclose(gzhash_file);
  gzclose(gzfile_index);
  
  return 0;
}


/*******************************************************/
/* main                                                */
/*******************************************************/
int main(int argc, char *argv[])
{
  // check params
  if (argc<2)
  {
    printf("Usage %s fbin_file [chunk_size [--skip_sort]]\n\n",argv[0]);
    return -1;
  }
  
  int chunk_size=10000;
  int skip_sort=0;
  
  if(argc==3){
      chunk_size=atoi(argv[2]);
  }
  
  if (argc==4){
      skip_sort=1;
  }
  
  int res=hash_index_file(argv[1],chunk_size, skip_sort);
  
  return res;
}

