
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


long long last_chunk_file(char *filename){
  
  // open file again to annotate chunk
  int file=open(filename,O_APPEND);

  //goto end of file
  long long pos=lseek(file,0,SEEK_END);
  if (pos==-1)  {fprintf(stderr,"error %d seeking file %s :%s\n",errno,filename,strerror(errno));return -1;}

  close(file);
  
  return pos;
}

// creates a hash from an index file with the desired chunk size. Chunk size can be 
// adjusted to fit a good compromise between access speed and used space.
int hash_index_file(char *filename, int chunk_size, int skip_sort)
{

  char hash_file_name[MAXFNAME];
  char indexname[MAXFNAME];
  char sorted_indexname[MAXFNAME];
  
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
      snprintf(cmd,10000,"sort_index %s",indexname);
      system(cmd);
  }

  // use sorted index
  snprintf(sorted_indexname,MAXFNAME,"%s.index.sorted",filename);
  
  // open hash file
  gzFile gzhash_file=gzopen(hash_file_name,"wb");

  // open sorted index file
  gzFile gzsorted_file_index=gzopen(sorted_indexname,"r");
  
  // open output index file
  // int file_index=open(indexname,flags,0644);
  gzFile gzfile_index=gzopen(indexname,"w");
  
  // write header
  gzprintf(gzfile_index,"UMACOMPRESSEDFORMAT 1 0 0 999999999999 999999999999\n");

  //reopen
  gzclose(gzfile_index);
  gzfile_index=gzopen(indexname,"ab");


  if (gzsorted_file_index==NULL) {
	fprintf(stderr,"error opening gzsorted_file_index :%s\n",gzerror(gzsorted_file_index,&error));
	return -2;
  }
  
  if (gzfile_index==NULL) {
	fprintf(stderr,"error opening gzfile_index :%s\n",gzerror(gzfile_index,&error));
	return -2;
  }
  
  if (gzhash_file==NULL) {
  	fprintf(stderr,"error opening gzhash_file :%s\n",gzerror(gzhash_file,&error));
  	return -2;
  }
  
  // podria leerse saltando linea 1, y luego leyendo 10000 lineas sin sscanf
  
  // repeat until EOF
  while ( gzgets(gzsorted_file_index,tmp,sizeof(tmp))!=Z_NULL ) {
    
    // parse string
    sscanf(tmp,"%s %lld %lld",sname,&gz_chunk,&beginH);

    if(strcmp(sname,"UMACOMPRESSEDFORMAT")!=0) // valid index line
    {
        
        // clear chunk_data if any
        // if (gz_chunk!=current_chunk){
        if ((count%chunk_size)==0){
            if (strcmp(min_name,"")!=0){
                // there is data to write
                res=gzprintf(gzhash_file,"%s %s %lld\n",min_name,max_name,current_chunk);
            }
            
            strcpy(min_name,"");
            strcpy(max_name,"");
            // current_chunk=gz_chunk;
            current_chunk = gztell(gzfile_index);

            //reopen new gzchunk
            gzclose(gzfile_index);
            current_chunk = last_chunk_file(indexname);
            gzfile_index=gzopen(indexname,"ab");

        }
         
        // write line to current gzchunk in index
        gzprintf(gzfile_index,tmp);
        
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
  gzclose(gzsorted_file_index);
  
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
    printf("Usage %s fqbin_file [chunk_size [--skip_sort]]\n\n",argv[0]);
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

