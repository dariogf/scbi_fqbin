#include "lib_fqbin.h"


int check_error(int error_condition, char *message, int return_value){
    if (error_condition) {
        fprintf(stderr,"Error %d; %s\nMSG:%s\n",errno ,message, strerror(errno)); 
        return return_value;
    }
}


int write_seq(struct file_data *file, char *seq_name, char *fasta, char *qual, char *extras)
{
	// compress data
	char metainfo[SEQ_METADATA];
  	int error=0;

	// Convert  to fastq

    // char qual[150000];
    // sprintf(qual,"");
    // char *sep = " ";
    // char *word, *brkt;
    // 
    // for (word = strtok_r(in_qual, sep, &brkt);
    //      word;
    //      word = strtok_r(NULL, sep, &brkt))
    // {
    //  sprintf(qual,"%s%c",qual,atoi(word)+33);
    //     // strcat(qual2,".");
    //     // printf("%s,%s,%c\n",word,qual2, atoi(word)+33);
    // 
    // }

    // printf("write seq\n");

	
  if (file->gzf_bin==NULL) {fprintf(stderr,"error with gzfile_bin, is NULL :%s\n",gzerror(file->gzf_bin,&error));return -2;}
  
  long fasta_len=strlen(fasta);

  // preprocess fasta string
  // if (fasta_len>0)
  // {
  //     // printf("IN:%s, %ld\n",fasta,strlen(fasta));
  //     // BWXform(fasta, fasta,1);
  //     // printf("\nOUT:\n");
  //     // write(1,fasta,fasta_len);
  //     // printf("\n");
  //     // fflush(1);
  // }
  
  
  // preprocess quality string
  	if (strlen(qual)>0){
	    
      // printf("========================\nQ:%s\n",qual);
      // process qual
      int i;
      
      // annotate first qual to compare with rest of quals
      char old=qual[0];
      
      // if qual is going to be flattened, take the flatten limit as old
      if ((file->flatten_qual>0) & (old>=file->flatten_qual))
      {
          old=file->flatten_qual;
      }
      
      int same_qual=1;
      
      // process qual string
      for( i = 0; i < strlen(qual); ++i)
      {

          //discretize up by discretize_qual value
          if (file->discretize_qual>1){
             // qual[i]=qual[i]-(qual[i] % 2)+2-1;
             qual[i]=(qual[i] / file->discretize_qual)*file->discretize_qual;
          }
          
          // printf("FL:%c>=%c\n",qual[i],file->flatten_qual);
          // trim high qualitys
          if ((file->flatten_qual>0) & (qual[i]>file->flatten_qual)){
              qual[i]=file->flatten_qual;
          }
          
          if (qual[i]!=old) {same_qual=0;}
      }
      
      if (same_qual)
      {
          // trim qual string
          sprintf(qual,"%c",old);
          qual[1]=0;
          // printf("EQUAL: %s,%s,%ld\n",seq_name, qual,strlen(qual));
      }
      
      // printf("\nQ:%s\n",qual);
	    
	}
    // printf("Calc fasta_len\n");
    // printf("Calculated fasta_len\n");

    snprintf(metainfo,SEQ_METADATA-1,"9999%s %ld %ld %ld", seq_name, fasta_len, strlen(qual), strlen(extras));
    snprintf(metainfo,SEQ_METADATA-1,"%4ld%s %ld %ld %ld", strlen(metainfo)-4, seq_name, fasta_len, strlen(qual), strlen(extras));
    // snprintf(metainfo,SEQ_METADATA-1,"9999%s %ld %ld", seq_name, strlen(fasta), strlen(extras));
    // snprintf(metainfo,SEQ_METADATA-1,"%4ld%s %ld %ld", strlen(metainfo)-4, seq_name, strlen(fasta), strlen(extras));

	// get begin pos of header
	long beginH=gztell(file->gzf_bin);

	// TODO check gztell
	if (beginH==-1) {fprintf(stderr,"error with pos of beginH of gzfile_bin :%s\n", gzerror(file->gzf_bin,&error)); return -2;}

	// write seq to bin file
    gzwrite(file->gzf_bin, metainfo, strlen(metainfo));

	// TODO check gzwrite
	long beginI=gztell(file->gzf_bin);
  
	if (beginI==-1) {fprintf(stderr,"error with pos of beginI of gzfile :%s\n",gzerror(file->gzf_bin,&error));return -2;}
	
    // printf("Before write fasta\n");
    
	int res=1;
	if (fasta_len>0) res=gzwrite(file->gzf_bin,fasta,fasta_len); //Z_FILTERED);
	
	if ( res==0 ) { fprintf(stderr,"Error when writting fasta\n");return -8;}
	long fastaS=gztell(file->gzf_bin)-beginI;

    // printf("After write fasta\n");
	
	
	if (strlen(qual)>0){
	    res=gzwrite(file->gzf_bin,qual,strlen(qual)); //Z_FILTERED);
	} 
	
	if ( res==0 ) { fprintf(stderr,"Error when writting qual\n");return -8;}
	long qualS=gztell(file->gzf_bin)-fastaS-beginI;
	
	if (strlen(extras)>0) res=gzwrite(file->gzf_bin,extras,strlen(extras)); //Z_FILTERED);
	
	if ( res==0 ) { fprintf(stderr,"Error when writting extras\n");return -8;}
	long extrasS=gztell(file->gzf_bin)-qualS-fastaS-beginI;
	
	//	add_sequence(&seql,seq_name,pos_chunk_gz,beginI,fastaS,qualS,extrasS);

	// Write index file
	if((file)->create_index)
	{
        char tmp[SEQ_METADATA];
        sprintf(tmp,"%s %lld %ld\n",seq_name,file->pos_chunk_gz,beginH);
        gzwrite(file->gzf_index,tmp,strlen(tmp));
    }

	(file->counter)++;
    //	if (counter > 2) fprintf(stderr,"Probando static counter para llamadas desde ruby, valor %d\n",counter);

  // create new chunk
  if (((file->counter)%10000)==0) {
        //      curr_time=time(NULL);
        // printf("10k seqs in:%f secs\n",difftime(curr_time,prev_time));
        // prev_time=curr_time;

		// close current chunk
	  	gzclose(file->gzf_bin);
		
		// open file again to annotate chunk
		int file_bin=open(file->name,O_APPEND);
		
		//goto end of file
		long long pos=lseek(file_bin,0,SEEK_END);
		if (pos==-1)  {fprintf(stderr,"error %d seeking file :%s\n",errno,strerror(errno));return -1;}
		
		// annotate chunk pos
		file->pos_chunk_gz=pos;
		
		close(file_bin);
		
		// open new gzfile 
		file->gzf_bin=gzopen(file->name,"ab");
		if (file->gzf_bin==NULL) {fprintf(stderr,"error opening gzfile :%s\n",gzerror(file->gzf_bin,&error));return -2;}
	}
	
	return 0;
}



/* Reads the metadata from the main file 
   It initializes the version variable 
*/
int read_bin_file_metadata(struct file_data *filed)
{
	char header[SEQ_METADATA];
	int fastaS,qualS,extrasS=0;
	int ver,subver;

    // printf("pos1b %ld\n",gztell(filed->gzf_bin));//lseek ((*filed), 0, SEEK_CUR))

	int res=read_seq_header(filed->gzf_bin, header, &fastaS, &qualS, &extrasS);
    // printf("pos2b %ld\n",gztell(filed->gzf_bin));//lseek ((*filed), 0, SEEK_CUR))
	

	if ( res!=0 ) {return -1;}
	if ( strlen(header)<20 ) {fprintf(stderr,"Too short sequence header:%s. lenght:%ld\n",header,strlen(header));return -1;}

	// 28UMACOMPRESSEDFORMAT_1_0 0 0 0
//	header[strlen(header)-2]=0;

	if (strncmp(header,"UMACOMPRESSEDFORMAT_",19)!=0) {fprintf(stderr,"Incorrect header in file, header:%s\n",header);return -1;}
	// TODO fill the file_data structure with the header data
	if (sscanf(header,"UMACOMPRESSEDFORMAT_%d_%d",&ver,&subver)!=2) return -1;
	//if (sscanf(header,"UMACOMPRESSEDFORMAT_%d_%d",&(filed->version),&(filed->subversion))!=2) return -1;
	filed->version=11;//ver;
	filed->subversion=subver;
//        fprintf(stderr,"file version:%d,%d\n",filed->version,filed->subversion);
        return 0;
}

/* Reads the metadata from the index file 
   It initializes the version and binary_search variable 
*/
int read_index_file_metadata(struct file_data *filed)
{
	char header[SEQ_METADATA];
	int fastaS,qualS,extrasS=0;


	int res=read_seq_header(filed->gzf_bin, header, &fastaS, &qualS, &extrasS);

	if ( strlen(header)<19 ) {return -1;}

	// 28UMACOMPRESSEDFORMAT 1 0 0 0 0
	header[strlen(header)-2]=0;

	if (strncmp(header,"UMACOMPRESSEDFORMAT",19)!=0) return -1;
	// TODO fill the file_data structure with the header data
	if (sscanf(header,"UMACOMPRESSEDFORMAT %d %d",&(filed->version),&(filed->subversion))!=2) {
		fprintf(stderr,"SEQ READ:Header incorrect when reading versions:%s.\n",header);
		return -1;
	}

        return 0;
}

/* reads the header of a sequence in the main file.
   the pointer to the file points to the fasta data after calling read_seq_header 
   returns 0 if ok
 	-1 if there is an error
	-2 if EOF
*/

int read_seq_header(gzFile *gzf_bin, char *seq_name,int *fastaS, int *qualS, int *extrasS)
{
  int header_size=4;
  char hsize[40];
  char tmp[1000];
  char sname[SEQ_METADATA];

  long pos=gzread(gzf_bin,hsize,header_size);
               
  // EOF found
  if ( pos==0 ) return -2;

  // Error reading file
  if ( pos==-1 ) {fprintf(stderr,"Incorrect sequence header. File may be corrupted\n");return -1;}

  hsize[pos]=0;
  sscanf(hsize,"%d",&header_size);
  pos=gzread(gzf_bin,tmp,header_size);

  if ( pos==0 ) return -2;
	
  if ( pos==-1 ) {fprintf(stderr,"Incorrect sequence header. File may be corrupted\n");return -1;}
	
  tmp[header_size]=0;
  int reads=sscanf(tmp,"%s %d %d %d",sname,fastaS,qualS,extrasS);
	
  if (reads!=4) {return -1;};
  	
  if (seq_name!=NULL) strncpy(seq_name,sname,SEQ_METADATA);

  return 0;
}

// check files before reading
// it initializes the previous variables, file_version and binary_search
// result :
// 0 : if both the bin and index files exists and are from the current version
// 1 : if both the bin and index files exists but are from another version
// 2 : if both files are missing
// 3 : if bin file is missing
// 4 : if index file is missing
int check_files()
{
  
// open the files, read and check the header
  return 0;
}

// returns the version of the opened file
int version(struct file_data *filed)
{
  if (filed->gzf_bin==NULL) return -1;
  return filed->version;
}

// returns the version of the opened file
int subversion(struct file_data *filed)
{
  if (filed->gzf_bin==NULL) return -1;
  return filed->subversion;
}

/* 
mode can be:
1 - random, for each read it begins to read from the beggining of index
2 - sequential, it keeps the position inside the index and main files.
*/
int initialize_sequential_reads(struct file_data ** filed, char *filename)
{
	char header[SEQ_METADATA];
	int fastaS,qualS,extrasS=0;
	
	if ( *filed  == NULL ) {*filed=malloc(sizeof(struct file_data));}
	
	
	(*filed)->gzf_bin=gzopen(filename,"r");
	
	int res=check_error((*filed)->gzf_bin==NULL,"Unable to open file",-1);
	
	strncpy((*filed)->name,filename,MAXFNAME);
	(*filed)->error=0;
			


	// reads the metadata
/*
	int res=read_seq_header(filed->gzf_bin, header, &fastaS, &qualS, &extrasS);

	if ( strlen(header)<19 ) {fprintf(stderr,"SEQ READ:Header incorrect:%s.\n",header);return -1;}

	// 28UMACOMPRESSEDFORMAT_1 0 0 0
	header[strlen(header)-2]=0;

	if (strncmp(header,"UMACOMPRESSEDFORMAT",19)!=0) return -1;
	// TODO fill the file_data structure with the header data
*/
    // printf("pos1 %ld\n",gztell((*filed)->gzf_bin));
	res= read_bin_file_metadata(*filed);
    // printf("pos2 %ld\n",gztell((*filed)->gzf_bin));
    
  // inspect_file_data_struct(filed);

	return res;
}


int read_data_sequential(struct file_data *filed,char **seq_name, char **fasta, char **qual, char **extras)
{
  int res=0;
  int error=0;
  int fastaS,qualS,extrasS=0;

  if ( *seq_name == NULL ) {*seq_name=(char *)malloc(SEQ_METADATA);strncpy(*seq_name,"",4);}

  res=read_seq_header(filed->gzf_bin, *seq_name, &fastaS, &qualS, &extrasS);
  // printf("FS:%d,QS:%d\n",fastaS,qualS);
  if (res==-2) // EOF
		return -9;
		
  if ( *fasta  == NULL ) {*fasta=(char *)malloc(fastaS+1);strncpy(*fasta,"",fastaS);}
  if ( *qual == NULL) {*qual=(char *)malloc(qualS+1);strncpy(*qual,"",qualS);}
  if (( *extras == NULL )&&(extrasS>0)) {*extras=(char *)malloc(extrasS+1);strncpy(*extras,"",extrasS);}

  long pos=gzread(filed->gzf_bin,*fasta,fastaS);
  
  // BWReverseXform(*fasta, *fasta, 1, fastaS);
  
  (*fasta)[fastaS]=0;
  pos=gzread(filed->gzf_bin,*qual,qualS);

  // if only one qual read, repeat it
  if((qualS==1) & (qualS!=fastaS))
  {
      char q=*qual[0];
      free(*qual);
      *qual=NULL;
      qualS=fastaS;
      if ( *qual == NULL) {*qual=(char *)malloc(qualS+1);strncpy(*qual,"",qualS);}
      memset(*qual,q,qualS);
  }
  
  (*qual)[qualS]=0;
  // printf("LLEGA\n:");
  if (extrasS>0) {pos=gzread(filed->gzf_bin,*extras,extrasS);(*extras)[extrasS]=0;}
  return 0;
  
}
int close_sequential_reads(struct file_data *file_d)
{
  gzclose(file_d->gzf_bin); 
  if ( (file_d)!= NULL ) {  free(file_d);}
}

int regenerate_index(char * filename){


/*    NO SE PUEDE USAR ESTA FUNCION PORQUE CUANDO SE LEE EL FICHERO, LA POSICION DE LECTURA NO TIENE PORQUE COINCIDIR CON LA DE ESCRITURA CUANDO SE CREÓ EL FBIN, DE MODO QUE AL GENERAR EL INDICE NO SABEMOS DE QUE CHUNK LEER. LA LIB ZLIB LEE CON EL GZREAD UN TROZO DE FICHERO DEL TAMAÑO QUE ELLA DECIDA, Y ASI NO COINCIDE LUEGO CON LOS BLOQUES QUE NOSOTROS ESCRIBIAMOS
  */  
  
    return -1;
    
    int res=0;
    
    long long bcount=0;
    long long pos=0;
    
    struct file_data *filed=NULL;
    int fastaS,qualS,extrasS=0;
    
    char header[SEQ_METADATA];
    char metainfo[SEQ_METADATA];
    
    char *seq_name = NULL;//[SEQ_METADATA];
    
    if ( (filed)  == NULL ) {filed=malloc(sizeof(struct file_data));}
    
    filed->pos_chunk_gz=0;
    filed->counter=0;

    // memset(metainfo,'a',20);
    // printf("%s\n",metainfo);
    // open file
    
    int file_bin=open(filename,O_RDONLY);
    
    // filed->gzf_bin=gzopen(filename,"rb");
    
    filed->gzf_bin = gzdopen(file_bin,"rb");
    
	res=check_error(filed->gzf_bin==NULL,"Unable to open file",-1);
    
	filed->error=0;
	strncpy(filed->name,filename,MAXFNAME);
    
    // read header
	res=read_bin_file_metadata(filed);
	
	
	int error=0;
    int seek_res=0;


    if ( seq_name == NULL ) {seq_name=(char *)malloc(SEQ_METADATA);strncpy(seq_name,"",4);}
    char basura[100000];
    
    // char ** basura;
    //     if ( *basura  == NULL ) {*basura=(char *)malloc(fastaS+qualS+extrasS+1); strncpy(*basura,"",fastaS+qualS+extrasS);}
    //     
    while (res==0){
        long beginH=gztell(filed->gzf_bin);

        // long long pos=lseek(file_bin,0,SEEK_CUR);
        // printf("BEF: %lld\n",pos);
    	
        res=read_seq_header(filed->gzf_bin, seq_name, &fastaS, &qualS, &extrasS);
        if (res==-2) // EOF
      		return -9;
      		
        bcount=bcount+4+fastaS+qualS+extrasS;
  		
        // printf("SEQ: %s, skip: %d, res:%d\n",seq_name,fastaS+qualS+extrasS,res);
        
        snprintf(metainfo,SEQ_METADATA-1,"9999%s %d %d %d", seq_name, fastaS, qualS, extrasS);
        // snprintf(metainfo,SEQ_METADATA-1,"%4ld%s %d %d %d", strlen(metainfo)-4, seq_name, fastaS,qualS,extrasS);
        
        // Write index file
    	char tmp[SEQ_METADATA];
    	sprintf(tmp,"%s %lld %ld\n",seq_name,filed->pos_chunk_gz,beginH);
        
        printf("%s %lld %ld\n",seq_name,filed->pos_chunk_gz,beginH);
        printf("%s\n",metainfo);
    	
        pos=lseek(file_bin,0,SEEK_CUR);
            // long long pos2=gztell(filed->gzf_bin);

            printf("Antes seek: %lld\n",pos);
        
        // printf("%s\n",seq_name);

        seek_res=gzseek(filed->gzf_bin,fastaS+qualS+extrasS,SEEK_CUR);
        
        printf("bcount:%lld\n",bcount);
        
        // long long pos4=lseek(file_bin,0,SEEK_CUR);
		
		    
		    pos=lseek(file_bin,0,SEEK_CUR);
		    printf("Despues seek: %lld\n",pos);
        // printf("AFT: %lld\n=============\n",pos4);
        
        // long pos3=gzread(filed->gzf_bin,&basura,fastaS+qualS+extrasS);
        
        (filed->counter)++;

		long long pos2=gztell(filed->gzf_bin);
	  	


        //  new chunk
        if (((filed->counter)%10000)==0) {
            
            printf("SEQ 10K:%s\n",seq_name);
    		// close current chunk
            pos=lseek(file_bin,0,SEEK_CUR);
            // long long pos2=gztell(filed->gzf_bin);

            printf("FINAL BLOCK POSSSSSSSSSSSSSS: %lld\n",pos);
            
    	  	gzclose(filed->gzf_bin);
    	  	
    		// open file again to annotate chunk
    		file_bin=open(filed->name,O_RDONLY);

    		//goto end of file
    		pos=lseek(file_bin,pos,SEEK_SET);
    		if (pos==-1)  {fprintf(stderr,"error %d seeking file :%s\n",errno,strerror(errno));return -1;}

            printf("FINAL BLOCK POSSSSSSSSSSSSSS: %lld\n",pos);

    		// annotate chunk pos
    		filed->pos_chunk_gz=pos;

    		close(file_bin);

    		// open new gzfile 
    		filed->gzf_bin=gzdopen(file_bin,"rb");
    		if (filed->gzf_bin==NULL) {fprintf(stderr,"error opening gzfile :%s\n",gzerror(filed->gzf_bin,&error));return -2;}
    	}
        
        
        
    }

    // if ( *fasta  == NULL ) {*fasta=(char *)malloc(fastaS+1);strncpy(*fasta,"",fastaS);}
    // if ( *qual == NULL) {*qual=(char *)malloc(qualS+1);strncpy(*qual,"",qualS);}
    // if (( *extras == NULL )&&(extrasS>0)) {*extras=(char *)malloc(extrasS+1);strncpy(*extras,"",extrasS);}


    // long pos=gzread(filed->gzf_bin,*fasta,fastaS);
    // (*fasta)[fastaS]=0;
    // pos=gzread(filed->gzf_bin,*qual,qualS);
    // (*qual)[qualS]=0;
    // if (extrasS>0) {pos=gzread(filed->gzf_bin,*extras,extrasS);(*extras)[extrasS]=0;}
    
	
    // close files
    gzclose(filed->gzf_bin);
    if ( (filed)!= NULL ) {free(filed);}

    return 0;

    
}

long long find_seq_in_hash(char *filename,char *sname)
{

  char hash_file_name[MAXFNAME];
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
  snprintf(hash_file_name,MAXFNAME,"%s.index.hash",filename);
  
  
  // open index and hash file
  gzFile gzhash_file=gzopen(hash_file_name,"r");
  
  if (gzhash_file==NULL) {
    // fprintf(stderr,"error opening gzhash_file :%s\n",gzerror(gzhash_file,&error));
    // no hash file found
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
        if((strcmp(sname,sname1)>=0) && (strcmp(sname,sname2)<=0))
        {   
            #if DEBUG
                printf("%s in [%s,%s]\n",sname,sname1,sname2);
            #endif
            res = gz_chunk;
            break;
        }else{
            // printf("%s NOT IN [%s,%s]\n",sname,sname1,sname2);
        }
        
    }
        
  }
  
  // close files
  gzclose(gzhash_file);
  
  return res;
}


long long find_seq_in_index(char *filename,char *sname, long long index_chunk, long long *gz_chunk, long long *beginH){
    
    long long chunk=-1;

    char file_name[MAXFNAME];
    // char indexname[MAXFNAME];
    int error;
    char sname1[MAXSEQNAME];// sequence name
    char sname2[MAXSEQNAME];// sequence name
    long long aux_beginH=-1,aux_gz_chunk=-1;
    char tmp[SEQ_METADATA];
    long long res=-1;
    
    *gz_chunk = 0;
    *beginH=0;

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
    // gzFile filegz=gzopen(file_name,"r");
    int file=open(file_name, O_RDONLY);
        
    if (file<0)
    {
      return -2;
    }
    
    if(index_chunk>0)
    {
        #if DEBUG
            printf("Seek to %lld\n",index_chunk);
        #endif
        // res=gzseek(filegz,index_chunk,SEEK_SET);
        res=lseek(file,index_chunk,SEEK_SET);
        
        // printf("Seeked\n");
    }
    
    gzFile filegz=gzdopen(file,"r");
    
    if (filegz==NULL) {
        // fprintf(stderr,"error opening gzhash_file :%s\n",gzerror(filegz,&error));
    	return -2;
    }
    
    // repeat until EOF
    while ( gzgets(filegz,tmp,sizeof(tmp))!=Z_NULL ) {

      // printf("%s\n",tmp);
      // parse string
      // int reads=sscanf(tmp,"%s %s %lld",sname1,sname2,&gz_chunk);
      int reads=sscanf(tmp,"%s %lld %lld",sname1,&aux_gz_chunk,&aux_beginH);
      

      if(reads==3) // valid index line
      {
          // 
          if(strcmp(sname,sname1)==0)
          {
              #if DEBUG
                  printf("%s IN %s\n",sname, tmp);
              #endif

              chunk = aux_gz_chunk;
              *gz_chunk = aux_gz_chunk;
              *beginH=aux_beginH;
              
              // beginH=gz_beginH;
              break;
          }else{
              #if DEBUG
               printf("NOT IN %s",tmp);
              #endif
               // break;
          }

      }

    }

    // close files
    gzclose(filegz);
    
    return chunk;
}



/* 
   read_seq reads from filename the sequence named seq_name and returns its
	fasta, qual and extras in those variables.
   It returns 0 if there are no errors, otherwise it returns:
   -2 : error opening index file (it doesn't exists)
   -3 : error reading index file
   -4 : error sequence not found in index file
   -5 : error opening file (it doesn't exists)
   -6 : error reading file 
   -7 : error sequence not found
   -8 : error uncompressing sequence
   -9 : EOF

*/

int read_seq(char *filename, char *seq_name, char **fasta, char **qual, char **extras)
{
/* Hacer grep en filename.index de seq_name */
/* Una vez encontrado leer su info (indice y offsets) */
/* leer de filename en sus offests el fasta qual y extras */
/* Descomprimirlo y devolverlo */

  char indexname[MAXFNAME];
  char sname[MAXSEQNAME];// sequence name
  // char *fasta_comp;      // compressed fasta
  // char *qual_comp;     // compressed qual
  // char *extras_comp;     // compressed extras
  long long beginH, gz_chunk=0;
  int fastaS, qualS, extrasS=0;
  char tmp[SEQ_METADATA];
  int res=0;
  int error=0;

  // int bufsize=MAXSEQLENGTH;
  // 
  //    // allocate memory for return data if necessary
  //   if ( *fasta  == NULL ) {*fasta=(char *)malloc(bufsize);strncpy(*fasta,"",bufsize);}
  //   if ( *qual == NULL) {*qual=(char *)malloc(bufsize);strncpy(*qual,"",bufsize);}
  //   if ( *extras == NULL ) {*extras=(char *)malloc(bufsize);strncpy(*extras,"",bufsize);}
  
	// calc index name
  // snprintf(indexname,MAXFNAME,"%s.index",filename);

  long long chunk=find_seq_in_hash(filename,seq_name);
  // printf("Chunk: %lld\n",chunk);
   
  if (chunk<0){
      chunk=0;
  }
  
  if ((res=find_seq_in_index(filename,seq_name,chunk,&gz_chunk,&beginH))<0){
    return res;
  };


	// open index file
    //   gzFile gzfile_index=gzopen(indexname,"r");
    //   if (gzfile_index==NULL) {
    // fprintf(stderr,"error opening gzfile_index :%s\n",gzerror(gzfile_index,&error));
    // return -2;
    //   }
    // 
    //   // Reads the index to this info, and the offset to its data
    //   int reads=3;
    //   while ( reads == 3 ) {
    //          
    //          // read a chunk of data from index with the size of tmp
    //         gzgets(gzfile_index,tmp,sizeof(tmp));
    //          reads=sscanf(tmp,"%s %lld %lld",sname,&gz_chunk,&beginH);
    //          
    //          
    // 
    //          
    //          if (( reads != 3 ) && ( reads!=EOF )) {
    //              fprintf(stderr,"Error scanning index: %d\n",reads);
    //              gzclose(gzfile_index);
    //              return -3;
    //          }
    //          
    //          // sequence was finally found, exit loop
    //         if ( strncmp(sname, seq_name,MAXSEQNAME)==0) reads=999; // to get out, seq found
    //   }
    // 
    // // close index file
    //   gzclose(gzfile_index);
    // 
	// maybe sequence was not found
	// fprintf(stderr,"Sequence not found\n");
  // if (reads==EOF) {return -4;}

	// We get here if sequence was found
	
	#if DEBUG
	    printf("Index found %lld. Seeking\n",gz_chunk);
	#endif
	// open bin file to extract data
  int dataf=open(filename, O_RDONLY);

	// seek to chunk pos
	// TODO- ¿como se salta el chunk?
  // res=lseek(dataf,gz_chunk,SEEK_SET);
  res=lseek(dataf,gz_chunk,SEEK_SET);

	// TODO check res
  gzFile gzfile_bin=gzdopen(dataf,"r");
  
	// seek to seq inside chunk
  res=gzseek(gzfile_bin,beginH,SEEK_SET);
	// TODO check res

    // printf("Seeked\n");

  // read sequence header 
  res=read_seq_header(gzfile_bin,NULL, &fastaS, &qualS, &extrasS);
  
  // int bufsize=MAXSEQLENGTH;

  // memset(*qual,q,qualS);
  // allocate memory for return data if necessary
  if ( *fasta  == NULL ) {*fasta=(char *)malloc(fastaS+1);(*fasta)[0]=0;}
  if ( *qual == NULL) {*qual=(char *)malloc(qualS+1);(*qual)[0]=0;}
  if ( *extras == NULL ) {*extras=(char *)malloc(extrasS+1);(*extras)[0]=0;}
  
  long pos=gzread(gzfile_bin,*fasta,fastaS);
  // printf("LEIDO:%ld\n",pos);
  // BWReverseXform(*fasta, *fasta, 1, fastaS);

  (*fasta)[fastaS]=0;
  pos=gzread(gzfile_bin,*qual,qualS);
  
  // if only one qual read, repeat it
  if((qualS==1) & (qualS!=fastaS))
  {
      char q=*qual[0];
      free(*qual);
      *qual=NULL;
      qualS=fastaS;
      if ( *qual == NULL) {*qual=(char *)malloc(qualS+1);strncpy(*qual,"",qualS);}
      memset(*qual,q,qualS);
  }
  
  (*qual)[qualS]=0;

  if (extrasS>0) {pos=gzread(gzfile_bin,*extras,extrasS); (*extras)[extrasS]=0;}
  gzclose(gzfile_bin);
  
  return 0;
}

void inspect_file_data_struct(struct file_data *file){
	
	printf("file name:%s\n",file->name);
	printf("file index_name:%s\n",file->index_name);
	printf("file version:%d\n",file->version);
	printf("file subversion:%d\n",file->subversion);
	printf("error:%d\n",file->error);
	/*
	if (file->bin_search==TRUE) printf("file binary search is possible\n");
	else printf("file binary search is not possible\n");
	*/
	
}

// initialize the state for doing writes
// two modes:
// 1 .- new files
// 2 .- add data to files, if they don't exist they are created
int initialize_writes(struct file_data ** file, char *output_name, int mode, int discretize_qual, int flatten_qual, int create_index)
{

// check if the files exists, in case it exists check if it has the
// correct metadata and if it is of the correct version
// in other case exits with an error
	// struct file_data *file = malloc(sizeof(struct write_file));
	if ( *file  == NULL ) {*file=malloc(sizeof(struct file_data));}
	
	(*file)->pos_chunk_gz=0;
	(*file)->discretize_qual=discretize_qual;
	(*file)->flatten_qual=flatten_qual;
	(*file)->create_index=create_index;
	
	int state=check_files(output_name);
	if (state==1) {
		fprintf(stderr,"File is from a different version\n");
		return -1;
	}
	if ((state!=2)&&(state!=0)) {
		fprintf(stderr,"Error %d when checking files\n",state);
		return -1;
	}
	
	// copy the name of the file
	strncpy((*file)->name,output_name,MAXFNAME);

	// open the compressed files
	int error=0;
	int flags=O_WRONLY|O_CREAT|O_TRUNC;
	if (mode==2) flags=O_RDWR;
	// printf("mode:%d\n",mode);

	//set index name
	snprintf((*file)->index_name,MAXFNAME,"%s.index",(*file)->name);
	
    int file_index=-1;
    
	if ((*file)->create_index){
    	//open index file
    	file_index=open((*file)->index_name,flags,0644);

    	if (file_index==-1) return -2;
    }
    
	// open bin file
	int file_bin=open((*file)->name,flags,0644);
	// printf("fd:%d\n",file_bin);
	if (file_bin==-1) {fprintf(stderr,"error opening file_bin for writting:%s\n",strerror(errno));return -2;}
	if (mode==2) {
		long long pos=lseek(file_index,0,SEEK_END);
 		if (pos==-1) {fprintf(stderr,"error going to end of index file %s\n",strerror(errno)); return -2;}
		pos=lseek(file_bin,0,SEEK_END);
 		if (pos==-1) {fprintf(stderr,"error going to end of bin file %s\n",strerror(errno)); return -2;}
		(*file)->pos_chunk_gz=pos;
	}
	
	if ((*file)->create_index){	    
    	// open zlib index file
    	(*file)->gzf_index=gzdopen(file_index,"wb");
    	if ((*file)->gzf_index==NULL) {
    		fprintf(stderr,"error opening gzfile_index for writting:%s\n",gzerror((*file)->gzf_index,&error));
    		return -2;
    	}
    }
    
	// open zlib bin file
	(*file)->gzf_bin=gzdopen(file_bin,"wb");
	if ((*file)->gzf_bin==NULL) {
		fprintf(stderr,"error opening gzfile for writting:%s\n",gzerror((*file)->gzf_bin,&error));
		return -2;
	}

	// initializes the files, writting the metadata
 	if (mode==1) {
	    char header[SEQ_METADATA];
		(*file)->version=VERSION;
		(*file)->subversion=SUBVERSION;
		(*file)->error=0;
		// TODO put correct size
		snprintf(header,SEQ_METADATA-1,"9999UMACOMPRESSEDFORMAT_%d_%d %d %d %d", (*file)->version,(*file)->subversion, 0, 0, 0);
		snprintf(header,SEQ_METADATA-1,"%4ldUMACOMPRESSEDFORMAT_%d_%d %d %d %d",  strlen(header)-4,(*file)->version,(*file)->subversion, 0, 0, 0);
//		snprintf(header,SEQ_METADATA-1,"%4ld%s %ld %ld %ld", strlen(metainfo)-4, (*file)->version,(*file)->subversion, 0, 0, 0);

//		sprintf(header,"  29UMACOMPRESSEDFORMAT_%d_%d 0 0 0\n",(*file)->version,(*file)->subversion);
		int res=gzwrite((*file)->gzf_bin,header,strlen(header));

        if((*file)->create_index)
        {
    		sprintf(header,"UMACOMPRESSEDFORMAT 1 0 0 999999999999 999999999999\n");
    		res=gzwrite((*file)->gzf_index,header,strlen(header));
        }
	}
	(*file)->counter=0;
	
    // prev_time=time(NULL);
	
	// printf("Init writes done\n");
	return 0;
}



int close_writes(struct file_data *file)
{
   gzclose(file->gzf_bin);
   if((file)->create_index)
   {
       gzclose(file->gzf_index);
   }
   
   if ( file  != NULL ) {  free(file);}
}

// open a file for reading and check error
int open_file(char *fname, FILE **file){
    
    *file=fopen(fname,"r");

    if (*file==NULL) { fprintf(stderr,"Error opening fasta file %s, result %d %s\n",fname,errno,strerror(errno));
        return -1;
    };
    
    return 0;
}

// open a file for reading and check error
int close_file(FILE *file){
    
    fclose(file);
    return 0;
}


// removes last \n from string if any
int chomp(char *str){
    
    // printf("LEN: %s, %ld\n",str,strlen(str));
    if (str[strlen(str)-1]=='\n'){
        str[strlen(str)-1]='\0';
        // printf("LEN2: %s, %ld\n",str,strlen(str));
    }
    
}

// split name in name and comments and remove @ or > from first char 
int split_name(char *fname, char *name, char *comments){
    
    char *name_part;
    char *comment_part;
    
    // remove first char (@, >, etc)
    memmove(fname, fname+1, strlen(fname));
    
    // split name by space
    name_part = strtok(fname, " ");

    // get remaining until end of line
    comment_part=strtok(NULL, "\n");

    // assign name and comments
    if(name_part)
    {
        strcpy(name,name_part);
    }else{
        strcpy(name,"");
    }
    
    if(comment_part)
    {
        strcpy(comments,comment_part);
    }else{
        strcpy(comments,"");
    }    
    
    return 1;
}

int check_mem(char **var,int size)
{
    if ( *var == NULL) {*var=(char *)malloc(size);strncpy(*var,"",size);}
}

// read next seq from fastq file
int get_next_seq_fastq(FILE *file, char **name, char **fasta, char **qual, char **comments){
    
    check_mem(name,MAXSEQNAME);
    check_mem(fasta,MAXSEQLENGTH);
    check_mem(qual,MAXSEQLENGTH);
    check_mem(comments,MAXSEQLENGTH);
    
    
    
    
    char fname[MAXSEQNAME];// sequence name
    char qname[MAXSEQNAME];// sequence name
    
    strcpy(*name,"");
    strcpy(*fasta,"");
    strcpy(*qual,"");
    strcpy(*comments,"");

    errno = 0;
    
    // read sequence name line ---------------------------------
    if (fgets( fname, MAXSEQNAME, file )==NULL) {return 0; };

    if (errno) {fprintf(stderr,"Error reading fastq, result %d %s\n",errno,strerror(errno));return INVALID_FASTQ_FORMAT;}
    chomp(fname);

    // check for @ at beginning
    if (fname[0]!='@'){
        fprintf(stderr,"ERROR: Invalid FASTQ format %s. Missing @ in name line.\n", fname);
        return INVALID_FASTQ_FORMAT;
    }
    
    // split name by space in name and comments
    split_name(fname,*name,*comments);
    
    // read fasta line ---------------------------------
    if (fgets( *fasta, MAXSEQLENGTH, file )==NULL) { return 0; };
    
    if (errno) {fprintf(stderr,"Error reading fastq, result %d %s\n",errno,strerror(errno));return INVALID_FASTQ_FORMAT;}
    chomp(*fasta);
    
    
    // read qual name line ---------------------------------
    if (fgets( qname, MAXSEQLENGTH, file )==NULL) { return 0; };
    
    if (errno) {fprintf(stderr,"Error reading fastq, result %d %s\n",errno,strerror(errno));return INVALID_FASTQ_FORMAT;}
    chomp(qname);
    
    // check for + sign at beginning of qual name line
    if (qname[0]!='+'){
        fprintf(stderr,"ERROR: Invalid FASTQ format. Missing + in qual line. %s.\n", *name);
        return INVALID_FASTQ_FORMAT;
    }

    // read qual line  ---------------------------------
    if (fgets( *qual, MAXSEQLENGTH, file )==NULL) { return 0; };
    if (errno) {fprintf(stderr,"Error reading fastq, result %d %s\n",errno,strerror(errno));return INVALID_FASTQ_FORMAT;}

    chomp(*qual);
    
    return 1;
}

// read next seq from fasta file
int get_next_seq_fasta(FILE *file, char *name, char *fasta, char *comments){
    
    char fname[MAXSEQNAME];// sequence name
    char *line;
    if ((line = malloc(MAXSEQLENGTH)) == NULL) {
        puts("Memory allocation error!");
        return EXIT_FAILURE;
    }
    
    // init vars
    strcpy(name,"");
    strcpy(fasta,"");
    strcpy(comments,"");
    
    // read sequence name line ---------------------------------
    if (fgets( fname, MAXSEQNAME, file )==NULL) { return 0; };
    if (errno) {fprintf(stderr,"Error reading fasta, result %d %s\n",errno,strerror(errno));return INVALID_FASTA_FORMAT;}
    chomp(fname);
    
    // check for @ at beginning
    if (fname[0]!='>'){
        fprintf(stderr,"ERROR: Invalid FASTA format %s. Missing @ in name line.\n", fname);
        return INVALID_FASTA_FORMAT;
    }
    
    // split name by space in name and comments
    split_name(fname,name,comments);
    
    // get current pos in file
    fpos_t pos=0;
    
    // read fasta line ---------------------------------
    char c=0;
    int len=0;
    int num_lines=0;
    while (1) {
        // fgetpos( file, &pos );
        
        // inspect first char
        c=fgetc(file);
        ungetc(c,file);

        if (c!='>'){
            // get following line
            if (fgets( line, MAXSEQLENGTH, file )==NULL) { break; };
            if (errno) {fprintf(stderr,"Error reading fasta, result %d %s\n",errno,strerror(errno));return INVALID_FASTA_FORMAT;}
            chomp(line);
            
            // append to fasta
            // strcat(fasta,line);
            
            len = len + sprintf(fasta+len,"%s",line);
            
            if (len>=MAXSEQLENGTH)
            {
                fprintf(stderr,"Error, maximun sequence size error (%d). You can recompile lib with a bigger MAXSEQLENGTH\n",MAXSEQLENGTH);return MAX_SEQ_SIZE_ERROR;
            }
            // fasta[strlen(fasta)]=0;
             // num_lines++;
            // printf("%d\n",num_lines);
            
             // if((num_lines%1000)==0)
             // {
             //     printf("%ld\n,",strlen(fasta));
             // }
        }else{ // name line
            // rewind file and exit
            // fsetpos(file, &pos);
            break;
        }
    }

    // printf("%ld\n,",strlen(fasta));
    
    free(line);
    
    return 1; 
}

// process a fastq file adding it to fbin file
int process_fastq(char *fname, char *efname, char *outname, int discretize_qual, int flatten_qual, int create_index)
{

    // allocate strings
    char *name;
    if ((name = malloc(MAXSEQNAME)) == NULL) {
        puts("Memory allocation error!");
        return EXIT_FAILURE;
    }
    
    char *fasta;
    if ((fasta = malloc(MAXSEQLENGTH)) == NULL) {
        puts("Memory allocation error!");
        return EXIT_FAILURE;
    }
    
    char *qual;
    if ((qual = malloc(MAXSEQLENGTH)) == NULL) {
        puts("Memory allocation error!");
        return EXIT_FAILURE;
    }
    
    char *comments;
    if ((comments = malloc(MAXSEQLENGTH)) == NULL) {
        puts("Memory allocation error!");
        return EXIT_FAILURE;
    }


    char *extras_name;
    if ((extras_name = malloc(MAXSEQLENGTH)) == NULL) {
        puts("Memory allocation error!");
        return EXIT_FAILURE;
    }

    
    char *extras;
    if ((extras = malloc(MAXSEQLENGTH)) == NULL) {
        puts("Memory allocation error!");
        return EXIT_FAILURE;
    }

    char *final_extras;
    if ((final_extras = malloc(MAXSEQLENGTH)) == NULL) {
        puts("Memory allocation error!");
        return EXIT_FAILURE;
    }
    
    
    // char name[MAXSEQLENGTH];
    // char fasta[MAXSEQLENGTH];
    // char qual[MAXSEQLENGTH];
    // char extras[MAXSEQLENGTH];
    static time_t curr_time=0;
    static time_t prev_time=0;
    
    prev_time=time(NULL);
    
    
    FILE *fastq_file=NULL;
    FILE *extras_file=NULL;
    
    int valid=0;
    int res=0;
    int r=0;
    
    // Open fasta and qual files
    if (strcmp(fname,"-")==0){
        fastq_file=stdin;
    }else{
        open_file(fname,&fastq_file);
    }

    if(efname!=NULL)
    {
        
      open_file(efname,&extras_file);
    }

    // open output file
    struct file_data *file=NULL;
    int error2=initialize_writes(&file, outname,1, discretize_qual, flatten_qual,create_index);

    // printf("Init writes\n");
    // read first extra entry
    if(extras_file!=NULL)
    {
      get_next_seq_fasta(extras_file,extras_name,extras,comments);
    }
    
    
    // for each sequence on fastq file
    while (valid=get_next_seq_fastq(fastq_file,&name,&fasta,&qual,&comments)){
        if(valid==1)
        {
            r++;

            // printf("======================\nNAME:%s\nSEQ :%s\nQUAL:%s\n", name,fasta,qual);
            // if(strlen(comments)>0)
            // {
            //     printf("COM :%s\n", comments);
            // }

            strcpy(final_extras,comments);
            

            // check if there are extras available
            if (strcmp(name,extras_name)==0){
                strcat(final_extras,extras);
                
                // read next extras
                if(extras_file!=NULL)
                {
                  get_next_seq_fasta(extras_file,extras_name,extras,comments);
                }
            }
            
            int error_wr=write_seq(file,name,fasta,qual,final_extras);
            if (error_wr!=0)  {res=error_wr; break;};
            // if (error_wr==0)  cnt++;
            
        }else{
            fprintf(stderr,"Invalid sequence found. Aborting import.");
            res=-1;
            break;
        }
        
        if ((r%10000)==0) {
            printf(".");
            fflush(stdout);
            // curr_time=time(NULL);
            //              printf("10k seqs in:%8.0f secs\n",difftime(curr_time,prev_time));
            //             prev_time=curr_time;
        }
        
    }

    curr_time=time(NULL);    
    printf("\nEnd fastq processing. %d seqs in %.0f s. Rate: %8.2f seqs/s\n",r,difftime(curr_time,prev_time),r/difftime(curr_time,prev_time));

    // free mem
    free(name);
    free(fasta);
    free(qual);
    free(comments);
    free(extras_name);
    free(extras);
    free(final_extras);
    
    // close files
    fclose(fastq_file);
    if(extras_file!=NULL)
    {
      fclose(extras_file);
    }
    close_writes(file);
    
    return res;
}


// process a fastq file adding it to fbin file
int process_fasta(char *fname, char *efname, char *outname, int discretize_qual, int flatten_qual, int create_index)
{

    // allocate strings
    char *name;
    if ((name = malloc(MAXSEQNAME)) == NULL) {
        puts("Memory allocation error!");
        return EXIT_FAILURE;
    }
    
    char *fasta;
    if ((fasta = malloc(MAXSEQLENGTH)) == NULL) {
        puts("Memory allocation error!");
        return EXIT_FAILURE;
    }
    
    char *qual;
    if ((qual = malloc(MAXSEQLENGTH)) == NULL) {
        puts("Memory allocation error!");
        return EXIT_FAILURE;
    }
    
    char *comments;
    if ((comments = malloc(MAXSEQLENGTH)) == NULL) {
        puts("Memory allocation error!");
        return EXIT_FAILURE;
    }


    char *extras_name;
    if ((extras_name = malloc(MAXSEQLENGTH)) == NULL) {
        puts("Memory allocation error!");
        return EXIT_FAILURE;
    }

    
    char *extras;
    if ((extras = malloc(MAXSEQLENGTH)) == NULL) {
        puts("Memory allocation error!");
        return EXIT_FAILURE;
    }

    char *final_extras;
    if ((final_extras = malloc(MAXSEQLENGTH)) == NULL) {
        puts("Memory allocation error!");
        return EXIT_FAILURE;
    }
    
    
    static time_t curr_time=0;
    static time_t prev_time=0;
    
    prev_time=time(NULL);
    
    
    FILE *fastq_file=NULL;
    FILE *extras_file=NULL;
    
    int valid=0;
    int res=0;
    int r=0;
    
    // Open fasta and qual files
    if (strcmp(fname,"-")==0){
        fastq_file=stdin;
    }else{
        open_file(fname,&fastq_file);
    }

    if(efname!=NULL)
    {
        
      open_file(efname,&extras_file);
    }

    // open output file
    struct file_data *file=NULL;
    int error2=initialize_writes(&file, outname,1, discretize_qual, flatten_qual,create_index);

    // printf("Init writes\n");
    // read first extra entry
    if(extras_file!=NULL)
    {
      get_next_seq_fasta(extras_file,extras_name,extras,comments);
    }
    
    strcpy(qual,"");
    qual[0]=0;
    
    // for each sequence on fastq file
    while (valid=get_next_seq_fasta(fastq_file,name,fasta,comments)){
        if(valid==1)
        {
            r++;

            // printf("======================\nNAME:%s\nSEQ :%s\nQUAL:%s\n", name,fasta,qual);
            // if(strlen(comments)>0)
            // {
            //     printf("COM :%s\n", comments);
            // }

            strcpy(final_extras,comments);
            

            // check if there are extras available
            if (strcmp(name,extras_name)==0){
                strcat(final_extras,extras);
                
                // read next extras
                if(extras_file!=NULL)
                {
                  get_next_seq_fasta(extras_file,extras_name,extras,comments);
                }
            }
            
            int error_wr=write_seq(file,name,fasta,qual,final_extras);
            if (error_wr!=0)  {res=error_wr; break;};
            // if (error_wr==0)  cnt++;
            
        }else{
            fprintf(stderr,"Invalid sequence found. Aborting import.");
            res=-1;
            break;
        }
        
        if ((r%10000)==0) {
            printf(".");
            fflush(stdout);
            // curr_time=time(NULL);
            //              printf("10k seqs in:%8.0f secs\n",difftime(curr_time,prev_time));
            //             prev_time=curr_time;
        }
        
    }

    curr_time=time(NULL);    
    printf("\nEnd fastq processing. %d seqs in %.0f s. Rate: %8.2f seqs/s\n",r,difftime(curr_time,prev_time),r/difftime(curr_time,prev_time));

    // free mem
    free(name);
    free(fasta);
    free(qual);
    free(comments);
    free(extras_name);
    free(extras);
    free(final_extras);
    
    // close files
    fclose(fastq_file);
    if(extras_file!=NULL)
    {
      fclose(extras_file);
    }
    close_writes(file);
    
    return res;
}

// int process_biofile(char *fname, char *qfname, char *efname, char *outname)
// {
// 
//         char sname[MAXSEQNAME];// sequence name       
//         char qname[MAXSEQNAME];// sequence name       
//         char ename[MAXSEQNAME];// sequence name       
//         char next_sname[MAXSEQNAME];// sequence name  
//         char next_qname[MAXSEQNAME];// sequence name  
//         char next_ename[MAXSEQNAME];// sequence name  
// 
//         char fasta[150000];
//         char qual[150000]; 
//         char extras[150000]; 
//         char extras_used[150000];
//         char next_fcomment[150000];
//         char next_qcomment[150000];
//         char next_ecomment[150000];
//         char tmp[150000];          
//         int extras_bool=TRUE;
// 
//         int cnt=1;
// 
//  sprintf(extras_used,"INITIALIZED");
// 
//         // Open fasta and qual files
//         FILE *file_fasta=fopen(fname,"r");
// 
//         if (file_fasta==NULL) { fprintf(stderr,"error opening fasta file %s, result %d %s\n",fname,errno,strerror(errno));return -2;};                                                                                    
// //      setvbuf(file_fasta,NULL,_IONBF,0);                                                                   
//         FILE *file_qual=fopen(qfname,"r");                                                               
//         if (file_qual==NULL) { fprintf(stderr,"error opening qual file %s, result %d %s\n",qfname,errno,strerror(errno));return -2;};                                                                                 
//         FILE *file_extras=fopen(efname,"r");                                                               
//         if (file_extras==NULL) {fprintf(stderr,"error opening extras file %s, result %d %s\n",efname,errno,strerror(errno)); extras_bool=FALSE;sprintf(extras,"");};
// 
// //      setvbuf(file_qual,NULL,_IONBF,0);                                                                    
//         int error=0;                                                                                         
//         int end=0;  //0 is false                                                                             
//         char *res;                                                                                           
// 
//         // reads the name of the sequence from both
// 
// //      fscanf(file_qual,">%9000s",qname);
// //      fscanf(file_fasta,">%9000s",sname);
// 
// 
//         res=fgets(tmp,150000,file_fasta);
//         if (res!=NULL) {                 
//                 sscanf(tmp,">%9000s",sname);
//                 strncpy(next_fcomment,tmp+strlen(sname)+2,150000);
//         }                                                         
// 
//         res=fgets(tmp,150000,file_qual);
//         if (res!=NULL) {                
//                 sscanf(tmp,">%9000s",qname);
//                 strncpy(next_qcomment,tmp+strlen(qname)+2,150000);
//         }                                                         
// 
//  if ( extras_bool ) {
//      res=fgets(tmp,150000,file_extras);
//      if (res!=NULL) {                
//          sscanf(tmp,">%9000s",ename);
//          strncpy(next_ecomment,tmp+strlen(ename)+2,150000);
//      } else sprintf(ename,"");
//  }
//  printf("extras seq:%s\n",ename);
// 
//         printf("file:%s q:%s seqname:%s qseqname%s efname:%s extras:%s\n",fname, qfname,sname,qname,efname,extras);
//         printf("next_fcomment:%s next_qcomment:%s\n",next_fcomment,next_qcomment);   
// 
//  struct file_data *file=NULL;
//         int error2=initialize_writes(&file, outname,1);
// 
// //      sprintf(next_fcomment,"");
// //      sprintf(next_qcomment,"");
// 
//         while (!end) {
//                 if ( strcmp(sname,qname)!=0 ) {error = -9; goto end;}
//      /*
//      if (extras_bool)
//                  if ( strcmp(sname,ename)!=0 ) {error = -9; goto end;}
//      */
//                 // load the qual and fasta                           
// 
//                 sprintf(fasta,"");
//                 sprintf(fasta,"%s",next_fcomment);
//                 sprintf(next_fcomment,"");        
//                 sprintf(tmp,"");                  
//      res=fasta;                        
//      while (( res!=NULL ) && (tmp[0]!='>' )) {
//          res=fgets(tmp,150000,file_fasta);
//          if ((tmp[0]!='>')&&(res!=NULL)) sprintf (fasta,"%s%s",fasta,tmp);
//          else if (res!=NULL) {sscanf(tmp,">%9000s",next_sname); strncpy(next_fcomment,tmp+strlen(next_sname)+2,sizeof(next_fcomment));}                                                                    
//      }                                                                                            
//                 if (res==NULL) end=1;                                                                        
// 
//                 sprintf(qual,"");
//                 sprintf(qual,"%s",next_qcomment);
//                 sprintf(next_qcomment,"");       
//                 res=qual;                        
//                 sprintf(tmp,"");                 
//                 while (( res!=NULL ) && (tmp[0]!='>' )) {
//                         res=fgets(tmp,150000,file_qual); 
//                         if ((tmp[0]!='>')&&(res!=NULL)) sprintf (qual,"%s%s",qual,tmp);
//                         else if (res!=NULL) {sscanf(tmp,">%9000s",next_qname); strncpy(next_qcomment,tmp+strlen(next_qname)+2,sizeof(next_qcomment));}                                                                    
//                 }
//                 if (res==NULL) end=1;                                                                        
// 
//      // If extra_used!=NULL then it means that it has been used and a new one must be read
//      if (extras_bool && (strcmp(extras_used,"")!=0)) {
//          sprintf(extras,"");
//          sprintf(extras,"%s",next_ecomment);
//          sprintf(next_ecomment,"");
//          res=extras;
//          sprintf(tmp,"");
//          while (( res!=NULL ) && (tmp[0]!='>' )) {
//              res=fgets(tmp,150000,file_extras);
//              if ((tmp[0]!='>')&&(res!=NULL)) sprintf (extras,"%s%s",extras,tmp);
//              else if (res!=NULL) {sscanf(tmp,">%9000s",next_ename); strncpy(next_ecomment,tmp+strlen(next_ename)+2,sizeof(next_ecomment));}                                        
//          }
//          //if (res==NULL) end=1; Extras file can be finished and processing will continue
//      }
// 
//      /* If the name of the name is equal to the name of the actual sequence then it will be used for writting */
//      if ( strcmp(sname,ename)==0 ) {
//         strcpy(extras_used,extras);
//             strcpy(ename,next_ename);                                                                    
//      } else sprintf(extras_used,"");
//      
//                              
//                 int error_wr=write_seq(file,sname, fasta,qual,extras_used);
//                 if (error_wr!=0)  { end=1;error=error_wr; };                                                 
//                 if (error_wr==0)  cnt++;                                                                     
//                 strcpy(sname,next_sname);                                                                    
//                 strcpy(qname,next_qname);                                                                    
// 
//         }
// 
//   // repeat until EOF or error
//   end:                      
//         fclose(file_fasta); 
//         fclose(file_qual);
//         
//         close_writes(file);
//         //fclose(file_index);
// //              print_seqs(seql);
//         return error;            
// }
// 



int init_dicts(char *d_fasta,char *d_qual,int size)
{
	char *dict_f="fasta.dic";
	char *dict_q="qual.dic";
        FILE *f_d_fasta=fopen(dict_f,"r");
	if (f_d_fasta==NULL) { fprintf(stderr,"error opening qual file %s, result %d %s\n",dict_f,errno,strerror(errno));return -2;};
	fread(d_fasta,size,1,f_d_fasta);
	fclose(f_d_fasta);

        FILE *f_d_qual=fopen(dict_q,"r");
	if (f_d_qual==NULL) { fprintf(stderr,"error opening qual file %s, result %d %s\n",dict_q,errno,strerror(errno));return -2;};
	fread(d_qual,size,1,f_d_fasta);
	fclose(f_d_qual);
}



int free_string(char **string){
  if (string!=NULL){
    free(*string);
    *string=NULL;
  }
}

