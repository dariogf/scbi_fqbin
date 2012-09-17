
#include <stdio.h>
#include <string.h>
#include <time.h>


#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <errno.h>

#include <zlib.h>
#include <stdlib.h>
#include "libreria_gz.h"

#define CHUNK 262144

// Maximum file name (including .idx)
#define MAXFNAME 512

// Maximum lenght of the name of a sequence
#define MAXSEQNAME 1024
#define DEBUG 0
#define FALSE 0 
#define TRUE 1 

char dict_fasta[65536];
char dict_qual[65536];

// Maximum size of the metadata of a sequence, including name, lenght of fasta, qual and extras. 
// It should be a maximum of 10000
#define SEQ_METADATA 10000

static time_t curr_time=0;
static time_t prev_time=0;

int write_seq(struct file_data *file, char *seq_name, char *fasta, char *qual, char *extras)
{
	// compress data
	char metainfo[SEQ_METADATA];
  	int error=0;

	
  if (file->gzf_bin==NULL) {fprintf(stderr,"error with gzfile_bin, is NULL :%s\n",gzerror(file->gzf_bin,&error));return -2;}
  
	snprintf(metainfo,SEQ_METADATA-1,"9999%s %ld %ld %ld", seq_name, strlen(fasta), strlen(qual), strlen(extras));
	snprintf(metainfo,SEQ_METADATA-1,"%4ld%s %ld %ld %ld", strlen(metainfo)-4, seq_name, strlen(fasta), strlen(qual), strlen(extras));

	// get begin pos of header
	long beginH=gztell(file->gzf_bin);

	// TODO check gztell
	if (beginH==-1) {fprintf(stderr,"error with pos of beginH of gzfile_bin :%s\n", gzerror(file->gzf_bin,&error)); return -2;}

	// write seq to bin file
	gzwrite(file->gzf_bin, metainfo, strlen(metainfo));

	// TODO check gzwrite
	long beginI=gztell(file->gzf_bin);
  
	if (beginI==-1) {fprintf(stderr,"error with pos of beginI of gzfile :%s\n",gzerror(file->gzf_bin,&error));return -2;}
	
	
	int res=1;
	if (strlen(fasta)>0) res=gzwrite(file->gzf_bin,fasta,strlen(fasta)); //Z_FILTERED);
	
	if ( res==0 ) { fprintf(stderr,"Error when writting fasta\n");return -8;}
	long fastaS=gztell(file->gzf_bin)-beginI;
	
	if (strlen(qual)>0) res=gzwrite(file->gzf_bin,qual,strlen(qual)); //Z_FILTERED);
	
	if ( res==0 ) { fprintf(stderr,"Error when writting qual\n");return -8;}
	long qualS=gztell(file->gzf_bin)-fastaS-beginI;
	
	if (strlen(extras)>0) res=gzwrite(file->gzf_bin,extras,strlen(extras)); //Z_FILTERED);
	
	if ( res==0 ) { fprintf(stderr,"Error when writting extras\n");return -8;}
	long extrasS=gztell(file->gzf_bin)-qualS-fastaS-beginI;
	

	//	add_sequence(&seql,seq_name,pos_chunk_gz,beginI,fastaS,qualS,extrasS);

	// Write index file
	char tmp[SEQ_METADATA];
	sprintf(tmp,"%s %lld %ld\n",seq_name,file->pos_chunk_gz,beginH);
	
	gzwrite(file->gzf_index,tmp,strlen(tmp));

	(file->counter)++;
//	if (counter > 2) fprintf(stderr,"Probando static counter para llamadas desde ruby, valor %d\n",counter);

// create new chunk
  if (((file->counter)%10000)==0) {
	  	curr_time=time(NULL);
		printf("time passed:%ld\n",curr_time-prev_time);
	  	prev_time=curr_time;

		// close current chunk
	  	gzclose(file->gzf_bin);
		
		// open file again
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

	int res=read_seq_header(filed->gzf_bin, header, &fastaS, &qualS, &extrasS);

	if ( res!=0 ) {fprintf(stderr,"SEQ READ incorrect:%d\n",res);return -1;}
	if ( strlen(header)<20 ) {fprintf(stderr,"SEQ READ:Header incorrect:%s. lenght:%ld\n",header,strlen(header));return -1;}

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

	if ( strlen(header)<19 ) {fprintf(stderr,"SEQ READ:Header incorrect:%s.\n",header);return -1;}

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

int read_seq_header(gzFile *gzf_bin, char *seq_name,int *fastaS, int *qualS, int*extrasS)
{
  int header_size=4;
  char hsize[40];
  char tmp[1000];
  char sname[SEQ_METADATA];

  long pos=gzread(gzf_bin,hsize,header_size);
               
  // EOF found
  if ( pos==0 ) return -2;

  // Error reading file
  if ( pos==-1 ) {fprintf(stderr,"error reading header\n");return -1;}

  hsize[pos]=0;
  sscanf(hsize,"%d",&header_size);
  pos=gzread(gzf_bin,tmp,header_size);

  if ( pos==0 ) return -2;
	
  if ( pos==-1 ) {fprintf(stderr,"error reading header\n");return -1;}
	
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
	int res= read_bin_file_metadata(*filed);
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
  if (res==-2) // EOF
		return -9;
		
  if ( *fasta  == NULL ) {*fasta=(char *)malloc(fastaS+1);strncpy(*fasta,"",fastaS);}
  if ( *qual == NULL) {*qual=(char *)malloc(qualS+1);strncpy(*qual,"",qualS);}
  if (( *extras == NULL )&&(extrasS>0)) {*extras=(char *)malloc(extrasS+1);strncpy(*extras,"",extrasS);}

  long pos=gzread(filed->gzf_bin,*fasta,fastaS);
  (*fasta)[fastaS]=0;
  pos=gzread(filed->gzf_bin,*qual,qualS);
  (*qual)[qualS]=0;
  if (extrasS>0) {pos=gzread(filed->gzf_bin,*extras,extrasS);(*extras)[extrasS]=0;}
  return 0;

}
int close_sequential_reads(struct file_data *file_d)
{
  gzclose(file_d->gzf_bin);
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

  int bufsize=150000;

	// allocate memory for return data if necessary
  if ( *fasta  == NULL ) {*fasta=(char *)malloc(bufsize);strncpy(*fasta,"",bufsize);}
  if ( *qual == NULL) {*qual=(char *)malloc(bufsize);strncpy(*qual,"",bufsize);}
  if ( *extras == NULL ) {*extras=(char *)malloc(bufsize);strncpy(*extras,"",bufsize);}

	// calc index name
  snprintf(indexname,MAXFNAME,"%s.index",filename);
  //FILE * filein=fopen(indexname,"r");

	// open index file
  gzFile gzfile_index=gzopen(indexname,"r");
  if (gzfile_index==NULL) {
	fprintf(stderr,"error opening gzfile_index :%s\n",gzerror(gzfile_index,&error));
	return -2;
  }

  // Reads the index to this info, and the offset to its data
  int reads=3;
  while ( reads == 3 ) {
				
				// read a chunk of data from index with the size of tmp
        gzgets(gzfile_index,tmp,sizeof(tmp));
				reads=sscanf(tmp,"%s %lld %lld",sname,&gz_chunk,&beginH);
				
				

				
				if (( reads != 3 ) && ( reads!=EOF )) {
					fprintf(stderr,"Error scanning index: %d\n",reads);
					gzclose(gzfile_index);
					return -3;
				}
				
				// sequence was finally found, exit loop
        if ( strncmp(sname, seq_name,MAXSEQNAME)==0) reads=999; // to get out, seq found
  }
	
	// close index file
  gzclose(gzfile_index);
	
	// maybe sequence was not found
	// fprintf(stderr,"Sequence not found\n");
  if (reads==EOF) {return -4;}

	// We get here if sequence was found
	
	// open bin file to extract data
  int dataf=open(filename, O_RDONLY);

	// seek to chunk pos
	// TODO- ¿como se salta el chunk?
  res=lseek(dataf,gz_chunk,SEEK_SET);

	// TODO check res
  gzFile gzfile_bin=gzdopen(dataf,"r");
  
	// seek to seq inside chunk
  res=gzseek(gzfile_bin,beginH,SEEK_SET);
	// TODO check res

//  fasta=malloc(fastaO+1);
//  qual=malloc(qualO+1);
//  extras=malloc(extrasO+1);
//  long pos=gzread(gzfile_bin,header,4); 
// read sequence header 

  res=read_seq_header(gzfile_bin,NULL, &fastaS, &qualS, &extrasS);

	

  long pos=gzread(gzfile_bin,*fasta,fastaS);

  (*fasta)[fastaS]=0;

  pos=gzread(gzfile_bin,*qual,qualS);
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
int initialize_writes(struct file_data ** file, char *output_name, int mode)
{

// check if the files exists, in case it exists check if it has the
// correct metadata and if it is of the correct version
// in other case exits with an error
	// struct file_data *file = malloc(sizeof(struct write_file));
	if ( *file  == NULL ) {*file=malloc(sizeof(struct file_data));}
	
	(*file)->pos_chunk_gz=0;
	
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
	
	//open index file
	int file_index=open((*file)->index_name,flags,0644);

	if (file_index==-1) return -2;
	
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
	
	// open zlib index file
	(*file)->gzf_index=gzdopen(file_index,"w");
	if ((*file)->gzf_index==NULL) {
		fprintf(stderr,"error opening gzfile_index for writting:%s\n",gzerror((*file)->gzf_index,&error));
		return -2;
	}
	
	// open zlib bin file
	(*file)->gzf_bin=gzdopen(file_bin,"w");
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

		sprintf(header,"UMACOMPRESSEDFORMAT 1 0 0 999999999999 999999999999\n");
		res=gzwrite((*file)->gzf_index,header,strlen(header));
	}
	(*file)->counter=0;
	
	// printf("Init writes done\n");
	return 0;
}



int close_writes(struct file_data *file)
{
   gzclose(file->gzf_bin);
   gzclose(file->gzf_index);
}


int process_biofile(char *fname, char *qfname, char *efname, char *outname)
{

        char sname[MAXSEQNAME];// sequence name       
        char qname[MAXSEQNAME];// sequence name       
        char ename[MAXSEQNAME];// sequence name       
        char next_sname[MAXSEQNAME];// sequence name  
        char next_qname[MAXSEQNAME];// sequence name  
        char next_ename[MAXSEQNAME];// sequence name  

        char fasta[150000];
        char qual[150000]; 
        char extras[150000]; 
	char extras_used[150000];
        char next_fcomment[150000];
        char next_qcomment[150000];
        char next_ecomment[150000];
        char tmp[150000];          
        int extras_bool=TRUE;

        int cnt=1;

	sprintf(extras_used,"INITIALIZED");

        // Open fasta and qual files
        FILE *file_fasta=fopen(fname,"r");

        if (file_fasta==NULL) { fprintf(stderr,"error opening fasta file %s, result %d %s\n",fname,errno,strerror(errno));return -2;};                                                                                    
//      setvbuf(file_fasta,NULL,_IONBF,0);                                                                   
        FILE *file_qual=fopen(qfname,"r");                                                               
        if (file_qual==NULL) { fprintf(stderr,"error opening qual file %s, result %d %s\n",qfname,errno,strerror(errno));return -2;};                                                                                 
        FILE *file_extras=fopen(efname,"r");                                                               
        if (file_extras==NULL) {fprintf(stderr,"error opening extras file %s, result %d %s\n",efname,errno,strerror(errno)); extras_bool=FALSE;sprintf(extras,"");};

//      setvbuf(file_qual,NULL,_IONBF,0);                                                                    
        int error=0;                                                                                         
        int end=0;  //0 is false                                                                             
        char *res;                                                                                           

        // reads the name of the sequence from both

//      fscanf(file_qual,">%9000s",qname);
//      fscanf(file_fasta,">%9000s",sname);


        res=fgets(tmp,150000,file_fasta);
        if (res!=NULL) {                 
                sscanf(tmp,">%9000s",sname);
                strncpy(next_fcomment,tmp+strlen(sname)+2,150000);
        }                                                         

        res=fgets(tmp,150000,file_qual);
        if (res!=NULL) {                
                sscanf(tmp,">%9000s",qname);
                strncpy(next_qcomment,tmp+strlen(qname)+2,150000);
        }                                                         

	if ( extras_bool ) {
		res=fgets(tmp,150000,file_extras);
		if (res!=NULL) {                
			sscanf(tmp,">%9000s",ename);
			strncpy(next_ecomment,tmp+strlen(ename)+2,150000);
		} else sprintf(ename,"");
	}
	printf("extras seq:%s\n",ename);

        printf("file:%s q:%s seqname:%s qseqname%s efname:%s extras:%s\n",fname, qfname,sname,qname,efname,extras);
        printf("next_fcomment:%s next_qcomment:%s\n",next_fcomment,next_qcomment);   

	struct file_data *file=NULL;
        int error2=initialize_writes(&file, outname,1);

//      sprintf(next_fcomment,"");
//      sprintf(next_qcomment,"");

        while (!end) {
                if ( strcmp(sname,qname)!=0 ) {error = -9; goto end;}
		/*
		if (extras_bool)
	                if ( strcmp(sname,ename)!=0 ) {error = -9; goto end;}
		*/
                // load the qual and fasta                           

                sprintf(fasta,"");
                sprintf(fasta,"%s",next_fcomment);
                sprintf(next_fcomment,"");        
                sprintf(tmp,"");                  
		res=fasta;                        
		while (( res!=NULL ) && (tmp[0]!='>' )) {
			res=fgets(tmp,150000,file_fasta);
			if ((tmp[0]!='>')&&(res!=NULL)) sprintf (fasta,"%s%s",fasta,tmp);
			else if (res!=NULL) {sscanf(tmp,">%9000s",next_sname); strncpy(next_fcomment,tmp+strlen(next_sname)+2,sizeof(next_fcomment));}                                                                    
		}                                                                                            
                if (res==NULL) end=1;                                                                        

                sprintf(qual,"");
                sprintf(qual,"%s",next_qcomment);
                sprintf(next_qcomment,"");       
                res=qual;                        
                sprintf(tmp,"");                 
                while (( res!=NULL ) && (tmp[0]!='>' )) {
                        res=fgets(tmp,150000,file_qual); 
                        if ((tmp[0]!='>')&&(res!=NULL)) sprintf (qual,"%s%s",qual,tmp);
                        else if (res!=NULL) {sscanf(tmp,">%9000s",next_qname); strncpy(next_qcomment,tmp+strlen(next_qname)+2,sizeof(next_qcomment));}                                                                    
                }
                if (res==NULL) end=1;                                                                        

		// If extra_used!=NULL then it means that it has been used and a new one must be read
		if (extras_bool && (strcmp(extras_used,"")!=0)) {
			sprintf(extras,"");
			sprintf(extras,"%s",next_ecomment);
			sprintf(next_ecomment,"");
			res=extras;
			sprintf(tmp,"");
			while (( res!=NULL ) && (tmp[0]!='>' )) {
				res=fgets(tmp,150000,file_extras);
				if ((tmp[0]!='>')&&(res!=NULL)) sprintf (extras,"%s%s",extras,tmp);
				else if (res!=NULL) {sscanf(tmp,">%9000s",next_ename); strncpy(next_ecomment,tmp+strlen(next_ename)+2,sizeof(next_ecomment));}                                        
			}
			//if (res==NULL) end=1; Extras file can be finished and processing will continue
		}

		/* If the name of the name is equal to the name of the actual sequence then it will be used for writting */
		if ( strcmp(sname,ename)==0 ) {
		   strcpy(extras_used,extras);
	           strcpy(ename,next_ename);                                                                    
		} else sprintf(extras_used,"");

                int error_wr=write_seq(file,sname, fasta,qual,extras_used);
                if (error_wr!=0)  { end=1;error=error_wr; };                                                 
                if (error_wr==0)  cnt++;                                                                     
                strcpy(sname,next_sname);                                                                    
                strcpy(qname,next_qname);                                                                    

        }

  // repeat until EOF or error
  end:                      
        fclose(file_fasta); 
        fclose(file_qual);  
        close_writes(file);
        //fclose(file_index);
//              print_seqs(seql);
        return error;            
}


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




