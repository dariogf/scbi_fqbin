#include <zlib.h>

#define VERSION 1
#define SUBVERSION 0

struct file_data {
	char name[10000];
	char index_name[10000];
	gzFile gzf_bin;
	// int file_bin;
	gzFile gzf_index;
	// int file_index;
	// char file_outname[10000];
	long long pos_chunk_gz;
	// Contains the version and subversion of this file
	int version;
	int subversion;
	// bin_search is true when a binary search can be used.
	int bin_search;
	// Counts the number of sequences written to the bin file, so it can
	// decide where to create a new gz chunk
	long long counter;
	// If there is an error it is stored here so it can be retrieved.
	int error;
};

// two modes:
// 1 .- new files
// 2 .- add data to files, if they don't exist they are created
int initialize_writes(struct file_data ** file, char *output_name, int mode);

/* write_seq writes a sequence to the files f_bin and its index to f_index 
   pos_chunk_gz is the offset of the beggining of the current gz chunk inside the file 
   seq_name is a pointer to the name of the sequence
   fasta, quanta and extras are pointers to strings, must be zero terminated.
   Returns 0 if all goes fine.
 */

void inspect_file_data_struct(struct file_data *file);

// int write_seq(gzFile *f_bin, FILE *f_index, long pos_chunk_gz, char *seq_name, char *fasta, char *quanta, char *extras);
int write_seq(struct file_data *file, char *seq_name, char *fasta, char *qual, char *extras);

int close_writes(struct file_data *file);

/*
   read_seq reads from filename the sequence named seq_name and returns its
   fasta, quanta and extras in those variables.
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
int read_seq(char *filename, char *seq_name, char **fasta, char **quanta, char **extras);

// For doing sequential reads of the whole file:
// int initialize_sequential_reads(struct file_data *filed, char *filename);
int initialize_sequential_reads(struct file_data ** filed, char *filename);

// return -9 on EOF
int read_data_sequential(struct file_data *filed, char **seq_name, char **fasta, char **qual, char **extras);
int close_sequential_reads(struct file_data *filed);


/* process_biofile reads from fname (and fname.quanta) and writes to outname
   (and outname.index) with the binary format 
   Returns 0 if all goes fine */
int process_biofile(char *fname,char *qfname, char *efname, char *outname);

/*
Format definition

Main file that contains chunks compressed in gz
For each sequence the information of that sequence is written with the format:
  28F143CJN01EBIJN 105 312 0

That is:
4 chars for the size of this header, excluding itself, that is, it is the size of
	the rest of the header
sequence name
fasta size
qual size
extras size

The First sequence can be a special sequence with metainfo for this file:
  30UMACOMPRESSEDFORMAT_version 0 0 0 
  27UMACOMPRESSEDFORMAT_1 0 0 0



Index file

Compressed using chunks

At the beggining a special sequence can be used to store metadata
like the number of fields, if a binary search can be used, etc.

That sequence will be:
UMACOMPRESSEDFORMAT version binary_search begin_of_sequential_index  0 0

If binary_search is yes then a metaindex follows to do a fast access to the
index data.
That will be the first sequence of each chunk and its offset inside the file.
(Or perphaps it can be put in another file....)


The rest of the index file will be indexes to the stored sequences, with
the following fields separated by spaces:

F143CJN01ETK00 0 471 

Sequence name
begin of the compressed chunk
offset inside the chunk of the header of that sequence.

*/






