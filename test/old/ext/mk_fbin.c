#include "libfbin.h"
#include <stdio.h>
#include <ctype.h>


/*******************************************************/
/* main                                                */
/*******************************************************/
int main(int argc, char *argv[])
{
  // check params
  if (argc!=5)
  {
    printf("Usage %s fasta_file qual_file extras_file output_file\n\n",argv[0]);
    return -1;
  }
 
  // process file
//  int res=process_biofile(argv[1],argv[2],argv[3],argv[4]);
  int res=process_biofile(argv[1],argv[2],argv[3] ,argv[4]);

  return res;
}

