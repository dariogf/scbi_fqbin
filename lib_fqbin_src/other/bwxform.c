/***************************************************************************
*                    Burrows-Wheeler Transform Library
*
*   File    : bwxform.c
*   Purpose : Provides prototypes for functions that apply and reverse the
*             Burrows-Wheeler transform (with or without move to front
*             coding/decoding).  The algorithms implemented are based upon
*             those described in "A Block-sorting Lossless Data Compression
*             Algorithm" by M. Burrows and D.J. Wheeler.
*   Author  : Michael Dipperstein
*   Date    : August 20, 2004
*
****************************************************************************
*   UPDATES
*
*   $Id: bwxform.c,v 1.6 2007/09/17 13:21:19 michael Exp $
*   $Log: bwxform.c,v $
*   Revision 1.6  2007/09/17 13:21:19  michael
*   Changes required for LGPL v3.
*
*   Revision 1.5  2005/11/03 15:01:46  michael
*   Speed up block sorting using the algorithm suggested by the
*   Burrows-Wheeler paper.  Radix sort all rotations by the first
*   two charcters before employing quicksort.
*
*   Revision 1.4  2005/05/02 13:33:41  michael
*   Allocate large arrays on heap instead of stack so that gcc builds code
*   that can handle larger blocks.
*
*   Update e-mail address
*
*   Revision 1.3  2004/08/27 01:24:16  michael
*   Write S[0] index (I) before transformed block to aviod having to
*   find I in a partial block.
*
*   Revision 1.2  2004/08/26 06:16:08  michael
*   Handle partial blocks without need to store block size.  Use size
*   returned by fread() to indicate smaller than standard block.
*
*   Revision 1.1.1.1  2004/08/23 04:34:18  michael
*   Burrows-Wheeler Transform
*
****************************************************************************
*
* bwxform: An ANSI C Burrows-Wheeler Transform/Reverse Transform Routines
* Copyright (C) 2004-2005, 2007 by
* Michael Dipperstein (mdipper@alumni.engr.ucsb.edu)
*
* This file is part of the BWT library.
*
* The BWT library is free software; you can redistribute it and/or modify
* it under the terms of the GNU Lesser General Public License as published
* by the Free Software Foundation; either version 3 of the License, or (at
* your option) any later version.
*
* The BWT library is distributed in the hope that it will be useful, but
* WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser
* General Public License for more details.
*
* You should have received a copy of the GNU Lesser General Public License
* along with this program.  If not, see <http://www.gnu.org/licenses/>.
*
***************************************************************************/

/***************************************************************************
*                             INCLUDED FILES
***************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include "bwxform.h"

/***************************************************************************
*                                CONSTANTS
***************************************************************************/
#define BLOCK_SIZE  4096        /* size of blocks */

#if BLOCK_SIZE > INT_MAX
#error BLOCK_SIZE must be <= INT_MAX and maximum size_t
#endif

/* NOTE: Need to find a way to check for maximum size_t */

/***************************************************************************
*                            TYPE DEFINITIONS
***************************************************************************/
unsigned char block[BLOCK_SIZE];        /* block being (un)transformed */
size_t blockSize;                       /* actual size of block */

/* counters and offsets used for radix sorting with characters */
unsigned int counters[256];
unsigned int offsetTable[256];

/***************************************************************************
*                                 MACROS
***************************************************************************/
/* wraps array index within array bounds (assumes value < 2 * limit) */
#define Wrap(value, limit)      (((value) < (limit)) ? (value) : ((value) - (limit)))

/***************************************************************************
*                               PROTOTYPES
***************************************************************************/
/* move to front functions */

void DoMTF(unsigned char *last, int length);
void UndoMTF(unsigned char *last, int length);

/***************************************************************************
*                                FUNCTIONS
***************************************************************************/

/***************************************************************************
*   Function   : ComparePresorted
*   Description: This comparison function is designed for use with qsort
*                and "block", a global array of "blockSize" unsigned chars.
*                It compares two strings in "block" starting at indices
*                s1 and s2 and ending at indices s1 - 1 and s2 - 1.
*                The strings are assumed to be presorted so that first two
*                characters are known to be matching.
*   Parameters : s1 - The starting index of a string in block
*                s2 - The starting index of a string in block
*   Effects    : NONE
*   Returned   : > 0 if string s1 > string s2
*                0 if string s1 == string s2
*                < 0 if string s1 < string s2
***************************************************************************/
int ComparePresorted(const void *s1, const void *s2)
{
    int offset1, offset2;
    int i;
    int result;

    offset1 = *((int *)s1);
    offset2 = *((int *)s2);

    /***********************************************************************
    * Compare 1 character at a time until there's difference or the end of
    * the block is reached.  Since we're only sorting strings that already
    * match at the first two characters, start with the third character.
    ***********************************************************************/
    for(i = 2; i < blockSize; i++)
    {
        result = (int)block[Wrap((offset1 + i), blockSize)] -
            (int)block[Wrap((offset2 + i), blockSize)];

        if (result != 0)
        {
            return result;
        }
    }

    /* strings are identical */
    return 0;
}

/***************************************************************************
*   Function   : BWXformFile
*   Description: This function performs a Burrows-Wheeler transformation
*                on a file (with optional move to front) and writes the
*                resulting data to the specified output file.  Comments in
*                this function indicate corresponding variables, labels,
*                and sections in "A Block-sorting Lossless Data Compression
*                Algorithm" by M. Burrows and D.J. Wheeler.
*   Parameters : inFile - Name of file to transform
*                outFile - Name of file to write transformed output to
*                mtf - Set to TRUE if move to front coding should be
*                      applied.
*   Effects    : A Burrows-Wheeler transformation (and possibly move to
*                front encoding) is applied to inFile.   The results of
*                the transformation are written to outFile.
*   Returned   : TRUE for success, otherwise FALSE.
***************************************************************************/
int BWXformFile(char *inFile, char *outFile, char mtf)
{
    int i, j, k;
    FILE *fpIn, *fpOut;
    unsigned int *rotationIdx;      /* index of first char in rotation */
    unsigned int *v;                /* index of radix sorted charaters */
    int s0Idx;                      /* index of S0 in rotations (I) */
    unsigned char *last;            /* last characters from sorted rotations */

    /***********************************************************************
    * BLOCK_SIZE arrays are allocated on the heap, because gcc generates
    * code that throws a Segmentation fault when the large arrays are
    * allocated on the stack.
    ***********************************************************************/
    rotationIdx = (unsigned int *)malloc(BLOCK_SIZE * sizeof(unsigned int));

    if (NULL == rotationIdx)
    {
        perror("Allocating array of rotation indices");
        return FALSE;
    }

    v = (unsigned int *)malloc(BLOCK_SIZE * sizeof(unsigned int));

    if (v == rotationIdx)
    {
        perror("Allocating array of sort indices");
        free(rotationIdx);
        return FALSE;
    }

    last = (unsigned char *)malloc(BLOCK_SIZE * sizeof(unsigned char));

    if (NULL == last)
    {
        perror("Allocating array of last characters");
        free(rotationIdx);
        free(v);
        return FALSE;
    }

    /* open input and output files */
    if ((fpIn = fopen(inFile, "rb")) == NULL)
    {
        perror(inFile);
        return FALSE;
    }

    if (outFile == NULL)
    {
        fpOut = stdout;
    }
    else
    {
        if ((fpOut = fopen(outFile, "wb")) == NULL)
        {
            fclose(fpIn);
            perror(outFile);
            return FALSE;
        }
    }

    while((blockSize = fread(block, sizeof(unsigned char), BLOCK_SIZE, fpIn))
        != 0)
    {
        /*******************************************************************
        * Sort the rotated strings in the block.  A radix sort is performed
        * on the first to characters of all the rotated strings (2nd
        * character then 1st).  All rotated strings with matching initial
        * characters are then quicksorted. - Q4..Q7
        *******************************************************************/

        /*** radix sort on second character in rotation ***/

        /* count number of characters for radix sort */
        memset(counters, 0, 256 * sizeof(int));
        for (i = 0; i < blockSize; i++)
        {
            counters[block[i]]++;
        }

        offsetTable[0] = 0;

        for(i = 1; i < 256; i++)
        {
            /* determine number of values before those sorted under i */
            offsetTable[i] = offsetTable[i - 1] + counters[i - 1];
        }

        /* sort on 2nd character */
        for (i = 0; i < blockSize - 1; i++)
        {
            j = block[i + 1];
            v[offsetTable[j]] = i;
            offsetTable[j] = offsetTable[j] + 1;
        }

        /* handle wrap around for string starting at end of block */
        j = block[0];
        v[offsetTable[j]] = i;
        offsetTable[0] = 0;

        /*** radix sort on first character in rotation ***/

        for(i = 1; i < 256; i++)
        {
            /* determine number of values before those sorted under i */
            offsetTable[i] = offsetTable[i - 1] + counters[i - 1];
        }

        for (i = 0; i < blockSize; i++)
        {
            j = v[i];
            j = block[j];
            rotationIdx[offsetTable[j]] = v[i];
            offsetTable[j] = offsetTable[j] + 1;
        }

        /*******************************************************************
        * now rotationIdx contains the sort order of all strings sorted
        * by their first 2 characters.  Use qsort to sort the strings
        * that have their first two characters matching.
        *******************************************************************/
        for (i = 0, k = 0; (i <= UCHAR_MAX) && (k < (blockSize - 1)); i++)
        {
            for (j = 0; (j <= UCHAR_MAX) && (k < (blockSize - 1)); j++)
            {
                int first = k;

                /* count strings starting with ij */
                while ((i == block[rotationIdx[k]]) &&
                    (j == block[Wrap(rotationIdx[k] + 1,  blockSize)]))
                {
                    k++;

                    if (k == blockSize)
                    {
                        /* we've searched the whole block */
                        break;
                    }
                }

                if (k - first > 1)
                {
                    /* there are at least 2 strings staring with ij, sort them */
                    qsort(&rotationIdx[first], k - first, sizeof(int),
                        ComparePresorted);
                }
            }
        }

        /* find last characters of rotations (L) - C2 */
        s0Idx = 0;
        for (i = 0; i < blockSize; i++)
        {
            if (rotationIdx[i] != 0)
            {
                last[i] = block[rotationIdx[i] - 1];
            }
            else
            {
                /* unrotated string 1st character is end of string */
                s0Idx = i;
                last[i] = block[blockSize - 1];
            }
        }

        if (mtf)
        {
            DoMTF(last, blockSize);
        }

        /* write index of end of unrotated string (I) */
        fwrite(&s0Idx, sizeof(int), 1, fpOut);

        /* write out last characters of rotations (L) */
        fwrite(last, sizeof(unsigned char), blockSize, fpOut);
    }

    /* clean up */
    free(rotationIdx);
    free(v);
    free(last);
    fclose(fpIn);
    fclose(fpOut);
    return TRUE;
}



int BWXform(char *inString, char *outString, int mtf)
{
    int i, j, k;
    FILE *fpIn, *fpOut;
    unsigned int *rotationIdx;      /* index of first char in rotation */
    unsigned int *v;                /* index of radix sorted charaters */
    int s0Idx;                      /* index of S0 in rotations (I) */
    unsigned char *last;            /* last characters from sorted rotations */

    /***********************************************************************
    * BLOCK_SIZE arrays are allocated on the heap, because gcc generates
    * code that throws a Segmentation fault when the large arrays are
    * allocated on the stack.
    ***********************************************************************/
    rotationIdx = (unsigned int *)malloc(BLOCK_SIZE * sizeof(unsigned int));

    if (NULL == rotationIdx)
    {
        perror("Allocating array of rotation indices");
        return FALSE;
    }

    v = (unsigned int *)malloc(BLOCK_SIZE * sizeof(unsigned int));

    if (v == rotationIdx)
    {
        perror("Allocating array of sort indices");
        free(rotationIdx);
        return FALSE;
    }

    last = (unsigned char *)malloc(BLOCK_SIZE * sizeof(unsigned char));

    if (NULL == last)
    {
        perror("Allocating array of last characters");
        free(rotationIdx);
        free(v);
        return FALSE;
    }

    strcpy(block,inString);
    blockSize=strlen(inString);
    block[blockSize]='\0';
    // printf("block:%s, SIZE: %ld\n",block,blockSize);
    
    // while((blockSize = fread(block, sizeof(unsigned char), BLOCK_SIZE, fpIn))
    //     != 0)
    // {
        /*******************************************************************
        * Sort the rotated strings in the block.  A radix sort is performed
        * on the first to characters of all the rotated strings (2nd
        * character then 1st).  All rotated strings with matching initial
        * characters are then quicksorted. - Q4..Q7
        *******************************************************************/

        /*** radix sort on second character in rotation ***/

        /* count number of characters for radix sort */
        memset(counters, 0, 256 * sizeof(int));
        for (i = 0; i < blockSize; i++)
        {
            counters[block[i]]++;
        }

        offsetTable[0] = 0;

        for(i = 1; i < 256; i++)
        {
            /* determine number of values before those sorted under i */
            offsetTable[i] = offsetTable[i - 1] + counters[i - 1];
        }

        /* sort on 2nd character */
        for (i = 0; i < blockSize - 1; i++)
        {
            j = block[i + 1];
            v[offsetTable[j]] = i;
            offsetTable[j] = offsetTable[j] + 1;
        }

        /* handle wrap around for string starting at end of block */
        j = block[0];
        v[offsetTable[j]] = i;
        offsetTable[0] = 0;

        /*** radix sort on first character in rotation ***/

        for(i = 1; i < 256; i++)
        {
            /* determine number of values before those sorted under i */
            offsetTable[i] = offsetTable[i - 1] + counters[i - 1];
        }

        for (i = 0; i < blockSize; i++)
        {
            j = v[i];
            j = block[j];
            rotationIdx[offsetTable[j]] = v[i];
            offsetTable[j] = offsetTable[j] + 1;
        }

        /*******************************************************************
        * now rotationIdx contains the sort order of all strings sorted
        * by their first 2 characters.  Use qsort to sort the strings
        * that have their first two characters matching.
        *******************************************************************/
        for (i = 0, k = 0; (i <= UCHAR_MAX) && (k < (blockSize - 1)); i++)
        {
            for (j = 0; (j <= UCHAR_MAX) && (k < (blockSize - 1)); j++)
            {
                int first = k;

                /* count strings starting with ij */
                while ((i == block[rotationIdx[k]]) &&
                    (j == block[Wrap(rotationIdx[k] + 1,  blockSize)]))
                {
                    k++;

                    if (k == blockSize)
                    {
                        /* we've searched the whole block */
                        break;
                    }
                }

                if (k - first > 1)
                {
                    /* there are at least 2 strings staring with ij, sort them */
                    qsort(&rotationIdx[first], k - first, sizeof(int),
                        ComparePresorted);
                }
            }
        }

        /* find last characters of rotations (L) - C2 */
        s0Idx = 0;
        for (i = 0; i < blockSize; i++)
        {
            if (rotationIdx[i] != 0)
            {
                last[i] = block[rotationIdx[i] - 1];
            }
            else
            {
                /* unrotated string 1st character is end of string */
                s0Idx = i;
                last[i] = block[blockSize - 1];
            }
        }
        // printf("ANTES:\n");
        // for (i = 0; i < blockSize; i++)
        // {
        //     printf("%2d,",last[i]);
        // }
        // printf("\nDESPUES:\n");
        
        if (mtf)
        {
            DoMTF(last, blockSize);
        }

        // for (i = 0; i < blockSize; i++)
        // {
        //     printf("%2d,",last[i]);
        // }
        // printf("\n");

        // /* write index of end of unrotated string (I) */
        // fwrite(&s0Idx, sizeof(int), 1, fpOut);
        // 
        // /* write out last characters of rotations (L) */
        // fwrite(last, sizeof(unsigned char), blockSize, fpOut);
    // } //FIN WHILE


        memcpy((char *)outString, (void *)last, sizeof(unsigned char) * blockSize);
        // strncpy(outString,last,blockSize+1);
        // blockSize=strlen(last);
        // printf("block2:%s, SIZE: %ld\n",last,strlen(last));


    /* clean up */
    free(rotationIdx);
    free(v);
    free(last);
    // fclose(fpIn);
    // fclose(fpOut);
    return TRUE;
}

/***************************************************************************
*   Function   : DoMTF
*   Description: This function performs move to front encoding on a block
*                on of data that has already had the Burrows-Wheeler
*                transformation applied to it.  Comments in this function
*                indicate corresponding variables, labels, and sections in
*                "A Block-sorting Lossless Data Compression Algorithm" by
*                M. Burrows and D.J. Wheeler.
*   Parameters : last - pointer an array of "last" characters from
*                       Burrows-Wheeler rotations (L)
*                length - the number of unsigned chars contained in last.
*   Effects    : Move to front encoding is applied on an array of last
*                characters.  The results of the encoding replace the data
*                that was stored in last.
*   Returned   : NONE
***************************************************************************/
void DoMTF(unsigned char *last, int length)
{
    unsigned char list[UCHAR_MAX + 1];      /* list of characters (Y) */
    unsigned char *encoded;                 /* mtf encoded block (R) */
    int i, j;

    /***********************************************************************
    * BLOCK_SIZE arrays are allocated on the heap, because gcc generates
    * code that throws a Segmentation fault when the large arrays are
    * allocated on the stack.
    ***********************************************************************/
    encoded = (unsigned char *)malloc(BLOCK_SIZE * sizeof(unsigned char));

    if (NULL == encoded)
    {
        perror("Allocating array to store MTF encoding");
        return;
    }

    /* start with alphabetically sorted list of characters */
    for(i = 0; i <= UCHAR_MAX; i++)
    {
        list[i] = (unsigned char)i;
    }

    /* move-to-front coding - M1 */
    for (i = 0; i < length; i++)
    {
        /*******************************************************************
        * Find the character in the list of characters.  I do a sequential
        * search because of move to front causes common characters to be
        * near the front of the list.
        *******************************************************************/
        for (j = 0; j <= UCHAR_MAX; j++)
        {
            if (list[j] == last[i])
            {
                /* we found the character */
                encoded[i] = j;
                break;
            }
        }

        /* now move the current character to the front of the list */
        for (; j > 0; j--)
        {
            list[j] = list[j - 1];
        }
        list[0] = last[i];
    }

    /* copy mtf encoded vector of last characters (R) to input */
    memcpy((void *)last, (void *)encoded, sizeof(unsigned char) * length);
    free(encoded);

    return;
}

/***************************************************************************
*   Function   : BWReverseXformFile
*   Description: This function reverses a Burrows-Wheeler transformation
*                on a file (with optional move to front) and writes the
*                resulting data to the specified output file.  Comments in
*                this function indicate corresponding variables, labels,
*                and sections in "A Block-sorting Lossless Data Compression
*                Algorithm" by M. Burrows and D.J. Wheeler.
*   Parameters : inFile - Name of file to reverse transform
*                outFile - Name of file to write reverse transformed
*                          output to
*                mtf - Set to TRUE if move to front decoding should be
*                      applied
*   Effects    : A Burrows-Wheeler reverse transformation (and possibly
*                move to front encoding) is applied to inFile.   The results
*                of the reverse transformation are written to outFile.
*   Returned   : TRUE for success, otherwise FALSE.
***************************************************************************/
int BWReverseXformFile(char *inFile, char *outFile, char mtf)
{
    FILE *fpIn, *fpOut;
    int i, j, sum;
    int count[UCHAR_MAX + 1];   /* count[i] = # of chars in block <= i */
    int *pred;                  /* pred[i] = # of times block[i] appears in
                                   block[0 .. i - 1] */
    unsigned char *unrotated;   /* original block */
    int s0Idx;                  /* index of S0 in rotations (I) */

    /***********************************************************************
    * BLOCK_SIZE arrays are allocated on the heap, because gcc generates
    * code that throws a Segmentation fault when the large arrays are
    * allocated on the stack.
    ***********************************************************************/
    pred = (int *)malloc(BLOCK_SIZE * sizeof(int));

    if (NULL == pred)
    {
        perror("Allocating array of matching predicessors");
        return FALSE;
    }

    unrotated = (unsigned char *)malloc(BLOCK_SIZE * sizeof(unsigned char));

    if (NULL == unrotated)
    {
        perror("Allocating array to store unrotated block");
        free(pred);
        return FALSE;
    }

    /* open input and output files */
    if ((fpIn = fopen(inFile, "rb")) == NULL)
    {
        perror(inFile);
        return FALSE;
    }

    if (outFile == NULL)
    {
        fpOut = stdout;
    }
    else
    {
        if ((fpOut = fopen(outFile, "wb")) == NULL)
        {
            fclose(fpIn);
            perror(outFile);
            return FALSE;
        }
    }

    while(fread(&s0Idx, sizeof(int), 1, fpIn) != 0)
    {
        blockSize = fread(block, sizeof(unsigned char), BLOCK_SIZE, fpIn);

        if(mtf)
        {
            UndoMTF(block, blockSize);
        }

        /* code based on pseudo code from section 4.2 (D1 and D2) follows */
        for(i = 0; i <= UCHAR_MAX; i++)
        {
            count[i] = 0;
        }

        /*******************************************************************
        * Set pred[i] to the number of times block[i] appears in the
        * substring block[0 .. i - 1].  As a useful side effect count[i]
        * will be the number of times character i appears in block.
        *******************************************************************/
        for (i = 0; i < blockSize; i++)
        {
            pred[i] = count[block[i]];
            count[block[i]]++;
        }

        /*******************************************************************
        * Finally, set count[i] to the number of characters in block
        * lexicographically less than i.
        *******************************************************************/
        sum = 0;
        for(i = 0; i <= UCHAR_MAX; i++)
        {
            j = count[i];
            count[i] = sum;
            sum += j;
        }

        /* construct the initial unrotated string (S[0]) */
        i = s0Idx;
        for(j = blockSize - 1; j >= 0; j--)
        {
            unrotated[j] = block[i];
            i = pred[i] + count[block[i]];
        }

        fwrite(unrotated, sizeof(unsigned char), blockSize, fpOut);
    }

    /* clean up */
    free(pred);
    free(unrotated);
    fclose(fpIn);
    fclose(fpOut);
    return TRUE;
}


int BWReverseXform(char *inString, char *outString, int mtf, long size)
{
    // FILE *fpIn, *fpOut;
    int i, j, sum;
    int count[UCHAR_MAX + 1];   /* count[i] = # of chars in block <= i */
    int *pred;                  /* pred[i] = # of times block[i] appears in
                                   block[0 .. i - 1] */
    unsigned char *unrotated;   /* original block */
    int s0Idx;                  /* index of S0 in rotations (I) */

    /***********************************************************************
    * BLOCK_SIZE arrays are allocated on the heap, because gcc generates
    * code that throws a Segmentation fault when the large arrays are
    * allocated on the stack.
    ***********************************************************************/
    pred = (int *)malloc(BLOCK_SIZE * sizeof(int));

    if (NULL == pred)
    {
        perror("Allocating array of matching predicessors");
        return FALSE;
    }

    unrotated = (unsigned char *)malloc(BLOCK_SIZE * sizeof(unsigned char));

    if (NULL == unrotated)
    {
        perror("Allocating array to store unrotated block");
        free(pred);
        return FALSE;
    }

    // while(fread(&s0Idx, sizeof(int), 1, fpIn) != 0)
    //     {
        // blockSize = fread(block, sizeof(unsigned char), BLOCK_SIZE, fpIn);

    blockSize=size;
    strncpy(block,inString,blockSize);
    
        if(mtf)
        {
            UndoMTF(block, blockSize);
        }

        /* code based on pseudo code from section 4.2 (D1 and D2) follows */
        for(i = 0; i <= UCHAR_MAX; i++)
        {
            count[i] = 0;
        }

        /*******************************************************************
        * Set pred[i] to the number of times block[i] appears in the
        * substring block[0 .. i - 1].  As a useful side effect count[i]
        * will be the number of times character i appears in block.
        *******************************************************************/
        for (i = 0; i < blockSize; i++)
        {
            pred[i] = count[block[i]];
            count[block[i]]++;
        }

        /*******************************************************************
        * Finally, set count[i] to the number of characters in block
        * lexicographically less than i.
        *******************************************************************/
        sum = 0;
        for(i = 0; i <= UCHAR_MAX; i++)
        {
            j = count[i];
            count[i] = sum;
            sum += j;
        }

        /* construct the initial unrotated string (S[0]) */
        i = s0Idx;
        for(j = blockSize - 1; j >= 0; j--)
        {
            unrotated[j] = block[i];
            i = pred[i] + count[block[i]];
        }

        // fwrite(unrotated, sizeof(unsigned char), blockSize, fpOut);
    // }

    strncpy(outString,unrotated,blockSize);

    /* clean up */
    free(pred);
    free(unrotated);
    // fclose(fpIn);
    //  fclose(fpOut);
    //  
    return TRUE;
}

/***************************************************************************
*   Function   : UndoMTF
*   Description: This function reverses move to front encoding on a block
*                on of data that has already had the Burrows-Wheeler
*                transformation applied to it.  Comments in this function
*                indicate corresponding variables, labels, and sections in
*                "A Block-sorting Lossless Data Compression Algorithm" by
*                M. Burrows and D.J. Wheeler.
*   Parameters : last - pointer an array of mtf encoded characters from
*                       Burrows-Wheeler rotations.
*                length - the number of unsigned chars contained in last.
*   Effects    : Move to front encoding is reversed on an array of last
*                characters.  The results of the reversal are stored in
*                the array last (L), providing an array of last characters
*                of sorted rotations.
*   Returned   : NONE
***************************************************************************/
void UndoMTF(unsigned char *last, int length)
{
    unsigned char list[UCHAR_MAX + 1];      /* list of characters (Y) */
    unsigned char *encoded;                 /* mtf encoded block (R) */
    int i, j;

    /***********************************************************************
    * BLOCK_SIZE arrays are allocated on the heap, because gcc generates
    * code that throws a Segmentation fault when the large arrays are
    * allocated on the stack.
    ***********************************************************************/
    encoded = (unsigned char *)malloc(BLOCK_SIZE * sizeof(unsigned char));

    if (NULL == encoded)
    {
        perror("Allocating array to store MTF encoding");
        return;
    }

    /* copy last into encoded */
    memcpy((void *)encoded, (void *)last, sizeof(unsigned char) * length);

    /* start with alphabetically sorted list of characters */
    for(i = 0; i <= UCHAR_MAX; i++)
    {
        list[i] = (unsigned char)i;
    }

    /* move-to-front decoding - W2 */
    for (i = 0; i < length; i++)
    {
        /* decode the character */
        last[i] = list[encoded[i]];

        /* now move the current character to the front of the list */
        for (j = encoded[i]; j > 0; j--)
        {
            list[j] = list[j - 1];
        }
        list[0] = last[i];
    }

    free(encoded);
    return;
}
