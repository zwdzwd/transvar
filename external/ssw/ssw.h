/*  ssw.h
 *  copied from SSW library
 *  PloS One, 2013 DOI: 10.1371/journal.pone.0082138
 *  
 *  Created by Mengyao Zhao on 6/22/10.
 *  Copyright 2010 Boston College. All rights reserved.
 *	Version 0.1.4
 *	Last revision by Mengyao Zhao on 01/30/13.
 */

#ifndef SSW_H
#define SSW_H

#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <emmintrin.h>


typedef struct {
  uint32_t* seq;
  int32_t length;
} cigar;

typedef struct _profile {
  __m128i* profile_byte;	/* 0: none */
  __m128i* profile_word;	/* 0: none */
  const int8_t* read;
  const int8_t* mat;
  int32_t readLen;
  int32_t n;
  uint8_t bias;
} s_profile;


/* Create the query profile using the query sequence.
   read	pointer to the query sequence; the query sequence needs to be numbers
   readLen	length of the query sequence
   mat	pointer to the substitution matrix; mat needs to be corresponding to the read sequence
   n the number of rows/columns in mat (mat has n*n elements)
   score_size	estimated Smith-Waterman score; if your estimated best alignment score is surely < 255 please set 0; if 
   your estimated best alignment score >= 255, please set 1; if you don't know, please set 2 

   example for parameter read and mat:
   If the query sequence is: ACGTATC, the sequence that read points to can be: 1234142
   Then if the penalty for match is 2 and for mismatch is -2, the substitution matrix of parameter mat will be:
    A  C  G  T  
    2 -2 -2 -2 A
   -2  2 -2 -2 C
   -2 -2  2 -2 G
   -2 -2 -2  2 T
   mat is the pointer to the array {2, -2, -2, -2, -2, 2, -2, -2, -2, -2, 2, -2, -2, -2, -2, 2} */
s_profile* ssw_init (const int8_t* read,
		     const int32_t readLen,
		     const int8_t* mat,
		     const int32_t n,
		     const int8_t score_size);

static inline void ssw_destroy (s_profile* p) {
  free(p->profile_byte);
  free(p->profile_word);
  free(p);
}

/* structure of the alignment result
   score1	the best alignment score
   score2	sub-optimal alignment score
   ref_begin1	0-based best alignment beginning position on reference;	ref_begin1 = -1 when the best alignment beginning 
   position is not available
   ref_end1	0-based best alignment ending position on reference
   read_begin1	0-based best alignment beginning position on read; read_begin1 = -1 when the best alignment beginning 
   position is not available
   read_end1	0-based best alignment ending position on read
   read_end2	0-based sub-optimal alignment ending position on read
   cigar	best alignment cigar; stored the same as that in BAM format, high 28 bits: length, low 4 bits: M/I/D (0/1/2); 
   cigar = 0 when the best alignment path is not available
   cigarLen	length of the cigar string; cigarLen = 0 when the best alignment path is not available
*/
typedef struct {
  uint16_t score1;
  uint16_t score2;
  int32_t ref_begin1;
  int32_t ref_end1;
  int32_t read_begin1;
  int32_t read_end1;
  int32_t ref_end2;
  uint32_t* cigar;
  int32_t cigarLen;
} s_align;

static inline void free_align(s_align* a) {
  free(a->cigar);
  free(a);
}

/* Striped Smith-Waterman alignment.
   ref	target sequence; the target sequence needs to be numbers and corresponding to the mat parameter of function ssw_init
   refLen length of the target sequence
   weight_gapO	the absolute value of gap open penalty  
   weight_gapE	the absolute value of gap extension penalty

   flag	bitwise FLAG; (from high to low)
   bit 5: when setted as 1, function ssw_align will return the best alignment beginning position;
   bit 6: when setted as 1, if (ref_end1 - ref_begin1 < filterd && read_end1 - read_begin1 < filterd), (whatever bit 5 is setted) the function will return the best alignment beginning position and cigar;
   bit 7: when setted as 1, if the best alignment score >= filters, (whatever bit 5 is setted) the function will return the best alignment beginning position and cigar;
   bit 8: when setted as 1, (whatever bit 5, 6 or 7 is setted) the function will always return the best alignment beginning position and cigar. When flag == 0, only the optimal and sub-optimal scores and the optimal alignment ending position will be returned.

   filters: score filter, when bit 7 of flag is setted as 1 and bit 8 is setted as 0, filters will be used.
   filterd: distance filter, when bit 6 of flag is setted as 1 and bit 8 is setted as 0, filterd will be used
   maskLen: The distance between the optimal and suboptimal alignment ending position >= maskLen. We suggest using readLen/2, if you don't have special concerns. Note: maskLen has to be >= 15, otherwise this function will NOT return the suboptimal alignment information. After locating the optimal alignment ending position, the suboptimal alignment score can be heuristically found by checking the second largest score in the array that contains the maximal score of each column of the SW matrix. In order to avoid picking the scores that belong to the alignments sharing the partial best alignment, SSW C library masks the reference loci nearby (mask length = maskLen) the best alignment ending position and locates the second largest score from the unmasked elements.
   
   Whatever the parameter flag is set, this function will at least return the optimal and sub-optimal alignment score, and the optimal alignment ending positions on target and query sequences.
   If both bit 6 and 7 of the flag are set while bit 8 is not, the function will return cigar only when both criteria are fulfilled.

   All returned positions are 0-based coordinate.
*/
s_align* ssw_align (const s_profile* prof, 
		    const int8_t* ref, 
		    int32_t refLen, 
		    const uint8_t weight_gapO, 
		    const uint8_t weight_gapE, 
		    const uint8_t flag,	
		    const uint16_t filters,
		    const int32_t filterd,
		    const int32_t maskLen);


/* user defined functions */

#define ALIT 20U

#define S_MATCH 1
#define S_MISMATCH 2
#define S_GAPOPEN 2
#define S_GAPEXT 1
#define S_GAPOPEN_LONGGAP 10
#define S_GAPEXT_LONGGAP 0
#define S_NOISE 20

/* tolerate two mismatches from seed extension */
#define MIST 2

s_align* _ssw_align(int8_t *query_seq, int query_len, int8_t *ref_seq, int ref_len);

#endif
