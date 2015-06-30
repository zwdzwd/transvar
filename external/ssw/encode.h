/* The MIT License
 * Copyright (c) 2015 by The University of Texas MD Anderson Cancer Center (kchen3@mdanderson.org)

 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 * 
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
 * BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
 * ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 * CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 * 
 * This file contains function for encoding/decoding/transfering
 * between different coding methods for DNA sequence.
 *
 * 1. nt256char, native
 * Each base occupies 8 bits.
 * A, C, G, T are just ASCII codes as respresented in char
 *
 * 2. nt256int8
 * Each base occupies 8 bits.
 * A, C, G, T are 0, 1, 2, 3. This is used in the SSW library
 *
 * 3. nt16
 * Each base occupies 4 bits
 * This is used in samtools/bwa.
 * A, C, G, T is 1 (0001), 2 (0010), 4 (0100) and 8 (1000)
 * 
 * 4. nt8
 * Each base occupies 3 bits
 * The downside is that the encoding and decoding would be
 * slower due to non-power-of-2 of 3. Some bases may be
 * astride two bytes. Would be faster if implemented as
 * nt4 + mask array
 * A is 0 (000)
 * C is 3 (001)
 * G is 6 (110)
 * T is 7 (111)
 * N is 5 (101) or 2 (010)
 * 4 (100) and 3 (011) is reserved for future, currently put N
 * 
 * 5. nt4
 * Each base occupies 2 bits.
 * A, C, G, T are 0, 2, 1, 3
 * This has the highest compression rate.
 * But this doesn't allow ambiguous characters such as 'N'.
 *
 * All codings are little-endian,
 * i.e., first base resides on the highest (left-most) bits
 * and last base resides on the lowest (right-most) bits
 * All bitarrays are left-aligned.
 *
 * e.g, ATGC
 * bit-pattern in nt16: [AT] [GC]
 * bit-pattern in nt4: [ATGC]
 * 
 */

#ifndef _ENCODE_H
#define _ENCODE_H

#include <limits.h>
#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#define LEN_NT256_TO_NT4(l) ((l)+3)>>2
#define LEN_NT256_TO_NT8(l) ((l)*3+7)>>3
#define LEN_NT256_TO_NT16(l) ((l)+1)>>1
#define LEN_NT16_TO_NT4(l)  ((l)+1)>>1

/**************************
 * operations on bitarray *
 **************************/

/* print bit array */
static inline void bitarr_show(const unsigned char *bitarr, size_t bitarr_len) {
  const unsigned char *byte;
  for ( byte = bitarr; bitarr_len--; ++byte ) {
    unsigned char mask;
    for ( mask = (1 << (CHAR_BIT - 1)); mask; mask >>= 1 ) {
      putchar(mask & *byte ? '1' : '0');
    }
    putchar(' ');
  }
  putchar('\n');
}

/* shift the bit array to the left by step */
static inline void
bitarr_shiftl(uint8_t *bitarr, size_t l, size_t step) {

  if (step >> 3) {
    size_t i;
    size_t b = step >> 3;
    for(i=0; i<l-b; ++i) bitarr[i] = bitarr[i+b];
    for(i=l-b; i<l; ++i) bitarr[i] = 0;
  }
  if (step &= 7) {
    uint8_t *byte;
    for (byte=bitarr; l--; byte++) {
      uint8_t cross = 0;
      if (l) cross = (byte[1] >> (8 - step)); /* skip the last byte which does not have neighbor on his right */
      *byte <<= step;
      *byte |= cross;
    }
  }
}

/* shift the bit array to the right by step */
static inline void
bitarr_shiftr(uint8_t *bitarr, size_t l, size_t step) {

  if (step >> 3) {
    size_t i;
    size_t b = step >> 3;
    for(i=l-b-1; i--;) bitarr[i+b] = bitarr[i];
    for(i=0; i<b; ++i) bitarr[i] = 0;
  }
  if (step &= 7) {
    uint8_t *byte;
    for (byte=bitarr+l-1; l--; byte--) {
      uint8_t cross = 0;
      if (byte != bitarr) cross = (byte[-1] << (8 - step)); /* skip the first byte which does not have neighbor on his left */
      *byte >>= step;
      *byte |= cross;
    }    
  }
}

static inline uint8_t*
bitarr_clone(const uint8_t *b, size_t l){
  uint8_t *c = (uint8_t*) malloc(l);
  memcpy(c, b, l);
  return c;
}

static inline void
bitarr_not(uint8_t *b, size_t l) {
  for (;l--;b++) *b = ~*b;
}

static inline uint8_t* bitarr_or(const uint8_t *bitarr1, const uint8_t *bitarr2, size_t bitarr_len) {
  uint8_t *out = (uint8_t*) malloc(bitarr_len);
  size_t i;
  for (i=0; i<bitarr_len; ++i) out[i] = bitarr1[i] | bitarr2[i];
  return out;
}

static inline void
bitarr_or_eq(uint8_t *dest, const uint8_t *src, size_t l) {
  size_t i;
  for (i=0; i<l; ++i) dest[i] |= src[i];
  return;
}

static inline void _bitarr_and(uint8_t *out, const uint8_t *b, const uint8_t *m, size_t l) {
  size_t i;
  for (i=0; i<l; ++i) out[i] = b[i] & m[i];
}

static inline uint8_t* bitarr_and(const uint8_t *b, const uint8_t *m, size_t l) {
  uint8_t *out = (uint8_t*) malloc(l);
  _bitarr_and(out, b, m, l);
  return out;
}

/* print bit array as char array */
static inline void bitarr_printchar(unsigned char *bitarr, size_t bitarr_len) {
  unsigned char *byte;
  for (byte=bitarr; bitarr_len--; byte++)
    printf("%c", *byte);
  printf("\n");
}

static inline void bitarr_reverse_ip(void *bitarr, size_t len) {
  unsigned char *p = bitarr;
  unsigned char *q = bitarr+len;
  for(--q; p<q; ++p, --q) {
    *p = *p ^ *q;
    *q = *p ^ *q;
    *p = *p ^ *q;
  }
}


/******************
 * nt4 operations *
 ******************
 *
 * 'A' <=> 0x00 00000000
 * 'C' <=> 0x40 01000000
 * 'G' <=> 0x80 10000000
 * 'T' <=> 0xC0 11000000
 * 
 * */

extern const uint8_t nt256char_to_nt4_table[128];
extern const int8_t  nt4_to_nt256int8_table[256];
extern const uint8_t nt4_to_nt16_table[256];
extern const uint8_t nt256char_to_nt16_table[128];
extern const uint8_t nt16_rev_table[256];
extern const uint8_t nt8_rev_table[256];
extern const uint8_t nt4_rev_table[256];
extern const int8_t  nt256char_to_nt256int8_table[128];
extern const int8_t  nt16_to_nt256int8_table[256];
extern const uint8_t nt256char_to_nt8_table[256];
extern const char    nt4_to_nt256char_table[256];
extern const char    nt8_to_nt256char_table[256];
extern const char    nt16_to_nt256char_table[256];
extern const char    nt256int8_to_nt256char_table[4];
extern const char    nt256char_rev_table[128];

static inline void nt256char_rev_ip(char *seq_nt256char, size_t len) {
  bitarr_reverse_ip(seq_nt256char, len);
  size_t i;
  for (i=0; i<len; ++i) {
    seq_nt256char[i] = nt256char_rev_table[(unsigned) seq_nt256char[i]];
  }
}

static inline void
_nt256char_encode_nt4(uint8_t *seq_nt4, const char *seq_nt256char, size_t len) {
  size_t l = LEN_NT256_TO_NT4(len);
  const char *p;
  seq_nt4[0] = 0;
  for(p=seq_nt256char+len-1; len--; p--) {
    bitarr_shiftr(seq_nt4, l, 2);
    seq_nt4[0] |= nt256char_to_nt4_table[(unsigned char)*p];
  }
}

static inline uint8_t*
nt256char_encode_nt4(const char *seq_nt256char, size_t len) {
  uint8_t *seq_nt4 = (uint8_t*) malloc(LEN_NT256_TO_NT4(len));
  _nt256char_encode_nt4(seq_nt4, seq_nt256char, len);
  return seq_nt4;
}

#define cstr_encode_nt4(cstr) nt256char_encode_nt4((cstr), strlen(cstr))


static inline void
_nt4_decode_nt256char(char *seq_nt256char, const uint8_t *seq_nt4, size_t len) {
  size_t i;
  for (i=0; i<len; i++) {
    seq_nt256char[i] = nt4_to_nt256char_table[(uint8_t)(seq_nt4[i>>2]<<((i&3)<<1))];
  }
}

static inline char*
nt4_decode_nt256char(const uint8_t *seq_nt4, size_t len) {
  char *seq_nt256char = (char*) malloc(len);
  _nt4_decode_nt256char(seq_nt256char, seq_nt4, len);
  return seq_nt256char;
}

static inline char*
nt4_decode_cstr(const uint8_t *seq_nt4, size_t len) {
  char *seq_nt256char = (char*) malloc(len+1);
  _nt4_decode_nt256char(seq_nt256char, seq_nt4, len);
  seq_nt256char[len] = 0;
  return seq_nt256char;
}


static inline void
_nt4_decode_nt256int8(int8_t *seq_nt256int8, const uint8_t *seq_nt4, size_t len) {
  size_t i;
  for (i=0; i<len; i++) {
    seq_nt256int8[i] = nt4_to_nt256int8_table[(uint8_t)(seq_nt4[i>>2]<<((i&3)<<1))];
  }
}

static inline int8_t*
nt4_decode_nt256int8(const uint8_t *seq_nt4, size_t len) {
  int8_t *seq_nt256int8 = (int8_t*) malloc(len);
  _nt4_decode_nt256int8(seq_nt256int8, seq_nt4, len);
  return seq_nt256int8;
}


static inline void
_nt4_decode_nt16(uint8_t *seq_nt16, const uint8_t *seq_nt4, size_t len) {

  memset(seq_nt16, 0, (len+1)>>1);
  size_t i;
  for (i=0; i<len; i++) {
    seq_nt16[i>>1] |= nt4_to_nt16_table[(uint8_t) (seq_nt4[i>>2]<<((i&3)<<1))] >> ((i&1)<<2);
  }
  return;
}

static inline uint8_t*
nt4_decode_nt16(uint8_t *seq_nt4, size_t len) {
  uint8_t *seq_nt16 = (uint8_t*) malloc(LEN_NT256_TO_NT16(len));
  _nt4_decode_nt16(seq_nt16, seq_nt4, len);
  return seq_nt16;
}

/* in-place reverse complementary */
static inline void
nt4_rev_ip(uint8_t *seq, size_t len) {
  size_t l=LEN_NT256_TO_NT4(len);
  bitarr_reverse_ip(seq, l);
  size_t i;
  for(i=0; i<l; i++) seq[i] = nt4_rev_table[seq[i]];

  bitarr_shiftl(seq, l, (l<<3)-(len<<1));
}

static inline void
_nt4_rev(uint8_t *seq_nt4_rev, const uint8_t *seq_nt4, size_t len) {
  size_t l=LEN_NT256_TO_NT4(len);
  size_t i;
  for (i=0; i<l; i++) seq_nt4_rev[i] = nt4_rev_table[seq_nt4[l-i-1]];
  bitarr_shiftl(seq_nt4_rev, l, (l<<3)-(len<<1));
}

static inline uint8_t*
nt4_rev(uint8_t *seq_nt4, size_t len) {
  uint8_t *seq_nt4_rev = (uint8_t*) malloc(LEN_NT256_TO_NT4(len));
  _nt4_rev(seq_nt4_rev, seq_nt4, len);
  return seq_nt4_rev;
}



/*******************
 * nt16 operations *
 *******************
 *
 * 'A' <=> 0x01 00000001
 * 'C' <=> 0x02 00000010
 * 'G' <=> 0x04 00000100
 * 'T' <=> 0x08 00001000
 * 
 * */

static inline void
_nt256char_encode_nt16(uint8_t *seq_nt16, const char *seq_nt256char, size_t len) {
  size_t i;
  memset(seq_nt16, 0, LEN_NT256_TO_NT16(len));
  for (i=0; i<len; i++) {
    seq_nt16[i>>1] |= nt256char_to_nt16_table[(unsigned) seq_nt256char[i]] >> ((i&1)<<2);
  }
}

static inline uint8_t*
nt256char_encode_nt16(const char *seq_nt256char, size_t len) {
  uint8_t *seq_nt16 = (uint8_t*) malloc(LEN_NT256_TO_NT16(len));
  _nt256char_encode_nt16(seq_nt16, seq_nt256char, len);
  return seq_nt16;
}

#define cstr_encode_nt16(cstr) nt256char_encode_nt16((cstr), strlen((cstr)))

static inline void
_nt16_decode_nt256char(void *seq_nt256char, const uint8_t *seq_nt16, size_t len) {
  size_t i;
  for (i=0; i<len; i++) {
    ((unsigned char*) seq_nt256char)[i] = nt16_to_nt256char_table[(uint8_t) (seq_nt16[i>>1] << ((i&1)<<2))];
  }
}

static inline char*
nt16_decode_nt256char(const uint8_t *seq_nt16, size_t len) {
  char *seq_nt256char = (char*) malloc(len);
  _nt16_decode_nt256char(seq_nt256char, seq_nt16, len);
  return seq_nt256char;
}

static inline char*
nt16_decode_cstr(const uint8_t *seq_nt16, size_t len) {
  char *seq_nt256char = (char*) malloc(len+1);
  _nt16_decode_nt256char(seq_nt256char, seq_nt16, len);
  seq_nt256char[len] = 0;
  return seq_nt256char;
}

static inline void
_nt16_decode_nt256int8(int8_t *seq_nt256int8, const uint8_t *seq_nt16, size_t len) {
  size_t i;
  for (i=0; i<len; ++i) {
    seq_nt256int8[i] = nt16_to_nt256int8_table[(uint8_t) (seq_nt16[i>>1] << ((i&1)<<2))];
  }
}

static inline int8_t*
nt16_decode_nt256int8(const uint8_t *seq_nt16, size_t len) {
  int8_t *seq_nt256int8 = (int8_t*) malloc(len);
  _nt16_decode_nt256int8(seq_nt256int8, seq_nt16, len);
  return seq_nt256int8;
}

static inline void
nt16_rev_ip(uint8_t *seq_nt16, size_t len) {
  size_t l=LEN_NT256_TO_NT16(len);
  bitarr_reverse_ip(seq_nt16, l);
  size_t i;
  for(i=0; i<l; i++) seq_nt16[i] = nt16_rev_table[seq_nt16[i]];
  bitarr_shiftl(seq_nt16, l, (l<<3)-(len<<2));
}

static inline void
_nt16_rev(uint8_t *seq_nt16_rev, const uint8_t *seq_nt16, size_t len) {
  size_t l=LEN_NT256_TO_NT16(len);
  size_t i;
  for (i=0; i<l; i++) seq_nt16_rev[i] = nt16_rev_table[seq_nt16[l-i-1]];
  bitarr_shiftl(seq_nt16_rev, l, (l<<3)-(len<<2));
}

static inline uint8_t*
nt16_rev(const uint8_t *seq_nt16, size_t len) {
  uint8_t *seq_nt16_rev = (uint8_t*) malloc(LEN_NT256_TO_NT16(len));
  _nt16_rev(seq_nt16_rev, seq_nt16, len);
  return seq_nt16_rev;
}


/******************
 * nt8 operations *
 ******************
 *
 * 'A' <=> 0x00 00000000
 * 'C' <=> 0x20 00100000
 * 'G' <=> 0xC0 11000000
 * 'T' <=> 0xE0 11100000
 * 'N' => 0xA0  10100000
 */

static inline void
_nt256char_encode_nt8(uint8_t *seq_nt8, const char *seq_nt256char, size_t len) {
  const char *p;
  size_t l = LEN_NT256_TO_NT8(len);
  memset(seq_nt8, 0, l);
  seq_nt8[0] = 0;
  for (p=seq_nt256char+len-1; len--; p--) {
    bitarr_shiftr(seq_nt8, l, 3);
    seq_nt8[0] |= nt256char_to_nt8_table[(unsigned char)*p];
  }
}

static inline uint8_t*
nt256char_encode_nt8(const char *seq_nt256char, size_t len) {
  uint8_t *seq_nt8 = malloc(LEN_NT256_TO_NT8(len));
  _nt256char_encode_nt8(seq_nt8, seq_nt256char, len);
  return seq_nt8;
}

#define cstr_encode_nt8(cstr) nt256char_encode_nt8((cstr), strlen(cstr))

static inline void
_nt8_decode_nt256char(char *seq_nt256char, const uint8_t *seq_nt8, size_t len) {
  size_t l = LEN_NT256_TO_NT8(len);
  uint8_t *tmp = (uint8_t*) malloc(l);
  memcpy(tmp, seq_nt8, l);
  
  size_t i;
  for (i=0; i<len; i++) {
    seq_nt256char[i] = nt8_to_nt256char_table[tmp[0]];
    bitarr_shiftl(tmp, l, 3);
  }
  free(tmp);
}

static inline char*
nt8_decode_nt256char(const uint8_t *seq_nt8, size_t len) {
  char *seq_nt256char = (char*) malloc(len);
  _nt8_decode_nt256char(seq_nt256char, seq_nt8, len);
  return seq_nt256char;
}

static inline char*
nt8_decode_cstr(const uint8_t *seq_nt8, size_t len) {
  char *seq_nt256char = (char*) malloc(len+1);
  _nt8_decode_nt256char(seq_nt256char, seq_nt8, len);
  seq_nt256char[len] = 0;
  return seq_nt256char;
}

static inline void
_nt8_rev(uint8_t *seq_nt8_rev, const uint8_t *seq_nt8, size_t len) {
  size_t l = LEN_NT256_TO_NT8(len);
  uint8_t *tmp = (uint8_t*) malloc(l);
  memcpy(tmp, seq_nt8, l);
  seq_nt8_rev[0] = 0;
  size_t i;
  for (i=0; i<len; i++) {
    bitarr_shiftr(seq_nt8_rev, l, 3);
    seq_nt8_rev[0] |= nt8_rev_table[tmp[0]];
    bitarr_shiftl(tmp, l, 3);
  }
  free(tmp);
}

static inline uint8_t*
nt8_rev(const uint8_t *seq_nt8, size_t len) {
  uint8_t *seq_nt8_rev = (uint8_t*) malloc(LEN_NT256_TO_NT8(len));
  _nt8_rev(seq_nt8_rev, seq_nt8, len);
  return seq_nt8_rev;
}

/************************
 * nt256int8 operations *
 ************************
 *
 * 'A' <=> 0x00 00000000
 * 'C' <=> 0x01 00000001
 * 'G' <=> 0x02 00000010
 * 'T' <=> 0x03 00000011
 */

static inline void
_nt256char_encode_nt256int8(int8_t *seq_nt256int8, const unsigned char *seq_nt256char, size_t len) {
  size_t i;
  for (i=0; i<len; i++) {
    seq_nt256int8[i] = nt256char_to_nt256int8_table[seq_nt256char[i]];
  }
}

static inline int8_t*
nt256char_encode_nt256int8(const unsigned char *seq_nt256char, size_t len) {
  int8_t *seq_nt256int8 = (int8_t*) malloc(len);
  _nt256char_encode_nt256int8(seq_nt256int8, seq_nt256char, len);
  return seq_nt256int8;
}

#define cstr_encode_nt256int8(cstr) nt256char_encode_nt256int8(((unsigned char*) cstr), strlen(cstr))

static inline void
_nt256int8_decode_nt256char(char *seq_nt256char, const int8_t *seq_nt256int8, size_t len) {
  size_t i;
  for(i=0; i<len; i++) {
    seq_nt256char[i] = nt256int8_to_nt256char_table[seq_nt256int8[i]];
  }
}

static inline char*
nt256int8_decode_nt256char(const int8_t *seq_nt256int8, size_t len) {
  char *seq_nt256char = (char*) malloc(len);
  _nt256int8_decode_nt256char(seq_nt256char, seq_nt256int8, len);
  return seq_nt256char;
}

static inline char*
nt256int8_decode_cstr(const int8_t *seq_nt256int8, size_t len) {
  char *seq_nt256char = (char*) malloc(len+1);
  _nt256int8_decode_nt256char(seq_nt256char, seq_nt256int8, len);
  seq_nt256char[len] = 0;
  return seq_nt256char;
}

static inline void
_nt256int8_rev(int8_t *seq_nt256int8_rev, const int8_t *seq_nt256int8,  size_t len) {

  for (seq_nt256int8_rev += len-1; len--; seq_nt256int8++, seq_nt256int8_rev--) {
    *seq_nt256int8_rev = 3 - *seq_nt256int8;
  }
}

static inline int8_t*
nt256int8_rev(const int8_t *seq_nt256int8,  size_t len) {
  
  int8_t *seq_nt256int8_rev = (int8_t*) malloc(len);
  _nt256int8_rev(seq_nt256int8_rev, seq_nt256int8, len);
  return seq_nt256int8_rev;

}

static inline void
nt256int8_rev_ip(int8_t *seq_nt256int8, size_t len) {
  bitarr_reverse_ip(seq_nt256int8, len);
  size_t i;
  for (i=0; i<len; ++i) seq_nt256int8[i] = 3 - seq_nt256int8[i];
}

/* #define CBITS 2 */
/* #define CLEN(len) LEN_NT256_TO_NT4(len) */
/* #define ENCODE(seq, len) nt256char_encode_nt4(seq, len) */
/* #define _ENCODE(code, seq, len) _nt256char_encode_nt4(code, seq, len) */
/* #define _REV(code_rev, code, len) _nt4_rev(code_rev, code, len) */
/* #define DECODE(code, len) nt4_decode_cstr(code, len) */

/* #define CBITS 3 */
/* #define CLEN(len) LEN_NT256_TO_NT8(len) */
/* #define ENCODE(seq, len) nt256char_encode_nt8(seq, len) */
/* #define _ENCODE(code, seq, len) _nt256char_encode_nt8(code, seq, len) */
/* #define _REV(code_rev, code, len) _nt8_rev(code_rev, code, len) */
/* #define DECODE(code, len) nt8_decode_cstr(code, len) */

#define CBITS 4
#define CLEN(len) LEN_NT256_TO_NT16(len)
#define ENCODE(seq, len) nt256char_encode_nt16(seq, len)
#define _ENCODE(code, seq, len) _nt256char_encode_nt16(code, seq, len)
#define _REV(code_rev, code, len) _nt16_rev(code_rev, code, len)
#define REV(code, len) nt16_rev(code, len)
#define DECODE(code, len) nt16_decode_cstr(code, len)
#define _DECODE(code_rev, code, len) _nt16_decode_nt256char(code_rev, code, len)
#define TRCODE_NT256INT8(seq, len) nt16_decode_nt256int8(seq, len);


#endif
