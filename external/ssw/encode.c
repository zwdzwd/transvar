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
 */
#include <stdio.h>
#include "encode.h"


const uint8_t nt256char_to_nt4_table[128] = {
  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,
  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,
  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,
  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,
  1,  0,  1, 64,  1,  1,  1,128,  1,  1,  1,  1,  1,  1,  1,  1,
  1,  1,  1,  1,192,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,
  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,
  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1
};


const char nt4_to_nt256char_table[256] = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT";


const int8_t nt4_to_nt256int8_table[256] = {
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
  1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
  1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
  1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
  2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,
  2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,
  2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,
  2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,
  3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,
  3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,
  3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,
  3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3
};

const uint8_t nt4_to_nt16_table[256] = {
  16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16,
  16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16,
  16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16,
  16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16,
  32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32,
  32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32,
  32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32,
  32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32,
  64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64,
  64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64,
  64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64,
  64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64,
  128,128,128,128,128,128,128,128,128,128,128,128,128,128,128,128,
  128,128,128,128,128,128,128,128,128,128,128,128,128,128,128,128,
  128,128,128,128,128,128,128,128,128,128,128,128,128,128,128,128,
  128,128,128,128,128,128,128,128,128,128,128,128,128,128,128,128
};

const uint8_t nt4_rev_table[256] = {
  255,191,127, 63,239,175,111, 47,223,159, 95, 31,207,143, 79, 15,
  251,187,123, 59,235,171,107, 43,219,155, 91, 27,203,139, 75, 11,
  247,183,119, 55,231,167,103, 39,215,151, 87, 23,199,135, 71,  7,
  243,179,115, 51,227,163, 99, 35,211,147, 83, 19,195,131, 67,  3,
  254,190,126, 62,238,174,110, 46,222,158, 94, 30,206,142, 78, 14,
  250,186,122, 58,234,170,106, 42,218,154, 90, 26,202,138, 74, 10,
  246,182,118, 54,230,166,102, 38,214,150, 86, 22,198,134, 70,  6,
  242,178,114, 50,226,162, 98, 34,210,146, 82, 18,194,130, 66,  2,
  253,189,125, 61,237,173,109, 45,221,157, 93, 29,205,141, 77, 13,
  249,185,121, 57,233,169,105, 41,217,153, 89, 25,201,137, 73,  9,
  245,181,117, 53,229,165,101, 37,213,149, 85, 21,197,133, 69,  5,
  241,177,113, 49,225,161, 97, 33,209,145, 81, 17,193,129, 65,  1,
  252,188,124, 60,236,172,108, 44,220,156, 92, 28,204,140, 76, 12,
  248,184,120, 56,232,168,104, 40,216,152, 88, 24,200,136, 72,  8,
  244,180,116, 52,228,164,100, 36,212,148, 84, 20,196,132, 68,  4,
  240,176,112, 48,224,160, 96, 32,208,144, 80, 16,192,128, 64,  0
};


/* similar to bam_nt16_table from samtools' bam_import.c
 * except that here uses upper four bits.
 * 
 * Note that this table contains several conversion
 * '1','2','3','4' => 1,2,4,8
 * 'A','C','G','T', 'N' => 1,2,4,8,15
 * */
const uint8_t nt256char_to_nt16_table[128] = {
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0, 16, 32, 64,128,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0, 16,  0, 32,  0,  0,  0, 64,  0,  0,  0,  0,  0,  0,240,  0,
  0,  0,  0,  0,128,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0
};

const char nt16_to_nt256char_table[256]="================AAAAAAAAAAAAAAAACCCCCCCCCCCCCCCC================GGGGGGGGGGGGGGGG================================================TTTTTTTTTTTTTTTT================================================================================================NNNNNNNNNNNNNNNN";

const int8_t nt16_to_nt256int8_table[256]={
  4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
  4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
  2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,
  4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
  4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
  4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
  3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,
  4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
  4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
  4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
  4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
  4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
  4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
  4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4
};

const uint8_t nt16_rev_table[256] = {
   0,128, 64,192, 32,160, 96,224, 16,144, 80,208, 48,176,112,240,
   8,136, 72,200, 40,168,104,232, 24,152, 88,216, 56,184,120,248,
   4,132, 68,196, 36,164,100,228, 20,148, 84,212, 52,180,116,244,
  12,140, 76,204, 44,172,108,236, 28,156, 92,220, 60,188,124,252,
   2,130, 66,194, 34,162, 98,226, 18,146, 82,210, 50,178,114,242,
  10,138, 74,202, 42,170,106,234, 26,154, 90,218, 58,186,122,250,
   6,134, 70,198, 38,166,102,230, 22,150, 86,214, 54,182,118,246,
  14,142, 78,206, 46,174,110,238, 30,158, 94,222, 62,190,126,254,
   1,129, 65,193, 33,161, 97,225, 17,145, 81,209, 49,177,113,241,
   9,137, 73,201, 41,169,105,233, 25,153, 89,217, 57,185,121,249,
   5,133, 69,197, 37,165,101,229, 21,149, 85,213, 53,181,117,245,
  13,141, 77,205, 45,173,109,237, 29,157, 93,221, 61,189,125,253,
   3,131, 67,195, 35,163, 99,227, 19,147, 83,211, 51,179,115,243,
  11,139, 75,203, 43,171,107,235, 27,155, 91,219, 59,187,123,251,
   7,135, 71,199, 39,167,103,231, 23,151, 87,215, 55,183,119,247,
  15,143, 79,207, 47,175,111,239, 31,159, 95,223, 63,191,127,255, 
};

const char nt256int8_to_nt256char_table[4] = "ACGT";

const int8_t nt256char_to_nt256int8_table[128] = {
  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
  4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
  4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
  4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
  4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4 
};

const uint8_t nt256char_to_nt8_table[256] = {
  160,160,160,160,160,160,160,160,160,160,160,160,160,160,160,160,
  160,160,160,160,160,160,160,160,160,160,160,160,160,160,160,160,
  160,160,160,160,160,160,160,160,160,160,160,160,160,160,160,160,
  160,160,160,160,160,160,160,160,160,160,160,160,160,160,160,160,
  160,  0,160, 32,160,160,160,192,160,160,160,160,160,160,160,160,
  160,160,160,160,224,160,160,160,160,160,160,160,160,160,160,160,
  160,160,160,160,160,160,160,160,160,160,160,160,160,160,160,160,
  160,160,160,160,160,160,160,160,160,160,160,160,160,160,160,160,
  160,160,160,160,160,160,160,160,160,160,160,160,160,160,160,160,
  160,160,160,160,160,160,160,160,160,160,160,160,160,160,160,160,
  160,160,160,160,160,160,160,160,160,160,160,160,160,160,160,160,
  160,160,160,160,160,160,160,160,160,160,160,160,160,160,160,160,
  160,160,160,160,160,160,160,160,160,160,160,160,160,160,160,160,
  160,160,160,160,160,160,160,160,160,160,160,160,160,160,160,160,
  160,160,160,160,160,160,160,160,160,160,160,160,160,160,160,160,
  160,160,160,160,160,160,160,160,160,160,160,160,160,160,160,160
};

const char nt8_to_nt256char_table[256] = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT";

const uint8_t nt8_rev_table[256] = {
  224,224,224,224,224,224,224,224,224,224,224,224,224,224,224,224,
  224,224,224,224,224,224,224,224,224,224,224,224,224,224,224,224,
  192,192,192,192,192,192,192,192,192,192,192,192,192,192,192,192,
  192,192,192,192,192,192,192,192,192,192,192,192,192,192,192,192,
  160,160,160,160,160,160,160,160,160,160,160,160,160,160,160,160,
  160,160,160,160,160,160,160,160,160,160,160,160,160,160,160,160,
  128,128,128,128,128,128,128,128,128,128,128,128,128,128,128,128,
  128,128,128,128,128,128,128,128,128,128,128,128,128,128,128,128,
  96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96,
  96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96,
  64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64,
  64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64,
  32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32,
  32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0
};

const char nt256char_rev_table[128] = "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTNGNNNCNNNNNNNNNNNNANNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN";

#ifdef MAIN_PRINT_CODING_TABLE
/* gcc encode.c -DMAIN_PRINT_CODING_TABLE */
int main() {
  int i;
  uint8_t j;
  printf("nt256char => nt4\n");
  for (i=0; i<128; i++) {
    if (i == 'A') printf("%3d,", 0);
    else if (i == 'C') printf("%3d,", 0x40);
    else if (i == 'G') printf("%3d,", 0x80);
    else if (i == 'T') printf("%3d,", 0xC0);
    else printf("%3d,", 1);
    if (((i+1) & 0xf) == 0) printf("\n");
  }

  printf("\n\nnt4 => nt256char\n");
  for (i=0; i<256; i++) {
    if ((i & 0xC0) == 0) putchar('A');
    else if ((i & 0xC0) == 0x40) putchar('C');
    else if ((i & 0xC0) == 0x80) putchar('G');
    else if ((i & 0xC0) == 0xC0) putchar('T');
    else putchar('N');
  }

  printf("\n\nnt4 => nt256int8\n");
  for (i=0; i<256; i++) {
    if ((i & 0xC0) == 0) printf("%d,", 0);
    else if ((i & 0xC0) == 0x40) printf("%d,", 1);
    else if ((i & 0xC0) == 0x80) printf("%d,", 2);
    else if ((i & 0xC0) == 0xC0) printf("%d,", 3);
    else printf("%d,", 4);
    if (((i+1) & 0xf) == 0) printf("\n");
  }

  puts("\n\nnt4 rev");
  for (i=0; i<256; i++) {
    j=i;
    j=((j & 0xF0)>>4 | (j & 0x0F)<<4);
    j=((j & 0xCC)>>2 | (j & 0x33)<<2);
    printf("%3u,", (uint8_t) ~j);
    if (((i+1) & 0xf) == 0) printf("\n");
  }

  printf("\n\nnt4 => nt16\n");
  for (i=0; i<256; i++) {
    if ((i & 0xC0) == 0) printf("%3d,", 0x10);
    else if ((i & 0xC0) == 0x40) printf("%3d,", 0x20);
    else if ((i & 0xC0) == 0x80) printf("%3d,", 0x40);
    else if ((i & 0xC0) == 0xC0) printf("%3d,", 0x80);
    else printf("%d,", 0);
    if (((i+1) & 0xf) == 0) printf("\n");
  }

  printf("\n\nnt256char => nt16\n");
  for (i=0; i<128; i++) {
    if (i == 'A') printf("%3d,", 0x10);
    else if (i == 'C') printf("%3d,", 0x20);
    else if (i == 'G') printf("%3d,", 0x40);
    else if (i == 'T') printf("%3d,", 0x80);
    else if (i == 'N') printf("%3d,", 0xF0);
    else if (i == '1') printf("%3d,", 0x10);
    else if (i == '2') printf("%3d,", 0x20);
    else if (i == '3') printf("%3d,", 0x40);
    else if (i == '4') printf("%3d,", 0x80);
    else printf("%3d,", 0);
    if (((i+1) & 0xf) == 0) printf("\n");
  }

  puts("\n\nnt16=>nt256char");
  for (i=0; i<256; i++) {
    if ((i & 0xF0) == 0x10) putchar('A');
    else if ((i & 0xF0) == 0x20) putchar('C');
    else if ((i & 0xF0) == 0x40) putchar('G');
    else if ((i & 0xF0) == 0x80) putchar('T');
    else if ((i & 0xF0) == 0xF0) putchar('N');
    else putchar('=');
  }
  puts("\n");

  puts("\n\nnt16 => nt256int");
  for (i=0; i<256; i++) {
    if ((i & 0xF0) == 0x10) printf("0,");
    else if ((i & 0xF0) == 0x20) printf("1,");
    else if ((i & 0xF0) == 0x40) printf("2,");
    else if ((i & 0xF0) == 0x80) printf("3,");
    else printf("4,");
    if (((i+1) & 0xf) == 0) printf("\n");
  }


  puts("\n\nnt16 rev");
  for (i=0; i<256; i++) {
    j=i;
    j=((j & 0xF0)>>4 | (j & 0x0F)<<4);
    j=((j & 0xCC)>>2 | (j & 0x33)<<2);
    j=((j & 0xAA)>>1 | (j & 0x55)<<1);
    printf("%3d,", j);
    if (((i+1) & 0xf) == 0) printf("\n");
  }

  puts("\n\nnt256car => nt8");
  for (i=0; i<256; i++) {
    if (i == 'A') printf("%3d,", 0x00);
    else if (i == 'C') printf("%3d,", 0x20);
    else if (i == 'G') printf("%3d,", 0xC0);
    else if (i == 'T') printf("%3d,", 0xE0);
    else printf("%3d,", 0xA0);
    if (((i+1) & 0xf) == 0) printf("\n");
  }

  puts("\n\nnt8=>nt256char");
  for (i=0; i<256; i++) {
    if ((i & 0xE0) == 0x00) putchar('A');
    else if ((i & 0xE0) == 0x20) putchar('C');
    else if ((i & 0xE0) == 0xC0) putchar('G');
    else if ((i & 0xE0) == 0xE0) putchar('T');
    else putchar('N');
  }
  puts("\n");

  puts("\n\nnt8 rev");
  for (i=0; i<256; i++) {
    if ((i & 0xE0) == 0x00) printf("%3d,", 0xE0);
    else if ((i & 0xE0) == 0x20) printf("%3d,", 0xC0);
    else if ((i & 0xE0) == 0xC0) printf("%3d,", 0x20);
    else if ((i & 0xE0) == 0xE0) printf("%3d,", 0x00);
    else  printf("%3d,", (~i & 0xE0));
    if (((i+1) & 0xf) == 0) printf("\n");
  }
  puts("\n");

  puts("\n\nnt256char rev");
  for (i=0; i<128; i++) {
    if (i=='A') putchar('T');
    else if (i=='T') putchar('A');
    else if (i=='C') putchar('G');
    else if (i=='G') putchar('C');
    else putchar('N');
  }
  puts("\n");
  

  return 0;
}

#endif

#ifdef MAIN_ENCODE

int main() {
  char *seq_nt256char = "AAAAGGAACTTAGTGAAATTTTAGGTAGAAAGTTGACAGGGTACCAGGAGATGATGTAAGGGACAAGCAGCCACACCCCATTCTTGAGGGGCTGAGGTGGAAGAGACAGGCCCGGAGGGGTGAGGCAGTCTTTACTCACCTGTAGATGTCTCGGGCCATCCCGAAGTCTCCAATCTTGGCCACTCTTCCAGGGCCTGGACAGGTCAAGAGGCAGTTTCTGGCAGCAATGTCTCTGGGAAGAAAGGAAATGCATTTCCTAATTTTATCCCTAGGAAGATGAGTGTACAACGGCCATCACTAG";
  /* char *seq_nt256char = "ATGCATGCAAG"; */
  size_t l = strlen(seq_nt256char);
  printf("seq_nt256char: %s\n", seq_nt256char);

  puts("\n\n==== test nt4 ====");
  uint8_t *seq_nt4 = cstr_encode_nt4(seq_nt256char);
  puts("encode to nt4:");
  bitarr_show(seq_nt4, LEN_NT256_TO_NT4(l));
  printf("seq_nt256char (decoded): %s\n", nt4_decode_cstr(seq_nt4, l));
  printf("seq_nt256char rev: %s\n", nt4_decode_cstr(nt4_rev(seq_nt4, l), l));
  nt4_rev_ip(seq_nt4, l);
  printf("seq_nt256char rev in-place: %s\n", nt4_decode_nt256char(seq_nt4, l));
  printf("transcode to nt16:\n");
  uint8_t *seq_nt16tr=nt4_decode_nt16(seq_nt4, l);
  bitarr_show(seq_nt16tr, LEN_NT256_TO_NT16(l));
  printf("decoded %s\n", nt16_decode_cstr(seq_nt16tr, l));

  puts("\n\n==== test nt16 =====");
  uint8_t *seq_nt16 = cstr_encode_nt16(seq_nt256char);
  puts("encode to nt16:");
  bitarr_show(seq_nt16, LEN_NT256_TO_NT16(l));
  printf("seq_nt256char (decoded): %s\n", nt16_decode_cstr(seq_nt16, l));
  puts("reverse nt16:");
  bitarr_show(nt16_rev(seq_nt16, l), LEN_NT256_TO_NT16(l));
  printf("reversed decoded: %s\n", nt16_decode_cstr(nt16_rev(seq_nt16, l), l));

  puts("\n\n===== test nt256int8 ====");
  int8_t *seq_nt256int8 = cstr_encode_nt256int8(seq_nt256char);
  puts("encode to nt256int8 and decode back:");
  puts(nt256int8_decode_cstr(seq_nt256int8, l));
  puts("reversed");
  puts(nt256int8_decode_cstr(nt256int8_rev(seq_nt256int8, l), l));
  puts("transcode from nt16:");
  puts(nt16_decode_cstr(seq_nt16, l));
  seq_nt256int8 = nt16_decode_nt256int8(seq_nt16,l);
  puts(nt256int8_decode_nt256char(seq_nt256int8, l));

  puts("\n\n===== test nt8 ====");
  int8_t *seq_nt8 = cstr_encode_nt8(seq_nt256char);
  bitarr_show(seq_nt8, LEN_NT256_TO_NT8(l));
  printf("decoded: %s\n", nt8_decode_cstr(seq_nt8, l));
  printf("reversed: %s\n", nt8_decode_cstr(nt8_rev(seq_nt8, l), l));
  
  return 0;
}

#endif
