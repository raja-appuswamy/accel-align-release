#pragma once

#define MOD ((1UL<<29)-1)
#define MOD ((1UL<<29)-1)
#define BATCH_SIZE 15000
//percent of the read length
#define MAX_INDEL 15
#define SC_MCH 2
#define SC_MIS 8
#define GAPO 12
#define GAPE 2
#define SC_AMBI 1 // score when one or both bases are "N"
#define END_BONUS 10
#define Z_DROP 100
#define BANDWIDTH 151

//modes
#define SHORT_READ_MODE 0

//Embedding related
#define EMBED_PAD 4
#define NUM_STR 3 //r
#define NUM_CHAR 5 // num of chars in strings
#define MAX_ELEN 1500 //the maximum length of random string
#define RBITS_PER_STRING (MAX_ELEN * NUM_CHAR)
#define TOTAL_RBITS (RBITS_PER_STRING * NUM_STR)
#define BITPOS(STR_ID, OFFSET, CHAR_ID) (STR_ID * RBITS_PER_STRING +\
        OFFSET * NUM_CHAR + CHAR_ID)

#define MAX_LEN 512

#define KSW_EZ_SCORE_ONLY  0x01 // don't record alignment path/cigar
#define KSW_EZ_RIGHT       0x02 // right-align gaps
#define KSW_EZ_GENERIC_SC  0x04 // without this flag: match/mismatch only; last symbol is a wildcard
#define KSW_EZ_APPROX_MAX  0x08 // approximate max; this is faster with sse
#define KSW_EZ_APPROX_DROP 0x10 // approximate Z-drop; faster with sse
#define KSW_EZ_EXTZ_ONLY   0x40 // only perform extension
#define KSW_EZ_REV_CIGAR   0x80 // reverse CIGAR in the output

#ifndef kroundup32
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif

#define MAX_OCC 1000