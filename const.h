#pragma once

#define MOD ((1UL<<29)-1)
#define MOD ((1UL<<29)-1)
#define BATCH_SIZE 15000
//percent of the read length
#define MAX_INDEL 5

//modes
#define SHORT_READ_MODE 0

//Embedding related
#define EMBED_PAD 4
#define NUM_STR 1 //r
#define NUM_CHAR 5 // num of chars in strings
#define MAX_ELEN 1500 //the maximum length of random string
#define RBITS_PER_STRING (MAX_ELEN * NUM_CHAR)
#define TOTAL_RBITS (RBITS_PER_STRING * NUM_STR)
#define BITPOS(STR_ID, OFFSET, CHAR_ID) (STR_ID * RBITS_PER_STRING +\
        OFFSET * NUM_CHAR + CHAR_ID)

// mapq
#define MIS_PENALTY -1

#define MAX_LEN 512