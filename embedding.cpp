#include "header.h"

#define EMBED_PAD 4
#define CGK2_EMBED 1

int Embedding::cgk2_embed(const char ** oridata, unsigned rlen, int threshold, int id,
        int strid, char *embeddedQ)
{
    int nmismatch = 0;
    int elen = efactor * rlen;
    assert(elen <= MAX_ELEN);

    int j = 0;
    if (id == 0) {
        for (unsigned i = 0; i < rlen; i++) {
            uint8_t s = oridata[id][i];
            char bit = hash_eb[BITPOS(strid, j, s)];
            if (!bit) {
                // here, jth embedded value should be s. for query (id = 0) we need
                // to generate the embedding. for id > 1, we need to verify.
                embeddedQ[j] = s;
                j++;
            } else {
                // here, jth and j+1th value are both s.
                embeddedQ[j + 1] = embeddedQ[j] = s;
                j += 2;
            }
        }

        //append the rest with EMBED_PAD
        //because the embedded candidate may be longer than j and need to count nmismatch with embeddedQ[j]
        for(; j < elen; j++){
            embeddedQ[j] = EMBED_PAD;
        }

    } else {
        for (unsigned i = 0; i < rlen; i++) {
            uint8_t s = oridata[id][i];
            char bit = hash_eb[BITPOS(strid, j, s)];
            if (!bit) {
                // here, jth embedded value should be s. for query (id = 0) we need
                // to generate the embedding. for id > 1, we need to verify.
                nmismatch += (embeddedQ[j] == s ? 0 : 1);

                if (nmismatch > threshold) {
                    nmismatch = elen;
                    goto end;
                }
                j++;
            } else {
                // here, jth and j+1th value are both s.
                nmismatch += (embeddedQ[j] == s ? 0 : 1);
                nmismatch += (embeddedQ[j + 1] == s ? 0 : 1);
                if (nmismatch > threshold) {
                    goto end;
                }

                j += 2;
            }
        }
    }

    assert(j <= elen);

    end:

    return nmismatch;
}

int Embedding::cgk_embed(const char ** oridata, unsigned rlen, int threshold, int id,
        int strid, char *embeddedQ)
{
    unsigned i = 0;
    int nmismatch = 0;
    int elen = 3 * rlen;
    assert(elen < MAX_ELEN);

    if (id == 0) {
        for (int j = 0; j < elen; j++)
        {
            // we either use original input or we pad
            uint8_t s = i < rlen ? oridata[id][i] : EMBED_PAD;

            // id 0 is reserved for the query itself. so we do not check it.
            embeddedQ[j] = s;
            char bit = hash_eb[BITPOS(strid, j, s)];
            i = i + (bit >= 1 ? 0 : 1);
        }
    } else {
        for (int j = 0; j < elen; j++)
        {
            // we either use original input or we pad
            uint8_t s = i < rlen ? oridata[id][i] : EMBED_PAD;
            nmismatch += (embeddedQ[j] == s ? 0 : 1);

            if (nmismatch > threshold) {
                nmismatch = elen;
                goto end;
            }
            char bit = hash_eb[BITPOS(strid, j, s)];
            i = i + (bit >= 1 ? 0 : 1);
        }
    }

end:

    return nmismatch;
}

int Embedding::embedstr(const char ** oridata, unsigned rlen, int threshold, int id,
        int strid, char *embeddedQ)
{
#ifdef CGK2_EMBED
    return cgk2_embed(oridata, rlen, threshold, id, strid, embeddedQ);
#else
    return cgk_embed(oridata, rlen, threshold, id, strid, embeddedQ);
#endif
}

void Embedding::embed_two_pairs(vector<Region> &candidate_regions_f1, vector<Region> &candidate_regions_r2,
                                unsigned rlen, int* thresholds, char embeddedQ_f1[], char embeddedQ_r2[],
                                const char **candidate_refs_f1, const char **candidate_refs_r2,
                                unsigned nregions_f1, unsigned nregions_r2, bool flag_f1[], bool flag_r2[],
                                int best_f1, int next_f1, int best_r2, int next_r2){

    if (best_f1 == nregions_f1 || best_r2 == nregions_r2) // not equal to the init value
        return;

    int elen = rlen * efactor;

    // embed first pair, threshold is elen (embed whole read)
    int nmismatch_f1 = embedstr(candidate_refs_f1, rlen, elen, best_f1+1, 0, embeddedQ_f1);
    candidate_regions_f1[best_f1].embed_dist = nmismatch_f1;
    int nmismatch_r2 = embedstr(candidate_refs_r2, rlen, elen, best_r2+1, 0, embeddedQ_r2);
    candidate_regions_r2[best_r2].embed_dist = nmismatch_r2;
    thresholds[0] = nmismatch_f1 + nmismatch_r2;
    // set flag to 0, no need to embed again later
    flag_f1[best_f1] = 0;
    flag_r2[best_r2] = 0;

    // embed the second pair
    if (next_f1 != nregions_f1 && next_r2 != nregions_r2) { // not equal to the init value
        int nmismatch;
        if (next_f1 == best_f1){
            assert(next_r2 != best_r2);
            nmismatch = embedstr(candidate_refs_r2, rlen, elen, next_r2+1, 0, embeddedQ_r2);
            candidate_regions_r2[next_r2].embed_dist = nmismatch;
            thresholds[1] = nmismatch_f1 + nmismatch;
        }else if (next_r2 == best_r2){
            assert(next_f1 != best_f1);
            nmismatch = embedstr(candidate_refs_f1, rlen, elen, next_f1+1, 0, embeddedQ_f1);
            candidate_regions_f1[next_f1].embed_dist = nmismatch;
            thresholds[1] = nmismatch_r2 + nmismatch;
        }else{
            nmismatch = embedstr(candidate_refs_f1, rlen, elen, next_f1+1, 0, embeddedQ_f1);
            candidate_regions_f1[next_f1].embed_dist = nmismatch;
            thresholds[1] = nmismatch;
            nmismatch = embedstr(candidate_refs_r2, rlen, elen, next_r2+1, 0, embeddedQ_r2);
            candidate_regions_r2[next_r2].embed_dist = nmismatch;
            thresholds[1] += nmismatch;
        }
        flag_f1[next_f1] = 0;
        flag_r2[next_r2] = 0;
    }

}

void Embedding::embeddata_pair(vector<Region> &candidate_regions, char embeddedQ[],
                          const char **candidate_refs, unsigned ncandidates, bool flag[],
                          unsigned rlen, int threshold)
{
    struct timeval start, end;
    gettimeofday(&start, NULL);

    int step = 1;

    for (unsigned i = 0; i < ncandidates; i += step) {
        if (!flag[i])
            continue;
        candidate_regions[i].embed_dist = embedstr(candidate_refs, rlen, threshold, i + 1, 0, embeddedQ);
    }

    gettimeofday(&end, NULL);
    embed_time += compute_elapsed(&start, &end);
}

void Embedding::embeddata_iterative_update(vector<Region> &candidate_regions,
                          const char **input, unsigned ninput, unsigned rlen,
                          int &best_threshold, int &next_threshold,
                          bool max_rnd, int &best_idx, int &next_idx)
{
    int step = 1;
    int elen = rlen * efactor;
    char embeddedQ[elen];

    // do first candidate
    struct timeval start, end;
    gettimeofday(&start, NULL);
    int nmismatch = 0;
    if (best_threshold == 0 && next_threshold == 0) {
        // if we already have 2 exact match (one for best, one for second best for mapq), look for exact matches only
        nmismatch = (memcmp(input[1], input[0], rlen) == 0 ? 0 : elen);
    } else {
        // embed Q
        embedstr(input, rlen, elen, 0, 0, embeddedQ);
        nmismatch = embedstr(input, rlen, next_threshold, 1, 0, embeddedQ);
    }
    candidate_regions[0].embed_dist = nmismatch;
    if (nmismatch < best_threshold) {
        next_threshold = best_threshold;
        best_threshold = nmismatch;
    } else if (nmismatch < next_threshold) {
        next_threshold = nmismatch;
    }
    int best_dist = nmismatch;
    int next_dist = numeric_limits<int>::max();
    best_idx = next_idx = 0;

//#if DBGPRINT
//    cout << "Candidate region 0 at pos " << candidate_regions[0].pos <<
//        " with cov " << candidate_regions[0].cov << " has nmismatch " <<
//        nmismatch << endl;
//#endif

    for (unsigned i = 2; i < ninput; i += step) {
        if (best_threshold == 0 && next_threshold == 0 ) {
            // if we already have 2 exact match (one for best, one for second best for mapq), look for exact matches only
            nmismatch = (memcmp(input[i], input[0], rlen) == 0 ? 0 : elen);
        } else {
            nmismatch = embedstr(input, rlen, next_threshold, i, 0, embeddedQ);
        }
        candidate_regions[i - 1].embed_dist = nmismatch;

//#if DBGPRINT
//        cout << "Candidate region " << i - 1 << " at pos " << candidate_regions[i - 1].pos <<
//            " with cov " << candidate_regions[0].cov << " has nmismatch " <<
//            nmismatch << endl;
//#endif

        if (nmismatch < best_threshold) {
            next_threshold = best_threshold;
            best_threshold = nmismatch;
        } else if (nmismatch < next_threshold) {
            next_threshold = nmismatch;
        }

        // set best and next idx so that we dont have to sort regions later
        // pick the minimal pos hit when several hits have same embed dist and the best one is not the hcov(0)
        if (nmismatch < best_dist ||
            (nmismatch == best_dist && best_idx!= 0 && candidate_regions[i-1].pos < candidate_regions[best_idx].pos)) {
            next_dist = best_dist;
            next_idx = best_idx;
            best_dist = nmismatch;
            best_idx = i - 1;
        } else if (nmismatch < next_dist) {
            next_dist = nmismatch;
            next_idx = i - 1;
        }
    }

    gettimeofday(&end, NULL);
    embed_time += compute_elapsed(&start, &end);
}

Embedding::Embedding()
{
#ifdef CGK2_EMBED
    efactor = 2;
#else
    efactor = 3;
#endif

    time_t seed = time(NULL);
    //srand(seed);
    cerr << "Embedding using random seed " << seed << endl;
    srand(1559063236);
    cerr << "Creating " << NUM_STR << " random string(s) of length: " <<
        MAX_ELEN << endl;
    embed_time = 0;

    /* initialize the hash structures. These are our random strings.
     * So we need to generate randomstring of length MAX_ELEN.
     * For each position in this string, we need one random value for
     * each of the characters in our alphabet A,C,G,T,N
     * hash_eb[x] loops over each random string. Within that, we first
     * store 0th bit of all chars first, then 1st bit of all chars next
     * in a flattended 2d array. So hash_eb[0][0-NUM_CHAR] is 0th bit
     * of all NUM_CHARS, [0[NUM_CHAR -- 2*NUM_CHAR] is bit 1 for all
     * NUM_CHARS and so on.
     */

    for (int j = 0; j < NUM_STR; j++)
        for (int t = 0; t < NUM_CHAR; t++)
            for (int d = 0; d < MAX_ELEN; d++)
                hash_eb[BITPOS(j, d, t)] = 1 - rand() % 2;
}

Embedding::Embedding(const char *fname)
{
#ifdef CGK2_EMBED
    efactor = 2;
#else
    efactor = 3;
#endif
    cerr << "Loading embedding from file " << fname << endl;
    ifstream input(fname);
    if (input) {
        char num_str, num_char;
        input.read((char *)&num_str, sizeof(int));
        input.read((char *)&num_char, sizeof(int));
        assert(num_str == NUM_STR && num_char == NUM_CHAR);

        cerr << "NUM_STR: " << NUM_STR << ", NUM_CHAR: " << NUM_CHAR <<
            " ,MAX_ELEN: " << MAX_ELEN << endl;

         // initialize the hash structures
        string str_hasheb;
        std::getline(input, str_hasheb);
        for (int i = 0; i < TOTAL_RBITS; i++)
            hash_eb[i] = str_hasheb[i];
    }
}

Embedding::~Embedding()
{
}

