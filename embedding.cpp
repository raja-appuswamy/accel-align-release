#include "header.h"

#define EMBED_PAD 4
#define CGK2_EMBED 1

int Embedding::cgk2_unmatched(const char *r, const char *ref,
                              const vector<uint32_t> &mch,
                              const unsigned rlen,
                              const unsigned kmer_step,
                              const int threshold,
                              const int strid) {
  unsigned elen = efactor * rlen;
  int nmismatch = 0;
  char embeddedQ[elen];

  // embed ref: cadidate pos
  unsigned i = 0, m_idx = 0, j = 0;
  while (i < rlen) {
    if (m_idx < mch.size() && i == mch[m_idx]) {
      i += kmer_step;
      m_idx++;
    } else {
      uint8_t s = ref[i];
      char bit = hash_eb[BITPOS(strid, j, s)];
      if (!bit) {
        embeddedQ[j] = s;
        j++;
      } else {
        embeddedQ[j + 1] = embeddedQ[j] = s;
        j += 2;
      }
      i++;
    }
  }

  //append the rest with EMBED_PAD
  //because the embedded candidate may be longer than j and need to count nmismatch with embeddedQ[j]
  for (; j < elen; j++) {
    embeddedQ[j] = EMBED_PAD;
  }

  //reset idx, embed read and cal mismatch
  i = j = m_idx = 0;
  while (i < rlen) {
    if (m_idx < mch.size() && i == mch[m_idx]) {
      i += kmer_step;
      m_idx++;
    } else {
      uint8_t s = r[i];
      char bit = hash_eb[BITPOS(strid, j, s)];
      if (!bit) {
        nmismatch += (embeddedQ[j] == s ? 0 : 1);
        j++;
      } else {
        nmismatch += (embeddedQ[j] == s ? 0 : 1);
        nmismatch += (embeddedQ[j + 1] == s ? 0 : 1);
        j += 2;
      }
      if (nmismatch > threshold) {
        nmismatch = elen;
        goto end;
      }
      i++;
    }
  }

  for (; j < elen; j++)
    nmismatch += (embeddedQ[j] == EMBED_PAD ? 0 : 1);

  assert(j == elen);
  end:
  return nmismatch;
}

void Embedding::cgk2_embedQ(const char *oridata, unsigned rlen, int strid, char *embeddedQ) {
  int j = 0;
  int elen = efactor * rlen;
  for (unsigned i = 0; i < rlen; i++) {
    uint8_t s = oridata[i];
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
  for (; j < elen; j++) {
    embeddedQ[j] = EMBED_PAD;
  }
  assert(j <= elen);
}

int Embedding::cgk2_embed_nmismatch(const char *oridata, unsigned rlen, int threshold, int strid, char *embeddedQ) {
  int nmismatch = 0;
  int j = 0;
  int elen = efactor * rlen;

  for (unsigned i = 0; i < rlen; i++) {
    uint8_t s = oridata[i];
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

  for (; j < elen; j++)
    nmismatch += (embeddedQ[j] == EMBED_PAD ? 0 : 1);

  assert(j == elen);
  end:
  return nmismatch;
}

int Embedding::cgk2_embed(const char **oridata, unsigned rlen, int threshold, int id,
                          int strid, char *embeddedQ) {
  int nmismatch = 0;

  if (id == 0) {
    cgk2_embedQ(oridata[id], rlen, strid, embeddedQ);
  } else {
    nmismatch = cgk2_embed_nmismatch(oridata[id], rlen, threshold, strid, embeddedQ);
  }

  return nmismatch;
}

int Embedding::cgk_embed(const char **oridata, unsigned rlen, int threshold, int id,
                         int strid, char *embeddedQ) {
  unsigned i = 0;
  int nmismatch = 0;
  int elen = 3 * rlen;
  assert(elen < MAX_ELEN);

  if (id == 0) {
    for (int j = 0; j < elen; j++) {
      // we either use original input or we pad
      uint8_t s = i < rlen ? oridata[id][i] : EMBED_PAD;

      // id 0 is reserved for the query itself. so we do not check it.
      embeddedQ[j] = s;
      char bit = hash_eb[BITPOS(strid, j, s)];
      i = i + (bit >= 1 ? 0 : 1);
    }
  } else {
    for (int j = 0; j < elen; j++) {
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

int Embedding::embedstr(const char **oridata, unsigned rlen, int threshold, int id,
                        int strid, char *embeddedQ) {
#ifdef CGK2_EMBED
  return cgk2_embed(oridata, rlen, threshold, id, strid, embeddedQ);
#else
  return cgk_embed(oridata, rlen, threshold, id, strid, embeddedQ);
#endif
}

void Embedding::embed_unmatch_iter(vector<Region> &candidate_regions, const char *ptr_ref, const char *r,
                                   const unsigned rlen, const unsigned kmer_step, int &best_threshold,
                                   int &next_threshold, unsigned &best_idx, unsigned &next_idx) {

  auto start = std::chrono::system_clock::now();

  int elen = rlen * efactor, nmismatch;

  for (unsigned i = 0; i < candidate_regions.size(); ++i) {
    Region &region = candidate_regions[i];
    region.embed_dist = elen;

    for (unsigned strid = 0; strid < NUM_STR; ++strid) {
      if (best_threshold <= 1 && next_threshold <= 1) {
        // if we already have 2 exact match/or dist 1 (one for best, one for second best for mapq), look for exact matches only
        nmismatch = (memcmp(r, ptr_ref + region.rs, rlen) == 0 ? 0 : elen);
      } else {
        nmismatch = cgk2_unmatched(r, ptr_ref + region.rs, region.matched_intervals,
                                   rlen, kmer_step, next_threshold, strid);
      }
      region.embed_dist = region.embed_dist < nmismatch ? region.embed_dist : nmismatch;

      // if embed_dist is 0/1, no need to embed again
      if (region.embed_dist == 0 || region.embed_dist == 1)
        break;
    }

    // set best and next idx so that we dont have to sort regions later
    // pick the minimal pos hit when several hits have same embed dist and the best one is not the hcov(0)
    if (region.embed_dist < best_threshold ||
        (region.embed_dist == best_threshold && best_idx != 0 && region.rs < candidate_regions[best_idx].rs)) {
      if (i != best_idx) {
        next_threshold = best_threshold;
        next_idx = best_idx;
      }
      best_threshold = region.embed_dist;
      best_idx = i;
    } else if (region.embed_dist < next_threshold) {
      if (i != best_idx) {
        next_threshold = region.embed_dist;
        next_idx = i;
      }
    }

  }

  auto end = std::chrono::system_clock::now();
  auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
  embed_time += elapsed.count();
}

void Embedding::embed_unmatch(vector<Region> &candidate_regions,
                              const char *ptr_ref,
                              const char *r,
                              const unsigned rlen,
                              const unsigned kmer_step,
                              bool flag_f1[]) {
  auto start = std::chrono::system_clock::now();

  int elen = rlen * efactor, nmismatch;

  for (unsigned i = 0; i < candidate_regions.size(); ++i) {
    if (!flag_f1[i]) {
      continue;
    }

    Region &region = candidate_regions[i];
    region.embed_dist = elen;

    for (unsigned strid = 0; strid < NUM_STR; ++strid) {
      nmismatch = cgk2_unmatched(r, ptr_ref + region.rs, region.matched_intervals, rlen, kmer_step, elen, strid);
      region.embed_dist = region.embed_dist < nmismatch ? region.embed_dist : nmismatch;

      // if embed_dist is 0/1, no need to embed again
      if (region.embed_dist == 0 || region.embed_dist == 1)
        break;
    }

  }

  auto end = std::chrono::system_clock::now();
  auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
  embed_time += elapsed.count();
}

void Embedding::embed_unmatch_pair(Read &mate1, Read &mate2,
                                   vector<Region> &candidate_regions_f1, vector<Region> &candidate_regions_r2,
                                   const char *ptr_ref, const char *r, const unsigned rlen, const unsigned kmer_step,
                                   bool flag_r2[], unsigned pairdis, int &best_threshold, int &next_threshold,
                                   unsigned &best_f1, unsigned &best_r2) {
  auto start = std::chrono::system_clock::now();

  int elen = rlen * efactor, nmismatch;
  best_threshold = next_threshold = elen;

  for (unsigned i = 0; i < candidate_regions_r2.size(); i++) {
    if (!flag_r2[i]) {
      continue;
    }

    Region &region = candidate_regions_r2[i];
    region.embed_dist = elen;

    Region tmp;
    tmp.rs = region.rs < pairdis ? 0 : region.rs - pairdis;
    auto start = lower_bound(candidate_regions_f1.begin(), candidate_regions_f1.end(), tmp,
                             [](const Region &left, const Region &right) {
                               return left.rs < right.rs;
                             }
    );

    tmp.rs = region.rs + pairdis;
    auto end = upper_bound(candidate_regions_f1.begin(), candidate_regions_f1.end(), tmp,
                           [](const Region &left, const Region &right) {
                             return left.rs < right.rs;
                           }
    );

    assert(start != end);

    for (unsigned strid = 0; strid < NUM_STR; ++strid) {
      if (best_threshold <= 1 && next_threshold <= 1) {
        // if we already have 2 exact match/or dist 1 (one for best, one for second best for mapq), look for exact matches only
        nmismatch = (memcmp(r, ptr_ref + region.rs, rlen) == 0 ? 0 : elen);
      } else {
        nmismatch = cgk2_unmatched(r, ptr_ref + region.rs, region.matched_intervals, rlen, kmer_step,
                                   min(int(region.embed_dist), next_threshold), strid);
      }
      region.embed_dist = region.embed_dist < nmismatch ? region.embed_dist : nmismatch;

      // if embed_dist is 0/1, no need to embed again
      if (region.embed_dist == 0 || region.embed_dist == 1)
        break;
    }

    for (auto itr = start; itr != end; ++itr) {
      int sum_dist = region.embed_dist + itr->embed_dist;

      if (sum_dist <= best_threshold) {
        next_threshold = best_threshold;
        best_threshold = sum_dist;
        best_f1 = itr - candidate_regions_f1.begin();
        best_r2 = i;
        mate1.secBest = mate1.best;
        mate2.secBest = mate2.best;
        mate1.best = mate1.best < itr->embed_dist ? mate1.best : itr->embed_dist;
        mate2.best = mate2.best < region.embed_dist ? mate2.best : region.embed_dist;
      } else if (sum_dist < next_threshold) {
        next_threshold = sum_dist;
      }
    }
  }

  auto end = std::chrono::system_clock::now();
  auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
  embed_time += elapsed.count();
}

Embedding::Embedding() {
#ifdef CGK2_EMBED
  efactor = 2;
#else
  efactor = 3;
#endif

  time_t seed = time(NULL);
  srand(seed);
  cerr << "Embedding using random seed " << seed << endl;
//  srand(1559063236);
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

Embedding::Embedding(const char *fname) {
#ifdef CGK2_EMBED
  efactor = 2;
#else
  efactor = 3;
#endif
  cerr << "Loading embedding from file " << fname << endl;
  ifstream input(fname);
  if (input) {
    char num_str, num_char;
    input.read((char *) &num_str, sizeof(int));
    input.read((char *) &num_char, sizeof(int));
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

Embedding::~Embedding() {
}

