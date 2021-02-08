#pragma once
#include "header.h"

using namespace std;

class Embedding {
public:
  ~Embedding();
  Embedding();
  Embedding(const char *fname);
  int embedstr(const char **oridata, unsigned rlen, int threshold, int id,
               int strid, char *embeddedQ);
  void embeddata_iterative_update(vector<Region> &candidate_regions,
                                  const char **input, unsigned ninput, unsigned rlen,
                                  int &best_threshold, int &next_threshold,
                                  bool max_rnd, int &best_idx, int &next_idx);
  void embeddata_pair(vector<Region> &candidate_regions, char embeddedQ[],
                      const char **candidate_refs, unsigned ncandidates, bool flag[],
                      unsigned rlen, int threshold);
  void embed_two_pairs(vector<Region> &candidate_regions_f1, vector<Region> &candidate_regions_r2,
                       unsigned rlen, int *thresholds, char embeddedQ_f1[], char embeddedQ_r2[],
                       const char **candidate_refs_f1, const char **candidate_refs_r2,
                       unsigned nregions_f1, unsigned nregions_r2, bool flag_f1[], bool flag_r2[],
                       int best_f1, int next_f1, int best_r2, int next_r2);
  //unsigned char **hash_eb;
  std::bitset<TOTAL_RBITS> hash_eb;
  float embed_time;
  int efactor;

private:

  int cgk2_embed(const char **oridata, unsigned rlen, int threshold, int id,
                 int strid, char *embeddedQ);
  int cgk_embed(const char **oridata, unsigned rlen, int threshold, int id,
                int strid, char *embeddedQ);

};
