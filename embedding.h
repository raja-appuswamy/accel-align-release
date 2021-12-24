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
  void cgk2_embedQ(const char *oridata, unsigned rlen, int strid, char *embeddedQ);
  int cgk2_embed_nmismatch(const char *oridata, unsigned rlen, int threshold, int strid, char *embeddedQ);
  void embed_unmatch_iter(vector<Region> &candidate_regions, const char *ptr_ref, const char *r, const unsigned rlen,
                          const unsigned kmer_len, int &best_threshold, int &next_threshold,
                          unsigned &best_idx, unsigned &next_idx);
  void embed_unmatch(vector<Region> &candidate_regions, const char *ptr_ref, const char *r, const unsigned rlen,
                     const unsigned kmer_step, bool flag_f1[]);
  void embed_unmatch_pair(Read &mate1, Read &mate2,
                          vector<Region> &candidate_regions_f1, vector<Region> &candidate_regions_r2,
                          const char *ptr_ref, const char *r, const unsigned rlen, const unsigned kmer_step,
                          bool flag_r2[], unsigned pairdis, int &best_threshold, int &next_threshold,
                          unsigned &best_f1, unsigned &best_r2);

  //unsigned char **hash_eb;
  std::bitset<TOTAL_RBITS> hash_eb;
  float embed_time;
  int efactor;

 private:

  int cgk2_embed(const char **oridata, unsigned rlen, int threshold, int id,
                 int strid, char *embeddedQ);
  int cgk_embed(const char **oridata, unsigned rlen, int threshold, int id,
                int strid, char *embeddedQ);
  int cgk2_unmatched(const char *ref,
                     const char *r,
                     const vector<uint32_t> &mch,
                     const unsigned rlen,
                     const unsigned kmer_len,
                     const int threshold,
                     const int strid);

};
