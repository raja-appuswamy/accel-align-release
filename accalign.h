#pragma once
class AccAlign {
private:
  std::string &ref;
  std::vector<std::string> &name;
  std::vector<uint32_t> &offset;
  Embedding *embedding;

  float input_io_time, parse_time;
  float seeding_time, hit_count_time, vpair_build_time;
  float embedding_time, swap_time;
  float sw_time, sam_time, sam_pre_time, sam_out_time;
  long vpair_sort_count;

  // output fields
  std::string sam_name;
  std::ofstream sam_stream;

  void cpu_root_fn(tbb::concurrent_bounded_queue<ReadCnt> *inputQ,
                   tbb::concurrent_bounded_queue<ReadCnt> *outputQ);
  void output_root_fn(tbb::concurrent_bounded_queue<ReadCnt> *outputQ,
                      tbb::concurrent_bounded_queue<ReadPair> *dataQ);
  void align_wrapper(int tid, int soff, int eoff,
                     Read *ptlread, Read *ptlread2,
                     tbb::concurrent_bounded_queue<ReadPair> *dataQ);
  void embed_wrapper(Read &R, bool ispe, vector<Region> &fregion,
                     vector<Region> &rregion, int &fbest, int &fnext, int &rbest,
                     int &rnext);
  void embed_wrapper_pair(Read &R1, Read &R2,
                          vector<Region> &candidate_regions_f1, vector<Region> &candidate_regions_r1,
                          vector<Region> &candidate_regions2_f2, vector<Region> &rcandidate_regions2,
                          bool flag_f1[], bool flag_r1[], bool flag_f2[], bool flag_r2[],
                          int &best_f1, int &next_f1, int &best_r1, int &next_r1,
                          int &best_f2, int &next_f2, int &best_r2, int &next_r2);
  void pghole_wrapper(Read &R, vector<Region> &fcandidate_regions,
                      vector<Region> &rcandidate_regions, int &fbest, int &rbest, int ori_slide);
  void pigeonhole_query(char *Q, size_t rlen, vector<Region> &candidate_regions,
                        char S, int err_threshold, int &best, int ori_slide);
  void pairdis_filter(vector<Region> &in_regions1, vector<Region> &in_regions2,
                      bool flag1[], bool flag2[],
                      int &best1, int &next1, int &best2, int &next2);
  void lsh_filter(char *Q, size_t rlen,
                  vector<Region> &candidate_regions,
                  int &best_threshold, int &next_threshold,
                  int &best_idx, int &next_idx);
  void set_mapq(Read &R, vector<Region> &fcandidate_regions,
                vector<Region> &rcandidate_regions, int fbest, int fnext,
                int rbest, int rnext);
  void mark_for_extension(Read &read, char S, Region &cregion);
  void save_region(Read &R, size_t rlen, Region &region,
                   Alignment &a);
  void score_region(Read &r, char *strand, Region &region,
                    Alignment &a);
  void sam_header(void);
  int get_mapq(int best, int secbest, bool hasSecbest, int rlen);

public:
  uint32_t *keyv, *posv;

  void open_output(std::string &out_file);
  void close_output();
  bool fastq(const char *F1, const char *F2, bool enable_gpu);
  void print_stats();
  void map_read(Read &R);
  void align_read(Read &R);
  void out_sam(string *sam);
  void snprintf_sam(Read &R, string *s);
  void snprintf_pair_sam(Read &R, string *s, Read &R2, string *s2);
  void map_paired_read(Read &mate1, Read &mate2);
  void wfa_align_read(Read &R);
  AccAlign(Reference &r);
  ~AccAlign();
};


