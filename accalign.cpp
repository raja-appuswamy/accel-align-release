#include "header.h"
#include "accalign.h"
#include "ksw2.h"

using namespace tbb::flow;
using namespace std;

unsigned kmer_len = 32;
int kmer_step = 1;
uint64_t mask;
unsigned pairdis = 1000;
string g_out, g_batch_file, g_embed_file;
char rcsymbol[6] = "TGCAN";
uint8_t code[256];
bool enable_extension = true, enable_wfa_extension = false, extend_all = false;

int g_ncpus = 1;
float delTime = 0, mapqTime = 0, keyvTime = 0, posvTime = 0, sortTime = 0;
int8_t mat[25];

void make_code(void) {
  for (size_t i = 0; i < 256; i++)
    code[i] = 4;

  code['A'] = code['a'] = 0;
  code['C'] = code['c'] = 1;
  code['G'] = code['g'] = 2;
  code['T'] = code['t'] = 3;

  // we set N's also to 0. Shouldn't matter as index construction logic
  // doesn't consider kmers with Ns anyway
  code['N'] = code['n'] = 0;
}

static void parse(char seq[], char fwd[], char rev[], char rev_str[]) {
  unsigned len = strlen(seq);

  for (size_t i = 0; i < len; i++) {
    uint8_t c = *(code + seq[i]);
    fwd[i] = c;
    rev[len - 1 - i] = c == 4 ? c : 3 - c;
    rev_str[len - 1 - i] = rcsymbol[c];
  }
  *(fwd + len) = '\0';
  *(rev + len) = '\0';
  *(rev_str + len) = '\0';
}

gzFile &operator>>(gzFile &in, Read &r) {
  char temp[MAX_LEN];

  if (gzgets(in, r.name, MAX_LEN) == NULL)
    return in;
  if (gzgets(in, r.seq, MAX_LEN) == NULL)
    return in;
  if (gzgets(in, temp, MAX_LEN) == NULL)
    return in;
  if (gzgets(in, r.qua, MAX_LEN) == NULL)
    return in;

  unsigned i = 0;
  while (i < strlen(r.name)) {
    if (isspace(r.name[i])) { // isspace(): \t, \n, \v, \f, \r
      memset(r.name + i, '\0', strlen(r.name) - i);
      break;
    }
    i++;
  }

  r.qua[strlen(r.qua) - 1] = '\0';
  r.seq[strlen(r.seq) - 1] = '\0';

  r.tid = r.pos = 0;
  r.as = numeric_limits<int32_t>::min();
  r.strand = '*';

  return in;
}

void print_usage() {
  cerr << "accalign [options] <ref.fa> [read1.fastq] [read2.fastq]\n";
  cerr << "\t Maximum read length supported is 512\n";
  cerr << "options:\n";
  cerr << "\t-t INT Number of cpu threads to use [all]\n";
  cerr << "\t-l INT Length of seed [32]\n";
  cerr << "\t-o Name of the output file \n";
  cerr << "\t-x Alignment-free mode\n";
  cerr << "\t-w Use WFA for extension. KSW used by default. \n";
  cerr << "\t-p Maximum distance allowed between the paired-end reads [1000]\n";
  cerr << "\t-d Disable embedding, extend all candidates from seeding (this mode is super slow, only for benchmark).\n";
}

void AccAlign::print_stats() {
//#if DBGPRINT
  cerr << "Breakdown:\n" <<
       "Input IO time:\t" << input_io_time / 1000000 << "\n" <<
       "Parse time:\t" << parse_time / 1000000 << "\n" <<
       "Seeding time: \t" << seeding_time / 1000000 << "\n" <<
       "\t lookup keyv time:\t" << keyvTime / 1000000 << "\n" <<
       "\t lookup posv time:\t" << posvTime / 1000000 << "\n" <<
       "\t sort hits time:\t" << sortTime / 1000000 << "\n" <<
       "\t Hit count time:\t" << hit_count_time / 1000000 << "\n" <<
       "\t Swap high cov time:\t" << swap_time / 1000000 << "\n" <<
       "\t Vpair build (only for pe):\t" << vpair_build_time / 1000000 << "\n" <<
       "Embedding time(total/actual):\t" << embedding->embed_time / 1000000 << "\n" <<
       "Extending time (+ build output string if ENABLE_GPU):\t" << sw_time / 1000000 << "\n" <<
       "Mark best region time:\t" << mapqTime / 1000000 << "\n" <<
       "SAM output time :\t" << sam_time / 1000000 << "\n" <<
       std::endl << endl;

  cerr << "Total pairs sorted: " << vpair_sort_count << endl;
//#endif
}

bool AccAlign::fastq(const char *F1, const char *F2, bool enable_gpu) {

  bool is_paired = false;

  gzFile in1 = gzopen(F1, "rt");
  if (in1 == Z_NULL)
    return false;

  gzFile in2 = Z_NULL;
  if (strlen(F2) > 0) {
    is_paired = true;
    in2 = gzopen(F2, "rt");
    if (in2 == Z_NULL)
      return false;
  }

  cerr << "Reading fastq file " << F1 << ", " << F2 << "\n";

  // start CPU and GPU master threads, they consume reads from inputQ
  // dataQ is to re-use
  tbb::concurrent_bounded_queue<ReadCnt> inputQ;
  tbb::concurrent_bounded_queue<ReadCnt> outputQ;
  tbb::concurrent_bounded_queue<ReadPair> dataQ;

  thread cpu_thread = thread(&AccAlign::cpu_root_fn, this, &inputQ, &outputQ);
  thread out_thread = thread(&AccAlign::output_root_fn, this, &outputQ, &dataQ);

  auto start = std::chrono::system_clock::now();

  int total_nreads = 0, nreads_per_vec = 0, vec_index = 0, vec_size = 50;

  int batch_size = BATCH_SIZE;
  if (is_paired)
    batch_size /= 2;
  Read *reads[vec_size];
  reads[vec_index] = new Read[batch_size];
  Read *reads2[vec_size];
  if (is_paired) {
    reads2[vec_index] = new Read[batch_size];
  }

  bool neof1 = (!gzeof(in1) && gzgetc(in1) != EOF);
  bool neof2 = (!is_paired || (!gzeof(in2) && gzgetc(in2) != EOF));
  while (vec_index < vec_size && neof1 && neof2) {
    Read &r = *(reads[vec_index] + nreads_per_vec);
    in1 >> r;

    if (!strlen(r.seq)) {
      break;
    }

    if (is_paired) {
      Read &r2 = *(reads2[vec_index] + nreads_per_vec);
      in2 >> r2;

      if (!strlen(r2.seq)) {
        break;
      }
    }
    neof1 = (!gzeof(in1) && gzgetc(in1) != EOF);
    neof2 = (!is_paired || (!gzeof(in2) && gzgetc(in2) != EOF));

    ++nreads_per_vec;

    if (nreads_per_vec == batch_size) {
      if (is_paired)
        inputQ.push(make_tuple(reads[vec_index], reads2[vec_index], batch_size));
      else
        inputQ.push(make_tuple(reads[vec_index], (Read *) NULL, batch_size));

      vec_index++;

      if (vec_index < vec_size) {
        reads[vec_index] = new Read[batch_size];
        if (is_paired)
          reads2[vec_index] = new Read[batch_size];
      }

      total_nreads += nreads_per_vec;
      nreads_per_vec = 0;
    }
  }

  ReadPair cur_vec = make_tuple((Read *) NULL, (Read *) NULL);

  // the nb of reads is less than vec_size *BATCH_SIZE, and there are some reads not pushed to inputQ
  if (nreads_per_vec && vec_index < vec_size) {
    // the remaining reads
    if (is_paired)
      inputQ.push(make_tuple(reads[vec_index], reads2[vec_index], nreads_per_vec));
    else
      inputQ.push(make_tuple(reads[vec_index], (Read *) NULL, nreads_per_vec));

    total_nreads += nreads_per_vec;
  } else {
    // still have reads not loaded
    dataQ.pop(cur_vec);

    while (neof1 && neof2) {
      Read &r = *(std::get<0>(cur_vec) + nreads_per_vec);
      in1 >> r;

      if (!strlen(r.seq)) {
        break;
      }

      if (is_paired) {
        Read &r2 = *(std::get<1>(cur_vec) + nreads_per_vec);
        in2 >> r2;

        if (!strlen(r2.seq)) {
          break;
        }
      }
      neof1 = (!gzeof(in1) && gzgetc(in1) != EOF);
      neof2 = (!is_paired || (!gzeof(in2) && gzgetc(in2) != EOF));

      ++nreads_per_vec;

      if (nreads_per_vec == batch_size) {
        inputQ.push(make_tuple(std::get<0>(cur_vec), std::get<1>(cur_vec), batch_size));
        dataQ.pop(cur_vec);
        total_nreads += nreads_per_vec;
        nreads_per_vec = 0;
      }
    }

    // the remaining reads
    if (nreads_per_vec) {
      total_nreads += nreads_per_vec;
      inputQ.push(make_tuple(std::get<0>(cur_vec), std::get<1>(cur_vec), nreads_per_vec));
    }
  }

  gzclose(in1);
  if (is_paired) {
    gzclose(in2);
  }

  auto end = std::chrono::system_clock::now();
  auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
  input_io_time += elapsed.count();

  cerr << "done reading " << total_nreads << " reads from fastq file " << F1 << ", " << F2 << " in " <<
       input_io_time / 1000000.0 << " secs\n";

  ReadCnt sentinel = make_tuple((Read *) NULL, (Read *) NULL, 0);
  inputQ.push(sentinel);

  int size = vec_index < vec_size ? vec_index : vec_size;
  start = std::chrono::system_clock::now();
  if (total_nreads % batch_size == 0) {
    //because the last popped cur_vec has not been pushed back
    size -= 1;
    delete[] std::get<0>(cur_vec);

    if (is_paired)
      delete[] std::get<1>(cur_vec);
  }
  for (int i = 0; i < size; i++) {
    dataQ.pop(cur_vec);
    delete[] std::get<0>(cur_vec);

    if (is_paired)
      delete[] std::get<1>(cur_vec);
  }
  end = std::chrono::system_clock::now();
  elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
  delTime += elapsed.count();

  cpu_thread.join();

  outputQ.push(sentinel);
  out_thread.join();

  cerr << "Processed " << total_nreads << " in total \n";
  return true;
}

void AccAlign::output_root_fn(tbb::concurrent_bounded_queue<ReadCnt> *outputQ,
                              tbb::concurrent_bounded_queue<ReadPair> *dataQ) {
  cerr << "Extension and output function starting.." << endl;

  unsigned nreads = 0;
  tbb::concurrent_bounded_queue<ReadCnt> *targetQ = outputQ;
  do {
    ReadCnt gpu_reads;
    targetQ->pop(gpu_reads);
    nreads = std::get<2>(gpu_reads);
    if (!nreads) {
      targetQ->push(gpu_reads);   //put sentinel back
      break;
    }
    align_wrapper(0, 0, nreads, std::get<0>(gpu_reads), std::get<1>(gpu_reads), dataQ);
  } while (1);

  cerr << "Extension and output function quitting...\n";
}

class Parallel_mapper {
  Read *all_reads1;
  Read *all_reads2;
  AccAlign *acc_obj;

 public:
  Parallel_mapper(Read *_all_reads1, Read *_all_reads2, AccAlign *_acc_obj) :
      all_reads1(_all_reads1), all_reads2(_all_reads2), acc_obj(_acc_obj) {}

  void operator()(const tbb::blocked_range<size_t> &r) const {
    if (!all_reads2) {
      for (size_t i = r.begin(); i != r.end(); ++i) {
        acc_obj->map_read(*(all_reads1 + i));
      }
    } else {
      for (size_t i = r.begin(); i != r.end(); ++i) {
        acc_obj->map_paired_read(*(all_reads1 + i), *(all_reads2 + i));
      }
    }
  }
};

void AccAlign::cpu_root_fn(tbb::concurrent_bounded_queue<ReadCnt> *inputQ,
                           tbb::concurrent_bounded_queue<ReadCnt> *outputQ) {
  cerr << "CPU Root function starting.." << endl;

  tbb::concurrent_bounded_queue<ReadCnt> *targetQ = inputQ;
  int nreads = 0, total = 0;
  do {
    ReadCnt cpu_readcnt;
    targetQ->pop(cpu_readcnt);
    nreads = std::get<2>(cpu_readcnt);
    total += nreads;
    if (nreads == 0) {
      inputQ->push(cpu_readcnt);    // push sentinel back
      break;
    }

    tbb::task_scheduler_init init(g_ncpus);
    tbb::parallel_for(tbb::blocked_range<size_t>(0, nreads),
                      Parallel_mapper(std::get<0>(cpu_readcnt), std::get<1>(cpu_readcnt), this)
    );

    outputQ->push(cpu_readcnt);
  } while (1);

  cerr << "Processed " << total << " reads in cpu \n";
  cerr << "CPU Root function quitting.." << endl;
}

void AccAlign::mark_for_extension(Read &read, char S, Region &cregion) {
  read.strand = S;
  int rlen = strlen(read.seq);

  cregion.re = cregion.rs + rlen < ref.size() ? cregion.rs + rlen :
               ref.size();

  char *strand = S == '+' ? read.fwd : read.rev;

  if (cregion.embed_dist && !enable_extension)
    rectify_start_pos(strand, cregion, rlen);

  read.best_region = cregion;
}

void AccAlign::pigeonhole_query_topcov(char *Q,
                                       size_t rlen,
                                       vector<Region> &candidate_regions,
                                       char S,
                                       int err_threshold,
                                       unsigned kmer_step,
                                       unsigned max_occ,
                                       unsigned &best,
                                       unsigned ori_slide) {
  int max_cov = 0;
  unsigned nkmers = (rlen - ori_slide - kmer_len) / kmer_step + 1;
  size_t ntotal_hits = 0;
  size_t b[nkmers], e[nkmers];
  unsigned kmer_idx = 0;
  unsigned ori_slide_bk = ori_slide;
  unsigned nseed_freq = 0;
  bool high_freq = false;

  // Take non-overlapping seeds and find all hits
  auto start = std::chrono::system_clock::now();
  for (size_t i = ori_slide; i + kmer_len <= rlen; i += kmer_step) {
    uint64_t k = 0;
    for (size_t j = i; j < i + kmer_len; j++)
      k = (k << 2) + *(Q + j);
    size_t hash = (k & mask) % MOD;
    b[kmer_idx] = keyv[hash];
    e[kmer_idx] = keyv[hash + 1];
    if (e[kmer_idx] - b[kmer_idx] >= max_occ)
      nseed_freq++;
//    if (e[kmer_idx] - b[kmer_idx] < max_occ) {
//      ntotal_hits += (e[kmer_idx] - b[kmer_idx]);
//    }
    kmer_idx++;
  }
  assert(kmer_idx == nkmers);
  auto end = std::chrono::system_clock::now();
  auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
  keyvTime += elapsed.count();

  if (nseed_freq > nkmers / 2)
    high_freq = true;

  for (size_t i = 0; i < nkmers; i++) {
    if ((!high_freq && e[i] - b[i] < max_occ) || high_freq)
      ntotal_hits += (e[i] - b[i]);
  }

  // if we have no hits, we are done
  if (!ntotal_hits)
    return;

  uint32_t top_pos[nkmers];
  int rel_off[nkmers];
  uint32_t MAX_POS = numeric_limits<uint32_t>::max();

  start = std::chrono::system_clock::now();
  // initialize top values with first values for each kmer.
  for (unsigned i = 0; i < nkmers; i++) {
    if (b[i] < e[i] && ((!high_freq && e[i] - b[i] < max_occ) || high_freq)) {
//    if (b[i] < e[i] && e[i] - b[i] < max_occ) {
      top_pos[i] = posv[b[i]];
      rel_off[i] = i * kmer_step;
      uint32_t shift_pos = rel_off[i] + ori_slide_bk;
      //TODO: for each chrome, happen to < the start pos
      if (top_pos[i] < shift_pos)
        top_pos[i] = 0; // there is insertion before this kmer
      else
        top_pos[i] -= shift_pos;
    } else {
      top_pos[i] = MAX_POS;
    }
  }
  end = std::chrono::system_clock::now();
  elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
  posvTime += elapsed.count();

  size_t nprocessed = 0;
  uint32_t last_pos = MAX_POS, last_qs = ori_slide_bk; //last query start pos
  int last_cov = 0;

  start = std::chrono::system_clock::now();

  vector<Region> unique_regions;
  unique_regions.reserve(ntotal_hits);
  size_t idx = 0;
  Region r;
  r.matched_intervals.reserve(nkmers);
//  Region unique_regions[ntotal_hits];
//  size_t idx = 0;
//  Region *r = unique_regions + idx;

  while (nprocessed < ntotal_hits) {
    //find min
    uint32_t *min_item = min_element(top_pos, top_pos + nkmers);
    uint32_t min_pos = *min_item;
    int min_kmer = min_item - top_pos;

    if ((!high_freq && e[min_kmer] - b[min_kmer] < max_occ) || high_freq) {
//    if (e[min_kmer] - b[min_kmer] < max_occ) {
      // kick off prefetch for next round
      __builtin_prefetch(posv + b[min_kmer] + 1);

      // if previous min element was same as current one, increment coverage.
      // otherwise, check if last min element's coverage was high enough to make it a candidate region

      if (min_pos == last_pos) {
        r.matched_intervals.push_back(last_qs);
        last_cov++;
      } else {
        if (nprocessed != 0) {
          r.cov = last_cov;
          r.rs = last_pos;
          r.matched_intervals.push_back(last_qs);
          r.qs = r.matched_intervals[0]; //the first match seed, so left extension could be accurate
          r.qe = r.qs + kmer_len;

          if (last_cov > max_cov)
            max_cov = last_cov;

          assert(r.rs != MAX_POS && r.rs < MAX_POS);
          unique_regions.push_back(move(r));
          ++idx;
//          r =  unique_regions + idx;
        }

        last_cov = 1;
      }
      last_qs = min_kmer * kmer_step + ori_slide_bk;
      last_pos = min_pos;
    }

    // add next element
    b[min_kmer]++;
    uint32_t next_pos = b[min_kmer] < e[min_kmer] ? posv[b[min_kmer]] : MAX_POS;
    if (next_pos != MAX_POS) {
      uint32_t shift_pos = rel_off[min_kmer] + ori_slide_bk;
      //TODO: for each chrome, happen to < the start pos
      if (next_pos < shift_pos)
        *min_item = 0; // there is insertion before this kmer
      else
        *min_item = next_pos - shift_pos;
    } else
      *min_item = MAX_POS;

    ++nprocessed;
  }

  // we will have the last few positions not processed. check here.
  if (last_pos != MAX_POS) {
    r.cov = last_cov;
    r.rs = last_pos;
    r.matched_intervals.push_back(last_qs);
    r.qs = r.matched_intervals[0]; //the first match seed, so left extension could be accurate
    r.qe = r.qs + kmer_len;

    if (last_cov > max_cov)
      max_cov = last_cov;

    assert(r.rs != MAX_POS && r.rs < MAX_POS);
    unique_regions.push_back(move(r));
    ++idx;
  }

  err_threshold = max(err_threshold, max_cov - 1);
  assert(idx <= ntotal_hits);
  for (size_t i = 0; i < idx; i++) {
    Region &r = unique_regions[i];
    if (r.cov >= err_threshold) {
      if (r.cov == max_cov)
        best = candidate_regions.size();

      candidate_regions.push_back(move(r));
    }
  }

  end = std::chrono::system_clock::now();
  elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
  hit_count_time += elapsed.count();
}

void AccAlign::pigeonhole_query_sort(char *Q,
                                     size_t rlen,
                                     vector<Region> &candidate_regions,
                                     char S,
                                     unsigned err_threshold,
                                     unsigned kmer_step,
                                     unsigned max_occ,
                                     unsigned &best,
                                     unsigned ori_slide) {
  unsigned max_cov = 0;
  unsigned nkmers = (rlen - ori_slide - kmer_len) / kmer_step + 1;
  size_t ntotal_hits = 0;
  size_t b[nkmers], e[nkmers];
  unsigned kmer_idx = 0;
  unsigned nseed_freq = 0;
  bool high_freq = false;

  // Take non-overlapping seeds and find all hits
  auto start = std::chrono::system_clock::now();
  for (size_t i = ori_slide; i + kmer_len <= rlen; i += kmer_step) {
    uint64_t k = 0;
    for (size_t j = i; j < i + kmer_len; j++)
      k = (k << 2) + *(Q + j);
    size_t hash = (k & mask) % MOD;
    b[kmer_idx] = keyv[hash];
    e[kmer_idx] = keyv[hash + 1];
    if (e[kmer_idx] - b[kmer_idx] >= max_occ)
      nseed_freq++;
//    if (e[kmer_idx] - b[kmer_idx] < max_occ) {
//      ntotal_hits += (e[kmer_idx] - b[kmer_idx]);
//    }
    kmer_idx++;
  }
  assert(kmer_idx == nkmers);
  auto end = std::chrono::system_clock::now();
  auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
  keyvTime += elapsed.count();

  if (nseed_freq > nkmers / 2)
    high_freq = true;

  for (size_t i = 0; i < nkmers; i++) {
    if ((!high_freq && e[i] - b[i] < max_occ) || high_freq)
      ntotal_hits += (e[i] - b[i]);
  }

  // if we have no hits, we are done
  if (!ntotal_hits)
    return;

  start = std::chrono::system_clock::now();
  // initialize top values with first values for each kmer.
  uint32_t MAX_POS = numeric_limits<uint32_t>::max();
  vector<Region> regions;
  regions.reserve(ntotal_hits);
  for (unsigned i = 0; i < nkmers; i++) {
    if (b[i] < e[i] && ((!high_freq && e[i] - b[i] < max_occ) || high_freq)) {
//    if (b[i] < e[i] && e[i] - b[i] < max_occ) {
      for (uint32_t j = b[i]; j < e[i]; j++) {
        Region r;
        r.rs = posv[j];
        r.qs = i * kmer_step + ori_slide;
        r.rs -= min(r.rs, r.qs);
        regions.push_back(r);
        // rs can't be samller than 0, if insertion before this kmer, set rs to 0 instead of -1
      }
    }
  }
  assert(regions.size() == ntotal_hits);
  end = std::chrono::system_clock::now();
  elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
  posvTime += elapsed.count();

  start = std::chrono::system_clock::now();

  sort(regions.begin(), regions.end(), Region());

  size_t nprocessed = 0, last_cov = 0;
  uint32_t last_pos = MAX_POS;

  while (nprocessed < ntotal_hits) {

    if (regions[nprocessed].rs == last_pos) {
      last_cov++;
    } else {
      if (last_cov >= err_threshold) {
        Region r;
        r.cov = last_cov;
        r.rs = last_pos;
        for (unsigned i = nprocessed - last_cov; i < nprocessed; i++)
          r.matched_intervals.push_back(regions[i].qs);
        r.qs = r.matched_intervals[0]; //the first match seed, so left extension could be accurate
        r.qe = r.qs + kmer_len;

        assert(r.rs < MAX_POS);

        if (last_cov >= max_cov) {
          max_cov = last_cov;
          best = candidate_regions.size();
        }
        candidate_regions.push_back(r);
      }
      last_cov = 1;
    }
    last_pos = regions[nprocessed].rs;

    ++nprocessed;
  }

  // we will have the last few positions not processed. check here.
  if (last_cov >= err_threshold && last_pos != MAX_POS) {
    Region r;
    r.cov = last_cov;
    r.rs = last_pos;
    for (unsigned i = nprocessed - last_cov; i < nprocessed; i++)
      r.matched_intervals.push_back(regions[i].qs);
    r.qs = r.matched_intervals[0]; //the first match seed, so left extension could be accurate
    r.qe = r.qs + kmer_len;
    assert(r.rs < MAX_POS);

    if (last_cov >= max_cov) {
      max_cov = last_cov;
      best = candidate_regions.size();
    }
    candidate_regions.push_back(r);
  }

  end = std::chrono::system_clock::now();
  elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
  posvTime += elapsed.count();
}

void AccAlign::pghole_wrapper(Read &R,
                              vector<Region> &fcandidate_regions,
                              vector<Region> &rcandidate_regions,
                              unsigned &fbest,
                              unsigned &rbest) {
  size_t rlen = strlen(R.seq);
  unsigned kmer_step = kmer_len, nfregions = 0, nrregions = 0;
  bool high_freq = false;
  unsigned ori_slide = 0;

  unsigned slide = kmer_len < rlen - kmer_len ? kmer_len : rlen - kmer_len;

  // MAX_OCC, cov >= 2
//  while (kmer_step > 0 && !nfregions && !nrregions) {
  while (ori_slide < slide && !nfregions && !nrregions) {

    unsigned nkmers = (rlen - ori_slide - kmer_len) / kmer_step + 1;

    if (nkmers < 4) {
      //nkmer 3, 2, 1, top 2 cov of cov >=2, is 3, 2, is as same as cov>=2
      // as cov2 is faster than top2, use cov2
      pigeonhole_query(R.fwd, rlen, fcandidate_regions, '+', fbest, ori_slide, 2, kmer_step, MAX_OCC, high_freq);
      pigeonhole_query(R.rev, rlen, rcandidate_regions, '-', rbest, ori_slide, 2, kmer_step, MAX_OCC, high_freq);
    } else {
      pigeonhole_query_topcov(R.fwd, rlen, fcandidate_regions, '+', 2, kmer_step, MAX_OCC, fbest, ori_slide);
      pigeonhole_query_topcov(R.rev, rlen, rcandidate_regions, '-', 2, kmer_step, MAX_OCC, rbest, ori_slide);
    }
    nfregions = fcandidate_regions.size();
    nrregions = rcandidate_regions.size();

    if (!nfregions && !nrregions) {
      pigeonhole_query(R.fwd, rlen, fcandidate_regions, '+', fbest, ori_slide, 1, kmer_step, MAX_OCC, high_freq);
      pigeonhole_query(R.rev, rlen, rcandidate_regions, '-', rbest, ori_slide, 1, kmer_step, MAX_OCC, high_freq);
      nfregions = fcandidate_regions.size();
      nrregions = rcandidate_regions.size();
    }

    R.kmer_step = kmer_step;
    ori_slide++;
//    kmer_step = kmer_step / 2;
  }

}

void AccAlign::pigeonhole_query(char *Q,
                                size_t rlen,
                                vector<Region> &candidate_regions,
                                char S,
                                unsigned &best,
                                unsigned ori_slide,
                                int err_threshold,
                                unsigned kmer_step,
                                unsigned max_occ,
                                bool &high_freq) {
  int max_cov = 0;
  unsigned nkmers = (rlen - ori_slide - kmer_len) / kmer_step + 1;
  size_t ntotal_hits = 0;
  size_t b[nkmers], e[nkmers];
  unsigned kmer_idx = 0;
  unsigned nseed_freq = 0;

  // Take non-overlapping seeds and find all hits
  auto start = std::chrono::system_clock::now();
  for (size_t i = ori_slide; i + kmer_len <= rlen; i += kmer_step) {
    uint64_t k = 0;
    for (size_t j = i; j < i + kmer_len; j++)
      k = (k << 2) + *(Q + j);
    size_t hash = (k & mask) % MOD;
    b[kmer_idx] = keyv[hash];
    e[kmer_idx] = keyv[hash + 1];
    if (e[kmer_idx] - b[kmer_idx] >= max_occ)
      nseed_freq++;
    kmer_idx++;
  }
  assert(kmer_idx == nkmers);
  auto end = std::chrono::system_clock::now();
  auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
  keyvTime += elapsed.count();

  if (nseed_freq > nkmers / 2)
    high_freq = true;

  for (size_t i = 0; i < nkmers; i++) {
    if ((!high_freq && e[i] - b[i] < max_occ) || high_freq)
      ntotal_hits += (e[i] - b[i]);
  }

  // if we have no hits, we are done
  if (!ntotal_hits)
    return;

  uint32_t top_pos[nkmers], MAX_POS = numeric_limits<uint32_t>::max();
  int rel_off[nkmers];

  start = std::chrono::system_clock::now();
  // initialize top values with first values for each kmer.
  for (unsigned i = 0; i < nkmers; i++) {
    if (b[i] < e[i] && ((!high_freq && e[i] - b[i] < max_occ) || high_freq)) {
      top_pos[i] = posv[b[i]];
      rel_off[i] = i * kmer_step;
      uint32_t shift_pos = rel_off[i] + ori_slide;
      top_pos[i] -= min(top_pos[i], shift_pos); //pos can't <0, e.g. insertion before this kmer, set 0 instead of -1
    } else {
      top_pos[i] = MAX_POS;
    }
  }
  end = std::chrono::system_clock::now();
  elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
  posvTime += elapsed.count();

  size_t nprocessed = 0;
  uint32_t last_pos = MAX_POS, last_qs = ori_slide; //last query start pos
  int last_cov = 0;

  start = std::chrono::system_clock::now();

  Region r;
  r.matched_intervals.reserve(nkmers);
  while (nprocessed < ntotal_hits) {
    //find min
    uint32_t *min_item = min_element(top_pos, top_pos + nkmers);
    uint32_t min_pos = *min_item;
    int min_kmer = min_item - top_pos;

    if ((!high_freq && e[min_kmer] - b[min_kmer] < max_occ) || high_freq) {
      // kick off prefetch for next round
      __builtin_prefetch(posv + b[min_kmer] + 1);

      // if previous min element was same as current one, increment coverage.
      // otherwise, check if last min element's coverage was high enough to make it a candidate region
      if (min_pos == last_pos) {
        r.matched_intervals.push_back(last_qs);
        last_cov++;
      } else {
        if (last_cov >= err_threshold) {
          r.cov = last_cov;
          r.rs = last_pos;
          r.matched_intervals.push_back(last_qs);
          r.qs = r.matched_intervals[0]; //let it be the first match seed, so the left extension could be accurate
          r.qe = r.qs + kmer_len;

          if (last_cov >= max_cov) {
            max_cov = last_cov;
            best = candidate_regions.size();
          }

          assert(r.rs < MAX_POS);
          candidate_regions.push_back(move(r));
        }

        last_cov = 1;
      }
      last_qs = min_kmer * kmer_step + ori_slide;
      last_pos = min_pos;
    }

    // add next element
    b[min_kmer]++;
    uint32_t next_pos = b[min_kmer] < e[min_kmer] ? posv[b[min_kmer]] : MAX_POS;
    if (next_pos != MAX_POS) {
      uint32_t shift_pos = rel_off[min_kmer] + ori_slide;
      *min_item = next_pos - min(next_pos, shift_pos);
      //pos can't <0, e.g. insertion before this kmer, set 0 instead of -1
    } else
      *min_item = MAX_POS;

    ++nprocessed;
  }

  // we will have the last few positions not processed. check here.
  if (last_pos != MAX_POS) {
    if (last_cov >= err_threshold) {
      r.cov = last_cov;
      r.rs = last_pos;
      r.matched_intervals.push_back(last_qs);
      r.qs = r.matched_intervals[0]; //let it be the first match seed, so the left extension could be accurate
      r.qe = r.qs + kmer_len;

      if (last_cov >= max_cov) {
        max_cov = last_cov;
        best = candidate_regions.size();
      }

      assert(r.rs < MAX_POS);
      candidate_regions.push_back(move(r));
    }
  }

  end = std::chrono::system_clock::now();
  elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
  hit_count_time += elapsed.count();
}

void AccAlign::pghole_wrapper_mates(Read &R,
                                    vector<Region> &fcandidate_regions,
                                    vector<Region> &rcandidate_regions,
                                    unsigned &fbest,
                                    unsigned &rbest,
                                    unsigned ori_slide,
                                    unsigned kmer_step, unsigned max_occ, bool &high_freq) {
  unsigned rlen = strlen(R.seq);

  // MAX_OCC, cov >= 2
  pigeonhole_query(R.fwd, rlen, fcandidate_regions, '+', fbest, ori_slide, 2, kmer_step, max_occ, high_freq);
  pigeonhole_query(R.rev, rlen, rcandidate_regions, '-', rbest, ori_slide, 2, kmer_step, max_occ, high_freq);
  R.kmer_step = kmer_step;
  unsigned nfregions = fcandidate_regions.size();
  unsigned nrregions = rcandidate_regions.size();

  if (!nfregions && !nrregions) {
    pigeonhole_query(R.fwd, rlen, fcandidate_regions, '+', fbest, ori_slide, 1, kmer_step, max_occ, high_freq);
    pigeonhole_query(R.rev, rlen, rcandidate_regions, '-', rbest, ori_slide, 1, kmer_step, max_occ, high_freq);
  }
}

//check by f1r1
void AccAlign::pghole_wrapper_pair(Read &mate1, Read &mate2,
                                   vector<Region> &region_f1, vector<Region> &region_r1,
                                   vector<Region> &region_f2, vector<Region> &region_r2,
                                   unsigned &best_f1, unsigned &best_r1, unsigned &best_f2, unsigned &best_r2,
                                   unsigned &next_f1, unsigned &next_r1, unsigned &next_f2, unsigned &next_r2,
                                   bool *&flag_f1, bool *&flag_r1, bool *&flag_f2, bool *&flag_r2,
                                   bool &has_f1r2, bool &has_r1f2) {
  int min_rlen = strlen(mate1.seq) < strlen(mate2.seq) ? strlen(mate1.seq) : strlen(mate2.seq);
  unsigned slide = kmer_len < min_rlen - kmer_len ? kmer_len : min_rlen - kmer_len;
  unsigned kmer_step1 = kmer_len, kmer_step2 = kmer_len;
  unsigned slide1 = 0, slide2 = 0;

  bool high_freq_1 = false, high_freq_2 = false; //read is from high repetitive region
  int mac_occ_1 = MAX_OCC, mac_occ_2 = MAX_OCC;

  while (slide1 < slide && slide2 < slide) {
//  while (kmer_step1 > 0 && kmer_step2 > 0) {
    if (has_f1r2 || has_r1f2)
      break;

    pghole_wrapper_mates(mate1, region_f1, region_r1, best_f1, best_r1, slide1, kmer_step1, mac_occ_1, high_freq_1);
    pghole_wrapper_mates(mate2, region_f2, region_r2, best_f2, best_r2, slide2, kmer_step2, mac_occ_2, high_freq_2);

    // filter based on pairdis
    flag_f1 = new bool[region_f1.size()]();
    flag_r1 = new bool[region_r1.size()]();
    flag_f2 = new bool[region_f2.size()]();
    flag_r2 = new bool[region_r2.size()]();
    has_f1r2 = pairdis_filter(region_f1, region_r2, flag_f1, flag_r2, best_f1, next_f1, best_r2, next_r2);
    has_r1f2 = pairdis_filter(region_r1, region_f2, flag_r1, flag_f2, best_r1, next_r1, best_f2, next_f2);

    if (!has_f1r2 && !has_r1f2) {
      region_f1.clear();
      region_f2.clear();
      region_r1.clear();
      region_r2.clear();
      delete[] flag_f1;
      delete[] flag_r1;
      delete[] flag_f2;
      delete[] flag_r2;
    }

    slide1++;
    slide2++;
//    kmer_step1 /= 2;
//    kmer_step2 /= 2;
  }
}

void AccAlign::embed_wrapper_pair(Read &R1, Read &R2,
                                  vector<Region> &candidate_regions_f1, vector<Region> &candidate_regions_r2,
                                  bool flag_f1[], bool flag_r2[], unsigned &best_f1, unsigned &best_r2,
                                  int &best_threshold, int &next_threshold, char strand) {
  const char *ptr_ref = ref.c_str();
  char *seq1, *seq2;
  if (strand == '+') {
    //f1r2
    seq1 = R1.fwd;
    seq2 = R2.rev;
  } else {
    //r1f2
    seq1 = R1.rev;
    seq2 = R2.fwd;
  }

  //embed r1
  embedding->embed_unmatch(candidate_regions_f1, ptr_ref, seq1, strlen(R1.seq), R1.kmer_step, flag_f1);
  //embed r2
  embedding->embed_unmatch_pair(R1, R2, candidate_regions_f1, candidate_regions_r2, ptr_ref, seq2, strlen(R2.seq),
                                R2.kmer_step, flag_r2, pairdis, best_threshold, next_threshold, best_f1, best_r2);
}

void AccAlign::embed_wrapper(Read &R, bool ispe,
                             vector<Region> &fcandidate_regions, vector<Region> &rcandidate_regions,
                             unsigned &fbest, unsigned &fnext, unsigned &rbest, unsigned &rnext,
                             int &best_threshold, int &next_threshold) {
  unsigned nfregions = fcandidate_regions.size();
  unsigned nrregions = rcandidate_regions.size();
  assert(nfregions + nrregions > 0); //at least one hit

  // for single-end mapping, we better have things sorted by coverage. for
  // paired end, it is possible for a region with high cov to be eliminated if
  // no pairs are found
  if (!ispe && nfregions > 1)
    assert(fcandidate_regions[0].cov >= fcandidate_regions[1].cov);
  if (!ispe && nrregions > 1)
    assert(rcandidate_regions[0].cov >= rcandidate_regions[1].cov);
  vpair_sort_count += nfregions + nrregions;

  // pass the one with the highest coverage in first
  bool fwd_first = true;
  // there is no fwd, or fwd cov < rev cov
  if (nfregions == 0 || (nrregions > 0 && fcandidate_regions[0].cov < rcandidate_regions[0].cov))
    fwd_first = false;

  // embed now, but only in the case where we have > 1 regions either for the
  // forward or for reverse or for both strands. If we have only 1 region
  // globally, there is no point embedding.
  size_t rlen = strlen(R.seq);
  fbest = fnext = rbest = rnext = 0;
  const char *ptr_ref = ref.c_str();

  if (fwd_first) {
    embedding->embed_unmatch_iter(fcandidate_regions, ptr_ref, R.fwd, rlen, R.kmer_step,
                                  best_threshold, next_threshold, fbest, fnext);
    if (nrregions)
      embedding->embed_unmatch_iter(rcandidate_regions, ptr_ref, R.rev, rlen, R.kmer_step,
                                    best_threshold, next_threshold, rbest, rnext);
  } else {
    embedding->embed_unmatch_iter(rcandidate_regions, ptr_ref, R.rev, rlen, R.kmer_step,
                                  best_threshold, next_threshold, rbest, rnext);
    if (nfregions)
      embedding->embed_unmatch_iter(fcandidate_regions, ptr_ref, R.fwd, rlen, R.kmer_step,
                                    best_threshold, next_threshold, fbest, fnext);
  }

}

int AccAlign::get_mapq(int best, int secBest) {
  int mapq;
  if (secBest == 0 && best == 0) {
    mapq = 40;
  } else {
    float x = (float) best / secBest;
    mapq = ceil(60 * (1 - x * x));
  }
  mapq = mapq < 60 ? mapq : 60;
  mapq = mapq < 0 ? 0 : mapq;
  return mapq;
}

void AccAlign::map_read(Read &R) {

  if (strlen(R.seq) < kmer_len){
    R.strand = '*';
    return;
  }

  auto start = std::chrono::system_clock::now();
  parse(R.seq, R.fwd, R.rev, R.rev_str);
  auto end = std::chrono::system_clock::now();
  auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
  parse_time += elapsed.count();

  start = std::chrono::system_clock::now();
  vector<Region> fcandidate_regions, rcandidate_regions;

  // first try pigenhole. if we generate any candidates, pass them through
  // lsh. If something comes out, we are done.
  // XXX: On experimentation, it was found that using pigeonhole filtering
  // produces wrong results and invalid mappings when errors are too large.
  unsigned fbest, rbest;
  pghole_wrapper(R, fcandidate_regions, rcandidate_regions, fbest, rbest);
  unsigned nfregions = fcandidate_regions.size();
  unsigned nrregions = rcandidate_regions.size();
  end = std::chrono::system_clock::now();
  elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
  seeding_time += elapsed.count();

  if (nfregions == 0 && nrregions == 0) {
    R.strand = '*';
    return;
  }

  if (extend_all) {
    int max_as = INT_MIN;
    char strand = '+';
    Region r;
    for (Region &region: fcandidate_regions) {
      char *s = R.fwd;
      Alignment a;
      score_region(R, s, region, a);
      if (region.score >= max_as) {
        r = region;
        max_as = region.score;
      }
    }

    for (Region &region: rcandidate_regions) {
      char *s = R.rev;
      Alignment a;
      score_region(R, s, region, a);
      if (region.score >= max_as) {
        r = region;
        strand = '-';
        max_as = region.score;
      }
    }

    R.best_region = r;
    R.strand = strand;
  } else {
    start = std::chrono::system_clock::now();
    // before doing embedding, move highest cov region to front
    if (nfregions > 1 && fbest != 0) {
      iter_swap(fcandidate_regions.begin() + fbest, fcandidate_regions.begin());
    }
    if (nrregions > 1 && rbest != 0) {
      iter_swap(rcandidate_regions.begin() + rbest, rcandidate_regions.begin());
    }
    end = std::chrono::system_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    swap_time += elapsed.count();
    seeding_time += elapsed.count();

    unsigned fnext, rnext;
    int best_threshold = strlen(R.seq) * embedding->efactor;
    int next_threshold = strlen(R.seq) * embedding->efactor;
    embed_wrapper(R, false, fcandidate_regions, rcandidate_regions, fbest, fnext, rbest, rnext,
                  best_threshold, next_threshold);

    start = std::chrono::system_clock::now();

    if (nfregions == 0) {
      if (nrregions == 0)
        return;
      mark_for_extension(R, '-', rcandidate_regions[rbest]);
    } else if (nrregions == 0) {
      mark_for_extension(R, '+', fcandidate_regions[fbest]);
    } else {
      // pick the candidate with smallest embed dist
      // if fwd/rev have same embed_dist, take the hcov one
      // if hcov one not the min dist, take the one with smaller pos (to be consistent with gpu)
      if (fcandidate_regions[fbest].embed_dist < rcandidate_regions[rbest].embed_dist) {
        mark_for_extension(R, '+', fcandidate_regions[fbest]);
      } else if (fcandidate_regions[fbest].embed_dist > rcandidate_regions[rbest].embed_dist) {
        mark_for_extension(R, '-', rcandidate_regions[rbest]);
      } else {
        if (fcandidate_regions[fbest].rs < rcandidate_regions[rbest].rs) {
          mark_for_extension(R, '+', fcandidate_regions[fbest]);
        } else {
          mark_for_extension(R, '-', rcandidate_regions[rbest]);
        }
      }
    }

    R.best = best_threshold;
    R.secBest = next_threshold;

    end = std::chrono::system_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    mapqTime += elapsed.count();
  }

}

void AccAlign::extend_pair(Read &mate1, Read &mate2,
                           vector<Region> &candidate_regions_f1, vector<Region> &candidate_regions_r2,
                           bool flag_f1[], bool flag_r2[], unsigned &best_f1, unsigned &best_r2,
                           int &best_threshold, int &next_threshold, char strand) {
  char *seq1, *seq2;
  if (strand == '+') { //f1r2
    seq1 = mate1.fwd;
    seq2 = mate2.rev;
  } else { //r1f2
    seq1 = mate1.rev;
    seq2 = mate2.fwd;
  }

  //extend mate1
  for (unsigned i = 0; i < candidate_regions_f1.size(); ++i) {
    if (!flag_f1[i]) {
      continue;
    }

    Region &region = candidate_regions_f1[i];
    Alignment a;
    int rs_bk = region.rs;
    score_region(mate1, seq1, region, a);
    region.rs = rs_bk; // sore_region will justify the start pos, convert back, because the later itr need it
  }

  //extend mate2
  for (unsigned i = 0; i < candidate_regions_r2.size(); ++i) {
    if (!flag_r2[i]) {
      continue;
    }

    Region &region = candidate_regions_r2[i];
    Alignment a;
    int rs_bk = region.rs;
    score_region(mate2, seq2, region, a);
    region.rs = rs_bk; // sore_region will justify the start pos, convert back, because the later itr need it

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
    for (auto itr = start; itr != end; ++itr) {
      int sum_as = region.score + itr->score;
      if (sum_as >= best_threshold) {
        best_f1 = itr - candidate_regions_f1.begin();
        best_r2 = i;
        best_threshold = sum_as;
      } else if (sum_as >= next_threshold) {
        next_threshold = sum_as;
      }
    }
  }

//  best_threshold = -best_threshold;  //get the neg of as, so that lower score means better, consist
//  next_threshold = -next_threshold;
}

void AccAlign::map_paired_read(Read &mate1, Read &mate2) {

  if (strlen(mate1.seq) < kmer_len || strlen(mate2.seq) < kmer_len){
    mate1.strand = '*';
    mate2.strand = '*';
    return;
  }

  auto start = std::chrono::system_clock::now();
  parse(mate1.seq, mate1.fwd, mate1.rev, mate1.rev_str);
  parse(mate2.seq, mate2.fwd, mate2.rev, mate2.rev_str);
  auto end = std::chrono::system_clock::now();
  auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
  parse_time += elapsed.count();

  start = std::chrono::system_clock::now();
  vector<Region> region_f1, region_r1, region_f2, region_r2;
  unsigned best_f1 = 0, best_r1 = 0, best_f2 = 0, best_r2 = 0;
  unsigned next_f1 = 0, next_r1 = 0, next_f2 = 0, next_r2 = 0;
  bool *flag_f1 = nullptr, *flag_r1 = nullptr, *flag_f2 = nullptr, *flag_r2 = nullptr;
  bool has_f1r2 = false, has_r1f2 = false;

  // lookup candidates
  unsigned min_rlen = strlen(mate1.seq) < strlen(mate2.seq) ? strlen(mate1.seq) : strlen(mate2.seq);
  if (min_rlen < kmer_len) {
    mate1.strand = '*';
    mate2.strand = '*';
    return;
  }

  pghole_wrapper_pair(mate1, mate2, region_f1, region_r1, region_f2, region_r2,
                      best_f1, best_r1, best_f2, best_r2, next_f1, next_r1, next_f2, next_r2,
                      flag_f1, flag_r1, flag_f2, flag_r2, has_f1r2, has_r1f2);
  end = std::chrono::system_clock::now();
  elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
  seeding_time += elapsed.count();

  if (!has_f1r2 && !has_r1f2) {
    mate1.strand = '*';
    mate2.strand = '*';
    return;
  }

  if (extend_all) {
    int best_f1r2, next_f1r2, best_r1f2, next_r1f2;
    best_f1r2 = next_f1r2 = best_r1f2 = next_r1f2 = INT_MIN;

    if (has_f1r2)
      extend_pair(mate1, mate2, region_f1, region_r2, flag_f1, flag_r2, best_f1, best_r2, best_f1r2, next_f1r2, '+');
    if (has_r1f2)
      extend_pair(mate1, mate2, region_r1, region_f2, flag_r1, flag_f2, best_r1, best_f2, best_r1f2, next_r1f2, '-');

    // if there is no candidates, the strand will remain *
    int secmin_dist;
    if (best_f1r2 > best_r1f2) {
      mark_for_extension(mate1, '+', region_f1[best_f1]);
      mark_for_extension(mate2, '-', region_r2[best_r2]);
      if (best_r1f2 > next_f1r2)
        secmin_dist = best_r1f2;
      else
        secmin_dist = next_f1r2;
    } else {
      mark_for_extension(mate2, '+', region_f2[best_f2]);
      mark_for_extension(mate1, '-', region_r1[best_r1]);
      if (best_f1r2 > next_r1f2)
        secmin_dist = best_f1r2;
      else
        secmin_dist = next_r1f2;
    }

    mate1.best = mate2.best = min(best_f1r2, best_r1f2);
    mate1.secBest = mate2.secBest = secmin_dist;
    return;
  }

  // now apply embedding filter on filtered regions.
  // But before, rearrange so that regions with high coverage are at the top
  //no need to swap, just embed the best and next first..
  int best_f1r2 = INT_MAX, next_f1r2 = INT_MAX, best_r1f2 = INT_MAX, next_r1f2 = INT_MAX;
  if (has_f1r2)
    embed_wrapper_pair(mate1, mate2, region_f1, region_r2, flag_f1, flag_r2,
                       best_f1, best_r2, best_f1r2, next_f1r2, '+');
  if (has_r1f2)
    embed_wrapper_pair(mate1, mate2, region_r1, region_f2, flag_r1, flag_f2,
                       best_r1, best_f2, best_r1f2, next_r1f2, '-');

  delete[] flag_f1;
  delete[] flag_r1;
  delete[] flag_f2;
  delete[] flag_r2;

  start = std::chrono::system_clock::now();
  if (best_f1r2 <= best_r1f2) {
    mark_for_extension(mate1, '+', region_f1[best_f1]);
    mark_for_extension(mate2, '-', region_r2[best_r2]);
  } else {
    mark_for_extension(mate2, '+', region_f2[best_f2]);
    mark_for_extension(mate1, '-', region_r1[best_r1]);
  }

  end = std::chrono::system_clock::now();
  elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
  mapqTime += elapsed.count();
}

void AccAlign::print_paired_sam(Read &R, Read &R2) {
  auto start = std::chrono::system_clock::now();

  stringstream ss;
  char strand1 = R.strand;
  char strand2 = R2.strand;

  //mate 1
  string rname = R.name;
  ss << rname.substr(0, rname.find_last_of("/")) << '\t';

  uint16_t flag = 0x1;
  if (strand1 == '*')
    flag |= 0x4;
  if (strand2 == '*')
    flag |= 0x8;
  if (!(flag & 0x4) && !(flag & 0x8))
    flag |= 0x2;
  if (strand1 == '-')
    flag |= 0x10;
  if (strand2 == '-')
    flag |= 0x20;
  flag |= 0x40;
  ss << flag;

  ss << '\t' << (strand1 == '*' ? "*" : name[R.tid]);
  ss << '\t' << (strand1 == '*' ? 0 : R.pos);
  ss << '\t' << (strand1 == '*' ? 0 : (int) R.mapq);
  ss << '\t' << (strand1 == '*' ? "*" : R.cigar) << '\t';

  if (R.strand == '*' || R2.strand == '*')
    ss << '*';
  else
    ss << (name[R.tid] == name[R2.tid] ? "=" : name[R2.tid].c_str());

  ss << '\t' << (strand2 == '*' ? 0 : R2.pos);

  int isize = 0;
  if (R.strand != '*' && R2.strand != '*') {
    if (R.pos > R2.pos)
      isize = R2.pos - R.pos - strlen(R.seq);
    else
      isize = R2.pos - R.pos + strlen(R2.seq);
  }
  ss << '\t' << isize;

  if (strand1 == '-') {
    std::reverse(R.qua, R.qua + strlen(R.qua));
    ss << '\t' << R.rev_str;
  } else {
    ss << '\t' << R.seq;
  }
  ss << '\t' << R.qua << "\tNM:i:" << R.nm << "\tAS:i:" << R.as << endl;


  //for mate2
  string rname2 = R2.name;
  ss << rname2.substr(0, rname2.find_last_of("/")) << '\t';

  flag = 0x1;
  if (strand2 == '*')
    flag |= 0x4;
  if (strand1 == '*')
    flag |= 0x8;
  if (!(flag & 0x4) && !(flag & 0x8))
    flag |= 0x2;
  if (strand2 == '-')
    flag |= 0x10;
  if (strand1 == '-')
    flag |= 0x20;
  flag |= 0x80;

  ss << flag;
  ss << '\t' << (strand2 == '*' ? "*" : name[R2.tid]);
  ss << '\t' << (strand2 == '*' ? 0 : R2.pos);
  ss << '\t' << (strand2 == '*' ? 0 : (int) R2.mapq);
  ss << '\t' << (strand2 == '*' ? "*" : R2.cigar) << '\t';

  if (R.strand == '*' || R2.strand == '*')
    ss << '*';
  else
    ss << (name[R.tid] == name[R2.tid] ? "=" : name[R.tid].c_str());

  ss << '\t' << (strand2 == '*' ? 0 : R.pos);
  ss << '\t' << -isize;
  if (strand2 == '-') {
    std::reverse(R2.qua, R2.qua + strlen(R2.qua));
    ss << '\t' << R2.rev_str;
  } else {
    ss << '\t' << R2.seq;
  }
  ss << '\t' << R2.qua << "\tNM:i:" << R2.nm << "\tAS:i:" << R2.as << endl;

  {
    std::lock_guard<std::mutex> guard(sam_mutex);
    if (sam_name.length()) {
      sam_stream << ss.str();
      sam_stream.flush();
    } else {
      cout << ss.str();
    }
  }

  auto end = std::chrono::system_clock::now();
  auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
  sam_time += elapsed.count();
}

void AccAlign::snprintf_pair_sam(Read &R, string *s, Read &R2, string *s2) {
  auto start = std::chrono::system_clock::now();

  // 60 is the approximate length for all int
  int size;
  if (!enable_extension) {
    size = 60;
  } else {
    size = 60 + strlen(R.seq); //assume the length of cigar will not longer than the read
  }
  char strand1 = R.strand;
  char strand2 = R2.strand;

  //mate 1
  string rname = R.name;
  string nn = rname.substr(0, rname.find_last_of("/"));

  uint16_t flag = 0x1;
  if (strand1 == '*')
    flag |= 0x4;
  if (strand2 == '*')
    flag |= 0x8;
  if (!(flag & 0x4) && !(flag & 0x8))
    flag |= 0x2;
  if (strand1 == '-')
    flag |= 0x10;
  if (strand2 == '-')
    flag |= 0x20;
  flag |= 0x40;

  int isize = 0;
  if (R.strand != '*' && R2.strand != '*') {
    if (R.pos > R2.pos)
      isize = R2.pos - R.pos - strlen(R.seq);
    else
      isize = R2.pos - R.pos + strlen(R2.seq);
  }

  string format = "%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s\tNM:i:%d\tAS:i:%d\n";
  if (R.strand == '+' && R2.strand != '*') {
    size += strlen(R.name) + name[R.tid].length() + name[R2.tid].length() + 2 * strlen(R.seq);
    char buf[size];
    snprintf(buf, size, format.c_str(), nn.c_str(), flag, name[R.tid].c_str(), R.pos,
             (int) R.mapq, R.cigar, name[R.tid] == name[R2.tid] ? "=" : name[R2.tid].c_str(),
             R2.pos, isize, R.seq, R.qua, R.nm, R.as);
    *s = buf;
  } else if (R.strand == '-' && R2.strand != '*') {
    size += strlen(R.name) + name[R.tid].length() + name[R2.tid].length() + 2 * strlen(R.seq);
    char buf[size];
    std::reverse(R.qua, R.qua + strlen(R.qua));
    snprintf(buf, size, format.c_str(), nn.c_str(), flag, name[R.tid].c_str(), R.pos,
             (int) R.mapq, R.cigar, name[R.tid] == name[R2.tid] ? "=" : name[R2.tid].c_str(),
             R2.pos, isize, R.rev_str, R.qua, R.nm, R.as);
    *s = buf;
  } else if (R.strand == '+' && R2.strand == '*') {
    size += strlen(R.name) + name[R.tid].length() + 2 * strlen(R.seq);
    char buf[size];
    snprintf(buf, size, format.c_str(),
             nn.c_str(), flag, name[R.tid].c_str(), R.pos, (int) R.mapq, R.cigar, "*", 0,
             isize, R.seq, R.qua, R.nm, R.as);
    *s = buf;
  } else if (R.strand == '-' && R2.strand == '*') {
    size += strlen(R.name) + name[R2.tid].length() + 2 * strlen(R.seq);
    char buf[size];
    std::reverse(R.qua, R.qua + strlen(R.qua));
    snprintf(buf, size, format.c_str(),
             nn.c_str(), flag, name[R.tid].c_str(), R.pos, (int) R.mapq, R.cigar, "*", 0,
             isize, R.rev_str, R.qua, R.nm, R.as);
    *s = buf;
  } else if (R.strand == '*' && R2.strand != '*') {
    size += strlen(R.name) + name[R2.tid].length() + 2 * strlen(R.seq);
    char buf[size];
    snprintf(buf, size, format.c_str(),
             nn.c_str(), flag, "*", 0, 0, "*", name[R2.tid].c_str(), R2.pos,
             isize, R.seq, R.qua, R.nm, R.as);
    *s = buf;
  } else if (R.strand == '*' && R2.strand == '*') {
    size += strlen(R.name) + 2 * strlen(R.seq);
    char buf[size];
    snprintf(buf, size, format.c_str(),
             nn.c_str(), flag, "*", 0, 0, "*", "*", 0,
             isize, R.seq, R.qua, R.nm, R.as);
    *s = buf;
  }

  //for mate2
  string rname2 = R2.name;
  string nn2 = rname2.substr(0, rname2.find_last_of("/"));

  flag = 0x1;
  if (strand2 == '*')
    flag |= 0x4;
  if (strand1 == '*')
    flag |= 0x8;
  if (!(flag & 0x4) && !(flag & 0x8))
    flag |= 0x2;
  if (strand2 == '-')
    flag |= 0x10;
  if (strand1 == '-')
    flag |= 0x20;
  flag |= 0x80;

  if (R2.strand == '+' && R.strand != '*') {
    size += strlen(R2.name) + name[R2.tid].length() + name[R.tid].length() + 2 * strlen(R2.seq);
    char buf[size];
    snprintf(buf, size, format.c_str(), nn2.c_str(), flag, name[R2.tid].c_str(), R2.pos,
             (int) R2.mapq, R2.cigar, name[R.tid] == name[R2.tid] ? "=" : name[R.tid].c_str(),
             R.pos, -isize, R2.seq, R2.qua, R2.nm, R2.as);
    *s2 = buf;
  } else if (R2.strand == '-' && R.strand != '*') {
    size += strlen(R2.name) + name[R2.tid].length() + name[R.tid].length() + 2 * strlen(R2.seq);
    char buf[size];
    std::reverse(R2.qua, R2.qua + strlen(R2.qua));
    snprintf(buf, size, format.c_str(), nn2.c_str(), flag, name[R2.tid].c_str(), R2.pos,
             (int) R2.mapq, R2.cigar, name[R.tid] == name[R2.tid] ? "=" : name[R.tid].c_str(),
             R.pos, -isize, R2.rev_str, R2.qua, R2.nm, R2.as);
    *s2 = buf;
  } else if (R2.strand == '+' && R.strand == '*') {
    size += strlen(R2.name) + name[R2.tid].length() + 2 * strlen(R2.seq);
    char buf[size];
    snprintf(buf, size, format.c_str(),
             nn2.c_str(), flag, name[R2.tid].c_str(), R2.pos, (int) R2.mapq, R2.cigar, "*", 0,
             -isize, R2.seq, R2.qua, R2.nm, R2.as);
    *s2 = buf;
  } else if (R2.strand == '-' && R.strand == '*') {
    size += strlen(R2.name) + name[R.tid].length() + 2 * strlen(R2.seq);
    char buf[size];
    std::reverse(R2.qua, R2.qua + strlen(R2.qua));
    snprintf(buf, size, format.c_str(),
             nn2.c_str(), flag, name[R2.tid].c_str(), R2.pos, (int) R2.mapq, R2.cigar, "*", 0,
             -isize, R2.rev_str, R2.qua, R2.nm, R2.as);
    *s2 = buf;
  } else if (R2.strand == '*' && R.strand != '*') {
    size += strlen(R2.name) + name[R.tid].length() + 2 * strlen(R2.seq);
    char buf[size];
    snprintf(buf, size, format.c_str(),
             nn2.c_str(), flag, "*", 0, 0, "*", name[R.tid].c_str(), R.pos,
             -isize, R2.seq, R2.qua, R2.nm, R2.as);
    *s2 = buf;
  } else if (R2.strand == '*' && R.strand == '*') {
    size += strlen(R2.name) + 2 * strlen(R2.seq);
    char buf[size];
    snprintf(buf, size, format.c_str(),
             nn2.c_str(), flag, "*", 0, 0, "*", "*", 0,
             -isize, R2.seq, R2.qua, R2.nm, R2.as);
    *s2 = buf;
  }

  auto end = std::chrono::system_clock::now();
  auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
  sam_pre_time += elapsed.count();
}

void AccAlign::print_sam(Read &R) {
  auto start = std::chrono::system_clock::now();

  stringstream ss;
  string rname = R.name;
  ss << rname.substr(0, rname.find_last_of("/")) << '\t';

  if (R.strand == '*') {
    ss << "4\t*\t0\t255\t*\t*\t0\t0\t";
  } else {
    ss << (R.strand == '+' ? 0 : 16) << "\t" <<
       name[R.tid] << "\t" <<
       R.pos << "\t" <<
       (int) R.mapq << "\t" <<
       R.cigar << "\t*\t0\t0\t";
  }

  if (R.strand == '-') {
    std::reverse(R.qua, R.qua + strlen(R.qua));
    ss << '\t' << R.rev_str;
  } else {
    ss << '\t' << R.seq;
  }
  ss << '\t' << R.qua << "\tNM:i:" << R.nm << "\tAS:i:" << R.as << endl;

  {
    std::lock_guard<std::mutex> guard(sam_mutex);
    if (sam_name.length()) {
      sam_stream << ss.str();
      sam_stream.flush();
    } else {
      cout << ss.str();
    }
  }

  auto end = std::chrono::system_clock::now();
  auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
  sam_pre_time += elapsed.count();
}

void AccAlign::snprintf_sam(Read &R, string *s) {
  auto start = std::chrono::system_clock::now();

  // 50 is the approximate length for all int
  int size;
  if (!enable_extension) {
    size = 50;
  } else {
    size = 50 + strlen(R.seq); //assume the length of cigar will not longer than the read
  }

  string format = "%s\t%d\t%s\t%d\t%d\t%s\t*\t0\t0\t%s\t%s\tNM:i:%d\tAS:i:%d\n";
  if (R.strand == '*') {
    size += strlen(R.name) + 2 * strlen(R.seq);
    char buf[size];
    snprintf(buf, size, format.c_str(),
             R.name, 0, "*", 0, 0, "*", R.seq, R.qua, 0, 0);
    *s = buf;
  } else {
    size += strlen(R.name) + name[R.tid].length() + 2 * strlen(R.seq);
    char buf[size];
    if (R.strand == '+') {
      snprintf(buf, size, format.c_str(),
               R.name, 0, name[R.tid].c_str(), R.pos, (int) R.mapq, R.cigar,
               R.seq, R.qua, R.nm, R.as);
    } else {
      std::reverse(R.qua, R.qua + strlen(R.qua));
      snprintf(buf, size, format.c_str(),
               R.name, 16, name[R.tid].c_str(), R.pos, (int) R.mapq, R.cigar,
               R.rev_str, R.qua, R.nm, R.as);
    }
    *s = buf;
  }

  auto end = std::chrono::system_clock::now();
  auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
  sam_pre_time += elapsed.count();
}

void AccAlign::out_sam(string *s) {
  auto start = std::chrono::system_clock::now();
  {
    if (sam_name.length()) {
      sam_stream << *s;
    } else {
      cout << *s;
    }
  }
  auto end = std::chrono::system_clock::now();
  auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
  sam_out_time += elapsed.count();
}

class Tbb_aligner {
  Read *all_reads;
  string *sams;
  AccAlign *acc_obj;

 public:
  Tbb_aligner(Read *_all_reads, string *_sams, AccAlign *_acc_obj) :
      all_reads(_all_reads), sams(_sams), acc_obj(_acc_obj) {}

  void operator()(const tbb::blocked_range<size_t> &r) const {
    for (size_t i = r.begin(); i != r.end(); ++i) {
      if ((all_reads + i)->strand != '*') {
        if (enable_wfa_extension)
          acc_obj->wfa_align_read(*(all_reads + i));
        else
          acc_obj->align_read(*(all_reads + i));
      }
      acc_obj->snprintf_sam(*(all_reads + i), sams + i);
    }
  }
};

class Tbb_aligner_paired {
  Read *all_reads;
  Read *all_reads2;
  string *sams;
  AccAlign *acc_obj;

 public:
  Tbb_aligner_paired(Read *_all_reads, Read *_all_reads2, string *_sams, AccAlign *_acc_obj) :
      all_reads(_all_reads), all_reads2(_all_reads2), sams(_sams), acc_obj(_acc_obj) {}

  void operator()(const tbb::blocked_range<size_t> &r) const {
    for (size_t i = r.begin(); i != r.end(); ++i) {
      if ((all_reads + i)->strand != '*') {
        if (enable_wfa_extension)
          acc_obj->wfa_align_read(*(all_reads + i));
        else
          acc_obj->align_read(*(all_reads + i));
      }

      if ((all_reads2 + i)->strand != '*') {
        if (enable_wfa_extension)
          acc_obj->wfa_align_read(*(all_reads2 + i));
        else
          acc_obj->align_read(*(all_reads2 + i));
      }

      acc_obj->snprintf_pair_sam(*(all_reads + i), sams + 2 * i, *(all_reads2 + i), sams + 2 * i + 1);
    }
  }
};

void AccAlign::align_wrapper(int tid, int soff, int eoff, Read *ptlread, Read *ptlread2,
                             tbb::concurrent_bounded_queue<ReadPair> *dataQ) {

  if (!ptlread2) {
    // single-end read alignment
    string sams[eoff];
    tbb::task_scheduler_init init(g_ncpus);
    tbb::parallel_for(tbb::blocked_range<size_t>(soff, eoff), Tbb_aligner(ptlread, sams, this));

    auto start = std::chrono::system_clock::now();
    for (int i = soff; i < eoff; i++) {
      out_sam(sams + i);
    }
    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    sam_time += elapsed.count();

    dataQ->push(make_tuple(ptlread, (Read *) NULL));
  }
//  else {
//    string sams[2 * eoff];
//    tbb::task_scheduler_init init(g_ncpus);
//    tbb::parallel_for(tbb::blocked_range<size_t>(soff, eoff), Tbb_aligner_paired(ptlread, ptlread2, sams, this));
//
//    auto start = std::chrono::system_clock::now();
//    for (int i = soff; i < 2 * eoff; i++) {
//      out_sam(sams + i);
//    }
//    auto end = std::chrono::system_clock::now();
//    auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
//    sam_time += elapsed.count();
//
//    dataQ->push(make_tuple(ptlread, ptlread2));
//  }
}

static void ksw_gen_simple_mat(int m, int8_t *mat, int8_t a, int8_t b, int8_t sc_ambi) {
  int i, j;
  a = a < 0 ? -a : a;
  b = b > 0 ? -b : b;
  sc_ambi = sc_ambi > 0 ? -sc_ambi : sc_ambi;
  for (i = 0; i < m - 1; ++i) {
    for (j = 0; j < m - 1; ++j)
      mat[i * m + j] = i == j ? a : b;
    mat[i * m + m - 1] = sc_ambi;
  }
  for (j = 0; j < m; ++j)
    mat[(m - 1) * m + j] = sc_ambi;
}

//rectify_start_pos at the end of map
void AccAlign::rectify_start_pos(char *strand, Region &region, unsigned rlen) {
  //embed first kmer of read
  int elen = kmer_len * embedding->efactor;
  char embeddedQ[elen];
  embedding->cgk2_embedQ(strand, kmer_len, 0, embeddedQ);

  //shift the start pos and find the best
  const char *ptr_ref = ref.c_str();

  // the pos without shift
  int threshold = embedding->cgk2_embed_nmismatch(ptr_ref + region.rs, kmer_len, elen, 0, embeddedQ);
  int shift = 0;

  float indel_len = ceil(MAX_INDEL * rlen / float(100));
  for (int i = -indel_len; i < indel_len; ++i) {
    if (i == 0 || region.rs + i < 0 || region.rs + i >= ref.size())
      continue;

    int nmismatch = embedding->cgk2_embed_nmismatch(ptr_ref + region.rs + i, kmer_len, threshold, 0, embeddedQ);
    if (nmismatch < threshold) {
      threshold = nmismatch;
      shift = i;
    }
  }

  region.rs = region.rs + shift;
}

static inline void mm_seq_rev(uint32_t len, char *seq) {
  uint32_t i;
  char t;
  for (i = 0; i < len >> 1; ++i)
    t = seq[i], seq[i] = seq[len - 1 - i], seq[len - 1 - i] = t;
}

static void append_cigar(Extension *&rp, uint32_t n_cigar, uint32_t *cigar) {
  if (n_cigar == 0) return;
  if (rp == 0) {
    uint32_t capacity = n_cigar + sizeof(Extension) / 4;
    kroundup32(capacity);
    rp = (Extension *) calloc(capacity, 4);
    rp->capacity = capacity;
  } else if (rp->n_cigar + n_cigar + sizeof(Extension) / 4 > rp->capacity) {
    rp->capacity = rp->n_cigar + n_cigar + sizeof(Extension) / 4;
    kroundup32(rp->capacity);
    rp = (Extension *) realloc(rp, rp->capacity * 4);
  }
  if (rp->n_cigar > 0 && (rp->cigar[rp->n_cigar - 1] & 0xf) == (cigar[0] & 0xf)) { // same CIGAR op at the boundary
    rp->cigar[rp->n_cigar - 1] += cigar[0] >> 4 << 4;
    if (n_cigar > 1) memcpy(rp->cigar + rp->n_cigar, cigar + 1, (n_cigar - 1) * 4);
    rp->n_cigar += n_cigar - 1;
  } else {
    memcpy(rp->cigar + rp->n_cigar, cigar, n_cigar * 4);
    rp->n_cigar += n_cigar;
  }
}

void AccAlign::score_region(Read &r, char *qseq, Region &region,
                            Alignment &a) {
  unsigned qlen = strlen(r.seq);

  if (!enable_extension) {
    region.score = qlen * SC_MCH;
    r.mapq = get_mapq(r.best, r.secBest);
  } else if (!extend_all && (region.embed_dist == 0 || region.embed_dist == 1)) {
    // if the region has a embed distance of 0, then its an exact match
    if (region.embed_dist == 0)
      region.score = qlen * SC_MCH;
    if (region.embed_dist == 1)
      region.score = (qlen - 1) * SC_MCH - SC_MIS;

    r.mapq = get_mapq(r.best, r.secBest);
  } else {
    Extension *extension = nullptr;
    uint32_t raw_rs = region.rs;

    uint32_t qs, qe, qs0, qe0, rs, re, rs0, re0, l;
    uint32_t beginclip = 0, endclip = 0;
    qs = region.qs;
    assert(region.qs == region.matched_intervals[0]);
    qe = qs + kmer_len;
    rs = region.rs + qs;
    re = region.rs + qe;
    l = qs;
    l += l * SC_MCH + END_BONUS > GAPO ? (l * SC_MCH + END_BONUS - GAPO) / GAPE : 0;
    qs0 = 0, qe0 = qlen;
    rs0 = rs > l ? rs - l : 0;

    l = qlen - qe;
    l += l * SC_MCH + END_BONUS > GAPO ? (l * SC_MCH + END_BONUS - GAPO) / GAPE : 0;
    re0 = re + l < offset.back() ? re + l: offset.back();

    //left extension
    if (qs > 0 && rs > 0) {
      ksw_extz_t ez_l;
      memset(&ez_l, 0, sizeof(ksw_extz_t));
      char tseq[rs - rs0];
      memcpy(tseq, ref.c_str() + rs0, rs - rs0);
      mm_seq_rev(rs - rs0, tseq);
      mm_seq_rev(qs - qs0, qseq);

      const uint8_t *_tseq = reinterpret_cast<const uint8_t *>(tseq);
      const uint8_t *_qseq = reinterpret_cast<const uint8_t *>(qseq);
      ksw_extz2_sse(0, qs - qs0, _qseq, rs - rs0, _tseq, 5, mat, GAPO, GAPE, BANDWIDTH,
                    Z_DROP, END_BONUS, KSW_EZ_EXTZ_ONLY | KSW_EZ_RIGHT | KSW_EZ_REV_CIGAR, &ez_l);

      if (ez_l.n_cigar > 0) {
        append_cigar(extension, ez_l.n_cigar, ez_l.cigar);
        extension->dp_score += ez_l.max;
      }
      beginclip = ez_l.reach_end ? 0 : qs - qs0 - ez_l.max_q - 1;
      region.rs = rs - (ez_l.reach_end ? ez_l.mqe_t + 1 : ez_l.max_t + 1);
      mm_seq_rev(qs - qs0, qseq);
      free(ez_l.cigar);
    }

    // matched seed
    uint32_t cigar_m[] = {kmer_len << 4};
    append_cigar(extension, 1, cigar_m);
    for (unsigned j = 0; j < kmer_len; ++j) {
      extension->dp_score += qseq[qs + j] == ref.c_str()[raw_rs + qs + j] ? SC_MCH : -SC_MIS;
    }

    // right extension
    if (qe < qe0 && re < re0) {
      ksw_extz_t ez_r;
      memset(&ez_r, 0, sizeof(ksw_extz_t));
      const uint8_t *_tseq = reinterpret_cast<const uint8_t *>(ref.c_str() + re);
      const uint8_t *_qseq = reinterpret_cast<const uint8_t *>(qseq + qe);
      ksw_extz2_sse(0, qe0 - qe, _qseq, re0 - re, _tseq, 5, mat, GAPO, GAPE, BANDWIDTH,
                    Z_DROP, END_BONUS, KSW_EZ_EXTZ_ONLY, &ez_r);

      if (ez_r.n_cigar > 0) {
        append_cigar(extension, ez_r.n_cigar, ez_r.cigar);
        extension->dp_score += ez_r.max;
      }
      endclip = ez_r.reach_end ? 0 : qe0 - qe - ez_r.max_q - 1;
      free(ez_r.cigar);
    }

    stringstream cigar_string;

    if (beginclip)
      cigar_string << beginclip << 'S';

    unsigned i = 0;
    int edit_mismatch = 0;
    unsigned ref_pos = region.rs, read_pos = beginclip;
    while (i < extension->n_cigar) {
      int count = extension->cigar[i] >> 4;
      char op = "MID"[extension->cigar[i] & 0xf];
      i++;
      while (i < extension->n_cigar && "MID"[extension->cigar[i] & 0xf] == op) {
        count += extension->cigar[i] >> 4;
        i++;
      }
      cigar_string << count << op;
      switch (op) {
        case 'M':
          for (int j = 0; j < count; j++, ref_pos++, read_pos++) {
            if (ref.c_str()[ref_pos] != qseq[read_pos])
              edit_mismatch++;
          }
          break;
        case 'D':edit_mismatch += count;
          ref_pos += count;
          break;
        case 'I':edit_mismatch += count;
          read_pos += count;
          break;
        default:assert(0);
      }
    }

    if (endclip)
      cigar_string << endclip << 'S';

    r.mapq = get_mapq(r.best, r.secBest);
    a.cigar_string = cigar_string.str();
    a.ref_begin = 0;
    region.score = extension->dp_score;
    a.mismatches = edit_mismatch;
    free(extension);
  }

}

//determine chromo id
int AccAlign::get_tid(Read &R) {
  int tid = R.pos / (offset[1] - offset[0]);
  if (R.pos >= offset[tid] && (tid == (int) name.size() - 1 || R.pos < offset[tid + 1])) {
    return tid;
  } else if (R.pos < offset[tid]) {
    while (tid >= 0) {
      if (R.pos >= offset[tid]) {
        return tid;
      }
      --tid;
    }
  } else {
    while (tid < (int) name.size()) {
      if (R.pos < offset[tid + 1]) {
        return tid;
      }++tid;
    }
  }

  return INT_MAX;
}

void AccAlign::save_region(Read &R, size_t rlen, Region &region,
                           Alignment &a) {
  if (!enable_extension || !region.embed_dist || region.embed_dist == 1) {
    R.pos = region.rs;
    sprintf(R.cigar, "%uM", (unsigned) rlen);
    R.nm = 0;
  } else {
    R.pos = region.rs;
    int cigar_len = a.cigar_string.size();
    strncpy(R.cigar, a.cigar_string.c_str(), cigar_len);
    R.cigar[cigar_len] = '\0';
    R.nm = a.mismatches;
  }
  R.as = region.score;

  R.tid = get_tid(R);

  if (R.tid + 1 < (int) name.size() && R.pos + rlen/2 > offset[R.tid + 1]){
    //reach the end of chromo, switch to next
    R.pos = 1;
    R.tid += 1;
  } else
    R.pos = R.pos - offset[R.tid] + 1;

  //cerr << "Saving region at pos " << R.pos << " as pos " << R.pos -
  //    offset[R.tid] + 1 << " for read " << R.name << endl;
}

void AccAlign::align_read(Read &R) {
  auto start = std::chrono::system_clock::now();

  if (R.strand == '*') {
    return;
  }

  Region region = R.best_region;
  char *s = R.strand == '+' ? R.fwd : R.rev;

  Alignment a;
  size_t rlen = strlen(R.seq);
  score_region(R, s, region, a);
  save_region(R, rlen, region, a);

  auto end = std::chrono::system_clock::now();
  auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
  sw_time += elapsed.count();
}

void AccAlign::wfa_align_read(Read &R) {
  auto start = std::chrono::system_clock::now();

  if (R.strand == '*')
    return;

  Region region = R.best_region;
  size_t rlen = strlen(R.seq);
  char *text = R.strand == '+' ? R.fwd : R.rev;
  const char *pattern = ref.c_str() + region.rs;

  if (enable_extension && region.embed_dist) {
    // Allocate MM
    mm_allocator_t *const mm_allocator = mm_allocator_new(BUFFER_SIZE_8M);
    // Set penalties
    affine_penalties_t affine_penalties = {
        .match = -SC_MCH,
        .mismatch = SC_MIS,
        .gap_opening = GAPO,
        .gap_extension = GAPE,
    };

    // Init Affine-WFA
    affine_wavefronts_t *affine_wavefronts = affine_wavefronts_new_complete(
        rlen, rlen, &affine_penalties, NULL, mm_allocator);
    // Align
    affine_wavefronts_align(affine_wavefronts, pattern, rlen, text, rlen);

    // Display alignment
    edit_cigar_t *edit_cigar = &affine_wavefronts->edit_cigar;
    std::stringstream cigar;
    char last_op = edit_cigar->operations[edit_cigar->begin_offset];

    int last_op_length = 1, nm = 0;
    int i;
    for (i = edit_cigar->begin_offset + 1; i < edit_cigar->end_offset; ++i) {
      if (edit_cigar->operations[i] != 'M')
        nm++;

      //convert to M if X
      last_op = last_op == 'X' ? 'M' : last_op;
      if (edit_cigar->operations[i] == last_op || (edit_cigar->operations[i] == 'X' && last_op == 'M')) {
        ++last_op_length;
      } else {
        cigar << last_op_length << last_op;
        last_op = edit_cigar->operations[i];
        last_op_length = 1;
      }
    }
    cigar << last_op_length << last_op;
    int cigar_len = cigar.str().length();
    strncpy(R.cigar, cigar.str().c_str(), cigar_len);
    R.cigar[cigar_len] = '\0';
    R.nm = nm;

    R.as = edit_cigar_score_gap_affine(edit_cigar, &affine_penalties);

    // Free
    affine_wavefronts_delete(affine_wavefronts);
    mm_allocator_delete(mm_allocator);
  } else {
    sprintf(R.cigar, "%uM", (unsigned) rlen);
    R.as = 0;
    R.nm = 0;
  }

  R.pos = region.rs;
  R.tid = 0;
  for (size_t j = 0; j < name.size(); j++) {
    if (offset[j + 1] > R.pos) {
      R.tid = j;
      break;
    }
  }
  R.pos = R.pos - offset[R.tid] + 1;

  auto end = std::chrono::system_clock::now();
  auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
  sw_time += elapsed.count();
}

bool AccAlign::pairdis_filter(vector<Region> &in_regions1, vector<Region> &in_regions2,
                              bool flag1[], bool flag2[],
                              unsigned &best1, unsigned &next1, unsigned &best2, unsigned &next2) {
  unsigned sum_best = 0, sum_next = 0;
  // init best/next, it's out of the index, could be used to check if the value has been set
  best1 = next1 = in_regions1.size();
  best2 = next2 = in_regions2.size();
  bool has_pair = false;

  for (unsigned i = 0; i < in_regions1.size(); i++) {
    Region tmp;
    tmp.rs = in_regions1[i].rs - pairdis;
    int start = lower_bound(in_regions2.begin(), in_regions2.end(), tmp,
                            [](const Region &left, const Region &right) {
                              return left.rs < right.rs;
                            }
    ) - in_regions2.begin();

    tmp.rs = in_regions1[i].rs + pairdis;
    int end = upper_bound(in_regions2.begin(), in_regions2.end(), tmp,
                          [](const Region &left, const Region &right) {
                            return left.rs < right.rs;
                          }
    ) - in_regions2.begin();

    for (int j = start; j < end; j++) {
      unsigned sum_cov = in_regions1[i].cov + in_regions2[j].cov;
      if (sum_cov > sum_best) {
        sum_next = sum_best;
        next1 = best1;
        next2 = best2;
        sum_best = sum_cov;
        best1 = i;
        best2 = j;
      } else if (sum_cov > sum_next) {
        sum_next = sum_cov;
        next1 = i;
        next2 = j;
      }
      has_pair = true;
      flag1[i] = 1;
      flag2[j] = 1;
    }
  }

  return has_pair;
}

void AccAlign::sam_header(void) {
  ostringstream so;

  so << "@HD\tVN:1.3\tSO:coordinate\n";
  for (size_t i = 0; i < name.size(); i++)
    so << "@SQ\tSN:" << name[i] << '\t' << "LN:" << offset[i + 1] - offset[i] << '\n';
  so << "@PG\tID:AccAlign\tPN:AccAlign\tVN:0.0\n";

  if (sam_name.length())
    sam_stream << so.str();
  else
    cout << so.str();
}

void AccAlign::open_output(string &out_file) {
  sam_name = out_file;

  if (out_file.length()) {
    cerr << "setting output file as " << out_file << endl;
    sam_stream.open(out_file);
  } else {
    assert(g_batch_file.length() == 0);
    setvbuf(stdout, NULL, _IOFBF, 16 * 1024 * 1024);
  }

  sam_header();
}

void AccAlign::close_output() {
  if (sam_name.length()) {
    cerr << "Closing output file " << sam_name << endl;
    sam_stream.close();
  }
}

AccAlign::AccAlign(Reference &r) :
    ref(r.ref), name(r.name),
    offset(r.offset),
    keyv(r.keyv), posv(r.posv) {

  input_io_time = parse_time = 0;
  seeding_time = hit_count_time = 0;
  vpair_build_time = 0;
  sw_time = sam_time = sam_pre_time = sam_out_time = 0;
  vpair_sort_count = 0;

  if (g_embed_file.size())
    embedding = new Embedding(g_embed_file.c_str());
  else
    embedding = new Embedding;
}

AccAlign::~AccAlign() {
  delete embedding;
}

struct tbb_map {
  AccAlign *accalign;

 public:
  tbb_map(AccAlign *obj) : accalign(obj) {}

  Read *operator()(Read *r) {
    accalign->map_read(*r);
    return r;
  }

  ReadPair operator()(ReadPair p) {
    Read *mate1 = std::get<0>(p);
    Read *mate2 = std::get<1>(p);
    accalign->map_paired_read(*mate1, *mate2);

    return p;
  }
};

struct tbb_align {
  AccAlign *accalign;

 public:
  tbb_align(AccAlign *obj) : accalign(obj) {}

  Read *operator()(Read *r) {
    accalign->align_read(*r);
    return r;
  }

  ReadPair operator()(ReadPair p) {
    Read *mate1 = std::get<0>(p);
    Read *mate2 = std::get<1>(p);
    accalign->align_read(*mate1);
    accalign->align_read(*mate2);

    int mapq_pe = mate1->mapq > mate2->mapq ? mate1->mapq : mate2->mapq;
    if (mate1->mapq < mapq_pe) mate1->mapq = (int) (.2f * mate1->mapq + .8f * mapq_pe + .499f);
    if (mate2->mapq < mapq_pe) mate2->mapq = (int) (.2f * mate2->mapq + .8f * mapq_pe + .499f);

    return p;
  }
};

struct tbb_score {
  AccAlign *accalign;

 public:
  tbb_score(AccAlign *obj) : accalign(obj) {}

  continue_msg operator()(Read *r) {
    accalign->print_sam(*r);
    delete r;
    return continue_msg();
  }

  continue_msg operator()(ReadPair p) {
    Read *mate1 = std::get<0>(p);
    Read *mate2 = std::get<1>(p);
    accalign->print_paired_sam(*mate1, *mate2);
    delete mate1;
    delete mate2;

    return continue_msg();
  }
};

bool AccAlign::tbb_fastq(const char *F1, const char *F2) {
  gzFile in1 = gzopen(F1, "rt");
  if (in1 == Z_NULL)
    return false;

  gzFile in2 = Z_NULL;
  char is_paired = false;
  if (strlen(F2) > 0) {
    is_paired = true;
    in2 = gzopen(F2, "rt");
    if (in2 == Z_NULL)
      return false;
  }

  //struct timeval start, end;

  cerr << "Reading fastq file " << F1 << ", " << F2 << "\n";

  // replace this broadcast node with source node
//  if (!is_paired) {
//    graph g;
//
//    source_node < Read * > input_node(g, [&](Read *&r) -> bool {
//      auto start = std::chrono::system_clock::now();
//
//      if (gzeof(in1) || (gzgetc(in1) == EOF))
//        return false;
//
//      r = new Read;
//      in1 >> *r;
//
//      auto end = std::chrono::system_clock::now();
//      auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
//      input_io_time += elapsed.count();
//
//      if (!strlen(r->seq))
//        return false;
//
//      return true;
//
//    }, false);
//
//    int max_objects = 10000000;
//    limiter_node < Read * > lnode(g, max_objects);
//    function_node < Read * , Read * > map_node(g, unlimited, tbb_map(this));
//    function_node < Read * , Read * > align_node(g, unlimited, tbb_align(this));
//    function_node < Read * , continue_msg > score_node(g, 1, tbb_score(this));
//
//    make_edge(score_node, lnode.decrement);
//    make_edge(align_node, score_node);
//    make_edge(map_node, align_node);
//    make_edge(lnode, map_node);
//    make_edge(input_node, map_node);
//    input_node.activate();
//    g.wait_for_all();
//  } else {
  if (is_paired) {
    graph g;
    source_node<ReadPair> input_node(g, [&](ReadPair &rp) -> bool {
      auto start = std::chrono::system_clock::now();

      bool end1 = gzgetc(in1) == EOF;
      bool end2 = gzgetc(in2) == EOF;
      if (gzeof(in1) || end1 || end2)
        return false;

      Read *r = new Read;
      in1 >> *r;
      if (!strlen(r->seq))
        return false;
      get<0>(rp) = r;

      Read *r2 = new Read;
      in2 >> *r2;
      if (!strlen(r2->seq))
        return false;
      get<1>(rp) = r2;

      auto end = std::chrono::system_clock::now();
      auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
      input_io_time += elapsed.count();

      return true;
    }, false);

    int max_objects = 1000000;
    limiter_node<ReadPair> lnode(g, max_objects);
    function_node<ReadPair, ReadPair> map_node(g, unlimited, tbb_map(this));
    function_node<ReadPair, ReadPair> align_node(g, unlimited, tbb_align(this));
    function_node<ReadPair, continue_msg> score_node(g, 1, tbb_score(this));

    make_edge(score_node, lnode.decrement);
    make_edge(align_node, score_node);
    make_edge(map_node, align_node);
    make_edge(lnode, map_node);
    make_edge(input_node, map_node);
    input_node.activate();
    g.wait_for_all();
  }

  gzclose(in1);
  gzclose(in2);

  return true;
}

int main(int ac, char **av) {
  if (ac < 3) {
    print_usage();
    return 0;
  }

  int opn = 1;
  int kmer_temp = 0;
  while (opn < ac) {
    bool flag = false;
    if (av[opn][0] == '-') {
      if (av[opn][1] == 't') {
        g_ncpus = atoi(av[opn + 1]);
        opn += 2;
        flag = true;
      } else if (av[opn][1] == 'l') {
        kmer_temp = atoi(av[opn + 1]);
        opn += 2;
        flag = true;
      } else if (av[opn][1] == 'o') {
        g_out = av[opn + 1];
        opn += 2;
        flag = true;
      } else if (av[opn][1] == 'e') {
        g_embed_file = av[opn + 1];
        opn += 2;
        flag = true;
      } else if (av[opn][1] == 'b') {
        g_batch_file = av[opn + 1];
        opn += 2;
        flag = true;
      } else if (av[opn][1] == 'p') {
        pairdis = atoi(av[opn + 1]);
        opn += 2;
        flag = true;
      } else if (av[opn][1] == 'x') {
        enable_extension = false;
        opn += 1;
        flag = true;
      } else if (av[opn][1] == 'w') {
        enable_wfa_extension = true;
        opn += 1;
        flag = true;
      } else if (av[opn][1] == 'd') {
        extend_all = true;
        opn += 1;
        flag = true;
      } else {
        print_usage();
      }
    }
    if (!flag)
      break;
  }
  if (kmer_temp != 0)
    kmer_len = kmer_temp;
  mask = kmer_len == 32 ? ~0 : (1ULL << (kmer_len * 2)) - 1;

  cerr << "Using " << g_ncpus << " cpus " << endl;
  cerr << "Using kmer length " << kmer_len << " and step size " << kmer_step << endl;

  tbb::task_scheduler_init init(g_ncpus);
  make_code();

  // load reference once
  Reference *r = new Reference(av[opn]);
  opn++;
  if (enable_extension && !enable_wfa_extension)
    ksw_gen_simple_mat(5, mat, SC_MCH, SC_MIS, SC_AMBI);

  size_t total_begin = time(NULL);

  auto start = std::chrono::system_clock::now();

  AccAlign f(*r);
  f.open_output(g_out);

  if (opn == ac - 1) {
    f.fastq(av[opn], "\0", false);
//    f.tbb_fastq(av[opn], "\0");
  } else if (opn == ac - 2) {
//    f.fastq(av[opn], av[opn + 1], false);
    f.tbb_fastq(av[opn], av[opn + 1]);
  } else {
    print_usage();
    return 0;
  }

  auto end = std::chrono::system_clock::now();
  auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
  cerr << "Time to align: " << elapsed.count() / 1000 << " secs\n";

  f.print_stats();
  f.close_output();

  delete r;

  cerr << "Total time: " << (time(NULL) - total_begin) << " secs\n";

  return 0;
}
