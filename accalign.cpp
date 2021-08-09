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
bool toExtend = true, useWFA = false, enable_softclip = false;

int g_ncpus = 1;
float delTime = 0, alignTime = 0, mapqTime, keyvTime = 0, posvTime = 0, sortTime = 0;

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

  int i = 0;
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
  cerr << "\t-s Take into account softclipping. No softclipping by default. \n";
}

void AccAlign::print_stats() {
#if DBGPRINT
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
       "Embedding time(total/actual):\t" << embedding_time / 1000000 << "/" << embedding->embed_time / 1000000 << "\n"
       <<
       "Extending time (+ build output string if ENABLE_GPU):\t" << sw_time / 1000000 << "\n" <<
       "Mapq && mark for extention time:\t" << mapqTime / 1000000 << "\n" <<
       "SAM output time :\t" << sam_time / 1000000 << "\n" <<
       std::endl << endl;

  cerr << "Total pairs sorted: " << vpair_sort_count << endl;
#endif
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
  int rlen = strlen(read.seq);

  cregion.beg = cregion.pos;
  cregion.end = cregion.pos + rlen < ref.size() ? cregion.pos + rlen :
                ref.size();

  char *strand = S == '+' ? read.fwd : read.rev;

  if (cregion.embed_dist)
    rectify_start_pos(strand, cregion, rlen);

  cregion.is_exact = cregion.is_aligned = false;
  read.best_region = cregion;
}

void AccAlign::lsh_filter(char *Q, size_t rlen,
                          vector<Region> &candidate_regions,
                          int &best_threshold, int &next_threshold,
                          unsigned &best_idx, unsigned &next_idx) {
  const char **candidate_refs;
  unsigned ncandidates = candidate_regions.size();

  if (!ncandidates)
    return;

  // alloc space
  candidate_refs = new const char *[ncandidates + 1];
  candidate_refs[0] = Q;

  // add all remaining candidates for embedding
  const char *ptr_ref = ref.c_str();
  for (unsigned i = 0; i < ncandidates; ++i) {
    Region &r = candidate_regions[i];
    candidate_refs[i + 1] = ptr_ref + r.pos;
  }

  // now do embedding
  embedding->embeddata_iterative_update(candidate_regions, candidate_refs, ncandidates + 1,
                                        rlen, best_threshold, next_threshold, true, best_idx, next_idx);

  delete[] candidate_refs;
}

void AccAlign::pigeonhole_query(char *Q, size_t rlen, vector<Region> &candidate_regions,
                                char S, int err_threshold, unsigned &best, unsigned ori_slide) {
  int max_cov = 0;
  unsigned kmer_window = kmer_len + kmer_step - 1;
  unsigned nwindows = (rlen - ori_slide) / kmer_window;
  unsigned nkmers = nwindows * kmer_step;
  size_t ntotal_hits = 0;
  size_t b[nkmers], e[nkmers];
  unsigned kmer_idx = 0;
  unsigned ori_slide_bk = ori_slide;

  // Take non-overlapping seeds and find all hits
  auto start = std::chrono::system_clock::now();
  for (; ori_slide + kmer_window <= rlen; ori_slide += kmer_window) {
    for (int step = 0; step < kmer_step; step++) {
      uint64_t k = 0;
      for (size_t j = ori_slide + step; j < ori_slide + kmer_len + step; j++)
        k = (k << 2) + *(Q + j);
      size_t hash = (k & mask) % MOD;
      b[kmer_idx] = keyv[hash];
      e[kmer_idx] = keyv[hash + 1];
      ntotal_hits += (e[kmer_idx] - b[kmer_idx]);
      kmer_idx++;
    }
  }
  assert(kmer_idx == nkmers);
  auto end = std::chrono::system_clock::now();
  auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
  keyvTime += elapsed.count();

  // if we have no hits, we are done
  if (!ntotal_hits)
    return;

  uint32_t top_pos[nkmers];
  int step_off[nkmers], rel_off[nkmers];
  uint32_t MAX_POS = numeric_limits<uint32_t>::max();

  start = std::chrono::system_clock::now();
  // initialize top values with first values for each kmer.
  for (unsigned i = 0; i < nkmers; i++) {
    if (b[i] < e[i]) {
      top_pos[i] = posv[b[i]];
      step_off[i] = i % kmer_step;
      rel_off[i] = (i / kmer_step) * kmer_window;
      int shift_pos = rel_off[i] + step_off[i] + ori_slide_bk;
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
  uint32_t last_pos = MAX_POS;
  int last_cov = 0;

  start = std::chrono::system_clock::now();
  while (nprocessed < ntotal_hits) {
    //find min
    uint32_t *min_item = min_element(top_pos, top_pos + nkmers);
    uint32_t min_pos = *min_item;
    int min_kmer = min_item - top_pos;

    // kick off prefetch for next round
    __builtin_prefetch(posv + b[min_kmer] + 1);

    // if previous min element was same as current one, increment coverage.
    // otherwise, check if last min element's coverage was high enough to
    // make it a candidate region
    if (min_pos == last_pos) {
      last_cov++;
    } else {
      if (last_cov >= err_threshold) {
        Region r;
        r.cov = last_cov;
        r.pos = last_pos;
        //cerr << "Adding " << last_pos << " with cov " << r.cov <<
        //    " as candidate for dir " << S << endl;
        assert(r.pos != MAX_POS && r.pos < MAX_POS);

        // add high coverage regions up front so that we dont have to
        // sort later.
        // NOTE: Inserting at vector front is crazy expensive. If we
        // simply change > max_cov to >=max_cov, timing will blow up.
        // Just keep it > max_cov as we are anyway interested in only
        // the region with max cov, and we need to embed all regions any
        // way.
        if (last_cov > max_cov) {
          max_cov = last_cov;
          best = candidate_regions.size();
        }
        candidate_regions.push_back(r);
      }
      last_cov = 1;
    }
    last_pos = min_pos;

    // add next element
    b[min_kmer]++;
    uint32_t next_pos = b[min_kmer] < e[min_kmer] ? posv[b[min_kmer]] : MAX_POS;
    if (next_pos != MAX_POS) {
      int shift_pos = rel_off[min_kmer] + step_off[min_kmer] + ori_slide_bk;
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
  if (last_cov >= err_threshold && last_pos != MAX_POS) {
    Region r;
    r.cov = last_cov;
    r.pos = last_pos;
    if (last_cov > max_cov) {
      max_cov = last_cov;
      best = candidate_regions.size();
    }
    candidate_regions.push_back(r);
  }
  end = std::chrono::system_clock::now();
  elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
  hit_count_time += elapsed.count();
}

void AccAlign::pghole_wrapper(Read &R,
                              vector<Region> &fcandidate_regions,
                              vector<Region> &rcandidate_regions,
                              unsigned &fbest,
                              unsigned &rbest,
                              unsigned ori_slide) {
  size_t rlen = strlen(R.seq);
  pigeonhole_query(R.fwd, rlen, fcandidate_regions, '+', 2, fbest, ori_slide);
  pigeonhole_query(R.rev, rlen, rcandidate_regions, '-', 2, rbest, ori_slide);
  unsigned nfregions = fcandidate_regions.size();
  unsigned nrregions = rcandidate_regions.size();

  // if we don't find any regions, try lowering the error threshold
  if (!nfregions && !nrregions) {
    pigeonhole_query(R.fwd, rlen, fcandidate_regions, '+', 1, fbest, ori_slide);
    pigeonhole_query(R.rev, rlen, rcandidate_regions, '-', 1, rbest, ori_slide);
  }
}

void AccAlign::embed_wrapper_pair(Read &R1, Read &R2,
                                  vector<Region> &candidate_regions_f1, vector<Region> &candidate_regions_r1,
                                  vector<Region> &candidate_regions_f2, vector<Region> &candidate_regions_r2,
                                  bool flag_f1[], bool flag_r1[], bool flag_f2[], bool flag_r2[],
                                  unsigned &best_f1, unsigned &next_f1, unsigned &best_r1, unsigned &next_r1,
                                  unsigned &best_f2, unsigned &next_f2, unsigned &best_r2, unsigned &next_r2) {

  unsigned nregions_f1 = candidate_regions_f1.size();
  unsigned nregions_r1 = candidate_regions_r1.size();
  unsigned nregions_f2 = candidate_regions_f2.size();
  unsigned nregions_r2 = candidate_regions_r2.size();

  if (best_f1 != nregions_f1 && next_f1 != nregions_f1 && best_r2 != nregions_r2 && next_r2 != nregions_r2)
    assert(candidate_regions_f1[best_f1].cov + candidate_regions_r2[best_r2].cov >=
        candidate_regions_f1[next_f1].cov + candidate_regions_r2[next_r2].cov);
  if (best_r1 != nregions_r1 && next_r1 != nregions_r1 && best_f2 != nregions_f2 && next_f2 != nregions_f2)
    assert(candidate_regions_r1[best_r1].cov + candidate_regions_f2[best_f2].cov >=
        candidate_regions_r1[next_r1].cov + candidate_regions_f2[next_f2].cov);

  //candidate_refs: the first one is the read, then the all the candidates' coresponding reference
  const char **candidate_refs_f1, **candidate_refs_r1, **candidate_refs_f2, **candidate_refs_r2;
  const char *ptr_ref = ref.c_str();

  candidate_refs_f1 = new const char *[nregions_f1 + 1];
  candidate_refs_f1[0] = R1.fwd;
  for (unsigned i = 0; i < nregions_f1; ++i) {
    candidate_refs_f1[i + 1] = ptr_ref + candidate_regions_f1[i].pos;
  }

  candidate_refs_r1 = new const char *[nregions_r1 + 1];
  candidate_refs_r1[0] = R1.rev;
  for (unsigned i = 0; i < nregions_r1; ++i) {
    candidate_refs_r1[i + 1] = ptr_ref + candidate_regions_r1[i].pos;
  }

  candidate_refs_f2 = new const char *[nregions_f2 + 1];
  candidate_refs_f2[0] = R2.fwd;
  for (unsigned i = 0; i < nregions_f2; ++i) {
    candidate_refs_f2[i + 1] = ptr_ref + candidate_regions_f2[i].pos;
  }

  candidate_refs_r2 = new const char *[nregions_r2 + 1];
  candidate_refs_r2[0] = R2.rev;
  for (unsigned i = 0; i < nregions_r2; ++i) {
    candidate_refs_r2[i + 1] = ptr_ref + candidate_regions_r2[i].pos;
  }

  //embedQ
  size_t rlen = strlen(R1.seq);
  int elen = rlen * embedding->efactor;
  char embeddedQ_f1[elen], embeddedQ_r1[elen], embeddedQ_f2[elen], embeddedQ_r2[elen];
  embedding->embedstr(candidate_refs_f1, rlen, elen, 0, 0, embeddedQ_f1);
  embedding->embedstr(candidate_refs_r1, rlen, elen, 0, 0, embeddedQ_r1);
  embedding->embedstr(candidate_refs_f2, rlen, elen, 0, 0, embeddedQ_f2);
  embedding->embedstr(candidate_refs_r2, rlen, elen, 0, 0, embeddedQ_r2);

  //embed best and next pairs
  int thresholds[4] = {elen, elen, elen, elen};
  embedding->embed_two_pairs(candidate_regions_f1, candidate_regions_r2, rlen, thresholds,
                             embeddedQ_f1, embeddedQ_r2, candidate_refs_f1, candidate_refs_r2,
                             nregions_f1, nregions_r2, flag_f1, flag_r2,
                             best_f1, next_f1, best_r2, next_r2);
  embedding->embed_two_pairs(candidate_regions_r1, candidate_regions_f2, rlen, thresholds + 2,
                             embeddedQ_r1, embeddedQ_f2, candidate_refs_r1, candidate_refs_f2,
                             nregions_r1, nregions_f2, flag_r1, flag_f2,
                             best_r1, next_r1, best_f2, next_f2);

  // threshold is the second minimal embed distance (because mapq needs second best candidate's embed dist)
  sort(thresholds, thresholds + 4);
  int threshold = thresholds[1];

  // embed rest candidates
  embedding->embeddata_pair(candidate_regions_f1, embeddedQ_f1, candidate_refs_f1, nregions_f1,
                            flag_f1, rlen, threshold);
  embedding->embeddata_pair(candidate_regions_r1, embeddedQ_r1, candidate_refs_r1, nregions_r1,
                            flag_r1, rlen, threshold);
  embedding->embeddata_pair(candidate_regions_f2, embeddedQ_f2, candidate_refs_f2, nregions_f2,
                            flag_f2, rlen, threshold);
  embedding->embeddata_pair(candidate_regions_r2, embeddedQ_r2, candidate_refs_r2, nregions_r2,
                            flag_r2, rlen, threshold);

  delete[] candidate_refs_f1;
  delete[] candidate_refs_r1;
  delete[] candidate_refs_f2;
  delete[] candidate_refs_r2;
}

void AccAlign::embed_wrapper(Read &R, bool ispe,
                             vector<Region> &fcandidate_regions, vector<Region> &rcandidate_regions,
                             unsigned &fbest, unsigned &fnext, unsigned &rbest, unsigned &rnext) {
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

  // pass the one with the highest coverage in first
  bool fwd_first = true;
  vpair_sort_count += nfregions + nrregions;
  if (nfregions > 1 && nrregions > 1) {
    if (fcandidate_regions[0].cov < rcandidate_regions[0].cov) {
      fwd_first = false;
    }
  }

  // embed now, but only in the case where we have > 1 regions either for the
  // forward or for reverse or for both strands. If we have only 1 region
  // globally, there is no point embedding.
  size_t rlen = strlen(R.seq);
  int best_threshold = rlen * embedding->efactor;
  int next_threshold = rlen * embedding->efactor;

  if ((nfregions == 0) || (nrregions == 0)) {
    fbest = fnext = rbest = rnext = 0;
    if (nfregions) {
      lsh_filter(R.fwd, rlen, fcandidate_regions, best_threshold, next_threshold, fbest, fnext);
    }
    if (nrregions) {
      lsh_filter(R.rev, rlen, rcandidate_regions, best_threshold, next_threshold, rbest, rnext);
    }
  } else {
    if (fwd_first) {
      lsh_filter(R.fwd, rlen, fcandidate_regions, best_threshold, next_threshold, fbest, fnext);
      lsh_filter(R.rev, rlen, rcandidate_regions, best_threshold, next_threshold, rbest, rnext);
    } else {
      lsh_filter(R.rev, rlen, rcandidate_regions, best_threshold, next_threshold, rbest, rnext);
      lsh_filter(R.fwd, rlen, fcandidate_regions, best_threshold, next_threshold, fbest, fnext);
    }
  }

}

int AccAlign::get_mapq(int best, int secbest, bool hasSecbest, int rlen) {
  best = MIS_PENALTY * best;
  secbest = MIS_PENALTY * secbest;

  int64_t scMin = -0.6 + -0.6 * rlen;
  int64_t diff = abs(scMin);
  int64_t bestOver = best - scMin;
  int ret;

  if (!hasSecbest) {
    if (bestOver >= diff * (double) 0.9f)
      ret = 42;
    else if (bestOver >= diff * (double) 0.8f)
      ret = 40;
    else if (bestOver >= diff * (double) 0.7f)
      ret = 24;
    else if (bestOver >= diff * (double) 0.6f)
      ret = 23;
    else if (bestOver >= diff * (double) 0.5f)
      ret = 8;
    else if (bestOver >= diff * (double) 0.4f)
      ret = 3;
    else
      ret = 0;
  } else {
    int64_t bestdiff = abs(abs(static_cast<long>(best)) - abs(static_cast<long>(secbest)));
    if (bestdiff >= diff * (double) 0.9f) {
      if (bestOver == diff) {
        ret = 39;
      } else {
        ret = 33;
      }
    } else if (bestdiff >= diff * (double) 0.8f) {
      if (bestOver == diff) {
        ret = 38;
      } else {
        ret = 27;
      }
    } else if (bestdiff >= diff * (double) 0.7f) {
      if (bestOver == diff) {
        ret = 37;
      } else {
        ret = 26;
      }
    } else if (bestdiff >= diff * (double) 0.6f) {
      if (bestOver == diff) {
        ret = 36;
      } else {
        ret = 22;
      }
    } else if (bestdiff >= diff * (double) 0.5f) {
      if (bestOver == diff) {
        ret = 35;
      } else if (bestOver >= diff * (double) 0.84f) {
        ret = 25;
      } else if (bestOver >= diff * (double) 0.68f) {
        ret = 16;
      } else {
        ret = 5;
      }
    } else if (bestdiff >= diff * (double) 0.4f) {
      if (bestOver == diff) {
        ret = 34;
      } else if (bestOver >= diff * (double) 0.84f) {
        ret = 21;
      } else if (bestOver >= diff * (double) 0.68f) {
        ret = 14;
      } else {
        ret = 4;
      }
    } else if (bestdiff >= diff * (double) 0.3f) {
      if (bestOver == diff) {
        ret = 32;
      } else if (bestOver >= diff * (double) 0.88f) {
        ret = 18;
      } else if (bestOver >= diff * (double) 0.67f) {
        ret = 15;
      } else {
        ret = 3;
      }
    } else if (bestdiff >= diff * (double) 0.2f) {
      if (bestOver == diff) {
        ret = 31;
      } else if (bestOver >= diff * (double) 0.88f) {
        ret = 17;
      } else if (bestOver >= diff * (double) 0.67f) {
        ret = 11;
      } else {
        ret = 0;
      }
    } else if (bestdiff >= diff * (double) 0.1f) {
      if (bestOver == diff) {
        ret = 30;
      } else if (bestOver >= diff * (double) 0.88f) {
        ret = 12;
      } else if (bestOver >= diff * (double) 0.67f) {
        ret = 7;
      } else {
        ret = 0;
      }
    } else if (bestdiff > 0) {
      if (bestOver >= diff * (double) 0.67f) {
        ret = 6;
      } else {
        ret = 2;
      }
    } else {
      assert(bestdiff == 0);
      if (bestOver >= diff * (double) 0.67f) {
        ret = 1;
      } else {
        ret = 0;
      }
    }
  }
  return ret;
}

void AccAlign::set_mapq(Read &R, vector<Region> &fcandidate_regions,
                        vector<Region> &rcandidate_regions, int fbest, int fnext,
                        int rbest, int rnext) {
  int len = embedding->efactor * strlen(R.seq);

  unsigned nfregions = fcandidate_regions.size();
  unsigned nrregions = rcandidate_regions.size();
  if (nfregions == 0) {
    if (nrregions == 0)
      return;
    int best_score = rcandidate_regions[rbest].embed_dist;
    if (nrregions > 1) {
      int next_score = rcandidate_regions[rnext].embed_dist;
      R.mapq = get_mapq(best_score, next_score, true, len);
    } else {
      //only 1 alignment
      R.mapq = get_mapq(best_score, 0, false, len);
    }
  } else if (nrregions == 0) {
    int best_score = fcandidate_regions[fbest].embed_dist;
    if (nfregions > 1) {
      int next_score = fcandidate_regions[fnext].embed_dist;
      R.mapq = get_mapq(best_score, next_score, true, len);
    } else {
      // only 1 alignment
      R.mapq = get_mapq(best_score, 0, false, len);
    }
  } else {
    // we have atleast 1 forward and 1 reverse candidate
    if (fcandidate_regions[fbest].embed_dist <
        rcandidate_regions[rbest].embed_dist) {
      int second_best = rcandidate_regions[rbest].embed_dist;
      if (nfregions > 1 && fcandidate_regions[fnext].embed_dist < second_best)
        second_best = fcandidate_regions[fnext].embed_dist;
      R.mapq = get_mapq(fcandidate_regions[fbest].embed_dist, second_best, true, len);
    } else {
      int second_best = fcandidate_regions[fbest].embed_dist;
      if (nrregions > 1 && rcandidate_regions[rnext].embed_dist < second_best)
        second_best = rcandidate_regions[rnext].embed_dist;
      R.mapq = get_mapq(rcandidate_regions[rbest].embed_dist, second_best, true, len);
    }
  }
}

void AccAlign::map_read(Read &R) {

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
  pghole_wrapper(R, fcandidate_regions, rcandidate_regions, fbest, rbest, 0);
  unsigned nfregions = fcandidate_regions.size();
  unsigned nrregions = rcandidate_regions.size();
  end = std::chrono::system_clock::now();
  elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
  seeding_time += elapsed.count();

  if (nfregions == 0 && nrregions == 0) {
    R.strand = '*';
    return;
  }

  start = std::chrono::system_clock::now();
  // before doing embedding, move highest cov region to front
  if (nfregions > 1 && fbest != 0) {
    iter_swap(fcandidate_regions.begin() + fbest,
              fcandidate_regions.begin());
  }
  if (nrregions > 1 && rbest != 0) {
    iter_swap(rcandidate_regions.begin() + rbest,
              rcandidate_regions.begin());
  }
  end = std::chrono::system_clock::now();
  elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
  swap_time += elapsed.count();
  seeding_time += elapsed.count();

  unsigned fnext, rnext;
  start = std::chrono::system_clock::now();
  embed_wrapper(R, false, fcandidate_regions, rcandidate_regions,
                fbest, fnext, rbest, rnext);
  end = std::chrono::system_clock::now();
  elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
  embedding_time += elapsed.count();

  start = std::chrono::system_clock::now();
  set_mapq(R, fcandidate_regions, rcandidate_regions, fbest, fnext,
           rbest, rnext);
  if (nfregions == 0) {
    if (nrregions == 0)
      return;
    mark_for_extension(R, '-', rcandidate_regions[rbest]);
    R.strand = '-';
  } else if (nrregions == 0) {
    mark_for_extension(R, '+', fcandidate_regions[fbest]);
    R.strand = '+';
  } else {
    // pick the candidate with smallest embed dist
    // if fwd/rev have same embed_dist, take the hcov one
    // if hcov one not the min dist, take the one with smaller pos (to be consistent with gpu)
    if (fcandidate_regions[fbest].embed_dist < rcandidate_regions[rbest].embed_dist) {
      mark_for_extension(R, '+', fcandidate_regions[fbest]);
      R.strand = '+';
    } else if (fcandidate_regions[fbest].embed_dist > rcandidate_regions[rbest].embed_dist) {
      mark_for_extension(R, '-', rcandidate_regions[rbest]);
      R.strand = '-';
    } else {
      if (fcandidate_regions[fbest].pos < rcandidate_regions[rbest].pos) {
        mark_for_extension(R, '+', fcandidate_regions[fbest]);
        R.strand = '+';
      } else {
        mark_for_extension(R, '-', rcandidate_regions[rbest]);
        R.strand = '-';
      }
    }
  }

  end = std::chrono::system_clock::now();
  elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
  mapqTime += elapsed.count();
}

void AccAlign::map_paired_read(Read &mate1, Read &mate2) {

  auto start = std::chrono::system_clock::now();
  parse(mate1.seq, mate1.fwd, mate1.rev, mate1.rev_str);
  parse(mate2.seq, mate2.fwd, mate2.rev, mate2.rev_str);
  auto end = std::chrono::system_clock::now();
  auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
  parse_time += elapsed.count();

  vector<Region> fcandidate_regions1, rcandidate_regions1;
  vector<Region> fcandidate_regions2, rcandidate_regions2;
  unsigned fbest1 = 0, rbest1 = 0, fbest2 = 0, rbest2 = 0;
  unsigned fnext1 = 0, rnext1 = 0, fnext2 = 0, rnext2 = 0;
  unsigned nfregions1, nrregions1, nfregions2, nrregions2;
  bool *flag_f1, *flag_r1, *flag_f2, *flag_r2;

  bool has_f1r2 = false, has_r1f2 = false;

  // lookup candidates
  for (unsigned i = 0; i < kmer_len; i++) {
    if (has_f1r2 || has_r1f2)
      break;

    start = std::chrono::system_clock::now();
    pghole_wrapper(mate1, fcandidate_regions1, rcandidate_regions1,
                   fbest1, rbest1, i);
    nfregions1 = fcandidate_regions1.size();
    nrregions1 = rcandidate_regions1.size();
    pghole_wrapper(mate2, fcandidate_regions2, rcandidate_regions2,
                   fbest2, rbest2, i);
    nfregions2 = fcandidate_regions2.size();
    nrregions2 = rcandidate_regions2.size();
    end = std::chrono::system_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    seeding_time += elapsed.count();

    //init the next, it's out of the index, could be used to check if the value has been set
    fnext1 = nfregions1, rnext1 = nrregions1, fnext2 = nfregions2, rnext2 = nrregions2;
    // reset best, because best/next depends on sum of cov, not only itself, will be set in pairdist_filter
    fbest1 = nfregions1, rbest1 = nrregions1, fbest2 = nfregions2, rbest2 = nrregions2;

    // filter based on pairdis
    // Note that after filtering, regions wont be sorted based on coverage. They
    // will instead be sorted based on position. pairdis filter needs this
    // ordering so that it can adopt some shortcuts.
    start = std::chrono::system_clock::now();

    flag_f1 = new bool[nfregions1]();
    flag_r1 = new bool[nrregions1]();
    flag_f2 = new bool[nfregions2]();
    flag_r2 = new bool[nrregions2]();
    pairdis_filter(fcandidate_regions1, rcandidate_regions2, flag_f1, flag_r2,
                   fbest1, fnext1, rbest2, rnext2);
    pairdis_filter(rcandidate_regions1, fcandidate_regions2, flag_r1, flag_f2,
                   rbest1, rnext1, fbest2, fnext2);
    end = std::chrono::system_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    vpair_build_time += elapsed.count();

    for (unsigned i = 0; i < nfregions1; i++) {
      if (flag_f1[i]) {
        has_f1r2 = true;
        break;
      }
    }

    if (!has_f1r2) {
      for (unsigned i = 0; i < nrregions1; i++) {
        if (flag_r1[i]) {
          has_r1f2 = true;
          break;
        }
      }
    }

    if (!has_f1r2 && !has_r1f2) {
      fcandidate_regions1.clear();
      fcandidate_regions2.clear();
      rcandidate_regions1.clear();
      rcandidate_regions2.clear();
      delete[] flag_f1;
      delete[] flag_r1;
      delete[] flag_f2;
      delete[] flag_r2;
    }
  }

  if (!has_f1r2 && !has_r1f2) {
    mate1.strand = '*';
    mate2.strand = '*';
    return;
  }

  // now apply embedding filter on filtered regions.
  // But before, rearrange so that regions with high coverage are at the top
  //no need to swap, just embed the best and next first..
  start = std::chrono::system_clock::now();
  embed_wrapper_pair(mate1, mate2,
                     fcandidate_regions1, rcandidate_regions1,
                     fcandidate_regions2, rcandidate_regions2,
                     flag_f1, flag_r1, flag_f2, flag_r2,
                     fbest1, fnext1, rbest1, rnext1,
                     fbest2, fnext2, rbest2, rnext2);
  delete[] flag_f1;
  delete[] flag_r1;
  delete[] flag_f2;
  delete[] flag_r2;
  end = std::chrono::system_clock::now();
  elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
  embedding_time += elapsed.count();

  // finally, pick regions for further extension. here, we make use of pairing
  // again by first checking overall, across mate1 and mate2, who has lowest
  // embed dist. then, we use the fwd or rev strand from that mate as the
  // deciding factor. if we chose fwd for that mate, we pick rev for other
  // mate, and vice versa.
  start = std::chrono::system_clock::now();

  int min_dist_f1r2 = INT_MAX, min_dist = INT_MAX, secmin_dist = INT_MAX;
  Region best_fregion, best_rregion;

  for (unsigned i = 0; i < nfregions1; i++) {
    Region tmp;
    tmp.pos = fcandidate_regions1[i].pos - pairdis;
    auto start = lower_bound(rcandidate_regions2.begin(), rcandidate_regions2.end(), tmp,
                             [](const Region &left, const Region &right) {
                               return left.pos < right.pos;
                             }
    );

    tmp.pos = fcandidate_regions1[i].pos + pairdis;
    auto end = upper_bound(rcandidate_regions2.begin(), rcandidate_regions2.end(), tmp,
                           [](const Region &left, const Region &right) {
                             return left.pos < right.pos;
                           }
    );

    for (auto itr = start; itr != end; ++itr) {
      int sum_embed_dist = fcandidate_regions1[i].embed_dist + itr->embed_dist;
      if (sum_embed_dist < min_dist) {
        secmin_dist = min_dist;
        min_dist = sum_embed_dist;
        best_fregion = fcandidate_regions1[i];
        best_rregion = *itr;
      } else if (sum_embed_dist < secmin_dist) {
        secmin_dist = sum_embed_dist;
      }
    }
  }
  min_dist_f1r2 = min_dist;

  for (unsigned i = 0; i < nrregions1; i++) {
    Region tmp;
    tmp.pos = rcandidate_regions1[i].pos - pairdis;
    auto start = lower_bound(fcandidate_regions2.begin(), fcandidate_regions2.end(), tmp,
                             [](const Region &left, const Region &right) {
                               return left.pos < right.pos;
                             }
    );

    tmp.pos = rcandidate_regions1[i].pos + pairdis;
    auto end = upper_bound(fcandidate_regions2.begin(), fcandidate_regions2.end(), tmp,
                           [](const Region &left, const Region &right) {
                             return left.pos < right.pos;
                           }
    );

    for (auto itr = start; itr != end; ++itr) {
      int sum_embed_dist = rcandidate_regions1[i].embed_dist + itr->embed_dist;
      if (sum_embed_dist < min_dist) {
        secmin_dist = min_dist;
        min_dist = sum_embed_dist;
        best_fregion = *itr;
        best_rregion = rcandidate_regions1[i];
      } else if (sum_embed_dist < secmin_dist) {
        secmin_dist = sum_embed_dist;
      }
    }
  }


  // if there is no candidates, the strand will remain *
  if (min_dist < INT_MAX) {
    if (min_dist_f1r2 == min_dist) {
      mark_for_extension(mate1, '+', best_fregion);
      mark_for_extension(mate2, '-', best_rregion);
      mate1.strand = '+';
      mate2.strand = '-';
    } else {
      mark_for_extension(mate2, '+', best_fregion);
      mark_for_extension(mate1, '-', best_rregion);
      mate1.strand = '-';
      mate2.strand = '+';
    }

    int mapq;
    size_t max_hamming = strlen(mate1.seq) * embedding->efactor * 2; //*2 because it's the sum of embed dist
    if (secmin_dist < INT_MAX)
      mapq = get_mapq(min_dist, secmin_dist, true, max_hamming);
    else
      mapq = get_mapq(min_dist, 0, false, max_hamming);
    mate1.mapq = mapq;
    mate2.mapq = mapq;
  }

  end = std::chrono::system_clock::now();
  elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
  mapqTime += elapsed.count();
}

void AccAlign::snprintf_pair_sam(Read &R, string *s, Read &R2, string *s2) {
  auto start = std::chrono::system_clock::now();

  // 60 is the approximate length for all int
  int size;
  if (!toExtend) {
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
    snprintf(buf,
             size,
             format.c_str(),
             nn.c_str(),
             flag,
             name[R.tid].c_str(),
             R.pos,
             (int) R.mapq,
             R.cigar,
             name[R.tid] == name[R2.tid] ? "=" : name[R2.tid].c_str(),
             R2.pos,
             isize,
             R.seq,
             R.qua,
             R.nm,
             R.as);
    *s = buf;
  } else if (R.strand == '-' && R2.strand != '*') {
    size += strlen(R.name) + name[R.tid].length() + name[R2.tid].length() + 2 * strlen(R.seq);
    char buf[size];
    std::reverse(R.qua, R.qua + strlen(R.qua));
    snprintf(buf,
             size,
             format.c_str(),
             nn.c_str(),
             flag,
             name[R.tid].c_str(),
             R.pos,
             (int) R.mapq,
             R.cigar,
             name[R.tid] == name[R2.tid] ? "=" : name[R2.tid].c_str(),
             R2.pos,
             isize,
             R.rev_str,
             R.qua,
             R.nm,
             R.as);
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
    snprintf(buf,
             size,
             format.c_str(),
             nn2.c_str(),
             flag,
             name[R2.tid].c_str(),
             R2.pos,
             (int) R2.mapq,
             R2.cigar,
             name[R.tid] == name[R2.tid] ? "=" : name[R.tid].c_str(),
             R.pos,
             -isize,
             R2.seq,
             R2.qua,
             R2.nm,
             R2.as);
    *s2 = buf;
  } else if (R2.strand == '-' && R.strand != '*') {
    size += strlen(R2.name) + name[R2.tid].length() + name[R.tid].length() + 2 * strlen(R2.seq);
    char buf[size];
    std::reverse(R2.qua, R2.qua + strlen(R2.qua));
    snprintf(buf,
             size,
             format.c_str(),
             nn2.c_str(),
             flag,
             name[R2.tid].c_str(),
             R2.pos,
             (int) R2.mapq,
             R2.cigar,
             name[R.tid] == name[R2.tid] ? "=" : name[R.tid].c_str(),
             R.pos,
             -isize,
             R2.rev_str,
             R2.qua,
             R2.nm,
             R2.as);
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

void AccAlign::snprintf_sam(Read &R, string *s) {
  auto start = std::chrono::system_clock::now();

  // 50 is the approximate length for all int
  int size;
  if (!toExtend) {
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
        if (useWFA)
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
        if (useWFA)
          acc_obj->wfa_align_read(*(all_reads + i));
        else
          acc_obj->align_read(*(all_reads + i));
      }

      if ((all_reads2 + i)->strand != '*') {
        if (useWFA)
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
    auto start = std::chrono::system_clock::now();
    tbb::task_scheduler_init init(g_ncpus);
    tbb::parallel_for(tbb::blocked_range<size_t>(soff, eoff), Tbb_aligner(ptlread, sams, this));
    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    alignTime += elapsed.count();

    start = std::chrono::system_clock::now();
    for (int i = soff; i < eoff; i++) {
      out_sam(sams + i);
    }
    end = std::chrono::system_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    sam_time += elapsed.count();

    dataQ->push(make_tuple(ptlread, (Read *) NULL));
  } else {
    string sams[2 * eoff];
    auto start = std::chrono::system_clock::now();
    tbb::task_scheduler_init init(g_ncpus);
    tbb::parallel_for(tbb::blocked_range<size_t>(soff, eoff), Tbb_aligner_paired(ptlread, ptlread2, sams, this));

    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    alignTime += elapsed.count();

    start = std::chrono::system_clock::now();
    for (int i = soff; i < 2 * eoff; i++) {
      out_sam(sams + i);
    }
    end = std::chrono::system_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    sam_time += elapsed.count();

    dataQ->push(make_tuple(ptlread, ptlread2));
  }
}

void ksw_align(const char *tseq, int tlen, const char *qseq, int qlen,
               int sc_mch, int sc_mis, int gapo, int gape, ksw_extz_t &ez) {
  int8_t a = sc_mch, b = sc_mis < 0 ? sc_mis : -sc_mis; // a>0 and b<0
  int8_t mat[25] = {a, b, b, b, 0, b, a, b, b, 0, b, b, a, b, 0, b, b, b, a, 0, 0, 0, 0, 0, 0};
  const uint8_t *ts = reinterpret_cast<const uint8_t *>(tseq);
  const uint8_t *qs = reinterpret_cast<const uint8_t *>(qseq);
  memset(&ez, 0, sizeof(ksw_extz_t));
  ksw_extz2_sse(0, qlen, qs, tlen, ts, 5, mat, gapo, gape, -1, -1, 0, 0, &ez);
}

void AccAlign::soft_clip(char *cigar, int &len, Read &R, Region &region) {
  int MAX_MATCH = ceil(0.2 * strlen(R.seq));
  const char *ptr_ref = ref.c_str() + region.beg;

  /*
   * soft clip in the start of read: cigar e.g. 9I18D????,
   * -> 'I' means this subseq only in the read, as it's in the beginning of the read could be regarded as soft clip
   * -> 'D' means this subseq only in the ref, if this follows 'I', should be removed and shift the start pos accordingly
   */
  assert(len >= 2);
  int i = 0;
  int first_val = 0;
  while (i < len && isdigit(cigar[i])) {
    first_val = first_val * 10 + (cigar[i] - '0');
    i++;
  }
  char first_char = cigar[i];
  i++;

  /*
   * limit the match length: e.g. 62M13I13D1M should be soft clip at the end not in the beginning
   */
  if (first_char != 'M' || first_val < MAX_MATCH) {
    int second_val = 0;
    while (i < len && isdigit(cigar[i])) {
      second_val = second_val * 10 + (cigar[i] - '0');
      i++;
    }
    char second_char = cigar[i];
    i++;
    int second_index = i;

    int third_val = 0;
    while (i < len && isdigit(cigar[i])) {
      third_val = third_val * 10 + (cigar[i] - '0');
      i++;
    }
    char third_char = cigar[i];
    i++;

    if ((first_char == 'M' && second_char == 'I' && third_char == 'D') || (first_char == 'I' && second_char == 'D')) {
      int m_val = first_char == 'M' ? first_val : 0;
      int i_val = first_char == 'I' ? first_val : second_val;
      int d_val = second_char == 'D' ? second_val : third_val;

      R.pos += d_val;

      string s = to_string(m_val + i_val);
      size_t s_len = s.length();
      strncpy(cigar, s.c_str(), s_len);
      cigar[s_len] = 'S';

      if (first_char == 'M')
        memcpy(cigar + s_len + 1, cigar + i, len + 1 - i);
      else
        memcpy(cigar + s_len + 1, cigar + second_index, len + 1 - second_index);

      /* justify the AS */
      region.score += GAPO + GAPE * (i_val - 1) + GAPO + GAPE * (d_val - 1);
      for (int j = 0; j < m_val; j++) {
        if (R.seq[j] == ptr_ref[j]) {
          region.score -= SC_MCH;
        } else
          region.score += SC_MIS;
      }

    }
    len = strlen(cigar);
  }


  /*
 * soft clip in the end of read: cigar e.g. ???9I18D1M,
 */
  assert(len >= 2);
  i = len - 1;
  int num_m = 0;
  if (cigar[i] == 'M') {
    i--;
    while (i >= 0 && isdigit(cigar[i])) {
      num_m += (cigar[i] - '0') * pow(10, len - 2 - i);
      i--;
    }
  }

  if (cigar[i] == 'D' && num_m < MAX_MATCH) {
    i--;
    int num_d = 0, d_bk = i;
    while (i >= 0 && isdigit(cigar[i])) {
      num_d += (cigar[i] - '0') * pow(10, d_bk - i);
      i--;
    }

    if (cigar[i] == 'I') {
      i--;
      int num_i = 0, i_bk = i;
      while (i >= 0 && isdigit(cigar[i])) {
        num_i += (cigar[i] - '0') * pow(10, i_bk - i);
        i--;
      }

      std::string s = std::to_string(num_i + num_m);
      char const *pchar = s.c_str();

      strncpy(cigar + i + 1, pchar, strlen(pchar));
      cigar[i + strlen(pchar) + 1] = 'S';
      cigar[i + strlen(pchar) + 2] = '\0';

      /* justify the AS */
      region.score += GAPO + GAPE * (num_i - 1) + GAPO + GAPE * (num_d - 1);
      for (int j = 0; j < num_m; j++) {
        if (R.seq[len - 1 - j] == ptr_ref[len - 1 - j]) {
          region.score -= SC_MCH;
        } else
          region.score += SC_MIS;
      }

    }
  }
  len = strlen(cigar);
}

//wield to find some cigar SRR098401.31169: 52M24I24D, so need to take care for the case ..?I?D and ..?D?I
void AccAlign::rectify_cigar(char *cigar, int len, Read &R, Region &region) {
  if (enable_softclip)
    soft_clip(cigar, len, R, region);

  // D at the end => the read has been total aligned, but ref not
  // remove the D and number of D
  if (cigar[len - 1] == 'D') {
    int end = len - 2;
    int num_d = 0;
    while (isdigit(cigar[end])) {
      num_d += (cigar[end] - '0') * pow(10, len - 2 - end);
      end--;
    }
    region.score += GAPO + GAPE * (num_d - 1);

    cigar[end + 1] = '\0';
  }


  // I at the end => the ref has been total aligned, but read not
  // change the I to M, aggregate the count
  if (cigar[len - 1] == 'I') {
    int num_I = 0;
    int i = len - 2;
    char op_before = 0;
    for (; i >= 0; i--) {
      if (!isdigit(cigar[i])) {
        op_before = cigar[i];
        break;
      }

      num_I += (cigar[i] - '0') * pow(10, len - 2 - i);
    }

    if (op_before == 'M') {
      int num_M = 0;
      int index_M = i;
      i--;
      for (; i >= 0; i--) {
        if (cigar[i] == 'I' || cigar[i] == 'D' || cigar[i] == 'S')
          break;
        num_M += (cigar[i] - '0') * pow(10, index_M - 1 - i);
      }
      num_I += num_M;
    }

    std::string s = std::to_string(num_I);
    char const *pchar = s.c_str();

    strncpy(cigar + i + 1, pchar, strlen(pchar));
    cigar[i + strlen(pchar) + 1] = 'M';
    cigar[i + strlen(pchar) + 2] = '\0';

    /* justify the AS */
    region.score += GAPO + GAPE * (num_I - 1);
  }
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
  int threshold = embedding->cgk2_embed_nmismatch(ptr_ref + region.beg, kmer_len, elen, 0, embeddedQ);
  int shift = 0;

  float indel_len = ceil(MAX_INDEL * rlen / float(100));
  for (int i = -indel_len; i < indel_len; ++i) {
    if (i == 0 || region.beg + i < 0 || region.beg + i >= ref.size())
      continue;

    int nmismatch = embedding->cgk2_embed_nmismatch(ptr_ref + region.beg + i, kmer_len, threshold, 0, embeddedQ);
    if (nmismatch < threshold) {
      threshold = nmismatch;
      shift = i;
    }
  }

  region.beg = region.beg + shift;
}

void AccAlign::score_region(Read &r, char *strand, Region &region,
                            Alignment &a) {
  unsigned len = strlen(r.seq);


  // if the region has a embed distance of 0, then its an exact match
  if (!region.embed_dist || !toExtend) {
    // we found an exact match
    region.is_exact = true;

    // XXX: the scoring here of setting it to len is based on the
    // assumption that our current ssw impl. gives a best score of 150
    region.score = len;
  } else {
    region.is_exact = false;
    const char *ptr_ref = ref.c_str() + region.beg;
    const char *ptr_read = strand;

    ksw_extz_t ez;
    ksw_align(ptr_ref, len, ptr_read, len, SC_MCH, SC_MIS, GAPO, GAPE, ez);

    stringstream cigar_string;
    int edit_mismatch = 0;
    unsigned ref_pos = 0, read_pos = 0;
    for (int i = 0; i < ez.n_cigar; i++) {
      int count = ez.cigar[i] >> 4;
      char op = "MID"[ez.cigar[i] & 0xf];
      cigar_string << count << op;
      switch (op) {
        case 'M':
          for (int j = 0; j < count; j++, ref_pos++, read_pos++) {
            if (ptr_ref[ref_pos] != ptr_read[read_pos])
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

    a.cigar_string = cigar_string.str();
    free(ez.cigar);
    a.ref_begin = 0;
    region.score = a.sw_score = ez.score;
    a.mismatches = edit_mismatch;
  }

  region.is_aligned = true;
}

void AccAlign::save_region(Read &R, size_t rlen, Region &region,
                           Alignment &a) {
  if (region.is_exact) {
    R.pos = region.pos;
    sprintf(R.cigar, "%uM", (unsigned) rlen);
    R.nm = 0;
  } else {
    R.pos = region.beg + a.ref_begin;
    int cigar_len = a.cigar_string.size();
    strncpy(R.cigar, a.cigar_string.c_str(), cigar_len);
    R.cigar[cigar_len] = '\0';
    rectify_cigar(R.cigar, strlen(R.cigar), R, region);
    R.nm = a.mismatches;
  }

  R.tid = 0;
  for (size_t j = 0; j < name.size(); j++) {
    if (offset[j + 1] > R.pos) {
      R.tid = j;
      break;
    }
  }
  //cerr << "Saving region at pos " << R.pos << " as pos " << R.pos -
  //    offset[R.tid] + 1 << " for read " << R.name << endl;

  R.pos = R.pos - offset[R.tid] + 1;
  R.as = region.score;

}

void AccAlign::align_read(Read &R) {
  auto start = std::chrono::system_clock::now();

  Region region = R.best_region;
  char *s = R.strand == '+' ? R.fwd : R.rev;

  Alignment a;
  size_t rlen = strlen(R.seq);
  if (!region.is_aligned)
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
  const char *pattern = ref.c_str() + region.beg;

  if (toExtend && region.embed_dist) {
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
    rectify_cigar(R.cigar, strlen(R.cigar), R, region);
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

  R.pos = region.pos;
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

void AccAlign::pairdis_filter(vector<Region> &in_regions1, vector<Region> &in_regions2,
                              bool flag1[], bool flag2[],
                              unsigned &best1, unsigned &next1, unsigned &best2, unsigned &next2) {
  unsigned sum_best = 0, sum_next = 0;

  for (unsigned i = 0; i < in_regions1.size(); i++) {
    Region tmp;
    tmp.pos = in_regions1[i].pos - pairdis;
    int start = lower_bound(in_regions2.begin(), in_regions2.end(), tmp,
                            [](const Region &left, const Region &right) {
                              return left.pos < right.pos;
                            }
    ) - in_regions2.begin();

    tmp.pos = in_regions1[i].pos + pairdis;
    int end = upper_bound(in_regions2.begin(), in_regions2.end(), tmp,
                          [](const Region &left, const Region &right) {
                            return left.pos < right.pos;
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
      flag1[i] = 1;
      flag2[j] = 1;
    }
  }
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
  embedding_time = 0;
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
        toExtend = false;
        opn += 1;
        flag = true;
      } else if (av[opn][1] == 'w') {
        useWFA = true;
        opn += 1;
        flag = true;
      } else if (av[opn][1] == 's') {
        enable_softclip = true;
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

  size_t total_begin = time(NULL);

  auto start = std::chrono::system_clock::now();

  AccAlign f(*r);
  f.open_output(g_out);

  if (opn == ac - 1) {
    f.fastq(av[opn], "\0", false);
  } else if (opn == ac - 2) {
    f.fastq(av[opn], av[opn + 1], false);
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
