#include "header.h"

using namespace std;
const unsigned mod = (1UL << 29) - 1;
const unsigned step = 1;
unsigned kmer;

struct Data {
  uint32_t key, pos;
  Data() : key(-1), pos(-1) {}
  Data(uint32_t k, uint32_t p) : key(k), pos(p) {}

  bool operator()(const Data &X, const Data &Y) const {
    return X.key == Y.key ? X.pos < Y.pos : X.key < Y.key;
  }

};

class Index {
 private:
  string ref;
 public:
  bool load_ref(const char *F);
  bool make_index(const char *F);
  void cal_key(size_t i, vector<Data> &data);
};

bool Index::load_ref(const char *F) {
  char code[256], buf[65536];
  for (size_t i = 0; i < 256; i++)
    code[i] = 4;
  code['A'] = code['a'] = 0;
  code['C'] = code['c'] = 1;
  code['G'] = code['g'] = 2;
  code['T'] = code['t'] = 3;
  cerr << "Loading ref\n";
  FILE *f = fopen(F, "rb");
  if (f == NULL)
    return false;
  fseek(f, 0, SEEK_END);
  ref.reserve(ftell(f) + 1);
  fclose(f);
  f = fopen(F, "rt");
  if (f == NULL)
    return false;
  while (fgets(buf, 65536, f) != NULL) {
    if (buf[0] == '>')
      continue;
    for (char *p = buf; *p; p++)
      if (*p >= 33)
        ref.push_back(*(code + *p));
  }
  fclose(f);
  cerr << "genome\t" << ref.size() << '\n';
  return true;
}

void Index::cal_key(size_t i, vector<Data> &data) {
  uint64_t h = 0;
  bool hasn = false;
  for (unsigned j = 0; j < kmer; j++) {
    if (ref[i + j] == 4) {
      hasn = true;
    }
    h = (h << 2) + ref[i + j];
  }
  if (!hasn) {
    data[i / step].key = h % mod;
    data[i / step].pos = i;
  }
}

class Tbb_cal_key {
  vector<Data> &data;
  Index *index_obj;

 public:
  Tbb_cal_key(vector<Data> &_data, Index *_index_obj) :
      data(_data), index_obj(_index_obj) {}

  void operator()(const tbb::blocked_range<size_t> &r) const {
    for (size_t i = r.begin(); i != r.end(); ++i) {
      index_obj->cal_key(i, data);
    }
  }
};

bool Index::make_index(const char *F) {
  size_t limit = ref.size() - kmer + 1;
  size_t vsz;
  if (step == 1)
    vsz = limit;
  else
    vsz = ref.size() / step + 1;

  vector<Data> data(vsz, Data());
  cerr << "hashing :limit = " << limit << ", vsz = " << vsz << endl;

  tbb::parallel_for(tbb::blocked_range<size_t>(0, limit), Tbb_cal_key(data, this));
  cerr << "hash\t" << data.size() << endl;

  //XXX: Parallel sort uses lots of memory. Need to fix this. In general, we
  //use 8 bytes per item. Its a waste.
  try {
    cerr << "Attempting parallel sorting\n";
    tbb::task_scheduler_init init(tbb::task_scheduler_init::automatic);
    tbb::parallel_sort(data.begin(), data.end(), Data());
  } catch (std::bad_alloc e) {
    cerr << "Fall back to serial sorting (low mem)\n";
    sort(data.begin(), data.end(), Data());
  }

  cerr << "writing\n";
  string fn = F;
  fn += ".hash";
  ofstream fo(fn.c_str(), ios::binary);

  // determine the number of valid entries based on first junk entry
  auto joff = std::lower_bound(data.begin(), data.end(), Data(-1, -1), Data());
  size_t eof = joff - data.begin();
  cerr << "Found " << eof << " valid entries out of " <<
       data.size() << " total\n";
  fo.write((char *) &eof, 4);

  // write out keys
  for (size_t i = eof; i < data.size(); i++)
    assert(data[i].key == (uint32_t) -1);
  try {
    cerr << "Fast writing posv (" << eof << ")\n";
    uint32_t *buf = new uint32_t[eof];
    for (size_t i = 0; i < eof; i++) {
      buf[i] = data[i].pos;
    }
    fo.write((char *) buf, eof * sizeof(uint32_t));
    delete[] buf;
  } catch (std::bad_alloc e) {
    cerr << "Fall back to slow writing posv due to low mem.\n";
    for (size_t i = 0; i < eof; i++) {
      fo.write((char *) &data[i].pos, 4);
    }
  }

  size_t last_key = 0, offset;
  try {
    cerr << "Fast writing keyv\n";
    size_t buf_idx = 0;
    uint32_t *buf = new uint32_t[mod + 1];
    for (size_t i = 0; i < eof;) {
      assert (data[i].key != (uint32_t) -1);
      size_t h = data[i].key, n;
      offset = i;
      for (size_t j = last_key; j <= h; j++) {
        buf[buf_idx] = offset;
        ++buf_idx;
      }
      last_key = h + 1;
      for (n = i + 1; n < eof && data[n].key == h; n++);
      i = n;
    }
    offset = eof;
    for (size_t j = last_key; j <= mod; j++) {
      buf[buf_idx] = offset;
      ++buf_idx;
    }
    assert(buf_idx == (mod + 1));
    fo.write((char *) buf, buf_idx * sizeof(uint32_t));
    delete[] buf;
  } catch (std::bad_alloc e) {
    cerr << "Fall back to slow writing keyv (low mem)\n";
    for (size_t i = 0; i < eof;) {
      assert (data[i].key != (uint32_t) -1);
      size_t h = data[i].key, n;
      offset = i;
      for (size_t j = last_key; j <= h; j++) {
        fo.write((char *) &offset, 4);
      }
      last_key = h + 1;
      for (n = i + 1; n < eof && data[n].key == h; n++);
      i = n;
    }
    offset = eof;
    for (size_t j = last_key; j <= mod; j++) {
      fo.write((char *) &offset, 4);
    }
  }

  cerr << "Indexing complete\n";
  fo.close();
  return true;
}

int main(int ac, char **av) {
  if (ac < 2) {
    cerr << "index [options] <ref.fa>\n";
    cerr << "options:\n";
    cerr << "\t-l INT length of seed [32]\n";
    return 0;
  }
  unsigned kmer_temp = 0;
  for (int it = 1; it < ac; it++) {
    if (strcmp(av[it], "-l") == 0)
      kmer_temp = atoi(av[it + 1]);
  }
  kmer = 32;
  if (kmer_temp != 0)
    kmer = kmer_temp;

  cerr << "Using kmer length " << kmer << " and step size " << step << endl;

  Index i;
  if (!i.load_ref(av[ac - 1]))
    return 0;
  if (!i.make_index(av[ac - 1]))
    return 0;
  return 0;
}
