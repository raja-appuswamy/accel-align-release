#include "header.h"

extern uint8_t code[256];

class RefParser {
  string &ref;

public:
  RefParser(string &ref) : ref(ref) {}

  void operator()(const tbb::blocked_range<size_t> &r) const {
    for (size_t i = r.begin(); i != r.end(); ++i)
      ref[i] = *(code + ref[i]);
  }
};

void Reference::load_index(const char *F) {
  string fn = string(F) + ".hash";

  cerr << "loading hashtable from " << fn << endl;
  ifstream fi;
  fi.open(fn.c_str(), ios::binary);
  if (!fi) {
    cerr << "Unable to open index file " << fn << endl;
    exit(0);
  }
  fi.read((char *) &nposv, 4);
  fi.close();
  nkeyv = MOD + 1;

  cerr << "Mapping keyv of size: " << nkeyv * 4 <<
       " and posv of size " << (size_t) nposv * 4 <<
       " from index file " << fn << endl;
  size_t posv_sz = (size_t) nposv * sizeof(uint32_t);
  size_t keyv_sz = (size_t) nkeyv * sizeof(uint32_t);
  int fd = open(fn.c_str(), O_RDONLY);
  char *base = reinterpret_cast<char *>(
      mmap(NULL, 4 + posv_sz + keyv_sz, PROT_READ, MAP_PRIVATE |
          MAP_POPULATE, fd, 0));
  assert(base != MAP_FAILED);
  posv = (uint32_t *) (base + 4);
  keyv = posv + nposv;
  cerr << "Mapping done" << endl;
  cerr << "done loading hashtable\n";

}

Reference::Reference(const char *F) {
  struct timeval start, end;
  gettimeofday(&start, NULL);

  // Start index load in parallel
  thread t(&Reference::load_index, this, F);

  size_t ref_size = 0;
  string bref(F);
  bref += ".bref";
  ifstream brefin(bref.c_str(), ios::in);
  if (brefin) {
    /* Experimentation revealed that using prebuild binary ref leads to
     * embeding procedure taking more latency on memory stalls. Might have
     * something to do with ref now being loaded in cache. So while support
     * is there so that we can test it on other servers, it is disabled.
     */
    cout << "WARNING: USING PREBUILD BINARY REF.\n";
    cout << "-----------------------------------\n";
    cout << "Loading index and binary ref from " << bref << endl;
    size_t count;
    brefin >> count;
    cout << "Done loading " << count << " names " << endl;
    for (size_t i = 0; i < count; ++i) {
      string tmp;
      brefin >> tmp;
      name.push_back(tmp);
    }
    brefin >> count;
    for (size_t i = 0; i < count; ++i) {
      uint32_t tmp;
      brefin >> tmp;
      offset.push_back(tmp);
    }
    cout << "Done loading " << count << " offsets " << endl;
    brefin >> ref_size;
    ref.resize(ref_size);
    brefin.read(&ref[0], ref_size);
    cout << "Done loading binary ref of sz " << ref_size << endl;
    brefin.close();
  } else {
    FILE *f = fopen(F, "rb");
    fseek(f, 0, SEEK_END);
    ref.reserve(ftell(f) + 1);
    ref = "";
    fclose(f);

    ifstream fi(F);
    if (!fi) {
      cerr << "fail to open " << F << '\n';
      return;
    }
    cerr << "loading references\n";
    string buf = "";
    while (!fi.eof()) {
      getline(fi, buf);
      if (!buf.size())
        continue;

      if (buf[0] == '>') {
        string tmp(buf.c_str() + 1);
        istringstream iss(tmp);
        vector<string> parsed((istream_iterator<string>(iss)),
                              istream_iterator<string>());
        //cerr<<"Loading chromosome " << parsed[0] << "\n";
        name.push_back(parsed[0]);
        offset.push_back(ref_size);
      } else {
        ref += buf;
        ref_size += buf.size();
      }
    }

    // process last chromosome
    cerr << "Parsing reference.\n";
    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, ref_size),
        RefParser(ref));
    fi.close();
    offset.push_back(ref_size);
    cerr << "Loaded reference size: " << ref_size << endl;

#if FLUSH_BINARY_REF
    // optionally write out the reference to a binary file
    cout << "Writing out the binary ref to " << bref << endl;
    ofstream refout(bref.c_str(), ios::out);
    refout << name.size() << endl;
    for (string &n : name)
        refout << n << endl;
    refout << offset.size() << endl;
    for (uint32_t off : offset)
        refout << off << endl;
    refout << ref_size << endl;
    refout << ref << endl;
    refout.close();
#endif
  }

  // wait for index load to finish
  t.join();

  gettimeofday(&end, NULL);
  cerr << "Setup reference in " <<
       compute_elapsed(&start, &end) / 1000000 << " secs\n";
}

Reference::~Reference() {
  size_t posv_sz = (size_t) nposv * sizeof(uint32_t);
  size_t keyv_sz = (size_t) nkeyv * sizeof(uint32_t);
  char *base = (char *) posv - 4;
  int r = munmap(base, posv_sz + keyv_sz);
  assert(r == 0);
}


