#include "header.h"

using namespace std;
const unsigned mod=(1UL<<29)-1;
const unsigned step = 1;
unsigned kmer;
bool low_mem = false;

struct Data{
    uint32_t key, pos;
    Data() : key(-1), pos(-1) {}
    Data(uint32_t k, uint32_t p) : key(k), pos(p) {}
    bool operator()(Data X, Data Y) {
        return X.key == Y.key ? X.pos < Y.pos : X.key < Y.key;
    }
};

class Index{
    private:
        string ref;
    public:
        bool load_ref(const char *F);
        bool make_index(const char *F);
};

bool Index::load_ref(const char *F){
    char code[256], buf[65536];
    for(size_t i=0; i<256; i++) code[i]=4;
    code['A']=code['a']=0; code['C']=code['c']=1; code['G']=code['g']=2; code['T']=code['t']=3;
    cerr << "Loading ref\n";
    FILE *f=fopen(F, "rb");
    if(f==NULL) return false;
    fseek(f, 0, SEEK_END);
    ref.reserve(ftell(f)+1);
    fclose(f);
    f=fopen(F, "rt");
    if(f==NULL) return false;
    while(fgets(buf, 65536, f)!=NULL){
        if(buf[0]=='>') continue;
        for(char *p=buf; *p; p++) if(*p>=33) ref.push_back(*(code+*p));
    }
    fclose(f);
    cerr<<"genome\t"<<ref.size()<<'\n';
    return true;
}

bool Index::make_index(const char *F) {
    size_t limit = ref.size() - kmer + 1;
    size_t vsz;
    if (step == 1)
        vsz = limit;
    else
        vsz = ref.size() / step + 1;

    vector<Data> data(vsz, Data());
    cerr<<"hashing :limit = " << limit << ", vsz = " << vsz << endl;
    size_t i, j;
    uint64_t h;
    bool hasn;

#pragma omp parallel for private(i,j,h,hasn)
    for(i = 0; i < limit; i += step){
        h = 0;
        hasn = false;
        for(j = 0; j < kmer; j++){
            if(ref[i+j] == 4) {
                hasn = true;
                //        break;
            }
            h = (h << 2) + ref[i + j];
        }
        if(!hasn) {
            data[i / step].key = h % mod;
            data[i / step].pos = i;
        }
    }
    cerr << "hash\t" << data.size( ) << endl;

    //XXX: Parallel sort uses lots of memory. Need to fix this. In general, we
    //use 8 bytes per item. Its a waste.
    cerr<<"sorting\n";
    if (low_mem)
        sort(data.begin(), data.end(), Data());
    else
        __gnu_parallel::sort(data.begin(), data.end(), Data());

    //map<unsigned, unsigned> m;
    cerr<<"writing\n";
    string fn=F; fn+=".hash";
    //string keys = F + string(".keys");
    ofstream fo(fn.c_str(), ios::binary);
    //ofstream ko(keys.c_str(), ios::binary);

    // determine the number of valid entries based on first junk entry
    auto joff = std::lower_bound(data.begin(), data.end(), Data(-1, -1), Data());
    size_t eof = joff - data.begin();
    cerr << "Found " << eof << " valid entries out of " <<
        data.size() << " total\n";
    fo.write((char*)&eof, 4);

    // write out keys
    cerr << "Writing posv (" << eof << ")\n";
    for (size_t i = eof; i < data.size(); i++)
        assert(data[i].key == (uint32_t)-1);
    if (!low_mem) {
        uint32_t *buf = new uint32_t[eof];
        for(size_t i = 0; i < eof; i++) {
            buf[i] = data[i].pos;
        }
        fo.write((char *)buf, eof * sizeof(uint32_t));
        delete[] buf;
    } else {
        for(size_t i = 0; i < eof; i++) {
            fo.write((char*)&data[i].pos, 4);
            //ko.write((char *)&data[i].key, 4);
        }
    }

    cerr << "Writing keyv\n";
    size_t last_key = 0, offset;
    if (!low_mem) {
        size_t buf_idx = 0;
        uint32_t *buf = new uint32_t[mod + 1];
        for(size_t i = 0; i < eof; ){
            assert (data[i].key != (uint32_t) -1);
            size_t h = data[i].key, n;
            offset = i;
            for(size_t j = last_key; j <= h; j++) {
                buf[buf_idx] = offset;
                ++buf_idx;
                //fo.write((char*)&offset, 4);
            }
            last_key = h + 1;
            for(n = i + 1; n < eof && data[n].key == h; n++);
            //m[n-i]++;
            i=n;
        }
        offset = eof;
        for(size_t j = last_key; j <= mod; j++) {
            buf[buf_idx] = offset;
            ++buf_idx;
            //fo.write((char*)&offset, 4);
        }
        assert(buf_idx == (mod + 1));
        fo.write((char *)buf, buf_idx * sizeof(uint32_t));
        delete[] buf;
    } else {
        for(size_t i = 0; i < eof; ){
            assert (data[i].key != (uint32_t) -1);
            size_t h = data[i].key, n;
            offset = i;
            for(size_t j = last_key; j <= h; j++) {
                fo.write((char*)&offset, 4);
            }
            last_key = h + 1;
            for(n = i + 1; n < eof && data[n].key == h; n++);
            //m[n-i]++;
            i=n;
        }
        offset = eof;
        for(size_t j = last_key; j <= mod; j++) {
            fo.write((char*)&offset, 4);
        }
    }

    fo.close();
    return true;
}

int main(int ac, char **av){
    if(ac<2){ cerr<<"index [options] <ref.fa>\n";
        cerr<<"options:\n";
        cerr<<"\t-l INT length of seed. [32]\n\n";
        cerr<<"\t-m Use low mem \n\n";
        return 0;
    }
    unsigned kmer_temp=0;
    for(int it=1;it<ac;it++)
    {
        if(strcmp(av[it],"-l")==0)
            kmer_temp=atoi(av[it+1]);
        if(strcmp(av[it],"-m")==0)
            low_mem = true;
    }
    kmer = 32;
    if(kmer_temp!=0)kmer=kmer_temp;

    cerr << "Using kmer length " << kmer << " and step size " << step << endl;

    Index i;
    if(!i.load_ref(av[ac-1])) return 0;
    if(!i.make_index(av[ac-1])) return 0;
    return 0;
}
