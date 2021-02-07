#pragma once

#include <cstdint>
#include <tuple>
#include <algorithm>
#include <parallel/algorithm>
#include <assert.h>
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <bitset>
#include <unordered_map>
#include <chrono>
#include <zlib.h>
#include <cmath>
#include <omp.h>
#include <string>
#include <thread>
#include <mutex>
#include <future>
#include <chrono>
#include <sys/time.h>
#include <sys/mman.h>
#include <sched.h>
#include <fcntl.h>

#include "tbb/tbb.h"
#include "tbb/flow_graph.h"
#include "tbb/task_scheduler_init.h"
#include "tbb/concurrent_queue.h"

#include "ssw.h"
#include "ssw_cpp.h"
#include "type.h"
#include "const.h"
#include "embedding.h"
#include "accalign.h"
#include "util.h"

#include "gap_affine/affine_wavefront_align.h"
