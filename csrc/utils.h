//
// Created by Yukio Fukuzawa on 26/07/18.
//

#ifndef CSPADE_ORIGINAL_UTILS_H
#define CSPADE_ORIGINAL_UTILS_H

#include <sstream>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <stdexcept>
#include <strstream>
#include <cstring>

enum {
    Pruning_No = 0, Pruning_L2 = 1, Pruning_Zero = 2, Pruning_Follow = 4
};

const long MBYTE = (1024*1024);
const int ITSZ = sizeof(int);

#define min(a, b) ((a) < (b) ? (a) : (b))

#ifndef INFINITY
#define INFINITY 2147483647
#endif

#ifndef bzero
#define bzero(b,len) (memset((b), '\0', (len)), (void) 0)
#endif

//join type
#define LJOIN 0
#define EJOIN 1
#define MJOIN 2

class EqGrNode;

using std::cout;
using std::cerr;
using std::endl;
using std::ostringstream;
using std::ifstream;
using std::vector;
using std::string;
using std::runtime_error;
using std::ostream;
using std::ofstream;
using std::istrstream;

namespace global {
    extern double FOLLOWTHRESH;
    extern long AVAILMEM;
    extern int DBASE_NUM_TRANS;
    extern int DBASE_MAXITEM;
    extern double DBASE_AVG_TRANS_SZ;
    extern double DBASE_AVG_CUST_SZ;
    extern int DBASE_TOT_TRANS;

    extern int num_partitions;
    extern int maxiter;
    extern int *backidx;            // in extl2.cc
    extern int *NumLargeItemset;    // in sequence.cc
    extern EqGrNode **eqgraph;      // -do-
    extern int pruning_type;

    extern int NUMCLASS;
    extern int min_gap;
    extern int max_gap;
    extern int max_seq_len;
    extern int max_iset_len;
    extern int prepruning;
    extern int postpruning;
    extern int L2pruning;

    extern double MINSUP_PER;
    extern int MINSUP_ABS;
    extern long MEMUSED;
}

namespace sequence {
    struct arg_t {
        string name;
        string binf;
        string dataf;
        string idxf;
        string conf;
        string it2f;
        string seqf;
        string classf;

        int num_partitions = 1;
        double min_support = 2;
        int use_ascending = -1;
        bool use_class = false;
        bool do_l2 = false;
        bool use_hash = false;
        int min_gap = 1;
        double maxmem = 128;
        bool recursive = false;
        int pruning_type = Pruning_No;
        int max_gap = 2147483647;
        bool use_window = false;
        int max_seq_len = 100;
        int max_iset_len = 100;

        bool twoseq = false;
        bool use_diff = false;
        bool use_newformat = true;
        bool no_minus_off = true;
    };

    extern arg_t cspade_args;
}

// Logger and mem logger
extern ostringstream logger;
extern ostringstream mined;
extern ostringstream memlog;
extern ostringstream summary;

//extern char *optarg;

extern int cmp2it(const void *a, const void *b);


struct result_t {
    int nsequences;
    bool success;
    string error;
    string mined;
    string logger;
    string summary;
    string memlog;
};


#endif //CSPADE_ORIGINAL_UTILS_H
