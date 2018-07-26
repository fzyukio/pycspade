//
// Created by Yukio Fukuzawa on 26/07/18.
//

#include <sstream>
#include "utils.h"

namespace global {
    int num_partitions = 1;
    int maxiter = 2;
    int min_gap = 1;
    int max_gap = INFINITY;
    int max_seq_len = 100;
    int max_iset_len = 100;
    EqGrNode **eqgraph;
    int L2pruning = 0;
    int prepruning = 0;
    int postpruning = 0;
    int pruning_type = Pruning_No;
    int *NumLargeItemset;

    int DBASE_NUM_TRANS;
    int DBASE_MAXITEM;
    double DBASE_AVG_TRANS_SZ;
    double DBASE_AVG_CUST_SZ;
    int DBASE_TOT_TRANS;
    double MINSUP_PER;
    int MINSUP_ABS;

    double FOLLOWTHRESH = 1.0;
    int NUMCLASS = 1;
    long MEMUSED = 0;
    long AVAILMEM = 128 * MBYTE;
}



ostringstream logger;
ostringstream mined;
ostringstream memlog;
ostringstream summary;

int cmp2it(const void *a, const void *b) {
    int *ary = (int *) a;
    int *bry = (int *) b;
    if (ary[0] < bry[0]) return -1;
    else if (ary[0] > bry[0]) return 1;
    else {
        if (ary[1] < bry[1]) return -1;
        else if (ary[1] > bry[1]) return 1;
        else return 0;
    }
}