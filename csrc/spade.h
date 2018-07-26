#ifndef SPADE_H
#define SPADE_H

#include <climits>
#include "Itemset.h"
#include "Eqclass.h"

namespace sequence {

    extern void get_tmpnewf_intersect(Itemset *&ljoin, Itemset *&ejoin,
                                      Itemset *&mjoin, int &lcnt, int &ecnt, int &mcnt,
                                      Itemset *it1, Itemset *it2, int iter);

    extern void fill_join(Itemset *join, Itemset *hdr1, Itemset *hdr2);

    extern void
    pre_pruning(Itemset *&iset, unsigned int ptempl, Itemset *it1, Itemset *it2, char use_seq);

    extern void post_pruning(Itemset *&iset, unsigned int templ);

    extern Itemset *prune_decision(Itemset *it1, Itemset *it2, unsigned int ptempl, int jflg);
}

#endif
