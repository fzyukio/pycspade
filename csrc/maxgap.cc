#include <cmath>
#include <iostream>

#include "spade.h"
#include "maxgap.h"
#include "partition.h"

ItBufMgr *IBM;

void ItBufMgr::get_ext_item(int it) {
    int supsz = partition_get_idxsup(it);
    //logger << "GETEXT " << it << " " << supsz << std::endl;
    int *newit = (int *) malloc(supsz * sizeof(int));
    partition_read_item(newit, it);
    _items[F1::fidx[it]]->set_support(supsz);
    _items[F1::fidx[it]]->ival()->set_size(supsz);
    _items[F1::fidx[it]]->ival()->set_array(newit);
}

void process_itemset(Itemset *iset, unsigned int templ, int iter, int min_support_all, bool use_maxgap, bool use_hash,
                     bool recursive, std::ostream& out) {
    int i, it2;
    int lsup, esup;
    unsigned int ntpl;
    Itemset *iset2, *ejoin, *ljoin;
    int it = (*iset)[0];
    if (maxiter < iter) maxiter = iter;
    if (eqgraph[it]) {
        for (i = 0; i < eqgraph[it]->num_elements(); i++) {
            it2 = eqgraph[it]->get_element(i);
            ntpl = templ;
            iset2 = IBM->get_item(it2);
            ljoin = NULL;
            ejoin = prune_decision(iset2, iset, ntpl, EJOIN, use_maxgap, use_hash, recursive);
            if (pruning_type > 1) {
                pre_pruning(ejoin, ntpl, iset, iset2, 0, out);
            }
            if (ejoin)
                get_tmpnewf_intersect(ljoin, ejoin, ljoin, lsup, esup, lsup,
                                      iset2, iset, iter, min_support_all);
            if (ejoin) fill_join(ejoin, iset, iset2);
            if (pruning_type > 1) post_pruning(ejoin, ntpl, out);
            if (ejoin) {
                NumLargeItemset[iter - 1]++;
                if (outputfreq) {
                    ejoin->print_seq(ntpl, out);
                }
                process_itemset(ejoin, ntpl, iter + 1, min_support_all, use_maxgap, use_hash, recursive, out);
                delete ejoin;
            }
        }

        for (i = 0; i < eqgraph[it]->seqnum_elements(); i++) {
            it2 = eqgraph[it]->seqget_element(i);
            ntpl = SETBIT(templ, 1, iter - 2);
            iset2 = IBM->get_item(it2);
            ejoin = NULL;
            ljoin = prune_decision(iset2, iset, ntpl, LJOIN, use_maxgap, use_hash, recursive);
            if (pruning_type > 1) pre_pruning(ljoin, ntpl, iset, iset2, 1, out);
            if (ljoin)
                get_tmpnewf_intersect(ljoin, ejoin, ejoin, lsup, esup, esup,
                                      iset2, iset, iter, min_support_all);
            if (ljoin) fill_join(ljoin, iset, iset2);
            if (pruning_type > 1) post_pruning(ljoin, ntpl, out);
            if (ljoin) {
                NumLargeItemset[iter - 1]++;
                if (outputfreq) ljoin->print_seq(ntpl, out);
                process_itemset(ljoin, ntpl, iter + 1, min_support_all, use_maxgap, use_hash, recursive, out);
                delete ljoin;
            }
        }
    }
}

void process_maxgap(Eqclass *L2, int min_support_all, bool use_maxgap, bool use_hash, bool recursive, std::ostream& out) {
    ListNodes<Itemset *> *hdr = L2->list()->head();
    for (; hdr; hdr = hdr->next()) {
        process_itemset(hdr->item(), L2->templ(), 3, min_support_all, use_maxgap, use_hash, recursive, out);
    }

    hdr = L2->list2()->head();
    for (; hdr; hdr = hdr->next()) {
        process_itemset(hdr->item(), L2->templ2(), 3, min_support_all, use_maxgap, use_hash, recursive, out);
    }
}
