#include <cerrno>
#include "Itemset.h"
#include "utils.h"

Itemset::Itemset(int it_sz, int ival_sz, int nclass) {
    theItemset = new Array(it_sz);
    if (theItemset == nullptr) {
        throw runtime_error("memory:: Itemset");
    }
    theIval = new Array(ival_sz);
    if (theIval == nullptr) {
        throw runtime_error("memory:: Ival");
    }

    theSupport = 0;
    clsSup = new int[nclass];
    for (int i = 0; i < nclass; i++) clsSup[i] = 0;

    global::MEMUSED += sizeof(Itemset) + nclass * sizeof(int);
}

Itemset::~Itemset() {
    if (theItemset) delete theItemset;
    if (theIval) delete theIval;
    theItemset = nullptr;
    theSupport = 0;
    global::MEMUSED -= sizeof(Itemset);
}

int Itemset::compare(Itemset &ar2) {
    int len;
    if (size() <= ar2.size()) len = size();
    else len = ar2.size();
    for (int i = 0; i < len; i++) {
        if ((*theItemset)[i] > (*ar2.theItemset)[i]) return 1;
        else if ((*theItemset)[i] < (*ar2.theItemset)[i]) return -1;
    }
    if (size() < ar2.size()) return -1;
    else if (size() > ar2.size()) return 1;
    else return 0;
}

//len must be less than length of both Itemsets
int Itemset::compare(Itemset &ar2, int len) {
    for (int i = 0; i < len; i++) {
        if ((*theItemset)[i] > (*ar2.theItemset)[i]) return 1;
        else if ((*theItemset)[i] < (*ar2.theItemset)[i]) return -1;
    }
    return 0;
}

int Itemset::compare(Array &ar2, int len) {
    for (int i = 0; i < len; i++) {
        if ((*theItemset)[i] > ar2[i]) return 1;
        else if ((*theItemset)[i] < ar2[i]) return -1;
    }
    return 0;
}

int Itemset::compare(Itemset &ar2, int len, unsigned int bvec) {
    int pos = 0;
    int it;
    for (int i = 0; i < len; i++) {
        while (!GETBIT(bvec, pos)) {
            pos++;
        }
        it = (*theItemset)[pos++];
        if (it > (*ar2.theItemset)[i]) return 1;
        else if (it < (*ar2.theItemset)[i]) return -1;
    }
    return 0;
}

int Itemset::subsequence(Itemset *ar) {
    int i, j;
    if (size() > ar->size()) return 0;
    int start = 0;
    for (i = 0; i < size(); i++) {
        for (j = start; j < ar->size(); j++) {
            if ((*theItemset)[i] == (*ar->theItemset)[j]) {
                start = j + 1;
                break;
            }
        }
        if (j >= ar->size()) return 0;
    }
    return 1;
}

ostream &operator<<(ostream &outputStream, Itemset &itemset) {
    outputStream << "ITEM: ";
    outputStream << *itemset.theItemset;
    outputStream << "(" << itemset.theSupport << ")";
    outputStream << "\n";
    return outputStream;
}

void Itemset::print_seq(int itempl) {
    int i;
    int sz = size();
    mined << (*theItemset)[0] << " ";

    for (i = 1; i < sz - 1; i++) {
        if (GETBIT(itempl, sz - 1 - i))
            mined << "-> ";
        mined << (*theItemset)[i] << " ";
    }
    if (GETBIT(itempl, sz - 1 - i))
        mined << "-> ";
    mined << (*theItemset)[sz - 1] << " ";
    mined << "-- " << theSupport;
    for (i = 0; i < global::NUMCLASS; i++)
        mined << " " << clsSup[i];
    mined << " ";
    mined << endl;
}
