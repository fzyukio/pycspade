#include <sys/types.h>
#include <unistd.h>
#include <cerrno>

#include "Eqclass.h"
#include "utils.h"

Eqclass::Eqclass(int iset_sz, int eqt) {
    Iset_size = iset_sz;
    Eqtype = eqt;
    theList = new Lists<Itemset *>;
    if (theList == nullptr) {
        throw runtime_error("memory :: Eqclass");
    }
    seqTemplate = seqTemplate2 = 0;
    theList2 = nullptr;
    if (Eqtype == EQCTYP1) {
        theList2 = new Lists<Itemset *>;
        if (theList2 == nullptr) {
            throw runtime_error("memory :: Eqclass");
        }
    }
    global::MEMUSED += sizeof(Eqclass);
}

Eqclass::~Eqclass() {
    if (theList) {
        theList->clear();
        delete theList;
    }
    theList = nullptr;
    if (theList2) {
        theList2->clear();
        delete theList2;
    }
    theList2 = nullptr;
    global::MEMUSED -= sizeof(Eqclass);
}

Itemset *Eqclass::uniqsorted(Itemset *it, CMP_FUNC func) {
    Itemset *rval;
    if (!(rval = theList->find(it, Itemset::Itemcompare))) {
        theList->sortedAscend(it, func);
    }
    return rval;
}

int Eqclass::subseq(Itemset *it) {
    ListNodes<Itemset *> *hd = theList->head();
    for (; hd; hd = hd->next()) {
        if (it->subsequence(hd->item())) {
            return 1;
        }
    }
    return 0;
}


EqGrNode::EqGrNode(int sz) {
    if (sz > 0) {
        theElements = new Array(sz);
        stheElements = new Array(sz);
        _set_sup = new Array *[global::NUMCLASS];
        _seq_sup = new Array *[global::NUMCLASS];
        for (int i = 0; i < global::NUMCLASS; i++) {
            _set_sup[i] = new Array(sz);
            _seq_sup[i] = new Array(sz);
        }
    } else {
        theElements = nullptr;
        stheElements = nullptr;
        _set_sup = nullptr;
        _seq_sup = nullptr;
    }

    freqArray = nullptr;
    freqArraySz = 0;
    theFlg = 0;
    global::MEMUSED += sizeof(EqGrNode);
}

EqGrNode::~EqGrNode() {
    if (theElements) delete theElements;
    if (stheElements) delete stheElements;
    if (_set_sup) {
        for (int i = 0; i < global::NUMCLASS; i++)
            delete _set_sup[i];
    }
    if (_seq_sup) {
        for (int i = 0; i < global::NUMCLASS; i++)
            delete _seq_sup[i];
    }
    if (freqArray) {
        for (int i = 0; i < freqArraySz; i++) delete freqArray[i];
        delete[] freqArray;
    }
    theElements = nullptr;
    theFlg = 0;
    global::MEMUSED -= sizeof(EqGrNode);
}

//assume that elements are sorted in descending order
int EqGrNode::bsearch(int min, int max, FreqIt **freqArray,
                      FreqIt &fit, int recursive) {
    int mid = (max + min) / 2;
    if (max < min) return -1;

    int res = freqArray[mid]->compare(&fit, recursive);
    if (res == 0) return mid;
    else if (res < 0) return bsearch(min, mid - 1, freqArray, fit, recursive);
    else return bsearch(mid + 1, max, freqArray, fit, recursive);
}

int EqGrNode::bsearch(int min, int max, int *itary, int it) {
    int mid = (max + min) / 2;
    if (max < min) return -1;

    if (it == itary[mid]) return mid;
    else if (it < itary[mid]) return bsearch(min, mid - 1, itary, it);
    else return bsearch(mid + 1, max, itary, it);
}


int EqGrNode::find_freqarray(FreqIt &fit, int recursive) {
    if (freqArraySz > 0)
        return bsearch(0, freqArraySz - 1, freqArray, fit, recursive);
    else return 0;
}


ostream &operator<<(ostream &outputStream, EqGrNode &EQ) {
    int i;
    if (EQ.theElements) {
        logger << "SET " << *EQ.theElements << endl;
        for (i = 0; i < global::NUMCLASS; i++)
            logger << "Sup" << i << " : " << *EQ._set_sup[i] << endl;
        logger << "Tot";
        for (i = 0; i < EQ.theElements->size(); i++)
            logger << " " << EQ.get_sup(i);
        logger << endl;
    }
    if (EQ.stheElements) {
        logger << "SEQ " << *EQ.stheElements << endl;
        for (i = 0; i < global::NUMCLASS; i++)
            logger << "SSup" << i << " : " << *EQ._seq_sup[i] << endl;
        logger << "Tot";
        for (i = 0; i < EQ.stheElements->size(); i++)
            logger << " " << EQ.get_seqsup(i);
        logger << endl;
    }

    return outputStream;
}


int FreqIt::compare(Itemset *fit, unsigned int itpl) {
    int i;

    //first compare seqsz, one with larger seqsz is smaller
    if (seqsz > fit->size()) return -1;
    else if (seqsz < fit->size()) return 1;

    //compare items & template bits
    if (seq[0] < (*fit)[0]) return -1;
    else if (seq[0] > (*fit)[0]) return 1;

    int bpos = seqsz - 1;
    int b1, b2;
    for (i = 1; i < seqsz; i++) {
        b1 = GETBIT(templ, bpos - i);
        b2 = GETBIT(itpl, bpos - i);
        if (b1 < b2) return -1;
        else if (b1 > b2) return 1;

        if (seq[i] < (*fit)[i]) return -1;
        else if (seq[i] > (*fit)[i]) return 1;
    }
    return 0;
}


int FreqIt::compare(FreqIt *fit, int recursive) {
    int i;

    //first compare seqsz, one with larger seqsz is smaller
    if (seqsz > fit->seqsz) return -1;
    else if (seqsz < fit->seqsz) return 1;

    //compare items & template bits
    if (seq[seqsz - 1] < fit->seq[fit->seqsz - 1]) return -1;
    else if (seq[seqsz - 1] > fit->seq[fit->seqsz - 1]) return 1;

    int bpos = 0;
    int b1, b2;
    for (i = seqsz - 2; i >= 0; i--, bpos++) {
        b1 = GETBIT(templ, bpos);
        b2 = GETBIT(fit->templ, bpos);
        if (b1 < b2) return -1;
        else if (b1 > b2) return 1;

        if (seq[i] < fit->seq[i]) return -1;
        else if (seq[i] > fit->seq[i]) return 1;
    }
    return 0;
}

//////F1
Array **F1::itsup = nullptr;
int *F1::backidx = nullptr;
int *F1::fidx = nullptr;
int F1::numfreq = 0;
 
