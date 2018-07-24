#include <cerrno>
#include <iostream>
#include <unistd.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <stdexcept>
#include <sstream>

#include "utils.h"
#include "extl2.h"

#define seqitcntbufsz 4086

int seqbuf[seqitcntbufsz];
int seqpos = 0;

unsigned int ***set_sup, ***seq_sup;
invdb *invDB;
int EXTBLKSZ;

invdb::invdb(int sz) {
    int i;
    numcust = sz;
    curit = (int **) malloc(numcust * sizeof(int *));
    curcnt = (int *) malloc(numcust * ITSZ);
    curcid = (int *) malloc(numcust * ITSZ);
    curitsz = (int *) malloc(numcust * ITSZ);
    int ttval = (int) (DBASE_AVG_CUST_SZ * DBASE_AVG_TRANS_SZ);
    for (i = 0; i < numcust; i++) {
        curitsz[i] = ttval;
        curit[i] = (int *) malloc(curitsz[i] * ITSZ);
        curcnt[i] = 0;
        curcid[i] = NOCLASS;
    }
}

invdb::~invdb() {
    int i;
    for (i = 0; i < numcust; i++) {
        free(curit[i]);
    }
    free(curit);
    free(curcnt);
    free(curcid);
    free(curitsz);
}

void invdb::incr(int sz) {
    int oldsz = numcust;
    numcust = sz;

    //logger << "INCR SIZE " << sz << " "<< oldsz << std::endl;

    curit = (int **) realloc(curit, numcust * sizeof(int *));
    curcnt = (int *) realloc(curcnt, numcust * ITSZ);
    curcid = (int *) realloc(curcid, numcust * ITSZ);
    curitsz = (int *) realloc(curitsz, numcust * ITSZ);
    if (curit == NULL || curcnt == NULL || curitsz == NULL || curcid == NULL) {
        throw std::runtime_error("REALLCO  curit");
    }

    int i;
    int ttval = (int) (DBASE_AVG_CUST_SZ * DBASE_AVG_TRANS_SZ);
    for (i = oldsz; i < numcust; i++) {
        curitsz[i] = ttval;
        curit[i] = (int *) malloc(curitsz[i] * ITSZ);
        curcnt[i] = 0;
        curcid[i] = NOCLASS;
    }
}

void invdb::incr_curit(int midx) {
    //logger << "OLD " << curitsz[midx] << std::endl;
    curitsz[midx] = (int) (2 * curitsz[midx]);
    //logger << "NEW " << curitsz[midx] << std::endl;
    curit[midx] = (int *) realloc(curit[midx], curitsz[midx] * ITSZ);
    if (curit[midx] == NULL) {
        throw std::runtime_error("REALLCO  curit");
    }
    //logger << "INCR " << midx << " " << curitsz[midx] << std::endl;
}

void print_idlist(int *ival, int supsz) {
    std::ostringstream& logger = *_logger;
    int i, cid, cnt;

    if (supsz > 0) {
        cid = ival[0];
        cnt = 0;
        for (i = 0; i < supsz;) {
            if (cid == ival[i]) {
                cnt++;
                i += 2;
            } else {
                logger << cid << " " << cnt << " ";
                cid = ival[i];
                cnt = 0;
            }
        }
        logger << cid << " " << cnt;
    }
}

int make_l1_pass(std::ostream& out) {
    int i, j;
    int *sup, supsz;
    int bsz = 100;

    F1::init();
    F1::backidx = (int *) malloc(bsz * ITSZ);
    if (F1::backidx == NULL) {
        throw std::invalid_argument("F1::BACKIDX NULL");
    }
    F1::fidx = new int[DBASE_MAXITEM];

    F1::numfreq = 0;
    int ivalsz = 100;
    int *ival = (int *) malloc(ivalsz * ITSZ);
    if (ival == NULL) {
        throw std::invalid_argument("IVAL NULL");
    }
    int tt = 0;
    for (i = 0; i < DBASE_MAXITEM; i++) {
        supsz = partition_get_idxsup(i);
        if (ivalsz < supsz) {
            ivalsz = supsz;
            ival = (int *) realloc(ival, ivalsz * ITSZ);
            if (ival == NULL) {
                throw std::invalid_argument("IVAL NULL");
            }
        }
        partition_read_item(ival, i);
        for (j = 0; j < NUMCLASS; j++) ClassInfo::TMPE[j] = 0;

        int cid = -1;
        for (j = 0; j < supsz; j += 2) {
            if (cid != ival[j])
                ClassInfo::TMPE[ClassInfo::getcls(ival[j])]++;
            cid = ival[j];
            //if (tt < cid) tt = cid;
        }

        char lflg = 0;
        F1::fidx[i] = -1;       // default init
        for (j = 0; j < NUMCLASS; j++) {
            if (ClassInfo::TMPE[j] >= ClassInfo::MINSUP[j]) {
                lflg = 1;
                if (F1::numfreq + 1 > bsz) {
                    bsz = 2 * bsz;
                    F1::backidx = (int *) realloc(F1::backidx, bsz * ITSZ);
                    if (F1::backidx == NULL) {
                        throw std::invalid_argument("F1::BACKIDX NULL");
                    }
                }
                F1::backidx[F1::numfreq] = i;
                F1::fidx[i] = F1::numfreq;
                F1::numfreq++;
                break;
            }
        }

        if (lflg) {
            if (outputfreq) out << i << " --";
            for (j = 0; j < NUMCLASS; j++) {
                F1::add_sup(ClassInfo::TMPE[j], j);
                if (outputfreq) out << " " << ClassInfo::TMPE[j];
            }
            if (outputfreq) out << " " << F1::get_sup(i) << " ";
            if (print_tidlist) print_idlist(ival, supsz);
            if (outputfreq) out << std::endl;
        }
    }

    F1::backidx = (int *) realloc(F1::backidx, F1::numfreq * ITSZ);
    if (F1::backidx == NULL && F1::numfreq != 0) {
        throw std::invalid_argument("F1::BACKIDX NULL");
    }

    free(ival);

    return F1::numfreq;
}

void suffix_add_item_eqgraph(char use_seq, int it1, int it2) {
    if (eqgraph[it2] == NULL) {
        eqgraph[it2] = new EqGrNode(2);
    }
    if (use_seq) eqgraph[it2]->seqadd_element(it1);
    else eqgraph[it2]->add_element(it1);
}


void process_cust_invert(int cid, int curcnt, int *curit) {
    int i, j, k, l;
    int nv1, nv2, diff;
    int it1, it2;

    //logger << "PROCESS " << cid << " " << curcnt << " -- ";
    //for (i=0; i < curcnt; ++i)
    //    logger << " " << curit[i];
    //logger << std::endl;

    for (i = 0; i < curcnt; i = nv1) {
        nv1 = i;
        it1 = curit[i];
        while (it1 == curit[nv1] && nv1 < curcnt) nv1 += 2;
        for (j = i; j < curcnt; j = nv2) {
            nv2 = j;
            it2 = curit[j];
            while (it2 == curit[nv2] && nv2 < curcnt) nv2 += 2;
            if (seq_sup[it1] && curit[i + 1] + min_gap <= curit[nv2 - 1]) {
                for (k = i, l = j; k < nv1 && l < nv2;) {
                    diff = curit[l + 1] - curit[k + 1];
                    if (diff < min_gap) l += 2;
                    else if (diff > max_gap) k += 2;
                    else {
                        seq_sup[it1][it2][ClassInfo::getcls(cid)]++;
                        break;
                    }
                }
            }

            if (j > i) {
                if (seq_sup[it2] && curit[j + 1] + min_gap <= curit[nv1 - 1]) {
                    for (k = j, l = i; k < nv2 && l < nv1;) {
                        diff = curit[l + 1] - curit[k + 1];
                        if (diff < min_gap) l += 2;
                        else if (diff > max_gap) k += 2;
                        else {
                            seq_sup[it2][it1][ClassInfo::getcls(cid)]++;
                            break;
                        }
                    }
                }

                if (set_sup[it1]) {
                    for (k = i, l = j; k < nv1 && l < nv2;) {
                        if (curit[k + 1] > curit[l + 1]) l += 2;
                        else if (curit[k + 1] < curit[l + 1]) k += 2;
                        else {
                            set_sup[it1][it2 - it1 - 1][ClassInfo::getcls(cid)]++;
                            break;
                        }
                    }
                }
            }
        }
    }
}

void process_invert(int pnum) {
    int i, k;
    int minv, maxv;
    partition_get_minmaxcustid(F1::backidx, F1::numfreq, pnum, minv, maxv);
    //logger << "MVAL " << minv << " " << maxv << std::endl;
    if (invDB->numcust < maxv - minv + 1)
        invDB->incr(maxv - minv + 1);

    int supsz;
    int ivalsz = 0;
    int *ival = NULL;
    for (i = 0; i < F1::numfreq; i++) {
        supsz = partition_get_lidxsup(pnum, F1::backidx[i]);
        if (ivalsz < supsz) {
            ivalsz = supsz;
            ival = (int *) realloc(ival, ivalsz * ITSZ);
            if (ival == NULL) {
                throw std::runtime_error("IVAL NULL");
            }
        }
        //logger << "READ " << i << " " << F1::backidx[i] << std::endl;
        partition_lclread_item(ival, pnum, F1::backidx[i]);

        int cid;
        int midx;
        for (int pos = 0; pos < supsz; pos += 2) {
            //if (cid != ival[pos]){
            cid = ival[pos];
            midx = cid - minv;

            if (invDB->curcnt[midx] + 2 > invDB->curitsz[midx]) {
                invDB->incr_curit(midx);
            }
            invDB->curcid[midx] = cid;
            invDB->curit[midx][invDB->curcnt[midx]++] = i;
            invDB->curit[midx][invDB->curcnt[midx]++] = ival[pos + 1];
        }
    }
    for (k = 0; k < maxv - minv + 1; k++) {
        if (invDB->curcnt[k] > 0) {
            process_cust_invert(invDB->curcid[k], invDB->curcnt[k], invDB->curit[k]);
        }
        invDB->curcnt[k] = 0;
        invDB->curcid[k] = NOCLASS;
    }
}


//return 1 to prune, else 0
char extl2_pre_pruning(int totsup, int it, int pit, char use_seq, unsigned int *clsup = NULL) {
    std::ostringstream& logger = *_logger;
    float conf, conf2;
    int itsup;
    if (pruning_type == NOPRUNING) return 0;
    if (use_seq) return 0;
    if (GETBIT(pruning_type, FOLLOWPRUNING - 1)) {
        itsup = F1::get_sup(it);
        conf = (1.0 * totsup) / itsup;
        conf2 = (1.0 * totsup) / F1::get_sup(pit);
        //logger << "prune " << conf << " " << totsup << " " << itsup << std::endl;
        if (conf >= FOLLOWTHRESH || conf2 >= FOLLOWTHRESH) {
            if (outputfreq) {
                logger << "PRUNE_EXT " << pit << (use_seq ? " -2 " : " ")
                     << it << " -1 " << totsup;
                for (int i = 0; i < NUMCLASS; i++)
                    logger << " " << clsup[i];
                logger << std::endl;
            }
            prepruning++;
            return 1;
        }
    }
    return 0;
}

void get_F2(int &l2cnt) {
    int i, j, k;
    int fcnt;
    char lflg;
    char use_seq;

    for (j = 0; j < F1::numfreq; j++) {
        //seqpos = 0;
        if (set_sup[j]) {
            use_seq = 0;
            for (k = j + 1; k < F1::numfreq; k++) {
                lflg = 0;
                for (i = 0; i < NUMCLASS; i++) {
                    fcnt = set_sup[j][k - j - 1][i];
                    if (fcnt >= ClassInfo::MINSUP[i]) {
                        //seqbuf[seqpos++] = F1::backidx[k];
                        lflg = 1;
                        break;
                    }
                }
                if (lflg) {
                    fcnt = 0;
                    for (i = 0; i < NUMCLASS; i++) {
                        fcnt += set_sup[j][k - j - 1][i];
                    }
                    if (!extl2_pre_pruning(fcnt, F1::backidx[k],
                                           F1::backidx[j], use_seq, set_sup[j][k - j - 1])) {
                        suffix_add_item_eqgraph(use_seq, F1::backidx[j], F1::backidx[k]);
                        //logger << F1::backidx[j] << " " << F1::backidx[k];
                        for (i = 0; i < NUMCLASS; i++) {
                            int ffcnt = set_sup[j][k - j - 1][i];
                            //logger << " " << fcnt;
                            eqgraph[F1::backidx[k]]->add_sup(ffcnt, i);
                        }
                        //logger << std::endl;
                        l2cnt++;
                    }
                    //if (outputfreq){
                    //   logger << F1::backidx[j] << " " << F1::backidx[k] << " -2";
                    //   for (i=0; i < NUMCLASS; i++)
                    //      logger << " " << set_sup[j][k-j-1][i];
                    //   logger << " " << fcnt << std::endl;
                    //}
                }
            }
        }
        //if (seqpos > 0) add_item_eqgraph(F1::backidx[j], 0);
        //seqpos = 0;
        if (seq_sup[j]) {
            use_seq = 1;
            for (k = 0; k < F1::numfreq; k++) {
                lflg = 0;
                for (i = 0; i < NUMCLASS; i++) {
                    fcnt = seq_sup[j][k][i];
                    if (fcnt >= ClassInfo::MINSUP[i]) {
                        //seqbuf[seqpos++] = F1::backidx[k];
                        lflg = 1;
                        break;
                    }
                }
                if (lflg) {
                    fcnt = 0;
                    for (i = 0; i < NUMCLASS; i++) {
                        fcnt += seq_sup[j][k][i];
                    }
                    if (!extl2_pre_pruning(fcnt, F1::backidx[k],
                                           F1::backidx[j], use_seq, seq_sup[j][k])) {
                        suffix_add_item_eqgraph(use_seq, F1::backidx[j], F1::backidx[k]);
                        l2cnt++;
                        //logger << F1::backidx[j] << " -1 " << F1::backidx[k] << " -2 ";
                        for (i = 0; i < NUMCLASS; i++) {
                            int ffcnt = seq_sup[j][k][i];
                            //logger << fcnt << " ";
                            eqgraph[F1::backidx[k]]->add_seqsup(ffcnt, i);
                        }
                        //logger << std::endl;
                    }
                    //if (outputfreq){
                    //  logger << F1::backidx[j] << " -1 " << F1::backidx[k] << " -2";
                    //   for (i=0; i < NUMCLASS; i++)
                    //      logger << " " << seq_sup[j][k][i];
                    //   logger << " " << fcnt << std::endl;
                    //}
                }
            }
        }
        //if (seqpos > 0) add_item_eqgraph(F1::backidx[j], 1);
    }
}


int make_l2_pass() {
    int i, j;

    int l2cnt = 0;
    int mem_used = 0;

    EXTBLKSZ = num_partitions + (DBASE_NUM_TRANS + num_partitions - 1) / num_partitions;
    int tsz = (int) (DBASE_AVG_CUST_SZ * DBASE_AVG_TRANS_SZ);
    invDB = new invdb(EXTBLKSZ);

    set_sup = new unsigned int **[F1::numfreq];        // support for 2-itemsets
    bzero((char *) set_sup, F1::numfreq * sizeof(unsigned int **));
    mem_used += F1::numfreq * sizeof(unsigned int **);
    seq_sup = new unsigned int **[F1::numfreq];        // support for 2-itemsets
    bzero((char *) seq_sup, F1::numfreq * sizeof(unsigned int **));
    mem_used += F1::numfreq * sizeof(unsigned int **);

    int low, high;

    int itsz = sizeof(unsigned int);
    for (low = 0; low < F1::numfreq; low = high) {
        for (high = low; high < F1::numfreq &&
                         mem_used + ((2 * F1::numfreq - high - 1) * itsz * NUMCLASS) <
                         AVAILMEM; high++) {
            if (max_iset_len > 1 && F1::numfreq - high - 1 > 0) {
                set_sup[high] = new unsigned int *[F1::numfreq - high - 1];
                for (i = 0; i < F1::numfreq - high - 1; i++) {
                    set_sup[high][i] = new unsigned int[NUMCLASS];
                    for (j = 0; j < NUMCLASS; j++) set_sup[high][i][j] = 0;
                }
                mem_used += (F1::numfreq - high - 1) * itsz * NUMCLASS;
            }
            if (max_seq_len > 1) {
                seq_sup[high] = new unsigned int *[F1::numfreq];
                for (i = 0; i < F1::numfreq; i++) {
                    seq_sup[high][i] = new unsigned int[NUMCLASS];
                    for (j = 0; j < NUMCLASS; j++) seq_sup[high][i][j] = 0;
                }
                mem_used += F1::numfreq * itsz * NUMCLASS;
            }
        }
        for (int p = 0; p < num_partitions; p++) {
            process_invert(p);
        }
        get_F2(l2cnt);
        // reclaim memory
        for (i = low; i < high; i++) {
            if (set_sup[i]) {
                for (j = 0; j < F1::numfreq - i - 1; j++)
                    delete[] set_sup[i][j];
                delete[] set_sup[i];
                mem_used -= (F1::numfreq - i - 1) * itsz * NUMCLASS;
            }
            set_sup[i] = NULL;
            if (seq_sup[i]) {
                for (j = 0; j < F1::numfreq; j++)
                    delete[] seq_sup[i][j];
                delete[] seq_sup[i];
                mem_used -= F1::numfreq * itsz * NUMCLASS;
            }
            seq_sup[i] = NULL;
        }
    }

    delete[] set_sup;
    delete[] seq_sup;
    delete invDB;
    return l2cnt;
}


void get_l2file(char *fname, char use_seq, int &l2cnt) {
    int *cntary;
    int fd = open(fname, O_RDONLY);
    if (fd < 1) {
        throw std::runtime_error("can't open l2 file");
    }
    int flen = lseek(fd, 0, SEEK_END);
    if (flen > 0) {
#ifndef DEC
        cntary = (int *) mmap((char *) NULL, flen, PROT_READ,
                              MAP_PRIVATE, fd, 0);
#else
        cntary = (int *) mmap((char *)NULL, flen, PROT_READ,
                               (MAP_FILE|MAP_VARIABLE|MAP_PRIVATE), fd, 0);
#endif
        if (cntary == (int *) -1) {
            throw std::runtime_error("MMAP ERROR:cntary");
        }

        // build eqgraph -- large 2-itemset relations
        int lim = flen / ITSZ;
        //int oldit = -1;
        //seqpos = 0;
        char lflg = 0;
        int i, j;
        for (i = 0; i < lim; i += 3) {
            lflg = 0;
            for (j = 0; j < NUMCLASS; j++) {
                if (cntary[i + 2] >= ClassInfo::MINSUP[j]) {
                    lflg = 1;
                    //logger << "FILELARGE " << cntary[i] << ((use_seq)?"->":" ")
                    //     << cntary[i+1] << " SUPP " << cntary[i+2] << std::endl;
                    break;
                }
            }
            if (lflg) {
                if (!extl2_pre_pruning(cntary[i + 2], cntary[i + 1],
                                       cntary[i], use_seq)) {
                    suffix_add_item_eqgraph(use_seq, cntary[i], cntary[i + 1]);
                    l2cnt++;
                    //assign sup to a single class, sice we don't know breakup
                    if (use_seq) eqgraph[cntary[i + 1]]->add_seqsup(cntary[i + 2], 0);
                    else eqgraph[cntary[i + 1]]->add_sup(cntary[i + 2], 0);
                    for (j = 1; j < NUMCLASS; j++)
                        if (use_seq) eqgraph[cntary[i + 1]]->add_seqsup(0, j);
                        else eqgraph[cntary[i + 1]]->add_sup(0, j);
                }
                //if (outputfreq){
                //  logger << cntary[i] << ((use_seq)?" -1 ":" ") << cntary[i+1]
                //        << " -2";
                //   logger << " " << cntary[i+2];
                //   for (j=1; j < NUMCLASS; j++) logger << " 0";
                //   logger << " " << cntary[i+2];
                //   logger << std::endl;
                //}
            }
        }
        //if (seqpos > 0) add_item_eqgraph(oldit, use_seq);
        munmap((caddr_t) cntary, flen);
    }
    close(fd);
}

int get_file_l2(char *it2f, char *seqf) {
    int l2cnt = 0;

    if (max_iset_len > 1) get_l2file(it2f, 0, l2cnt);
    if (max_seq_len > 1) get_l2file(seqf, 1, l2cnt);
    std::ostringstream os;
    os << "L2 : " << l2cnt;

    throw std::runtime_error(os.str());
    return l2cnt;
}
















