#include <cerrno>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>
#include <fcntl.h>
#include <unistd.h>
#include <sys/time.h>
#include "Eqclass.h"
#include "extl2.h"
#include "maxgap.h"
#include "utils.h"
#include "makebin.h"
#include "getconf.h"
#include "exttpose.h"


namespace sequence {
    int FreqArrayPos = 0;
    FreqIt **FreqArray;
    int FreqArraySz = 100;
    int num_intersect = 0;
    bool use_maxgap = false;

    arg_t cspade_args;

    Array *interval, *interval2, *interval3;

    double L2ISECTTIME = 0, EXTL1TIME = 0, EXTL2TIME = 0;
    int maxitemsup;

    void process_cluster1(Eqclass *cluster, Lists<Eqclass *> *LargeL, int iter);

    void add_freq(Itemset *it, int templ) {
        auto *freq = new FreqIt(it->itemset()->array(), it->size(), templ);
        if (FreqArrayPos + 1 >= FreqArraySz) {
            FreqArraySz = (int) (1.5 * FreqArraySz);
            FreqArray = (FreqIt **) realloc(FreqArray, FreqArraySz * sizeof(FreqIt *));
            if (FreqArray == nullptr) {
                throw runtime_error("no mmeory fro FREqArray ");
            }
        }
        FreqArray[FreqArrayPos++] = freq;
    }


    arg_t &populate_names(arg_t &_args) {
        _args.binf = _args.name + ".data";
        _args.dataf = _args.name + ".tpose";
        _args.idxf = _args.name + ".idx";
        _args.conf = _args.name + ".conf";
        _args.it2f = _args.name + ".2it";
        _args.seqf = _args.name + ".2seq";
        _args.classf = _args.name + ".class";
        return _args;
    }

    void populate_global() {

        if (use_maxgap)
            cspade_args.use_hash = 0;

        if (cspade_args.max_gap < cspade_args.min_gap)
            cspade_args.max_gap = cspade_args.min_gap;

        global::num_partitions = cspade_args.num_partitions;
        global::AVAILMEM = (long) cspade_args.maxmem * MBYTE;
        global::max_gap = cspade_args.max_gap;
        global::min_gap = cspade_args.min_gap;
        global::max_seq_len = cspade_args.max_seq_len;
        global::max_iset_len = cspade_args.max_iset_len;
        global::pruning_type = cspade_args.pruning_type;
        global::max_seq_len = cspade_args.max_seq_len;
        global::max_iset_len = cspade_args.max_iset_len;

        int c = open(cspade_args.conf.c_str(), O_RDONLY);
        if (c < 0) {
            throw runtime_error("ERROR: invalid conf file\n");
        }
        read(c, &global::DBASE_NUM_TRANS, ITSZ);

        if (cspade_args.min_support < 1) {
            global::MINSUP_PER = cspade_args.min_support;
            global::MINSUP_ABS = static_cast<int>(lround(global::MINSUP_PER * global::DBASE_NUM_TRANS));
        } else {
            global::MINSUP_ABS = static_cast<int>(cspade_args.min_support);
            global::MINSUP_PER = static_cast<double >(global::MINSUP_ABS) / global::DBASE_NUM_TRANS;
        }

        logger << "min_support_all " << global::MINSUP_ABS << " out of " << global::DBASE_NUM_TRANS
               << " sequences" << endl;
        read(c, &global::DBASE_MAXITEM, ITSZ);
        read(c, &global::DBASE_AVG_CUST_SZ, sizeof(double));
        read(c, &global::DBASE_AVG_TRANS_SZ, sizeof(double));
        read(c, &global::DBASE_TOT_TRANS, ITSZ);
        close(c);
    }

    void get_2newf_intersect(Itemset *ljoin, Itemset *ejoin,
                             int *it1, int *it2, int sup1, int sup2) {
        int i, j, k, l;
        int nval1, nval2, diff;
        int lflge;

        num_intersect++;

        int icid, jcid;
        for (i = 0, j = 0; i < sup1 && j < sup2;) {
            icid = it1[i];
            jcid = it2[j];
            if (icid > jcid) {
                j += 2;
            } else if (icid < jcid) {
                i += 2;
            } else {
                nval1 = i;
                nval2 = j;
                while (it1[i] == it1[nval1] && nval1 < sup1) nval1 += 2;
                while (it2[j] == it2[nval2] && nval2 < sup2) nval2 += 2;
                if (ljoin && it1[i + 1] + global::min_gap <= it2[nval2 - 1]) {
                    //add tid
                    lflge = 0;
                    for (k = i, l = j; k < nval1 && l < nval2;) {
                        diff = it2[l + 1] - it1[k + 1];
                        if (diff < global::min_gap) l += 2;
                        else if (diff > global::max_gap) k += 2;
                        else {
                            ljoin->ival()->optadd(icid);
                            ljoin->ival()->optadd(it1[k + 1]);
                            lflge = 1;
                            k += 2;
                        }
                    }
                    if (lflge) {
                        ljoin->increment_support();
                        ljoin->increment_cls_support(ClassInfo::getcls(icid));
                    }
                }

                if (ejoin) {
                    lflge = 0;
                    for (k = i, l = j; k < nval1 && l < nval2;) {
                        if (it1[k + 1] < it2[l + 1]) k += 2;
                        else if (it1[k + 1] > it2[l + 1]) l += 2;
                        else {
                            ejoin->ival()->optadd(icid);
                            ejoin->ival()->optadd(it2[l + 1]);
                            lflge = 1;
                            k += 2;
                            l += 2;
                        }
                    }
                    if (lflge) {
                        ejoin->increment_support();
                        ejoin->increment_cls_support(ClassInfo::getcls(icid));
                    }
                }
                i = nval1;
                j = nval2;
            }
        }
    }

    void make_itemset(Itemset *it, Array *ary, int cnt, int *clscnt) {
        int i;
        for (i = 0; i < ary->size(); i++) {
            it->ival()->optadd((*ary)[i]);
        }
        it->set_support(cnt);
        for (i = 0; i < global::NUMCLASS; i++) it->set_cls_support(clscnt[i], i);
    }

    void get_tmpnewf_intersect(Itemset *&ljoin, Itemset *&ejoin, Itemset *&mjoin,
                               int &lcnt, int &ecnt, int &mcnt,
                               Itemset *it1, Itemset *it2, int iter) {
        int i, j, k, l;
        int nval1, nval2, diff;
        int lflge;
        Array *lary, *eary, *mary;

        num_intersect++;

        lary = interval;
        eary = interval2;
        mary = interval3;
        lary->reset();
        eary->reset();
        mary->reset();

        lcnt = ecnt = mcnt = 0;
        for (i = 0; i < global::NUMCLASS; i++) {
            ClassInfo::TMPL[i] = 0;
            ClassInfo::TMPE[i] = 0;
            ClassInfo::TMPM[i] = 0;
        }

        int dc1 = it1->support() - global::MINSUP_ABS;
        int dc2 = it2->support() - global::MINSUP_ABS;
        int df1 = 0;
        int df2 = 0;
        int icid, jcid;
        for (i = 0, j = 0; i < it1->ivalsize() && j < it2->ivalsize();) {
            if (df1 > dc1 || df2 > dc2) break;
            icid = it1->ival(i);
            jcid = it2->ival(j);
            if (icid > jcid) {
                //df must be incremented only once per customer
                while (jcid == it2->ival(j) && j < it2->ivalsize()) j += 2;
                df2++;
            } else if (icid < jcid) {
                while (icid == it1->ival(i) && i < it1->ivalsize()) i += 2;
                df1++;
            } else {
                nval1 = i;
                nval2 = j;
                while (it1->ival(i) == it1->ival(nval1) && nval1 < it1->ivalsize())
                    nval1 += 2;
                while (it2->ival(j) == it2->ival(nval2) && nval2 < it2->ivalsize())
                    nval2 += 2;

                if (ljoin && it1->ival(i + 1) + global::min_gap <= it2->ival(nval2 - 1)) {
                    lflge = 0;
                    for (k = i, l = j; k < nval1 && l < nval2;) {
                        diff = it2->ival(l + 1) - it1->ival(k + 1);
                        if (diff < global::min_gap) l += 2;
                        else if (diff > global::max_gap) k += 2;
                        else {
                            lary->optadd(icid);
                            lary->optadd(it1->ival(k + 1));
                            lflge = 1;
                            k += 2;
                        }
                    }
                    if (lflge) {
                        lcnt++;
                        ClassInfo::TMPL[ClassInfo::getcls(icid)]++;
                    }
                }

                if (ejoin) {
                    lflge = 0;
                    for (k = i, l = j; k < nval1 && l < nval2;) {
                        if (it1->ival(k + 1) < it2->ival(l + 1)) k += 2;
                        else if (it1->ival(k + 1) > it2->ival(l + 1)) l += 2;
                        else {
                            eary->optadd(icid);
                            eary->optadd(it2->ival(l + 1));
                            lflge = 1;
                            k += 2;
                            l += 2;
                        }
                    }
                    if (lflge) {
                        ecnt++;
                        ClassInfo::TMPE[ClassInfo::getcls(icid)]++;
                    }
                }

                if (mjoin && it2->ival(j + 1) + global::min_gap <= it1->ival(nval1 - 1)) {
                    lflge = 0;
                    for (k = i, l = j; k < nval1 && l < nval2;) {
                        diff = it1->ival(k + 1) - it2->ival(l + 1);
                        if (diff < global::min_gap) k += 2;
                        else if (diff > global::max_gap) l += 2;
                        else {
                            mary->optadd(icid);
                            mary->optadd(it2->ival(l + 1));
                            lflge = 1;
                            l += 2;
                        }
                    }
                    if (lflge) {
                        mcnt++;
                        ClassInfo::TMPM[ClassInfo::getcls(icid)]++;
                    }
                }
                i = nval1;
                j = nval2;
            }
        }
        if (ljoin) {
            ljoin = nullptr;
            for (i = 0; i < global::NUMCLASS; i++) {
                if (ClassInfo::TMPL[i] >= ClassInfo::MINSUP[i]) {
                    ljoin = new Itemset(iter, lary->size());
                    make_itemset(ljoin, lary, lcnt, ClassInfo::TMPL);
                    break;
                }
            }
        }
        if (ejoin) {
            ejoin = nullptr;
            for (i = 0; i < global::NUMCLASS; i++) {
                if (ClassInfo::TMPE[i] >= ClassInfo::MINSUP[i]) {
                    ejoin = new Itemset(iter, eary->size());
                    make_itemset(ejoin, eary, ecnt, ClassInfo::TMPE);
                    break;
                }
            }
        }
        if (mjoin) {
            mjoin = nullptr;
            for (i = 0; i < global::NUMCLASS; i++) {
                if (ClassInfo::TMPM[i] >= ClassInfo::MINSUP[i]) {
                    mjoin = new Itemset(iter, mary->size());
                    make_itemset(mjoin, mary, mcnt, ClassInfo::TMPM);
                    break;
                }
            }
        }
    }

    void pre_pruning(Itemset *&join, unsigned int ptempl,
                     Itemset *clas, Itemset *prefix, char use_seq) {
        double conf, conf2;
        int i, res, cit, pit;
        if (join == nullptr) return;
        pit = (*prefix)[0];
        int bitval = 0;
        int nsz = clas->size() - 2;
        if (GETBIT(global::pruning_type, Pruning_Follow - 1)) {
            for (i = 0; i <= nsz + 1 && !bitval; i++) {
                cit = (*clas)[i];
                if (use_seq) {
                    return; //TURN OFF FOR SEQUENCES

                    res = global::eqgraph[cit]->seqfind(pit);
                    if (res != -1) {
                        conf = (global::eqgraph[cit]->get_seqsup(res) * 1.0) / F1::get_sup(cit);
                        if (conf >= global::FOLLOWTHRESH) {
                            mined << "PRUNE_PRE " << pit << " -1 ";
                            clas->print_seq(SETBIT(ptempl, 1, nsz + 1));
                            global::prepruning++;
                            join = nullptr;
                            break;
                        }
                    }
                } else {
                    res = global::eqgraph[cit]->find(pit);
                    if (res != -1) {
                        conf = (global::eqgraph[cit]->get_sup(res) * 1.0) / F1::get_sup(cit);
                        conf2 = (global::eqgraph[cit]->get_sup(res) * 1.0) / F1::get_sup(pit);
                        if (conf >= global::FOLLOWTHRESH || conf2 >= global::FOLLOWTHRESH) {
                            mined << "PRUNE_PRE " << pit << " ";
                            clas->print_seq(SETBIT(ptempl, 1, nsz + 1));
                            global::prepruning++;
                            join = nullptr;
                            break;
                        }
                    }
                }

                if (nsz - i >= 0) bitval = GETBIT(ptempl, nsz - i);
            }
        }
    }

    void post_pruning(Itemset *&iset, unsigned int templ) {
        int i;
        int remsup;
        double remdb;
        if (iset == nullptr || global::NUMCLASS <= 1) return;

        if (GETBIT(global::pruning_type, Pruning_Zero - 1)) {
            for (i = 0; i < global::NUMCLASS; i++) {
                remsup = iset->support() - iset->cls_support(i);
                remdb = ClassInfo::getcnt() - ClassInfo::getcnt(i);
                if (remsup / remdb <= Pruning_Zero) {
                    mined << "PRUNE_POST ";
                    iset->print_seq(templ);
                    global::postpruning++;
                    delete iset;
                    iset = nullptr;
                    break;
                }
            }
        }
    }


    void fill_seq_template(Eqclass *EQ, Eqclass *parent, int LR) {
        if (LR == 1) {
            EQ->set_templ(SETBIT(parent->templ(), 1, EQ->templ_sz() - 1));
            EQ->set_templ2(parent->templ());
        } else if (LR == 2) {
            EQ->set_templ(SETBIT(parent->templ2(), 1, EQ->templ_sz() - 1));
            EQ->set_templ2(parent->templ2());
        }
    }

    int get_valid_el(int it, char *ibvec, char *sbvec) {
        int i, j;
        int i1, i2;
        int rval = 0;

        if (global::pruning_type == Pruning_Zero) {
            for (i = 0; i < global::eqgraph[it]->seqnum_elements(); i++) sbvec[i] = 1;
            for (i = 0; i < global::eqgraph[it]->num_elements(); i++) ibvec[i] = 1;
            rval = 1;
            return rval;
        }

        for (i = 0; i < global::eqgraph[it]->seqnum_elements(); i++) {
            sbvec[i] = 0;
        }
        for (i = 0; i < global::eqgraph[it]->num_elements(); i++) {
            ibvec[i] = 0;
        }

        for (i = 0; i < global::eqgraph[it]->seqnum_elements(); i++) {
            i1 = global::eqgraph[it]->seqget_element(i);
            for (j = i; j < global::eqgraph[it]->seqnum_elements(); j++) {
                i2 = global::eqgraph[it]->seqget_element(j);
                if (global::eqgraph[i2] && global::eqgraph[i2]->seqfind(i1) != -1) {
                    sbvec[i] = 1;
                    sbvec[j] = 1;
                    rval = 1;
                }
                if (j > i) {
                    if ((global::eqgraph[i2] && global::eqgraph[i2]->find(i1) != -1) ||
                        (global::eqgraph[i1] && global::eqgraph[i1]->seqfind(i2) != -1)) {
                        sbvec[i] = 1;
                        sbvec[j] = 1;
                        rval = 1;
                    }
                }
            }
        }


        for (i = 0; i < global::eqgraph[it]->num_elements(); i++) {
            i1 = global::eqgraph[it]->get_element(i);
            for (j = i + 1; j < global::eqgraph[it]->num_elements(); j++) {
                i2 = global::eqgraph[it]->get_element(j);
                if (global::eqgraph[i2] && global::eqgraph[i2]->find(i1) != -1) {
                    ibvec[i] = 1;
                    ibvec[j] = 1;
                    rval = 1;
                }
            }
            for (j = 0; j < global::eqgraph[it]->seqnum_elements(); j++) {
                i2 = global::eqgraph[it]->seqget_element(j);
                if (global::eqgraph[i1] && global::eqgraph[i1]->seqfind(i2) != -1) {
                    ibvec[i] = 1;
                    sbvec[j] = 1;
                    rval = 1;
                }
            }
        }

        for (i = 0; i < global::eqgraph[it]->seqnum_elements(); i++)
            if (!sbvec[i]) {
                global::L2pruning++;
                mined << "PRUNE_L2 " << it << " -1 " << global::eqgraph[it]->seqget_element(i)
                      << " " << global::eqgraph[it]->get_seqsup(i) << endl;
            }

        for (i = 0; i < global::eqgraph[it]->num_elements(); i++)
            if (!ibvec[i]) {
                global::L2pruning++;
                mined << "PRUNE_L2 " << it << " " << global::eqgraph[it]->get_element(i)
                      << " " << global::eqgraph[it]->get_sup(i) << endl;
            }
        return rval;
    }

//construct the next set of eqclasses from external disk
    Eqclass *get_ext_eqclass(int it) {
        double t1, t2;
        seconds(t1);
        int i, k, it2, supsz, supsz2;
        Itemset *ljoin = nullptr;
        Itemset *ejoin = nullptr;

        char *ibvec = nullptr, *sbvec = nullptr;
        if (!use_maxgap) {
            ibvec = sbvec = nullptr;
            if (global::eqgraph[it]->num_elements() > 0)
                ibvec = new char[global::eqgraph[it]->num_elements()];
            if (global::eqgraph[it]->seqnum_elements() > 0)
                sbvec = new char[global::eqgraph[it]->seqnum_elements()];

            if (!get_valid_el(it, ibvec, sbvec)) return nullptr;
        }

        auto *L2 = new Eqclass(1, EQCTYP1);
        if (L2 == nullptr) {
            throw runtime_error("memory exceeded : ext_class ");
        }
        //init seq pattern templates
        L2->set_templ(1);
        L2->set_templ2(0);

        interval->reset();
        interval2->reset();

        supsz = partition_get_idxsup(it);
        partition_read_item(interval->array(), it);

        int tmpit;
        for (i = 0, k = 0; i < global::eqgraph[it]->num_elements() ||
                           k < global::eqgraph[it]->seqnum_elements();) {
            ljoin = nullptr;
            ejoin = nullptr;

            it2 = global::DBASE_MAXITEM + 1;
            tmpit = global::DBASE_MAXITEM + 1;
            if (i < global::eqgraph[it]->num_elements() && (use_maxgap || ibvec[i]))
                it2 = global::eqgraph[it]->get_element(i);
            if (k < global::eqgraph[it]->seqnum_elements() && (use_maxgap || sbvec[k]))
                tmpit = global::eqgraph[it]->seqget_element(k);

            if (it2 == tmpit) {
                ejoin = (Itemset *) 1;
                ljoin = (Itemset *) 1;
                k++;
                i++;
                if (it2 == global::DBASE_MAXITEM + 1) continue;
            } else if (it2 < tmpit) {
                ejoin = (Itemset *) 1;
                i++;
            } else {
                ljoin = (Itemset *) 1;
                k++;
                it2 = tmpit;
            }
            supsz2 = partition_get_idxsup(it2);

            partition_read_item(interval2->array(), it2);

            if (ejoin) {
                ejoin = new Itemset(2, min(supsz, supsz2));
                if (ejoin == nullptr) {
                    throw runtime_error("memory exceeded");
                }
            } else ejoin = nullptr;
            if (ljoin) {
                ljoin = new Itemset(2, supsz2);
                if (ljoin == nullptr) {
                    throw runtime_error("memory exceeded");
                }
            } else ljoin = nullptr;
            get_2newf_intersect(ljoin, ejoin, interval2->array(), interval->array(),
                                supsz2, supsz);

            if (ljoin) {
                ljoin->add_item(it2);
                ljoin->add_item(it);
            }
            if (global::pruning_type > 1) post_pruning(ljoin, L2->templ());
            if (ljoin) {
                ljoin->reallocival();
                L2->prepend(ljoin);
                ljoin->print_seq(L2->templ());
            }

            if (ejoin) {
                ejoin->add_item(it2);
                ejoin->add_item(it);
            }
            if (global::pruning_type > 1) post_pruning(ejoin, L2->templ2());
            if (ejoin) {
                ejoin->reallocival();
                L2->prepend2(ejoin);
                ejoin->print_seq(L2->templ2());
            }
        }

        seconds(t2);
        L2ISECTTIME += t2 - t1;
        return L2;
    }

    void fill_join(Itemset *join, Itemset *hdr1, Itemset *hdr2) {
        int i;

        join->add_item((*hdr2)[0]);
        for (i = 0; i < hdr1->size(); i++) {
            join->add_item((*hdr1)[i]);
        }
    }

    Itemset *prune_decision(Itemset *it1, Itemset *it2, unsigned int ptempl, int jflg) {
        int i, j, k;

        //prune if seq pat exceeds the max seq len or iset len
        int bit, seqlen = 0, isetlen = 0, maxisetlen = 0;
        for (i = 0; i < it2->size(); i++) {
            bit = GETBIT(ptempl, i);
            if (bit) {
                seqlen++;
                if (maxisetlen < isetlen) maxisetlen = isetlen;
                isetlen = 0;
            } else isetlen++;
        }
        if (maxisetlen < isetlen) maxisetlen = isetlen;
        seqlen++;
        maxisetlen++;

        if (seqlen > global::max_seq_len) return nullptr;
        if (maxisetlen > global::max_iset_len) return nullptr;


        //global::max_gap destroys the downward closure property, so we cannot prune
        if (use_maxgap) return (Itemset *) 1;

        int l1 = (*it1)[0];
        int l2 = (*it2)[0];
        int nsz;
        if (cspade_args.use_hash && (it2->size() > 2)) {
            if (cspade_args.recursive) return (Itemset *) 1;

            unsigned int ttpl;
            FreqIt fit(it2->size(), 0);

            //skip the last two subsets (or omit the first two elements)
            nsz = it2->size() - 2;

            for (i = nsz + 1; i >= 1; i--) {
                k = 0;
                ttpl = 0;
                //form new subset template
                if (i == nsz + 1) ttpl = (ptempl >> 1);
                else {
                    for (j = 0; j < i; j++) {
                        bit = GETBIT(ptempl, nsz - j + 1);
                        ttpl = SETBIT(ttpl, bit, nsz - j);
                    }
                    bit = GETBIT(ptempl, nsz - j + 1);
                    bit = bit || GETBIT(ptempl, nsz - j);
                    ttpl = SETBIT(ttpl, bit, nsz - j);
                    j += 2;
                    for (; j < nsz + 2; j++) {
                        bit = GETBIT(ptempl, nsz - j + 1);
                        ttpl = SETBIT(ttpl, bit, nsz - j + 1);
                    }
                }
                //form new subset by omitting the i-th item
                fit.seq[k++] = l1;
                fit.seq[k++] = l2;
                for (j = 1; j < nsz + 2; j++) {
                    if (j != i) {
                        fit.seq[k++] = (*it2)[j];
                    }
                }
                fit.templ = ttpl;

                //???? Does this work for suffix classes
                if (fit.seq[fit.size() - 1] == (*it1)[it1->size() - 1] && !cspade_args.recursive) {
                    //elements should be in current class
                    if (FreqArrayPos > 0) {
                        if (!EqGrNode::bsearch(0, FreqArrayPos - 1, FreqArray, fit, cspade_args.recursive)) {
                            return nullptr;
                        }
                    } else return nullptr;
                } else if (fit.seq[fit.size() - 1] > (*it1)[it1->size() - 1]) {
                    // class must already have been processed, otherwise we can't prune
                    if (!global::eqgraph[fit.seq[fit.size() - 1]]->find_freqarray(fit, cspade_args.recursive)) {
                        return nullptr;
                    }
                }
            }
        } else {// if (it1->size() == 2){
            bit = 0;
            nsz = it2->size() - 2;
            for (i = 0; i <= nsz + 1 && !bit; i++) {
                l2 = (*it2)[i];
                if (global::eqgraph[l2]) {
                    if (jflg == LJOIN || jflg == MJOIN) {
                        if (global::eqgraph[l2]->seqfind(l1) == -1)
                            return nullptr;
                    } else {
                        if (global::eqgraph[l2]->find(l1) == -1)
                            return nullptr;
                    }
                } else return nullptr;
                if (nsz - i >= 0) bit = GETBIT(ptempl, nsz - i);
            }
        }
        return (Itemset *) 1;
    }


    void insert_freqarray(Lists<Eqclass *> *LargeL) {
        //insert frequent itemsets into hash table
        ListNodes<Eqclass *> *chd;
        ListNodes<Itemset *> *hdr1, *hdr2;
        Eqclass *cluster;

        chd = LargeL->head();
        for (; chd; chd = chd->next()) {
            cluster = chd->item();
            hdr1 = cluster->list()->head();
            for (; hdr1; hdr1 = hdr1->next()) {
                add_freq(hdr1->item(), cluster->templ());
            }
            hdr2 = cluster->list2()->head();
            for (; hdr2; hdr2 = hdr2->next()) {
                add_freq(hdr2->item(), cluster->templ2());
            }
        }
    }

    void
    process_cluster_list1(ListNodes<Itemset *> *hdr1, Lists<Itemset *> *cluster2, Lists<Eqclass *> *LargeL, int iter,
                          int eqtype, Eqclass *parent) {
        ListNodes<Itemset *> *hdr2;
        auto *EQ = new Eqclass(iter - 1, eqtype);
        if (EQ == nullptr) {
            throw runtime_error("memory exceeded");
        }
        fill_seq_template(EQ, parent, 2);
        Itemset *ljoin, *ejoin, *mjoin;
        int lsup, esup, msup;

        hdr2 = cluster2->head();
        for (; hdr2; hdr2 = hdr2->next()) {
            ljoin = prune_decision(hdr2->item(), hdr1->item(), EQ->templ(), LJOIN);
            ejoin = nullptr;
            mjoin = nullptr;
            lsup = esup = msup = 0;
            if (global::pruning_type > 1)
                pre_pruning(ljoin, EQ->templ(), hdr1->item(), hdr2->item(), 1);
            if (ljoin)
                get_tmpnewf_intersect(ljoin, ejoin, mjoin, lsup, esup, msup,
                                      hdr2->item(), hdr1->item(), iter);
            if (ljoin) fill_join(ljoin, hdr1->item(), hdr2->item());
            if (global::pruning_type > 1) post_pruning(ljoin, EQ->templ());
            if (ljoin) {
                global::NumLargeItemset[iter - 1]++;
                ljoin->print_seq(EQ->templ());
                EQ->append(ljoin);
            }
        }

        hdr2 = hdr1->next();
        for (; hdr2 != nullptr; hdr2 = hdr2->next()) {
            ejoin = prune_decision(hdr2->item(), hdr1->item(), EQ->templ2(), EJOIN);
            ljoin = nullptr;
            mjoin = nullptr;
            lsup = esup = msup = 0;
            if (global::pruning_type > 1)
                pre_pruning(ejoin, EQ->templ2(), hdr1->item(), hdr2->item(), 0);
            if (ejoin)
                get_tmpnewf_intersect(ljoin, ejoin, mjoin, lsup, esup, msup,
                                      hdr2->item(), hdr1->item(), iter);
            if (ejoin) fill_join(ejoin, hdr1->item(), hdr2->item());
            if (global::pruning_type > 1) post_pruning(ejoin, EQ->templ2());
            if (ejoin) {
                global::NumLargeItemset[iter - 1]++;
                ejoin->print_seq(EQ->templ2());
                EQ->append2(ejoin);
            }
        }

        if (EQ) {
            if ((EQ->list()->size() > 0) || (EQ->list2()->size() > 0)) {
                if (cspade_args.recursive) {
                    process_cluster1(EQ, nullptr, iter + 1);
                    delete EQ;
                } else LargeL->append(EQ);
            } else {
                delete EQ;
                EQ = nullptr;
            }
        }
    }

    void process_cluster_list2(ListNodes<Itemset *> *hdr1, int i, Eqclass **EQ, Lists<Eqclass *> *LargeL, int iter) {
        int j;

        ListNodes<Itemset *> *hdr2;
        Itemset *ljoin, *ejoin, *mjoin;
        int lsup, esup, msup;

        //join with sequences
        hdr2 = hdr1;
        for (j = i; hdr2; j++, hdr2 = hdr2->next()) {
            ljoin = prune_decision(hdr1->item(), hdr2->item(), EQ[j]->templ(), LJOIN);
            if (hdr2 == hdr1) {
                ejoin = mjoin = nullptr;
            } else {
                ejoin = prune_decision(hdr2->item(), hdr1->item(), EQ[i]->templ2(), EJOIN);
                mjoin = prune_decision(hdr2->item(), hdr1->item(), EQ[i]->templ(), MJOIN);
            }
            lsup = esup = msup = 0;
            if (global::pruning_type > 1) {
                pre_pruning(ejoin, EQ[i]->templ2(), hdr1->item(), hdr2->item(), 0);
                pre_pruning(ljoin, EQ[j]->templ(), hdr2->item(), hdr1->item(), 1);
                pre_pruning(mjoin, EQ[i]->templ(), hdr1->item(), hdr2->item(), 1);
            }

            if (ljoin || ejoin || mjoin)
                get_tmpnewf_intersect(ljoin, ejoin, mjoin, lsup, esup, msup,
                                      hdr1->item(), hdr2->item(), iter);
            if (ljoin) fill_join(ljoin, hdr2->item(), hdr1->item());
            if (global::pruning_type > 1) post_pruning(ljoin, EQ[j]->templ());
            if (ljoin) {
                global::NumLargeItemset[iter - 1]++;
                ljoin->print_seq(EQ[j]->templ());
                EQ[j]->append(ljoin);
            }

            if (ejoin) fill_join(ejoin, hdr1->item(), hdr2->item());
            if (global::pruning_type > 1) post_pruning(ejoin, EQ[i]->templ2());
            if (ejoin) {
                global::NumLargeItemset[iter - 1]++;
                ejoin->print_seq(EQ[i]->templ2());
                EQ[i]->append2(ejoin);
            }

            if (mjoin) fill_join(mjoin, hdr1->item(), hdr2->item());
            if (global::pruning_type > 1) post_pruning(mjoin, EQ[i]->templ());
            if (mjoin) {
                global::NumLargeItemset[iter - 1]++;
                mjoin->print_seq(EQ[i]->templ());
                EQ[i]->append(mjoin);
            }
        }
        if ((EQ[i]->list()->size() > 0) || (EQ[i]->list2()->size() > 0)) {
            if (cspade_args.recursive) {
                process_cluster1(EQ[i], nullptr, iter + 1);
                delete EQ[i];
                EQ[i] = nullptr;
            } else LargeL->append(EQ[i]);
        } else {
            delete EQ[i];
            EQ[i] = nullptr;
        }

    }


    void process_cluster1(Eqclass *cluster, Lists<Eqclass *> *LargeL, int iter) {
        Eqclass **EQ = nullptr;
        ListNodes<Itemset *> *hdr1, *hdr2;
        int i;

        if (cluster->list()->head()) {
            EQ = new Eqclass *[cluster->list()->size()];
            if (EQ == nullptr) {
                throw runtime_error("memory exceeded");
            }
            for (i = 0; i < cluster->list()->size(); i++) {
                EQ[i] = new Eqclass(iter - 1, EQCTYP1);
                if (EQ[i] == nullptr) {
                    throw runtime_error("memory exceeded");
                }
                fill_seq_template(EQ[i], cluster, 1);
            }
        }

        hdr1 = cluster->list()->head();
        for (i = 0; hdr1; hdr1 = hdr1->next(), i++) {
            process_cluster_list2(hdr1, i, EQ, LargeL, iter);
        }

        delete[] EQ;

        hdr2 = cluster->list2()->head();
        for (; hdr2; hdr2 = hdr2->next()) {
            process_cluster_list1(hdr2, cluster->list(), LargeL, iter, EQCTYP1, cluster);
        }

        if (global::maxiter < iter) global::maxiter = iter;

    }


    void find_large(Eqclass *cluster, int it) {
        Lists<Eqclass *> *LargeL, *Candidate;
        ListNodes<Eqclass *> *chd;
        int iter;
        int LargelistSum = 0;
        int more;

        more = 1;
        Candidate = new Lists<Eqclass *>;
        Candidate->append(cluster);
        for (iter = 3; more; iter++) {
            LargeL = new Lists<Eqclass *>;
            chd = Candidate->head();
            for (; chd; chd = chd->next()) {
                process_cluster1(chd->item(), LargeL, iter);
                //reclaim memory for this class immediately
                delete chd->item();
                chd->set_item(nullptr);
            }
            Candidate->clear();
            delete Candidate;

            if (cspade_args.use_hash) insert_freqarray(LargeL);
            chd = LargeL->head();
            LargelistSum = 0;
            for (; chd; chd = chd->next()) {
                LargelistSum += chd->item()->list()->size();
                if (chd->item()->list2())
                    LargelistSum += chd->item()->list2()->size();
            }
            more = (LargelistSum > 0);

            Candidate = LargeL;
            memlog << it << " " << global::MEMUSED << endl;

            if (!more) {
                LargeL->clear();
                delete LargeL;
            }
        }
    }


    void process_class(int it) {

        //from 2-itemsets from ext disk
        Eqclass *large2it = get_ext_eqclass(it);
        if (large2it == nullptr) return;

        mined << "PROCESS " << it << endl;

        memlog << it << " " << global::MEMUSED << endl;
        if (use_maxgap) {
            process_maxgap(large2it);
        } else {
            if (cspade_args.recursive) {
                process_cluster1(large2it, nullptr, 3);
                delete large2it;
            } else find_large(large2it, it);
        }
    }

    void newSeq() {
        int i, j;

        if (cspade_args.use_hash)
            FreqArray = (FreqIt **) malloc(FreqArraySz * sizeof(FreqIt *));
        //form large itemsets for each eqclass
        if (cspade_args.use_ascending != -2) {
            if (cspade_args.use_ascending == -1) {
                for (i = 0; i < global::DBASE_MAXITEM; i++)
                    if (global::eqgraph[i]) {
                        memlog << i << " " << global::MEMUSED << endl;
                        process_class(i);
                        memlog << i << " " << global::MEMUSED << endl;
                    }
            } else if (global::eqgraph[cspade_args.use_ascending])
                process_class(cspade_args.use_ascending);
        } else {
            for (i = global::DBASE_MAXITEM - 1; i >= 0; i--) {
                if (global::eqgraph[i]) {
                    memlog << i << " " << global::MEMUSED << endl;
                    if (cspade_args.use_hash) FreqArrayPos = 0;
                    process_class(i);
                    if (cspade_args.use_hash) {
                        if (FreqArrayPos > 0) {
                            auto **fit = new FreqIt *[FreqArrayPos];
                            for (j = 0; j < FreqArrayPos; j++) {
                                fit[j] = FreqArray[j];
                            }
                            global::eqgraph[i]->set_freqarray(fit, FreqArrayPos);
                        }
                    }
                    memlog << i << " " << global::MEMUSED << endl;
                }
            }
        }
    }


    void read_files() {
        int i;

        global::NumLargeItemset = new int[(int) global::DBASE_AVG_TRANS_SZ * 30];
        bzero((char *) global::NumLargeItemset, sizeof(int) * ((int) global::DBASE_AVG_TRANS_SZ * 30));

        global::eqgraph = new EqGrNode *[global::DBASE_MAXITEM];
        bzero((char *) global::eqgraph, global::DBASE_MAXITEM * sizeof(EqGrNode *));

        double t1, t2;
        seconds(t1);
        global::NumLargeItemset[0] = make_l1_pass();
        seconds(t2);
        EXTL1TIME = t2 - t1;

        if (!cspade_args.do_l2) {
            global::NumLargeItemset[1] = make_l2_pass();
            seconds(t1);
            EXTL2TIME = t1 - t2;
        } else {
            global::NumLargeItemset[1] = get_file_l2(cspade_args.it2f, cspade_args.seqf);
            seconds(t1);
            EXTL2TIME = t1 - t2;
        }

        for (i = 0; i < global::DBASE_MAXITEM; i++) {
            if (global::eqgraph[i]) {
                if (global::eqgraph[i]->num_elements() > 0)
                    global::eqgraph[i]->elements()->compact();
                if (global::eqgraph[i]->seqnum_elements() > 0)
                    global::eqgraph[i]->seqelements()->compact();
            }
        }

        maxitemsup = 0;
        int sup;
        for (i = 0; i < global::DBASE_MAXITEM; i++) {
            sup = partition_get_idxsup(i);
            if (maxitemsup < sup) maxitemsup = sup;
        }
        interval = new Array(maxitemsup);
        interval2 = new Array(maxitemsup);
        interval3 = new Array(maxitemsup);
    }

    void mine() {
        int i;
        double ts, te;
        double t1, t2;
        seconds(ts);

        partition_alloc(cspade_args.dataf, cspade_args.idxf);
        ClassInfo cls(cspade_args.use_class, cspade_args.classf);
        read_files();


        if (use_maxgap)
            IBM = new ItBufMgr(global::NumLargeItemset[0]);
        seconds(t1);
        newSeq();
        seconds(t2);
        double FKtime = t2 - t1;
        if (use_maxgap)
            delete IBM;
        seconds(te);

        summary << "SPADE ";
        if (cspade_args.use_hash)
            summary << "USEHASH ";
        summary << cspade_args.dataf << ' ' << global::MINSUP_PER << ' ' << global::MINSUP_ABS << ' '
                << num_intersect << ' ' << L2ISECTTIME
                << ' '
                << global::pruning_type << ' ' << global::L2pruning << ' ' << global::prepruning << ' '
                << global::postpruning;
        summary << "0 ";
        if (use_maxgap)
            summary << global::max_gap;
        else
            summary << "-1 ";
        summary << global::min_gap << ' ' << global::max_iset_len << ' ' << global::max_seq_len << " :";
        for (i = 0; i < global::maxiter; i++) {
            summary << global::NumLargeItemset[i] << ' ';
        }
        summary << ": " << EXTL1TIME << ' ' << EXTL2TIME << ' ' << FKtime << ' ' << te - ts << endl;

        partition_dealloc();

        delete interval;
        delete interval2;
        delete interval3;
        for (i = 0; i < global::DBASE_MAXITEM; i++) {
            if (global::eqgraph[i]) delete global::eqgraph[i];
        }
        delete[] global::eqgraph;

        memlog << global::MEMUSED << endl;
    }

    void print_args() {
        cout << "num_partitions = " << cspade_args.num_partitions << endl;
        cout << "min_support = " << cspade_args.min_support << endl;
        cout << "use_ascending = " << cspade_args.use_ascending << endl;
        cout << "use_class = " << cspade_args.use_class << endl;
        cout << "do_l2 = " << cspade_args.do_l2 << endl;
        cout << "use_hash = " << cspade_args.use_hash << endl;
        cout << "min_gap = " << cspade_args.min_gap << endl;
        cout << "maxmem = " << cspade_args.maxmem << endl;
        cout << "recursive = " << cspade_args.recursive << endl;
        cout << "pruning_type = " << cspade_args.pruning_type << endl;
        cout << "max_gap = " << cspade_args.max_gap << endl;
        cout << "use_window = " << cspade_args.num_partitions << endl;
        cout << "max_seq_len = " << cspade_args.max_seq_len << endl;
        cout << "max_iset_len = " << cspade_args.max_iset_len << endl;
        cout << "twoseq = " << cspade_args.twoseq << endl;
        cout << "use_diff = " << cspade_args.use_diff << endl;
        cout << "use_newformat = " << cspade_args.use_newformat << endl;
        cout << "no_minus_off = " << cspade_args.no_minus_off << endl;


        cout << "global::num_partitions =" << global::num_partitions << endl;
        cout << "global::AVAILMEM =" << global::AVAILMEM << endl;
        cout << "global::MINSUP_PER =" << global::MINSUP_PER << endl;
        cout << "global::MINSUP_ABS = " << global::MINSUP_ABS << endl;
        cout << "global::max_gap =" << global::max_gap << endl;
        cout << "global::min_gap =" << global::min_gap << endl;
        cout << "global::max_seq_len =" << global::max_seq_len << endl;
        cout << "global::max_iset_len =" << global::max_iset_len << endl;
        cout << "global::pruning_type =" << global::pruning_type << endl;
        cout << "global::max_seq_len =" << global::max_seq_len << endl;
        cout << "global::max_iset_len =" << global::max_iset_len << endl;
    }

}


result_t cspade(const string &asciifile, sequence::arg_t &_args) {
    sequence::cspade_args = _args;
    result_t result;

    try {
        convert_bin(asciifile);
        create_conf(false);
        sequence::populate_global();
        exttpose();
        sequence::mine();
        result.success = true;
        result.mined = mined.str();
        result.logger = logger.str();
        result.summary = summary.str();
        result.memlog = memlog.str();
        result.nsequences = global::DBASE_NUM_TRANS;
    }
    catch (std::exception &error) {
        result.success = false;
        result.error = error.what();
    }

    // Clear the stringstream otherwise the next run will get duplicated result
    mined.str(std::string());
    mined.clear();
    logger.str(std::string());
    logger.clear();
    summary.str(std::string());
    summary.clear();
    memlog.str(std::string());
    memlog.clear();
    return result;
}

int main(int argc, char **argv) {
    sequence::arg_t _args;
    _args.name = "testdata/bb";
    _args.num_partitions = 1;
    _args.min_support = 2;
    _args.do_l2 = true;
    _args.recursive = true;
    _args.max_iset_len = 3;
    _args.max_seq_len = 4;

    _args = sequence::populate_names(_args);

    auto outcome = cspade("testdata/zaki.txt", _args);
    cout << outcome.mined;
//    cout << outcome.logger;
//    cout << outcome.summary;
//    cout << outcome.memlog;
//    cout << outcome.nsequences;
}