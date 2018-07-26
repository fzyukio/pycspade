#include <cerrno>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/time.h>
#include <cmath>

#include "utils.h"
#include "calcdb.h"
#include "ArrayT.h"
#include "exttpose.h"

void sort_get_l2(int &l2cnt, int fd, ofstream &ofd, int *backidx, int *freqidx,
                 int numfreq, int *offsets, unsigned char *cntary, char use_seq, int minsupport, bool twoseq) {
    //write 2-itemsets counts to file

    int i, j, k, fcnt;
    int lit;
    long sortflen;
    int *sortary;
    int itbuf[3];

    sortflen = lseek(fd, 0, SEEK_CUR);
    if (sortflen < 0) {
        throw runtime_error("SEEK SEQ");
    }
    //logger << "SORT " << sortflen << endl;
    if (sortflen > 0) {
#ifdef DEC
        sortary = (int *) mmap((char *)nullptr, sortflen,
                               (PROT_READ|PROT_WRITE),
                               (MAP_FILE|MAP_VARIABLE|MAP_PRIVATE), fd, 0);
#else
        sortary = (int *) mmap((char *) nullptr, sortflen,
                               (PROT_READ | PROT_WRITE),
                               MAP_PRIVATE, fd, 0);
#endif
        if (sortary == (int *) -1) {
            throw runtime_error("SEQFd MMAP ERROR");
        }

        qsort(sortary, (sortflen / sizeof(int)) / 2, 2 * sizeof(int), cmp2it);
    }

    int numel = sortflen / sizeof(int);
    i = 0;
    fcnt = 0;
    for (j = 0; j < numfreq; j++) {
        if (use_seq) k = 0;
        else k = j + 1;
        for (; k < numfreq; k++) {
            fcnt = 0;
            if (sortflen > 0 && i < numel) {
                while (i < numel &&
                       j == freqidx[sortary[i]] && k == freqidx[sortary[i + 1]]) {
                    fcnt += 256;
                    i += 2;
                }
            }
            if (use_seq) fcnt += (int) cntary[j * numfreq + k];
            else {
                lit = j;
                lit = (offsets[lit] - lit - 1);
                fcnt += (int) cntary[lit + k];
            }

            if (fcnt >= minsupport) {
                if (twoseq) {
                    ofd.write((char *) &fcnt, ITSZ);
                } else {
                    itbuf[0] = backidx[j];
                    itbuf[1] = backidx[k];
                    itbuf[2] = fcnt;
                    ofd.write((char *) itbuf, 3 * ITSZ);
                }
                //logger << backidx[j] << ((use_seq)?" -> ":" ")
                //     << backidx[k] << " SUPP " << fcnt << endl;
                l2cnt++;
            }
        }
    }
    if (sortflen > 0) munmap((caddr_t) sortary, sortflen);
}


void process_cust(int *fidx, int fcnt, int numfreq, int *backidx,
                  ArrayT **extary, unsigned char *seq2, unsigned char *itcnt2, char *ocust,
                  int *offsets, int seqfd, int isetfd, bool use_seq) {
    int j, k, lit;
    int ii1, ii2;

    for (k = 0; k < fcnt; k++) {
        for (j = k; j < fcnt; j++) {
            if (use_seq && extary[fidx[j]]->size() > 0) {
                lit = extary[fidx[j]]->size() - 1;
                if ((*extary[fidx[k]])[0] < (*extary[fidx[j]])[lit]) {
                    if ((++seq2[fidx[k] * numfreq + fidx[j]]) == 0) {
                        write(seqfd, (char *) &backidx[fidx[k]], ITSZ);
                        write(seqfd, (char *) &backidx[fidx[j]], ITSZ);
                    }
                }
            }
            if (j > k) {
                if (fidx[k] < fidx[j]) {
                    ii1 = fidx[k];
                    ii2 = fidx[j];
                } else {
                    ii2 = fidx[k];
                    ii1 = fidx[j];
                }
                lit = offsets[ii1] - ii1 - 1;
                if (ocust[lit + ii2] == 1) {
                    if ((++itcnt2[lit + ii2]) == 0) {
                        write(isetfd, (char *) &backidx[ii1], ITSZ);
                        write(isetfd, (char *) &backidx[ii2], ITSZ);
                        //itcnt2[lit+ii2] = 0;
                    }
                    ocust[lit + ii2] = 0;
                }

                if (extary[fidx[k]]->size() > 0) {
                    lit = extary[fidx[k]]->size() - 1;
                    if ((*extary[fidx[j]])[0] < (*extary[fidx[k]])[lit]) {
                        if ((++seq2[fidx[j] * numfreq + fidx[k]]) == 0) {
                            write(seqfd, (char *) &backidx[fidx[j]], ITSZ);
                            write(seqfd, (char *) &backidx[fidx[k]], ITSZ);
                        }
                    }
                }
            }
        }
        extary[fidx[k]]->reset();
    }
}

void do_invert_db(CalcDb *DCB, int pblk, ArrayT **extary, int numfreq, int *freqidx, int *backidx, int *fidx,
                  int mincustid, int maxcustid, bool use_seq) {
    using sequence::cspade_args;

    int numitem, tid, custid;
    int *buf;
    char tmpnam[300];
    int i, j, k;
    int fd;
    int idx;

    DCB->get_first_blk();
    DCB->get_next_trans(buf, numitem, tid, custid);
    int ocid;// = -1;
    for (int p = 0; p < cspade_args.num_partitions; p++) {
        if (cspade_args.num_partitions > 1) sprintf(tmpnam, "%s.P%d", cspade_args.dataf.c_str(), p);
        else sprintf(tmpnam, "%s", cspade_args.dataf.c_str());

        // Creating tpose file
        if ((fd = open(tmpnam, (O_WRONLY | O_CREAT | O_TRUNC), 0666)) < 0) {
            string error_message = "can't open tpose file: " + string(tmpnam);
            throw runtime_error(error_message);
        }

        for (i = 0; i < numfreq; i++) {
            extary[i]->reset();
        }
        //count 2-itemsets
        int plb = p * pblk + mincustid;
        int pub = plb + pblk;
        if (pub >= maxcustid) pub = maxcustid + 1;
        logger << "BOUNDS " << plb << " " << pub << endl;
        int fcnt;
        for (; !DCB->eof() && custid < pub;) {
            fcnt = 0;
            ocid = custid;
            while (!DCB->eof() && ocid == custid && custid < pub) {
                if (cspade_args.use_diff) {
                    //add this tid to all items not in the trans
                    k = 0;
                    for (j = 0; j < numitem; j++) {
                        if (freqidx[buf[j]] == -1) continue;

                        while (backidx[k] < buf[j]) {
                            idx = k;
                            if (!cspade_args.use_newformat)
                                extary[idx]->add(fd, tid, use_seq, p);
                            else extary[idx]->add(fd, tid, use_seq, p, custid);
                            //}
                            k++;
                        }
                        k++; //skip over buf[j]
                    }
                    for (; k < numfreq; k++) {
                        idx = k;
                        if (!cspade_args.use_newformat)
                            extary[idx]->add(fd, tid, use_seq, p);
                        else extary[idx]->add(fd, tid, use_seq, p, custid);
                        //}
                    }
                } else {
                    // add this tid to all items in the trans
                    for (j = 0; j < numitem; j++) {
                        idx = freqidx[buf[j]];
                        if (idx != -1) {
                            if (!cspade_args.use_newformat) {
                                if (use_seq && extary[idx]->flg() == 0) {
                                    fidx[fcnt] = idx;
                                    fcnt++;
                                    extary[idx]->setflg(1);
                                    extary[idx]->add(fd, tid, use_seq, p, custid);
                                } else {
                                    extary[idx]->add(fd, tid, use_seq, p);
                                }
                            } else {
                                extary[idx]->add(fd, tid, use_seq, p, custid);
                            }
                        }
                    }
                }

                DCB->get_next_trans(buf, numitem, tid, custid);
            }
            if (!cspade_args.use_newformat && use_seq) {
                for (k = 0; k < fcnt; k++) {
                    extary[fidx[k]]->setlastpos();
                    extary[fidx[k]]->setflg(0);
                }
                fcnt = 0;
            }
        }

        for (i = 0; i < numfreq; i++) {
            extary[i]->flushbuf(fd, use_seq, p);
        }
        close(fd);
    }
    logger << "WROTE INVERT " << endl;
}

void tpose(bool use_seq) {
    using sequence::cspade_args;

    int i, j, l;
    int idx;
    int custid, tid, numitem, fcnt;
    ofstream ofd;
    double t1, t2;
    int sumsup = 0, sumdiff = 0;

    auto *DCB = new CalcDb(cspade_args.binf);

    auto *itcnt = new int[global::DBASE_MAXITEM];
    auto *ocnt = new int[global::DBASE_MAXITEM];
    auto *itlen = new int[global::DBASE_MAXITEM];
    bzero((char *) itcnt, ((global::DBASE_MAXITEM) * ITSZ));
    for (i = 0; i < global::DBASE_MAXITEM; i++) ocnt[i] = -1;
    bzero((char *) itlen, ((global::DBASE_MAXITEM) * ITSZ));

    //count 1 items
    int *buf;
    DCB->get_first_blk();
    DCB->get_next_trans(buf, numitem, tid, custid);
    int mincustid = custid;
    while (!DCB->eof()) {
        //logger << custid << " " << tid << " " << numitem;
        for (j = 0; j < numitem; j++) {
            //logger << " " << buf[j] << flush;
            itlen[buf[j]]++;
            if (use_seq && ocnt[buf[j]] != custid) {
                itcnt[buf[j]]++;
                ocnt[buf[j]] = custid;
            }
            //if (buf[j] == 17) logger << " " << tid;
        }
        //logger << endl;
        DCB->get_next_trans(buf, numitem, tid, custid);
    }
    //logger << endl;
    int maxcustid = custid;
    logger << "MINMAX " << mincustid << " " << maxcustid << endl;

    int *freqidx = new int[global::DBASE_MAXITEM];
    int numfreq = 0;
    for (i = 0; i < global::DBASE_MAXITEM; i++) {
        if (use_seq) {
            if (itcnt[i] >= global::MINSUP_ABS) {
                logger << i << " SUPP " << itcnt[i] << endl;
                freqidx[i] = numfreq;
                numfreq++;
            } else freqidx[i] = -1;
        } else {
            if (itlen[i] >= global::MINSUP_ABS) {
                freqidx[i] = numfreq;
                numfreq++;
                sumsup += itlen[i];
                sumdiff += (global::DBASE_NUM_TRANS - itlen[i]);
            } else freqidx[i] = -1;
        }
        //if (i == 17) logger << " 17 SUP " << itlen[17] << endl;
    }
    int *backidx = new int[numfreq];
    numfreq = 0;
    for (i = 0; i < global::DBASE_MAXITEM; i++) {
        if (use_seq) {
            if (itcnt[i] >= global::MINSUP_ABS)
                backidx[numfreq++] = i;
        } else {
            if (itlen[i] >= global::MINSUP_ABS)
                backidx[numfreq++] = i;
        }
    }

    logger << "numfreq " << numfreq << " :  " << " SUMSUP SUMDIFF = " << sumsup << " " << sumdiff << endl;

    if (numfreq == 0)
        return;

    int extarysz = global::AVAILMEM / numfreq;
    extarysz /= sizeof(int);
    logger << "EXTRARYSZ " << extarysz << endl;
    if (extarysz < 2) extarysz = 2;
    ArrayT **extary = new ArrayT *[numfreq];
    for (i = 0; i < numfreq; i++) {
        extary[i] = new ArrayT(extarysz, cspade_args.num_partitions);
    }

    char tmpnam[300];
    int plb, pub, pblk;
    pblk = (int) ceil(((double) (maxcustid - mincustid + 1)) / cspade_args.num_partitions);
    if (cspade_args.num_partitions > 1) {
        DCB->get_first_blk();
        DCB->get_next_trans(buf, numitem, tid, custid);
    }
    for (j = 0; j < cspade_args.num_partitions; j++) {
        //construct offsets for 1-itemsets
        if (cspade_args.num_partitions > 1) {
            sprintf(tmpnam, "%s.P%d", cspade_args.idxf.c_str(), j);
            plb = j * pblk + mincustid;
            pub = plb + pblk;
            if (pub > maxcustid) pub = maxcustid + 1;
            bzero((char *) itcnt, ((global::DBASE_MAXITEM) * ITSZ));
            for (i = 0; i < global::DBASE_MAXITEM; i++) ocnt[i] = -1;
            bzero((char *) itlen, ((global::DBASE_MAXITEM) * ITSZ));
            for (; !DCB->eof() && custid < pub;) {
                for (i = 0; i < numitem; i++) {
                    itlen[buf[i]]++;
                    if (use_seq && ocnt[buf[i]] != custid) {
                        itcnt[buf[i]]++;
                        ocnt[buf[i]] = custid;
                    }
                }
                DCB->get_next_trans(buf, numitem, tid, custid);
            }
        } else sprintf(tmpnam, "%s", cspade_args.idxf.c_str());
        logger << "OPENED " << tmpnam << endl;
        ofd.open(tmpnam);
        if (!ofd) {
            string error_message = "can't open idx file: " + string(tmpnam);
            throw runtime_error(error_message);
        }

        int file_offset = 0;
        int null = -1;
        for (i = 0; i < global::DBASE_MAXITEM; i++) {
            if (freqidx[i] != -1) {
                ofd.write((char *) &file_offset, ITSZ);
                extary[freqidx[i]]->set_offset(file_offset, j);
                if (use_seq) {
                    if (cspade_args.use_newformat) file_offset += (2 * itlen[i]);
                    else file_offset += (2 * itcnt[i] + itlen[i]);
                } else {
                    if (cspade_args.use_diff) file_offset += (global::DBASE_NUM_TRANS - itlen[i]);
                    else file_offset += itlen[i];
                }
            } else if (cspade_args.no_minus_off) {
                ofd.write((char *) &file_offset, ITSZ);
            } else ofd.write((char *) &null, ITSZ);
        }
        logger << "OFF " << i << " " << file_offset << endl;
        ofd.write((char *) &file_offset, ITSZ);
        ofd.close();
    }

    delete[] ocnt;
    delete[] itlen;
    delete[] itcnt;

    logger << "Wrote Offt " << endl;

    int *fidx = new int[numfreq];

    int ocid = -1;
    if (cspade_args.do_l2) {
        int seqfd, isetfd;
        if (use_seq) {
            if ((seqfd = open("tmpseq", (O_RDWR | O_CREAT | O_TRUNC), 0666)) < 0) {
                string error_message = "can't open tmpseq file: ";
                throw runtime_error(error_message);
            }
        }

        if ((isetfd = open("tmpiset", (O_RDWR | O_CREAT | O_TRUNC), 0666)) < 0) {
            string error_message = "can't open tmpseq file: ";
            throw runtime_error(error_message);
        }

        unsigned char *seq2;
        if (use_seq) {
            seq2 = new unsigned char[numfreq * numfreq];
            bzero((char *) seq2, numfreq * numfreq * sizeof(unsigned char));
        }

        auto *itcnt2 = new unsigned char[(numfreq * (numfreq - 1) / 2)];
        bzero((char *) itcnt2, (numfreq * (numfreq - 1) / 2) * sizeof(unsigned char));
        //for (i=0; i < numfreq*(numfreq-1)/2; i++) itcnt2[i] = 0;
        auto *ocust = new char[(numfreq * (numfreq - 1) / 2)];
        bzero((char *) ocust, (numfreq * (numfreq - 1) / 2) * sizeof(char));
        auto *offsets = new int[numfreq];
        int offt = 0;
        for (i = numfreq - 1; i >= 0; i--) {
            offsets[numfreq - i - 1] = offt;
            offt += i;
        }

        ocid = -1;
        int lit;
        //count 2-itemsets
        DCB->get_first_blk();
        DCB->get_next_trans(buf, numitem, tid, custid);
        while (!DCB->eof()) {
            fcnt = 0;
            ocid = custid;
            while (!DCB->eof() && ocid == custid) {
                for (j = 0; j < numitem; j++) {
                    idx = freqidx[buf[j]];
                    if (idx != -1) {
                        if (use_seq) {
                            if (extary[idx]->size() == 0) {
                                fidx[fcnt] = idx;
                                fcnt++;
                                extary[idx]->setitem(0, tid);
                                extary[idx]->setitem(1, tid);
                                extary[idx]->setsize(2);
                            } else {
                                extary[idx]->setitem(1, tid);
                            }

                            lit = offsets[idx] - idx - 1;
                            for (l = j + 1; l < numitem; l++) {
                                if (freqidx[buf[l]] != -1) {
                                    ocust[lit + freqidx[buf[l]]] = 1;
                                }
                            }
                        } else {
                            lit = offsets[idx] - idx - 1;
                            for (l = j + 1; l < numitem; l++) {
                                if (freqidx[buf[l]] != -1) {
                                    if ((++itcnt2[lit + freqidx[buf[l]]]) == 0) {
                                        write(isetfd, (char *) &buf[j], ITSZ);
                                        write(isetfd, (char *) &buf[l], ITSZ);
                                    }
                                }
                            }
                        }
                    }
                }
                DCB->get_next_trans(buf, numitem, tid, custid);
            }

            if (use_seq) {
                process_cust(fidx, fcnt, numfreq, backidx, extary, seq2, itcnt2,
                             ocust, offsets, seqfd, isetfd, use_seq);
            }
        }
        delete[] ocust;
        logger << "2-IT " << " " << endl;

        //write 2-itemsets counts to file
        int l2cnt = 0;
        if (use_seq) {
            ofd.open(cspade_args.seqf.c_str());
            if (ofd.fail()) {
                string error_message = "can't open seq file: " + cspade_args.seqf;
                throw runtime_error(error_message);
            }
            sort_get_l2(l2cnt, seqfd, ofd, backidx, freqidx,
                        numfreq, offsets, seq2, 1, global::MINSUP_ABS, cspade_args.twoseq);

            ofd.close();
            logger << "SEQ2 cnt " << l2cnt << endl;
        }
        int seqs = l2cnt;

        ofd.open(cspade_args.it2f.c_str());
        //if ((fd = open(it2fn, (O_WRONLY|O_CREAT|O_TRUNC), 0666)) < 0){
        if (ofd.fail()) {
            string error_message = "can't open it2 file: " + cspade_args.it2f;
            throw runtime_error(error_message);
        }
        sort_get_l2(l2cnt, isetfd, ofd, backidx, freqidx,
                    numfreq, offsets, itcnt2, 0, global::MINSUP_ABS, cspade_args.twoseq);
        ofd.close();
        logger << "SORT " << l2cnt << "  " << endl;

        if (use_seq) unlink("tmpseq");
        unlink("tmpiset");
        delete[] offsets;
        delete[] itcnt2;
        if (use_seq) delete[] seq2;
    }

    do_invert_db(DCB, pblk, extary, numfreq, freqidx, backidx, fidx, mincustid, maxcustid, use_seq);

    delete[] freqidx;
    delete[] backidx;

    delete DCB;
}


/**
 *
 * @return
 */
vector<string> exttpose() {
    using sequence::cspade_args;
    bool use_seq = true;
    if (cspade_args.twoseq) {
        use_seq = false;
    }

    logger << "CONF " << global::DBASE_NUM_TRANS << " " << global::DBASE_MAXITEM << " " <<
           global::DBASE_AVG_TRANS_SZ << " " << global::DBASE_AVG_CUST_SZ << endl;

    if (cspade_args.use_diff) {
        use_seq = false;
        cspade_args.num_partitions = 1;
        logger << "SEQ TURNED OFF and PARTITIONS = 1\n";
    }

    tpose(use_seq);

    vector<string> retval;
    retval.push_back(cspade_args.binf);
    retval.push_back(cspade_args.idxf);
    retval.push_back(cspade_args.it2f);
    retval.push_back(cspade_args.seqf);

    return retval;
}