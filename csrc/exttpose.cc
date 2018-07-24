#include <cerrno>
#include <iostream>
#include <fstream>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/time.h>
#include <cmath>
#include <sstream>

#include "calcdb.h"
#include "ArrayT.h"
#include "exttpose.h"
#include "utils.h"


void sort_get_l2(int &l2cnt, int fd, std::ofstream &ofd, int *backidx, int *freqidx,
                 int numfreq, int *offsets, unsigned char *cntary, char use_seq, int MINSUPPORT, bool twoseq) {
    //write 2-itemsets counts to file

    int i, j, k, fcnt;
    int lit;
    long sortflen;
    int *sortary;
    int itbuf[3];

    sortflen = lseek(fd, 0, SEEK_CUR);
    if (sortflen < 0) {
        throw std::runtime_error("SEEK SEQ");
    }
    //std::cout << "SORT " << sortflen << std::endl;
    if (sortflen > 0) {
#ifdef DEC
        sortary = (int *) mmap((char *)NULL, sortflen,
                               (PROT_READ|PROT_WRITE),
                               (MAP_FILE|MAP_VARIABLE|MAP_PRIVATE), fd, 0);
#else
        sortary = (int *) mmap((char *) NULL, sortflen,
                               (PROT_READ | PROT_WRITE),
                               MAP_PRIVATE, fd, 0);
#endif
        if (sortary == (int *) -1) {
            throw std::runtime_error("SEQFd MMAP ERROR");
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

            if (fcnt >= MINSUPPORT) {
                if (twoseq) {
                    ofd.write((char *) &fcnt, ITSZ);
                } else {
                    itbuf[0] = backidx[j];
                    itbuf[1] = backidx[k];
                    itbuf[2] = fcnt;
                    ofd.write((char *) itbuf, 3 * ITSZ);
                }
                //std::cout << backidx[j] << ((use_seq)?" -> ":" ")
                //     << backidx[k] << " SUPP " << fcnt << std::endl;
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
                  int mincustid, int maxcustid, int num_partitions, char *output, bool use_diff, bool use_newformat,
                  bool use_seq) {
    double t1, t2;
    int numitem, tid, custid;
    int *buf;
    char tmpnam[300];
    int i, j, k;
    int fd;
    int idx;

    DCB->get_first_blk();
    DCB->get_next_trans(buf, numitem, tid, custid);
    int ocid;// = -1;
    for (int p = 0; p < num_partitions; p++) {
        if (num_partitions > 1) sprintf(tmpnam, "%s.P%d", output, p);
        else sprintf(tmpnam, "%s", output);
        if ((fd = open(tmpnam, (O_WRONLY | O_CREAT | O_TRUNC), 0666)) < 0) {
            throw std::runtime_error("Can't open out file");
        }

        for (i = 0; i < numfreq; i++) {
            extary[i]->reset();
        }
        //count 2-itemsets
        int plb = p * pblk + mincustid;
        int pub = plb + pblk;
        if (pub >= maxcustid) pub = maxcustid + 1;
        std::cout << "BOUNDS " << plb << " " << pub << std::endl;
        int fcnt;
        for (; !DCB->eof() && custid < pub;) {
            fcnt = 0;
            ocid = custid;
            //std::cout << "TID " << custid << " " << tid << " " << numitem << std::endl;
            while (!DCB->eof() && ocid == custid && custid < pub) {
                //for (k=0; k < numitem; k++){

                // }

                if (use_diff) {
                    //add this tid to all items not in the trans
                    k = 0;
                    for (j = 0; j < numitem; j++) {
                        if (freqidx[buf[j]] == -1) continue;

                        while (backidx[k] < buf[j]) {
                            //if ((idx = freqidx[backidx[k]]) != -1){
                            idx = k;
                            if (!use_newformat)
                                extary[idx]->add(fd, tid, use_seq, p);
                            else extary[idx]->add(fd, tid, use_seq, p, custid);
                            //}
                            k++;
                        }
                        k++; //skip over buf[j]
                    }
                    for (; k < numfreq; k++) {
                        //if ((idx = freqidx[backidx[k]]) != -1){
                        idx = k;
                        if (!use_newformat)
                            extary[idx]->add(fd, tid, use_seq, p);
                        else extary[idx]->add(fd, tid, use_seq, p, custid);
                        //}
                    }
                } else {
                    // add this tid to all items in the trans
                    for (j = 0; j < numitem; j++) {
                        idx = freqidx[buf[j]];
                        if (idx != -1) {
                            if (!use_newformat) {
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
            if (!use_newformat && use_seq) {
                for (k = 0; k < fcnt; k++) {
                    extary[fidx[k]]->setlastpos();
                    extary[fidx[k]]->setflg(0);
                }
                fcnt = 0;
            }
        }

        for (i = 0; i < numfreq; i++) {
            //std::cout << "FLUSH " << i << " " << extary[i]->lastPos << " " <<
            //   extary[i]->theSize << std::endl;
            extary[i]->flushbuf(fd, use_seq, p);
        }
        close(fd);
    }
    std::cout << "WROTE INVERT " << std::endl;
}

void tpose(char *input, char *idxfn, char *it2fn, char *seqfn, char *output, int MINSUPPORT, bool twoseq,
           int DBASE_MAXITEM, int DBASE_NUM_TRANS, long AMEM, int num_partitions, bool do_invert, bool use_newformat,
           bool use_diff, bool no_minus_off, bool do_l2, bool use_seq) {
    int i, j, l;
    int idx;
    int custid, tid, numitem, fcnt;
    std::ofstream ofd;
    double t1, t2;
    int sumsup = 0, sumdiff = 0;

    CalcDb *DCB = new CalcDb(input);

    int *itcnt = new int[DBASE_MAXITEM];
    int *ocnt = new int[DBASE_MAXITEM];
    int *itlen = new int[DBASE_MAXITEM];
    bzero((char *) itcnt, ((DBASE_MAXITEM) * ITSZ));
    //bzero((char *)ocnt, ((DBASE_MAXITEM)*ITSZ));
    for (i = 0; i < DBASE_MAXITEM; i++) ocnt[i] = -1;
    bzero((char *) itlen, ((DBASE_MAXITEM) * ITSZ));

    //count 1 items
    int *buf;
    DCB->get_first_blk();
    DCB->get_next_trans(buf, numitem, tid, custid);
    int mincustid = custid;
    while (!DCB->eof()) {
        //std::cout << custid << " " << tid << " " << numitem;
        for (j = 0; j < numitem; j++) {
            //std::cout << " " << buf[j] << std::flush;
            itlen[buf[j]]++;
            if (use_seq && ocnt[buf[j]] != custid) {
                itcnt[buf[j]]++;
                ocnt[buf[j]] = custid;
            }
            //if (buf[j] == 17) std::cout << " " << tid;
        }
        //std::cout << std::endl;
        DCB->get_next_trans(buf, numitem, tid, custid);
    }
    //std::cout << std::endl;
    int maxcustid = custid;
    std::cout << "MINMAX " << mincustid << " " << maxcustid << std::endl;

    int *freqidx = new int[DBASE_MAXITEM];
    int numfreq = 0;
    for (i = 0; i < DBASE_MAXITEM; i++) {
        if (use_seq) {
            if (itcnt[i] >= MINSUPPORT) {
                std::cout << i << " SUPP " << itcnt[i] << std::endl;
                freqidx[i] = numfreq;
                numfreq++;
            } else freqidx[i] = -1;
        } else {
            if (itlen[i] >= MINSUPPORT) {
                freqidx[i] = numfreq;
                numfreq++;
                sumsup += itlen[i];
                sumdiff += (DBASE_NUM_TRANS - itlen[i]);
            } else freqidx[i] = -1;
        }
        //if (i == 17) std::cout << " 17 SUP " << itlen[17] << std::endl;
    }
    int *backidx = new int[numfreq];
    numfreq = 0;
    for (i = 0; i < DBASE_MAXITEM; i++) {
        if (use_seq) {
            if (itcnt[i] >= MINSUPPORT)
                backidx[numfreq++] = i;
        } else {
            if (itlen[i] >= MINSUPPORT)
                backidx[numfreq++] = i;
        }
    }

    std::cout << "numfreq " << numfreq << " :  " << " SUMSUP SUMDIFF = " << sumsup << " " << sumdiff << std::endl;

    if (numfreq == 0) return;

    int extarysz = AMEM / numfreq;
    extarysz /= sizeof(int);
    std::cout << "EXTRARYSZ " << extarysz << std::endl;
    if (extarysz < 2) extarysz = 2;
    ArrayT **extary = new ArrayT *[numfreq];
    for (i = 0; i < numfreq; i++) {
        extary[i] = new ArrayT(extarysz, num_partitions);
    }


    char tmpnam[300];
    int plb, pub, pblk;
    pblk = (int) ceil(((double) (maxcustid - mincustid + 1)) / num_partitions);
    if (do_invert) {
        if (num_partitions > 1) {
            DCB->get_first_blk();
            DCB->get_next_trans(buf, numitem, tid, custid);
        }
        for (j = 0; j < num_partitions; j++) {
            //construct offsets for 1-itemsets
            if (num_partitions > 1) {
                sprintf(tmpnam, "%s.P%d", idxfn, j);
                plb = j * pblk + mincustid;
                pub = plb + pblk;
                if (pub > maxcustid) pub = maxcustid + 1;
                bzero((char *) itcnt, ((DBASE_MAXITEM) * ITSZ));
                //bzero((char *)ocnt, ((DBASE_MAXITEM)*ITSZ));
                for (i = 0; i < DBASE_MAXITEM; i++) ocnt[i] = -1;
                bzero((char *) itlen, ((DBASE_MAXITEM) * ITSZ));
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
            } else sprintf(tmpnam, "%s", idxfn);
            //std::cout << "100 VAL " << itcnt[100] << std::endl;
            std::cout << "OPENED " << tmpnam << std::endl;
            ofd.open(tmpnam);
            if (!ofd) {
                throw std::runtime_error("Can't open out file");
            }

            int file_offset = 0;
            int null = -1;
            for (i = 0; i < DBASE_MAXITEM; i++) {
                //if (i == 17) std::cout << "LIDX " << i << " " << itlen[i] << std::endl;
                if (freqidx[i] != -1) {
                    ofd.write((char *) &file_offset, ITSZ);
                    extary[freqidx[i]]->set_offset(file_offset, j);
                    if (use_seq) {
                        if (use_newformat) file_offset += (2 * itlen[i]);
                        else file_offset += (2 * itcnt[i] + itlen[i]);
                    } else {
                        if (use_diff) file_offset += (DBASE_NUM_TRANS - itlen[i]);
                        else file_offset += itlen[i];
                    }
                } else if (no_minus_off) {
                    ofd.write((char *) &file_offset, ITSZ);
                } else ofd.write((char *) &null, ITSZ);
                //std::cout << "OFF " << i <<" " << file_offset << std::endl;
            }
            std::cout << "OFF " << i << " " << file_offset << std::endl;
            ofd.write((char *) &file_offset, ITSZ);
            ofd.close();
        }
    }

    delete[] ocnt;
    delete[] itlen;
    delete[] itcnt;

    std::cout << "Wrote Offt " << std::endl;

    int *fidx = new int[numfreq];
    if (fidx == NULL) {
        throw std::runtime_error("Can't alloc fidx");
        exit(errno);
    }

    int ocid = -1;
    if (do_l2) {
        int seqfd, isetfd;
        if (use_seq) {
            if ((seqfd = open("tmpseq", (O_RDWR | O_CREAT | O_TRUNC), 0666)) < 0) {
                throw std::runtime_error("Can't open out file");
            }
        }

        if ((isetfd = open("tmpiset", (O_RDWR | O_CREAT | O_TRUNC), 0666)) < 0) {
            throw std::runtime_error("Can't open out file");
        }

        unsigned char *seq2;
        if (use_seq) {
            seq2 = new unsigned char[numfreq * numfreq];
            if (seq2 == NULL) {
                throw std::runtime_error("SEQ MMAP ERROR");
            }
            //for (i=0; i < numfreq*numfreq; i++) seq2[i] = 0;
            bzero((char *) seq2, numfreq * numfreq * sizeof(unsigned char));
        }

        unsigned char *itcnt2 = new unsigned char[(numfreq * (numfreq - 1) / 2)];
        if (itcnt2 == NULL) {
            throw std::runtime_error("ITCNT MMAP ERROR");
        }
        bzero((char *) itcnt2, (numfreq * (numfreq - 1) / 2) * sizeof(unsigned char));
        //for (i=0; i < numfreq*(numfreq-1)/2; i++) itcnt2[i] = 0;
        char *ocust = new char[(numfreq * (numfreq - 1) / 2)];
        if (ocust == NULL) {
            throw std::runtime_error("OCUSt MMAP ERROR");
        }
        bzero((char *) ocust, (numfreq * (numfreq - 1) / 2) * sizeof(char));
        int *offsets = new int[numfreq];
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
                                //extary[idx]->add(isetfd,tid,use_seq,0);
                                //extary[idx]->add(isetfd,tid,use_seq,0);
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
        std::cout << "2-IT " << " " << std::endl;

        //write 2-itemsets counts to file
        int l2cnt = 0;
        if (use_seq) {
            ofd.open(seqfn);
            if (ofd.fail()) {
                throw std::runtime_error("Can't open seq file");
            }
            sort_get_l2(l2cnt, seqfd, ofd, backidx, freqidx,
                        numfreq, offsets, seq2, 1, MINSUPPORT, twoseq);

            ofd.close();
            std::cout << "SEQ2 cnt " << l2cnt << std::endl;
        }
        int seqs = l2cnt;

        ofd.open(it2fn);
        //if ((fd = open(it2fn, (O_WRONLY|O_CREAT|O_TRUNC), 0666)) < 0){
        if (ofd.fail()) {
            throw std::runtime_error("Can't open it2 file");
        }
        sort_get_l2(l2cnt, isetfd, ofd, backidx, freqidx,
                    numfreq, offsets, itcnt2, 0, MINSUPPORT, twoseq);
        ofd.close();
        std::cout << "SORT " << l2cnt << "  " << std::endl;

        if (use_seq) unlink("tmpseq");
        unlink("tmpiset");
        delete[] offsets;
        delete[] itcnt2;
        if (use_seq) delete[] seq2;
    }

    if (do_invert) {
        do_invert_db(DCB, pblk, extary, numfreq, freqidx, backidx, fidx, mincustid, maxcustid, num_partitions, output,
                     use_diff, use_newformat, use_seq);
    }

    delete[] freqidx;
    delete[] backidx;

    delete DCB;
}


/**
 *
 * @return
 */
std::string _exttpose(const std::string& dbname, int num_partitions, double min_support, bool twoseq, bool use_diff, bool do_l2,
              bool do_invert, bool use_newformat, int maxmem, bool no_minus_off) {
    char input[300];       //input file name
    char output[300];      //output file name
    char idxfn[300];
    char inconfn[300];
    char it2fn[300];
    char seqfn[300];
    bool use_seq = true;

    int DBASE_NUM_TRANS; //tot trans for assoc, num cust for sequences
    int DBASE_MAXITEM;   //number of items
    float DBASE_AVG_TRANS_SZ; //avg trans size
    float DBASE_AVG_CUST_SZ = 0; //avg cust size for sequences
    int DBASE_TOT_TRANS; //tot trans for sequences


    sprintf(input, "%s.data", dbname.c_str());
    sprintf(inconfn, "%s.conf", dbname.c_str());
    sprintf(output, "%s.tpose", dbname.c_str());
    sprintf(idxfn, "%s.idx", dbname.c_str());
    sprintf(it2fn, "%s.2it", dbname.c_str());
    sprintf(seqfn, "%s.2seq", dbname.c_str());

    std::cout << "input = " << input << std::endl;
    std::cout << "inconfn = " << inconfn << std::endl;
    std::cout << "output = " << output << std::endl;
    std::cout << "idxfn = " << idxfn << std::endl;
    std::cout << "it2fn = " << it2fn << std::endl;
    std::cout << "seqfn = " << seqfn << std::endl;

    std::cout << "num_partitions = " << num_partitions << std::endl;
    std::cout << "min_support = " << min_support << std::endl;
    std::cout << "twoseq = " << twoseq << std::endl;
    std::cout << "use_diff = " << use_diff << std::endl;
    std::cout << "do_l2 = " << do_l2 << std::endl;
    std::cout << "do_invert = " << do_invert << std::endl;
    std::cout << "use_newformat = " << use_newformat << std::endl;
    std::cout << "maxmem = " << maxmem << std::endl;
    std::cout << "no_minus_off = " << no_minus_off << std::endl;

    if (twoseq) {
        use_seq = false;
    }

    long AMEM = maxmem * MEG;

    int c = open(inconfn, O_RDONLY);
    if (c < 0) {
        throw std::runtime_error("ERROR: invalid conf file");
    }
    if (use_seq) {
        read(c, (char *) &DBASE_NUM_TRANS, ITSZ);
        read(c, (char *) &DBASE_MAXITEM, ITSZ);
        read(c, (char *) &DBASE_AVG_CUST_SZ, sizeof(float));
        read(c, (char *) &DBASE_AVG_TRANS_SZ, sizeof(float));
        read(c, (char *) &DBASE_TOT_TRANS, ITSZ);
    } else {
        read(c, (char *) &DBASE_NUM_TRANS, ITSZ);
        read(c, (char *) &DBASE_MAXITEM, ITSZ);
        read(c, (char *) &DBASE_AVG_TRANS_SZ, sizeof(float));
    }
    std::cout << "CONF " << DBASE_NUM_TRANS << " " << DBASE_MAXITEM << " " <<
         DBASE_AVG_TRANS_SZ << " " << DBASE_AVG_CUST_SZ << std::endl;

    close(c);

    if (use_diff) {
        use_seq = 0;
        num_partitions = 1;
        std::cout << "SEQ TURNED OFF and PARTITIONS = 1\n";
    }
    int MINSUPPORT = (int) (min_support * DBASE_NUM_TRANS + 0.5);

    //ensure that support is at least 2
    if (!twoseq && MINSUPPORT < 1) MINSUPPORT = 1;
    std::cout << "MINSUPPORT " << MINSUPPORT << " " << DBASE_NUM_TRANS << std::endl;

    tpose(input, idxfn, it2fn, seqfn, output, MINSUPPORT, twoseq, DBASE_MAXITEM, DBASE_NUM_TRANS, AMEM, num_partitions,
          do_invert, use_newformat, use_diff, no_minus_off, do_l2, use_seq);

    std::ostringstream os;
    os << output << '\t' << idxfn << '\t' << it2fn << '\t' << seqfn;

    return os.str();
}