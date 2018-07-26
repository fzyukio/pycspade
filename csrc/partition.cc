#include <cmath>
#include <cerrno>
#include <iostream>
#include <unistd.h>
#include <fcntl.h>
#include <sys/mman.h>

#include "partition.h"
#include "utils.h"

struct timeval tp;

int *DATAFD, *IDXFD, *IDXFLEN, **ITEMIDX;

void partition_alloc(const string& dataf, const string& idxf) {
    DATAFD = new int[global::num_partitions];
    IDXFD = new int[global::num_partitions];
    IDXFLEN = new int[global::num_partitions];
    ITEMIDX = new int *[global::num_partitions];
    char tmpnam[300];
    for (int i = 0; i < global::num_partitions; i++) {
        if (global::num_partitions > 1) sprintf(tmpnam, "%s.P%d", dataf.c_str(), i);
        else sprintf(tmpnam, "%s", dataf.c_str());
        DATAFD[i] = open(tmpnam, O_RDONLY);
        if (DATAFD[i] < 0) {
            string error_message = "can't open data file: " + string(tmpnam);
            throw runtime_error(error_message);
        }

        if (global::num_partitions > 1) sprintf(tmpnam, "%s.P%d", idxf.c_str(), i);
        else sprintf(tmpnam, "%s", idxf.c_str());
        IDXFD[i] = open(tmpnam, O_RDONLY);
        if (IDXFD[i] < 0) {
            string error_message = "can't open idxa file: " + string(tmpnam);
            throw runtime_error(error_message);
        }
        IDXFLEN[i] = lseek(IDXFD[i], 0, SEEK_END);
        lseek(IDXFD[i], 0, SEEK_SET);
#ifndef DEC
        ITEMIDX[i] = (int *) mmap((char *) nullptr, IDXFLEN[i], PROT_READ,
                                  MAP_PRIVATE, IDXFD[i], 0);
#else
        ITEMIDX[i] = (int *) mmap((char *)nullptr, IDXFLEN[i], PROT_READ,
                                (MAP_FILE|MAP_VARIABLE|MAP_PRIVATE),
                                IDXFD[i], 0);
#endif
        if (ITEMIDX[i] == (int *) -1) {
            throw runtime_error("MMAP ERROR:item_idx");
        }
    }
}

void partition_dealloc() {
    for (int i = 0; i < global::num_partitions; i++) {
        close(DATAFD[i]);
        close(IDXFD[i]);
        munmap((caddr_t) ITEMIDX[i], IDXFLEN[i]);
    }
    delete[] DATAFD;
    delete[] IDXFD;
    delete[] IDXFLEN;
    delete[] ITEMIDX;
}

int partition_get_blk_sz(int p) {
    return lseek(DATAFD[p], 0, SEEK_END);
}

int partition_get_max_blksz() {
    int max = 0;
    int flen;
    for (int i = 0; i < global::num_partitions; i++) {
        flen = lseek(DATAFD[i], 0, SEEK_END);
        if (max < flen) max = flen;
    }
    return max;
}

void partition_get_blk(int *MAINBUF, int p) {
    int flen = lseek(DATAFD[p], 0, SEEK_END);
    logger << "FILESZ " << flen << endl;
    lseek(DATAFD[p], 0, SEEK_SET);
    if (read(DATAFD[p], (char *) MAINBUF, flen) < 0) {
        throw runtime_error("read item1");
    }
}

int partition_get_idxsup(int it) {
    int supsz = 0;
    for (int i = 0; i < global::num_partitions; i++) {
        supsz += ITEMIDX[i][it + 1] - ITEMIDX[i][it];
    }
    return supsz;
}

int partition_get_lidxsup(int idx, int it) {
    return (ITEMIDX[idx][it + 1] - ITEMIDX[idx][it]);
}

int partition_get_idx(int idx, int it) {
    return ITEMIDX[idx][it];
}

int *partition_idx(int idx) {
    return ITEMIDX[idx];
}

void partition_read_item(int *ival, int it) {
    int ipos = 0;
    int supsz;
    for (int i = 0; i < global::num_partitions; i++) {
        supsz = ITEMIDX[i][it + 1] - ITEMIDX[i][it];
        if (supsz > 0) {
            lseek(DATAFD[i], ITEMIDX[i][it] * sizeof(int), SEEK_SET);
            if (read(DATAFD[i], (char *) &ival[ipos], supsz * sizeof(int)) < 0) {
                throw runtime_error("read item1");
            }
            ipos += supsz;
        }
    }
}

void partition_lclread_item(int *ival, int pnum, int it) {
    int supsz;
    supsz = ITEMIDX[pnum][it + 1] - ITEMIDX[pnum][it];
    if (supsz > 0) {
        lseek(DATAFD[pnum], ITEMIDX[pnum][it] * sizeof(int), SEEK_SET);
        if (read(DATAFD[pnum], (char *) ival, supsz * sizeof(int)) < 0) {
            throw runtime_error("read item1");
        }
    }
}


void partition_get_minmaxcustid(int *backidx, int numit, int pnum,
                                int &minv, int &maxv) {
    int custid, it, i, supsz;
    minv = INT_MAX;
    maxv = 0;
    for (i = 0; i < numit; i++) {
        it = backidx[i];
        supsz = ITEMIDX[pnum][it + 1] - ITEMIDX[pnum][it];
        if (supsz > 0) {
            lseek(DATAFD[pnum], ITEMIDX[pnum][it] * sizeof(int), SEEK_SET);
            read(DATAFD[pnum], (char *) &custid, sizeof(int));
            if (minv > custid) minv = custid;
            lseek(DATAFD[pnum], (supsz - 3) * sizeof(int), SEEK_CUR);
            read(DATAFD[pnum], (char *) &custid, sizeof(int));
            if (maxv < custid) maxv = custid;
        }
    }
}


//public
int *ClassInfo::CLASSCNT = nullptr;
int *ClassInfo::MINSUP = nullptr;
int *ClassInfo::TMPE = nullptr;
int *ClassInfo::TMPM = nullptr;
int *ClassInfo::TMPL = nullptr;

//private                                
int ClassInfo::fd = -1;
int *ClassInfo::classes = nullptr;
int *ClassInfo::clsaddr = nullptr;

ClassInfo::ClassInfo(bool use_class, const string& classf) {
    int i, numtrans, maxval;
    if (use_class) {
        fd = open(classf.c_str(), O_RDONLY);
        if (fd < 0) {
            throw runtime_error("ERROR: InvalidClassFile");
        }

        long fdlen = lseek(fd, 0, SEEK_END);
        clsaddr = (int *) mmap((char *) nullptr, fdlen, PROT_READ, MAP_PRIVATE, fd, 0);
        if (clsaddr == (int *) -1) {
            throw runtime_error("MMAP ERROR:classfile_idx");
        }
        // first entry contains num classes
        global::NUMCLASS = clsaddr[0];
        //input is global::NUMCLASS followed by <cid, class> pairs
        numtrans = (fdlen / sizeof(int) - 1) / 2;
        maxval = clsaddr[numtrans * 2 - 1] + 1;
        classes = new int[maxval];
        for (i = 0; i < maxval; i++) classes[i] = NOCLASS;
        for (i = 1; i < fdlen / sizeof(int); i += 2) {
            classes[clsaddr[i]] = clsaddr[i + 1];
        }
    }

    CLASSCNT = new int[global::NUMCLASS];
    TMPE = new int[global::NUMCLASS];
    TMPM = new int[global::NUMCLASS];
    TMPL = new int[global::NUMCLASS];
    MINSUP = new int[global::NUMCLASS];

    for (i = 0; i < global::NUMCLASS; i++)
        CLASSCNT[i] = 0;

    if (use_class) {
        // class frequency
        for (i = 0; i < maxval; i++)
            if (classes[i] != NOCLASS)
                CLASSCNT[classes[i]]++;
    } else CLASSCNT[0] = global::DBASE_NUM_TRANS;

    for (i = 0; i < global::NUMCLASS; i++) {
        MINSUP[i] = (int) (global::MINSUP_PER * CLASSCNT[i] + 0.5);
        if (MINSUP[i] < 1) MINSUP[i] = 1;
    }
}

ClassInfo::~ClassInfo() {
    delete[] CLASSCNT;
    delete[] MINSUP;
    delete[] TMPE;
    delete[] TMPM;
    delete[] TMPL;
    if (fd != -1) {
        long fdlen = lseek(fd, 0, SEEK_END);
        munmap((caddr_t) clsaddr, fdlen);
        close(fd);
    }
}

