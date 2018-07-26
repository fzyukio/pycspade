#ifndef __PARTITION_H_
#define __PARTITION_H_

#include "spade.h"
#include <sys/time.h>

extern struct timeval tp;
#define seconds(tm) gettimeofday(&tp,(struct timezone *)0);\
tm=tp.tv_sec+tp.tv_usec/1000000.0

extern void partition_alloc(const string& dataf, const string& idxf);

extern void partition_dealloc();

extern void partition_get_blk(int *MAINBUF, int p);

extern int partition_get_blk_sz(int p);

extern int partition_get_max_blksz();

extern int partition_get_idxsup(int it);

extern int partition_get_lidxsup(int idx, int it);

extern int partition_get_idx(int idx, int it);

extern int *partition_idx(int idx);

extern void partition_read_item(int *ival, int it);

extern void partition_lclread_item(int *ival, int pnum, int it);

extern void partition_get_minmaxcustid(int *backidx, int numit, int pnum,
                                       int &minv, int &maxv);

const int NOCLASS = -1;

class ClassInfo {
private:
    static int fd;
    static int *clsaddr;
    static int *classes;
public:
    static int *CLASSCNT;
    static int *MINSUP;
    static int *TMPE;            // temporary variables to keep support
    static int *TMPM;            // counts during intersections
    static int *TMPL;


    ClassInfo(bool use_class, const string& classf);

    ~ClassInfo();

    static int getcnt(int cls = -1) {
        if (cls == -1) {
            int sum = 0;
            for (int i = 0; i < global::NUMCLASS; i++)
                sum += CLASSCNT[i];
            return sum;
        } else return CLASSCNT[cls];
    }

    static int getcls(int idx) {
        if (fd == -1)
            return 0;
        else return classes[idx];
    }
};

#endif// __PARTITION_H_
