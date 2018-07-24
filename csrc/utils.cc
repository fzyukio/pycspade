#include "utils.h"
std::ostringstream *_logger = new std::ostringstream();

int cmp2it(const void *a, const void *b) {
    int *ary = (int *) a;
    int *bry = (int *) b;
    if (ary[0] < bry[0]) return -1;
    else if (ary[0] > bry[0]) return 1;
    else {
        if (ary[1] < bry[1]) return -1;
        else if (ary[1] > bry[1]) return 1;
        else return 0;
    }
}

