#ifndef __ARRAY_T_H
#define __ARRAY_T_H

#include <cerrno>
#include <fcntl.h>
#include <unistd.h>
#include <cstdlib>
#include <sys/types.h>

#include "utils.h"

using namespace std;

class ArrayT {
protected:
    int *theArray;
    char theFlg;
    int lastPos;
    unsigned int theSize;
    unsigned int totSize;
    long *offset;
public:

    ArrayT(int sz, int npart = 1);

    ~ArrayT();

    int operator[](unsigned int index) {
        return theArray[index];
    };

    char flg() {
        return theFlg;
    }

    void setflg(char flg) {
        theFlg = flg;
    }

    int lastpos() {
        return lastPos;
    }

    //to be used ony for use_seq
    void setlastpos() {
        theArray[lastPos + 1] = theSize - lastPos - 2;
        lastPos = theSize;
    }

    long get_offset(int pos = 0) {
        return offset[pos];
    }

    void set_offset(long off, int pos = 0) {
        offset[pos] = off;
    }

    int totsize() {
        return totSize;
    }

    void reset() {
        theSize = 0;
        lastPos = 0;
        theFlg = 0;
    }

    int *array() {
        return theArray;
    }

    int size() {
        return theSize;
    }

    void setsize(int size) {
        theSize = size;
    }

    void setitem(int pos, int item) {
        theArray[pos] = item;
    }

    void additem(int item) {
        theArray[theSize] = item;
        theSize++;
    }

    void flushbuf(int fd, int use_seq, int pos = 0) {
        lseek(fd, offset[pos] * sizeof(int), SEEK_SET);
        int wblk = theSize;
        if (wblk > 0) {
            int res = ::write(fd, (char *) theArray, wblk * sizeof(int));
            if (res < wblk * sizeof(int)) {
                throw runtime_error("Error writing");
            }
            offset[pos] += wblk;
        }
        theSize = 0;
    }

    void add(int fd, int item, int use_seq, int pos, int custid = -1) {
        if (use_seq) {
            if (theSize + 2 > totSize) {
                flushbuf(fd, use_seq, pos);
            }
            theArray[theSize++] = custid;
        } else {
            if (theSize + 1 > totSize) {
                flushbuf(fd, use_seq, pos);
            }
        }
        theArray[theSize++] = item;
    }
};

#endif //__ARRAY_T_H


