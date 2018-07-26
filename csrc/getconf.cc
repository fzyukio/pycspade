#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fcntl.h>
#include <cmath>

#include "calcdb.h"
#include "getconf.h"

#include "utils.h"

//using namespace std;

string create_conf(bool assoc) {
    using sequence::cspade_args;

    int DBASE_NUM_TRANS = 0;
    int DBASE_MAXITEM = 0;
    int DBASE_NUM_CUST = 0;
    int DBASE_MINTRANS = 0;
    int DBASE_MAXTRANS = 0;
    double DBASE_AVG_TRANS_SZ = 0;
    double DBASE_AVG_CUST_SZ = 0;

    int i;

    int custid, tid, nitem;
    int *buf;
    int oldcustid = -1;
    int oldtcnt = 0;
    int tsizesum = 0;
    int tcustsum = 0;
    int tsizesq = 0;
    int maxnitem = 0;

    auto *DCB = new CalcDb(cspade_args.binf);
    DCB->get_first_blk();
    DCB->get_next_trans(buf, nitem, tid, custid);
    DBASE_MINTRANS = custid;
    while (!DCB->eof()) {
        DBASE_MAXTRANS = custid;
        if (!assoc) {
            if (oldcustid != custid) {
                tcustsum += DBASE_NUM_TRANS - oldtcnt;
                oldtcnt = DBASE_NUM_TRANS;
                DBASE_NUM_CUST++;
                oldcustid = custid;
            }
        }
        DBASE_NUM_TRANS++;
        tsizesum += nitem;
        if (nitem > maxnitem) maxnitem = nitem;

        tsizesq += (nitem * nitem);
        for (i = 0; i < nitem; i++)
            if (buf[i] > DBASE_MAXITEM) DBASE_MAXITEM = buf[i];
        DCB->get_next_trans(buf, nitem, tid, custid);
    }
    tcustsum += DBASE_NUM_TRANS - oldtcnt;
    DBASE_MAXITEM++;

    if (!assoc) DBASE_AVG_CUST_SZ = (1.0 * tcustsum) / DBASE_NUM_CUST;
    DBASE_AVG_TRANS_SZ = (1.0 * tsizesum) / DBASE_NUM_TRANS;
    double trans_sq_avg = (1.0 * tsizesq) / DBASE_NUM_TRANS;
    double stddev = sqrt(trans_sq_avg -
                         (DBASE_AVG_TRANS_SZ * DBASE_AVG_TRANS_SZ));


    //write config info to new file
    int conffd;
    if ((conffd = open(cspade_args.conf.c_str(), (O_WRONLY | O_CREAT), 0666)) < 0) {
        string error_message = "can't open conf file: " + cspade_args.conf;
        throw runtime_error(error_message);
    }
    if (assoc) {
        write(conffd, (char *) &DBASE_NUM_TRANS, ITSZ);
        write(conffd, (char *) &DBASE_MAXITEM, ITSZ);
        write(conffd, (char *) &DBASE_AVG_TRANS_SZ, sizeof(double));
        write(conffd, (char *) &DBASE_MINTRANS, ITSZ);
        write(conffd, (char *) &DBASE_MAXTRANS, ITSZ);
    } else {
        write(conffd, (char *) &DBASE_NUM_CUST, ITSZ);
        write(conffd, (char *) &DBASE_MAXITEM, ITSZ);
        write(conffd, (char *) &DBASE_AVG_CUST_SZ, sizeof(double));
        write(conffd, (char *) &DBASE_AVG_TRANS_SZ, sizeof(double));
        write(conffd, (char *) &DBASE_NUM_TRANS, ITSZ);
        write(conffd, (char *) &DBASE_MINTRANS, ITSZ);
        write(conffd, (char *) &DBASE_MAXTRANS, ITSZ);
    }

    close(conffd);
    delete DCB;
    return string(cspade_args.conf);
}








