#ifndef __MAKEBIN_H
#define __MAKEBIN_H

#include "utils.h"

//using namespace std;

extern void convbin(char *inBuf, std::streamsize inSize, ofstream &fout);

extern void convert_bin(const string& ifname);

#endif //__MAKEBIN_H