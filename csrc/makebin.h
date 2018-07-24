#ifndef __MAKEBIN_H
#define __MAKEBIN_H

#include <fstream>

//using namespace std;

extern void convbin(char *inBuf, int inSize, std::ofstream &fout);

extern void convert_bin(const std::string& ifname, const std::string& ofname);

#endif //__MAKEBIN_H