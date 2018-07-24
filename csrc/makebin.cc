#include <cerrno>
#include <iostream>
#include <fstream>
#include <strstream>
#include "makebin.h"

#include "utils.h"
const int lineSize = 8192;
const int wdSize = 256;

//using namespace std;


void convbin(char *inBuf, int inSize, std::ofstream &fout) {
    char inStr[wdSize];
    std::istrstream ist(inBuf, inSize);
    int it;
    while (ist >> inStr) {
        it = atoi(inStr);
        fout.write((char *) &it, ITSZ);
    }
}

void convert_bin(const std::string& ifname, const std::string& ofname) {
    std::ifstream fin(ifname.c_str());
    std::ofstream fout(ofname.c_str());
    char inBuf[lineSize];
    int inSize;
    if (!fin) {
        throw std::runtime_error("cannot open in file");
    }
    if (!fout) {
        throw std::runtime_error("cannot open out file");
    }

    while (fin.getline(inBuf, lineSize)) {
        inSize = fin.gcount();
        //std::cout << "IN SIZE " << inSize << std::endl;
        convbin(inBuf, inSize, fout);
    }
}
