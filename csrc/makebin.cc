#include <cerrno>
#include "makebin.h"

const int lineSize = 8192;
const int wdSize = 256;


void convbin(char *inBuf, std::streamsize inSize, ofstream &fout) {
    char inStr[wdSize];
    istrstream ist(inBuf, inSize);
    int it;
    while (ist >> inStr) {
        it = atoi(inStr);
        fout.write((char *) &it, ITSZ);
    }
}

void convert_bin(const string& ifname) {
    using sequence::cspade_args;

    ifstream fin(ifname.c_str());
    ofstream fout(cspade_args.binf.c_str());
    char inBuf[lineSize];
    std::streamsize inSize;
    if (!fin) {
        string error_message = "can't open ascii file: " + ifname;
        throw runtime_error(error_message);
    }
    if (!fout) {
        string error_message = "can't open binary file: " + cspade_args.binf;
        throw runtime_error(error_message);
    }

    while (fin.getline(inBuf, lineSize)) {
        inSize = fin.gcount();
        convbin(inBuf, inSize, fout);
    }
}
